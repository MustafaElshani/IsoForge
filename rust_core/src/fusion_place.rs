// Read parts the aligner left unaligned, placed on annotated exons.
//
// - A piece of `min_len` or more unaligned read bases beyond an alignment's end in an annotated exon is read outward
//   and placed step by step: at an exon boundary on the exons that start there, else anywhere in an annotated exon.
// - An end is searched only when another read's alignment ends within `window` bp of it with unaligned bases beyond.
// - A step on the aligned gene's own next exon extends the alignment; one on another gene is one more part of the read.
// - A poly(A) tail, a simple repeat, recurring words and a repeat are not placed, nor is a piece that lies in the
//   genome beyond the alignment or a part whose gene the read's bases before it continue into.
// - Parts under the mapping-quality floor and short soft clips are placed by class, over the reads that share an end.

use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering::Relaxed};

use rayon::prelude::*;
use rust_htslib::bam::{self, record::Cigar, Read};
use rustc_hash::{FxHashMap, FxHashSet};

use coitrees::{COITree, Interval, IntervalTree};

use crate::fusion::ChimericSegment;
use crate::genome_ref::GenomeRef;
use crate::post_processing::PostDB;

#[derive(Debug, Clone)]
pub struct PlaceOptions {
    /// (dataset, BAM) pairs, the datasets named as in the chimeric segments.
    pub bams: Vec<(String, PathBuf)>,
    /// Unaligned read bases looked at, and placed bases a part needs.
    pub min_len: usize,
    /// Bases of a piece compared, from the alignment outward.
    pub max_len: usize,
    /// How far short of an annotated exon boundary an alignment may end (bp); it runs past one by `free` bases at most.
    pub window: i64,
    /// Bases free at the start of the comparison.
    pub free: usize,
    /// Score (match 2, mismatch -3, inserted or deleted base -4) a part placed at an exon boundary needs.
    pub min_score: i32,
    /// The score a part entering an exon anywhere needs.
    pub min_score_free: i32,
    /// Score per placed base a part needs.
    pub min_score_per_base: f64,
    pub margin: i32,
    pub min_mapq: u8,
    /// Read bases two consecutive parts may share or leave between them for the caller to join them (`fusion.rs`).
    pub max_query_gap: usize,
}

impl Default for PlaceOptions {
    fn default() -> Self {
        Self { bams: Vec::new(), min_len: 25, max_len: 150, window: 10, free: 3, min_score: 32, min_score_free: 40, min_score_per_base: 0.9, margin: 8, min_mapq: 20, max_query_gap: 30 }
    }
}

const K: usize = 10;
/// Bases of a target, from its boundary, whose k-mers are indexed.
const SEED: usize = 48;
/// A part the aligner could not place among copies (mapping quality under the floor) is named with the other reads'
/// parts at the same junction end: every copy each part fits is scored, the scores add over the reads, and the copy
/// that leads by `CLASS_MARGIN` with `CLASS_MIN_READS` reads fitting it names them all.
const CLASS_MARGIN: i32 = 30;
const CLASS_MIN_READS: usize = 2;
/// The margin the leading gene needs over the runner-up for each read in the class.
const CLASS_MARGIN_PER_READ: i32 = 10;

/// Short parts placed as a class: primary alignments ending in an annotated exon with a soft clip of `SHORT_CLIP_MIN`
/// to `SHORT_CLIP_MAX` bases; ends within `SHORT_WIN` bp on one side form a class. Its pieces (the first `SHORT_QMAX`
/// bases) name a partner when `SHORT_FIT` pieces and `SHORT_SHARE` of them fit its target (`SHORT_FIT_PER_BASE` a
/// base), the pooled score leads another gene's target by `SHORT_MARGIN` and the pieces are `SHORT_MIN` bases or more.
pub const SHORT_MIN: usize = 14;
const SHORT_CLIP_MIN: usize = 8;
const SHORT_CLIP_MAX: usize = 40;
const SHORT_WIN: i64 = 5;
const SHORT_QMAX: usize = 26;
const SHORT_FIT: usize = 3;
const SHORT_SHARE: f64 = 0.8;
const SHORT_FIT_PER_BASE: f64 = 1.4;
const SHORT_MARGIN: i32 = 24;
/// Pieces of a class scored for its decision, evenly spaced over it; every piece is then aligned to the target chosen.
const SHORT_SAMPLE: usize = 24;
/// Candidate targets a class is scored on: those with most seed hits, after targets of one gene alike over their
/// first `SHORT_HEAD` bases are taken once.
const SHORT_CANDIDATES: usize = 4000;
const SHORT_HEAD: usize = 40;
/// Words in more targets than this seed nothing.
const SHORT_BUCKET_MAX: usize = 1500;

/// A soft clip of `SHORT_CLIP_MIN` to `SHORT_CLIP_MAX` bases at an end of a primary alignment: the end's side on the
/// genome and position, the clip's bases on the record and their place [ca, cb), the read (name hash, the record's
/// position), and the strand of the read's transcript when the read gives it (`read_strand`: 1 plus, -1 minus, 0).
struct ShortEnd { high: bool, pos: i64, clip: Vec<u8>, ca: usize, cb: usize, name: u64, rec_pos: i64, strand: i8, traced: bool }

/// A class's piece placed: the target, the alignment (score, piece bases, target bases), the piece's bases on the
/// record, the strand the class was read on, and the class's size and margin.
struct ShortHit { target: u32, aln: (i32, usize, usize, usize, usize), ra: usize, rb: usize, plus: bool, pieces: u32, excess: i32 }

/// The read-name hash the short classes key reads by.
fn name_hash(name: &[u8]) -> u64 {
    use std::hash::{Hash, Hasher};
    let mut h = rustc_hash::FxHasher::default();
    name.hash(&mut h);
    h.finish()
}

/// The canonical 10-mers of a sequence (the lower of each word's code and its reverse complement's).
fn canon10(s: &[u8]) -> Vec<u32> {
    let mut out = Vec::new();
    if s.len() < 10 { return out; }
    for w in s.windows(10) {
        let (mut f, mut r, mut ok) = (0u32, 0u32, true);
        for (i, &b) in w.iter().enumerate() {
            let Some(c) = code(b) else { ok = false; break };
            f = (f << 2) | c;
            r |= (3 - c) << (2 * i as u32);
        }
        if ok { out.push(f.min(r)); }
    }
    out
}

/// A 3' piece cut where a poly(A) tail starts: three A or more, at most one other base, then four A or more.
fn short_tail_cut(s: &[u8]) -> usize {
    let n = s.len();
    let mut i = 0;
    while i < n {
        if s[i] == b'A' {
            let mut j = i;
            while j < n && s[j] == b'A' { j += 1; }
            if j - i >= 7 { return i; }
            if j - i >= 3 && j < n {
                let mut k = j + 1;
                while k < n && s[k] == b'A' { k += 1; }
                if k - (j + 1) >= 4 { return i; }
            }
            i = j;
        } else { i += 1; }
    }
    n
}

/// A short piece of low complexity: three in four of its bases one base, two distinct bases at most, or a tandem
/// repeat, with no more distinct three-base words than half its words ((GGC)n, (CAG)n).
fn short_low_complexity(s: &[u8]) -> bool {
    if s.len() < 8 { return true; }
    let mut c = [0usize; 256];
    for &b in s { c[b as usize] += 1; }
    let distinct = c.iter().filter(|&&x| x > 0).count();
    let words: FxHashSet<&[u8]> = s.windows(3).collect();
    4 * c.iter().max().copied().unwrap_or(0) >= 3 * s.len() || distinct <= 2 || 2 * words.len() <= s.len() - 2
}

/// The strand of a read's transcript by its poly(A) tail: ten A's in the twelve bases after the alignment's high end
/// (plus, 1), or ten T's in the twelve before its low end (minus, -1); 0 with neither or both. `sq` is the record's
/// sequence (the reference's orientation), `lead` and `trail` its soft clips.
fn tail_strand(sq: &[u8], lead: usize, trail: usize) -> i8 {
    let n = sq.len();
    if lead + trail > n { return 0; }
    let run = |w: &[u8], b: u8| w.iter().filter(|&&x| x.to_ascii_uppercase() == b).count() >= 10;
    let a = trail >= 12 && run(&sq[n - trail..n - trail + 12], b'A');
    let t = lead >= 12 && run(&sq[lead - 12..lead], b'T');
    match (a, t) { (true, false) => 1, (false, true) => -1, _ => 0 }
}

/// The strand of a read's transcript: the `ts` tag (`+` when the transcript runs as the read aligns), else the poly(A)
/// tail (`tail_strand`).
fn read_strand(rec: &bam::Record, sq: &[u8], lead: usize, trail: usize) -> i8 {
    if let Ok(bam::record::Aux::Char(c)) = rec.aux(b"ts") {
        if c == b'+' || c == b'-' { return if (c == b'+') != rec.is_reverse() { 1 } else { -1 }; }
    }
    tail_strand(sq, lead, trail)
}

/// A low-MAPQ part of a read, on the primary record: the part in `known`, its bases on the record and their place,
/// the record's strand and length, the read's name, the end of the primary alignment it lies beyond, whether it lies
/// to the right of the alignment on the record, and the aligned bases next to it (`BACK` at most), for the copy test.
struct LowPart { known: usize, seq: Vec<u8>, ra: usize, rev: bool, len: usize, name: String, pos: i64, right: bool, inner: Vec<u8> }

/// An exon a sequence fits, from `ExonStore::candidates`.
struct Cand { gene: u32, chrom: u32, lo: i64, hi: i64, score: i32, q0: usize, q1: usize, plus: bool }

/// An annotated exon boundary read outward: `forward` from an exon start 3'-ward (transcript orientation), else from
/// an exon end 5'-ward (the transcript's sequence reversed). `blocks` are the genomic intervals read, in reading
/// order.
struct Target {
    gene: u32,
    chrom: u32,
    plus: bool,
    forward: bool,
    blocks: Vec<(i64, i64)>,
    seq: Vec<u8>,
    /// read from a transcript's own first base 3'-ward, or from its last 5'-ward (in some annotated model)
    tx_end: bool,
}

impl Target {
    /// Reading a block runs up the genome when the transcript does and the target is read 3'-ward, or neither.
    fn ascending(&self) -> bool { self.forward == self.plus }

    /// Genomic blocks (ascending) of the target's bases [t0, t1).
    fn cover(&self, t0: usize, t1: usize) -> Vec<(i64, i64)> {
        let mut out = Vec::new();
        let mut at = 0usize;
        for &(lo, hi) in &self.blocks {
            let n = (hi - lo + 1) as usize;
            let (a, b) = (t0.max(at), t1.min(at + n));
            if a < b {
                let (x, y) = ((a - at) as i64, (b - at) as i64);
                out.push(if self.ascending() { (lo + x, lo + y - 1) } else { (hi - y + 1, hi - x) });
            }
            at += n;
        }
        out.sort_unstable();
        out
    }
}

fn code(b: u8) -> Option<u32> {
    match b { b'A' => Some(0), b'C' => Some(1), b'G' => Some(2), b'T' => Some(3), _ => None }
}

/// Each k-mer of a sequence (k <= 15) with its position, those holding a base other than A, C, G, T left out.
fn kmers_of(s: &[u8], k: usize) -> impl Iterator<Item = (usize, u32)> + '_ {
    let mask = (1u32 << (2 * k)) - 1;
    let (mut x, mut valid) = (0u32, 0usize);
    s.iter().enumerate().filter_map(move |(i, &b)| {
        match code(b) {
            Some(c) => { x = ((x << 2) | c) & mask; valid += 1; }
            None => { valid = 0; }
        }
        (valid >= k).then(|| (i + 1 - k, x))
    })
}

fn kmers(s: &[u8]) -> impl Iterator<Item = (usize, u32)> + '_ { kmers_of(s, K) }

fn complement(b: u8) -> u8 {
    match b { b'A' => b'T', b'C' => b'G', b'G' => b'C', b'T' => b'A', x => x }
}

/// Whether the bases either side of a seed (24 each side, on the seed's diagonal or one beside it) match more often
/// than chance: a quarter of the `w` compared bases and three standard deviations more. `pair(j, shift)` gives the
/// piece's base at `j` and the other sequence's on the diagonal `shift` beside the seed's, None outside them.
fn flanks_agree(seed: usize, k: usize, n: usize, pair: impl Fn(usize, i64) -> Option<(u8, u8)>) -> bool {
    N_FLANK.fetch_add(1, Relaxed);
    // first the eight bases either side on the seed's own diagonal: of twelve or more compared, seven must match
    let (mut w, mut m) = (0u32, 0u32);
    for j in (seed.saturating_sub(8)..seed).chain(seed + k..(seed + k + 8).min(n)) {
        if let Some((a, b)) = pair(j, 0) { w += 1; if a == b { m += 1; } }
    }
    if w >= 12 && m < 7 { return false; }
    let side = |from: usize, to: usize| -> (u32, u32) {
        (-1..=1).map(|shift| {
            let (mut w, mut m) = (0u32, 0u32);
            for j in from..to {
                if let Some((a, b)) = pair(j, shift) { w += 1; if a == b { m += 1; } }
            }
            (m, w)
        }).max().map_or((0, 0), |(m, w)| (w, m))
    };
    let (wl, ml) = side(seed.saturating_sub(24), seed);
    let (wr, mr) = side(seed + k, (seed + k + 24).min(n));
    let (w, m) = (wl + wr, ml + mr);
    w < 8 || m as f64 >= 0.25 * w as f64 + 3.0 * (0.1875 * w as f64).sqrt()
}

/// The best score of a prefix of `q` (up to `free` leading bases skipped) aligned anywhere in `win`, both read in
/// the same direction, and the prefix's length.
#[cfg(test)]
fn own_locus_score(q: &[u8], win: &[u8], free: usize) -> (i32, usize) { locus_score(q, win, free) }

/// The same, with `free` leading bases of `q` skipped at no cost: with all of them free, the best local alignment of
/// any part of `q`.
fn locus_score(q: &[u8], win: &[u8], free: usize) -> (i32, usize) {
    const NEG: i32 = i32::MIN / 4;
    let n = q.len();
    let mut best = (0, 0);
    let mut prev = vec![NEG; n + 1];
    let mut cur = vec![NEG; n + 1];
    for i in 0..=free.min(n) { prev[i] = 0; }
    for &w in win {
        cur[0] = 0;
        for i in 1..=n {
            let d = prev[i - 1] + if q[i - 1] == w { MATCH } else { MISMATCH };
            let mut x = d.max(prev[i] + GAP).max(cur[i - 1] + GAP);
            if i <= free { x = x.max(0); }
            cur[i] = x;
            if x > best.0 { best = (x, i); }
        }
        std::mem::swap(&mut prev, &mut cur);
    }
    best
}

/// A 12-mer and its reverse complement as one word (the smaller code): a sequence read on either strand.
fn canonical12(s: &[u8]) -> impl Iterator<Item = u32> + '_ {
    kmers_of(s, 12).map(|(_, k)| {
        let mut rc = 0u32;
        let mut x = k;
        for _ in 0..12 { rc = (rc << 2) | (3 - (x & 3)); x >>= 2; }
        k.min(rc)
    })
}

/// 12-mers that recur in the unaligned pieces of a chromosome's reads, counted over the pieces' first `ARTEFACT_SPAN`
/// bases at the alignment (each word once a piece); a word in 0.5 % of the pieces (20 at least) is one.
#[derive(Default)]
struct Artefacts { counts: FxHashMap<u32, u32>, pieces: u32, words: FxHashSet<u32> }

const ARTEFACT_SPAN: usize = 60;

impl Artefacts {
    fn add(&mut self, piece: &[u8], scratch: &mut Vec<u32>) {
        scratch.clear();
        scratch.extend(canonical12(piece));
        scratch.sort_unstable(); scratch.dedup();
        for &k in scratch.iter() { *self.counts.entry(k).or_insert(0) += 1; }
        self.pieces += 1;
    }

    fn finish(&mut self) {
        let need = ((self.pieces as f64 * 0.005).ceil() as u32).max(20);
        self.words = self.counts.iter().filter(|&(_, &n)| n >= need).map(|(&k, _)| k).collect();
        self.counts = FxHashMap::default();
    }

    /// Whether a piece (its bases at the alignment) carries two of the recurring words.
    fn holds(&self, piece: &[u8]) -> bool {
        !self.words.is_empty() && canonical12(&piece[..piece.len().min(ARTEFACT_SPAN)]).filter(|k| self.words.contains(k)).take(2).count() >= 2
    }
}

/// The bases of a piece before its poly(A) tail (read as A or T): up to the first run of ten A, or T, or the first 24
/// bases of which seven in ten are.
fn before_tail(s: &[u8]) -> usize {
    let mut cut = s.len();
    for base in [b'A', b'T'] {
        let mut run = 0usize;
        for (i, &b) in s.iter().enumerate() {
            run = if b == base { run + 1 } else { 0 };
            if run >= 10 { cut = cut.min(i + 1 - run); break; }
        }
        let mut n = s.iter().take(24).filter(|&&b| b == base).count();
        for i in 0..s.len().saturating_sub(23) {
            if 10 * n >= 7 * 24 { cut = cut.min(i); break; }
            if i + 24 < s.len() { n += usize::from(s[i + 24] == base); n -= usize::from(s[i] == base); }
        }
    }
    cut
}

/// Whether a sequence is of low complexity (a poly(A) tail, a simple repeat), judged on its first 48 bases: seven in
/// ten are A, or T; or twelve A, or T, stand in a row; or fewer than a third of the three-base words are distinct.
fn low_complexity(s: &[u8]) -> bool {
    let s = &s[..s.len().min(48)];
    if s.len() < 12 { return false; }
    for base in [b'A', b'T'] {
        if 10 * s.iter().filter(|&&b| b == base).count() >= 7 * s.len() { return true; }
        let (mut run, mut longest) = (0usize, 0usize);
        for &b in s { run = if b == base { run + 1 } else { 0 }; longest = longest.max(run); }
        if longest >= 12 { return true; }
    }
    let mut seen = [false; 64];
    let (mut words, mut distinct) = (0u32, 0u32);
    for (_, k) in kmers_of(s, 3) {
        words += 1;
        if !seen[k as usize] { seen[k as usize] = true; distinct += 1; }
    }
    words >= 10 && 3 * distinct < words
}

/// Scores of the comparison: a match, a mismatch, a base inserted or deleted.
const MATCH: i32 = 2;
const MISMATCH: i32 = -3;
const GAP: i32 = -4;

/// The best-scoring alignment of a prefix of `q` with a prefix of `t`, up to `free` leading bases of one of them
/// skipped at no cost: (score, bases of q skipped, end in q, bases of t skipped, end in t). The alignment stays within
/// `band` of the diagonal; a row 30 below the best score ends it.
fn extend(q: &[u8], t: &[u8], free: usize, band: usize) -> (i32, usize, usize, usize, usize) {
    const NEG: i32 = i32::MIN / 4;
    let (m, n) = (q.len(), t.len());
    // origin: >= 0 bases of q skipped, < 0 bases of t skipped
    let mut prev: Vec<(i32, i32)> = (0..=n).map(|j| if j <= free { (0, -(j as i32)) } else { (NEG, 0) }).collect();
    let mut cur = vec![(NEG, 0i32); n + 1];
    let mut best = (0i32, 0usize, 0usize, 0i32);
    for i in 1..=m {
        let (lo, hi) = (i.saturating_sub(band).max(1), (i + band).min(n));
        if lo > hi { break; }
        cur[lo - 1] = if lo == 1 && i <= free { (0, i as i32) } else { (NEG, 0) };
        let mut row = NEG;
        for j in lo..=hi {
            let d = (prev[j - 1].0 + if q[i - 1] == t[j - 1] { MATCH } else { MISMATCH }, prev[j - 1].1);
            let u = (prev[j].0 + GAP, prev[j].1);
            let l = (cur[j - 1].0 + GAP, cur[j - 1].1);
            let x = if d.0 >= u.0 && d.0 >= l.0 { d } else if u.0 >= l.0 { u } else { l };
            cur[j] = x;
            row = row.max(x.0);
            if x.0 > best.0 { best = (x.0, i, j, x.1); }
        }
        if hi < n { cur[hi + 1] = (NEG, 0); }
        if row < best.0 - 30 { break; }
        std::mem::swap(&mut prev, &mut cur);
    }
    let (q0, t0) = if best.3 >= 0 { (best.3 as usize, 0) } else { (0, (-best.3) as usize) };
    (best.0, q0, best.1, t0, best.2)
}

/// Targets of one reading direction by the 10-mers of their first bases: each entry a target and the k-mer's position
/// on it.
struct KmerIndex {
    offsets: Vec<u32>,
    entries: Vec<u32>,
}

const POS_BITS: u32 = 6;

impl KmerIndex {
    fn build<'a>(targets: &'a [Target], forward: bool) -> Self {
        let mut counts = vec![0u32; (1usize << (2 * K)) + 1];
        let of = |t: &'a Target| kmers(&t.seq[..t.seq.len().min(SEED)]);
        for t in targets.iter().filter(|t| t.forward == forward) {
            for (_, k) in of(t) { counts[k as usize + 1] += 1; }
        }
        for i in 1..counts.len() { counts[i] += counts[i - 1]; }
        let mut entries = vec![0u32; *counts.last().unwrap() as usize];
        let mut fill = counts.clone();
        for (ti, t) in targets.iter().enumerate().filter(|(_, t)| t.forward == forward) {
            for (pos, k) in of(t) { entries[fill[k as usize] as usize] = ((ti as u32) << POS_BITS) | pos as u32; fill[k as usize] += 1; }
        }
        Self { offsets: counts, entries }
    }

    fn get(&self, k: u32) -> &[u32] {
        &self.entries[self.offsets[k as usize] as usize..self.offsets[k as usize + 1] as usize]
    }
}

/// Work counters, logged once per run: pieces looked at, searches at exon boundaries and inside exons, their seed
/// hits and extensions, own-locus alignments.
static N_PIECES: AtomicU64 = AtomicU64::new(0);
static N_PLACE: AtomicU64 = AtomicU64::new(0);
static N_LOCATE: AtomicU64 = AtomicU64::new(0);
static N_HITS: AtomicU64 = AtomicU64::new(0);
static N_EXTEND: AtomicU64 = AtomicU64::new(0);
static N_OWN: AtomicU64 = AtomicU64::new(0);
static N_ARTEFACT: AtomicU64 = AtomicU64::new(0);
/// Low-MAPQ parts named by their class.
static N_CLASS: AtomicU64 = AtomicU64::new(0);
/// Short parts placed as a class: classes, and pieces placed.
static N_SHORT_CLASS: AtomicU64 = AtomicU64::new(0);
static N_SHORT: AtomicU64 = AtomicU64::new(0);
/// Nanoseconds in: chain (all searches), own-locus alignments, seeding, extensions.
static T_CHAIN: AtomicU64 = AtomicU64::new(0);
static T_OWN: AtomicU64 = AtomicU64::new(0);
static T_SEED: AtomicU64 = AtomicU64::new(0);
static T_EXT: AtomicU64 = AtomicU64::new(0);
static N_FLANK: AtomicU64 = AtomicU64::new(0);

thread_local! { static DEBUG: std::cell::Cell<bool> = const { std::cell::Cell::new(false) }; }

/// Bases of the read before a placed part looked at for the part's own gene.
const BACK: usize = 80;
/// How far beyond an alignment's end a piece is looked for in the genome, as an exon of the aligned gene.
const OWN_REACH: i64 = 200_000;
const KF: usize = 12;
/// Bases of a piece, from its start, whose k-mers seed the search inside exons.
const SEED_FREE: usize = 60;

/// Every gene's annotated exons, merged, as genomic sequence, with the places of each 12-mer.
struct ExonStore {
    /// (gene, chromosome, low, high, offset in `seq`), by offset
    exons: Vec<(u32, u32, i64, i64, u32)>,
    seq: Vec<u8>,
    offsets: Vec<u32>,
    entries: Vec<u32>,
    /// gene -> plus strand
    plus: FxHashMap<u32, bool>,
    /// chromosome -> exon intervals with the exon's index
    trees: FxHashMap<String, COITree<u32, u32>>,
}

/// A piece entering an exon: its bases [q0, q1) on the genome's [lo, hi] of the exon, and whether they reach the
/// exon's end in the reading direction.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Located { exon: u32, score: i32, q0: usize, q1: usize, lo: i64, hi: i64, at_end: bool }

impl ExonStore {
    fn build(models: &[(u32, u32, bool, Vec<(i64, i64)>)], chroms: &[String], genome: &Path) -> Self {
        let mut per_gene: FxHashMap<u32, (u32, Vec<(i64, i64)>)> = FxHashMap::default();
        let mut plus: FxHashMap<u32, bool> = FxHashMap::default();
        for (gene, chrom, p, exons) in models {
            per_gene.entry(*gene).or_insert_with(|| (*chrom, Vec::new())).1.extend(exons.iter().copied());
            plus.insert(*gene, *p);
        }
        let mut merged: Vec<(u32, u32, i64, i64)> = Vec::new();
        let mut genes: Vec<u32> = per_gene.keys().copied().collect();
        genes.sort_unstable();
        for g in genes {
            let (chrom, mut ex) = per_gene.remove(&g).unwrap();
            ex.sort_unstable();
            let mut cur: Option<(i64, i64)> = None;
            for (a, b) in ex {
                match cur {
                    Some((lo, hi)) if a <= hi + 1 => cur = Some((lo, hi.max(b))),
                    Some((lo, hi)) => { merged.push((g, chrom, lo, hi)); cur = Some((a, b)); }
                    None => cur = Some((a, b)),
                }
            }
            if let Some((lo, hi)) = cur { merged.push((g, chrom, lo, hi)); }
        }
        let seqs: Vec<Vec<u8>> = merged.par_iter().map_init(
            || GenomeRef::open(genome.to_path_buf()).ok(),
            |g, &(_, c, lo, hi)| g.as_ref().and_then(|g| g.fetch_seq(&chroms[c as usize], lo - 1, hi).ok()).unwrap_or_default()).collect();
        let mut seq: Vec<u8> = Vec::with_capacity(seqs.iter().map(|x| x.len() + 1).sum());
        let mut exons = Vec::with_capacity(merged.len());
        let mut iv: FxHashMap<String, Vec<Interval<u32>>> = FxHashMap::default();
        for (&(g, c, lo, hi), x) in merged.iter().zip(&seqs) {
            if x.len() as i64 != hi - lo + 1 || seq.len() + x.len() + 1 >= u32::MAX as usize { continue; }
            iv.entry(chroms[c as usize].clone()).or_default().push(Interval::new(lo as i32, hi as i32, exons.len() as u32));
            exons.push((g, c, lo, hi, seq.len() as u32));
            seq.extend_from_slice(x);
            seq.push(b'N');
        }
        let mut counts = vec![0u32; (1usize << (2 * KF)) + 1];
        for (_, k) in kmers_of(&seq, KF) { counts[k as usize + 1] += 1; }
        for i in 1..counts.len() { counts[i] += counts[i - 1]; }
        let mut entries = vec![0u32; *counts.last().unwrap() as usize];
        let mut fill = counts.clone();
        for (pos, k) in kmers_of(&seq, KF) { entries[fill[k as usize] as usize] = pos as u32; fill[k as usize] += 1; }
        let trees = iv.into_iter().map(|(k, v)| (k, COITree::new(&v))).collect();
        Self { exons, seq, offsets: counts, entries, plus, trees }
    }

    fn get(&self, k: u32) -> &[u32] {
        &self.entries[self.offsets[k as usize] as usize..self.offsets[k as usize + 1] as usize]
    }

    /// The genes (with strand) whose exons hold a position.
    fn genes_at(&self, chrom: &str, pos: i64) -> Vec<(u32, bool)> {
        let mut out = Vec::new();
        if let Some(t) = self.trees.get(chrom) {
            t.query(pos as i32, pos as i32, |iv| {
                let g = self.exons[iv.metadata as usize].0;
                let x = (g, self.plus[&g]);
                if !out.contains(&x) { out.push(x); }
            });
        }
        out
    }

    /// Every exon a sequence (in the record's orientation) fits on either strand, scoring `min_score` at least: the
    /// seeding and extension of `locate` for a sequence read 3'-ward, every extension kept, no decision taken.
    fn candidates(&self, q: &[u8], o: &PlaceOptions, hits: &mut Vec<u64>) -> Vec<Cand> {
        const BIAS: i64 = 1 << 20;
        let n = q.len();
        let free = 2 * o.free;
        let mut out = Vec::new();
        for pp in [true, false] {
            let up = pp;
            let x: Vec<u8> = if pp { q.to_vec() } else { q.iter().rev().map(|&b| complement(b)).collect() };
            let m = n.min(SEED_FREE);
            let from = if up { 0 } else { n - m };
            hits.clear();
            for (i, k) in kmers_of(&x[from..from + m], KF) {
                let e = self.get(k);
                if e.len() > 400 { continue; }
                for &s in e { hits.push(((s as i64 - (from + i) as i64 + BIAS) as u64) << 8 | ((from + i) as u64 & 255)); }
            }
            hits.sort_unstable();
            let mut groups: Vec<(usize, i64, usize)> = Vec::new();
            let mut i = 0;
            while i < hits.len() {
                let mut k = i + 1;
                while k < hits.len() && (hits[k] >> 8) - (hits[k - 1] >> 8) <= 3 { k += 1; }
                let mut places: Vec<u64> = hits[i..k].iter().map(|h| h & 255).collect();
                places.sort_unstable(); places.dedup();
                if places.len() >= 2 {
                    let diag = (hits[i] >> 8) as i64 - BIAS;
                    let at = (diag + (hits[i] & 255) as i64).max(0) as u32;
                    let e = self.exons.partition_point(|x| x.4 <= at).saturating_sub(1);
                    let (gene, _, lo, hi, eo) = self.exons[e];
                    let (eo, elen) = (eo as i64, hi - lo + 1);
                    if self.plus[&gene] == pp && (at as i64) < eo + elen { groups.push((places.len(), diag, e)); }
                }
                i = k;
            }
            groups.sort_unstable_by(|a, b| b.0.cmp(&a.0).then(a.1.cmp(&b.1)));
            for &(_, diag, e) in groups.iter().take(24) {
                let (gene, chrom, lo, hi, eo) = self.exons[e];
                let (eo, elen) = (eo as i64, hi - lo + 1);
                let band = free + 8;
                let (t, top): (Vec<u8>, i64) = if up {
                    let a = (diag - o.free as i64).max(eo);
                    let b = (diag + (n + band) as i64).min(eo + elen);
                    if b <= a { continue; }
                    (self.seq[a as usize..b as usize].to_vec(), a)
                } else {
                    let b = (diag + n as i64 - 1 + o.free as i64).min(eo + elen - 1);
                    let a = (diag - band as i64).max(eo);
                    if b < a { continue; }
                    (self.seq[a as usize..=b as usize].iter().rev().map(|&c| complement(c)).collect(), b)
                };
                let (score, q0, q1, t0, t1) = extend(q, &t, free, band);
                if score >= o.min_score && t1 > t0 {
                    let (s0, s1) = if up { (top + t0 as i64, top + t1 as i64 - 1) } else { (top - t1 as i64 + 1, top - t0 as i64) };
                    out.push(Cand { gene, chrom, lo: lo + (s0 - eo), hi: lo + (s1 - eo), score, q0, q1, plus: pp });
                }
            }
        }
        out
    }

    /// Where a piece (read outward, its first base next to the alignment) enters an exon: the piece is the
    /// transcript read 3'-ward (`forward`) or 5'-ward, so on a plus-strand gene it runs up the genome when read
    /// 3'-ward, on a minus-strand gene when read 5'-ward. Whether the exon is of a `local` gene; or why none.
    fn locate(&self, q: &[u8], forward: bool, local: &[u32], only: Option<u32>, o: &PlaceOptions, hits: &mut Vec<u64>) -> Result<(Located, bool), &'static str> {
        const BIAS: i64 = 1 << 20;
        let n = q.len();
        let free = 2 * o.free;
        N_LOCATE.fetch_add(1, Relaxed);
        let mut found: Vec<Located> = Vec::new();
        let (mut words, mut common) = (0u32, 0u32);
        for pp in [true, false] {
            let up = forward == pp;
            // the piece as the genome's plus strand reads it
            let x: Vec<u8> = match (up, pp) {
                (true, true) => q.to_vec(),
                (true, false) => q.iter().map(|&b| complement(b)).collect(),
                (false, true) => q.iter().rev().copied().collect(),
                (false, false) => q.iter().rev().map(|&b| complement(b)).collect(),
            };
            let m = n.min(SEED_FREE);
            let from = if up { 0 } else { n - m };
            let t_seed = std::time::Instant::now();
            hits.clear();
            for (i, k) in kmers_of(&x[from..from + m], KF) {
                let e = self.get(k);
                words += 1;
                // a k-mer at more than 400 places is common: no seed
                if e.len() > 400 { common += 1; continue; }
                for &s in e { hits.push(((s as i64 - (from + i) as i64 + BIAS) as u64) << 8 | ((from + i) as u64 & 255)); }
            }
            hits.sort_unstable();
            N_HITS.fetch_add(hits.len() as u64, Relaxed);
            T_SEED.fetch_add(t_seed.elapsed().as_nanos() as u64, Relaxed);
            if DEBUG.with(|d| d.get()) { log::info!("[Locate debug] partner strand {} up {} x {} hits {}", if pp { '+' } else { '-' }, up, String::from_utf8_lossy(&x[..x.len().min(40)]), hits.len()); }
            // seed groups: hits on one diagonal (within 3 of each other), at two places of the piece or more, in an exon
            // of the strand, whose flanks agree; the groups with most seeds are extended
            let mut groups: Vec<(usize, i64, usize)> = Vec::new();
            let mut i = 0;
            while i < hits.len() {
                let mut k = i + 1;
                while k < hits.len() && (hits[k] >> 8) - (hits[k - 1] >> 8) <= 3 { k += 1; }
                let mut places: Vec<u64> = hits[i..k].iter().map(|h| h & 255).collect();
                places.sort_unstable(); places.dedup();
                if places.len() >= 2 {
                    let diag = (hits[i] >> 8) as i64 - BIAS;
                    let at = (diag + (hits[i] & 255) as i64).max(0) as u32;
                    let e = self.exons.partition_point(|x| x.4 <= at).saturating_sub(1);
                    let (gene, _, lo, hi, eo) = self.exons[e];
                    let (eo, elen) = (eo as i64, hi - lo + 1);
                    if self.plus[&gene] == pp && only.map_or(true, |g| g == gene) && (at as i64) < eo + elen
                        && flanks_agree((hits[i] & 255) as usize, KF, n, |j, shift| {
                            let sj = diag + j as i64 + shift;
                            (sj >= eo && sj < eo + elen).then(|| (x[j], self.seq[sj as usize]))
                        }) { groups.push((places.len(), diag, e)); }
                }
                i = k;
            }
            groups.sort_unstable_by(|a, b| b.0.cmp(&a.0).then(a.1.cmp(&b.1)));
            for &(_, diag, e) in groups.iter().take(12) {
                let (_, _, lo, hi, eo) = self.exons[e];
                let (eo, elen) = (eo as i64, hi - lo + 1);
                let band = free + 8;
                let (t, top): (Vec<u8>, i64) = if up {
                    let a = (diag - o.free as i64).max(eo);
                    let b = (diag + (n + band) as i64).min(eo + elen);
                    if b <= a { continue; }
                    (self.seq[a as usize..b as usize].iter().map(|&c| if pp { c } else { complement(c) }).collect(), a)
                } else {
                    let b = (diag + n as i64 - 1 + o.free as i64).min(eo + elen - 1);
                    let a = (diag - band as i64).max(eo);
                    if b < a { continue; }
                    (self.seq[a as usize..=b as usize].iter().rev().map(|&c| if pp { c } else { complement(c) }).collect(), b)
                };
                N_EXTEND.fetch_add(1, Relaxed);
                let t_ext = std::time::Instant::now();
                let (score, q0, q1, t0, t1) = extend(q, &t, free, band);
                T_EXT.fetch_add(t_ext.elapsed().as_nanos() as u64, Relaxed);
                if DEBUG.with(|d| d.get()) { log::info!("[Locate debug] exon {}:{}-{} diag {} t {} score {} q {}-{} t {}-{}", self.exons[e].1, lo, hi, diag, String::from_utf8_lossy(&t[..t.len().min(50)]), score, q0, q1, t0, t1); }
                if score >= o.min_score_free.min(o.min_score) - 2 * o.margin && t1 > t0 {
                    // the target's bases [t0, t1) from `top`, up or down the exon's sequence
                    let (s0, s1) = if up { (top + t0 as i64, top + t1 as i64 - 1) } else { (top - t1 as i64 + 1, top - t0 as i64) };
                    let at_end = if up { eo + elen - 1 - s1 <= 3 } else { s0 - eo <= 3 };
                    found.push(Located { exon: e as u32, score, q0, q1, lo: lo + (s0 - eo), hi: lo + (s1 - eo), at_end });
                }
            }
        }
        // a piece a fifth of whose k-mers (on either strand) are common is a repeat
        if only.is_none() && common >= 4 && 5 * common >= words { return Err("repeat"); }
        let Some(&b) = found.iter().max_by_key(|x| x.score) else { return Err("no exon") };
        let placed = b.q1 - b.q0;
        // within one gene, the score of a part at an exon boundary is enough
        if only.is_some() { return if b.score >= o.min_score && placed >= o.min_len { Ok((b, false)) } else { Err("weak") }; }
        if b.score < o.min_score_free || placed < o.min_len || (b.score as f64) < o.min_score_per_base * placed as f64 { return Err("weak"); }
        let (gene, chrom, _, _, _) = self.exons[b.exon as usize];
        // a piece that fits three other places six tenths as well is a repeat
        let mut alike: Vec<(u32, i64)> = found.iter().filter(|x| 10 * x.score >= 6 * b.score)
            .filter(|x| { let c = self.exons[x.exon as usize].1; c != chrom || x.hi < b.lo || b.hi < x.lo }).map(|x| (self.exons[x.exon as usize].1, x.lo)).collect();
        alike.sort_unstable(); alike.dedup();
        if alike.len() >= 3 { return Err("repeat"); }
        // a rival: another place scoring within `margin` (the same place under an overlapping gene is none)
        let elsewhere = |x: &&Located| { let c = self.exons[x.exon as usize].1; x.score >= b.score - o.margin && (c != chrom || x.hi < b.lo || b.hi < x.lo) };
        if found.iter().filter(elsewhere).any(|x| self.exons[x.exon as usize].0 != gene) { return Err("two genes"); }
        let own = local.contains(&gene) || found.iter().any(|x| x.score >= b.score - 2 * o.margin && local.contains(&self.exons[x.exon as usize].0));
        Ok((b, own))
    }
}

/// One placed part of a piece: the piece's bases [off + q0, off + q1) on `blocks` (ascending) of a gene, whether the
/// gene is the one the piece left (its own next exon), and whether the part ends at an exon's end.
#[derive(Debug, Clone)]
struct Step { off: usize, q0: usize, q1: usize, score: i32, gene: u32, chrom: u32, plus: bool, blocks: Vec<(i64, i64)>, own: bool, at_end: bool, anchored: bool }

/// A piece placed on a target: its bases [q0, q1) on the target's bases [t0, t1).
#[derive(Debug, Clone, Copy, PartialEq)]
struct Placed { target: u32, score: i32, q0: usize, q1: usize, t0: usize, t1: usize }

struct ExonIndex {
    chroms: Vec<String>,
    targets: Vec<Target>,
    forward: KmerIndex,
    backward: KmerIndex,
    /// chromosome -> position -> (gene, plus strand, exon start in transcript orientation) of the exons bounded there
    bounds: FxHashMap<String, FxHashMap<i64, Vec<(u32, bool, bool)>>>,
    store: ExonStore,
}

impl ExonIndex {
    /// The target a piece (read outward from the boundary) continues into, and whether it is an exon of a `local`
    /// gene (the aligned gene's own next exon); or why none: no target scores enough, or two genes' targets at
    /// different places score alike.
    fn place(&self, q: &[u8], forward: bool, local: &[u32], only: Option<u32>, o: &PlaceOptions, hits: &mut Vec<u64>) -> Result<(Placed, bool), (&'static str, Option<Placed>)> {
        let index = if forward { &self.forward } else { &self.backward };
        // seeds: k-mers of the piece's first bases at about the same place on the target
        let drift = o.free as i64 + 4;
        N_PLACE.fetch_add(1, Relaxed);
        let t_seed = std::time::Instant::now();
        hits.clear();
        let (mut words, mut common) = (0u32, 0u32);
        for (qpos, k) in kmers(&q[..q.len().min(SEED + o.free)]) {
            let bucket = index.get(k);
            words += 1;
            // a k-mer of more than 1500 targets is common: no seed
            if bucket.len() > 1500 { common += 1; continue; }
            for &e in bucket {
                let diag = qpos as i64 - (e & ((1 << POS_BITS) - 1)) as i64;
                if diag.abs() <= drift { hits.push(((e >> POS_BITS) as u64) << 16 | ((diag + 128) as u64) << 8 | qpos as u64); }
            }
        }
        hits.sort_unstable();
        N_HITS.fetch_add(hits.len() as u64, Relaxed);
        T_SEED.fetch_add(t_seed.elapsed().as_nanos() as u64, Relaxed);
        // a piece a fifth of whose k-mers are common is a repeat
        if only.is_none() && common >= 4 && 5 * common >= words { return Err(("repeat", None)); }
        let mut best: Option<Placed> = None;
        let mut rivals: Vec<(u32, i32)> = Vec::new();
        // targets with two seeds on one diagonal (within 3), at different places of the piece, whose flanks agree; the
        // ones with most seeds are extended
        let mut seeded: Vec<(usize, u32)> = Vec::new();
        let mut i = 0;
        while i < hits.len() {
            let target = (hits[i] >> 16) as u32;
            let mut k = i;
            while k < hits.len() && (hits[k] >> 16) as u32 == target { k += 1; }
            let t = &self.targets[target as usize];
            if only.map_or(true, |g| g == t.gene)
                && hits[i..k].windows(2).filter(|w| ((w[1] >> 8) & 255) - ((w[0] >> 8) & 255) <= 3 && (w[1] & 255) != (w[0] & 255)).any(|w| {
                    let (diag, qpos) = (((w[0] >> 8) & 255) as i64 - 128, (w[0] & 255) as usize);
                    flanks_agree(qpos, K, q.len(), |j, shift| { let tj = j as i64 - diag + shift; (tj >= 0 && (tj as usize) < t.seq.len()).then(|| (q[j], t.seq[tj as usize])) })
                }) { seeded.push((k - i, target)); }
            i = k;
        }
        seeded.sort_unstable_by(|a, b| b.0.cmp(&a.0).then(a.1.cmp(&b.1)));
        for &(_, target) in seeded.iter().take(12) {
            let t = &self.targets[target as usize];
            N_EXTEND.fetch_add(1, Relaxed);
            let t_ext = std::time::Instant::now();
            let (score, q0, q1, t0, t1) = extend(q, &t.seq, o.free, o.free + 8);
            T_EXT.fetch_add(t_ext.elapsed().as_nanos() as u64, Relaxed);
            if score >= o.min_score - 2 * o.margin {
                rivals.push((target, score));
                if best.map_or(true, |b| score > b.score) { best = Some(Placed { target, score, q0, q1, t0, t1 }); }
            }
        }
        let Some(b) = best else { return Err(("no exon", None)) };
        let placed = b.q1 - b.q0;
        if b.score < o.min_score || placed < o.min_len || (b.score as f64) < o.min_score_per_base * placed as f64 { return Err(("weak", Some(b))); }
        let t = &self.targets[b.target as usize];
        if only.is_some() { return Ok((b, false)); }
        // a piece that fits three other places six tenths as well is a repeat
        let apart = |x: &Target| x.chrom != t.chrom || x.blocks[0].1 < t.blocks[0].0 || t.blocks[0].1 < x.blocks[0].0;
        let mut alike: Vec<(u32, i64)> = rivals.iter().filter(|r| 10 * r.1 >= 6 * b.score).map(|r| &self.targets[r.0 as usize]).filter(|x| apart(x)).map(|x| (x.chrom, x.blocks[0].0)).collect();
        alike.sort_unstable(); alike.dedup();
        if alike.len() >= 3 { return Err(("repeat", Some(b))); }
        // a rival: a target scoring within `margin` of the best one, of another gene and at another place
        let same_place = |x: &Target| x.chrom == t.chrom && x.blocks[0].0 <= t.blocks[0].1 && t.blocks[0].0 <= x.blocks[0].1;
        let close = |&&(x, sc): &&(u32, i32)| sc >= b.score - o.margin && x != b.target;
        // the aligned gene's own exon is taken when it scores within twice the margin of another gene's
        if let Some(&(x, sc)) = rivals.iter().filter(|r| r.1 >= (b.score - 2 * o.margin).max(o.min_score) && local.contains(&self.targets[r.0 as usize].gene)).max_by_key(|r| r.1) {
            let own = if x == b.target { b } else {
                let t = &self.targets[x as usize];
                let (score, q0, q1, t0, t1) = extend(q, &t.seq, o.free, o.free + 8);
                debug_assert_eq!(score, sc);
                Placed { target: x, score, q0, q1, t0, t1 }
            };
            return Ok((own, true));
        }
        if rivals.iter().filter(close).any(|&(x, _)| { let x = &self.targets[x as usize]; x.gene != t.gene && !same_place(x) }) {
            return Err(("two genes", Some(b)));
        }
        Ok((b, false))
    }

    /// Whether the read's bases before a part placed on `gene` are that gene's too: read back from the part's start
    /// (`before`, the bases next to it first), they continue into the gene's transcript the other way, at an exon
    /// boundary or inside the exon.
    fn continues_back(&self, before: &[u8], forward: bool, gene: u32, o: &PlaceOptions, hits: &mut Vec<u64>) -> bool {
        before.len() >= o.min_len
            && (self.place(before, !forward, &[], Some(gene), o, hits).is_ok() || self.store.locate(before, !forward, &[], Some(gene), o, hits).is_ok())
    }

    /// Whether the genome beside a placed part continues the read backwards: the read's bases before the part
    /// (`before`, the bases next to it first) against the genome's bases before the part's first, on its strand.
    fn genome_continues_back(&self, g: &GenomeRef, before: &[u8], forward: bool, st: &Step, o: &PlaceOptions) -> bool {
        if before.len() < o.min_len || st.blocks.is_empty() { return false; }
        let up = forward == st.plus;
        let n = before.len() as i64 + 2 * o.free as i64 + 8;
        let chrom = &self.chroms[st.chrom as usize];
        // read backwards, the bases before a part read up the genome run down from its first base, and the other way
        let t: Vec<u8> = if up {
            let first = st.blocks[0].0;
            let Ok(mut t) = g.fetch_seq(chrom, (first - 1 - n).max(0), first - 1) else { return false };
            t.reverse();
            t
        } else {
            let last = st.blocks[st.blocks.len() - 1].1;
            let Ok(t) = g.fetch_seq(chrom, last, last + n) else { return false };
            t
        };
        let t: Vec<u8> = if st.plus { t } else { t.iter().map(|&b| complement(b)).collect() };
        let (score, q0, q1, _, _) = extend(before, &t, 2 * o.free, 2 * o.free + 8);
        if DEBUG.with(|d| d.get()) { log::info!("[Locate debug] genome beside the part on {}: the read's bases before it {} against {} score {} over {}", chrom, String::from_utf8_lossy(&before[..before.len().min(40)]), String::from_utf8_lossy(&t[..t.len().min(40)]), score, q1 - q0); }
        score >= o.min_score && q1 - q0 >= o.min_len
    }

    /// Whether a piece (in the reference's orientation) lies in the genome within `reach` bp beyond the alignment's end
    /// (`high` or low), outside the annotated exons of genes of its strand (`plus`) other than `local`: its first 60
    /// bases score 28 and within twice `margin` of its `score` on the other gene, or any part 38 and within `margin`.
    fn in_own_locus(&self, g: &GenomeRef, chrom: &str, boundary: i64, high: bool, plus: bool, piece: &[u8], local: &[u32], score: i32, reach: i64, o: &PlaceOptions) -> bool {
        let n = piece.len().min(60);
        if n < o.min_len { return false; }
        N_OWN.fetch_add(1, Relaxed);
        let t_own = std::time::Instant::now();
        // the piece's bases next to the boundary, read outward from it
        let q: Vec<u8> = if high { piece[..n].to_vec() } else { piece[piece.len() - n..].iter().rev().copied().collect() };
        let (lo, hi) = if high { (boundary + 1, boundary + reach) } else { ((boundary - reach).max(1), boundary - 1) };
        let Ok(mut win) = g.fetch_seq(chrom, lo - 1, hi) else { return false };
        // annotated exons of other genes of the same strand there are masked
        if let Some(t) = self.store.trees.get(chrom) {
            t.query(lo as i32, hi as i32, |iv| {
                let gene = self.store.exons[iv.metadata as usize].0;
                if !local.contains(&gene) && self.store.plus[&gene] == plus {
                    let (a, b) = ((iv.first as i64).max(lo) - lo, (iv.last as i64).min(hi) - lo);
                    let last = (b as usize).min(win.len() - 1);
                    for x in win[a as usize..=last].iter_mut() { *x = b'N'; }
                }
            });
        }
        // read outward too: down the genome beyond the low end
        if !high { win.reverse(); }
        // places sharing two 6-mers with the piece on one diagonal, each aligned without further seeds: the piece from
        // its start (and its rest from where that ends), or any part of it
        let mut at: Vec<Vec<u16>> = vec![Vec::new(); 1 << 12];
        // (a 6-mer of one or two distinct bases seeds nothing)
        for (i, k) in kmers_of(&q, 6) {
            let w = &q[i..i + 6];
            if w.iter().filter(|&&b| b != w[0]).map(|&b| b).collect::<FxHashSet<u8>>().len() >= 2 { at[k as usize].push(i as u16); }
        }
        let mut diags: Vec<i64> = Vec::new();
        for (p, k) in kmers_of(&win, 6) { for &i in &at[k as usize] { diags.push(p as i64 - i as i64); } }
        diags.sort_unstable();
        // the diagonals with most seeds (40 at most)
        let mut groups: Vec<(usize, i64, i64)> = Vec::new();
        let mut i = 0;
        while i < diags.len() {
            let mut k = i + 1;
            while k < diags.len() && diags[k] - diags[k - 1] <= 3 { k += 1; }
            if k - i >= 2 { groups.push((k - i, diags[i], diags[k - 1])); }
            i = k;
        }
        groups.sort_unstable_by(|a, b| b.0.cmp(&a.0).then(a.1.cmp(&b.1)));
        groups.truncate(40);
        let (mut own, mut end, mut local_best) = (0, 0, 0);
        let mut two_parts = false;
        for &(_, d0, d1) in &groups {
            let a = (d0 - 8).max(0) as usize;
            let b = ((d1 + n as i64 + 8).max(0) as usize).min(win.len());
            if a < b {
                let (x, e) = locus_score(&q, &win[a..b], o.free);
                if x > own { own = x; end = e; }
                local_best = local_best.max(locus_score(&q, &win[a..b], n).0);
            }
        }
        if own > 0 && end + 15 <= n {
            // the rest of the piece, anywhere in the stretch by the same seeds
            let rest = &q[end..];
            let mut best_rest = 0;
            for &(_, d0, d1) in &groups {
                let a = (d0 + end as i64 - 8).max(0) as usize;
                let b = ((d1 + n as i64 + 8).max(0) as usize).min(win.len());
                if a < b { best_rest = best_rest.max(locus_score(rest, &win[a..b], o.free).0); }
            }
            // two parts count only when the first scores 24 and the rest 28
            two_parts = own >= 24 && best_rest >= 28;
            if two_parts { own += best_rest; }
        }
        let local = local_best;
        if DEBUG.with(|d| d.get()) { log::info!("[Locate debug] own locus {}:{}-{} piece {} scores {} (first part to {}), any part {}, the other gene {}", chrom, lo, hi, String::from_utf8_lossy(&q), own, end, local, score); }
        T_OWN.fetch_add(t_own.elapsed().as_nanos() as u64, Relaxed);
        let _ = two_parts;
        own >= (score - 2 * o.margin).max(28) || local >= (score - o.margin).max(38)
    }

    /// The best pooled score of a short class's pieces (reference orientation, read outward) at one place of the genome
    /// within `reach` bp beyond its alignment end, outside the annotated exons of other genes of the strand `plus` and
    /// of `partner`. The places tried are the 16 where most pieces share an 8-mer of their first `SHORT_MIN` bases.
    fn short_genome_score(&self, g: &GenomeRef, chrom: &str, boundary: i64, high: bool, plus: bool, pieces: &[Vec<u8>], local: &[u32], partner: u32, reach: i64, free: usize) -> i32 {
        let (lo, hi) = if high { (boundary + 1, boundary + reach) } else { ((boundary - reach).max(1), boundary - 1) };
        if hi < lo { return 0; }
        let Ok(mut win) = g.fetch_seq(chrom, lo - 1, hi) else { return 0 };
        if win.is_empty() { return 0; }
        if let Some(t) = self.store.trees.get(chrom) {
            t.query(lo as i32, hi as i32, |iv| {
                let gene = self.store.exons[iv.metadata as usize].0;
                if !local.contains(&gene) && (self.store.plus[&gene] == plus || gene == partner) {
                    let (a, b) = ((iv.first as i64).max(lo) - lo, (iv.last as i64).min(hi) - lo);
                    if a <= b && (a as usize) < win.len() {
                        let last = (b as usize).min(win.len() - 1);
                        for x in win[a as usize..=last].iter_mut() { *x = b'N'; }
                    }
                }
            });
        }
        for x in win.iter_mut() { *x = x.to_ascii_uppercase(); }
        // read outward: down the genome beyond the low end
        if !high { win.reverse(); }
        let mut words: FxHashMap<u32, Vec<(u16, u16)>> = FxHashMap::default();
        for (j, q) in pieces.iter().enumerate() {
            for (i, k) in kmers_of(&q[..q.len().min(SHORT_MIN)], 8) { words.entry(k).or_default().push((j as u16, i as u16)); }
        }
        let mut at: FxHashMap<i64, FxHashSet<u16>> = FxHashMap::default();
        for (p, k) in kmers_of(&win, 8) {
            if let Some(v) = words.get(&k) { for &(j, i) in v { at.entry(p as i64 - i as i64).or_default().insert(j); } }
        }
        let mut places: Vec<(usize, i64)> = at.into_iter().map(|(a, s)| (s.len(), a)).collect();
        places.sort_unstable_by(|x, y| y.0.cmp(&x.0).then(x.1.cmp(&y.1)));
        places.truncate(16);
        places.iter().map(|&(_, a)| {
            let a0 = a.max(0) as usize;
            let s0 = a0.saturating_sub(free);
            pieces.iter().map(|q| {
                let e = (a0 + q.len() + free + 8).min(win.len());
                if s0 >= e { 0 } else { extend(q, &win[s0..e], 2 * free, 2 * free + 8).0 }
            }).sum::<i32>()
        }).max().unwrap_or(0)
    }

    /// A piece placed step by step: at an exon boundary (`at_boundary` for the first step, then wherever the step
    /// before ended at its exon's end) on the exons that start there, else, or when none fits, where it enters an
    /// exon; the genes to tell from being the gene of the step before. A step on the gene the piece left that is not
    /// its next exon from the boundary ends the chain.
    fn chain(&self, q: &[u8], forward: bool, local: &[u32], at_boundary: bool, o: &PlaceOptions, hits: &mut Vec<u64>) -> Vec<Step> {
        let mut steps: Vec<Step> = Vec::new();
        let mut off = 0usize;
        let mut local = local.to_vec();
        let mut anchored = at_boundary;
        while steps.len() < 4 && q.len() >= off + o.min_len {
            let sub = &q[off..q.len().min(off + o.max_len)];
            let mut step: Option<Step> = None;
            if anchored {
                if let Ok((p, own)) = self.place(sub, forward, &local, None, o, hits) {
                    let t = &self.targets[p.target as usize];
                    let mut at = 0i64;
                    let mut at_end = false;
                    for &(lo, hi) in &t.blocks {
                        at += hi - lo + 1;
                        if (at - p.t1 as i64).abs() <= 3 { at_end = true; }
                    }
                    step = Some(Step { off, q0: p.q0, q1: p.q1, score: p.score, gene: t.gene, chrom: t.chrom, plus: t.plus, blocks: t.cover(p.t0, p.t1), own, at_end, anchored: true });
                }
            }
            if step.is_none() {
                if let Ok((l, own)) = self.store.locate(sub, forward, &local, None, o, hits) {
                    if own { break; }
                    let (gene, chrom, _, _, _) = self.store.exons[l.exon as usize];
                    step = Some(Step { off, q0: l.q0, q1: l.q1, score: l.score, gene, chrom, plus: self.store.plus[&gene], blocks: vec![(l.lo, l.hi)], own: false, at_end: l.at_end, anchored: false });
                }
            }
            let Some(st) = step else { break };
            if low_complexity(&sub[st.q0..st.q1]) { break; }
            off += st.q1;
            local = vec![st.gene];
            anchored = st.at_end;
            let blocks_empty = st.blocks.is_empty();
            steps.push(st);
            if blocks_empty { break; }
        }
        steps
    }
}

/// Annotated transcripts' exons in transcript order as genomic (low, high), with gene (index in `db.genes`).
fn annotated_models(db: &PostDB) -> (Vec<String>, Vec<(u32, u32, bool, Vec<(i64, i64)>)>) {
    let gidx: FxHashMap<i64, u32> = db.genes.iter().enumerate().map(|(i, g)| (g.gene_id, i as u32)).collect();
    let edge_ends: FxHashMap<i64, (i64, i64)> = db.edges.iter().map(|e| (e.edge_id, (e.v1, e.v2))).collect();
    let mut chroms: Vec<String> = Vec::new();
    let mut chrom_id: FxHashMap<String, u32> = FxHashMap::default();
    let mut out = Vec::new();
    for t in db.transcripts.iter().filter(|t| db.annot.known_tx_ids.contains(&t.transcript_id)) {
        let Some(&g) = gidx.get(&t.gene_id) else { continue };
        let Some(pos) = crate::quant_em::model_positions(db, t, &edge_ends) else { continue };
        if pos.len() < 2 { continue; }
        let plus = db.genes[g as usize].strand == "+";
        let c = *chrom_id.entry(t.chromosome.clone()).or_insert_with(|| { chroms.push(t.chromosome.clone()); chroms.len() as u32 - 1 });
        out.push((g, c, plus, pos.chunks(2).map(|p| (p[0].min(p[1]), p[0].max(p[1]))).collect()));
    }
    (chroms, out)
}

/// Every exon boundary of the models read outward for `len` bases: target blocks without sequence, each once.
fn target_specs(models: &[(u32, u32, bool, Vec<(i64, i64)>)], len: usize) -> Vec<Target> {
    let mut seen: FxHashMap<(u32, bool, bool, Vec<(i64, i64)>), usize> = FxHashMap::default();
    let mut out: Vec<Target> = Vec::new();
    for (gene, chrom, plus, exons) in models {
        for k in 0..exons.len() {
            for forward in [true, false] {
                // the model's own first exon read 3'-ward, or its last read 5'-ward: a transcript end
                let tx_end = (forward && k == 0) || (!forward && k + 1 == exons.len());
                let mut blocks = Vec::new();
                let mut rem = len as i64;
                let mut j = k as i64;
                while rem > 0 && j >= 0 && (j as usize) < exons.len() {
                    let (lo, hi) = exons[j as usize];
                    let take = rem.min(hi - lo + 1);
                    // reading up the genome starts at the exon's low end
                    blocks.push(if forward == *plus { (lo, lo + take - 1) } else { (hi - take + 1, hi) });
                    rem -= take;
                    j += if forward { 1 } else { -1 };
                }
                match seen.get(&(*chrom, *plus, forward, blocks.clone())) {
                    Some(&i) => { out[i].tx_end |= tx_end; }
                    None => {
                        seen.insert((*chrom, *plus, forward, blocks.clone()), out.len());
                        out.push(Target { gene: *gene, chrom: *chrom, plus: *plus, forward, blocks, seq: Vec::new(), tx_end });
                    }
                }
            }
        }
    }
    out
}

fn build_index(db: &PostDB, genome: &Path, o: &PlaceOptions) -> ExonIndex {
    let (chroms, models) = annotated_models(db);
    let mut targets = target_specs(&models, o.max_len + 2 * o.free + 8);
    let seqs: Vec<Vec<u8>> = targets.par_iter().map_init(
        || GenomeRef::open(genome.to_path_buf()).ok(),
        |g, t| {
            let Some(g) = g.as_ref() else { return Vec::new() };
            let mut s = Vec::new();
            for &(lo, hi) in &t.blocks {
                let Ok(mut b) = g.fetch_seq(&chroms[t.chrom as usize], lo - 1, hi) else { return Vec::new() };
                // the block in reading order: reversed when read down the genome, complemented for a minus-strand gene
                if !t.ascending() { b.reverse(); }
                if !t.plus { for x in b.iter_mut() { *x = complement(*x); } }
                s.extend(b);
            }
            s
        }).collect();
    for (t, s) in targets.iter_mut().zip(seqs) { t.seq = s; }
    let mut bounds: FxHashMap<String, FxHashMap<i64, Vec<(u32, bool, bool)>>> = FxHashMap::default();
    for (gene, chrom, plus, exons) in &models {
        let m = bounds.entry(chroms[*chrom as usize].clone()).or_default();
        for &(lo, hi) in exons {
            let (start, end) = if *plus { (lo, hi) } else { (hi, lo) };
            for (p, is_start) in [(start, true), (end, false)] {
                let v = m.entry(p).or_default();
                if !v.contains(&(*gene, *plus, is_start)) { v.push((*gene, *plus, is_start)); }
            }
        }
    }
    let (forward, backward) = (KmerIndex::build(&targets, true), KmerIndex::build(&targets, false));
    let store = ExonStore::build(&models, &chroms, genome);
    ExonIndex { chroms, targets, forward, backward, bounds, store }
}

/// The piece beyond an alignment's end, read outward from it along the transcript of a gene on `plus`: beyond the
/// alignment's high end (`high`) of a plus gene, or its low end of a minus gene, the piece runs 3'-ward (returned
/// true) and is compared with exon starts; otherwise 5'-ward, with exon ends. `piece` is in the reference's
/// orientation.
fn outward(piece: &[u8], high: bool, plus: bool) -> (Vec<u8>, bool) {
    match (high, plus) {
        (true, true) => (piece.to_vec(), true),
        (true, false) => (piece.iter().map(|&b| complement(b)).collect(), false),
        (false, true) => (piece.iter().rev().copied().collect(), false),
        (false, false) => (piece.iter().rev().map(|&b| complement(b)).collect(), true),
    }
}

/// A read's alignment extended over its gene's own exons beyond one end: the exons' blocks (ascending) join the
/// alignment across the intron from the boundary, and the read's bases placed on them join its part of the read.
/// False (nothing changed) when the blocks do not lie beyond the boundary.
fn extend_segment(base: &mut ChimericSegment, high: bool, boundary: i64, blocks: &[(i64, i64)], query: (i32, i32)) -> bool {
    let (Some(&(lo, _)), Some(&(_, hi))) = (blocks.first(), blocks.last()) else { return false };
    let mut add = Vec::with_capacity(2 * blocks.len());
    if high {
        if lo <= boundary + 1 { return false; }
        add.push(boundary + 1); add.push(lo - 1);
        for x in blocks.windows(2) { add.push(x[0].1 + 1); add.push(x[1].0 - 1); }
        base.introns.extend(add);
        base.ref_end = hi;
    } else {
        if hi >= boundary - 1 { return false; }
        for x in blocks.windows(2) { add.push(x[0].1 + 1); add.push(x[1].0 - 1); }
        add.push(hi + 1); add.push(boundary - 1);
        add.extend(base.introns.iter().copied());
        base.introns = add;
        base.ref_start = lo;
    }
    base.query_start = base.query_start.min(query.0);
    base.query_end = base.query_end.max(query.1);
    true
}

/// Each aligned base of a record: its position on the read as sequenced, its reference position (1-based; -1 for an
/// inserted base) and whether it matches the genome.
fn aligned_bases(rec: &bam::Record, genome: &GenomeRef, chrom: &str) -> Vec<(i32, i64, bool)> {
    let seq = rec.seq().as_bytes();
    let cig = rec.cigar();
    let lead_h = match cig.first() { Some(Cigar::HardClip(n)) => *n as usize, _ => 0 };
    let tail_h = if cig.len() > 1 { match cig.last() { Some(Cigar::HardClip(n)) => *n as usize, _ => 0 } } else { 0 };
    let full = lead_h + seq.len() + tail_h;
    let rev = rec.is_reverse();
    let to_seq = |q: usize| -> i32 { let p = lead_h + q; if rev { (full - 1 - p) as i32 } else { p as i32 } };
    let mut out = Vec::with_capacity(seq.len());
    let (mut q, mut r) = (0usize, rec.pos() + 1);
    for op in cig.iter() {
        match *op {
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                let n = n as usize;
                let g = genome.fetch_seq(chrom, r - 1, r - 1 + n as i64).unwrap_or_default();
                for k in 0..n {
                    let m = g.get(k).is_some_and(|&b| b == seq[q + k].to_ascii_uppercase());
                    out.push((to_seq(q + k), r + k as i64, m));
                }
                q += n; r += n as i64;
            }
            Cigar::Ins(n) => { for k in 0..n as usize { out.push((to_seq(q + k), -1, false)); } q += n as usize; }
            Cigar::SoftClip(n) => { q += n as usize; }
            Cigar::Del(n) | Cigar::RefSkip(n) => { r += n as i64; }
            Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }
    out
}

/// The record of a read's supplementary alignment starting at `pos` (1-based) on `chrom`.
fn fetch_supplementary(f: &mut bam::IndexedReader, chrom: &str, pos: i64, name: &[u8]) -> Option<bam::Record> {
    f.fetch((chrom, pos - 1, pos)).ok()?;
    let mut r = bam::Record::new();
    while let Some(Ok(())) = f.read(&mut r) {
        if r.is_supplementary() && !r.is_secondary() && r.pos() + 1 == pos && r.qname() == name { return Some(r); }
    }
    None
}

/// A segment cut down to its aligned bases at the read positions `keep` holds; false when fewer than `min` are left.
fn trim_segment(seg: &mut ChimericSegment, bases: &[(i32, i64, bool)], keep: impl Fn(i32) -> bool, min: usize) -> bool {
    let (mut lo_q, mut hi_q, mut lo_r, mut hi_r, mut n) = (i32::MAX, i32::MIN, i64::MAX, i64::MIN, 0usize);
    for &(q, r, _) in bases {
        if !keep(q) { continue; }
        n += 1; lo_q = lo_q.min(q); hi_q = hi_q.max(q);
        if r >= 0 { lo_r = lo_r.min(r); hi_r = hi_r.max(r); }
    }
    if n < min || lo_r == i64::MAX { return false; }
    seg.query_start = lo_q; seg.query_end = hi_q + 1;
    seg.ref_start = lo_r; seg.ref_end = hi_r;
    let mut introns = Vec::new();
    for p in seg.introns.chunks(2) { if p[0] > lo_r && p[1] < hi_r { introns.push(p[0]); introns.push(p[1]); } }
    seg.introns = introns;
    true
}

/// Places the junction of two consecutive parts of a read that overlap by more than `max_gap` bases, or lie more than
/// that apart: an overlap is cut where each part keeps the bases it matches best, unless both match them; a gap is
/// closed where the parts, extended along their own genome, meet. `seqd`: the read as sequenced. Returns what was done.
#[allow(clippy::too_many_arguments)]
fn refine_pair(p: &mut ChimericSegment, pb: Option<&[(i32, i64, bool)]>, q: &mut ChimericSegment, qb: Option<&[(i32, i64, bool)]>,
               seqd: &[u8], genome: &GenomeRef, max_gap: i32, min_keep: usize) -> Option<&'static str> {
    let (pqe, qqs) = (p.query_end, q.query_start);
    if qqs + max_gap < pqe {
        // overlap: the bases from Q's start to P's end
        let (pb, qb) = (pb?, qb?);
        let l = (pqe - qqs) as usize;
        let (mut mp, mut mq) = (vec![0u32; l], vec![0u32; l]);
        for &(x, _, m) in pb { if x >= qqs && x < pqe && m { mp[(x - qqs) as usize] = 1; } }
        for &(x, _, m) in qb { if x >= qqs && x < pqe && m { mq[(x - qqs) as usize] = 1; } }
        let (sp, sq) = (mp.iter().sum::<u32>() as f64 / l as f64, mq.iter().sum::<u32>() as f64 / l as f64);
        // a cut is made only when one part matches 0.85 of the shared bases or more and the other 0.65 at most
        if sp.max(sq) < 0.85 || sp.min(sq) > 0.65 { return None; }
        let (mut best, mut cut, mut acc_p, total_q) = (i64::MIN, 0usize, 0i64, mq.iter().sum::<u32>() as i64);
        let mut acc_q = 0i64;
        for c in 0..=l {
            let sc = acc_p + (total_q - acc_q);
            if sc > best { best = sc; cut = c; }
            if c < l { acc_p += mp[c] as i64; acc_q += mq[c] as i64; }
        }
        let cut = qqs + cut as i32;
        let (p0, q0) = (p.clone(), q.clone());
        if !trim_segment(p, pb, |x| x < cut, min_keep) || !trim_segment(q, qb, |x| x >= cut, min_keep) { *p = p0; *q = q0; return None; }
        return Some("overlap cut");
    }
    let g = qqs - pqe;
    if g <= max_gap || g > 300 { return None; }
    let gap: Vec<u8> = seqd[pqe as usize..qqs as usize].iter().map(|x| x.to_ascii_uppercase()).collect();
    let n = (g + 20) as i64;
    let rc = |v: Vec<u8>| -> Vec<u8> { v.iter().rev().map(|&b| complement(b)).collect() };
    // P's genome beyond its end, read along the sequenced read; Q's before its start, read backwards from it
    let tp = if !p.reverse { genome.fetch_seq(&p.chromosome, p.ref_end, p.ref_end + n).ok()? }
             else { rc(genome.fetch_seq(&p.chromosome, (p.ref_start - 1 - n).max(0), p.ref_start - 1).ok()?) };
    let tq = if !q.reverse { let mut v = genome.fetch_seq(&q.chromosome, (q.ref_start - 1 - n).max(0), q.ref_start - 1).ok()?; v.reverse(); v }
             else { genome.fetch_seq(&q.chromosome, q.ref_end, q.ref_end + n).ok()?.iter().map(|&b| complement(b)).collect() };
    let mut gap_rev = gap.clone(); gap_rev.reverse();
    let (sp, _, qp, _, tp1) = extend(&gap, &tp, 3, 16);
    let (sq, _, qq, _, tq1) = extend(&gap_rev, &tq, 3, 16);
    let qp = if qp >= 10 && sp >= qp as i32 { qp } else { 0 };
    let qq = if qq >= 10 && sq >= qq as i32 { qq } else { 0 };
    if (qp + qq) as i32 + 10 < g { return None; }
    // both extensions covering the whole gap: no cut
    if qp as i32 + 10 >= g && qq as i32 + 10 >= g { return None; }
    let over = (qp + qq) as i32 - g;
    let cut = pqe + qp as i32 - over.max(0) / 2;
    let take_p = (cut - pqe).max(0) as usize;
    let take_q = (qqs - cut).max(0) as usize;
    let dp = if qp > 0 { (tp1 * take_p / qp) as i64 } else { 0 };
    let dq = if qq > 0 { (tq1 * take_q / qq) as i64 } else { 0 };
    p.query_end = cut;
    if !p.reverse { p.ref_end += dp } else { p.ref_start -= dp }
    q.query_start = cut;
    if !q.reverse { q.ref_start -= dq } else { q.ref_end += dq }
    Some("gap closed")
}

/// A read part to replace: (dataset, read, its place on the read).
pub type PartKey = (String, String, i32, i32);

/// The parts of reads the aligner left unaligned that continue into annotated exons, as aligned parts: a part on
/// another gene's exons is one more part of the read (for a read with no other alignment, its alignment too); a part
/// on the aligned gene's own next exons extends the alignment, which is returned in place of the part in `known` it
/// replaces. `known` are the reads' parts already in hand.
pub fn place(db: &PostDB, known: &[ChimericSegment], o: &PlaceOptions, genome_path: &Path) -> (Vec<ChimericSegment>, Vec<PartKey>) {
    if o.bams.is_empty() { return (Vec::new(), Vec::new()); }
    let index = build_index(db, genome_path, o);
    log::info!("[Fusion] Placing unaligned read parts: {} annotated exon boundaries read outward, {} merged exons ({:.0} Mb)",
               index.targets.len(), index.store.exons.len(), index.store.seq.len() as f64 / 1e6);
    // each read's parts in hand, along the read; and, per chromosome, the ends of those parts with `min_len` unaligned
    // bases beyond them (a supplementary alignment's end is searched when another part ends there too)
    let mut parts: FxHashMap<(&str, &str), Vec<usize>> = FxHashMap::default();
    for (i, s) in known.iter().enumerate() { parts.entry((s.dataset.as_str(), s.read_name.as_str())).or_default().push(i); }
    let mut kends: FxHashMap<&str, (Vec<i64>, Vec<i64>)> = FxHashMap::default();
    for v in parts.values_mut() {
        v.sort_by_key(|&i| (known[i].query_start, known[i].query_end));
        for k in 0..v.len() {
            let s = &known[v[k]];
            let right = if k + 1 < v.len() { known[v[k + 1]].query_start } else { s.read_length } - s.query_end;
            let left = s.query_start - if k > 0 { known[v[k - 1]].query_end } else { 0 };
            let (high_gap, low_gap) = if s.reverse { (left, right) } else { (right, left) };
            let e = kends.entry(s.chromosome.as_str()).or_default();
            if high_gap >= o.min_len as i32 { e.0.push(s.ref_end); }
            if low_gap >= o.min_len as i32 { e.1.push(s.ref_start); }
        }
    }
    for e in kends.values_mut() { e.0.sort_unstable(); e.1.sort_unstable(); }
    let mut jobs: Vec<(usize, String)> = Vec::new();
    for (bi, (_, p)) in o.bams.iter().enumerate() {
        let Ok(r) = bam::IndexedReader::from_path(p) else { continue };
        for name in r.header().target_names() {
            let c = String::from_utf8_lossy(name).to_string();
            if index.bounds.contains_key(&c) || index.store.trees.contains_key(&c) { jobs.push((bi, c)); }
        }
    }
    let w = o.window;
    let started = std::time::Instant::now();
    // ISOFORGE_PLACE_CHROMS=<chrom,chrom> limits the search to those chromosomes.
    if let Ok(only) = std::env::var("ISOFORGE_PLACE_CHROMS") {
        let only: Vec<&str> = only.split(',').collect();
        jobs.retain(|(_, c)| only.contains(&c.as_str()));
    }
    // ISOFORGE_PLACE_TRACE=<read name prefix> logs what became of those reads' unaligned pieces
    let trace: Option<String> = std::env::var("ISOFORGE_PLACE_TRACE").ok();
    // ISOFORGE_REFINE=1 places the junction of overlapping or separated parts on the bases (`refine_pair`)
    let refine = std::env::var("ISOFORGE_REFINE").is_ok();
    let genes: Vec<String> = db.genes.iter().map(|g| db.annot.gene_name_map.get(&g.gene_id).cloned().unwrap_or_else(|| g.gene_id.to_string())).collect();
    let found: Vec<(Vec<ChimericSegment>, Vec<PartKey>)> = jobs.par_iter().map(|(bi, chrom)| {
        let (mut out, mut replaced) = (Vec::new(), Vec::new());
        let (ds, path) = &o.bams[*bi];
        let at = index.bounds.get(chrom);
        let Ok(mut reader) = bam::IndexedReader::from_path(path) else { return (out, replaced) };
        let mut fetcher = bam::IndexedReader::from_path(path).ok();
        let names: Vec<String> = reader.header().target_names().iter().map(|n| String::from_utf8_lossy(n).to_string()).collect();
        if reader.fetch(chrom.as_str()).is_err() { return (out, replaced); }
        let mut rec = bam::Record::new();
        let mut hits: Vec<u64> = Vec::new();
        let genome = GenomeRef::open(genome_path.to_path_buf()).ok();
        // A first pass over the chromosome gathers the alignment ends with unaligned bases beyond (an end is looked at
        // only when another read's ends within `window` bp of it) and counts the words of the pieces (see `Artefacts`).
        let (mut ends_high, mut ends_low): (Vec<i64>, Vec<i64>) = (Vec::new(), Vec::new());
        let mut artefacts = Artefacts::default();
        let mut scratch: Vec<u32> = Vec::new();
        // the short clips at alignment ends, for the short classes
        let mut short_ends: Vec<ShortEnd> = Vec::new();
        while let Some(Ok(())) = reader.read(&mut rec) {
            if rec.is_unmapped() || rec.is_secondary() || rec.is_supplementary() || rec.mapq() < o.min_mapq { continue; }
            let cig = rec.cigar();
            let lead = match cig.first() { Some(Cigar::SoftClip(n)) => *n as usize, _ => 0 };
            let trail = match cig.last() { Some(Cigar::SoftClip(n)) => *n as usize, _ => 0 };
            if (SHORT_CLIP_MIN..=SHORT_CLIP_MAX).contains(&trail) || (SHORT_CLIP_MIN..=SHORT_CLIP_MAX).contains(&lead) {
                let sq = rec.seq().as_bytes();
                let n = sq.len();
                if n >= lead + trail {
                    let (h, rp) = (name_hash(rec.qname()), rec.pos());
                    let strand = read_strand(&rec, &sq, lead, trail);
                    let traced = trace.as_deref().is_some_and(|t| rec.qname().starts_with(t.as_bytes()));
                    if (SHORT_CLIP_MIN..=SHORT_CLIP_MAX).contains(&trail) {
                        short_ends.push(ShortEnd { high: true, pos: cig.end_pos(), clip: sq[n - trail..].iter().map(|b| b.to_ascii_uppercase()).collect(), ca: n - trail, cb: n, name: h, rec_pos: rp, strand, traced });
                    }
                    if (SHORT_CLIP_MIN..=SHORT_CLIP_MAX).contains(&lead) {
                        short_ends.push(ShortEnd { high: false, pos: rec.pos() + 1, clip: sq[..lead].iter().map(|b| b.to_ascii_uppercase()).collect(), ca: 0, cb: lead, name: h, rec_pos: rp, strand, traced });
                    }
                }
            }
            if lead < o.min_len && trail < o.min_len { continue; }
            if lead >= o.min_len { ends_low.push(rec.pos() + 1); }
            if trail >= o.min_len { ends_high.push(cig.end_pos()); }
            let seq = rec.seq().as_bytes();
            if seq.len() < lead + trail { continue; }
            if lead >= o.min_len {
                let x: Vec<u8> = seq[lead.saturating_sub(ARTEFACT_SPAN)..lead].iter().map(|b| b.to_ascii_uppercase()).collect();
                artefacts.add(&x, &mut scratch);
            }
            if trail >= o.min_len {
                let a = seq.len() - trail;
                let x: Vec<u8> = seq[a..(a + ARTEFACT_SPAN).min(seq.len())].iter().map(|b| b.to_ascii_uppercase()).collect();
                artefacts.add(&x, &mut scratch);
            }
        }
        artefacts.finish();
        ends_high.sort_unstable(); ends_low.sort_unstable();
        let recurs = |v: &Vec<i64>, p: i64| v.partition_point(|&x| x <= p + w) - v.partition_point(|&x| x < p - w) >= 2;
        // the ends whose piece was searched, and each placed part's end: for the pieces a placed part competes with
        let mut searched: Vec<(bool, i64)> = Vec::new();
        let mut placed_at: Vec<(usize, bool, i64)> = Vec::new();
        if reader.fetch(chrom.as_str()).is_err() { return (out, replaced); }
        // under the trace, a traced read's alignment that is passed over says why
        let skipped = |rec: &bam::Record, why: &str| {
            if let Some(t) = trace.as_deref() {
                let name = String::from_utf8_lossy(rec.qname());
                if name.starts_with(t) { log::info!("[Place trace] {} passed over at {}:{}: {}", name, chrom, rec.pos() + 1, why); }
            }
        };
        let near_in = |at: Option<&FxHashMap<i64, Vec<(u32, bool, bool)>>>, p: i64| -> Vec<(u32, bool, bool, i64)> {
            (-w..=w).filter_map(|d| at.and_then(|at| at.get(&(p + d))).map(|v| v.iter().map(move |&(g, pl, st)| (g, pl, st, p + d)))).flatten().collect()
        };
        let count_near = |v: &[i64], p: i64| v.partition_point(|&x| x <= p + w) - v.partition_point(|&x| x < p - w);
        // One end of one of the read's parts, with the piece beyond it on the record: `right` is the piece's side on
        // the record, `high` the end's side on the genome (they differ for a part aligned on the other strand from the
        // primary), `same` whether the part and the primary are aligned on one strand.
        struct EndJob { part: Option<usize>, chrom: String, pos: i64, high: bool, right: bool, same: bool, part_rev: bool, a: usize, b: usize,
                        bounds: Vec<(u32, bool, bool, i64)>, inside: Vec<(u32, bool)> }
        // the low-MAPQ parts by the end of the primary alignment they lie beyond: (high end, position)
        let mut low_classes: FxHashMap<(bool, i64), Vec<LowPart>> = FxHashMap::default();
        // the lone reads whose alignment the main pass wrote out (with a part placed on another gene)
        let mut based: FxHashSet<(u64, i64)> = FxHashSet::default();
        while let Some(Ok(())) = reader.read(&mut rec) {
            if rec.is_unmapped() || rec.is_secondary() || rec.is_supplementary() { continue; }
            let cig = rec.cigar();
            let lead = match cig.first() { Some(Cigar::SoftClip(n)) => *n as usize, _ => 0 };
            let trail = match cig.last() { Some(Cigar::SoftClip(n)) => *n as usize, _ => 0 };
            let qname = std::str::from_utf8(rec.qname()).unwrap_or("");
            let my_parts: &[usize] = parts.get(&(ds.as_str(), qname)).map_or(&[], |v| &v[..]);
            // a part the aligner could not place among copies, beyond an end of this confident alignment: gathered
            // with the other reads' at that end, to be named together after the pass
            if rec.mapq() >= o.min_mapq && my_parts.iter().any(|&i| known[i].supplementary && known[i].mapq < o.min_mapq && !known[i].placed) {
                let sq = rec.seq().as_bytes();
                let (ln, rv) = (sq.len(), rec.is_reverse());
                for &i in my_parts {
                    let s = &known[i];
                    if !s.supplementary || s.mapq >= o.min_mapq || s.placed { continue; }
                    let (a, b) = { let (a, b) = (s.query_start.max(0) as usize, (s.query_end.max(0) as usize).min(ln)); if rv { (ln - b, ln - a) } else { (a, b) } };
                    if b < a + o.min_len { continue; }
                    let right = a + 10 >= ln.saturating_sub(trail);
                    let left = b <= lead + 10;
                    if right == left { continue; }
                    let high = right != rv;
                    let pos = if high { cig.end_pos() } else { rec.pos() + 1 };
                    let inner: Vec<u8> = if right { sq[a.saturating_sub(BACK)..a].iter().map(|x| x.to_ascii_uppercase()).collect() }
                        else { sq[b..(b + BACK).min(ln)].iter().map(|x| x.to_ascii_uppercase()).collect() };
                    low_classes.entry((high, pos / 10)).or_default().push(LowPart {
                        known: i, seq: sq[a..b].iter().map(|x| x.to_ascii_uppercase()).collect(), ra: a, rev: rv, len: ln, name: qname.to_string(), pos, right, inner });
                }
            }
            if lead < o.min_len && trail < o.min_len { continue; }
            if rec.mapq() < o.min_mapq { skipped(&rec, &format!("mapping quality {} with {} and {} unaligned bases", rec.mapq(), lead, trail)); continue; }
            let (first, last) = (rec.pos() + 1, cig.end_pos());
            let (try_high, try_low) = (trail >= o.min_len && recurs(&ends_high, last), lead >= o.min_len && recurs(&ends_low, first));
            if !try_high && !try_low && my_parts.is_empty() { skipped(&rec, &format!("no other read's alignment ends there with unaligned bases beyond ({} and {} unaligned)", lead, trail)); continue; }
            let (b_high, b_low) = (if try_high { near_in(at, last) } else { Vec::new() }, if try_low { near_in(at, first) } else { Vec::new() });
            let (g_high, g_low) = (if try_high { index.store.genes_at(chrom, last) } else { Vec::new() },
                                   if try_low { index.store.genes_at(chrom, first) } else { Vec::new() });
            if b_high.is_empty() && b_low.is_empty() && g_high.is_empty() && g_low.is_empty() && my_parts.is_empty() {
                skipped(&rec, &format!("the alignment does not end in an annotated exon ({} and {} unaligned)", lead, trail));
                continue;
            }
            let name = qname.to_string();
            let traced = trace.as_deref().is_some_and(|t| name.starts_with(t));
            DEBUG.with(|d| d.set(traced && std::env::var("ISOFORGE_LOCATE_DEBUG").is_ok()));
            let seq = rec.seq().as_bytes();
            let len = seq.len();
            if len < lead + trail || len == 0 { continue; }
            let rev = rec.is_reverse();
            // every part of the read in the record's orientation, the primary among them
            let on_record = |a: i32, b: i32| -> (usize, usize) {
                let (a, b) = (a.max(0) as usize, (b.max(0) as usize).min(len));
                if rev { (len - b, len - a) } else { (a, b) }
            };
            let mut ranges: Vec<(usize, usize)> = my_parts.iter().map(|&i| on_record(known[i].query_start, known[i].query_end)).collect();
            ranges.push((lead, len - trail));
            let next_start = |p: usize| ranges.iter().filter(|x| x.0 + 10 >= p && x.1 > p + 10).map(|x| x.0.max(p)).min().unwrap_or(len);
            let prev_end = |p: usize| ranges.iter().filter(|x| x.1 <= p + 10 && x.0 + 10 < p).map(|x| x.1.min(p)).max().unwrap_or(0);
            let mut jobs: Vec<EndJob> = Vec::new();
            for (high, bounds, inside) in [(true, b_high, g_high), (false, b_low, g_low)] {
                if bounds.is_empty() && inside.is_empty() { continue; }
                let (a, b) = if high { (len - trail, next_start(len - trail)) } else { (prev_end(lead), lead) };
                jobs.push(EndJob { part: None, chrom: chrom.clone(), pos: if high { last } else { first }, high, right: high, same: true, part_rev: rev, a, b, bounds, inside });
            }
            // the ends of the read's supplementary alignments: their pieces are in this record's sequence
            for (k, &i) in my_parts.iter().enumerate() {
                let s = &known[i];
                if !s.supplementary || s.mapq < o.min_mapq || s.placed { continue; }
                let (sa, sb) = ranges[k];
                let same = s.reverse == rev;
                let ke = kends.get(s.chromosome.as_str());
                let at_s = index.bounds.get(&s.chromosome);
                for right in [true, false] {
                    let high = right == same;
                    let (a, b) = if right { (sb, next_start(sb)) } else { (prev_end(sa), sa) };
                    if b < a + o.min_len { continue; }
                    // the bases beyond this end are unaligned only when no other part of the read covers them
                    let edge = if right { sb } else { sa };
                    if ranges.iter().any(|x| x.0 + 10 < edge && x.1 > edge + 10) { continue; }
                    let pos = if high { s.ref_end } else { s.ref_start };
                    let recurring = ke.is_some_and(|e| count_near(if high { &e.0 } else { &e.1 }, pos) >= 2);
                    if !recurring { continue; }
                    let bounds = near_in(at_s, pos);
                    let inside = index.store.genes_at(&s.chromosome, pos);
                    if bounds.is_empty() && inside.is_empty() { continue; }
                    jobs.push(EndJob { part: Some(i), chrom: s.chromosome.clone(), pos, high, right, same, part_rev: s.reverse, a, b, bounds, inside });
                }
            }
            if jobs.is_empty() { continue; }
            let Some(mut base) = ChimericSegment::of_record(&rec, &names, ds) else { continue };
            let was = (base.query_start, base.query_end);
            let (mut placed_any, mut extended) = (false, false);
            // the parts' alignments extended over their own exons, by part
            let mut ext: FxHashMap<usize, ChimericSegment> = FxHashMap::default();
            // pieces placed on another gene, on the record: a gap searched from both of its flanks is placed once
            let mut placed_ranges: Vec<(usize, usize)> = Vec::new();
            for job in jobs {
                let EndJob { part, chrom: pchrom, pos, high, right, same, part_rev, a, b, bounds, inside } = job;
                let pchrom = pchrom.as_str();
                let at_p = index.bounds.get(pchrom);
                if b < a + o.min_len {
                    if traced { log::info!("[Place trace] {} {} end: {} unaligned bases before the read's next part", name, if high { "high" } else { "low" }, b.saturating_sub(a)); }
                    continue;
                }
                if placed_ranges.iter().any(|&(x, y)| y.min(b) > x.max(a) + 10) { continue; }
                let mut local: Vec<u32> = bounds.iter().map(|x| x.0).chain(inside.iter().map(|x| x.0)).collect();
                local.sort_unstable(); local.dedup();
                // (steps, gene strand, piece on the record, boundary)
                let mut best: Option<(Vec<Step>, bool, usize, usize, i64)> = None;
                let mut counted = false;
                for plus in [true, false] {
                    // beyond an exon's 3' boundary on its gene's strand, the piece runs 3'-ward; beyond its 5' one, 5'-ward
                    let forward = high == plus;
                    // the read's bases from the boundary: the alignment's end is `d` short of it, or past it; without a
                    // boundary near, from the alignment's end inside an exon of a gene of this strand
                    let mut shifts: Vec<(i64, bool)> = bounds.iter().filter(|x| x.1 == plus && x.2 != forward)
                        .map(|x| (if high { x.3 - pos } else { pos - x.3 }, true)).collect();
                    // an alignment may run past the boundary by `free` bases at most
                    shifts.retain(|x| x.0 >= -(o.free as i64));
                    shifts.sort_unstable(); shifts.dedup();
                    if shifts.is_empty() && inside.iter().any(|x| x.1 == plus) {
                        // Without a boundary near: a boundary up to 60 bp on is taken when the piece's first bases
                        // are the genome's up to it (aligned to their end, scoring one a base or more).
                        if let (Some(at), Some(g)) = (at_p, genome.as_ref()) {
                            for d in (w + 1)..=60 {
                                let p = if high { pos + d } else { pos - d };
                                if !at.get(&p).is_some_and(|v| v.iter().any(|x| x.1 == plus && x.2 != forward)) { continue; }
                                let (lo, hi) = if high { (pos + 1, p) } else { (p, pos - 1) };
                                let Ok(mut gs) = g.fetch_seq(pchrom, lo - 1, hi) else { continue };
                                if !high { gs.reverse(); }
                                let dd = d as usize;
                                // the read's bases next to the alignment, read outward, in the part's reference orientation
                                let mut mine: Vec<u8> = if right { seq.get(a..a + dd).unwrap_or(&[]).iter().map(|x| x.to_ascii_uppercase()).collect() }
                                    else { seq.get(b.saturating_sub(dd)..b).unwrap_or(&[]).iter().rev().map(|x| x.to_ascii_uppercase()).collect() };
                                if !same { for x in mine.iter_mut() { *x = complement(*x); } }
                                if mine.len() == dd && gs.len() == dd {
                                    let (sc, _, q1, _, _) = extend(&mine, &gs, o.free, o.free + 8);
                                    if q1 + o.free >= dd && sc >= d as i32 { shifts.push((d, true)); break; }
                                }
                            }
                        }
                        if shifts.is_empty() { shifts.push((0, false)); }
                    }
                    for (d, at_boundary) in shifts {
                        let (pa, pb) = if right { (a as i64 + d, b as i64) } else { (a as i64, b as i64 - d) };
                        if pa < 0 || pb > len as i64 || pb < pa + o.min_len as i64 { continue; }
                        let (pa, pb) = (pa as usize, pb as usize);
                        let piece: Vec<u8> = seq[pa..pb].iter().map(|x| x.to_ascii_uppercase()).collect();
                        // the piece in the part's reference orientation
                        let piece_ref: Vec<u8> = if same { piece.clone() } else { piece.iter().rev().map(|&x| complement(x)).collect() };
                        let (mut q, _) = outward(&piece_ref, high, plus);
                        // the piece ends where its poly(A) tail starts; a short or low-complexity piece is left out
                        q.truncate(before_tail(&q));
                        if q.len() < o.min_len || low_complexity(&q) {
                            if traced { log::info!("[Place trace] {} {} end: a tail or a simple repeat ({} of {} bases before it)", name, if high { "high" } else { "low" }, q.len(), piece.len()); }
                            continue;
                        }
                        // a piece holding words that recur across the reads' pieces is left out
                        let at_end: &[u8] = if right { &piece[..piece.len().min(ARTEFACT_SPAN)] } else { &piece[piece.len().saturating_sub(ARTEFACT_SPAN)..] };
                        if artefacts.holds(at_end) {
                            N_ARTEFACT.fetch_add(1, Relaxed);
                            if traced { log::info!("[Place trace] {} {} end: recurring words (an adapter) in the piece's {} bases", name, if high { "high" } else { "low" }, piece.len()); }
                            continue;
                        }
                        N_PIECES.fetch_add(1, Relaxed);
                        if !counted && part.is_none() { searched.push((high, pos)); counted = true; }
                        let t0 = std::time::Instant::now();
                        let mut steps = index.chain(&q, forward, &local, at_boundary, o, &mut hits);
                        T_CHAIN.fetch_add(t0.elapsed().as_nanos() as u64, Relaxed);
                        // the first part on another gene: when the read's bases before it continue into that gene too,
                        // or along the genome beside the part, nothing is placed
                        let mut copy_of = None;
                        let mut own_locus = false;
                        if let Some(st) = steps.iter().find(|x| !x.own) {
                            let inner: &[u8] = if right { &seq[pa.saturating_sub(BACK)..pa] } else { &seq[pb..(pb + BACK).min(len)] };
                            let inner: Vec<u8> = if same { inner.iter().map(|x| x.to_ascii_uppercase()).collect() }
                                else { inner.iter().rev().map(|x| complement(x.to_ascii_uppercase())).collect() };
                            let mut before = outward(&inner, high, plus).0;
                            before.extend_from_slice(&q[..st.off + st.q0]);
                            before.reverse();
                            before.truncate(BACK);
                            if index.continues_back(&before, forward, st.gene, o, &mut hits)
                                || genome.as_ref().is_some_and(|g| index.genome_continues_back(g, &before, forward, st, o)) { copy_of = Some(st.gene); }
                            // nor when the part lies in the genome beyond the alignment (`in_own_locus`)
                            else {
                                // beyond the boundary, or beyond the gene's own exons the piece ran through first
                                let mut edge = if high { pos + d } else { pos - d };
                                for x in steps.iter().take_while(|x| x.own) {
                                    if let (Some(a), Some(b)) = (x.blocks.first(), x.blocks.last()) { edge = if high { edge.max(b.1) } else { edge.min(a.0) }; }
                                }
                                let rest: &[u8] = if right { &piece[st.off.min(piece.len())..] } else { &piece[..piece.len().saturating_sub(st.off)] };
                                let rest: Vec<u8> = if same { rest.to_vec() } else { rest.iter().rev().map(|&x| complement(x)).collect() };
                                if genome.as_ref().is_some_and(|g| index.in_own_locus(g, pchrom, edge, high, plus, &rest, &local, st.score, OWN_REACH, o)) { own_locus = true; }
                            }
                        }
                        if traced {
                            let shown: Vec<String> = steps.iter().map(|st| format!("{}{}{} score {} piece {}-{} on {}:{:?}", if st.own { "own exon of " } else { "" },
                                genes[st.gene as usize], if st.anchored { "" } else { " (inside an exon)" }, st.score, st.off + st.q0, st.off + st.q1,
                                index.chroms[st.chrom as usize], st.blocks)).collect();
                            let miss = if steps.is_empty() { format!("{:?} / {:?}", index.place(&q[..q.len().min(o.max_len)], forward, &local, None, o, &mut hits).err().map(|e| (e.0, e.1.map(|p| p.score))),
                                index.store.locate(&q[..q.len().min(o.max_len)], forward, &local, None, o, &mut hits).err()) } else { String::new() };
                            if own_locus { log::info!("[Place trace] {} the piece lies in the genome beyond the alignment too: an exon of the aligned gene", name); }
                            if let Some(g) = copy_of { log::info!("[Place trace] {} the read's bases before the part placed on {} continue there too: the read lies on a copy of one gene", name, genes[g as usize]); }
                            log::info!("[Place trace] {} {} end at {}:{}{} gene strand {} shift {} ({}) piece {} bases ({}): {:?} {}", name, if high { "high" } else { "low" },
                                       pchrom, pos, if part.is_some() { " (supplementary alignment)" } else { "" }, if plus { '+' } else { '-' }, d,
                                       if at_boundary { "exon boundary" } else { "inside an exon" }, piece.len(), String::from_utf8_lossy(&q[..q.len().min(40)]), shown, miss);
                        }
                        if copy_of.is_some() || own_locus { steps.clear(); }
                        if steps.is_empty() { continue; }
                        // A chain that starts on the gene's own exon is taken unless one starting on another gene
                        // scores `2 * margin` more on its first part.
                        let rank = |v: &Vec<Step>| (if v[0].own { v[0].score + 2 * o.margin } else { v[0].score }, v.iter().map(|x| x.score).sum::<i32>());
                        if best.as_ref().map_or(true, |x| rank(&steps) > rank(&x.0)) {
                            best = Some((steps, plus, pa, pb, if high { pos + d } else { pos - d }));
                        }
                    }
                }
                let Some((steps, plus, pa, pb, boundary)) = best else { continue };
                // the read runs along the aligned gene's transcript or against it, and along a placed gene's alike
                let sense = part_rev != plus;
                let mut edge = boundary;
                let mut on_base = true;
                // the part's own alignment, to extend over its gene's exons
                let own_seg: &mut ChimericSegment = match part {
                    None => &mut base,
                    Some(i) => ext.entry(i).or_insert_with(|| known[i].clone()),
                };
                let mut own_extended = false;
                for st in steps {
                    if st.blocks.is_empty() { break; }
                    // on the record, the placed bases run outward from the boundary
                    let (ra, rb) = if right { (pa + st.off + st.q0, pa + st.off + st.q1) } else { (pb - st.off - st.q1, pb - st.off - st.q0) };
                    let (qs, qe) = if rev { ((len - rb) as i32, (len - ra) as i32) } else { (ra as i32, rb as i32) };
                    if st.own {
                        // the aligned gene's own exons extend the alignment; none is taken after another gene's part
                        if !on_base || index.chroms[st.chrom as usize] != *pchrom || !extend_segment(own_seg, high, edge, &st.blocks, (qs, qe)) { break; }
                        edge = if high { st.blocks[st.blocks.len() - 1].1 } else { st.blocks[0].0 };
                        own_extended = true;
                        continue;
                    }
                    on_base = false;
                    let mut introns = Vec::new();
                    for x in st.blocks.windows(2) { introns.push(x[0].1 + 1); introns.push(x[1].0 - 1); }
                    out.push(ChimericSegment {
                        read_name: name.clone(), dataset: ds.clone(), chromosome: index.chroms[st.chrom as usize].clone(),
                        reverse: sense != st.plus, supplementary: true, mapq: 60, ref_start: st.blocks[0].0, ref_end: st.blocks[st.blocks.len() - 1].1,
                        query_start: qs, query_end: qe, read_length: len as i32, introns, placed: true, short_class: false, pieces: 0,
                        excess: st.score - if st.anchored { o.min_score } else { o.min_score_free },
                    });
                    match part {
                        None => placed_at.push((out.len() - 1, high, pos)),
                        // the pieces at a supplementary alignment's end: the parts in hand ending there with unaligned bases beyond
                        Some(_) => { let k = out.len() - 1; out[k].pieces = kends.get(pchrom).map_or(1, |e| count_near(if high { &e.0 } else { &e.1 }, pos).max(1) as u32); }
                    }
                    placed_ranges.push((ra.min(rb), ra.max(rb)));
                    placed_any = true;
                }
                if own_extended && part.is_none() { extended = true; }
            }
            // consecutive parts that overlap on the read, or stand apart, have their junction placed on the bases
            if refine && !my_parts.is_empty() {
                if let (Some(g), Some(f)) = (genome.as_ref(), fetcher.as_mut()) {
                    // the parts in sequenced order by their current extents; the primary alignment (its record in hand,
                    // its segment `base`) is among them, as -1
                    let mut order: Vec<(i32, i32, i64)> = my_parts.iter().map(|&i| {
                        if !known[i].supplementary { (base.query_start, base.query_end, -1) }
                        else { let s = ext.get(&i).unwrap_or(&known[i]); (s.query_start, s.query_end, i as i64) } }).collect();
                    if !order.iter().any(|x| x.2 < 0) { order.push((base.query_start, base.query_end, -1)); }
                    order.sort_unstable();
                    let seqd: Vec<u8> = if rev { seq.iter().rev().map(|&b| complement(b.to_ascii_uppercase())).collect() } else { seq.iter().map(|b| b.to_ascii_uppercase()).collect() };
                    let mut bases_cache: FxHashMap<i64, Option<Vec<(i32, i64, bool)>>> = FxHashMap::default();
                    for k in 0..order.len().saturating_sub(1) {
                        let (pi, qi) = (order[k].2, order[k + 1].2);
                        let (pqe, qqs) = (order[k].1, order[k + 1].0);
                        let apart = qqs - pqe;
                        if apart.abs() <= o.max_query_gap as i32 || apart > 300 { continue; }
                        // a gap holding a placed part is skipped
                        if apart > 0 {
                            let (ga, gb) = if rev { ((len as i32 - qqs) as usize, (len as i32 - pqe) as usize) } else { (pqe as usize, qqs as usize) };
                            if placed_ranges.iter().any(|&(x, y)| y.min(gb) > x.max(ga) + 10) { continue; }
                        }
                        let mut take = |i: i64| -> ChimericSegment { if i < 0 { base.clone() } else { ext.remove(&(i as usize)).unwrap_or_else(|| known[i as usize].clone()) } };
                        let (mut sp, mut sq) = (take(pi), take(qi));
                        let need_bases = apart < 0;
                        let mut bases_of = |i: i64| -> Option<Vec<(i32, i64, bool)>> {
                            if !need_bases { return None; }
                            bases_cache.entry(i).or_insert_with(|| {
                                if i < 0 { Some(aligned_bases(&rec, g, chrom)) } else {
                                    let s = &known[i as usize];
                                    fetch_supplementary(f, &s.chromosome, s.ref_start, rec.qname()).map(|r| aligned_bases(&r, g, &s.chromosome))
                                }
                            }).clone()
                        };
                        let (bp, bq) = (bases_of(pi), bases_of(qi));
                        let done = refine_pair(&mut sp, bp.as_deref(), &mut sq, bq.as_deref(), &seqd, g, o.max_query_gap as i32, o.min_len);
                        if traced { log::info!("[Place trace] {} parts {}-{} and {}-{} on the read, {} bases {}: {}", name, order[k].0, pqe, qqs, order[k + 1].1,
                                               apart.abs(), if apart < 0 { "shared" } else { "between" }, done.unwrap_or("left as aligned")); }
                        if done.is_some() {
                            order[k].1 = sp.query_end; order[k + 1].0 = sq.query_start;
                            if pi < 0 { extended = true; }
                            if qi < 0 { extended = true; }
                        }
                        for (i, sg) in [(pi, sp), (qi, sq)] { if i < 0 { base = sg; } else { ext.insert(i as usize, sg); } }
                    }
                }
            }
            // a read with no other part: its alignment, when a part was placed on another gene; a read with other
            // parts: its alignment in place of the part in hand, when its gene's own exons extended it
            if my_parts.is_empty() {
                if placed_any { based.insert((name_hash(rec.qname()), rec.pos())); out.push(base); }
            } else if extended {
                replaced.push((ds.clone(), name.clone(), was.0, was.1));
                out.push(base);
            }
            for (i, seg) in ext {
                let k = &known[i];
                if seg.introns == k.introns && seg.query_start == k.query_start && seg.query_end == k.query_end && seg.ref_start == k.ref_start && seg.ref_end == k.ref_end { continue; }
                replaced.push((ds.clone(), name.clone(), known[i].query_start, known[i].query_end));
                out.push(seg);
            }
        }
        // Low-MAPQ parts named by their class: for the reads whose confident alignment ends at one place, every copy
        // each part fits is scored (`candidates`), the best score per gene adds over the reads, and the gene that leads
        // by `CLASS_MARGIN` with `CLASS_MIN_READS` reads fitting it (not a gene of the confident alignment) names them.
        for ((high, _), class) in low_classes {
            if class.len() < CLASS_MIN_READS { continue; }
            let own: Vec<u32> = index.store.genes_at(chrom, class[0].pos).iter().map(|x| x.0).collect();
            let mut per_read: Vec<(Vec<Cand>, FxHashMap<u32, (i32, usize)>)> = Vec::new();
            let mut totals: FxHashMap<u32, (i32, u32, i64, i64)> = FxHashMap::default();
            for p in &class {
                let c = index.store.candidates(&p.seq, o, &mut hits);
                let mut best: FxHashMap<u32, (i32, usize)> = FxHashMap::default();
                for (k, x) in c.iter().enumerate() {
                    let e = best.entry(x.gene).or_insert((x.score, k));
                    if x.score > e.0 { *e = (x.score, k); }
                }
                for (g, (sc, k)) in &best {
                    let x = &c[*k];
                    let e = totals.entry(*g).or_insert((0, x.chrom, x.lo, x.hi));
                    e.0 += sc;
                }
                per_read.push((c, best));
            }
            let mut ranked: Vec<(i32, u32)> = totals.iter().map(|(g, t)| (t.0, *g)).collect();
            ranked.sort_unstable_by(|a, b| b.cmp(a));
            let Some(&(top, g_star)) = ranked.first() else { continue };
            let star = totals[&g_star];
            let second = ranked.iter().skip(1).find(|x| { let t = totals[&x.1]; !(t.1 == star.1 && t.2 <= star.3 && star.2 <= t.3) }).map_or(0, |x| x.0);
            let backing = per_read.iter().filter(|(_, m)| m.get(&g_star).is_some_and(|x| x.0 >= o.min_score)).count();
            let traced_any = trace.as_deref().is_some_and(|t| class.iter().any(|p| p.name.starts_with(t)));
            // whether the confident alignment's own gene fits the parts within `CLASS_MARGIN` of the leader
            let own_near = own.iter().any(|g| totals.get(g).is_some_and(|t| t.0 >= top - CLASS_MARGIN));
            let need = CLASS_MARGIN.max(CLASS_MARGIN_PER_READ * class.len() as i32);
            // the other genes elsewhere that the parts fit six tenths as well as the leader: three make a repeat
            let alike = ranked.iter().filter(|x| x.1 != g_star && 10 * x.0 >= 6 * top && { let t = totals[&x.1]; !(t.1 == star.1 && t.2 <= star.3 && star.2 <= t.3) }).count();
            let why = if own.contains(&g_star) || own_near { Some("the confident alignment's own gene fits as well") } else if alike >= 3 { Some("a repeat: three other genes fit it six tenths as well") }
                else if top - second < need { Some("another gene scores alike over the reads") }
                else if backing < CLASS_MIN_READS { Some("too few reads fit it") } else { None };
            if traced_any {
                for p in &class {
                    if trace.as_deref().is_some_and(|t| p.name.starts_with(t)) {
                        log::info!("[Place trace] {} {} end at {}:{}: a part below the mapping-quality floor, {} reads at this end, best {} {} over the reads, runner-up {}: {}", p.name,
                                   if high { "high" } else { "low" }, chrom, p.pos, class.len(), genes[g_star as usize], top, second, why.unwrap_or("named by the class"));
                    }
                }
            }
            if why.is_some() { continue; }
            for (p, (c, m)) in class.iter().zip(&per_read) {
                let Some(&(sc, k)) = m.get(&g_star) else { continue };
                if sc < o.min_score { continue; }
                let x = &c[k];
                // the read's bases before the part, in the named gene's transcript orientation, nearest the part
                // first, taken when the transcript runs from the primary into the part: when they continue into that
                // gene too (its exons, or the genome beside the part), the part is not named
                let before: Option<Vec<u8>> = match (x.plus, p.right) {
                    (true, true) => Some(p.inner.iter().rev().copied().collect()),
                    (false, false) => Some(p.inner.iter().map(|&b| complement(b)).collect()),
                    _ => None,
                };
                if let Some(before) = before {
                    let st = Step { off: 0, q0: x.q0, q1: x.q1, score: sc, gene: x.gene, chrom: x.chrom, plus: x.plus, blocks: vec![(x.lo, x.hi)], own: false, at_end: false, anchored: false };
                    if index.continues_back(&before, true, x.gene, o, &mut hits) || genome.as_ref().is_some_and(|g| index.genome_continues_back(g, &before, true, &st, o)) {
                        if trace.as_deref().is_some_and(|t| p.name.starts_with(t)) { log::info!("[Place trace] {} the read's bases before the part named {} by its class continue there too: the read lies on a copy of one gene", p.name, genes[x.gene as usize]); }
                        continue;
                    }
                }
                // nor when the part is found in the genome beside the confident alignment (within `OWN_REACH`, outside
                // other genes' exons)
                if genome.as_ref().is_some_and(|g| index.in_own_locus(g, chrom, p.pos, high, !p.rev, &p.seq, &own, sc, OWN_REACH, o)) {
                    if trace.as_deref().is_some_and(|t| p.name.starts_with(t)) { log::info!("[Place trace] {} the part named {} by its class lies in the genome beside the confident alignment too: the aligned gene's own", p.name, genes[x.gene as usize]); }
                    continue;
                }
                let (ra, rb) = (p.ra + x.q0, p.ra + x.q1);
                let (qs, qe) = if p.rev { ((p.len - rb) as i32, (p.len - ra) as i32) } else { (ra as i32, rb as i32) };
                out.push(ChimericSegment {
                    read_name: p.name.clone(), dataset: ds.clone(), chromosome: index.chroms[x.chrom as usize].clone(), reverse: if x.plus { p.rev } else { !p.rev },
                    supplementary: true, mapq: 60, ref_start: x.lo, ref_end: x.hi, query_start: qs, query_end: qe, read_length: p.len as i32, introns: Vec::new(),
                    placed: true, short_class: false, pieces: class.len() as u32, excess: top - second - need + o.margin,
                });
                replaced.push((ds.clone(), p.name.clone(), known[p.known].query_start, known[p.known].query_end));
                N_CLASS.fetch_add(1, Relaxed);
            }
        }
        // Short parts placed as a class (see `SHORT_MIN`): the reads whose alignment ends at one place with a short
        // clip beyond it, their pieces scored together on the partner side's exon boundaries and on the transcripts'
        // own ends; a read whose clip is `min_len` or longer has had its piece searched above and is not placed here.
        let mut short_hits: FxHashMap<(u64, i64), ShortHit> = FxHashMap::default();
        if !short_ends.is_empty() {
            // words that recur in the short clips of the chromosome
            let mut words: FxHashMap<u32, u32> = FxHashMap::default();
            for e in &short_ends {
                let mut ws = canon10(&e.clip);
                ws.sort_unstable(); ws.dedup();
                for w in ws { *words.entry(w).or_insert(0) += 1; }
            }
            let floor = (20.0f64).max(0.002 * short_ends.len() as f64) as u32;
            let common: FxHashSet<u32> = words.into_iter().filter(|&(_, n)| n >= floor).map(|(w, _)| w).collect();
            short_ends.sort_by_key(|e| (e.high, e.pos));
            let drift = o.free as i64 + 4;
            let mut i = 0;
            while i < short_ends.len() {
                let mut k = i + 1;
                while k < short_ends.len() && short_ends[k].high == short_ends[i].high && short_ends[k].pos - short_ends[k - 1].pos <= SHORT_WIN { k += 1; }
                let group = &short_ends[i..k];
                i = k;
                let mut names: Vec<u64> = group.iter().map(|e| e.name).collect();
                names.sort_unstable(); names.dedup();
                if names.len() < SHORT_FIT { continue; }
                let high = group[0].high;
                let mut ps: Vec<i64> = group.iter().map(|e| e.pos).collect();
                ps.sort_unstable();
                let pos = ps[ps.len() / 2];
                let gs = index.store.genes_at(chrom, pos);
                if gs.is_empty() { continue; }
                let own: Vec<u32> = gs.iter().map(|x| x.0).collect();
                let mut strands: Vec<bool> = gs.iter().map(|x| x.1).collect();
                strands.sort_unstable(); strands.dedup();
                // the strand the reads give: two of them at least, three for each that says the other
                let (np, nm) = group.iter().fold((0usize, 0usize), |(p, m), e| match e.strand { 1 => (p + 1, m), -1 => (p, m + 1), _ => (p, m) });
                let read_plus = if np >= 2 && np >= 3 * nm { Some(true) } else if nm >= 2 && nm >= 3 * np { Some(false) } else { None };
                for plus in strands {
                    if read_plus.is_some_and(|r| r != plus) {
                        if group.iter().any(|e| e.traced) { log::info!("[Place trace] short class at {}:{} ({} end) not read on the {} strand: {} of its reads are plus, {} minus", chrom, pos, if high { "high" } else { "low" }, if plus { "+" } else { "-" }, np, nm); }
                        continue;
                    }
                    // the partner follows the aligned part in transcript order
                    let three = high == plus;
                    // (outward piece, piece from the read's end, bases of the clip kept, the end)
                    let mut pieces: Vec<(Vec<u8>, Vec<u8>, usize, &ShortEnd)> = Vec::new();
                    for e in group {
                        let t: Vec<u8> = if plus { e.clip.clone() } else { e.clip.iter().rev().map(|&b| complement(b)).collect() };
                        let l = if three { short_tail_cut(&t) } else { t.len() };
                        if l < SHORT_CLIP_MIN { continue; }
                        let tc = &t[..l];
                        let outward: Vec<u8> = if three { tc.to_vec() } else { tc.iter().rev().copied().collect() };
                        if short_low_complexity(&outward) || canon10(tc).iter().any(|w| common.contains(w)) { continue; }
                        let readend: Vec<u8> = if three { tc.iter().rev().copied().collect() } else { tc.to_vec() };
                        pieces.push((outward[..l.min(SHORT_QMAX)].to_vec(), readend[..l.min(SHORT_QMAX)].to_vec(), l, e));
                    }
                    let n = pieces.len();
                    let traced = group.iter().any(|e| e.traced);
                    if n < SHORT_FIT {
                        if traced { log::info!("[Place trace] short class at {}:{} ({} end, {} strand): {} reads, {} pieces left after the tail, complexity and adapter tests", chrom, pos, if high { "high" } else { "low" }, if plus { "+" } else { "-" }, names.len(), n); }
                        continue;
                    }
                    // (margin, target, from the read's end, pooled score)
                    let mut chosen: Vec<(i32, u32, bool, i32)> = Vec::new();
                    // the pieces the decision is taken on
                    let sample: Vec<usize> = if n <= SHORT_SAMPLE { (0..n).collect() } else { (0..SHORT_SAMPLE).map(|k| k * n / SHORT_SAMPLE).collect() };
                    for from_end in [false, true] {
                        let qs: Vec<&Vec<u8>> = sample.iter().map(|&k| if from_end { &pieces[k].1 } else { &pieces[k].0 }).collect();
                        let mut lens: Vec<usize> = pieces.iter().map(|p| if from_end { p.1.len() } else { p.0.len() }).collect();
                        lens.sort_unstable();
                        if lens[lens.len() / 2] < SHORT_MIN { continue; }
                        // read outward: the partner side's boundaries; from the read's end: the transcripts' own ends
                        let fwd = if from_end { !three } else { three };
                        let kidx = if fwd { &index.forward } else { &index.backward };
                        let mut hitc: FxHashMap<u32, u32> = FxHashMap::default();
                        for q in &qs {
                            for (qpos, km) in kmers(&q[..q.len().min(SEED + o.free)]) {
                                let bucket = kidx.get(km);
                                if bucket.len() > SHORT_BUCKET_MAX { continue; }
                                for &e in bucket {
                                    let diag = qpos as i64 - (e & ((1 << POS_BITS) - 1)) as i64;
                                    if diag.abs() <= drift { *hitc.entry(e >> POS_BITS).or_insert(0) += 1; }
                                }
                            }
                        }
                        // one target per gene and first bases (the smallest index), with the seed hits of all of them;
                        // the aligned gene's own targets compete too
                        let mut byhead: FxHashMap<(u32, &[u8]), (u32, u32)> = FxHashMap::default();
                        for (&c, &h) in &hitc {
                            let t = &index.targets[c as usize];
                            if from_end && !t.tx_end { continue; }
                            let e = byhead.entry((t.gene, &t.seq[..t.seq.len().min(SHORT_HEAD)])).or_insert((c, 0));
                            e.0 = e.0.min(c); e.1 += h;
                        }
                        let mut cand: Vec<(u32, u32)> = byhead.into_values().collect();
                        let seeded = cand.len();
                        if cand.is_empty() {
                            if traced { log::info!("[Place trace] short class at {}:{} ({} end, {} strand) read {}: {} pieces, no candidate target", chrom, pos, if high { "high" } else { "low" }, if plus { "+" } else { "-" }, if from_end { "from the read's end" } else { "outward" }, n); }
                            continue;
                        }
                        cand.sort_unstable_by_key(|&(c, h)| (std::cmp::Reverse(h), c));
                        cand.truncate(SHORT_CANDIDATES);
                        let mut cand: Vec<u32> = cand.into_iter().map(|x| x.0).collect();
                        cand.sort_unstable();
                        let alns: Vec<Vec<(i32, usize, usize, usize, usize)>> = cand.iter().map(|&c| {
                            qs.iter().map(|q| extend(q, &index.targets[c as usize].seq, o.free, o.free + 8)).collect()
                        }).collect();
                        let tot: Vec<i32> = alns.iter().map(|v| v.iter().map(|x| x.0).sum()).collect();
                        let top = (0..cand.len()).max_by_key(|&j| (tot[j], std::cmp::Reverse(cand[j]))).unwrap();
                        let tt = &index.targets[cand[top] as usize];
                        let head = &tt.seq[..tt.seq.len().min(16)];
                        let second = (0..cand.len()).filter(|&j| {
                            let u = &index.targets[cand[j] as usize];
                            u.gene != tt.gene && &u.seq[..u.seq.len().min(16)] != head
                        }).map(|j| tot[j]).max().unwrap_or(0);
                        let fit = qs.iter().zip(&alns[top]).filter(|(q, a)| a.0 as f64 >= SHORT_FIT_PER_BASE * q.len() as f64).count();
                        let margin = tot[top] - second;
                        // whether the leading target is the aligned gene's own
                        let own_top = own.contains(&tt.gene);
                        let ok = !own_top && fit >= SHORT_FIT && fit as f64 >= SHORT_SHARE * qs.len() as f64 && margin >= SHORT_MARGIN;
                        if traced {
                            log::info!("[Place trace] short class at {}:{} ({} end, {} strand) read {}: {} pieces ({} scored, median {} bases), {} candidates ({} scored), top {} {}:{:?} pooled {}, runner-up {}, {} fitting: {}",
                                       chrom, pos, if high { "high" } else { "low" }, if plus { "+" } else { "-" }, if from_end { "from the read's end" } else { "outward" },
                                       n, qs.len(), lens[lens.len() / 2], seeded, cand.len(), genes[tt.gene as usize], index.chroms[tt.chrom as usize], tt.blocks.first(), tot[top], second, fit,
                                       if ok { "accepted" } else if own_top { "rejected, the aligned gene's own exon" } else { "rejected" });
                        }
                        if ok { chosen.push((margin, cand[top], from_end, tot[top])); }
                    }
                    // both readings naming different genes: no call
                    if chosen.len() == 2 && index.targets[chosen[0].1 as usize].gene != index.targets[chosen[1].1 as usize].gene { continue; }
                    let Some((margin, target, from_end, pooled)) = chosen.into_iter().max_by_key(|x| x.0) else { continue };
                    // the genome beyond the alignment end competes like another gene: the pooled score must lead its
                    // score by `SHORT_MARGIN`
                    if let Some(g) = genome.as_ref() {
                        let refs: Vec<Vec<u8>> = sample.iter().map(|&k| {
                            let (l, c) = (pieces[k].2.min(pieces[k].3.clip.len()), &pieces[k].3.clip);
                            let kept = if plus { &c[..l] } else { &c[c.len() - l..] };
                            let mut v: Vec<u8> = if high { kept.to_vec() } else { kept.iter().rev().copied().collect() };
                            v.truncate(SHORT_QMAX);
                            v
                        }).collect();
                        let partner = index.targets[target as usize].gene;
                        let at = index.short_genome_score(g, chrom, pos, high, plus, &refs, &own, partner, OWN_REACH, o.free);
                        if traced { log::info!("[Place trace] short class at {}:{} ({} end, {} strand): pooled {} on {}, {} on the genome beyond the alignment: {}", chrom, pos, if high { "high" } else { "low" }, if plus { "+" } else { "-" }, pooled, genes[partner as usize], at, if pooled - at >= SHORT_MARGIN { "placed" } else { "the aligned gene's own sequence, no call" }); }
                        if pooled - at < SHORT_MARGIN { continue; }
                    }
                    N_SHORT_CLASS.fetch_add(1, Relaxed);
                    // every piece of the class on the target chosen
                    let tseq = &index.targets[target as usize].seq;
                    let alns: Vec<(i32, usize, usize, usize, usize)> = pieces.iter().map(|p| extend(if from_end { &p.1 } else { &p.0 }, tseq, o.free, o.free + 8)).collect();
                    for ((q, _, l, e), a) in pieces.iter().map(|p| (if from_end { &p.1 } else { &p.0 }, &p.1, p.2, p.3)).zip(&alns) {
                        if (a.0 as f64) < SHORT_FIT_PER_BASE * q.len() as f64 || e.cb - e.ca >= o.min_len { continue; }
                        // the placed bases on the record: piece index -> clip index (in transcript order) -> record
                        let tidx = |x: usize| if from_end == three { l - 1 - x } else { x };
                        let rp = |j: usize| if plus { e.ca + j } else { e.cb - 1 - j };
                        let (p1, p2) = (rp(tidx(a.1)), rp(tidx(a.2.saturating_sub(1).max(a.1))));
                        let (ra, rb) = (p1.min(p2), p1.max(p2) + 1);
                        short_hits.entry((e.name, e.rec_pos)).or_insert(ShortHit { target, aln: *a, ra, rb, plus, pieces: n as u32, excess: margin - SHORT_MARGIN + o.margin });
                    }
                }
            }
        }
        // the reads of the placed short classes, read again: each gets its placed part, and its alignment when it has
        // no other part
        if !short_hits.is_empty() && reader.fetch(chrom.as_str()).is_ok() {
            while let Some(Ok(())) = reader.read(&mut rec) {
                if rec.is_unmapped() || rec.is_secondary() || rec.is_supplementary() { continue; }
                let key = (name_hash(rec.qname()), rec.pos());
                let Some(h) = short_hits.remove(&key) else { continue };
                let len = rec.seq().len();
                let rev = rec.is_reverse();
                let t = &index.targets[h.target as usize];
                let (score, _, _, t0, t1) = h.aln;
                let blocks = t.cover(t0, t1);
                if blocks.is_empty() || h.rb > len { continue; }
                let name = String::from_utf8_lossy(rec.qname()).to_string();
                let sense = rev != h.plus;
                let (qs, qe) = if rev { ((len - h.rb) as i32, (len - h.ra) as i32) } else { (h.ra as i32, h.rb as i32) };
                let mut introns = Vec::new();
                for x in blocks.windows(2) { introns.push(x[0].1 + 1); introns.push(x[1].0 - 1); }
                if trace.as_deref().is_some_and(|tr| name.starts_with(tr)) {
                    log::info!("[Place trace] {} short part of {} bases placed with its class ({} pieces) on {} {}:{:?} score {}", name, h.rb - h.ra, h.pieces,
                               genes[t.gene as usize], index.chroms[t.chrom as usize], blocks, score);
                }
                out.push(ChimericSegment {
                    read_name: name.clone(), dataset: ds.clone(), chromosome: index.chroms[t.chrom as usize].clone(), reverse: sense != t.plus,
                    supplementary: true, mapq: 60, ref_start: blocks[0].0, ref_end: blocks[blocks.len() - 1].1, query_start: qs, query_end: qe,
                    read_length: len as i32, introns, placed: true, short_class: true, pieces: h.pieces, excess: h.excess,
                });
                N_SHORT.fetch_add(1, Relaxed);
                let qname = std::str::from_utf8(rec.qname()).unwrap_or("");
                if parts.get(&(ds.as_str(), qname)).map_or(true, |v| v.is_empty()) && !based.contains(&key) {
                    if let Some(base) = ChimericSegment::of_record(&rec, &names, ds) { out.push(base); }
                }
            }
        }
        searched.sort_unstable();
        for (k, high, end) in placed_at {
            let n = searched.partition_point(|&x| x <= (high, end + w)) - searched.partition_point(|&x| x < (high, end - w));
            out[k].pieces = n as u32;
        }
        (out, replaced)
    }).collect();
    let (mut out, mut replaced) = (Vec::new(), Vec::new());
    for (o, r) in found { out.extend(o); replaced.extend(r); }
    log::info!("[Fusion] Placing: {} pieces of recurring words (adapters) left out", N_ARTEFACT.load(Relaxed));
    log::info!("[Fusion] Placing: {} pieces read, {} searches at exon boundaries and {} inside exons ({} seed hits, {} extensions), {} own-locus alignments, {:.0} s (thread seconds: searches {:.0}, of which seeding {:.0} and extensions {:.0}, {} flank tests; own locus {:.0})",
               N_PIECES.load(Relaxed), N_PLACE.load(Relaxed), N_LOCATE.load(Relaxed), N_HITS.load(Relaxed), N_EXTEND.load(Relaxed), N_OWN.load(Relaxed), started.elapsed().as_secs_f64(),
               T_CHAIN.load(Relaxed) as f64 / 1e9, T_SEED.load(Relaxed) as f64 / 1e9, T_EXT.load(Relaxed) as f64 / 1e9, N_FLANK.load(Relaxed), T_OWN.load(Relaxed) as f64 / 1e9);
    log::info!("[Fusion] {} unaligned read parts placed on another gene's annotated exons ({} of them parts below the mapping-quality floor named by their class, {} short parts placed with their class in {} classes), {} alignments of split reads extended over their gene's own exons",
               out.iter().filter(|s| s.placed).count(), N_CLASS.load(Relaxed), N_SHORT.load(Relaxed), N_SHORT_CLASS.load(Relaxed), replaced.len());
    (out, replaced)
}

#[cfg(test)]
mod tests {
    use super::*;

    const EXON: &[u8] = b"ACGTTGCAAGGCTTAACGGATCCATGCATGCAAGTCCGATCGATCGGCTAGCTAGGATCCGGTTAACC";
    const OTHER: &[u8] = b"TTGACCGGTATAGCGCGATTACGGCATCGATTTAGCCGGATATCGCGCTAATGCGCATTAGCGGCTA";

    fn target(gene: u32, plus: bool, forward: bool, blocks: Vec<(i64, i64)>, seq: &[u8]) -> Target {
        Target { gene, chrom: 0, plus, forward, blocks, seq: seq.to_vec(), tx_end: false }
    }

    fn store(exons: Vec<(u32, bool, i64, &[u8])>) -> ExonStore {
        let mut s = ExonStore { exons: Vec::new(), seq: Vec::new(), offsets: Vec::new(), entries: Vec::new(), plus: FxHashMap::default(), trees: FxHashMap::default() };
        for (gene, plus, lo, x) in exons {
            s.exons.push((gene, 0, lo, lo + x.len() as i64 - 1, s.seq.len() as u32));
            s.seq.extend_from_slice(x); s.seq.push(b'N');
            s.plus.insert(gene, plus);
        }
        let mut counts = vec![0u32; (1usize << (2 * KF)) + 1];
        for (_, k) in kmers_of(&s.seq, KF) { counts[k as usize + 1] += 1; }
        for i in 1..counts.len() { counts[i] += counts[i - 1]; }
        let mut entries = vec![0u32; *counts.last().unwrap() as usize];
        let mut fill = counts.clone();
        for (pos, k) in kmers_of(&s.seq, KF) { entries[fill[k as usize] as usize] = pos as u32; fill[k as usize] += 1; }
        s.offsets = counts; s.entries = entries;
        s
    }

    fn index(targets: Vec<Target>) -> ExonIndex {
        let (forward, backward) = (KmerIndex::build(&targets, true), KmerIndex::build(&targets, false));
        ExonIndex { chroms: vec!["chr1".into()], targets, forward, backward, bounds: FxHashMap::default(), store: store(vec![]) }
    }

    #[test]
    fn a_piece_extends_from_the_boundary_and_stops_where_the_exon_does() {
        // 40 bases of the exon with one substitution, then unrelated sequence
        let mut q = EXON[..40].to_vec();
        q[17] = if q[17] == b'A' { b'C' } else { b'A' };
        q.extend_from_slice(b"TTTTTTTTTTTTTTTTTTTT");
        let (score, q0, q1, t0, t1) = extend(&q, EXON, 6, 14);
        assert_eq!((q0, t0), (0, 0));
        assert_eq!((q1, t1), (40, 40));
        assert_eq!(score, 39 * 2 - 3);
        // the alignment stopped three bases short of the boundary: the piece's first three bases are skipped
        let mut q = b"GGG".to_vec();
        q.extend_from_slice(&EXON[..30]);
        let (score, q0, q1, t0, _) = extend(&q, EXON, 6, 14);
        assert_eq!((score, q0, q1, t0), (60, 3, 33, 0));
        // it ran three bases past it: the target's first three are
        let (score, q0, _, t0, t1) = extend(&EXON[3..33], EXON, 6, 14);
        assert_eq!((score, q0, t0, t1), (60, 0, 3, 33));
    }

    #[test]
    fn target_bases_map_to_genomic_blocks_in_either_direction() {
        // read 3'-ward on the plus strand: up the genome, across an intron
        let t = target(0, true, true, vec![(101, 110), (201, 230)], b"");
        assert_eq!(t.cover(0, 15), vec![(101, 110), (201, 205)]);
        assert_eq!(t.cover(3, 8), vec![(104, 108)]);
        // read 3'-ward on the minus strand: down the genome from the exon's high end
        let t = target(0, false, true, vec![(291, 300), (150, 180)], b"");
        assert_eq!(t.cover(0, 15), vec![(176, 180), (291, 300)]);
        // read 5'-ward on the plus strand: down from the exon's end
        let t = target(0, true, false, vec![(291, 300), (150, 180)], b"");
        assert_eq!(t.cover(2, 12), vec![(179, 180), (291, 298)]);
    }

    #[test]
    fn a_piece_is_placed_on_one_gene_s_exon_only() {
        let o = PlaceOptions::default();
        let mut hits = Vec::new();
        let ix = index(vec![target(1, true, true, vec![(1001, 1067)], EXON), target(2, true, true, vec![(5001, 5067)], OTHER)]);
        let mut q = EXON[..45].to_vec();
        q[20] = if q[20] == b'G' { b'T' } else { b'G' };
        let (p, own) = ix.place(&q, true, &[], None, &o, &mut hits).unwrap();
        assert_eq!((p.target, p.q0, p.q1, p.t0, p.t1, own), (0, 0, 45, 0, 45, false));
        // the aligned gene's own next exon is not another gene
        assert!(ix.place(&q, true, &[1], None, &o, &mut hits).unwrap().1);
        // read 5'-ward there is no such target
        assert!(ix.place(&q, false, &[], None, &o, &mut hits).is_err());
        // two genes sharing the exon's sequence at different places: not placed
        let ix = index(vec![target(1, true, true, vec![(1001, 1067)], EXON), target(2, true, true, vec![(5001, 5067)], EXON)]);
        assert_eq!(ix.place(&q, true, &[], None, &o, &mut hits).unwrap_err().0, "two genes");
        // two overlapping genes sharing the exon itself: placed, and the aligned gene's own exon when either is
        let ix = index(vec![target(1, true, true, vec![(1001, 1067)], EXON), target(2, true, true, vec![(1001, 1060)], EXON)]);
        assert!(!ix.place(&q, true, &[], None, &o, &mut hits).unwrap().1);
        assert!(ix.place(&q, true, &[2], None, &o, &mut hits).unwrap().1);
        // too short a piece
        let ix = index(vec![target(1, true, true, vec![(1001, 1067)], EXON)]);
        assert!(ix.place(&EXON[..15], true, &[], None, &o, &mut hits).is_err());
    }

    #[test]
    fn a_piece_runs_through_the_gene_s_own_short_exon_into_another_gene() {
        let o = PlaceOptions::default();
        let mut hits = Vec::new();
        // the aligned gene's next exon is 40 bases, then the read continues into gene 2
        let ix = index(vec![target(1, true, true, vec![(1001, 1040)], &EXON[..40]), target(2, true, true, vec![(5001, 5067)], OTHER)]);
        let mut q = EXON[..40].to_vec();
        q.extend_from_slice(&OTHER[..50]);
        let steps = ix.chain(&q, true, &[1], true, &o, &mut hits);
        assert_eq!(steps.len(), 2);
        assert_eq!((steps[0].gene, steps[0].off, steps[0].own, steps[0].at_end), (1, 0, true, true));
        assert_eq!((steps[1].gene, steps[1].off, steps[1].own), (2, 40, false));
        assert_eq!((steps[1].q0, steps[1].q1, steps[1].blocks.clone()), (0, 50, vec![(5001, 5050)]));
    }

    #[test]
    fn a_piece_enters_an_exon_anywhere_on_the_strand_that_continues_the_transcript() {
        let o = PlaceOptions::default();
        let mut hits = Vec::new();
        let revcomp = |x: &[u8]| -> Vec<u8> { x.iter().rev().map(|&b| complement(b)).collect() };
        // gene 1 on the plus strand at 1001..: read 3'-ward, the piece runs up the genome from the exon's 11th base
        let st = store(vec![(1, true, 1001, EXON), (2, false, 5001, OTHER)]);
        let mut q = EXON[10..55].to_vec();
        q[20] = if q[20] == b'G' { b'T' } else { b'G' };
        let (l, own) = st.locate(&q, true, &[], None, &o, &mut hits).unwrap();
        assert_eq!((l.lo, l.hi, l.q0, l.q1, l.at_end, own), (1011, 1055, 0, 45, false, false));
        // read 5'-ward the same piece is the transcript backwards: reversed, it runs down the genome
        let back: Vec<u8> = EXON[10..55].iter().rev().copied().collect();
        let (l, _) = st.locate(&back, false, &[], None, &o, &mut hits).unwrap();
        assert_eq!((l.lo, l.hi, l.at_end), (1011, 1055, false));
        // a plus-strand gene's sequence does not continue a transcript read the other way
        assert!(st.locate(&q, false, &[], None, &o, &mut hits).is_err());
        // gene 2 on the minus strand: its transcript is the reverse complement; read 3'-ward the piece runs down the
        // genome, here to the exon's low end
        let t2 = revcomp(OTHER);
        let (l, _) = st.locate(&t2[20..], true, &[], None, &o, &mut hits).unwrap();
        assert_eq!((l.lo, l.hi, l.at_end), (5001, 5001 + OTHER.len() as i64 - 21, true));
        // the aligned gene's own exon is told apart
        assert!(st.locate(&q, true, &[1], None, &o, &mut hits).unwrap().1);
    }

    #[test]
    fn a_gene_s_own_exons_extend_the_alignment() {
        let seg = |s: i64, e: i64, introns: Vec<i64>| ChimericSegment {
            read_name: "r".into(), dataset: "d".into(), chromosome: "chr1".into(), reverse: false, supplementary: false, mapq: 60,
            ref_start: s, ref_end: e, query_start: 100, query_end: 400, read_length: 600, introns, placed: false, short_class: false, pieces: 0, excess: 0 };
        // beyond the high end: the alignment stopped two bases short of the boundary at 2000
        let mut b = seg(1000, 1998, vec![1201, 1499]);
        assert!(extend_segment(&mut b, true, 2000, &[(3001, 3040)], (400, 440)));
        assert_eq!((b.ref_start, b.ref_end, b.query_start, b.query_end), (1000, 3040, 100, 440));
        assert_eq!(b.introns, vec![1201, 1499, 2001, 3000]);
        // beyond the low end
        let mut b = seg(1000, 1998, vec![1201, 1499]);
        assert!(extend_segment(&mut b, false, 1000, &[(501, 540), (701, 730)], (30, 100)));
        assert_eq!((b.ref_start, b.ref_end, b.query_start), (501, 1998, 30));
        assert_eq!(b.introns, vec![541, 700, 731, 999, 1201, 1499]);
        // an exon inside the alignment is not beyond its end
        let mut b = seg(1000, 1998, vec![]);
        assert!(!extend_segment(&mut b, true, 1998, &[(1500, 1540)], (400, 440)));
    }

    #[test]
    fn a_read_whose_bases_before_the_part_are_the_part_s_gene_s_is_that_gene_s() {
        let o = PlaceOptions::default();
        let mut hits = Vec::new();
        let mut ix = index(vec![]);
        ix.store = store(vec![(1, true, 1001, EXON)]);
        // a part placed from the exon's 41st base 3'-ward; the 40 bases before it, read back, are the same exon's
        let before: Vec<u8> = EXON[..40].iter().rev().copied().collect();
        assert!(ix.continues_back(&before, true, 1, &o, &mut hits));
        // another gene's bases are not, nor the same bases asked of another gene
        let other: Vec<u8> = OTHER[..40].iter().rev().copied().collect();
        assert!(!ix.continues_back(&other, true, 1, &o, &mut hits));
        assert!(!ix.continues_back(&before, true, 2, &o, &mut hits));
    }

    #[test]
    fn a_piece_is_found_anywhere_in_a_stretch_of_genome_without_seeds() {
        let mut win = OTHER.to_vec();
        win.extend_from_slice(EXON);
        win.extend_from_slice(OTHER);
        // every ninth base wrong: no seed of nine bases, and still the piece's place
        let mut q = EXON[5..41].to_vec();
        for i in (4..q.len()).step_by(9) { q[i] = complement(q[i]); }
        assert!(own_locus_score(&q, &win, 3).0 >= 36 * 2 - 4 * 5);
        assert!(own_locus_score(&q, OTHER, 3).0 < 28);
        // a tail and a simple repeat are of low complexity, an exon is not
        assert!(low_complexity(b"TTTTTTTTTTATTTTTTGTTTTTTTTTTCTTTTTTTT"));
        assert!(low_complexity(b"CTTTCTTTCTTTCTTTCTTCTTTCTTTCTCTTTCTTTCTTT"));
        assert!(!low_complexity(EXON) && !low_complexity(OTHER));
        assert!(!low_complexity(b"ACGGGCGTGGGCGGCCGGAGCGGCGGCGGCGGCGGCGGGC"));
        assert!(low_complexity(b"TATCAACCAAAAAAAAAAAAAAAATAAAAAAAAATTAAAA"));
        assert!(low_complexity(b"GCATTGACCGTAAAAAAAAAAAAAGCATTGCAGGCTAGCATC"));
        // an adapter is in the pieces of read after read, on either strand; an exon in a few
        let mut art = Artefacts::default();
        let mut scratch = Vec::new();
        let rc: Vec<u8> = OTHER.iter().rev().map(|&b| complement(b)).collect();
        for i in 0..400 { art.add(if i % 2 == 0 { &OTHER[..50] } else { &rc[rc.len() - 50..] }, &mut scratch); }
        for _ in 0..3 { art.add(&EXON[..50], &mut scratch); }
        art.finish();
        assert!(art.holds(&OTHER[5..45]) && art.holds(&rc[20..60]) && !art.holds(EXON));
        // a piece ends where its tail starts
        assert_eq!(before_tail(b"GATTGACACGAGATAAGCGTGGTTTTTTTTTTGATTGGGT"), 22);
        assert_eq!(before_tail(EXON), EXON.len());
        assert!(before_tail(b"GCATCGGATCCATTTATTTTCTTTTTATTTTTTATTTCTTTTTAGC") <= 12);
    }

    #[test]
    fn a_piece_is_read_outward_along_the_gene() {
        let piece = b"AACCG";
        // beyond the high end of a plus gene's exon: as it is, 3'-ward
        assert_eq!(outward(piece, true, true), (b"AACCG".to_vec(), true));
        // beyond the low end of a plus gene's exon: backwards from the alignment, 5'-ward
        assert_eq!(outward(piece, false, true), (b"GCCAA".to_vec(), false));
        // a minus gene's transcript runs down the genome: beyond the low end is 3'-ward, on the other strand
        assert_eq!(outward(piece, false, false), (b"CGGTT".to_vec(), true));
        assert_eq!(outward(piece, true, false), (b"TTGGC".to_vec(), false));
    }

    #[test]
    fn short_pieces_are_cut_at_their_tail_and_judged_by_complexity() {
        // a run of seven A or more, or three A, one other base and four A, is where the tail starts
        assert_eq!(short_tail_cut(b"ACGTAAAAAAAAACC"), 4);
        assert_eq!(short_tail_cut(b"CCTAAAGAAAAT"), 3);
        // shorter runs are sequence: nothing is cut
        assert_eq!(short_tail_cut(b"CCAAGAAAATT"), 11);
        assert_eq!(short_tail_cut(b"GAAAGTAAAC"), 10);
        // three in four bases one base, two distinct bases, or under 8 bases: low complexity
        assert!(short_low_complexity(b"AAAAAAAC"));
        assert!(short_low_complexity(b"ACACACACAC"));
        assert!(short_low_complexity(b"ACGTAC"));
        assert!(!short_low_complexity(b"ACGTTGCAAG"));
        // a tandem repeat of three bases: a GGC repeat after a run of G, a CAG repeat
        assert!(short_low_complexity(b"GGGGGCGGCGGCGGCGGCA"));
        assert!(short_low_complexity(b"CAGCAGCAGCAGCAG"));
        // pieces of 14 and 19 bases without a repeat are sequence
        assert!(!short_low_complexity(b"GCCCTTGATGAGCC"));
        assert!(!short_low_complexity(b"GGTGGACGCTACACAGTGG"));
    }

    #[test]
    fn a_read_s_tail_gives_its_strand() {
        let aligned = b"CGTAGCTAGGCTAGCATCGATCGGACTAGCATCGACT";
        let with = |lead: &[u8], trail: &[u8]| -> (Vec<u8>, usize, usize) {
            ([lead, &aligned[..], trail].concat(), lead.len(), trail.len())
        };
        // A's after the high end: a plus-strand transcript; T's before the low end: minus
        let (s, l, t) = with(b"ACTGG", b"AAAAAAAAAAAAAAAAGTACTCTGCG");
        assert_eq!(tail_strand(&s, l, t), 1);
        let (s, l, t) = with(b"CGCAGAGTACTTTTTTTTTTTTTTTT", b"");
        assert_eq!(tail_strand(&s, l, t), -1);
        // two other bases in the twelve leave a tail; neither end, or both, say nothing
        let (s, l, t) = with(b"", b"AAAAGAAAACAAAAAA");
        assert_eq!(tail_strand(&s, l, t), 1);
        let (s, l, t) = with(b"GCATGCATGCATGCAT", b"CGATCGATCGATCGAT");
        assert_eq!(tail_strand(&s, l, t), 0);
        let (s, l, t) = with(b"TTTTTTTTTTTTTTTT", b"AAAAAAAAAAAAAAAA");
        assert_eq!(tail_strand(&s, l, t), 0);
    }

    #[test]
    fn canonical_words_are_the_same_on_both_strands() {
        let s = b"ACGTTGCAAGTCCATG";
        let rc: Vec<u8> = s.iter().rev().map(|&b| complement(b)).collect();
        let (mut a, mut b) = (canon10(s), canon10(&rc));
        a.sort_unstable(); b.sort_unstable();
        assert_eq!(a, b);
        assert!(canon10(b"ACGTN").is_empty());
    }
}
