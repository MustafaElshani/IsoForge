// Fusion reads the aligner did not split.
//
// - Each candidate junction's sequence is taken from the genome: the last bases of the 5' part and the first bases of
//   the 3' part, in transcript orientation.
// - Primary unsplit reads whose alignment ends within `window` bp of a breakpoint with a soft clip of `min_clip` bases
//   or more on the junction side are read from the BAMs.
// - A read supports the junction when its clip's first, or last, bases align to the other part's sequence at the
//   junction with at most `max_mismatch` of them edited. A clip that is mostly A or T supports none.

use std::path::{Path, PathBuf};

use rayon::prelude::*;
use rust_htslib::bam::{self, record::Aux, record::Cigar, Read};

use crate::genome_ref::GenomeRef;

#[derive(Debug, Clone)]
pub struct RescueOptions {
    /// (dataset, BAM) pairs, the datasets named as in the chimeric segments.
    pub bams: Vec<(String, PathBuf)>,
    pub min_clip: usize,
    /// Bases of each part's junction sequence compared.
    pub seq_len: usize,
    pub max_mismatch: f64,
    pub window: i64,
    pub min_mapq: u8,
    /// A breakpoint shared by more seeds than this is not looked at.
    pub max_targets: usize,
    /// Reads rescued per seed at most, per BAM chromosome.
    pub max_per_seed: usize,
}

impl Default for RescueOptions {
    fn default() -> Self {
        Self { bams: Vec::new(), min_clip: 30, seq_len: 100, max_mismatch: 0.2, window: 10, min_mapq: 20, max_targets: 20, max_per_seed: 50 }
    }
}

/// A candidate junction: its breakpoints (1-based, the 5' part's last base and the 3' part's first), the parts'
/// chromosomes, strands and exon blocks (1-based inclusive, genomic order).
#[derive(Debug, Clone)]
pub struct Seed {
    pub chrom5: String,
    pub chrom3: String,
    pub plus5: bool,
    pub plus3: bool,
    pub bp5: i64,
    pub bp3: i64,
    pub blocks5: Vec<(i64, i64)>,
    pub blocks3: Vec<(i64, i64)>,
    /// A seed of the second rank is looked at after the others and does not count towards their `max_targets`.
    pub second: bool,
}

/// An unsplit read supporting seed `seed`: aligned to the 5' partner (`on5`) or the 3' one, with its aligned blocks.
#[derive(Debug, Clone)]
pub struct Rescued {
    pub seed: usize,
    pub dataset: String,
    pub name: String,
    pub on5: bool,
    pub blocks: Vec<(i64, i64)>,
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter().rev().map(|&b| match b { b'A' => b'T', b'C' => b'G', b'G' => b'C', b'T' => b'A', x => x }).collect()
}

/// Up to `n` bases of a part adjacent to its junction end, as genomic sequence in ascending order: from the high end
/// of `blocks` downward (`junction_high`), or from the low end upward.
fn part_bases(g: &GenomeRef, chrom: &str, blocks: &[(i64, i64)], junction_high: bool, n: usize) -> Option<Vec<u8>> {
    let mut pieces: Vec<Vec<u8>> = Vec::new();
    let mut have = 0usize;
    let order: Vec<&(i64, i64)> = if junction_high { blocks.iter().rev().collect() } else { blocks.iter().collect() };
    for &&(s, e) in &order {
        if have >= n { break; }
        let len = (e - s + 1).max(0) as usize;
        let take = len.min(n - have) as i64;
        let (a, b) = if junction_high { (e - take + 1, e) } else { (s, s + take - 1) };
        let seq = g.fetch_seq(chrom, a - 1, b).ok()?;
        have += seq.len();
        pieces.push(seq);
    }
    if junction_high { pieces.reverse(); }
    Some(pieces.concat())
}

/// The seed's junction sequence in transcript orientation: the 5' part's last bases and the 3' part's first.
fn junction_sides(g: &GenomeRef, s: &Seed, n: usize) -> Option<(Vec<u8>, Vec<u8>)> {
    // the 5' part ends at its high genomic end on +, its low end on -; the 3' part starts at its low end on +
    let g5 = part_bases(g, &s.chrom5, &s.blocks5, s.plus5, n)?;
    let g3 = part_bases(g, &s.chrom3, &s.blocks3, !s.plus3, n)?;
    let tail5 = if s.plus5 { g5 } else { revcomp(&g5) };
    let head3 = if s.plus3 { g3 } else { revcomp(&g3) };
    Some((tail5, head3))
}

/// Edit distance of `q` aligned whole against `t`, starting within the first `free` bases of `t` and ending anywhere.
fn fit_distance(q: &[u8], t: &[u8], free: usize) -> usize {
    let (m, n) = (q.len(), t.len());
    if m == 0 { return 0; }
    // starting past the first `free` bases costs the bases skipped
    let mut prev: Vec<usize> = (0..=n).map(|j| j.saturating_sub(free)).collect();
    let mut cur = vec![0usize; n + 1];
    for i in 1..=m {
        cur[0] = i;
        for j in 1..=n {
            let sub = prev[j - 1] + usize::from(q[i - 1] != t[j - 1]);
            cur[j] = sub.min(prev[j] + 1).min(cur[j - 1] + 1);
        }
        std::mem::swap(&mut prev, &mut cur);
    }
    *prev.iter().min().unwrap()
}

/// Whether a clip (transcript orientation) continues into the other part: the clip's first bases against the 3'
/// part's head (`forward`), or its last bases against the 5' part's tail.
fn clip_matches(clip: &[u8], other: &[u8], forward: bool, o: &RescueOptions) -> bool {
    let l = clip.len().min(o.seq_len).min(other.len());
    if l < o.min_clip { return false; }
    let at = clip.iter().filter(|&&b| b == b'A').count().max(clip.iter().filter(|&&b| b == b'T').count());
    if at as f64 > 0.7 * clip.len() as f64 { return false; }
    let w = o.window as usize;
    let (q, t): (Vec<u8>, Vec<u8>) = if forward {
        (clip[..l].to_vec(), other[..(l + 2 * w).min(other.len())].to_vec())
    } else {
        let q: Vec<u8> = clip[clip.len() - l..].iter().rev().copied().collect();
        let t: Vec<u8> = other[other.len().saturating_sub(l + 2 * w)..].iter().rev().copied().collect();
        (q, t)
    };
    fit_distance(&q, &t, 2 * w) as f64 <= o.max_mismatch * l as f64
}

/// Aligned blocks (1-based inclusive) and soft clips (leading, trailing) of a record.
fn alignment(record: &bam::Record) -> (Vec<(i64, i64)>, usize, usize) {
    let mut blocks = Vec::new();
    let mut pos = record.pos() + 1;
    let mut start = pos;
    let (mut lead, mut trail, mut seen) = (0usize, 0usize, false);
    for op in record.cigar().iter() {
        match *op {
            Cigar::SoftClip(n) => { if seen { trail += n as usize } else { lead += n as usize } }
            Cigar::HardClip(_) | Cigar::Pad(_) | Cigar::Ins(_) => {}
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) | Cigar::Del(n) => { seen = true; pos += n as i64; }
            Cigar::RefSkip(n) => { blocks.push((start, pos - 1)); pos += n as i64; start = pos; }
        }
    }
    blocks.push((start, pos - 1));
    (blocks, lead, trail)
}

/// One side of a seed to look for (indexed by chromosome, alignment end and breakpoint): a read's clip there must
/// continue into `other` (`forward`: the clip's first bases against `other`'s first; else its last against
/// `other`'s last), the clip read in transcript orientation after reverse-complementing when `minus`.
struct Target {
    seed: usize,
    second: bool,
    on5: bool,
    forward: bool,
    minus: bool,
    other: std::sync::Arc<Vec<u8>>,
    /// 12-mers of `other` near the junction, tested before the alignment.
    kmers: std::sync::Arc<rustc_hash::FxHashSet<u32>>,
}

const K: usize = 12;

fn kmer_codes(s: &[u8]) -> impl Iterator<Item = u32> + '_ {
    s.windows(K).filter_map(|w| {
        let mut x = 0u32;
        for &b in w {
            x = (x << 2) | match b { b'A' => 0, b'C' => 1, b'G' => 2, b'T' => 3, _ => return None };
        }
        Some(x)
    })
}

/// The unsplit reads supporting each seed, from every BAM: one pass over each chromosome of each BAM, every primary
/// unsplit record's clipped ends looked up among the seeds' breakpoints; a read supports at most one seed, and a
/// breakpoint shared by more than `max_targets` seeds is not looked at.
pub fn rescue(seeds: &[Seed], o: &RescueOptions, genome: &Path) -> Vec<Rescued> {
    use rustc_hash::{FxHashMap, FxHashSet};
    use std::sync::Arc;
    if o.bams.is_empty() || seeds.is_empty() { return Vec::new(); }
    let sides: Vec<Option<(Vec<u8>, Vec<u8>)>> = seeds.par_iter().map_init(
        || GenomeRef::open(genome.to_path_buf()).ok(),
        |g, s| g.as_ref().and_then(|g| junction_sides(g, s, o.seq_len + 2 * o.window as usize))).collect();
    let kmer_set = |seq: &[u8], forward: bool| -> FxHashSet<u32> {
        let n = (o.seq_len / 2 + 2 * o.window as usize).min(seq.len());
        let part = if forward { &seq[..n] } else { &seq[seq.len() - n..] };
        kmer_codes(part).collect()
    };
    // chromosome -> (at_high, position) -> targets
    let mut index: FxHashMap<String, FxHashMap<(bool, i64), Vec<Target>>> = FxHashMap::default();
    for (si, (s, sd)) in seeds.iter().zip(sides).enumerate() {
        let Some((tail5, head3)) = sd else { continue };
        let (k5, k3) = (Arc::new(kmer_set(&tail5, false)), Arc::new(kmer_set(&head3, true)));
        let (tail5, head3) = (Arc::new(tail5), Arc::new(head3));
        // reads on the 5' part: the junction at its high end on +, its low end on -; the clip continues into the 3'
        // part's head. Reads on the 3' part: the junction at its low end on +, its high end on -; the clip leads out
        // of the 5' part's tail.
        index.entry(s.chrom5.clone()).or_default().entry((s.plus5, s.bp5)).or_default().push(Target {
            seed: si, second: s.second, on5: true, forward: true, minus: !s.plus5, other: head3.clone(), kmers: k3.clone() });
        index.entry(s.chrom3.clone()).or_default().entry((!s.plus3, s.bp3)).or_default().push(Target {
            seed: si, second: s.second, on5: false, forward: false, minus: !s.plus3, other: tail5.clone(), kmers: k5.clone() });
    }
    for m in index.values_mut() {
        m.retain(|_, v| {
            if v.iter().filter(|t| !t.second).count() > o.max_targets { return false; }
            // over the limit only with the seeds of the second rank: those are dropped
            if v.len() > o.max_targets { v.retain(|t| !t.second); }
            !v.is_empty()
        });
    }
    let mut jobs: Vec<(usize, String)> = Vec::new();
    for (bi, (_, p)) in o.bams.iter().enumerate() {
        let Ok(r) = bam::IndexedReader::from_path(p) else { continue };
        for name in r.header().target_names() {
            let c = String::from_utf8_lossy(name).to_string();
            if index.get(&c).is_some_and(|m| !m.is_empty()) { jobs.push((bi, c)); }
        }
    }
    let found: Vec<Vec<Rescued>> = jobs.par_iter().map(|(bi, chrom)| {
        let mut out = Vec::new();
        let (ds, path) = &o.bams[*bi];
        let at = &index[chrom];
        let Ok(mut reader) = bam::IndexedReader::from_path(path) else { return out };
        if reader.fetch(chrom.as_str()).is_err() { return out; }
        let mut per_seed: FxHashMap<usize, usize> = FxHashMap::default();
        let mut rec = bam::Record::new();
        while let Some(Ok(())) = reader.read(&mut rec) {
            if rec.is_unmapped() || rec.is_secondary() || rec.is_supplementary() || rec.mapq() < o.min_mapq { continue; }
            let cig = rec.cigar();
            let lead = match cig.first() { Some(Cigar::SoftClip(n)) => *n as usize, _ => 0 };
            let trail = match cig.last() { Some(Cigar::SoftClip(n)) => *n as usize, _ => 0 };
            if lead < o.min_clip && trail < o.min_clip { continue; }
            let first = rec.pos() + 1;
            let last = cig.end_pos();
            let near = |high: bool, end: i64| (-o.window..=o.window).any(|d| at.contains_key(&(high, end + d)));
            let (try_high, try_low) = (trail >= o.min_clip && near(true, last), lead >= o.min_clip && near(false, first));
            if !try_high && !try_low { continue; }
            if matches!(rec.aux(b"SA"), Ok(Aux::String(_))) { continue; }
            let seq = rec.seq().as_bytes();
            let mut done = false;
            for (high, end, clip, go) in [(true, last, trail, try_high), (false, first, lead, try_low)] {
                if done || !go { continue; }
                let piece: Vec<u8> = (if high { &seq[seq.len() - clip..] } else { &seq[..clip] }).iter().map(|b| b.to_ascii_uppercase()).collect();
                let rc = revcomp(&piece);
                // the seeds of the first rank at every breakpoint within the window, then those of the second
                for (rank, d) in [false, true].into_iter().flat_map(|r| (-o.window..=o.window).map(move |d| (r, d))) {
                    let Some(ts) = at.get(&(high, end + d)) else { continue };
                    for t in ts.iter().filter(|t| t.second == rank) {
                        if per_seed.get(&t.seed).is_some_and(|&n| n >= o.max_per_seed) { continue; }
                        let clip_t: &[u8] = if t.minus { &rc } else { &piece };
                        let n = (o.seq_len / 2).min(clip_t.len());
                        let probe = if t.forward { &clip_t[..n] } else { &clip_t[clip_t.len() - n..] };
                        if kmer_codes(probe).filter(|k| t.kmers.contains(k)).take(2).count() < 2 { continue; }
                        if clip_matches(clip_t, &t.other, t.forward, o) {
                            let (blocks, _, _) = alignment(&rec);
                            out.push(Rescued { seed: t.seed, dataset: ds.clone(), name: String::from_utf8_lossy(rec.qname()).to_string(),
                                               on5: t.on5, blocks });
                            *per_seed.entry(t.seed).or_insert(0) += 1;
                            done = true;
                            break;
                        }
                    }
                    if done { break; }
                }
            }
        }
        out
    }).collect();
    let mut seen = FxHashSet::default();
    found.into_iter().flatten().filter(|r| seen.insert((r.dataset.clone(), r.name.clone()))).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_clip_is_fitted_against_the_other_part_with_a_free_start() {
        let other = b"ACGTTGCAAGGCTTAACGGATCCATGCATGCAAGTCCGATCGATCGGCTAGCTAGGATCC".to_vec();
        let o = RescueOptions { min_clip: 20, seq_len: 40, window: 5, ..Default::default() };
        // the clip starts three bases into the other part (the alignment ran past the breakpoint), one substitution
        let mut clip = other[3..43].to_vec();
        clip[10] = if clip[10] == b'A' { b'C' } else { b'A' };
        assert!(clip_matches(&clip, &other, true, &o));
        // an unrelated clip does not fit
        let junk = b"GGGGCCCCTTTTAAAAGGGGCCCCTTTTAAAAGGGGCCCC".to_vec();
        assert!(!clip_matches(&junk, &other, true, &o));
        // a poly(A) clip is a tail
        assert!(!clip_matches(&vec![b'A'; 40], &vec![b'A'; 60], true, &o));
        // the last bases of a clip against the 5' part's tail
        let tail = other.clone();
        let clip = tail[tail.len() - 35..tail.len() - 2].to_vec();
        assert!(clip_matches(&clip, &tail, false, &o));
    }

    #[test]
    fn edit_distance_allows_a_free_start_only_at_the_beginning() {
        assert_eq!(fit_distance(b"ACGT", b"TTACGT", 2), 0);
        assert_eq!(fit_distance(b"ACGT", b"TTTTACGT", 2), 2);
    }
}
