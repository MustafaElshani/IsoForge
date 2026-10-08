//! Terminal exons an aligner left in a soft clip.
//!
//! When an end of the alignment lies within `BOUNDARY_SLACK` bp of an annotated intron's flanking exon boundary,
//! the read's bases beyond that boundary (the clip, and any bases the aligner ran into the intron) are aligned,
//! anchored at the intron, to the genome beyond it; a match of at least `MIN_RESCUED` exon bases at about 83 %
//! identity adds the intron and the matched exon to the read. Only annotated introns on the read's strand are
//! tried, one exon at each end.
use crate::data_structures::SharedRefs;
use crate::genome_ref::GenomeRef;
use crate::types::{EdgeType, Strand};
use rust_htslib::bam;
use rust_htslib::bam::record::Cigar;
use std::collections::HashMap;
use std::sync::atomic::{AtomicU64, Ordering};

/// A read end this far (bp) from an annotated exon boundary is tried.
pub const BOUNDARY_SLACK: i64 = 5;
/// Bases of the exon beyond the intron a clip must match to be placed there.
pub const MIN_RESCUED: usize = 12;
/// Read bases compared beyond the alignment.
const MAX_OUTWARD: usize = 400;
const MATCH: i32 = 1;
const MISMATCH: i32 = -2;
const GAP: i32 = -2;

/// Reads with an exon rescued at their low and at their high genomic end, over the run.
pub static RESCUED_LOW: AtomicU64 = AtomicU64::new(0);
pub static RESCUED_HIGH: AtomicU64 = AtomicU64::new(0);

#[derive(Default, Debug)]
pub struct TerminalExonIndex {
    /// (chromosome, strand) -> first exon base right of an annotated intron (last intron base + 1) -> the
    /// intron's first bases.
    right_of: HashMap<(String, Strand), HashMap<i64, Vec<i64>>>,
    /// (chromosome, strand) -> last exon base left of an annotated intron (first intron base - 1) -> the
    /// intron's last bases.
    left_of: HashMap<(String, Strand), HashMap<i64, Vec<i64>>>,
}

/// What a rescue changed: clip bases left beyond the alignment at each genomic end, whether an exon was added at
/// each end, and, when an exon was added at the 5' end, whether the bases beyond it carry the cap signature.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Rescued {
    pub leading_clip: u32,
    pub trailing_clip: u32,
    pub low: bool,
    pub high: bool,
    pub capped: Option<bool>,
}

impl TerminalExonIndex {
    pub fn from_reference(refs: &SharedRefs) -> Self {
        let pos: HashMap<_, _> = refs.location_dict.chromosomes.values()
            .flat_map(|c| c.positions.iter().map(|&(p, v)| (v, p))).collect();
        let chrom_of: HashMap<_, _> = refs.location_dict.chromosomes.iter()
            .flat_map(|(chrom, c)| c.positions.iter().map(move |&(_, v)| (v, chrom.clone()))).collect();
        let mut index = Self::default();
        for e in refs.edge_dict.values() {
            if e.edge_type != EdgeType::Intron { continue; }
            let (Some(&p1), Some(&p2), Some(chrom)) = (pos.get(&e.v1), pos.get(&e.v2), chrom_of.get(&e.v1)) else { continue };
            index.insert(chrom, e.strand, p1.min(p2) + 1, p1.max(p2) - 1);
        }
        index
    }

    /// An annotated intron, first and last intron base (1-based).
    pub fn insert(&mut self, chrom: &str, strand: Strand, first: i64, last: i64) {
        if first > last { return; }
        let key = (chrom.to_string(), strand);
        let r = self.right_of.entry(key.clone()).or_default().entry(last + 1).or_default();
        if !r.contains(&first) { r.push(first); }
        let l = self.left_of.entry(key).or_default().entry(first - 1).or_default();
        if !l.contains(&last) { l.push(last); }
    }

    /// Places the read's clipped bases at either end of its alignment on the annotated exon beyond the intron at
    /// that end, when they match it. `introns` (first, last intron base, flat, genomic order), `start` and `end`
    /// (1-based, inclusive) are updated in place.
    pub fn rescue(&self, record: &bam::Record, chrom: &str, strand: Strand, genome: &GenomeRef,
                  introns: &mut Vec<i64>, start: &mut i64, end: &mut i64) -> Rescued {
        let cigar = record.cigar();
        let (lead, trail) = crate::bam_reader::extract_soft_clips(&cigar);
        let mut out = Rescued { leading_clip: lead, trailing_clip: trail, low: false, high: false, capped: None };
        let key = (chrom.to_string(), strand);
        let seq = record.seq().as_bytes();
        let aligned = |op: &Cigar| match *op { Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => Some(n as i64), _ => None };
        let first_block = cigar.iter().find(|op| !matches!(op, Cigar::SoftClip(_) | Cigar::HardClip(_))).and_then(aligned).unwrap_or(0);
        let last_block = cigar.iter().rev().find(|op| !matches!(op, Cigar::SoftClip(_) | Cigar::HardClip(_))).and_then(aligned).unwrap_or(0);
        let five_is_low = strand == Strand::Plus;

        // low genomic end: the exon left of an annotated intron whose right exon starts near `start`
        if let Some(map) = self.right_of.get(&key) {
            let mut best: Option<(i32, i64, i64, usize, i64, i64)> = None; // (score, first, last, exon bases, read used, overhang)
            for b in (*start - BOUNDARY_SLACK)..=(*start + BOUNDARY_SLACK) {
                let Some(firsts) = map.get(&b) else { continue };
                let delta = *start - b; // > 0: the alignment starts inside the right exon
                let over = (-delta).max(0);
                if over > first_block || (lead as i64 + over) < MIN_RESCUED as i64 { continue; }
                let q_end = lead as usize + over as usize;
                let n = q_end.min(MAX_OUTWARD);
                let read_out: Vec<u8> = seq[q_end - n..q_end].iter().rev().copied().collect();
                let head: Vec<u8> = if delta > 0 {
                    match genome.fetch_seq(chrom, b - 1, *start - 1) { Ok(s) => s.into_iter().rev().collect(), Err(_) => continue }
                } else { Vec::new() };
                for &first in firsts {
                    if introns.first().is_some_and(|&f| f <= b) { continue; }
                    let e = first - 1;
                    let k = n as i64 + 20;
                    let Ok(up) = genome.fetch_seq(chrom, (e - k).max(0), e) else { continue };
                    let mut reference = head.clone();
                    reference.extend(up.iter().rev());
                    let (score, read_used, ref_used) = anchored(&read_out, &reference);
                    let exon = ref_used.saturating_sub(head.len());
                    if accept(score, ref_used, exon) && best.map_or(true, |x| score > x.0) {
                        best = Some((score, first, b - 1, exon, read_used as i64, over));
                    }
                }
            }
            if let Some((_, first, last, exon, read_used, over)) = best {
                // the rescued exon ends at first - 1; the aligned bases that ran into the intron are part of it
                introns.insert(0, last);
                introns.insert(0, first);
                *start = first - exon as i64;
                out.leading_clip = (lead as i64 - (read_used - over)).max(0) as u32;
                out.low = true;
                RESCUED_LOW.fetch_add(1, Ordering::Relaxed);
                if five_is_low {
                    // the read bases beyond the rescued exon, next to it
                    let rest = out.leading_clip as usize;
                    let cap = rest >= 2 && seq[rest - 2..rest].eq_ignore_ascii_case(b"GG")
                        && genome.fetch_seq(chrom, *start - 3, *start - 1).map_or(true, |g| !g.eq_ignore_ascii_case(b"GG"));
                    out.capped = Some(cap);
                }
            }
        }

        // high genomic end: the exon right of an annotated intron whose left exon ends near `end`
        if let Some(map) = self.left_of.get(&key) {
            let mut best: Option<(i32, i64, i64, usize, i64, i64)> = None; // (score, first, last, exon bases, read used, overhang)
            let len = seq.len();
            for b in (*end - BOUNDARY_SLACK)..=(*end + BOUNDARY_SLACK) {
                let Some(lasts) = map.get(&b) else { continue };
                let delta = b - *end; // > 0: the alignment ends before the left exon's end
                let over = (-delta).max(0);
                if over > last_block || (trail as i64 + over) < MIN_RESCUED as i64 { continue; }
                let q_start = len - trail as usize - over as usize;
                let n = (len - q_start).min(MAX_OUTWARD);
                let read_out: Vec<u8> = seq[q_start..q_start + n].to_vec();
                let head: Vec<u8> = if delta > 0 {
                    match genome.fetch_seq(chrom, *end, b) { Ok(s) => s, Err(_) => continue }
                } else { Vec::new() };
                for &last in lasts {
                    if introns.last().is_some_and(|&l| l >= b) { continue; }
                    let s = last + 1;
                    let k = n as i64 + 20;
                    let Ok(down) = genome.fetch_seq(chrom, s - 1, s - 1 + k) else { continue };
                    let mut reference = head.clone();
                    reference.extend(down.iter());
                    let (score, read_used, ref_used) = anchored(&read_out, &reference);
                    let exon = ref_used.saturating_sub(head.len());
                    if accept(score, ref_used, exon) && best.map_or(true, |x| score > x.0) {
                        best = Some((score, b + 1, last, exon, read_used as i64, over));
                    }
                }
            }
            if let Some((_, first, last, exon, read_used, over)) = best {
                introns.push(first);
                introns.push(last);
                *end = last + exon as i64;
                out.trailing_clip = (trail as i64 - (read_used - over)).max(0) as u32;
                out.high = true;
                RESCUED_HIGH.fetch_add(1, Ordering::Relaxed);
                if !five_is_low {
                    // bases beyond the rescued exon, next to it; a minus-strand cap reads CC on the genome strand
                    let rest = out.trailing_clip as usize;
                    let cap = rest >= 2 && seq[len - rest..len - rest + 2].eq_ignore_ascii_case(b"CC")
                        && genome.fetch_seq(chrom, *end, *end + 2).map_or(true, |g| !g.eq_ignore_ascii_case(b"CC"));
                    out.capped = Some(cap);
                }
            }
        }
        out
    }
}

/// A rescue needs `MIN_RESCUED` bases of the exon beyond the intron and a score of at least half the reference
/// bases aligned (about 83 % identity with mismatches only).
fn accept(score: i32, ref_used: usize, exon: usize) -> bool {
    exon >= MIN_RESCUED && score >= 10 && 2 * score as usize >= ref_used
}

/// Alignment of `read` and `reference` anchored at their first bases, free at their ends: the best score over all
/// end points and the read and reference bases it uses (ties go to the longer reference).
pub fn anchored(read: &[u8], reference: &[u8]) -> (i32, usize, usize) {
    let (n, m) = (read.len(), reference.len());
    let mut prev: Vec<i32> = (0..=m as i32).map(|j| j * GAP).collect();
    let mut cur = vec![0i32; m + 1];
    let mut best = (0i32, 0usize, 0usize);
    for i in 1..=n {
        cur[0] = i as i32 * GAP;
        let a = read[i - 1].to_ascii_uppercase();
        for j in 1..=m {
            let s = if a == reference[j - 1].to_ascii_uppercase() && a != b'N' { MATCH } else { MISMATCH };
            cur[j] = (prev[j - 1] + s).max(prev[j] + GAP).max(cur[j - 1] + GAP);
            if cur[j] > best.0 || (cur[j] == best.0 && j > best.2) { best = (cur[j], i, j); }
        }
        std::mem::swap(&mut prev, &mut cur);
    }
    best
}

#[cfg(test)]
mod tests {
    use super::*;
    use bam::record::CigarString;

    /// A 1,000-bp pseudo-random genome; intron 131-400 (exon 1 is 106-130, 25 bp) and intron 601-800.
    fn fixture() -> (tempfile::TempDir, GenomeRef, TerminalExonIndex) {
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("genome.fa");
        let mut x: u64 = 7;
        let g: Vec<u8> = (0..1000).map(|_| { x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407); b"ACGT"[(x >> 62) as usize] }).collect();
        std::fs::write(&fasta, format!(">chr1\n{}\n", String::from_utf8(g).unwrap())).unwrap();
        std::fs::write(dir.path().join("genome.fa.fai"), "chr1\t1000\t6\t1000\t1001\n").unwrap();
        let genome = GenomeRef::open(fasta).unwrap();
        let mut index = TerminalExonIndex::default();
        index.insert("chr1", Strand::Plus, 131, 400);
        index.insert("chr1", Strand::Plus, 601, 800);
        (dir, genome, index)
    }

    fn record(seq: &[u8], cigar: Vec<Cigar>, pos: i64) -> bam::Record {
        let mut r = bam::Record::new();
        r.set(b"r", Some(&CigarString(cigar)), seq, &vec![30; seq.len()]);
        r.set_pos(pos);
        r
    }

    #[test]
    fn a_short_first_exon_in_the_clip_is_placed_across_its_intron() {
        let (_d, genome, index) = fixture();
        let mut seq = genome.fetch_seq("chr1", 105, 130).unwrap(); // exon 1, 1-based 106-130
        seq.extend(genome.fetch_seq("chr1", 400, 500).unwrap()); // exon 2 from 401
        let r = record(&seq, vec![Cigar::SoftClip(25), Cigar::Match(100)], 400);
        let (mut introns, mut start, mut end) = (vec![], 401, 500);
        let out = index.rescue(&r, "chr1", Strand::Plus, &genome, &mut introns, &mut start, &mut end);
        assert_eq!((introns, start, end, out.leading_clip, out.low, out.high), (vec![131, 400], 106, 500, 0, true, false));
        // the aligner ran 3 bases into the intron: they belong to exon 1
        let r = record(&seq, vec![Cigar::SoftClip(22), Cigar::Match(103)], 397);
        let (mut introns, mut start, mut end) = (vec![], 398, 500);
        let out = index.rescue(&r, "chr1", Strand::Plus, &genome, &mut introns, &mut start, &mut end);
        assert_eq!((introns, start, out.leading_clip), (vec![131, 400], 106, 0));
    }

    #[test]
    fn a_last_exon_before_the_tail_is_placed_and_the_tail_stays_clipped() {
        let (_d, genome, index) = fixture();
        let mut seq = genome.fetch_seq("chr1", 500, 600).unwrap(); // 1-based 501-600
        seq.extend(genome.fetch_seq("chr1", 800, 860).unwrap()); // exon from 801, 60 bp
        seq.extend([b'A'; 20]);
        let r = record(&seq, vec![Cigar::Match(100), Cigar::SoftClip(80)], 500);
        let (mut introns, mut start, mut end) = (vec![], 501, 600);
        let out = index.rescue(&r, "chr1", Strand::Plus, &genome, &mut introns, &mut start, &mut end);
        assert_eq!((introns, end, out.trailing_clip, out.high), (vec![601, 800], 860, 20, true));
    }

    #[test]
    fn a_clip_matching_nothing_is_left_alone() {
        let (_d, genome, index) = fixture();
        let mut seq = b"GATTACAGATTACAGATTACAGATTACA".to_vec();
        seq.extend(genome.fetch_seq("chr1", 400, 500).unwrap());
        let r = record(&seq, vec![Cigar::SoftClip(28), Cigar::Match(100)], 400);
        let (mut introns, mut start, mut end) = (vec![], 401, 500);
        let out = index.rescue(&r, "chr1", Strand::Plus, &genome, &mut introns, &mut start, &mut end);
        assert!(introns.is_empty() && start == 401 && !out.low && out.leading_clip == 28);
    }

    #[test]
    fn anchored_alignment_scores_matches_and_stops_where_they_end() {
        assert_eq!(anchored(b"ACGTACGTAC", b"ACGTACGTACTTTT"), (10, 10, 10));
        let (score, _, used) = anchored(b"ACGTTCGTACGGGGGG", b"ACGTACGTACCCCCCC");
        assert_eq!((score, used), (7, 10));
    }
}
