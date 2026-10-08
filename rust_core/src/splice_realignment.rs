//! Sequence-based realignment of splice sites onto nearby annotated introns.
//!
//! When a read's intron is unannotated but shares one site with an annotated
//! intron whose other site lies within `window` bp, the read's bases across the
//! junction are scored against both spliced genomic templates. The read takes
//! the annotated intron unless its observed intron fits better by at least
//! `MIN_OBSERVED_ADVANTAGE` edits.
use crate::data_structures::SharedRefs;
use crate::genome_ref::GenomeRef;
use crate::microexon_repair::{edit_distance, query_boundary};
use crate::types::{EdgeType, Strand};
use rust_htslib::bam;
use std::collections::{BTreeMap, HashMap, HashSet};

/// Exonic bases scored on each side, beyond the outermost of the two sites.
const FLANK: i64 = 20;
/// Edits by which the observed intron must beat the annotated one to be kept.
const MIN_OBSERVED_ADVANTAGE: usize = 2;
/// Edits by which the read's bases must favour its unannotated intron over every annotated one nearby for its sites
/// to be settled against snapping (`Realigned::kept`).
const SETTLE_ADVANTAGE: usize = 3;
/// Candidate introns considered per read intron; beyond this the locus is left alone.
const MAX_CANDIDATES: usize = 8;

#[derive(Default, Debug)]
pub struct SpliceRealignIndex {
    /// Annotated introns as (first intron base, last intron base), 1-based,
    /// keyed by (chromosome, strand); the same convention as a read's intron list.
    exact: HashMap<(String, Strand), HashSet<(i64, i64)>>,
    by_start: HashMap<(String, Strand), BTreeMap<i64, Vec<i64>>>,
    by_end: HashMap<(String, Strand), BTreeMap<i64, Vec<i64>>>,
    window: i64,
}

impl SpliceRealignIndex {
    pub fn from_reference(refs: &SharedRefs, window: i64) -> Self {
        let pos: HashMap<_, _> = refs
            .location_dict
            .chromosomes
            .values()
            .flat_map(|c| c.positions.iter().map(|&(p, v)| (v, p)))
            .collect();
        let chrom_of: HashMap<_, _> = refs
            .location_dict
            .chromosomes
            .iter()
            .flat_map(|(chrom, c)| c.positions.iter().map(move |&(_, v)| (v, chrom.clone())))
            .collect();
        let mut index = Self { window, ..Self::default() };
        for e in refs.edge_dict.values() {
            if e.edge_type != EdgeType::Intron {
                continue;
            }
            let (Some(&p1), Some(&p2), Some(chrom)) = (pos.get(&e.v1), pos.get(&e.v2), chrom_of.get(&e.v1)) else {
                continue;
            };
            index.insert(chrom, e.strand, p1.min(p2) + 1, p1.max(p2) - 1);
        }
        index
    }

    fn insert(&mut self, chrom: &str, strand: Strand, first: i64, last: i64) {
        if first > last {
            return;
        }
        let key = (chrom.to_string(), strand);
        if self.exact.entry(key.clone()).or_default().insert((first, last)) {
            self.by_start.entry(key.clone()).or_default().entry(first).or_default().push(last);
            self.by_end.entry(key).or_default().entry(last).or_default().push(first);
        }
    }

    /// Annotated introns other than (first, last) that share one of its sites and
    /// whose other site lies within the window.
    fn alternatives(&self, key: &(String, Strand), first: i64, last: i64) -> Vec<(i64, i64)> {
        let mut out = Vec::new();
        if let Some(ends) = self.by_start.get(key).and_then(|m| m.get(&first)) {
            out.extend(ends.iter().filter(|&&l| l != last && (l - last).abs() <= self.window).map(|&l| (first, l)));
        }
        if let Some(starts) = self.by_end.get(key).and_then(|m| m.get(&last)) {
            out.extend(starts.iter().filter(|&&f| f != first && (f - first).abs() <= self.window).map(|&f| (f, last)));
        }
        out
    }

    /// Edit distance of the read's bases across the junction to the observed
    /// intron's spliced template and to each alternative's, sorted best first.
    fn score(
        record: &bam::Record,
        chrom: &str,
        genome: &GenomeRef,
        observed: (i64, i64),
        alternatives: &[(i64, i64)],
    ) -> Option<(usize, Vec<(usize, (i64, i64))>, usize)> {
        // 0-based half-open exon flanks around every intron considered.
        let lo = alternatives.iter().map(|c| c.0).min()?.min(observed.0) - 1 - FLANK;
        let hi = alternatives.iter().map(|c| c.1).max()?.max(observed.1) + FLANK;
        if lo < record.pos() || hi > record.cigar().end_pos() {
            return None;
        }
        let (qlo, qhi) = (query_boundary(record, lo)?, query_boundary(record, hi)?);
        if qhi <= qlo || qhi - qlo > 4 * (hi - lo) as usize {
            return None;
        }
        let sequence = record.seq().as_bytes();
        let query = &sequence[qlo..qhi];
        let template = |f: i64, l: i64| -> Option<Vec<u8>> {
            let mut t = genome.fetch_seq(chrom, lo, f - 1).ok()?;
            t.extend(genome.fetch_seq(chrom, l, hi).ok()?);
            Some(t)
        };
        let observed_score = edit_distance(query, &template(observed.0, observed.1)?);
        let mut scored: Vec<(usize, (i64, i64))> = alternatives
            .iter()
            .filter_map(|&c| template(c.0, c.1).map(|t| (edit_distance(query, &t), c)))
            .collect();
        scored.sort_unstable();
        Some((observed_score, scored, query.len()))
    }

    /// Realigns unannotated introns onto nearby annotated ones and flags annotated
    /// introns whose read sequence cannot tell them from a nearby annotated
    /// alternative. `ambiguous[k]` refers to intron k of the returned list.
    pub fn realign(
        &self,
        record: &bam::Record,
        chrom: &str,
        strand: Strand,
        genome: &GenomeRef,
        introns: &[i64],
    ) -> Realigned {
        let n = introns.len() / 2;
        let unchanged = Realigned { introns: introns.to_vec(), moved: 0, ambiguous: vec![false; n], kept: vec![false; n] };
        let key = (chrom.to_string(), strand);
        let Some(exact) = self.exact.get(&key) else { return unchanged };
        if self.window <= 0 {
            return unchanged;
        }
        let mut output = introns.to_vec();
        let mut ambiguous = vec![false; n];
        let mut kept = vec![false; n];
        let mut moved = 0;
        for k in 0..n {
            let observed = (introns[2 * k], introns[2 * k + 1]);
            let alternatives = self.alternatives(&key, observed.0, observed.1);
            if alternatives.is_empty() || alternatives.len() > MAX_CANDIDATES {
                continue;
            }
            let Some((observed_score, scored, query_len)) = Self::score(record, chrom, genome, observed, &alternatives) else {
                continue;
            };
            let (best, site) = scored[0];
            if exact.contains(&observed) {
                // Two annotated sites a few bases apart: unless the read's bases
                // favour the observed one, the site is left to chain assignment.
                ambiguous[k] = best < observed_score + MIN_OBSERVED_ADVANTAGE;
                continue;
            }
            // A tie between two annotated introns leaves the read as aligned.
            if scored.get(1).is_some_and(|s| s.0 == best) {
                continue;
            }
            if best < observed_score + MIN_OBSERVED_ADVANTAGE && best * 4 <= query_len {
                output[2 * k] = site.0;
                output[2 * k + 1] = site.1;
                moved += 1;
                ambiguous[k] = scored.get(1).is_some_and(|s| s.0 < best + MIN_OBSERVED_ADVANTAGE);
            } else if best >= observed_score + SETTLE_ADVANTAGE {
                // the read's own bases favour its unannotated intron over every annotated one nearby
                kept[k] = true;
            }
        }
        // A realigned intron may not cross its neighbours.
        if output.windows(2).any(|w| w[0] > w[1]) {
            unchanged
        } else {
            Realigned { introns: output, moved, ambiguous, kept }
        }
    }
}

pub struct Realigned {
    pub introns: Vec<i64>,
    pub moved: usize,
    /// Per intron of `introns`: its read sequence fits a nearby annotated alternative within one edit.
    pub ambiguous: Vec<bool>,
    /// Per intron of `introns`: unannotated, and the read's bases favour it over every annotated intron nearby by at
    /// least `SETTLE_ADVANTAGE` edits.
    pub kept: Vec<bool>,
}

/// Bit mask over a read's internal splice sites in transcript order (bit i is
/// `positions[i + 1]`), set for both sites of every ambiguous intron. `ambiguous`
/// is in genomic order; on the minus strand the sites run the other way.
pub fn ambiguous_site_mask(ambiguous: &[bool], strand: Strand) -> u64 {
    let n_sites = 2 * ambiguous.len();
    let mut mask = 0u64;
    for (k, _) in ambiguous.iter().enumerate().filter(|(_, &a)| a) {
        for genomic in [2 * k, 2 * k + 1] {
            let i = if strand == Strand::Minus { n_sites - 1 - genomic } else { genomic };
            if i < 64 {
                mask |= 1u64 << i;
            }
        }
    }
    mask
}

#[cfg(test)]
mod tests {
    use super::*;
    use bam::record::{Cigar, CigarString};

    /// Genome: exon A (0..100), intron with an annotated acceptor and a cryptic
    /// AG six bases earlier, exon B. Reads carry exon A then exon B's bases from
    /// either the annotated or the cryptic acceptor.
    fn fixture() -> (tempfile::TempDir, GenomeRef, SpliceRealignIndex) {
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("genome.fa");
        let mut g: Vec<u8> = Vec::new();
        let unit = b"ACGTTGCAAGCTTACGGATC";
        for i in 0..400 {
            g.push(unit[(i * 7 + i / 3) % unit.len()]);
        }
        g[100..102].copy_from_slice(b"GT"); // donor, intron starts at 1-based 101
        g[192..194].copy_from_slice(b"AG"); // cryptic acceptor, intron would end at 1-based 194
        g[198..200].copy_from_slice(b"AG"); // annotated acceptor, intron ends at 1-based 200
        std::fs::write(&fasta, format!(">chr1\n{}\n", String::from_utf8(g).unwrap())).unwrap();
        std::fs::write(dir.path().join("genome.fa.fai"), "chr1\t400\t6\t400\t401\n").unwrap();
        let genome = GenomeRef::open(fasta).unwrap();
        let mut index = SpliceRealignIndex { window: 30, ..Default::default() };
        index.insert("chr1", Strand::Plus, 101, 200);
        (dir, genome, index)
    }

    /// A read made of exon A and exon B's bases starting at 0-based `exon_b`,
    /// aligned with its junction ending at 1-based intron base `intron_last`.
    fn read(genome: &GenomeRef, exon_b: i64, intron_last: i64) -> bam::Record {
        let mut seq = genome.fetch_seq("chr1", 0, 100).unwrap();
        seq.extend(genome.fetch_seq("chr1", exon_b, exon_b + 100).unwrap());
        let mut r = bam::Record::new();
        r.set(b"r", Some(&CigarString(vec![Cigar::Match(100), Cigar::RefSkip((intron_last - 100) as u32), Cigar::Match(100)])), &seq, &vec![30; seq.len()]);
        r.set_pos(0);
        r
    }

    #[test]
    fn a_read_carrying_the_annotated_exon_bases_is_moved_to_the_annotated_intron() {
        let (_d, genome, index) = fixture();
        // Exon B's bases from the annotated acceptor, aligned 6 bp early onto
        // the cryptic AG.
        let r = read(&genome, 200, 194);
        let out = index.realign(&r, "chr1", Strand::Plus, &genome, &[101, 194]);
        assert_eq!((out.introns, out.moved), (vec![101, 200], 1));
    }

    #[test]
    fn a_read_whose_bases_support_the_unannotated_site_keeps_it() {
        let (_d, genome, index) = fixture();
        // Exon B's bases really begin at the cryptic acceptor: six genomic
        // bases the annotated template lacks.
        let r = read(&genome, 194, 194);
        let out = index.realign(&r, "chr1", Strand::Plus, &genome, &[101, 194]);
        assert_eq!((out.introns, out.moved), (vec![101, 194], 0));
    }

    #[test]
    fn an_annotated_intron_is_ambiguous_only_when_the_read_cannot_tell_it_from_its_neighbour() {
        let (_d, genome, mut index) = fixture();
        // Both acceptors annotated (two isoforms), six bases apart.
        index.insert("chr1", Strand::Plus, 101, 194);
        // Bases of the 200 acceptor aligned on the 194 one: the neighbour fits better.
        let out = index.realign(&read(&genome, 200, 194), "chr1", Strand::Plus, &genome, &[101, 194]);
        assert_eq!((out.introns, out.moved, out.ambiguous), (vec![101, 194], 0, vec![true]));
        // Bases of the 194 acceptor: the read supports its own site.
        let out = index.realign(&read(&genome, 194, 194), "chr1", Strand::Plus, &genome, &[101, 194]);
        assert_eq!(out.ambiguous, vec![false]);
    }

    #[test]
    fn ambiguous_sites_are_indexed_in_transcript_order() {
        // Second of two introns ambiguous: genomic sites 2 and 3.
        assert_eq!(ambiguous_site_mask(&[false, true], Strand::Plus), 0b1100);
        // On the minus strand the transcript runs right to left.
        assert_eq!(ambiguous_site_mask(&[false, true], Strand::Minus), 0b0011);
    }

    #[test]
    fn annotated_introns_are_left_alone() {
        let (_d, genome, index) = fixture();
        let r = read(&genome, 200, 200);
        let out = index.realign(&r, "chr1", Strand::Plus, &genome, &[101, 200]);
        assert_eq!((out.introns, out.moved, out.ambiguous), (vec![101, 200], 0, vec![false]));
    }
}
