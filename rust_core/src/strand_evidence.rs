//! Strand of reads that carry no strand evidence of their own.
//!
//! This module measures, from the unspliced reads that carry a poly(A) tail, whether
//! a dataset's alignment orientation carries strand information. When it does not, a
//! read without evidence of its own takes the strand of the tailed reads covering it,
//! failing that of the annotated gene there, and is set aside when neither exists.

use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering::Relaxed};

use rust_htslib::bam::{self, Read as BamRead};
use rustc_hash::FxHashMap;

use crate::bam_reader::{extract_soft_clips, strand_from_polya};
use crate::data_structures::SharedRefs;
use crate::errors::{IsoForgeError, IsoForgeResult};
use crate::types::{RunConfig, Strand};

/// Width of the genomic bins that hold the strand votes.
const BIN: i64 = 250;
/// Tailed unspliced reads seen before a library's orientation is called.
const SURVEY_MIN_TAILED: u64 = 2_000;
/// Tail/orientation agreement above which alignment orientation is trusted.
const ORIENTED_AGREEMENT: f64 = 0.9;

/// Per-dataset strand lookup for reads without strand evidence of their own.
pub struct StrandIndex {
    /// Alignment orientation of unspliced reads carries strand information.
    pub oriented: bool,
    pub tailed_unspliced: u64,
    pub agreement: f64,
    /// Tailed reads covering each bin; only built for an unoriented library.
    votes: FxHashMap<String, [Vec<u32>; 2]>,
    /// Annotated gene bodies per bin.
    genes: FxHashMap<String, [Vec<u8>; 2]>,
    pub from_locus: AtomicU64,
    pub from_gene: AtomicU64,
    pub unresolved: AtomicU64,
}

impl std::fmt::Debug for StrandIndex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("StrandIndex")
            .field("oriented", &self.oriented)
            .field("tailed_unspliced", &self.tailed_unspliced)
            .field("agreement", &self.agreement)
            .finish()
    }
}

fn idx(strand: Strand) -> usize {
    match strand {
        Strand::Plus => 0,
        Strand::Minus => 1,
    }
}

fn bins(start: i64, end: i64) -> (usize, usize) {
    ((start.min(end).max(1) / BIN) as usize, (start.max(end) / BIN) as usize)
}

impl StrandIndex {
    /// Survey a BAM's orientation and, for an unoriented library, index its tailed reads.
    pub fn build(bam_path: &Path, config: &RunConfig, refs: &SharedRefs) -> IsoForgeResult<Self> {
        let mut reader = bam::Reader::from_path(bam_path)
            .map_err(|e| IsoForgeError::BamError(format!("{}: {}", bam_path.display(), e)))?;
        reader
            .set_threads(config.threads.clamp(1, 4))
            .map_err(|e| IsoForgeError::HtsError(e.to_string()))?;
        let header = reader.header().clone();
        let names: Vec<String> = header
            .target_names()
            .iter()
            .map(|n| String::from_utf8_lossy(n).to_string())
            .collect();
        let lengths: Vec<u64> = (0..names.len()).map(|i| header.target_len(i as u32).unwrap_or(0)).collect();

        // Votes accumulate as difference arrays and are summed once the pass is over.
        let mut diffs: Vec<Option<[Vec<i32>; 2]>> = (0..names.len()).map(|_| None).collect();
        let (mut tailed, mut agree) = (0u64, 0u64);
        let mut oriented: Option<bool> = None;
        let mut record = bam::Record::new();
        while let Some(r) = reader.read(&mut record) {
            r.map_err(|e| IsoForgeError::HtsError(format!("Failed to read record: {}", e)))?;
            let flag = record.flags();
            if (flag != 0 && flag != 16) || record.mapq() < config.min_mapq || record.tid() < 0 {
                continue;
            }
            let cigar = record.cigar();
            let spliced = cigar.iter().any(|op| matches!(op, bam::record::Cigar::RefSkip(_)));
            let (lead, trail) = extract_soft_clips(&cigar);
            let aln = if record.is_reverse() { Strand::Minus } else { Strand::Plus };
            let tail = strand_from_polya(&record.seq().as_bytes(), lead, trail);
            let Some(vote) = tail else { continue };
            if !spliced {
                tailed += 1;
                agree += (vote == aln) as u64;
                if oriented.is_none() && tailed >= SURVEY_MIN_TAILED {
                    let a = agree as f64 / tailed as f64;
                    oriented = Some(a >= ORIENTED_AGREEMENT);
                    if oriented == Some(true) {
                        break;
                    }
                }
            }
            let tid = record.tid() as usize;
            let d = diffs[tid].get_or_insert_with(|| {
                let n = (lengths[tid] / BIN as u64) as usize + 2;
                [vec![0; n], vec![0; n]]
            });
            let (lo, hi) = bins(record.pos() + 1, cigar.end_pos());
            let d = &mut d[idx(vote)];
            let hi = hi.min(d.len() - 2);
            d[lo] += 1;
            d[hi + 1] -= 1;
        }
        // Too few tailed reads to judge: keep the alignment orientation.
        let oriented = oriented.unwrap_or(true);
        let agreement = if tailed > 0 { agree as f64 / tailed as f64 } else { 0.0 };

        let mut votes = FxHashMap::default();
        if !oriented {
            for (tid, d) in diffs.into_iter().enumerate() {
                let Some(d) = d else { continue };
                let sum = |v: Vec<i32>| {
                    let mut acc = 0i64;
                    v.into_iter().map(|x| { acc += x as i64; acc.max(0) as u32 }).collect::<Vec<u32>>()
                };
                let [p, m] = d;
                votes.insert(names[tid].clone(), [sum(p), sum(m)]);
            }
        }

        let mut genes: FxHashMap<String, [Vec<u8>; 2]> = FxHashMap::default();
        if !oriented {
            for g in refs.gene_tree.records() {
                let Some(tid) = names.iter().position(|n| *n == g.chromosome) else { continue };
                let n = (lengths[tid] / BIN as u64) as usize + 2;
                let e = genes.entry(g.chromosome.clone()).or_insert_with(|| [vec![0; n], vec![0; n]]);
                let (lo, hi) = bins(g.start, g.end);
                for b in lo..=hi.min(n - 1) {
                    e[idx(g.strand)][b] = 1;
                }
            }
        }

        Ok(Self {
            oriented,
            tailed_unspliced: tailed,
            agreement,
            votes,
            genes,
            from_locus: AtomicU64::new(0),
            from_gene: AtomicU64::new(0),
            unresolved: AtomicU64::new(0),
        })
    }

    /// Strand for a read without evidence of its own at `start..end`: the tailed
    /// reads covering it when they agree, else the annotated gene there, else none.
    pub fn resolve(&self, chrom: &str, start: i64, end: i64) -> Option<Strand> {
        let (lo, hi) = bins(start, end);
        if let Some(v) = self.votes.get(chrom) {
            let sum = |a: &Vec<u32>| a[lo.min(a.len() - 1)..=hi.min(a.len() - 1)].iter().map(|&x| x as u64).sum::<u64>();
            let (p, m) = (sum(&v[0]), sum(&v[1]));
            if p + m > 0 {
                // Both strands hold at least a third of the votes here: the read is ambiguous.
                if p.min(m) * 3 >= p + m {
                    self.unresolved.fetch_add(1, Relaxed);
                    return None;
                }
                self.from_locus.fetch_add(1, Relaxed);
                return Some(if p > m { Strand::Plus } else { Strand::Minus });
            }
        }
        if let Some(g) = self.genes.get(chrom) {
            let any = |a: &Vec<u8>| a[lo.min(a.len() - 1)..=hi.min(a.len() - 1)].iter().any(|&x| x > 0);
            match (any(&g[0]), any(&g[1])) {
                (true, false) => {
                    self.from_gene.fetch_add(1, Relaxed);
                    return Some(Strand::Plus);
                }
                (false, true) => {
                    self.from_gene.fetch_add(1, Relaxed);
                    return Some(Strand::Minus);
                }
                _ => {}
            }
        }
        self.unresolved.fetch_add(1, Relaxed);
        None
    }

    pub fn log_survey(&self, dataset: &str) {
        if self.oriented {
            eprintln!(
                "[IsoForge] {}: alignment orientation agrees with poly(A) tails on {:.1}% of {} tailed unspliced reads; orientation kept for reads without their own evidence",
                dataset, 100.0 * self.agreement, self.tailed_unspliced
            );
        } else {
            eprintln!(
                "[IsoForge] {}: alignment orientation agrees with poly(A) tails on only {:.1}% of {} tailed unspliced reads; reads without a tail or clean splice motifs take the strand of the tailed reads or annotated gene covering them",
                dataset, 100.0 * self.agreement, self.tailed_unspliced
            );
        }
    }

    /// Clears the resolution counts, for an index read again by a later pass over the same alignments.
    pub fn reset_counts(&self) {
        for c in [&self.from_locus, &self.from_gene, &self.unresolved] { c.store(0, Relaxed); }
    }

    pub fn log_resolution(&self, dataset: &str) {
        if self.oriented {
            return;
        }
        eprintln!(
            "[IsoForge] {}: reads without their own strand evidence stranded from covering tailed reads: {}, from the annotated gene: {}, set aside as ambiguous or unsupported: {}",
            dataset,
            self.from_locus.load(Relaxed),
            self.from_gene.load(Relaxed),
            self.unresolved.load(Relaxed)
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn index(votes: [Vec<u32>; 2], genes: [Vec<u8>; 2]) -> StrandIndex {
        let mut v = FxHashMap::default();
        v.insert("chr1".to_string(), votes);
        let mut g = FxHashMap::default();
        g.insert("chr1".to_string(), genes);
        StrandIndex {
            oriented: false,
            tailed_unspliced: 0,
            agreement: 0.0,
            votes: v,
            genes: g,
            from_locus: AtomicU64::new(0),
            from_gene: AtomicU64::new(0),
            unresolved: AtomicU64::new(0),
        }
    }

    #[test]
    fn covering_reads_decide_the_strand() {
        let ix = index([vec![0, 30, 30, 0], vec![0, 1, 0, 0]], [vec![0; 4], vec![0; 4]]);
        assert_eq!(ix.resolve("chr1", 300, 600), Some(Strand::Plus));
        assert_eq!(ix.from_locus.load(Relaxed), 1);
    }

    #[test]
    fn split_support_is_ambiguous() {
        let ix = index([vec![0, 10, 0, 0], vec![0, 6, 0, 0]], [vec![0; 4], vec![0; 4]]);
        assert_eq!(ix.resolve("chr1", 260, 400), None);
        assert_eq!(ix.unresolved.load(Relaxed), 1);
    }

    #[test]
    fn annotated_gene_is_the_fallback() {
        let ix = index([vec![0; 4], vec![0; 4]], [vec![0; 4], vec![0, 1, 1, 0]]);
        assert_eq!(ix.resolve("chr1", 300, 600), Some(Strand::Minus));
        assert_eq!(ix.from_gene.load(Relaxed), 1);
        assert_eq!(ix.resolve("chr1", 800, 900), None);
        assert_eq!(ix.resolve("chr2", 300, 600), None);
    }

    #[test]
    fn genes_on_both_strands_do_not_decide() {
        let ix = index([vec![0; 4], vec![0; 4]], [vec![0, 1, 0, 0], vec![0, 1, 0, 0]]);
        assert_eq!(ix.resolve("chr1", 300, 400), None);
    }
}
