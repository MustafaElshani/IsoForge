//! Small RNAs: annotated one-exon transcripts of small-RNA gene types, counted from the reads that match them.
//!
//! - An unspliced primary alignment on a small RNA's strand is one of its reads when it starts no more than
//!   `SLACK_5P` bp before the annotated start, ends no more than `SLACK_3P` bp past the annotated end and has at
//!   least half of its aligned span on the small RNA.
//! - No end variant is split from a small RNA, and no other read counts for it.
//! - A small RNA inside an exon of another gene on its strand takes only the reads that start at its own start.

use std::path::Path;

use rustc_hash::{FxHashMap, FxHashSet};
use serde::{Deserialize, Serialize};

use crate::consistent_assignment::containing_exon;
use crate::data_structures::{MonoexonRecord, SharedRefs};
use crate::errors::{IsoForgeError, IsoForgeResult};
use crate::types::{GeneId, RunConfig, Strand, TranscriptId, TranscriptModel};

/// Bases a read may start before the annotated start of the small RNA holding it.
pub const SLACK_5P: i64 = 15;
/// Bases a read may end past its annotated end.
pub const SLACK_3P: i64 = 100;
/// Bases inside the annotated start within which a read still starts at the start of a small RNA nested in an exon
/// of another gene.
pub const ANCHOR_INSIDE: i64 = 5;

/// An annotated small RNA: its span (1-based, closed), strand, and whether it lies inside an exon of another gene on
/// that strand.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct SmallRna {
    pub lo: i64,
    pub hi: i64,
    pub minus: bool,
    pub nested: bool,
}

impl SmallRna {
    /// Whether the alignment `a..=b` (either order) on the small RNA's strand is one of its reads.
    pub fn holds(&self, a: i64, b: i64) -> bool {
        let (lo, hi) = (a.min(b), a.max(b));
        let (before, past) = if self.minus { (hi - self.hi, self.lo - lo) } else { (self.lo - lo, hi - self.hi) };
        let shared = hi.min(self.hi) - lo.max(self.lo) + 1;
        before <= SLACK_5P && past <= SLACK_3P && 2 * shared >= hi - lo + 1 && (!self.nested || before >= -ANCHOR_INSIDE)
    }

    fn strand(&self) -> Strand {
        if self.minus { Strand::Minus } else { Strand::Plus }
    }
}

/// The run database's table of small-RNA transcripts; absent when the run had the rule off.
pub const IDS_FILE: &str = "small_rna_transcripts.json";

/// The small-RNA transcripts of a run, from its database: empty when the rule was off.
pub fn load(run_db: &Path) -> IsoForgeResult<FxHashMap<TranscriptId, SmallRna>> {
    let path = run_db.join(IDS_FILE);
    if !path.exists() {
        return Ok(FxHashMap::default());
    }
    let rows: Vec<(TranscriptId, SmallRna)> = serde_json::from_str(&std::fs::read_to_string(&path)?).map_err(invalid)?;
    Ok(rows.into_iter().collect())
}

fn invalid(e: serde_json::Error) -> IsoForgeError {
    IsoForgeError::IoError(std::io::Error::new(std::io::ErrorKind::InvalidData, e))
}

/// The strands on which a small RNA holds an interval.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Holding {
    pub plus: bool,
    pub minus: bool,
}

impl Holding {
    pub fn any(&self) -> bool {
        self.plus || self.minus
    }

    pub fn on(&self, strand: Strand) -> bool {
        match strand {
            Strand::Plus => self.plus,
            Strand::Minus => self.minus,
        }
    }
}

/// The annotated small-RNA transcripts, by chromosome.
#[derive(Debug, Default)]
pub struct SmallRnaIndex {
    /// Sorted by start.
    by_chrom: FxHashMap<String, Vec<SmallRna>>,
    by_id: FxHashMap<TranscriptId, SmallRna>,
    /// Longest span among the entries, which bounds the search.
    max_span: i64,
}

impl SmallRnaIndex {
    /// From (chromosome, start, end, strand, transcript, nested in an exon of another gene) records, coordinates
    /// 1-based and closed.
    pub fn from_records<'a>(records: impl IntoIterator<Item = (&'a str, i64, i64, Strand, TranscriptId, bool)>) -> Self {
        let mut index = Self::default();
        for (chrom, a, b, strand, transcript_id, nested) in records {
            let entry = SmallRna { lo: a.min(b), hi: a.max(b), minus: strand == Strand::Minus, nested };
            index.by_chrom.entry(chrom.to_string()).or_default().push(entry);
            index.by_id.insert(transcript_id, entry);
            index.max_span = index.max_span.max(entry.hi - entry.lo + 1);
        }
        for entries in index.by_chrom.values_mut() {
            entries.sort_unstable_by_key(|e| (e.lo, e.hi));
        }
        index
    }

    /// The one-exon transcripts of genes whose type is in `config.small_rna_types`, no longer than
    /// `config.small_rna_max_length`; each is marked when it lies inside an exon of an annotated transcript of
    /// another gene on its strand that is no small RNA itself.
    pub fn from_reference(
        refs: &SharedRefs,
        gene_strand: &FxHashMap<GeneId, Strand>,
        gene_types: &FxHashMap<GeneId, String>,
        config: &RunConfig,
    ) -> Self {
        let types: FxHashSet<&str> = config.small_rna_types.iter().map(|t| t.as_str()).collect();
        let small: Vec<(&TranscriptModel, Strand)> = refs.transcript_dict.values().filter_map(|t| {
            if t.n_exons != 1 || (t.start_pos - t.end_pos).abs() + 1 > config.small_rna_max_length {
                return None;
            }
            if !gene_types.get(&t.gene_id).is_some_and(|ty| types.contains(ty.as_str())) {
                return None;
            }
            Some((t, *gene_strand.get(&t.gene_id)?))
        }).collect();
        let ids: FxHashSet<TranscriptId> = small.iter().map(|(t, _)| t.transcript_id).collect();
        Self::from_records(small.iter().map(|&(t, strand)| {
            let (lo, hi) = (t.start_pos.min(t.end_pos), t.start_pos.max(t.end_pos));
            let nested = refs.transcript_tree.query_overlap(&t.chromosome, lo, hi).iter().any(|host| {
                if host.strand != strand || host.gene_id == t.gene_id || ids.contains(&host.transcript_id) {
                    return false;
                }
                match refs.parent_chains.as_ref().and_then(|chains| chains.get(host.transcript_id)) {
                    Some(chain) => containing_exon(chain, lo, hi, 0).is_some(),
                    // a one-exon host
                    None => host.min_pos <= lo && hi <= host.max_pos,
                }
            });
            (t.chromosome.as_str(), t.start_pos, t.end_pos, strand, t.transcript_id, nested)
        }))
    }

    pub fn len(&self) -> usize {
        self.by_id.len()
    }

    pub fn is_empty(&self) -> bool {
        self.by_id.is_empty()
    }

    /// How many of the small RNAs lie inside an exon of another gene.
    pub fn n_nested(&self) -> usize {
        self.by_id.values().filter(|s| s.nested).count()
    }

    pub fn is_small_rna(&self, transcript_id: TranscriptId) -> bool {
        self.by_id.contains_key(&transcript_id)
    }

    /// Writes the small RNAs into the run database, where post-processing reads them (`load`).
    pub fn write(&self, run_db: &Path) -> IsoForgeResult<()> {
        let mut rows: Vec<(TranscriptId, SmallRna)> = self.by_id.iter().map(|(&t, &s)| (t, s)).collect();
        rows.sort_unstable_by_key(|r| r.0);
        std::fs::write(run_db.join(IDS_FILE), serde_json::to_string(&rows).map_err(invalid)?)?;
        Ok(())
    }

    /// The strands on which a small RNA holds the alignment `lo..=hi`.
    pub fn holding(&self, chrom: &str, lo: i64, hi: i64) -> Holding {
        let mut held = Holding::default();
        let Some(entries) = self.by_chrom.get(chrom) else { return held };
        // the entries starting from lo - max_span + 1 to hi
        let first = entries.partition_point(|e| e.lo < lo - self.max_span + 1);
        for e in &entries[first..] {
            if e.lo > hi {
                break;
            }
            if e.holds(lo, hi) {
                match e.strand() {
                    Strand::Plus => held.plus = true,
                    Strand::Minus => held.minus = true,
                }
            }
        }
        held
    }

    /// The annotated one-exon transcripts a small-RNA read at `lo..=hi` may be assigned to: the small RNAs holding
    /// it, or the small RNAs among the candidates when none holds it.
    pub fn holders(&self, candidates: Vec<MonoexonRecord>, lo: i64, hi: i64) -> Vec<MonoexonRecord> {
        let small: Vec<MonoexonRecord> = candidates.into_iter().filter(|r| self.is_small_rna(r.transcript_id)).collect();
        let holding: Vec<MonoexonRecord> = small.iter()
            .filter(|r| self.by_id.get(&r.transcript_id).is_some_and(|s| s.holds(lo, hi)))
            .cloned()
            .collect();
        // a read that no small RNA holds stays with the small RNAs among its candidates
        if holding.is_empty() { small } else { holding }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn index() -> SmallRnaIndex {
        // a small RNA on the plus strand, one on the minus strand of another chromosome, two nested entries
        SmallRnaIndex::from_records([
            ("chr21", 8212572, 8212724, Strand::Plus, 1, false),
            ("chr1", 5000, 5129, Strand::Minus, 2, false),
            ("chr2", 1000, 1300, Strand::Plus, 3, false),
            ("chr2", 1100, 1180, Strand::Plus, 4, false),
            // a small RNA inside an exon of another gene, on either strand
            ("chr3", 3000, 3128, Strand::Plus, 5, true),
            ("chr3", 7000, 7128, Strand::Minus, 6, true),
        ])
    }

    #[test]
    fn an_alignment_inside_a_small_rna_is_held_on_its_strand() {
        let ix = index();
        // a read starting one base upstream of the annotated start and ending a few past the end
        let held = ix.holding("chr21", 8212571, 8212730);
        assert!(held.on(Strand::Plus) && !held.on(Strand::Minus));
        assert!(ix.holding("chr1", 5010, 5100).on(Strand::Minus));
        assert!(!ix.holding("chr1", 5010, 5100).on(Strand::Plus));
    }

    #[test]
    fn a_read_reaching_beyond_the_allowances_is_not_held() {
        let ix = index();
        // a read starting far upstream of the small RNA, and reads just past the allowance at either end
        assert!(!ix.holding("chr21", 8210554, 8212695).any());
        assert!(!ix.holding("chr21", 8212572 - SLACK_5P - 1, 8212700).any());
        assert!(!ix.holding("chr21", 8212600, 8212724 + SLACK_3P + 1).any());
        assert!(ix.holding("chr21", 8212572 - SLACK_5P, 8212724 + SLACK_3P).any());
        assert!(!ix.holding("chr5", 8212600, 8212700).any());
    }

    #[test]
    fn the_wider_allowance_is_at_the_three_prime_end_of_either_strand() {
        let ix = index();
        // the minus-strand small RNA at 5000..=5129 starts at 5129: 60 bp past its end is below 5000
        assert!(ix.holding("chr1", 4940, 5129).on(Strand::Minus));
        assert!(!ix.holding("chr1", 5000, 5189).any());
        // a read on the plus strand ending 65 bp past the annotated end
        assert!(ix.holding("chr21", 8212572, 8212724 + 65).on(Strand::Plus));
    }

    #[test]
    fn a_read_with_most_of_its_span_off_the_small_rna_is_not_held() {
        let ix = index();
        // starts in the last 30 bases of the 153-nt entry and runs 90 past it: within the allowances, mostly off it
        assert!(!ix.holding("chr21", 8212695, 8212724 + 90).any());
        // a fragment inside counts in full
        assert!(ix.holding("chr21", 8212600, 8212700).on(Strand::Plus));
    }

    #[test]
    fn a_small_rna_nested_in_another_genes_exon_takes_only_reads_starting_at_its_start() {
        let ix = index();
        assert_eq!(ix.n_nested(), 2);
        // from its start, or a few bases either side of it
        assert!(ix.holding("chr3", 3000, 3128).on(Strand::Plus));
        assert!(ix.holding("chr3", 3000 - SLACK_5P, 3160).on(Strand::Plus));
        assert!(ix.holding("chr3", 3000 + ANCHOR_INSIDE, 3128).on(Strand::Plus));
        // a fragment of the host starting further in, which an unnested small RNA would take
        assert!(!ix.holding("chr3", 3000 + ANCHOR_INSIDE + 1, 3128).any());
        assert!(!ix.holding("chr3", 3040, 3200).any());
        assert!(ix.holding("chr21", 8212600, 8212700).on(Strand::Plus));
        // on the minus strand the start is the higher coordinate
        assert!(ix.holding("chr3", 6950, 7128).on(Strand::Minus));
        assert!(ix.holding("chr3", 6950, 7128 - ANCHOR_INSIDE).on(Strand::Minus));
        assert!(!ix.holding("chr3", 7000, 7100).any());
    }

    #[test]
    fn the_small_rnas_round_trip_through_the_run_database() {
        let dir = tempfile::tempdir().unwrap();
        assert!(load(dir.path()).unwrap().is_empty());
        index().write(dir.path()).unwrap();
        let small = load(dir.path()).unwrap();
        assert_eq!(small.len(), 6);
        assert!(!small.contains_key(&9));
        assert_eq!(small[&6], SmallRna { lo: 7000, hi: 7128, minus: true, nested: true });
        // a loaded small RNA holds the same reads, with ends in either order
        assert!(small[&1].holds(8212571, 8212730) && small[&1].holds(8212730, 8212571));
        assert!(!small[&1].holds(8210554, 8212695));
        assert!(small[&6].holds(7128, 6950) && !small[&6].holds(7100, 7000));
    }

    #[test]
    fn nested_small_rnas_are_both_found() {
        let ix = index();
        assert!(ix.holding("chr2", 1105, 1175).any());
        // inside the longer one only
        assert!(ix.holding("chr2", 1005, 1290).any());
        assert_eq!(ix.len(), 6);
    }

    #[test]
    fn a_small_rna_read_goes_to_the_small_rnas_holding_it() {
        let ix = index();
        let rec = |transcript_id, start, end| MonoexonRecord {
            gene_id: transcript_id, transcript_id, chromosome: "chr2".to_string(), strand: Strand::Plus,
            start, end, start_vertex: 0, end_vertex: 0, exon_id: 0,
        };
        // transcript 9 is an annotated one-exon transcript that is no small RNA
        let cands = || vec![rec(3, 1000, 1300), rec(4, 1100, 1180), rec(9, 900, 2000)];
        let ids = |v: Vec<MonoexonRecord>| v.into_iter().map(|r| r.transcript_id).collect::<Vec<_>>();
        // inside both small RNAs, or inside the longer one only; never the other annotated transcript
        assert_eq!(ids(ix.holders(cands(), 1105, 1175)), vec![3, 4]);
        assert_eq!(ids(ix.holders(cands(), 1005, 1290)), vec![3]);
        // a small-RNA read whose corrected ends no small RNA holds stays with the small RNAs it overlaps
        assert_eq!(ids(ix.holders(cands(), 1005, 1900)), vec![3, 4]);
    }
}
