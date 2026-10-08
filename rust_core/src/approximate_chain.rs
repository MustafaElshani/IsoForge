//! Chain-level approximate assignment of a read to an annotated transcript.
//!
//! A read is assigned to the annotated transcript whose internal splice
//! positions all lie within `tolerance` of the read's. A read splice site on an
//! annotated site of its role (donor or acceptor) is locked and may not move
//! onto a different annotated site. Candidates are ranked by total offset, then
//! the largest single offset, then agreement of the read's own 5'/3' ends with
//! the model; a read is declined only when two models tie on every criterion.

use crate::data_structures::SharedRefs;
use crate::types::{EdgeType, GeneId, Strand, TranscriptId};
use rustc_hash::FxHashMap;

/// An annotated transcript flattened into everything the classifier needs.
#[derive(Debug, Clone, Default)]
pub struct ApproxModel {
    pub transcript_id: TranscriptId,
    pub gene_id: GeneId,
    /// Annotated TSS and TES in transcript orientation.
    pub start_pos: i64,
    pub end_pos: i64,
    pub start_vertex: i64,
    pub end_vertex: i64,
    pub start_exon: i64,
    pub end_exon: i64,
    /// Internal edges in annotation order (introns and exons alternating).
    pub internal_edges: Vec<i64>,
    /// Internal splice vertices, 5'→3'.
    pub internal_vertices: Vec<i64>,
}

/// One annotated transcript, plus the coordinates a read is compared against.
#[derive(Debug, Clone)]
struct ChainEntry {
    internal: Vec<i64>,
    model: ApproxModel,
}

/// Annotated splice sites of one (chromosome, strand), by role, sorted.
#[derive(Debug, Clone, Default)]
struct SiteSets {
    donors: Vec<i64>,
    acceptors: Vec<i64>,
}

/// Annotated chains bucketed by chromosome, strand and exon count, sorted
/// within a bucket by first internal position.
#[derive(Debug, Default)]
pub struct ApproxChainIndex {
    buckets: FxHashMap<(String, Strand, i32), Vec<ChainEntry>>,
    sites: FxHashMap<(String, Strand), SiteSets>,
    /// Whether a query site on an annotated site of its role is locked.
    lock_annotated: bool,
}

/// What an approximate match resolved to.
pub struct ApproxMatch<'a> {
    pub model: &'a ApproxModel,
    /// Sum of absolute offsets across internal splice positions.
    pub total_offset: i64,
    /// Largest single offset.
    pub max_offset: i64,
    /// Junctions that had to move to reach the model.
    pub n_shifted: usize,
    /// Next-best candidate and its total offset.
    pub runner_up: Option<(TranscriptId, i64)>,
    /// Per-junction offsets to the chosen model, 5'→3'.
    pub offsets: Vec<i64>,
}

/// The resolution of a query, used for assignment and for the trace.
pub struct Resolution<'a> {
    pub matched: Option<ApproxMatch<'a>>,
    /// The two models that could not be told apart, when the query declined.
    pub tie: Option<(TranscriptId, TranscriptId)>,
    /// Which query junctions were locked to their annotated site.
    pub locked: Vec<bool>,
}

pub const CHAIN_TRACE_HEADER: &str = "stage\tchrom\tstrand\tn_exons\traw_internal\tlocked\toutcome\ttranscript_id\tgene_id\toffsets\tmax_offset\tn_shifted\trunner_up\trunner_up_total\tmargin";

impl<'a> Resolution<'a> {
    /// One TSV row describing this decision.
    pub fn trace_row(&self, stage: &str, chrom: &str, strand: Strand, positions: &[i64]) -> String {
        let join = |v: &[i64]| v.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",");
        let raw = join(&positions[1..positions.len() - 1]);
        let locked: String = self.locked.iter().map(|&l| if l { '1' } else { '0' }).collect();
        let strand_c = if strand == Strand::Plus { '+' } else { '-' };
        let n_exons = positions.len() / 2;
        match (&self.matched, self.tie) {
            (Some(m), _) => {
                let (ru, ru_total, margin) = match m.runner_up {
                    Some((t, tot)) => (t.to_string(), tot.to_string(), (tot - m.total_offset).to_string()),
                    None => (".".into(), ".".into(), ".".into()),
                };
                format!(
                    "{stage}\t{chrom}\t{strand_c}\t{n_exons}\t{raw}\t{locked}\tassigned\t{}\t{}\t{}\t{}\t{}\t{ru}\t{ru_total}\t{margin}",
                    m.model.transcript_id, m.model.gene_id, join(&m.offsets), m.max_offset, m.n_shifted
                )
            }
            (None, Some((a, b))) => format!(
                "{stage}\t{chrom}\t{strand_c}\t{n_exons}\t{raw}\t{locked}\tdeclined_tie\t{a}\t.\t.\t.\t.\t{b}\t.\t0"
            ),
            (None, None) => format!(
                "{stage}\t{chrom}\t{strand_c}\t{n_exons}\t{raw}\t{locked}\tdeclined_none\t.\t.\t.\t.\t.\t.\t.\t."
            ),
        }
    }
}

impl ApproxChainIndex {
    /// Build from the reference annotation (multi-exon models only, strand from
    /// the intron edges). Annotated donors and acceptors come from intron edges
    /// only; `v1` of an intron edge is its donor on both strands.
    pub fn from_reference(refs: &SharedRefs) -> Self {
        let (positions, _) = reference_maps(refs);
        let chrom_of: FxHashMap<i64, &str> = refs
            .location_dict
            .chromosomes
            .iter()
            .flat_map(|(chrom, c)| c.positions.iter().map(move |&(_, v)| (v, chrom.as_str())))
            .collect();
        let mut index = Self { lock_annotated: refs.config.approx_chain_lock_annotated, ..Self::default() };
        for e in refs.edge_dict.values() {
            if e.edge_type != EdgeType::Intron {
                continue;
            }
            let (Some(&d), Some(&a), Some(&chrom)) = (positions.get(&e.v1), positions.get(&e.v2), chrom_of.get(&e.v1))
            else {
                continue;
            };
            let s = index.sites.entry((chrom.to_string(), e.strand)).or_default();
            s.donors.push(d);
            s.acceptors.push(a);
        }
        for s in index.sites.values_mut() {
            s.donors.sort_unstable();
            s.donors.dedup();
            s.acceptors.sort_unstable();
            s.acceptors.dedup();
        }
        for tx in refs.transcript_dict.values() {
            if tx.n_exons < 2 {
                continue;
            }
            let Some((internal, internal_vertices, internal_edges, strand)) = refs.hypergraph.internal(tx.transcript_id)
            else {
                continue;
            };
            if internal.len() != (tx.n_exons as usize - 1) * 2 {
                continue;
            }
            let (start_pos, end_pos) = match strand {
                Strand::Plus => (tx.start_pos, tx.end_pos),
                Strand::Minus => (tx.end_pos, tx.start_pos),
            };
            index
                .buckets
                .entry((tx.chromosome.clone(), strand, tx.n_exons))
                .or_default()
                .push(ChainEntry {
                    internal,
                    model: ApproxModel {
                        transcript_id: tx.transcript_id,
                        gene_id: tx.gene_id,
                        start_pos,
                        end_pos,
                        start_vertex: tx.start_vertex,
                        end_vertex: tx.end_vertex,
                        start_exon: tx.start_exon,
                        end_exon: tx.end_exon,
                        internal_edges,
                        internal_vertices,
                    },
                });
        }
        for entries in index.buckets.values_mut() {
            entries.sort_by_key(|e| (e.internal[0], e.model.transcript_id));
        }
        index
    }

    pub fn is_empty(&self) -> bool {
        self.buckets.is_empty()
    }

    /// Best annotated transcript whose every internal splice position lies
    /// within `tolerance` of the read's, or None when nothing fits or two
    /// models are indistinguishable.
    ///
    /// `positions` is the read's full coordinate vector, 5'→3', so the internal
    /// splice positions are everything but the first and last entry.
    pub fn find(&self, chrom: &str, strand: Strand, positions: &[i64], tolerance: i64) -> Option<ApproxMatch<'_>> {
        self.resolve(chrom, strand, positions, tolerance).matched
    }

    /// `find`, plus why the query resolved the way it did.
    pub fn resolve(&self, chrom: &str, strand: Strand, positions: &[i64], tolerance: i64) -> Resolution<'_> {
        self.resolve_with_unlocked(chrom, strand, positions, tolerance, 0, 0)
    }

    /// `resolve`, with internal sites that are never locked (`unlocked`) and
    /// novel sites the splice-site consensus settled, which are always locked
    /// (`settled`); bit i of either is positions[i + 1].
    pub fn resolve_with_unlocked(&self, chrom: &str, strand: Strand, positions: &[i64], tolerance: i64, unlocked: u64,
                                 settled: u64) -> Resolution<'_> {
        let none = Resolution { matched: None, tie: None, locked: Vec::new() };
        let n_exons = positions.len() as i32 / 2;
        if n_exons < 2 || tolerance <= 0 {
            return none;
        }
        let query = &positions[1..positions.len() - 1];

        // Internal positions alternate donor, acceptor in transcript order on
        // both strands, so the role of a query site is its index parity.
        let sites = self.sites.get(&(chrom.to_string(), strand));
        let locked: Vec<bool> = query
            .iter()
            .enumerate()
            .map(|(i, q)| {
                (i < 64 && settled & (1u64 << i) != 0)
                    || self.lock_annotated
                        && (i >= 64 || unlocked & (1u64 << i) == 0)
                        && sites.is_some_and(|s| {
                            let set = if i % 2 == 0 { &s.donors } else { &s.acceptors };
                            set.binary_search(q).is_ok()
                        })
            })
            .collect();

        let Some(entries) = self.buckets.get(&(chrom.to_string(), strand, n_exons)) else {
            return Resolution { locked, ..none };
        };

        // The scan starts at the first chain whose first internal position is
        // at least `first - tolerance`.
        let first = query[0];
        let lo = entries.partition_point(|e| e.internal[0] < first - tolerance);
        let read_5p = positions[0];
        let read_3p = *positions.last().unwrap();

        // (total, worst, ends, transcript, entry index)
        let mut best: Option<(i64, i64, i64, TranscriptId, usize)> = None;
        let mut second: Option<(i64, TranscriptId)> = None;
        let mut tied_with: Option<TranscriptId> = None;
        for (offset, e) in entries[lo..].iter().enumerate() {
            if e.internal[0] > first + tolerance {
                break;
            }
            let mut total = 0i64;
            let mut worst = 0i64;
            let mut fits = true;
            for (i, (a, b)) in query.iter().zip(&e.internal).enumerate() {
                let d = (a - b).abs();
                if d > tolerance || (locked[i] && d != 0) {
                    fits = false;
                    break;
                }
                total += d;
                worst = worst.max(d);
            }
            if !fits {
                continue;
            }
            // Agreement of the read's own ends with this model, for chain-level ties.
            let ends = (read_5p - e.model.start_pos).abs() + (read_3p - e.model.end_pos).abs();
            let key = (total, worst, ends, e.model.transcript_id, lo + offset);
            match best {
                None => best = Some(key),
                Some(cur) => {
                    if key.0 < cur.0
                        || (key.0 == cur.0 && key.1 < cur.1)
                        || (key.0 == cur.0 && key.1 == cur.1 && key.2 < cur.2)
                    {
                        second = Some((cur.0, cur.3));
                        best = Some(key);
                        tied_with = None;
                    } else {
                        if (key.0, key.1, key.2) == (cur.0, cur.1, cur.2) {
                            // Indistinguishable on chain offset and on endpoints.
                            tied_with = Some(key.3);
                        }
                        if second.is_none_or(|s| key.0 < s.0) {
                            second = Some((key.0, key.3));
                        }
                    }
                }
            }
        }
        let Some((total, worst, _, tid, idx)) = best else {
            return Resolution { locked, ..none };
        };
        if let Some(other) = tied_with {
            return Resolution { matched: None, tie: Some((tid, other)), locked };
        }
        let e = &entries[idx];
        let offsets: Vec<i64> = query.iter().zip(&e.internal).map(|(a, b)| (a - b).abs()).collect();
        let n_shifted = offsets.iter().filter(|&&d| d > 0).count();
        Resolution {
            matched: Some(ApproxMatch {
                model: &e.model,
                total_offset: total,
                max_offset: worst,
                n_shifted,
                runner_up: second.map(|(total, id)| (id, total)),
                offsets,
            }),
            tie: None,
            locked,
        }
    }
}

/// Vertex → coordinate and edge → (v1, v2, strand) lookups over the reference.
pub(crate) fn reference_maps(
    refs: &SharedRefs,
) -> (FxHashMap<i64, i64>, FxHashMap<i64, (i64, i64, Strand)>) {
    let positions = refs
        .location_dict
        .chromosomes
        .values()
        .flat_map(|c| c.positions.iter().map(|&(p, v)| (v, p)))
        .collect();
    let edges = refs
        .edge_dict
        .values()
        .map(|e| (e.edge_id, (e.v1, e.v2, e.strand)))
        .collect();
    (positions, edges)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn entry(id: i64, internal: &[i64], start: i64, end: i64) -> ChainEntry {
        ChainEntry {
            internal: internal.to_vec(),
            model: ApproxModel {
                transcript_id: id,
                gene_id: 1,
                start_pos: start,
                end_pos: end,
                ..ApproxModel::default()
            },
        }
    }

    fn index(entries: Vec<ChainEntry>) -> ApproxChainIndex {
        let mut i = ApproxChainIndex { lock_annotated: true, ..ApproxChainIndex::default() };
        let mut e = entries;
        e.sort_by_key(|x| (x.internal[0], x.model.transcript_id));
        i.buckets.insert(("chr1".into(), Strand::Plus, 3), e);
        i
    }

    fn with_sites(mut ix: ApproxChainIndex, donors: &[i64], acceptors: &[i64]) -> ApproxChainIndex {
        let mut s = SiteSets { donors: donors.to_vec(), acceptors: acceptors.to_vec() };
        s.donors.sort_unstable();
        s.acceptors.sort_unstable();
        ix.sites.insert(("chr1".into(), Strand::Plus), s);
        ix
    }

    #[test]
    fn matches_a_chain_shifted_beyond_the_snap_radius() {
        let ix = index(vec![entry(7, &[200, 300, 400, 500], 100, 600)]);
        // Every junction 40 bp off: beyond snapping, still one clear chain.
        let read = [100, 240, 340, 440, 540, 600];
        let m = ix.find("chr1", Strand::Plus, &read, 50).unwrap();
        assert_eq!(m.model.transcript_id, 7);
        assert_eq!(m.max_offset, 40);
        assert_eq!(m.total_offset, 160);
        assert_eq!(m.n_shifted, 4);
        assert_eq!(m.offsets, vec![40, 40, 40, 40]);
        assert!(m.runner_up.is_none());
        // The same read is refused at a tolerance that cannot cover the shift.
        assert!(ix.find("chr1", Strand::Plus, &read, 25).is_none());
    }

    #[test]
    fn one_junction_outside_tolerance_rejects_the_whole_chain() {
        let ix = index(vec![entry(7, &[200, 300, 400, 500], 100, 600)]);
        let read = [100, 205, 305, 405, 900, 1000];
        assert!(ix.find("chr1", Strand::Plus, &read, 50).is_none());
    }

    #[test]
    fn wrong_strand_or_exon_count_never_matches() {
        let ix = index(vec![entry(7, &[200, 300, 400, 500], 100, 600)]);
        let read = [100, 205, 305, 405, 505, 600];
        assert!(ix.find("chr1", Strand::Minus, &read, 50).is_none());
        let two_exon = [100, 205, 305, 600];
        assert!(ix.find("chr1", Strand::Plus, &two_exon, 50).is_none());
    }

    /// The read is not thrown away when two models are both in range: the
    /// nearer chain wins, and if the chains tie the read's own ends decide.
    #[test]
    fn ambiguity_is_resolved_rather_than_discarded() {
        let ix = index(vec![
            entry(7, &[200, 300, 400, 500], 100, 600),
            entry(9, &[210, 310, 410, 510], 100, 600),
        ]);
        let closer = [100, 202, 302, 402, 502, 600];
        let m = ix.find("chr1", Strand::Plus, &closer, 50).unwrap();
        assert_eq!(m.model.transcript_id, 7);
        assert_eq!(m.runner_up, Some((9, 32)));
        let other = [100, 208, 308, 408, 508, 600];
        let m = ix.find("chr1", Strand::Plus, &other, 50).unwrap();
        assert_eq!(m.model.transcript_id, 9);
        assert_eq!(m.runner_up, Some((7, 32)));

        // Identical chains differing only in annotated ends: the read's own
        // 3' end resolves it, so the read is still assigned.
        let ix2 = index(vec![
            entry(7, &[200, 300, 400, 500], 100, 600),
            entry(9, &[200, 300, 400, 500], 100, 900),
        ]);
        let short = [100, 200, 300, 400, 500, 610];
        assert_eq!(ix2.find("chr1", Strand::Plus, &short, 50).unwrap().model.transcript_id, 7);
        let long = [100, 200, 300, 400, 500, 890];
        assert_eq!(ix2.find("chr1", Strand::Plus, &long, 50).unwrap().model.transcript_id, 9);
    }

    /// Only a tie on every criterion declines.
    #[test]
    fn indistinguishable_models_decline() {
        let ix = index(vec![
            entry(7, &[200, 300, 400, 500], 100, 600),
            entry(9, &[200, 300, 400, 500], 100, 600),
        ]);
        let read = [100, 205, 305, 405, 505, 600];
        let r = ix.resolve("chr1", Strand::Plus, &read, 50);
        assert!(r.matched.is_none());
        assert_eq!(r.tie, Some((7, 9)));
        assert!(r.trace_row("test", "chr1", Strand::Plus, &read).contains("\tdeclined_tie\t7\t"));
    }

    /// A read donor sitting exactly on an annotated donor is locked: the chain
    /// may not move it 12 bp onto a different annotated donor. The same
    /// coordinate, when it is not an annotated donor, is free to shift.
    #[test]
    fn locked_annotated_site_refuses_a_shift_to_another_annotated_site() {
        let base = || index(vec![entry(7, &[200, 300, 400, 500], 100, 600)]);
        let read = [100, 212, 300, 400, 500, 600];
        let locked = with_sites(base(), &[200, 212, 400], &[300, 500]);
        let r = locked.resolve("chr1", Strand::Plus, &read, 50);
        assert!(r.matched.is_none() && r.tie.is_none());
        assert_eq!(r.locked, vec![true, true, true, true]);

        let unlocked = with_sites(base(), &[200, 400], &[300, 500]);
        let m = unlocked.find("chr1", Strand::Plus, &read, 50).unwrap();
        assert_eq!(m.model.transcript_id, 7);
        assert_eq!(m.max_offset, 12);
        assert_eq!(m.n_shifted, 1);

        // A site whose read sequence could not tell the two annotated donors
        // apart is not locked: the whole chain places the read on transcript 7.
        let ambiguous = with_sites(base(), &[200, 212, 400], &[300, 500]);
        let r = ambiguous.resolve_with_unlocked("chr1", Strand::Plus, &read, 50, 0b0001, 0);
        assert_eq!(r.locked, vec![false, true, true, true]);
        assert_eq!(r.matched.unwrap().model.transcript_id, 7);

        // With the lock off the annotated site does not bind.
        let mut off = with_sites(base(), &[200, 212, 400], &[300, 500]);
        off.lock_annotated = false;
        assert_eq!(off.find("chr1", Strand::Plus, &read, 50).unwrap().model.transcript_id, 7);
    }

    /// Roles are not interchangeable: an annotated acceptor at a coordinate
    /// does not lock a query site that is a donor.
    #[test]
    fn annotated_acceptor_does_not_lock_a_donor_query() {
        let ix = with_sites(
            index(vec![entry(7, &[200, 300, 400, 500], 100, 600)]),
            &[200, 400],
            &[212, 300, 500],
        );
        let read = [100, 212, 300, 400, 500, 600];
        let r = ix.resolve("chr1", Strand::Plus, &read, 50);
        assert_eq!(r.locked, vec![false, true, true, true]);
        assert_eq!(r.matched.unwrap().model.transcript_id, 7);
    }

    /// A novel donor the splice-site consensus settled is locked like an
    /// annotated one: the chain does not move it 12 bp onto the annotated donor.
    #[test]
    fn a_settled_novel_site_is_locked() {
        let ix = with_sites(index(vec![entry(7, &[200, 300, 400, 500], 100, 600)]), &[200, 400], &[300, 500]);
        let read = [100, 212, 300, 400, 500, 600];
        assert_eq!(ix.resolve_with_unlocked("chr1", Strand::Plus, &read, 50, 0, 0).matched.unwrap().model.transcript_id, 7);
        let r = ix.resolve_with_unlocked("chr1", Strand::Plus, &read, 50, 0, 0b0001);
        assert_eq!(r.locked, vec![true, true, true, true]);
        assert!(r.matched.is_none() && r.tie.is_none());
    }

    /// A read whose junctions all sit on annotated sites in a combination no
    /// transcript has is a NIC: every site is locked, nothing fits.
    #[test]
    fn all_annotated_novel_combination_is_never_absorbed() {
        let ix = with_sites(
            index(vec![
                entry(7, &[200, 300, 400, 500], 100, 600),
                entry(9, &[200, 330, 420, 500], 100, 600),
            ]),
            &[200, 400, 420],
            &[300, 330, 500],
        );
        // A's acceptor (300) with B's donor (420).
        let read = [100, 200, 300, 420, 500, 600];
        let r = ix.resolve("chr1", Strand::Plus, &read, 50);
        assert!(r.matched.is_none() && r.tie.is_none());
        assert!(r.trace_row("test", "chr1", Strand::Plus, &read).contains("\t1111\tdeclined_none\t"));
    }
}
