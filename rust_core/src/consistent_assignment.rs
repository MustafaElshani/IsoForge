//! Consistency-based attribution of truncated reads to annotated transcripts.
//!
//! A truncated (ISM) read is consistent with an annotated transcript when its
//! internal edges occur as a contiguous slice of the parent's chain and neither
//! unspliced terminal exon runs into a parent intron by more than
//! [`ISM_END_SLACK`]. A read consistent with exactly one annotated transcript
//! may take that transcript's ID and count towards its abundance; one
//! consistent with several confirms nothing. A truncated read never anchors
//! the parent's TSS or TES with its truncated end.

use crate::data_structures::SharedRefs;
use crate::types::{GeneId, Strand, TranscriptId, TranscriptModel};
use rustc_hash::FxHashMap;

/// How far an unspliced terminal exon of a truncated read may run past the
/// parent's neighbouring splice site and still be consistent.
pub const ISM_END_SLACK: i64 = 25;

/// One annotated multi-exon transcript in the form the consistency test needs.
#[derive(Debug, Clone)]
pub struct ParentChain {
    pub gene_id: GeneId,
    pub strand: Strand,
    /// Internal edges in transcript order: intron, exon, intron, …, intron.
    pub edges: Vec<i64>,
    /// Internal splice coordinates in transcript order: intron j occupies
    /// `coords[2j]` (donor) and `coords[2j + 1]` (acceptor).
    pub coords: Vec<i64>,
    /// Annotated 5' end (TSS) and 3' end (TES) in transcript orientation —
    /// `TranscriptModel` stores genomic left/right instead.
    pub start_pos: i64,
    pub end_pos: i64,
}

impl ParentChain {
    pub(crate) fn n_introns(&self) -> usize {
        (self.edges.len() + 1) / 2
    }

    /// Signed distance from `from` to `to` in the direction of transcription.
    fn downstream(&self, from: i64, to: i64) -> i64 {
        match self.strand {
            Strand::Plus => to - from,
            Strand::Minus => from - to,
        }
    }
}

/// Every annotated multi-exon transcript, keyed by id.
#[derive(Debug, Default)]
pub struct ParentIndex {
    chains: FxHashMap<TranscriptId, ParentChain>,
}

/// Outcome of testing a read against its candidate parents.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Parents {
    /// Exactly one annotated transcript is consistent with the read.
    Unique(Consistency),
    /// Several are; the read cannot tell them apart.
    Ambiguous(usize),
    /// None is, even though each contains the read's edges: a terminal exon
    /// runs into a parent intron.
    NoParent,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Consistency {
    pub parent: TranscriptId,
    /// The read's first intron is the parent's first intron.
    pub intact_5p: bool,
    /// The read's last intron is the parent's last intron.
    pub intact_3p: bool,
}

impl ParentIndex {
    pub fn from_reference(refs: &SharedRefs) -> Self {
        let mut index = Self::default();
        for tx in refs.transcript_dict.values() {
            if tx.n_exons < 2 {
                continue;
            }
            let Some((coords, _, path, strand)) = refs.hypergraph.internal(tx.transcript_id) else {
                continue;
            };
            if coords.len() != (tx.n_exons as usize - 1) * 2 {
                continue;
            }
            let (start_pos, end_pos) = match strand {
                Strand::Plus => (tx.start_pos, tx.end_pos),
                Strand::Minus => (tx.end_pos, tx.start_pos),
            };
            index.chains.insert(
                tx.transcript_id,
                ParentChain { gene_id: tx.gene_id, strand, edges: path, coords, start_pos, end_pos },
            );
        }
        index
    }

    pub fn len(&self) -> usize {
        self.chains.len()
    }

    pub fn is_empty(&self) -> bool {
        self.chains.is_empty()
    }

    pub fn get(&self, tid: TranscriptId) -> Option<&ParentChain> {
        self.chains.get(&tid)
    }
}

/// How close an unspliced read's 3' end must lie to a parent's annotated 3' end to
/// choose that parent among several containing exons.
const UNSPLICED_END_WINDOW: i64 = 50;

/// Parent of an unspliced read.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UnsplicedParent {
    /// Exactly one annotated transcript has an exon containing the read.
    Unique(Consistency),
    /// Several do; the read is credited to the one given (the lowest-id parent whose 3' end lies
    /// within `UNSPLICED_END_WINDOW` of the read's, else the lowest id) but cannot place any end.
    Ambiguous(Consistency, usize),
    /// No annotated exon contains the read.
    None,
}

/// An annotated single-exon transcript overlapping an unspliced read: its genomic
/// span and its ends in transcript orientation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MonoexonCandidate {
    pub transcript_id: TranscriptId,
    pub lo: i64,
    pub hi: i64,
    pub five: i64,
    pub three: i64,
}

/// The annotated single-exon transcript an unspliced read belongs to, as (index into `cands`, consistency, whether it
/// is the only candidate): among those containing the read, the one whose ends lie nearest the read's, otherwise the
/// largest overlap. The 3' end is intact, the 5' end up to [`UNSPLICED_END_WINDOW`] inside the annotated start.
pub fn monoexon_parent(cands: &[MonoexonCandidate], read_5p: i64, read_3p: i64, strand: Strand) -> Option<(usize, Consistency, bool)> {
    let (lo, hi) = (read_5p.min(read_3p), read_5p.max(read_3p));
    let containing: Vec<usize> = (0..cands.len())
        .filter(|&i| cands[i].lo - ISM_END_SLACK <= lo && hi <= cands[i].hi + ISM_END_SLACK)
        .collect();
    let (best, unique) = if containing.is_empty() {
        let b = (0..cands.len()).max_by_key(|&i| {
            let c = &cands[i];
            ((hi.min(c.hi) - lo.max(c.lo)).max(0), c.hi - c.lo)
        })?;
        (b, cands.len() == 1)
    } else {
        let b = *containing.iter().min_by_key(|&&i| {
            let c = &cands[i];
            ((read_5p - c.five).abs() + (read_3p - c.three).abs(), -(c.hi - c.lo))
        })?;
        (b, containing.len() == 1)
    };
    let c = &cands[best];
    let inside = match strand {
        Strand::Plus => read_5p - c.five,
        Strand::Minus => c.five - read_5p,
    };
    Some((best, Consistency { parent: c.transcript_id, intact_5p: inside <= UNSPLICED_END_WINDOW, intact_3p: true }, unique))
}

/// Exon of `parent` containing the unspliced read at genomic `lo..=hi`, as
/// (exon index in transcript order, exon count), allowing `slack` bp of alignment
/// noise at the exon boundaries.
pub fn containing_exon(parent: &ParentChain, lo: i64, hi: i64, slack: i64) -> Option<(usize, usize)> {
    let n_introns = parent.n_introns();
    for e in 0..=n_introns {
        let a = if e == 0 { parent.start_pos } else { parent.coords[2 * e - 1] };
        let b = if e == n_introns { parent.end_pos } else { parent.coords[2 * e] };
        if lo >= a.min(b) - slack && hi <= a.max(b) + slack {
            return Some((e, n_introns + 1));
        }
    }
    None
}

/// Which of `candidates` an unspliced read at `lo..=hi` (3' end at `read_3p`) is a
/// fragment of: the annotated spliced transcript with an exon containing it.
pub fn unspliced_parent(
    index: &ParentIndex,
    candidates: &[TranscriptId],
    lo: i64,
    hi: i64,
    read_3p: i64,
) -> UnsplicedParent {
    let mut hits: Vec<(TranscriptId, usize, usize)> = candidates
        .iter()
        .filter_map(|&tid| index.get(tid).and_then(|p| containing_exon(p, lo, hi, ISM_END_SLACK).map(|(e, n)| (tid, e, n))))
        .collect();
    if hits.is_empty() {
        return UnsplicedParent::None;
    }
    hits.sort_unstable();
    hits.dedup();
    let consistency = |&(tid, e, n): &(TranscriptId, usize, usize)| Consistency { parent: tid, intact_5p: e == 0, intact_3p: e + 1 == n };
    if hits.len() == 1 {
        return UnsplicedParent::Unique(consistency(&hits[0]));
    }
    let at_end: Vec<&(TranscriptId, usize, usize)> = hits
        .iter()
        .filter(|(tid, _, _)| index.get(*tid).is_some_and(|p| (p.end_pos - read_3p).abs() <= UNSPLICED_END_WINDOW))
        .collect();
    let chosen = at_end.first().copied().unwrap_or(&hits[0]);
    UnsplicedParent::Ambiguous(consistency(chosen), hits.len())
}

/// Whether a read with internal edges `read_edges` (transcript order) and
/// terminal coordinates `read_5p`/`read_3p` is consistent with `parent`.
/// Returns which of the read's ends coincide with the parent's chain ends.
pub fn consistent(parent: &ParentChain, read_edges: &[i64], read_5p: i64, read_3p: i64) -> Option<(bool, bool)> {
    let n = read_edges.len();
    if n == 0 || n > parent.edges.len() {
        return None;
    }
    let k = parent.edges.windows(n).position(|w| w == read_edges)?;
    // A slice that starts at an exon edge of the parent (odd k) is not a match.
    if k % 2 != 0 {
        return None;
    }
    let first = k / 2;
    let last = (k + n - 1) / 2;
    let intact_5p = first == 0;
    if !intact_5p {
        // The read's 5' exon must lie inside the parent exon that precedes
        // its first intron, not reach back into the parent's intron.
        let prev_acceptor = parent.coords[2 * first - 1];
        if parent.downstream(prev_acceptor, read_5p) < -ISM_END_SLACK {
            return None;
        }
    }
    let intact_3p = last + 1 == parent.n_introns();
    if !intact_3p {
        let next_donor = parent.coords[2 * (last + 1)];
        if parent.downstream(next_donor, read_3p) > ISM_END_SLACK {
            return None;
        }
    }
    Some((intact_5p, intact_3p))
}

/// Test a read against every candidate whose edge set contains the read's.
pub fn unique_parent(
    index: &ParentIndex,
    candidates: &[&TranscriptModel],
    read_edges: &[i64],
    read_5p: i64,
    read_3p: i64,
) -> Parents {
    let mut found: Option<Consistency> = None;
    let mut n = 0usize;
    for m in candidates {
        let Some(parent) = index.get(m.transcript_id) else {
            continue;
        };
        if let Some((intact_5p, intact_3p)) = consistent(parent, read_edges, read_5p, read_3p) {
            n += 1;
            if found.is_none() {
                found = Some(Consistency { parent: m.transcript_id, intact_5p, intact_3p });
            }
        }
    }
    match (n, found) {
        (1, Some(c)) => Parents::Unique(c),
        (0, _) => Parents::NoParent,
        _ => Parents::Ambiguous(n),
    }
}

/// Whether an FSM read of `short` could equally be a truncated read of a
/// longer annotated transcript that contains `short`'s chain. True when the
/// read's terminal on the nested side lies inside the longer transcript's
/// exon and more than `window` bp from `short`'s annotated end.
pub fn nested_in_longer(
    index: &ParentIndex,
    short: TranscriptId,
    short_n_exons: i32,
    candidates: &[&TranscriptModel],
    read_5p: i64,
    read_3p: i64,
    window: i64,
) -> bool {
    let Some(s) = index.get(short) else {
        return false;
    };
    for m in candidates {
        if m.n_exons <= short_n_exons || m.transcript_id == short {
            continue;
        }
        let Some(l) = index.get(m.transcript_id) else {
            continue;
        };
        let n = s.edges.len();
        let Some(k) = l.edges.windows(n).position(|w| w == s.edges.as_slice()) else {
            continue;
        };
        let first = k / 2;
        let last = (k + n - 1) / 2;
        if first > 0 {
            let prev_acceptor = l.coords[2 * first - 1];
            let inside_longer_exon = l.downstream(prev_acceptor, read_5p) >= -ISM_END_SLACK;
            let far_from_short_start = l.downstream(s.start_pos, read_5p).abs() > window;
            if inside_longer_exon && far_from_short_start {
                return true;
            }
        }
        if last + 1 < l.n_introns() {
            let next_donor = l.coords[2 * (last + 1)];
            let inside_longer_exon = l.downstream(next_donor, read_3p) <= ISM_END_SLACK;
            let far_from_short_end = l.downstream(s.end_pos, read_3p).abs() > window;
            if inside_longer_exon && far_from_short_end {
                return true;
            }
        }
    }
    false
}

#[cfg(test)]
mod tests {
    fn plus_parent(tid: TranscriptId) -> (TranscriptId, ParentChain) {
        // exons 100-200, 300-400, 500-600 on the plus strand
        (tid, ParentChain { gene_id: 7, strand: Strand::Plus, edges: vec![1, 2, 3], coords: vec![201, 299, 401, 499], start_pos: 100, end_pos: 600 })
    }

    #[test]
    fn monoexon_reads_find_their_transcript() {
        let cand = |id: TranscriptId, lo: i64, hi: i64, strand: Strand| {
            let (five, three) = if strand == Strand::Minus { (hi, lo) } else { (lo, hi) };
            MonoexonCandidate { transcript_id: id, lo, hi, five, three }
        };
        for strand in [Strand::Plus, Strand::Minus] {
            // A short transcript nested in a longer one: a read spanning the short one goes to it.
            let c = [cand(1, 1000, 1600, strand), cand(2, 1000, 5000, strand)];
            let (r5, r3) = if strand == Strand::Plus { (1003, 1598) } else { (1598, 1003) };
            let (i, k, unique) = monoexon_parent(&c, r5, r3, strand).unwrap();
            assert_eq!((c[i].transcript_id, k.intact_5p, k.intact_3p, unique), (1, true, true, false), "{strand:?}");
            // A read inside only the long one, starting 2 kb in: a 5' fragment of it.
            let (r5, r3) = if strand == Strand::Plus { (3000, 5000) } else { (3000, 1000) };
            let (i, k, unique) = monoexon_parent(&c, r5, r3, strand).unwrap();
            assert_eq!((c[i].transcript_id, k.intact_5p, k.intact_3p, unique), (2, false, true, true), "{strand:?}");
            // Starting upstream of the annotated start leaves the 5' end intact.
            let (r5, r3) = if strand == Strand::Plus { (900, 4000) } else { (5100, 2000) };
            let (i, k, _) = monoexon_parent(&c, r5, r3, strand).unwrap();
            assert_eq!((c[i].transcript_id, k.intact_5p), (2, true), "{strand:?}");
        }
        // Partial overlaps pick the larger overlap on either strand.
        let c = [cand(1, 1000, 2000, Strand::Minus), cand(2, 1800, 6000, Strand::Minus)];
        let (i, _, _) = monoexon_parent(&c, 1900, 900, Strand::Minus).unwrap();
        assert_eq!(c[i].transcript_id, 1);
    }

    #[test]
    fn unspliced_reads_find_their_exon() {
        let (_, p) = plus_parent(1);
        assert_eq!(containing_exon(&p, 120, 180, 5), Some((0, 3)));
        assert_eq!(containing_exon(&p, 310, 395, 5), Some((1, 3)));
        assert_eq!(containing_exon(&p, 520, 603, 5), Some((2, 3)));
        assert_eq!(containing_exon(&p, 180, 320, 5), None);
        let mut index = ParentIndex::default();
        let (tid, chain) = plus_parent(1);
        index.chains.insert(tid, chain);
        assert_eq!(unspliced_parent(&index, &[1], 520, 600, 600), UnsplicedParent::Unique(Consistency { parent: 1, intact_5p: false, intact_3p: true }));
        assert_eq!(unspliced_parent(&index, &[1], 180, 320, 320), UnsplicedParent::None);
        // A second parent sharing the last exon but ending later: the read's 3' end picks the first.
        let (tid2, mut chain2) = plus_parent(2);
        chain2.end_pos = 900;
        index.chains.insert(tid2, chain2);
        assert_eq!(unspliced_parent(&index, &[1, 2], 520, 600, 600), UnsplicedParent::Ambiguous(Consistency { parent: 1, intact_5p: false, intact_3p: true }, 2));
    }

    use super::*;

    /// Four exons on the plus strand: exons [100,200] [300,400] [500,600] [700,800].
    fn plus() -> ParentChain {
        ParentChain {
            gene_id: 1,
            strand: Strand::Plus,
            edges: vec![10, 11, 12, 13, 14],
            coords: vec![200, 300, 400, 500, 600, 700],
            start_pos: 100,
            end_pos: 800,
        }
    }

    /// The same transcript on the minus strand, coordinates in transcript order.
    fn minus() -> ParentChain {
        ParentChain {
            gene_id: 1,
            strand: Strand::Minus,
            edges: vec![10, 11, 12, 13, 14],
            coords: vec![700, 600, 500, 400, 300, 200],
            start_pos: 800,
            end_pos: 100,
        }
    }

    #[test]
    fn five_prime_truncation_inside_the_exon_is_consistent() {
        assert_eq!(consistent(&plus(), &[12, 13, 14], 350, 800), Some((false, true)));
        // Up to the slack past the acceptor is still an alignment-soft end…
        assert_eq!(consistent(&plus(), &[12, 13, 14], 300 - ISM_END_SLACK, 800), Some((false, true)));
        // …beyond it the read reaches into the parent's intron.
        assert_eq!(consistent(&plus(), &[12, 13, 14], 250, 800), None);
        // Minus strand: the read's 5' end is the higher coordinate.
        assert_eq!(consistent(&minus(), &[12, 13, 14], 450, 100), Some((false, true)));
        assert_eq!(consistent(&minus(), &[12, 13, 14], 640, 100), None);
    }

    #[test]
    fn three_prime_truncation_mirrors_it() {
        assert_eq!(consistent(&plus(), &[10, 11, 12], 100, 550), Some((true, false)));
        assert_eq!(consistent(&plus(), &[10, 11, 12], 100, 640), None);
        assert_eq!(consistent(&minus(), &[10, 11, 12], 800, 350), Some((true, false)));
        assert_eq!(consistent(&minus(), &[10, 11, 12], 800, 260), None);
        // Truncated at both ends.
        assert_eq!(consistent(&plus(), &[12], 350, 550), Some((false, false)));
    }

    #[test]
    fn non_subchains_never_match() {
        assert_eq!(consistent(&plus(), &[10, 11, 14], 100, 800), None);
        assert_eq!(consistent(&plus(), &[11, 12, 13], 100, 800), None);
        assert_eq!(consistent(&plus(), &[], 100, 800), None);
    }

    fn model(id: i64, n_exons: i32) -> TranscriptModel {
        TranscriptModel {
            transcript_id: id,
            gene_id: 1,
            start_exon: 0,
            jn_path: None,
            end_exon: 0,
            start_vertex: 0,
            end_vertex: 0,
            n_exons,
            chromosome: "chr1".into(),
            start_pos: 0,
            end_pos: 0,
        }
    }

    #[test]
    fn unique_and_ambiguous_parents() {
        let mut ix = ParentIndex::default();
        ix.chains.insert(7, plus());
        // A second parent sharing exons 2-4 but with a different first exon.
        let mut other = plus();
        other.edges = vec![20, 11, 12, 13, 14];
        other.coords = vec![150, 300, 400, 500, 600, 700];
        ix.chains.insert(9, other);
        let m7 = model(7, 4);
        let m9 = model(9, 4);
        // Exons 2-4 fit both parents: ambiguous.
        assert_eq!(unique_parent(&ix, &[&m7, &m9], &[12, 13, 14], 350, 800), Parents::Ambiguous(2));
        // Exons 1-3 fit only parent 7.
        assert_eq!(
            unique_parent(&ix, &[&m7, &m9], &[10, 11, 12], 100, 550),
            Parents::Unique(Consistency { parent: 7, intact_5p: true, intact_3p: false })
        );
        // Edges fit, but the read runs 60 bp into the parent's intron.
        assert_eq!(unique_parent(&ix, &[&m7], &[12, 13, 14], 240, 800), Parents::NoParent);
    }

    #[test]
    fn nested_short_isoform_is_ambiguous_only_when_the_read_looks_like_the_long_one() {
        let mut ix = ParentIndex::default();
        ix.chains.insert(7, plus());
        // S: exons 2-4 of L with its own TSS 60 bp inside L's exon 2.
        ix.chains.insert(
            3,
            ParentChain {
                gene_id: 1,
                strand: Strand::Plus,
                edges: vec![12, 13, 14],
                coords: vec![400, 500, 600, 700],
                start_pos: 360,
                end_pos: 800,
            },
        );
        let l = model(7, 4);
        // Starts 30 bp from S's TSS: an S read.
        assert!(!nested_in_longer(&ix, 3, 3, &[&l], 390, 800, 100));
        // Starts 160 bp past S's TSS, inside L's exon 2: a fragment of L.
        assert!(nested_in_longer(&ix, 3, 3, &[&l], 520, 800, 100));
        // A longer exon 2 for L.
        let mut long = plus();
        long.coords = vec![200, 300, 900, 1000, 1100, 1200];
        ix.chains.insert(7, long);
        ix.get(3).unwrap();
        ix.chains.get_mut(&3).unwrap().coords = vec![900, 1000, 1100, 1200];
        assert!(nested_in_longer(&ix, 3, 3, &[&l], 600, 1300, 100));
        // Starting upstream of L's acceptor by more than the slack is not inside L's exon.
        assert!(!nested_in_longer(&ix, 3, 3, &[&l], 250, 1300, 100));
        // Without a longer candidate nothing is nested.
        assert!(!nested_in_longer(&ix, 3, 3, &[], 600, 1300, 100));
    }
}
