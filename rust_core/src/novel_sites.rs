//! Consensus of novel splice sites: positions with no annotated site of their role within `annotated_zone` are
//! weighed by their reads and by the strength of the splice site on the read's strand.
//! - A position that is no splice site on its strand moves onto the strongest site within `window` bp.
//! - A weaker site moves onto a stronger one holding at least `1/ratio` times its reads, unless it is itself a
//!   strong site holding `protect_fraction` of that support.
//! - Between `snap_zone` and `annotated_zone` of an annotated site, a strong site at least as strong as the
//!   annotated one is settled where it lies.

use crate::bam_reader::ParsedRead;
use crate::data_structures::SharedRefs;
use crate::genome_ref::GenomeRef;
use crate::splice_scores::SpliceScores;
use crate::types::{EdgeType, Strand};
use rustc_hash::FxHashMap;
use std::io::Write;

/// Annotated donors and acceptors per chromosome and strand, sorted.
#[derive(Debug, Default)]
pub struct AnnotatedSites {
    /// Chromosome → (donors, acceptors) on the plus and the minus strand.
    sites: FxHashMap<String, [(Vec<i64>, Vec<i64>); 2]>,
}

impl AnnotatedSites {
    /// From the intron edges of the reference; `v1` of an intron edge is its donor on both strands.
    pub fn from_reference(refs: &SharedRefs) -> Self {
        let mut position: FxHashMap<i64, (i64, &str)> = FxHashMap::default();
        for (chrom, c) in refs.location_dict.chromosomes.iter() {
            for &(p, v) in c.positions.iter() {
                position.insert(v, (p, chrom.as_str()));
            }
        }
        let mut out = AnnotatedSites::default();
        for e in refs.edge_dict.values() {
            if e.edge_type != EdgeType::Intron {
                continue;
            }
            let (Some(&(d, chrom)), Some(&(a, _))) = (position.get(&e.v1), position.get(&e.v2)) else { continue };
            out.add(chrom, e.strand, d, a);
        }
        for s in out.sites.values_mut().flat_map(|c| c.iter_mut()) {
            s.0.sort_unstable();
            s.0.dedup();
            s.1.sort_unstable();
            s.1.dedup();
        }
        out
    }

    fn add(&mut self, chrom: &str, strand: Strand, donor: i64, acceptor: i64) {
        let s = &mut self.sites.entry(chrom.to_string()).or_default()[usize::from(strand == Strand::Minus)];
        s.0.push(donor);
        s.1.push(acceptor);
    }

    /// The annotated site of this strand and role nearest to `pos` at a distance in `lo..=hi`, with the
    /// distance; the lower position on a tie.
    pub fn nearest(&self, chrom: &str, strand: Strand, donor: bool, pos: i64, lo: i64, hi: i64) -> Option<(i64, i64)> {
        let s = &self.sites.get(chrom)?[usize::from(strand == Strand::Minus)];
        let v = if donor { &s.0 } else { &s.1 };
        let mut best: Option<(i64, i64)> = None;
        for &x in &v[v.partition_point(|&x| x < pos - hi)..] {
            if x > pos + hi {
                break;
            }
            let d = (x - pos).abs();
            if d >= lo && best.is_none_or(|(_, b)| d < b) {
                best = Some((x, d));
            }
        }
        best
    }

    /// Sites from (chromosome, strand, donor, acceptor) introns, sorted as `from_reference` leaves them.
    #[cfg(test)]
    pub(crate) fn from_introns(introns: &[(&str, Strand, i64, i64)]) -> Self {
        let mut out = AnnotatedSites::default();
        for &(chrom, strand, d, a) in introns {
            out.add(chrom, strand, d, a);
        }
        for s in out.sites.values_mut().flat_map(|c| c.iter_mut()) {
            s.0.sort_unstable();
            s.1.sort_unstable();
        }
        out
    }
}

/// Strand and role of a splice position, as one key.
pub fn code(strand: Strand, donor: bool) -> u8 {
    (if strand == Strand::Minus { 2 } else { 0 }) + u8::from(!donor)
}

/// (raw position, strand-and-role code) → (consensus position, settled), per chromosome: every position the
/// consensus moves, and every settled centre mapping to itself. A settled site is kept by vertex matching.
pub type SiteMap = FxHashMap<(i64, u8), (i64, bool)>;

#[derive(Debug, Clone)]
pub struct NovelSiteParams {
    /// Farthest a position moves onto a stronger one.
    pub window: i64,
    /// Positions with an annotated site of their role this close are left to annotation snapping.
    pub annotated_zone: i64,
    /// Positions this close to an annotated site of their role always snap onto it.
    pub snap_zone: i64,
    /// A position farther than `snap_zone` from an annotated site of its role, but within `annotated_zone`, is
    /// settled with two reads when at least this strong and no weaker than that site plus `near_settle_margin`.
    pub near_settle_strength: f64,
    pub near_settle_margin: f64,
    /// A position moves only onto one holding at least `1/ratio` times its reads.
    pub ratio: f64,
    /// A site at least this strong, holding `protect_fraction` of its neighbour's reads, stays.
    pub protect_score: i8,
    pub protect_fraction: f64,
    /// Weight of one strength unit against one log2 unit of read support.
    pub score_weight: f64,
    /// Strength counted for a position that is no splice site on its strand.
    pub missing_score: f64,
    /// A position that is no splice site on its strand moves only onto a site holding this share of its reads.
    pub target_fraction: f64,
    /// A centre this strong (with two reads) is settled: vertex matching keeps it and the positions moved onto it.
    pub settle_strength: f64,
}

impl Default for NovelSiteParams {
    fn default() -> Self {
        NovelSiteParams { window: 20, annotated_zone: 25, snap_zone: 6, near_settle_strength: 11.0, near_settle_margin: 0.0,
            ratio: 0.5, protect_score: 8, protect_fraction: 0.2, score_weight: 0.25, missing_score: -6.0,
            target_fraction: 0.1, settle_strength: 6.0 }
    }
}

/// Strength of the splice site at IsoForge position `pos` on the read's strand: its score in
/// `scores`, 3 for the minor GC and AT–AC classes, and, without scores, 8 for GT–AG; None when the
/// genome carries no splice dinucleotide of the role there.
pub fn site_strength(genome: Option<&GenomeRef>, scores: Option<&SpliceScores>, chrom: &str, strand: Strand, donor: bool, pos: i64) -> Option<f64> {
    if let Some(sc) = scores.and_then(|s| s.score(chrom, strand, donor, pos)) {
        return Some(sc as f64);
    }
    let g = genome?;
    // Intron side of the site: after the exon end at a left boundary, before the exon start at a right one.
    let left = (strand == Strand::Plus) == donor;
    let raw = if left { g.fetch_seq(chrom, pos, pos + 2).ok()? } else { g.fetch_seq(chrom, pos - 3, pos - 1).ok()? };
    if raw.len() < 2 {
        return None;
    }
    let d = if strand == Strand::Minus {
        let c = |b: u8| match b { b'A' => b'T', b'C' => b'G', b'G' => b'C', b'T' => b'A', x => x };
        [c(raw[1]), c(raw[0])]
    } else {
        [raw[0], raw[1]]
    };
    match (donor, &d) {
        (true, b"GT") | (false, b"AG") => if scores.is_some() { None } else { Some(8.0) },
        (true, b"GC") | (true, b"AT") | (false, b"AC") => Some(3.0),
        _ => None,
    }
}

/// The consensus of the novel splice positions of `reads`.
pub fn novel_site_consensus(
    reads: &[ParsedRead],
    annotated: &AnnotatedSites,
    genome: Option<&GenomeRef>,
    scores: Option<&SpliceScores>,
    p: &NovelSiteParams,
) -> FxHashMap<String, SiteMap> {
    consensus_of(reads.iter().map(|r| (r.chromosome.as_str(), r.strand, r.positions.as_slice())), annotated,
        &|chrom, strand, donor, pos| site_strength(genome, scores, chrom, strand, donor, pos), genome.is_some() || scores.is_some(), p)
}

fn consensus_of<'a>(
    reads: impl Iterator<Item = (&'a str, Strand, &'a [i64])>,
    annotated: &AnnotatedSites,
    strength: &dyn Fn(&str, Strand, bool, i64) -> Option<f64>,
    weighed: bool,
    p: &NovelSiteParams,
) -> FxHashMap<String, SiteMap> {
    let mut support: FxHashMap<(&str, u8), FxHashMap<i64, f64>> = FxHashMap::default();
    // Positions between snapping reach and `annotated_zone` of an annotated site of their role: reads and that site.
    let mut near: FxHashMap<(&str, u8), FxHashMap<i64, (f64, i64)>> = FxHashMap::default();
    for (chrom, strand, positions) in reads {
        let n = positions.len();
        if n < 4 {
            continue;
        }
        for i in 1..n - 1 {
            let donor = (i - 1) % 2 == 0;
            let pos = positions[i];
            match annotated.nearest(chrom, strand, donor, pos, 0, p.annotated_zone) {
                None => *support.entry((chrom, code(strand, donor))).or_default().entry(pos).or_default() += 1.0,
                Some((site, d)) if weighed && d > p.snap_zone => {
                    near.entry((chrom, code(strand, donor))).or_default().entry(pos).or_insert((0.0, site)).0 += 1.0;
                }
                Some(_) => {}
            }
        }
    }
    let mut out: FxHashMap<String, SiteMap> = FxHashMap::default();
    // ISOFORGE_NOVEL_SITE_TRACE=<file>: every novel position with its reads, strength, score and fate.
    let mut trace = std::env::var_os("ISOFORGE_NOVEL_SITE_TRACE").and_then(|p| {
        std::fs::OpenOptions::new().create(true).append(true).open(p).ok().map(std::io::BufWriter::new)
    });
    for ((chrom, c), sites) in support {
        let strand = if c >= 2 { Strand::Minus } else { Strand::Plus };
        let donor = c % 2 == 0;
        let mut v: Vec<(i64, f64, f64, Option<f64>)> = sites
            .into_iter()
            .map(|(pos, w)| {
                let st = if weighed { strength(chrom, strand, donor, pos) } else { Some(0.0) };
                (pos, w, (1.0 + w).log2() + p.score_weight * st.unwrap_or(p.missing_score), st)
            })
            .collect();
        if v.len() < 2 {
            continue;
        }
        v.sort_unstable_by_key(|x| x.0);
        let mut order: Vec<usize> = (0..v.len()).collect();
        order.sort_by(|&a, &b| v[b].2.total_cmp(&v[a].2).then(v[b].1.total_cmp(&v[a].1)).then(v[a].0.cmp(&v[b].0)));
        let mut state = vec![0u8; v.len()]; // 0 free, 1 centre, 2 moved
        let mut target = vec![0i64; v.len()];
        let mut moved_onto_i: Vec<i64> = Vec::new();
        let map = out.entry(chrom.to_string()).or_default();
        for &i in &order {
            if state[i] != 0 {
                continue;
            }
            state[i] = 1;
            let (pi, wi, si, sti) = v[i];
            let lo = v.partition_point(|x| x.0 < pi - p.window);
            for j in lo..v.len() {
                let (pj, wj, sj, stj) = v[j];
                if pj > pi + p.window {
                    break;
                }
                if state[j] != 0 {
                    continue;
                }
                let moves = match (stj, sti) {
                    // No splice site on its strand: onto the strongest site near it that carries reads.
                    (None, Some(_)) => wi >= 2.0 && wi >= p.target_fraction * wj,
                    (None, None) => false,
                    (Some(x), _) => sj < si && wj <= p.ratio * wi
                        && !(x >= p.protect_score as f64 && wj >= p.protect_fraction * wi && wj >= 3.0),
                };
                if moves {
                    state[j] = 2;
                    target[j] = pi;
                    moved_onto_i.push(pj);
                }
            }
            let settled = wi >= 2.0 && (!weighed || sti.is_some_and(|x| x >= p.settle_strength));
            for pj in moved_onto_i.drain(..) {
                map.insert((pj, c), (pi, settled));
            }
            if settled {
                map.insert((pi, c), (pi, true));
            }
        }
        if let Some(t) = trace.as_mut() {
            for (k, &(pos, w, sc, st)) in v.iter().enumerate() {
                let _ = writeln!(t, "{chrom}\t{}\t{}\t{pos}\t{w}\t{}\t{sc:.2}\t{}\t{}",
                    if strand == Strand::Minus { '-' } else { '+' }, if donor { "donor" } else { "acceptor" },
                    st.map_or("NA".to_string(), |x| format!("{x}")), ["free", "centre", "moved"][state[k] as usize],
                    if state[k] == 2 { target[k].to_string() } else { ".".to_string() });
            }
        }
    }
    for ((chrom, c), sites) in near {
        let strand = if c >= 2 { Strand::Minus } else { Strand::Plus };
        let donor = c % 2 == 0;
        for (pos, (w, site)) in sites {
            let st = strength(chrom, strand, donor, pos);
            let annotated_st = strength(chrom, strand, donor, site).unwrap_or(p.missing_score);
            let settled = w >= 2.0
                && st.is_some_and(|x| x >= p.near_settle_strength && x >= annotated_st + p.near_settle_margin);
            if settled {
                out.entry(chrom.to_string()).or_default().insert((pos, c), (pos, true));
            }
            if let Some(t) = trace.as_mut() {
                let _ = writeln!(t, "{chrom}\t{}\t{}\t{pos}\t{w}\t{}\t{annotated_st}\t{}\t{site}",
                    if strand == Strand::Minus { '-' } else { '+' }, if donor { "donor" } else { "acceptor" },
                    st.map_or("NA".to_string(), |x| format!("{x}")), if settled { "near_settled" } else { "near_snapped" });
            }
        }
    }
    out.retain(|_, m| !m.is_empty());
    out
}

/// A read's positions with its novel splice sites moved onto their consensus (None when nothing
/// moves), and the internal sites the consensus settled (bit i = positions[i + 1]). A move that would
/// leave an exon or intron without length is skipped.
pub fn apply(positions: &[i64], strand: Strand, map: &SiteMap) -> (Option<Vec<i64>>, u64) {
    apply_keeping(positions, strand, map, 0)
}

/// As `apply`, leaving the sites in `keep` (bit i = positions[i + 1]) where they lie.
pub fn apply_keeping(positions: &[i64], strand: Strand, map: &SiteMap, keep: u64) -> (Option<Vec<i64>>, u64) {
    let n = positions.len();
    if n < 4 {
        return (None, 0);
    }
    let mut out: Option<Vec<i64>> = None;
    let mut settled = 0u64;
    for i in 1..n - 1 {
        if i - 1 < 64 && keep & (1u64 << (i - 1)) != 0 { continue; }
        let donor = (i - 1) % 2 == 0;
        let Some(&(t, settle)) = map.get(&(positions[i], code(strand, donor))) else { continue };
        let cur = out.as_deref().unwrap_or(positions);
        let (prev, next) = (cur[i - 1], cur[i + 1]);
        let ordered = if strand == Strand::Minus { prev > t && t > next } else { prev < t && t < next };
        if !ordered {
            continue;
        }
        if t != positions[i] {
            out.get_or_insert_with(|| positions.to_vec())[i] = t;
        }
        if settle && i - 1 < 64 {
            settled |= 1u64 << (i - 1);
        }
    }
    (out, settled)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn read(chrom: &str, strand: Strand, positions: &[i64]) -> (String, Strand, Vec<i64>) {
        (chrom.to_string(), strand, positions.to_vec())
    }

    fn novel_site_consensus(
        reads: &[(String, Strand, Vec<i64>)],
        annotated: &AnnotatedSites,
        scores: Option<&SpliceScores>,
        p: &NovelSiteParams,
    ) -> FxHashMap<String, SiteMap> {
        // Without a genome, a site the scores do not list counts as no splice site.
        let strength = |chrom: &str, strand: Strand, donor: bool, pos: i64| {
            scores.and_then(|s| s.score(chrom, strand, donor, pos)).map(f64::from)
        };
        consensus_of(reads.iter().map(|(c, s, v)| (c.as_str(), *s, v.as_slice())), annotated, &strength, scores.is_some(), p)
    }

    #[test]
    fn scattered_novel_site_moves_onto_its_strong_mode() {
        let mut reads: Vec<(String, Strand, Vec<i64>)> = (0..10).map(|_| read("chr1", Strand::Plus, &[1, 100, 200, 300])).collect();
        reads.extend((0..3).map(|_| read("chr1", Strand::Plus, &[1, 108, 200, 300])));
        reads.push(read("chr1", Strand::Plus, &[1, 91, 200, 300]));
        let m = novel_site_consensus(&reads, &AnnotatedSites::default(), None, &NovelSiteParams::default());
        let chr1 = &m["chr1"];
        assert_eq!(chr1.get(&(108, code(Strand::Plus, true))).map(|x| x.0), Some(100));
        assert_eq!(chr1.get(&(91, code(Strand::Plus, true))).map(|x| x.0), Some(100));
        assert_eq!(chr1.get(&(100, code(Strand::Plus, true))), Some(&(100, true)));
        // The acceptor at 200 is one position with every read: nothing to move.
        assert_eq!(chr1.get(&(200, code(Strand::Plus, false))), None);
        assert_eq!(apply(&[1, 108, 200, 300], Strand::Plus, chr1), (Some(vec![1, 100, 200, 300]), 0b1));
        // The centre itself is settled without moving.
        assert_eq!(apply(&[1, 100, 200, 300], Strand::Plus, chr1), (None, 0b1));
    }

    #[test]
    fn a_strong_alternative_site_keeps_its_reads() {
        let tsv = "chr1\t100\t+\tD\t11\nchr1\t103\t+\tD\t10\nchr1\t106\t+\tD\t-2\n";
        let scores = SpliceScores::from_reader(tsv.as_bytes()).unwrap();
        let mut reads: Vec<(String, Strand, Vec<i64>)> = (0..10).map(|_| read("chr1", Strand::Plus, &[1, 100, 200, 300])).collect();
        reads.extend((0..4).map(|_| read("chr1", Strand::Plus, &[1, 103, 200, 300])));
        reads.extend((0..4).map(|_| read("chr1", Strand::Plus, &[1, 106, 200, 300])));
        let m = novel_site_consensus(&reads, &AnnotatedSites::default(), Some(&scores), &NovelSiteParams::default());
        let chr1 = &m["chr1"];
        assert_eq!(chr1.get(&(103, code(Strand::Plus, true))), Some(&(103, true)));
        assert_eq!(chr1.get(&(106, code(Strand::Plus, true))).map(|x| x.0), Some(100));
    }

    #[test]
    fn splice_strength_can_outweigh_a_small_majority() {
        // 5 reads at a weak site, 4 at a strong one 6 bp away: the weak site moves onto the strong one.
        let tsv = "chr1\t100\t+\tD\t12\nchr1\t106\t+\tD\t-3\n";
        let scores = SpliceScores::from_reader(tsv.as_bytes()).unwrap();
        let mut reads: Vec<(String, Strand, Vec<i64>)> = (0..4).map(|_| read("chr1", Strand::Plus, &[1, 100, 200, 300])).collect();
        reads.extend((0..5).map(|_| read("chr1", Strand::Plus, &[1, 106, 200, 300])));
        let p = NovelSiteParams { ratio: 1.5, ..NovelSiteParams::default() };
        let m = novel_site_consensus(&reads, &AnnotatedSites::default(), Some(&scores), &p);
        assert_eq!(m["chr1"].get(&(106, code(Strand::Plus, true))).map(|x| x.0), Some(100));
    }

    #[test]
    fn a_site_off_its_strand_moves_onto_the_splice_site_even_with_more_reads() {
        // 9 reads at a position that is no splice site on its strand, 6 at a splice site 7 bp away.
        let tsv = "chr1\t100\t+\tD\t12\n";
        let scores = SpliceScores::from_reader(tsv.as_bytes()).unwrap();
        let mut reads: Vec<(String, Strand, Vec<i64>)> = (0..6).map(|_| read("chr1", Strand::Plus, &[1, 100, 200, 300])).collect();
        reads.extend((0..9).map(|_| read("chr1", Strand::Plus, &[1, 107, 200, 300])));
        let m = novel_site_consensus(&reads, &AnnotatedSites::default(), Some(&scores), &NovelSiteParams::default());
        assert_eq!(m["chr1"].get(&(107, code(Strand::Plus, true))).map(|x| x.0), Some(100));
    }

    #[test]
    fn positions_near_annotated_sites_are_left_to_snapping() {
        let ann = AnnotatedSites::from_introns(&[("chr1", Strand::Plus, 110, 150)]);
        let mut reads: Vec<(String, Strand, Vec<i64>)> = (0..10).map(|_| read("chr1", Strand::Plus, &[1, 100, 200, 300])).collect();
        reads.extend((0..3).map(|_| read("chr1", Strand::Plus, &[1, 104, 200, 300])));
        let m = novel_site_consensus(&reads, &ann, None, &NovelSiteParams::default());
        assert!(m.get("chr1").map_or(true, |c| !c.contains_key(&(104, code(Strand::Plus, true)))));
    }

    #[test]
    fn a_strong_site_near_an_annotated_one_is_settled_when_at_least_as_strong() {
        // Annotated donors at 110 (score 9), 300 (score 12) and 500 (score 5); read donors 10 bp from each,
        // scored 13, 11 and 10; one 5 bp from 110 inside snapping reach; and one with a single read.
        let tsv = "chr1\t100\t+\tD\t13\nchr1\t110\t+\tD\t9\nchr1\t105\t+\tD\t13\nchr1\t90\t+\tD\t13\n\
                   chr1\t290\t+\tD\t11\nchr1\t300\t+\tD\t12\nchr1\t490\t+\tD\t10\nchr1\t500\t+\tD\t5\n";
        let scores = SpliceScores::from_reader(tsv.as_bytes()).unwrap();
        let ann = AnnotatedSites::from_introns(&[("chr1", Strand::Plus, 110, 200), ("chr1", Strand::Plus, 300, 400),
                                                 ("chr1", Strand::Plus, 500, 600)]);
        let mut reads: Vec<(String, Strand, Vec<i64>)> = Vec::new();
        for donor in [100, 105, 290, 490] {
            reads.extend((0..3).map(|_| read("chr1", Strand::Plus, &[1, donor, 650, 700])));
        }
        reads.push(read("chr1", Strand::Plus, &[1, 90, 650, 700]));
        let m = novel_site_consensus(&reads, &ann, Some(&scores), &NovelSiteParams::default());
        let chr1 = &m["chr1"];
        let d = code(Strand::Plus, true);
        assert_eq!(chr1.get(&(100, d)), Some(&(100, true)));
        assert_eq!(chr1.get(&(105, d)), None, "within snapping reach");
        assert_eq!(chr1.get(&(290, d)), None, "weaker than the annotated site");
        assert_eq!(chr1.get(&(490, d)), None, "below the settling strength");
        assert_eq!(chr1.get(&(90, d)), None, "one read");
        assert_eq!(apply(&[1, 100, 650, 700], Strand::Plus, chr1), (None, 0b1));
    }

    #[test]
    fn minus_strand_roles_follow_transcript_order() {
        let mut reads: Vec<(String, Strand, Vec<i64>)> = (0..10).map(|_| read("chr1", Strand::Minus, &[400, 300, 200, 100])).collect();
        reads.extend((0..3).map(|_| read("chr1", Strand::Minus, &[400, 293, 200, 100])));
        let m = novel_site_consensus(&reads, &AnnotatedSites::default(), None, &NovelSiteParams::default());
        assert_eq!(m["chr1"].get(&(293, code(Strand::Minus, true))).map(|x| x.0), Some(300));
        assert_eq!(apply(&[400, 293, 200, 100], Strand::Minus, &m["chr1"]).0, Some(vec![400, 300, 200, 100]));
    }

    #[test]
    fn a_move_that_would_cross_the_next_site_is_skipped() {
        let mut map = SiteMap::default();
        map.insert((110, code(Strand::Plus, false)), (125, true));
        assert_eq!(apply(&[1, 100, 110, 120, 200, 300], Strand::Plus, &map), (None, 0));
    }
}
