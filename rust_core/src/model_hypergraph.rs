// The run's models as a hypergraph, as the filter's rules see them.
//
// Every model, annotated or built from reads, is a hyperedge: the chain of its vertices 5'->3' (PostDB::chains, read
// from the hyperedge tables). Between consecutive splice vertices runs an intron link, identified by the edge the run
// database stores for it, which also holds the link's two vertices. A model's intron chain is the sequence of its
// links, and two models share an intron when they share a link. `ModelHypergraph` holds each model's link chain and
// each link's vertices.

use rustc_hash::FxHashMap;

use crate::post_processing::{PostDB, TranscriptRecord};
use crate::types::Strand;

/// An intron link: the chromosome, strand and the positions of its two vertices (donor, acceptor: the last base of
/// the exon before, the first of the exon after).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Link {
    pub chromosome: String,
    pub strand: Strand,
    pub donor: i64,
    pub acceptor: i64,
}

/// The models' hyperedges in rule terms: each model's link chain 5'->3', and each link's vertices.
pub struct ModelHypergraph {
    links: FxHashMap<i64, Vec<i64>>,
    link: FxHashMap<i64, Link>,
}

/// The intron links of a record, 5'->3': the edges at the even places of its internal path.
pub fn record_links(t: &TranscriptRecord) -> impl Iterator<Item = i64> + '_ {
    t.jn_path.as_deref().unwrap_or("").split(',').step_by(2).filter_map(|s| s.parse::<i64>().ok())
}

impl ModelHypergraph {
    /// Every record's links, and each link's vertices from the link itself (the intron edge: its 5' vertex, the donor,
    /// then the acceptor).
    pub fn new(db: &PostDB) -> Self {
        let edges: FxHashMap<i64, (i64, i64, &str)> = db.edges.iter().map(|e| (e.edge_id, (e.v1, e.v2, e.strand.as_str()))).collect();
        let mut links: FxHashMap<i64, Vec<i64>> = FxHashMap::default();
        let mut link: FxHashMap<i64, Link> = FxHashMap::default();
        for t in &db.transcripts {
            let ids: Vec<i64> = record_links(t).collect();
            for &id in &ids {
                if link.contains_key(&id) { continue; }
                let Some(&(v1, v2, strand)) = edges.get(&id) else { continue };
                let (Some(&d), Some(&a)) = (db.loc_map.get(&v1), db.loc_map.get(&v2)) else { continue };
                let strand = if strand == "-" { Strand::Minus } else { Strand::Plus };
                link.insert(id, Link { chromosome: t.chromosome.clone(), strand, donor: d, acceptor: a });
            }
            links.insert(t.transcript_id, ids);
        }
        Self { links, link }
    }

    /// The model's intron links, 5'->3' (none for a single-exon model or an unknown one).
    pub fn links(&self, tid: i64) -> &[i64] {
        self.links.get(&tid).map_or(&[], |v| v.as_slice())
    }

    pub fn is_spliced(&self, tid: i64) -> bool {
        !self.links(tid).is_empty()
    }

    /// The vertices of an intron link.
    pub fn link(&self, id: i64) -> Option<&Link> {
        self.link.get(&id)
    }
}

/// Whether `shorter` is a contiguous run of `longer` (and shorter than it): one hyperedge's link chain inside another's.
pub fn contained(shorter: &[i64], longer: &[i64]) -> bool {
    !shorter.is_empty() && longer.len() > shorter.len() && longer.windows(shorter.len()).any(|w| w == shorter)
}

/// The genomic exons of a chain (exon boundaries 5'->3', two per exon), sorted.
pub fn genomic_exons(chain: &[i64]) -> Vec<(i64, i64)> {
    let mut ex: Vec<(i64, i64)> = chain.chunks(2).map(|p| (p[0].min(p[1]), p[0].max(p[1]))).collect();
    ex.sort_unstable();
    ex
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn containment_needs_a_contiguous_shorter_run() {
        assert!(contained(&[2, 3], &[1, 2, 3, 4]));
        assert!(!contained(&[2, 4], &[1, 2, 3, 4]));
        assert!(!contained(&[1, 2], &[1, 2]));
        assert!(!contained(&[], &[1, 2]));
    }

    #[test]
    fn exons_come_out_sorted_on_both_strands() {
        assert_eq!(genomic_exons(&[100, 200, 300, 400]), vec![(100, 200), (300, 400)]);
        assert_eq!(genomic_exons(&[400, 300, 200, 100]), vec![(100, 200), (300, 400)]);
    }
}
