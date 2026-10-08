// Edge matching: search, creation and matching of splice and transcript edges.

use crate::data_structures::*;
use crate::types::*;

/// Search for an existing edge: the reference dict first, then the concurrent
/// novel overlay.
pub fn search_for_edge(
    v1: VertexId,
    v2: VertexId,
    edge_type: EdgeType,
    edge_dict: &EdgeDict,
    novel_edge_dict: Option<&ConcurrentEdgeDict>,
) -> Option<EdgeId> {
    let key = EdgeKey {
        v1,
        v2,
        edge_type,
    };
    // Reference dict
    if let Some(e) = edge_dict.get(&key) {
        return Some(e.edge_id);
    }
    // Novel overlay
    if let Some(ned) = novel_edge_dict {
        if let Some(entry) = ned.get(&key) {
            return Some(entry.edge_id);
        }
    }
    None
}

/// Create a novel edge in the per-thread buffer and the concurrent overlay;
/// threads creating the same edge get the same ID. The overlay keeps the strand
/// of the first read in file order that uses the edge (`novel_buffer.read_idx`),
/// which the edge takes when the window's results are merged.
pub fn create_or_get_edge(
    v1: VertexId,
    v2: VertexId,
    edge_type: EdgeType,
    strand: Strand,
    counters: &Counters,
    novel_buffer: &mut NovelBuffer,
    novel_edge_dict: &ConcurrentEdgeDict,
) -> EdgeId {
    let key = EdgeKey { v1, v2, edge_type };
    let read = novel_buffer.read_idx;

    let mut is_new = false;
    let eid = {
        let mut entry = novel_edge_dict.entry(key).or_insert_with(|| {
            is_new = true;
            NovelEdgeEntry { edge_id: counters.edge.increment(), first: read, strand }
        });
        if read < entry.first {
            entry.first = read;
            entry.strand = strand;
        }
        entry.edge_id
    };

    if is_new {
        let edge = Edge {
            edge_id: eid,
            v1,
            v2,
            edge_type,
            strand,
        };
        novel_buffer.edges.push(edge);
    }

    eid
}

/// Match an edge or create it. Returns (edge_id, novelty) where novelty is 0
/// for known, 1 for novel.
pub fn match_or_create_edge(
    v1: VertexId,
    v2: VertexId,
    edge_type: EdgeType,
    strand: Strand,
    edge_dict: &EdgeDict,
    counters: &Counters,
    novel_buffer: &mut NovelBuffer,
    novel_edge_dict: &ConcurrentEdgeDict,
) -> (EdgeId, u8) {
    // Novelty follows reference membership; a reused novel edge stays novel.
    if let Some(eid) = search_for_edge(v1, v2, edge_type, edge_dict, None) {
        (eid, 0)
    } else {
        let eid = create_or_get_edge(v1, v2, edge_type, strand, counters, novel_buffer, novel_edge_dict);
        (eid, 1)
    }
}

#[cfg(test)]
mod provenance_tests {
    use super::*;
    #[test]
    fn observed_novel_edge_stays_novel_when_reused() {
        let refs = EdgeDict::default();
        let overlay = ConcurrentEdgeDict::default();
        let counters = Counters::new(0, 0, 0, 0, 0, 0);
        let mut buffer = NovelBuffer::new();
        let first = match_or_create_edge(1, 2, EdgeType::Intron, Strand::Plus,
            &refs, &counters, &mut buffer, &overlay);
        let second = match_or_create_edge(1, 2, EdgeType::Intron, Strand::Plus,
            &refs, &counters, &mut buffer, &overlay);
        assert_eq!(first, second);
        assert_eq!(second.1, 1);
        assert_eq!(buffer.edges.len(), 1);
    }

    #[test]
    fn shared_edge_takes_the_strand_of_its_first_read_in_file_order() {
        // one vertex pair used by reads of both strands, reached by the later read first
        let refs = EdgeDict::default();
        let overlay = ConcurrentEdgeDict::default();
        let counters = Counters::new(0, 0, 0, 0, 0, 0);
        let (mut late, mut early) = (NovelBuffer::new(), NovelBuffer::new());
        late.read_idx = 7;
        early.read_idx = 3;
        let a = match_or_create_edge(5, 5, EdgeType::Exon, Strand::Plus, &refs, &counters, &mut late, &overlay);
        let b = match_or_create_edge(5, 5, EdgeType::Exon, Strand::Minus, &refs, &counters, &mut early, &overlay);
        assert_eq!(a.0, b.0);
        let key = EdgeKey { v1: 5, v2: 5, edge_type: EdgeType::Exon };
        assert_eq!(overlay.get(&key).unwrap().strand, Strand::Minus);
        // a still later read of the other strand does not change it
        late.read_idx = 9;
        match_or_create_edge(5, 5, EdgeType::Exon, Strand::Plus, &refs, &counters, &mut late, &overlay);
        let entry = *overlay.get(&key).unwrap();
        assert_eq!((entry.first, entry.strand), (3, Strand::Minus));
    }
}

/// Match all splice edges (internal introns and exons between splice sites).
pub fn match_all_splice_edges(
    vertices: &[VertexId],
    strand: Strand,
    edge_dict: &EdgeDict,
    counters: &Counters,
    novel_buffer: &mut NovelBuffer,
    novel_edge_dict: &ConcurrentEdgeDict,
) -> (Vec<EdgeId>, Vec<u8>) {
    if vertices.len() < 2 {
        return (Vec::new(), Vec::new());
    }
    let mut edge_matches = Vec::with_capacity(vertices.len() - 1);
    let mut novelty = Vec::with_capacity(vertices.len() - 1);

    for i in 0..vertices.len() - 1 {
        // Alternate between intron and exon edges
        let edge_type = if i % 2 == 0 {
            EdgeType::Intron
        } else {
            EdgeType::Exon
        };

        let (edge_id, is_novel) = match_or_create_edge(
            vertices[i],
            vertices[i + 1],
            edge_type,
            strand,
            edge_dict,
            counters,
            novel_buffer,
            novel_edge_dict,
        );
        edge_matches.push(edge_id);
        novelty.push(is_novel);
    }

    (edge_matches, novelty)
}

/// Match all transcript edges, starting with an exon (monoexon case).
pub fn match_all_transcript_edges(
    vertices: &[VertexId],
    strand: Strand,
    edge_dict: &EdgeDict,
    counters: &Counters,
    novel_buffer: &mut NovelBuffer,
    novel_edge_dict: &ConcurrentEdgeDict,
) -> (Vec<EdgeId>, Vec<u8>) {
    if vertices.len() < 2 {
        return (Vec::new(), Vec::new());
    }
    let mut edge_matches = Vec::with_capacity(vertices.len() - 1);
    let mut novelty = Vec::with_capacity(vertices.len() - 1);

    for i in 0..vertices.len() - 1 {
        // Monoexon: start with exon, then alternate
        let edge_type = if i % 2 == 0 {
            EdgeType::Exon
        } else {
            EdgeType::Intron
        };

        let (edge_id, is_novel) = match_or_create_edge(
            vertices[i],
            vertices[i + 1],
            edge_type,
            strand,
            edge_dict,
            counters,
            novel_buffer,
            novel_edge_dict,
        );
        edge_matches.push(edge_id);
        novelty.push(is_novel);
    }

    (edge_matches, novelty)
}
