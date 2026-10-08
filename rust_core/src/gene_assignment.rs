// Gene assignment: vertex concordance, interval-tree overlap and tiebreaking.

use rustc_hash::FxHashMap;

use crate::data_structures::*;
use crate::types::*;

/// Result of gene assignment
#[derive(Debug, Clone)]
pub enum GeneAssignment {
    /// Single gene found
    Found(GeneId),
    /// Multiple candidate genes, need tiebreaking
    Ambiguous(Vec<GeneId>),
    /// Potential fusion transcript
    Fusion,
    /// No gene found
    NotFound,
}

/// Find a gene match based on vertex concordance.
pub fn find_gene_match_on_vertex_basis(
    vertex_ids: &[VertexId],
    strand: Strand,
    vertex_2_gene: &Vertex2Gene,
) -> GeneAssignment {
    let mut gene_tally: FxHashMap<GeneId, usize> = FxHashMap::default();
    let mut n_gene_matches: Vec<usize> = Vec::with_capacity(vertex_ids.len());

    for &vertex in vertex_ids {
        if let Some(gene_set) = vertex_2_gene.get(&vertex) {
            // Filter to same strand
            let matches: Vec<GeneId> = gene_set
                .iter()
                .filter(|(_, s)| *s == strand)
                .map(|(g, _)| *g)
                .collect();

            n_gene_matches.push(matches.len());

            for gene_id in &matches {
                *gene_tally.entry(*gene_id).or_insert(0) += 1;
            }
        } else {
            n_gene_matches.push(0);
        }
    }

    // No shared splice junctions
    if gene_tally.is_empty() {
        return GeneAssignment::NotFound;
    }

    let max_per_vertex = n_gene_matches.iter().copied().max().unwrap_or(0);

    // Potential fusion: multiple genes but no shared splice sites
    if max_per_vertex <= 1 && gene_tally.len() > 1 {
        return GeneAssignment::Fusion;
    }

    // Multiple genes with overlapping vertices -> need tiebreaking
    if gene_tally.len() > 1 {
        let gene_ids: Vec<GeneId> = gene_tally.keys().copied().collect();
        return GeneAssignment::Ambiguous(gene_ids);
    }

    // Single gene: pick the one with highest tally
    let gene_id = gene_tally
        .into_iter()
        .max_by_key(|(_, count)| *count)
        .map(|(id, _)| id)
        .unwrap();

    GeneAssignment::Found(gene_id)
}

/// Search for an overlapping gene using the transcript interval tree.
pub fn search_for_overlap_with_gene(
    chromosome: &str,
    start: i64,
    end: i64,
    strand: Strand,
    _gene_tree: &GeneIntervalTree,
    transcript_tree: &TranscriptIntervalTree,
    gene_ids: Option<&[GeneId]>,
) -> Option<(GeneId, Strand)> {
    let min_start = start.min(end);
    let max_end = start.max(end);

    // Query transcript tree for overlap matches
    let matches = transcript_tree.query_overlap(chromosome, min_start, max_end);

    if matches.is_empty() {
        return None;
    }

    // Filter by gene_ids if provided
    let filtered: Vec<&TranscriptRecord> = if let Some(ids) = gene_ids {
        let id_set: rustc_hash::FxHashSet<GeneId> = ids.iter().copied().collect();
        matches.iter().filter(|m| id_set.contains(&m.gene_id)).collect()
    } else {
        matches.iter().collect()
    };

    if filtered.is_empty() {
        return None;
    }

    // Prefer same-strand matches
    let same_strand: Vec<&TranscriptRecord> = filtered
        .iter()
        .filter(|m| m.strand == strand)
        .copied()
        .collect();

    let candidates = if !same_strand.is_empty() {
        same_strand
    } else {
        filtered
    };

    // Find best match by distance to 5'/3' ends
    get_best_match(&candidates, min_start, max_end)
}

/// Best gene match by the smallest total distance of read ends to transcript ends.
fn get_best_match(
    matches: &[&TranscriptRecord],
    read_start: i64,
    read_end: i64,
) -> Option<(GeneId, Strand)> {
    let mut min_dist = i64::MAX;
    let mut best: Option<(GeneId, Strand)> = None;

    for m in matches {
        let start_dist = (m.min_pos - read_start).abs();
        let end_dist = (m.max_pos - read_end).abs();
        let dist = start_dist + end_dist;

        if dist < min_dist {
            min_dist = dist;
            best = Some((m.gene_id, m.strand));
        }
    }

    best
}

/// Full gene assignment.
pub fn assign_gene(
    vertex_ids: &[VertexId],
    strand: Strand,
    vertex_2_gene: &Vertex2Gene,
    chromosome: &str,
    start: i64,
    end: i64,
    gene_tree: &GeneIntervalTree,
    transcript_tree: &TranscriptIntervalTree,
) -> (GeneAssignment, bool) {
    let result = find_gene_match_on_vertex_basis(vertex_ids, strand, vertex_2_gene);

    match result {
        GeneAssignment::Found(id) => (GeneAssignment::Found(id), false),
        GeneAssignment::Fusion => (GeneAssignment::NotFound, true),
        GeneAssignment::Ambiguous(gene_ids) => {
            // Tiebreak by overlap distance
            match search_for_overlap_with_gene(
                chromosome,
                start,
                end,
                strand,
                gene_tree,
                transcript_tree,
                Some(&gene_ids),
            ) {
                Some((gene_id, _)) => (GeneAssignment::Found(gene_id), false),
                None => (GeneAssignment::NotFound, false),
            }
        }
        GeneAssignment::NotFound => (GeneAssignment::NotFound, false),
    }
}

/// Update vertex_2_gene for a newly assigned gene.
pub fn update_vertex_2_gene(
    gene_id: GeneId,
    vertex_ids: &[VertexId],
    strand: Strand,
    novel_buffer: &mut NovelBuffer,
) {
    for &vertex in vertex_ids {
        novel_buffer
            .vertex_2_gene
            .push((vertex, gene_id, strand));
    }
}
