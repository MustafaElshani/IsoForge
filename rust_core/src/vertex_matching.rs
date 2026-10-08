// Vertex matching: exact, nearby and extended splice-vertex search,
// permissive end matching, and novel vertex creation.

use crate::data_structures::*;
use crate::transcript_utils::compute_delta;
use crate::types::*;

/// Result of a vertex match attempt
#[derive(Debug, Clone)]
pub struct VertexMatch {
    pub vertex_id: VertexId,
    pub delta: Option<i64>,
    pub is_novel: bool,
}

/// Every novel vertex is visible (callers outside the windowed annotation).
pub const ALL_VISIBLE: VertexId = VertexId::MAX;

/// Search for a vertex at an exact position: the reference dict first, then
/// the concurrent novel overlay. Novel vertices above `visible_max` are not
/// seen.
pub fn search_vertex_at_pos(
    location_dict: &LocationDict,
    chromosome: &str,
    position: i64,
    novel_location_dict: Option<&ConcurrentLocationDict>,
    visible_max: VertexId,
) -> Option<VertexId> {
    // Reference dict
    if let Some(vid) = location_dict.get_vertex(chromosome, position) {
        return Some(vid);
    }
    // Novel overlay
    if let Some(nld) = novel_location_dict {
        if let Some(entry) = nld.get(&(chromosome.to_string(), position)) {
            if *entry <= visible_max {
                return Some(*entry);
            }
        }
    }
    None
}

/// Search for a splice vertex within ±max_dist bp, first in the reference
/// annotation, then in the novel overlay.
///
/// Returns the vertex ID of the nearest match, or None if nothing is within
/// max_dist bp.
pub fn search_nearby_splice_vertex(
    location_dict: &LocationDict,
    chromosome: &str,
    position: i64,
    max_dist: i64,
    novel_location_dict: &ConcurrentLocationDict,
    visible_max: VertexId,
) -> Option<VertexId> {
    search_nearby_splice_vertex_shift(location_dict, chromosome, position, max_dist, novel_location_dict, visible_max).map(|(v, _)| v)
}

/// As `search_nearby_splice_vertex`, also returning how far the read's
/// position was moved onto the chosen vertex.
pub fn search_nearby_splice_vertex_shift(
    location_dict: &LocationDict,
    chromosome: &str,
    position: i64,
    max_dist: i64,
    novel_location_dict: &ConcurrentLocationDict,
    visible_max: VertexId,
) -> Option<(VertexId, i64)> {
    search_nearby_splice_vertex_shift_ref(location_dict, chromosome, position, max_dist, novel_location_dict, visible_max).map(|(v, s, _)| (v, s))
}

/// As `search_nearby_splice_vertex_shift`, also reporting whether the hit is
/// a reference vertex (true) or one from the novel overlay (false).
pub fn search_nearby_splice_vertex_shift_ref(
    location_dict: &LocationDict,
    chromosome: &str,
    position: i64,
    max_dist: i64,
    novel_location_dict: &ConcurrentLocationDict,
    visible_max: VertexId,
) -> Option<(VertexId, i64, bool)> {
    // 1. Check reference annotation: nearest within ±max_dist
    if let Some((pos, vid)) = location_dict.find_nearest(
        chromosome,
        position,
        position - max_dist,
        position + max_dist,
    ) {
        return Some((vid, (pos - position).abs(), true));
    }

    // 2. Novel overlay: scan ±max_dist
    let chrom = chromosome.to_string();
    let mut best: Option<(VertexId, i64)> = None; // (vertex_id, abs_distance)

    for offset in 1..=max_dist {
        // Check position - offset
        if let Some(entry) = novel_location_dict.get(&(chrom.clone(), position - offset)).filter(|e| **e <= visible_max) {
            match best {
                None => best = Some((*entry, offset)),
                Some((_, best_dist)) if offset < best_dist => best = Some((*entry, offset)),
                _ => {}
            }
        }
        // Check position + offset
        if let Some(entry) = novel_location_dict.get(&(chrom.clone(), position + offset)).filter(|e| **e <= visible_max) {
            match best {
                None => best = Some((*entry, offset)),
                Some((_, best_dist)) if offset < best_dist => best = Some((*entry, offset)),
                _ => {}
            }
        }
        // Stop at the first distance with a hit
        if best.is_some() {
            break;
        }
    }

    best.map(|(v, s)| (v, s, false))
}

/// Extended reference-only splice vertex search.
///
/// Searches only the reference annotation in (max_splice_dist,
/// max_extended_snap] bp from the position. A position with short-read
/// support is kept (returns None); otherwise the read snaps to the nearest
/// reference vertex.
pub fn search_extended_ref_splice_vertex(
    location_dict: &LocationDict,
    chromosome: &str,
    position: i64,
    max_splice_dist: i64,
    max_extended_snap: i64,
    sr_splice_positions: Option<&rustc_hash::FxHashSet<(String, i64)>>,
) -> Option<VertexId> {
    search_extended_ref_splice_vertex_shift(location_dict, chromosome, position, max_splice_dist, max_extended_snap, sr_splice_positions).map(|(v, _)| v)
}

/// As `search_extended_ref_splice_vertex`, also returning the snap distance.
pub fn search_extended_ref_splice_vertex_shift(
    location_dict: &LocationDict,
    chromosome: &str,
    position: i64,
    max_splice_dist: i64,
    max_extended_snap: i64,
    sr_splice_positions: Option<&rustc_hash::FxHashSet<(String, i64)>>,
) -> Option<(VertexId, i64)> {
    // Only the extended zone (max_splice_dist, max_extended_snap]
    let inner = max_splice_dist + 1;
    let outer = max_extended_snap;

    // Binary search reference dict for nearest vertex in extended zone
    if let Some((_ref_pos, vid)) = location_dict.find_nearest(
        chromosome,
        position,
        position - outer,
        position + outer,
    ) {
        let dist = (position - _ref_pos).abs();
        // Only snap beyond max_splice_dist
        if dist >= inner && dist <= outer {
            // Positions with short-read support are kept
            if let Some(sr_positions) = sr_splice_positions {
                let key = (chromosome.to_string(), position);
                if sr_positions.contains(&key) {
                    return None;
                }
            }
            // Snap to the reference vertex
            return Some((vid, dist));
        }
    }
    None
}

/// Create a novel vertex and return its ID, recording it in the per-thread
/// buffer and the concurrent overlay; threads creating a vertex at the same
/// position get the same ID.
pub fn create_or_get_vertex(
    chromosome: &str,
    position: i64,
    counters: &Counters,
    config: &RunConfig,
    novel_buffer: &mut NovelBuffer,
    novel_location_dict: &ConcurrentLocationDict,
) -> (VertexId, bool) {
    let key = (chromosome.to_string(), position);

    // Check-or-insert into the concurrent dict; an existing ID is reused.
    let mut is_new = false;
    let vid = *novel_location_dict.entry(key).or_insert_with(|| {
        is_new = true;
        counters.vertex.increment()
    });

    if is_new {
        // Newly created: record it in this thread's buffer
        let location = Location {
            location_id: vid,
            genome_build: config.build.clone(),
            chromosome: chromosome.to_string(),
            position,
        };
        novel_buffer.locations.push(location);
    }

    (vid, is_new)
}

/// Create a vertex without concurrent deduplication.
pub fn create_vertex(
    chromosome: &str,
    position: i64,
    counters: &Counters,
    config: &RunConfig,
) -> (VertexId, Location) {
    let new_id = counters.vertex.increment();
    let location = Location {
        location_id: new_id,
        genome_build: config.build.clone(),
        chromosome: chromosome.to_string(),
        position,
    };
    (new_id, location)
}

/// Permissive vertex search: nearest annotated vertex within the cutoff,
/// checking the priority direction first at each distance. Only the
/// read-only annotation dict is searched.
pub fn permissive_vertex_search(
    location_dict: &LocationDict,
    chromosome: &str,
    position: i64,
    strand: Strand,
    sj_pos: i64,
    pos_type: &str,
    config: &RunConfig,
) -> (Option<VertexId>, Option<i64>) {
    // Try strict match first
    if let Some(vid) = location_dict.get_vertex(chromosome, position) {
        return (Some(vid), Some(0));
    }

    let max_dist = match pos_type {
        "start" => config.cutoff_5p,
        "end" => config.cutoff_3p,
        _ => return (None, None),
    };

    // Search window and priority direction
    let (window_start, window_end, direction_priority) =
        if (strand == Strand::Plus && pos_type == "start")
            || (strand == Strand::Minus && pos_type == "end")
        {
            (position - max_dist, sj_pos, -1i32)
        } else {
            (sj_pos, position + max_dist, 1i32)
        };

    // Use binary search on sorted annotation positions
    match location_dict.find_nearest_with_priority(
        chromosome,
        position,
        window_start,
        window_end,
        direction_priority,
    ) {
        Some((found_pos, vid)) => {
            let dist = compute_delta(found_pos, position, strand);
            (Some(vid), Some(dist))
        }
        None => (None, None),
    }
}

/// Permissive match with gene priority: a start/end vertex of the same gene
/// first, then a general permissive match.
pub fn permissive_match_with_gene_priority(
    location_dict: &LocationDict,
    chromosome: &str,
    position: i64,
    strand: Strand,
    sj_pos: i64,
    pos_type: &str,
    gene_id: GeneId,
    gene_endpoints: &GeneEndpointDict,
    config: &RunConfig,
) -> (Option<VertexId>, Option<i64>, bool) {
    // gene_in_match: whether the matched vertex belongs to this gene
    // Try exact match first
    if let Some(vid) = location_dict.get_vertex(chromosome, position) {
        if let Some(gene_locs) = gene_endpoints.get(&gene_id) {
            if gene_locs.contains_key(&position) {
                return (Some(vid), Some(0), true);
            }
        }
        return (Some(vid), Some(0), false);
    }

    // Try gene-specific matching
    if let Some(gene_locs) = gene_endpoints.get(&gene_id) {
        let max_dist = match pos_type {
            "start" => config.cutoff_5p,
            "end" => config.cutoff_3p,
            _ => return (None, None, false),
        };

        let (search_window_start, search_window_end) =
            if (strand == Strand::Plus && pos_type == "start")
                || (strand == Strand::Minus && pos_type == "end")
            {
                (position - max_dist, sj_pos)
            } else {
                (sj_pos, position + max_dist)
            };

        let mut min_abs_dist = max_dist + 1;
        let mut best_dist = None;
        let mut closest_vertex = None;

        for (&known_pos, &vertex_id) in gene_locs {
            if known_pos < search_window_start || known_pos > search_window_end {
                continue;
            }
            let curr_dist = compute_delta(known_pos, position, strand);
            if curr_dist.abs() < min_abs_dist {
                best_dist = Some(curr_dist);
                min_abs_dist = curr_dist.abs();
                closest_vertex = Some(vertex_id);
            }
        }

        if min_abs_dist <= max_dist {
            return (closest_vertex, best_dist, true);
        }
    }

    // Fall back to general permissive match (annotation-only)
    let (vid, dist) = permissive_vertex_search(
        location_dict,
        chromosome,
        position,
        strand,
        sj_pos,
        pos_type,
        config,
    );
    (vid, dist, false)
}

/// How one internal splice position resolves against the vertices a read may see.
pub enum SiteResolution {
    /// An existing vertex takes it: (vertex, novelty flag, bases the position moved).
    Vertex(VertexId, u8, i64),
    /// No vertex takes it: one is made (or found) at this position; bases the position moved.
    Create(i64, i64),
}

/// Resolve the `k`-th internal splice position of a read (see `match_splice_vertices`), without creating anything.
pub fn resolve_splice_site(
    location_dict: &LocationDict,
    chromosome: &str,
    position: i64,
    k: usize,
    config: &RunConfig,
    novel_location_dict: &ConcurrentLocationDict,
    sr_splice_positions: Option<&rustc_hash::FxHashSet<(String, i64)>>,
    consensus: Option<&rustc_hash::FxHashMap<i64, (i64, u32)>>,
    settled_sites: u64,
    visible_max: VertexId,
) -> SiteResolution {
    // 1. Exact match: check reference dict + novel overlay
    if let Some(vid) = search_vertex_at_pos(location_dict, chromosome, position, Some(novel_location_dict), visible_max) {
        return SiteResolution::Vertex(vid, u8::from(location_dict.get_vertex(chromosome, position) != Some(vid)), 0);
    }

    // A novel site settled by the splice-site consensus keeps its position.
    if k < 64 && settled_sites & (1u64 << k) != 0 {
        return SiteResolution::Create(position, 0);
    }

    // 2. Nearby match: snap to an existing vertex within ±max_splice_dist
    //    (reference and novel overlay). Any vertex takes the site, whatever its role.
    if config.max_splice_dist > 0 {
        if let Some((vid, shift, is_reference)) = search_nearby_splice_vertex_shift_ref(
            location_dict,
            chromosome,
            position,
            config.max_splice_dist,
            novel_location_dict,
            visible_max,
        ) {
            // A novel vertex off this window's pileup mode does not capture a read
            // on the mode when the mode has >= 3 reads; the mode vertex takes it instead.
            let strong_mode = !is_reference && shift > 0
                && consensus.and_then(|c| c.get(&position)).is_some_and(|&(m, n)| m == position && n >= 3);
            if !strong_mode {
                return SiteResolution::Vertex(vid, u8::from(!is_reference), shift);
            }
        }
    }

    // 3. Extended reference-only snap, skipped for positions with short-read
    //    support
    if config.max_extended_snap > config.max_splice_dist {
        if let Some((vid, shift)) = search_extended_ref_splice_vertex_shift(
            location_dict,
            chromosome,
            position,
            config.max_splice_dist,
            config.max_extended_snap,
            sr_splice_positions,
        ) {
            return SiteResolution::Vertex(vid, 0, shift);
        }
    }

    // 4. A novel vertex at the window's consensus position for this site when one was computed
    let target = consensus.and_then(|c| c.get(&position)).map(|&(m, _)| m).unwrap_or(position);
    SiteResolution::Create(target, (target - position).abs())
}

/// The positions at which a read's internal splice sites need a novel vertex, seen against the vertices visible
/// before its window (`visible_max`): (position, whether the splice-site consensus settled it there).
pub fn proposed_splice_vertices(
    location_dict: &LocationDict,
    chromosome: &str,
    positions: &[i64],
    config: &RunConfig,
    novel_location_dict: &ConcurrentLocationDict,
    sr_splice_positions: Option<&rustc_hash::FxHashSet<(String, i64)>>,
    consensus: Option<&rustc_hash::FxHashMap<i64, (i64, u32)>>,
    settled_sites: u64,
    visible_max: VertexId,
    out: &mut Vec<(i64, bool)>,
) {
    if positions.len() < 4 {
        return;
    }
    for i in 1..positions.len() - 1 {
        if let SiteResolution::Create(target, _) = resolve_splice_site(location_dict, chromosome, positions[i], i - 1, config,
            novel_location_dict, sr_splice_positions, consensus, settled_sites, visible_max) {
            out.push((target, i - 1 < 64 && settled_sites & (1u64 << (i - 1)) != 0));
        }
    }
}

/// Create a window's novel splice vertices from its reads' proposals (`proposed_splice_vertices`, with each one's
/// chromosome): by chromosome and position, each position once, a settled site first among equals. Without a window
/// consensus (`thin`, the snapping distance) an unsettled position within that distance of the vertex made before it
/// joins that vertex instead of becoming one.
pub fn create_window_vertices(
    proposed: &mut Vec<(&str, i64, bool)>,
    thin: Option<i64>,
    counters: &Counters,
    config: &RunConfig,
    novel_buffer: &mut NovelBuffer,
    novel_location_dict: &ConcurrentLocationDict,
) {
    proposed.sort_unstable_by(|a, b| (a.0, a.1, !a.2).cmp(&(b.0, b.1, !b.2)));
    proposed.dedup_by(|b, a| a.0 == b.0 && a.1 == b.1);
    let mut last: Option<(&str, i64)> = None;
    for &(chrom, pos, forced) in proposed.iter() {
        if !forced && thin.is_some_and(|dist| last.is_some_and(|(c, p)| c == chrom && pos - p <= dist)) {
            continue;
        }
        create_or_get_vertex(chrom, pos, counters, config, novel_buffer, novel_location_dict);
        last = Some((chrom, pos));
    }
}

/// Match internal splice vertices.
///
/// Each internal position is matched exactly (reference, then novel overlay),
/// then snapped to a vertex within ±max_splice_dist, then snapped to a
/// reference vertex in the extended zone unless short reads support it;
/// otherwise a novel vertex is created in the concurrent overlay. Novel
/// vertices above `visible_max` are not seen (`search_vertex_at_pos`).
///
/// Also returns the largest distance any internal position was moved (0 when
/// every site matched exactly or was created novel).
pub fn match_splice_vertices(
    location_dict: &LocationDict,
    chromosome: &str,
    positions: &[i64],
    _strand: Strand,
    counters: &Counters,
    config: &RunConfig,
    novel_buffer: &mut NovelBuffer,
    novel_location_dict: &ConcurrentLocationDict,
    sr_splice_positions: Option<&rustc_hash::FxHashSet<(String, i64)>>,
    consensus: Option<&rustc_hash::FxHashMap<i64, (i64, u32)>>,
    settled_sites: u64,
    visible_max: VertexId,
) -> (Vec<VertexId>, Vec<u8>, i32) {
    let mut vertex_matches = Vec::with_capacity(positions.len() - 2);
    let mut novelty = Vec::with_capacity(positions.len() - 2);
    let mut max_shift: i64 = 0;

    // Skip first and last positions (those are transcript ends)
    for i in 1..positions.len() - 1 {
        match resolve_splice_site(location_dict, chromosome, positions[i], i - 1, config, novel_location_dict,
                                  sr_splice_positions, consensus, settled_sites, visible_max) {
            SiteResolution::Vertex(vid, flag, shift) => {
                max_shift = max_shift.max(shift);
                vertex_matches.push(vid);
                novelty.push(flag);
            }
            SiteResolution::Create(target, shift) => {
                // concurrent creators get the same ID
                max_shift = max_shift.max(shift);
                let (vid, _is_new) = create_or_get_vertex(chromosome, target, counters, config, novel_buffer, novel_location_dict);
                vertex_matches.push(vid);
                novelty.push(1);
            }
        }
    }

    (vertex_matches, novelty, max_shift.min(i32::MAX as i64) as i32)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn observed_vertices_remain_novel_on_exact_and_nearby_reuse() {
        let locations = LocationDict::new();
        let counters = Counters::new(0, 0, 0, 0, 0, 0);
        let config = RunConfig::default();
        let mut novel = NovelBuffer::new();
        let overlay = ConcurrentLocationDict::new();
        for positions in [[1, 100, 201, 300], [1, 100, 201, 300], [1, 102, 203, 300]] {
            let (_, flags, _) = match_splice_vertices(&locations, "chr1", &positions, Strand::Plus,
                &counters, &config, &mut novel, &overlay, None, None, 0, ALL_VISIBLE);
            assert_eq!(flags, vec![1, 1]);
        }
        assert_eq!(novel.locations.len(), 2);
    }

    #[test]
    fn test_acceptor_snaps_to_annotated_vertex_nearby() {
        let mut locations = LocationDict::new();
        locations.insert("chr12".to_string(), 6_967_892, 1);
        locations.insert("chr12".to_string(), 6_967_248, 2);
        locations.finalize();
        let positions = [6_967_991, 6_967_892, 6_967_238, 6_967_139];
        let counters = Counters::new(0, 0, 2, 0, 0, 0);
        let config = RunConfig::default();
        let mut novel = NovelBuffer::new();
        let overlay = ConcurrentLocationDict::new();
        let (vertices, _, _) = match_splice_vertices(
            &locations, "chr12", &positions, Strand::Minus, &counters, &config,
            &mut novel, &overlay, None, None, 0, ALL_VISIBLE,
        );
        assert_eq!(vertices, vec![1, 2]);
        assert_eq!(positions, [6_967_991, 6_967_892, 6_967_238, 6_967_139]);
        assert!(novel.locations.is_empty());
    }

    /// One window as the engine runs it: the reads' proposals against the vertices from before the window, the
    /// window's vertices made in coordinate order, then every read matched against that set. Returns each read's
    /// vertices and the positions of the vertices made.
    fn window(locations: &LocationDict, reads: &[[i64; 4]], order: &[usize], settled: &[u64]) -> (Vec<Vec<VertexId>>, Vec<i64>) {
        let config = RunConfig::default();
        let counters = Counters::new(0, 0, 1, 0, 0, 0);
        let overlay = ConcurrentLocationDict::new();
        let before = counters.vertex.get();
        let mut proposed: Vec<(&str, i64, bool)> = Vec::new();
        let mut sites = Vec::new();
        for &i in order {
            sites.clear();
            proposed_splice_vertices(locations, "chr1", &reads[i], &config, &overlay, None, None, settled[i], before, &mut sites);
            proposed.extend(sites.iter().map(|&(p, forced)| ("chr1", p, forced)));
        }
        let mut made = NovelBuffer::new();
        create_window_vertices(&mut proposed, Some(config.max_splice_dist), &counters, &config, &mut made, &overlay);
        let visible = counters.vertex.get();
        let mut later = NovelBuffer::new();
        let mut vertices = vec![Vec::new(); reads.len()];
        for &i in order {
            vertices[i] = match_splice_vertices(locations, "chr1", &reads[i], Strand::Plus, &counters, &config, &mut later,
                                                &overlay, None, None, settled[i], visible).0;
        }
        assert!(later.locations.is_empty(), "matching a window's reads makes no splice vertex of its own");
        (vertices, made.locations.iter().map(|l| l.position).collect())
    }

    #[test]
    fn a_window_resolves_its_reads_the_same_in_any_order() {
        let mut locations = LocationDict::new();
        locations.insert("chr1".to_string(), 5_000, 1);
        locations.finalize();
        // novel sites a few bases apart, and one beside an annotated vertex
        let reads = [[1, 103, 300, 400], [1, 100, 302, 400], [1, 101, 300, 400], [1, 100, 5_003, 6_000]];
        let unsettled = [0u64; 4];
        let forward = window(&locations, &reads, &[0, 1, 2, 3], &unsettled);
        assert_eq!(forward, window(&locations, &reads, &[3, 2, 1, 0], &unsettled));
        assert_eq!(forward, window(&locations, &reads, &[2, 0, 3, 1], &unsettled));
        // the lowest position of each group becomes the vertex; the positions within 6 bp of it join it
        assert_eq!(forward.1, vec![100, 300]);
        assert_eq!(forward.0, vec![vec![2, 3], vec![2, 3], vec![2, 3], vec![2, 1]]);
    }

    #[test]
    fn matching_reads_one_after_another_depends_on_their_order() {
        // with every novel vertex visible at once, the first read's position wins
        let locations = LocationDict::new();
        let config = RunConfig::default();
        let reads = [[1, 103, 300, 400], [1, 100, 302, 400]];
        let made = |order: [usize; 2]| {
            let counters = Counters::new(0, 0, 0, 0, 0, 0);
            let (mut novel, overlay) = (NovelBuffer::new(), ConcurrentLocationDict::new());
            for i in order {
                match_splice_vertices(&locations, "chr1", &reads[i], Strand::Plus, &counters, &config, &mut novel, &overlay,
                                      None, None, 0, ALL_VISIBLE);
            }
            novel.locations.iter().map(|l| l.position).collect::<Vec<_>>()
        };
        assert_eq!(made([0, 1]), vec![103, 300]);
        assert_eq!(made([1, 0]), vec![100, 302]);
    }

    #[test]
    fn window_vertices_thin_by_distance_and_keep_settled_sites() {
        let config = RunConfig::default();
        let positions = |proposed: &mut Vec<(&str, i64, bool)>, thin| {
            let counters = Counters::new(0, 0, 0, 0, 0, 0);
            let (mut made, overlay) = (NovelBuffer::new(), ConcurrentLocationDict::new());
            create_window_vertices(proposed, thin, &counters, &config, &mut made, &overlay);
            made.locations.iter().map(|l| (l.chromosome.clone(), l.position, l.location_id)).collect::<Vec<_>>()
        };
        let at = |v: &[(&str, i64)]| v.iter().enumerate().map(|(i, &(c, p))| (c.to_string(), p, i as i64 + 1)).collect::<Vec<_>>();
        // the distance is taken from the last vertex made, not from the last position seen; IDs follow the coordinate
        let mut p = vec![("chr2", 50, false), ("chr1", 110, false), ("chr1", 100, false), ("chr1", 105, false), ("chr1", 105, false)];
        assert_eq!(positions(&mut p.clone(), Some(6)), at(&[("chr1", 100), ("chr1", 110), ("chr2", 50)]));
        // with a window consensus every target is made
        assert_eq!(positions(&mut p, None), at(&[("chr1", 100), ("chr1", 105), ("chr1", 110), ("chr2", 50)]));
        // a site the splice-site consensus settled is made where it lies, whichever read proposed it first
        let mut q = vec![("chr1", 103, false), ("chr1", 100, false), ("chr1", 103, true)];
        assert_eq!(positions(&mut q, Some(6)), at(&[("chr1", 100), ("chr1", 103)]));
    }

    #[test]
    fn vertices_above_the_visible_bound_are_not_seen() {
        let locations = LocationDict::new();
        let counters = Counters::new(0, 0, 0, 0, 0, 0);
        let config = RunConfig::default();
        let (mut novel, overlay) = (NovelBuffer::new(), ConcurrentLocationDict::new());
        let (first, _) = create_or_get_vertex("chr1", 100, &counters, &config, &mut novel, &overlay);
        let (second, _) = create_or_get_vertex("chr1", 200, &counters, &config, &mut novel, &overlay);
        assert_eq!(search_vertex_at_pos(&locations, "chr1", 200, Some(&overlay), first), None);
        assert_eq!(search_vertex_at_pos(&locations, "chr1", 200, Some(&overlay), second), Some(second));
        assert_eq!(search_nearby_splice_vertex(&locations, "chr1", 203, 6, &overlay, first), None);
        assert_eq!(search_nearby_splice_vertex(&locations, "chr1", 203, 6, &overlay, ALL_VISIBLE), Some(second));
        assert_eq!(search_nearby_splice_vertex(&locations, "chr1", 103, 6, &overlay, first), Some(first));
    }
}

/// Match all transcript vertices, with permissive matching at the ends.
pub fn match_all_transcript_vertices(
    location_dict: &LocationDict,
    chromosome: &str,
    positions: &[i64],
    strand: Strand,
    counters: &Counters,
    config: &RunConfig,
    novel_buffer: &mut NovelBuffer,
    novel_location_dict: &ConcurrentLocationDict,
    sr_splice_positions: Option<&rustc_hash::FxHashSet<(String, i64)>>,
) -> (Vec<VertexId>, Vec<u8>, Option<i64>, Option<i64>) {
    let mut vertex_matches = Vec::with_capacity(positions.len());
    let mut novelty = Vec::with_capacity(positions.len());
    let mut diff_5p = None;
    let mut diff_3p = None;
    let last = positions.len() - 1;

    for i in 0..positions.len() {
        let position = positions[i];

        if i == 0 {
            // Start position: permissive match
            let sj_pos = positions[i + 1];
            let (vid, delta) = permissive_vertex_search(
                location_dict,
                chromosome,
                position,
                strand,
                sj_pos,
                "start",
                config,
            );
            diff_5p = delta;
            match vid {
                Some(v) => {
                    vertex_matches.push(v);
                    novelty.push(0);
                }
                None => {
                    let (vid, _) = create_or_get_vertex(
                        chromosome, position, counters, config, novel_buffer, novel_location_dict,
                    );
                    vertex_matches.push(vid);
                    novelty.push(1);
                }
            }
        } else if i == last {
            // End position: permissive match
            let sj_pos = positions[i - 1];
            let (vid, delta) = permissive_vertex_search(
                location_dict,
                chromosome,
                position,
                strand,
                sj_pos,
                "end",
                config,
            );
            diff_3p = delta;
            match vid {
                Some(v) => {
                    vertex_matches.push(v);
                    novelty.push(0);
                }
                None => {
                    let (vid, _) = create_or_get_vertex(
                        chromosome, position, counters, config, novel_buffer, novel_location_dict,
                    );
                    vertex_matches.push(vid);
                    novelty.push(1);
                }
            }
        } else {
            // Internal position: exact match first, then nearby snap
            match search_vertex_at_pos(location_dict, chromosome, position, Some(novel_location_dict), ALL_VISIBLE) {
                Some(vid) => {
                    vertex_matches.push(vid);
                    novelty.push(0);
                }
                None => {
                    // Nearby snap
                    if config.max_splice_dist > 0 {
                        if let Some(vid) = search_nearby_splice_vertex(
                            location_dict, chromosome, position,
                            config.max_splice_dist, novel_location_dict, ALL_VISIBLE,
                        ) {
                            vertex_matches.push(vid);
                            novelty.push(0);
                            continue;
                        }
                    }
                    // Extended reference-only snap
                    if config.max_extended_snap > config.max_splice_dist {
                        if let Some(vid) = search_extended_ref_splice_vertex(
                            location_dict, chromosome, position,
                            config.max_splice_dist, config.max_extended_snap,
                            sr_splice_positions,
                        ) {
                            vertex_matches.push(vid);
                            novelty.push(0);
                            continue;
                        }
                    }
                    let (vid, _) = create_or_get_vertex(
                        chromosome, position, counters, config, novel_buffer, novel_location_dict,
                    );
                    vertex_matches.push(vid);
                    novelty.push(1);
                }
            }
        }
    }

    (vertex_matches, novelty, diff_5p, diff_3p)
}

/// Match monoexon vertices (start and end, permissive).
pub fn match_monoexon_vertices(
    location_dict: &LocationDict,
    chromosome: &str,
    positions: &[i64],
    strand: Strand,
    counters: &Counters,
    config: &RunConfig,
    novel_buffer: &mut NovelBuffer,
    novel_location_dict: &ConcurrentLocationDict,
    sr_splice_positions: Option<&rustc_hash::FxHashSet<(String, i64)>>,
) -> (Vec<VertexId>, Vec<u8>, Option<i64>, Option<i64>) {
    assert!(positions.len() == 2, "Monoexon should have exactly 2 positions");
    match_all_transcript_vertices(location_dict, chromosome, positions, strand, counters, config, novel_buffer, novel_location_dict, sr_splice_positions)
}

/// Process the 5' end of a transcript.
pub fn process_5p(
    location_dict: &LocationDict,
    edge_dict: &EdgeDict,
    chromosome: &str,
    positions: &[i64],
    strand: Strand,
    vertex_ids: &[VertexId],
    gene_id: GeneId,
    gene_starts: &GeneEndpointDict,
    counters: &Counters,
    config: &RunConfig,
    novel_buffer: &mut NovelBuffer,
    novel_location_dict: &ConcurrentLocationDict,
    novel_edge_dict: &ConcurrentEdgeDict,
) -> (VertexId, EdgeId, u8, bool, Option<i64>) {
    let (start_vertex, diff_5p, known_start) = permissive_match_with_gene_priority(
        location_dict,
        chromosome,
        positions[0],
        strand,
        positions[1],
        "start",
        gene_id,
        gene_starts,
        config,
    );

    let start_vertex = match start_vertex {
        Some(v) => v,
        None => {
            let (vid, _) = create_or_get_vertex(
                chromosome, positions[0], counters, config, novel_buffer, novel_location_dict,
            );
            vid
        }
    };

    // Create the start exon edge
    let (start_exon, start_novelty) = super::edge_matching::match_or_create_edge(
        start_vertex,
        vertex_ids[0],
        EdgeType::Exon,
        strand,
        edge_dict,
        counters,
        novel_buffer,
        novel_edge_dict,
    );

    (start_vertex, start_exon, start_novelty, known_start, diff_5p)
}

/// Process the 3' end of a transcript.
pub fn process_3p(
    location_dict: &LocationDict,
    edge_dict: &EdgeDict,
    chromosome: &str,
    positions: &[i64],
    strand: Strand,
    vertex_ids: &[VertexId],
    gene_id: GeneId,
    gene_ends: &GeneEndpointDict,
    counters: &Counters,
    config: &RunConfig,
    novel_buffer: &mut NovelBuffer,
    novel_location_dict: &ConcurrentLocationDict,
    novel_edge_dict: &ConcurrentEdgeDict,
) -> (VertexId, EdgeId, u8, bool, Option<i64>) {
    let last = positions.len() - 1;
    let (end_vertex, diff_3p, known_end) = permissive_match_with_gene_priority(
        location_dict,
        chromosome,
        positions[last],
        strand,
        positions[last - 1],
        "end",
        gene_id,
        gene_ends,
        config,
    );

    let end_vertex = match end_vertex {
        Some(v) => v,
        None => {
            let (vid, _) = create_or_get_vertex(
                chromosome, positions[last], counters, config, novel_buffer, novel_location_dict,
            );
            vid
        }
    };

    // Create the end exon edge
    let last_vid = vertex_ids.len() - 1;
    let (end_exon, end_novelty) = super::edge_matching::match_or_create_edge(
        vertex_ids[last_vid],
        end_vertex,
        EdgeType::Exon,
        strand,
        edge_dict,
        counters,
        novel_buffer,
        novel_edge_dict,
    );

    (end_vertex, end_exon, end_novelty, known_end, diff_3p)
}
