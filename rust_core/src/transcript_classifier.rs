// Transcript classification: FSM, ISM, NIC, NNC, Antisense, Intergenic,
// Genomic, Fusion. FSM/ISM candidates are the annotated hyperedges that contain the read's chain.

use dashmap::DashMap;

use crate::data_structures::*;
use crate::edge_matching;
use crate::gene_assignment::{self, GeneAssignment};
use crate::transcript_utils::{self, compute_delta, construct_names};
use crate::types::*;
use crate::vertex_matching;

/// Shared concurrent registry of novel transcripts (sorted edge path →
/// transcript_id), so reads with identical exon structures share one ID.
pub type NovelTranscriptRegistry = DashMap<TranscriptPathKey, TranscriptId>;

/// Classify and annotate a multi-exon transcript.
pub fn identify_transcript(
    chrom: &str,
    positions: &[i64],
    strand: Strand,
    refs: &SharedRefs,
    novel_buffer: &mut NovelBuffer,
    novel_tx_registry: &NovelTranscriptRegistry,
    consensus: Option<&rustc_hash::FxHashMap<i64, (i64, u32)>>,
    ambiguous_sites: u64,
    settled_sites: u64,
    small_rna_read: bool,
) -> ReadAnnotation {
    let n_exons = positions.len() as i32 / 2;

    // ── Monoexon dispatch ────────────────────────────────────────────
    if n_exons <= 1 {
        return identify_monoexon_transcript(chrom, positions, strand, refs, novel_buffer, novel_tx_registry, small_rna_read);
    }

    let mut gene_novelty: Vec<NoveltyAnnotation> = Vec::new();
    let mut transcript_novelty: Vec<NoveltyAnnotation> = Vec::new();
    let mut hint_flags: u8 = 0;
    let mut ism_parent: Option<TranscriptId> = None;

    // Step 1: Match splice junction vertices
    let (vertex_ids, v_novelty, max_splice_shift) = vertex_matching::match_splice_vertices(
        &refs.location_dict,
        chrom,
        positions,
        strand,
        &refs.counters,
        &refs.config,
        novel_buffer,
        &refs.novel_location_dict,
        refs.sr_splice_positions.as_ref(),
        consensus,
        settled_sites,
        refs.visible_vertex_max.load(std::sync::atomic::Ordering::Relaxed),
    );

    // Step 2: Match splice edges
    let (edge_ids, e_novelty) = edge_matching::match_all_splice_edges(
        &vertex_ids,
        strand,
        &refs.edge_dict,
        &refs.counters,
        novel_buffer,
        &refs.novel_edge_dict,
    );

    // Step 3: Check novelty patterns
    let all_sjs_known = transcript_utils::check_all_sjs_known(&e_novelty);
    let _all_exons_known = transcript_utils::check_all_exons_known(&e_novelty);
    let all_exons_novel = transcript_utils::check_all_exons_novel(&e_novelty);
    let splice_vertices_known = v_novelty.iter().all(|&v| v == 0);

    let mut gene_id: Option<GeneId> = None;
    let mut transcript_id: Option<TranscriptId> = None;
    let mut start_end_info: Option<StartEndInfo> = None;
    let mut fusion = false;

    // ─── FSM / ISM search ────────────────────────────────────────────────
    if all_sjs_known {
        // Candidates: the annotated hyperedges holding every edge of the read's chain
        let all_matches = refs.hypergraph.containing(&edge_ids);

        if !all_matches.is_empty() {
            // Try FSM first
            let fsm_result = process_fsm(
                chrom, positions, strand, &edge_ids, &vertex_ids,
                &all_matches, &refs, novel_buffer, novel_tx_registry,
            );

            if let Some((gid, tid, nov, sei, fsm_flags)) = fsm_result {
                gene_id = Some(gid);
                transcript_id = Some(tid);
                transcript_novelty = nov;
                start_end_info = Some(sei);
                hint_flags |= HINT_IS_KNOWN | fsm_flags;
            } else {
                // Try ISM
                let ism_result = process_ism(
                    chrom, positions, strand, &edge_ids, &vertex_ids,
                    &all_matches, refs, novel_buffer, novel_tx_registry,
                );
                if let Some((gid, tid, nov, sei, ism_flags, parent)) = ism_result {
                    gene_id = Some(gid);
                    transcript_id = Some(tid);
                    transcript_novelty = nov;
                    start_end_info = Some(sei);
                    hint_flags |= ism_flags;
                    ism_parent = parent;
                }
            }
        }

    }

    // ─── Chain-first assignment ──────────────────────────────────────────
    // A read whose junctions all lie within tolerance of one annotated
    // transcript, with junctions on annotated sites held in place, is assigned
    // to it before the NIC path is tried.
    if gene_id.is_none() && refs.config.chain_first && refs.config.approx_chain_tolerance > 0 {
        if let Some((gid, tid, nov, sei)) =
            process_approximate_fsm(chrom, positions, strand, refs, novel_buffer, "chain_first", ambiguous_sites, settled_sites)
        {
            gene_id = Some(gid);
            transcript_id = Some(tid);
            transcript_novelty = nov;
            start_end_info = Some(sei);
            hint_flags |= HINT_IS_KNOWN;
        }
    }

    // ─── NIC (known vertices, new connections) ───────────────────────────
    if gene_id.is_none() && (all_sjs_known || splice_vertices_known) {
        let nic_result = process_nic(
            chrom, positions, strand, &edge_ids, &vertex_ids,
            refs, novel_buffer, novel_tx_registry,
        );
        if let Some((gid, tid, nov, sei, fus)) = nic_result {
            gene_id = Some(gid);
            transcript_id = Some(tid);
            transcript_novelty = nov;
            start_end_info = Some(sei);
            fusion = fus;
        } else {
            // NIC returned None — check if assign_gene detected fusion
            let (_, fus) = gene_assignment::assign_gene(
                &vertex_ids, strand, &refs.vertex_2_gene,
                chrom, positions[0], *positions.last().unwrap(),
                &refs.gene_tree, &refs.transcript_tree,
            );
            fusion = fus;
        }
    }

    // ─── Approximate FSM rescue ──────────────────────────────────────────
    // Without chain-first, the chain test runs after NIC has failed and before
    // an NNC model is created.
    if gene_id.is_none() && !fusion && !refs.config.chain_first && refs.config.approx_chain_tolerance > 0 {
        if let Some((gid, tid, nov, sei)) =
            process_approximate_fsm(chrom, positions, strand, refs, novel_buffer, "post_nic", ambiguous_sites, settled_sites)
        {
            gene_id = Some(gid);
            transcript_id = Some(tid);
            transcript_novelty = nov;
            start_end_info = Some(sei);
            hint_flags |= HINT_IS_KNOWN;
        }
    }

    // Antisense (known splice vertices on opposite strand)
    if splice_vertices_known && gene_id.is_none() && !fusion {
        let anti_result = process_spliced_antisense(
            chrom, positions, strand, &edge_ids, &vertex_ids,
            refs, novel_buffer, novel_tx_registry,
        );
        if let Some((gid, tid, g_nov, t_nov, sei)) = anti_result {
            gene_id = Some(gid);
            transcript_id = Some(tid);
            gene_novelty = g_nov;
            transcript_novelty = t_nov;
            start_end_info = Some(sei);
        }
    }
    // NNC (novel splice donors/acceptors, at least one known exon)
    else if !splice_vertices_known && !fusion && !all_exons_novel && gene_id.is_none() {
        let nnc_result = process_nnc(
            chrom, positions, strand, &edge_ids, &vertex_ids,
            refs, novel_buffer, novel_tx_registry,
        );
        if let Some((gid, tid, nov, sei, fus)) = nnc_result {
            gene_id = Some(gid);
            transcript_id = Some(tid);
            transcript_novelty = nov;
            start_end_info = Some(sei);
            fusion = fus;
        } else {
            // NNC returned None — check if assign_gene detected fusion
            let (_, fus) = gene_assignment::assign_gene(
                &vertex_ids, strand, &refs.vertex_2_gene,
                chrom, positions[0], *positions.last().unwrap(),
                &refs.gene_tree, &refs.transcript_tree,
            );
            fusion = fus;
        }
    }

    // Remaining cases: intergenic, genomic, fusion
    if gene_id.is_none() {
        let remaining_result = process_remaining(
            chrom, positions, strand, &edge_ids, &vertex_ids,
            refs, novel_buffer, novel_tx_registry, fusion,
        );
        let (gid, tid, g_nov, t_nov, sei) = remaining_result;
        gene_id = Some(gid);
        transcript_id = Some(tid);
        gene_novelty = g_nov;
        transcript_novelty = t_nov;
        start_end_info = Some(sei);
    }

    // ─── Finalize ────────────────────────────────────────────────────────
    let sei = start_end_info.unwrap_or_default();
    let gene_id = gene_id.unwrap();
    let transcript_id = transcript_id.unwrap();

    // Update vertex_2_gene
    gene_assignment::update_vertex_2_gene(gene_id, &sei.vertex_ids, strand, novel_buffer);

    // Add names for novel entities
    let (gene_name, transcript_name) =
        construct_names(gene_id, transcript_id, &refs.config.id_prefix, refs.config.n_places);

    if !gene_novelty.is_empty() {
        gene_novelty.push(make_annotation(gene_id, &refs.config, "gene_status", "NOVEL"));
        gene_novelty.push(make_annotation(gene_id, &refs.config, "gene_name", &gene_name));
        gene_novelty.push(make_annotation(gene_id, &refs.config, "gene_id", &gene_name));
    }
    if !transcript_novelty.is_empty() {
        transcript_novelty.push(make_annotation(transcript_id, &refs.config, "transcript_status", "NOVEL"));
        transcript_novelty.push(make_annotation(transcript_id, &refs.config, "transcript_name", &transcript_name));
        transcript_novelty.push(make_annotation(transcript_id, &refs.config, "transcript_id", &transcript_name));
    }

    // Exon novelty
    let mut exon_novelty = Vec::new();
    let full_e_novelty = build_full_edge_novelty(&e_novelty, &sei);
    let exon_edges = collect_exon_edges(&sei.edge_ids, &full_e_novelty);
    for (exon_id, is_novel) in exon_edges {
        if is_novel {
            exon_novelty.push(make_annotation(exon_id, &refs.config, "exon_status", "NOVEL"));
        }
    }

    ReadAnnotation {
        gene_id,
        transcript_id,
        gene_novelty,
        transcript_novelty,
        exon_novelty,
        start_vertex: sei.start_vertex,
        end_vertex: sei.end_vertex,
        start_exon: sei.start_exon,
        end_exon: sei.end_exon,
        start_delta: sei.diff_5p,
        end_delta: sei.diff_3p,
        // These are filled in by the caller
        read_id: String::new(),
        read_idx: READ_IDX_NONE,
        dataset: empty_arc_str(),
        location: String::new(),
        strand,
        read_length: 0,
        start_pos: positions[0],
        end_pos: *positions.last().unwrap(),
        n_exons,
        max_splice_shift,
        introns: positions.get(1..positions.len().saturating_sub(1)).map_or_else(Vec::new, |s| s.to_vec()),
        cap_signature: false,
        unaligned_3p: false,
        fraction_as: None,
        polya_tail_len: None,
        polya_tail_softclip: None,
        custom_label: None,
        allelic_label: None,
        start_support: None,
        end_support: None,
        internal_priming: None,
        end_site_info: None,
        hint_flags,
        ism_parent,
    }
}

/// Process FSM (Full Splice Match).
///
/// A read matching all internal splice junctions of an annotated transcript
/// is assigned to it. Terminals beyond the 5'/3' cutoffs are recorded as
/// novel ends on the annotated ID, or with `fsm_terminal_new_id` create a
/// novel FSM terminal variant keyed on the internal edges.
fn process_fsm(
    chrom: &str,
    positions: &[i64],
    strand: Strand,
    edge_ids: &[EdgeId],
    vertex_ids: &[VertexId],
    all_matches: &[&TranscriptModel],
    refs: &SharedRefs,
    novel_buffer: &mut NovelBuffer,
    novel_tx_registry: &NovelTranscriptRegistry,
) -> Option<(GeneId, TranscriptId, Vec<NoveltyAnnotation>, StartEndInfo, u8)> {
    let n_exons = positions.len() as i32 / 2;

    // FSM matches must have same number of exons
    let fsm_matches: Vec<&&TranscriptModel> = all_matches
        .iter()
        .filter(|m| m.n_exons == n_exons)
        .collect();

    if fsm_matches.is_empty() {
        return None;
    }

    // Among several FSM candidates, choose the one whose annotated last exon end
    // is closest to the read's.
    let read_3p = *positions.last().unwrap();
    let transcript_match = if fsm_matches.len() > 1 {
        fsm_matches.iter().min_by_key(|m| {
            let delta = compute_delta(m.end_pos, read_3p, strand);
            delta.abs()
        }).unwrap()
    } else {
        fsm_matches[0]
    };
    let gene_id = transcript_match.gene_id;
    let known_transcript_id = transcript_match.transcript_id;

    // ── Nested-isoform ambiguity ────────────────────────────────────────
    // A read whose end lies inside a longer containing isoform's exon, far from
    // this short isoform's annotated end, is flagged `HINT_ISM_AMBIGUOUS`.
    let mut fsm_flags: u8 = 0;
    if refs.config.nested_fsm_ambiguous {
        if let Some(index) = refs.parent_chains.as_ref() {
            if crate::consistent_assignment::nested_in_longer(
                index, known_transcript_id, n_exons, all_matches,
                positions[0], *positions.last().unwrap(), refs.config.nested_fsm_end_window,
            ) {
                fsm_flags |= HINT_ISM_AMBIGUOUS;
            }
        }
    }

    let curr_5p_diff = compute_delta(transcript_match.start_pos, positions[0], strand);
    let curr_3p_diff = compute_delta(transcript_match.end_pos, *positions.last().unwrap(), strand);

    let mut sei = StartEndInfo::default();

    // 5' end
    let novel_5p = curr_5p_diff.abs() > refs.config.cutoff_5p;
    if !novel_5p {
        sei.start_vertex = transcript_match.start_vertex;
        sei.start_exon = transcript_match.start_exon;
        sei.diff_5p = Some(curr_5p_diff);
        sei.start_novelty = 0;
    } else {
        let (sv, se, sn, _, d5) = vertex_matching::process_5p(
            &refs.location_dict, &refs.edge_dict, chrom, positions, strand,
            vertex_ids, gene_id, &refs.gene_starts, &refs.counters, &refs.config, novel_buffer,
            &refs.novel_location_dict, &refs.novel_edge_dict,
        );
        sei.start_vertex = sv;
        sei.start_exon = se;
        sei.start_novelty = sn;
        sei.diff_5p = d5;
    }

    // 3' end
    let novel_3p = curr_3p_diff.abs() > refs.config.cutoff_3p;
    if !novel_3p {
        sei.end_vertex = transcript_match.end_vertex;
        sei.end_exon = transcript_match.end_exon;
        sei.diff_3p = Some(curr_3p_diff);
        sei.end_novelty = 0;
    } else {
        let (ev, ee, en, _, d3) = vertex_matching::process_3p(
            &refs.location_dict, &refs.edge_dict, chrom, positions, strand,
            vertex_ids, gene_id, &refs.gene_ends, &refs.counters, &refs.config, novel_buffer,
            &refs.novel_location_dict, &refs.novel_edge_dict,
        );
        sei.end_vertex = ev;
        sei.end_exon = ee;
        sei.end_novelty = en;
        sei.diff_3p = d3;
    }

    // Build full edge and vertex lists
    let mut full_edges = vec![sei.start_exon];
    full_edges.extend_from_slice(edge_ids);
    full_edges.push(sei.end_exon);

    let mut full_vertices = vec![sei.start_vertex];
    full_vertices.extend_from_slice(vertex_ids);
    full_vertices.push(sei.end_vertex);

    sei.edge_ids = full_edges.clone();
    sei.vertex_ids = full_vertices.clone();

    // ── Decide: canonical FSM vs FSM terminal variant ────────────────
    if novel_5p || novel_3p {
        // The annotated ID is kept unless fsm_terminal_new_id is set.
        let tid = if refs.config.fsm_terminal_new_id {
            create_novel_transcript(
                strand, chrom, positions[0], *positions.last().unwrap(),
                gene_id, &full_edges, &full_vertices,
                &refs.counters, novel_buffer, novel_tx_registry,
            )
            .transcript_id
        } else {
            known_transcript_id
        };

        let mut novelty = Vec::new();
        novelty.push(make_annotation(tid, &refs.config, "FSM_transcript", "TRUE"));
        novelty.push(make_annotation(tid, &refs.config, "FSM_to_IDs",
            &known_transcript_id.to_string()));
        if novel_5p {
            novelty.push(make_annotation(tid, &refs.config, "FSM_novel_5p", "TRUE"));
        }
        if novel_3p {
            novelty.push(make_annotation(tid, &refs.config, "FSM_novel_3p", "TRUE"));
        }

        Some((gene_id, tid, novelty, sei, fsm_flags))
    } else {
        // Both terminals within cutoff — canonical FSM assignment
        Some((gene_id, known_transcript_id, Vec::new(), sei, fsm_flags))
    }
}

/// Process ISM (Incomplete Splice Match)
/// Returns (gene_id, transcript_id, novelty, start_end_info, ism_hint_flags, ism_parent)
fn process_ism(
    chrom: &str,
    positions: &[i64],
    strand: Strand,
    edge_ids: &[EdgeId],
    vertex_ids: &[VertexId],
    all_matches: &[&TranscriptModel],
    refs: &SharedRefs,
    novel_buffer: &mut NovelBuffer,
    novel_tx_registry: &NovelTranscriptRegistry,
) -> Option<(GeneId, TranscriptId, Vec<NoveltyAnnotation>, StartEndInfo, u8, Option<TranscriptId>)> {
    let n_exons = positions.len() as i32 / 2;

    // Choose gene: pick from matching genes, tiebreak if needed
    let gene_ids: Vec<GeneId> = all_matches.iter().map(|m| m.gene_id).collect();
    let unique_genes: Vec<GeneId> = {
        let mut g = gene_ids.clone();
        g.sort_unstable();
        g.dedup();
        g
    };

    let gene_id = if unique_genes.len() > 1 {
        match gene_assignment::search_for_overlap_with_gene(
            chrom, positions[0], *positions.last().unwrap(), strand,
            &refs.gene_tree, &refs.transcript_tree, Some(&unique_genes),
        ) {
            Some((gid, _)) => gid,
            None => return None,
        }
    } else {
        unique_genes[0]
    };

    // Process ends
    let mut sei = StartEndInfo::default();
    if n_exons > 1 {
        let (sv, se, sn, _, d5) = vertex_matching::process_5p(
            &refs.location_dict, &refs.edge_dict, chrom, positions, strand,
            vertex_ids, gene_id, &refs.gene_starts, &refs.counters, &refs.config, novel_buffer,
            &refs.novel_location_dict, &refs.novel_edge_dict,
        );
        let (ev, ee, en, _, d3) = vertex_matching::process_3p(
            &refs.location_dict, &refs.edge_dict, chrom, positions, strand,
            vertex_ids, gene_id, &refs.gene_ends, &refs.counters, &refs.config, novel_buffer,
            &refs.novel_location_dict, &refs.novel_edge_dict,
        );

        sei.start_vertex = sv;
        sei.start_exon = se;
        sei.start_novelty = sn;
        sei.diff_5p = d5;
        sei.end_vertex = ev;
        sei.end_exon = ee;
        sei.end_novelty = en;
        sei.diff_3p = d3;
    }

    let mut full_edges = vec![sei.start_exon];
    full_edges.extend_from_slice(edge_ids);
    full_edges.push(sei.end_exon);

    let mut full_vertices = vec![sei.start_vertex];
    full_vertices.extend_from_slice(vertex_ids);
    full_vertices.push(sei.end_vertex);

    sei.edge_ids = full_edges.clone();
    sei.vertex_ids = full_vertices.clone();

    // Create novel transcript
    let transcript_model = create_novel_transcript(
        strand, chrom, positions[0], *positions.last().unwrap(),
        gene_id, &full_edges, &full_vertices, &refs.counters, novel_buffer, novel_tx_registry,
    );
    let transcript_id = transcript_model.transcript_id;

    // Build novelty annotations
    let mut novelty = Vec::new();
    novelty.push(make_annotation(transcript_id, &refs.config, "ISM_transcript", "TRUE"));

    // Collect ISM IDs
    let ism_ids: Vec<String> = all_matches.iter().map(|m| m.transcript_id.to_string()).collect();
    novelty.push(make_annotation(transcript_id, &refs.config, "ISM_to_IDs", &ism_ids.join(",")));

    // Check prefix/suffix
    let edge_str = if n_exons > 1 {
        edge_ids.iter().map(|e| e.to_string()).collect::<Vec<_>>().join(",")
    } else {
        full_edges[0].to_string()
    };

    let mut prefix_ids = Vec::new();
    let mut suffix_ids = Vec::new();
    for m in all_matches {
        if let Some(ref jn) = m.jn_path {
            if jn.starts_with(&edge_str) {
                prefix_ids.push(m.transcript_id.to_string());
            }
            if jn.ends_with(&edge_str) {
                suffix_ids.push(m.transcript_id.to_string());
            }
        }
    }

    let mut ism_flags: u8 = HINT_IS_ISM;
    if !prefix_ids.is_empty() {
        novelty.push(make_annotation(transcript_id, &refs.config, "ISM-prefix_transcript", "TRUE"));
        novelty.push(make_annotation(transcript_id, &refs.config, "ISM-prefix_to_IDs", &prefix_ids.join(",")));
        ism_flags |= HINT_IS_ISM_PREFIX;
    }
    if !suffix_ids.is_empty() {
        novelty.push(make_annotation(transcript_id, &refs.config, "ISM-suffix_transcript", "TRUE"));
        novelty.push(make_annotation(transcript_id, &refs.config, "ISM-suffix_to_IDs", &suffix_ids.join(",")));
        ism_flags |= HINT_IS_ISM_SUFFIX;
    }

    // ── Consistency with the annotated parents ──────────────────────────
    // Records whether the truncated read is consistent with exactly one
    // annotated transcript; the engine later decides whether it moves onto it.
    let mut ism_parent = None;
    if let Some(index) = refs.parent_chains.as_ref() {
        use crate::consistent_assignment::{unique_parent, Parents};
        match unique_parent(index, all_matches, edge_ids, positions[0], *positions.last().unwrap()) {
            Parents::Unique(c) => {
                ism_flags |= HINT_ISM_UNIQUE;
                if c.intact_5p {
                    ism_flags |= HINT_ISM_5P_INTACT;
                }
                if c.intact_3p {
                    ism_flags |= HINT_ISM_3P_INTACT;
                }
                ism_parent = Some(c.parent);
                novelty.push(make_annotation(transcript_id, &refs.config, "ISM_unique_parent", &c.parent.to_string()));
            }
            Parents::Ambiguous(n) => {
                ism_flags |= HINT_ISM_AMBIGUOUS;
                novelty.push(make_annotation(transcript_id, &refs.config, "ISM_consistent_parents", &n.to_string()));
            }
            Parents::NoParent => {}
        }
    }

    Some((gene_id, transcript_id, novelty, sei, ism_flags, ism_parent))
}

/// Process NIC (Novel In Catalog)
fn process_nic(
    chrom: &str,
    positions: &[i64],
    strand: Strand,
    edge_ids: &[EdgeId],
    vertex_ids: &[VertexId],
    refs: &SharedRefs,
    novel_buffer: &mut NovelBuffer,
    novel_tx_registry: &NovelTranscriptRegistry,
) -> Option<(GeneId, TranscriptId, Vec<NoveltyAnnotation>, StartEndInfo, bool)> {
    let (assignment, fusion) = gene_assignment::assign_gene(
        vertex_ids, strand, &refs.vertex_2_gene,
        chrom, positions[0], *positions.last().unwrap(),
        &refs.gene_tree, &refs.transcript_tree,
    );

    let gene_id = match assignment {
        GeneAssignment::Found(id) => id,
        // With a fusion, return None so the caller falls through to
        // process_remaining with fusion=true
        _ => return None,
    };

    let (sv, se, sn, _, d5) = vertex_matching::process_5p(
        &refs.location_dict, &refs.edge_dict, chrom, positions, strand,
        vertex_ids, gene_id, &refs.gene_starts, &refs.counters, &refs.config, novel_buffer,
        &refs.novel_location_dict, &refs.novel_edge_dict,
    );
    let (ev, ee, en, _, d3) = vertex_matching::process_3p(
        &refs.location_dict, &refs.edge_dict, chrom, positions, strand,
        vertex_ids, gene_id, &refs.gene_ends, &refs.counters, &refs.config, novel_buffer,
        &refs.novel_location_dict, &refs.novel_edge_dict,
    );

    let mut full_edges = vec![se];
    full_edges.extend_from_slice(edge_ids);
    full_edges.push(ee);

    let mut full_vertices = vec![sv];
    full_vertices.extend_from_slice(vertex_ids);
    full_vertices.push(ev);

    let sei = StartEndInfo {
        start_vertex: sv, end_vertex: ev,
        start_exon: se, end_exon: ee,
        diff_5p: d5, diff_3p: d3,
        start_novelty: sn, end_novelty: en,
        vertex_ids: full_vertices.clone(),
        edge_ids: full_edges.clone(),
    };

    let transcript_model = create_novel_transcript(
        strand, chrom, positions[0], *positions.last().unwrap(),
        gene_id, &full_edges, &full_vertices, &refs.counters, novel_buffer, novel_tx_registry,
    );

    let novelty = vec![
        make_annotation(transcript_model.transcript_id, &refs.config, "NIC_transcript", "TRUE"),
    ];

    Some((gene_id, transcript_model.transcript_id, novelty, sei, fusion))
}

/// Process NNC (Novel Not in Catalog)
fn process_nnc(
    chrom: &str,
    positions: &[i64],
    strand: Strand,
    edge_ids: &[EdgeId],
    vertex_ids: &[VertexId],
    refs: &SharedRefs,
    novel_buffer: &mut NovelBuffer,
    novel_tx_registry: &NovelTranscriptRegistry,
) -> Option<(GeneId, TranscriptId, Vec<NoveltyAnnotation>, StartEndInfo, bool)> {
    let (assignment, fusion) = gene_assignment::assign_gene(
        vertex_ids, strand, &refs.vertex_2_gene,
        chrom, positions[0], *positions.last().unwrap(),
        &refs.gene_tree, &refs.transcript_tree,
    );

    let gene_id = match assignment {
        GeneAssignment::Found(id) => id,
        // With a fusion, return None so the caller falls through to
        // process_remaining with fusion=true
        _ => return None,
    };

    let (sv, se, sn, _, d5) = vertex_matching::process_5p(
        &refs.location_dict, &refs.edge_dict, chrom, positions, strand,
        vertex_ids, gene_id, &refs.gene_starts, &refs.counters, &refs.config, novel_buffer,
        &refs.novel_location_dict, &refs.novel_edge_dict,
    );
    let (ev, ee, en, _, d3) = vertex_matching::process_3p(
        &refs.location_dict, &refs.edge_dict, chrom, positions, strand,
        vertex_ids, gene_id, &refs.gene_ends, &refs.counters, &refs.config, novel_buffer,
        &refs.novel_location_dict, &refs.novel_edge_dict,
    );

    let mut full_edges = vec![se];
    full_edges.extend_from_slice(edge_ids);
    full_edges.push(ee);

    let mut full_vertices = vec![sv];
    full_vertices.extend_from_slice(vertex_ids);
    full_vertices.push(ev);

    let sei = StartEndInfo {
        start_vertex: sv, end_vertex: ev,
        start_exon: se, end_exon: ee,
        diff_5p: d5, diff_3p: d3,
        start_novelty: sn, end_novelty: en,
        vertex_ids: full_vertices.clone(),
        edge_ids: full_edges.clone(),
    };

    let transcript_model = create_novel_transcript(
        strand, chrom, positions[0], *positions.last().unwrap(),
        gene_id, &full_edges, &full_vertices, &refs.counters, novel_buffer, novel_tx_registry,
    );

    let novelty = vec![
        make_annotation(transcript_model.transcript_id, &refs.config, "NNC_transcript", "TRUE"),
    ];

    Some((gene_id, transcript_model.transcript_id, novelty, sei, fusion))
}

/// Process spliced antisense transcript
fn process_spliced_antisense(
    chrom: &str,
    positions: &[i64],
    strand: Strand,
    edge_ids: &[EdgeId],
    vertex_ids: &[VertexId],
    refs: &SharedRefs,
    novel_buffer: &mut NovelBuffer,
    novel_tx_registry: &NovelTranscriptRegistry,
) -> Option<(GeneId, TranscriptId, Vec<NoveltyAnnotation>, Vec<NoveltyAnnotation>, StartEndInfo)> {
    let anti_strand = strand.opposite();

    let anti_gene_result =
        gene_assignment::find_gene_match_on_vertex_basis(vertex_ids, anti_strand, &refs.vertex_2_gene);

    let anti_gene_id = match anti_gene_result {
        GeneAssignment::Found(id) => id,
        GeneAssignment::Ambiguous(ids) => {
            match gene_assignment::search_for_overlap_with_gene(
                chrom, positions[0], *positions.last().unwrap(), strand,
                &refs.gene_tree, &refs.transcript_tree, Some(&ids),
            ) {
                Some((id, _)) => id,
                None => return None,
            }
        }
        _ => return None,
    };

    // Process ends using the antisense gene's endpoints (swapped)
    let (sv, se, sn, _, d5) = vertex_matching::process_5p(
        &refs.location_dict, &refs.edge_dict, chrom, positions, strand,
        vertex_ids, anti_gene_id, &refs.gene_ends, &refs.counters, &refs.config, novel_buffer,
        &refs.novel_location_dict, &refs.novel_edge_dict,
    );
    let (ev, ee, en, _, d3) = vertex_matching::process_3p(
        &refs.location_dict, &refs.edge_dict, chrom, positions, strand,
        vertex_ids, anti_gene_id, &refs.gene_starts, &refs.counters, &refs.config, novel_buffer,
        &refs.novel_location_dict, &refs.novel_edge_dict,
    );

    let mut full_edges = vec![se];
    full_edges.extend_from_slice(edge_ids);
    full_edges.push(ee);

    let mut full_vertices = vec![sv];
    full_vertices.extend_from_slice(vertex_ids);
    full_vertices.push(ev);

    let sei = StartEndInfo {
        start_vertex: sv, end_vertex: ev,
        start_exon: se, end_exon: ee,
        diff_5p: d5, diff_3p: d3,
        start_novelty: sn, end_novelty: en,
        vertex_ids: full_vertices.clone(),
        edge_ids: full_edges.clone(),
    };

    // Create new gene
    let gene_id = get_or_create_novel_gene(
        NovelGeneKey::Antisense(anti_gene_id, strand),
        chrom, positions[0], *positions.last().unwrap(), strand,
        &refs.counters, &refs.novel_gene_dict,
    );

    // Create transcript
    let transcript_model = create_novel_transcript(
        strand, chrom, positions[0], *positions.last().unwrap(),
        gene_id, &full_edges, &full_vertices, &refs.counters, novel_buffer, novel_tx_registry,
    );

    let gene_novelty = vec![
        make_annotation(gene_id, &refs.config, "antisense_gene", "TRUE"),
        make_annotation(gene_id, &refs.config, "gene_antisense_to_IDs", &anti_gene_id.to_string()),
    ];
    let transcript_novelty = vec![
        make_annotation(transcript_model.transcript_id, &refs.config, "antisense_transcript", "TRUE"),
    ];

    Some((gene_id, transcript_model.transcript_id, gene_novelty, transcript_novelty, sei))
}

/// Identify and annotate a monoexon transcript.
///
///   1. Match start+end vertices permissively
///   2. Query the monoexon interval tree for annotated overlapping transcripts
///   3. If one overlaps on the same strand, use that gene/transcript
///   4. Otherwise create a novel gene/transcript or assign to an overlapping gene
fn identify_monoexon_transcript(
    chrom: &str,
    positions: &[i64],
    strand: Strand,
    refs: &SharedRefs,
    novel_buffer: &mut NovelBuffer,
    novel_tx_registry: &NovelTranscriptRegistry,
    small_rna_read: bool,
) -> ReadAnnotation {
    let mut gene_novelty: Vec<NoveltyAnnotation> = Vec::new();
    let mut transcript_novelty: Vec<NoveltyAnnotation> = Vec::new();
    let mut hint_flags: u8 = 0;
    let mut ism_parent: Option<TranscriptId> = None;
    let n_exons = 1i32;

    let read_start = positions[0];
    let read_end = positions[1];

    // Step 1: Match start/end vertices permissively
    let (vertex_ids, _v_novelty, diff_5p, diff_3p) = vertex_matching::match_monoexon_vertices(
        &refs.location_dict,
        chrom,
        positions,
        strand,
        &refs.counters,
        &refs.config,
        novel_buffer,
        &refs.novel_location_dict,
        refs.sr_splice_positions.as_ref(),
    );

    let start_vertex = vertex_ids[0];
    let end_vertex = vertex_ids[1];

    // Step 2: Create the single exon edge
    let (exon_edge, exon_novelty_flag) = edge_matching::match_or_create_edge(
        start_vertex,
        end_vertex,
        EdgeType::Exon,
        strand,
        &refs.edge_dict,
        &refs.counters,
        novel_buffer,
        &refs.novel_edge_dict,
    );

    let full_vertices = vec![start_vertex, end_vertex];
    let full_edges = vec![exon_edge];

    let sei = StartEndInfo {
        start_vertex,
        end_vertex,
        start_exon: exon_edge,
        end_exon: exon_edge,
        diff_5p,
        diff_3p,
        start_novelty: exon_novelty_flag,
        end_novelty: exon_novelty_flag,
        vertex_ids: full_vertices.clone(),
        edge_ids: full_edges.clone(),
    };

    // Step 3: Query monoexon tree for known overlapping transcripts on same strand
    let overlaps = refs.monoexon_tree.query(chrom, strand, read_start, read_end);
    // A read of an annotated small RNA goes to the small RNAs holding it (small_rna.rs).
    let overlaps = match &refs.small_rna {
        Some(index) if small_rna_read => index.holders(overlaps, read_start.min(read_end), read_start.max(read_end)),
        _ => overlaps,
    };

    let (gene_id, transcript_id) = if !overlaps.is_empty() {
        use crate::consistent_assignment::{monoexon_parent, MonoexonCandidate};
        let cands: Vec<MonoexonCandidate> = overlaps
            .iter()
            .map(|r| {
                let (lo, hi) = (r.start.min(r.end), r.start.max(r.end));
                let (five, three) = if strand == Strand::Minus { (hi, lo) } else { (lo, hi) };
                MonoexonCandidate { transcript_id: r.transcript_id, lo, hi, five, three }
            })
            .collect();
        let (i, c, unique) = monoexon_parent(&cands, read_start, read_end, strand).expect("overlaps is not empty");
        let best = &overlaps[i];
        if c.intact_5p {
            hint_flags |= HINT_IS_KNOWN;
        } else {
            // A read starting inside the transcript is a 5' fragment of it: it counts
            // for the transcript and places no start of its own.
            hint_flags |= HINT_IS_ISM | HINT_IS_ISM_SUFFIX | HINT_ISM_3P_INTACT;
            hint_flags |= if unique { HINT_ISM_UNIQUE } else { HINT_ISM_AMBIGUOUS };
            ism_parent = Some(best.transcript_id);
        }
        (best.gene_id, best.transcript_id)
    } else {
        // No monoexon overlap — check gene overlap
        let overlap_gene = gene_assignment::search_for_overlap_with_gene(
            chrom, read_start, read_end, strand,
            &refs.gene_tree, &refs.transcript_tree, None,
        );

        match overlap_gene {
            None => {
                // Intergenic
                let gid = get_or_create_novel_gene(
                    locus_key(chrom, read_start, read_end, strand, refs.config.novel_gene_bin),
                    chrom, read_start, read_end, strand,
                    &refs.counters, &refs.novel_gene_dict,
                );
                gene_novelty.push(make_annotation(gid, &refs.config, "intergenic_novel", "TRUE"));

                let tm = create_novel_transcript(
                    strand, chrom, read_start, read_end,
                    gid, &full_edges, &full_vertices, &refs.counters, novel_buffer, novel_tx_registry,
                );
                transcript_novelty.push(make_annotation(tm.transcript_id, &refs.config, "intergenic_transcript", "TRUE"));
                (gid, tm.transcript_id)
            }
            Some((gid, match_strand)) if match_strand != strand => {
                // Antisense
                let new_gid = get_or_create_novel_gene(
                    NovelGeneKey::Antisense(gid, strand),
                    chrom, read_start, read_end, strand,
                    &refs.counters, &refs.novel_gene_dict,
                );
                gene_novelty.push(make_annotation(new_gid, &refs.config, "antisense_gene", "TRUE"));

                let tm = create_novel_transcript(
                    strand, chrom, read_start, read_end,
                    new_gid, &full_edges, &full_vertices, &refs.counters, novel_buffer, novel_tx_registry,
                );
                transcript_novelty.push(make_annotation(tm.transcript_id, &refs.config, "antisense_transcript", "TRUE"));
                (new_gid, tm.transcript_id)
            }
            Some((gid, _)) => {
                // An unspliced read inside one exon of an annotated spliced transcript
                // is a fragment of that transcript: it counts for the parent and can
                // place no end of its own.
                use crate::consistent_assignment::{unspliced_parent, UnsplicedParent};
                let fragment = refs.parent_chains.as_ref().map(|index| {
                    let candidates: Vec<TranscriptId> = refs
                        .transcript_tree
                        .query_overlap(chrom, read_start, read_end)
                        .into_iter()
                        .filter(|t| t.strand == strand)
                        .map(|t| t.transcript_id)
                        .collect();
                    unspliced_parent(index, &candidates, read_start.min(read_end), read_start.max(read_end), read_end)
                });
                match fragment {
                    Some(UnsplicedParent::Unique(c)) | Some(UnsplicedParent::Ambiguous(c, _)) => {
                        hint_flags |= HINT_IS_ISM;
                        if !c.intact_5p { hint_flags |= HINT_IS_ISM_SUFFIX; }
                        if !c.intact_3p { hint_flags |= HINT_IS_ISM_PREFIX; }
                        if c.intact_5p { hint_flags |= HINT_ISM_5P_INTACT; }
                        if c.intact_3p { hint_flags |= HINT_ISM_3P_INTACT; }
                        hint_flags |= if matches!(fragment, Some(UnsplicedParent::Unique(_))) { HINT_ISM_UNIQUE } else { HINT_ISM_AMBIGUOUS };
                        ism_parent = Some(c.parent);
                        let parent_gene = refs.parent_chains.as_ref().and_then(|i| i.get(c.parent)).map_or(gid, |p| p.gene_id);
                        (parent_gene, c.parent)
                    }
                    _ => {
                        // Genomic — known gene, novel monoexon transcript
                        let tm = create_novel_transcript(
                            strand, chrom, read_start, read_end,
                            gid, &full_edges, &full_vertices, &refs.counters, novel_buffer, novel_tx_registry,
                        );
                        transcript_novelty.push(make_annotation(tm.transcript_id, &refs.config, "genomic_transcript", "TRUE"));
                        (gid, tm.transcript_id)
                    }
                }
            }
        }
    };

    // Exon novelty
    let exon_novelty = if exon_novelty_flag == 1 {
        vec![make_annotation(exon_edge, &refs.config, "exon_status", "NOVEL")]
    } else {
        Vec::new()
    };

    // Gene/transcript names for novel entities
    let (gene_name, transcript_name) =
        transcript_utils::construct_names(gene_id, transcript_id, &refs.config.id_prefix, refs.config.n_places);

    if !gene_novelty.is_empty() {
        gene_novelty.push(make_annotation(gene_id, &refs.config, "gene_status", "NOVEL"));
        gene_novelty.push(make_annotation(gene_id, &refs.config, "gene_name", &gene_name));
        gene_novelty.push(make_annotation(gene_id, &refs.config, "gene_id", &gene_name));
    }
    if !transcript_novelty.is_empty() {
        transcript_novelty.push(make_annotation(transcript_id, &refs.config, "transcript_status", "NOVEL"));
        transcript_novelty.push(make_annotation(transcript_id, &refs.config, "transcript_name", &transcript_name));
        transcript_novelty.push(make_annotation(transcript_id, &refs.config, "transcript_id", &transcript_name));
    }

    gene_assignment::update_vertex_2_gene(gene_id, &sei.vertex_ids, strand, novel_buffer);

    ReadAnnotation {
        gene_id,
        transcript_id,
        gene_novelty,
        transcript_novelty,
        exon_novelty,
        start_vertex: sei.start_vertex,
        end_vertex: sei.end_vertex,
        start_exon: sei.start_exon,
        end_exon: sei.end_exon,
        start_delta: sei.diff_5p,
        end_delta: sei.diff_3p,
        read_id: String::new(),
        read_idx: READ_IDX_NONE,
        dataset: empty_arc_str(),
        location: String::new(),
        strand,
        read_length: 0,
        start_pos: positions[0],
        end_pos: *positions.last().unwrap(),
        n_exons,
        max_splice_shift: 0,
        introns: positions.get(1..positions.len().saturating_sub(1)).map_or_else(Vec::new, |s| s.to_vec()),
        cap_signature: false,
        unaligned_3p: false,
        fraction_as: None,
        polya_tail_len: None,
        polya_tail_softclip: None,
        custom_label: None,
        allelic_label: None,
        start_support: None,
        end_support: None,
        internal_priming: None,
        end_site_info: None,
        hint_flags,
        ism_parent,
    }
}

/// Process remaining cases: intergenic, genomic, fusion
fn process_remaining(
    chrom: &str,
    positions: &[i64],
    strand: Strand,
    edge_ids: &[EdgeId],
    vertex_ids: &[VertexId],
    refs: &SharedRefs,
    novel_buffer: &mut NovelBuffer,
    novel_tx_registry: &NovelTranscriptRegistry,
    fusion: bool,
) -> (GeneId, TranscriptId, Vec<NoveltyAnnotation>, Vec<NoveltyAnnotation>, StartEndInfo) {
    let mut gene_novelty = Vec::new();
    let mut transcript_novelty = Vec::new();

    let overlap_gene = if !refs.config.create_novel_spliced_genes || !fusion {
        gene_assignment::search_for_overlap_with_gene(
            chrom, positions[0], *positions.last().unwrap(), strand,
            &refs.gene_tree, &refs.transcript_tree, None,
        )
    } else {
        None
    };

    // Process ends; without an overlapping gene the gene ID passed is 0
    let dummy_gene_id = overlap_gene.map(|(id, _)| id).unwrap_or(0);

    let (sv, se, sn, _, d5) = vertex_matching::process_5p(
        &refs.location_dict, &refs.edge_dict, chrom, positions, strand,
        vertex_ids, dummy_gene_id, &refs.gene_starts, &refs.counters, &refs.config, novel_buffer,
        &refs.novel_location_dict, &refs.novel_edge_dict,
    );
    let (ev, ee, en, _, d3) = vertex_matching::process_3p(
        &refs.location_dict, &refs.edge_dict, chrom, positions, strand,
        vertex_ids, dummy_gene_id, &refs.gene_ends, &refs.counters, &refs.config, novel_buffer,
        &refs.novel_location_dict, &refs.novel_edge_dict,
    );

    let mut full_edges = vec![se];
    full_edges.extend_from_slice(edge_ids);
    full_edges.push(ee);

    let mut full_vertices = vec![sv];
    full_vertices.extend_from_slice(vertex_ids);
    full_vertices.push(ev);

    let sei = StartEndInfo {
        start_vertex: sv, end_vertex: ev,
        start_exon: se, end_exon: ee,
        diff_5p: d5, diff_3p: d3,
        start_novelty: sn, end_novelty: en,
        vertex_ids: full_vertices.clone(),
        edge_ids: full_edges.clone(),
    };

    let (gene_id, transcript_id) = match overlap_gene {
        None => {
            let (t_nov, g_nov) = if fusion {
                ("fusion_transcript", "fusion_novel")
            } else {
                ("intergenic_transcript", "intergenic_novel")
            };

            let gid = get_or_create_novel_gene(
                locus_key(chrom, positions[0], *positions.last().unwrap(), strand, refs.config.novel_gene_bin),
                chrom, positions[0], *positions.last().unwrap(), strand,
                &refs.counters, &refs.novel_gene_dict,
            );
            gene_novelty.push(make_annotation(gid, &refs.config, g_nov, "TRUE"));

            let tm = create_novel_transcript(
                strand, chrom, positions[0], *positions.last().unwrap(),
                gid, &full_edges, &full_vertices, &refs.counters, novel_buffer, novel_tx_registry,
            );
            transcript_novelty.push(make_annotation(tm.transcript_id, &refs.config, t_nov, "TRUE"));
            (gid, tm.transcript_id)
        }
        Some((overlap_gid, match_strand)) if match_strand != strand => {
            let gid = get_or_create_novel_gene(
                NovelGeneKey::Antisense(overlap_gid, strand),
                chrom, positions[0], *positions.last().unwrap(), strand,
                &refs.counters, &refs.novel_gene_dict,
            );
            gene_novelty.push(make_annotation(gid, &refs.config, "antisense_gene", "TRUE"));
            gene_novelty.push(make_annotation(gid, &refs.config, "gene_antisense_to_IDs", &overlap_gid.to_string()));

            let tm = create_novel_transcript(
                strand, chrom, positions[0], *positions.last().unwrap(),
                gid, &full_edges, &full_vertices, &refs.counters, novel_buffer, novel_tx_registry,
            );
            transcript_novelty.push(make_annotation(tm.transcript_id, &refs.config, "antisense_transcript", "TRUE"));
            (gid, tm.transcript_id)
        }
        Some((gid, _)) => {
            let tm = create_novel_transcript(
                strand, chrom, positions[0], *positions.last().unwrap(),
                gid, &full_edges, &full_vertices, &refs.counters, novel_buffer, novel_tx_registry,
            );
            transcript_novelty.push(make_annotation(tm.transcript_id, &refs.config, "genomic_transcript", "TRUE"));
            (gid, tm.transcript_id)
        }
    };

    (gene_id, transcript_id, gene_novelty, transcript_novelty, sei)
}

// ─── Helper functions ────────────────────────────────────────────────────────

fn make_annotation(id: i64, config: &RunConfig, attribute: &str, value: &str) -> NoveltyAnnotation {
    NoveltyAnnotation {
        id,
        annot_name: config.annot_pool.annot_name(),
        source: config.annot_pool.source(),
        attribute: config.annot_pool.get(attribute),
        value: config.annot_pool.get(value),
    }
}

/// Key an intergenic/fusion read onto a shared novel locus, binned by the
/// read's leftmost coordinate.
fn locus_key(chrom: &str, start: i64, end: i64, strand: Strand, bin: i64) -> NovelGeneKey {
    let anchor = start.min(end);
    let binned = if bin > 0 { anchor.div_euclid(bin) } else { anchor };
    NovelGeneKey::Locus(chrom.to_string(), strand, binned)
}

/// Resolve the novel gene for this read, creating one only if the locus is
/// new; reads on the same key share the gene and widen its span.
fn get_or_create_novel_gene(
    key: NovelGeneKey,
    chrom: &str,
    start: i64,
    end: i64,
    strand: Strand,
    counters: &Counters,
    novel_gene_dict: &ConcurrentGeneDict,
) -> GeneId {
    let lo = start.min(end);
    let hi = start.max(end);

    match novel_gene_dict.entry(key) {
        dashmap::mapref::entry::Entry::Occupied(mut occupied) => {
            let entry = occupied.get_mut();
            entry.start = entry.start.min(lo);
            entry.end = entry.end.max(hi);
            entry.gene_id
        }
        dashmap::mapref::entry::Entry::Vacant(vacant) => {
            let gene_id = counters.gene.increment();
            vacant.insert(NovelGeneEntry {
                gene_id,
                chromosome: chrom.to_string(),
                strand,
                start: lo,
                end: hi,
            });
            gene_id
        }
    }
}

fn create_novel_transcript(
    _strand: Strand,
    chrom: &str,
    start: i64,
    end: i64,
    gene_id: GeneId,
    edge_ids: &[EdgeId],
    vertex_ids: &[VertexId],
    counters: &Counters,
    novel_buffer: &mut NovelBuffer,
    novel_tx_registry: &NovelTranscriptRegistry,
) -> TranscriptModel {
    // ── Dedup key ──────────────────────────────────────────
    // Multi-exon transcripts (3+ edges) are keyed on internal edges only, so
    // reads with the same internal splice structure share one model whatever
    // their terminal positions; shorter paths use the full edge path.
    let dedup_key = if edge_ids.len() > 2 {
        // Internal edges: skip first exon edge and last exon edge
        TranscriptPathKey::new(edge_ids[1..edge_ids.len() - 1].to_vec())
    } else {
        // Shorter paths: use full path
        TranscriptPathKey::new(edge_ids.to_vec())
    };

    // Full path key for the transcript model storage
    let full_path_key = TranscriptPathKey::new(edge_ids.to_vec());

    // ── Reuse the ID of an already registered internal structure ─────────────
    let tid = if let Some(existing) = novel_tx_registry.get(&dedup_key) {
        *existing
    } else {
        let new_id = counters.transcript.increment();
        let actual_id = *novel_tx_registry.entry(dedup_key.clone()).or_insert(new_id);
        if actual_id == new_id {
            // First registration of this structure
            let jn = if edge_ids.len() > 1 {
                let internal: Vec<String> = edge_ids[1..edge_ids.len() - 1]
                    .iter().map(|e| e.to_string()).collect();
                if internal.is_empty() { None } else { Some(internal.join(",")) }
            } else { None };
            novel_buffer.transcripts.push((full_path_key.clone(), TranscriptModel {
                transcript_id: new_id,
                gene_id,
                start_exon: edge_ids[0],
                jn_path: jn,
                end_exon: *edge_ids.last().unwrap(),
                start_vertex: vertex_ids[0],
                end_vertex: *vertex_ids.last().unwrap(),
                n_exons: ((edge_ids.len() + 1) / 2) as i32,
                chromosome: chrom.to_string(),
                start_pos: start,
                end_pos: end,
            }));
        }
        actual_id
    };

    let jn_path = if edge_ids.len() > 1 {
        let internal: Vec<String> = edge_ids[1..edge_ids.len() - 1]
            .iter()
            .map(|e| e.to_string())
            .collect();
        if internal.is_empty() { None } else { Some(internal.join(",")) }
    } else {
        None
    };

    let n_exons = ((edge_ids.len() + 1) / 2) as i32;

    TranscriptModel {
        transcript_id: tid,
        gene_id,
        start_exon: edge_ids[0],
        jn_path,
        end_exon: *edge_ids.last().unwrap(),
        start_vertex: vertex_ids[0],
        end_vertex: *vertex_ids.last().unwrap(),
        n_exons,
        chromosome: chrom.to_string(),
        start_pos: start,
        end_pos: end,
    }
}

fn build_full_edge_novelty(splice_novelty: &[u8], sei: &StartEndInfo) -> Vec<u8> {
    let mut full = Vec::with_capacity(splice_novelty.len() + 2);
    full.push(sei.start_novelty);
    full.extend_from_slice(splice_novelty);
    full.push(sei.end_novelty);
    full
}

fn collect_exon_edges(edge_ids: &[EdgeId], novelty: &[u8]) -> Vec<(EdgeId, bool)> {
    edge_ids
        .iter()
        .zip(novelty.iter())
        .step_by(2) // Exons are at even indices (0, 2, 4, ...)
        .map(|(&eid, &nov)| (eid, nov == 1))
        .collect()
}

/// Assign a read to an annotated transcript whose chain it matches within a
/// tolerance, when exact edge-ID matching has failed.
///
/// The read takes the annotated chain; its terminals go through the same
/// cutoff logic as `process_fsm`.
fn process_approximate_fsm(
    chrom: &str,
    positions: &[i64],
    strand: Strand,
    refs: &SharedRefs,
    novel_buffer: &mut NovelBuffer,
    stage: &str,
    ambiguous_sites: u64,
    settled_sites: u64,
) -> Option<(GeneId, TranscriptId, Vec<NoveltyAnnotation>, StartEndInfo)> {
    let index = refs.approx_chain.as_ref()?;
    let resolution = index.resolve_with_unlocked(chrom, strand, positions, refs.config.approx_chain_tolerance,
                                                 ambiguous_sites, settled_sites);
    if refs.config.chain_trace {
        novel_buffer.chain_traces.push(resolution.trace_row(stage, chrom, strand, positions));
    }
    let hit = resolution.matched?;
    // The index holds the flattened annotated chain.
    let model = hit.model;
    let internal_edges = &model.internal_edges;
    let internal_vertices = &model.internal_vertices;

    let gene_id = model.gene_id;
    let curr_5p_diff = compute_delta(model.start_pos, positions[0], strand);
    let curr_3p_diff = compute_delta(model.end_pos, *positions.last().unwrap(), strand);
    let mut sei = StartEndInfo::default();

    let novel_5p = curr_5p_diff.abs() > refs.config.cutoff_5p;
    if !novel_5p {
        sei.start_vertex = model.start_vertex;
        sei.start_exon = model.start_exon;
        sei.diff_5p = Some(curr_5p_diff);
        sei.start_novelty = 0;
    } else {
        let (sv, se, sn, _, d5) = vertex_matching::process_5p(
            &refs.location_dict, &refs.edge_dict, chrom, positions, strand,
            internal_vertices, gene_id, &refs.gene_starts, &refs.counters, &refs.config,
            novel_buffer, &refs.novel_location_dict, &refs.novel_edge_dict,
        );
        sei.start_vertex = sv;
        sei.start_exon = se;
        sei.start_novelty = sn;
        sei.diff_5p = d5;
    }

    let novel_3p = curr_3p_diff.abs() > refs.config.cutoff_3p;
    if !novel_3p {
        sei.end_vertex = model.end_vertex;
        sei.end_exon = model.end_exon;
        sei.diff_3p = Some(curr_3p_diff);
        sei.end_novelty = 0;
    } else {
        let (ev, ee, en, _, d3) = vertex_matching::process_3p(
            &refs.location_dict, &refs.edge_dict, chrom, positions, strand,
            internal_vertices, gene_id, &refs.gene_ends, &refs.counters, &refs.config,
            novel_buffer, &refs.novel_location_dict, &refs.novel_edge_dict,
        );
        sei.end_vertex = ev;
        sei.end_exon = ee;
        sei.end_novelty = en;
        sei.diff_3p = d3;
    }

    let mut full_edges = vec![sei.start_exon];
    full_edges.extend_from_slice(internal_edges);
    full_edges.push(sei.end_exon);
    let mut full_vertices = vec![sei.start_vertex];
    full_vertices.extend_from_slice(internal_vertices);
    full_vertices.push(sei.end_vertex);
    sei.edge_ids = full_edges;
    sei.vertex_ids = full_vertices;

    // Record the assigned model and the largest junction shift.
    let tid = model.transcript_id;
    let mut novelty = vec![
        make_annotation(tid, &refs.config, "approximate_FSM", "TRUE"),
        make_annotation(tid, &refs.config, "approximate_FSM_to_IDs", &tid.to_string()),
        make_annotation(tid, &refs.config, "approximate_FSM_max_offset",
                        &hit.max_offset.to_string()),
    ];
    if novel_5p {
        novelty.push(make_annotation(tid, &refs.config, "FSM_novel_5p", "TRUE"));
    }
    if novel_3p {
        novelty.push(make_annotation(tid, &refs.config, "FSM_novel_3p", "TRUE"));
    }
    Some((gene_id, tid, novelty, sei))
}
