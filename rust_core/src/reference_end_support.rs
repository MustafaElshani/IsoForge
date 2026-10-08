//! Adds clusters at annotated transcript ends that have read support but no called peak nearby.
use crate::annotation_engine::ObsCore;
use crate::cage_peaks::CagePeakIndex;
use crate::data_structures::SharedRefs;
use crate::errors::IsoForgeResult;
use crate::genome_ref::GenomeRef;
use crate::tss_tes_clustering::{
    ClusterIndex, FeatureAnnotationTree, ReadFilterHints, TssTesEngine,
};
use crate::types::{ClusterType, EndSiteCluster, Strand, HINT_IS_KNOWN, HINT_IS_ISM, HINT_ISM_5P_INTACT, HINT_ISM_3P_INTACT};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::io::Write;

type Site = (String, Strand, u8, i64);

/// Adds a cluster at each annotated transcript end that has no called peak
/// within 25 bp and enough reads of the transcript ending within
/// `reference_end_window` of it. Writes the added clusters to `output`.
pub fn supplement(
    engine: &mut TssTesEngine,
    refs: &SharedRefs,
    cores: &[ObsCore],
    hints: Option<&ReadFilterHints>,
    cage: Option<&CagePeakIndex>,
    features: &FeatureAnnotationTree,
    genome: Option<&GenomeRef>,
    output: &std::path::Path,
) -> IsoForgeResult<()> {
    let positions: HashMap<_, _> = refs
        .location_dict
        .chromosomes
        .values()
        .flat_map(|c| c.positions.iter().map(|&(p, v)| (v, p)))
        .collect();
    let transcripts: HashMap<_, _> = refs
        .transcript_dict
        .values()
        .map(|t| (t.transcript_id, t))
        .collect();
    let mut sites: BTreeMap<Site, Vec<&ObsCore>> = BTreeMap::new();
    for c in cores {
        if c.hint_flags & HINT_IS_KNOWN == 0
            || c.is_internal_priming
            || c.fraction_as.unwrap_or(0.0) > refs.config.max_frac_as
        {
            continue;
        }
        // The tail gate follows the run's poly(A) policy.
        if !tail_gate(engine.config.require_polya_tail, c.polya_tail_len, c.polya_tail_softclip) {
            continue;
        }
        // A full-length read counts for both annotated ends or for neither; a
        // truncated read only for the end it kept intact, when the run allows it.
        let ism = c.hint_flags & HINT_IS_ISM != 0;
        if ism && !engine.config.ism_intact_end_support {
            continue;
        }
        let (mut ok_5p, mut ok_3p) = (true, true);
        if let Some(h) = hints {
            if h.is_degraded_5p.contains(c.read_idx) || h.exclude_tss_evidence.contains(c.read_idx) {
                ok_5p = false;
            }
            if h.is_degraded_3p.contains(c.read_idx) || h.exclude_tes_evidence.contains(c.read_idx) {
                ok_3p = false;
            }
        }
        if cage.is_some_and(|idx| {
            !idx.overlaps_with_window(
                &c.location,
                c.strand,
                c.start_pos,
                engine.config.cage_window,
            )
        }) {
            ok_5p = false;
        }
        if ism {
            ok_5p &= c.hint_flags & HINT_ISM_5P_INTACT != 0;
            ok_3p &= c.hint_flags & HINT_ISM_3P_INTACT != 0;
        } else if !(ok_5p && ok_3p) {
            continue;
        }
        let Some(t) = transcripts.get(&c.transcript_id) else {
            continue;
        };
        let (Some(&start), Some(&end)) =
            (positions.get(&t.start_vertex), positions.get(&t.end_vertex))
        else {
            continue;
        };
        let window = engine.config.reference_end_window;
        let near_5p = (c.start_pos - start).abs() <= window;
        let near_3p = (c.end_pos - end).abs() <= window;
        if !ism && !(near_5p && near_3p) {
            continue;
        }
        if ok_5p && near_5p {
            sites
                .entry((c.location.to_string(), c.strand, 0, start))
                .or_default()
                .push(c);
        }
        if ok_3p && near_3p {
            sites
                .entry((c.location.to_string(), c.strand, 1, end))
                .or_default()
                .push(c);
        }
    }
    // A site needs two reads, or one with `known_singleton_ends`.
    let min_reads = if engine.config.known_singleton_ends { 1 } else { 2 };
    let mut proposals = Vec::new();
    for (key, mut rs) in sites {
        rs.sort_by_key(|r| (if key.2 == 0 { r.start_pos } else { r.end_pos }, r.read_idx));
        let raw: Vec<_> = rs
            .iter()
            .map(|r| if key.2 == 0 { r.start_pos } else { r.end_pos })
            .collect();
        let Some((lo, hi)) = supported_window(&raw, min_reads) else {
            continue;
        };
        let selected = rs[lo..hi].to_vec();
        proposals.push((key, selected, raw[lo], raw[hi - 1], raw[(lo + hi - 1) / 2]));
    }
    proposals.sort_by_key(|(key, rs, _, _, _)| (std::cmp::Reverse(rs.len()), key.clone()));
    let mut occupied: BTreeMap<(String, Strand, u8), BTreeSet<i64>> = BTreeMap::new();
    for (kind, clusters) in [(0, &engine.tss_clusters), (1, &engine.tes_clusters)] {
        for c in clusters {
            occupied
                .entry((c.chromosome.clone(), c.strand, kind))
                .or_default()
                .insert(c.peak_site);
        }
    }
    let mut next_id = engine
        .tss_clusters
        .iter()
        .chain(&engine.tes_clusters)
        .map(|c| c.cluster_id)
        .max()
        .unwrap_or(0);
    let mut audit = std::io::BufWriter::new(std::fs::File::create(output)?);
    writeln!(
        audit,
        "cluster_id\tchrom\tstrand\ttype\treference_site\tobserved_peak\tindependent_reads\tfull_length_reads\tintact_end_reads"
    )?;
    let n_datasets = engine.counts.dataset_read_counts.len() as u32;
    for ((chrom, strand, kind, reference), rs, start, end, peak) in proposals {
        let occupied = occupied.entry((chrom.clone(), strand, kind)).or_default();
        if occupied
            .range(reference - 25..=reference + 25)
            .next()
            .is_some()
        {
            continue;
        }
        let mut datasets = HashMap::new();
        for r in &rs {
            *datasets.entry(r.dataset.to_string()).or_insert(0usize) += 1;
        }
        let supporting = datasets.values().filter(|&&n| n >= min_reads).count() as u32;
        if supporting == 0 {
            continue;
        }
        if n_datasets >= engine.config.min_replication_samples {
            let minimum = (n_datasets as f64 * engine.config.min_replication_rate).ceil() as u32;
            if supporting < minimum
                && datasets.values().copied().max().unwrap_or(0)
                    < engine.config.min_condition_reads as usize
            {
                continue;
            }
        }
        let kind = if kind == 0 {
            ClusterType::TSS
        } else {
            ClusterType::TES
        };
        let (feature, gene_id) = features.annotate(&chrom, peak, strand, kind);
        let (_, pas, distance) = if kind == ClusterType::TES {
            genome
                .map(|g| g.scan_pas_upstream(&chrom, peak, strand))
                .unwrap_or((false, None, None))
        } else {
            (false, None, None)
        };
        next_id += 1;
        let count = rs.len() as u64;
        let intact = rs.iter().filter(|r| r.hint_flags & HINT_IS_ISM != 0).count() as u64;
        // A single-read cluster is recorded in `singleton_anchored`.
        if count < 2 {
            engine.singleton_anchored.insert(next_id);
        }
        let cluster = EndSiteCluster {
            cluster_id: next_id,
            cluster_type: kind,
            chromosome: chrom.clone(),
            start,
            end,
            peak_site: peak,
            strand,
            read_count: count,
            feature,
            gene_id,
            tpm: count as f64 / cores.len().max(1) as f64 * 1e6,
            gene_count: 0,
            usage: 0.0,
            frac_a: if kind == ClusterType::TES {
                Some(rs.iter().map(|r| r.fraction_as.unwrap_or(0.0)).sum::<f64>() / count as f64)
            } else {
                None
            },
            pas_signal: pas,
            pas_distance: distance,
            n_datasets_supporting: supporting,
        };
        writeln!(
            audit,
            "{next_id}\t{chrom}\t{}\t{}\t{reference}\t{peak}\t{count}\t{}\t{intact}",
            strand.as_str(),
            if kind == ClusterType::TSS {
                "TSS"
            } else {
                "TES"
            },
            count - intact
        )?;
        occupied.insert(peak);
        if kind == ClusterType::TSS {
            engine.tss_clusters.push(cluster);
        } else {
            engine.tes_clusters.push(cluster);
        }
    }
    engine.compute_usage_stats();
    engine.tss_index = Some(ClusterIndex::build(&engine.tss_clusters));
    engine.tes_index = Some(ClusterIndex::build(&engine.tes_clusters));
    Ok(())
}

/// Whether a read satisfies the poly(A) requirement in force for this run.
/// Evidence is a basecaller tail length or a detected soft-clipped tail.
fn tail_gate(require: bool, tail_len: Option<i32>, softclip: Option<i32>) -> bool {
    if !require {
        return true;
    }
    tail_len.is_some_and(|n| n > 0) || softclip.is_some_and(|n| n > 0)
}

/// Half-open index range of the largest run of sorted positions spanning at
/// most 20 bp, or None when it holds fewer than `min_reads`. Ties resolve to
/// the first window.
fn supported_window(sorted: &[i64], min_reads: usize) -> Option<(usize, usize)> {
    let mut left = 0;
    let mut best = (0, 0);
    for right in 0..sorted.len() {
        while sorted[right] - sorted[left] > 20 {
            left += 1;
        }
        if right + 1 - left > best.1 - best.0 {
            best = (left, right + 1);
        }
    }
    (best.1 - best.0 >= min_reads.max(1)).then_some(best)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tail_gate_follows_the_run_policy() {
        // Required: needs either evidence source, and a zero-length tag is absence.
        assert!(tail_gate(true, Some(20), None));
        assert!(tail_gate(true, None, Some(15)));
        assert!(!tail_gate(true, Some(0), None));
        assert!(!tail_gate(true, None, None));
        // Relaxed: never blocks, regardless of evidence.
        assert!(tail_gate(false, None, None));
        assert!(tail_gate(false, Some(0), Some(0)));
    }

    #[test]
    fn no_singletons_or_transitive_end_support() {
        assert_eq!(supported_window(&[], 2), None);
        assert_eq!(supported_window(&[100], 2), None);
        assert_eq!(supported_window(&[100, 105], 2), Some((0, 2)));
        assert_eq!(supported_window(&[100, 130], 2), None);
        assert_eq!(supported_window(&[100, 120, 140], 2), Some((0, 2)));
    }

    /// With a floor of one read an empty site gives None and positions more
    /// than 20 bp apart are not joined into one window.
    #[test]
    fn known_singleton_floor_admits_one_read_but_not_a_staircase() {
        assert_eq!(supported_window(&[], 1), None);
        assert_eq!(supported_window(&[100], 1), Some((0, 1)));
        assert_eq!(supported_window(&[100, 105], 1), Some((0, 2)));
        // 100/120/140 spans 40 bp: the largest window holds two reads.
        assert_eq!(supported_window(&[100, 120, 140], 1), Some((0, 2)));
    }
}
