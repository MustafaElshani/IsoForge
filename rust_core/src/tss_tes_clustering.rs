// TSS/TES clustering engine.
//
//   1. COUNT: accumulate 5' (TSS) and 3' (TES) positions per chromosome/strand
//   2. CLUSTER: Gaussian-smoothed peak calling
//   3. ANNOTATE: assign each cluster to UTR/exon/intron/gene features
//   4. PAS SCAN: scan upstream of TES clusters for poly(A) signal hexamers
//   5. JOIN: assign each read's 5'/3' ends to the nearest TSS/TES cluster
//   6. FILTER: mark reads failing the evidence gates as unsupported

// HashMap and HashSet with deterministic hashing.
use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};
use std::path::Path;
use std::sync::Arc;

use coitrees::{COITree, Interval, IntervalTree};
use rayon::prelude::*;

use crate::bam_reader::ClusteringInput;
use crate::cage_peaks::CagePeakIndex;
use crate::errors::*;
use crate::genome_ref::GenomeRef;
use crate::types::*;

// ═══════════════════════════════════════════════════════════════════════════════
// STEP 1: COUNTING — Accumulate 5'/3' end positions
// ═══════════════════════════════════════════════════════════════════════════════

/// Position counts for a single chromosome/strand.
/// Key: genomic position (1-based), Value: (total_count, per_dataset_counts)
type PositionCounts = HashMap<i64, (u64, HashMap<String, u64>)>;

/// With ISOFORGE_TSS_GATE_TRACE set to a path, every TSS peak the gate judges is written there with its window
/// statistics, the 5' ends in the 200 bp downstream and upstream of its window and whether it was kept.
static TSS_GATE_TRACE: std::sync::Mutex<Vec<String>> = std::sync::Mutex::new(Vec::new());
const TSS_TRACE_FLANK: i64 = 200;
/// Downstream bases whose 5' ends the TSS rescue weighs a rejected peak against.
const TSS_RESCUE_FLANK: i64 = 200;

/// Read-level filter hints for novelty-aware counting, built after transcript
/// assignment and keyed by `ReadIdx`.
#[derive(Debug, Clone, Default)]
pub struct ReadFilterHints {
    /// Reads flagged as internal priming artifacts
    pub is_internal_priming: ReadIdxSet,
    /// Reads whose transcript-oriented 5′ portion is missing (ISM-suffix).
    pub exclude_tss_evidence: ReadIdxSet,
    /// Reads whose transcript-oriented 3′ portion is missing (ISM-prefix).
    pub exclude_tes_evidence: ReadIdxSet,
    /// Reads whose 5′ end is degraded (deviates from per-transcript modal position)
    pub is_degraded_5p: ReadIdxSet,
    /// Reads whose 3′ end is degraded (deviates from per-transcript modal position)
    pub is_degraded_3p: ReadIdxSet,
    /// Reads classified as Known/FSM; exempt from weak-cluster revalidation.
    pub is_known_read: ReadIdxSet,
}

impl ReadFilterHints {
    pub fn mark_truncation(&mut self, read_idx: ReadIdx, flags: u8) {
        use crate::types::{HINT_IS_ISM, HINT_IS_ISM_PREFIX, HINT_IS_ISM_SUFFIX};
        if flags & HINT_IS_ISM == 0 { return; }
        if flags & HINT_IS_ISM_PREFIX != 0 { self.exclude_tes_evidence.insert(read_idx); }
        if flags & HINT_IS_ISM_SUFFIX != 0 { self.exclude_tss_evidence.insert(read_idx); }
    }
}


/// Accumulated end-site counts across all reads.
#[derive(Debug, Default)]
pub struct EndSiteCounts {
    /// TSS counts: (chrom, strand) → position → (count, dataset_counts)
    pub tss: HashMap<(String, Strand), PositionCounts>,
    /// TES counts: (chrom, strand) → position → (count, dataset_counts)
    pub tes: HashMap<(String, Strand), PositionCounts>,
    /// Per-read fraction-A values keyed by 3' position for cluster-level aggregation
    pub frac_a_by_pos: HashMap<(String, Strand, i64), Vec<f64>>,
    /// Total reads per dataset
    pub dataset_read_counts: HashMap<String, u64>,
}

impl EndSiteCounts {
    pub fn new() -> Self {
        Self::default()
    }

    /// Accumulate counts from a batch of annotated reads: the 5' end goes to the
    /// TSS counts and the 3' end to the TES counts, unless one of the filters
    /// below excludes the read from that end.
    pub fn add_reads<R: ClusteringInput>(
        &mut self,
        reads: &[R],
        config: &ClusterConfig,
        hints: Option<&ReadFilterHints>,
        cage_index: Option<&CagePeakIndex>,
    ) {
        let mut n_ip_filtered = 0u64;
        let mut n_tss_truncated_filtered = 0u64;
        let mut n_tes_truncated_filtered = 0u64;
        let mut n_degraded_5p_filtered = 0u64;
        let mut n_degraded_3p_filtered = 0u64;
        let mut n_no_polya_filtered = 0u64;
        let mut n_no_cage_filtered = 0u64;
        let mut n_no_tail_start_filtered = 0u64;
        let mut n_unaligned_3p_filtered = 0u64;

        for read in reads {
            let pos_5p = match read.pos_5p() {
                Some(p) => p,
                None => continue,
            };
            let pos_3p = match read.pos_3p() {
                Some(p) => p,
                None => continue,
            };

            // Check internal priming status
            let is_ip = if config.filter_internal_priming {
                // Hints from annotation first
                if let Some(h) = hints {
                    h.is_internal_priming.contains(read.read_idx())
                } else {
                    // Otherwise the read's own internal_priming field
                    read.is_ip()
                }
            } else {
                false
            };

            // ISM-prefix/suffix status (from hints)
            let exclude_tss_evidence = if config.filter_truncated_reads {
                hints.map_or(false, |h| h.exclude_tss_evidence.contains(read.read_idx()))
            } else {
                false
            };
            let exclude_tes_evidence = if config.filter_truncated_reads {
                hints.map_or(false, |h| h.exclude_tes_evidence.contains(read.read_idx()))
            } else {
                false
            };

            // Degradation status
            let is_degraded_5p = if config.filter_degraded_reads {
                hints.map_or(false, |h| h.is_degraded_5p.contains(read.read_idx()))
            } else {
                false
            };
            let is_degraded_3p = if config.filter_degraded_reads {
                hints.map_or(false, |h| h.is_degraded_3p.contains(read.read_idx()))
            } else {
                false
            };

            let key_tss = (read.chromosome().to_owned(), read.strand_val());
            let key_tes = (read.chromosome().to_owned(), read.strand_val());

            // Reads without a poly(A) tail (pt:i:0 counts as none) are excluded from TES counting
            let has_polya = read.polya_len().map_or(false, |v| v > 0)
                         || read.polya_softclip().map_or(false, |v| v > 0);
            let no_polya_filter = config.require_polya_tail && !has_polya;
            // A spliced read without its tail places no start when starts come from
            // full-length reads.
            let no_tail_start = config.tss_from_full_length && read.n_exons() > 1 && !has_polya;

            // Reads whose 5' end misses CAGE peaks are excluded from TSS counting
            let no_cage_filter = if let Some(cage) = cage_index {
                !cage.overlaps_with_window(read.chromosome(), read.strand_val(), pos_5p, config.cage_window)
            } else {
                false // no CAGE index → no CAGE filter
            };

            // TSS count: skipped for 5'-truncated (ISM-suffix), 5'-degraded or
            // CAGE-unsupported reads.
            
            if !exclude_tss_evidence && !is_degraded_5p && !no_cage_filter && !no_tail_start {
                let tss_entry = self.tss.entry(key_tss).or_default();
                let (total, per_ds) = tss_entry.entry(pos_5p).or_insert_with(|| (0, HashMap::default()));
                *total += 1;
                *per_ds.entry(read.dataset_name().to_owned()).or_insert(0) += 1;
            } else {
                if exclude_tss_evidence { n_tss_truncated_filtered += 1; }
                if is_degraded_5p { n_degraded_5p_filtered += 1; }
                if no_cage_filter { n_no_cage_filtered += 1; }
                if no_tail_start { n_no_tail_start_filtered += 1; }
            }

            // TES count: skipped for internally primed, 3'-truncated (ISM-prefix),
            // 3'-degraded or tail-less reads.
            
            // A read whose alignment stops before its 3' end (transcript sequence before the tail) places no TES;
            // it still joins a cluster other reads make.
            let unaligned_3p = read.unaligned_3p();
            if !is_ip && !exclude_tes_evidence && !is_degraded_3p && !no_polya_filter && !unaligned_3p {
                let tes_entry = self.tes.entry(key_tes).or_default();
                let (total, per_ds) = tes_entry.entry(pos_3p).or_insert_with(|| (0, HashMap::default()));
                *total += 1;
                *per_ds.entry(read.dataset_name().to_owned()).or_insert(0) += 1;
            } else {
                if is_ip { n_ip_filtered += 1; }
                if exclude_tes_evidence { n_tes_truncated_filtered += 1; }
                if is_degraded_3p { n_degraded_3p_filtered += 1; }
                if no_polya_filter { n_no_polya_filtered += 1; }
                if unaligned_3p { n_unaligned_3p_filtered += 1; }
            }

            // Track fraction-A at 3' positions for cluster-level aggregation
            if let Some(frac_a) = read.frac_a() {
                self.frac_a_by_pos
                    .entry((read.chromosome().to_owned(), read.strand_val(), pos_3p))
                    .or_default()
                    .push(frac_a);
            }

            // Dataset read counts
            *self.dataset_read_counts.entry(read.dataset_name().to_owned()).or_insert(0) += 1;
        }

        if n_ip_filtered > 0 || n_tss_truncated_filtered > 0 || n_tes_truncated_filtered > 0
            || n_unaligned_3p_filtered > 0 || n_degraded_5p_filtered > 0 || n_degraded_3p_filtered > 0
            || n_no_polya_filtered > 0 || n_no_cage_filtered > 0 || n_no_tail_start_filtered > 0
        {
            eprintln!(
                "[IsoForge]   Filtered from counting: {} internal-priming (TES), {} 5prime-truncated (TSS), {} 3prime-truncated (TES), {} degraded-5p (TSS), {} degraded-3p (TES), {} no-polyA (TES), {} no-CAGE (TSS), {} spliced without tail (TSS), {} unaligned before the tail (TES)",
                n_ip_filtered, n_tss_truncated_filtered, n_tes_truncated_filtered,
                n_degraded_5p_filtered, n_degraded_3p_filtered,
                n_no_polya_filtered, n_no_cage_filtered, n_no_tail_start_filtered, n_unaligned_3p_filtered
            );
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// STEP 2: CLUSTERING — Gaussian-smoothed peak calling
// ═══════════════════════════════════════════════════════════════════════════════

/// Compact start seeds: windows of at most 6 bp holding at least three starts
/// (or `minimum`), whose density is at least four times that of the
/// surrounding 50 bp.
///
/// Returns: Vec of (start, end, peak_position, total_count, per_dataset_counts)
fn compact_start_candidates(counts: &PositionCounts, minimum: u64)
    -> Vec<(i64, i64, i64, u64, HashMap<String, u64>)> {
    let mut sites: Vec<_> = counts.keys().copied().collect(); sites.sort_unstable();
    let mut prefix = vec![0u64];
    for p in &sites { prefix.push(prefix.last().unwrap()+counts[p].0); }
    let mut proposals = Vec::new();
    for (i, &start) in sites.iter().enumerate() {
        let end_idx = sites.partition_point(|&p| p <= start+6);
        let total = prefix[end_idx]-prefix[i];
        if total < minimum.max(3) { continue; }
        let end = sites[end_idx-1];
        let lo = sites.partition_point(|&p| p < start-50);
        let hi = sites.partition_point(|&p| p <= end+50);
        let outside = prefix[hi]-prefix[lo]-total;
        if total*100 < 4*outside*((end-start+1) as u64) { continue; }
        // Weighted median, independent of iteration order.
        let mid = (i..end_idx).find(|&j| 2*(prefix[j+1]-prefix[i]) >= total).unwrap();
        let mut ds = HashMap::default();
        for p in &sites[i..end_idx] {
            for (d, n) in &counts[p].1 { *ds.entry(d.clone()).or_default() += n; }
        }
        proposals.push((start, end, sites[mid], total, ds));
    }
    proposals.sort_by_key(|p| (std::cmp::Reverse(p.3), p.0, p.2));
    let mut result: Vec<(i64, i64, i64, u64, HashMap<String, u64>)> = Vec::new();
    for p in proposals {
        if !result.iter().any(|q| p.0 <= q.1+6 && q.0 <= p.1+6) { result.push(p); }
    }
    result.sort_by_key(|p| p.0);
    result
}

/// Attrition through the peak-calling gates.
#[derive(Default, Debug, Clone, Copy)]
pub struct PeakGateStats {
    /// Local maximum whose positions were entirely claimed by a taller peak.
    pub suppressed_by_neighbour: u64,
    pub below_min_cluster_reads: u64,
    pub rejected_entropy: u64,
    pub rejected_peak_ratio: u64,
    pub rejected_width: u64,
    /// Rejected peaks kept by the rescue (a gene's strongest start, above every kept TSS peak of the gene).
    pub rescued: u64,
}

impl PeakGateStats {
    fn merge(&mut self, o: &PeakGateStats) {
        self.suppressed_by_neighbour += o.suppressed_by_neighbour;
        self.below_min_cluster_reads += o.below_min_cluster_reads;
        self.rejected_entropy += o.rejected_entropy;
        self.rejected_peak_ratio += o.rejected_peak_ratio;
        self.rejected_width += o.rejected_width;
        self.rescued += o.rescued;
    }
    fn total(&self) -> u64 {
        self.suppressed_by_neighbour + self.below_min_cluster_reads
            + self.rejected_entropy + self.rejected_peak_ratio + self.rejected_width
    }
}

/// Call peaks from position counts using Gaussian smoothing: each local maximum of the smoothed counts, tallest
/// first, is extended outward until the smoothed count drops below extent_cutoff or ratio_cutoff × peak height.
///
/// Returns: Vec of (start, end, peak_position, total_count, per_dataset_counts)
fn call_peaks(
    counts: &PositionCounts,
    config: &ClusterConfig,
    min_position_reads: u64,
    stats: &mut PeakGateStats,
) -> Vec<(i64, i64, i64, u64, HashMap<String, u64>)> {
    if counts.is_empty() {
        return Vec::new();
    }
    let pad = config.cluster_window as i64;
    let lo = *counts.keys().min().unwrap() - pad;
    let hi = *counts.keys().max().unwrap() + pad;
    // Positions with fewer reads than the threshold are left out of the smoothing input.
    let mut sites: Vec<(i64, u64)> = counts
        .iter()
        .filter(|(_, v)| v.0 >= min_position_reads)
        .map(|(&p, v)| (p, v.0))
        .collect();
    sites.sort_unstable();

    let mut clusters = if config.extent_cutoff > 0 || config.ratio_cutoff > 0.0 {
        // Runs of positions more than twice `reach` apart are called on arrays of their own.
        let reach = pad + config.gaussian_std.ceil() as usize as i64 + 1;
        let mut out = Vec::new();
        let mut first = 0;
        while first < sites.len() {
            let mut last = first;
            while last + 1 < sites.len() && sites[last + 1].0 - sites[last].0 <= 2 * reach {
                last += 1;
            }
            let span = ((sites[first].0 - reach).max(lo), (sites[last].0 + reach).min(hi));
            out.extend(peaks_on_span(counts, config, &sites[first..=last], span, stats));
            first = last + 1;
        }
        out
    } else {
        // With both extent thresholds at zero, all positions are called on one array.
        peaks_on_span(counts, config, &sites, (lo, hi), stats)
    };

    // Sort by position
    clusters.sort_by_key(|c| c.0);
    clusters
}

/// Peaks of the counted `sites` (sorted by position) on the bases `span.0..=span.1`, found as
/// `call_peaks` describes, in the order they are claimed.
fn peaks_on_span(
    counts: &PositionCounts,
    config: &ClusterConfig,
    sites: &[(i64, u64)],
    span: (i64, i64),
    stats: &mut PeakGateStats,
) -> Vec<(i64, i64, i64, u64, HashMap<String, u64>)> {
    let window = config.cluster_window;
    let arr_len = (span.1 - span.0 + 1) as usize;
    let mut coverage = vec![0.0f64; arr_len];
    for &(pos, count) in sites {
        coverage[(pos - span.0) as usize] = count as f64;
    }

    // Apply Gaussian smoothing
    let smoothed = gaussian_smooth(&coverage, config.gaussian_std, window);

    // Find local maxima (peaks)
    let mut peaks: Vec<(usize, f64)> = Vec::new();
    for i in 1..smoothed.len().saturating_sub(1) {
        if smoothed[i] > smoothed[i - 1] && smoothed[i] >= smoothed[i + 1] && smoothed[i] > 0.0 {
            peaks.push((i, smoothed[i]));
        }
    }

    // Sort peaks by height (descending) for greedy assignment
    peaks.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

    // For each peak, extend to find cluster boundaries
    let mut used = vec![false; arr_len];
    let mut clusters: Vec<(i64, i64, i64, u64, HashMap<String, u64>)> = Vec::new();

    for (peak_idx, peak_height) in &peaks {
        if used[*peak_idx] {
            continue;
        }

        let extent_threshold = (config.extent_cutoff as f64).max(config.ratio_cutoff * peak_height);

        // Minimum cluster extent of ±ceil(σ) around the peak.
        let min_half_width = config.gaussian_std.ceil() as usize;
        let mut left = peak_idx.saturating_sub(min_half_width);
        let mut right = (*peak_idx + min_half_width).min(arr_len - 1);

        // Further extend if smoothed values support wider clusters
        while left > 0 && smoothed[left - 1] >= extent_threshold && !used[left - 1] {
            left -= 1;
        }
        while right + 1 < arr_len && smoothed[right + 1] >= extent_threshold && !used[right + 1] {
            right += 1;
        }

        // Don't include positions already used by another (higher) peak
        while left <= *peak_idx && used[left] {
            left += 1;
        }
        while right >= *peak_idx && used[right] {
            if right == 0 { break; }
            right -= 1;
        }

        // Skip if entirely consumed by higher peaks
        if left > right {
            stats.suppressed_by_neighbour += 1;
            continue;
        }

        // Mark positions as used
        for i in left..=right {
            used[i] = true;
        }

        // Convert back to genomic coordinates
        let genomic_start = span.0 + left as i64;
        let genomic_end = span.0 + right as i64;
        let genomic_peak = span.0 + *peak_idx as i64;

        // Sum counts in this cluster region
        let mut total_count: u64 = 0;
        let mut ds_counts: HashMap<String, u64> = HashMap::default();
        for pos in genomic_start..=genomic_end {
            if let Some((count, per_ds)) = counts.get(&pos) {
                total_count += count;
                for (ds, c) in per_ds {
                    *ds_counts.entry(ds.clone()).or_insert(0) += c;
                }
            }
        }

        if total_count > 0 {
            clusters.push((genomic_start, genomic_end, genomic_peak, total_count, ds_counts));
        }
    }
    clusters
}

/// Move each cluster's peak to its most-used base within the cluster extent.
/// Ties go to the base nearest the smoothed maximum, then to the leftmost.
fn snap_peaks_to_raw_mode(
    clusters: &mut [(i64, i64, i64, u64, HashMap<String, u64>)],
    counts: &PositionCounts,
) {
    for c in clusters.iter_mut() {
        let smoothed = c.2;
        let mut best = smoothed;
        let mut best_n = counts.get(&smoothed).map_or(0, |x| x.0);
        for pos in c.0..=c.1 {
            let n = counts.get(&pos).map_or(0, |x| x.0);
            if n > best_n || (n == best_n && n > 0 && (pos - smoothed).abs() < (best - smoothed).abs()) {
                best = pos;
                best_n = n;
            }
        }
        c.2 = best;
    }
}

/// Compute TSS cluster stats within a narrow window around the peak.
/// Returns (win_start, win_end, total_count, per_dataset_counts, peak_count, entropy, width).
fn compute_tss_window_stats(
    pos_counts: &PositionCounts,
    start: i64,
    end: i64,
    peak: i64,
    window: i64,
) -> (i64, i64, u64, HashMap<String, u64>, u64, f64, i64) {
    let win_start = if window > 0 { (peak - window).max(start) } else { start };
    let win_end = if window > 0 { (peak + window).min(end) } else { end };

    let mut total: u64 = 0;
    let mut per_ds: HashMap<String, u64> = HashMap::default();
    let mut counts: Vec<u64> = Vec::new();
    let mut min_pos: Option<i64> = None;
    let mut max_pos: Option<i64> = None;

    for pos in win_start..=win_end {
        if let Some((count, ds_counts)) = pos_counts.get(&pos) {
            if *count > 0 {
                total += *count;
                counts.push(*count);
                min_pos = Some(min_pos.map_or(pos, |v| v.min(pos)));
                max_pos = Some(max_pos.map_or(pos, |v| v.max(pos)));
                for (ds, c) in ds_counts {
                    *per_ds.entry(ds.clone()).or_insert(0) += *c;
                }
            }
        }
    }

    let peak_count = pos_counts.get(&peak).map(|(c, _)| *c).unwrap_or(0);

    // Normalized entropy (0..1); 0 = single sharp peak, 1 = uniform spread.
    let entropy = if total > 0 && counts.len() > 1 {
        let total_f = total as f64;
        let mut ent = 0.0;
        for c in &counts {
            let p = *c as f64 / total_f;
            ent -= p * p.ln();
        }
        ent / (counts.len() as f64).ln()
    } else {
        0.0
    };

    let width = if let (Some(lo), Some(hi)) = (min_pos, max_pos) {
        hi - lo
    } else {
        0
    };

    (win_start, win_end, total, per_ds, peak_count, entropy, width)
}

/// Apply 1D Gaussian smoothing with the given standard deviation and kernel
/// half-width.
fn gaussian_smooth(signal: &[f64], std: f64, half_window: usize) -> Vec<f64> {
    let n = signal.len();
    let mut result = vec![0.0; n];

    // Build Gaussian kernel
    let kernel_size = 2 * half_window + 1;
    let mut kernel = vec![0.0; kernel_size];
    let mut kernel_sum = 0.0;
    for i in 0..kernel_size {
        let x = i as f64 - half_window as f64;
        let val = (-0.5 * (x / std).powi(2)).exp();
        kernel[i] = val;
        kernel_sum += val;
    }
    // Normalize
    for k in &mut kernel {
        *k /= kernel_sum;
    }

    // Convolve
    for i in 0..n {
        let mut sum = 0.0;
        for j in 0..kernel_size {
            let src = i as i64 + j as i64 - half_window as i64;
            if src >= 0 && (src as usize) < n {
                sum += signal[src as usize] * kernel[j];
            }
        }
        result[i] = sum;
    }

    result
}

// ═══════════════════════════════════════════════════════════════════════════════
// STEP 2b: POST-HOC CLUSTER MERGING — Merge nearby peaks
// ═══════════════════════════════════════════════════════════════════════════════

/// Merge clusters whose peaks are within `merge_distance` bp of each other.
///
/// Clusters are sorted by peak position and consecutive ones merged; the
/// merged cluster takes the peak of the higher-count cluster and sums counts.
fn merge_nearby_clusters(
    clusters: &mut Vec<(i64, i64, i64, u64, HashMap<String, u64>)>,
    merge_distance: i64,
) {
    if clusters.len() < 2 || merge_distance <= 0 {
        return;
    }

    // Sort by peak position
    clusters.sort_by_key(|c| c.2);

    let mut merged: Vec<(i64, i64, i64, u64, HashMap<String, u64>)> = Vec::new();
    merged.push(clusters[0].clone());

    for i in 1..clusters.len() {
        let last = merged.last().unwrap().clone();
        let curr = &clusters[i];

        if (curr.2 - last.2).abs() <= merge_distance {
            // Merge: expand boundaries, pick peak of higher-count cluster, sum counts
            let m = merged.last_mut().unwrap();
            m.0 = m.0.min(curr.0); // start = min
            m.1 = m.1.max(curr.1); // end = max
            if curr.3 > m.3 {
                m.2 = curr.2; // peak = peak of higher-count cluster
            }
            m.3 += curr.3; // total_count = sum
            for (ds, c) in &curr.4 {
                *m.4.entry(ds.clone()).or_insert(0) += c;
            }
        } else {
            merged.push(curr.clone());
        }
    }

    *clusters = merged;
}

// ═══════════════════════════════════════════════════════════════════════════════
// STEP 2c: SATELLITE CLUSTER ABSORPTION — Remove gene-level minor satellites
// ═══════════════════════════════════════════════════════════════════════════════

/// Remove satellite clusters: those with less than `max_usage` of their gene's
/// reads lying within `max_distance` bp of a larger cluster in the same gene.
///
/// Returns the number of clusters removed.
fn absorb_satellite_clusters(
    clusters: &mut Vec<EndSiteCluster>,
    max_distance: i64,
    max_usage: f64,
) -> usize {
    if clusters.len() < 2 {
        return 0;
    }

    // Group clusters by (gene_id, chromosome, strand)
    let mut gene_groups: HashMap<(Option<String>, String, Strand), Vec<usize>> = HashMap::default();
    for (idx, c) in clusters.iter().enumerate() {
        gene_groups
            .entry((c.gene_id.clone(), c.chromosome.clone(), c.strand))
            .or_default()
            .push(idx);
    }

    let mut to_remove: Vec<usize> = Vec::new();

    for (_, indices) in &gene_groups {
        if indices.len() < 2 {
            continue;
        }

        // For each minor cluster, check if it's near a larger cluster
        for &idx in indices {
            let c = &clusters[idx];
            if c.usage >= max_usage {
                continue; // not a satellite
            }

            // Is there a larger cluster within max_distance?
            let is_satellite = indices.iter().any(|&other_idx| {
                if other_idx == idx {
                    return false;
                }
                let other = &clusters[other_idx];
                other.read_count > c.read_count
                    && (other.peak_site - c.peak_site).abs() <= max_distance
            });

            if is_satellite {
                to_remove.push(idx);
            }
        }
    }

    if to_remove.is_empty() {
        return 0;
    }

    let n = to_remove.len();
    to_remove.sort_unstable();
    to_remove.dedup();
    // Remove in reverse order to preserve indices
    for &idx in to_remove.iter().rev() {
        clusters.remove(idx);
    }
    n
}

// ═══════════════════════════════════════════════════════════════════════════════
// STEP 3: FEATURE ANNOTATION — Assign clusters to genomic features
// ═══════════════════════════════════════════════════════════════════════════════

/// Interval tree for UTR and gene feature annotation, built from
/// utr_regions.parquet and gene records.
pub struct FeatureAnnotationTree {
    /// Per-(chrom, strand) interval trees of (feature, gene_id) intervals
    trees: HashMap<(String, Strand), COITree<(GenomicFeature, String), usize>>,
}

/// Entry in the feature annotation tree.
#[allow(dead_code)]
#[derive(Debug, Clone)]
struct FeatureEntry {
    feature: GenomicFeature,
    gene_id: String,
}

impl FeatureAnnotationTree {
    /// Build from UTR regions and gene records: one interval per UTR region with
    /// its feature and one per gene body.
    pub fn build(utr_regions: &[UtrRegion], genes: &[Gene]) -> Self {
        let mut intervals: HashMap<(String, Strand), Vec<Interval<(GenomicFeature, String)>>> =
            HashMap::default();

        // Add UTR regions
        for utr in utr_regions {
            let key = (utr.chromosome.clone(), utr.strand);
            intervals
                .entry(key)
                .or_default()
                .push(Interval::new(
                    utr.start as i32,
                    utr.end as i32,
                    (utr.feature.clone(), utr.gene_id.clone()),
                ));
        }

        // Add gene body regions (lowest priority)
        for gene in genes {
            let key = (gene.chromosome.clone(), gene.strand);
            intervals
                .entry(key)
                .or_default()
                .push(Interval::new(
                    gene.start as i32,
                    gene.end as i32,
                    (GenomicFeature::Gene, format!("{}", gene.gene_id)),
                ));
        }

        // Build COITrees
        let mut trees = HashMap::default();
        for (key, ivs) in &intervals {
            trees.insert(key.clone(), COITree::new(ivs));
        }

        Self { trees }
    }

    /// Annotate a cluster position with the highest-priority overlapping feature.
    ///
    /// Returns (GenomicFeature, gene_id) for the highest-priority feature.
    pub fn annotate(
        &self,
        chrom: &str,
        position: i64,
        strand: Strand,
        cluster_type: ClusterType,
    ) -> (GenomicFeature, Option<String>) {
        let key = (chrom.to_string(), strand);
        if let Some(tree) = self.trees.get(&key) {
            let pos = position as i32;
            let mut best_feature = GenomicFeature::Intergenic;
            let mut best_gene_id = None;
            let mut best_priority = u8::MAX;

            tree.query(pos, pos + 1, |interval| {
                let (ref feature, ref gene_id) = interval.metadata;
                let priority = feature.priority(cluster_type);
                if priority < best_priority {
                    best_priority = priority;
                    best_feature = feature.clone();
                    best_gene_id = Some(gene_id.clone());
                }
            });

            (best_feature, best_gene_id)
        } else {
            (GenomicFeature::Intergenic, None)
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// STEP 4: PAS HEXAMER SCANNING — Search upstream of TES clusters
// ═══════════════════════════════════════════════════════════════════════════════

/// The poly(A) signal of a TES cluster peak: the strongest hexamer in the PAS_SCAN_WINDOW bp upstream that is
/// positioned to direct cleavage there (see [`positioned_pas`]), else a U-rich downstream element reported as
/// [`DSE_LABEL`]. Returns the motif and its distance from the hexamer's first base to the peak.
pub fn scan_pas_hexamers(
    genome: &GenomeRef,
    chrom: &str,
    peak_pos: i64,
    strand: Strand,
) -> Option<(String, i64)> {
    match genome.scan_pas_upstream(chrom, peak_pos, strand) {
        (true, Some(motif), Some(gap)) if motif != DSE_LABEL => Some((motif, gap + 6)),
        (true, Some(motif), _) => Some((motif, 0)),
        _ => None,
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// STEP 5: JOIN — Assign reads to nearest TSS/TES clusters
// ═══════════════════════════════════════════════════════════════════════════════

/// Index of clusters for nearest-cluster lookup.
pub struct ClusterIndex {
    /// Per-chrom/strand sorted list of (peak_position, cluster_id, cluster)
    index: HashMap<(String, Strand), Vec<(i64, i64, EndSiteCluster)>>,
}

impl ClusterIndex {
    /// Build from a list of clusters.
    pub fn build(clusters: &[EndSiteCluster]) -> Self {
        let mut index: HashMap<(String, Strand), Vec<(i64, i64, EndSiteCluster)>> = HashMap::default();

        for cluster in clusters {
            index
                .entry((cluster.chromosome.clone(), cluster.strand))
                .or_default()
                .push((cluster.peak_site, cluster.cluster_id, cluster.clone()));
        }

        // Sort each list by position for binary search
        for list in index.values_mut() {
            list.sort_by_key(|&(pos, _, _)| pos);
        }

        Self { index }
    }

    /// Find the nearest cluster to a given position within max_distance.
    ///
    /// Returns the nearest cluster, or None when no peak lies within max_distance.
    pub fn find_nearest(
        &self,
        chrom: &str,
        position: i64,
        strand: Strand,
        max_distance: i64,
    ) -> Option<&EndSiteCluster> {
        let key = (chrom.to_string(), strand);
        let list = self.index.get(&key)?;

        if list.is_empty() {
            return None;
        }

        // Binary search for the closest position
        let idx = list.partition_point(|&(pos, _, _)| pos < position);

        let mut best: Option<(i64, usize)> = None; // (distance, index)

        // Check idx and idx-1
        if idx < list.len() {
            let dist = (list[idx].0 - position).abs();
            if dist <= max_distance {
                best = Some((dist, idx));
            }
        }
        if idx > 0 {
            let dist = (list[idx - 1].0 - position).abs();
            if dist <= max_distance {
                if best.is_none() || dist < best.unwrap().0 {
                    best = Some((dist, idx - 1));
                }
            }
        }

        best.map(|(_, i)| &list[i].2)
    }
}

/// Join a read to TSS/TES clusters: the 5' end to the nearest TSS cluster and
/// the 3' end to the nearest TES cluster, both within max_distance. A read that
/// fails a support check is unsupported, with the check in `rejection_reason`.
pub fn join_read_to_clusters<R: ClusteringInput>(
    read: &R,
    tss_index: &ClusterIndex,
    tes_index: &ClusterIndex,
    max_distance: i64,
    hints: Option<&ReadFilterHints>,
    cage_index: Option<&CagePeakIndex>,
    config: &ClusterConfig,
) -> ReadEndSiteInfo {
    let pos_5p = match read.pos_5p() {
        Some(p) => p,
        None => return ReadEndSiteInfo {
            rejection_reason: Some("insufficient_positions".to_string()),
            ..ReadEndSiteInfo::default()
        },
    };
    let pos_3p = match read.pos_3p() {
        Some(p) => p,
        None => return ReadEndSiteInfo {
            rejection_reason: Some("insufficient_positions".to_string()),
            ..ReadEndSiteInfo::default()
        },
    };

    let tss_match = tss_index.find_nearest(read.chromosome(), pos_5p, read.strand_val(), max_distance);
    let tes_match = tes_index.find_nearest(read.chromosome(), pos_3p, read.strand_val(), max_distance);

    let tss_cluster_id = tss_match.map(|c| c.cluster_id);
    let tes_cluster_id = tes_match.map(|c| c.cluster_id);

    let polya_tail_len = read.polya_len().or(read.polya_softclip());

    // ── Support status ──────
    // A read needs a TSS or a TES cluster (OR gate).
    let cluster_ok = tss_cluster_id.is_some() || tes_cluster_id.is_some();
    let mut rejection_reason: Option<String> = None;

    // All reads, Known included, must pass the poly(A), degradation and cluster checks.

    // Degradation hints
    if let Some(h) = hints {
        let is_degraded_5p = h.is_degraded_5p.contains(read.read_idx());
        let is_degraded_3p = h.is_degraded_3p.contains(read.read_idx());
        if is_degraded_5p || is_degraded_3p {
            let reasons = if is_degraded_5p && is_degraded_3p {
                "degraded_5p_3p"
            } else if is_degraded_5p {
                "degraded_5p"
            } else {
                "degraded_3p"
            };
            rejection_reason = Some(reasons.to_string());
        }
    }

    // Poly(A) tail requirement; pt:i:0 counts as no tail
    let has_polya = polya_tail_len.map_or(false, |v| v > 0);
    if rejection_reason.is_none() && config.require_polya_tail && !has_polya {
        rejection_reason = Some("no_polya_tail".to_string());
    }

    // With poly(A) tails required, TES cluster support is also required.
    if rejection_reason.is_none() && config.require_polya_tail && tes_cluster_id.is_none() {
        rejection_reason = Some("no_tes_cluster".to_string());
    }

    // CAGE peak overlap (when a CAGE index is provided)
    if rejection_reason.is_none() {
        if let Some(cage) = cage_index {
            if !cage.overlaps_with_window(read.chromosome(), read.strand_val(), pos_5p, config.cage_window) {
                rejection_reason = Some("no_cage_overlap".to_string());
            }
        }
    }

    // Cluster support (TSS or TES)
    if rejection_reason.is_none() && !cluster_ok {
        if tss_cluster_id.is_none() && tes_cluster_id.is_none() {
            rejection_reason = Some("no_tss_tes_cluster".to_string());
        } else if tss_cluster_id.is_none() {
            rejection_reason = Some("no_tss_cluster".to_string());
        } else {
            rejection_reason = Some("no_tes_cluster".to_string());
        }
    }

    let is_supported = rejection_reason.is_none();

    ReadEndSiteInfo {
        tss_cluster_id,
        tes_cluster_id,
        tss_site: tss_match.map(|c| c.peak_site),
        polya_site: tes_match.map(|c| c.peak_site),
        tss_feature: tss_match.map(|c| c.feature.clone()),
        polya_feature: tes_match.map(|c| c.feature.clone()),
        polya_signal: tes_match.and_then(|c| c.pas_signal.clone()),
        polya_tail_len,
        is_supported,
        rejection_reason,
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// MAIN ENGINE — Orchestrates the full TSS/TES pipeline
// ═══════════════════════════════════════════════════════════════════════════════

/// TSS/TES clustering engine: accumulate_counts() per BAM, run_clustering()
/// for peak calling, then join_reads() to assign reads to clusters.
pub struct TssTesEngine {
    pub config: ClusterConfig,
    pub counts: EndSiteCounts,
    pub tss_clusters: Vec<EndSiteCluster>,
    pub tes_clusters: Vec<EndSiteCluster>,
    pub tss_index: Option<ClusterIndex>,
    pub tes_index: Option<ClusterIndex>,
    /// Single-read clusters placed on an annotated end by reference_end_support; their revalidation floor is one read.
    pub singleton_anchored: HashSet<i64>,
    cluster_counter: i64,
}

impl TssTesEngine {
    pub fn new(config: ClusterConfig) -> Self {
        Self {
            config,
            counts: EndSiteCounts::new(),
            tss_clusters: Vec::new(),
            tes_clusters: Vec::new(),
            tss_index: None,
            tes_index: None,
            singleton_anchored: HashSet::default(),
            cluster_counter: 0,
        }
    }

    /// Accumulate end-site counts from a batch of reads, with optional
    /// ReadFilterHints and CagePeakIndex.
    pub fn accumulate_counts<R: ClusteringInput>(
        &mut self,
        reads: &[R],
        hints: Option<&ReadFilterHints>,
        cage_index: Option<&CagePeakIndex>,
    ) {
        self.counts.add_reads(reads, &self.config, hints, cage_index);
    }

    /// Run Gaussian peak calling on TSS and TES positions, then annotate
    /// clusters with genomic features and, given a genome, PAS signals.
    pub fn run_clustering(
        &mut self,
        feature_tree: &FeatureAnnotationTree,
        genome: Option<&GenomeRef>,
    ) {
        use rayon::prelude::*;

        let total_reads: u64 = self.counts.dataset_read_counts.values().sum();
        let n_datasets = self.counts.dataset_read_counts.len() as u32;

        // ── TSS clustering (parallelized per chromosome/strand) ──────────
        eprintln!("[IsoForge] TSS clustering (σ={}, merge_dist={}bp, min_reads={}, min_pos_reads={})...",
                 self.config.gaussian_std, self.config.merge_distance, self.config.min_cluster_reads,
                 self.config.min_tss_position_reads);

        // Sorted so cluster IDs are assigned in genome order.
        let mut tss_keys: Vec<_> = self.counts.tss.keys().cloned().collect();
        tss_keys.sort();
        let tss_results: Vec<(Vec<EndSiteCluster>, PeakGateStats)> = tss_keys
            .par_iter()
            .map(|(chrom, strand)| {
                let pos_counts = &self.counts.tss[&(chrom.clone(), *strand)];
                let mut stats = PeakGateStats::default();
                let mut raw_clusters = call_peaks(pos_counts, &self.config, self.config.min_tss_position_reads, &mut stats);
                if self.config.tss_peak_raw_mode {
                    snap_peaks_to_raw_mode(&mut raw_clusters, pos_counts);
                }
                merge_nearby_clusters(&mut raw_clusters, self.config.merge_distance);
                let mut compact = HashSet::default();
                if self.config.compact_tss {
                    for c in compact_start_candidates(pos_counts, self.config.min_cluster_reads) {
                        if raw_clusters.iter().any(|r| c.0 <= r.1 && r.0 <= c.1) { continue; }
                        // A candidate whose first and last base lie in different annotated genes is skipped.
                        let (_, first_gene) = feature_tree.annotate(chrom, c.0, *strand, ClusterType::TSS);
                        let (_, last_gene) = feature_tree.annotate(chrom, c.1, *strand, ClusterType::TSS);
                        if first_gene != last_gene { continue; }
                        compact.insert((c.0, c.1, c.2));
                        raw_clusters.push(c);
                    }
                    raw_clusters.sort_by_key(|c| c.0);
                }

                let mut clusters = Vec::new();
                // peaks the shape gates rejected, and the reads of each gene's strongest kept TSS peak, for the rescue
                let mut rejected: Vec<(i64, i64, i64, u64, HashMap<String, u64>)> = Vec::new();
                let mut kept_max: HashMap<String, u64> = HashMap::default();
                for (start, end, peak, _count, _ds_counts) in raw_clusters {
                    let (win_start, win_end, win_count, win_ds_counts, peak_count, entropy, width) =
                        compute_tss_window_stats(pos_counts, start, end, peak, self.config.tss_peak_window);

                    if win_count < self.config.min_cluster_reads {
                        stats.below_min_cluster_reads += 1;
                        continue;
                    }

                    let peak_ratio = peak_count as f64 / win_count.max(1) as f64;
                    if std::env::var_os("ISOFORGE_TSS_GATE_TRACE").is_some() {
                        let flank = |a: i64, b: i64| -> u64 { (a..=b).filter_map(|p| pos_counts.get(&p)).map(|(c, _)| *c).sum() };
                        let (down, up) = if *strand == Strand::Minus {
                            (flank(win_start - TSS_TRACE_FLANK, win_start - 1), flank(win_end + 1, win_end + TSS_TRACE_FLANK))
                        } else {
                            (flank(win_end + 1, win_end + TSS_TRACE_FLANK), flank(win_start - TSS_TRACE_FLANK, win_start - 1))
                        };
                        let status = if compact.contains(&(start, end, peak)) { "compact".to_string() } else {
                            let mut r = Vec::new();
                            if entropy > self.config.tss_peak_max_entropy { r.push("entropy"); }
                            if peak_ratio < self.config.tss_peak_min_ratio { r.push("peak_ratio"); }
                            if width > self.config.tss_peak_max_width { r.push("width"); }
                            if r.is_empty() { "kept".to_string() } else { r.join(",") }
                        };
                        let (feature, gene) = feature_tree.annotate(chrom, peak, *strand, ClusterType::TSS);
                        TSS_GATE_TRACE.lock().unwrap().push(format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{}\t{}\t{}\t{}\t{}\t{}",
                            chrom, if *strand == Strand::Minus { "-" } else { "+" }, win_start, win_end, peak, win_count, peak_count,
                            entropy, peak_ratio, width, down, up, format!("{:?}", feature), gene.unwrap_or_default(), status));
                    }
                    if !compact.contains(&(start, end, peak)) {
                        // Each gate is counted independently; a cluster may fail more than one.
                        let bad_entropy = entropy > self.config.tss_peak_max_entropy;
                        let bad_ratio = peak_ratio < self.config.tss_peak_min_ratio;
                        let bad_width = width > self.config.tss_peak_max_width;
                        if bad_entropy { stats.rejected_entropy += 1; }
                        if bad_ratio { stats.rejected_peak_ratio += 1; }
                        if bad_width { stats.rejected_width += 1; }
                        if bad_entropy || bad_ratio || bad_width {
                            if !bad_width { rejected.push((win_start, win_end, peak, win_count, win_ds_counts)); }
                            continue;
                        }
                    }

                    let (feature, gene_id) = feature_tree.annotate(chrom, peak, *strand, ClusterType::TSS);
                    if let Some(g) = &gene_id { let m = kept_max.entry(g.clone()).or_default(); *m = (*m).max(win_count); }
                    let tpm = if total_reads > 0 {
                        (win_count as f64 / total_reads as f64) * 1_000_000.0
                    } else {
                        0.0
                    };

                    clusters.push(EndSiteCluster {
                        cluster_id: 0, // assigned below
                        cluster_type: ClusterType::TSS,
                        chromosome: chrom.clone(),
                        start: win_start,
                        end: win_end,
                        peak_site: peak,
                        strand: *strand,
                        read_count: win_count,
                        feature,
                        gene_id,
                        tpm,
                        gene_count: 0,
                        usage: 0.0,
                        frac_a: None,
                        pas_signal: None,
                        pas_distance: None,
                        n_datasets_supporting: win_ds_counts.len() as u32,
                    });
                }
                // Rescue: keeps a gene's largest rejected peak that passed the width gate, holds more reads than every
                // TSS peak the gene kept and has `tss_peak_rescue_enrichment` times the 5'-end density downstream.
                if self.config.tss_peak_rescue_min_reads > 0 {
                    rejected.sort_by(|a, b| b.3.cmp(&a.3).then(a.2.cmp(&b.2)));
                    let mut rescued_genes: HashSet<String> = HashSet::default();
                    for (win_start, win_end, peak, win_count, win_ds_counts) in rejected {
                        if win_count < self.config.tss_peak_rescue_min_reads { break; }
                        let (a, b) = if *strand == Strand::Minus { (win_start - TSS_RESCUE_FLANK, win_start - 1) } else { (win_end + 1, win_end + TSS_RESCUE_FLANK) };
                        let flank: u64 = (a..=b).filter_map(|p| pos_counts.get(&p)).map(|(c, _)| *c).sum();
                        let density = win_count as f64 / (win_end - win_start + 1).max(1) as f64;
                        if density < self.config.tss_peak_rescue_enrichment * (flank + 1) as f64 / TSS_RESCUE_FLANK as f64 { continue; }
                        let (feature, gene_id) = feature_tree.annotate(chrom, peak, *strand, ClusterType::TSS);
                        let Some(g) = gene_id.clone() else { continue };
                        if win_count <= kept_max.get(&g).copied().unwrap_or(0) || !rescued_genes.insert(g) { continue; }
                        stats.rescued += 1;
                        let tpm = if total_reads > 0 { (win_count as f64 / total_reads as f64) * 1_000_000.0 } else { 0.0 };
                        clusters.push(EndSiteCluster {
                            cluster_id: 0, cluster_type: ClusterType::TSS, chromosome: chrom.clone(), start: win_start, end: win_end,
                            peak_site: peak, strand: *strand, read_count: win_count, feature, gene_id, tpm, gene_count: 0, usage: 0.0,
                            frac_a: None, pas_signal: None, pas_distance: None, n_datasets_supporting: win_ds_counts.len() as u32,
                        });
                    }
                    clusters.sort_by_key(|c| (c.start, c.peak_site));
                }
                (clusters, stats)
            })
            .collect();

        let mut tss_gate_stats = PeakGateStats::default();
        for (_, s) in &tss_results { tss_gate_stats.merge(s); }
        if let Some(path) = std::env::var_os("ISOFORGE_TSS_GATE_TRACE") {
            let mut rows = std::mem::take(&mut *TSS_GATE_TRACE.lock().unwrap());
            rows.sort();
            let header = "chrom\tstrand\twin_start\twin_end\tpeak\twin_count\tpeak_count\tentropy\tpeak_ratio\twidth\tdown_flank\tup_flank\tfeature\tgene\tstatus\n";
            if let Err(e) = std::fs::write(&path, header.to_string() + &rows.join("\n") + "\n") {
                eprintln!("[IsoForge]   TSS gate trace not written: {e}");
            } else {
                eprintln!("[IsoForge]   TSS gate trace: {} peaks -> {:?}", rows.len(), path);
            }
        }
        if tss_gate_stats.total() > 0 {
            eprintln!(
                "[IsoForge]   TSS peaks discarded: {} suppressed by a taller neighbour, \
{} below min_cluster_reads, {} entropy, {} peak-ratio, {} width; {} of the rejected kept as their gene's strongest start",
                tss_gate_stats.suppressed_by_neighbour, tss_gate_stats.below_min_cluster_reads,
                tss_gate_stats.rejected_entropy, tss_gate_stats.rejected_peak_ratio,
                tss_gate_stats.rejected_width, tss_gate_stats.rescued);
        }

        // Flatten and assign cluster IDs
        for (clusters, _) in tss_results {
            for mut c in clusters {
                self.cluster_counter += 1;
                c.cluster_id = self.cluster_counter;
                self.tss_clusters.push(c);
            }
        }
        eprintln!("[IsoForge]   Found {} TSS clusters", self.tss_clusters.len());

        // ── TES clustering (parallelized per chromosome/strand) ──────────
        eprintln!("[IsoForge] TES/polyA clustering (σ={}, merge_dist={}bp, min_reads={})...",
                 self.config.gaussian_std, self.config.merge_distance, self.config.min_cluster_reads);

        let mut tes_keys: Vec<_> = self.counts.tes.keys().cloned().collect();
        tes_keys.sort();
        let frac_a_ref = &self.counts.frac_a_by_pos;
        let tes_results: Vec<(Vec<EndSiteCluster>, PeakGateStats)> = tes_keys
            .par_iter()
            .map(|(chrom, strand)| {
                let pos_counts = &self.counts.tes[&(chrom.clone(), *strand)];
                let mut stats = PeakGateStats::default();
                let mut raw_clusters = call_peaks(pos_counts, &self.config, 1, &mut stats);
                merge_nearby_clusters(&mut raw_clusters, self.config.merge_distance);

                let mut clusters = Vec::new();
                for (start, end, peak, count, ds_counts) in raw_clusters {
                    if count < self.config.min_cluster_reads {
                        stats.below_min_cluster_reads += 1;
                        continue;
                    }

                    let (feature, gene_id) = feature_tree.annotate(chrom, peak, *strand, ClusterType::TES);
                    let tpm = if total_reads > 0 {
                        (count as f64 / total_reads as f64) * 1_000_000.0
                    } else {
                        0.0
                    };

                    // Compute mean fraction-A for this cluster's reads
                    let frac_a = {
                        let mut all_fa: Vec<f64> = Vec::new();
                        for pos in start..=end {
                            if let Some(fas) = frac_a_ref.get(&(chrom.clone(), *strand, pos)) {
                                all_fa.extend(fas);
                            }
                        }
                        if all_fa.is_empty() { None }
                        else { Some(all_fa.iter().sum::<f64>() / all_fa.len() as f64) }
                    };

                    clusters.push(EndSiteCluster {
                        cluster_id: 0, // assigned below
                        cluster_type: ClusterType::TES,
                        chromosome: chrom.clone(),
                        start,
                        end,
                        peak_site: peak,
                        strand: *strand,
                        read_count: count,
                        feature,
                        gene_id,
                        tpm,
                        gene_count: 0,
                        usage: 0.0,
                        frac_a,
                        pas_signal: None,
                        pas_distance: None,
                        n_datasets_supporting: ds_counts.len() as u32,
                    });
                }
                (clusters, stats)
            })
            .collect();

        let mut tes_gate_stats = PeakGateStats::default();
        for (_, s) in &tes_results { tes_gate_stats.merge(s); }
        if tes_gate_stats.total() > 0 {
            eprintln!(
                "[IsoForge]   TES peaks discarded: {} suppressed by a taller neighbour, \
{} below min_cluster_reads",
                tes_gate_stats.suppressed_by_neighbour, tes_gate_stats.below_min_cluster_reads);
        }

        // Flatten, assign cluster IDs, and scan for PAS hexamers
        for (clusters, _) in tes_results {
            for mut c in clusters {
                self.cluster_counter += 1;
                c.cluster_id = self.cluster_counter;

                // PAS hexamer scanning (requires a genome reference)
                if let Some(ref genome_ref) = genome {
                    if let Some((motif, dist)) = scan_pas_hexamers(genome_ref, &c.chromosome, c.peak_site, c.strand) {
                        c.pas_signal = Some(motif);
                        c.pas_distance = Some(dist);
                    }
                }

                self.tes_clusters.push(c);
            }
        }
        eprintln!("[IsoForge]   Found {} TES clusters", self.tes_clusters.len());

        // ── Compute gene-level usage statistics ──────────────────────────
        self.compute_usage_stats();

        // ── Satellite cluster absorption ─────────────────────────────────
        // Within each gene, minor clusters (usage < satellite_max_usage) within
        // satellite_absorption_distance of a larger cluster are removed.
        {
            let max_dist = self.config.satellite_absorption_distance;
            let max_usage = self.config.satellite_max_usage;
            if max_dist > 0 && max_usage > 0.0 {
                let n_tes = absorb_satellite_clusters(
                    &mut self.tes_clusters, max_dist, max_usage,
                );
                let n_tss = absorb_satellite_clusters(
                    &mut self.tss_clusters, max_dist, max_usage,
                );
                if n_tes > 0 || n_tss > 0 {
                    eprintln!(
                        "[IsoForge]   Satellite absorption (≤{}bp, <{:.0}% usage): \
                         removed {} TES, {} TSS clusters",
                        max_dist, max_usage * 100.0, n_tes, n_tss
                    );
                    // Recompute usage stats after removal
                    self.compute_usage_stats();
                }
            }
        }

        // ── Apply replication filter ─────────────────────────────────────
        // With at least min_replication_samples datasets, a cluster survives if it
        // appears in >= min_supporting datasets or has >= min_condition_reads in total.
        if n_datasets >= self.config.min_replication_samples {
            let min_supporting = (n_datasets as f64 * self.config.min_replication_rate).ceil() as u32;
            let min_cond = self.config.min_condition_reads;
            let tss_before = self.tss_clusters.len();
            let tes_before = self.tes_clusters.len();

            self.tss_clusters.retain(|c| c.n_datasets_supporting >= min_supporting || c.read_count >= min_cond);
            // TES: replication filter, and clusters with fraction-A > 0.7 and no PAS
            // hexamer are dropped.
            self.tes_clusters.retain(|c| {
                let passes_replication = c.n_datasets_supporting >= min_supporting || c.read_count >= min_cond;
                let is_internal_priming = c.frac_a.map_or(false, |fa| fa > 0.7) && c.pas_signal.is_none();
                passes_replication && !is_internal_priming
            });

            eprintln!(
                "[IsoForge]   Replication filter (min {} datasets OR {} reads): TSS {} → {}, TES {} → {}",
                min_supporting,
                min_cond,
                tss_before,
                self.tss_clusters.len(),
                tes_before,
                self.tes_clusters.len()
            );
        }

        // ── Build search indices ─────────────────────────────────────────
        self.tss_index = Some(ClusterIndex::build(&self.tss_clusters));
        self.tes_index = Some(ClusterIndex::build(&self.tes_clusters));

        eprintln!(
            "[IsoForge] Clustering complete: {} TSS, {} TES clusters",
            self.tss_clusters.len(),
            self.tes_clusters.len()
        );
    }

    /// Compute per-gene usage statistics for each cluster.
    pub(crate) fn compute_usage_stats(&mut self) {
        // Accumulate total reads per gene for each cluster type
        let mut tss_gene_totals: HashMap<String, u64> = HashMap::default();
        let mut tes_gene_totals: HashMap<String, u64> = HashMap::default();

        for c in &self.tss_clusters {
            if let Some(ref gid) = c.gene_id {
                *tss_gene_totals.entry(gid.clone()).or_insert(0) += c.read_count;
            }
        }
        for c in &self.tes_clusters {
            if let Some(ref gid) = c.gene_id {
                *tes_gene_totals.entry(gid.clone()).or_insert(0) += c.read_count;
            }
        }

        // Set gene_count and usage for each cluster
        for c in &mut self.tss_clusters {
            if let Some(ref gid) = c.gene_id {
                let gene_total = *tss_gene_totals.get(gid).unwrap_or(&0);
                c.gene_count = gene_total;
                c.usage = if gene_total > 0 {
                    c.read_count as f64 / gene_total as f64
                } else {
                    0.0
                };
            }
        }
        for c in &mut self.tes_clusters {
            if let Some(ref gid) = c.gene_id {
                let gene_total = *tes_gene_totals.get(gid).unwrap_or(&0);
                c.gene_count = gene_total;
                c.usage = if gene_total > 0 {
                    c.read_count as f64 / gene_total as f64
                } else {
                    0.0
                };
            }
        }
    }

    /// Join all reads to their nearest TSS/TES clusters, applying the poly(A),
    /// CAGE and degradation gates. Returns a vector indexed by `ReadIdx`, sized
    /// `n_reads`; reads outside clustering keep the default (unsupported) entry.
    pub fn join_reads<R: ClusteringInput>(
        &self,
        all_reads: &[R],
        n_reads: usize,
        hints: Option<&ReadFilterHints>,
        cage_index: Option<&CagePeakIndex>,
    ) -> Vec<ReadEndSiteInfo> {
        let tss_idx = self.tss_index.as_ref().expect("Must run clustering first");
        let tes_idx = self.tes_index.as_ref().expect("Must run clustering first");
        let max_dist = self.config.tss_tes_distance;
        let config = &self.config;

        // One record per read, in the order of `all_reads`.
        let results: Vec<ReadEndSiteInfo> = all_reads
            .par_iter()
            .map(|read| join_read_to_clusters(read, tss_idx, tes_idx, max_dist, hints, cage_index, config))
            .collect();

        // Collect rejection stats
        let mut n_rejected_polya = 0u64;
        let mut n_rejected_cage = 0u64;
        let mut n_rejected_degraded = 0u64;
        let mut n_rejected_no_cluster = 0u64;
        let mut n_supported = 0u64;

        for info in &results {
            if info.is_supported {
                n_supported += 1;
            } else if let Some(ref reason) = info.rejection_reason {
                match reason.as_str() {
                    "no_polya_tail" => n_rejected_polya += 1,
                    "no_cage_overlap" => n_rejected_cage += 1,
                    r if r.starts_with("degraded") => n_rejected_degraded += 1,
                    _ => n_rejected_no_cluster += 1,
                }
            }
        }

        let total = results.len() as u64;
        eprintln!(
            "[IsoForge] Read support summary: {}/{} supported ({:.1}%)",
            n_supported, total, 100.0 * n_supported as f64 / total.max(1) as f64
        );
        eprintln!(
            "[IsoForge]   Rejected: {} no-polyA, {} no-CAGE, {} degraded, {} no-cluster",
            n_rejected_polya, n_rejected_cage, n_rejected_degraded, n_rejected_no_cluster
        );

        // Records are indexed by read index: when the reads already come in that order the records
        // are returned as they are, otherwise each is placed at its read's index.
        if results.len() == n_reads && all_reads.iter().enumerate().all(|(i, r)| r.read_idx() as usize == i) {
            return results;
        }
        let mut dense = vec![ReadEndSiteInfo::default(); n_reads];
        for (read, info) in all_reads.iter().zip(results) {
            let idx = read.read_idx() as usize;
            if idx < dense.len() {
                dense[idx] = info;
            }
        }
        dense
    }

    /// Post-filter cluster revalidation.
    ///
    /// Counts supported reads per TSS/TES cluster; clusters below
    /// min_cluster_reads have support revoked from their reads, repeated until
    /// no further reads are revoked.
    pub fn revalidate_clusters(
        &self,
        join_results: &mut [ReadEndSiteInfo],
        hints: Option<&ReadFilterHints>,
    ) {
        let min_reads = self.config.min_cluster_reads;
        let min_reads_known = self.config.min_cluster_reads_known;
        let mut total_revoked = 0u64;
        let mut total_zombie_tss = 0usize;
        let mut total_zombie_tes = 0usize;
        let mut iteration = 0u32;

        loop {
            iteration += 1;
            let mut revoked_this_round = 0u64;

            // Count supported reads per cluster
            let mut tss_counts: HashMap<i64, u64> = HashMap::default();
            let mut tes_counts: HashMap<i64, u64> = HashMap::default();

            for info in join_results.iter() {
                if info.is_supported {
                    if let Some(cid) = info.tss_cluster_id {
                        *tss_counts.entry(cid).or_insert(0) += 1;
                    }
                    if let Some(cid) = info.tes_cluster_id {
                        *tes_counts.entry(cid).or_insert(0) += 1;
                    }
                }
            }

            // Hard zombie: below min_reads_known (1 for clusters in `singleton_anchored`); all reads revoked.
            // Weak zombie: between min_reads_known and min_reads, only non-Known reads revoked.
            let hard_floor = |cid: i64| -> u64 {
                if self.singleton_anchored.contains(&cid) { 1 } else { min_reads_known }
            };
            let zombie_tss: HashSet<i64> = self.tss_clusters.iter()
                .filter(|c| tss_counts.get(&c.cluster_id).copied().unwrap_or(0) < hard_floor(c.cluster_id))
                .map(|c| c.cluster_id)
                .collect();
            let zombie_tes: HashSet<i64> = self.tes_clusters.iter()
                .filter(|c| tes_counts.get(&c.cluster_id).copied().unwrap_or(0) < hard_floor(c.cluster_id))
                .map(|c| c.cluster_id)
                .collect();
            // Weak zombies: clusters between min_reads_known and min_reads
            let weak_zombie_tss: HashSet<i64> = self.tss_clusters.iter()
                .filter(|c| {
                    let n = tss_counts.get(&c.cluster_id).copied().unwrap_or(0);
                    n >= min_reads_known && n < min_reads
                })
                .map(|c| c.cluster_id)
                .collect();
            let weak_zombie_tes: HashSet<i64> = self.tes_clusters.iter()
                .filter(|c| {
                    let n = tes_counts.get(&c.cluster_id).copied().unwrap_or(0);
                    n >= min_reads_known && n < min_reads
                })
                .map(|c| c.cluster_id)
                .collect();

            if iteration == 1 {
                total_zombie_tss = zombie_tss.len() + weak_zombie_tss.len();
                total_zombie_tes = zombie_tes.len() + weak_zombie_tes.len();
            }

            // Revoke support for reads joined to zombie/weak-zombie clusters
            for (read_idx, info) in join_results.iter_mut().enumerate() {
                if !info.is_supported { continue; }

                // Known/FSM reads are revoked only by hard-zombie clusters.
                let is_known = hints.map_or(false, |h| h.is_known_read.contains(read_idx as ReadIdx));

                let tss_hard = info.tss_cluster_id.map_or(false, |cid| zombie_tss.contains(&cid));
                let tes_hard = info.tes_cluster_id.map_or(false, |cid| zombie_tes.contains(&cid));
                let tss_weak = !is_known && info.tss_cluster_id.map_or(false, |cid| weak_zombie_tss.contains(&cid));
                let tes_weak = !is_known && info.tes_cluster_id.map_or(false, |cid| weak_zombie_tes.contains(&cid));

                if tss_hard || tes_hard || tss_weak || tes_weak {
                    info.is_supported = false;
                    info.rejection_reason = Some(
                        if tss_hard && tes_hard {
                            "zombie_tss_tes_cluster".to_string()
                        } else if tss_hard {
                            "zombie_tss_cluster".to_string()
                        } else if tes_hard {
                            "zombie_tes_cluster".to_string()
                        } else if tss_weak {
                            "weak_zombie_tss_cluster".to_string()
                        } else {
                            "weak_zombie_tes_cluster".to_string()
                        }
                    );
                    revoked_this_round += 1;
                }
            }

            total_revoked += revoked_this_round;

            // Converged: no new zombie reads this round
            if revoked_this_round == 0 {
                break;
            }

            // Iteration limit
            if iteration > 10 {
                eprintln!("[IsoForge] WARNING: Cluster revalidation did not converge after {} iterations", iteration);
                break;
            }
        }

        if total_revoked > 0 || total_zombie_tss > 0 || total_zombie_tes > 0 {
            eprintln!(
                "[IsoForge] Cluster revalidation: {} zombie TSS clusters, {} zombie TES clusters, \
                 {} reads revoked in {} iterations (min {}/{} reads: novel/Known)",
                total_zombie_tss, total_zombie_tes,
                total_revoked, iteration, min_reads, min_reads_known
            );
        }
    }

    /// Get summary statistics for logging.
    pub fn summary(&self) -> ClusteringSummary {
        let tss_supported = self.tss_clusters.iter().map(|c| c.read_count).sum();
        let tes_supported = self.tes_clusters.iter().map(|c| c.read_count).sum();
        let tes_with_pas = self.tes_clusters.iter().filter(|c| c.pas_signal.is_some()).count();

        ClusteringSummary {
            n_tss_clusters: self.tss_clusters.len(),
            n_tes_clusters: self.tes_clusters.len(),
            total_tss_reads: tss_supported,
            total_tes_reads: tes_supported,
            tes_with_pas_signal: tes_with_pas,
            tes_without_pas_signal: self.tes_clusters.len() - tes_with_pas,
        }
    }
}

/// Summary statistics from clustering.
#[derive(Debug, Clone)]
pub struct ClusteringSummary {
    pub n_tss_clusters: usize,
    pub n_tes_clusters: usize,
    pub total_tss_reads: u64,
    pub total_tes_reads: u64,
    pub tes_with_pas_signal: usize,
    pub tes_without_pas_signal: usize,
}

// ═══════════════════════════════════════════════════════════════════════════════
// PARQUET I/O — Write cluster files
// ═══════════════════════════════════════════════════════════════════════════════

use arrow::array::*;
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::properties::WriterProperties;

fn parquet_props() -> WriterProperties {
    WriterProperties::builder()
        .set_compression(Compression::ZSTD(Default::default()))
        .build()
}

/// Write TSS or TES clusters to Parquet.
pub fn write_clusters(
    path: &Path,
    clusters: &[EndSiteCluster],
) -> IsoForgeResult<()> {
    let mut cluster_ids = Int64Builder::new();
    let mut cluster_types = StringBuilder::new();
    let mut chroms = StringBuilder::new();
    let mut starts = Int64Builder::new();
    let mut ends = Int64Builder::new();
    let mut peaks = Int64Builder::new();
    let mut strands = StringBuilder::new();
    let mut counts = UInt64Builder::new();
    let mut features = StringBuilder::new();
    let mut gene_ids = StringBuilder::new();
    let mut tpms = Float64Builder::new();
    let mut gene_counts = UInt64Builder::new();
    let mut usages = Float64Builder::new();
    let mut frac_as = Float64Builder::new();
    let mut pas_signals = StringBuilder::new();
    let mut pas_distances = Int64Builder::new();
    let mut n_datasets = UInt32Builder::new();

    for c in clusters {
        cluster_ids.append_value(c.cluster_id);
        cluster_types.append_value(c.cluster_type.to_string());
        chroms.append_value(&c.chromosome);
        starts.append_value(c.start);
        ends.append_value(c.end);
        peaks.append_value(c.peak_site);
        strands.append_value(c.strand.as_str());
        counts.append_value(c.read_count);
        features.append_value(c.feature.as_str());
        gene_ids.append_value(c.gene_id.as_deref().unwrap_or(""));
        tpms.append_value(c.tpm);
        gene_counts.append_value(c.gene_count);
        usages.append_value(c.usage);
        match c.frac_a {
            Some(fa) => frac_as.append_value(fa),
            None => frac_as.append_null(),
        }
        match &c.pas_signal {
            Some(s) => pas_signals.append_value(s),
            None => pas_signals.append_null(),
        }
        match c.pas_distance {
            Some(d) => pas_distances.append_value(d),
            None => pas_distances.append_null(),
        }
        n_datasets.append_value(c.n_datasets_supporting);
    }

    let schema = Arc::new(Schema::new(vec![
        Field::new("cluster_id", DataType::Int64, false),
        Field::new("cluster_type", DataType::Utf8, false),
        Field::new("chromosome", DataType::Utf8, false),
        Field::new("start", DataType::Int64, false),
        Field::new("end", DataType::Int64, false),
        Field::new("peak_site", DataType::Int64, false),
        Field::new("strand", DataType::Utf8, false),
        Field::new("read_count", DataType::UInt64, false),
        Field::new("feature", DataType::Utf8, false),
        Field::new("gene_id", DataType::Utf8, true),
        Field::new("tpm", DataType::Float64, false),
        Field::new("gene_count", DataType::UInt64, false),
        Field::new("usage", DataType::Float64, false),
        Field::new("frac_a", DataType::Float64, true),
        Field::new("pas_signal", DataType::Utf8, true),
        Field::new("pas_distance", DataType::Int64, true),
        Field::new("n_datasets_supporting", DataType::UInt32, false),
    ]));

    let batch = RecordBatch::try_new(
        schema.clone(),
        vec![
            Arc::new(cluster_ids.finish()),
            Arc::new(cluster_types.finish()),
            Arc::new(chroms.finish()),
            Arc::new(starts.finish()),
            Arc::new(ends.finish()),
            Arc::new(peaks.finish()),
            Arc::new(strands.finish()),
            Arc::new(counts.finish()),
            Arc::new(features.finish()),
            Arc::new(gene_ids.finish()),
            Arc::new(tpms.finish()),
            Arc::new(gene_counts.finish()),
            Arc::new(usages.finish()),
            Arc::new(frac_as.finish()),
            Arc::new(pas_signals.finish()),
            Arc::new(pas_distances.finish()),
            Arc::new(n_datasets.finish()),
        ],
    )?;

    let file = std::fs::File::create(path)?;
    let mut writer = ArrowWriter::try_new(file, schema, Some(parquet_props()))?;
    writer.write(&batch)?;
    writer.close()?;
    Ok(())
}

/// Write a BED cluster file.
///
/// 14-column BED for TES (polyA): Chrom Start End Site Count Strand Feature Gene TPM GeneCount Usage FracA Signal AnnotatedSite
/// 12-column BED for TSS:         Chrom Start End Site Count Strand Feature Gene TPM GeneCount Usage AnnotatedSite
pub fn write_clusters_bed(
    path: &Path,
    clusters: &[EndSiteCluster],
) -> IsoForgeResult<()> {
    use std::io::Write;
    let mut file = std::io::BufWriter::new(std::fs::File::create(path)?);

    for c in clusters {
        let gene_id = c.gene_id.as_deref().unwrap_or(".");

        match c.cluster_type {
            ClusterType::TES => {
                let frac_a_str = c.frac_a.map(|f| format!("{:.4}", f)).unwrap_or_else(|| ".".to_string());
                let signal = c.pas_signal.as_deref().unwrap_or(".");
                writeln!(
                    file,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{}\t{:.4}\t{}\t{}\t.",
                    c.chromosome,
                    c.start - 1, // BED is 0-based
                    c.end,
                    c.peak_site,
                    c.read_count,
                    c.strand.as_str(),
                    c.feature.as_str(),
                    gene_id,
                    c.tpm,
                    c.gene_count,
                    c.usage,
                    frac_a_str,
                    signal,
                ).map_err(|e| IsoForgeError::IoError(e))?;
            }
            ClusterType::TSS => {
                writeln!(
                    file,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{}\t{:.4}\t.",
                    c.chromosome,
                    c.start - 1,
                    c.end,
                    c.peak_site,
                    c.read_count,
                    c.strand.as_str(),
                    c.feature.as_str(),
                    gene_id,
                    c.tpm,
                    c.gene_count,
                    c.usage,
                ).map_err(|e| IsoForgeError::IoError(e))?;
            }
        }
    }

    Ok(())
}

#[derive(Debug, Clone, Default)]
pub struct GroupClusterStat {
    pub read_count: u64,
    pub tpm: f64,
    pub n_datasets_supporting: u32,
}

pub type GroupClusterStats = HashMap<(String, i64), GroupClusterStat>;

pub fn infer_treatment_groups(datasets: &[String]) -> HashMap<String, Vec<String>> {
    let mut groups: HashMap<String, Vec<String>> = HashMap::default();
    for ds in datasets {
        let treatment = ds
            .split_once('_')
            .map(|(prefix, _)| prefix)
            .unwrap_or(ds.as_str())
            .to_string();
        groups.entry(treatment).or_default().push(ds.clone());
    }
    groups
}

pub fn compute_per_group_cluster_stats(
    tss_clusters: &[EndSiteCluster],
    tes_clusters: &[EndSiteCluster],
    group_tss_counts: &HashMap<(String, i64), u64>,
    group_tes_counts: &HashMap<(String, i64), u64>,
    group_total_reads: &HashMap<String, u64>,
    group_dataset_support: &HashMap<(String, i64), HashSet<String>>,
) -> (GroupClusterStats, GroupClusterStats) {
    let mut tss_stats: GroupClusterStats = HashMap::default();
    let mut tes_stats: GroupClusterStats = HashMap::default();

    for cluster in tss_clusters {
        for ((group_name, cluster_id), &count) in group_tss_counts {
            if *cluster_id != cluster.cluster_id {
                continue;
            }
            let total = group_total_reads.get(group_name).copied().unwrap_or(0);
            let tpm = if total > 0 {
                (count as f64 / total as f64) * 1_000_000.0
            } else {
                0.0
            };
            let n_datasets_supporting = group_dataset_support
                .get(&(group_name.clone(), *cluster_id))
                .map(|set| set.len() as u32)
                .unwrap_or(0);
            tss_stats.insert(
                (group_name.clone(), *cluster_id),
                GroupClusterStat {
                    read_count: count,
                    tpm,
                    n_datasets_supporting,
                },
            );
        }
    }

    for cluster in tes_clusters {
        for ((group_name, cluster_id), &count) in group_tes_counts {
            if *cluster_id != cluster.cluster_id {
                continue;
            }
            let total = group_total_reads.get(group_name).copied().unwrap_or(0);
            let tpm = if total > 0 {
                (count as f64 / total as f64) * 1_000_000.0
            } else {
                0.0
            };
            let n_datasets_supporting = group_dataset_support
                .get(&(group_name.clone(), *cluster_id))
                .map(|set| set.len() as u32)
                .unwrap_or(0);
            tes_stats.insert(
                (group_name.clone(), *cluster_id),
                GroupClusterStat {
                    read_count: count,
                    tpm,
                    n_datasets_supporting,
                },
            );
        }
    }

    (tss_stats, tes_stats)
}

pub fn write_group_tes_bed(
    path: &Path,
    clusters: &[EndSiteCluster],
    stats: &GroupClusterStats,
    group_name: &str,
) -> IsoForgeResult<()> {
    use std::io::Write;
    let mut file = std::io::BufWriter::new(std::fs::File::create(path)?);

    for c in clusters.iter().filter(|c| c.cluster_type == ClusterType::TES) {
        let Some(s) = stats.get(&(group_name.to_string(), c.cluster_id)) else {
            continue;
        };
        if s.read_count == 0 {
            continue;
        }

        let gene_id = c.gene_id.as_deref().unwrap_or(".");
        let frac_a_str = c.frac_a.map(|f| format!("{:.4}", f)).unwrap_or_else(|| ".".to_string());
        let signal = c.pas_signal.as_deref().unwrap_or(".");
        writeln!(
            file,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{}\t{:.4}\t{}\t{}\t.",
            c.chromosome,
            c.start - 1,
            c.end,
            c.peak_site,
            s.read_count,
            c.strand.as_str(),
            c.feature.as_str(),
            gene_id,
            s.tpm,
            c.gene_count,
            c.usage,
            frac_a_str,
            signal,
        )?;
    }

    Ok(())
}

pub fn write_group_tss_bed(
    path: &Path,
    clusters: &[EndSiteCluster],
    stats: &GroupClusterStats,
    group_name: &str,
) -> IsoForgeResult<()> {
    use std::io::Write;
    let mut file = std::io::BufWriter::new(std::fs::File::create(path)?);

    for c in clusters.iter().filter(|c| c.cluster_type == ClusterType::TSS) {
        let Some(s) = stats.get(&(group_name.to_string(), c.cluster_id)) else {
            continue;
        };
        if s.read_count == 0 {
            continue;
        }

        let gene_id = c.gene_id.as_deref().unwrap_or(".");
        writeln!(
            file,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{}\t{:.4}\t.",
            c.chromosome,
            c.start - 1,
            c.end,
            c.peak_site,
            s.read_count,
            c.strand.as_str(),
            c.feature.as_str(),
            gene_id,
            s.tpm,
            c.gene_count,
            c.usage,
        )?;
    }

    Ok(())
}

pub fn load_chrom_sizes(path: &Path) -> IsoForgeResult<Vec<(String, i64)>> {
    use std::io::{BufRead, BufReader};
    let f = std::fs::File::open(path)?;
    let r = BufReader::new(f);
    let mut out = Vec::new();

    for line in r.lines() {
        let line = line?;
        let mut parts = line.split('\t');
        let Some(chrom) = parts.next() else { continue; };
        let Some(size_str) = parts.next() else { continue; };
        if let Ok(size) = size_str.parse::<i64>() {
            out.push((chrom.to_string(), size));
        }
    }

    Ok(out)
}

pub fn write_bedgraph(
    path: &Path,
    counts: &HashMap<(String, i64), u64>,
    _chrom_sizes: &[(String, i64)],
) -> IsoForgeResult<()> {
    use std::io::Write;
    let mut rows: Vec<(&String, i64, u64)> = counts
        .iter()
        .map(|((chrom, pos), &count)| (chrom, *pos, count))
        .collect();
    rows.sort_by(|a, b| a.0.cmp(b.0).then_with(|| a.1.cmp(&b.1)));

    let mut file = std::io::BufWriter::new(std::fs::File::create(path)?);
    for (chrom, pos, count) in rows {
        if pos <= 0 {
            continue;
        }
        let start0 = pos - 1;
        let end = pos;
        writeln!(file, "{}\t{}\t{}\t{}", chrom, start0, end, count)?;
    }

    Ok(())
}

pub fn bedgraph_to_bigwig(
    bedgraph_path: &Path,
    chrom_sizes_path: &Path,
    bw_path: &Path,
) -> IsoForgeResult<bool> {
    use std::process::Command;

    let status = Command::new("bedGraphToBigWig")
        .arg(bedgraph_path)
        .arg(chrom_sizes_path)
        .arg(bw_path)
        .status();

    match status {
        Ok(s) => Ok(s.success()),
        Err(_) => Ok(false),
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// UTR REGION LOADING — Read utr_regions.parquet for feature annotation
// ═══════════════════════════════════════════════════════════════════════════════

/// Load UTR regions from the Parquet file written by gtf_to_parquet.
pub fn load_utr_regions(path: &Path) -> IsoForgeResult<Vec<UtrRegion>> {
    if !path.exists() {
        eprintln!("[IsoForge] Warning: utr_regions.parquet not found at {:?}, feature annotation will use gene-level only", path);
        return Ok(Vec::new());
    }

    let table = parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder::try_new(
        std::fs::File::open(path)?,
    )?
    .build()?;

    let mut regions = Vec::new();
    for batch_result in table {
        let batch = batch_result?;
        let chroms = batch.column(0).as_any().downcast_ref::<StringArray>().unwrap();
        let starts = batch.column(1).as_any().downcast_ref::<Int64Array>().unwrap();
        let ends = batch.column(2).as_any().downcast_ref::<Int64Array>().unwrap();
        let strands = batch.column(3).as_any().downcast_ref::<StringArray>().unwrap();
        let features = batch.column(4).as_any().downcast_ref::<StringArray>().unwrap();
        let gene_ids_col = batch.column(5).as_any().downcast_ref::<StringArray>().unwrap();
        let transcript_ids_col = batch.column(6).as_any().downcast_ref::<StringArray>().unwrap();

        for i in 0..batch.num_rows() {
            let strand = match strands.value(i) {
                "+" => Strand::Plus,
                "-" => Strand::Minus,
                _ => continue,
            };

            regions.push(UtrRegion {
                chromosome: chroms.value(i).to_string(),
                start: starts.value(i),
                end: ends.value(i),
                strand,
                feature: GenomicFeature::from_str(features.value(i)),
                gene_id: gene_ids_col.value(i).to_string(),
                transcript_id: transcript_ids_col.value(i).to_string(),
            });
        }
    }

    eprintln!("[IsoForge] Loaded {} UTR regions for feature annotation", regions.len());
    Ok(regions)
}

// ═══════════════════════════════════════════════════════════════════════════════
// TESTS
// ═══════════════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn truncation_excludes_only_the_missing_end_on_both_strands() {
        use crate::annotation_engine::ObsCore;
        use std::sync::Arc;
        for strand in [Strand::Plus, Strand::Minus] {
            for missing_five in [false, true] {
                let read = ObsCore {
                    gene_id: 1, transcript_id: 1, start_pos: if strand == Strand::Plus { 100 } else { 900 },
                    end_pos: if strand == Strand::Plus { 900 } else { 100 },
                    dataset: Arc::from("sample"), location: Arc::from("chr1"), fraction_as: Some(0.0),
                    polya_tail_len: Some(30), polya_tail_softclip: None, strand, n_exons: 2,
                    read_idx: 0, hint_flags: 0, is_internal_priming: false, max_splice_shift: 0,
                    cap_signature: false,
                    unaligned_3p: false,
                };
                let mut hints = ReadFilterHints::default();
                hints.mark_truncation(0, crate::types::HINT_IS_ISM | if missing_five {
                    crate::types::HINT_IS_ISM_SUFFIX
                } else { crate::types::HINT_IS_ISM_PREFIX });
                let mut counts = EndSiteCounts::default();
                counts.add_reads(&[read], &ClusterConfig::default(), Some(&hints), None);
                assert_eq!(counts.tss.is_empty(), missing_five);
                assert_eq!(counts.tes.is_empty(), !missing_five);
            }
        }
    }

    #[test]
    fn test_gaussian_smooth() {
        let signal = vec![0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0];
        let smoothed = gaussian_smooth(&signal, 1.0, 3);
        // Peak should be at index 3, smoothed values around it
        assert!(smoothed[3] > smoothed[2]);
        assert!(smoothed[3] > smoothed[4]);
        // Smoothed peak should be lower than original
        assert!(smoothed[3] < 10.0);
        // Values should be symmetric around peak
        assert!((smoothed[2] - smoothed[4]).abs() < 1e-10);
    }

    #[test]
    fn test_call_peaks_simple() {
        let config = ClusterConfig::default();
        let mut counts: PositionCounts = HashMap::default();

        // Create two clear peaks at positions 100 and 200
        for offset in -5i64..=5 {
            let count = (10 - offset.abs()) as u64;
            let mut ds = HashMap::default();
            ds.insert("ds1".to_string(), count);
            counts.insert(100 + offset as i64, (count, ds.clone()));
            counts.insert(200 + offset as i64, (count, ds));
        }

        let peaks = call_peaks(&counts, &config, 1, &mut PeakGateStats::default());
        assert!(peaks.len() >= 2, "Should find at least 2 peaks, got {}", peaks.len());
    }

    #[test]
    fn test_pas_hexamers_count() {
        assert_eq!(PAS_HEXAMERS.len(), 24, "Should have 24 curated PAS hexamers");
    }

    #[test]
    fn test_cluster_index_nearest() {
        let clusters = vec![
            EndSiteCluster {
                cluster_id: 1,
                cluster_type: ClusterType::TES,
                chromosome: "chr1".to_string(),
                start: 95,
                end: 105,
                peak_site: 100,
                strand: Strand::Plus,
                read_count: 50,
                feature: GenomicFeature::ThreePrimeUtr,
                gene_id: Some("GENE1".to_string()),
                tpm: 10.0,
                gene_count: 50,
                usage: 1.0,
                frac_a: Some(0.3),
                pas_signal: Some("AATAAA".to_string()),
                pas_distance: Some(20),
                n_datasets_supporting: 3,
            },
        ];

        let index = ClusterIndex::build(&clusters);

        // Within distance
        let found = index.find_nearest("chr1", 110, Strand::Plus, 50);
        assert!(found.is_some());
        assert_eq!(found.unwrap().cluster_id, 1);

        // Too far
        let not_found = index.find_nearest("chr1", 200, Strand::Plus, 50);
        assert!(not_found.is_none());

        // Wrong strand
        let wrong_strand = index.find_nearest("chr1", 100, Strand::Minus, 50);
        assert!(wrong_strand.is_none());
    }
}

#[cfg(test)]
mod compact_start_tests {
    use super::*;
    #[test]
    fn a_genes_strongest_start_spread_over_a_few_bases_is_rescued() {
        let region = |start, end, gene: &str| UtrRegion { chromosome: "chr1".into(), start, end, strand: Strand::Plus,
            feature: GenomicFeature::FivePrimeUtr, gene_id: gene.into(), transcript_id: format!("T{gene}") };
        let tree = FeatureAnnotationTree::build(&[region(1000, 1100, "G1"), region(5000, 5600, "G3"), region(9000, 9600, "G2")], &[]);
        let mut e = TssTesEngine::new(ClusterConfig::default());
        {
            let m = e.counts.tss.entry(("chr1".to_string(), Strand::Plus)).or_default();
            let mut put = |p: i64, n: u64| { m.insert(p, (n, [("d".to_string(), n)].into_iter().collect())); };
            // G1: its only start, 6 reads at each of 21 bases: fails the entropy and peak-ratio gates
            for p in 1010..=1030 { put(p, 6); }
            // G3: a sharp kept peak and a spread one with more reads 400 bp downstream, which is kept too
            put(5050, 20);
            for p in 5450..=5470 { put(p, 12); }
            // G2: a sharp kept peak and a spread one with fewer reads 400 bp downstream, which stays rejected
            put(9050, 200);
            for p in 9450..=9470 { put(p, 6); }
        }
        e.counts.dataset_read_counts.insert("d".into(), 1000);
        e.run_clustering(&tree, None);
        let peaks: Vec<i64> = e.tss_clusters.iter().map(|c| c.peak_site).collect();
        assert!(peaks.iter().any(|p| (1010..=1030).contains(p)), "{peaks:?}");
        assert!(peaks.contains(&5050) && peaks.iter().any(|p| (5450..=5470).contains(p)), "{peaks:?}");
        assert!(peaks.contains(&9050) && !peaks.iter().any(|p| (9450..=9470).contains(p)), "{peaks:?}");
        // disabled: the spread start is rejected
        let mut off = TssTesEngine::new(ClusterConfig { tss_peak_rescue_min_reads: 0, ..ClusterConfig::default() });
        off.counts = std::mem::take(&mut e.counts);
        off.run_clustering(&tree, None);
        assert!(!off.tss_clusters.iter().any(|c| (1010..=1030).contains(&c.peak_site)));
    }

    #[test]
    fn raw_mode_peak_keeps_a_dispersed_promoter() {
        // Read starts per base where smoothing places the maximum on a minor base.
        let mut pos = Vec::new();
        for (p, n) in [(66, 6), (68, 211), (69, 29), (70, 120), (71, 89), (72, 152),
                       (73, 38), (74, 3), (75, 9), (76, 72), (77, 17), (78, 2)] {
            pos.extend(std::iter::repeat(1000 + p).take(n));
        }
        let c = counts(&pos);
        let cfg = ClusterConfig::default();
        let mut stats = PeakGateStats::default();
        let mut raw = call_peaks(&c, &cfg, cfg.min_tss_position_reads, &mut stats);
        assert_eq!(raw.len(), 1, "{:?}", raw.iter().map(|r| (r.0, r.1, r.2)).collect::<Vec<_>>());
        let ratio = |r: &(i64, i64, i64, u64, HashMap<String, u64>)| {
            let (_, _, total, _, at_peak, _, _) = compute_tss_window_stats(&c, r.0, r.1, r.2, cfg.tss_peak_window);
            at_peak as f64 / total as f64
        };
        let smoothed = raw[0].2;
        assert_ne!(smoothed, 1068, "the fixture must reproduce the smoothing shift");
        assert!(ratio(&raw[0]) < cfg.tss_peak_min_ratio, "smoothed peak should fail the ratio gate");
        snap_peaks_to_raw_mode(&mut raw, &c);
        assert_eq!(raw[0].2, 1068);
        assert!(ratio(&raw[0]) >= cfg.tss_peak_min_ratio);
    }

    #[test]
    fn raw_mode_leaves_a_sharp_peak_unchanged() {
        let mut pos = Vec::new();
        for (p, n) in [(34, 3), (35, 16), (36, 4), (40, 5)] {
            pos.extend(std::iter::repeat(1000 + p).take(n));
        }
        let c = counts(&pos);
        let cfg = ClusterConfig::default();
        let mut stats = PeakGateStats::default();
        let mut raw = call_peaks(&c, &cfg, cfg.min_tss_position_reads, &mut stats);
        let before: Vec<i64> = raw.iter().map(|r| r.2).collect();
        snap_peaks_to_raw_mode(&mut raw, &c);
        assert_eq!(before, raw.iter().map(|r| r.2).collect::<Vec<_>>());
    }

    /// Peaks over one array spanning every counted position.
    fn dense_peaks(c: &PositionCounts, cfg: &ClusterConfig, min_reads: u64, stats: &mut PeakGateStats)
        -> Vec<(i64, i64, i64, u64, HashMap<String, u64>)> {
        let pad = cfg.cluster_window as i64;
        let span = (*c.keys().min().unwrap() - pad, *c.keys().max().unwrap() + pad);
        let mut sites: Vec<(i64, u64)> = c.iter().filter(|(_, v)| v.0 >= min_reads).map(|(&p, v)| (p, v.0)).collect();
        sites.sort_unstable();
        let mut out = peaks_on_span(c, cfg, &sites, span, stats);
        out.sort_by_key(|x| x.0);
        out
    }

    #[test]
    fn peaks_called_per_run_match_one_chromosome_wide_array() {
        // Clustered ends with gaps on both sides of the split distance, over a range of kernel,
        // window, threshold and position-floor settings.
        let mut state = 0x2545_f491_4f6c_dd1du64;
        let mut next = move |n: u64| { state ^= state << 13; state ^= state >> 7; state ^= state << 17; state % n };
        for trial in 0..300 {
            let mut pos = Vec::new();
            let mut at = 1_000i64;
            for _ in 0..(1 + next(12)) {
                at += [1, 5, 30, 44, 45, 46, 60, 400, 5_000][next(9) as usize];
                let width = 1 + next(15) as i64;
                for _ in 0..(1 + next(60)) {
                    pos.push(at + next(width as u64) as i64);
                }
            }
            let c = counts(&pos);
            let cfg = ClusterConfig {
                cluster_window: [20, 5, 0, 30][trial % 4],
                gaussian_std: [1.0, 2.5, 0.7, 4.0][(trial / 4) % 4],
                extent_cutoff: [10, 1, 0, 0][(trial / 16) % 4],
                ratio_cutoff: [0.05, 0.3, 0.05, 0.0][(trial / 16) % 4],
                ..ClusterConfig::default()
            };
            for min_reads in [1, 2, 3] {
                let (mut a, mut b) = (PeakGateStats::default(), PeakGateStats::default());
                assert_eq!(call_peaks(&c, &cfg, min_reads, &mut a), dense_peaks(&c, &cfg, min_reads, &mut b),
                           "trial {trial}, floor {min_reads}");
                assert_eq!(a.suppressed_by_neighbour, b.suppressed_by_neighbour);
            }
        }
    }

    fn counts(positions: &[i64]) -> PositionCounts {
        let mut result = PositionCounts::default();
        for p in positions {
            let r = result.entry(*p).or_insert_with(|| (0, HashMap::default()));
            r.0 += 1; *r.1.entry("sample".into()).or_default() += 1;
        }
        result
    }
    #[test]
    fn compact_multiread_starts_without_singleton_peaks() {
        assert!(compact_start_candidates(&counts(&[100]), 3).is_empty());
        assert!(compact_start_candidates(&counts(&[100, 102]), 3).is_empty());
        let peaks = compact_start_candidates(&counts(&[100, 102, 105]), 3);
        assert_eq!(peaks.len(), 1); assert_eq!(peaks[0].3, 3);
        assert_eq!(peaks[0].2, 102);
    }
    #[test]
    fn diffuse_staircase_does_not_seed_compact_peaks() {
        let p: Vec<_> = (0..300).collect();
        assert!(compact_start_candidates(&counts(&p), 3).is_empty());
        let p: Vec<_> = (0..30).map(|p| p*10).collect();
        assert!(compact_start_candidates(&counts(&p), 3).is_empty());
    }
}
