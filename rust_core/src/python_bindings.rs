// PyO3 bindings exposing the IsoForge engine and post-processing steps to Python.

use std::path::PathBuf;

use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

use crate::annotation_engine::IsoForgeEngine as RustEngine;
use crate::post_processing::{self, PostDB, FilterConfig};
use crate::types::{AnnotStringPool, ClusterConfig, RunConfig};

/// Python-visible IsoForgeEngine wrapper.
#[pyclass]
pub struct IsoForgeEngine {
    inner: RustEngine,
}

#[pymethods]
impl IsoForgeEngine {
    /// Create engine from Parquet annotation directory. `annot_dir` holds the Parquet reference files and
    /// `output_dir` takes the output. `config` is a dict of configuration parameters; optional keys: genome_fasta,
    /// sj_files (list of STAR SJ.out.tab paths).
    #[new]
    fn new(annot_dir: &str, config: &Bound<'_, PyDict>, output_dir: &str) -> PyResult<Self> {
        let rust_config = dict_to_config(config)?;

        // Configure rayon global thread pool to use the requested thread count
        let n_threads = rust_config.threads.max(1);
        rayon::ThreadPoolBuilder::new()
            .num_threads(n_threads)
            .build_global()
            .ok(); // Ignore error if already initialized

        // Initialize logger based on verbosity level
        init_logger(rust_config.verbosity);

        let mut engine = RustEngine::from_parquet(
            &PathBuf::from(annot_dir),
            rust_config,
            PathBuf::from(output_dir),
        )
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;

        // Set genome reference if provided
        if let Some(val) = config.get_item("genome_fasta")? {
            let genome_path: String = val.extract()?;
            engine = engine.with_genome(PathBuf::from(genome_path));
        }

        // Set splice junction database if provided
        if let Some(val) = config.get_item("sj_files")? {
            let sj_paths: Vec<String> = val.extract()?;
            let sj_pathbufs: Vec<PathBuf> = sj_paths.into_iter().map(PathBuf::from).collect();
            engine = engine.with_sj_db(&sj_pathbufs)
                .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
        }

        // Set annotation directory for UTR regions
        engine = engine.with_annot_dir(PathBuf::from(annot_dir));

        // Set clustering config if enabled
        if let Some(cluster_cfg) = dict_to_cluster_config(config)? {
            engine = engine.with_clustering(cluster_cfg);
        }

        Ok(IsoForgeEngine { inner: engine })
    }

    /// Run the IsoForge pipeline on one or more BAM files. `bam_configs` is a list of (bam_path, dataset_name)
    /// tuples. Returns a dict with the status of the run.
    fn run(&self, py: Python<'_>, bam_configs: &Bound<'_, PyList>) -> PyResult<PyObject> {
        let mut configs = Vec::new();

        for item in bam_configs.iter() {
            let tuple = item.extract::<(String, String)>()?;
            configs.push((PathBuf::from(tuple.0), tuple.1));
        }

        // The GIL is released while the Rust pipeline runs.
        py.allow_threads(|| {
            self.inner
                .run(&configs)
                .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))
        })?;

        Python::with_gil(|py| {
            let dict = PyDict::new_bound(py);
            dict.set_item("status", "success")?;
            Ok(dict.into())
        })
    }

}

/// Build the reference index from a GTF file.
#[pyfunction]
#[pyo3(signature = (gtf_path, genome_build, annot_name, output_dir, min_length=0))]
fn init_from_gtf(
    gtf_path: &str,
    genome_build: &str,
    annot_name: &str,
    output_dir: &str,
    min_length: i64,
) -> PyResult<()> {
    crate::gtf_to_parquet::initialize_from_gtf(
        &PathBuf::from(gtf_path),
        genome_build,
        annot_name,
        &PathBuf::from(output_dir),
        min_length,
    )
    .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))
}

/// Get version string.
#[pyfunction]
fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

// ═══════════════════════════════════════════════════════════════════════════════
// POST-PROCESSING PyO3 BINDINGS
// ═══════════════════════════════════════════════════════════════════════════════

/// Load the annotation and results directories into a PostDB handle, which is
/// passed to the other post_* functions and freed with post_free_db().
#[pyfunction]
fn post_load_db(py: Python<'_>, annot_dir: &str, results_dir: &str) -> PyResult<u64> {
    py.allow_threads(|| {
        let db = PostDB::load(annot_dir, results_dir)
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
        let boxed = Box::new(db);
        let ptr = Box::into_raw(boxed) as u64;
        Ok(ptr)
    })
}

/// Keep the annotated form of annotated records that take read-derived ends (see
/// `post_processing::split_annotated_forms`). Call before filtering. Returns the number
/// of records split.
#[pyfunction]
#[pyo3(signature = (handle, known_end_snap=25, known_end_max_extension=25, known_end_extension_min_fraction=0.5))]
fn post_split_annotated_forms(py: Python<'_>, handle: u64, known_end_snap: i64, known_end_max_extension: i64,
                              known_end_extension_min_fraction: f64) -> usize {
    py.allow_threads(|| {
        let db = unsafe { &mut *(handle as *mut PostDB) };
        post_processing::split_annotated_forms(db, known_end_snap, known_end_max_extension, known_end_extension_min_fraction)
    })
}

/// Report under its annotated ID a transcript that passed only as an end variant of its
/// own annotated chain (see `post_processing::reclaim_annotated_records`). Call after
/// filtering; rewrites the whitelist. Returns the number of records reclaimed.
#[pyfunction]
fn post_reclaim_annotated_records(py: Python<'_>, handle: u64, whitelist_file: &str) -> PyResult<usize> {
    py.allow_threads(|| {
        let db = unsafe { &mut *(handle as *mut PostDB) };
        post_processing::reclaim_annotated_records(db, whitelist_file)
    }).map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))
}

/// Free a PostDB handle. Must be called after all post-processing is done.
#[pyfunction]
fn post_free_db(handle: u64) {
    if handle != 0 {
        unsafe {
            let _ = Box::from_raw(handle as *mut PostDB);
        }
    }
}

/// Summarize the datasets (step 1): writes the summary TSV and returns its path.
#[pyfunction]
fn post_summarize(py: Python<'_>, handle: u64, output_prefix: &str) -> PyResult<String> {
    py.allow_threads(|| {
        let db = unsafe { &*(handle as *const PostDB) };
        post_processing::run_summarize(db, output_prefix)
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))
    })
}

/// Filter the transcripts (step 2): writes the whitelist to `output_file` and returns its path.
#[pyfunction]
#[pyo3(signature = (handle, output_file, max_frac_a=0.35, min_count=5, min_datasets=0, allow_genomic=false, filter_known=false, min_count_known=0, min_tss_support=0.0, tss_upstream_rescue=false, evidence_tss_rescue=false, junction_genome=None, junction_sj_files=None, junction_scores=None, min_count_supported=0, require_phased_novel_ends=false, phased_ends_spliced=false, junction_prune_radius=50, junction_prune_ratio=0.7, min_novel_locus_fraction=0.005, frac_a_annotated_end_window=30, novel_exact_support=3, novel_exact_tol=0, ism_evidence_rescue=false, ism_rescue_min_reads=5, ism_rescue_min_distance=0, min_count_nic=5, weak_source_tier=false, sub_path_min_fraction=0.3, sub_path_tail_evidence=false, sibling_edits=true, sibling_ratio_terminal=1.0, sibling_ratio_microexon=0.5, sibling_ratio_shift=0.1, sibling_ratio_intron=0.05, known_end_evidence=2, known_end_evidence_window=50, min_count_end_variant=3, end_variant_min_datasets=0, cap_tss_min_fraction=0.15, cap_tss_min_reads=3, quant_requires_end_support=true, start_model=false, start_min_llr=20.0, admission_model=None, admission_add=0.6, admission_drop=0.5, admission_by_class=Vec::new(), admission_guard_tailed=0.6, admission_guard_tes=0.6, admission_guard_ip=0.1, admission_guard_tss=0.8, admission_guard_tss_all=0.6, admission_min_datasets=2, known_phased_reads=0, pool_datasets=false))]
fn post_filter(
    py: Python<'_>,
    handle: u64,
    output_file: &str,
    max_frac_a: f64,
    min_count: i64,
    min_datasets: i64,
    allow_genomic: bool,
    filter_known: bool,
    min_count_known: i64,
    min_tss_support: f64,
    tss_upstream_rescue: bool,
    evidence_tss_rescue: bool,
    junction_genome: Option<String>,
    junction_sj_files: Option<Vec<String>>,
    junction_scores: Option<String>,
    min_count_supported: i64,
    require_phased_novel_ends: bool,
    phased_ends_spliced: bool,
    junction_prune_radius: i64,
    junction_prune_ratio: f64,
    min_novel_locus_fraction: f64,
    frac_a_annotated_end_window: i64,
    novel_exact_support: i64,
    novel_exact_tol: i64,
    ism_evidence_rescue: bool,
    ism_rescue_min_reads: i64,
    ism_rescue_min_distance: i64,
    min_count_nic: i64,
    weak_source_tier: bool,
    sub_path_min_fraction: f64,
    sub_path_tail_evidence: bool,
    sibling_edits: bool,
    sibling_ratio_terminal: f64,
    sibling_ratio_microexon: f64,
    sibling_ratio_shift: f64,
    sibling_ratio_intron: f64,
    known_end_evidence: i64,
    known_end_evidence_window: i64,
    min_count_end_variant: i64,
    end_variant_min_datasets: i64,
    cap_tss_min_fraction: f64,
    cap_tss_min_reads: i64,
    quant_requires_end_support: bool,
    start_model: bool,
    start_min_llr: f64,
    admission_model: Option<String>,
    admission_add: f64,
    admission_drop: f64,
    admission_by_class: Vec<(String, f64, f64)>,
    admission_guard_tailed: f64,
    admission_guard_tes: f64,
    admission_guard_ip: f64,
    admission_guard_tss: f64,
    admission_guard_tss_all: f64,
    admission_min_datasets: i64,
    known_phased_reads: i64,
    pool_datasets: bool,
) -> PyResult<String> {
    py.allow_threads(|| {
        let db = unsafe { &*(handle as *const PostDB) };
        let config = FilterConfig {
            quant_requires_end_support,
            start_model,
            start_min_llr,
            max_frac_a,
            min_count,
            min_datasets,
            allow_genomic,
            filter_known,
            min_count_known,
            min_tss_support,
            tss_upstream_rescue,
            evidence_tss_rescue,
            junction_genome: junction_genome.clone(),
            junction_scores: junction_scores.clone(),
            junction_sj_files: junction_sj_files.clone().unwrap_or_default(),
            min_count_supported,
            require_phased_novel_ends,
            phased_ends_spliced,
            junction_prune_radius,
            junction_prune_ratio,
            min_novel_locus_fraction,
            frac_a_annotated_end_window,
            novel_exact_support,
            novel_exact_tol,
            ism_evidence_rescue,
            ism_rescue_min_reads,
            ism_rescue_min_distance,
            min_count_nic,
            weak_source_tier,
            sub_path_min_fraction,
            sub_path_tail_evidence,
            sibling_edits,
            sibling_ratio_terminal,
            sibling_ratio_microexon,
            sibling_ratio_shift,
            sibling_ratio_intron,
            known_end_evidence,
            known_end_evidence_window,
            min_count_end_variant,
            end_variant_min_datasets,
            cap_tss_min_fraction,
            cap_tss_min_reads,
            admission_model: admission_model.clone(),
            admission_add,
            admission_drop,
            admission_by_class,
            admission_guard_tailed,
            admission_guard_tes,
            admission_guard_ip,
            admission_guard_tss,
            admission_guard_tss_all,
            admission_min_datasets,
            known_phased_reads,
            pool_datasets,
        };
        post_processing::run_filter(db, output_file, &config)
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))
    })
}

/// Write the abundance table (step 3). Returns its path.
#[pyfunction]
#[pyo3(signature = (handle, output_prefix, genome_build, whitelist_file=None, quant_requires_end_support=true, max_frac_a=0.35, frac_a_annotated_end_window=30, quant_reassign=false, quant_partial_reads=false, quant_hypergraph=true))]
fn post_abundance(
    py: Python<'_>,
    handle: u64,
    output_prefix: &str,
    genome_build: &str,
    whitelist_file: Option<&str>,
    quant_requires_end_support: bool,
    max_frac_a: f64,
    frac_a_annotated_end_window: i64,
    quant_reassign: bool,
    quant_partial_reads: bool,
    quant_hypergraph: bool,
) -> PyResult<String> {
    py.allow_threads(|| {
        let db = unsafe { &*(handle as *const PostDB) };
        post_processing::run_abundance(db, output_prefix, genome_build, whitelist_file, quant_requires_end_support, max_frac_a, frac_a_annotated_end_window, quant_reassign, quant_partial_reads, quant_hypergraph)
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))
    })
}

/// Align read files (FASTA or FASTQ, gzipped or not) to a genome with minimap2, each into a BAM, from one index;
/// with `sort` the BAM is coordinate-sorted and indexed, otherwise unsorted. `jobs` holds (reads, output BAM) pairs.
/// Returns per file (reads, unmapped, primary, supplementary, secondary records).
#[pyfunction]
#[pyo3(signature = (genome, jobs, preset="splice", k=None, forward_only=false, secondary=false, md=true, splice_scores=None, splice_score_scale=0.7, junc_bed=None, copy_comment=false, index_out=None, threads=8, sort=true, tmp_dir=None, sort_mem_gb=4.0))]
fn align_reads(
    py: Python<'_>,
    genome: &str,
    jobs: Vec<(String, String)>,
    preset: &str,
    k: Option<i32>,
    forward_only: bool,
    secondary: bool,
    md: bool,
    splice_scores: Option<String>,
    splice_score_scale: f32,
    junc_bed: Option<String>,
    copy_comment: bool,
    index_out: Option<String>,
    threads: usize,
    sort: bool,
    tmp_dir: Option<String>,
    sort_mem_gb: f64,
) -> PyResult<Vec<(u64, u64, u64, u64, u64)>> {
    py.allow_threads(|| {
        let run = || -> crate::errors::IsoForgeResult<Vec<(u64, u64, u64, u64, u64)>> {
            let opts = crate::align::AlignOptions {
                preset: preset.to_string(), k, forward_only, secondary, md,
                splice_scores: splice_scores.map(std::path::PathBuf::from), splice_score_scale,
                junc_bed: junc_bed.map(std::path::PathBuf::from),
                copy_comment, index_out: index_out.map(std::path::PathBuf::from),
            };
            let pool = rayon::ThreadPoolBuilder::new().num_threads(threads.max(1)).build()
                .map_err(|e| crate::errors::IsoForgeError::ConfigError(e.to_string()))?;
            pool.install(|| {
                let t0 = std::time::Instant::now();
                let aligner = crate::align::Aligner::new(std::path::Path::new(genome), &opts, threads)?;
                eprintln!("[IsoForge] minimap2 index of {} ready in {:.1}s", genome, t0.elapsed().as_secs_f64());
                let mut out = Vec::with_capacity(jobs.len());
                for (reads, bam) in &jobs {
                    let t = std::time::Instant::now();
                    let bam_path = std::path::Path::new(bam);
                    let st = if sort {
                        let tmp = tmp_dir.as_ref().map(std::path::PathBuf::from)
                            .unwrap_or_else(|| bam_path.parent().map(|p| p.to_path_buf()).unwrap_or_default());
                        aligner.align_to_sorted_bam(std::path::Path::new(reads), bam_path, &tmp, threads,
                                                    (sort_mem_gb.max(0.01) * 1e9) as usize)?
                    } else {
                        aligner.align_to_bam(std::path::Path::new(reads), bam_path, threads)?
                    };
                    eprintln!("[IsoForge] aligned {} reads of {} in {:.1}s: {} primary, {} supplementary, {} unmapped -> {}",
                              st.reads, reads, t.elapsed().as_secs_f64(), st.primary, st.supplementary, st.unmapped, bam);
                    out.push((st.reads, st.unmapped, st.primary, st.supplementary, st.secondary));
                }
                Ok(out)
            })
        };
        run().map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))
    })
}

/// Call fusions from the chimeric-read segments kept during annotation.
/// Returns (TSV path, candidates listed, candidates passing).
#[pyfunction]
#[pyo3(signature = (handle, output_prefix, min_reads=2, min_reads_intraexonic=4, min_mapq=20, readthrough_distance=1000000, min_fraction=0.01, max_partners=3, single_dataset_min_reads=5, foldback_min_reads=5, foldback_gap=1000, hub_dominance=20.0, hub_dominance_min_reads=10, genome_fasta=None, rescue_bams=None, pseudogene_partners=true, unnamed_partners=true, min_reads_unnamed=5, min_reads_unnamed_intraexonic=10))]
fn post_fusions(
    py: Python<'_>,
    handle: u64,
    output_prefix: &str,
    min_reads: usize,
    min_reads_intraexonic: usize,
    min_mapq: u8,
    readthrough_distance: i64,
    min_fraction: f64,
    max_partners: usize,
    single_dataset_min_reads: usize,
    foldback_min_reads: usize,
    foldback_gap: i64,
    hub_dominance: f64,
    hub_dominance_min_reads: usize,
    genome_fasta: Option<&str>,
    rescue_bams: Option<Vec<(String, String)>>,
    pseudogene_partners: bool,
    unnamed_partners: bool,
    min_reads_unnamed: usize,
    min_reads_unnamed_intraexonic: usize,
) -> PyResult<(String, usize, usize)> {
    py.allow_threads(|| {
        let db = unsafe { &*(handle as *const PostDB) };
        let run = || -> crate::errors::IsoForgeResult<(String, usize, usize)> {
            let path = db.observed_path.parent().map(|p| p.join("chimeric_segments.parquet")).unwrap_or_default();
            let mut segs = crate::fusion::load_segments(&path)?;
            crate::fusion::name_pseudoautosomal(&mut segs, min_mapq);
            let opts = crate::fusion::FusionOptions { min_reads, min_reads_intraexonic, min_mapq, readthrough_distance,
                                                      min_fraction, max_partners, single_dataset_min_reads, foldback_min_reads, foldback_gap, hub_dominance,
                                                      hub_dominance_min_reads, pseudogene_partners, unnamed_partners, min_reads_unnamed, min_reads_unnamed_intraexonic,
                                                      genome: genome_fasta.map(std::path::PathBuf::from),
                                                      rescue: rescue_bams.filter(|b| !b.is_empty()).map(|b| crate::fusion_rescue::RescueOptions {
                                                          bams: b.into_iter().map(|(d, p)| (d, std::path::PathBuf::from(p))).collect(),
                                                          min_mapq, ..Default::default() }),
                                                      ..Default::default() };
            // read parts the aligner left unaligned, placed on annotated exons, join the reads' aligned parts
            if let (Some(ro), Some(genome)) = (opts.rescue.as_ref(), opts.genome.as_ref()) {
                let po = crate::fusion_place::PlaceOptions { bams: ro.bams.clone(), min_mapq, ..Default::default() };
                let (placed, replaced) = crate::fusion_place::place(db, &segs, &po, genome);
                if !replaced.is_empty() {
                    let mut gone: rustc_hash::FxHashMap<&str, Vec<(&str, i32, i32)>> = rustc_hash::FxHashMap::default();
                    for (d, n, a, b) in &replaced { gone.entry(n.as_str()).or_default().push((d.as_str(), *a, *b)); }
                    segs.retain(|s| !gone.get(s.read_name.as_str()).is_some_and(|v| v.contains(&(s.dataset.as_str(), s.query_start, s.query_end))));
                }
                segs.extend(placed);
            }
            let (calls, multi, genes, datasets) = crate::fusion::call_fusions(db, &segs, opts);
            let listed = calls.iter().filter(|c| c.reads >= min_reads).count();
            let (tsv, n_pass) = crate::fusion::write_fusions(output_prefix, &calls, &genes, &datasets, min_reads)?;
            let (_, n_multi) = crate::fusion::write_multi_fusions(output_prefix, &multi, &calls, &genes, &datasets)?;
            crate::fusion::write_fusion_models(output_prefix, &calls, &multi, &genes, min_reads)?;
            crate::fusion::write_fusion_parquet(output_prefix, &calls, &multi, &genes, &datasets, &segs, min_reads)?;
            crate::fusion::write_fusion_gtf(output_prefix, &calls, &multi, &genes, min_reads)?;
            let n_pass = n_pass + n_multi;
            Ok((tsv, listed, n_pass))
        };
        run().map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))
    })
}

/// Fit every read's own intron chain to every candidate model that holds it (the read hypergraph) and write the
/// fits to `out_path`. Returns (read classes, fits, classes with a fit); (0, 0, 0) when the `introns` column is
/// absent.
#[pyfunction]
#[pyo3(signature = (handle, out_path, tol=10, slack=20, end_slack=100, end_bin=20))]
fn post_hyper_fits(py: Python<'_>, handle: u64, out_path: &str, tol: i64, slack: i64, end_slack: i64, end_bin: i64) -> PyResult<(usize, usize, usize)> {
    py.allow_threads(|| {
        let db = unsafe { &*(handle as *const PostDB) };
        let run = || -> crate::errors::IsoForgeResult<(usize, usize, usize)> {
            let Some(chains) = crate::hyper_core::load_read_chains(&db.observed_path)? else { return Ok((0, 0, 0)) };
            // the fits alone: offsets from the chains' own ends, no sites
            let opts = crate::hyper_core::HyperOptions { tol, slack, end_slack, end_bin, read_ends: false };
            let core = crate::hyper_core::HyperCore::build(db, &chains, &crate::hyper_core::Sites::default(), opts);
            let n = core.write_fits(std::path::Path::new(out_path), &chains)?;
            let with_fit = (0..core.offsets.len() - 1).filter(|&i| core.offsets[i + 1] > core.offsets[i]).count();
            Ok((core.offsets.len() - 1, n, with_fit))
        };
        run().map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))
    })
}

/// Read hypergraph of the reported models (<prefix>_isoforge_hypergraph.tsv) and splicing coordination between
/// pairs of alternative segments (<prefix>_isoforge_coordination.tsv).
/// Returns (coordination TSV path, hyperedges, segment pairs tested, pairs with q < 0.05).
#[pyfunction]
#[pyo3(signature = (handle, output_prefix, whitelist_file=None, min_reads=20.0))]
fn post_hypergraph(
    py: Python<'_>,
    handle: u64,
    output_prefix: &str,
    whitelist_file: Option<&str>,
    min_reads: f64,
) -> PyResult<(String, usize, usize, usize)> {
    py.allow_threads(|| {
        let db = unsafe { &*(handle as *const PostDB) };
        let run = || -> crate::errors::IsoForgeResult<(String, usize, usize, usize)> {
            let reported: Vec<i64> = match whitelist_file.filter(|p| std::path::Path::new(p).exists()) {
                Some(p) => { let mut v: Vec<i64> = post_processing::load_whitelist(p)?.into_iter().collect(); v.sort_unstable(); v }
                None => db.transcripts.iter().map(|t| t.transcript_id).collect(),
            };
            let hg = crate::hypergraph::ReadHypergraph::build(db, &reported);
            let coord = hg.coordination(min_reads);
            let (_, cpath) = crate::hypergraph::write_hypergraph(db, &hg, output_prefix, &coord)?;
            let n_sig = coord.iter().filter(|c| c.q < 0.05).count();
            Ok((cpath, hg.hyperedges.len(), coord.len(), n_sig))
        };
        run().map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))
    })
}

/// Write the GTF (step 4). Returns its path.
#[pyfunction]
#[pyo3(signature = (handle, output_prefix, genome_build, whitelist_file=None, observed_only=false, known_end_snap=25, known_end_max_extension=25, known_end_extension_min_fraction=0.5, known_ends_as_annotated=false))]
fn post_gtf(
    py: Python<'_>,
    handle: u64,
    output_prefix: &str,
    genome_build: &str,
    whitelist_file: Option<&str>,
    observed_only: bool,
    known_end_snap: i64,
    known_end_max_extension: i64,
    known_end_extension_min_fraction: f64,
    known_ends_as_annotated: bool,
) -> PyResult<String> {
    py.allow_threads(|| {
        let db = unsafe { &*(handle as *const PostDB) };
        post_processing::run_gtf(db, output_prefix, genome_build, whitelist_file, observed_only, known_end_snap,
                                 known_end_max_extension, known_end_extension_min_fraction, known_ends_as_annotated)
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))
    })
}

/// Write the read annotations (step 5). Returns (tsv_path, parquet_path).
#[pyfunction]
#[pyo3(signature = (handle, output_prefix, genome_build, whitelist_file=None, split_unsupported=false, max_frac_a=0.35, frac_a_annotated_end_window=30, emit_tsv=true))]
fn post_read_annot(
    py: Python<'_>,
    handle: u64,
    output_prefix: &str,
    genome_build: &str,
    whitelist_file: Option<&str>,
    split_unsupported: bool,
    max_frac_a: f64,
    frac_a_annotated_end_window: i64,
    emit_tsv: bool,
) -> PyResult<(String, String)> {
    py.allow_threads(|| {
        let db = unsafe { &*(handle as *const PostDB) };
        post_processing::run_read_annot(db, output_prefix, genome_build, whitelist_file, split_unsupported, max_frac_a, frac_a_annotated_end_window, emit_tsv)
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))
    })
}

/// Prune the cluster BED files by the whitelist (step 6).
#[pyfunction]
fn post_prune_clusters(py: Python<'_>, handle: u64, whitelist_file: &str) -> PyResult<()> {
    py.allow_threads(|| {
        let db = unsafe { &*(handle as *const PostDB) };
        post_processing::run_prune_clusters(db, whitelist_file)
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))
    })
}

/// Python module definition.
#[pymodule]
fn isoforge_core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Initialize logger with default INFO level
    init_logger(1);

    m.add_class::<IsoForgeEngine>()?;
    m.add_function(wrap_pyfunction!(init_from_gtf, m)?)?;
    m.add_function(wrap_pyfunction!(version, m)?)?;
    m.add_function(wrap_pyfunction!(post_load_db, m)?)?;
    m.add_function(wrap_pyfunction!(post_hyper_fits, m)?)?;
    m.add_function(wrap_pyfunction!(post_free_db, m)?)?;
    m.add_function(wrap_pyfunction!(post_split_annotated_forms, m)?)?;
    m.add_function(wrap_pyfunction!(post_reclaim_annotated_records, m)?)?;
    m.add_function(wrap_pyfunction!(post_summarize, m)?)?;
    m.add_function(wrap_pyfunction!(post_filter, m)?)?;
    m.add_function(wrap_pyfunction!(post_abundance, m)?)?;
    m.add_function(wrap_pyfunction!(align_reads, m)?)?;
    m.add_function(wrap_pyfunction!(post_fusions, m)?)?;
    m.add_function(wrap_pyfunction!(post_hypergraph, m)?)?;
    m.add_function(wrap_pyfunction!(post_gtf, m)?)?;
    m.add_function(wrap_pyfunction!(post_read_annot, m)?)?;
    m.add_function(wrap_pyfunction!(post_prune_clusters, m)?)?;
    Ok(())
}

// ─── Config conversion ──────────────────────────────────────────────────────

fn dict_to_config(dict: &Bound<'_, PyDict>) -> PyResult<RunConfig> {
    let mut config = RunConfig::default();
    if let Some(v) = dict.get_item("repair_microexons")? { config.repair_microexons = v.extract()?; }
    if let Some(v) = dict.get_item("terminal_exon_rescue")? { config.terminal_exon_rescue = v.extract()?; }
    if let Some(v) = dict.get_item("sequence_settled_sites")? { config.sequence_settled_sites = v.extract()?; }
    if let Some(v) = dict.get_item("splice_realign_window")? { config.splice_realign_window = v.extract()?; }
    if let Some(v) = dict.get_item("approx_chain_tolerance")? { config.approx_chain_tolerance = v.extract()?; }
    if let Some(v) = dict.get_item("chain_first")? { config.chain_first = v.extract()?; }
    if let Some(v) = dict.get_item("approx_chain_lock_annotated")? { config.approx_chain_lock_annotated = v.extract()?; }
    if let Some(v) = dict.get_item("assign_unique_ism")? {
        let s: String = v.extract()?;
        config.assign_unique_ism = crate::types::AssignUniqueIsm::parse(&s)
            .ok_or_else(|| pyo3::exceptions::PyValueError::new_err(format!("unknown assign_unique_ism mode: {s}")))?;
    }
    if let Some(v) = dict.get_item("ism_reassign_max_reads")? { config.ism_reassign_max_reads = v.extract()?; }
    if let Some(v) = dict.get_item("nested_fsm_ambiguous")? { config.nested_fsm_ambiguous = v.extract()?; }
    if let Some(v) = dict.get_item("nested_fsm_end_window")? { config.nested_fsm_end_window = v.extract()?; }

    if let Some(val) = dict.get_item("threads")? {
        config.threads = val.extract::<usize>()?;
    }
    if let Some(val) = dict.get_item("min_mapq")? {
        config.min_mapq = val.extract::<u8>()?;
    }
    if let Some(val) = dict.get_item("quant_min_mapq")? {
        config.quant_min_mapq = val.extract::<u8>()?;
    }
    if let Some(v) = dict.get_item("small_rna")? { config.small_rna = v.extract()?; }
    if let Some(v) = dict.get_item("small_rna_types")? { config.small_rna_types = v.extract()?; }
    if let Some(v) = dict.get_item("small_rna_min_read_length")? { config.small_rna_min_read_length = v.extract()?; }
    if let Some(v) = dict.get_item("small_rna_max_length")? { config.small_rna_max_length = v.extract()?; }
    if let Some(val) = dict.get_item("detect_fusions")? {
        config.detect_fusions = val.extract::<bool>()?;
    }
    if let Some(val) = dict.get_item("fusion_spliced_joins")? {
        config.fusion_spliced_joins = val.extract::<bool>()?;
    }
    if let Some(val) = dict.get_item("junction_hmm")? {
        config.junction_hmm = val.extract::<bool>()?;
    }
    if let Some(val) = dict.get_item("junction_hmm_window")? {
        config.junction_hmm_window = val.extract::<i64>()?;
    }
    if let Some(val) = dict.get_item("junction_hmm_min_support")? {
        config.junction_hmm_min_support = val.extract::<u32>()?;
    }
    if let Some(val) = dict.get_item("junction_hmm_model")? {
        config.junction_hmm_model = val.extract::<Option<String>>()?;
    }
    if let Some(val) = dict.get_item("junction_indel_min")? {
        config.junction_indel_min = val.extract::<u32>()?;
    }
    if let Some(val) = dict.get_item("junction_indel_max")? {
        config.junction_indel_max = val.extract::<u32>()?;
    }
    if let Some(val) = dict.get_item("batch_size")? {
        config.batch_size = val.extract::<usize>()?;
    }
    if let Some(val) = dict.get_item("min_coverage")? {
        config.min_coverage = val.extract::<f64>()?;
    }
    if let Some(val) = dict.get_item("min_length")? {
        config.min_length = val.extract::<i64>()?;
    }
    if let Some(val) = dict.get_item("coverage_excludes_polya")? {
        config.coverage_excludes_polya = val.extract::<bool>()?;
    }
    if let Some(val) = dict.get_item("novel_site_window")? {
        config.novel_site_window = val.extract::<i64>()?;
    }
    if let Some(val) = dict.get_item("splice_scores")? {
        config.splice_scores = val.extract::<Option<String>>()?;
    }
    if let Some(val) = dict.get_item("junction_consensus_window")? {
        config.junction_consensus_window = val.extract::<i64>()?;
    }
    if let Some(val) = dict.get_item("min_identity")? {
        config.min_identity = val.extract::<f64>()?;
    }
    if let Some(val) = dict.get_item("max_frac_as")? {
        config.max_frac_as = val.extract::<f64>()?;
    }
    if let Some(val) = dict.get_item("max_5p_dist")? {
        config.cutoff_5p = val.extract::<i64>()?;
    }
    if let Some(val) = dict.get_item("max_3p_dist")? {
        config.cutoff_3p = val.extract::<i64>()?;
    }
    if let Some(val) = dict.get_item("max_splice_dist")? {
        config.max_splice_dist = val.extract::<i64>()?;
    }
    if let Some(val) = dict.get_item("sj_correction_window")? {
        config.sj_correction_window = val.extract::<Option<i64>>()?;
        if config.sj_correction_window.is_some_and(|window| window < 0) {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "sj_correction_window must be nonnegative or None",
            ));
        }
    }
    if let Some(val) = dict.get_item("read_window")? {
        config.read_window = val.extract::<usize>()?.max(1);
    }
    if let Some(val) = dict.get_item("expected_reads")? {
        config.expected_reads = val.extract::<usize>()?;
    }
    if let Some(val) = dict.get_item("id_prefix")? {
        config.id_prefix = val.extract::<String>()?;
    }
    if let Some(val) = dict.get_item("n_places")? {
        config.n_places = val.extract::<usize>()?;
    }
    if let Some(val) = dict.get_item("create_novel_spliced_genes")? {
        config.create_novel_spliced_genes = val.extract::<bool>()?;
    }
    if let Some(val) = dict.get_item("use_cb_tag")? {
        config.use_cb_tag = val.extract::<bool>()?;
    }
    if let Some(val) = dict.get_item("verbosity")? {
        config.verbosity = val.extract::<u8>()?;
    }
    if let Some(val) = dict.get_item("correct_sj")? {
        config.correct_sj = val.extract::<bool>()?;
    }
    if let Some(val) = dict.get_item("correct_microindels")? {
        config.correct_microindels = val.extract::<bool>()?;
    }
    if let Some(val) = dict.get_item("max_indel_size")? {
        config.max_indel_size = val.extract::<u32>()?;
    }
    if let Some(val) = dict.get_item("max_extended_snap")? {
        config.max_extended_snap = val.extract::<i64>()?;
    }
    if let Some(val) = dict.get_item("novel_gene_bin")? {
        config.novel_gene_bin = val.extract::<i64>()?;
    }
    if let Some(val) = dict.get_item("fsm_terminal_new_id")? {
        config.fsm_terminal_new_id = val.extract::<bool>()?;
    }
    if let Some(val) = dict.get_item("strict_full_length")? {
        config.strict_full_length = val.extract::<bool>()?;
    }
    if let Some(val) = dict.get_item("infer_strand_from_polya")? {
        config.infer_strand_from_polya = val.extract::<bool>()?;
    }

    // Rebuild annotation string pool from the (potentially updated) id_prefix
    config.annot_pool = AnnotStringPool::new(&config.id_prefix);

    Ok(config)
}

/// Extract optional TSS/TES clustering configuration from Python dict. Clustering is enabled when "clustering" is
/// True or any clustering-specific key is present; absent keys keep the ClusterConfig defaults.
fn dict_to_cluster_config(dict: &Bound<'_, PyDict>) -> PyResult<Option<ClusterConfig>> {
    // Check if clustering is explicitly enabled
    let enabled = if let Some(val) = dict.get_item("clustering")? {
        val.extract::<bool>()?
    } else {
        // Check if any clustering-specific key is present
        dict.get_item("cluster_window")?.is_some()
            || dict.get_item("cluster_extent_cutoff")?.is_some()
            || dict.get_item("tss_tes_distance")?.is_some()
    };

    if !enabled {
        return Ok(None);
    }

    let mut cfg = ClusterConfig::default();
    if let Some(v) = dict.get_item("reference_end_support")? { cfg.reference_end_support = v.extract()?; }
    if let Some(v) = dict.get_item("known_singleton_ends")? { cfg.known_singleton_ends = v.extract()?; }
    if let Some(v) = dict.get_item("reference_end_window")? { cfg.reference_end_window = v.extract()?; }
    if let Some(v) = dict.get_item("ism_intact_end_support")? { cfg.ism_intact_end_support = v.extract()?; }
    if let Some(v) = dict.get_item("known_ends_not_degraded")? { cfg.known_ends_not_degraded = v.extract()?; }
    if let Some(v) = dict.get_item("compact_tss")? { cfg.compact_tss = v.extract()?; }
    if let Some(v) = dict.get_item("tss_peak_raw_mode")? { cfg.tss_peak_raw_mode = v.extract()?; }
    if let Some(v) = dict.get_item("annotated_id_follows_ends")? { cfg.annotated_id_follows_ends = v.extract()?; }
    if let Some(v) = dict.get_item("end_contradiction_distance")? { cfg.end_contradiction_distance = v.extract()?; }
    if let Some(v) = dict.get_item("annotated_tes_polya_window")? { cfg.annotated_tes_polya_window = v.extract()?; }
    if let Some(v) = dict.get_item("exclude_upstream_unclustered_5p")? { cfg.exclude_upstream_unclustered_5p = v.extract()?; }
    if let Some(v) = dict.get_item("split_capped_tss")? { cfg.split_capped_tss = v.extract()?; }
    if let Some(v) = dict.get_item("split_tss_min_distance")? { cfg.split_tss_min_distance = v.extract()?; }
    if let Some(v) = dict.get_item("split_tss_min_reads")? { cfg.split_tss_min_reads = v.extract()?; }
    if let Some(v) = dict.get_item("split_tss_min_share")? { cfg.split_tss_min_share = v.extract()?; }
    if let Some(v) = dict.get_item("split_tss_cap_ratio")? { cfg.split_tss_cap_ratio = v.extract()?; }

    if let Some(val) = dict.get_item("cluster_window")? {
        cfg.cluster_window = val.extract::<usize>()?;
    }
    if let Some(val) = dict.get_item("cluster_extent_cutoff")? {
        cfg.extent_cutoff = val.extract::<u64>()?;
    }
    if let Some(val) = dict.get_item("cluster_ratio_cutoff")? {
        cfg.ratio_cutoff = val.extract::<f64>()?;
    }
    if let Some(val) = dict.get_item("tss_tes_distance")? {
        cfg.tss_tes_distance = val.extract::<i64>()?;
    }
    if let Some(val) = dict.get_item("min_replication_rate")? {
        cfg.min_replication_rate = val.extract::<f64>()?;
    }
    if let Some(val) = dict.get_item("gaussian_std")? {
        cfg.gaussian_std = val.extract::<f64>()?;
    }
    if let Some(val) = dict.get_item("merge_distance")? {
        cfg.merge_distance = val.extract::<i64>()?;
    }
    if let Some(val) = dict.get_item("min_cluster_reads")? {
        cfg.min_cluster_reads = val.extract::<u64>()?;
    }
    if let Some(val) = dict.get_item("min_condition_reads")? {
        cfg.min_condition_reads = val.extract::<u64>()?;
    }
    if let Some(val) = dict.get_item("min_tss_position_reads")? {
        cfg.min_tss_position_reads = val.extract::<u64>()?;
    }
    if let Some(val) = dict.get_item("min_cluster_reads_known")? {
        cfg.min_cluster_reads_known = val.extract::<u64>()?;
    }
    if let Some(val) = dict.get_item("tss_peak_window")? {
        cfg.tss_peak_window = val.extract::<i64>()?;
    }
    if let Some(val) = dict.get_item("tss_peak_max_entropy")? {
        cfg.tss_peak_max_entropy = val.extract::<f64>()?;
    }
    if let Some(val) = dict.get_item("tss_peak_rescue_min_reads")? { cfg.tss_peak_rescue_min_reads = val.extract::<u64>()?; }
    if let Some(val) = dict.get_item("tss_peak_rescue_enrichment")? { cfg.tss_peak_rescue_enrichment = val.extract::<f64>()?; }
    if let Some(val) = dict.get_item("tss_peak_min_ratio")? {
        cfg.tss_peak_min_ratio = val.extract::<f64>()?;
    }
    if let Some(val) = dict.get_item("tss_peak_max_width")? {
        cfg.tss_peak_max_width = val.extract::<i64>()?;
    }
    if let Some(val) = dict.get_item("filter_internal_priming")? {
        cfg.filter_internal_priming = val.extract::<bool>()?;
    }
    if let Some(val) = dict.get_item("filter_truncated_reads")? {
        cfg.filter_truncated_reads = val.extract::<bool>()?;
    }
    if let Some(val) = dict.get_item("write_end_tracks")? {
        cfg.write_end_tracks = val.extract::<bool>()?;
    }
    if let Some(val) = dict.get_item("filter_degraded_reads")? {
        cfg.filter_degraded_reads = val.extract::<bool>()?;
    }
    if let Some(val) = dict.get_item("degradation_tolerance")? {
        cfg.degradation_tolerance = val.extract::<i64>()?;
    }
    if let Some(val) = dict.get_item("cage_peaks_file")? {
        cfg.cage_peaks_file = Some(val.extract::<String>()?);
    }
    if let Some(val) = dict.get_item("cage_window")? {
        cfg.cage_window = val.extract::<i64>()?;
    }
    if let Some(val) = dict.get_item("require_polya_tail")? {
        cfg.require_polya_tail = val.extract::<bool>()?;
    }

    // ── Satellite cluster absorption ─────────────────────────────────────────
    if let Some(val) = dict.get_item("satellite_absorption_distance")? {
        cfg.satellite_absorption_distance = val.extract::<i64>()?;
    }
    if let Some(val) = dict.get_item("satellite_max_usage")? {
        cfg.satellite_max_usage = val.extract::<f64>()?;
    }

    Ok(Some(cfg))
}

/// Sets the log level: 0 = errors and warnings, 1 = info, 2+ = debug. The logger is created on the first call.
fn init_logger(verbosity: u8) {
    use std::sync::Once;
    static INIT: Once = Once::new();

    INIT.call_once(|| {
        env_logger::Builder::new()
            .filter_level(log::LevelFilter::Debug)
            .format_timestamp_secs()
            .format_target(false)
            .init();
    });
    log::set_max_level(match verbosity {
        0 => log::LevelFilter::Warn,
        1 => log::LevelFilter::Info,
        _ => log::LevelFilter::Debug,
    });
}
