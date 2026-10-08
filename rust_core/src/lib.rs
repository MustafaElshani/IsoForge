// IsoForge core: Rust engine for long-read transcript annotation,
// quantification and TSS/TES analysis.

// ── Global allocator: jemalloc on Linux ──────────────────────────────────────
#[cfg(target_os = "linux")]
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

// Freed pages are returned to the OS immediately; overridable at runtime via
// _RJEM_MALLOC_CONF.
#[cfg(target_os = "linux")]
#[allow(non_upper_case_globals)]
#[export_name = "_rjem_malloc_conf"]
pub static MALLOC_CONF: &[u8] = b"dirty_decay_ms:0,muzzy_decay_ms:0\0";

pub mod types;
pub mod data_structures;
pub mod annotation_hypergraph;
pub mod bam_reader;
pub mod strand_evidence;
pub mod approximate_chain;
pub mod consistent_assignment;
pub mod microexon_repair;
pub mod splice_realignment;
pub mod terminal_rescue;
pub mod small_rna;
pub mod splice_scores;
pub mod novel_sites;
pub mod reference_end_support;
pub mod genome_ref;
pub mod transcript_utils;
pub mod vertex_matching;
pub mod edge_matching;
pub mod transcript_classifier;
pub mod gene_assignment;
pub mod parquet_io;
pub mod gtf_to_parquet;
pub mod tss_tes_clustering;
pub mod cage_peaks;
pub mod annotation_engine;
pub mod errors;
pub mod post_processing;
pub mod quant_em;
pub mod align;
pub mod junction_hmm;
pub mod junction_groups;
pub mod fusion;
pub mod fusion_rescue;
pub mod fusion_place;
pub mod hypergraph;
pub mod hyper_core;
pub mod model_hypergraph;
pub mod admission;
pub mod python_bindings;

pub use annotation_engine::IsoForgeEngine;
pub use types::*;
