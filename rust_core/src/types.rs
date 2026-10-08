// Core types used throughout IsoForge.

use rustc_hash::FxHashSet;
use serde::{Deserialize, Serialize};
use std::fmt;
use std::sync::Arc;

// ─── ID types ─────────────────────────────────────────────────────────────────

pub type GeneId = i64;
pub type TranscriptId = i64;
pub type VertexId = i64;
pub type EdgeId = i64;
pub type ObservedId = i64;
pub type DatasetId = i64;

// ─── Strand ───────────────────────────────────────────────────────────────────

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub enum Strand {
    Plus,
    Minus,
}

impl Strand {
    pub fn from_char(c: char) -> Option<Self> {
        match c {
            '+' => Some(Strand::Plus),
            '-' => Some(Strand::Minus),
            _ => None,
        }
    }

    pub fn from_str(s: &str) -> Option<Self> {
        match s {
            "+" => Some(Strand::Plus),
            "-" => Some(Strand::Minus),
            _ => None,
        }
    }

    pub fn opposite(&self) -> Self {
        match self {
            Strand::Plus => Strand::Minus,
            Strand::Minus => Strand::Plus,
        }
    }

    pub fn as_str(&self) -> &'static str {
        match self {
            Strand::Plus => "+",
            Strand::Minus => "-",
        }
    }
}

impl fmt::Display for Strand {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.as_str())
    }
}

// ─── Edge type ────────────────────────────────────────────────────────────────

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum EdgeType {
    Exon,
    Intron,
}

impl EdgeType {
    pub fn from_str(s: &str) -> Option<Self> {
        match s {
            "exon" => Some(EdgeType::Exon),
            "intron" => Some(EdgeType::Intron),
            _ => None,
        }
    }

    pub fn as_str(&self) -> &'static str {
        match self {
            EdgeType::Exon => "exon",
            EdgeType::Intron => "intron",
        }
    }
}

// ─── Novelty Category ─────────────────────────────────────────────────────────

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum NoveltyCategory {
    Known,
    FSM, // Full Splice Match
    ISM, // Incomplete Splice Match
    NIC, // Novel In Catalog
    NNC, // Novel Not in Catalog
    Antisense,
    Intergenic,
    Genomic,
    Fusion,
}

impl NoveltyCategory {
    pub fn as_str(&self) -> &'static str {
        match self {
            NoveltyCategory::Known => "Known",
            NoveltyCategory::FSM => "FSM",
            NoveltyCategory::ISM => "ISM",
            NoveltyCategory::NIC => "NIC",
            NoveltyCategory::NNC => "NNC",
            NoveltyCategory::Antisense => "Antisense",
            NoveltyCategory::Intergenic => "Intergenic",
            NoveltyCategory::Genomic => "Genomic",
            NoveltyCategory::Fusion => "Fusion",
        }
    }
}

// ─── ISM Subtype ──────────────────────────────────────────────────────────────

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum IsmSubtype {
    Prefix,
    Suffix,
    Internal,
    Both, // prefix and suffix
}

// ─── Core Data Records ───────────────────────────────────────────────────────

/// A genomic location (vertex) in the annotation graph
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Location {
    pub location_id: VertexId,
    pub genome_build: String,
    pub chromosome: String,
    pub position: i64,
}

/// An edge (exon or intron) connecting two vertices
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Edge {
    pub edge_id: EdgeId,
    pub v1: VertexId,
    pub v2: VertexId,
    pub edge_type: EdgeType,
    pub strand: Strand,
}

/// A known transcript model
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TranscriptModel {
    pub transcript_id: TranscriptId,
    pub gene_id: GeneId,
    pub start_exon: EdgeId,
    pub jn_path: Option<String>,    // comma-separated internal edge IDs
    pub end_exon: EdgeId,
    pub start_vertex: VertexId,
    pub end_vertex: VertexId,
    pub n_exons: i32,
    pub chromosome: String,
    pub start_pos: i64,
    pub end_pos: i64,
}

/// A gene record
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Gene {
    pub gene_id: GeneId,
    pub strand: Strand,
    pub chromosome: String,
    pub start: i64,
    pub end: i64,
}

/// An observed read record (output row)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ObservedRecord {
    pub obs_id: ObservedId,
    pub gene_id: GeneId,
    pub transcript_id: TranscriptId,
    pub read_name: String,
    pub dataset: String,
    pub start_vertex: VertexId,
    pub end_vertex: VertexId,
    pub start_exon: EdgeId,
    pub end_exon: EdgeId,
    pub start_delta: Option<i64>,
    pub end_delta: Option<i64>,
    pub read_length: i64,
    pub fraction_as: Option<f64>,
    pub custom_label: Option<String>,
    pub allelic_label: Option<String>,
    pub start_support: Option<String>,
    pub end_support: Option<String>,
}

/// Annotation novelty entry (gene, transcript, or exon); string fields are
/// interned `Arc<str>`.
#[derive(Debug, Clone)]
pub struct NoveltyAnnotation {
    pub id: i64,
    pub annot_name: Arc<str>,
    pub source: Arc<str>,
    pub attribute: Arc<str>,
    pub value: Arc<str>,
}

/// QC metrics for a single read
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QcMetrics {
    pub dataset: String,
    pub read_id: String,
    pub passed_qc: bool,
    pub primary_mapped: bool,
    pub read_length: i64,
    pub fraction_aligned: Option<f64>,
    pub identity: Option<f64>,
}

// ═════════════════════════════════════════════════════════════════════════════
// DENSE READ INDEX
// ═════════════════════════════════════════════════════════════════════════════

/// Dense per-read index, assigned in BAM order at parse time and stable across
/// all input BAMs; the key that joins reads between annotation, clustering and
/// filtering.
pub type ReadIdx = u32;

/// Sentinel for "no index assigned" (reads dropped before indexing).
pub const READ_IDX_NONE: ReadIdx = u32::MAX;

/// Bitset over dense read indices. Grows on demand; `contains` on an
/// out-of-range index returns false.
#[derive(Debug, Clone, Default)]
pub struct ReadIdxSet {
    bits: Vec<u64>,
    count: usize,
}

impl ReadIdxSet {
    pub fn new() -> Self {
        Self::default()
    }

    /// Pre-allocate for `n` reads.
    pub fn with_capacity(n: usize) -> Self {
        Self { bits: vec![0u64; n.div_ceil(64)], count: 0 }
    }

    #[inline]
    pub fn insert(&mut self, idx: ReadIdx) -> bool {
        if idx == READ_IDX_NONE {
            return false;
        }
        let (word, mask) = (idx as usize / 64, 1u64 << (idx as usize % 64));
        if word >= self.bits.len() {
            self.bits.resize(word + 1, 0);
        }
        let was_set = self.bits[word] & mask != 0;
        self.bits[word] |= mask;
        if !was_set {
            self.count += 1;
        }
        !was_set
    }

    #[inline]
    pub fn contains(&self, idx: ReadIdx) -> bool {
        if idx == READ_IDX_NONE {
            return false;
        }
        let (word, mask) = (idx as usize / 64, 1u64 << (idx as usize % 64));
        self.bits.get(word).is_some_and(|w| w & mask != 0)
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.count
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.count == 0
    }

    /// Union another set into this one.
    pub fn union_with(&mut self, other: &ReadIdxSet) {
        if other.bits.len() > self.bits.len() {
            self.bits.resize(other.bits.len(), 0);
        }
        let mut count = 0usize;
        for (dst, src) in self.bits.iter_mut().zip(other.bits.iter()) {
            *dst |= *src;
        }
        for w in &self.bits {
            count += w.count_ones() as usize;
        }
        self.count = count;
    }
}

// ─── Hint flags for ISM / Known classification ───────────────────────────────
pub const HINT_IS_ISM: u8        = 0x01;
pub const HINT_IS_ISM_PREFIX: u8 = 0x02;
pub const HINT_IS_ISM_SUFFIX: u8 = 0x04;
pub const HINT_IS_KNOWN: u8      = 0x08;
/// ISM read consistent with exactly one annotated transcript (see
/// `consistent_assignment`); `HINT_ISM_5P_INTACT`/`3P_INTACT` say which of
/// its ends coincide with that parent's chain ends.
pub const HINT_ISM_UNIQUE: u8    = 0x10;
/// Consistent with several annotated transcripts, so it confirms none. On an
/// FSM read (no `HINT_IS_ISM`) it marks a short isoform's read that could
/// equally be a fragment of a longer isoform containing it.
pub const HINT_ISM_AMBIGUOUS: u8 = 0x20;
pub const HINT_ISM_5P_INTACT: u8 = 0x40;
pub const HINT_ISM_3P_INTACT: u8 = 0x80;

/// Whether the flags mark an FSM read that could be a fragment of a longer
/// isoform containing it: `HINT_ISM_AMBIGUOUS` set without `HINT_IS_ISM`.
pub fn is_nested_ambiguous(flags: u8) -> bool {
    flags & (HINT_ISM_AMBIGUOUS | HINT_IS_ISM) == HINT_ISM_AMBIGUOUS
}

/// When a truncated (ISM) read consistent with exactly one annotated
/// transcript is attributed to that transcript instead of its own model.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "snake_case")]
pub enum AssignUniqueIsm {
    #[default]
    Off,
    /// Every unique-parent ISM read.
    Always,
    /// Only reads of ISM models with at most `ism_reassign_max_reads` reads.
    BelowNovelFloor,
    /// Only when the parent has no full-length read of its own.
    NoFsmParent,
}

impl AssignUniqueIsm {
    pub fn parse(s: &str) -> Option<Self> {
        match s {
            "off" => Some(Self::Off),
            "always" => Some(Self::Always),
            "below_novel_floor" => Some(Self::BelowNovelFloor),
            "no_fsm_parent" => Some(Self::NoFsmParent),
            _ => None,
        }
    }
}

/// Result of annotating a single read
#[derive(Debug, Clone)]
pub struct ReadAnnotation {
    pub gene_id: GeneId,
    pub transcript_id: TranscriptId,
    pub gene_novelty: Vec<NoveltyAnnotation>,
    pub transcript_novelty: Vec<NoveltyAnnotation>,
    pub exon_novelty: Vec<NoveltyAnnotation>,
    pub start_vertex: VertexId,
    pub end_vertex: VertexId,
    pub start_exon: EdgeId,
    pub end_exon: EdgeId,
    pub start_delta: Option<i64>,
    pub end_delta: Option<i64>,
    // Read metadata
    pub read_id: String,
    /// Dense run-stable index for this read; the join key used by clustering.
    pub read_idx: ReadIdx,
    /// Shared per-dataset name.
    pub dataset: Arc<str>,
    pub location: String,
    pub strand: Strand,
    pub read_length: i64,
    /// Actual BAM alignment position at 5' end (positions[0])
    pub start_pos: i64,
    /// Actual BAM alignment position at 3' end (positions[last])
    pub end_pos: i64,
    pub n_exons: i32,
    /// Largest distance (bp) an internal splice site of this read was moved
    /// by vertex snapping; 0 = every site exact or created novel.
    pub max_splice_shift: i32,
    /// The read's own internal exon boundaries in 5'->3' order (each intron's
    /// donor, then its acceptor), as corrected before they were matched to the
    /// index's sites: the read's hyperedge, whichever model it is counted for.
    pub introns: Vec<i64>,
    /// The read's 5' end carries untemplated G's in its soft clip.
    pub cap_signature: bool,
    /// The alignment ends before the transcript's 3' end (see `ParsedRead::unaligned_3p`).
    pub unaligned_3p: bool,
    pub fraction_as: Option<f64>,
    /// PolyA tail length from the basecaller `pt:i` tag (None if absent)
    pub polya_tail_len: Option<i64>,
    /// PolyA tail length derived from soft-clipped bases.
    pub polya_tail_softclip: Option<i64>,
    pub custom_label: Option<String>,
    pub allelic_label: Option<String>,
    /// TSS peak position; formatted only when written out.
    pub start_support: Option<i64>,
    /// TES/polyA peak position; formatted only when written out.
    pub end_support: Option<i64>,
    // ── PAS-aware internal priming ───────────────────────────────────────
    pub internal_priming: Option<InternalPrimingResult>,
    // ── TSS/TES cluster join result ──────────────────────────────────────
    pub end_site_info: Option<ReadEndSiteInfo>,
    // ── Classification hint flags ──
    pub hint_flags: u8,
    /// The one annotated transcript a truncated read is consistent with.
    pub ism_parent: Option<TranscriptId>,
}
#[derive(Debug, Clone, Default)]
pub struct StartEndInfo {
    pub start_vertex: VertexId,
    pub end_vertex: VertexId,
    pub start_exon: EdgeId,
    pub end_exon: EdgeId,
    pub diff_5p: Option<i64>,
    pub diff_3p: Option<i64>,
    pub start_novelty: u8, // 0 = known, 1 = novel
    pub end_novelty: u8,
    pub vertex_ids: Vec<VertexId>,
    pub edge_ids: Vec<EdgeId>,
}

/// Edge key for lookups
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct EdgeKey {
    pub v1: VertexId,
    pub v2: VertexId,
    pub edge_type: EdgeType,
}

/// Pre-populated `Arc<str>` pool for NoveltyAnnotation strings; read-only
/// after construction.
#[derive(Debug, Clone)]
pub struct AnnotStringPool {
    annot_name: Arc<str>,
    source: Arc<str>,
    pool: std::collections::HashMap<&'static str, Arc<str>>,
}

impl Default for AnnotStringPool {
    fn default() -> Self {
        // The id_prefix of RunConfig::default().
        Self::new("ISOFOR")
    }
}

impl AnnotStringPool {
    pub fn new(annot_name: &str) -> Self {
        let mut pool = std::collections::HashMap::new();
        for s in &[
            "TRUE", "FALSE", "NOVEL", "KNOWN",
            "ISM_transcript", "NIC_transcript", "NNC_transcript",
            "FSM_transcript", "ISM-prefix_transcript", "ISM-suffix_transcript",
            "antisense_transcript", "antisense_gene",
            "genomic_transcript", "intergenic_transcript", "intergenic_novel",
            "gene_status", "transcript_status", "exon_status",
            "gene_name", "gene_id", "transcript_name", "transcript_id",
            "ISM_to_IDs", "ISM-prefix_to_IDs", "ISM-suffix_to_IDs",
            "FSM_to_IDs", "FSM_novel_5p", "FSM_novel_3p",
            "gene_antisense_to_IDs", "fusion_transcript", "TSS_TES_variant",
        ] {
            pool.insert(*s, Arc::from(*s));
        }
        Self {
            annot_name: Arc::from(annot_name),
            source: Arc::from("IsoForge"),
            pool,
        }
    }

    /// Get a cached `Arc<str>` for a known string, or allocate a new one.
    #[inline]
    pub fn get(&self, s: &str) -> Arc<str> {
        self.pool.get(s).map(Arc::clone).unwrap_or_else(|| Arc::from(s))
    }

    /// Get the shared annot_name (id_prefix).
    #[inline]
    pub fn annot_name(&self) -> Arc<str> {
        Arc::clone(&self.annot_name)
    }

    /// Get the shared source string ("IsoForge").
    #[inline]
    pub fn source(&self) -> Arc<str> {
        Arc::clone(&self.source)
    }
}

/// Run configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RunConfig {
    /// Sequence-supported recovery of annotated microexons.
    #[serde(default)]
    pub repair_microexons: bool,
    /// Realign a read's unannotated intron onto an annotated intron sharing one
    /// site when the other lies within this many bp and the read's own bases do
    /// not favour the observed intron by at least two edits (0 disables).
    #[serde(default = "default_splice_realign_window")]
    pub splice_realign_window: i64,
    /// Place a terminal exon the aligner left in a soft clip back on the read when the clip matches the exon beyond
    /// the annotated intron at that end (see `terminal_rescue`).
    #[serde(default = "default_true")]
    pub terminal_exon_rescue: bool,
    /// Keep a read's unannotated splice site where it lies when the read's bases favour its intron over every
    /// annotated intron nearby by at least two edits: the site consensus and vertex snapping leave it in place.
    #[serde(default)]
    pub sequence_settled_sites: bool,
    /// Tolerance (bp) for chain-level approximate assignment to an annotated
    /// transcript; 0 disables it.
    #[serde(default)]
    pub approx_chain_tolerance: i64,
    /// Try chain-level assignment before the NIC path, so a read whose only
    /// novelty is a splice site within tolerance of one annotated transcript is
    /// attributed to it instead of minting a model. Needs a tolerance > 0.
    #[serde(default)]
    pub chain_first: bool,
    /// A read splice site sitting exactly on an annotated site of its role may
    /// not be moved onto a different annotated site.
    #[serde(default = "default_true")]
    pub approx_chain_lock_annotated: bool,
    /// Write one row per chain-level decision to chain_assignment.tsv.
    #[serde(default)]
    pub chain_trace: bool,
    /// Attribute a truncated read consistent with exactly one annotated
    /// transcript to that transcript.
    #[serde(default)]
    pub assign_unique_ism: AssignUniqueIsm,
    /// `below_novel_floor`: dissolve only ISM models with at most this many reads in total.
    #[serde(default = "default_ism_reassign_max_reads")]
    pub ism_reassign_max_reads: usize,
    /// An FSM read of a short isoform whose end lies inside a longer isoform's
    /// exon and far from the short isoform's annotated end confirms neither.
    #[serde(default)]
    pub nested_fsm_ambiguous: bool,
    /// How far (bp) from the short isoform's annotated end a read must lie to
    /// count as ambiguous under `nested_fsm_ambiguous`.
    #[serde(default = "default_nested_fsm_end_window")]
    pub nested_fsm_end_window: i64,
    pub build: String,
    #[serde(default = "default_min_mapq")]
    pub min_mapq: u8,
    /// Reads with MAPQ in [quant_min_mapq, min_mapq) are parsed but build no
    /// model; they are kept for quantification only and counted towards the
    /// reported models they fit. `quant_min_mapq >= min_mapq` disables this.
    #[serde(default = "default_quant_min_mapq")]
    pub quant_min_mapq: u8,
    /// Annotated small RNAs counted, not assembled (small_rna.rs): unspliced reads lying inside the one-exon
    /// transcripts of genes of `small_rna_types` are kept from `small_rna_min_read_length` bases and whatever their
    /// mapping quality, count for those transcripts only and build no model.
    #[serde(default)]
    pub small_rna: bool,
    #[serde(default = "default_small_rna_types")]
    pub small_rna_types: Vec<String>,
    #[serde(default = "default_small_rna_min_read_length")]
    pub small_rna_min_read_length: i64,
    #[serde(default = "default_small_rna_max_length")]
    pub small_rna_max_length: i64,
    /// Keep the segments of chimeric reads (records with an SA tag) for fusion calling.
    #[serde(default)]
    pub detect_fusions: bool,
    /// Junction groups decided by pair-HMM likelihoods of all their reads (junction_groups.rs): every dataset is read
    /// twice before annotation, to count introns and to score the reads of contested groups.
    #[serde(default)]
    pub junction_hmm: bool,
    /// bp at each site within which introns form one group.
    #[serde(default = "default_junction_hmm_window")]
    pub junction_hmm_window: i64,
    /// Reads an intron needs to be a candidate site of its group (annotated introns are candidates without).
    #[serde(default = "default_junction_hmm_min_support")]
    pub junction_hmm_min_support: u32,
    /// Pair-HMM parameters to use instead of training them (a `<prefix>_junction_hmm.json` written by a run).
    #[serde(default)]
    pub junction_hmm_model: Option<String>,
    /// With detect_fusions, a spliced alignment is also cut at an intron that joins two annotated genes, and its
    /// parts are kept as segments for fusion calling.
    #[serde(default)]
    pub fusion_spliced_joins: bool,
    /// Introns with an indel of this many bases or more next to them are placed where the read needs no
    /// indel (see `place_indel_adjacent_introns`); 0 disables.
    #[serde(default)]
    pub junction_indel_min: u32,
    #[serde(default = "default_junction_indel_max")]
    pub junction_indel_max: u32,
    pub min_coverage: f64,
    /// Longest soft clip, other than a poly(A)/poly(T) tail, tolerated at a read end.
    /// Shorter clips are left out of the coverage measure; a read with a longer one is dropped.
    #[serde(default = "default_max_terminal_clip")]
    pub max_terminal_clip: u32,
    pub min_identity: f64,
    pub min_length: i64,
    pub cutoff_5p: i64,
    pub cutoff_3p: i64,
    pub use_cb_tag: bool,
    pub create_novel_spliced_genes: bool,
    pub id_prefix: String,
    pub n_places: usize,
    pub threads: usize,
    pub tmp_dir: String,
    pub out_prefix: String,
    pub batch_size: usize,
    /// Reads parsed and annotated per window before being dropped.
    #[serde(default = "default_read_window")]
    pub read_window: usize,
    /// Optional total read count, used to size per-read vectors once.
    #[serde(default)]
    pub expected_reads: usize,
    /// Width (bp) of the genomic bin that groups intergenic reads onto a shared
    /// novel gene. 0 mints a separate gene per distinct start position.
    #[serde(default = "default_novel_gene_bin")]
    pub novel_gene_bin: i64,
    /// Give an FSM read whose terminals exceed the cutoffs a new transcript ID
    /// instead of refining the annotated transcript's ends.
    #[serde(default)]
    pub fsm_terminal_new_id: bool,
    /// Require every read to pass end validation before it counts towards a
    /// transcript's depth; otherwise all reads of a validated transcript count.
    #[serde(default)]
    pub strict_full_length: bool,
    /// Take the transcript strand from a polyA tail in the soft clips instead of
    /// the read's alignment orientation.
    #[serde(default = "default_infer_strand_from_polya")]
    pub infer_strand_from_polya: bool,
    pub max_frac_as: f64,
    pub max_splice_dist: i64,
    #[serde(default)]
    pub sj_correction_window: Option<i64>,
    pub verbosity: u8,
    /// Correct non-canonical splice junctions during BAM reading.
    pub correct_sj: bool,
    /// Correct microindels (≤ max_indel_size bp) in CIGAR strings before
    /// extracting splice junctions.
    pub correct_microindels: bool,
    /// Maximum size (bp) of insertions/deletions to correct.
    pub max_indel_size: u32,
    /// Maximum distance (bp) for extended reference-only splice vertex snapping;
    /// novel vertices within it are snapped to the reference unless short-read
    /// evidence validates them. 0 disables.
    pub max_extended_snap: i64,
    /// Place a novel splice vertex at the most-supported raw position within
    /// this many bp among the window's reads instead of the first read's
    /// position (0 disables).
    pub junction_consensus_window: i64,
    /// Novel splice positions (no annotated site of their role within
    /// `max_extended_snap`) move onto a stronger position within this many bp
    /// holding at least twice their reads (0 disables; see `novel_sites`).
    #[serde(default = "default_novel_site_window")]
    pub novel_site_window: i64,
    /// Splice-site scores (TSV, optionally gzipped) weighing the strength of
    /// each novel splice position.
    #[serde(default)]
    pub splice_scores: Option<String>,
    /// Exclude a poly(A)/poly(T) soft clip from the coverage denominator.
    pub coverage_excludes_polya: bool,
    /// String pool for NoveltyAnnotation; not serialized, rebuilt from id_prefix.
    #[serde(skip)]
    pub annot_pool: AnnotStringPool,
}

impl Default for RunConfig {
    fn default() -> Self {
        let id_prefix = "ISOFOR".to_string();
        let annot_pool = AnnotStringPool::new(&id_prefix);
        Self {
            repair_microexons: false,
            splice_realign_window: 30,
            terminal_exon_rescue: true,
            sequence_settled_sites: false,
            approx_chain_tolerance: 50,
            chain_first: true,
            approx_chain_lock_annotated: true,
            chain_trace: false,
            assign_unique_ism: AssignUniqueIsm::BelowNovelFloor,
            ism_reassign_max_reads: default_ism_reassign_max_reads(),
            nested_fsm_ambiguous: false,
            nested_fsm_end_window: default_nested_fsm_end_window(),
            build: String::new(),
            min_mapq: default_min_mapq(),
            quant_min_mapq: default_quant_min_mapq(),
            small_rna: false,
            small_rna_types: default_small_rna_types(),
            small_rna_min_read_length: default_small_rna_min_read_length(),
            small_rna_max_length: default_small_rna_max_length(),
            detect_fusions: false,
            fusion_spliced_joins: false,
            junction_hmm: false,
            junction_hmm_window: default_junction_hmm_window(),
            junction_hmm_min_support: default_junction_hmm_min_support(),
            junction_hmm_model: None,
            junction_indel_min: 0,
            junction_indel_max: default_junction_indel_max(),
            min_coverage: 0.9,
            max_terminal_clip: default_max_terminal_clip(),
            min_identity: 0.8,
            min_length: 200,
            cutoff_5p: 500,
            cutoff_3p: 200,
            use_cb_tag: false,
            create_novel_spliced_genes: false,
            id_prefix,
            n_places: 9,
            threads: 4,
            tmp_dir: "isoforge_tmp/".to_string(),
            out_prefix: "isoforge".to_string(),
            batch_size: 1000,
            read_window: default_read_window(),
            expected_reads: 0,
            novel_gene_bin: default_novel_gene_bin(),
            fsm_terminal_new_id: false,
            strict_full_length: false,
            infer_strand_from_polya: default_infer_strand_from_polya(),
            max_frac_as: 0.5,
            max_splice_dist: 6,
            sj_correction_window: None,
            verbosity: 1,
            correct_sj: true,
            correct_microindels: true,
            max_indel_size: 5,
            max_extended_snap: 25,
            junction_consensus_window: 6,
            novel_site_window: default_novel_site_window(),
            splice_scores: None,
            coverage_excludes_polya: true,
            annot_pool,
        }
    }
}

/// Shared empty string for placeholder `Arc<str>` fields.
pub fn empty_arc_str() -> Arc<str> {
    static EMPTY: std::sync::OnceLock<Arc<str>> = std::sync::OnceLock::new();
    Arc::clone(EMPTY.get_or_init(|| Arc::from("")))
}

fn default_novel_site_window() -> i64 {
    20
}

fn default_read_window() -> usize {
    500_000
}

fn default_novel_gene_bin() -> i64 {
    1_000
}

fn default_infer_strand_from_polya() -> bool {
    true
}

fn default_min_mapq() -> u8 {
    16
}

fn default_quant_min_mapq() -> u8 {
    u8::MAX
}

/// Gene types of the small-RNA rule.
pub fn default_small_rna_types() -> Vec<String> {
    ["rRNA", "snRNA", "snoRNA", "scaRNA", "misc_RNA", "vault_RNA", "Y_RNA", "Mt_tRNA", "ribozyme", "sRNA", "scRNA"]
        .iter().map(|s| s.to_string()).collect()
}

fn default_small_rna_min_read_length() -> i64 { 50 }

fn default_small_rna_max_length() -> i64 { 500 }

fn default_junction_indel_max() -> u32 {
    12
}

fn default_split_tss_min_distance() -> i64 { 20 }
fn default_split_tss_min_reads() -> usize { 50 }
fn default_split_tss_min_share() -> f64 { 0.1 }
fn default_tss_peak_rescue_min_reads() -> u64 { 10 }
fn default_junction_hmm_window() -> i64 { 5 }
fn default_junction_hmm_min_support() -> u32 { 2 }
fn default_tss_peak_rescue_enrichment() -> f64 { 5.0 }
fn default_split_tss_cap_ratio() -> f64 { 0.7 }

fn default_max_terminal_clip() -> u32 {
    300
}

/// Transcript dict key: sorted edge IDs.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct TranscriptPathKey(pub Vec<EdgeId>);

impl TranscriptPathKey {
    pub fn new(mut edges: Vec<EdgeId>) -> Self {
        edges.sort_unstable();
        Self(edges)
    }

    pub fn from_sorted(edges: Vec<EdgeId>) -> Self {
        Self(edges)
    }

    /// Check if self is a subset of other
    pub fn is_subset_of(&self, other: &TranscriptPathKey) -> bool {
        if self.0.len() > other.0.len() {
            return false;
        }
        let mut j = 0;
        for &edge in &self.0 {
            while j < other.0.len() && other.0[j] < edge {
                j += 1;
            }
            if j >= other.0.len() || other.0[j] != edge {
                return false;
            }
            j += 1;
        }
        true
    }

    pub fn contains(&self, edge: EdgeId) -> bool {
        self.0.binary_search(&edge).is_ok()
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn iter(&self) -> impl Iterator<Item = &EdgeId> {
        self.0.iter()
    }

    pub fn as_frozenset(&self) -> FxHashSet<EdgeId> {
        self.0.iter().copied().collect()
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// TSS / TES CLUSTERING TYPES
// ═══════════════════════════════════════════════════════════════════════════════

/// Poly(A) signal (PAS) hexamers (24), scanned upstream of TES clusters to
/// distinguish poly(A) sites from internal priming.
pub const PAS_HEXAMERS: &[&str] = &[
    "AATAAA", "ATTAAA", "AGTAAA", "TATAAA",
    "CATAAA", "GATAAA", "AATATA", "AATACA",
    "AATAGA", "AATGAA", "ACTAAA",
    "AATAAT", "AATAAG", "ATTACA", "ATTATA", "AACAAG",
    "AAGAAA",
    "AACAAA",
    "AAAACA", "AATTAA", "ATGAAA", "ATAAAA",
    "AATCAA", "AATAAC",
];

/// Number of bases upstream of a TES cluster peak to scan for PAS hexamers.
pub const PAS_SCAN_WINDOW: i64 = 100;

/// The strongest poly(A) signals among [`PAS_HEXAMERS`].
pub const CANONICAL_PAS: &[&str] = &["AATAAA", "ATTAAA"];

/// Furthest a canonical hexamer may lie upstream of a cleavage site and still direct
/// it, as bases between the hexamer and the 3' end.
pub const PAS_MAX_GAP: i64 = 40;

/// Furthest a variant hexamer may lie upstream and still direct cleavage.
pub const PAS_VARIANT_MAX_GAP: i64 = 24;

/// Nearest a hexamer may lie to the 3' end and still count.
pub const PAS_MIN_GAP: i64 = 10;

/// Bases between a poly(A) signal and the cleavage site it most often directs.
pub const PAS_TYPICAL_GAP: i64 = 18;

/// The poly(A) signal directing cleavage at the 3' end of `upstream` (sequence in transcript orientation that ends at
/// the 3' end), as (motif, bases between it and the 3' end). A canonical hexamer counts from [`PAS_MIN_GAP`] to
/// [`PAS_MAX_GAP`] nt upstream, a variant only to [`PAS_VARIANT_MAX_GAP`]; a canonical hexamer wins over a variant,
/// then the one nearest [`PAS_TYPICAL_GAP`].
pub fn positioned_pas(upstream: &[u8]) -> Option<(&'static str, i64)> {
    let n = upstream.len();
    let mut best: Option<((bool, i64), &'static str, i64)> = None;
    for i in 0..n.saturating_sub(5) {
        let gap = (n - i - 6) as i64;
        if !(PAS_MIN_GAP..=PAS_MAX_GAP).contains(&gap) {
            continue;
        }
        let hex = &upstream[i..i + 6];
        let Some(&motif) = PAS_HEXAMERS.iter().find(|p| p.as_bytes().eq_ignore_ascii_case(hex)) else { continue };
        if gap > PAS_VARIANT_MAX_GAP && !CANONICAL_PAS.contains(&motif) {
            continue;
        }
        let key = (CANONICAL_PAS.contains(&motif), -(gap - PAS_TYPICAL_GAP).abs());
        if best.map_or(true, |(k, _, _)| key > k) {
            best = Some((key, motif, gap));
        }
    }
    best.map(|(_, motif, gap)| (motif, gap))
}

/// Bases downstream of a cleavage site searched for a U-rich element.
pub const DSE_WINDOW: i64 = 40;

/// Fraction of U (T in the genome) in [`DSE_WINDOW`] that marks a downstream element.
pub const DSE_MIN_U: f64 = 0.4;

/// Label reported for a site whose only poly(A) signal is a downstream element.
pub const DSE_LABEL: &str = "DSE";

/// Whether `downstream` (sequence in transcript orientation starting just after the
/// 3' end) carries a U-rich downstream element.
pub fn u_rich_downstream(downstream: &[u8]) -> bool {
    !downstream.is_empty()
        && downstream.iter().filter(|b| b.eq_ignore_ascii_case(&b'T')).count() as f64 >= DSE_MIN_U * downstream.len() as f64
}

#[cfg(test)]
mod pas_tests {
    use super::*;

    #[test]
    fn downstream_element_is_u_rich() {
        assert!(u_rich_downstream(b"TTTTGTTTTGTTTTCTTTTGCCCCCCCCCCCCCCCCCCCC"));
        assert!(!u_rich_downstream(b"TTGTTTCTGTTTTACGTGTTCCCCCCCCCCCCCCCCCCCC"));
        assert!(!u_rich_downstream(b"CCCCCCCCCCCCCCCCCCCCGGGGGGGGGGAAAAAAAAAA"));
        assert!(!u_rich_downstream(b""));
    }

    fn upstream(hexamers: &[(&str, usize)]) -> Vec<u8> {
        // 100 nt of C, with each hexamer placed `gap` bases before the 3' end.
        let mut s = vec![b'C'; 100];
        for (h, gap) in hexamers {
            let i = 100 - gap - 6;
            s[i..i + 6].copy_from_slice(h.as_bytes());
        }
        s
    }

    #[test]
    fn positioned_signal_wins_over_nearer_or_farther_hexamers() {
        assert_eq!(positioned_pas(&upstream(&[("AATAAA", 18)])), Some(("AATAAA", 18)));
        // A hexamer beyond 40 nt is not this site's signal.
        assert_eq!(positioned_pas(&upstream(&[("AATAAA", 60)])), None);
        // A canonical hexamer beats a variant nearer the typical position.
        assert_eq!(positioned_pas(&upstream(&[("AATAAA", 25), ("CATAAA", 17)])), Some(("AATAAA", 25)));
        // Among variants, the one nearest the typical position.
        assert_eq!(positioned_pas(&upstream(&[("CATAAA", 12), ("GATAAA", 20)])), Some(("GATAAA", 20)));
        assert_eq!(positioned_pas(&upstream(&[])), None);
    }

    #[test]
    fn a_hexamer_counts_where_its_strength_can_place_cleavage() {
        // A hexamer within 9 nt of the end does not count.
        assert_eq!(positioned_pas(&upstream(&[("AAGAAA", 2)])), None);
        assert_eq!(positioned_pas(&upstream(&[("AATAAA", 9)])), None);
        assert_eq!(positioned_pas(&upstream(&[("AATAAA", 10)])), Some(("AATAAA", 10)));
        // A canonical hexamer directs cleavage from further upstream than a variant.
        assert_eq!(positioned_pas(&upstream(&[("AATAAA", 35)])), Some(("AATAAA", 35)));
        assert_eq!(positioned_pas(&upstream(&[("AGTAAA", 35)])), None);
        assert_eq!(positioned_pas(&upstream(&[("AGTAAA", 24)])), Some(("AGTAAA", 24)));
        // A variant too close or too far leaves the site to a canonical hexamer in range.
        assert_eq!(positioned_pas(&upstream(&[("AAGAAA", 2), ("ATTAAA", 33)])), Some(("ATTAAA", 33)));
    }
}

/// Default parameters for TSS/TES clustering.
pub const DEFAULT_CLUSTER_WINDOW: usize = 20;
pub const DEFAULT_CLUSTER_EXTENT_CUTOFF: u64 = 10;
pub const DEFAULT_CLUSTER_RATIO_CUTOFF: f64 = 0.05;
pub const DEFAULT_GAUSSIAN_STD: f64 = 1.0;
pub const DEFAULT_TSS_TES_DISTANCE: i64 = 50;
pub const DEFAULT_MIN_REPLICATION_RATE: f64 = 0.5;
pub const DEFAULT_MERGE_DISTANCE: i64 = 25;
pub const DEFAULT_DEGRADATION_TOLERANCE: i64 = 20;
pub const DEFAULT_MIN_CLUSTER_READS: u64 = 3;
pub const DEFAULT_MIN_CONDITION_READS: u64 = 5;
/// Minimum reads at a single genomic position before it contributes to TSS
/// Gaussian smoothing.
pub const DEFAULT_MIN_TSS_POSITION_READS: u64 = 2;
/// Minimum post-filter reads per cluster before revalidation revokes Known/FSM reads.
pub const DEFAULT_MIN_CLUSTER_READS_KNOWN: u64 = 2;
/// Containing transcripts must carry this many times an annotated transcript's
/// reads before reads that never reach its distinguishing terminal sequence are
/// refused its ID (see `split_transcripts_by_clusters`).
pub const DISTINGUISHING_PARENT_RATIO: usize = 100;
/// Reads a non-dominant TSS/TES group needs before it becomes its own transcript variant.
pub const MIN_END_VARIANT_READS: usize = 2;
fn default_reference_end_window() -> i64 { 25 }
fn default_splice_realign_window() -> i64 {
    30
}

fn default_end_contradiction_distance() -> i64 { 100 }
fn default_annotated_tes_polya_window() -> i64 { 30 }
fn default_true() -> bool { true }
fn default_ism_reassign_max_reads() -> usize { 9 }
fn default_nested_fsm_end_window() -> i64 { 100 }
/// TSS peak gating defaults (window, entropy, peak ratio, width).
pub const DEFAULT_TSS_PEAK_WINDOW: i64 = 10;        // +/- bp around peak to keep
pub const DEFAULT_TSS_PEAK_MAX_ENTROPY: f64 = 0.85;  // normalized entropy (0..1)
pub const DEFAULT_TSS_PEAK_MIN_RATIO: f64 = 0.2;    // peak_count / total_count
pub const DEFAULT_TSS_PEAK_MAX_WIDTH: i64 = 30;     // max span with support (bp)

/// PolyA tail detection parameters.
pub const POLYA_TAIL_MIN_LEN: usize = 10;
pub const POLYA_TAIL_MIN_PERCENT_A: f64 = 0.9;
pub const POLYA_TAIL_SCORE_MATCH: i32 = 1;     // +1 per A (or T on minus)
pub const POLYA_TAIL_SCORE_MISMATCH: i32 = -5;  // -5 per non-A
pub const POLYA_TAIL_SCORE_ABORT: i32 = -20;     // stop if score drops below this

/// Type of end site cluster.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum ClusterType {
    /// Transcript End Site (poly(A) site)
    TES,
    /// Transcript Start Site
    TSS,
}

impl fmt::Display for ClusterType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ClusterType::TES => write!(f, "TES"),
            ClusterType::TSS => write!(f, "TSS"),
        }
    }
}

/// Genomic feature where a cluster falls relative to annotation.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum GenomicFeature {
    ThreePrimeUtr,
    FivePrimeUtr,
    Exon,
    Intron,
    Gene,       // within gene but not a specific feature
    Intergenic, // outside any gene
}

impl GenomicFeature {
    pub fn as_str(&self) -> &'static str {
        match self {
            GenomicFeature::ThreePrimeUtr => "three_prime_utr",
            GenomicFeature::FivePrimeUtr => "five_prime_utr",
            GenomicFeature::Exon => "exon",
            GenomicFeature::Intron => "intron",
            GenomicFeature::Gene => "gene",
            GenomicFeature::Intergenic => "intergenic",
        }
    }

    pub fn from_str(s: &str) -> Self {
        match s {
            "three_prime_utr" => GenomicFeature::ThreePrimeUtr,
            "five_prime_utr" => GenomicFeature::FivePrimeUtr,
            "exon" => GenomicFeature::Exon,
            "intron" => GenomicFeature::Intron,
            "gene" => GenomicFeature::Gene,
            _ => GenomicFeature::Intergenic,
        }
    }

    /// Priority ordering for feature annotation (lower = higher priority).
    /// When a cluster overlaps multiple features, the highest priority wins.
    /// For TES: three_prime_utr > exon > intron > five_prime_utr > gene > intergenic
    /// For TSS: five_prime_utr > exon > intron > three_prime_utr > gene > intergenic
    pub fn priority(&self, cluster_type: ClusterType) -> u8 {
        match cluster_type {
            ClusterType::TES => match self {
                GenomicFeature::ThreePrimeUtr => 0,
                GenomicFeature::Exon => 1,
                GenomicFeature::Intron => 2,
                GenomicFeature::FivePrimeUtr => 3,
                GenomicFeature::Gene => 4,
                GenomicFeature::Intergenic => 5,
            },
            ClusterType::TSS => match self {
                GenomicFeature::FivePrimeUtr => 0,
                GenomicFeature::Exon => 1,
                GenomicFeature::Intron => 2,
                GenomicFeature::ThreePrimeUtr => 3,
                GenomicFeature::Gene => 4,
                GenomicFeature::Intergenic => 5,
            },
        }
    }
}

/// A TSS or TES cluster identified by peak-calling.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EndSiteCluster {
    pub cluster_id: i64,
    pub cluster_type: ClusterType,
    pub chromosome: String,
    pub start: i64,            // cluster interval start (1-based)
    pub end: i64,              // cluster interval end (1-based, inclusive)
    pub peak_site: i64,        // position of the peak within the cluster
    pub strand: Strand,
    pub read_count: u64,       // total reads supporting this cluster
    pub feature: GenomicFeature,
    pub gene_id: Option<String>,  // gene associated with this cluster
    pub tpm: f64,              // transcripts per million
    pub gene_count: u64,       // total reads for the gene at this site type
    pub usage: f64,            // fraction of gene's reads at this cluster

    // TES-specific fields (None for TSS clusters)
    pub frac_a: Option<f64>,        // mean fraction-A of member reads
    pub pas_signal: Option<String>, // PAS hexamer found upstream
    pub pas_distance: Option<i64>,  // distance from peak to PAS hexamer

    // Replication info
    pub n_datasets_supporting: u32, // number of datasets with reads in this cluster
}

/// UTR region from the annotation, used for cluster feature annotation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct UtrRegion {
    pub chromosome: String,
    pub start: i64,
    pub end: i64,
    pub strand: Strand,
    pub feature: GenomicFeature, // FivePrimeUtr or ThreePrimeUtr
    pub gene_id: String,
    pub transcript_id: String,
}

/// TSS/TES clustering configuration.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ClusterConfig {
    /// Compact multi-read TSS seeding.
    #[serde(default)]
    pub compact_tss: bool,
    /// Two-read support for existing annotated ends only.
    #[serde(default)]
    pub reference_end_support: bool,
    /// How far (bp) a read's own ends may sit from the annotated ends and still
    /// count as full-length for reference end support.
    #[serde(default = "default_reference_end_window")]
    pub reference_end_window: i64,
    /// Let a truncated read attributed to its unique parent support the
    /// parent's annotated end on the side the read kept intact. Requires
    /// reference_end_support.
    #[serde(default)]
    pub ism_intact_end_support: bool,
    /// A full-length read of an annotated transcript whose end lies within
    /// `reference_end_window` of the annotated end is never called degraded on
    /// that side.
    #[serde(default)]
    pub known_ends_not_degraded: bool,
    /// Let one full-length exact-FSM read confirm an already-annotated end;
    /// requires reference_end_support and never applies to novel ends.
    #[serde(default)]
    pub known_singleton_ends: bool,
    /// Judge each TSS cluster at its most-used base rather than at the Gaussian
    /// maximum, and report that base as the TSS.
    #[serde(default)]
    pub tss_peak_raw_mode: bool,
    /// When the reads of an annotated transcript split into TSS/TES groups, the
    /// annotated ID stays with the group whose cluster peaks match the annotated
    /// TSS and TES rather than with the largest group.
    #[serde(default)]
    pub annotated_id_follows_ends: bool,
    /// How far (bp) a read group's end cluster must lie outside an annotated end
    /// before the group is treated as a different molecule.
    #[serde(default = "default_end_contradiction_distance")]
    pub end_contradiction_distance: i64,
    /// A read ending within this many bp of any annotated TES on its strand is
    /// never treated as internally primed (0 disables).
    #[serde(default = "default_annotated_tes_polya_window")]
    pub annotated_tes_polya_window: i64,
    /// A read whose 5' end joined no TSS cluster and starts upstream of its
    /// transcript's main TSS is marked end-unsupported and kept out of the
    /// transcript's count.
    #[serde(default)]
    pub exclude_upstream_unclustered_5p: bool,
    /// A start peak within `tss_tes_distance` of its model's dominant start, on the same TES cluster, stays an end
    /// variant of its own with at least `split_tss_min_distance` bp between them, `split_tss_min_reads` reads,
    /// `split_tss_min_share` of the model's reads and `split_tss_cap_ratio` times the dominant's cap-signature rate.
    /// A library without cap signatures splits nothing.
    #[serde(default)]
    pub split_capped_tss: bool,
    #[serde(default = "default_split_tss_min_distance")]
    pub split_tss_min_distance: i64,
    #[serde(default = "default_split_tss_min_reads")]
    pub split_tss_min_reads: usize,
    #[serde(default = "default_split_tss_min_share")]
    pub split_tss_min_share: f64,
    #[serde(default = "default_split_tss_cap_ratio")]
    pub split_tss_cap_ratio: f64,
    pub cluster_window: usize,
    pub extent_cutoff: u64,
    pub ratio_cutoff: f64,
    pub gaussian_std: f64,
    pub tss_tes_distance: i64,
    pub min_replication_rate: f64,
    pub min_replication_samples: u32,
    /// Merge nearby clusters whose peaks are ≤ merge_distance bp apart
    pub merge_distance: i64,
    /// Minimum total reads per cluster to keep
    pub min_cluster_reads: u64,
    /// Minimum total reads that keep a cluster below the replication rate
    pub min_condition_reads: u64,
    /// Exclude internally primed reads from TES counting
    pub filter_internal_priming: bool,
    /// Exclude ISM-suffix reads from TSS and ISM-prefix reads from TES counting
    pub filter_truncated_reads: bool,
    /// Exclude degraded reads (5'/3' end deviates from per-transcript mode) from TSS/TES counting
    pub filter_degraded_reads: bool,
    /// Write per-sample read-end coverage tracks (bedGraph, and BigWig when bedGraphToBigWig is on PATH).
    pub write_end_tracks: bool,
    /// Max bp deviation from per-transcript modal 5'/3' position to consider a read non-degraded
    pub degradation_tolerance: i64,
    /// Enable TSS/TES clustering and joining
    pub enabled: bool,
    /// Optional path to CAGE peak BED file (plain or .gz) for 5′ end validation.
    /// When provided, reads whose 5′ end does not overlap a CAGE peak are marked
    /// unsupported and excluded from output.
    pub cage_peaks_file: Option<String>,
    /// Window (bp) around a CAGE peak for overlap matching
    pub cage_window: i64,
    /// Require poly(A) tail evidence (BAM pt:i tag or soft-clip detection) for a
    /// read to be considered full-length.
    pub require_polya_tail: bool,
    /// Starts are placed by full-length reads: a spliced read's 5' end counts toward
    /// a TSS only when its poly(A) tail is present. Decided per run from the tail rate,
    /// never from the command line.
    #[serde(default)]
    pub tss_from_full_length: bool,
    /// Minimum reads at the same genomic position for that position to contribute
    /// to the TSS Gaussian smoothing signal; applied to TSS clustering only.
    pub min_tss_position_reads: u64,
    /// Minimum post-filter supported reads per cluster before revalidation
    /// revokes Known/FSM reads.
    pub min_cluster_reads_known: u64,

    // ── TSS peak shape gating (entropy/ratio/width) ─────────────────────────
    /// TSS cluster window half-width (bp) around the peak for counting.
    pub tss_peak_window: i64,
    /// Max normalized entropy of the TSS start distribution.
    pub tss_peak_max_entropy: f64,
    /// A peak the entropy or peak-ratio gate rejects is kept, the strongest one per gene, when it holds more reads
    /// than every TSS peak its gene kept (in a gene without one, always), at least this many reads, and its reads per
    /// bp are at least `tss_peak_rescue_enrichment` times those of the 200 bp downstream of it. 0 disables.
    #[serde(default = "default_tss_peak_rescue_min_reads")]
    pub tss_peak_rescue_min_reads: u64,
    #[serde(default = "default_tss_peak_rescue_enrichment")]
    pub tss_peak_rescue_enrichment: f64,
    /// Min peak ratio (peak_count / total_count) for a TSS cluster.
    pub tss_peak_min_ratio: f64,
    /// Max supported width (bp) of a TSS cluster within the window.
    pub tss_peak_max_width: i64,

    // ── TES satellite cluster absorption ─────────────────────────────────────
    /// Maximum distance (bp) from a larger cluster's peak for a minor cluster to
    /// be absorbed as a satellite.
    pub satellite_absorption_distance: i64,
    /// Maximum gene-level usage (fraction) for a cluster to be absorbed as a
    /// satellite.
    pub satellite_max_usage: f64,
}

impl Default for ClusterConfig {
    fn default() -> Self {
        Self {
            compact_tss: false,
            reference_end_support: false,
            reference_end_window: default_reference_end_window(),
            ism_intact_end_support: false,
            known_ends_not_degraded: true,
            known_singleton_ends: false,
            tss_peak_raw_mode: true,
            annotated_id_follows_ends: true,
            end_contradiction_distance: default_end_contradiction_distance(),
            annotated_tes_polya_window: default_annotated_tes_polya_window(),
            exclude_upstream_unclustered_5p: true,
            split_capped_tss: false,
            split_tss_min_distance: default_split_tss_min_distance(),
            split_tss_min_reads: default_split_tss_min_reads(),
            split_tss_min_share: default_split_tss_min_share(),
            split_tss_cap_ratio: default_split_tss_cap_ratio(),
            cluster_window: DEFAULT_CLUSTER_WINDOW,
            extent_cutoff: DEFAULT_CLUSTER_EXTENT_CUTOFF,
            ratio_cutoff: DEFAULT_CLUSTER_RATIO_CUTOFF,
            gaussian_std: DEFAULT_GAUSSIAN_STD,
            tss_tes_distance: DEFAULT_TSS_TES_DISTANCE,
            min_replication_rate: DEFAULT_MIN_REPLICATION_RATE,
            min_replication_samples: 2,
            merge_distance: DEFAULT_MERGE_DISTANCE,
            min_cluster_reads: DEFAULT_MIN_CLUSTER_READS,
            min_condition_reads: DEFAULT_MIN_CONDITION_READS,
            filter_internal_priming: true,
            filter_truncated_reads: true,
            filter_degraded_reads: true,
            write_end_tracks: true,
            degradation_tolerance: DEFAULT_DEGRADATION_TOLERANCE,
            enabled: false,
            cage_peaks_file: None,
            cage_window: 50,
            require_polya_tail: true,
            tss_from_full_length: false,
            min_tss_position_reads: DEFAULT_MIN_TSS_POSITION_READS,
            min_cluster_reads_known: DEFAULT_MIN_CLUSTER_READS_KNOWN,
            tss_peak_window: DEFAULT_TSS_PEAK_WINDOW,
            tss_peak_max_entropy: DEFAULT_TSS_PEAK_MAX_ENTROPY,
            tss_peak_rescue_min_reads: default_tss_peak_rescue_min_reads(),
            tss_peak_rescue_enrichment: default_tss_peak_rescue_enrichment(),
            tss_peak_min_ratio: DEFAULT_TSS_PEAK_MIN_RATIO,
            tss_peak_max_width: DEFAULT_TSS_PEAK_MAX_WIDTH,
            satellite_absorption_distance: 100,
            satellite_max_usage: 0.01,
        }
    }
}

/// Per-read TSS/TES join result (attached to ReadAnnotation after clustering).
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct ReadEndSiteInfo {
    /// TSS cluster ID this read's 5′ end was joined to (None if no match)
    pub tss_cluster_id: Option<i64>,
    /// TES cluster ID this read's 3′ end was joined to (None if no match)
    pub tes_cluster_id: Option<i64>,
    /// TSS peak position
    pub tss_site: Option<i64>,
    /// TES (polyA) peak position
    pub polya_site: Option<i64>,
    /// Feature at the TSS cluster
    pub tss_feature: Option<GenomicFeature>,
    /// Feature at the TES cluster
    pub polya_feature: Option<GenomicFeature>,
    /// PAS hexamer found near the TES cluster
    pub polya_signal: Option<String>,
    /// PolyA tail length from BAM pt:i tag or soft-clip detection
    pub polya_tail_len: Option<i64>,
    /// Whether the read passes the configured gates; a read may pass on its TES
    /// cluster alone.
    pub is_supported: bool,
    /// Reason for rejection (None if supported)
    pub rejection_reason: Option<String>,
}

/// Internal priming assessment: a read is internally primed when its
/// downstream fraction-A exceeds the threshold and no poly(A) signal is
/// positioned to direct cleavage at its 3' end.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct InternalPrimingResult {
    pub fraction_a: f64,
    pub has_pas_signal: bool,
    pub pas_motif: Option<String>,
    pub pas_distance: Option<i64>,
    /// True if this read is likely an internal priming artifact
    /// (high frac-A AND no PAS signal)
    pub is_internal_priming: bool,
}
