// Core data structures. Reference annotation data is shared read-only across
// threads via Arc; novel vertices and edges live in concurrent DashMaps so
// every thread sees the others' novel entities.

use std::sync::atomic::{AtomicI64, Ordering};
use std::sync::Arc;

use coitrees::{COITree, Interval, IntervalTree};
use dashmap::DashMap;
use parking_lot::RwLock;
use rustc_hash::{FxHashMap, FxHashSet};
use crate::types::*;

// ─── Atomic counters ─────────────────────────────────────────────────────────

/// Lock-free atomic counter for ID generation.
#[derive(Debug)]
pub struct AtomicCounter {
    value: AtomicI64,
}

impl AtomicCounter {
    pub fn new(initial: i64) -> Self {
        Self {
            value: AtomicI64::new(initial),
        }
    }

    /// Atomically increment and return the new value.
    pub fn increment(&self) -> i64 {
        self.value.fetch_add(1, Ordering::SeqCst) + 1
    }

    pub fn get(&self) -> i64 {
        self.value.load(Ordering::SeqCst)
    }
}

/// All IsoForge counters grouped together
#[derive(Debug)]
pub struct Counters {
    pub gene: AtomicCounter,
    pub transcript: AtomicCounter,
    pub vertex: AtomicCounter,
    pub edge: AtomicCounter,
    pub observed: AtomicCounter,
    pub dataset: AtomicCounter,
}

impl Counters {
    pub fn new(
        genes: i64,
        transcripts: i64,
        vertices: i64,
        edges: i64,
        observed: i64,
        datasets: i64,
    ) -> Self {
        Self {
            gene: AtomicCounter::new(genes),
            transcript: AtomicCounter::new(transcripts),
            vertex: AtomicCounter::new(vertices),
            edge: AtomicCounter::new(edges),
            observed: AtomicCounter::new(observed),
            dataset: AtomicCounter::new(datasets),
        }
    }
}

// ─── Location Dict ───────────────────────────────────────────────────────────
// Sorted positions per chromosome for exact lookups and range queries.

/// Positions for a single chromosome, sorted for binary search.
#[derive(Debug, Clone)]
pub struct ChromosomePositions {
    /// Sorted array of (position, vertex_id) pairs
    pub positions: Vec<(i64, VertexId)>,
}

impl ChromosomePositions {
    pub fn new() -> Self {
        Self {
            positions: Vec::new(),
        }
    }

    pub fn from_pairs(mut pairs: Vec<(i64, VertexId)>) -> Self {
        pairs.sort_unstable_by_key(|&(pos, _)| pos);
        Self { positions: pairs }
    }

    /// Finalize: sort and deduplicate
    pub fn finalize(&mut self) {
        self.positions.sort_unstable_by_key(|&(pos, _)| pos);
        self.positions.dedup_by_key(|pair| pair.0);
    }

    /// Exact position lookup: O(log n)
    pub fn get(&self, position: i64) -> Option<VertexId> {
        self.positions
            .binary_search_by_key(&position, |&(pos, _)| pos)
            .ok()
            .map(|idx| self.positions[idx].1)
    }

    /// Find the nearest vertex within a range: O(log n + k).
    pub fn find_nearest_in_range(
        &self,
        position: i64,
        window_start: i64,
        window_end: i64,
    ) -> Option<(i64, VertexId)> {
        // Find the insertion point
        let start_idx = self
            .positions
            .partition_point(|&(pos, _)| pos < window_start);
        let end_idx = self
            .positions
            .partition_point(|&(pos, _)| pos <= window_end);

        let mut best: Option<(i64, VertexId, i64)> = None; // (pos, id, abs_dist)

        for &(pos, vid) in &self.positions[start_idx..end_idx] {
            let dist = (pos - position).abs();
            match &best {
                None => best = Some((pos, vid, dist)),
                Some((_, _, best_dist)) if dist < *best_dist => {
                    best = Some((pos, vid, dist));
                }
                _ => {}
            }
        }

        best.map(|(pos, vid, _)| (pos, vid))
    }

    /// Find nearest vertex with directional priority (for permissive matching).
    /// priority_direction: -1 means prefer smaller positions first,
    ///                     +1 means prefer larger positions first.
    pub fn find_nearest_with_priority(
        &self,
        position: i64,
        window_start: i64,
        window_end: i64,
        priority_direction: i32,
    ) -> Option<(i64, VertexId)> {
        let start_idx = self
            .positions
            .partition_point(|&(pos, _)| pos < window_start);
        let end_idx = self
            .positions
            .partition_point(|&(pos, _)| pos <= window_end);

        if start_idx >= end_idx {
            return None;
        }

        // Candidates are sorted by absolute distance; at equal distance the
        // priority direction comes first.
        let mut candidates: Vec<(i64, VertexId, i64)> = self.positions[start_idx..end_idx]
            .iter()
            .map(|&(pos, vid)| {
                let diff = pos - position;
                let abs_dist = diff.abs();
                // aligned: 0 for an exact match or a vertex in the priority direction (smaller positions
                // for -1, larger for +1), 1 for a vertex in the opposite direction
                let aligned = if diff == 0 {
                    0 // exact match always preferred
                } else if (diff < 0 && priority_direction < 0) || (diff > 0 && priority_direction > 0) {
                    0 // in priority direction
                } else {
                    1 // opposite direction
                };
                (pos, vid, abs_dist * 2 + aligned as i64) // composite sort key
            })
            .collect();

        // Sort by composite key: smallest absolute distance first,
        // then priority direction preferred at ties
        candidates.sort_unstable_by_key(|&(_, _, key)| key);

        candidates.first().map(|&(pos, vid, _)| (pos, vid))
    }

    pub fn insert(&mut self, position: i64, vertex_id: VertexId) {
        match self
            .positions
            .binary_search_by_key(&position, |&(pos, _)| pos)
        {
            Ok(idx) => self.positions[idx].1 = vertex_id,
            Err(idx) => self.positions.insert(idx, (position, vertex_id)),
        }
    }

    pub fn len(&self) -> usize {
        self.positions.len()
    }

    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }
}

/// The complete location dictionary, shared read-only across threads.
#[derive(Debug, Clone)]
pub struct LocationDict {
    pub chromosomes: FxHashMap<String, ChromosomePositions>,
}

impl LocationDict {
    pub fn new() -> Self {
        Self {
            chromosomes: FxHashMap::default(),
        }
    }

    pub fn get_vertex(&self, chromosome: &str, position: i64) -> Option<VertexId> {
        self.chromosomes.get(chromosome)?.get(position)
    }

    pub fn find_nearest(
        &self,
        chromosome: &str,
        position: i64,
        window_start: i64,
        window_end: i64,
    ) -> Option<(i64, VertexId)> {
        self.chromosomes
            .get(chromosome)?
            .find_nearest_in_range(position, window_start, window_end)
    }

    pub fn find_nearest_with_priority(
        &self,
        chromosome: &str,
        position: i64,
        window_start: i64,
        window_end: i64,
        priority_direction: i32,
    ) -> Option<(i64, VertexId)> {
        self.chromosomes.get(chromosome)?.find_nearest_with_priority(
            position,
            window_start,
            window_end,
            priority_direction,
        )
    }

    pub fn insert(&mut self, chromosome: String, position: i64, vertex_id: VertexId) {
        self.chromosomes
            .entry(chromosome)
            .or_insert_with(ChromosomePositions::new)
            .insert(position, vertex_id);
    }

    pub fn finalize(&mut self) {
        for chrom in self.chromosomes.values_mut() {
            chrom.finalize();
        }
    }

    pub fn total_vertices(&self) -> usize {
        self.chromosomes.values().map(|c| c.len()).sum()
    }
}

// ─── Edge Dict ───────────────────────────────────────────────────────────────

pub type EdgeDict = FxHashMap<EdgeKey, Edge>;

// ─── Transcript Dict ─────────────────────────────────────────────────────────
// Transcripts keyed by sorted edge path; containment lookups go through the annotation hypergraph
// (annotation_hypergraph.rs).

/// Main transcript dictionary
pub type TranscriptDict = FxHashMap<TranscriptPathKey, TranscriptModel>;

// ─── Vertex to Gene mapping ─────────────────────────────────────────────────
// vertex_id -> set of (gene_id, strand)

pub type Vertex2Gene = FxHashMap<VertexId, FxHashSet<(GeneId, Strand)>>;

// ─── Gene start/end positions ────────────────────────────────────────────────
// gene_id -> position -> vertex_id, used for permissive matching with gene priority

pub type GeneEndpointDict = FxHashMap<GeneId, FxHashMap<i64, VertexId>>;

// ─── Interval tree for monoexon overlap queries ──────────────────────────────

/// A monoexonic transcript record for interval tree queries
#[derive(Debug, Clone)]
pub struct MonoexonRecord {
    pub gene_id: GeneId,
    pub transcript_id: TranscriptId,
    pub chromosome: String,
    pub start: i64,
    pub end: i64,
    pub strand: Strand,
    pub start_vertex: VertexId,
    pub end_vertex: VertexId,
    pub exon_id: EdgeId,
}

/// Interval tree for overlap queries on monoexonic transcripts.
pub struct MonoexonIntervalTree {
    /// chromosome+strand -> interval tree
    trees: FxHashMap<String, COITree<MonoexonRecord, usize>>,
    /// Buffer for building trees
    intervals: FxHashMap<String, Vec<Interval<MonoexonRecord>>>,
    finalized: bool,
}

impl std::fmt::Debug for MonoexonIntervalTree {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("MonoexonIntervalTree")
            .field("finalized", &self.finalized)
            .finish()
    }
}

impl MonoexonIntervalTree {
    pub fn new() -> Self {
        Self {
            trees: FxHashMap::default(),
            intervals: FxHashMap::default(),
            finalized: false,
        }
    }

    pub fn insert(&mut self, record: MonoexonRecord) {
        let key = format!("{}_{}", record.chromosome, record.strand);
        let start = record.start.min(record.end) as i32;
        let end = record.start.max(record.end) as i32;
        self.intervals
            .entry(key)
            .or_insert_with(Vec::new)
            .push(Interval::new(start, end, record));
        self.finalized = false;
    }

    pub fn finalize(&mut self) {
        self.trees.clear();
        for (key, intervals) in self.intervals.iter() {
            self.trees.insert(key.clone(), COITree::new(intervals));
        }
        self.finalized = true;
    }

    pub fn query(
        &self,
        chromosome: &str,
        strand: Strand,
        start: i64,
        end: i64,
    ) -> Vec<MonoexonRecord> {
        let key = format!("{}_{}", chromosome, strand);
        let min_pos = start.min(end) as i32;
        let max_pos = start.max(end) as i32;

        match self.trees.get(&key) {
            Some(tree) => {
                let mut results = Vec::new();
                tree.query(min_pos, max_pos, |interval| {
                    results.push(interval.metadata.clone());
                });
                results
            }
            None => Vec::new(),
        }
    }
}

/// A gene record of the interval tree for overlap-based gene assignment
#[derive(Debug, Clone)]
pub struct GeneRecord {
    pub gene_id: GeneId,
    pub chromosome: String,
    pub start: i64,
    pub end: i64,
    pub strand: Strand,
}

pub struct GeneIntervalTree {
    trees: FxHashMap<String, COITree<GeneRecord, usize>>,
    intervals: FxHashMap<String, Vec<Interval<GeneRecord>>>,
    finalized: bool,
}

impl std::fmt::Debug for GeneIntervalTree {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("GeneIntervalTree")
            .field("finalized", &self.finalized)
            .finish()
    }
}

impl GeneIntervalTree {
    /// Every annotated gene, for building per-strand lookups.
    pub fn records(&self) -> impl Iterator<Item = &GeneRecord> {
        self.intervals.values().flatten().map(|iv| &iv.metadata)
    }

    pub fn new() -> Self {
        Self {
            trees: FxHashMap::default(),
            intervals: FxHashMap::default(),
            finalized: false,
        }
    }

    pub fn insert(&mut self, record: GeneRecord) {
        let key = record.chromosome.clone();
        let start = record.start.min(record.end) as i32;
        let end = record.start.max(record.end) as i32;
        self.intervals
            .entry(key)
            .or_insert_with(Vec::new)
            .push(Interval::new(start, end, record));
        self.finalized = false;
    }

    pub fn finalize(&mut self) {
        self.trees.clear();
        for (key, intervals) in self.intervals.iter() {
            self.trees.insert(key.clone(), COITree::new(intervals));
        }
        self.finalized = true;
    }

    pub fn query_overlap(
        &self,
        chromosome: &str,
        start: i64,
        end: i64,
    ) -> Vec<GeneRecord> {
        let min_pos = start.min(end) as i32;
        let max_pos = start.max(end) as i32;

        match self.trees.get(chromosome) {
            Some(tree) => {
                let mut results = Vec::new();
                tree.query(min_pos, max_pos, |interval| {
                    results.push(interval.metadata.clone());
                });
                results
            }
            None => Vec::new(),
        }
    }
}

// ─── Transcript interval tree for overlap queries ────────────────────────────

#[derive(Debug, Clone)]
pub struct TranscriptRecord {
    pub gene_id: GeneId,
    pub transcript_id: TranscriptId,
    pub chromosome: String,
    pub strand: Strand,
    pub min_pos: i64,
    pub max_pos: i64,
}

pub struct TranscriptIntervalTree {
    trees: FxHashMap<String, COITree<TranscriptRecord, usize>>,
    intervals: FxHashMap<String, Vec<Interval<TranscriptRecord>>>,
    finalized: bool,
}

impl std::fmt::Debug for TranscriptIntervalTree {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("TranscriptIntervalTree")
            .field("finalized", &self.finalized)
            .finish()
    }
}

impl TranscriptIntervalTree {
    /// Every annotated transcript, for building per-strand lookups.
    pub fn records(&self) -> impl Iterator<Item = &TranscriptRecord> {
        self.intervals.values().flatten().map(|iv| &iv.metadata)
    }

    pub fn new() -> Self {
        Self {
            trees: FxHashMap::default(),
            intervals: FxHashMap::default(),
            finalized: false,
        }
    }

    pub fn insert(&mut self, record: TranscriptRecord) {
        let key = record.chromosome.clone();
        let start = record.min_pos as i32;
        let end = record.max_pos as i32;
        self.intervals
            .entry(key)
            .or_insert_with(Vec::new)
            .push(Interval::new(start, end, record));
        self.finalized = false;
    }

    pub fn finalize(&mut self) {
        self.trees.clear();
        for (key, intervals) in self.intervals.iter() {
            self.trees.insert(key.clone(), COITree::new(intervals));
        }
        self.finalized = true;
    }

    pub fn query_overlap(
        &self,
        chromosome: &str,
        start: i64,
        end: i64,
    ) -> Vec<TranscriptRecord> {
        let min_pos = start.min(end) as i32;
        let max_pos = start.max(end) as i32;

        match self.trees.get(chromosome) {
            Some(tree) => {
                let mut results = Vec::new();
                tree.query(min_pos, max_pos, |interval| {
                    results.push(interval.metadata.clone());
                });
                results
            }
            None => Vec::new(),
        }
    }
}

// ─── Thread-local novel entity buffer ────────────────────────────────────────
// Each rayon thread collects novel entities here, then they're merged after.

#[derive(Debug, Default)]
pub struct NovelBuffer {
    pub locations: Vec<Location>,
    pub edges: Vec<Edge>,
    pub transcripts: Vec<(TranscriptPathKey, TranscriptModel)>,
    pub vertex_2_gene: Vec<(VertexId, GeneId, Strand)>,
    pub observed: Vec<ObservedRecord>,
    pub gene_annotations: Vec<NoveltyAnnotation>,
    pub transcript_annotations: Vec<NoveltyAnnotation>,
    pub exon_annotations: Vec<NoveltyAnnotation>,
    pub qc_records: Vec<QcMetrics>,
    /// One TSV row per chain-level assignment decision (ISOFORGE_CHAIN_TRACE).
    pub chain_traces: Vec<String>,
    /// The read being matched (its index in file order); `READ_IDX_NONE` outside the pass over the reads.
    pub read_idx: ReadIdx,
}

impl NovelBuffer {
    pub fn new() -> Self {
        Self { read_idx: READ_IDX_NONE, ..Self::default() }
    }

    pub fn clear(&mut self) {
        self.locations.clear();
        self.edges.clear();
        self.transcripts.clear();
        self.vertex_2_gene.clear();
        self.observed.clear();
        self.gene_annotations.clear();
        self.transcript_annotations.clear();
        self.exon_annotations.clear();
        self.qc_records.clear();
    }
}

// ─── Concurrent overlay dicts for novel entity sharing ───────────────────────
// Lookups check the read-only reference dict first, then the concurrent
// overlay; novel entities are inserted into the overlay so all threads see them.

/// Concurrent location dictionary for novel vertices: (chromosome, position)
/// → vertex_id, so threads share one ID per novel splice site.
pub type ConcurrentLocationDict = DashMap<(String, i64), VertexId>;

/// A novel edge in the concurrent dictionary: its ID, and the strand of the first read in file order that uses it
/// (`first`, that read's index). An edge is keyed by its vertices alone, so reads of both strands may use one (a
/// terminal exon of no length, or one whose end was matched past its splice site).
#[derive(Debug, Clone, Copy)]
pub struct NovelEdgeEntry {
    pub edge_id: EdgeId,
    pub first: ReadIdx,
    pub strand: Strand,
}

/// Concurrent edge dictionary for novel edges: EdgeKey → its ID and strand, so threads
/// share one ID per vertex pair.
pub type ConcurrentEdgeDict = DashMap<EdgeKey, NovelEdgeEntry>;

/// Identity of a novel gene locus, so reads at the same novel locus share
/// one gene.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum NovelGeneKey {
    /// Antisense to a known gene, keyed on that gene and the read's strand.
    Antisense(GeneId, Strand),
    /// Intergenic or fusion locus, keyed on (chromosome, strand, coordinate bin).
    Locus(String, Strand, i64),
}

/// A novel gene plus the union of the spans of every read assigned to it.
#[derive(Debug, Clone)]
pub struct NovelGeneEntry {
    pub gene_id: GeneId,
    pub chromosome: String,
    pub strand: Strand,
    pub start: i64,
    pub end: i64,
}

/// Concurrent registry of novel genes, the gene-level analogue of
/// `NovelTranscriptRegistry`.
pub type ConcurrentGeneDict = DashMap<NovelGeneKey, NovelGeneEntry>;

// ─── Shared reference collection ─────────────────────────────────────────────
// Shared read-only (Arc) across threads; novel vertices/edges go into the
// concurrent overlays.

/// Highest ID in the reference annotation for each entity space.
///
/// Novel entities are allocated above these, so an ID greater than the bound
/// that is not in the final catalogue belongs to an entity that was discarded.
#[derive(Debug, Clone, Copy, Default)]
pub struct RefIdBounds {
    pub gene: i64,
    pub transcript: i64,
    pub vertex: i64,
    pub edge: i64,
}

pub struct SharedRefs {
    pub location_dict: LocationDict,
    pub edge_dict: EdgeDict,
    pub transcript_dict: TranscriptDict,
    /// The annotated transcripts as hyperedges: their vertex chains and the containment index over them.
    pub hypergraph: crate::annotation_hypergraph::AnnotationHypergraph,
    pub vertex_2_gene: Vertex2Gene,
    pub gene_starts: GeneEndpointDict,
    pub gene_ends: GeneEndpointDict,
    pub gene_tree: GeneIntervalTree,
    pub transcript_tree: TranscriptIntervalTree,
    pub monoexon_tree: MonoexonIntervalTree,
    pub counters: Counters,
    /// Reference/novel ID boundary, captured before any novel entity is minted.
    pub ref_max: RefIdBounds,
    pub config: RunConfig,
    /// Concurrent overlay for novel vertices — shared across all threads.
    pub novel_location_dict: ConcurrentLocationDict,
    /// Concurrent overlay for novel edges — shared across all threads.
    pub novel_edge_dict: ConcurrentEdgeDict,
    /// Concurrent registry of novel genes — shared across all threads.
    pub novel_gene_dict: ConcurrentGeneDict,
    /// Optional set of short-read splice positions (STAR SJ.out.tab); positions
    /// in it are exempt from extended snapping.
    pub sr_splice_positions: Option<rustc_hash::FxHashSet<(String, i64)>>,
    /// Annotated donors and acceptors by strand, for the novel splice-site consensus.
    pub splice_sites: crate::novel_sites::AnnotatedSites,
    /// Annotated chains for approximate assignment; None when disabled.
    pub approx_chain: Option<crate::approximate_chain::ApproxChainIndex>,
    /// Annotated chains for consistency-based assignment; None when disabled.
    pub parent_chains: Option<crate::consistent_assignment::ParentIndex>,
    /// Annotated small RNAs counted, not assembled (small_rna.rs); None unless `small_rna` is set.
    pub small_rna: Option<Arc<crate::small_rna::SmallRnaIndex>>,
    /// The highest novel vertex ID splice-site matching may see (vertex_matching.rs): the engine sets it to the
    /// vertex counter before a window's reads are matched in parallel.
    pub visible_vertex_max: std::sync::atomic::AtomicI64,
}

/// Shared references behind Arc; mutable novel data is per thread.
pub struct ThreadContext {
    pub refs: Arc<RwLock<SharedRefs>>,
    pub novel: NovelBuffer,
}

impl ThreadContext {
    pub fn new(refs: Arc<RwLock<SharedRefs>>) -> Self {
        Self {
            refs,
            novel: NovelBuffer::new(),
        }
    }
}
