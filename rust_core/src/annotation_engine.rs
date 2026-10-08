// Annotation engine: loads the reference once, annotates the reads of each BAM
// in parallel batches, clusters TSS/TES end sites and writes Parquet outputs.

// HashMap and HashSet with deterministic hashing.
use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use rayon::prelude::*;

use crate::bam_reader::{BamReaderWithHeader, IntegratedBamReader, ParsedRead, ReaderPipeline};
use crate::cage_peaks::CagePeakIndex;
use crate::bam_reader::ClusteringInput;
use crate::data_structures::*;
use crate::errors::*;
use crate::genome_ref::{GenomeRef, SpliceJunctionDb};
use crate::parquet_io;
use crate::tss_tes_clustering::{self, TssTesEngine};
use crate::transcript_classifier::{self, NovelTranscriptRegistry};
use crate::types::*;

/// Main IsoForge engine.
pub struct IsoForgeEngine {
    pub refs: Arc<SharedRefs>,
    pub config: RunConfig,
    pub output_dir: PathBuf,
    /// Directory containing annotation Parquet files (for UTR regions).
    pub annot_dir: Option<PathBuf>,
    /// Optional genome FASTA path for integrated labeling + SJ correction.
    pub genome_path: Option<PathBuf>,
    /// Optional splice junction database paths (STAR SJ.out.tab files).
    pub sj_db: Option<Arc<SpliceJunctionDb>>,
    /// Shared concurrent registry for novel transcript deduplication.
    /// Ensures reads with identical exon structures share a single transcript ID.
    pub novel_tx_registry: NovelTranscriptRegistry,
    /// TSS/TES clustering configuration.
    pub cluster_config: Option<ClusterConfig>,
    /// Set for the run when the fraction of spliced reads with a poly(A) tail lies between
    /// `FULL_LENGTH_START_RATE` and `FULL_LENGTH_START_MAX_RATE`: starts are then placed by full-length reads.
    pub full_length_starts: std::sync::atomic::AtomicBool,
    /// Set for the run when the fraction of spliced reads with a poly(A) tail reaches `FULL_LENGTH_START_RATE`:
    /// an unspliced group with a start of its own then stays with its transcript.
    pub tails_present: std::sync::atomic::AtomicBool,
    /// Segments of chimeric reads from all datasets, for fusion calling.
    pub chimeric: std::sync::Mutex<Vec<crate::fusion::ChimericSegment>>,
    /// Each dataset's strand index, built on the first pass over its alignments and reused by later ones.
    strand_indexes: std::sync::Mutex<rustc_hash::FxHashMap<PathBuf, Arc<crate::strand_evidence::StrandIndex>>>,
}

/// Fraction of spliced reads carrying a poly(A) tail above which a run places starts
/// from full-length reads.
const FULL_LENGTH_START_RATE: f64 = 0.1;
/// Tailed reads a group needs before its own start becomes an alternative-start variant.
const MIN_FULL_LENGTH_START_READS: usize = 3;

/// Of the end variants made so far (cluster keys, largest first), the one a read group
/// missing one end cluster belongs to: the first sharing its other cluster, where that
/// cluster is not the dominant group's.
fn shared_end_variant(
    key: (Option<i64>, Option<i64>),
    dominant: (Option<i64>, Option<i64>),
    created: impl Iterator<Item = (Option<i64>, Option<i64>)>,
) -> Option<usize> {
    let mut created = created;
    match key {
        (None, Some(e)) if Some(e) != dominant.1 => created.position(|k| k.1 == Some(e)),
        (Some(t), None) if Some(t) != dominant.0 => created.position(|k| k.0 == Some(t)),
        _ => None,
    }
}
/// Fraction of spliced reads carrying a poly(A) tail up to which a run places starts
/// from full-length reads.
const FULL_LENGTH_START_MAX_RATE: f64 = 0.9;
/// Window around an annotated start within which a read's 5' end counts as reaching it,
/// for the run summary.
const ANNOTATED_START_WINDOW: i64 = 25;
/// Half-width of the window around an annotated transcript start or end that counts as
/// its promoter or poly(A) feature.
const ANNOTATED_END_FEATURE_WINDOW: i64 = 25;
/// Reads per dataset ending at a 3' cluster that confirm an annotated transcript no
/// read reaches the start of.
const THREE_PRIME_CONFIRM_READS: usize = 3;

/// Per-read state kept in memory after annotation; the full `ReadAnnotation`
/// is spilled to Parquet. Indexed by `ReadIdx`.
#[derive(Clone)]
pub struct ObsCore {
    pub gene_id: GeneId,
    pub transcript_id: TranscriptId,
    pub start_pos: i64,
    pub end_pos: i64,
    pub dataset: Arc<str>,
    pub location: Arc<str>,
    pub fraction_as: Option<f64>,
    /// Basecaller poly(A) tail length.
    pub polya_tail_len: Option<i32>,
    pub polya_tail_softclip: Option<i32>,
    pub strand: Strand,
    pub n_exons: i32,
    /// Largest snap distance of an internal splice site (see `ReadAnnotation`).
    pub max_splice_shift: i32,
    pub read_idx: ReadIdx,
    pub hint_flags: u8,
    pub is_internal_priming: bool,
    /// Untemplated G's at the 5' end: a capped start (see `cap_signature` in the reader).
    pub cap_signature: bool,
    /// The alignment ends before the transcript's 3' end (see `ParsedRead::unaligned_3p`).
    pub unaligned_3p: bool,
}

impl ClusteringInput for ObsCore {
    #[inline] fn read_idx(&self) -> ReadIdx { self.read_idx }
    #[inline] fn chromosome(&self) -> &str { &self.location }
    #[inline] fn strand_val(&self) -> Strand { self.strand }
    #[inline] fn pos_5p(&self) -> Option<i64> { Some(self.start_pos) }
    #[inline] fn pos_3p(&self) -> Option<i64> { Some(self.end_pos) }
    #[inline] fn dataset_name(&self) -> &str { &self.dataset }
    #[inline] fn polya_len(&self) -> Option<i64> { self.polya_tail_len.map(|v| v as i64) }
    #[inline] fn polya_softclip(&self) -> Option<i64> { self.polya_tail_softclip.map(|v| v as i64) }
    #[inline] fn is_ip(&self) -> bool { self.is_internal_priming }
    #[inline] fn frac_a(&self) -> Option<f64> { self.fraction_as }
    #[inline] fn n_exons(&self) -> i32 { self.n_exons }
    #[inline] fn unaligned_3p(&self) -> bool { self.unaligned_3p }
}

/// Per-chromosome map from every raw internal splice position of these reads
/// to (most-supported raw position within `window` bp, its read count).
/// Ties go to the position nearest the query, then the smaller coordinate.
pub fn junction_consensus(reads: &[ParsedRead], window: i64) -> rustc_hash::FxHashMap<String, rustc_hash::FxHashMap<i64, (i64, u32)>> {
    use rustc_hash::FxHashMap;
    let mut counts: FxHashMap<&str, FxHashMap<i64, u32>> = FxHashMap::default();
    for r in reads {
        let n = r.positions.len();
        if n < 4 { continue; }
        let c = counts.entry(r.chromosome.as_str()).or_default();
        for &p in &r.positions[1..n - 1] { *c.entry(p).or_default() += 1; }
    }
    let mut out: FxHashMap<String, FxHashMap<i64, (i64, u32)>> = FxHashMap::default();
    for (chrom, c) in counts {
        let mut v: Vec<(i64, u32)> = c.into_iter().collect();
        v.sort_unstable();
        let mut map = FxHashMap::default();
        let mut lo = 0usize;
        for i in 0..v.len() {
            let (p, n) = v[i];
            while v[lo].0 < p - window { lo += 1; }
            let mut best = (n, 0i64, p);
            for &(q, m) in &v[lo..] {
                if q > p + window { break; }
                let cand = (m, -(q - p).abs(), q);
                if cand.0 > best.0 || (cand.0 == best.0 && (cand.1 > best.1 || (cand.1 == best.1 && cand.2 < best.2))) { best = cand; }
            }
            map.insert(p, (best.2, best.0));
        }
        if !map.is_empty() { out.insert(chrom.to_string(), map); }
    }
    out
}

/// Everything the novel splice-site consensus needs beyond the reads and the annotated sites.
pub struct NovelSiteContext {
    pub genome: Option<crate::genome_ref::GenomeRef>,
    pub scores: Option<crate::splice_scores::SpliceScores>,
    pub params: crate::novel_sites::NovelSiteParams,
}

/// Results from processing all reads.
#[derive(Default)]
pub struct IsoForgeResults {
    pub observed: Vec<ReadAnnotation>,
    pub gene_annotations: Vec<NoveltyAnnotation>,
    pub transcript_annotations: Vec<NoveltyAnnotation>,
    pub exon_annotations: Vec<NoveltyAnnotation>,
    pub abundance: HashMap<(GeneId, TranscriptId, Arc<str>), i64>,
    pub novel_genes: Vec<Gene>,
    pub novel_transcripts: Vec<TranscriptModel>,
    pub novel_edges: Vec<Edge>,
    pub novel_locations: Vec<Location>,
    /// End variants of annotated transcripts: (variant, annotated parent, whether the
    /// variant's ends are compatible with the parent's annotated ends).
    pub end_variant_parents: Vec<(TranscriptId, TranscriptId, bool)>,
    pub qc_log: Vec<QcMetrics>,
    // Summary stats
    pub total_reads: usize,
    pub reads_processed: usize,
    pub reads_failed_qc: usize,
}

impl IsoForgeEngine {
    /// Create engine with pre-loaded references.
    pub fn new(refs: SharedRefs, config: RunConfig, output_dir: PathBuf) -> Self {
        Self {
            refs: Arc::new(refs),
            config,
            output_dir,
            annot_dir: None,
            genome_path: None,
            sj_db: None,
            novel_tx_registry: NovelTranscriptRegistry::new(),
            cluster_config: None,
            full_length_starts: std::sync::atomic::AtomicBool::new(false),
            tails_present: std::sync::atomic::AtomicBool::new(false),
            chimeric: std::sync::Mutex::new(Vec::new()),
            strand_indexes: std::sync::Mutex::new(rustc_hash::FxHashMap::default()),
        }
    }

    /// `<prefix>_<name>` next to the run's outputs (the database directory is `<prefix>_isoforge_db`).
    fn beside_outputs(&self, name: &str) -> PathBuf {
        let db_name = self.output_dir.file_name().and_then(|n| n.to_str()).unwrap_or("");
        let prefix = db_name.strip_suffix("_isoforge_db").unwrap_or(db_name);
        self.output_dir.parent().unwrap_or(&self.output_dir).join(format!("{prefix}_{name}"))
    }

    /// Set genome reference path for integrated processing.
    pub fn with_genome(mut self, genome_path: PathBuf) -> Self {
        self.genome_path = Some(genome_path);
        self
    }

    /// Set annotation directory (for UTR regions used by clustering).
    pub fn with_annot_dir(mut self, annot_dir: PathBuf) -> Self {
        self.annot_dir = Some(annot_dir);
        self
    }

    /// Enable TSS/TES clustering with the given configuration.
    pub fn with_clustering(mut self, cluster_config: ClusterConfig) -> Self {
        self.cluster_config = Some(cluster_config);
        self
    }

    /// Set splice junction database from STAR SJ.out.tab files.
    ///
    /// When `max_extended_snap > max_splice_dist`, also collects every splice
    /// position (±2 bp) so vertex matching can gate extended snapping.
    pub fn with_sj_db(mut self, sj_paths: &[PathBuf]) -> IsoForgeResult<Self> {
        if !sj_paths.is_empty() {
            let db = SpliceJunctionDb::from_star_sj_tabs(sj_paths)?;
            log::info!("[IsoForge] Loaded {} known splice junctions from {} files",
                     db.len(), sj_paths.len());

            // Short-read splice positions for extended-snap gating
            if self.config.max_extended_snap > self.config.max_splice_dist {
                let positions = db.splice_positions(2); // ±2 bp fuzzy
                log::info!(
                    "[IsoForge] Built SR splice position index: {} positions (±2 bp fuzzy) for Layer 3 gating",
                    positions.len()
                );
                if let Some(refs) = Arc::get_mut(&mut self.refs) {
                    refs.sr_splice_positions = Some(positions);
                }
            }

            self.sj_db = Some(Arc::new(db));
        }
        Ok(self)
    }

    /// The splice-site scores for the novel splice-site consensus; None when
    /// `novel_site_window` is 0 or no scores are given.
    fn novel_site_context(&self) -> IsoForgeResult<Option<NovelSiteContext>> {
        let Some(path) = self.config.splice_scores.as_ref().filter(|_| self.config.novel_site_window > 0) else {
            return Ok(None);
        };
        let t = std::time::Instant::now();
        let s = crate::splice_scores::SpliceScores::load(Path::new(path))?;
        log::info!("[IsoForge] Splice-site scores: {} sites from {} in {:.1} s", s.len(), path, t.elapsed().as_secs_f64());
        let scores = Some(s);
        let params = crate::novel_sites::NovelSiteParams {
            window: self.config.novel_site_window,
            annotated_zone: self.config.max_extended_snap,
            snap_zone: self.config.max_splice_dist,
            ..crate::novel_sites::NovelSiteParams::default()
        };
        let genome = match &self.genome_path {
            Some(path) => Some(crate::genome_ref::GenomeRef::open(path)?),
            None => None,
        };
        Ok(Some(NovelSiteContext { genome, scores, params }))
    }

    /// Load references from Parquet files in the given directory.
    pub fn from_parquet(annot_dir: &Path, config: RunConfig, output_dir: PathBuf) -> IsoForgeResult<Self> {
        let mut config = config;
        config.chain_trace = std::env::var_os("ISOFORGE_CHAIN_TRACE").is_some();
        let location_dict = parquet_io::load_location_dict(&annot_dir.join("locations.parquet"))?;
        let edge_dict = parquet_io::load_edge_dict(&annot_dir.join("edges.parquet"))?;
        let transcript_dict = parquet_io::load_transcript_dict(&annot_dir.join("transcripts.parquet"))?;
        let vertex_2_gene = parquet_io::load_vertex_2_gene(&annot_dir.join("vertex_2_gene.parquet"))?;
        let gene_starts = parquet_io::load_gene_endpoints(&annot_dir.join("gene_starts.parquet"))?;
        let gene_ends = parquet_io::load_gene_endpoints(&annot_dir.join("gene_ends.parquet"))?;

        // The annotation as a hypergraph: one hyperedge (vertex chain) per transcript, with the containment index
        // FSM/ISM matching uses; chains from the index's hyperedges.parquet when it has one
        let hypergraph = crate::annotation_hypergraph::AnnotationHypergraph::load(annot_dir, &transcript_dict, &location_dict, &edge_dict)?;
        eprintln!("[IsoForge] Annotation hypergraph: {} hyperedges ({} chains from hyperedges.parquet)",
                  hypergraph.len(), hypergraph.from_file);

        // Load genes and build interval trees
        let genes_path = annot_dir.join("genes.parquet");
        let mut gene_tree = GeneIntervalTree::new();
        let mut transcript_tree = TranscriptIntervalTree::new();
        let mut monoexon_tree = MonoexonIntervalTree::new();

        // a transcript's strand is its gene's
        let mut gene_strand: rustc_hash::FxHashMap<GeneId, Strand> = rustc_hash::FxHashMap::default();
        if genes_path.exists() {
            let genes = parquet_io::load_genes(&genes_path)?;
            for gene in &genes {
                gene_strand.insert(gene.gene_id, gene.strand);
                gene_tree.insert(GeneRecord {
                    gene_id: gene.gene_id,
                    chromosome: gene.chromosome.clone(),
                    start: gene.start,
                    end: gene.end,
                    strand: gene.strand,
                });
            }
            gene_tree.finalize();
            log::info!("[IsoForge] Built gene interval tree ({} genes)", genes.len());
        }

        // Build transcript interval tree and monoexon tree from transcript dict
        let mut n_multi = 0usize;
        let mut n_mono = 0usize;
        for model in transcript_dict.values() {
            let strand = gene_strand.get(&model.gene_id).copied()
                .unwrap_or_else(|| strand_at_start(&vertex_2_gene, model.start_vertex, model.gene_id));

            let min_pos = model.start_pos.min(model.end_pos);
            let max_pos = model.start_pos.max(model.end_pos);

            transcript_tree.insert(TranscriptRecord {
                gene_id: model.gene_id,
                transcript_id: model.transcript_id,
                chromosome: model.chromosome.clone(),
                strand,
                min_pos,
                max_pos,
            });

            if model.n_exons == 1 {
                monoexon_tree.insert(MonoexonRecord {
                    gene_id: model.gene_id,
                    transcript_id: model.transcript_id,
                    chromosome: model.chromosome.clone(),
                    strand,
                    start: min_pos,
                    end: max_pos,
                    start_vertex: model.start_vertex,
                    end_vertex: model.end_vertex,
                    exon_id: model.start_exon,
                });
                n_mono += 1;
            } else {
                n_multi += 1;
            }
        }
        transcript_tree.finalize();
        monoexon_tree.finalize();
        eprintln!(
            "[IsoForge] Built transcript trees ({} multi-exon, {} mono-exon)",
            n_multi, n_mono
        );

        // Initialize counters — prefer counters.json if present, else derive from max IDs
        let counters_path = annot_dir.join("counters.json");
        let counters = if counters_path.exists() {
            let (g, t, v, e, o, d) = parquet_io::load_counters(&counters_path)?;
            eprintln!(
                "[IsoForge] Loaded counters: genes={}, transcripts={}, vertices={}, edges={}",
                g, t, v, e
            );
            Counters::new(g, t, v, e, o, d)
        } else {
            let max_gene = transcript_dict.values().map(|t| t.gene_id).max().unwrap_or(0);
            let max_transcript = transcript_dict.values().map(|t| t.transcript_id).max().unwrap_or(0);
            let max_vertex = location_dict.chromosomes.values()
                .flat_map(|c| c.positions.iter().map(|(_, id)| *id))
                .max()
                .unwrap_or(0);
            let max_edge = edge_dict.values().map(|e| e.edge_id).max().unwrap_or(0);
            Counters::new(max_gene, max_transcript, max_vertex, max_edge, 0, 0)
        };

        // Snapshot before annotation starts minting novel entities above these.
        let ref_max = RefIdBounds {
            gene: counters.gene.get(),
            transcript: counters.transcript.get(),
            vertex: counters.vertex.get(),
            edge: counters.edge.get(),
        };

        let refs = SharedRefs {
            location_dict,
            edge_dict,
            transcript_dict,
            hypergraph,
            vertex_2_gene,
            gene_starts,
            gene_ends,
            gene_tree,
            transcript_tree,
            monoexon_tree,
            counters,
            ref_max,
            config: config.clone(),
            novel_location_dict: ConcurrentLocationDict::new(),
            novel_edge_dict: ConcurrentEdgeDict::new(),
            novel_gene_dict: ConcurrentGeneDict::new(),
            sr_splice_positions: None,
            splice_sites: crate::novel_sites::AnnotatedSites::default(),
            approx_chain: None,
            parent_chains: None,
            small_rna: None,
            visible_vertex_max: std::sync::atomic::AtomicI64::new(crate::vertex_matching::ALL_VISIBLE),
        };
        let refs = SharedRefs { splice_sites: crate::novel_sites::AnnotatedSites::from_reference(&refs), ..refs };
        let refs = if config.approx_chain_tolerance > 0 {
            let index = crate::approximate_chain::ApproxChainIndex::from_reference(&refs);
            eprintln!("[IsoForge] Approximate chain assignment enabled (tolerance {} bp, {}, annotated sites {})",
                      config.approx_chain_tolerance,
                      if config.chain_first { "chain-first" } else { "after NIC" },
                      if config.approx_chain_lock_annotated { "locked" } else { "unlocked" });
            SharedRefs { approx_chain: Some(index), ..refs }
        } else {
            refs
        };

        // The parent index: every annotated multi-exon chain.
        let refs = {
            let index = crate::consistent_assignment::ParentIndex::from_reference(&refs);
            eprintln!(
                "[IsoForge] Consistency-based assignment: {} annotated multi-exon chains indexed \
                 (unique-parent ISM {:?}, nested FSM ambiguity {})",
                index.len(),
                config.assign_unique_ism,
                if config.nested_fsm_ambiguous { "on" } else { "off" }
            );
            SharedRefs { parent_chains: Some(index), ..refs }
        };

        // Annotated small RNAs counted, not assembled (small_rna.rs)
        let refs = if config.small_rna {
            let gene_types = parquet_io::load_gene_types(&annot_dir.join("gene_annotations.parquet"))?;
            let index = crate::small_rna::SmallRnaIndex::from_reference(&refs, &gene_strand, &gene_types, &config);
            eprintln!("[IsoForge] Small-RNA rule: {} annotated one-exon transcripts of {} (up to {} nt), {} of them inside \
                       an exon of another gene; their reads are kept from {} bases at any mapping quality, count for \
                       them and build no model",
                      index.len(), config.small_rna_types.join(", "), config.small_rna_max_length, index.n_nested(),
                      config.small_rna_min_read_length);
            SharedRefs { small_rna: Some(Arc::new(index)), ..refs }
        } else {
            refs
        };

        Ok(Self::new(refs, config, output_dir))
    }

    /// Process reads in parallel using rayon.
    fn process_reads_parallel(
        &self,
        reads: &[ParsedRead],
        dataset_name: &Arc<str>,
        batch_size: usize,
        novel_ctx: Option<&NovelSiteContext>,
    ) -> IsoForgeResult<IsoForgeResults> {
        let n_threads = rayon::current_num_threads();
        log::info!("[IsoForge] Processing with {} threads, batch_size={}", n_threads, batch_size);

        // Partition reads into batches
        let batches: Vec<&[ParsedRead]> = reads.chunks(batch_size).collect();

        // Junction consensus for this window: maps each raw internal splice position to the most-supported raw
        // position within the consensus window, per chromosome.
        let consensus: Option<rustc_hash::FxHashMap<String, rustc_hash::FxHashMap<i64, (i64, u32)>>> =
            if self.config.junction_consensus_window > 0 {
                Some(junction_consensus(reads, self.config.junction_consensus_window))
            } else { None };

        // Novel splice positions join their consensus site.
        let novel_sites = novel_ctx.map(|c| crate::novel_sites::novel_site_consensus(reads, &self.refs.splice_sites, c.genome.as_ref(), c.scores.as_ref(), &c.params));

        // Creates the window's novel splice vertices before any read is matched: the internal sites that no earlier
        // vertex takes, in coordinate order. The pass below sees these and the earlier vertices only.
        let mut pre_buffer = NovelBuffer::new();
        {
            use std::sync::atomic::Ordering;
            let before = self.refs.counters.vertex.get();
            self.refs.visible_vertex_max.store(before, Ordering::Relaxed);
            let mut proposed: Vec<(&str, i64, bool)> = batches
                .par_iter()
                .map(|batch| {
                    let mut out: Vec<(&str, i64, bool)> = Vec::new();
                    let mut sites: Vec<(i64, bool)> = Vec::new();
                    for read in batch.iter() {
                        if read.positions.len() < 4 { continue; }
                        let chrom = read.chromosome.as_str();
                        let cons = consensus.as_ref().and_then(|c| c.get(chrom));
                        let (positions, settled) = Self::read_sites(read, novel_sites.as_ref().and_then(|m| m.get(chrom)), &self.config);
                        sites.clear();
                        crate::vertex_matching::proposed_splice_vertices(
                            &self.refs.location_dict, chrom, &positions, &self.refs.config, &self.refs.novel_location_dict,
                            self.refs.sr_splice_positions.as_ref(), cons, settled, before, &mut sites);
                        out.extend(sites.iter().map(|&(p, forced)| (chrom, p, forced)));
                    }
                    out
                })
                .flatten()
                .collect();
            // without a window consensus, raw positions within the snapping distance of one already taken join it
            let thin = consensus.is_none().then_some(self.config.max_splice_dist);
            crate::vertex_matching::create_window_vertices(&mut proposed, thin, &self.refs.counters, &self.refs.config,
                                                           &mut pre_buffer, &self.refs.novel_location_dict);
            self.refs.visible_vertex_max.store(self.refs.counters.vertex.get(), Ordering::Relaxed);
        }

        // Process each batch in parallel
        let batch_results: Vec<BatchResult> = batches
            .par_iter()
            .map(|batch| {
                let mut novel_buffer = NovelBuffer::new();
                let mut annotations = Vec::new();
                let mut qc_log = Vec::new();

                for read in batch.iter() {
                    let cons = consensus.as_ref().and_then(|c| c.get(read.chromosome.as_str()));
                    let sites = novel_sites.as_ref().and_then(|m| m.get(read.chromosome.as_str()));
                    novel_buffer.read_idx = read.read_idx;
                    match self.annotate_read(read, dataset_name, &mut novel_buffer, cons, sites) {
                        Ok(annot) => annotations.push(annot),
                        Err(IsoForgeError::QcError(_msg)) => {
                            qc_log.push(QcMetrics {
                                dataset: dataset_name.to_string(),
                                read_id: read.read_id.clone(),
                                passed_qc: false,
                                primary_mapped: true,
                                read_length: read.read_length,
                                fraction_aligned: None,
                                identity: None,
                            });
                        }
                        Err(_e) => {
                            qc_log.push(QcMetrics {
                                dataset: dataset_name.to_string(),
                                read_id: read.read_id.clone(),
                                passed_qc: false,
                                primary_mapped: true,
                                read_length: read.read_length,
                                fraction_aligned: None,
                                identity: None,
                            });
                        }
                    }
                }

                BatchResult {
                    annotations,
                    novel_buffer,
                    qc_log,
                }
            })
            .collect();

        // Merge results
        let mut results = self.merge_results(batch_results, dataset_name)?;
        results.novel_locations.extend(pre_buffer.locations);
        Ok(results)
    }

    /// A read's exon boundaries as annotation matches them (novel sites moved onto their consensus) and which of
    /// its internal sites are settled where they lie.
    fn read_sites<'a>(read: &'a ParsedRead, novel_sites: Option<&crate::novel_sites::SiteMap>, config: &RunConfig)
                      -> (std::borrow::Cow<'a, [i64]>, u64) {
        // sites the read's own bases favour over nearby annotated ones stay where they lie
        let keep = if config.sequence_settled_sites { read.sequence_settled } else { 0 };
        let (moved, settled) = novel_sites.map_or((None, 0), |m| crate::novel_sites::apply_keeping(&read.positions, read.strand, m, keep));
        let positions = match moved {
            Some(v) => std::borrow::Cow::Owned(v),
            None => std::borrow::Cow::Borrowed(read.positions.as_slice()),
        };
        (positions, settled | keep)
    }

    /// Annotate a single read.
    fn annotate_read(
        &self,
        read: &ParsedRead,
        dataset_name: &Arc<str>,
        novel_buffer: &mut NovelBuffer,
        consensus: Option<&rustc_hash::FxHashMap<i64, (i64, u32)>>,
        novel_sites: Option<&crate::novel_sites::SiteMap>,
    ) -> IsoForgeResult<ReadAnnotation> {
        // Reads that fail the coverage or identity checks are filtered in the reader.

        let strand = read.strand;
        let chrom = &read.chromosome;

        // Multi-exon vs mono-exon dispatch
        // positions are already computed in ParsedRead
        let (positions, settled) = Self::read_sites(read, novel_sites, &self.config);
        let positions: &[i64] = &positions;
        let mut annot = transcript_classifier::identify_transcript(
            chrom,
            positions,
            strand,
            &self.refs,
            novel_buffer,
            &self.novel_tx_registry,
            consensus,
            read.ambiguous_sites,
            settled,
            read.small_rna,
        );

        // Fill in read-level metadata
        annot.read_id = read.read_id.clone();
        annot.read_idx = read.read_idx;
        annot.dataset = Arc::clone(dataset_name);
        annot.location = chrom.clone();
        annot.strand = strand;
        annot.read_length = read.read_length;
        annot.fraction_as = read.fraction_as;
        annot.polya_tail_len = read.polya_tail_len;
        annot.cap_signature = read.cap_signature;
        annot.unaligned_3p = read.unaligned_3p;
        annot.polya_tail_softclip = read.polya_tail_softclip;
        annot.internal_priming = read.internal_priming.clone();

        Ok(annot)
    }

    /// Append chain-level decision rows to `<db>/chain_assignment.tsv`,
    /// writing the header when the file is new.
    fn append_chain_traces(&self, rows: &[String]) -> IsoForgeResult<()> {
        use std::io::Write;
        let path = self.output_dir.join("chain_assignment.tsv");
        let new = !path.exists();
        let mut f = std::fs::OpenOptions::new().create(true).append(true).open(&path)?;
        if new {
            writeln!(f, "{}", crate::approximate_chain::CHAIN_TRACE_HEADER)?;
        }
        for r in rows {
            writeln!(f, "{r}")?;
        }
        Ok(())
    }

    fn merge_results(
        &self,
        batch_results: Vec<BatchResult>,
        dataset_name: &Arc<str>,
    ) -> IsoForgeResult<IsoForgeResults> {
        let mut results = IsoForgeResults::default();
        let mut seen_novel_tx_keys = rustc_hash::FxHashSet::default();

        for br in batch_results {
            results.reads_processed += br.annotations.len();
            results.reads_failed_qc += br.qc_log.len();

            for mut annot in br.annotations {
                // Update abundance
                let key = (annot.gene_id, annot.transcript_id, Arc::clone(dataset_name));
                *results.abundance.entry(key).or_insert(0) += 1;

                // Moves the read's novelty annotations into the results.
                results.gene_annotations.extend(std::mem::take(&mut annot.gene_novelty));
                results.transcript_annotations.extend(std::mem::take(&mut annot.transcript_novelty));
                results.exon_annotations.extend(std::mem::take(&mut annot.exon_novelty));

                results.observed.push(annot);
            }

            results.novel_edges.extend(br.novel_buffer.edges);
            results.novel_locations.extend(br.novel_buffer.locations);
            if !br.novel_buffer.chain_traces.is_empty() {
                self.append_chain_traces(&br.novel_buffer.chain_traces)?;
            }
            for (key, model) in br.novel_buffer.transcripts {
                // Only keep the first occurrence of each path key
                if !seen_novel_tx_keys.contains(&key) {
                    seen_novel_tx_keys.insert(key);
                    results.novel_transcripts.push(model);
                }
            }

            results.qc_log.extend(br.qc_log);
        }

        // A novel edge takes the strand the dictionary kept for it: that of the first read in file order that uses it.
        for e in results.novel_edges.iter_mut() {
            if let Some(entry) = self.refs.novel_edge_dict.get(&EdgeKey { v1: e.v1, v2: e.v2, edge_type: e.edge_type }) {
                e.strand = entry.strand;
            }
        }

        // A novel model takes its ends and gene from the first read of the window that carries it, in file order.
        if !results.novel_transcripts.is_empty() {
            let mut first: rustc_hash::FxHashMap<TranscriptId, usize> =
                results.novel_transcripts.iter().map(|m| (m.transcript_id, usize::MAX)).collect();
            for (i, a) in results.observed.iter().enumerate() {
                if let Some(slot) = first.get_mut(&a.transcript_id) {
                    if *slot == usize::MAX { *slot = i; }
                }
            }
            for m in results.novel_transcripts.iter_mut() {
                let i = first[&m.transcript_id];
                if i == usize::MAX { continue; }
                let a = &results.observed[i];
                m.gene_id = a.gene_id;
                m.start_exon = a.start_exon;
                m.end_exon = a.end_exon;
                m.start_vertex = a.start_vertex;
                m.end_vertex = a.end_vertex;
                m.start_pos = a.start_pos;
                m.end_pos = a.end_pos;
            }
        }

        results.total_reads = results.reads_processed + results.reads_failed_qc;

        eprintln!(
            "[IsoForge] Processed {}/{} reads ({} failed QC, {} novel genes, {} novel transcripts)",
            results.reads_processed,
            results.total_reads,
            results.reads_failed_qc,
            self.refs.novel_gene_dict.len(),
            results.novel_transcripts.len(),
        );

        if std::env::var_os("ISOFORGE_DEBUG_MEM").is_some() {
            let (antisense, loci) = self.refs.novel_gene_dict.iter().fold(
                (0usize, 0usize),
                |(a, l), e| match e.key() {
                    NovelGeneKey::Antisense(..) => (a + 1, l),
                    NovelGeneKey::Locus(..) => (a, l + 1),
                },
            );
            eprintln!(
                "[IsoForge][mem] vertices={} edges={} transcripts={} genes={} \
                 (antisense={} loci={})",
                self.refs.novel_location_dict.len(),
                self.refs.novel_edge_dict.len(),
                self.novel_tx_registry.len(),
                self.refs.novel_gene_dict.len(),
                antisense,
                loci,
            );
        }

        Ok(results)
    }

    /// Run full pipeline: process BAM files, write output.
    ///
    ///   1. Read all BAMs → Vec<ParsedRead> per dataset
    ///   2. Annotate all reads (transcript identification)
    ///   3. If clustering enabled: TSS/TES clustering on all reads
    ///   4. Write all outputs (observed, annotations, abundance, clusters)
    pub fn run(
        &self,
        bam_configs: &[(PathBuf, String)], // (bam_path, dataset_name)
    ) -> IsoForgeResult<()> {
        // ── Phase 1+2: per-BAM read → annotate → compact ────────
        // Each BAM is read and annotated in windows and reduced to ObsCore; full
        // ParsedReads are dropped and novelty annotations are streamed to Parquet.

        // Ensure output directory exists before opening streaming writer
        std::fs::create_dir_all(&self.output_dir)?;
        let annot_path = self.output_dir.join("annotations.parquet");
        let mut annot_writer = parquet_io::StreamingAnnotationWriter::new(&annot_path)?;

        let mut all_results = IsoForgeResults::default();
        // Per-read state kept after annotation; full ReadAnnotations are spilled to disk.
        let mut cores: Vec<ObsCore> = Vec::new();
        // (core index, unique annotated parent) for truncated reads.
        let mut ism_parents: Vec<(usize, TranscriptId)> = Vec::new();
        // Parallel to `cores`: whether the read survives the full-length filter.
        let mut keep: Vec<bool> = Vec::new();
        let mut chrom_cache: rustc_hash::FxHashMap<String, Arc<str>> =
            rustc_hash::FxHashMap::default();
        let spill_path = self.output_dir.join("observed.spill.parquet");
        let mut obs_writer = parquet_io::StreamingObservedWriter::new(&spill_path)?;
        // Reads kept for quantification only (MAPQ in [quant_min_mapq, min_mapq)).
        let mut quant_writer = if self.config.quant_min_mapq < self.config.min_mapq {
            Some(parquet_io::QuantReadWriter::new(&self.output_dir.join("quant_reads.parquet"))?)
        } else {
            None
        };
        // Distinct (entity, attribute, value) triples already written.
        let mut seen_annotations: rustc_hash::FxHashSet<(i64, Arc<str>, Arc<str>)> =
            rustc_hash::FxHashSet::default();
        // Half-open range of `cores` belonging to each dataset.
        let mut dataset_ranges: Vec<(String, std::ops::Range<usize>)> =
            Vec::with_capacity(bam_configs.len());
        let mut total_reads_count: usize = 0;

        // Capacity is reserved from the BAM index read count.
        let expected_reads = if self.config.expected_reads > 0 {
            self.config.expected_reads
        } else {
            bam_configs.iter().map(|(p, _)| count_indexed_reads(p)).sum()
        };
        if expected_reads > 0 {
            cores.reserve_exact(expected_reads);
            eprintln!("[IsoForge] Preallocating for {} reads", expected_reads);
        }

        let novel_ctx = self.novel_site_context()?;

        // Junction groups (junction_groups.rs): every dataset's introns counted, then the reads of contested groups
        // scored and each group's mixture fitted, before any read is annotated.
        let junctions = match (&self.genome_path, self.config.junction_hmm) {
            (Some(genome_path), true) => {
                let t0 = std::time::Instant::now();
                let mut jg = crate::junction_groups::JunctionGroups::new(
                    &self.refs, self.config.junction_hmm_window, self.config.junction_hmm_min_support);
                if let Some(path) = &self.config.junction_hmm_model {
                    let text = std::fs::read_to_string(path)
                        .map_err(|e| IsoForgeError::ConfigError(format!("--junction_hmm_model {path}: {e}")))?;
                    let doc: serde_json::Value = serde_json::from_str(&text)
                        .map_err(|e| IsoForgeError::ConfigError(format!("--junction_hmm_model {path}: {e}")))?;
                    let params: crate::junction_hmm::HmmParams = serde_json::from_value(doc.get("params").cloned().unwrap_or(doc))
                        .map_err(|e| IsoForgeError::ConfigError(format!("--junction_hmm_model {path}: {e}")))?;
                    jg = jg.with_params(params);
                }
                let jg = Arc::new(jg);
                let genome = crate::genome_ref::GenomeRef::open(genome_path.clone())?;
                for (bam_path, dataset_name) in bam_configs {
                    self.read_bam_windows(bam_path, dataset_name, self.config.read_window, Some(jg.clone()), |_| Ok(()))?;
                }
                let t_count = t0.elapsed().as_secs_f64();
                jg.build_groups(&genome);
                for (bam_path, dataset_name) in bam_configs {
                    self.read_bam_windows(bam_path, dataset_name, self.config.read_window, Some(jg.clone()), |_| Ok(()))?;
                }
                let t_score = t0.elapsed().as_secs_f64();
                jg.fit();
                eprintln!("[IsoForge] Junction groups ready in {:.1}s (counting {:.1}s, scoring {:.1}s, fitting {:.1}s)",
                          t0.elapsed().as_secs_f64(), t_count, t_score - t_count, t0.elapsed().as_secs_f64() - t_score);
                // the strand counters are reset after the pre-passes
                for c in [&crate::bam_reader::STRAND_FROM_MOTIF, &crate::bam_reader::STRAND_MOTIF_FLIPPED,
                          &crate::bam_reader::STRAND_INFERRED, &crate::bam_reader::STRAND_FLIPPED] {
                    c.store(0, std::sync::atomic::Ordering::Relaxed);
                }
                Some(jg)
            }
            _ => None,
        };

        for (bam_path, dataset_name) in bam_configs {
            log::info!("[IsoForge] Reading dataset: {} ({:?})", dataset_name, bam_path);

            let batch_size = self.config.batch_size;
            let window_size = self.config.read_window;
            let ds_arc: Arc<str> = Arc::from(dataset_name.as_str());
            let ds_start = cores.len();
            let mut ds_reads: usize = 0;

            self.read_bam_windows(bam_path, dataset_name, window_size, junctions.clone(), |mut reads| {
                // Quantification-only reads are written out and leave the window before the reads are indexed.
                if let Some(w) = quant_writer.as_mut() {
                    if reads.iter().any(|r| r.quant_only) {
                        let (quant, rest): (Vec<ParsedRead>, Vec<ParsedRead>) = reads.into_iter().partition(|r| r.quant_only);
                        w.flush(&quant)?;
                        reads = rest;
                    }
                }
                // Dense indices continue across windows and across datasets.
                for (i, read) in reads.iter_mut().enumerate() {
                    read.read_idx = (total_reads_count + i) as ReadIdx;
                }

                let mut results = self.process_reads_parallel(&reads, &ds_arc, batch_size, novel_ctx.as_ref())?;

                // The full annotations go to disk; only ObsCore stays resident.
                obs_writer.flush(&results.observed)?;
                cores.reserve(results.observed.len());
                for annot in results.observed.iter() {
                    let location = match chrom_cache.get(annot.location.as_str()) {
                        Some(c) => Arc::clone(c),
                        None => {
                            let c: Arc<str> = Arc::from(annot.location.as_str());
                            chrom_cache.insert(annot.location.clone(), Arc::clone(&c));
                            c
                        }
                    };
                    if let Some(parent) = annot.ism_parent {
                        ism_parents.push((cores.len(), parent));
                    }
                    cores.push(ObsCore {
                        gene_id: annot.gene_id,
                        transcript_id: annot.transcript_id,
                        start_pos: annot.start_pos,
                        end_pos: annot.end_pos,
                        dataset: Arc::clone(&annot.dataset),
                        location,
                        fraction_as: annot.fraction_as,
                        polya_tail_len: annot.polya_tail_len.map(|v| v as i32),
                        polya_tail_softclip: annot.polya_tail_softclip.map(|v| v as i32),
                        read_idx: annot.read_idx,
                        strand: annot.strand,
                        n_exons: annot.n_exons,
                        max_splice_shift: annot.max_splice_shift,
                        hint_flags: annot.hint_flags,
                        is_internal_priming: annot
                            .internal_priming
                            .as_ref()
                            .map_or(false, |ip| ip.is_internal_priming),
                        cap_signature: annot.cap_signature,
                        unaligned_3p: annot.unaligned_3p,
                    });
                }
                results.observed.clear();
                results.observed.shrink_to_fit();

                total_reads_count += reads.len();
                ds_reads += reads.len();

                // Full ParsedReads die here; ObsCore carries what survives.
                drop(reads);

                merge_into(&mut all_results, results);

                // Novelty annotations are streamed to disk, keeping only the first
                // occurrence of each (entity, attribute, value).
                let mut batch_annotations: Vec<NoveltyAnnotation> = Vec::new();
                batch_annotations.append(&mut all_results.gene_annotations);
                batch_annotations.append(&mut all_results.transcript_annotations);
                batch_annotations.append(&mut all_results.exon_annotations);
                batch_annotations.retain(|a| {
                    seen_annotations.insert((a.id, Arc::clone(&a.attribute), Arc::clone(&a.value)))
                });
                annot_writer.flush(&batch_annotations)?;

                if std::env::var_os("ISOFORGE_DEBUG_MEM").is_some() {
                    debug_heap_report(
                        &cores,
                        &all_results,
                        seen_annotations.len(),
                        &self.refs,
                    );
                }

                Ok(())
            })?;

            eprintln!("[IsoForge] Loaded {} reads from {:?}", ds_reads, bam_path);
            dataset_ranges.push((dataset_name.clone(), ds_start..cores.len()));

            eprintln!(
                "[IsoForge] Memory: processed + compacted dataset {}, annotations flushed to disk",
                dataset_name
            );
        }
        if let Some(jg) = &junctions {
            jg.log_decisions();
            let (model, table) = (self.beside_outputs("junction_hmm.json"), self.beside_outputs("junction_groups.tsv"));
            match jg.write_outputs(&model, &table) {
                Ok(()) => eprintln!("[IsoForge] Junction groups: pair-HMM written to {}, groups to {}", model.display(), table.display()),
                Err(e) => eprintln!("[IsoForge] Junction groups: outputs not written: {e}"),
            }
        }

        if self.config.detect_fusions {
            let segs = std::mem::take(&mut *self.chimeric.lock().unwrap());
            let n = crate::fusion::write_segments(&self.output_dir.join("chimeric_segments.parquet"), &segs)?;
            eprintln!("[IsoForge] {} alignment segments of chimeric reads kept for fusion calling", n);
        }
        if let Some(w) = quant_writer {
            let n = w.close()?;
            eprintln!("[IsoForge] {} reads with MAPQ {}-{} kept for quantification only", n, self.config.quant_min_mapq, self.config.min_mapq.saturating_sub(1));
        }
        // Finalize streaming annotation writer
        let n_annotations = annot_writer.close()?;
        // Finalize the spill of observed reads
        let n_spilled = obs_writer.close()?;
        log::info!("[IsoForge] Spilled {} observed reads to {:?}", n_spilled, spill_path);

        // Report how often the poly(A) tail reversed the alignment strand.
        if self.config.infer_strand_from_polya {
            use std::sync::atomic::Ordering::Relaxed;
            let motif = crate::bam_reader::STRAND_FROM_MOTIF.load(Relaxed);
            let motif_flipped = crate::bam_reader::STRAND_MOTIF_FLIPPED.load(Relaxed);
            if motif > 0 {
                eprintln!(
                    "[IsoForge] Strand from splice motifs: {} spliced reads, {} corrected against alignment orientation ({:.1}%)",
                    motif, motif_flipped, 100.0 * motif_flipped as f64 / motif as f64,
                );
            }
            let inferred = crate::bam_reader::STRAND_INFERRED.load(Relaxed);
            let flipped = crate::bam_reader::STRAND_FLIPPED.load(Relaxed);
            if inferred > 0 {
                eprintln!(
                    "[IsoForge] Strand from polyA tail: {} reads resolved ({:.1}% of {}), \
                     {} corrected against alignment orientation ({:.1}%)",
                    inferred,
                    100.0 * inferred as f64 / total_reads_count.max(1) as f64,
                    total_reads_count,
                    flipped,
                    100.0 * flipped as f64 / inferred as f64,
                );
            }
        }
        eprintln!(
            "[IsoForge] Memory: converted {} reads to compact clustering format, \
             {} annotations streamed to disk",
            total_reads_count,
            n_annotations
        );

        // ── Unique-parent ISM assignment ─────────────────────────────────
        // Runs after all datasets are read.
        if self.config.assign_unique_ism != AssignUniqueIsm::Off && !ism_parents.is_empty() {
            let (moved, dissolved) = self.reassign_unique_ism(&mut cores, &ism_parents);
            eprintln!(
                "[IsoForge] Unique-parent ISM assignment ({:?}): {} of {} unique-parent reads \
                 moved onto their annotated transcript, {} ISM models dissolved",
                self.config.assign_unique_ism, moved, ism_parents.len(), dissolved
            );
            // Rebuild abundance from the reassigned reads.
            all_results.abundance.clear();
            for core in &cores {
                let key = (core.gene_id, core.transcript_id, Arc::clone(&core.dataset));
                *all_results.abundance.entry(key).or_insert(0) += 1;
            }
        }
        drop(ism_parents);

        // ── Phase 3: TSS/TES clustering ──────────────────────────────────
        // Filled in by clustering; empty when clustering is disabled.
        let mut join_results: Vec<ReadEndSiteInfo> = Vec::new();
        // Transcripts merged away by splitting, each with the transcript it merged into.
        let mut tx_merges: HashMap<TranscriptId, TranscriptId> = HashMap::default();
        if let Some(ref cluster_cfg) = self.cluster_config {
            // Load CAGE peak index if configured
            let cage_index = if let Some(ref cage_path) = cluster_cfg.cage_peaks_file {
                let path = std::path::Path::new(cage_path);
                match CagePeakIndex::from_bed(path) {
                    Ok(idx) => {
                        eprintln!("[IsoForge] CAGE peaks loaded: {} peaks from {:?}", idx.n_peaks, cage_path);
                        Some(idx)
                    }
                    Err(e) => {
                        eprintln!("[IsoForge] WARNING: Failed to load CAGE peaks from {:?}: {}", cage_path, e);
                        None
                    }
                }
            } else {
                eprintln!("[IsoForge] CAGE peaks: not provided (5′ validation by cluster-only)");
                None
            };

            // Without poly(A) evidence in the BAM, the poly(A) requirement is
            // disabled for the run.
            let relaxed_cfg;
            let cluster_cfg = if cluster_cfg.require_polya_tail {
                let n_polya = cores
                    .iter()
                    .filter(|r| {
                        r.polya_tail_len.map_or(false, |v| v > 0)
                            || r.polya_tail_softclip.map_or(false, |v| v > 0)
                    })
                    .count();

                if n_polya == 0 {
                    eprintln!(
                        "[IsoForge] PolyA tail: no pt:i tags and no soft-clip polyA detected in \
                         any of {} reads — polyA requirement disabled for this run; 3′ ends will \
                         be validated by cluster support alone.",
                        total_reads_count
                    );
                    relaxed_cfg = ClusterConfig { require_polya_tail: false, ..cluster_cfg.clone() };
                    &relaxed_cfg
                } else {
                    eprintln!(
                        "[IsoForge] PolyA tail: REQUIRED ({}/{} reads carry polyA evidence; \
                         reads without polyA → rejected)",
                        n_polya, total_reads_count
                    );
                    cluster_cfg
                }
            } else {
                eprintln!("[IsoForge] PolyA tail: optional (no enforcement)");
                cluster_cfg
            };
            // Starts from full-length reads, when the library carries tails.
            let (n_spliced, n_spliced_tailed) = cores.iter().filter(|r| r.n_exons > 1).fold((0usize, 0usize), |(n, k), r| {
                let tailed = r.polya_tail_len.map_or(false, |v| v > 0) || r.polya_tail_softclip.map_or(false, |v| v > 0);
                (n + 1, k + tailed as usize)
            });
            let tails_present = n_spliced > 0 && n_spliced_tailed as f64 >= FULL_LENGTH_START_RATE * n_spliced as f64;
            self.tails_present.store(tails_present, std::sync::atomic::Ordering::Relaxed);
            // How often the 5' ends of tailed and tail-less spliced reads reach an annotated start,
            // for the run summary.
            let mut starts: HashMap<(Arc<str>, Strand), Vec<i64>> = HashMap::default();
            for tr in self.refs.transcript_tree.records() {
                let s = if tr.strand == Strand::Plus { tr.min_pos } else { tr.max_pos };
                starts.entry((Arc::from(tr.chromosome.as_str()), tr.strand)).or_default().push(s);
            }
            for v in starts.values_mut() { v.sort_unstable(); }
            let reaches_start = |r: &ObsCore| {
                starts.get(&(r.location.clone(), r.strand)).is_some_and(|v| {
                    let i = v.partition_point(|&x| x < r.start_pos - ANNOTATED_START_WINDOW);
                    i < v.len() && v[i] <= r.start_pos + ANNOTATED_START_WINDOW
                })
            };
            let (mut tailed_n, mut tailed_at, mut untailed_n, mut untailed_at) = (0usize, 0usize, 0usize, 0usize);
            for r in cores.iter().filter(|r| r.n_exons > 1) {
                let tailed = r.polya_tail_len.map_or(false, |v| v > 0) || r.polya_tail_softclip.map_or(false, |v| v > 0);
                let at = reaches_start(r);
                if tailed { tailed_n += 1; tailed_at += at as usize; } else { untailed_n += 1; untailed_at += at as usize; }
            }
            let (tailed_rate, untailed_rate) = (tailed_at as f64 / tailed_n.max(1) as f64, untailed_at as f64 / untailed_n.max(1) as f64);
            let full_length_starts = tails_present
                && (n_spliced_tailed as f64) <= FULL_LENGTH_START_MAX_RATE * n_spliced as f64;
            self.full_length_starts.store(full_length_starts, std::sync::atomic::Ordering::Relaxed);
            eprintln!(
                "[IsoForge] Starts from full-length reads: {} — {} of {} spliced reads carry a poly(A) tail ({:.1}%); 5' ends at an annotated start: {:.1}% of tailed, {:.1}% of tail-less spliced reads",
                if full_length_starts { "on (a spliced read places a TSS only with its tail present)" } else if tails_present { "off, tails on nearly every read (every read's 5' end counts)" } else { "off, no tails (every read's 5' end counts)" },
                n_spliced_tailed, n_spliced, 100.0 * n_spliced_tailed as f64 / n_spliced.max(1) as f64, 100.0 * tailed_rate, 100.0 * untailed_rate
            );
            let full_length_cfg = ClusterConfig { tss_from_full_length: full_length_starts, ..cluster_cfg.clone() };
            let cluster_cfg = &full_length_cfg;

            // Build novelty hints from Phase 2 results for filtering
            let hints = self.build_filter_hints(&cores, cluster_cfg);
            let (join, mut tss_clusters, mut tes_clusters) = self.run_clustering(
                &cores, &dataset_ranges, total_reads_count, cluster_cfg,
                Some(&hints), cage_index.as_ref(),
            )?;
            join_results = join;

            // Counts the reads joined to a TSS and to a TES cluster.
            let mut n_tss_joined = 0u64;
            let mut n_tes_joined = 0u64;
            for idx in 0..cores.len() {
                if let Some(info) = join_results.get(idx) {
                    if info.tss_site.is_some() {
                        n_tss_joined += 1;
                    }
                    if info.polya_site.is_some() {
                        n_tes_joined += 1;
                    }
                }
            }
            let total = cores.len() as u64;
            eprintln!("[IsoForge] Joined reads: {}/{} TSS ({:.1}%), {}/{} TES ({:.1}%)",
                     n_tss_joined, total, 100.0 * n_tss_joined as f64 / total.max(1) as f64,
                     n_tes_joined, total, 100.0 * n_tes_joined as f64 / total.max(1) as f64);

            // ── Split end-validation from existence ──────────────────────
            // Reads that pass the evidence gates decide whether a transcript exists;
            // once it does, every read assigned to it counts towards its depth.
            // Transcripts supported only by unvalidated reads are dropped.
            let validated: Vec<bool> = (0..cores.len())
                .map(|i| join_results.get(i).map_or(false, |info| info.is_supported))
                .collect();

            keep = if self.config.strict_full_length {
                validated.clone()
            } else {
                // A read that may be a fragment of a longer nested isoform does not
                // establish the short isoform on its own, but counts once it is anchored.
                let mut anchored: rustc_hash::FxHashSet<TranscriptId> = cores
                    .iter()
                    .enumerate()
                    .filter(|(i, c)| validated[*i] && !is_nested_ambiguous(c.hint_flags))
                    .map(|(_, c)| c.transcript_id)
                    .collect();
                let confirmed = self.three_prime_confirmed(&cores, &mut join_results, &anchored, cluster_cfg);
                anchored.extend(confirmed);
                cores.iter().map(|c| anchored.contains(&c.transcript_id)).collect()
            };

            let n_before = cores.len();
            let n_valid = validated.iter().filter(|k| **k).count();
            let n_after = keep.iter().filter(|k| **k).count();
            eprintln!(
                "[IsoForge] Full-length filter: {} reads end-validated ({:.1}%), \
                 {} retained for quantification ({:.1}%), {} rejected",
                n_valid, 100.0 * n_valid as f64 / n_before.max(1) as f64,
                n_after, 100.0 * n_after as f64 / n_before.max(1) as f64,
                n_before - n_after
            );

            // Recompute abundance after filtering
            all_results.abundance.clear();
            for (i, core) in cores.iter().enumerate() {
                if !keep[i] {
                    continue;
                }
                let key = (core.gene_id, core.transcript_id, Arc::clone(&core.dataset));
                *all_results.abundance.entry(key).or_insert(0) += 1;
            }

            log_end_consistency("before the cluster split", &all_results, &self.refs);
            // ── Transcript splitting by TSS/TES clusters ─────────────────
            let n_variants = self.split_transcripts_by_clusters(
                &mut cores, &keep, &mut join_results, &mut all_results, &mut tx_merges,
            );
            if n_variants > 0 {
                eprintln!(
                    "[IsoForge] Created {} new TSS/TES transcript variants",
                    n_variants
                );
            }


            // ── Post-filter cluster pruning ──────────────────────────────
            // Keep only clusters that still support at least one surviving read.
            {
                // Collect surviving TSS and TES cluster IDs from reads that passed all filters
                let mut surviving_tss_ids = rustc_hash::FxHashSet::default();
                let mut surviving_tes_ids = rustc_hash::FxHashSet::default();
                for idx in 0..cores.len() {
                    if !keep[idx] {
                        continue;
                    }
                    if let Some(info) = join_results.get(idx) {
                        if let Some(tss_id) = info.tss_cluster_id {
                            surviving_tss_ids.insert(tss_id);
                        }
                        if let Some(tes_id) = info.tes_cluster_id {
                            surviving_tes_ids.insert(tes_id);
                        }
                    }
                }

                let tss_before = tss_clusters.len();
                let tes_before = tes_clusters.len();

                tss_clusters.retain(|c| surviving_tss_ids.contains(&c.cluster_id));
                tes_clusters.retain(|c| surviving_tes_ids.contains(&c.cluster_id));

                let tss_pruned = tss_before - tss_clusters.len();
                let tes_pruned = tes_before - tes_clusters.len();
                eprintln!(
                    "[IsoForge] Post-filter cluster pruning: TSS {}/{} retained ({} pruned), \
                     TES {}/{} retained ({} pruned)",
                    tss_clusters.len(), tss_before, tss_pruned,
                    tes_clusters.len(), tes_before, tes_pruned,
                );
            }

            // ── Write cluster outputs (post-filter) ─────────────────────────
            // Cluster files go into a sibling _TSS_TES_Clusters directory.
            let cluster_dir = {
                let db_name = self.output_dir.file_name()
                    .and_then(|n| n.to_str())
                    .unwrap_or("");
                let prefix = db_name.strip_suffix("_isoforge_db").unwrap_or(db_name);
                self.output_dir.parent()
                    .unwrap_or(&self.output_dir)
                    .join(format!("{}_TSS_TES_Clusters", prefix))
            };
            {
                std::fs::create_dir_all(&cluster_dir)?;

                let mut all_clusters = Vec::new();
                all_clusters.extend_from_slice(&tss_clusters);
                all_clusters.extend_from_slice(&tes_clusters);

                let clusters_path = cluster_dir.join("tss_tes_clusters.parquet");
                tss_tes_clustering::write_clusters(&clusters_path, &all_clusters)?;
                eprintln!("[IsoForge] Wrote {} TSS/TES clusters to {:?}",
                         all_clusters.len(), clusters_path);

                let tes_bed_path = cluster_dir.join("tes_clusters.bed");
                tss_tes_clustering::write_clusters_bed(&tes_bed_path, &tes_clusters)?;
                let tss_bed_path = cluster_dir.join("tss_clusters.bed");
                tss_tes_clustering::write_clusters_bed(&tss_bed_path, &tss_clusters)?;
                eprintln!("[IsoForge] Wrote cluster BED files (post-filter)");
            }

            // ── Per-sample & per-treatment cluster output ────────────────
            // Per-group BED and bedGraph files for clusters with surviving reads.
            {
                eprintln!("[IsoForge] ── Per-sample/treatment cluster output ──");

                // Collect datasets
                let mut datasets: Vec<String> = cores.iter()
                    .enumerate()
                    .filter(|(i, _)| keep[*i])
                    .map(|(_, c)| c.dataset.to_string())
                    .collect::<HashSet<_>>()
                    .into_iter()
                    .collect();
                datasets.sort();
                eprintln!("[IsoForge] Datasets: {:?}", datasets);

                // Infer treatment groups
                let treatment_groups = tss_tes_clustering::infer_treatment_groups(&datasets);
                eprintln!("[IsoForge] Treatment groups:");
                for (tname, members) in &treatment_groups {
                    eprintln!("  {} → {:?}", tname, members);
                }

                // Accumulate per-sample cluster counts
                let mut sample_tes_counts: HashMap<(String, i64), u64> = HashMap::default();
                let mut sample_tss_counts: HashMap<(String, i64), u64> = HashMap::default();
                let mut sample_total_reads: HashMap<String, u64> = HashMap::default();
                // Per-sample per-cluster dataset membership (n_datasets in treatment groups)
                let mut sample_ds_cluster_tes: HashMap<(String, i64), HashSet<String>> = HashMap::default();
                let mut sample_ds_cluster_tss: HashMap<(String, i64), HashSet<String>> = HashMap::default();

                for (idx, core) in cores.iter().enumerate() {
                    if !keep[idx] {
                        continue;
                    }
                    let ds = core.dataset.to_string();
                    *sample_total_reads.entry(ds.clone()).or_insert(0) += 1;

                    if let Some(info) = join_results.get(idx) {
                        if let Some(tes_id) = info.tes_cluster_id {
                            *sample_tes_counts.entry((ds.clone(), tes_id)).or_insert(0) += 1;
                            sample_ds_cluster_tes.entry((ds.clone(), tes_id))
                                .or_default().insert(ds.clone());
                        }
                        if let Some(tss_id) = info.tss_cluster_id {
                            *sample_tss_counts.entry((ds.clone(), tss_id)).or_insert(0) += 1;
                            sample_ds_cluster_tss.entry((ds.clone(), tss_id))
                                .or_default().insert(ds.clone());
                        }
                    }
                }

                // Per-sample stats
                let (sample_tss_stats, sample_tes_stats) =
                    tss_tes_clustering::compute_per_group_cluster_stats(
                        &tss_clusters, &tes_clusters,
                        &sample_tss_counts, &sample_tes_counts,
                        &sample_total_reads,
                        &{
                            let mut combined = sample_ds_cluster_tes.clone();
                            combined.extend(sample_ds_cluster_tss.clone());
                            combined
                        },
                    );

                // Per-sample BED files
                let sample_dir = cluster_dir.join("sample");
                std::fs::create_dir_all(&sample_dir)?;
                for ds in &datasets {
                    let tes_path = sample_dir.join(format!("{}_tes.bed", ds));
                    tss_tes_clustering::write_group_tes_bed(
                        &tes_path, &tes_clusters, &sample_tes_stats, ds,
                    )?;
                    let tss_path = sample_dir.join(format!("{}_tss.bed", ds));
                    tss_tes_clustering::write_group_tss_bed(
                        &tss_path, &tss_clusters, &sample_tss_stats, ds,
                    )?;
                }
                eprintln!("[IsoForge] Wrote {} per-sample BED files to {:?}",
                         datasets.len() * 2, sample_dir);

                // Per-treatment stats and BED files
                let dataset_dir = cluster_dir.join("dataset");
                std::fs::create_dir_all(&dataset_dir)?;

                if treatment_groups.len() >= 2 {
                    let mut treat_tes_counts: HashMap<(String, i64), u64> = HashMap::default();
                    let mut treat_tss_counts: HashMap<(String, i64), u64> = HashMap::default();
                    let mut treat_total_reads: HashMap<String, u64> = HashMap::default();
                    let mut treat_ds_cluster: HashMap<(String, i64), HashSet<String>> = HashMap::default();

                    for (tname, members) in &treatment_groups {
                        for ds in members {
                            // Aggregate sample counts into treatment counts
                            if let Some(&total) = sample_total_reads.get(ds) {
                                *treat_total_reads.entry(tname.clone()).or_insert(0) += total;
                            }
                            for ((sds, cid), cnt) in &sample_tes_counts {
                                if sds == ds {
                                    *treat_tes_counts.entry((tname.clone(), *cid)).or_insert(0) += cnt;
                                    treat_ds_cluster.entry((tname.clone(), *cid))
                                        .or_default().insert(ds.clone());
                                }
                            }
                            for ((sds, cid), cnt) in &sample_tss_counts {
                                if sds == ds {
                                    *treat_tss_counts.entry((tname.clone(), *cid)).or_insert(0) += cnt;
                                    treat_ds_cluster.entry((tname.clone(), *cid))
                                        .or_default().insert(ds.clone());
                                }
                            }
                        }
                    }

                    let treat_names: Vec<String> = treatment_groups.iter().map(|(t, _)| t.clone()).collect();
                    let (treat_tss_stats, treat_tes_stats) =
                        tss_tes_clustering::compute_per_group_cluster_stats(
                            &tss_clusters, &tes_clusters,
                            &treat_tss_counts, &treat_tes_counts,
                            &treat_total_reads,
                            &treat_ds_cluster,
                        );

                    for tname in &treat_names {
                        let tes_path = dataset_dir.join(format!("{}_tes.bed", tname));
                        tss_tes_clustering::write_group_tes_bed(
                            &tes_path, &tes_clusters, &treat_tes_stats, tname,
                        )?;
                        let tss_path = dataset_dir.join(format!("{}_tss.bed", tname));
                        tss_tes_clustering::write_group_tss_bed(
                            &tss_path, &tss_clusters, &treat_tss_stats, tname,
                        )?;
                    }
                    eprintln!("[IsoForge] Wrote {} per-treatment BED files to {:?}",
                             treat_names.len() * 2, dataset_dir);
                } else {
                    eprintln!("[IsoForge] Only 1 treatment group detected, skipping per-treatment output");
                }

                // bedGraph count files (per-sample and merged, strand-separated), converted
                // to BigWig when bedGraphToBigWig is on PATH. Skipped with --no_end_tracks.
                let write_end_tracks = self.cluster_config.as_ref().map_or(true, |c| c.write_end_tracks);
                if write_end_tracks {
                let counts_dir = cluster_dir.join("counts");
                std::fs::create_dir_all(&counts_dir)?;

                // Load chromosome sizes from genome .fai
                let chrom_sizes = if let Some(ref genome_path) = self.genome_path {
                    let fai_path = genome_path.with_extension(
                        format!("{}.fai", genome_path.extension().unwrap_or_default().to_str().unwrap_or(""))
                    );
                    if fai_path.exists() {
                        tss_tes_clustering::load_chrom_sizes(&fai_path)?
                    } else {
                        // Try appending .fai directly
                        let fai2 = PathBuf::from(format!("{}.fai", genome_path.display()));
                        if fai2.exists() {
                            tss_tes_clustering::load_chrom_sizes(&fai2)?
                        } else {
                            eprintln!("[IsoForge] Warning: No .fai index found, skipping bedGraph output");
                            Vec::new()
                        }
                    }
                } else {
                    eprintln!("[IsoForge] Warning: No genome path set, skipping bedGraph output");
                    Vec::new()
                };

                if !chrom_sizes.is_empty() {
                    // Write chrom.sizes file for bedGraphToBigWig
                    let chrom_sizes_path = counts_dir.join("chrom.sizes");
                    {
                        use std::io::Write;
                        let mut f = std::io::BufWriter::new(std::fs::File::create(&chrom_sizes_path)?);
                        for (name, size) in &chrom_sizes {
                            writeln!(f, "{}\t{}", name, size).map_err(|e| IsoForgeError::IoError(e))?;
                        }
                    }

                    let mut bigwig_available = false;

                    // Per-sample bedGraph (strand-separated, TES + TSS)
                    let all_datasets_plus_all: Vec<Option<&str>> = datasets.iter()
                        .map(|d| Some(d.as_str()))
                        .chain(std::iter::once(None)) // None = "all" merged
                        .collect();

                    // Accumulate strand-aware end-site positions
                    let mut tes_strand_positions: HashMap<(String, String, String, i64), u64> = HashMap::default();
                    let mut tss_strand_positions: HashMap<(String, String, String, i64), u64> = HashMap::default();

                    for (idx, core) in cores.iter().enumerate() {
                        if !keep[idx] {
                            continue;
                        }
                        let ds = core.dataset.to_string();
                        let strand_str = core.strand.as_str().to_string();

                        if let Some(info) = join_results.get(idx) {
                            if let Some(polya_site) = info.polya_site {
                                *tes_strand_positions.entry((
                                    ds.clone(), core.location.to_string(), strand_str.clone(), polya_site
                                )).or_insert(0) += 1;
                            }
                            if let Some(tss_site) = info.tss_site {
                                *tss_strand_positions.entry((
                                    ds.clone(), core.location.to_string(), strand_str.clone(), tss_site
                                )).or_insert(0) += 1;
                            }
                        }
                    }

                    // Write per-sample + merged bedGraph files
                    for ds_opt in &all_datasets_plus_all {
                        let label = ds_opt.unwrap_or("all");

                        for (end_type, strand_positions) in &[
                            ("polyA", &tes_strand_positions),
                            ("tss", &tss_strand_positions),
                        ] {
                            for (strand_label, strand_char) in &[("pos", "+"), ("neg", "-")] {
                                let mut counts: HashMap<(String, i64), u64> = HashMap::default();
                                for ((ds, chrom, strand, pos), &cnt) in *strand_positions {
                                    let ds_match = match ds_opt {
                                        Some(target_ds) => ds.as_str() == *target_ds,
                                        None => true,
                                    };
                                    if ds_match && strand.as_str() == *strand_char {
                                        *counts.entry((chrom.clone(), *pos)).or_insert(0) += cnt;
                                    }
                                }

                                if counts.is_empty() { continue; }

                                let bg_name = format!("{}_{}_counts_{}.bedGraph", label, end_type, strand_label);
                                let bg_path = counts_dir.join(&bg_name);
                                tss_tes_clustering::write_bedgraph(&bg_path, &counts, &chrom_sizes)?;

                                // Try BigWig conversion
                                let bw_name = bg_name.replace(".bedGraph", ".bw");
                                let bw_path = counts_dir.join(&bw_name);
                                match tss_tes_clustering::bedgraph_to_bigwig(&bg_path, &chrom_sizes_path, &bw_path) {
                                    Ok(true) => { bigwig_available = true; },
                                    _ => {},
                                }
                            }
                        }
                    }

                    let format_label = if bigwig_available { "BigWig" } else { "bedGraph" };
                    eprintln!("[IsoForge] Wrote per-sample + merged {} count files to {:?}",
                             format_label, counts_dir);
                }
                }
            }
        }

        // ── Phase 4: Write outputs ───────────────────────────────────────

        // Reads keep everything when clustering never ran.
        if keep.len() != cores.len() {
            keep = vec![true; cores.len()];
        }

        // A novel transcript whose reads were all removed by the full-length filter is not reported.
        let dropped_models: Vec<TranscriptModel>;
        let dropped_tx: rustc_hash::FxHashSet<TranscriptId> = {
            let supported: rustc_hash::FxHashSet<TranscriptId> = cores
                .iter()
                .enumerate()
                .filter(|(i, _)| keep[*i])
                .map(|(_, c)| c.transcript_id)
                .collect();
            let dropped: rustc_hash::FxHashSet<TranscriptId> = all_results
                .novel_transcripts
                .iter()
                .map(|t| t.transcript_id)
                .filter(|id| !supported.contains(id))
                .collect();
            dropped_models = all_results.novel_transcripts.iter().filter(|t| !supported.contains(&t.transcript_id)).cloned().collect();
            all_results.novel_transcripts.retain(|t| supported.contains(&t.transcript_id));
            if !dropped.is_empty() {
                eprintln!(
                    "[IsoForge] Dropped {} novel transcripts with no surviving read support \
                     ({} retained)",
                    dropped.len(),
                    all_results.novel_transcripts.len()
                );
            }
            dropped
        };

        // each model's first read, in file order (cores are in file order)
        let mut first_read: rustc_hash::FxHashMap<TranscriptId, ReadIdx> = rustc_hash::FxHashMap::default();
        for c in cores.iter() {
            let e = first_read.entry(c.transcript_id).or_insert(c.read_idx);
            if c.read_idx < *e { *e = c.read_idx; }
        }

        log_end_consistency("before end anchoring", &all_results, &self.refs);
        let n_anchored =
            self.anchor_novel_transcript_ends(&mut all_results, &cores, &keep, &join_results, &first_read);
        log_end_consistency("after end anchoring", &all_results, &self.refs);
        if n_anchored > 0 {
            eprintln!(
                "[IsoForge] Anchored {} novel transcript ends to validated TSS/TES peaks",
                n_anchored
            );
        }

        // Novel genes come from the registry: one per locus, with the union of the spans of its reads.
        all_results.novel_genes = {
            // sorted by locus, then by the registry key (what the gene is antisense to, or its locus bin)
            let mut genes: Vec<((u8, i64), Gene)> = self.refs.novel_gene_dict.iter().map(|e| {
                let rank = match e.key() {
                    NovelGeneKey::Antisense(g, _) => (0u8, *g),
                    NovelGeneKey::Locus(_, _, bin) => (1u8, *bin),
                };
                (rank, Gene { gene_id: e.gene_id, strand: e.strand, chromosome: e.chromosome.clone(), start: e.start, end: e.end })
            }).collect();
            genes.sort_by(|a, b| (&a.1.chromosome, a.1.strand, a.1.start, a.1.end, a.0).cmp(&(&b.1.chromosome, b.1.strand, b.1.start, b.1.end, b.0)));
            genes.into_iter().map(|(_, g)| g).collect()
        };

        let remap = renumber_novel_entities(&mut all_results, &self.refs.ref_max, &dropped_models, &first_read);
        drop(first_read);
        if std::env::var_os("ISOFORGE_END_CHECK").is_some() {
            let stray = cores.iter().enumerate().filter(|(i, c)| keep[*i] && c.transcript_id > self.refs.ref_max.transcript
                && !remap.transcripts.contains_key(&c.transcript_id)).count();
            eprintln!("[IsoForge] Renumbering: {stray} kept reads on a novel transcript outside the renumbering");
        }
        log_end_consistency("after renumbering", &all_results, &self.refs);
        remap.rewrite_annotation_ids(
            &annot_path, &dropped_tx, &tx_merges, &self.refs.ref_max,
            &self.config.id_prefix, self.config.n_places,
        )?;
        {
            let final_ids: rustc_hash::FxHashSet<TranscriptId> =
                all_results.novel_transcripts.iter().map(|m| m.transcript_id).collect();
            let rows: Vec<(TranscriptId, TranscriptId, bool)> = all_results.end_variant_parents.iter()
                .map(|&(v, parent, ok)| (remap.transcript(v), parent, ok))
                .filter(|(v, _, _)| final_ids.contains(v))
                .collect();
            let mut rows = rows;
            rows.sort_unstable();
            parquet_io::write_end_variants(&self.output_dir.join("end_variants.parquet"), &rows)?;
            log::info!("[IsoForge] Wrote {} end variants of annotated transcripts", rows.len());
        }
        // writes the list of annotated small RNAs for post-processing, or removes the list when the rule is off
        match &self.refs.small_rna {
            Some(index) => index.write(&self.output_dir)?,
            None => { let _ = std::fs::remove_file(self.output_dir.join(crate::small_rna::IDS_FILE)); }
        }

        if std::env::var_os("ISOFORGE_READ_DECISIONS").is_some() {
            parquet_io::write_read_decisions(&spill_path, &self.output_dir.join("read_decisions.tsv"),
                &cores, &join_results, &keep, |tid| remap.transcript(tid))?;
        }

        // Replay the spill, filling in the cluster columns and renumbered IDs.
        let observed_path = self.output_dir.join("observed.parquet");
        let empty_info = ReadEndSiteInfo::default();
        let n_observed = parquet_io::finalize_observed(&spill_path, &observed_path, &|idx, ids| {
            if !keep.get(idx).copied().unwrap_or(false) {
                return None;
            }
            let core = &cores[idx];
            let info = join_results.get(idx).unwrap_or(&empty_info);
            Some(parquet_io::ObservedPatch {
                gene_id: remap.gene(core.gene_id),
                transcript_id: remap.transcript(core.transcript_id),
                start_vertex: remap.vertex(ids.start_vertex),
                end_vertex: remap.vertex(ids.end_vertex),
                start_exon: remap.edge(ids.start_exon),
                end_exon: remap.edge(ids.end_exon),
                start_support: info.tss_site,
                end_support: info.polya_site,
                tss_site: info.tss_site,
                tss_feature: info.tss_feature.as_ref().map(|f| f.as_str().to_string()),
                polya_site: info.polya_site,
                polya_feature: info.polya_feature.as_ref().map(|f| f.as_str().to_string()),
                polya_signal: info.polya_signal.clone(),
                cluster_supported: Some(info.is_supported),
            })
        })?;
        let _ = std::fs::remove_file(&spill_path);
        log::info!("[IsoForge] Wrote {} observed reads to {:?}", n_observed, observed_path);
        log::info!("[IsoForge] Wrote observed reads to {:?}", observed_path);

        // Annotations were streamed to disk during Phase 1+2.
        log::info!("[IsoForge] Annotations already written to {:?}", annot_path);

        let abundance_path = self.output_dir.join("abundance.parquet");
        parquet_io::write_abundance(&abundance_path, &all_results.abundance)?;
        log::info!("[IsoForge] Wrote abundance to {:?}", abundance_path);

        let genes_path = self.output_dir.join("novel_genes.parquet");
        // in ID order
        all_results.novel_genes.sort_by_key(|g| g.gene_id);
        parquet_io::write_novel_genes(&genes_path, &all_results.novel_genes)?;
        log::info!("[IsoForge] Wrote {} novel genes", all_results.novel_genes.len());

        let transcripts_path = self.output_dir.join("novel_transcripts.parquet");
        parquet_io::write_novel_transcripts(&transcripts_path, &all_results.novel_transcripts)?;
        eprintln!(
            "[IsoForge] Wrote {} novel transcripts",
            all_results.novel_transcripts.len()
        );

        // The models as hyperedges: each novel transcript's vertex chain over the annotated and novel vertices
        {
            let (mut pos, mut ends) = crate::annotation_hypergraph::vertex_and_edge_maps(&self.refs.location_dict, &self.refs.edge_dict);
            pos.extend(all_results.novel_locations.iter().map(|l| (l.location_id, l.position)));
            ends.extend(all_results.novel_edges.iter().map(|e| (e.edge_id, (e.v1, e.v2, e.strand))));
            let rows: Vec<_> = all_results.novel_transcripts.iter()
                .filter_map(|m| crate::annotation_hypergraph::derive_row(m, &pos, &ends)).collect();
            if rows.len() < all_results.novel_transcripts.len() {
                log::warn!("[IsoForge] {} novel transcripts without a complete vertex chain",
                           all_results.novel_transcripts.len() - rows.len());
            }
            crate::annotation_hypergraph::write_hyperedges(&self.output_dir.join("novel_hyperedges.parquet"), &rows)?;
            eprintln!("[IsoForge] Wrote {} novel hyperedges", rows.len());
        }

        // Deduplicate and write novel edges
        let mut seen_edge_ids = rustc_hash::FxHashSet::default();
        let mut unique_edges: Vec<Edge> = Vec::new();
        for edge in all_results.novel_edges {
            if seen_edge_ids.insert(edge.edge_id) {
                unique_edges.push(edge);
            }
        }
        let edges_path = self.output_dir.join("novel_edges.parquet");
        parquet_io::write_novel_edges(&edges_path, &unique_edges)?;
        eprintln!("[IsoForge] Wrote {} novel edges", unique_edges.len());

        // Deduplicate and write novel locations
        let mut seen_loc_ids = rustc_hash::FxHashSet::default();
        let mut unique_locs: Vec<Location> = Vec::new();
        for loc in all_results.novel_locations {
            if seen_loc_ids.insert(loc.location_id) {
                unique_locs.push(loc);
            }
        }
        let locs_path = self.output_dir.join("novel_locations.parquet");
        parquet_io::write_novel_locations(&locs_path, &unique_locs)?;
        eprintln!("[IsoForge] Wrote {} novel locations", unique_locs.len());

        log::info!("[IsoForge] Done. {} total reads processed.", all_results.total_reads);

        Ok(())
    }

    /// Anchor novel transcript ends to their validated TSS/TES cluster peaks.
    ///
    /// Each end moves to the most-voted peak its supporting reads joined, or to their most-voted
    /// boundary where they joined none. A vertex and terminal exon edge are created at the
    /// anchored position when none exists. Returns the number of ends moved.
    fn anchor_novel_transcript_ends(
        &self,
        results: &mut IsoForgeResults,
        cores: &[ObsCore],
        keep: &[bool],
        join_results: &[ReadEndSiteInfo],
        first_read: &rustc_hash::FxHashMap<TranscriptId, ReadIdx>,
    ) -> usize {
        let mut votes = AnchorVotes::default();
        votes.collect(cores, keep, join_results);

        // Annotated and novel edges by ID, for the terminal exons that are re-pointed.
        let mut edge_by_id: rustc_hash::FxHashMap<EdgeId, Edge> = self
            .refs
            .edge_dict
            .values()
            .map(|e| (e.edge_id, e.clone()))
            .collect();
        for e in results.novel_edges.iter() {
            edge_by_id.insert(e.edge_id, e.clone());
        }

        let mut new_locations: Vec<Location> = Vec::new();
        let mut new_edges: Vec<Edge> = Vec::new();
        let mut vertex_at: rustc_hash::FxHashMap<(String, i64), VertexId> = results
            .novel_locations
            .iter()
            .map(|l| ((l.chromosome.clone(), l.position), l.location_id))
            .collect();
        let mut edge_at: rustc_hash::FxHashMap<(VertexId, VertexId), EdgeId> = edge_by_id
            .values()
            .filter(|e| e.edge_type == EdgeType::Exon)
            .map(|e| ((e.v1, e.v2), e.edge_id))
            .collect();

        let genome_build = self.config.build.clone();
        let mut n_anchored = 0usize;

        // The models in the order of their first reads in the files: a new terminal exon that the ends of two models
        // move onto is made with the strand of the earlier.
        let mut order: Vec<usize> = (0..results.novel_transcripts.len()).collect();
        order.sort_by_key(|&i| {
            let t = &results.novel_transcripts[i];
            (first_read.get(&t.transcript_id).copied().unwrap_or(ReadIdx::MAX), t.start_pos, t.end_pos, t.n_exons)
        });

        for i in order {
            let tx = &mut results.novel_transcripts[i];
            for is_start in [true, false] {
                let anchor = if is_start {
                    votes.modal_tss(tx.transcript_id)
                } else {
                    votes.modal_tes(tx.transcript_id)
                };
                let Some(anchor) = anchor else { continue };

                // The model's terminal exon edge at this end.
                let term_edge_id = if is_start { tx.start_exon } else { tx.end_exon };
                let Some(term_edge) = edge_by_id.get(&term_edge_id).cloned() else { continue };
                // A single-exon model's one edge is both its first and its last
                // exon: both are re-pointed.
                let single_exon = tx.start_exon == tx.end_exon;
                // The end of the terminal exon that stays put: edges run 5'->3', so the first exon's 3' vertex and
                // the last exon's 5' vertex.
                let inner_vertex = if is_start { term_edge.v2 } else { term_edge.v1 };

                let key = (tx.chromosome.clone(), anchor);
                let new_vertex = match vertex_at.get(&key) {
                    Some(&v) => v,
                    None => match self.refs.location_dict.get_vertex(&tx.chromosome, anchor) {
                        Some(v) => {
                            vertex_at.insert(key, v);
                            v
                        }
                        None => {
                            let v = self.refs.counters.vertex.increment();
                            new_locations.push(Location {
                                location_id: v,
                                genome_build: genome_build.clone(),
                                chromosome: tx.chromosome.clone(),
                                position: anchor,
                            });
                            vertex_at.insert(key, v);
                            v
                        }
                    },
                };

                let (v1, v2) = if is_start { (new_vertex, inner_vertex) } else { (inner_vertex, new_vertex) };
                let new_edge = match edge_at.get(&(v1, v2)) {
                    Some(&e) => e,
                    None => {
                        let e = self.refs.counters.edge.increment();
                        new_edges.push(Edge {
                            edge_id: e,
                            v1,
                            v2,
                            edge_type: EdgeType::Exon,
                            strand: term_edge.strand,
                        });
                        edge_at.insert((v1, v2), e);
                        e
                    }
                };

                if is_start {
                    if tx.start_pos != anchor {
                        n_anchored += 1;
                    }
                    tx.start_vertex = new_vertex;
                    tx.start_exon = new_edge;
                    tx.start_pos = anchor;
                } else {
                    if tx.end_pos != anchor {
                        n_anchored += 1;
                    }
                    tx.end_vertex = new_vertex;
                    tx.end_exon = new_edge;
                    tx.end_pos = anchor;
                }
                if single_exon {
                    tx.start_exon = new_edge;
                    tx.end_exon = new_edge;
                }
            }
        }

        results.novel_locations.extend(new_locations);
        results.novel_edges.extend(new_edges);
        n_anchored
    }

    /// Annotated transcripts, not already anchored, whose 3' end is confirmed although no read reaches their start:
    /// reads that fit the transcript alone end at one of its 3' clusters, at least [`THREE_PRIME_CONFIRM_READS`] per
    /// dataset in two datasets (or in the only one). Those reads are then marked end-validated.
    fn three_prime_confirmed(
        &self,
        cores: &[ObsCore],
        join_results: &mut [ReadEndSiteInfo],
        anchored: &rustc_hash::FxHashSet<TranscriptId>,
        cfg: &ClusterConfig,
    ) -> rustc_hash::FxHashSet<TranscriptId> {
        const FITS_ALONE: u8 = HINT_IS_ISM | HINT_ISM_UNIQUE | HINT_ISM_3P_INTACT;
        let mut per: rustc_hash::FxHashMap<TranscriptId, rustc_hash::FxHashMap<&str, usize>> = rustc_hash::FxHashMap::default();
        let mut eligible: Vec<usize> = Vec::new();
        let mut datasets: rustc_hash::FxHashSet<&str> = rustc_hash::FxHashSet::default();
        for (i, c) in cores.iter().enumerate() {
            datasets.insert(&c.dataset);
            if c.transcript_id > self.refs.ref_max.transcript || anchored.contains(&c.transcript_id) || c.is_internal_priming {
                continue;
            }
            let Some(info) = join_results.get(i) else { continue };
            if info.rejection_reason.as_deref() != Some("degraded_5p") || info.tes_cluster_id.is_none() {
                continue;
            }
            let fits = c.hint_flags & HINT_IS_KNOWN != 0 || c.hint_flags & FITS_ALONE == FITS_ALONE;
            let tailed = c.polya_len().or(c.polya_softclip()).map_or(false, |v| v > 0);
            if !fits || is_nested_ambiguous(c.hint_flags) || (cfg.require_polya_tail && !tailed) {
                continue;
            }
            *per.entry(c.transcript_id).or_default().entry(&c.dataset).or_default() += 1;
            eligible.push(i);
        }
        let need = datasets.len().min(2);
        let confirmed: rustc_hash::FxHashSet<TranscriptId> = per
            .into_iter()
            .filter(|(_, ds)| ds.values().filter(|&&n| n >= THREE_PRIME_CONFIRM_READS).count() >= need)
            .map(|(tid, _)| tid)
            .collect();
        let mut n_reads = 0usize;
        for i in eligible {
            if confirmed.contains(&cores[i].transcript_id) {
                let info = &mut join_results[i];
                info.is_supported = true;
                info.rejection_reason = None;
                n_reads += 1;
            }
        }
        eprintln!("[IsoForge] 3' ends confirm {} annotated transcripts that no read reaches the start of ({} reads)", confirmed.len(), n_reads);
        confirmed
    }

    /// A truncated read that was moved onto its unique annotated parent.
    fn is_reassigned_ism(&self, c: &ObsCore) -> bool {
        c.hint_flags & HINT_IS_ISM != 0 && c.transcript_id <= self.refs.ref_max.transcript
    }

    /// Move unique-parent ISM reads onto their annotated transcript according
    /// to `assign_unique_ism`. Returns (reads moved, ISM models left empty).
    fn reassign_unique_ism(&self, cores: &mut [ObsCore], ism_parents: &[(usize, TranscriptId)]) -> (usize, usize) {
        let mode = self.config.assign_unique_ism;
        let index = self.refs.parent_chains.as_ref().expect("parent index built when assignment is on");
        let mut model_reads: rustc_hash::FxHashMap<TranscriptId, usize> = rustc_hash::FxHashMap::default();
        let mut fsm_parents: rustc_hash::FxHashSet<TranscriptId> = rustc_hash::FxHashSet::default();
        for c in cores.iter() {
            if c.hint_flags & HINT_IS_ISM != 0 {
                *model_reads.entry(c.transcript_id).or_insert(0) += 1;
            } else if c.hint_flags & HINT_IS_KNOWN != 0 {
                fsm_parents.insert(c.transcript_id);
            }
        }
        let mut touched: rustc_hash::FxHashSet<TranscriptId> = rustc_hash::FxHashSet::default();
        let mut moved = 0usize;
        for &(i, parent) in ism_parents {
            let c = &mut cores[i];
            let eligible = match mode {
                AssignUniqueIsm::Off => false,
                AssignUniqueIsm::Always => true,
                AssignUniqueIsm::BelowNovelFloor => {
                    model_reads.get(&c.transcript_id).copied().unwrap_or(0) <= self.config.ism_reassign_max_reads
                }
                AssignUniqueIsm::NoFsmParent => !fsm_parents.contains(&parent),
            };
            if !eligible {
                continue;
            }
            let Some(chain) = index.get(parent) else {
                continue;
            };
            if c.transcript_id == parent {
                continue; // already credited to the parent at classification
            }
            touched.insert(c.transcript_id);
            c.transcript_id = parent;
            c.gene_id = chain.gene_id;
            c.hint_flags |= HINT_IS_KNOWN;
            moved += 1;
        }
        let mut remaining: rustc_hash::FxHashSet<TranscriptId> = rustc_hash::FxHashSet::default();
        for c in cores.iter() {
            if touched.contains(&c.transcript_id) {
                remaining.insert(c.transcript_id);
            }
        }
        (moved, touched.len() - remaining.len())
    }

    /// Build ReadFilterHints from annotation results.
    ///
    /// Flags internally primed reads, ISM-prefix reads (truncated 3' end,
    /// excluded from TES counting), ISM-suffix reads (truncated 5' end, excluded
    /// from TSS counting), and degraded reads whose 5' or 3' end lies more than
    /// `degradation_tolerance` bp from their transcript's modal position.
    fn build_filter_hints(
        &self,
        cores: &[ObsCore],
        config: &ClusterConfig,
    ) -> tss_tes_clustering::ReadFilterHints {
        let mut hints = tss_tes_clustering::ReadFilterHints::default();

        // Annotated TES positions per (chromosome, strand); reads ending near one
        // are never flagged as internally primed.
        let window = config.annotated_tes_polya_window;
        let mut annotated_tes: HashMap<(&str, Strand), Vec<i64>> = HashMap::default();
        if window > 0 {
            for t in self.refs.transcript_dict.values() {
                // the 3' end vertex of the transcript's hyperedge
                if let Some(row) = self.refs.hypergraph.get(t.transcript_id).and_then(|h| h.row.as_ref()) {
                    annotated_tes.entry((t.chromosome.as_str(), row.strand)).or_default().push(*row.positions.last().unwrap());
                }
            }
            for v in annotated_tes.values_mut() { v.sort_unstable(); v.dedup(); }
        }
        let mut n_primed_at_annotated_tes = 0usize;

        for (idx, annot) in cores.iter().enumerate() {
            let read_idx = idx as ReadIdx;

            // Internal priming
            if annot.is_internal_priming {
                let at_annotated_tes = annotated_tes.get(&(&*annot.location, annot.strand)).is_some_and(|ends| {
                    let i = ends.partition_point(|&x| x < annot.end_pos - window);
                    ends.get(i).is_some_and(|&x| x <= annot.end_pos + window)
                });
                if at_annotated_tes {
                    n_primed_at_annotated_tes += 1;
                } else {
                    hints.is_internal_priming.insert(read_idx);
                }
            }

            // ISM subtype from hint_flags
            hints.mark_truncation(read_idx, annot.hint_flags);

            // Known/FSM reads, including truncated reads moved onto their unique
            // parent, bypass weak-cluster revocation.
            if annot.hint_flags & HINT_IS_KNOWN != 0 {
                hints.is_known_read.insert(read_idx);
            }
        }

        if self.config.terminal_exon_rescue {
            eprintln!("[IsoForge] Terminal exons placed from soft clips: {} reads at the low genomic end, {} at the high end",
                      crate::terminal_rescue::RESCUED_LOW.load(std::sync::atomic::Ordering::Relaxed),
                      crate::terminal_rescue::RESCUED_HIGH.load(std::sync::atomic::Ordering::Relaxed));
        }
        eprintln!(
            "[IsoForge] Filter hints: {} internal-priming ({} more at an annotated TES kept), {} 5prime-truncated (TSS-excluded), {} 3prime-truncated (TES-excluded), {} Known/FSM (zombie-bypass)",
            hints.is_internal_priming.len(), n_primed_at_annotated_tes,
            hints.exclude_tss_evidence.len(),
            hints.exclude_tes_evidence.len(),
            hints.is_known_read.len(),
        );

        // ── Degradation-aware read filtering ────────────────────────────
        // Reads whose 5' or 3' end lies more than `degradation_tolerance` from
        // their transcript's modal position are flagged as degraded.
        if config.filter_degraded_reads {
            let tolerance = config.degradation_tolerance;


            // Group (read_idx, pos_5p, pos_3p) by transcript_id.
            let mut tx_reads: HashMap<TranscriptId, Vec<(ReadIdx, i64, i64)>> = HashMap::default();
            for (idx, annot) in cores.iter().enumerate() {
                // A parent's modal ends are defined by its full-length reads.
                if self.is_reassigned_ism(annot) {
                    continue;
                }
                tx_reads
                    .entry(annot.transcript_id)
                    .or_default()
                    .push((idx as ReadIdx, annot.start_pos, annot.end_pos));
            }

            // Per transcript, pass 1 flags reads beyond tolerance from the modal 5'/3'
            // position; pass 2 rescues outliers that form coherent sub-clusters
            // (≥3 reads within 2 bp) as alternative ends.
            let mut n_degraded_5p = 0u64;
            let mut n_degraded_3p = 0u64;
            let mut n_rescued_5p = 0u64;
            let mut n_rescued_3p = 0u64;
            let mut n_subclusters_5p = 0u64;
            let mut n_subclusters_3p = 0u64;
            let mut n_tx_processed = 0u64;
            let mut n_tx_no_mode = 0u64;

            const SUBCLUSTER_GAP: i64 = 2;       // max bp gap between reads in a sub-cluster
            const MIN_SUBCLUSTER_READS: usize = 3; // min reads to rescue a sub-cluster

            // Annotated (TSS, TES) per reference transcript, in transcript
            // orientation, when a full-length read near its annotated end is
            // exempt from the read-derived mode.
            let annotated_ends: HashMap<TranscriptId, (i64, i64)> = if config.known_ends_not_degraded {
                let positions: rustc_hash::FxHashMap<i64, i64> = self.refs.location_dict.chromosomes.values()
                    .flat_map(|c| c.positions.iter().map(|&(p, v)| (v, p))).collect();
                self.refs.transcript_dict.values()
                    .filter_map(|t| Some((t.transcript_id, (*positions.get(&t.start_vertex)?, *positions.get(&t.end_vertex)?))))
                    .collect()
            } else {
                HashMap::default()
            };
            let end_window = config.reference_end_window;
            let mut n_exempt_5p = 0u64;
            let mut n_exempt_3p = 0u64;

            for (tx_id, reads) in &tx_reads {
                // Transcripts with fewer than 3 reads are skipped
                if reads.len() < 3 {
                    continue;
                }
                n_tx_processed += 1;

                // Find modal 5' position (most frequent exact coordinate)
                let mut pos_5p_counts: HashMap<i64, usize> = HashMap::default();
                let mut pos_3p_counts: HashMap<i64, usize> = HashMap::default();
                for (_, p5, p3) in reads {
                    *pos_5p_counts.entry(*p5).or_insert(0) += 1;
                    *pos_3p_counts.entry(*p3).or_insert(0) += 1;
                }

                // A mode requires at least two reads at the same coordinate.
                if pos_5p_counts.values().all(|&n| n < 2) && pos_3p_counts.values().all(|&n| n < 2) {
                    n_tx_no_mode += 1;
                    continue;
                }

                // Mode = position with highest count (ties broken by leftmost for 5', rightmost for 3')
                let mode_5p = pos_5p_counts
                    .iter()
                    .max_by(|a, b| a.1.cmp(b.1).then_with(|| b.0.cmp(a.0)))
                    .map(|(pos, _)| *pos)
                    .unwrap();
                let mode_3p = pos_3p_counts
                    .iter()
                    .max_by(|a, b| a.1.cmp(b.1).then_with(|| a.0.cmp(b.0)))
                    .map(|(pos, _)| *pos)
                    .unwrap();

                // ── Pass 1: Identify outlier reads beyond tolerance ──
                let mut outliers_5p: Vec<(ReadIdx, i64)> = Vec::new();
                let mut outliers_3p: Vec<(ReadIdx, i64)> = Vec::new();

                // A full-length read within the window of an annotated end is not
                // degraded on that side.
                let annotated = annotated_ends.get(tx_id);
                for (read_idx, p5, p3) in reads {
                    if (*p5 - mode_5p).abs() > tolerance {
                        if annotated.is_some_and(|(tss, _)| (*p5 - tss).abs() <= end_window) {
                            n_exempt_5p += 1;
                        } else {
                            outliers_5p.push((*read_idx, *p5));
                        }
                    }
                    if (*p3 - mode_3p).abs() > tolerance {
                        if annotated.is_some_and(|(_, tes)| (*p3 - tes).abs() <= end_window) {
                            n_exempt_3p += 1;
                        } else {
                            outliers_3p.push((*read_idx, *p3));
                        }
                    }
                }

                // ── Pass 2: rescue coherent sub-clusters (≥3 reads within 2 bp) ──
                let (rescued_set_5p, sc5, nr5) = rescue_coherent_subclusters(
                    &mut outliers_5p, SUBCLUSTER_GAP, MIN_SUBCLUSTER_READS,
                );
                let (rescued_set_3p, sc3, nr3) = rescue_coherent_subclusters(
                    &mut outliers_3p, SUBCLUSTER_GAP, MIN_SUBCLUSTER_READS,
                );
                n_subclusters_5p += sc5;
                n_subclusters_3p += sc3;
                n_rescued_5p += nr5;
                n_rescued_3p += nr3;

                // Flag the outlier reads that are in no coherent sub-cluster
                for (read_idx, _) in &outliers_5p {
                    if !rescued_set_5p.contains(read_idx) {
                        hints.is_degraded_5p.insert(*read_idx);
                        n_degraded_5p += 1;
                    }
                }
                for (read_idx, _) in &outliers_3p {
                    if !rescued_set_3p.contains(read_idx) {
                        hints.is_degraded_3p.insert(*read_idx);
                        n_degraded_3p += 1;
                    }
                }
            }

            eprintln!(
                "[IsoForge] Degradation filter (tolerance={}bp): {} transcripts scanned, \
                 {} skipped for having no modal end, \
                 {} reads flagged 5'-degraded, {} reads flagged 3'-degraded",
                tolerance, n_tx_processed, n_tx_no_mode, n_degraded_5p, n_degraded_3p
            );
            eprintln!(
                "[IsoForge]   Sub-cluster rescue (gap≤{}bp, min≥{} reads): \
                 {} 5'-subclusters rescued {} reads, {} 3'-subclusters rescued {} reads",
                SUBCLUSTER_GAP, MIN_SUBCLUSTER_READS,
                n_subclusters_5p, n_rescued_5p, n_subclusters_3p, n_rescued_3p
            );
            if config.known_ends_not_degraded {
                eprintln!(
                    "[IsoForge]   Annotated-end exemption (window {} bp): {} reads kept at 5', {} at 3'",
                    end_window, n_exempt_5p, n_exempt_3p
                );
            }

            // ── 5' dispersion filter ────────────────────────
            // Stage 2: for novel transcripts with dispersed 5' positions (IQR above
            //   threshold), reads whose 5' ends lie inside the body of a more abundant
            //   transcript of the same gene are flagged as 5'-degraded. Annotated
            //   transcripts are exempt.
            // Stage 3: ISM-suffix reads are grouped by gene; if too few lie within the
            //   window around the mode, all of them are flagged.

            const IQR_THRESHOLD: f64 = 50.0;          // Stage 2: min IQR (bp) to trigger
            const BODY_SLACK: i64 = 50;                // Stage 2: slack from parent ends
            const BODY_FRACTION: f64 = 0.5;            // Stage 2: min fraction inside body
            const MIN_READS_BODY: usize = 5;           // Stage 2: min reads to evaluate
            const DISPERSION_WINDOW: i64 = 5;          // Stage 3: ±bp around mode
            const DISPERSION_THRESHOLD: f64 = 0.30;    // Stage 3: min concentration
            const MIN_ISM_SUFFIX_READS: usize = 5;     // Stage 3: min reads to judge

            // ── Stage 2: Cross-transcript body-overlap (novel tx only) ──

            // Build transcript_id → gene_id mapping
            let mut tx_gene_id: HashMap<TranscriptId, GeneId> = HashMap::default();
            for annot in cores.iter() {
                tx_gene_id.entry(annot.transcript_id).or_insert(annot.gene_id);
            }

            // Annotated transcript IDs from the reference, exempt from Stage 2.
            let known_tx_ids: HashSet<TranscriptId> = self.refs.transcript_dict
                .values()
                .map(|m| m.transcript_id)
                .collect();

            // Per-transcript stats: (gene_id, mode_5p, mode_3p, iqr_5p, reads with 5' positions)
            let mut disp_stats: HashMap<TranscriptId, (GeneId, i64, i64, f64, Vec<(ReadIdx, i64)>)> =
                HashMap::default();
            let mut gene_tx_list: HashMap<GeneId, Vec<TranscriptId>> = HashMap::default();

            for (tx_id, reads) in &tx_reads {
                if reads.len() < MIN_READS_BODY { continue; }

                let gene_id = match tx_gene_id.get(tx_id) {
                    Some(&gid) => gid,
                    None => continue,
                };

                // Compute mode_5p and mode_3p
                let mut p5_counts: HashMap<i64, usize> = HashMap::default();
                let mut p3_counts: HashMap<i64, usize> = HashMap::default();
                for (_, p5, p3) in reads {
                    *p5_counts.entry(*p5).or_insert(0) += 1;
                    *p3_counts.entry(*p3).or_insert(0) += 1;
                }
                let m5 = p5_counts.iter()
                    .max_by(|a, b| a.1.cmp(b.1).then_with(|| b.0.cmp(a.0)))
                    .map(|(pos, _)| *pos).unwrap();
                let m3 = p3_counts.iter()
                    .max_by(|a, b| a.1.cmp(b.1).then_with(|| a.0.cmp(b.0)))
                    .map(|(pos, _)| *pos).unwrap();

                // Compute IQR of 5' positions
                let mut positions_5p: Vec<i64> = reads.iter().map(|(_, p5, _)| *p5).collect();
                positions_5p.sort_unstable();
                let n = positions_5p.len();
                let q1 = positions_5p[n / 4];
                let q3 = positions_5p[3 * n / 4];
                let iqr = (q3 - q1) as f64;

                let read_5p: Vec<(ReadIdx, i64)> = reads.iter()
                    .map(|(rid, p5, _)| (*rid, *p5)).collect();

                disp_stats.insert(*tx_id, (gene_id, m5, m3, iqr, read_5p));
                gene_tx_list.entry(gene_id).or_default().push(*tx_id);
            }

            let mut n_stage2_flagged = 0u64;
            let mut n_stage2_tx = 0u64;

            for (_gene_id, tx_ids) in &gene_tx_list {
                // Sort transcripts descending by read count (most abundant first)
                let mut sorted: Vec<_> = tx_ids.iter()
                    .filter_map(|tid| disp_stats.get(tid).map(|s| (*tid, s)))
                    .collect();
                sorted.sort_by(|a, b| b.1.4.len().cmp(&a.1.4.len()));

                for i in 0..sorted.len() {
                    let (tx_id, &(_, _mode_5p, _mode_3p, iqr, ref reads)) = sorted[i];

                    // Annotated transcripts are exempt from Stage 2
                    if known_tx_ids.contains(&tx_id) { continue; }

                    // IQR screening — skip transcripts with tight 5' distributions
                    if iqr <= IQR_THRESHOLD { continue; }

                    // Check if reads' 5' positions fall inside a more-abundant
                    // transcript's body (between its TSS+slack and TES-slack)
                    let mut inside_body: rustc_hash::FxHashSet<ReadIdx> = rustc_hash::FxHashSet::default();

                    for j in 0..i {
                        let (_, &(_, p_m5, p_m3, _, ref p_reads)) = sorted[j];
                        // Parent must have strictly more reads
                        if p_reads.len() <= reads.len() { continue; }

                        let p_low = p_m5.min(p_m3);
                        let p_high = p_m5.max(p_m3);
                        // Parents spanning less than 2 * BODY_SLACK + 100 bp are skipped
                        if p_high - p_low < 2 * BODY_SLACK + 100 { continue; }

                        for (read_idx, pos_5p) in reads {
                            if *pos_5p > p_low + BODY_SLACK && *pos_5p < p_high - BODY_SLACK {
                                inside_body.insert(*read_idx);
                            }
                        }
                    }

                    let body_frac = inside_body.len() as f64 / reads.len() as f64;

                    if body_frac >= BODY_FRACTION {
                        // Stage 2 triggered — flag reads inside parent body
                        for read_idx in &inside_body {
                            hints.is_degraded_5p.insert(*read_idx);
                        }
                        n_stage2_flagged += inside_body.len() as u64;
                        n_stage2_tx += 1;
                    }
                }
            }

            if n_stage2_flagged > 0 {
                eprintln!(
                    "[IsoForge] 5' Body-overlap filter (IQR>{:.0}bp, novel only): \
                     {} tx, {} reads flagged",
                    IQR_THRESHOLD, n_stage2_tx, n_stage2_flagged,
                );
            }

            // ── Stage 3: ISM-suffix 5' dispersion filter ────────────────
            // ISM-suffix reads are grouped by gene; if fewer than 30% lie within
            // ±5 bp of the mode, all of the gene's ISM-suffix reads are flagged.

            let mut ism_suffix_genes: HashMap<GeneId, Vec<(ReadIdx, i64)>> = HashMap::default();

            for (idx, annot) in cores.iter().enumerate() {
                let is_ism_suffix = (annot.hint_flags & HINT_IS_ISM != 0)
                    && (annot.hint_flags & HINT_IS_ISM_SUFFIX != 0);
                if !is_ism_suffix { continue; }

                let bam_5p = annot.start_pos;
                if bam_5p == 0 { continue; }

                ism_suffix_genes
                    .entry(annot.gene_id)
                    .or_default()
                    .push((idx as ReadIdx, bam_5p));
            }

            let mut n_stage3_flagged = 0u64;
            let mut n_stage3_genes = 0u64;

            for (_gene_id, reads) in &ism_suffix_genes {
                if reads.len() < MIN_ISM_SUFFIX_READS { continue; }

                let mut pos_counts: HashMap<i64, usize> = HashMap::default();
                for (_, p) in reads {
                    *pos_counts.entry(*p).or_insert(0) += 1;
                }
                let mode_pos = pos_counts
                    .iter()
                    .max_by_key(|(_, &c)| c)
                    .map(|(&pos, _)| pos)
                    .unwrap();

                let near_mode = reads.iter()
                    .filter(|(_, p)| (*p - mode_pos).abs() <= DISPERSION_WINDOW)
                    .count();
                let concentration = near_mode as f64 / reads.len() as f64;

                if concentration < DISPERSION_THRESHOLD {
                    n_stage3_genes += 1;
                    for (read_idx, _) in reads {
                        hints.is_degraded_5p.insert(*read_idx);
                        n_stage3_flagged += 1;
                    }
                }
            }

            if n_stage2_flagged > 0 || n_stage3_flagged > 0 {
                eprintln!(
                    "[IsoForge] 5' Dispersion filter: \
                     Stage 2 (body-overlap): {} tx, {} reads | \
                     Stage 3 (ISM-Suffix): {} genes, {} reads",
                    n_stage2_tx, n_stage2_flagged,
                    n_stage3_genes, n_stage3_flagged,
                );
            }

            eprintln!(
                "[IsoForge] Total degradation hints: {} is_degraded_5p, {} is_degraded_3p",
                hints.is_degraded_5p.len(), hints.is_degraded_3p.len()
            );
        }


        hints
    }

    /// Split transcripts by TSS/TES cluster assignment.
    ///
    /// Reads of one transcript that join different TSS or TES clusters form
    /// separate variants with their own IDs, provided the group's exon structure
    /// is consistent with the parent. Pass 1 counts reads per
    /// (tx, tss_cid, tes_cid, dataset); the decision step allocates IDs; pass 2
    /// relabels reads in place.
    ///
    /// Returns the number of new transcript variants created.
    fn split_transcripts_by_clusters(
        &self,
        cores: &mut [ObsCore],
        keep: &[bool],
        join_results: &mut [ReadEndSiteInfo],
        results: &mut IsoForgeResults,
        tx_merges: &mut HashMap<TranscriptId, TranscriptId>,
    ) -> usize {
        let (exclude_upstream, follow_ends, join_distance, contradiction, polya_required) = self.cluster_config.as_ref()
            .map(|c| (c.exclude_upstream_unclustered_5p, c.annotated_id_follows_ends, c.tss_tes_distance,
                      c.end_contradiction_distance, c.require_polya_tail))
            .unwrap_or((false, false, 50, 100, true));

        // ── Upstream reads with no TSS cluster ──────────────────────────────
        // A read whose 5' end joined no TSS cluster and lies upstream of its
        // transcript's TSS is marked end-unsupported, removed from the count and
        // kept out of the grouping.
        let mut upstream_unclustered: HashSet<usize> = HashSet::default();
        if exclude_upstream {
            // Upstream is measured from the further upstream of the main TSS cluster
            // and, for annotated transcripts, the annotated TSS.
            let with_reads: HashSet<TranscriptId> = cores.iter().enumerate()
                .filter(|(i, _)| keep.get(*i).copied().unwrap_or(false))
                .map(|(_, c)| c.transcript_id).collect();
            let positions: rustc_hash::FxHashMap<i64, i64> = self.refs.location_dict.chromosomes.values()
                .flat_map(|c| c.positions.iter().map(|&(p, v)| (v, p))).collect();
            let annotated_tss: HashMap<TranscriptId, i64> = self.refs.transcript_dict.values()
                .filter(|t| with_reads.contains(&t.transcript_id))
                .filter_map(|t| Some((t.transcript_id, *positions.get(&t.start_vertex)?)))
                .collect();
            drop(positions);
            let mut tss_votes: HashMap<TranscriptId, HashMap<i64, usize>> = HashMap::default();
            for (idx, annot) in cores.iter().enumerate() {
                if !keep.get(idx).copied().unwrap_or(false) || self.is_reassigned_ism(annot) { continue; }
                if let Some(site) = join_results.get(idx).and_then(|e| e.tss_site) {
                    *tss_votes.entry(annot.transcript_id).or_default().entry(site).or_insert(0) += 1;
                }
            }
            let main_tss: HashMap<TranscriptId, i64> = tss_votes.into_iter()
                .filter_map(|(t, votes)| votes.into_iter()
                    .max_by(|a, b| a.1.cmp(&b.1).then(b.0.cmp(&a.0)))
                    .map(|(site, _)| (t, site)))
                .collect();
            for (idx, annot) in cores.iter().enumerate() {
                if !keep.get(idx).copied().unwrap_or(false) || self.is_reassigned_ism(annot) { continue; }
                if join_results.get(idx).map_or(true, |e| e.tss_cluster_id.is_some()) { continue; }
                let Some(&main) = main_tss.get(&annot.transcript_id) else { continue };
                let minus = matches!(annot.strand, crate::types::Strand::Minus);
                let tss = match annotated_tss.get(&annot.transcript_id) {
                    Some(&a) if minus => main.max(a),
                    Some(&a) => main.min(a),
                    None => main,
                };
                let five = if minus { annot.start_pos.max(annot.end_pos) } else { annot.start_pos.min(annot.end_pos) };
                let upstream = if minus { five - tss } else { tss - five };
                if upstream > join_distance {
                    upstream_unclustered.insert(idx);
                }
            }
            for &idx in &upstream_unclustered {
                if let Some(info) = join_results.get_mut(idx) {
                    info.is_supported = false;
                    info.rejection_reason = Some("upstream_unclustered_5p".to_string());
                }
                let annot = &cores[idx];
                let key = (annot.gene_id, annot.transcript_id, Arc::clone(&annot.dataset));
                if let Some(n) = results.abundance.get_mut(&key) {
                    *n -= 1;
                    if *n <= 0 { results.abundance.remove(&key); }
                }
            }
            eprintln!(
                "[IsoForge] Upstream unclustered 5' ends: {} reads start more than {} bp upstream \
                 of their transcript's reference TSS (main cluster or, if further upstream, the annotated TSS) \
                 and joined no TSS cluster; kept out of its count",
                upstream_unclustered.len(), join_distance
            );
        }

        // ── Pass 1: count reads per (tx_id, tss_cid, tes_cid, dataset) ──────
        struct GroupInfo {
            total: usize,
            ds_counts: HashMap<Arc<str>, i64>,
            /// Reads per exon count, for the structure guard.
            n_exons_counts: HashMap<i32, usize>,
            /// Representative TSS cluster peak position, used for proximity absorption.
            tss_site: Option<i64>,
            /// Representative TES cluster peak position.
            polya_site: Option<i64>,
            /// Reads with a poly(A) tail, for the full-length start guard.
            tailed: usize,
            /// Reads with a cap signature, for keeping a capped start near the dominant one apart.
            capped: usize,
        }

        let mut tx_groups: HashMap<
            TranscriptId,
            HashMap<(Option<i64>, Option<i64>), GroupInfo>,
        > = HashMap::default();
        let small_rna = self.refs.small_rna.as_deref();

        for (idx, annot) in cores.iter().enumerate() {
            if !keep.get(idx).copied().unwrap_or(false) {
                continue;
            }
            // A truncated read never seeds an end variant of its parent.
            if self.is_reassigned_ism(annot) {
                continue;
            }
            // Nor is an end variant split from an annotated small RNA (small_rna.rs).
            if small_rna.is_some_and(|index| index.is_small_rna(annot.transcript_id)) {
                continue;
            }
            if upstream_unclustered.contains(&idx) {
                continue;
            }
            let info = join_results.get(idx);
            let tss_cid = info.and_then(|e| e.tss_cluster_id);
            let tes_cid = info.and_then(|e| e.tes_cluster_id);

            let group = tx_groups
                .entry(annot.transcript_id)
                .or_default()
                .entry((tss_cid, tes_cid))
                .or_insert_with(|| GroupInfo {
                    total: 0,
                    ds_counts: HashMap::default(),
                    n_exons_counts: HashMap::default(),
                    tss_site: None,
                    polya_site: None,
                    tailed: 0,
                    capped: 0,
                });
            group.total += 1;
            if annot.cap_signature {
                group.capped += 1;
            }
            if annot.polya_tail_len.map_or(false, |v| v > 0) || annot.polya_tail_softclip.map_or(false, |v| v > 0) {
                group.tailed += 1;
            }
            *group.ds_counts.entry(annot.dataset.clone()).or_insert(0) += 1;
            *group.n_exons_counts
                .entry(annot.n_exons)
                .or_insert(0) += 1;
            // Representative peak position (identical for all reads in a cluster)
            if group.tss_site.is_none() {
                group.tss_site = info.and_then(|e| e.tss_site);
            }
            if group.polya_site.is_none() {
                group.polya_site = info.and_then(|e| e.polya_site);
            }
        }

        // ── Cross-parent splice merge ───────────────────────────────────────
        // Novel transcripts of the same gene with identical internal splice
        // junctions but different terminal exons are merged under one parent.
        {
            // tx_id → (gene_id, jn_path), novel transcripts only.
            let known_tx_ids: HashSet<TranscriptId> = self
                .refs
                .transcript_dict
                .values()
                .map(|m| m.transcript_id)
                .collect();

            let mut tx_splice_sig: HashMap<TranscriptId, (GeneId, Option<String>)> =
                HashMap::default();
            for m in results.novel_transcripts.iter() {
                // Multi-exon novel transcripts only.
                if tx_groups.contains_key(&m.transcript_id)
                    && !known_tx_ids.contains(&m.transcript_id)
                    && m.jn_path.is_some()
                    && m.n_exons > 1
                {
                    tx_splice_sig
                        .insert(m.transcript_id, (m.gene_id, m.jn_path.clone()));
                }
            }

            // A novel single-exon transcript is keyed on the TSS and TES clusters of its largest
            // group that has both.
            for m in results.novel_transcripts.iter() {
                if m.n_exons != 1 || known_tx_ids.contains(&m.transcript_id) { continue; }
                let Some(groups) = tx_groups.get(&m.transcript_id) else { continue };
                let Some(((tss, tes), _)) = groups.iter()
                    .filter(|(k, _)| k.0.is_some() && k.1.is_some())
                    .max_by_key(|(_, v)| v.total) else { continue };
                tx_splice_sig.insert(m.transcript_id,
                    (m.gene_id, Some(format!("mono:{}:{}", tss.unwrap(), tes.unwrap()))));
            }

            // Group tx_ids by (gene_id, jn_path)
            let mut splice_groups: HashMap<(GeneId, Option<String>), Vec<TranscriptId>> =
                HashMap::default();
            for (tx_id, (gene_id, jn_path)) in &tx_splice_sig {
                splice_groups
                    .entry((*gene_id, jn_path.clone()))
                    .or_default()
                    .push(*tx_id);
            }

            // Merge minor tx_ids of each group into the one with the most reads.
            let mut n_cross_parent_merged = 0usize;
            let mut n_cross_parent_reads = 0usize;
            let mut remap_tx: HashMap<TranscriptId, TranscriptId> = HashMap::default();
            // Equal read totals go to the model with the lower ends.
            let tx_ends: HashMap<TranscriptId, (i64, i64)> = results.novel_transcripts.iter()
                .map(|m| (m.transcript_id, (m.start_pos, m.end_pos))).collect();

            for (_sig, tx_ids) in &splice_groups {
                if tx_ids.len() <= 1 {
                    continue;
                }

                // Dominant tx_id: most reads across all its cluster groups
                let mut tx_totals: Vec<(TranscriptId, usize)> = tx_ids
                    .iter()
                    .map(|&tid| {
                        let total: usize = tx_groups
                            .get(&tid)
                            .map(|g| g.values().map(|v| v.total).sum())
                            .unwrap_or(0);
                        (tid, total)
                    })
                    .collect();
                tx_totals.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| {
                    (tx_ends.get(&a.0), a.0).cmp(&(tx_ends.get(&b.0), b.0))
                }));

                let dominant_tx = tx_totals[0].0;

                // Move non-dominant groups into the dominant tx_id
                for &(minor_tx, minor_total) in &tx_totals[1..] {
                    if minor_total == 0 {
                        continue;
                    }

                    if let Some(minor_groups) = tx_groups.remove(&minor_tx) {
                        let dom_groups = tx_groups.entry(dominant_tx).or_default();
                        for (cluster_key, minor_info) in minor_groups {
                            let dom_entry = dom_groups
                                .entry(cluster_key)
                                .or_insert_with(|| GroupInfo {
                                    total: 0,
                                    ds_counts: HashMap::default(),
                                    n_exons_counts: HashMap::default(),
                                    tss_site: minor_info.tss_site,
                                    polya_site: minor_info.polya_site,
                                tailed: 0,
                                capped: 0,
                                });
                            dom_entry.total += minor_info.total;
                            for (ds, cnt) in minor_info.ds_counts {
                                *dom_entry.ds_counts.entry(ds).or_insert(0) += cnt;
                            }
                            dom_entry.tailed += minor_info.tailed;
                            dom_entry.capped += minor_info.capped;
                            for (ne, cnt) in minor_info.n_exons_counts {
                                *dom_entry.n_exons_counts.entry(ne).or_insert(0) += cnt;
                            }
                        }
                        remap_tx.insert(minor_tx, dominant_tx);
                        tx_merges.insert(minor_tx, dominant_tx);
                        n_cross_parent_merged += 1;
                        n_cross_parent_reads += minor_total;
                    }
                }
            }

            // Remap reads from minor to dominant transcript_id
            if !remap_tx.is_empty() {
                for annot in cores.iter_mut() {
                    if let Some(&new_tx) = remap_tx.get(&annot.transcript_id) {
                        annot.transcript_id = new_tx;
                    }
                }

                // Remap abundance entries
                let old_keys: Vec<_> = results
                    .abundance
                    .keys()
                    .filter(|(_, tid, _)| remap_tx.contains_key(tid))
                    .cloned()
                    .collect();
                for (gid, old_tid, ds) in old_keys {
                    if let Some(&new_tid) = remap_tx.get(&old_tid) {
                        let count = results
                            .abundance
                            .remove(&(gid, old_tid, ds.clone()))
                            .unwrap_or(0);
                        *results.abundance.entry((gid, new_tid, ds)).or_insert(0) +=
                            count;
                    }
                }

                // Remove remapped novel transcripts from the list
                let remap_set: HashSet<TranscriptId> =
                    remap_tx.keys().copied().collect();
                results
                    .novel_transcripts
                    .retain(|m| !remap_set.contains(&m.transcript_id));

                // Remove novelty annotations for remapped transcripts
                results
                    .transcript_annotations
                    .retain(|a| !remap_set.contains(&a.id));

                eprintln!(
                    "[IsoForge] Cross-parent splice merge: {} minor transcript(s) \
                     ({} reads) merged into dominant parents (same gene + internal \
                     splice junctions, different first/last exon endpoints)",
                    n_cross_parent_merged, n_cross_parent_reads
                );
            }
        }

        // Annotated (TSS, TES, gene) of each reference transcript that has reads
        // here, in transcript orientation, for deciding which group keeps its ID.
        let annotated_ends: HashMap<TranscriptId, (i64, i64, GeneId)> = if follow_ends {
            let positions: rustc_hash::FxHashMap<i64, i64> = self.refs.location_dict.chromosomes.values()
                .flat_map(|c| c.positions.iter().map(|&(p, v)| (v, p))).collect();
            self.refs.transcript_dict.values()
                .filter(|t| tx_groups.contains_key(&t.transcript_id))
                .filter_map(|t| Some((t.transcript_id, (*positions.get(&t.start_vertex)?, *positions.get(&t.end_vertex)?, t.gene_id))))
                .collect()
        } else {
            HashMap::default()
        };
        // Distinguishing terminal sequence. When every other annotated transcript of
        // the gene that contains T's whole intron chain has an intron just before T's
        // first donor, and T's first exon extends into that intron, a read group must
        // start beyond that exon start (by more than the join distance) to carry T's
        // ID; likewise at the 3' end. Applies only when the containing transcripts
        // carry at least DISTINGUISHING_PARENT_RATIO times T's reads. Maps T to its
        // (5' boundary, 3' boundary) in genomic coordinates.
        let distinguishing: HashMap<TranscriptId, (Option<i64>, Option<i64>)> = if follow_ends {
            let genes: HashSet<GeneId> = annotated_ends.values().map(|e| e.2).collect();
            let mut by_gene: HashMap<GeneId, Vec<(TranscriptId, Vec<i64>, Strand)>> = HashMap::default();
            for m in self.refs.transcript_dict.values() {
                if !genes.contains(&m.gene_id) { continue; }
                if let Some((chain, _, _, strand)) = self.refs.hypergraph.internal(m.transcript_id) {
                    by_gene.entry(m.gene_id).or_default().push((m.transcript_id, chain, strand));
                }
            }
            let reads_of = |t: &TranscriptId| -> usize {
                tx_groups.get(t).map_or(0, |g| g.values().map(|v| v.total).sum())
            };
            let mut out = HashMap::default();
            for (&tx, &(a_tss, a_tes, gene)) in &annotated_ends {
                let Some(peers) = by_gene.get(&gene) else { continue };
                let Some((_, own, strand)) = peers.iter().find(|p| p.0 == tx) else { continue };
                let up = |x: i64, y: i64| if *strand == Strand::Minus { x - y } else { y - x }; // x upstream of y by
                let (mut b5, mut b3, mut five_ok, mut three_ok) = (None::<i64>, None::<i64>, true, true);
                let (mut five_seen, mut three_seen) = (false, false);
                let mut containing_reads = 0usize;
                for (other, chain, st) in peers {
                    if *other == tx || st != strand || chain.len() <= own.len() { continue; }
                    let Some(j) = (0..=(chain.len() - own.len()) / 2).map(|j| 2 * j).find(|&j| chain[j..j + own.len()] == own[..]) else { continue };
                    containing_reads += reads_of(other);
                    if j > 0 {
                        five_seen = true;
                        let exon_start = chain[j - 1];
                        if up(a_tss, exon_start) > 0 {
                            if b5.map_or(true, |b| up(exon_start, b) > 0) { b5 = Some(exon_start); }
                        } else {
                            five_ok = false;
                        }
                    }
                    if j + own.len() < chain.len() {
                        three_seen = true;
                        let exon_end = chain[j + own.len()];
                        if up(exon_end, a_tes) > 0 {
                            if b3.map_or(true, |b| up(b, exon_end) > 0) { b3 = Some(exon_end); }
                        } else {
                            three_ok = false;
                        }
                    }
                }
                if containing_reads < DISTINGUISHING_PARENT_RATIO * reads_of(&tx).max(1) { continue; }
                let b5 = if five_seen && five_ok { b5 } else { None };
                let b3 = if three_seen && three_ok { b3 } else { None };
                if b5.is_some() || b3.is_some() {
                    out.insert(tx, (b5, b3));
                }
            }
            out
        } else {
            HashMap::default()
        };
        // A read group is compatible with an annotated transcript unless a TSS
        // cluster lies upstream of the annotated TSS, or a TES cluster downstream of
        // the annotated TES, by more than end_contradiction_distance. An end with no
        // cluster does not contradict.
        let compatible = |tx: &TranscriptId, tss: Option<i64>, tes: Option<i64>| -> bool {
            let Some(&(a_tss, a_tes, _)) = annotated_ends.get(tx) else { return true };
            let minus = a_tss > a_tes;
            let (b5, b3) = distinguishing.get(tx).copied().unwrap_or((None, None));
            if let Some(s) = tss {
                let upstream = if minus { s - a_tss } else { a_tss - s };
                if upstream > contradiction { return false; }
                // The group must reach into T's distinguishing 5' sequence.
                if b5.is_some_and(|b| (if minus { s - b } else { b - s }) <= join_distance) { return false; }
            }
            if let Some(e) = tes {
                let downstream = if minus { a_tes - e } else { e - a_tes };
                if downstream > contradiction { return false; }
                // When the poly(A) tail is required, a TES cluster upstream of the annotated TES by more
                // than end_contradiction_distance contradicts as well.
                if polya_required && downstream < -contradiction { return false; }
                if b3.is_some_and(|b| (if minus { b - e } else { e - b }) <= join_distance) { return false; }
            }
            true
        };
        let mut id_swaps: Vec<(TranscriptId, TranscriptId)> = Vec::new();
        let mut variant_parents: Vec<(TranscriptId, TranscriptId, bool)> = Vec::new();
        // (parent, dominant group key, reads) of annotated transcripts none of whose
        // groups is compatible: the annotated ID is carried by no group.
        let mut id_demotes: Vec<(TranscriptId, (Option<i64>, Option<i64>), usize)> = Vec::new();
        // Group keys and read totals per annotated transcript, for merge targets.
        let mut annotated_group_keys: HashMap<TranscriptId, Vec<((Option<i64>, Option<i64>), usize, bool)>> = HashMap::default();

        // ── Decision: allocate new IDs for non-dominant groups ───────────────
        // combo_to_new_tx: (tx_id, tss_cid, tes_cid) → new_tx_id
        let mut combo_to_new_tx: HashMap<
            (TranscriptId, Option<i64>, Option<i64>),
            TranscriptId,
        > = HashMap::default();

        // Collect (new_tx_id, parent_tx_id, n_reads, ds_counts)
        struct VariantInfo {
            new_tx_id: TranscriptId,
            parent_tx_id: TranscriptId,
            n_reads: usize,
            ds_counts: HashMap<Arc<str>, i64>,
        }
        let mut variants: Vec<VariantInfo> = Vec::new();
        let mut n_rejected_structure = 0usize;
        let mut n_absorbed_none_cluster = 0usize;
        let mut n_absorbed_minor_cluster = 0usize;
        let mut n_absorbed_annotated_ends = 0usize;
        let mut n_absorbed_singleton_variant = 0usize;
        let mut n_absorbed_untailed_start = 0usize;
        let mut n_joined_shared_end = 0usize;

        // Cluster peaks within this distance of the dominant's peaks are the same
        // endpoint for proximity absorption.
        let proximity_threshold = self.cluster_config
            .as_ref()
            .map(|c| c.tss_tes_distance)
            .unwrap_or(50);
        let split_tss = self.cluster_config.as_ref().filter(|c| c.split_capped_tss).cloned();
        let mut n_split_capped_tss = 0usize;

        for (tx_id, cluster_groups) in &tx_groups {
            // Dominant group (most reads) keeps original transcript_id
            let mut groups_vec: Vec<_> = cluster_groups.iter().collect();
            // Most reads first; a tie goes to more tailed reads, then to the ends nearer the
            // annotated ends, then to the cluster order.
            let off_annotation = |g: &GroupInfo| -> i64 {
                annotated_ends.get(tx_id).map_or(0, |&(a_tss, a_tes, _)| {
                    g.tss_site.map_or(1 << 40, |s| (s - a_tss).abs()) + g.polya_site.map_or(1 << 40, |s| (s - a_tes).abs())
                })
            };
            groups_vec.sort_by(|a, b| b.1.total.cmp(&a.1.total)
                .then(b.1.tailed.cmp(&a.1.tailed))
                .then(off_annotation(a.1).cmp(&off_annotation(b.1)))
                .then(a.0.cmp(b.0)));

            // Get the dominant group's cluster IDs and peak positions
            let (dom_tss, dom_tes) = groups_vec[0].0;
            let dom_tss_site = groups_vec[0].1.tss_site;
            let dom_polya_site = groups_vec[0].1.polya_site;

            if annotated_ends.contains_key(tx_id) {
                annotated_group_keys.insert(*tx_id, groups_vec.iter().enumerate()
                    .map(|(i, (k, v))| (**k, v.total, i == 0)).collect());
            }
            if cluster_groups.len() <= 1 {
                // A single TSS/TES combination is not split, but its reads may still
                // contradict the annotated ends.
                if annotated_ends.contains_key(tx_id) && !compatible(tx_id, dom_tss_site, dom_polya_site) {
                    id_demotes.push((*tx_id, (*dom_tss, *dom_tes), groups_vec[0].1.total));
                }
                continue;
            }
            let dom_n_exons = groups_vec[0].1.n_exons_counts.iter()
                .max_by_key(|&(_, &count)| count)
                .map(|(&n, _)| n);
            let tx_total: usize = groups_vec.iter().map(|(_, g)| g.total).sum();

            let mut created_here: Vec<(TranscriptId, Option<i64>, Option<i64>, usize)> = Vec::new();
            // Cluster keys of the variants made here, with their index in `variants` and `created_here`.
            let mut created_keys: Vec<((Option<i64>, Option<i64>), TranscriptId, usize, usize)> = Vec::new();
            // Groups with both end clusters are decided first, then those lacking one; the
            // largest first within each.
            let mut rest: Vec<_> = groups_vec.iter().skip(1).collect();
            rest.sort_by_key(|((t, e), info)| (t.is_none() || e.is_none(), std::cmp::Reverse(info.total), *t, *e));
            for ((tss_cid, tes_cid), info) in rest {
                // A group missing one end cluster whose other cluster is an end variant's,
                // not the dominant group's, belongs to that variant.
                let shared = shared_end_variant((*tss_cid, *tes_cid), (*dom_tss, *dom_tes),
                                                created_keys.iter().map(|k| k.0))
                    .map(|i| created_keys[i]);
                if let Some((_, new_tx_id, vi, ci)) = shared {
                    combo_to_new_tx.insert((*tx_id, *tss_cid, *tes_cid), new_tx_id);
                    variants[vi].n_reads += info.total;
                    for (ds, n) in &info.ds_counts {
                        *variants[vi].ds_counts.entry(ds.clone()).or_insert(0) += n;
                    }
                    created_here[ci].3 += info.total;
                    n_joined_shared_end += info.total;
                    continue;
                }
                // A single-read group is absorbed into the dominant group rather than split.
                if info.total < MIN_END_VARIANT_READS {
                    n_absorbed_singleton_variant += info.total;
                    continue;
                }

                // ── None-cluster absorption ────────────────────────────
                // A group whose cluster IDs are None where the dominant group's are Some,
                // and equal where both are Some, is absorbed into the dominant group:
                //   (None, TES_Y), (TSS_X, None), (None, None) into (TSS_X, TES_Y).
                // (TSS_Z, TES_Y) and (None, TES_Z) remain variants.
                let tss_compatible = tss_cid.is_none() || *tss_cid == *dom_tss;
                let tes_compatible = tes_cid.is_none() || *tes_cid == *dom_tes;
                if tss_compatible && tes_compatible {
                    n_absorbed_none_cluster += info.total;
                    log::debug!(
                        "[IsoForge] Absorbed {} reads into dominant group of tx {} \
                         (None-cluster absorption: group=({:?},{:?}), dominant=({:?},{:?}))",
                        info.total, tx_id, tss_cid, tes_cid, dom_tss, dom_tes
                    );
                    continue;
                }

                // ── Minor-cluster proximity absorption ──────────────────────
                // A group with the dominant's exon count whose cluster peaks lie within the
                // proximity threshold of the dominant's is absorbed.
                {
                    let tss_close = match (info.tss_site, dom_tss_site) {
                        (Some(a), Some(b)) => (a - b).abs() <= proximity_threshold,
                        _ => true, // missing peak → treat as compatible
                    };
                    let tes_close = match (info.polya_site, dom_polya_site) {
                        (Some(a), Some(b)) => (a - b).abs() <= proximity_threshold,
                        _ => true,
                    };
                    let group_dom_n_exons = info.n_exons_counts.iter()
                        .max_by_key(|&(_, &count)| count)
                        .map(|(&n, _)| n);

                    // A capped start of its own near the dominant start stays apart.
                    let distinct_start = split_tss.as_ref().is_some_and(|c| {
                        let (Some(a), Some(b)) = (info.tss_site, dom_tss_site) else { return false };
                        let dom = &groups_vec[0].1;
                        let dom_rate = dom.capped as f64 / dom.total.max(1) as f64;
                        (a - b).abs() >= c.split_tss_min_distance
                            && *tes_cid == *dom_tes && tes_cid.is_some()
                            && info.total >= c.split_tss_min_reads
                            && info.total as f64 >= c.split_tss_min_share * tx_total as f64
                            && dom_rate >= 0.05
                            && info.capped as f64 / info.total as f64 >= c.split_tss_cap_ratio * dom_rate
                    });
                    if distinct_start { n_split_capped_tss += info.total; }

                    if tss_close && tes_close
                        && dom_n_exons.is_some()
                        && group_dom_n_exons == dom_n_exons
                        && !distinct_start
                    {
                        n_absorbed_minor_cluster += info.total;
                        log::debug!(
                            "[IsoForge] Minor-cluster absorption: {} reads of tx {} \
                             (group peaks=({:?},{:?}), dom peaks=({:?},{:?}), dist_thresh={}bp)",
                            info.total, tx_id, info.tss_site, info.polya_site,
                            dom_tss_site, dom_polya_site, proximity_threshold
                        );
                        continue;
                    }

                    // ── Annotated-end absorption ─────────────────────────────
                    // Groups of an annotated transcript whose end peaks all lie within the join
                    // distance of the annotated TSS and TES are merged. When the largest group
                    // does not match the annotation, the matching group stays a variant.
                    if let Some(&(ann_tss, ann_tes, _)) = annotated_ends.get(tx_id) {
                        let at_annotation = |tss: Option<i64>, tes: Option<i64>| {
                            tss.map_or(true, |s| (s - ann_tss).abs() <= join_distance)
                                && tes.map_or(true, |s| (s - ann_tes).abs() <= join_distance)
                        };
                        if dom_n_exons.is_some()
                            && group_dom_n_exons == dom_n_exons
                            && at_annotation(info.tss_site, info.polya_site)
                            && at_annotation(dom_tss_site, dom_polya_site)
                            && !distinct_start
                        {
                            n_absorbed_annotated_ends += info.total;
                            continue;
                        }
                    }
                }

                // ── Structure guard ──────────────────────────────────────────
                // The dominant exon count must account for ≥60% of the group's reads.
                let dominant_structure_count = info.n_exons_counts
                    .values()
                    .max()
                    .copied()
                    .unwrap_or(0);
                let structure_consistency = dominant_structure_count as f64 / info.total as f64;
                if structure_consistency < 0.6 {
                    n_rejected_structure += 1;
                    log::debug!(
                        "[IsoForge] Structure guard rejected variant of tx {} ({} reads, {:.0}% n_exons consistency, {} distinct n_exons values)",
                        tx_id, info.total, structure_consistency * 100.0,
                        info.n_exons_counts.len()
                    );
                    continue;
                }

                // ── Full-length start guard ─────────────────────────────────
                // A spliced group with a start of its own is absorbed when starts come from full-length reads and it
                // has fewer than MIN_FULL_LENGTH_START_READS tailed reads; an unspliced one when the run has tails.
                let group_unspliced = info.n_exons_counts.iter().max_by_key(|&(_, &count)| count).map(|(&n, _)| n) == Some(1);
                let own_start = *tss_cid != *dom_tss;
                if own_start
                    && ((group_unspliced && self.tails_present.load(std::sync::atomic::Ordering::Relaxed))
                        || (!group_unspliced && self.full_length_starts.load(std::sync::atomic::Ordering::Relaxed) && info.tailed < MIN_FULL_LENGTH_START_READS))
                {
                    n_absorbed_untailed_start += info.total;
                    continue;
                }
                let new_tx_id = self.refs.counters.transcript.increment();
                combo_to_new_tx.insert((*tx_id, *tss_cid, *tes_cid), new_tx_id);
                created_keys.push(((*tss_cid, *tes_cid), new_tx_id, variants.len(), created_here.len()));
                variants.push(VariantInfo {
                    new_tx_id,
                    parent_tx_id: *tx_id,
                    n_reads: info.total,
                    ds_counts: info.ds_counts.clone(),
                });
                created_here.push((new_tx_id, info.tss_site, info.polya_site, info.total));
            }

            // The annotated ID goes to the group whose cluster peaks best match the
            // annotated TSS (weight 2) and TES (weight 1) within the join distance; if a
            // variant outscores the largest group, the two exchange IDs.
            if let Some(&(ann_tss, ann_tes, _)) = annotated_ends.get(tx_id) {
                let score = |tss: Option<i64>, tes: Option<i64>| -> u8 {
                    2 * u8::from(tss.is_some_and(|s| (s - ann_tss).abs() <= join_distance))
                        + u8::from(tes.is_some_and(|s| (s - ann_tes).abs() <= join_distance))
                };
                let dom_ok = compatible(tx_id, dom_tss_site, dom_polya_site);
                let dom_score = score(dom_tss_site, dom_polya_site);
                let best = created_here.iter()
                    .filter(|v| compatible(tx_id, v.1, v.2))
                    .max_by(|a, b| score(a.1, a.2).cmp(&score(b.1, b.2)).then(a.3.cmp(&b.3)));
                let swapped = match best {
                    Some(b) if !dom_ok || score(b.1, b.2) > dom_score => { id_swaps.push((*tx_id, b.0)); Some(b.0) }
                    None if !dom_ok => { id_demotes.push((*tx_id, (*dom_tss, *dom_tes), groups_vec[0].1.total)); None }
                    _ => None,
                };
                // A swap leaves the new ID on the largest group's reads.
                for v in &created_here {
                    let ok = if swapped == Some(v.0) { dom_ok } else { compatible(tx_id, v.1, v.2) };
                    variant_parents.push((v.0, *tx_id, ok));
                }
            }
        }

        if n_rejected_structure > 0 {
            eprintln!(
                "[IsoForge] Structure guard: rejected {} inconsistent variant groups",
                n_rejected_structure
            );
        }

        if n_joined_shared_end > 0 {
            eprintln!(
                "[IsoForge] Shared-end joins: {} reads without a cluster at one end joined the end variant \
                 sharing their cluster at the other",
                n_joined_shared_end
            );
        }

        if n_absorbed_none_cluster > 0 {
            eprintln!(
                "[IsoForge] None-cluster absorption: {} reads absorbed into dominant groups \
                 (same structure, missing TSS/TES cluster assignment)",
                n_absorbed_none_cluster
            );
        }

        if n_absorbed_annotated_ends > 0 {
            eprintln!(
                "[IsoForge] Annotated-end absorption: {} reads absorbed into annotated transcripts \
                 (same structure, every end peak within {}bp of the annotated TSS/TES)",
                n_absorbed_annotated_ends, join_distance
            );
        }

        if n_absorbed_minor_cluster > 0 {
            eprintln!(
                "[IsoForge] Minor-cluster proximity absorption: {} reads absorbed into dominant groups \
                 (same structure, adjacent TSS/TES cluster peaks within {}bp)",
                n_absorbed_minor_cluster, proximity_threshold
            );
        }
        if n_split_capped_tss > 0 {
            eprintln!(
                "[IsoForge] Capped starts near the dominant start kept apart: {} reads (a distinct capped TSS peak within {}bp)",
                n_split_capped_tss, proximity_threshold
            );
        }

        if n_absorbed_singleton_variant > 0 {
            eprintln!(
                "[IsoForge] Singleton end-variant absorption: {} reads kept with their dominant \
                 group (one read is not an alternative end pair)",
                n_absorbed_singleton_variant
            );
        }

        // Free the grouping map before pass 2
        drop(tx_groups);

        // ── Annotated IDs no read group is compatible with ────────────────
        // Their reads join a group of another annotated transcript of the same gene
        // on the same TSS and TES clusters whose intron chain is within the
        // splice-snapping distance at every junction; otherwise they become an end
        // variant with a new ID.
        let mut demote_merge: Vec<(TranscriptId, TranscriptId)> = Vec::new();
        let mut demote_new: Vec<(TranscriptId, TranscriptId)> = Vec::new();
        if !id_demotes.is_empty() {
            let models_by_id: HashMap<TranscriptId, &TranscriptModel> = self.refs.transcript_dict.values()
                .filter(|m| annotated_ends.contains_key(&m.transcript_id))
                .map(|m| (m.transcript_id, m)).collect();
            let mut by_gene: HashMap<GeneId, Vec<TranscriptId>> = HashMap::default();
            for (&t, e) in &annotated_ends { by_gene.entry(e.2).or_default().push(t); }
            let mut chains: HashMap<TranscriptId, Option<Vec<i64>>> = HashMap::default();
            let mut chain_of = |tx: TranscriptId| -> Option<Vec<i64>> {
                chains.entry(tx).or_insert_with(|| models_by_id.get(&tx)
                    .and_then(|m| self.refs.hypergraph.internal(m.transcript_id))
                    .map(|c| c.0)).clone()
            };
            let demoted: HashSet<TranscriptId> = id_demotes.iter().map(|d| d.0).collect();
            let swapped: HashMap<TranscriptId, TranscriptId> = id_swaps.iter().copied().collect();
            let tol = self.config.max_splice_dist.max(0);
            for &(parent, key, _) in &id_demotes {
                let gene = annotated_ends.get(&parent).map(|e| e.2);
                let target = if key.0.is_some() && key.1.is_some() {
                    let own = chain_of(parent);
                    let peers: Vec<TranscriptId> = gene.and_then(|g| by_gene.get(&g)).cloned().unwrap_or_default();
                    peers.iter()
                        .filter(|q| **q != parent && !demoted.contains(q))
                        .filter_map(|q| {
                            let keys = annotated_group_keys.get(q)?;
                            let &(_, total, is_dom) = keys.iter().find(|(k, _, _)| *k == key)?;
                            let (a, b) = (own.as_ref()?, chain_of(*q)?);
                            if a.len() != b.len() || a.iter().zip(&b).any(|(x, y)| (x - y).abs() > tol) { return None; }
                            // Final ID of that group once variants and swaps are applied.
                            let final_id = if !is_dom {
                                let v = *combo_to_new_tx.get(&(*q, key.0, key.1))?;
                                if swapped.get(q) == Some(&v) { *q } else { v }
                            } else {
                                swapped.get(q).copied().unwrap_or(*q)
                            };
                            Some((total, final_id))
                        })
                        .max()
                        .map(|(_, id)| id)
                } else {
                    None
                };
                match target {
                    Some(id) => demote_merge.push((parent, id)),
                    None => demote_new.push((parent, self.refs.counters.transcript.increment())),
                }
            }
            eprintln!(
                "[IsoForge] Annotated IDs contradicted by their reads: {} transcript(s); {} joined a group of the same \
                 gene on the same TSS and TES clusters with a splice chain within {} bp, {} became end variants",
                id_demotes.len(), demote_merge.len(), tol, demote_new.len()
            );
        }
        variant_parents.extend(demote_new.iter().map(|&(parent, new_id)| (new_id, parent, false)));
        results.end_variant_parents.extend(variant_parents);

        let n_new_variants = variants.len();
        if n_new_variants == 0 && demote_new.is_empty() && demote_merge.is_empty() {
            return 0;
        }

        eprintln!(
            "[IsoForge] Creating {} TSS/TES transcript variants...",
            n_new_variants
        );

        // ── Create transcript models for new variants ────────────────────
        // Parent models come from the annotation dict, then from novel transcripts.
        let parent_ids: HashSet<TranscriptId> =
            variants.iter().map(|v| v.parent_tx_id).collect();

        let mut parent_models: HashMap<TranscriptId, &TranscriptModel> = HashMap::default();
        for m in self.refs.transcript_dict.values() {
            if parent_ids.contains(&m.transcript_id) {
                parent_models.insert(m.transcript_id, m);
            }
        }
        // Novel transcripts
        for m in results.novel_transcripts.iter() {
            if parent_ids.contains(&m.transcript_id) {
                parent_models.insert(m.transcript_id, m);
            }
        }

        // The variants' models: copies of their parents under the new IDs.
        let mut new_models: Vec<TranscriptModel> = Vec::with_capacity(n_new_variants);
        for v in &variants {
            if let Some(&orig) = parent_models.get(&v.parent_tx_id) {
                let mut variant_model = orig.clone();
                variant_model.transcript_id = v.new_tx_id;
                new_models.push(variant_model);
            }
        }
        drop(parent_models);
        let demote_parents: HashSet<TranscriptId> = demote_new.iter().map(|d| d.0).collect();
        let demote_models: HashMap<TranscriptId, &TranscriptModel> = if demote_parents.is_empty() {
            HashMap::default()
        } else {
            self.refs.transcript_dict.values().filter(|m| demote_parents.contains(&m.transcript_id))
                .map(|m| (m.transcript_id, m)).collect()
        };
        for &(parent, new_id) in &demote_new {
            if let Some(orig) = demote_models.get(&parent) {
                let mut model = (*orig).clone();
                model.transcript_id = new_id;
                new_models.push(model);
            }
        }
        results.novel_transcripts.extend(new_models);

        // ── Update abundance counts ──────────────────────────────────────
        // Index tx_id → [(gid, ds)] built once.
        let mut tx_abundance_keys: HashMap<TranscriptId, Vec<(GeneId, Arc<str>)>> =
            HashMap::default();
        for (gid, tid, ds) in results.abundance.keys() {
            if parent_ids.contains(tid) {
                tx_abundance_keys
                    .entry(*tid)
                    .or_default()
                    .push((*gid, ds.clone()));
            }
        }

        for v in &variants {
            if let Some(keys) = tx_abundance_keys.get(&v.parent_tx_id) {
                for (gid, ds) in keys {
                    if let Some(&variant_count) = v.ds_counts.get(ds) {
                        if variant_count > 0 {
                            // Subtract from original
                            let orig_key = (*gid, v.parent_tx_id, ds.clone());
                            if let Some(orig_count) = results.abundance.get_mut(&orig_key) {
                                *orig_count -= variant_count;
                            }
                            // Add to variant
                            let new_key = (*gid, v.new_tx_id, ds.clone());
                            *results.abundance.entry(new_key).or_insert(0) += variant_count;
                        }
                    }
                }
            }

            log::debug!(
                "[IsoForge] Created TSS/TES variant {p}T{:0width$} from {p}T{:0width$} ({} reads)",
                v.new_tx_id,
                v.parent_tx_id,
                v.n_reads,
                p = self.config.id_prefix,
                width = self.config.n_places,
            );
        }

        log::info!(
            "[IsoForge] Created {} TSS/TES cluster variants",
            variants.len()
        );

        // ── Pass 2: relabel reads whose (tx_id, tss_cid, tes_cid) has a new ID ──
        for (idx, annot) in cores.iter_mut().enumerate() {
            if !keep.get(idx).copied().unwrap_or(false) || upstream_unclustered.contains(&idx) {
                continue;
            }
            let info = join_results.get(idx);
            let tss_cid = info.and_then(|e| e.tss_cluster_id);
            let tes_cid = info.and_then(|e| e.tes_cluster_id);

            if let Some(&new_tx_id) =
                combo_to_new_tx.get(&(annot.transcript_id, tss_cid, tes_cid))
            {
                annot.transcript_id = new_tx_id;
            }
        }

        // ── Exchange IDs where the annotated ends sit on a minority group ──
        if !id_swaps.is_empty() {
            let swap_map: HashMap<TranscriptId, TranscriptId> =
                id_swaps.iter().flat_map(|&(p, v)| [(p, v), (v, p)]).collect();
            for (idx, annot) in cores.iter_mut().enumerate() {
                if !keep.get(idx).copied().unwrap_or(false) || upstream_unclustered.contains(&idx) {
                    continue;
                }
                if let Some(&t) = swap_map.get(&annot.transcript_id) {
                    annot.transcript_id = t;
                }
            }
            let keys: Vec<_> = results.abundance.keys()
                .filter(|(_, t, _)| swap_map.contains_key(t)).cloned().collect();
            let mut moved = Vec::with_capacity(keys.len());
            for key in keys {
                if let Some(n) = results.abundance.remove(&key) {
                    moved.push(((key.0, swap_map[&key.1], key.2), n));
                }
            }
            for (key, n) in moved {
                *results.abundance.entry(key).or_insert(0) += n;
            }
            eprintln!(
                "[IsoForge] Annotated IDs follow annotated ends: {} transcript(s) now keep their ID on \
                 the reads matching the annotated TSS/TES; the larger group starting elsewhere is the variant",
                id_swaps.len()
            );
        }

        // ── Relabel reads of annotated IDs their reads contradict ─────────
        // Reads that reach the transcript's distinguishing terminal sequence stay
        // on the annotated ID; the rest move.
        if !demote_merge.is_empty() || !demote_new.is_empty() {
            let relabel: HashMap<TranscriptId, TranscriptId> =
                demote_merge.iter().chain(demote_new.iter()).copied().collect();
            let reaches = |annot: &ObsCore| -> bool {
                let (Some(&(a_tss, a_tes, _)), Some(&(b5, b3))) = (annotated_ends.get(&annot.transcript_id), distinguishing.get(&annot.transcript_id)) else { return false };
                let minus = a_tss > a_tes;
                let (lo, hi) = (annot.start_pos.min(annot.end_pos), annot.start_pos.max(annot.end_pos));
                let (five, three) = if minus { (hi, lo) } else { (lo, hi) };
                b5.is_some_and(|b| (if minus { five - b } else { b - five }) > join_distance)
                    || b3.is_some_and(|b| (if minus { b - three } else { three - b }) > join_distance)
            };
            let mut stayed: HashMap<(GeneId, TranscriptId, Arc<str>), i64> = HashMap::default();
            for (idx, annot) in cores.iter_mut().enumerate() {
                if !keep.get(idx).copied().unwrap_or(false) || upstream_unclustered.contains(&idx) {
                    continue;
                }
                if let Some(&t) = relabel.get(&annot.transcript_id) {
                    if reaches(annot) {
                        *stayed.entry((annot.gene_id, annot.transcript_id, Arc::clone(&annot.dataset))).or_insert(0) += 1;
                    } else {
                        annot.transcript_id = t;
                    }
                }
            }
            let keys: Vec<_> = results.abundance.keys()
                .filter(|(_, t, _)| relabel.contains_key(t)).cloned().collect();
            let mut moved = Vec::with_capacity(keys.len());
            for key in keys {
                let keep_here = stayed.get(&key).copied().unwrap_or(0);
                if let Some(n) = results.abundance.get_mut(&key) {
                    let go = (*n - keep_here).max(0);
                    *n -= go;
                    moved.push(((key.0, relabel[&key.1], key.2.clone()), go));
                    if *n <= 0 { results.abundance.remove(&key); }
                }
            }
            for (key, n) in moved {
                if n > 0 { *results.abundance.entry(key).or_insert(0) += n; }
            }
            let n_stayed: i64 = stayed.values().sum();
            if n_stayed > 0 {
                eprintln!(
                    "[IsoForge] Refused annotated IDs: {} reads that reach the transcript's distinguishing terminal sequence stay on it",
                    n_stayed
                );
            }
        }

        if n_absorbed_untailed_start > 0 {
            eprintln!(
                "[IsoForge] End variants: {} reads of alternative-start groups without {} tailed reads absorbed into the dominant model",
                n_absorbed_untailed_start, MIN_FULL_LENGTH_START_READS
            );
        }
        n_new_variants + demote_new.len()
    }

    /// Parse a BAM and pass its reads to `on_window` in windows of at most
    /// `window_size`. Without a genome reference the whole file is one window.
    fn read_bam_windows<F>(
        &self,
        bam_path: &Path,
        dataset_name: &str,
        window_size: usize,
        junctions: Option<Arc<crate::junction_groups::JunctionGroups>>,
        mut on_window: F,
    ) -> IsoForgeResult<()>
    where
        F: FnMut(Vec<ParsedRead>) -> IsoForgeResult<()>,
    {
        if let Some(ref genome_path) = self.genome_path {
            let strand_index = if self.config.infer_strand_from_polya {
                let cached = self.strand_indexes.lock().unwrap().get(bam_path).cloned();
                match cached {
                    Some(ix) => { ix.reset_counts(); Some(ix) }
                    None => {
                        let ix = Arc::new(crate::strand_evidence::StrandIndex::build(bam_path, &self.config, &self.refs)?);
                        ix.log_survey(dataset_name);
                        self.strand_indexes.lock().unwrap().insert(bam_path.to_path_buf(), ix.clone());
                        Some(ix)
                    }
                }
            } else {
                None
            };
            let pipeline = ReaderPipeline {
                strand_index: strand_index.clone(),
                microexons: self.config.repair_microexons.then(|| Arc::new(crate::microexon_repair::MicroexonIndex::from_reference(&self.refs))),
                splice_realign: (self.config.splice_realign_window > 0).then(|| Arc::new(crate::splice_realignment::SpliceRealignIndex::from_reference(&self.refs, self.config.splice_realign_window))),
                terminal_rescue: self.config.terminal_exon_rescue.then(|| Arc::new(crate::terminal_rescue::TerminalExonIndex::from_reference(&self.refs))),
                junction_groups: junctions.clone(),
                small_rna: self.refs.small_rna.clone(),
                genome_path: Some(genome_path.clone()),
                sj_db: self.sj_db.clone(),
                correct_sj: self.config.correct_sj,
                max_sj_offset: self.config.sj_correction_window.unwrap_or(self.config.max_splice_dist),
                frac_a_range: 20,
                compute_md_if_missing: true,
                correct_microindels: self.config.correct_microindels,
                max_len_indel: self.config.max_indel_size,
            };

            let mut reader = IntegratedBamReader::new(bam_path.to_str().unwrap_or(""), pipeline)?;
            let prepass = junctions.as_ref().is_some_and(|j| j.is_prepass());
            reader.collect_chimeric = self.config.detect_fusions && !prepass;
            if self.config.detect_fusions && self.config.fusion_spliced_joins {
                reader.gene_spans = Some(Arc::new(crate::fusion::GeneSpans::from_tree(&self.refs.gene_tree)));
            }
            reader.process_windows(&self.config, dataset_name, window_size, on_window)?;
            if let Some(ix) = strand_index.filter(|_| !prepass) {
                ix.log_resolution(dataset_name);
            }
            if self.config.detect_fusions && !prepass {
                let segs = reader.take_chimeric();
                self.chimeric.lock().unwrap().extend(segs);
            }
            Ok(())
        } else {
            let reader = BamReaderWithHeader::new(bam_path.to_str().unwrap_or(""))?;
            on_window(reader.read_all_with_chroms(&self.config)?)
        }
    }

    /// Run TSS/TES clustering on all reads from all datasets.
    ///
    /// Steps:
    ///   1. Load UTR regions from annotation directory
    ///   2. Accumulate end-site counts from all reads (with novelty-aware filtering)
    ///   3. Call peaks and form clusters
    ///   4. Annotate clusters with genomic features (UTR, gene body)
    ///   5. Scan for PAS hexamers near TES clusters
    ///   6. Join reads to nearest clusters
    ///   7. Revalidate clusters against the joined reads
    fn run_clustering(
        &self,
        cores: &[ObsCore],
        dataset_ranges: &[(String, std::ops::Range<usize>)],
        n_reads: usize,
        cluster_cfg: &ClusterConfig,
        hints: Option<&tss_tes_clustering::ReadFilterHints>,
        cage_index: Option<&CagePeakIndex>,
    ) -> IsoForgeResult<(Vec<ReadEndSiteInfo>, Vec<EndSiteCluster>, Vec<EndSiteCluster>)> {
        eprintln!("[IsoForge] ── TSS/TES Clustering ──");

        // Load UTR regions
        let mut all_feature_regions = if let Some(ref annot_dir) = self.annot_dir {
            let utr_path = annot_dir.join("utr_regions.parquet");
            if utr_path.exists() {
                tss_tes_clustering::load_utr_regions(&utr_path)?
            } else {
                eprintln!("[IsoForge] Warning: utr_regions.parquet not found, feature annotation disabled");
                Vec::new()
            }
        } else {
            eprintln!("[IsoForge] Warning: annotation directory not set, feature annotation disabled");
            Vec::new()
        };

        // Load exon/intron regions from the reference annotation
        if let Some(ref annot_dir) = self.annot_dir {
            let exon_intron_path = annot_dir.join("exon_intron_regions.parquet");
            if exon_intron_path.exists() {
                let exon_intron = tss_tes_clustering::load_utr_regions(&exon_intron_path)?;
                eprintln!("[IsoForge] Loaded {} exon/intron regions for feature annotation", exon_intron.len());
                all_feature_regions.extend(exon_intron);
            }
        }

        // Load gene records for gene-body feature annotation
        let genes = if let Some(ref annot_dir) = self.annot_dir {
            let genes_path = annot_dir.join("genes.parquet");
            if genes_path.exists() {
                crate::parquet_io::load_genes(&genes_path)?
            } else {
                Vec::new()
            }
        } else {
            Vec::new()
        };

        // Open genome for PAS scanning
        let genome = if let Some(ref genome_path) = self.genome_path {
            Some(GenomeRef::open(genome_path)?)
        } else {
            None
        };

        // Each annotated transcript start and end contributes a UTR-labelled window of
        // ANNOTATED_END_FEATURE_WINDOW bp either side.
        {
            let gene_strand: HashMap<GeneId, Strand> = genes.iter().map(|g| (g.gene_id, g.strand)).collect();
            let w = ANNOTATED_END_FEATURE_WINDOW;
            for m in self.refs.transcript_dict.values() {
                let Some(&strand) = gene_strand.get(&m.gene_id) else { continue };
                let (tss, tes) = if strand == Strand::Plus { (m.start_pos, m.end_pos) } else { (m.end_pos, m.start_pos) };
                for (pos, feature) in [(tss, GenomicFeature::FivePrimeUtr), (tes, GenomicFeature::ThreePrimeUtr)] {
                    all_feature_regions.push(UtrRegion {
                        chromosome: m.chromosome.clone(),
                        start: pos - w,
                        end: pos + w,
                        strand,
                        feature,
                        gene_id: m.gene_id.to_string(),
                        transcript_id: m.transcript_id.to_string(),
                    });
                }
            }
        }

        // Build feature annotation tree from UTR + exon/intron regions + gene bodies
        let feature_tree = tss_tes_clustering::FeatureAnnotationTree::build(&all_feature_regions, &genes);

        // Create engine and accumulate counts
        let mut engine = TssTesEngine::new(cluster_cfg.clone());

        eprintln!("[IsoForge] Accumulating end-site counts from {} reads across {} datasets",
                 cores.len(), dataset_ranges.len());

        // Datasets accumulate into one set of counts in turn.
        for (_dataset_name, range) in dataset_ranges.iter() {
            engine.counts.add_reads(&cores[range.clone()], cluster_cfg, hints, cage_index);
        }

        // Run clustering
        eprintln!("[IsoForge] Calling TSS/TES peaks...");
        engine.run_clustering(&feature_tree, genome.as_ref());

        if cluster_cfg.reference_end_support {
            crate::reference_end_support::supplement(&mut engine, &self.refs, cores, hints,
                cage_index, &feature_tree, genome.as_ref(), &self.output_dir.join("reference_end_support.tsv"))?;
        }

        // Join reads to clusters
        eprintln!("[IsoForge] Joining reads to clusters...");
        let mut join_results = engine.join_reads(cores, n_reads, hints, cage_index);

        // ── Post-filter cluster revalidation ─────────────────────────────
        // Clusters that lost too many supported reads have their remaining reads
        // revoked, iterated to convergence; Known/FSM reads are exempt.
        engine.revalidate_clusters(&mut join_results, hints);

        // Clusters are written in run() after transcript-level filtering.
        let tss_clusters = engine.tss_clusters;
        let tes_clusters = engine.tes_clusters;

        eprintln!("[IsoForge] Returning {} TSS + {} TES clusters for post-filter writing",
                 tss_clusters.len(), tes_clusters.len());

        Ok((join_results, tss_clusters, tes_clusters))
    }
}

/// Per-batch intermediate result.
struct BatchResult {
    annotations: Vec<ReadAnnotation>,
    novel_buffer: NovelBuffer,
    qc_log: Vec<QcMetrics>,
}

/// Build positions array from start, end, and intron coordinates.
/// Positions format: [exon1_start, exon1_end, exon2_start, exon2_end, ...]
#[allow(dead_code)]
fn build_positions_from_introns(
    read_start: i64,
    read_end: i64,
    introns: &[(i64, i64)],
) -> Vec<i64> {
    let mut positions = Vec::with_capacity(2 + introns.len() * 2);

    positions.push(read_start);

    for (_i, (intron_start, intron_end)) in introns.iter().enumerate() {
        positions.push(intron_start - 1);   // exon end
        positions.push(intron_end + 1);     // next exon start
    }

    positions.push(read_end);
    positions
}

/// Merge one IsoForgeResults into another.
fn merge_into(target: &mut IsoForgeResults, source: IsoForgeResults) {
    target.observed.extend(source.observed);
    target.gene_annotations.extend(source.gene_annotations);
    target.transcript_annotations.extend(source.transcript_annotations);
    target.exon_annotations.extend(source.exon_annotations);
    target.novel_genes.extend(source.novel_genes);
    target.novel_transcripts.extend(source.novel_transcripts);
    target.novel_edges.extend(source.novel_edges);
    target.novel_locations.extend(source.novel_locations);
    target.qc_log.extend(source.qc_log);
    target.total_reads += source.total_reads;
    target.reads_processed += source.reads_processed;
    target.reads_failed_qc += source.reads_failed_qc;

    for (key, count) in source.abundance {
        *target.abundance.entry(key).or_insert(0) += count;
    }
}

/// Per-transcript TSS/TES peak tallies used to anchor novel transcript ends,
/// collected after splitting and filtering. Cluster peaks take precedence;
/// read boundaries are used where no peak was called.
#[derive(Default)]
pub struct AnchorVotes {
    tss: HashMap<TranscriptId, HashMap<i64, usize>>,
    tes: HashMap<TranscriptId, HashMap<i64, usize>>,
    /// Read alignment boundaries, used where a cluster peak is unavailable.
    tss_read: HashMap<TranscriptId, HashMap<i64, usize>>,
    tes_read: HashMap<TranscriptId, HashMap<i64, usize>>,
}

impl AnchorVotes {
    fn collect(&mut self, cores: &[ObsCore], keep: &[bool], join_results: &[ReadEndSiteInfo]) {
        for (idx, annot) in cores.iter().enumerate() {
            if !keep.get(idx).copied().unwrap_or(false) {
                continue;
            }
            let Some(info) = join_results.get(idx) else { continue };
            if !info.is_supported {
                continue;
            }
            if let Some(p) = info.tss_site {
                *self.tss.entry(annot.transcript_id).or_default().entry(p).or_insert(0) += 1;
            }
            if let Some(p) = info.polya_site {
                *self.tes.entry(annot.transcript_id).or_default().entry(p).or_insert(0) += 1;
            }
            *self.tss_read.entry(annot.transcript_id).or_default()
                .entry(annot.start_pos).or_insert(0) += 1;
            *self.tes_read.entry(annot.transcript_id).or_default()
                .entry(annot.end_pos).or_insert(0) += 1;
        }
    }

    /// Most-voted position; ties resolve to the smaller coordinate.
    fn modal(m: Option<&HashMap<i64, usize>>) -> Option<i64> {
        m.and_then(|m| {
            m.iter()
                .max_by(|a, b| a.1.cmp(b.1).then_with(|| b.0.cmp(a.0)))
                .map(|(pos, _)| *pos)
        })
    }

    fn modal_tss(&self, tx: TranscriptId) -> Option<i64> {
        Self::modal(self.tss.get(&tx)).or_else(|| Self::modal(self.tss_read.get(&tx)))
    }

    fn modal_tes(&self, tx: TranscriptId) -> Option<i64> {
        Self::modal(self.tes.get(&tx)).or_else(|| Self::modal(self.tes_read.get(&tx)))
    }
}

/// Read count from the BAM index, or 0 without an index; used to size buffers.
fn count_indexed_reads(bam_path: &Path) -> usize {
    let Ok(mut reader) = rust_htslib::bam::IndexedReader::from_path(bam_path) else {
        return 0;
    };
    match reader.index_stats() {
        Ok(stats) => stats.iter().map(|(_, _, mapped, unmapped)| mapped + unmapped).sum::<u64>()
            as usize,
        Err(_) => 0,
    }
}

/// Which entity an annotation row describes. Gene, transcript and edge IDs are
/// separate counters; the row's attribute tells them apart.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum IdSpace {
    Gene,
    Transcript,
    Edge,
}

fn annotation_id_space(attribute: &str) -> Option<IdSpace> {
    if attribute == "exon_status" {
        return Some(IdSpace::Edge);
    }
    if attribute.starts_with("gene_") || attribute == "intergenic_novel" || attribute == "antisense_gene" {
        return Some(IdSpace::Gene);
    }
    if attribute.starts_with("transcript_")
        || attribute.starts_with("FSM")
        || attribute.starts_with("ISM")
        || attribute.starts_with("NIC")
        || attribute.starts_with("NNC")
        || attribute.ends_with("_transcript")
        || attribute == "TSS_TES_variant"
    {
        return Some(IdSpace::Transcript);
    }
    None
}

/// The canonical renumbering of novel entity IDs (`renumber_novel_entities`): old ID -> final ID per ID space. An ID
/// that is not in a map is a reference ID and stays.
#[derive(Default)]
pub struct NovelRemap {
    vertices: rustc_hash::FxHashMap<VertexId, VertexId>,
    edges: rustc_hash::FxHashMap<EdgeId, EdgeId>,
    genes: rustc_hash::FxHashMap<GeneId, GeneId>,
    transcripts: rustc_hash::FxHashMap<TranscriptId, TranscriptId>,
}

impl NovelRemap {
    #[inline]
    fn vertex(&self, id: VertexId) -> VertexId {
        self.vertices.get(&id).copied().unwrap_or(id)
    }
    #[inline]
    fn edge(&self, id: EdgeId) -> EdgeId {
        self.edges.get(&id).copied().unwrap_or(id)
    }
    #[inline]
    fn gene(&self, id: GeneId) -> GeneId {
        self.genes.get(&id).copied().unwrap_or(id)
    }
    #[inline]
    fn transcript(&self, id: TranscriptId) -> TranscriptId {
        self.transcripts.get(&id).copied().unwrap_or(id)
    }

    fn jn_path(&self, path: &Option<String>) -> Option<String> {
        path.as_ref().map(|p| {
            p.split(',')
                .map(|tok| match tok.trim().parse::<EdgeId>() {
                    Ok(id) => self.edge(id).to_string(),
                    Err(_) => tok.to_string(),
                })
                .collect::<Vec<_>>()
                .join(",")
        })
    }

    /// Rewrite the `id` column of the streamed annotations file after
    /// renumbering; rows of `dropped` transcripts are removed.
    fn rewrite_annotation_ids(
        &self,
        path: &Path,
        dropped: &rustc_hash::FxHashSet<TranscriptId>,
        tx_merges: &HashMap<TranscriptId, TranscriptId>,
        ref_max: &RefIdBounds,
        id_prefix: &str,
        n_places: usize,
    ) -> IsoForgeResult<()> {
        if !path.exists() {
            return Ok(());
        }
        // Rows whose novel ID is not in the remap are dropped.
        parquet_io::remap_annotation_ids_by_attribute(
            path,
            &|attribute, id| match annotation_id_space(attribute) {
                Some(IdSpace::Gene) => match self.genes.get(&id) {
                    Some(&g) => Some(g),
                    None if id > ref_max.gene => None,
                    None => Some(id),
                },
                Some(IdSpace::Transcript) => {
                    // Merged transcripts are resolved to their survivor before renumbering.
                    let id = tx_merges.get(&id).copied().unwrap_or(id);
                    match self.transcripts.get(&id) {
                        Some(&t) => Some(t),
                        None if dropped.contains(&id) || id > ref_max.transcript => None,
                        None => Some(id),
                    }
                }
                Some(IdSpace::Edge) => match self.edges.get(&id) {
                    Some(&e) => Some(e),
                    None if id > ref_max.edge => None,
                    None => Some(id),
                },
                // Unrecognised attribute: leave the row untouched.
                None => Some(id),
            },
            // Generated names embed the numeric ID and are rebuilt; reference names are unchanged.
            &|attribute, new_id, old_value| {
                let marker = match attribute {
                    "gene_id" | "gene_name" => 'G',
                    "transcript_id" | "transcript_name" => 'T',
                    _ => return None,
                };
                let generated = format!("{}{}", id_prefix, marker);
                if !old_value.starts_with(&generated) {
                    return None;
                }
                Some(format!("{}{:0width$}", generated, new_id, width = n_places))
            },
        )
    }
}

/// The strand of an annotated transcript read off the genes at its start vertex, for an index without a gene table:
/// its own gene's entry, else the entry of the lowest gene ID.
fn strand_at_start(vertex_2_gene: &Vertex2Gene, start_vertex: VertexId, gene_id: GeneId) -> Strand {
    vertex_2_gene
        .get(&start_vertex)
        .and_then(|set| set.iter().find(|(g, _)| *g == gene_id).or_else(|| set.iter().min_by_key(|(g, _)| *g)))
        .map(|(_, s)| *s)
        .unwrap_or(Strand::Plus)
}

/// Number `items` from `base + 1` in the order of their keys; items with equal keys keep the order they come in,
/// and an ID listed twice is numbered once.
fn assign_dense_ids<K: Ord>(items: &[(i64, K)], base: i64) -> rustc_hash::FxHashMap<i64, i64> {
    let mut order: Vec<&(i64, K)> = items.iter().collect();
    order.sort_by(|a, b| a.1.cmp(&b.1));

    let mut out = rustc_hash::FxHashMap::default();
    for (id, _) in order {
        if !out.contains_key(id) {
            let next = base + 1 + out.len() as i64;
            out.insert(*id, next);
        }
    }
    out
}

/// Prints the heap bytes held by each container that survives a window.
fn debug_heap_report(
    cores: &[ObsCore],
    results: &IsoForgeResults,
    seen_annotations: usize,
    refs: &SharedRefs,
) {
    const MB: f64 = 1024.0 * 1024.0;
    let str_bytes = |s: &String| std::mem::size_of::<String>() + s.capacity();

    let cores_b = cores.len() * std::mem::size_of::<ObsCore>();
    let tx_b: usize = results
        .novel_transcripts
        .iter()
        .map(|t| {
            std::mem::size_of::<TranscriptModel>()
                + t.chromosome.capacity()
                + t.jn_path.as_ref().map_or(0, |s| s.capacity())
        })
        .sum();
    let edges_b = results.novel_edges.len() * std::mem::size_of::<Edge>();
    let locs_b: usize = results
        .novel_locations
        .iter()
        .map(|l| {
            std::mem::size_of::<Location>()
                + l.chromosome.capacity()
                + l.genome_build.capacity()
        })
        .sum();
    let abund_b = results.abundance.len()
        * (std::mem::size_of::<(GeneId, TranscriptId, Arc<str>)>() + std::mem::size_of::<i64>());
    let annot_b = seen_annotations * std::mem::size_of::<(i64, Arc<str>, Arc<str>)>();
    let vtx_dict_b = refs.novel_location_dict.len()
        * (std::mem::size_of::<(String, i64)>() + std::mem::size_of::<VertexId>() + 24);
    let edge_dict_b = refs.novel_edge_dict.len()
        * (std::mem::size_of::<EdgeKey>() + std::mem::size_of::<NovelEdgeEntry>());
    let gene_dict_b: usize = refs
        .novel_gene_dict
        .iter()
        .map(|e| {
            std::mem::size_of::<NovelGeneKey>()
                + std::mem::size_of::<NovelGeneEntry>()
                + str_bytes(&e.chromosome)
        })
        .sum();

    let tracked = cores_b + tx_b + edges_b + locs_b + abund_b + annot_b
        + vtx_dict_b + edge_dict_b + gene_dict_b;

    eprintln!(
        "[IsoForge][heap] cores={:.0}M novel_tx={:.0}M novel_edges={:.0}M \
         novel_locs={:.0}M abundance={:.0}M annots={:.0}M vtx_dict={:.0}M \
         edge_dict={:.0}M gene_dict={:.0}M | tracked={:.1}G{}",
        cores_b as f64 / MB,
        tx_b as f64 / MB,
        edges_b as f64 / MB,
        locs_b as f64 / MB,
        abund_b as f64 / MB,
        annot_b as f64 / MB,
        vtx_dict_b as f64 / MB,
        edge_dict_b as f64 / MB,
        gene_dict_b as f64 / MB,
        tracked as f64 / (MB * 1024.0),
        jemalloc_summary(),
    );
}

#[cfg(target_os = "linux")]
fn jemalloc_summary() -> String {
    use tikv_jemalloc_ctl::{epoch, stats};
    const G: f64 = 1024.0 * 1024.0 * 1024.0;
    if epoch::advance().is_err() {
        return String::new();
    }
    match (stats::allocated::read(), stats::resident::read()) {
        (Ok(a), Ok(r)) => format!(" live={:.1}G rss={:.1}G", a as f64 / G, r as f64 / G),
        _ => String::new(),
    }
}

#[cfg(not(target_os = "linux"))]
fn jemalloc_summary() -> String {
    String::new()
}

/// Models whose terminal exon edges do not end at their terminal vertices (debug; ISOFORGE_END_CHECK).
fn log_end_consistency(stage: &str, results: &IsoForgeResults, refs: &SharedRefs) {
    if std::env::var_os("ISOFORGE_END_CHECK").is_none() { return; }
    let mut ends: rustc_hash::FxHashMap<EdgeId, (VertexId, VertexId)> = refs.edge_dict.values().map(|e| (e.edge_id, (e.v1, e.v2))).collect();
    ends.extend(results.novel_edges.iter().map(|e| (e.edge_id, (e.v1, e.v2))));
    let (mut bad5, mut bad3, mut missing) = (0usize, 0usize, 0usize);
    for t in &results.novel_transcripts {
        match (ends.get(&t.start_exon), ends.get(&t.end_exon)) {
            (Some(s), Some(e)) => {
                if s.0 != t.start_vertex { bad5 += 1; }
                if e.1 != t.end_vertex { bad3 += 1; }
            }
            _ => missing += 1,
        }
    }
    eprintln!("[IsoForge] End check {stage}: {} novel models, start edge off its vertex {bad5}, end edge off its vertex {bad3}, edges missing {missing}",
              results.novel_transcripts.len());
}

/// Canonical renumbering of novel entity IDs. Every ID space is numbered densely above the reference's highest ID, in
/// an order of the entities themselves (vertices by position, edges by their vertices, genes by locus, transcripts by
/// structure). Dropped novel transcripts (`dropped`) are numbered after the kept ones. The tables come out in ID order.
fn renumber_novel_entities(results: &mut IsoForgeResults, ref_max: &RefIdBounds, dropped: &[TranscriptModel],
                           first_read: &rustc_hash::FxHashMap<TranscriptId, ReadIdx>) -> NovelRemap {
    let mut remap = NovelRemap::default();

    // Vertices, keyed by genomic position.
    let loc_keys: Vec<_> = results
        .novel_locations
        .iter()
        .map(|l| (l.location_id, (l.chromosome.clone(), l.position)))
        .collect();
    remap.vertices = assign_dense_ids(&loc_keys, ref_max.vertex);
    for l in results.novel_locations.iter_mut() {
        l.location_id = remap.vertex(l.location_id);
    }

    // Edges, keyed by their (already renumbered) endpoints.
    let edge_keys: Vec<_> = results
        .novel_edges
        .iter()
        .map(|e| {
            (
                e.edge_id,
                (remap.vertex(e.v1), remap.vertex(e.v2), e.edge_type as u8, e.strand),
            )
        })
        .collect();
    remap.edges = assign_dense_ids(&edge_keys, ref_max.edge);
    for e in results.novel_edges.iter_mut() {
        e.edge_id = remap.edge(e.edge_id);
        e.v1 = remap.vertex(e.v1);
        e.v2 = remap.vertex(e.v2);
    }

    // Genes, keyed by locus; genes sharing a locus keep the order they come in (by registry key).
    let gene_keys: Vec<_> = results
        .novel_genes
        .iter()
        .map(|g| (g.gene_id, (g.chromosome.clone(), g.strand, g.start, g.end)))
        .collect();
    remap.genes = assign_dense_ids(&gene_keys, ref_max.gene);
    for g in results.novel_genes.iter_mut() {
        g.gene_id = remap.gene(g.gene_id);
    }

    // Transcripts, keyed by locus then exon chain: the kept ones, then the dropped ones.
    let tx_key = |t: &TranscriptModel| (
        t.transcript_id,
        (
            t.chromosome.clone(),
            t.start_pos,
            t.end_pos,
            t.n_exons,
            remap.edge(t.start_exon),
            remap.jn_path(&t.jn_path),
            remap.edge(t.end_exon),
            // then the model's first read, for end variants that share all of the above
            first_read.get(&t.transcript_id).copied().unwrap_or(ReadIdx::MAX),
        ),
    );
    let tx_keys: Vec<_> = results.novel_transcripts.iter().map(tx_key).collect();
    let dropped_keys: Vec<_> = dropped.iter().map(tx_key).collect();
    if std::env::var_os("ISOFORGE_END_CHECK").is_some() {
        let mut k: Vec<_> = tx_keys.iter().chain(dropped_keys.iter()).map(|x| &x.1).collect();
        k.sort();
        let dup = k.windows(2).filter(|w| w[0] == w[1]).count();
        eprintln!("[IsoForge] Renumbering: {} novel transcripts kept, {} dropped, {dup} with the key of another", tx_keys.len(), dropped_keys.len());
    }
    let mut transcripts = assign_dense_ids(&tx_keys, ref_max.transcript);
    let kept = transcripts.len() as i64;
    transcripts.extend(assign_dense_ids(&dropped_keys, ref_max.transcript + kept));
    remap.transcripts = transcripts;
    for t in results.novel_transcripts.iter_mut() {
        t.transcript_id = remap.transcript(t.transcript_id);
        t.gene_id = remap.gene(t.gene_id);
        t.start_exon = remap.edge(t.start_exon);
        t.end_exon = remap.edge(t.end_exon);
        t.jn_path = remap.jn_path(&t.jn_path);
        t.start_vertex = remap.vertex(t.start_vertex);
        t.end_vertex = remap.vertex(t.end_vertex);
    }

    results.abundance = results
        .abundance
        .drain()
        .map(|((g, t, ds), c)| ((remap.gene(g), remap.transcript(t), ds), c))
        .collect();

    results.novel_locations.sort_by_key(|l| l.location_id);
    results.novel_edges.sort_by_key(|e| e.edge_id);
    results.novel_genes.sort_by_key(|g| g.gene_id);
    results.novel_transcripts.sort_by_key(|t| t.transcript_id);

    remap
}

/// Rescue coherent sub-clusters among outlier reads.
///
/// Outlier positions are sorted and grouped while consecutive reads lie
/// within `max_gap` bp; groups with at least `min_reads` reads are not
/// flagged as degraded.
///
/// Returns (rescued read indices, number of sub-clusters, number of rescued reads).
fn rescue_coherent_subclusters(
    outliers: &mut [(ReadIdx, i64)],
    max_gap: i64,
    min_reads: usize,
) -> (rustc_hash::FxHashSet<ReadIdx>, u64, u64) {
    let mut rescued = rustc_hash::FxHashSet::default();

    if outliers.is_empty() {
        return (rescued, 0, 0);
    }

    // Sort by genomic position
    outliers.sort_by_key(|(_, pos)| *pos);

    // Greedy grouping: a new cluster starts when the gap exceeds max_gap
    let mut n_subclusters = 0u64;
    let mut n_rescued = 0u64;
    let mut cluster_start = 0usize;

    for i in 1..=outliers.len() {
        let gap_exceeded = i == outliers.len()
            || (outliers[i].1 - outliers[i - 1].1).abs() > max_gap;

        if gap_exceeded {
            let cluster = &outliers[cluster_start..i];
            if cluster.len() >= min_reads {
                n_subclusters += 1;
                for (read_idx, _) in cluster {
                    rescued.insert(*read_idx);
                    n_rescued += 1;
                }
            }
            cluster_start = i;
        }
    }

    (rescued, n_subclusters, n_rescued)
}
#[cfg(test)]
mod shared_end_tests {
    use super::shared_end_variant;

    #[test]
    fn a_group_missing_one_cluster_joins_the_variant_sharing_the_other() {
        let dominant = (Some(735), Some(195));
        let made = [(Some(735), Some(514)), (Some(54), Some(195)), (Some(54), Some(514))];
        // No start cluster, the poly(A) cluster of a variant: that variant, the largest first.
        assert_eq!(shared_end_variant((None, Some(514)), dominant, made.iter().copied()), Some(0));
        // No poly(A) cluster, the start cluster of a variant.
        assert_eq!(shared_end_variant((Some(54), None), dominant, made.iter().copied()), Some(1));
        // The dominant group's own cluster, both clusters, or neither: not a variant's.
        assert_eq!(shared_end_variant((None, Some(195)), dominant, made.iter().copied()), None);
        assert_eq!(shared_end_variant((Some(735), None), dominant, made.iter().copied()), None);
        assert_eq!(shared_end_variant((Some(9), Some(514)), dominant, made.iter().copied()), None);
        assert_eq!(shared_end_variant((None, None), dominant, made.iter().copied()), None);
        assert_eq!(shared_end_variant((None, Some(777)), dominant, made.iter().copied()), None);
    }
}

#[cfg(test)]
mod strand_tests {
    use super::*;

    #[test]
    fn a_transcript_takes_the_strand_of_its_own_gene_at_a_shared_vertex() {
        // one vertex, a gene on each strand; whatever order the index's rows come in
        for rows in [[(7, Strand::Plus), (9, Strand::Minus)], [(9, Strand::Minus), (7, Strand::Plus)]] {
            let mut v2g = Vertex2Gene::default();
            for r in rows { v2g.entry(100).or_default().insert(r); }
            assert_eq!(strand_at_start(&v2g, 100, 7), Strand::Plus);
            assert_eq!(strand_at_start(&v2g, 100, 9), Strand::Minus);
            // a gene the vertex does not list: the lowest gene's strand, in either order
            assert_eq!(strand_at_start(&v2g, 100, 8), Strand::Plus);
        }
        assert_eq!(strand_at_start(&Vertex2Gene::default(), 100, 7), Strand::Plus);
    }
}

#[cfg(test)]
mod renumbering_tests {
    use super::assign_dense_ids;

    #[test]
    fn ids_follow_the_entities_not_the_order_they_were_made_in() {
        // the same three vertices from two runs: other provisional IDs, another order of arrival
        let one = [(907, ("chr1", 500)), (901, ("chr1", 100)), (950, ("chr2", 7))];
        let other = [(11, ("chr2", 7)), (12, ("chr1", 500)), (13, ("chr1", 100))];
        let (a, b) = (assign_dense_ids(&one, 1000), assign_dense_ids(&other, 1000));
        assert_eq!((a[&901], a[&907], a[&950]), (1001, 1002, 1003));
        assert_eq!((b[&13], b[&12], b[&11]), (1001, 1002, 1003));
    }

    #[test]
    fn an_id_listed_twice_is_numbered_once_and_equal_keys_keep_their_order() {
        let items = [(5, 20), (7, 10), (5, 20), (9, 10)];
        let ids = assign_dense_ids(&items, 0);
        assert_eq!(ids.len(), 3);
        assert_eq!((ids[&7], ids[&9], ids[&5]), (1, 2, 3));
    }
}
