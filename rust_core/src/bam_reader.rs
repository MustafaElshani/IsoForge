// BAM/SAM reading with integrated preprocessing: QC filtering, microindel and
// splice junction correction, fraction-A and PAS assessment, poly(A) tail and
// cap detection, all in one streaming pass.

use rust_htslib::bam::{self, record::Aux, Read as BamRead};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use crate::errors::{IsoForgeError, IsoForgeResult};
use crate::genome_ref::{GenomeRef, SpliceJunctionDb};
use crate::types::*;

/// Parsed read from a BAM record with all fields needed for annotation and
/// TSS/TES clustering.
#[derive(Debug, Clone)]
pub struct ParsedRead {
    pub read_id: String,
    /// Dense run-stable index, assigned in BAM order after parsing.
    pub read_idx: ReadIdx,
    pub dataset: String,
    pub chromosome: String,
    pub strand: Strand,
    pub start: i64,     // 1-based leftmost
    pub end: i64,       // 1-based rightmost
    pub read_length: i64,
    pub positions: Vec<i64>, // [start, splice_sites..., end] in 5'->3' order
    /// Untemplated G's in the 5' soft clip (see `cap_signature`).
    pub cap_signature: bool,
    /// At least `UNALIGNED_3P_MIN` bases of the 3' soft clip precede its poly(A) run (see `unaligned_before_tail`).
    pub unaligned_3p: bool,
    /// Internal splice sites (bit i = positions[i + 1]) whose read sequence fits a
    /// nearby annotated alternative intron within one edit.
    pub ambiguous_sites: u64,
    /// Splice sites of unannotated introns the read's own bases favour over annotated ones nearby (bit i: splice
    /// site i in transcript order).
    pub sequence_settled: u64,
    pub n_exons: i32,
    pub cigar_string: String,
    pub flag: u16,
    // Internal priming label (from the genome, or from the fA tag)
    pub fraction_as: Option<f64>,
    // Custom IsoForge tags
    pub custom_label: Option<String>,
    pub allelic_label: Option<String>,
    pub start_support: Option<String>,
    pub end_support: Option<String>,
    // Splice junction correction stats
    pub sj_corrected: u32,       // number of splice junctions corrected
    pub sj_noncanonical: u32,    // remaining non-canonical after correction
    pub sj_canonical: u32,       // number of canonical splice junctions
    pub microindels_ins: u32,    // number of small insertions absorbed
    pub microindels_del: u32,    // number of small deletions absorbed
    // QC metrics computed during parsing
    pub alignment_coverage: f64,
    pub alignment_identity: Option<f64>,
    // ── Soft-clip and polyA tail fields (for TSS/TES clustering) ──────────
    /// Number of soft-clipped bases at the 5' end of the alignment.
    pub soft_clip_5p: u32,
    /// Number of soft-clipped bases at the 3' end of the alignment.
    pub soft_clip_3p: u32,
    /// PolyA tail length from the basecaller pt:i tag; soft-clip scoring is the fallback.
    pub polya_tail_len: Option<i64>,
    /// PolyA tail length detected from soft-clipped bases
    /// (+1 A/T, -5 mismatch, stop at -20).
    pub polya_tail_softclip: Option<i64>,
    // ── PAS-aware internal priming assessment ─────────────────────────────
    /// Internal priming result with PAS evidence.
    pub internal_priming: Option<InternalPrimingResult>,
    /// MAPQ below `min_mapq` but at least `quant_min_mapq`: the read builds no
    /// model and is kept for quantification only.
    pub quant_only: bool,
    /// The unspliced read is held by an annotated small RNA on its strand (see `small_rna.rs`).
    pub small_rna: bool,
}

// ═════════════════════════════════════════════════════════════════════════════
// CLUSTERING INPUT — per-read data for TSS/TES clustering
// ═════════════════════════════════════════════════════════════════════════════

/// Fields needed by TSS/TES clustering, implemented by `ParsedRead` and
/// `ObsCore`.
pub trait ClusteringInput: Sync + Send {
    fn read_idx(&self) -> ReadIdx;
    fn chromosome(&self) -> &str;
    fn strand_val(&self) -> Strand;
    fn pos_5p(&self) -> Option<i64>;
    fn pos_3p(&self) -> Option<i64>;
    fn dataset_name(&self) -> &str;
    fn polya_len(&self) -> Option<i64>;
    fn polya_softclip(&self) -> Option<i64>;
    fn is_ip(&self) -> bool;
    fn frac_a(&self) -> Option<f64>;
    fn n_exons(&self) -> i32;
    /// See `ParsedRead::unaligned_3p`.
    fn unaligned_3p(&self) -> bool;
}

impl ClusteringInput for ParsedRead {
    #[inline] fn read_idx(&self) -> ReadIdx { self.read_idx }
    #[inline] fn chromosome(&self) -> &str { &self.chromosome }
    #[inline] fn strand_val(&self) -> Strand { self.strand }
    #[inline] fn pos_5p(&self) -> Option<i64> { self.positions.first().copied() }
    #[inline] fn pos_3p(&self) -> Option<i64> { self.positions.last().copied() }
    #[inline] fn dataset_name(&self) -> &str { &self.dataset }
    #[inline] fn polya_len(&self) -> Option<i64> { self.polya_tail_len }
    #[inline] fn polya_softclip(&self) -> Option<i64> { self.polya_tail_softclip }
    #[inline] fn is_ip(&self) -> bool {
        self.internal_priming.as_ref().map_or(false, |ip| ip.is_internal_priming)
    }
    #[inline] fn frac_a(&self) -> Option<f64> { self.fraction_as }
    #[inline] fn n_exons(&self) -> i32 { self.n_exons }
    #[inline] fn unaligned_3p(&self) -> bool { self.unaligned_3p }
}

/// Quality check result
#[allow(dead_code)]
pub struct QcResult {
    pub passed: bool,
    pub metrics: QcMetrics,
}

/// Configuration for the integrated reader pipeline.
#[derive(Debug, Clone)]
pub struct ReaderPipeline {
    pub microexons: Option<Arc<crate::microexon_repair::MicroexonIndex>>,
    /// Annotated introns for sequence-based realignment of shifted splice sites.
    pub splice_realign: Option<Arc<crate::splice_realignment::SpliceRealignIndex>>,
    /// Terminal exons placed from soft clips (see `terminal_rescue`).
    pub terminal_rescue: Option<Arc<crate::terminal_rescue::TerminalExonIndex>>,
    /// Junction groups decided by the reads' pair-HMM likelihoods (see `junction_groups`); while they are counted or
    /// scored, parsing stops after this step.
    pub junction_groups: Option<Arc<crate::junction_groups::JunctionGroups>>,
    /// Annotated small RNAs whose reads pass the length and mapping-quality gates (see `small_rna`).
    pub small_rna: Option<Arc<crate::small_rna::SmallRnaIndex>>,
    /// Reference genome for fraction-A and splice motif checking; without it the
    /// fA tag is read from the BAM record.
    pub genome_path: Option<PathBuf>,

    /// Known splice junction database, preferred targets for splice junction correction.
    pub sj_db: Option<Arc<SpliceJunctionDb>>,

    /// Strand lookup for tail-less unspliced reads of an unoriented library.
    pub strand_index: Option<Arc<crate::strand_evidence::StrandIndex>>,

    /// Whether to correct non-canonical splice junctions.
    pub correct_sj: bool,

    /// Maximum base-pair offset for splice junction correction.
    pub max_sj_offset: i64,

    /// Range size for fraction-A computation.
    pub frac_a_range: i64,

    /// Whether to compute MD tag from genome if not present in BAM.
    pub compute_md_if_missing: bool,

    /// Correct microindels in the CIGAR string before extracting intron positions.
    pub correct_microindels: bool,

    /// Maximum size (bp) of insertions and deletions to correct.
    pub max_len_indel: u32,
}

impl Default for ReaderPipeline {
    fn default() -> Self {
        Self {
            microexons: None,
            splice_realign: None,
            terminal_rescue: None,
            junction_groups: None,
            small_rna: None,
            genome_path: None,
            sj_db: None,
            strand_index: None,
            correct_sj: true,
            max_sj_offset: 5,
            frac_a_range: 20,
            compute_md_if_missing: true,
            correct_microindels: true,
            max_len_indel: 5,
        }
    }
}

/// BAM reader with integrated preprocessing: streams BAM records and in a single pass performs QC filtering (flag,
/// coverage, identity, length), intron extraction, splice junction motif checking and correction, fraction-A
/// computation for internal priming, and MD tag computation when missing.
pub struct IntegratedBamReader {
    path: String,
    ref_names: Vec<String>,
    genome: Option<GenomeRef>,
    pipeline: ReaderPipeline,
    /// Segments of chimeric reads (records with an SA tag), when `collect_chimeric` is set.
    chimeric: std::sync::Mutex<Vec<crate::fusion::ChimericSegment>>,
    pub collect_chimeric: bool,
    /// Gene bodies: with them a spliced record is also cut at an intron joining two genes for fusion calling.
    pub gene_spans: Option<Arc<crate::fusion::GeneSpans>>,
}

impl IntegratedBamReader {
    /// Open a BAM/SAM file with optional genome reference for integrated processing.
    pub fn new(path: &str, pipeline: ReaderPipeline) -> IsoForgeResult<Self> {
        if !Path::new(path).exists() {
            return Err(IsoForgeError::BamError(format!(
                "BAM/SAM file does not exist: {}",
                path
            )));
        }

        // Open BAM to read header for reference names
        let reader = bam::Reader::from_path(path)
            .map_err(|e| IsoForgeError::HtsError(format!("Failed to open BAM: {}", e)))?;

        let header_view = reader.header();
        let ref_names: Vec<String> = (0..header_view.target_count())
            .map(|tid| String::from_utf8_lossy(header_view.tid2name(tid)).to_string())
            .collect();

        // Open genome reference if path is provided
        let genome = if let Some(ref genome_path) = pipeline.genome_path {
            Some(GenomeRef::open(genome_path)?)
        } else {
            None
        };

        Ok(Self {
            path: path.to_string(),
            ref_names,
            genome,
            pipeline,
            chimeric: std::sync::Mutex::new(Vec::new()),
            collect_chimeric: false,
            gene_spans: None,
        })
    }

    /// The chimeric segments collected by `process_windows`.
    pub fn take_chimeric(&self) -> Vec<crate::fusion::ChimericSegment> {
        std::mem::take(&mut *self.chimeric.lock().unwrap())
    }

    /// Parse the BAM and pass the reads back in windows of at most
    /// `window_size`. A reader thread reads records into batches and rayon
    /// workers process them in parallel, each with its own GenomeRef handle.
    pub fn process_windows<F>(
        &self,
        config: &RunConfig,
        dataset_name: &str,
        window_size: usize,
        mut on_window: F,
    ) -> IsoForgeResult<()>
    where
        F: FnMut(Vec<ParsedRead>) -> IsoForgeResult<()>,
    {
        use rayon::prelude::*;
        use std::sync::atomic::{AtomicU64, Ordering};

        let n_threads = rayon::current_num_threads();
        let process_batch_size = 5000; // records per parallel batch

        // Open per-thread genome references
        let thread_genomes: Vec<Option<GenomeRef>> = if self.genome.is_some() {
            (0..n_threads)
                .map(|_| self.genome.as_ref().and_then(|g| g.reopen().ok()))
                .collect()
        } else {
            (0..n_threads).map(|_| None).collect()
        };

        // Per-thread readers share ref_names and pipeline config but have their own
        // genome handles; each thread locks only its own slot.
        struct ThreadReader {
            ref_names: Vec<String>,
            genome: Option<GenomeRef>,
            pipeline: ReaderPipeline,
        }

        let thread_readers: Vec<std::sync::Mutex<ThreadReader>> = thread_genomes
            .into_iter()
            .map(|genome| std::sync::Mutex::new(ThreadReader {
                ref_names: self.ref_names.clone(),
                genome,
                pipeline: self.pipeline.clone(),
            }))
            .collect();

        let mut reader = bam::Reader::from_path(&self.path)
            .map_err(|e| IsoForgeError::HtsError(format!("Failed to open BAM: {}", e)))?;

        // htslib decompression threads
        let decomp_threads = (config.threads / 4).max(1).min(4);
        reader
            .set_threads(decomp_threads)
            .map_err(|e| IsoForgeError::HtsError(format!("Failed to set threads: {}", e)))?;

        let stats_total = AtomicU64::new(0);
        let stats_passed = AtomicU64::new(0);
        let stats_small_rna = AtomicU64::new(0);
        let stats_filtered = AtomicU64::new(0);
        let stats_errors = AtomicU64::new(0);
        let stats_sj_corrected = AtomicU64::new(0);
        let stats_sj_noncanonical = AtomicU64::new(0);
        let stats_mi_ins = AtomicU64::new(0);
        let stats_mi_del = AtomicU64::new(0);

        let mut window: Vec<ParsedRead> = Vec::new();

        // A reader thread decodes the next batches while this one parses the current. Batches arrive in file order.
        let (tx, rx) = std::sync::mpsc::sync_channel::<Result<Vec<bam::Record>, String>>(3);
        std::thread::scope(|scope| -> IsoForgeResult<()> {
            scope.spawn(move || {
                let mut record = bam::Record::new();
                loop {
                    let mut batch: Vec<bam::Record> = Vec::with_capacity(process_batch_size);
                    let mut end = false;
                    while batch.len() < process_batch_size {
                        match reader.read(&mut record) {
                            Some(Ok(())) => batch.push(record.clone()),
                            Some(Err(e)) => { let _ = tx.send(Err(format!("Failed to read record: {}", e))); return; }
                            None => { end = true; break; }
                        }
                    }
                    if !batch.is_empty() && tx.send(Ok(batch)).is_err() { return; }
                    if end { return; }
                }
            });
            for received in rx {
                let record_batch = received.map_err(IsoForgeError::HtsError)?;
                let batch_count = record_batch.len();

                stats_total.fetch_add(batch_count as u64, Ordering::Relaxed);

                // Segments of chimeric reads, for fusion calling.
                if self.collect_chimeric {
                    let segs: Vec<crate::fusion::ChimericSegment> = record_batch
                        .par_iter()
                        .flat_map_iter(|rec| match &self.gene_spans {
                            Some(gs) => crate::fusion::ChimericSegment::split_at_gene_joins(rec, &self.ref_names, dataset_name, gs),
                            None => crate::fusion::ChimericSegment::from_record(rec, &self.ref_names, dataset_name).into_iter().collect(),
                        })
                        .collect();
                    if !segs.is_empty() {
                        self.chimeric.lock().unwrap().extend(segs);
                    }
                }

                // ── Process records in parallel with rayon ──
                let batch_results: Vec<Option<ParsedRead>> = record_batch
                    .par_iter()
                    .map(|rec| {
                        // Get this thread's genome reference (each thread uses its own slot)
                        let tid = rayon::current_thread_index().unwrap_or(0);
                        let tr = thread_readers[tid % thread_readers.len()].lock().unwrap();

                        // Per-thread reader context
                        let result = process_record_with_refs(
                            rec,
                            &tr.ref_names,
                            tr.genome.as_ref(),
                            &tr.pipeline,
                            config,
                            dataset_name,
                        );

                        match result {
                            Ok(Some(mut parsed)) => {
                                stats_passed.fetch_add(1, Ordering::Relaxed);
                                if parsed.small_rna { stats_small_rna.fetch_add(1, Ordering::Relaxed); }
                                stats_sj_corrected.fetch_add(parsed.sj_corrected as u64, Ordering::Relaxed);
                                stats_sj_noncanonical.fetch_add(parsed.sj_noncanonical as u64, Ordering::Relaxed);
                                stats_mi_ins.fetch_add(parsed.microindels_ins as u64, Ordering::Relaxed);
                                stats_mi_del.fetch_add(parsed.microindels_del as u64, Ordering::Relaxed);
                                parsed.cigar_string = String::new();
                                Some(parsed)
                            }
                            Ok(None) => {
                                stats_filtered.fetch_add(1, Ordering::Relaxed);
                                None
                            }
                            Err(_) => {
                                stats_errors.fetch_add(1, Ordering::Relaxed);
                                None
                            }
                        }
                    })
                    .collect();

                // Collect results
                window.extend(batch_results.into_iter().flatten());
                if window.len() >= window_size {
                    on_window(std::mem::take(&mut window))?;
                }
            }
            Ok(())
        })?;

        if !window.is_empty() {
            on_window(window)?;
        }

        let total = stats_total.load(Ordering::Relaxed);
        let passed = stats_passed.load(Ordering::Relaxed);
        let filtered = stats_filtered.load(Ordering::Relaxed);
        let errors = stats_errors.load(Ordering::Relaxed);
        let sj_corr = stats_sj_corrected.load(Ordering::Relaxed);
        let sj_nc = stats_sj_noncanonical.load(Ordering::Relaxed);
        let mi_i = stats_mi_ins.load(Ordering::Relaxed);
        let mi_d = stats_mi_del.load(Ordering::Relaxed);

        // the junction groups' pre-passes print no summary
        if self.pipeline.junction_groups.as_ref().is_some_and(|j| j.is_prepass()) {
            return Ok(());
        }
        eprintln!(
            "[IsoForge] Processed {} records: {} passed QC, {} filtered, {} errors ({} threads)",
            total, passed, filtered, errors, n_threads
        );
        if self.pipeline.small_rna.is_some() {
            eprintln!("[IsoForge] Small RNAs: {} reads of annotated small RNAs", stats_small_rna.load(Ordering::Relaxed));
        }
        if self.genome.is_some() {
            eprintln!(
                "[IsoForge] Splice junctions: {} corrected, {} non-canonical remaining",
                sj_corr, sj_nc
            );
        }
        if self.pipeline.correct_microindels && (mi_i > 0 || mi_d > 0) {
            eprintln!(
                "[IsoForge] Microindel correction: {} insertions, {} deletions absorbed (≤{}bp)",
                mi_i, mi_d, self.pipeline.max_len_indel
            );
        }

        Ok(())
    }

    /// Stream reads with a callback.
    #[allow(dead_code)]
    pub fn stream_reads<F>(
        &self,
        config: &RunConfig,
        dataset_name: &str,
        mut callback: F,
    ) -> IsoForgeResult<ProcessingStats>
    where
        F: FnMut(ParsedRead) -> IsoForgeResult<()>,
    {
        let mut reader = bam::Reader::from_path(&self.path)
            .map_err(|e| IsoForgeError::HtsError(format!("Failed to open BAM: {}", e)))?;

        reader
            .set_threads(config.threads.max(1))
            .map_err(|e| IsoForgeError::HtsError(format!("Failed to set threads: {}", e)))?;

        let mut record = bam::Record::new();
        let mut stats = ProcessingStats::default();

        while let Some(result) = reader.read(&mut record) {
            result.map_err(|e| IsoForgeError::HtsError(format!("Failed to read record: {}", e)))?;
            stats.total_records += 1;

            match self.process_record(&record, config, dataset_name) {
                Ok(Some(parsed)) => {
                    stats.passed_qc += 1;
                    stats.sj_corrected += parsed.sj_corrected as u64;
                    callback(parsed)?;
                }
                Ok(None) => {
                    stats.filtered += 1;
                }
                Err(_) => {
                    stats.errors += 1;
                }
            }
        }

        Ok(stats)
    }

    /// Read reads in batches of at most `batch_size`.
    #[allow(dead_code)]
    pub fn read_batches(
        &self,
        config: &RunConfig,
        dataset_name: &str,
        batch_size: usize,
    ) -> IsoForgeResult<Vec<Vec<ParsedRead>>> {
        let mut reader = bam::Reader::from_path(&self.path)
            .map_err(|e| IsoForgeError::HtsError(format!("Failed to open BAM: {}", e)))?;

        reader
            .set_threads(config.threads.max(1))
            .map_err(|e| IsoForgeError::HtsError(format!("Failed to set threads: {}", e)))?;

        let mut batches = Vec::new();
        let mut current_batch = Vec::with_capacity(batch_size);
        let mut record = bam::Record::new();

        while let Some(result) = reader.read(&mut record) {
            result.map_err(|e| IsoForgeError::HtsError(format!("Failed to read record: {}", e)))?;

            if let Ok(Some(parsed)) = self.process_record(&record, config, dataset_name) {
                current_batch.push(parsed);
                if current_batch.len() >= batch_size {
                    batches.push(std::mem::replace(
                        &mut current_batch,
                        Vec::with_capacity(batch_size),
                    ));
                }
            }
        }

        if !current_batch.is_empty() {
            batches.push(current_batch);
        }

        Ok(batches)
    }

    /// Process a single BAM record: filter by flag, MAPQ, read length, alignment coverage and identity; correct
    /// microindels and extract introns from the CIGAR; check and correct splice junction motifs; build the positions
    /// array; compute fraction-A and internal priming; read soft clips, the poly(A) tail and custom SAM tags.
    fn process_record(
        &self,
        record: &bam::Record,
        config: &RunConfig,
        dataset_name: &str,
    ) -> IsoForgeResult<Option<ParsedRead>> {
        // ── Step 1: Flag filtering ───────────────────────────────────────
        if record.is_unmapped() || record.is_secondary() || record.is_supplementary() {
            return Ok(None);
        }
        let flag = record.flags();
        if flag != 0 && flag != 16 {
            return Ok(None);
        }

        // ── Step 1b: MAPQ filtering ──────────────────────────────────────
        let mapq = record.mapq();
        if mapq < config.min_mapq {
            return Ok(None);
        }

        let tid = record.tid();
        if tid < 0 || (tid as usize) >= self.ref_names.len() {
            return Ok(None);
        }

        let chromosome = self.ref_names[tid as usize].clone();
        let mut strand = if record.is_reverse() {
            Strand::Minus
        } else {
            Strand::Plus
        };

        let read_id = String::from_utf8_lossy(record.qname()).to_string();
        let dataset = if !dataset_name.is_empty() {
            dataset_name.to_string()
        } else if config.use_cb_tag {
            get_string_tag(record, b"CB").unwrap_or_default()
        } else {
            get_string_tag(record, b"RG").unwrap_or_default()
        };

        // ── Step 2: Length filtering ─────────────────────────────────────
        let sam_start = record.pos() + 1; // 0-based to 1-based
        let sam_end = record.cigar().end_pos(); // already 1-based end
        let read_length = record.seq_len() as i64;

        if read_length < config.min_length {
            return Ok(None);
        }

        // ── Step 3: Coverage from CIGAR ──────────────────────────────────
        let cigar = record.cigar();
        let Some(excluded) = coverage_exclusion(&record.seq().as_bytes(), &cigar, config) else { return Ok(None) };
        let coverage = compute_alignment_coverage_excluding(&cigar, excluded);
        if coverage < config.min_coverage {
            return Ok(None);
        }

        let Some(resolved) = resolve_strand(record.qname(), strand, &record.seq().as_bytes(), &cigar, &chromosome, sam_start, sam_end,
                                            self.genome.as_ref(), self.pipeline.strand_index.as_deref(), config) else { return Ok(None) };
        strand = resolved;

        // ── Step 4: Identity from MD tag ─────────────────────────────────
        let mut md_tag = get_string_tag(record, b"MD");

        // Compute the MD tag from the genome when missing
        if md_tag.is_none() && self.genome.is_some() && self.pipeline.compute_md_if_missing {
            if let Some(ref genome) = self.genome {
                let cigar_ops = extract_cigar_ops(&cigar);
                let seq = record.seq().as_bytes();
                if let Ok(computed_md) =
                    genome.compute_md_tag(&chromosome, sam_start, &cigar_ops, &seq)
                {
                    md_tag = Some(computed_md);
                }
            }
        }

        let identity = md_tag
            .as_ref()
            .map(|md| compute_alignment_identity(md, record.seq_len()));

        if let Some(id) = identity {
            if id < config.min_identity {
                return Ok(None);
            }
        }

        // ── Step 5: Microindel correction + intron extraction ──────────
        // Small insertions/deletions (≤ max_len_indel bp) are absorbed into adjacent
        // match operations before introns are extracted.
        let (intron_list, mi_ins, mi_del) = if self.pipeline.correct_microindels {
            let (cleaned_cigar, n_ins, n_del) = correct_microindels_in_cigar(
                &cigar,
                self.pipeline.max_len_indel,
            );
            let introns = get_introns_from_corrected_cigar(&cleaned_cigar, sam_start);
            (introns, n_ins, n_del)
        } else {
            (get_introns_from_cigar(&cigar, sam_start), 0, 0)
        };

        let intron_list = match (&self.pipeline.microexons, &self.genome) {
            (Some(index), Some(genome)) => index.repair(record, &chromosome, strand, genome, &intron_list).0,
            _ => intron_list,
        };
        let capped = cap_signature(record, strand, self.genome.as_ref(), &chromosome);
        let (intron_list, n_realigned, ambiguous_sites, sequence_settled) = match (&self.pipeline.splice_realign, &self.genome) {
            (Some(index), Some(genome)) => {
                let r = index.realign(record, &chromosome, strand, genome, &intron_list);
                (r.introns, r.moved, crate::splice_realignment::ambiguous_site_mask(&r.ambiguous, strand),
                 crate::splice_realignment::ambiguous_site_mask(&r.kept, strand))
            }
            _ => (intron_list, 0, 0, 0),
        };

        // ── Step 6: Splice junction correction ───────────────────────────
        let (corrected_introns, sj_corrected, sj_noncanonical, sj_canonical) =
            if self.genome.is_some() && !intron_list.is_empty() {
                self.correct_and_check_introns(
                    &chromosome,
                    &intron_list,
                    strand,
                )?
            } else {
                (intron_list.clone(), 0, 0, intron_list.len() as u32 / 2)
            };

        // ── Step 6b: Terminal exons left in a soft clip ──────────────────
        let (mut corrected_introns, mut sam_start, mut sam_end) = (corrected_introns, sam_start, sam_end);
        let (mut capped, mut ambiguous_sites, mut sequence_settled) = (capped, ambiguous_sites, sequence_settled);
        let (leading_clip, trailing_clip) = match (&self.pipeline.terminal_rescue, &self.genome) {
            (Some(index), Some(genome)) => {
                let r = index.rescue(record, &chromosome, strand, genome, &mut corrected_introns, &mut sam_start, &mut sam_end);
                shift_ambiguous(&mut ambiguous_sites, strand, r.low, r.high);
                shift_ambiguous(&mut sequence_settled, strand, r.low, r.high);
                if let Some(c) = r.capped { capped = c; }
                (r.leading_clip, r.trailing_clip)
            }
            _ => extract_soft_clips(&cigar),
        };

        // ── Step 7: Build positions array ────────────────────────────────
        let splice_sites: Vec<i64> = corrected_introns
            .iter()
            .enumerate()
            .map(|(i, &x)| {
                if i % 2 == 1 {
                    x + 1 // intron end → exon start
                } else {
                    x - 1 // intron start → exon end
                }
            })
            .collect();

        let mut positions = Vec::with_capacity(splice_sites.len() + 2);
        positions.push(sam_start);
        positions.extend_from_slice(&splice_sites);
        positions.push(sam_end);

        if strand == Strand::Minus {
            positions.reverse();
        }

        let n_exons = (positions.len() / 2) as i32;

        // ── Step 8: Fraction-A + PAS-aware internal priming ──────────────
        let fraction_as = self.compute_fraction_a_for_read(
            &chromosome,
            sam_start,
            sam_end,
            strand,
            record,
        );

        // PAS-aware internal priming assessment
        let internal_priming = self.assess_internal_priming_for_read(
            &chromosome,
            sam_start,
            sam_end,
            strand,
            record,
            config.max_frac_as,
        );

        // ── Step 9: Soft-clip extraction + polyA tail detection ───────────
        // (clips as left after Step 6b)
        // Convert alignment-oriented clips to strand-aware 5'/3'
        let (soft_clip_5p, soft_clip_3p) = match strand {
            Strand::Plus => (leading_clip, trailing_clip),
            Strand::Minus => (trailing_clip, leading_clip),
        };

        // Poly(A) tail length from the basecaller pt:i tag
        let polya_tail_len = get_int_tag(record, b"pt");

        // Poly(A) tail from soft-clipped bases (fallback)
        let polya_tail_softclip = detect_polya_tail_from_softclip(
            &record.seq().as_bytes(),
            strand,
            leading_clip,
            trailing_clip,
        );
        let unaligned_3p = unaligned_before_tail(&record.seq().as_bytes(), strand, leading_clip, trailing_clip)
            >= UNALIGNED_3P_MIN;

        // ── Step 10: Parse custom SAM tags ───────────────────────────────
        let cigar_string = cigar_to_string(&cigar);

        Ok(Some(ParsedRead {
            read_id,
            read_idx: READ_IDX_NONE,
            dataset,
            chromosome,
            strand,
            start: sam_start,
            end: sam_end,
            read_length,
            positions,
            n_exons,
            cigar_string,
            flag,
            fraction_as,
            custom_label: get_string_tag(record, b"lC"),
            allelic_label: get_string_tag(record, b"lA"),
            start_support: get_string_tag(record, b"tS"),
            end_support: get_string_tag(record, b"tE"),
            sj_corrected: sj_corrected + n_realigned as u32,
            ambiguous_sites,
            sequence_settled,
            cap_signature: capped,
            unaligned_3p,
            sj_noncanonical,
            sj_canonical,
            microindels_ins: mi_ins,
            microindels_del: mi_del,
            alignment_coverage: coverage,
            alignment_identity: identity,
            soft_clip_5p,
            soft_clip_3p,
            polya_tail_len,
            polya_tail_softclip,
            internal_priming,
            quant_only: false,
            small_rna: false,
        }))
    }

    /// Check and optionally correct splice junctions using the genome reference: a non-canonical intron moves to a
    /// known junction within ±max_sj_offset when the database holds one, else to the nearest canonical motif in that
    /// range. Returns: (corrected_introns, n_corrected, n_remaining_noncanonical, n_canonical)
    fn correct_and_check_introns(
        &self,
        chrom: &str,
        introns: &[i64],
        strand: Strand,
    ) -> IsoForgeResult<(Vec<i64>, u32, u32, u32)> {
        let genome = self.genome.as_ref().unwrap();
        let mut corrected = introns.to_vec();
        let mut n_corrected: u32 = 0;
        let mut n_noncanonical: u32 = 0;
        let mut n_canonical: u32 = 0;

        // Process pairs (intron_start, intron_end)
        for i in (0..corrected.len()).step_by(2) {
            if i + 1 >= corrected.len() {
                break;
            }

            let intron_start = corrected[i];
            let intron_end = corrected[i + 1];

            // Check current motif
            let motif_result = genome.check_splice_motif(chrom, intron_start, intron_end, strand)?;

            if motif_result.is_canonical {
                n_canonical += 1;
                continue;
            }

            // Non-canonical — try correction if enabled
            if self.pipeline.correct_sj {
                // First check SJ database for a known junction nearby
                if let Some(ref sj_db) = self.pipeline.sj_db {
                    if let Some((new_start, new_end)) = self.find_known_sj_nearby(
                        sj_db,
                        chrom,
                        intron_start,
                        intron_end,
                        self.pipeline.max_sj_offset,
                    ) {
                        corrected[i] = new_start;
                        corrected[i + 1] = new_end;
                        n_corrected += 1;
                        n_canonical += 1;
                        continue;
                    }
                }

                // Fall back to genome-based correction
                if let Ok(Some((new_start, new_end))) = genome.correct_splice_junction(
                    chrom,
                    intron_start,
                    intron_end,
                    strand,
                    self.pipeline.max_sj_offset,
                ) {
                    corrected[i] = new_start;
                    corrected[i + 1] = new_end;
                    n_corrected += 1;
                    n_canonical += 1;
                    continue;
                }
            }

            // Couldn't correct
            n_noncanonical += 1;
        }

        Ok((corrected, n_corrected, n_noncanonical, n_canonical))
    }

    /// Search for a known splice junction within ±offset of the given position.
    fn find_known_sj_nearby(
        &self,
        sj_db: &SpliceJunctionDb,
        chrom: &str,
        intron_start: i64,
        intron_end: i64,
        max_offset: i64,
    ) -> Option<(i64, i64)> {
        let mut best: Option<(i64, i64, i64)> = None;

        for d_off in -max_offset..=max_offset {
            for a_off in -max_offset..=max_offset {
                let new_start = intron_start + d_off;
                let new_end = intron_end + a_off;
                if new_start <= 0 || new_end <= 0 || new_start >= new_end {
                    continue;
                }

                if let Some(info) = sj_db.get(chrom, new_start, new_end) {
                    let total_shift = d_off.abs() + a_off.abs();
                    let score = info.unique_reads + info.multi_reads;
                    match &best {
                        None => best = Some((new_start, new_end, total_shift)),
                        Some((_, _, prev_shift)) => {
                            if total_shift < *prev_shift
                                || (total_shift == *prev_shift && score > 0)
                            {
                                best = Some((new_start, new_end, total_shift));
                            }
                        }
                    }
                }
            }
        }

        best.map(|(s, e, _)| (s, e))
    }

    /// Compute fraction-A for internal priming detection, from the genome when
    /// loaded, otherwise from the fA tag.
    fn compute_fraction_a_for_read(
        &self,
        chrom: &str,
        sam_start: i64,
        sam_end: i64,
        strand: Strand,
        record: &bam::Record,
    ) -> Option<f64> {
        // Compute from the genome
        if let Some(ref genome) = self.genome {
            // transcript_end: rightmost aligned position on the + strand, leftmost on the - strand
            let transcript_end = match strand {
                Strand::Plus => sam_end,
                Strand::Minus => sam_start,
            };

            return genome.compute_fraction_a(
                chrom,
                transcript_end,
                strand,
                self.pipeline.frac_a_range,
            );
        }

        // Otherwise the fA tag
        get_float_tag(record, b"fA")
    }

    /// PAS-aware internal priming assessment: fraction-A and a PAS hexamer scan
    /// with a genome, otherwise fraction-A from the fA tag only.
    fn assess_internal_priming_for_read(
        &self,
        chrom: &str,
        sam_start: i64,
        sam_end: i64,
        strand: Strand,
        record: &bam::Record,
        frac_a_threshold: f64,
    ) -> Option<InternalPrimingResult> {
        // transcript_end: rightmost aligned position on the + strand, leftmost on the - strand
        let transcript_end = match strand {
            Strand::Plus => sam_end,
            Strand::Minus => sam_start,
        };

        if let Some(ref genome) = self.genome {
            // PAS-aware assessment
            Some(genome.assess_internal_priming(
                chrom,
                transcript_end,
                strand,
                self.pipeline.frac_a_range,
                frac_a_threshold,
            ))
        } else {
            // fA tag only (no PAS information)
            let fraction_a = get_float_tag(record, b"fA")?;
            Some(InternalPrimingResult {
                fraction_a,
                has_pas_signal: false,    // unknown without genome
                pas_motif: None,
                pas_distance: None,
                is_internal_priming: fraction_a > frac_a_threshold,
            })
        }
    }
}

/// Record processing for rayon worker threads: parses one BAM record with the reference names, genome handle and
/// pipeline of the calling thread.
fn process_record_with_refs(
    record: &bam::Record,
    ref_names: &[String],
    genome: Option<&GenomeRef>,
    pipeline: &ReaderPipeline,
    config: &RunConfig,
    dataset_name: &str,
) -> IsoForgeResult<Option<ParsedRead>> {
    // ── Step 1: Flag filtering ───────────────────────────────────────
    if record.is_unmapped() || record.is_secondary() || record.is_supplementary() {
        return Ok(None);
    }
    let flag = record.flags();
    if flag != 0 && flag != 16 {
        return Ok(None);
    }

    let tid = record.tid();
    if tid < 0 || (tid as usize) >= ref_names.len() {
        return Ok(None);
    }

    let chromosome = ref_names[tid as usize].clone();

    // An unspliced alignment lying inside an annotated small RNA passes the mapping-quality and length gates below
    // (small_rna.rs); whether the read is on the small RNA's strand is settled once its strand is resolved.
    let held = match &pipeline.small_rna {
        Some(index) if !record.cigar().iter().any(|op| matches!(op, bam::record::Cigar::RefSkip(_))) =>
            index.holding(&chromosome, record.pos() + 1, record.cigar().end_pos()),
        _ => crate::small_rna::Holding::default(),
    };

    // ── Step 1b: MAPQ filtering ──────────────────────────────────────
    // Reads in [quant_min_mapq, min_mapq) are parsed for quantification only.
    let mapq = record.mapq();
    let below_mapq = mapq < config.min_mapq;
    let no_quant = config.quant_min_mapq >= config.min_mapq || mapq < config.quant_min_mapq;
    if below_mapq && no_quant && !held.any() {
        return Ok(None);
    }
    // The junction groups' counting pass needs only spliced reads, their scoring pass only reads with an intron near
    // a contested site.
    if let Some(jg) = &pipeline.junction_groups {
        let skip = match jg.phase() {
            crate::junction_groups::Phase::Count => !record.cigar().iter().any(|op| matches!(op, bam::record::Cigar::RefSkip(_))),
            crate::junction_groups::Phase::Score => !jg.near_contested(&chromosome, record),
            crate::junction_groups::Phase::Decide => false,
        };
        if skip { return Ok(None); }
    }
    let mut strand = if record.is_reverse() {
        Strand::Minus
    } else {
        Strand::Plus
    };

    let read_id = String::from_utf8_lossy(record.qname()).to_string();
    let dataset = if !dataset_name.is_empty() {
        dataset_name.to_string()
    } else if config.use_cb_tag {
        get_string_tag(record, b"CB").unwrap_or_default()
    } else {
        get_string_tag(record, b"RG").unwrap_or_default()
    };

    // ── Step 2: Length filtering ─────────────────────────────────────
    let sam_start = record.pos() + 1;
    let sam_end = record.cigar().end_pos();
    let read_length = record.seq_len() as i64;

    let min_length = if held.any() { config.min_length.min(config.small_rna_min_read_length) } else { config.min_length };
    if read_length < min_length {
        return Ok(None);
    }

    // ── Step 3: Coverage from CIGAR ──────────────────────────────────
    let cigar = record.cigar();
    let Some(excluded) = coverage_exclusion(&record.seq().as_bytes(), &cigar, config) else { return Ok(None) };
    let coverage = compute_alignment_coverage_excluding(&cigar, excluded);
    if coverage < config.min_coverage {
        return Ok(None);
    }

    let Some(resolved) = resolve_strand(record.qname(), strand, &record.seq().as_bytes(), &cigar, &chromosome, sam_start, sam_end,
                                        genome, pipeline.strand_index.as_deref(), config) else { return Ok(None) };
    strand = resolved;

    // A read held by a small RNA of the other strand only is no small-RNA read: the gates it was let through apply.
    let held = held.on(strand);
    if !held && ((below_mapq && no_quant) || read_length < config.min_length) {
        return Ok(None);
    }

    // ── Step 4: Identity from MD tag ─────────────────────────────────
    let mut md_tag = get_string_tag(record, b"MD");

    if md_tag.is_none() && genome.is_some() && pipeline.compute_md_if_missing {
        if let Some(genome_ref) = genome {
            let cigar_ops = extract_cigar_ops(&cigar);
            let seq = record.seq().as_bytes();
            if let Ok(computed_md) =
                genome_ref.compute_md_tag(&chromosome, sam_start, &cigar_ops, &seq)
            {
                md_tag = Some(computed_md);
            }
        }
    }

    let identity = md_tag
        .as_ref()
        .map(|md| compute_alignment_identity(md, record.seq_len()));

    if let Some(id) = identity {
        if id < config.min_identity {
            return Ok(None);
        }
    }

    // ── Step 5: Microindel correction + intron extraction ──────────
    let (intron_list, mi_ins, mi_del) = if pipeline.correct_microindels {
        let (cleaned_cigar, n_ins, n_del) = correct_microindels_in_cigar(
            &cigar,
            pipeline.max_len_indel,
        );
        let introns = get_introns_from_corrected_cigar(&cleaned_cigar, sam_start);
        (introns, n_ins, n_del)
    } else {
        (get_introns_from_cigar(&cigar, sam_start), 0, 0)
    };

    let mut intron_list = intron_list;
    let n_indel_placed = match genome {
        Some(genome) => place_indel_adjacent_introns(record, &chromosome, strand, genome, &mut intron_list,
                                                      config.junction_indel_min, config.junction_indel_max),
        None => 0,
    };
    let intron_list = match (&pipeline.microexons, genome) {
        (Some(index), Some(genome)) => index.repair(record, &chromosome, strand, genome, &intron_list).0,
        _ => intron_list,
    };
    // Junction groups first: an intron in a contested group takes its group's decision; counting and scoring stop
    // the parse here.
    let (intron_list, grouped) = match (&pipeline.junction_groups, genome) {
        (Some(jg), Some(genome)) => {
            let (introns, handled) = jg.apply(record, &dataset, &chromosome, genome, &intron_list);
            if jg.is_prepass() { return Ok(None); }
            (introns, handled)
        }
        _ => { let n = intron_list.len() / 2; (intron_list, vec![false; n]) }
    };
    let capped = cap_signature(record, strand, genome, &chromosome);
    let (intron_list, n_realigned, ambiguous_sites, sequence_settled) = match (&pipeline.splice_realign, genome) {
        (Some(index), Some(genome)) => {
            let mut r = index.realign(record, &chromosome, strand, genome, &intron_list);
            // an intron the groups decided keeps the groups' coordinates; the realignment's ambiguity flag stays unless
            // the realignment had moved it
            if grouped.iter().any(|&g| g) && r.introns.len() == intron_list.len() {
                for (k, _) in grouped.iter().enumerate().filter(|(_, &g)| g) {
                    if (r.introns[2 * k], r.introns[2 * k + 1]) != (intron_list[2 * k], intron_list[2 * k + 1]) {
                        r.moved -= 1;
                        r.ambiguous[k] = false;
                    }
                    r.introns[2 * k] = intron_list[2 * k];
                    r.introns[2 * k + 1] = intron_list[2 * k + 1];
                    r.kept[k] = false;
                }
                if r.introns.windows(2).any(|w| w[0] > w[1]) { r.introns = intron_list.clone(); }
            }
            (r.introns, r.moved, crate::splice_realignment::ambiguous_site_mask(&r.ambiguous, strand),
             crate::splice_realignment::ambiguous_site_mask(&r.kept, strand))
        }
        _ => (intron_list, 0, 0, 0),
    };

    // ── Step 6: Splice junction correction ───────────────────────────
    let (corrected_introns, sj_corrected, sj_noncanonical, sj_canonical) =
        if genome.is_some() && !intron_list.is_empty() {
            correct_and_check_introns_standalone(
                genome.unwrap(),
                pipeline.sj_db.as_deref(),
                &chromosome,
                &intron_list,
                strand,
                pipeline.correct_sj,
                pipeline.max_sj_offset,
            )?
        } else {
            (intron_list.clone(), 0, 0, intron_list.len() as u32 / 2)
        };

    // ── Step 6b: Terminal exons left in a soft clip ──────────────────
    let (mut corrected_introns, mut sam_start, mut sam_end) = (corrected_introns, sam_start, sam_end);
    let (mut capped, mut ambiguous_sites, mut sequence_settled) = (capped, ambiguous_sites, sequence_settled);
    let (leading_clip, trailing_clip) = match (&pipeline.terminal_rescue, genome) {
        (Some(index), Some(genome)) => {
            let r = index.rescue(record, &chromosome, strand, genome, &mut corrected_introns, &mut sam_start, &mut sam_end);
            shift_ambiguous(&mut ambiguous_sites, strand, r.low, r.high);
                shift_ambiguous(&mut sequence_settled, strand, r.low, r.high);
            if let Some(c) = r.capped { capped = c; }
            (r.leading_clip, r.trailing_clip)
        }
        _ => extract_soft_clips(&cigar),
    };

    // ── Step 7: Build positions array ────────────────────────────────
    let splice_sites: Vec<i64> = corrected_introns
        .iter()
        .enumerate()
        .map(|(i, &x)| {
            if i % 2 == 1 { x + 1 } else { x - 1 }
        })
        .collect();

    let mut positions = Vec::with_capacity(splice_sites.len() + 2);
    positions.push(sam_start);
    positions.extend_from_slice(&splice_sites);
    positions.push(sam_end);

    if strand == Strand::Minus {
        positions.reverse();
    }

    let n_exons = (positions.len() / 2) as i32;

    // ── Step 8: Fraction-A + PAS-aware internal priming ──────────────
    let transcript_end = match strand {
        Strand::Plus => sam_end,
        Strand::Minus => sam_start,
    };

    let fraction_as = if let Some(genome_ref) = genome {
        genome_ref.compute_fraction_a(
            &chromosome,
            transcript_end,
            strand,
            pipeline.frac_a_range,
        )
    } else {
        get_float_tag(record, b"fA")
    };

    let internal_priming = if let Some(genome_ref) = genome {
        Some(genome_ref.assess_internal_priming(
            &chromosome,
            transcript_end,
            strand,
            pipeline.frac_a_range,
            config.max_frac_as,
        ))
    } else {
        get_float_tag(record, b"fA").map(|fraction_a| {
            InternalPrimingResult {
                fraction_a,
                has_pas_signal: false,
                pas_motif: None,
                pas_distance: None,
                is_internal_priming: fraction_a > config.max_frac_as,
            }
        })
    };

    // ── Step 9: Soft-clip extraction + polyA tail detection ───────────
    // (clips as left after Step 6b)
    let (soft_clip_5p, soft_clip_3p) = match strand {
        Strand::Plus => (leading_clip, trailing_clip),
        Strand::Minus => (trailing_clip, leading_clip),
    };

    let polya_tail_len = get_int_tag(record, b"pt");
    let polya_tail_softclip = detect_polya_tail_from_softclip(
        &record.seq().as_bytes(),
        strand,
        leading_clip,
        trailing_clip,
    );
    let unaligned_3p = unaligned_before_tail(&record.seq().as_bytes(), strand, leading_clip, trailing_clip)
        >= UNALIGNED_3P_MIN;

    // ── Step 10: Parse custom SAM tags ───────────────────────────────
    let cigar_string = cigar_to_string(&cigar);

    // A terminal exon placed from a soft clip makes the read spliced, and it is no small-RNA read then.
    let small_rna = held && n_exons == 1;
    if held && !small_rna && ((below_mapq && no_quant) || read_length < config.min_length) {
        return Ok(None);
    }
    let quant_only = below_mapq && !small_rna;

    Ok(Some(ParsedRead {
        read_id,
        read_idx: READ_IDX_NONE,
        dataset,
        chromosome,
        strand,
        start: sam_start,
        end: sam_end,
        read_length,
        positions,
        n_exons,
        cigar_string,
        flag,
        fraction_as,
        custom_label: get_string_tag(record, b"lC"),
        allelic_label: get_string_tag(record, b"lA"),
        start_support: get_string_tag(record, b"tS"),
        end_support: get_string_tag(record, b"tE"),
        sj_corrected: sj_corrected + n_realigned as u32 + n_indel_placed,
        ambiguous_sites,
        sequence_settled,
        cap_signature: capped,
        unaligned_3p,
        sj_noncanonical,
        sj_canonical,
        microindels_ins: mi_ins,
        microindels_del: mi_del,
        alignment_coverage: coverage,
        alignment_identity: identity,
        soft_clip_5p,
        soft_clip_3p,
        polya_tail_len,
        polya_tail_softclip,
        internal_priming,
        quant_only,
        small_rna,
    }))
}

/// True when the soft clip at the read's 5' end (in transcript orientation) ends in GG next to the first aligned
/// base and the genome does not read GG there. The clip may be any length.
fn cap_signature(record: &bam::Record, strand: Strand, genome: Option<&GenomeRef>, chrom: &str) -> bool {
    use rust_htslib::bam::record::Cigar::SoftClip;
    let cigar = record.cigar();
    let seq = record.seq().as_bytes();
    let clip: Vec<u8> = match (strand, cigar.first(), cigar.last()) {
        (Strand::Plus, Some(SoftClip(n)), _) => seq[..*n as usize].to_vec(),
        (Strand::Minus, _, Some(SoftClip(n))) => seq[seq.len() - *n as usize..].iter().rev().map(|&b| complement(b)).collect(),
        _ => return false,
    };
    // Only the two clip bases adjacent to the alignment are read.
    if clip.len() < 2 || !clip[clip.len() - 2..].eq_ignore_ascii_case(b"GG") {
        return false;
    }
    let Some(genome) = genome else { return true };
    let before = match strand {
        Strand::Plus => genome.fetch_seq(chrom, record.pos() - 2, record.pos()).ok(),
        Strand::Minus => genome.fetch_seq(chrom, cigar.end_pos(), cigar.end_pos() + 2).ok()
            .map(|s| s.iter().rev().map(|&b| complement(b)).collect()),
    };
    before.as_deref().map_or(true, |b| !b.eq_ignore_ascii_case(b"GG"))
}

fn complement(b: u8) -> u8 {
    match b.to_ascii_uppercase() { b'A' => b'T', b'C' => b'G', b'G' => b'C', b'T' => b'A', other => other }
}

/// Splice junction correction for worker threads.
fn correct_and_check_introns_standalone(
    genome: &GenomeRef,
    sj_db: Option<&SpliceJunctionDb>,
    chrom: &str,
    introns: &[i64],
    strand: Strand,
    correct_sj: bool,
    max_sj_offset: i64,
) -> IsoForgeResult<(Vec<i64>, u32, u32, u32)> {
    let mut corrected = introns.to_vec();
    let mut n_corrected: u32 = 0;
    let mut n_noncanonical: u32 = 0;
    let mut n_canonical: u32 = 0;

    for i in (0..corrected.len()).step_by(2) {
        if i + 1 >= corrected.len() {
            break;
        }

        let intron_start = corrected[i];
        let intron_end = corrected[i + 1];

        let motif_result = genome.check_splice_motif(chrom, intron_start, intron_end, strand)?;

        if motif_result.is_canonical {
            n_canonical += 1;
            continue;
        }

        if correct_sj {
            // Check SJ database first
            if let Some(db) = sj_db {
                if let Some((new_start, new_end)) = find_known_sj_nearby_standalone(
                    db, chrom, intron_start, intron_end, max_sj_offset,
                ) {
                    corrected[i] = new_start;
                    corrected[i + 1] = new_end;
                    n_corrected += 1;
                    n_canonical += 1;
                    continue;
                }
            }

            // Genome-based correction fallback
            if let Ok(Some((new_start, new_end))) = genome.correct_splice_junction(
                chrom, intron_start, intron_end, strand, max_sj_offset,
            ) {
                corrected[i] = new_start;
                corrected[i + 1] = new_end;
                n_corrected += 1;
                n_canonical += 1;
                continue;
            }
        }

        n_noncanonical += 1;
    }

    Ok((corrected, n_corrected, n_noncanonical, n_canonical))
}

/// Splice junction database search for worker threads.
fn find_known_sj_nearby_standalone(
    sj_db: &SpliceJunctionDb,
    chrom: &str,
    intron_start: i64,
    intron_end: i64,
    max_offset: i64,
) -> Option<(i64, i64)> {
    let mut best: Option<(i64, i64, i64)> = None;

    for d_off in -max_offset..=max_offset {
        for a_off in -max_offset..=max_offset {
            let new_start = intron_start + d_off;
            let new_end = intron_end + a_off;
            if new_start <= 0 || new_end <= 0 || new_start >= new_end {
                continue;
            }

            if let Some(info) = sj_db.get(chrom, new_start, new_end) {
                let total_shift = d_off.abs() + a_off.abs();
                let score = info.unique_reads + info.multi_reads;
                match &best {
                    None => best = Some((new_start, new_end, total_shift)),
                    Some((_, _, prev_shift)) => {
                        if total_shift < *prev_shift
                            || (total_shift == *prev_shift && score > 0)
                        {
                            best = Some((new_start, new_end, total_shift));
                        }
                    }
                }
            }
        }
    }

    best.map(|(s, e, _)| (s, e))
}

/// Processing statistics for a BAM file.
#[derive(Debug, Default, Clone)]
pub struct ProcessingStats {
    pub total_records: u64,
    pub passed_qc: u64,
    pub filtered: u64,
    pub errors: u64,
    pub sj_corrected: u64,
    pub sj_noncanonical: u64,
    pub microindels_insertions: u64,
    pub microindels_deletions: u64,
}

// ─── Basic reader (no integrated processing) ─────────────────────────────────

/// BamReader without integrated processing.
pub struct BamReader {
    path: String,
}

impl BamReader {
    pub fn new(path: &str) -> IsoForgeResult<Self> {
        if !Path::new(path).exists() {
            return Err(IsoForgeError::BamError(format!(
                "BAM file does not exist: {}",
                path
            )));
        }
        Ok(Self {
            path: path.to_string(),
        })
    }

    #[allow(dead_code)]
    pub fn read_all(&self, config: &RunConfig) -> IsoForgeResult<Vec<ParsedRead>> {
        let mut reader = bam::Reader::from_path(&self.path)
            .map_err(|e| IsoForgeError::HtsError(format!("Failed to open BAM: {}", e)))?;
        reader
            .set_threads(config.threads.max(1))
            .map_err(|e| IsoForgeError::HtsError(format!("Failed to set threads: {}", e)))?;

        let mut reads = Vec::new();
        let mut record = bam::Record::new();
        while let Some(result) = reader.read(&mut record) {
            result.map_err(|e| IsoForgeError::HtsError(format!("Failed to read record: {}", e)))?;
            if let Some(parsed) = self.parse_record(&record, config)? {
                reads.push(parsed);
            }
        }
        Ok(reads)
    }

    fn parse_record(
        &self,
        record: &bam::Record,
        config: &RunConfig,
    ) -> IsoForgeResult<Option<ParsedRead>> {
        if record.is_unmapped() || record.is_secondary() || record.is_supplementary() {
            return Ok(None);
        }
        let flag = record.flags();
        if flag != 0 && flag != 16 {
            return Ok(None);
        }

        // MAPQ filtering
        let mapq = record.mapq();
        if mapq < config.min_mapq {
            return Ok(None);
        }

        let read_id = String::from_utf8_lossy(record.qname()).to_string();
        let dataset = if config.use_cb_tag {
            get_string_tag(record, b"CB").unwrap_or_default()
        } else {
            get_string_tag(record, b"RG").unwrap_or_default()
        };

        let tid = record.tid();
        if tid < 0 {
            return Ok(None);
        }

        let strand = if record.is_reverse() {
            Strand::Minus
        } else {
            Strand::Plus
        };

        let sam_start = record.pos() + 1;
        let sam_end = record.cigar().end_pos();
        let read_length = record.seq_len() as i64;

        if read_length < config.min_length {
            return Ok(None);
        }

        let cigar = record.cigar();
        let Some(excluded) = coverage_exclusion(&record.seq().as_bytes(), &cigar, config) else { return Ok(None) };
        let coverage = compute_alignment_coverage_excluding(&cigar, excluded);
        if coverage < config.min_coverage {
            return Ok(None);
        }

        let md_tag = get_string_tag(record, b"MD");
        let identity = md_tag
            .as_ref()
            .map(|md| compute_alignment_identity(md, record.seq_len()));
        if let Some(id) = identity {
            if id < config.min_identity {
                return Ok(None);
            }
        }

        let intron_list = get_introns_from_cigar(&cigar, sam_start);
        let splice_sites: Vec<i64> = intron_list
            .iter()
            .enumerate()
            .map(|(i, &x)| if i % 2 == 1 { x + 1 } else { x - 1 })
            .collect();

        let mut positions = Vec::with_capacity(splice_sites.len() + 2);
        positions.push(sam_start);
        positions.extend_from_slice(&splice_sites);
        positions.push(sam_end);

        if strand == Strand::Minus {
            positions.reverse();
        }

        let n_exons = (positions.len() / 2) as i32;
        let chromosome = format!("tid_{}", tid);
        let cigar_string = cigar_to_string(&cigar);

        Ok(Some(ParsedRead {
            read_id,
            read_idx: READ_IDX_NONE,
            dataset,
            chromosome,
            strand,
            start: sam_start,
            end: sam_end,
            read_length,
            positions,
            n_exons,
            cigar_string,
            flag,
            fraction_as: get_float_tag(record, b"fA"),
            custom_label: get_string_tag(record, b"lC"),
            allelic_label: get_string_tag(record, b"lA"),
            start_support: get_string_tag(record, b"tS"),
            end_support: get_string_tag(record, b"tE"),
            sj_corrected: 0,
            ambiguous_sites: 0,
            sequence_settled: 0,
            cap_signature: false,
            unaligned_3p: false,
            sj_noncanonical: 0,
            sj_canonical: intron_list.len() as u32 / 2,
            microindels_ins: 0,
            microindels_del: 0,
            alignment_coverage: coverage,
            alignment_identity: identity,
            soft_clip_5p: 0,
            soft_clip_3p: 0,
            polya_tail_len: get_int_tag(record, b"pt"),
            polya_tail_softclip: None,
            internal_priming: None,
            quant_only: false,
            small_rna: false,
        }))
    }
}

/// BamReader with header for chromosome name resolution.
pub struct BamReaderWithHeader {
    reader: BamReader,
    ref_names: Vec<String>,
}

impl BamReaderWithHeader {
    pub fn new(path: &str) -> IsoForgeResult<Self> {
        let bam_reader = BamReader::new(path)?;
        let reader = bam::Reader::from_path(path)
            .map_err(|e| IsoForgeError::HtsError(format!("Failed to open BAM: {}", e)))?;
        let header_view = reader.header();
        let ref_names: Vec<String> = (0..header_view.target_count())
            .map(|tid| String::from_utf8_lossy(header_view.tid2name(tid)).to_string())
            .collect();
        Ok(Self {
            reader: bam_reader,
            ref_names,
        })
    }

    #[allow(dead_code)]
    pub fn ref_name(&self, tid: u32) -> &str {
        &self.ref_names[tid as usize]
    }

    pub fn read_all_with_chroms(&self, config: &RunConfig) -> IsoForgeResult<Vec<ParsedRead>> {
        let mut reader = bam::Reader::from_path(&self.reader.path)
            .map_err(|e| IsoForgeError::HtsError(format!("Failed to open BAM: {}", e)))?;
        reader
            .set_threads(config.threads.max(1))
            .map_err(|e| IsoForgeError::HtsError(format!("Failed to set threads: {}", e)))?;

        let mut reads = Vec::new();
        let mut record = bam::Record::new();
        while let Some(result) = reader.read(&mut record) {
            result.map_err(|e| IsoForgeError::HtsError(format!("Failed to read record: {}", e)))?;
            if record.is_unmapped() || record.is_secondary() || record.is_supplementary() {
                continue;
            }
            let flag = record.flags();
            if flag != 0 && flag != 16 {
                continue;
            }
            // MAPQ filtering
            if record.mapq() < config.min_mapq {
                continue;
            }
            let tid = record.tid();
            if tid < 0 || (tid as usize) >= self.ref_names.len() {
                continue;
            }

            let read_id = String::from_utf8_lossy(record.qname()).to_string();
            let dataset = if config.use_cb_tag {
                get_string_tag(&record, b"CB").unwrap_or_default()
            } else {
                get_string_tag(&record, b"RG").unwrap_or_default()
            };

            let strand = if record.is_reverse() {
                Strand::Minus
            } else {
                Strand::Plus
            };

            let sam_start = record.pos() + 1;
            let sam_end = record.cigar().end_pos();
            let read_length = record.seq_len() as i64;
            if read_length < config.min_length {
                continue;
            }

            let cigar = record.cigar();
            let Some(excluded) = coverage_exclusion(&record.seq().as_bytes(), &cigar, config) else { continue };
            let coverage = compute_alignment_coverage_excluding(&cigar, excluded);
            if coverage < config.min_coverage {
                continue;
            }

            let md_tag = get_string_tag(&record, b"MD");
            let identity = md_tag
                .as_ref()
                .map(|md| compute_alignment_identity(md, record.seq_len()));
            if let Some(id) = identity {
                if id < config.min_identity {
                    continue;
                }
            }

            let intron_list = get_introns_from_cigar(&cigar, sam_start);
            let splice_sites: Vec<i64> = intron_list
                .iter()
                .enumerate()
                .map(|(i, &x)| if i % 2 == 1 { x + 1 } else { x - 1 })
                .collect();

            let mut positions = Vec::with_capacity(splice_sites.len() + 2);
            positions.push(sam_start);
            positions.extend_from_slice(&splice_sites);
            positions.push(sam_end);

            if strand == Strand::Minus {
                positions.reverse();
            }

            let n_exons = (positions.len() / 2) as i32;
            let chromosome = self.ref_names[tid as usize].clone();
            let cigar_string = cigar_to_string(&cigar);

            reads.push(ParsedRead {
                read_id,
                read_idx: READ_IDX_NONE,
                dataset,
                chromosome,
                strand,
                start: sam_start,
                end: sam_end,
                read_length,
                positions,
                n_exons,
                cigar_string,
                flag,
                fraction_as: get_float_tag(&record, b"fA"),
                custom_label: get_string_tag(&record, b"lC"),
                allelic_label: get_string_tag(&record, b"lA"),
                start_support: get_string_tag(&record, b"tS"),
                end_support: get_string_tag(&record, b"tE"),
                sj_corrected: 0,
                ambiguous_sites: 0,
                sequence_settled: 0,
                cap_signature: false,
                unaligned_3p: false,
                sj_noncanonical: 0,
                sj_canonical: intron_list.len() as u32 / 2,
                microindels_ins: 0,
                microindels_del: 0,
                alignment_coverage: coverage,
                alignment_identity: identity,
                soft_clip_5p: 0,
                soft_clip_3p: 0,
                polya_tail_len: get_int_tag(&record, b"pt"),
                polya_tail_softclip: None,
                internal_priming: None,
                quant_only: false,
                small_rna: false,
            });
        }

        Ok(reads)
    }
}

// ─── Helper functions ────────────────────────────────────────────────────────

fn get_string_tag(record: &bam::Record, tag: &[u8]) -> Option<String> {
    match record.aux(tag) {
        Ok(Aux::String(s)) => Some(s.to_string()),
        _ => None,
    }
}

fn get_float_tag(record: &bam::Record, tag: &[u8]) -> Option<f64> {
    match record.aux(tag) {
        Ok(Aux::Float(f)) => Some(f as f64),
        Ok(Aux::Double(d)) => Some(d),
        _ => None,
    }
}

/// Read an integer tag from a BAM record (supports i8, u8, i16, u16, i32, u32).
fn get_int_tag(record: &bam::Record, tag: &[u8]) -> Option<i64> {
    match record.aux(tag) {
        Ok(Aux::I8(v)) => Some(v as i64),
        Ok(Aux::U8(v)) => Some(v as i64),
        Ok(Aux::I16(v)) => Some(v as i64),
        Ok(Aux::U16(v)) => Some(v as i64),
        Ok(Aux::I32(v)) => Some(v as i64),
        Ok(Aux::U32(v)) => Some(v as i64),
        _ => None,
    }
}

/// Splice-site flags follow transcript order: an intron added at the read's 5' end shifts them by its two sites.
fn shift_ambiguous(mask: &mut u64, strand: Strand, low: bool, high: bool) {
    if (strand == Strand::Plus && low) || (strand == Strand::Minus && high) { *mask <<= 2; }
}

/// Extract soft-clip lengths from the CIGAR string.
///
/// Returns (leading_soft_clip, trailing_soft_clip) in alignment orientation.
pub(crate) fn extract_soft_clips(cigar: &bam::record::CigarStringView) -> (u32, u32) {
    use rust_htslib::bam::record::Cigar::*;
    let ops: Vec<_> = cigar.iter().collect();
    if ops.is_empty() {
        return (0, 0);
    }

    let leading = match ops.first() {
        Some(SoftClip(len)) => *len,
        _ => 0,
    };

    let trailing = if ops.len() > 1 {
        match ops.last() {
            Some(SoftClip(len)) => *len,
            _ => 0,
        }
    } else {
        0
    };

    (leading, trailing)
}

/// Longest scoring run of `target_base` at the start of `clip_seq`.
///
/// +1 per match, -5 per mismatch, abort below -20; a hit needs >=10bp at >=90%.
fn score_polya_run(clip_seq: &[u8], target_base: u8) -> Option<i64> {
    use crate::types::*;

    if clip_seq.len() < POLYA_TAIL_MIN_LEN {
        return None;
    }
    let mut score: i32 = 0;
    let mut max_score: i32 = 0;
    let mut tail_len: usize = 0;
    let mut n_target: usize = 0;

    for (i, &base) in clip_seq.iter().enumerate() {
        if base.to_ascii_uppercase() == target_base {
            score += POLYA_TAIL_SCORE_MATCH;
            n_target += 1;
        } else {
            score += POLYA_TAIL_SCORE_MISMATCH;
        }
        if score > max_score {
            max_score = score;
            tail_len = i + 1;
        }
        if score < POLYA_TAIL_SCORE_ABORT {
            break;
        }
    }

    if tail_len >= POLYA_TAIL_MIN_LEN
        && (n_target as f64 / tail_len as f64) >= POLYA_TAIL_MIN_PERCENT_A
    {
        Some(tail_len as i64)
    } else {
        None
    }
}

/// Transcript strand implied by a polyA tail in the soft clips: plus for a trailing run of A, minus for a leading
/// run of T (SEQ is genome-forward). Returns None unless exactly one end carries a tail.
pub(crate) fn strand_from_polya(
    read_seq: &[u8],
    leading_clip: u32,
    trailing_clip: u32,
) -> Option<Strand> {
    let plus = if trailing_clip > 0 && (trailing_clip as usize) <= read_seq.len() {
        let start = read_seq.len() - trailing_clip as usize;
        score_polya_run(&read_seq[start..], b'A')
    } else {
        None
    };
    let minus = if leading_clip > 0 && (leading_clip as usize) <= read_seq.len() {
        // The leading clip is read outward from the alignment.
        let outward: Vec<u8> = read_seq[..leading_clip as usize].iter().rev().copied().collect();
        score_polya_run(&outward, b'T')
    } else {
        None
    };
    match (plus, minus) {
        (Some(_), None) => Some(Strand::Plus),
        (None, Some(_)) => Some(Strand::Minus),
        _ => None,
    }
}

/// Reads whose strand came from splice motifs, from a poly(A) tail, and how often
/// each overrode the alignment orientation, for the run summary.
pub static STRAND_FROM_MOTIF: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
pub static STRAND_MOTIF_FLIPPED: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
pub static STRAND_INFERRED: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
pub static STRAND_FLIPPED: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

/// Strand of a read: from a poly(A) tail in the soft clips, else from splice motifs when every junction is canonical
/// on one strand, else from the alignment orientation or, in an unoriented library, from the tailed reads and
/// annotation covering it (None when they give none). The alignment strand when `infer_strand_from_polya` is off.
fn resolve_strand(
    qname: &[u8],
    aln_strand: Strand,
    read_seq: &[u8],
    cigar: &bam::record::CigarStringView,
    chrom: &str,
    start: i64,
    end: i64,
    genome: Option<&GenomeRef>,
    index: Option<&crate::strand_evidence::StrandIndex>,
    config: &RunConfig,
) -> Option<Strand> {
    use std::sync::atomic::Ordering::Relaxed;
    if !config.infer_strand_from_polya {
        return Some(aln_strand);
    }
    static TRACE: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    let trace = *TRACE.get_or_init(|| std::env::var_os("ISOFORGE_STRAND_TRACE").is_some());
    let introns = get_introns_from_cigar(cigar, start);
    let motif = genome.and_then(|g| motif_strand(g, chrom, &introns));
    let (leading_clip, trailing_clip) = extract_soft_clips(cigar);
    let tail = strand_from_polya(read_seq, leading_clip, trailing_clip);
    let (resolved, source) = if let Some(s) = tail {
        STRAND_INFERRED.fetch_add(1, Relaxed);
        if s != aln_strand {
            STRAND_FLIPPED.fetch_add(1, Relaxed);
        }
        (Some(s), "tail")
    } else if let Some(s) = motif {
        STRAND_FROM_MOTIF.fetch_add(1, Relaxed);
        if s != aln_strand {
            STRAND_MOTIF_FLIPPED.fetch_add(1, Relaxed);
        }
        (Some(s), "motif")
    } else {
        match index {
            Some(ix) if !ix.oriented => (ix.resolve(chrom, start, end), "locus"),
            _ => (Some(aln_strand), "alignment"),
        }
    };
    if trace {
        let ch = |s: Strand| if s == Strand::Plus { '+' } else { '-' };
        eprintln!(
            "[strand] {}\t{}:{}-{}\tintrons={}\taln={}\tmotif={:?}\ttail={:?}\t{}\t{:?}",
            String::from_utf8_lossy(qname), chrom, start, end, introns.len() / 2,
            ch(aln_strand), motif.map(ch), tail.map(ch), source, resolved.map(ch)
        );
    }
    resolved
}

/// Strand on which every intron of a read carries a canonical splice motif. None when any junction is canonical on
/// neither strand or on both, when the junctions disagree, and for a read without introns.
fn motif_strand(genome: &GenomeRef, chrom: &str, introns: &[i64]) -> Option<Strand> {
    let (mut plus, mut minus, mut n) = (0u32, 0u32, 0u32);
    for pair in introns.chunks_exact(2) {
        n += 1;
        match genome.splice_motif_strand(chrom, pair[0], pair[1]) {
            Ok(Some(Strand::Plus)) => plus += 1,
            Ok(Some(Strand::Minus)) => minus += 1,
            _ => {}
        }
    }
    if n > 0 && plus == n {
        Some(Strand::Plus)
    } else if n > 0 && minus == n {
        Some(Strand::Minus)
    } else {
        None
    }
}

/// Detect a poly(A) tail in the soft-clipped bases. Scoring: +1 per A (+ strand) or T (- strand), -5 per other base,
/// stop when the cumulative score drops below -20; a tail needs ≥10 bp at ≥90% A. Returns the tail length, or None
/// if no tail is found.
fn detect_polya_tail_from_softclip(
    read_seq: &[u8],
    strand: Strand,
    leading_clip: u32,
    trailing_clip: u32,
) -> Option<i64> {
    use crate::types::*;

    // On + strand, polyA tail is at the 3' end (trailing soft-clip)
    // On - strand, polyA tail is at the 5' end (leading soft-clip) as poly-T
    let (clip_len, clip_start, target_base) = match strand {
        Strand::Plus => {
            if trailing_clip == 0 {
                return None;
            }
            let start = read_seq.len().saturating_sub(trailing_clip as usize);
            (trailing_clip, start, b'A')
        }
        Strand::Minus => {
            if leading_clip == 0 {
                return None;
            }
            (leading_clip, 0usize, b'T')
        }
    };

    if (clip_len as usize) < POLYA_TAIL_MIN_LEN {
        return None;
    }

    // Score the soft-clipped region
    let clip_end = clip_start + clip_len as usize;
    let clip_seq = if clip_end <= read_seq.len() {
        &read_seq[clip_start..clip_end]
    } else {
        return None;
    };

    match strand {
        Strand::Plus => score_polya_run(clip_seq, target_base),
        Strand::Minus => {
            // Read the leading clip outward from the alignment, where the tail lies.
            let outward: Vec<u8> = clip_seq.iter().rev().copied().collect();
            score_polya_run(&outward, target_base)
        }
    }
}

/// Unaligned bases before the tail that mean the alignment ends before the transcript's 3' end.
pub(crate) const UNALIGNED_3P_MIN: u32 = 30;

/// Bases of the 3' soft clip, read outward from the alignment, before its poly(A) run starts (a window of
/// `POLYA_TAIL_MIN_LEN` bases at least `POLYA_TAIL_MIN_PERCENT_A` A, T on a minus-strand read); 0 when the clip holds
/// no run.
pub(crate) fn unaligned_before_tail(read_seq: &[u8], strand: Strand, leading_clip: u32, trailing_clip: u32) -> u32 {
    use crate::types::{POLYA_TAIL_MIN_LEN, POLYA_TAIL_MIN_PERCENT_A};
    let (outward, target): (Vec<u8>, u8) = match strand {
        Strand::Plus => (read_seq[read_seq.len().saturating_sub(trailing_clip as usize)..].to_vec(), b'A'),
        Strand::Minus => (read_seq[..(leading_clip as usize).min(read_seq.len())].iter().rev().copied().collect(), b'T'),
    };
    let w = POLYA_TAIL_MIN_LEN;
    if outward.len() < w { return 0; }
    let need = (w as f64 * POLYA_TAIL_MIN_PERCENT_A).ceil() as usize;
    let hit = |b: &u8| b.to_ascii_uppercase() == target;
    let mut n = outward[..w].iter().filter(|b| hit(b)).count();
    for i in 0..=outward.len() - w {
        if i > 0 {
            n -= hit(&outward[i - 1]) as usize;
            n += hit(&outward[i + w - 1]) as usize;
        }
        if n >= need { return i as u32; }
    }
    0
}

/// Extract CIGAR operations as (length, op_code) pairs.
fn extract_cigar_ops(cigar: &bam::record::CigarStringView) -> Vec<(u32, u8)> {
    cigar
        .iter()
        .map(|op| {
            use rust_htslib::bam::record::Cigar::*;
            match op {
                Match(l) => (*l, 0),
                Ins(l) => (*l, 1),
                Del(l) => (*l, 2),
                RefSkip(l) => (*l, 3),
                SoftClip(l) => (*l, 4),
                HardClip(l) => (*l, 5),
                Pad(l) => (*l, 6),
                Equal(l) => (*l, 7),
                Diff(l) => (*l, 8),
            }
        })
        .collect()
}

/// Correct microindels in a CIGAR string: insertions ≤ max_len are removed (the bases are dropped from the read),
/// deletions ≤ max_len are converted to match (reference bases fill in), and adjacent M/=/X runs are then merged
/// into a single M. Returns the operations and the numbers of insertions and deletions absorbed.
fn correct_microindels_in_cigar(
    cigar: &bam::record::CigarStringView,
    max_len: u32,
) -> (Vec<rust_htslib::bam::record::Cigar>, u32, u32) {
    use rust_htslib::bam::record::Cigar::*;

    let mut n_ins: u32 = 0;
    let mut n_del: u32 = 0;

    // First pass: absorb small I/D into adjacent matches
    let mut ops: Vec<rust_htslib::bam::record::Cigar> = Vec::with_capacity(cigar.len());

    for op in cigar.iter() {
        match op {
            // Small insertion: removed.
            Ins(len) if *len <= max_len => {
                n_ins += 1;
            }
            // Small deletion: converted to match over the same reference positions.
            Del(len) if *len <= max_len => {
                ops.push(Match(*len));
                n_del += 1;
            }
            // Everything else passes through unchanged
            other => {
                ops.push(other.clone());
            }
        }
    }

    // Second pass: merge adjacent M/=/X blocks into single M
    let mut merged: Vec<rust_htslib::bam::record::Cigar> = Vec::with_capacity(ops.len());
    for op in &ops {
        match op {
            Match(len) | Equal(len) | Diff(len) => {
                if let Some(Match(prev_len)) = merged.last_mut() {
                    *prev_len += len;
                } else {
                    merged.push(Match(*len));
                }
            }
            other => {
                merged.push(other.clone());
            }
        }
    }

    (merged, n_ins, n_del)
}

/// Introns placed where the read needs no indel against them: an intron takes in an adjacent deletion, or gives up
/// the bases an adjacent insertion repeats (indels of `min_len..=max_len` bases), if a slide that keeps the spliced
/// sequence puts it on GT-AG, else GC-AG, else AT-AC. Updates `introns` (1-based) in place; returns how many moved.
fn place_indel_adjacent_introns(
    record: &bam::Record,
    chrom: &str,
    strand: Strand,
    genome: &GenomeRef,
    introns: &mut [i64],
    min_len: u32,
    max_len: u32,
) -> u32 {
    use rust_htslib::bam::record::Cigar::*;
    if min_len == 0 || introns.is_empty() { return 0; }
    const MAX_SLIDE: i64 = 12;
    let ops: Vec<_> = record.cigar().iter().copied().collect();
    let seq = record.seq().as_bytes();
    // reference and query position at the start of each op
    let (mut r, mut q) = (record.pos() + 1, 0usize);
    let mut at = Vec::with_capacity(ops.len());
    for op in &ops {
        at.push((r, q));
        match op {
            Match(n) | Equal(n) | Diff(n) => { r += *n as i64; q += *n as usize; }
            Del(n) | RefSkip(n) => { r += *n as i64; }
            Ins(n) | SoftClip(n) => { q += *n as usize; }
            HardClip(_) | Pad(_) => {}
        }
    }
    let fetch = |a: i64, b: i64| -> Vec<u8> { genome.fetch_seq(chrom, a - 1, b).unwrap_or_default() }; // 1-based inclusive
    let rank = |u: i64, v: i64| -> Option<u8> {
        let (d, a) = (fetch(u, u + 1), fetch(v - 1, v));
        if d.len() < 2 || a.len() < 2 { return None; }
        match crate::genome_ref::classify_splice_motif_class([d[0], d[1]], [a[0], a[1]], strand) {
            crate::genome_ref::SpliceMotifClass::GtAg => Some(0),
            crate::genome_ref::SpliceMotifClass::GcAg => Some(1),
            crate::genome_ref::SpliceMotifClass::AtAc => Some(2),
            _ => None,
        }
    };
    // best canonical placement of [u, v] within its slide class
    let slide = |u: i64, v: i64| -> Option<(i64, i64)> {
        let left = fetch(u - MAX_SLIDE, u + MAX_SLIDE - 1);   // positions u-12 .. u+11
        let right = fetch(v - MAX_SLIDE + 1, v + MAX_SLIDE);  // positions v-11 .. v+12
        if left.len() as i64 != 2 * MAX_SLIDE || right.len() as i64 != 2 * MAX_SLIDE { return None; }
        let l = |p: i64| left[(p - (u - MAX_SLIDE)) as usize];
        let rr = |p: i64| right[(p - (v - MAX_SLIDE + 1)) as usize];
        let mut best: Option<((u8, i64, i64), i64)> = None;
        for k in -MAX_SLIDE..=MAX_SLIDE {
            let equivalent = if k > 0 { (0..k).all(|i| l(u + i) == rr(v + 1 + i)) }
                else if k < 0 { (k..0).all(|i| l(u + i) == rr(v + 1 + i)) }
                else { true };
            if !equivalent { continue; }
            let Some(m) = rank(u + k, v + k) else { continue };
            let key = (m, k.abs(), -k);
            if best.map_or(true, |(b, _)| key < b) { best = Some((key, k)); }
        }
        best.map(|(_, k)| (u + k, v + k))
    };
    let n = introns.len() / 2;
    let mut moved = 0;
    for (i, op) in ops.iter().enumerate() {
        let RefSkip(len) = op else { continue };
        let (u, v) = (at[i].0, at[i].0 + *len as i64 - 1);
        let Some(k) = (0..n).find(|&k| introns[2 * k] == u && introns[2 * k + 1] == v) else { continue };
        let mut placed = None;
        for j in [i.wrapping_sub(1), i + 1] {
            let Some(adj) = ops.get(j) else { continue };
            let cands: Vec<(i64, i64)> = match adj {
                Del(l) if (min_len..=max_len).contains(l) => {
                    let l = *l as i64;
                    vec![if j + 1 == i { (u - l, v) } else { (u, v + l) }]
                }
                Ins(l) if (min_len..=max_len).contains(l) => {
                    let (l, qs) = (*l as i64, at[j].1);
                    let ins = &seq[qs..qs + l as usize];
                    let mut c = Vec::new();
                    if ins.eq_ignore_ascii_case(&fetch(u, u + l - 1)) { c.push((u + l, v)); }
                    if ins.eq_ignore_ascii_case(&fetch(v - l + 1, v)) { c.push((u, v - l)); }
                    c
                }
                _ => continue,
            };
            placed = cands.into_iter().find_map(|(cu, cv)| slide(cu, cv));
            if placed.is_some() { break; }
        }
        let Some((pu, pv)) = placed else { continue };
        // at least one exonic base must remain on either side, next to the read's ends and neighbouring introns
        let lo = if k == 0 { record.pos() + 1 } else { introns[2 * k - 1] + 1 };
        let hi = if k + 1 == n { record.cigar().end_pos() } else { introns[2 * k + 2] - 1 };
        if pu <= lo || pv >= hi || (pu, pv) == (u, v) { continue; }
        introns[2 * k] = pu;
        introns[2 * k + 1] = pv;
        moved += 1;
    }
    moved
}

/// Extract introns (N operations) from corrected CIGAR operations.
fn get_introns_from_corrected_cigar(
    ops: &[rust_htslib::bam::record::Cigar],
    start: i64,
) -> Vec<i64> {
    use rust_htslib::bam::record::Cigar::*;
    let mut introns = Vec::new();
    let mut pos = start;

    for op in ops {
        match op {
            Match(len) | Equal(len) | Diff(len) | Del(len) => {
                pos += *len as i64;
            }
            RefSkip(len) => {
                let intron_start = pos;
                let intron_end = pos + *len as i64 - 1;
                introns.push(intron_start);
                introns.push(intron_end);
                pos += *len as i64;
            }
            Ins(_) | SoftClip(_) | HardClip(_) | Pad(_) => {}
        }
    }

    introns
}

/// Extract introns from CIGAR operations (N operations).
fn get_introns_from_cigar(cigar: &bam::record::CigarStringView, start: i64) -> Vec<i64> {
    let mut introns = Vec::new();
    let mut pos = start;

    for op in cigar.iter() {
        use rust_htslib::bam::record::Cigar::*;
        match op {
            Match(len) | Equal(len) | Diff(len) | Del(len) => {
                pos += *len as i64;
            }
            RefSkip(len) => {
                let intron_start = pos;
                let intron_end = pos + *len as i64 - 1;
                introns.push(intron_start);
                introns.push(intron_end);
                pos += *len as i64;
            }
            Ins(_) | SoftClip(_) | HardClip(_) | Pad(_) => {}
        }
    }

    introns
}

/// Soft-clipped bases left out of the coverage measure: a poly(A)/poly(T) tail and, at each end, other clipped
/// sequence up to `max_terminal_clip`. None when a non-tail clip is longer than that.
fn coverage_exclusion(seq: &[u8], cigar: &bam::record::CigarStringView, config: &RunConfig) -> Option<u32> {
    let (lead, trail) = extract_soft_clips(cigar);
    // Tail length in each clip, read outward from the alignment.
    let (tail_lead, tail_trail) = if config.coverage_excludes_polya {
        (
            detect_polya_tail_from_softclip(seq, Strand::Minus, lead, trail).unwrap_or(0) as u32,
            detect_polya_tail_from_softclip(seq, Strand::Plus, lead, trail).unwrap_or(0) as u32,
        )
    } else {
        (0, 0)
    };
    if lead - tail_lead > config.max_terminal_clip || trail - tail_trail > config.max_terminal_clip {
        return None;
    }
    Some(lead + trail)
}

/// Coverage with `excluded_clip` soft-clipped bases left out of both the denominator and the unaligned count.
fn compute_alignment_coverage_excluding(cigar: &bam::record::CigarStringView, excluded_clip: u32) -> f64 {
    let mut total_bases: f64 = -(excluded_clip as f64);
    let mut unaligned_bases: f64 = -(excluded_clip as f64);

    for op in cigar.iter() {
        use rust_htslib::bam::record::Cigar::*;
        match op {
            Match(len) | Equal(len) | Diff(len) | Ins(len) | Del(len) => {
                total_bases += *len as f64;
            }
            SoftClip(len) => {
                total_bases += *len as f64;
                unaligned_bases += *len as f64;
            }
            HardClip(len) => {
                total_bases += *len as f64;
                unaligned_bases += *len as f64;
            }
            RefSkip(_) | Pad(_) => {}
        }
    }

    if total_bases <= 0.0 {
        return 0.0;
    }
    (total_bases - unaligned_bases) / total_bases
}

/// Compute alignment identity from MD tag.
fn compute_alignment_identity(md_tag: &str, seq_len: usize) -> f64 {
    let mut matches: f64 = 0.0;
    let mut total_bases = seq_len as f64;
    let mut current_num = String::new();
    let mut in_deletion = false;

    for ch in md_tag.chars() {
        if ch.is_ascii_digit() {
            current_num.push(ch);
            in_deletion = false;
        } else if ch == '^' {
            if !current_num.is_empty() {
                matches += current_num.parse::<f64>().unwrap_or(0.0);
                current_num.clear();
            }
            in_deletion = true;
        } else {
            if !current_num.is_empty() {
                matches += current_num.parse::<f64>().unwrap_or(0.0);
                current_num.clear();
            }
            if in_deletion {
                total_bases += 1.0;
            }
            in_deletion = false;
        }
    }

    if !current_num.is_empty() {
        matches += current_num.parse::<f64>().unwrap_or(0.0);
    }

    if total_bases == 0.0 {
        return 0.0;
    }
    matches / total_bases
}

/// Convert CIGAR view to string representation.
fn cigar_to_string(cigar: &bam::record::CigarStringView) -> String {
    let mut s = String::new();
    for op in cigar.iter() {
        use rust_htslib::bam::record::Cigar::*;
        let (len, ch) = match op {
            Match(l) => (l, 'M'),
            Ins(l) => (l, 'I'),
            Del(l) => (l, 'D'),
            RefSkip(l) => (l, 'N'),
            SoftClip(l) => (l, 'S'),
            HardClip(l) => (l, 'H'),
            Pad(l) => (l, 'P'),
            Equal(l) => (l, '='),
            Diff(l) => (l, 'X'),
        };
        s.push_str(&format!("{}{}", len, ch));
    }
    s
}

#[cfg(test)]
mod unaligned_tail_tests {
    use super::*;

    #[test]
    fn transcript_sequence_before_the_tail_counts_as_unaligned() {
        // + strand: trailing clip read outward
        let tail_at_once = b"ACGTACGTACGTAAAAAAAAAAAAAAAAAAAA";
        assert_eq!(unaligned_before_tail(tail_at_once, Strand::Plus, 0, 20), 0);
        let exon_then_tail = [&b"ACGT"[..], &[b'C'; 40][..], &[b'A'; 25][..]].concat();
        // the first window of 10 bases with 9 A's starts after 39 bases
        assert_eq!(unaligned_before_tail(&exon_then_tail, Strand::Plus, 0, 65), 39);
        // - strand: leading clip, the tail as T's farthest from the alignment
        let minus = [&[b'T'; 25][..], &[b'G'; 35][..], &b"ACGTACGT"[..]].concat();
        assert_eq!(unaligned_before_tail(&minus, Strand::Minus, 60, 0), 34);
        // a clip without a tail gives 0
        assert_eq!(unaligned_before_tail(&[b'G'; 50], Strand::Plus, 0, 50), 0);
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn tails_are_read_outward_from_the_alignment() {
        // Genome-forward SEQ of a minus-strand read: adapter, then poly(T) against the alignment.
        let adapter = b"ACTTGCCTGTCGCTCTATCTTCGGAAGCAGTGGTATCAACGCAGAGTAC";
        let mut seq = adapter.to_vec();
        seq.extend_from_slice(&[b'T'; 30]);
        let lead = seq.len() as u32;
        seq.extend_from_slice(b"GATTACAGATTACAGATTACAGATTACA");
        assert_eq!(strand_from_polya(&seq, lead, 0), Some(Strand::Minus));
        assert_eq!(detect_polya_tail_from_softclip(&seq, Strand::Minus, lead, 0), Some(30));
        // Plus-strand read: poly(A) against the alignment, then adapter.
        let mut seq = b"GATTACAGATTACAGATTACAGATTACA".to_vec();
        let aligned = seq.len();
        seq.extend_from_slice(&[b'A'; 30]);
        seq.extend_from_slice(adapter);
        let trail = (seq.len() - aligned) as u32;
        assert_eq!(strand_from_polya(&seq, 0, trail), Some(Strand::Plus));
        assert!(detect_polya_tail_from_softclip(&seq, Strand::Plus, 0, trail) >= Some(30));
    }

    use super::*;

    #[test]
    fn test_microindels_preserve_intron_coordinates() {
        use bam::record::{Cigar::*, CigarString};

        let cigar = CigarString(vec![
            SoftClip(12), Match(40), Ins(3), Match(60), Del(2),
            RefSkip(100), Match(30), Del(5), Match(40), Ins(8),
            RefSkip(80), Match(50), HardClip(10),
        ]).into_view(100);
        let expected = get_introns_from_cigar(&cigar, 101);
        assert_eq!(expected, vec![203, 302, 378, 457]);
        for max_len in [0, 1, 3, 5, 10] {
            let (corrected, _, _) = correct_microindels_in_cigar(&cigar, max_len);
            assert_eq!(get_introns_from_corrected_cigar(&corrected, 101), expected);
        }
    }

    #[test]
    fn test_noncanonical_reads_preserve_ends_and_polya_in_both_readers() {
        use bam::record::{Cigar::*, CigarString};

        let directory = tempfile::tempdir().unwrap();
        let fasta = directory.path().join("genome.fa");
        std::fs::write(&fasta, format!(">chr1\n{}\n", "C".repeat(1000))).unwrap();
        std::fs::write(directory.path().join("genome.fa.fai"),
            "chr1\t1000\t6\t1000\t1001\n").unwrap();
        let reader = IntegratedBamReader {
            path: String::new(),
            ref_names: vec!["chr1".to_string()],
            genome: Some(GenomeRef::open(&fasta).unwrap()),
            pipeline: ReaderPipeline::default(),
            chimeric: std::sync::Mutex::new(Vec::new()),
            collect_chimeric: false,
            gene_spans: None,
        };
        let config = RunConfig { min_length: 1, ..RunConfig::default() };

        for strand in [Strand::Plus, Strand::Minus] {
            let (ops, sequence, expected) = match strand {
                Strand::Plus => (
                    vec![Match(100), RefSkip(100), Match(100), SoftClip(20)],
                    format!("{}{}", "C".repeat(200), "A".repeat(20)),
                    vec![101, 200, 301, 400],
                ),
                Strand::Minus => (
                    vec![SoftClip(20), Match(100), RefSkip(100), Match(100)],
                    format!("{}{}", "T".repeat(20), "C".repeat(200)),
                    vec![400, 301, 200, 101],
                ),
            };
            let mut record = bam::Record::new();
            record.set(b"noncanonical", Some(&CigarString(ops)), sequence.as_bytes(), &[30; 220]);
            record.set_tid(0);
            record.set_pos(100);
            record.set_mapq(60);
            record.set_flags(0);
            let serial = reader.process_record(&record, &config, "sample").unwrap().unwrap();
            let parallel = process_record_with_refs(
                &record, &reader.ref_names, reader.genome.as_ref(),
                &reader.pipeline, &config, "sample",
            ).unwrap().unwrap();
            for parsed in [serial, parallel] {
                assert_eq!(parsed.positions, expected);
                assert_eq!(parsed.strand, strand);
                assert_eq!((parsed.start, parsed.end), (101, 400));
                assert_eq!(parsed.sj_noncanonical, 1);
                assert_eq!(parsed.sj_corrected, 0);
                assert_eq!(parsed.sj_canonical, 0);
                assert_eq!(parsed.polya_tail_softclip, Some(20));
            }
        }
    }

    #[test]
    fn a_read_of_an_annotated_small_rna_passes_the_length_and_mapping_quality_gates() {
        use bam::record::{Cigar::*, CigarString};

        let directory = tempfile::tempdir().unwrap();
        let fasta = directory.path().join("genome.fa");
        std::fs::write(&fasta, format!(">chr1\n{}\n", "C".repeat(1000))).unwrap();
        std::fs::write(directory.path().join("genome.fa.fai"), "chr1\t1000\t6\t1000\t1001\n").unwrap();
        let genome = GenomeRef::open(&fasta).unwrap();
        let names = vec!["chr1".to_string()];
        // a 153-nt small RNA on the plus strand at 300..=452
        let index = crate::small_rna::SmallRnaIndex::from_records([("chr1", 300, 452, Strand::Plus, 7, false)]);
        let with_rule = ReaderPipeline { small_rna: Some(Arc::new(index)), ..ReaderPipeline::default() };
        let without = ReaderPipeline::default();
        let config = RunConfig::default();
        assert!(config.min_length == 200 && config.min_mapq == 16);

        let record = |pos: i64, ops: Vec<bam::record::Cigar>, sequence: String, mapq: u8, flags: u16| {
            let mut r = bam::Record::new();
            r.set(b"read", Some(&CigarString(ops)), sequence.as_bytes(), &vec![30; sequence.len()]);
            r.set_tid(0);
            r.set_pos(pos);
            r.set_mapq(mapq);
            r.set_flags(flags);
            r
        };
        let parse = |r: &bam::Record, p: &ReaderPipeline| process_record_with_refs(r, &names, Some(&genome), p, &config, "sample").unwrap();

        // 160 bases at 300..=459 with mapping quality 0: a mature small-RNA read of a multi-copy gene
        let mature = record(299, vec![Match(160)], "C".repeat(160), 0, 0);
        assert!(parse(&mature, &without).is_none());
        let kept = parse(&mature, &with_rule).unwrap();
        assert!(kept.small_rna && !kept.quant_only);
        assert_eq!((kept.start, kept.end, kept.read_length), (300, 459, 160));

        // the same alignment on the other strand, and one elsewhere, stay out
        assert!(parse(&record(299, vec![Match(160)], "C".repeat(160), 0, 16), &with_rule).is_none());
        assert!(parse(&record(600, vec![Match(160)], "C".repeat(160), 0, 0), &with_rule).is_none());

        // a longer read over the small RNA passes the usual gates as any read and is no small-RNA read
        let precursor = parse(&record(99, vec![Match(360)], "C".repeat(360), 60, 0), &with_rule).unwrap();
        assert!(!precursor.small_rna);

        // a 3'-extended, polyadenylated molecule of the small RNA: 60 bases past its end, then the tail
        let tailed = record(299, vec![Match(213), SoftClip(40)], format!("{}{}", "C".repeat(213), "A".repeat(40)), 60, 0);
        let kept = parse(&tailed, &with_rule).unwrap();
        assert!(kept.small_rna);
        assert_eq!(kept.polya_tail_softclip, Some(40));
        assert!(!parse(&tailed, &without).unwrap().small_rna);
    }

    #[test]
    fn an_indel_against_an_intron_is_placed_where_the_read_needs_none() {
        use bam::record::{Cigar::*, CigarString};
        // A tandem acceptor (NAGNAG): exon 1 ends ...GGACAG at 100, the intron runs GTAGGG..CGCAGCAG over 101-206,
        // exon 2 starts CTCT at 207. The read is aligned to the upstream AG (203) with the CAG at 98-100 deleted.
        let exon1 = format!("{}GGACAG", "C".repeat(94));
        let intron = format!("GTAGGG{}CCCTCGCAGCAG", "T".repeat(88));
        let exon2 = format!("CTCTACCAG{}", "C".repeat(141));
        let chrom = format!("{}{}{}", exon1, intron, exon2);
        let directory = tempfile::tempdir().unwrap();
        let fasta = directory.path().join("genome.fa");
        std::fs::write(&fasta, format!(">chr1\n{}\n", chrom)).unwrap();
        std::fs::write(directory.path().join("genome.fa.fai"), format!("chr1\t{}\t6\t{}\t{}\n", chrom.len(), chrom.len(), chrom.len() + 1)).unwrap();
        let genome = GenomeRef::open(&fasta).unwrap();
        let read = format!("{}{}", &chrom[..97], &chrom[203..306]);
        let parse = |ops: Vec<bam::record::Cigar>, seq: &str, min_len: u32| {
            let mut record = bam::Record::new();
            record.set(b"r", Some(&CigarString(ops)), seq.as_bytes(), &vec![30; seq.len()]);
            record.set_tid(0);
            record.set_pos(0);
            let mut introns = get_introns_from_cigar(&record.cigar(), 1);
            let moved = place_indel_adjacent_introns(&record, "chr1", Strand::Plus, &genome, &mut introns, min_len, 12);
            (introns, moved)
        };
        // the deletion joins the intron, which slides three bases to GT-AG: exon 2 starts at 207
        assert_eq!(parse(vec![Match(97), Del(3), RefSkip(103), Match(103)], &read, 3), (vec![101, 206], 1));
        // below the minimum indel length the aligner's intron stands
        assert_eq!(parse(vec![Match(97), Del(3), RefSkip(103), Match(103)], &read, 4), (vec![101, 203], 0));
        // an alignment that already needs no indel is left alone
        let exact = format!("{}{}", &chrom[..100], &chrom[206..306]);
        assert_eq!(parse(vec![Match(100), RefSkip(106), Match(100)], &exact, 3), (vec![101, 206], 0));
    }

    #[test]
    fn test_alignment_identity() {
        let identity = compute_alignment_identity("100", 100);
        assert!((identity - 1.0).abs() < 0.001);

        let identity = compute_alignment_identity("50A50", 101);
        assert!((identity - 100.0 / 101.0).abs() < 0.001);
    }

    #[test]
    fn test_pipeline_config_defaults() {
        let pipeline = ReaderPipeline::default();
        assert_eq!(pipeline.frac_a_range, 20);
        assert_eq!(pipeline.max_sj_offset, 5);
        assert!(pipeline.correct_sj);
    }
}
