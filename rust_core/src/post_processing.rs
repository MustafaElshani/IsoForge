// Post-processing: reads the annotation and result Parquet files and writes
// the summary, transcript whitelist, abundance matrix, GTF, per-read
// annotations and pruned cluster files.

use std::collections::{HashMap, HashSet, BTreeMap};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write, BufWriter};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use arrow::array::*;
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};

use crate::errors::*;
use crate::parquet_io;

// ═══════════════════════════════════════════════════════════════════════════════
// DATA LOADING
// ═══════════════════════════════════════════════════════════════════════════════

/// Column-oriented records loaded from Parquet.

/// A single gene record (from genes.parquet + novel_genes.parquet)
#[derive(Debug, Clone)]
pub struct GeneRecord {
    pub gene_id: i64,
    pub gene_name: Option<String>,
    pub chromosome: String,
    pub start: i64,
    pub end: i64,
    pub strand: String,
}

/// A single transcript record (from transcripts.parquet + novel_transcripts.parquet)
#[derive(Debug, Clone)]
pub struct TranscriptRecord {
    pub transcript_id: i64,
    pub gene_id: i64,
    pub transcript_name: Option<String>,
    pub start_exon: i64,
    pub jn_path: Option<String>,
    pub end_exon: i64,
    pub start_vertex: i64,
    pub end_vertex: i64,
    pub n_exons: i32,
    pub chromosome: String,
    pub start_pos: i64,
    pub end_pos: i64,
}

/// A single edge record
#[derive(Debug, Clone)]
pub struct EdgeRecord {
    pub edge_id: i64,
    pub v1: i64,
    pub v2: i64,
    pub edge_type: String,
    pub strand: String,
}

/// A single location record
#[derive(Debug, Clone)]
pub struct LocationRecord {
    pub location_id: i64,
    pub chromosome: String,
    pub position: i64,
}

/// Observed read — one row of observed.parquet. Repeated strings are shared (interned) and read
/// names are not held: the read-annotation writer reads them back from the file in row order.
#[derive(Debug, Clone, Default)]
pub struct ObservedRead {
    pub dataset: Arc<str>,
    pub gene_id: i64,
    pub transcript_id: i64,
    pub chromosome: Arc<str>,
    pub start_pos: i64,
    pub end_pos: i64,
    pub strand: Arc<str>,
    pub read_length: i64,
    pub n_exons: i32,
    pub start_vertex: i64,
    pub end_vertex: i64,
    pub start_exon: i64,
    pub end_exon: i64,
    pub start_delta: Option<i64>,
    pub end_delta: Option<i64>,
    pub fraction_as: Option<f64>,
    pub custom_label: Option<Arc<str>>,
    pub allelic_label: Option<Arc<str>>,
    pub start_support: Option<Arc<str>>,
    pub end_support: Option<Arc<str>>,
    // PAS fields
    pub has_pas_signal: Option<bool>,
    pub pas_motif: Option<Arc<str>>,
    pub pas_distance: Option<i64>,
    pub is_internal_priming: Option<bool>,
    // Cluster fields
    pub tss_site: Option<i64>,
    pub tss_feature: Option<Arc<str>>,
    pub polya_site: Option<i64>,
    pub polya_feature: Option<Arc<str>>,
    pub polya_signal: Option<Arc<str>>,
    pub cluster_supported: Option<bool>,
    /// PolyA tail length from the basecaller `pt:i` tag (None if the tag is absent)
    pub polya_length: Option<i64>,
    /// Largest vertex-snap distance of an internal splice site; None when the column is absent
    pub max_splice_shift: Option<i32>,
    /// Untemplated G's at the 5' end (a capped 5' end); None when the column is absent
    pub cap_signature: Option<bool>,
    /// Poly(A) tail found in the soft clip, in bases; None when none was found or the
    /// column is absent
    pub polya_softclip: Option<i64>,
    /// A fragment of the model it is credited to (`mark_fragments`): bit 0 when its 5' part is missing, bit 1 when
    /// its 3' part is; 0 for a read that holds every exon of its model.
    pub fragment: u8,
}

impl ObservedRead {
    /// End-validated and holding every exon of its model (no fragment of it).
    pub fn validated(&self) -> bool { self.cluster_supported == Some(true) && self.fragment == 0 }
}

const MISSING_5P: u8 = 1;
const MISSING_3P: u8 = 2;

/// Marks each read credited to a model that holds more exons than the read does (a fragment of it) with the part it
/// lacks: the 5' part when it starts past the model's first exon, the 3' part when it ends before its last, both
/// when neither holds. Returns the number of reads marked.
fn mark_fragments(observed: &mut [ObservedRead], transcripts: &[TranscriptRecord], tx_map: &FxHashMap<i64, usize>,
                  chains: &FxHashMap<i64, Vec<i64>>) -> usize {
    let mut n = 0;
    for r in observed.iter_mut() {
        r.fragment = 0;
        let Some(&i) = tx_map.get(&r.transcript_id) else { continue };
        if r.n_exons >= transcripts[i].n_exons { continue; }
        let mut f = 0;
        if let Some(c) = chains.get(&r.transcript_id).filter(|c| c.len() >= 4) {
            let (donor, acceptor) = (c[1], c[c.len() - 2]);
            let plus = &*r.strand != "-";
            let (lo, hi) = (r.start_pos.min(r.end_pos), r.start_pos.max(r.end_pos));
            if if plus { lo > donor } else { hi < donor } { f |= MISSING_5P; }
            if if plus { hi < acceptor } else { lo > acceptor } { f |= MISSING_3P; }
        }
        r.fragment = if f == 0 { MISSING_5P | MISSING_3P } else { f };
        n += 1;
    }
    n
}

/// One shared copy of each distinct string.
#[derive(Default)]
struct Interner(FxHashMap<Box<str>, Arc<str>>);

impl Interner {
    fn get(&mut self, s: &str) -> Arc<str> {
        if let Some(a) = self.0.get(s) {
            return a.clone();
        }
        let a: Arc<str> = Arc::from(s);
        self.0.insert(s.into(), a.clone());
        a
    }
    fn opt(&mut self, col: &StringArray, i: usize) -> Option<Arc<str>> {
        (!col.is_null(i)).then(|| self.get(col.value(i)))
    }
}

/// The read names of observed.parquet in row order, taken a run of rows at a time.
struct ReadNames {
    reader: parquet::arrow::arrow_reader::ParquetRecordBatchReader,
    batch: Option<StringArray>,
    pos: usize,
}

impl ReadNames {
    fn open(path: &Path) -> IsoForgeResult<Self> {
        let builder = ParquetRecordBatchReaderBuilder::try_new(File::open(path)?)?;
        let mask = parquet::arrow::ProjectionMask::roots(builder.parquet_schema(), [0]);
        Ok(ReadNames { reader: builder.with_projection(mask).build()?, batch: None, pos: 0 })
    }

    fn take(&mut self, n: usize) -> IsoForgeResult<Vec<String>> {
        let mut out = Vec::with_capacity(n);
        while out.len() < n {
            if self.batch.as_ref().is_none_or(|b| self.pos >= b.len()) {
                let next = self.reader.next().ok_or_else(|| {
                    IsoForgeError::IntegrityError("observed.parquet holds fewer reads than were loaded".into())
                })??;
                let names = next.column(0).as_any().downcast_ref::<StringArray>().ok_or_else(|| {
                    IsoForgeError::IntegrityError("observed.parquet read_name is not a string column".into())
                })?;
                self.batch = Some(names.clone());
                self.pos = 0;
            }
            let b = self.batch.as_ref().unwrap();
            let k = (n - out.len()).min(b.len() - self.pos);
            out.extend((self.pos..self.pos + k).map(|i| b.value(i).to_string()));
            self.pos += k;
        }
        Ok(out)
    }
}

/// Abundance row
#[derive(Debug, Clone)]
pub struct AbundanceRecord {
    pub gene_id: i64,
    pub transcript_id: i64,
    pub dataset: String,
    pub count: i64,
}

/// Counters from counters.json
#[derive(Debug, Clone)]
pub struct Counters {
    pub genes: i64,
    pub transcripts: i64,
    pub id_prefix: String,
    data: HashMap<String, serde_json::Value>,
}

/// Annotation lookup tables holding the distinct (id → value) pairs used by
/// post-processing.
pub struct AnnotationDicts {
    // ── Gene annotation lookups ──────────────────────────────────────────
    /// gene_id attribute → formatted ID  (annot_name-filtered)
    pub gene_id_map: FxHashMap<i64, String>,
    /// gene_name attribute  (annot_name-filtered)
    pub gene_name_map: FxHashMap<i64, String>,
    /// gene_type attribute  (annot_name-filtered)
    pub gene_type_map: FxHashMap<i64, String>,
    /// gene_status attribute (all sources)
    pub gene_status_map: FxHashMap<i64, String>,
    /// Derived gene novelty: gene_id → novelty label  (Known / Antisense / …)
    pub gene_novelty: FxHashMap<i64, String>,

    // ── Transcript annotation lookups ────────────────────────────────────
    /// transcript_id attribute → formatted ID  (annot_name-filtered)
    pub tx_id_map: FxHashMap<i64, String>,
    /// transcript_name attribute  (annot_name-filtered)
    pub tx_name_map: FxHashMap<i64, String>,
    /// transcript_type attribute  (annot_name-filtered)
    pub tx_type_map: FxHashMap<i64, String>,
    /// transcript_status attribute (all sources)
    pub tx_status_map: FxHashMap<i64, String>,
    /// Derived transcript novelty: tx_id → novelty label
    pub tx_novelty: FxHashMap<i64, String>,

    // ── Special sets ─────────────────────────────────────────────────────
    pub ism_prefix: FxHashSet<i64>,
    pub ism_suffix: FxHashSet<i64>,
    pub known_tx_ids: FxHashSet<i64>,
    /// Annotated records split off by `split_annotated_forms`, written at their annotated ends.
    pub annotated_forms: FxHashSet<i64>,
    /// Annotated records that took the ID back from their passing end variant
    /// (`reclaim_annotated_records`), written where the variant was: at its peaks.
    pub reclaimed: FxHashSet<i64>,
    /// End variants of annotated transcripts: variant -> (annotated parent, whether the
    /// variant's ends are compatible with the parent's annotated ends), from the run's
    /// end_variants.parquet and from `split_annotated_forms`.
    pub variant_parent: FxHashMap<i64, (i64, bool)>,
    /// Annotated small RNAs counted, not assembled (small_rna.rs): no end variant is split from them, and only the
    /// reads they hold count for them.
    pub small_rna: FxHashMap<i64, crate::small_rna::SmallRna>,
    /// The read-derived 5' and 3' ends of a variant made by `split_annotated_forms`, where they differ from the
    /// annotated ones. Its chain stays the annotated record's; `PostDB::chain_ends` takes its ends from here.
    pub variant_ends: FxHashMap<i64, (Option<i64>, Option<i64>)>,
}

impl AnnotationDicts {
    fn new() -> Self {
        AnnotationDicts {
            gene_id_map: FxHashMap::default(),
            gene_name_map: FxHashMap::default(),
            gene_type_map: FxHashMap::default(),
            gene_status_map: FxHashMap::default(),
            gene_novelty: FxHashMap::default(),
            tx_id_map: FxHashMap::default(),
            tx_name_map: FxHashMap::default(),
            tx_type_map: FxHashMap::default(),
            tx_status_map: FxHashMap::default(),
            tx_novelty: FxHashMap::default(),
            ism_prefix: FxHashSet::default(),
            ism_suffix: FxHashSet::default(),
            known_tx_ids: FxHashSet::default(),
            annotated_forms: FxHashSet::default(),
            reclaimed: FxHashSet::default(),
            variant_parent: FxHashMap::default(),
            small_rna: FxHashMap::default(),
            variant_ends: FxHashMap::default(),
        }
    }
}

/// Discover the annotation source name from an annotation Parquet file.
///
/// Returns the first `annot_name` value that is not in the exclude set.
fn discover_annot_name(path: &Path, exclude: &HashSet<&str>) -> Option<String> {
    if !path.exists() {
        return None;
    }
    let file = match File::open(path) {
        Ok(f) => f,
        Err(_) => return None,
    };
    let builder = match ParquetRecordBatchReaderBuilder::try_new(file) {
        Ok(b) => b,
        Err(_) => return None,
    };
    let reader = match builder.build() {
        Ok(r) => r,
        Err(_) => return None,
    };
    for batch in reader {
        let batch = match batch {
            Ok(b) => b,
            Err(_) => continue,
        };
        let annot_names = batch.column(1).as_any().downcast_ref::<StringArray>()?;
        for i in 0..batch.num_rows() {
            let name = annot_names.value(i);
            if !exclude.contains(name) {
                return Some(name.to_string());
            }
        }
    }
    None
}

/// Stream through an annotation Parquet file and populate `dicts`.
///
/// Parquet schema: (id: i64, annot_name: Utf8, source: Utf8, attribute: Utf8, value: Utf8)
///
/// `is_gene` selects whether to populate gene-level dicts or transcript-level
/// dicts.  `annot_name_filter` is the annotation source used for name/id/type
/// lookups (gene_id, gene_name, etc.).
fn stream_annot_into_dicts(
    path: &Path,
    dicts: &mut AnnotationDicts,
    annot_name_filter: &str,
    is_gene: bool,
) -> IsoForgeResult<()> {
    if !path.exists() {
        return Ok(());
    }
    let file = File::open(path)?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file)?;
    let reader = builder.build()?;

    // Novelty maps
    let gene_nov: FxHashMap<&str, &str> = GENE_NOVELTY_MAP.iter().copied().collect();
    let tx_nov: FxHashMap<&str, &str> = TX_NOVELTY_MAP.iter().copied().collect();

    for batch_result in reader {
        let batch = batch_result?;
        let ids = batch.column(0).as_any().downcast_ref::<Int64Array>().unwrap();
        let annot_names = batch.column(1).as_any().downcast_ref::<StringArray>().unwrap();
        // column 2 = source (unused)
        let attributes = batch.column(3).as_any().downcast_ref::<StringArray>().unwrap();
        let values = batch.column(4).as_any().downcast_ref::<StringArray>().unwrap();

        for i in 0..batch.num_rows() {
            let id = ids.value(i);
            let an = annot_names.value(i);
            let attr = attributes.value(i);
            let val = values.value(i);

            if is_gene {
                // Annot-name-filtered gene dicts
                if an == annot_name_filter {
                    match attr {
                        "gene_id" => { dicts.gene_id_map.entry(id).or_insert_with(|| val.to_string()); }
                        "gene_name" => { dicts.gene_name_map.entry(id).or_insert_with(|| val.to_string()); }
                        "gene_type" => { dicts.gene_type_map.entry(id).or_insert_with(|| val.to_string()); }
                        _ => {}
                    }
                }
                // Unfiltered gene dicts
                if attr == "gene_status" {
                    dicts.gene_status_map.entry(id).or_insert_with(|| val.to_string());
                }
                // Gene novelty
                if let Some(&label) = gene_nov.get(attr) {
                    let match_val = if attr.ends_with("_status") { "KNOWN" } else { "TRUE" };
                    if val == match_val {
                        dicts.gene_novelty.entry(id).or_insert_with(|| label.to_string());
                    }
                }
            } else {
                // Annot-name-filtered transcript dicts
                if an == annot_name_filter {
                    match attr {
                        "transcript_id" => { dicts.tx_id_map.entry(id).or_insert_with(|| val.to_string()); }
                        "transcript_name" => { dicts.tx_name_map.entry(id).or_insert_with(|| val.to_string()); }
                        "transcript_type" => { dicts.tx_type_map.entry(id).or_insert_with(|| val.to_string()); }
                        _ => {}
                    }
                }
                // Unfiltered transcript dicts
                if attr == "transcript_status" {
                    dicts.tx_status_map.entry(id).or_insert_with(|| val.to_string());
                    if val == "KNOWN" {
                        dicts.known_tx_ids.insert(id);
                    }
                }
                // ISM sets
                if attr == "ISM-prefix_transcript" && val == "TRUE" {
                    dicts.ism_prefix.insert(id);
                }
                if attr == "ISM-suffix_transcript" && val == "TRUE" {
                    dicts.ism_suffix.insert(id);
                }
                // Transcript novelty
                if let Some(&label) = tx_nov.get(attr) {
                    let match_val = if attr.ends_with("_status") { "KNOWN" } else { "TRUE" };
                    if val == match_val {
                        dicts.tx_novelty.entry(id).or_insert_with(|| label.to_string());
                    }
                }
            }
        }
    }
    Ok(())
}

/// Stream run-time annotations and split them into gene / transcript dicts.
///
/// Gene-relevant attributes go into the gene dicts; the rest go into
/// transcript dicts (excluding `exon_status`).
fn stream_run_annotations_into_dicts(
    path: &Path,
    dicts: &mut AnnotationDicts,
    annot_name_filter: &str,
) -> IsoForgeResult<()> {
    if !path.exists() {
        return Ok(());
    }
    let file = File::open(path)?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file)?;
    let reader = builder.build()?;

    let gene_attrs: HashSet<&str> = [
        "gene_status", "gene_id", "gene_name", "gene_type",
        "antisense_gene", "intergenic_novel", "fusion_novel",
        "gene_antisense_to_IDs",
    ].iter().copied().collect();
    let exon_attrs: HashSet<&str> = ["exon_status"].iter().copied().collect();

    let gene_nov: FxHashMap<&str, &str> = GENE_NOVELTY_MAP.iter().copied().collect();
    let tx_nov: FxHashMap<&str, &str> = TX_NOVELTY_MAP.iter().copied().collect();

    for batch_result in reader {
        let batch = batch_result?;
        let ids = batch.column(0).as_any().downcast_ref::<Int64Array>().unwrap();
        let annot_names = batch.column(1).as_any().downcast_ref::<StringArray>().unwrap();
        let attributes = batch.column(3).as_any().downcast_ref::<StringArray>().unwrap();
        let values = batch.column(4).as_any().downcast_ref::<StringArray>().unwrap();

        for i in 0..batch.num_rows() {
            let id = ids.value(i);
            let an = annot_names.value(i);
            let attr = attributes.value(i);
            let val = values.value(i);

            if gene_attrs.contains(attr) {
                // Gene-relevant attributes
                if an == annot_name_filter {
                    match attr {
                        "gene_id" => { dicts.gene_id_map.entry(id).or_insert_with(|| val.to_string()); }
                        "gene_name" => { dicts.gene_name_map.entry(id).or_insert_with(|| val.to_string()); }
                        "gene_type" => { dicts.gene_type_map.entry(id).or_insert_with(|| val.to_string()); }
                        _ => {}
                    }
                }
                if attr == "gene_status" {
                    dicts.gene_status_map.entry(id).or_insert_with(|| val.to_string());
                }
                if let Some(&label) = gene_nov.get(attr) {
                    let match_val = if attr.ends_with("_status") { "KNOWN" } else { "TRUE" };
                    if val == match_val {
                        dicts.gene_novelty.entry(id).or_insert_with(|| label.to_string());
                    }
                }
            } else if !exon_attrs.contains(attr) {
                // Transcript-relevant attributes
                if an == annot_name_filter {
                    match attr {
                        "transcript_id" => { dicts.tx_id_map.entry(id).or_insert_with(|| val.to_string()); }
                        "transcript_name" => { dicts.tx_name_map.entry(id).or_insert_with(|| val.to_string()); }
                        "transcript_type" => { dicts.tx_type_map.entry(id).or_insert_with(|| val.to_string()); }
                        _ => {}
                    }
                }
                if attr == "transcript_status" {
                    dicts.tx_status_map.entry(id).or_insert_with(|| val.to_string());
                    if val == "KNOWN" {
                        dicts.known_tx_ids.insert(id);
                    }
                }
                if attr == "ISM-prefix_transcript" && val == "TRUE" {
                    dicts.ism_prefix.insert(id);
                }
                if attr == "ISM-suffix_transcript" && val == "TRUE" {
                    dicts.ism_suffix.insert(id);
                }
                if let Some(&label) = tx_nov.get(attr) {
                    let match_val = if attr.ends_with("_status") { "KNOWN" } else { "TRUE" };
                    if val == match_val {
                        dicts.tx_novelty.entry(id).or_insert_with(|| label.to_string());
                    }
                }
            }
        }
    }
    Ok(())
}

// ─── Parquet Loaders ────────────────────────────────────────────────────────

fn load_genes_from_parquet(path: &Path) -> IsoForgeResult<Vec<GeneRecord>> {
    if !path.exists() {
        return Ok(Vec::new());
    }
    let batches = parquet_io::read_parquet_batches(path)?;
    let mut result = Vec::new();
    for batch in &batches {
        let ids = batch.column(0).as_any().downcast_ref::<Int64Array>().unwrap();
        let names = batch.column(1).as_any().downcast_ref::<StringArray>().unwrap();
        let chroms = batch.column(2).as_any().downcast_ref::<StringArray>().unwrap();
        let starts = batch.column(3).as_any().downcast_ref::<Int64Array>().unwrap();
        let ends = batch.column(4).as_any().downcast_ref::<Int64Array>().unwrap();
        let strands = batch.column(5).as_any().downcast_ref::<StringArray>().unwrap();
        for i in 0..batch.num_rows() {
            result.push(GeneRecord {
                gene_id: ids.value(i),
                gene_name: if names.is_null(i) { None } else { Some(names.value(i).to_string()) },
                chromosome: chroms.value(i).to_string(),
                start: starts.value(i),
                end: ends.value(i),
                strand: strands.value(i).to_string(),
            });
        }
    }
    Ok(result)
}

fn load_transcripts_from_parquet(path: &Path) -> IsoForgeResult<Vec<TranscriptRecord>> {
    if !path.exists() {
        return Ok(Vec::new());
    }
    let batches = parquet_io::read_parquet_batches(path)?;
    let mut result = Vec::new();
    for batch in &batches {
        let t_ids = batch.column(0).as_any().downcast_ref::<Int64Array>().unwrap();
        let g_ids = batch.column(1).as_any().downcast_ref::<Int64Array>().unwrap();
        let t_names = batch.column(2).as_any().downcast_ref::<StringArray>().unwrap();
        let start_exons = batch.column(3).as_any().downcast_ref::<Int64Array>().unwrap();
        let jn_paths = batch.column(4).as_any().downcast_ref::<StringArray>().unwrap();
        let end_exons = batch.column(5).as_any().downcast_ref::<Int64Array>().unwrap();
        let start_verts = batch.column(6).as_any().downcast_ref::<Int64Array>().unwrap();
        let end_verts = batch.column(7).as_any().downcast_ref::<Int64Array>().unwrap();
        let n_exons_col = batch.column(8).as_any().downcast_ref::<Int32Array>().unwrap();
        let chroms = batch.column(9).as_any().downcast_ref::<StringArray>().unwrap();
        let starts = batch.column(10).as_any().downcast_ref::<Int64Array>().unwrap();
        let ends = batch.column(11).as_any().downcast_ref::<Int64Array>().unwrap();
        for i in 0..batch.num_rows() {
            result.push(TranscriptRecord {
                transcript_id: t_ids.value(i),
                gene_id: g_ids.value(i),
                transcript_name: if t_names.is_null(i) { None } else { Some(t_names.value(i).to_string()) },
                start_exon: start_exons.value(i),
                jn_path: if jn_paths.is_null(i) { None } else { Some(jn_paths.value(i).to_string()) },
                end_exon: end_exons.value(i),
                start_vertex: start_verts.value(i),
                end_vertex: end_verts.value(i),
                n_exons: n_exons_col.value(i),
                chromosome: chroms.value(i).to_string(),
                start_pos: starts.value(i),
                end_pos: ends.value(i),
            });
        }
    }
    Ok(result)
}

fn load_edges_from_parquet(path: &Path) -> IsoForgeResult<Vec<EdgeRecord>> {
    if !path.exists() {
        return Ok(Vec::new());
    }
    let batches = parquet_io::read_parquet_batches(path)?;
    let mut result = Vec::new();
    for batch in &batches {
        let ids = batch.column(0).as_any().downcast_ref::<Int64Array>().unwrap();
        let v1s = batch.column(1).as_any().downcast_ref::<Int64Array>().unwrap();
        let v2s = batch.column(2).as_any().downcast_ref::<Int64Array>().unwrap();
        let types = batch.column(3).as_any().downcast_ref::<StringArray>().unwrap();
        let strands = batch.column(4).as_any().downcast_ref::<StringArray>().unwrap();
        for i in 0..batch.num_rows() {
            result.push(EdgeRecord {
                edge_id: ids.value(i),
                v1: v1s.value(i),
                v2: v2s.value(i),
                edge_type: types.value(i).to_string(),
                strand: strands.value(i).to_string(),
            });
        }
    }
    Ok(result)
}

fn load_locations_from_parquet(path: &Path) -> IsoForgeResult<Vec<LocationRecord>> {
    if !path.exists() {
        return Ok(Vec::new());
    }
    let batches = parquet_io::read_parquet_batches(path)?;
    let mut result = Vec::new();
    for batch in &batches {
        let ids = batch.column(0).as_any().downcast_ref::<Int64Array>().unwrap();
        let chroms = batch.column(1).as_any().downcast_ref::<StringArray>().unwrap();
        let positions = batch.column(2).as_any().downcast_ref::<Int64Array>().unwrap();
        for i in 0..batch.num_rows() {
            result.push(LocationRecord {
                location_id: ids.value(i),
                chromosome: chroms.value(i).to_string(),
                position: positions.value(i),
            });
        }
    }
    Ok(result)
}


/// Distance in bp downstream of every start of its model beyond which a read without a TSS cluster is demoted.
const UNCLUSTERED_START_WINDOW: i64 = 100;

/// Demotes each end-validated read without a TSS cluster that starts more than `UNCLUSTERED_START_WINDOW` downstream
/// of every start of its model (the TSS peaks its clustered reads reach or, when none reaches one, the start peaks
/// of its own reads, `read_start_peaks`); a model without either keeps its reads. Returns the number of reads demoted.
pub(crate) fn demote_unclustered_degraded_starts(observed: &mut [ObservedRead]) -> usize {
    let orient = |r: &ObservedRead| if &*r.strand == "-" { -1 } else { 1 };
    // positions in transcript orientation (x = position on +, -position on -)
    let five = |r: &ObservedRead| { let s = orient(r); (s * r.start_pos).min(s * r.end_pos) };
    let mut peaks: FxHashMap<i64, Vec<i64>> = FxHashMap::default();
    let mut starts: FxHashMap<i64, Vec<i64>> = FxHashMap::default();
    for r in observed.iter().filter(|r| r.cluster_supported == Some(true)) {
        match r.tss_site {
            Some(t) => {
                let p = orient(r) * t;
                let v = peaks.entry(r.transcript_id).or_default();
                if !v.contains(&p) { v.push(p); }
            }
            None => starts.entry(r.transcript_id).or_default().push(five(r)),
        }
    }
    for (tid, mut x) in starts {
        if peaks.contains_key(&tid) { continue; }
        let p = read_start_peaks(&mut x);
        if !p.is_empty() { peaks.insert(tid, p); }
    }
    let mut n = 0;
    for r in observed.iter_mut() {
        if r.cluster_supported != Some(true) || r.tss_site.is_some() { continue; }
        let Some(v) = peaks.get(&r.transcript_id) else { continue };
        let x5 = five(r);
        // downstream of every peak by more than the window
        if v.iter().all(|&p| x5 - p > UNCLUSTERED_START_WINDOW) {
            r.cluster_supported = Some(false);
            n += 1;
        }
    }
    n
}

const START_PEAK_WIDTH: i64 = 20;
const START_PEAK_MIN_SHARE: f64 = 0.1;
const START_PEAK_MIN_READS: usize = 3;

/// Start peaks of a model's reads (5' ends in transcript orientation): windows of `START_PEAK_WIDTH` + 1 bp holding
/// at least `START_PEAK_MIN_SHARE` of the reads and `START_PEAK_MIN_READS`, densest first (the most upstream of
/// equal ones), each at its window's median start.
fn read_start_peaks(x: &mut Vec<i64>) -> Vec<i64> {
    x.sort_unstable();
    let need = ((x.len() as f64 * START_PEAK_MIN_SHARE).ceil() as usize).max(START_PEAK_MIN_READS);
    let mut peaks = Vec::new();
    loop {
        let (mut best, mut at, mut j) = (0usize, 0usize, 0usize);
        for i in 0..x.len() {
            if j < i + 1 { j = i + 1; }
            while j < x.len() && x[j] - x[i] <= START_PEAK_WIDTH { j += 1; }
            if j - i > best { best = j - i; at = i; }
        }
        if best < need { break; }
        peaks.push(x[at + best / 2]);
        x.drain(at..at + best);
    }
    peaks
}

fn load_observed_from_parquet(path: &Path) -> IsoForgeResult<Vec<ObservedRead>> {
    if !path.exists() {
        return Ok(Vec::new());
    }
    // Decodes one batch at a time into rows sized from the file's row count.
    let builder = ParquetRecordBatchReaderBuilder::try_new(File::open(path)?)?;
    let mut result = Vec::with_capacity(builder.metadata().file_metadata().num_rows() as usize);
    let mut strings = Interner::default();
    for batch in builder.build()? {
        let batch = batch?;
        let datasets = batch.column(1).as_any().downcast_ref::<StringArray>().unwrap();
        let gene_ids = batch.column(2).as_any().downcast_ref::<Int64Array>().unwrap();
        let transcript_ids = batch.column(3).as_any().downcast_ref::<Int64Array>().unwrap();
        let chroms = batch.column(4).as_any().downcast_ref::<StringArray>().unwrap();
        let start_positions = batch.column(5).as_any().downcast_ref::<Int64Array>().unwrap();
        let end_positions = batch.column(6).as_any().downcast_ref::<Int64Array>().unwrap();
        let strands = batch.column(7).as_any().downcast_ref::<StringArray>().unwrap();
        let read_lengths = batch.column(8).as_any().downcast_ref::<Int64Array>().unwrap();
        let n_exons_col = batch.column(9).as_any().downcast_ref::<Int32Array>().unwrap();
        let start_vertices = batch.column(10).as_any().downcast_ref::<Int64Array>().unwrap();
        let end_vertices = batch.column(11).as_any().downcast_ref::<Int64Array>().unwrap();
        let start_exons = batch.column(12).as_any().downcast_ref::<Int64Array>().unwrap();
        let end_exons = batch.column(13).as_any().downcast_ref::<Int64Array>().unwrap();
        let start_deltas = batch.column(14).as_any().downcast_ref::<Int64Array>().unwrap();
        let end_deltas = batch.column(15).as_any().downcast_ref::<Int64Array>().unwrap();
        let fractions_as = batch.column(16).as_any().downcast_ref::<Float64Array>().unwrap();
        let custom_labels = batch.column(17).as_any().downcast_ref::<StringArray>().unwrap();
        let allelic_labels = batch.column(18).as_any().downcast_ref::<StringArray>().unwrap();
        let start_supports = batch.column(19).as_any().downcast_ref::<StringArray>().unwrap();
        let end_supports = batch.column(20).as_any().downcast_ref::<StringArray>().unwrap();
        let has_pas_signals = batch.column(21).as_any().downcast_ref::<BooleanArray>().unwrap();
        let pas_motifs = batch.column(22).as_any().downcast_ref::<StringArray>().unwrap();
        let pas_distances = batch.column(23).as_any().downcast_ref::<Int64Array>().unwrap();
        let is_internal_primings = batch.column(24).as_any().downcast_ref::<BooleanArray>().unwrap();
        let tss_sites = batch.column(25).as_any().downcast_ref::<Int64Array>().unwrap();
        let tss_features = batch.column(26).as_any().downcast_ref::<StringArray>().unwrap();
        let polya_sites = batch.column(27).as_any().downcast_ref::<Int64Array>().unwrap();
        let polya_features = batch.column(28).as_any().downcast_ref::<StringArray>().unwrap();
        let polya_signals = batch.column(29).as_any().downcast_ref::<StringArray>().unwrap();
        let cluster_supporteds = batch.column(30).as_any().downcast_ref::<BooleanArray>().unwrap();
        // Optional column
        let polya_lengths = if batch.num_columns() > 31 {
            batch.column(31).as_any().downcast_ref::<Int64Array>()
        } else {
            None
        };
        let max_splice_shifts = if batch.num_columns() > 32 {
            batch.column(32).as_any().downcast_ref::<Int32Array>()
        } else {
            None
        };
        let cap_signatures = if batch.num_columns() > 33 {
            batch.column(33).as_any().downcast_ref::<BooleanArray>()
        } else {
            None
        };
        let polya_softclips = if batch.num_columns() > 34 {
            batch.column(34).as_any().downcast_ref::<Int64Array>()
        } else {
            None
        };

        for i in 0..batch.num_rows() {
            result.push(ObservedRead {
                dataset: strings.get(datasets.value(i)),
                gene_id: gene_ids.value(i),
                transcript_id: transcript_ids.value(i),
                chromosome: strings.get(chroms.value(i)),
                start_pos: start_positions.value(i),
                end_pos: end_positions.value(i),
                strand: strings.get(strands.value(i)),
                read_length: read_lengths.value(i),
                n_exons: n_exons_col.value(i),
                start_vertex: start_vertices.value(i),
                end_vertex: end_vertices.value(i),
                start_exon: start_exons.value(i),
                end_exon: end_exons.value(i),
                start_delta: if start_deltas.is_null(i) { None } else { Some(start_deltas.value(i)) },
                end_delta: if end_deltas.is_null(i) { None } else { Some(end_deltas.value(i)) },
                fraction_as: if fractions_as.is_null(i) { None } else { Some(fractions_as.value(i)) },
                custom_label: strings.opt(custom_labels, i),
                allelic_label: strings.opt(allelic_labels, i),
                start_support: strings.opt(start_supports, i),
                end_support: strings.opt(end_supports, i),
                has_pas_signal: if has_pas_signals.is_null(i) { None } else { Some(has_pas_signals.value(i)) },
                pas_motif: strings.opt(pas_motifs, i),
                pas_distance: if pas_distances.is_null(i) { None } else { Some(pas_distances.value(i)) },
                is_internal_priming: if is_internal_primings.is_null(i) { None } else { Some(is_internal_primings.value(i)) },
                tss_site: if tss_sites.is_null(i) { None } else { Some(tss_sites.value(i)) },
                tss_feature: strings.opt(tss_features, i),
                polya_site: if polya_sites.is_null(i) { None } else { Some(polya_sites.value(i)) },
                polya_feature: strings.opt(polya_features, i),
                polya_signal: strings.opt(polya_signals, i),
                cluster_supported: if cluster_supporteds.is_null(i) { None } else { Some(cluster_supporteds.value(i)) },
                polya_length: polya_lengths.and_then(|a| if a.is_null(i) { None } else { Some(a.value(i)) }),
                max_splice_shift: max_splice_shifts.and_then(|a| if a.is_null(i) { None } else { Some(a.value(i)) }),
                cap_signature: cap_signatures.and_then(|a| if a.is_null(i) { None } else { Some(a.value(i)) }),
                polya_softclip: polya_softclips.and_then(|a| if a.is_null(i) { None } else { Some(a.value(i)) }),
                fragment: 0,
            });
        }
    }
    Ok(result)
}

fn load_abundance_from_parquet(path: &Path) -> IsoForgeResult<Vec<AbundanceRecord>> {
    if !path.exists() {
        return Ok(Vec::new());
    }
    let batches = parquet_io::read_parquet_batches(path)?;
    let mut result = Vec::new();
    for batch in &batches {
        let gene_ids = batch.column(0).as_any().downcast_ref::<Int64Array>().unwrap();
        let transcript_ids = batch.column(1).as_any().downcast_ref::<Int64Array>().unwrap();
        let datasets = batch.column(2).as_any().downcast_ref::<StringArray>().unwrap();
        let counts = batch.column(3).as_any().downcast_ref::<Int64Array>().unwrap();
        for i in 0..batch.num_rows() {
            result.push(AbundanceRecord {
                gene_id: gene_ids.value(i),
                transcript_id: transcript_ids.value(i),
                dataset: datasets.value(i).to_string(),
                count: counts.value(i),
            });
        }
    }
    Ok(result)
}

/// variant -> (annotated parent, compatible) from end_variants.parquet; empty when the file is absent.
fn load_end_variants(path: &Path) -> IsoForgeResult<FxHashMap<i64, (i64, bool)>> {
    let mut out = FxHashMap::default();
    if !path.exists() {
        return Ok(out);
    }
    for batch in &parquet_io::read_parquet_batches(path)? {
        let ids = batch.column(0).as_any().downcast_ref::<Int64Array>().unwrap();
        let parents = batch.column(1).as_any().downcast_ref::<Int64Array>().unwrap();
        let ok = batch.column(2).as_any().downcast_ref::<BooleanArray>().unwrap();
        for i in 0..batch.num_rows() {
            out.insert(ids.value(i), (parents.value(i), ok.value(i)));
        }
    }
    Ok(out)
}

fn load_counters(annot_dir: &Path, results_dir: &Path) -> Counters {
    let mut data = HashMap::new();
    let annot_path = annot_dir.join("counters.json");
    if annot_path.exists() {
        if let Ok(text) = fs::read_to_string(&annot_path) {
            if let Ok(val) = serde_json::from_str::<serde_json::Value>(&text) {
                if let Some(obj) = val.as_object() {
                    for (k, v) in obj {
                        data.insert(k.clone(), v.clone());
                    }
                }
            }
        }
    }
    let run_path = results_dir.join("counters.json");
    if run_path.exists() {
        if let Ok(text) = fs::read_to_string(&run_path) {
            if let Ok(val) = serde_json::from_str::<serde_json::Value>(&text) {
                if let Some(obj) = val.as_object() {
                    for (k, v) in obj {
                        data.insert(k.clone(), v.clone());
                    }
                }
            }
        }
    }
    let genes = data.get("genes").and_then(|v| v.as_i64()).unwrap_or(0);
    let transcripts = data.get("transcripts").and_then(|v| v.as_i64()).unwrap_or(0);
    let id_prefix = data.get("id_prefix")
        .and_then(|v| v.as_str())
        .unwrap_or("ISOFOR")
        .to_string();
    Counters { genes, transcripts, id_prefix, data }
}

// ═══════════════════════════════════════════════════════════════════════════════
// PostDB — shared pre-loaded data for all post-processing steps
// ═══════════════════════════════════════════════════════════════════════════════

/// Data shared by all post-processing steps, loaded once. Annotation rows are
/// held as lookup tables (`AnnotationDicts`).
pub struct PostDB {
    pub genes: Vec<GeneRecord>,         // annot + novel, deduped
    pub transcripts: Vec<TranscriptRecord>, // annot + novel, deduped
    pub edges: Vec<EdgeRecord>,         // annot + novel, deduped
    pub locations: Vec<LocationRecord>, // annot + novel, deduped
    pub annot: AnnotationDicts,         // pre-built annotation lookup tables
    pub observed: Vec<ObservedRead>,
    /// observed.parquet, whose rows `observed` holds in order; read names are read back from it.
    pub observed_path: PathBuf,
    pub abundance: Vec<AbundanceRecord>,
    pub counters: Counters,
    // Derived indices
    pub datasets: Vec<String>,
    pub annot_name: String,
    pub id_prefix: String,
    pub n_places: usize,
    // lookup maps
    pub loc_map: FxHashMap<i64, i64>,       // location_id → position
    pub gene_map: FxHashMap<i64, usize>,    // gene_id → index in genes
    pub tx_map: FxHashMap<i64, usize>,      // transcript_id → index in transcripts
    /// Every model as a hyperedge: the positions of its vertices 5'->3' (two per exon), from the index's
    /// hyperedges.parquet and the run's novel_hyperedges.parquet where they agree with the record, else from its edges.
    pub chains: FxHashMap<i64, Vec<i64>>,
    /// How many of `chains` came from the hyperedge tables.
    pub chains_from_tables: usize,
}

/// The edges of a transcript record in order: first exon, the internal edges, last exon.
pub fn record_path(t: &TranscriptRecord) -> Vec<i64> {
    let mut path = vec![t.start_exon];
    if let Some(p) = t.jn_path.as_deref() {
        path.extend(p.split(',').filter_map(|s| s.parse::<i64>().ok()));
    }
    if path.len() > 1 || t.end_exon != t.start_exon {
        path.push(t.end_exon);
    }
    path
}

/// Exon boundaries 5'->3' from an edge path (its exon edges, at even places), None when a vertex is unknown.
pub fn chain_from_edges(path: &[i64], edge_ends: &FxHashMap<i64, (i64, i64)>, loc_map: &FxHashMap<i64, i64>) -> Option<Vec<i64>> {
    let mut pos = Vec::with_capacity(path.len() + 1);
    for e in path.iter().step_by(2) {
        let &(v1, v2) = edge_ends.get(e)?;
        pos.push(*loc_map.get(&v1)?);
        pos.push(*loc_map.get(&v2)?);
    }
    Some(pos)
}

impl PostDB {
    /// Load all data from the annotation and results Parquet directories.
    ///
    /// Annotation rows are streamed batch by batch into `AnnotationDicts`.
    pub fn load(annot_dir: &str, results_dir: &str) -> IsoForgeResult<Self> {
        let annot = Path::new(annot_dir);
        let results = Path::new(results_dir);

        log::info!("[PostDB] Loading annotation data from {:?}", annot);
        log::info!("[PostDB] Loading results data from {:?}", results);

        // Load genes (annot + novel, dedup)
        let mut genes = load_genes_from_parquet(&annot.join("genes.parquet"))?;
        let novel_genes = load_genes_from_parquet(&results.join("novel_genes.parquet"))?;
        let mut seen_gids: FxHashSet<i64> = genes.iter().map(|g| g.gene_id).collect();
        for g in novel_genes {
            if seen_gids.insert(g.gene_id) {
                genes.push(g);
            }
        }
        log::info!("[PostDB] {} genes loaded", genes.len());

        // Load transcripts (annot + novel, dedup)
        let mut transcripts = load_transcripts_from_parquet(&annot.join("transcripts.parquet"))?;
        let novel_txs = load_transcripts_from_parquet(&results.join("novel_transcripts.parquet"))?;
        let mut seen_tids: FxHashSet<i64> = transcripts.iter().map(|t| t.transcript_id).collect();
        for t in novel_txs {
            if seen_tids.insert(t.transcript_id) {
                transcripts.push(t);
            }
        }
        log::info!("[PostDB] {} transcripts loaded", transcripts.len());

        // Load edges (annot + novel, dedup)
        let mut edges = load_edges_from_parquet(&annot.join("edges.parquet"))?;
        let novel_edges = load_edges_from_parquet(&results.join("novel_edges.parquet"))?;
        let mut seen_eids: FxHashSet<i64> = edges.iter().map(|e| e.edge_id).collect();
        for e in novel_edges {
            if seen_eids.insert(e.edge_id) {
                edges.push(e);
            }
        }
        log::info!("[PostDB] {} edges loaded", edges.len());

        // Load locations (annot + novel, dedup)
        let mut locations = load_locations_from_parquet(&annot.join("locations.parquet"))?;
        let novel_locs = load_locations_from_parquet(&results.join("novel_locations.parquet"))?;
        let mut seen_lids: FxHashSet<i64> = locations.iter().map(|l| l.location_id).collect();
        for l in novel_locs {
            if seen_lids.insert(l.location_id) {
                locations.push(l);
            }
        }
        log::info!("[PostDB] {} locations loaded", locations.len());

        // Load counters
        let counters = load_counters(annot, results);
        let id_prefix = counters.id_prefix.clone();

        // ── Stream annotation data into lookup tables ────────────────────
        let exclude: HashSet<&str> = ["IsoForge", &id_prefix].into_iter().collect();
        let annot_name = discover_annot_name(
            &annot.join("gene_annotations.parquet"), &exclude,
        ).unwrap_or_else(|| "unknown".to_string());
        log::info!("[PostDB] Detected annotation source: {}", annot_name);

        let mut annot_dicts = AnnotationDicts::new();

        // Gene annotations (annot dir)
        stream_annot_into_dicts(
            &annot.join("gene_annotations.parquet"),
            &mut annot_dicts, &annot_name, true,
        )?;

        // Transcript annotations (annot dir)
        stream_annot_into_dicts(
            &annot.join("transcript_annotations.parquet"),
            &mut annot_dicts, &annot_name, false,
        )?;

        // Run-time annotations (results dir) — split gene/tx automatically
        stream_run_annotations_into_dicts(
            &results.join("annotations.parquet"),
            &mut annot_dicts, &annot_name,
        )?;

        log::info!("[PostDB] Annotation dicts built: gene_id={}, gene_name={}, \
                    tx_id={}, tx_name={}, gene_novelty={}, tx_novelty={}, \
                    ism_prefix={}, ism_suffix={}, known_tx={}",
                   annot_dicts.gene_id_map.len(), annot_dicts.gene_name_map.len(),
                   annot_dicts.tx_id_map.len(), annot_dicts.tx_name_map.len(),
                   annot_dicts.gene_novelty.len(), annot_dicts.tx_novelty.len(),
                   annot_dicts.ism_prefix.len(), annot_dicts.ism_suffix.len(),
                   annot_dicts.known_tx_ids.len());

        annot_dicts.variant_parent = load_end_variants(&results.join("end_variants.parquet"))?;
        log::info!("[PostDB] {} end variants of annotated transcripts", annot_dicts.variant_parent.len());
        annot_dicts.small_rna = crate::small_rna::load(results)?;
        if !annot_dicts.small_rna.is_empty() {
            log::info!("[PostDB] {} annotated small RNAs counted, not assembled", annot_dicts.small_rna.len());
        }

        // Load observed reads
        let observed_path = results.join("observed.parquet");
        let mut observed = load_observed_from_parquet(&observed_path)?;
        log::info!("[PostDB] {} observed reads loaded", observed.len());
        let n_demoted = demote_unclustered_degraded_starts(&mut observed);
        log::info!("[PostDB] Reads without a TSS cluster starting more than {} bp inside every start of their model, no longer end-validated: {}",
                   UNCLUSTERED_START_WINDOW, n_demoted);

        // Load abundance
        let abundance = load_abundance_from_parquet(&results.join("abundance.parquet"))?;
        log::info!("[PostDB] {} abundance records loaded", abundance.len());

        // Derive datasets (sorted unique)
        let mut ds_set: BTreeMap<String, ()> = BTreeMap::new();
        for r in &observed {
            ds_set.entry(r.dataset.to_string()).or_default();
        }
        let datasets: Vec<String> = ds_set.into_keys().collect();
        log::info!("[PostDB] {} datasets: {:?}", datasets.len(), datasets);

        let n_places = {
            let max_id = counters.genes.max(counters.transcripts);
            std::cmp::max(9, format!("{}", max_id).len() + 1)
        };

        // Build lookup maps
        let loc_map: FxHashMap<i64, i64> = locations.iter()
            .map(|l| (l.location_id, l.position))
            .collect();

        let gene_map: FxHashMap<i64, usize> = genes.iter()
            .enumerate()
            .map(|(i, g)| (g.gene_id, i))
            .collect();

        let tx_map: FxHashMap<i64, usize> = transcripts.iter()
            .enumerate()
            .map(|(i, t)| (t.transcript_id, i))
            .collect();

        // The hypergraph: every model's vertex chain
        let mut rows: FxHashMap<i64, crate::annotation_hypergraph::HyperedgeRow> = FxHashMap::default();
        for path in [annot.join("hyperedges.parquet"), results.join("novel_hyperedges.parquet")] {
            if path.exists() {
                rows.extend(crate::annotation_hypergraph::read_hyperedges(&path)?.into_iter().map(|r| (r.transcript_id, r)));
            }
        }
        let edge_ends: FxHashMap<i64, (i64, i64)> = edges.iter().map(|e| (e.edge_id, (e.v1, e.v2))).collect();
        let mut chains: FxHashMap<i64, Vec<i64>> = FxHashMap::default();
        let mut chains_from_tables = 0usize;
        for t in &transcripts {
            let path = record_path(t);
            let row = rows.remove(&t.transcript_id).filter(|r| r.edges == path && r.vertices.first() == Some(&t.start_vertex)
                && r.vertices.last() == Some(&t.end_vertex));
            let chain = match row {
                Some(r) => { chains_from_tables += 1; Some(r.positions) }
                None => chain_from_edges(&path, &edge_ends, &loc_map),
            };
            if let Some(c) = chain { chains.insert(t.transcript_id, c); }
        }
        log::info!("[PostDB] Hypergraph: {} model chains ({} from the hyperedge tables)", chains.len(), chains_from_tables);
        let n_fragments = mark_fragments(&mut observed, &transcripts, &tx_map, &chains);
        log::info!("[PostDB] Reads holding fewer exons than the model they are credited to (fragments: they place no end and are not whole reads of it): {}", n_fragments);

        log::info!("[PostDB] Database loaded successfully");

        Ok(PostDB {
            genes,
            transcripts,
            edges,
            locations,
            annot: annot_dicts,
            observed,
            observed_path,
            abundance,
            counters,
            datasets,
            annot_name,
            id_prefix,
            n_places,
            loc_map,
            gene_map,
            tx_map,
            chains,
            chains_from_tables,
        })
    }

    // ── Helper methods ───────────────────────────────────────────────────

    /// The 5' and 3' end of a model's chain; for a variant split from an annotated record, the ends its reads place.
    pub fn chain_ends(&self, tid: i64) -> Option<(i64, i64)> {
        let c = self.chains.get(&tid)?;
        let (five, three) = (*c.first()?, *c.last()?);
        let (r5, r3) = self.annot.variant_ends.get(&tid).copied().unwrap_or((None, None));
        Some((r5.unwrap_or(five), r3.unwrap_or(three)))
    }

    /// Whether `dataset` is one of the datasets being processed.
    pub fn in_datasets(&self, dataset: &str) -> bool {
        self.datasets.iter().any(|d| d == dataset)
    }

    /// Get exon edge IDs from a transcript record.
    fn get_exon_edge_ids(tx: &TranscriptRecord) -> Vec<i64> {
        let mut edges = vec![tx.start_exon];
        if let Some(ref jn) = tx.jn_path {
            if !jn.is_empty() {
                let parts: Vec<&str> = jn.split(',').collect();
                for (i, p) in parts.iter().enumerate() {
                    if i % 2 == 1 {  // exon edges at odd indices
                        if let Ok(eid) = p.parse::<i64>() {
                            edges.push(eid);
                        }
                    }
                }
            }
        }
        if tx.end_exon != *edges.last().unwrap_or(&0) {
            edges.push(tx.end_exon);
        }
        edges
    }

    /// Format gene ID for output
    fn format_gene_id(&self, gid: i64, gene_id_map: &FxHashMap<i64, String>) -> String {
        match gene_id_map.get(&gid) {
            Some(name) => name.clone(),
            None => format!("{}G{:0width$}", self.id_prefix, gid, width = self.n_places),
        }
    }

    /// Format transcript ID for output
    fn format_transcript_id(&self, tid: i64, tx_id_map: &FxHashMap<i64, String>) -> String {
        match tx_id_map.get(&tid) {
            Some(name) => name.clone(),
            None => format!("{}T{:0width$}", self.id_prefix, tid, width = self.n_places),
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// GENE/TRANSCRIPT NOVELTY MAPS (used by multiple steps)
// ═══════════════════════════════════════════════════════════════════════════════

const GENE_NOVELTY_MAP: &[(&str, &str)] = &[
    ("gene_status", "Known"),
    ("antisense_gene", "Antisense"),
    ("intergenic_novel", "Intergenic"),
    ("fusion_novel", "Fusion"),
];

const TX_NOVELTY_MAP: &[(&str, &str)] = &[
    ("transcript_status", "Known"),
    ("ISM_transcript", "ISM"),
    ("NIC_transcript", "NIC"),
    ("NNC_transcript", "NNC"),
    ("antisense_transcript", "Antisense"),
    ("intergenic_transcript", "Intergenic"),
    ("genomic_transcript", "Genomic"),
    ("fusion_transcript", "Fusion"),
];

// ═══════════════════════════════════════════════════════════════════════════════
// 1. SUMMARIZE DATASETS
// ═══════════════════════════════════════════════════════════════════════════════

pub fn run_summarize(db: &PostDB, output_prefix: &str) -> IsoForgeResult<String> {
    let gene_novelty = &db.annot.gene_novelty;
    let tx_novelty = &db.annot.tx_novelty;
    let ism_prefix = &db.annot.ism_prefix;
    let ism_suffix = &db.annot.ism_suffix;

    let outpath = format!("{}_isoforge_summary.tsv", output_prefix);
    let mut out = BufWriter::new(File::create(&outpath)?);

    // Header
    writeln!(out, "dataset\treads_annotated\tknown_genes\tantisense_genes\tother_novel_genes\t\
                   known_transcripts\tnovel_transcripts\tISMs\tprefix_ISMs\tsuffix_ISMs\t\
                   NICs\tNNCs\tantisense_transcripts\tgenomic_transcripts")?;

    for ds in &db.datasets {
        // Collect unique genes and transcripts for this dataset
        let mut unique_genes = FxHashSet::default();
        let mut unique_txs = FxHashSet::default();
        let mut n_reads = 0i64;

        for r in &db.observed {
            if *r.dataset == **ds {
                n_reads += 1;
                unique_genes.insert(r.gene_id);
                unique_txs.insert(r.transcript_id);
            }
        }

        let mut known_genes = 0i64;
        let mut antisense_genes = 0i64;
        for &gid in &unique_genes {
            match gene_novelty.get(&gid).map(|s| s.as_str()) {
                Some("Known") => known_genes += 1,
                Some("Antisense") => antisense_genes += 1,
                _ => {}
            }
        }
        let other_novel_genes = unique_genes.len() as i64 - known_genes - antisense_genes;

        let tx_set: FxHashSet<i64> = unique_txs.iter().copied().collect();
        let mut known_tx = 0i64;
        let mut ism_tx = 0i64;
        let mut nic_tx = 0i64;
        let mut nnc_tx = 0i64;
        let mut antisense_tx = 0i64;
        let mut genomic_tx = 0i64;
        for &tid in &unique_txs {
            match tx_novelty.get(&tid).map(|s| s.as_str()) {
                Some("Known") => known_tx += 1,
                Some("ISM") => ism_tx += 1,
                Some("NIC") => nic_tx += 1,
                Some("NNC") => nnc_tx += 1,
                Some("Antisense") => antisense_tx += 1,
                Some("Genomic") => genomic_tx += 1,
                _ => {}
            }
        }
        let prefix_ism = ism_prefix.intersection(&tx_set).count() as i64;
        let suffix_ism = ism_suffix.intersection(&tx_set).count() as i64;
        let novel_tx = unique_txs.len() as i64 - known_tx;

        writeln!(out, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                 ds, n_reads, known_genes, antisense_genes, other_novel_genes,
                 known_tx, novel_tx, ism_tx, prefix_ism, suffix_ism,
                 nic_tx, nnc_tx, antisense_tx, genomic_tx)?;
    }

    log::info!("[Post] Summary written to {}", outpath);
    Ok(outpath)
}

// ═══════════════════════════════════════════════════════════════════════════════
// 2. FILTER TRANSCRIPTS
// ═══════════════════════════════════════════════════════════════════════════════

/// Filter configuration
pub struct FilterConfig {
    /// Abundance counts only reads with end support, so a model needs one such read to be reported.
    pub quant_requires_end_support: bool,
    /// Remove novel spliced models starting on a parent's acceptor, and exempt those whose start is a spike of at
    /// least `start_min_llr` (log-likelihood ratio against the smear) from the phased-end and sub-path rules.
    pub start_model: bool,
    pub start_min_llr: f64,
    pub max_frac_a: f64,
    pub min_count: i64,
    pub min_datasets: i64,
    pub allow_genomic: bool,
    pub filter_known: bool,
    pub min_count_known: i64,
    pub min_tss_support: f64,
    /// A spliced novel model failing the TSS-feature rule passes when it starts more than `TSS_UPSTREAM_MIN` bp
    /// upstream of every reported exon ending at its first donor (there is one).
    pub tss_upstream_rescue: bool,
    /// A novel model failing the TSS-feature rule passes when at least
    /// `cap_tss_min_reads` of its reads with a TSS cluster, and this fraction of
    /// them, carry the cap signature (0 disables).
    pub cap_tss_min_fraction: f64,
    pub cap_tss_min_reads: i64,
    pub evidence_tss_rescue: bool,
    pub junction_genome: Option<String>,
    pub junction_sj_files: Vec<String>,
    /// Splice-site scores (`--splice_scores`): with them, an unannotated splice site of a novel model needs a score.
    pub junction_scores: Option<String>,
    pub min_count_supported: i64,
    pub require_phased_novel_ends: bool,
    /// Hold spliced novel models to `require_phased_novel_ends` too; without it only unspliced ones are.
    pub phased_ends_spliced: bool,
    /// Prune a novel model carrying an unannotated intron that is a minor
    /// variant of a better-supported intron within this many bp on both
    /// sites; 0 disables.
    pub junction_prune_radius: i64,
    /// Support ratio (this intron / alternative) below which it is pruned.
    pub junction_prune_ratio: f64,
    /// A novel model also needs at least this fraction of its gene's
    /// best-supported model, per dataset; 0 disables.
    pub min_novel_locus_fraction: f64,
    /// A read of an annotated transcript whose 3' end lies within this many bp
    /// of the annotated TES is fraction-A eligible whatever its downstream A
    /// content. 0 disables.
    pub frac_a_annotated_end_window: i64,
    /// Reads per dataset whose every internal splice site sits within
    /// `novel_exact_tol` bp of the model's vertices, required of a novel
    /// multi-exon model (0 disables).
    pub novel_exact_support: i64,
    pub novel_exact_tol: i64,
    /// Let ISM models use the evidence-based TSS rescue route.
    pub ism_evidence_rescue: bool,
    /// Compact start reads an ISM needs on that route.
    pub ism_rescue_min_reads: i64,
    /// An ISM's start must lie at least this far from every containing
    /// annotated parent's start (0 = no distance test).
    pub ism_rescue_min_distance: i64,
    /// Lower read floor for NIC models (every junction annotated); 0 = use min_count.
    pub min_count_nic: i64,
    /// A spliced model whose junctions are all annotated (NIC) passes the read floor with `WEAK_SOURCE_MIN_READS`
    /// reads when the reported model of its gene nearest to it (one whose chain contains its chain, else the one with
    /// fewest junctions not shared) has at most `WEAK_SOURCE_MAX_READS`.
    pub weak_source_tier: bool,
    /// A novel multi-exon model whose intron path is a contiguous sub-path of
    /// another reported model of the same gene needs at least this fraction of
    /// that model's reads (0 disables).
    pub sub_path_min_fraction: f64,
    /// A sub-path that lacks the longer model's 3' junctions passes the sub-path rule when its reads carry poly(A)
    /// tails at `SUB_PATH_TAIL_RATIO` of the longer model's rate. A path at the longer model's 3' end stays under
    /// the rule.
    pub sub_path_tail_evidence: bool,
    /// Drop a novel multi-exon model that differs from another reported model of its gene by one edit (see
    /// `SiblingEdit`), below that model's reads times the edit's ratio.
    pub sibling_edits: bool,
    /// The read ratios of the four edits: a lost terminal exon, a skipped microexon, a moved junction, an extra
    /// intron. An edit with ratio 0 drops nothing.
    pub sibling_ratio_terminal: f64,
    pub sibling_ratio_microexon: f64,
    pub sibling_ratio_shift: f64,
    pub sibling_ratio_intron: f64,
    /// An annotated transcript needs this many reads anchoring each of its
    /// annotated ends: a 5' end in a validated TSS cluster or within
    /// `known_end_evidence_window` bp of the annotated TSS, and a 3' end in a
    /// validated poly(A) cluster or within that window of the annotated TES.
    /// Enforced only for transcripts whose reads carry basecaller poly(A) tags.
    /// 0 disables.
    pub known_end_evidence: i64,
    pub known_end_evidence_window: i64,
    /// Lower per-dataset read floor for an end variant of an annotated intron
    /// chain (a model whose splice structure is entirely annotated and only whose
    /// terminal differs). Its reads must anchor a validated, promoter-annotated
    /// TSS in at least `end_variant_min_datasets` datasets. 0 disables.
    pub min_count_end_variant: i64,
    /// 0 means two thirds of the datasets, rounded up: 4 of 6 replicates, 1 of 1.
    pub end_variant_min_datasets: i64,
    /// Learned admission on the read hypergraph (`admission.rs`): the tree-ensemble JSON; None keeps the rules' report.
    pub admission_model: Option<String>,
    /// Rejected candidates at or above `admission_add` join the report (new intron chains only, one per chain, with
    /// counted reads); reported candidates below `admission_drop` leave it.
    pub admission_add: f64,
    pub admission_drop: f64,
    /// Thresholds of single novelty classes, (class, add, drop), in place of the two above.
    pub admission_by_class: Vec<(String, f64, f64)>,
    /// 3' evidence a joining candidate needs unless it is a novel combination of annotated junctions: shares of its
    /// reads with a poly(A) tail, in a TES cluster (at least), and internally primed (at most).
    pub admission_guard_tailed: f64,
    pub admission_guard_tes: f64,
    pub admission_guard_ip: f64,
    /// TSS-clustered share a joining candidate outside annotated genes needs; it must also be spliced.
    pub admission_guard_tss: f64,
    /// TSS-clustered share every joining candidate needs (0: none).
    pub admission_guard_tss_all: f64,
    /// Datasets a joining candidate's counted reads must come from, in a run of that many datasets or more.
    pub admission_min_datasets: i64,
    /// An annotated transcript also passes the known-evidence filter with this many reads whose two ends lie in
    /// validated clusters (a start cluster and a poly(A) cluster; not internally primed), in any datasets
    /// (0 disables).
    pub known_phased_reads: i64,
    /// The rules for novel models count reads per dataset. With this, the datasets of a run count as one library:
    /// the read floors, exact-junction and end-evidence counts and the locus fraction are taken over the run.
    pub pool_datasets: bool,
}

impl Default for FilterConfig {
    fn default() -> Self {
        FilterConfig {
            quant_requires_end_support: true,
            start_model: false,
            start_min_llr: 20.0,
            max_frac_a: 0.35,
            min_count: 5,
            min_datasets: 0,
            allow_genomic: false,
            filter_known: false,
            min_count_known: 0,
            min_tss_support: 0.0,
            tss_upstream_rescue: false,
            cap_tss_min_fraction: 0.15,
            cap_tss_min_reads: 3,
            evidence_tss_rescue: false,
            junction_genome: None,
            junction_scores: None,
            junction_sj_files: Vec::new(),
            min_count_supported: 0,
            require_phased_novel_ends: false,
            phased_ends_spliced: false,
            junction_prune_radius: 50,
            junction_prune_ratio: 0.7,
            min_novel_locus_fraction: 0.005,
            frac_a_annotated_end_window: 30,
            novel_exact_support: 3,
            novel_exact_tol: 0,
            ism_evidence_rescue: false,
            ism_rescue_min_reads: 5,
            ism_rescue_min_distance: 0,
            min_count_nic: 5,
            weak_source_tier: false,
            sub_path_min_fraction: 0.3,
            sub_path_tail_evidence: false,
            sibling_edits: true,
            sibling_ratio_terminal: 1.0,
            sibling_ratio_microexon: 0.5,
            sibling_ratio_shift: 0.1,
            sibling_ratio_intron: 0.05,
            known_end_evidence: 2,
            known_end_evidence_window: 50,
            min_count_end_variant: 3,
            end_variant_min_datasets: 0,
            admission_model: None,
            admission_add: 0.6,
            admission_drop: 0.5,
            admission_by_class: vec![("NIC".to_string(), 0.5, 0.1), ("NNC".to_string(), f64::INFINITY, 0.7)],
            admission_guard_tailed: 0.6,
            admission_guard_tes: 0.6,
            admission_guard_ip: 0.1,
            admission_guard_tss: 0.8,
            admission_guard_tss_all: 0.6,
            admission_min_datasets: 2,
            known_phased_reads: 0,
            pool_datasets: false,
        }
    }
}

/// The dataset a per-dataset count of a novel-model rule is kept under: the read's own, or one key for the whole run
/// when the datasets are pooled (`FilterConfig::pool_datasets`).
fn dataset_key(pool: bool, dataset: &str) -> &str {
    if pool { "" } else { dataset }
}

/// Return candidates with concentrated, independently phased end evidence,
/// excluding those contained in a longer expressed path with the same 3' introns.
fn evidence_supported_starts(db: &PostDB, reads: &[&ObservedRead], allow_nonreference: bool, max_span: i64, allow_ism: bool, ism_min_reads: i64, ism_min_distance: i64) -> FxHashSet<i64> {
    let mut groups: FxHashMap<i64, Vec<&ObservedRead>> = FxHashMap::default();
    for &r in reads { groups.entry(r.transcript_id).or_default().push(r); }
    let paths: FxHashMap<i64, Vec<i64>> = db.transcripts.iter().map(|t| {
        let path: Vec<i64> = t.jn_path.as_deref().unwrap_or("").split(',')
            .step_by(2).filter_map(|s| s.parse().ok()).collect();
        (t.transcript_id, path)
    }).collect();
    // Every intron must occur in the reference annotation on this route.
    let reference_introns: FxHashSet<i64> = paths.iter()
        .filter(|(tid, _)| db.annot.known_tx_ids.contains(tid))
        .flat_map(|(_, path)| path.iter().copied()).collect();
    let mut parents: FxHashMap<i64, Vec<i64>> = FxHashMap::default();
    for (&tid, rs) in &groups {
        if rs.len() >= 5 { parents.entry(rs[0].gene_id).or_default().push(tid); }
    }
    let mut rescued = FxHashSet::default();
    for (&tid, rs) in &groups {
        if rs.len() < 5 || rs[0].n_exons < 2 { continue; }
        let is_ism = db.annot.tx_novelty.get(&tid).is_some_and(|s| s.starts_with("ISM"));
        if is_ism && !allow_ism { continue; }
        let Some(path) = paths.get(&tid).filter(|p| !p.is_empty()) else { continue; };
        if !allow_nonreference && !path.iter().all(|e| reference_introns.contains(e)) { continue; }
        let shadow = parents.get(&rs[0].gene_id).is_some_and(|ps| ps.iter().any(|p| {
            let parent = &paths[p];
            parent.len() > path.len() && parent.ends_with(path)
        }));
        // With `allow_ism`, an ISM qualifies through its own compact, phased start.
        if shadow && !(allow_ism && is_ism) { continue; }
        let mut positions: Vec<i64> = rs.iter().filter(|r|
            r.cluster_supported == Some(true) && r.tss_site.is_some() && r.polya_site.is_some()
            && r.is_internal_priming != Some(true))
            .map(|r| r.start_pos).collect();
        positions.sort_unstable();
        let need = if is_ism { ism_min_reads.max(1) as usize } else { 3 };
        if is_ism && ism_min_distance > 0 && !positions.is_empty() {
            // The ISM's start must lie farther than `ism_min_distance` from the start
            // of every annotated transcript that contains its path.
            let start = positions[positions.len() / 2];
            let too_close = db.transcripts.iter().any(|p| {
                p.gene_id == rs[0].gene_id && db.annot.known_tx_ids.contains(&p.transcript_id)
                    && paths.get(&p.transcript_id).is_some_and(|pp| pp.len() > path.len() && pp.windows(path.len()).any(|w| w == path.as_slice()))
                    && db.loc_map.get(&p.start_vertex).is_some_and(|&ps| (ps - start).abs() < ism_min_distance)
            });
            if too_close { continue; }
        }
        if compact_support(&positions, max_span) >= need && compact_support(&positions, max_span)*10 >= rs.len()*7 {
            rescued.insert(tid);
        }
    }
    rescued
}

fn compact_support(sorted: &[i64], span: i64) -> usize {
    let mut left = 0; let mut best = 0;
    for right in 0..sorted.len() {
        while sorted[right]-sorted[left] > span { left += 1; }
        best = best.max(right-left+1);
    }
    best
}

fn supported_short_read_junction(info: &crate::genome_ref::SjInfo, strand: crate::types::Strand) -> bool {
    info.unique_reads >= 3 && info.strand == if strand == crate::types::Strand::Plus { 1 } else { 2 }
}

/// Models among `candidates` with an unannotated intron that no short-read junction supports and whose ends are not
/// GT-AG or GC-AG or, with splice-site scores, one of whose splice sites that no annotated intron uses has no score
/// (chromosomes without scores are exempt). Returns the models, and those among them that fail on the score.
fn unsupported_junction_models(db: &PostDB, mh: &crate::model_hypergraph::ModelHypergraph, candidates: &FxHashSet<i64>,
    genome: &crate::genome_ref::GenomeRef, sj: &[crate::genome_ref::SpliceJunctionDb], scores: Option<&crate::splice_scores::SpliceScores>)
    -> IsoForgeResult<(FxHashSet<i64>, FxHashSet<i64>)> {
    use crate::genome_ref::SpliceMotifClass;
    use crate::types::Strand;
    // the intron links of annotated hyperedges, and their splice sites by role
    let reference: FxHashSet<i64> = db.transcripts.iter()
        .filter(|t| db.annot.known_tx_ids.contains(&t.transcript_id))
        .flat_map(|t| mh.links(t.transcript_id).iter().copied()).collect();
    let mut sites: FxHashSet<(&str, bool, bool, i64)> = FxHashSet::default();
    if scores.is_some() {
        for &id in &reference {
            let Some(l) = mh.link(id) else { continue };
            let plus = l.strand == Strand::Plus;
            sites.insert((l.chromosome.as_str(), plus, true, l.donor));
            sites.insert((l.chromosome.as_str(), plus, false, l.acceptor));
        }
    }
    // per unannotated link: 0 supported, 1 no intron by its ends, 2 a new site without a score
    let mut verdict: FxHashMap<i64, u8> = FxHashMap::default();
    let (mut rejected, mut unscored) = (FxHashSet::default(), FxHashSet::default());
    for t in db.transcripts.iter().filter(|t| candidates.contains(&t.transcript_id)) {
        for &id in mh.links(t.transcript_id) {
            if reference.contains(&id) { continue; }
            let v = if let Some(&v) = verdict.get(&id) { v } else {
                let l = mh.link(id).ok_or_else(|| IsoForgeError::IntegrityError(format!("Missing intron link {id}")))?;
                let (a, b) = (l.donor, l.acceptor);
                let start = a.min(b)+1; let end = a.max(b)-1;
                let strand = l.strand;
                let short_reads = start < end && sj.iter().any(|db| db.get(&t.chromosome,start,end).is_some_and(|info| supported_short_read_junction(info,strand)));
                let v = if short_reads { 0 }
                    else if !(start < end && matches!(genome.check_splice_motif(&t.chromosome,start,end,strand)?.motif_class, SpliceMotifClass::GtAg | SpliceMotifClass::GcAg)) { 1 }
                    else if scores.is_some_and(|sc| sc.has_chrom(&t.chromosome) && [(true, l.donor), (false, l.acceptor)].into_iter().any(|(donor, pos)|
                        !sites.contains(&(t.chromosome.as_str(), strand == Strand::Plus, donor, pos)) && sc.score(&t.chromosome, strand, donor, pos).is_none())) { 2 }
                    else { 0 };
                verdict.insert(id, v); v
            };
            if v != 0 {
                rejected.insert(t.transcript_id);
                if v == 2 { unscored.insert(t.transcript_id); }
                break;
            }
        }
    }
    Ok((rejected, unscored))
}

/// The single edit turning a reported model into a novel one: a first or last exon lost (the novel model's chain is
/// the reported one's without its first or last junction, and its end sits within 10 bp of the exon that junction
/// led into), a microexon (<= 30 bp) skipped, one junction placed up to 20 bp away, or an extra intron inside an exon.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
enum SiblingEdit { FirstLost, LastLost, MicroexonSkip, Shifted, ExtraIntron }

impl SiblingEdit {
    /// The novel model is dropped below this fraction of the reported model's reads.
    fn max_ratio(self, c: &FilterConfig) -> f64 {
        match self { SiblingEdit::FirstLost | SiblingEdit::LastLost => c.sibling_ratio_terminal, SiblingEdit::MicroexonSkip => c.sibling_ratio_microexon,
                     SiblingEdit::Shifted => c.sibling_ratio_shift, SiblingEdit::ExtraIntron => c.sibling_ratio_intron }
    }
}

/// The one intron `m` has that `a` lacks (exon end, next exon start), when there is exactly one.
fn extra_intron(m: &[(i64, i64)], a: &[(i64, i64)]) -> Option<(i64, i64)> {
    let introns = |ex: &[(i64, i64)]| -> Vec<(i64, i64)> { ex.windows(2).map(|w| (w[0].1, w[1].0)).collect() };
    let ai = introns(a);
    let mut extra = introns(m).into_iter().filter(|x| !ai.contains(x));
    let first = extra.next();
    if extra.next().is_some() { None } else { first }
}

/// `m` and `a`: genomic exons, sorted; `plus`: strand.
fn sibling_edit(m: &[(i64, i64)], a: &[(i64, i64)], plus: bool) -> Option<SiblingEdit> {
    let junctions = |ex: &[(i64, i64)]| -> Vec<(i64, i64)> { ex.windows(2).map(|w| (w[0].1, w[1].0)).collect() };
    let (mc, ac) = (junctions(m), junctions(a));
    let (low_lost, high_lost) = if plus { (SiblingEdit::FirstLost, SiblingEdit::LastLost) } else { (SiblingEdit::LastLost, SiblingEdit::FirstLost) };
    if mc.len() + 1 == ac.len() {
        if mc[..] == ac[1..] && (m[0].0 - ac[0].1).abs() <= 10 { return Some(low_lost); }
        if mc[..] == ac[..ac.len() - 1] && (m[m.len() - 1].1 - ac[ac.len() - 1].0).abs() <= 10 { return Some(high_lost); }
        for i in 0..ac.len() - 1 {
            if mc[..i] == ac[..i] && mc[i] == (ac[i].0, ac[i + 1].1) && mc[i + 1..] == ac[i + 2..] {
                return (ac[i + 1].0 - ac[i].1 + 1 <= 30).then_some(SiblingEdit::MicroexonSkip);
            }
        }
    }
    if mc.len() == ac.len() + 1 {
        // an exon the reported model skips is no edit of this kind
        for i in 0..mc.len() - 1 {
            if ac[..i] == mc[..i] && ac[i] == (mc[i].0, mc[i + 1].1) && ac[i + 1..] == mc[i + 2..] { return None; }
        }
        for i in 0..mc.len() {
            if mc[..i] == ac[..i] && mc[i + 1..] == ac[i..] { return Some(SiblingEdit::ExtraIntron); }
        }
    }
    if mc.len() == ac.len() {
        let diff: Vec<usize> = (0..mc.len()).filter(|&i| mc[i] != ac[i]).collect();
        if let [i] = diff[..] {
            if (mc[i].0 - ac[i].0).abs() <= 20 && (mc[i].1 - ac[i].1).abs() <= 20 { return Some(SiblingEdit::Shifted); }
        }
    }
    None
}

/// Evidence for the start of a novel spliced model: `llr` is the log-likelihood ratio of a spike at the start
/// against the smear of truncated molecules alone; `on_parent_acceptor` marks a start on the acceptor of the
/// junction just upstream in a longer chain that contains this one.
struct StartEvidence {
    llr: f64,
    on_parent_acceptor: bool,
}

/// How far upstream of a reported model's exon with the same donor a first exon must start to be its own under
/// `tss_upstream_rescue`, in bp.
const TSS_UPSTREAM_MIN: i64 = 50;

/// Reads a NIC model needs under `weak_source_tier`, and the most its nearest reported model may have.
const WEAK_SOURCE_MIN_READS: i64 = 3;
const WEAK_SOURCE_MAX_READS: i64 = 5;

/// The nearest reported model of `tid` among `models` of its gene (intron paths, transcript order): one whose path
/// contains its path first, then the one with the fewest junctions not shared, and among those the one with most
/// reads. None when the gene has no other.
fn nearest_relative(tid: i64, mine: &[i64], models: &[(i64, Vec<i64>)], reads: &FxHashMap<i64, i64>) -> Option<i64> {
    let set: FxHashSet<i64> = mine.iter().copied().collect();
    models.iter().filter(|(o, _)| *o != tid).map(|(o, c)| {
        let contains = c.len() > mine.len() && c.windows(mine.len()).any(|w| w == mine);
        let other: FxHashSet<i64> = c.iter().copied().collect();
        ((!contains, set.symmetric_difference(&other).count(), -reads.get(o).copied().unwrap_or(0), *o), *o)
    }).min().map(|(_, o)| o)
}

/// Share of the longer model's poly(A) tail rate a 3'-short sub-path's reads need under `sub_path_tail_evidence`.
const SUB_PATH_TAIL_RATIO: f64 = 0.75;

/// Half-width of the window a start's reads are fitted over, and the width a start's own reads lie within.
const START_WINDOW: i64 = 150;
const START_CORE: i64 = 10;

/// Log-likelihood ratio of p x Laplace(mu, b) + (1 - p) x Uniform[lo, hi] against Uniform[lo, hi] for the
/// positions `x` (all within [lo, hi]), with p fitted by EM.
fn spike_llr(x: &[i64], mu: i64, lo: i64, hi: i64, b: f64) -> f64 {
    if x.is_empty() { return 0.0; }
    let uni = 1.0 / (hi - lo + 1) as f64;
    let lap: Vec<f64> = x.iter().map(|&v| (-((v - mu).abs() as f64) / b).exp() / (2.0 * b)).collect();
    let mut p = 0.5f64;
    for _ in 0..200 {
        let np = lap.iter().map(|&l| p * l / (p * l + (1.0 - p) * uni)).sum::<f64>() / lap.len() as f64;
        let done = (np - p).abs() < 1e-7;
        p = np;
        if done { break; }
    }
    lap.iter().map(|&l| (p * l + (1.0 - p) * uni).ln()).sum::<f64>() - lap.len() as f64 * uni.ln()
}

/// Spike-versus-smear evidence for the starts of novel spliced models.
///
/// The 5' ends of the reads on every model with the same intron chain are placed along the transcript (x grows
/// towards the 3' end on both strands). Around a model's start they are fitted, by EM on the spike weight p, as
/// p x Laplace(start, b) + (1 - p) x a flat smear over the window, with b the mean distance of annotated models'
/// reads to the annotated start. A start on the acceptor of the upstream junction of a longer chain containing
/// this one is flagged.
fn start_evidence(db: &PostDB, mh: &crate::model_hypergraph::ModelHypergraph, filt: &[&ObservedRead], candidates: &FxHashSet<i64>)
    -> FxHashMap<i64, StartEvidence> {
    let known = &db.annot.known_tx_ids;
    // Per model: its chain id, strand sign (+1 / -1), start in x and first donor in x.
    let mut chain_id: FxHashMap<Vec<i64>, usize> = FxHashMap::default();
    let mut chains: Vec<Vec<i64>> = Vec::new();
    let mut model: FxHashMap<i64, (usize, i64, i64, i64)> = FxHashMap::default();
    for t in &db.transcripts {
        let path: Vec<i64> = mh.links(t.transcript_id).to_vec();
        if path.is_empty() { continue; }
        let Some(first) = mh.link(path[0]) else { continue };
        let sign = if first.strand == crate::types::Strand::Minus { -1 } else { 1 };
        let (Some(&start), donor) = (db.loc_map.get(&t.start_vertex), first.donor) else { continue };
        let next = chains.len();
        let id = *chain_id.entry(path.clone()).or_insert_with(|| { chains.push(path); next });
        model.insert(t.transcript_id, (id, sign, sign * start, sign * donor));
    }
    // 5' ends of the reads of each chain, in x.
    let mut xs: Vec<Vec<i64>> = vec![Vec::new(); chains.len()];
    let mut offsets: Vec<f64> = Vec::new();
    for r in filt {
        let Some(&(id, sign, start_x, _)) = model.get(&r.transcript_id) else { continue };
        let five = if sign > 0 { r.start_pos.min(r.end_pos) } else { r.start_pos.max(r.end_pos) };
        let x = sign * five;
        xs[id].push(x);
        if known.contains(&r.transcript_id) && (x - start_x).abs() <= START_CORE {
            offsets.push((x - start_x).abs() as f64);
        }
    }
    for v in xs.iter_mut() { v.sort_unstable(); }
    let b = if offsets.is_empty() { 2.0 } else { (offsets.iter().sum::<f64>() / offsets.len() as f64).max(0.5) };
    log::info!("[Post] Start model: spike width {:.2} bp from {} reads at annotated starts", b, offsets.len());
    // Where each intron occurs after the first position of a chain, to find chains extending one upstream.
    let mut occurs: FxHashMap<i64, Vec<(usize, usize)>> = FxHashMap::default();
    for (id, c) in chains.iter().enumerate() {
        for (j, &e) in c.iter().enumerate().skip(1) { occurs.entry(e).or_default().push((id, j)); }
    }
    let mut out = FxHashMap::default();
    for &tid in candidates {
        let Some(&(id, _, mu, donor)) = model.get(&tid) else { continue };
        let c = &chains[id];
        let (lo, hi) = (mu - START_WINDOW, (mu + START_WINDOW).min(donor - 1).max(mu + 2));
        let v = &xs[id];
        let llr = spike_llr(&v[v.partition_point(|&x| x < lo)..v.partition_point(|&x| x <= hi)], mu, lo, hi, b);
        let on_parent_acceptor = occurs.get(&c[0]).is_some_and(|occ| occ.iter().any(|&(pid, j)| {
            let p = &chains[pid];
            p.len() >= j + c.len() && p[j..j + c.len()] == c[..]
                && mh.link(p[j - 1]).map(|l| l.acceptor)
                    .is_some_and(|acc| (acc - mu.abs()).abs() <= 3)
        }));
        out.insert(tid, StartEvidence { llr, on_parent_acceptor });
    }
    out
}

/// Topology-aware pruning of unannotated introns.
///
/// An unannotated intron in a novel model is a bulge when an alternative
/// intron (annotated, or novel and carried by a passing model with at least
/// three reads) lies within `radius` bp on both its donor and its acceptor
/// and outscores it:
///   score = log2(alt_reads / reads) + 1·annotated(alt) + motif(this, alt)
/// where motif is +1 when only the alternative is canonical and −1 when only
/// this intron is. The model is rejected when the best alternative scores at
/// least log2(1/ratio) and beats the runner-up by 0.5; ties are recorded,
/// not pruned. Annotated introns are never pruned. With `trace_path`, every
/// decision is written there with its alternatives.
fn bulge_junction_models(db: &PostDB, mh: &crate::model_hypergraph::ModelHypergraph, candidates: &FxHashSet<i64>,
    genome: Option<&crate::genome_ref::GenomeRef>, radius: i64, ratio: f64, trace_path: Option<&str>) -> IsoForgeResult<FxHashSet<i64>> {
    use crate::types::Strand;
    let mut model_reads: FxHashMap<i64, u64> = FxHashMap::default();
    for r in db.observed.iter().filter(|r| db.in_datasets(&r.dataset)) {
        *model_reads.entry(r.transcript_id).or_default() += 1;
    }
    let mut reference: FxHashSet<i64> = FxHashSet::default();
    let mut passing_introns: FxHashSet<i64> = FxHashSet::default();
    let mut support: FxHashMap<i64, u64> = FxHashMap::default();
    // Every intron in any observed model gets its coordinates once.
    let mut coord: FxHashMap<i64, (String, Strand, i64, i64)> = FxHashMap::default();
    for t in &db.transcripts {
        let known = db.annot.known_tx_ids.contains(&t.transcript_id);
        let n = model_reads.get(&t.transcript_id).copied().unwrap_or(0);
        for &id in mh.links(t.transcript_id) {
            if known { reference.insert(id); }
            if candidates.contains(&t.transcript_id) { passing_introns.insert(id); }
            *support.entry(id).or_default() += n;
            if !coord.contains_key(&id) {
                let l = mh.link(id).ok_or_else(|| IsoForgeError::IntegrityError(format!("Missing intron link {id}")))?;
                let strand: Strand = l.strand;
                // donor, then acceptor: the link's 5' vertex first on both strands
                coord.insert(id, (t.chromosome.clone(), strand, l.donor, l.acceptor));
            }
        }
    }
    // Neighbour index: per (chrom, strand), introns sorted by donor.
    let mut index: FxHashMap<(String, Strand), Vec<(i64, i64, i64)>> = FxHashMap::default();
    for (&id, (chrom, strand, d, a)) in &coord {
        if reference.contains(&id) || (passing_introns.contains(&id) && support.get(&id).copied().unwrap_or(0) >= 3) {
            index.entry((chrom.clone(), *strand)).or_default().push((*d, *a, id));
        }
    }
    for v in index.values_mut() { v.sort_unstable(); }
    let mut canonical_cache: FxHashMap<i64, bool> = FxHashMap::default();
    let mut canonical = |id: i64| -> IsoForgeResult<Option<bool>> {
        let Some(g) = genome else { return Ok(None) };
        if let Some(&c) = canonical_cache.get(&id) { return Ok(Some(c)); }
        let (chrom, strand, d, a) = &coord[&id];
        let (start, end) = (d.min(a) + 1, d.max(a) - 1);
        let c = start < end && g.check_splice_motif(chrom, start, end, *strand)?.is_canonical;
        canonical_cache.insert(id, c);
        Ok(Some(c))
    };
    let threshold = (1.0 / ratio).log2();
    let mut trace = match trace_path {
        Some(path) => Some(BufWriter::new(File::create(path)?)),
        None => None,
    };
    if let Some(trace) = trace.as_mut() {
        writeln!(trace, "transcript_id\tchrom\tstrand\tintron_id\tdonor\tacceptor\treads\tcanonical\talt_intron_id\talt_donor_shift\talt_acceptor_shift\talt_annotated\talt_reads\talt_canonical\tscore\trunner_up_score\tdecision")?;
    }
    let mut rejected = FxHashSet::default();
    for t in db.transcripts.iter().filter(|t| candidates.contains(&t.transcript_id) && !db.annot.known_tx_ids.contains(&t.transcript_id)) {
        for &id in mh.links(t.transcript_id) {
            if reference.contains(&id) { continue; }
            let (chrom, strand, d, a) = coord[&id].clone();
            let reads = support.get(&id).copied().unwrap_or(0);
            let mine = canonical(id)?;
            let Some(list) = index.get(&(chrom.clone(), strand)) else { continue };
            let lo = list.partition_point(|x| x.0 < d - radius);
            let mut scores: Vec<(f64, i64)> = Vec::new();
            for &(d2, a2, id2) in &list[lo..] {
                if d2 > d + radius { break; }
                if id2 == id || (a2 - a).abs() > radius { continue; }
                let annotated = reference.contains(&id2);
                let alt_reads = support.get(&id2).copied().unwrap_or(0);
                let mut score = (alt_reads.max(1) as f64 / reads.max(1) as f64).log2() + if annotated { 1.0 } else { 0.0 };
                if let (Some(m), Some(c)) = (mine, canonical(id2)?) { score += c as i8 as f64 - m as i8 as f64; }
                scores.push((score, id2));
            }
            if scores.is_empty() { continue; }
            scores.sort_by(|x, y| y.0.partial_cmp(&x.0).unwrap_or(std::cmp::Ordering::Equal));
            let (best, alt) = scores[0];
            let runner = scores.get(1).map(|x| x.0);
            let prune = best >= threshold && runner.map_or(true, |r| best - r >= 0.5);
            let (_, _, d2, a2) = &coord[&alt];
            let fmt = |c: Option<bool>| c.map_or(".".to_string(), |v| v.to_string());
            if let Some(trace) = trace.as_mut() {
                writeln!(trace, "{}\t{chrom}\t{}\t{id}\t{d}\t{a}\t{reads}\t{}\t{alt}\t{}\t{}\t{}\t{}\t{}\t{best:.3}\t{}\t{}",
                    t.transcript_id, strand.as_str(), fmt(mine), d2 - d, a2 - a, reference.contains(&alt),
                    support.get(&alt).copied().unwrap_or(0), fmt(canonical(alt)?),
                    runner.map_or(".".to_string(), |r| format!("{r:.3}")),
                    if prune { "pruned" } else if best >= threshold { "kept_tie" } else { "kept" })?;
            }
            if prune { rejected.insert(t.transcript_id); break; }
        }
    }
    Ok(rejected)
}

/// Which reads may support a model and count towards abundance.
///
/// A read is internally primed when the 20 nt downstream of its 3' end are
/// more than `max_frac_a` A and no polyadenylation signal hexamer lies within
/// 100 nt upstream. A read ending within `annotated_end_window` bp of any
/// annotated TES on its chromosome and strand is never excluded. Internally
/// primed reads are excluded from both model support and quantification.
pub struct InternalPrimingGate {
    max_frac_a: f64,
    window: i64,
    annotated_tes: FxHashMap<(String, String), Vec<i64>>,
}

impl InternalPrimingGate {
    pub fn new(db: &PostDB, max_frac_a: f64, annotated_end_window: i64) -> Self {
        let mut annotated_tes: FxHashMap<(String, String), Vec<i64>> = FxHashMap::default();
        if annotated_end_window > 0 {
            let strand_of: FxHashMap<i64, &str> = db.genes.iter().map(|g| (g.gene_id, g.strand.as_str())).collect();
            for t in db.transcripts.iter().filter(|t| db.annot.known_tx_ids.contains(&t.transcript_id)) {
                if let (Some(&tes), Some(strand)) = (db.loc_map.get(&t.end_vertex), strand_of.get(&t.gene_id)) {
                    annotated_tes.entry((t.chromosome.clone(), strand.to_string())).or_default().push(tes);
                }
            }
            for v in annotated_tes.values_mut() { v.sort_unstable(); v.dedup(); }
        }
        Self { max_frac_a, window: annotated_end_window, annotated_tes }
    }

    fn near_annotated_tes(&self, chromosome: &str, strand: &str, end_pos: i64) -> bool {
        let Some(ends) = self.annotated_tes.get(&(chromosome.to_string(), strand.to_string())) else { return false };
        let i = ends.partition_point(|&x| x < end_pos - self.window);
        ends.get(i).is_some_and(|&x| x <= end_pos + self.window)
    }

    pub fn eligible(&self, r: &ObservedRead) -> bool {
        self.eligible_fields(&r.chromosome, &r.strand, r.end_pos, r.fraction_as, r.has_pas_signal)
    }

    /// `eligible` for a read given by its fields: `end_pos` is its 3' end.
    pub fn eligible_fields(&self, chromosome: &str, strand: &str, end_pos: i64, fraction_as: Option<f64>, has_pas_signal: Option<bool>) -> bool {
        // Tolerance for floating-point fractions.
        fraction_as.unwrap_or(0.0) <= self.max_frac_a + 1e-6
            || has_pas_signal == Some(true)
            || (self.window > 0 && self.near_annotated_tes(chromosome, strand, end_pos))
    }
}

pub fn run_filter(db: &PostDB, output_file: &str, config: &FilterConfig) -> IsoForgeResult<String> {
    let filter_trace = std::env::var_os("ISOFORGE_FILTER_TRACE").is_some();
    if config.min_count_supported != 0 && config.min_count_supported < 5 {
        return Err(IsoForgeError::ConfigError("min_count_supported must be zero or at least five".into()));
    }
    let min_datasets = if config.min_datasets <= 0 {
        db.datasets.len() as i64
    } else {
        config.min_datasets
    };

    // Known transcript IDs
    let known_ids = &db.annot.known_tx_ids;
    // The models as hyperedges: their intron links and the links' vertices, which the rules below ask about
    let mh = crate::model_hypergraph::ModelHypergraph::new(db);

    let mut whitelist: FxHashSet<i64> = FxHashSet::default();

    if !config.filter_known {
        // Only include observed known transcripts
        let observed_txs: FxHashSet<i64> = db.observed.iter()
            .map(|r| r.transcript_id)
            .collect();
        whitelist = known_ids.intersection(&observed_txs).copied().collect();
    }

    // Remove internally primed reads (see `InternalPrimingGate`).
    let gate = InternalPrimingGate::new(db, config.max_frac_a, config.frac_a_annotated_end_window);
    let frac_a_eligible = |r: &ObservedRead| -> bool { gate.eligible(r) };
    let filt: Vec<&ObservedRead> = db.observed.iter()
        .filter(|r| db.in_datasets(&r.dataset) && frac_a_eligible(r))
        .collect();

    // Known evidence filter
    // annotated transcripts that pass it on whole molecules alone (`known_phased_reads`)
    let mut known_by_phased: FxHashSet<i64> = FxHashSet::default();
    if !config.filter_known && config.min_count_known > 0 && !whitelist.is_empty() {
        // Count per (known_tid, dataset)
        let mut known_counts: FxHashMap<(i64, &str), i64> = FxHashMap::default();
        for r in &filt {
            if whitelist.contains(&r.transcript_id) {
                *known_counts.entry((r.transcript_id, &r.dataset)).or_default() += 1;
            }
        }
        // Group by tid: count datasets with >= min_count_known
        let mut known_ds_counts: FxHashMap<i64, i64> = FxHashMap::default();
        for (&(tid, _ds), &count) in &known_counts {
            if count >= config.min_count_known {
                *known_ds_counts.entry(tid).or_default() += 1;
            }
        }
        let mut known_with_evidence: FxHashSet<i64> = known_ds_counts.iter()
            .filter(|(_, &n_ds)| n_ds >= min_datasets)
            .map(|(&tid, _)| tid)
            .collect();
        // Whole molecules: reads with both ends in validated clusters, over every dataset.
        if config.known_phased_reads > 0 {
            let mut whole: FxHashMap<i64, i64> = FxHashMap::default();
            for r in &filt {
                if whitelist.contains(&r.transcript_id) && r.cluster_supported == Some(true) && r.tss_site.is_some()
                    && r.polya_site.is_some() && r.is_internal_priming != Some(true) {
                    *whole.entry(r.transcript_id).or_default() += 1;
                }
            }
            for (&tid, &n) in &whole {
                if n >= config.known_phased_reads && known_with_evidence.insert(tid) { known_by_phased.insert(tid); }
            }
        }
        let known_evidence_filtered = whitelist.len() - known_with_evidence.len();
        whitelist = known_with_evidence;
        log::info!("[Post] Known evidence filter (≥{} reads in ≥{} datasets{}): {} pass, {} removed",
                   config.min_count_known, min_datasets,
                   if config.known_phased_reads > 0 { format!(", or ≥{} reads with both ends in clusters: {}", config.known_phased_reads, known_by_phased.len()) } else { String::new() },
                   whitelist.len(), known_evidence_filtered);
    }

    // Annotated transcripts need reads anchoring both annotated ends; applied
    // only where the basecaller reported poly(A) tails.
    let mut known_end_rejected: FxHashSet<i64> = FxHashSet::default();
    if config.known_end_evidence > 0 && !config.filter_known && !whitelist.is_empty() {
        let w = config.known_end_evidence_window;
        let mut ann_ends: FxHashMap<i64, (i64, i64)> = FxHashMap::default();
        for t in &db.transcripts {
            if !whitelist.contains(&t.transcript_id) { continue; }
            if let (Some(&s), Some(&e)) = (db.loc_map.get(&t.start_vertex), db.loc_map.get(&t.end_vertex)) {
                ann_ends.insert(t.transcript_id, (s, e));
            }
        }
        // per transcript: (reads carrying a poly(A) tag, 5' anchored, 3' anchored)
        let mut ev: FxHashMap<i64, (i64, i64, i64)> = FxHashMap::default();
        for r in &filt {
            let Some(&(tss, tes)) = ann_ends.get(&r.transcript_id) else { continue };
            let e = ev.entry(r.transcript_id).or_insert((0, 0, 0));
            if r.polya_length.is_some() { e.0 += 1; }
            // A 5' anchor is a read starting at the annotated TSS or one in a TSS
            // cluster located there. A 3' anchor is a read in a validated poly(A)
            // cluster or within the window of the annotated TES.
            if (r.start_pos - tss).abs() <= w || r.tss_site.is_some_and(|s| (s - tss).abs() <= w) { e.1 += 1; }
            if (r.end_pos - tes).abs() <= w || r.polya_site.is_some() { e.2 += 1; }
        }
        // Reads with both ends in validated clusters, counted per transcript;
        // transcripts carrying them are exempt from this rule.
        let mut phased: FxHashMap<i64, i64> = FxHashMap::default();
        for r in &filt {
            if r.cluster_supported == Some(true) && r.tss_site.is_some() && r.polya_site.is_some()
                && r.is_internal_priming != Some(true) {
                *phased.entry(r.transcript_id).or_default() += 1;
            }
        }
        for &tid in &whitelist {
            let (tagged, five, three) = ev.get(&tid).copied().unwrap_or((0, 0, 0));
            if phased.get(&tid).copied().unwrap_or(0) >= config.known_end_evidence { continue; }
            // skips transcripts that passed on whole molecules
            if known_by_phased.contains(&tid) { continue; }
            // The rule does not apply to transcripts without basecaller tail evidence.
            if tagged == 0 { continue; }
            if five < config.known_end_evidence || three < config.known_end_evidence {
                known_end_rejected.insert(tid);
            }
        }
        whitelist.retain(|tid| !known_end_rejected.contains(tid));
        log::info!("[Post] Annotated end evidence (≥{} reads anchoring each end, window {} bp): {} models removed",
                   config.known_end_evidence, w, known_end_rejected.len());
    }

    // Novelty lookup
    let tx_novelty = &db.annot.tx_novelty;

    // Novelty categories to exclude
    let mut exclude_cats: HashSet<&str> = HashSet::new();
    if !config.allow_genomic {
        exclude_cats.insert("Genomic");
    }

    // Filter novel reads
    let novel_reads: Vec<&&ObservedRead> = if config.filter_known {
        filt.iter().collect()
    } else {
        filt.iter().filter(|r| !known_ids.contains(&r.transcript_id)).collect()
    };

    // Count per (transcript, dataset); with pooled datasets, per transcript over the run.
    let pool = config.pool_datasets && db.datasets.len() > 1;
    let ds_keys: Vec<&str> = if pool { vec![""] } else { db.datasets.iter().map(|d| d.as_str()).collect() };
    let mut counts: FxHashMap<(i64, &str), i64> = FxHashMap::default();
    for r in &novel_reads {
        *counts.entry((r.transcript_id, dataset_key(pool, &r.dataset))).or_default() += 1;
    }
    if pool { log::info!("[Post] Novel-model rules count reads over the run's {} datasets (pooled)", db.datasets.len()); }

    // All novel transcripts observed
    let _all_novel_observed: FxHashSet<i64> = novel_reads.iter()
        .map(|r| r.transcript_id)
        .collect();

    let supported_starts = if config.evidence_tss_rescue || config.min_count_supported > 0 {
        evidence_supported_starts(db, &filt, false, 10, false, 5, 0)
    } else { FxHashSet::default() };

    // The weak-source tier: NIC models below the floor whose nearest reported relative is weak.
    let weak_source: FxHashSet<i64> = if config.weak_source_tier {
        let mut reads_of: FxHashMap<i64, i64> = FxHashMap::default();
        for r in &filt { *reads_of.entry(r.transcript_id).or_default() += 1; }
        let floor_passing: FxHashSet<i64> = counts.iter().filter(|(_, &c)| c >= config.min_count).map(|(&(t, _), _)| t).collect();
        let chain = |t: &TranscriptRecord| -> Vec<i64> { mh.links(t.transcript_id).to_vec() };
        // the reported models so far, by gene: observed annotated ones and novel ones above the floor
        let mut by_gene: FxHashMap<i64, Vec<(i64, Vec<i64>)>> = FxHashMap::default();
        for t in &db.transcripts {
            if !(whitelist.contains(&t.transcript_id) || floor_passing.contains(&t.transcript_id)) { continue; }
            let c = chain(t);
            if !c.is_empty() { by_gene.entry(t.gene_id).or_default().push((t.transcript_id, c)); }
        }
        let mut out = FxHashSet::default();
        for (&(tid, _), &count) in &counts {
            if count < WEAK_SOURCE_MIN_READS || count >= config.min_count || out.contains(&tid) { continue; }
            if !tx_novelty.get(&tid).is_some_and(|c| c == "NIC") { continue; }
            let Some(&i) = db.tx_map.get(&tid) else { continue };
            let t = &db.transcripts[i];
            let mine = chain(t);
            if mine.is_empty() { continue; }
            if let Some(o) = by_gene.get(&t.gene_id).and_then(|v| nearest_relative(tid, &mine, v, &reads_of)) {
                if reads_of.get(&o).copied().unwrap_or(0) <= WEAK_SOURCE_MAX_READS { out.insert(tid); }
            }
        }
        log::info!("[Post] Weak-source tier: {} NIC models with {}-{} reads whose nearest reported model has at most {}",
                   out.len(), WEAK_SOURCE_MIN_READS, config.min_count - 1, WEAK_SOURCE_MAX_READS);
        out
    } else { FxHashSet::default() };

    // minCount filter
    let mut counts_pass_min: FxHashMap<(i64, &str), i64> = FxHashMap::default();
    for (&(tid, ds), &count) in &counts {
        let nic_tier = config.min_count_nic > 0 && count >= config.min_count_nic
            && tx_novelty.get(&tid).is_some_and(|c| c == "NIC");
        let weak_tier = count >= WEAK_SOURCE_MIN_READS && weak_source.contains(&tid);
        if count >= config.min_count || nic_tier || weak_tier || (config.min_count_supported > 0
            && count >= config.min_count_supported && supported_starts.contains(&tid)) {
            counts_pass_min.insert((tid, ds), count);
        }
    }

    // Count datasets per transcript
    let mut ds_counts: FxHashMap<i64, i64> = FxHashMap::default();
    for &(tid, _) in counts_pass_min.keys() {
        *ds_counts.entry(tid).or_default() += 1;
    }

    let mut passing_txs: FxHashSet<i64> = ds_counts.iter()
        .filter(|(_, &n)| n >= min_datasets)
        .map(|(&tid, _)| tid)
        .collect();

    // Starts from full-length reads: on when the share of spliced reads with a poly(A) tail lies between
    // `FULL_LENGTH_START_RATE` and `FULL_LENGTH_START_MAX_RATE`; a start is then established by tailed reads.
    let starts_from_tails = {
        // Annotated starts per chromosome and strand, from the transcripts in the index.
        let mut starts: FxHashMap<(String, String), Vec<i64>> = FxHashMap::default();
        for tx in &db.transcripts {
            if !known_ids.contains(&tx.transcript_id) { continue; }
            let Some(strand) = db.gene_map.get(&tx.gene_id).map(|&i| normalize_strand(&db.genes[i].strand)) else { continue };
            let (Some(&a), Some(&b)) = (db.loc_map.get(&tx.start_vertex), db.loc_map.get(&tx.end_vertex)) else { continue };
            let s = if strand == "-" { a.max(b) } else { a.min(b) };
            starts.entry((tx.chromosome.clone(), strand)).or_default().push(s);
        }
        for v in starts.values_mut() { v.sort_unstable(); }
        let reaches = |r: &ObservedRead| starts.get(&(r.chromosome.to_string(), r.strand.to_string())).is_some_and(|v| {
            let i = v.partition_point(|&x| x < r.start_pos - ANNOTATED_START_WINDOW);
            i < v.len() && v[i] <= r.start_pos + ANNOTATED_START_WINDOW
        });
        let (mut tn, mut ta, mut un, mut ua) = (0usize, 0usize, 0usize, 0usize);
        for r in filt.iter().filter(|r| r.n_exons > 1) {
            let tailed = r.polya_length.or(r.polya_softclip).map_or(false, |v| v > 0);
            let at = reaches(r);
            if tailed { tn += 1; ta += at as usize; } else { un += 1; ua += at as usize; }
        }
        let n = tn + un;
        let (tr, ur) = (ta as f64 / tn.max(1) as f64, ua as f64 / un.max(1) as f64);
        let on = n > 0 && tn as f64 >= FULL_LENGTH_START_RATE * n as f64 && tn as f64 <= FULL_LENGTH_START_MAX_RATE * n as f64;
        log::info!("[Post] Starts from full-length reads: {} (5' ends at an annotated start: {:.1}% of tailed, {:.1}% of tail-less spliced reads)",
                   if on { "on" } else { "off" }, 100.0 * tr, 100.0 * ur);
        on
    };
    fn tailed_start(r: &ObservedRead) -> bool {
        r.cluster_supported == Some(true) && r.tss_site.is_some() && r.polya_length.or(r.polya_softclip).map_or(false, |v| v > 0)
    }
    // End variants of an annotated intron chain pass on a lower read floor when
    // a validated, promoter-annotated TSS is carried by enough reads in enough
    // datasets.
    let mut end_variant_admitted: FxHashSet<i64> = FxHashSet::default();
    if config.min_count_end_variant > 0 {
        let mut promoter_counts: FxHashMap<(i64, &str), i64> = FxHashMap::default();
        for r in &filt {
            let anchored = if starts_from_tails {
                tailed_start(r)
            } else {
                r.cluster_supported == Some(true) && r.tss_site.is_some()
                    && matches!(r.tss_feature.as_deref(), Some("five_prime_utr" | "intergenic"))
            };
            if anchored {
                *promoter_counts.entry((r.transcript_id, &*r.dataset)).or_default() += 1;
            }
        }
        let annotated_chains: FxHashSet<Vec<i64>> = db.transcripts.iter()
            .filter(|t| known_ids.contains(&t.transcript_id))
            .map(|t| mh.links(t.transcript_id).to_vec())
            .filter(|c| !c.is_empty())
            .collect();
        let need_datasets = if config.end_variant_min_datasets > 0 {
            config.end_variant_min_datasets
        } else {
            ((2 * db.datasets.len() as i64 + 2) / 3).max(1)
        };
        let candidates: FxHashSet<i64> = novel_reads.iter().map(|r| r.transcript_id).collect();
        for tid in candidates {
            if passing_txs.contains(&tid) || known_ids.contains(&tid) { continue; }
            let Some(&i) = db.tx_map.get(&tid) else { continue };
            let chain: Vec<i64> = mh.links(db.transcripts[i].transcript_id).to_vec();
            if chain.is_empty() || !annotated_chains.contains(&chain) { continue; }
            let n_ds = db.datasets.iter().filter(|ds| {
                promoter_counts.get(&(tid, ds.as_str())).copied().unwrap_or(0) >= config.min_count_end_variant
            }).count() as i64;
            if n_ds >= need_datasets {
                end_variant_admitted.insert(tid);
            }
        }
        passing_txs.extend(end_variant_admitted.iter());
        log::info!("[Post] End-variant tier (>= {} promoter-anchored reads in >= {} datasets, annotated intron chain): {} models admitted",
                   config.min_count_end_variant, need_datasets, end_variant_admitted.len());
    }

    // Filter by novelty category
    let mut to_remove = Vec::new();
    for &tid in &passing_txs {
        if let Some(cat) = tx_novelty.get(&tid) {
            if exclude_cats.contains(cat.as_str()) {
                to_remove.push(tid);
            }
        }
    }
    for tid in to_remove {
        passing_txs.remove(&tid);
    }

    // Exact junction support: a novel multi-exon model needs, in enough
    // datasets, reads whose internal splice sites all lie within
    // `novel_exact_tol` bp of the model's vertices.
    let mut exact_counts: FxHashMap<i64, i64> = FxHashMap::default();
    let mut per_ds: FxHashMap<(i64, &str), i64> = FxHashMap::default();
    for r in &filt {
        if r.max_splice_shift.is_some_and(|s| (s as i64) <= config.novel_exact_tol) {
            *per_ds.entry((r.transcript_id, dataset_key(pool, &r.dataset))).or_default() += 1;
            *exact_counts.entry(r.transcript_id).or_default() += 1;
        }
    }
    let mut exact_rejected = FxHashSet::default();
    if config.novel_exact_support > 0 {
        let mut ds_ok: FxHashMap<i64, i64> = FxHashMap::default();
        for (&(tid, _), &n) in &per_ds { if n >= config.novel_exact_support { *ds_ok.entry(tid).or_default() += 1; } }
        for &tid in &passing_txs {
            if known_ids.contains(&tid) { continue; }
            let multi = db.tx_map.get(&tid).is_some_and(|&i| db.transcripts[i].n_exons > 1);
            if multi && ds_ok.get(&tid).copied().unwrap_or(0) < min_datasets { exact_rejected.insert(tid); }
        }
        passing_txs.retain(|tid| !exact_rejected.contains(tid));
        log::info!("[Post] Novel exact junction support (>= {} reads within {} bp): {} models removed",
            config.novel_exact_support, config.novel_exact_tol, exact_rejected.len());
    }

    // A novel single-exon model needs end evidence: in enough datasets, reads whose start lies on a promoter-like
    // feature or whose end carries a polyadenylation signal. Same read floor and dataset count as the junction rule.
    let mut mono_rejected = FxHashSet::default();
    if config.novel_exact_support > 0 {
        let mut per_ds: FxHashMap<(i64, &str), i64> = FxHashMap::default();
        for r in &filt {
            let start_ok = matches!(r.tss_feature.as_deref(), Some("five_prime_utr" | "intergenic"));
            let end_ok = r.polya_signal.as_deref().is_some_and(|s| !s.is_empty())
                || r.has_pas_signal == Some(true);
            if start_ok || end_ok {
                *per_ds.entry((r.transcript_id, dataset_key(pool, &r.dataset))).or_default() += 1;
            }
        }
        let mut ds_ok: FxHashMap<i64, i64> = FxHashMap::default();
        for (&(tid, _), &n) in &per_ds { if n >= config.novel_exact_support { *ds_ok.entry(tid).or_default() += 1; } }
        for &tid in &passing_txs {
            if known_ids.contains(&tid) { continue; }
            let mono = db.tx_map.get(&tid).is_some_and(|&i| db.transcripts[i].n_exons == 1);
            if mono && ds_ok.get(&tid).copied().unwrap_or(0) < min_datasets { mono_rejected.insert(tid); }
        }
        passing_txs.retain(|tid| !mono_rejected.contains(tid));
        log::info!("[Post] Novel single-exon end evidence (>= {} reads on a promoter feature or with a poly(A) signal): {} models removed",
            config.novel_exact_support, mono_rejected.len());
    }

    // Removes a novel model with an annotated transcript's structure (its intron chain, or for a single-exon model
    // an overlapping annotated single-exon transcript) whose 3' cluster with most reads lies more than NEW_END_WINDOW
    // bp from every annotated 3' end of it and has a poly(A) signal on fewer than half of its reads.
    let mut new_end_rejected = FxHashSet::default();
    if starts_from_tails {
        let strand_of = |t: &TranscriptRecord| db.gene_map.get(&t.gene_id).map(|&i| normalize_strand(&db.genes[i].strand));
        let span = |t: &TranscriptRecord| Some((*db.loc_map.get(&t.start_vertex)?, *db.loc_map.get(&t.end_vertex)?));
        let mut chain_ends: FxHashMap<Vec<i64>, Vec<i64>> = FxHashMap::default();
        let mut mono_spans: FxHashMap<(String, String), Vec<(i64, i64, i64)>> = FxHashMap::default();
        for t in db.transcripts.iter().filter(|t| known_ids.contains(&t.transcript_id)) {
            let (Some(s), Some((a, b))) = (strand_of(t), span(t)) else { continue };
            let three = if s == "-" { a.min(b) } else { a.max(b) };
            let chain: Vec<i64> = mh.links(t.transcript_id).to_vec();
            if chain.is_empty() {
                mono_spans.entry((t.chromosome.clone(), s)).or_default().push((a.min(b), a.max(b), three));
            } else {
                chain_ends.entry(chain).or_default().push(three);
            }
        }
        // Per model and 3' cluster: reads, and reads whose cluster carries a signal.
        let mut sites: FxHashMap<i64, FxHashMap<i64, (usize, usize)>> = FxHashMap::default();
        for r in &filt {
            let Some(p) = r.polya_site else { continue };
            if known_ids.contains(&r.transcript_id) { continue; }
            let e = sites.entry(r.transcript_id).or_default().entry(p).or_default();
            e.0 += 1;
            if r.polya_signal.as_deref().is_some_and(|s| !s.is_empty()) { e.1 += 1; }
        }
        for &tid in &passing_txs {
            if known_ids.contains(&tid) { continue; }
            let (Some(&i), Some(per_site)) = (db.tx_map.get(&tid), sites.get(&tid)) else { continue };
            let t = &db.transcripts[i];
            let (Some(s), Some((a, b))) = (strand_of(t), span(t)) else { continue };
            let Some((&site, &(n, with_signal))) = per_site.iter().max_by_key(|(&p, &(n, _))| (n, -p)) else { continue };
            let chain: Vec<i64> = mh.links(t.transcript_id).to_vec();
            let annotated: Vec<i64> = if chain.is_empty() {
                let (lo, hi) = (a.min(b), a.max(b));
                mono_spans.get(&(t.chromosome.clone(), s)).map_or_else(Vec::new, |v| {
                    v.iter().filter(|(l, h, _)| *l <= hi && *h >= lo).map(|x| x.2).collect()
                })
            } else {
                chain_ends.get(&chain).cloned().unwrap_or_default()
            };
            if !annotated.is_empty() && annotated.iter().all(|&e| (e - site).abs() > NEW_END_WINDOW) && 2 * with_signal < n {
                new_end_rejected.insert(tid);
            }
        }
        passing_txs.retain(|tid| !new_end_rejected.contains(tid));
        log::info!("[Post] New 3' ends of annotated structures without a positioned poly(A) signal: {} models removed", new_end_rejected.len());
    }

    let rescue_ids = if config.evidence_tss_rescue {
        if config.junction_genome.is_some() {
            // Novel junctions may use this endpoint route; junction validation below
            // applies separately.
            evidence_supported_starts(db, &filt, true, 20, config.ism_evidence_rescue, config.ism_rescue_min_reads, config.ism_rescue_min_distance)
        } else { supported_starts.clone() }
    } else { FxHashSet::default() };
    let mut tss_rejected = FxHashSet::default();
    let mut tss_rescued = FxHashSet::default();
    let mut tss_upstream: FxHashSet<i64> = FxHashSet::default();
    let mut cap_rescued: FxHashSet<i64> = FxHashSet::default();
    // Novel models whose clustered-start reads carry untemplated 5' G's at the
    // configured rate (see `cap_signature` in the BAM reader) pass the TSS
    // feature rule and are exempt from the sub-path rule. Empty with `starts_from_tails`.
    let capped_starts: FxHashSet<i64> = if config.cap_tss_min_fraction > 0.0 && !starts_from_tails {
        let mut per_tx: FxHashMap<i64, (i64, i64)> = FxHashMap::default();
        for r in &filt {
            if known_ids.contains(&r.transcript_id) || r.tss_site.is_none() || r.cluster_supported != Some(true) { continue; }
            let e = per_tx.entry(r.transcript_id).or_default();
            e.0 += 1;
            if r.cap_signature == Some(true) { e.1 += 1; }
        }
        per_tx.into_iter()
            .filter(|&(_, (n, cap))| cap >= config.cap_tss_min_reads && cap as f64 >= config.cap_tss_min_fraction * n as f64)
            .map(|(tid, _)| tid)
            .collect()
    } else {
        FxHashSet::default()
    };

    // With the start model, a novel spliced model starting on the acceptor of the junction just upstream in a
    // longer chain containing its own is removed, and a model whose start spike reaches `start_min_llr`
    // (`start_evidence`) is exempt from the phased-end and sub-path rules.
    let mut start_rejected: FxHashSet<i64> = FxHashSet::default();
    let mut start_spiked: FxHashSet<i64> = FxHashSet::default();
    let is_multi = |tid: &i64| db.tx_map.get(tid).is_some_and(|&i| db.transcripts[i].n_exons > 1);
    if config.start_model {
        let candidates: FxHashSet<i64> = passing_txs.iter().copied()
            .filter(|t| !known_ids.contains(t) && is_multi(t))
            .collect();
        let evidence = start_evidence(db, &mh, &filt, &candidates);
        if filter_trace {
            let mut w = BufWriter::new(File::create(format!("{output_file}.start_model.tsv"))?);
            writeln!(w, "transcript_id\tllr\ton_parent_acceptor")?;
            let mut ids: Vec<_> = evidence.keys().copied().collect(); ids.sort_unstable();
            for t in ids { let e = &evidence[&t]; writeln!(w, "{t}\t{:.3}\t{}", e.llr, e.on_parent_acceptor)?; }
        }
        for (&tid, e) in &evidence {
            if e.on_parent_acceptor { start_rejected.insert(tid); }
            else if e.llr >= config.start_min_llr { start_spiked.insert(tid); }
        }
        passing_txs.retain(|t| !start_rejected.contains(t));
        log::info!("[Post] Start model: {} novel spliced models removed on a parent's acceptor, {} with a start spike \
                    (log LR >= {}) exempt from the phased-end and sub-path rules ({} checked)",
                   start_rejected.len(), start_spiked.len(), config.start_min_llr, candidates.len());
    }

    // Novel transcript TSS quality filter
    if config.min_tss_support > 0.0 {
        let novel_in_passing: FxHashSet<i64> = passing_txs.iter()
            .filter(|&&tid| {
                tx_novelty.get(&tid).map_or(false, |cat| cat != "Known")
            })
            .copied()
            .collect();

        if !novel_in_passing.is_empty() {
            let genuine_tss_features: HashSet<&str> = ["five_prime_utr", "intergenic"].into_iter().collect();

            // Count annotated and genuine reads per novel transcript
            let mut per_tx_annotated: FxHashMap<i64, i64> = FxHashMap::default();
            let mut per_tx_genuine: FxHashMap<i64, i64> = FxHashMap::default();
            // Tailed reads at the start, per transcript and dataset
            let mut per_tx_tailed: FxHashMap<(i64, &str), i64> = FxHashMap::default();

            for r in &filt {
                if novel_in_passing.contains(&r.transcript_id) {
                    if let Some(ref feat) = r.tss_feature {
                        *per_tx_annotated.entry(r.transcript_id).or_default() += 1;
                        if genuine_tss_features.contains(&**feat) {
                            *per_tx_genuine.entry(r.transcript_id).or_default() += 1;
                        }
                    }
                    if tailed_start(r) {
                        *per_tx_tailed.entry((r.transcript_id, dataset_key(pool, &r.dataset))).or_default() += 1;
                    }
                }
            }

            // the exons of the reported spliced models, by gene, for the own-promoter rescue
            let edge_ends: FxHashMap<i64, (i64, i64)> = db.edges.iter().map(|e| (e.edge_id, (e.v1, e.v2))).collect();
            let edge_plus: FxHashMap<i64, bool> = db.edges.iter().map(|e| (e.edge_id, e.strand == "+")).collect();
            let exons_of = |t: &TranscriptRecord| -> Option<(bool, Vec<(i64, i64)>)> {
                let plus = *edge_plus.get(&t.start_exon)?;
                let pos = crate::quant_em::model_positions(db, t, &edge_ends)?;
                let mut ex: Vec<(i64, i64)> = pos.chunks(2).map(|p| (p[0].min(p[1]), p[0].max(p[1]))).collect();
                ex.sort_unstable();
                Some((plus, ex))
            };
            let mut gene_exons: FxHashMap<i64, Vec<(i64, bool, Vec<(i64, i64)>)>> = FxHashMap::default();
            if config.tss_upstream_rescue {
                for t in &db.transcripts {
                    if !(passing_txs.contains(&t.transcript_id) || whitelist.contains(&t.transcript_id)) || t.n_exons < 2 { continue; }
                    if let Some((plus, ex)) = exons_of(t) { gene_exons.entry(t.gene_id).or_default().push((t.transcript_id, plus, ex)); }
                }
            }
            let own_first_exon = |tid: i64| -> bool {
                let Some(&i) = db.tx_map.get(&tid) else { return false };
                let t = &db.transcripts[i];
                let Some((plus, ex)) = exons_of(t) else { return false };
                if ex.len() < 2 { return false; }
                // the first exon in transcript order: its 5' end and its donor
                let (start, donor) = if plus { (ex[0].0, ex[0].1) } else { (ex[ex.len() - 1].1, ex[ex.len() - 1].0) };
                let mut nearest: Option<i64> = None;   // the least upstream offset of a reported exon at the donor
                for (o, oplus, oex) in gene_exons.get(&t.gene_id).map(|v| v.as_slice()).unwrap_or(&[]) {
                    if *o == tid || *oplus != plus { continue; }
                    for &(a, b) in oex {
                        let off = if plus && b == donor { Some(start - a) } else if !plus && a == donor { Some(b - start) } else { None };
                        if let Some(off) = off { nearest = Some(nearest.map_or(off, |n| n.max(off))); }
                    }
                }
                nearest.is_some_and(|off| off < -TSS_UPSTREAM_MIN)
            };
            let mut upstream_rescued = 0usize;
            let mut excluded_by_tss = 0;
            for &tid in &novel_in_passing {
                let multi = db.tx_map.get(&tid).is_some_and(|&i| db.transcripts[i].n_exons > 1);
                let fails = if starts_from_tails && multi {
                    // A spliced novel model's start is established by tailed reads.
                    let n_ds = ds_keys.iter()
                        .filter(|ds| per_tx_tailed.get(&(tid, **ds)).copied().unwrap_or(0) >= config.novel_exact_support)
                        .count() as i64;
                    n_ds < min_datasets
                } else {
                    let n_annotated = *per_tx_annotated.get(&tid).unwrap_or(&0);
                    if n_annotated == 0 {
                        continue; // no read with a TSS feature: the model passes
                    }
                    let n_genuine = *per_tx_genuine.get(&tid).unwrap_or(&0);
                    (n_genuine as f64 / n_annotated as f64) < config.min_tss_support
                };
                if fails {
                    // Capped starts pass wherever they lie.
                    if capped_starts.contains(&tid) {
                        cap_rescued.insert(tid);
                    } else if rescue_ids.contains(&tid) {
                        tss_rescued.insert(tid);
                    } else if config.tss_upstream_rescue && multi && own_first_exon(tid) {
                        tss_upstream.insert(tid);
                        upstream_rescued += 1;
                    } else {
                        passing_txs.remove(&tid);
                        tss_rejected.insert(tid);
                        excluded_by_tss += 1;
                    }
                }
            }
            log::info!("[Post] Novel TSS quality filter (≥{:.0}%): {} removed, {} kept on capped starts, {} on a first exon of \
                        their own ({} novel checked)",
                       config.min_tss_support * 100.0, excluded_by_tss, cap_rescued.len(), upstream_rescued, novel_in_passing.len());
        }
    }

    let mut phase_rejected = FxHashSet::default();
    if config.require_phased_novel_ends {
        let mut phased: FxHashMap<(i64, &str), usize> = FxHashMap::default();
        for r in &filt {
            if r.cluster_supported == Some(true) && r.tss_site.is_some() && r.polya_site.is_some()
                && r.is_internal_priming != Some(true) {
                *phased.entry((r.transcript_id, dataset_key(pool, &r.dataset))).or_default() += 1;
            }
        }
        let mut phased_datasets: FxHashMap<i64, i64> = FxHashMap::default();
        for ((tid, _), n) in phased { if n >= 3 { *phased_datasets.entry(tid).or_default() += 1; } }
        let spliced = |tid: i64| db.tx_map.contains_key(&tid) && mh.is_spliced(tid);
        for &tid in &passing_txs {
            if start_spiked.contains(&tid) { continue; }
            if !config.phased_ends_spliced && spliced(tid) { continue; }
            if !known_ids.contains(&tid) && phased_datasets.get(&tid).copied().unwrap_or(0) < min_datasets {
                phase_rejected.insert(tid);
            }
        }
        passing_txs.retain(|tid| !phase_rejected.contains(tid));
        log::info!("[Post] Novel phased end evidence ({}): {} models removed",
                   if config.phased_ends_spliced { "all novel models" } else { "unspliced novel models" }, phase_rejected.len());
    }
    let mut junction_rejected = FxHashSet::default();
    // the genome, short-read junctions and splice-site scores of the junction rule, kept for learned admission:
    // a candidate it would add is held to the same rule
    let mut junction_ctx: Option<(crate::genome_ref::GenomeRef, Vec<crate::genome_ref::SpliceJunctionDb>, Option<crate::splice_scores::SpliceScores>)> = None;
    if let Some(path) = &config.junction_genome {
        let genome = crate::genome_ref::GenomeRef::open(path)?;
        let sj_paths: Vec<_> = config.junction_sj_files.iter().map(std::path::PathBuf::from).collect();
        // Each file is loaded as its own junction table and checked separately.
        let sj: Vec<_> = sj_paths.iter().map(crate::genome_ref::SpliceJunctionDb::from_star_sj_tab)
            .collect::<IsoForgeResult<_>>()?;
        let scores = match &config.junction_scores {
            Some(p) => Some(crate::splice_scores::SpliceScores::load(std::path::Path::new(p))?),
            None => None,
        };
        let (rejected, unscored) = unsupported_junction_models(db, &mh, &passing_txs, &genome, &sj, scores.as_ref())?;
        junction_rejected = rejected;
        passing_txs.retain(|tid| !junction_rejected.contains(tid));
        log::info!("[Post] Unsupported novel junctions: {} models removed{}", junction_rejected.len(),
                   if scores.is_some() { format!(" ({} for a new splice site without a splice-site score)", unscored.len()) } else { String::new() });
        junction_ctx = Some((genome, sj, scores));
    }
    let mut prune_rejected = FxHashSet::default();
    if config.junction_prune_radius > 0 {
        let genome = match &config.junction_genome {
            Some(path) => Some(crate::genome_ref::GenomeRef::open(path)?),
            None => None,
        };
        let trace = filter_trace.then(|| format!("{output_file}.junction_pruning.tsv"));
        prune_rejected = bulge_junction_models(db, &mh, &passing_txs, genome.as_ref(), config.junction_prune_radius,
            config.junction_prune_ratio, trace.as_deref())?;
        passing_txs.retain(|tid| !prune_rejected.contains(tid));
        log::info!("[Post] Junction pruning (radius {} bp, ratio {}): {} models removed",
            config.junction_prune_radius, config.junction_prune_ratio, prune_rejected.len());
    }
    let mut subpath_rejected = FxHashSet::default();
    if config.sub_path_min_fraction > 0.0 {
        // each model's share of reads with a poly(A) tail
        let tail_rate: FxHashMap<i64, f64> = if config.sub_path_tail_evidence {
            let mut n: FxHashMap<i64, (u64, u64)> = FxHashMap::default();
            for r in &filt {
                let e = n.entry(r.transcript_id).or_default();
                e.0 += 1;
                if r.polya_softclip.unwrap_or(0) > 0 || r.polya_length.unwrap_or(0) > 0 { e.1 += 1; }
            }
            n.into_iter().map(|(t, (a, b))| (t, b as f64 / a.max(1) as f64)).collect()
        } else { FxHashMap::default() };
        let tailed_end = |mine: i64, longer: i64| -> bool {
            let (a, b) = (tail_rate.get(&mine).copied().unwrap_or(0.0), tail_rate.get(&longer).copied().unwrap_or(0.0));
            b > 0.0 && a >= SUB_PATH_TAIL_RATIO * b
        };
        // Read totals of every passing model, and their intron paths.
        let mut reads_of: FxHashMap<i64, i64> = FxHashMap::default();
        for r in &filt { if passing_txs.contains(&r.transcript_id) || whitelist.contains(&r.transcript_id) { *reads_of.entry(r.transcript_id).or_default() += 1; } }
        let mut by_gene: FxHashMap<i64, Vec<(i64, Vec<i64>)>> = FxHashMap::default();
        for t in &db.transcripts {
            if !(passing_txs.contains(&t.transcript_id) || whitelist.contains(&t.transcript_id)) { continue; }
            let path: Vec<i64> = mh.links(t.transcript_id).to_vec();
            if !path.is_empty() { by_gene.entry(t.gene_id).or_default().push((t.transcript_id, path)); }
        }
        for &tid in &passing_txs {
            if known_ids.contains(&tid) || capped_starts.contains(&tid) || start_spiked.contains(&tid) { continue; }
            let Some(&i) = db.tx_map.get(&tid) else { continue };
            let t = &db.transcripts[i]; let path: Vec<i64> = mh.links(t.transcript_id).to_vec();
            if path.is_empty() { continue; }
            let mine = reads_of.get(&tid).copied().unwrap_or(0) as f64;
            let parent_max = by_gene.get(&t.gene_id).map(|v| v.iter()
                .filter(|(o, p)| *o != tid && crate::model_hypergraph::contained(&path, p)
                    && !(config.sub_path_tail_evidence && !p.ends_with(&path) && tailed_end(tid, *o)))
                .map(|(o, _)| reads_of.get(o).copied().unwrap_or(0)).max().unwrap_or(0)).unwrap_or(0) as f64;
            if parent_max > 0.0 && mine < config.sub_path_min_fraction * parent_max { subpath_rejected.insert(tid); }
        }
        passing_txs.retain(|tid| !subpath_rejected.contains(tid));
        log::info!("[Post] Sub-path fraction {}: {} models removed", config.sub_path_min_fraction, subpath_rejected.len());
    }
    let mut sibling_rejected = FxHashSet::default();
    let mut sibling_parent: Vec<(i64, i64, SiblingEdit)> = Vec::new();
    if config.sibling_edits {
        let edge_ends: FxHashMap<i64, (i64, i64)> = db.edges.iter().map(|e| (e.edge_id, (e.v1, e.v2))).collect();
        let edge_plus: FxHashMap<i64, bool> = db.edges.iter().map(|e| (e.edge_id, e.strand == "+")).collect();
        let mut reads_of: FxHashMap<i64, i64> = FxHashMap::default();
        for r in &filt { if passing_txs.contains(&r.transcript_id) || whitelist.contains(&r.transcript_id) { *reads_of.entry(r.transcript_id).or_default() += 1; } }
        // Genomic exons of every reported multi-exon model, by gene.
        let mut by_gene: FxHashMap<i64, Vec<(i64, bool, Vec<(i64, i64)>, &str)>> = FxHashMap::default();
        // Introns of annotated transcripts (exon end, next exon start): an extra intron among them is no edit.
        let mut annotated_introns: FxHashSet<(&str, i64, i64)> = FxHashSet::default();
        for t in db.transcripts.iter().filter(|t| db.annot.known_tx_ids.contains(&t.transcript_id)) {
            let Some(pos) = crate::quant_em::model_positions(db, t, &edge_ends) else { continue };
            let ex = crate::model_hypergraph::genomic_exons(&pos);
            for w in ex.windows(2) { annotated_introns.insert((t.chromosome.as_str(), w[0].1, w[1].0)); }
        }
        for t in &db.transcripts {
            if !(passing_txs.contains(&t.transcript_id) || whitelist.contains(&t.transcript_id)) { continue; }
            let (Some(&plus), Some(pos)) = (edge_plus.get(&t.start_exon), crate::quant_em::model_positions(db, t, &edge_ends)) else { continue };
            let ex = crate::model_hypergraph::genomic_exons(&pos);
            if ex.len() < 2 { continue; }
            by_gene.entry(t.gene_id).or_default().push((t.transcript_id, plus, ex, t.chromosome.as_str()));
        }
        let mut by_edit: FxHashMap<SiblingEdit, usize> = FxHashMap::default();
        for models in by_gene.values() {
            for (tid, plus, ex, chrom) in models {
                if known_ids.contains(tid) || !passing_txs.contains(tid) { continue; }
                let mine = reads_of.get(tid).copied().unwrap_or(0) as f64;
                let hit = models.iter().filter(|(o, oplus, _, _)| o != tid && oplus == plus).find_map(|(o, _, oex, _)| {
                    let e = sibling_edit(ex, oex, *plus)?;
                    if e == SiblingEdit::ExtraIntron
                        && extra_intron(ex, oex).is_some_and(|(a, b)| annotated_introns.contains(&(*chrom, a, b))) { return None; }
                    // a model with a capped start is not held to a lost first exon
                    if e == SiblingEdit::FirstLost && capped_starts.contains(tid) { return None; }
                    // a lost terminal exon of an annotated transcript drops the model whatever its reads
                    let terminal = matches!(e, SiblingEdit::FirstLost | SiblingEdit::LastLost);
                    let ratio = e.max_ratio(config);
                    (terminal && ratio > 0.0 && known_ids.contains(o) || mine < ratio * reads_of.get(o).copied().unwrap_or(0) as f64)
                        .then_some((e, *o))
                });
                if let Some((e, o)) = hit {
                    sibling_rejected.insert(*tid); sibling_parent.push((*tid, o, e)); *by_edit.entry(e).or_default() += 1;
                }
            }
        }
        passing_txs.retain(|tid| !sibling_rejected.contains(tid));
        let mut detail: Vec<String> = by_edit.iter().map(|(e, n)| format!("{} {:?}", n, e)).collect();
        detail.sort();
        log::info!("[Post] One-edit siblings of reported models: {} models removed ({})", sibling_rejected.len(), detail.join(", "));
    }
    let mut locus_rejected = FxHashSet::default();
    if config.min_novel_locus_fraction > 0.0 {
        // A novel model must reach a fraction of its gene's best-supported
        // model in enough datasets; known models are never held to it.
        let mut counts: FxHashMap<(i64, &str), i64> = FxHashMap::default();
        for r in &filt { *counts.entry((r.transcript_id, dataset_key(pool, &r.dataset))).or_default() += 1; }
        let gene_of = |tid: i64| db.tx_map.get(&tid).map(|&i| db.transcripts[i].gene_id);
        let mut gene_max: FxHashMap<(i64, &str), i64> = FxHashMap::default();
        for (&(tid, ds), &n) in &counts {
            if let Some(g) = gene_of(tid) { let m = gene_max.entry((g, ds)).or_default(); *m = (*m).max(n); }
        }
        for &tid in &passing_txs {
            if known_ids.contains(&tid) { continue; }
            let Some(g) = gene_of(tid) else { continue };
            let ok = ds_keys.iter().filter(|ds| {
                let need = (config.min_novel_locus_fraction * gene_max.get(&(g, **ds)).copied().unwrap_or(0) as f64).ceil() as i64;
                counts.get(&(tid, **ds)).copied().unwrap_or(0) >= need.max(1)
            }).count() as i64;
            if ok < min_datasets { locus_rejected.insert(tid); }
        }
        passing_txs.retain(|tid| !locus_rejected.contains(tid));
        log::info!("[Post] Novel locus fraction {}: {} models removed", config.min_novel_locus_fraction, locus_rejected.len());
    }
    whitelist.extend(passing_txs.iter());

    // When counting requires end support, a model is reported only with a read that abundance counts: end-validated
    // and not internally primed, and for a novel model a whole read (`WholeReads`).
    let mut uncounted: FxHashSet<i64> = FxHashSet::default();
    if config.quant_requires_end_support {
        let novel = |t: i64| tx_novelty.get(&t).is_none_or(|x| x != "Known");
        let whole = WholeReads::new(db, &build_tss_tes_info(db, &whitelist));
        let counted: FxHashSet<i64> = db.observed.iter()
            .filter(|r| db.in_datasets(&r.dataset) && r.cluster_supported == Some(true) && gate.eligible(r)
                    && (!novel(r.transcript_id) || whole.holds(db, r)))
            .map(|r| r.transcript_id)
            .collect();
        uncounted = whitelist.iter().copied().filter(|t| !counted.contains(t)).collect();
        whitelist.retain(|t| counted.contains(t));
        log::info!("[Post] Models without a counted read (none with end support; novel models: no whole read): {} removed", uncounted.len());
    }

    // One row per observed model, including rejected candidates: the rules' read counts and decision. Written with
    // ISOFORGE_FILTER_TRACE, and the input of learned admission.
    let admission_on = config.admission_model.is_some() || std::env::var_os("ISOFORGE_ADMISSION_DUMP").is_some();
    let mut ledger_rows: FxHashMap<i64, crate::admission::LedgerRow> = FxHashMap::default();
    if filter_trace || admission_on {
        let mut all_counts: FxHashMap<i64, (u64, u64, u64, u64, u64)> = FxHashMap::default();
        for r in db.observed.iter().filter(|r| db.in_datasets(&r.dataset)) {
            let n = all_counts.entry(r.transcript_id).or_default();
            n.0 += 1;
            if frac_a_eligible(r) { n.1 += 1; }
            if r.cluster_supported == Some(true) {
                n.2 += 1;
                if r.tss_site.is_some() && r.polya_site.is_some() && r.is_internal_priming != Some(true)
                    && frac_a_eligible(r) { n.4 += 1; }
            }
            if matches!(r.tss_feature.as_deref(), Some("five_prime_utr" | "intergenic")) { n.3 += 1; }
        }
        for (&tid, &n) in &all_counts {
            let decision = if exact_rejected.contains(&tid) { "reject_novel_exact_support" }
                else if phase_rejected.contains(&tid) { "reject_unphased_novel_ends" }
                else if junction_rejected.contains(&tid) { "reject_unsupported_novel_junction" }
                else if prune_rejected.contains(&tid) { "reject_junction_prune" }
                else if start_rejected.contains(&tid) { "reject_start_model" }
                else if subpath_rejected.contains(&tid) { "reject_sub_path_fraction" }
                else if sibling_rejected.contains(&tid) { "reject_sibling_edit" }
                else if locus_rejected.contains(&tid) { "reject_novel_locus_fraction" }
                else if known_end_rejected.contains(&tid) { "reject_known_end_evidence" }
                else if new_end_rejected.contains(&tid) { "reject_new_end_without_signal" }
                else if uncounted.contains(&tid) { "reject_no_counted_reads" }
                else if end_variant_admitted.contains(&tid) && whitelist.contains(&tid) { "pass_end_variant_tier" }
                else if tss_rescued.contains(&tid) { "pass_evidence_tss_rescue" }
                else if tss_upstream.contains(&tid) && whitelist.contains(&tid) { "pass_own_first_exon" }
                else if cap_rescued.contains(&tid) && whitelist.contains(&tid) { "pass_cap_tss" }
                else if whitelist.contains(&tid) { "pass" }
                else if tss_rejected.contains(&tid) { "reject_tss_feature" }
                else if n.1 == 0 { "reject_fraction_a" }
                else if tx_novelty.get(&tid).is_some_and(|x| exclude_cats.contains(x.as_str())) { "reject_category" }
                else { "reject_depth_or_dataset_support" };
            ledger_rows.insert(tid, crate::admission::LedgerRow {
                reads: n.0, frac_a_eligible: n.1, passes_read_gates: n.2, utr_or_intergenic_tss: n.3, phased_end_reads: n.4,
                exact_junction_reads: exact_counts.get(&tid).copied().unwrap_or(0).max(0) as u64, decision,
                novelty: tx_novelty.get(&tid).cloned().unwrap_or_else(|| "Unknown".to_string()),
            });
        }
    }
    if filter_trace {
        let mut ids: Vec<_> = ledger_rows.keys().copied().collect(); ids.sort_unstable();
        let mut ledger = BufWriter::new(File::create(format!("{output_file}.decisions.tsv"))?);
        writeln!(ledger, "transcript_id\tnovelty\treads\tfrac_a_eligible\tpasses_read_gates\tutr_or_intergenic_tss\tphased_end_reads\texact_junction_reads\tdecision")?;
        for tid in ids {
            let l = &ledger_rows[&tid];
            writeln!(ledger, "{tid}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", l.novelty, l.reads, l.frac_a_eligible, l.passes_read_gates,
                l.utr_or_intergenic_tss, l.phased_end_reads, l.exact_junction_reads, l.decision)?;
        }
    }

    // Learned admission on the read hypergraph (see admission.rs): after the rules, before the report is written.
    // With ISOFORGE_ADMISSION_DUMP the features are computed and written even without a model.
    let dump_path = std::env::var_os("ISOFORGE_ADMISSION_DUMP").map(|p| p.to_string_lossy().to_string());
    if config.admission_model.is_some() || dump_path.is_some() {
        let t0 = std::time::Instant::now();
        let model = match &config.admission_model { Some(path) => Some(crate::admission::AdmissionModel::load(path)?), None => None };
        if ledger_rows.is_empty() {
            log::warn!("[Post] Admission: no ledger rows; the rules' report is kept");
        } else {
            match crate::hyper_core::load_read_chains(&db.observed_path)? {
                None => log::warn!("[Post] Admission: the run database holds no read chains; the rules' report is kept"),
                Some(chains) if chains.chain_of.len() == db.observed.len() => {
                    // The features take read ends against the chains' own ends, not the sites the reads place.
                    let opts = crate::hyper_core::HyperOptions { read_ends: false, ..Default::default() };
                    let core = crate::hyper_core::HyperCore::build(db, &chains, &crate::hyper_core::Sites::default(), opts);
                    let t_core = t0.elapsed().as_secs_f64();
                    let feats = crate::admission::features(db, &core, &chains, &ledger_rows, &whitelist, dump_path.is_some());
                    let t_feat = t0.elapsed().as_secs_f64() - t_core;
                    // per candidate with reads that abundance counts: the datasets those reads come from
                    let counted: FxHashMap<i64, u32> = {
                        let pairs: FxHashSet<(i64, &str)> = db.observed.iter()
                            .filter(|r| db.in_datasets(&r.dataset) && r.validated() && gate.eligible(r))
                            .map(|r| (r.transcript_id, &*r.dataset)).collect();
                        let mut n: FxHashMap<i64, u32> = FxHashMap::default();
                        for (tid, _) in pairs { *n.entry(tid).or_default() += 1; }
                        n
                    };
                    // for a dump, the junction rule's verdict on every candidate (a column of the table)
                    let junction_bad: Option<FxHashSet<i64>> = match (&dump_path, &junction_ctx) {
                        (Some(_), Some((genome, sj, scores))) => {
                            let all: FxHashSet<i64> = feats.iter().map(|f| f.tid).collect();
                            Some(unsupported_junction_models(db, &mh, &all, genome, sj, scores.as_ref())?.0)
                        }
                        _ => None,
                    };
                    match &model {
                        None => {
                            if let Some(path) = &dump_path { crate::admission::dump(path, &feats, &[], &counted, junction_bad.as_ref())?; }
                            log::info!("[Post] Admission features of {} candidates written (no model) in {:.1}s", feats.len(), t0.elapsed().as_secs_f64());
                        }
                        Some(model) => {
                            let prm = crate::admission::AdmissionParams { add: config.admission_add, drop: config.admission_drop,
                                by_class: config.admission_by_class.clone(),
                                guard_tailed: config.admission_guard_tailed, guard_tes: config.admission_guard_tes, guard_ip: config.admission_guard_ip,
                                guard_tss: config.admission_guard_tss, guard_tss_all: config.admission_guard_tss_all,
                                min_datasets: config.admission_min_datasets.clamp(1, db.datasets.len().max(1) as i64) as u32 };
                            // a candidate the junction rule rejects cannot join: the selection is repeated without the
                            // rejected candidates until none of the chosen is rejected
                            let mut no_junction: FxHashSet<i64> = FxHashSet::default();
                            let (mut add, mut drop, mut probs) = crate::admission::select(model, &feats, &counted, &prm, &no_junction);
                            if let Some((genome, sj, scores)) = &junction_ctx {
                                loop {
                                    let bad = unsupported_junction_models(db, &mh, &add, genome, sj, scores.as_ref())?.0;
                                    if bad.is_empty() { break; }
                                    no_junction.extend(bad);
                                    (add, drop, probs) = crate::admission::select(model, &feats, &counted, &prm, &no_junction);
                                }
                            }
                            if let Some(path) = &dump_path { crate::admission::dump(path, &feats, &probs, &counted, junction_bad.as_ref())?; }
                            let mut log_rows = BufWriter::new(File::create(format!("{output_file}.admission.tsv"))?);
                            writeln!(log_rows, "transcript_id\tp\trules\tadmission")?;
                            for (tid, p) in &probs {
                                let rules = if whitelist.contains(tid) { "reported" } else { "rejected" };
                                let act = if add.contains(tid) { "added" } else if drop.contains(tid) { "dropped" } else if whitelist.contains(tid) { "kept" } else { "rejected" };
                                writeln!(log_rows, "{tid}\t{p:.6}\t{rules}\t{act}")?;
                            }
                            whitelist.retain(|t| !drop.contains(t));
                            whitelist.extend(add.iter().copied());
                            let t_all = t0.elapsed().as_secs_f64();
                            if !no_junction.is_empty() {
                                log::info!("[Post] Admission: {} probable candidates not added, an intron of theirs fails the junction rule", no_junction.len());
                            }
                            let own: Vec<String> = config.admission_by_class.iter()
                                .map(|c| format!("{} {} / {}", c.0, if c.1.is_finite() { c.1.to_string() } else { "never".to_string() }, c.2)).collect();
                            log::info!("[Post] Admission on the read hypergraph: {} candidates scored, {} added (p >= {}), {} dropped (p < {}){} in {:.1}s (fits {:.1}s, features {:.1}s, trees {:.1}s)",
                                       feats.len(), add.len(), config.admission_add, drop.len(), config.admission_drop,
                                       if own.is_empty() { String::new() } else { format!("; classes with their own thresholds: {}", own.join(", ")) },
                                       t_all, t_core, t_feat, t_all - t_core - t_feat);
                        }
                    }
                }
                Some(_) => log::warn!("[Post] Admission: read chains do not match the observed reads; the rules' report is kept"),
            }
        }
    }

    // Gene IDs for each transcript
    let mut result_pairs: Vec<(i64, i64)> = Vec::with_capacity(whitelist.len());
    for &tid in &whitelist {
        let gid = db.tx_map.get(&tid)
            .map(|&idx| db.transcripts[idx].gene_id)
            .unwrap_or(0);
        result_pairs.push((gid, tid));
    }
    result_pairs.sort();

    // Write CSV (no header)
    let mut out = BufWriter::new(File::create(output_file)?);
    for (gid, tid) in &result_pairs {
        writeln!(out, "{},{}", gid, tid)?;
    }

    log::info!("[Post] Filter: {} transcripts written to {}", result_pairs.len(), output_file);
    // Writes the removed one-edit siblings, each with the model it is an edit of and the edit.
    let mut sib = BufWriter::new(File::create(format!("{output_file}.siblings.tsv"))?);
    writeln!(sib, "transcript_ID\tparent_transcript_ID\tedit")?;
    sibling_parent.sort_unstable_by_key(|x| x.0);
    for (t, o, e) in &sibling_parent { writeln!(sib, "{t}\t{o}\t{e:?}")?; }
    Ok(output_file.to_string())
}

// ═══════════════════════════════════════════════════════════════════════════════
// 3. CREATE ABUNDANCE
// ═══════════════════════════════════════════════════════════════════════════════

pub fn run_abundance(
    db: &PostDB,
    output_prefix: &str,
    _genome_build: &str,
    whitelist_file: Option<&str>,
    quant_requires_end_support: bool,
    max_frac_a: f64,
    frac_a_annotated_end_window: i64,
    quant_reassign: bool,
    quant_partial_reads: bool,
    quant_hypergraph: bool,
) -> IsoForgeResult<String> {
    // Load whitelist
    let (mut whitelist, suffix) = if let Some(wl_path) = whitelist_file {
        if Path::new(wl_path).exists() {
            let wl = load_whitelist(wl_path)?;
            (Some(wl), "_filtered")
        } else {
            (None, "")
        }
    } else {
        (None, "")
    };

    // Pre-built annotation dicts
    let gene_id_map = &db.annot.gene_id_map;
    let gene_name_map = &db.annot.gene_name_map;
    let tx_id_map = &db.annot.tx_id_map;
    let tx_name_map = &db.annot.tx_name_map;
    let gene_novelty = &db.annot.gene_novelty;
    let tx_novelty = &db.annot.tx_novelty;
    let ism_prefix = &db.annot.ism_prefix;
    let ism_suffix = &db.annot.ism_suffix;

    // TSS/TES info for length computation
    let tx_ids_needed: FxHashSet<i64> = if let Some(ref wl) = whitelist {
        wl.iter().copied().collect()
    } else {
        db.transcripts.iter().map(|t| t.transcript_id).collect()
    };
    let tss_tes_info = build_tss_tes_info(db, &tx_ids_needed);
    let tx_lengths = compute_transcript_lengths(db, Some(&tss_tes_info));

    // Build abundance pivot: tid → { dataset → count }
    let mut ab_pivot: FxHashMap<i64, FxHashMap<&str, i64>> = FxHashMap::default();
    let gate = InternalPrimingGate::new(db, max_frac_a, frac_a_annotated_end_window);
    if quant_requires_end_support {
        // Reads whose ends joined no TSS/TES cluster are not counted when
        // `quant_requires_end_support` is set; transcript models are unaffected.
        let (mut dropped, mut primed): (i64, i64) = (0, 0);
        for r in &db.observed {
            if !db.in_datasets(&r.dataset) { continue; }
            if let Some(ref wl) = whitelist {
                if !wl.contains(&r.transcript_id) { continue; }
            }
            if r.cluster_supported != Some(true) { dropped += 1; continue; }
            if !gate.eligible(r) { primed += 1; continue; }
            *ab_pivot.entry(r.transcript_id).or_default()
                .entry(&*r.dataset).or_default() += 1;
        }
        log::info!("[Post] Quantification: {} reads without end support and {} internally primed reads excluded from abundance", dropped, primed);
    } else {
        for rec in &db.abundance {
            if !db.datasets.iter().any(|d| d == &rec.dataset) { continue; }
            if let Some(ref wl) = whitelist {
                if !wl.contains(&rec.transcript_id) { continue; }
            }
            *ab_pivot.entry(rec.transcript_id).or_default()
                .entry(rec.dataset.as_str()).or_default() += rec.count;
        }
    }

    // Over the read hypergraph the whole reads of a reported model count for it and the other eligible reads are
    // shared by EM among the reported models they fit. A model left with less than half a read leaves the report and
    // the whitelist is rewritten; when that moves a site of a model that stays, the counts are made again.
    let (mut tss_tes_info, mut tx_lengths) = (tss_tes_info, tx_lengths);
    let hyper: Option<FxHashMap<i64, FxHashMap<String, f64>>> = if quant_hypergraph {
        let quant_path = db.observed_path.parent().map(|p| p.join("quant_reads.parquet"));
        let quant_reads = match quant_path { Some(p) => crate::quant_em::load_quant_reads(&p)?, None => Vec::new() };
        let mut counts = None;
        let mut gone: Vec<(i64, i64)> = Vec::new();
        for pass in 1..=3 {
            let reported: FxHashSet<i64> = db.transcripts.iter().map(|t| t.transcript_id)
                .filter(|t| whitelist.as_ref().map_or(true, |wl| wl.contains(t))).collect();
            let whole = WholeReads::new(db, &tss_tes_info);
            counts = crate::hyper_core::hypergraph_counts(db, &reported, |r| gate.eligible(r), |r| whole.holds(db, r), &quant_reads,
                |q| gate.eligible_fields(&q.chromosome, &q.strand, *q.positions.last().unwrap(), q.fraction_as, q.has_pas_signal),
                &whole.sites, crate::hyper_core::HyperOptions::default())?;
            let (Some(c), Some(wl)) = (counts.as_ref(), whitelist.as_mut()) else { break };
            let unread: FxHashSet<i64> = wl.iter().copied()
                .filter(|t| c.get(t).map_or(0.0, |m| m.values().sum::<f64>()) < 0.5).collect();
            log::info!("[Post] Hypergraph quantification: {} reported models left with less than half a read, removed from the report", unread.len());
            if unread.is_empty() { break; }
            wl.retain(|t| !unread.contains(t));
            gone.extend(unread.iter().map(|&t| (t, ab_pivot.get(&t).map_or(0, |m| m.values().sum::<i64>()))));
            tss_tes_info = build_tss_tes_info(db, wl);
            tx_lengths = compute_transcript_lengths(db, Some(&tss_tes_info));
            let after = hyper_sites(db, &tss_tes_info);
            let moved = wl.iter().filter(|t| after.get(t) != whole.sites.get(t)).count();
            if moved == 0 { break; }
            if pass < 3 { log::info!("[Post] Hypergraph quantification: the removal moves a site of {} models that stay; counting again", moved); }
        }
        if let (false, Some(wl), Some(path)) = (gone.is_empty(), whitelist.as_ref(), whitelist_file) {
            let mut pairs: Vec<(i64, i64)> = wl.iter()
                .map(|&t| (db.tx_map.get(&t).map(|&i| db.transcripts[i].gene_id).unwrap_or(0), t)).collect();
            pairs.sort();
            let mut out = BufWriter::new(File::create(path)?);
            for (gid, tid) in &pairs { writeln!(out, "{},{}", gid, tid)?; }
            gone.sort();
            let mut out = BufWriter::new(File::create(format!("{path}.unread.tsv"))?);
            writeln!(out, "transcript_id\treads")?;
            for (t, n) in gone { writeln!(out, "{}\t{}", t, n)?; }
        }
        counts
    } else { None };
    let shared: Option<FxHashMap<i64, FxHashMap<String, f64>>> = if hyper.is_some() { hyper } else if quant_reassign && quant_requires_end_support {
        let reported: Vec<i64> = db.transcripts.iter().map(|t| t.transcript_id)
            .filter(|t| whitelist.as_ref().map_or(true, |wl| wl.contains(t))).collect();
        let direct: FxHashMap<i64, FxHashMap<String, f64>> = ab_pivot.iter()
            .map(|(&t, m)| (t, m.iter().map(|(d, &c)| (d.to_string(), c as f64)).collect())).collect();
        let quant_path = db.observed_path.parent().map(|p| p.join("quant_reads.parquet"));
        let quant_reads = match quant_path { Some(p) => crate::quant_em::load_quant_reads(&p)?, None => Vec::new() };
        let siblings = whitelist_file.map(|w| crate::quant_em::load_sibling_parents(&format!("{w}.siblings.tsv"))).unwrap_or_default();
        let (counts, st) = crate::quant_em::reassigned_counts(
            db, &reported, &direct, &siblings,
            |r| r.cluster_supported == Some(true) && gate.eligible(r),
            |r| gate.eligible(r),
            &quant_reads,
            |q| gate.eligible_fields(&q.chromosome, &q.strand, *q.positions.last().unwrap(), q.fraction_as, q.has_pas_signal),
            crate::quant_em::QuantOptions { partial_reads: quant_partial_reads, ..Default::default() },
        );
        log::info!("[Post] Quantification: reads fitting reported models: {} of {} of unreported models, {} of {} of reported models without end support, {} of {} quantification-only ({} fit several, shared by EM)",
                   st.rejected_fitted, st.rejected_reads, st.unsupported_fitted, st.unsupported_reads, st.quant_fitted, st.quant_reads, st.ambiguous);
        Some(counts)
    } else {
        None
    };

    // Build output file
    let outpath = format!("{}_isoforge_abundance{}.tsv", output_prefix, suffix);
    let mut out = BufWriter::new(File::create(&outpath)?);

    // Header
    write!(out, "gene_ID\ttranscript_ID\tannot_gene_id\tannot_transcript_id\t\
                 annot_gene_name\tannot_transcript_name\tn_exons\tlength\t\
                 gene_novelty\ttranscript_novelty\tISM_subtype")?;
    for ds in &db.datasets {
        write!(out, "\t{}", ds)?;
    }
    writeln!(out)?;

    // Iterate transcripts
    let mut n_written = 0;
    for tx in &db.transcripts {
        let tid = tx.transcript_id;
        if let Some(ref wl) = whitelist {
            if !wl.contains(&tid) { continue; }
        }

        let gid = tx.gene_id;
        let a_gid = db.format_gene_id(gid, &gene_id_map);
        let a_gname = gene_name_map.get(&gid).cloned().unwrap_or_else(|| a_gid.clone());
        let a_tid = db.format_transcript_id(tid, &tx_id_map);
        let a_tname = tx_name_map.get(&tid).cloned().unwrap_or_else(|| a_tid.clone());

        let g_nov = gene_novelty.get(&gid).cloned().unwrap_or_else(|| "Novel".to_string());
        let t_nov = tx_novelty.get(&tid).cloned().unwrap_or_else(|| "Novel".to_string());

        let is_pfx = ism_prefix.contains(&tid);
        let is_sfx = ism_suffix.contains(&tid);
        let ism_sub = if is_pfx && is_sfx { "Both" }
            else if is_pfx { "Prefix" }
            else if is_sfx { "Suffix" }
            else { "None" };

        let length = tx_lengths.get(&tid).copied().unwrap_or(0);

        write!(out, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
               gid, tid, a_gid, a_tid, a_gname, a_tname,
               tx.n_exons, length, g_nov, t_nov, ism_sub)?;

        if let Some(ref shared) = shared {
            let counts = shared.get(&tid);
            for ds in &db.datasets {
                let c = counts.and_then(|m| m.get(ds.as_str())).copied().unwrap_or(0.0);
                write!(out, "\t{}", format_count(c))?;
            }
        } else {
            let counts = ab_pivot.get(&tid);
            for ds in &db.datasets {
                let c = counts.and_then(|m| m.get(ds.as_str())).copied().unwrap_or(0);
                write!(out, "\t{}", c)?;
            }
        }
        writeln!(out)?;
        n_written += 1;
    }

    log::info!("[Post] Abundance: {} transcripts × {} datasets → {}",
               n_written, db.datasets.len(), outpath);
    Ok(outpath)
}

/// A read count: whole numbers as integers, shared reads to two decimals.
fn format_count(c: f64) -> String {
    let r = c.round();
    if (c - r).abs() < 0.005 { format!("{}", r as i64) } else { format!("{:.2}", c) }
}

// ═══════════════════════════════════════════════════════════════════════════════
// 4. CREATE GTF
// ═══════════════════════════════════════════════════════════════════════════════

pub fn run_gtf(
    db: &PostDB,
    output_prefix: &str,
    _genome_build: &str,
    whitelist_file: Option<&str>,
    observed_only: bool,
    known_end_snap: i64,
    known_end_max_extension: i64,
    known_end_extension_min_fraction: f64,
    known_ends_as_annotated: bool,
) -> IsoForgeResult<String> {
    // Load whitelist
    let whitelist = if let Some(wl_path) = whitelist_file {
        if Path::new(wl_path).exists() {
            Some(load_whitelist(wl_path)?)
        } else { None }
    } else { None };

    // Or use observed_only
    let observed_txs: FxHashSet<i64> = if observed_only && whitelist.is_none() {
        db.observed.iter().map(|r| r.transcript_id).collect()
    } else {
        FxHashSet::default()
    };

    // Pre-built annotation dicts
    let gene_id_map = &db.annot.gene_id_map;
    let gene_name_map = &db.annot.gene_name_map;
    let gene_type_map = &db.annot.gene_type_map;
    let gene_status_map = &db.annot.gene_status_map;
    let tx_id_map = &db.annot.tx_id_map;
    let tx_name_map = &db.annot.tx_name_map;
    let tx_type_map = &db.annot.tx_type_map;
    let tx_status_map = &db.annot.tx_status_map;

    // Build exon coordinate lookup: edge_id → (start, end) where start < end
    let mut exon_coords: FxHashMap<i64, (i64, i64)> = FxHashMap::default();
    for e in &db.edges {
        if e.edge_type == "exon" {
            let p1 = db.loc_map.get(&e.v1).copied().unwrap_or(0);
            let p2 = db.loc_map.get(&e.v2).copied().unwrap_or(0);
            exon_coords.insert(e.edge_id, (p1.min(p2), p1.max(p2)));
        }
    }

    // Group transcripts by gene
    let mut gene_txs: BTreeMap<i64, Vec<&TranscriptRecord>> = BTreeMap::new();
    for tx in &db.transcripts {
        if let Some(ref wl) = whitelist {
            if !wl.contains(&tx.transcript_id) { continue; }
        } else if observed_only && !observed_txs.contains(&tx.transcript_id) {
            continue;
        }
        gene_txs.entry(tx.gene_id).or_default().push(tx);
    }

    // TSS/TES info
    let tx_ids_needed: FxHashSet<i64> = gene_txs.values()
        .flat_map(|txs| txs.iter().map(|t| t.transcript_id))
        .collect();
    let tss_tes_info = build_tss_tes_info(db, &tx_ids_needed);

    let outpath = format!("{}_isoforge.gtf", output_prefix);
    let mut out = BufWriter::new(File::create(&outpath)?);
    let mut lines_written = 0usize;

    for (&gid, txs) in &gene_txs {
        let g_idx = match db.gene_map.get(&gid) {
            Some(&idx) => idx,
            None => continue,
        };
        let g = &db.genes[g_idx];

        let strand = normalize_strand(&g.strand);
        let a_gid = db.format_gene_id(gid, &gene_id_map);
        let a_gname = gene_name_map.get(&gid).cloned().unwrap_or_else(|| a_gid.clone());
        let g_type = gene_type_map.get(&gid).cloned().unwrap_or_else(|| "lncRNA".to_string());
        let g_status = gene_status_map.get(&gid).cloned().unwrap_or_else(|| "NOVEL".to_string());

        // Gene line: gene-level attributes only, no transcript_id.
        writeln!(out, "{}\tIsoForge\tgene\t{}\t{}\t.\t{}\t.\tgene_id \"{}\"; \
                       gene_name \"{}\"; gene_type \"{}\"; gene_status \"{}\"; isoforge_gene \"{}\";",
                 g.chromosome, g.start, g.end, strand, a_gid, a_gname, g_type, g_status, gid)?;
        lines_written += 1;

        for tx in txs {
            let tid = tx.transcript_id;
            let a_tid = db.format_transcript_id(tid, &tx_id_map);
            let a_tname = tx_name_map.get(&tid).cloned().unwrap_or_else(|| a_tid.clone());
            let t_type = tx_type_map.get(&tid).cloned().unwrap_or_else(|| "lncRNA".to_string());
            let t_status = tx_status_map.get(&tid).cloned().unwrap_or_else(|| "NOVEL".to_string());

            // Transcript span from exon coordinates, which are (min, max) on both strands.
            let edge_list = PostDB::get_exon_edge_ids(tx);
            let (mut tx_start, mut tx_end) = {
                let mut gmin = i64::MAX;
                let mut gmax = i64::MIN;
                for &eid in &edge_list {
                    if let Some(&(es, ee)) = exon_coords.get(&eid) {
                        gmin = gmin.min(es);
                        gmax = gmax.max(ee);
                    }
                }
                if gmin == i64::MAX {
                    // Fallback: use stored positions, normalised
                    (tx.start_pos.min(tx.end_pos), tx.start_pos.max(tx.end_pos))
                } else {
                    (gmin, gmax)
                }
            };

            // The transcript line's span covers the exons after terminal clipping.
            let tti = tss_tes_info.get(&tid);
            let n_exons_total = edge_list.len();
            // An annotated transcript keeps its annotated end when the evidence
            // peak lies within `known_end_snap` bp of it; the peak stays in the
            // attributes. Novel models always take their peaks.
            let is_known = db.annot.known_tx_ids.contains(&tid);
            // `outward` is +1 where a larger coordinate lies outside the
            // annotated span and -1 where a smaller one does, so `beyond`
            // measures extension past the annotation in either orientation.
            let keep_annotated = |site: i64, annotated: i64, outward: i64, fraction: f64, feature: Option<&str>| -> bool {
                if !is_known { return false; }
                // A record carrying an annotated ID keeps the annotated coordinates; read
                // evidence stays in tss_site / polyA_site. So does an annotated small RNA, which is counted, not
                // assembled (small_rna.rs).
                if known_ends_as_annotated || db.annot.annotated_forms.contains(&tid) || db.annot.small_rna.contains_key(&tid) { return true; }
                if db.annot.reclaimed.contains(&tid) { return false; }
                keeps_annotated_end(site, annotated, outward, fraction, feature,
                                    known_end_snap, known_end_max_extension, known_end_extension_min_fraction)
            };
            {
                let mut clipped_min = i64::MAX;
                let mut clipped_max = i64::MIN;
                for (idx, &eid) in edge_list.iter().enumerate() {
                    if let Some(&(mut es, mut ee)) = exon_coords.get(&eid) {
                        let exon_num = idx + 1;
                        // Apply the same terminal clipping used for exon lines
                        if let Some(info) = tti {
                            if strand == "+" {
                                if exon_num == 1 {
                                    if let Some(tss) = info.tss_site {
                                        if tss <= ee && !keep_annotated(tss, es, -1, info.tss_fraction, info.tss_feature.as_deref()) { es = tss; }
                                    }
                                }
                                if exon_num == n_exons_total {
                                    if let Some(polya) = info.polya_site {
                                        if polya >= es && !keep_annotated(polya, ee, 1, info.polya_fraction, None) { ee = polya; }
                                    }
                                }
                            } else {
                                if exon_num == 1 {
                                    if let Some(tss) = info.tss_site {
                                        if tss >= es && !keep_annotated(tss, ee, 1, info.tss_fraction, info.tss_feature.as_deref()) { ee = tss; }
                                    }
                                }
                                if exon_num == n_exons_total {
                                    if let Some(polya) = info.polya_site {
                                        if polya <= ee && !keep_annotated(polya, es, -1, info.polya_fraction, None) { es = polya; }
                                    }
                                }
                            }
                        }
                        clipped_min = clipped_min.min(es);
                        clipped_max = clipped_max.max(ee);
                    }
                }
                if clipped_min != i64::MAX {
                    tx_start = clipped_min;
                    tx_end = clipped_max;
                }
            }

            // Transcript attributes
            let mut t_attrs = format!(
                "gene_id \"{}\"; transcript_id \"{}\"; gene_name \"{}\"; transcript_name \"{}\"; \
                 transcript_type \"{}\"; transcript_status \"{}\"; isoforge_gene \"{}\"; isoforge_transcript \"{}\";",
                a_gid, a_tid, a_gname, a_tname, t_type, t_status, gid, tid
            );

            // An end variant names the annotated transcript whose intron chain it carries.
            if let Some(&(parent, _)) = db.annot.variant_parent.get(&tid) {
                if parent != tid {
                    if let Some(name) = db.annot.tx_id_map.get(&parent) {
                        t_attrs.push_str(&format!(" variant_of \"{}\";", name));
                    }
                }
            }

            // Append TSS/polyA cluster attributes
            if let Some(info) = tti {
                if let Some(tss) = info.tss_site {
                    t_attrs.push_str(&format!(" tss_site \"{}\";", tss));
                    t_attrs.push_str(&format!(" tss_Feature \"{}\";",
                                              info.tss_feature.as_deref().unwrap_or("exon")));
                }
                if let Some(polya) = info.polya_site {
                    t_attrs.push_str(&format!(" polyA_site \"{}\";", polya));
                    t_attrs.push_str(&format!(" polyA_Feature \"{}\";",
                                              info.polya_feature.as_deref().unwrap_or("exon")));
                    if let Some(ref sig) = info.polya_signal {
                        t_attrs.push_str(&format!(" polyA_signal \"{}\";", sig));
                    }
                }
            }

            writeln!(out, "{}\tIsoForge\ttranscript\t{}\t{}\t.\t{}\t.\t{}",
                     g.chromosome, tx_start, tx_end, strand, t_attrs)?;
            lines_written += 1;

            // Exon lines
            let n_exons = edge_list.len();
            for (exon_num_0, &eid) in edge_list.iter().enumerate() {
                let exon_num = exon_num_0 + 1;
                let coords = match exon_coords.get(&eid) {
                    Some(c) => *c,
                    None => continue,
                };
                let mut e_start = coords.0;
                let mut e_end = coords.1;

                // Set terminal exon boundaries to TSS/polyA (extend or clip)
                if let Some(info) = tti {
                    if strand == "+" {
                        if exon_num == 1 {
                            if let Some(tss) = info.tss_site {
                                if tss <= e_end && !keep_annotated(tss, coords.0, -1, info.tss_fraction, info.tss_feature.as_deref()) {
                                    e_start = tss;
                                }
                            }
                        }
                        if exon_num == n_exons {
                            if let Some(polya) = info.polya_site {
                                if polya >= e_start && !keep_annotated(polya, coords.1, 1, info.polya_fraction, None) {
                                    e_end = polya;
                                }
                            }
                        }
                    } else {
                        if exon_num == 1 {
                            if let Some(tss) = info.tss_site {
                                if tss >= e_start && !keep_annotated(tss, coords.1, 1, info.tss_fraction, info.tss_feature.as_deref()) {
                                    e_end = tss;
                                }
                            }
                        }
                        if exon_num == n_exons {
                            if let Some(polya) = info.polya_site {
                                if polya <= e_end && !keep_annotated(polya, coords.0, -1, info.polya_fraction, None) {
                                    e_start = polya;
                                }
                            }
                        }
                    }
                }

                writeln!(out, "{}\tIsoForge\texon\t{}\t{}\t.\t{}\t.\t\
                              gene_id \"{}\"; transcript_id \"{}\"; gene_name \"{}\"; \
                              transcript_name \"{}\"; exon_number \"{}\"; \
                              isoforge_gene \"{}\"; isoforge_transcript \"{}\"; isoforge_exon \"{}\";",
                         g.chromosome, e_start, e_end, strand,
                         a_gid, a_tid, a_gname, a_tname, exon_num,
                         gid, tid, eid)?;
                lines_written += 1;
            }
        }
    }

    log::info!("[Post] GTF: {} lines ({} genes, {} transcripts) → {}",
               lines_written, gene_txs.len(),
               gene_txs.values().map(|v| v.len()).sum::<usize>(), outpath);
    Ok(outpath)
}

// ═══════════════════════════════════════════════════════════════════════════════
// 5. GET READ ANNOTATIONS
// ═══════════════════════════════════════════════════════════════════════════════

/// Pre-formatted fields for one observed read, computed in parallel and
/// written sequentially.
struct ReadAnnotRow {
    a_gid: String,
    a_tid: String,
    a_gname: String,
    a_tname: String,
    g_nov: String,
    t_nov: String,
    ism_sub: String,
    read_start: i32,
    read_end: i32,
    strand: String,
    n_exons: i32,
}

/// Write per-read annotations.
///
/// With `split_unsupported`, reads that are not validated and reads excluded
/// by the internal-priming gate (see `InternalPrimingGate`) are written to
/// `_read_annot_unsupported.tsv` instead of the main table. Both files carry
/// identical columns.
pub fn run_read_annot(
    db: &PostDB,
    output_prefix: &str,
    genome_build: &str,
    whitelist_file: Option<&str>,
    split_unsupported: bool,
    max_frac_a: f64,
    frac_a_annotated_end_window: i64,
    emit_tsv: bool,
) -> IsoForgeResult<(String, String)> {
    let gate = InternalPrimingGate::new(db, max_frac_a, frac_a_annotated_end_window);
    // Annotation dicts
    let gene_id_map = &db.annot.gene_id_map;
    let gene_name_map = &db.annot.gene_name_map;
    let tx_id_map = &db.annot.tx_id_map;
    let tx_name_map = &db.annot.tx_name_map;
    let gene_novelty = &db.annot.gene_novelty;
    let tx_novelty = &db.annot.tx_novelty;
    let ism_prefix = &db.annot.ism_prefix;
    let ism_suffix = &db.annot.ism_suffix;

    // Gene strand map
    let gene_strand: FxHashMap<i64, String> = db.genes.iter()
        .map(|g| (g.gene_id, normalize_strand(&g.strand)))
        .collect();

    // Transcript n_exons map
    let tx_nexons: FxHashMap<i64, i32> = db.transcripts.iter()
        .map(|t| (t.transcript_id, t.n_exons))
        .collect();

    // Dataset set for O(1) membership
    let dataset_set: FxHashSet<&str> = db.datasets.iter().map(|s| s.as_str()).collect();

    // Load whitelist
    let whitelist_tids = if let Some(wl_path) = whitelist_file {
        if Path::new(wl_path).exists() {
            Some(load_whitelist(wl_path)?)
        } else { None }
    } else { None };
    // The validated reads of a reported model are its whole reads, the ones the counts keep with it; a read of a
    // model outside the report keeps its end validation.
    let whole = whitelist_tids.as_ref().map(|wl| WholeReads::new(db, &build_tss_tes_info(db, wl)));
    let validated = |r: &ObservedRead| match (&whole, &whitelist_tids) {
        // an annotated small RNA's validated reads are the reads it holds, as the counts take them (small_rna.rs)
        _ if db.annot.small_rna.contains_key(&r.transcript_id) =>
            db.annot.small_rna[&r.transcript_id].holds(r.start_pos, r.end_pos) && gate.eligible(r),
        (Some(w), Some(wl)) if wl.contains(&r.transcript_id) => w.holds(db, r) && gate.eligible(r),
        _ => r.validated(),
    };

    // ── Write TSV (whitelist-filtered) ────────────────────────────────────
    let outpath_tsv = format!("{}_isoforge_read_annot.tsv", output_prefix);
    let mut out_tsv = if emit_tsv { Some(BufWriter::new(File::create(&outpath_tsv)?)) } else { None };

    // Header
    if let Some(ref mut out_tsv) = out_tsv {
    writeln!(out_tsv, "read_name\tdataset\tgenome_build\tchromosome\t\
                       read_start\tread_end\tstrand\tn_exons\tread_length\t\
                       gene_ID\ttranscript_ID\t\
                       annot_gene_id\tannot_transcript_id\t\
                       annot_gene_name\tannot_transcript_name\t\
                       gene_novelty\ttranscript_novelty\tISM_subtype\t\
                       fraction_As\tcustom_label\tallelic_label\t\
                       start_support\tend_support\t\
                       tss_site\ttss_feature\tpolya_site\tpolya_feature\tpolya_signal\t\
                       cluster_supported\tpolya_length")?;
    }

    let mut tsv_count = 0usize;

    // Second table, opened only when the split is requested.
    let outpath_unsup = format!("{}_isoforge_read_annot_unsupported.tsv", output_prefix);
    let mut out_unsup = if split_unsupported && emit_tsv {
        let mut w = BufWriter::new(File::create(&outpath_unsup)?);
        writeln!(w, "read_name\tdataset\tgenome_build\tchromosome\t\
                     read_start\tread_end\tstrand\tn_exons\tread_length\t\
                     gene_ID\ttranscript_ID\t\
                     annot_gene_id\tannot_transcript_id\t\
                     annot_gene_name\tannot_transcript_name\t\
                     gene_novelty\ttranscript_novelty\tISM_subtype\t\
                     fraction_As\tcustom_label\tallelic_label\t\
                     start_support\tend_support\t\
                     tss_site\ttss_feature\tpolya_site\tpolya_feature\tpolya_signal\t\
                     cluster_supported\tpolya_length")?;
        Some(w)
    } else {
        None
    };
    let mut unsup_count = 0usize;

    // ── Parquet schema & writer ──────────────────────────────────────────
    let outpath_pq = format!("{}_isoforge_read_annot.parquet", output_prefix);
    let pq_schema = Arc::new(Schema::new(vec![
        Field::new("read_name", DataType::Utf8, false),
        Field::new("dataset", DataType::Utf8, false),
        Field::new("genome_build", DataType::Utf8, false),
        Field::new("chromosome", DataType::Utf8, false),
        Field::new("read_start", DataType::Int32, false),
        Field::new("read_end", DataType::Int32, false),
        Field::new("strand", DataType::Utf8, false),
        Field::new("n_exons", DataType::Int32, false),
        Field::new("read_length", DataType::Int32, false),
        Field::new("gene_ID", DataType::Int32, false),
        Field::new("transcript_ID", DataType::Int32, false),
        Field::new("annot_gene_id", DataType::Utf8, false),
        Field::new("annot_transcript_id", DataType::Utf8, false),
        Field::new("annot_gene_name", DataType::Utf8, false),
        Field::new("annot_transcript_name", DataType::Utf8, false),
        Field::new("gene_novelty", DataType::Utf8, false),
        Field::new("transcript_novelty", DataType::Utf8, false),
        Field::new("ISM_subtype", DataType::Utf8, false),
        Field::new("fraction_As", DataType::Float32, true),
        Field::new("custom_label", DataType::Utf8, true),
        Field::new("allelic_label", DataType::Utf8, true),
        Field::new("start_support", DataType::Utf8, true),
        Field::new("end_support", DataType::Utf8, true),
        Field::new("tss_site", DataType::Int64, true),
        Field::new("tss_feature", DataType::Utf8, true),
        Field::new("polya_site", DataType::Int64, true),
        Field::new("polya_feature", DataType::Utf8, true),
        Field::new("polya_signal", DataType::Utf8, true),
        Field::new("cluster_supported", DataType::Boolean, true),
        Field::new("polya_length", DataType::Int64, true),
    ]));

    let pq_file = File::create(&outpath_pq)?;
    let mut pq_writer = ArrowWriter::try_new(pq_file, pq_schema.clone(), Some(parquet_io::writer_props()))?;

    // ── Chunked parallel processing ──────────────────────────────────────
    // Reads are processed in CHUNK_SIZE batches; within each chunk TSV lines and
    // Arrow batches are built in parallel and written sequentially.
    const CHUNK_SIZE: usize = 1 << 18;

    let mut chunks_written = 0usize;

    // Format one read's annotation fields; returns
    // (Option<tsv_line_bytes>, parquet_row_data).
    #[inline]
    fn format_read<'a>(
        r: &'a ObservedRead,
        read_name: &str,
        genome_build: &str,
        gene_id_map: &FxHashMap<i64, String>,
        gene_name_map: &FxHashMap<i64, String>,
        tx_id_map: &FxHashMap<i64, String>,
        tx_name_map: &FxHashMap<i64, String>,
        gene_novelty: &FxHashMap<i64, String>,
        tx_novelty: &FxHashMap<i64, String>,
        ism_prefix: &FxHashSet<i64>,
        ism_suffix: &FxHashSet<i64>,
        gene_strand: &FxHashMap<i64, String>,
        tx_nexons: &FxHashMap<i64, i32>,
        loc_map: &FxHashMap<i64, i64>,
        id_prefix: &str,
        n_places: usize,
        whitelist_tids: &Option<FxHashSet<i64>>,
        emit_tsv: bool,
        validated: bool,
    ) -> (Option<Vec<u8>>, ReadAnnotRow) {
        let gid = r.gene_id;
        let tid = r.transcript_id;

        let a_gid = match gene_id_map.get(&gid) {
            Some(name) => name.clone(),
            None => format!("{}G{:0width$}", id_prefix, gid, width = n_places),
        };
        let a_gname = gene_name_map.get(&gid).cloned().unwrap_or_else(|| a_gid.clone());
        let a_tid = match tx_id_map.get(&tid) {
            Some(name) => name.clone(),
            None => format!("{}T{:0width$}", id_prefix, tid, width = n_places),
        };
        let a_tname = tx_name_map.get(&tid).cloned().unwrap_or_else(|| a_tid.clone());

        let g_nov = gene_novelty.get(&gid).map(|s| s.as_str()).unwrap_or("Novel");
        let t_nov = tx_novelty.get(&tid).map(|s| s.as_str()).unwrap_or("Novel");

        let is_pfx = ism_prefix.contains(&tid);
        let is_sfx = ism_suffix.contains(&tid);
        let ism_sub: &str = if is_pfx && is_sfx { "Both" }
            else if is_pfx { "Prefix" }
            else if is_sfx { "Suffix" }
            else { "None" };

        let sv_pos = loc_map.get(&r.start_vertex).copied().unwrap_or(0);
        let ev_pos = loc_map.get(&r.end_vertex).copied().unwrap_or(0);
        let strand = gene_strand.get(&gid).map(|s| s.as_str()).unwrap_or("+");

        // Use actual BAM alignment positions when available; fall back to vertex positions.
        // read_start is the 5' end and read_end the 3' end, so read_start > read_end on the - strand.
        let (a, b) = if r.start_pos != 0 || r.end_pos != 0 {
            (r.start_pos, r.end_pos)
        } else {
            (sv_pos, ev_pos)
        };
        let (read_start, read_end) = if strand == "-" {
            (a.max(b), a.min(b))
        } else {
            (a.min(b), a.max(b))
        };
        let n_exons = tx_nexons.get(&tid).copied().unwrap_or(r.n_exons);

        // TSV line (only if whitelisted)
        let write_tsv = match whitelist_tids {
            Some(wl) => wl.contains(&tid),
            None => true,
        };
        let tsv_line = if write_tsv && emit_tsv {
            let frac_as_str = match r.fraction_as {
                Some(f) => format!("{}", f),
                None => String::new(),
            };
            let custom = r.custom_label.as_deref().unwrap_or("");
            let allelic = r.allelic_label.as_deref().unwrap_or("");
            let start_sup = r.start_support.as_deref().unwrap_or("");
            let end_sup = r.end_support.as_deref().unwrap_or("");
            let tss_s = r.tss_site.map(|v| v.to_string()).unwrap_or_default();
            let tss_f = r.tss_feature.as_deref().unwrap_or("");
            let polya_s = r.polya_site.map(|v| v.to_string()).unwrap_or_default();
            let polya_f = r.polya_feature.as_deref().unwrap_or("");
            let polya_sig = r.polya_signal.as_deref().unwrap_or("");
            let cs = r.cluster_supported.map(|_| if validated { "True" } else { "False" }).unwrap_or("");
            let polya_len = r.polya_length.or(r.polya_softclip).map(|v| v.to_string()).unwrap_or_default();

            let mut buf = Vec::with_capacity(300);
            use std::io::Write as _;
            let _ = writeln!(buf, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                     read_name, r.dataset, genome_build, r.chromosome,
                     read_start, read_end, strand, n_exons, r.read_length,
                     gid, tid, a_gid, a_tid, a_gname, a_tname,
                     g_nov, t_nov, ism_sub, frac_as_str,
                     custom, allelic, start_sup, end_sup,
                     tss_s, tss_f, polya_s, polya_f, polya_sig, cs, polya_len);
            Some(buf)
        } else {
            None
        };

        let row = ReadAnnotRow {
            a_gid,
            a_tid,
            a_gname,
            a_tname,
            g_nov: g_nov.to_string(),
            t_nov: t_nov.to_string(),
            ism_sub: ism_sub.to_string(),
            read_start: read_start as i32,
            read_end: read_end as i32,
            strand: strand.to_string(),
            n_exons,
        };

        (tsv_line, row)
    }

    // Read names come back from observed.parquet, whose rows db.observed holds in order.
    let mut names = if db.observed.is_empty() { None } else { Some(ReadNames::open(&db.observed_path)?) };
    for chunk in db.observed.chunks(CHUNK_SIZE) {
        let chunk_names = names.as_mut().map(|n| n.take(chunk.len())).transpose()?.unwrap_or_default();
        // Filter to known datasets and process in parallel
        let results: Vec<(usize, Option<Vec<u8>>, ReadAnnotRow, &ObservedRead, bool)> = chunk.par_iter()
            .enumerate()
            .filter(|(_, r)| dataset_set.contains(&*r.dataset))
            .map(|(idx, r)| {
                let ok = validated(r);
                let (tsv, row) = format_read(
                    r, &chunk_names[idx], genome_build,
                    gene_id_map, gene_name_map, tx_id_map, tx_name_map,
                    gene_novelty, tx_novelty, ism_prefix, ism_suffix,
                    &gene_strand, &tx_nexons, &db.loc_map,
                    &db.id_prefix, db.n_places, &whitelist_tids, emit_tsv, ok,
                );
                (idx, tsv, row, r, ok)
            })
            .collect();

        // Build Arrow arrays + write TSV — sequential (maintains row order)
        let n = results.len();
        let mut pq_read_names = StringBuilder::with_capacity(n, n * 40);
        let mut pq_datasets = StringBuilder::with_capacity(n, n * 20);
        let mut pq_genome_builds = StringBuilder::with_capacity(n, n * 8);
        let mut pq_chroms = StringBuilder::with_capacity(n, n * 6);
        let mut pq_read_starts = Int32Builder::with_capacity(n);
        let mut pq_read_ends = Int32Builder::with_capacity(n);
        let mut pq_strands = StringBuilder::with_capacity(n, n);
        let mut pq_n_exons = Int32Builder::with_capacity(n);
        let mut pq_read_lengths = Int32Builder::with_capacity(n);
        let mut pq_gene_ids = Int32Builder::with_capacity(n);
        let mut pq_transcript_ids = Int32Builder::with_capacity(n);
        let mut pq_annot_gene_ids = StringBuilder::with_capacity(n, n * 20);
        let mut pq_annot_tx_ids = StringBuilder::with_capacity(n, n * 20);
        let mut pq_annot_gene_names = StringBuilder::with_capacity(n, n * 20);
        let mut pq_annot_tx_names = StringBuilder::with_capacity(n, n * 20);
        let mut pq_gene_novelties = StringBuilder::with_capacity(n, n * 10);
        let mut pq_tx_novelties = StringBuilder::with_capacity(n, n * 10);
        let mut pq_ism_subtypes = StringBuilder::with_capacity(n, n * 6);
        let mut pq_fraction_as = Float32Builder::with_capacity(n);
        let mut pq_custom_labels = StringBuilder::new();
        let mut pq_allelic_labels = StringBuilder::new();
        let mut pq_start_supports = StringBuilder::new();
        let mut pq_end_supports = StringBuilder::new();
        let mut pq_tss_sites = Int64Builder::with_capacity(n);
        let mut pq_tss_features = StringBuilder::new();
        let mut pq_polya_sites = Int64Builder::with_capacity(n);
        let mut pq_polya_features = StringBuilder::new();
        let mut pq_polya_signals = StringBuilder::new();
        let mut pq_cluster_supporteds = BooleanBuilder::with_capacity(n);
        let mut pq_polya_lengths = Int64Builder::with_capacity(n);

        for (idx, tsv_line, row, r, ok) in &results {
            // With splitting, reads that failed end validation or are internally
            // primed go to the unsupported table.
            if let Some(ref line) = tsv_line {
                match out_unsup {
                    Some(ref mut w) if !*ok || !gate.eligible(r) => {
                        w.write_all(line)?;
                        unsup_count += 1;
                    }
                    _ => {
                        if let Some(ref mut w) = out_tsv {
                            w.write_all(line)?;
                            tsv_count += 1;
                        }
                    }
                }
            }

            // Build Parquet arrays
            pq_read_names.append_value(&chunk_names[*idx]);
            pq_datasets.append_value(&r.dataset);
            pq_genome_builds.append_value(genome_build);
            pq_chroms.append_value(&r.chromosome);
            pq_read_starts.append_value(row.read_start);
            pq_read_ends.append_value(row.read_end);
            pq_strands.append_value(&row.strand);
            pq_n_exons.append_value(row.n_exons);
            pq_read_lengths.append_value(r.read_length as i32);
            pq_gene_ids.append_value(r.gene_id as i32);
            pq_transcript_ids.append_value(r.transcript_id as i32);
            pq_annot_gene_ids.append_value(&row.a_gid);
            pq_annot_tx_ids.append_value(&row.a_tid);
            pq_annot_gene_names.append_value(&row.a_gname);
            pq_annot_tx_names.append_value(&row.a_tname);
            pq_gene_novelties.append_value(&row.g_nov);
            pq_tx_novelties.append_value(&row.t_nov);
            pq_ism_subtypes.append_value(&row.ism_sub);
            match r.fraction_as {
                Some(f) => pq_fraction_as.append_value(f as f32),
                None => pq_fraction_as.append_null(),
            }
            match &r.custom_label {
                Some(l) => pq_custom_labels.append_value(l),
                None => pq_custom_labels.append_null(),
            }
            match &r.allelic_label {
                Some(l) => pq_allelic_labels.append_value(l),
                None => pq_allelic_labels.append_null(),
            }
            match &r.start_support {
                Some(s) => pq_start_supports.append_value(s),
                None => pq_start_supports.append_null(),
            }
            match &r.end_support {
                Some(s) => pq_end_supports.append_value(s),
                None => pq_end_supports.append_null(),
            }
            match r.tss_site {
                Some(v) => pq_tss_sites.append_value(v),
                None => pq_tss_sites.append_null(),
            }
            match &r.tss_feature {
                Some(f) => pq_tss_features.append_value(f),
                None => pq_tss_features.append_null(),
            }
            match r.polya_site {
                Some(v) => pq_polya_sites.append_value(v),
                None => pq_polya_sites.append_null(),
            }
            match &r.polya_feature {
                Some(f) => pq_polya_features.append_value(f),
                None => pq_polya_features.append_null(),
            }
            match &r.polya_signal {
                Some(s) => pq_polya_signals.append_value(s),
                None => pq_polya_signals.append_null(),
            }
            match r.cluster_supported {
                Some(_) => pq_cluster_supporteds.append_value(*ok),
                None => pq_cluster_supporteds.append_null(),
            }
            match r.polya_length.or(r.polya_softclip) {
                Some(v) => pq_polya_lengths.append_value(v),
                None => pq_polya_lengths.append_null(),
            }
        }

        // Flush chunk as one Parquet row group
        if !results.is_empty() {
            let batch = RecordBatch::try_new(
                pq_schema.clone(),
                vec![
                    Arc::new(pq_read_names.finish()),
                    Arc::new(pq_datasets.finish()),
                    Arc::new(pq_genome_builds.finish()),
                    Arc::new(pq_chroms.finish()),
                    Arc::new(pq_read_starts.finish()),
                    Arc::new(pq_read_ends.finish()),
                    Arc::new(pq_strands.finish()),
                    Arc::new(pq_n_exons.finish()),
                    Arc::new(pq_read_lengths.finish()),
                    Arc::new(pq_gene_ids.finish()),
                    Arc::new(pq_transcript_ids.finish()),
                    Arc::new(pq_annot_gene_ids.finish()),
                    Arc::new(pq_annot_tx_ids.finish()),
                    Arc::new(pq_annot_gene_names.finish()),
                    Arc::new(pq_annot_tx_names.finish()),
                    Arc::new(pq_gene_novelties.finish()),
                    Arc::new(pq_tx_novelties.finish()),
                    Arc::new(pq_ism_subtypes.finish()),
                    Arc::new(pq_fraction_as.finish()),
                    Arc::new(pq_custom_labels.finish()),
                    Arc::new(pq_allelic_labels.finish()),
                    Arc::new(pq_start_supports.finish()),
                    Arc::new(pq_end_supports.finish()),
                    Arc::new(pq_tss_sites.finish()),
                    Arc::new(pq_tss_features.finish()),
                    Arc::new(pq_polya_sites.finish()),
                    Arc::new(pq_polya_features.finish()),
                    Arc::new(pq_polya_signals.finish()),
                    Arc::new(pq_cluster_supporteds.finish()),
                    Arc::new(pq_polya_lengths.finish()),
                ],
            )?;
            pq_writer.write(&batch)?;
            chunks_written += 1;
        }
    }

    drop(out_tsv);
    if let Some(w) = out_unsup {
        drop(w);
        log::info!("[Post] End-unvalidated reads written separately: {} reads -> {}",
                   unsup_count, outpath_unsup);
    }

    pq_writer.close()?;

    log::info!("[Post] Read annotations Parquet: {} chunks written", chunks_written);

    let tsv_size = if emit_tsv { fs::metadata(&outpath_tsv).map(|m| m.len()).unwrap_or(0) } else { 0 };
    let pq_size = fs::metadata(&outpath_pq).map(|m| m.len()).unwrap_or(0);
    if emit_tsv {
        log::info!("[Post] Read annotations: {} reads, TSV={:.1}MB ({}), Parquet={:.1}MB",
                   db.observed.len(), tsv_size as f64 / 1e6, tsv_count,
                   pq_size as f64 / 1e6);
    } else {
        log::info!("[Post] Read annotations: {} reads, Parquet={:.1}MB (TSV not written)",
                   db.observed.len(), pq_size as f64 / 1e6);
    }

    Ok((if emit_tsv { outpath_tsv } else { String::new() }, outpath_pq))
}

// ═══════════════════════════════════════════════════════════════════════════════
// 6. PRUNE CLUSTER BEDS BY WHITELIST
// ═══════════════════════════════════════════════════════════════════════════════

pub fn run_prune_clusters(db: &PostDB, whitelist_file: &str) -> IsoForgeResult<()> {
    let wl = load_whitelist(whitelist_file)?;

    // Collect surviving cluster peak positions from whitelist reads
    let mut surviving_polya: FxHashSet<i64> = FxHashSet::default();
    let mut surviving_tss: FxHashSet<i64> = FxHashSet::default();

    // (a fragment keeps no cluster at the end it lacks)
    for r in &db.observed {
        if wl.contains(&r.transcript_id) {
            if let Some(polya) = r.polya_site.filter(|_| r.fragment & MISSING_3P == 0) {
                surviving_polya.insert(polya);
            }
            if let Some(tss) = r.tss_site.filter(|_| r.fragment & MISSING_5P == 0) {
                surviving_tss.insert(tss);
            }
        }
    }

    // Locate cluster directory
    let _results_dir = Path::new(&db.counters.data.get("_results_dir")
        .and_then(|v| v.as_str())
        .unwrap_or(""))
        .to_path_buf();

    // Infer from whitelist path
    let wl_path = Path::new(whitelist_file);
    let wl_dir = wl_path.parent().unwrap_or(Path::new("."));

    // Try to find the cluster directory based on naming convention
    let db_entries: Vec<_> = fs::read_dir(wl_dir)
        .ok()
        .map(|rd| rd.filter_map(|e| e.ok()).collect())
        .unwrap_or_default();

    let mut cluster_dir: Option<PathBuf> = None;
    for entry in &db_entries {
        let name = entry.file_name();
        let name_str = name.to_string_lossy();
        if name_str.ends_with("_TSS_TES_Clusters") && entry.path().is_dir() {
            cluster_dir = Some(entry.path());
            break;
        }
    }

    let cluster_dir = match cluster_dir {
        Some(d) => d,
        None => {
            log::info!("[Post] Cluster directory not found, skipping pruning");
            return Ok(());
        }
    };

    // Filter BED files
    let tes_pruned = filter_bed_file(&cluster_dir.join("tes_clusters.bed"), &surviving_polya)?;
    let tss_pruned = filter_bed_file(&cluster_dir.join("tss_clusters.bed"), &surviving_tss)?;

    // Filter per-sample BED files
    let mut _sample_tes_pruned = 0usize;
    let mut _sample_tss_pruned = 0usize;
    let sample_dir = cluster_dir.join("sample");
    if sample_dir.is_dir() {
        for entry in fs::read_dir(&sample_dir)? {
            let entry = entry?;
            let name = entry.file_name();
            let name_str = name.to_string_lossy();
            if name_str.ends_with("_tes.bed") {
                _sample_tes_pruned += filter_bed_file(&entry.path(), &surviving_polya)?;
            } else if name_str.ends_with("_tss.bed") {
                _sample_tss_pruned += filter_bed_file(&entry.path(), &surviving_tss)?;
            }
        }
    }

    // Filter per-treatment BED files
    let treat_dir = cluster_dir.join("treatment");
    if treat_dir.is_dir() {
        for entry in fs::read_dir(&treat_dir)? {
            let entry = entry?;
            let name = entry.file_name();
            let name_str = name.to_string_lossy();
            if name_str.ends_with("_tes.bed") {
                filter_bed_file(&entry.path(), &surviving_polya)?;
            } else if name_str.ends_with("_tss.bed") {
                filter_bed_file(&entry.path(), &surviving_tss)?;
            }
        }
    }

    // Filter cluster Parquet
    let pq_path = cluster_dir.join("tss_tes_clusters.parquet");
    if pq_path.exists() {
        filter_cluster_parquet(&pq_path, &surviving_polya, &surviving_tss)?;
    }

    log::info!("[Post] Cluster BED pruning: TES={}, TSS={} pruned from main BED. \
                Surviving: {} TES + {} TSS peaks",
               tes_pruned, tss_pruned, surviving_polya.len(), surviving_tss.len());

    Ok(())
}


// ═══════════════════════════════════════════════════════════════════════════════
// HELPER FUNCTIONS
// ═══════════════════════════════════════════════════════════════════════════════

/// Spliced-read tail rate from which starts are established by tailed reads.
const FULL_LENGTH_START_RATE: f64 = 0.1;
/// Spliced-read tail rate above which starts are not established by tailed reads.
const FULL_LENGTH_START_MAX_RATE: f64 = 0.9;
const ANNOTATED_START_WINDOW: i64 = 25;
/// A known transcript keeps its annotated start when its read-derived 5' end lies more than this far inside the
/// annotated span, off any promoter feature.
const KNOWN_END_INSIDE_TOLERANCE: i64 = 100;
/// Reads within ANNOTATED_FORM_WINDOW bp of an annotated end its record does not report
/// that keep the annotated form as a model of its own.
const ANNOTATED_FORM_READS: usize = 3;
const ANNOTATED_FORM_WINDOW: i64 = 50;
/// Distance from every annotated 3' end beyond which a novel model's 3' end is a new
/// poly(A) site.
const NEW_END_WINDOW: i64 = 50;

/// TSS/TES info per transcript
#[derive(Debug, Clone)]
pub struct TssTesInfo {
    pub tss_site: Option<i64>,
    pub tss_feature: Option<String>,
    pub polya_site: Option<i64>,
    pub polya_feature: Option<String>,
    pub polya_signal: Option<String>,
    /// Fraction of this transcript's end-assigned reads carried by the
    /// selected site.
    pub tss_fraction: f64,
    pub polya_fraction: f64,
    /// The site is an end cluster its reads joined (not the commonest read end of a model none of whose reads
    /// joined one).
    pub tss_clustered: bool,
    pub polya_clustered: bool,
}

/// Whether an annotated record keeps its annotated end rather than the read-derived
/// `site`. `outward` is +1 where a larger coordinate lies outside the annotated span and
/// -1 where a smaller one does. The annotated end stands when the peak lies within
/// `snap` bp of it; when a 5' peak (`feature` given) lies more than
/// KNOWN_END_INSIDE_TOLERANCE bp inside the annotated start on a feature that is not a
/// promoter; and when a peak extends the transcript by more than `max_extension` bp on
/// less than `min_fraction` of its end-assigned reads.
fn keeps_annotated_end(site: i64, annotated: i64, outward: i64, fraction: f64, feature: Option<&str>,
                       snap: i64, max_extension: i64, min_fraction: f64) -> bool {
    if snap > 0 && (site - annotated).abs() <= snap { return true; }
    if let Some(f) = feature {
        let inside = (annotated - site) * outward;
        if inside > KNOWN_END_INSIDE_TOLERANCE && !matches!(f, "five_prime_utr" | "intergenic") {
            return true;
        }
    }
    if max_extension > 0 {
        let beyond = (site - annotated) * outward;
        if beyond > max_extension && fraction < min_fraction {
            return true;
        }
    }
    false
}

/// Keep the annotated form of an annotated transcript whose record takes a read-derived
/// end more than ANNOTATED_FORM_WINDOW bp from the annotation. When at least
/// ANNOTATED_FORM_READS of its end-validated reads end within that window of the annotated
/// end, and some end-validated reads end elsewhere, the reads there keep the annotated ID
/// and are written at the annotated ends, and the others become an end variant under a
/// new ID with the read-derived ends and the annotated record's standing in the filters.
/// Both forms are then reported. Returns the number of records split.
pub fn split_annotated_forms(db: &mut PostDB, snap: i64, max_extension: i64, min_fraction: f64) -> usize {
    // an annotated small RNA keeps all its reads under its annotated ID (small_rna.rs)
    let known: FxHashSet<i64> = db.observed.iter().map(|r| r.transcript_id)
        .filter(|t| db.annot.known_tx_ids.contains(t) && !db.annot.small_rna.contains_key(t)).collect();
    let info = build_tss_tes_info(db, &known);
    let mut splits: Vec<(i64, i64, i64, bool, bool)> = Vec::new(); // (tid, tss, tes, 5' differs, 3' differs)
    let mut read_ends: FxHashMap<i64, (Option<i64>, Option<i64>)> = FxHashMap::default();
    for &tid in &known {
        let (Some(&i), Some(ti)) = (db.tx_map.get(&tid), info.get(&tid)) else { continue };
        let t = &db.transcripts[i];
        let Some(strand) = db.gene_map.get(&t.gene_id).map(|&g| normalize_strand(&db.genes[g].strand)) else { continue };
        let (Some(&a), Some(&b)) = (db.loc_map.get(&t.start_vertex), db.loc_map.get(&t.end_vertex)) else { continue };
        let (lo, hi) = (a.min(b), a.max(b));
        let minus = strand == "-";
        let (tss, tes) = if minus { (hi, lo) } else { (lo, hi) };
        // The two forms are distinct only when the read-derived end lies beyond the window
        // in which reads count as reaching the annotated end.
        let differs5 = ti.tss_site.is_some_and(|x| (x - tss).abs() > ANNOTATED_FORM_WINDOW
            && !keeps_annotated_end(x, tss, if minus { 1 } else { -1 }, ti.tss_fraction, ti.tss_feature.as_deref(), snap, max_extension, min_fraction));
        let differs3 = ti.polya_site.is_some_and(|x| (x - tes).abs() > ANNOTATED_FORM_WINDOW
            && !keeps_annotated_end(x, tes, if minus { -1 } else { 1 }, ti.polya_fraction, None, snap, max_extension, min_fraction));
        if differs5 || differs3 {
            splits.push((tid, tss, tes, differs5, differs3));
            read_ends.insert(tid, (ti.tss_site.filter(|_| differs5), ti.polya_site.filter(|_| differs3)));
        }
    }
    let at_annotated = |r: &ObservedRead, tss: i64, tes: i64, d5: bool, d3: bool| {
        (d5 && (r.start_pos - tss).abs() <= ANNOTATED_FORM_WINDOW) || (d3 && (r.end_pos - tes).abs() <= ANNOTATED_FORM_WINDOW)
    };
    // End-validated reads at the annotated end, and elsewhere: both forms must keep reads
    // that count.
    let mut support: FxHashMap<i64, (usize, usize)> = FxHashMap::default();
    let plan: FxHashMap<i64, (i64, i64, bool, bool)> = splits.iter().map(|&(t, a, b, c, d)| (t, (a, b, c, d))).collect();
    for r in db.observed.iter().filter(|r| r.cluster_supported == Some(true)) {
        if let Some(&(tss, tes, d5, d3)) = plan.get(&r.transcript_id) {
            let e = support.entry(r.transcript_id).or_default();
            if at_annotated(r, tss, tes, d5, d3) { e.0 += 1 } else { e.1 += 1 }
        }
    }
    let mut next = db.transcripts.iter().map(|t| t.transcript_id).max().unwrap_or(0).max(db.counters.transcripts);
    let mut variant_of: FxHashMap<i64, i64> = FxHashMap::default();
    for &(tid, ..) in &splits {
        let (at, elsewhere) = support.get(&tid).copied().unwrap_or((0, 0));
        if at < ANNOTATED_FORM_READS || elsewhere == 0 { continue; }
        next += 1;
        let mut rec = db.transcripts[db.tx_map[&tid]].clone();
        rec.transcript_id = next;
        rec.transcript_name = None;
        db.tx_map.insert(next, db.transcripts.len());
        db.transcripts.push(rec);
        // the variant keeps the record's chain, and ends where its reads place them (`variant_ends`)
        if let Some(c) = db.chains.get(&tid).cloned() {
            let n = c.len();
            if n >= 2 {
                let (plus, (r5, r3)) = (c[0] <= c[n - 1], read_ends.get(&tid).copied().unwrap_or((None, None)));
                db.annot.variant_ends.insert(next, (r5.filter(|&x| if plus { x <= c[1] } else { x >= c[1] }),
                                                    r3.filter(|&x| if plus { x >= c[n - 2] } else { x <= c[n - 2] })));
            }
            db.chains.insert(next, c);
        }
        let a = &mut db.annot;
        if let Some(ty) = a.tx_type_map.get(&tid).cloned() { a.tx_type_map.insert(next, ty); }
        a.tx_status_map.insert(next, "NOVEL".to_string());
        a.tx_novelty.insert(next, "Novel".to_string());
        a.known_tx_ids.insert(next);
        a.annotated_forms.insert(tid);
        a.variant_parent.insert(next, (tid, true));
        variant_of.insert(tid, next);
    }
    db.counters.transcripts = next;
    for r in db.observed.iter_mut() {
        if let (Some(&v), Some(&(tss, tes, d5, d3))) = (variant_of.get(&r.transcript_id), plan.get(&r.transcript_id)) {
            if !at_annotated(r, tss, tes, d5, d3) { r.transcript_id = v; }
        }
    }
    // Recount the abundance of every split record from its reads.
    let touched: FxHashSet<i64> = variant_of.iter().flat_map(|(&a, &b)| [a, b]).collect();
    db.abundance.retain(|a| !touched.contains(&a.transcript_id));
    let mut counts: BTreeMap<(i64, i64, String), i64> = BTreeMap::new();
    for r in db.observed.iter().filter(|r| touched.contains(&r.transcript_id)) {
        *counts.entry((r.gene_id, r.transcript_id, r.dataset.to_string())).or_default() += 1;
    }
    db.abundance.extend(counts.into_iter().map(|((gene_id, transcript_id, dataset), count)| AbundanceRecord { gene_id, transcript_id, dataset, count }));
    // demotes degraded starts again after reads moved between the two forms
    let n_demoted = demote_unclustered_degraded_starts(&mut db.observed);
    log::info!("[Post] After the split, reads starting more than {} bp inside every start of their model no longer end-validated: {}",
               UNCLUSTERED_START_WINDOW, n_demoted);
    log::info!("[Post] Annotated forms kept beside read-derived ends (>= {} reads within {} bp of the annotated end): {} records split",
               ANNOTATED_FORM_READS, ANNOTATED_FORM_WINDOW, variant_of.len());
    variant_of.len()
}

/// Report under its annotated ID a transcript that passed the filters only as an end
/// variant of its own annotated transcript: when the annotated record does not pass, the
/// largest passing variant marked compatible with the annotated ends takes the ID. Its
/// reads and those of the annotated record are counted under it, and it is written where
/// the variant was, at its peaks. An incompatible variant keeps its novel ID. Rewrites the
/// whitelist and returns the number of records reclaimed.
pub fn reclaim_annotated_records(db: &mut PostDB, whitelist_file: &str) -> IsoForgeResult<usize> {
    if db.annot.variant_parent.is_empty() {
        log::info!("[Post] No end-variant table in this run database: annotated IDs are not reclaimed");
        return Ok(0);
    }
    let mut wl = load_whitelist(whitelist_file)?;
    let mut reads: FxHashMap<i64, usize> = FxHashMap::default();
    for r in &db.observed { *reads.entry(r.transcript_id).or_default() += 1; }
    let mut best: FxHashMap<i64, (usize, i64)> = FxHashMap::default(); // annotated -> (reads, variant)
    for &tid in &wl {
        let Some(&(k, true)) = db.annot.variant_parent.get(&tid) else { continue };
        if wl.contains(&k) || !db.tx_map.contains_key(&k) { continue; }
        let n = reads.get(&tid).copied().unwrap_or(0);
        let e = best.entry(k).or_insert((n, tid));
        if n > e.0 || (n == e.0 && tid < e.1) { *e = (n, tid); }
    }
    let fold: FxHashMap<i64, (i64, i64)> = best.iter()
        .filter_map(|(&k, &(_, v))| db.tx_map.get(&k).map(|&i| (v, (k, db.transcripts[i].gene_id))))
        .collect();
    let touched: FxHashSet<i64> = fold.iter().flat_map(|(&v, &(k, _))| [v, k]).collect();
    for r in db.observed.iter_mut() {
        if let Some(&(k, gene)) = fold.get(&r.transcript_id) {
            r.transcript_id = k;
            r.gene_id = gene;
        }
    }
    db.abundance.retain(|a| !touched.contains(&a.transcript_id));
    let mut counts: BTreeMap<(i64, i64, String), i64> = BTreeMap::new();
    for r in db.observed.iter().filter(|r| touched.contains(&r.transcript_id)) {
        *counts.entry((r.gene_id, r.transcript_id, r.dataset.to_string())).or_default() += 1;
    }
    db.abundance.extend(counts.into_iter().map(|((gene_id, transcript_id, dataset), count)| AbundanceRecord { gene_id, transcript_id, dataset, count }));
    let n_demoted = demote_unclustered_degraded_starts(&mut db.observed);
    log::info!("[Post] After the reclaim, reads starting more than {} bp inside every start of their model no longer end-validated: {}",
               UNCLUSTERED_START_WINDOW, n_demoted);
    for (&v, &(k, _)) in &fold {
        wl.remove(&v);
        wl.insert(k);
        db.annot.reclaimed.insert(k);
    }
    let mut pairs: Vec<(i64, i64)> = wl.iter()
        .map(|&t| (db.tx_map.get(&t).map(|&i| db.transcripts[i].gene_id).unwrap_or(0), t))
        .collect();
    pairs.sort();
    let mut out = BufWriter::new(File::create(whitelist_file)?);
    for (gid, tid) in &pairs {
        writeln!(out, "{},{}", gid, tid)?;
    }
    let incompatible = wl.iter().filter(|t| db.annot.variant_parent.get(t).is_some_and(|&(k, ok)| !ok && !wl.contains(&k))).count();
    log::info!("[Post] Annotated records reported in place of their only passing, compatible end variant: {}; \
                end variants reported as isoforms of their own, their annotated transcript not passing: {}",
               fold.len(), incompatible);
    Ok(fold.len())
}

/// The end cluster a transcript is placed at, with the share of the counted reads it
/// holds: the cluster most of its supported reads carry, or without any, most of all
/// its reads; on a tie the outermost, the higher coordinate when `outermost_high`.
fn end_site(supported: Option<&FxHashMap<i64, i64>>, all: Option<&FxHashMap<i64, i64>>, outermost_high: bool)
    -> Option<(i64, f64)>
{
    let sites = supported.filter(|m| !m.is_empty()).or(all)?;
    let total: i64 = sites.values().sum();
    sites.iter()
        .max_by(|a, b| a.1.cmp(b.1).then(if outermost_high { a.0.cmp(b.0) } else { b.0.cmp(a.0) }))
        .map(|(&site, &count)| (site, count as f64 / total.max(1) as f64))
}

/// The outermost cluster holding at least a tenth of a transcript's reads, leaving out the
/// clusters in `taken`, with its share of the counted reads; the higher coordinate is the
/// outermost when `outermost_high`.
fn outer_site(all: Option<&FxHashMap<i64, i64>>, taken: &FxHashSet<i64>, outermost_high: bool) -> Option<(i64, f64)> {
    let sites = all?;
    let total: i64 = sites.values().sum();
    let min_support = ((total as f64) * 0.1).max(1.0) as i64;
    let qualified = sites.iter().filter(|(s, &c)| c >= min_support && !taken.contains(s));
    let best = if outermost_high { qualified.max_by_key(|(s, _)| **s) } else { qualified.min_by_key(|(s, _)| **s) };
    best.map(|(&s, &c)| (s, c as f64 / total.max(1) as f64))
}

/// The TSS and poly(A) site the reads of each model place, as the read hypergraph takes them: the end clusters its
/// reads joined, in the model's terminal exon or beyond its end. No site: the commonest read end of a model without
/// a cluster, and a 5' cluster more than KNOWN_END_INSIDE_TOLERANCE bp inside an annotated start off any promoter.
fn hyper_sites(db: &PostDB, info: &FxHashMap<i64, TssTesInfo>) -> crate::hyper_core::Sites {
    info.iter().map(|(&t, i)| {
        let (mut tss, mut tes) = (i.tss_site.filter(|_| i.tss_clustered), i.polya_site.filter(|_| i.polya_clustered));
        if let (Some(c), Some((five, _))) = (db.chains.get(&t).filter(|c| c.len() >= 2), db.chain_ends(t)) {
            let (n, plus) = (c.len(), c[0] <= c[c.len() - 1]);
            // a site lies in the model's terminal exon or beyond its end
            tss = tss.filter(|&s| if plus { s <= c[1] } else { s >= c[1] });
            tes = tes.filter(|&s| if plus { s >= c[n - 2] } else { s <= c[n - 2] });
            if let Some(site) = tss {
                let inside = if plus { site - five } else { five - site };
                let promoter = matches!(i.tss_feature.as_deref(), Some("five_prime_utr" | "intergenic"));
                if db.annot.known_tx_ids.contains(&t) && !db.annot.reclaimed.contains(&t) && inside > KNOWN_END_INSIDE_TOLERANCE && !promoter {
                    tss = None;
                }
            }
        }
        (t, (tss, tes))
    }).collect()
}

/// Which reads are whole reads of the reported model they are counted for: end-validated, carrying every exon of the
/// model, and with both ends at the model's (`at_model_ends`).
struct WholeReads {
    sites: crate::hyper_core::Sites,
}

impl WholeReads {
    fn new(db: &PostDB, info: &FxHashMap<i64, TssTesInfo>) -> Self { Self { sites: hyper_sites(db, info) } }

    fn holds(&self, db: &PostDB, r: &ObservedRead) -> bool {
        r.validated() && db.chain_ends(r.transcript_id).is_some_and(|ends| at_model_ends(r, ends, self.sites.get(&r.transcript_id)))
    }
}

/// Whether a read's ends are its model's: each within `END_WINDOW` bp of the end of the model's chain (`ends`, 5' and
/// 3') or of the site the model's reads place there.
fn at_model_ends(r: &ObservedRead, (c5, c3): (i64, i64), sites: Option<&(Option<i64>, Option<i64>)>) -> bool {
    let w = crate::hyper_core::END_WINDOW as i64;
    let (s5, s3) = sites.copied().unwrap_or((None, None));
    let near = |x: i64, end: i64, site: Option<i64>| (x - end).abs() <= w || site.is_some_and(|s| (x - s).abs() <= w);
    let (lo, hi) = (r.start_pos.min(r.end_pos), r.start_pos.max(r.end_pos));
    let (five, three) = if &*r.strand == "-" { (hi, lo) } else { (lo, hi) };
    near(five, c5, s5) && near(three, c3, s3)
}

/// Build per-transcript TSS/TES sites from observed reads.
fn build_tss_tes_info(db: &PostDB, tx_ids: &FxHashSet<i64>) -> FxHashMap<i64, TssTesInfo> {
    let mut result: FxHashMap<i64, TssTesInfo> = FxHashMap::default();

    // Collect TSS counts: tid → (tss_site → count)
    let mut tss_counts: FxHashMap<i64, FxHashMap<i64, i64>> = FxHashMap::default();
    let mut tss_features: FxHashMap<(i64, i64), String> = FxHashMap::default(); // (tid, site) → feature
    let mut polya_counts: FxHashMap<i64, FxHashMap<i64, i64>> = FxHashMap::default();
    let mut polya_features: FxHashMap<(i64, i64), String> = FxHashMap::default();
    let mut polya_signals: FxHashMap<(i64, i64), String> = FxHashMap::default();
    // The same counts over supported reads only, which place a transcript's ends
    // whenever it has any.
    let mut tss_supported: FxHashMap<i64, FxHashMap<i64, i64>> = FxHashMap::default();
    let mut polya_supported: FxHashMap<i64, FxHashMap<i64, i64>> = FxHashMap::default();

    // Fallback: collect actual BAM alignment positions (start_pos = 5' end, end_pos = 3' end)
    let mut read_5p_counts: FxHashMap<i64, FxHashMap<i64, i64>> = FxHashMap::default();
    let mut read_3p_counts: FxHashMap<i64, FxHashMap<i64, i64>> = FxHashMap::default();

    for r in &db.observed {
        if !tx_ids.contains(&r.transcript_id) { continue; }

        let supported = r.cluster_supported == Some(true);
        // a fragment does not count at the end it lacks
        let (has5, has3) = (r.fragment & MISSING_5P == 0, r.fragment & MISSING_3P == 0);
        if let Some(tss) = r.tss_site.filter(|_| has5) {
            *tss_counts.entry(r.transcript_id).or_default()
                .entry(tss).or_default() += 1;
            if supported {
                *tss_supported.entry(r.transcript_id).or_default().entry(tss).or_default() += 1;
            }
            if let Some(ref f) = r.tss_feature {
                tss_features.entry((r.transcript_id, tss)).or_insert_with(|| f.to_string());
            }
        }

        if let Some(polya) = r.polya_site.filter(|_| has3) {
            *polya_counts.entry(r.transcript_id).or_default()
                .entry(polya).or_default() += 1;
            if supported {
                *polya_supported.entry(r.transcript_id).or_default().entry(polya).or_default() += 1;
            }
            if let Some(ref f) = r.polya_feature {
                polya_features.entry((r.transcript_id, polya)).or_insert_with(|| f.to_string());
            }
            if let Some(ref s) = r.polya_signal {
                polya_signals.entry((r.transcript_id, polya)).or_insert_with(|| s.to_string());
            }
        }

        // Actual BAM positions for the fallback (zero positions skipped)
        if r.start_pos != 0 && has5 {
            *read_5p_counts.entry(r.transcript_id).or_default()
                .entry(r.start_pos).or_default() += 1;
        }
        if r.end_pos != 0 && has3 {
            *read_3p_counts.entry(r.transcript_id).or_default()
                .entry(r.end_pos).or_default() += 1;
        }
    }

    // Build strand lookup for outermost-cluster selection
    let gene_strand: FxHashMap<i64, bool> = db.genes.iter()
        .map(|g| (g.gene_id, normalize_strand(&g.strand) == "+"))
        .collect();
    let tx_gene: FxHashMap<i64, i64> = db.transcripts.iter()
        .map(|tx| (tx.transcript_id, tx.gene_id))
        .collect();

    // A novel model, an end variant among them, sits at the clusters most of its supported
    // reads carry. An annotated record takes the outermost cluster holding at least a tenth
    // of all its reads, tailed or not, but not a cluster that one of its own end variants
    // stands for where that variant differs from it.
    let all_tids: FxHashSet<i64> = tss_counts.keys()
        .chain(polya_counts.keys())
        .chain(read_5p_counts.keys())
        .chain(read_3p_counts.keys())
        .copied().collect();
    let strand_of = |tid: i64| tx_gene.get(&tid).and_then(|gid| gene_strand.get(gid)).copied().unwrap_or(true);
    let modal: FxHashMap<i64, (Option<i64>, Option<i64>)> = all_tids.iter()
        .filter(|t| tx_ids.contains(t))
        .map(|&t| {
            let plus = strand_of(t);
            (t, (end_site(tss_supported.get(&t), tss_counts.get(&t), !plus).map(|x| x.0),
                 end_site(polya_supported.get(&t), polya_counts.get(&t), plus).map(|x| x.0)))
        })
        .collect();
    let mut taken: FxHashMap<i64, (FxHashSet<i64>, FxHashSet<i64>)> = FxHashMap::default();
    for (&v, &(parent, _)) in &db.annot.variant_parent {
        if v == parent { continue; }
        let (Some(&(v5, v3)), Some(&(p5, p3))) = (modal.get(&v), modal.get(&parent)) else { continue };
        let e = taken.entry(parent).or_default();
        if let Some(x) = v5.filter(|&x| Some(x) != p5) { e.0.insert(x); }
        if let Some(x) = v3.filter(|&x| Some(x) != p3) { e.1.insert(x); }
    }
    let none = (FxHashSet::default(), FxHashSet::default());
    for tid in all_tids {
        if !tx_ids.contains(&tid) { continue; }

        let is_plus = strand_of(tid);
        let known = db.annot.known_tx_ids.contains(&tid);
        let (taken5, taken3) = taken.get(&tid).unwrap_or(&none);

        let mut info = TssTesInfo {
            tss_site: None, tss_feature: None,
            polya_site: None, polya_feature: None, polya_signal: None,
            tss_fraction: 0.0, polya_fraction: 0.0,
            tss_clustered: false, polya_clustered: false,
        };

        // Outermost 5': min coord on +strand, max coord on -strand
        let tss = if known { outer_site(tss_counts.get(&tid), taken5, !is_plus) } else { None }
            .or_else(|| end_site(tss_supported.get(&tid), tss_counts.get(&tid), !is_plus));
        if let Some((site, fraction)) = tss {
            info.tss_site = Some(site);
            info.tss_clustered = true;
            info.tss_fraction = fraction;
            info.tss_feature = tss_features.get(&(tid, site)).cloned()
                .or(Some("exon".to_string()));
        }

        // Fallback: when no TSS cluster, use modal actual BAM 5' position
        if info.tss_site.is_none() {
            if let Some(positions) = read_5p_counts.get(&tid) {
                if let Some((&max_pos, &max_count)) = positions.iter().max_by_key(|(_, &c)| c) {
                    info.tss_site = Some(max_pos);
                    info.tss_fraction = max_count as f64 / positions.values().sum::<i64>().max(1) as f64;
                    info.tss_feature = Some("exon".to_string());
                }
            }
        }

        // Outermost 3': max coord on +strand, min coord on -strand
        let polya = if known { outer_site(polya_counts.get(&tid), taken3, is_plus) } else { None }
            .or_else(|| end_site(polya_supported.get(&tid), polya_counts.get(&tid), is_plus));
        if let Some((site, fraction)) = polya {
            info.polya_site = Some(site);
            info.polya_clustered = true;
            info.polya_fraction = fraction;
            info.polya_feature = polya_features.get(&(tid, site)).cloned()
                .or(Some("exon".to_string()));
            info.polya_signal = polya_signals.get(&(tid, site)).cloned();
        }

        // Fallback: when no polyA cluster, use modal actual BAM 3' position
        if info.polya_site.is_none() {
            if let Some(positions) = read_3p_counts.get(&tid) {
                if let Some((&max_pos, &max_count)) = positions.iter().max_by_key(|(_, &c)| c) {
                    info.polya_site = Some(max_pos);
                    info.polya_fraction = max_count as f64 / positions.values().sum::<i64>().max(1) as f64;
                    info.polya_feature = Some("exon".to_string());
                }
            }
        }

        result.insert(tid, info);
    }

    log::info!("[Post] TSS/TES info: {} transcripts with cluster assignments", result.len());
    result
}

/// Compute transcript lengths from exon edge coordinates, with TSS/TES clipping.
fn compute_transcript_lengths(db: &PostDB, tss_tes_info: Option<&FxHashMap<i64, TssTesInfo>>)
    -> FxHashMap<i64, i64>
{
    // Build exon edge lengths
    let mut edge_len: FxHashMap<i64, i64> = FxHashMap::default();
    let mut edge_coords: FxHashMap<i64, (i64, i64)> = FxHashMap::default();
    for e in &db.edges {
        if e.edge_type == "exon" {
            let p1 = db.loc_map.get(&e.v1).copied().unwrap_or(0);
            let p2 = db.loc_map.get(&e.v2).copied().unwrap_or(0);
            let len = (p2 - p1).abs() + 1;
            edge_len.insert(e.edge_id, len);
            edge_coords.insert(e.edge_id, (p1.min(p2), p1.max(p2)));
        }
    }

    // Gene strand map
    let gene_strand: FxHashMap<i64, String> = db.genes.iter()
        .map(|g| (g.gene_id, normalize_strand(&g.strand)))
        .collect();

    let mut lengths: FxHashMap<i64, i64> = FxHashMap::default();

    for tx in &db.transcripts {
        let se = tx.start_exon;
        let ee = tx.end_exon;

        let mut total = edge_len.get(&se).copied().unwrap_or(0);

        if let Some(ref jn) = tx.jn_path {
            if !jn.is_empty() {
                let parts: Vec<&str> = jn.split(',').collect();
                for j in (1..parts.len()).step_by(2) {
                    if let Ok(eid) = parts[j].parse::<i64>() {
                        total += edge_len.get(&eid).copied().unwrap_or(0);
                    }
                }
            }
        }

        if ee != se {
            total += edge_len.get(&ee).copied().unwrap_or(0);
        }

        // Clip terminal exons to actual TSS/polyA sites
        if let Some(tti_map) = tss_tes_info {
            if let Some(tti) = tti_map.get(&tx.transcript_id) {
                let is_plus = gene_strand.get(&tx.gene_id)
                    .map(|s| s == "+").unwrap_or(true);

                let edge_list = PostDB::get_exon_edge_ids(tx);
                let first_eid = edge_list[0];
                let last_eid = *edge_list.last().unwrap();

                let first_coords = edge_coords.get(&first_eid);
                let last_coords = edge_coords.get(&last_eid);

                let mut adjustment = 0i64;

                if is_plus {
                    if let (Some(tss), Some(&(e_s, e_e))) = (tti.tss_site, first_coords) {
                        if e_s < tss && tss <= e_e {
                            adjustment += tss - e_s;
                        }
                    }
                    if let (Some(polya), Some(&(e_s, e_e))) = (tti.polya_site, last_coords) {
                        if e_s <= polya && polya < e_e {
                            adjustment += e_e - polya;
                        }
                    }
                } else {
                    if let (Some(tss), Some(&(e_s, e_e))) = (tti.tss_site, first_coords) {
                        if e_s <= tss && tss < e_e {
                            adjustment += e_e - tss;
                        }
                    }
                    if let (Some(polya), Some(&(e_s, e_e))) = (tti.polya_site, last_coords) {
                        if e_s < polya && polya <= e_e {
                            adjustment += polya - e_s;
                        }
                    }
                }

                if adjustment > 0 {
                    total = (total - adjustment).max(1);
                }
            }
        }

        lengths.insert(tx.transcript_id, total);
    }

    lengths
}

/// Load whitelist file (CSV: gene_ID, transcript_ID — no header)
pub fn load_whitelist(path: &str) -> IsoForgeResult<FxHashSet<i64>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut wl = FxHashSet::default();
    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split(',').collect();
        if parts.len() >= 2 {
            if let Ok(tid) = parts[1].trim().parse::<i64>() {
                wl.insert(tid);
            }
        }
    }
    Ok(wl)
}

/// Normalize strand string: "Plus" → "+", "Minus" → "-"
fn normalize_strand(s: &str) -> String {
    match s {
        "Plus" | "+" => "+".to_string(),
        "Minus" | "-" => "-".to_string(),
        _ => "+".to_string(),
    }
}

/// Filter a BED file, keeping only lines whose column-4 peak is in surviving_peaks.
fn filter_bed_file(path: &Path, surviving_peaks: &FxHashSet<i64>) -> IsoForgeResult<usize> {
    if !path.exists() {
        return Ok(0);
    }
    let content = fs::read_to_string(path)?;
    let lines: Vec<&str> = content.lines().collect();
    let mut kept = Vec::new();
    let mut pruned = 0usize;

    for line in &lines {
        if line.is_empty() { continue; }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() > 3 {
            if let Ok(peak) = parts[3].parse::<i64>() {
                if surviving_peaks.contains(&peak) {
                    kept.push(*line);
                } else {
                    pruned += 1;
                }
                continue;
            }
        }
        kept.push(*line); // keep headers/malformed lines
    }

    let mut out = File::create(path)?;
    for line in &kept {
        writeln!(out, "{}", line)?;
    }

    Ok(pruned)
}

/// Filter cluster Parquet by surviving peak positions.
fn filter_cluster_parquet(path: &Path, surviving_polya: &FxHashSet<i64>, surviving_tss: &FxHashSet<i64>)
    -> IsoForgeResult<()>
{
    let batches = parquet_io::read_parquet_batches(path)?;
    if batches.is_empty() { return Ok(()); }

    let schema = batches[0].schema();

    // Find column indices
    let type_idx = schema.index_of("cluster_type").ok();
    let peak_idx = schema.index_of("peak_site").ok();

    if type_idx.is_none() || peak_idx.is_none() {
        return Ok(()); // schema doesn't match, skip
    }
    let type_idx = type_idx.unwrap();
    let peak_idx = peak_idx.unwrap();

    // Filter rows and rebuild batches
    let mut kept_indices = Vec::new();
    let mut row_offset = 0usize;

    for batch in &batches {
        let types = batch.column(type_idx).as_any().downcast_ref::<StringArray>().unwrap();
        let peaks = batch.column(peak_idx).as_any().downcast_ref::<Int64Array>().unwrap();

        for i in 0..batch.num_rows() {
            let ct = types.value(i);
            let peak = peaks.value(i);
            let keep = match ct {
                "TES" => surviving_polya.contains(&peak),
                "TSS" => surviving_tss.contains(&peak),
                _ => true,
            };
            if keep {
                kept_indices.push((row_offset, i));
            }
        }
        row_offset += 1;
    }

    // Concatenate all batches and apply a boolean filter
    let all_batch = arrow::compute::concat_batches(&schema, &batches)?;
    let n = all_batch.num_rows();
    let mut filter_bools = BooleanBuilder::with_capacity(n);

    let types = all_batch.column(type_idx).as_any().downcast_ref::<StringArray>().unwrap();
    let peaks = all_batch.column(peak_idx).as_any().downcast_ref::<Int64Array>().unwrap();

    for i in 0..n {
        let ct = types.value(i);
        let peak = peaks.value(i);
        let keep = match ct {
            "TES" => surviving_polya.contains(&peak),
            "TSS" => surviving_tss.contains(&peak),
            _ => true,
        };
        filter_bools.append_value(keep);
    }

    let filter_arr = filter_bools.finish();
    let filtered = arrow::compute::filter_record_batch(&all_batch, &filter_arr)?;

    // Write back
    let file = File::create(path)?;
    let mut writer = ArrowWriter::try_new(file, schema, Some(parquet_io::writer_props()))?;
    writer.write(&filtered)?;
    writer.close()?;

    Ok(())
}

#[cfg(test)]
mod start_peak_tests {
    use super::*;

    #[test]
    fn start_peaks_come_from_dense_starts_not_scattered_ones() {
        // 60 reads over 12 bp at 1000, 20 at a second start 300 bp downstream, 20 scattered degraded starts
        let mut x: Vec<i64> = (0..60).map(|i| 1000 + i % 12).collect();
        x.extend((0..20).map(|i| 1300 + i % 5));
        x.extend((0..20).map(|i| 1500 + 97 * i));
        let p = read_start_peaks(&mut x);
        assert_eq!(p.len(), 2);
        assert!((1000..1012).contains(&p[0]) && (1300..1305).contains(&p[1]), "{:?}", p);
        // scattered starts only: no peak
        let mut y: Vec<i64> = (0..40).map(|i| 1000 + 53 * i).collect();
        assert!(read_start_peaks(&mut y).is_empty());
    }
}

#[cfg(test)]
mod fragment_tests {
    use super::*;

    #[test]
    fn a_read_with_fewer_exons_than_its_model_is_a_fragment_of_it() {
        let tx = |id: i64, n: i32| TranscriptRecord { transcript_id: id, gene_id: 1, transcript_name: None, start_exon: 0, jn_path: None,
            end_exon: 0, start_vertex: 0, end_vertex: 0, n_exons: n, chromosome: "chr1".into(), start_pos: 0, end_pos: 0 };
        let transcripts = vec![tx(1, 3), tx(2, 3)];
        let tx_map: FxHashMap<i64, usize> = [(1, 0), (2, 1)].into_iter().collect();
        // model 1 on the plus strand (exons 100-200, 300-400, 500-600), model 2 the same on the minus strand
        let chains: FxHashMap<i64, Vec<i64>> = [(1, vec![100, 200, 300, 400, 500, 600]), (2, vec![600, 500, 400, 300, 200, 100])].into_iter().collect();
        let read = |tid: i64, strand: &str, n: i32, five: i64, three: i64| ObservedRead {
            transcript_id: tid, strand: Arc::from(strand), n_exons: n, start_pos: five, end_pos: three, cluster_supported: Some(true), ..Default::default() };
        let mut reads = vec![
            read(1, "+", 3, 110, 590),   // every exon
            read(1, "+", 1, 520, 600),   // unspliced, inside the last exon
            read(1, "+", 2, 310, 598),   // the last two exons
            read(1, "+", 2, 100, 390),   // the first two
            read(1, "+", 1, 320, 380),   // inside the middle exon
            read(2, "-", 1, 180, 100),   // minus strand, inside the last exon
            read(2, "-", 2, 600, 310),   // minus strand, the first two exons
        ];
        assert_eq!(mark_fragments(&mut reads, &transcripts, &tx_map, &chains), 6);
        let marks: Vec<u8> = reads.iter().map(|r| r.fragment).collect();
        assert_eq!(marks, vec![0, MISSING_5P, MISSING_5P, MISSING_3P, MISSING_5P | MISSING_3P, MISSING_5P, MISSING_3P]);
        assert!(reads[0].validated() && reads.iter().skip(1).all(|r| !r.validated()));
    }

    #[test]
    fn a_whole_read_ends_at_the_chains_ends_or_at_the_sites_the_model_takes() {
        let read = |five: i64, three: i64| ObservedRead { strand: Arc::from("+"), start_pos: five, end_pos: three, ..Default::default() };
        let chain = (1000, 2300);
        // at the chain's own ends
        assert!(at_model_ends(&read(1010, 2290), chain, None));
        // 69 bp inside the start: only with the site the model's reads place there
        assert!(!at_model_ends(&read(1069, 2290), chain, None));
        assert!(at_model_ends(&read(1069, 2290), chain, Some(&(Some(1069), None))));
        // a proximal poly(A) site of the model, and a read that reaches neither 3' end
        assert!(at_model_ends(&read(1000, 2100), chain, Some(&(None, Some(2110)))));
        assert!(!at_model_ends(&read(1000, 2100), chain, Some(&(Some(1069), None))));
        // minus strand: the 5' end is the higher coordinate, whichever way the row holds it
        let minus = |a: i64, b: i64| ObservedRead { strand: Arc::from("-"), start_pos: a, end_pos: b, ..Default::default() };
        assert!(at_model_ends(&minus(2295, 1005), (2300, 1000), None));
        assert!(at_model_ends(&minus(1005, 2295), (2300, 1000), None));
        assert!(!at_model_ends(&minus(2100, 1005), (2300, 1000), None));
    }
}

#[cfg(test)]
mod evidence_tests {

    #[test]
    fn annotated_end_rule_is_symmetric_across_strands() {
        // Within the snap distance the annotated end stands, on either side.
        assert!(keeps_annotated_end(1010, 1000, 1, 1.0, None, 25, 25, 0.5));
        assert!(keeps_annotated_end(990, 1000, -1, 1.0, None, 25, 25, 0.5));
        // A proximal 3' peak 300 bp inside the annotated end is taken.
        assert!(!keeps_annotated_end(700, 1000, 1, 1.0, None, 25, 25, 0.5));
        assert!(!keeps_annotated_end(1300, 1000, -1, 1.0, None, 25, 25, 0.5));
        // An extension carried by a minority of reads is not.
        assert!(keeps_annotated_end(1200, 1000, 1, 0.2, None, 25, 25, 0.5));
        assert!(!keeps_annotated_end(1200, 1000, 1, 0.8, None, 25, 25, 0.5));
        // A 5' peak deep inside the start, off a promoter feature, is not taken.
        assert!(keeps_annotated_end(1300, 1000, -1, 1.0, Some("exon"), 25, 25, 0.5));
        assert!(!keeps_annotated_end(1300, 1000, -1, 1.0, Some("five_prime_utr"), 25, 25, 0.5));
    }

    use super::*;
    #[test]
    fn concentration_is_bounded_not_transitive() {
        assert_eq!(compact_support(&[],10),0);
        assert_eq!(compact_support(&[1],10),1);
        assert_eq!(compact_support(&[100,102,104,106,109],10),5);
        assert_eq!(compact_support(&[100,110,120,130,140],10),2);
        assert_eq!(compact_support(&[100,100,100,200,210],10),3);
    }
}

#[cfg(test)]
mod junction_evidence_tests {
    use super::*;
    use crate::types::Strand;
    #[test]
    fn independent_junction_evidence_requires_three_unique_reads_and_correct_strand() {
        let mut info = crate::genome_ref::SjInfo { strand:1, motif:0, annotated:false, unique_reads:1, multi_reads:100 };
        assert!(!supported_short_read_junction(&info, Strand::Plus));
        info.unique_reads = 3;
        assert!(supported_short_read_junction(&info, Strand::Plus));
        assert!(!supported_short_read_junction(&info, Strand::Minus));
        info.strand = 0;
        assert!(!supported_short_read_junction(&info, Strand::Plus));
        info.strand = 2;
        assert!(supported_short_read_junction(&info, Strand::Minus));
    }
}

#[cfg(test)]
mod end_site_tests {
    use super::*;
    fn counts(v: &[(i64, i64)]) -> FxHashMap<i64, i64> { v.iter().copied().collect() }

    #[test]
    fn a_model_sits_at_the_cluster_its_supported_reads_carry() {
        // 17 supported reads at 195; 7 unsupported reads merged in from 514, further out.
        let sup = counts(&[(195, 17)]);
        let all = counts(&[(195, 25), (514, 7)]);
        assert_eq!(end_site(Some(&sup), Some(&all), true).map(|s| s.0), Some(195));
        assert_eq!(end_site(Some(&sup), Some(&all), true).map(|s| s.1), Some(1.0));
        // Without supported reads, the cluster most of all reads carry.
        assert_eq!(end_site(None, Some(&all), true).map(|s| s.0), Some(195));
        assert_eq!(end_site(Some(&counts(&[])), Some(&all), false).map(|s| s.0), Some(195));
        // A tie goes to the outermost cluster on either strand.
        let tie = counts(&[(100, 3), (200, 3)]);
        assert_eq!(end_site(Some(&tie), None, true).map(|s| s.0), Some(200));
        assert_eq!(end_site(Some(&tie), None, false).map(|s| s.0), Some(100));
        assert_eq!(end_site(None, None, true), None);
    }

    #[test]
    fn an_annotated_record_takes_its_outermost_cluster_but_not_its_variants() {
        // 96 reads at 629, 41 at the annotated end 583, 18 at 955 (minus strand: lower is outermost).
        let all = counts(&[(629, 96), (583, 41), (955, 18)]);
        assert_eq!(outer_site(Some(&all), &FxHashSet::default(), false).map(|s| s.0), Some(583));
        // A cluster holding less than a tenth of the reads does not count.
        let thin = counts(&[(629, 96), (583, 5)]);
        assert_eq!(outer_site(Some(&thin), &FxHashSet::default(), false).map(|s| s.0), Some(629));
        // The cluster an end variant of the record stands for is left to that variant.
        let gng = counts(&[(195, 25), (514, 7)]);
        let taken: FxHashSet<i64> = [514].into_iter().collect();
        assert_eq!(outer_site(Some(&gng), &FxHashSet::default(), true).map(|s| s.0), Some(514));
        assert_eq!(outer_site(Some(&gng), &taken, true).map(|s| s.0), Some(195));
        let only: FxHashSet<i64> = [195, 514].into_iter().collect();
        assert_eq!(outer_site(Some(&gng), &only, true), None);
    }
}

#[cfg(test)]
mod end_variant_table_tests {
    use super::*;
    #[test]
    fn end_variant_table_round_trips_and_is_optional() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("end_variants.parquet");
        assert!(load_end_variants(&path).unwrap().is_empty());
        parquet_io::write_end_variants(&path, &[(900, 12, true), (901, 12, false)]).unwrap();
        let v = load_end_variants(&path).unwrap();
        assert_eq!(v.get(&900), Some(&(12, true)));
        assert_eq!(v.get(&901), Some(&(12, false)));
        assert_eq!(v.len(), 2);
    }
}

#[cfg(test)]
mod start_model_tests {
    use super::*;

    #[test]
    fn a_pile_up_at_the_start_is_a_spike_and_a_smear_is_not() {
        // 20 reads within a base of the start, against 20 spread evenly over the window.
        let mut spike: Vec<i64> = (0..20).map(|i| 1000 + (i % 3) - 1).collect();
        spike.sort_unstable();
        let smear: Vec<i64> = (0..20).map(|i| 850 + i * 15).collect();
        let (lo, hi) = (850, 1150);
        assert!(spike_llr(&spike, 1000, lo, hi, 1.0) > 50.0);
        assert!(spike_llr(&smear, 1000, lo, hi, 1.0) < 2.0);
        // The same pile-up 40 bp from the start tested is no spike there.
        let off: Vec<i64> = spike.iter().map(|x| x + 40).collect();
        assert!(spike_llr(&off, 1000, lo, hi, 1.0) < 2.0);
        assert_eq!(spike_llr(&[], 1000, lo, hi, 1.0), 0.0);
    }
}

#[cfg(test)]
mod sibling_tests {
    use super::*;

    // a reported model on +: exons 100-150, 300-320 (a 21 bp microexon), 500-600, 800-900
    const A: [(i64, i64); 4] = [(100, 150), (300, 320), (500, 600), (800, 900)];

    #[test]
    fn the_extra_intron_of_a_sibling() {
        assert_eq!(extra_intron(&[(100, 150), (160, 200), (300, 400)], &[(100, 200), (300, 400)]), Some((150, 160)));
        assert_eq!(extra_intron(&[(100, 150), (160, 200), (300, 350), (360, 400)], &[(100, 200), (300, 400)]), None);
    }

    #[test]
    fn edits_of_one_step() {
        // first exon soft-clipped: starts at the second exon, chain without the first junction
        assert_eq!(sibling_edit(&[(302, 320), (500, 600), (800, 900)], &A, true), Some(SiblingEdit::FirstLost));
        // the same on -: the genomic first junction is the 3' one
        assert_eq!(sibling_edit(&[(302, 320), (500, 600), (800, 900)], &A, false), Some(SiblingEdit::LastLost));
        // a start 40 bp into the second exon is no lost exon
        assert_eq!(sibling_edit(&[(340, 360), (500, 600), (800, 900)], &[(100, 150), (300, 360), (500, 600), (800, 900)], true), None);
        // microexon skipped
        assert_eq!(sibling_edit(&[(100, 150), (500, 600), (800, 900)], &A, true), Some(SiblingEdit::MicroexonSkip));
        // a longer exon skipped: no rule
        assert_eq!(sibling_edit(&[(100, 150), (300, 320), (800, 900)], &A, true), None);
        // one junction 6 bp away
        assert_eq!(sibling_edit(&[(100, 150), (300, 320), (506, 600), (800, 900)], &A, true), Some(SiblingEdit::Shifted));
        // an extra intron inside the third exon
        assert_eq!(sibling_edit(&[(100, 150), (300, 320), (500, 540), (561, 600), (800, 900)], &A, true), Some(SiblingEdit::ExtraIntron));
        // an extra exon inside an intron of A is an inclusion, not an extra intron
        assert_eq!(sibling_edit(&[(100, 150), (300, 320), (500, 600), (700, 720), (800, 900)], &A, true), None);
        assert_eq!(sibling_edit(&A, &A, true), None);
    }

    #[test]
    fn edit_ratios_come_from_the_settings() {
        let mut c = FilterConfig::default();
        assert_eq!([SiblingEdit::FirstLost, SiblingEdit::LastLost, SiblingEdit::MicroexonSkip, SiblingEdit::Shifted, SiblingEdit::ExtraIntron]
            .map(|e| e.max_ratio(&c)), [1.0, 1.0, 0.5, 0.1, 0.05]);
        // a ratio of 0 leaves that edit out and the others as they were
        c.sibling_ratio_shift = 0.0;
        assert_eq!(SiblingEdit::Shifted.max_ratio(&c), 0.0);
        assert_eq!(SiblingEdit::MicroexonSkip.max_ratio(&c), 0.5);
    }
}
