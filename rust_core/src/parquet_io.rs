// Parquet/Arrow I/O: reads the reference annotation and writes run results
// (observed reads, abundance, annotations, novel entities).

use std::collections::HashMap;
use std::fs::File;
use std::path::Path;
use std::sync::Arc;

use arrow::array::*;
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::arrow::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::properties::WriterProperties;

use crate::data_structures::*;
use crate::errors::*;
use crate::types::*;
use rustc_hash::FxHashMap;

// ─── Schemas ─────────────────────────────────────────────────────────────────

pub(crate) fn gene_schema() -> Schema {
    Schema::new(vec![
        Field::new("gene_id", DataType::Int64, false),
        Field::new("gene_name", DataType::Utf8, true),
        Field::new("chromosome", DataType::Utf8, false),
        Field::new("start", DataType::Int64, false),
        Field::new("end", DataType::Int64, false),
        Field::new("strand", DataType::Utf8, false),
    ])
}

pub(crate) fn transcript_schema() -> Schema {
    Schema::new(vec![
        Field::new("transcript_id", DataType::Int64, false),
        Field::new("gene_id", DataType::Int64, false),
        Field::new("transcript_name", DataType::Utf8, true),
        Field::new("start_exon", DataType::Int64, false),
        Field::new("jn_path", DataType::Utf8, true),
        Field::new("end_exon", DataType::Int64, false),
        Field::new("start_vertex", DataType::Int64, false),
        Field::new("end_vertex", DataType::Int64, false),
        Field::new("n_exons", DataType::Int32, false),
        Field::new("chromosome", DataType::Utf8, false),
        Field::new("start_pos", DataType::Int64, false),
        Field::new("end_pos", DataType::Int64, false),
    ])
}

#[allow(dead_code)]
pub(crate) fn location_schema() -> Schema {
    Schema::new(vec![
        Field::new("location_id", DataType::Int64, false),
        Field::new("chromosome", DataType::Utf8, false),
        Field::new("position", DataType::Int64, false),
    ])
}

#[allow(dead_code)]
pub(crate) fn edge_schema() -> Schema {
    Schema::new(vec![
        Field::new("edge_id", DataType::Int64, false),
        Field::new("v1", DataType::Int64, false),
        Field::new("v2", DataType::Int64, false),
        Field::new("edge_type", DataType::Utf8, false),
        Field::new("strand", DataType::Utf8, false),
    ])
}

pub(crate) fn observed_schema() -> Schema {
    Schema::new(vec![
        Field::new("read_name", DataType::Utf8, false),
        Field::new("dataset", DataType::Utf8, false),
        Field::new("gene_id", DataType::Int64, false),
        Field::new("transcript_id", DataType::Int64, false),
        Field::new("chromosome", DataType::Utf8, false),
        Field::new("start_pos", DataType::Int64, false),
        Field::new("end_pos", DataType::Int64, false),
        Field::new("strand", DataType::Utf8, false),
        Field::new("read_length", DataType::Int64, false),
        Field::new("n_exons", DataType::Int32, false),
        Field::new("start_vertex", DataType::Int64, false),
        Field::new("end_vertex", DataType::Int64, false),
        Field::new("start_exon", DataType::Int64, false),
        Field::new("end_exon", DataType::Int64, false),
        Field::new("start_delta", DataType::Int64, true),
        Field::new("end_delta", DataType::Int64, true),
        Field::new("fraction_as", DataType::Float64, true),
        Field::new("custom_label", DataType::Utf8, true),
        Field::new("allelic_label", DataType::Utf8, true),
        Field::new("start_support", DataType::Utf8, true),
        Field::new("end_support", DataType::Utf8, true),
        // PAS-aware internal priming
        Field::new("has_pas_signal", DataType::Boolean, true),
        Field::new("pas_motif", DataType::Utf8, true),
        Field::new("pas_distance", DataType::Int64, true),
        Field::new("is_internal_priming", DataType::Boolean, true),
        // TSS/TES cluster join results
        Field::new("tss_site", DataType::Int64, true),
        Field::new("tss_feature", DataType::Utf8, true),
        Field::new("polya_site", DataType::Int64, true),
        Field::new("polya_feature", DataType::Utf8, true),
        Field::new("polya_signal", DataType::Utf8, true),
        Field::new("cluster_supported", DataType::Boolean, true),
        // PolyA tail length from the basecaller pt:i tag
        Field::new("polya_length", DataType::Int64, true),
        // Largest vertex-snap distance of an internal splice site (bp)
        Field::new("max_splice_shift", DataType::Int32, true),
        // Untemplated G's at the 5' end: a capped 5' end
        Field::new("cap_signature", DataType::Boolean, true),
        // Poly(A) tail found in the soft clip, in bases; null when none was found
        Field::new("polya_softclip", DataType::Int64, true),
        // The read's own internal exon boundaries, 5'->3', comma-separated (its hyperedge); null when unspliced
        Field::new("introns", DataType::Utf8, true),
    ])
}

pub(crate) fn annotation_schema() -> Schema {
    Schema::new(vec![
        Field::new("id", DataType::Int64, false),
        Field::new("annot_name", DataType::Utf8, false),
        Field::new("source", DataType::Utf8, false),
        Field::new("attribute", DataType::Utf8, false),
        Field::new("value", DataType::Utf8, false),
    ])
}

pub(crate) fn abundance_schema() -> Schema {
    Schema::new(vec![
        Field::new("gene_id", DataType::Int64, false),
        Field::new("transcript_id", DataType::Int64, false),
        Field::new("dataset", DataType::Utf8, false),
        Field::new("count", DataType::Int64, false),
    ])
}

// ─── Reader ──────────────────────────────────────────────────────────────────

/// Read a Parquet file into record batches.
pub(crate) fn read_parquet_batches(path: &Path) -> IsoForgeResult<Vec<RecordBatch>> {
    let file = File::open(path)?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file)?;
    let reader = builder.build()?;

    let mut batches = Vec::new();
    for batch_result in reader {
        batches.push(batch_result?);
    }
    Ok(batches)
}

/// Load location dict from Parquet.
pub fn load_location_dict(path: &Path) -> IsoForgeResult<LocationDict> {
    let batches = read_parquet_batches(path)?;
    let mut dict = LocationDict::new();

    for batch in &batches {
        let ids = batch
            .column(0)
            .as_any()
            .downcast_ref::<Int64Array>()
            .ok_or(IsoForgeError::InvalidPosType("expected Int64Array for location_id".to_string()))?;
        let chroms = batch
            .column(1)
            .as_any()
            .downcast_ref::<StringArray>()
            .ok_or(IsoForgeError::InvalidPosType("expected StringArray for chromosome".to_string()))?;
        let positions = batch
            .column(2)
            .as_any()
            .downcast_ref::<Int64Array>()
            .ok_or(IsoForgeError::InvalidPosType("expected Int64Array for position".to_string()))?;

        for i in 0..batch.num_rows() {
            let chrom = chroms.value(i);
            let pos = positions.value(i);
            let id = ids.value(i);
            dict.insert(chrom.to_string(), pos, id);
        }
    }

    dict.finalize();
    Ok(dict)
}

/// Load edge dict from Parquet.
pub fn load_edge_dict(path: &Path) -> IsoForgeResult<EdgeDict> {
    let batches = read_parquet_batches(path)?;
    let mut dict = EdgeDict::default();

    for batch in &batches {
        let ids = batch.column(0).as_any().downcast_ref::<Int64Array>().unwrap();
        let v1s = batch.column(1).as_any().downcast_ref::<Int64Array>().unwrap();
        let v2s = batch.column(2).as_any().downcast_ref::<Int64Array>().unwrap();
        let types = batch.column(3).as_any().downcast_ref::<StringArray>().unwrap();
        let strands = batch.column(4).as_any().downcast_ref::<StringArray>().unwrap();

        for i in 0..batch.num_rows() {
            let edge_type = match types.value(i) {
                "exon" => EdgeType::Exon,
                _ => EdgeType::Intron,
            };
            let strand = match strands.value(i) {
                "+" => Strand::Plus,
                _ => Strand::Minus,
            };
            let key = EdgeKey {
                v1: v1s.value(i),
                v2: v2s.value(i),
                edge_type,
            };
            let edge = Edge {
                edge_id: ids.value(i),
                v1: v1s.value(i),
                v2: v2s.value(i),
                edge_type,
                strand,
            };
            dict.insert(key, edge);
        }
    }
    Ok(dict)
}

/// Load transcript dict from Parquet.
pub fn load_transcript_dict(path: &Path) -> IsoForgeResult<TranscriptDict> {
    let batches = read_parquet_batches(path)?;
    let mut dict = TranscriptDict::default();

    for batch in &batches {
        let t_ids = batch.column(0).as_any().downcast_ref::<Int64Array>().unwrap();
        let g_ids = batch.column(1).as_any().downcast_ref::<Int64Array>().unwrap();
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
            let jn = if jn_paths.is_null(i) {
                None
            } else {
                Some(jn_paths.value(i).to_string())
            };

            let model = TranscriptModel {
                transcript_id: t_ids.value(i),
                gene_id: g_ids.value(i),
                start_exon: start_exons.value(i),
                jn_path: jn.clone(),
                end_exon: end_exons.value(i),
                start_vertex: start_verts.value(i),
                end_vertex: end_verts.value(i),
                n_exons: n_exons_col.value(i),
                chromosome: chroms.value(i).to_string(),
                start_pos: starts.value(i),
                end_pos: ends.value(i),
            };

            // Build path key from edges
            let mut all_edges = vec![model.start_exon];
            if let Some(ref jn_str) = jn {
                for edge_str in jn_str.split(',') {
                    if let Ok(eid) = edge_str.parse::<i64>() {
                        all_edges.push(eid);
                    }
                }
            }
            all_edges.push(model.end_exon);

            let key = TranscriptPathKey::new(all_edges);
            dict.insert(key, model);
        }
    }
    Ok(dict)
}

/// Load vertex to gene mapping from Parquet.
pub fn load_vertex_2_gene(path: &Path) -> IsoForgeResult<Vertex2Gene> {
    let batches = read_parquet_batches(path)?;
    let mut v2g = Vertex2Gene::default();

    for batch in &batches {
        let verts = batch.column(0).as_any().downcast_ref::<Int64Array>().unwrap();
        let genes = batch.column(1).as_any().downcast_ref::<Int64Array>().unwrap();
        let strands = batch.column(2).as_any().downcast_ref::<StringArray>().unwrap();

        for i in 0..batch.num_rows() {
            let strand = match strands.value(i) {
                "+" => Strand::Plus,
                _ => Strand::Minus,
            };
            v2g.entry(verts.value(i))
                .or_default()
                .insert((genes.value(i), strand));
        }
    }
    Ok(v2g)
}

/// Load gene start/end dicts from Parquet.
pub fn load_gene_endpoints(path: &Path) -> IsoForgeResult<GeneEndpointDict> {
    let batches = read_parquet_batches(path)?;
    let mut dict = GeneEndpointDict::default();

    for batch in &batches {
        let genes = batch.column(0).as_any().downcast_ref::<Int64Array>().unwrap();
        let positions = batch.column(1).as_any().downcast_ref::<Int64Array>().unwrap();
        let vertex_ids = batch.column(2).as_any().downcast_ref::<Int64Array>().unwrap();

        for i in 0..batch.num_rows() {
            dict.entry(genes.value(i))
                .or_insert_with(FxHashMap::default)
                .insert(positions.value(i), vertex_ids.value(i));
        }
    }
    Ok(dict)
}

// ─── Writer ──────────────────────────────────────────────────────────────────

pub(crate) fn writer_props() -> WriterProperties {
    WriterProperties::builder()
        .set_compression(Compression::ZSTD(Default::default()))
        .build()
}

/// Write observed reads to Parquet in chunks.
pub fn write_observed(path: &Path, records: &[ReadAnnotation]) -> IsoForgeResult<()> {
    let schema = Arc::new(observed_schema());

    let file = File::create(path)?;
    let mut writer = ArrowWriter::try_new(file, schema.clone(), Some(writer_props()))?;

    // Chunks of 1M records
    const CHUNK_SIZE: usize = 1_000_000;

    for chunk in records.chunks(CHUNK_SIZE) {
        let batch = build_observed_batch(&schema, chunk)?;
        writer.write(&batch)?;
    }

    writer.close()?;
    Ok(())
}

/// Build a single RecordBatch from a slice of ReadAnnotation.
fn build_observed_batch(
    schema: &Arc<Schema>,
    records: &[ReadAnnotation],
) -> IsoForgeResult<RecordBatch> {
    let n = records.len();
    let mut read_names = StringBuilder::with_capacity(n, n * 40);
    let mut datasets = StringBuilder::with_capacity(n, n * 20);
    let mut gene_ids = Int64Builder::with_capacity(n);
    let mut transcript_ids = Int64Builder::with_capacity(n);
    let mut chroms = StringBuilder::with_capacity(n, n * 6);
    let mut start_positions = Int64Builder::with_capacity(n);
    let mut end_positions = Int64Builder::with_capacity(n);
    let mut strands = StringBuilder::with_capacity(n, n);
    let mut read_lengths = Int64Builder::with_capacity(n);
    let mut n_exons_col = Int32Builder::with_capacity(n);
    let mut start_vertices = Int64Builder::with_capacity(n);
    let mut end_vertices = Int64Builder::with_capacity(n);
    let mut start_exons = Int64Builder::with_capacity(n);
    let mut end_exons = Int64Builder::with_capacity(n);
    let mut start_deltas = Int64Builder::with_capacity(n);
    let mut end_deltas = Int64Builder::with_capacity(n);
    let mut fractions_as = Float64Builder::with_capacity(n);
    let mut custom_labels = StringBuilder::new();
    let mut allelic_labels = StringBuilder::new();
    let mut start_supports = StringBuilder::new();
    let mut end_supports = StringBuilder::new();
    // PAS-aware internal priming columns
    let mut has_pas_signals = BooleanBuilder::with_capacity(n);
    let mut pas_motifs = StringBuilder::new();
    let mut pas_distances = Int64Builder::with_capacity(n);
    let mut is_internal_primings = BooleanBuilder::with_capacity(n);
    // TSS/TES cluster join columns
    let mut tss_sites = Int64Builder::with_capacity(n);
    let mut tss_features = StringBuilder::new();
    let mut polya_sites = Int64Builder::with_capacity(n);
    let mut polya_features = StringBuilder::new();
    let mut polya_signals = StringBuilder::new();
    let mut cluster_supporteds = BooleanBuilder::with_capacity(n);
    let mut polya_lengths = Int64Builder::with_capacity(n);
    let mut max_splice_shifts = Int32Builder::with_capacity(n);
    let mut cap_signatures = BooleanBuilder::with_capacity(n);
    let mut polya_softclips = Int64Builder::with_capacity(n);
    let mut introns_col = StringBuilder::with_capacity(n, n * 32);
    let mut buf = String::new();

    for r in records {
        read_names.append_value(&r.read_id);
        datasets.append_value(&r.dataset);
        gene_ids.append_value(r.gene_id);
        transcript_ids.append_value(r.transcript_id);
        chroms.append_value(&r.location);
        start_positions.append_value(r.start_pos);
        end_positions.append_value(r.end_pos);
        strands.append_value(match r.strand {
            Strand::Plus => "+",
            Strand::Minus => "-",
        });
        read_lengths.append_value(r.read_length as i64);
        n_exons_col.append_value(r.n_exons);
        start_vertices.append_value(r.start_vertex);
        end_vertices.append_value(r.end_vertex);
        start_exons.append_value(r.start_exon);
        end_exons.append_value(r.end_exon);

        match r.start_delta {
            Some(d) => start_deltas.append_value(d),
            None => start_deltas.append_null(),
        };
        match r.end_delta {
            Some(d) => end_deltas.append_value(d),
            None => end_deltas.append_null(),
        };
        match r.fraction_as {
            Some(f) => fractions_as.append_value(f),
            None => fractions_as.append_null(),
        };
        match &r.custom_label {
            Some(l) => custom_labels.append_value(l),
            None => custom_labels.append_null(),
        };
        match &r.allelic_label {
            Some(l) => allelic_labels.append_value(l),
            None => allelic_labels.append_null(),
        };
        match &r.start_support {
            Some(v) => start_supports.append_value(v.to_string()),
            None => start_supports.append_null(),
        };
        match &r.end_support {
            Some(v) => end_supports.append_value(v.to_string()),
            None => end_supports.append_null(),
        };
        // PAS-aware internal priming fields
        match &r.internal_priming {
            Some(ip) => {
                has_pas_signals.append_value(ip.has_pas_signal);
                match &ip.pas_motif {
                    Some(m) => pas_motifs.append_value(m),
                    None => pas_motifs.append_null(),
                };
                match ip.pas_distance {
                    Some(d) => pas_distances.append_value(d),
                    None => pas_distances.append_null(),
                };
                is_internal_primings.append_value(ip.is_internal_priming);
            }
            None => {
                has_pas_signals.append_null();
                pas_motifs.append_null();
                pas_distances.append_null();
                is_internal_primings.append_null();
            }
        };
        // TSS/TES cluster join fields
        match &r.end_site_info {
            Some(info) => {
                match info.tss_site {
                    Some(s) => tss_sites.append_value(s),
                    None => tss_sites.append_null(),
                };
                match &info.tss_feature {
                    Some(f) => tss_features.append_value(f.as_str()),
                    None => tss_features.append_null(),
                };
                match info.polya_site {
                    Some(s) => polya_sites.append_value(s),
                    None => polya_sites.append_null(),
                };
                match &info.polya_feature {
                    Some(f) => polya_features.append_value(f.as_str()),
                    None => polya_features.append_null(),
                };
                match &info.polya_signal {
                    Some(s) => polya_signals.append_value(s),
                    None => polya_signals.append_null(),
                };
                cluster_supporteds.append_value(info.is_supported);
            }
            None => {
                tss_sites.append_null();
                tss_features.append_null();
                polya_sites.append_null();
                polya_features.append_null();
                polya_signals.append_null();
                cluster_supporteds.append_null();
            }
        };

        // pt:i tag from the basecaller only; left null when the BAM lacks the tag
        match r.polya_tail_len {
            Some(v) => polya_lengths.append_value(v),
            None => polya_lengths.append_null(),
        };
        max_splice_shifts.append_value(r.max_splice_shift);
        cap_signatures.append_value(r.cap_signature);
        match r.polya_tail_softclip {
            Some(v) => polya_softclips.append_value(v),
            None => polya_softclips.append_null(),
        };
        if r.introns.is_empty() {
            introns_col.append_null();
        } else {
            use std::fmt::Write as _;
            buf.clear();
            for (k, p) in r.introns.iter().enumerate() {
                if k > 0 { buf.push(','); }
                let _ = write!(buf, "{p}");
            }
            introns_col.append_value(&buf);
        }
    }

    let batch = RecordBatch::try_new(
        schema.clone(),
        vec![
            Arc::new(read_names.finish()),
            Arc::new(datasets.finish()),
            Arc::new(gene_ids.finish()),
            Arc::new(transcript_ids.finish()),
            Arc::new(chroms.finish()),
            Arc::new(start_positions.finish()),
            Arc::new(end_positions.finish()),
            Arc::new(strands.finish()),
            Arc::new(read_lengths.finish()),
            Arc::new(n_exons_col.finish()),
            Arc::new(start_vertices.finish()),
            Arc::new(end_vertices.finish()),
            Arc::new(start_exons.finish()),
            Arc::new(end_exons.finish()),
            Arc::new(start_deltas.finish()),
            Arc::new(end_deltas.finish()),
            Arc::new(fractions_as.finish()),
            Arc::new(custom_labels.finish()),
            Arc::new(allelic_labels.finish()),
            Arc::new(start_supports.finish()),
            Arc::new(end_supports.finish()),
            // PAS columns
            Arc::new(has_pas_signals.finish()),
            Arc::new(pas_motifs.finish()),
            Arc::new(pas_distances.finish()),
            Arc::new(is_internal_primings.finish()),
            // TSS/TES cluster join columns
            Arc::new(tss_sites.finish()),
            Arc::new(tss_features.finish()),
            Arc::new(polya_sites.finish()),
            Arc::new(polya_features.finish()),
            Arc::new(polya_signals.finish()),
            Arc::new(cluster_supporteds.finish()),
            Arc::new(polya_lengths.finish()),
            Arc::new(max_splice_shifts.finish()),
            Arc::new(cap_signatures.finish()),
            Arc::new(polya_softclips.finish()),
            Arc::new(introns_col.finish()),
        ],
    )?;

    Ok(batch)
}

/// Rewrite the `id` column of an existing annotations Parquet file in place. `remap` returns `None` for rows that
/// should be dropped. `rename` regenerates the `value` column for rows whose value is a generated ID string.
pub fn remap_annotation_ids_by_attribute(
    path: &Path,
    remap: &dyn Fn(&str, i64) -> Option<i64>,
    rename: &dyn Fn(&str, i64, &str) -> Option<String>,
) -> IsoForgeResult<()> {
    let schema = Arc::new(annotation_schema());
    let batches = read_parquet_batches(path)?;

    let tmp_path = path.with_extension("parquet.remap");
    {
        let file = File::create(&tmp_path)?;
        let mut writer = ArrowWriter::try_new(file, schema.clone(), Some(writer_props()))?;

        for batch in &batches {
            let ids = batch
                .column(0)
                .as_any()
                .downcast_ref::<Int64Array>()
                .ok_or(IsoForgeError::InvalidPosType(
                    "expected Int64Array for annotation id".to_string(),
                ))?;
            let attrs = batch.column(3).as_any().downcast_ref::<StringArray>().ok_or(
                IsoForgeError::InvalidPosType("expected StringArray for attribute".to_string()),
            )?;
            let values = batch.column(4).as_any().downcast_ref::<StringArray>().ok_or(
                IsoForgeError::InvalidPosType("expected StringArray for value".to_string()),
            )?;

            let mut keep = Vec::with_capacity(ids.len());
            let mut remapped = Int64Builder::with_capacity(ids.len());
            let mut renamed = StringBuilder::new();
            for i in 0..ids.len() {
                if ids.is_null(i) {
                    keep.push(true);
                    remapped.append_null();
                    renamed.append_value(values.value(i));
                    continue;
                }
                let mapped = remap(attrs.value(i), ids.value(i));
                keep.push(mapped.is_some());
                // Dropped rows keep a placeholder until the mask is applied.
                let new_id = mapped.unwrap_or(ids.value(i));
                remapped.append_value(new_id);
                let old_value = values.value(i);
                match rename(attrs.value(i), new_id, old_value) {
                    Some(v) => renamed.append_value(v),
                    None => renamed.append_value(old_value),
                }
            }

            let mut columns = batch.columns().to_vec();
            columns[0] = Arc::new(remapped.finish());
            columns[4] = Arc::new(renamed.finish());

            let out = if keep.iter().all(|k| *k) {
                RecordBatch::try_new(schema.clone(), columns)?
            } else {
                let mask = BooleanArray::from(keep);
                let full = RecordBatch::try_new(schema.clone(), columns)?;
                arrow::compute::filter_record_batch(&full, &mask)?
            };
            if out.num_rows() > 0 {
                writer.write(&out)?;
            }
        }
        writer.close()?;
    }

    std::fs::rename(&tmp_path, path)?;
    Ok(())
}

// ─── Streaming observed writer + replay ──────────────────────────────────────
// Observed reads are spilled to Parquet window by window with the cluster
// columns left null, then replayed once clustering has filled them in.

/// Streaming writer for observed reads; each window is flushed to disk.
pub struct StreamingObservedWriter {
    writer: ArrowWriter<File>,
    schema: Arc<Schema>,
    n_written: usize,
}

impl StreamingObservedWriter {
    pub fn new(path: &Path) -> IsoForgeResult<Self> {
        let schema = Arc::new(observed_schema());
        let file = File::create(path)?;
        let writer = ArrowWriter::try_new(file, schema.clone(), Some(writer_props()))?;
        Ok(Self { writer, schema, n_written: 0 })
    }

    pub fn flush(&mut self, records: &[ReadAnnotation]) -> IsoForgeResult<()> {
        const CHUNK_SIZE: usize = 1_000_000;
        for chunk in records.chunks(CHUNK_SIZE) {
            let batch = build_observed_batch(&self.schema, chunk)?;
            self.writer.write(&batch)?;
            self.n_written += chunk.len();
        }
        Ok(())
    }

    pub fn close(self) -> IsoForgeResult<usize> {
        self.writer.close()?;
        Ok(self.n_written)
    }
}

/// Reads below the MAPQ threshold that are kept for quantification only: they build no model and are counted in
/// post-processing towards the reported models they fit. `positions` is the parsed read's exon boundaries in 5'->3'
/// order, comma-joined.
pub(crate) fn quant_read_schema() -> Schema {
    Schema::new(vec![
        Field::new("dataset", DataType::Utf8, false),
        Field::new("chromosome", DataType::Utf8, false),
        Field::new("strand", DataType::Utf8, false),
        Field::new("positions", DataType::Utf8, false),
        Field::new("fraction_as", DataType::Float64, true),
        Field::new("has_pas_signal", DataType::Boolean, true),
    ])
}

pub struct QuantReadWriter {
    writer: ArrowWriter<File>,
    schema: Arc<Schema>,
    n_written: usize,
}

impl QuantReadWriter {
    pub fn new(path: &Path) -> IsoForgeResult<Self> {
        let schema = Arc::new(quant_read_schema());
        let file = File::create(path)?;
        let writer = ArrowWriter::try_new(file, schema.clone(), Some(writer_props()))?;
        Ok(Self { writer, schema, n_written: 0 })
    }

    pub fn flush(&mut self, reads: &[crate::bam_reader::ParsedRead]) -> IsoForgeResult<()> {
        if reads.is_empty() {
            return Ok(());
        }
        let mut dataset = StringBuilder::new();
        let mut chromosome = StringBuilder::new();
        let mut strand = StringBuilder::new();
        let mut positions = StringBuilder::new();
        let mut fraction_as = Float64Builder::new();
        let mut has_pas = BooleanBuilder::new();
        for r in reads {
            dataset.append_value(&r.dataset);
            chromosome.append_value(&r.chromosome);
            strand.append_value(r.strand.to_string());
            let p: Vec<String> = r.positions.iter().map(|x| x.to_string()).collect();
            positions.append_value(p.join(","));
            match r.fraction_as { Some(v) => fraction_as.append_value(v), None => fraction_as.append_null() }
            match &r.internal_priming { Some(ip) => has_pas.append_value(ip.has_pas_signal), None => has_pas.append_null() }
        }
        let batch = RecordBatch::try_new(self.schema.clone(), vec![
            Arc::new(dataset.finish()), Arc::new(chromosome.finish()), Arc::new(strand.finish()),
            Arc::new(positions.finish()), Arc::new(fraction_as.finish()), Arc::new(has_pas.finish()),
        ])?;
        self.writer.write(&batch)?;
        self.n_written += reads.len();
        Ok(())
    }

    pub fn close(self) -> IsoForgeResult<usize> {
        self.writer.close()?;
        Ok(self.n_written)
    }
}

/// Graph IDs of a spilled row, so the caller can renumber them.
#[derive(Debug, Clone, Copy)]
pub struct ObservedIds {
    pub start_vertex: i64,
    pub end_vertex: i64,
    pub start_exon: i64,
    pub end_exon: i64,
}

/// Replacement values for a spilled row. Returning `None` from the patch
/// callback drops the read from the output.
#[derive(Debug, Default, Clone)]
pub struct ObservedPatch {
    pub gene_id: i64,
    pub transcript_id: i64,
    pub start_vertex: i64,
    pub end_vertex: i64,
    pub start_exon: i64,
    pub end_exon: i64,
    pub start_support: Option<i64>,
    pub end_support: Option<i64>,
    pub tss_site: Option<i64>,
    pub tss_feature: Option<String>,
    pub polya_site: Option<i64>,
    pub polya_feature: Option<String>,
    pub polya_signal: Option<String>,
    pub cluster_supported: Option<bool>,
}

/// Stream the spill file into the final observed output, applying `patch` to
/// each row in read order, one record batch at a time.
pub fn finalize_observed(
    spill: &Path,
    out: &Path,
    patch: &dyn Fn(usize, ObservedIds) -> Option<ObservedPatch>,
) -> IsoForgeResult<usize> {
    use arrow::array::{BooleanArray, Int64Array, StringArray};

    let schema = Arc::new(observed_schema());
    let file = File::open(spill)?;
    let reader = ParquetRecordBatchReaderBuilder::try_new(file)?.build()?;

    let out_file = File::create(out)?;
    let mut writer = ArrowWriter::try_new(out_file, schema.clone(), Some(writer_props()))?;

    let col_i64 = |b: &RecordBatch, i: usize| -> Int64Array {
        b.column(i).as_any().downcast_ref::<Int64Array>().unwrap().clone()
    };

    let mut read_idx = 0usize;
    let mut n_written = 0usize;

    for batch in reader {
        let batch = batch?;
        let n = batch.num_rows();
        let (sv, ev, sx, ex) = (
            col_i64(&batch, 10),
            col_i64(&batch, 11),
            col_i64(&batch, 12),
            col_i64(&batch, 13),
        );

        let mut keep = Vec::with_capacity(n);
        let mut patches: Vec<ObservedPatch> = Vec::with_capacity(n);
        for i in 0..n {
            let ids = ObservedIds {
                start_vertex: sv.value(i),
                end_vertex: ev.value(i),
                start_exon: sx.value(i),
                end_exon: ex.value(i),
            };
            match patch(read_idx, ids) {
                Some(p) => {
                    keep.push(true);
                    patches.push(p);
                }
                None => {
                    keep.push(false);
                    patches.push(ObservedPatch::default());
                }
            }
            read_idx += 1;
        }

        let mut columns = batch.columns().to_vec();
        columns[2] = Arc::new(Int64Array::from_iter_values(patches.iter().map(|p| p.gene_id)));
        columns[3] =
            Arc::new(Int64Array::from_iter_values(patches.iter().map(|p| p.transcript_id)));
        columns[10] =
            Arc::new(Int64Array::from_iter_values(patches.iter().map(|p| p.start_vertex)));
        columns[11] = Arc::new(Int64Array::from_iter_values(patches.iter().map(|p| p.end_vertex)));
        columns[12] =
            Arc::new(Int64Array::from_iter_values(patches.iter().map(|p| p.start_exon)));
        columns[13] = Arc::new(Int64Array::from_iter_values(patches.iter().map(|p| p.end_exon)));
        columns[19] = Arc::new(StringArray::from_iter(
            patches.iter().map(|p| p.start_support.map(|v| v.to_string())),
        ));
        columns[20] = Arc::new(StringArray::from_iter(
            patches.iter().map(|p| p.end_support.map(|v| v.to_string())),
        ));
        columns[25] = Arc::new(Int64Array::from_iter(patches.iter().map(|p| p.tss_site)));
        columns[26] =
            Arc::new(StringArray::from_iter(patches.iter().map(|p| p.tss_feature.clone())));
        columns[27] = Arc::new(Int64Array::from_iter(patches.iter().map(|p| p.polya_site)));
        columns[28] =
            Arc::new(StringArray::from_iter(patches.iter().map(|p| p.polya_feature.clone())));
        columns[29] =
            Arc::new(StringArray::from_iter(patches.iter().map(|p| p.polya_signal.clone())));
        columns[30] =
            Arc::new(BooleanArray::from_iter(patches.iter().map(|p| p.cluster_supported)));

        let full = RecordBatch::try_new(schema.clone(), columns)?;
        let out_batch = if keep.iter().all(|k| *k) {
            full
        } else {
            arrow::compute::filter_record_batch(&full, &BooleanArray::from(keep))?
        };
        if out_batch.num_rows() > 0 {
            n_written += out_batch.num_rows();
            writer.write(&out_batch)?;
        }
    }

    writer.close()?;
    Ok(n_written)
}

/// Write novelty annotations to Parquet in chunks.
pub fn write_annotations(path: &Path, annotations: &[NoveltyAnnotation]) -> IsoForgeResult<()> {
    let schema = Arc::new(annotation_schema());

    let file = File::create(path)?;
    let mut writer = ArrowWriter::try_new(file, schema.clone(), Some(writer_props()))?;

    const CHUNK_SIZE: usize = 2_000_000;

    for chunk in annotations.chunks(CHUNK_SIZE) {
        let n = chunk.len();
        let mut ids = Int64Builder::with_capacity(n);
        let mut annot_names = StringBuilder::with_capacity(n, n * 10);
        let mut sources = StringBuilder::with_capacity(n, n * 10);
        let mut attributes = StringBuilder::with_capacity(n, n * 15);
        let mut values = StringBuilder::with_capacity(n, n * 10);

        for a in chunk {
            ids.append_value(a.id);
            annot_names.append_value(&a.annot_name);
            sources.append_value(&a.source);
            attributes.append_value(&a.attribute);
            values.append_value(&a.value);
        }

        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(ids.finish()),
                Arc::new(annot_names.finish()),
                Arc::new(sources.finish()),
                Arc::new(attributes.finish()),
                Arc::new(values.finish()),
            ],
        )?;

        writer.write(&batch)?;
    }

    writer.close()?;
    Ok(())
}

// ─── Streaming annotation writer ─────────────────────────────────────────────

/// Streaming writer for novelty annotations; the ArrowWriter stays open across flush() calls.
pub struct StreamingAnnotationWriter {
    writer: ArrowWriter<File>,
    schema: Arc<Schema>,
    n_written: usize,
}

impl StreamingAnnotationWriter {
    /// Open a new streaming annotation writer.
    pub fn new(path: &Path) -> IsoForgeResult<Self> {
        let schema = Arc::new(annotation_schema());
        let file = File::create(path)?;
        let writer = ArrowWriter::try_new(file, schema.clone(), Some(writer_props()))?;
        Ok(Self { writer, schema, n_written: 0 })
    }

    /// Write a batch of annotations to disk in chunks.
    pub fn flush(&mut self, annotations: &[NoveltyAnnotation]) -> IsoForgeResult<()> {
        const CHUNK_SIZE: usize = 2_000_000;
        for chunk in annotations.chunks(CHUNK_SIZE) {
            let n = chunk.len();
            let mut ids = Int64Builder::with_capacity(n);
            let mut annot_names = StringBuilder::with_capacity(n, n * 10);
            let mut sources = StringBuilder::with_capacity(n, n * 10);
            let mut attributes = StringBuilder::with_capacity(n, n * 15);
            let mut values = StringBuilder::with_capacity(n, n * 10);

            for a in chunk {
                ids.append_value(a.id);
                annot_names.append_value(&a.annot_name);
                sources.append_value(&a.source);
                attributes.append_value(&a.attribute);
                values.append_value(&a.value);
            }

            let batch = RecordBatch::try_new(
                self.schema.clone(),
                vec![
                    Arc::new(ids.finish()),
                    Arc::new(annot_names.finish()),
                    Arc::new(sources.finish()),
                    Arc::new(attributes.finish()),
                    Arc::new(values.finish()),
                ],
            )?;

            self.writer.write(&batch)?;
        }
        self.n_written += annotations.len();
        Ok(())
    }

    /// Close the writer and finalize the Parquet file.
    pub fn close(self) -> IsoForgeResult<usize> {
        let total = self.n_written;
        self.writer.close()?;
        Ok(total)
    }
}

/// End variants of annotated transcripts: the variant, its annotated parent, and whether
/// the variant's ends are compatible with the parent's annotated ends.
pub fn write_end_variants(path: &Path, rows: &[(i64, i64, bool)]) -> IsoForgeResult<()> {
    let schema = Arc::new(Schema::new(vec![
        Field::new("transcript_id", DataType::Int64, false),
        Field::new("parent_id", DataType::Int64, false),
        Field::new("compatible", DataType::Boolean, false),
    ]));
    let mut sorted = rows.to_vec();
    sorted.sort_unstable();
    let batch = RecordBatch::try_new(
        schema.clone(),
        vec![
            Arc::new(Int64Array::from_iter_values(sorted.iter().map(|r| r.0))),
            Arc::new(Int64Array::from_iter_values(sorted.iter().map(|r| r.1))),
            Arc::new(BooleanArray::from(sorted.iter().map(|r| r.2).collect::<Vec<bool>>())),
        ],
    )?;
    let file = File::create(path)?;
    let mut writer = ArrowWriter::try_new(file, schema, None)?;
    writer.write(&batch)?;
    writer.close()?;
    Ok(())
}

/// Write abundance counts to Parquet.
pub fn write_abundance<S: std::hash::BuildHasher>(
    path: &Path,
    counts: &HashMap<(GeneId, TranscriptId, Arc<str>), i64, S>,
) -> IsoForgeResult<()> {
    let schema = Arc::new(abundance_schema());

    let mut gene_ids = Int64Builder::new();
    let mut transcript_ids = Int64Builder::new();
    let mut datasets = StringBuilder::new();
    let mut count_vals = Int64Builder::new();

    // Sorted by gene, transcript and dataset; zero counts are dropped.
    let mut rows: Vec<_> = counts.iter().filter(|(_, &count)| count != 0).collect();
    rows.sort_by(|(a, _), (b, _)| a.cmp(b));

    for ((gid, tid, ds), count) in rows {
        gene_ids.append_value(*gid);
        transcript_ids.append_value(*tid);
        datasets.append_value(ds);
        count_vals.append_value(*count);
    }

    let batch = RecordBatch::try_new(
        schema.clone(),
        vec![
            Arc::new(gene_ids.finish()),
            Arc::new(transcript_ids.finish()),
            Arc::new(datasets.finish()),
            Arc::new(count_vals.finish()),
        ],
    )?;

    let file = File::create(path)?;
    let mut writer = ArrowWriter::try_new(file, schema, Some(writer_props()))?;
    writer.write(&batch)?;
    writer.close()?;
    Ok(())
}

/// Write novel genes to Parquet (for updating reference).
pub fn write_novel_genes(path: &Path, genes: &[Gene]) -> IsoForgeResult<()> {
    let schema = Arc::new(gene_schema());

    let mut gene_ids = Int64Builder::new();
    let mut gene_names = StringBuilder::new();
    let mut chroms = StringBuilder::new();
    let mut starts = Int64Builder::new();
    let mut ends = Int64Builder::new();
    let mut strands = StringBuilder::new();

    for g in genes {
        gene_ids.append_value(g.gene_id);
        gene_names.append_null();
        chroms.append_value(&g.chromosome);
        starts.append_value(g.start);
        ends.append_value(g.end);
        strands.append_value(match g.strand {
            Strand::Plus => "+",
            Strand::Minus => "-",
        });
    }

    let batch = RecordBatch::try_new(
        schema.clone(),
        vec![
            Arc::new(gene_ids.finish()),
            Arc::new(gene_names.finish()),
            Arc::new(chroms.finish()),
            Arc::new(starts.finish()),
            Arc::new(ends.finish()),
            Arc::new(strands.finish()),
        ],
    )?;

    let file = File::create(path)?;
    let mut writer = ArrowWriter::try_new(file, schema, Some(writer_props()))?;
    writer.write(&batch)?;
    writer.close()?;
    Ok(())
}

/// Write novel transcripts to Parquet in chunks.
pub fn write_novel_transcripts(path: &Path, transcripts: &[TranscriptModel]) -> IsoForgeResult<()> {
    let schema = Arc::new(transcript_schema());

    let file = File::create(path)?;
    let mut writer = ArrowWriter::try_new(file, schema.clone(), Some(writer_props()))?;

    const CHUNK_SIZE: usize = 1_000_000;

    for chunk in transcripts.chunks(CHUNK_SIZE) {
        let n = chunk.len();
        let mut t_ids = Int64Builder::with_capacity(n);
        let mut g_ids = Int64Builder::with_capacity(n);
        let mut t_names = StringBuilder::new();
        let mut start_exons = Int64Builder::with_capacity(n);
        let mut jn_paths = StringBuilder::new();
        let mut end_exons = Int64Builder::with_capacity(n);
        let mut start_verts = Int64Builder::with_capacity(n);
        let mut end_verts = Int64Builder::with_capacity(n);
        let mut n_exons_col = Int32Builder::with_capacity(n);
        let mut chroms = StringBuilder::with_capacity(n, n * 6);
        let mut starts = Int64Builder::with_capacity(n);
        let mut ends = Int64Builder::with_capacity(n);

        for t in chunk {
            t_ids.append_value(t.transcript_id);
            g_ids.append_value(t.gene_id);
            t_names.append_null();
            start_exons.append_value(t.start_exon);
            match &t.jn_path {
                Some(jn) => jn_paths.append_value(jn),
                None => jn_paths.append_null(),
            };
            end_exons.append_value(t.end_exon);
            start_verts.append_value(t.start_vertex);
            end_verts.append_value(t.end_vertex);
            n_exons_col.append_value(t.n_exons);
            chroms.append_value(&t.chromosome);
            starts.append_value(t.start_pos);
            ends.append_value(t.end_pos);
        }

        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(t_ids.finish()),
                Arc::new(g_ids.finish()),
                Arc::new(t_names.finish()),
                Arc::new(start_exons.finish()),
                Arc::new(jn_paths.finish()),
                Arc::new(end_exons.finish()),
                Arc::new(start_verts.finish()),
                Arc::new(end_verts.finish()),
                Arc::new(n_exons_col.finish()),
                Arc::new(chroms.finish()),
                Arc::new(starts.finish()),
                Arc::new(ends.finish()),
            ],
        )?;

        writer.write(&batch)?;
    }

    writer.close()?;
    Ok(())
}

/// Write novel edges to Parquet.
pub fn write_novel_edges(path: &Path, edges: &[Edge]) -> IsoForgeResult<()> {
    let schema = Arc::new(edge_schema());

    let mut edge_ids = Int64Builder::with_capacity(edges.len());
    let mut v1s = Int64Builder::with_capacity(edges.len());
    let mut v2s = Int64Builder::with_capacity(edges.len());
    let mut types = StringBuilder::new();
    let mut strands = StringBuilder::new();

    for e in edges {
        edge_ids.append_value(e.edge_id);
        v1s.append_value(e.v1);
        v2s.append_value(e.v2);
        types.append_value(e.edge_type.as_str());
        strands.append_value(match e.strand {
            Strand::Plus => "+",
            Strand::Minus => "-",
        });
    }

    let batch = RecordBatch::try_new(
        schema.clone(),
        vec![
            Arc::new(edge_ids.finish()),
            Arc::new(v1s.finish()),
            Arc::new(v2s.finish()),
            Arc::new(types.finish()),
            Arc::new(strands.finish()),
        ],
    )?;

    let file = File::create(path)?;
    let mut writer = ArrowWriter::try_new(file, schema, Some(writer_props()))?;
    writer.write(&batch)?;
    writer.close()?;
    Ok(())
}

/// Write novel locations to Parquet.
pub fn write_novel_locations(path: &Path, locations: &[Location]) -> IsoForgeResult<()> {
    let schema = Arc::new(location_schema());

    let mut loc_ids = Int64Builder::with_capacity(locations.len());
    let mut chroms = StringBuilder::new();
    let mut positions = Int64Builder::with_capacity(locations.len());

    for l in locations {
        loc_ids.append_value(l.location_id);
        chroms.append_value(&l.chromosome);
        positions.append_value(l.position);
    }

    let batch = RecordBatch::try_new(
        schema.clone(),
        vec![
            Arc::new(loc_ids.finish()),
            Arc::new(chroms.finish()),
            Arc::new(positions.finish()),
        ],
    )?;

    let file = File::create(path)?;
    let mut writer = ArrowWriter::try_new(file, schema, Some(writer_props()))?;
    writer.write(&batch)?;
    writer.close()?;
    Ok(())
}

// ─── Gene loader (for building interval trees) ───────────────────────────────

/// Gene types (the GTF's `gene_type`, or `gene_biotype`) by gene, from the index's gene_annotations.parquet.
pub fn load_gene_types(path: &Path) -> IsoForgeResult<rustc_hash::FxHashMap<GeneId, String>> {
    let mut types = rustc_hash::FxHashMap::default();
    for batch in &read_parquet_batches(path)? {
        let ids = batch.column(0).as_any().downcast_ref::<Int64Array>().unwrap();
        let attributes = batch.column(3).as_any().downcast_ref::<StringArray>().unwrap();
        let values = batch.column(4).as_any().downcast_ref::<StringArray>().unwrap();
        for i in 0..batch.num_rows() {
            if matches!(attributes.value(i), "gene_type" | "gene_biotype") {
                types.entry(ids.value(i)).or_insert_with(|| values.value(i).to_string());
            }
        }
    }
    Ok(types)
}

/// Load genes from Parquet for building interval trees.
pub fn load_genes(path: &Path) -> IsoForgeResult<Vec<Gene>> {
    let batches = read_parquet_batches(path)?;
    let mut genes = Vec::new();

    for batch in &batches {
        let ids = batch.column(0).as_any().downcast_ref::<Int64Array>().unwrap();
        let chroms = batch.column(2).as_any().downcast_ref::<StringArray>().unwrap();
        let starts = batch.column(3).as_any().downcast_ref::<Int64Array>().unwrap();
        let ends = batch.column(4).as_any().downcast_ref::<Int64Array>().unwrap();
        let strands_col = batch.column(5).as_any().downcast_ref::<StringArray>().unwrap();

        for i in 0..batch.num_rows() {
            let strand = match strands_col.value(i) {
                "+" => Strand::Plus,
                _ => Strand::Minus,
            };
            genes.push(Gene {
                gene_id: ids.value(i),
                strand,
                chromosome: chroms.value(i).to_string(),
                start: starts.value(i),
                end: ends.value(i),
            });
        }
    }
    Ok(genes)
}

/// Load counters from JSON file.
pub fn load_counters(path: &Path) -> IsoForgeResult<(i64, i64, i64, i64, i64, i64)> {
    let data = std::fs::read_to_string(path)?;
    let map: HashMap<String, i64> = serde_json::from_str(&data)
        .map_err(|e| IsoForgeError::IoError(std::io::Error::new(std::io::ErrorKind::Other, e.to_string())))?;

    Ok((
        *map.get("genes").unwrap_or(&0),
        *map.get("transcripts").unwrap_or(&0),
        *map.get("vertex").unwrap_or(&0),
        *map.get("edge").unwrap_or(&0),
        *map.get("observed").unwrap_or(&0),
        *map.get("dataset").unwrap_or(&0),
    ))
}

/// Optional audit of every assigned read, including reads dropped before export.
/// Original IDs are run-local; final IDs are given only for retained reads.
pub fn write_read_decisions<F>(
    spill: &Path, output: &Path, cores: &[crate::annotation_engine::ObsCore],
    joins: &[crate::types::ReadEndSiteInfo], keep: &[bool], final_id: F,
) -> IsoForgeResult<()>
where F: Fn(i64) -> i64 {
    use std::io::Write;
    let mut out = std::io::BufWriter::new(File::create(output)?);
    writeln!(out, "read_idx\tread_name\tdataset\tchrom\tstrand\tstart\tend\tn_exons\toriginal_transcript_id\tfinal_transcript_id\tretained_for_quantification\tpasses_read_gates\tend_rejection_reason\ttss\ttes\ttss_feature\tfrac_a\tinternal_priming\thint_flags\tmax_splice_shift")?;
    let reader = ParquetRecordBatchReaderBuilder::try_new(File::open(spill)?)?
        .with_batch_size(65536).build()?;
    let mut idx = 0;
    let empty = crate::types::ReadEndSiteInfo::default();
    for batch in reader {
        let batch = batch?;
        let names = batch.column_by_name("read_name").unwrap().as_any().downcast_ref::<StringArray>().unwrap();
        for row in 0..batch.num_rows() {
            let c = &cores[idx]; let info = joins.get(idx).unwrap_or(&empty);
            let retained = keep.get(idx).copied().unwrap_or(false);
            let tid = if retained { final_id(c.transcript_id).to_string() } else { String::new() };
            writeln!(out, "{idx}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{tid}\t{retained}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                names.value(row), c.dataset, c.location, c.strand.as_str(), c.start_pos, c.end_pos,
                c.n_exons, c.transcript_id, info.is_supported,
                info.rejection_reason.as_deref().unwrap_or(if joins.is_empty() { "clustering_disabled" } else { "" }),
                info.tss_site.map(|x|x.to_string()).unwrap_or_default(),
                info.polya_site.map(|x|x.to_string()).unwrap_or_default(),
                info.tss_feature.as_ref().map(|x|x.as_str()).unwrap_or(""),
                c.fraction_as.map(|x|x.to_string()).unwrap_or_default(), c.is_internal_priming, c.hint_flags, c.max_splice_shift)?;
            idx += 1;
        }
    }
    Ok(())
}
