// GTF → Parquet index builder.
//
// Reads a GTF annotation and writes the reference graph to output_dir, classifying generic "UTR" entries as
// five_prime_utr or three_prime_utr from CDS boundaries. Files written: locations, edges, transcripts, hyperedges
// (each transcript's vertex chain 5'->3'), genes, vertex_2_gene, gene_starts, gene_ends, gene_annotations,
// transcript_annotations, utr_regions and exon_intron_regions as .parquet, and counters.json.

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::sync::Arc;

use arrow::array::*;
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::properties::WriterProperties;

use crate::errors::*;

// ─── GTF data structures ─────────────────────────────────────────────────────

#[derive(Debug)]
struct GtfGene {
    identifier: String,
    chromosome: String,
    start: i64,
    end: i64,
    strand: String,
    source: String,
    annotations: HashMap<String, String>,
}

#[derive(Debug)]
struct GtfTranscript {
    identifier: String,
    gene_id: String,
    chromosome: String,
    start: i64,
    end: i64,
    strand: String,
    source: String,
    annotations: HashMap<String, String>,
    exons: Vec<GtfExon>,
}

#[derive(Debug, Clone)]
#[allow(dead_code)]
struct GtfExon {
    chromosome: String,
    start: i64,
    end: i64,
    strand: String,
    source: String,
    annotations: HashMap<String, String>,
}

/// A UTR entry from the GTF (may be generic "UTR" or already classified).
#[derive(Debug, Clone)]
struct GtfUtr {
    chromosome: String,
    start: i64,
    end: i64,
    strand: String,
    transcript_id: String,
    gene_id: String,
    /// Original feature type from GTF: "UTR", "five_prime_utr", or "three_prime_utr"
    original_feature: String,
    /// Classified feature
    classified_feature: Option<String>,
}

/// A CDS entry from the GTF (used to classify generic UTR entries).
#[derive(Debug, Clone)]
struct GtfCds {
    start: i64,
    end: i64,
    transcript_id: String,
}

// ─── Main entry point ────────────────────────────────────────────────────────

/// Build the IsoForge index from a GTF file: reads `gtf_path` and writes the Parquet files to `output_dir`.
/// `genome_build` and `annot_name` name the genome build and the annotation source; transcripts shorter than
/// `min_length` are left out (0 = no filter).
pub fn initialize_from_gtf(
    gtf_path: &Path,
    genome_build: &str,
    annot_name: &str,
    output_dir: &Path,
    min_length: i64,
) -> IsoForgeResult<()> {
    std::fs::create_dir_all(output_dir)?;

    log::info!("[IsoForge] Reading GTF: {:?}", gtf_path);
    let (genes, transcripts, utrs, cds_entries) = read_gtf_file(gtf_path)?;
    eprintln!(
        "[IsoForge]   Parsed {} genes, {} transcripts, {} UTR entries, {} CDS entries",
        genes.len(),
        transcripts.len(),
        utrs.len(),
        cds_entries.len(),
    );

    // Filter by length
    let transcripts: Vec<GtfTranscript> = if min_length > 0 {
        transcripts
            .into_iter()
            .filter(|t| {
                let len: i64 = t.exons.iter().map(|e| e.end - e.start + 1).sum();
                len >= min_length
            })
            .collect()
    } else {
        transcripts
    };

    // Collect gene IDs referenced by surviving transcripts
    let active_gene_ids: HashSet<&str> = transcripts.iter().map(|t| t.gene_id.as_str()).collect();
    let genes: Vec<GtfGene> = genes
        .into_iter()
        .filter(|g| active_gene_ids.contains(g.identifier.as_str()))
        .collect();

    eprintln!(
        "[IsoForge]   After filtering: {} genes, {} transcripts",
        genes.len(),
        transcripts.len()
    );

    // Build the graph (vertices, edges, transcript records)
    let graph = build_graph(&genes, &transcripts, genome_build);
    eprintln!(
        "[IsoForge]   Graph: {} vertices, {} edges",
        graph.locations.len(),
        graph.edges.len()
    );

    // Write all Parquet files
    write_locations(output_dir, &graph)?;
    log::info!("[IsoForge]   ✓ locations.parquet");

    write_edges(output_dir, &graph)?;
    log::info!("[IsoForge]   ✓ edges.parquet");

    write_genes_parquet(output_dir, &genes, &graph)?;
    log::info!("[IsoForge]   ✓ genes.parquet");

    write_transcripts_parquet(output_dir, &graph)?;
    log::info!("[IsoForge]   ✓ transcripts.parquet");

    crate::annotation_hypergraph::write_hyperedges(&output_dir.join("hyperedges.parquet"), &graph.hyperedges)?;
    log::info!("[IsoForge]   ✓ hyperedges.parquet ({} hyperedges)", graph.hyperedges.len());

    write_vertex_2_gene(output_dir, &graph)?;
    log::info!("[IsoForge]   ✓ vertex_2_gene.parquet");

    write_gene_endpoints(output_dir, &graph)?;
    log::info!("[IsoForge]   ✓ gene_starts.parquet, gene_ends.parquet");

    write_counters(output_dir, &graph)?;
    log::info!("[IsoForge]   ✓ counters.json");

    write_gene_annotations(output_dir, &genes, annot_name)?;
    log::info!("[IsoForge]   ✓ gene_annotations.parquet");

    write_transcript_annotations(output_dir, &transcripts, &graph, annot_name)?;
    log::info!("[IsoForge]   ✓ transcript_annotations.parquet");

    // Classify UTRs and write utr_regions.parquet
    let classified_utrs = classify_utrs(utrs, &cds_entries);
    let n_5p = classified_utrs.iter().filter(|u| u.classified_feature.as_deref() == Some("five_prime_utr")).count();
    let n_3p = classified_utrs.iter().filter(|u| u.classified_feature.as_deref() == Some("three_prime_utr")).count();
    write_utr_regions(output_dir, &classified_utrs)?;
    eprintln!(
        "[IsoForge]   UTR classification: {} five_prime_utr, {} three_prime_utr",
        n_5p, n_3p
    );
    log::info!("[IsoForge]   ✓ utr_regions.parquet");

    // Exon and intron regions for feature annotation
    let (n_exon, n_intron) = write_exon_intron_regions(output_dir, &genes, &transcripts)?;
    eprintln!(
        "[IsoForge]   Exon/intron regions: {} exon, {} intron intervals",
        n_exon, n_intron
    );
    log::info!("[IsoForge]   ✓ exon_intron_regions.parquet");

    log::info!("[IsoForge] GTF → Parquet initialization complete.");
    Ok(())
}

// ─── GTF parsing ─────────────────────────────────────────────────────────────

fn read_gtf_file(path: &Path) -> IsoForgeResult<(Vec<GtfGene>, Vec<GtfTranscript>, Vec<GtfUtr>, Vec<GtfCds>)> {
    let file = File::open(path)?;
    let reader = BufReader::with_capacity(1 << 20, file); // 1 MB buffer

    let mut genes: Vec<GtfGene> = Vec::new();
    let mut transcripts: HashMap<String, GtfTranscript> = HashMap::new();
    let mut gene_order: Vec<String> = Vec::new();
    let mut transcript_order: Vec<String> = Vec::new();
    let mut utrs: Vec<GtfUtr> = Vec::new();
    let mut cds_entries: Vec<GtfCds> = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            continue;
        }

        let entry_type = fields[2];
        let chrom = fields[0].to_string();
        let source = fields[1].to_string();
        let start = fields[3].parse::<i64>().unwrap_or(0);
        let end = fields[4].parse::<i64>().unwrap_or(0);
        let strand = fields[6].to_string();
        let attrs = parse_gtf_attributes(fields[8]);

        match entry_type {
            "gene" => {
                let gene_id = attrs.get("gene_id").cloned().unwrap_or_default();
                gene_order.push(gene_id.clone());
                genes.push(GtfGene {
                    identifier: gene_id,
                    chromosome: chrom,
                    start,
                    end,
                    strand,
                    source,
                    annotations: attrs,
                });
            }
            "transcript" => {
                let transcript_id = attrs.get("transcript_id").cloned().unwrap_or_default();
                let gene_id = attrs.get("gene_id").cloned().unwrap_or_default();
                transcript_order.push(transcript_id.clone());
                transcripts.insert(
                    transcript_id.clone(),
                    GtfTranscript {
                        identifier: transcript_id,
                        gene_id,
                        chromosome: chrom,
                        start,
                        end,
                        strand,
                        source,
                        annotations: attrs,
                        exons: Vec::new(),
                    },
                );
            }
            "exon" => {
                let transcript_id = attrs.get("transcript_id").cloned().unwrap_or_default();
                let exon = GtfExon {
                    chromosome: chrom,
                    start,
                    end,
                    strand,
                    source,
                    annotations: attrs,
                };
                if let Some(t) = transcripts.get_mut(&transcript_id) {
                    // Insert in sorted order by start position
                    let pos = t.exons.iter().position(|e| exon.start < e.start);
                    match pos {
                        Some(i) => t.exons.insert(i, exon),
                        None => t.exons.push(exon),
                    }
                }
            }
            // Parse UTR entries (generic "UTR" or pre-classified)
            "UTR" | "five_prime_utr" | "three_prime_utr" => {
                let transcript_id = attrs.get("transcript_id").cloned().unwrap_or_default();
                let gene_id = attrs.get("gene_id").cloned().unwrap_or_default();
                let classified = if entry_type != "UTR" {
                    Some(entry_type.to_string())
                } else {
                    None  // Classified later from CDS boundaries
                };
                utrs.push(GtfUtr {
                    chromosome: chrom,
                    start,
                    end,
                    strand,
                    transcript_id,
                    gene_id,
                    original_feature: entry_type.to_string(),
                    classified_feature: classified,
                });
            }
            // Parse CDS entries (needed for UTR classification)
            "CDS" => {
                let transcript_id = attrs.get("transcript_id").cloned().unwrap_or_default();
                cds_entries.push(GtfCds {
                    start,
                    end,
                    transcript_id,
                });
            }
            _ => {} // Skip other entry types
        }
    }

    // Return transcripts in insertion order
    let ordered_transcripts: Vec<GtfTranscript> = transcript_order
        .into_iter()
        .filter_map(|id| transcripts.remove(&id))
        .collect();

    Ok((genes, ordered_transcripts, utrs, cds_entries))
}

fn parse_gtf_attributes(attr_str: &str) -> HashMap<String, String> {
    let mut attrs = HashMap::new();
    for entry in attr_str.split(';') {
        let entry = entry.trim();
        if entry.is_empty() {
            continue;
        }
        let parts: Vec<&str> = entry.splitn(2, ' ').collect();
        if parts.len() == 2 {
            let key = parts[0].trim().trim_matches('"').to_string();
            let val = parts[1].trim().trim_matches('"').to_string();
            attrs.insert(key, val);
        }
    }
    attrs
}

// ─── Graph construction ──────────────────────────────────────────────────────

/// All the data needed to write Parquet files.
struct AnnotationGraph {
    // Locations: (location_id, genome_build, chromosome, position)
    locations: Vec<(i64, String, String, i64)>,
    // Edges: (edge_id, v1, v2, edge_type, strand)
    edges: Vec<(i64, i64, i64, String, String)>,
    // Transcript records: (tid, gid, tname, start_exon, jn_path, end_exon,
    //                      start_vertex, end_vertex, n_exons, chrom, start_pos, end_pos)
    transcript_records: Vec<TranscriptRecord>,
    // Each transcript as a hyperedge: its vertices 5'->3' with their positions and the edges between them
    hyperedges: Vec<crate::annotation_hypergraph::HyperedgeRow>,
    // Vertex→Gene: (vertex_id, gene_id, strand)
    vertex_2_gene: Vec<(i64, i64, String)>,
    // Gene start endpoints: (gene_id, position, vertex_id)
    gene_starts: Vec<(i64, i64, i64)>,
    // Gene end endpoints: (gene_id, position, vertex_id)
    gene_ends: Vec<(i64, i64, i64)>,
    // Gene ID map: native_id → db_id
    gene_id_map: HashMap<String, i64>,
    // Transcript ID map: native_id → db_id
    transcript_id_map: HashMap<String, i64>,
    // Counters
    gene_counter: i64,
    transcript_counter: i64,
    vertex_counter: i64,
    edge_counter: i64,
}

struct TranscriptRecord {
    transcript_id: i64,
    gene_id: i64,
    transcript_name: String,
    start_exon: i64,
    jn_path: Option<String>,
    end_exon: i64,
    start_vertex: i64,
    end_vertex: i64,
    n_exons: i32,
    chromosome: String,
    start_pos: i64,
    end_pos: i64,
}

fn build_graph(
    genes: &[GtfGene],
    transcripts: &[GtfTranscript],
    genome_build: &str,
) -> AnnotationGraph {
    let mut graph = AnnotationGraph {
        locations: Vec::new(),
        edges: Vec::new(),
        transcript_records: Vec::new(),
        hyperedges: Vec::new(),
        vertex_2_gene: Vec::new(),
        gene_starts: Vec::new(),
        gene_ends: Vec::new(),
        gene_id_map: HashMap::new(),
        transcript_id_map: HashMap::new(),
        gene_counter: 0,
        transcript_counter: 0,
        vertex_counter: 0,
        edge_counter: 0,
    };

    // Vertex dedup: (genome_build, chromosome, position) → vertex_id
    let mut vertex_lookup: HashMap<(String, String, i64), i64> = HashMap::new();
    // Track which genes own each vertex: vertex_id → set of gene_ids
    let mut vertex_genes: HashMap<i64, HashSet<i64>> = HashMap::new();
    // Edge dedup: (v1, v2, edge_type, strand) → edge_id
    let mut edge_lookup: HashMap<(i64, i64, String, String), i64> = HashMap::new();

    // 1. Assign gene IDs
    for gene in genes {
        graph.gene_counter += 1;
        graph.gene_id_map.insert(gene.identifier.clone(), graph.gene_counter);
    }

    // 2. Process transcripts → create vertices and edges
    for transcript in transcripts {
        graph.transcript_counter += 1;
        let db_tid = graph.transcript_counter;
        graph.transcript_id_map.insert(transcript.identifier.clone(), db_tid);

        let db_gid = *graph.gene_id_map.get(&transcript.gene_id).unwrap_or(&0);
        let strand = &transcript.strand;

        // Create vertices for all exon endpoints
        let mut t_vertices: Vec<i64> = Vec::new();
        let mut t_positions: Vec<i64> = Vec::new();

        for exon in &transcript.exons {
            t_positions.push(exon.start);
            t_positions.push(exon.end);
            let v_left = get_or_create_vertex(
                &mut graph,
                &mut vertex_lookup,
                genome_build,
                &exon.chromosome,
                exon.start,
            );
            t_vertices.push(v_left);

            let v_right = get_or_create_vertex(
                &mut graph,
                &mut vertex_lookup,
                genome_build,
                &exon.chromosome,
                exon.end,
            );
            t_vertices.push(v_right);

            // Track vertex→gene
            vertex_genes.entry(v_left).or_default().insert(db_gid);
            vertex_genes.entry(v_right).or_default().insert(db_gid);
        }

        // If minus strand, reverse vertex order (IsoForge convention)
        if strand == "-" {
            t_vertices.reverse();
            t_positions.reverse();
        }

        let start_vertex = *t_vertices.first().unwrap_or(&0);
        let end_vertex = *t_vertices.last().unwrap_or(&0);
        let n_exons = transcript.exons.len() as i32;

        // Create edges between consecutive vertices
        let mut t_edges: Vec<i64> = Vec::new();
        let mut prev_type = "none";

        for i in 0..t_vertices.len().saturating_sub(1) {
            let v1 = t_vertices[i];
            let v2 = t_vertices[i + 1];
            let edge_type = if prev_type == "none" || prev_type == "intron" {
                "exon"
            } else {
                "intron"
            };

            let eid = get_or_create_edge(
                &mut graph,
                &mut edge_lookup,
                v1,
                v2,
                edge_type,
                strand,
            );
            t_edges.push(eid);
            prev_type = edge_type;
        }

        let start_exon = *t_edges.first().unwrap_or(&0);
        let end_exon = *t_edges.last().unwrap_or(&0);

        let jn_path = if t_edges.len() > 2 {
            Some(
                t_edges[1..t_edges.len() - 1]
                    .iter()
                    .map(|e| e.to_string())
                    .collect::<Vec<_>>()
                    .join(","),
            )
        } else {
            None
        };

        // Get chromosome and positions from the transcript entry
        let chrom = &transcript.chromosome;
        let start_pos = transcript.start;
        let end_pos = transcript.end;

        let t_name = transcript
            .annotations
            .get("transcript_name")
            .cloned()
            .unwrap_or_else(|| format!("transcript_{}", db_tid));

        graph.transcript_records.push(TranscriptRecord {
            transcript_id: db_tid,
            gene_id: db_gid,
            transcript_name: t_name,
            start_exon,
            jn_path,
            end_exon,
            start_vertex,
            end_vertex,
            n_exons,
            chromosome: chrom.clone(),
            start_pos,
            end_pos,
        });

        graph.hyperedges.push(crate::annotation_hypergraph::HyperedgeRow {
            transcript_id: db_tid,
            gene_id: db_gid,
            chromosome: chrom.clone(),
            strand: crate::types::Strand::from_str(strand).unwrap_or(crate::types::Strand::Plus),
            vertices: t_vertices.clone(),
            positions: t_positions,
            edges: t_edges.clone(),
        });

        // Gene endpoints: start_vertex and end_vertex
        graph.gene_starts.push((db_gid, start_pos, start_vertex));
        graph.gene_ends.push((db_gid, end_pos, end_vertex));
    }

    // 3. Build vertex_2_gene from accumulated data
    let gid_strand: HashMap<i64, String> = genes
        .iter()
        .filter_map(|g| {
            graph.gene_id_map.get(&g.identifier).map(|&gid| (gid, g.strand.clone()))
        })
        .collect();

    for (vid, gids) in &vertex_genes {
        for &gid in gids {
            let strand = gid_strand
                .get(&gid)
                .cloned()
                .unwrap_or_else(|| "+".to_string());
            graph.vertex_2_gene.push((*vid, gid, strand));
        }
    }
    // sorted by vertex and gene
    graph.vertex_2_gene.sort();

    // Deduplicate gene endpoints
    let mut seen_starts: HashSet<(i64, i64)> = HashSet::new();
    graph.gene_starts.retain(|(gid, pos, _)| seen_starts.insert((*gid, *pos)));
    let mut seen_ends: HashSet<(i64, i64)> = HashSet::new();
    graph.gene_ends.retain(|(gid, pos, _)| seen_ends.insert((*gid, *pos)));

    graph
}

fn get_or_create_vertex(
    graph: &mut AnnotationGraph,
    lookup: &mut HashMap<(String, String, i64), i64>,
    genome_build: &str,
    chromosome: &str,
    position: i64,
) -> i64 {
    let key = (genome_build.to_string(), chromosome.to_string(), position);
    if let Some(&vid) = lookup.get(&key) {
        return vid;
    }

    graph.vertex_counter += 1;
    let vid = graph.vertex_counter;
    lookup.insert(key, vid);
    graph.locations.push((
        vid,
        genome_build.to_string(),
        chromosome.to_string(),
        position,
    ));
    vid
}

fn get_or_create_edge(
    graph: &mut AnnotationGraph,
    lookup: &mut HashMap<(i64, i64, String, String), i64>,
    v1: i64,
    v2: i64,
    edge_type: &str,
    strand: &str,
) -> i64 {
    let key = (v1, v2, edge_type.to_string(), strand.to_string());
    if let Some(&eid) = lookup.get(&key) {
        return eid;
    }

    graph.edge_counter += 1;
    let eid = graph.edge_counter;
    lookup.insert(key, eid);
    graph.edges.push((
        eid,
        v1,
        v2,
        edge_type.to_string(),
        strand.to_string(),
    ));
    eid
}

// ─── Parquet writers ─────────────────────────────────────────────────────────

fn props() -> WriterProperties {
    WriterProperties::builder()
        .set_compression(Compression::ZSTD(Default::default()))
        .build()
}

fn write_to_parquet(path: &Path, schema: Arc<Schema>, batch: &RecordBatch) -> IsoForgeResult<()> {
    let file = File::create(path)?;
    let mut writer = ArrowWriter::try_new(file, schema, Some(props()))?;
    writer.write(batch)?;
    writer.close()?;
    Ok(())
}

fn write_locations(output_dir: &Path, graph: &AnnotationGraph) -> IsoForgeResult<()> {
    let mut ids = Int64Builder::new();
    let mut chroms = StringBuilder::new();
    let mut positions = Int64Builder::new();

    for (lid, _build, chrom, pos) in &graph.locations {
        ids.append_value(*lid);
        chroms.append_value(chrom);
        positions.append_value(*pos);
    }

    let schema = Arc::new(Schema::new(vec![
        Field::new("location_id", DataType::Int64, false),
        Field::new("chromosome", DataType::Utf8, false),
        Field::new("position", DataType::Int64, false),
    ]));

    let batch = RecordBatch::try_new(
        schema.clone(),
        vec![
            Arc::new(ids.finish()),
            Arc::new(chroms.finish()),
            Arc::new(positions.finish()),
        ],
    )?;

    write_to_parquet(&output_dir.join("locations.parquet"), schema, &batch)
}

fn write_edges(output_dir: &Path, graph: &AnnotationGraph) -> IsoForgeResult<()> {
    let mut edge_ids = Int64Builder::new();
    let mut v1s = Int64Builder::new();
    let mut v2s = Int64Builder::new();
    let mut types = StringBuilder::new();
    let mut strands = StringBuilder::new();

    for (eid, v1, v2, etype, strand) in &graph.edges {
        edge_ids.append_value(*eid);
        v1s.append_value(*v1);
        v2s.append_value(*v2);
        types.append_value(etype);
        strands.append_value(strand);
    }

    let schema = Arc::new(Schema::new(vec![
        Field::new("edge_id", DataType::Int64, false),
        Field::new("v1", DataType::Int64, false),
        Field::new("v2", DataType::Int64, false),
        Field::new("edge_type", DataType::Utf8, false),
        Field::new("strand", DataType::Utf8, false),
    ]));

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

    write_to_parquet(&output_dir.join("edges.parquet"), schema, &batch)
}

fn write_genes_parquet(
    output_dir: &Path,
    genes: &[GtfGene],
    graph: &AnnotationGraph,
) -> IsoForgeResult<()> {
    let mut g_ids = Int64Builder::new();
    let mut g_names = StringBuilder::new();
    let mut g_chroms = StringBuilder::new();
    let mut g_starts = Int64Builder::new();
    let mut g_ends = Int64Builder::new();
    let mut g_strands = StringBuilder::new();

    for gene in genes {
        let db_id = *graph.gene_id_map.get(&gene.identifier).unwrap_or(&0);
        let name = gene
            .annotations
            .get("gene_name")
            .cloned()
            .unwrap_or_else(|| gene.identifier.clone());
        g_ids.append_value(db_id);
        g_names.append_value(&name);
        g_chroms.append_value(&gene.chromosome);
        g_starts.append_value(gene.start);
        g_ends.append_value(gene.end);
        g_strands.append_value(&gene.strand);
    }

    let schema = Arc::new(Schema::new(vec![
        Field::new("gene_id", DataType::Int64, false),
        Field::new("gene_name", DataType::Utf8, true),
        Field::new("chromosome", DataType::Utf8, false),
        Field::new("start", DataType::Int64, false),
        Field::new("end", DataType::Int64, false),
        Field::new("strand", DataType::Utf8, false),
    ]));

    let batch = RecordBatch::try_new(
        schema.clone(),
        vec![
            Arc::new(g_ids.finish()),
            Arc::new(g_names.finish()),
            Arc::new(g_chroms.finish()),
            Arc::new(g_starts.finish()),
            Arc::new(g_ends.finish()),
            Arc::new(g_strands.finish()),
        ],
    )?;

    write_to_parquet(&output_dir.join("genes.parquet"), schema, &batch)
}

fn write_transcripts_parquet(output_dir: &Path, graph: &AnnotationGraph) -> IsoForgeResult<()> {
    let mut t_ids = Int64Builder::new();
    let mut g_ids = Int64Builder::new();
    let mut t_names = StringBuilder::new();
    let mut start_exons = Int64Builder::new();
    let mut jn_paths = StringBuilder::new();
    let mut end_exons = Int64Builder::new();
    let mut start_verts = Int64Builder::new();
    let mut end_verts = Int64Builder::new();
    let mut n_exons_col = Int32Builder::new();
    let mut chroms = StringBuilder::new();
    let mut starts = Int64Builder::new();
    let mut ends = Int64Builder::new();

    for t in &graph.transcript_records {
        t_ids.append_value(t.transcript_id);
        g_ids.append_value(t.gene_id);
        t_names.append_value(&t.transcript_name);
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

    let schema = Arc::new(Schema::new(vec![
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
    ]));

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

    write_to_parquet(&output_dir.join("transcripts.parquet"), schema, &batch)
}

fn write_vertex_2_gene(output_dir: &Path, graph: &AnnotationGraph) -> IsoForgeResult<()> {
    let mut verts = Int64Builder::new();
    let mut genes = Int64Builder::new();
    let mut strands = StringBuilder::new();

    for (vid, gid, strand) in &graph.vertex_2_gene {
        verts.append_value(*vid);
        genes.append_value(*gid);
        strands.append_value(strand);
    }

    let schema = Arc::new(Schema::new(vec![
        Field::new("vertex_id", DataType::Int64, false),
        Field::new("gene_id", DataType::Int64, false),
        Field::new("strand", DataType::Utf8, false),
    ]));

    let batch = RecordBatch::try_new(
        schema.clone(),
        vec![
            Arc::new(verts.finish()),
            Arc::new(genes.finish()),
            Arc::new(strands.finish()),
        ],
    )?;

    write_to_parquet(&output_dir.join("vertex_2_gene.parquet"), schema, &batch)
}

fn write_gene_endpoints(output_dir: &Path, graph: &AnnotationGraph) -> IsoForgeResult<()> {
    let schema = Arc::new(Schema::new(vec![
        Field::new("gene_id", DataType::Int64, false),
        Field::new("position", DataType::Int64, false),
        Field::new("vertex_id", DataType::Int64, false),
    ]));

    // Starts
    {
        let mut g_ids = Int64Builder::new();
        let mut positions = Int64Builder::new();
        let mut v_ids = Int64Builder::new();
        for (gid, pos, vid) in &graph.gene_starts {
            g_ids.append_value(*gid);
            positions.append_value(*pos);
            v_ids.append_value(*vid);
        }
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(g_ids.finish()),
                Arc::new(positions.finish()),
                Arc::new(v_ids.finish()),
            ],
        )?;
        write_to_parquet(&output_dir.join("gene_starts.parquet"), schema.clone(), &batch)?;
    }

    // Ends
    {
        let mut g_ids = Int64Builder::new();
        let mut positions = Int64Builder::new();
        let mut v_ids = Int64Builder::new();
        for (gid, pos, vid) in &graph.gene_ends {
            g_ids.append_value(*gid);
            positions.append_value(*pos);
            v_ids.append_value(*vid);
        }
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(g_ids.finish()),
                Arc::new(positions.finish()),
                Arc::new(v_ids.finish()),
            ],
        )?;
        write_to_parquet(&output_dir.join("gene_ends.parquet"), schema, &batch)?;
    }

    Ok(())
}

fn write_counters(output_dir: &Path, graph: &AnnotationGraph) -> IsoForgeResult<()> {
    let mut counters = std::collections::BTreeMap::new();
    counters.insert("genes".to_string(), graph.gene_counter);
    counters.insert("transcripts".to_string(), graph.transcript_counter);
    counters.insert("vertex".to_string(), graph.vertex_counter);
    counters.insert("edge".to_string(), graph.edge_counter);
    counters.insert("observed".to_string(), 0i64);
    counters.insert("dataset".to_string(), 0i64);
    counters.insert("genome_build".to_string(), 1i64);

    let json = serde_json::to_string_pretty(&counters)
        .map_err(|e| IsoForgeError::IoError(std::io::Error::new(std::io::ErrorKind::Other, e.to_string())))?;
    std::fs::write(output_dir.join("counters.json"), json)?;
    Ok(())
}

fn write_gene_annotations(
    output_dir: &Path,
    genes: &[GtfGene],
    annot_name: &str,
) -> IsoForgeResult<()> {
    let mut ids = Int64Builder::new();
    let mut annot_names = StringBuilder::new();
    let mut sources = StringBuilder::new();
    let mut attributes = StringBuilder::new();
    let mut values_col = StringBuilder::new();

    // Enumerate using the same order as graph building (gene_counter = idx + 1)
    for (idx, gene) in genes.iter().enumerate() {
        let db_id = (idx + 1) as i64;

        // Add gene_status = KNOWN
        let mut attrs = gene.annotations.clone();
        if !attrs.contains_key("gene_status") {
            attrs.insert("gene_status".to_string(), "KNOWN".to_string());
        }

        // attributes in name order, not the map's
        let mut attrs: Vec<(&String, &String)> = attrs.iter().collect();
        attrs.sort();
        for (key, val) in attrs {
            ids.append_value(db_id);
            annot_names.append_value(annot_name);
            sources.append_value(&gene.source);
            attributes.append_value(key);
            values_col.append_value(val);
        }
    }

    let schema = Arc::new(Schema::new(vec![
        Field::new("id", DataType::Int64, false),
        Field::new("annot_name", DataType::Utf8, false),
        Field::new("source", DataType::Utf8, false),
        Field::new("attribute", DataType::Utf8, false),
        Field::new("value", DataType::Utf8, false),
    ]));

    let batch = RecordBatch::try_new(
        schema.clone(),
        vec![
            Arc::new(ids.finish()),
            Arc::new(annot_names.finish()),
            Arc::new(sources.finish()),
            Arc::new(attributes.finish()),
            Arc::new(values_col.finish()),
        ],
    )?;

    write_to_parquet(
        &output_dir.join("gene_annotations.parquet"),
        schema,
        &batch,
    )
}

fn write_transcript_annotations(
    output_dir: &Path,
    transcripts: &[GtfTranscript],
    graph: &AnnotationGraph,
    annot_name: &str,
) -> IsoForgeResult<()> {
    let mut ids = Int64Builder::new();
    let mut annot_names = StringBuilder::new();
    let mut sources = StringBuilder::new();
    let mut attributes = StringBuilder::new();
    let mut values_col = StringBuilder::new();

    for transcript in transcripts {
        let db_id = *graph
            .transcript_id_map
            .get(&transcript.identifier)
            .unwrap_or(&0);

        let mut attrs = transcript.annotations.clone();
        // Skip gene-level attributes
        attrs.remove("gene_id");
        attrs.remove("gene_name");
        attrs.remove("gene_type");
        attrs.remove("gene_status");

        if !attrs.contains_key("transcript_status") {
            attrs.insert("transcript_status".to_string(), "KNOWN".to_string());
        }

        let mut attrs: Vec<(&String, &String)> = attrs.iter().collect();
        attrs.sort();
        for (key, val) in attrs {
            ids.append_value(db_id);
            annot_names.append_value(annot_name);
            sources.append_value(&transcript.source);
            attributes.append_value(key);
            values_col.append_value(val);
        }
    }

    let schema = Arc::new(Schema::new(vec![
        Field::new("id", DataType::Int64, false),
        Field::new("annot_name", DataType::Utf8, false),
        Field::new("source", DataType::Utf8, false),
        Field::new("attribute", DataType::Utf8, false),
        Field::new("value", DataType::Utf8, false),
    ]));

    let batch = RecordBatch::try_new(
        schema.clone(),
        vec![
            Arc::new(ids.finish()),
            Arc::new(annot_names.finish()),
            Arc::new(sources.finish()),
            Arc::new(attributes.finish()),
            Arc::new(values_col.finish()),
        ],
    )?;

    write_to_parquet(
        &output_dir.join("transcript_annotations.parquet"),
        schema,
        &batch,
    )
}

// ─── UTR classification ──────────────────────────────────────────────────────

/// Classify generic "UTR" entries as "five_prime_utr" or "three_prime_utr" from the transcript's CDS span: on + strand
/// a UTR before the CDS start is 5' and after the CDS end is 3'; on - strand the reverse. Already classified entries
/// pass through unchanged; UTRs of transcripts without CDS are marked as exon.
fn classify_utrs(mut utrs: Vec<GtfUtr>, cds_entries: &[GtfCds]) -> Vec<GtfUtr> {
    // Build transcript_id → (min_cds_start, max_cds_end) lookup
    let mut cds_bounds: HashMap<String, (i64, i64)> = HashMap::new();
    for cds in cds_entries {
        let entry = cds_bounds
            .entry(cds.transcript_id.clone())
            .or_insert((i64::MAX, i64::MIN));
        entry.0 = entry.0.min(cds.start);
        entry.1 = entry.1.max(cds.end);
    }

    for utr in utrs.iter_mut() {
        // Already classified UTRs
        if utr.classified_feature.is_some() {
            continue;
        }

        // Look up CDS boundaries for this transcript
        if let Some(&(cds_start, cds_end)) = cds_bounds.get(&utr.transcript_id) {
            let is_plus = utr.strand == "+";

            if utr.end <= cds_start {
                // UTR is entirely upstream of CDS
                utr.classified_feature = Some(if is_plus {
                    "five_prime_utr".to_string()
                } else {
                    "three_prime_utr".to_string()
                });
            } else if utr.start >= cds_end {
                // UTR is entirely downstream of CDS
                utr.classified_feature = Some(if is_plus {
                    "three_prime_utr".to_string()
                } else {
                    "five_prime_utr".to_string()
                });
            } else {
                // UTR overlapping the CDS boundary: classified by its midpoint
                let utr_mid = (utr.start + utr.end) / 2;
                let cds_mid = (cds_start + cds_end) / 2;
                if (is_plus && utr_mid < cds_mid) || (!is_plus && utr_mid > cds_mid) {
                    utr.classified_feature = Some("five_prime_utr".to_string());
                } else {
                    utr.classified_feature = Some("three_prime_utr".to_string());
                }
            }
        } else {
            // No CDS for this transcript: mark as exon
            utr.classified_feature = Some("exon".to_string());
        }
    }

    utrs
}

/// Write classified UTR regions to Parquet for use by TSS/TES clustering.
///
/// Output columns: chromosome, start, end, strand, feature, gene_id, transcript_id
fn write_utr_regions(output_dir: &Path, utrs: &[GtfUtr]) -> IsoForgeResult<()> {
    let mut chroms = StringBuilder::new();
    let mut starts = Int64Builder::new();
    let mut ends = Int64Builder::new();
    let mut strands = StringBuilder::new();
    let mut features = StringBuilder::new();
    let mut gene_ids = StringBuilder::new();
    let mut transcript_ids = StringBuilder::new();

    for utr in utrs {
        let feature = utr
            .classified_feature
            .as_deref()
            .unwrap_or(&utr.original_feature);

        // Only write classified UTR entries (skip entries that couldn't be classified)
        if feature == "UTR" {
            continue;
        }

        chroms.append_value(&utr.chromosome);
        starts.append_value(utr.start);
        ends.append_value(utr.end);
        strands.append_value(&utr.strand);
        features.append_value(feature);
        gene_ids.append_value(&utr.gene_id);
        transcript_ids.append_value(&utr.transcript_id);
    }

    let schema = Arc::new(Schema::new(vec![
        Field::new("chromosome", DataType::Utf8, false),
        Field::new("start", DataType::Int64, false),
        Field::new("end", DataType::Int64, false),
        Field::new("strand", DataType::Utf8, false),
        Field::new("feature", DataType::Utf8, false),
        Field::new("gene_id", DataType::Utf8, false),
        Field::new("transcript_id", DataType::Utf8, false),
    ]));

    let batch = RecordBatch::try_new(
        schema.clone(),
        vec![
            Arc::new(chroms.finish()),
            Arc::new(starts.finish()),
            Arc::new(ends.finish()),
            Arc::new(strands.finish()),
            Arc::new(features.finish()),
            Arc::new(gene_ids.finish()),
            Arc::new(transcript_ids.finish()),
        ],
    )?;

    write_to_parquet(&output_dir.join("utr_regions.parquet"), schema, &batch)
}

/// Write exon and intron regions (introns are the gaps between consecutive
/// exons of each transcript) to exon_intron_regions.parquet, with the same
/// schema as utr_regions.parquet.
fn write_exon_intron_regions(
    output_dir: &Path,
    genes: &[GtfGene],
    transcripts: &[GtfTranscript],
) -> IsoForgeResult<(usize, usize)> {
    let mut chroms = StringBuilder::new();
    let mut starts = Int64Builder::new();
    let mut ends = Int64Builder::new();
    let mut strands = StringBuilder::new();
    let mut features = StringBuilder::new();
    let mut gene_ids = StringBuilder::new();
    let mut transcript_ids = StringBuilder::new();

    // Transcript gene_id → gene identifier
    let _gene_map: HashMap<&str, &str> = genes
        .iter()
        .map(|g| (g.identifier.as_str(), g.identifier.as_str()))
        .collect();

    let mut n_exon = 0usize;
    let mut n_intron = 0usize;

    for tx in transcripts {
        if tx.exons.is_empty() {
            continue;
        }

        // Write exon intervals
        for exon in &tx.exons {
            chroms.append_value(&exon.chromosome);
            starts.append_value(exon.start);
            ends.append_value(exon.end);
            strands.append_value(&exon.strand);
            features.append_value("exon");
            gene_ids.append_value(&tx.gene_id);
            transcript_ids.append_value(&tx.identifier);
            n_exon += 1;
        }

        // Compute and write intron intervals (gaps between consecutive exons)
        let mut sorted_exons: Vec<(i64, i64)> = tx.exons.iter().map(|e| (e.start, e.end)).collect();
        sorted_exons.sort_by_key(|&(s, _)| s);

        for i in 0..sorted_exons.len() - 1 {
            let intron_start = sorted_exons[i].1 + 1;
            let intron_end = sorted_exons[i + 1].0 - 1;
            if intron_end >= intron_start {
                chroms.append_value(&tx.chromosome);
                starts.append_value(intron_start);
                ends.append_value(intron_end);
                strands.append_value(&tx.strand);
                features.append_value("intron");
                gene_ids.append_value(&tx.gene_id);
                transcript_ids.append_value(&tx.identifier);
                n_intron += 1;
            }
        }
    }

    let schema = Arc::new(Schema::new(vec![
        Field::new("chromosome", DataType::Utf8, false),
        Field::new("start", DataType::Int64, false),
        Field::new("end", DataType::Int64, false),
        Field::new("strand", DataType::Utf8, false),
        Field::new("feature", DataType::Utf8, false),
        Field::new("gene_id", DataType::Utf8, false),
        Field::new("transcript_id", DataType::Utf8, false),
    ]));

    let batch = RecordBatch::try_new(
        schema.clone(),
        vec![
            Arc::new(chroms.finish()),
            Arc::new(starts.finish()),
            Arc::new(ends.finish()),
            Arc::new(strands.finish()),
            Arc::new(features.finish()),
            Arc::new(gene_ids.finish()),
            Arc::new(transcript_ids.finish()),
        ],
    )?;

    write_to_parquet(
        &output_dir.join("exon_intron_regions.parquet"),
        schema,
        &batch,
    )?;

    Ok((n_exon, n_intron))
}
