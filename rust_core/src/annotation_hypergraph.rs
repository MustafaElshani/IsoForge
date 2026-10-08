// The annotation as a hypergraph.
//
// Vertices are exon boundaries (locations.parquet): a transcript's 5' end, its splice sites, its 3' end. Every
// annotated transcript is a hyperedge, the chain of its vertices 5'->3' (two per exon), and every model built from
// reads is one too. The annotation's hyperedges are stored in hyperedges.parquet, a run's models' in
// novel_hyperedges.parquet. A read is matched by the annotated hyperedges that hold every edge of its chain
// (`containing`).

use std::fs::File;
use std::path::Path;
use std::sync::Arc;

use arrow::array::{Array, Int32Array, Int32Builder, Int64Array, Int64Builder, ListArray, ListBuilder, StringArray, StringBuilder};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::arrow::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::properties::WriterProperties;
use rustc_hash::{FxHashMap, FxHashSet};

use crate::data_structures::{EdgeDict, LocationDict, TranscriptDict};
use crate::errors::{IsoForgeError, IsoForgeResult};
use crate::types::*;

/// One hyperedge: a transcript's vertices 5'->3' (two per exon), their positions, and the edges between consecutive
/// vertices (exon, intron, exon, ..., exon; one exon edge for a single-exon transcript).
#[derive(Debug, Clone, PartialEq)]
pub struct HyperedgeRow {
    pub transcript_id: TranscriptId,
    pub gene_id: GeneId,
    pub chromosome: String,
    pub strand: Strand,
    pub vertices: Vec<VertexId>,
    pub positions: Vec<i64>,
    pub edges: Vec<EdgeId>,
}

impl HyperedgeRow {
    pub fn n_exons(&self) -> usize {
        self.vertices.len() / 2
    }

    /// The splice sites (every vertex but the two ends), as (positions, vertices, the edges between them), 5'->3';
    /// None for a single-exon transcript.
    pub fn internal(&self) -> Option<(Vec<i64>, Vec<VertexId>, Vec<EdgeId>)> {
        let n = self.vertices.len();
        if n < 4 {
            return None;
        }
        Some((self.positions[1..n - 1].to_vec(), self.vertices[1..n - 1].to_vec(), self.edges[1..self.edges.len() - 1].to_vec()))
    }
}

/// The edges of a transcript record in order: first exon, the internal edges of `jn_path`, last exon.
pub fn edge_path(model: &TranscriptModel) -> Vec<EdgeId> {
    let mut path = vec![model.start_exon];
    if let Some(jn) = model.jn_path.as_deref() {
        path.extend(jn.split(',').filter_map(|s| s.parse::<i64>().ok()));
    }
    if model.n_exons > 1 || model.end_exon != model.start_exon {
        path.push(model.end_exon);
    }
    path
}

/// A transcript record's hyperedge, from the vertex ends of its intron edges (`edge_ends`: edge -> (v1, v2, strand),
/// oriented 5'->3') and the vertices' positions; None when an edge or a position is missing.
pub fn derive_row(model: &TranscriptModel, positions: &FxHashMap<VertexId, i64>,
                  edge_ends: &FxHashMap<EdgeId, (VertexId, VertexId, Strand)>) -> Option<HyperedgeRow> {
    let edges = edge_path(model);
    let mut vertices = vec![model.start_vertex];
    let strand;
    if model.n_exons > 1 {
        let jn: Vec<EdgeId> = model.jn_path.as_deref()?.split(',').filter_map(|s| s.parse().ok()).collect();
        if jn.is_empty() {
            return None;
        }
        let mut st = None;
        for id in jn.iter().step_by(2) {
            let &(v1, v2, s) = edge_ends.get(id)?;
            st.get_or_insert(s);
            vertices.push(v1);
            vertices.push(v2);
        }
        strand = st?;
    } else {
        strand = edge_ends.get(&model.start_exon)?.2;
    }
    vertices.push(model.end_vertex);
    let positions = vertices.iter().map(|v| positions.get(v).copied()).collect::<Option<Vec<i64>>>()?;
    Some(HyperedgeRow { transcript_id: model.transcript_id, gene_id: model.gene_id, chromosome: model.chromosome.clone(),
                        strand, vertices, positions, edges })
}

/// Vertex -> position and edge -> (v1, v2, strand) maps of an annotation's dictionaries.
pub fn vertex_and_edge_maps(locations: &LocationDict, edges: &EdgeDict)
                            -> (FxHashMap<VertexId, i64>, FxHashMap<EdgeId, (VertexId, VertexId, Strand)>) {
    let pos = locations.chromosomes.values().flat_map(|c| c.positions.iter().map(|&(p, v)| (v, p))).collect();
    let ends = edges.values().map(|e| (e.edge_id, (e.v1, e.v2, e.strand))).collect();
    (pos, ends)
}

/// An annotated transcript as a hyperedge: the record annotate keeps and, when its chain could be built, the chain.
#[derive(Debug, Clone)]
pub struct Hyperedge {
    pub model: TranscriptModel,
    pub row: Option<HyperedgeRow>,
}

/// The annotation's hyperedges, one per distinct edge path (as the transcript dictionary keeps them), in
/// transcript-ID order, with the edge -> hyperedges postings that answer containment.
#[derive(Debug, Clone, Default)]
pub struct AnnotationHypergraph {
    hyperedges: Vec<Hyperedge>,
    by_tid: FxHashMap<TranscriptId, u32>,
    postings: FxHashMap<EdgeId, Vec<u32>>,
    /// Chains read from hyperedges.parquet (the others were derived from the edge paths).
    pub from_file: usize,
}

impl AnnotationHypergraph {
    /// The hyperedges of the transcript dictionary; chains from `rows` (hyperedges.parquet) where they agree with
    /// the record's edge path, derived from the dictionaries otherwise.
    pub fn build(dict: &TranscriptDict, locations: &LocationDict, edges: &EdgeDict,
                 mut rows: FxHashMap<TranscriptId, HyperedgeRow>) -> Self {
        let (pos, ends) = vertex_and_edge_maps(locations, edges);
        let mut models: Vec<&TranscriptModel> = dict.values().collect();
        models.sort_unstable_by_key(|m| m.transcript_id);
        let mut g = Self::default();
        for m in models {
            let path = edge_path(m);
            let row = match rows.remove(&m.transcript_id) {
                Some(r) if r.edges == path && r.vertices.first() == Some(&m.start_vertex)
                    && r.vertices.last() == Some(&m.end_vertex) => { g.from_file += 1; Some(r) }
                _ => derive_row(m, &pos, &ends),
            };
            let k = g.hyperedges.len() as u32;
            g.by_tid.insert(m.transcript_id, k);
            let distinct: FxHashSet<EdgeId> = path.iter().copied().collect();
            let mut distinct: Vec<EdgeId> = distinct.into_iter().collect();
            distinct.sort_unstable();
            for e in distinct {
                g.postings.entry(e).or_default().push(k);
            }
            g.hyperedges.push(Hyperedge { model: m.clone(), row });
        }
        g
    }

    /// From an index directory: hyperedges.parquet when the index has it.
    pub fn load(annot_dir: &Path, dict: &TranscriptDict, locations: &LocationDict, edges: &EdgeDict) -> IsoForgeResult<Self> {
        let path = annot_dir.join("hyperedges.parquet");
        let rows = if path.exists() { read_hyperedges(&path)?.into_iter().map(|r| (r.transcript_id, r)).collect() }
                   else { FxHashMap::default() };
        Ok(Self::build(dict, locations, edges, rows))
    }

    pub fn len(&self) -> usize {
        self.hyperedges.len()
    }

    pub fn is_empty(&self) -> bool {
        self.hyperedges.is_empty()
    }

    /// Hyperedges in transcript-ID order.
    pub fn iter(&self) -> impl Iterator<Item = &Hyperedge> {
        self.hyperedges.iter()
    }

    pub fn get(&self, tid: TranscriptId) -> Option<&Hyperedge> {
        self.by_tid.get(&tid).map(|&k| &self.hyperedges[k as usize])
    }

    /// The splice sites of an annotated transcript (positions, vertices, edges between them) and its strand.
    pub fn internal(&self, tid: TranscriptId) -> Option<(Vec<i64>, Vec<VertexId>, Vec<EdgeId>, Strand)> {
        let row = self.get(tid)?.row.as_ref()?;
        let (p, v, e) = row.internal()?;
        Some((p, v, e, row.strand))
    }

    /// The annotated transcripts whose hyperedge holds every edge of `query` (a read's chain), in transcript-ID order.
    pub fn containing(&self, query: &[EdgeId]) -> Vec<&TranscriptModel> {
        if query.is_empty() {
            return Vec::new();
        }
        let mut lists: Vec<&Vec<u32>> = Vec::with_capacity(query.len());
        for e in query {
            match self.postings.get(e) {
                Some(l) => lists.push(l),
                None => return Vec::new(),
            }
        }
        lists.sort_unstable_by_key(|l| l.len());
        let mut cur: Vec<u32> = lists[0].clone();
        for l in &lists[1..] {
            cur = intersect_sorted(&cur, l);
            if cur.is_empty() {
                return Vec::new();
            }
        }
        cur.into_iter().map(|k| &self.hyperedges[k as usize].model).collect()
    }
}

fn intersect_sorted(a: &[u32], b: &[u32]) -> Vec<u32> {
    let (mut i, mut j) = (0, 0);
    let mut out = Vec::with_capacity(a.len().min(b.len()));
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                out.push(a[i]);
                i += 1;
                j += 1;
            }
        }
    }
    out
}

fn schema() -> Arc<Schema> {
    let list = || DataType::List(Arc::new(Field::new("item", DataType::Int64, true)));
    Arc::new(Schema::new(vec![
        Field::new("transcript_id", DataType::Int64, false),
        Field::new("gene_id", DataType::Int64, false),
        Field::new("chromosome", DataType::Utf8, false),
        Field::new("strand", DataType::Utf8, false),
        Field::new("n_exons", DataType::Int32, false),
        Field::new("vertices", list(), false),
        Field::new("positions", list(), false),
        Field::new("edges", list(), false),
    ]))
}

/// Writes hyperedges, in the order given.
pub fn write_hyperedges(path: &Path, rows: &[HyperedgeRow]) -> IsoForgeResult<()> {
    let (mut tid, mut gid, mut chrom, mut strand, mut nex) =
        (Int64Builder::new(), Int64Builder::new(), StringBuilder::new(), StringBuilder::new(), Int32Builder::new());
    let (mut verts, mut poss, mut edges) =
        (ListBuilder::new(Int64Builder::new()), ListBuilder::new(Int64Builder::new()), ListBuilder::new(Int64Builder::new()));
    for r in rows {
        tid.append_value(r.transcript_id);
        gid.append_value(r.gene_id);
        chrom.append_value(&r.chromosome);
        strand.append_value(r.strand.as_str());
        nex.append_value(r.n_exons() as i32);
        for v in &r.vertices { verts.values().append_value(*v); }
        verts.append(true);
        for p in &r.positions { poss.values().append_value(*p); }
        poss.append(true);
        for e in &r.edges { edges.values().append_value(*e); }
        edges.append(true);
    }
    let schema = schema();
    let batch = RecordBatch::try_new(schema.clone(), vec![
        Arc::new(tid.finish()), Arc::new(gid.finish()), Arc::new(chrom.finish()), Arc::new(strand.finish()),
        Arc::new(nex.finish()), Arc::new(verts.finish()), Arc::new(poss.finish()), Arc::new(edges.finish()),
    ])?;
    let props = WriterProperties::builder().set_compression(Compression::ZSTD(Default::default())).build();
    let mut w = ArrowWriter::try_new(File::create(path)?, schema, Some(props))?;
    w.write(&batch)?;
    w.close()?;
    Ok(())
}

fn list_values(a: &ListArray, i: usize) -> Vec<i64> {
    let v = a.value(i);
    let v = v.as_any().downcast_ref::<Int64Array>().expect("int64 list");
    (0..v.len()).map(|k| v.value(k)).collect()
}

/// Reads hyperedges written by `write_hyperedges`.
pub fn read_hyperedges(path: &Path) -> IsoForgeResult<Vec<HyperedgeRow>> {
    let reader = ParquetRecordBatchReaderBuilder::try_new(File::open(path)?)?.build()?;
    let mut out = Vec::new();
    for batch in reader {
        let batch = batch?;
        let col = |name: &str| batch.column_by_name(name)
            .ok_or_else(|| IsoForgeError::ConfigError(format!("{}: no column {name}", path.display())));
        let tid = col("transcript_id")?.as_any().downcast_ref::<Int64Array>().expect("int64").clone();
        let gid = col("gene_id")?.as_any().downcast_ref::<Int64Array>().expect("int64").clone();
        let chrom = col("chromosome")?.as_any().downcast_ref::<StringArray>().expect("utf8").clone();
        let strand = col("strand")?.as_any().downcast_ref::<StringArray>().expect("utf8").clone();
        let _ = col("n_exons")?.as_any().downcast_ref::<Int32Array>().expect("int32");
        let verts = col("vertices")?.as_any().downcast_ref::<ListArray>().expect("list").clone();
        let poss = col("positions")?.as_any().downcast_ref::<ListArray>().expect("list").clone();
        let edges = col("edges")?.as_any().downcast_ref::<ListArray>().expect("list").clone();
        for i in 0..batch.num_rows() {
            out.push(HyperedgeRow {
                transcript_id: tid.value(i), gene_id: gid.value(i), chromosome: chrom.value(i).to_string(),
                strand: Strand::from_str(strand.value(i)).unwrap_or(Strand::Plus),
                vertices: list_values(&verts, i), positions: list_values(&poss, i), edges: list_values(&edges, i),
            });
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn model(tid: i64, start_exon: i64, jn: Option<&str>, end_exon: i64, sv: i64, ev: i64, n_exons: i32) -> TranscriptModel {
        TranscriptModel { transcript_id: tid, gene_id: 1, start_exon, jn_path: jn.map(String::from), end_exon,
                          start_vertex: sv, end_vertex: ev, n_exons, chromosome: "chr1".into(), start_pos: 0, end_pos: 0 }
    }

    fn maps() -> (FxHashMap<VertexId, i64>, FxHashMap<EdgeId, (VertexId, VertexId, Strand)>) {
        // vertices 1..6 at 100, 200, 300, 400, 500, 600; a three-exon chain 1-2 | 3-4 | 5-6
        let pos = (1..=6).map(|v| (v, v * 100)).collect();
        let ends = [(10, (1, 2)), (11, (2, 3)), (12, (3, 4)), (13, (4, 5)), (14, (5, 6))]
            .into_iter().map(|(e, (a, b))| (e, (a, b, Strand::Plus))).collect();
        (pos, ends)
    }

    #[test]
    fn chains_follow_the_edge_path() {
        let (pos, ends) = maps();
        let r = derive_row(&model(7, 10, Some("11,12,13"), 14, 1, 6, 3), &pos, &ends).unwrap();
        assert_eq!(r.vertices, vec![1, 2, 3, 4, 5, 6]);
        assert_eq!(r.positions, vec![100, 200, 300, 400, 500, 600]);
        assert_eq!(r.edges, vec![10, 11, 12, 13, 14]);
        assert_eq!(r.internal().unwrap(), (vec![200, 300, 400, 500], vec![2, 3, 4, 5], vec![11, 12, 13]));
        let mono = derive_row(&model(8, 10, None, 10, 1, 2, 1), &pos, &ends).unwrap();
        assert_eq!((mono.vertices, mono.edges), (vec![1, 2], vec![10]));
        assert!(derive_row(&model(9, 10, Some("11,99,13"), 14, 1, 6, 3), &pos, &ends).is_some(), "exon edges are not looked up");
        assert!(derive_row(&model(9, 10, Some("98,12,13"), 14, 1, 6, 3), &pos, &ends).is_none());
    }

    #[test]
    fn containment_is_an_intersection_in_transcript_order() {
        let mut g = AnnotationHypergraph::default();
        for (k, (tid, path)) in [(5, vec![1, 2, 3]), (9, vec![2, 3, 4]), (12, vec![1, 2, 3, 4])].into_iter().enumerate() {
            g.by_tid.insert(tid, k as u32);
            for e in &path { g.postings.entry(*e).or_default().push(k as u32); }
            g.hyperedges.push(Hyperedge { model: model(tid, path[0], None, *path.last().unwrap(), 0, 0, 2), row: None });
        }
        let ids = |q: &[i64]| g.containing(q).iter().map(|m| m.transcript_id).collect::<Vec<_>>();
        assert_eq!(ids(&[2, 3]), vec![5, 9, 12]);
        assert_eq!(ids(&[1, 2]), vec![5, 12]);
        assert_eq!(ids(&[3, 4]), vec![9, 12]);
        assert!(ids(&[1, 4, 7]).is_empty());
    }

    #[test]
    fn hyperedges_round_trip_through_parquet() {
        let rows = vec![
            HyperedgeRow { transcript_id: 3, gene_id: 2, chromosome: "chr2".into(), strand: Strand::Minus,
                           vertices: vec![9, 8, 7, 6], positions: vec![900, 800, 700, 600], edges: vec![21, 22, 23] },
            HyperedgeRow { transcript_id: 4, gene_id: 2, chromosome: "chr2".into(), strand: Strand::Plus,
                           vertices: vec![1, 2], positions: vec![100, 200], edges: vec![30] },
        ];
        let dir = std::env::temp_dir().join(format!("isoforge_hyperedges_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("hyperedges.parquet");
        write_hyperedges(&p, &rows).unwrap();
        assert_eq!(read_hyperedges(&p).unwrap(), rows);
        let _ = std::fs::remove_dir_all(&dir);
    }
}
