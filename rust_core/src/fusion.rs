// Fusion transcripts from chimeric alignments.
//
// - A read's segments are sorted along the read; each pair of consecutive segments gives a junction between a 5' and
//   a 3' partner, and a breakpoint within `snap_window` of an annotated exon boundary of its gene moves onto it.
// - Junctions are clustered by breakpoint into calls; each call gets a category (inter-chromosomal,
//   intra-chromosomal, read-through, fold-back), a filter status and the isoforms of its reads.
// - Reads with the same ordered path of two or more calls make one multi-segment fusion.
// - The calls are written as TSV, BEDPE, GTF and Parquet files.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::sync::Arc;

use arrow::array::*;
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use coitrees::{COITree, Interval, IntervalTree};
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::arrow::ArrowWriter;
use rust_htslib::bam;
use rust_htslib::bam::record::{Aux, Cigar};
use rustc_hash::{FxHashMap, FxHashSet};

use crate::errors::IsoForgeResult;
use crate::post_processing::PostDB;

/// One aligned part of a chimeric read.
#[derive(Debug, Clone)]
pub struct ChimericSegment {
    pub read_name: String,
    pub dataset: String,
    pub chromosome: String,
    pub reverse: bool,
    pub supplementary: bool,
    pub mapq: u8,
    /// 1-based, inclusive.
    pub ref_start: i64,
    pub ref_end: i64,
    /// The covered part of the read in its sequenced orientation, 0-based half-open.
    pub query_start: i32,
    pub query_end: i32,
    pub read_length: i32,
    /// Introns as [first, last] intronic base pairs, flattened.
    pub introns: Vec<i64>,
    /// Not an alignment: a part of the read the aligner left unaligned, placed on annotated exons (`fusion_place`).
    pub placed: bool,
    /// A placed part under `min_placed` bases, named with the other reads' pieces at the same alignment end
    /// (`fusion_place`, short parts placed as a class; `SHORT_MIN` bases at least).
    pub short_class: bool,
    /// For a placed part: the reads with an unaligned piece at the same end of their alignment (this one among them),
    /// and its score above the score a placed part needs.
    pub pieces: u32,
    pub excess: i32,
}

impl ChimericSegment {
    /// The segment of a primary or supplementary record that carries an SA tag.
    pub fn from_record(record: &bam::Record, ref_names: &[String], dataset: &str) -> Option<Self> {
        if !matches!(record.aux(b"SA"), Ok(Aux::String(_))) { return None; }
        Self::of_record(record, ref_names, dataset)
    }

    /// A primary or supplementary record's alignment as a segment.
    pub fn of_record(record: &bam::Record, ref_names: &[String], dataset: &str) -> Option<Self> {
        if record.is_unmapped() || record.is_secondary() { return None; }
        let tid = record.tid();
        if tid < 0 || tid as usize >= ref_names.len() { return None; }
        let cigar = record.cigar();
        let (mut lead, mut trail, mut qlen) = (0i32, 0i32, 0i32);
        let mut seen_aligned = false;
        let mut introns = Vec::new();
        let mut pos = record.pos() + 1;
        for op in cigar.iter() {
            match *op {
                Cigar::SoftClip(n) | Cigar::HardClip(n) => {
                    if seen_aligned { trail += n as i32 } else { lead += n as i32 }
                }
                Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => { seen_aligned = true; qlen += n as i32; pos += n as i64; }
                Cigar::Ins(n) => { seen_aligned = true; qlen += n as i32; }
                Cigar::Del(n) => { pos += n as i64; }
                Cigar::RefSkip(n) => { introns.push(pos); introns.push(pos + n as i64 - 1); pos += n as i64; }
                Cigar::Pad(_) => {}
            }
        }
        let read_length = lead + qlen + trail;
        // Query coordinates on the read as sequenced: a reverse alignment's CIGAR runs along the reverse complement.
        let (qs, qe) = if record.is_reverse() { (trail, trail + qlen) } else { (lead, lead + qlen) };
        let dataset = if !dataset.is_empty() {
            dataset.to_string()
        } else {
            match record.aux(b"RG") { Ok(Aux::String(s)) => s.to_string(), _ => String::new() }
        };
        Some(Self {
            read_name: String::from_utf8_lossy(record.qname()).to_string(),
            dataset,
            chromosome: ref_names[tid as usize].clone(),
            reverse: record.is_reverse(),
            supplementary: record.is_supplementary(),
            mapq: record.mapq(),
            ref_start: record.pos() + 1,
            ref_end: cigar.end_pos(),
            query_start: qs,
            query_end: qe,
            read_length,
            introns,
            placed: false,
            short_class: false,
            pieces: 0,
            excess: 0,
        })
    }

    /// The parts of a record's alignment between introns that join two genes (see `GeneSpans::joins_genes`). A record
    /// with an SA tag gives its parts even when nothing splits it; any other record only when such an intron does.
    pub fn split_at_gene_joins(record: &bam::Record, ref_names: &[String], dataset: &str, spans: &GeneSpans) -> Vec<Self> {
        if record.is_unmapped() || record.is_secondary() { return Vec::new(); }
        let has_sa = matches!(record.aux(b"SA"), Ok(Aux::String(_)));
        if !has_sa && record.is_supplementary() { return Vec::new(); }
        let tid = record.tid();
        if tid < 0 || tid as usize >= ref_names.len() { return Vec::new(); }
        let chrom = &ref_names[tid as usize];
        let cigar = record.cigar();
        let (mut lead, mut trail, mut q) = (0i32, 0i32, 0i32);
        let mut seen_aligned = false;
        let mut pos = record.pos() + 1;
        // parts in alignment orientation: (ref_start, ref_end, query_start, query_end, introns)
        let mut parts: Vec<(i64, i64, i32, i32, Vec<i64>)> = Vec::new();
        let (mut rs, mut qs, mut introns) = (pos, 0i32, Vec::new());
        for op in cigar.iter() {
            match *op {
                Cigar::SoftClip(n) | Cigar::HardClip(n) => {
                    if seen_aligned { trail += n as i32 } else { lead += n as i32 }
                }
                Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => { seen_aligned = true; q += n as i32; pos += n as i64; }
                Cigar::Ins(n) => { seen_aligned = true; q += n as i32; }
                Cigar::Del(n) => { pos += n as i64; }
                Cigar::RefSkip(n) => {
                    let (first, last) = (pos, pos + n as i64 - 1);
                    if spans.joins_genes(chrom, first, last) {
                        parts.push((rs, first - 1, qs, q, std::mem::take(&mut introns)));
                        rs = last + 1; qs = q;
                    } else {
                        introns.push(first); introns.push(last);
                    }
                    pos += n as i64;
                }
                Cigar::Pad(_) => {}
            }
        }
        parts.push((rs, pos - 1, qs, q, introns));
        if parts.len() < 2 && !has_sa { return Vec::new(); }
        let read_length = lead + q + trail;
        let dataset = if !dataset.is_empty() { dataset.to_string() }
            else { match record.aux(b"RG") { Ok(Aux::String(s)) => s.to_string(), _ => String::new() } };
        let name = String::from_utf8_lossy(record.qname()).to_string();
        parts.into_iter().map(|(ref_start, ref_end, a, b, introns)| {
            // query coordinates on the read as sequenced: a reverse alignment's CIGAR runs along the reverse complement
            let (qa, qb) = (lead + a, lead + b);
            let (query_start, query_end) = if record.is_reverse() { (read_length - qb, read_length - qa) } else { (qa, qb) };
            Self { read_name: name.clone(), dataset: dataset.clone(), chromosome: chrom.clone(), reverse: record.is_reverse(),
                   supplementary: record.is_supplementary(), mapq: record.mapq(), ref_start, ref_end, query_start, query_end,
                   read_length, introns, placed: false, short_class: false, pieces: 0, excess: 0 }
        }).collect()
    }

    /// Aligned reference blocks (exonic parts), 1-based inclusive.
    fn blocks(&self) -> Vec<(i64, i64)> {
        let mut out = Vec::with_capacity(self.introns.len() / 2 + 1);
        let mut s = self.ref_start;
        for p in self.introns.chunks(2) {
            out.push((s, p[0] - 1));
            s = p[1] + 1;
        }
        out.push((s, self.ref_end));
        out
    }
}

fn segment_schema() -> Schema {
    Schema::new(vec![
        Field::new("read_name", DataType::Utf8, false),
        Field::new("dataset", DataType::Utf8, false),
        Field::new("chromosome", DataType::Utf8, false),
        Field::new("reverse", DataType::Boolean, false),
        Field::new("supplementary", DataType::Boolean, false),
        Field::new("mapq", DataType::Int32, false),
        Field::new("ref_start", DataType::Int64, false),
        Field::new("ref_end", DataType::Int64, false),
        Field::new("query_start", DataType::Int32, false),
        Field::new("query_end", DataType::Int32, false),
        Field::new("read_length", DataType::Int32, false),
        Field::new("introns", DataType::Utf8, false),
    ])
}

pub fn write_segments(path: &Path, segs: &[ChimericSegment]) -> IsoForgeResult<usize> {
    let schema = Arc::new(segment_schema());
    let mut w = ArrowWriter::try_new(File::create(path)?, schema.clone(), Some(crate::parquet_io::writer_props()))?;
    for chunk in segs.chunks(1 << 18) {
        let introns: Vec<String> = chunk.iter()
            .map(|s| s.introns.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",")).collect();
        let batch = RecordBatch::try_new(schema.clone(), vec![
            Arc::new(StringArray::from_iter_values(chunk.iter().map(|s| s.read_name.as_str()))),
            Arc::new(StringArray::from_iter_values(chunk.iter().map(|s| s.dataset.as_str()))),
            Arc::new(StringArray::from_iter_values(chunk.iter().map(|s| s.chromosome.as_str()))),
            Arc::new(BooleanArray::from(chunk.iter().map(|s| s.reverse).collect::<Vec<_>>())),
            Arc::new(BooleanArray::from(chunk.iter().map(|s| s.supplementary).collect::<Vec<_>>())),
            Arc::new(Int32Array::from_iter_values(chunk.iter().map(|s| s.mapq as i32))),
            Arc::new(Int64Array::from_iter_values(chunk.iter().map(|s| s.ref_start))),
            Arc::new(Int64Array::from_iter_values(chunk.iter().map(|s| s.ref_end))),
            Arc::new(Int32Array::from_iter_values(chunk.iter().map(|s| s.query_start))),
            Arc::new(Int32Array::from_iter_values(chunk.iter().map(|s| s.query_end))),
            Arc::new(Int32Array::from_iter_values(chunk.iter().map(|s| s.read_length))),
            Arc::new(StringArray::from_iter_values(introns.iter().map(|s| s.as_str()))),
        ])?;
        w.write(&batch)?;
    }
    w.close()?;
    Ok(segs.len())
}

/// A part in a pseudoautosomal region of chrX or chrY with mapping quality below `min_mapq` is placed on chrX (PAR1
/// at the same coordinates, PAR2 shifted) with mapping quality `min_mapq`. The regions are in GRCh38 coordinates.
pub fn name_pseudoautosomal(segs: &mut [ChimericSegment], min_mapq: u8) {
    const PAR1: (i64, i64) = (10_001, 2_781_479);
    const PAR2_Y: (i64, i64) = (56_887_903, 57_217_415);
    const PAR2_X: (i64, i64) = (155_701_383, 156_030_895);
    let within = |r: (i64, i64), s: &ChimericSegment| s.ref_start >= r.0 && s.ref_end <= r.1;
    let mut n = 0;
    for s in segs.iter_mut() {
        if s.mapq >= min_mapq { continue; }
        let fix = match s.chromosome.as_str() {
            "chrY" if within(PAR1, s) => Some(0),
            "chrY" if within(PAR2_Y, s) => Some(PAR2_X.0 - PAR2_Y.0),
            "chrX" if within(PAR1, s) || within(PAR2_X, s) => Some(0),
            _ => None,
        };
        let Some(shift) = fix else { continue };
        if shift != 0 { s.ref_start += shift; s.ref_end += shift; for x in s.introns.iter_mut() { *x += shift; } }
        if s.chromosome == "chrY" { s.chromosome = "chrX".to_string(); }
        s.mapq = min_mapq;
        n += 1;
    }
    if n > 0 { log::info!("[Fusion] {} parts in the pseudoautosomal regions named by their chrX copy", n); }
}

pub fn load_segments(path: &Path) -> IsoForgeResult<Vec<ChimericSegment>> {
    let mut out = Vec::new();
    if !path.exists() { return Ok(out); }
    let reader = ParquetRecordBatchReaderBuilder::try_new(File::open(path)?)?.build()?;
    for batch in reader {
        let b = batch?;
        let s = |i: usize| b.column(i).as_any().downcast_ref::<StringArray>().unwrap().clone();
        let bo = |i: usize| b.column(i).as_any().downcast_ref::<BooleanArray>().unwrap().clone();
        let i32c = |i: usize| b.column(i).as_any().downcast_ref::<Int32Array>().unwrap().clone();
        let i64c = |i: usize| b.column(i).as_any().downcast_ref::<Int64Array>().unwrap().clone();
        let (name, ds, ch, rev, sup, mq, rs, re, qs, qe, rl, intr) =
            (s(0), s(1), s(2), bo(3), bo(4), i32c(5), i64c(6), i64c(7), i32c(8), i32c(9), i32c(10), s(11));
        for i in 0..b.num_rows() {
            out.push(ChimericSegment {
                read_name: name.value(i).to_string(), dataset: ds.value(i).to_string(), chromosome: ch.value(i).to_string(),
                reverse: rev.value(i), supplementary: sup.value(i), mapq: mq.value(i).clamp(0, 255) as u8,
                ref_start: rs.value(i), ref_end: re.value(i), query_start: qs.value(i), query_end: qe.value(i),
                read_length: rl.value(i),
                introns: intr.value(i).split(',').filter_map(|x| x.parse().ok()).collect(),
                placed: false,
                short_class: false,
                pieces: 0,
                excess: 0,
            });
        }
    }
    Ok(out)
}

/// Annotated gene bodies per chromosome, for finding introns that join two genes while the BAM is parsed.
pub struct GeneSpans {
    trees: FxHashMap<String, coitrees::COITree<u32, u32>>,
}

impl GeneSpans {
    pub fn from_tree(tree: &crate::data_structures::GeneIntervalTree) -> Self {
        use coitrees::{COITree, Interval, IntervalTree};
        let mut by_chrom: FxHashMap<String, Vec<Interval<u32>>> = FxHashMap::default();
        for (i, g) in tree.records().enumerate() {
            by_chrom.entry(g.chromosome.clone()).or_default()
                .push(Interval::new(g.start.min(g.end) as i32, g.start.max(g.end) as i32, i as u32));
        }
        Self { trees: by_chrom.into_iter().map(|(c, v)| (c, COITree::new(&v))).collect() }
    }

    /// Genes whose body holds `pos`.
    fn at(&self, chrom: &str, pos: i64) -> Vec<u32> {
        use coitrees::IntervalTree;
        let mut out = Vec::new();
        if let Some(t) = self.trees.get(chrom) { t.query(pos as i32, pos as i32, |iv| out.push(iv.metadata.clone())); }
        out
    }

    /// Whether the two sides of an intron [first, last] each lie in an annotated gene and no gene holds both.
    pub fn joins_genes(&self, chrom: &str, first: i64, last: i64) -> bool {
        let a = self.at(chrom, first - 1);
        if a.is_empty() { return false; }
        let b = self.at(chrom, last + 1);
        !b.is_empty() && !a.iter().any(|x| b.contains(x))
    }
}

#[derive(Debug, Clone)]
pub struct FusionOptions {
    pub min_reads: usize,
    /// Reads needed when neither breakpoint lies on an annotated exon boundary.
    pub min_reads_intraexonic: usize,
    pub min_mapq: u8,
    /// Aligned read bases each segment needs (a part placed on annotated exons, `min_placed`).
    pub min_anchor: i32,
    pub min_placed: i32,
    /// Score above the score a placed part needs that one read of a call made of placed parts alone must have.
    pub placed_margin: i32,
    /// Largest gap or overlap between consecutive segments on the read.
    pub max_query_gap: i32,
    pub cluster_window: i64,
    pub boundary_tol: i64,
    pub readthrough_distance: i64,
    /// A 5' breakpoint this close to an annotated end of its gene marks a concatemer.
    pub concatemer_window: i64,
    /// Breakpoints within this distance of an annotated exon boundary of their gene move onto it.
    pub snap_window: i64,
    /// Opposite-strand partners on one chromosome with breakpoints this close are a fold-back.
    pub foldback_distance: i64,
    /// With both breakpoints on annotated exon boundaries and at least this many reads, a call is not a fold-back.
    pub foldback_min_reads: usize,
    /// Opposite-strand parts further apart than this on the genome are not a fold-back.
    pub foldback_gap: i64,
    /// Supporting reads as a fraction of the less-expressed partner's reads a fusion needs.
    pub min_fraction: f64,
    /// A partner joined to more than this many genes among the candidates is a chimera hub.
    pub max_partners: usize,
    /// A call joining a hub passes with this many times the median reads of the hub's other partner pairs and
    /// `hub_dominance_min_reads` reads.
    pub hub_dominance: f64,
    pub hub_dominance_min_reads: usize,
    /// With two or more datasets, a call supported in only one needs this many reads.
    pub single_dataset_min_reads: usize,
    /// Genome FASTA, for the paralog test: a passing call is a paralog when `paralog_min_cover` of one partner's part
    /// (the genome's sequence along its blocks) lies in the other partner's gene span (`paralog_flank` bp either side)
    /// with `paralog_min_identity`.
    pub genome: Option<std::path::PathBuf>,
    /// A call with a pseudogene partner goes through the same tests as any other (the paralog test among them);
    /// when false it is listed as `pseudogene`.
    pub pseudogene_partners: bool,
    /// BAMs to re-read for unsplit reads through each candidate junction (see `fusion_rescue`); needs `genome`.
    pub rescue: Option<crate::fusion_rescue::RescueOptions>,
    pub paralog_min_cover: f64,
    pub paralog_min_identity: f64,
    pub paralog_flank: i64,
    /// List junctions with one part that no annotated gene names in the read's sense (`unnamed_partner`).
    pub unnamed_partners: bool,
    /// Reads such a call needs when the gene's part ends at an annotated exon boundary, and when it does not.
    pub min_reads_unnamed: usize,
    pub min_reads_unnamed_intraexonic: usize,
}

impl Default for FusionOptions {
    fn default() -> Self {
        Self { min_reads: 2, min_reads_intraexonic: 4, min_mapq: 20, min_anchor: 50, min_placed: 25, placed_margin: 8, max_query_gap: 30,
               cluster_window: 10, boundary_tol: 3, readthrough_distance: 1_000_000, concatemer_window: 50, snap_window: 10, foldback_distance: 100_000,
               foldback_min_reads: 5, foldback_gap: 1000, min_fraction: 0.01, max_partners: 3, hub_dominance: 20.0, hub_dominance_min_reads: 10, single_dataset_min_reads: 5,
               genome: None, pseudogene_partners: true, rescue: None, paralog_min_cover: 0.75, paralog_min_identity: 0.9, paralog_flank: 1000,
               unnamed_partners: true, min_reads_unnamed: 5, min_reads_unnamed_intraexonic: 10 }
    }
}

/// Annotated genes and exon boundaries, for placing breakpoints.
struct GeneIndex {
    /// chromosome -> exon intervals with (gene index)
    exons: FxHashMap<String, COITree<u32, u32>>,
    /// (gene index) -> (id, name, chromosome, plus strand, start, end)
    genes: Vec<(String, String, String, bool, i64, i64)>,
    /// (gene index, position): exon ends in transcript orientation (donors, and annotated 3' ends)
    donors: FxHashSet<(u32, i64)>,
    /// (gene index, position): exon starts in transcript orientation (acceptors, and annotated 5' ends)
    acceptors: FxHashSet<(u32, i64)>,
    /// gene index -> annotated 3' ends
    tes: FxHashMap<u32, Vec<i64>>,
    /// (chromosome, donor, acceptor): annotated introns as exon end and next exon start in transcript order
    introns: FxHashSet<(String, i64, i64)>,
    /// gene index -> protein-coding
    coding: Vec<bool>,
    /// (chromosome, plus strand) -> protein-coding gene spans, by start
    coding_spans: FxHashMap<(String, bool), Vec<(i64, i64)>>,
    /// gene index -> a pseudogene
    pseudo: Vec<bool>,
    /// gene index -> reads annotated to the gene
    expr: Vec<u64>,
    /// gene index -> annotated introns (donor, acceptor) and exons (start, end), in transcript orientation
    gene_introns: FxHashMap<u32, Vec<(i64, i64)>>,
    gene_exons: FxHashMap<u32, Vec<(i64, i64)>>,
    /// chromosome -> gene bodies with (gene index)
    bodies: FxHashMap<String, COITree<u32, u32>>,
    /// gene index -> an annotated gene (not one made from this run's reads, nor a locus)
    known: Vec<bool>,
}

impl GeneIndex {
    fn build(db: &PostDB) -> Self {
        let gidx: FxHashMap<i64, u32> = db.genes.iter().enumerate().map(|(i, g)| (g.gene_id, i as u32)).collect();
        let genes: Vec<(String, String, String, bool, i64, i64)> = db.genes.iter().map(|g| {
            let id = db.annot.gene_id_map.get(&g.gene_id).cloned().unwrap_or_else(|| format!("{}", g.gene_id));
            let name = db.annot.gene_name_map.get(&g.gene_id).cloned().or_else(|| g.gene_name.clone()).unwrap_or_else(|| id.clone());
            (id, name, g.chromosome.clone(), g.strand == "+", g.start.min(g.end), g.start.max(g.end))
        }).collect();
        let edge_ends: FxHashMap<i64, (i64, i64)> = db.edges.iter().map(|e| (e.edge_id, (e.v1, e.v2))).collect();
        let mut per_gene: FxHashMap<u32, (String, Vec<(i64, i64)>)> = FxHashMap::default();
        let mut donors = FxHashSet::default();
        let mut acceptors = FxHashSet::default();
        let mut tes: FxHashMap<u32, Vec<i64>> = FxHashMap::default();
        let mut introns: FxHashSet<(String, i64, i64)> = FxHashSet::default();
        let mut gene_introns_set: FxHashMap<u32, FxHashSet<(i64, i64)>> = FxHashMap::default();
        let mut gene_exons_set: FxHashMap<u32, FxHashSet<(i64, i64)>> = FxHashMap::default();
        for t in db.transcripts.iter().filter(|t| db.annot.known_tx_ids.contains(&t.transcript_id)) {
            let Some(&gi) = gidx.get(&t.gene_id) else { continue };
            let Some(pos) = crate::quant_em::model_positions(db, t, &edge_ends) else { continue };
            let tx_introns: Vec<(i64, i64)> = (1..pos.len().saturating_sub(1)).step_by(2).map(|k| (pos[k], pos[k + 1])).collect();
            for k in (1..pos.len().saturating_sub(1)).step_by(2) {
                introns.insert((t.chromosome.clone(), pos[k], pos[k + 1]));
            }
            gene_introns_set.entry(gi).or_default().extend(tx_introns.iter().copied());
            gene_exons_set.entry(gi).or_default().extend(pos.chunks(2).map(|p| (p[0], p[1])));
            // pos: exon boundaries in 5'->3' order, [start, end] per exon
            let e = per_gene.entry(gi).or_insert_with(|| (t.chromosome.clone(), Vec::new()));
            for p in pos.chunks(2) {
                e.1.push((p[0].min(p[1]), p[0].max(p[1])));
                acceptors.insert((gi, p[0]));
                donors.insert((gi, p[1]));
            }
            if let Some(&last) = pos.last() { tes.entry(gi).or_default().push(last); }
        }
        // each gene's exons merged
        let mut iv: FxHashMap<String, Vec<Interval<u32>>> = FxHashMap::default();
        let mut keys: Vec<u32> = per_gene.keys().copied().collect();
        keys.sort_unstable();
        for gi in keys {
            let (chrom, mut ex) = per_gene.remove(&gi).unwrap();
            ex.sort_unstable();
            let mut merged: Vec<(i64, i64)> = Vec::with_capacity(ex.len());
            for (a, b) in ex {
                match merged.last_mut() {
                    Some(m) if a <= m.1 + 1 => m.1 = m.1.max(b),
                    _ => merged.push((a, b)),
                }
            }
            let v = iv.entry(chrom).or_default();
            for (a, b) in merged { v.push(Interval::new(a as i32, b as i32, gi)); }
        }
        let exons = iv.into_iter().map(|(k, v)| (k, COITree::new(&v))).collect();
        let coding: Vec<bool> = db.genes.iter().map(|g| db.annot.gene_type_map.get(&g.gene_id).is_some_and(|t| t == "protein_coding")).collect();
        let mut coding_spans: FxHashMap<(String, bool), Vec<(i64, i64)>> = FxHashMap::default();
        for (g, _) in genes.iter().zip(&coding).filter(|(_, &c)| c) {
            coding_spans.entry((g.2.clone(), g.3)).or_default().push((g.4, g.5));
        }
        for v in coding_spans.values_mut() { v.sort_unstable(); }
        let pseudo = db.genes.iter().map(|g| db.annot.gene_type_map.get(&g.gene_id).is_some_and(|t| t.contains("pseudogene"))).collect();
        let mut expr = vec![0u64; db.genes.len()];
        for r in &db.observed {
            if let Some(&g) = gidx.get(&r.gene_id) { expr[g as usize] += 1; }
        }
        let gene_introns = gene_introns_set.into_iter().map(|(g, v)| (g, v.into_iter().collect())).collect();
        let gene_exons = gene_exons_set.into_iter().map(|(g, v)| (g, v.into_iter().collect())).collect();
        let mut biv: FxHashMap<String, Vec<Interval<u32>>> = FxHashMap::default();
        for (i, g) in genes.iter().enumerate() { biv.entry(g.2.clone()).or_default().push(Interval::new(g.4 as i32, g.5 as i32, i as u32)); }
        let bodies = biv.into_iter().map(|(k, v)| (k, COITree::new(&v))).collect();
        let known = db.genes.iter().map(|g| db.annot.gene_id_map.contains_key(&g.gene_id)).collect();
        Self { exons, genes, donors, acceptors, tes, introns, coding, coding_spans, pseudo, expr, gene_introns, gene_exons, bodies, known }
    }

    /// Genes whose body holds `pos`, by index.
    fn body_genes(&self, chrom: &str, pos: i64) -> Vec<u32> {
        let mut out = Vec::new();
        if let Some(t) = self.bodies.get(chrom) { t.query(pos as i32, pos as i32, |iv| out.push(iv.metadata)); }
        out.sort_unstable();
        out
    }

    /// Adds partners that are loci, not genes, to the index (`unnamed_partner`): no exons and no exon boundaries.
    fn add_loci(&mut self, loci: Vec<(String, String, String, bool, i64, i64)>) {
        for g in loci {
            self.genes.push(g);
            self.coding.push(false);
            self.pseudo.push(false);
            self.known.push(false);
            self.expr.push(u64::MAX);
        }
    }

    /// Genes whose exons overlap the segment's aligned bases, with the overlap, most first.
    fn genes_of(&self, s: &ChimericSegment) -> Vec<(u32, i64)> {
        let Some(tree) = self.exons.get(&s.chromosome) else { return Vec::new() };
        let mut cover: FxHashMap<u32, i64> = FxHashMap::default();
        for (a, b) in s.blocks() {
            tree.query(a as i32, b as i32, |iv| {
                let ov = (b.min(iv.last as i64) - a.max(iv.first as i64) + 1).max(0);
                *cover.entry(iv.metadata).or_insert(0) += ov;
            });
        }
        let mut v: Vec<(u32, i64)> = cover.into_iter().filter(|&(_, c)| c > 0).collect();
        v.sort_by_key(|&(g, c)| (std::cmp::Reverse(c), g));
        v
    }

    /// The nearest annotated exon boundary of gene g within `tol` of p, if any.
    fn snap(set: &FxHashSet<(u32, i64)>, g: u32, p: i64, tol: i64) -> Option<i64> {
        (0..=tol).flat_map(|d| [p - d, p + d]).find(|&q| set.contains(&(g, q)))
    }

    /// The segment's own splice sites (its blocks' inner boundaries) that are annotated exon boundaries of gene g.
    fn sites_matched(&self, g: u32, s: &ChimericSegment) -> u32 {
        let b = s.blocks();
        let mut n = 0;
        for w in b.windows(2) {
            for p in [w[0].1, w[1].0] {
                if self.donors.contains(&(g, p)) || self.acceptors.contains(&(g, p)) { n += 1; }
            }
        }
        n
    }

    fn near(set: &FxHashSet<(u32, i64)>, g: u32, p: i64, tol: i64) -> bool {
        (-tol..=tol).any(|d| set.contains(&(g, p + d)))
    }

    /// Whether a protein-coding gene of the strand lies wholly between two positions.
    fn coding_between(&self, chrom: &str, plus: bool, lo: i64, hi: i64) -> bool {
        self.coding_spans.get(&(chrom.to_string(), plus)).is_some_and(|v| span_between(v, lo, hi))
    }

    fn annotated_intron(&self, chrom: &str, donor: i64, acceptor: i64, tol: i64) -> bool {
        (-tol..=tol).any(|a| (-tol..=tol).any(|b| self.introns.contains(&(chrom.to_string(), donor + a, acceptor + b))))
    }
}

/// One reported fusion.
#[derive(Debug, Clone)]
pub struct FusionCall {
    pub gene5: u32,
    pub gene3: u32,
    pub bp5: i64,
    pub bp3: i64,
    pub reads: usize,
    /// Of `reads`, those the aligner split at this junction (both parts alignments); the others carry a part placed
    /// on annotated exons, or an unsplit clip that continues into the other partner.
    pub reads_split: usize,
    pub reads_by_dataset: FxHashMap<String, usize>,
    pub boundary5: bool,
    pub boundary3: bool,
    pub category: &'static str,
    pub filter: &'static str,
    /// Up to five supporting reads, for inspection.
    pub example_reads: Vec<String>,
    /// Every supporting read (dataset, name).
    pub read_names: Vec<(String, String)>,
    /// Each partner's part of the fusion transcript: the aligned blocks (genomic order) of the supporting reads'
    /// most common exon chain on that side, ending (5' partner) or starting (3' partner) at the breakpoint, and the
    /// reads carrying that chain.
    pub blocks5: Vec<(i64, i64)>,
    pub blocks3: Vec<(i64, i64)>,
    pub chain_reads5: usize,
    pub chain_reads3: usize,
    /// The call's transcripts (`call_isoforms`) and each supporting read's, with whether it fits it, in `read_names`
    /// order.
    pub isoforms: Vec<FusionIsoform>,
    pub read_isoform: Vec<Option<(u32, bool)>>,
    /// One part is not named by an annotated gene in the read's sense (`unnamed_partner`).
    pub unnamed: bool,
}

/// Reads a call needs: `min_reads` when both breakpoints are annotated exon boundaries, one more when one is,
/// `min_reads_intraexonic` when neither; a call with an unnamed part, `min_reads_unnamed` with a boundary and
/// `min_reads_unnamed_intraexonic` without.
fn reads_needed(opts: &FusionOptions, boundary5: bool, boundary3: bool, unnamed: bool) -> usize {
    if unnamed { if boundary5 || boundary3 { opts.min_reads_unnamed } else { opts.min_reads_unnamed_intraexonic } }
    else if boundary5 && boundary3 { opts.min_reads }
    else if boundary5 || boundary3 { opts.min_reads + 1 }
    else { opts.min_reads_intraexonic }
}

/// The partners of a junction that no pair of annotated genes names in one orientation: (gene of `a`, gene of `b`,
/// whether the read is sense). One part runs along a gene whose exons it covers; the other is named by a gene whose
/// body holds its breakpoint on the transcript's strand, else by a locus: `antisense_<gene>` or an intergenic bin.
#[allow(clippy::too_many_arguments)]
fn unnamed_partner(gi: &GeneIndex, a: &ChimericSegment, b: &ChimericSegment, ca: &[(u32, i64)], cb: &[(u32, i64)], opts: &FusionOptions,
                   locus_of: &mut FxHashMap<String, u32>, loci: &mut Vec<(String, String, String, bool, i64, i64)>) -> Option<(u32, u32, bool)> {
    const BIN: i64 = 100_000;
    let mito = |c: &str| matches!(c, "chrM" | "MT" | "chrMT" | "M");
    if mito(&a.chromosome) || mito(&b.chromosome) { return None; }
    // (boundary, exon bases), the read's sense, whether the gene is the 5' partner, the gene
    let mut best: Option<((u8, i64), bool, bool, u32)> = None;
    for sense in [true, false] {
        let (s5, c5, s3, c3) = if sense { (a, ca, b, cb) } else { (b, cb, a, ca) };
        let (t5, t3) = (s5.reverse != sense, s3.reverse != sense);        // the transcript's strand at each part: plus?
        let bp5 = if t5 { s5.ref_end } else { s5.ref_start };
        let bp3 = if t3 { s3.ref_start } else { s3.ref_end };
        for &(g, cov) in c5.iter().take(8) {
            if gi.genes[g as usize].3 != t5 { continue; }
            let key = (GeneIndex::near(&gi.donors, g, bp5, opts.snap_window) as u8, cov);
            if best.map_or(true, |(k, ..)| key > k) { best = Some((key, sense, true, g)); }
        }
        for &(g, cov) in c3.iter().take(8) {
            if gi.genes[g as usize].3 != t3 { continue; }
            let key = (GeneIndex::near(&gi.acceptors, g, bp3, opts.snap_window) as u8, cov);
            if best.map_or(true, |(k, ..)| key > k) { best = Some((key, sense, false, g)); }
        }
    }
    let (_, sense, named5, g) = best?;
    let (s5, s3) = if sense { (a, b) } else { (b, a) };
    let (t5, t3) = (s5.reverse != sense, s3.reverse != sense);
    let bp5 = if t5 { s5.ref_end } else { s5.ref_start };
    let bp3 = if t3 { s3.ref_start } else { s3.ref_end };
    // the shape of a read-through or of a fold-back is not named
    if s5.chromosome == s3.chromosome {
        if t5 == t3 && ((t5 && bp3 > bp5 && bp3 - bp5 <= opts.readthrough_distance) || (!t5 && bp5 > bp3 && bp5 - bp3 <= opts.readthrough_distance)) { return None; }
        if t5 != t3 && (bp5 - bp3).abs() <= opts.foldback_distance { return None; }
    }
    let (u, cu, tu, bpu) = if named5 { (s3, if sense { cb } else { ca }, t3, bp3) } else { (s5, if sense { ca } else { cb }, t5, bp5) };
    if gi.pseudo[g as usize] || cu.iter().any(|&(h, _)| gi.pseudo[h as usize]) { return None; }
    // annotated genes whose body holds the part's end at the junction
    let at: Vec<u32> = gi.body_genes(&u.chromosome, bpu).into_iter().filter(|&h| gi.known[h as usize]).collect();
    if at.iter().any(|&h| gi.pseudo[h as usize]) { return None; }
    // a gene of the transcript's strand whose body holds the part: protein-coding first, then the shortest
    let pick = |plus: bool| -> Option<u32> {
        at.iter().copied().filter(|&h| gi.genes[h as usize].3 == plus)
            .min_by_key(|&h| (!gi.coding[h as usize], gi.genes[h as usize].5 - gi.genes[h as usize].4, h))
    };
    let partner = if let Some(h) = pick(tu) {
        if h == g { return None; }
        h
    } else {
        let host = cu.first().map(|&(h, _)| h).or_else(|| pick(!tu));
        let (key, entry) = match host {
            Some(h) => {
                let x = &gi.genes[h as usize];
                (format!("as:{h}"), (format!("antisense_{}", x.0), format!("antisense_{}", x.1), x.2.clone(), tu, x.4, x.5))
            }
            None => {
                let lo = bpu / BIN * BIN;
                let name = format!("intergenic_{}_{:.1}Mb", u.chromosome, lo as f64 / 1e6);
                (format!("ig:{}:{}:{}", u.chromosome, tu, lo), (name.clone(), name, u.chromosome.clone(), tu, lo.max(1), lo + BIN))
            }
        };
        let n = (gi.genes.len() + loci.len()) as u32;
        *locus_of.entry(key).or_insert_with(|| { loci.push(entry); n })
    };
    let (g5, g3) = if named5 { (g, partner) } else { (partner, g) };
    Some(if sense { (g5, g3, true) } else { (g3, g5, false) })
}

struct Junction {
    side: (u32, bool, u32, bool), gene5: u32, gene3: u32, bp5: i64, bp3: i64, read: usize, dataset: usize,
    /// Creation index, position of the segment pair along the read, and whether the read is sense.
    jid: usize, pos: u32, sense: bool,
    /// Aligned blocks (genomic order) of the 5' and 3' partner's segments.
    blocks5: Vec<(i64, i64)>, blocks3: Vec<(i64, i64)>,
    /// Both parts are alignments (none placed on annotated exons, none an unsplit read's clip).
    split: bool,
    /// The reads with an unaligned piece where this read's placed part starts (0 without a placed part), and the
    /// placed part's score above the score it needs.
    pieces: u32,
    excess: i32,
    /// One part is not named by an annotated gene in the read's sense (`unnamed_partner`).
    unnamed: bool,
}

/// The most common chain of introns among a side's segments (blocks), its reads, and blocks spanning the chain's
/// widest extent; the side's junction end is set to the breakpoint (`junction_high`: the junction is at the high
/// genomic end of these blocks).
fn consensus_blocks(sides: &[&Vec<(i64, i64)>], breakpoint: i64, junction_high: bool) -> (Vec<(i64, i64)>, usize) {
    let introns = |b: &Vec<(i64, i64)>| -> Vec<(i64, i64)> { b.windows(2).map(|w| (w[0].1, w[1].0)).collect() };
    let mut count: FxHashMap<Vec<(i64, i64)>, (usize, i64, i64)> = FxHashMap::default();
    for b in sides {
        if b.is_empty() { continue; }
        let e = count.entry(introns(b)).or_insert((0, i64::MAX, i64::MIN));
        e.0 += 1; e.1 = e.1.min(b[0].0); e.2 = e.2.max(b[b.len() - 1].1);
    }
    let Some((chain, &(n, lo, hi))) = count.iter()
        .max_by(|a, b| a.1 .0.cmp(&b.1 .0).then(b.0.len().cmp(&a.0.len())).then(b.0.cmp(a.0))) else { return (Vec::new(), 0) };
    let mut blocks = Vec::with_capacity(chain.len() + 1);
    let mut start = lo;
    for &(a, b) in chain { blocks.push((start, a)); start = b; }
    blocks.push((start, hi));
    set_junction_end(&mut blocks, breakpoint, junction_high);
    (blocks, n)
}

/// Whether one of the spans (sorted by start) lies wholly between `lo` and `hi`.
fn span_between(spans: &[(i64, i64)], lo: i64, hi: i64) -> bool {
    let i = spans.partition_point(|&(s, _)| s <= lo);
    spans[i..].iter().take_while(|&&(s, _)| s < hi).any(|&(_, e)| e < hi)
}

/// Bases between two parts on the genome, negative when their spans overlap.
fn parts_apart(a: &[(i64, i64)], b: &[(i64, i64)]) -> i64 {
    match (a.first(), a.last(), b.first(), b.last()) {
        (Some(a0), Some(a1), Some(b0), Some(b1)) => a0.0.max(b0.0) - a1.1.min(b1.1),
        _ => 0,
    }
}

fn set_junction_end(blocks: &mut [(i64, i64)], breakpoint: i64, high: bool) {
    if high { if let Some(l) = blocks.last_mut() { if breakpoint >= l.0 { l.1 = breakpoint; } } }
    else if let Some(f) = blocks.first_mut() { if breakpoint <= f.1 { f.0 = breakpoint; } }
}

/// A part of a fusion transcript from its segments: the consensus chain (`consensus_blocks`) with the ends where the
/// transcript enters (`bp_in`) and leaves (`bp_out`) the part set to those breakpoints; the transcript runs from
/// low to high genomic positions on a plus-strand gene, so it leaves at the high end there.
fn part_blocks(sides: &[&Vec<(i64, i64)>], bp_in: Option<i64>, bp_out: Option<i64>, plus: bool) -> (Vec<(i64, i64)>, usize) {
    let (mut blocks, n) = match (bp_in, bp_out) {
        (_, Some(o)) => consensus_blocks(sides, o, plus),
        (Some(i), None) => consensus_blocks(sides, i, !plus),
        (None, None) => consensus_blocks(sides, i64::MIN, true),
    };
    if let (Some(i), Some(_)) = (bp_in, bp_out) { set_junction_end(&mut blocks, i, !plus); }
    (blocks, n)
}

// Fusion isoforms. A call's reads are grouped by the parts they cover and their introns corrected; within a group,
// chains are chosen by the reads that fit them as 5'-truncated copies, with one isoform per poly(A) site of a chain.

/// Bases an intron's ends may differ by and still be the same intron.
const ISOFORM_INTRON_TOL: i64 = 5;
/// An intron shorter than this is dropped when few of the pieces spanning it carry it.
const ISOFORM_SHORT_INTRON: i64 = 60;
/// Share of the pieces spanning a short intron that must carry it for it to stand.
const ISOFORM_INTRON_SHARE: f64 = 0.2;
/// Reads, and share of its group's reads, a chain needs to be an isoform of its own.
const ISOFORM_MIN_READS: usize = 2;
const ISOFORM_MIN_SHARE: f64 = 0.05;
/// Read 3' ends within this many bp of the next are one poly(A) site; a site with this share of the chain's reads is
/// an isoform of its own.
const ISOFORM_END_GAP: i64 = 50;
const ISOFORM_END_SHARE: f64 = 0.1;
/// Read 3' ends in no site, within this many bp of the next and holding `ISOFORM_SPREAD_SHARE` of the chain's reads,
/// are a site of their own.
const ISOFORM_SPREAD_GAP: i64 = 500;
const ISOFORM_SPREAD_SHARE: f64 = 0.05;
/// Width of a start peak's window, and the share of the reads reaching the first exon a peak needs to be the 5' end.
const ISOFORM_START_WINDOW: i64 = 20;
const ISOFORM_START_SHARE: f64 = 0.1;
/// Bases a truncated read may reach into an intron of the chain it is placed in.
const ISOFORM_END_TOL: i64 = 10;

/// A transcript of a fusion call: the first part it covers (0-based, in the call's transcript order) and its exon
/// blocks on each part it covers from there (genomic order), and the reads assigned to it.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct FusionIsoform {
    pub first_part: usize,
    pub blocks: Vec<Vec<(i64, i64)>>,
    pub reads: usize,
    /// IsoForge's ID of the fusion transcript (`fusion_isoform_id`), set once the calls are final.
    pub id: String,
}

impl FusionIsoform {
    pub fn last_part(&self) -> usize { self.first_part + self.blocks.len() - 1 }
}

/// The ID of the `n`th fusion transcript of a run: ISOFORF and nine digits.
pub fn fusion_isoform_id(n: usize) -> String { format!("ISOFORF{:09}", n) }

/// One read of a fusion call: its pieces on consecutive parts from `first`, each as aligned blocks in genomic order.
struct CallRead { read: usize, first: usize, pieces: Vec<Vec<(i64, i64)>> }

/// Introns of a piece as (exon end, next exon start), genomic order.
fn introns_of(b: &[(i64, i64)]) -> Vec<(i64, i64)> {
    b.windows(2).map(|w| (w[0].1, w[1].0)).collect()
}

/// Intron correction for one part of a call: each intron of the part's pieces maps to the most used intron within
/// `ISOFORM_INTRON_TOL` bp at both ends; a short intron few pieces carry, and an intron of any length only one of
/// three or more pieces spanning it carries, is dropped (the blocks either side join).
struct IntronFix { map: FxHashMap<(i64, i64), (i64, i64)>, drop: FxHashSet<(i64, i64)> }

impl IntronFix {
    fn build(pieces: &[&Vec<(i64, i64)>]) -> Self {
        let mut count: FxHashMap<(i64, i64), usize> = FxHashMap::default();
        for b in pieces {
            let mut seen = FxHashSet::default();
            for x in introns_of(b) { if seen.insert(x) { *count.entry(x).or_insert(0) += 1; } }
        }
        let mut order: Vec<((i64, i64), usize)> = count.iter().map(|(&k, &v)| (k, v)).collect();
        order.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
        let mut by_donor: Vec<(i64, i64)> = count.keys().copied().collect();
        by_donor.sort_unstable();
        let mut map: FxHashMap<(i64, i64), (i64, i64)> = FxHashMap::default();
        for &(rep, _) in &order {
            if map.contains_key(&rep) { continue; }
            let lo = by_donor.partition_point(|x| x.0 < rep.0 - ISOFORM_INTRON_TOL);
            for &x in by_donor[lo..].iter().take_while(|x| x.0 <= rep.0 + ISOFORM_INTRON_TOL) {
                if !map.contains_key(&x) && (x.1 - rep.1).abs() <= ISOFORM_INTRON_TOL { map.insert(x, rep); }
            }
        }
        let mut support: FxHashMap<(i64, i64), usize> = FxHashMap::default();
        for b in pieces {
            let mut seen = FxHashSet::default();
            for x in introns_of(b) { let r = map[&x]; if seen.insert(r) { *support.entry(r).or_insert(0) += 1; } }
        }
        let mut drop = FxHashSet::default();
        for (&r, &s) in &support {
            let spanning = pieces.iter().filter(|b| !b.is_empty() && b[0].0 <= r.0 && b[b.len() - 1].1 >= r.1).count();
            let short = r.1 - r.0 - 1 < ISOFORM_SHORT_INTRON;
            if (short && (s < ISOFORM_MIN_READS || (s as f64) < ISOFORM_INTRON_SHARE * spanning as f64)) || (s == 1 && spanning >= 3) {
                drop.insert(r);
            }
        }
        Self { map, drop }
    }

    /// The piece's blocks with its introns corrected.
    fn apply(&self, b: &[(i64, i64)]) -> Vec<(i64, i64)> {
        if b.is_empty() { return Vec::new(); }
        let mut blocks = Vec::with_capacity(b.len());
        let mut start = b[0].0;
        for x in introns_of(b) {
            let r = self.map.get(&x).copied().unwrap_or(x);
            if self.drop.contains(&r) { continue; }
            if r.0 >= start { blocks.push((start, r.0)); }
            start = r.1;
        }
        let end = b[b.len() - 1].1;
        if end >= start { blocks.push((start, end)); } else if let Some(l) = blocks.last_mut() { l.1 = l.1.max(end); }
        blocks
    }
}

/// A read's corrected piece on one part in transcript order: its introns 5' to 3', and its 5' and 3' ends.
struct OrientedPiece { introns: Vec<(i64, i64)>, five: i64, three: i64 }

fn orient(blocks: &[(i64, i64)], plus: bool) -> OrientedPiece {
    let mut introns = introns_of(blocks);
    let (lo, hi) = (blocks[0].0, blocks[blocks.len() - 1].1);
    if plus { OrientedPiece { introns, five: lo, three: hi } }
    else { introns.reverse(); OrientedPiece { introns, five: hi, three: lo } }
}

/// Whether position `a` lies upstream of `b` (5' of it) on a part's strand.
fn upstream(a: i64, b: i64, plus: bool) -> bool { if plus { a < b } else { a > b } }

/// Whether a read (pieces on the parts `first..`) fits chain `c` (introns per part in transcript order, same parts)
/// as a 5'-truncated copy: on its first part its introns are the chain's last ones and it starts inside the exon they
/// follow; on every later part the introns are the chain's.
fn fits(read: &[OrientedPiece], c: &[Vec<(i64, i64)>], plus: &[bool], first: usize) -> bool {
    let n = read.len();
    for (k, (r, ci)) in read.iter().zip(c).enumerate() {
        let p = plus[first + k];
        if n > 1 && k == 0 {
            if r.introns.len() > ci.len() || ci[ci.len() - r.introns.len()..] != r.introns[..] { return false; }
            // a read starting in a later exon must not reach back into the intron before it
            if r.introns.len() < ci.len() {
                // (exon end, next exon start) in genomic order: the acceptor is the high end on +, the low end on -
                let t = ci[ci.len() - r.introns.len() - 1];
                let acceptor = if p { t.1 } else { t.0 };
                let slack = if p { acceptor - ISOFORM_END_TOL } else { acceptor + ISOFORM_END_TOL };
                if upstream(r.five, slack, p) { return false; }
            }
        } else if r.introns != *ci {
            return false;
        }
    }
    true
}

/// The densest `ISOFORM_START_WINDOW` window of positions, and its most frequent position (a tight site's 3' end).
fn start_peak(mut v: Vec<i64>) -> Option<i64> {
    if v.is_empty() { return None; }
    v.sort_unstable();
    let (mut best, mut at, mut j) = (0usize, 0usize, 0usize);
    for i in 0..v.len() {
        while v[i] - v[j] > ISOFORM_START_WINDOW { j += 1; }
        if i + 1 - j > best { best = i + 1 - j; at = j; }
    }
    let window = &v[at..at + best];
    let mut mode = (window[0], 0usize);
    let mut k = 0;
    while k < window.len() {
        let e = k + window[k..].iter().take_while(|&&x| x == window[k]).count();
        if e - k > mode.1 { mode = (window[k], e - k); }
        k = e;
    }
    Some(mode.0)
}

/// An isoform's 5' end from the starts of its reads that reach its first exon: the most upstream start peak (the most
/// frequent start in an `ISOFORM_START_WINDOW` window) holding `ISOFORM_START_SHARE` of them, or where all but the most
/// upstream tenth start when that is further upstream or there is no peak.
fn five_end(mut v: Vec<i64>, plus: bool) -> Option<i64> {
    if v.is_empty() { return None; }
    if plus { v.sort_unstable(); } else { v.sort_unstable_by(|a, b| b.cmp(a)); }
    let need = ISOFORM_MIN_READS.max((ISOFORM_START_SHARE * v.len() as f64).ceil() as usize);
    let mut j = 0;
    for i in 0..v.len() {
        if j < i { j = i; }
        while j < v.len() && (v[j] - v[i]).abs() <= ISOFORM_START_WINDOW { j += 1; }
        if j - i >= need {
            let w = &v[i..j];
            let (mut best, mut k) = ((w[0], 0usize), 0usize);
            while k < w.len() {
                let e = k + w[k..].iter().take_while(|&&x| x == w[k]).count();
                if e - k > best.1 { best = (w[k], e - k); }
                k = e;
            }
            let q = v[v.len() / 10];
            return Some(if upstream(q, best.0, plus) { q } else { best.0 });
        }
    }
    Some(v[v.len() / 10])
}

/// Isoforms of one fusion call from its reads (see the section comment). `plus` holds each part's strand, `bp_in`
/// and `bp_out` where the call's transcript enters and leaves each part. Returns the isoforms, reads covering all
/// parts first and then by reads, and each read's isoform with whether it fits it as a copy (false: placed with the
/// isoform sharing most of its introns).
fn call_isoforms(reads: &[CallRead], plus: &[bool], bp_in: &[Option<i64>], bp_out: &[Option<i64>])
                 -> (Vec<FusionIsoform>, FxHashMap<usize, (usize, bool)>) {
    let n_parts = plus.len();
    // intron correction per part, over every read of the call
    let fixes: Vec<IntronFix> = (0..n_parts).map(|p| {
        let pieces: Vec<&Vec<(i64, i64)>> = reads.iter()
            .filter(|r| p >= r.first && p < r.first + r.pieces.len())
            .map(|r| &r.pieces[p - r.first]).filter(|b| !b.is_empty()).collect();
        IntronFix::build(&pieces)
    }).collect();
    let oriented: Vec<Option<Vec<OrientedPiece>>> = reads.iter().map(|r| {
        let v: Vec<OrientedPiece> = r.pieces.iter().enumerate().filter_map(|(k, b)| {
            let fixed = fixes[r.first + k].apply(b);
            (!fixed.is_empty()).then(|| orient(&fixed, plus[r.first + k]))
        }).collect();
        (v.len() == r.pieces.len() && v.len() >= 2).then_some(v)
    }).collect();
    // groups by the parts covered
    let mut groups: FxHashMap<(usize, usize), Vec<usize>> = FxHashMap::default();
    for (i, r) in reads.iter().enumerate() {
        if oriented[i].is_some() { groups.entry((r.first, r.first + r.pieces.len() - 1)).or_default().push(i); }
    }
    let mut keys: Vec<(usize, usize)> = groups.keys().copied().collect();
    keys.sort_by_key(|&(f, l)| (std::cmp::Reverse(l - f), f));
    let mut isoforms: Vec<(FusionIsoform, Vec<usize>)> = Vec::new();
    let mut nearest: FxHashSet<usize> = FxHashSet::default();
    for (f, l) in keys {
        let members = &groups[&(f, l)];
        let chain_of = |i: usize| -> Vec<Vec<(i64, i64)>> { oriented[i].as_ref().unwrap().iter().map(|p| p.introns.clone()).collect() };
        // candidate chains: those at least two reads carry (every distinct chain when none does). A chain another
        // contains is a candidate too.
        let mut exact: FxHashMap<Vec<Vec<(i64, i64)>>, usize> = FxHashMap::default();
        for &i in members { *exact.entry(chain_of(i)).or_insert(0) += 1; }
        let mut maximal: Vec<&Vec<Vec<(i64, i64)>>> = exact.iter().filter(|(_, &n)| n >= ISOFORM_MIN_READS).map(|(c, _)| c).collect();
        if maximal.is_empty() { maximal = exact.keys().collect(); }
        maximal.sort();
        let need = ISOFORM_MIN_READS.max((ISOFORM_MIN_SHARE * members.len() as f64).ceil() as usize);
        let mut left: Vec<usize> = members.clone();
        let mut chosen: Vec<(&Vec<Vec<(i64, i64)>>, Vec<usize>)> = Vec::new();
        while !left.is_empty() {
            let best = maximal.iter()
                .map(|c| (*c, left.iter().copied().filter(|&i| fits(oriented[i].as_ref().unwrap(), c, plus, f)).collect::<Vec<_>>()))
                .max_by(|a, b| a.1.len().cmp(&b.1.len())
                    .then(a.0.iter().map(|x| x.len()).sum::<usize>().cmp(&b.0.iter().map(|x| x.len()).sum::<usize>()))
                    .then(b.0.cmp(a.0)));
            let Some((c, fitting)) = best else { break };
            if fitting.is_empty() || (fitting.len() < need && !chosen.is_empty()) { break; }
            let taken: FxHashSet<usize> = fitting.iter().copied().collect();
            left.retain(|i| !taken.contains(i));
            chosen.push((c, fitting));
        }
        // reads fitting no isoform join the one sharing most of their introns
        for &i in &left { nearest.insert(i); }
        for i in left {
            let ci = chain_of(i);
            let score = |c: &Vec<Vec<(i64, i64)>>| -> usize {
                c.iter().zip(&ci).map(|(a, b)| a.iter().filter(|x| b.contains(x)).count()).sum()
            };
            if let Some(k) = (0..chosen.len()).max_by(|&a, &b| score(chosen[a].0).cmp(&score(chosen[b].0))
                .then(chosen[a].1.len().cmp(&chosen[b].1.len()))) {
                chosen[k].1.push(i);
            }
        }
        for (c, mut fitting) in chosen {
            fitting.sort_unstable();
            let lp = l - f;
            let reach_last: Vec<usize> = fitting.iter().copied().filter(|&i| oriented[i].as_ref().unwrap()[lp].introns.len() == c[lp].len()).collect();
            let reach_first: Vec<usize> = fitting.iter().copied().filter(|&i| oriented[i].as_ref().unwrap()[0].introns.len() == c[0].len()).collect();
            // poly(A) sites among the reads that reach the last exon
            let pl = plus[l];
            let mut ends: Vec<(i64, usize)> = reach_last.iter().map(|&i| (oriented[i].as_ref().unwrap()[lp].three, i)).collect();
            ends.sort_unstable();
            let mut clusters: Vec<Vec<(i64, usize)>> = Vec::new();
            for e in ends {
                match clusters.last_mut() {
                    Some(cl) if e.0 - cl[cl.len() - 1].0 <= ISOFORM_END_GAP => cl.push(e),
                    _ => clusters.push(vec![e]),
                }
            }
            let site_need = ISOFORM_MIN_READS.max((ISOFORM_END_SHARE * reach_last.len() as f64).ceil() as usize);
            let mut sites: Vec<Vec<(i64, usize)>> = clusters.iter().filter(|cl| cl.len() >= site_need).cloned().collect();
            // ends in no site, scattered over a region with enough reads, are a site of their own
            let spread_need = ISOFORM_MIN_READS.max((ISOFORM_SPREAD_SHARE * reach_last.len() as f64).ceil() as usize);
            let mut spread: Vec<Vec<(i64, usize)>> = Vec::new();
            for cl in clusters.iter().filter(|cl| cl.len() < site_need) {
                for &e in cl {
                    match spread.last_mut() {
                        Some(g) if e.0 - g[g.len() - 1].0 <= ISOFORM_SPREAD_GAP
                            && !sites.iter().any(|st| st[0].0 > g[g.len() - 1].0 && st[0].0 < e.0) => g.push(e),
                        _ => spread.push(vec![e]),
                    }
                }
            }
            sites.extend(spread.into_iter().filter(|g| g.len() >= spread_need));
            sites.sort_by_key(|s| std::cmp::Reverse(s.len()));
            // a site's 3' end: its most frequent end when its ends are tight (the middle 80 % within
            // `ISOFORM_END_GAP`), else where all but the most downstream tenth of its reads end
            let site_end = |v: &[(i64, usize)]| -> i64 {
                let mut e: Vec<i64> = v.iter().map(|x| x.0).collect();
                if pl { e.sort_unstable(); } else { e.sort_unstable_by(|a, b| b.cmp(a)); }
                let q = |f: f64| e[((e.len() as f64 * f).ceil() as usize).clamp(1, e.len()) - 1];
                if (q(0.9) - q(0.1)).abs() <= ISOFORM_END_GAP { start_peak(e.clone()).unwrap_or(q(0.5)) } else { q(0.9) }
            };
            // one isoform per site; reads of small end clusters go to the nearest site, reads short of the last exon
            // to the site with most reads
            let mut split: Vec<(Option<i64>, Vec<usize>)> = if sites.len() >= 2 {
                // reads go to the site whose span holds their end, else the nearest one
                let spans: Vec<(i64, i64)> = sites.iter().map(|s| (s[0].0, s[s.len() - 1].0)).collect();
                let peaks: Vec<i64> = sites.iter().map(|s| site_end(s)).collect();
                let mut by_site: Vec<Vec<usize>> = vec![Vec::new(); peaks.len()];
                let reached: FxHashMap<usize, i64> = reach_last.iter().map(|&i| (i, oriented[i].as_ref().unwrap()[lp].three)).collect();
                for &i in &fitting {
                    let k = match reached.get(&i) {
                        Some(&e) => (0..spans.len()).min_by_key(|&k| {
                            let (lo, hi) = spans[k];
                            (if e < lo { lo - e } else if e > hi { e - hi } else { 0 }, k)
                        }).unwrap(),
                        None => 0,
                    };
                    by_site[k].push(i);
                }
                peaks.into_iter().map(Some).zip(by_site).collect()
            } else {
                // no site: where all but the most downstream tenth of the reads end
                let end = sites.first().map(|s| site_end(s)).or_else(|| {
                    let mut v: Vec<i64> = fitting.iter().map(|&i| oriented[i].as_ref().unwrap()[lp].three).collect();
                    if pl { v.sort_unstable_by(|a, b| b.cmp(a)); } else { v.sort_unstable(); }
                    v.get(v.len() / 10).copied()
                });
                vec![(end, fitting.clone())]
            };
            for (end3, rs) in split.drain(..) {
                let first_reach: Vec<i64> = rs.iter().copied().filter(|i| reach_first.contains(i))
                    .map(|i| oriented[i].as_ref().unwrap()[0].five).collect();
                let pf = plus[f];
                let five = five_end(first_reach, pf).or_else(|| {
                    five_end(rs.iter().map(|&i| oriented[i].as_ref().unwrap()[0].five).collect(), pf)
                });
                let (Some(five), Some(three)) = (five, end3) else { continue };
                let blocks: Vec<Vec<(i64, i64)>> = (f..=l).map(|p| {
                    let outer5 = if p == f { five } else { bp_in[p].unwrap_or(five) };
                    let outer3 = if p == l { three } else { bp_out[p].unwrap_or(three) };
                    let (lo, hi) = if plus[p] { (outer5, outer3) } else { (outer3, outer5) };
                    let mut introns = c[p - f].clone();
                    introns.sort_unstable();
                    let mut b = Vec::with_capacity(introns.len() + 1);
                    let mut start = lo;
                    for (d, a) in introns {
                        if d >= start && a <= hi { b.push((start, d)); start = a; }
                    }
                    if hi >= start { b.push((start, hi)); }
                    b
                }).collect();
                if blocks.iter().any(|b| b.is_empty()) { continue; }
                isoforms.push((FusionIsoform { first_part: f, blocks, reads: rs.len(), id: String::new() }, rs));
            }
        }
    }
    isoforms.sort_by(|a, b| {
        let full = |x: &FusionIsoform| x.first_part == 0 && x.last_part() == n_parts - 1;
        full(&b.0).cmp(&full(&a.0)).then(b.0.reads.cmp(&a.0.reads)).then(a.0.first_part.cmp(&b.0.first_part))
    });
    let mut of = FxHashMap::default();
    for (k, (_, rs)) in isoforms.iter().enumerate() { for &i in rs { of.insert(reads[i].read, (k, !nearest.contains(&i))); } }
    (isoforms.into_iter().map(|x| x.0).collect(), of)
}

/// A fusion of three or more segments: one molecule joining several loci in order, the hyperedge of which the
/// pairwise calls in `calls` are the consecutive junctions.
#[derive(Debug, Clone)]
pub struct MultiFusionCall {
    /// Pairwise calls (indices into the call list) in transcript order.
    pub calls: Vec<usize>,
    /// Reads spanning every junction.
    pub reads: usize,
    /// Reads of the constituent pairwise calls: the spanning ones and those that cover only part of the path.
    pub compatible_reads: usize,
    pub reads_by_dataset: FxHashMap<String, usize>,
    pub filter: &'static str,
    pub example_reads: Vec<String>,
    /// Every read spanning the whole path (dataset, name).
    pub read_names: Vec<(String, String)>,
    /// Each part (gene, in transcript order) of the fusion transcript from the reads spanning the path: the aligned
    /// blocks of their most common exon chain in it, entering and leaving at the junctions' breakpoints, and the
    /// reads carrying that chain.
    pub blocks: Vec<Vec<(i64, i64)>>,
    pub chain_reads: Vec<usize>,
    /// The transcripts of the path's reads (`call_isoforms`): those spanning it and those over part of it; each
    /// spanning read's isoform in `read_names` order, and the reads over part of the path with theirs.
    pub isoforms: Vec<FusionIsoform>,
    pub read_isoform: Vec<Option<(u32, bool)>>,
    pub partial_reads: Vec<((String, String), Option<(u32, bool)>)>,
}

/// Whether `blocks` of gene `from` (the genome's sequence along them) lie in gene `to` as well: they are aligned
/// (`splice:hq`) to `to`'s span with `paralog_flank` bp either side, their own place masked when the genes overlap,
/// and `paralog_min_cover` of them must align with `paralog_min_identity`.
fn part_in_gene(g: &crate::genome_ref::GenomeRef, gi: &GeneIndex, from: u32, blocks: &[(i64, i64)], to: u32, opts: &FusionOptions) -> bool {
    let chrom_from = &gi.genes[from as usize].2;
    let q: Vec<u8> = blocks.iter().flat_map(|&(a, b)| g.fetch_seq(chrom_from, a - 1, b).unwrap_or_default()).collect();
    if q.len() < 50 { return false; }
    let (_, _, chrom, _, s, e) = &gi.genes[to as usize];
    let lo = (s - opts.paralog_flank).max(1);
    let mut target = g.fetch_seq(chrom, lo - 1, e + opts.paralog_flank).unwrap_or_default();
    if chrom == chrom_from {
        for &(a, b) in blocks {
            for p in a.max(lo)..=b.min(lo + target.len() as i64 - 1) { target[(p - lo) as usize] = b'N'; }
        }
    }
    let Some(ix) = crate::align::SequenceIndex::new(&[target], "splice:hq") else { return false };
    let (covered, matches) = ix.cover(&q);
    covered as f64 >= opts.paralog_min_cover * q.len() as f64 && matches as f64 >= opts.paralog_min_identity * covered as f64
}

/// The passing calls one of whose partners' parts lies in the other partner's gene as well. A part is the genome's
/// sequence along the call's blocks on its partner (`part_in_gene`).
fn paralog_calls(calls: &[FusionCall], gi: &GeneIndex, genome: &Path, opts: &FusionOptions) -> Vec<usize> {
    use rayon::prelude::*;
    let todo: Vec<usize> = (0..calls.len()).filter(|&i| calls[i].filter == "PASS").collect();
    todo.par_iter().map_init(|| crate::genome_ref::GenomeRef::open(genome.to_path_buf()).ok(), |g, &i| {
        let Some(g) = g.as_ref() else { return None };
        let c = &calls[i];
        (part_in_gene(g, gi, c.gene5, &c.blocks5, c.gene3, opts) || part_in_gene(g, gi, c.gene3, &c.blocks3, c.gene5, opts))
            .then_some(i)
    }).flatten().collect()
}

/// The passing calls that copy a stronger call with a paralog in place of one partner: they share the other partner
/// and its breakpoint (`cluster_window`), and the differing partner's part lies in the stronger call's gene. `calls`
/// are in order of reads, most first.
fn paralog_partner_calls(calls: &[FusionCall], gi: &GeneIndex, genome: &Path, opts: &FusionOptions) -> Vec<usize> {
    use rayon::prelude::*;
    let todo: Vec<usize> = (0..calls.len()).filter(|&i| calls[i].filter == "PASS").collect();
    todo.par_iter().map_init(|| crate::genome_ref::GenomeRef::open(genome.to_path_buf()).ok(), |g, &i| {
        let Some(g) = g.as_ref() else { return None };
        let c = &calls[i];
        // a call between genes is compared with calls between genes, a call with an unnamed part with its like
        let copy = (0..i).filter(|&j| matches!(calls[j].filter, "PASS" | "secondary") && calls[j].unnamed == c.unnamed).any(|j| {
            let m = &calls[j];
            (m.gene3 == c.gene3 && m.gene5 != c.gene5 && (m.bp3 - c.bp3).abs() <= opts.cluster_window
                && part_in_gene(g, gi, c.gene5, &c.blocks5, m.gene5, opts))
            || (m.gene5 == c.gene5 && m.gene3 != c.gene3 && (m.bp5 - c.bp5).abs() <= opts.cluster_window
                && part_in_gene(g, gi, c.gene3, &c.blocks3, m.gene3, opts))
        });
        copy.then_some(i)
    }).flatten().collect()
}

/// A candidate junction for re-reading unsplit reads: a cluster of junctions by breakpoints, with its labels.
struct JunctionSeed {
    seed: crate::fusion_rescue::Seed,
    side: (u32, bool, u32, bool),
    gene5: u32,
    gene3: u32,
    unnamed: bool,
}

/// Clusters of junctions with breakpoints within `cluster_window` of the first of the cluster, with the most common
/// breakpoints and genes, and the consensus blocks of each part. The junctions two genes' exons name make the seeds
/// of the first rank, those with an unnamed part (`unnamed`) the second.
fn junction_seeds(junctions: &[Junction], gi: &GeneIndex, opts: &FusionOptions, unnamed: bool) -> Vec<JunctionSeed> {
    let mut idx: Vec<usize> = (0..junctions.len()).filter(|&i| junctions[i].unnamed == unnamed).collect();
    idx.sort_by_key(|&i| { let j = &junctions[i]; (j.side, j.bp5, j.bp3) });
    let mut out = Vec::new();
    let mut i = 0;
    while i < idx.len() {
        let first = &junctions[idx[i]];
        let mut k = i + 1;
        while k < idx.len() {
            let j = &junctions[idx[k]];
            if j.side != first.side || (j.bp5 - first.bp5).abs() > opts.cluster_window || (j.bp3 - first.bp3).abs() > opts.cluster_window { break; }
            k += 1;
        }
        let members: Vec<&Junction> = idx[i..k].iter().map(|&m| &junctions[m]).collect();
        // a cluster makes a seed when the aligner itself split one of its reads at least
        if !members.iter().any(|m| m.split) { i = k; continue; }
        let mode = |f: &dyn Fn(&Junction) -> (i64, i64)| -> (i64, i64) {
            let mut c: FxHashMap<(i64, i64), usize> = FxHashMap::default();
            for m in &members { *c.entry(f(m)).or_insert(0) += 1; }
            c.into_iter().max_by_key(|&(p, n)| (n, std::cmp::Reverse(p))).unwrap().0
        };
        let (bp5, bp3) = mode(&|j: &Junction| (j.bp5, j.bp3));
        let (g5, g3) = mode(&|j: &Junction| (j.gene5 as i64, j.gene3 as i64));
        let (g5, g3) = (g5 as u32, g3 as u32);
        let (x5, x3) = (&gi.genes[g5 as usize], &gi.genes[g3 as usize]);
        let side5: Vec<&Vec<(i64, i64)>> = members.iter().map(|m| &m.blocks5).collect();
        let side3: Vec<&Vec<(i64, i64)>> = members.iter().map(|m| &m.blocks3).collect();
        let (blocks5, _) = consensus_blocks(&side5, bp5, x5.3);
        let (blocks3, _) = consensus_blocks(&side3, bp3, !x3.3);
        if !blocks5.is_empty() && !blocks3.is_empty() {
            out.push(JunctionSeed {
                seed: crate::fusion_rescue::Seed { chrom5: x5.2.clone(), chrom3: x3.2.clone(), plus5: x5.3, plus3: x3.3, bp5, bp3, blocks5, blocks3, second: unnamed },
                side: first.side, gene5: g5, gene3: g3, unnamed,
            });
        }
        i = k;
    }
    out
}

/// Fusions from chimeric segments, sorted by supporting reads.
pub fn call_fusions(db: &PostDB, segs: &[ChimericSegment], opts: FusionOptions) -> (Vec<FusionCall>, Vec<MultiFusionCall>, Vec<String>, Vec<String>) {
    let mut gi = GeneIndex::build(db);
    // partners that are loci (`unnamed_partner`): their gene indexes follow the annotated genes'
    let mut locus_of: FxHashMap<String, u32> = FxHashMap::default();
    let mut loci: Vec<(String, String, String, bool, i64, i64)> = Vec::new();
    let mut n_unnamed = 0usize;
    let pseudo_partners = opts.pseudogene_partners;
    let mut by_read: FxHashMap<(&str, &str), Vec<&ChimericSegment>> = FxHashMap::default();
    for s in segs { by_read.entry((s.dataset.as_str(), s.read_name.as_str())).or_default().push(s); }
    let mut datasets: Vec<String> = segs.iter().map(|s| s.dataset.clone()).collect::<FxHashSet<_>>().into_iter().collect();
    datasets.sort();
    let ds_of: FxHashMap<&str, usize> = datasets.iter().enumerate().map(|(i, d)| (d.as_str(), i)).collect();
    let mut keys: Vec<_> = by_read.keys().copied().collect();
    keys.sort_unstable();
    let mut junctions: Vec<Junction> = Vec::new();
    // junctions with an unnamed part get their chromosomes' numbers after the loop (the numbers order the clusters)
    let mut unnumbered: Vec<(usize, String, String)> = Vec::new();
    let mut chrom_ids: FxHashMap<String, u32> = FxHashMap::default();
    let mut chrom_id = |c: &str| -> u32 {
        let n = chrom_ids.len() as u32;
        *chrom_ids.entry(c.to_string()).or_insert(n)
    };
    for (ri, key) in keys.iter().enumerate() {
        let mut v = by_read[key].clone();
        if v.len() < 2 { continue; }
        v.sort_by_key(|s| (s.query_start, s.query_end));
        let cands: Vec<Vec<(u32, i64)>> = v.iter().map(|s| gi.genes_of(s)).collect();
        let genes: Vec<Option<u32>> = cands.iter().map(|c| c.first().map(|&(g, _)| g)).collect();
        // a segment between two segments of one other gene makes no junction
        let embedded: Vec<bool> = (0..v.len()).map(|i| {
            i > 0 && i + 1 < v.len() && genes[i - 1].is_some() && genes[i - 1] == genes[i + 1] && genes[i] != genes[i - 1]
        }).collect();
        for i in 0..v.len() - 1 {
            let (a, b) = (v[i], v[i + 1]);
            if embedded[i] || embedded[i + 1] { continue; }
            let gap = b.query_start - a.query_end;
            if gap.abs() > opts.max_query_gap { continue; }
            if a.mapq < opts.min_mapq || b.mapq < opts.min_mapq { continue; }
            let short = |s: &ChimericSegment| s.query_end - s.query_start < if s.short_class { crate::fusion_place::SHORT_MIN as i32 } else if s.placed { opts.min_placed } else { opts.min_anchor };
            if short(a) || short(b) { continue; }
            // One gene holding both segments in the read's orientation makes no junction, unless the segments are out
            // of its order, each on another gene, those two not overlapping and both on the strand opposite to it.
            let overlapping = |x: u32, y: u32| {
                let (gx, gy) = (&gi.genes[x as usize], &gi.genes[y as usize]);
                gx.2 == gy.2 && gx.4 <= gy.5 && gy.4 <= gx.5
            };
            let intragenic = cands[i].iter().any(|&(g, _)| {
                let pg = gi.genes[g as usize].3;
                if !cands[i + 1].iter().any(|&(h, _)| h == g) { return false; }
                let sense = a.reverse != pg;
                if sense != (b.reverse != pg) { return false; }
                let (s5, s3) = if sense { (a, b) } else { (b, a) };
                let colinear = if pg { s5.ref_end <= s3.ref_start + opts.snap_window } else { s5.ref_start + opts.snap_window >= s3.ref_end };
                if colinear { return true; }
                let (Some(&(ta, _)), Some(&(tb, _))) = (cands[i].iter().find(|x| x.0 != g), cands[i + 1].iter().find(|x| x.0 != g)) else { return true };
                overlapping(ta, tb) || gi.genes[ta as usize].3 == pg || gi.genes[tb as usize].3 == pg
            });
            if intragenic { continue; }
            // The partners: of the genes the segments overlap, the pair whose strands agree on the read's orientation,
            // by annotated exon boundaries at the breakpoints, then annotated splice sites, coding genes, exon bases.
            let mut best: Option<((u8, u32, u8, i64), (u32, u32, bool))> = None;
            for &(ga, ca) in cands[i].iter().take(8) {
                for &(gb, cb) in cands[i + 1].iter().take(8) {
                    if ga == gb { continue; }
                    let (pa, pb) = (gi.genes[ga as usize].3, gi.genes[gb as usize].3);
                    let (sa, sb) = (a.reverse != pa, b.reverse != pb);
                    if sa != sb { continue; }
                    let (s5, g5, p5, s3, g3, p3) = if sa { (a, ga, pa, b, gb, pb) } else { (b, gb, pb, a, ga, pa) };
                    let bp5 = if p5 { s5.ref_end } else { s5.ref_start };
                    let bp3 = if p3 { s3.ref_start } else { s3.ref_end };
                    let nb = GeneIndex::near(&gi.donors, g5, bp5, opts.snap_window) as u8
                           + GeneIndex::near(&gi.acceptors, g3, bp3, opts.snap_window) as u8;
                    let pc = gi.coding[ga as usize] as u8 + gi.coding[gb as usize] as u8;
                    // then the genes whose exons the parts follow (their splice sites among the gene's annotated ones)
                    let sites = gi.sites_matched(ga, a) + gi.sites_matched(gb, b);
                    let key = (nb, sites, pc, ca + cb);
                    if best.map_or(true, |(k, _)| key > k) { best = Some((key, (ga, gb, sa))); }
                }
            }
            let (ga, gb, sa, unnamed) = match best {
                Some((_, x)) => (x.0, x.1, x.2, false),
                None => {
                    // only parts the aligner aligned are given an unnamed partner
                    if !opts.unnamed_partners || a.placed || b.placed { continue; }
                    let Some(x) = unnamed_partner(&gi, a, b, &cands[i], &cands[i + 1], &opts, &mut locus_of, &mut loci) else { continue };
                    n_unnamed += 1;
                    (x.0, x.1, x.2, true)
                }
            };
            let plus_of = |g: u32| if (g as usize) < gi.genes.len() { gi.genes[g as usize].3 } else { loci[g as usize - gi.genes.len()].3 };
            let (pa, pb) = (plus_of(ga), plus_of(gb));
            let (s5, g5, p5, s3, g3, p3) = if sa { (a, ga, pa, b, gb, pb) } else { (b, gb, pb, a, ga, pa) };
            let bp5 = if p5 { s5.ref_end } else { s5.ref_start };
            let bp3 = if p3 { s3.ref_start } else { s3.ref_end };
            // a breakpoint moves onto an annotated boundary within `snap_window`
            let bp5 = GeneIndex::snap(&gi.donors, g5, bp5, opts.snap_window).unwrap_or(bp5);
            let bp3 = GeneIndex::snap(&gi.acceptors, g3, bp3, opts.snap_window).unwrap_or(bp3);
            // an annotated intron makes no junction
            if s5.chromosome == s3.chromosome && p5 == p3 && gi.annotated_intron(&s5.chromosome, bp5, bp3, opts.boundary_tol) { continue; }
            let jid = junctions.len();
            let side = if unnamed {
                unnumbered.push((jid, s5.chromosome.clone(), s3.chromosome.clone()));
                (u32::MAX, p5, u32::MAX, p3)
            } else {
                (chrom_id(&s5.chromosome), p5, chrom_id(&s3.chromosome), p3)
            };
            junctions.push(Junction { side, gene5: g5, gene3: g3, bp5, bp3, read: ri, dataset: ds_of[s5.dataset.as_str()],
                                      jid, pos: i as u32, sense: sa, blocks5: s5.blocks(), blocks3: s3.blocks(), split: !a.placed && !b.placed, pieces: a.pieces.max(b.pieces), excess: a.excess.max(b.excess), unnamed });
        }
    }
    for (jid, c5, c3) in unnumbered {
        let j = &mut junctions[jid];
        j.side = (chrom_id(&c5), j.side.1, chrom_id(&c3), j.side.3);
    }
    if opts.unnamed_partners {
        log::info!("[Fusion] {} junctions with a part no annotated gene names in the read's sense, at {} loci", n_unnamed, loci.len());
    }
    // genes come first in the index, loci after them
    let n_genes = gi.genes.len() as u32;
    gi.add_loci(std::mem::take(&mut loci));
    let gi = gi;
    // reads the aligner did not split: a soft clip at a candidate's breakpoint that continues into the other part
    let mut read_keys: Vec<(String, String)> = keys.iter().map(|&(d, n)| (d.to_string(), n.to_string())).collect();
    if let (Some(ro), Some(genome)) = (opts.rescue.as_ref(), opts.genome.as_ref()) {
        let mut seeds = junction_seeds(&junctions, &gi, &opts, false);
        if opts.unnamed_partners { seeds.extend(junction_seeds(&junctions, &gi, &opts, true)); }
        let plain: Vec<crate::fusion_rescue::Seed> = seeds.iter().map(|s| s.seed.clone()).collect();
        // a read already split, or with a part placed on annotated exons, is skipped
        let found: Vec<_> = crate::fusion_rescue::rescue(&plain, ro, genome).into_iter()
            .filter(|r| !by_read.contains_key(&(r.dataset.as_str(), r.name.as_str()))).collect();
        let n_found = found.len();
        for r in found {
            let Some(&di) = ds_of.get(r.dataset.as_str()) else { continue };
            let s = &seeds[r.seed];
            let ri = read_keys.len();
            read_keys.push((r.dataset, r.name));
            let jid = junctions.len();
            let (blocks5, blocks3) = if r.on5 { (r.blocks, s.seed.blocks3.clone()) } else { (s.seed.blocks5.clone(), r.blocks) };
            junctions.push(Junction { side: s.side, gene5: s.gene5, gene3: s.gene3, bp5: s.seed.bp5, bp3: s.seed.bp3, read: ri,
                                      dataset: di, jid, pos: 0, sense: true, blocks5, blocks3, split: false, pieces: 0, excess: 0, unnamed: s.unnamed });
        }
        log::info!("[Fusion] {} unsplit reads continue into the other part at {} candidate junctions", n_found, seeds.len());
    }
    // cluster by breakpoints (the label is the pair most reads chose)
    let t_cluster = std::time::Instant::now();
    let mut call_of = vec![usize::MAX; junctions.len()];
    junctions.sort_by_key(|j| (j.side, j.bp5, j.bp3, j.read));
    let mut calls = Vec::new();
    let mut i = 0;
    while i < junctions.len() {
        let mut k = i;
        while k < junctions.len() && junctions[k].side == junctions[i].side { k += 1; }
        // single linkage on bp5, then split on bp3
        let group = &junctions[i..k];
        // the group is sorted by bp5, so the junctions within `cluster_window` of a member's bp5 are one run of it
        let bp5s: Vec<i64> = group.iter().map(|j| j.bp5).collect();
        let mut used = vec![false; group.len()];
        for s in 0..group.len() {
            if used[s] { continue; }
            let mut members = vec![s]; used[s] = true;
            let mut idx = 0;
            while idx < members.len() {
                let m = &group[members[idx]];
                let lo = bp5s.partition_point(|&b| b < m.bp5 - opts.cluster_window);
                let hi = bp5s.partition_point(|&b| b <= m.bp5 + opts.cluster_window);
                for t in lo..hi {
                    if !used[t] && (group[t].bp3 - m.bp3).abs() <= opts.cluster_window {
                        used[t] = true; members.push(t);
                    }
                }
                idx += 1;
            }
            let mut reads: std::collections::BTreeSet<usize> = std::collections::BTreeSet::new();
            let mut split_reads: FxHashSet<usize> = FxHashSet::default();
            let mut by_ds: FxHashMap<String, FxHashSet<usize>> = FxHashMap::default();
            let mut modes: FxHashMap<(i64, i64), usize> = FxHashMap::default();
            // per pair of partners: the junctions two genes' exons name, and all of them
            let mut labels: FxHashMap<(u32, u32), (usize, usize)> = FxHashMap::default();
            let mut named_reads: FxHashSet<usize> = FxHashSet::default();
            // Reads with a placed part count when they are a tenth at least of the reads with an unaligned piece at
            // that end and, with no read the aligner split, one is placed `placed_margin` above the score needed.
            let placed_reads: FxHashSet<usize> = members.iter().map(|&m| &group[m]).filter(|j| j.pieces > 0).map(|j| j.read).collect();
            let any_split = members.iter().any(|&m| group[m].split);
            let pieces = members.iter().map(|&m| group[m].pieces).max().unwrap_or(0);
            let excess = members.iter().map(|&m| &group[m]).filter(|j| j.pieces > 0).map(|j| j.excess).max().unwrap_or(0);
            let count_placed = 10 * placed_reads.len() >= pieces as usize && (any_split || excess >= opts.placed_margin);
            for &m in &members {
                let j = &group[m];
                call_of[j.jid] = calls.len();
                *modes.entry((j.bp5, j.bp3)).or_insert(0) += 1;
                let l = labels.entry((j.gene5, j.gene3)).or_insert((0, 0));
                l.1 += 1;
                if !j.unnamed { l.0 += 1; }
                if j.pieces > 0 && !count_placed { continue; }
                reads.insert(j.read);
                if !j.unnamed { named_reads.insert(j.read); }
                if j.split { split_reads.insert(j.read); }
                by_ds.entry(datasets[j.dataset].clone()).or_default().insert(j.read);
            }
            let (&(bp5, bp3), _) = modes.iter().max_by_key(|&(&p, &n)| (n, std::cmp::Reverse(p))).unwrap();
            // the label: the pair most reads chose among those two genes' exons name; without one, among all
            let (&(g5, g3), _) = labels.iter().max_by_key(|&(&g, &(named, n))| (named > 0, if named > 0 { named } else { n }, std::cmp::Reverse(g))).unwrap();
            let (ref x5, ref x3) = (&gi.genes[g5 as usize], &gi.genes[g3 as usize]);
            let boundary5 = GeneIndex::near(&gi.donors, g5, bp5, opts.boundary_tol);
            let boundary3 = GeneIndex::near(&gi.acceptors, g3, bp3, opts.boundary_tol);
            let side5: Vec<&Vec<(i64, i64)>> = members.iter().map(|&m| &group[m].blocks5).collect();
            let side3: Vec<&Vec<(i64, i64)>> = members.iter().map(|&m| &group[m].blocks3).collect();
            // the 5' partner's part ends at the breakpoint: its high end on +, its low end on -
            let (blocks5, chain_reads5) = consensus_blocks(&side5, bp5, x5.3);
            let (blocks3, chain_reads3) = consensus_blocks(&side3, bp3, !x3.3);
            let apart = parts_apart(&blocks5, &blocks3);
            // read-through: same strand, 3' partner downstream within `readthrough_distance`, no coding gene between
            let category = if x5.2 != x3.2 {
                "inter-chromosomal"
            } else if x5.3 != x3.3 && (bp5 - bp3).abs() <= opts.foldback_distance
                      && !(boundary5 && boundary3 && reads.len() >= opts.foldback_min_reads)
                      && apart <= opts.foldback_gap {
                "fold-back"
            } else if x5.3 == x3.3 && ((x5.3 && bp3 > bp5 && bp3 - bp5 <= opts.readthrough_distance)
                                        || (!x5.3 && bp5 > bp3 && bp5 - bp3 <= opts.readthrough_distance))
                      && !gi.coding_between(&x5.2, x5.3, bp5.min(bp3), bp5.max(bp3)) {
                "read-through"
            } else {
                "intra-chromosomal"
            };
            // concatemer: a 5' breakpoint within `concatemer_window` of an annotated 3' end of its gene and at no donor
            let concatemer = !boundary5 && gi.tes.get(&g5).is_some_and(|v| v.iter().any(|&t| (t - bp5).abs() <= opts.concatemer_window));
            // a call is two genes' when the reads their exons name reach those calls' floor by themselves; the reads of
            // the same junction that no gene names in their sense (`unnamed_partner`) then count with them. Otherwise the
            // call is held to the floor of calls with an unnamed part, on all its reads.
            let unnamed = members.iter().any(|&m| group[m].unnamed) && named_reads.len() < reads_needed(&opts, boundary5, boundary3, false);
            let need = reads_needed(&opts, boundary5, boundary3, unnamed);
            let mito = |c: &str| matches!(c, "chrM" | "MT" | "chrMT" | "M");
            let n = reads.len() as f64;
            let low_expr = gi.expr[g5 as usize].min(gi.expr[g3 as usize]) as f64;
            let fraction = n / (low_expr + n);
            let filter = if concatemer { "concatemer" }
                else if !pseudo_partners && (gi.pseudo[g5 as usize] || gi.pseudo[g3 as usize]) { "pseudogene" }
                else if mito(&x5.2) || mito(&x3.2) { "mitochondrial" }
                else if category == "read-through" || category == "fold-back" { category }
                else if reads.len() < need { "low_support" }
                else if fraction < opts.min_fraction { "background" }
                else { "PASS" };
            let example_reads = reads.iter().take(5).map(|&r| read_keys[r].1.clone()).collect();
            let read_names = reads.iter().map(|&r| (read_keys[r].0.clone(), read_keys[r].1.clone())).collect();
            calls.push(FusionCall { gene5: g5, gene3: g3, bp5, bp3, reads: reads.len(), reads_split: split_reads.len(), read_names, blocks5, blocks3, chain_reads5, chain_reads3,
                                    isoforms: Vec::new(), read_isoform: Vec::new(),
                                    reads_by_dataset: by_ds.into_iter().map(|(d, r)| (d, r.len())).collect(),
                                    boundary5, boundary3, category, filter, example_reads, unnamed });
        }
        i = k;
    }
    let s_cluster = t_cluster.elapsed().as_secs_f64();
    let t_merge = std::time::Instant::now();
    // The status of a call from its own fields, for calls changed below.
    let status_of = |c: &FusionCall| -> &'static str {
        let (x5, x3) = (&gi.genes[c.gene5 as usize], &gi.genes[c.gene3 as usize]);
        let concatemer = !c.boundary5 && gi.tes.get(&c.gene5).is_some_and(|v| v.iter().any(|&t| (t - c.bp5).abs() <= opts.concatemer_window));
        let need = reads_needed(&opts, c.boundary5, c.boundary3, c.unnamed);
        let mito = |ch: &str| matches!(ch, "chrM" | "MT" | "chrMT" | "M");
        let n = c.reads as f64;
        let low_expr = gi.expr[c.gene5 as usize].min(gi.expr[c.gene3 as usize]) as f64;
        if concatemer { "concatemer" }
            else if !pseudo_partners && (gi.pseudo[c.gene5 as usize] || gi.pseudo[c.gene3 as usize]) { "pseudogene" }
            else if mito(&x5.2) || mito(&x3.2) { "mitochondrial" }
            else if c.category == "read-through" || c.category == "fold-back" { c.category }
            else if c.reads < need { "low_support" }
            else if n / (low_expr + n) < opts.min_fraction { "background" }
            else { "PASS" }
    };
    // A cluster of a gene pair with one breakpoint at a donor (acceptor) joins the pair's cluster that shares its other
    // breakpoint and has this one in the next (previous) exon, at most `max_query_gap` + `boundary_tol` bases in.
    {
        let in_exon = |g: u32, first: i64, p: i64, from_start: bool| -> bool {
            gi.gene_exons.get(&g).is_some_and(|ex| ex.iter().any(|&(st, en)| {
                (if from_start { st } else { en } - first).abs() <= opts.boundary_tol && p >= st.min(en) && p <= st.max(en)
                    && (p - first).abs() + 1 <= opts.max_query_gap as i64 + opts.boundary_tol
            }))
        };
        let exon_after = |g: u32, donor: i64, p: i64| gi.gene_introns.get(&g)
            .is_some_and(|v| v.iter().any(|&(d, acc)| (d - donor).abs() <= opts.boundary_tol && in_exon(g, acc, p, true)));
        let exon_before = |g: u32, acceptor: i64, p: i64| gi.gene_introns.get(&g)
            .is_some_and(|v| v.iter().any(|&(d, acc)| (acc - acceptor).abs() <= opts.boundary_tol && in_exon(g, d, p, false)));
        let mut target: Vec<Option<usize>> = vec![None; calls.len()];
        // only the calls of one gene pair are compared, each pair's in the order of the calls
        let mut of_pair: FxHashMap<(u32, u32, bool), Vec<usize>> = FxHashMap::default();
        for (x, c) in calls.iter().enumerate() { of_pair.entry((c.gene5, c.gene3, c.unnamed)).or_default().push(x); }
        for x in 0..calls.len() {
            let cx = &calls[x];
            if cx.reads_split == 0 { continue; }
            for &y in &of_pair[&(cx.gene5, cx.gene3, cx.unnamed)] {
                if x == y { continue; }
                let cy = &calls[y];
                let same5 = (cx.bp5 - cy.bp5).abs() <= opts.cluster_window;
                let same3 = (cx.bp3 - cy.bp3).abs() <= opts.cluster_window;
                if cy.reads_split == 0 { continue; }
                if (same3 && !same5 && cx.boundary5 && exon_after(cx.gene5, cx.bp5, cy.bp5))
                    || (same5 && !same3 && cx.boundary3 && exon_before(cx.gene3, cx.bp3, cy.bp3)) { target[x] = Some(y); break; }
            }
        }
        // holders[u]: the calls whose junctions now belong to call u (each call's at first, then those that moved in)
        let mut holders: Vec<Vec<usize>> = (0..calls.len()).map(|v| vec![v]).collect();
        for x in 0..calls.len() {
            let Some(mut y) = target[x] else { continue };
            let mut seen = vec![x];
            while let Some(z) = target[y] { if seen.contains(&z) { break; } seen.push(y); y = z; }
            let moved: Vec<(String, String)> = calls[x].read_names.clone();
            let split = calls[x].reads_split;
            for (d, n) in moved {
                if calls[y].read_names.iter().any(|r| r.0 == d && r.1 == n) { continue; }
                *calls[y].reads_by_dataset.entry(d.clone()).or_insert(0) += 1;
                calls[y].read_names.push((d, n));
            }
            calls[y].reads = calls[y].read_names.len();
            calls[y].reads_split = (calls[y].reads_split + split).min(calls[y].reads);
            calls[y].filter = status_of(&calls[y]);
            calls[x].filter = "secondary";
            let moved_calls = std::mem::take(&mut holders[x]);
            holders[y].extend(moved_calls);
        }
        let mut now_of = vec![0usize; calls.len()];
        for (u, hs) in holders.iter().enumerate() { for &v in hs { now_of[v] = u; } }
        for c in call_of.iter_mut() { if *c != usize::MAX { *c = now_of[*c]; } }
    }
    let s_merge = t_merge.elapsed().as_secs_f64();
    let t_mid = std::time::Instant::now();
    let mut order: Vec<usize> = (0..calls.len()).collect();
    order.sort_by(|&x, &y| { let (a, b) = (&calls[x], &calls[y]);
        b.reads.cmp(&a.reads).then(a.gene5.cmp(&b.gene5)).then(a.gene3.cmp(&b.gene3)).then(a.bp5.cmp(&b.bp5)) });
    let mut new_index = vec![0usize; calls.len()];
    for (k, &o) in order.iter().enumerate() { new_index[o] = k; }
    let mut calls: Vec<FusionCall> = order.iter().map(|&o| calls[o].clone()).collect();
    for c in call_of.iter_mut() { if *c != usize::MAX { *c = new_index[*c]; } }
    // Paralogs: one partner's part of the reads lies in the other partner's gene too.
    if let Some(genome) = &opts.genome {
        for i in paralog_calls(&calls, &gi, genome, &opts) { calls[i].filter = "paralog"; }
    }
    // Chimera hubs: a partner joined to more than `max_partners` genes among the candidates with `min_reads` reads,
    // structural artefacts and opposite-strand neighbours within `foldback_distance` left out.
    let structural = |f: &str| matches!(f, "concatemer" | "pseudogene" | "mitochondrial" | "fold-back" | "paralog");
    let local_inversion = |c: &FusionCall| {
        let (x5, x3) = (&gi.genes[c.gene5 as usize], &gi.genes[c.gene3 as usize]);
        x5.2 == x3.2 && x5.3 != x3.3 && (c.bp5 - c.bp3).abs() <= opts.foldback_distance
    };
    let counted = |c: &FusionCall| c.reads >= opts.min_reads && !structural(c.filter) && !local_inversion(c);
    // A call with an unnamed part (`unnamed_partner`) does not count among a gene's partners. A gene's unnamed
    // partners are counted apart.
    let is_locus = |g: u32| g >= n_genes;
    let mut partners: FxHashMap<u32, FxHashSet<u32>> = FxHashMap::default();
    let mut unnamed_of: FxHashMap<u32, Vec<(u32, usize)>> = FxHashMap::default();
    for c in calls.iter().filter(|c| counted(c)) {
        if !c.unnamed {
            partners.entry(c.gene5).or_default().insert(c.gene3);
            partners.entry(c.gene3).or_default().insert(c.gene5);
        } else {
            for (g, other) in [(c.gene5, c.gene3), (c.gene3, c.gene5)] {
                if is_locus(g) { continue; }
                let v = unnamed_of.entry(g).or_default();
                match v.iter_mut().find(|x| x.0 == other) { Some(x) => x.1 += c.reads, None => v.push((other, c.reads)) }
            }
        }
    }
    // A call dominates a hub with `hub_dominance` times the median reads of the hub's other pairs. A call with both
    // breakpoints on annotated exon boundaries counts the reads of the pair's other such calls with its own.
    let mut pair_spliced: FxHashMap<(u32, u32), usize> = FxHashMap::default();
    for c in calls.iter().filter(|c| counted(c) && !c.unnamed && c.boundary5 && c.boundary3) { *pair_spliced.entry((c.gene5, c.gene3)).or_insert(0) += c.reads; }
    let mut pair_reads: FxHashMap<(u32, u32), usize> = FxHashMap::default();
    for c in calls.iter().filter(|c| counted(c) && !c.unnamed) { *pair_reads.entry((c.gene5, c.gene3)).or_insert(0) += c.reads; }
    let mut partner_pairs: FxHashMap<u32, Vec<((u32, u32), usize)>> = FxHashMap::default();
    for (&p, &n) in &pair_reads {
        partner_pairs.entry(p.0).or_default().push((p, n));
        partner_pairs.entry(p.1).or_default().push((p, n));
    }
    let dominates = |c: &FusionCall, g: u32| -> bool {
        let p = (c.gene5, c.gene3);
        let reads = if !c.unnamed && c.boundary5 && c.boundary3 { pair_spliced.get(&p).copied().unwrap_or(c.reads) } else { c.reads };
        let mut others: Vec<usize> = partner_pairs.get(&g).map_or(Vec::new(), |v| v.iter().filter(|x| x.0 != p).map(|x| x.1).collect());
        if others.is_empty() { return true; }
        others.sort_unstable();
        let n = others.len();
        let median = if n % 2 == 1 { others[n / 2] as f64 } else { (others[n / 2 - 1] + others[n / 2]) as f64 / 2.0 };
        reads >= opts.hub_dominance_min_reads && reads as f64 >= opts.hub_dominance * median
    };
    for i in 0..calls.len() {
        if calls[i].filter != "PASS" { continue; }
        let deg = |g: u32| partners.get(&g).map_or(0, |v| v.len());
        let (g5, g3) = (calls[i].gene5, calls[i].gene3);
        let hub5 = deg(g5) > opts.max_partners;
        let hub3 = deg(g3) > opts.max_partners;
        if (hub5 && !dominates(&calls[i], g5)) || (hub3 && !dominates(&calls[i], g3)) { calls[i].filter = "promiscuous"; continue; }
        if calls[i].unnamed {
            // a gene joined to more than `max_partners` unnamed partners: the call needs `hub_dominance` times the
            // median reads of the others
            let mut scattered = false;
            for (g, other) in [(g5, g3), (g3, g5)] {
                if is_locus(g) { continue; }
                let v = unnamed_of.get(&g).map(|v| v.as_slice()).unwrap_or(&[]);
                if v.len() <= opts.max_partners { continue; }
                let mine = v.iter().find(|x| x.0 == other).map_or(calls[i].reads, |x| x.1);
                let mut others: Vec<usize> = v.iter().filter(|x| x.0 != other).map(|x| x.1).collect();
                others.sort_unstable();
                let n = others.len();
                let median = if n % 2 == 1 { others[n / 2] as f64 } else { (others[n / 2 - 1] + others[n / 2]) as f64 / 2.0 };
                if !(mine >= opts.hub_dominance_min_reads && mine as f64 >= opts.hub_dominance * median) { scattered = true; }
            }
            if scattered { calls[i].filter = "promiscuous"; continue; }
        }
        let c = &mut calls[i];
        let n_ds = c.reads_by_dataset.values().filter(|&&n| n > 0).count();
        if datasets.len() >= 2 && n_ds == 1 && c.reads < opts.single_dataset_min_reads { c.filter = "single_sample"; }
    }
    // A breakpoint of a gene pair that shares one end with a call of four times its reads is listed as secondary.
    for i in 0..calls.len() {
        if calls[i].filter != "PASS" { continue; }
        let minor = (0..i).any(|j| {
            let (m, c) = (&calls[j], &calls[i]);
            m.gene5 == c.gene5 && m.gene3 == c.gene3 && m.unnamed == c.unnamed && m.reads >= 4 * c.reads
                && ((m.bp5 - c.bp5).abs() <= opts.cluster_window || (m.bp3 - c.bp3).abs() <= opts.cluster_window)
        });
        if minor { calls[i].filter = "secondary"; }
    }
    // A copy of a stronger call with a paralog for one partner is listed as secondary too.
    if let Some(genome) = &opts.genome {
        for i in paralog_partner_calls(&calls, &gi, genome, &opts) { calls[i].filter = "secondary"; }
    }
    // Each pairwise call's transcripts from its reads, and the call's parts from the one with most reads.
    let mut junctions_of: Vec<Vec<usize>> = vec![Vec::new(); calls.len()];
    for (x, j) in junctions.iter().enumerate() {
        if call_of[j.jid] != usize::MAX { junctions_of[call_of[j.jid]].push(x); }
    }
    for (ci, c) in calls.iter_mut().enumerate() {
        let mut seen = FxHashSet::default();
        let reads: Vec<CallRead> = junctions_of[ci].iter().filter(|&&x| seen.insert(junctions[x].read))
            .map(|&x| CallRead { read: junctions[x].read, first: 0, pieces: vec![junctions[x].blocks5.clone(), junctions[x].blocks3.clone()] })
            .collect();
        let plus = [gi.genes[c.gene5 as usize].3, gi.genes[c.gene3 as usize].3];
        let (isoforms, of) = call_isoforms(&reads, &plus, &[None, Some(c.bp3)], &[Some(c.bp5), None]);
        let mut order: Vec<usize> = reads.iter().map(|r| r.read).collect();
        order.sort_unstable();
        c.read_isoform = order.iter().map(|r| of.get(r).map(|&(k, fit)| (k as u32, fit))).collect();
        if let Some(top) = isoforms.first().filter(|t| t.first_part == 0 && t.blocks.len() == 2) {
            c.blocks5 = top.blocks[0].clone();
            c.blocks3 = top.blocks[1].clone();
            c.chain_reads5 = top.reads;
            c.chain_reads3 = top.reads;
        }
        c.isoforms = isoforms;
    }
    let s_mid = t_mid.elapsed().as_secs_f64();
    let t_paths = std::time::Instant::now();
    // Multi-segment fusions: a read whose consecutive segment pairs all made junctions is one hyperedge over their
    // loci, in transcript order (read order for a sense read, reversed for an antisense one).
    let mut per_read: FxHashMap<usize, Vec<(u32, bool, usize, usize)>> = FxHashMap::default();
    // (a junction with an unnamed part ends a path)
    for (x, j) in junctions.iter().enumerate() {
        if call_of[j.jid] != usize::MAX && !calls[call_of[j.jid]].unnamed { per_read.entry(j.read).or_default().push((j.pos, j.sense, call_of[j.jid], x)); }
    }
    // path -> spanning reads, those per dataset, and each spanning read's junctions in path order
    type PathReads = (std::collections::BTreeSet<usize>, FxHashMap<usize, FxHashSet<usize>>, Vec<(usize, Vec<usize>)>);
    let mut paths: FxHashMap<Vec<usize>, PathReads> = FxHashMap::default();
    let mut read_ds: FxHashMap<usize, usize> = FxHashMap::default();
    for j in &junctions { read_ds.insert(j.read, j.dataset); }
    for (&r, v) in per_read.iter_mut() {
        if v.len() < 2 { continue; }
        v.sort_unstable();
        let mut start = 0;
        for e in 1..=v.len() {
            let run_ends = e == v.len() || v[e].0 != v[e - 1].0 + 1 || v[e].1 != v[e - 1].1;
            if !run_ends { continue; }
            if e - start >= 2 {
                let mut path: Vec<usize> = v[start..e].iter().map(|x| x.2).collect();
                let mut js: Vec<usize> = v[start..e].iter().map(|x| x.3).collect();
                if !v[start].1 { path.reverse(); js.reverse(); }
                let entry = paths.entry(path).or_default();
                entry.0.insert(r);
                entry.1.entry(read_ds[&r]).or_default().insert(r);
                entry.2.push((r, js));
            }
            start = e;
        }
    }
    let reads_of_call: Vec<Vec<usize>> = junctions_of.iter().map(|xs| {
        let mut v: Vec<usize> = xs.iter().map(|&x| junctions[x].read).collect();
        v.sort_unstable();
        v.dedup();
        v
    }).collect();
    let artefact = |f: &str| matches!(f, "concatemer" | "pseudogene" | "mitochondrial" | "fold-back" | "paralog" | "background");
    // reads a junction needs by itself, as for a pairwise call: fewer when its breakpoints are exon boundaries
    let need_of = |c: &FusionCall| reads_needed(&opts, c.boundary5, c.boundary3, c.unnamed);
    // A concatemer-status junction continues a path with a passing junction when it has the reads a pairwise call
    // needs, reaches `min_fraction`, and has `min_reads` split reads or its 3' breakpoint on an exon boundary.
    let continues: Vec<bool> = calls.iter().map(|c| {
        let n = c.reads as f64;
        let low_expr = gi.expr[c.gene5 as usize].min(gi.expr[c.gene3 as usize]) as f64;
        c.filter == "concatemer" && c.reads >= need_of(c) && (c.boundary3 || c.reads_split >= opts.min_reads)
            && n / (low_expr + n) >= opts.min_fraction
    }).collect();
    let mut multi: Vec<MultiFusionCall> = Vec::new();
    let mut keys_sorted: Vec<&Vec<usize>> = paths.keys().collect();
    keys_sorted.sort();
    for path in keys_sorted {
        let (reads, by_ds, spans) = &paths[path];
        let compatible: FxHashSet<usize> = path.iter().flat_map(|&c| reads_of_call[c].iter().copied()).collect();
        let n_ds = by_ds.values().filter(|r| !r.is_empty()).count();
        // with a read over the whole path, the path is supported when each junction has the reads a pairwise call needs
        let each_junction = !reads.is_empty() && path.iter().all(|&c| calls[c].reads >= need_of(&calls[c]));
        let need = path.iter().map(|&c| need_of(&calls[c])).max().unwrap_or(opts.min_reads);
        let carried = |c: usize| continues[c] && path.iter().any(|&d| d != c && calls[d].filter == "PASS");
        let filter = if let Some(c) = path.iter().find(|&&c| artefact(calls[c].filter) && !carried(c)) { calls[*c].filter }
            else if reads.len() < need && !each_junction { "low_support" }
            else if datasets.len() >= 2 && n_ds == 1 && reads.len() < opts.single_dataset_min_reads { "single_sample" }
            else { "PASS" };
        // part p is the 5' segment of junction p; the last part the 3' segment of the last junction
        // (a read carrying the path twice counts once)
        let mut spans: Vec<&(usize, Vec<usize>)> = spans.iter().collect();
        spans.sort_unstable();
        spans.dedup_by_key(|x| x.0);
        let k = path.len();
        let (mut blocks, mut chain_reads) = (Vec::with_capacity(k + 1), Vec::with_capacity(k + 1));
        for p in 0..=k {
            let sides: Vec<&Vec<(i64, i64)>> = spans.iter()
                .map(|(_, js)| if p < k { &junctions[js[p]].blocks5 } else { &junctions[js[k - 1]].blocks3 }).collect();
            let g = if p == 0 { calls[path[0]].gene5 } else { calls[path[p - 1]].gene3 };
            let bp_in = (p > 0).then(|| calls[path[p - 1]].bp3);
            let bp_out = (p < k).then(|| calls[path[p]].bp5);
            let (b, n) = part_blocks(&sides, bp_in, bp_out, gi.genes[g as usize].3);
            blocks.push(b);
            chain_reads.push(n);
        }
        // the path's transcripts: from the reads spanning it, and from those over a run of its junctions
        let part_gene = |p: usize| if p == 0 { calls[path[0]].gene5 } else { calls[path[p - 1]].gene3 };
        let plus: Vec<bool> = (0..=k).map(|p| gi.genes[part_gene(p) as usize].3).collect();
        let bp_in: Vec<Option<i64>> = (0..=k).map(|p| (p > 0).then(|| calls[path[p - 1]].bp3)).collect();
        let bp_out: Vec<Option<i64>> = (0..=k).map(|p| (p < k).then(|| calls[path[p]].bp5)).collect();
        let pieces_of = |xs: &[usize]| -> Vec<Vec<(i64, i64)>> {
            std::iter::once(junctions[xs[0]].blocks5.clone()).chain(xs.iter().map(|&x| junctions[x].blocks3.clone())).collect()
        };
        let mut call_reads: Vec<CallRead> = spans.iter().map(|(r, js)| CallRead { read: *r, first: 0, pieces: pieces_of(js) }).collect();
        let position: FxHashMap<usize, usize> = path.iter().enumerate().map(|(q, &c)| (c, q)).collect();
        let mut partial: Vec<usize> = path.iter().flat_map(|&c| reads_of_call[c].iter().copied()).filter(|r| !reads.contains(r)).collect();
        partial.sort_unstable();
        partial.dedup();
        for &r in &partial {
            let mut at: Vec<(usize, usize)> = per_read.get(&r).map_or(Vec::new(), |v| v.iter()
                .filter_map(|e| position.get(&e.2).map(|&q| (q, e.3))).collect());
            at.sort_unstable();
            at.dedup_by_key(|e| e.0);
            // the longest run of consecutive junctions of the path
            let (mut best, mut s0) = ((0usize, 0usize), 0usize);
            for e in 1..=at.len() {
                if e == at.len() || at[e].0 != at[e - 1].0 + 1 {
                    if e - s0 > best.1 - best.0 { best = (s0, e); }
                    s0 = e;
                }
            }
            if best.1 > best.0 {
                let xs: Vec<usize> = at[best.0..best.1].iter().map(|e| e.1).collect();
                call_reads.push(CallRead { read: r, first: at[best.0].0, pieces: pieces_of(&xs) });
            }
        }
        let (isoforms, of) = call_isoforms(&call_reads, &plus, &bp_in, &bp_out);
        if let Some(top) = isoforms.first().filter(|t| t.first_part == 0 && t.blocks.len() == k + 1) {
            blocks = top.blocks.clone();
            chain_reads = vec![top.reads; k + 1];
        }
        let iso = |r: &usize| of.get(r).map(|&(x, fit)| (x as u32, fit));
        multi.push(MultiFusionCall {
            calls: path.clone(), reads: reads.len(), compatible_reads: compatible.len(),
            reads_by_dataset: by_ds.iter().map(|(&d, r)| (datasets[d].clone(), r.len())).collect(),
            filter, example_reads: reads.iter().take(5).map(|&r| read_keys[r].1.clone()).collect(),
            read_names: reads.iter().map(|&r| (read_keys[r].0.clone(), read_keys[r].1.clone())).collect(),
            blocks, chain_reads,
            read_isoform: reads.iter().map(iso).collect(),
            partial_reads: partial.iter().map(|r| ((read_keys[*r].0.clone(), read_keys[*r].1.clone()), iso(r))).collect(),
            isoforms,
        });
    }
    multi.sort_by(|a, b| b.reads.cmp(&a.reads).then(a.calls.cmp(&b.calls)));
    log::info!("[Fusion] {} junctions: {} candidates by breakpoints in {:.1} s, merged across short exons in {:.1} s, paralogs, hubs and transcripts in {:.1} s, {} multi-segment paths in {:.1} s",
               junctions.len(), calls.len(), s_cluster, s_merge, s_mid, multi.len(), t_paths.elapsed().as_secs_f64());
    // A passing path through the same chain of genes as one with more reads is listed as its isoform.
    let mut chains_seen: FxHashSet<Vec<(u32, u32)>> = FxHashSet::default();
    for m in multi.iter_mut().filter(|m| m.filter == "PASS") {
        let chain: Vec<(u32, u32)> = m.calls.iter().map(|&c| (calls[c].gene5, calls[c].gene3)).collect();
        if !chains_seen.insert(chain) { m.filter = "isoform"; }
    }
    // The junctions of a passing multi-segment fusion, or of an isoform of one, are reported under it.
    for m in multi.iter().filter(|m| m.filter == "PASS" || m.filter == "isoform") {
        for &c in &m.calls {
            if !artefact(calls[c].filter) || continues[c] { calls[c].filter = "in_multi_segment"; }
        }
    }
    let names = gi.genes.iter().map(|g| format!("{}\t{}\t{}\t{}\t{}\t{}", g.0, g.1, g.2, if g.3 { "+" } else { "-" }, g.4, g.5)).collect();
    // Each listed call's transcripts get their IDs in the order the calls are written: pairwise calls with `min_reads`
    // reads, then multi-segment calls.
    let mut n = 0usize;
    for c in calls.iter_mut().filter(|c| c.reads >= opts.min_reads) {
        for iso in c.isoforms.iter_mut() { n += 1; iso.id = fusion_isoform_id(n); }
    }
    for m in multi.iter_mut() {
        for iso in m.isoforms.iter_mut() { n += 1; iso.id = fusion_isoform_id(n); }
    }
    (calls, multi, names, datasets)
}

/// Writes the Parquet tables of the listed calls (pairwise with `min_reads_listed` reads, and multi-segment):
/// <prefix>_isoforge_fusions.parquet and the `fusion_parts`, `fusion_reads` and `fusion_isoforms` tables beside it.
/// Genomic coordinates are 1-based inclusive; query coordinates 0-based half-open on the read as sequenced.
pub fn write_fusion_parquet(prefix: &str, calls: &[FusionCall], multi: &[MultiFusionCall], genes: &[String], datasets: &[String],
                            segs: &[ChimericSegment], min_reads_listed: usize) -> IsoForgeResult<(String, String, String)> {
    // gene strings: id, name, chromosome, strand, start, end
    let gene = |g: u32| -> Vec<&str> { genes[g as usize].split('\t').collect() };
    let span = |g: u32| -> (i64, i64) { let x = gene(g); (x[4].parse().unwrap_or(0), x[5].parse().unwrap_or(0)) };
    let fid = |i: usize| format!("F{:06}", i + 1);
    struct Part { gene: u32, blocks: Vec<(i64, i64)>, jin: Option<(i64, bool)>, jout: Option<(i64, bool)>, chain_reads: usize }
    struct Row<'a> {
        id: String, kind: &'static str, filter: &'static str, category: &'static str, parts: Vec<Part>, reads: usize,
        /// Of `reads`, those the aligner split (a multi-segment call's: those with no part placed on annotated exons).
        split: usize,
        compatible: usize, constituents: Vec<String>, per_dataset: Vec<i32>,
        /// (read, spans every junction, its isoform and whether it fits it)
        read_names: Vec<(&'a (String, String), bool, Option<(u32, bool)>)>,
        isoforms: &'a [FusionIsoform],
    }
    let per_ds = |m: &FxHashMap<String, usize>| -> Vec<i32> { datasets.iter().map(|d| m.get(d).copied().unwrap_or(0) as i32).collect() };
    let mut by_read: FxHashMap<(&str, &str), Vec<&ChimericSegment>> = FxHashMap::default();
    for sg in segs { by_read.entry((sg.dataset.as_str(), sg.read_name.as_str())).or_default().push(sg); }
    for v in by_read.values_mut() { v.sort_by_key(|x| (x.query_start, x.query_end, x.ref_start)); }
    let mut rows: Vec<Row> = Vec::new();
    for (i, c) in calls.iter().enumerate().filter(|(_, c)| c.reads >= min_reads_listed) {
        rows.push(Row {
            id: fid(i), kind: "pairwise", filter: c.filter, category: c.category,
            parts: vec![
                Part { gene: c.gene5, blocks: c.blocks5.clone(), jin: None, jout: Some((c.bp5, c.boundary5)), chain_reads: c.chain_reads5 },
                Part { gene: c.gene3, blocks: c.blocks3.clone(), jin: Some((c.bp3, c.boundary3)), jout: None, chain_reads: c.chain_reads3 },
            ],
            reads: c.reads, split: c.reads_split, compatible: c.reads, constituents: Vec::new(), per_dataset: per_ds(&c.reads_by_dataset),
            read_names: c.read_names.iter().enumerate().map(|(k, r)| (r, true, c.read_isoform.get(k).copied().flatten())).collect(),
            isoforms: &c.isoforms,
        });
    }
    for (i, mc) in multi.iter().enumerate() {
        let k = mc.calls.len();
        let parts = (0..=k).map(|p| {
            let jin = (p > 0).then(|| { let c = &calls[mc.calls[p - 1]]; (c.bp3, c.boundary3) });
            let jout = (p < k).then(|| { let c = &calls[mc.calls[p]]; (c.bp5, c.boundary5) });
            let g = if p == 0 { calls[mc.calls[0]].gene5 } else { calls[mc.calls[p - 1]].gene3 };
            Part { gene: g, blocks: mc.blocks[p].clone(), jin, jout, chain_reads: mc.chain_reads[p] }
        }).collect();
        let mut read_names: Vec<(&(String, String), bool, Option<(u32, bool)>)> = mc.read_names.iter().enumerate()
            .map(|(k, r)| (r, true, mc.read_isoform.get(k).copied().flatten())).collect();
        read_names.extend(mc.partial_reads.iter().map(|(r, iso)| (r, false, *iso)));
        rows.push(Row {
            id: format!("M{:06}", i + 1), kind: "multi_segment", filter: mc.filter, category: "multi-segment", parts,
            reads: mc.reads, compatible: mc.compatible_reads, constituents: mc.calls.iter().map(|&c| fid(c)).collect(),
            split: mc.read_names.iter().filter(|(ds, name)| by_read.get(&(ds.as_str(), name.as_str())).map_or(true, |v| v.iter().all(|x| !x.placed))).count(),
            per_dataset: per_ds(&mc.reads_by_dataset), read_names, isoforms: &mc.isoforms,
        });
    }
    let list_str = |v: Vec<Vec<&str>>| -> ListArray {
        let mut b = ListBuilder::new(StringBuilder::new());
        for x in v { for y in x { b.values().append_value(y); } b.append(true); }
        b.finish()
    };
    let list_i64 = |v: Vec<Vec<i64>>| -> ListArray {
        ListArray::from_iter_primitive::<arrow::datatypes::Int64Type, _, _>(v.into_iter().map(|x| Some(x.into_iter().map(Some).collect::<Vec<_>>())))
    };
    let write = |path: &str, fields: Vec<(String, ArrayRef, bool)>| -> IsoForgeResult<()> {
        let schema = Arc::new(Schema::new(fields.iter().map(|(n, a, null)| Field::new(n.as_str(), a.data_type().clone(), *null)).collect::<Vec<_>>()));
        let batch = RecordBatch::try_new(schema.clone(), fields.into_iter().map(|(_, a, _)| a).collect())?;
        let mut w = ArrowWriter::try_new(File::create(path)?, schema, Some(crate::parquet_io::writer_props()))?;
        w.write(&batch)?;
        w.close()?;
        Ok(())
    };
    let col = |n: &str, a: ArrayRef| (n.to_string(), a, false);
    let nullable = |n: &str, a: ArrayRef| (n.to_string(), a, true);
    // calls
    let cpath = format!("{}_isoforge_fusions.parquet", prefix);
    let mut cf = vec![
        col("fusion_id", Arc::new(StringArray::from_iter_values(rows.iter().map(|r| r.id.as_str())))),
        col("kind", Arc::new(StringArray::from_iter_values(rows.iter().map(|r| r.kind)))),
        col("filter", Arc::new(StringArray::from_iter_values(rows.iter().map(|r| r.filter)))),
        col("category", Arc::new(StringArray::from_iter_values(rows.iter().map(|r| r.category)))),
        col("n_parts", Arc::new(Int32Array::from_iter_values(rows.iter().map(|r| r.parts.len() as i32)))),
        col("gene_ids", Arc::new(list_str(rows.iter().map(|r| r.parts.iter().map(|p| gene(p.gene)[0]).collect()).collect()))),
        col("gene_names", Arc::new(list_str(rows.iter().map(|r| r.parts.iter().map(|p| gene(p.gene)[1]).collect()).collect()))),
        col("reads", Arc::new(Int32Array::from_iter_values(rows.iter().map(|r| r.reads as i32)))),
        col("reads_split", Arc::new(Int32Array::from_iter_values(rows.iter().map(|r| r.split as i32)))),
        col("reads_compatible", Arc::new(Int32Array::from_iter_values(rows.iter().map(|r| r.compatible as i32)))),
        col("constituent_ids", Arc::new(list_str(rows.iter().map(|r| r.constituents.iter().map(|x| x.as_str()).collect()).collect()))),
        col("n_isoforms", Arc::new(Int32Array::from_iter_values(rows.iter().map(|r| r.isoforms.len() as i32)))),
        nullable("top_isoform", Arc::new(StringArray::from(rows.iter().map(|r| r.isoforms.first().map(|x| x.id.as_str())).collect::<Vec<_>>()))),
    ];
    for (d, ds) in datasets.iter().enumerate() {
        cf.push(col(&format!("reads_{}", ds), Arc::new(Int32Array::from_iter_values(rows.iter().map(|r| r.per_dataset[d])))));
    }
    write(&cpath, cf)?;
    // parts
    let ppath = format!("{}_isoforge_fusion_parts.parquet", prefix);
    let prow: Vec<(&str, usize, &Part)> = rows.iter().flat_map(|r| r.parts.iter().enumerate().map(move |(k, p)| (r.id.as_str(), k + 1, p))).collect();
    let pg: Vec<Vec<&str>> = prow.iter().map(|(_, _, p)| gene(p.gene)).collect();
    write(&ppath, vec![
        col("fusion_id", Arc::new(StringArray::from_iter_values(prow.iter().map(|x| x.0)))),
        col("part", Arc::new(Int32Array::from_iter_values(prow.iter().map(|x| x.1 as i32)))),
        col("gene_id", Arc::new(StringArray::from_iter_values(pg.iter().map(|g| g[0])))),
        col("gene_name", Arc::new(StringArray::from_iter_values(pg.iter().map(|g| g[1])))),
        col("chromosome", Arc::new(StringArray::from_iter_values(pg.iter().map(|g| g[2])))),
        col("strand", Arc::new(StringArray::from_iter_values(pg.iter().map(|g| g[3])))),
        col("gene_start", Arc::new(Int64Array::from_iter_values(prow.iter().map(|x| span(x.2.gene).0)))),
        col("gene_end", Arc::new(Int64Array::from_iter_values(prow.iter().map(|x| span(x.2.gene).1)))),
        col("block_starts", Arc::new(list_i64(prow.iter().map(|x| x.2.blocks.iter().map(|b| b.0).collect()).collect()))),
        col("block_ends", Arc::new(list_i64(prow.iter().map(|x| x.2.blocks.iter().map(|b| b.1).collect()).collect()))),
        nullable("junction_in", Arc::new(Int64Array::from(prow.iter().map(|x| x.2.jin.map(|j| j.0)).collect::<Vec<_>>()))),
        nullable("junction_out", Arc::new(Int64Array::from(prow.iter().map(|x| x.2.jout.map(|j| j.0)).collect::<Vec<_>>()))),
        nullable("boundary_in", Arc::new(BooleanArray::from(prow.iter().map(|x| x.2.jin.map(|j| j.1)).collect::<Vec<_>>()))),
        nullable("boundary_out", Arc::new(BooleanArray::from(prow.iter().map(|x| x.2.jout.map(|j| j.1)).collect::<Vec<_>>()))),
        col("chain_reads", Arc::new(Int32Array::from_iter_values(prow.iter().map(|x| x.2.chain_reads as i32)))),
    ])?;
    // reads: every aligned piece of every supporting read, on the part whose locus (gene span and blocks) it overlaps
    // most; a piece on none of them (a further partner, a paralog placement) has no part
    let (mut r_id, mut r_ds, mut r_name, mut r_span, mut r_seg, mut r_part) = (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new());
    let (mut r_iso, mut r_fit): (Vec<Option<i32>>, Vec<Option<bool>>) = (Vec::new(), Vec::new());
    let (mut r_chrom, mut r_rev, mut r_sup, mut r_mapq, mut r_placed) = (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new());
    let mut r_iso_id: Vec<Option<&str>> = Vec::new();
    let (mut r_rs, mut r_re, mut r_bs, mut r_be, mut r_qs, mut r_qe, mut r_len) = (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new());
    for row in &rows {
        let loci: Vec<(&str, i64, i64)> = row.parts.iter().map(|p| {
            let (s, e) = span(p.gene);
            let lo = p.blocks.first().map_or(s, |b| b.0.min(s));
            let hi = p.blocks.last().map_or(e, |b| b.1.max(e));
            (gene(p.gene)[2], lo, hi)
        }).collect();
        for &((ds, name), spans_all, iso) in &row.read_names {
            let Some(pieces) = by_read.get(&(ds.as_str(), name.as_str())) else { continue };
            for (k, sg) in pieces.iter().enumerate() {
                let part = loci.iter().enumerate()
                    .filter(|(_, (c, lo, hi))| *c == sg.chromosome && sg.ref_end >= *lo && sg.ref_start <= *hi)
                    .max_by_key(|(x, (_, lo, hi))| (sg.ref_end.min(*hi) - sg.ref_start.max(*lo), std::cmp::Reverse(*x)))
                    .map(|(x, _)| x as i32 + 1);
                let b = sg.blocks();
                r_id.push(row.id.as_str()); r_ds.push(ds.as_str()); r_name.push(name.as_str()); r_span.push(spans_all);
                r_iso.push(iso.map(|(k, _)| k as i32 + 1));
                r_iso_id.push(iso.and_then(|(k, _)| row.isoforms.get(k as usize)).map(|x| x.id.as_str()));
                r_fit.push(iso.map(|(_, f)| f));
                r_seg.push(k as i32 + 1); r_part.push(part);
                r_chrom.push(sg.chromosome.as_str()); r_rev.push(sg.reverse); r_sup.push(sg.supplementary); r_mapq.push(sg.mapq as i32); r_placed.push(sg.placed);
                r_rs.push(sg.ref_start); r_re.push(sg.ref_end);
                r_bs.push(b.iter().map(|x| x.0).collect::<Vec<_>>()); r_be.push(b.iter().map(|x| x.1).collect::<Vec<_>>());
                r_qs.push(sg.query_start); r_qe.push(sg.query_end); r_len.push(sg.read_length);
            }
        }
    }
    let rpath = format!("{}_isoforge_fusion_reads.parquet", prefix);
    write(&rpath, vec![
        col("fusion_id", Arc::new(StringArray::from(r_id))),
        col("dataset", Arc::new(StringArray::from(r_ds))),
        col("read_name", Arc::new(StringArray::from(r_name))),
        col("spans_all", Arc::new(BooleanArray::from(r_span))),
        nullable("isoform", Arc::new(Int32Array::from(r_iso))),
        nullable("isoform_id", Arc::new(StringArray::from(r_iso_id))),
        nullable("fits_isoform", Arc::new(BooleanArray::from(r_fit))),
        col("segment", Arc::new(Int32Array::from(r_seg))),
        nullable("part", Arc::new(Int32Array::from(r_part))),
        col("chromosome", Arc::new(StringArray::from(r_chrom))),
        col("reverse", Arc::new(BooleanArray::from(r_rev))),
        col("supplementary", Arc::new(BooleanArray::from(r_sup))),
        col("mapq", Arc::new(Int32Array::from(r_mapq))),
        col("ref_start", Arc::new(Int64Array::from(r_rs))),
        col("ref_end", Arc::new(Int64Array::from(r_re))),
        col("block_starts", Arc::new(list_i64(r_bs))),
        col("block_ends", Arc::new(list_i64(r_be))),
        col("query_start", Arc::new(Int32Array::from(r_qs))),
        col("query_end", Arc::new(Int32Array::from(r_qe))),
        col("read_length", Arc::new(Int32Array::from(r_len))),
        col("placed", Arc::new(BooleanArray::from(r_placed))),
    ])?;
    // isoforms: one row per call, isoform and part it covers, with the part's exon blocks and the breakpoints where
    // the isoform enters and leaves it (none at its own 5' and 3' ends, which come from its reads)
    let ipath = format!("{}_isoforge_fusion_isoforms.parquet", prefix);
    let (mut i_id, mut i_iso, mut i_full, mut i_first, mut i_last, mut i_reads, mut i_part) =
        (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new());
    let mut i_iso_id: Vec<&str> = Vec::new();
    let (mut i_g, mut i_bs, mut i_be, mut i_in, mut i_out) = (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new());
    for row in &rows {
        for (k, iso) in row.isoforms.iter().enumerate() {
            let full = iso.first_part == 0 && iso.last_part() == row.parts.len() - 1;
            for (q, b) in iso.blocks.iter().enumerate() {
                let p = iso.first_part + q;
                let part = &row.parts[p];
                i_id.push(row.id.as_str()); i_iso.push(k as i32 + 1); i_iso_id.push(iso.id.as_str()); i_full.push(full);
                i_first.push(iso.first_part as i32 + 1); i_last.push(iso.last_part() as i32 + 1); i_reads.push(iso.reads as i32);
                i_part.push(p as i32 + 1); i_g.push(gene(part.gene));
                i_bs.push(b.iter().map(|x| x.0).collect::<Vec<_>>()); i_be.push(b.iter().map(|x| x.1).collect::<Vec<_>>());
                i_in.push(if p > iso.first_part { part.jin.map(|j| j.0) } else { None });
                i_out.push(if p < iso.last_part() { part.jout.map(|j| j.0) } else { None });
            }
        }
    }
    write(&ipath, vec![
        col("fusion_id", Arc::new(StringArray::from(i_id))),
        col("isoform", Arc::new(Int32Array::from(i_iso))),
        col("isoform_id", Arc::new(StringArray::from(i_iso_id))),
        col("complete", Arc::new(BooleanArray::from(i_full))),
        col("first_part", Arc::new(Int32Array::from(i_first))),
        col("last_part", Arc::new(Int32Array::from(i_last))),
        col("reads", Arc::new(Int32Array::from(i_reads))),
        col("part", Arc::new(Int32Array::from(i_part))),
        col("gene_id", Arc::new(StringArray::from_iter_values(i_g.iter().map(|g| g[0])))),
        col("gene_name", Arc::new(StringArray::from_iter_values(i_g.iter().map(|g| g[1])))),
        col("chromosome", Arc::new(StringArray::from_iter_values(i_g.iter().map(|g| g[2])))),
        col("strand", Arc::new(StringArray::from_iter_values(i_g.iter().map(|g| g[3])))),
        col("block_starts", Arc::new(list_i64(i_bs))),
        col("block_ends", Arc::new(list_i64(i_be))),
        nullable("junction_in", Arc::new(Int64Array::from(i_in))),
        nullable("junction_out", Arc::new(Int64Array::from(i_out))),
    ])?;
    Ok((cpath, ppath, rpath))
}

/// Writes <prefix>_isoforge_fusion_models.tsv, one row per part of each listed fusion (pairwise calls with
/// `min_reads_listed` reads, and multi-segment calls) with its gene, blocks, breakpoints, reads and status, and
/// <prefix>_isoforge_fusion_reads.tsv, every supporting read of each call.
pub fn write_fusion_models(prefix: &str, calls: &[FusionCall], multi: &[MultiFusionCall], genes: &[String],
                           min_reads_listed: usize) -> IsoForgeResult<(String, String)> {
    let mpath = format!("{}_isoforge_fusion_models.tsv", prefix);
    let rpath = format!("{}_isoforge_fusion_reads.tsv", prefix);
    let mut m = BufWriter::new(File::create(&mpath)?);
    let mut r = BufWriter::new(File::create(&rpath)?);
    writeln!(m, "fusion_id\tpart\tgene_id\tgene_name\tchromosome\tstrand\tblocks\tjunction_breakpoints\treads\tchain_reads\tfilter")?;
    writeln!(r, "fusion_id\tdataset\tread_name")?;
    let fmt = |b: &[(i64, i64)]| b.iter().map(|(a, e)| format!("{}-{}", a, e)).collect::<Vec<_>>().join(",");
    let gene = |g: u32| -> Vec<String> { genes[g as usize].split('\t').map(|x| x.to_string()).collect() };
    for (i, c) in calls.iter().enumerate().filter(|(_, c)| c.reads >= min_reads_listed) {
        let id = format!("F{:06}", i + 1);
        for (part, g, blocks, bp, n) in [(1, c.gene5, &c.blocks5, c.bp5, c.chain_reads5), (2, c.gene3, &c.blocks3, c.bp3, c.chain_reads3)] {
            let x = gene(g);
            writeln!(m, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", id, part, x[0], x[1], x[2], x[3], fmt(blocks), bp, c.reads, n, c.filter)?;
        }
        for (d, n) in &c.read_names { writeln!(r, "{}\t{}\t{}", id, d, n)?; }
    }
    for (i, mc) in multi.iter().enumerate() {
        let id = format!("M{:06}", i + 1);
        let k = mc.calls.len();
        for part in 0..=k {
            // the breakpoints where the transcript enters and leaves the part
            let (g, junction) = if part == 0 {
                let c = &calls[mc.calls[0]]; (c.gene5, c.bp5.to_string())
            } else if part == k {
                let c = &calls[mc.calls[k - 1]]; (c.gene3, c.bp3.to_string())
            } else {
                let (a, b) = (&calls[mc.calls[part - 1]], &calls[mc.calls[part]]);
                (a.gene3, format!("{};{}", a.bp3, b.bp5))
            };
            let x = gene(g);
            writeln!(m, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", id, part + 1, x[0], x[1], x[2], x[3], fmt(&mc.blocks[part]), junction,
                     mc.reads, mc.chain_reads[part], mc.filter)?;
        }
        for (d, n) in &mc.read_names { writeln!(r, "{}\t{}\t{}", id, d, n)?; }
    }
    Ok((mpath, rpath))
}

/// Writes <prefix>_isoforge_fusion_isoforms.gtf: every transcript of each listed call, one transcript record per part
/// it covers, with its exons. gene_id is the call (F... or M...); transcript_id the transcript's ID and its part,
/// fusion_transcript_id the transcript's ID alone; exons are numbered along the whole fusion transcript, 5' to 3'.
pub fn write_fusion_gtf(prefix: &str, calls: &[FusionCall], multi: &[MultiFusionCall], genes: &[String],
                        min_reads_listed: usize) -> IsoForgeResult<String> {
    let path = format!("{}_isoforge_fusion_isoforms.gtf", prefix);
    let mut out = BufWriter::new(File::create(&path)?);
    writeln!(out, "##description: IsoForge fusion transcripts; one transcript record per part, exons numbered along the fusion")?;
    let gene = |g: u32| -> Vec<&str> { genes[g as usize].split('\t').collect() };
    let mut write_call = |id: &str, part_genes: &[u32], status: &str, isoforms: &[FusionIsoform]| -> IsoForgeResult<()> {
        let names: Vec<&str> = part_genes.iter().map(|&g| gene(g)[1]).collect();
        let fusion_genes = names.join("::");
        for iso in isoforms {
            let complete = iso.first_part == 0 && iso.last_part() + 1 == part_genes.len();
            let mut exon_number = 0usize;
            for (q, blocks) in iso.blocks.iter().enumerate() {
                let p = iso.first_part + q;
                let g = gene(part_genes[p]);
                let (chrom, strand) = (g[2], g[3]);
                let (Some(first), Some(last)) = (blocks.first(), blocks.last()) else { continue };
                let attrs = format!("gene_id \"{id}\"; transcript_id \"{}.{}\"; fusion_transcript_id \"{}\"; fusion_genes \"{fusion_genes}\"; \
                                     part \"{}\"; parts \"{}\"; partner_gene_id \"{}\"; partner_gene_name \"{}\"; reads \"{}\"; \
                                     complete \"{complete}\"; status \"{status}\";",
                                    iso.id, p + 1, iso.id, p + 1, part_genes.len(), g[0], g[1], iso.reads);
                writeln!(out, "{chrom}\tIsoForge\ttranscript\t{}\t{}\t.\t{strand}\t.\t{attrs}", first.0, last.1)?;
                let order: Vec<&(i64, i64)> = if strand == "-" { blocks.iter().rev().collect() } else { blocks.iter().collect() };
                for b in order {
                    exon_number += 1;
                    writeln!(out, "{chrom}\tIsoForge\texon\t{}\t{}\t.\t{strand}\t.\t{attrs} exon_number \"{exon_number}\";", b.0, b.1)?;
                }
            }
        }
        Ok(())
    };
    for (i, c) in calls.iter().enumerate().filter(|(_, c)| c.reads >= min_reads_listed) {
        write_call(&format!("F{:06}", i + 1), &[c.gene5, c.gene3], c.filter, &c.isoforms)?;
    }
    for (i, m) in multi.iter().enumerate() {
        let mut g: Vec<u32> = vec![calls[m.calls[0]].gene5];
        g.extend(m.calls.iter().map(|&c| calls[c].gene3));
        write_call(&format!("M{:06}", i + 1), &g, m.filter, &m.isoforms)?;
    }
    Ok(path)
}

/// Writes <prefix>_isoforge_fusions.tsv (all candidates with a filter status) and <prefix>_isoforge_fusions.bedpe (PASS).
pub fn write_fusions(prefix: &str, calls: &[FusionCall], genes: &[String], datasets: &[String], min_reads_listed: usize)
                     -> IsoForgeResult<(String, usize)> {
    let tsv = format!("{}_isoforge_fusions.tsv", prefix);
    let mut out = BufWriter::new(File::create(&tsv)?);
    write!(out, "fusion_id\tgene_id_5p\tgene_name_5p\tchrom_5p\tstrand_5p\tbreakpoint_5p\texon_boundary_5p\t\
                 gene_id_3p\tgene_name_3p\tchrom_3p\tstrand_3p\tbreakpoint_3p\texon_boundary_3p\tcategory\treads")?;
    for d in datasets { write!(out, "\treads_{}", d)?; }
    writeln!(out, "\tfilter\texample_reads\tisoforms\ttop_isoform\treads_split")?;
    let bedpe = format!("{}_isoforge_fusions.bedpe", prefix);
    let mut bp = BufWriter::new(File::create(&bedpe)?);
    let mut n_pass = 0;
    for (i, c) in calls.iter().enumerate().filter(|(_, c)| c.reads >= min_reads_listed) {
        let g5: Vec<&str> = genes[c.gene5 as usize].split('\t').collect();
        let g3: Vec<&str> = genes[c.gene3 as usize].split('\t').collect();
        write!(out, "F{:06}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
               i + 1, g5[0], g5[1], g5[2], g5[3], c.bp5, c.boundary5, g3[0], g3[1], g3[2], g3[3], c.bp3, c.boundary3, c.category, c.reads)?;
        for d in datasets { write!(out, "\t{}", c.reads_by_dataset.get(d).copied().unwrap_or(0))?; }
        writeln!(out, "\t{}\t{}\t{}\t{}\t{}", c.filter, c.example_reads.join(","), c.isoforms.len(),
                 c.isoforms.first().map(|x| x.id.as_str()).unwrap_or(""), c.reads_split)?;
        if c.filter == "PASS" || c.filter == "in_multi_segment" {
            if c.filter == "PASS" { n_pass += 1; }
            writeln!(bp, "{}\t{}\t{}\t{}\t{}\t{}\tF{:06}_{}--{}\t{}\t{}\t{}", g5[2], c.bp5 - 1, c.bp5, g3[2], c.bp3 - 1, c.bp3,
                     i + 1, g5[1], g3[1], c.reads, g5[3], g3[3])?;
        }
    }
    Ok((tsv, n_pass))
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::{Cigar::*, CigarString};

    #[test]
    fn fusion_transcripts_are_named_like_the_run_s_models() {
        assert_eq!(fusion_isoform_id(1), "ISOFORF000000001");
        assert_eq!(fusion_isoform_id(123456), "ISOFORF000123456");
    }

    fn record(ops: Vec<Cigar>, reverse: bool, supplementary: bool) -> bam::Record {
        let len: u32 = ops.iter().map(|op| match op { Match(n) | Ins(n) | SoftClip(n) => *n, _ => 0 }).sum();
        let mut r = bam::Record::new();
        r.set(b"read1", Some(&CigarString(ops)), &vec![b'A'; len as usize], &vec![30; len as usize]);
        r.set_tid(0);
        r.set_pos(999);
        r.set_mapq(60);
        let mut flags = 0u16;
        if reverse { flags |= 0x10; }
        if supplementary { flags |= 0x800; }
        r.set_flags(flags);
        r.push_aux(b"SA", bam::record::Aux::String("chr2,100,+,50M,60,0;")).unwrap();
        r
    }

    #[test]
    fn read_through_does_not_skip_a_coding_gene() {
        let spans = [(100, 200), (1_000, 5_000), (8_000, 9_000)];
        // partners at 500 and 7,000: the gene at 1,000-5,000 lies between them
        assert!(span_between(&spans, 500, 7_000));
        // a breakpoint inside it: the gene is a partner, not between
        assert!(!span_between(&spans, 4_000, 7_000));
        assert!(!span_between(&spans, 500, 4_000));
        assert!(!span_between(&spans, 5_500, 7_900));
        assert!(!span_between(&[], 0, 10));
    }

    #[test]
    fn a_folded_molecule_retraces_itself_and_an_inversion_does_not() {
        // the second part copies the first on the other strand: their spans overlap
        assert!(parts_apart(&[(100, 200), (500, 700)], &[(520, 700)]) < 0);
        // it starts a few bases on, within the fold
        assert_eq!(parts_apart(&[(100, 200), (500, 700)], &[(760, 900)]), 60);
        // an inverted neighbour 20 kb away, on either side
        assert_eq!(parts_apart(&[(100, 200), (500, 700)], &[(20_700, 21_000), (22_000, 22_100)]), 20_000);
        assert_eq!(parts_apart(&[(20_700, 21_000)], &[(100, 200), (500, 700)]), 20_000);
        assert_eq!(parts_apart(&[], &[(1, 2)]), 0);
    }

    #[test]
    fn segments_place_their_part_of_the_read_in_sequenced_orientation() {
        let names = vec!["chr1".to_string()];
        // forward primary: 20 clipped, 100 aligned over an intron, 30 clipped
        let s = ChimericSegment::from_record(&record(vec![SoftClip(20), Match(60), RefSkip(500), Match(40), SoftClip(30)], false, false), &names, "d").unwrap();
        assert_eq!((s.query_start, s.query_end, s.read_length), (20, 120, 150));
        assert_eq!((s.ref_start, s.ref_end), (1000, 1599));
        assert_eq!(s.introns, vec![1060, 1559]);
        assert_eq!(s.blocks(), vec![(1000, 1059), (1560, 1599)]);
        // reverse, hard-clipped supplementary: the CIGAR runs along the reverse complement, so its trailing
        // clip is the start of the read as sequenced
        let s = ChimericSegment::from_record(&record(vec![HardClip(70), Match(50), HardClip(30)], true, true), &names, "d").unwrap();
        assert_eq!((s.query_start, s.query_end, s.read_length), (30, 80, 150));
        assert!(s.supplementary && s.reverse);
        // no SA tag: not a chimeric read
        let mut plain = record(vec![Match(50)], false, false);
        plain.remove_aux(b"SA").unwrap();
        assert!(ChimericSegment::from_record(&plain, &names, "d").is_none());
    }

    #[test]
    fn spliced_reads_are_cut_only_where_an_intron_joins_two_genes() {
        use crate::data_structures::{GeneIntervalTree, GeneRecord};
        use crate::types::Strand;
        let mut tree = GeneIntervalTree::new();
        // gene A 900-1600, gene B 3000-5000, gene C 900-5000 on chr2 only
        tree.insert(GeneRecord { gene_id: 1, chromosome: "chr1".into(), start: 900, end: 1600, strand: Strand::Plus });
        tree.insert(GeneRecord { gene_id: 2, chromosome: "chr1".into(), start: 3000, end: 5000, strand: Strand::Plus });
        tree.insert(GeneRecord { gene_id: 3, chromosome: "chr2".into(), start: 900, end: 5000, strand: Strand::Plus });
        tree.finalize();
        let spans = GeneSpans::from_tree(&tree);
        let names = vec!["chr1".to_string()];
        let mut plain = record(vec![SoftClip(10), Match(60), RefSkip(500), Match(40), RefSkip(1900), Match(80)], false, false);
        plain.remove_aux(b"SA").unwrap();
        // 1000-1059 | intron 1060-1559 (inside A) | 1560-1599 | intron 1600-3499 (A to B) | 3500-3579
        let parts = ChimericSegment::split_at_gene_joins(&plain, &names, "d", &spans);
        assert_eq!(parts.len(), 2);
        assert_eq!((parts[0].ref_start, parts[0].ref_end, parts[0].query_start, parts[0].query_end), (1000, 1599, 10, 110));
        assert_eq!(parts[0].introns, vec![1060, 1559]);
        assert_eq!((parts[1].ref_start, parts[1].ref_end, parts[1].query_start, parts[1].query_end), (3500, 3579, 110, 190));
        // one gene holding both ends: no cut, and without an SA tag no segment at all
        let names2 = vec!["chr2".to_string()];
        assert!(ChimericSegment::split_at_gene_joins(&plain, &names2, "d", &spans).is_empty());
        // an SA record is kept whole when nothing cuts it
        let sa = record(vec![Match(50)], false, false);
        assert_eq!(ChimericSegment::split_at_gene_joins(&sa, &names, "d", &spans).len(), 1);
    }

    #[test]
    fn fusion_isoforms_follow_the_reads_ends_and_introns() {
        let r = |read: usize, five: Vec<(i64, i64)>, three: Vec<(i64, i64)>| CallRead { read, first: 0, pieces: vec![five, three] };
        let mut reads = Vec::new();
        // full-length reads, three ending at each of two poly(A) sites; one with a 19-bp false intron
        for i in 0..3 { reads.push(r(i, vec![(200, 300), (400, 1000)], vec![(5000, 5100), (5200, 5600)])); }
        for i in 3..6 { reads.push(r(i, vec![(200 + (i as i64 - 3) * 2, 300), (400, 1000)], vec![(5000, 5100), (5200, 5900)])); }
        reads.push(r(6, vec![(200, 300), (400, 700), (720, 1000)], vec![(5000, 5100), (5200, 5602)]));
        // 5'-truncated reads starting in the second exon
        for i in 7..11 { reads.push(r(i, vec![(450 + i as i64, 1000)], vec![(5000, 5100), (5200, 5600 - (i as i64 % 2) * 2)])); }
        // a read reaching back into the first intron fits no chain as a truncated copy
        reads.push(r(11, vec![(350, 1000)], vec![(5000, 5100), (5200, 5600)]));
        let (iso, of) = call_isoforms(&reads, &[true, true], &[None, Some(5000)], &[Some(1000), None]);
        assert_eq!(iso.len(), 2, "{iso:?}");
        assert_eq!(iso[0].blocks, vec![vec![(200, 300), (400, 1000)], vec![(5000, 5100), (5200, 5600)]]);
        assert_eq!(iso[0].reads, 9);
        assert_eq!(iso[1].blocks[1], vec![(5000, 5100), (5200, 5900)]);
        assert_eq!(iso[1].reads, 3);
        assert_eq!(of.len(), 12);
        assert!((7..11).all(|i| of[&i] == (0, true)) && of[&6] == (0, true) && (3..6).all(|i| of[&i] == (1, true)));
        // the read reaching back into the first intron is placed, not fitted
        assert_eq!(of[&11], (0, false));
    }

    #[test]
    fn fusion_isoforms_read_a_minus_strand_partner_3_prime_to_5_prime() {
        // 3' partner on -: entered at 9000, exons [8500, 9000] then [8000, 8400] in transcript order
        let r = |read: usize, three: Vec<(i64, i64)>| CallRead { read, first: 0, pieces: vec![vec![(100, 500)], three] };
        let mut reads: Vec<CallRead> = (0..4).map(|i| r(i, vec![(8000, 8400), (8500, 9000)])).collect();
        reads.push(r(4, vec![(8600, 9000)]));       // ends in the first exon: a 3'-truncated copy
        let (iso, of) = call_isoforms(&reads, &[true, false], &[None, Some(9000)], &[Some(500), None]);
        assert_eq!(iso.len(), 1, "{iso:?}");
        assert_eq!(iso[0].blocks[1], vec![(8000, 8400), (8500, 9000)]);
        assert_eq!(iso[0].reads, 5);
        assert_eq!(of[&4], (0, false));
        // a read ending in the first exon of the 3' part ends at another poly(A) site: not a copy of the chain
        let chain = vec![vec![], vec![(8400, 8500)]];
        let (five, a) = (orient(&[(100, 500)], true), orient(&[(8600, 9000)], false));
        assert!(!fits(&[five, a], &chain, &[true, false], 0));
        // 5' partner on -: exons [700, 900] then [300, 500]; a read starting in the second exon is a truncated copy,
        // one starting 100 bp into the intron is not
        let chain = vec![vec![(500, 700)], vec![(8400, 8500)]];
        let three = || orient(&[(8000, 8400), (8500, 9000)], false);
        assert!(fits(&[orient(&[(300, 480)], false), three()], &chain, &[false, false], 0));
        assert!(!fits(&[orient(&[(300, 600)], false), three()], &chain, &[false, false], 0));
        assert!(fits(&[orient(&[(300, 500), (700, 880)], false), three()], &chain, &[false, false], 0));
    }

    #[test]
    fn a_middle_part_enters_and_leaves_at_its_breakpoints_on_either_strand() {
        let (a, b, c) = (vec![(100, 200), (300, 400)], vec![(100, 200), (300, 400)], vec![(90, 200), (300, 410)]);
        let sides = vec![&a, &b, &c];
        // plus strand: enters at the low end (95), leaves at the high end (405)
        assert_eq!(part_blocks(&sides, Some(95), Some(405), true), (vec![(95, 200), (300, 405)], 3));
        // minus strand: enters at the high end, leaves at the low end
        assert_eq!(part_blocks(&sides, Some(405), Some(95), false), (vec![(95, 200), (300, 405)], 3));
        // a first part (5' partner on +) ends at its breakpoint, a last part (3' partner on +) starts at it
        assert_eq!(part_blocks(&sides, None, Some(405), true).0, vec![(90, 200), (300, 405)]);
        assert_eq!(part_blocks(&sides, Some(95), None, true).0, vec![(95, 200), (300, 410)]);
    }

    #[test]
    fn fusion_tables_place_each_read_piece_on_its_part() {
        use arrow::array::{BooleanArray, Int32Array, StringArray};
        let genes = vec!["GA\tA\tchr1\t+\t1000\t2000".to_string(), "GB\tB\tchr1\t-\t5000\t6000".to_string(),
                         "GC\tC\tchr2\t+\t100\t900".to_string()];
        let seg = |name: &str, chrom: &str, rs: i64, re: i64, qs: i32| ChimericSegment {
            read_name: name.into(), dataset: "d1".into(), chromosome: chrom.into(), reverse: false, supplementary: qs > 0,
            mapq: 60, ref_start: rs, ref_end: re, query_start: qs, query_end: qs + 100, read_length: 300, introns: Vec::new(), placed: false, short_class: false, pieces: 0, excess: 0 };
        // r1 spans A -> B -> C; r2 carries only B -> C (5'-truncated) and a piece elsewhere
        let segs = vec![seg("r1", "chr1", 1500, 1600, 0), seg("r1", "chr1", 5500, 5600, 100), seg("r1", "chr2", 200, 300, 200),
                        seg("r2", "chr1", 5500, 5600, 0), seg("r2", "chr2", 200, 300, 100), seg("r2", "chr9", 1, 50, 200)];
        let call = |g5: u32, g3: u32, bp5: i64, bp3: i64, reads: &[&str]| FusionCall {
            gene5: g5, gene3: g3, bp5, bp3, reads: reads.len(), reads_split: reads.len(), reads_by_dataset: [("d1".to_string(), reads.len())].into_iter().collect(),
            boundary5: true, boundary3: false, category: "intra-chromosomal", filter: "in_multi_segment", example_reads: Vec::new(),
            read_names: reads.iter().map(|r| ("d1".to_string(), r.to_string())).collect(),
            blocks5: Vec::new(), blocks3: Vec::new(), chain_reads5: 0, chain_reads3: 0, isoforms: Vec::new(), read_isoform: Vec::new(), unnamed: false };
        let calls = vec![call(1, 2, 5500, 200, &["r1", "r2"]), call(0, 1, 1600, 5600, &["r1"])];
        let multi = vec![MultiFusionCall {
            calls: vec![1, 0], reads: 1, compatible_reads: 2, reads_by_dataset: [("d1".to_string(), 1)].into_iter().collect(),
            filter: "PASS", example_reads: vec!["r1".into()], read_names: vec![("d1".into(), "r1".into())],
            blocks: vec![vec![(1500, 1600)], vec![(5500, 5600)], vec![(200, 300)]], chain_reads: vec![1, 1, 1],
            isoforms: vec![FusionIsoform { first_part: 0, blocks: vec![vec![(1500, 1600)], vec![(5500, 5600)], vec![(200, 300)]], reads: 1,
                                           id: fusion_isoform_id(1) },
                           FusionIsoform { first_part: 1, blocks: vec![vec![(5500, 5600)], vec![(200, 300)]], reads: 1,
                                           id: fusion_isoform_id(2) }],
            read_isoform: vec![Some((0, true))], partial_reads: vec![(("d1".into(), "r2".into()), Some((1, true)))] }];
        let dir = std::env::temp_dir().join(format!("isoforge_fusion_tables_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let prefix = dir.join("t").to_string_lossy().to_string();
        let (cpath, ppath, rpath) = write_fusion_parquet(&prefix, &calls, &multi, &genes, &["d1".to_string()], &segs, 1).unwrap();
        let read = |p: &str| {
            let r = ParquetRecordBatchReaderBuilder::try_new(File::open(p).unwrap()).unwrap().build().unwrap();
            arrow::compute::concat_batches(&r.schema(), &r.map(|b| b.unwrap()).collect::<Vec<_>>()).unwrap()
        };
        let (c, p, r) = (read(&cpath), read(&ppath), read(&rpath));
        let s = |b: &RecordBatch, n: &str| b.column_by_name(n).unwrap().as_any().downcast_ref::<StringArray>().unwrap().iter().map(|x| x.unwrap().to_string()).collect::<Vec<_>>();
        let i = |b: &RecordBatch, n: &str| b.column_by_name(n).unwrap().as_any().downcast_ref::<Int32Array>().unwrap().iter().collect::<Vec<_>>();
        assert_eq!(s(&c, "fusion_id"), vec!["F000001", "F000002", "M000001"]);
        assert_eq!(i(&c, "n_parts"), vec![Some(2), Some(2), Some(3)]);
        assert_eq!(i(&c, "reads_d1"), vec![Some(2), Some(1), Some(1)]);
        // the multi-segment call's parts in transcript order, with the gene spans
        let m: Vec<usize> = (0..p.num_rows()).filter(|&k| s(&p, "fusion_id")[k] == "M000001").collect();
        assert_eq!(m.iter().map(|&k| s(&p, "gene_name")[k].clone()).collect::<Vec<_>>(), vec!["A", "B", "C"]);
        let gs = p.column_by_name("gene_start").unwrap().as_any().downcast_ref::<Int64Array>().unwrap();
        assert_eq!(m.iter().map(|&k| gs.value(k)).collect::<Vec<_>>(), vec![1000, 5000, 100]);
        let jin = p.column_by_name("junction_in").unwrap().as_any().downcast_ref::<Int64Array>().unwrap();
        assert!(jin.is_null(m[0]));
        assert_eq!((jin.value(m[1]), jin.value(m[2])), (5600, 200));
        // every piece of each read, on the part it overlaps; r2 carries only part of the multi-segment path
        let mr: Vec<usize> = (0..r.num_rows()).filter(|&k| s(&r, "fusion_id")[k] == "M000001").collect();
        let names = s(&r, "read_name");
        let parts = i(&r, "part");
        let spans = r.column_by_name("spans_all").unwrap().as_any().downcast_ref::<BooleanArray>().unwrap();
        assert_eq!(mr.iter().map(|&k| (names[k].as_str(), parts[k], spans.value(k))).collect::<Vec<_>>(),
                   vec![("r1", Some(1), true), ("r1", Some(2), true), ("r1", Some(3), true),
                        ("r2", Some(2), false), ("r2", Some(3), false), ("r2", None, false)]);
        // a pairwise call's reads keep their other pieces, without a part
        let f1: Vec<Option<i32>> = (0..r.num_rows()).filter(|&k| s(&r, "fusion_id")[k] == "F000001" && names[k] == "r1").map(|k| parts[k]).collect();
        assert_eq!(f1, vec![None, Some(1), Some(2)]);
        // each read carries its isoform; the partial one has an isoform over parts 2-3 of its own
        let iso = i(&r, "isoform");
        assert_eq!(mr.iter().map(|&k| (names[k].as_str(), iso[k])).collect::<Vec<_>>(),
                   vec![("r1", Some(1)), ("r1", Some(1)), ("r1", Some(1)), ("r2", Some(2)), ("r2", Some(2)), ("r2", Some(2))]);
        let t = read(&format!("{}_isoforge_fusion_isoforms.parquet", prefix));
        let rows: Vec<usize> = (0..t.num_rows()).filter(|&k| s(&t, "fusion_id")[k] == "M000001").collect();
        assert_eq!(rows.iter().map(|&k| (i(&t, "isoform")[k], i(&t, "part")[k], i(&t, "first_part")[k])).collect::<Vec<_>>(),
                   vec![(Some(1), Some(1), Some(1)), (Some(1), Some(2), Some(1)), (Some(1), Some(3), Some(1)),
                        (Some(2), Some(2), Some(2)), (Some(2), Some(3), Some(2))]);
        let tin = t.column_by_name("junction_in").unwrap().as_any().downcast_ref::<Int64Array>().unwrap();
        // isoform 2 enters part 2 at its reads' 5' end, not at the path's junction
        assert!(tin.is_null(rows[3]) && tin.value(rows[4]) == 200);
        std::fs::remove_dir_all(&dir).ok();
    }
}

/// Writes <prefix>_isoforge_fusions_multi.tsv: fusions of three or more segments, the genes and breakpoints in
/// transcript order, the reads spanning every junction and those compatible with part of the path.
pub fn write_multi_fusions(prefix: &str, multi: &[MultiFusionCall], calls: &[FusionCall], genes: &[String], datasets: &[String])
                           -> IsoForgeResult<(String, usize)> {
    let tsv = format!("{}_isoforge_fusions_multi.tsv", prefix);
    let mut out = BufWriter::new(File::create(&tsv)?);
    write!(out, "fusion_id\tn_genes\tgene_ids\tgene_names\tchromosomes\tstrands\tbreakpoints\texon_boundaries\treads_spanning\treads_compatible")?;
    for d in datasets { write!(out, "\treads_{}", d)?; }
    writeln!(out, "\tfilter\texample_reads\tisoforms\ttop_isoform")?;
    let mut n_pass = 0;
    for (i, m) in multi.iter().enumerate() {
        let mut gl: Vec<u32> = vec![calls[m.calls[0]].gene5];
        gl.extend(m.calls.iter().map(|&c| calls[c].gene3));
        let f = |k: usize| -> Vec<String> { gl.iter().map(|&g| genes[g as usize].split('\t').nth(k).unwrap_or("").to_string()).collect() };
        let (ids, names, chroms, strands) = (f(0), f(1), f(2), f(3));
        let bps: Vec<String> = m.calls.iter().map(|&c| {
            let x = &calls[c];
            let c5 = genes[x.gene5 as usize].split('\t').nth(2).unwrap_or("");
            let c3 = genes[x.gene3 as usize].split('\t').nth(2).unwrap_or("");
            format!("{}:{}>{}:{}", c5, x.bp5, c3, x.bp3)
        }).collect();
        let eb: Vec<&str> = m.calls.iter().map(|&c| {
            let x = &calls[c];
            match (x.boundary5, x.boundary3) { (true, true) => "both", (true, false) => "5p", (false, true) => "3p", _ => "neither" }
        }).collect();
        write!(out, "M{:06}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", i + 1, gl.len(), ids.join("::"), names.join("::"),
               chroms.join("::"), strands.join("::"), bps.join(";"), eb.join(";"), m.reads, m.compatible_reads)?;
        for d in datasets { write!(out, "\t{}", m.reads_by_dataset.get(d).copied().unwrap_or(0))?; }
        writeln!(out, "\t{}\t{}\t{}\t{}", m.filter, m.example_reads.join(","), m.isoforms.len(),
                 m.isoforms.first().map(|x| x.id.as_str()).unwrap_or(""))?;
        if m.filter == "PASS" { n_pass += 1; }
    }
    Ok((tsv, n_pass))
}
