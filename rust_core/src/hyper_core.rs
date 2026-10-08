// The read hypergraph.
//
// Each spliced read's intron chain (the `introns` column of observed.parquet) is a hyperedge over splice sites and
// each candidate model, annotated or novel, a hyperpath. A read fits a model when its introns are consecutive introns
// of the model, exactly or each within `tol` bp at both sites, and its 5' and 3' ends lie in the model's exons before
// and after them (`slack` bp beyond an inner exon's boundary, `end_slack` beyond the model's first and last exon); an
// unspliced read fits the models holding it inside one exon. A fit records the read's 5' and 3' ends as exonic
// distances from the model's ends.

use std::fs::File;
use std::path::Path;
use std::sync::Arc;

use arrow::array::{Array, Int32Builder, Int64Builder, StringArray, UInt32Builder};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use coitrees::{COITree, Interval, IntervalTree};
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::arrow::{ArrowWriter, ProjectionMask};
use rayon::prelude::*;
use rustc_hash::FxHashMap;

use crate::errors::IsoForgeResult;
use crate::post_processing::PostDB;

/// Each read's own intron chain, interned.
pub struct ReadChains {
    /// Chain of each read of `db.observed`, in row order (0: unspliced).
    pub chain_of: Vec<u32>,
    /// Distinct chains: introns as (donor, acceptor) in 5'->3' order.
    pub chains: Vec<Box<[(i64, i64)]>>,
}

/// Reads the `introns` column of observed.parquet; None when the column is absent.
pub fn load_read_chains(observed: &Path) -> IsoForgeResult<Option<ReadChains>> {
    if !observed.exists() { return Ok(None); }
    let builder = ParquetRecordBatchReaderBuilder::try_new(File::open(observed)?)?;
    let Ok(col) = builder.schema().index_of("introns") else { return Ok(None) };
    let mask = ProjectionMask::roots(builder.parquet_schema(), [col]);
    let reader = builder.with_projection(mask).build()?;
    let mut index: FxHashMap<Box<[(i64, i64)]>, u32> = FxHashMap::default();
    let mut chains: Vec<Box<[(i64, i64)]>> = vec![Box::new([])];
    index.insert(Box::new([]), 0);
    let mut chain_of = Vec::new();
    let mut v: Vec<i64> = Vec::new();
    for batch in reader {
        let batch = batch?;
        let a = batch.column(0).as_any().downcast_ref::<StringArray>().unwrap();
        for i in 0..batch.num_rows() {
            if a.is_null(i) { chain_of.push(0); continue; }
            v.clear();
            v.extend(a.value(i).split(',').filter_map(|x| x.parse::<i64>().ok()));
            let c: Box<[(i64, i64)]> = v.chunks_exact(2).map(|p| (p[0], p[1])).collect();
            let n = chains.len() as u32;
            let id = match index.get(&c) {
                Some(&id) => id,
                None => { index.insert(c.clone(), n); chains.push(c); n }
            };
            chain_of.push(id);
        }
    }
    Ok(Some(ReadChains { chain_of, chains }))
}

/// A candidate model in transcript order.
pub struct Candidate {
    pub tid: i64,
    pub gene: i64,
    pub chrom: Arc<str>,
    pub plus: bool,
    /// Exons as (5' boundary, 3' boundary), 5' to 3'.
    pub exons: Vec<(i64, i64)>,
    /// Introns as (donor, acceptor), 5' to 3'.
    pub introns: Vec<(i64, i64)>,
    /// Exonic distance from the model's 5' end to each exon's 5' boundary.
    pub before: Vec<i64>,
    pub len: i64,
    /// Exonic offset of the site its reads place at each end, inward from the chain's end (negative: beyond it).
    pub site5: Option<i32>,
    pub site3: Option<i32>,
    /// A read beyond the chain's end and every site is not at the model's end.
    pub bounded: bool,
    /// Reads the model is counted for.
    pub reads: u32,
}

/// Where a read's end lies against a model's.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EndFit { At, Short, Beyond }

fn end_fit(off: i32, site: Option<i32>, bounded: bool) -> EndFit {
    if !bounded { return if off <= END_WINDOW { EndFit::At } else { EndFit::Short }; }
    if off.abs() <= END_WINDOW || site.is_some_and(|s| (off - s).abs() <= END_WINDOW) { return EndFit::At; }
    if off < site.map_or(0, |s| s.min(0)) { EndFit::Beyond } else { EndFit::Short }
}

impl Candidate {
    pub fn fit5(&self, f: &Fit) -> EndFit { end_fit(f.off5, self.site5, self.bounded) }
    pub fn fit3(&self, f: &Fit) -> EndFit { end_fit(f.off3, self.site3, self.bounded) }
    /// A read of the whole model: an exact fit with both ends at the model's.
    pub fn whole(&self, f: &Fit) -> bool { f.shift == 0 && self.fit5(f) == EndFit::At && self.fit3(f) == EndFit::At }
    /// Exonic length between the sites its reads place, where they lie inside the chain's ends.
    pub fn span(&self) -> i64 {
        (self.len - self.site5.map_or(0, |s| s.max(0)) as i64 - self.site3.map_or(0, |s| s.max(0)) as i64).max(1)
    }
}

/// Bases within which a read's end is at a model's end.
pub const END_WINDOW: i32 = 50;
/// What a read costs a model that does not hold one of its ends: a tailed read ending away from the model's 3' end,
/// a read starting upstream of the model's start.
const OFF_END: f64 = 0.02;
/// What a fit costs whose splice sites are near the model's, not on them.
const SHIFTED: f64 = 0.02;

/// The genomic TSS and poly(A) site the reads of each model place (post_processing's `build_tss_tes_info`).
pub type Sites = FxHashMap<i64, (Option<i64>, Option<i64>)>;

/// One read fitting one candidate.
#[derive(Clone, Copy, Debug)]
pub struct Fit {
    pub cand: u32,
    /// Largest distance between one of the read's splice sites and the model's (0: exact).
    pub shift: u16,
    pub off5: i32,
    pub off3: i32,
}

#[derive(Clone, Copy, Debug)]
pub struct HyperOptions {
    pub tol: i64,
    pub slack: i64,
    pub end_slack: i64,
    pub end_bin: i64,
    /// A model's ends are also the sites its reads place, and a read beyond them is not at them.
    pub read_ends: bool,
}

impl Default for HyperOptions {
    fn default() -> Self { Self { tol: 10, slack: 20, end_slack: 100, end_bin: 20, read_ends: true } }
}

/// A read class: reads with one chain, ends within `end_bin` bp and the same end evidence — one weighted hyperedge.
pub struct ReadClass {
    pub chrom: Arc<str>,
    pub plus: bool,
    pub chain: u32,
    /// Mean 5' and 3' positions of its reads.
    pub five: i64,
    pub three: i64,
    /// A poly(A) tail (basecaller length or soft clip over 10 bases) and a capped 5' end.
    pub tail: bool,
    pub cap: bool,
    pub weight: f64,
}

/// Candidate introns by donor and candidate exons, per chromosome and strand, for fitting read chains.
pub struct FitIndex {
    sites: FxHashMap<(Arc<str>, bool), Vec<(i64, u32, u32)>>,
    trees: FxHashMap<(Arc<str>, bool), COITree<(u32, u32), u32>>,
    opts: HyperOptions,
}

pub struct HyperCore {
    pub index: FitIndex,
    pub cands: Vec<Candidate>,
    pub classes: Vec<ReadClass>,
    /// Class of each read of `db.observed`.
    pub class_of: Vec<u32>,
    /// Fits of class i: fits[offsets[i]..offsets[i + 1]].
    pub offsets: Vec<u32>,
    pub fits: Vec<Fit>,
}

fn exon_len(e: (i64, i64)) -> i64 { (e.1 - e.0).abs() + 1 }

impl FitIndex {
    /// The candidates holding a read chain (introns 5'->3') with its 5' and 3' ends.
    pub fn fit(&self, cands: &[Candidate], chrom: &Arc<str>, plus: bool, ch: &[(i64, i64)], five: i64, three: i64) -> Vec<Fit> {
        let key = (chrom.clone(), plus);
        let mut out: Vec<Fit> = Vec::new();
        let push = |f: Fit, out: &mut Vec<Fit>| {
            match out.iter_mut().find(|x| x.cand == f.cand) {
                Some(x) => if f.shift < x.shift { *x = f },
                None => out.push(f),
            }
        };
        if ch.is_empty() {
            let Some(tree) = self.trees.get(&key) else { return out };
            let (lo, hi) = (five.min(three), five.max(three));
            tree.query(lo as i32, hi as i32, |iv| {
                let (ci, j) = iv.metadata;
                let c = &cands[ci as usize];
                let e = c.exons[j as usize];
                let (elo, ehi) = (e.0.min(e.1), e.0.max(e.1));
                let s5 = if j == 0 { self.opts.end_slack } else { self.opts.slack };
                let s3 = if j as usize + 1 == c.exons.len() { self.opts.end_slack } else { self.opts.slack };
                let (lo_s, hi_s) = if plus { (s5, s3) } else { (s3, s5) };
                if lo < elo - lo_s || hi > ehi + hi_s { return; }
                let off5 = c.before[j as usize] + if plus { five - e.0 } else { e.0 - five };
                let off3 = c.len - c.before[j as usize] - exon_len(e) + if plus { e.1 - three } else { three - e.1 };
                push(Fit { cand: ci, shift: 0, off5: off5 as i32, off3: off3 as i32 }, &mut out);
            });
            return out;
        }
        let Some(v) = self.sites.get(&key) else { return out };
        let k = ch.len();
        let d0 = ch[0].0;
        let start = v.partition_point(|&(d, _, _)| d < d0 - self.opts.tol);
        for &(d, ci, j) in &v[start..] {
            if d > d0 + self.opts.tol { break; }
            let c = &cands[ci as usize];
            let j = j as usize;
            if j + k > c.introns.len() { continue; }
            let mut shift = 0i64;
            let mut ok = true;
            for (q, &(rd, ra)) in ch.iter().enumerate() {
                let (md, ma) = c.introns[j + q];
                let s = (md - rd).abs().max((ma - ra).abs());
                if s > self.opts.tol { ok = false; break; }
                shift = shift.max(s);
            }
            if !ok { continue; }
            let (e5, e3) = (c.exons[j], c.exons[j + k]);
            let s5 = if j == 0 { self.opts.end_slack } else { self.opts.slack };
            let s3 = if j + k + 1 == c.exons.len() { self.opts.end_slack } else { self.opts.slack };
            let inside = if plus {
                e5.0 - s5 <= five && five <= e5.1 && e3.0 <= three && three <= e3.1 + s3
            } else {
                e5.1 <= five && five <= e5.0 + s5 && e3.1 - s3 <= three && three <= e3.0
            };
            if !inside { continue; }
            let off5 = c.before[j] + if plus { five - e5.0 } else { e5.0 - five };
            let off3 = c.len - c.before[j + k] - exon_len(e3) + if plus { e3.1 - three } else { three - e3.1 };
            push(Fit { cand: ci, shift: shift.min(u16::MAX as i64) as u16, off5: off5 as i32, off3: off3 as i32 }, &mut out);
        }
        out
    }
}

impl HyperCore {
    /// Candidates are the models any read is counted for; `sites` gives the ends their reads place.
    pub fn build(db: &PostDB, chains: &ReadChains, sites: &Sites, opts: HyperOptions) -> Self {
        let mut counted: FxHashMap<i64, u32> = FxHashMap::default();
        for r in &db.observed { *counted.entry(r.transcript_id).or_insert(0) += 1; }
        let edge_ends: FxHashMap<i64, (i64, i64)> = db.edges.iter().map(|e| (e.edge_id, (e.v1, e.v2))).collect();
        let edge_strand: FxHashMap<i64, bool> = db.edges.iter().map(|e| (e.edge_id, e.strand == "+")).collect();
        let mut cands = Vec::new();
        let mut chroms: FxHashMap<String, Arc<str>> = FxHashMap::default();
        for t in &db.transcripts {
            let Some(&reads) = counted.get(&t.transcript_id) else { continue };
            let Some(&plus) = edge_strand.get(&t.start_exon) else { continue };
            let Some(mut pos) = crate::quant_em::model_positions(db, t, &edge_ends) else { continue };
            // a variant split from an annotated record ends where its reads place it
            if let (Some(&(r5, r3)), true) = (db.annot.variant_ends.get(&t.transcript_id), pos.len() >= 2) {
                let n = pos.len();
                if let Some(x) = r5 { pos[0] = x; }
                if let Some(x) = r3 { pos[n - 1] = x; }
            }
            let exons: Vec<(i64, i64)> = pos.chunks(2).map(|p| (p[0], p[1])).collect();
            if exons.is_empty() { continue; }
            let introns = exons.windows(2).map(|w| (w[0].1, w[1].0)).collect();
            let mut before = Vec::with_capacity(exons.len());
            let mut acc = 0;
            for &e in &exons { before.push(acc); acc += exon_len(e); }
            let chrom = chroms.entry(t.chromosome.clone()).or_insert_with(|| Arc::from(t.chromosome.as_str())).clone();
            // a site counts where it lies in the terminal exon or beyond the chain's end
            let (first, last) = (exons[0], exons[exons.len() - 1]);
            let inward = |site: i64, end: i64, five: bool| if plus == five { site - end } else { end - site };
            let (mut site5, mut site3) = (None, None);
            if opts.read_ends {
                if let Some(&(tss, tes)) = sites.get(&t.transcript_id) {
                    site5 = tss.map(|s| inward(s, first.0, true)).filter(|&o| o < exon_len(first)).map(|o| o as i32);
                    site3 = tes.map(|s| inward(s, last.1, false)).filter(|&o| o < exon_len(last)).map(|o| o as i32);
                }
            }
            cands.push(Candidate { tid: t.transcript_id, gene: t.gene_id, chrom, plus, exons, introns, before, len: acc,
                                   site5, site3, bounded: opts.read_ends, reads });
        }
        // (chromosome, strand) -> introns by donor; and exon intervals for unspliced reads
        let mut sites: FxHashMap<(Arc<str>, bool), Vec<(i64, u32, u32)>> = FxHashMap::default();
        let mut exon_iv: FxHashMap<(Arc<str>, bool), Vec<Interval<(u32, u32)>>> = FxHashMap::default();
        for (ci, c) in cands.iter().enumerate() {
            let key = (c.chrom.clone(), c.plus);
            let v = sites.entry(key.clone()).or_default();
            for (j, &(d, _)) in c.introns.iter().enumerate() { v.push((d, ci as u32, j as u32)); }
            let iv = exon_iv.entry(key).or_default();
            for (j, &e) in c.exons.iter().enumerate() {
                iv.push(Interval::new(e.0.min(e.1) as i32, e.0.max(e.1) as i32, (ci as u32, j as u32)));
            }
        }
        for v in sites.values_mut() { v.sort_unstable(); }
        let trees: FxHashMap<(Arc<str>, bool), COITree<(u32, u32), u32>> =
            exon_iv.into_iter().map(|(k, v)| (k, COITree::new(&v))).collect();

        // read classes
        let mut class_index: FxHashMap<(Arc<str>, bool, u32, i64, i64, bool, bool), u32> = FxHashMap::default();
        let mut classes: Vec<ReadClass> = Vec::new();
        let mut sums: Vec<(i64, i64)> = Vec::new();
        let mut class_of = Vec::with_capacity(db.observed.len());
        for (i, r) in db.observed.iter().enumerate() {
            let plus = &*r.strand == "+";
            let chain = chains.chain_of.get(i).copied().unwrap_or(0);
            let tail = r.polya_length.is_some_and(|x| x > 10) || r.polya_softclip.is_some_and(|x| x > 10);
            let cap = r.cap_signature.unwrap_or(false);
            let key = (r.chromosome.clone(), plus, chain, r.start_pos.div_euclid(opts.end_bin), r.end_pos.div_euclid(opts.end_bin), tail, cap);
            let n = classes.len() as u32;
            let c = *class_index.entry(key).or_insert_with(|| {
                classes.push(ReadClass { chrom: r.chromosome.clone(), plus, chain, five: 0, three: 0, tail, cap, weight: 0.0 });
                sums.push((0, 0));
                n
            });
            let k = c as usize;
            classes[k].weight += 1.0;
            sums[k].0 += r.start_pos; sums[k].1 += r.end_pos;
            class_of.push(c);
        }
        drop(class_index);
        for (c, (a, b)) in classes.iter_mut().zip(sums) {
            c.five = (a as f64 / c.weight).round() as i64;
            c.three = (b as f64 / c.weight).round() as i64;
        }

        let index = FitIndex { sites, trees, opts };
        let per_class: Vec<Vec<Fit>> = classes.par_iter()
            .map(|cl| index.fit(&cands, &cl.chrom, cl.plus, &chains.chains[cl.chain as usize], cl.five, cl.three)).collect();
        let mut offsets = Vec::with_capacity(per_class.len() + 1);
        let mut fits = Vec::with_capacity(per_class.iter().map(|v| v.len()).sum());
        offsets.push(0u32);
        for v in per_class {
            fits.extend(v);
            offsets.push(fits.len() as u32);
        }
        Self { cands, classes, class_of, offsets, fits, index }
    }

    pub fn fits_of(&self, class: usize) -> &[Fit] {
        &self.fits[self.offsets[class] as usize..self.offsets[class + 1] as usize]
    }

    /// Writes every fit (class, model, largest shift, 5' and 3' offsets) to a Parquet file, the classes to
    /// `<path>.classes.parquet`, the candidates to `<path>.models.parquet` and each read's class to
    /// `<path>.reads.parquet`.
    pub fn write_fits(&self, path: &Path, chains: &ReadChains) -> IsoForgeResult<usize> {
        {
            let schema = Arc::new(Schema::new(vec![
                Field::new("class", DataType::UInt32, false),
                Field::new("chrom", DataType::Utf8, false),
                Field::new("plus", DataType::Boolean, false),
                Field::new("n_introns", DataType::Int32, false),
                Field::new("chain", DataType::UInt32, false),
                Field::new("five", DataType::Int64, false),
                Field::new("three", DataType::Int64, false),
                Field::new("tail", DataType::Boolean, false),
                Field::new("cap", DataType::Boolean, false),
                Field::new("weight", DataType::Float64, false),
            ]));
            let (mut cb, mut chb, mut pb, mut nb, mut kb, mut fb, mut tb, mut tlb, mut capb, mut wb) = (
                UInt32Builder::new(), arrow::array::StringBuilder::new(), arrow::array::BooleanBuilder::new(), Int32Builder::new(),
                UInt32Builder::new(), Int64Builder::new(), Int64Builder::new(), arrow::array::BooleanBuilder::new(),
                arrow::array::BooleanBuilder::new(), arrow::array::Float64Builder::new());
            for (i, c) in self.classes.iter().enumerate() {
                cb.append_value(i as u32); chb.append_value(&*c.chrom); pb.append_value(c.plus);
                nb.append_value(chains.chains[c.chain as usize].len() as i32); kb.append_value(c.chain);
                fb.append_value(c.five); tb.append_value(c.three); tlb.append_value(c.tail); capb.append_value(c.cap);
                wb.append_value(c.weight);
            }
            let batch = RecordBatch::try_new(schema.clone(), vec![
                Arc::new(cb.finish()), Arc::new(chb.finish()), Arc::new(pb.finish()), Arc::new(nb.finish()), Arc::new(kb.finish()),
                Arc::new(fb.finish()), Arc::new(tb.finish()), Arc::new(tlb.finish()), Arc::new(capb.finish()), Arc::new(wb.finish()),
            ])?;
            let mut w = ArrowWriter::try_new(File::create(path.with_extension("classes.parquet"))?, schema, None)?;
            w.write(&batch)?;
            w.close()?;
            let schema = Arc::new(Schema::new(vec![Field::new("class", DataType::UInt32, false)]));
            let mut b = UInt32Builder::new();
            for &c in &self.class_of { b.append_value(c); }
            let batch = RecordBatch::try_new(schema.clone(), vec![Arc::new(b.finish())])?;
            let mut w = ArrowWriter::try_new(File::create(path.with_extension("reads.parquet"))?, schema, None)?;
            w.write(&batch)?;
            w.close()?;
        }
        let schema = Arc::new(Schema::new(vec![
            Field::new("class", DataType::UInt32, false),
            Field::new("model", DataType::Int64, false),
            Field::new("shift", DataType::Int32, false),
            Field::new("off5", DataType::Int32, false),
            Field::new("off3", DataType::Int32, false),
        ]));
        let mut w = ArrowWriter::try_new(File::create(path)?, schema.clone(), None)?;
        let n_classes = self.offsets.len() - 1;
        const CHUNK: usize = 1 << 18;
        let mut start = 0;
        while start < n_classes {
            let end = (start + CHUNK).min(n_classes);
            let (mut rb, mut mb, mut sb, mut b5, mut b3) =
                (UInt32Builder::new(), Int64Builder::new(), Int32Builder::new(), Int32Builder::new(), Int32Builder::new());
            for i in start..end {
                for f in self.fits_of(i) {
                    rb.append_value(i as u32);
                    mb.append_value(self.cands[f.cand as usize].tid);
                    sb.append_value(f.shift as i32);
                    b5.append_value(f.off5);
                    b3.append_value(f.off3);
                }
            }
            let batch = RecordBatch::try_new(schema.clone(), vec![
                Arc::new(rb.finish()), Arc::new(mb.finish()), Arc::new(sb.finish()), Arc::new(b5.finish()), Arc::new(b3.finish()),
            ])?;
            w.write(&batch)?;
            start = end;
        }
        w.close()?;
        // and the candidates: model, gene, exons, exonic length, reads counted for it
        let schema = Arc::new(Schema::new(vec![
            Field::new("model", DataType::Int64, false),
            Field::new("gene", DataType::Int64, false),
            Field::new("n_exons", DataType::Int32, false),
            Field::new("length", DataType::Int64, false),
            Field::new("reads", DataType::UInt32, false),
        ]));
        let (mut tb, mut gb, mut nb, mut lb, mut rb) =
            (Int64Builder::new(), Int64Builder::new(), Int32Builder::new(), Int64Builder::new(), UInt32Builder::new());
        for c in &self.cands {
            tb.append_value(c.tid); gb.append_value(c.gene); nb.append_value(c.exons.len() as i32);
            lb.append_value(c.len); rb.append_value(c.reads);
        }
        let batch = RecordBatch::try_new(schema.clone(), vec![
            Arc::new(tb.finish()), Arc::new(gb.finish()), Arc::new(nb.finish()), Arc::new(lb.finish()), Arc::new(rb.finish()),
        ])?;
        let mpath = path.with_extension("models.parquet");
        let mut w = ArrowWriter::try_new(File::create(&mpath)?, schema, None)?;
        w.write(&batch)?;
        w.close()?;
        Ok(self.fits.len())
    }
}

/// How reads' ends fall on the models they fit, measured on the run from exact fits of tailed, spliced read classes
/// that end at the model's 3' end: the share whose 5' end is at the model's start, by model length (`p_full`), and
/// the lengths the others cover (`dens`, per base, in bins of `COV_BIN`).
pub struct EndModel {
    p_full: [f64; 8],
    dens: Vec<f64>,
}

const LEN_BINS: [i64; 8] = [500, 1000, 1500, 2000, 3000, 4000, 6000, i64::MAX];
const COV_BIN: i64 = 100;
const COV_BINS: usize = 200;

fn len_bin(l: i64) -> usize { LEN_BINS.iter().position(|&b| l < b).unwrap_or(LEN_BINS.len() - 1) }

impl EndModel {
    /// Measured on the fits to the candidates `on` selects.
    pub fn measure(core: &HyperCore, chains: &ReadChains, on: impl Fn(usize) -> bool) -> Self {
        let mut full = [(0.0f64, 0.0f64); 8];
        let mut hist = vec![1.0f64; COV_BINS];
        for (ci, cl) in core.classes.iter().enumerate() {
            if !cl.tail || chains.chains[cl.chain as usize].is_empty() { continue; }
            for f in core.fits_of(ci) {
                if !on(f.cand as usize) || f.shift != 0 { continue; }
                let c = &core.cands[f.cand as usize];
                if c.fit3(f) != EndFit::At { continue; }
                match c.fit5(f) {
                    EndFit::At => { let b = len_bin(c.span()); full[b].1 += cl.weight; full[b].0 += cl.weight; }
                    EndFit::Short => {
                        full[len_bin(c.span())].1 += cl.weight;
                        let cov = (c.len - f.off5 as i64 - f.off3 as i64).max(0);
                        hist[((cov / COV_BIN) as usize).min(COV_BINS - 1)] += cl.weight;
                    }
                    EndFit::Beyond => {}
                }
            }
        }
        let mut p_full = [0.5f64; 8];
        for (p, &(a, n)) in p_full.iter_mut().zip(&full) { if n > 0.0 { *p = (a / n).clamp(0.02, 0.98); } }
        let tot: f64 = hist.iter().sum();
        Self { p_full, dens: hist.iter().map(|&h| h / tot / COV_BIN as f64).collect() }
    }

    /// Likelihood of a read class (tailed or not) under one model it fits.
    pub fn lik(&self, core: &HyperCore, tail: bool, f: &Fit) -> f64 {
        let c = &core.cands[f.cand as usize];
        let pf = self.p_full[len_bin(c.span())];
        let at = pf / (END_WINDOW as f64 + 100.0);
        let l5 = match c.fit5(f) {
            EndFit::At => at,
            EndFit::Short => {
                let cov = (c.len - f.off5 as i64 - f.off3 as i64).max(0);
                (1.0 - pf) * self.dens[((cov / COV_BIN) as usize).min(COV_BINS - 1)]
            }
            EndFit::Beyond => at * OFF_END,
        };
        let l3 = if tail && c.fit3(f) != EndFit::At { OFF_END } else { 1.0 };
        l5 * l3 * if f.shift == 0 { 1.0 } else { SHIFTED }
    }
}

/// Counts per reported model and dataset. A read is counted for its own reported model when `stays` holds, its class
/// fits no more than one reported model or it is a whole read of that model (`Candidate::whole`); the other eligible
/// reads are shared by EM among the reported models their class fits. None when the reads' chains are absent.
pub fn hypergraph_counts(db: &PostDB, reported: &rustc_hash::FxHashSet<i64>, eligible: impl Fn(&crate::post_processing::ObservedRead) -> bool,
                         stays: impl Fn(&crate::post_processing::ObservedRead) -> bool,
                         quant_reads: &[crate::quant_em::QuantRead], quant_eligible: impl Fn(&crate::quant_em::QuantRead) -> bool,
                         sites: &Sites, opts: HyperOptions) -> IsoForgeResult<Option<FxHashMap<i64, FxHashMap<String, f64>>>> {
    let Some(chains) = load_read_chains(&db.observed_path)? else { return Ok(None) };
    if chains.chain_of.len() != db.observed.len() { return Ok(None); }
    let core = HyperCore::build(db, &chains, sites, opts);
    let rep_cand: Vec<bool> = core.cands.iter().map(|c| reported.contains(&c.tid)).collect();
    // Annotated small RNAs are counted for the reads they hold and are no fit for a read class.
    let small_cand: Vec<bool> = core.cands.iter().map(|c| db.annot.small_rna.contains_key(&c.tid)).collect();
    let rep_fits = |c: usize| core.fits_of(c).iter().filter(|f| rep_cand[f.cand as usize] && !small_cand[f.cand as usize]);
    let ends = EndModel::measure(&core, &chains, |c| rep_cand[c]);
    let lik = |tail: bool, f: &Fit| ends.lik(&core, tail, f);
    let slot: FxHashMap<i64, usize> = core.cands.iter().enumerate().filter(|(i, _)| rep_cand[*i]).map(|(i, c)| (c.tid, i)).collect();
    let ds_of: FxHashMap<&str, usize> = db.datasets.iter().enumerate().map(|(i, d)| (d.as_str(), i)).collect();
    let nd = db.datasets.len();
    let mut fixed: FxHashMap<usize, Vec<f64>> = FxHashMap::default();
    let mut shared: FxHashMap<usize, Vec<f64>> = FxHashMap::default();
    // ISOFORGE_HQ_TRACE=<model,model,...> logs what became of those models' reads
    let trace: Vec<i64> = std::env::var("ISOFORGE_HQ_TRACE").ok().map(|v| v.split(',').filter_map(|x| x.parse().ok()).collect()).unwrap_or_default();
    let mut traced: FxHashMap<i64, [u32; 7]> = FxHashMap::default();
    let mut traced_classes: FxHashMap<i64, Vec<usize>> = FxHashMap::default();
    // reads kept with their model: [validated whole reads, for fitting no other model, whole by chain and ends]
    let mut kept = [0u64; 3];
    for (i, r) in db.observed.iter().enumerate() {
        let Some(&di) = ds_of.get(&*r.dataset) else { continue };
        if trace.contains(&r.transcript_id) {
            let t = traced.entry(r.transcript_id).or_default();
            t[0] += 1;
            if eligible(r) {
                t[1] += 1;
                let c = core.class_of[i] as usize;
                let n = rep_fits(c).count();
                if n == 0 { t[2] += 1 } else if n == 1 { t[3] += 1 } else { t[4] += 1 }
                if let Some(&m) = slot.get(&r.transcript_id) {
                    if rep_fits(c).any(|f| f.cand as usize == m && core.cands[m].whole(f)) { t[5] += 1 }
                    if stays(r) { t[6] += 1 }
                }
                let v = traced_classes.entry(r.transcript_id).or_default();
                if !v.contains(&c) && v.len() < 4 { v.push(c) }
            }
        }
        if !eligible(r) { continue; }
        let c = core.class_of[i] as usize;
        let fits: Vec<&Fit> = rep_fits(c).collect();
        // a read assigned to a small RNA that does not hold it is counted like a removed model's read
        let small_rna = db.annot.small_rna.get(&r.transcript_id);
        let held = small_rna.is_some_and(|s| s.holds(r.start_pos, r.end_pos));
        let own = slot.get(&r.transcript_id).copied().filter(|_| small_rna.is_none() || held);
        let to = match own {
            Some(m) if held => Some(m),
            Some(m) => {
                if stays(r) { kept[0] += 1; Some(m) }
                else if fits.len() <= 1 { kept[1] += 1; Some(m) }
                else if fits.iter().any(|f| f.cand as usize == m && core.cands[m].whole(f)) { kept[2] += 1; Some(m) }
                else { None }
            }
            None => {
                if fits.is_empty() { continue; }
                if fits.len() == 1 { Some(fits[0].cand as usize) } else { None }
            }
        };
        match to {
            Some(m) => fixed.entry(m).or_insert_with(|| vec![0.0; nd])[di] += 1.0,
            None => shared.entry(c).or_insert_with(|| vec![0.0; nd])[di] += 1.0,
        }
    }
    // quantification-only reads (below the mapping quality that builds models): fitted like a removed model's read
    let mut extra: Vec<(Vec<f64>, Vec<(usize, f64)>)> = Vec::new();
    let mut chain_buf: Vec<(i64, i64)> = Vec::new();
    for q in quant_reads {
        let Some(&di) = ds_of.get(&*q.dataset) else { continue };
        if !quant_eligible(q) || q.positions.len() < 2 { continue; }
        let p = &q.positions;
        chain_buf.clear();
        chain_buf.extend(p[1..p.len() - 1].chunks_exact(2).map(|x| (x[0], x[1])));
        let plus = &*q.strand == "+";
        let fits: Vec<Fit> = core.index.fit(&core.cands, &q.chromosome, plus, &chain_buf, p[0], p[p.len() - 1])
            .into_iter().filter(|f| rep_cand[f.cand as usize] && !small_cand[f.cand as usize]).collect();
        match fits.len() {
            0 => {}
            1 => fixed.entry(fits[0].cand as usize).or_insert_with(|| vec![0.0; nd])[di] += 1.0,
            _ => {
                let mut w = vec![0.0; nd];
                w[di] = 1.0;
                extra.push((w, fits.iter().map(|f| (f.cand as usize, lik(false, f))).collect()));
            }
        }
    }
    // EM over the shared classes, abundances started from the reads that stay
    let mut sh: Vec<(Vec<f64>, Vec<(usize, f64)>)> = shared.iter().map(|(&c, w)| {
        let cl = &core.classes[c];
        (w.clone(), rep_fits(c).map(|f| (f.cand as usize, lik(cl.tail, f))).collect())
    }).collect();
    let n_shared_classes = sh.len();
    sh.extend(extra);
    let n_quant_shared = sh.len() - n_shared_classes;
    let stay: FxHashMap<usize, f64> = slot.values().map(|&m| (m, fixed.get(&m).map_or(0.0, |v| v.iter().sum::<f64>()))).collect();
    let weights: Vec<f64> = sh.iter().map(|(w, _)| w.iter().sum()).collect();
    let theta = em_abundance(&stay, &weights, &sh.iter().map(|(_, fl)| &fl[..]).collect::<Vec<_>>());
    let mut out: FxHashMap<i64, FxHashMap<String, f64>> = FxHashMap::default();
    for (&m, v) in &fixed {
        let e = out.entry(core.cands[m].tid).or_default();
        for (d, &x) in v.iter().enumerate() { if x > 0.0 { *e.entry(db.datasets[d].clone()).or_insert(0.0) += x; } }
    }
    for (w, fl) in &sh {
        let tot: f64 = fl.iter().map(|&(m, l)| theta[&m] * l).sum();
        if tot <= 0.0 { continue; }
        for &(m, l) in fl {
            let share = theta[&m] * l / tot;
            let e = out.entry(core.cands[m].tid).or_default();
            for (d, &x) in w.iter().enumerate() { if x > 0.0 { *e.entry(db.datasets[d].clone()).or_insert(0.0) += x * share; } }
        }
    }
    for &t in &trace {
        let cand = core.cands.iter().position(|c| c.tid == t);
        log::info!("[HQ trace] model {}: candidate {:?}, reported {}, slot {:?}, reads [all, eligible, fit 0, fit 1, fit several, whole by fit, validated whole] {:?}, theta {:?}, out {:?}",
                   t, cand, reported.contains(&t), slot.get(&t), traced.get(&t), cand.and_then(|m| theta.get(&m)), out.get(&t));
        if let Some(m) = cand { let c = &core.cands[m]; log::info!("[HQ trace]   model exons {:?} len {} sites {:?} {:?}", c.exons, c.len, c.site5, c.site3); }
        for &c in traced_classes.get(&t).map(|v| &v[..]).unwrap_or(&[]) {
            let cl = &core.classes[c];
            log::info!("[HQ trace]   class {} weight {} five {} three {} tail {} chain {:?}", c, cl.weight, cl.five, cl.three, cl.tail, chains.chains[cl.chain as usize]);
            for f in core.fits_of(c) {
                log::info!("[HQ trace]     fit model {} reported {} shift {} off5 {} off3 {} lik {:.3e} theta {:?}", core.cands[f.cand as usize].tid,
                           rep_cand[f.cand as usize], f.shift, f.off5, f.off3, lik(cl.tail, f), theta.get(&(f.cand as usize)));
            }
        }
    }
    log::info!("[Post] Hypergraph quantification: {} read classes, {} shared among several reported models, {} quantification-only reads shared",
               core.classes.len(), n_shared_classes, n_quant_shared);
    log::info!("[Post] Hypergraph quantification: reads kept with their reported model: {} whole reads with validated ends, {} that fit no other model, {} whole by their chain and ends",
               kept[0], kept[1], kept[2]);
    Ok(Some(out))
}

/// Abundances by EM: `stay` holds each model's reads that are its own, `shared[i]` the models a shared read class fits
/// with the class's likelihood under each, `weights[i]` its reads. Each round shares every class among its models by
/// abundance times likelihood.
fn em_abundance(stay: &FxHashMap<usize, f64>, weights: &[f64], shared: &[&[(usize, f64)]]) -> FxHashMap<usize, f64> {
    let mut theta: FxHashMap<usize, f64> = stay.iter().map(|(&m, &n)| (m, n + 0.5)).collect();
    for _ in 0..200 {
        let mut add: FxHashMap<usize, f64> = FxHashMap::default();
        for (&w, fl) in weights.iter().zip(shared) {
            let tot: f64 = fl.iter().map(|&(m, l)| theta[&m] * l).sum();
            if tot <= 0.0 { continue; }
            for &(m, l) in fl.iter() { *add.entry(m).or_insert(0.0) += w * theta[&m] * l / tot; }
        }
        let mut delta = 0.0;
        let mut total = 0.0;
        for (&m, th) in theta.iter_mut() {
            let new = stay[&m] + add.get(&m).copied().unwrap_or(0.0);
            delta += (new - *th).abs();
            total += new;
            *th = new.max(1e-9);
        }
        if delta < 1e-4 * total.max(1.0) { break; }
    }
    theta
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shared_reads_follow_the_likelihood_and_then_the_abundance() {
        // model 0 keeps 100 reads of its own, model 1 none; a class of 50 reads fits both equally: it goes to model 0
        let stay: FxHashMap<usize, f64> = [(0, 100.0), (1, 0.0)].into_iter().collect();
        let th = em_abundance(&stay, &[50.0], &[&[(0, 1e-3), (1, 1e-3)]]);
        assert!(th[&1] < 0.5 && (th[&0] - 150.0).abs() < 0.5, "{th:?}");
        // the same class a hundred times likelier under model 1, which also keeps reads: model 1 takes it
        let stay: FxHashMap<usize, f64> = [(0, 100.0), (1, 20.0)].into_iter().collect();
        let th = em_abundance(&stay, &[50.0], &[&[(0, 1e-5), (1, 1e-3)]]);
        assert!(th[&1] > 65.0 && (th[&0] + th[&1] - 170.0).abs() < 1e-6, "{th:?}");
        // two models with nothing of their own and one equal class: split, no read lost
        let stay: FxHashMap<usize, f64> = [(0, 0.0), (1, 0.0)].into_iter().collect();
        let th = em_abundance(&stay, &[10.0], &[&[(0, 1e-3), (1, 1e-3)]]);
        assert!((th[&0] - 5.0).abs() < 1e-6 && (th[&1] - 5.0).abs() < 1e-6, "{th:?}");
    }

    #[test]
    fn a_reads_end_is_at_the_chains_end_or_at_the_site_the_reads_place() {
        // the chain's own end alone, unbounded: up to the window inside it, however far beyond
        assert_eq!(end_fit(-90, None, false), EndFit::At);
        assert_eq!(end_fit(51, None, false), EndFit::Short);
        // bounded: past the window on the outside the read is beyond the model
        assert_eq!(end_fit(-50, None, true), EndFit::At);
        assert_eq!(end_fit(-51, None, true), EndFit::Beyond);
        // reads place the start 69 bp inside the annotated one: both are ends of the model
        assert_eq!(end_fit(69, Some(69), true), EndFit::At);
        assert_eq!(end_fit(0, Some(69), true), EndFit::At);
        assert_eq!(end_fit(130, Some(69), true), EndFit::Short);
        // a site 1.9 kb inside: a read between the two ends is short of the outer one
        assert_eq!(end_fit(900, Some(1898), true), EndFit::Short);
        assert_eq!(end_fit(1900, Some(1898), true), EndFit::At);
        // a site upstream of the chain's end
        assert_eq!(end_fit(-300, Some(-300), true), EndFit::At);
        assert_eq!(end_fit(-150, Some(-300), true), EndFit::Short);
        assert_eq!(end_fit(-400, Some(-300), true), EndFit::Beyond);
    }

    fn cand(len: i64, site5: Option<i32>, bounded: bool) -> Candidate {
        Candidate { tid: 0, gene: 0, chrom: Arc::from("chr1"), plus: true, exons: vec![(1, len)], introns: Vec::new(), before: vec![0], len,
                    site5, site3: None, bounded, reads: 0 }
    }

    #[test]
    fn a_read_at_the_site_of_an_annotated_model_is_not_given_to_a_shorter_variant() {
        // an annotated model whose reads start 69 bp inside its annotated start, and an end variant 65 bp further in
        let core = |bounded: bool| HyperCore {
            index: FitIndex { sites: FxHashMap::default(), trees: FxHashMap::default(), opts: HyperOptions::default() },
            cands: vec![cand(337, if bounded { Some(69) } else { None }, bounded), cand(272, None, bounded)],
            classes: Vec::new(), class_of: Vec::new(), offsets: vec![0], fits: Vec::new(),
        };
        let ends = EndModel { p_full: [0.9; 8], dens: vec![1e-3; COV_BINS] };
        let (own, variant) = (Fit { cand: 0, shift: 0, off5: 69, off3: 0 }, Fit { cand: 1, shift: 0, off5: -65, off3: 0 });
        // by the chains' ends alone the read is truncated for its own model and whole for the variant
        let c = core(false);
        assert!(!c.cands[0].whole(&own) && c.cands[1].whole(&variant));
        assert!(ends.lik(&c, true, &variant) > 10.0 * ends.lik(&c, true, &own));
        // with the site its reads place it is a whole read of its own model, and starts upstream of the variant
        let c = core(true);
        assert!(c.cands[0].whole(&own) && !c.cands[1].whole(&variant));
        assert!(ends.lik(&c, true, &own) > 10.0 * ends.lik(&c, true, &variant));
        assert_eq!(c.cands[0].span(), 337 - 69);
        // a read ending past the model's 3' end is beyond it, and less likely with a tail than without
        let past = Fit { cand: 1, shift: 0, off5: 0, off3: -80 };
        assert_eq!(c.cands[1].fit3(&past), EndFit::Beyond);
        assert!(ends.lik(&c, true, &past) < 0.1 * ends.lik(&c, false, &past));
    }

    #[test]
    fn exon_lengths_count_both_ends() {
        assert_eq!(exon_len((100, 199)), 100);
        assert_eq!(exon_len((199, 100)), 100);
    }
}
