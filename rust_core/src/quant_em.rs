// Quantification: reads of models the filters removed, and reads kept for quantification only (MAPQ below the
// model-building threshold), are counted towards the reported models they fit, shared by EM.
//
// - A read fits a reported model when each of its introns lies within `tol` bp of one of the model's at both sites,
//   the matched introns are contiguous and in order in the model, and its ends stay inside the flanking exons (within
//   `slack` bp, `end_slack` beyond the model's own ends); an unspliced read, when one exon holds it (same slacks).
// - Without `partial_reads`, a read must also span the model, both ends within `end_window` bp of the model's own.
// - Reads that fit several models are shared in proportion to the models' abundance.

use std::path::Path;
use std::sync::Arc;

use arrow::array::{Array, BooleanArray, Float64Array, StringArray};
use coitrees::{COITree, Interval, IntervalTree};
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use rustc_hash::FxHashMap;

use crate::errors::IsoForgeResult;
use crate::post_processing::{PostDB, TranscriptRecord};

#[derive(Debug, Clone, Copy)]
pub struct QuantOptions {
    pub tol: i64,
    pub slack: i64,
    pub end_slack: i64,
    pub max_iter: usize,
    /// Also count reads that do not span the model they fit; off, both ends must be within `end_window` of the model's.
    pub partial_reads: bool,
    pub end_window: i64,
}

impl Default for QuantOptions {
    fn default() -> Self {
        Self { tol: 10, slack: 20, end_slack: 100, max_iter: 200, partial_reads: false, end_window: 50 }
    }
}

/// A read kept for quantification only, as written by `QuantReadWriter`.
pub struct QuantRead {
    pub dataset: Arc<str>,
    pub chromosome: Arc<str>,
    pub strand: Arc<str>,
    /// Exon boundaries in 5'->3' order.
    pub positions: Vec<i64>,
    pub fraction_as: Option<f64>,
    pub has_pas_signal: Option<bool>,
}

pub fn load_quant_reads(path: &Path) -> IsoForgeResult<Vec<QuantRead>> {
    let mut out = Vec::new();
    if !path.exists() {
        return Ok(out);
    }
    let mut strings: FxHashMap<String, Arc<str>> = FxHashMap::default();
    let mut intern = |s: &str| -> Arc<str> {
        if let Some(a) = strings.get(s) { return Arc::clone(a); }
        let a: Arc<str> = Arc::from(s);
        strings.insert(s.to_string(), Arc::clone(&a));
        a
    };
    let reader = ParquetRecordBatchReaderBuilder::try_new(std::fs::File::open(path)?)?.build()?;
    for batch in reader {
        let batch = batch?;
        let col = |name: &str| batch.column(batch.schema().index_of(name).unwrap()).clone();
        let ds = col("dataset"); let ds = ds.as_any().downcast_ref::<StringArray>().unwrap();
        let ch = col("chromosome"); let ch = ch.as_any().downcast_ref::<StringArray>().unwrap();
        let st = col("strand"); let st = st.as_any().downcast_ref::<StringArray>().unwrap();
        let po = col("positions"); let po = po.as_any().downcast_ref::<StringArray>().unwrap();
        let fa = col("fraction_as"); let fa = fa.as_any().downcast_ref::<Float64Array>().unwrap();
        let pa = col("has_pas_signal"); let pa = pa.as_any().downcast_ref::<BooleanArray>().unwrap();
        for i in 0..batch.num_rows() {
            let positions: Vec<i64> = po.value(i).split(',').filter_map(|x| x.parse().ok()).collect();
            if positions.len() < 2 { continue; }
            out.push(QuantRead {
                dataset: intern(ds.value(i)),
                chromosome: intern(ch.value(i)),
                strand: intern(st.value(i)),
                positions,
                fraction_as: if fa.is_null(i) { None } else { Some(fa.value(i)) },
                has_pas_signal: if pa.is_null(i) { None } else { Some(pa.value(i)) },
            });
        }
    }
    Ok(out)
}

/// A reported model in transcript-oriented coordinates (x = position on +, -position on -), so that
/// x increases 5' -> 3' on both strands.
struct Model {
    exons: Vec<(i64, i64)>,
    introns: Vec<(i64, i64)>,
}

/// The introns of a chain in oriented coordinates (donor, acceptor), from exon boundaries in 5'->3' order.
fn oriented_introns(positions: &[i64], plus: bool) -> Vec<(i64, i64)> {
    let s = if plus { 1 } else { -1 };
    positions[1..positions.len() - 1].chunks(2).map(|p| (s * p[0], s * p[1])).collect()
}

/// Reported models indexed for compatibility queries.
pub struct ModelIndex {
    ids: Vec<i64>,
    models: Vec<Model>,
    /// (chromosome, plus strand) -> distinct introns sorted by donor, acceptor
    introns: FxHashMap<(Arc<str>, bool), Vec<(i64, i64)>>,
    /// (chromosome, plus, donor, acceptor) -> (model, intron index)
    at: FxHashMap<(Arc<str>, bool, i64, i64), Vec<(u32, u32)>>,
    /// (chromosome, plus) -> exons as oriented intervals, metadata (model, exon index)
    exons: FxHashMap<(Arc<str>, bool), COITree<(u32, u32), u32>>,
    opts: QuantOptions,
}

/// Exon boundaries of a model in 5'->3' order: its hyperedge's vertex chain (PostDB::chains), or from its edge path
/// when the database holds no chain for it; None when a vertex is unknown.
pub fn model_positions(db: &PostDB, t: &TranscriptRecord, edge_ends: &FxHashMap<i64, (i64, i64)>) -> Option<Vec<i64>> {
    if let Some(c) = db.chains.get(&t.transcript_id) {
        return Some(c.clone());
    }
    crate::post_processing::chain_from_edges(&crate::post_processing::record_path(t), edge_ends, &db.loc_map)
}

impl ModelIndex {
    pub fn build(db: &PostDB, reported: &[i64], opts: QuantOptions) -> Self {
        let edge_ends: FxHashMap<i64, (i64, i64)> = db.edges.iter().map(|e| (e.edge_id, (e.v1, e.v2))).collect();
        let edge_strand: FxHashMap<i64, bool> = db.edges.iter().map(|e| (e.edge_id, e.strand == "+")).collect();
        let mut ids = Vec::new();
        let mut models = Vec::new();
        let mut introns: FxHashMap<(Arc<str>, bool), Vec<(i64, i64)>> = FxHashMap::default();
        let mut at: FxHashMap<(Arc<str>, bool, i64, i64), Vec<(u32, u32)>> = FxHashMap::default();
        let mut ex_iv: FxHashMap<(Arc<str>, bool), Vec<Interval<(u32, u32)>>> = FxHashMap::default();
        let mut chroms: FxHashMap<String, Arc<str>> = FxHashMap::default();
        for &tid in reported {
            let Some(&ti) = db.tx_map.get(&tid) else { continue };
            let t = &db.transcripts[ti];
            let Some(&plus) = edge_strand.get(&t.start_exon) else { continue };
            let Some(pos) = model_positions(db, t, &edge_ends) else { continue };
            let s = if plus { 1 } else { -1 };
            let exons: Vec<(i64, i64)> = pos.chunks(2).map(|p| { let (a, b) = (s * p[0], s * p[1]); (a.min(b), a.max(b)) }).collect();
            let ints: Vec<(i64, i64)> = exons.windows(2).map(|w| (w[0].1, w[1].0)).collect();
            let chrom = chroms.entry(t.chromosome.clone()).or_insert_with(|| Arc::from(t.chromosome.as_str())).clone();
            let m = models.len() as u32;
            for (j, &(d, a)) in ints.iter().enumerate() {
                introns.entry((Arc::clone(&chrom), plus)).or_default().push((d, a));
                at.entry((Arc::clone(&chrom), plus, d, a)).or_default().push((m, j as u32));
            }
            let n = exons.len();
            for (j, &(lo, hi)) in exons.iter().enumerate() {
                let lo_s = if j == 0 { opts.end_slack } else { opts.slack };
                let hi_s = if j + 1 == n { opts.end_slack } else { opts.slack };
                ex_iv.entry((Arc::clone(&chrom), plus)).or_default()
                    .push(Interval::new((lo - lo_s) as i32, (hi + hi_s) as i32, (m, j as u32)));
            }
            ids.push(tid);
            models.push(Model { exons, introns: ints });
        }
        for v in introns.values_mut() { v.sort_unstable(); v.dedup(); }
        let exons = ex_iv.into_iter().map(|(k, v)| (k, COITree::new(&v))).collect();
        Self { ids, models, introns, at, exons, opts }
    }

    pub fn transcript_id(&self, m: u32) -> i64 { self.ids[m as usize] }

    /// The reported intron nearest to (d, a) with both sites within `tol`.
    fn snap(&self, key: &(Arc<str>, bool), d: i64, a: i64) -> Option<(i64, i64)> {
        let v = self.introns.get(key)?;
        let tol = self.opts.tol;
        let i0 = v.partition_point(|&(x, _)| x < d - tol);
        let mut best: Option<(i64, (i64, i64))> = None;
        for &(x, y) in &v[i0..] {
            if x > d + tol { break; }
            if (y - a).abs() > tol { continue; }
            let dist = (x - d).abs() + (y - a).abs();
            if best.map_or(true, |(b, _)| dist < b) { best = Some((dist, (x, y))); }
        }
        best.map(|(_, iv)| iv)
    }

    /// Candidate (model, first matched intron) pairs for an oriented intron chain, before end checks.
    pub fn chain_candidates(&self, chrom: &Arc<str>, plus: bool, chain: &[(i64, i64)]) -> Vec<(u32, u32)> {
        let key = (Arc::clone(chrom), plus);
        let mut snapped = Vec::with_capacity(chain.len());
        for &(d, a) in chain {
            match self.snap(&key, d, a) {
                Some(iv) => snapped.push(iv),
                None => return Vec::new(),
            }
        }
        let Some(first) = self.at.get(&(Arc::clone(chrom), plus, snapped[0].0, snapped[0].1)) else { return Vec::new() };
        first.iter().copied().filter(|&(m, j)| {
            let ints = &self.models[m as usize].introns;
            let j = j as usize;
            j + snapped.len() <= ints.len() && ints[j..j + snapped.len()] == snapped[..]
        }).collect()
    }

    /// Whether a read with oriented ends x5 <= x3 and `k` introns fits model m at intron offset j: with
    /// `partial_reads`, inside the flanking exons; otherwise spanning the model, both ends within `end_window`
    /// of the model's own.
    pub fn ends_fit(&self, m: u32, j: u32, k: usize, x5: i64, x3: i64) -> bool {
        let ex = &self.models[m as usize].exons;
        let (j, last) = (j as usize, j as usize + k);
        if !self.opts.partial_reads {
            let w = self.opts.end_window;
            return j == 0 && last + 1 == ex.len() && (x5 - ex[0].0).abs() <= w && (x3 - ex[last].1).abs() <= w;
        }
        let lo_s = if j == 0 { self.opts.end_slack } else { self.opts.slack };
        let hi_s = if last + 1 == ex.len() { self.opts.end_slack } else { self.opts.slack };
        x5 >= ex[j].0 - lo_s && x3 <= ex[last].1 + hi_s
    }

    /// Models holding the oriented span [x5, x3] inside one exon (within the slacks); without `partial_reads`,
    /// only single-exon models the read spans, both ends within `end_window`.
    pub fn exon_candidates(&self, chrom: &Arc<str>, plus: bool, x5: i64, x3: i64, out: &mut Vec<u32>) {
        let Some(tree) = self.exons.get(&(Arc::clone(chrom), plus)) else { return };
        let (partial, w) = (self.opts.partial_reads, self.opts.end_window);
        tree.query(x5 as i32, x3 as i32, |iv| {
            if (iv.first as i64) > x5 || (iv.last as i64) < x3 { return; }
            let m = iv.metadata.0;
            if !partial {
                let ex = &self.models[m as usize].exons;
                if ex.len() != 1 || (x5 - ex[0].0).abs() > w || (x3 - ex[0].1).abs() > w { return; }
            }
            out.push(m);
        });
    }

    /// Compatible models of a read given as exon boundaries in 5'->3' order.
    pub fn read_candidates(&self, chrom: &Arc<str>, plus: bool, positions: &[i64], out: &mut Vec<u32>) {
        out.clear();
        let s = if plus { 1 } else { -1 };
        let (x5, x3) = { let (a, b) = (s * positions[0], s * positions[positions.len() - 1]); (a.min(b), a.max(b)) };
        if positions.len() <= 2 {
            self.exon_candidates(chrom, plus, x5, x3, out);
        } else {
            let chain = oriented_introns(positions, plus);
            for (m, j) in self.chain_candidates(chrom, plus, &chain) {
                if self.ends_fit(m, j, chain.len(), x5, x3) { out.push(m); }
            }
        }
        out.sort_unstable();
        out.dedup();
    }
}

/// Reads shared among reported models: `base[m]` reads counted for model m directly, and classes of reads
/// that fit several models (a sorted model list and its read count). Returns base + EM share per model.
pub fn share_by_em(base: &[f64], classes: &FxHashMap<Vec<u32>, f64>, max_iter: usize) -> Vec<f64> {
    let mut extra = vec![0.0f64; base.len()];
    let mut multi: Vec<(&Vec<u32>, f64)> = Vec::new();
    for (c, &n) in classes {
        if c.len() == 1 { extra[c[0] as usize] += n; } else { multi.push((c, n)); }
    }
    // Classes in the order of their model lists.
    multi.sort_unstable_by(|a, b| a.0.cmp(b.0));
    let fixed = extra.clone();
    let mut theta: Vec<f64> = base.iter().zip(&fixed).map(|(b, f)| b + f).collect();
    let mut share = vec![0.0f64; base.len()];
    for _ in 0..max_iter {
        let mut new = vec![0.0f64; base.len()];
        for &(c, n) in &multi {
            let tot: f64 = c.iter().map(|&m| theta[m as usize]).sum();
            if tot > 0.0 {
                for &m in c { new[m as usize] += n * theta[m as usize] / tot; }
            } else {
                for &m in c { new[m as usize] += n / c.len() as f64; }
            }
        }
        let diff: f64 = new.iter().zip(&share).map(|(a, b)| (a - b).abs()).sum();
        share = new;
        for i in 0..base.len() { theta[i] = base[i] + fixed[i] + share[i]; }
        if diff < 1e-3 { break; }
    }
    theta
}

/// What `reassigned_counts` did, for the log.
#[derive(Debug, Default)]
pub struct QuantStats {
    pub rejected_reads: usize,
    pub rejected_fitted: usize,
    pub unsupported_reads: usize,
    pub unsupported_fitted: usize,
    pub quant_reads: usize,
    pub quant_fitted: usize,
    pub ambiguous: usize,
}

/// Removed one-edit siblings and their parent models, from the filter's `<filter csv>.siblings.tsv` (empty when
/// the file is absent).
pub fn load_sibling_parents(path: &str) -> FxHashMap<i64, i64> {
    let Ok(text) = std::fs::read_to_string(path) else { return FxHashMap::default() };
    text.lines().skip(1).filter_map(|l| {
        let mut f = l.split('\t');
        Some((f.next()?.parse().ok()?, f.next()?.parse().ok()?))
    }).collect()
}

/// Counts per reported model and dataset: `direct` (reads counted for their own reported model) plus the reads of
/// unreported models and the quantification-only reads, each shared among the reported models it fits. `may_share`
/// is the read gate for sharing, `quant_counts` the one for quantification-only reads.
pub fn reassigned_counts(
    db: &PostDB,
    reported: &[i64],
    direct: &FxHashMap<i64, FxHashMap<String, f64>>,
    sibling_parent: &FxHashMap<i64, i64>,
    counted_directly: impl Fn(&crate::post_processing::ObservedRead) -> bool,
    may_share: impl Fn(&crate::post_processing::ObservedRead) -> bool,
    quant_reads: &[QuantRead],
    quant_counts: impl Fn(&QuantRead) -> bool,
    opts: QuantOptions,
) -> (FxHashMap<i64, FxHashMap<String, f64>>, QuantStats) {
    let index = ModelIndex::build(db, reported, opts);
    let n = index.ids.len();
    let slot: FxHashMap<i64, u32> = index.ids.iter().enumerate().map(|(i, &t)| (t, i as u32)).collect();
    let ds_of: FxHashMap<&str, usize> = db.datasets.iter().enumerate().map(|(i, d)| (d.as_str(), i)).collect();
    let nd = db.datasets.len();
    let mut base = vec![vec![0.0f64; n]; nd];
    for (tid, by_ds) in direct {
        let Some(&m) = slot.get(tid) else { continue };
        for (d, &c) in by_ds {
            if let Some(&di) = ds_of.get(d.as_str()) { base[di][m as usize] += c; }
        }
    }
    let mut classes: Vec<FxHashMap<Vec<u32>, f64>> = vec![FxHashMap::default(); nd];
    let mut stats = QuantStats::default();

    // Reads of unreported models: the model's chain is matched once, the ends per read.
    let edge_ends: FxHashMap<i64, (i64, i64)> = db.edges.iter().map(|e| (e.edge_id, (e.v1, e.v2))).collect();
    let edge_strand: FxHashMap<i64, bool> = db.edges.iter().map(|e| (e.edge_id, e.strand == "+")).collect();
    struct Chain { chrom: Arc<str>, plus: bool, k: usize, cands: Vec<(u32, u32)> }
    let mut chains: FxHashMap<i64, Option<Chain>> = FxHashMap::default();
    let mut buf: Vec<u32> = Vec::new();
    for r in &db.observed {
        let on_reported = slot.contains_key(&r.transcript_id);
        // Without partial reads: a read of a reported model counts only for it, and only with end support
        // (as `direct` has it); a read of a removed model is shared only with end support too.
        if on_reported && (counted_directly(r) || !opts.partial_reads) { continue; }
        let Some(&di) = ds_of.get(&*r.dataset) else { continue };
        if !may_share(r) { continue; }
        if !opts.partial_reads && !counted_directly(r) { continue; }
        if on_reported { stats.unsupported_reads += 1; } else { stats.rejected_reads += 1; }
        // A read of a one-edit sibling removed by the filter counts for the sibling's parent model.
        if let Some(m) = sibling_parent.get(&r.transcript_id).and_then(|p| slot.get(p)) {
            stats.rejected_fitted += 1;
            *classes[di].entry(vec![*m]).or_insert(0.0) += 1.0;
            continue;
        }
        let entry = chains.entry(r.transcript_id).or_insert_with(|| {
            let &ti = db.tx_map.get(&r.transcript_id)?;
            let t = &db.transcripts[ti];
            let &plus = edge_strand.get(&t.start_exon)?;
            let pos = model_positions(db, t, &edge_ends)?;
            let chrom: Arc<str> = Arc::from(t.chromosome.as_str());
            let chain = oriented_introns(&pos, plus);
            let cands = if chain.is_empty() { Vec::new() } else { index.chain_candidates(&chrom, plus, &chain) };
            Some(Chain { chrom, plus, k: chain.len(), cands })
        });
        let Some(c) = entry.as_ref() else { continue };
        let s = if c.plus { 1 } else { -1 };
        let (x5, x3) = { let (a, b) = (s * r.start_pos, s * r.end_pos); (a.min(b), a.max(b)) };
        buf.clear();
        if c.k == 0 {
            index.exon_candidates(&c.chrom, c.plus, x5, x3, &mut buf);
        } else {
            buf.extend(c.cands.iter().filter(|&&(m, j)| index.ends_fit(m, j, c.k, x5, x3)).map(|&(m, _)| m));
        }
        buf.sort_unstable();
        buf.dedup();
        if buf.is_empty() { continue; }
        if on_reported { stats.unsupported_fitted += 1; } else { stats.rejected_fitted += 1; }
        if buf.len() > 1 { stats.ambiguous += 1; }
        *classes[di].entry(buf.clone()).or_insert(0.0) += 1.0;
    }

    // Quantification-only reads, matched on their own parsed chains.
    for q in quant_reads {
        let Some(&di) = ds_of.get(&*q.dataset) else { continue };
        if !quant_counts(q) { continue; }
        stats.quant_reads += 1;
        index.read_candidates(&q.chromosome, &*q.strand == "+", &q.positions, &mut buf);
        if buf.is_empty() { continue; }
        stats.quant_fitted += 1;
        if buf.len() > 1 { stats.ambiguous += 1; }
        *classes[di].entry(buf.clone()).or_insert(0.0) += 1.0;
    }

    let mut out: FxHashMap<i64, FxHashMap<String, f64>> = FxHashMap::default();
    for di in 0..nd {
        let total = share_by_em(&base[di], &classes[di], opts.max_iter);
        for (m, &v) in total.iter().enumerate() {
            if v > 0.0 {
                out.entry(index.ids[m]).or_default().insert(db.datasets[di].clone(), v);
            }
        }
    }
    (out, stats)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn em_shares_ambiguous_reads_by_abundance_and_conserves_them() {
        // model 0 has 90 reads of its own, model 1 has 10; 50 reads fit both, 5 fit only model 1
        let base = vec![90.0, 10.0];
        let mut classes = FxHashMap::default();
        classes.insert(vec![0u32, 1], 50.0);
        classes.insert(vec![1u32], 5.0);
        let n = share_by_em(&base, &classes, 500);
        assert!((n.iter().sum::<f64>() - 155.0).abs() < 1e-6);
        // the fixed point: theta0 / theta1 = 90 / 15 for the shared reads
        assert!((n[0] - (90.0 + 50.0 * 90.0 / 105.0)).abs() < 0.05, "{:?}", n);
    }

    #[test]
    fn introns_are_oriented_five_to_three_on_both_strands() {
        // + strand: exons [100,200], [300,400]; - strand: the same exons read from 400 down to 100
        assert_eq!(oriented_introns(&[100, 200, 300, 400], true), vec![(200, 300)]);
        assert_eq!(oriented_introns(&[400, 300, 200, 100], false), vec![(-300, -200)]);
    }
}
