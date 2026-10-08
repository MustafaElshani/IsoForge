// Learned admission on the read hypergraph.
//
// - `features`: for every candidate model (annotated transcripts and novel paths with reads), its structure, the
//   rules' read counts and decision, the reads carrying its introns, the read classes that fit it, its relation to the
//   nearest reported model of its gene, the ends of its reads, its gene and the candidates sharing its intron chain.
// - `AdmissionModel`: a tree ensemble read from a JSON file, which gives each candidate a probability.
// - `select`: reported candidates below `drop` leave the report; unreported candidates at or above `add` join it,
//   one per intron chain, when they pass the guards. A class of candidates may have thresholds of its own (`by_class`).

use std::fs;
use std::io::Write;

use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use serde::Deserialize;

use crate::errors::{IsoForgeError, IsoForgeResult};
use crate::hyper_core::{EndModel, Fit, HyperCore, ReadChains};
use crate::post_processing::PostDB;

const EM_ROUNDS: usize = 100;

/// Names of the numeric features, in the order of `CandFeatures::num`.
pub const NUM: [&str; 43] = [
    "n_exons", "length", "reads", "frac_a_eligible", "passes_read_gates", "utr_or_intergenic_tss",
    "phased_end_reads", "exact_junction_reads", "j_min_support", "j_min_novel_support", "n_novel_junctions",
    "w_fit", "w_exact", "w_full", "w_unique", "fanout", "em_reads", "em_frac", "full_share", "exact_share",
    "unique_share", "rel_share", "rel_reads", "rel_diff", "n_containers", "tss_cl", "tes_cl", "supported", "ip",
    "tailed", "shifted", "read_len", "five_iqr", "three_iqr", "gene_reads", "gene_max", "gene_cands", "gene_share",
    "gene_rank", "chain_n", "chain_rank", "chain_share", "chain_reported",
];

#[derive(Deserialize)]
struct TreeJson { feature: Vec<i64>, threshold: Vec<f64>, left: Vec<i64>, right: Vec<i64>, missing_left: Vec<i64>, leaf: Vec<i64>, value: Vec<f64> }

#[derive(Deserialize)]
struct ModelJson { format: String, features: Vec<String>, baseline: f64, trees: Vec<TreeJson> }

struct Tree { feature: Vec<u32>, threshold: Vec<f64>, left: Vec<u32>, right: Vec<u32>, missing_left: Vec<bool>, leaf: Vec<bool>, value: Vec<f64> }

/// A gradient-boosted tree ensemble with numeric splits (categories one-hot), read from a JSON file.
pub struct AdmissionModel {
    pub features: Vec<String>,
    baseline: f64,
    trees: Vec<Tree>,
}

impl AdmissionModel {
    pub fn load(path: &str) -> IsoForgeResult<Self> {
        let text = fs::read_to_string(path).map_err(|e| IsoForgeError::ConfigError(format!("admission model {path}: {e}")))?;
        let m: ModelJson = serde_json::from_str(&text).map_err(|e| IsoForgeError::ConfigError(format!("admission model {path}: {e}")))?;
        if m.format != "isoforge-admission-trees-1" {
            return Err(IsoForgeError::ConfigError(format!("admission model {path}: unknown format {}", m.format)));
        }
        let trees = m.trees.into_iter().map(|t| Tree {
            feature: t.feature.iter().map(|&x| x.max(0) as u32).collect(),
            threshold: t.threshold,
            left: t.left.iter().map(|&x| x.max(0) as u32).collect(),
            right: t.right.iter().map(|&x| x.max(0) as u32).collect(),
            missing_left: t.missing_left.iter().map(|&x| x != 0).collect(),
            leaf: t.leaf.iter().map(|&x| x != 0).collect(),
            value: t.value,
        }).collect();
        Ok(Self { features: m.features, baseline: m.baseline, trees })
    }

    /// Probability for one feature vector (in `self.features` order); NaN is a missing value.
    pub fn predict(&self, x: &[f64]) -> f64 {
        let mut raw = self.baseline;
        for t in &self.trees {
            let mut n = 0usize;
            loop {
                if t.leaf[n] { raw += t.value[n]; break; }
                let v = x[t.feature[n] as usize];
                let left = if v.is_nan() { t.missing_left[n] } else { v <= t.threshold[n] };
                n = if left { t.left[n] } else { t.right[n] } as usize;
            }
        }
        1.0 / (1.0 + (-raw).exp())
    }
}

/// What the rules' ledger holds for one candidate (the columns of ISOFORGE_FILTER_TRACE's decisions table).
#[derive(Clone)]
pub struct LedgerRow {
    pub reads: u64,
    pub frac_a_eligible: u64,
    pub passes_read_gates: u64,
    pub utr_or_intergenic_tss: u64,
    pub phased_end_reads: u64,
    pub exact_junction_reads: u64,
    pub decision: &'static str,
    pub novelty: String,
}

/// One candidate's features: numeric (NUM order) and the three categorical ones.
pub struct CandFeatures {
    pub tid: i64,
    pub num: [f64; 43],
    pub novelty: String,
    pub decision: &'static str,
    pub rel_kind: &'static str,
    pub chain_key: Vec<(i64, i64)>,
    pub chrom: String,
    pub plus: bool,
    pub reported: bool,
    /// exons 5'->3' as (5' boundary, 3' boundary); empty unless kept for a dump
    pub exons: Vec<(i64, i64)>,
}

impl CandFeatures {
    pub fn get(&self, name: &str) -> f64 {
        if let Some(i) = NUM.iter().position(|&n| n == name) { return self.num[i]; }
        if let Some(v) = name.strip_prefix("novelty=") { return (self.novelty == v) as u8 as f64; }
        if let Some(v) = name.strip_prefix("decision=") { return (self.decision == v) as u8 as f64; }
        if let Some(v) = name.strip_prefix("rel_kind=") { return (self.rel_kind == v) as u8 as f64; }
        f64::NAN
    }
}

fn quantile(sorted: &[f64], q: f64) -> f64 {
    if sorted.is_empty() { return f64::NAN; }
    let h = (sorted.len() - 1) as f64 * q;
    let lo = h.floor() as usize;
    let hi = (lo + 1).min(sorted.len() - 1);
    sorted[lo] + (h - lo as f64) * (sorted[hi] - sorted[lo])
}

fn median(sorted: &[f64]) -> f64 { quantile(sorted, 0.5) }

/// Introns two sorted, duplicate-free chains share.
fn sorted_common(a: &[(i64, i64)], b: &[(i64, i64)]) -> usize {
    let (mut i, mut j, mut n) = (0, 0, 0);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => { n += 1; i += 1; j += 1; }
        }
    }
    n
}

/// The relation of chain `a` to the reported chain `b`.
fn edit_kind(a: &[(i64, i64)], b: &[(i64, i64)], plus: bool) -> &'static str {
    if a.len() < b.len() {
        for i in 0..=(b.len() - a.len()) {
            if &b[i..i + a.len()] == a {
                return if plus {
                    if i + a.len() == b.len() { "sub5" } else if i == 0 { "sub3" } else { "submid" }
                } else if i + a.len() == b.len() { "sub3" } else if i == 0 { "sub5" } else { "submid" };
            }
        }
    }
    let sa: FxHashSet<(i64, i64)> = a.iter().copied().collect();
    let sb: FxHashSet<(i64, i64)> = b.iter().copied().collect();
    let only_a: Vec<(i64, i64)> = sa.difference(&sb).copied().collect();
    let only_b: Vec<(i64, i64)> = sb.difference(&sa).copied().collect();
    if only_a.len() == 1 && only_b.len() == 1 {
        let ((x1, y1), (x2, y2)) = (only_a[0], only_b[0]);
        if x1 == x2 || y1 == y2 {
            return if (x1 - x2).abs() + (y1 - y2).abs() <= 30 { "shift" } else { "altsite" };
        }
        return "swap";
    }
    if only_a.len() == 1 && only_b.len() == 2 { return "skip"; }
    if only_a.len() == 2 && only_b.len() == 1 { return "extra"; }
    if only_a.is_empty() && !only_b.is_empty() { return "subset"; }
    if only_b.is_empty() && !only_a.is_empty() { return "superset"; }
    "other"
}

/// Features of every candidate of `core` that the ledger lists, sorted by transcript ID.
pub fn features(db: &PostDB, core: &HyperCore, chains: &ReadChains, ledger: &FxHashMap<i64, LedgerRow>,
                reported: &FxHashSet<i64>, with_exons: bool) -> Vec<CandFeatures> {
    // candidates in transcript-ID order (the ledger's order, which ties follow)
    let mut order: Vec<usize> = (0..core.cands.len()).filter(|&i| ledger.contains_key(&core.cands[i].tid)).collect();
    order.sort_by_key(|&i| core.cands[i].tid);
    let slot: FxHashMap<i64, usize> = core.cands.iter().enumerate().map(|(i, c)| (c.tid, i)).collect();
    // genomic introns (lo, hi) of each candidate
    let key_of = |a: i64, b: i64| (a.min(b), a.max(b));
    let cand_junctions = |ci: usize| -> Vec<(i64, i64)> {
        core.cands[ci].introns.iter().map(|&(d, a)| key_of(d, a)).collect()
    };
    // guide introns: those of annotated transcripts
    let edge_ends: FxHashMap<i64, (i64, i64)> = db.edges.iter().map(|e| (e.edge_id, (e.v1, e.v2))).collect();
    let edge_strand: FxHashMap<i64, bool> = db.edges.iter().map(|e| (e.edge_id, e.strand == "+")).collect();
    let mut guide: FxHashSet<(String, bool, i64, i64)> = FxHashSet::default();
    for t in &db.transcripts {
        if !db.annot.known_tx_ids.contains(&t.transcript_id) { continue; }
        let Some(&plus) = edge_strand.get(&t.start_exon) else { continue };
        let Some(pos) = crate::quant_em::model_positions(db, t, &edge_ends) else { continue };
        let ex: Vec<(i64, i64)> = pos.chunks(2).map(|p| (p[0], p[1])).collect();
        for w in ex.windows(2) {
            let (lo, hi) = key_of(w[0].1, w[1].0);
            guide.insert((t.chromosome.clone(), plus, lo, hi));
        }
    }
    let t_feat = std::time::Instant::now();
    // reads carrying each intron, over every read of the run
    let mut support: FxHashMap<(&str, bool, i64, i64), f64> = FxHashMap::default();
    for (i, r) in db.observed.iter().enumerate() {
        let c = chains.chain_of.get(i).copied().unwrap_or(0) as usize;
        if c == 0 { continue; }
        let plus = &*r.strand == "+";
        for &(d, a) in chains.chains[c].iter() {
            let (lo, hi) = key_of(d, a);
            *support.entry((&*r.chromosome, plus, lo, hi)).or_insert(0.0) += 1.0;
        }
    }
    log::debug!("[Post] Admission features: guide and junction support {:.1}s", t_feat.elapsed().as_secs_f64());
    // read-hypergraph fits: weight fitting each candidate, exactly, full-length, uniquely; fan-out
    let nc = core.cands.len();
    let (mut w_fit, mut w_exact, mut w_full, mut w_unique, mut fan) = (vec![0.0; nc], vec![0.0; nc], vec![0.0; nc], vec![0.0; nc], vec![0.0; nc]);
    for (k, cl) in core.classes.iter().enumerate() {
        let fs = core.fits_of(k);
        let n = fs.len() as f64;
        for f in fs {
            let m = f.cand as usize;
            w_fit[m] += cl.weight;
            fan[m] += cl.weight * n;
            if f.shift == 0 { w_exact[m] += cl.weight; }
            if core.cands[m].whole(f) { w_full[m] += cl.weight; }
            if fs.len() == 1 { w_unique[m] += cl.weight; }
        }
    }
    // the counting's likelihood, measured on fits to reported models, then an EM over every candidate
    let rep_cand: Vec<bool> = core.cands.iter().map(|c| reported.contains(&c.tid)).collect();
    let ends = EndModel::measure(core, chains, |c| rep_cand[c]);
    let lik = |cl_tail: bool, f: &Fit| ends.lik(core, cl_tail, f);
    // only classes with fits and candidates with fits take part
    let mut theta = vec![1.0f64; nc];
    let has_fit: Vec<bool> = { let mut v = vec![false; nc]; for f in &core.fits { v[f.cand as usize] = true; } v };
    // The EM splits into sets of candidates linked by the classes that fit them; each set runs on its own (in
    // parallel), with its classes in class order.
    let mut parent: Vec<u32> = (0..nc as u32).collect();
    fn find(p: &mut [u32], mut x: u32) -> u32 {
        while p[x as usize] != x { p[x as usize] = p[p[x as usize] as usize]; x = p[x as usize]; }
        x
    }
    for k in 0..core.classes.len() {
        let fs = core.fits_of(k);
        if fs.len() < 2 { continue; }
        let a = find(&mut parent, fs[0].cand);
        for f in &fs[1..] {
            let b = find(&mut parent, f.cand);
            if a != b { parent[b as usize] = a; }
        }
    }
    let mut set_of_root: FxHashMap<u32, usize> = FxHashMap::default();
    let mut sets: Vec<Vec<u32>> = Vec::new();
    for k in 0..core.classes.len() {
        let fs = core.fits_of(k);
        if fs.is_empty() { continue; }
        let r = find(&mut parent, fs[0].cand);
        let s = match set_of_root.get(&r) { Some(&s) => s, None => { sets.push(Vec::new()); set_of_root.insert(r, sets.len() - 1); sets.len() - 1 } };
        sets[s].push(k as u32);
    }
    let mut local = vec![u32::MAX; nc];
    let mut set_cands: Vec<Vec<u32>> = vec![Vec::new(); sets.len()];
    for (s, ks) in sets.iter().enumerate() {
        for &k in ks {
            for f in core.fits_of(k as usize) {
                if local[f.cand as usize] == u32::MAX { local[f.cand as usize] = set_cands[s].len() as u32; set_cands[s].push(f.cand); }
            }
        }
    }
    let solved: Vec<Vec<f64>> = sets.par_iter().enumerate().map(|(s, ks)| {
        let n = set_cands[s].len();
        // the set's fits class by class, each class's in the order of the candidates' transcript IDs, with each fit's
        // set-local candidate and likelihood
        let mut bounds: Vec<u32> = Vec::with_capacity(ks.len() + 1);
        let (mut slot, mut lv): (Vec<u32>, Vec<f64>) = (Vec::new(), Vec::new());
        bounds.push(0);
        for &k in ks {
            let k = k as usize;
            let mut idx: Vec<u32> = (core.offsets[k]..core.offsets[k + 1]).collect();
            idx.sort_by_key(|&q| core.cands[core.fits[q as usize].cand as usize].tid);
            let tail = core.classes[k].tail;
            for q in idx {
                let f = &core.fits[q as usize];
                slot.push(local[f.cand as usize]);
                lv.push(lik(tail, f));
            }
            bounds.push(slot.len() as u32);
        }
        let (mut th, mut acc) = (vec![1.0f64; n], vec![0.0f64; n]);
        for _ in 0..EM_ROUNDS {
            acc.iter_mut().for_each(|a| *a = 0.0);
            for (j, &k) in ks.iter().enumerate() {
                let (lo, hi) = (bounds[j] as usize, bounds[j + 1] as usize);
                let mut sum = 0.0;
                for q in lo..hi { sum += lv[q] * th[slot[q] as usize]; }
                if sum <= 0.0 { sum = 1e-300; }
                let qv = core.classes[k as usize].weight / sum;
                for q in lo..hi { acc[slot[q] as usize] += lv[q] * qv; }
            }
            for i in 0..n { th[i] = (th[i] * acc[i]).max(1e-12); }
        }
        th
    }).collect();
    for (s, th) in solved.into_iter().enumerate() {
        for (i, &m) in set_cands[s].iter().enumerate() { theta[m as usize] = th[i]; }
    }
    log::debug!("[Post] Admission EM: {} sets of candidates, the largest {} classes, done at {:.1}s",
                sets.len(), sets.iter().map(|v| v.len()).max().unwrap_or(0), t_feat.elapsed().as_secs_f64());
    // gene and chain groups over the ledger's candidates
    let reads_of = |ci: usize| ledger[&core.cands[ci].tid].reads as f64;
    let mut gene_sum: FxHashMap<i64, f64> = FxHashMap::default();
    let mut gene_max: FxHashMap<i64, f64> = FxHashMap::default();
    let mut gene_n: FxHashMap<i64, f64> = FxHashMap::default();
    let mut gene_reads_list: FxHashMap<i64, Vec<f64>> = FxHashMap::default();
    type ChainKey = (String, bool, Vec<(i64, i64)>);
    let chain_key_of = |ci: usize| -> Option<ChainKey> {
        let c = &core.cands[ci];
        if c.exons.len() < 2 { return None; }
        let mut js = cand_junctions(ci);
        js.sort_unstable();
        Some((c.chrom.to_string(), c.plus, js))
    };
    let mut chain_sum: FxHashMap<ChainKey, f64> = FxHashMap::default();
    let mut chain_n: FxHashMap<ChainKey, f64> = FxHashMap::default();
    let mut chain_rep: FxHashMap<ChainKey, bool> = FxHashMap::default();
    let mut chain_members: FxHashMap<ChainKey, Vec<(f64, i64)>> = FxHashMap::default();
    for &ci in &order {
        let c = &core.cands[ci];
        let r = reads_of(ci);
        *gene_sum.entry(c.gene).or_insert(0.0) += r;
        let gm = gene_max.entry(c.gene).or_insert(f64::MIN);
        *gm = gm.max(r);
        *gene_n.entry(c.gene).or_insert(0.0) += 1.0;
        gene_reads_list.entry(c.gene).or_default().push(r);
        if let Some(k) = chain_key_of(ci) {
            *chain_sum.entry(k.clone()).or_insert(0.0) += r;
            *chain_n.entry(k.clone()).or_insert(0.0) += 1.0;
            *chain_rep.entry(k.clone()).or_insert(false) |= rep_cand[ci];
            chain_members.entry(k).or_default().push((r, c.tid));
        }
    }
    // chain rank: by reads, descending; ties in transcript-ID order
    let mut chain_rank: FxHashMap<i64, f64> = FxHashMap::default();
    for v in chain_members.values_mut() {
        v.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap().then(a.1.cmp(&b.1)));
        for (i, &(_, t)) in v.iter().enumerate() { chain_rank.insert(t, (i + 1) as f64); }
    }
    // the reported multi-exon models of each gene, in transcript-ID order
    let mut rep_by_gene: FxHashMap<i64, Vec<(i64, Vec<(i64, i64)>, f64)>> = FxHashMap::default();
    let mut rep_sorted: Vec<i64> = reported.iter().copied().collect();
    rep_sorted.sort_unstable();
    for t in rep_sorted {
        let Some(&ci) = slot.get(&t) else { continue };
        if !ledger.contains_key(&t) || core.cands[ci].exons.len() < 2 { continue; }
        let mut ch = cand_junctions(ci);
        ch.sort_unstable();
        rep_by_gene.entry(core.cands[ci].gene).or_default().push((t, ch, reads_of(ci)));
    }
    // reads of each candidate (all of them, every dataset) for the end features
    let mut by_tid: FxHashMap<i64, Vec<u32>> = FxHashMap::default();
    for (i, r) in db.observed.iter().enumerate() {
        if ledger.contains_key(&r.transcript_id) { by_tid.entry(r.transcript_id).or_default().push(i as u32); }
    }
    let t_rows = std::time::Instant::now();
    let out: Vec<CandFeatures> = order.par_iter().map(|&ci| {
        let c = &core.cands[ci];
        let lg = &ledger[&c.tid];
        let mut x = [f64::NAN; 43];
        let set = |x: &mut [f64; 43], name: &str, v: f64| { x[NUM.iter().position(|&n| n == name).unwrap()] = v; };
        set(&mut x, "n_exons", c.exons.len() as f64);
        set(&mut x, "length", c.len as f64);
        set(&mut x, "reads", lg.reads as f64);
        set(&mut x, "frac_a_eligible", lg.frac_a_eligible as f64);
        set(&mut x, "passes_read_gates", lg.passes_read_gates as f64);
        set(&mut x, "utr_or_intergenic_tss", lg.utr_or_intergenic_tss as f64);
        set(&mut x, "phased_end_reads", lg.phased_end_reads as f64);
        set(&mut x, "exact_junction_reads", lg.exact_junction_reads as f64);
        let js = cand_junctions(ci);
        if js.is_empty() {
            set(&mut x, "n_novel_junctions", 0.0);
        } else {
            let sup: Vec<f64> = js.iter().map(|&(lo, hi)| support.get(&(&*c.chrom, c.plus, lo, hi)).copied().unwrap_or(0.0)).collect();
            set(&mut x, "j_min_support", sup.iter().copied().fold(f64::INFINITY, f64::min));
            let nov: Vec<f64> = js.iter().zip(&sup).filter(|((lo, hi), _)| !guide.contains(&(c.chrom.to_string(), c.plus, *lo, *hi)))
                .map(|(_, &s)| s).collect();
            if !nov.is_empty() { set(&mut x, "j_min_novel_support", nov.iter().copied().fold(f64::INFINITY, f64::min)); }
            set(&mut x, "n_novel_junctions", nov.len() as f64);
        }
        set(&mut x, "w_fit", w_fit[ci]);
        set(&mut x, "w_exact", w_exact[ci]);
        set(&mut x, "w_full", w_full[ci]);
        set(&mut x, "w_unique", w_unique[ci]);
        let has = w_fit[ci] > 0.0;
        set(&mut x, "fanout", if has { fan[ci] / w_fit[ci] } else { f64::NAN });
        let em = if has_fit[ci] { theta[ci] } else { 0.0 };
        set(&mut x, "em_reads", em);
        let denom = w_fit[ci].max(1e-9);
        set(&mut x, "em_frac", if has { em / denom } else { f64::NAN });
        set(&mut x, "full_share", w_full[ci] / denom);
        set(&mut x, "exact_share", w_exact[ci] / denom);
        set(&mut x, "unique_share", w_unique[ci] / denom);
        // relation to the nearest reported model of the gene
        let mut rel_kind = "single";
        let mut n_cont = 0.0;
        if c.exons.len() >= 2 {
            let mut ch = js.clone();
            ch.sort_unstable();
            let mut best: Option<((u8, usize, f64), i64, &Vec<(i64, i64)>, f64)> = None;
            if let Some(list) = rep_by_gene.get(&c.gene) {
                for (p, pc, pr) in list {
                    if *p == c.tid { continue; }
                    let sub = pc.len() > ch.len() && (0..=(pc.len() - ch.len())).any(|i| pc[i..i + ch.len()] == ch[..]);
                    if sub { n_cont += 1.0; }
                    let diff = ch.len() + pc.len() - 2 * sorted_common(&ch, pc);
                    let key = (if sub { 0u8 } else { 1u8 }, diff, -pr);
                    let better = match &best {
                        None => true,
                        Some((bk, _, _, _)) => key.0 < bk.0 || (key.0 == bk.0 && (key.1 < bk.1 || (key.1 == bk.1 && key.2 < bk.2))),
                    };
                    if better { best = Some((key, *p, pc, *pr)); }
                }
            }
            match best {
                None => { rel_kind = "none"; set(&mut x, "rel_reads", 0.0); }
                Some((key, _, pc, pr)) => {
                    // edit kind of the two chains, given in genomic order; 5' and 3' follow the strand
                    rel_kind = edit_kind(&ch, pc, c.plus);
                    let r = lg.reads as f64;
                    set(&mut x, "rel_share", if r + pr > 0.0 { r / (r + pr) } else { f64::NAN });
                    set(&mut x, "rel_reads", pr);
                    set(&mut x, "rel_diff", key.1 as f64);
                }
            }
        }
        set(&mut x, "n_containers", n_cont);
        // ends of its reads
        if let Some(rows) = by_tid.get(&c.tid) {
            let n = rows.len() as f64;
            let (mut tss, mut tes, mut sup, mut ip, mut tail, mut shf) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
            let (mut five, mut three, mut rl) = (Vec::with_capacity(rows.len()), Vec::with_capacity(rows.len()), Vec::with_capacity(rows.len()));
            for &i in rows {
                let r = &db.observed[i as usize];
                tss += r.tss_site.is_some() as u8 as f64;
                tes += r.polya_site.is_some() as u8 as f64;
                sup += (r.cluster_supported == Some(true)) as u8 as f64;
                ip += (r.is_internal_priming == Some(true)) as u8 as f64;
                tail += (r.polya_length.unwrap_or(0) > 10 || r.polya_softclip.unwrap_or(0) > 10) as u8 as f64;
                shf += (r.max_splice_shift.unwrap_or(0) > 0) as u8 as f64;
                let plus = &*r.strand == "+";
                five.push(if plus { r.start_pos as f64 } else { -(r.end_pos as f64) });
                three.push(if plus { r.end_pos as f64 } else { -(r.start_pos as f64) });
                rl.push(r.read_length as f64);
            }
            for v in [&mut five, &mut three, &mut rl] { v.sort_by(|a, b| a.partial_cmp(b).unwrap()); }
            set(&mut x, "tss_cl", tss / n);
            set(&mut x, "tes_cl", tes / n);
            set(&mut x, "supported", sup / n);
            set(&mut x, "ip", ip / n);
            set(&mut x, "tailed", tail / n);
            set(&mut x, "shifted", shf / n);
            set(&mut x, "read_len", median(&rl));
            set(&mut x, "five_iqr", quantile(&five, 0.75) - quantile(&five, 0.25));
            set(&mut x, "three_iqr", quantile(&three, 0.75) - quantile(&three, 0.25));
        }
        // gene context
        let gs = gene_sum[&c.gene];
        let r = lg.reads as f64;
        set(&mut x, "gene_reads", gs);
        set(&mut x, "gene_max", gene_max[&c.gene]);
        set(&mut x, "gene_cands", gene_n[&c.gene]);
        set(&mut x, "gene_share", r / gs.max(1.0));
        set(&mut x, "gene_rank", 1.0 + gene_reads_list[&c.gene].iter().filter(|&&v| v > r).count() as f64);
        // candidates sharing its intron chain
        let key = chain_key_of(ci);
        match &key {
            Some(k) => {
                set(&mut x, "chain_n", chain_n[k]);
                set(&mut x, "chain_rank", chain_rank[&c.tid]);
                set(&mut x, "chain_share", r / chain_sum[k].max(1.0));
                set(&mut x, "chain_reported", (chain_rep[k] && !rep_cand[ci]) as u8 as f64);
            }
            None => {
                set(&mut x, "chain_n", 1.0);
                set(&mut x, "chain_rank", 1.0);
                set(&mut x, "chain_share", r / r.max(1.0));
                set(&mut x, "chain_reported", 0.0);
            }
        }
        CandFeatures {
            tid: c.tid, num: x, novelty: lg.novelty.clone(), decision: lg.decision, rel_kind,
            chain_key: key.map(|k| k.2).unwrap_or_default(), chrom: c.chrom.to_string(), plus: c.plus,
            reported: rep_cand[ci], exons: if with_exons { c.exons.clone() } else { Vec::new() },
        }
    }).collect();
    log::debug!("[Post] Admission features: candidate rows {:.1}s", t_rows.elapsed().as_secs_f64());
    out
}

/// Novelty classes of candidates inside annotated genes.
const GENIC: [&str; 4] = ["Known", "ISM", "NIC", "NNC"];

/// Thresholds and the guards of the selection.
pub struct AdmissionParams {
    pub add: f64,
    pub drop: f64,
    /// Thresholds of single novelty classes, (class, add, drop), in place of `add` and `drop`.
    pub by_class: Vec<(String, f64, f64)>,
    pub guard_tailed: f64,
    pub guard_tes: f64,
    pub guard_ip: f64,
    /// TSS-clustered share a joining candidate outside annotated genes needs (it must also be spliced).
    pub guard_tss: f64,
    /// TSS-clustered share every joining candidate needs (0: none).
    pub guard_tss_all: f64,
    /// Datasets a joining candidate's counted reads must come from (the caller caps it at the run's datasets).
    pub min_datasets: u32,
}

/// The candidates to add to and to drop from the report, with each candidate's probability. `counted`: per candidate
/// with reads that abundance counts, the number of datasets those reads come from. `barred`: candidates that cannot
/// join whatever their probability (the junction rule rejects an intron of theirs).
pub fn select(model: &AdmissionModel, feats: &[CandFeatures], counted: &FxHashMap<i64, u32>, prm: &AdmissionParams,
              barred: &FxHashSet<i64>) -> (FxHashSet<i64>, FxHashSet<i64>, Vec<(i64, f64)>) {
    // each model column resolved once: a numeric feature's slot, or a category value (one-hot)
    enum Col { Num(usize), Novelty(String), Decision(String), RelKind(String), Missing }
    let cols: Vec<Col> = model.features.iter().map(|n| {
        if let Some(i) = NUM.iter().position(|&m| m == n.as_str()) { Col::Num(i) }
        else if let Some(v) = n.strip_prefix("novelty=") { Col::Novelty(v.to_string()) }
        else if let Some(v) = n.strip_prefix("decision=") { Col::Decision(v.to_string()) }
        else if let Some(v) = n.strip_prefix("rel_kind=") { Col::RelKind(v.to_string()) }
        else { Col::Missing }
    }).collect();
    let probs: Vec<f64> = feats.par_iter().map(|f| {
        let x: Vec<f64> = cols.iter().map(|c| match c {
            Col::Num(i) => f.num[*i],
            Col::Novelty(v) => (f.novelty == *v) as u8 as f64,
            Col::Decision(v) => (f.decision == v.as_str()) as u8 as f64,
            Col::RelKind(v) => (f.rel_kind == v.as_str()) as u8 as f64,
            Col::Missing => f64::NAN,
        }).collect();
        model.predict(&x)
    }).collect();
    let mut drop = FxHashSet::default();
    let mut best_new: FxHashMap<(String, bool, Vec<(i64, i64)>, i64), (f64, i64)> = FxHashMap::default();
    let get = |f: &CandFeatures, n: &str| { let v = f.get(n); if v.is_nan() { 0.0 } else { v } };
    for (f, &p) in feats.iter().zip(&probs) {
        let (t_add, t_drop) = prm.by_class.iter().find(|c| c.0 == f.novelty).map_or((prm.add, prm.drop), |c| (c.1, c.2));
        if f.reported {
            if p < t_drop { drop.insert(f.tid); }
            continue;
        }
        if p < t_add || get(f, "chain_reported") > 0.0 || barred.contains(&f.tid) { continue; }
        if counted.get(&f.tid).copied().unwrap_or(0) < prm.min_datasets.max(1) { continue; }
        // its reads' TSS-clustered share must reach `guard_tss_all`; unless it is NIC, their tailed and TES-clustered
        // shares must reach `guard_tailed` and `guard_tes` and their internally primed share stay within `guard_ip`
        if get(f, "tss_cl") < prm.guard_tss_all { continue; }
        if f.novelty != "NIC" && !(get(f, "tailed") >= prm.guard_tailed && get(f, "tes_cl") >= prm.guard_tes && get(f, "ip") <= prm.guard_ip) {
            continue;
        }
        // outside annotated genes (antisense, intergenic, genomic) only spliced candidates with a clustered 5' end join
        if !GENIC.contains(&f.novelty.as_str()) && !(get(f, "n_exons") >= 2.0 && get(f, "tss_cl") >= prm.guard_tss) { continue; }
        // one per intron chain (mono-exon candidates are their own chain)
        let key = (f.chrom.clone(), f.plus, f.chain_key.clone(), if f.chain_key.is_empty() { f.tid } else { 0 });
        let e = best_new.entry(key).or_insert((p, f.tid));
        if p > e.0 { *e = (p, f.tid); }
    }
    let add: FxHashSet<i64> = best_new.values().map(|&(_, t)| t).collect();
    (add, drop, feats.iter().map(|f| f.tid).zip(probs).collect())
}

/// One row per candidate with its features, structure and probability (ISOFORGE_ADMISSION_DUMP). Probabilities are
/// NaN when no model is given; counted is the number of datasets with reads of the candidate that abundance counts.
pub fn dump(path: &str, feats: &[CandFeatures], probs: &[(i64, f64)], counted: &FxHashMap<i64, u32>,
            junction_bad: Option<&FxHashSet<i64>>) -> std::io::Result<()> {
    let mut w = std::io::BufWriter::new(fs::File::create(path)?);
    write!(w, "transcript_id\tchrom\tstrand\texons\tnovelty\tdecision\trel_kind\treported\tcounted")?;
    for n in NUM { write!(w, "\t{n}")?; }
    writeln!(w, "\tp\tjunction_ok")?;
    for (k, f) in feats.iter().enumerate() {
        let ex: Vec<String> = f.exons.iter().map(|&(a, b)| format!("{}-{}", a.min(b), a.max(b))).collect();
        write!(w, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", f.tid, f.chrom, if f.plus { '+' } else { '-' }, ex.join(","),
               f.novelty, f.decision, f.rel_kind, f.reported as u8, counted.get(&f.tid).copied().unwrap_or(0))?;
        for v in f.num { if v.is_nan() { write!(w, "\tnan")? } else { write!(w, "\t{v}")? } }
        let p = probs.get(k).map(|x| x.1).unwrap_or(f64::NAN);
        if p.is_nan() { write!(w, "\tnan")? } else { write!(w, "\t{p}")? }
        // whether the junction rule accepts every intron of the candidate (nan: the rule is off)
        match junction_bad { Some(b) => writeln!(w, "\t{}", !b.contains(&f.tid) as u8)?, None => writeln!(w, "\tnan")? }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn stump() -> AdmissionModel {
        // one split on feature 0 at 0.5 (missing goes left): left leaf -1, right leaf +2; baseline 0.1
        AdmissionModel {
            features: vec!["x".into()],
            baseline: 0.1,
            trees: vec![Tree { feature: vec![0, 0, 0], threshold: vec![0.5, 0.0, 0.0], left: vec![1, 0, 0], right: vec![2, 0, 0],
                               missing_left: vec![true, false, false], leaf: vec![false, true, true], value: vec![0.0, -1.0, 2.0] }],
        }
    }

    #[test]
    fn trees_follow_thresholds_and_missing_values() {
        let m = stump();
        let s = |r: f64| 1.0 / (1.0 + (-r).exp());
        assert!((m.predict(&[0.2]) - s(-0.9)).abs() < 1e-15);
        assert!((m.predict(&[0.5]) - s(-0.9)).abs() < 1e-15, "x <= threshold goes left");
        assert!((m.predict(&[0.7]) - s(2.1)).abs() < 1e-15);
        assert!((m.predict(&[f64::NAN]) - s(-0.9)).abs() < 1e-15, "missing goes left here");
    }

    fn cand(tid: i64, novelty: &str, reported: bool, vals: &[(&str, f64)]) -> CandFeatures {
        let mut num = [f64::NAN; 43];
        for (n, v) in vals { num[NUM.iter().position(|m| m == n).unwrap()] = *v; }
        CandFeatures { tid, num, novelty: novelty.into(), decision: "pass", rel_kind: "none", chain_key: vec![(tid * 100, tid * 100 + 50)],
                       chrom: "chr1".into(), plus: true, reported, exons: Vec::new() }
    }

    #[test]
    fn additions_need_clustered_starts_and_reads_of_two_datasets() {
        // probable (p 0.89) with a read, improbable (0.29) without
        let mut m = stump();
        m.features = vec!["reads".into()];
        let good = [("reads", 3.0), ("tss_cl", 1.0), ("tailed", 1.0), ("tes_cl", 1.0), ("ip", 0.0)];
        let feats = vec![
            cand(1, "NIC", false, &good),
            cand(2, "NIC", false, &[("reads", 3.0), ("tss_cl", 0.2)]),      // its reads start in no TSS cluster
            cand(3, "NIC", false, &good),                                    // reads of one dataset only
            cand(4, "NNC", false, &[("reads", 3.0), ("tss_cl", 1.0), ("tailed", 0.0), ("tes_cl", 1.0), ("ip", 0.0)]),  // no tails
            cand(5, "Known", true, &[("reads", 0.0)]),                       // reported, improbable
            cand(6, "NNC", false, &good),
        ];
        let counted: FxHashMap<i64, u32> = [(1, 2), (2, 2), (3, 1), (4, 3), (6, 2)].into_iter().collect();
        let prm = |tss_all, min_datasets| AdmissionParams { add: 0.7, drop: 0.3, by_class: Vec::new(), guard_tailed: 0.6, guard_tes: 0.6, guard_ip: 0.1,
                                                            guard_tss: 0.8, guard_tss_all: tss_all, min_datasets };
        let ids = |s: &FxHashSet<i64>| { let mut v: Vec<i64> = s.iter().copied().collect(); v.sort_unstable(); v };
        let (add, drop, _) = select(&m, &feats, &counted, &prm(0.6, 2), &FxHashSet::default());
        assert_eq!(ids(&add), vec![1, 6]);
        assert_eq!(ids(&drop), vec![5]);
        assert_eq!(ids(&select(&m, &feats, &counted, &prm(0.6, 1), &FxHashSet::default()).0), vec![1, 3, 6]);
        // neither guard: every probable candidate with counted reads and, unless NIC, tailed 3' ends
        assert_eq!(ids(&select(&m, &feats, &counted, &prm(0.0, 1), &FxHashSet::default()).0), vec![1, 2, 3, 6]);
        // a class with thresholds of its own: NNC candidates (p 0.89) need 0.95 to join, annotated models are never vetoed
        let mut own = prm(0.6, 2);
        own.by_class = vec![("NNC".to_string(), 0.95, 0.5), ("Known".to_string(), 0.7, 0.0)];
        let (add, drop, _) = select(&m, &feats, &counted, &own, &FxHashSet::default());
        assert_eq!(ids(&add), vec![1]);
        assert!(drop.is_empty());
        // a barred candidate does not join
        let barred: FxHashSet<i64> = [1].into_iter().collect();
        assert_eq!(ids(&select(&m, &feats, &counted, &prm(0.6, 2), &barred).0), vec![6]);
        // an infinite addition threshold: no candidate of the class joins, whatever its probability
        let mut never = prm(0.6, 2);
        never.by_class = vec![("NNC".to_string(), f64::INFINITY, 0.7)];
        assert_eq!(ids(&select(&m, &feats, &counted, &never, &FxHashSet::default()).0), vec![1]);
    }

    #[test]
    fn thresholds_read_back_exactly() {
        // floats written as shortest round-trip decimals are read back by serde_json bit for bit
        let mut x: u64 = 0x9E37_79B9_7F4A_7C15;
        let vals: Vec<f64> = (0..20000).map(|_| {
            x ^= x << 13; x ^= x >> 7; x ^= x << 17;
            (x >> 11) as f64 / (1u64 << 53) as f64 * 1e3
        }).collect();
        let text = format!("[{}]", vals.iter().map(|v| format!("{v}")).collect::<Vec<_>>().join(","));
        let back: Vec<f64> = serde_json::from_str(&text).unwrap();
        assert!(vals.iter().zip(&back).all(|(a, b)| a.to_bits() == b.to_bits()));
    }

    #[test]
    fn quantiles_interpolate_linearly() {
        let v = [1.0, 2.0, 4.0, 8.0];
        assert!((quantile(&v, 0.25) - 1.75).abs() < 1e-12);
        assert!((quantile(&v, 0.75) - 5.0).abs() < 1e-12);
        assert!((median(&v) - 3.0).abs() < 1e-12);
        assert_eq!(quantile(&[3.0], 0.75), 3.0);
    }

    #[test]
    fn edit_kinds() {
        let b = vec![(10, 20), (30, 40), (50, 60)];
        assert_eq!(edit_kind(&[(30, 40), (50, 60)], &b, true), "sub5");
        assert_eq!(edit_kind(&[(10, 20), (30, 40)], &b, true), "sub3");
        assert_eq!(edit_kind(&[(10, 20), (30, 40)], &b, false), "sub5");
        assert_eq!(edit_kind(&[(10, 20), (33, 40), (50, 60)], &b, true), "shift");
        assert_eq!(edit_kind(&[(10, 40), (50, 60)], &b, true), "skip");
        assert_eq!(edit_kind(&[(10, 20), (30, 34), (36, 40), (50, 60)], &b, true), "extra");
    }
}
