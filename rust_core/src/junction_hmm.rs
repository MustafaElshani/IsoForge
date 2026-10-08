// Junction posteriors from a pair-HMM.
//
// - A candidate junction's spliced template is the exon bases before its donor joined to those after its acceptor; a
//   read's likelihood under it is the forward sum over all alignments of the read's bases to the template.
// - Gap opening depends on the length of the homopolymer at the template position.
// - A junction group's reads decide together: a mixture over its candidates is fitted by EM from the reads' likelihoods
//   and the candidates' priors, sites the reads do not need are dropped, and each read takes its most probable site.
// - Emission and transition probabilities are estimated by Baum-Welch, with Dirichlet pseudocounts from the defaults.

use serde::{Deserialize, Serialize};

/// Homopolymer lengths the gap-opening probabilities are kept for; longer runs share the last class.
pub const HP_CLASSES: usize = 6;

/// Pair-HMM parameters (probabilities, not logs).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct HmmParams {
    /// A template base read as another base.
    pub mismatch: f64,
    /// Opening an insertion (read bases with no template base) or a deletion (template bases with no read base) after
    /// a match, by the length of the homopolymer at the template position (index 0: no run, 1 base).
    pub open_ins: [f64; HP_CLASSES],
    pub open_del: [f64; HP_CLASSES],
    /// Extending an insertion or a deletion.
    pub extend_ins: f64,
    pub extend_del: f64,
}

impl Default for HmmParams {
    /// Default parameters; the gap-opening probabilities rise with the homopolymer class.
    fn default() -> Self {
        let ramp = |base: f64| { let mut v = [0.0; HP_CLASSES]; for (h, x) in v.iter_mut().enumerate() { *x = (base * (1.0 + 0.6 * h as f64)).min(0.35); } v };
        Self { mismatch: 0.03, open_ins: ramp(0.015), open_del: ramp(0.02), extend_ins: 0.3, extend_del: 0.3 }
    }
}

/// The homopolymer class of each template position: the length of the run of equal bases it lies in, capped.
pub fn homopolymer_classes(t: &[u8]) -> Vec<usize> {
    let mut out = vec![0usize; t.len()];
    let mut i = 0;
    while i < t.len() {
        let j = i + t[i..].iter().take_while(|&&b| b.eq_ignore_ascii_case(&t[i])).count();
        for x in &mut out[i..j] { *x = (j - i - 1).min(HP_CLASSES - 1); }
        i = j;
    }
    out
}

#[inline]
pub(crate) fn logsum(a: f64, b: f64) -> f64 {
    if a == f64::NEG_INFINITY { return b; }
    if b == f64::NEG_INFINITY { return a; }
    let (hi, lo) = if a > b { (a, b) } else { (b, a) };
    hi + (lo - hi).exp().ln_1p()
}

/// Log parameters, precomputed.
struct LogParams {
    same: f64,
    diff: f64,
    background: f64,
    m_to_i: [f64; HP_CLASSES],
    m_to_d: [f64; HP_CLASSES],
    m_to_m: [f64; HP_CLASSES],
    i_to_i: f64,
    i_to_m: f64,
    d_to_d: f64,
    d_to_m: f64,
}

impl LogParams {
    fn from_prob(p: &ProbParams) -> Self {
        let ln = |a: [f64; HP_CLASSES]| a.map(|x| x.ln());
        Self {
            same: p.same.ln(), diff: p.diff.ln(), background: p.background.ln(),
            m_to_i: ln(p.m_to_i), m_to_d: ln(p.m_to_d), m_to_m: ln(p.m_to_m),
            i_to_i: p.i_to_i.ln(), i_to_m: p.i_to_m.ln(), d_to_d: p.d_to_d.ln(), d_to_m: p.d_to_m.ln(),
        }
    }

    fn new(p: &HmmParams) -> Self {
        let mut m_to_i = [0.0; HP_CLASSES];
        let mut m_to_d = [0.0; HP_CLASSES];
        let mut m_to_m = [0.0; HP_CLASSES];
        for h in 0..HP_CLASSES {
            m_to_i[h] = p.open_ins[h].ln();
            m_to_d[h] = p.open_del[h].ln();
            m_to_m[h] = (1.0 - p.open_ins[h] - p.open_del[h]).max(1e-9).ln();
        }
        Self {
            same: (1.0 - p.mismatch).ln(), diff: (p.mismatch / 3.0).ln(), background: 0.25f64.ln(),
            m_to_i, m_to_d, m_to_m,
            i_to_i: p.extend_ins.ln(), i_to_m: (1.0 - p.extend_ins).ln(),
            d_to_d: p.extend_del.ln(), d_to_m: (1.0 - p.extend_del).ln(),
        }
    }
}

/// log P(read | template): the pair-HMM forward algorithm, both sequences aligned end to end.
pub fn log_likelihood(read: &[u8], template: &[u8], params: &HmmParams) -> f64 {
    forward(read, template, &LogParams::new(params), &homopolymer_classes(template)).0
}

/// log P(read | template), the forward sum restricted to cells within `band` of the diagonal from (0, 0) to the end.
pub fn log_likelihood_banded(read: &[u8], template: &[u8], params: &HmmParams, band: usize) -> f64 {
    forward_banded(read, template, &homopolymer_classes(template), &ProbParams::new(params), band)
}

/// The pair-HMM's probabilities, for the banded forward sum, which runs in probability space with each row scaled.
#[derive(Debug, Clone)]
pub struct ProbParams {
    same: f64,
    diff: f64,
    background: f64,
    m_to_i: [f64; HP_CLASSES],
    m_to_d: [f64; HP_CLASSES],
    m_to_m: [f64; HP_CLASSES],
    i_to_i: f64,
    i_to_m: f64,
    d_to_d: f64,
    d_to_m: f64,
}

impl ProbParams {
    pub fn new(p: &HmmParams) -> Self {
        let mut m_to_m = [0.0; HP_CLASSES];
        for h in 0..HP_CLASSES { m_to_m[h] = (1.0 - p.open_ins[h] - p.open_del[h]).max(1e-9); }
        Self {
            same: 1.0 - p.mismatch, diff: p.mismatch / 3.0, background: 0.25,
            m_to_i: p.open_ins, m_to_d: p.open_del, m_to_m,
            i_to_i: p.extend_ins, i_to_m: 1.0 - p.extend_ins,
            d_to_d: p.extend_del, d_to_m: 1.0 - p.extend_del,
        }
    }
}

/// `log_likelihood_banded` with the template's homopolymer classes and the probabilities computed once: two rows per
/// state, each row scaled to sum 1 and the scales' logs added.
pub fn forward_banded(read: &[u8], template: &[u8], hp: &[usize], pp: &ProbParams, band: usize) -> f64 {
    let (n, m) = (read.len(), template.len());
    if n == 0 || m == 0 {
        return forward(read, template, &LogParams::from_prob(pp), hp).0;
    }
    let w = m + 1;
    let (mut pm, mut pi, mut pd) = (vec![0.0f64; w], vec![0.0f64; w], vec![0.0f64; w]);
    let (mut cm, mut ci, mut cd) = (vec![0.0f64; w], vec![0.0f64; w], vec![0.0f64; w]);
    let class = |j: usize| hp[j.min(m - 1)];
    let b = band as i64 + (n as i64 - m as i64).abs();
    let range = |i: usize| {
        let centre = (i * m / n) as i64;
        ((centre - b).max(0) as usize, ((centre + b) as usize).min(m))
    };
    pm[0] = 1.0;
    for j in 1..=range(0).1 {
        let h = class(j - 1);
        pd[j] = pm[j - 1] * pp.m_to_d[h] + pd[j - 1] * pp.d_to_d;
    }
    let mut log_scale = 0.0;
    for i in 1..=n {
        let (j0, j1) = range(i);
        cm.iter_mut().for_each(|x| *x = 0.0);
        ci.iter_mut().for_each(|x| *x = 0.0);
        cd.iter_mut().for_each(|x| *x = 0.0);
        let r = read[i - 1];
        for j in j0..=j1 {
            if j > 0 {
                let e = if r.eq_ignore_ascii_case(&template[j - 1]) { pp.same } else { pp.diff };
                let h = class(j - 1);
                cm[j] = e * (pm[j - 1] * pp.m_to_m[h] + pi[j - 1] * pp.i_to_m + pd[j - 1] * pp.d_to_m);
            }
            let h = class(j);
            ci[j] = pp.background * (pm[j] * pp.m_to_i[h] + pi[j] * pp.i_to_i);
            if j > 0 {
                let h = class(j - 1);
                cd[j] = cm[j - 1] * pp.m_to_d[h] + cd[j - 1] * pp.d_to_d;
            }
        }
        let s: f64 = (j0..=j1).map(|j| cm[j] + ci[j] + cd[j]).sum();
        if !(s > 0.0) { return f64::NEG_INFINITY; }
        for j in j0..=j1 { cm[j] /= s; ci[j] /= s; cd[j] /= s; }
        log_scale += s.ln();
        std::mem::swap(&mut pm, &mut cm);
        std::mem::swap(&mut pi, &mut ci);
        std::mem::swap(&mut pd, &mut cd);
    }
    log_scale + (pm[m] + pi[m] + pd[m]).ln()
}

/// Forward matrices (flattened, (n+1) x (m+1) each, for M, I, D) and the total log-likelihood.
fn forward(read: &[u8], t: &[u8], lp: &LogParams, hp: &[usize]) -> (f64, Vec<f64>, Vec<f64>, Vec<f64>) {
    let (n, m) = (read.len(), t.len());
    let w = m + 1;
    let ninf = f64::NEG_INFINITY;
    let (mut fm, mut fi, mut fd) = (vec![ninf; (n + 1) * w], vec![ninf; (n + 1) * w], vec![ninf; (n + 1) * w]);
    fm[0] = 0.0;
    // the class of the template position a transition out of cell (i, j) looks at: the next template base
    let class = |j: usize| hp[j.min(m.saturating_sub(1))];
    for i in 0..=n {
        for j in 0..=m {
            if i == 0 && j == 0 { continue; }
            let k = i * w + j;
            if i > 0 && j > 0 {
                let e = if read[i - 1].eq_ignore_ascii_case(&t[j - 1]) { lp.same } else { lp.diff };
                let p = (i - 1) * w + (j - 1);
                let h = class(j - 1);
                fm[k] = e + logsum(logsum(fm[p] + lp.m_to_m[h], fi[p] + lp.i_to_m), fd[p] + lp.d_to_m);
            }
            if i > 0 {
                let p = (i - 1) * w + j;
                let h = class(j);
                fi[k] = lp.background + logsum(fm[p] + lp.m_to_i[h], fi[p] + lp.i_to_i);
            }
            if j > 0 {
                let p = i * w + (j - 1);
                let h = class(j - 1);
                fd[k] = logsum(fm[p] + lp.m_to_d[h], fd[p] + lp.d_to_d);
            }
        }
    }
    let end = n * w + m;
    (logsum(logsum(fm[end], fi[end]), fd[end]), fm, fi, fd)
}

/// Backward matrices matching `forward`.
fn backward(read: &[u8], t: &[u8], lp: &LogParams, hp: &[usize]) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let (n, m) = (read.len(), t.len());
    let w = m + 1;
    let ninf = f64::NEG_INFINITY;
    let (mut bm, mut bi, mut bd) = (vec![ninf; (n + 1) * w], vec![ninf; (n + 1) * w], vec![ninf; (n + 1) * w]);
    let end = n * w + m;
    bm[end] = 0.0; bi[end] = 0.0; bd[end] = 0.0;
    let class = |j: usize| hp[j.min(m.saturating_sub(1))];
    for i in (0..=n).rev() {
        for j in (0..=m).rev() {
            let k = i * w + j;
            if k == end { continue; }
            // out of (i, j) into M(i+1, j+1), I(i+1, j), D(i, j+1)
            let em = if i < n && j < m { if read[i].eq_ignore_ascii_case(&t[j]) { lp.same } else { lp.diff } } else { ninf };
            let to_m = if i < n && j < m { em + bm[(i + 1) * w + j + 1] } else { ninf };
            let to_i = if i < n { lp.background + bi[(i + 1) * w + j] } else { ninf };
            let to_d = if j < m { bd[i * w + j + 1] } else { ninf };
            let hm = class(j);
            bm[k] = logsum(logsum(lp.m_to_m[hm] + to_m, lp.m_to_i[class(j)] + to_i), lp.m_to_d[class(j)] + to_d);
            bi[k] = logsum(lp.i_to_m + to_m, lp.i_to_i + to_i);
            bd[k] = logsum(lp.d_to_m + to_m, lp.d_to_d + to_d);
        }
    }
    (bm, bi, bd)
}

/// Expected counts for Baum-Welch, per homopolymer class for the gap opens.
#[derive(Debug, Clone, Default)]
pub struct HmmCounts {
    same: f64,
    diff: f64,
    m_to_i: [f64; HP_CLASSES],
    m_to_d: [f64; HP_CLASSES],
    m_to_m: [f64; HP_CLASSES],
    i_to_i: f64,
    i_to_m: f64,
    d_to_d: f64,
    d_to_m: f64,
    pairs: usize,
}

impl HmmCounts {
    /// Adds one read-template pair's expected transition and emission counts under `params`.
    pub fn add(&mut self, read: &[u8], t: &[u8], params: &HmmParams) {
        if read.is_empty() || t.is_empty() { return; }
        let lp = LogParams::new(params);
        let hp = homopolymer_classes(t);
        let (total, fm, fi, fd) = forward(read, t, &lp, &hp);
        if !total.is_finite() { return; }
        let (bm, bi, bd) = backward(read, t, &lp, &hp);
        let (n, m) = (read.len(), t.len());
        let w = m + 1;
        let class = |j: usize| hp[j.min(m.saturating_sub(1))];
        let post = |x: f64| (x - total).exp();
        for i in 0..=n {
            for j in 0..=m {
                let k = i * w + j;
                if i > 0 && j > 0 {
                    let same = read[i - 1].eq_ignore_ascii_case(&t[j - 1]);
                    let g = post(fm[k] + bm[k]);
                    if same { self.same += g } else { self.diff += g }
                }
                // transitions out of (i, j)
                if i < n && j < m {
                    let em = if read[i].eq_ignore_ascii_case(&t[j]) { lp.same } else { lp.diff };
                    let next = em + bm[(i + 1) * w + j + 1];
                    let h = class(j);
                    self.m_to_m[h] += post(fm[k] + lp.m_to_m[h] + next);
                    self.i_to_m += post(fi[k] + lp.i_to_m + next);
                    self.d_to_m += post(fd[k] + lp.d_to_m + next);
                }
                if i < n {
                    let next = lp.background + bi[(i + 1) * w + j];
                    self.m_to_i[class(j)] += post(fm[k] + lp.m_to_i[class(j)] + next);
                    self.i_to_i += post(fi[k] + lp.i_to_i + next);
                }
                if j < m {
                    let next = bd[i * w + j + 1];
                    self.m_to_d[class(j)] += post(fm[k] + lp.m_to_d[class(j)] + next);
                    self.d_to_d += post(fd[k] + lp.d_to_d + next);
                }
            }
        }
        self.pairs += 1;
    }

    pub fn pairs(&self) -> usize { self.pairs }

    /// New parameters from the counts, each probability shrunk towards `prior` by `strength` pseudo-observations
    /// (a Dirichlet prior with the prior's probabilities as its mean).
    pub fn estimate(&self, prior: &HmmParams, strength: f64) -> HmmParams {
        let beta = |x: f64, total: f64, p: f64| ((x + strength * p) / (total + strength)).clamp(1e-6, 0.5);
        let mut open_ins = [0.0; HP_CLASSES];
        let mut open_del = [0.0; HP_CLASSES];
        for h in 0..HP_CLASSES {
            let out = self.m_to_i[h] + self.m_to_d[h] + self.m_to_m[h];
            open_ins[h] = beta(self.m_to_i[h], out, prior.open_ins[h]);
            open_del[h] = beta(self.m_to_d[h], out, prior.open_del[h]);
        }
        HmmParams {
            mismatch: beta(self.diff, self.same + self.diff, prior.mismatch),
            open_ins, open_del,
            extend_ins: beta(self.i_to_i, self.i_to_i + self.i_to_m, prior.extend_ins).min(0.9),
            extend_del: beta(self.d_to_d, self.d_to_d + self.d_to_m, prior.extend_del).min(0.9),
        }
    }
}

/// Baum-Welch: `rounds` of expectation over `pairs` (read bases, spliced template) and re-estimation, from `start`,
/// shrunk towards the defaults by `strength` pseudo-observations per probability.
pub fn train(pairs: &[(Vec<u8>, Vec<u8>)], start: &HmmParams, rounds: usize, strength: f64) -> HmmParams {
    use rayon::prelude::*;
    let mut p = start.clone();
    for _ in 0..rounds {
        // expected counts per chunk of 64 pairs, summed in chunk order
        let parts: Vec<HmmCounts> = pairs.par_chunks(64).map(|ch| {
            let mut c = HmmCounts::default();
            for (r, t) in ch { c.add(r, t, &p); }
            c
        }).collect();
        let counts = parts.into_iter().fold(HmmCounts::default(), |mut a, b| {
            a.same += b.same; a.diff += b.diff; a.i_to_i += b.i_to_i; a.i_to_m += b.i_to_m; a.d_to_d += b.d_to_d;
            a.d_to_m += b.d_to_m; a.pairs += b.pairs;
            for h in 0..HP_CLASSES { a.m_to_i[h] += b.m_to_i[h]; a.m_to_d[h] += b.m_to_d[h]; a.m_to_m[h] += b.m_to_m[h]; }
            a
        });
        if counts.pairs == 0 { break; }
        p = counts.estimate(&HmmParams::default(), strength);
    }
    p
}

/// The prior weight of a candidate junction from its splice motif (GT-AG, GC-AG, AT-AC, other) and whether it is
/// annotated: pseudo-counts of the mixture's Dirichlet prior.
pub fn site_prior(donor_motif: &[u8], acceptor_motif: &[u8], annotated: bool) -> f64 {
    let up = |s: &[u8]| s.iter().map(|b| b.to_ascii_uppercase()).collect::<Vec<u8>>();
    let (d, a) = (up(donor_motif), up(acceptor_motif));
    let motif = match (d.as_slice(), a.as_slice()) {
        (b"GT", b"AG") => 1.0,
        (b"GC", b"AG") => 0.1,
        (b"AT", b"AC") => 0.02,
        _ => 0.002,
    };
    motif * if annotated { 5.0 } else { 1.0 }
}

/// A junction group's decision: the mixture weight of each candidate (zero for dropped ones) and each read's site.
#[derive(Debug, Clone, PartialEq)]
pub struct GroupFit {
    pub weights: Vec<f64>,
    pub assignment: Vec<usize>,
    /// log-likelihood of the reads under the fitted mixture
    pub log_likelihood: f64,
}

/// EM for the mixture weights over the `alive` candidates, the priors as pseudo-counts; the weights and the reads'
/// log-likelihood under them.
fn fit_weights(ll: &[Vec<f64>], prior: &[f64], alive: &[bool]) -> (Vec<f64>, f64) {
    let m = prior.len();
    let k_alive = alive.iter().filter(|&&a| a).count().max(1);
    let mut weights: Vec<f64> = (0..m).map(|c| if alive[c] { 1.0 / k_alive as f64 } else { 0.0 }).collect();
    let mut expected = vec![0.0; m];
    let log_weights = |w: &[f64]| -> Vec<f64> { w.iter().map(|&x| if x > 0.0 { x.ln() } else { f64::NEG_INFINITY }).collect() };
    for _ in 0..100 {
        expected.iter_mut().for_each(|x| *x = 0.0);
        let lw = log_weights(&weights);
        for row in ll {
            let z = (0..m).map(|c| lw[c] + row[c]).fold(f64::NEG_INFINITY, logsum);
            if !z.is_finite() { continue; }
            for c in 0..m { if alive[c] { expected[c] += (lw[c] + row[c] - z).exp(); } }
        }
        let total: f64 = (0..m).filter(|&c| alive[c]).map(|c| expected[c] + prior[c]).sum();
        let mut shift = 0.0;
        for c in 0..m {
            let w = if alive[c] { (expected[c] + prior[c]) / total } else { 0.0 };
            shift += (w - weights[c]).abs();
            weights[c] = w;
        }
        if shift < 1e-7 { break; }
    }
    let lw = log_weights(&weights);
    let total = ll.iter().map(|row| (0..m).map(|c| lw[c] + row[c]).fold(f64::NEG_INFINITY, logsum)).filter(|z| z.is_finite()).sum();
    (weights, total)
}

/// Fits the mixture over a group's candidates: `ll[r][c]` is read r's log-likelihood under candidate c, `prior[c]` the
/// candidate's pseudo-count. EM for the weights; then the candidate whose removal costs the reads least is dropped
/// while that cost is below `max_cost` nats, unless `keep[c]` is set. Each read then takes its most probable candidate.
pub fn fit_group(ll: &[Vec<f64>], prior: &[f64], keep: &[bool], max_cost: f64) -> GroupFit {
    let m = prior.len();
    let mut alive = vec![true; m];
    let (mut weights, mut log_likelihood) = fit_weights(ll, prior, &alive);
    loop {
        if alive.iter().filter(|&&a| a).count() <= 1 { break; }
        let mut best: Option<(f64, usize, Vec<f64>, f64)> = None;
        let candidates: Vec<usize> = (0..m).filter(|&c| alive[c] && !keep[c]).collect();
        for c in candidates {
            alive[c] = false;
            let (w, total) = fit_weights(ll, prior, &alive);
            alive[c] = true;
            // the log-likelihood the reads lose without the candidate, plus the log of its prior when below one
            let cost = log_likelihood - total + prior[c].max(1e-9).ln().min(0.0);
            if cost < max_cost && best.as_ref().is_none_or(|b| cost < b.0) { best = Some((cost, c, w, total)); }
        }
        let Some((_, c, w, total)) = best else { break };
        alive[c] = false;
        weights = w;
        log_likelihood = total;
    }
    let lw: Vec<f64> = weights.iter().map(|&x| if x > 0.0 { x.ln() } else { f64::NEG_INFINITY }).collect();
    let assignment = ll.iter().map(|row| {
        (0..m).max_by(|&a, &b| (lw[a] + row[a]).partial_cmp(&(lw[b] + row[b])).unwrap()).unwrap_or(0)
    }).collect();
    GroupFit { weights, assignment, log_likelihood }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A deterministic sequence generator for the fixtures.
    fn rng(seed: u64) -> impl FnMut() -> u64 {
        let mut s = seed;
        move || { s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407); s >> 33 }
    }

    /// `t` read with independent substitutions, insertions and deletions at the given rates.
    fn noisy(t: &[u8], sub: f64, ins: f64, del: f64, next: &mut impl FnMut() -> u64) -> Vec<u8> {
        let mut out = Vec::new();
        for &b in t {
            if (next() % 1_000_000) as f64 / 1e6 < ins { out.push(b"ACGT"[(next() % 4) as usize]); }
            if (next() % 1_000_000) as f64 / 1e6 < del { continue; }
            if (next() % 1_000_000) as f64 / 1e6 < sub {
                let at = b"ACGT".iter().position(|&x| x == b).unwrap() as u64;
                out.push(b"ACGT"[((at + 1 + next() % 3) % 4) as usize]);
            } else {
                out.push(b);
            }
        }
        out
    }

    #[test]
    fn the_true_template_is_the_most_likely() {
        let mut next = rng(7);
        let t: Vec<u8> = (0..60).map(|_| b"ACGT"[(next() % 4) as usize]).collect();
        let r = noisy(&t, 0.03, 0.02, 0.03, &mut next);
        let mut shifted = t.clone();
        shifted.drain(28..31);                               // three bases fewer, as a site 3 bp further in
        shifted.extend_from_slice(b"ACG");
        let p = HmmParams::default();
        assert!(log_likelihood(&r, &t, &p) > log_likelihood(&r, &shifted, &p) + 2.0);
        // identical sequences: likelihood near the no-error path
        let exact = log_likelihood(&t, &t, &p);
        assert!(exact > 60.0 * (1.0 - 0.03f64).ln() + 60.0 * (1.0 - 0.02 - 0.03f64).ln() - 1.0, "{exact}");
    }

    #[test]
    fn the_band_changes_nothing_for_reads_near_the_diagonal() {
        let mut next = rng(3);
        let p = HmmParams::default();
        for _ in 0..20 {
            let t: Vec<u8> = (0..60).map(|_| b"ACGT"[(next() % 4) as usize]).collect();
            let r = noisy(&t, 0.03, 0.02, 0.03, &mut next);
            let (full, banded) = (log_likelihood(&r, &t, &p), log_likelihood_banded(&r, &t, &p, 12));
            assert!((full - banded).abs() < 1e-6 * full.abs().max(1.0) + 1e-3, "{full} {banded}");
        }
    }

    #[test]
    fn the_scaled_forward_sum_equals_the_log_space_one() {
        let mut next = rng(5);
        let p = HmmParams::default();
        for _ in 0..50 {
            let t: Vec<u8> = (0..70).map(|_| b"ACGT"[(next() % 4) as usize]).collect();
            let r = noisy(&t, 0.05, 0.03, 0.04, &mut next);
            let (full, banded) = (log_likelihood(&r, &t, &p), log_likelihood_banded(&r, &t, &p, 40));
            assert!((full - banded).abs() < 1e-9 * full.abs().max(1.0), "{full} {banded}");
        }
    }

    #[test]
    fn homopolymer_classes_count_the_run() {
        assert_eq!(homopolymer_classes(b"ACCCGTTTTTTTTA"), vec![0, 2, 2, 2, 0, 5, 5, 5, 5, 5, 5, 5, 5, 0]);
    }

    #[test]
    fn baum_welch_recovers_the_error_rates() {
        let mut next = rng(11);
        let pairs: Vec<(Vec<u8>, Vec<u8>)> = (0..300).map(|_| {
            let t: Vec<u8> = (0..60).map(|_| b"ACGT"[(next() % 4) as usize]).collect();
            (noisy(&t, 0.06, 0.01, 0.01, &mut next), t)
        }).collect();
        let fitted = train(&pairs, &HmmParams::default(), 5, 10.0);
        assert!((fitted.mismatch - 0.06).abs() < 0.015, "{fitted:?}");
        assert!(fitted.open_del[0] < 0.02 && fitted.open_ins[0] < 0.02, "{fitted:?}");
    }

    #[test]
    fn a_group_takes_back_reads_the_aligner_moved_but_keeps_two_real_sites() {
        // candidate 0: the true site; candidate 1: a site 3 bp away. Reads of site 0 favour it mildly each.
        let one_site: Vec<Vec<f64>> = (0..20).map(|i| if i % 3 == 0 { vec![-50.0, -49.5] } else { vec![-50.0, -52.0] }).collect();
        let fit = fit_group(&one_site, &[1.0, 1.0], &[false, false], 3.0);
        assert_eq!(fit.weights[1], 0.0, "{fit:?}");
        assert!(fit.assignment.iter().all(|&c| c == 0));
        // two real sites: each read favours its own by 10 nats
        let two: Vec<Vec<f64>> = (0..20).map(|i| if i < 8 { vec![-60.0, -50.0] } else { vec![-50.0, -60.0] }).collect();
        let fit = fit_group(&two, &[1.0, 1.0], &[false, false], 3.0);
        assert!(fit.weights[0] > 0.4 && fit.weights[1] > 0.3, "{fit:?}");
        assert_eq!(fit.assignment.iter().filter(|&&c| c == 1).count(), 8);
        // two sites the reads cannot tell apart: the one with the weaker motif goes
        let tie: Vec<Vec<f64>> = (0..30).map(|_| vec![-50.0, -50.0]).collect();
        let fit = fit_group(&tie, &[0.002, 1.0], &[false, false], 3.0);
        assert_eq!(fit.weights[0], 0.0, "{fit:?}");
        let fit = fit_group(&tie, &[1.0, 0.002], &[false, false], 3.0);
        assert_eq!(fit.weights[1], 0.0, "{fit:?}");
        // a site to keep stays, though its reads do not need it
        let fit = fit_group(&one_site, &[1.0, 1.0], &[false, true], 3.0);
        assert!(fit.weights[1] > 0.0, "{fit:?}");
    }

    #[test]
    fn canonical_annotated_sites_weigh_most() {
        assert!(site_prior(b"GT", b"AG", true) > site_prior(b"GT", b"AG", false));
        assert!(site_prior(b"GT", b"AG", false) > site_prior(b"GC", b"AG", false));
        assert!(site_prior(b"GC", b"AG", false) > site_prior(b"AT", b"AC", false));
        assert!(site_prior(b"at", b"ac", false) > site_prior(b"CT", b"AC", false));
    }
}
