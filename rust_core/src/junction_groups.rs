// Junction groups: where the reads of one junction put its sites a few bases apart, the group decides from all their
// bases (junction_hmm). A group is the introns with at least `min_support` reads, and the annotated introns, whose
// donors and acceptors are both within `window` bp of another member's; it is contested when two of its introns have
// the support and the second holds `CONTEST_SHARE` of the group's reads. Before annotation every dataset is read twice:
// once to count the introns the spliced reads carry and to sample reads at annotated introns for training the pair-HMM,
// once to score the reads of contested groups against every candidate site of their group. Each group's mixture is
// then fitted, its annotated sites kept, and every scored read decided; annotation reads the decisions.

use std::collections::BTreeMap;
use std::sync::atomic::{AtomicU64, AtomicU8, Ordering};
use std::sync::{Mutex, RwLock};

use rust_htslib::bam;
use rustc_hash::{FxHashMap, FxHashSet};

use crate::data_structures::SharedRefs;
use crate::genome_ref::GenomeRef;
use crate::junction_hmm::{fit_group, forward_banded, homopolymer_classes, site_prior, train, HmmParams, ProbParams};
use crate::microexon_repair::query_boundary;
use crate::types::EdgeType;

/// Exonic bases scored on each side, beyond the outermost candidate site.
const FLANK: i64 = 20;
/// Candidates a group may hold; the best-supported ones are kept.
const MAX_CANDIDATES: usize = 8;
/// Share of a group's reads its second-best intron needs for the group to be contested.
const CONTEST_SHARE: f64 = 0.05;
/// Log-likelihood (nats) the group's reads must lose without a site for the site to stay.
const DROP_COST: f64 = 3.0;
/// Posterior below which a read at a kept site moves to the kept site most probable for it.
const MOVE_POSTERIOR: f64 = 0.01;
/// Cells either side of the diagonal the forward sum keeps.
const BAND: usize = 12;
/// Read-template pairs sampled for training the pair-HMM, and one read in `TRAIN_EVERY` at an annotated intron taken;
/// of more, the pairs with the smallest sampling keys are kept.
const TRAIN_PAIRS: usize = 20_000;
const TRAIN_EVERY: u64 = 37;
const SHARDS: usize = 64;
/// Bases by which the parse may move an aligned intron before the groups see it (the junction-indel rule, microindel
/// absorption): the scoring pass parses only reads with an aligned intron boundary this close to a contested site.
const PARSE_SLACK: i64 = 30;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum Phase { Count = 0, Score = 1, Decide = 2 }

/// A contested group: its candidate introns (first, last intron base, 1-based), their priors, and the spliced
/// templates of every candidate over the group's window [lo, hi) (0-based, `FLANK` exonic bases outside the
/// outermost sites).
struct Group {
    chrom: u32,
    candidates: Vec<(i64, i64)>,
    priors: Vec<f64>,
    annotated: Vec<bool>,
    lo: i64,
    hi: i64,
    templates: Vec<Vec<u8>>,
    /// each template's homopolymer classes
    hp: Vec<Vec<usize>>,
}

/// A fitted group, for the junction table: its candidates' reads as aligned and as decided.
struct GroupSummary { group: u32, chrom: u32, observed: Vec<u32>, decided: Vec<u32>, weights: Vec<f64>, other: u32 }

/// One scored read intron: the read (dataset and name hashed), its intron, its group, which candidate it is (if any),
/// and its log-likelihood under each candidate.
struct Scored { read: u64, first: i64, last: i64, group: u32, observed: Option<u8>, ll: Vec<f32>, name: Option<String> }

pub struct JunctionGroups {
    phase: AtomicU8,
    window: i64,
    min_support: u32,
    chroms: RwLock<FxHashMap<String, u32>>,
    /// annotated introns as (chromosome, first intron base, last intron base), 1-based
    annotated: FxHashSet<(String, i64, i64)>,
    counts: Vec<Mutex<FxHashMap<(u32, i64, i64), u32>>>,
    /// sampled (key, read bases, template): a max-heap on the key holding the `TRAIN_PAIRS` smallest
    training: Mutex<std::collections::BinaryHeap<(u64, Vec<u8>, Vec<u8>)>>,
    /// chromosome -> first intron base -> (last intron base, group), for members of contested groups
    by_first: RwLock<FxHashMap<u32, BTreeMap<i64, Vec<(i64, u32)>>>>,
    groups: RwLock<Vec<Group>>,
    /// chromosome -> every donor and acceptor of the contested groups' candidates, sorted
    sites: RwLock<FxHashMap<u32, Vec<i64>>>,
    params: RwLock<HmmParams>,
    /// parameters given by the user (not trained), and the read-template pairs the trained ones came from
    preloaded: bool,
    trained_on: AtomicU64,
    summaries: Mutex<Vec<GroupSummary>>,
    scores: Vec<Mutex<Vec<Scored>>>,
    /// (read, first, last) -> the intron the read takes
    decisions: RwLock<FxHashMap<(u64, i64, i64), (i64, i64)>>,
    scored: AtomicU64,
    decided: AtomicU64,
    moved: AtomicU64,
    /// With ISOFORGE_JUNCTION_TRACE=<file>: every scored read intron, its group and its decision, written there.
    trace: bool,
}

fn shard_of<T: std::hash::Hash>(x: &T) -> usize {
    use std::hash::Hasher;
    let mut h = rustc_hash::FxHasher::default();
    x.hash(&mut h);
    (h.finish() as usize) % SHARDS
}

fn read_key(dataset: &str, qname: &[u8]) -> u64 {
    use std::hash::Hasher;
    let mut h = rustc_hash::FxHasher::default();
    h.write(dataset.as_bytes());
    h.write_u8(0);
    h.write(qname);
    h.finish()
}

impl JunctionGroups {
    /// Groups over the annotated introns of `refs`; `window` bp at each site, members need `min_support` reads.
    pub fn new(refs: &SharedRefs, window: i64, min_support: u32) -> Self {
        let pos: FxHashMap<_, _> = refs.location_dict.chromosomes.values()
            .flat_map(|c| c.positions.iter().map(|&(p, v)| (v, p))).collect();
        let chrom_of: FxHashMap<_, _> = refs.location_dict.chromosomes.iter()
            .flat_map(|(chrom, c)| c.positions.iter().map(move |&(_, v)| (v, chrom.clone()))).collect();
        let mut annotated = FxHashSet::default();
        for e in refs.edge_dict.values().filter(|e| e.edge_type == EdgeType::Intron) {
            if let (Some(&p1), Some(&p2), Some(chrom)) = (pos.get(&e.v1), pos.get(&e.v2), chrom_of.get(&e.v1)) {
                annotated.insert((chrom.clone(), p1.min(p2) + 1, p1.max(p2) - 1));
            }
        }
        Self {
            phase: AtomicU8::new(Phase::Count as u8), window, min_support,
            chroms: RwLock::new(FxHashMap::default()), annotated,
            counts: (0..SHARDS).map(|_| Mutex::new(FxHashMap::default())).collect(),
            training: Mutex::new(std::collections::BinaryHeap::new()),
            by_first: RwLock::new(FxHashMap::default()),
            groups: RwLock::new(Vec::new()),
            sites: RwLock::new(FxHashMap::default()),
            params: RwLock::new(HmmParams::default()),
            preloaded: false, trained_on: AtomicU64::new(0), summaries: Mutex::new(Vec::new()),
            scores: (0..SHARDS).map(|_| Mutex::new(Vec::new())).collect(),
            decisions: RwLock::new(FxHashMap::default()),
            scored: AtomicU64::new(0), decided: AtomicU64::new(0), moved: AtomicU64::new(0),
            trace: std::env::var_os("ISOFORGE_JUNCTION_TRACE").is_some(),
        }
    }

    /// Uses `params` instead of training the pair-HMM on the run's reads.
    pub fn with_params(mut self, params: HmmParams) -> Self {
        *self.params.get_mut().unwrap() = params;
        self.preloaded = true;
        self
    }

    pub fn phase(&self) -> Phase {
        match self.phase.load(Ordering::Relaxed) { 0 => Phase::Count, 1 => Phase::Score, _ => Phase::Decide }
    }
    fn set_phase(&self, p: Phase) { self.phase.store(p as u8, Ordering::Relaxed); }
    /// Whether reads are only being counted or scored (the parse stops after this step).
    pub fn is_prepass(&self) -> bool { self.phase() != Phase::Decide }

    fn chrom_id(&self, chrom: &str) -> u32 {
        if let Some(&c) = self.chroms.read().unwrap().get(chrom) { return c; }
        let mut w = self.chroms.write().unwrap();
        let n = w.len() as u32;
        *w.entry(chrom.to_string()).or_insert(n)
    }

    /// Whether the read is needed by the scoring pass: an aligned intron with a boundary within `window` +
    /// `PARSE_SLACK` bp of a contested site. The other reads leave the scoring pass before they are parsed.
    pub fn near_contested(&self, chrom: &str, record: &bam::Record) -> bool {
        let Some(&c) = self.chroms.read().unwrap().get(chrom) else { return false };
        let sites = self.sites.read().unwrap();
        let Some(s) = sites.get(&c) else { return false };
        let reach = self.window + PARSE_SLACK;
        let near = |b: i64| { let i = s.partition_point(|&x| x < b - reach); i < s.len() && s[i] <= b + reach };
        let mut pos = record.pos() + 1;
        for op in record.cigar().iter() {
            use rust_htslib::bam::record::Cigar::*;
            match op {
                Match(n) | Equal(n) | Diff(n) | Del(n) => pos += *n as i64,
                RefSkip(n) => {
                    if near(pos) || near(pos + *n as i64 - 1) { return true; }
                    pos += *n as i64;
                }
                _ => {}
            }
        }
        false
    }

    /// The contested group a read intron belongs to: a member within `window` bp at both sites.
    fn group_of(&self, chrom: u32, first: i64, last: i64) -> Option<u32> {
        let map = self.by_first.read().unwrap();
        map.get(&chrom)?.range(first - self.window..=first + self.window)
            .flat_map(|(_, v)| v.iter())
            .find(|&&(l, _)| (l - last).abs() <= self.window)
            .map(|&(_, g)| g)
    }

    /// The read's bases over [lo, hi), when its alignment covers the window and holds no intron there but intron `k`.
    fn read_window(record: &bam::Record, introns: &[i64], k: usize, lo: i64, hi: i64) -> Option<Vec<u8>> {
        if lo < record.pos() || hi > record.cigar().end_pos() { return None; }
        for j in 0..introns.len() / 2 {
            if j != k && introns[2 * j + 1] >= lo && introns[2 * j] <= hi { return None; }
        }
        let (qlo, qhi) = (query_boundary(record, lo)?, query_boundary(record, hi)?);
        if qhi <= qlo || qhi - qlo > 4 * (hi - lo) as usize { return None; }
        Some(record.seq().as_bytes()[qlo..qhi].to_vec())
    }

    fn template(genome: &GenomeRef, chrom: &str, lo: i64, hi: i64, c: (i64, i64)) -> Option<Vec<u8>> {
        let mut t = genome.fetch_seq(chrom, lo, c.0 - 1).ok()?;
        t.extend(genome.fetch_seq(chrom, c.1, hi).ok()?);
        Some(t)
    }

    /// The parse step: counts (Count), scores (Score) or applies the decisions (Decide) for the read's introns. `introns`
    /// holds first and last intron bases, 1-based, in genomic order. Returns them, moved where the group decided, and
    /// which introns a group decided (those the per-read realignment then leaves alone).
    pub fn apply(&self, record: &bam::Record, dataset: &str, chrom: &str, genome: &GenomeRef, introns: &[i64]) -> (Vec<i64>, Vec<bool>) {
        let n = introns.len() / 2;
        let mut out = introns.to_vec();
        let mut handled = vec![false; n];
        let cid = self.chrom_id(chrom);
        match self.phase() {
            Phase::Count => {
                for k in 0..n {
                    let (f, l) = (introns[2 * k], introns[2 * k + 1]);
                    let key = (cid, f, l);
                    *self.counts[shard_of(&key)].lock().unwrap().entry(key).or_insert(0) += 1;
                    let key = {
                        use std::hash::Hasher;
                        let mut h = rustc_hash::FxHasher::default();
                        h.write(record.qname());
                        h.write_i64(f);
                        h.finish()
                    };
                    let room = |tr: &std::collections::BinaryHeap<(u64, Vec<u8>, Vec<u8>)>| {
                        tr.len() < TRAIN_PAIRS || tr.peek().is_some_and(|top| key < top.0)
                    };
                    if key % TRAIN_EVERY == 0 && self.annotated.contains(&(chrom.to_string(), f, l)) && room(&self.training.lock().unwrap()) {
                        let (lo, hi) = (f - 1 - FLANK, l + FLANK);
                        if let (Some(q), Some(t)) = (Self::read_window(record, introns, k, lo, hi), Self::template(genome, chrom, lo, hi, (f, l))) {
                            let mut tr = self.training.lock().unwrap();
                            tr.push((key, q, t));
                            if tr.len() > TRAIN_PAIRS { tr.pop(); }
                        }
                    }
                }
            }
            Phase::Score => {
                let pp = ProbParams::new(&self.params.read().unwrap());
                let rk = read_key(dataset, record.qname());
                for k in 0..n {
                    let (f, l) = (introns[2 * k], introns[2 * k + 1]);
                    let Some(g) = self.group_of(cid, f, l) else { continue };
                    let groups = self.groups.read().unwrap();
                    let gr = &groups[g as usize];
                    let Some(q) = Self::read_window(record, introns, k, gr.lo, gr.hi) else { continue };
                    let ll: Vec<f32> = gr.templates.iter().zip(&gr.hp).map(|(t, h)| forward_banded(&q, t, h, &pp, BAND) as f32).collect();
                    let observed = gr.candidates.iter().position(|&c| c == (f, l)).map(|x| x as u8);
                    drop(groups);
                    let name = self.trace.then(|| String::from_utf8_lossy(record.qname()).to_string());
                    self.scores[g as usize % SHARDS].lock().unwrap().push(Scored { read: rk, first: f, last: l, group: g, observed, ll, name });
                    self.scored.fetch_add(1, Ordering::Relaxed);
                }
            }
            Phase::Decide => {
                let rk = read_key(dataset, record.qname());
                let d = self.decisions.read().unwrap();
                for k in 0..n {
                    let (f, l) = (introns[2 * k], introns[2 * k + 1]);
                    let Some(&(nf, nl)) = d.get(&(rk, f, l)) else { continue };
                    out[2 * k] = nf;
                    out[2 * k + 1] = nl;
                    handled[k] = true;
                    self.decided.fetch_add(1, Ordering::Relaxed);
                    if (nf, nl) != (f, l) { self.moved.fetch_add(1, Ordering::Relaxed); }
                }
            }
        }
        // a moved intron may not cross its neighbours, nor leave the read's alignment
        if out.windows(2).any(|w| w[0] > w[1]) || out.first().is_some_and(|&x| x <= record.pos() + 1)
            || out.last().is_some_and(|&x| x >= record.cigar().end_pos()) {
            return (introns.to_vec(), vec![false; n]);
        }
        (out, handled)
    }

    /// After counting: the contested groups with their candidates, priors and templates, and the pair-HMM's parameters
    /// trained on the sampled reads at annotated introns. Then scoring begins.
    pub fn build_groups(&self, genome: &GenomeRef) -> usize {
        let mut counts: FxHashMap<(u32, i64, i64), u32> = FxHashMap::default();
        for s in &self.counts { counts.extend(s.lock().unwrap().drain()); }
        let chroms = self.chroms.read().unwrap().clone();
        let name_of: FxHashMap<u32, String> = chroms.iter().map(|(n, &c)| (c, n.clone())).collect();
        // members: supported or annotated introns, per chromosome
        let mut per_chrom: FxHashMap<u32, Vec<(i64, i64, u32, bool)>> = FxHashMap::default();
        for (&(c, f, l), &n) in &counts {
            let ann = self.annotated.contains(&(name_of[&c].clone(), f, l));
            if n >= self.min_support || ann { per_chrom.entry(c).or_default().push((f, l, n, ann)); }
        }
        // unobserved annotated introns next to observed ones compete too
        for (chrom, f, l) in &self.annotated {
            let Some(&c) = chroms.get(chrom) else { continue };
            if counts.contains_key(&(c, *f, *l)) { continue; }
            if let Some(v) = per_chrom.get_mut(&c) { v.push((*f, *l, 0, true)); }
        }
        let dinuc = |chrom: &str, a: i64| genome.fetch_seq(chrom, a, a + 2).unwrap_or_default();
        let revcomp = |s: &[u8]| s.iter().rev().map(|b| match b.to_ascii_uppercase() { b'A' => b'T', b'C' => b'G', b'G' => b'C', b'T' => b'A', x => x }).collect::<Vec<u8>>();
        let mut by_first: FxHashMap<u32, BTreeMap<i64, Vec<(i64, u32)>>> = FxHashMap::default();
        let mut groups: Vec<Group> = Vec::new();
        let mut keys: Vec<u32> = per_chrom.keys().copied().collect();
        keys.sort_unstable_by(|a, b| name_of[a].cmp(&name_of[b]));
        for c in keys {
            let mut v = per_chrom.remove(&c).unwrap();
            let chrom = &name_of[&c];
            v.sort_unstable();
            let mut parent: Vec<usize> = (0..v.len()).collect();
            fn find(p: &mut [usize], x: usize) -> usize { let mut r = x; while p[r] != r { r = p[r]; } let mut y = x; while p[y] != r { let n = p[y]; p[y] = r; y = n; } r }
            for i in 0..v.len() {
                let mut j = i + 1;
                while j < v.len() && v[j].0 - v[i].0 <= self.window {
                    if (v[j].1 - v[i].1).abs() <= self.window { let (a, b) = (find(&mut parent, i), find(&mut parent, j)); if a != b { parent[a] = b; } }
                    j += 1;
                }
            }
            let mut comps: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
            for i in 0..v.len() { let r = find(&mut parent, i); comps.entry(r).or_default().push(i); }
            for (_, mut members) in comps {
                // contested: two supported introns, the second with a share of the group's reads
                let mut support: Vec<u32> = members.iter().map(|&i| v[i].2).filter(|&n| n >= self.min_support).collect();
                support.sort_unstable_by(|a, b| b.cmp(a));
                let total: u32 = members.iter().map(|&i| v[i].2).sum();
                if support.len() < 2 || (support[1] as f64) < CONTEST_SHARE * total as f64 { continue; }
                members.sort_by_key(|&i| (std::cmp::Reverse(v[i].2), std::cmp::Reverse(v[i].3), v[i].0, v[i].1));
                members.truncate(MAX_CANDIDATES);
                members.sort_by_key(|&i| (v[i].0, v[i].1));
                let candidates: Vec<(i64, i64)> = members.iter().map(|&i| (v[i].0, v[i].1)).collect();
                let lo = candidates.iter().map(|x| x.0).min().unwrap() - 1 - FLANK;
                let hi = candidates.iter().map(|x| x.1).max().unwrap() + FLANK;
                let templates: Option<Vec<Vec<u8>>> = candidates.iter().map(|&x| Self::template(genome, chrom, lo, hi, x)).collect();
                let Some(templates) = templates else { continue };
                let priors: Vec<f64> = members.iter().map(|&i| {
                    let (f, l, _, ann) = v[i];
                    let (d, a) = (dinuc(chrom, f - 1), dinuc(chrom, l - 2));
                    site_prior(&d, &a, ann).max(site_prior(&revcomp(&a), &revcomp(&d), ann))
                }).collect();
                let annotated: Vec<bool> = members.iter().map(|&i| v[i].3).collect();
                let g = groups.len() as u32;
                for &(f, l) in &candidates { by_first.entry(c).or_default().entry(f).or_default().push((l, g)); }
                let hp = templates.iter().map(|t| homopolymer_classes(t)).collect();
                groups.push(Group { chrom: c, candidates, priors, annotated, lo, hi, templates, hp });
            }
        }
        let n_groups = groups.len();
        let mut sites: FxHashMap<u32, Vec<i64>> = FxHashMap::default();
        for (&c, m) in &by_first {
            let v = sites.entry(c).or_default();
            for (&f, ls) in m { v.push(f); v.extend(ls.iter().map(|&(l, _)| l)); }
            v.sort_unstable();
            v.dedup();
        }
        *self.sites.write().unwrap() = sites;
        *self.by_first.write().unwrap() = by_first;
        *self.groups.write().unwrap() = groups;
        let pairs: Vec<(Vec<u8>, Vec<u8>)> = std::mem::take(&mut *self.training.lock().unwrap())
            .into_sorted_vec().into_iter().map(|(_, q, t)| (q, t)).collect();
        if !self.preloaded && pairs.len() >= 200 {
            *self.params.write().unwrap() = train(&pairs, &HmmParams::default(), 5, 20.0);
            self.trained_on.store(pairs.len() as u64, Ordering::Relaxed);
        }
        self.set_phase(Phase::Score);
        let p = self.params.read().unwrap();
        let source = if self.preloaded { "given".to_string() } else { format!("from {} reads at annotated introns", pairs.len()) };
        eprintln!("[IsoForge] Junction groups: {} contested of {} distinct introns; pair-HMM {}: mismatch {:.3}, \
                   insertion/deletion opening {:.3}/{:.3} (no run) to {:.3}/{:.3} (runs of 6+)",
                  n_groups, counts.len(), source, p.mismatch, p.open_ins[0], p.open_del[0],
                  p.open_ins[crate::junction_hmm::HP_CLASSES - 1], p.open_del[crate::junction_hmm::HP_CLASSES - 1]);
        n_groups
    }

    /// After scoring: each contested group's mixture and each scored read's decision. Then annotation reads them.
    pub fn fit(&self) {
        let groups = self.groups.read().unwrap();
        let mut per_group: Vec<Vec<Scored>> = (0..groups.len()).map(|_| Vec::new()).collect();
        for s in &self.scores { for x in s.lock().unwrap().drain(..) { per_group[x.group as usize].push(x); } }
        // each group's rows sorted by read, then intron
        for rows in per_group.iter_mut() { rows.sort_unstable_by_key(|r| (r.read, r.first, r.last)); }
        use rayon::prelude::*;
        type Fitted = (Vec<((u64, i64, i64), (i64, i64))>, usize, Vec<String>, Option<GroupSummary>);
        let results: Vec<Fitted> = per_group.par_iter().enumerate().map(|(g, rows)| {
            if rows.is_empty() { return (Vec::new(), 0, Vec::new(), None); }
            let gr = &groups[g];
            let ll: Vec<Vec<f64>> = rows.iter().map(|r| r.ll.iter().map(|&x| x as f64).collect()).collect();
            let w = fit_group(&ll, &gr.priors, &gr.annotated, DROP_COST).weights;
            let live: Vec<usize> = (0..w.len()).filter(|&c| w[c] > 0.0).collect();
            let dropped = w.len() - live.len();
            let decisions = rows.iter().zip(&ll).map(|(r, l)| {
                let best = live.iter().copied().max_by(|&a, &b| (w[a].ln() + l[a]).partial_cmp(&(w[b].ln() + l[b])).unwrap()).unwrap_or(0);
                let pick = match r.observed.map(|o| o as usize).filter(|o| w[*o] > 0.0) {
                    // at a kept site: move only when the site is improbable for the read under the mixture
                    Some(o) => {
                        let z = live.iter().map(|&c| w[c].ln() + l[c]).fold(f64::NEG_INFINITY, crate::junction_hmm::logsum);
                        if best != o && (w[o].ln() + l[o] - z).exp() < MOVE_POSTERIOR { best } else { o }
                    }
                    // at a dropped site or none: the most probable kept site
                    None => best,
                };
                ((r.read, r.first, r.last), gr.candidates[pick])
            }).collect::<Vec<_>>();
            let trace: Vec<String> = if self.trace {
                let cands = gr.candidates.iter().map(|c| format!("{}-{}", c.0, c.1)).collect::<Vec<_>>().join(",");
                let ws = w.iter().map(|x| format!("{x:.4}")).collect::<Vec<_>>().join(",");
                let pr = gr.priors.iter().map(|x| format!("{x:.4}")).collect::<Vec<_>>().join(",");
                rows.iter().zip(&decisions).map(|(r, (_, d))| format!("{}\t{}\t{}\t{}\t{}-{}\t{}\t{}\t{}\t{}", r.name.as_deref().unwrap_or(""), g,
                    r.first, r.last, d.0, d.1, cands, ws, pr, r.ll.iter().map(|x| format!("{x:.2}")).collect::<Vec<_>>().join(","))).collect()
            } else { Vec::new() };
            let m = gr.candidates.len();
            let (mut observed, mut decided, mut other) = (vec![0u32; m], vec![0u32; m], 0u32);
            for (r, (_, d)) in rows.iter().zip(&decisions) {
                match r.observed { Some(o) => observed[o as usize] += 1, None => other += 1 }
                if let Some(i) = gr.candidates.iter().position(|c| c == d) { decided[i] += 1; }
            }
            let summary = GroupSummary { group: g as u32, chrom: gr.chrom, observed, decided, weights: w.clone(), other };
            (decisions, dropped, trace, Some(summary))
        }).collect();
        let mut map: FxHashMap<(u64, i64, i64), (i64, i64)> = FxHashMap::default();
        let (mut dropped, mut to_move) = (0usize, 0usize);
        let with_reads = per_group.iter().filter(|r| !r.is_empty()).count();
        let mut trace_rows: Vec<String> = Vec::new();
        let mut summaries = Vec::new();
        for (d, x, t, s) in results {
            dropped += x;
            for (k, v) in d { if (k.1, k.2) != v { to_move += 1; } map.insert(k, v); }
            trace_rows.extend(t);
            summaries.extend(s);
        }
        *self.summaries.lock().unwrap() = summaries;
        if let Some(path) = std::env::var_os("ISOFORGE_JUNCTION_TRACE") {
            let text = format!("read\tgroup\tfirst\tlast\tdecided\tcandidates\tweights\tpriors\tll\n{}\n", trace_rows.join("\n"));
            if let Err(e) = std::fs::write(&path, text) { eprintln!("[IsoForge] junction trace not written: {e}"); }
        }
        eprintln!("[IsoForge] Junction groups fitted: {} with reads, {} read introns scored, {} candidate sites dropped, \
                   {} read introns to move", with_reads, self.scored.load(Ordering::Relaxed), dropped, to_move);
        *self.decisions.write().unwrap() = map;
        drop(groups);
        self.set_phase(Phase::Decide);
    }

    /// The pair-HMM parameters (JSON) and one line per candidate site of each fitted group: its annotation, prior,
    /// mixture weight, and reads as aligned and as decided.
    pub fn write_outputs(&self, model: &std::path::Path, table: &std::path::Path) -> std::io::Result<()> {
        let params = self.params.read().unwrap().clone();
        let doc = serde_json::json!({
            "params": params,
            "source": if self.preloaded { "given" } else { "trained" },
            "training_pairs": self.trained_on.load(Ordering::Relaxed),
        });
        std::fs::write(model, serde_json::to_string_pretty(&doc).unwrap_or_default() + "\n")?;
        let names: FxHashMap<u32, String> = self.chroms.read().unwrap().iter().map(|(n, &c)| (c, n.clone())).collect();
        let groups = self.groups.read().unwrap();
        let mut out = String::from("group\tchromosome\tintron_start\tintron_end\tannotated\tprior\tweight\treads_aligned\treads_decided\treads_elsewhere\n");
        let mut summaries = self.summaries.lock().unwrap();
        summaries.sort_by_key(|s| s.group);
        for s in summaries.iter() {
            let gr = &groups[s.group as usize];
            for (i, c) in gr.candidates.iter().enumerate() {
                out.push_str(&format!("{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{}\t{}\t{}\n", s.group, names.get(&s.chrom).map(|x| x.as_str()).unwrap_or(""),
                                      c.0, c.1, gr.annotated[i], gr.priors[i], s.weights[i], s.observed[i], s.decided[i], s.other));
            }
        }
        std::fs::write(table, out)
    }

    pub fn log_decisions(&self) {
        eprintln!("[IsoForge] Junction groups: {} read introns decided, {} moved to another site of their group",
                  self.decided.load(Ordering::Relaxed), self.moved.load(Ordering::Relaxed));
    }
}

impl std::fmt::Debug for JunctionGroups {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "JunctionGroups({:?}, {} groups)", self.phase(), self.groups.read().map(|g| g.len()).unwrap_or(0))
    }
}
