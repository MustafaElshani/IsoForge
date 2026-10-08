// The read hypergraph: exon segments as nodes, each distinct read class as one directed, weighted hyperedge.
//
// - Per gene, the exon boundaries of its reported models cut exonic sequence into segments (nodes), numbered 5' to 3'.
//   A model is the ordered list of the segments its exons cover; a read, the part of its model between its own ends.
// - Touching segments that every model includes or skips together form one exon unit; a unit is alternative when
//   some model includes it and another skips it while covering exons on both sides.
// - Splicing coordination: for a pair of alternative units, the reads spanning both are tabulated by inclusion and
//   the 2x2 table is tested for dependence (Fisher's exact test, two-sided; Benjamini-Hochberg over all pairs).

use std::fs::File;
use std::io::{BufWriter, Write};

use rustc_hash::FxHashMap;

use crate::errors::IsoForgeResult;
use crate::post_processing::PostDB;

/// One gene's segments in transcript-oriented coordinates (x = position on +, -position on -).
struct GeneGraph {
    /// Segment intervals [lo, hi], sorted 5' to 3'.
    nodes: Vec<(i64, i64)>,
    /// Model -> the segment indices its exons cover, in order.
    models: FxHashMap<i64, Vec<u32>>,
    /// Segment -> exon unit; units are numbered 5' to 3'.
    unit_of: Vec<u32>,
    /// Unit -> its interval and the exon (run of touching units) it belongs to.
    units: Vec<((i64, i64), u32)>,
    /// Unit -> some model's exon starts at its 5' boundary; some model's exon ends at its 3' boundary.
    exon_start: Vec<bool>,
    exon_end: Vec<bool>,
    /// Alternative units.
    alternative: Vec<u32>,
}

/// A read class: gene, the ordered segments it covers, and its weight per dataset.
pub struct Hyperedge {
    pub gene: i64,
    pub nodes: Vec<u32>,
    pub weights: Vec<f64>,
}

pub struct ReadHypergraph {
    genes: FxHashMap<i64, GeneGraph>,
    pub hyperedges: Vec<Hyperedge>,
    pub datasets: Vec<String>,
}

fn oriented(a: i64, b: i64, plus: bool) -> (i64, i64) {
    let s = if plus { 1 } else { -1 };
    let (x, y) = (s * a, s * b);
    (x.min(y), x.max(y))
}

impl GeneGraph {
    fn new(ms: &[(i64, Vec<(i64, i64)>)]) -> Self {
        // boundaries: every exon start and end + 1, so segments are half-open cuts of exonic sequence
        let mut cuts: Vec<i64> = ms.iter().flat_map(|(_, ex)| ex.iter().flat_map(|&(a, b)| [a, b + 1])).collect();
        cuts.sort_unstable();
        cuts.dedup();
        // a segment exists where some model's exon covers it
        let mut nodes = Vec::new();
        for w in cuts.windows(2) {
            let (lo, hi) = (w[0], w[1] - 1);
            if ms.iter().any(|(_, ex)| ex.iter().any(|&(a, b)| a <= lo && hi <= b)) { nodes.push((lo, hi)); }
        }
        let mut models = FxHashMap::default();
        let mut order: Vec<i64> = Vec::with_capacity(ms.len());
        for (tid, ex) in ms {
            let covered: Vec<u32> = nodes.iter().enumerate()
                .filter(|(_, &(lo, hi))| ex.iter().any(|&(a, b)| a <= lo && hi <= b))
                .map(|(i, _)| i as u32).collect();
            models.insert(*tid, covered);
            order.push(*tid);
        }
        // inclusion pattern of each segment over the models; touching segments with one pattern are one unit
        let pattern = |n: u32| -> Vec<bool> { order.iter().map(|t| models[t].binary_search(&n).is_ok()).collect() };
        let mut unit_of = Vec::with_capacity(nodes.len());
        let mut units: Vec<((i64, i64), u32)> = Vec::new();
        let mut prev: Option<Vec<bool>> = None;
        for (i, &(lo, hi)) in nodes.iter().enumerate() {
            let p = pattern(i as u32);
            let touches = i > 0 && nodes[i - 1].1 + 1 == lo;
            if touches && prev.as_ref() == Some(&p) {
                units.last_mut().unwrap().0 .1 = hi;
            } else {
                let exon = match units.last() { Some(&(_, e)) if touches => e, Some(&(_, e)) => e + 1, None => 0 };
                units.push(((lo, hi), exon));
            }
            unit_of.push(units.len() as u32 - 1);
            prev = Some(p);
        }
        let unit_lists: Vec<Vec<u32>> = order.iter().map(|t| {
            let mut u: Vec<u32> = models[t].iter().map(|&n| unit_of[n as usize]).collect();
            u.dedup();
            u
        }).collect();
        let alternative = (0..units.len() as u32).filter(|&v| {
            let inc = unit_lists.iter().any(|m| m.contains(&v));
            let skip = unit_lists.iter().any(|m| !m.contains(&v) && m.first().is_some_and(|&x| x < v) && m.last().is_some_and(|&x| x > v));
            inc && skip
        }).collect();
        let starts: rustc_hash::FxHashSet<i64> = ms.iter().flat_map(|(_, ex)| ex.iter().map(|&(a, _)| a)).collect();
        let ends: rustc_hash::FxHashSet<i64> = ms.iter().flat_map(|(_, ex)| ex.iter().map(|&(_, b)| b)).collect();
        let exon_start = units.iter().map(|&((lo, _), _)| starts.contains(&lo)).collect();
        let exon_end = units.iter().map(|&((_, hi), _)| ends.contains(&hi)).collect();
        Self { nodes, models, unit_of, units, exon_start, exon_end, alternative }
    }

    /// Whether two alternative units (u before v) are paired: units of different runs always; units of one run only
    /// when they are two exons of their own (an exon ends at u, an exon starts at v) that a model joins by retaining
    /// what lies between.
    fn paired(&self, u: u32, v: u32) -> bool {
        self.units[u as usize].1 != self.units[v as usize].1 || (self.exon_end[u as usize] && self.exon_start[v as usize])
    }
}

impl ReadHypergraph {
    /// Builds the hypergraph of the reported models `reported` and the reads on them.
    pub fn build(db: &PostDB, reported: &[i64]) -> Self {
        let edge_ends: FxHashMap<i64, (i64, i64)> = db.edges.iter().map(|e| (e.edge_id, (e.v1, e.v2))).collect();
        let edge_strand: FxHashMap<i64, bool> = db.edges.iter().map(|e| (e.edge_id, e.strand == "+")).collect();
        // exons of each reported model, oriented
        let mut by_gene: FxHashMap<i64, Vec<(i64, Vec<(i64, i64)>)>> = FxHashMap::default();
        let mut plus_of: FxHashMap<i64, bool> = FxHashMap::default();
        for &tid in reported {
            let Some(&ti) = db.tx_map.get(&tid) else { continue };
            let t = &db.transcripts[ti];
            let Some(&plus) = edge_strand.get(&t.start_exon) else { continue };
            let Some(pos) = crate::quant_em::model_positions(db, t, &edge_ends) else { continue };
            let mut exons: Vec<(i64, i64)> = pos.chunks(2).map(|p| oriented(p[0], p[1], plus)).collect();
            exons.sort_unstable();
            by_gene.entry(t.gene_id).or_default().push((tid, exons));
            plus_of.insert(tid, plus);
        }
        let genes: FxHashMap<i64, GeneGraph> = by_gene.into_iter().map(|(g, mut ms)| {
            ms.sort_by_key(|m| m.0);
            (g, GeneGraph::new(&ms))
        }).collect();
        // reads: the part of their model's segment list between their ends
        let datasets = db.datasets.clone();
        let ds_of: FxHashMap<&str, usize> = datasets.iter().enumerate().map(|(i, d)| (d.as_str(), i)).collect();
        let mut classes: FxHashMap<(i64, Vec<u32>), Vec<f64>> = FxHashMap::default();
        for r in &db.observed {
            let (Some(&plus), Some(&di)) = (plus_of.get(&r.transcript_id), ds_of.get(&*r.dataset)) else { continue };
            let Some(gg) = genes.get(&r.gene_id) else { continue };
            let Some(path) = gg.models.get(&r.transcript_id) else { continue };
            let (x5, x3) = oriented(r.start_pos, r.end_pos, plus);
            let nodes: Vec<u32> = path.iter().copied()
                .filter(|&n| { let (lo, hi) = gg.nodes[n as usize]; hi >= x5 && lo <= x3 }).collect();
            if nodes.is_empty() { continue; }
            let w = classes.entry((r.gene_id, nodes)).or_insert_with(|| vec![0.0; datasets.len()]);
            w[di] += 1.0;
        }
        let mut hyperedges: Vec<Hyperedge> = classes.into_iter().map(|((gene, nodes), weights)| Hyperedge { gene, nodes, weights }).collect();
        hyperedges.sort_by(|a, b| a.gene.cmp(&b.gene).then(a.nodes.cmp(&b.nodes)));
        Self { genes, hyperedges, datasets }
    }

    pub fn n_nodes(&self) -> usize { self.genes.values().map(|g| g.nodes.len()).sum() }

    /// Coordination of pairs of alternative units on different exons: per gene and pair, the reads that span both,
    /// by inclusion. Pairs need `min_reads` such reads and each unit seen both included and skipped among them.
    pub fn coordination(&self, min_reads: f64) -> Vec<Coordination> {
        let mut by_gene: FxHashMap<i64, Vec<&Hyperedge>> = FxHashMap::default();
        for h in &self.hyperedges { by_gene.entry(h.gene).or_default().push(h); }
        let mut out = Vec::new();
        let mut gene_ids: Vec<&i64> = self.genes.keys().collect();
        gene_ids.sort();
        for g in gene_ids {
            let gg = &self.genes[g];
            let alt = &gg.alternative;
            if alt.len() < 2 { continue; }
            let Some(hs) = by_gene.get(g) else { continue };
            // per hyperedge: first and last unit, inclusion of each alternative unit, reads
            let reads: Vec<(u32, u32, Vec<bool>, f64)> = hs.iter().map(|h| {
                let units: Vec<u32> = h.nodes.iter().map(|&n| gg.unit_of[n as usize]).collect();
                let inc = alt.iter().map(|u| units.contains(u)).collect();
                (units[0], *units.last().unwrap(), inc, h.weights.iter().sum())
            }).collect();
            for (ai, &u) in alt.iter().enumerate() {
                for (bi, &v) in alt.iter().enumerate().skip(ai + 1) {
                    if !gg.paired(u, v) { continue; }
                    let mut t = [[0.0f64; 2]; 2];
                    for (first, last, inc, w) in &reads {
                        if *first >= u || *last <= v { continue; }
                        t[inc[ai] as usize][inc[bi] as usize] += w;
                    }
                    let n: f64 = t.iter().flatten().sum();
                    let (u_in, u_out) = (t[1][0] + t[1][1], t[0][0] + t[0][1]);
                    let (v_in, v_out) = (t[0][1] + t[1][1], t[0][0] + t[1][0]);
                    if n < min_reads || u_in == 0.0 || u_out == 0.0 || v_in == 0.0 || v_out == 0.0 { continue; }
                    let log_or = ((t[1][1] + 0.5) * (t[0][0] + 0.5) / ((t[1][0] + 0.5) * (t[0][1] + 0.5))).ln();
                    let c = |x: f64| x.round() as u64;
                    let p = fisher_two_sided([[c(t[0][0]), c(t[0][1])], [c(t[1][0]), c(t[1][1])]]);
                    out.push(Coordination { gene: *g, u_span: gg.units[u as usize].0, v_span: gg.units[v as usize].0,
                                            table: t, log_or, p, q: 1.0 });
                }
            }
        }
        // Benjamini-Hochberg across all tested pairs
        let mut idx: Vec<usize> = (0..out.len()).collect();
        idx.sort_by(|&a, &b| out[a].p.partial_cmp(&out[b].p).unwrap());
        let m = out.len() as f64;
        let mut prev = 1.0f64;
        for (rank, &i) in idx.iter().enumerate().rev() {
            let q = (out[i].p * m / (rank as f64 + 1.0)).min(prev);
            out[i].q = q;
            prev = q;
        }
        out
    }
}

/// One pair of alternative units: the reads spanning both by inclusion (table[a included][b included]).
pub struct Coordination {
    pub gene: i64,
    pub u_span: (i64, i64),
    pub v_span: (i64, i64),
    pub table: [[f64; 2]; 2],
    pub log_or: f64,
    pub p: f64,
    pub q: f64,
}

fn ln_factorial(n: u64) -> f64 {
    // exact sum for small n, Stirling series beyond
    if n < 256 { (1..=n).map(|k| (k as f64).ln()).sum() }
    else {
        let x = n as f64;
        x * x.ln() - x + 0.5 * (2.0 * std::f64::consts::PI * x).ln() + 1.0 / (12.0 * x) - 1.0 / (360.0 * x.powi(3))
    }
}

/// Two-sided Fisher's exact test of a 2x2 table: the probability of the tables with the observed margins that are no
/// more likely than the observed one (relative tolerance 1e-7). Each tail is summed outward from where it starts,
/// until a term falls below 1e-17 of the sum.
pub fn fisher_two_sided(t: [[u64; 2]; 2]) -> f64 {
    let (a, b, c, d) = (t[0][0], t[0][1], t[1][0], t[1][1]);
    let (r1, r2, c1) = (a + b, c + d, a + c);
    let n = r1 + r2;
    if n == 0 { return 1.0; }
    let lo = c1.saturating_sub(r2);
    let hi = r1.min(c1);
    let base = ln_factorial(r1) + ln_factorial(r2) + ln_factorial(c1) + ln_factorial(n - c1) - ln_factorial(n);
    let lp = |x: u64| base - ln_factorial(x) - ln_factorial(r1 - x) - ln_factorial(c1 - x) - ln_factorial(r2 + x - c1);
    let p0 = lp(a);
    let eps = 1e-7;
    let mode = ((((c1 + 1) as f64) * ((r1 + 1) as f64) / ((n + 2) as f64)).floor() as u64).clamp(lo, hi);
    if lp(mode) <= p0 + eps { return 1.0; }
    // sum of exp(lp(x) - p0) from `start` outward (away from the mode)
    let tail = |start: u64, up: bool| -> f64 {
        let (mut s, mut x) = (0.0, start);
        loop {
            let v = (lp(x) - p0).exp();
            s += v;
            if v < 1e-17 * s { break; }
            if up { if x == hi { break; } x += 1; } else { if x == lo { break; } x -= 1; }
        }
        s
    };
    let mut sum;
    if a > mode {
        sum = tail(a, true);
        // other side: the largest x below the mode with lp(x) <= p0 (lp increases towards the mode)
        if lp(lo) <= p0 + eps {
            let (mut l, mut h) = (lo, mode);
            while l < h { let m = (l + h + 1) / 2; if lp(m) <= p0 + eps { l = m } else { h = m - 1 } }
            sum += tail(l, false);
        }
    } else {
        sum = tail(a, false);
        // other side: the smallest x above the mode with lp(x) <= p0 (lp decreases away from the mode)
        if lp(hi) <= p0 + eps {
            let (mut l, mut h) = (mode, hi);
            while l < h { let m = (l + h) / 2; if lp(m) <= p0 + eps { h = m } else { l = m + 1 } }
            sum += tail(l, true);
        }
    }
    (sum * p0.exp()).min(1.0)
}

/// Writes <prefix>_isoforge_hypergraph.tsv (hyperedges: gene, segments, weights) and
/// <prefix>_isoforge_coordination.tsv (pairs of alternative units with their table and test).
pub fn write_hypergraph(db: &PostDB, hg: &ReadHypergraph, prefix: &str, coordination: &[Coordination]) -> IsoForgeResult<(String, String)> {
    let gname = |g: i64| db.annot.gene_name_map.get(&g).cloned().unwrap_or_else(|| g.to_string());
    let chrom_of: FxHashMap<i64, (String, String)> = db.genes.iter().map(|g| (g.gene_id, (g.chromosome.clone(), g.strand.clone()))).collect();
    let hpath = format!("{}_isoforge_hypergraph.tsv", prefix);
    let mut h = BufWriter::new(File::create(&hpath)?);
    write!(h, "gene\tgene_name\tsegments")?;
    for d in &hg.datasets { write!(h, "\t{}", d)?; }
    writeln!(h)?;
    for e in &hg.hyperedges {
        let nodes: Vec<String> = e.nodes.iter().map(|n| n.to_string()).collect();
        write!(h, "{}\t{}\t{}", e.gene, gname(e.gene), nodes.join(","))?;
        for w in &e.weights { write!(h, "\t{}", w)?; }
        writeln!(h)?;
    }
    let cpath = format!("{}_isoforge_coordination.tsv", prefix);
    let mut c = BufWriter::new(File::create(&cpath)?);
    writeln!(c, "gene\tgene_name\tchromosome\tstrand\tunit_a\tunit_b\tboth_skipped\ta_only\tb_only\tboth_included\tlog_odds_ratio\tp_value\tq_value")?;
    let span = |s: (i64, i64), strand: &str| -> String {
        let (a, b) = if strand == "+" { (s.0, s.1) } else { (-s.1, -s.0) };
        format!("{}-{}", a, b)
    };
    for x in coordination {
        let (chrom, strand) = chrom_of.get(&x.gene).cloned().unwrap_or_default();
        writeln!(c, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{:.3e}\t{:.3e}", x.gene, gname(x.gene), chrom, strand,
                 span(x.u_span, &strand), span(x.v_span, &strand), x.table[0][0], x.table[1][0], x.table[0][1], x.table[1][1],
                 x.log_or, x.p, x.q)?;
    }
    Ok((hpath, cpath))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fisher_two_sided_p_values() {
        // reference P values of the two-sided test
        let close = |t: [[u64; 2]; 2], want: f64| {
            let got = fisher_two_sided(t);
            assert!((got - want).abs() <= 1e-9 + 1e-6 * want, "{:?}: {} vs {}", t, got, want);
        };
        close([[3, 1], [1, 3]], 0.48571428571428565);
        close([[100, 200], [300, 50]], 2.3546829814954165e-44);
        close([[0, 10], [10, 0]], 1.082508822446903e-05);
        close([[5, 5], [5, 5]], 1.0);
        close([[7, 3], [2, 9]], 0.029973122852379817);
        close([[0, 0], [5, 7]], 1.0);
        assert!(fisher_two_sided([[40000, 20000], [19000, 41000]]) < 1e-300);
    }

    #[test]
    fn units_merge_touching_segments_and_pair_only_across_exons() {
        // exons 1-3 of four models, + strand: exon 2 is a cassette; exon 3 has an alternative 3' splice site (a piece
        // 500-509 used by two models), so units: [100-199], [300-349], [500-509], [510-599]
        let ms = vec![
            (1, vec![(100, 199), (300, 349), (500, 599)]),
            (2, vec![(100, 199), (510, 599)]),
            (3, vec![(100, 199), (300, 349), (510, 599)]),
            (4, vec![(100, 199), (500, 599)]),
        ];
        let g = GeneGraph::new(&ms);
        assert_eq!(g.units.iter().map(|u| u.0).collect::<Vec<_>>(), vec![(100, 199), (300, 349), (500, 509), (510, 599)]);
        // exons: 0, 1, 2 (500-509 touches 510-599)
        assert_eq!(g.units.iter().map(|u| u.1).collect::<Vec<_>>(), vec![0, 1, 2, 2]);
        // alternative: the cassette (unit 1) and the 3' splice-site piece (unit 2)
        assert_eq!(g.alternative, vec![1, 2]);
        // the cassette pairs with the piece of exon 3; the two pieces of exon 3 are one exon
        assert!(g.paired(1, 2));
        assert!(!g.paired(2, 3));
    }

    #[test]
    fn exons_joined_by_a_retained_intron_are_still_paired() {
        // two mutually exclusive exons (300-349, 500-549) and a model that retains the intron between them: one run
        // of touching units [300-349], [350-499], [500-549]
        let ms = vec![
            (1, vec![(100, 199), (300, 349), (700, 799)]),
            (2, vec![(100, 199), (500, 549), (700, 799)]),
            (3, vec![(100, 199), (300, 549), (700, 799)]),
        ];
        let g = GeneGraph::new(&ms);
        assert_eq!(g.units.iter().map(|u| u.0).collect::<Vec<_>>(), vec![(100, 199), (300, 349), (350, 499), (500, 549), (700, 799)]);
        assert_eq!(g.units.iter().map(|u| u.1).collect::<Vec<_>>(), vec![0, 1, 1, 1, 2]);
        assert_eq!(g.alternative, vec![1, 2, 3]);
        // the two exons are paired; neither pairs with the retained intron beside it
        assert!(g.paired(1, 3));
        assert!(!g.paired(1, 2) && !g.paired(2, 3));
        // reads: 40 with the first exon, 50 with the second, 3 retaining the intron, 2 skipping both
        let edge = |nodes: Vec<u32>, w: f64| Hyperedge { gene: 7, nodes, weights: vec![w] };
        let h = ReadHypergraph {
            genes: [(7, g)].into_iter().collect(),
            hyperedges: vec![edge(vec![0, 1, 4], 40.0), edge(vec![0, 3, 4], 50.0), edge(vec![0, 1, 2, 3, 4], 3.0), edge(vec![0, 4], 2.0)],
            datasets: vec!["d".into()],
        };
        let c = h.coordination(20.0);
        let pair = c.iter().find(|x| x.u_span == (300, 349) && x.v_span == (500, 549)).expect("the two exons are tested");
        assert_eq!(pair.table, [[2.0, 50.0], [40.0, 3.0]]);
        assert!(pair.log_or < -4.0 && pair.p < 1e-10);
        // an exon and the intron retained beside it are never a pair
        assert!(!c.iter().any(|x| x.u_span == (350, 499) || x.v_span == (350, 499)));
    }
}
