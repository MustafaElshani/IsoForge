//! Bounded sequence-based repair of annotated internal microexons.
//! Inclusion must fit the local read better than the observed skipped alignment.
use crate::data_structures::SharedRefs;
use crate::genome_ref::GenomeRef;
use crate::types::Strand;
use rust_htslib::bam;
use std::collections::{BTreeMap, BTreeSet};

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct Candidate {
    left: i64,
    start: i64,
    end: i64,
    right: i64,
}

#[derive(Default, Debug)]
pub struct MicroexonIndex {
    // All coordinates 0-based half-open; left/right bound the combined skip.
    loci: BTreeMap<(String, Strand), BTreeMap<i64, Vec<Candidate>>>,
}

impl MicroexonIndex {
    pub fn from_reference(refs: &SharedRefs) -> Self {
        let mut seen = BTreeSet::new();
        let mut index = Self::default();
        for tx in refs.transcript_dict.values() {
            // the transcript's hyperedge: vertices 5'->3', two per exon; intron j runs from vertex 2j+1 to 2j+2
            let Some(row) = refs.hypergraph.get(tx.transcript_id).and_then(|h| h.row.as_ref()) else {
                continue;
            };
            let p = &row.positions;
            for j in 0..row.n_exons().saturating_sub(2) {
                // intron j, internal exon j+1, intron j+1
                let (a, b, c) = ((p[2 * j + 1], p[2 * j + 2]), (p[2 * j + 2], p[2 * j + 3]), (p[2 * j + 3], p[2 * j + 4]));
                let bounds = [a.0, a.1, c.0, c.1];
                let candidate = Candidate {
                    left: *bounds.iter().min().unwrap(),
                    start: b.0.min(b.1) - 1,
                    end: b.0.max(b.1),
                    right: *bounds.iter().max().unwrap() - 1,
                };
                if candidate.end - candidate.start > 30
                    || candidate.start <= candidate.left
                    || candidate.end >= candidate.right
                {
                    continue;
                }
                if !seen.insert((tx.chromosome.clone(), row.strand, candidate)) {
                    continue;
                }
                index
                    .loci
                    .entry((tx.chromosome.clone(), row.strand))
                    .or_default()
                    .entry(candidate.left)
                    .or_default()
                    .push(candidate);
            }
        }
        index
    }

    pub fn repair(
        &self,
        record: &bam::Record,
        chrom: &str,
        strand: Strand,
        genome: &GenomeRef,
        introns: &[i64],
    ) -> (Vec<i64>, usize) {
        let Some(locus) = self.loci.get(&(chrom.to_string(), strand)) else {
            return (introns.to_vec(), 0);
        };
        let mut output = Vec::with_capacity(introns.len());
        let sequence = record.seq().as_bytes();
        let mut repaired = 0;
        for pair in introns.chunks_exact(2) {
            let left = pair[0] - 1;
            let right = pair[1];
            let lo = (left - 24).max(record.pos());
            let hi = (right + 24).min(record.cigar().end_pos());
            let mut scored = Vec::new();
            let mut best_skipped = usize::MAX;
            let mut query_length = 0;
            // Candidates whose bounds lie within 15 bp of the intron's; with more than 32 the intron is left as it is.
            let local: Vec<_> = locus
                .range(left - 15..=left + 15)
                .flat_map(|(_, cs)| cs.iter())
                .filter(|c| (c.right - right).abs() <= 15)
                .take(33)
                .collect();
            if local.len() > 32 {
                output.extend(pair);
                continue;
            }
            if left - lo >= 12 && hi - right >= 12 {
                if let (Some(qlo), Some(qhi)) =
                    (query_boundary(record, lo), query_boundary(record, hi))
                {
                    if qhi > qlo && qhi - qlo <= 128 {
                        let query = &sequence[qlo..qhi];
                        query_length = query.len();
                        if let (Ok(mut skip), Ok(tail)) = (
                            genome.fetch_seq(chrom, lo, left),
                            genome.fetch_seq(chrom, right, hi),
                        ) {
                            skip.extend(tail);
                            best_skipped = edit_distance(query, &skip);
                            for c in local {
                                let c = *c;
                                if (c.right - right).abs() > 15
                                    || c.left - lo < 8
                                    || hi - c.right < 8
                                {
                                    continue;
                                }
                                if let (Ok(mut included), Ok(exon), Ok(tail)) = (
                                    genome.fetch_seq(chrom, lo, c.left),
                                    genome.fetch_seq(chrom, c.start, c.end),
                                    genome.fetch_seq(chrom, c.right, hi),
                                ) {
                                    // Inclusion is compared with skipping at the
                                    // same corrected flanks.
                                    let mut adjusted_skip = included.clone();
                                    adjusted_skip.extend(&tail);
                                    best_skipped =
                                        best_skipped.min(edit_distance(query, &adjusted_skip));
                                    included.extend(exon);
                                    included.extend(tail);
                                    let score = edit_distance(query, &included);
                                    // All competitors are kept until uniqueness
                                    // is assessed.
                                    scored.push((score, c));
                                }
                            }
                        }
                    }
                }
            }
            scored.sort_unstable();
            if let Some(&(score, c)) = scored.first() {
                if score + 2 <= best_skipped
                    && score * 4 <= query_length
                    && scored.get(1).is_none_or(|x| x.0 >= score + 2)
                {
                    output.extend([c.left + 1, c.start, c.end + 1, c.right]);
                    repaired += 1;
                    continue;
                }
            }
            output.extend(pair);
        }
        // The local proposal cannot cross another observed intron.
        if output.windows(2).any(|w| w[0] >= w[1]) {
            (introns.to_vec(), 0)
        } else {
            (output, repaired)
        }
    }
}

pub(crate) fn query_boundary(record: &bam::Record, target: i64) -> Option<usize> {
    use bam::record::Cigar::*;
    let mut r = record.pos();
    let mut q = 0usize;
    for op in record.cigar().iter() {
        match *op {
            Match(n) | Equal(n) | Diff(n) => {
                if target >= r && target <= r + n as i64 {
                    return Some(q + (target - r) as usize);
                }
                r += n as i64;
                q += n as usize;
            }
            Ins(n) | SoftClip(n) => q += n as usize,
            Del(n) | RefSkip(n) => r += n as i64,
            _ => (),
        }
    }
    None
}

pub(crate) fn edit_distance(a: &[u8], b: &[u8]) -> usize {
    let mut previous: Vec<usize> = (0..=b.len()).collect();
    let mut current = vec![0; b.len() + 1];
    for (i, x) in a.iter().enumerate() {
        current[0] = i + 1;
        for (j, y) in b.iter().enumerate() {
            current[j + 1] = (previous[j + 1] + 1)
                .min(current[j] + 1)
                .min(previous[j] + usize::from(!x.eq_ignore_ascii_case(y)));
        }
        std::mem::swap(&mut previous, &mut current);
    }
    previous[b.len()]
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn sequence_evidence_distinguishes_inclusion_from_skipping() {
        let a = b"ACGTACGTACGT";
        let b = b"TGACTGACTGAC";
        let mut inclusion = a.to_vec();
        inclusion.extend(b"GGCC");
        inclusion.extend(b);
        let mut skipping = a.to_vec();
        skipping.extend(b);
        assert_eq!(edit_distance(&inclusion, &inclusion), 0);
        assert_eq!(edit_distance(&skipping, &inclusion), 4);
        assert_eq!(edit_distance(&inclusion, &skipping), 4);
    }
}

#[cfg(test)]
mod alignment_tests {
    use super::*;
    use bam::record::{Cigar, CigarString};
    #[test]
    fn repairs_sequence_supported_four_base_exon_but_preserves_true_skipping() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("genome.fa");
        let mut bases = vec![b'C'; 400];
        bases[200..204].copy_from_slice(b"GGCC");
        std::fs::write(
            &fasta,
            format!(">chr1\n{}\n", String::from_utf8(bases).unwrap()),
        )
        .unwrap();
        std::fs::write(dir.path().join("genome.fa.fai"), "chr1\t400\t6\t400\t401\n").unwrap();
        let genome = GenomeRef::open(fasta).unwrap();
        for strand in [Strand::Plus, Strand::Minus] {
            let mut index = MicroexonIndex::default();
            index
                .loci
                .entry(("chr1".into(), strand))
                .or_default()
                .entry(100)
                .or_default()
                .push(Candidate {
                    left: 100,
                    start: 200,
                    end: 204,
                    right: 300,
                });
            for included in [true, false] {
                let mut record = bam::Record::new();
                record.set_pos(0);
                let mut seq = vec![b'C'; 100];
                let mut cigar = vec![Cigar::Match(100)];
                if included {
                    seq.extend(b"GGCC");
                    cigar.push(Cigar::Ins(4));
                }
                seq.extend(vec![b'C'; 100]);
                cigar.extend([Cigar::RefSkip(200), Cigar::Match(100)]);
                record.set(
                    b"read",
                    Some(&CigarString(cigar)),
                    &seq,
                    &vec![30; seq.len()],
                );
                record.set_flags(if strand == Strand::Minus { 16 } else { 0 });
                let (introns, n) = index.repair(&record, "chr1", strand, &genome, &[101, 300]);
                if included {
                    assert_eq!(introns, vec![101, 200, 205, 300]);
                    assert_eq!(n, 1);
                } else {
                    assert_eq!(introns, vec![101, 300]);
                    assert_eq!(n, 0);
                }
            }
        }
    }
}

#[cfg(test)]
mod boundary_controls {
    use super::*;
    use bam::record::{Cigar, CigarString};
    #[test]
    fn boundary_improvement_alone_does_not_establish_an_exon() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("genome.fa");
        let mut sequence = vec![b'C'; 400];
        sequence[100..104].copy_from_slice(b"GGGG");
        std::fs::write(
            &fasta,
            format!(">chr1\n{}\n", String::from_utf8(sequence).unwrap()),
        )
        .unwrap();
        std::fs::write(dir.path().join("genome.fa.fai"), "chr1\t400\t6\t400\t401\n").unwrap();
        let genome = GenomeRef::open(fasta).unwrap();
        for strand in [Strand::Plus, Strand::Minus] {
            let mut index = MicroexonIndex::default();
            index
                .loci
                .entry(("chr1".into(), strand))
                .or_default()
                .entry(100)
                .or_default()
                .push(Candidate {
                    left: 100,
                    start: 200,
                    end: 203,
                    right: 300,
                });
            let mut record = bam::Record::new();
            record.set_pos(0);
            record.set(
                b"boundary_error",
                Some(&CigarString(vec![
                    Cigar::Match(102),
                    Cigar::Del(2),
                    Cigar::RefSkip(196),
                    Cigar::Match(100),
                ])),
                &vec![b'C'; 202],
                &vec![30; 202],
            );
            let (introns, n) = index.repair(&record, "chr1", strand, &genome, &[105, 300]);
            assert_eq!(introns, vec![105, 300]);
            assert_eq!(n, 0);
        }
    }
}
