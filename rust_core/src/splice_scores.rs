//! Splice-site scores: one score per site, read from a TSV table (`chrom pos
//! strand D|A score`, optionally gzipped). The table places a site at the base
//! before the intron for the intron's left boundary and at the intron's last
//! base for its right boundary; `score` takes an IsoForge splice position
//! (1-based exon boundary) and its role.

use crate::errors::{IsoForgeError, IsoForgeResult};
use crate::types::Strand;
use rustc_hash::FxHashMap;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

/// Sorted positions and their scores for one (strand, role).
#[derive(Debug, Default)]
struct Track {
    positions: Vec<u32>,
    scores: Vec<i8>,
}

/// Scores per chromosome, four tracks each: (+, donor), (+, acceptor), (-, donor), (-, acceptor).
#[derive(Debug, Default)]
pub struct SpliceScores {
    chroms: FxHashMap<String, [Track; 4]>,
    n_sites: usize,
}

fn track_index(strand: Strand, donor: bool) -> usize {
    (if strand == Strand::Minus { 2 } else { 0 }) + usize::from(!donor)
}

impl SpliceScores {
    pub fn load(path: &Path) -> IsoForgeResult<Self> {
        let file = std::fs::File::open(path)
            .map_err(|e| IsoForgeError::ConfigError(format!("splice scores {}: {e}", path.display())))?;
        let reader: Box<dyn Read> = if path.extension().is_some_and(|x| x == "gz") {
            Box::new(flate2::read::MultiGzDecoder::new(file))
        } else {
            Box::new(file)
        };
        Self::from_reader(BufReader::with_capacity(1 << 20, reader))
    }

    pub fn from_reader<R: BufRead>(reader: R) -> IsoForgeResult<Self> {
        let mut out = SpliceScores::default();
        let mut unsorted = false;
        for line in reader.lines() {
            let line = line?;
            let mut f = line.split('\t');
            let (Some(chrom), Some(pos), Some(strand), Some(kind), Some(score)) = (f.next(), f.next(), f.next(), f.next(), f.next())
            else {
                continue;
            };
            let (Ok(pos), Ok(score)) = (pos.parse::<u32>(), score.trim().parse::<i32>()) else { continue };
            let strand = match strand { "+" => Strand::Plus, "-" => Strand::Minus, _ => continue };
            let donor = match kind { "D" => true, "A" => false, _ => continue };
            if !out.chroms.contains_key(chrom) {
                out.chroms.insert(chrom.to_string(), Default::default());
            }
            let t = &mut out.chroms.get_mut(chrom).expect("inserted")[track_index(strand, donor)];
            unsorted |= t.positions.last().is_some_and(|&p| p > pos);
            t.positions.push(pos);
            t.scores.push(score.clamp(i8::MIN as i32, i8::MAX as i32) as i8);
            out.n_sites += 1;
        }
        if unsorted {
            for tracks in out.chroms.values_mut() {
                for t in tracks.iter_mut() {
                    let mut idx: Vec<usize> = (0..t.positions.len()).collect();
                    idx.sort_unstable_by_key(|&i| t.positions[i]);
                    t.positions = idx.iter().map(|&i| t.positions[i]).collect();
                    t.scores = idx.iter().map(|&i| t.scores[i]).collect();
                }
            }
        }
        Ok(out)
    }

    pub fn len(&self) -> usize {
        self.n_sites
    }

    pub fn is_empty(&self) -> bool {
        self.n_sites == 0
    }

    /// Whether the scores cover this chromosome.
    pub fn has_chrom(&self, chrom: &str) -> bool {
        self.chroms.contains_key(chrom)
    }

    /// Score of the splice site at IsoForge position `pos` (the exon's last base
    /// before a donor on +, the exon's first base after an acceptor on +, and the
    /// mirror on -), or None when the table has no site there.
    pub fn score(&self, chrom: &str, strand: Strand, donor: bool, pos: i64) -> Option<i8> {
        let tracks = self.chroms.get(chrom)?;
        // Left boundary of the intron (a donor on +, an acceptor on -) is the exon's
        // last base; the right boundary is one base before the next exon.
        let left = (strand == Strand::Plus) == donor;
        let key = if left { pos } else { pos - 1 };
        let key = u32::try_from(key).ok()?;
        let t = &tracks[track_index(strand, donor)];
        t.positions.binary_search(&key).ok().map(|i| t.scores[i])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const TSV: &str = "chr1\t100\t+\tD\t11\nchr1\t199\t+\tA\t9\nchr1\t299\t-\tD\t12\nchr1\t200\t-\tA\t-3\n";

    #[test]
    fn scores_line_up_with_isoforge_positions_on_both_strands() {
        let s = SpliceScores::from_reader(TSV.as_bytes()).unwrap();
        assert_eq!(s.len(), 4);
        // + strand intron after exon end 100, next exon from 200
        assert_eq!(s.score("chr1", Strand::Plus, true, 100), Some(11));
        assert_eq!(s.score("chr1", Strand::Plus, false, 200), Some(9));
        // - strand intron between exon ends 200 (left) and 300 (right): donor at 300, acceptor at 200
        assert_eq!(s.score("chr1", Strand::Minus, true, 300), Some(12));
        assert_eq!(s.score("chr1", Strand::Minus, false, 200), Some(-3));
        assert_eq!(s.score("chr1", Strand::Plus, true, 101), None);
        assert_eq!(s.score("chr2", Strand::Plus, true, 100), None);
    }
}
