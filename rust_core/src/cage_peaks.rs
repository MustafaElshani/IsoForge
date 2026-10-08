// CAGE peak index for 5′ end validation: per-(chrom, strand) interval trees
// built from a BED or BED.gz file (columns chrom, start, end and optional
// strand).

use std::collections::HashMap;
use std::io::{BufRead, BufReader};
use std::path::Path;

use coitrees::{COITree, Interval, IntervalTree};

use crate::errors::*;
use crate::types::Strand;

/// CAGE peak interval tree index for overlap queries.
pub struct CagePeakIndex {
    /// Per-(chrom, strand) interval trees
    trees: HashMap<(String, Strand), COITree<(), u32>>,
    /// Total number of peaks loaded
    pub n_peaks: usize,
}

impl CagePeakIndex {
    /// Load CAGE peaks from a BED file (plain or gzipped).
    ///
    /// Accepts both `.bed` and `.bed.gz` formats. If no strand column (col 5)
    /// is present, peaks are indexed on BOTH strands.
    pub fn from_bed(path: &Path) -> IsoForgeResult<Self> {
        let file = std::fs::File::open(path)
            .map_err(|e| IsoForgeError::IoError(e))?;

        let is_gzipped = path
            .to_str()
            .map(|s| s.ends_with(".gz"))
            .unwrap_or(false);

        let reader: Box<dyn BufRead> = if is_gzipped {
            use flate2::read::GzDecoder;
            Box::new(BufReader::new(GzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };

        let mut intervals: HashMap<(String, Strand), Vec<Interval<()>>> = HashMap::new();
        let mut n_peaks = 0usize;

        for line_result in reader.lines() {
            let line = line_result.map_err(|e| IsoForgeError::IoError(e))?;
            let line = line.trim();

            // Skip headers, comments, track lines, browser lines
            if line.is_empty()
                || line.starts_with('#')
                || line.starts_with("track")
                || line.starts_with("browser")
            {
                continue;
            }

            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 3 {
                continue;
            }

            let chrom = fields[0].to_string();
            let start: i32 = match fields[1].parse() {
                Ok(v) => v,
                Err(_) => continue,
            };
            let end: i32 = match fields[2].parse() {
                Ok(v) => v,
                Err(_) => continue,
            };

            // BED is 0-based half-open: [start, end)
            // coitrees uses inclusive intervals, so use [start, end-1]
            let interval = Interval::new(start, end.saturating_sub(1), ());

            // Parse strand if available (column 5)
            let strands: Vec<Strand> = if fields.len() > 5 {
                match fields[5] {
                    "+" => vec![Strand::Plus],
                    "-" => vec![Strand::Minus],
                    _ => vec![Strand::Plus, Strand::Minus], // unstranded → both
                }
            } else {
                vec![Strand::Plus, Strand::Minus] // no strand → both
            };

            for strand in strands {
                intervals
                    .entry((chrom.clone(), strand))
                    .or_default()
                    .push(interval.clone());
            }
            n_peaks += 1;
        }

        // Build interval trees
        let mut trees = HashMap::new();
        for (key, ivs) in &intervals {
            trees.insert(key.clone(), COITree::new(&*ivs));
        }

        eprintln!(
            "[IsoForge] Loaded {} CAGE peaks from {:?} ({} chrom/strand groups)",
            n_peaks,
            path,
            trees.len()
        );

        Ok(CagePeakIndex { trees, n_peaks })
    }

    /// Check if a genomic position overlaps any CAGE peak.
    ///
    /// Position is 1-based (genomic coordinate). Internally converts to
    /// 0-based for the BED-origin interval tree.
    pub fn overlaps(&self, chrom: &str, strand: Strand, position: i64) -> bool {
        let key = (chrom.to_string(), strand);
        if let Some(tree) = self.trees.get(&key) {
            let pos_0based = (position - 1) as i32; // convert 1-based → 0-based
            tree.query_count(pos_0based, pos_0based) > 0
        } else {
            false
        }
    }

    /// Check if a position is within `window` bp of any CAGE peak.
    pub fn overlaps_with_window(
        &self,
        chrom: &str,
        strand: Strand,
        position: i64,
        window: i64,
    ) -> bool {
        let key = (chrom.to_string(), strand);
        if let Some(tree) = self.trees.get(&key) {
            let pos_0based = (position - 1) as i32;
            let start = pos_0based - window as i32;
            let end = pos_0based + window as i32;
            tree.query_count(start, end) > 0
        } else {
            false
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_cage_peaks_from_bed() {
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, "chr1\t100\t200\tpeak1\t0\t+").unwrap();
        writeln!(tmp, "chr1\t500\t600\tpeak2\t0\t-").unwrap();
        writeln!(tmp, "chrX\t1000\t1100\tpeak3\t0\t+").unwrap();
        tmp.flush().unwrap();

        let index = CagePeakIndex::from_bed(tmp.path()).unwrap();
        assert_eq!(index.n_peaks, 3);

        // chr1:+ position 150 (1-based) → 0-based 149, within [100,199]
        assert!(index.overlaps("chr1", Strand::Plus, 150));
        // chr1:+ position 50 → outside
        assert!(!index.overlaps("chr1", Strand::Plus, 50));
        // chr1:- position 550 (1-based) → 0-based 549, within [500,599]
        assert!(index.overlaps("chr1", Strand::Minus, 550));
        // Wrong strand
        assert!(!index.overlaps("chr1", Strand::Plus, 550));
    }

    #[test]
    fn test_cage_peaks_with_window() {
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(tmp, "chr1\t100\t200\tpeak1\t0\t+").unwrap();
        tmp.flush().unwrap();

        let index = CagePeakIndex::from_bed(tmp.path()).unwrap();

        // Position 95 (1-based) = 0-based 94, not within [100,199]
        assert!(!index.overlaps("chr1", Strand::Plus, 95));
        // But with 10bp window: [84, 104] overlaps [100,199]
        assert!(index.overlaps_with_window("chr1", Strand::Plus, 95, 10));
    }
}
