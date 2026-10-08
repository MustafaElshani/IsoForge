// Indexed genome FASTA access for fraction-A, PAS hexamer scanning, splice
// motif checks and MD tag computation. Each thread uses its own reader.

use rust_htslib::faidx;
use std::collections::HashMap;
use std::path::{Path, PathBuf};

use crate::errors::{IsoForgeError, IsoForgeResult};
use crate::types::{positioned_pas, u_rich_downstream, InternalPrimingResult, Strand, DSE_LABEL, DSE_WINDOW, PAS_SCAN_WINDOW};

/// Genome reference FASTA reader (htslib faidx); `reopen` gives each thread
/// its own file handle.
pub struct GenomeRef {
    path: PathBuf,
    reader: faidx::Reader,
}

// Safety: each GenomeRef owns its own faidx handle.
unsafe impl Send for GenomeRef {}
// Safety: each thread uses its own GenomeRef.
unsafe impl Sync for GenomeRef {}

impl GenomeRef {
    /// Open an indexed FASTA file (.fai must exist alongside).
    pub fn open<P: AsRef<Path>>(fasta_path: P) -> IsoForgeResult<Self> {
        let path = fasta_path.as_ref().to_path_buf();
        let fai_path = format!("{}.fai", path.display());
        if !Path::new(&fai_path).exists() {
            return Err(IsoForgeError::IoError(std::io::Error::new(
                std::io::ErrorKind::NotFound,
                format!(
                    "FASTA index not found: {}. Run `samtools faidx {}`",
                    fai_path,
                    path.display()
                ),
            )));
        }
        let reader = faidx::Reader::from_path(&path)
            .map_err(|e| IsoForgeError::HtsError(format!("Failed to open FASTA: {}", e)))?;
        Ok(Self { path, reader })
    }

    /// Create a new reader handle for the same file (for use in another thread).
    pub fn reopen(&self) -> IsoForgeResult<Self> {
        Self::open(&self.path)
    }

    /// Fetch sequence for a 0-based half-open interval [start, end).
    /// Returns uppercase ASCII.
    pub fn fetch_seq(&self, chrom: &str, start: i64, end: i64) -> IsoForgeResult<Vec<u8>> {
        if start >= end || start < 0 {
            return Ok(Vec::new());
        }
        // faidx::Reader::fetch_seq takes 0-based start and INCLUSIVE end
        let seq = self
            .reader
            .fetch_seq(chrom, start as usize, (end - 1) as usize)
            .map_err(|e| {
                IsoForgeError::HtsError(format!(
                    "Failed to fetch {}:{}-{}: {}",
                    chrom, start, end, e
                ))
            })?;
        // Uppercase
        Ok(seq.iter().map(|b| b.to_ascii_uppercase()).collect())
    }

    /// Fetch sequence as a String.
    pub fn fetch_seq_string(&self, chrom: &str, start: i64, end: i64) -> IsoForgeResult<String> {
        let seq = self.fetch_seq(chrom, start, end)?;
        Ok(String::from_utf8_lossy(&seq).to_string())
    }

    /// Fraction of A in the `range_size` genomic bases after a transcript's 3' end, on the transcript's strand;
    /// `transcript_end` is the 1-based last aligned base in transcript orientation.
    /// Returns the fraction in [0.0, 1.0], or None if the sequence cannot be fetched.
    pub fn compute_fraction_a(
        &self,
        chrom: &str,
        transcript_end: i64,
        strand: Strand,
        range_size: i64,
    ) -> Option<f64> {
        if range_size <= 0 {
            return Some(0.0);
        }

        // Genomic window (1-based)
        let (window_start_1based, window_end_1based) = match strand {
            Strand::Plus => {
                // After the 3' end on + strand
                let start = transcript_end + 1;
                let end = transcript_end + range_size;
                (start, end)
            }
            Strand::Minus => {
                // Before the 3' end on - strand (which is the leftmost position)
                let end = transcript_end - 1;
                let start = transcript_end - range_size;
                (start, end)
            }
        };

        // Clamp to positive coordinates
        let start_0 = (window_start_1based - 1).max(0);
        let end_0 = window_end_1based; // Convert to half-open [start, end)

        // Fetch sequence
        let seq = self.fetch_seq(chrom, start_0, end_0).ok()?;
        if seq.is_empty() {
            return Some(0.0);
        }

        // Count A's (T's on the minus strand)
        let target = match strand {
            Strand::Plus => b'A',
            Strand::Minus => b'T', // reverse complement of A
        };

        let count = seq.iter().filter(|&&b| b == target).count();
        Some(count as f64 / seq.len() as f64)
    }

    /// PAS-aware internal priming assessment. A read is internally primed when its downstream fraction-A exceeds the
    /// threshold and no poly(A) signal is positioned to direct cleavage at its 3' end (see [`positioned_pas`]).
    pub fn assess_internal_priming(
        &self,
        chrom: &str,
        transcript_end: i64,
        strand: Strand,
        range_size: i64,
        frac_a_threshold: f64,
    ) -> InternalPrimingResult {
        // Fraction-A
        let fraction_a = self
            .compute_fraction_a(chrom, transcript_end, strand, range_size)
            .unwrap_or(0.0);

        // PAS hexamers upstream of the 3' end
        let (has_pas, pas_motif, pas_distance) =
            self.scan_pas_upstream(chrom, transcript_end, strand);

        // Internally primed = high fraction-A without a PAS signal
        let is_internal_priming = fraction_a > frac_a_threshold && !has_pas;

        InternalPrimingResult {
            fraction_a,
            has_pas_signal: has_pas,
            pas_motif,
            pas_distance,
            is_internal_priming,
        }
    }

    /// The poly(A) signal upstream (in transcription direction) of a transcript 3' end: the strongest hexamer placed
    /// to direct cleavage there (see [`positioned_pas`]), else a U-rich downstream element.
    /// Returns: (found, motif, distance in bases upstream of the 3' end)
    pub fn scan_pas_upstream(
        &self,
        chrom: &str,
        transcript_end: i64,
        strand: Strand,
    ) -> (bool, Option<String>, Option<i64>) {
        let scan_window = PAS_SCAN_WINDOW;

        // The scan_window bases ending at the 3' end base, in 0-based half-open coordinates;
        // `transcript_end` is the 1-based last aligned base in transcript orientation.
        let (region_start, region_end) = match strand {
            Strand::Plus => ((transcript_end - scan_window).max(0), transcript_end),
            Strand::Minus => (transcript_end - 1, transcript_end - 1 + scan_window),
        };

        // Fetch the genomic sequence
        let seq = match self.fetch_seq(chrom, region_start, region_end) {
            Ok(s) if !s.is_empty() => s,
            _ => return (false, None, None),
        };

        // Minus strand is searched on the reverse complement
        let search_seq = if strand == Strand::Minus {
            reverse_complement(&seq)
        } else {
            seq
        };

        // Only a hexamer positioned to direct cleavage at this end is its signal; failing
        // one, a U-rich downstream element is.
        if let Some((motif, gap)) = positioned_pas(&search_seq) {
            return (true, Some(motif.to_string()), Some(gap));
        }
        if self.downstream_element(chrom, transcript_end, strand) {
            return (true, Some(DSE_LABEL.to_string()), Some(0));
        }
        (false, None, None)
    }

    /// Whether a U-rich downstream element follows the 3' end (see [`u_rich_downstream`]).
    pub fn downstream_element(&self, chrom: &str, transcript_end: i64, strand: Strand) -> bool {
        // The DSE_WINDOW bases after the 3' end base, 0-based half-open.
        let (a, b) = match strand {
            Strand::Plus => (transcript_end, transcript_end + DSE_WINDOW),
            Strand::Minus => ((transcript_end - 1 - DSE_WINDOW).max(0), transcript_end - 1),
        };
        match self.fetch_seq(chrom, a, b) {
            Ok(seq) => u_rich_downstream(&if strand == Strand::Minus { reverse_complement(&seq) } else { seq }),
            Err(_) => false,
        }
    }

    /// Check whether a splice junction has a canonical dinucleotide motif (GT-AG, GC-AG or AT-AC on the + strand).
    /// `intron_start` and `intron_end` are the 1-based first and last intronic bases; `strand` selects the motifs
    /// checked. Returns is_canonical, donor_motif, acceptor_motif and motif_class.
    pub fn check_splice_motif(
        &self,
        chrom: &str,
        intron_start: i64,
        intron_end: i64,
        strand: Strand,
    ) -> IsoForgeResult<SpliceMotifResult> {
        // Fetch donor dinucleotide (first 2 bases of intron)
        let donor_start = intron_start - 1; // 0-based
        let donor = self.fetch_seq(chrom, donor_start, donor_start + 2)?;

        // Acceptor dinucleotide: 0-based (intron_end - 2, intron_end - 1)
        let acceptor_start = intron_end - 2; // 0-based start of acceptor dinucleotide
        let acceptor = self.fetch_seq(chrom, acceptor_start, acceptor_start + 2)?;

        if donor.len() < 2 || acceptor.len() < 2 {
            return Ok(SpliceMotifResult {
                is_canonical: false,
                donor_motif: [b'N', b'N'],
                acceptor_motif: [b'N', b'N'],
                motif_class: SpliceMotifClass::Unknown,
            });
        }

        let d = [donor[0], donor[1]];
        let a = [acceptor[0], acceptor[1]];

        let (is_canonical, motif_class) = classify_splice_motif(d, a, strand);

        Ok(SpliceMotifResult {
            is_canonical,
            donor_motif: d,
            acceptor_motif: a,
            motif_class,
        })
    }

    /// Strand on which an intron's donor and acceptor dinucleotides form a canonical
    /// motif, when exactly one strand does.
    pub fn splice_motif_strand(
        &self,
        chrom: &str,
        intron_start: i64,
        intron_end: i64,
    ) -> IsoForgeResult<Option<Strand>> {
        let donor = self.fetch_seq(chrom, intron_start - 1, intron_start + 1)?;
        let acceptor = self.fetch_seq(chrom, intron_end - 2, intron_end)?;
        if donor.len() < 2 || acceptor.len() < 2 {
            return Ok(None);
        }
        let d = [donor[0], donor[1]];
        let a = [acceptor[0], acceptor[1]];
        let plus = classify_splice_motif(d, a, Strand::Plus).0;
        let minus = classify_splice_motif(d, a, Strand::Minus).0;
        Ok(match (plus, minus) {
            (true, false) => Some(Strand::Plus),
            (false, true) => Some(Strand::Minus),
            _ => None,
        })
    }

    /// Try to correct a non-canonical splice junction by shifting it within a window. Searches ±max_offset positions
    /// around the original junction for the nearest position that produces a canonical motif.
    /// Returns: Some((new_intron_start, new_intron_end)) if correctable, None otherwise.
    pub fn correct_splice_junction(
        &self,
        chrom: &str,
        intron_start: i64,
        intron_end: i64,
        strand: Strand,
        max_offset: i64,
    ) -> IsoForgeResult<Option<(i64, i64)>> {
        // First check if already canonical
        let current = self.check_splice_motif(chrom, intron_start, intron_end, strand)?;
        if current.is_canonical {
            return Ok(None); // No correction needed
        }

        // Try offsets, preferring smaller shifts
        let mut best: Option<(i64, i64, i64)> = None; // (new_start, new_end, total_shift)

        for donor_offset in -max_offset..=max_offset {
            for acceptor_offset in -max_offset..=max_offset {
                let new_start = intron_start + donor_offset;
                let new_end = intron_end + acceptor_offset;

                // Skip invalid coordinates
                if new_start <= 0 || new_end <= 0 || new_start >= new_end {
                    continue;
                }

                // Skip if this is the original position
                if donor_offset == 0 && acceptor_offset == 0 {
                    continue;
                }

                if let Ok(result) = self.check_splice_motif(chrom, new_start, new_end, strand) {
                    if result.is_canonical {
                        let total_shift = donor_offset.abs() + acceptor_offset.abs();
                        if best.is_none() || total_shift < best.unwrap().2 {
                            best = Some((new_start, new_end, total_shift));
                        }
                    }
                }
            }
        }

        Ok(best.map(|(s, e, _)| (s, e)))
    }

    /// Compute the MD tag from CIGAR, read sequence and reference, for alignments
    /// that lack one.
    pub fn compute_md_tag(
        &self,
        chrom: &str,
        ref_start: i64,       // 1-based start
        cigar_ops: &[(u32, u8)], // (length, op_code) pairs
        read_seq: &[u8],
    ) -> IsoForgeResult<String> {
        let mut md = String::new();
        let mut ref_pos = ref_start - 1; // 0-based
        let mut read_pos: usize = 0;
        let mut matches = 0u64;

        for &(len, op) in cigar_ops {
            let len = len as i64;
            match op {
                0 | 7 | 8 => {
                    // M, =, X — consumes both ref and read
                    let ref_seq = self.fetch_seq(chrom, ref_pos, ref_pos + len)?;
                    for i in 0..len as usize {
                        if read_pos + i < read_seq.len()
                            && i < ref_seq.len()
                            && read_seq[read_pos + i].to_ascii_uppercase()
                                == ref_seq[i].to_ascii_uppercase()
                        {
                            matches += 1;
                        } else {
                            // Mismatch
                            md.push_str(&matches.to_string());
                            matches = 0;
                            if i < ref_seq.len() {
                                md.push(ref_seq[i] as char);
                            } else {
                                md.push('N');
                            }
                        }
                    }
                    ref_pos += len;
                    read_pos += len as usize;
                }
                1 | 4 => {
                    // I, S — consumes read only
                    read_pos += len as usize;
                }
                2 => {
                    // D — consumes ref only (deletion from read)
                    md.push_str(&matches.to_string());
                    matches = 0;
                    md.push('^');
                    let del_seq = self.fetch_seq(chrom, ref_pos, ref_pos + len)?;
                    for &b in &del_seq {
                        md.push(b as char);
                    }
                    ref_pos += len;
                }
                3 => {
                    // N — intron skip, consumes ref only
                    ref_pos += len;
                }
                5 | 6 => {
                    // H, P — no consume
                }
                _ => {}
            }
        }

        // Final match count
        md.push_str(&matches.to_string());

        Ok(md)
    }
}

// ─── Splice motif types ──────────────────────────────────────────────────────

/// Result of checking a splice junction motif.
#[derive(Debug, Clone)]
pub struct SpliceMotifResult {
    pub is_canonical: bool,
    pub donor_motif: [u8; 2],
    pub acceptor_motif: [u8; 2],
    pub motif_class: SpliceMotifClass,
}

/// Classification of splice site motif.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpliceMotifClass {
    /// GT-AG (CT-AC on the minus strand)
    GtAg,
    /// GC-AG (CT-GC on the minus strand)
    GcAg,
    /// AT-AC (GT-AT on the minus strand)
    AtAc,
    /// Non-canonical
    NonCanonical,
    /// Could not determine (sequence fetch failed)
    Unknown,
}

/// `classify_splice_motif`'s class alone, for callers outside this module.
pub fn classify_splice_motif_class(donor: [u8; 2], acceptor: [u8; 2], strand: Strand) -> SpliceMotifClass {
    classify_splice_motif(donor, acceptor, strand).1
}

/// Classify a splice motif pair as canonical or non-canonical. On + strand, canonical motifs are GT-AG, GC-AG and
/// AT-AC; on - strand, the reverse complements CT-AC, CT-GC and GT-AT.
fn classify_splice_motif(donor: [u8; 2], acceptor: [u8; 2], strand: Strand) -> (bool, SpliceMotifClass) {
    match strand {
        Strand::Plus => {
            if donor == [b'G', b'T'] && acceptor == [b'A', b'G'] {
                (true, SpliceMotifClass::GtAg)
            } else if donor == [b'G', b'C'] && acceptor == [b'A', b'G'] {
                (true, SpliceMotifClass::GcAg)
            } else if donor == [b'A', b'T'] && acceptor == [b'A', b'C'] {
                (true, SpliceMotifClass::AtAc)
            } else {
                (false, SpliceMotifClass::NonCanonical)
            }
        }
        Strand::Minus => {
            // Reverse complement: read the motif from the other end
            if donor == [b'C', b'T'] && acceptor == [b'A', b'C'] {
                (true, SpliceMotifClass::GtAg)
            } else if donor == [b'C', b'T'] && acceptor == [b'G', b'C'] {
                (true, SpliceMotifClass::GcAg)
            } else if donor == [b'G', b'T'] && acceptor == [b'A', b'T'] {
                (true, SpliceMotifClass::AtAc)
            } else {
                (false, SpliceMotifClass::NonCanonical)
            }
        }
    }
}

/// Reverse complement a DNA sequence.
fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'a' => b't',
            b't' => b'a',
            b'c' => b'g',
            b'g' => b'c',
            _ => b'N',
        })
        .collect()
}

// ─── Splice junction database ────────────────────────────────────────────────

/// Known splice junctions from short-read data, used to validate and
/// prioritise splice junction corrections.
#[derive(Debug, Default)]
pub struct SpliceJunctionDb {
    /// Junction record by (chromosome, intron start, intron end).
    junctions: HashMap<(String, i64, i64), SjInfo>,
}

/// Information about a known splice junction from short-read alignment.
#[derive(Debug, Clone)]
pub struct SjInfo {
    pub strand: i8,        // 0=undefined, 1=+, 2=-
    pub motif: i8,         // 0=non-canonical, 1=GT/AG, 2=CT/AC, 3=GC/AG, 4=CT/GC, 5=AT/AC, 6=GT/AT
    pub annotated: bool,   // from GTF annotation
    pub unique_reads: i64, // unique mapping reads
    pub multi_reads: i64,  // multi-mapping reads
}

impl SpliceJunctionDb {
    /// Load from a STAR SJ.out.tab file (1-based start and end).
    ///
    /// Format: chrom start end strand motif annotated unique_reads multi_reads overhang
    pub fn from_star_sj_tab<P: AsRef<Path>>(path: P) -> IsoForgeResult<Self> {
        let mut db = Self::default();
        let content = std::fs::read_to_string(path.as_ref())
            .map_err(|e| IsoForgeError::IoError(e))?;

        for line in content.lines() {
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 9 {
                continue;
            }

            let chrom = fields[0].to_string();
            let start: i64 = fields[1].parse().unwrap_or(0);
            let end: i64 = fields[2].parse().unwrap_or(0);
            let strand: i8 = fields[3].parse().unwrap_or(0);
            let motif: i8 = fields[4].parse().unwrap_or(0);
            let annotated: bool = fields[5] == "1";
            let unique_reads: i64 = fields[6].parse().unwrap_or(0);
            let multi_reads: i64 = fields[7].parse().unwrap_or(0);

            db.junctions.insert(
                (chrom, start, end),
                SjInfo {
                    strand,
                    motif,
                    annotated,
                    unique_reads,
                    multi_reads,
                },
            );
        }

        Ok(db)
    }

    /// Load and merge multiple STAR SJ.out.tab files.
    pub fn from_star_sj_tabs(paths: &[PathBuf]) -> IsoForgeResult<Self> {
        let mut merged = Self::default();
        for path in paths {
            let db = Self::from_star_sj_tab(path)?;
            for (key, info) in db.junctions {
                let entry = merged.junctions.entry(key).or_insert(SjInfo {
                    strand: info.strand,
                    motif: info.motif,
                    annotated: info.annotated,
                    unique_reads: 0,
                    multi_reads: 0,
                });
                entry.unique_reads += info.unique_reads;
                entry.multi_reads += info.multi_reads;
                if info.annotated {
                    entry.annotated = true;
                }
            }
        }
        Ok(merged)
    }

    /// Check if a splice junction is in the database.
    pub fn contains(&self, chrom: &str, start: i64, end: i64) -> bool {
        self.junctions.contains_key(&(chrom.to_string(), start, end))
    }

    /// Look up a splice junction.
    pub fn get(&self, chrom: &str, start: i64, end: i64) -> Option<&SjInfo> {
        self.junctions.get(&(chrom.to_string(), start, end))
    }

    /// Number of junctions in the database.
    pub fn len(&self) -> usize {
        self.junctions.len()
    }

    pub fn is_empty(&self) -> bool {
        self.junctions.is_empty()
    }

    /// All splice positions (donor starts and acceptor ends), each widened by
    /// ±`fuzzy` bp, as (chromosome, position) pairs. Positions in this set have
    /// short-read support and are exempt from extended snapping.
    pub fn splice_positions(&self, fuzzy: i64) -> rustc_hash::FxHashSet<(String, i64)> {
        let mut positions = rustc_hash::FxHashSet::default();
        for (chrom, start, end) in self.junctions.keys() {
            for offset in -fuzzy..=fuzzy {
                positions.insert((chrom.clone(), start + offset));
                positions.insert((chrom.clone(), end + offset));
            }
        }
        positions
    }
}

// ─── Tests ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_classify_splice_motif_plus() {
        // GT-AG on + strand
        let (is_canon, class) = classify_splice_motif([b'G', b'T'], [b'A', b'G'], Strand::Plus);
        assert!(is_canon);
        assert_eq!(class, SpliceMotifClass::GtAg);

        // GC-AG on + strand
        let (is_canon, class) = classify_splice_motif([b'G', b'C'], [b'A', b'G'], Strand::Plus);
        assert!(is_canon);
        assert_eq!(class, SpliceMotifClass::GcAg);

        // Non-canonical
        let (is_canon, _class) = classify_splice_motif([b'A', b'A'], [b'T', b'T'], Strand::Plus);
        assert!(!is_canon);
    }

    /// A genome of C, 500 bases on chr1, with each sequence of `sets` written at its 0-based offset.
    fn genome_with(sets: &[(usize, &[u8])]) -> (tempfile::TempDir, GenomeRef) {
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("g.fa");
        let mut g = vec![b'C'; 500];
        for (at, seq) in sets {
            g[*at..*at + seq.len()].copy_from_slice(seq);
        }
        std::fs::write(&fasta, format!(">chr1\n{}\n", String::from_utf8(g).unwrap())).unwrap();
        std::fs::write(dir.path().join("g.fa.fai"), "chr1\t500\t6\t500\t501\n").unwrap();
        (dir, GenomeRef::open(fasta).unwrap())
    }

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement(b"AATAAA"), b"TTTATT");
        assert_eq!(reverse_complement(b"ATTAAA"), b"TTTAAT");
    }

    #[test]
    fn poly_a_signal_windows_are_the_same_on_both_strands() {
        // + strand: 3' end at 1-based 200 (0-based 199); AATAAA at 0-based 176..182, 18 bases
        // (182..=199) from its last base to the 3' end base inclusive.
        // - strand: 3' end at 1-based 301 (0-based 300); the same hexamer read on the reverse
        // strand, TTTATT at 0-based 318..324.
        let (_d, g) = genome_with(&[(176, b"AATAAA"), (318, b"TTTATT")]);
        assert_eq!(g.scan_pas_upstream("chr1", 200, Strand::Plus), (true, Some("AATAAA".into()), Some(18)));
        assert_eq!(g.scan_pas_upstream("chr1", 301, Strand::Minus), (true, Some("AATAAA".into()), Some(18)));
        // A downstream element of exactly 16 U in 40 bases: + strand after 1-based 200 is
        // 0-based 200..240; - strand after 1-based 301 is 0-based 260..300, read reversed.
        let (_d2, g2) = genome_with(&[(200, b"TTTTTTTTTTTTTTTT"), (260, b"AAAAAAAAAAAAAAAA")]);
        assert!(g2.downstream_element("chr1", 200, Strand::Plus));
        assert!(g2.downstream_element("chr1", 301, Strand::Minus));
        assert!(!g2.downstream_element("chr1", 302, Strand::Minus));
    }

    #[test]
    fn test_classify_splice_motif_minus() {
        // CT-AC on - strand (= GT-AG on + strand)
        let (is_canon, class) = classify_splice_motif([b'C', b'T'], [b'A', b'C'], Strand::Minus);
        assert!(is_canon);
        assert_eq!(class, SpliceMotifClass::GtAg);
    }
}

