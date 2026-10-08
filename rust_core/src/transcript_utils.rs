// Transcript utility functions.

use crate::types::Strand;

/// Signed distance between two positions in the direction of the strand.
pub fn compute_delta(orig_pos: i64, new_pos: i64, strand: Strand) -> i64 {
    let abs_dist = (orig_pos - new_pos).abs();
    match strand {
        Strand::Plus => {
            if new_pos < orig_pos {
                -abs_dist
            } else {
                abs_dist
            }
        }
        Strand::Minus => {
            if new_pos < orig_pos {
                abs_dist
            } else {
                -abs_dist
            }
        }
    }
}

/// Check whether all SJs (introns) in the novelty list are known.
/// Input novelty list should not include first or last exon.
/// Introns are at even indices (0, 2, 4, ...) in the splice-edge novelty list.
pub fn check_all_sjs_known(novelty: &[u8]) -> bool {
    if novelty.len() == 1 {
        return novelty[0] == 0;
    }
    // Introns are at even indices
    novelty.iter().step_by(2).all(|&n| n == 0)
}

/// Check whether all exons in the novelty list are known.
/// Exons are at odd indices (1, 3, 5, ...) in the splice-edge novelty list.
pub fn check_all_exons_known(novelty: &[u8]) -> bool {
    if novelty.len() == 1 {
        return novelty[0] == 0;
    }
    // Exons are at odd indices
    novelty.iter().skip(1).step_by(2).all(|&n| n == 0)
}

/// Check whether all exons in the novelty list are novel.
pub fn check_all_exons_novel(novelty: &[u8]) -> bool {
    if novelty.len() == 1 {
        return false; // No exons to analyze
    }
    // Exons are at odd indices
    novelty.iter().skip(1).step_by(2).all(|&n| n == 1)
}

/// Construct gene and transcript names from IDs
pub fn construct_names(
    gene_id: i64,
    transcript_id: i64,
    prefix: &str,
    n_places: usize,
) -> (String, String) {
    let gene_name = format!("{}G{:0>width$}", prefix, gene_id, width = n_places);
    let transcript_name = format!("{}T{:0>width$}", prefix, transcript_id, width = n_places);
    (gene_name, transcript_name)
}

/// Get overlap between two intervals (inclusive endpoints)
pub fn get_overlap(a: (i64, i64), b: (i64, i64)) -> (i64, Option<f64>) {
    let overlap = 0i64.max(a.1.min(b.1) - a.0.max(b.0) + 1);
    let ref_len = (b.1 - b.0).abs();
    let perc_overlap = if ref_len != 0 {
        Some((overlap as f64 / ref_len as f64) * 100.0)
    } else {
        None
    };
    (overlap, perc_overlap)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compute_delta_plus_strand() {
        assert_eq!(compute_delta(100, 90, Strand::Plus), -10);
        assert_eq!(compute_delta(100, 110, Strand::Plus), 10);
        assert_eq!(compute_delta(100, 100, Strand::Plus), 0);
    }

    #[test]
    fn test_compute_delta_minus_strand() {
        assert_eq!(compute_delta(100, 90, Strand::Minus), 10);
        assert_eq!(compute_delta(100, 110, Strand::Minus), -10);
        assert_eq!(compute_delta(100, 100, Strand::Minus), 0);
    }

    #[test]
    fn test_all_sjs_known() {
        assert!(check_all_sjs_known(&[0]));
        assert!(!check_all_sjs_known(&[1]));
        assert!(check_all_sjs_known(&[0, 0, 0])); // intron, exon, intron
        assert!(!check_all_sjs_known(&[1, 0, 0]));
    }

    #[test]
    fn test_all_exons_known() {
        assert!(check_all_exons_known(&[0]));
        assert!(check_all_exons_known(&[0, 0, 0])); // intron, exon, intron
        assert!(!check_all_exons_known(&[0, 1, 0]));
    }

    #[test]
    fn test_construct_names() {
        let (gene, transcript) = construct_names(42, 100, "ISOFOR", 9);
        assert_eq!(gene, "ISOFORG000000042");
        assert_eq!(transcript, "ISOFORT000000100");
    }

    #[test]
    fn test_overlap() {
        let (overlap, _) = get_overlap((10, 20), (15, 25));
        assert_eq!(overlap, 6);

        let (overlap, _) = get_overlap((10, 20), (30, 40));
        assert_eq!(overlap, 0);
    }
}
