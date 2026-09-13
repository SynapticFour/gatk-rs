//! GATK `ReadPosRankSumTest` via `MannWhitneyU` (list-based parity gate).
//!
//! Java `RankSumTest.annotate` (GATK 4.4.0.0): finite z (including 0.0) →
//! `String.format("%.3f", zScore)` map entry; `Double.isNaN(z)` or empty
//! lists → `Collections.emptyMap()` (no INFO key). 6R.172: `Option<f64>`
//! restores that undefined-vs-zero distinction. Do not treat `0.0` as omit.

use crate::mann_whitney_u::{MannWhitneyU, TestType};

/// Mann-Whitney z-score for REF vs ALT read positions within reads.
///
/// `None` means the statistic is undefined (empty REF, empty ALT, or NaN).
/// `Some(0.0)` is a legitimate finite result and must remain emittable.
pub fn read_pos_rank_sum(ref_positions: &[f64], alt_positions: &[f64]) -> Option<f64> {
    rank_sum_z(alt_positions, ref_positions)
}

fn rank_sum_z(alt: &[f64], reference: &[f64]) -> Option<f64> {
    if alt.is_empty() || reference.is_empty() {
        return None;
    }
    let result = MannWhitneyU::default().test(alt, reference, TestType::FirstDominates);
    if result.z.is_nan() {
        None
    } else {
        Some(result.z)
    }
}

#[cfg(test)]
mod forensic_6r172 {
    use super::read_pos_rank_sum;

    /// 6R.172 matrix: undefined is `None`, not numeric zero.
    #[test]
    fn forensic_6r172_semantic_matrix_production_fn() {
        let a = read_pos_rank_sum(&[10.0, 20.0], &[30.0, 40.0]);
        let b = read_pos_rank_sum(&[10.0, 20.0], &[10.0, 20.0]);
        let c = read_pos_rank_sum(&[], &[51.0, 61.0, 70.0]);
        let d = read_pos_rank_sum(&[10.0, 20.0], &[]);
        let e = read_pos_rank_sum(&[], &[]);
        eprintln!("6R172\tA_nonzero\t{a:?}");
        eprintln!("6R172\tB_valid_zero\t{b:?}");
        eprintln!("6R172\tC_no_ref\t{c:?}");
        eprintln!("6R172\tD_no_alt\t{d:?}");
        eprintln!("6R172\tE_both_empty\t{e:?}");
        assert!(
            a.is_some_and(|z| z.is_finite() && z != 0.0),
            "case A finite non-zero, got {a:?}"
        );
        assert_eq!(b, Some(0.0), "case B legitimate zero must remain Some(0.0)");
        assert_eq!(
            c, None,
            "case C target-equivalent must be None, not Some(0.0)"
        );
        assert_eq!(d, None, "case D must be None");
        assert_eq!(e, None, "case E must be None");
        assert_ne!(
            c, b,
            "undefined must not collapse onto the legitimate-zero representation"
        );
    }
}
