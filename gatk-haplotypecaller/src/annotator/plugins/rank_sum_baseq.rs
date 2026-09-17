//! GATK `BaseQualityRankSumTest` via `MannWhitneyU` (production).
//!
//! Java `RankSumTest.annotate` (GATK 4.4.0.0): finite z (including 0.0) →
//! `String.format("%.3f", zScore)` map entry; `Double.isNaN(z)` or empty
//! lists → `Collections.emptyMap()` (no INFO key). Reuse 6R.172 `Option<f64>`:
//! do not treat `0.0` as omit.

use crate::mann_whitney_u::{MannWhitneyU, TestType};

/// Mann-Whitney z-score for REF vs ALT base qualities.
///
/// `None` means the statistic is undefined (empty REF, empty ALT, or NaN).
/// `Some(0.0)` is a legitimate finite result and must remain emittable.
pub fn base_quality_rank_sum(ref_quals: &[f64], alt_quals: &[f64]) -> Option<f64> {
    rank_sum_z(alt_quals, ref_quals)
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
mod forensic_6r193_zero {
    use super::base_quality_rank_sum;

    #[test]
    fn forensic_6r193_finite_zero_is_some() {
        assert_eq!(
            base_quality_rank_sum(&[30.0], &[30.0]),
            Some(0.0),
            "equal singleton lists are a legitimate finite zero"
        );
        assert_eq!(base_quality_rank_sum(&[], &[30.0]), None);
        assert_eq!(base_quality_rank_sum(&[30.0], &[]), None);
        assert_eq!(base_quality_rank_sum(&[], &[]), None);
    }
}
