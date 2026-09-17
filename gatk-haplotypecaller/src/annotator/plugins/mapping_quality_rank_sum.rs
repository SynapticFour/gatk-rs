//! GATK `MappingQualityRankSumTest` via `MannWhitneyU` (production).
//!
//! Java `RankSumTest.annotate` (GATK 4.4.0.0): finite z (including 0.0) →
//! `String.format("%.3f", zScore)` map entry; `Double.isNaN(z)` or empty
//! lists → `Collections.emptyMap()` (no INFO key). Reuse 6R.172 `Option<f64>`:
//! do not treat `0.0` as omit.

use crate::mann_whitney_u::{MannWhitneyU, TestType};

/// Mann-Whitney z-score for REF vs ALT mapping qualities.
///
/// `None` means the statistic is undefined (empty REF, empty ALT, or NaN).
/// `Some(0.0)` is a legitimate finite result and must remain emittable.
pub fn mapping_quality_rank_sum(ref_mqs: &[f64], alt_mqs: &[f64]) -> Option<f64> {
    rank_sum_z(alt_mqs, ref_mqs)
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
mod forensic_6r194_zero {
    use super::mapping_quality_rank_sum;

    #[test]
    fn forensic_6r194_finite_zero_is_some() {
        assert_eq!(
            mapping_quality_rank_sum(&[60.0], &[60.0]),
            Some(0.0),
            "equal singleton lists are a legitimate finite zero"
        );
        assert_eq!(mapping_quality_rank_sum(&[], &[60.0]), None);
        assert_eq!(mapping_quality_rank_sum(&[60.0], &[]), None);
        assert_eq!(mapping_quality_rank_sum(&[], &[]), None);
    }

    #[test]
    fn forensic_6r194_unequal_singletons_match_java_millirounding() {
        let z = mapping_quality_rank_sum(&[60.0], &[40.0]).expect("finite");
        assert!(
            (z * 1000.0).round() == -674.0,
            "ALT-lower singleton pair is Java MQRankSum=-0.674, got {z}"
        );
        assert!(z.is_finite() && z != 0.0);
    }
}
