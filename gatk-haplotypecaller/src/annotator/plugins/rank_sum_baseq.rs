//! GATK `BaseQualityRankSumTest` via `MannWhitneyU` (production).

use crate::mann_whitney_u::{MannWhitneyU, TestType};

/// Mann-Whitney z-score for REF vs ALT base qualities (GATK `RankSumTest` / `BaseQualityRankSumTest`).
pub fn base_quality_rank_sum(ref_quals: &[u8], alt_quals: &[u8]) -> f64 {
    if ref_quals.is_empty() || alt_quals.is_empty() {
        return 0.0;
    }
    let ref_f: Vec<f64> = ref_quals.iter().map(|&q| f64::from(q)).collect();
    let alt_f: Vec<f64> = alt_quals.iter().map(|&q| f64::from(q)).collect();
    let result = MannWhitneyU::default().test(&alt_f, &ref_f, TestType::FirstDominates);
    if result.z.is_nan() {
        0.0
    } else {
        result.z
    }
}
