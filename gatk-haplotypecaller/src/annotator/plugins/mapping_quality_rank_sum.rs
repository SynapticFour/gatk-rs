//! GATK `MappingQualityRankSumTest` via `MannWhitneyU` (list-based parity gate).

use crate::mann_whitney_u::{MannWhitneyU, TestType};

/// Mann-Whitney z-score for REF vs ALT mapping qualities.
pub fn mapping_quality_rank_sum(ref_mqs: &[u8], alt_mqs: &[u8]) -> f64 {
    if ref_mqs.is_empty() || alt_mqs.is_empty() {
        return 0.0;
    }
    let ref_f: Vec<f64> = ref_mqs.iter().map(|&q| f64::from(q)).collect();
    let alt_f: Vec<f64> = alt_mqs.iter().map(|&q| f64::from(q)).collect();
    let result = MannWhitneyU::default().test(&alt_f, &ref_f, TestType::FirstDominates);
    if result.z.is_nan() {
        0.0
    } else {
        result.z
    }
}
