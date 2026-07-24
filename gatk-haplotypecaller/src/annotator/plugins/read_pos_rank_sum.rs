//! GATK `ReadPosRankSumTest` via `MannWhitneyU` (list-based parity gate).

use crate::mann_whitney_u::{MannWhitneyU, TestType};

/// Mann-Whitney z-score for REF vs ALT read positions within reads.
pub fn read_pos_rank_sum(ref_positions: &[f64], alt_positions: &[f64]) -> f64 {
    rank_sum_z(alt_positions, ref_positions)
}

fn rank_sum_z(alt: &[f64], reference: &[f64]) -> f64 {
    if alt.is_empty() || reference.is_empty() {
        return 0.0;
    }
    let result = MannWhitneyU::default().test(alt, reference, TestType::FirstDominates);
    if result.z.is_nan() {
        0.0
    } else {
        result.z
    }
}
