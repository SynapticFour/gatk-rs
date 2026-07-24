//! GATK `PairHMMLikelihoodCalculationEngine` base-quality capping.

use gatk_core::BaseQuality;

/// GATK `QualityUtils.MIN_USABLE_Q_SCORE`.
pub const MIN_USABLE_Q_SCORE: u8 = 6;

/// GATK `setToFixedValueIfTooLow`: values below `min_qual` become `fixed_qual`.
#[inline]
pub fn set_to_fixed_value_if_too_low(current: u8, min_qual: u8, fixed_qual: u8) -> u8 {
    if current < min_qual {
        fixed_qual
    } else {
        current
    }
}

/// Cap base qualities by mapping quality, then enforce `base_quality_score_threshold`.
/// Each byte is first normalized through [`BaseQuality::new`] (SAM max 93) so PairHMM
/// never sees uncapped Phred values at this boundary.
pub fn cap_read_base_qualities(
    quals: &mut [u8],
    mapq: u8,
    base_quality_score_threshold: u8,
    disable_cap_read_qualities_to_mapq: bool,
) {
    for q in quals.iter_mut() {
        let mut v = BaseQuality::new(*q).value();
        if !disable_cap_read_qualities_to_mapq {
            v = v.min(mapq);
        }
        *q = set_to_fixed_value_if_too_low(v, base_quality_score_threshold, MIN_USABLE_Q_SCORE);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn caps_low_bases_to_min_usable() {
        let mut q = vec![5u8, 18, 30];
        cap_read_base_qualities(&mut q, 60, 18, true);
        assert_eq!(q, vec![6, 18, 30]);
    }

    #[test]
    fn caps_by_mapq_when_enabled() {
        let mut q = vec![40u8, 40];
        cap_read_base_qualities(&mut q, 25, 6, false);
        assert_eq!(q, vec![25, 25]);
    }
}
