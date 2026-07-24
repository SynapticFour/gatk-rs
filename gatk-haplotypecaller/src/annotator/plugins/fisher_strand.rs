//! GATK `FisherStrand` — Fisher's exact test on ref/alt × forward/reverse counts (production).

use crate::activity_scoring::log10_sum_log10;
use statrs::distribution::{Discrete, Hypergeometric};

/// GATK `FisherStrand.MIN_PVALUE`.
pub const MIN_PVALUE: f64 = 1e-320;
/// GATK `FisherStrand.TARGET_TABLE_SIZE`.
const TARGET_TABLE_SIZE: f64 = 200.0;
/// GATK `FisherExactTest.REL_ERR`.
const REL_ERR: f64 = 1.0 - 1e-7;
const LOG10_E: f64 = std::f64::consts::LOG10_E;

/// Phred-scaled FS annotation value (GATK `makeValueObjectForAnnotation` without the `%.3f` round-trip).
pub fn fisher_strand_statistic(ref_fw: u32, ref_rv: u32, alt_fw: u32, alt_rv: u32) -> f64 {
    let table = [
        [ref_fw as i32, ref_rv as i32],
        [alt_fw as i32, alt_rv as i32],
    ];
    let p = p_value_for_contingency_table(&table);
    phred_scale_error_rate(p.max(MIN_PVALUE))
}

/// GATK `FisherStrand.pValueForContingencyTable`.
pub fn p_value_for_contingency_table(original_table: &[[i32; 2]; 2]) -> f64 {
    let normalized = normalize_contingency_table(original_table);
    fisher_exact_two_sided_p_value(&normalized)
}

fn normalize_contingency_table(table: &[[i32; 2]; 2]) -> [[i32; 2]; 2] {
    let sum = table[0][0]
        .checked_add(table[0][1])
        .and_then(|s| s.checked_add(table[1][0]))
        .and_then(|s| s.checked_add(table[1][1]))
        .unwrap_or(i32::MAX);
    if (sum as f64) <= TARGET_TABLE_SIZE * 2.0 {
        return *table;
    }
    let norm_factor = (sum as f64) / TARGET_TABLE_SIZE;
    [
        [
            (table[0][0] as f64 / norm_factor) as i32,
            (table[0][1] as f64 / norm_factor) as i32,
        ],
        [
            (table[1][0] as f64 / norm_factor) as i32,
            (table[1][1] as f64 / norm_factor) as i32,
        ],
    ]
}

/// GATK `FisherExactTest.twoSidedPValue` (R `fisher.test` two-sided).
fn fisher_exact_two_sided_p_value(table: &[[i32; 2]; 2]) -> f64 {
    let m = table[0][0] + table[0][1];
    let n = table[1][0] + table[1][1];
    let k = table[0][0] + table[1][0];
    let lo = 0.max(k - n);
    let hi = k.min(m);
    if lo >= hi {
        return 1.0;
    }
    let population = (m + n) as u64;
    let successes = m as u64;
    let sample = k as u64;
    let dist = match Hypergeometric::new(population, successes, sample) {
        Ok(d) => d,
        Err(_) => return 1.0,
    };
    let observed = table[0][0];
    let mut log_ds = Vec::with_capacity((hi - lo + 1) as usize);
    for x in lo..=hi {
        let lp = dist.ln_pmf(x as u64);
        log_ds.push(lp);
    }
    let obs_idx = (observed - lo) as usize;
    let threshold = log_ds[obs_idx] * REL_ERR;
    let log10_ds: Vec<f64> = log_ds
        .iter()
        .copied()
        .filter(|ln| *ln <= threshold)
        .map(|ln| ln * LOG10_E)
        .collect();
    if log10_ds.is_empty() {
        return 1.0;
    }
    let log10_sum = log10_sum_log10(&log10_ds);
    (10_f64.powf(log10_sum)).min(1.0)
}

/// GATK `QualityUtils.phredScaleErrorRate`.
fn phred_scale_error_rate(error_rate: f64) -> f64 {
    let log10_err = error_rate.log10();
    phred_scale_log10_error_rate(log10_err)
}

/// GATK `QualityUtils.phredScaleLog10ErrorRate` (`MIN_LOG10_SCALED_QUAL` = log10 of Java `Double.MIN_VALUE`).
fn phred_scale_log10_error_rate(error_rate_log10: f64) -> f64 {
    const MIN_LOG10_SCALED_QUAL: f64 = -323.306_103_180_901_6;
    (-10.0 * error_rate_log10.max(MIN_LOG10_SCALED_QUAL)).abs()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn table(ref_fw: i32, ref_rv: i32, alt_fw: i32, alt_rv: i32) -> [[i32; 2]; 2] {
        [[ref_fw, ref_rv], [alt_fw, alt_rv]]
    }

    #[test]
    fn gatk_unit_test_vectors() {
        let cases: &[(i32, i32, i32, i32, f64)] = &[
            (9, 11, 12, 10, 0.7578618),
            (12, 10, 9, 11, 0.7578618),
            (9, 10, 12, 10, 0.7578618),
            (9, 9, 12, 10, 1.0),
            (9, 13, 12, 10, 0.5466948),
            (12, 10, 9, 13, 0.5466948),
            (9, 12, 11, 9, 0.5377362),
            (0, 0, 0, 3, 1.0),
            (9, 0, 0, 0, 1.0),
            (0, 0, 0, 0, 1.0),
            (100000, 100000, 100000, 100000, 1.0),
            (0, 0, 100000, 100000, 1.0),
            (66, 14, 64, 4, 0.04243330),
            (137, 159, 9, 23, 0.06088506),
        ];
        for &(rf, rr, af, ar, expected) in cases {
            let p = p_value_for_contingency_table(&table(rf, rr, af, ar));
            assert!(
                (p - expected).abs() < 1e-6,
                "p-value mismatch for [{rf},{rr};{af},{ar}]: got {p}, expected {expected}"
            );
        }
    }

    #[test]
    fn parity_fixture_strand_bias() {
        let fs = fisher_strand_statistic(10, 2, 2, 10);
        assert!(
            fs > 10.0,
            "strong strand bias should yield high FS, got {fs}"
        );
    }
}
