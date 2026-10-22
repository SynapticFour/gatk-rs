//! GATK `FisherStrand` — Fisher's exact test on ref/alt × forward/reverse counts (production).

use crate::activity_scoring::log10_sum_log10;

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
    // Former `Hypergeometric::new` rejection: successes or draws above the population.
    if successes > population || sample > population {
        return 1.0;
    }
    let observed = table[0][0];
    let mut log_ds = Vec::with_capacity((hi - lo + 1) as usize);
    for x in lo..=hi {
        let lp = hypergeometric_log_probability(m + n, m, k, x);
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

/// Apache Commons Math 3.5 `HypergeometricDistribution.logProbability`.
///
/// `N` is the population, `K` the number of successes, `n` the sample size.
fn hypergeometric_log_probability(population: i32, successes: i32, sample: i32, x: i32) -> f64 {
    let p = f64::from(sample) / f64::from(population);
    let q = f64::from(population - sample) / f64::from(population);
    let p1 = log_binomial_probability(x, successes, p, q);
    let p2 = log_binomial_probability(sample - x, population - successes, p, q);
    let p3 = log_binomial_probability(sample, population, p, q);
    p1 + p2 - p3
}

/// Apache Commons Math 3.5 `SaddlePointExpansion.logBinomialProbability`.
fn log_binomial_probability(x: i32, n: i32, p: f64, q: f64) -> f64 {
    if x == 0 {
        if p < 0.1 {
            -get_deviance_part(f64::from(n), f64::from(n) * q) - f64::from(n) * p
        } else {
            f64::from(n) * q.ln()
        }
    } else if x == n {
        if q < 0.1 {
            -get_deviance_part(f64::from(n), f64::from(n) * p) - f64::from(n) * q
        } else {
            f64::from(n) * p.ln()
        }
    } else {
        let ret = get_stirling_error(f64::from(n))
            - get_stirling_error(f64::from(x))
            - get_stirling_error(f64::from(n - x))
            - get_deviance_part(f64::from(x), f64::from(n) * p)
            - get_deviance_part(f64::from(n - x), f64::from(n) * q);
        let f = (std::f64::consts::PI * 2.0 * f64::from(x) * f64::from(n - x)) / f64::from(n);
        -0.5 * f.ln() + ret
    }
}

/// Apache Commons Math 3.5 `SaddlePointExpansion.getStirlingError`.
fn get_stirling_error(z: f64) -> f64 {
    /// Exact Stirling expansion error for `z = 0, 0.5, ..., 15`.
    const EXACT_STIRLING_ERRORS: [f64; 31] = [
        0.0,
        0.1534264097200273452913848,
        0.0810614667953272582196702,
        0.0548141210519176538961390,
        0.0413406959554092940938221,
        0.03316287351993628748511048,
        0.02767792568499833914878929,
        0.02374616365629749597132920,
        0.02079067210376509311152277,
        0.01848845053267318523077934,
        0.01664469118982119216319487,
        0.01513497322191737887351255,
        0.01387612882307074799874573,
        0.01281046524292022692424986,
        0.01189670994589177009505572,
        0.01110455975820691732662991,
        0.010411265261972096497478567,
        0.009799416126158803298389475,
        0.009255462182712732917728637,
        0.008768700134139385462952823,
        0.008330563433362871256469318,
        0.007934114564314020547248100,
        0.007573675487951840794972024,
        0.007244554301320383179543912,
        0.006942840107209529865664152,
        0.006665247032707682442354394,
        0.006408994188004207068439631,
        0.006171712263039457647532867,
        0.005951370112758847735624416,
        0.005746216513010115682023589,
        0.005554733551962801371038690,
    ];
    if z < 15.0 {
        let z2 = 2.0 * z;
        if z2.floor() == z2 {
            EXACT_STIRLING_ERRORS[z2 as usize]
        } else {
            // `Gamma.logGamma(z + 1) - (z + 0.5) * ln(z) + z - 0.5 * ln(2π)`.
            // Integer arguments from `log_binomial_probability` take the table arm.
            f64::NAN
        }
    } else {
        let z2 = z * z;
        (0.083333333333333333333
            - (0.00277777777777777777778
                - (0.00079365079365079365079365
                    - (0.000595238095238095238095238 - 0.0008417508417508417508417508 / z2) / z2)
                    / z2)
                / z2)
            / z
    }
}

/// Apache Commons Math 3.5 `SaddlePointExpansion.getDeviancePart`.
fn get_deviance_part(x: f64, mu: f64) -> f64 {
    if (x - mu).abs() < 0.1 * (x + mu) {
        let d = x - mu;
        let mut v = d / (x + mu);
        let mut s1 = v * d;
        let mut s = f64::NAN;
        let mut ej = 2.0 * x * v;
        v *= v;
        let mut j = 1i32;
        while s1 != s {
            s = s1;
            ej *= v;
            s1 = s + ej / f64::from((j * 2) + 1);
            j += 1;
        }
        s1
    } else {
        x * (x / mu).ln() + mu - x
    }
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

    /// Java `HypergeometricDistribution.logProbability` bits for `N=102, K=83, n=70`, `x=51..=70`.
    /// Recorded from the GATK 4.4.0.0 Commons Math 3.5 jar in 7R.68/7R.69 (`doubleToRawLongBits`).
    #[test]
    fn java_hypergeometric_ln_bits_for_support_51_to_70() {
        const JAVA_LN_BITS: [i64; 20] = [
            -4602646795984669864,
            -4605383283164480118,
            -4607324206578857635,
            -4609906603009685541,
            -4611587024396473687,
            -4613314082875347134,
            -4613764818702010310,
            -4612939276175630970,
            -4611292139482635911,
            -4609680706417973933,
            -4607484346707785347,
            -4605940911827207960,
            -4604247471226293456,
            -4602459861411729866,
            -4601290475787714290,
            -4599940580755583262,
            -4598388092215394602,
            -4597386137601643653,
            -4596338127202424441,
            -4595056324810417985,
        ];
        for (i, x) in (51i32..=70).enumerate() {
            let lp = hypergeometric_log_probability(102, 83, 70, x);
            let java_bits = JAVA_LN_BITS[i] as u64;
            assert_eq!(
                lp.to_bits(),
                java_bits,
                "x={x} rust_bits={} java_bits={java_bits} rust={lp:?} java={:?}",
                lp.to_bits(),
                f64::from_bits(java_bits)
            );
            let log10 = lp * LOG10_E;
            let java_log10 = f64::from_bits(java_bits) * LOG10_E;
            assert_eq!(
                log10.to_bits(),
                java_log10.to_bits(),
                "log10 x={x} rust_bits={} java_bits={}",
                log10.to_bits(),
                java_log10.to_bits()
            );
        }
    }

    /// Java two-sided p and raw FS for contingency table `60 23 / 10 9`.
    #[test]
    fn java_fisher_table_60_23_10_9_bits() {
        const JAVA_P_BITS: u64 = 4592425646947568419;
        const JAVA_FS_BITS: u64 = 4621637562915352156;
        let p = p_value_for_contingency_table(&table(60, 23, 10, 9));
        let fs = fisher_strand_statistic(60, 23, 10, 9);
        assert_eq!(
            p.to_bits(),
            JAVA_P_BITS,
            "p rust_bits={} java_bits={JAVA_P_BITS} rust={p:?} java={:?}",
            p.to_bits(),
            f64::from_bits(JAVA_P_BITS)
        );
        assert_eq!(
            fs.to_bits(),
            JAVA_FS_BITS,
            "FS rust_bits={} java_bits={JAVA_FS_BITS} rust={fs:?} java={:?}",
            fs.to_bits(),
            f64::from_bits(JAVA_FS_BITS)
        );
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
