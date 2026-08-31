//! GATK `AlleleFrequencyCalculator` / `AFCalculationResult` slice.

use crate::activity_scoring::{
    log10_sum_log10, HaplotypeCallerActivityScoringParams,
    DEFAULT_HETEROZYGOSITY_STANDARD_DEVIATION, DEFAULT_SNP_HETEROZYGOSITY,
};
use gatk_common::GatkResult;

/// MLE allele frequency result (biallelic diploid parity v1).
/// # Invariants
/// `af` is the estimated alt allele frequency; log10 posteriors are in log10 probability space.
/// `em_iterations` counts EM steps used to converge the MLE.
/// # Ownership
/// Owned scalar bundle returned from AF calculators.
/// # Mutation
/// Immutable result snapshot.
/// # Biological assumptions
/// Biallelic diploid site with SNP-like Dirichlet prior (parity v1 slice).
/// # Java equivalence
/// GATK `AFCalculationResult` / `AlleleFrequencyCalculator` output.
#[derive(Debug, Clone, PartialEq)]
pub struct AfCalculationResult {
    pub alt_allele_count: i32,
    pub af: f64,
    pub log10_posterior_no_variant: f64,
    pub log10_p_no_alt: f64,
    pub em_iterations: usize,
}

/// Dirichlet / heterozygosity knobs for biallelic MLE allele-frequency calculation.
/// # Invariants
/// `snp_heterozygosity` and `heterozygosity_standard_deviation` feed the same prior family as
/// HC activity scoring defaults when constructed via [`Default`].
/// Pseudocounts are non-negative Dirichlet weights for ref vs SNP alleles.
/// # Ownership
/// [`Copy`] value type; safe to pass by value into calculators.
/// # Mutation
/// Field updates are caller-owned; calculators read a snapshot per call.
/// # Biological assumptions
/// Biallelic diploid site with SNP-like heterozygosity prior (not indel-specific AF model).
/// # Java equivalence
/// GATK 4.4 `AlleleFrequencyCalculator` / `AFCalculationResult` configuration slice.
#[derive(Debug, Clone, Copy)]
pub struct AfCalculatorConfig {
    /// Java config mirror; EM loop currently consumes derived pseudocounts.
    #[allow(dead_code)]
    pub snp_heterozygosity: f64,
    /// Java config mirror; EM loop currently consumes derived pseudocounts.
    #[allow(dead_code)]
    pub heterozygosity_standard_deviation: f64,
    pub ref_pseudocount: f64,
    pub snp_pseudocount: f64,
}

impl Default for AfCalculatorConfig {
    fn default() -> Self {
        let params = HaplotypeCallerActivityScoringParams::default();
        Self {
            snp_heterozygosity: DEFAULT_SNP_HETEROZYGOSITY,
            heterozygosity_standard_deviation: DEFAULT_HETEROZYGOSITY_STANDARD_DEVIATION,
            ref_pseudocount: params.ref_pseudocount(),
            snp_pseudocount: params.snp_pseudocount(),
        }
    }
}

fn log10_dirichlet_mean_weights(pseudocounts: &[f64]) -> Vec<f64> {
    let sum: f64 = pseudocounts.iter().sum();
    pseudocounts
        .iter()
        .map(|c| (c / sum).max(1e-300).log10())
        .collect()
}

fn diploid_genotype_allele_counts(genotype_index: usize) -> [usize; 2] {
    match genotype_index {
        0 => [2, 0],
        1 => [1, 1],
        _ => [0, 2],
    }
}

fn log10_normalized_genotype_posteriors_biallelic(
    log10_likelihoods: &[f64],
    log10_allele_frequencies: &[f64; 2],
) -> [f64; 3] {
    let mut log10_post = [f64::NEG_INFINITY; 3];
    for (gi, ll) in log10_likelihoods.iter().enumerate().take(3) {
        let [ref_c, alt_c] = diploid_genotype_allele_counts(gi);
        let log10_combo =
            crate::activity_scoring::log_binomial_coefficient_natural(2, alt_c as u32)
                * std::f64::consts::LOG10_E;
        log10_post[gi] = *ll
            + log10_combo
            + (ref_c as f64) * log10_allele_frequencies[0]
            + (alt_c as f64) * log10_allele_frequencies[1];
    }
    let log10_sum = log10_sum_log10(&log10_post);
    let mut out = log10_post;
    for x in &mut out {
        *x -= log10_sum;
    }
    out
}

fn effective_allele_counts_biallelic(
    samples_log10_likelihoods: &[&[f64]],
    log10_allele_frequencies: &[f64; 2],
) -> [f64; 2] {
    let mut log10_counts = [f64::NEG_INFINITY; 2];
    for gl in samples_log10_likelihoods {
        if gl.len() < 3 {
            continue;
        }
        let post = log10_normalized_genotype_posteriors_biallelic(gl, log10_allele_frequencies);
        for (gi, &log_p) in post.iter().enumerate() {
            let [ref_c, alt_c] = diploid_genotype_allele_counts(gi);
            if ref_c > 0 {
                log10_counts[0] =
                    log10_sum_log10(&[log10_counts[0], log_p + (ref_c as f64).log10()]);
            }
            if alt_c > 0 {
                log10_counts[1] =
                    log10_sum_log10(&[log10_counts[1], log_p + (alt_c as f64).log10()]);
            }
        }
    }
    [10_f64.powf(log10_counts[0]), 10_f64.powf(log10_counts[1])]
}

/// Biallelic diploid EM AF (GATK `AlleleFrequencyCalculator#calculate` core loop).
pub fn calculate_biallelic_af_em(
    samples_log10_likelihoods: &[&[f64]],
    config: &AfCalculatorConfig,
) -> GatkResult<AfCalculationResult> {
    const THRESHOLD: f64 = 0.1;
    let flat = -(2.0_f64).log10();
    let mut log10_af = [flat, flat];
    let mut allele_counts = [0.0_f64, 0.0];
    let mut iterations = 0usize;
    // GATK 4.4 `AlleleFrequencyCalculator.calculate` (SHA 2dbc0258): after each
    // `effectiveAlleleCounts`, always `log10AlleleFrequencies = Dirichlet(prior+counts).log10MeanWeights()`
    // then test the count delta. Breaking *before* that update used the previous AF for
    // `log10PNoVariant` and produced QUAL 78.583 vs Java 78.32 on PL 90,6,0.
    let mut allele_counts_maximum_difference = f64::INFINITY;
    while allele_counts_maximum_difference > THRESHOLD {
        iterations += 1;
        let new_counts = effective_allele_counts_biallelic(samples_log10_likelihoods, &log10_af);
        allele_counts_maximum_difference = new_counts
            .iter()
            .zip(allele_counts.iter())
            .map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max);
        allele_counts = new_counts;
        let posterior_pseudo = [
            config.ref_pseudocount + allele_counts[0],
            config.snp_pseudocount + allele_counts[1],
        ];
        let means = log10_dirichlet_mean_weights(&posterior_pseudo);
        log10_af = [means[0], means[1]];
        if iterations > 100 {
            break;
        }
    }

    let mut log10_p_no_variant = 0.0_f64;
    for gl in samples_log10_likelihoods {
        if gl.len() >= 3 {
            let post = log10_normalized_genotype_posteriors_biallelic(gl, &log10_af);
            log10_p_no_variant += post[0];
        }
    }
    log10_p_no_variant = log10_p_no_variant.min(0.0);

    let total_alleles = allele_counts[0] + allele_counts[1];
    let af = if total_alleles > 0.0 {
        allele_counts[1] / total_alleles
    } else {
        0.0
    };
    Ok(AfCalculationResult {
        alt_allele_count: allele_counts[1].round() as i32,
        af,
        log10_posterior_no_variant: log10_p_no_variant,
        log10_p_no_alt: log10_p_no_variant,
        em_iterations: iterations,
    })
}

/// Multi-allelic diploid AF EM result (GenotypeGVCFs joint genotyping).
/// # Invariants
/// `allele_frequencies.len == n_alleles` (REF + ALTs); sums to ~1.
/// `log10_posterior_no_variant` is Σ_s log10 P(hom-ref | data_s, π̂).
/// # Java equivalence
/// GATK `AlleleFrequencyCalculator` multi-allelic path used by `GenotypeGVCFsEngine`.
#[derive(Debug, Clone, PartialEq)]
pub struct MultiallelicAfResult {
    pub allele_frequencies: Vec<f64>,
    pub mle_allele_counts: Vec<f64>,
    pub log10_posterior_no_variant: f64,
    pub em_iterations: usize,
}

fn diploid_genotype_pairs(n_alleles: usize) -> Vec<(usize, usize)> {
    let mut pairs = Vec::new();
    for j in 0..n_alleles {
        for i in 0..=j {
            pairs.push((i, j));
        }
    }
    pairs
}

fn log10_normalized_genotype_posteriors_multi(
    log10_likelihoods: &[f64],
    log10_allele_frequencies: &[f64],
    pairs: &[(usize, usize)],
) -> Vec<f64> {
    let n = pairs.len().min(log10_likelihoods.len());
    let mut log10_post = vec![f64::NEG_INFINITY; n];
    for (gi, &(a, b)) in pairs.iter().enumerate().take(n) {
        let log10_hwe = if a == b {
            2.0 * log10_allele_frequencies[a]
        } else {
            (2.0_f64).log10() + log10_allele_frequencies[a] + log10_allele_frequencies[b]
        };
        log10_post[gi] = log10_likelihoods[gi] + log10_hwe;
    }
    let log10_sum = log10_sum_log10(&log10_post);
    for x in &mut log10_post {
        *x -= log10_sum;
    }
    log10_post
}

fn effective_allele_counts_multi(
    samples_log10_likelihoods: &[&[f64]],
    log10_allele_frequencies: &[f64],
    pairs: &[(usize, usize)],
) -> Vec<f64> {
    let n_alleles = log10_allele_frequencies.len();
    let mut log10_counts = vec![f64::NEG_INFINITY; n_alleles];
    for gl in samples_log10_likelihoods {
        if gl.len() < pairs.len() {
            continue;
        }
        let post = log10_normalized_genotype_posteriors_multi(gl, log10_allele_frequencies, pairs);
        for (gi, &(a, b)) in pairs.iter().enumerate() {
            let log_p = post[gi];
            log10_counts[a] = log10_sum_log10(&[log10_counts[a], log_p]);
            log10_counts[b] = log10_sum_log10(&[log10_counts[b], log_p]);
        }
    }
    log10_counts
        .iter()
        .map(|&c| if c.is_finite() { 10_f64.powf(c) } else { 0.0 })
        .collect()
}

/// Multi-allelic diploid EM AF for GenotypeGVCFs (Dirichlet / HWE mean-field).
pub fn calculate_multiallelic_af_em(
    samples_log10_likelihoods: &[&[f64]],
    n_alleles: usize,
    config: &AfCalculatorConfig,
) -> GatkResult<MultiallelicAfResult> {
    if n_alleles < 2 {
        return Err(gatk_common::GatkError::argument(
            "multiallelic AF requires REF + at least one ALT",
        ));
    }
    if n_alleles == 2 {
        let bi = calculate_biallelic_af_em(samples_log10_likelihoods, config)?;
        let alt_c = bi.alt_allele_count as f64;
        let ref_c = if bi.af > 1e-12 {
            alt_c * (1.0 - bi.af) / bi.af
        } else {
            0.0
        };
        return Ok(MultiallelicAfResult {
            allele_frequencies: vec![1.0 - bi.af, bi.af],
            mle_allele_counts: vec![ref_c, alt_c],
            log10_posterior_no_variant: bi.log10_posterior_no_variant,
            em_iterations: bi.em_iterations,
        });
    }

    let pairs = diploid_genotype_pairs(n_alleles);
    let flat = -(n_alleles as f64).log10();
    let mut log10_af = vec![flat; n_alleles];
    let mut allele_counts = vec![0.0_f64; n_alleles];
    let mut iterations = 0usize;
    const THRESHOLD: f64 = 0.1;
    loop {
        iterations += 1;
        let new_counts =
            effective_allele_counts_multi(samples_log10_likelihoods, &log10_af, &pairs);
        let max_diff = new_counts
            .iter()
            .zip(allele_counts.iter())
            .map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max);
        allele_counts = new_counts;
        if max_diff <= THRESHOLD || iterations > 100 {
            break;
        }
        let mut pseudo: Vec<f64> = allele_counts
            .iter()
            .enumerate()
            .map(|(i, &c)| {
                if i == 0 {
                    config.ref_pseudocount + c
                } else {
                    config.snp_pseudocount + c
                }
            })
            .collect();
        // Keep pseudocounts positive.
        for p in &mut pseudo {
            if *p < 1e-12 {
                *p = 1e-12;
            }
        }
        log10_af = log10_dirichlet_mean_weights(&pseudo);
    }

    let mut log10_p_no_variant = 0.0_f64;
    for gl in samples_log10_likelihoods {
        if gl.len() >= pairs.len() {
            let post = log10_normalized_genotype_posteriors_multi(gl, &log10_af, &pairs);
            if let Some(&p0) = post.first() {
                log10_p_no_variant += p0;
            }
        }
    }
    log10_p_no_variant = log10_p_no_variant.min(0.0);

    let total: f64 = allele_counts.iter().sum();
    let allele_frequencies = if total > 0.0 {
        allele_counts.iter().map(|c| c / total).collect()
    } else {
        vec![1.0 / n_alleles as f64; n_alleles]
    };

    Ok(MultiallelicAfResult {
        allele_frequencies,
        mle_allele_counts: allele_counts,
        log10_posterior_no_variant: log10_p_no_variant,
        em_iterations: iterations,
    })
}

/// Site QUAL (phred) from multi-sample AF posterior P(no variant).
pub fn qual_from_log10_p_no_variant(log10_p_no_variant: f64) -> f64 {
    let clamped = log10_p_no_variant.min(0.0);
    (-10.0 * clamped).max(0.0)
}

#[cfg(test)]
mod six_r40_af_calculator_tests {
    use super::*;

    /// Java `AlleleFrequencyCalculator.calculate` loop (pre-6R.40 production used a
    /// break-before-Dirichlet order). Kept as an independent copy so the dump test can
    /// assert production now matches this order.
    fn calculate_biallelic_af_em_java_loop_order(
        samples_log10_likelihoods: &[&[f64]],
        config: &AfCalculatorConfig,
    ) -> AfCalculationResult {
        const THRESHOLD: f64 = 0.1;
        let flat = -(2.0_f64).log10();
        let mut log10_af = [flat, flat];
        let mut allele_counts = [0.0_f64, 0.0];
        let mut iterations = 0usize;
        let mut allele_counts_maximum_difference = f64::INFINITY;
        while allele_counts_maximum_difference > THRESHOLD {
            iterations += 1;
            let new_counts =
                effective_allele_counts_biallelic(samples_log10_likelihoods, &log10_af);
            allele_counts_maximum_difference = new_counts
                .iter()
                .zip(allele_counts.iter())
                .map(|(a, b)| (a - b).abs())
                .fold(0.0_f64, f64::max);
            allele_counts = new_counts;
            let posterior_pseudo = [
                config.ref_pseudocount + allele_counts[0],
                config.snp_pseudocount + allele_counts[1],
            ];
            let means = log10_dirichlet_mean_weights(&posterior_pseudo);
            log10_af = [means[0], means[1]];
            if iterations > 100 {
                break;
            }
        }
        let mut log10_p_no_variant = 0.0_f64;
        for gl in samples_log10_likelihoods {
            if gl.len() >= 3 {
                let post = log10_normalized_genotype_posteriors_biallelic(gl, &log10_af);
                log10_p_no_variant += post[0];
            }
        }
        log10_p_no_variant = log10_p_no_variant.min(0.0);
        let total_alleles = allele_counts[0] + allele_counts[1];
        let af = if total_alleles > 0.0 {
            allele_counts[1] / total_alleles
        } else {
            0.0
        };
        AfCalculationResult {
            alt_allele_count: allele_counts[1].round() as i32,
            af,
            log10_posterior_no_variant: log10_p_no_variant,
            log10_p_no_alt: log10_p_no_variant,
            em_iterations: iterations,
        }
    }

    fn gls_from_pl_90_6_0() -> [f64; 3] {
        [-9.0, -0.6, 0.0]
    }

    #[test]
    fn six_r40_af_em_on_canonical_pl_dumps_mle_and_qual() {
        let gl = gls_from_pl_90_6_0();
        let cfg = AfCalculatorConfig::default();
        let rust = calculate_biallelic_af_em(&[&gl], &cfg).expect("af");
        let java_loop = calculate_biallelic_af_em_java_loop_order(&[&gl], &cfg);
        let rust_qual = -10.0 * rust.log10_posterior_no_variant;
        let java_loop_qual = -10.0 * java_loop.log10_posterior_no_variant;
        eprintln!(
            "AF_EM rust: mleac={} af={:.6} log10PNoVar={:.8} QUAL={:.6} iters={}",
            rust.alt_allele_count,
            rust.af,
            rust.log10_posterior_no_variant,
            rust_qual,
            rust.em_iterations
        );
        eprintln!(
            "AF_EM java_loop_order: mleac={} af={:.6} log10PNoVar={:.8} QUAL={:.6} iters={}",
            java_loop.alt_allele_count,
            java_loop.af,
            java_loop.log10_posterior_no_variant,
            java_loop_qual,
            java_loop.em_iterations
        );
        eprintln!("Java VCF QUAL=78.32 MLEAC=1; called GT index would give MLEAC=2");
        assert_eq!(rust.alt_allele_count, 1);
        assert!((rust_qual - 78.32).abs() < 0.02, "QUAL rust={rust_qual}");
        assert!((java_loop_qual - 78.32).abs() < 0.02);
        assert!((rust_qual - java_loop_qual).abs() < 1e-6);
    }

    #[test]
    fn six_r40_mleac_from_called_gt_is_not_af_mle() {
        // Java composeCallAttributes: MLEAC = round(EM effective alt count), not the called GT.
        // Production annotate_hc_variant_site now uses AF MLEAC; AC remains called GT (=2).
        let gl = gls_from_pl_90_6_0();
        let af = calculate_biallelic_af_em(&[&gl], &AfCalculatorConfig::default()).expect("af");
        let called_hom_alt_ac = 2i32;
        eprintln!(
            "MLEAC_CONTRACT called_GT_AC={} AF_round_alt={} differ={}",
            called_hom_alt_ac,
            af.alt_allele_count,
            called_hom_alt_ac != af.alt_allele_count
        );
        assert_eq!(
            called_hom_alt_ac, 2,
            "1/1 diploid called AC is 2 (Java AC=2, not MLEAC)"
        );
        // Document: AC from called GT is 2; AF MLE alt count is 1 (Java MLEAC).
        assert_ne!(
            af.alt_allele_count, called_hom_alt_ac,
            "on this PL, AF MLE alt count must not equal called hom-alt AC (Java MLEAC=1)"
        );
        assert_eq!(af.alt_allele_count, 1);
    }
}
