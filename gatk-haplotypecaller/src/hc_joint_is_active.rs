//! Multisample HaplotypeCaller `isActive` (`HaplotypeCallerEngine`: joint `MinimalGenotypingEngine#calculateGenotypes`,
//! AlleleFrequencyCalculator EM, QUAL → [`QualityUtils#qualToProb`](https://github.com/broadinstitute/gatk)).
//! Per-sample genotype likelihoods are quantized like `htsjdk` `GenotypeBuilder#PL(double[])` before AF/QUAL, matching VC round-tripping in GATK.
//! Assumes stratified nonempty per-sample piles (Java `AlignmentContext.splitContextBySampleName`).

use crate::activity_profile::ActivityProfileState;
use crate::activity_scoring::{
    calc_ref_vs_any_log10_genotype_likelihoods,
    genotype_log10_likelihoods_after_java_genotype_pl_roundtrip, log10_sum_log10,
    HaplotypeCallerActivityScoringParams, PileupObservation,
};
use crate::minimal_genotyping::{
    haplotype_caller_activity_profile_state_minimal_genotyping, max_element_index,
};
const JAVA_EM_THRESHOLD: f64 = 0.1;
/// `AFCalculationResult.EPSILON`.
const AFC_RESULT_EPSILON: f64 = 1.0e-10;

/// Branch threshold for `NaturalLogUtils.log1mexp` (= `Math.log(0.5)` = `-LN_2`).
const LOG1MEXP_THRESHOLD: f64 = -std::f64::consts::LN_2;

const LOG10_TWO: f64 = std::f64::consts::LOG10_2;

#[inline]
fn log1mexp(a: f64) -> f64 {
    if a > 0.0 {
        f64::NAN
    } else if a == 0.0 {
        f64::NEG_INFINITY
    } else if a < LOG1MEXP_THRESHOLD {
        (1.0 - a.exp()).ln()
    } else {
        (-a.exp_m1()).ln()
    }
}

/// GATK `MathUtils.log10OneMinusPow10`.
pub fn log10_one_minus_pow10(a: f64) -> f64 {
    if a > 0.0 {
        f64::NAN
    } else if a == 0.0 {
        f64::NEG_INFINITY
    } else {
        let b = a * std::f64::consts::LN_10;
        log1mexp(b) * std::f64::consts::LOG10_E
    }
}

#[inline]
fn qual_to_error_prob_log10(phred_threshold: f64) -> f64 {
    phred_threshold * -0.1
}

#[inline]
fn qual_to_prob_from_nonnegative_qual(qual: f64) -> f64 {
    debug_assert!(qual >= 0.0, "qual={qual}");
    1.0 - 10_f64.powf(-qual / 10.0)
}

fn dirichlet_log10_mean_weights(alpha: &[f64]) -> Vec<f64> {
    let sum: f64 = alpha.iter().sum();
    alpha.iter().map(|a| (a / sum).log10()).collect()
}

/// Heterozygous genotype log₁₀ combination multiplier (`MathUtils.log10Factorial(ploidy) − ∑log10(count!)`).
#[inline]
fn log10_combo_diploid_heterozygous() -> f64 {
    LOG10_TWO
}

fn normalize_log10_in_place(vals: &mut [f64]) {
    let s = log10_sum_log10(vals);
    for x in vals.iter_mut() {
        *x -= s;
    }
}

fn genotype_log10_posteriors_diploid_biallelic_sample(
    log10_gl: &[f64; 3],
    log10_af: &[f64; 2],
) -> [f64; 3] {
    let combos = [0.0_f64, log10_combo_diploid_heterozygous(), 0.0];
    let mut v = [
        combos[0] + log10_gl[0] + 2.0 * log10_af[0],
        combos[1] + log10_gl[1] + log10_af[0] + log10_af[1],
        combos[2] + log10_gl[2] + 2.0 * log10_af[1],
    ];
    normalize_log10_in_place(&mut v);
    v
}

fn effective_diploid_biallelic_allele_counts(
    per_sample_gl: &[Vec<f64>],
    log10_af: &[f64; 2],
) -> Option<[f64; 2]> {
    let mut acc = [f64::NEG_INFINITY, f64::NEG_INFINITY];
    for gl in per_sample_gl {
        let arr: &[f64; 3] = gl.as_slice().try_into().ok()?;
        let post = genotype_log10_posteriors_diploid_biallelic_sample(arr, log10_af);
        // `post[g]` is log₁₀ P(g|θ); allele copy mass is P(g)×count ⇒ log₁₀(P×c)=post[g]+log₁₀(c).
        // `log10_sum_log10(post[g], log10(c))` would mean log₁₀(10^post+10^log c) and does **not** match Java `AF`.
        acc[0] = log10_sum_log10(&[acc[0], post[0] + LOG10_TWO]);
        acc[1] = log10_sum_log10(&[acc[1], post[2] + LOG10_TWO]);
        acc[0] = log10_sum_log10(&[acc[0], post[1]]);
        acc[1] = log10_sum_log10(&[acc[1], post[1]]);
    }
    Some([10_f64.powf(acc[0]), 10_f64.powf(acc[1])])
}

/// Allele-frequency EM + `calculateGenotypes` QUAL path (biallelic diploid, no spanning deletion).
pub fn haplotype_caller_joint_multisample_is_active_activity_state(
    contig: impl Into<String>,
    pos: u64,
    nonempty_sample_piles: &[&[PileupObservation]],
    scoring: &HaplotypeCallerActivityScoringParams,
) -> ActivityProfileState {
    let contig = contig.into();
    if nonempty_sample_piles.len() < 2 {
        let merged: Vec<_> = nonempty_sample_piles
            .iter()
            .flat_map(|p| p.iter().copied())
            .collect();
        return haplotype_caller_activity_profile_state_minimal_genotyping(
            contig, pos, &merged, scoring,
        );
    }

    let hq_soft_clip_running_mean =
        crate::activity_scoring::hq_soft_clip_running_mean_rcm_stratified(
            nonempty_sample_piles,
            scoring,
        );

    let mut per_sample_gl = Vec::with_capacity(nonempty_sample_piles.len());
    for p in nonempty_sample_piles {
        per_sample_gl.push(calc_ref_vs_any_log10_genotype_likelihoods(
            scoring.sample_ploidy.as_u32(),
            p,
            scoring,
        ));
    }

    let per_sample_java_pl: Vec<Vec<f64>> = per_sample_gl
        .iter()
        .map(|gl| genotype_log10_likelihoods_after_java_genotype_pl_roundtrip(gl))
        .collect();

    let active_prob =
        if scoring.sample_ploidy.get() == 2 && per_sample_java_pl.iter().all(|gl| gl.len() == 3) {
            joint_active_prob_via_af_em_and_qual(&per_sample_java_pl, scoring).unwrap_or(0.0)
        } else {
            crate::activity_scoring::calculate_multi_sample_any_non_ref_posterior(
                &per_sample_java_pl,
                true,
                scoring,
            )
        };

    let original = per_sample_gl
        .iter()
        .map(|gl| {
            let max_i = max_element_index(gl);
            gl[max_i] - gl[0]
        })
        .fold(f64::NEG_INFINITY, f64::max);

    activity_state_meta(
        contig,
        pos,
        active_prob,
        original,
        hq_soft_clip_running_mean,
    )
}

fn activity_state_meta(
    contig: String,
    pos: u64,
    active_prob: f64,
    original: f64,
    hq_soft_clip_running_mean: f64,
) -> ActivityProfileState {
    let evidence = if hq_soft_clip_running_mean
        > crate::activity_scoring::AVERAGE_HQ_SOFTCLIPS_HQ_BASES_THRESHOLD
    {
        crate::activity_profile::ActivityEvidence::HighQualitySoftClips {
            clip_bases: hq_soft_clip_running_mean as u32,
        }
    } else {
        crate::activity_profile::ActivityEvidence::Plain
    };

    ActivityProfileState {
        contig: std::sync::Arc::from(contig),
        pos,
        active_prob,
        original_active_prob: original,
        evidence,
    }
}

/// `per_sample_gl` must match `htsjdk GenotypeLikelihoods` after VC construction (`PL(double[])` quantization).
fn joint_active_prob_via_af_em_and_qual(
    per_sample_gl: &[Vec<f64>],
    scoring: &HaplotypeCallerActivityScoringParams,
) -> Option<f64> {
    // `HaplotypeCallerEngine` uses `FAKE_ALT = <FAKE_ALT>` (symbolic, `Allele#length` 0) vs `FAKE_REF = N`
    // (length 1). `AlleleFrequencyCalculator` therefore applies `indelPseudocount` to the alt, not `snpPseudocount`.
    let ref_pc = scoring.ref_pseudocount();
    let alt_pc_prior = scoring.indel_heterozygosity * ref_pc;
    let prior_counts = [ref_pc, alt_pc_prior];

    let n_allele = 2_usize;
    let flat_freq = -(n_allele as f64).log10();

    let mut log10_af = [flat_freq, flat_freq];
    let mut allele_counts = [0.0_f64; 2];

    for _ in 0..256 {
        let new_counts = effective_diploid_biallelic_allele_counts(per_sample_gl, &log10_af)?;
        let max_diff = (new_counts[0] - allele_counts[0])
            .abs()
            .max((new_counts[1] - allele_counts[1]).abs());
        allele_counts = new_counts;
        let posterior_pc = [
            prior_counts[0] + allele_counts[0],
            prior_counts[1] + allele_counts[1],
        ];
        log10_af = {
            let w = dirichlet_log10_mean_weights(&posterior_pc);
            [w[0], w[1]]
        };
        if max_diff <= JAVA_EM_THRESHOLD {
            break;
        }
    }

    let mut log10_p_no_variant = 0.0_f64;
    for gl in per_sample_gl {
        let arr: &[f64; 3] = gl.as_slice().try_into().ok()?;
        let post = genotype_log10_posteriors_diploid_biallelic_sample(arr, &log10_af);
        // GATK uses `+=` over per-sample log10 P(hom-ref | θ), i.e. log10 ∏_s p_s.
        log10_p_no_variant += post[0];
    }

    let call_conf = scoring.active_region_standard_confidence_for_calling;
    let alt_plausible =
        log10_p_no_variant + AFC_RESULT_EPSILON < qual_to_error_prob_log10(call_conf);
    let site_is_monomorphic = !alt_plausible;

    let log10_vc_confidence = vc_log10_p_error_assignment(
        site_is_monomorphic,
        scoring.annotate_all_sites_with_pls,
        log10_p_no_variant,
    );
    let phred_scaled_qual = (-10.0_f64) * log10_vc_confidence;

    let emits = gc_passes_emit_threshold(
        phred_scaled_qual,
        site_is_monomorphic,
        scoring.emit_all_confident_sites,
        scoring.active_region_standard_confidence_for_calling,
    ) && supplemental_emit_filters();

    if !emits {
        return Some(0.0);
    }

    Some(qual_to_prob_from_nonnegative_qual(phred_scaled_qual))
}

/// Remaining unconditional `calculateGenotypes` veto paths are no-ops for HC IsActive probes (no spanning del, allele cap).
#[inline]
fn supplemental_emit_filters() -> bool {
    true
}

/// `GenotypingEngine#calculateGenotypes` veto that requires `calculateGenotypes` continuation (false here).
/// Mirrors `emitAllActiveSites` / GVCF quirks absent in active-region evaluator.
#[inline]
fn vc_log10_p_error_assignment(
    site_is_monomorphic: bool,
    annotate_all_pls: bool,
    log10_p_posterior_no_variant: f64,
) -> f64 {
    if !site_is_monomorphic || annotate_all_pls {
        log10_p_posterior_no_variant + 0.0
    } else {
        log10_one_minus_pow10(log10_p_posterior_no_variant) + 0.0
    }
}

fn gc_passes_emit_threshold(
    phred_scaled_confidence: f64,
    site_is_monomorphic: bool,
    emit_all_confident_sites: bool,
    standard_confidence_for_calling: f64,
) -> bool {
    let mode_all_confident = emit_all_confident_sites;
    (mode_all_confident || !site_is_monomorphic)
        && phred_scaled_confidence >= standard_confidence_for_calling
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::activity_scoring::{
        calc_ref_vs_any_log10_genotype_likelihoods,
        genotype_log10_likelihoods_after_java_genotype_pl_roundtrip, is_alt_before_assembly,
        PileupObservation,
    };

    #[test]
    fn frozen_joint_snp10_prob_matches_java_hc_probe() {
        let params = HaplotypeCallerActivityScoringParams::default();
        // Ref fixture `parity/fixtures/reference.fa` chr1 ACGT repeating; base at 10 (1-based) is A (index 9).
        let ra = vec![PileupObservation {
            read_base: b'T',
            qual: 37,
            is_deletion: false,
            is_alt: is_alt_before_assembly(b'T', b'A', false, false, false, false, false, false),
            is_next_to_soft_clip: false,
            read_hq_soft_clip_base_count: 0,
        }];
        let rb = vec![PileupObservation {
            read_base: b'A',
            qual: 37,
            is_deletion: false,
            is_alt: false,
            is_next_to_soft_clip: false,
            read_hq_soft_clip_base_count: 0,
        }];
        let gl_a = genotype_log10_likelihoods_after_java_genotype_pl_roundtrip(
            &calc_ref_vs_any_log10_genotype_likelihoods(2, &ra, &params),
        );
        let gl_b = genotype_log10_likelihoods_after_java_genotype_pl_roundtrip(
            &calc_ref_vs_any_log10_genotype_likelihoods(2, &rb, &params),
        );
        let p = joint_active_prob_via_af_em_and_qual(&[gl_a, gl_b], &params)
            .expect("vc path returns Some");
        assert!(
            (p - 0.99942045).abs() < 1e-5,
            "expected Java raw-activity joint probe (~0.99942045), got {p}"
        );
    }

    #[test]
    fn inactive_all_ref_joint_pair() {
        let params = HaplotypeCallerActivityScoringParams::default();
        let obs = vec![PileupObservation {
            read_base: b'A',
            qual: 30,
            is_deletion: false,
            is_alt: false,
            is_next_to_soft_clip: false,
            read_hq_soft_clip_base_count: 0,
        }];
        let ga = genotype_log10_likelihoods_after_java_genotype_pl_roundtrip(
            &calc_ref_vs_any_log10_genotype_likelihoods(2, &obs, &params),
        );
        let gb = genotype_log10_likelihoods_after_java_genotype_pl_roundtrip(
            &calc_ref_vs_any_log10_genotype_likelihoods(2, &obs, &params),
        );
        let p = joint_active_prob_via_af_em_and_qual(&[ga, gb], &params).expect("Some");
        assert_eq!(p, 0.0);
    }
}
