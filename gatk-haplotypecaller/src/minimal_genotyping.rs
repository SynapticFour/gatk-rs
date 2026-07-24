//! GATK `MinimalGenotypingEngine` path for HC `isActive` (GAP-C-01).
//! Java delegates to `AlleleFrequencyCalculator#calculateSingleSampleBiallelicNonRefPosterior`
//! with `returnZeroIfRefIsMax = true` on PL-rounded likelihoods from uncapped RCM output.
//! Hom-ref capping (`getGenotypeLikelihoodsCappedByHomRefLikelihood`) is for gVCF emission only, not `isActive`.

use crate::activity_profile::ActivityProfileState;
use crate::activity_scoring::{
    calc_ref_vs_any_log10_genotype_likelihoods,
    calculate_single_sample_biallelic_non_ref_posterior,
    genotype_log10_likelihoods_after_java_genotype_pl_roundtrip,
    HaplotypeCallerActivityScoringParams, PileupObservation,
    AVERAGE_HQ_SOFTCLIPS_HQ_BASES_THRESHOLD,
};

/// GATK `RefVsAnyResult#getGenotypeLikelihoodsCappedByHomRefLikelihood`.
pub fn cap_genotype_likelihoods_by_hom_ref(gl: &[f64]) -> Vec<f64> {
    if gl.is_empty() {
        return Vec::new();
    }
    let hom = gl[0];
    gl.iter().map(|&g| g.min(hom)).collect()
}

/// GATK `GenotypingEngine#calculateSingleSampleRefVsAnyActiveStateProfileValue`.
#[inline]
pub fn calculate_single_sample_ref_vs_any_active_state_profile_value(
    log10_genotype_likelihoods: &[f64],
    params: &HaplotypeCallerActivityScoringParams,
) -> f64 {
    calculate_single_sample_biallelic_non_ref_posterior(log10_genotype_likelihoods, true, params)
}

/// Full HC single-sample activity state via the MinimalGenotypingEngine + ReferenceConfidenceModel path.
pub fn haplotype_caller_activity_profile_state_minimal_genotyping(
    contig: impl Into<String>,
    pos: u64,
    pileup: &[PileupObservation],
    params: &HaplotypeCallerActivityScoringParams,
) -> ActivityProfileState {
    let contig = contig.into();
    if pileup.is_empty() {
        return ActivityProfileState::new(contig, pos, 0.0);
    }

    let hq_soft_clip_running_mean =
        crate::activity_scoring::hq_soft_clip_running_mean_rcm_path(pileup, params);

    let gl_raw =
        calc_ref_vs_any_log10_genotype_likelihoods(params.sample_ploidy.as_u32(), pileup, params);
    // GATK `HaplotypeCallerEngine#isActive`: `GenotypeBuilder.PL(uncapped RCM GL)` then posterior on
    // `getLikelihoods` — no hom-ref cap before active-state scoring.
    let gl_pl = genotype_log10_likelihoods_after_java_genotype_pl_roundtrip(&gl_raw);
    let is_active_prob =
        calculate_single_sample_ref_vs_any_active_state_profile_value(&gl_pl, params);

    let max_i = max_element_index(&gl_pl);
    let original = gl_pl[max_i] - gl_pl[0];

    let evidence = if hq_soft_clip_running_mean > AVERAGE_HQ_SOFTCLIPS_HQ_BASES_THRESHOLD {
        crate::activity_profile::ActivityEvidence::HighQualitySoftClips {
            clip_bases: hq_soft_clip_running_mean as u32,
        }
    } else {
        crate::activity_profile::ActivityEvidence::Plain
    };

    ActivityProfileState {
        contig: std::sync::Arc::from(contig),
        pos,
        active_prob: is_active_prob,
        original_active_prob: original,
        evidence,
    }
}

pub(crate) fn max_element_index(values: &[f64]) -> usize {
    debug_assert!(!values.is_empty());
    let mut m = 0usize;
    for i in 1..values.len() {
        if values[i] > values[m] {
            m = i;
        }
    }
    m
}
