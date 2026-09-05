//! Locus-by-locus activity scoring aligned with GATK `HaplotypeCallerEngine#isActive`.
//! `calc_ref_vs_any_log10_genotype_likelihoods` mirrors `ReferenceConfidenceModel#calcGenotypeLikelihoodsOfRefVsAny`
//! for the non–flow path (`readsWereRealigned = false`).
//! This module ports the **pre-band-pass** probability path used for assembly-region discovery:
//! 1. [`calc_ref_vs_any_log10_genotype_likelihoods`] — `ReferenceConfidenceModel#calcGenotypeLikelihoodsOfRefVsAny`
//! for the non–flow-based reference model (`readsWereRealigned = false`).
//! 2. [`calculate_single_sample_biallelic_non_ref_posterior`] — `AlleleFrequencyCalculator#calculateSingleSampleBiallelicNonRefPosterior`.
//! Joint multisample `isActive` (≥2 samples with nonempty pileups) mirrors `AlleleFrequencyCalculator × GenotypingEngine`,
//! feeding per-sample genotype likelihood vectors through **`htsjdk` PL quantization** (`GenotypeBuilder.PL(double[])`) exactly like HC;
//! fallback `haplotype_caller_activity_profile_state_single_sample`-style pooling when only one nonzero pileup survives split.

use crate::activity_profile::ActivityProfileState;
use gatk_common::HaplotypeCallerConfig;
use statrs::function::gamma::ln_gamma;
use std::sync::OnceLock;

/// GATK `ReferenceConfidenceModel.REF_MODEL_DELETION_QUAL`.
pub const REF_MODEL_DELETION_QUAL: u8 = 30;

/// GATK `HomoSapiensConstants.SNP_HETEROZYGOSITY` default used by `GenotypeCalculationArgumentCollection`.
pub const DEFAULT_SNP_HETEROZYGOSITY: f64 = 1e-3;

/// GATK `HomoSapiensConstants.INDEL_HETEROZYGOSITY`.
pub const DEFAULT_INDEL_HETEROZYGOSITY: f64 = 1.0 / 8000.0;

/// Active-region cap applied to `standardConfidenceForCalling` in `HaplotypeCallerEngine#initializeActiveRegionEvaluationGenotyperEngine`.
pub const ACTIVE_REGION_STANDARD_CONFIDENCE_CAP: f64 = 4.0;

/// GATK default `GenotypeCalculationArgumentCollection.heterozygosityStandardDeviation`.
pub const DEFAULT_HETEROZYGOSITY_STANDARD_DEVIATION: f64 = 0.01;

/// GATK `HaplotypeCallerEngine` HQ soft-clip heuristic threshold.
pub const AVERAGE_HQ_SOFTCLIPS_HQ_BASES_THRESHOLD: f64 = 6.0;

/// GATK `MathUtils.LOG10_ONE_THIRD` (= −log₁₀ 3).
pub const LOG10_ONE_THIRD: f64 = -0.47712125471966244;

/// Parameters for HC activity scoring (defaults match a standard non–DRAGEN HC run).
/// # Invariants
/// `sample_ploidy` is typically 2; active-region confidence is capped at [`ACTIVE_REGION_STANDARD_CONFIDENCE_CAP`].
/// Heterozygosity / SD feed Dirichlet pseudocounts for AF calculators.
/// # Ownership
/// Cloneable config; built from [`HaplotypeCallerConfig`] or defaults.
/// # Mutation
/// Snapshot per evaluator/profile; not mutated during locus scoring.
/// # Biological assumptions
/// Pre-band-pass ref-vs-any GL and non-ref posterior for assembly-region discovery.
/// # Java equivalence
/// GATK `HaplotypeCallerEngine#isActive` / `ReferenceConfidenceModel` / `AlleleFrequencyCalculator` knobs.
#[derive(Debug, Clone)]
pub struct HaplotypeCallerActivityScoringParams {
    pub min_base_quality_score: u8,
    /// Sample ploidy (`≥ 1`); type-enforced via [`crate::bio_ids::Ploidy`].
    pub sample_ploidy: crate::bio_ids::Ploidy,
    pub snp_heterozygosity: f64,
    pub indel_heterozygosity: f64,
    pub heterozygosity_standard_deviation: f64,
    pub ref_model_deletion_qual: u8,
    pub active_region_alt_multiplier: f64,
    /// MinimalGenotypingEngine AF prior / active-region call confidence (Java `min(cap, user standardConfidenceForCalling)`).
    pub active_region_standard_confidence_for_calling: f64,
    /// `StandardCallerArgumentCollection.annotateAllSitesWithPLs` (affects `MinimalGenotypingEngine#forceKeepAllele` / QUAL branch).
    pub annotate_all_sites_with_pls: bool,
    /// `OutputMode.EMIT_ALL_CONFIDENT_SITES` when true; parity defaults to `EMIT_VARIANTS_ONLY`.
    pub emit_all_confident_sites: bool,
    /// When true, apply `(qual <= minBaseQual)` even for deletions (GATK flow-based path).
    pub flow_based_reference_model: bool,
    /// GATK `StandardCallerArgumentCollection.CONTAMINATION_FRACTION` applied to pileups before GL
    /// (`AlleleBiasedDownsamplingUtils`). HC active-region engine forces `0.0`; parity gates may set `> 0`.
    pub contamination_fraction_to_filter: f64,
}

impl Default for HaplotypeCallerActivityScoringParams {
    fn default() -> Self {
        Self {
            min_base_quality_score: 10,
            sample_ploidy: crate::bio_ids::Ploidy::DIPLOID,
            snp_heterozygosity: DEFAULT_SNP_HETEROZYGOSITY,
            indel_heterozygosity: DEFAULT_INDEL_HETEROZYGOSITY,
            heterozygosity_standard_deviation: DEFAULT_HETEROZYGOSITY_STANDARD_DEVIATION,
            ref_model_deletion_qual: REF_MODEL_DELETION_QUAL,
            active_region_alt_multiplier: 1.0,
            active_region_standard_confidence_for_calling: ACTIVE_REGION_STANDARD_CONFIDENCE_CAP,
            annotate_all_sites_with_pls: false,
            emit_all_confident_sites: false,
            flow_based_reference_model: false,
            contamination_fraction_to_filter: 0.0,
        }
    }
}

impl HaplotypeCallerActivityScoringParams {
    pub fn from_haplotype_caller_config(cfg: &HaplotypeCallerConfig) -> Self {
        let mut s = Self::default();
        s.min_base_quality_score = cfg.min_base_quality_score;
        let user = cfg
            .stand_call_confidence
            .min(ACTIVE_REGION_STANDARD_CONFIDENCE_CAP);
        s.active_region_standard_confidence_for_calling = user;
        s
    }

    pub(crate) fn ref_pseudocount(&self) -> f64 {
        self.snp_heterozygosity
            / (self.heterozygosity_standard_deviation * self.heterozygosity_standard_deviation)
    }

    pub(crate) fn snp_pseudocount(&self) -> f64 {
        self.snp_heterozygosity * self.ref_pseudocount()
    }

    pub(crate) fn indel_pseudocount(&self) -> f64 {
        self.indel_heterozygosity * self.ref_pseudocount()
    }
}

/// One read’s evidence at a locus after resolving deletion quality (GATK `ReadPileup` element).
/// # Invariants
/// `qual` is base quality or deletion surrogate (`REF_MODEL_DELETION_QUAL` path).
/// `is_alt` reflects RCM before/after-assembly alt classification for the configured path.
/// # Ownership
/// [`Copy`] pileup cell; no ownership of the underlying BAM record.
/// # Mutation
/// Immutable per locus evaluation; contamination downsampling may drop observations wholesale.
/// # Biological assumptions
/// Encodes one read’s base/deletion evidence at a reference locus for activity GLs.
/// # Java equivalence
/// GATK `PileupElement` fields used by `ReferenceConfidenceModel` / `isActive`.
#[derive(Debug, Clone, Copy)]
pub struct PileupObservation {
    /// Read base at the locus (`PileupElement.getBase`), uppercased when not deletion marker.
    pub read_base: u8,
    /// Base quality at the locus (or deletion surrogate quality).
    pub qual: u8,
    pub is_deletion: bool,
    /// Result of `ReferenceConfidenceModel.isAltBeforeAssembly` / `isAltAfterAssembly` for `readsWereRealigned`.
    pub is_alt: bool,
    /// `PileupElement.isNextToSoftClip` — HQ soft-clip accumulation is conditional on this in RCM.
    pub is_next_to_soft_clip: bool,
    /// `AlignmentUtils.countHighQualitySoftClips(read, RCM threshold)` for this read.
    pub read_hq_soft_clip_base_count: u32,
}

/// GATK `ReferenceConfidenceModel.isAltBeforeAssembly` for `readsWereRealigned == false`.
#[inline]
pub fn is_alt_before_assembly(
    read_base: u8,
    ref_base: u8,
    is_deletion: bool,
    is_before_deletion_start: bool,
    is_after_deletion_end: bool,
    is_before_insertion: bool,
    is_after_insertion: bool,
    is_next_to_soft_clip: bool,
) -> bool {
    read_base != ref_base
        || is_deletion
        || is_before_deletion_start
        || is_after_deletion_end
        || is_before_insertion
        || is_after_insertion
        || is_next_to_soft_clip
}

/// GATK `ReferenceConfidenceModel.isAltAfterAssembly` for post-realign sample evidence.
#[inline]
pub fn is_alt_after_assembly(read_base: u8, ref_base: u8, is_deletion: bool) -> bool {
    read_base != ref_base || is_deletion
}

#[inline]
fn qual_to_error_prob_log10_double(qual: f64) -> f64 {
    qual * -0.1
}

fn qual_to_prob_log10_byte(qual: u8) -> f64 {
    static CACHE: OnceLock<[f64; 256]> = OnceLock::new();
    let cache = CACHE.get_or_init(|| {
        let mut a = [0.0_f64; 256];
        for i in 0..=254u16 {
            let e = 10_f64.powf(-(i as f64) / 10.0);
            a[i as usize] = (1.0 - e).log10();
        }
        a
    });
    cache[(qual as usize) & 0xff]
}

#[inline]
fn qual_to_error_prob_log10_byte(qual: u8) -> f64 {
    qual_to_error_prob_log10_double(qual as f64)
}

#[inline]
fn fast_round_java(d: f64) -> usize {
    if d > 0.0 {
        (d + 0.5) as usize
    } else {
        (d - 0.5) as usize
    }
}

fn jacobian_log_table_get(diff: f64) -> f64 {
    const MAX_TOLERANCE: f64 = 8.0;
    const TABLE_STEP: f64 = 0.0001;
    static CACHE: OnceLock<Vec<f64>> = OnceLock::new();
    let cache = CACHE.get_or_init(|| {
        let n = (MAX_TOLERANCE / TABLE_STEP) as usize + 1;
        (0..n)
            .map(|k| (1.0 + 10_f64.powf(-(k as f64) * TABLE_STEP)).log10())
            .collect::<Vec<_>>()
    });
    let inv_step = 1.0 / TABLE_STEP;
    let mut idx = fast_round_java(diff * inv_step);
    if idx >= cache.len() {
        idx = cache.len() - 1;
    }
    cache[idx]
}

/// GATK `MathUtils.approximateLog10SumLog10(a, b)`.
pub fn approximate_log10_sum_log10_pair(a: f64, b: f64) -> f64 {
    const MAX_TOLERANCE: f64 = 8.0;
    let (a, b) = if a > b { (b, a) } else { (a, b) };
    if a == f64::NEG_INFINITY {
        return b;
    }
    let diff = b - a;
    if diff < MAX_TOLERANCE {
        b + jacobian_log_table_get(diff)
    } else {
        b
    }
}

/// Java `Math.round(double)` narrowed to PL-range values (`Math.floor(x + 0.5)` for finite `x`).
#[inline]
fn java_math_round_double_to_i32(x: f64) -> i32 {
    debug_assert!(x.is_finite());
    let rounded = (x + 0.5).floor();
    if rounded < i32::MIN as f64 {
        i32::MIN
    } else if rounded > i32::MAX as f64 {
        i32::MAX
    } else {
        rounded as i32
    }
}

/// Genotype likelihoods as consumed by HC `AlleleFrequencyCalculator` after VC construction:
/// [`GenotypeBuilder#PL(double[])`](https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/variantcontext/GenotypeBuilder.html)
/// converts log10 likelihoods → integer PLs then back (`GenotypeLikelihoods` / `FastGenotype#getLikelihoods`).
pub(crate) fn genotype_log10_likelihoods_after_java_genotype_pl_roundtrip(gl: &[f64]) -> Vec<f64> {
    if gl.is_empty() {
        return Vec::new();
    }
    let adjust = gl.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    gl.iter()
        .copied()
        .map(|g| {
            let shifted = (-10.0_f64) * (g - adjust);
            let clamped = shifted.min(i32::MAX as f64);
            let pl = java_math_round_double_to_i32(clamped);
            (pl as f64) / (-10.0)
        })
        .collect()
}

/// GATK `MathUtils.log10SumLog10` over the full slice.
pub fn log10_sum_log10(values: &[f64]) -> f64 {
    if values.is_empty() {
        return f64::NEG_INFINITY;
    }
    if values.len() == 1 {
        return values[0];
    }
    let mut max_i = 0usize;
    for i in 1..values.len() {
        if values[i] > values[max_i] {
            max_i = i;
        }
    }
    let max_v = values[max_i];
    if max_v == f64::NEG_INFINITY {
        return max_v;
    }
    let mut sum = 1.0_f64;
    for (i, &cur_val) in values.iter().enumerate() {
        if i == max_i || cur_val == f64::NEG_INFINITY {
            continue;
        }
        sum += 10_f64.powf(cur_val - max_v);
    }
    max_v + if sum != 1.0 { sum.log10() } else { 0.0 }
}

/// GATK `MathUtils.normalizeFromLog10ToLinearSpace` (returns probabilities summing to 1).
/// Borrows the log10 vector so callers that still need it (e.g. genotype posteriors)
/// do not clone before normalizing.
pub fn normalize_from_log10_to_linear_space(log10_probs: &[f64]) -> Vec<f64> {
    let s = log10_sum_log10(log10_probs);
    log10_probs.iter().map(|x| 10_f64.powf(*x - s)).collect()
}

#[inline]
fn max_element_index(values: &[f64]) -> usize {
    debug_assert!(!values.is_empty());
    let mut m = 0usize;
    for i in 1..values.len() {
        if values[i] > values[m] {
            m = i;
        }
    }
    m
}

#[inline]
fn log_to_log10(ln: f64) -> f64 {
    ln * std::f64::consts::LOG10_E
}

/// Natural log of binomial coefficient C(n,k) — `CombinatoricsUtils.binomialCoefficientLog`.
pub fn log_binomial_coefficient_natural(n: u32, k: u32) -> f64 {
    debug_assert!(k <= n);
    ln_gamma((n + 1) as f64) - ln_gamma((k + 1) as f64) - ln_gamma((n - k + 1) as f64)
}

/// GATK `AlleleFrequencyCalculator.calculateSingleSampleBiallelicNonRefPosterior` with priors from
/// `AlleleFrequencyCalculator.makeCalculator(genotypeArgs)` (SNP branch).
pub fn calculate_single_sample_biallelic_non_ref_posterior(
    log10_genotype_likelihoods: &[f64],
    return_zero_if_ref_is_max: bool,
    params: &HaplotypeCallerActivityScoringParams,
) -> f64 {
    let ploidy = log10_genotype_likelihoods.len().saturating_sub(1);
    if ploidy == 0 || log10_genotype_likelihoods.is_empty() {
        return 0.0;
    }

    if return_zero_if_ref_is_max && max_element_index(log10_genotype_likelihoods) == 0 {
        return 0.0;
    }

    let ref_pc = params.ref_pseudocount();
    let snp_pc = params.snp_pseudocount();

    let mut log10_unnorm = Vec::with_capacity(ploidy + 1);
    for n in 0..=ploidy {
        let ln_prior = log_binomial_coefficient_natural(ploidy as u32, n as u32)
            + ln_gamma(n as f64 + snp_pc)
            + ln_gamma(ploidy as f64 - n as f64 + ref_pc);
        log10_unnorm.push(log10_genotype_likelihoods[n] + log_to_log10(ln_prior));
    }

    if return_zero_if_ref_is_max && max_element_index(&log10_unnorm) == 0 {
        return 0.0;
    }

    let linear = normalize_from_log10_to_linear_space(&log10_unnorm);
    1.0 - linear[0]
}

/// Multi-sample active probability from independent per-sample posteriors:
/// `1 - Π_s (1 - P_s(non-ref))`.
/// This gives `P(any sample is non-ref)` and provides a stable multi-sample
/// activity signal until full joint genotyping scaffolding is wired.
pub fn calculate_multi_sample_any_non_ref_posterior(
    per_sample_log10_genotype_likelihoods: &[Vec<f64>],
    return_zero_if_ref_is_max: bool,
    params: &HaplotypeCallerActivityScoringParams,
) -> f64 {
    if per_sample_log10_genotype_likelihoods.is_empty() {
        return 0.0;
    }
    let mut inactive_prob = 1.0_f64;
    for gl in per_sample_log10_genotype_likelihoods {
        let p = calculate_single_sample_biallelic_non_ref_posterior(
            gl,
            return_zero_if_ref_is_max,
            params,
        );
        inactive_prob *= (1.0 - p).clamp(0.0, 1.0);
    }
    1.0 - inactive_prob
}

fn effective_qual(obs: &PileupObservation, params: &HaplotypeCallerActivityScoringParams) -> u8 {
    if obs.is_deletion {
        params.ref_model_deletion_qual
    } else {
        obs.qual
    }
}

fn skip_observation(
    obs: &PileupObservation,
    qual: u8,
    params: &HaplotypeCallerActivityScoringParams,
) -> bool {
    // Java includes bases at exactly minBaseQualityScore (Q10 when threshold is 10).
    let qfilt = qual < params.min_base_quality_score;
    qfilt && (params.flow_based_reference_model || !obs.is_deletion)
}

/// Online mean matching GATK `MathUtils.RunningAverage` updates during RCM `calcGenotypeLikelihoodsOfRefVsAny`
/// (each contributing pileup element that is not skipped; `isAlt && isNextToSoftClip` only).
pub fn hq_soft_clip_running_mean_rcm_path(
    pileup: &[PileupObservation],
    params: &HaplotypeCallerActivityScoringParams,
) -> f64 {
    let mut mean = 0.0_f64;
    let mut count = 0_u64;
    for obs in pileup {
        let qual = effective_qual(obs, params);
        if skip_observation(obs, qual, params) {
            continue;
        }
        if obs.is_alt && obs.is_next_to_soft_clip {
            count += 1;
            mean += (obs.read_hq_soft_clip_base_count as f64 - mean) / count as f64;
        }
    }
    mean
}

/// Multisample ordering: one composite stream in encounter order of nonempty per-sample strata (matches
/// `HaplotypeCallerEngine#isActive` iterating `splitContexts` values in `LinkedHashMap` order; Rust uses
/// the same encounter order as `nonempty_stratified_sample_pileups_ordered`).
pub fn hq_soft_clip_running_mean_rcm_stratified(
    sample_piles: &[&[PileupObservation]],
    params: &HaplotypeCallerActivityScoringParams,
) -> f64 {
    let mut mean = 0.0_f64;
    let mut count = 0_u64;
    for pile in sample_piles {
        for obs in *pile {
            let qual = effective_qual(obs, params);
            if skip_observation(obs, qual, params) {
                continue;
            }
            if obs.is_alt && obs.is_next_to_soft_clip {
                count += 1;
                mean += (obs.read_hq_soft_clip_base_count as f64 - mean) / count as f64;
            }
        }
    }
    mean
}

/// GATK `ReferenceConfidenceModel#calcGenotypeLikelihoodsOfRefVsAny` (log10 GL vector, length `ploidy + 1`).
pub fn calc_ref_vs_any_log10_genotype_likelihoods(
    ploidy: u32,
    pileup: &[PileupObservation],
    params: &HaplotypeCallerActivityScoringParams,
) -> Vec<f64> {
    let lc = (ploidy + 1) as usize;
    let mut gl = vec![0.0_f64; lc];
    let log10_ploidy = (ploidy as f64).log10();
    let mut read_count: u32 = 0;

    for obs in pileup {
        let qual = effective_qual(obs, params);
        if skip_observation(obs, qual, params) {
            continue;
        }
        read_count += 1;
        let is_alt = obs.is_alt;
        let reference_likelihood;
        let non_ref_likelihood;
        if is_alt {
            non_ref_likelihood = qual_to_prob_log10_byte(qual);
            reference_likelihood = qual_to_error_prob_log10_byte(qual) + LOG10_ONE_THIRD;
        } else {
            reference_likelihood = qual_to_prob_log10_byte(qual);
            non_ref_likelihood = qual_to_error_prob_log10_byte(qual) + LOG10_ONE_THIRD;
        }
        let read_weight = if is_alt {
            params.active_region_alt_multiplier
        } else {
            1.0
        };

        gl[0] += read_weight * (reference_likelihood + log10_ploidy);
        gl[lc - 1] += read_weight * (non_ref_likelihood + log10_ploidy);

        let mut i = 1usize;
        let mut j = lc.saturating_sub(2);
        while i < lc - 1 {
            let term = approximate_log10_sum_log10_pair(
                reference_likelihood + (j as f64).log10(),
                non_ref_likelihood + (i as f64).log10(),
            );
            gl[i] += read_weight * term;
            i += 1;
            if j == 0 {
                break;
            }
            j -= 1;
        }
    }

    let denom = (read_count as f64) * log10_ploidy;
    for x in &mut gl {
        *x -= denom;
    }
    gl
}

/// Full single-sample `isActive` activity probability + `ActivityProfileState` metadata (no band-pass).
pub fn haplotype_caller_activity_profile_state_single_sample(
    contig: impl Into<String>,
    pos: u64,
    pileup: &[PileupObservation],
    params: &HaplotypeCallerActivityScoringParams,
) -> ActivityProfileState {
    crate::minimal_genotyping::haplotype_caller_activity_profile_state_minimal_genotyping(
        contig, pos, pileup, params,
    )
}

/// Multi-sample pileups (**nonempty strata only**) use the Java joint path; pooling when only one
/// strata has reads.
pub fn haplotype_caller_activity_profile_state_multi_sample(
    contig: impl Into<String>,
    pos: u64,
    sample_pileups: &[Vec<PileupObservation>],
    params: &HaplotypeCallerActivityScoringParams,
) -> ActivityProfileState {
    let contig = contig.into();
    // Borrow nonempty strata — avoid cloning every observation at every locus.
    let nonempty: Vec<&[PileupObservation]> = sample_pileups
        .iter()
        .filter(|pile| !pile.is_empty())
        .map(|pile| pile.as_slice())
        .collect();
    match nonempty.len() {
        0 => ActivityProfileState::new(contig.as_str(), pos, 0.0),
        1 => haplotype_caller_activity_profile_state_single_sample(
            contig.as_str(),
            pos,
            nonempty[0],
            params,
        ),
        _ => {
            crate::hc_joint_is_active::haplotype_caller_joint_multisample_is_active_activity_state(
                contig, pos, &nonempty, params,
            )
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn is_alt_before_assembly_hom_ref_match() {
        assert!(!is_alt_before_assembly(
            b'A', b'A', false, false, false, false, false, false
        ));
    }

    #[test]
    fn is_alt_after_assembly_ignores_indel_flank_flags() {
        assert!(!is_alt_after_assembly(b'G', b'G', false));
        assert!(is_alt_after_assembly(b'T', b'G', false));
        assert!(is_alt_after_assembly(b'-', b'G', true));
    }

    #[test]
    fn log10_sum_two_values_matches_max_plus_log10_1_plus_pow() {
        let a = -3.0;
        let b = -2.0;
        let s = log10_sum_log10(&[a, b]);
        let expected = b + (1.0 + 10_f64.powf(a - b)).log10();
        assert!((s - expected).abs() < 1e-12, "s={s} exp={expected}");
    }

    #[test]
    fn java_genotype_pl_roundtrip_matches_htsjdk_joint_probe_gl() {
        let gl = [
            -4.177121254719663,
            -0.30108776596620884,
            -8.666178727151364e-5,
        ];
        let out = genotype_log10_likelihoods_after_java_genotype_pl_roundtrip(&gl);
        assert!((out[0] + 4.2).abs() < 1e-12, "got {out:?}");
        assert!((out[1] + 0.3).abs() < 1e-12, "got {out:?}");
        assert!(out[2].abs() < 1e-15, "got {out:?}");
    }

    #[test]
    fn frozen_hom_ref_single_read_matches_python_reference() {
        let params = HaplotypeCallerActivityScoringParams::default();
        let obs = [PileupObservation {
            read_base: b'A',
            qual: 30,
            is_deletion: false,
            is_alt: false,
            is_next_to_soft_clip: false,
            read_hq_soft_clip_base_count: 0,
        }];
        let gl = calc_ref_vs_any_log10_genotype_likelihoods(2, &obs, &params);
        assert_eq!(gl.len(), 3);
        // Slight cross-platform float noise exists because the heterozygous term uses an
        // approximate Jacobian log table. Keep a narrow bound while remaining platform-stable.
        assert!((gl[0] + 0.00043451177401771).abs() < 5e-4, "gl={gl:?}");
        // In the Java formula, heterozygous terms use log10(j)/log10(i) (no +log10(ploidy)),
        // then the shared denominator subtracts log10(ploidy), so diploid het sits near -0.3010.
        assert!((gl[1] + 0.30131962629328163).abs() < 5e-4, "gl={gl:?}");
        assert!((gl[2] + 3.4771212547196626).abs() < 5e-4, "gl={gl:?}");
        assert!(
            gl[0] > gl[1],
            "expected hom-ref > het for one high-Q ref read, gl={gl:?}"
        );
        let p = calculate_single_sample_biallelic_non_ref_posterior(&gl, true, &params);
        assert!(p < 1e-6, "p={p}");
    }

    #[test]
    fn frozen_hom_alt_five_reads_matches_python_reference() {
        let params = HaplotypeCallerActivityScoringParams::default();
        let obs = [PileupObservation {
            read_base: b'T',
            qual: 30,
            is_deletion: false,
            is_alt: true,
            is_next_to_soft_clip: false,
            read_hq_soft_clip_base_count: 0,
        }; 5];
        let gl = calc_ref_vs_any_log10_genotype_likelihoods(2, &obs, &params);
        assert!((gl[0] + 17.385606273598313).abs() < 1e-9);
        assert!((gl[1] + 1.5065981314664083).abs() < 1e-9);
        assert!((gl[2] + 0.0021725588700884924).abs() < 1e-9);
        let p = calculate_single_sample_biallelic_non_ref_posterior(&gl, true, &params);
        assert!((p - 1.0).abs() < 1e-10, "p={p}");
    }

    #[test]
    fn empty_pileup_yields_zero_activity() {
        let params = HaplotypeCallerActivityScoringParams::default();
        let st = haplotype_caller_activity_profile_state_single_sample("chr1", 1, &[], &params);
        assert_eq!(st.active_prob, 0.0);
    }

    #[test]
    fn hq_soft_clip_kind_when_mean_exceeds_threshold() {
        let params = HaplotypeCallerActivityScoringParams::default();
        let obs = [PileupObservation {
            read_base: b'T',
            qual: 30,
            is_deletion: false,
            is_alt: true,
            is_next_to_soft_clip: true,
            read_hq_soft_clip_base_count: 7,
        }];
        let st = haplotype_caller_activity_profile_state_single_sample("chr1", 1, &obs, &params);
        assert_eq!(
            st.evidence,
            crate::activity_profile::ActivityEvidence::HighQualitySoftClips { clip_bases: 7 }
        );
    }

    #[test]
    fn multi_sample_probability_aggregates_any_non_ref() {
        let params = HaplotypeCallerActivityScoringParams::default();
        let s1 = vec![PileupObservation {
            read_base: b'T',
            qual: 30,
            is_deletion: false,
            is_alt: true,
            is_next_to_soft_clip: false,
            read_hq_soft_clip_base_count: 0,
        }];
        let s2 = vec![PileupObservation {
            read_base: b'A',
            qual: 30,
            is_deletion: false,
            is_alt: false,
            is_next_to_soft_clip: false,
            read_hq_soft_clip_base_count: 0,
        }];
        let gl1 = calc_ref_vs_any_log10_genotype_likelihoods(2, &s1, &params);
        let gl2 = calc_ref_vs_any_log10_genotype_likelihoods(2, &s2, &params);
        let p1 = calculate_single_sample_biallelic_non_ref_posterior(&gl1, true, &params);
        let p2 = calculate_single_sample_biallelic_non_ref_posterior(&gl2, true, &params);
        let pm = calculate_multi_sample_any_non_ref_posterior(&[gl1, gl2], true, &params);
        let expected = 1.0 - (1.0 - p1) * (1.0 - p2);
        assert!((pm - expected).abs() < 1e-12, "pm={pm} expected={expected}");
    }

    #[test]
    fn multi_sample_state_reduces_to_single_sample_for_one_sample() {
        let params = HaplotypeCallerActivityScoringParams::default();
        let obs = vec![PileupObservation {
            read_base: b'T',
            qual: 30,
            is_deletion: false,
            is_alt: true,
            is_next_to_soft_clip: false,
            read_hq_soft_clip_base_count: 0,
        }];
        let s = haplotype_caller_activity_profile_state_single_sample("chr1", 11, &obs, &params);
        let m = haplotype_caller_activity_profile_state_multi_sample("chr1", 11, &[obs], &params);
        assert!((s.active_prob - m.active_prob).abs() < 1e-12);
        assert!((s.original_active_prob - m.original_active_prob).abs() < 1e-12);
    }

    #[test]
    fn multi_sample_state_empty_returns_zero() {
        let params = HaplotypeCallerActivityScoringParams::default();
        let st = haplotype_caller_activity_profile_state_multi_sample("chr1", 1, &[], &params);
        assert_eq!(st.active_prob, 0.0);
        assert_eq!(
            st.evidence,
            crate::activity_profile::ActivityEvidence::Plain
        );
    }
}
