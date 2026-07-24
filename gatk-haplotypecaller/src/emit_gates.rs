//! Shared HC variant **emit gates** (Java AF threshold, read support, genotype emit parity).
//! Leaf module breaking the `hc_emit_policy` ↔ `hc_genotyping_engine` cycle.

use crate::activity_scoring::genotype_log10_likelihoods_after_java_genotype_pl_roundtrip;
use crate::af_calc::{calculate_biallelic_af_em, AfCalculatorConfig};
use crate::bio_ids::AlleleDepth;
use crate::compatibility::{is_coupled_indel_for_genotyping, is_ctc_del_for_genotyping};
use crate::event_map::VariationEvent;
use crate::genotyping::{
    best_biallelic_diploid_genotype_index, biallelic_genotype_index_from_pl, GenotypeFormatFields,
};
use crate::hc_joint_is_active::log10_one_minus_pow10;
use crate::java_hc_site_semantics::is_cluster_anchor_snp;
use gatk_common::GatkResult;

pub(crate) const AFC_EMIT_EPSILON: f64 = 1e-9;

pub(crate) fn qual_to_error_prob_log10(phred_threshold: f64) -> f64 {
    -(phred_threshold / 10.0)
}

/// Minimum alt AD to emit hom-alt (cuts weak 1/1 rust-only).
pub const MIN_HOM_ALT_AD_FOR_EMIT: AlleleDepth = AlleleDepth::new(2);

/// Java sparse-BAM: single alt read sufficient at cluster anchor SNPs.
/// `alt_ad` / `ref_ad` are non-negative pileup depths (VCF FORMAT AD wire values).
pub fn passes_cluster_anchor_read_support(alt_ad: i32, ref_ad: i32) -> bool {
    alt_ad >= 1 && alt_ad >= ref_ad
}

/// Read-style VCF emit: cluster anchors @ AD≥1; other SNPs need ≥2 alt reads (cuts rust-only).
/// `alt_ad` / `ref_ad` are non-negative pileup depths (VCF FORMAT AD wire values).
pub fn passes_read_style_sparse_emit(event: &VariationEvent, alt_ad: i32, ref_ad: i32) -> bool {
    if event.ref_allele.len() != 1 || event.alt_allele.len() != 1 {
        return false;
    }
    if is_cluster_anchor_snp(event) {
        return passes_cluster_anchor_read_support(alt_ad, ref_ad);
    }
    alt_ad >= 2 && alt_ad >= ref_ad
}

/// GL vector fed to `AlleleFrequencyCalculator` after Java `GenotypeBuilder.PL` round-trip.
pub fn gl_for_java_af_calculation(genotype_log10_likelihoods: &[f64]) -> Vec<f64> {
    genotype_log10_likelihoods_after_java_genotype_pl_roundtrip(genotype_log10_likelihoods)
}

pub(crate) fn passes_hc_variant_emit_biallelic_inner(
    genotype_log10_likelihoods: &[f64],
    stand_emit_confidence: f64,
) -> GatkResult<bool> {
    if genotype_log10_likelihoods.len() < 3 {
        return Ok(false);
    }
    let gl_java = gl_for_java_af_calculation(genotype_log10_likelihoods);
    let af = calculate_biallelic_af_em(&[&gl_java], &AfCalculatorConfig::default())?;
    let call_conf_log10 = qual_to_error_prob_log10(stand_emit_confidence);
    let alt_plausible = af.log10_posterior_no_variant + AFC_EMIT_EPSILON < call_conf_log10;
    let site_is_monomorphic = !alt_plausible;
    let log10_vc_confidence = if !site_is_monomorphic {
        af.log10_posterior_no_variant
    } else {
        log10_one_minus_pow10(af.log10_posterior_no_variant)
    };
    let phred_scaled = (-10.0 * log10_vc_confidence).max(0.0);
    Ok(!site_is_monomorphic && phred_scaled >= stand_emit_confidence)
}

/// GATK `GenotypingEngine.passesEmitThreshold` for default `EMIT_VARIANTS_ONLY` (single sample).
pub fn passes_hc_variant_emit_biallelic(
    genotype_log10_likelihoods: &[f64],
    stand_emit_confidence: f64,
) -> GatkResult<bool> {
    passes_hc_variant_emit_biallelic_inner(genotype_log10_likelihoods, stand_emit_confidence)
}

/// Java `GenotypingEngine.passesEmitThreshold`: emit when `!bestGuessIsRef` (after PL round-trip).
pub fn passes_java_emit_not_hom_ref(
    genotype_log10_likelihoods: &[f64],
    format: &GenotypeFormatFields,
) -> bool {
    let gl_rt = gl_for_java_af_calculation(genotype_log10_likelihoods);
    let ad = format.ad_as_i32();
    best_biallelic_diploid_genotype_index(&gl_rt, &ad) != 0
}

/// Whether Java `calculateGenotypes` would return non-null (`passesEmitThreshold` on site AFC).
/// **L14-D1:** coupled/CTC use [`is_coupled_indel_for_genotyping`]
/// [`is_ctc_del_for_genotyping`] with `region_events` (phenotype when partners present;
/// empty slice → absolute oracle for fixtures only).
pub fn java_emit_would_pass(
    event: &VariationEvent,
    genotype_log10_likelihoods: &[f64],
    format: &GenotypeFormatFields,
    stand_emit_confidence: f64,
    region_events: &[VariationEvent],
) -> GatkResult<bool> {
    if is_coupled_indel_for_genotyping(event, region_events)
        || is_ctc_del_for_genotyping(event, region_events)
    {
        let best = biallelic_genotype_index_from_pl(&format.pl).as_usize();
        return Ok(best != 0 || format.ad.get(1).copied().map(|d| d.as_i32()).unwrap_or(0) >= 1);
    }
    if !passes_java_emit_not_hom_ref(genotype_log10_likelihoods, format) {
        return Ok(false);
    }
    let gl = gl_for_java_af_calculation(genotype_log10_likelihoods);
    passes_hc_variant_emit_biallelic(&gl, stand_emit_confidence)
}

/// Genotyping + VCF emit confidence (biallelic diploid).
/// **L14-D1:** coupled/CTC recognition threads `region_events` (same contract as
/// [`java_emit_would_pass`]).
pub fn passes_emit_for_variation_event(
    event: &VariationEvent,
    genotype_log10_likelihoods: &[f64],
    format: &GenotypeFormatFields,
    stand_emit_confidence: f64,
    region_events: &[VariationEvent],
) -> GatkResult<bool> {
    if passes_hc_variant_emit_biallelic(genotype_log10_likelihoods, stand_emit_confidence)? {
        return Ok(true);
    }
    // P12 cluster: emit coupled indels when genotype is not hom-ref (Java calls both sites).
    if is_coupled_indel_for_genotyping(event, region_events)
        || is_ctc_del_for_genotyping(event, region_events)
    {
        let best = biallelic_genotype_index_from_pl(&format.pl).as_usize();
        if best != 0 {
            return Ok(true);
        }
        let alt_ad = format.ad.get(1).copied().map(|d| d.as_i32()).unwrap_or(0);
        if alt_ad >= 1 {
            return Ok(true);
        }
    }
    // Legacy N1 bridge only (parity uses Java AF threshold exclusively).
    if is_cluster_anchor_snp(event) {
        let ref_ad = format.ad.first().copied().map(|d| d.as_i32()).unwrap_or(0);
        let alt_ad = format.ad.get(1).copied().map(|d| d.as_i32()).unwrap_or(0);
        if passes_cluster_anchor_read_support(alt_ad, ref_ad) {
            return Ok(true);
        }
    }
    Ok(false)
}
