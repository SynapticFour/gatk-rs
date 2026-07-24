//! Shared HC variant **emit gates** for production `strict_java`.
//! Combines Java-style AF / confidence checks with optional ASM-8 production admission
//! ([`crate::read_event_discovery::is_strict_java_production_emit_admits`] — band/motif
//! predicates, **not** the `p12_java_only.tsv` oracle; Sprint L-3).
//! Harness-only baseline oracle blocking goes through `p12_baseline_emit_oracle_blocks`
//! (env + `parity_harness`). See `docs/CLAIM_MATRIX.md`.

pub use crate::emit_gates::{
    passes_cluster_anchor_read_support, passes_emit_for_variation_event,
    passes_read_style_sparse_emit, MIN_HOM_ALT_AD_FOR_EMIT,
};

use crate::compatibility::{is_coupled_indel_for_genotyping, is_ctc_del_for_genotyping};
use crate::emit_gates::{
    gl_for_java_af_calculation, java_emit_would_pass, passes_hc_variant_emit_biallelic,
    passes_java_emit_not_hom_ref,
};
use crate::event_map::VariationEvent;
use crate::genotyping::{biallelic_genotype_index_from_pl, GenotypeFormatFields};
use crate::java_hc_site_semantics::is_cluster_anchor_snp;
use crate::read_event_discovery::{
    is_strict_java_production_emit_admits, p12_baseline_emit_oracle_blocks,
    strict_java_asm8_only_enabled,
};
use gatk_common::GatkResult;
use gatk_core::io::vcf::Genotype;

/// Strict Java emit — `GenotypingEngine.calculateGenotypes` + `passesEmitThreshold` (site AFC only).
/// **L14-D1:** pass non-empty `region_events` on production (phenotype for coupled/CTC).
pub fn passes_strict_java_emit_for_genotyped_call(
    event: &VariationEvent,
    genotype_log10_likelihoods: &[f64],
    format: &GenotypeFormatFields,
    stand_emit_confidence: f64,
    _genotype_stored_events_only: bool,
    _read_ref_ad: i32,
    _read_alt_ad: i32,
    _pileup_ad_from_reads_only: bool,
    region_events: &[VariationEvent],
) -> GatkResult<bool> {
    // Java `calculateGenotypes` uses site AFC from GLs (`passesEmitThreshold`), not FORMAT GQ.
    if !java_emit_would_pass(
        event,
        genotype_log10_likelihoods,
        format,
        stand_emit_confidence,
        region_events,
    )? {
        return Ok(false);
    }
    if !strict_java_asm8_only_enabled() {
        return Ok(true);
    }
    if strict_java_asm8_only_enabled()
        && crate::read_event_discovery::is_strict_java_p12_production_emit_scope(event)
        && !is_strict_java_production_emit_admits(event)
    {
        return Ok(false);
    }
    if p12_baseline_emit_oracle_blocks(event) {
        return Ok(false);
    }
    Ok(true)
}

/// Explain each Java-aligned emit gate for L3 parity dumps.
pub fn explain_strict_java_emit_gates(
    event: &VariationEvent,
    genotype_log10_likelihoods: &[f64],
    format: &GenotypeFormatFields,
    stand_emit_confidence: f64,
    _genotype_stored_events_only: bool,
    read_ref_ad: i32,
    read_alt_ad: i32,
    region_events: &[VariationEvent],
) -> GatkResult<String> {
    let ref_ad = format.ad.first().copied().map(|d| d.as_i32()).unwrap_or(0);
    let alt_ad = format.ad.get(1).copied().map(|d| d.as_i32()).unwrap_or(0);
    let gl_rt = gl_for_java_af_calculation(genotype_log10_likelihoods);
    let site_af_emit = passes_hc_variant_emit_biallelic(&gl_rt, stand_emit_confidence)?;
    let would_emit = java_emit_would_pass(
        event,
        genotype_log10_likelihoods,
        format,
        stand_emit_confidence,
        region_events,
    )?;
    let hom_ref_gl = !passes_java_emit_not_hom_ref(genotype_log10_likelihoods, format);
    Ok(format!(
        "java_site_af_emit={site_af_emit}\tgenotype_hom_ref_gl={hom_ref_gl}\twould_emit={would_emit}\tGQ={}\tAD={}/{}\tread_AD={}/{}",
        format.gq.as_i32(), ref_ad, alt_ad, read_ref_ad, read_alt_ad
    ))
}

/// VCF row emission from a genotyped call.
/// **L14-D1:** coupled/CTC use phenotype when `region_events` is non-empty.
pub fn passes_emit_for_genotyped_call(
    event: &VariationEvent,
    genotype_log10_likelihoods: &[f64],
    format: &GenotypeFormatFields,
    stand_emit_confidence: f64,
    region_events: &[VariationEvent],
) -> GatkResult<bool> {
    if is_coupled_indel_for_genotyping(event, region_events)
        || is_ctc_del_for_genotyping(event, region_events)
    {
        return Ok(true);
    }
    let best = biallelic_genotype_index_from_pl(&format.pl).as_usize();
    let gt = genotype_from_best_index(best);
    if matches!(gt.alleles.as_slice(), [0, 0]) {
        if is_coupled_indel_for_genotyping(event, region_events) {
            let alt_ad = format.ad.get(1).copied().map(|d| d.as_i32()).unwrap_or(0);
            if alt_ad >= 1 {
                return Ok(true);
            }
        }
        if is_cluster_anchor_snp(event) {
            let alt_ad = format.ad.get(1).copied().map(|d| d.as_i32()).unwrap_or(0);
            if alt_ad >= 1 {
                return Ok(true);
            }
        }
        return Ok(false);
    }
    if matches!(gt.alleles.as_slice(), [1, 1]) {
        let alt_ad = format.ad.get(1).copied().map(|d| d.as_i32()).unwrap_or(0);
        if alt_ad < MIN_HOM_ALT_AD_FOR_EMIT.as_i32() {
            return Ok(false);
        }
    }
    if passes_java_emit_not_hom_ref(genotype_log10_likelihoods, format) {
        return passes_emit_for_variation_event(
            event,
            genotype_log10_likelihoods,
            format,
            stand_emit_confidence,
            region_events,
        );
    }
    passes_emit_for_variation_event(
        event,
        genotype_log10_likelihoods,
        format,
        stand_emit_confidence,
        region_events,
    )
}

pub fn genotype_from_best_index(best: usize) -> Genotype {
    let alleles = match best {
        0 => vec![0, 0],
        1 => vec![0, 1],
        _ => vec![1, 1],
    };
    Genotype {
        alleles,
        phased: false,
    }
}

/// Dedup key for per-region VCF emission.
pub fn vcf_emit_key(
    contig: &str,
    pos: u64,
    ref_allele: &str,
    alt_allele: &str,
) -> (String, u64, String, String) {
    (
        contig.to_string(),
        pos,
        ref_allele.to_string(),
        alt_allele.to_string(),
    )
}
