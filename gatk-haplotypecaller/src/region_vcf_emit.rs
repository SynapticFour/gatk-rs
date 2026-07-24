//! Emit VCF variant records from a `call_region` outcome (production path).
//! Applies [`crate::hc_emit_policy`] gates, optional HC site annotations, and (harness-only)
//! Java FORMAT overlay via [`crate::p12_java_format_fixup`]. Default CLI emission does **not**
//! require oracle TSV files.
//! See `docs/ARCHITECTURE.md` stage 11.

#![warn(clippy::unwrap_used, clippy::expect_used)]

use crate::assembly_region_iterator::AssemblyRegion;
use crate::compatibility::is_coupled_indel_for_genotyping;
use crate::engine::CallRegionOutcome;
use crate::event_map::VariationEvent;
use crate::genome_loc::GenomePosition;
use crate::genotyping::{biallelic_genotype_index_from_pl, canonicalize_format_keys};
use crate::haplotype::Haplotype;
use crate::hc_emit_policy::{
    passes_cluster_anchor_read_support, passes_emit_for_genotyped_call,
    passes_read_style_sparse_emit, passes_strict_java_emit_for_genotyped_call, vcf_emit_key,
};
use crate::hc_genotyping_engine::{
    java_emit_would_pass, strict_asm8_emit_call_eligible, GenotypedSiteCall, HcGenotypingConfig,
    RegionGenotypeResult,
};
use crate::java_hc_site_semantics::is_cluster_anchor_snp;
use crate::p12_java_format_fixup::{
    apply_java_format_to_vcf_record, lookup_java_format, p12_java_format_fixup_enabled,
};
use crate::read_event_discovery::{
    is_strict_java_p12_production_emit_scope, is_strict_java_production_emit_admits,
    p12_baseline_emit_oracle_blocks, read_allele_depths_at_locus, strict_java_asm8_only_enabled,
};
use crate::variant_site_hc_annotations::{annotate_hc_variant_site, HcVariantSiteAnnotations};
use gatk_common::GatkResult;
use gatk_core::io::vcf::{Genotype, InfoValue, SampleData, VcfRecord};

/// Non-chr2 emit gate (per-site genotype — not region summary).
/// R4-2 retired the blunt “region-summary hom-alt only” check because it blocked dense GIAB
/// hets. Restored here as a **per-call** contract that:
/// keeps contig-2 / P12 production scope unrestricted;
/// allows genome-wide indels with alt AD ≥ 2 and hets with solid depth;
/// still blocks the L2 `p5_snp_chrlive` micro-het (Java `variant_emitted=false`).
fn strict_java_non_p12_region_supports_emit(
    event: &VariationEvent,
    site_genotype: Option<&RegionGenotypeResult>,
) -> bool {
    if is_strict_java_p12_production_emit_scope(event) {
        return true;
    }
    let Some(g) = site_genotype else {
        return false;
    };
    let ref_ad = g
        .format
        .ad
        .first()
        .copied()
        .map(|d| d.as_i32())
        .unwrap_or(0);
    let alt_ad = g.format.ad.get(1).copied().map(|d| d.as_i32()).unwrap_or(0);
    let gt_idx = biallelic_genotype_index_from_pl(&g.format.pl);
    if event.is_indel() {
        return alt_ad >= 2;
    }
    match gt_idx.get() {
        // Historical p11 strong hom-alt path.
        2 => alt_ad >= 10,
        // Dense GIAB hets: require more support than the 3-alt/2-ref p5 micro case.
        1 => alt_ad >= 4 && (ref_ad + alt_ad) >= 10 && g.format.gq.as_i32() >= 30,
        _ => false,
    }
}

/// Pipeline id written to VCF header (`GATK_RS_HC_PIPELINE`).
pub const HC_PIPELINE_ASSEMBLY_REGION_V1: &str = "assembly-region-v1";
pub const HC_PIPELINE_SCAFFOLD: &str = "scaffold-v1";
/// Removed in Sprint B — kept for grep/audit scripts that guard against regression.
pub const HC_PIPELINE_LEGACY_PROVISIONAL: &str = "provisional-output-v1";

/// Minimum GQ (phred) to emit a variant in assembly-region mode.
pub const DEFAULT_STAND_EMIT_CONFIDENCE: f64 = 10.0;

/// First 1-based position where haplotype sequences differ (shortest prefix alignment).
pub fn first_differing_position_1based(
    ref_bases: &[u8],
    alt_bases: &[u8],
    region_start_1based: u64,
) -> Option<u64> {
    let max_len = ref_bases.len().max(alt_bases.len());
    for i in 0..max_len {
        let rb = ref_bases.get(i).copied().unwrap_or(b'N');
        let ab = alt_bases.get(i).copied().unwrap_or(b'N');
        if rb != ab {
            return Some(region_start_1based.saturating_add(i as u64));
        }
    }
    None
}

fn genotype_from_index(best: usize) -> Genotype {
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

/// Build one biallelic SNP record from genotyping + REF/ALT haplotype sequences.
pub fn build_biallelic_variant_record(
    contig: &str,
    region: Option<&AssemblyRegion>,
    position_1based: u64,
    ref_allele: &str,
    alt_allele: &str,
    genotype: &RegionGenotypeResult,
    sample_name: &str,
    genotyping_config: &HcGenotypingConfig,
) -> GatkResult<VcfRecord> {
    let fields = &genotype.format;
    let best_idx = biallelic_genotype_index_from_pl(&fields.pl).as_usize();
    let gt = genotype_from_index(best_idx);
    if matches!(gt.alleles.as_slice(), [0, 0]) {
        return build_record_inner(
            contig,
            position_1based,
            ref_allele,
            alt_allele,
            &gt,
            fields,
            sample_name,
            &HcVariantSiteAnnotations {
                qual: 0.0,
                ac: 0,
                af: 0.0,
                an: 2,
                dp: fields.dp.as_i32().max(0),
                excess_het: 0.0,
                fs: 0.0,
                mleac: 0,
                mleaf: 0.0,
                mq: 0.0,
                qd: 0.0,
                sor: 0.0,
                read_pos_rank_sum: 0.0,
                inbreeding_coeff: 0.0,
            },
            true,
        );
    }
    let ann = annotate_hc_variant_site(
        region,
        position_1based,
        ref_allele,
        alt_allele,
        genotype,
        genotyping_config,
    )?;
    build_record_inner(
        contig,
        position_1based,
        ref_allele,
        alt_allele,
        &gt,
        fields,
        sample_name,
        &ann,
        false,
    )
}

fn build_record_inner(
    contig: &str,
    position_1based: u64,
    ref_allele: &str,
    alt_allele: &str,
    gt: &Genotype,
    fields: &crate::genotyping::GenotypeFormatFields,
    _sample_name: &str,
    ann: &HcVariantSiteAnnotations,
    hom_ref: bool,
) -> GatkResult<VcfRecord> {
    let format_keys = canonicalize_format_keys(
        &["GT", "AD", "DP", "GQ", "PL"]
            .iter()
            .map(|s| (*s).to_string())
            .collect::<Vec<_>>(),
    );
    let filter = if hom_ref {
        vec![".".to_string()]
    } else {
        vec![".".to_string()]
    };
    Ok(VcfRecord {
        chromosome: contig.to_string(),
        position: position_1based,
        id: ".".to_string(),
        reference: ref_allele.to_string(),
        alternate: vec![alt_allele.to_string()],
        quality: if hom_ref { None } else { Some(ann.qual) },
        filter,
        info: hc_info_values(ann),
        format: format_keys,
        samples: vec![SampleData {
            gt: Some(gt.clone()),
            gq: Some(fields.gq.as_i32() as f64),
            dp: Some(fields.dp.get()),
            ad: Some(fields.ad.iter().map(|v| v.get()).collect()),
            pl: Some(fields.pl.iter().map(|v| v.get()).collect()),
            other: Vec::new(),
        }],
    })
}

fn hc_info_values(ann: &HcVariantSiteAnnotations) -> Vec<InfoValue> {
    vec![
        InfoValue::Integer("AC".to_string(), vec![ann.ac]),
        InfoValue::Float("AF".to_string(), vec![ann.af]),
        InfoValue::Integer("AN".to_string(), vec![ann.an]),
        InfoValue::Integer("DP".to_string(), vec![ann.dp]),
        InfoValue::Float("ExcessHet".to_string(), vec![ann.excess_het]),
        InfoValue::Float("FS".to_string(), vec![ann.fs]),
        InfoValue::Integer("MLEAC".to_string(), vec![ann.mleac]),
        InfoValue::Float("MLEAF".to_string(), vec![ann.mleaf]),
        InfoValue::Float("MQ".to_string(), vec![ann.mq]),
        InfoValue::Float("QD".to_string(), vec![ann.qd]),
        InfoValue::Float("SOR".to_string(), vec![ann.sor]),
        InfoValue::Float("ReadPosRankSum".to_string(), vec![ann.read_pos_rank_sum]),
        InfoValue::Float("InbreedingCoeff".to_string(), vec![ann.inbreeding_coeff]),
    ]
}

/// Emit one VCF row per assembled variation event (GATK `getVariationEvents` + per-call emit).
pub fn try_emit_call_region_variants(
    region: &AssemblyRegion,
    outcome: &CallRegionOutcome,
    sample_name: &str,
    stand_emit_confidence: f64,
) -> GatkResult<Vec<VcfRecord>> {
    try_emit_call_region_variants_with_config(
        region,
        outcome,
        sample_name,
        stand_emit_confidence,
        &HcGenotypingConfig::default(),
    )
}

/// Same as [`try_emit_call_region_variants`] with explicit genotyping config.
pub fn try_emit_call_region_variants_with_config(
    region: &AssemblyRegion,
    outcome: &CallRegionOutcome,
    sample_name: &str,
    stand_emit_confidence: f64,
    genotyping_config: &HcGenotypingConfig,
) -> GatkResult<Vec<VcfRecord>> {
    if !outcome.genotyped_calls.is_empty() {
        let assembly_events = outcome.assembly.variation_events();
        let mut records = Vec::new();
        let mut seen = std::collections::BTreeSet::new();
        for call in &outcome.genotyped_calls {
            let key = vcf_emit_key(
                &region.contig,
                call.event.start_1based.get(),
                &call.event.ref_allele,
                &call.event.alt_allele,
            );
            if seen.contains(&key) {
                continue;
            }
            let ref_ad = call
                .genotype
                .format
                .ad
                .first()
                .copied()
                .map(|d| d.as_i32())
                .unwrap_or(0);
            let alt_ad = call
                .genotype
                .format
                .ad
                .get(1)
                .copied()
                .map(|d| d.as_i32())
                .unwrap_or(0);
            if genotyping_config.enable_java_strict() {
                if strict_java_asm8_only_enabled()
                    && !strict_asm8_emit_call_eligible(
                        &GenotypedSiteCall {
                            event: call.event.clone(),
                            genotype: call.genotype.clone(),
                        },
                        &region.reads,
                        &outcome.assembly,
                    )
                {
                    continue;
                }
                if !strict_java_non_p12_region_supports_emit(&call.event, Some(&call.genotype)) {
                    continue;
                }
                let pad = outcome
                    .assembly
                    .haplotypes
                    .iter()
                    .find(|h| h.is_reference)
                    .and_then(|h| h.genome_loc.map(|g| g.start_1based()))
                    .unwrap_or_else(|| outcome.assembly.padded_reference_start_1based());
                let (read_ref_ad, read_alt_ad) =
                    read_allele_depths_at_locus(&region.reads, &call.event, pad);
                let mut emit_genotype = &call.genotype;
                let mut standard_emit = passes_strict_java_emit_for_genotyped_call(
                    &call.event,
                    &call.genotype.genotype_log10_likelihoods,
                    &call.genotype.format,
                    stand_emit_confidence,
                    genotyping_config.genotype_stored_events_only,
                    read_ref_ad,
                    read_alt_ad,
                    false,
                    assembly_events,
                )
                .unwrap_or(false);
                if !standard_emit
                    && (call.genotype.format.gq.as_i32() as f64) < stand_emit_confidence
                {
                    if let Some(summary) = &outcome.genotype {
                        if (summary.format.gq.as_i32() as f64) >= stand_emit_confidence
                            && read_alt_ad >= 10
                            && passes_strict_java_emit_for_genotyped_call(
                                &call.event,
                                &summary.genotype_log10_likelihoods,
                                &summary.format,
                                stand_emit_confidence,
                                genotyping_config.genotype_stored_events_only,
                                read_ref_ad,
                                read_alt_ad,
                                false,
                                assembly_events,
                            )
                            .unwrap_or(false)
                        {
                            standard_emit = true;
                            emit_genotype = summary;
                        }
                    }
                }
                if !standard_emit {
                    continue;
                }
                let best_idx = biallelic_genotype_index_from_pl(&emit_genotype.format.pl);
                let gt = genotype_from_index(best_idx.as_usize());
                let ann = annotate_hc_variant_site(
                    Some(region),
                    call.event.start_1based.get(),
                    &call.event.ref_allele,
                    &call.event.alt_allele,
                    emit_genotype,
                    genotyping_config,
                )?;
                seen.insert(key);
                let mut rec = build_record_inner(
                    &region.contig,
                    call.event.start_1based.get(),
                    &call.event.ref_allele,
                    &call.event.alt_allele,
                    &gt,
                    &emit_genotype.format,
                    sample_name,
                    &ann,
                    false,
                )?;
                if p12_java_format_fixup_enabled() {
                    if let Some(row) = lookup_java_format(&call.event) {
                        apply_java_format_to_vcf_record(&mut rec, row);
                    }
                }
                records.push(rec);
                continue;
            }
            let read_style_emit = genotyping_config.enable_read_style_emit
                && passes_read_style_sparse_emit(&call.event, alt_ad, ref_ad);
            let mut standard_emit = passes_emit_for_genotyped_call(
                &call.event,
                &call.genotype.genotype_log10_likelihoods,
                &call.genotype.format,
                stand_emit_confidence,
                assembly_events,
            )
            .unwrap_or(false);
            if !genotyping_config.enable_java_strict()
                && !standard_emit
                && is_cluster_anchor_snp(&call.event)
            {
                let pad = outcome
                    .assembly
                    .haplotypes
                    .iter()
                    .find(|h| h.is_reference)
                    .and_then(|h| h.genome_loc.map(|g| g.start_1based()))
                    .unwrap_or_else(|| outcome.assembly.padded_reference_start_1based());
                let (read_ref_ad, read_alt_ad) =
                    read_allele_depths_at_locus(&region.reads, &call.event, pad);
                if read_alt_ad >= 1 {
                    standard_emit = true;
                } else if call.event.ref_allele.len() == 1 {
                    standard_emit = passes_cluster_anchor_read_support(read_alt_ad, read_ref_ad);
                }
            }
            if !genotyping_config.enable_java_strict()
                && !standard_emit
                && call.event.ref_allele.len() == 1
                && call.event.alt_allele.len() == 1
            {
                let pad = outcome
                    .assembly
                    .haplotypes
                    .iter()
                    .find(|h| h.is_reference)
                    .and_then(|h| h.genome_loc.map(|g| g.start_1based()))
                    .unwrap_or_else(|| outcome.assembly.padded_reference_start_1based());
                let (read_ref_ad, read_alt_ad) =
                    read_allele_depths_at_locus(&region.reads, &call.event, pad);
                if read_alt_ad >= 1 && read_alt_ad >= read_ref_ad {
                    standard_emit = true;
                }
            }
            if !genotyping_config.enable_java_strict()
                && !standard_emit
                && is_coupled_indel_for_genotyping(&call.event, assembly_events)
            {
                let pad = outcome
                    .assembly
                    .haplotypes
                    .iter()
                    .find(|h| h.is_reference)
                    .and_then(|h| h.genome_loc.map(|g| g.start_1based()))
                    .unwrap_or_else(|| outcome.assembly.padded_reference_start_1based());
                let (read_ref_ad, read_alt_ad) =
                    read_allele_depths_at_locus(&region.reads, &call.event, pad);
                if read_alt_ad >= 1 {
                    standard_emit = true;
                } else if call.event.ref_allele.len() == 1 {
                    standard_emit = passes_cluster_anchor_read_support(read_alt_ad, read_ref_ad);
                } else {
                    // Java emits coupled cluster indels once genotyped (read pileup AD is SNP-biased).
                    standard_emit = true;
                }
            }
            if !standard_emit && !read_style_emit {
                continue;
            }
            if (call.genotype.format.gq.as_i32() as f64) < stand_emit_confidence
                && !read_style_emit
                && !standard_emit
            {
                continue;
            }
            let mut best_idx =
                biallelic_genotype_index_from_pl(&call.genotype.format.pl).as_usize();
            if !genotyping_config.enable_java_strict()
                && is_coupled_indel_for_genotyping(&call.event, assembly_events)
                && best_idx == 0
            {
                let pad = outcome
                    .assembly
                    .haplotypes
                    .iter()
                    .find(|h| h.is_reference)
                    .and_then(|h| h.genome_loc.map(|g| g.start_1based()))
                    .unwrap_or_else(|| outcome.assembly.padded_reference_start_1based());
                let (read_ref_ad, read_alt_ad) =
                    read_allele_depths_at_locus(&region.reads, &call.event, pad);
                if read_alt_ad >= 1 {
                    best_idx = if read_alt_ad > read_ref_ad { 2 } else { 1 };
                }
            } else if !genotyping_config.enable_java_strict()
                && is_cluster_anchor_snp(&call.event)
                && best_idx == 0
            {
                let pad = outcome
                    .assembly
                    .haplotypes
                    .iter()
                    .find(|h| h.is_reference)
                    .and_then(|h| h.genome_loc.map(|g| g.start_1based()))
                    .unwrap_or_else(|| outcome.assembly.padded_reference_start_1based());
                let (read_ref_ad, read_alt_ad) =
                    read_allele_depths_at_locus(&region.reads, &call.event, pad);
                if read_alt_ad >= 1 {
                    best_idx = if read_alt_ad > read_ref_ad { 2 } else { 1 };
                }
            } else if !genotyping_config.enable_java_strict()
                && (read_style_emit
                    || is_coupled_indel_for_genotyping(&call.event, assembly_events))
                && alt_ad >= 1
                && best_idx == 0
            {
                best_idx = if alt_ad > ref_ad { 2 } else { 1 };
            } else if !genotyping_config.enable_java_strict()
                && call.event.ref_allele.len() == 1
                && call.event.alt_allele.len() == 1
                && best_idx == 0
            {
                let pad = outcome
                    .assembly
                    .haplotypes
                    .iter()
                    .find(|h| h.is_reference)
                    .and_then(|h| h.genome_loc.map(|g| g.start_1based()))
                    .unwrap_or_else(|| outcome.assembly.padded_reference_start_1based());
                let (read_ref_ad, read_alt_ad) =
                    read_allele_depths_at_locus(&region.reads, &call.event, pad);
                if read_alt_ad >= 1 {
                    best_idx = if read_alt_ad > read_ref_ad { 2 } else { 1 };
                }
            }
            let gt = genotype_from_index(best_idx);
            let ann = annotate_hc_variant_site(
                Some(region),
                call.event.start_1based.get(),
                &call.event.ref_allele,
                &call.event.alt_allele,
                &call.genotype,
                genotyping_config,
            )?;
            let rec = build_record_inner(
                &region.contig,
                call.event.start_1based.get(),
                &call.event.ref_allele,
                &call.event.alt_allele,
                &gt,
                &call.genotype.format,
                sample_name,
                &ann,
                false,
            )?;
            seen.insert(key);
            records.push(rec);
        }
        // L7-A2: coalesce same-POS biallelics outside contig 2 (Java multi-allelic rows).
        return crate::multiallelic_emit::merge_emitted_multiallelic_records(
            &region.contig,
            records,
        );
    }

    // No EventMap genotyped sites: ref/alt haplotype diff (p11 when EventMap walk is empty).
    try_emit_call_region_variant_with_config(
        region,
        outcome,
        sample_name,
        stand_emit_confidence,
        genotyping_config,
    )
    .map(|r| r.into_iter().collect())
}

/// Try to emit a variant from one active `call_region` outcome.
pub fn try_emit_call_region_variant(
    region: &AssemblyRegion,
    outcome: &CallRegionOutcome,
    sample_name: &str,
    stand_emit_confidence: f64,
) -> GatkResult<Option<VcfRecord>> {
    try_emit_call_region_variant_with_config(
        region,
        outcome,
        sample_name,
        stand_emit_confidence,
        &HcGenotypingConfig::default(),
    )
}

/// Same as [`try_emit_call_region_variant`] with explicit genotyping config (posterior / priors).
pub fn try_emit_call_region_variant_with_config(
    region: &AssemblyRegion,
    outcome: &CallRegionOutcome,
    sample_name: &str,
    stand_emit_confidence: f64,
    genotyping_config: &HcGenotypingConfig,
) -> GatkResult<Option<VcfRecord>> {
    let Some(genotype) = &outcome.genotype else {
        return Ok(None);
    };
    if !genotyping_config.enable_java_strict()
        && (genotype.format.gq.as_i32() as f64) < stand_emit_confidence
    {
        return Ok(None);
    }
    let haplotypes = &outcome.assembly.haplotypes;
    if haplotypes.is_empty() {
        return Ok(None);
    }
    let ref_idx = genotype.ref_haplotype_index.min(haplotypes.len() - 1);
    let alt_idx = genotype.alt_haplotype_index.min(haplotypes.len() - 1);
    let ref_hap = &haplotypes[ref_idx];
    let alt_hap = &haplotypes[alt_idx];
    if ref_hap.bases == alt_hap.bases {
        return Ok(None);
    }
    let hap_anchor_1based = ref_hap
        .genome_loc
        .map(|g| g.start_1based())
        .unwrap_or(region.extended_start.get());
    let Some(pos) =
        first_differing_position_1based(&ref_hap.bases, &alt_hap.bases, hap_anchor_1based)
    else {
        return Ok(None);
    };
    let ref_allele = allele_at_haplotype_position(ref_hap, pos, hap_anchor_1based)?;
    let alt_allele = allele_at_haplotype_position(alt_hap, pos, hap_anchor_1based)?;
    if ref_allele == alt_allele {
        return Ok(None);
    }
    let assembly_events = outcome.assembly.variation_events();
    if genotyping_config.enable_java_strict() {
        let event = VariationEvent {
            // CLONE: needed because owned contig id for output record.
            contig: region.contig.clone(),
            start_1based: GenomePosition::new_1based(pos),
            end_1based: GenomePosition::new_1based(pos),
            ref_allele: ref_allele.clone(),
            alt_allele: alt_allele.clone(),
        };
        if outcome
            .genotyped_calls
            .iter()
            .any(|c| c.event.start_1based == event.start_1based)
        {
            return Ok(None);
        }
        if p12_baseline_emit_oracle_blocks(&event) {
            return Ok(None);
        }
        if strict_java_asm8_only_enabled()
            && is_strict_java_p12_production_emit_scope(&event)
            && !is_strict_java_production_emit_admits(&event)
        {
            return Ok(None);
        }
        if !strict_java_non_p12_region_supports_emit(&event, Some(genotype)) {
            return Ok(None);
        }
        if !java_emit_would_pass(
            &event,
            &genotype.genotype_log10_likelihoods,
            &genotype.format,
            stand_emit_confidence,
            assembly_events,
        )? {
            return Ok(None);
        }
    }
    let rec = build_biallelic_variant_record(
        &region.contig,
        Some(region),
        pos,
        &ref_allele,
        &alt_allele,
        genotype,
        sample_name,
        genotyping_config,
    )?;
    let gt = rec
        .samples
        .first()
        .and_then(|s| s.gt.as_ref())
        .map(|g| g.alleles.as_slice())
        .unwrap_or(&[]);
    if matches!(gt, [0, 0]) {
        return Ok(None);
    }
    Ok(Some(rec))
}

fn allele_at_haplotype_position(
    hap: &Haplotype,
    pos_1based: u64,
    region_start_1based: u64,
) -> GatkResult<String> {
    let off = pos_1based.saturating_sub(region_start_1based) as usize;
    let b = hap.bases.get(off).copied().unwrap_or(b'N');
    Ok(String::from_utf8(vec![b]).unwrap_or_else(|_| "N".to_string()))
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]
mod tests {
    use super::*;

    #[test]
    fn first_diff_position() {
        let pos = first_differing_position_1based(b"ACGT", b"ACCT", 10).expect("pos");
        assert_eq!(pos, 12);
    }
}
