//! GATK HC variant-site INFO / QUAL for assembly-region VCF emission (J-D01).

use crate::activity_scoring::genotype_log10_likelihoods_after_java_genotype_pl_roundtrip;
use crate::af_calc::{
    biallelic_alt_pseudocount, calculate_biallelic_af_em,
    calculate_biallelic_af_em_with_alt_pseudocount, diploid_af_log10_prob_only_ref_allele_exists,
    AfCalculatorConfig,
};
use crate::annotator::plugins::{
    excess_het, fisher_strand, mapping_quality_rank_sum, qual_by_depth, rank_sum_baseq,
    read_pos_rank_sum, strand_odds_ratio,
};
use crate::assembly_region_iterator::AssemblyRegion;
use crate::event_map::{
    build_per_haplotype_variation_events, merged_alleles_for_genotyping, overlapping_events,
    variation_events_at_position_from_cache, VariationEvent,
};
use crate::fragment_overlap::{
    pre_realign_alignment, read_base_at_ref_coord_1based, read_base_quality_at_ref_coord_1based,
};
use crate::genotyping::GenotypeFormatFields;
use crate::haplotype::Haplotype;
use crate::hc_allele_mapping::{create_allele_mapper_with_events, replace_span_del_events};
use crate::hc_genotyping_engine::{
    java_create_allele_mapper_pools, marginalize_rows_to_biallelic_alleles,
    region_likelihoods_to_rows, HcGenotypingConfig, RegionGenotypeResult,
};
use crate::read_model::MAPPING_QUALITY_UNAVAILABLE;
use crate::read_realignment::LOG_10_INFORMATIVE_THRESHOLD;
use crate::region_read_likelihood::RegionReadLikelihood;
use crate::shared_bam::SharedBamRecord;
use gatk_common::GatkResult;
use gatk_core::io::vcf::Genotype;
use rust_htslib::bam::record::{Cigar, CigarString};
use rust_htslib::bam::Record;
use std::collections::BTreeSet;

const FLAG_REVERSE: u16 = 0x10;
/// GATK `StrandBiasTest.ARRAY_DIM` — `FisherStrand.MIN_COUNT`.
pub const FISHER_STRAND_MIN_COUNT: u32 = 2;
/// GATK `StrandOddsRatio.MIN_COUNT`.
pub const STRAND_ODDS_RATIO_MIN_COUNT: u32 = 0;

/// HC INFO + QUAL slice aligned with GATK `HaplotypeCaller` default annotators on biallelic sites.
/// # Invariants
/// Fields mirror GATK default INFO keys for emitted biallelic variant sites (AC/AF/AN/DP/FS/SOR/etc.).
/// `qual` uses AF-calculator path with Java-observed posterior convention for emitted variants.
/// # Ownership
/// Owns scalar annotation values; genotype and region data are consumed at construction time.
/// # Mutation
/// Immutable annotation bundle produced by [`annotate_hc_variant_site`].
/// # Biological assumptions
/// Single-sample diploid biallelic SNP/indels with standard short-read INFO semantics.
/// # Java equivalence
/// GATK `HaplotypeCaller` default variant annotators (J-D01 INFO/QUAL slice).
#[derive(Debug, Clone, PartialEq)]
pub struct HcVariantSiteAnnotations {
    pub qual: f64,
    pub ac: i32,
    pub af: f64,
    pub an: i32,
    pub dp: i32,
    pub excess_het: f64,
    pub fs: f64,
    pub mleac: i32,
    pub mleaf: f64,
    pub mq: f64,
    pub qd: f64,
    pub sor: f64,
    /// `None` = Java `RankSumTest` emptyMap (undefined / NaN). `Some(0.0)` emits.
    pub read_pos_rank_sum: Option<f64>,
    /// 6R.193: `BaseQRankSum`. Same Option contract as [`Self::read_pos_rank_sum`].
    pub base_q_rank_sum: Option<f64>,
    /// 6R.194: `MQRankSum`. Same Option contract as [`Self::read_pos_rank_sum`].
    pub mq_rank_sum: Option<f64>,
    pub inbreeding_coeff: f64,
}

/// Phred-scaled site QUAL from GATK `GenotypingEngine` + `AlleleFrequencyCalculator`
/// (`use-posteriors-to-calculate-qual` false, default HC).
pub fn qual_from_af_calculation(genotype_log10_likelihoods: &[f64]) -> GatkResult<f64> {
    let af = calculate_biallelic_af_em(
        &[genotype_log10_likelihoods],
        &AfCalculatorConfig::default(),
    )?;
    // Observed on p11: Java uses `log10ProbOnlyRefAlleleExists` for emitted variant QUAL.
    let log10_conf = af.log10_posterior_no_variant + 0.0;
    Ok((-10.0 * log10_conf) + 0.0)
}

/// Site QUAL from GATK `AlleleFrequencyCalculator.calculate` on a diploid merged VC.
/// `alleles` is REF first. When `*` is present, P(no variant) includes REF/SPAN_DEL genotypes.
pub fn qual_from_merged_diploid_af_calculate(
    log10_likelihoods: &[f64],
    alleles: &[&str],
) -> GatkResult<f64> {
    let log10_pe = diploid_af_log10_prob_only_ref_allele_exists(
        log10_likelihoods,
        alleles,
        &AfCalculatorConfig::default(),
    )?;
    Ok((-10.0 * log10_pe.min(0.0)) + 0.0)
}

/// GATK `QualByDepth.getDepth` for a single variant genotype.
pub fn qd_depth_for_variant(gt: &Genotype, fields: &GenotypeFormatFields) -> i32 {
    if !is_het_or_hom_var(gt) {
        return 0;
    }
    if !fields.ad.is_empty() {
        let total_ad: i32 = fields.ad.iter().map(|d| d.as_i32()).sum();
        if total_ad != 0 {
            let alt_ad = total_ad.saturating_sub(fields.ad[0].as_i32().max(0));
            let ad_restricted = if alt_ad > 1 { total_ad } else { 0 };
            if ad_restricted > 0 {
                return ad_restricted;
            }
            return total_ad;
        }
    }
    fields.dp.as_i32().max(0)
}

fn is_het_or_hom_var(gt: &Genotype) -> bool {
    matches!(gt.alleles.as_slice(), [0, 1] | [1, 0] | [1, 1])
}

/// Post-filter `AlleleLikelihoods` slice for Java `StrandBiasTest.getContingencyTable`.
///
/// # Java equivalence
/// GATK 4.4 `FisherStrand` / `StrandOddsRatio` `calculateAnnotationFromLikelihoods`
/// consume the same allele-marginalized evidence object passed to
/// `VariantAnnotatorEngine.annotateContext` (after `filterPoorlyModeledEvidence`).
/// MQ uses [`Self::reads`] as Java `likelihoods.sampleEvidence` (6R.167);
/// it does **not** apply the informative best-allele filter used by FS/SOR.
pub struct HcStrandBiasLikelihoods<'a> {
    pub reads: &'a [SharedBamRecord],
    pub likelihoods: &'a [RegionReadLikelihood],
    pub haplotypes: &'a [Haplotype],
    pub contig: &'a str,
    pub ref_bytes: &'a [u8],
    pub pad_start_1based: u64,
    pub full_ref_bytes: &'a [u8],
    pub full_pad_1based: u64,
    pub max_mnp_distance: usize,
    pub emit_spanning_dels: bool,
}

fn apply_strand_min_count(table: (u32, u32, u32, u32), min_count: u32) -> (u32, u32, u32, u32) {
    let total = table.0 + table.1 + table.2 + table.3;
    if total > min_count {
        table
    } else {
        (0, 0, 0, 0)
    }
}

/// GATK `StrandBiasTest.getContingencyTable` (REF/ALT × forward/reverse).
///
/// Per-read: best allele over the locus candidate set (`searchBestAllele`,
/// reference priority 1.0), then `isInformative` (`confidence > 0.2`).
/// Only an informative best allele equal to the emitted REF or ALT is counted.
/// Sample table is copied only when total counts `> min_count`.
pub fn strand_bias_contingency_table(
    evidence: &HcStrandBiasLikelihoods<'_>,
    loc_1based: u64,
    ref_allele: &str,
    alt_allele: &str,
    min_count: u32,
) -> (u32, u32, u32, u32) {
    apply_strand_min_count(
        strand_bias_sample_counts(evidence, loc_1based, ref_allele, alt_allele, None),
        min_count,
    )
}

fn strand_bias_sample_counts(
    evidence: &HcStrandBiasLikelihoods<'_>,
    loc_1based: u64,
    ref_allele: &str,
    alt_allele: &str,
    annotation_output: PretrimAnnotationAlleles<'_>,
) -> (u32, u32, u32, u32) {
    let mut ref_fw = 0u32;
    let mut ref_rv = 0u32;
    let mut alt_fw = 0u32;
    let mut alt_rv = 0u32;
    for decision in annotation_allele_decisions(
        evidence,
        loc_1based,
        ref_allele,
        alt_allele,
        annotation_output,
    ) {
        if !decision.informative {
            continue;
        }
        let Some(rec) = evidence.reads.get(decision.read_index) else {
            continue;
        };
        let reverse = rec.flags() & FLAG_REVERSE != 0;
        match decision.class {
            AnnotationAlleleClass::EmittedRef => {
                if reverse {
                    ref_rv += 1;
                } else {
                    ref_fw += 1;
                }
            }
            AnnotationAlleleClass::EmittedAlt => {
                if reverse {
                    alt_rv += 1;
                } else {
                    alt_fw += 1;
                }
            }
            // A best allele that is not the emitted REF or ALT is ignored, not reassigned.
            AnnotationAlleleClass::Other => {}
        }
    }
    (ref_fw, ref_rv, alt_fw, alt_rv)
}

/// GATK `RankSumTest.isUsableRead`: MQ != 0 and MQ != 255.
fn java_ranksum_mq_usable(rec: &Record) -> bool {
    let mq = rec.mapq();
    mq != 0 && mq != MAPPING_QUALITY_UNAVAILABLE
}

/// GATK `ReadPosRankSumTest.getReadPosition`: `min(left, right)` including hard clips.
fn java_read_pos_rank_element_on_cigar(
    alignment_pos0: i64,
    cig: &CigarString,
    seq_len: i64,
    vc_start_1based: i32,
) -> Option<f64> {
    let mut leading_hard = 0i64;
    let mut trailing_hard = 0i64;
    if let Some(Cigar::HardClip(n)) = cig.iter().next() {
        leading_hard = i64::from(*n);
    }
    if let Some(Cigar::HardClip(n)) = cig.iter().last() {
        trailing_hard = i64::from(*n);
    }
    let vc = i64::from(vc_start_1based);
    let mut ref_pos = alignment_pos0 + 1;
    let mut q = 0usize;
    for c in cig.iter() {
        match c {
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                let n = i64::from(*n);
                if vc >= ref_pos && vc < ref_pos + n {
                    let idx = q + (vc - ref_pos) as usize;
                    let left = leading_hard + idx as i64;
                    let right = seq_len - 1 - idx as i64 + trailing_hard;
                    return Some(left.min(right) as f64);
                }
                ref_pos += n;
                q += n as usize;
            }
            Cigar::Del(n) | Cigar::RefSkip(n) => {
                ref_pos += i64::from(*n);
            }
            Cigar::Ins(n) | Cigar::SoftClip(n) => {
                q += *n as usize;
            }
            Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }
    None
}

/// Inclusive end of the emitted reference allele.
///
/// GATK `VariantContext.getEnd()` is `start + reference.length() - 1`.
fn java_variant_context_end_1based(vc_start_1based: u64, ref_allele: &str) -> i32 {
    let end = vc_start_1based.saturating_add(ref_allele.len().saturating_sub(1) as u64);
    end as i32
}

/// `read.getStart() == vc.getEnd() + 1`.
fn read_starts_on_base_after_variant_end(rec: &Record, vc_end_1based: i32) -> bool {
    rec.pos() + 1 == i64::from(vc_end_1based) + 1
}

/// First attached operator that is not soft or hard clipping is `I`.
///
/// GATK `CigarOperator.isClipping` is `S` or `H` only.
fn first_non_clipping_cigar_operator_is_insertion(cigar: &CigarString) -> bool {
    cigar
        .iter()
        .find(|op| !matches!(op, Cigar::SoftClip(_) | Cigar::HardClip(_)))
        .is_some_and(|op| matches!(op, Cigar::Ins(_)))
}

/// GATK `ReadPosRankSumTest.getReadPosition`.
///
/// A leading insertion whose read start is the reference base immediately
/// after `vc.getEnd()` returns 0 before the coordinate walk. Otherwise
/// `min(left, right)` on the attached CIGAR, then the pre-realign POS/CIGAR
/// when the attached CIGAR does not cover `vc.getStart()` (6R.203).
fn java_read_pos_rank_element(
    rec: &Record,
    vc_start_1based: i32,
    vc_end_1based: i32,
) -> Option<f64> {
    let attached = CigarString::from(rec.cigar().iter().copied().collect::<Vec<_>>());
    if read_starts_on_base_after_variant_end(rec, vc_end_1based)
        && first_non_clipping_cigar_operator_is_insertion(&attached)
    {
        return Some(0.0);
    }
    if let Some(v) = java_read_pos_rank_element_on_cigar(
        rec.pos(),
        &attached,
        rec.seq_len() as i64,
        vc_start_1based,
    ) {
        return Some(v);
    }
    let (pos0, cig) = pre_realign_alignment(rec)?;
    java_read_pos_rank_element_on_cigar(pos0, &cig, rec.seq_len() as i64, vc_start_1based)
}

/// GATK `RankSumTest.fillQualsFromLikelihood` REF/ALT lists.
///
/// Observable Java 4.4 contract on the annotation `AlleleLikelihoods`:
/// remarginalize `evidence.likelihoods` to site REF/ALT,
/// `bestAllelesBreakingTies` (REF wins ties), `isInformative` (gap > 0.2),
/// `isUsableRead`, then `getElementForRead`. Does **not** walk `region.reads`.
///
/// 6R.191: ReadPosRankSum uses this membership plus `getReadPosition`.
/// 6R.193: BaseQRankSum uses the same membership plus
/// `ReadUtils.getReadBaseQualityAtReferenceCoordinate`.
/// 6R.194: MQRankSum uses the same membership plus `read.getMappingQuality()`.
/// Mann-Whitney / undefined-vs-zero stay in the RankSum plugins.
fn fill_ranksum_quals_from_likelihood<F>(
    evidence: &HcStrandBiasLikelihoods<'_>,
    loc_1based: u64,
    ref_allele: &str,
    alt_allele: &str,
    annotation_output: PretrimAnnotationAlleles<'_>,
    mut element_for_read: F,
) -> (Vec<f64>, Vec<f64>)
where
    F: FnMut(&Record, i32) -> Option<f64>,
{
    let mut ref_vals = Vec::new();
    let mut alt_vals = Vec::new();
    for decision in annotation_allele_decisions(
        evidence,
        loc_1based,
        ref_allele,
        alt_allele,
        annotation_output,
    ) {
        if !decision.informative {
            continue;
        }
        if !matches!(
            decision.class,
            AnnotationAlleleClass::EmittedRef | AnnotationAlleleClass::EmittedAlt
        ) {
            continue;
        }
        let Some(rec) = evidence.reads.get(decision.read_index) else {
            continue;
        };
        if !java_ranksum_mq_usable(rec) {
            continue;
        }
        let Some(v) = element_for_read(rec, loc_1based as i32) else {
            continue;
        };
        match decision.class {
            AnnotationAlleleClass::EmittedRef => ref_vals.push(v),
            AnnotationAlleleClass::EmittedAlt => alt_vals.push(v),
            AnnotationAlleleClass::Other => {}
        }
    }
    (ref_vals, alt_vals)
}

/// Best allele for one read after annotation-boundary selection.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum AnnotationAlleleClass {
    EmittedRef,
    EmittedAlt,
    /// Best allele is a candidate that this record does not emit.
    Other,
}

struct AnnotationAlleleDecision {
    read_index: usize,
    /// Index into `[reference, *alts]` for the selected best allele.
    best_index: usize,
    informative: bool,
    class: AnnotationAlleleClass,
}

fn annotation_allele_class(
    bases: &str,
    ref_allele: &str,
    alt_allele: &str,
) -> AnnotationAlleleClass {
    if bases == ref_allele {
        AnnotationAlleleClass::EmittedRef
    } else if bases == alt_allele {
        AnnotationAlleleClass::EmittedAlt
    } else {
        AnnotationAlleleClass::Other
    }
}

/// Select the best allele over every candidate at `loc`, then label it against the
/// emitted REF/ALT pair.
///
/// Java `AlleleLikelihoods.searchBestAllele` walks the merged allele list in order.
/// A strictly greater likelihood replaces the best allele and demotes the previous
/// best to second. An equal likelihood does not replace an earlier allele. When
/// `best - second < 0.2`, alleles within `0.2` of that best likelihood are walked
/// again. Reference priority is `1.0` and every other priority is `0`; only a
/// strictly higher priority replaces, and an allele more than `0.2` below the
/// best likelihood cannot take the tie. Confidence is `0` when the two
/// likelihoods are equal, otherwise `best - second` after that reassignment.
/// The read is informative only when confidence is strictly greater than
/// [`LOG_10_INFORMATIVE_THRESHOLD`] (0.2).
///
/// A site whose events do not recover a multi-allele list containing both emitted
/// alleles keeps the previous two-allele projection.
/// Pre-trim reference and kept output alts. `None` keeps the emitted-allele path.
type PretrimAnnotationAlleles<'a> = Option<(&'a str, &'a [String])>;

fn annotation_allele_decisions(
    evidence: &HcStrandBiasLikelihoods<'_>,
    loc_1based: u64,
    ref_allele: &str,
    alt_allele: &str,
    annotation_output: PretrimAnnotationAlleles<'_>,
) -> Vec<AnnotationAlleleDecision> {
    if evidence.likelihoods.is_empty() || evidence.haplotypes.is_empty() {
        return Vec::new();
    }
    if let Some((pre_ref, kept_alts)) =
        pretrim_output_when_different(ref_allele, alt_allele, annotation_output)
    {
        if let Some(decisions) =
            pretrim_annotation_decisions(evidence, loc_1based, pre_ref, kept_alts)
        {
            return decisions;
        }
    }
    if let Some(decisions) =
        multi_allele_annotation_decisions(evidence, loc_1based, ref_allele, alt_allele)
    {
        return decisions;
    }
    biallelic_annotation_decisions(evidence, loc_1based, ref_allele, alt_allele)
}

/// Use the retained pre-trim pair only when it is not the emitted representation.
fn pretrim_output_when_different<'a>(
    emitted_ref: &str,
    emitted_alt: &str,
    annotation_output: PretrimAnnotationAlleles<'a>,
) -> PretrimAnnotationAlleles<'a> {
    let (pre_ref, kept_alts) = annotation_output?;
    let same = pre_ref == emitted_ref && kept_alts.len() == 1 && kept_alts[0] == emitted_alt;
    if same {
        None
    } else {
        Some((pre_ref, kept_alts))
    }
}

/// Java `createAlleleMapper` on the pre-trim call: search the full candidate set,
/// including spanning `*`, and count REF/ALT only for the kept pre-trim alleles.
fn pretrim_annotation_decisions(
    evidence: &HcStrandBiasLikelihoods<'_>,
    loc_1based: u64,
    pre_ref: &str,
    kept_alts: &[String],
) -> Option<Vec<AnnotationAlleleDecision>> {
    let cache = build_per_haplotype_variation_events(
        evidence.haplotypes,
        evidence.full_ref_bytes,
        evidence.full_pad_1based,
        evidence.max_mnp_distance,
        evidence.contig,
    );
    let raw =
        variation_events_at_position_from_cache(&cache, loc_1based, evidence.emit_spanning_dels);
    let events = if evidence.emit_spanning_dels {
        replace_span_del_events(
            &raw,
            loc_1based,
            evidence.full_pad_1based,
            evidence.full_ref_bytes,
        )
    } else {
        raw
    };
    let (long_ref, mut alts) = merged_alleles_for_genotyping(&events, loc_1based)?;
    if long_ref != pre_ref || kept_alts.is_empty() {
        return None;
    }
    if kept_alts
        .iter()
        .any(|kept| !alts.iter().any(|alt| alt == kept))
    {
        return None;
    }
    if evidence.emit_spanning_dels && !alts.iter().any(|alt| alt == "*") {
        let spans = (0..evidence.haplotypes.len()).any(|hap| {
            overlapping_events(cache.events_for(hap), loc_1based)
                .iter()
                .any(|event| event.start_1based.get() != loc_1based)
        });
        if spans {
            alts.push("*".to_string());
        }
    }
    let (pools, _, _) = java_create_allele_mapper_pools(
        evidence.haplotypes.len(),
        &cache,
        loc_1based,
        &long_ref,
        &alts,
        evidence.emit_spanning_dels,
    );
    if pools.len() != 1 + alts.len() {
        return None;
    }
    let mut allele_bases = Vec::with_capacity(pools.len());
    allele_bases.push(long_ref.as_str());
    for alt in &alts {
        allele_bases.push(alt.as_str());
    }
    let rows = region_likelihoods_to_rows(evidence.likelihoods, evidence.haplotypes.len());
    let mut priorities = vec![0.0; allele_bases.len()];
    priorities[0] = 1.0;
    let mut out = Vec::with_capacity(rows.len());
    for row in &rows {
        let lls: Vec<f64> = pools
            .iter()
            .map(|pool| max_log10_over_haplotype_indices(&row.haplotype_log10_likelihoods, pool))
            .collect();
        let (best_i, _second_i, best_l, second_l) = search_best_allele(&lls, &priorities);
        let confidence = if best_l == second_l {
            0.0
        } else {
            best_l - second_l
        };
        let bases = allele_bases.get(best_i).copied().unwrap_or("");
        out.push(AnnotationAlleleDecision {
            read_index: row.read_index,
            best_index: best_i,
            informative: confidence > LOG_10_INFORMATIVE_THRESHOLD,
            class: annotation_allele_class_kept(bases, pre_ref, kept_alts),
        });
    }
    Some(out)
}

fn annotation_allele_class_kept(
    bases: &str,
    pre_ref: &str,
    kept_alts: &[String],
) -> AnnotationAlleleClass {
    if bases == pre_ref {
        AnnotationAlleleClass::EmittedRef
    } else if kept_alts.iter().any(|alt| alt == bases) {
        AnnotationAlleleClass::EmittedAlt
    } else {
        AnnotationAlleleClass::Other
    }
}

/// Candidate alleles in haplotype-encounter order, reference first.
#[cfg(test)]
fn annotation_candidate_alleles(
    evidence: &HcStrandBiasLikelihoods<'_>,
    loc_1based: u64,
) -> Option<(String, Vec<String>)> {
    let cache = build_per_haplotype_variation_events(
        evidence.haplotypes,
        evidence.full_ref_bytes,
        evidence.full_pad_1based,
        evidence.max_mnp_distance,
        evidence.contig,
    );
    let events = variation_events_at_position_from_cache(&cache, loc_1based, false);
    let (long_ref, alts) = merged_alleles_for_genotyping(&events, loc_1based)?;
    if alts.len() < 2 {
        return None;
    }
    Some((long_ref, alts))
}

fn multi_allele_annotation_decisions(
    evidence: &HcStrandBiasLikelihoods<'_>,
    loc_1based: u64,
    ref_allele: &str,
    alt_allele: &str,
) -> Option<Vec<AnnotationAlleleDecision>> {
    let cache = build_per_haplotype_variation_events(
        evidence.haplotypes,
        evidence.full_ref_bytes,
        evidence.full_pad_1based,
        evidence.max_mnp_distance,
        evidence.contig,
    );
    let events = variation_events_at_position_from_cache(&cache, loc_1based, false);
    let (long_ref, alts) = merged_alleles_for_genotyping(&events, loc_1based)?;
    if alts.len() < 2 || long_ref != ref_allele || !alts.iter().any(|alt| alt == alt_allele) {
        return None;
    }
    let (pools, _, _) = java_create_allele_mapper_pools(
        evidence.haplotypes.len(),
        &cache,
        loc_1based,
        &long_ref,
        &alts,
        evidence.emit_spanning_dels,
    );
    if pools.len() != 1 + alts.len() {
        return None;
    }
    let mut allele_bases = Vec::with_capacity(pools.len());
    allele_bases.push(long_ref.as_str());
    for alt in &alts {
        allele_bases.push(alt.as_str());
    }
    let rows = region_likelihoods_to_rows(evidence.likelihoods, evidence.haplotypes.len());
    let mut priorities = vec![0.0; allele_bases.len()];
    priorities[0] = 1.0;
    let mut out = Vec::with_capacity(rows.len());
    for row in &rows {
        let lls: Vec<f64> = pools
            .iter()
            .map(|pool| max_log10_over_haplotype_indices(&row.haplotype_log10_likelihoods, pool))
            .collect();
        let (best_i, _second_i, best_l, second_l) = search_best_allele(&lls, &priorities);
        let confidence = if best_l == second_l {
            0.0
        } else {
            best_l - second_l
        };
        let bases = allele_bases.get(best_i).copied().unwrap_or("");
        out.push(AnnotationAlleleDecision {
            read_index: row.read_index,
            best_index: best_i,
            informative: confidence > LOG_10_INFORMATIVE_THRESHOLD,
            class: annotation_allele_class(bases, ref_allele, alt_allele),
        });
    }
    Some(out)
}

/// Previous emitted-pair projection, used when the locus has a single ALT.
fn biallelic_annotation_decisions(
    evidence: &HcStrandBiasLikelihoods<'_>,
    loc_1based: u64,
    ref_allele: &str,
    alt_allele: &str,
) -> Vec<AnnotationAlleleDecision> {
    let event = VariationEvent::from_alleles(evidence.contig, loc_1based, ref_allele, alt_allele);
    let hap_cache = build_per_haplotype_variation_events(
        evidence.haplotypes,
        evidence.full_ref_bytes,
        evidence.full_pad_1based,
        evidence.max_mnp_distance,
        evidence.contig,
    );
    let mapping = create_allele_mapper_with_events(
        &event,
        loc_1based,
        evidence.haplotypes,
        evidence.pad_start_1based,
        evidence.ref_bytes,
        evidence.max_mnp_distance,
        evidence.emit_spanning_dels,
        Some(&hap_cache),
    );
    let rows = region_likelihoods_to_rows(evidence.likelihoods, evidence.haplotypes.len());
    let marg = marginalize_rows_to_biallelic_alleles(
        &rows,
        &mapping.ref_haplotype_indices,
        &mapping.alt_haplotype_indices,
    );
    let mut out = Vec::with_capacity(marg.len());
    for row in &marg {
        let lls = &row.haplotype_log10_likelihoods;
        let ll_ref = lls.first().copied().unwrap_or(f64::NEG_INFINITY);
        let ll_alt = lls.get(1).copied().unwrap_or(f64::NEG_INFINITY);
        if !ll_ref.is_finite() && !ll_alt.is_finite() {
            out.push(AnnotationAlleleDecision {
                read_index: row.read_index,
                best_index: 0,
                informative: false,
                class: AnnotationAlleleClass::Other,
            });
            continue;
        }
        let (best_is_ref, best, second) = if ll_ref > ll_alt {
            (true, ll_ref, ll_alt)
        } else if ll_alt > ll_ref {
            (false, ll_alt, ll_ref)
        } else {
            (true, ll_ref, ll_alt)
        };
        let gap = if second.is_finite() {
            best - second
        } else {
            f64::INFINITY
        };
        out.push(AnnotationAlleleDecision {
            read_index: row.read_index,
            best_index: usize::from(!best_is_ref),
            informative: gap > LOG_10_INFORMATIVE_THRESHOLD,
            class: if best_is_ref {
                AnnotationAlleleClass::EmittedRef
            } else {
                AnnotationAlleleClass::EmittedAlt
            },
        });
    }
    out
}

/// GATK `AlleleLikelihoods.searchBestAllele` plus the default tie priority.
fn search_best_allele(lls: &[f64], priorities: &[f64]) -> (usize, usize, f64, f64) {
    let n = lls.len();
    if n == 0 {
        return (0, 0, f64::NEG_INFINITY, f64::NEG_INFINITY);
    }
    let mut best_i = 0usize;
    let mut second_i = 0usize;
    let mut best_l = lls[0];
    let mut second_l = f64::NEG_INFINITY;
    for i in 1..n {
        let likelihood = lls[i];
        if likelihood > best_l {
            second_i = best_i;
            second_l = best_l;
            best_i = i;
            best_l = likelihood;
        } else if likelihood > second_l {
            second_i = i;
            second_l = likelihood;
        }
    }
    if best_l - second_l < LOG_10_INFORMATIVE_THRESHOLD {
        // The window stays on the pre-tie maximum. A high-priority allele farther
        // than 0.2 below that maximum does not inherit the tie.
        let window_best = best_l;
        let mut tied_best = best_i;
        let mut tied_second = second_i;
        let mut best_p = priorities.get(best_i).copied().unwrap_or(0.0);
        let mut second_p = priorities.get(second_i).copied().unwrap_or(0.0);
        for i in 0..n {
            let likelihood = lls[i];
            if i == tied_best || window_best - likelihood > LOG_10_INFORMATIVE_THRESHOLD {
                continue;
            }
            let priority = priorities.get(i).copied().unwrap_or(0.0);
            if priority > best_p {
                second_p = best_p;
                tied_second = tied_best;
                best_p = priority;
                tied_best = i;
            } else if priority > second_p {
                second_p = priority;
                tied_second = i;
            }
        }
        best_i = tied_best;
        second_i = tied_second;
    }
    let best_l = lls[best_i];
    let second_l = if second_i != best_i {
        lls[second_i]
    } else {
        f64::NEG_INFINITY
    };
    (best_i, second_i, best_l, second_l)
}

fn max_log10_over_haplotype_indices(row: &[f64], haps: &[crate::bio_ids::HaplotypeIndex]) -> f64 {
    let mut best = f64::NEG_INFINITY;
    for hap in haps {
        if let Some(&value) = row.get(hap.get()) {
            if value > best {
                best = value;
            }
        }
    }
    best
}

/// 6R.191: ReadPosRankSum `fillQualsFromLikelihood` lists (`getReadPosition`).
pub fn read_pos_rank_sum_quals_from_likelihoods(
    evidence: &HcStrandBiasLikelihoods<'_>,
    loc_1based: u64,
    ref_allele: &str,
    alt_allele: &str,
) -> (Vec<f64>, Vec<f64>) {
    read_pos_rank_sum_quals(evidence, loc_1based, ref_allele, alt_allele, None)
}

fn read_pos_rank_sum_quals(
    evidence: &HcStrandBiasLikelihoods<'_>,
    loc_1based: u64,
    ref_allele: &str,
    alt_allele: &str,
    annotation_output: PretrimAnnotationAlleles<'_>,
) -> (Vec<f64>, Vec<f64>) {
    let vc_end_1based = java_variant_context_end_1based(loc_1based, ref_allele);
    fill_ranksum_quals_from_likelihood(
        evidence,
        loc_1based,
        ref_allele,
        alt_allele,
        annotation_output,
        |rec, vc_start_1based| java_read_pos_rank_element(rec, vc_start_1based, vc_end_1based),
    )
}

/// 6R.193: BaseQRankSum `fillQualsFromLikelihood` lists
/// (`getReadBaseQualityAtReferenceCoordinate`). Same membership as ReadPos.
pub fn baseq_rank_sum_quals_from_likelihoods(
    evidence: &HcStrandBiasLikelihoods<'_>,
    loc_1based: u64,
    ref_allele: &str,
    alt_allele: &str,
) -> (Vec<f64>, Vec<f64>) {
    baseq_rank_sum_quals(evidence, loc_1based, ref_allele, alt_allele, None)
}

fn baseq_rank_sum_quals(
    evidence: &HcStrandBiasLikelihoods<'_>,
    loc_1based: u64,
    ref_allele: &str,
    alt_allele: &str,
    annotation_output: PretrimAnnotationAlleles<'_>,
) -> (Vec<f64>, Vec<f64>) {
    fill_ranksum_quals_from_likelihood(
        evidence,
        loc_1based,
        ref_allele,
        alt_allele,
        annotation_output,
        java_baseq_rank_element,
    )
}

/// 6R.194: MQRankSum `fillQualsFromLikelihood` lists (`read.getMappingQuality()`).
/// Same membership as ReadPos / BaseQ. Does **not** use RMS `sampleEvidence`.
pub fn mq_rank_sum_quals_from_likelihoods(
    evidence: &HcStrandBiasLikelihoods<'_>,
    loc_1based: u64,
    ref_allele: &str,
    alt_allele: &str,
) -> (Vec<f64>, Vec<f64>) {
    mq_rank_sum_quals(evidence, loc_1based, ref_allele, alt_allele, None)
}

fn mq_rank_sum_quals(
    evidence: &HcStrandBiasLikelihoods<'_>,
    loc_1based: u64,
    ref_allele: &str,
    alt_allele: &str,
    annotation_output: PretrimAnnotationAlleles<'_>,
) -> (Vec<f64>, Vec<f64>) {
    fill_ranksum_quals_from_likelihood(
        evidence,
        loc_1based,
        ref_allele,
        alt_allele,
        annotation_output,
        java_mq_rank_element,
    )
}

/// GATK `BaseQualityRankSumTest.getElementForRead`.
///
/// `ReadUtils.getReadBaseQualityAtReferenceCoordinate` on the attached
/// alignment. An empty lookup contributes no observation. Java does not
/// reread a pre-realignment alignment.
fn java_baseq_rank_element(rec: &Record, vc_start_1based: i32) -> Option<f64> {
    read_base_quality_at_ref_coord_1based(rec, vc_start_1based).map(f64::from)
}

/// GATK `MappingQualityRankSumTest.getElementForRead`: `read.getMappingQuality()`.
fn java_mq_rank_element(rec: &Record, _vc_start_1based: i32) -> Option<f64> {
    Some(f64::from(rec.mapq()))
}

/// GATK `RMSMappingQuality.calculateRawData` membership:
/// unique reads remaining in `likelihoods.sampleEvidence` after
/// `filterPoorlyModeledEvidence`, with `MQ != QualityUtils.MAPPING_QUALITY_UNAVAILABLE`.
///
/// This is **not** the full `genotyping_reads` overlap list, **not** an ALT pileup,
/// and **not** the informative-best-allele subset used by FS/SOR.
/// MQ=0 is included; MQ=255 is skipped. Duplicates / secondary / supplementary
/// are already absent from this post-filter evidence list.
pub fn rms_mapping_quality_sample_reads<'a>(
    reads: &'a [SharedBamRecord],
    likelihoods: &[RegionReadLikelihood],
) -> Vec<&'a SharedBamRecord> {
    let mut seen = BTreeSet::new();
    let mut out = Vec::new();
    for cell in likelihoods {
        let idx = cell.read_index.get();
        if !seen.insert(idx) {
            continue;
        }
        if let Some(rec) = reads.get(idx) {
            out.push(rec);
        }
    }
    out
}

pub fn rms_mapping_quality_sample_mapqs(
    reads: &[SharedBamRecord],
    likelihoods: &[RegionReadLikelihood],
) -> Vec<u8> {
    rms_mapping_quality_sample_reads(reads, likelihoods)
        .into_iter()
        .map(|rec| rec.mapq())
        .filter(|&mq| mq != MAPPING_QUALITY_UNAVAILABLE)
        .collect()
}

/// GATK `RMSMappingQuality.calculateRawData` squares + count, then
/// `makeFinalizedAnnotationString`: `Math.sqrt(sumOfSquaredMQs / (double) numOfReads)`.
///
/// Empty eligible list returns `None` (caller keeps the pre-existing MQ=0.0 emit).
pub fn rms_mapping_quality_raw(mqs: &[u8]) -> Option<(usize, u64, f64)> {
    if mqs.is_empty() {
        return None;
    }
    let n = mqs.len();
    let sum_sq: u64 = mqs
        .iter()
        .map(|&m| {
            let mq = u64::from(m);
            mq * mq
        })
        .sum();
    let rms = (sum_sq as f64 / n as f64).sqrt();
    Some((n, sum_sq, rms))
}

/// Java `RMSMappingQuality.OUTPUT_PRECISION` (`String.format("%.2f", rms)`).
fn java_finalize_mq(rms: f64) -> f64 {
    (rms * 100.0).round() / 100.0
}

/// GATK 4.4 `Coverage.annotate`: `likelihoods.evidenceCount()` on the
/// already-canonical annotation AlleleLikelihoods passed to
/// `VariantAnnotatorEngine.annotateContext`.
///
/// Java `evidenceCount()` is remaining-list cardinality
/// (`evidenceBySampleIndex.stream().mapToInt(List::size).sum()`). Loc-loop
/// `retainEvidence` / overlap is **object construction**, not a second
/// Coverage filter (6R.201). Unique `read_index` is the attached
/// unique-evidence cardinality. This is **not** FORMAT/`DepthPerSampleHC`
/// (informative `bestAlleles` only) and **not** the unfiltered region
/// `sampleEvidence` used by MQ.
///
/// Default HC `prepareReadAlleleLikelihoodsForAnnotation` also
/// `addEvidence(overlappingFilteredReads, 0)`. Reads already remaining in the
/// post-filter matrix are not double-counted. INFO DP may exceed FORMAT DP.
pub fn coverage_evidence_count(
    reads: &[SharedBamRecord],
    likelihoods: &[RegionReadLikelihood],
    start_1based: u64,
    end_1based: u64,
    margin: i32,
) -> i32 {
    let mut seen = BTreeSet::new();
    for cell in likelihoods {
        seen.insert(cell.read_index.get());
    }
    // 7R.6: mate-on-other-contig reads are omitted from PairHMM and have no
    // likelihood row. Java Coverage still counts the overlapping ones.
    crate::engine::engine_observe::with_coverage_only_failed_mates(|mates| {
        for (i, rec) in mates.iter().enumerate() {
            if !crate::hc_genotyping_engine::java_alignment_read_overlaps_interval(
                rec,
                start_1based,
                end_1based,
                margin,
            ) {
                continue;
            }
            seen.insert(reads.len().saturating_add(i));
        }
    });
    seen.len() as i32
}

/// 6R.168: Java RMS on [`rms_mapping_quality_sample_mapqs`], then `%.2f`.
///
/// 7R.9: Java `sampleEvidence` after `addEvidence` also contains overlapping
/// mate-filtered reads. `with_coverage_only_failed_mates` is exactly the reads
/// the 6R.311 retain removed from the scored list, so they have no likelihood
/// index and are not already in `mqs`.
fn mq_rms_of_sample_evidence(
    reads: &[SharedBamRecord],
    likelihoods: &[RegionReadLikelihood],
    start_1based: u64,
    end_1based: u64,
    margin: i32,
) -> f64 {
    let mut mqs = rms_mapping_quality_sample_mapqs(reads, likelihoods);
    crate::engine::engine_observe::with_coverage_only_failed_mates(|mates| {
        for rec in mates {
            if rec.mapq() == MAPPING_QUALITY_UNAVAILABLE {
                continue;
            }
            if crate::hc_genotyping_engine::java_alignment_read_overlaps_interval(
                rec,
                start_1based,
                end_1based,
                margin,
            ) {
                mqs.push(rec.mapq());
            }
        }
    });
    match rms_mapping_quality_raw(&mqs) {
        Some((_, _, rms)) => java_finalize_mq(rms),
        None => 0.0,
    }
}

pub fn annotate_hc_variant_site(
    region: Option<&AssemblyRegion>,
    position_1based: u64,
    ref_allele: &str,
    alt_allele: &str,
    genotype: &RegionGenotypeResult,
    config: &HcGenotypingConfig,
    qual_log10_p_error: Option<f64>,
    strand_likelihoods: Option<&HcStrandBiasLikelihoods<'_>>,
    annotation_output: Option<(&str, &[String])>,
) -> GatkResult<HcVariantSiteAnnotations> {
    let gl_for_qual = genotype_log10_likelihoods_after_java_genotype_pl_roundtrip(
        &genotype.genotype_log10_likelihoods,
    );
    // 6R.196: Java AF alt Dirichlet weight is SNP iff `alt.length()==ref.length()`.
    let af_cfg = AfCalculatorConfig::default();
    let alt_pc = biallelic_alt_pseudocount(ref_allele, alt_allele, &af_cfg);
    let af_result =
        calculate_biallelic_af_em_with_alt_pseudocount(&[&gl_for_qual], &af_cfg, alt_pc)?;
    let qual = match qual_log10_p_error {
        Some(log10_pe) => (-10.0 * log10_pe.min(0.0)) + 0.0,
        None => (-10.0 * af_result.log10_posterior_no_variant) + 0.0,
    };
    let best_idx =
        crate::genotyping::biallelic_genotype_index_from_pl(&genotype.format.pl).as_usize();
    let gt = genotype_from_index(best_idx);
    // 6R.174: INFO DP is Java `Coverage.evidenceCount`, not FORMAT/`DepthPerSampleHC`.
    // 6R.180: `strand_likelihoods` is the per-variant genotyping AlleleLikelihoods when
    // present; formulas below are unchanged.
    let var_end = VariationEvent::vcf_end_1based(position_1based, ref_allele);
    let dp = match strand_likelihoods {
        Some(ev) => coverage_evidence_count(
            ev.reads,
            ev.likelihoods,
            position_1based,
            var_end,
            config.informative_read_overlap_margin,
        ),
        None => genotype.format.dp.as_i32().max(0),
    };
    let qd_depth = qd_depth_for_variant(&gt, &genotype.format);
    // Java AC/AF/AN from called genotypes; MLEAC/MLEAF from AFCalculationResult
    // (`composeCallAttributes`: round(EM alt count), MLEAF = MLEAC / AN).
    let (ac, af, an) = called_ac_af_an(best_idx);
    let mleac = af_result.alt_allele_count;
    let mleaf = if an > 0 {
        (mleac as f64 / an as f64).min(1.0)
    } else {
        0.0
    };
    // 6R.191: ReadPosRankSum uses Java `fillQualsFromLikelihood` on the same
    // per-variant annotation `AlleleLikelihoods` as DP/MQ/SOR (`strand_likelihoods`).
    // 6R.193: BaseQRankSum uses that same object and membership; only the
    // element is base quality. Do not walk `region.reads` pileup.
    // 6R.194: MQRankSum uses that same object and membership; element is
    // `read.getMappingQuality()`. Do not use RMS `sampleEvidence`.
    let _ = region;
    let (ref_positions, alt_positions, ref_baseq, alt_baseq, ref_mapq, alt_mapq) =
        match strand_likelihoods {
            Some(ev) => {
                let (rp_ref, rp_alt) = match annotation_output {
                    None => read_pos_rank_sum_quals_from_likelihoods(
                        ev,
                        position_1based,
                        ref_allele,
                        alt_allele,
                    ),
                    Some(output) => read_pos_rank_sum_quals(
                        ev,
                        position_1based,
                        ref_allele,
                        alt_allele,
                        Some(output),
                    ),
                };
                let (bq_ref, bq_alt) = match annotation_output {
                    None => baseq_rank_sum_quals_from_likelihoods(
                        ev,
                        position_1based,
                        ref_allele,
                        alt_allele,
                    ),
                    Some(output) => baseq_rank_sum_quals(
                        ev,
                        position_1based,
                        ref_allele,
                        alt_allele,
                        Some(output),
                    ),
                };
                let (mq_ref, mq_alt) = match annotation_output {
                    None => mq_rank_sum_quals_from_likelihoods(
                        ev,
                        position_1based,
                        ref_allele,
                        alt_allele,
                    ),
                    Some(output) => {
                        mq_rank_sum_quals(ev, position_1based, ref_allele, alt_allele, Some(output))
                    }
                };
                (rp_ref, rp_alt, bq_ref, bq_alt, mq_ref, mq_alt)
            }
            None => (
                Vec::new(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
            ),
        };
    // 6R.166: FS/SOR use Java `getContingencyTable` on post-filter allele likelihoods.
    // 6R.167: MQ membership is Java `sampleEvidence` (MQ != 255), including
    // overlapping coverage-only failed mates (7R.9).
    // 6R.168: MQ aggregation is Java `makeFinalizedAnnotationString` RMS + `%.2f`.
    let (fs_rf, fs_rr, fs_af, fs_ar, sor_rf, sor_rr, sor_af, sor_ar, mq) = match strand_likelihoods
    {
        Some(ev) => {
            let sample = strand_bias_sample_counts(
                ev,
                position_1based,
                ref_allele,
                alt_allele,
                annotation_output,
            );
            let fs = apply_strand_min_count(sample, FISHER_STRAND_MIN_COUNT);
            let sor = apply_strand_min_count(sample, STRAND_ODDS_RATIO_MIN_COUNT);
            (
                fs.0,
                fs.1,
                fs.2,
                fs.3,
                sor.0,
                sor.1,
                sor.2,
                sor.3,
                mq_rms_of_sample_evidence(
                    ev.reads,
                    ev.likelihoods,
                    position_1based,
                    var_end,
                    config.informative_read_overlap_margin,
                ),
            )
        }
        None => (0, 0, 0, 0, 0, 0, 0, 0, 0.0),
    };
    let fs = fisher_strand::fisher_strand_statistic(fs_rf, fs_rr, fs_af, fs_ar);
    let sor = strand_odds_ratio::strand_odds_ratio(sor_rf, sor_rr, sor_af, sor_ar);
    let rp = read_pos_rank_sum::read_pos_rank_sum(&ref_positions, &alt_positions);
    let bq = rank_sum_baseq::base_quality_rank_sum(&ref_baseq, &alt_baseq);
    let mqrs = mapping_quality_rank_sum::mapping_quality_rank_sum(&ref_mapq, &alt_mapq);
    // Raw QUAL/depth only. `fixTooHighQD` consumes the process-global Java RNG and
    // must run later in genomic emit order (see `apply_fix_too_high_qd_to_vcf_records`).
    let qd = qual_by_depth::raw_qual_by_depth(qual, qd_depth);
    let (ref_n, het_n, hom_alt_n) = genotype_counts_from_index(best_idx); // usize diploid index
    let excess_het_phred = excess_het::excess_heterozygosity_phred(ref_n, het_n, hom_alt_n);
    let inbreeding_coeff = if ref_n + het_n + hom_alt_n > 0 {
        1.0 - (het_n as f64) / (ref_n + het_n + hom_alt_n) as f64
    } else {
        0.0
    };
    Ok(HcVariantSiteAnnotations {
        qual,
        ac,
        af,
        an,
        dp,
        excess_het: excess_het_phred,
        fs,
        mleac,
        mleaf,
        mq,
        qd,
        sor,
        read_pos_rank_sum: rp,
        base_q_rank_sum: bq,
        mq_rank_sum: mqrs,
        inbreeding_coeff,
    })
}

#[allow(dead_code)] // 6R.170 pileup RankSum; production moved in 6R.191.
fn read_offset_evidence_at_site(
    region: Option<&AssemblyRegion>,
    position_1based: u64,
    ref_allele: &str,
    alt_allele: &str,
) -> (Vec<f64>, Vec<f64>) {
    let mut ref_pos = Vec::new();
    let mut alt_pos = Vec::new();
    let ref_b = ref_allele.as_bytes().first().copied().unwrap_or(b'N');
    let alt_b = alt_allele.as_bytes().first().copied().unwrap_or(b'N');
    let pos = position_1based as i32;
    let Some(region) = region else {
        return (ref_pos, alt_pos);
    };
    for rec in &region.reads {
        let read_start = rec.pos() + 1;
        let Some(base) = read_base_at_ref_coord_1based(rec, pos) else {
            continue;
        };
        let offset = (position_1based as i64 - read_start) as f64;
        if base.eq_ignore_ascii_case(&alt_b) {
            alt_pos.push(offset);
        } else if base.eq_ignore_ascii_case(&ref_b) {
            ref_pos.push(offset);
        }
    }
    (ref_pos, alt_pos)
}

fn called_ac_af_an(best: usize) -> (i32, f64, i32) {
    let an = 2;
    match best {
        0 => (0, 0.0, an),
        1 => (1, 0.5, an),
        _ => (2, 1.0, an),
    }
}

#[allow(dead_code)] // 6R.38/40: former MLEAC path (called GT); Java uses AF MLE instead.
fn mle_alleles_from_genotype_index(best: usize) -> (i32, f64, i32, i32, f64) {
    let (ac, af, an) = called_ac_af_an(best);
    (ac, af, an, ac, af)
}

fn genotype_counts_from_index(best: usize) -> (u32, u32, u32) {
    match best {
        0 => (1, 0, 0),
        1 => (0, 1, 0),
        _ => (0, 0, 1),
    }
}

#[allow(dead_code)] // 6R.165 ALT-pileup MQ walk; production MQ moved in 6R.167.
fn read_strand_evidence_at_site(
    region: Option<&AssemblyRegion>,
    position_1based: u64,
    ref_allele: &str,
    alt_allele: &str,
) -> (u32, u32, u32, u32, u64, u32) {
    let mut ref_fw = 0u32;
    let mut ref_rv = 0u32;
    let mut alt_fw = 0u32;
    let mut alt_rv = 0u32;
    let mut mq_sum = 0u64;
    let mut mq_n = 0u32;
    let ref_b = ref_allele.as_bytes().first().copied().unwrap_or(b'N');
    let alt_b = alt_allele.as_bytes().first().copied().unwrap_or(b'N');
    let pos = position_1based as i32;
    let Some(region) = region else {
        return (ref_fw, ref_rv, alt_fw, alt_rv, mq_sum, mq_n);
    };
    for rec in &region.reads {
        let Some(base) = read_base_at_ref_coord_1based(rec, pos) else {
            continue;
        };
        let reverse = rec.flags() & FLAG_REVERSE != 0;
        let supports_alt = base.eq_ignore_ascii_case(&alt_b);
        let supports_ref = base.eq_ignore_ascii_case(&ref_b);
        if supports_alt {
            if reverse {
                alt_rv += 1;
            } else {
                alt_fw += 1;
            }
            mq_sum += u64::from(rec.mapq());
            mq_n += 1;
        } else if supports_ref {
            if reverse {
                ref_rv += 1;
            } else {
                ref_fw += 1;
            }
        }
    }
    (ref_fw, ref_rv, alt_fw, alt_rv, mq_sum, mq_n)
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

#[cfg(test)]
#[path = "tests/variant_site_hc_annotations_tests.rs"]
mod tests;
