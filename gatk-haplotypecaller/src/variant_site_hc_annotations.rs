//! GATK HC variant-site INFO / QUAL for assembly-region VCF emission (J-D01).

use crate::activity_scoring::genotype_log10_likelihoods_after_java_genotype_pl_roundtrip;
use crate::af_calc::{
    calculate_biallelic_af_em, diploid_af_log10_prob_only_ref_allele_exists, AfCalculatorConfig,
};
use crate::annotator::plugins::{
    excess_het, fisher_strand, qual_by_depth, read_pos_rank_sum, strand_odds_ratio,
};
use crate::assembly_region_iterator::AssemblyRegion;
use crate::event_map::{build_per_haplotype_variation_events, VariationEvent};
use crate::fragment_overlap::read_base_at_ref_coord_1based;
use crate::genotyping::GenotypeFormatFields;
use crate::haplotype::Haplotype;
use crate::hc_allele_mapping::create_allele_mapper_with_events;
use crate::hc_genotyping_engine::{
    java_alignment_read_overlaps_interval, marginalize_rows_to_biallelic_alleles,
    region_likelihoods_to_rows, HcGenotypingConfig, RegionGenotypeResult,
};
use crate::read_model::MAPPING_QUALITY_UNAVAILABLE;
use crate::read_realignment::LOG_10_INFORMATIVE_THRESHOLD;
use crate::region_read_likelihood::RegionReadLikelihood;
use crate::shared_bam::SharedBamRecord;
use gatk_common::GatkResult;
use gatk_core::io::vcf::Genotype;
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
/// Per-read: remarg to site REF/ALT haplotype pools, `bestAllelesBreakingTies`
/// (REF wins likelihood ties), then `isInformative` (`confidence > 0.2`).
/// Best alleles that are neither REF nor an ALT of this site are skipped.
/// Sample table is copied only when total counts `> min_count`.
pub fn strand_bias_contingency_table(
    evidence: &HcStrandBiasLikelihoods<'_>,
    loc_1based: u64,
    ref_allele: &str,
    alt_allele: &str,
    min_count: u32,
) -> (u32, u32, u32, u32) {
    apply_strand_min_count(
        strand_bias_sample_counts(evidence, loc_1based, ref_allele, alt_allele),
        min_count,
    )
}

fn strand_bias_sample_counts(
    evidence: &HcStrandBiasLikelihoods<'_>,
    loc_1based: u64,
    ref_allele: &str,
    alt_allele: &str,
) -> (u32, u32, u32, u32) {
    if evidence.likelihoods.is_empty() || evidence.haplotypes.is_empty() {
        return (0, 0, 0, 0);
    }
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
    let mut ref_fw = 0u32;
    let mut ref_rv = 0u32;
    let mut alt_fw = 0u32;
    let mut alt_rv = 0u32;
    for row in &marg {
        let Some(rec) = evidence.reads.get(row.read_index) else {
            continue;
        };
        let lls = &row.haplotype_log10_likelihoods;
        let ll_ref = lls.first().copied().unwrap_or(f64::NEG_INFINITY);
        let ll_alt = lls.get(1).copied().unwrap_or(f64::NEG_INFINITY);
        if !ll_ref.is_finite() && !ll_alt.is_finite() {
            continue;
        }
        let (best_is_ref, best, second) = if ll_ref > ll_alt {
            (true, ll_ref, ll_alt)
        } else if ll_alt > ll_ref {
            (false, ll_alt, ll_ref)
        } else {
            // Java `bestAllelesBreakingTies`: REF priority 1.0 vs ALT 0 on a tie.
            (true, ll_ref, ll_alt)
        };
        let gap = if second.is_finite() {
            best - second
        } else {
            f64::INFINITY
        };
        // Java `BestAllele.isInformative`: strict `confidence > 0.2`. Near-ties do not count.
        if gap <= LOG_10_INFORMATIVE_THRESHOLD {
            continue;
        }
        let reverse = rec.flags() & FLAG_REVERSE != 0;
        if best_is_ref {
            if reverse {
                ref_rv += 1;
            } else {
                ref_fw += 1;
            }
        } else if reverse {
            alt_rv += 1;
        } else {
            alt_fw += 1;
        }
    }
    (ref_fw, ref_rv, alt_fw, alt_rv)
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
/// allele-level object passed to `VariantAnnotatorEngine.annotateContext`.
///
/// That object is genotyping `retainEvidence` (overlap with
/// `variantCallingRelevantOverlap` = merged VC ± `informativeReadOverlapMargin`).
/// This is **not** FORMAT/`DepthPerSampleHC` (informative `bestAlleles` only)
/// and **not** the unfiltered region `sampleEvidence` used by MQ.
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
        let idx = cell.read_index.get();
        if !seen.insert(idx) {
            continue;
        }
        let Some(rec) = reads.get(idx) else {
            seen.remove(&idx);
            continue;
        };
        if !java_alignment_read_overlaps_interval(rec, start_1based, end_1based, margin) {
            seen.remove(&idx);
        }
    }
    seen.len() as i32
}

/// 6R.168: Java RMS on [`rms_mapping_quality_sample_mapqs`], then `%.2f`.
fn mq_rms_of_sample_evidence(
    reads: &[SharedBamRecord],
    likelihoods: &[RegionReadLikelihood],
) -> f64 {
    let mqs = rms_mapping_quality_sample_mapqs(reads, likelihoods);
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
) -> GatkResult<HcVariantSiteAnnotations> {
    let gl_for_qual = genotype_log10_likelihoods_after_java_genotype_pl_roundtrip(
        &genotype.genotype_log10_likelihoods,
    );
    let af_result = calculate_biallelic_af_em(&[&gl_for_qual], &AfCalculatorConfig::default())?;
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
    let (ref_positions, alt_positions) =
        read_offset_evidence_at_site(region, position_1based, ref_allele, alt_allele);
    // 6R.166: FS/SOR use Java `getContingencyTable` on post-filter allele likelihoods.
    // 6R.167: MQ membership is Java `sampleEvidence` (MQ != 255). Do not change it.
    // 6R.168: MQ aggregation is Java `makeFinalizedAnnotationString` RMS + `%.2f`.
    let (fs_rf, fs_rr, fs_af, fs_ar, sor_rf, sor_rr, sor_af, sor_ar, mq) = match strand_likelihoods
    {
        Some(ev) => {
            let sample = strand_bias_sample_counts(ev, position_1based, ref_allele, alt_allele);
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
                mq_rms_of_sample_evidence(ev.reads, ev.likelihoods),
            )
        }
        None => (0, 0, 0, 0, 0, 0, 0, 0, 0.0),
    };
    let fs = fisher_strand::fisher_strand_statistic(fs_rf, fs_rr, fs_af, fs_ar);
    let sor = strand_odds_ratio::strand_odds_ratio(sor_rf, sor_rr, sor_af, sor_ar);
    let rp = read_pos_rank_sum::read_pos_rank_sum(&ref_positions, &alt_positions);
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
        inbreeding_coeff,
    })
}

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
mod tests {
    use super::*;

    #[test]
    fn p11_pl_yields_java_af_qual() {
        let pl = [2238_i32, 151, 0];
        let gl: Vec<f64> = pl.iter().map(|p| -((p - pl[2]) as f64) / 10.0).collect();
        let qual = qual_from_af_calculation(&gl).expect("qual");
        assert!((qual - 2224.06).abs() < 0.15, "qual={qual}");
    }

    #[test]
    fn six_r64_java_emitted_pl_does_not_reproduce_java_qual() {
        // Java VCF PL=542,0,1353 QUAL=510.06. Rust AFCalculator on those PLs is ~534.64,
        // so Java QUAL is not AF(emitted biallelic PL). GenotypingEngine computes QUAL
        // from AF on the pre-subset merged VC (6 GLs), then copies it through subset+trim.
        let gl = [-54.2, 0.0, -135.3];
        let qual = qual_from_af_calculation(&gl).expect("qual");
        assert!(
            (qual - 534.64).abs() < 0.15,
            "Rust AF on Java emitted PL: got {qual}"
        );
        assert!(
            (qual - 510.06).abs() > 10.0,
            "must not confuse AF(emitted PL) with Java site QUAL 510.06, got {qual}"
        );
    }

    #[test]
    fn six_r64_rust_pl_298_0_1103_qual_matches_emitted() {
        // Rust 6R.63 VCF: PL=298,0,1103 QUAL=290.64.
        let gl = [-29.8, 0.0, -110.3];
        let qual = qual_from_af_calculation(&gl).expect("qual");
        assert!(
            (qual - 290.64).abs() < 0.15,
            "Rust PL→QUAL: got {qual}, expected ~290.64"
        );
    }

    #[test]
    fn six_r64_rust_merged_6gl_qual_is_not_java_510() {
        // AF on merged 6-state GLs (Java QUAL source stage). If this were ~510, QUAL
        // would be a post-subset AF-order issue. It is not — 6-state inputs already differ.
        let gl6 = [-29.8, -33.7, -162.0, 0.0, -105.8, -110.3];
        let af = crate::af_calc::calculate_multiallelic_af_em(
            &[&gl6],
            3,
            &AfCalculatorConfig::default(),
        )
        .expect("af6");
        let qual = (-10.0 * af.log10_posterior_no_variant) + 0.0;
        // Measured: same as AF on the subsetted 3-GLs (290.64). Java lifecycle F
        // (AF before unused-ALT subset) would not produce 510.06 from these 6-GLs.
        assert!(
            (qual - 290.64).abs() < 0.15,
            "Rust 6-GL AF QUAL must match subsetted QUAL ~290.64, not Java 510.06; got {qual}"
        );
    }

    #[test]
    fn six_r40_pl_90_6_0_qual_mleac_match_java_af() {
        let gl = [-9.0, -0.6, 0.0];
        let af = calculate_biallelic_af_em(&[&gl], &AfCalculatorConfig::default()).expect("af");
        let qual = (-10.0 * af.log10_posterior_no_variant) + 0.0;
        assert_eq!(af.alt_allele_count, 1);
        assert!((qual - 78.32).abs() < 0.02, "qual={qual}");
        let mleaf = (af.alt_allele_count as f64 / 2.0).min(1.0);
        assert!((mleaf - 0.5).abs() < 1e-9);
        let _qd = crate::annotator::plugins::qual_by_depth::hold_process_qd_rng_for_test();
        crate::annotator::plugins::qual_by_depth::reset_gatk_qual_by_depth_rng();
        let qd = crate::annotator::plugins::qual_by_depth::qual_by_depth(qual, 2);
        assert!(
            (qd - 25.36).abs() < 0.005,
            "QD first GATK-seed gaussian, qd={qd}"
        );
    }
}
