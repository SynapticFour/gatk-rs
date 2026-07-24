//! GATK HC variant-site INFO / QUAL for assembly-region VCF emission (J-D01).

use crate::activity_scoring::genotype_log10_likelihoods_after_java_genotype_pl_roundtrip;
use crate::af_calc::{calculate_biallelic_af_em, AfCalculatorConfig};
use crate::annotator::plugins::{
    excess_het, fisher_strand, qual_by_depth, read_pos_rank_sum, strand_odds_ratio,
};
use crate::assembly_region_iterator::AssemblyRegion;
use crate::fragment_overlap::read_base_at_ref_coord_1based;
use crate::genotyping::GenotypeFormatFields;
use crate::hc_genotyping_engine::{HcGenotypingConfig, RegionGenotypeResult};
use gatk_common::GatkResult;
use gatk_core::io::vcf::Genotype;

const FLAG_REVERSE: u16 = 0x10;

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
    pub read_pos_rank_sum: f64,
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

pub fn annotate_hc_variant_site(
    region: Option<&AssemblyRegion>,
    position_1based: u64,
    ref_allele: &str,
    alt_allele: &str,
    genotype: &RegionGenotypeResult,
    _config: &HcGenotypingConfig,
) -> GatkResult<HcVariantSiteAnnotations> {
    let gl_for_qual = genotype_log10_likelihoods_after_java_genotype_pl_roundtrip(
        &genotype.genotype_log10_likelihoods,
    );
    let qual = qual_from_af_calculation(&gl_for_qual)?;
    let best_idx =
        crate::genotyping::biallelic_genotype_index_from_pl(&genotype.format.pl).as_usize();
    let gt = genotype_from_index(best_idx);
    let dp = genotype.format.dp.as_i32().max(0);
    let qd_depth = qd_depth_for_variant(&gt, &genotype.format);
    let (ac, af, an, mleac, mleaf) = mle_alleles_from_genotype_index(best_idx);
    let (ref_fw, ref_rv, alt_fw, alt_rv, mq_sum, mq_n) =
        read_strand_evidence_at_site(region, position_1based, ref_allele, alt_allele);
    let (ref_positions, alt_positions) =
        read_offset_evidence_at_site(region, position_1based, ref_allele, alt_allele);
    let fs = fisher_strand::fisher_strand_statistic(ref_fw, ref_rv, alt_fw, alt_rv);
    let sor = strand_odds_ratio::strand_odds_ratio(ref_fw, ref_rv, alt_fw, alt_rv);
    let rp = read_pos_rank_sum::read_pos_rank_sum(&ref_positions, &alt_positions);
    let qd = qual_by_depth::qual_by_depth(qual, qd_depth);
    let mq = if mq_n > 0 {
        mq_sum as f64 / mq_n as f64
    } else {
        0.0
    };
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

fn mle_alleles_from_genotype_index(best: usize) -> (i32, f64, i32, i32, f64) {
    let an = 2;
    match best {
        0 => (0, 0.0, an, 0, 0.0),
        1 => (1, 0.5, an, 1, 0.5),
        _ => (2, 1.0, an, 2, 1.0),
    }
}

fn genotype_counts_from_index(best: usize) -> (u32, u32, u32) {
    match best {
        0 => (1, 0, 0),
        1 => (0, 1, 0),
        _ => (0, 0, 1),
    }
}

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
        assert!((qual - 2224.06).abs() < 0.1, "qual={qual}");
    }
}
