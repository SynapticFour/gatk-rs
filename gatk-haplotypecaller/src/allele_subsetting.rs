//! GATK haplotype-score allele trimming (`whichAllelesToKeepBasedonHapScores`).

use crate::genotype_limits::which_alleles_to_keep_by_haplotype_scores_with_ref;
use crate::genotyping::HaplotypeLikelihoodAggregation;
use crate::haplotype::Haplotype;
use crate::hc_genotyping_engine::subset_biallelic_haplotype_indices;
use gatk_common::GatkResult;

/// Pick alleles to retain when genotype enumeration would explode (G.3 + G-D05).
pub fn subset_alleles_for_genotyping(
    haplotypes: &[Haplotype],
    aggregation: &HaplotypeLikelihoodAggregation,
    max_allele_count: usize,
) -> GatkResult<Vec<usize>> {
    if haplotypes.len() <= max_allele_count {
        return Ok((0..haplotypes.len()).collect());
    }
    // One score per haplotype/allele — borrow as `&[f64]` rows (no `Vec<Vec<f64>>`).
    let scores: Vec<&[f64]> = aggregation
        .haplotype_log10_sums
        .iter()
        .map(std::slice::from_ref)
        .collect();
    let is_ref: Vec<bool> = haplotypes.iter().map(|h| h.is_reference).collect();
    Ok(which_alleles_to_keep_by_haplotype_scores_with_ref(
        &scores,
        max_allele_count,
        Some(&is_ref),
    ))
}

/// Biallelic ref+top-alt subsetting (existing G.1.1) exposed for dumps.
pub fn biallelic_ref_alt_indices(
    aggregation: &HaplotypeLikelihoodAggregation,
    haplotypes: &[Haplotype],
) -> (usize, usize) {
    subset_biallelic_haplotype_indices(aggregation, haplotypes)
}
