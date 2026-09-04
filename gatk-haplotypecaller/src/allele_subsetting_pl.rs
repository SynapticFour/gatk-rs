//! GATK `AlleleSubsettingUtils#subsetAlleles` PL/AD slice (G-D05 / `g-subset-pl`).

use gatk_common::{GatkError, GatkResult};

/// GATK `MathUtils.scaleLogSpaceArrayForNumericalStability`.
fn scale_log_space_array_for_numerical_stability(gl: &[f64]) -> Vec<f64> {
    if gl.is_empty() {
        return Vec::new();
    }
    let max_v = gl.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    gl.iter().map(|g| g - max_v).collect()
}

/// Diploid genotype pairs in GATK order: for each `j`, all `(i, j)` with `i <= j`.
fn diploid_genotype_pairs(allele_count: usize) -> Vec<(usize, usize)> {
    let mut pairs = Vec::new();
    for j in 0..allele_count {
        for i in 0..=j {
            pairs.push((i, j));
        }
    }
    pairs
}

fn num_genotype_likelihoods(allele_count: usize, ploidy: usize) -> usize {
    debug_assert!(ploidy == 2, "parity v1 supports diploid only");
    allele_count.saturating_mul(allele_count + 1) / 2
}

/// GATK `AlleleSubsettingUtils.subsettedPLIndices` for diploid genotypes.
fn subsetted_pl_indices(
    ploidy: usize,
    original_allele_count: usize,
    keep_allele_indices: &[usize],
) -> Vec<usize> {
    debug_assert_eq!(ploidy, 2);
    let old_pairs = diploid_genotype_pairs(original_allele_count);
    let new_pairs = diploid_genotype_pairs(keep_allele_indices.len());
    let old_index: std::collections::HashMap<(usize, usize), usize> = old_pairs
        .iter()
        .enumerate()
        .map(|(idx, &(i, j))| ((i, j), idx))
        .collect();
    new_pairs
        .iter()
        .map(|&(ni, nj)| {
            let oi = keep_allele_indices[ni];
            let oj = keep_allele_indices[nj];
            let (oi, oj) = if oi <= oj { (oi, oj) } else { (oj, oi) };
            *old_index.get(&(oi, oj)).expect("old genotype index")
        })
        .collect()
}

fn subset_ad(original_ad: &[i32], keep_allele_indices: &[usize]) -> Vec<i32> {
    keep_allele_indices
        .iter()
        .map(|&i| original_ad[i])
        .collect()
}

/// GATK `calculateOutputAlleleSubset` keep-set for a single-sample diploid call.
///
/// Java keeps ALT `a` when `AFCalculationResult.passesThreshold(a, stand-call-conf)`
/// (`GenotypingEngine.calculateOutputAlleleSubset`, 4.4.0.0). Default HC
/// `forceKeepAllele` is false (not GVCF). When one genotype posterior dominates
/// (USE_PLS_TO_ASSIGN), an ALT absent from that genotype has log10 P(absent) ≈ 0
/// and fails the threshold — equivalent to keeping REF plus ALTs present in GT,
/// in original allele order.
pub fn output_allele_keep_indices_from_assigned_gt(
    n_alleles: usize,
    gt_allele_indices: &[i32],
) -> Vec<usize> {
    let mut keep = Vec::with_capacity(n_alleles);
    keep.push(0);
    for i in 1..n_alleles {
        if gt_allele_indices
            .iter()
            .any(|&g| g >= 0 && (g as usize) == i)
        {
            keep.push(i);
        }
    }
    keep
}

/// Diploid unused-ALT subset after merged genotyping (6R.62).
///
/// PLs are remapped with [`subset_log10_genotype_likelihoods`] (`AlleleSubsettingUtils.subsettedPLIndices`
/// + `scaleLogSpaceArrayForNumericalStability`). AD is sliced by keep indices.
/// GLs are **not** recalculated from reads.
#[derive(Debug, Clone, PartialEq)]
pub struct UnusedAltSubsetResult {
    pub alt_alleles: Vec<String>,
    pub log10_gls: Vec<f64>,
    pub ad: Vec<i32>,
}

pub fn subset_unused_alts_after_merged_genotyping(
    alt_alleles: &[String],
    gt_allele_indices: &[i32],
    log10_gls: &[f64],
    ad: &[i32],
) -> GatkResult<UnusedAltSubsetResult> {
    let n_alleles = 1 + alt_alleles.len();
    let keep = output_allele_keep_indices_from_assigned_gt(n_alleles, gt_allele_indices);
    if keep.len() == n_alleles {
        return Ok(UnusedAltSubsetResult {
            alt_alleles: alt_alleles.to_vec(),
            log10_gls: log10_gls.to_vec(),
            ad: ad.to_vec(),
        });
    }
    let new_gl = subset_log10_genotype_likelihoods(log10_gls, n_alleles, &keep)?;
    let new_ad = subset_ad(ad, &keep);
    let new_alts: Vec<String> = keep
        .iter()
        .skip(1)
        .map(|&i| alt_alleles[i - 1].clone())
        .collect();
    Ok(UnusedAltSubsetResult {
        alt_alleles: new_alts,
        log10_gls: new_gl,
        ad: new_ad,
    })
}

/// Subset log10 genotype likelihoods from `original_allele_count` to `keep_allele_indices`.
pub fn subset_log10_genotype_likelihoods(
    original_gl: &[f64],
    original_allele_count: usize,
    keep_allele_indices: &[usize],
) -> GatkResult<Vec<f64>> {
    let ploidy = 2;
    let expected = num_genotype_likelihoods(original_allele_count, ploidy);
    if original_gl.len() != expected {
        return Err(GatkError::argument(format!(
            "subset PL: expected {} log10 likelihoods for {original_allele_count} alleles, got {}",
            expected,
            original_gl.len()
        )));
    }
    let map = subsetted_pl_indices(ploidy, original_allele_count, keep_allele_indices);
    let picked: Vec<f64> = map.iter().map(|&idx| original_gl[idx]).collect();
    Ok(scale_log_space_array_for_numerical_stability(&picked))
}

/// Log10 GL to integer PL (GATK `GenotypeBuilder` / `FastGenotype` round-trip).
fn log10_gl_to_int_pl(gl: &[f64]) -> Vec<i32> {
    if gl.is_empty() {
        return Vec::new();
    }
    let max_v = gl.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let mut pls: Vec<i32> = gl
        .iter()
        .map(|g| {
            let shifted = (-10.0_f64) * (g - max_v);
            shifted.min(i32::MAX as f64).round() as i32
        })
        .collect();
    let min_pl = pls.iter().copied().min().unwrap_or(0);
    for pl in &mut pls {
        *pl -= min_pl;
    }
    pls
}

const SAC_STRANDS: usize = 2;

fn subset_sac(original_sac: &[i32], keep_allele_indices: &[usize]) -> GatkResult<Vec<i32>> {
    if original_sac.len() % SAC_STRANDS != 0 {
        return Err(GatkError::argument(format!(
            "SAC length {} not divisible by {SAC_STRANDS}",
            original_sac.len()
        )));
    }
    let mut out = Vec::with_capacity(keep_allele_indices.len() * SAC_STRANDS);
    for &allele_idx in keep_allele_indices {
        out.push(original_sac[SAC_STRANDS * allele_idx]);
        out.push(original_sac[SAC_STRANDS * allele_idx + 1]);
    }
    Ok(out)
}

/// One row of `g-subset-pl` / `g-subset-vc` fixture output.
/// # Invariants
/// `pl` / `ad` lengths match diploid genotype count / allele count **after** subsetting.
/// `allele_count_before` ≥ `allele_count_after`.
/// # Ownership
/// Owns PL/AD(/SAC) vectors for parity fixture comparison.
/// # Mutation
/// Immutable subset result.
/// # Biological assumptions
/// Allele subsetting remaps genotype likelihoods when dropping unused ALTs.
/// # Java equivalence
/// GATK allele-subsetting utilities (`AlleleSubsettingUtils` PL/AD/SAC remap).
#[derive(Debug, Clone, PartialEq)]
pub struct SubsetAllelesPlResult {
    pub allele_count_before: usize,
    pub allele_count_after: usize,
    pub pl: Vec<i32>,
    pub ad: Vec<i32>,
    pub gq: Option<i32>,
    pub sac: Option<Vec<i32>>,
}

/// Trim ACG to AC with het ref/C log10 PL fixture (GATK unit-test values).
pub fn subset_trim_acg_to_ac(log10_pl: &[f64], ad: &[i32]) -> GatkResult<SubsetAllelesPlResult> {
    let gl = subset_log10_genotype_likelihoods(log10_pl, 3, &[0, 1])?;
    let pl = log10_gl_to_int_pl(&gl);
    let ad_out = subset_ad(ad, &[0, 1]);
    Ok(SubsetAllelesPlResult {
        allele_count_before: 3,
        allele_count_after: 2,
        pl,
        ad: ad_out,
        gq: Some(200),
        sac: None,
    })
}

/// ACG to AG (GATK `AlleleSubsettingUtilsUnitTest` hetRefG3AllelesPL row).
pub fn subset_trim_acg_to_ag(log10_pl: &[f64], ad: &[i32]) -> GatkResult<SubsetAllelesPlResult> {
    let gl = subset_log10_genotype_likelihoods(log10_pl, 3, &[0, 2])?;
    let pl = log10_gl_to_int_pl(&gl);
    let ad_out = subset_ad(ad, &[0, 2]);
    Ok(SubsetAllelesPlResult {
        allele_count_before: 3,
        allele_count_after: 2,
        pl,
        ad: ad_out,
        gq: Some(200),
        sac: None,
    })
}

/// ACG to AC with PL + AD + SAC (GATK `AlleleSubsettingUtilsUnitTest` het AC row).
pub fn subset_het_ac_sac(
    log10_pl: &[f64],
    ad: &[i32],
    sac: &[i32],
) -> GatkResult<SubsetAllelesPlResult> {
    let gl = subset_log10_genotype_likelihoods(log10_pl, 3, &[0, 1])?;
    let pl = log10_gl_to_int_pl(&gl);
    let ad_out = subset_ad(ad, &[0, 1]);
    let sac_out = subset_sac(sac, &[0, 1])?;
    Ok(SubsetAllelesPlResult {
        allele_count_before: 3,
        allele_count_after: 2,
        pl,
        ad: ad_out,
        gq: Some(200),
        sac: Some(sac_out),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn trim_acg_to_ac_matches_gatk_unit_test() {
        let pl = [-20.0, 0.0, -20.0, -30.0, -40.0, -60.0];
        let ad = [14, 7, 1];
        let r = subset_trim_acg_to_ac(&pl, &ad).unwrap();
        assert_eq!(r.allele_count_after, 2);
        assert_eq!(r.pl, vec![200, 0, 200]);
        assert_eq!(r.ad, vec![14, 7]);
    }

    #[test]
    fn het_ac_sac_matches_gatk_unit_test() {
        let pl = [-20.0, 0.0, -20.0, -30.0, -40.0, -60.0];
        let ad = [14, 7, 1];
        let sac = [10, 9, 10, 9, 1, 1];
        let r = subset_het_ac_sac(&pl, &ad, &sac).unwrap();
        assert_eq!(r.pl, vec![200, 0, 200]);
        assert_eq!(r.ad, vec![14, 7]);
        assert_eq!(r.gq, Some(200));
        assert_eq!(r.sac, Some(vec![10, 9, 10, 9]));
    }

    #[test]
    fn unused_alt_subset_remaps_genotype_after_merged_genotyping() {
        // [TG, T, CG] GT=0/2 → drop unused deletion T, remap 0/2 → 0/1 on [TG, CG].
        let alts = vec!["T".to_string(), "CG".to_string()];
        let gls = vec![-29.8, -33.7, -162.0, 0.0, -105.8, -110.3];
        let ad = vec![28, 2, 10];
        let r = subset_unused_alts_after_merged_genotyping(&alts, &[0, 2], &gls, &ad).unwrap();
        assert_eq!(r.alt_alleles, vec!["CG".to_string()]);
        assert_eq!(r.ad, vec![28, 10]);
        assert_eq!(r.log10_gls.len(), 3);
        let pl = log10_gl_to_int_pl(&r.log10_gls);
        assert_eq!(pl, vec![298, 0, 1103]);
        assert_ne!(pl, vec![90, 30, 60, 30, 0, 60]);
        let best = pl
            .iter()
            .enumerate()
            .min_by_key(|(_, p)| *p)
            .map(|(i, _)| i)
            .unwrap();
        assert_eq!(best, 1, "USE_PLS_TO_ASSIGN on subsetted PLs is 0/1");
    }

    #[test]
    fn unused_alt_subset_keeps_second_alt_when_gt_uses_first() {
        let alts = vec!["T".to_string(), "CG".to_string()];
        let gls = vec![0.0, -1.0, -20.0, -30.0, -40.0, -60.0];
        let ad = vec![10, 8, 1];
        // GT=0/1 uses deletion T; CG is unused.
        let r = subset_unused_alts_after_merged_genotyping(&alts, &[0, 1], &gls, &ad).unwrap();
        assert_eq!(r.alt_alleles, vec!["T".to_string()]);
        assert_eq!(r.ad, vec![10, 8]);
        assert_eq!(r.log10_gls.len(), 3);
    }

    #[test]
    fn unused_alt_subset_keeps_both_when_gt_is_1_2() {
        let alts = vec!["T".to_string(), "CG".to_string()];
        let gls = vec![-9.0, -3.0, -6.0, -3.0, 0.0, -6.0];
        let ad = vec![5, 10, 12];
        let r = subset_unused_alts_after_merged_genotyping(&alts, &[1, 2], &gls, &ad).unwrap();
        assert_eq!(r.alt_alleles, vec!["T".to_string(), "CG".to_string()]);
        assert_eq!(r.ad, vec![5, 10, 12]);
        assert_eq!(r.log10_gls, gls);
    }

    #[test]
    fn unused_alt_subset_is_not_locus_specific() {
        let alts = vec!["A".to_string(), "CC".to_string()];
        let gls = vec![-20.0, -30.0, -40.0, 0.0, -10.0, -15.0];
        let ad = vec![12, 1, 9];
        let r = subset_unused_alts_after_merged_genotyping(&alts, &[0, 2], &gls, &ad).unwrap();
        assert_eq!(r.alt_alleles, vec!["CC".to_string()]);
        assert_eq!(r.ad, vec![12, 9]);
    }
}
