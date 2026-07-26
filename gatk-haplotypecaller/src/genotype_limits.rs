//! GATK `GenotypeLikelihoodCalculators.computeMaxAcceptableAlleleCount`.

use gatk_common::{GatkError, GatkResult};

fn log10_factorial(n: u32) -> f64 {
    if n <= 1 {
        return 0.0;
    }
    (2..=n).map(|i| (i as f64).log10()).sum()
}

fn log10_binomial(n: u32, k: u32) -> f64 {
    log10_factorial(n) - log10_factorial(k) - log10_factorial(n - k)
}

/// GATK `computeMaxAcceptableAlleleCount(ploidy, maxGenotypeCount)`.
pub fn compute_max_acceptable_allele_count(
    ploidy: u32,
    max_genotype_count: u32,
) -> GatkResult<u32> {
    if ploidy == 0 {
        return Err(GatkError::argument("ploidy must be >= 1"));
    }
    if max_genotype_count == 0 {
        return Err(GatkError::argument("max_genotype_count must be >= 1"));
    }
    if ploidy == 1 {
        return Ok(max_genotype_count);
    }
    let log10_max = (max_genotype_count as f64).log10();
    let x = 10_f64.powf((log10_factorial(ploidy) + log10_max) / (ploidy as f64));
    let lower = ((x.floor() as i32) - (ploidy as i32) - 1).max(2) as u32;
    let upper = x.ceil().max(2.0) as u32;
    for a in (lower..=upper).rev() {
        let log10_gt = log10_binomial(ploidy + a - 1, a - 1);
        if log10_max >= log10_gt {
            return Ok(a);
        }
    }
    Ok(2)
}

/// Haplotype score precedence for allele trimming (GATK `whichAllelesToKeepBasedonHapScores` lite).
/// # Invariants
/// `best_score` ≥ `second_best_score` for the allele’s haplotype score list.
/// `allele_index` indexes the allele list being trimmed.
/// # Ownership
/// Owned score triple used while ranking alleles to keep.
/// # Mutation
/// Ephemeral ranking record; discarded after allele selection.
/// # Biological assumptions
/// Alleles with stronger supporting haplotype scores are preferred when capping allele count.
/// # Java equivalence
/// GATK `HaplotypeCallerGenotypingEngine.whichAllelesToKeepBasedonHapScores` ranking cell.
#[derive(Debug, Clone, PartialEq)]
pub struct AlleleHaplotypeScore {
    pub allele_index: usize,
    pub best_score: f64,
    pub second_best_score: f64,
}

/// Best and second-best scores without allocating a sorted clone (ascending-sort last two).
#[inline]
fn top_two_scores(scores: &[f64]) -> (f64, f64) {
    let mut best = f64::NEG_INFINITY;
    let mut second = f64::NEG_INFINITY;
    for &s in scores {
        if s > best {
            second = best;
            best = s;
        } else if s > second {
            second = s;
        }
    }
    (best, second)
}

/// GATK `HaplotypeCallerGenotypingEngine.whichAllelesToKeepBasedonHapScores`.
#[allow(dead_code)] // thin wrapper over `_with_ref`; kept for GATK-shaped call sites
pub fn which_alleles_to_keep_by_haplotype_scores<S: AsRef<[f64]>>(
    scores_per_allele: &[S],
    desired_allele_count: usize,
) -> Vec<usize> {
    which_alleles_to_keep_by_haplotype_scores_with_ref(
        scores_per_allele,
        desired_allele_count,
        None,
    )
}

/// Same as [`which_alleles_to_keep_by_haplotype_scores`] with optional per-allele reference flags (GATK ref-priority tie-break).
pub fn which_alleles_to_keep_by_haplotype_scores_with_ref<S: AsRef<[f64]>>(
    scores_per_allele: &[S],
    desired_allele_count: usize,
    is_reference: Option<&[bool]>,
) -> Vec<usize> {
    if scores_per_allele.len() <= desired_allele_count {
        return (0..scores_per_allele.len()).collect();
    }
    let mut ranked: Vec<AlleleHaplotypeScore> = scores_per_allele
        .iter()
        .enumerate()
        .map(|(i, scores)| {
            // Top-2 without sorting a cloned score vector (same as sort-then-take-last-two).
            let (best, second) = top_two_scores(scores.as_ref());
            AlleleHaplotypeScore {
                allele_index: i,
                best_score: best,
                second_best_score: second,
            }
        })
        .collect();
    ranked.sort_by(|a, b| compare_haplotype_allele_scores(a, b, is_reference));
    let keep = desired_allele_count.min(ranked.len());
    let mut out: Vec<usize> = ranked.iter().take(keep).map(|a| a.allele_index).collect();
    out.sort_unstable();
    out
}

fn compare_haplotype_allele_scores(
    a: &AlleleHaplotypeScore,
    b: &AlleleHaplotypeScore,
    is_reference: Option<&[bool]>,
) -> std::cmp::Ordering {
    if let Some(flags) = is_reference {
        let a_ref = flags.get(a.allele_index).copied().unwrap_or(false);
        let b_ref = flags.get(b.allele_index).copied().unwrap_or(false);
        if a_ref && !b_ref {
            return std::cmp::Ordering::Less;
        }
        if !a_ref && b_ref {
            return std::cmp::Ordering::Greater;
        }
    }
    b.best_score
        .total_cmp(&a.best_score)
        .then(b.second_best_score.total_cmp(&a.second_best_score))
        .then(a.allele_index.cmp(&b.allele_index))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn trim_three_fixture_matches_gatk() {
        let scores = vec![vec![0.0], vec![-0.5], vec![-1.0], vec![0.1]];
        let is_ref = [true, false, false, false];
        let kept = which_alleles_to_keep_by_haplotype_scores_with_ref(&scores, 2, Some(&is_ref));
        assert_eq!(kept, vec![0, 3]);
    }
}
