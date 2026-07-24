//! GATK `ExcessHet` exact test (I-D03 production parity).

const MIN_NEEDED_VALUE: f64 = 1.0e-16;
/// GATK `ExcessHet.PHRED_SCALED_MIN_P_VALUE` (-10 * log10(1e-16)).
const PHRED_SCALED_MIN_P_VALUE: f64 = 160.0;

/// Right-sided exact test p-value (GATK `ExcessHet.exactTest`).
pub fn exact_test(het_count: i32, ref_count: i32, hom_count: i32) -> f64 {
    assert!(het_count >= 0 && ref_count >= 0 && hom_count >= 0);
    let het_count = het_count as usize;
    let ref_count = ref_count as usize;
    let hom_count = hom_count as usize;

    let (obs_hom_r, obs_hom_c) = if ref_count < hom_count {
        (ref_count, hom_count)
    } else {
        (hom_count, ref_count)
    };

    let rare_copies = 2 * obs_hom_r + het_count;
    let n = het_count + obs_hom_c + obs_hom_r;
    if n == 0 {
        return 1.0;
    }

    let mut probs = vec![0.0_f64; rare_copies + 1];
    let mut mid = ((rare_copies as f64) * (2.0 * n as f64 - rare_copies as f64)
        / (2.0 * n as f64 - 1.0))
        .floor() as usize;
    if (mid % 2) != (rare_copies % 2) {
        mid += 1;
    }

    probs[mid] = 1.0;
    let mut mysum = 1.0;

    let mut curr_hets = mid;
    let mut curr_hom_r = (rare_copies - mid) / 2;
    let mut curr_hom_c = n - curr_hets - curr_hom_r;

    while curr_hets >= 2 {
        let potential_prob = probs[curr_hets] * curr_hets as f64 * (curr_hets - 1) as f64
            / (4.0 * (curr_hom_r + 1) as f64 * (curr_hom_c + 1) as f64);
        if potential_prob < MIN_NEEDED_VALUE {
            break;
        }
        probs[curr_hets - 2] = potential_prob;
        mysum += probs[curr_hets - 2];
        curr_hets -= 2;
        curr_hom_r += 1;
        curr_hom_c += 1;
    }

    curr_hets = mid;
    curr_hom_r = (rare_copies - mid) / 2;
    curr_hom_c = n - curr_hets - curr_hom_r;

    while curr_hets <= rare_copies.saturating_sub(2) {
        let potential_prob = probs[curr_hets] * 4.0 * curr_hom_r as f64 * curr_hom_c as f64
            / ((curr_hets + 2) as f64 * (curr_hets + 1) as f64);
        if potential_prob < MIN_NEEDED_VALUE {
            break;
        }
        probs[curr_hets + 2] = potential_prob;
        mysum += probs[curr_hets + 2];
        curr_hets += 2;
        curr_hom_r = curr_hom_r.saturating_sub(1);
        curr_hom_c = curr_hom_c.saturating_sub(1);
    }

    let right_pval = probs[het_count] / mysum;
    if het_count == rare_copies {
        return right_pval.clamp(0.0, 1.0);
    }
    let tail: f64 = probs[(het_count + 1)..].iter().sum();
    (right_pval + tail / mysum).clamp(0.0, 1.0)
}

/// Phred-scaled excess het (GATK `ExcessHet.calculateEH(GenotypeCounts, sampleCount)`).
pub fn excess_heterozygosity_phred(ref_count: u32, het_count: u32, hom_count: u32) -> f64 {
    let sample_count = ref_count + het_count + hom_count;
    if sample_count == 0 {
        return 0.0;
    }
    let pval = exact_test(het_count as i32, ref_count as i32, hom_count as i32);
    if pval < 1.0e-60 {
        return PHRED_SCALED_MIN_P_VALUE;
    }
    let phred = -10.0 * pval.log10();
    if phred == 0.0 {
        0.0
    } else {
        phred
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `i1-excess-het` / `balanced` row — matches Java `ExcessHet.calculateEH(0, 5, 5)`.
    #[test]
    fn balanced_ref_het_hom_matches_java() {
        let eh = excess_heterozygosity_phred(0, 5, 5);
        assert!(
            (eh - 2.838932).abs() < 1e-5,
            "expected Java golden 2.838932, got {eh}"
        );
    }

    /// All hom-ref samples → p≈1, phred≈0 (GATK `testPositiveZeroPhredScore` style counts).
    #[test]
    fn all_hom_ref_zero_phred() {
        let eh = excess_heterozygosity_phred(100, 0, 0);
        assert!(eh.abs() < 1e-6, "expected ~0 phred, got {eh}");
    }
}
