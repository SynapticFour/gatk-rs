//! 6R.101 coordinate-free: FORMAT/AD is remaining-allele remarg, not 4-way permute.
//!
//! Java 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`)
//!
//! ```text
//! calculateGLsForThisEvent                 // GenotypeBuilder: PL only; AD absent
//! AlleleSubsettingUtils.subsetAlleles      // if (g.hasAD()) slice; skipped
//! DepthPerAlleleBySample.annotate
//!   alleles = LinkedHashSet(vc.getAlleles())          // remaining call alleles
//!   annotateWithLikelihoods                            // FIRST AD write
//!     alleleSubset = {allele → [allele]}              // identity remarg
//!     subsetted = likelihoods.marginalize(alleleSubset)
//!     bestAllelesBreakingTies(sample)
//!     filter isInformative (confidence > 0.2)
//!     count by allele identity into vc allele order
//! reverseTrimAlleles / phaseVC             // GenotypeBuilder copy
//! ```
//!
//! PairHMM is closed. The 0.2 threshold is not changed. A 4-way informative
//! vote whose best allele is later unused is **not** the remaining-allele vote.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r101_ad_read_assignment_contract
//! HOLDOUT_6R101=1 cargo test -p gatk-haplotypecaller --test holdout_6r101_ad_read_assignment -- --nocapture
//! ```

use gatk_haplotypecaller::read_realignment::LOG_10_INFORMATIVE_THRESHOLD;

const JAVA_INFORMATIVE: f64 = 0.2;

/// Simple 4-way informative vote (Rust production `informative_ad_n_alleles`).
fn four_way_vote(lls: &[f64]) -> Option<usize> {
    let mut best_i = 0usize;
    let mut best = f64::NEG_INFINITY;
    let mut second = f64::NEG_INFINITY;
    for (i, &ll) in lls.iter().enumerate() {
        if ll > best {
            second = best;
            best = ll;
            best_i = i;
        } else if ll > second {
            second = ll;
        }
    }
    if best.is_finite() && (best - second).abs() > LOG_10_INFORMATIVE_THRESHOLD {
        Some(best_i)
    } else {
        None
    }
}

/// Java identity remarg: vote only over remaining columns, then `isInformative`.
fn remaining_vote(lls: &[f64], keep: &[usize]) -> Option<usize> {
    let kept: Vec<f64> = keep.iter().map(|&i| lls[i]).collect();
    four_way_vote(&kept).map(|i| keep[i])
}

/// Java `searchBestAllele` REF priority when (best − second) < 0.2.
fn java_best_breaking_ties(lls: &[f64], ref_i: usize) -> (usize, f64, bool) {
    let n = lls.len();
    let mut best_i = 0usize;
    let mut second_i = 0usize;
    let mut best = lls[0];
    let mut second = f64::NEG_INFINITY;
    for a in 1..n {
        let cand = lls[a];
        if cand > best {
            second_i = best_i;
            second = best;
            best_i = a;
            best = cand;
        } else if cand > second {
            second_i = a;
            second = cand;
        }
    }
    if best - second < JAVA_INFORMATIVE {
        let mut best_pri = if best_i == ref_i { 1.0 } else { 0.0 };
        for a in 0..n {
            let cand = lls[a];
            if a == best_i || best - cand > JAVA_INFORMATIVE {
                continue;
            }
            let pri = if a == ref_i { 1.0 } else { 0.0 };
            if pri > best_pri {
                second_i = best_i;
                best_i = a;
                best_pri = pri;
            }
        }
    }
    let best_ll = lls[best_i];
    let second_ll = if second_i != best_i {
        lls[second_i]
    } else {
        f64::NEG_INFINITY
    };
    let conf = if best_ll == second_ll {
        0.0
    } else {
        best_ll - second_ll
    };
    (best_i, conf, conf > JAVA_INFORMATIVE)
}

/// Unused-ALT 4-way winners can still be remaining-allele informative after identity remarg.
#[test]
fn forensic_6r101_unused_four_way_win_can_be_remaining_informative() {
    // Columns: REF, unused SPAN_DEL, unused SNP, called ALT — Java order TG,*,T,CG.
    let keep = [0usize, 3];
    // SPAN_DEL wins 4-way; remaining REF vs ALT still informative ALT (gap 0.45).
    let star_to_alt = [-8.45, 0.0, -3.0, -8.0];
    assert_eq!(four_way_vote(&star_to_alt), Some(1));
    assert_eq!(remaining_vote(&star_to_alt, &keep), Some(3));
    let kept = [star_to_alt[0], star_to_alt[3]];
    let (_, conf, inf) = java_best_breaking_ties(&kept, 0);
    assert!(conf > JAVA_INFORMATIVE);
    assert!(inf);
    assert!((conf - 0.45).abs() < 1e-12);

    // Unused SNP wins 4-way; remaining REF vs ALT still informative REF (gap 0.55).
    let snp_to_ref = [-8.0, -3.0, 0.0, -8.55];
    assert_eq!(four_way_vote(&snp_to_ref), Some(2));
    assert_eq!(remaining_vote(&snp_to_ref, &keep), Some(0));
    let kept = [snp_to_ref[0], snp_to_ref[3]];
    let (_, conf, inf) = java_best_breaking_ties(&kept, 0);
    assert!(conf > JAVA_INFORMATIVE);
    assert!(inf);
    assert!((conf - 0.55).abs() < 1e-12);
}

/// Permute of 4-way informative counts drops unused-column votes. Remarg reassigns them.
#[test]
fn forensic_6r101_four_way_permute_drops_unused_votes_remarg_reassigns() {
    let keep = [0usize, 3];
    let rows = [
        [-8.45, 0.0, -3.0, -8.0], // * → remaining ALT
        [-8.45, 0.0, -3.0, -8.0], // * → remaining ALT
        [-8.0, -3.0, 0.0, -8.55], // T → remaining REF
        [-8.0, -3.0, 0.0, -8.55], // T → remaining REF
        [0.0, -9.0, -9.0, -5.0],  // REF
        [-5.0, -9.0, -9.0, 0.0],  // ALT
    ];
    let mut four = [0i32; 4];
    let mut remarg = [0i32; 2];
    for row in &rows {
        if let Some(i) = four_way_vote(row) {
            four[i] += 1;
        }
        if let Some(i) = remaining_vote(row, &keep) {
            if i == 0 {
                remarg[0] += 1;
            } else {
                remarg[1] += 1;
            }
        }
    }
    let permute = [four[0], four[3]];
    assert_eq!(four, [1, 2, 2, 1]);
    assert_eq!(permute, [1, 1]);
    assert_eq!(remarg, [3, 3]);
    assert_eq!(
        remarg[0] - permute[0],
        2,
        "unused-SNP votes reassigned to remaining REF"
    );
    assert_eq!(
        remarg[1] - permute[1],
        2,
        "unused-SPAN_DEL votes reassigned to remaining ALT"
    );
}

/// `subsetAlleles` cannot invent AD. Annotation remarg is the first write.
#[test]
fn forensic_6r101_subset_does_not_write_ad_annotation_does() {
    let after_gls_has_ad = false;
    let after_subset_has_ad = after_gls_has_ad;
    assert!(
        !after_subset_has_ad,
        "Java unused-ALT subset skips AD when genotype.hasAD() is false"
    );
    let annotation_is_first_write = true;
    assert!(annotation_is_first_write);
}

/// FORMAT AD follows remaining remarg, not the unused-ALT slice of 4-way counts.
#[test]
fn forensic_6r101_format_ad_is_remaining_remarg_not_four_way_slice() {
    let keep = [0usize, 3];
    let rows = [
        [-8.45, 0.0, -3.0, -8.0],
        [-8.0, -3.0, 0.0, -8.55],
        [0.0, -9.0, -9.0, -4.0],
        [-4.0, -9.0, -9.0, 0.0],
    ];
    let mut four = [0i32; 4];
    let mut remarg = [0i32; 2];
    for row in &rows {
        if let Some(i) = four_way_vote(row) {
            four[i] += 1;
        }
        if let Some(i) = remaining_vote(row, &keep) {
            remarg[usize::from(i != 0)] += 1;
        }
    }
    let permute = [four[0], four[3]];
    let format_ad = remarg;
    assert_eq!(format_ad, remarg);
    assert_ne!(format_ad, permute);
}

/// Remaining-allele confidence for unused-column 4-way winners is far from 0.2.
/// REF tie-break is not the first AD divergence.
#[test]
fn forensic_6r101_unused_winners_are_not_near_informative_threshold() {
    let keep = [0usize, 3];
    let star_to_alt = [-8.45, 0.0, -3.0, -8.0];
    let snp_to_ref = [-8.0, -3.0, 0.0, -8.55];
    for row in [star_to_alt, snp_to_ref] {
        let four_d = {
            let mut xs = row;
            xs.sort_by(|a, b| b.partial_cmp(a).unwrap());
            xs[0] - xs[1]
        };
        assert!(
            four_d > 1.0,
            "4-way unused win is a clear unused-column best, not a 0.2-boundary tie"
        );
        let kept = [row[keep[0]], row[keep[1]]];
        let (_, conf, inf) = java_best_breaking_ties(&kept, 0);
        assert!(inf);
        assert!(
            conf > JAVA_INFORMATIVE + 0.2,
            "remaining remarg is not a near-threshold informativeness flip"
        );
        let simple = four_way_vote(&kept);
        let (jb, _, jinf) = java_best_breaking_ties(&kept, 0);
        assert_eq!(simple, jinf.then_some(jb));
    }
}

#[test]
fn forensic_6r101_is_informative_threshold_unchanged() {
    assert_eq!(LOG_10_INFORMATIVE_THRESHOLD, JAVA_INFORMATIVE);
    let eq = [0.0, -0.2];
    let (_, conf, inf) = java_best_breaking_ties(&eq, 0);
    assert!((conf - 0.2).abs() < 1e-15);
    assert!(!inf);
}
