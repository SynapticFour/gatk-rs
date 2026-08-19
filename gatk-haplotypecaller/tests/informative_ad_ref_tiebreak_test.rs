//! Property: informative AD near-ties are uninformative (Java DepthPerAlleleBySample).
//! GATK `BestAllele.isInformative` requires confidence > LOG_10_INFORMATIVE_THRESHOLD.
//! Tie-break toward REF chooses the allele label only; uninformative votes still do not count.

use gatk_haplotypecaller::genotyping::ReadLikelihoodRow;
use gatk_haplotypecaller::hc_genotyping_engine::{
    biallelic_allele_depths_from_rows_min_conf, InformativeAd,
};

#[test]
fn informative_ad_near_tie_is_uninformative() {
    // Within threshold: ALT slightly ahead — Java drops the vote (does not count as REF).
    let rows = vec![ReadLikelihoodRow {
        read_index: 0,
        read_id: "r0".into(),
        haplotype_log10_likelihoods: vec![-1.0, -1.01],
    }];
    let ad = biallelic_allele_depths_from_rows_min_conf(&rows, 0, 1, Some(0.1));
    assert_eq!(ad, vec![0, 0]);
    assert_eq!(
        InformativeAd::from_marginalized_rows(&rows, 0, 1, Some(0.1)),
        InformativeAd {
            ref_depth: 0,
            alt_depth: 0
        }
    );
}

#[test]
fn informative_ad_clear_alt_winner_counts_alt() {
    let rows = vec![ReadLikelihoodRow {
        read_index: 0,
        read_id: "r0".into(),
        haplotype_log10_likelihoods: vec![-3.0, -0.5],
    }];
    let ad = biallelic_allele_depths_from_rows_min_conf(&rows, 0, 1, Some(0.1));
    assert_eq!(ad, vec![0, 1]);
}

#[test]
fn informative_ad_clear_ref_winner_counts_ref() {
    let rows = vec![ReadLikelihoodRow {
        read_index: 0,
        read_id: "r0".into(),
        haplotype_log10_likelihoods: vec![-0.5, -3.0],
    }];
    let ad = biallelic_allele_depths_from_rows_min_conf(&rows, 0, 1, Some(0.1));
    assert_eq!(ad, vec![1, 0]);
}
