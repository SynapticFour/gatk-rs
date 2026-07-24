//! synthetic PairHMM rows → biallelic GL marginalization (one site).

use gatk_haplotypecaller::bio_ids::HaplotypeIndex;
use gatk_haplotypecaller::genotyping::ReadLikelihoodRow;
use gatk_haplotypecaller::hc_genotyping_engine::{
    biallelic_genotype_log10_likelihoods_gatk, marginalize_rows_to_biallelic_alleles,
};

#[test]
fn marginalize_pools_ref_and_alt_hap_columns() {
    let rows = vec![
        ReadLikelihoodRow {
            read_id: "r0".into(),
            haplotype_log10_likelihoods: vec![-0.1, -5.0, -0.5],
        },
        ReadLikelihoodRow {
            read_id: "r1".into(),
            haplotype_log10_likelihoods: vec![-0.2, -4.8, -0.4],
        },
    ];
    let marg = marginalize_rows_to_biallelic_alleles(
        &rows,
        &[HaplotypeIndex::new(0)],
        &[HaplotypeIndex::new(1), HaplotypeIndex::new(2)],
    );
    assert_eq!(marg.len(), 2);
    assert_eq!(marg[0].haplotype_log10_likelihoods.len(), 2);
    assert!(marg[0].haplotype_log10_likelihoods[0] > -0.25);
    assert!(marg[0].haplotype_log10_likelihoods[1] >= -0.5);
}

#[test]
fn biallelic_gl_fixture_heterozygous_with_read_supporting_alt() {
    let rows = vec![
        ReadLikelihoodRow {
            read_id: "r0".into(),
            haplotype_log10_likelihoods: vec![-2.0, -0.1],
        },
        ReadLikelihoodRow {
            read_id: "r1".into(),
            haplotype_log10_likelihoods: vec![-2.0, -0.1],
        },
    ];
    let gls = biallelic_genotype_log10_likelihoods_gatk(&rows, 0, 1);
    assert_eq!(gls.len(), 3);
    assert!(gls[1] > gls[0], "0/1 should beat 0/0: {:?}", gls);
}
