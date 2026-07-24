use gatk_haplotypecaller::{
    pairhmm_log10_likelihood, pairhmm_log10_likelihoods_vectorized, PairHmmInput, PairHmmParams,
};

#[test]
fn scalar_and_vectorized_pairhmm_are_equivalent() {
    let params = PairHmmParams::default();
    let read = "ACGTACGTAC";
    let quals = vec![32; read.len()];
    let haps = vec![
        "ACGTACGTAC".to_string(),
        "ACGTTCGTAC".to_string(),
        "ACGTACGTTC".to_string(),
        "ACGTACGGAC".to_string(),
    ];
    let vec_out = pairhmm_log10_likelihoods_vectorized(read, &quals, 60, &haps, &params)
        .expect("vectorized output");
    assert_eq!(vec_out.len(), haps.len());

    for (idx, hap) in haps.iter().enumerate() {
        let scalar = pairhmm_log10_likelihood(
            &PairHmmInput {
                read_bases: read.to_string(),
                read_base_quals: quals.clone(),
                read_mapping_quality: 60,
                haplotype_bases: hap.clone(),
            },
            &params,
        )
        .expect("scalar output");
        assert!(
            (scalar - vec_out[idx]).abs() <= 1e-12,
            "idx={idx} scalar={scalar} vectorized={}",
            vec_out[idx]
        );
    }
}
