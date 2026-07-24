use gatk_haplotypecaller::{
    pairhmm_log10_likelihood, pairhmm_log10_likelihoods_vectorized, PairHmmInput, PairHmmParams,
};

#[test]
fn step83_long_indel_and_low_quality_tail_boundaries() {
    let params = PairHmmParams::default();
    let read = format!("{}{}", "ACGT".repeat(40), "T".repeat(20));
    let hap = "ACGT".repeat(30);
    let mut quals = vec![30; 160];
    quals.extend(vec![2; 20]); // low-quality tail
    let input = PairHmmInput {
        read_bases: read,
        read_base_quals: quals,
        read_mapping_quality: 45,
        haplotype_bases: hap,
    };
    let ll = pairhmm_log10_likelihood(&input, &params).expect("likelihood");
    assert!(ll.is_finite());
}

#[test]
fn step84_clipping_adapter_n_base_artifact_corpus() {
    let params = PairHmmParams::default();
    let read = "NNNNACGTACGTNNNN";
    let quals = vec![5, 5, 5, 5, 30, 30, 30, 30, 30, 30, 30, 30, 5, 5, 5, 5];
    let haps = vec![
        "ACGTACGTACGTACGT".to_string(),
        "NNNNACGTACGTNNNN".to_string(),
        "ACGTTCGTACGTACGT".to_string(),
    ];
    let out = pairhmm_log10_likelihoods_vectorized(read, &quals, 20, &haps, &params).expect("vec");
    assert_eq!(out.len(), haps.len());
    assert!(out.iter().all(|v| v.is_finite()));
}
