use gatk_haplotypecaller::{pairhmm_log10_likelihood, PairHmmInput, PairHmmParams};

#[test]
fn step87_rejects_invalid_probability_params() {
    let input = PairHmmInput {
        read_bases: "ACGT".to_string(),
        read_base_quals: vec![30; 4],
        read_mapping_quality: 60,
        haplotype_bases: "ACGT".to_string(),
    };
    let bad = PairHmmParams {
        gap_open_prob: 1.2,
        gap_extend_prob: 0.1,
        insertion_emission_prob: 0.25,
    };
    assert!(pairhmm_log10_likelihood(&input, &bad).is_err());
}

#[test]
fn step87_rejects_malformed_read_quality_lengths() {
    let input = PairHmmInput {
        read_bases: "ACGT".to_string(),
        read_base_quals: vec![30; 3],
        read_mapping_quality: 60,
        haplotype_bases: "ACGT".to_string(),
    };
    let ok = PairHmmParams::default();
    assert!(pairhmm_log10_likelihood(&input, &ok).is_err());
}

#[test]
fn step87_handles_extreme_quality_edges_without_nan() {
    let input = PairHmmInput {
        read_bases: "NNNNNNNN".to_string(),
        read_base_quals: vec![0, 0, 0, 0, 60, 60, 60, 60],
        read_mapping_quality: 0,
        haplotype_bases: "ACGTACGT".to_string(),
    };
    let ll = pairhmm_log10_likelihood(&input, &PairHmmParams::default()).expect("likelihood");
    assert!(ll.is_finite());
    assert!(!ll.is_nan());
}
