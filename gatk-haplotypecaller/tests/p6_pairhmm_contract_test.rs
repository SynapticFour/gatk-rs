use gatk_haplotypecaller::{pairhmm_log10_likelihood, PairHmmInput, PairHmmParams};

fn mk_input(read: &str, quals: &[u8], mapq: u8, hap: &str) -> PairHmmInput {
    PairHmmInput {
        read_bases: read.to_string(),
        read_base_quals: quals.to_vec(),
        read_mapping_quality: mapq,
        haplotype_bases: hap.to_string(),
    }
}

#[test]
fn contract_step77_state_machine_prefers_match_path() {
    let params = PairHmmParams::default();
    let exact = mk_input("ACTGACTG", &[30; 8], 60, "ACTGACTG");
    let off = mk_input("ACTGACTG", &[30; 8], 60, "ACTAACTG");
    let ll_exact = pairhmm_log10_likelihood(&exact, &params).expect("exact");
    let ll_off = pairhmm_log10_likelihood(&off, &params).expect("off");
    assert!(ll_exact > ll_off);
}

#[test]
fn contract_step78_base_and_map_quality_cap_affects_likelihood() {
    let params = PairHmmParams::default();
    let hi = mk_input("AAAA", &[35, 35, 35, 35], 60, "AAAA");
    let lo = mk_input("AAAA", &[35, 35, 35, 35], 5, "AAAA");
    let ll_hi = pairhmm_log10_likelihood(&hi, &params).expect("high");
    let ll_lo = pairhmm_log10_likelihood(&lo, &params).expect("low");
    assert!(ll_hi > ll_lo);
}

#[test]
fn contract_step79_long_indel_tail_stays_finite() {
    let params = PairHmmParams::default();
    let read = "ACGT".repeat(64);
    let hap = format!("{}{}", "ACGT".repeat(48), "T".repeat(64));
    let quals = vec![20_u8; read.len()];
    let input = mk_input(&read, &quals, 40, &hap);
    let ll = pairhmm_log10_likelihood(&input, &params).expect("likelihood");
    assert!(
        ll.is_finite(),
        "likelihood must remain finite for long tails"
    );
}
