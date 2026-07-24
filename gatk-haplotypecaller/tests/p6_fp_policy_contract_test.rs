use gatk_haplotypecaller::{pairhmm_fp_eq, PairHmmFpPolicy};

#[test]
fn step86_fp_policy_abs_and_rel_epsilon_contracts() {
    let p = PairHmmFpPolicy::default();
    assert!(pairhmm_fp_eq(-10.0, -10.0 + 1e-11, p));
    assert!(pairhmm_fp_eq(
        -1_000_000.0,
        -1_000_000.005,
        PairHmmFpPolicy {
            abs_epsilon: 1e-10,
            rel_epsilon: 1e-8
        }
    ));
    assert!(!pairhmm_fp_eq(-10.0, -9.8, p));
}
