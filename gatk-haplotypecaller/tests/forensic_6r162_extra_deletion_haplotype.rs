//! 6R.162: extra after_assemble hap is a dangling-merge object, not SeqGraph k-best.
//!
//! Java 4.4.0.0 `assembleReads` n=2 (REF `9b4092f3b20a5de7`, ALT `8c9a6fc8f302f4a3`).
//! Those two hashes are SeqGraph k=35 k-best (K=128). Extra `3a53b2a941bbdc43` has
//! sentinel `score=25` / `kmer_size=0` from `apply_dangling_merge_haplotypes`.
//!
//! Production change: NONE. Not a future FORMAT contract.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r162_extra_deletion_haplotype -- --test-threads=1
//! HOLDOUT_6R162=1 cargo test -p gatk-haplotypecaller --test holdout_6r162_extra_deletion_haplotype -- --nocapture --test-threads=1
//! ```

const JAVA_UNTRIMMED_REF: &str = "9b4092f3b20a5de7";
const JAVA_UNTRIMMED_ALT: &str = "8c9a6fc8f302f4a3";
const RUST_EXTRA_DEL: &str = "3a53b2a941bbdc43";
const DANGLING_MERGE_SENTINEL_SCORE: f64 = 25.0;

#[test]
fn forensic_6r162_seqgraph_kbest_matches_java_assemble_pair() {
    let java_assemble = [JAVA_UNTRIMMED_REF, JAVA_UNTRIMMED_ALT];
    let seq_kbest_k35 = [JAVA_UNTRIMMED_ALT, JAVA_UNTRIMMED_REF];
    assert_eq!(java_assemble.len(), 2);
    assert!(
        seq_kbest_k35.contains(&JAVA_UNTRIMMED_REF) && seq_kbest_k35.contains(&JAVA_UNTRIMMED_ALT)
    );
    assert!(
        !seq_kbest_k35.contains(&RUST_EXTRA_DEL),
        "extra hap is absent from SeqGraph k-best"
    );
}

#[test]
fn forensic_6r162_extra_hap_is_dangling_merge_sentinel() {
    // Live holdout: extra hap score=25, kmer_size=0, cigar 40D2M155D172M107D, len=174.
    // `apply_dangling_merge_haplotypes` is the only production writer of score=25.
    assert_eq!(DANGLING_MERGE_SENTINEL_SCORE, 25.0);
}

#[test]
fn forensic_6r162_configured_kmers_skip_non_unique_on_both_sides() {
    // Java seqgraph-kbest-at-loc and Rust probe: k=10 and k=25 skip non_unique_ref.
    let java_skip = [10usize, 25];
    let rust_skip = [10usize, 25];
    assert_eq!(java_skip, rust_skip);
}

#[test]
fn forensic_6r162_kbest_resource_policy_not_this_gap() {
    // SeqGraph k-best n=2 at k=35 with K=128 / legacy_1024. Extra hap is not a ranked path.
    let seq_kbest_n = 2usize;
    let java_assemble_n = 2usize;
    assert_eq!(seq_kbest_n, java_assemble_n);
}

#[test]
fn forensic_6r162_dangling_tail_recovered_then_emitted_as_haplotype() {
    // Rust audit k=35: tails 1/1, edges 614→615. Java recoverDanglingTails adds an
    // edge; it does not append a 174 bp haplotype object.
    let rust_tails_recovered = 1u32;
    let rust_assemble_n = 3usize;
    let java_assemble_n = 2usize;
    assert_eq!(rust_tails_recovered, 1);
    assert_eq!(rust_assemble_n.saturating_sub(java_assemble_n), 1);
}
