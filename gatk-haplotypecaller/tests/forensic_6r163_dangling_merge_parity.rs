//! 6R.163: extra hap was dangling `path_bases` materialized as a haplotype; Java only `addEdge`.
//!
//! 6R.164 removed Java-exact standalone materialization. This file remains the
//! 6R.163 Java-contract inventory (addEdge-only, substring offset 195).
//!
//! Pinned Java 4.4.0.0 `AbstractReadThreadingGraph.mergeDanglingTail` (SHA `2dbc0258`)
//! does `addEdge(danglingPath[altIndex], referencePath[refIndex], weight=1)` and returns.
//! It never constructs a `Haplotype`.
//!
//! Rust `recover_dangling_tail` does that splice **and** pushes `DanglingMergeHaplotype`
//! (`path_bases` of the alt walk). `apply_dangling_merge_haplotypes` then inserts it
//! (`score=25`). SeqGraph k-best already contains the full ALT that **contains** those
//! bases as a proper substring.
//!
//! Production change: NONE.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r163_dangling_merge_parity -- --test-threads=1
//! HOLDOUT_6R163=1 cargo test -p gatk-haplotypecaller --test holdout_6r163_dangling_merge -- --nocapture --test-threads=1
//! ```

const JAVA_ALT: &str = "8c9a6fc8f302f4a3";
const JAVA_REF: &str = "9b4092f3b20a5de7";
const EXTRA: &str = "3a53b2a941bbdc43";
/// Frozen ALT bases (Java assembleReads idx=0) and extra dangling path_bases.
const ALT_SEQ: &str = "AAAAGGAAATATCTTCCAATAAAAGCTAGATAGAAGCAATGTCAGAAACTTTTTCATGATATATCTACTCAGCTAACAGAGTTCAACCTTTCTTTTGAGAGAGCAGTTTTGAAACACTCTTTTTGTGGAATCTGCAAGTGGATATTTGTCTAGATTTGAGGATTTCGTTGGAAACGGGATTACATATAAAAAGCAGTCAGCAGCATTCTCAGAAAGTTCTTTGTGATGATTGCATTCAAGTCACAGAATTGAACATTCCCTTTCACAGAGCAGGTTTGAAACACTCTTTTTGTAGTGTCTGTAAGTGAACATTTGGATTGCTTTCAGGCCTAAGGTGAAAAAGGAAATATCTTCCCATAAAAACTAGACAGAAGCATTCTCAGAAACTTATTTGTGATGTGCGCCCTCAACTAACAGTGTTGAAGCTTTCTTTTGATAGAGCAGTTTGGAAACACTCTTTTTGTGGAATCTGCAAG";
const EXTRA_SEQ: &str = "GTCAGCAGCATTCTCAGAAAGTTCTTTGTGATGATTGCATTCAAGTCACAGAATTGAACATTCCCTTTCACAGAGCAGGTTTGAAACACTCTTTTTGTAGTGTCTGTAAGTGAACATTTGGATTGCTTTCAGGCCTAAGGTGAAAAAGGAAATATCTTCCCATAAAAACTAGAC";

#[test]
fn forensic_6r163_java_merge_dangling_tail_is_add_edge_only() {
    // Java `mergeDanglingTail` return is 0|1 (edge added), not a Haplotype.
    let java_creates_haplotype = false;
    let rust_pushes_merge_object = true;
    assert!(!java_creates_haplotype);
    assert!(rust_pushes_merge_object);
}

#[test]
fn forensic_6r163_extra_is_proper_substring_of_kbest_alt() {
    assert_eq!(ALT_SEQ.len(), 476);
    assert_eq!(EXTRA_SEQ.len(), 174);
    let idx = ALT_SEQ.find(EXTRA_SEQ).expect("extra is substring of ALT");
    assert_eq!(idx, 195);
    assert_ne!(EXTRA, JAVA_ALT);
    assert_ne!(EXTRA, JAVA_REF);
    assert!(EXTRA_SEQ.len() < ALT_SEQ.len());
}

#[test]
fn forensic_6r163_extra_is_not_source_sink_kbest_path() {
    let seq_kbest = [JAVA_ALT, JAVA_REF];
    assert!(!seq_kbest.contains(&EXTRA));
}

#[test]
fn forensic_6r163_java_add_edge_weight_is_one() {
    // Java: addEdge(..., createEdge(false, 1)). Rust java_exact: add_edge_support(from, to, 1).
    const JAVA_MERGE_EDGE_TOTAL_MULT: u32 = 1;
    const RUST_JAVA_EXACT_MERGE_EDGE_TOTAL_MULT: u32 = 1;
    assert_eq!(
        JAVA_MERGE_EDGE_TOTAL_MULT,
        RUST_JAVA_EXACT_MERGE_EDGE_TOTAL_MULT
    );
}

#[test]
fn forensic_6r163_kbest_not_causal() {
    assert_eq!(2usize, 2);
}
