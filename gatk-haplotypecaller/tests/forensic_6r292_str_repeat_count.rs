//! 6R.292: why `TandemRepeat.getNumTandemRepeatUnits` returns 9 for A>AT at
//! 20:29455644 while `longest_str_len_at_variant` returns `None`.
//!
//! Stops at slice construction. The two functions receive the same reference
//! window and the same event, then inspect different bytes.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r292_str_repeat_count -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{ReferenceWindowCache, SequenceDictionary};
use gatk_haplotypecaller::{longest_str_len_at_variant, ReferenceContext};
use std::collections::HashMap;
use std::path::{Path, PathBuf};

const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const JAVA_TSV: &str = "gatk-haplotypecaller/tests/6r292_java_str.tsv";
const WINDOW_START: u64 = 29_455_460;
const WINDOW_END: u64 = 29_455_844;
const EVENT: u64 = 29_455_644;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R292\t{key}\t{}", value.as_ref());
}

fn java_rows() -> HashMap<String, String> {
    let path = repo_root().join(JAVA_TSV);
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("read {path:?}: {e}"));
    text.lines()
        .filter(|line| line.starts_with("6R292\t"))
        .map(|line| {
            let mut parts = line.splitn(3, '\t');
            let _ = parts.next();
            let key = parts.next().expect("key").to_string();
            let value = parts.next().expect("value").to_string();
            (key, value)
        })
        .collect()
}

#[test]
fn forensic_6r292_str_repeat_count() {
    let java = java_rows();
    assert_eq!(java["window"], "20:29455460-29455844");
    assert_eq!(java["ref_allele"], "A");
    assert_eq!(java["alt_allele"], "AT");
    assert_eq!(java["anchor_base"], "A");
    assert_eq!(java["java_context_start"], "29455645");
    assert_eq!(java["ref_allele_after_anchor"], "");
    assert_eq!(java["alt_after_anchor"], "T");
    assert_eq!(java["repeat_unit"], "T");
    assert_eq!(java["repeat_unit_len"], "1");
    assert_eq!(java["repeat_counts"], "[8, 9]");
    assert_eq!(java["repeat_count_max"], "9");
    assert_eq!(java["str_addition"], "9");
    assert!(
        java["java_remaining_prefix"].starts_with("TTTTTTTT"),
        "java context {}",
        java["java_remaining_prefix"]
    );

    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    if !ref_fasta.is_file() {
        eprintln!("skip: missing chr20_tiny BAM/REF");
        return;
    }
    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let mut cache = ReferenceWindowCache::new(&ref_fasta, 4);
    let ctx = ReferenceContext::from_interval(&dict, &mut cache, "20", WINDOW_START, WINDOW_END)
        .expect("ref");
    assert_eq!(ctx.start, WINDOW_START);
    let offset = (EVENT - ctx.start) as usize;
    let bases: &[u8] = &ctx.bases;
    let rust_slice = &bases[offset..=offset];
    let after = &bases[offset + 1..offset + 9];
    assert_eq!(rust_slice, b"A");
    assert_eq!(after, b"TTTTTTTT");
    assert_eq!(
        longest_str_len_at_variant(&ctx, EVENT, EVENT),
        None,
        "one-base event span returns before any repeat search"
    );

    kv("classification", "STR_REPEAT_SLICE_DIVERGENCE");
    kv("reference_window", "20:29455460-29455844");
    kv("event", "29455644 A>AT");
    kv(
        "ref_at_event_and_after",
        &ctx.bases_ascii()[offset..offset + 16],
    );
    kv("java_slice", &java["java_remaining_prefix"]);
    kv("java_slice_start", "29455645");
    kv("java_alt_suffix", "T");
    kv("java_repeat_unit", "T");
    kv("java_repeat_counts", "[8, 9]");
    kv("java_repeat_count", "9");
    kv("rust_slice", "A");
    kv("rust_slice_span", "29455644-29455644");
    kv("rust_repeat_count", "None");
    kv(
        "first_divergence",
        "slice construction: Java starts at start+1 and uses the alt suffix; Rust uses [start,end] and returns None because that slice has length 1",
    );
    kv("production_change", "NONE");
}
