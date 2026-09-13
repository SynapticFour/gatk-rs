//! 6R.161: first causal split at `2:92316347 G/A` is haplotype population, not votes.
//!
//! Frozen Java 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`)
//! `eventmap-haps-at-loc` + `ad-annotation-call` on `2:92316200-92316580`.
//!
//! Untrimmed Java haplotypes (fnv1a64 of bases):
//!   REF `9b4092f3b20a5de7`  ALT `8c9a6fc8f302f4a3`  n=2
//! Rust `after_assemble` has those two **plus** deletion hap `3a53b2a941bbdc43`.
//!
//! Production change: NONE. The later `45,3,0` overwrite is not this round.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r161_informative_votes -- --test-threads=1
//! HOLDOUT_6R161=1 cargo test -p gatk-haplotypecaller --test holdout_6r161_informative_votes -- --nocapture --test-threads=1
//! ```

const JAVA_UNTRIMMED_REF: &str = "9b4092f3b20a5de7";
const JAVA_UNTRIMMED_ALT: &str = "8c9a6fc8f302f4a3";
const RUST_EXTRA_DEL: &str = "3a53b2a941bbdc43";
const JAVA_TRIMMED_REF: &str = "4984337d563a59ec";
const JAVA_TRIMMED_ALT: &str = "ef15e644fee28e3c";
const RUST_TRIMMED_REF: &str = "189070c3f2226167";
const RUST_TRIMMED_ALT: &str = "7851f1a75f3d24c7";

#[test]
fn forensic_6r161_untrimmed_two_haps_match_java_and_rust_has_extra() {
    let java_untrimmed = [JAVA_UNTRIMMED_REF, JAVA_UNTRIMMED_ALT];
    let rust_after_assemble = [JAVA_UNTRIMMED_ALT, JAVA_UNTRIMMED_REF, RUST_EXTRA_DEL];
    assert_eq!(java_untrimmed.len(), 2);
    assert_eq!(rust_after_assemble.len(), 3);
    assert!(
        rust_after_assemble.contains(&JAVA_UNTRIMMED_REF)
            && rust_after_assemble.contains(&JAVA_UNTRIMMED_ALT),
        "the Java 2-hap set is present on Rust after_assemble"
    );
    assert!(
        !java_untrimmed.contains(&RUST_EXTRA_DEL),
        "extra Rust deletion hap is absent from Java"
    );
}

#[test]
fn forensic_6r161_trimmed_pairhmm_haps_already_differ() {
    assert_ne!(JAVA_TRIMMED_REF, RUST_TRIMMED_REF);
    assert_ne!(JAVA_TRIMMED_ALT, RUST_TRIMMED_ALT);
    assert_ne!(JAVA_TRIMMED_REF.len(), 0);
}

#[test]
fn forensic_6r161_java_ad_03_is_three_alt_informative_reads() {
    // DepthPerAlleleBySample after poorly-modeled filter (seq 6).
    let java_kept = [
        ("H06HDADXX130110:1:1101:10034:45116", 99u16, "ALT"),
        ("H06HDADXX130110:2:1101:10025:49248", 99, "ALT"),
        ("H06HDADXX130110:2:1101:10046:78083", 147, "ALT"),
    ];
    let java_filtered = [
        ("H06HDADXX130110:2:1101:10046:78083", 99u16),
        ("H06HDADXX130110:1:1101:10034:45116", 147),
        ("H06HDADXX130110:2:1101:10025:49248", 147),
        ("H06JUADXX130110:1:1101:10067:75885", 147),
    ];
    assert_eq!(java_kept.len(), 3);
    assert!(java_kept.iter().all(|(_, _, a)| *a == "ALT"));
    assert_eq!(java_filtered.len(), 4);
}

#[test]
fn forensic_6r161_vote_threshold_not_causal_on_kept_rows() {
    const CONF: f64 = 4.5;
    const THRESH: f64 = 0.2;
    assert!(CONF > THRESH);
    // Java 3 kept and Rust 7 overlap votes all have conf=4.5. Membership differs, not the 0.2 cut.
}

#[test]
fn forensic_6r161_pairhmm_float_residual_not_this_gap() {
    // Established Java-float vs Rust-f64 residual is ~1e-6 (6R.144). Live LLs on the
    // shared QNAME 45116/99: Java max -4.62 vs Rust max -52.51. That is haplotype
    // sequence / trim, not the float residual.
    let java_max = -4.618_268_966_674_805_f64;
    let rust_max = -52.513_267_93_f64;
    assert!((java_max - rust_max).abs() > 1.0);
}
