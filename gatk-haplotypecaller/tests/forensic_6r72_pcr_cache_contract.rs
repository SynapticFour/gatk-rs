//! 6R.72 coordinate-free: Java 4.4 CONSERVATIVE PCR cache vs production Rust.
//!
//! SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`
//! `PairHMMLikelihoodCalculationEngine.getErrorModelAdjustedQual` /
//! `initializePCRErrorModel`.
//!
//! Repeat finder is frozen (6R.71). 6R.73: production GOP source is BI/BD.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r72_pcr_cache_contract
//! ```

use gatk_haplotypecaller::pairhmm_log10::GATK_PARITY_DEFAULT_INS_QUAL;
use gatk_haplotypecaller::pcr_error_model::{
    apply_pcr_error_model, error_model_adjusted_qual, find_tandem_repeat_units,
    tandem_repeat_units, PcrErrorModel, MAX_REPEAT_LENGTH,
};

/// GATK 4.4 `MathUtils.fastRound`: `(int)(d+0.5)` if `d>0`, else `(int)(d-0.5)`.
fn java_fast_round(d: f64) -> i32 {
    if d > 0.0 {
        (d + 0.5) as i32
    } else {
        (d - 0.5) as i32
    }
}

fn java_raw(r: usize) -> f64 {
    40.0 - (r as f64 / (3.0 * std::f64::consts::PI)).exp() + 1.0
}

fn java_44_cache(r: usize) -> u8 {
    java_fast_round(java_raw(r)).max(10) as u8
}

fn rust_cache(r: usize) -> u8 {
    error_model_adjusted_qual(r, PcrErrorModel::Conservative.rate_factor().unwrap())
}

/// Pinned Java CONSERVATIVE `pcrIndelErrorModelCache` for r=0..=40
/// (`INITIAL_QSCORE=40`, rate 3.0, min 10, `fastRound`).
const JAVA_CACHE_0_TO_40: [u8; 41] = [
    40, 40, 40, 40, 39, 39, 39, 39, 39, 38, 38, 38, 37, 37, 37, 36, 36, 35, 34, 33, 33, 32, 31, 30,
    28, 27, 25, 23, 21, 19, 17, 14, 11, 10, 10, 10, 10, 10, 10, 10, 10,
];

#[test]
fn conservative_rate_and_min_match_java_4_4() {
    assert_eq!(PcrErrorModel::Conservative.rate_factor(), Some(3.0));
    assert_eq!(PcrErrorModel::Aggressive.rate_factor(), Some(2.0));
    assert_eq!(MAX_REPEAT_LENGTH, 20);
}

#[test]
fn java_cache_table_0_to_40() {
    for r in 0..=40 {
        let raw = java_raw(r);
        let rounded = java_fast_round(raw);
        let clamped = rounded.max(10) as u8;
        assert_eq!(clamped, JAVA_CACHE_0_TO_40[r], "Java cache[{r}]");
        assert_eq!(java_44_cache(r), JAVA_CACHE_0_TO_40[r]);
        eprintln!(
            "r={r} raw={raw:.12} rounded={rounded} clamp={clamped} rust={}",
            rust_cache(r)
        );
    }
    assert_eq!(java_44_cache(0), 40);
    assert_eq!(java_44_cache(1), 40);
    assert_eq!(java_44_cache(2), 40);
    assert_eq!(java_44_cache(3), 40);
    assert_eq!(java_44_cache(4), 39);
    assert_eq!(java_44_cache(5), 39);
    assert_eq!(java_44_cache(10), 38);
    assert_eq!(java_44_cache(20), 33);
    assert_eq!(java_44_cache(32), 11);
    assert_eq!(java_44_cache(33), 10);
    assert_eq!(java_44_cache(34), 10);
}

#[test]
fn production_cache_matches_java_on_finder_and_clamp_indices() {
    for &r in &[0, 1, 2, 3, 4, 5, 10, 20, 33, 34] {
        assert_eq!(rust_cache(r), java_44_cache(r), "cache[{r}]");
    }
    // Live cache is 0..=20; 33/34 are the formula clamp, not a finder index.
    for r in 0..=MAX_REPEAT_LENGTH {
        assert_eq!(rust_cache(r), JAVA_CACHE_0_TO_40[r]);
    }
}

#[test]
fn truncating_the_raw_float_would_change_cache_1() {
    // r=1 raw ≈ 39.888. Truncation → 39; fastRound → 40.
    let raw = java_raw(1);
    assert!(raw > 39.0 && raw < 40.0, "{raw}");
    let truncated = raw as i32;
    assert_eq!(truncated, 39);
    assert_eq!(java_fast_round(raw), 40);
    assert_eq!(rust_cache(1), 40);
}

#[test]
fn min_clamp_transitions_at_repeat_33() {
    assert!(java_fast_round(java_raw(32)) > 10);
    assert!(java_fast_round(java_raw(33)) < 10);
    assert_eq!(java_44_cache(32), 11);
    assert_eq!(java_44_cache(33), 10);
    assert_eq!(rust_cache(33), 10);
    assert_eq!(rust_cache(34), 10);
    // Finder clips unit count to 20, so applyPCR never indexes 33.
    let long_a = vec![b'A'; 40];
    assert_eq!(tandem_repeat_units(&long_a, 0), 20);
}

#[test]
fn cache_index_equals_repeat_length_and_zero_is_not_selected() {
    assert_eq!(find_tandem_repeat_units(b"T", 0), (b"T".to_vec(), 1));
    assert_eq!(rust_cache(1), 40);
    let seq = b"TAAGAAAA";
    assert_eq!(find_tandem_repeat_units(seq, 0), (b"A".to_vec(), 2));
    assert_eq!(rust_cache(2), 40);
    assert_ne!(
        tandem_repeat_units(seq, 0),
        0,
        "canonical cell does not use cache[0]"
    );
}

#[test]
fn canonical_taagaaaa_after_cache_port() {
    let seq = b"TAAGAAAA";
    let (unit, n) = find_tandem_repeat_units(seq, 0);
    assert_eq!(unit, b"A");
    assert_eq!(n, 2);
    let bi = 44u8;
    let fill_q45 = GATK_PARITY_DEFAULT_INS_QUAL;
    assert_eq!(fill_q45, 45);
    assert_eq!(rust_cache(2), 40);
    assert_eq!(java_44_cache(2), 40);
    assert_eq!(bi.min(rust_cache(2)), 40);
    assert_eq!(fill_q45.min(rust_cache(2)), 40);
    // Pre-PCR sources still differ: BI 44 vs production fill 45.
    assert_ne!(bi, fill_q45);
}

#[test]
fn q45_fill_is_masked_by_q40_cache_on_canonical_repeat() {
    let seq = b"TAAGAAAA";
    let mut from_q45 = vec![GATK_PARITY_DEFAULT_INS_QUAL; seq.len()];
    let mut from_bi = vec![44u8; seq.len()];
    let mut d1 = from_q45.clone();
    let mut d2 = from_bi.clone();
    apply_pcr_error_model(seq, &mut from_q45, &mut d1, PcrErrorModel::Conservative);
    apply_pcr_error_model(seq, &mut from_bi, &mut d2, PcrErrorModel::Conservative);
    assert_eq!(from_q45[0], 40, "Q45 PCR-capped to cache[2]=40");
    assert_eq!(from_bi[0], 40, "BI 44 PCR-capped to cache[2]=40");
    assert_eq!(
        from_q45[0], from_bi[0],
        "final GOP masks the 45 vs 44 source"
    );
    // Last base is never PCR-written: source difference remains visible.
    assert_eq!(from_q45[seq.len() - 1], 45);
    assert_eq!(from_bi[seq.len() - 1], 44);
}

#[test]
fn anti_masking_repeat_4_selects_q39_not_q40() {
    // Java cache[3]=40, cache[4]=39. GOP 44 makes the min visible.
    assert_eq!(java_44_cache(3), 40);
    assert_eq!(java_44_cache(4), 39);
    assert_eq!(rust_cache(4), 39);
    let seq = b"AAAATAAAA";
    assert_eq!(tandem_repeat_units(seq, 4), 4);
    let mut ins = vec![44u8; seq.len()];
    let mut del = ins.clone();
    apply_pcr_error_model(seq, &mut ins, &mut del, PcrErrorModel::Conservative);
    assert_eq!(ins[4], 39);
    assert_ne!(ins[4], 40);
}

#[test]
fn min_selects_cache_when_input_gop_is_above() {
    for &(r, expect) in &[(1u8, 40u8), (2, 40), (3, 40), (4, 39), (5, 39), (10, 38)] {
        let input = 44u8;
        assert!(input > expect);
        assert_eq!(input.min(rust_cache(r as usize)), expect);
        assert_eq!(input.min(java_44_cache(r as usize)), expect);
    }
}
