//! 6R.98 coordinate-free: Java `--native-pair-hmm-use-double-precision`
//! result-buffer contract vs default float store.
//!
//! Java 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`) + GKL **0.8.8**.
//! Do not inspect PairHMM DP / SIMD arithmetic. No production precision change.
//!
//! Default (`useDoublePrecision=false`):
//! ```text
//! result_float = g_compute_full_prob_float(...)
//! result_final = (double)(log10f(result_float) - g_ctxf.LOG10_INITIAL_CONSTANT)
//! ```
//!
//! Double (`useDoublePrecision=true`): GKL sets `g_use_double`, so
//! `result_float = 0.0f` which is `< MIN_ACCEPTED` (`1e-28f`), forcing:
//! ```text
//! result_double = g_compute_full_prob_double(...)
//! result_final = log10(result_double) - g_ctxd.LOG10_INITIAL_CONSTANT
//! javaResults[i] = result_final   // still jdoubleArray / double[]
//! ```
//! The switch changes the **native numerical backend**, not only the Java array type.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r98_java_double_precision_result_contract
//! HOLDOUT_6R98=1 cargo test -p gatk-haplotypecaller --test holdout_6r98_java_double_precision -- --nocapture
//! ```

use std::collections::HashMap;

const JAVA_FLOAT_DUMP: &str = include_str!("6r96_java_seq6_post_kernel.tsv");
const JAVA_DOUBLE_DUMP: &str = include_str!("6r98_java_seq6_double_post_kernel.tsv");

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Classification {
    JavaFloatResultPath,
    JavaDoubleMatchesRust,
    JavaDoubleStillDiverges,
    PrecisionSwitchNotActive,
    NoProvenRelationship,
}

fn is_exact_f32_widened(x: f64) -> bool {
    (f64::from(x as f32)).to_bits() == x.to_bits()
}

fn ordered_f64_bits(x: f64) -> i64 {
    let bits = x.to_bits() as i64;
    if bits < 0 {
        i64::MIN - bits
    } else {
        bits
    }
}

fn ulp_distance(a: f64, b: f64) -> u64 {
    ordered_f64_bits(a).abs_diff(ordered_f64_bits(b))
}

/// GKL 0.8.8: `g_use_double ? 0.0f : g_compute_full_prob_float`. Zero is always
/// `< MIN_ACCEPTED` (1e-28f), so the double flag always takes the double-retry store.
fn gkl_use_double_forces_double_backend(use_double: bool, min_accepted: f32) -> bool {
    let result_float = if use_double { 0.0f32 } else { 1.0f32 };
    use_double && result_float < min_accepted
}

fn classify(
    j_float_all_f32_wide: bool,
    j_double_all_f32_wide: bool,
    j_float_eq_j_double: bool,
    j_double_eq_rust: bool,
) -> Classification {
    if j_double_all_f32_wide && j_float_eq_j_double {
        Classification::PrecisionSwitchNotActive
    } else if j_double_eq_rust {
        Classification::JavaDoubleMatchesRust
    } else if j_float_all_f32_wide && !j_double_all_f32_wide && !j_float_eq_j_double {
        Classification::JavaDoubleStillDiverges
    } else if j_float_all_f32_wide && j_double_eq_rust {
        Classification::JavaFloatResultPath
    } else {
        Classification::NoProvenRelationship
    }
}

fn parse_kv(s: &str) -> HashMap<&str, &str> {
    let mut d = HashMap::new();
    for part in s.split('\t') {
        if let Some((k, v)) = part.split_once('=') {
            d.insert(k, v);
        }
    }
    d
}

fn dump_cell_bits(text: &str, prefix: &str) -> Vec<u64> {
    let mut bits = Vec::new();
    for line in text.lines() {
        let p: Vec<&str> = line.split('\t').collect();
        if p.len() < 4 || p[0] != prefix {
            continue;
        }
        if !p[2].starts_with("rowbits_") {
            continue;
        }
        let val = p[3..].join("\t");
        let kv = parse_kv(&val);
        for s in kv.get("bits").copied().unwrap_or("").split(',') {
            if s.is_empty() {
                continue;
            }
            bits.push(u64::from_str_radix(s, 16).unwrap_or(0));
        }
    }
    bits
}

fn dump_buffer_kv(text: &str, prefix: &str) -> HashMap<String, String> {
    let mut d = HashMap::new();
    for line in text.lines() {
        let p: Vec<&str> = line.split('\t').collect();
        if p.len() < 4 || p[0] != prefix {
            continue;
        }
        d.insert(p[2].to_string(), p[3].to_string());
    }
    d
}

#[test]
fn forensic_6r98_gkl_double_flag_forces_double_backend() {
    const MIN_ACCEPTED: f32 = 1e-28;
    assert!(gkl_use_double_forces_double_backend(true, MIN_ACCEPTED));
    assert!(!gkl_use_double_forces_double_backend(false, MIN_ACCEPTED));
}

#[test]
fn forensic_6r98_default_store_is_f32_widen_double_store_is_not() {
    let float_store = f64::from(-16.542_921_f32);
    assert!(is_exact_f32_widened(float_store));
    let double_like: f64 = -16.542_921_340_320_96;
    assert!(!is_exact_f32_widened(double_like));
}

#[test]
fn forensic_6r98_java_float_dump_is_f32_wide_double_dump_is_not() {
    let float_bits = dump_cell_bits(JAVA_FLOAT_DUMP, "6R96");
    let double_bits = dump_cell_bits(JAVA_DOUBLE_DUMP, "6R98");
    assert_eq!(float_bits.len(), 153 * 70);
    assert_eq!(double_bits.len(), 153 * 70);
    let float_wide = float_bits
        .iter()
        .filter(|&&b| is_exact_f32_widened(f64::from_bits(b)))
        .count();
    let double_wide = double_bits
        .iter()
        .filter(|&&b| is_exact_f32_widened(f64::from_bits(b)))
        .count();
    assert_eq!(float_wide, float_bits.len());
    assert_eq!(double_wide, 0);
    let exact_eq = float_bits
        .iter()
        .zip(double_bits.iter())
        .filter(|(a, b)| a == b)
        .count();
    assert_eq!(exact_eq, 0);
}

#[test]
fn forensic_6r98_precision_switch_is_active_on_full_matrix() {
    let kv = dump_buffer_kv(JAVA_DOUBLE_DUMP, "6R98");
    assert_eq!(
        kv.get("use_double_precision").map(String::as_str),
        Some("true")
    );
    assert_eq!(
        kv.get("kernel_buffer_f32_wide").map(String::as_str),
        Some("0")
    );
    assert_eq!(kv.get("matrix_f32_wide").map(String::as_str), Some("0"));
    assert_eq!(
        kv.get("buffer_matrix_sorted_bits_equal")
            .map(String::as_str),
        Some("true")
    );
    assert_eq!(
        kv.get("pairhmm_class").map(String::as_str),
        Some("VectorLoglessPairHMM")
    );
}

#[test]
fn forensic_6r98_result_buffer_type_is_still_double_array() {
    let kv = dump_buffer_kv(JAVA_DOUBLE_DUMP, "6R98");
    assert_eq!(
        kv.get("kernel_buffer_type").map(String::as_str),
        Some("double[] mLogLikelihoodArray")
    );
}

#[test]
fn forensic_6r98_ulp_zero_iff_bits_equal() {
    let a: f64 = -16.542_921_340_320_96;
    assert_eq!(ulp_distance(a, a), 0);
    assert!(ulp_distance(a, f64::from(a as f32)) > 0);
}

#[test]
fn forensic_6r98_classify_switch_inactive_when_double_still_f32_wide() {
    assert_eq!(
        classify(true, true, true, false),
        Classification::PrecisionSwitchNotActive
    );
}

#[test]
fn forensic_6r98_classify_match_when_double_equals_rust() {
    assert_eq!(
        classify(true, false, false, true),
        Classification::JavaDoubleMatchesRust
    );
}

#[test]
fn forensic_6r98_classify_diverge_when_switch_active_and_rust_differs() {
    assert_eq!(
        classify(true, false, false, false),
        Classification::JavaDoubleStillDiverges
    );
}
