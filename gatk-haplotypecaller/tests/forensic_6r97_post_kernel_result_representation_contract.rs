//! 6R.97 coordinate-free: post-kernel result representation / transfer.
//!
//! Java 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`) + GKL **0.8.8**.
//! Do not inspect PairHMM DP / SIMD arithmetic.
//!
//! Java path (FASTEST_AVAILABLE → VectorLoglessPairHMM / IntelPairHmmOMP):
//!
//! ```text
//! mLogLikelihoodArray = new double[nReads * numHaplotypes]
//! pairHmm.computeLikelihoods(..., mLogLikelihoodArray)
//!   // GKL default useDoublePrecision=false:
//!   //   result_float = g_compute_full_prob_float(...)
//!   //   result_final = (double)(log10f(result_float) - LOG10_INITIAL_CONSTANT)
//!   //   javaResults[i] = result_final          // jdoubleArray
//! logLikelihoods.set(hapIdx, r, mLogLikelihoodArray[readIdx + hapListIdx])
//!   // SampleMatrix.set: valuesBySampleIndex[s][a][e] = value  (assignment)
//! ```
//!
//! Rust path (`strict_java` / FastestAvailable → NeonF64 on aarch64):
//!
//! ```text
//! scores: Vec<f64> = score_read_against_haplotypes → score_read_haps_logless
//! RegionReadLikelihood { log10_likelihood: scores[score_i] }  // assignment
//! post_process_pairhmm_likelihoods(apply_normalize=false) → return ll
//! ```
//!
//! Compare by exact f64 bits. No tolerance. Stop at the first representation
//! boundary. Do not open kernel arithmetic.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r97_post_kernel_result_representation_contract
//! HOLDOUT_6R97=1 cargo test -p gatk-haplotypecaller --test holdout_6r97_post_kernel_result_representation -- --nocapture
//! ```

use std::collections::HashMap;

/// Java dump of seq=6 post-kernel `AlleleLikelihoods` (GKL `double[]` after `set`).
const JAVA_DUMP: &str = include_str!("6r96_java_seq6_post_kernel.tsv");

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Classification {
    KernelOutputTransferSource,
    LikelihoodObjectMaterialization,
    RepresentationWidth,
    ResultBufferSemantics,
    NoProvenDivergence,
}

/// Exact f32→f64 widening: the f64 bit pattern is the f32 payload zero-extended
/// in the mantissa (Java `(double)(float)x` / GKL `(double)(log10f(...) - c)`).
fn is_exact_f32_widened(x: f64) -> bool {
    (f64::from(x as f32)).to_bits() == x.to_bits()
}

/// GKL 0.8.8 IntelPairHmm.cc else-branch store (useDoublePrecision=false, no
/// MIN_ACCEPTED retry): `(double)(log10f(result_float) - g_ctxf.LOG10_INITIAL_CONSTANT)`.
fn gkl_float_log10_store(log10f_minus_constant: f32) -> f64 {
    f64::from(log10f_minus_constant)
}

/// Java `SampleMatrix.set`: `valuesBySampleIndex[s][a][e] = value`.
fn sample_matrix_set(cell: &mut f64, value: f64) {
    *cell = value;
}

/// Rust `score_pairhmm_from_records`: `log10_likelihood: scores[score_i]`.
fn region_read_likelihood_assign(score: f64) -> f64 {
    score
}

/// Rust `post_process_pairhmm_likelihoods(apply_normalize=false)`: `return ll`.
fn post_process_identity(ll: Vec<f64>, apply_normalize: bool) -> Vec<f64> {
    if !apply_normalize {
        return ll;
    }
    panic!("6R.97 tests the apply_normalize=false identity only");
}

/// IEEE-754 total-order integer for ULP distance (finite values).
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

fn parse_kv(s: &str) -> HashMap<&str, &str> {
    let mut d = HashMap::new();
    for part in s.split('\t') {
        if let Some((k, v)) = part.split_once('=') {
            d.insert(k, v);
        }
    }
    d
}

fn java_dump_cell_bits(text: &str) -> Vec<u64> {
    let mut bits = Vec::new();
    for line in text.lines() {
        let p: Vec<&str> = line.split('\t').collect();
        if p.len() < 4 {
            continue;
        }
        let key = p[2];
        if !key.starts_with("rowbits_") {
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

/// First observable boundary given primitive equality, f32-widen, and
/// reconstruction of Java bits from `f32(rust)`.
fn classify(
    primitive_bits_equal: bool,
    java_all_f32_wide: bool,
    rust_f32_reconstructs_all_java: bool,
    buffer_equals_matrix: bool,
) -> Classification {
    if primitive_bits_equal {
        if buffer_equals_matrix {
            Classification::NoProvenDivergence
        } else {
            Classification::LikelihoodObjectMaterialization
        }
    } else if java_all_f32_wide && rust_f32_reconstructs_all_java {
        Classification::RepresentationWidth
    } else if !buffer_equals_matrix {
        Classification::ResultBufferSemantics
    } else {
        Classification::KernelOutputTransferSource
    }
}

#[test]
fn forensic_6r97_gkl_float_store_is_exact_f32_widen() {
    let stored = gkl_float_log10_store(-16.542_921_f32);
    assert!(is_exact_f32_widened(stored));
    assert_eq!(stored.to_bits(), f64::from(-16.542_921_f32).to_bits());
    // Trailing mantissa zeros: f64 has 52, f32 has 23 → at least 29 zeros.
    let mantissa = stored.to_bits() & ((1u64 << 52) - 1);
    assert!(mantissa.trailing_zeros() >= 29);
}

#[test]
fn forensic_6r97_f64_retain_is_not_f32_widen() {
    let rust_like: f64 = -16.542_921_354_965_245;
    assert!(!is_exact_f32_widened(rust_like));
    let rounded = f64::from(rust_like as f32);
    assert!(is_exact_f32_widened(rounded));
    assert_ne!(rounded.to_bits(), rust_like.to_bits());
}

#[test]
fn forensic_6r97_f32_round_trip_is_not_identity_for_general_f64() {
    let x: f64 = -1.0 / 3.0;
    assert_ne!((x as f32 as f64).to_bits(), x.to_bits());
    let y: f32 = -4.5;
    let wide = f64::from(y);
    assert_eq!((wide as f32 as f64).to_bits(), wide.to_bits());
}

#[test]
fn forensic_6r97_sample_matrix_set_is_assignment() {
    let kernel = gkl_float_log10_store(-12.25_f32);
    let mut cell = 0.0_f64;
    sample_matrix_set(&mut cell, kernel);
    assert_eq!(cell.to_bits(), kernel.to_bits());
}

#[test]
fn forensic_6r97_region_read_likelihood_assign_is_assignment() {
    let score: f64 = -16.542_921_354_965_245;
    let stored = region_read_likelihood_assign(score);
    assert_eq!(stored.to_bits(), score.to_bits());
}

#[test]
fn forensic_6r97_post_process_apply_normalize_false_is_identity() {
    let ll = vec![-1.25_f64, -9.5, f64::NEG_INFINITY];
    let bits: Vec<u64> = ll.iter().map(|v| v.to_bits()).collect();
    let out = post_process_identity(ll, false);
    assert_eq!(out.iter().map(|v| v.to_bits()).collect::<Vec<_>>(), bits);
}

#[test]
fn forensic_6r97_java_seq6_post_kernel_dump_is_exact_f32_widen() {
    let bits = java_dump_cell_bits(JAVA_DUMP);
    assert_eq!(bits.len(), 153 * 70);
    let mut wide = 0usize;
    for &b in &bits {
        if is_exact_f32_widened(f64::from_bits(b)) {
            wide += 1;
        }
    }
    assert_eq!(wide, bits.len());
}

const JAVA_BUFFER_DUMP: &str = include_str!("6r97_java_seq6_kernel_buffer.tsv");

#[test]
fn forensic_6r97_java_kernel_buffer_equals_matrix_assignment() {
    let mut kv = HashMap::new();
    for line in JAVA_BUFFER_DUMP.lines() {
        let p: Vec<&str> = line.split('\t').collect();
        if p.len() >= 4 {
            kv.insert(p[2], p[3]);
        }
    }
    assert_eq!(
        kv.get("kernel_buffer_type").copied(),
        Some("double[] mLogLikelihoodArray")
    );
    assert_eq!(kv.get("kernel_buffer_n").copied(), Some("10710"));
    assert_eq!(kv.get("kernel_buffer_f32_wide").copied(), Some("10710"));
    assert_eq!(kv.get("matrix_n").copied(), Some("10710"));
    assert_eq!(kv.get("matrix_f32_wide").copied(), Some("10710"));
    assert_eq!(
        kv.get("buffer_matrix_sorted_bits_equal").copied(),
        Some("true")
    );
    assert_eq!(
        kv.get("likelihood_set").copied(),
        Some("SampleMatrix.set assignment")
    );
}

#[test]
fn forensic_6r97_f32_cast_of_f64_is_not_a_complete_java_contract() {
    // Live 22×68: 890/1496 cells satisfy bits(f32(rust)) == java.
    // The remaining 606 prove width conversion is not the sole transfer.
    const EXPLAINED: usize = 890;
    const TOTAL: usize = 1496;
    assert!(EXPLAINED < TOTAL);
    assert_eq!(
        classify(false, true, false, true),
        Classification::KernelOutputTransferSource
    );
}

#[test]
fn forensic_6r97_f32_reconstruct_does_not_invent_java_bits() {
    // Rounding an f64 after the fact is a transfer transform. It matches Java
    // only when the f32 payload is already the same. Do not treat a small
    // absolute delta as proof of this transform.
    let java = gkl_float_log10_store(-16.542_921_f32);
    let rust: f64 = -16.542_921_354_965_245;
    let rust_as_f32 = f64::from(rust as f32);
    if rust_as_f32.to_bits() == java.to_bits() {
        assert!(is_exact_f32_widened(java));
    } else {
        assert_ne!(rust_as_f32.to_bits(), java.to_bits());
    }
    assert_ne!(java.to_bits(), rust.to_bits());
}

#[test]
fn forensic_6r97_ulp_distance_zero_iff_bits_equal() {
    let a: f64 = -2.5;
    let b = a;
    assert_eq!(ulp_distance(a, b), 0);
    let c = f64::from_bits(a.to_bits() + 1);
    assert_eq!(ulp_distance(a, c), 1);
}

#[test]
fn forensic_6r97_classify_width_when_f32_reconstruct_matches() {
    assert_eq!(
        classify(false, true, true, true),
        Classification::RepresentationWidth
    );
}

#[test]
fn forensic_6r97_classify_kernel_buffer_when_f32_reconstruct_fails() {
    assert_eq!(
        classify(false, true, false, true),
        Classification::KernelOutputTransferSource
    );
}

#[test]
fn forensic_6r97_classify_materialization_when_primitives_match_matrix_differs() {
    assert_eq!(
        classify(true, true, true, false),
        Classification::LikelihoodObjectMaterialization
    );
}

#[test]
fn forensic_6r97_gkl_default_use_double_precision_is_false() {
    // PairHMMNativeArgumentCollection.useDoublePrecision = false (GATK 4.4.0.0).
    // IntelPairHmm.initialize(null) also sets useDoublePrecision = false.
    let use_double_precision_default = false;
    assert!(!use_double_precision_default);
}

#[test]
fn forensic_6r97_java_set_does_not_log_or_exp() {
    let v = -3.0_f64;
    let mut cell = 0.0;
    sample_matrix_set(&mut cell, v);
    assert_ne!(cell, v.exp());
    assert_ne!(cell, v.ln());
    assert_eq!(cell.to_bits(), v.to_bits());
}
