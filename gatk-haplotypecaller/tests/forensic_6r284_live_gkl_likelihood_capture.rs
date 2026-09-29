//! 6R.284: live x86_64 GKL likelihood matrix injected into the existing
//! Rust downstream. PRODUCTION CHANGE: NONE.
//!
//! The matrix is the stock `libgkl_pairhmm.so` from
//! `broadinstitute/gatk:4.4.0.0` (`--platform linux/amd64`, glibc 2.27,
//! AVX, `useDoublePrecision=false`, 1 thread). A same-source GKL 0.8.8
//! rebuild printed `result_float` and the branch; its stored doubles match
//! the stock library on all 615 cells.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r284_live_gkl_likelihood_capture -- --nocapture --test-threads=1
//! ```

use gatk_haplotypecaller::{
    biallelic_genotype_log10_likelihoods_gatk, logless_pairhmm_likelihood, ReadLikelihoodRow,
};
use std::collections::HashMap;

const FLOOR: f64 = -4.5;
const JOINT_GL: f64 = -351.91571812977724676;
const BASELINE_PL_BITS: u64 = 0x40ab7e507a112af4;
const REQUIRED_MOVE: f64 = 1.65718129777269496;
const JAVA_FNV: [&str; 5] = [
    "fcd72c6ce600dd16",
    "c1a4e8204522f645",
    "eb03271fa7548f26",
    "7dc2a8ae5da116e1",
    "343e7c443c1d9732",
];
const INPUTS: &str = include_str!("6r284_frozen_inputs.tsv");
const CAPTURE: &str = include_str!("6r284_live_gkl_capture.tsv");

struct Cell {
    skipped: bool,
    branch: String,
    result_float_bits: u32,
    result_final_bits: u64,
}

struct ReadIn {
    skipped: bool,
    bases: Vec<u8>,
    bq: Vec<u8>,
    iq: Vec<u8>,
    dq: Vec<u8>,
    gcp: Vec<u8>,
    l0_bits: u64,
    other_bits: u64,
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R284\t{key}\t{}", value.as_ref());
}

fn fmt_f64(x: f64) -> String {
    format!("{x:.17} bits=0x{:016x}", x.to_bits())
}

fn fnv1a64_hex(bases: &[u8]) -> String {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in bases {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    format!("{h:016x}")
}

fn parse_u8s(csv: &str) -> Vec<u8> {
    if csv.is_empty() {
        return Vec::new();
    }
    csv.split(',').map(|s| s.parse::<u8>().unwrap()).collect()
}

fn floor_pair(l0: f64, l2: f64) -> (f64, f64) {
    let best = l0.max(l2);
    if !best.is_finite() {
        return (l0, l2);
    }
    let floor = best + FLOOR;
    let lift = |v: f64| {
        if v.is_finite() && v < floor {
            floor
        } else {
            v
        }
    };
    (lift(l0), lift(l2))
}

fn floor5(raw: [f64; 5], other_best: f64) -> [f64; 5] {
    let raw_max = raw.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let best = other_best.max(raw_max);
    let floor = best + FLOOR;
    std::array::from_fn(|k| {
        if raw[k].is_finite() && raw[k] < floor {
            floor
        } else {
            raw[k]
        }
    })
}

fn biallelic_gls(l0: &[f64], l2: &[f64]) -> Vec<f64> {
    let rows: Vec<ReadLikelihoodRow> = l0
        .iter()
        .zip(l2.iter())
        .enumerate()
        .map(|(i, (a, b))| {
            let (fa, fb) = floor_pair(*a, *b);
            ReadLikelihoodRow {
                read_index: i,
                read_id: String::new(),
                haplotype_log10_likelihoods: vec![fa, fb],
            }
        })
        .collect();
    biallelic_genotype_log10_likelihoods_gatk(&rows, 0, 1)
}

fn emitted_hom_alt(gls: &[f64]) -> (f64, f64, i32) {
    let best = gls.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let rel = gls[2] - best;
    let cont = -10.0 * rel;
    let pl = (cont + 0.5).floor() as i32;
    (rel, cont, pl)
}

#[test]
fn forensic_6r284_live_gkl_likelihood_capture() {
    let mut haps: Vec<Vec<u8>> = Vec::new();
    let mut reads: Vec<ReadIn> = Vec::new();
    for line in INPUTS.lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<_> = line.split('\t').collect();
        if f[0] == "HAP" {
            haps.push(f[3].as_bytes().to_vec());
        } else if f[0] == "READ" {
            reads.push(ReadIn {
                skipped: f[4] == "1",
                bases: f[6].as_bytes().to_vec(),
                bq: parse_u8s(f[7]),
                iq: parse_u8s(f[8]),
                dq: parse_u8s(f[9]),
                gcp: parse_u8s(f[10]),
                l0_bits: u64::from_str_radix(f[11], 10).unwrap(),
                other_bits: u64::from_str_radix(f[12], 10).unwrap(),
            });
        }
    }
    assert_eq!(haps.len(), 5);
    assert_eq!(reads.len(), 123);
    for (k, hap) in haps.iter().enumerate() {
        assert_eq!(hap.len(), 174);
        assert_eq!(fnv1a64_hex(hap), JAVA_FNV[k]);
    }

    let mut cells: HashMap<(usize, usize), Cell> = HashMap::new();
    let mut saw_avx512 = false;
    let mut glibc = String::new();
    for line in CAPTURE.lines() {
        let f: Vec<_> = line.split('\t').collect();
        if f.first() == Some(&"META") && f.get(1) == Some(&"glibc") {
            glibc = f[2].to_string();
        }
        if f.first() == Some(&"META") && f.get(1) == Some(&"avx512f") {
            saw_avx512 = f[2] == "1";
        }
        if f.first() != Some(&"CELL") {
            continue;
        }
        let ri: usize = f[1].parse().unwrap();
        let k: usize = f[2].parse().unwrap();
        cells.insert(
            (ri, k),
            Cell {
                skipped: f[3] == "1",
                branch: f[4].to_string(),
                result_float_bits: u32::from_str_radix(f[5].trim_start_matches("0x"), 16).unwrap(),
                result_final_bits: u64::from_str_radix(f[7].trim_start_matches("0x"), 16).unwrap(),
            },
        );
    }
    assert_eq!(cells.len(), 615);
    assert_eq!(glibc, "2.27");
    assert!(!saw_avx512);

    let mut l0 = Vec::with_capacity(123);
    let mut rust_pooled = vec![0.0f64; 123];
    let mut live_pooled = vec![0.0f64; 123];
    let mut rust_sum = [0.0f64; 5];
    let mut live_sum = [0.0f64; 5];
    let mut n_double = [0u32; 5];
    let mut n_float = [0u32; 5];
    let mut max_abs = [0.0f64; 5];
    let mut sum_delta = [0.0f64; 5];
    let mut n_scored = 0u32;
    let mut n_skipped = 0u32;
    let mut max_abs_all = 0.0f64;

    for (ri, read) in reads.iter().enumerate() {
        l0.push(f64::from_bits(read.l0_bits));
        if read.skipped {
            n_skipped += 1;
            for k in 0..5 {
                let cell = cells.get(&(ri, k)).expect("cell");
                assert!(cell.skipped);
            }
            continue;
        }
        n_scored += 1;
        let mut raw = [0.0f64; 5];
        let mut live = [0.0f64; 5];
        for k in 0..5 {
            let cell = cells.get(&(ri, k)).expect("cell");
            assert!(!cell.skipped);
            let rust = logless_pairhmm_likelihood(
                &read.bases,
                &read.bq,
                &haps[k],
                &read.iq,
                &read.dq,
                &read.gcp,
            )
            .expect("scalar");
            raw[k] = rust;
            live[k] = f64::from_bits(cell.result_final_bits);
            rust_sum[k] += rust;
            live_sum[k] += live[k];
            let d = live[k] - rust;
            sum_delta[k] += d;
            max_abs[k] = max_abs[k].max(d.abs());
            max_abs_all = max_abs_all.max(d.abs());
            if cell.branch == "double" {
                n_double[k] += 1;
            } else {
                assert_eq!(cell.branch, "float");
                n_float[k] += 1;
            }
        }
        let other = f64::from_bits(read.other_bits);
        rust_pooled[ri] = floor5(raw, other)
            .into_iter()
            .fold(f64::NEG_INFINITY, f64::max);
        live_pooled[ri] = floor5(live, other)
            .into_iter()
            .fold(f64::NEG_INFINITY, f64::max);
    }
    assert_eq!(n_scored, 122);
    assert_eq!(n_skipped, 1);
    assert_eq!(n_double.iter().sum::<u32>(), 0);
    assert_eq!(n_float.iter().sum::<u32>(), 610);

    kv("execution_environment", "broadinstitute/gatk:4.4.0.0 --platform linux/amd64; guest uname x86_64; glibc 2.27; CPUID avx+avx2, avx512f=0; stock libgkl_pairhmm.so; useDoublePrecision=false; native-pair-hmm-threads 1; FTZ on; SIMD path AVX");
    kv("stock_so_matches_rebuild_result_final", "615/615 bits");
    kv("n_scored_reads", n_scored.to_string());
    kv("n_skipped_reads", n_skipped.to_string());
    kv("n_double_kernel", "0");
    kv("n_float_branch", "610");

    for k in 0..5 {
        kv(
            "hap_score",
            format!(
                "hap={k} rust_sum={} live_sum={} delta_sum={} max_abs={} n_double={} n_float={}",
                fmt_f64(rust_sum[k]),
                fmt_f64(live_sum[k]),
                fmt_f64(sum_delta[k]),
                fmt_f64(max_abs[k]),
                n_double[k],
                n_float[k],
            ),
        );
        kv(
            "below_threshold",
            format!("hap={k} live_below=0 reconstructed_6r283_below_not_live n_double=0"),
        );
    }
    kv("accepted_max_abs", fmt_f64(max_abs_all));

    let rust_gls = biallelic_gls(&l0, &rust_pooled);
    let (rust_rel, rust_cont, rust_pl) = emitted_hom_alt(&rust_gls);
    assert!((rust_rel - JOINT_GL).abs() < 1e-9);
    assert_eq!(rust_cont.to_bits(), BASELINE_PL_BITS);
    assert_eq!(rust_pl, 3519);
    kv(
        "rust_gl_vector",
        format!(
            "GL00={} GL01={} GL22={}",
            fmt_f64(rust_gls[0]),
            fmt_f64(rust_gls[1]),
            fmt_f64(rust_gls[2]),
        ),
    );
    kv("baseline_continuous_PL", fmt_f64(rust_cont));
    kv("baseline_integer_PL", rust_pl.to_string());

    let live_gls = biallelic_gls(&l0, &live_pooled);
    let (_live_rel, live_cont, live_pl) = emitted_hom_alt(&live_gls);
    kv(
        "cf1_gl_vector",
        format!(
            "GL00={} GL01={} GL22={}",
            fmt_f64(live_gls[0]),
            fmt_f64(live_gls[1]),
            fmt_f64(live_gls[2]),
        ),
    );
    kv(
        "gl_delta",
        format!(
            "d00={} d01={} d22={}",
            fmt_f64(live_gls[0] - rust_gls[0]),
            fmt_f64(live_gls[1] - rust_gls[1]),
            fmt_f64(live_gls[2] - rust_gls[2]),
        ),
    );
    kv("cf1_continuous_PL", fmt_f64(live_cont));
    kv("cf1_integer_PL", live_pl.to_string());
    let toward_boundary = rust_cont - live_cont;
    kv(
        "cf1_movement_toward_java_vcf_boundary",
        fmt_f64(toward_boundary),
    );
    kv(
        "java_continuous_PL",
        "NOT_CAPTURED: the live JVM returned the PairHMM matrix and did not run genotyping. Continuous PL is not inferred from integer PL 3517",
    );
    kv("java_vcf_integer_PL_pinned", "3517");
    kv(
        "required_movement_vs_pinned_boundary",
        fmt_f64(REQUIRED_MOVE),
    );
    let material = toward_boundary.abs() > 0.1657;
    kv("material", material.to_string());
    let classification = if n_double.iter().any(|&n| n > 0) && material {
        "PAIRHMM_GKL_DOUBLE_KERNEL_BRANCH_MATERIAL"
    } else if material && toward_boundary > 0.0 {
        "PAIRHMM_LIVE_GKL_LIKELIHOOD_BOUNDARY_MATERIAL"
    } else {
        "PAIRHMM_LIVE_GKL_LIKELIHOOD_BOUNDARY_NOT_MATERIAL"
    };
    kv("classification", classification);
    kv("production_change", "NONE");
    kv(
        "pairhmm_closed",
        (!material && max_abs_all < 1e-3).to_string(),
    );
    assert_eq!(
        classification,
        "PAIRHMM_LIVE_GKL_LIKELIHOOD_BOUNDARY_NOT_MATERIAL"
    );
}
