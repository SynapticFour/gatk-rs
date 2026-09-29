//! 6R.302: design check for the Java-compatible heterozygote combine.
//! Production code is not modified. The candidate evaluation calls
//! `approximate_log10_sum_log10_pair`, which is the primitive the
//! diploid caller would use.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r302_jacobian_design -- --nocapture --test-threads=1
//! ```

use gatk_haplotypecaller::approximate_log10_sum_log10_pair;
use std::path::{Path, PathBuf};

const ORACLE: &str = "gatk-haplotypecaller/tests/6r301_jacobian_oracle.tsv";
const ALLELE_TSV: &str = "gatk-haplotypecaller/tests/6r297_java_hap_ll.tsv";

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R302\t{key}\t{}", value.as_ref());
}

fn fmt_f(x: f64) -> String {
    format!("{x:.17}")
}

fn cell_bits(cell: &str) -> f64 {
    let hex = cell.split("bits=").nth(1).unwrap().trim();
    f64::from_bits(u64::from_str_radix(hex.trim_start_matches("0x"), 16).unwrap())
}

fn runtime_table(k: i32) -> f64 {
    let exponent = (-(k as i32) as f64) * 1.0e-4;
    (1.0 + 10.0_f64.powf(exponent)).log10()
}

/// Java `MathUtils.fastRound` for the documented probe values.
fn java_fast_round(d: f64) -> i32 {
    if d > 0.0 {
        (d + 0.5) as i32
    } else {
        (d - 0.5) as i32
    }
}

fn java_round_non_negative(x: f64) -> i32 {
    debug_assert!(x >= 0.0);
    (x + 0.5).floor() as i32
}

#[test]
fn forensic_6r302_jacobian_design() {
    let oracle = std::fs::read_to_string(repo_root().join(ORACLE)).unwrap();
    let mut mismatches = Vec::new();
    let mut n_pairs = 0usize;
    for line in oracle.lines() {
        let p: Vec<_> = line.split('\t').collect();
        if p.len() < 11 || p[0] != "6R301" || p[1] != "pair" || p[8] != "table" {
            continue;
        }
        n_pairs += 1;
        let index: i32 = p[9].parse().unwrap();
        let java_cell = cell_bits(p[10]);
        let rust_cell = runtime_table(index);
        if rust_cell.to_bits() != java_cell.to_bits() {
            let a = cell_bits(p[5]);
            let b = cell_bits(p[6]);
            let java_result = cell_bits(p[7]);
            let helper = approximate_log10_sum_log10_pair(a, b);
            let smaller = if a > b { b } else { a };
            let larger = if a > b { a } else { b };
            mismatches.push((
                p[2].to_string(),
                p[3].to_string(),
                p[4].to_string(),
                index,
                larger - smaller,
                java_cell,
                rust_cell,
                java_result,
                helper,
            ));
        }
    }
    assert_eq!(n_pairs, 366);
    assert_eq!(mismatches.len(), 11);
    for (gt, qname, flags, index, diff, java_cell, rust_cell, java_result, helper) in &mismatches {
        assert_eq!(helper.to_bits(), java_result.to_bits());
        let ulps = rust_cell.to_bits().abs_diff(java_cell.to_bits());
        kv(
            "table_mismatch",
            format!(
                "genotype={gt} read={qname} flags={flags} index={index} diff={} java_bits=0x{:016x} rust_bits=0x{:016x} java={} rust={} ulps={ulps} per_read_result_delta={}",
                fmt_f(*diff),
                java_cell.to_bits(),
                rust_cell.to_bits(),
                fmt_f(*java_cell),
                fmt_f(*rust_cell),
                fmt_f(helper - java_result)
            ),
        );
    }

    let alleles = std::fs::read_to_string(repo_root().join(ALLELE_TSV)).unwrap();
    let mut java_gl = [0.0; 6];
    let mut java_pl = [0.0; 6];
    for line in alleles.lines() {
        let p: Vec<_> = line.split('\t').collect();
        if p.len() >= 8 && p[0] == "6R297" && p[1] == "java_gl_vector" {
            java_gl = [2, 3, 4, 5, 6, 7].map(|i| cell_bits(p[i]));
        }
        if p.len() >= 8 && p[0] == "6R297" && p[1] == "java_continuous_pl_vector" {
            java_pl = [2, 3, 4, 5, 6, 7].map(|i| cell_bits(p[i]));
        }
    }
    let log10_2 = 2.0_f64.log10();
    let mut gl = [0.0; 6];
    let mut n = 0usize;
    let mut non_finite = 0usize;
    for line in alleles.lines() {
        let p: Vec<_> = line.split('\t').collect();
        if p.len() < 7 || p[0] != "6R297" || p[1] != "allele_ll" {
            continue;
        }
        let a = [cell_bits(p[4]), cell_bits(p[5]), cell_bits(p[6])];
        if a.iter().any(|v| !v.is_finite()) {
            non_finite += 1;
        }
        n += 1;
        gl[0] += a[0] + log10_2;
        gl[2] += a[1] + log10_2;
        gl[5] += a[2] + log10_2;
        gl[1] += approximate_log10_sum_log10_pair(a[0], a[1]);
        gl[3] += approximate_log10_sum_log10_pair(a[0], a[2]);
        gl[4] += approximate_log10_sum_log10_pair(a[1], a[2]);
    }
    assert_eq!(n, 122);
    assert_eq!(non_finite, 0);
    let den = n as f64 * log10_2;
    for g in &mut gl {
        *g -= den;
    }
    for g in 0..6 {
        assert_eq!(gl[g].to_bits(), java_gl[g].to_bits());
    }
    let best = gl.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let mut pl = [0.0; 6];
    let mut ipl = [0i32; 6];
    for i in 0..6 {
        pl[i] = -10.0 * (gl[i] - best);
        ipl[i] = java_round_non_negative(pl[i]);
        assert_eq!(pl[i].to_bits(), java_pl[i].to_bits());
    }
    assert_eq!([ipl[0], ipl[3], ipl[5]], [570, 0, 3517]);

    let probes = [
        (0.0, 0i32),
        (1.5, 2),
        (f64::from_bits(0x3ff7ffffffffffff), 1),
        (f64::from_bits(0x3ff8000000000001), 2),
        (1e-15, 0),
        (-1.2, -1),
        (-1.5, -2),
        (f64::from_bits(0xbff8000000000001), -2),
        (f64::from_bits(0xbff7ffffffffffff), -1),
    ];
    for (d, expected) in probes {
        assert_eq!(java_fast_round(d), expected);
    }
    let specials = [
        (f64::NEG_INFINITY, -1.5, 0xbff8000000000000u64),
        (-1.5, f64::NEG_INFINITY, 0xbff8000000000000),
        (f64::NEG_INFINITY, f64::NEG_INFINITY, 0xfff0000000000000),
        (-2.0, -1.0, 0xbfeeace93f5e7523),
        (-3.25, -3.25, 0xc007977d95ec10c0),
        (f64::NAN, -1.0, 0xbff0000000000000),
        (-1.0, f64::NAN, 0x7ff8000000000000),
        (f64::INFINITY, -1.0, 0x7ff0000000000000),
        (-1.0, f64::INFINITY, 0x7ff0000000000000),
        (f64::INFINITY, f64::INFINITY, 0x7ff0000000000000),
        (0.0, f64::from_bits(0x401fffffffffffff), 0x4020000000254e3c),
        (0.0, 8.0, 0x4020000000000000),
        (0.0, f64::from_bits(0x4020000000000001), 0x4020000000000001),
    ];
    for (a, b, bits) in specials {
        assert_eq!(approximate_log10_sum_log10_pair(a, b).to_bits(), bits);
    }

    kv("table_mismatch_count", mismatches.len().to_string());
    kv("candidate_gl", gl.map(fmt_f).join(","));
    kv("candidate_pl", pl.map(fmt_f).join(","));
    kv(
        "candidate_integer_pl",
        format!(
            "{},{},{},{},{},{}",
            ipl[0], ipl[1], ipl[2], ipl[3], ipl[4], ipl[5]
        ),
    );
    kv(
        "candidate_emitted_pl",
        format!("{},{},{}", ipl[0], ipl[3], ipl[5]),
    );
    kv("hom_alt_cf", fmt_f(pl[5]));
    kv(
        "classification",
        "JAVA_COMPATIBLE_JACOBIAN_DESIGN_CONFIRMED",
    );
    kv("production_change", "NONE");
}
