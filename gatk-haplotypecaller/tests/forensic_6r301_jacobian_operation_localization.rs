//! 6R.301: heterozygote combine operation, diagnostic only.
//! Production `log10_sum_log10` is not modified.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r301_jacobian_operation_localization -- --nocapture --test-threads=1
//! ```

use gatk_haplotypecaller::{approximate_log10_sum_log10_pair, log10_sum_log10};
use std::path::{Path, PathBuf};
use std::sync::OnceLock;

const ORACLE: &str = "gatk-haplotypecaller/tests/6r301_jacobian_oracle.tsv";
const ALLELE_TSV: &str = "gatk-haplotypecaller/tests/6r297_java_hap_ll.tsv";

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R301\t{key}\t{}", value.as_ref());
}

fn fmt_f(x: f64) -> String {
    format!("{x:.17}")
}

fn cell_bits(cell: &str) -> f64 {
    let hex = cell.split("bits=").nth(1).unwrap().trim();
    f64::from_bits(u64::from_str_radix(hex.trim_start_matches("0x"), 16).unwrap())
}

/// Literal transcription of GATK 4.4.0.0 `JacobianLogTable` static init.
/// `Math.log10(1.0 + Math.pow(10.0, (double)(-k) * 1.0E-4))`, k in `0..80001`.
fn literal_cache() -> &'static [f64] {
    static CACHE: OnceLock<Vec<f64>> = OnceLock::new();
    CACHE.get_or_init(|| {
        (0..80_001)
            .map(|k| {
                let exponent = (-(k as i32) as f64) * 1.0e-4;
                (1.0 + 10.0_f64.powf(exponent)).log10()
            })
            .collect()
    })
}

/// GATK 4.4.0.0 `MathUtils.fastRound` line 452. `d2i` of `d ± 0.5`.
fn literal_fast_round(d: f64) -> i32 {
    if d > 0.0 {
        (d + 0.5) as i32
    } else {
        (d - 0.5) as i32
    }
}

/// GATK 4.4.0.0 `JacobianLogTable.get` lines 442–444. No clamp.
fn literal_jacobian_get(difference: f64) -> f64 {
    let index = literal_fast_round(difference * 10_000.0);
    literal_cache()[index as usize]
}

/// GATK 4.4.0.0 `MathUtils.approximateLog10SumLog10(double, double)` lines 490–502.
fn literal_c(a: f64, b: f64) -> f64 {
    if a > b {
        return literal_c(b, a);
    }
    if a == f64::NEG_INFINITY {
        return b;
    }
    let diff = b - a;
    let addend = if diff < 8.0 {
        literal_jacobian_get(diff)
    } else {
        0.0
    };
    b + addend
}

fn production_a(a: f64, b: f64) -> f64 {
    log10_sum_log10(&[a, b])
}

struct Pair {
    genotype: String,
    qname: String,
    flags: u16,
    a: f64,
    b: f64,
    java: f64,
    branch: String,
    index: Option<i32>,
    table: Option<f64>,
}

struct GenotypeStat {
    max_abs: f64,
    sum: f64,
    n_diff_ac: usize,
    n_diff_bc: usize,
    n: usize,
    qname: String,
    flags: u16,
    a: f64,
    b: f64,
    out_a: f64,
    out_b: f64,
    out_c: f64,
    diff_inputs: f64,
}

fn empty_stat() -> GenotypeStat {
    GenotypeStat {
        max_abs: -1.0,
        sum: 0.0,
        n_diff_ac: 0,
        n_diff_bc: 0,
        n: 0,
        qname: String::new(),
        flags: 0,
        a: 0.0,
        b: 0.0,
        out_a: 0.0,
        out_b: 0.0,
        out_c: 0.0,
        diff_inputs: 0.0,
    }
}

#[test]
fn forensic_6r301_jacobian_operation_localization() {
    let text = std::fs::read_to_string(repo_root().join(ORACLE)).unwrap();
    let mut cache_len = 0usize;
    let mut cache_xor = 0u64;
    let mut sample_k = Vec::new();
    let mut pairs = Vec::new();
    let mut thresholds = Vec::new();
    for line in text.lines() {
        let p: Vec<_> = line.split('\t').collect();
        if p.first().copied() != Some("6R301") {
            continue;
        }
        match p[1] {
            "cache_len" => cache_len = p[2].parse().unwrap(),
            "cache_xor" => {
                cache_xor = u64::from_str_radix(p[2].trim_start_matches("0x"), 16).unwrap()
            }
            "cache_k" => sample_k.push((p[2].parse::<usize>().unwrap(), cell_bits(p[3]))),
            "pair" => {
                let index = if p[9] == "NA" {
                    None
                } else {
                    Some(p[9].parse().unwrap())
                };
                let table = if p.len() > 10 {
                    Some(cell_bits(p[10]))
                } else {
                    None
                };
                pairs.push(Pair {
                    genotype: p[2].to_string(),
                    qname: p[3].to_string(),
                    flags: p[4].parse().unwrap(),
                    a: cell_bits(p[5]),
                    b: cell_bits(p[6]),
                    java: cell_bits(p[7]),
                    branch: p[8].to_string(),
                    index,
                    table,
                });
            }
            "threshold" | "threshold_swapped" => {
                let index = if p[6] == "NA" {
                    None
                } else {
                    Some(p[6].parse::<i32>().unwrap())
                };
                thresholds.push((
                    p[1].to_string(),
                    p[2].to_string(),
                    cell_bits(p[3]),
                    cell_bits(p[4]),
                    p[5].to_string(),
                    index,
                ));
            }
            _ => {}
        }
    }
    assert_eq!(cache_len, 80_001);
    assert_eq!(pairs.len(), 122 * 3);
    let mut xor = 0u64;
    for value in literal_cache() {
        xor ^= value.to_bits();
    }
    let table_xor_matches = xor == cache_xor;
    let mut sample_table_mismatch = 0usize;
    for (k, java) in &sample_k {
        let rust = literal_cache()[*k];
        if rust.to_bits() != java.to_bits() {
            sample_table_mismatch += 1;
            kv(
                "table_sample_mismatch",
                format!(
                    "k={k} java={} rust={} delta={}",
                    fmt_f(*java),
                    fmt_f(rust),
                    fmt_f(rust - *java)
                ),
            );
        }
    }

    let mut b_vs_c = 0usize;
    let mut c_vs_java = 0usize;
    let mut table_index_mismatch = 0usize;
    let mut used_table_mismatch = 0usize;
    let mut by_gt = [
        ("T/TTTG", empty_stat()),
        ("T/TGTTTG", empty_stat()),
        ("TTTG/TGTTTG", empty_stat()),
    ];
    for pair in &pairs {
        let out_a = production_a(pair.a, pair.b);
        let out_b = approximate_log10_sum_log10_pair(pair.a, pair.b);
        let out_c = literal_c(pair.a, pair.b);
        if out_b.to_bits() != out_c.to_bits() {
            b_vs_c += 1;
        }
        if out_c.to_bits() != pair.java.to_bits() {
            c_vs_java += 1;
        }
        if let Some(index) = pair.index {
            let smaller = if pair.a > pair.b { pair.b } else { pair.a };
            let larger = if pair.a > pair.b { pair.a } else { pair.b };
            let diff = larger - smaller;
            let rust_index = literal_fast_round(diff * 10_000.0);
            if rust_index != index {
                table_index_mismatch += 1;
            }
            if let Some(table) = pair.table {
                if literal_cache()[index as usize].to_bits() != table.to_bits() {
                    used_table_mismatch += 1;
                }
            }
        }
        let stat = by_gt
            .iter_mut()
            .find(|(name, _)| *name == pair.genotype)
            .unwrap();
        let d = out_a - out_c;
        stat.1.n += 1;
        stat.1.sum += d;
        if out_a.to_bits() != out_c.to_bits() {
            stat.1.n_diff_ac += 1;
        }
        if out_b.to_bits() != out_c.to_bits() {
            stat.1.n_diff_bc += 1;
        }
        if d.abs() > stat.1.max_abs {
            let smaller = if pair.a > pair.b { pair.b } else { pair.a };
            let larger = if pair.a > pair.b { pair.a } else { pair.b };
            stat.1.max_abs = d.abs();
            stat.1.qname = pair.qname.clone();
            stat.1.flags = pair.flags;
            stat.1.a = pair.a;
            stat.1.b = pair.b;
            stat.1.out_a = out_a;
            stat.1.out_b = out_b;
            stat.1.out_c = out_c;
            stat.1.diff_inputs = larger - smaller;
        }
    }
    assert_eq!(b_vs_c, 0);
    assert_eq!(c_vs_java, 0);
    assert_eq!(table_index_mismatch, 0);

    let mut threshold_same_branch = 0usize;
    let mut threshold_b_eq_c = 0usize;
    let mut threshold_a_eq_c = 0usize;
    let mut threshold_n = 0usize;
    for (kind, label, diff, java, branch, index) in &thresholds {
        let (a, b) = if kind == "threshold" {
            (0.0, *diff)
        } else {
            (*diff, 0.0)
        };
        let out_a = production_a(a, b);
        let out_b = approximate_log10_sum_log10_pair(a, b);
        let out_c = literal_c(a, b);
        threshold_n += 1;
        if out_b.to_bits() == out_c.to_bits() && out_c.to_bits() == java.to_bits() {
            threshold_b_eq_c += 1;
        }
        if out_a.to_bits() == out_c.to_bits() {
            threshold_a_eq_c += 1;
        }
        let smaller = if a > b { b } else { a };
        let larger = if a > b { a } else { b };
        let rust_branch = if smaller == f64::NEG_INFINITY {
            "neg_inf"
        } else if larger - smaller < 8.0 {
            "table"
        } else {
            "ignore"
        };
        if rust_branch == branch {
            threshold_same_branch += 1;
        }
        if let Some(index) = index {
            assert_eq!(literal_fast_round((larger - smaller) * 10_000.0), *index);
        }
        kv(
            "threshold_row",
            format!(
                "{kind} diff={label} java_branch={branch} literal_branch={rust_branch} A={} B={} C={} java={}",
                fmt_f(out_a),
                fmt_f(out_b),
                fmt_f(out_c),
                fmt_f(*java)
            ),
        );
    }
    assert_eq!(threshold_same_branch, threshold_n);
    assert_eq!(threshold_b_eq_c, threshold_n);

    let mut gl_a = [0.0; 6];
    let mut gl_c = [0.0; 6];
    let log10_2 = 2.0_f64.log10();
    let mut n_reads = 0usize;
    for line in std::fs::read_to_string(repo_root().join(ALLELE_TSV))
        .unwrap()
        .lines()
    {
        let p: Vec<_> = line.split('\t').collect();
        if p.len() < 7 || p[0] != "6R297" || p[1] != "allele_ll" {
            continue;
        }
        let alleles = [cell_bits(p[4]), cell_bits(p[5]), cell_bits(p[6])];
        n_reads += 1;
        let hom = |x: f64| x + log10_2;
        gl_a[0] += hom(alleles[0]);
        gl_c[0] += hom(alleles[0]);
        gl_a[2] += hom(alleles[1]);
        gl_c[2] += hom(alleles[1]);
        gl_a[5] += hom(alleles[2]);
        gl_c[5] += hom(alleles[2]);
        gl_a[1] += production_a(alleles[0], alleles[1]);
        gl_c[1] += literal_c(alleles[0], alleles[1]);
        gl_a[3] += production_a(alleles[0], alleles[2]);
        gl_c[3] += literal_c(alleles[0], alleles[2]);
        gl_a[4] += production_a(alleles[1], alleles[2]);
        gl_c[4] += literal_c(alleles[1], alleles[2]);
    }
    assert_eq!(n_reads, 122);
    let den = n_reads as f64 * log10_2;
    for g in 0..6 {
        gl_a[g] -= den;
        gl_c[g] -= den;
    }
    assert_eq!(gl_a[0].to_bits(), gl_c[0].to_bits());
    assert_eq!(gl_a[2].to_bits(), gl_c[2].to_bits());
    assert_eq!(gl_a[5].to_bits(), gl_c[5].to_bits());

    let java_gl = java_vector("java_gl_vector");
    let java_pl = java_vector("java_continuous_pl_vector");
    for g in 0..6 {
        assert_eq!(gl_c[g].to_bits(), java_gl[g].to_bits());
    }
    let pl = |gl: &[f64; 6]| {
        let best = gl.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let mut out = [0.0; 6];
        for i in 0..6 {
            out[i] = -10.0 * (gl[i] - best);
        }
        out
    };
    let pl_a = pl(&gl_a);
    let pl_c = pl(&gl_c);
    for g in 0..6 {
        assert_eq!(pl_c[g].to_bits(), java_pl[g].to_bits());
    }

    let max = pairs
        .iter()
        .find(|p| {
            p.genotype == "T/TTTG"
                && p.qname == "HWI-D00360:8:H88U0ADXX:1:2101:19764:71575"
                && p.flags == 83
        })
        .unwrap();
    let swapped_a = production_a(max.b, max.a);
    let forward_a = production_a(max.a, max.b);
    kv(
        "argument_order",
        format!(
            "production selects max by strict >; forward={} swapped={} bits_equal={}",
            fmt_f(forward_a),
            fmt_f(swapped_a),
            forward_a.to_bits() == swapped_a.to_bits()
        ),
    );
    kv("helper_inv_step", fmt_f(1.0 / 0.0001));
    for (name, stat) in &by_gt {
        kv(
            &format!("stat_{}", name.replace('/', "_")),
            format!(
                "max_abs={} sum={} n_diff_ac={} n_diff_bc={} n={} read={} flags={} a={} b={} input_diff={} A={} B={} C={}",
                fmt_f(stat.max_abs),
                fmt_f(stat.sum),
                stat.n_diff_ac,
                stat.n_diff_bc,
                stat.n,
                stat.qname,
                stat.flags,
                fmt_f(stat.a),
                fmt_f(stat.b),
                fmt_f(stat.diff_inputs),
                fmt_f(stat.out_a),
                fmt_f(stat.out_b),
                fmt_f(stat.out_c)
            ),
        );
    }
    kv("table_xor_matches", table_xor_matches.to_string());
    kv("sample_table_mismatch", sample_table_mismatch.to_string());
    kv("used_table_mismatch", used_table_mismatch.to_string());
    kv("b_minus_c_rows", b_vs_c.to_string());
    kv("c_minus_java_rows", c_vs_java.to_string());
    kv(
        "threshold_a_eq_c",
        format!("{threshold_a_eq_c}/{threshold_n}"),
    );
    kv("java_gl", java_gl.map(fmt_f).join(","));
    kv("baseline_gl", gl_a.map(fmt_f).join(","));
    kv("cf_gl", gl_c.map(fmt_f).join(","));
    kv("java_pl", java_pl.map(fmt_f).join(","));
    kv("baseline_pl", pl_a.map(fmt_f).join(","));
    kv("cf_pl", pl_c.map(fmt_f).join(","));
    kv("hom_alt_java", fmt_f(java_pl[5]));
    kv("hom_alt_baseline", fmt_f(pl_a[5]));
    kv("hom_alt_cf", fmt_f(pl_c[5]));
    kv(
        "classification",
        "JAVA_JACOBIAN_TABLE_VS_RUST_ANALYTIC_COMBINE",
    );
    kv("production_change", "NONE");
}

fn java_vector(key: &str) -> [f64; 6] {
    let text = std::fs::read_to_string(repo_root().join(ALLELE_TSV)).unwrap();
    for line in text.lines() {
        let p: Vec<_> = line.split('\t').collect();
        if p.len() >= 8 && p[0] == "6R297" && p[1] == key {
            return [
                cell_bits(p[2]),
                cell_bits(p[3]),
                cell_bits(p[4]),
                cell_bits(p[5]),
                cell_bits(p[6]),
                cell_bits(p[7]),
            ];
        }
    }
    panic!("missing {key}");
}
