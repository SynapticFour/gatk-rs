//! 6R.303: production heterozygote combine uses the Java Jacobian primitive.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r303_jacobian_production_fix -- --nocapture --test-threads=1
//! ```

use gatk_haplotypecaller::hc_genotyping_engine::{
    biallelic_genotype_log10_likelihoods_gatk, diploid_genotype_log10_likelihoods_from_allele_rows,
};
use gatk_haplotypecaller::{
    approximate_log10_sum_log10_pair, emit_genotype_format_fields, log10_sum_log10,
    ReadLikelihoodRow,
};
use std::path::{Path, PathBuf};

const ALLELE_TSV: &str = "gatk-haplotypecaller/tests/6r297_java_hap_ll.tsv";

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R303\t{key}\t{}", value.as_ref());
}

fn fmt_f(x: f64) -> String {
    format!("{x:.17}")
}

fn cell_bits(cell: &str) -> f64 {
    let hex = cell.split("bits=").nth(1).unwrap().trim();
    f64::from_bits(u64::from_str_radix(hex.trim_start_matches("0x"), 16).unwrap())
}

fn row(values: [f64; 2]) -> ReadLikelihoodRow {
    ReadLikelihoodRow {
        read_index: 0,
        read_id: String::new(),
        haplotype_log10_likelihoods: values.to_vec(),
    }
}

#[test]
fn forensic_6r303_jacobian_production_fix() {
    let log10_2 = 2.0_f64.log10();
    assert_eq!(
        approximate_log10_sum_log10_pair(-3.25, -3.25).to_bits(),
        (-3.25 + log10_2).to_bits()
    );
    let small = approximate_log10_sum_log10_pair(0.0, 1e-9);
    assert_ne!(small.to_bits(), log10_sum_log10(&[0.0, 1e-9]).to_bits());
    let below = f64::from_bits(0x401fffffffffffff);
    let above = f64::from_bits(0x4020000000000001);
    assert_eq!(
        approximate_log10_sum_log10_pair(0.0, below).to_bits(),
        0x4020000000254e3c
    );
    assert_eq!(
        approximate_log10_sum_log10_pair(0.0, 8.0).to_bits(),
        8.0_f64.to_bits()
    );
    assert_eq!(
        approximate_log10_sum_log10_pair(0.0, above).to_bits(),
        above.to_bits()
    );
    assert_eq!(
        approximate_log10_sum_log10_pair(-2.0, -1.0).to_bits(),
        approximate_log10_sum_log10_pair(-1.0, -2.0).to_bits()
    );
    assert_eq!(
        approximate_log10_sum_log10_pair(f64::NEG_INFINITY, -1.5).to_bits(),
        (-1.5_f64).to_bits()
    );
    assert_eq!(
        approximate_log10_sum_log10_pair(-1.5, f64::NEG_INFINITY).to_bits(),
        (-1.5_f64).to_bits()
    );
    assert_eq!(
        approximate_log10_sum_log10_pair(f64::NEG_INFINITY, f64::NEG_INFINITY).to_bits(),
        f64::NEG_INFINITY.to_bits()
    );
    let analytic = log10_sum_log10(&[-1.0, -2.0]);
    let expected = -1.0 + (1.0 + 10_f64.powf(-2.0 - -1.0)).log10();
    assert_eq!(analytic.to_bits(), expected.to_bits());

    let one = biallelic_genotype_log10_likelihoods_gatk(&[row([-0.4, -1.7])], 0, 1);
    assert_eq!(one[0].to_bits(), (-0.4_f64).to_bits());
    assert_eq!(one[2].to_bits(), (-1.7_f64).to_bits());
    assert_eq!(
        one[1].to_bits(),
        (approximate_log10_sum_log10_pair(-0.4, -1.7) - log10_2).to_bits()
    );

    let text = std::fs::read_to_string(repo_root().join(ALLELE_TSV)).unwrap();
    let mut rows = Vec::new();
    let mut java_gl = [0.0; 6];
    let mut java_pl = [0.0; 6];
    for line in text.lines() {
        let p: Vec<_> = line.split('\t').collect();
        if p.len() >= 8 && p[0] == "6R297" && p[1] == "java_gl_vector" {
            java_gl = [2, 3, 4, 5, 6, 7].map(|i| cell_bits(p[i]));
        }
        if p.len() >= 8 && p[0] == "6R297" && p[1] == "java_continuous_pl_vector" {
            java_pl = [2, 3, 4, 5, 6, 7].map(|i| cell_bits(p[i]));
        }
        if p.len() < 7 || p[0] != "6R297" || p[1] != "allele_ll" {
            continue;
        }
        let alleles = [cell_bits(p[4]), cell_bits(p[5]), cell_bits(p[6])];
        assert!(alleles.iter().all(|v| v.is_finite()));
        rows.push(ReadLikelihoodRow {
            read_index: rows.len(),
            read_id: p[2].to_string(),
            haplotype_log10_likelihoods: alleles.to_vec(),
        });
    }
    assert_eq!(rows.len(), 122);
    let gl = diploid_genotype_log10_likelihoods_from_allele_rows(&rows, 3);
    assert_eq!(gl.len(), 6);
    for i in 0..6 {
        assert_eq!(gl[i].to_bits(), java_gl[i].to_bits());
    }
    let best = gl.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let mut pl = [0.0; 6];
    for i in 0..6 {
        pl[i] = -10.0 * (gl[i] - best);
        assert_eq!(pl[i].to_bits(), java_pl[i].to_bits());
    }
    let format = emit_genotype_format_fields(&gl, &[1, 1, 1]).unwrap();
    let ipl = format.pl_as_i32();
    assert_eq!(ipl, vec![570, 148, 3485, 0, 2762, 3517]);
    kv(
        "homozygote_gl",
        format!("{},{},{}", fmt_f(gl[0]), fmt_f(gl[2]), fmt_f(gl[5])),
    );
    kv(
        "heterozygote_gl",
        format!("{},{},{}", fmt_f(gl[1]), fmt_f(gl[3]), fmt_f(gl[4])),
    );
    kv(
        "final_gl",
        gl.iter().map(|v| fmt_f(*v)).collect::<Vec<_>>().join(","),
    );
    kv(
        "final_pl",
        pl.iter().map(|v| fmt_f(*v)).collect::<Vec<_>>().join(","),
    );
    kv(
        "integer_pl",
        ipl.iter()
            .map(|v| v.to_string())
            .collect::<Vec<_>>()
            .join(","),
    );
    kv(
        "emitted_subset_pl",
        format!("{},{},{}", ipl[0], ipl[3], ipl[5]),
    );
    kv(
        "classification",
        "JAVA_COMPATIBLE_JACOBIAN_PRODUCTION_FIX_VALIDATED",
    );
}
