//! 6R.300: per-read genotype-likelihood contributions from the canonical
//! LOGLESS allele matrix. Production code is unchanged.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r300_genotype_likelihood_localization -- --nocapture --test-threads=1
//! ```

use gatk_haplotypecaller::{approximate_log10_sum_log10_pair, log10_sum_log10};
use std::path::{Path, PathBuf};

const ALLELE_TSV: &str = "gatk-haplotypecaller/tests/6r297_java_hap_ll.tsv";
const JAVA_READ_TSV: &str = "gatk-haplotypecaller/tests/6r300_per_read_gl.tsv";

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R300\t{key}\t{}", value.as_ref());
}

fn fmt_f(x: f64) -> String {
    format!("{x:.17}")
}

fn cell_bits(cell: &str) -> f64 {
    let hex = cell.split("bits=").nth(1).unwrap().trim();
    f64::from_bits(u64::from_str_radix(hex.trim_start_matches("0x"), 16).unwrap())
}

struct ReadGl {
    qname: String,
    flags: u16,
    alleles: [f64; 3],
    java: [f64; 6],
}

fn load() -> Vec<ReadGl> {
    let root = repo_root();
    let alleles = std::fs::read_to_string(root.join(ALLELE_TSV)).unwrap();
    let java = std::fs::read_to_string(root.join(JAVA_READ_TSV)).unwrap();
    let mut rows = Vec::new();
    for line in alleles.lines() {
        let p: Vec<_> = line.split('\t').collect();
        if p.len() < 7 || p[0] != "6R297" || p[1] != "allele_ll" {
            continue;
        }
        rows.push(ReadGl {
            qname: p[2].to_string(),
            flags: p[3].parse().unwrap(),
            alleles: [cell_bits(p[4]), cell_bits(p[5]), cell_bits(p[6])],
            java: [0.0; 6],
        });
    }
    let mut j = 0usize;
    for line in java.lines() {
        let p: Vec<_> = line.split('\t').collect();
        if p.len() < 9 || p[0] != "6R300" || p[1] != "read_gl" {
            continue;
        }
        assert_eq!(p[2], rows[j].qname);
        assert_eq!(p[3].parse::<u16>().unwrap(), rows[j].flags);
        for g in 0..6 {
            rows[j].java[g] = cell_bits(p[4 + g]);
        }
        j += 1;
    }
    assert_eq!(j, rows.len());
    rows
}

fn rust_exact(alleles: [f64; 3]) -> [f64; 6] {
    let log10_2 = 2.0_f64.log10();
    let het = |a: f64, b: f64| log10_sum_log10(&[a, b]);
    [
        alleles[0] + log10_2,
        het(alleles[0], alleles[1]),
        alleles[1] + log10_2,
        het(alleles[0], alleles[2]),
        het(alleles[1], alleles[2]),
        alleles[2] + log10_2,
    ]
}

fn rust_approx(alleles: [f64; 3]) -> [f64; 6] {
    let log10_2 = 2.0_f64.log10();
    let het = |a: f64, b: f64| approximate_log10_sum_log10_pair(a, b);
    [
        alleles[0] + log10_2,
        het(alleles[0], alleles[1]),
        alleles[1] + log10_2,
        het(alleles[0], alleles[2]),
        het(alleles[1], alleles[2]),
        alleles[2] + log10_2,
    ]
}

struct MaxDelta {
    abs: f64,
    signed: f64,
    qname: String,
    flags: u16,
    java: f64,
    rust: f64,
    n: usize,
}

fn max_delta(rows: &[ReadGl], g: usize, rust_of: impl Fn([f64; 3]) -> [f64; 6]) -> MaxDelta {
    let mut best = MaxDelta {
        abs: -1.0,
        signed: 0.0,
        qname: String::new(),
        flags: 0,
        java: 0.0,
        rust: 0.0,
        n: 0,
    };
    for row in rows {
        let rust = rust_of(row.alleles)[g];
        let d = rust - row.java[g];
        best.n += 1;
        if d.abs() > best.abs {
            best.abs = d.abs();
            best.signed = d;
            best.qname = row.qname.clone();
            best.flags = row.flags;
            best.java = row.java[g];
            best.rust = rust;
        }
    }
    best
}

#[test]
fn forensic_6r300_genotype_likelihood_localization() {
    let rows = load();
    assert_eq!(rows.len(), 122);
    let names = [
        "T/T",
        "T/TTTG",
        "TTTG/TTTG",
        "T/TGTTTG",
        "TTTG/TGTTTG",
        "TGTTTG/TGTTTG",
    ];
    let mut exact = Vec::new();
    let mut approx_mismatch = 0usize;
    for g in 0..6 {
        exact.push(max_delta(&rows, g, rust_exact));
        for row in &rows {
            if rust_approx(row.alleles)[g].to_bits() != row.java[g].to_bits() {
                approx_mismatch += 1;
            }
        }
    }
    assert_eq!(exact[0].abs, 0.0);
    assert_eq!(exact[2].abs, 0.0);
    assert_eq!(exact[5].abs, 0.0);
    assert!(exact[1].abs > exact[0].abs);
    assert_eq!(approx_mismatch, 0);
    let classification = "GENOTYPE_LIKELIHOOD_PER_READ_OPERATION_DIVERGENCE";
    kv("java_path", "HaplotypeCallerGenotypingEngine.calculateGLsForThisEvent:544 -> IndependentSampleGenotypesModel.calculateLikelihoods:66 -> GenotypeLikelihoodCalculator.genotypeLikelihoods:233 -> twoComponentGenotypeLikelihoodByRead:374");
    kv("rust_path", "hc_genotyping_engine/mod.rs diploid_genotype_log10_likelihoods_from_allele_rows:469-472 -> activity_scoring.rs log10_sum_log10:270");
    kv("java_operation", "MathUtils.approximateLog10SumLog10(double,double) line 492, JacobianLogTable for heterozygotes; homozygote is allele likelihood + log10(2)");
    kv("rust_operation", "exact log10(10^a + 10^b) via log10_sum_log10 for heterozygotes; homozygote is allele likelihood + log10(2)");
    kv(
        "genotype_order",
        "0/0,0/1,1/1,0/2,1/2,2/2 = T/T,T/TTTG,TTTG/TTTG,T/TGTTTG,TTTG/TGTTTG,TGTTTG/TGTTTG",
    );
    kv("per_read_dimensions", "122 reads x 6 genotypes");
    for (g, stat) in exact.iter().enumerate() {
        kv(
            &format!("max_delta_{}", names[g].replace('/', "_")),
            format!(
                "abs={} signed={} read={} flags={} java={} rust={}",
                fmt_f(stat.abs),
                fmt_f(stat.signed),
                stat.qname,
                stat.flags,
                fmt_f(stat.java),
                fmt_f(stat.rust)
            ),
        );
        kv(&format!("n_reads_{g}"), stat.n.to_string());
        let mut sum = 0.0;
        let mut nonzero = 0usize;
        for row in &rows {
            let d = rust_exact(row.alleles)[g] - row.java[g];
            if d != 0.0 {
                nonzero += 1;
            }
            sum += d;
        }
        kv(&format!("sum_delta_{g}"), fmt_f(sum));
        kv(&format!("nonzero_reads_{g}"), nonzero.to_string());
    }
    kv("approx_pair_mismatches", approx_mismatch.to_string());
    kv(
        "first_operation",
        "heterozygote per-read combine: log10_sum_log10 versus MathUtils.approximateLog10SumLog10",
    );
    kv("classification", classification);
    kv("production_change", "NONE");
}
