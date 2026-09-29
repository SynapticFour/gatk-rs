//! 6R.297: direct post-normalization comparison of all 24 haplotype
//! columns. Java columns come from `6r297_java_hap_ll.tsv`, captured at
//! `assignGenotypeLikelihoods` before marginalization.
//!
//! The host cannot execute AVX. This Java matrix is LOGLESS_CACHING.
//! The five-column AVX file `6r285_java_genotyping.tsv` is unchanged
//! and is used only to name the 122 genotyping reads and to check
//! whether columns 19–23 reproduce the 6R.295 AVX comparison.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r297_complete_java_hap_ll_capture -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::{
    begin_likelihood_pipeline_observe, call_disposition, flatten_assembly_regions,
    set_forensic_6r294_padded_span, take_likelihood_pipeline_cells,
    traverse_assembly_region_walker, AssemblyRegionCallDisposition, CallRegionArgs,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
};
use std::collections::HashMap;
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const JAVA297: &str = "gatk-haplotypecaller/tests/6r297_java_hap_ll.tsv";
const JAVA285: &str = "gatk-haplotypecaller/tests/6r285_java_genotyping.tsv";
const LOCUS: u64 = 29_455_649;
const JAVA_SPAN: (u64, u64) = (29_455_560, 29_455_728);
const FLOOR_DELTA: f64 = -4.5;
const TARGET_QNAME: &str = "HWI-D00360:5:H814YADXX:2:2207:19511:63503";
const TARGET_FLAGS: u16 = 163;
const AVX_ALLELE_T: f64 = f64::from_bits(0xc0029ba800000000);

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R297\t{key}\t{}", value.as_ref());
}

fn fmt_f(x: f64) -> String {
    format!("{x:.17}")
}

fn cell_bits(cell: &str) -> f64 {
    let hex = cell.split("bits=").nth(1).unwrap().trim();
    f64::from_bits(u64::from_str_radix(hex.trim_start_matches("0x"), 16).unwrap())
}

struct JavaHap {
    evidence_index: usize,
    values: Vec<f64>,
}

fn field(text: &str, tag: &str, key: &str) -> String {
    for line in text.lines() {
        let mut p = line.split('\t');
        if p.next() == Some(tag) && p.next() == Some(key) {
            return p.collect::<Vec<_>>().join("\t");
        }
    }
    panic!("missing {tag} {key}");
}

fn java297_haps(text: &str) -> HashMap<(String, u16), JavaHap> {
    let mut out = HashMap::new();
    for line in text.lines() {
        let mut p = line.split('\t');
        if p.next() != Some("6R297") || p.next() != Some("hap_ll") {
            continue;
        }
        let qname = p.next().unwrap().to_string();
        let flags: u16 = p.next().unwrap().parse().unwrap();
        let evidence_index: usize = p.next().unwrap().parse().unwrap();
        let values: Vec<f64> = p.map(cell_bits).collect();
        out.insert(
            (qname, flags),
            JavaHap {
                evidence_index,
                values,
            },
        );
    }
    out
}

fn java_allele_t(text: &str, tag: &str) -> HashMap<(String, u16), f64> {
    let mut out = HashMap::new();
    for line in text.lines() {
        let mut p = line.split('\t');
        if p.next() != Some(tag) || p.next() != Some("allele_ll") {
            continue;
        }
        let qname = p.next().unwrap().to_string();
        let flags: u16 = p.next().unwrap().parse().unwrap();
        let t = cell_bits(p.next().unwrap());
        out.insert((qname, flags), t);
    }
    out
}

fn avx_frozen(text: &str) -> HashMap<(String, u16), Vec<f64>> {
    let mut out = HashMap::new();
    for line in text.lines() {
        let mut p = line.split('\t');
        if p.next() != Some("6R285") || p.next() != Some("hap_ll") {
            continue;
        }
        let qname = p.next().unwrap().to_string();
        let flags: u16 = p.next().unwrap().parse().unwrap();
        out.insert((qname, flags), p.map(cell_bits).collect());
    }
    out
}

fn floor_row(raw: &[f64]) -> Vec<f64> {
    let best = raw.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let floor = best + FLOOR_DELTA;
    raw.iter()
        .map(|v| {
            if v.is_finite() && *v < floor {
                floor
            } else {
                *v
            }
        })
        .collect()
}

struct Pool {
    n: usize,
    n_diff: usize,
    max_abs: f64,
    where_max: String,
    sum: f64,
    sum_abs: f64,
}

fn pool_stats(pairs: impl IntoIterator<Item = (f64, f64, String)>) -> Pool {
    let mut s = Pool {
        n: 0,
        n_diff: 0,
        max_abs: 0.0,
        where_max: String::new(),
        sum: 0.0,
        sum_abs: 0.0,
    };
    for (rust, java, label) in pairs {
        let d = rust - java;
        s.n += 1;
        if d != 0.0 {
            s.n_diff += 1;
        }
        s.sum += d;
        s.sum_abs += d.abs();
        if d.abs() > s.max_abs {
            s.max_abs = d.abs();
            s.where_max = format!("{label} delta={}", fmt_f(d));
        }
    }
    s
}

#[test]
fn forensic_6r297_complete_java_hap_ll_capture() {
    let root = repo_root();
    let text297 = std::fs::read_to_string(root.join(JAVA297)).expect("6r297 tsv");
    let text285 = std::fs::read_to_string(root.join(JAVA285)).expect("6r285 tsv");
    assert_eq!(field(&text297, "6R297", "pairhmm_impl"), "LOGLESS_CACHING");
    assert_eq!(
        field(&text297, "6R297", "matrix_stage"),
        "assignGenotypeLikelihoods_post_normalize"
    );
    assert_eq!(field(&text297, "6R297", "haplotype_count"), "24");
    assert_eq!(
        field(&text297, "6R297", "frozen_hap_index"),
        "[19, 20, 21, 22, 23]"
    );
    let java = java297_haps(&text297);
    let java_allele = java_allele_t(&text297, "6R297");
    let avx_allele = java_allele_t(&text285, "6R285");
    let avx_hap = avx_frozen(&text285);
    let allele_reads: Vec<(String, u16)> = {
        let mut v: Vec<_> = avx_allele.keys().cloned().collect();
        v.sort();
        v
    };
    assert_eq!(allele_reads.len(), 122);
    let missing: Vec<_> = allele_reads
        .iter()
        .filter(|k| !java.contains_key(*k))
        .map(|(q, f)| format!("{q} flags={f}"))
        .collect();
    assert!(
        missing.is_empty(),
        "missing Java rows: {}",
        missing.join(",")
    );
    assert!(java.values().all(|h| h.values.len() == 24));

    gatk_haplotypecaller::hc_genotyping_engine::set_forensic_6r287_global_hap_floor(None, None);
    gatk_haplotypecaller::hc_genotyping_engine::set_forensic_6r288_haplotype_substitute(
        None, None, 0, None,
    );
    gatk_haplotypecaller::likelihood_engine::set_forensic_6r289_target(None, 0);
    gatk_haplotypecaller::likelihood_engine::set_forensic_6r289_substitute(None, None);
    set_forensic_6r294_padded_span(Some(JAVA_SPAN.0), Some(JAVA_SPAN.1));
    begin_likelihood_pipeline_observe();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing chr20_tiny BAM/REF");
        return;
    }
    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, INTERVAL).expect("interval");
    let walk = traverse_assembly_region_walker(
        &dict,
        &specs,
        &ref_fasta,
        &bam,
        &ReadFilterParams::gatk_standard_hc(),
        &WalkerTraversalConfig::gatk_haplotype_caller_production(100),
    )
    .expect("walk");
    let regions = flatten_assembly_regions(&walk);
    let region = regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= LOCUS
                && r.end.get() >= LOCUS
        })
        .expect("ActiveFull");
    let _outcome = HaplotypeCallerEngine::call_region(
        region,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");
    set_forensic_6r294_padded_span(None, None);
    let cells = take_likelihood_pipeline_cells();
    let mut raw: HashMap<(String, u16), Vec<f64>> = HashMap::new();
    for cell in cells.iter().filter(|c| c.stage == "post_kernel") {
        let row = raw
            .entry((cell.qname.clone(), cell.flags))
            .or_insert_with(|| vec![f64::NAN; 24]);
        if cell.hap_index < 24 && row[cell.hap_index].is_nan() {
            row[cell.hap_index] = cell.log10_likelihood;
        }
    }

    let mut pairs_a = Vec::new();
    let mut pairs_b = Vec::new();
    let mut pairs_c = Vec::new();
    let mut pairs_avx = Vec::new();
    let mut matched = 0usize;
    for key in &allele_reads {
        let java_row = &java[key];
        let raw_row = raw
            .get(key)
            .unwrap_or_else(|| panic!("rust missing {}", key.0));
        assert!(raw_row.iter().all(|v| v.is_finite()));
        let floored = floor_row(raw_row);
        matched += 1;
        for h in 0..24 {
            let label = format!("read={} flags={} hap={h}", key.0, key.1);
            let pair = (floored[h], java_row.values[h], label);
            if h < 14 {
                pairs_a.push(pair);
            } else if h < 19 {
                pairs_b.push(pair);
            } else {
                pairs_c.push(pair);
            }
        }
        if let Some(avx) = avx_hap.get(key) {
            for k in 0..5 {
                pairs_avx.push((
                    floored[19 + k],
                    avx[k],
                    format!("read={} flags={} hap={}", key.0, key.1, 19 + k),
                ));
            }
        }
    }
    assert_eq!(matched, 122);
    let a = pool_stats(pairs_a);
    let b = pool_stats(pairs_b);
    let c = pool_stats(pairs_c);
    let all = Pool {
        n: a.n + b.n + c.n,
        n_diff: a.n_diff + b.n_diff + c.n_diff,
        max_abs: a.max_abs.max(b.max_abs).max(c.max_abs),
        where_max: if a.max_abs >= b.max_abs && a.max_abs >= c.max_abs {
            a.where_max.clone()
        } else if b.max_abs >= c.max_abs {
            b.where_max.clone()
        } else {
            c.where_max.clone()
        },
        sum: a.sum + b.sum + c.sum,
        sum_abs: a.sum_abs + b.sum_abs + c.sum_abs,
    };
    let avx_c = pool_stats(pairs_avx);
    assert_eq!(a.n, 122 * 14);
    assert_eq!(b.n, 122 * 5);
    assert_eq!(c.n, 122 * 5);
    assert_eq!(all.n, 2928);

    let java_row = &java[&(TARGET_QNAME.to_string(), TARGET_FLAGS)];
    let java_max = java_row.values[..14]
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let java_argmax: Vec<usize> = (0..14)
        .filter(|&h| java_row.values[h] == java_max)
        .collect();
    let same_run_allele = java_allele[&(TARGET_QNAME.to_string(), TARGET_FLAGS)];
    let avx_allele_t = avx_allele[&(TARGET_QNAME.to_string(), TARGET_FLAGS)];
    assert_eq!(avx_allele_t.to_bits(), AVX_ALLELE_T.to_bits());
    let allele_equals_hap_max = same_run_allele.to_bits() == java_max.to_bits();
    let avx_allele_equals_hap_max = avx_allele_t.to_bits() == java_max.to_bits();
    let c19_reproduces = avx_c.max_abs == c.max_abs && avx_c.sum.to_bits() == c.sum.to_bits();

    kv("checkpoint", "assignGenotypeLikelihoods_post_normalize");
    kv("pairhmm_impl", "LOGLESS_CACHING");
    kv("avx_native", "NOT_EXECUTED_HOST_ARM64");
    kv("java_evidence_rows", java.len().to_string());
    kv("rust_rows_matched", matched.to_string());
    kv("missing_reads", "NONE");
    kv("pool_0_13_cells", a.n.to_string());
    kv("pool_0_13_differing", a.n_diff.to_string());
    kv("pool_0_13_max_abs", fmt_f(a.max_abs));
    kv("pool_0_13_where", &a.where_max);
    kv("pool_0_13_sum", fmt_f(a.sum));
    kv("pool_0_13_sum_abs", fmt_f(a.sum_abs));
    kv("pool_14_18_cells", b.n.to_string());
    kv("pool_14_18_differing", b.n_diff.to_string());
    kv("pool_14_18_max_abs", fmt_f(b.max_abs));
    kv("pool_14_18_where", &b.where_max);
    kv("pool_14_18_sum", fmt_f(b.sum));
    kv("pool_14_18_sum_abs", fmt_f(b.sum_abs));
    kv("pool_19_23_cells", c.n.to_string());
    kv("pool_19_23_differing", c.n_diff.to_string());
    kv("pool_19_23_max_abs", fmt_f(c.max_abs));
    kv("pool_19_23_where", &c.where_max);
    kv("pool_19_23_sum", fmt_f(c.sum));
    kv("pool_19_23_sum_abs", fmt_f(c.sum_abs));
    kv("avx_19_23_max_abs", fmt_f(avx_c.max_abs));
    kv("avx_19_23_sum", fmt_f(avx_c.sum));
    kv(
        "columns_19_23_consistency",
        if c19_reproduces {
            "REPRODUCES_AVX_6R295"
        } else {
            "DOES_NOT_REPRODUCE_AVX_6R295"
        },
    );
    kv("matrix_cells", all.n.to_string());
    kv("matrix_differing", all.n_diff.to_string());
    kv("matrix_max_abs", fmt_f(all.max_abs));
    kv("matrix_where", &all.where_max);
    kv("matrix_sum", fmt_f(all.sum));
    kv("matrix_sum_abs", fmt_f(all.sum_abs));
    kv(
        "target_java_evidence_index",
        java_row.evidence_index.to_string(),
    );
    kv(
        "target_java_T_argmax",
        java_argmax
            .iter()
            .map(|h| h.to_string())
            .collect::<Vec<_>>()
            .join(","),
    );
    kv("target_java_T_max", fmt_f(java_max));
    kv("target_same_run_allele_T", fmt_f(same_run_allele));
    kv(
        "target_allele_equals_hap_max",
        if allele_equals_hap_max { "YES" } else { "NO" },
    );
    kv("target_avx_allele_T", fmt_f(avx_allele_t));
    kv(
        "target_avx_allele_equals_this_hap_max",
        if avx_allele_equals_hap_max {
            "YES"
        } else {
            "NO"
        },
    );
    assert_eq!(all.n_diff, 0);
    assert!(allele_equals_hap_max);
    assert!(!avx_allele_equals_hap_max);
    assert!(!c19_reproduces);
    assert_eq!(fmt_f(avx_c.max_abs), "0.00000238641348460");
    kv("first_divergence", "NONE_IN_CAPTURED_24_COLUMN_MATRIX");
    kv(
        "classification",
        "LOGLESS_24_COLUMN_MATRIX_IDENTICAL_AVX_0_18_UNCAPTURED",
    );
    kv("production_change", "NONE");
}
