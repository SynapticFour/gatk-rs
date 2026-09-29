//! 6R.296: provenance of the 6R.295 matrices. No new numerical counterfactual.
//!
//! The Java file used by 6R.295 is `6r285_java_genotyping.tsv`.
//! `6r295_java_matrix.tsv` is a different kernel and is not a source.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r296_pairhmm_capture_provenance -- --nocapture --test-threads=1
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
const JAVA_TSV: &str = "gatk-haplotypecaller/tests/6r285_java_genotyping.tsv";
const LOCUS: u64 = 29_455_649;
const JAVA_SPAN: (u64, u64) = (29_455_560, 29_455_728);
const FLOOR_DELTA: f64 = -4.5;
const POOLS: [(usize, usize); 3] = [(0, 14), (14, 19), (19, 24)];
const TARGET_QNAME: &str = "HWI-D00360:5:H814YADXX:2:2207:19511:63503";
const TARGET_FLAGS: u16 = 163;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R296\t{key}\t{}", value.as_ref());
}

fn fmt_f(x: f64) -> String {
    format!("{x:.17}")
}

fn cell_bits(cell: &str) -> f64 {
    let hex = cell.split("bits=").nth(1).unwrap().trim();
    f64::from_bits(u64::from_str_radix(hex.trim_start_matches("0x"), 16).unwrap())
}

struct JavaRow {
    qname: String,
    flags: u16,
    values: Vec<f64>,
}

fn java_rows(text: &str, kind: &str) -> Vec<JavaRow> {
    let mut out = Vec::new();
    for line in text.lines() {
        let mut p = line.split('\t');
        if p.next() != Some("6R285") || p.next() != Some(kind) {
            continue;
        }
        let qname = p.next().unwrap().to_string();
        let flags: u16 = p.next().unwrap().parse().unwrap();
        out.push(JavaRow {
            qname,
            flags,
            values: p.map(cell_bits).collect(),
        });
    }
    out
}

fn field(text: &str, key: &str) -> String {
    for line in text.lines() {
        let mut p = line.split('\t');
        if p.next() == Some("6R285") && p.next() == Some(key) {
            return p.collect::<Vec<_>>().join("\t");
        }
    }
    panic!("missing {key}");
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

fn pool_max(row: &[f64]) -> [f64; 3] {
    std::array::from_fn(|a| {
        let (lo, hi) = POOLS[a];
        row[lo..hi]
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max)
    })
}

#[test]
fn forensic_6r296_pairhmm_capture_provenance() {
    let text = std::fs::read_to_string(repo_root().join(JAVA_TSV)).expect("java tsv");
    assert!(
        !text.contains("6r295_java_matrix"),
        "the AVX capture must not point at the LOGLESS dump"
    );
    let frozen = field(&text, "frozen_hap_index");
    let n_haps = field(&text, "n_haplotypes");
    let allele_columns = field(&text, "allele_columns");
    let n_evidence = field(&text, "n_evidence");
    let hap_ll = java_rows(&text, "hap_ll");
    let allele_ll = java_rows(&text, "allele_ll");
    assert_eq!(frozen, "[19, 20, 21, 22, 23]");
    assert_eq!(n_haps, "24");
    assert_eq!(allele_columns, "T,TTTG,TGTTTG");
    assert_eq!(n_evidence, "122");
    assert_eq!(allele_ll.len(), 122);
    assert!(hap_ll.iter().all(|r| r.values.len() == 5));
    assert!(allele_ll.iter().all(|r| r.values.len() == 3));
    let java_hap_keys: HashMap<_, _> = hap_ll
        .iter()
        .map(|r| ((r.qname.clone(), r.flags), r))
        .collect();
    let allele_with_hap = allele_ll
        .iter()
        .filter(|r| java_hap_keys.contains_key(&(r.qname.clone(), r.flags)))
        .count();
    assert_eq!(allele_with_hap, 122);

    gatk_haplotypecaller::hc_genotyping_engine::set_forensic_6r287_global_hap_floor(None, None);
    gatk_haplotypecaller::hc_genotyping_engine::set_forensic_6r288_haplotype_substitute(
        None, None, 0, None,
    );
    gatk_haplotypecaller::likelihood_engine::set_forensic_6r289_target(None, 0);
    gatk_haplotypecaller::likelihood_engine::set_forensic_6r289_substitute(None, None);
    set_forensic_6r294_padded_span(Some(JAVA_SPAN.0), Some(JAVA_SPAN.1));
    begin_likelihood_pipeline_observe();
    let root = repo_root();
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
    let mut post_kernel_reads = 0usize;
    for cell in cells.iter().filter(|c| c.stage == "post_kernel") {
        let row = raw
            .entry((cell.qname.clone(), cell.flags))
            .or_insert_with(|| {
                post_kernel_reads += 1;
                vec![f64::NAN; 24]
            });
        if cell.hap_index < 24 && row[cell.hap_index].is_nan() {
            row[cell.hap_index] = cell.log10_likelihood;
        }
    }
    assert!(raw.values().all(|row| row.iter().all(|v| v.is_finite())));
    assert!(raw.values().all(|row| row.len() == 24));

    let mut n_unfloored = 0usize;
    let mut n_post_norm = 0usize;
    let mut allele_sum = 0.0;
    let mut allele_cells = 0usize;
    for j in &allele_ll {
        let raw_row = raw.get(&(j.qname.clone(), j.flags)).expect("rust row");
        let floored = floor_row(raw_row);
        let java_hap = java_hap_keys.get(&(j.qname.clone(), j.flags)).unwrap();
        for k in 0..5 {
            n_post_norm += 1;
            if raw_row[19 + k] >= floored[19 + k] {
                n_unfloored += 1;
                let _ = java_hap.values[k];
            }
        }
        let rust_al = pool_max(&floored);
        for a in 0..3 {
            allele_sum += rust_al[a] - j.values[a];
            allele_cells += 1;
        }
    }
    assert_eq!(n_post_norm, 610);
    assert_eq!(n_unfloored, 73);
    assert_eq!(allele_cells, 366);
    assert!((allele_sum - 0.00005005406882219).abs() < 1e-16);

    let raw_row = raw
        .get(&(TARGET_QNAME.to_string(), TARGET_FLAGS))
        .expect("target rust row");
    let floored = floor_row(raw_row);
    let rust_argmax = (0..14)
        .max_by(|&i, &j| floored[i].partial_cmp(&floored[j]).unwrap())
        .unwrap();
    let java_allele = allele_ll
        .iter()
        .find(|r| r.qname == TARGET_QNAME && r.flags == TARGET_FLAGS)
        .expect("target allele");
    let java_hap = java_hap_keys
        .get(&(TARGET_QNAME.to_string(), TARGET_FLAGS))
        .expect("target hap");
    assert_eq!(rust_argmax, 5);
    assert!((floored[5] - java_allele.values[0]).abs() > 0.0);

    kv(
        "rust_source",
        "capture_scored_likelihood_pipeline stage=post_kernel",
    );
    kv(
        "rust_stage",
        "raw PairHMM return, before post_process_pairhmm_likelihoods",
    );
    kv(
        "rust_coverage",
        format!("{post_kernel_reads} reads x 24 haplotypes"),
    );
    kv("java_hap_source", "6r285_java_genotyping.tsv hap_ll");
    kv(
        "java_hap_producer",
        "Downstream285.dumpFrozenColumns at assignGenotypeLikelihoods",
    );
    kv(
        "java_hap_coverage",
        format!(
            "{} reads x 5 stored columns; frozen_hap_index {frozen}; n_haplotypes {n_haps}",
            hap_ll.len()
        ),
    );
    kv("java_haps_0_18", "NOT_CAPTURED");
    kv(
        "unfloored_73",
        "122 allele reads x columns 19-23, kept only where Rust raw was not raised by the test floor; Java side is post-normalize hap_ll",
    );
    kv(
        "post_norm_610",
        "122 allele reads x columns 19-23; Rust side is the test floor of post_kernel; Java side is post-normalize hap_ll",
    );
    kv(
        "target_rust_raw_0_13",
        raw_row[..14]
            .iter()
            .map(|v| fmt_f(*v))
            .collect::<Vec<_>>()
            .join(","),
    );
    kv(
        "target_rust_floored_0_13",
        floored[..14]
            .iter()
            .map(|v| fmt_f(*v))
            .collect::<Vec<_>>()
            .join(","),
    );
    kv("target_rust_T_argmax", rust_argmax.to_string());
    kv("target_rust_T_max", fmt_f(floored[rust_argmax]));
    kv("target_java_haps_0_13", "NOT_CAPTURED");
    kv("target_java_T_argmax", "NOT_CAPTURED");
    kv("target_java_allele_T", fmt_f(java_allele.values[0]));
    kv(
        "target_java_hap_19_23",
        java_hap
            .values
            .iter()
            .map(|v| fmt_f(*v))
            .collect::<Vec<_>>()
            .join(","),
    );
    kv("allele_direct_cells", allele_cells.to_string());
    kv("allele_sum_recomputed", fmt_f(allele_sum));
    kv(
        "allele_sum_source",
        "Rust test pool-max of floored post_kernel minus 6r285 allele_ll; 366 cells; Rust side reconstructed, Java side stored",
    );
    kv("complete_matrix_comparison_valid", "NO");
    kv("classification", "CAPTURE_SCOPE_EVIDENCE_GAP");
    kv("production_change", "NONE");
}
