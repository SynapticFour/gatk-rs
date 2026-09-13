//! 6R.177: InbreedingCoeff emit/suppress at `2:92305634 G/T` (proof of the arrow).
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! 6R.178 closed the emit gate (`n_genotypes >= 10`). This file still
//! documents Java `MIN_SAMPLES=10` and the remaining stacked formula split
//! (`1 - het/n` vs Java HWE `1 - het/(2pq n)`).
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r177_inbreeding_coeff_emit_predicate -- --nocapture --test-threads=1
//! HOLDOUT_6R177=1 cargo test -p gatk-haplotypecaller --test holdout_6r177_inbreeding -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

const INTERVAL: &str = "2:92305500-92305850";
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const JAVA_VCF_REL: &str = "parity/reports/6r43/p12_snp_cluster/java.vcf";
const TARGET: u64 = 92_305_634;
const MERGED_REF: &str = "G";
const MERGED_ALT: &str = "T";
const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
/// Java `InbreedingCoeff.MIN_SAMPLES` (SHA `2dbc0258`).
const JAVA_MIN_SAMPLES: usize = 10;
const LIVE_PL: [i32; 3] = [90, 6, 0];

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R177\t{key}\t{}", value.as_ref());
}

fn info_i32(info: &[InfoValue], key: &str) -> Option<i32> {
    for v in info {
        if let InfoValue::Integer(k, xs) = v {
            if k == key {
                return xs.first().copied();
            }
        }
    }
    None
}

fn info_f64(info: &[InfoValue], key: &str) -> Option<f64> {
    for v in info {
        if let InfoValue::Float(k, xs) = v {
            if k == key {
                return xs.first().copied();
            }
        }
    }
    None
}

fn info_has(info: &[InfoValue], key: &str) -> bool {
    info.iter().any(|v| match v {
        InfoValue::Flag(k)
        | InfoValue::Integer(k, _)
        | InfoValue::Float(k, _)
        | InfoValue::String(k, _)
        | InfoValue::Character(k, _) => k == key,
    })
}

fn parse_info_map(info: &str) -> BTreeMap<String, String> {
    let mut out = BTreeMap::new();
    for piece in info.split(';') {
        if piece.is_empty() {
            continue;
        }
        match piece.split_once('=') {
            Some((k, v)) => {
                out.insert(k.to_string(), v.to_string());
            }
            None => {
                out.insert(piece.to_string(), String::new());
            }
        }
    }
    out
}

fn java_target_record(path: &Path) -> (bool, BTreeMap<String, String>, String, usize) {
    let text = std::fs::read_to_string(path).expect("java.vcf");
    let mut header_declares_ic = false;
    let mut sample_n = 0usize;
    for line in text.lines() {
        if line.starts_with("##INFO=<ID=InbreedingCoeff") {
            header_declares_ic = true;
        }
        if line.starts_with("#CHROM") {
            sample_n = line.split('\t').count().saturating_sub(9);
            continue;
        }
        if line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 8 {
            continue;
        }
        if f[0] == "2" && f[1] == "92305634" && f[3] == MERGED_REF && f[4] == MERGED_ALT {
            return (
                header_declares_ic,
                parse_info_map(f[7]),
                line.to_string(),
                sample_n,
            );
        }
    }
    panic!("java.vcf missing 2:92305634 G/T");
}

/// Java `InbreedingCoeff.annotate` emit/suppress (SHA `2dbc0258`).
///
/// Gate 1: `genotypes.size() < MIN_SAMPLES` or `!vc.isVariant()` → emptyMap.
/// Gate 2: `calculateIC` sampleCount (called diploid with PL or GQ) `< MIN_SAMPLES` → emptyMap.
fn java_emits_inbreeding(
    n_vc_genotypes: usize,
    n_called_diploid_usable: usize,
    is_variant: bool,
) -> bool {
    if n_vc_genotypes < JAVA_MIN_SAMPLES || !is_variant {
        return false;
    }
    n_called_diploid_usable >= JAVA_MIN_SAMPLES
}

/// Production Rust formula in `annotate_hc_variant_site` (`genotype_counts_from_index`).
fn rust_ic_from_gt_index(best: usize) -> f64 {
    let (ref_n, het_n, hom_alt_n) = match best {
        0 => (1u32, 0, 0),
        1 => (0, 1, 0),
        _ => (0, 0, 1),
    };
    rust_ic_from_hard_counts(ref_n, het_n, hom_alt_n)
}

fn rust_ic_from_hard_counts(ref_n: u32, het_n: u32, hom_alt_n: u32) -> f64 {
    let n = ref_n + het_n + hom_alt_n;
    if n > 0 {
        1.0 - (het_n as f64) / n as f64
    } else {
        0.0
    }
}

/// 6R.178: `hc_info_values` inserts iff `n_genotypes >= MIN_SAMPLES`.
fn rust_emits_inbreeding(n_samples: usize) -> bool {
    n_samples >= JAVA_MIN_SAMPLES
}

/// GATK `MathUtils.normalizeFromLog10ToLinearSpace` on PL → [AA, AB, BB].
fn pl_to_normalized(pl: [i32; 3]) -> [f64; 3] {
    let log10 = [
        -(pl[0] as f64) / 10.0,
        -(pl[1] as f64) / 10.0,
        -(pl[2] as f64) / 10.0,
    ];
    let max = log10.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let lin = [
        10f64.powf(log10[0] - max),
        10f64.powf(log10[1] - max),
        10f64.powf(log10[2] - max),
    ];
    let s = lin[0] + lin[1] + lin[2];
    [lin[0] / s, lin[1] / s, lin[2] / s]
}

/// Java `InbreedingCoeff.calculateIC` with `ROUND_GENOTYPE_COUNTS = false`.
fn java_calculate_ic(normalized: &[[f64; 3]]) -> (usize, f64) {
    let mut ref_c = 0.0;
    let mut het_c = 0.0;
    let mut hom_c = 0.0;
    for n in normalized {
        ref_c += n[0];
        het_c += n[1];
        hom_c += n[2];
    }
    let sample_count = normalized.len();
    let p = (2.0 * ref_c + het_c) / (2.0 * (ref_c + het_c + hom_c));
    let q = 1.0 - p;
    let expected_hets = 2.0 * p * q * sample_count as f64;
    let f = 1.0 - (het_c / expected_hets);
    (sample_count, f)
}

#[test]
fn forensic_6r177_source_contract_no_production_change() {
    kv("java_pin", JAVA_PIN);
    kv(
        "java_path",
        "InbreedingCoeff.annotate → getFounderGenotypes → size < MIN_SAMPLES=10 → emptyMap",
    );
    kv(
        "rust_path",
        "annotate_hc_variant_site genotype_counts_from_index → 1-het/n; 6R.178 hc_info_values n>=10",
    );
    let emit = include_str!("../src/region_vcf_emit.rs");
    let ann = include_str!("../src/variant_site_hc_annotations.rs");
    let production_emit = emit.split("#[cfg(test)]").next().expect("production emit");
    assert!(
        production_emit.contains("\"InbreedingCoeff\".to_string()")
            && production_emit.contains("vec![ann.inbreeding_coeff]"),
        "production still has an InbreedingCoeff insert"
    );
    assert!(
        production_emit.contains("INBREEDING_COEFF_MIN_SAMPLES")
            && production_emit.contains("n_genotypes >= INBREEDING_COEFF_MIN_SAMPLES"),
        "6R.178 gates InbreedingCoeff on Java MIN_SAMPLES=10"
    );
    assert!(
        ann.contains("1.0 - (het_n as f64) / (ref_n + het_n + hom_alt_n) as f64"),
        "production formula is 1 - het/n, not Java HWE F"
    );
    assert!(
        !ann.contains("92305634") && !emit.contains("92305634"),
        "no locus-specific InbreedingCoeff patch"
    );
    assert_eq!(JAVA_MIN_SAMPLES, 10);
}

#[test]
fn forensic_6r177_same_input_emit_matrix() {
    kv("java_pin", JAVA_PIN);
    kv("min_samples", JAVA_MIN_SAMPLES.to_string());
    kv(
        "question",
        "n=1 omit is MIN_SAMPLES; n>=10 still uses 1-het/n vs Java HWE (stacked)",
    );

    // A: 1 sample, valid 1/1 (live site).
    let a_java = java_emits_inbreeding(1, 1, true);
    let a_rust_val = rust_ic_from_gt_index(2);
    let a_rust_emit = rust_emits_inbreeding(1);
    kv(
        "case_A",
        format!("1 sample 1/1 java_emit={a_java} rust_emit={a_rust_emit} rust_val={a_rust_val}"),
    );
    assert!(!a_java, "Java A: n=1 < 10 → emptyMap");
    assert!(!a_rust_emit, "6R.178: n=1 omits InbreedingCoeff");
    assert!(
        (a_rust_val - 1.0).abs() < 1e-12,
        "formula on 1/1 is still 1.0"
    );

    // A-het: 1 sample 0/1 still emits 0.0 in Rust — key exists because of insert, not F=1.
    let a_het = rust_ic_from_gt_index(1);
    kv(
        "case_A_het",
        format!(
            "1 sample 0/1 rust_val={a_het} rust_emit={}",
            rust_emits_inbreeding(1)
        ),
    );
    assert!((a_het - 0.0).abs() < 1e-12);
    assert!(
        !rust_emits_inbreeding(1),
        "6R.178 omits on n=1 even when the formula would be 0.0"
    );

    // B: 2 samples, valid genotypes.
    let b_java = java_emits_inbreeding(2, 2, true);
    kv(
        "case_B",
        format!(
            "2 samples java_emit={b_java} rust_emit={} rust_hard_2hom={:.4}",
            rust_emits_inbreeding(2),
            rust_ic_from_hard_counts(0, 0, 2)
        ),
    );
    assert!(!b_java);

    // C: 9 samples.
    let c_java = java_emits_inbreeding(9, 9, true);
    kv(
        "case_C",
        format!(
            "9 samples java_emit={c_java} rust_emit={}",
            rust_emits_inbreeding(9)
        ),
    );
    assert!(!c_java);

    // D: 10 samples, all called 1/1 with live PL — Java would emit HWE F, not omit.
    let d_java = java_emits_inbreeding(10, 10, true);
    let one = pl_to_normalized(LIVE_PL);
    let ten = vec![one; 10];
    let (d_n, d_f) = java_calculate_ic(&ten);
    kv(
        "case_D",
        format!(
            "10 samples 1/1 PL=90,6,0 java_emit={d_java} java_F={d_f:.4} sampleCount={d_n} rust_hard={:.4}",
            rust_ic_from_hard_counts(0, 0, 10)
        ),
    );
    assert!(d_java, "Java D: size=10 and called=10 → emit F");
    assert_eq!(d_n, 10);
    assert!(
        (d_f - 1.0).abs() > 0.05,
        "Java HWE F on 10× live PL is not Rust 1-het/n=1.0 (stacked later formula)"
    );

    // E: 10 VC genotypes, only 1 called diploid with PL; 9 no-call.
    let e_java = java_emits_inbreeding(10, 1, true);
    kv(
        "case_E",
        format!("10 genotypes / 1 called java_emit={e_java} (second MIN_SAMPLES gate)"),
    );
    assert!(!e_java, "Java E: size=10 but called=1 < 10 → emptyMap");

    assert!(
        !a_java && !a_rust_emit,
        "6R.178: n=1 emit predicate matches Java emptyMap; formula split remains at n>=10"
    );
}

#[test]
fn forensic_6r177_live_target_one_sample_hc() {
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    let java_vcf = root.join(JAVA_VCF_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing P12 ref/BAM");
        return;
    }

    kv("java_pin", JAVA_PIN);
    kv("target", "2:92305634 G/T");
    kv(
        "production_change",
        "6R.178: hc_info_values insert iff n_genotypes >= 10",
    );

    if java_vcf.is_file() {
        let (hdr_ic, jinfo, jline, jsamples) = java_target_record(&java_vcf);
        kv("java_record", jline);
        kv("java_sample_count", jsamples.to_string());
        kv("java_header_declares_ic", hdr_ic.to_string());
        kv(
            "java_ic_key",
            if jinfo.contains_key("InbreedingCoeff") {
                jinfo["InbreedingCoeff"].clone()
            } else {
                "OMITTED".to_string()
            },
        );
        kv("java_mq", jinfo.get("MQ").cloned().unwrap_or_default());
        assert!(
            hdr_ic,
            "Java default HC declares InbreedingCoeff in the header"
        );
        assert_eq!(jsamples, 1, "Java VCF is a one-sample HC call");
        assert!(
            !jinfo.contains_key("InbreedingCoeff"),
            "Java record omits InbreedingCoeff"
        );
        assert_eq!(jinfo.get("DP").map(String::as_str), Some("3"));
        assert_eq!(jinfo.get("SOR").map(String::as_str), Some("0.693"));
        assert_eq!(jinfo.get("MQ").map(String::as_str), Some("41.96"));
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
    let covering = regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= TARGET
                && r.end.get() >= TARGET
        })
        .expect("covering");

    let outcome = HaplotypeCallerEngine::call_region(
        covering,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("outcome");
    let call = outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(TARGET)
                && c.event.ref_allele == MERGED_REF
                && c.event.alt_allele == MERGED_ALT
        })
        .expect("genotyped G/T");
    let fmt = &call.genotype.format;
    assert_eq!(fmt.ad_as_i32(), vec![0, 2]);
    assert_eq!(fmt.pl_as_i32(), vec![90, 6, 0]);
    assert_eq!(fmt.dp.as_i32(), 2);
    assert_eq!(fmt.gq.as_i32(), 6);

    let emitted =
        try_emit_call_region_variants(covering, &outcome, "NA12878", DEFAULT_STAND_EMIT_CONFIDENCE)
            .expect("emit");
    let rec = emitted
        .iter()
        .find(|r| r.position == TARGET)
        .expect("emitted");
    assert_eq!(
        rec.samples.len(),
        1,
        "Rust emit is one-sample, not a grouping mismatch"
    );
    let sample = rec.samples.first().expect("sample");
    kv(
        "rust_sample",
        format!(
            "n={} name=NA12878 GT={:?} AD={:?} DP={:?} GQ={:?} PL={:?}",
            rec.samples.len(),
            sample.gt.as_ref().map(|g| g.to_string()),
            sample.ad,
            sample.dp,
            sample.gq,
            sample.pl
        ),
    );
    kv(
        "rust_ic_inputs",
        "best_idx=2 (1/1) → genotype_counts (ref=0, het=0, hom_alt=1) → 1-0/1=1.0",
    );
    assert_eq!(info_i32(&rec.info, "DP"), Some(3));
    assert_eq!(sample.dp.map(|d| d as i32), Some(2));
    let sor = info_f64(&rec.info, "SOR").unwrap_or(-1.0);
    assert!(
        (sor - 0.6931471805599453).abs() < 1e-6 || (sor - 0.693).abs() < 0.002,
        "6R.176 SOR stays 0.693, got {sor}"
    );
    assert!((rec.quality.unwrap_or(0.0) - 78.32).abs() < 0.02);
    let ic = info_f64(&rec.info, "InbreedingCoeff");
    kv("rust_ic", format!("{ic:?}"));
    assert!(
        !info_has(&rec.info, "InbreedingCoeff"),
        "6R.178: n=1 omits InbreedingCoeff, got {ic:?}"
    );
    let mq = info_f64(&rec.info, "MQ").unwrap_or(-1.0);
    kv(
        "mq_out_of_scope",
        format!("Java=41.96 Rust={mq} — not investigated"),
    );
    assert!((mq - 41.96).abs() < 0.005, "MQ recorded, not this arrow");
    assert!(!java_emits_inbreeding(1, 1, true));
    kv(
        "first_arrow",
        "6R.178 closed n<10 omit; remaining stacked formula is Java HWE vs 1-het/n",
    );
}
