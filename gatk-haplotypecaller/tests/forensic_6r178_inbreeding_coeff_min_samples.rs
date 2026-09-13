//! 6R.178: InbreedingCoeff `MIN_SAMPLES=10` emission gate at `2:92305634 G/T`.
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! PRODUCTION CHANGE: `hc_info_values` inserts InbreedingCoeff only when
//! `n_genotypes >= 10`. Formula in `annotate_hc_variant_site` is unchanged
//! (`1 - het/n`, not Java HWE `1 - het/(2pq n)`).
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r178_inbreeding_coeff_min_samples -- --nocapture --test-threads=1
//! HOLDOUT_6R178=1 cargo test -p gatk-haplotypecaller --test holdout_6r178_inbreeding -- --nocapture --test-threads=1
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
const JAVA_MIN_SAMPLES: usize = 10;
const RUST_FORMULA: &str = "1.0 - (het_n as f64) / (ref_n + het_n + hom_alt_n) as f64";

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R178\t{key}\t{}", value.as_ref());
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

fn rust_emits_inbreeding(n_genotypes: usize) -> bool {
    n_genotypes >= JAVA_MIN_SAMPLES
}

#[test]
fn forensic_6r178_source_contract_min_samples_gate_only() {
    kv("java_pin", JAVA_PIN);
    kv(
        "java_path",
        "InbreedingCoeff.annotate → genotypes.size() < MIN_SAMPLES=10 → emptyMap",
    );
    kv(
        "rust_path",
        "hc_info_values(ann, samples.len()) insert iff n_genotypes >= INBREEDING_COEFF_MIN_SAMPLES",
    );
    let emit = include_str!("../src/region_vcf_emit.rs");
    let ann = include_str!("../src/variant_site_hc_annotations.rs");
    let production_emit = emit.split("#[cfg(test)]").next().expect("production emit");
    assert!(
        production_emit.contains("const INBREEDING_COEFF_MIN_SAMPLES: usize = 10"),
        "production must name Java MIN_SAMPLES=10"
    );
    assert!(
        production_emit.contains("if n_genotypes >= INBREEDING_COEFF_MIN_SAMPLES"),
        "insert gate is genotype/sample count, not a value or GT special case"
    );
    assert!(
        production_emit.contains("hc_info_values(ann, samples.len())"),
        "cardinality is the emitted genotype collection (6R.177 population)"
    );
    assert!(
        ann.contains(RUST_FORMULA),
        "6R.178 must not change the 1-het/n formula"
    );
    assert!(
        !ann.contains("2.0 * p * q")
            && !ann.contains("expected_hets")
            && !ann.contains("INBREEDING_COEFF_MIN_SAMPLES"),
        "annotate_hc_variant_site must not grow Java HWE F or the emit gate"
    );
    assert!(
        !production_emit.contains("92305634")
            && !production_emit.contains("H06JUADXX")
            && !ann.contains("92305634"),
        "no locus/QNAME-specific InbreedingCoeff patch"
    );
    assert!(
        !production_emit.contains("inbreeding_coeff == 1")
            && !production_emit.contains("inbreeding_coeff == 0")
            && !production_emit.contains("ann.inbreeding_coeff == 1.0"),
        "must not suppress by value==1 or value==0"
    );
    assert!(
        production_emit.contains("\"InbreedingCoeff\"")
            && production_emit.contains(
                "Inbreeding coefficient as estimated from the genotype likelihoods per-sample"
            ),
        "header declaration of InbreedingCoeff must remain"
    );
    assert_eq!(JAVA_MIN_SAMPLES, 10);
    assert!(!rust_emits_inbreeding(1));
    assert!(!rust_emits_inbreeding(9));
    assert!(rust_emits_inbreeding(10));
}

#[test]
fn forensic_6r178_live_target_omits_inbreeding_coeff() {
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
        "hc_info_values insert InbreedingCoeff iff n_genotypes >= 10",
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
        assert!(hdr_ic, "Java header still declares InbreedingCoeff");
        assert_eq!(jsamples, 1);
        assert!(
            !jinfo.contains_key("InbreedingCoeff"),
            "Java record omits InbreedingCoeff"
        );
        assert_eq!(jinfo.get("DP").map(String::as_str), Some("3"));
        assert_eq!(jinfo.get("SOR").map(String::as_str), Some("0.693"));
        assert_eq!(jinfo.get("FS").map(String::as_str), Some("0.000"));
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
    assert_eq!(rec.samples.len(), 1, "same one-sample population as 6R.177");
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
    assert_eq!(
        sample.gt.as_ref().map(|g| g.to_string()).as_deref(),
        Some("1/1")
    );
    assert_eq!(sample.ad.as_deref(), Some(&[0u32, 2][..]));
    assert_eq!(sample.dp.map(|d| d as i32), Some(2));
    assert_eq!(sample.gq.map(|g| g as i32), Some(6));
    assert_eq!(sample.pl.as_deref(), Some(&[90u32, 6, 0][..]));
    assert_eq!(info_i32(&rec.info, "DP"), Some(3));
    let sor = info_f64(&rec.info, "SOR").unwrap_or(-1.0);
    let fs = info_f64(&rec.info, "FS").unwrap_or(-1.0);
    let mq = info_f64(&rec.info, "MQ").unwrap_or(-1.0);
    assert!(
        (sor - 0.6931471805599453).abs() < 1e-6 || (sor - 0.693).abs() < 0.002,
        "6R.176 SOR stays 0.693, got {sor}"
    );
    assert!(fs < 0.02, "FS stays 0, got {fs}");
    assert!((mq - 41.96).abs() < 0.005, "MQ stays 41.96, got {mq}");
    assert!((rec.quality.unwrap_or(0.0) - 78.32).abs() < 0.02);
    assert!(
        !info_has(&rec.info, "InbreedingCoeff"),
        "n=1 < 10 must omit InbreedingCoeff entirely, got {:?}",
        info_f64(&rec.info, "InbreedingCoeff")
    );
    kv("rust_ic", "OMITTED");
    kv(
        "formula",
        "UNCHANGED — annotate_hc_variant_site still 1-het/n; Java HWE is later",
    );
    kv(
        "next",
        "Java calculateIC F=1-het/(2pq n) vs Rust 1-het/n (only when n>=10)",
    );
}
