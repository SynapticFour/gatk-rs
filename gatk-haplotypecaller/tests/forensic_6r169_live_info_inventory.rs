//! 6R.169: live INFO inventory at `2:92316347 G/A` (proof-only).
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`, GKL 0.8.8.
//! Production change: NONE. FORMAT/QUAL/FS/SOR/MQ stay closed through 6R.168.
//!
//! Java default HC groups: `StandardAnnotation` + `StandardHCAnnotation`
//! (`HaplotypeCaller.getDefaultVariantAnnotationGroups`).
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r169_live_info_inventory -- --test-threads=1 --nocapture
//! HOLDOUT_6R169=1 cargo test -p gatk-haplotypecaller --test holdout_6r169_live_info -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::fragment_overlap::read_base_at_ref_coord_1based;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

const INTERVAL: &str = "2:92316200-92316580";
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const JAVA_VCF_REL: &str = "parity/reports/6r43/p12_mid_a/java.vcf";
const TARGET: u64 = 92_316_347;
const MERGED_REF: &str = "G";
const MERGED_ALT: &str = "A";
/// Java `InbreedingCoeff.MIN_SAMPLES`.
const JAVA_INBREEDING_MIN_SAMPLES: usize = 10;
/// Java `RankSumTest` + `MannWhitneyU.test`: empty REF or ALT → NaN → omit key.

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R169\t{key}\t{}", value.as_ref());
}

fn info_key(v: &InfoValue) -> &str {
    match v {
        InfoValue::Flag(k)
        | InfoValue::Integer(k, _)
        | InfoValue::Float(k, _)
        | InfoValue::String(k, _)
        | InfoValue::Character(k, _) => k.as_str(),
    }
}

fn info_display(v: &InfoValue) -> String {
    match v {
        InfoValue::Flag(k) => k.clone(),
        InfoValue::Integer(k, xs) => format!(
            "{k}={}",
            xs.iter()
                .map(|x| x.to_string())
                .collect::<Vec<_>>()
                .join(",")
        ),
        InfoValue::Float(k, xs) => format!(
            "{k}={}",
            xs.iter()
                .map(|x| format!("{x}"))
                .collect::<Vec<_>>()
                .join(",")
        ),
        InfoValue::String(k, xs) => format!("{k}={}", xs.join(",")),
        InfoValue::Character(k, xs) => format!(
            "{k}={}",
            xs.iter()
                .map(|c| c.to_string())
                .collect::<Vec<_>>()
                .join(",")
        ),
    }
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

fn java_target_record(path: &Path) -> (BTreeSet<String>, BTreeMap<String, String>, String) {
    let text = std::fs::read_to_string(path).expect("java.vcf");
    let mut header_keys = BTreeSet::new();
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("##INFO=<ID=") {
            let id = rest.split(',').next().unwrap_or("");
            header_keys.insert(id.to_string());
        }
        if line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 8 {
            continue;
        }
        if f[0] == "2" && f[1] == "92316347" && f[3] == MERGED_REF && f[4] == MERGED_ALT {
            return (header_keys, parse_info_map(f[7]), line.to_string());
        }
    }
    panic!("java.vcf missing 2:92316347 G/A");
}

fn rust_info_map(info: &[InfoValue]) -> BTreeMap<String, String> {
    let mut out = BTreeMap::new();
    for v in info {
        let k = info_key(v).to_string();
        let rendered = info_display(v);
        let val = rendered
            .split_once('=')
            .map(|(_, v)| v.to_string())
            .unwrap_or_default();
        out.insert(k, val);
    }
    out
}

fn numeric_close(java: &str, rust: &str) -> bool {
    match (java.parse::<f64>(), rust.parse::<f64>()) {
        (Ok(j), Ok(r)) => {
            if j == r {
                return true;
            }
            // Java annotators format FS/SOR %.3f, QD/MQ %.2f, ExcessHet %.4f, AF %.2f, MLEAF %.3f.
            let places = java
                .split('.')
                .nth(1)
                .map(|frac| frac.len())
                .unwrap_or(0)
                .min(6);
            let scale = 10f64.powi(places as i32);
            if places == 0 {
                (j - r).abs() < 1e-12
            } else {
                (j * scale).round() == (r * scale).round() || (j - r).abs() < 0.5 / scale
            }
        }
        _ => java == rust,
    }
}

#[test]
fn forensic_6r169_production_dispatch_unchanged() {
    let emit_src = include_str!("../src/region_vcf_emit.rs");
    assert!(
        emit_src.contains("if let Some(z) = ann.read_pos_rank_sum")
            && emit_src.contains("n_genotypes >= INBREEDING_COEFF_MIN_SAMPLES"),
        "6R.172: ReadPosRankSum is optional Some(z); 6R.178 gates InbreedingCoeff on n>=10"
    );
    let ann_src = include_str!("../src/variant_site_hc_annotations.rs");
    assert!(
        ann_src.contains("read_pos_rank_sum::read_pos_rank_sum")
            && ann_src.contains("inbreeding_coeff"),
        "production annotators must still run"
    );
    assert!(
        !ann_src.contains("92316347") && !emit_src.contains("92316347"),
        "no locus-specific INFO patch"
    );
    assert_eq!(JAVA_INBREEDING_MIN_SAMPLES, 10);
}

#[test]
fn forensic_6r169_live_info_inventory() {
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    let java_vcf = root.join(JAVA_VCF_REL);
    if !ref_fasta.is_file() || !bam.is_file() || !java_vcf.is_file() {
        eprintln!("skip: missing P12 ref/BAM/java.vcf");
        return;
    }

    kv(
        "java_pin",
        "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77",
    );
    kv("variant", format!("2:{TARGET} G/A"));
    kv("production_change", "NONE");
    kv(
        "java_default_groups",
        "StandardAnnotation + StandardHCAnnotation",
    );

    let (java_header, java_info, java_line) = java_target_record(&java_vcf);
    kv("java_header_info", format!("{java_header:?}"));
    kv("java_record", java_line);
    let java_keys: BTreeSet<String> = java_info.keys().cloned().collect();
    kv("java_keys", format!("{java_keys:?}"));

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
        .expect("ActiveFull covering target");

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
        .expect("genotyped G/A");
    let fmt = &call.genotype.format;
    assert_eq!(fmt.ad_as_i32(), vec![0, 3]);
    assert_eq!(fmt.pl_as_i32(), vec![135, 9, 0]);
    assert_eq!(fmt.dp.as_i32(), 3);
    assert_eq!(fmt.gq.as_i32(), 9);

    let emitted =
        try_emit_call_region_variants(covering, &outcome, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .expect("emit");
    let rec = emitted
        .iter()
        .find(|r| {
            r.position == TARGET
                && r.reference == MERGED_REF
                && r.alternate.first().map(String::as_str) == Some(MERGED_ALT)
        })
        .expect("emitted G/A");
    assert!((rec.quality.unwrap_or(0.0) - 121.84).abs() < 0.02);

    let rust_info = rust_info_map(&rec.info);
    let rust_keys: BTreeSet<String> = rust_info.keys().cloned().collect();
    kv("rust_keys", format!("{rust_keys:?}"));
    for v in &rec.info {
        kv("rust_info", info_display(v));
    }

    let common: BTreeSet<_> = java_keys.intersection(&rust_keys).cloned().collect();
    let java_only: BTreeSet<_> = java_keys.difference(&rust_keys).cloned().collect();
    let rust_only: BTreeSet<_> = rust_keys.difference(&java_keys).cloned().collect();
    kv("common", format!("{common:?}"));
    kv("java_only", format!("{java_only:?}"));
    kv("rust_only", format!("{rust_only:?}"));

    let mut common_different = Vec::new();
    let mut common_match = Vec::new();
    for k in &common {
        let j = java_info.get(k).map(String::as_str).unwrap_or("");
        let r = rust_info.get(k).map(String::as_str).unwrap_or("");
        if numeric_close(j, r) {
            common_match.push(k.clone());
            kv("common_match", format!("{k}\tjava={j}\trust={r}"));
        } else {
            common_different.push(k.clone());
            kv("common_different", format!("{k}\tjava={j}\trust={r}"));
        }
    }

    assert!(
        java_only.is_empty(),
        "unexpected Java-only INFO keys: {java_only:?}"
    );
    assert_eq!(
        rust_only,
        ["ReadPosRankSum"]
            .into_iter()
            .map(str::to_string)
            .collect::<BTreeSet<_>>(),
        "live Rust-only INFO must be exactly ReadPosRankSum, got {rust_only:?}"
    );
    assert!(
        common_different.is_empty(),
        "unexpected common-key value splits beyond Java print precision: {common_different:?}"
    );

    // Closed 6R.164–6R.168 fields.
    assert_eq!(java_info.get("FS").map(String::as_str), Some("0.000"));
    assert_eq!(java_info.get("MQ").map(String::as_str), Some("40.25"));
    assert_eq!(java_info.get("SOR").map(String::as_str), Some("1.179"));
    let rust_fs: f64 = rust_info
        .get("FS")
        .and_then(|s| s.parse().ok())
        .unwrap_or(-1.0);
    let rust_mq: f64 = rust_info
        .get("MQ")
        .and_then(|s| s.parse().ok())
        .unwrap_or(-1.0);
    let rust_sor: f64 = rust_info
        .get("SOR")
        .and_then(|s| s.parse().ok())
        .unwrap_or(-1.0);
    assert!(rust_fs < 0.02, "FS must stay Java 0, got {rust_fs}");
    assert!(
        (rust_mq - 40.25).abs() < 1e-12,
        "MQ must stay 40.25, got {rust_mq}"
    );
    assert!(
        (rust_sor - 1.1786549963416462).abs() < 1e-9 || (rust_sor - 1.179).abs() < 0.001,
        "SOR must stay Java 1.179, got {rust_sor}"
    );

    // Java header declares rank-sum / IC keys; this record omits them.
    for k in [
        "ReadPosRankSum",
        "InbreedingCoeff",
        "BaseQRankSum",
        "MQRankSum",
    ] {
        assert!(
            java_header.contains(k),
            "Java default HC header must declare {k}"
        );
        assert!(
            !java_keys.contains(k),
            "Java record at this site must omit {k}"
        );
    }

    let mut pileup_ref = 0usize;
    let mut pileup_alt = 0usize;
    let ref_b = MERGED_REF.as_bytes()[0];
    let alt_b = MERGED_ALT.as_bytes()[0];
    for rec in &covering.reads {
        let Some(base) = read_base_at_ref_coord_1based(rec, TARGET as i32) else {
            continue;
        };
        if base.eq_ignore_ascii_case(&alt_b) {
            pileup_alt += 1;
        } else if base.eq_ignore_ascii_case(&ref_b) {
            pileup_ref += 1;
        }
    }
    kv(
        "rust_pileup_ranksum_n",
        format!("ref={pileup_ref} alt={pileup_alt}"),
    );
    kv(
        "java_ranksum_evidence",
        "informative likelihoods REF=0 ALT=3 (6R.166 table [0,0;2,1]) → MannWhitneyU NaN → omit",
    );
    kv(
        "java_inbreeding",
        format!("enabled StandardAnnotation; MIN_SAMPLES={JAVA_INBREEDING_MIN_SAMPLES}; n_samples=1 → emptyMap"),
    );
    kv(
        "rust_read_pos_rank_sum",
        rust_info
            .get("ReadPosRankSum")
            .cloned()
            .unwrap_or_else(|| "ABSENT".into()),
    );
    kv(
        "rust_inbreeding_coeff",
        rust_info
            .get("InbreedingCoeff")
            .cloned()
            .unwrap_or_else(|| "ABSENT".into()),
    );

    assert!(
        pileup_ref > 0 && pileup_alt > 0,
        "Rust pileup at this site still has both REF and ALT bases (6R.165 [0,2;3,2] class)"
    );
    assert!(
        rust_info.get("InbreedingCoeff").is_none(),
        "6R.178: n=1 omits InbreedingCoeff"
    );
}
