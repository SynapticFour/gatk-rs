//! 6R.197: proof-only fresh VCF reconnaissance after 6R.196.
//! After 6R.239, `20:29455379 G/A` FORMAT matches Java in the covering dump.
//! First remaining genuine common-site split in that dump is
//! `20:29455649 T/TGTTTG` (PL/INFO). Rust-only `20:29455314 G>C` is 6R.240.
//! PRODUCTION CHANGE: NONE.
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! PRODUCTION CHANGE: NONE.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r197_fresh_info_reconnaissance -- --nocapture --test-threads=1
//! HOLDOUT_6R197=1 cargo test -p gatk-haplotypecaller --test holdout_6r197_fresh_info -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};

const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
const CURRENT_REGIONS: &[&str] = &[
    "ctrl_mid_b",
    "p12_snp_cluster",
    "p12_indel_mix",
    "p12_mid_a",
    "p12_post",
    "p12_het_tail",
    "p12_desert",
    "chr20_tiny",
];
const STALE_MISSING_BAM: &[&str] = &["chr20_w47", "chr21_w10"];
const QUAL_TOL: f64 = 0.05;
const FIRST_REMAINING: (&str, u64, &str, &str) = ("20", 29_455_649, "T", "TGTTTG");
const CLOSED_GA: (&str, u64, &str, &str) = ("20", 29_455_379, "G", "A");
const CLOSED_CHR20_PL: (&str, u64, &str, &str) = ("20", 29_455_015, "G", "T");
const CLOSED_WEAK: (&str, u64, &str, &str) = ("2", 92_325_268, "C", "T");
const CLOSED_HET_TAIL: (&str, u64, &str, &str) = ("2", 92_325_193, "C", "T");
const CLOSED_HET_SIB: (&str, u64, &str, &str) = ("2", 92_325_205, "G", "A");
const CLOSED_POST: (&str, u64, &str, &str) = ("2", 92_318_199, "C", "T");
const CLOSED_MID_B: (&str, u64, &str, &str) = ("2", 92_317_399, "C", "A");
const CLOSED_HOM_ALT: (&str, u64, &str, &str) = ("2", 92_316_296, "A", "T");
const CLOSED_ONE_READ: (&str, u64, &str, &str) = ("2", 92_316_416, "C", "A");
const CLOSED_CA: (&str, u64, &str, &str) = ("2", 92_307_403, "C", "A");
const CLOSED_GT: (&str, u64, &str, &str) = ("2", 92_305_634, "G", "T");
const CLOSED_AG: (&str, u64, &str, &str) = ("2", 92_305_635, "A", "G");
const CLOSED_AC: (&str, u64, &str, &str) = ("2", 92_305_716, "A", "C");
const CLOSED_INDEL: (&str, u64, &str, &str) = ("2", 92_307_324, "TTC", "T");
const CLOSED_TG: (&str, u64, &str, &str) = ("2", 92_307_333, "T", "G");
const CLOSED_QUAL: (&str, u64, &str, &str) = ("2", 92_307_359, "CT", "C");
const CLOSED_SNP_PRIOR: (&str, u64, &str, &str) = ("2", 92_307_364, "T", "C");
const CLOSED_MID: (&str, u64, &str, &str) = ("2", 92_316_347, "G", "A");
const INFO_KEYS: &[&str] = &[
    "AC",
    "AF",
    "AN",
    "DP",
    "ExcessHet",
    "FS",
    "InbreedingCoeff",
    "MLEAC",
    "MLEAF",
    "MQ",
    "QD",
    "SOR",
    "ReadPosRankSum",
    "BaseQRankSum",
    "MQRankSum",
];
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const GAP_INTERVAL: &str = "2:92305500-92305850";
const TG_INTERVAL: &str = "2:92307200-92307550";

#[derive(Clone, Debug)]
struct Rec {
    chrom: String,
    pos: u64,
    ref_a: String,
    alt: String,
    qual: String,
    filter: String,
    info: BTreeMap<String, String>,
    gt: String,
    ad: String,
    dp: String,
    gq: String,
    pl: String,
}

impl Rec {
    fn key(&self) -> (String, u64, String, String) {
        (
            self.chrom.clone(),
            self.pos,
            self.ref_a.clone(),
            self.alt.clone(),
        )
    }

    fn site(&self) -> String {
        format!("{}:{} {}/{}", self.chrom, self.pos, self.ref_a, self.alt)
    }
}

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R197\t{key}\t{}", value.as_ref());
}

fn parse_info(s: &str) -> BTreeMap<String, String> {
    let mut out = BTreeMap::new();
    if s.is_empty() || s == "." {
        return out;
    }
    for part in s.split(';') {
        if part.is_empty() {
            continue;
        }
        if let Some((k, v)) = part.split_once('=') {
            out.insert(k.to_string(), v.to_string());
        } else {
            out.insert(part.to_string(), "true".to_string());
        }
    }
    out
}

fn parse_vcf(path: &Path) -> Vec<Rec> {
    let text = fs::read_to_string(path).unwrap_or_default();
    let mut out = Vec::new();
    for line in text.lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 10 {
            continue;
        }
        let mut sample = BTreeMap::new();
        for (k, v) in f[8].split(':').zip(f[9].split(':')) {
            sample.insert(k.to_string(), v.to_string());
        }
        let alt0 = f[4].split(',').next().unwrap_or(f[4]);
        out.push(Rec {
            chrom: f[0].to_string(),
            pos: f[1].parse().expect("pos"),
            ref_a: f[3].to_string(),
            alt: alt0.to_string(),
            qual: f[5].to_string(),
            filter: f[6].to_string(),
            info: parse_info(f[7]),
            gt: sample.get("GT").cloned().unwrap_or_default(),
            ad: sample.get("AD").cloned().unwrap_or_default(),
            dp: sample.get("DP").cloned().unwrap_or_default(),
            gq: sample.get("GQ").cloned().unwrap_or_default(),
            pl: sample.get("PL").cloned().unwrap_or_default(),
        });
    }
    out
}

fn format_equal(j: &Rec, r: &Rec) -> bool {
    j.gt == r.gt
        && j.ad == r.ad
        && j.dp == r.dp
        && j.gq == r.gq
        && j.pl == r.pl
        && j.filter == r.filter
}

fn first_format_field(j: &Rec, r: &Rec) -> Option<&'static str> {
    if j.gt != r.gt {
        Some("GT")
    } else if j.ad != r.ad {
        Some("AD")
    } else if j.dp != r.dp {
        Some("DP")
    } else if j.gq != r.gq {
        Some("GQ")
    } else if j.pl != r.pl {
        Some("PL")
    } else if j.filter != r.filter {
        Some("FILTER")
    } else {
        None
    }
}

fn qual_close(a: &str, b: &str) -> bool {
    match (a.parse::<f64>(), b.parse::<f64>()) {
        (Ok(x), Ok(y)) => (x - y).abs() <= QUAL_TOL,
        _ => a == b,
    }
}

fn info_print_close(key: &str, jv: &str, rv: &str) -> bool {
    if jv == rv {
        return true;
    }
    let (Ok(j), Ok(r)) = (jv.parse::<f64>(), rv.parse::<f64>()) else {
        return false;
    };
    match key {
        "FS" | "SOR" | "ExcessHet" | "ReadPosRankSum" | "BaseQRankSum" | "MQRankSum" => {
            (j * 1000.0).round() == (r * 1000.0).round()
        }
        "MQ" | "QD" | "AF" | "MLEAF" => {
            (j * 100.0).round() == (r * 100.0).round() || (j - r).abs() < 0.005
        }
        _ => false,
    }
}

fn load_region(root: &Path, id: &str) -> (Vec<Rec>, Vec<Rec>) {
    let dir = root.join("parity/reports/6r43").join(id);
    (
        parse_vcf(&dir.join("java.vcf")),
        parse_vcf(&dir.join("rust.vcf")),
    )
}

fn rec_key(spec: (&str, u64, &str, &str)) -> (String, u64, String, String) {
    (
        spec.0.to_string(),
        spec.1,
        spec.2.to_string(),
        spec.3.to_string(),
    )
}

fn genuine_info_value_diffs(j: &Rec, r: &Rec) -> Vec<(String, String, String)> {
    let mut out = Vec::new();
    for (k, jv) in &j.info {
        let Some(rv) = r.info.get(k) else {
            continue;
        };
        if jv == rv || info_print_close(k, jv, rv) {
            continue;
        }
        out.push((k.clone(), jv.clone(), rv.clone()));
    }
    out
}

fn rust_only_info_keys(j: &Rec, r: &Rec) -> Vec<String> {
    r.info
        .keys()
        .filter(|k| !j.info.contains_key(*k))
        .cloned()
        .collect()
}

fn java_only_info_keys(j: &Rec, r: &Rec) -> Vec<String> {
    j.info
        .keys()
        .filter(|k| !r.info.contains_key(*k))
        .cloned()
        .collect()
}

fn unique_likelihood_indices(
    likelihoods: &[gatk_haplotypecaller::RegionReadLikelihood],
) -> BTreeSet<usize> {
    likelihoods.iter().map(|c| c.read_index.get()).collect()
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

#[test]
fn forensic_6r197_no_production_change() {
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");
    kv(
        "fresh_pass",
        "post-6R.196 production rust.vcf vs frozen Java 6R.43/6R.160 corpus",
    );
    kv(
        "classification",
        "G — UPSTREAM OBJECT / LIFECYCLE (provisional; causal arrow not yet proven)",
    );

    let af = include_str!("../src/af_calc.rs");
    let emit = include_str!("../src/region_vcf_emit.rs");
    let ann = include_str!("../src/variant_site_hc_annotations.rs");
    assert!(
        af.contains("biallelic_alt_pseudocount")
            && af.contains("alt_allele.len() == ref_allele.len()")
            && ann.contains("calculate_biallelic_af_em_with_alt_pseudocount"),
        "6R.196 length-based AF alt prior stays"
    );
    assert!(
        emit.contains("\"BaseQRankSum\"")
            && emit.contains("\"MQRankSum\"")
            && emit.contains("\"ReadPosRankSum\""),
        "RankSum headers/emit from 6R.191–6R.194 stay; reconnaissance does not unwire them"
    );
    assert!(
        !af.contains("6R.197") && !emit.contains("6R.197") && !ann.contains("6R.197"),
        "no 6R.197 production patch"
    );
}

#[test]
fn forensic_6r197_live_vcf_inventory() {
    let root = repo_root();
    let reports = root.join("parity/reports/6r43");
    let has_java = CURRENT_REGIONS
        .iter()
        .any(|id| reports.join(id).join("java.vcf").is_file());
    if !has_java {
        eprintln!("skip: missing parity/reports/6r43 Java VCFs");
        return;
    }
    let mut java_map = BTreeMap::new();
    let mut rust_map = BTreeMap::new();
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");
    kv(
        "fresh_pass",
        "HOLDOUT_6R43 rust.vcf after 6R.212 vs frozen Java 4.4.0.0",
    );
    kv("excluded", STALE_MISSING_BAM.join(","));
    for id in CURRENT_REGIONS {
        let (jv, rv) = load_region(&root, id);
        kv(
            "region",
            format!("{id}\tjava={}\trust={}", jv.len(), rv.len()),
        );
        for rec in jv {
            java_map.insert(rec.key(), rec);
        }
        for rec in rv {
            rust_map.insert(rec.key(), rec);
        }
    }
    for id in STALE_MISSING_BAM {
        let (jv, rv) = load_region(&root, id);
        kv(
            "excluded_stale",
            format!("{id}\tjava={}\trust={}\tBAM_missing", jv.len(), rv.len()),
        );
    }

    let java_keys: Vec<_> = java_map.keys().cloned().collect();
    let rust_keys: Vec<_> = rust_map.keys().cloned().collect();
    let java_only: Vec<_> = java_keys
        .iter()
        .filter(|k| !rust_map.contains_key(*k))
        .cloned()
        .collect();
    let rust_only: Vec<_> = rust_keys
        .iter()
        .filter(|k| !java_map.contains_key(*k))
        .cloned()
        .collect();
    let common: Vec<_> = java_keys
        .iter()
        .filter(|k| rust_map.contains_key(*k))
        .cloned()
        .collect();

    let mut fmt_n = 0usize;
    let mut qual_n = 0usize;
    let mut info_n = 0usize;
    let mut chr2_fmt = 0usize;
    let mut chr2_qual = 0usize;
    let mut extra_info_keys: BTreeSet<String> = BTreeSet::new();
    let mut first_fmt: Option<(String, &'static str)> = None;
    let mut first_qual: Option<String> = None;
    let mut first_info: Option<(Rec, Rec, Vec<(String, String, String)>)> = None;
    for k in &common {
        let j = &java_map[k];
        let r = &rust_map[k];
        for key in j.info.keys().chain(r.info.keys()) {
            if !INFO_KEYS.contains(&key.as_str()) {
                extra_info_keys.insert(key.clone());
            }
        }
        if let Some(fd) = first_format_field(j, r) {
            fmt_n += 1;
            if j.chrom == "2" {
                chr2_fmt += 1;
            }
            if first_fmt.is_none() {
                first_fmt = Some((j.site(), fd));
            }
        }
        if !qual_close(&j.qual, &r.qual) {
            qual_n += 1;
            if j.chrom == "2" {
                chr2_qual += 1;
            }
            if first_qual.is_none() {
                first_qual = Some(j.site());
            }
        }
        let val = genuine_info_value_diffs(j, r);
        let rust_only_keys = rust_only_info_keys(j, r);
        let java_only_keys = java_only_info_keys(j, r);
        if !val.is_empty() || !rust_only_keys.is_empty() || !java_only_keys.is_empty() {
            info_n += 1;
            if first_info.is_none() {
                first_info = Some((j.clone(), r.clone(), val));
            }
        }
    }
    kv("java_records", java_map.len().to_string());
    kv("rust_records", rust_map.len().to_string());
    kv("java_only", java_only.len().to_string());
    kv("rust_only", rust_only.len().to_string());
    kv("common", common.len().to_string());
    kv("common_format_diff", fmt_n.to_string());
    kv("common_qual_diff", qual_n.to_string());
    kv("common_info_diff_genuine", info_n.to_string());
    kv("chr2_format_diff", chr2_fmt.to_string());
    kv("chr2_qual_diff", chr2_qual.to_string());
    kv("extra_info_keys", format!("{extra_info_keys:?}"));
    if let Some((site, fd)) = &first_fmt {
        kv("first_format_diff", format!("{site} field={fd}"));
    }
    if let Some(site) = &first_qual {
        kv("first_qual_diff", site);
    }
    if let Some(k) = rust_only.first() {
        kv("first_rust_only", rust_map[k].site());
    }
    assert!(java_only.is_empty());
    assert_eq!(rust_only.len(), 4);
    assert_eq!(
        rust_map[rust_only.first().unwrap()].site(),
        "20:29455314 G/C"
    );
    assert!(
        rust_only.iter().all(|k| k.0 == "20"),
        "rust-only records stay on chr20_tiny, not chr2"
    );
    assert_eq!(chr2_fmt, 0, "chr2 FORMAT must stay closed after 6R.196");
    assert_eq!(chr2_qual, 0, "chr2 QUAL must stay closed within 0.05");
    assert_eq!(java_map.len(), 126);
    assert_eq!(rust_map.len(), 130);
    assert_eq!(common.len(), 126);
    assert!(
        first_fmt
            .as_ref()
            .is_some_and(|(site, fd)| site.starts_with("20:29455649") && *fd == "PL"),
        "6R.239 closed AD at 20:29455379; first remaining FORMAT split is PL at 20:29455649: {first_fmt:?}"
    );
    assert!(
        extra_info_keys.is_empty(),
        "no unexpected INFO keys beyond the compared set: {extra_info_keys:?}"
    );

    for (label, spec, dp, mq, sor) in [
        ("closed_92305634", CLOSED_GT, "3", 41.96, 0.693),
        ("closed_92305635", CLOSED_AG, "3", 41.96, 0.693),
        ("closed_92305716", CLOSED_AC, "4", 35.15, 1.179),
        ("closed_92307324", CLOSED_INDEL, "1", 44.0, 1.609),
        ("closed_92307333", CLOSED_TG, "1", 44.0, 1.609),
        ("closed_92316347", CLOSED_MID, "3", 40.25, 1.179),
        ("closed_92316296", CLOSED_HOM_ALT, "2", 47.00, 2.303),
        ("closed_92316416", CLOSED_ONE_READ, "1", 21.00, 1.609),
        ("closed_92317399", CLOSED_MID_B, "2", 27.00, 0.693),
        ("closed_92318199", CLOSED_POST, "1", 24.00, 1.609),
        ("closed_92325193", CLOSED_HET_TAIL, "3", 28.03, 0.223),
        ("closed_92325205", CLOSED_HET_SIB, "3", 28.03, 0.223),
        ("closed_92325268", CLOSED_WEAK, "3", 28.03, 1.179),
    ] {
        let k = rec_key(spec);
        let j = &java_map[&k];
        let r = &rust_map[&k];
        assert!(format_equal(j, r), "{label} FORMAT");
        assert!(qual_close(&j.qual, &r.qual), "{label} QUAL");
        assert_eq!(
            j.info.get("DP").map(String::as_str),
            Some(dp),
            "{label} java DP"
        );
        assert_eq!(
            r.info.get("DP").map(String::as_str),
            Some(dp),
            "{label} rust DP"
        );
        let jmq = j
            .info
            .get("MQ")
            .and_then(|s| s.parse::<f64>().ok())
            .unwrap_or(-1.0);
        let rmq = r
            .info
            .get("MQ")
            .and_then(|s| s.parse::<f64>().ok())
            .unwrap_or(-1.0);
        let jsor = j
            .info
            .get("SOR")
            .and_then(|s| s.parse::<f64>().ok())
            .unwrap_or(-1.0);
        let rsor = r
            .info
            .get("SOR")
            .and_then(|s| s.parse::<f64>().ok())
            .unwrap_or(-1.0);
        assert!((jmq - mq).abs() < 0.005, "{label} java MQ={jmq}");
        assert!((rmq - mq).abs() < 0.005, "{label} rust MQ={rmq}");
        assert!((jsor - sor).abs() < 0.002, "{label} java SOR={jsor}");
        assert!((rsor - sor).abs() < 0.002, "{label} rust SOR={rsor}");
        kv(label, format!("{} FORMAT/QUAL/DP/MQ/SOR closed", j.site()));
    }

    let ag = &rust_map[&rec_key(CLOSED_AG)];
    assert_eq!(ag.qual, "78.32");
    let indel = &rust_map[&rec_key(CLOSED_INDEL)];
    let jindel = &java_map[&rec_key(CLOSED_INDEL)];
    assert_eq!(indel.qual, "35.44");
    assert_eq!(jindel.qual, "35.44");
    let pin = rec_key(CLOSED_QUAL);
    let jp = &java_map[&pin];
    let rp = &rust_map[&pin];
    assert!(format_equal(jp, rp), "6R.196 pin FORMAT stays closed");
    assert_eq!(jp.qual, "31.60");
    assert_eq!(rp.qual, "31.60");
    assert_eq!(jp.info.get("QD").map(String::as_str), Some("15.80"));
    assert!(info_print_close(
        "QD",
        jp.info.get("QD").expect("java qd"),
        rp.info.get("QD").expect("rust qd"),
    ));
    assert!(info_print_close(
        "BaseQRankSum",
        jp.info.get("BaseQRankSum").expect("java bq"),
        rp.info.get("BaseQRankSum").expect("rust bq"),
    ));
    assert!(info_print_close(
        "MQRankSum",
        jp.info.get("MQRankSum").expect("java mqrs"),
        rp.info.get("MQRankSum").expect("rust mqrs"),
    ));
    assert!(info_print_close(
        "ReadPosRankSum",
        jp.info.get("ReadPosRankSum").expect("java rprs"),
        rp.info.get("ReadPosRankSum").expect("rust rprs"),
    ));
    assert!(
        genuine_info_value_diffs(jp, rp).is_empty()
            && rust_only_info_keys(jp, rp).is_empty()
            && java_only_info_keys(jp, rp).is_empty(),
        "6R.196 pin 2:92307359 stays closed; leftover={:?} rust_only={:?} java_only={:?}",
        genuine_info_value_diffs(jp, rp),
        rust_only_info_keys(jp, rp),
        java_only_info_keys(jp, rp)
    );
    kv(
        "closed_92307359",
        "2:92307359 CT/C QUAL=31.60 QD=15.80 RankSums match",
    );
    let snp = rec_key(CLOSED_SNP_PRIOR);
    let js = &java_map[&snp];
    let rs = &rust_map[&snp];
    assert!(format_equal(js, rs));
    assert_eq!(js.qual, "31.64");
    assert_eq!(rs.qual, "31.64");
    kv(
        "closed_92307364",
        "SNP QUAL 31.64 stays (same PL, SNP prior)",
    );

    let ca = rec_key(CLOSED_CA);
    let jca = &java_map[&ca];
    let rca = &rust_map[&ca];
    assert!(format_equal(jca, rca), "6R.203 pin FORMAT stays closed");
    assert_eq!(jca.qual, "154.64");
    assert_eq!(rca.qual, "154.64");
    assert_eq!(jca.gt, "0/1");
    assert_eq!(jca.ad, "2,4");
    assert_eq!(jca.dp, "6");
    assert_eq!(jca.gq, "72");
    assert_eq!(jca.pl, "162,0,72");
    assert_eq!(jca.info.get("DP").map(String::as_str), Some("6"));
    assert_eq!(rca.info.get("DP").map(String::as_str), Some("6"));
    assert!(info_print_close(
        "BaseQRankSum",
        jca.info.get("BaseQRankSum").expect("java bq"),
        rca.info.get("BaseQRankSum").expect("rust bq"),
    ));
    assert!(info_print_close(
        "ReadPosRankSum",
        jca.info.get("ReadPosRankSum").expect("java rprs"),
        rca.info.get("ReadPosRankSum").expect("rust rprs"),
    ));
    assert!(info_print_close(
        "MQRankSum",
        jca.info.get("MQRankSum").expect("java mqrs"),
        rca.info.get("MQRankSum").expect("rust mqrs"),
    ));
    assert!(
        genuine_info_value_diffs(jca, rca).is_empty()
            && rust_only_info_keys(jca, rca).is_empty()
            && java_only_info_keys(jca, rca).is_empty(),
        "6R.203 pin 2:92307403 stays closed; leftover={:?} rust_only={:?} java_only={:?}",
        genuine_info_value_diffs(jca, rca),
        rust_only_info_keys(jca, rca),
        java_only_info_keys(jca, rca)
    );
    kv(
        "closed_92307403",
        "2:92307403 C/A QUAL=154.64 DP=6 BaseQ=-1.834 ReadPos=1.282 MQRankSum=1.834",
    );

    let (j0, r0, vals) = first_info.expect("a remaining INFO split exists");
    kv("first_genuine_common", j0.site());
    kv("first_info_vals", format!("{vals:?}"));
    kv(
        "first_info_rust_only",
        format!("{:?}", rust_only_info_keys(&j0, &r0)),
    );
    kv(
        "first_info_java_only",
        format!("{:?}", java_only_info_keys(&j0, &r0)),
    );
    assert_eq!(
        (
            j0.chrom.as_str(),
            j0.pos,
            j0.ref_a.as_str(),
            j0.alt.as_str()
        ),
        FIRST_REMAINING,
        "first remaining common-site INFO split after 6R.239 G/A closure is 20:29455649 T/TGTTTG, got {}",
        j0.site()
    );
    assert!(
        vals.iter()
            .any(|(k, jv, rv)| k == "DP" && jv == "123" && rv == "230"),
        "INFO DP at 20:29455649 is 123 vs 230, got {vals:?}"
    );
    assert!(java_only_info_keys(&j0, &r0).is_empty());
    assert!(rust_only_info_keys(&j0, &r0).is_empty());
    let closed_ga = rec_key(CLOSED_GA);
    let jga = &java_map[&closed_ga];
    let rga = &rust_map[&closed_ga];
    assert!(format_equal(jga, rga), "6R.239 closed 20:29455379 FORMAT");
    let closed_pl = rec_key(CLOSED_CHR20_PL);
    let jc = &java_map[&closed_pl];
    let rc = &rust_map[&closed_pl];
    assert!(format_equal(jc, rc), "6R.218 closed 20:29455015 FORMAT");
    assert!(qual_close(&jc.qual, &rc.qual));
    assert!(
        genuine_info_value_diffs(jc, rc).is_empty(),
        "6R.218 closed 20:29455015 INFO; leftover={:?}",
        genuine_info_value_diffs(jc, rc)
    );
    kv(
        "primary_arrow_field",
        "6R.239 closed 20:29455379 G/A FORMAT; dump first remaining is PL/INFO at 20:29455649 — record only",
    );
    kv(
        "classification",
        "recorded only; 6R.240 rust-only covering emit is 20:29455314 G>C (not this round)",
    );
}

#[test]
fn forensic_6r197_target_live_emit() {
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing P12 ref/BAM");
        return;
    }
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");
    kv("closed_target", "2:92307403 C/A");

    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let gap_specs = parse_intervals_cli_string(&dict, GAP_INTERVAL).expect("gap");
    let tg_specs = parse_intervals_cli_string(&dict, TG_INTERVAL).expect("tg");
    let gap_walk = traverse_assembly_region_walker(
        &dict,
        &gap_specs,
        &ref_fasta,
        &bam,
        &ReadFilterParams::gatk_standard_hc(),
        &WalkerTraversalConfig::gatk_haplotype_caller_production(100),
    )
    .expect("gap walk");
    let tg_walk = traverse_assembly_region_walker(
        &dict,
        &tg_specs,
        &ref_fasta,
        &bam,
        &ReadFilterParams::gatk_standard_hc(),
        &WalkerTraversalConfig::gatk_haplotype_caller_production(100),
    )
    .expect("tg walk");
    let gap_regions = flatten_assembly_regions(&gap_walk);
    let tg_regions = flatten_assembly_regions(&tg_walk);
    let covering_ag = gap_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_AG.1
                && r.end.get() >= CLOSED_AG.1
        })
        .expect("covering A/G");
    let covering_for = |pos: u64| {
        tg_regions
            .iter()
            .find(|r| {
                matches!(
                    call_disposition(r),
                    AssemblyRegionCallDisposition::ActiveFull
                ) && r.start.get() <= pos
                    && r.end.get() >= pos
            })
            .unwrap_or_else(|| panic!("covering {pos}"))
    };
    let covering_qual = covering_for(CLOSED_QUAL.1);
    let covering_indel = covering_for(CLOSED_INDEL.1);
    let covering_tg = covering_for(CLOSED_CA.1);
    let ag_outcome = HaplotypeCallerEngine::call_region(
        covering_ag,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("ag call")
    .expect("ag outcome");
    let qual_outcome = HaplotypeCallerEngine::call_region(
        covering_qual,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("qual call")
    .expect("qual outcome");
    let indel_outcome = HaplotypeCallerEngine::call_region(
        covering_indel,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("indel call")
    .expect("indel outcome");
    let tg_outcome = HaplotypeCallerEngine::call_region(
        covering_tg,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("tg call")
    .expect("tg outcome");

    let closed_ag = ag_outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(CLOSED_AG.1)
                && c.event.ref_allele == CLOSED_AG.2
                && c.event.alt_allele == CLOSED_AG.3
        })
        .expect("A/G");
    let closed_qual = qual_outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(CLOSED_QUAL.1)
                && c.event.ref_allele == CLOSED_QUAL.2
                && c.event.alt_allele == CLOSED_QUAL.3
        })
        .expect("CT/C");
    let target = tg_outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(CLOSED_CA.1)
                && c.event.ref_allele == CLOSED_CA.2
                && c.event.alt_allele == CLOSED_CA.3
        })
        .expect("C/A");
    let closed_indel = indel_outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(CLOSED_INDEL.1)
                && c.event.ref_allele == CLOSED_INDEL.2
                && c.event.alt_allele == CLOSED_INDEL.3
        })
        .expect("TTC/T");

    let ag_n = unique_likelihood_indices(&closed_ag.annotation_likelihoods);
    let indel_n = unique_likelihood_indices(&closed_indel.annotation_likelihoods);
    let qual_n = unique_likelihood_indices(&closed_qual.annotation_likelihoods);
    let target_n = unique_likelihood_indices(&target.annotation_likelihoods);
    kv("closed_ag_ann_n", ag_n.len().to_string());
    kv("closed_indel_ann_n", indel_n.len().to_string());
    kv("closed_92307359_ann_n", qual_n.len().to_string());
    kv("target_ann_n", target_n.len().to_string());
    assert_eq!(ag_n.len(), 3, "6R.189 SiteScore n=3 stays");
    assert_eq!(indel_n.len(), 1, "6R.180 TTC/T n=1 stays");
    assert_eq!(
        qual_n.len(),
        0,
        "6R.192: CT/C has empty annotation_likelihoods and emit-fallback"
    );
    assert_eq!(target.genotype.format.ad_as_i32(), vec![2, 4]);
    assert_eq!(target.genotype.format.dp.as_i32(), 6);
    assert_eq!(target.genotype.format.pl_as_i32(), vec![162, 0, 72]);
    assert_eq!(
        target_n.len(),
        6,
        "6R.199: 2:92307403 attaches stored unique annotation_likelihoods n=6"
    );

    let ag_emitted = try_emit_call_region_variants(
        covering_ag,
        &ag_outcome,
        "NA12878",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .expect("ag emit");
    let qual_emitted = try_emit_call_region_variants(
        covering_qual,
        &qual_outcome,
        "NA12878",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .expect("qual emit");
    let tg_emitted = try_emit_call_region_variants(
        covering_tg,
        &tg_outcome,
        "NA12878",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .expect("tg emit");
    let ag_rec = ag_emitted
        .iter()
        .find(|r| r.position == CLOSED_AG.1 && r.reference == CLOSED_AG.2)
        .expect("A/G emit");
    let pin_rec = qual_emitted
        .iter()
        .find(|r| r.position == CLOSED_QUAL.1 && r.reference == CLOSED_QUAL.2)
        .expect("CT/C emit");
    let target_rec = tg_emitted
        .iter()
        .find(|r| r.position == CLOSED_CA.1 && r.reference == CLOSED_CA.2)
        .expect("C/A emit");
    assert_eq!(info_i32(&ag_rec.info, "DP"), Some(3));
    assert_eq!(info_i32(&pin_rec.info, "DP"), Some(2));
    assert!(info_has(&pin_rec.info, "BaseQRankSum"));
    assert!(info_has(&pin_rec.info, "MQRankSum"));
    assert!(info_has(&pin_rec.info, "ReadPosRankSum"));
    let pin_q = pin_rec.quality.expect("pin QUAL");
    assert!(
        (pin_q - 31.60).abs() < 0.005,
        "6R.196 QUAL stays 31.60, got {pin_q}"
    );
    assert_eq!(info_i32(&target_rec.info, "DP"), Some(6));
    let bq = info_f64(&target_rec.info, "BaseQRankSum").expect("6R.203 BaseQ");
    assert_eq!((bq * 1000.0).round(), -1834.0);
    let rprs = info_f64(&target_rec.info, "ReadPosRankSum").expect("6R.203 ReadPos");
    assert_eq!((rprs * 1000.0).round(), 1282.0);
    let mqrs = info_f64(&target_rec.info, "MQRankSum").expect("MQRankSum still emits");
    assert_eq!((mqrs * 1000.0).round(), 1834.0);
    let sample = target_rec.samples.first().expect("sample");
    assert_eq!(
        sample.gt.as_ref().map(|g| g.to_string()).as_deref(),
        Some("0/1")
    );
    assert_eq!(sample.ad.as_deref(), Some(&[2u32, 4][..]));
    kv("live_info_dp", "6");
    kv("live_java_info_dp", "6");
    kv(
        "next_arrow",
        "6R.203 closed BaseQRankSum/ReadPosRankSum at 2:92307403 via pre-realign getElementForRead",
    );
}
