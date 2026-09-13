//! 6R.173: live VCF reconnaissance vs frozen 6R.43 Java oracles (proof-only).
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! Rust VCFs: production `run_haplotype_caller` via `HOLDOUT_6R43=1`.
//! PRODUCTION CHANGE: NONE. 6R.174/176/178 closed `2:92305634 G/T`
//! (INFO DP, SOR, InbreedingCoeff omit). ReadPosRankSum at `2:92316347`
//! remains a known rust-only extra.
//!
//! ```text
//! HOLDOUT_6R43=1 cargo test -p gatk-haplotypecaller --test holdout_6r43_test -- --test-threads=1
//! cargo test -p gatk-haplotypecaller --test forensic_6r173_live_vcf_reconnaissance -- --nocapture --test-threads=1
//! HOLDOUT_6R173=1 cargo test -p gatk-haplotypecaller --test holdout_6r173_live_vcf_reconnaissance -- --nocapture --test-threads=1
//! ```

use std::collections::BTreeMap;
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

const KNOWN_RP: (&str, u64, &str, &str) = ("2", 92_316_347, "G", "A");
const CLOSED_6R178: (&str, u64, &str, &str) = ("2", 92_305_634, "G", "T");
/// First remaining INFO numeric split after 6R.174/176/178 closed 92305634.
const NEXT_INFO_SPLIT: (&str, u64, &str, &str) = ("2", 92_307_324, "TTC", "T");
const QUAL_TOL: f64 = 0.05;

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
    eprintln!("6R173\t{key}\t{}", value.as_ref());
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
        "FS" | "SOR" | "ExcessHet" => (j * 1000.0).round() == (r * 1000.0).round(),
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

fn unknown_rust_only_info(j: &Rec, r: &Rec) -> Vec<String> {
    r.info
        .keys()
        .filter(|k| !j.info.contains_key(*k) && *k != "ReadPosRankSum" && *k != "InbreedingCoeff")
        .cloned()
        .collect()
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

#[test]
fn forensic_6r173_live_vcf_reconnaissance() {
    let root = repo_root();
    let mut java_map = BTreeMap::new();
    let mut rust_map = BTreeMap::new();
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");
    kv("kbest", "legacy_1024 (production)");
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
    let mut chr2_fmt = 0usize;
    let mut ic_n = 0usize;
    let mut rp_rust_only_n = 0usize;
    let mut rp_java_only_n = 0usize;
    let mut bq_java_n = 0usize;
    let mut first_fmt: Option<(String, &'static str)> = None;
    for k in &common {
        let j = &java_map[k];
        let r = &rust_map[k];
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
        }
        if r.info.contains_key("InbreedingCoeff") && !j.info.contains_key("InbreedingCoeff") {
            ic_n += 1;
        }
        if r.info.contains_key("ReadPosRankSum") && !j.info.contains_key("ReadPosRankSum") {
            rp_rust_only_n += 1;
        }
        if j.info.contains_key("ReadPosRankSum") && !r.info.contains_key("ReadPosRankSum") {
            rp_java_only_n += 1;
        }
        if j.info.contains_key("BaseQRankSum") && !r.info.contains_key("BaseQRankSum") {
            bq_java_n += 1;
        }
    }

    kv("java_records", java_map.len().to_string());
    kv("rust_records", rust_map.len().to_string());
    kv("java_only", java_only.len().to_string());
    kv("rust_only", rust_only.len().to_string());
    kv("common", common.len().to_string());
    kv("common_format_diff", fmt_n.to_string());
    kv("common_qual_material", qual_n.to_string());
    kv("chr2_format_diff", chr2_fmt.to_string());
    kv("inbreeding_rust_only", ic_n.to_string());
    kv("readpos_rust_only", rp_rust_only_n.to_string());
    kv("readpos_java_only", rp_java_only_n.to_string());
    kv("baseq_java_only", bq_java_n.to_string());
    if let Some((site, fd)) = &first_fmt {
        kv("first_format_diff", format!("{site} field={fd}"));
    }
    if let Some(k) = rust_only.first() {
        kv("first_rust_only", rust_map[k].site());
    }

    let tgt = (
        KNOWN_RP.0.to_string(),
        KNOWN_RP.1,
        KNOWN_RP.2.to_string(),
        KNOWN_RP.3.to_string(),
    );
    let jt = &java_map[&tgt];
    let rt = &rust_map[&tgt];
    kv(
        "known_6r170_site",
        format!(
            "{} FORMAT_match={} QUAL_close={} rust_RP={:?} java_RP={:?} rust_IC={}",
            jt.site(),
            format_equal(jt, rt),
            qual_close(&jt.qual, &rt.qual),
            rt.info.get("ReadPosRankSum"),
            jt.info.get("ReadPosRankSum"),
            rt.info.contains_key("InbreedingCoeff")
        ),
    );
    assert!(
        format_equal(jt, rt),
        "6R.164 FORMAT at target must stay closed"
    );
    assert!(qual_close(&jt.qual, &rt.qual));
    assert_eq!(jt.gt, "1/1");
    assert_eq!(jt.ad, "0,3");
    assert_eq!(rt.info.get("ReadPosRankSum").map(String::as_str), Some("0"));
    assert!(!jt.info.contains_key("ReadPosRankSum"));
    kv(
        "known_6r170",
        "ReadPosRankSum rust=0 java=absent — not 6R.174",
    );
    kv(
        "known_6r169",
        "6R.178: InbreedingCoeff omitted when n<10 (rust.vcf refreshed by HOLDOUT_6R43)",
    );

    let closed_key = (
        CLOSED_6R178.0.to_string(),
        CLOSED_6R178.1,
        CLOSED_6R178.2.to_string(),
        CLOSED_6R178.3.to_string(),
    );
    let jc = &java_map[&closed_key];
    let rc = &rust_map[&closed_key];
    assert!(
        format_equal(jc, rc),
        "6R.164 FORMAT at 92305634 stays closed"
    );
    assert!(qual_close(&jc.qual, &rc.qual));
    assert_eq!(jc.info.get("DP").map(String::as_str), Some("3"));
    assert_eq!(rc.info.get("DP").map(String::as_str), Some("3"));
    assert_eq!(jc.info.get("SOR").map(String::as_str), Some("0.693"));
    assert!(
        !jc.info.contains_key("InbreedingCoeff") && !rc.info.contains_key("InbreedingCoeff"),
        "6R.178: both omit InbreedingCoeff at one-sample 92305634"
    );
    kv(
        "closed_6r178_target",
        "2:92305634 G/T FORMAT/QUAL/INFO DP/SOR/IC omit closed",
    );

    let mut first_genuine: Option<(Rec, Rec, Vec<(String, String, String)>)> = None;
    for k in &common {
        let j = &java_map[k];
        let r = &rust_map[k];
        if first_format_field(j, r).is_some() {
            continue;
        }
        if !qual_close(&j.qual, &r.qual) {
            continue;
        }
        let unk = unknown_rust_only_info(j, r);
        let val = genuine_info_value_diffs(j, r);
        if !unk.is_empty() {
            first_genuine = Some((j.clone(), r.clone(), val));
            kv("first_genuine_kind", "INFO_KEY_UNKNOWN");
            break;
        }
        if !val.is_empty() {
            kv("first_genuine_kind", format!("INFO_VALUE {}", val[0].0));
            first_genuine = Some((j.clone(), r.clone(), val));
            break;
        }
    }
    let (j0, r0, vals) = first_genuine.expect("expected a remaining INFO numeric split");
    kv("first_genuine", j0.site());
    kv(
        "first_genuine_format",
        format!(
            "GT={} AD={} DP={} GQ={} PL={} QUAL j={} r={}",
            j0.gt, j0.ad, j0.dp, j0.gq, j0.pl, j0.qual, r0.qual
        ),
    );
    kv("first_genuine_info_vals", format!("{vals:?}"));
    assert_eq!(j0.chrom, NEXT_INFO_SPLIT.0);
    assert_eq!(j0.pos, NEXT_INFO_SPLIT.1);
    assert_eq!(j0.ref_a, NEXT_INFO_SPLIT.2);
    assert_eq!(j0.alt, NEXT_INFO_SPLIT.3);
    assert!(format_equal(&j0, &r0));
    assert_eq!(vals[0].0, "DP");
    kv(
        "first_arrow_hypothesis",
        "2:92305634 closed by 6R.174/176/178. Remaining INFO DP/MQ/SOR at 2:92307324 TTC/T is later stacked, not 6R.178.",
    );

    let ann = include_str!("../src/variant_site_hc_annotations.rs");
    assert!(
        ann.contains("coverage_evidence_count") && ann.contains("Coverage.annotate"),
        "6R.174 INFO DP uses Java Coverage.evidenceCount"
    );
    assert!(
        !ann.contains("92305634") && !ann.contains("92316347"),
        "no coordinate-specific DP/RankSum rule"
    );

    assert!(java_only.is_empty());
    assert_eq!(rust_only.len(), 3);
    assert_eq!(chr2_fmt, 0);
    assert_eq!(ic_n, 0, "6R.178: rust.vcf must omit InbreedingCoeff at n=1");
    assert_eq!(java_map.len(), 126);
    assert_eq!(rust_map.len(), 129);
}
