//! 6R.192: proof-only fresh INFO/VCF reconnaissance after 6R.191.
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! PRODUCTION CHANGE: NONE.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r192_fresh_info_reconnaissance -- --nocapture --test-threads=1
//! HOLDOUT_6R192=1 cargo test -p gatk-haplotypecaller --test holdout_6r192_fresh_info -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::variant_site_hc_annotations::{
    read_pos_rank_sum_quals_from_likelihoods, HcStrandBiasLikelihoods,
};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    HaplotypeCallerEngine, HcGenotypingConfig, ReadFilterParams, RegionReadLikelihood,
    WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
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
const FIRST_REMAINING: (&str, u64, &str, &str) = ("2", 92_307_359, "CT", "C");
const CANDIDATE_BAIT: (&str, u64) = ("2", 92_305_759);
const CLOSED_GT: (&str, u64, &str, &str) = ("2", 92_305_634, "G", "T");
const CLOSED_AG: (&str, u64, &str, &str) = ("2", 92_305_635, "A", "G");
const CLOSED_AC: (&str, u64, &str, &str) = ("2", 92_305_716, "A", "C");
const CLOSED_INDEL: (&str, u64, &str, &str) = ("2", 92_307_324, "TTC", "T");
const CLOSED_TG: (&str, u64, &str, &str) = ("2", 92_307_333, "T", "G");
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
    eprintln!("6R192\t{key}\t{}", value.as_ref());
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

fn unique_likelihood_indices(likelihoods: &[RegionReadLikelihood]) -> BTreeSet<usize> {
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
fn forensic_6r192_no_production_change() {
    kv("java_pin", JAVA_PIN);
    kv(
        "production_change",
        "NONE at 6R.192; BaseQRankSum emit closed by 6R.193; MQRankSum closed by 6R.194",
    );
    kv(
        "fresh_pass",
        "6R.192 proved Java-only BaseQRankSum/MQRankSum at 2:92307359. 6R.193 closed BaseQRankSum. 6R.194 closed MQRankSum.",
    );
    kv(
        "classification",
        "F — WRONG EMISSION PREDICATE (BaseQ closed by 6R.193; MQ closed by 6R.194)",
    );

    let emit = include_str!("../src/region_vcf_emit.rs");
    let ann = include_str!("../src/variant_site_hc_annotations.rs");
    let pipe = include_str!("../src/hc_genotyping_engine/genotype_site_pipeline.rs");
    let early = include_str!("../src/hc_genotyping_engine/genotype_site_early_template.rs");
    let info_fn = emit
        .split("fn hc_info_values")
        .nth(1)
        .expect("hc_info_values");
    let info_body = info_fn
        .split("fn try_emit_call_region_variants")
        .next()
        .expect("body");
    assert!(
        info_body.contains("BaseQRankSum") && info_body.contains("InfoValue::Float(\"MQRankSum\""),
        "6R.193 emits BaseQRankSum; 6R.194 emits MQRankSum"
    );
    assert!(
        emit.contains("\"BaseQRankSum\"") && emit.contains("\"MQRankSum\""),
        "header still declares BaseQRankSum/MQRankSum; declaration is not emission"
    );
    assert!(
        ann.contains("read_pos_rank_sum_quals_from_likelihoods")
            && pipe.contains("let annotation = annotation_likelihoods_from_stored_haplotypes")
            && early.contains("if is_cluster_tg_snp(&event)"),
        "6R.191 ReadPos bind, 6R.189 SiteScore attach, 6R.186 cluster-TG stay"
    );
    assert!(
        !ann.contains("92307359")
            && !emit.contains("92307359")
            && !pipe.contains("92307359")
            && !ann.contains("92305759"),
        "no locus-specific 6R.192 production patch"
    );
}

#[test]
fn forensic_6r192_live_vcf_inventory() {
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
    kv(
        "production_change",
        "NONE at 6R.192; BaseQRankSum emit closed by 6R.193; MQRankSum closed by 6R.194",
    );
    kv(
        "fresh_pass",
        "post-6R.191 production rust.vcf vs frozen Java 6R.43/6R.160 corpus",
    );
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
    assert_eq!(rust_only.len(), 3);
    assert_eq!(chr2_fmt, 0, "chr2 FORMAT must stay closed after 6R.191");
    assert_eq!(chr2_qual, 0, "chr2 QUAL must stay closed within 0.05");
    assert_eq!(java_map.len(), 126);
    assert_eq!(rust_map.len(), 129);
    assert!(
        first_fmt
            .as_ref()
            .is_some_and(|(site, fd)| site.starts_with("20:29455379") && *fd == "AD"),
        "first FORMAT split after 6R.218 is AD at 20:29455379: {first_fmt:?}"
    );
    assert!(
        extra_info_keys.is_empty(),
        "no unexpected INFO keys beyond the compared set: {extra_info_keys:?}"
    );
    assert!(
        !java_map
            .keys()
            .any(|k| k.0 == CANDIDATE_BAIT.0 && k.1 == CANDIDATE_BAIT.1)
            && !rust_map
                .keys()
                .any(|k| k.0 == CANDIDATE_BAIT.0 && k.1 == CANDIDATE_BAIT.1),
        "bait locus 2:92305759 is not in this corpus"
    );

    for (label, spec, dp, mq, sor) in [
        ("closed_92305634", CLOSED_GT, "3", 41.96, 0.693),
        ("closed_92305635", CLOSED_AG, "3", 41.96, 0.693),
        ("closed_92305716", CLOSED_AC, "4", 35.15, 1.179),
        ("closed_92307324", CLOSED_INDEL, "1", 44.0, 1.609),
        ("closed_92307333", CLOSED_TG, "1", 44.0, 1.609),
        ("closed_92316347", CLOSED_MID, "3", 40.25, 1.179),
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
        assert!(!j.info.contains_key("BaseQRankSum") && !r.info.contains_key("BaseQRankSum"));
        assert!(!j.info.contains_key("MQRankSum") && !r.info.contains_key("MQRankSum"));
        kv(label, format!("{} FORMAT/QUAL/DP/MQ/SOR closed", j.site()));
    }
    let ag = &rust_map[&rec_key(CLOSED_AG)];
    assert_eq!(ag.gt, "1/1");
    assert_eq!(ag.ad, "0,2");
    assert_eq!(ag.dp, "2");
    assert_eq!(ag.gq, "6");
    assert_eq!(ag.pl, "90,6,0");
    assert_eq!(ag.qual, "78.32");
    let ac = &rust_map[&rec_key(CLOSED_AC)];
    let jac = &java_map[&rec_key(CLOSED_AC)];
    assert_eq!(ac.gt, "1/1");
    assert_eq!(ac.ad, "0,3");
    assert_eq!(ac.dp, "3");
    assert_eq!(ac.gq, "9");
    assert_eq!(ac.pl, "130,9,0");
    assert_eq!(ac.qual, "116.84");
    assert!(!jac.info.contains_key("ReadPosRankSum") && !ac.info.contains_key("ReadPosRankSum"));

    let pin = rec_key(FIRST_REMAINING);
    let jp = &java_map[&pin];
    let rp = &rust_map[&pin];
    assert!(format_equal(jp, rp), "6R.192 pin FORMAT stays closed");
    assert!(
        qual_close(&jp.qual, &rp.qual),
        "6R.196 closed QUAL at the 6R.192 pin"
    );
    assert_eq!(jp.qual, "31.60");
    assert_eq!(rp.qual, "31.60");
    assert!(
        genuine_info_value_diffs(jp, rp).is_empty()
            && rust_only_info_keys(jp, rp).is_empty()
            && java_only_info_keys(jp, rp).is_empty(),
        "6R.193–6R.196 closed the 6R.192 pin 2:92307359 CT/C; leftover={:?} rust_only={:?} java_only={:?}",
        genuine_info_value_diffs(jp, rp),
        rust_only_info_keys(jp, rp),
        java_only_info_keys(jp, rp)
    );
    assert_eq!(jp.info.get("QD").map(String::as_str), Some("15.80"));
    assert!(info_print_close(
        "QD",
        jp.info.get("QD").expect("java qd"),
        rp.info.get("QD").expect("rust qd"),
    ));
    assert!(jp.info.contains_key("BaseQRankSum") && rp.info.contains_key("BaseQRankSum"));
    assert!(jp.info.contains_key("MQRankSum") && rp.info.contains_key("MQRankSum"));
    kv(
        "closed_92307359",
        "2:92307359 CT/C FORMAT/QUAL/QD/RankSums closed by 6R.193–6R.196",
    );

    if let Some((j0, r0, vals)) = first_info {
        kv("first_info_after_pin", j0.site());
        kv("first_info_vals", format!("{vals:?}"));
        kv(
            "first_info_rust_only",
            format!("{:?}", rust_only_info_keys(&j0, &r0)),
        );
        kv(
            "first_info_java_only",
            format!("{:?}", java_only_info_keys(&j0, &r0)),
        );
        assert_ne!(
            (
                j0.chrom.as_str(),
                j0.pos,
                j0.ref_a.as_str(),
                j0.alt.as_str()
            ),
            FIRST_REMAINING,
            "6R.192 pin must not remain the first INFO split"
        );
    }
    kv(
        "primary_arrow_field",
        "6R.192 pin 2:92307359 closed by 6R.193/6R.194 RankSum emit, 6R.195 QD-follows-QUAL, 6R.196 indel AF prior",
    );
    kv(
        "classification",
        "F — WRONG EMISSION PREDICATE (BaseQ/MQ closed); QUAL prior closed by 6R.196",
    );
}

#[test]
fn forensic_6r192_first_arrow_emission_predicate() {
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing P12 ref/BAM");
        return;
    }
    kv("java_pin", JAVA_PIN);
    kv(
        "production_change",
        "NONE at 6R.192; BaseQRankSum emit closed by 6R.193; MQRankSum closed by 6R.194",
    );
    kv("target", "2:92307359 CT/C");
    kv(
        "java_path",
        "StandardAnnotation → BaseQualityRankSumTest/MappingQualityRankSumTest → RankSumTest.fillQualsFromLikelihood(annotation AlleleLikelihoods) → finite z including 0.0 → INFO insert",
    );
    kv(
        "rust_path",
        "6R.194: hc_info_values inserts MQRankSum from the same fillQuals object as ReadPos/BaseQ.",
    );

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
    let covering_ac = gap_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_AC.1
                && r.end.get() >= CLOSED_AC.1
        })
        .expect("covering A/C");
    let covering_tg = tg_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= FIRST_REMAINING.1
                && r.end.get() >= FIRST_REMAINING.1
        })
        .expect("tg covering");
    let ag_outcome = HaplotypeCallerEngine::call_region(
        covering_ag,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("ag call")
    .expect("ag outcome");
    let ac_outcome = HaplotypeCallerEngine::call_region(
        covering_ac,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("ac call")
    .expect("ac outcome");
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
    let closed_ac = ac_outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(CLOSED_AC.1)
                && c.event.ref_allele == CLOSED_AC.2
                && c.event.alt_allele == CLOSED_AC.3
        })
        .expect("A/C");
    let closed_tg = tg_outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(CLOSED_TG.1)
                && c.event.ref_allele == CLOSED_TG.2
                && c.event.alt_allele == CLOSED_TG.3
        })
        .expect("T/G");
    let target = tg_outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(FIRST_REMAINING.1)
                && c.event.ref_allele == FIRST_REMAINING.2
                && c.event.alt_allele == FIRST_REMAINING.3
        })
        .expect("CT/C");

    let ag_n = unique_likelihood_indices(&closed_ag.annotation_likelihoods).len();
    let ac_n = unique_likelihood_indices(&closed_ac.annotation_likelihoods).len();
    let tg_n = unique_likelihood_indices(&closed_tg.annotation_likelihoods).len();
    let target_n = unique_likelihood_indices(&target.annotation_likelihoods).len();
    kv("closed_92305635_annotation_n", ag_n.to_string());
    kv("closed_92305716_annotation_n", ac_n.to_string());
    kv("closed_92307333_annotation_n", tg_n.to_string());
    kv("target_annotation_n", target_n.to_string());
    assert_eq!(ag_n, 3, "6R.189 n=3 must not change");
    assert_eq!(ac_n, 4, "6R.191 annotation object must not change");
    assert_eq!(tg_n, 1, "6R.186 n=1 must not change");
    assert_eq!(
        target_n, 0,
        "this indel does not take the 6R.189 SNP loc-loop attach; emit falls back (DP/MQ/SOR/ReadPos still match)"
    );
    assert_eq!(closed_ag.genotype.format.ad_as_i32(), vec![0, 2]);
    assert_eq!(closed_ac.genotype.format.ad_as_i32(), vec![0, 3]);
    assert_eq!(target.genotype.format.ad_as_i32(), vec![1, 1]);
    assert_eq!(target.genotype.format.dp.as_i32(), 2);

    let (full_ref, full_pad) = tg_outcome.assembly.event_map_reference();
    let apply_pad = tg_outcome
        .assembly
        .haplotypes
        .iter()
        .find(|h| h.is_reference)
        .and_then(|h| h.genome_loc.as_ref().map(|g| g.start.get()))
        .unwrap_or(full_pad);
    let cfg = HcGenotypingConfig::default();
    let likelihoods = if target.annotation_likelihoods.is_empty() {
        tg_outcome.read_likelihoods.as_slice()
    } else {
        target.annotation_likelihoods.as_slice()
    };
    let evidence = HcStrandBiasLikelihoods {
        reads: &tg_outcome.genotyping_reads,
        likelihoods,
        haplotypes: &tg_outcome.assembly.haplotypes,
        contig: &covering_tg.contig,
        ref_bytes: tg_outcome.assembly.reference_bases(),
        pad_start_1based: apply_pad,
        full_ref_bytes: full_ref,
        full_pad_1based: full_pad,
        max_mnp_distance: tg_outcome.assembly.max_mnp_distance(),
        emit_spanning_dels: !cfg.disable_spanning_event_genotyping,
    };
    let (refs, alts) = read_pos_rank_sum_quals_from_likelihoods(
        &evidence,
        FIRST_REMAINING.1,
        FIRST_REMAINING.2,
        FIRST_REMAINING.3,
    );
    kv(
        "ranksum_lists",
        format!("REF={} ALT={}", refs.len(), alts.len()),
    );
    assert!(
        !refs.is_empty() && !alts.is_empty(),
        "Java-equivalent RankSum membership on the emit object (empty-annotation fallback) is defined on this het"
    );

    let ag_emitted = try_emit_call_region_variants(
        covering_ag,
        &ag_outcome,
        "NA12878",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .expect("ag emit");
    let ac_emitted = try_emit_call_region_variants(
        covering_ac,
        &ac_outcome,
        "NA12878",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .expect("ac emit");
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
    let ac_rec = ac_emitted
        .iter()
        .find(|r| r.position == CLOSED_AC.1 && r.reference == CLOSED_AC.2)
        .expect("A/C emit");
    let tg_rec = tg_emitted
        .iter()
        .find(|r| r.position == CLOSED_TG.1 && r.reference == CLOSED_TG.2)
        .expect("T/G emit");
    let target_rec = tg_emitted
        .iter()
        .find(|r| r.position == FIRST_REMAINING.1 && r.reference == FIRST_REMAINING.2)
        .expect("CT/C emit");
    assert_eq!(info_i32(&ag_rec.info, "DP"), Some(3));
    assert!((info_f64(&ag_rec.info, "MQ").unwrap_or(-1.0) - 41.96).abs() < 0.005);
    assert!(!info_has(&ac_rec.info, "ReadPosRankSum"));
    assert_eq!(info_i32(&tg_rec.info, "DP"), Some(1));
    let sample = target_rec.samples.first().expect("sample");
    assert_eq!(
        sample.gt.as_ref().map(|g| g.to_string()).as_deref(),
        Some("0/1")
    );
    assert_eq!(sample.ad.as_deref(), Some(&[1u32, 1][..]));
    assert_eq!(sample.dp.map(|d| d as i32), Some(2));
    assert!((target_rec.quality.unwrap_or(0.0) - 31.64).abs() < 0.05);
    assert!(
        info_has(&target_rec.info, "ReadPosRankSum"),
        "ReadPosRankSum is defined and already emitted; this is not another ReadPos source bug"
    );
    assert!(
        info_has(&target_rec.info, "BaseQRankSum"),
        "6R.193 closed BaseQRankSum emit at the 6R.192 first remaining site"
    );
    let bq = info_f64(&target_rec.info, "BaseQRankSum").unwrap_or(-1.0);
    assert!(bq.abs() < 0.0005, "finite zero must emit, got {bq}");
    let mqrs = info_f64(&target_rec.info, "MQRankSum").expect("6R.194 MQRankSum");
    assert_eq!(
        (mqrs * 1000.0).round(),
        -674.0,
        "6R.194 closed MQRankSum emit at the 6R.192 first remaining site, got {mqrs}"
    );
    kv(
        "first_divergent_arrow",
        "6R.192 pin 2:92307359 closed by 6R.193/6R.194 RankSum emit and 6R.196 QUAL prior",
    );
    kv(
        "classification",
        "F — WRONG EMISSION PREDICATE (BaseQ closed by 6R.193; MQ closed by 6R.194)",
    );
}
