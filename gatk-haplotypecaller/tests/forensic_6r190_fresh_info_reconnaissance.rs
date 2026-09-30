//! 6R.190: proof-only fresh INFO/VCF reconnaissance after 6R.189.
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! PRODUCTION CHANGE: NONE.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r190_fresh_info_reconnaissance -- --nocapture --test-threads=1
//! HOLDOUT_6R190=1 cargo test -p gatk-haplotypecaller --test holdout_6r190_fresh_info -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::{build_per_haplotype_variation_events, VariationEvent};
use gatk_haplotypecaller::fragment_overlap::read_base_at_ref_coord_1based;
use gatk_haplotypecaller::hc_allele_mapping::create_allele_mapper_with_events;
use gatk_haplotypecaller::read_model::MAPPING_QUALITY_UNAVAILABLE;
use gatk_haplotypecaller::read_realignment::LOG_10_INFORMATIVE_THRESHOLD;
use gatk_haplotypecaller::variant_site_hc_annotations::HcStrandBiasLikelihoods;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, marginalize_rows_to_biallelic_alleles,
    region_likelihoods_to_rows, traverse_assembly_region_walker, try_emit_call_region_variants,
    AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition, HaplotypeCallerEngine,
    HcGenotypingConfig, ReadFilterParams, RegionReadLikelihood, WalkerTraversalConfig,
    DEFAULT_STAND_EMIT_CONFIDENCE,
};
use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::Record;
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
const FIRST_REMAINING: (&str, u64, &str, &str) = ("2", 92_305_716, "A", "C");
const CLOSED_GT: (&str, u64, &str, &str) = ("2", 92_305_634, "G", "T");
const CLOSED_AG: (&str, u64, &str, &str) = ("2", 92_305_635, "A", "G");
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
const FLAG_REVERSE: u16 = 0x10;

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
    eprintln!("6R190\t{key}\t{}", value.as_ref());
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

fn classify_info_key(j: &Rec, r: &Rec, key: &str) -> &'static str {
    match (j.info.get(key), r.info.get(key)) {
        (None, None) => "absent both (header is not emission)",
        (None, Some(_)) => "Rust-only key",
        (Some(_), None) => "Java-only key",
        (Some(jv), Some(rv)) if jv == rv => "same key / same value",
        (Some(jv), Some(rv)) if info_print_close(key, jv, rv) => "formatting-only difference",
        (Some(_), Some(_)) => "same key / different value",
    }
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

fn qname(rec: &Record) -> String {
    String::from_utf8_lossy(rec.qname()).into_owned()
}

fn strand_label(rec: &Record) -> &'static str {
    if rec.flags() & FLAG_REVERSE != 0 {
        "rev"
    } else {
        "fwd"
    }
}

/// Java `RankSumTest.isUsableRead`: MQ != 0 and MQ != 255.
fn java_ranksum_mq_usable(rec: &Record) -> bool {
    let mq = rec.mapq();
    mq != 0 && mq != MAPPING_QUALITY_UNAVAILABLE
}

/// Java `ReadPosRankSumTest.getReadPosition`: `min(left, right)` including hard clips.
fn java_read_pos_rank_element(rec: &Record, vc_start_1based: i32) -> Option<f64> {
    let cig = rec.cigar();
    let mut leading_hard = 0i64;
    let mut trailing_hard = 0i64;
    if let Some(Cigar::HardClip(n)) = cig.iter().next() {
        leading_hard = i64::from(*n);
    }
    if let Some(Cigar::HardClip(n)) = cig.iter().last() {
        trailing_hard = i64::from(*n);
    }
    let vc = i64::from(vc_start_1based);
    let mut ref_pos = rec.pos() + 1;
    let mut q = 0usize;
    for c in cig.iter() {
        match c {
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                let n = i64::from(*n);
                if vc >= ref_pos && vc < ref_pos + n {
                    let idx = q + (vc - ref_pos) as usize;
                    let seq_len = rec.seq_len() as i64;
                    let left = leading_hard + idx as i64;
                    let right = seq_len - 1 - idx as i64 + trailing_hard;
                    return Some(left.min(right) as f64);
                }
                ref_pos += n;
                q += n as usize;
            }
            Cigar::Del(n) | Cigar::RefSkip(n) => {
                ref_pos += i64::from(*n);
            }
            Cigar::Ins(n) | Cigar::SoftClip(n) => {
                q += *n as usize;
            }
            Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }
    None
}

#[test]
fn forensic_6r190_no_production_change() {
    kv("java_pin", JAVA_PIN);
    kv(
        "production_change",
        "NONE at 6R.190; RankSum source closed by 6R.191",
    );
    kv(
        "fresh_pass",
        "6R.190 is a fresh post-6R.189 reconnaissance pass. 6R.191 closed the proven pileup source.",
    );
    kv(
        "classification",
        "A — WRONG SOURCE OBJECT (closed by 6R.191)",
    );

    let early = include_str!("../src/hc_genotyping_engine/genotype_site_early_template.rs");
    let pipe = include_str!("../src/hc_genotyping_engine/genotype_site_pipeline.rs");
    let emit = include_str!("../src/region_vcf_emit.rs");
    let ann = include_str!("../src/variant_site_hc_annotations.rs");
    let engine = include_str!("../src/engine.rs");
    let rp = include_str!("../src/annotator/plugins/read_pos_rank_sum.rs");

    assert!(
        pipe.contains("let annotation = annotation_likelihoods_from_stored_haplotypes")
            && pipe.contains("with_annotation_likelihoods(annotation)"),
        "6R.189 SiteScore stored-hap attach must stay"
    );
    assert!(
        early.contains("fn annotation_likelihoods_from_stored_haplotypes")
            && early.contains("if is_cluster_tg_snp(&event)"),
        "6R.186 cluster-TG helper must stay"
    );
    assert!(
        ann.contains("read_pos_rank_sum_quals_from_likelihoods")
            && ann.contains(
                "let rp = read_pos_rank_sum::read_pos_rank_sum(&ref_positions, &alt_positions)"
            ),
        "6R.191: ReadPosRankSum lists come from annotation likelihoods; formula stays 6R.172"
    );
    assert!(
        ann.contains("fn read_offset_evidence_at_site"),
        "6R.170 pileup helper remains as retired code, not the live bind"
    );
    assert!(
        emit.contains("if let Some(z) = ann.read_pos_rank_sum"),
        "6R.172 Some-only INFO insert stays"
    );
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
        rp.contains("if alt.is_empty() || reference.is_empty()") && rp.contains("return None"),
        "empty RankSum lists remain None (Java emptyMap)"
    );
    assert!(
        engine.contains("fn apply_java_order_normalize_and_filter"),
        "6R.185 stored-matrix lifecycle stays closed"
    );
    assert!(
        ann.contains("6R.174: INFO DP is Java `Coverage.evidenceCount`")
            && ann.contains("mq_rms_of_sample_evidence")
            && ann.contains("strand_bias_sample_counts"),
        "DP/MQ/SOR formulas must stay closed"
    );
    assert!(
        !ann.contains("92305716") && !emit.contains("92305716") && !pipe.contains("92305716"),
        "no locus-specific 6R.190 production patch in the annotation/emit path"
    );
}

#[test]
fn forensic_6r190_live_vcf_inventory() {
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
        "post-6R.189 production rust.vcf vs frozen Java 6R.43/6R.160 corpus",
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
    assert_eq!(rust_only.len(), 4);
    assert_eq!(
        rust_map[rust_only.first().unwrap()].site(),
        "20:29455314 G/C"
    );
    assert_eq!(chr2_fmt, 0, "chr2 FORMAT must stay closed after 6R.189");
    assert_eq!(chr2_qual, 0, "chr2 QUAL must stay closed within 0.05");
    assert_eq!(java_map.len(), 126);
    assert_eq!(rust_map.len(), 130);
    assert!(
        first_fmt
            .as_ref()
            .is_some_and(|(site, fd)| site.starts_with("20:29455649") && *fd == "PL"),
        "6R.239 closed AD at 20:29455379; first remaining FORMAT split is PL at 20:29455649: {first_fmt:?}"
    );

    kv(
        "6r187_stale_frontier",
        "6R.187 first remaining was 2:92305635 INFO DP vs then-current rust.vcf. 6R.189 closed that site. This pass does not skip rust-only ReadPosRankSum.",
    );
    for (label, spec, dp, mq, sor) in [
        ("closed_92305634", CLOSED_GT, "3", 41.96, 0.693),
        ("closed_92305635", CLOSED_AG, "3", 41.96, 0.693),
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
        kv(label, format!("{} FORMAT/QUAL/DP/MQ/SOR closed", j.site()));
    }
    let ag = &rust_map[&rec_key(CLOSED_AG)];
    let jag = &java_map[&rec_key(CLOSED_AG)];
    assert_eq!(ag.gt, "1/1");
    assert_eq!(ag.ad, "0,2");
    assert_eq!(ag.dp, "2");
    assert_eq!(ag.gq, "6");
    assert_eq!(ag.pl, "90,6,0");
    assert_eq!(ag.qual, "78.32");
    assert_eq!(jag.qual, "78.32");
    assert!(!jag.info.contains_key("ReadPosRankSum"));
    kv(
        "closed_92305635_format",
        "GT=1/1 AD=0,2 DP=2 GQ=6 PL=90,6,0 QUAL=78.32 INFO DP=3 MQ=41.96 SOR=0.693",
    );

    let site716 = rec_key(FIRST_REMAINING);
    let j716 = &java_map[&site716];
    let r716 = &rust_map[&site716];
    assert!(
        format_equal(j716, r716),
        "6R.190 target FORMAT stays closed"
    );
    assert!(qual_close(&j716.qual, &r716.qual));
    assert!(
        !j716.info.contains_key("ReadPosRankSum") && !r716.info.contains_key("ReadPosRankSum"),
        "6R.191 closed rust-only ReadPosRankSum at 2:92305716"
    );
    kv(
        "closed_92305716_by_6r191",
        "FORMAT/QUAL/DP stay; ReadPosRankSum absent both",
    );

    let (j0, r0, vals) = first_info.expect("expected a remaining INFO split after 6R.191");
    kv("post_6r191_first_info", j0.site());
    kv("post_6r191_vals", format!("{vals:?}"));
    kv(
        "post_6r191_rust_only",
        format!("{:?}", rust_only_info_keys(&j0, &r0)),
    );
    kv(
        "post_6r191_java_only",
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
        "6R.191 closed 2:92305716; this inventory must not still treat it as first remaining"
    );
    kv(
        "primary_arrow_field",
        "6R.190 rust-only ReadPosRankSum at 2:92305716 closed by 6R.191",
    );
}

#[test]
fn forensic_6r190_read_pos_rank_sum_source_arrow() {
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing P12 ref/BAM");
        return;
    }
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");
    kv("target", "2:92305716 A/C");
    kv(
        "java_path",
        "StandardAnnotation → ReadPosRankSumTest → RankSumTest.fillQualsFromLikelihood(annotation AlleleLikelihoods) → bestAllelesBreakingTies + isInformative + isUsableRead → empty REF on hom-alt → MannWhitneyU NaN → emptyMap",
    );
    kv(
        "rust_path",
        "pre-6R.191: annotate_hc_variant_site → read_offset_evidence_at_site(region.reads pileup). 6R.191 binds fillQualsFromLikelihood(annotation_likelihoods).",
    );

    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, GAP_INTERVAL).expect("interval");
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
    let covering_closed = regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_AG.1
                && r.end.get() >= CLOSED_AG.1
        })
        .expect("covering closed A/G");
    let covering = regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= FIRST_REMAINING.1
                && r.end.get() >= FIRST_REMAINING.1
        })
        .expect("covering");
    let closed_outcome = HaplotypeCallerEngine::call_region(
        covering_closed,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("closed call")
    .expect("closed outcome");
    let outcome = HaplotypeCallerEngine::call_region(
        covering,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("outcome");

    let closed_ag = closed_outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(CLOSED_AG.1)
                && c.event.ref_allele == CLOSED_AG.2
                && c.event.alt_allele == CLOSED_AG.3
        })
        .expect("closed A/G");
    let target = outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(FIRST_REMAINING.1)
                && c.event.ref_allele == FIRST_REMAINING.2
                && c.event.alt_allele == FIRST_REMAINING.3
        })
        .expect("target A/C");
    let closed_ann = unique_likelihood_indices(&closed_ag.annotation_likelihoods);
    let target_ann = unique_likelihood_indices(&target.annotation_likelihoods);
    kv("closed_92305635_annotation_n", closed_ann.len().to_string());
    kv("target_92305716_annotation_n", target_ann.len().to_string());
    assert_eq!(closed_ann.len(), 3, "6R.189 loc-loop n=3 must not regress");
    assert_eq!(closed_ag.genotype.format.ad_as_i32(), vec![0, 2]);
    assert_eq!(closed_ag.genotype.format.dp.as_i32(), 2);
    assert_eq!(closed_ag.genotype.format.pl_as_i32(), vec![90, 6, 0]);
    assert_eq!(target.genotype.format.ad_as_i32(), vec![0, 3]);
    assert_eq!(target.genotype.format.dp.as_i32(), 3);
    assert_eq!(target.genotype.format.pl_as_i32(), vec![130, 9, 0]);
    assert!(
        !target.annotation_likelihoods.is_empty(),
        "6R.189 SiteScore attach is present; RankSum still ignores it"
    );

    let (full_ref, full_pad) = outcome.assembly.event_map_reference();
    let apply_pad = outcome
        .assembly
        .haplotypes
        .iter()
        .find(|h| h.is_reference)
        .and_then(|h| h.genome_loc.as_ref().map(|g| g.start.get()))
        .unwrap_or(full_pad);
    let cfg = HcGenotypingConfig::default();
    let evidence = HcStrandBiasLikelihoods {
        reads: &outcome.genotyping_reads,
        likelihoods: &target.annotation_likelihoods,
        haplotypes: &outcome.assembly.haplotypes,
        contig: &covering.contig,
        ref_bytes: outcome.assembly.reference_bases(),
        pad_start_1based: apply_pad,
        full_ref_bytes: full_ref,
        full_pad_1based: full_pad,
        max_mnp_distance: outcome.assembly.max_mnp_distance(),
        emit_spanning_dels: !cfg.disable_spanning_event_genotyping,
    };
    kv(
        "java_ranksum_source_object",
        "attached annotation_likelihoods (loc-loop stored-hap after retainEvidence), independently of assuming it is correct for every annotator",
    );
    kv(
        "rust_ranksum_source_object",
        "6R.190 measured covering.reads pileup vs annotation_likelihoods; 6R.191 binds the latter",
    );
    let annotation_ids: BTreeSet<String> = target_ann
        .iter()
        .filter_map(|idx| evidence.reads.get(*idx))
        .map(|rec| format!("{} FLAG={}", qname(rec), rec.flags()))
        .collect();

    let hap_cache = build_per_haplotype_variation_events(
        evidence.haplotypes,
        evidence.full_ref_bytes,
        evidence.full_pad_1based,
        evidence.max_mnp_distance,
        evidence.contig,
    );
    let mapping = create_allele_mapper_with_events(
        &VariationEvent::from_alleles(
            evidence.contig,
            FIRST_REMAINING.1,
            FIRST_REMAINING.2,
            FIRST_REMAINING.3,
        ),
        FIRST_REMAINING.1,
        evidence.haplotypes,
        evidence.pad_start_1based,
        evidence.ref_bytes,
        evidence.max_mnp_distance,
        evidence.emit_spanning_dels,
        Some(&hap_cache),
    );
    let rows = region_likelihoods_to_rows(evidence.likelihoods, evidence.haplotypes.len());
    let marg = marginalize_rows_to_biallelic_alleles(
        &rows,
        &mapping.ref_haplotype_indices,
        &mapping.alt_haplotype_indices,
    );

    let mut java_ref = Vec::new();
    let mut java_alt = Vec::new();
    for row in &marg {
        let Some(rec) = evidence.reads.get(row.read_index) else {
            continue;
        };
        let qname = qname(rec);
        let lls = &row.haplotype_log10_likelihoods;
        let ll_ref = lls.first().copied().unwrap_or(f64::NEG_INFINITY);
        let ll_alt = lls.get(1).copied().unwrap_or(f64::NEG_INFINITY);
        if !ll_ref.is_finite() && !ll_alt.is_finite() {
            continue;
        }
        let (best_is_ref, best, second) = if ll_ref > ll_alt {
            (true, ll_ref, ll_alt)
        } else if ll_alt > ll_ref {
            (false, ll_alt, ll_ref)
        } else {
            (true, ll_ref, ll_alt)
        };
        let gap = if second.is_finite() {
            best - second
        } else {
            f64::INFINITY
        };
        let informative = gap > LOG_10_INFORMATIVE_THRESHOLD;
        let usable = java_ranksum_mq_usable(rec);
        let pos = java_read_pos_rank_element(rec, FIRST_REMAINING.1 as i32);
        kv(
            "likelihood_read",
            format!(
                "{qname} FLAG={} MAPQ={} strand={} CIGAR={} best={} informative={informative} mq_usable={usable} pos={pos:?} gap={gap:.4}",
                rec.flags(),
                rec.mapq(),
                strand_label(rec),
                rec.cigar(),
                if best_is_ref { "REF" } else { "ALT" }
            ),
        );
        if !informative || !usable {
            continue;
        }
        let Some(p) = pos else {
            continue;
        };
        if best_is_ref {
            java_ref.push(p);
        } else {
            java_alt.push(p);
        }
    }
    kv(
        "java_equivalent_counts",
        format!("REF={} ALT={}", java_ref.len(), java_alt.len()),
    );
    kv("java_equivalent_ref_pos", format!("{java_ref:?}"));
    kv("java_equivalent_alt_pos", format!("{java_alt:?}"));
    let java_omits = java_ref.is_empty() || java_alt.is_empty();
    kv("java_equivalent_omits", java_omits.to_string());
    assert!(
        java_ref.is_empty(),
        "Java fillQualsFromLikelihood REF list must be empty on this hom-alt site"
    );
    assert!(
        !java_alt.is_empty(),
        "Java-equivalent ALT RankSum evidence must be non-empty"
    );
    assert!(java_omits, "Java empty REF → undefined RankSum → omit");

    let mut rust_ref = Vec::new();
    let mut rust_alt = Vec::new();
    let ref_b = FIRST_REMAINING.2.as_bytes()[0];
    let alt_b = FIRST_REMAINING.3.as_bytes()[0];
    for rec in &covering.reads {
        let Some(base) = read_base_at_ref_coord_1based(rec, FIRST_REMAINING.1 as i32) else {
            continue;
        };
        let offset = (FIRST_REMAINING.1 as i64 - (rec.pos() + 1)) as f64;
        let allele = if base.eq_ignore_ascii_case(&alt_b) {
            rust_alt.push(offset);
            "ALT"
        } else if base.eq_ignore_ascii_case(&ref_b) {
            rust_ref.push(offset);
            "REF"
        } else {
            continue;
        };
        kv(
            "rust_pileup_read",
            format!(
                "{} FLAG={} MAPQ={} strand={} CIGAR={} allele={allele} pileup_offset={offset} in_annotation={}",
                qname(rec),
                rec.flags(),
                rec.mapq(),
                strand_label(rec),
                rec.cigar(),
                annotation_ids.contains(&format!("{} FLAG={}", qname(rec), rec.flags()))
            ),
        );
    }
    kv(
        "rust_pileup_counts",
        format!("REF={} ALT={}", rust_ref.len(), rust_alt.len()),
    );
    let rust_emits = !rust_ref.is_empty() && !rust_alt.is_empty();
    kv("rust_pileup_emits", rust_emits.to_string());
    assert!(
        rust_emits,
        "Rust pileup still fills both REF and ALT lists → Some(z) → INFO insert"
    );
    assert_ne!(
        (rust_ref.is_empty(), rust_alt.is_empty()),
        (java_ref.is_empty(), java_alt.is_empty()),
        "first arrow is RankSum source object / evidence class, not Mann-Whitney math"
    );

    let closed_emitted = try_emit_call_region_variants(
        covering_closed,
        &closed_outcome,
        "NA12878",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .expect("closed emit");
    let emitted =
        try_emit_call_region_variants(covering, &outcome, "NA12878", DEFAULT_STAND_EMIT_CONFIDENCE)
            .expect("emit");
    let closed_rec = closed_emitted
        .iter()
        .find(|r| r.position == CLOSED_AG.1 && r.reference == CLOSED_AG.2)
        .expect("closed emit");
    let target_rec = emitted
        .iter()
        .find(|r| r.position == FIRST_REMAINING.1 && r.reference == FIRST_REMAINING.2)
        .expect("target emit");
    let closed_sample = closed_rec.samples.first().expect("sample");
    let target_sample = target_rec.samples.first().expect("sample");
    assert_eq!(
        closed_sample.gt.as_ref().map(|g| g.to_string()).as_deref(),
        Some("1/1")
    );
    assert_eq!(closed_sample.ad.as_deref(), Some(&[0u32, 2][..]));
    assert_eq!(closed_sample.dp.map(|d| d as i32), Some(2));
    assert!((closed_rec.quality.unwrap_or(0.0) - 78.32).abs() < 0.05);
    assert_eq!(info_i32(&closed_rec.info, "DP"), Some(3));
    let closed_mq = info_f64(&closed_rec.info, "MQ").unwrap_or(-1.0);
    let closed_sor = info_f64(&closed_rec.info, "SOR").unwrap_or(-1.0);
    assert!((closed_mq - 41.96).abs() < 0.005);
    assert!((closed_sor - 0.693).abs() < 0.002);
    assert_eq!(info_i32(&target_rec.info, "DP"), Some(4));
    assert_eq!(target_sample.ad.as_deref(), Some(&[0u32, 3][..]));
    assert!(
        !info_has(&target_rec.info, "ReadPosRankSum"),
        "6R.191 closed the 6R.190 rust-only ReadPosRankSum emit"
    );
    kv(
        "emit_target",
        format!(
            "DP={} MQ={:?} SOR={:?} ReadPosRankSum={:?}",
            info_i32(&target_rec.info, "DP").unwrap_or(-1),
            info_f64(&target_rec.info, "MQ"),
            info_f64(&target_rec.info, "SOR"),
            info_f64(&target_rec.info, "ReadPosRankSum")
        ),
    );
    kv(
        "first_divergent_arrow",
        "6R.190 arrow: Java fillQualsFromLikelihood(annotation AlleleLikelihoods) vs then-Rust region.reads pileup; 6R.191 closed the source",
    );
    kv("classification", "A — WRONG SOURCE OBJECT");
}
