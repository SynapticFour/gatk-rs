//! 6R.187: proof-only fresh INFO reconnaissance after 6R.186.
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! PRODUCTION CHANGE: NONE.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r187_fresh_info_reconnaissance -- --nocapture --test-threads=1
//! HOLDOUT_6R187=1 cargo test -p gatk-haplotypecaller --test holdout_6r187_fresh_info -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::{build_per_haplotype_variation_events, VariationEvent};
use gatk_haplotypecaller::hc_allele_mapping::create_allele_mapper_with_events;
use gatk_haplotypecaller::hc_genotyping_engine::java_alignment_read_overlaps_interval;
use gatk_haplotypecaller::read_realignment::LOG_10_INFORMATIVE_THRESHOLD;
use gatk_haplotypecaller::variant_site_hc_annotations::{
    coverage_evidence_count, rms_mapping_quality_raw, rms_mapping_quality_sample_mapqs,
};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, marginalize_rows_to_biallelic_alleles,
    region_likelihoods_to_rows, traverse_assembly_region_walker, try_emit_call_region_variants,
    AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition, HaplotypeCallerEngine,
    ReadFilterParams, RegionReadLikelihood, WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
};
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
const FIRST_REMAINING: (&str, u64, &str, &str) = ("2", 92_305_635, "A", "G");
const CLOSED_GT: (&str, u64, &str, &str) = ("2", 92_305_634, "G", "T");
const CLOSED_INDEL: (&str, u64, &str, &str) = ("2", 92_307_324, "TTC", "T");
const CLOSED_TG: (&str, u64, &str, &str) = ("2", 92_307_333, "T", "G");
const CLOSED_MID: (&str, u64, &str, &str) = ("2", 92_316_347, "G", "A");
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const GAP_INTERVAL: &str = "2:92305500-92305850";
const TG_INTERVAL: &str = "2:92307200-92307550";
const MARGIN: i32 = 2;
const FLAG_REVERSE: u16 = 0x10;
const JAVA_MQ44_QNAME: &str = "H06JUADXX130110:1:1101:10052:88682";

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
    eprintln!("6R187\t{key}\t{}", value.as_ref());
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

fn dump_read(prefix: &str, rec: &Record, best: &str, informative: bool) {
    eprintln!(
        "6R187\t{prefix}\tQNAME={} FLAG={} MAPQ={} CIGAR={} start={} end={} strand={} best={} informative={}",
        qname(rec),
        rec.flags(),
        rec.mapq(),
        rec.cigar(),
        rec.pos() + 1,
        gatk_haplotypecaller::read_unclip::alignment_end_1based(rec),
        if rec.flags() & FLAG_REVERSE != 0 {
            "rev"
        } else {
            "fwd"
        },
        best,
        informative
    );
}

#[test]
fn forensic_6r187_no_production_change() {
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");
    kv(
        "fresh_pass",
        "6R.187 is a fresh post-6R.186 reconnaissance pass. It does not assume that the remaining annotation differences share the 6R.181/6R.186 cause.",
    );
    let early = include_str!("../src/hc_genotyping_engine/genotype_site_early_template.rs");
    let fin = include_str!("../src/hc_genotyping_engine/genotype_finalize.rs");
    let pipe = include_str!("../src/hc_genotyping_engine/genotype_site_pipeline.rs");
    let ge = include_str!("../src/hc_genotyping_engine/mod.rs");
    let emit = include_str!("../src/region_vcf_emit.rs");
    let ann = include_str!("../src/variant_site_hc_annotations.rs");
    let engine = include_str!("../src/engine.rs");
    let ll = include_str!("../src/engine_likelihoods.rs");
    assert!(
        early.contains("fn annotation_likelihoods_from_stored_haplotypes")
            && early.contains("marginalize_rows_to_biallelic_alleles")
            && early.contains("likelihood_subset_for_event"),
        "6R.186 helper must stay"
    );
    let helper = early
        .split("fn annotation_likelihoods_from_stored_haplotypes")
        .nth(1)
        .expect("helper");
    assert!(
        !helper.contains("92305635") && !helper.contains("92307333"),
        "6R.186 helper has no new 6R.187 locus rule"
    );
    let fin_fn = fin
        .split("fn finish_strict_java_shaped_site_call")
        .nth(1)
        .expect("finish fn");
    let fin_body = fin_fn
        .split("fn finalize_strict_java_variation_genotype")
        .next()
        .expect("body");
    assert!(
        !fin_body.contains("with_annotation_likelihoods"),
        "shared finish helper must not globally attach annotation_likelihoods"
    );
    assert!(
        pipe.contains("with_annotation_likelihoods")
            && ge.contains("fn per_variant_annotation_likelihoods"),
        "6R.180 SiteScore path stays"
    );
    assert!(
        emit.contains("if call.annotation_likelihoods.is_empty()")
            && emit.contains("likelihoods: &call.annotation_likelihoods"),
        "emit still consumes the attached object when present"
    );
    assert!(
        engine.contains("fn apply_java_order_normalize_and_filter")
            && engine.contains("refresh_region_read_likelihoods"),
        "6R.185 stored-matrix lifecycle stays closed"
    );
    assert!(
        ann.contains("6R.174: INFO DP is Java `Coverage.evidenceCount`")
            && ann.contains("mq_rms_of_sample_evidence")
            && ann.contains("strand_bias_sample_counts"),
        "DP/MQ/SOR formulas must stay closed"
    );
    assert!(
        ll.contains("score_pairhmm_from_records_java_mate_contig"),
        "6R.176 mate-contig gate must stay closed"
    );
    assert!(
        !ann.contains("92305635") && !emit.contains("92305635") && !engine.contains("92305635"),
        "no locus-specific 6R.187 production patch"
    );
}

#[test]
fn forensic_6r187_live_vcf_inventory() {
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
        "6R.187 is a fresh post-6R.186 reconnaissance pass. It does not assume that the remaining annotation differences share the 6R.181/6R.186 cause.",
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
    let mut chr2_fmt = 0usize;
    let mut first_fmt: Option<(String, &'static str)> = None;
    let mut first_qual: Option<String> = None;
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
            if first_qual.is_none() {
                first_qual = Some(j.site());
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
    kv("chr2_format_diff", chr2_fmt.to_string());
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
    assert_eq!(
        chr2_fmt, 0,
        "chr2 FORMAT must stay closed after 6R.164–6R.186"
    );
    assert_eq!(java_map.len(), 126);
    assert_eq!(rust_map.len(), 129);

    kv(
        "6r187_snapshot",
        "6R.187 first remaining was 2:92305635 INFO DP vs then-current rust.vcf. 6R.189 closed that site; current 6r43 rust.vcf is post-6R.189. First remaining after that frontier is 6R.190.",
    );
    for (label, spec, dp, mq, sor) in [
        ("closed_92305634", CLOSED_GT, "3", 41.96, 0.693),
        ("closed_92305635", FIRST_REMAINING, "3", 41.96, 0.693),
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
        assert!(
            !j.info.contains_key("InbreedingCoeff") && !r.info.contains_key("InbreedingCoeff"),
            "{label} InbreedingCoeff absent"
        );
        kv(label, format!("{} FORMAT/QUAL/DP/MQ/SOR closed", j.site()));
    }
    let tg = &rust_map[&rec_key(CLOSED_TG)];
    kv(
        "closed_92307333_readpos",
        format!("java_absent rust={:?}", tg.info.get("ReadPosRankSum")),
    );

    let mut first_fmt_before_info: Option<(String, &'static str)> = None;
    let mut first_info: Option<(Rec, Rec, Vec<(String, String, String)>)> = None;
    for k in &common {
        let j = &java_map[k];
        let r = &rust_map[k];
        if let Some(fd) = first_format_field(j, r) {
            first_fmt_before_info = Some((j.site(), fd));
            break;
        }
        if !qual_close(&j.qual, &r.qual) {
            first_fmt_before_info = Some((j.site(), "QUAL"));
            break;
        }
        let val = genuine_info_value_diffs(j, r);
        let rust_only_keys = rust_only_info_keys(j, r);
        let java_only_keys = java_only_info_keys(j, r);
        if !val.is_empty() || !java_only_keys.is_empty() {
            kv(
                "first_info_keys",
                format!("value={val:?} rust_only={rust_only_keys:?} java_only={java_only_keys:?}"),
            );
            first_info = Some((j.clone(), r.clone(), val));
            break;
        }
        // Rust-only ReadPosRankSum is inventoried but is not selected before a
        // numeric INFO split at an earlier locus. Track it only when it is first.
        if rust_only_keys
            .iter()
            .any(|k| k != "ReadPosRankSum" && k != "InbreedingCoeff")
        {
            first_info = Some((j.clone(), r.clone(), val));
            break;
        }
    }
    assert!(
        first_fmt_before_info.is_none(),
        "FORMAT/QUAL must stay closed on chr2 before later-contig FORMAT: {first_fmt_before_info:?}"
    );
    // 6R.187 originally pinned first remaining to 2:92305635 INFO DP=1 vs 3.
    // 6R.189 closed that site in current rust.vcf. This walk still skips
    // rust-only ReadPosRankSum (6R.187 methodology). The fresh first
    // remaining, including rust-only RankSum extras, is 6R.190.
    let (j0, r0, vals) = first_info.expect("expected a remaining INFO split after 6R.189");
    kv("post_6r189_skip_rp_first_info", j0.site());
    kv("post_6r189_skip_rp_vals", format!("{vals:?}"));
    kv(
        "post_6r189_skip_rp_rust_only",
        format!("{:?}", rust_only_info_keys(&j0, &r0)),
    );
    kv(
        "post_6r189_skip_rp_java_only",
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
        "6R.189 closed 2:92305635; this inventory must not still treat it as first remaining"
    );
    kv(
        "primary_arrow_field",
        "6R.187 DP site closed by 6R.189; 6R.190 owns the fresh first remaining",
    );
}

#[test]
fn forensic_6r187_target_annotation_objects() {
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing P12 ref/BAM");
        return;
    }
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");
    kv("target", "2:92305635 A/G");
    kv(
        "fresh_pass",
        "6R.187 is a fresh post-6R.186 reconnaissance pass. It does not assume that the remaining annotation differences share the 6R.181/6R.186 cause.",
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
    let outcome = HaplotypeCallerEngine::call_region(
        covering,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("outcome");

    let stored = unique_likelihood_indices(&outcome.read_likelihoods);
    kv("stored_hap_n", stored.len().to_string());

    let closed = outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(CLOSED_GT.1)
                && c.event.ref_allele == CLOSED_GT.2
                && c.event.alt_allele == CLOSED_GT.3
        })
        .expect("closed G/T");
    let target = outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(FIRST_REMAINING.1)
                && c.event.ref_allele == FIRST_REMAINING.2
                && c.event.alt_allele == FIRST_REMAINING.3
        })
        .expect("target A/G");

    let closed_ann = unique_likelihood_indices(&closed.annotation_likelihoods);
    let target_ann = unique_likelihood_indices(&target.annotation_likelihoods);
    kv(
        "closed_92305634_annotation_n",
        format!(
            "{} empty={}",
            closed_ann.len(),
            closed.annotation_likelihoods.is_empty()
        ),
    );
    kv(
        "target_92305635_annotation_n",
        format!(
            "{} empty={}",
            target_ann.len(),
            target.annotation_likelihoods.is_empty()
        ),
    );
    assert_eq!(closed.genotype.format.ad_as_i32(), vec![0, 2]);
    assert_eq!(closed.genotype.format.dp.as_i32(), 2);
    assert_eq!(target.genotype.format.ad_as_i32(), vec![0, 2]);
    assert_eq!(target.genotype.format.dp.as_i32(), 2);
    assert_eq!(target.genotype.format.pl_as_i32(), vec![90, 6, 0]);
    assert_eq!(target.genotype.format.gq.as_i32(), 6);

    let mut stored_retain_closed = BTreeSet::new();
    let mut stored_retain_target = BTreeSet::new();
    for idx in &stored {
        let Some(rec) = outcome.genotyping_reads.get(*idx) else {
            continue;
        };
        if java_alignment_read_overlaps_interval(rec, CLOSED_GT.1, CLOSED_GT.1, MARGIN) {
            stored_retain_closed.insert(*idx);
        }
        if java_alignment_read_overlaps_interval(rec, FIRST_REMAINING.1, FIRST_REMAINING.1, MARGIN)
        {
            stored_retain_target.insert(*idx);
        }
    }
    kv(
        "stored_retainEvidence_92305634_n",
        stored_retain_closed.len().to_string(),
    );
    kv(
        "stored_retainEvidence_92305635_n",
        stored_retain_target.len().to_string(),
    );
    assert_eq!(stored_retain_closed.len(), 3);
    assert_eq!(
        stored_retain_target.len(),
        3,
        "Java loc-loop retainEvidence on the stored hap matrix is n=3 at both adjacent SNPs"
    );

    for idx in &stored_retain_target {
        let rec = &outcome.genotyping_reads[*idx];
        dump_read("JAVA_RETAIN", rec, "?", false);
    }
    for idx in &target_ann {
        let rec = &outcome.genotyping_reads[*idx];
        dump_read("RUST_ATTACHED", rec, "?", false);
    }
    for idx in stored_retain_target.difference(&target_ann) {
        let rec = &outcome.genotyping_reads[*idx];
        dump_read("MISSING_FROM_RUST_OBJECT", rec, "?", false);
    }

    let (full_ref, full_pad) = outcome.assembly.event_map_reference();
    let apply_pad = outcome
        .assembly
        .haplotypes
        .iter()
        .find(|h| h.is_reference)
        .and_then(|h| h.genome_loc.as_ref().map(|g| g.start.get()))
        .unwrap_or(full_pad);
    let hap_cache = build_per_haplotype_variation_events(
        &outcome.assembly.haplotypes,
        full_ref,
        full_pad,
        outcome.assembly.max_mnp_distance(),
        &covering.contig,
    );
    let mapping = create_allele_mapper_with_events(
        &VariationEvent::from_alleles(
            &covering.contig,
            FIRST_REMAINING.1,
            FIRST_REMAINING.2,
            FIRST_REMAINING.3,
        ),
        FIRST_REMAINING.1,
        &outcome.assembly.haplotypes,
        apply_pad,
        outcome.assembly.reference_bases(),
        outcome.assembly.max_mnp_distance(),
        true,
        Some(&hap_cache),
    );
    let rows =
        region_likelihoods_to_rows(&outcome.read_likelihoods, outcome.assembly.haplotypes.len());
    let marg = marginalize_rows_to_biallelic_alleles(
        &rows,
        &mapping.ref_haplotype_indices,
        &mapping.alt_haplotype_indices,
    );
    for idx in stored_retain_target.union(&target_ann) {
        let rec = &outcome.genotyping_reads[*idx];
        let Some(row) = marg.iter().find(|r| r.read_index == *idx) else {
            dump_read("NO_MARG_ROW", rec, "none", false);
            continue;
        };
        let ll_ref = row
            .haplotype_log10_likelihoods
            .first()
            .copied()
            .unwrap_or(f64::NEG_INFINITY);
        let ll_alt = row
            .haplotype_log10_likelihoods
            .get(1)
            .copied()
            .unwrap_or(f64::NEG_INFINITY);
        let (best, gap) = if ll_alt > ll_ref {
            ("G", ll_alt - ll_ref)
        } else {
            ("A", ll_ref - ll_alt)
        };
        let informative = gap > LOG_10_INFORMATIVE_THRESHOLD;
        dump_read("MARG", rec, best, informative);
        kv(
            "marg_ll",
            format!(
                "{} A*={ll_ref:.6} G={ll_alt:.6} gap={gap:.6} informative={informative}",
                qname(rec)
            ),
        );
    }

    let closed_cov = if closed.annotation_likelihoods.is_empty() {
        coverage_evidence_count(
            &outcome.genotyping_reads,
            &outcome.read_likelihoods,
            CLOSED_GT.1,
            CLOSED_GT.1,
            MARGIN,
        )
    } else {
        coverage_evidence_count(
            &outcome.genotyping_reads,
            &closed.annotation_likelihoods,
            CLOSED_GT.1,
            CLOSED_GT.1,
            MARGIN,
        )
    };
    let target_cov = if target.annotation_likelihoods.is_empty() {
        coverage_evidence_count(
            &outcome.genotyping_reads,
            &outcome.read_likelihoods,
            FIRST_REMAINING.1,
            FIRST_REMAINING.1,
            MARGIN,
        )
    } else {
        coverage_evidence_count(
            &outcome.genotyping_reads,
            &target.annotation_likelihoods,
            FIRST_REMAINING.1,
            FIRST_REMAINING.1,
            MARGIN,
        )
    };
    kv("closed_coverage_evidence_count", closed_cov.to_string());
    kv("target_coverage_evidence_count", target_cov.to_string());
    assert_eq!(closed_cov, 3);
    assert_eq!(target_cov, 3, "6R.189 loc-loop Coverage.evidenceCount n=3");

    let target_mqs = if target.annotation_likelihoods.is_empty() {
        rms_mapping_quality_sample_mapqs(&outcome.genotyping_reads, &outcome.read_likelihoods)
    } else {
        rms_mapping_quality_sample_mapqs(&outcome.genotyping_reads, &target.annotation_likelihoods)
    };
    kv("target_mq_members", format!("{target_mqs:?}"));
    let target_rms = rms_mapping_quality_raw(&target_mqs)
        .map(|t| t.2)
        .unwrap_or(-1.0);
    kv("target_mq_rms", format!("{target_rms:.4}"));
    assert!(
        (target_rms - 41.96).abs() < 0.005,
        "6R.189 loc-loop MAPQ members produce Java MQ=41.96"
    );

    assert!(
        !target.annotation_likelihoods.is_empty(),
        "92305635 has an attached annotation object (not the empty 6R.181 fallback)"
    );
    assert_eq!(target_ann.len(), 3, "6R.189 attached loc-loop n=3");
    assert!(
        closed.annotation_likelihoods.is_empty() || closed_ann.len() == 3,
        "closed neighbor is empty-fallback or n=3"
    );
    assert_eq!(
        target_ann, stored_retain_target,
        "6R.189: attached object is Java stored-hap retainEvidence n=3"
    );

    let emitted =
        try_emit_call_region_variants(covering, &outcome, "NA12878", DEFAULT_STAND_EMIT_CONFIDENCE)
            .expect("emit");
    let closed_rec = emitted
        .iter()
        .find(|r| r.position == CLOSED_GT.1 && r.reference == CLOSED_GT.2)
        .expect("closed emit");
    let target_rec = emitted
        .iter()
        .find(|r| r.position == FIRST_REMAINING.1 && r.reference == FIRST_REMAINING.2)
        .expect("target emit");
    assert_eq!(info_i32(&closed_rec.info, "DP"), Some(3));
    assert_eq!(info_i32(&target_rec.info, "DP"), Some(3));
    let closed_mq = info_f64(&closed_rec.info, "MQ").unwrap_or(-1.0);
    let target_mq = info_f64(&target_rec.info, "MQ").unwrap_or(-1.0);
    let closed_sor = info_f64(&closed_rec.info, "SOR").unwrap_or(-1.0);
    let target_sor = info_f64(&target_rec.info, "SOR").unwrap_or(-1.0);
    kv(
        "emitted",
        format!(
            "closed DP=3 MQ={closed_mq:.2} SOR={closed_sor:.3}; target DP=3 MQ={target_mq:.2} SOR={target_sor:.3}"
        ),
    );
    assert!((closed_mq - 41.96).abs() < 0.005);
    assert!((target_mq - 41.96).abs() < 0.005);
    assert!((closed_sor - 0.693).abs() < 0.002);
    assert!((target_sor - 0.693).abs() < 0.002);
    assert!(!info_has(&target_rec.info, "InbreedingCoeff"));

    kv(
        "classification",
        "A closed by 6R.189 — SiteScore annotation is stored-hap retainEvidence n=3",
    );
    kv(
        "not_6r186",
        "cluster-TG 6R.186 was missing construction on an empty object; this site has a non-empty object with the wrong membership",
    );

    let tg_specs = parse_intervals_cli_string(&dict, TG_INTERVAL).expect("tg");
    let tg_walk = traverse_assembly_region_walker(
        &dict,
        &tg_specs,
        &ref_fasta,
        &bam,
        &ReadFilterParams::gatk_standard_hc(),
        &WalkerTraversalConfig::gatk_haplotype_caller_production(100),
    )
    .expect("tg walk");
    let tg_regions = flatten_assembly_regions(&tg_walk);
    let tg_covering = tg_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_TG.1
                && r.end.get() >= CLOSED_TG.1
        })
        .expect("tg covering");
    let tg_outcome = HaplotypeCallerEngine::call_region(
        tg_covering,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("tg call")
    .expect("tg outcome");
    let tg_call = tg_outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(CLOSED_TG.1)
                && c.event.ref_allele == CLOSED_TG.2
                && c.event.alt_allele == CLOSED_TG.3
        })
        .expect("T/G");
    let tg_ann = unique_likelihood_indices(&tg_call.annotation_likelihoods);
    kv("closed_92307333_annotation_n", tg_ann.len().to_string());
    assert_eq!(tg_ann.len(), 1, "6R.186 annotation_likelihoods n=1 stays");
    let keep_idx = *tg_ann.iter().next().expect("n=1");
    let rec = tg_outcome.genotyping_reads.get(keep_idx).expect("read");
    assert_eq!(qname(rec), JAVA_MQ44_QNAME);
    assert_eq!(rec.mapq(), 44);
    let tg_emitted = try_emit_call_region_variants(
        tg_covering,
        &tg_outcome,
        "NA12878",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .expect("tg emit");
    let tg_rec = tg_emitted
        .iter()
        .find(|r| r.position == CLOSED_TG.1 && r.reference == CLOSED_TG.2)
        .expect("tg record");
    assert_eq!(info_i32(&tg_rec.info, "DP"), Some(1));
    let tg_mq = info_f64(&tg_rec.info, "MQ").unwrap_or(-1.0);
    let tg_sor = info_f64(&tg_rec.info, "SOR").unwrap_or(-1.0);
    assert!((tg_mq - 44.0).abs() < 0.005);
    assert!((tg_sor - 1.609).abs() < 0.002);
    kv(
        "closed_92307333_info",
        format!(
            "DP=1 MQ={tg_mq} SOR={tg_sor} ReadPosRankSum={}",
            info_has(&tg_rec.info, "ReadPosRankSum")
        ),
    );
}
