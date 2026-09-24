//! 6R.242: proof-only fresh VCF inventory after 6R.241.
//! 6R.241 closed covering `20:29455314 G>C` emit. Do not re-investigate it.
//! PRODUCTION CHANGE: NONE.
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r242_fresh_parity_inventory -- --nocapture --test-threads=1
//! HOLDOUT_6R242=1 cargo test -p gatk-haplotypecaller --test holdout_6r242_fresh_parity_inventory -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::emit_gates::java_emit_would_pass;
use gatk_haplotypecaller::genotyping::best_pl_index;
use gatk_haplotypecaller::hc_genotyping_engine::{
    java_emit_af_decision, HcGenotypingConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use gatk_haplotypecaller::read_threading_assembler::DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH;
use gatk_haplotypecaller::variant_site_hc_annotations::coverage_evidence_count;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig, DEFAULT_STAND_CALL_CONF,
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
const CLOSED_GC: (&str, u64, &str, &str) = ("20", 29_455_314, "G", "C");
const FIRST_COMMON: (&str, u64, &str, &str) = ("20", 29_455_649, "T", "TGTTTG");
const FIRST_RUST_ONLY: (&str, u64, &str, &str) = ("20", 29_456_196, "A", "T");
const FIRST_JAVA_ONLY: (&str, u64, &str, &str) = ("20", 29_455_902, "G", "A");
const CLOSED_GA: (&str, u64, &str, &str) = ("20", 29_455_379, "G", "A");
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
const CHR20_INTERVAL: &str = "20:29455000-29456500";
const CHR20_BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const JAVA_VCF_REL: &str = "parity/reports/6r43/chr20_tiny/java.vcf";

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
    eprintln!("6R242\t{key}\t{}", value.as_ref());
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

#[test]
fn forensic_6r242_dump_inventory() {
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");
    kv("k", DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH.to_string());
    assert_eq!(DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH, 128);
    assert_eq!(
        HcGenotypingConfig::strict_java().stand_emit_confidence,
        DEFAULT_STAND_CALL_CONF
    );

    let root = repo_root();
    let reports = root.join("parity/reports/6r43");
    for id in STALE_MISSING_BAM {
        kv(
            "excluded_stale_bam",
            format!("{id} rust.vcf not used (BAM unavailable)"),
        );
    }
    let mut java_map = BTreeMap::new();
    let mut rust_map = BTreeMap::new();
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
        assert!(
            reports.join(id).join("java.vcf").is_file(),
            "missing java.vcf for {id}"
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
        }
        let val = genuine_info_value_diffs(j, r);
        let rok = rust_only_info_keys(j, r);
        let jok = java_only_info_keys(j, r);
        if !val.is_empty() || !rok.is_empty() || !jok.is_empty() {
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
    for k in &java_only {
        kv("java_only_site", java_map[k].site());
    }
    for k in &rust_only {
        let r = &rust_map[k];
        kv(
            "rust_only_site",
            format!(
                "{}\tGT={}\tAD={}\tDP={}\tGQ={}\tPL={}\tQUAL={}\tINFO_DP={}",
                r.site(),
                r.gt,
                r.ad,
                r.dp,
                r.gq,
                r.pl,
                r.qual,
                r.info.get("DP").cloned().unwrap_or_default()
            ),
        );
    }
    if let Some((site, fd)) = &first_fmt {
        kv("first_format_diff", format!("{site} field={fd}"));
    }

    assert!(
        !java_map.contains_key(&rec_key(CLOSED_GC)),
        "Java still omits 20:29455314 G>C"
    );
    assert!(
        !rust_map.contains_key(&rec_key(CLOSED_GC)),
        "6R.241: Rust must omit 20:29455314 G>C"
    );
    assert_eq!(chr2_fmt, 0);
    assert_eq!(chr2_qual, 0);
    assert_eq!(java_map.len(), 126);
    assert_eq!(rust_map.len(), 126);
    assert_eq!(java_only.len(), 2);
    assert_eq!(rust_only.len(), 2);
    assert_eq!(common.len(), 124);
    assert_eq!(java_map[&java_only[0]].site(), "20:29455902 G/A");
    assert_eq!(rust_map[&rust_only[0]].site(), "20:29456196 A/T");
    assert_eq!(
        (
            java_only[0].0.as_str(),
            java_only[0].1,
            java_only[0].2.as_str(),
            java_only[0].3.as_str()
        ),
        FIRST_JAVA_ONLY
    );
    assert_eq!(
        (
            rust_only[0].0.as_str(),
            rust_only[0].1,
            rust_only[0].2.as_str(),
            rust_only[0].3.as_str()
        ),
        FIRST_RUST_ONLY
    );
    let ro0 = &rust_map[&rust_only[0]];
    assert_eq!(ro0.gt, "0/1");
    assert_eq!(ro0.ad, "37,8");
    assert_eq!(ro0.dp, "45");
    assert_eq!(ro0.gq, "99");
    assert_eq!(ro0.pl, "101,0,1291");
    assert_eq!(ro0.qual, "93.64");
    assert!(extra_info_keys.is_empty());

    let ga = rec_key(CLOSED_GA);
    assert!(format_equal(&java_map[&ga], &rust_map[&ga]));

    let (j0, r0, vals) = first_info.expect("remaining INFO split");
    kv("first_genuine_common", j0.site());
    kv("first_info_vals", format!("{vals:?}"));
    kv(
        "first_common_java",
        format!(
            "GT={} AD={} DP={} GQ={} PL={} QUAL={} INFO_DP={}",
            j0.gt,
            j0.ad,
            j0.dp,
            j0.gq,
            j0.pl,
            j0.qual,
            j0.info.get("DP").cloned().unwrap_or_default()
        ),
    );
    kv(
        "first_common_rust",
        format!(
            "GT={} AD={} DP={} GQ={} PL={} QUAL={} INFO_DP={}",
            r0.gt,
            r0.ad,
            r0.dp,
            r0.gq,
            r0.pl,
            r0.qual,
            r0.info.get("DP").cloned().unwrap_or_default()
        ),
    );
    assert_eq!(
        (
            j0.chrom.as_str(),
            j0.pos,
            j0.ref_a.as_str(),
            j0.alt.as_str()
        ),
        FIRST_COMMON
    );
    assert_eq!(j0.gt, r0.gt);
    assert_eq!(j0.ad, r0.ad);
    assert_eq!(j0.dp, r0.dp);
    assert_eq!(j0.gq, r0.gq);
    assert!(qual_close(&j0.qual, &r0.qual));
    assert_eq!(j0.pl, "570,0,3517");
    assert_eq!(r0.pl, "570,0,3518");
    assert!(
        vals.iter()
            .any(|(k, jv, rv)| k == "DP" && jv == "123" && rv == "230"),
        "first genuine remaining field at common site is INFO DP 123 vs 230, got {vals:?}"
    );
    assert!(
        first_fmt
            .as_ref()
            .is_some_and(|(site, fd)| site.starts_with("20:29455649") && *fd == "PL"),
        "first FORMAT numeric split is PL ±1 at the same site: {first_fmt:?}"
    );
}

#[test]
fn forensic_6r242_first_common_live_layers() {
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");
    let root = repo_root();
    let java_vcf = root.join(JAVA_VCF_REL);
    let java_recs = parse_vcf(&java_vcf);
    let neigh: Vec<_> = java_recs
        .iter()
        .filter(|r| r.pos >= 29_455_621 && r.pos <= 29_455_674)
        .collect();
    kv(
        "java_neighborhood_29455621_29455674",
        neigh
            .iter()
            .map(|r| {
                format!(
                    "{} QUAL={} GT={} AD={} PL={} INFO_DP={}",
                    r.site(),
                    r.qual,
                    r.gt,
                    r.ad,
                    r.pl,
                    r.info.get("DP").cloned().unwrap_or_default()
                )
            })
            .collect::<Vec<_>>()
            .join(" | "),
    );
    let j_target = java_recs
        .iter()
        .find(|r| r.pos == FIRST_COMMON.1 && r.ref_a == FIRST_COMMON.2 && r.alt == FIRST_COMMON.3)
        .expect("Java VCF has T/TGTTTG");
    kv(
        "java_target",
        format!(
            "{} QUAL={} GT={} AD={} FORMAT_DP={} GQ={} PL={} INFO_DP={}",
            j_target.site(),
            j_target.qual,
            j_target.gt,
            j_target.ad,
            j_target.dp,
            j_target.gq,
            j_target.pl,
            j_target.info.get("DP").cloned().unwrap_or_default()
        ),
    );

    let rust_only_neigh: Vec<_> = java_recs
        .iter()
        .filter(|r| r.pos >= 29_456_168 && r.pos <= 29_456_261)
        .collect();
    kv(
        "java_neighborhood_around_first_rust_only",
        rust_only_neigh
            .iter()
            .map(|r| {
                format!(
                    "{} QUAL={} GT={} AD={} PL={}",
                    r.site(),
                    r.qual,
                    r.gt,
                    r.ad,
                    r.pl
                )
            })
            .collect::<Vec<_>>()
            .join(" | "),
    );
    assert!(
        !java_recs.iter().any(|r| r.pos == FIRST_RUST_ONLY.1
            && r.ref_a == FIRST_RUST_ONLY.2
            && r.alt == FIRST_RUST_ONLY.3),
        "Java covering VCF omits first rust-only 20:29456196 A/T"
    );

    let ref_fasta = root.join(REF_REL);
    let bam = root.join(CHR20_BAM_REL);
    assert!(ref_fasta.is_file() && bam.is_file());
    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, CHR20_INTERVAL).expect("interval");
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
            ) && r.start.get() <= FIRST_COMMON.1
                && r.end.get() >= FIRST_COMMON.1
        })
        .expect("ActiveFull containing 29455649");
    kv(
        "live_active_region",
        format!(
            "{}:{}-{} extended={}-{} reads={}",
            region.contig,
            region.start.get(),
            region.end.get(),
            region.extended_start.get(),
            region.extended_end.get(),
            region.reads.len()
        ),
    );

    let outcome = HaplotypeCallerEngine::call_region(
        region,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");
    kv("hap_n", outcome.assembly.haplotypes.len().to_string());
    let events_at: Vec<_> = outcome
        .assembly
        .variation_events()
        .iter()
        .filter(|e| e.start_1based.get() == FIRST_COMMON.1)
        .collect();
    kv(
        "eventmap_at_loc",
        format!(
            "n={}\t{}",
            events_at.len(),
            events_at
                .iter()
                .map(|e| format!("{}/{}", e.ref_allele, e.alt_allele))
                .collect::<Vec<_>>()
                .join(",")
        ),
    );
    assert!(
        events_at
            .iter()
            .any(|e| e.ref_allele == FIRST_COMMON.2 && e.alt_allele == FIRST_COMMON.3),
        "EventMap contains T/TGTTTG"
    );
    let alt_haps = outcome
        .assembly
        .haplotypes
        .iter()
        .filter(|h| !h.is_reference)
        .count();
    kv("alt_hap_n", alt_haps.to_string());

    let call = outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based.get() == FIRST_COMMON.1
                && c.event.ref_allele == FIRST_COMMON.2
                && c.event.alt_allele == FIRST_COMMON.3
        })
        .expect("genotyped T/TGTTTG");
    let pl = call.genotype.format.pl_as_i32();
    kv(
        "live_genotype",
        format!(
            "GT_idx={} AD={:?} PL={:?} GQ={} extra_alt={} post_merge={} ann_n={}",
            best_pl_index(&call.genotype.format.pl),
            call.genotype.format.ad_as_i32(),
            pl,
            call.genotype.format.gq.as_i32(),
            call.extra_alt_alleles.len(),
            call.post_merge_unused_alt_subset,
            call.annotation_likelihoods.len()
        ),
    );
    assert_eq!(best_pl_index(&call.genotype.format.pl), 1);
    assert_eq!(call.genotype.format.ad_as_i32(), vec![88, 22]);
    assert_eq!(pl[0], 570);
    assert_eq!(pl[1], 0);
    assert!(
        (pl[2] - 3517).abs() <= 1,
        "PL hom-alt is Java 3517 ±1 representation, got {}",
        pl[2]
    );
    let gls = &call.genotype.genotype_log10_likelihoods;
    let af = java_emit_af_decision(gls, DEFAULT_STAND_CALL_CONF).expect("af");
    kv(
        "emit_pred",
        format!(
            "stand30_pass={} phred={:.2} mono={} not_empty_ann={}",
            af.passes_emit,
            af.phred_scaled,
            af.site_is_monomorphic,
            !call.annotation_likelihoods.is_empty()
        ),
    );
    assert!(af.passes_emit, "site already passes Java calling threshold");
    assert!(java_emit_would_pass(
        &call.event,
        gls,
        &call.genotype.format,
        DEFAULT_STAND_CALL_CONF,
        &[],
    )
    .unwrap());
    let cov = coverage_evidence_count(
        &outcome.genotyping_reads,
        &call.annotation_likelihoods,
        FIRST_COMMON.1,
        FIRST_COMMON.1,
        2,
    );
    kv(
        "annotation_coverage_evidence_count",
        format!(
            "ann_n={} coverage={}",
            call.annotation_likelihoods.len(),
            cov
        ),
    );

    let recs =
        try_emit_call_region_variants(region, &outcome, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .unwrap_or_default();
    let emitted = recs.iter().find(|r| {
        r.position == FIRST_COMMON.1
            && r.reference == FIRST_COMMON.2
            && r.alternate.iter().any(|a| a == FIRST_COMMON.3)
    });
    assert!(emitted.is_some(), "both engines emit T/TGTTTG");
    let er = emitted.unwrap();
    kv(
        "live_emit_info_dp",
        info_i32(&er.info, "DP")
            .map(|d| d.to_string())
            .unwrap_or_default(),
    );
    kv(
        "layers",
        "EventMap has T/TGTTTG; GT/AD match Java; both emit; PL ±1 representation; INFO DP 123 vs 230 is the first genuine remaining field",
    );
    kv(
        "classification",
        "common-site INFO (Coverage/annotation evidence) after matching EventMap/allele-set/GT/AD/emission; not EMISSION_PREDICATE",
    );
}
