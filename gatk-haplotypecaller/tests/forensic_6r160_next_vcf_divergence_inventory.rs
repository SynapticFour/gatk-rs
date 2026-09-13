//! 6R.160 reconnaissance: inventory the next genuine Java-vs-Rust VCF divergence.
//!
//! Does **not** pin a future production FORMAT contract. It records the current
//! emitted-callset difference against frozen Java 4.4.0.0 VCFs.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r160_next_vcf_divergence_inventory -- --nocapture --test-threads=1
//! ```

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";

/// 6R.43 regions whose Rust VCF was regenerated this campaign (BAM present).
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

/// BAMs missing; rust.vcf is stale and must not enter first-divergence selection.
const STALE_MISSING_BAM: &[&str] = &["chr20_w47", "chr21_w10"];

const FIRST_CONTIG: &str = "2";
const FIRST_POS: u64 = 92_316_347;
const FIRST_REF: &str = "G";
const FIRST_ALT: &str = "A";
const JAVA_AD: &str = "0,3";
const JAVA_DP: &str = "3";
const JAVA_GQ: &str = "9";
const JAVA_PL: &str = "135,9,0";
const JAVA_GT: &str = "1/1";
const RUST_AD: &str = "0,1";
const RUST_DP: &str = "1";
const RUST_GQ: &str = "3";
const RUST_PL: &str = "45,3,0";
const RUST_GT: &str = "1/1";

const CANONICAL_POS: u64 = 29_456_196;

#[derive(Clone, Debug)]
struct Rec {
    chrom: String,
    pos: u64,
    ref_a: String,
    alt: String,
    qual: String,
    filter: String,
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
            gt: sample.get("GT").cloned().unwrap_or_default(),
            ad: sample.get("AD").cloned().unwrap_or_default(),
            dp: sample.get("DP").cloned().unwrap_or_default(),
            gq: sample.get("GQ").cloned().unwrap_or_default(),
            pl: sample.get("PL").cloned().unwrap_or_default(),
        });
    }
    out
}

const QUAL_TOL: f64 = 0.05;

fn format_fields_equal(j: &Rec, r: &Rec) -> bool {
    j.gt == r.gt
        && j.ad == r.ad
        && j.dp == r.dp
        && j.gq == r.gq
        && j.pl == r.pl
        && j.filter == r.filter
}

fn qual_close(a: &str, b: &str) -> bool {
    match (a.parse::<f64>(), b.parse::<f64>()) {
        (Ok(x), Ok(y)) => (x - y).abs() <= QUAL_TOL,
        _ => a == b,
    }
}

/// Reverse-trim / allele-string split. Identical CHROM/POS/REF/ALT is not representation-only.
fn is_representation_only(j: &Rec, r: &Rec) -> bool {
    j.chrom != r.chrom || j.pos != r.pos || j.ref_a != r.ref_a || j.alt != r.alt
}

fn load_region(root: &Path, id: &str) -> (Vec<Rec>, Vec<Rec>) {
    let dir = root.join("parity/reports/6r43").join(id);
    (
        parse_vcf(&dir.join("java.vcf")),
        parse_vcf(&dir.join("rust.vcf")),
    )
}

#[test]
fn forensic_6r160_next_vcf_divergence_inventory() {
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
    for id in CURRENT_REGIONS {
        let (jv, rv) = load_region(&root, id);
        eprintln!("6R160\tregion\t{id}\tjava={}\trust={}", jv.len(), rv.len());
        for rec in jv {
            java_map.insert(rec.key(), rec);
        }
        for rec in rv {
            rust_map.insert(rec.key(), rec);
        }
    }
    for id in STALE_MISSING_BAM {
        let (jv, rv) = load_region(&root, id);
        eprintln!(
            "6R160\texcluded_stale\t{id}\tjava={}\trust={}\treason=BAM_missing_rust_vcf_not_regenerated",
            jv.len(),
            rv.len()
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
    let mut common_identical = 0usize;
    let mut common_different = Vec::new();
    for k in &common {
        let j = &java_map[k];
        let r = &rust_map[k];
        if !format_fields_equal(j, r) {
            common_different.push(k.clone());
        } else if qual_close(&j.qual, &r.qual) {
            common_identical += 1;
        } else {
            eprintln!(
                "6R160\tqual_beyond_tol\t{}\tjava={}\trust={}",
                j.site(),
                j.qual,
                r.qual
            );
            common_different.push(k.clone());
        }
    }

    eprintln!("6R160\tjava_pin\t{JAVA_PIN}");
    eprintln!("6R160\tjava_records\t{}", java_map.len());
    eprintln!("6R160\trust_records\t{}", rust_map.len());
    eprintln!("6R160\tjava_only\t{}", java_only.len());
    eprintln!("6R160\trust_only\t{}", rust_only.len());
    eprintln!("6R160\tcommon\t{}", common.len());
    eprintln!("6R160\tcommon_identical\t{common_identical}");
    eprintln!("6R160\tcommon_different\t{}", common_different.len());
    if let Some(k) = java_only.first() {
        eprintln!("6R160\tfirst_java_only\t{}", java_map[k].site());
    }
    if let Some(k) = rust_only.first() {
        eprintln!("6R160\tfirst_rust_only\t{}", rust_map[k].site());
    }
    if let Some(k) = common_different.first() {
        let j = &java_map[k];
        let r = &rust_map[k];
        eprintln!(
            "6R160\tfirst_common_different\t{}\tjava_GT={} AD={} DP={} GQ={} PL={} QUAL={}\trust_GT={} AD={} DP={} GQ={} PL={} QUAL={}",
            j.site(),
            j.gt, j.ad, j.dp, j.gq, j.pl, j.qual,
            r.gt, r.ad, r.dp, r.gq, r.pl, r.qual
        );
        eprintln!(
            "6R160\trepresentation_only\t{}",
            is_representation_only(j, r)
        );
        eprintln!("6R160\tclassification\tTRUE_CALLSET_DIFFERENCE");
    }

    let first = common_different
        .first()
        .expect("need a common-different site");
    assert_eq!(first.0, FIRST_CONTIG);
    assert_eq!(first.1, FIRST_POS);
    assert_eq!(first.2, FIRST_REF);
    assert_eq!(first.3, FIRST_ALT);
    let j = &java_map[first];
    let r = &rust_map[first];
    assert!(
        !is_representation_only(j, r),
        "same alleles; not reverse-trim"
    );
    assert_eq!(j.gt, JAVA_GT);
    assert_eq!(j.ad, JAVA_AD);
    assert_eq!(j.dp, JAVA_DP);
    assert_eq!(j.gq, JAVA_GQ);
    assert_eq!(j.pl, JAVA_PL);
    assert_eq!(r.gt, RUST_GT);
    assert_eq!(r.ad, RUST_AD);
    assert_eq!(r.dp, RUST_DP);
    assert_eq!(r.gq, RUST_GQ);
    assert_eq!(r.pl, RUST_PL);

    let first_exclusive = java_only
        .iter()
        .map(|k| (k.0.as_str(), k.1))
        .chain(rust_only.iter().map(|k| (k.0.as_str(), k.1)))
        .min();
    if let Some((c, p)) = first_exclusive {
        assert!(
            (c, p) > (FIRST_CONTIG, FIRST_POS),
            "exclusive record {c}:{p} is earlier than the selected common-different site"
        );
    }

    let java_canon = java_map
        .keys()
        .any(|k| k.0 == "20" && k.1 == CANONICAL_POS && k.2 == "A" && k.3 == "T");
    let rust_canon = rust_map
        .keys()
        .any(|k| k.0 == "20" && k.1 == CANONICAL_POS && k.2 == "A" && k.3 == "T");
    eprintln!(
        "6R160\tcanonical_29456196_java_vcf\t{}",
        if java_canon { "present" } else { "absent" }
    );
    eprintln!(
        "6R160\tcanonical_29456196_production_rust_vcf\t{}",
        if rust_canon { "present" } else { "absent" }
    );
    assert!(
        !java_canon,
        "frozen Java 6R.43 chr20_tiny must omit 20:29456196"
    );
    // Production run_haplotype_caller still emits this locus (legacy_1024 / walker
    // regions). 6R.159 closed the diagnostic k-best + Java-bounds path, not this VCF.
    eprintln!(
        "6R160\tcanonical_note\t6R.159 diagnostic path both-absent; production rust.vcf still emits 0/1; later than 2:92316347"
    );

    eprintln!("6R160\taf_10_vs_30_causal_for_target\tNO\tboth_emit_hom_alt");
    eprintln!("6R160\tproduction_change\tNONE");
}
