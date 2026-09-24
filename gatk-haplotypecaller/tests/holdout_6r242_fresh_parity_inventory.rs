//! 6R.242 live: after 6R.241, covering G>C is omitted. First remaining
//! genuine VCF split is common-site INFO at `20:29455649 T/TGTTTG`.
//! PRODUCTION CHANGE: NONE. Skipped unless `HOLDOUT_6R242=1`.
//!
//! ```text
//! HOLDOUT_6R242=1 cargo test -p gatk-haplotypecaller --test holdout_6r242_fresh_parity_inventory -- --nocapture --test-threads=1
//! ```

use gatk_haplotypecaller::hc_genotyping_engine::HcGenotypingConfig;
use gatk_haplotypecaller::read_threading_assembler::DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH;
use gatk_haplotypecaller::DEFAULT_STAND_CALL_CONF;
use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

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
const QUAL_TOL: f64 = 0.05;
const INFO_PRINT: &[&str] = &[
    "FS",
    "SOR",
    "ExcessHet",
    "ReadPosRankSum",
    "BaseQRankSum",
    "MQRankSum",
    "MQ",
    "QD",
    "AF",
    "MLEAF",
];

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R242\t{key}\t{}", value.as_ref());
}

fn parse_info(s: &str) -> BTreeMap<String, String> {
    let mut out = BTreeMap::new();
    if s.is_empty() || s == "." {
        return out;
    }
    for part in s.split(';') {
        if let Some((k, v)) = part.split_once('=') {
            out.insert(k.to_string(), v.to_string());
        } else if !part.is_empty() {
            out.insert(part.to_string(), "true".to_string());
        }
    }
    out
}

fn keys(
    path: &Path,
) -> BTreeMap<
    (String, u64, String, String),
    (
        String,
        BTreeMap<String, String>,
        String,
        String,
        String,
        String,
        String,
        String,
    ),
> {
    let mut out = BTreeMap::new();
    for line in fs::read_to_string(path).unwrap_or_default().lines() {
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
        let alt0 = f[4].split(',').next().unwrap_or(f[4]).to_string();
        out.insert(
            (
                f[0].to_string(),
                f[1].parse().unwrap(),
                f[3].to_string(),
                alt0,
            ),
            (
                f[5].to_string(),
                parse_info(f[7]),
                sample.get("GT").cloned().unwrap_or_default(),
                sample.get("AD").cloned().unwrap_or_default(),
                sample.get("DP").cloned().unwrap_or_default(),
                sample.get("GQ").cloned().unwrap_or_default(),
                sample.get("PL").cloned().unwrap_or_default(),
                f[6].to_string(),
            ),
        );
    }
    out
}

fn info_print_close(key: &str, jv: &str, rv: &str) -> bool {
    if jv == rv {
        return true;
    }
    let (Ok(j), Ok(r)) = (jv.parse::<f64>(), rv.parse::<f64>()) else {
        return false;
    };
    if matches!(
        key,
        "FS" | "SOR" | "ExcessHet" | "ReadPosRankSum" | "BaseQRankSum" | "MQRankSum"
    ) {
        (j * 1000.0).round() == (r * 1000.0).round()
    } else if matches!(key, "MQ" | "QD" | "AF" | "MLEAF") {
        (j * 100.0).round() == (r * 100.0).round() || (j - r).abs() < 0.005
    } else {
        false
    }
}

fn genuine_info(j: &BTreeMap<String, String>, r: &BTreeMap<String, String>) -> bool {
    for (k, jv) in j {
        if let Some(rv) = r.get(k) {
            if jv != rv && !info_print_close(k, jv, rv) {
                return true;
            }
        } else {
            return true;
        }
    }
    r.keys().any(|k| !j.contains_key(k))
}

#[test]
fn holdout_6r242_fresh_parity_inventory() {
    if std::env::var("HOLDOUT_6R242").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R242=1");
        return;
    }
    assert_eq!(DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH, 128);
    assert_eq!(
        HcGenotypingConfig::strict_java().stand_emit_confidence,
        DEFAULT_STAND_CALL_CONF
    );
    let root = repo_root();
    let mut java = BTreeMap::new();
    let mut rust = BTreeMap::new();
    for id in CURRENT_REGIONS {
        let dir = root.join("parity/reports/6r43").join(id);
        java.extend(keys(&dir.join("java.vcf")));
        rust.extend(keys(&dir.join("rust.vcf")));
    }
    let jk: Vec<_> = java.keys().cloned().collect();
    let java_only: Vec<_> = jk
        .iter()
        .filter(|k| !rust.contains_key(*k))
        .cloned()
        .collect();
    let rust_only: Vec<_> = rust
        .keys()
        .filter(|k| !java.contains_key(*k))
        .cloned()
        .collect();
    let common: Vec<_> = jk
        .iter()
        .filter(|k| rust.contains_key(*k))
        .cloned()
        .collect();
    let mut first_info = None;
    for k in &common {
        let (jq, jinfo, jgt, jad, jdp, jgq, jpl, jf) = &java[k];
        let (rq, rinfo, rgt, rad, rdp, rgq, rpl, rf) = &rust[k];
        if genuine_info(jinfo, rinfo) {
            first_info = Some(k.clone());
            kv(
                "first_common",
                format!(
                    "{}:{} {}/{} java_PL={jpl} rust_PL={rpl} java_INFO_DP={} rust_INFO_DP={} GT {jgt}/{rgt} AD {jad}/{rad}",
                    k.0,
                    k.1,
                    k.2,
                    k.3,
                    jinfo.get("DP").cloned().unwrap_or_default(),
                    rinfo.get("DP").cloned().unwrap_or_default()
                ),
            );
            let _ = (jq, jdp, jgq, jf, rq, rdp, rgq, rf, rad, rgt);
            break;
        }
    }
    kv("java_records", java.len().to_string());
    kv("rust_records", rust.len().to_string());
    kv("java_only", java_only.len().to_string());
    kv("rust_only", rust_only.len().to_string());
    assert!(!java.contains_key(&("20".into(), 29_455_314, "G".into(), "C".into())));
    assert!(!rust.contains_key(&("20".into(), 29_455_314, "G".into(), "C".into())));
    assert_eq!(java.len(), 126);
    assert_eq!(rust.len(), 126);
    assert_eq!(java_only.len(), 2);
    assert_eq!(rust_only.len(), 2);
    assert_eq!(
        first_info,
        Some(("20".into(), 29_455_649, "T".into(), "TGTTTG".into()))
    );
    assert_eq!(
        rust_only[0],
        ("20".into(), 29_456_196, "A".into(), "T".into())
    );
    assert_eq!(
        java_only[0],
        ("20".into(), 29_455_902, "G".into(), "A".into())
    );
    kv(
        "classification",
        "common-site INFO DP at 20:29455649; not covering G>C emit",
    );
    kv("production_change", "NONE");
    let _ = QUAL_TOL;
    let _ = INFO_PRINT;
}
