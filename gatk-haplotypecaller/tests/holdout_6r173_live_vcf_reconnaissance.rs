//! 6R.173 live: `2:92305634 G/T` reconnaissance. INFO DP was 3 vs 2; 6R.174 closed it.
//! Skipped unless `HOLDOUT_6R173=1`.
//!
//! ```text
//! HOLDOUT_6R173=1 cargo test -p gatk-haplotypecaller --test holdout_6r173_live_vcf_reconnaissance -- --nocapture --test-threads=1
//! ```

use std::path::{Path, PathBuf};

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R173\t{key}\t{}", value.as_ref());
}

fn info_dp(line: &str) -> Option<&str> {
    let f: Vec<&str> = line.split('\t').collect();
    if f.len() < 8 {
        return None;
    }
    f[7].split(';').find_map(|p| p.strip_prefix("DP="))
}

fn sample_field(line: &str, name: &str) -> String {
    let f: Vec<&str> = line.split('\t').collect();
    if f.len() < 10 {
        return String::new();
    }
    f[8].split(':')
        .zip(f[9].split(':'))
        .find(|(k, _)| *k == name)
        .map(|(_, v)| v.to_string())
        .unwrap_or_default()
}

fn find_record(path: &Path, pos: u64, reff: &str, alt: &str) -> Option<String> {
    let text = std::fs::read_to_string(path).ok()?;
    for line in text.lines() {
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 5 {
            continue;
        }
        if f[0] == "2"
            && f[1].parse::<u64>().ok() == Some(pos)
            && f[3] == reff
            && f[4].split(',').next() == Some(alt)
        {
            return Some(line.to_string());
        }
    }
    None
}

#[test]
fn holdout_6r173_live_vcf_reconnaissance() {
    if std::env::var("HOLDOUT_6R173").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R173=1");
        return;
    }
    let root = repo_root();
    let java = root.join("parity/reports/6r43/p12_snp_cluster/java.vcf");
    let rust = root.join("parity/reports/6r43/p12_snp_cluster/rust.vcf");
    if !java.is_file() || !rust.is_file() {
        eprintln!("skip: missing 6r43 VCFs");
        return;
    }
    kv(
        "java_pin",
        "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77",
    );
    kv("production_change", "NONE");
    kv("target", "2:92305634 G/T");
    let j = find_record(&java, 92_305_634, "G", "T").expect("java record");
    let r = find_record(&rust, 92_305_634, "G", "T").expect("rust record");
    kv("java_info_dp", info_dp(&j).unwrap_or(""));
    kv("rust_info_dp", info_dp(&r).unwrap_or(""));
    kv(
        "format",
        format!(
            "J GT={} AD={} DP={} R GT={} AD={} DP={}",
            sample_field(&j, "GT"),
            sample_field(&j, "AD"),
            sample_field(&j, "DP"),
            sample_field(&r, "GT"),
            sample_field(&r, "AD"),
            sample_field(&r, "DP")
        ),
    );
    assert_eq!(sample_field(&j, "GT"), "1/1");
    assert_eq!(sample_field(&j, "AD"), "0,2");
    assert_eq!(sample_field(&j, "DP"), sample_field(&r, "DP"));
    assert_eq!(info_dp(&j), Some("3"));
    // 6R.174: live INFO DP matches Java Coverage. Pre-fix this was 2 (FORMAT DP copy).
    assert_eq!(info_dp(&r), Some("3"));
}
