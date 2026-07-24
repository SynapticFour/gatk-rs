//! Direct Java↔Rust VCF comparison (allele + GT + FORMAT), independent of truth.

use crate::types::DirectCompare;
use anyhow::{Context, Result};
use flate2::read::GzDecoder;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct AlleleKey {
    chrom: String,
    pos: u64,
    ref_a: String,
    alt: String,
}

#[derive(Debug, Clone)]
struct SiteRec {
    gt: String,
    /// Remaining FORMAT field values (excluding GT), joined for equality.
    format_rest: String,
}

fn open_vcf(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let name = path
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();
    if name.ends_with(".gz") || name.ends_with(".bgz") {
        Ok(Box::new(BufReader::new(GzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

fn canon_chrom(c: &str) -> String {
    let t = c.trim();
    if let Some(rest) = t.strip_prefix("chr") {
        rest.to_string()
    } else {
        t.to_string()
    }
}

fn normalize_gt(gt: &str) -> String {
    // Collapse phased/unphased separators for identity; keep allele indices.
    gt.replace('|', "/")
}

fn parse_sites(path: &Path) -> Result<HashMap<AlleleKey, SiteRec>> {
    let reader = open_vcf(path)?;
    let mut out: HashMap<AlleleKey, SiteRec> = HashMap::new();
    for line in reader.lines() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 5 {
            continue;
        }
        let chrom = canon_chrom(cols[0]);
        let pos: u64 = match cols[1].parse() {
            Ok(p) => p,
            Err(_) => continue,
        };
        let ref_a = cols[3].to_ascii_uppercase();
        let alts = cols[4];
        if alts == "." {
            continue;
        }
        let format = cols.get(8).copied().unwrap_or("");
        let sample = cols.get(9).copied().unwrap_or("");
        let (gt, format_rest) = extract_gt_and_rest(format, sample);
        for alt in alts.split(',') {
            let alt = alt.to_ascii_uppercase();
            if alt == "*" || alt.is_empty() {
                continue;
            }
            let key = AlleleKey {
                chrom: chrom.clone(),
                pos,
                ref_a: ref_a.clone(),
                alt,
            };
            out.insert(
                key,
                SiteRec {
                    gt: normalize_gt(&gt),
                    format_rest: format_rest.clone(),
                },
            );
        }
    }
    Ok(out)
}

fn extract_gt_and_rest(format: &str, sample: &str) -> (String, String) {
    if format.is_empty() || sample.is_empty() {
        return (String::new(), String::new());
    }
    let keys: Vec<&str> = format.split(':').collect();
    let vals: Vec<&str> = sample.split(':').collect();
    let mut gt = String::new();
    let mut rest_keys = Vec::new();
    let mut rest_vals = Vec::new();
    for (i, k) in keys.iter().enumerate() {
        let v = vals.get(i).copied().unwrap_or(".");
        if *k == "GT" {
            gt = v.to_string();
        } else {
            rest_keys.push(*k);
            rest_vals.push(v);
        }
    }
    let format_rest = rest_keys
        .into_iter()
        .zip(rest_vals)
        .map(|(k, v)| format!("{k}={v}"))
        .collect::<Vec<_>>()
        .join(";");
    (gt, format_rest)
}

/// Compare two callset VCFs (exact FORMAT string equality).
pub fn compare_callsets(java_vcf: &Path, rust_vcf: &Path) -> Result<DirectCompare> {
    compare_callsets_with_ad_tol(java_vcf, rust_vcf, 0)
}

/// Compare callsets; AD/DP numeric fields may differ by up to `ad_tol` without counting
/// as a FORMAT mismatch. Other FORMAT keys still require exact equality.
pub fn compare_callsets_with_ad_tol(
    java_vcf: &Path,
    rust_vcf: &Path,
    ad_tol: u32,
) -> Result<DirectCompare> {
    let java = parse_sites(java_vcf)?;
    let rust = parse_sites(rust_vcf)?;

    let mut identical = 0u64;
    let mut gt_mismatch = 0u64;
    let mut format_mismatch = 0u64;
    let mut java_only = 0u64;
    let mut rust_only = 0u64;

    for (key, j) in &java {
        match rust.get(key) {
            None => java_only += 1,
            Some(r) => {
                if j.gt == r.gt {
                    if format_rest_equal(&j.format_rest, &r.format_rest, ad_tol) {
                        identical += 1;
                    } else {
                        format_mismatch += 1;
                    }
                } else {
                    gt_mismatch += 1;
                }
            }
        }
    }
    for key in rust.keys() {
        if !java.contains_key(key) {
            rust_only += 1;
        }
    }

    Ok(DirectCompare {
        java_sites: java.len() as u64,
        rust_sites: rust.len() as u64,
        identical_sites: identical,
        allele_match_gt_mismatch: gt_mismatch,
        format_mismatch_same_gt: format_mismatch,
        java_only,
        rust_only,
    })
}

fn format_rest_equal(a: &str, b: &str, ad_tol: u32) -> bool {
    if a == b {
        return true;
    }
    if ad_tol == 0 {
        return false;
    }
    let map_a = parse_format_map(a);
    let map_b = parse_format_map(b);
    let mut keys: std::collections::BTreeSet<&str> = map_a.keys().copied().collect();
    keys.extend(map_b.keys().copied());
    for k in keys {
        let va = map_a.get(k).copied().unwrap_or(".");
        let vb = map_b.get(k).copied().unwrap_or(".");
        if k == "AD" {
            if !int_list_within_tol(va, vb, ad_tol) {
                return false;
            }
        } else if k == "DP" {
            if !int_within_tol(va, vb, ad_tol) {
                return false;
            }
        } else if va != vb {
            return false;
        }
    }
    true
}

fn parse_format_map(s: &str) -> std::collections::HashMap<&str, &str> {
    let mut m = std::collections::HashMap::new();
    for part in s.split(';') {
        if part.is_empty() {
            continue;
        }
        if let Some((k, v)) = part.split_once('=') {
            m.insert(k, v);
        }
    }
    m
}

fn int_within_tol(a: &str, b: &str, tol: u32) -> bool {
    match (a.parse::<i64>(), b.parse::<i64>()) {
        (Ok(x), Ok(y)) => (x - y).unsigned_abs() <= u64::from(tol),
        _ => a == b,
    }
}

fn int_list_within_tol(a: &str, b: &str, tol: u32) -> bool {
    let pa: Vec<&str> = a.split(',').collect();
    let pb: Vec<&str> = b.split(',').collect();
    if pa.len() != pb.len() {
        return false;
    }
    pa.iter()
        .zip(pb.iter())
        .all(|(x, y)| int_within_tol(x, y, tol))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::tempdir;

    #[test]
    fn detects_identical_and_format_drift() {
        let dir = tempdir().unwrap();
        let java = dir.path().join("j.vcf");
        let rust = dir.path().join("r.vcf");
        let header =
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n";
        let mut jf = File::create(&java).unwrap();
        write!(
            jf,
            "{header}1\t100\t.\tA\tG\t.\tPASS\t.\tGT:AD:DP\t0/1:10,10:20\n1\t200\t.\tC\tT\t.\tPASS\t.\tGT:DP\t1/1:30\n"
        )
        .unwrap();
        let mut rf = File::create(&rust).unwrap();
        write!(
            rf,
            "{header}1\t100\t.\tA\tG\t.\tPASS\t.\tGT:AD:DP\t0/1:9,11:20\n1\t200\t.\tC\tT\t.\tPASS\t.\tGT:DP\t1/1:30\n1\t300\t.\tG\tA\t.\tPASS\t.\tGT\t0/1\n"
        )
        .unwrap();
        let c = compare_callsets(&java, &rust).unwrap();
        assert_eq!(c.identical_sites, 1); // pos 200
        assert_eq!(c.format_mismatch_same_gt, 1); // pos 100
        assert_eq!(c.rust_only, 1);
        assert_eq!(c.java_only, 0);

        let c_tol = compare_callsets_with_ad_tol(&java, &rust, 1).unwrap();
        assert_eq!(c_tol.identical_sites, 2); // AD within ±1
        assert_eq!(c_tol.format_mismatch_same_gt, 0);
        assert_eq!(c_tol.rust_only, 1);
    }
}
