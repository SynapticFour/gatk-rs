//! Illumina hap.py wrapper and summary.csv parser.

use super::{prf_from_counts, run_checked, EvalInput};
use crate::types::{Prf, TruthMetrics};
use anyhow::{bail, Context, Result};
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

/// Run hap.py for one query callset.
pub fn run_happy(
    happy_bin: &Path,
    input: &EvalInput<'_>,
    query_label: &str,
) -> Result<Vec<TruthMetrics>> {
    if let Some(parent) = input.out_prefix.parent() {
        std::fs::create_dir_all(parent)?;
    }

    let mut cmd = Command::new(happy_bin);
    cmd.arg(input.truth_vcf)
        .arg(input.query_vcf)
        .arg("-r")
        .arg(input.reference)
        .arg("-f")
        .arg(input.confident_bed)
        .arg("-o")
        .arg(input.out_prefix)
        .arg("--threads")
        .arg(input.threads.to_string())
        .arg("--verbose");

    // Stratification: write a TSV for hap.py --stratification when beds are provided.
    let strat_tsv = write_stratification_tsv(input)?;
    if let Some(ref tsv) = strat_tsv {
        cmd.arg("--stratification").arg(tsv);
    }

    run_checked(&mut cmd, "hap.py")?;

    let summary = PathBuf::from(format!("{}.summary.csv", input.out_prefix.display()));
    if !summary.is_file() {
        // Some builds write without the intermediate dot form; try extended location.
        let alt = input.out_prefix.with_extension("summary.csv");
        if alt.is_file() {
            return parse_happy_summary_csv(&alt, query_label);
        }
        bail!(
            "hap.py did not produce summary.csv at {}",
            summary.display()
        );
    }
    parse_happy_summary_csv(&summary, query_label)
}

fn write_stratification_tsv(input: &EvalInput<'_>) -> Result<Option<PathBuf>> {
    if input.stratification.is_empty() {
        return Ok(None);
    }
    let path = PathBuf::from(format!("{}.stratification.tsv", input.out_prefix.display()));
    let mut f = File::create(&path)?;
    // hap.py stratification TSV: region_name <tab> bed_path
    for (name, bed) in input.stratification {
        if !bed.is_file() {
            bail!(
                "stratification BED not found for '{name}': {}",
                bed.display()
            );
        }
        writeln!(f, "{name}\t{}", bed.display())?;
    }
    Ok(Some(path))
}

/// Parse Illumina hap.py `*.summary.csv`.
pub fn parse_happy_summary_csv(path: &Path, query_label: &str) -> Result<Vec<TruthMetrics>> {
    let file = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let mut rdr = csv::ReaderBuilder::new()
        .flexible(true)
        .from_reader(BufReader::new(file));
    let headers = rdr
        .headers()
        .context("read hap.py summary headers")?
        .iter()
        .map(|s| s.to_string())
        .collect::<Vec<_>>();

    let idx =
        |name: &str| -> Option<usize> { headers.iter().position(|h| h.eq_ignore_ascii_case(name)) };

    let i_type = idx("Type").context("summary.csv missing Type")?;
    let i_filter = idx("Filter");
    let i_subtype = idx("Subtype")
        .or_else(|| idx("Subset"))
        .or_else(|| idx("QQ"));
    // Stratification column names vary; try common ones.
    let i_strat = idx("Subset")
        .or_else(|| idx("Stratifier"))
        .or_else(|| idx("Region"));

    let i_tp = idx("TRUTH.TP").context("summary.csv missing TRUTH.TP")?;
    let i_fn = idx("TRUTH.FN").context("summary.csv missing TRUTH.FN")?;
    let i_fp = idx("QUERY.FP").context("summary.csv missing QUERY.FP")?;
    let i_prec = idx("METRIC.Precision");
    let i_rec = idx("METRIC.Recall");
    let i_f1 = idx("METRIC.F1_Score");

    // Accumulate by (stratum, type)
    let mut snp: BTreeMap<String, Prf> = BTreeMap::new();
    let mut indel: BTreeMap<String, Prf> = BTreeMap::new();
    let mut all: BTreeMap<String, Prf> = BTreeMap::new();

    for rec in rdr.records() {
        let rec = rec?;
        let typ = rec.get(i_type).unwrap_or("").trim().to_ascii_uppercase();
        if let Some(i) = i_filter {
            let filt = rec.get(i).unwrap_or("").trim().to_ascii_uppercase();
            // Prefer PASS / ALL rows; skip others when Filter column exists.
            if !filt.is_empty() && filt != "PASS" && filt != "ALL" && filt != "*" {
                continue;
            }
        }
        // Skip per-genotype subtypes when Subtype is not * /.
        if let Some(i) = i_subtype {
            let sub = rec.get(i).unwrap_or("").trim();
            if !sub.is_empty() && sub != "*" && sub != "." && !sub.eq_ignore_ascii_case("ALL") {
                // Still allow stratification Subset column separately.
                if i_strat != Some(i) {
                    continue;
                }
            }
        }

        let stratum = i_strat
            .and_then(|i| rec.get(i))
            .map(|s| s.trim().to_string())
            .filter(|s| !s.is_empty())
            .unwrap_or_else(|| "*".to_string());

        let tp: u64 = parse_u64(rec.get(i_tp));
        let fn_: u64 = parse_u64(rec.get(i_fn));
        let fp: u64 = parse_u64(rec.get(i_fp));
        let mut prf = prf_from_counts(tp, fn_, fp);
        if let Some(i) = i_prec {
            if let Some(v) = parse_f64(rec.get(i)) {
                prf.precision = v;
            }
        }
        if let Some(i) = i_rec {
            if let Some(v) = parse_f64(rec.get(i)) {
                prf.recall = v;
            }
        }
        if let Some(i) = i_f1 {
            if let Some(v) = parse_f64(rec.get(i)) {
                prf.f1 = v;
            }
        }

        match typ.as_str() {
            "SNP" => {
                snp.insert(stratum, prf);
            }
            "INDEL" => {
                indel.insert(stratum, prf);
            }
            "ALL" | "*" => {
                all.insert(stratum, prf);
            }
            _ => {}
        }
    }

    let mut strata: BTreeMap<String, ()> = BTreeMap::new();
    for k in snp.keys().chain(indel.keys()).chain(all.keys()) {
        strata.insert(k.clone(), ());
    }
    if strata.is_empty() {
        strata.insert("*".into(), ());
    }

    let mut out = Vec::new();
    for stratum in strata.keys() {
        out.push(TruthMetrics {
            query_label: query_label.to_string(),
            stratum: stratum.clone(),
            snp: snp.get(stratum).cloned().unwrap_or_else(Prf::zero),
            indel: indel.get(stratum).cloned().unwrap_or_else(Prf::zero),
            all: all.get(stratum).cloned(),
        });
    }
    Ok(out)
}

fn parse_u64(s: Option<&str>) -> u64 {
    s.and_then(|x| x.trim().parse().ok()).unwrap_or(0)
}

fn parse_f64(s: Option<&str>) -> Option<f64> {
    s.and_then(|x| {
        let t = x.trim();
        if t.is_empty() || t == "." || t.eq_ignore_ascii_case("nan") {
            None
        } else {
            t.parse().ok()
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::tempdir;

    #[test]
    fn parses_minimal_summary() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("x.summary.csv");
        let mut f = File::create(&path).unwrap();
        writeln!(
            f,
            "Type,Filter,TRUTH.TP,TRUTH.FN,QUERY.FP,METRIC.Precision,METRIC.Recall,METRIC.F1_Score"
        )
        .unwrap();
        writeln!(f, "SNP,PASS,10,2,1,0.909,0.833,0.869").unwrap();
        writeln!(f, "INDEL,PASS,5,1,0,1.0,0.833,0.909").unwrap();
        let m = parse_happy_summary_csv(&path, "rust").unwrap();
        assert_eq!(m.len(), 1);
        assert!((m[0].snp.f1 - 0.869).abs() < 1e-6);
        assert!((m[0].indel.f1 - 0.909).abs() < 1e-6);
    }
}
