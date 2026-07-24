//! RTG Tools `vcfeval` fallback engine.

use super::{prf_from_counts, run_checked, EvalInput};
use crate::types::{Prf, TruthMetrics};
use anyhow::{bail, Context, Result};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

/// Run RTG vcfeval. Stratification BEDs are evaluated as separate vcfeval passes
/// (intersection with confident regions left to the caller-provided BED when possible).
pub fn run_vcfeval(
    rtg_bin: &Path,
    input: &EvalInput<'_>,
    query_label: &str,
) -> Result<Vec<TruthMetrics>> {
    let out_dir = PathBuf::from(format!("{}.vcfeval", input.out_prefix.display()));
    if out_dir.exists() {
        fs::remove_dir_all(&out_dir)?;
    }

    // Ensure SDF for reference (cached next to out_prefix).
    let sdf = PathBuf::from(format!("{}.sdf", input.reference.display()));
    if !sdf.exists() {
        let mut fmt = Command::new(rtg_bin);
        fmt.arg("format").arg("-o").arg(&sdf).arg(input.reference);
        run_checked(&mut fmt, "rtg format")?;
    }

    let mut metrics = Vec::new();
    metrics.push(run_one(
        rtg_bin,
        input,
        &out_dir,
        input.confident_bed,
        "*",
        query_label,
        &sdf,
    )?);

    for (name, bed) in input.stratification {
        let strat_out = PathBuf::from(format!("{}.strat_{name}", out_dir.display()));
        if strat_out.exists() {
            fs::remove_dir_all(&strat_out)?;
        }
        metrics.push(run_one(
            rtg_bin,
            input,
            &strat_out,
            bed,
            name,
            query_label,
            &sdf,
        )?);
    }
    Ok(metrics)
}

fn run_one(
    rtg_bin: &Path,
    input: &EvalInput<'_>,
    out_dir: &Path,
    bed: &Path,
    stratum: &str,
    query_label: &str,
    sdf: &Path,
) -> Result<TruthMetrics> {
    let mut cmd = Command::new(rtg_bin);
    cmd.arg("vcfeval")
        .arg("-b")
        .arg(input.truth_vcf)
        .arg("-c")
        .arg(input.query_vcf)
        .arg("-t")
        .arg(sdf)
        .arg("-e")
        .arg(bed)
        .arg("-o")
        .arg(out_dir)
        .arg("--threads")
        .arg(input.threads.to_string());
    run_checked(&mut cmd, "rtg vcfeval")?;

    let summary = out_dir.join("summary.txt");
    if !summary.is_file() {
        bail!("rtg vcfeval missing summary.txt in {}", out_dir.display());
    }
    parse_vcfeval_summary(&summary, query_label, stratum)
}

/// Parse RTG `summary.txt` (SNP / INDEL / Total rows).
pub fn parse_vcfeval_summary(
    path: &Path,
    query_label: &str,
    stratum: &str,
) -> Result<TruthMetrics> {
    let text = fs::read_to_string(path).with_context(|| format!("read {}", path.display()))?;
    let mut snp = Prf::zero();
    let mut indel = Prf::zero();
    let mut all = None;

    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('=') || line.starts_with("Threshold") {
            continue;
        }
        // Typical: None SNP 10 2 1... precision recall f-measure
        let cols: Vec<&str> = line.split_whitespace().collect();
        if cols.len() < 6 {
            continue;
        }
        // Find a type token.
        let type_idx = cols.iter().position(|c| {
            matches!(
                c.to_ascii_uppercase().as_str(),
                "SNP" | "INDEL" | "INDELS" | "TOTAL" | "ALL"
            )
        });
        let Some(ti) = type_idx else {
            continue;
        };
        let typ = cols[ti].to_ascii_uppercase();
        // After type: baseline/TP-style counts vary by RTG version.
        // Prefer trailing precision recall f-measure floats.
        let floats: Vec<f64> = cols.iter().filter_map(|c| c.parse::<f64>().ok()).collect();
        if floats.len() < 3 {
            continue;
        }
        let n = floats.len();
        let precision = floats[n - 3];
        let recall = floats[n - 2];
        let f1 = floats[n - 1];
        // Best-effort TP/FN/FP from integer columns near the type token.
        let ints: Vec<u64> = cols
            .iter()
            .skip(ti + 1)
            .filter_map(|c| c.parse::<u64>().ok())
            .collect();
        let (tp, fn_, fp) = match ints.as_slice() {
            [tp, fn_, fp, ..] => (*tp, *fn_, *fp),
            _ => (0, 0, 0),
        };
        let mut prf = prf_from_counts(tp, fn_, fp);
        prf.precision = precision;
        prf.recall = recall;
        prf.f1 = f1;
        match typ.as_str() {
            "SNP" => snp = prf,
            "INDEL" | "INDELS" => indel = prf,
            "TOTAL" | "ALL" => all = Some(prf),
            _ => {}
        }
    }

    Ok(TruthMetrics {
        query_label: query_label.to_string(),
        stratum: stratum.to_string(),
        snp,
        indel,
        all,
    })
}
