//! Markdown + JSON human reports.

use crate::cli::ReportArgs;
use crate::run;
use crate::types::{EquivResults, TruthMetrics};
use anyhow::Result;
use std::fs;
use std::path::Path;

pub fn report(args: ReportArgs) -> Result<i32> {
    let mut results = run::load_results(&args.results_dir)?;
    if let Some(thr) = args.f1_delta_threshold {
        for d in &mut results.f1_deltas {
            d.within_threshold = d.abs_delta <= thr;
        }
        results.manifest.f1_delta_threshold = thr;
        results.max_abs_delta = results
            .f1_deltas
            .iter()
            .map(|d| d.abs_delta)
            .fold(0.0_f64, f64::max);
        results.gate_passed = results.f1_deltas.iter().all(|d| d.within_threshold);
    }
    write_reports(&args.results_dir, &results)?;
    // Also refresh results.json with recomputed gate if threshold overridden.
    let json = serde_json::to_string_pretty(&results)?;
    fs::write(args.results_dir.join("results.json"), json)?;
    Ok(if results.gate_passed { 0 } else { 1 })
}

pub fn write_reports(out_dir: &Path, results: &EquivResults) -> Result<()> {
    let md = render_markdown(results);
    fs::write(out_dir.join("REPORT.md"), md)?;
    // Compact JSON already in results.json; also emit report.json alias.
    let json = serde_json::to_string_pretty(results)?;
    fs::write(out_dir.join("report.json"), json)?;
    Ok(())
}

fn render_markdown(r: &EquivResults) -> String {
    let mut s = String::new();
    s.push_str("# gatk-rs-equiv report\n\n");
    s.push_str(&format!(
        "- **Created:** {}\n- **Engine:** `{}`\n- **Interval:** {}\n- **F1 Δ threshold:** {:.4}\n- **Gate:** {}\n- **Max |ΔF1|:** {:.4}\n\n",
        r.manifest.created_utc,
        r.manifest.engine,
        r.manifest
            .interval
            .as_deref()
            .unwrap_or("(full BAM / BED)"),
        r.manifest.f1_delta_threshold,
        if r.gate_passed {
            "**PASS**"
        } else {
            "**FAIL**"
        },
        r.max_abs_delta
    ));

    s.push_str("## Equivalence metric (Rust F1 − Java F1)\n\n");
    s.push_str(
        "Absolute F1 against truth can be high for both tools while they still disagree. \
The gate uses **|ΔF1|** only.\n\n",
    );
    s.push_str(
        "| Stratum | Class | Java F1 | Rust F1 | Δ (Rust−Java) | |Δ| | Within threshold |\n",
    );
    s.push_str(
        "|---------|-------|--------:|--------:|--------------:|----:|:----------------:|\n",
    );
    for d in &r.f1_deltas {
        s.push_str(&format!(
            "| {} | {} | {:.4} | {:.4} | {:+.4} | {:.4} | {} |\n",
            d.stratum,
            d.class,
            d.java_f1,
            d.rust_f1,
            d.delta,
            d.abs_delta,
            if d.within_threshold { "yes" } else { "NO" }
        ));
    }
    s.push('\n');

    s.push_str("## Truth evaluation — Java\n\n");
    s.push_str(&metrics_table(&r.java_vs_truth));
    s.push_str("\n## Truth evaluation — Rust\n\n");
    s.push_str(&metrics_table(&r.rust_vs_truth));

    s.push_str("\n## Direct Rust vs Java (no truth)\n\n");
    let c = &r.direct_compare;
    s.push_str(&format!(
        "| Metric | Count |\n|--------|------:|\n\
         | Java sites (CHROM+POS+REF+ALT) | {} |\n\
         | Rust sites | {} |\n\
         | **Exact identical** (POS+REF+ALT+GT) | {} |\n\
         | Allele match, GT mismatch | {} |\n\
         | Same GT, FORMAT fields differ | {} |\n\
         | Java only | {} |\n\
         | Rust only | {} |\n\n",
        c.java_sites,
        c.rust_sites,
        c.identical_sites,
        c.allele_match_gt_mismatch,
        c.format_mismatch_same_gt,
        c.java_only,
        c.rust_only
    ));

    s.push_str("## Inputs\n\n");
    s.push_str(&format!(
        "- Reference: `{}`\n- BAM: `{}`\n- Truth VCF: `{}`\n- Confident BED: `{}`\n- Java VCF: `{}`\n- Rust VCF: `{}`\n",
        r.manifest.reference.display(),
        r.manifest.bam.display(),
        r.manifest.truth_vcf.display(),
        r.manifest.confident_regions.display(),
        r.manifest.java_vcf.display(),
        r.manifest.rust_vcf.display()
    ));
    if !r.manifest.stratification_beds.is_empty() {
        s.push_str("\n### Stratification BEDs\n\n");
        for (name, path) in &r.manifest.stratification_beds {
            s.push_str(&format!("- `{name}`: `{}`\n", path.display()));
        }
    }

    if !r.notes.is_empty() {
        s.push_str("\n## Notes\n\n");
        for n in &r.notes {
            s.push_str(&format!("- {n}\n"));
        }
    }

    s.push_str(
        "\n## Limits\n\n\
         This tool checks **output equivalence** via community truth engines and a direct site table. \
It does **not** measure runtime, memory, or full FORMAT/QUAL bitwise identity. \
See `gatk-rs-equiv/README.md`.\n",
    );
    s
}

fn metrics_table(rows: &[TruthMetrics]) -> String {
    let mut s = String::new();
    s.push_str("| Stratum | Class | Precision | Recall | F1 | TP | FN | FP |\n");
    s.push_str("|---------|-------|----------:|-------:|---:|---:|---:|---:|\n");
    for m in rows {
        for (class, p) in [("SNP", &m.snp), ("INDEL", &m.indel)] {
            s.push_str(&format!(
                "| {} | {class} | {:.4} | {:.4} | {:.4} | {} | {} | {} |\n",
                m.stratum, p.precision, p.recall, p.f1, p.truth_tp, p.truth_fn, p.query_fp
            ));
        }
    }
    s
}
