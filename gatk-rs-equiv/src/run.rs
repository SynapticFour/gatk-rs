//! Orchestrate callers + engines + gate.

use crate::cli::RunArgs;
use crate::compare::compare_callsets;
use crate::engine::{self, EnginePaths, EvalInput};
use crate::hc::{self, JavaGatk};
use crate::report;
use crate::types::{EquivResults, EquivRunManifest, F1Delta, TruthMetrics};
use anyhow::{bail, Context, Result};
use chrono::Utc;
use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

const MANIFEST_JSON: &str = "manifest.json";
const RESULTS_JSON: &str = "results.json";

pub fn run(args: RunArgs) -> Result<i32> {
    validate_inputs(&args)?;
    fs::create_dir_all(&args.out)?;
    if !args.skip_disk_check {
        crate::resources::require_free_gb(&args.out, args.min_free_gb)?;
    }
    let threads = crate::resources::clamp_threads(args.threads);
    if threads != args.threads {
        eprintln!(
            "[gatk-rs-equiv] clamping --threads {} → {threads} (max {})",
            args.threads,
            crate::resources::MAX_THREADS
        );
    }
    crate::resources::apply_process_env(threads);

    let java = JavaGatk::resolve(args.java_gatk_jar.clone(), args.java_gatk_bin.clone())?;
    let java_vcf = args.out.join("java.vcf");
    let rust_vcf = args.out.join("rust.vcf");

    if !(args.reuse_vcfs && java_vcf.is_file() && rust_vcf.is_file()) {
        hc::run_java_hc(
            &java,
            &args.reference,
            &args.bam,
            &java_vcf,
            args.interval.as_deref(),
            threads,
        )?;
        hc::run_rust_hc(
            &args.rust_binary,
            &args.reference,
            &args.bam,
            &rust_vcf,
            args.interval.as_deref(),
            threads,
        )?;
    } else {
        eprintln!("[gatk-rs-equiv] reusing existing java.vcf / rust.vcf");
    }

    let (engine_kind, engine_bin) = engine::select_engine(
        &args.engine,
        &EnginePaths {
            happy_bin: args.happy_bin.clone(),
            rtg_bin: args.rtg_bin.clone(),
        },
    )?;

    let stratification = parse_strat_beds(&args.stratification_bed)?;
    let eval_dir = args.out.join("eval");
    fs::create_dir_all(&eval_dir)?;
    let java_prefix = eval_dir.join("java");
    let rust_prefix = eval_dir.join("rust");

    let java_vs_truth = engine::evaluate(
        engine_kind,
        &engine_bin,
        &EvalInput {
            truth_vcf: &args.truth_vcf,
            query_vcf: &java_vcf,
            reference: &args.reference,
            confident_bed: &args.confident_regions,
            out_prefix: &java_prefix,
            threads,
            stratification: &stratification,
        },
        "java",
    )?;
    let rust_vs_truth = engine::evaluate(
        engine_kind,
        &engine_bin,
        &EvalInput {
            truth_vcf: &args.truth_vcf,
            query_vcf: &rust_vcf,
            reference: &args.reference,
            confident_bed: &args.confident_regions,
            out_prefix: &rust_prefix,
            threads,
            stratification: &stratification,
        },
        "rust",
    )?;

    let direct_compare = compare_callsets(&java_vcf, &rust_vcf)?;
    let f1_deltas = compute_deltas(&java_vs_truth, &rust_vs_truth, args.f1_delta_threshold);
    let max_abs_delta = f1_deltas
        .iter()
        .map(|d| d.abs_delta)
        .fold(0.0_f64, f64::max);
    let gate_passed = f1_deltas.iter().all(|d| d.within_threshold);

    let manifest = EquivRunManifest {
        created_utc: Utc::now().to_rfc3339(),
        tool_version: env!("CARGO_PKG_VERSION").to_string(),
        engine: engine_kind.as_str().to_string(),
        reference: args.reference.clone(),
        bam: args.bam.clone(),
        truth_vcf: args.truth_vcf.clone(),
        confident_regions: args.confident_regions.clone(),
        interval: args.interval.clone(),
        f1_delta_threshold: args.f1_delta_threshold,
        java_vcf: java_vcf.clone(),
        rust_vcf: rust_vcf.clone(),
        java_happy_prefix: java_prefix,
        rust_happy_prefix: rust_prefix,
        stratification_beds: stratification,
    };

    let mut notes = vec![
        "Primary equivalence metric is Rust−Java F1 delta (not absolute F1).".into(),
        "Evaluation uses community engines (hap.py or RTG vcfeval), not in-house VCF heuristics."
            .into(),
    ];
    if !gate_passed {
        notes.push(format!(
            "Gate FAILED: max |ΔF1|={max_abs_delta:.4} exceeds threshold {}",
            args.f1_delta_threshold
        ));
    }

    let results = EquivResults {
        manifest,
        java_vs_truth,
        rust_vs_truth,
        f1_deltas,
        direct_compare,
        gate_passed,
        max_abs_delta,
        notes,
    };

    write_json(&args.out.join(MANIFEST_JSON), &results.manifest)?;
    write_json(&args.out.join(RESULTS_JSON), &results)?;

    if args.write_report {
        report::write_reports(&args.out, &results)?;
    }

    eprintln!(
        "[gatk-rs-equiv] done: gate_passed={} max_|ΔF1|={:.4} threshold={}",
        results.gate_passed, results.max_abs_delta, args.f1_delta_threshold
    );
    Ok(if results.gate_passed { 0 } else { 1 })
}

fn validate_inputs(args: &RunArgs) -> Result<()> {
    for (label, p) in [
        ("--reference", &args.reference),
        ("--bam", &args.bam),
        ("--truth-vcf", &args.truth_vcf),
        ("--confident-regions", &args.confident_regions),
        ("--rust-binary", &args.rust_binary),
    ] {
        if !p.exists() {
            bail!("{label} not found: {}", p.display());
        }
    }
    if args.f1_delta_threshold < 0.0 {
        bail!("--f1-delta-threshold must be >= 0");
    }
    Ok(())
}

fn parse_strat_beds(specs: &[String]) -> Result<BTreeMap<String, PathBuf>> {
    let mut out = BTreeMap::new();
    for spec in specs {
        let (name, path) = spec
            .split_once('=')
            .with_context(|| format!("--stratification-bed expected NAME=PATH, got {spec}"))?;
        let path = PathBuf::from(path);
        if !path.is_file() {
            bail!("stratification BED '{name}' not found: {}", path.display());
        }
        out.insert(name.to_string(), path);
    }
    Ok(out)
}

fn compute_deltas(java: &[TruthMetrics], rust: &[TruthMetrics], threshold: f64) -> Vec<F1Delta> {
    let mut out = Vec::new();
    for j in java {
        let r = rust.iter().find(|x| x.stratum == j.stratum);
        let Some(r) = r else { continue };
        for (class, jf, rf) in [
            ("SNP", j.snp.f1, r.snp.f1),
            ("INDEL", j.indel.f1, r.indel.f1),
        ] {
            let delta = rf - jf;
            let abs_delta = delta.abs();
            out.push(F1Delta {
                stratum: j.stratum.clone(),
                class: class.to_string(),
                java_f1: jf,
                rust_f1: rf,
                delta,
                abs_delta,
                within_threshold: abs_delta <= threshold,
            });
        }
    }
    out
}

fn write_json<T: serde::Serialize>(path: &Path, value: &T) -> Result<()> {
    let s = serde_json::to_string_pretty(value)?;
    fs::write(path, s)?;
    Ok(())
}

/// Load results.json for `report` subcommand.
pub fn load_results(results_dir: &Path) -> Result<EquivResults> {
    let path = results_dir.join(RESULTS_JSON);
    let text = fs::read_to_string(&path)
        .with_context(|| format!("missing {}; run `gatk-rs-equiv run` first", path.display()))?;
    Ok(serde_json::from_str(&text)?)
}
