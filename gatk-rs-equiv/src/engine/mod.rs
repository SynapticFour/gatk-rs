//! External equivalence engines: hap.py (preferred) and RTG vcfeval (fallback).

mod happy;
mod vcfeval;

use crate::cli::EngineChoice;
use crate::types::{Prf, TruthMetrics};
use anyhow::{bail, Context, Result};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::Command;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EngineKind {
    Happy,
    Vcfeval,
}

impl EngineKind {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Happy => "hap.py",
            Self::Vcfeval => "rtg-vcfeval",
        }
    }
}

pub struct EnginePaths {
    pub happy_bin: Option<PathBuf>,
    pub rtg_bin: Option<PathBuf>,
}

fn which(name: &str) -> Option<PathBuf> {
    std::env::var_os("PATH").and_then(|paths| {
        for dir in std::env::split_paths(&paths) {
            let candidate = dir.join(name);
            if candidate.is_file() {
                return Some(candidate);
            }
        }
        None
    })
}

fn resolve_happy(explicit: Option<&Path>) -> Option<PathBuf> {
    if let Some(p) = explicit {
        return Some(p.to_path_buf());
    }
    if let Ok(p) = std::env::var("HAPPY_BIN") {
        return Some(PathBuf::from(p));
    }
    which("hap.py").or_else(|| which("happy"))
}

fn resolve_rtg(explicit: Option<&Path>) -> Option<PathBuf> {
    if let Some(p) = explicit {
        return Some(p.to_path_buf());
    }
    if let Ok(p) = std::env::var("RTG_BIN") {
        return Some(PathBuf::from(p));
    }
    which("rtg")
}

pub fn select_engine(choice: &EngineChoice, paths: &EnginePaths) -> Result<(EngineKind, PathBuf)> {
    let happy = resolve_happy(paths.happy_bin.as_deref());
    let rtg = resolve_rtg(paths.rtg_bin.as_deref());
    match choice {
        EngineChoice::Happy => {
            let p = happy.context("hap.py not found (set --happy-bin or HAPPY_BIN)")?;
            Ok((EngineKind::Happy, p))
        }
        EngineChoice::Vcfeval => {
            let p = rtg.context("rtg not found (set --rtg-bin or RTG_BIN)")?;
            Ok((EngineKind::Vcfeval, p))
        }
        EngineChoice::Auto => {
            if let Some(p) = happy {
                Ok((EngineKind::Happy, p))
            } else if let Some(p) = rtg {
                eprintln!("[gatk-rs-equiv] hap.py not found; falling back to RTG vcfeval");
                Ok((EngineKind::Vcfeval, p))
            } else {
                bail!(
                    "Neither hap.py nor rtg found. Install Illumina hap.py \
                     (https://github.com/Illumina/hap.py) or RTG Tools, \
                     or use Dockerfile.equiv."
                );
            }
        }
    }
}

pub struct EvalInput<'a> {
    pub truth_vcf: &'a Path,
    pub query_vcf: &'a Path,
    pub reference: &'a Path,
    pub confident_bed: &'a Path,
    pub out_prefix: &'a Path,
    pub threads: u32,
    /// Optional named stratification BEDs (GIAB strata, etc.).
    pub stratification: &'a BTreeMap<String, PathBuf>,
}

/// Run the selected engine; return metrics including `*` (all) and named strata when available.
pub fn evaluate(
    kind: EngineKind,
    engine_bin: &Path,
    input: &EvalInput<'_>,
    query_label: &str,
) -> Result<Vec<TruthMetrics>> {
    match kind {
        EngineKind::Happy => happy::run_happy(engine_bin, input, query_label),
        EngineKind::Vcfeval => vcfeval::run_vcfeval(engine_bin, input, query_label),
    }
}

pub(crate) fn run_checked(cmd: &mut Command, label: &str) -> Result<()> {
    eprintln!("[gatk-rs-equiv] {label}: {cmd:?}");
    let status = cmd
        .status()
        .with_context(|| format!("failed to spawn {label}"))?;
    if !status.success() {
        bail!("{label} failed with {status}");
    }
    Ok(())
}

pub(crate) fn prf_from_counts(tp: u64, fn_: u64, fp: u64) -> Prf {
    let precision = if tp + fp == 0 {
        0.0
    } else {
        tp as f64 / (tp + fp) as f64
    };
    let recall = if tp + fn_ == 0 {
        0.0
    } else {
        tp as f64 / (tp + fn_) as f64
    };
    let f1 = if precision + recall == 0.0 {
        0.0
    } else {
        2.0 * precision * recall / (precision + recall)
    };
    Prf {
        precision,
        recall,
        f1,
        truth_tp: tp,
        truth_fn: fn_,
        query_fp: fp,
    }
}
