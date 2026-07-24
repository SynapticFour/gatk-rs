//! Shared result types for equivalence runs and reports.

use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::path::PathBuf;

/// Precision / recall / F1 for one variant class.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct Prf {
    pub precision: f64,
    pub recall: f64,
    pub f1: f64,
    pub truth_tp: u64,
    pub truth_fn: u64,
    pub query_fp: u64,
}

impl Prf {
    pub fn zero() -> Self {
        Self {
            precision: 0.0,
            recall: 0.0,
            f1: 0.0,
            truth_tp: 0,
            truth_fn: 0,
            query_fp: 0,
        }
    }
}

/// Metrics for one query callset against truth (optionally one stratum).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TruthMetrics {
    pub query_label: String,
    pub stratum: String,
    pub snp: Prf,
    pub indel: Prf,
    /// Combined / ALL row when the engine provides it.
    pub all: Option<Prf>,
}

/// Direct Rust↔Java site comparison (independent of truth / hap.py).
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct DirectCompare {
    pub java_sites: u64,
    pub rust_sites: u64,
    /// Exact match on CHROM+POS+REF+ALT+GT.
    pub identical_sites: u64,
    /// Same CHROM+POS+REF+ALT but GT differs.
    pub allele_match_gt_mismatch: u64,
    /// Same CHROM+POS+REF+ALT+GT but other FORMAT fields differ.
    pub format_mismatch_same_gt: u64,
    /// In Java only (allele key).
    pub java_only: u64,
    /// In Rust only (allele key).
    pub rust_only: u64,
}

/// F1 delta: rust_f1 - java_f1 (negative means Rust behind Java).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct F1Delta {
    pub stratum: String,
    pub class: String,
    pub java_f1: f64,
    pub rust_f1: f64,
    pub delta: f64,
    pub abs_delta: f64,
    pub within_threshold: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EquivRunManifest {
    pub created_utc: String,
    pub tool_version: String,
    pub engine: String,
    pub reference: PathBuf,
    pub bam: PathBuf,
    pub truth_vcf: PathBuf,
    pub confident_regions: PathBuf,
    pub interval: Option<String>,
    pub f1_delta_threshold: f64,
    pub java_vcf: PathBuf,
    pub rust_vcf: PathBuf,
    pub java_happy_prefix: PathBuf,
    pub rust_happy_prefix: PathBuf,
    pub stratification_beds: BTreeMap<String, PathBuf>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EquivResults {
    pub manifest: EquivRunManifest,
    pub java_vs_truth: Vec<TruthMetrics>,
    pub rust_vs_truth: Vec<TruthMetrics>,
    pub f1_deltas: Vec<F1Delta>,
    pub direct_compare: DirectCompare,
    pub gate_passed: bool,
    pub max_abs_delta: f64,
    pub notes: Vec<String>,
}
