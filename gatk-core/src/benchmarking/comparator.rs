use super::*;
use serde::{Deserialize, Serialize};

/// Roll-up statistics across many suite-to-suite benchmark comparisons.
/// # Invariants
/// Rates are fractions in `[0.0, 1.0]` when `total_comparisons > 0`.
/// `meets_all_targets` uses 2× speedup and 2× memory reduction heuristics.
/// # Ownership
/// Plain scalars; serde clone.
/// # Mutation
/// Immutable output of [`SuiteComparator::compare_suites`].
/// # Biological assumptions
/// None (infrastructure).
/// # Java equivalence
/// None / Rust-native.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OverallStats {
    pub total_comparisons: usize,
    pub average_speedup: f64,
    pub average_memory_reduction: f64,
    pub average_overall_score: f64,
    pub bitwise_identical_rate: f64,
    pub functional_equivalent_rate: f64,
    pub meets_all_targets: bool,
}

/// Comparison of two [`BenchmarkSuite`] instances with aggregated stats.
/// # Invariants
/// Pairs GATK and GATK-RS results by zip index; empty suites yield zeroed stats.
/// # Ownership
/// Owns [`OverallStats`]; cheap clone.
/// # Mutation
/// Immutable comparison product.
/// # Biological assumptions
/// None (infrastructure).
/// # Java equivalence
/// None / Rust-native.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SuiteComparison {
    pub overall_stats: OverallStats,
}

/// Stateless comparator for benchmark suite pairs.
/// # Invariants
/// Unit struct; no configuration state.
/// # Ownership
/// Zero-sized; use [`SuiteComparator::new`] for value.
/// # Mutation
/// N/A.
/// # Biological assumptions
/// None (infrastructure).
/// # Java equivalence
/// None / Rust-native.
pub struct SuiteComparator;

impl Default for SuiteComparator {
    fn default() -> Self {
        Self::new()
    }
}

impl SuiteComparator {
    pub fn new() -> Self {
        Self
    }

    pub fn compare_suites(
        &self,
        gatk: &BenchmarkSuite,
        gatk_rs: &BenchmarkSuite,
    ) -> SuiteComparison {
        let pairs = gatk
            .gatk_results
            .iter()
            .zip(gatk_rs.gatk_rs_results.iter())
            .collect::<Vec<_>>();
        let total = pairs.len();
        if total == 0 {
            return SuiteComparison {
                overall_stats: OverallStats {
                    total_comparisons: 0,
                    average_speedup: 0.0,
                    average_memory_reduction: 0.0,
                    average_overall_score: 0.0,
                    bitwise_identical_rate: 0.0,
                    functional_equivalent_rate: 0.0,
                    meets_all_targets: false,
                },
            };
        }

        let mut speedups = Vec::new();
        let mut mem_red = Vec::new();
        let mut scores = Vec::new();
        let mut functional = 0usize;
        for (a, b) in pairs {
            let c = ComparisonResult::new(a, b);
            speedups.push(c.speedup);
            mem_red.push(c.memory_reduction);
            scores.push((c.speedup.min(4.0) / 4.0 + c.memory_reduction.min(4.0) / 4.0) * 50.0);
            if c.both_succeeded {
                functional += 1;
            }
        }
        let avg = |v: &[f64]| v.iter().sum::<f64>() / v.len() as f64;
        SuiteComparison {
            overall_stats: OverallStats {
                total_comparisons: total,
                average_speedup: avg(&speedups),
                average_memory_reduction: avg(&mem_red),
                average_overall_score: avg(&scores),
                bitwise_identical_rate: 0.0,
                functional_equivalent_rate: functional as f64 / total as f64,
                meets_all_targets: avg(&speedups) >= 2.0 && avg(&mem_red) >= 2.0,
            },
        }
    }
}
