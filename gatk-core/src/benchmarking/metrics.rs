use super::*;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::time::{Duration, Instant};

/// Typed benchmark metric value (duration, memory, count, or float).
/// # Invariants
/// Each variant carries a single scalar measurement; units are documented by collector key names.
/// # Ownership
/// `Copy`/`Clone` enum; no heap for numeric variants.
/// # Mutation
/// N/A (enum value).
/// # Biological assumptions
/// None (infrastructure).
/// # Java equivalence
/// None / Rust-native.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum MetricValue {
    DurationMs(u128),
    Memory(u64),
    Count(u64),
    Float(f64),
}

/// Mutable collector for named benchmark metrics and optional timing spans.
/// # Invariants
/// At most one active timer from [`MetricsCollector::start_timing`] until `stop_timing`.
/// Metric names are unique keys in the internal map (last write wins).
/// # Ownership
/// Owns metric map; exclusive mutable access required while timing.
/// # Mutation
/// All methods take `&mut self`; not thread-safe.
/// # Biological assumptions
/// None (infrastructure).
/// # Java equivalence
/// None / Rust-native.
pub struct MetricsCollector {
    started_at: Option<Instant>,
    metrics: HashMap<String, MetricValue>,
}

impl Default for MetricsCollector {
    fn default() -> Self {
        Self::new()
    }
}

impl MetricsCollector {
    pub fn new() -> Self {
        Self {
            started_at: None,
            metrics: HashMap::new(),
        }
    }

    pub fn start_timing(&mut self) {
        self.started_at = Some(Instant::now());
    }

    pub fn stop_timing(&mut self, name: &str) -> Duration {
        let elapsed = self
            .started_at
            .take()
            .map(|s| s.elapsed())
            .unwrap_or_default();
        self.metrics.insert(
            name.to_string(),
            MetricValue::DurationMs(elapsed.as_millis()),
        );
        elapsed
    }

    pub fn add_metric(&mut self, name: &str, value: MetricValue) {
        self.metrics.insert(name.to_string(), value);
    }
}

/// Summary statistics (mean/min/max/std dev) for a metric sample set.
/// # Invariants
/// Empty input yields default zeros via internal `stats` helper.
/// # Ownership
/// Plain scalars; serde-friendly clone.
/// # Mutation
/// Immutable output of analysis.
/// # Biological assumptions
/// None (infrastructure).
/// # Java equivalence
/// None / Rust-native.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct StatisticSet {
    pub mean: f64,
    pub min: f64,
    pub max: f64,
    pub std_dev: f64,
}

/// Aggregated performance statistics across multiple [`BenchmarkResult`] rows.
/// # Invariants
/// `success_rate` in `[0.0, 1.0]` when derived from analyzer.
/// # Ownership
/// Owns nested [`StatisticSet`] values; clone for reports.
/// # Mutation
/// Immutable product of [`PerformanceAnalyzer::analyze`].
/// # Biological assumptions
/// None (infrastructure).
/// # Java equivalence
/// None / Rust-native.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct PerformanceStats {
    pub execution_time_stats: StatisticSet,
    pub memory_usage_stats: StatisticSet,
    pub success_rate: f64,
}

/// Accumulates [`BenchmarkResult`] samples for statistical analysis.
/// # Invariants
/// Results appended in invocation order; empty set yields default stats.
/// # Ownership
/// Owns result vector; clone duplicates all stored runs.
/// # Mutation
/// `add_result` mutates; `analyze` is read-only.
/// # Biological assumptions
/// None (infrastructure).
/// # Java equivalence
/// None / Rust-native.
pub struct PerformanceAnalyzer {
    results: Vec<BenchmarkResult>,
}

impl Default for PerformanceAnalyzer {
    fn default() -> Self {
        Self::new()
    }
}

impl PerformanceAnalyzer {
    pub fn new() -> Self {
        Self {
            results: Vec::new(),
        }
    }

    pub fn add_result(&mut self, result: BenchmarkResult) {
        self.results.push(result);
    }

    pub fn analyze(&self) -> PerformanceStats {
        let times: Vec<f64> = self
            .results
            .iter()
            .map(|r| r.execution_time_seconds())
            .collect();
        let mems: Vec<f64> = self.results.iter().map(|r| r.peak_memory_mb()).collect();
        let success = if self.results.is_empty() {
            0.0
        } else {
            self.results.iter().filter(|r| r.success).count() as f64 / self.results.len() as f64
        };

        PerformanceStats {
            execution_time_stats: stats(&times),
            memory_usage_stats: stats(&mems),
            success_rate: success,
        }
    }
}

fn stats(values: &[f64]) -> StatisticSet {
    if values.is_empty() {
        return StatisticSet::default();
    }
    let mean = values.iter().sum::<f64>() / values.len() as f64;
    let min = values.iter().copied().fold(f64::INFINITY, f64::min);
    let max = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let var = values.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / values.len() as f64;
    StatisticSet {
        mean,
        min,
        max,
        std_dev: var.sqrt(),
    }
}

/// Boolean flag indicating regression on a single metric axis.
/// # Invariants
/// `detected` is set by [`RegressionDetector::detect_regression`].
/// # Ownership
/// Small cloneable struct.
/// # Mutation
/// Immutable after detection.
/// # Biological assumptions
/// None (infrastructure).
/// # Java equivalence
/// None / Rust-native.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RegressionFlag {
    pub detected: bool,
}

/// Regression summary comparing current stats to a baseline.
/// # Invariants
/// `overall_regression` is true if either time or memory regressed.
/// # Ownership
/// Owns nested flags; serde clone.
/// # Mutation
/// Immutable report object.
/// # Biological assumptions
/// None (infrastructure).
/// # Java equivalence
/// None / Rust-native.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RegressionReport {
    pub overall_regression: bool,
    pub execution_time_regression: RegressionFlag,
    pub memory_regression: RegressionFlag,
}

/// Detects performance regressions against a stored baseline within a threshold.
/// # Invariants
/// Regression when current mean exceeds baseline mean × `(1 + threshold_percent/100)`.
/// # Ownership
/// Owns baseline [`PerformanceStats`] copy.
/// # Mutation
/// Stateless detection on `&self`.
/// # Biological assumptions
/// None (infrastructure).
/// # Java equivalence
/// None / Rust-native.
pub struct RegressionDetector {
    baseline: PerformanceStats,
    threshold_percent: f64,
}

impl RegressionDetector {
    pub fn new(baseline: PerformanceStats, threshold_percent: f64) -> Self {
        Self {
            baseline,
            threshold_percent,
        }
    }

    pub fn detect_regression(&self, current: &PerformanceStats) -> RegressionReport {
        let t = 1.0 + self.threshold_percent / 100.0;
        let exec_reg =
            current.execution_time_stats.mean > self.baseline.execution_time_stats.mean * t;
        let mem_reg = current.memory_usage_stats.mean > self.baseline.memory_usage_stats.mean * t;
        RegressionReport {
            overall_regression: exec_reg || mem_reg,
            execution_time_regression: RegressionFlag { detected: exec_reg },
            memory_regression: RegressionFlag { detected: mem_reg },
        }
    }
}

/// Configurable performance SLO thresholds for benchmark validation.
/// # Invariants
/// Defaults: 300s max time, 4096 MiB max memory, 95% min success rate.
/// # Ownership
/// Small owned struct; clone to customize.
/// # Mutation
/// Public fields editable before passing to validator.
/// # Biological assumptions
/// None (infrastructure).
/// # Java equivalence
/// None / Rust-native.
#[derive(Debug, Clone)]
pub struct PerformanceTargets {
    pub max_execution_time_seconds: f64,
    pub max_memory_mb: f64,
    pub min_success_rate: f64,
}

impl Default for PerformanceTargets {
    fn default() -> Self {
        Self {
            max_execution_time_seconds: 300.0,
            max_memory_mb: 4096.0,
            min_success_rate: 0.95,
        }
    }
}

/// Pass/fail outcome for one performance target with actual vs target values.
/// # Invariants
/// `met` reflects comparator logic in [`PerformanceTargetValidator::validate`].
/// # Ownership
/// Plain scalars; serde-friendly.
/// # Mutation
/// Immutable result row.
/// # Biological assumptions
/// None (infrastructure).
/// # Java equivalence
/// None / Rust-native.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TargetResult {
    pub met: bool,
    pub actual: f64,
    pub target: f64,
}

/// Full validation report against [`PerformanceTargets`].
/// # Invariants
/// `all_targets_met` is conjunction of individual [`TargetResult::met`] flags.
/// # Ownership
/// Owns three [`TargetResult`] entries.
/// # Mutation
/// Immutable validation output.
/// # Biological assumptions
/// None (infrastructure).
/// # Java equivalence
/// None / Rust-native.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PerformanceValidation {
    pub all_targets_met: bool,
    pub execution_time_target: TargetResult,
    pub memory_target: TargetResult,
    pub success_rate_target: TargetResult,
}

/// Validates [`PerformanceStats`] against [`PerformanceTargets`].
/// # Invariants
/// Uses mean execution time and mean memory from stats; success rate compared directly.
/// # Ownership
/// Owns targets configuration.
/// # Mutation
/// Validation is read-only on `&self`.
/// # Biological assumptions
/// None (infrastructure).
/// # Java equivalence
/// None / Rust-native.
pub struct PerformanceTargetValidator {
    targets: PerformanceTargets,
}

impl PerformanceTargetValidator {
    pub fn new(targets: PerformanceTargets) -> Self {
        Self { targets }
    }

    pub fn validate(&self, stats: &PerformanceStats) -> PerformanceValidation {
        let exec = TargetResult {
            met: stats.execution_time_stats.mean <= self.targets.max_execution_time_seconds,
            actual: stats.execution_time_stats.mean,
            target: self.targets.max_execution_time_seconds,
        };
        let mem = TargetResult {
            met: stats.memory_usage_stats.mean <= self.targets.max_memory_mb,
            actual: stats.memory_usage_stats.mean,
            target: self.targets.max_memory_mb,
        };
        let success = TargetResult {
            met: stats.success_rate >= self.targets.min_success_rate,
            actual: stats.success_rate,
            target: self.targets.min_success_rate,
        };
        PerformanceValidation {
            all_targets_met: exec.met && mem.met && success.met,
            execution_time_target: exec,
            memory_target: mem,
            success_rate_target: success,
        }
    }
}
