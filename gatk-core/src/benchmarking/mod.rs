//! Benchmarking framework for GATK-RS performance comparison with original GATK

pub mod comparator;
pub mod datasets;
pub mod metrics;
pub mod reporter;
pub mod runner;

pub use comparator::*;
pub use datasets::*;
pub use metrics::*;
pub use reporter::*;
pub use runner::*;

use serde::{Deserialize, Serialize};
use std::time::Duration;

/// Benchmark run configuration (iterations, datasets, profiling flags).
/// # Invariants
/// `iterations` and `warmup_iterations` are counts, not durations.
/// `gatk_version` pins the Java baseline for comparisons.
/// # Ownership
/// Owns path/version strings; clone for suite copies.
/// # Mutation
/// Typically immutable for a suite after construction; fields are public for serde.
/// # Biological assumptions
/// None (infrastructure); `dataset` names benchmark fixtures, not biological types.
/// # Java equivalence
/// None / Rust-native benchmarking harness.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BenchmarkConfig {
    /// Number of iterations to run
    pub iterations: usize,
    /// Warmup iterations
    pub warmup_iterations: usize,
    /// Memory monitoring enabled
    pub memory_monitoring: bool,
    /// CPU profiling enabled
    pub cpu_profiling: bool,
    /// Output directory for results
    pub output_dir: String,
    /// Dataset to use for benchmarking
    pub dataset: String,
    /// GATK version to compare against
    pub gatk_version: String,
}

impl Default for BenchmarkConfig {
    fn default() -> Self {
        Self {
            iterations: 5,
            warmup_iterations: 2,
            memory_monitoring: true,
            cpu_profiling: false,
            output_dir: "benchmark_results".to_string(),
            dataset: "small".to_string(),
            gatk_version: "4.4.0.0".to_string(),
        }
    }
}

/// Single tool execution measurement from one benchmark iteration.
/// # Invariants
/// Successful runs have `success == true` and `error_message == None`.
/// Timings and memory are best-effort OS-level measurements.
/// # Ownership
/// Owns tool/command strings; clone for suite aggregation.
/// # Mutation
/// Built via constructors; fields public for serde roundtrips.
/// # Biological assumptions
/// `variant_count` summarizes caller output; not validated against truth set here.
/// # Java equivalence
/// None / Rust-native (compares against external GATK process metrics).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BenchmarkResult {
    /// Tool name (GATK-RS or GATK)
    pub tool: String,
    /// Command that was run
    pub command: String,
    /// Execution time
    pub execution_time: Duration,
    /// Peak memory usage in bytes
    pub peak_memory: u64,
    /// CPU usage percentage
    pub cpu_usage: f64,
    /// Number of variants produced
    pub variant_count: u64,
    /// Output file size in bytes
    pub output_size: u64,
    /// Success status
    pub success: bool,
    /// Error message if failed
    pub error_message: Option<String>,
}

impl BenchmarkResult {
    /// Create a new successful benchmark result
    pub fn success(
        tool: String,
        command: String,
        execution_time: Duration,
        peak_memory: u64,
        cpu_usage: f64,
        variant_count: u64,
        output_size: u64,
    ) -> Self {
        Self {
            tool,
            command,
            execution_time,
            peak_memory,
            cpu_usage,
            variant_count,
            output_size,
            success: true,
            error_message: None,
        }
    }

    /// Create a new failed benchmark result
    pub fn failure(tool: String, command: String, error_message: String) -> Self {
        Self {
            tool,
            command,
            execution_time: Duration::from_secs(0),
            peak_memory: 0,
            cpu_usage: 0.0,
            variant_count: 0,
            output_size: 0,
            success: false,
            error_message: Some(error_message),
        }
    }

    /// Get execution time in seconds
    pub fn execution_time_seconds(&self) -> f64 {
        self.execution_time.as_secs_f64()
    }

    /// Get peak memory in MB
    pub fn peak_memory_mb(&self) -> f64 {
        self.peak_memory as f64 / 1024.0 / 1024.0
    }

    /// Get output size in MB
    pub fn output_size_mb(&self) -> f64 {
        self.output_size as f64 / 1024.0 / 1024.0
    }
}

/// Aggregated speed/memory/output comparison between GATK and GATK-RS runs.
/// # Invariants
/// Ratios use GATK as baseline denominator where applicable.
/// `outputs_identical` defaults false until a comparator sets it.
/// # Ownership
/// Plain scalars; cheap clone.
/// # Mutation
/// Immutable result of [`ComparisonResult::new`].
/// # Biological assumptions
/// None (performance metadata); variant diff is count-based, not genotype concordance.
/// # Java equivalence
/// None / Rust-native.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ComparisonResult {
    /// Speedup factor (GATK time / GATK-RS time)
    pub speedup: f64,
    /// Memory reduction factor (GATK memory / GATK-RS memory)
    pub memory_reduction: f64,
    /// Output size difference percentage
    pub output_size_diff: f64,
    /// Variant count difference percentage
    pub variant_count_diff: f64,
    /// Both tools succeeded
    pub both_succeeded: bool,
    /// Outputs are identical (bitwise)
    pub outputs_identical: bool,
}

impl ComparisonResult {
    /// Create a new comparison result
    pub fn new(gatk_result: &BenchmarkResult, gatk_rs_result: &BenchmarkResult) -> Self {
        let speedup = if gatk_rs_result.execution_time_seconds() > 0.0 {
            gatk_result.execution_time_seconds() / gatk_rs_result.execution_time_seconds()
        } else {
            0.0
        };

        let memory_reduction = if gatk_rs_result.peak_memory_mb() > 0.0 {
            gatk_result.peak_memory_mb() / gatk_rs_result.peak_memory_mb()
        } else {
            0.0
        };

        let output_size_diff = if gatk_result.output_size_mb() > 0.0 {
            (gatk_rs_result.output_size_mb() - gatk_result.output_size_mb())
                / gatk_result.output_size_mb()
                * 100.0
        } else {
            0.0
        };

        let variant_count_diff = if gatk_result.variant_count > 0 {
            (gatk_rs_result.variant_count as f64 - gatk_result.variant_count as f64)
                / gatk_result.variant_count as f64
                * 100.0
        } else {
            0.0
        };

        Self {
            speedup,
            memory_reduction,
            output_size_diff,
            variant_count_diff,
            both_succeeded: gatk_result.success && gatk_rs_result.success,
            outputs_identical: false, // Will be set by comparator
        }
    }

    /// Check if GATK-RS meets performance targets
    pub fn meets_targets(&self) -> bool {
        self.speedup >= 2.0 && self.memory_reduction >= 2.0
    }

    /// Get performance summary
    pub fn summary(&self) -> String {
        format!(
            "Speedup: {:.2}x, Memory Reduction: {:.2}x, Output Size Diff: {:.1}%, Variant Count Diff: {:.1}%",
            self.speedup,
            self.memory_reduction,
            self.output_size_diff,
            self.variant_count_diff
        )
    }
}

/// Collection of paired GATK / GATK-RS benchmark runs plus metadata.
/// # Invariants
/// `gatk_results` and `gatk_rs_results` are paired by index for comparisons.
/// # Ownership
/// Owns vectors and nested config/metadata; clone or serialize for persistence.
/// # Mutation
/// Append results then call [`BenchmarkSuite::generate_comparisons`].
/// # Biological assumptions
/// None (infrastructure).
/// # Java equivalence
/// None / Rust-native.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BenchmarkSuite {
    /// Configuration used for benchmarks
    pub config: BenchmarkConfig,
    /// GATK benchmark results
    pub gatk_results: Vec<BenchmarkResult>,
    /// GATK-RS benchmark results
    pub gatk_rs_results: Vec<BenchmarkResult>,
    /// Comparison results
    pub comparisons: Vec<ComparisonResult>,
    /// Suite metadata
    pub metadata: SuiteMetadata,
}

/// Provenance and environment metadata for a [`BenchmarkSuite`].
/// # Invariants
/// `created_at` is UTC; git/rust/system fields best-effort at suite creation.
/// # Ownership
/// Owns strings and nested info structs.
/// # Mutation
/// Built once via [`SuiteMetadata::new`]; typically not mutated.
/// # Biological assumptions
/// None (infrastructure).
/// # Java equivalence
/// None / Rust-native.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SuiteMetadata {
    /// Timestamp when suite was created
    pub created_at: chrono::DateTime<chrono::Utc>,
    /// Git commit hash
    pub git_hash: String,
    /// Rust version
    pub rust_version: String,
    /// System information
    pub system_info: SystemInfo,
    /// GATK version used
    pub gatk_version: String,
    /// Dataset information
    pub dataset_info: DatasetInfo,
}

/// Host hardware/OS snapshot captured for benchmark reproducibility.
/// # Invariants
/// Memory fields are bytes; core counts from `num_cpus` / procfs when available.
/// # Ownership
/// Owns OS/CPU description strings.
/// # Mutation
/// Immutable after capture.
/// # Biological assumptions
/// None (infrastructure).
/// # Java equivalence
/// None / Rust-native.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SystemInfo {
    /// OS name and version
    pub os: String,
    /// CPU information
    pub cpu: String,
    /// Number of CPU cores
    pub cpu_cores: usize,
    /// Total memory in bytes
    pub total_memory: u64,
    /// Available memory in bytes
    pub available_memory: u64,
}

/// Description of a benchmark input dataset (reads + reference sizing).
/// # Invariants
/// `category`/`complexity` are taxonomy labels, not computed at runtime for unknown names.
/// # Ownership
/// Owns name/description strings and optional URL/checksum.
/// # Mutation
/// Serde-friendly public fields; normally loaded from fixtures.
/// # Biological assumptions
/// Summarizes sequencing experiment scale (read count, reference size).
/// # Java equivalence
/// None / Rust-native.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DatasetInfo {
    /// Dataset name
    pub name: String,
    /// Dataset size in bytes
    pub size: u64,
    /// Number of reads
    pub read_count: usize,
    /// Reference genome size
    pub reference_size: u64,
    /// Description
    pub description: String,
    /// Optional source URL
    pub url: Option<String>,
    /// Optional checksum
    pub checksum: Option<String>,
    /// Dataset creation timestamp
    pub created_at: chrono::DateTime<chrono::Utc>,
    /// Dataset category
    pub category: DatasetCategory,
    /// Dataset complexity
    pub complexity: DatasetComplexity,
}

/// Benchmark dataset provenance category.
/// # Invariants
/// Exhaustive enum; unknown filesystem datasets map to `Unknown` at discovery time.
/// # Ownership
/// `Copy` enum.
/// # Mutation
/// N/A.
/// # Biological assumptions
/// None (infrastructure taxonomy).
/// # Java equivalence
/// None / Rust-native.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum DatasetCategory {
    Synthetic,
    Public,
    Internal,
    Unknown,
}

/// Relative computational/biological difficulty label for a dataset.
/// # Invariants
/// Ordinal label only; no automatic scoring attached.
/// # Ownership
/// `Copy` enum.
/// # Mutation
/// N/A.
/// # Biological assumptions
/// Rough proxy for variant/read density in fixtures.
/// # Java equivalence
/// None / Rust-native.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum DatasetComplexity {
    Low,
    Medium,
    High,
}

impl BenchmarkSuite {
    /// Create a new benchmark suite
    pub fn new(config: BenchmarkConfig) -> Self {
        Self {
            metadata: SuiteMetadata::new(&config),
            config,
            gatk_results: Vec::new(),
            gatk_rs_results: Vec::new(),
            comparisons: Vec::new(),
        }
    }

    /// Add a GATK benchmark result
    pub fn add_gatk_result(&mut self, result: BenchmarkResult) {
        self.gatk_results.push(result);
    }

    /// Add a GATK-RS benchmark result
    pub fn add_gatk_rs_result(&mut self, result: BenchmarkResult) {
        self.gatk_rs_results.push(result);
    }

    /// Generate comparisons between results
    pub fn generate_comparisons(&mut self) {
        self.comparisons.clear();

        for (gatk_result, gatk_rs_result) in
            self.gatk_results.iter().zip(self.gatk_rs_results.iter())
        {
            let comparison = ComparisonResult::new(gatk_result, gatk_rs_result);
            self.comparisons.push(comparison);
        }
    }

    /// Get average speedup across all benchmarks
    pub fn average_speedup(&self) -> f64 {
        if self.comparisons.is_empty() {
            return 0.0;
        }

        let total: f64 = self.comparisons.iter().map(|c| c.speedup).sum();
        total / self.comparisons.len() as f64
    }

    /// Get average memory reduction across all benchmarks
    pub fn average_memory_reduction(&self) -> f64 {
        if self.comparisons.is_empty() {
            return 0.0;
        }

        let total: f64 = self.comparisons.iter().map(|c| c.memory_reduction).sum();
        total / self.comparisons.len() as f64
    }

    /// Check if all benchmarks meet performance targets
    pub fn meets_all_targets(&self) -> bool {
        self.comparisons.iter().all(|c| c.meets_targets())
    }

    /// Save benchmark suite to file
    pub fn save_to_file<P: AsRef<std::path::Path>>(&self, path: P) -> gatk_common::GatkResult<()> {
        let json = serde_json::to_string_pretty(self).map_err(|e| {
            gatk_common::GatkError::io("Failed to serialize benchmark suite", e.into())
        })?;

        std::fs::write(path, json)
            .map_err(|e| gatk_common::GatkError::io("Failed to write benchmark suite", e))?;

        Ok(())
    }

    /// Load benchmark suite from file
    pub fn load_from_file<P: AsRef<std::path::Path>>(path: P) -> gatk_common::GatkResult<Self> {
        let json = std::fs::read_to_string(path)
            .map_err(|e| gatk_common::GatkError::io("Failed to read benchmark suite", e))?;

        serde_json::from_str(&json).map_err(|e| {
            gatk_common::GatkError::io("Failed to deserialize benchmark suite", e.into())
        })
    }
}

impl SuiteMetadata {
    /// Create new suite metadata
    pub fn new(config: &BenchmarkConfig) -> Self {
        Self {
            created_at: chrono::Utc::now(),
            git_hash: get_git_hash(),
            rust_version: get_rust_version(),
            system_info: get_system_info(),
            gatk_version: config.gatk_version.clone(),
            dataset_info: get_dataset_info(&config.dataset),
        }
    }
}

/// Get current git hash
fn get_git_hash() -> String {
    std::process::Command::new("git")
        .args(["rev-parse", "HEAD"])
        .output()
        .ok()
        .and_then(|output| String::from_utf8(output.stdout).ok())
        .unwrap_or_else(|| "unknown".to_string())
}

/// Get current Rust version
fn get_rust_version() -> String {
    std::process::Command::new("rustc")
        .args(["--version"])
        .output()
        .ok()
        .and_then(|output| String::from_utf8(output.stdout).ok())
        .unwrap_or_else(|| "unknown".to_string())
}

/// Get system information
fn get_system_info() -> SystemInfo {
    let os = std::env::consts::OS.to_string();
    let cpu_cores = num_cpus::get();

    // Get CPU info (Linux specific)
    let cpu = if cfg!(target_os = "linux") {
        std::fs::read_to_string("/proc/cpuinfo")
            .ok()
            .and_then(|content| {
                content
                    .lines()
                    .find(|line| line.starts_with("model name"))
                    .and_then(|line| line.split(':').nth(1))
                    .map(|s| s.trim().to_string())
            })
            .unwrap_or_else(|| "unknown".to_string())
    } else {
        "unknown".to_string()
    };

    // Get memory info (Linux specific)
    let (total_memory, available_memory) = if cfg!(target_os = "linux") {
        if let Ok(meminfo) = std::fs::read_to_string("/proc/meminfo") {
            let mut total = 0u64;
            let mut available = 0u64;

            for line in meminfo.lines() {
                if line.starts_with("MemTotal:") {
                    if let Some(value) = extract_memory_value(line) {
                        total = value * 1024; // Convert KB to bytes
                    }
                } else if line.starts_with("MemAvailable:") {
                    if let Some(value) = extract_memory_value(line) {
                        available = value * 1024; // Convert KB to bytes
                    }
                }
            }

            (total, available)
        } else {
            (0, 0)
        }
    } else {
        (0, 0)
    };

    SystemInfo {
        os,
        cpu,
        cpu_cores,
        total_memory,
        available_memory,
    }
}

/// Extract memory value from /proc/meminfo line
fn extract_memory_value(line: &str) -> Option<u64> {
    line.split_whitespace()
        .nth(1)
        .and_then(|s| s.parse::<u64>().ok())
}

/// Get dataset information
fn get_dataset_info(dataset_name: &str) -> DatasetInfo {
    // This would be expanded to load actual dataset information
    match dataset_name {
        "small" => DatasetInfo {
            name: "small".to_string(),
            size: 100 * 1024 * 1024, // 100MB
            read_count: 1_000_000,
            reference_size: 3 * 1024 * 1024 * 1024, // 3GB
            description: "Small test dataset for quick benchmarking".to_string(),
            url: None,
            checksum: None,
            created_at: chrono::Utc::now(),
            category: DatasetCategory::Synthetic,
            complexity: DatasetComplexity::Low,
        },
        "medium" => DatasetInfo {
            name: "medium".to_string(),
            size: 1024 * 1024 * 1024, // 1GB
            read_count: 10_000_000,
            reference_size: 3 * 1024 * 1024 * 1024, // 3GB
            description: "Medium dataset for typical benchmarking".to_string(),
            url: None,
            checksum: None,
            created_at: chrono::Utc::now(),
            category: DatasetCategory::Synthetic,
            complexity: DatasetComplexity::Medium,
        },
        "large" => DatasetInfo {
            name: "large".to_string(),
            size: 10 * 1024 * 1024 * 1024, // 10GB
            read_count: 100_000_000,
            reference_size: 3 * 1024 * 1024 * 1024, // 3GB
            description: "Large dataset for stress testing".to_string(),
            url: None,
            checksum: None,
            created_at: chrono::Utc::now(),
            category: DatasetCategory::Synthetic,
            complexity: DatasetComplexity::High,
        },
        _ => DatasetInfo {
            name: dataset_name.to_string(),
            size: 0,
            read_count: 0,
            reference_size: 0,
            description: "Unknown dataset".to_string(),
            url: None,
            checksum: None,
            created_at: chrono::Utc::now(),
            category: DatasetCategory::Unknown,
            complexity: DatasetComplexity::Low,
        },
    }
}
