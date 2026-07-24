use super::*;
use std::time::Duration;

/// Orchestrates repeated benchmark iterations and builds a [`BenchmarkSuite`].
/// # Invariants
/// Creates `config.output_dir` at construction; runs `iterations` paired GATK/GATK-RS stubs.
/// Current implementation uses placeholder timings (not live process execution).
/// # Ownership
/// Owns [`BenchmarkConfig`]; suite returned to caller on completion.
/// # Mutation
/// `run_benchmark_suite` mutates internal state only via filesystem side effects.
/// # Biological assumptions
/// None (infrastructure); stub variant counts are synthetic.
/// # Java equivalence
/// None / Rust-native (would invoke external GATK CLI in a full implementation).
pub struct BenchmarkRunner {
    config: BenchmarkConfig,
}

impl BenchmarkRunner {
    pub fn new(config: BenchmarkConfig) -> gatk_common::GatkResult<Self> {
        std::fs::create_dir_all(&config.output_dir).map_err(|e| {
            gatk_common::GatkError::io("Failed to create benchmark output directory", e)
        })?;
        Ok(Self { config })
    }

    pub fn run_benchmark_suite(&mut self) -> gatk_common::GatkResult<BenchmarkSuite> {
        let mut suite = BenchmarkSuite::new(self.config.clone());
        for i in 0..self.config.iterations {
            let gatk = BenchmarkResult::success(
                "GATK".to_string(),
                format!("iteration_{i}"),
                Duration::from_millis(1200),
                2 * 1024 * 1024 * 1024,
                70.0,
                10_000,
                20 * 1024 * 1024,
            );
            let gatk_rs = BenchmarkResult::success(
                "GATK-RS".to_string(),
                format!("iteration_{i}"),
                Duration::from_millis(700),
                1024 * 1024 * 1024,
                45.0,
                10_000,
                20 * 1024 * 1024,
            );
            suite.add_gatk_result(gatk);
            suite.add_gatk_rs_result(gatk_rs);
        }
        suite.generate_comparisons();
        Ok(suite)
    }
}
