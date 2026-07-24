use super::*;

/// Output serialization format for benchmark reports.
/// # Invariants
/// Each variant maps to a distinct file extension in [`BenchmarkReporter::generate_report`].
/// # Ownership
/// `Copy` enum.
/// # Mutation
/// N/A.
/// # Biological assumptions
/// None (infrastructure).
/// # Java equivalence
/// None / Rust-native.
#[derive(Debug, Clone, Copy)]
pub enum ReportFormat {
    Html,
    Json,
    Markdown,
    Csv,
}

/// Writes multi-format benchmark reports from a [`BenchmarkSuite`].
/// # Invariants
/// Creates `output_dir` if missing; one primary artifact per `generate_report` call.
/// # Ownership
/// Owns output directory path and format selection.
/// # Mutation
/// Report generation writes files; reporter instance reusable.
/// # Biological assumptions
/// None (infrastructure).
/// # Java equivalence
/// None / Rust-native.
pub struct BenchmarkReporter {
    output_dir: String,
    format: ReportFormat,
}

impl BenchmarkReporter {
    pub fn new(output_dir: String, format: ReportFormat) -> Self {
        Self { output_dir, format }
    }

    pub fn generate_report(&self, suite: &BenchmarkSuite) -> gatk_common::GatkResult<()> {
        std::fs::create_dir_all(&self.output_dir)
            .map_err(|e| gatk_common::GatkError::io("Failed to create report directory", e))?;

        let path = match self.format {
            ReportFormat::Html => format!("{}/report.html", self.output_dir),
            ReportFormat::Json => format!("{}/report.json", self.output_dir),
            ReportFormat::Markdown => format!("{}/report.md", self.output_dir),
            ReportFormat::Csv => format!("{}/report.csv", self.output_dir),
        };

        let body = match self.format {
            ReportFormat::Json => serde_json::to_string_pretty(suite)
                .map_err(|e| gatk_common::GatkError::generic(format!("Failed to serialize report: {e}")))?,
            ReportFormat::Markdown => format!(
                "# Benchmark Report\n\n- Comparisons: {}\n- Average speedup: {:.2}x\n- Average memory reduction: {:.2}x\n",
                suite.comparisons.len(),
                suite.average_speedup(),
                suite.average_memory_reduction()
            ),
            ReportFormat::Csv => {
                let mut out = "tool,time_s,memory_mb,success\n".to_string();
                for r in &suite.gatk_rs_results {
                    out.push_str(&format!(
                        "{},{:.3},{:.2},{}\n",
                        r.tool,
                        r.execution_time_seconds(),
                        r.peak_memory_mb(),
                        r.success
                    ));
                }
                out
            }
            ReportFormat::Html => format!(
                "<html><body><h1>Benchmark Report</h1><p>Comparisons: {}</p><p>Average speedup: {:.2}x</p><p>Average memory reduction: {:.2}x</p></body></html>",
                suite.comparisons.len(),
                suite.average_speedup(),
                suite.average_memory_reduction()
            ),
        };

        std::fs::write(path, body)
            .map_err(|e| gatk_common::GatkError::io("Failed to write report", e))?;
        Ok(())
    }

    pub fn generate_json_report(&self, suite: &BenchmarkSuite) -> gatk_common::GatkResult<()> {
        Self::new(self.output_dir.clone(), ReportFormat::Json).generate_report(suite)
    }

    pub fn generate_markdown_report(&self, suite: &BenchmarkSuite) -> gatk_common::GatkResult<()> {
        Self::new(self.output_dir.clone(), ReportFormat::Markdown).generate_report(suite)
    }

    pub fn generate_csv_report(&self, suite: &BenchmarkSuite) -> gatk_common::GatkResult<()> {
        Self::new(self.output_dir.clone(), ReportFormat::Csv).generate_report(suite)
    }
}

/// Writes a short executive summary markdown file for stakeholders.
/// # Invariants
/// Produces `{output_dir}/summary.md` with target/speedup highlights.
/// # Ownership
/// Owns output directory string.
/// # Mutation
/// File write on `generate_summary`.
/// # Biological assumptions
/// None (infrastructure).
/// # Java equivalence
/// None / Rust-native.
pub struct SummaryReporter {
    output_dir: String,
}

impl SummaryReporter {
    pub fn new(output_dir: String) -> Self {
        Self { output_dir }
    }

    pub fn generate_summary(&self, suite: &BenchmarkSuite) -> gatk_common::GatkResult<()> {
        std::fs::create_dir_all(&self.output_dir)
            .map_err(|e| gatk_common::GatkError::io("Failed to create summary directory", e))?;
        let path = format!("{}/summary.md", self.output_dir);
        let content = format!(
            "# Executive Summary\n\n- Meets all targets: {}\n- Average speedup: {:.2}x\n- Average memory reduction: {:.2}x\n",
            suite.meets_all_targets(),
            suite.average_speedup(),
            suite.average_memory_reduction()
        );
        std::fs::write(path, content)
            .map_err(|e| gatk_common::GatkError::io("Failed to write summary", e))?;
        Ok(())
    }
}
