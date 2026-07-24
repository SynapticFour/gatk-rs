//! clap CLI for `gatk-rs-equiv`.

use clap::{Parser, Subcommand, ValueEnum};
use std::path::PathBuf;

#[derive(Debug, Clone, ValueEnum)]
pub enum EngineChoice {
    /// Prefer hap.py; fall back to RTG vcfeval if missing.
    Auto,
    /// Require Illumina hap.py (https://github.com/Illumina/hap.py).
    Happy,
    /// Require RTG Tools `vcfeval`.
    Vcfeval,
}

#[derive(Parser, Debug)]
#[command(
    name = "gatk-rs-equiv",
    version,
    about = "Scientific HaplotypeCaller equivalence: gatk-rs vs GATK4 via hap.py / RTG vcfeval",
    long_about = "Runs Java GATK4 HaplotypeCaller and gatk-rs HaplotypeCaller on identical inputs, \
evaluates both callsets with hap.py (or RTG vcfeval) against the same truth set, \
and reports Rust−Java F1 deltas plus a direct site comparison.\n\n\
Independent community tool — not affiliated with the Broad Institute. \
See gatk-rs-equiv/README.md for scope limits."
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Command,
}

#[derive(Subcommand, Debug)]
pub enum Command {
    /// Run both callers + equivalence engines; write results under --out.
    Run(RunArgs),
    /// Render Markdown + JSON report from a previous --out directory.
    Report(ReportArgs),
    /// Differential fuzzer: synthetic BAM → Java vs Rust → shrink → fixture (+ optional gh issue).
    DifferentialFuzz(DiffFuzzCliArgs),
}

#[derive(Parser, Debug)]
pub struct DiffFuzzCliArgs {
    /// Number of random scenarios to evaluate (M4 default: small).
    #[arg(long, default_value_t = 8)]
    pub iterations: u32,

    /// Work directory for intermediate BAM/VCF per iteration.
    #[arg(long = "out", default_value = "target/diff-fuzz")]
    pub out: PathBuf,

    /// Path to the gatk-rs binary.
    #[arg(long = "rust-binary", default_value = "target/release/gatk-rs")]
    pub rust_binary: PathBuf,

    /// Path to GATK4 package JAR (e.g. gatk-package-4.4.0.0-local.jar).
    #[arg(long = "java-gatk-jar")]
    pub java_gatk_jar: Option<PathBuf>,

    /// Alternative: `gatk` launcher script/binary on PATH or absolute path.
    #[arg(long = "java-gatk-bin")]
    pub java_gatk_bin: Option<PathBuf>,

    /// Where minimized regression fixtures are written.
    #[arg(
        long = "fixture-root",
        default_value = "gatk-haplotypecaller/tests/fixtures/regressions"
    )]
    pub fixture_root: PathBuf,

    /// Max |AD/DP| difference still treated as FORMAT-equal (default 0 = exact).
    #[arg(long = "format-ad-tol", default_value_t = 0)]
    pub format_ad_tol: u32,

    /// Create a GitHub issue with label `parity-divergence` for each fixture (`gh` required).
    #[arg(long = "open-github-issue", default_value_t = false)]
    pub open_github_issue: bool,

    /// Max shrink evaluation steps per divergence.
    #[arg(long = "shrink-steps", default_value_t = 24)]
    pub shrink_steps: usize,

    /// Fixed campaign seed (optional; default: wall-clock).
    #[arg(long = "seed")]
    pub seed: Option<u64>,

    /// Abort if fewer than this many GiB are free near --out.
    #[arg(long = "min-free-gb", default_value_t = crate::resources::DEFAULT_MIN_FREE_GB)]
    pub min_free_gb: u64,

    /// Replay a previously written fixture directory (skips campaign).
    #[arg(long = "replay-fixture")]
    pub replay_fixture: Option<PathBuf>,
}

#[derive(Parser, Debug)]
pub struct RunArgs {
    /// Path to GATK4 package JAR (e.g. gatk-package-4.4.0.0-local.jar).
    #[arg(long = "java-gatk-jar")]
    pub java_gatk_jar: Option<PathBuf>,

    /// Alternative: `gatk` launcher script/binary on PATH or absolute path.
    #[arg(long = "java-gatk-bin")]
    pub java_gatk_bin: Option<PathBuf>,

    /// Path to the gatk-rs binary.
    #[arg(long = "rust-binary", required = true)]
    pub rust_binary: PathBuf,

    /// Reference FASTA (indexed.fai required by both callers / hap.py).
    #[arg(long = "reference", required = true)]
    pub reference: PathBuf,

    /// Input BAM/CRAM (indexed).
    #[arg(long = "bam", required = true)]
    pub bam: PathBuf,

    /// Truth VCF (e.g. GIAB benchmark).
    #[arg(long = "truth-vcf", required = true)]
    pub truth_vcf: PathBuf,

    /// High-confidence / evaluation regions BED.
    #[arg(long = "confident-regions", required = true)]
    pub confident_regions: PathBuf,

    /// Output directory for VCFs, engine outputs, and reports.
    #[arg(long = "out", required = true)]
    pub out: PathBuf,

    /// Optional calling / evaluation interval (GATK `-L` syntax, e.g. 20:10000000-10050000).
    #[arg(long = "interval")]
    pub interval: Option<String>,

    /// Max allowed |rust_f1 − java_f1| for gate pass (default 0.02).
    #[arg(long = "f1-delta-threshold", default_value = "0.02")]
    pub f1_delta_threshold: f64,

    /// Equivalence engine selection.
    #[arg(long = "engine", value_enum, default_value_t = EngineChoice::Auto)]
    pub engine: EngineChoice,

    /// Path to hap.py executable (default: search PATH / HAPPY_BIN).
    #[arg(long = "happy-bin")]
    pub happy_bin: Option<PathBuf>,

    /// Path to `rtg` launcher (default: search PATH / RTG_BIN).
    #[arg(long = "rtg-bin")]
    pub rtg_bin: Option<PathBuf>,

    /// Extra stratification BED: NAME=path (repeatable). Name appears in reports.
    #[arg(long = "stratification-bed", value_name = "NAME=PATH")]
    pub stratification_bed: Vec<String>,

    /// Threads for hap.py / callers (clamped to 1..=4; default 2 for 16GB hosts).
    #[arg(long = "threads", default_value_t = crate::resources::DEFAULT_THREADS)]
    pub threads: u32,

    /// Abort if fewer than this many GiB are free near --out (default 8).
    #[arg(long = "min-free-gb", default_value_t = crate::resources::DEFAULT_MIN_FREE_GB)]
    pub min_free_gb: u64,

    /// Skip disk free-space check (not recommended on laptop-class hosts).
    #[arg(long = "skip-disk-check", default_value_t = false)]
    pub skip_disk_check: bool,

    /// Skip re-running callers if java.vcf / rust.vcf already exist under --out.
    #[arg(long = "reuse-vcfs", default_value_t = false)]
    pub reuse_vcfs: bool,

    /// Write report after run (default true).
    #[arg(long = "write-report", default_value_t = true)]
    pub write_report: bool,
}

#[derive(Parser, Debug)]
pub struct ReportArgs {
    /// Results directory produced by `gatk-rs-equiv run`.
    #[arg(long = "results-dir", required = true)]
    pub results_dir: PathBuf,

    /// Override F1 delta threshold when re-evaluating the gate (optional).
    #[arg(long = "f1-delta-threshold")]
    pub f1_delta_threshold: Option<f64>,
}
