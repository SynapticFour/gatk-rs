//! Configuration management for GATK-RS
//! This module handles configuration parsing and validation for GATK tools.
#![allow(clippy::result_large_err)]

use crate::error::{GatkError, GatkResult};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Main configuration structure for GATK-RS
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GatkConfig {
    /// Java options (for compatibility, parsed but ignored in Rust)
    pub java_options: Option<JavaOptions>,

    /// Spark options (for compatibility, parsed but ignored in Rust)
    pub spark_options: Option<SparkOptions>,

    /// Tool-specific configuration
    pub tool_config: ToolConfig,

    /// Global configuration
    pub global_config: GlobalConfig,
}

/// Java options configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct JavaOptions {
    /// Memory allocation (e.g., "-Xmx8G")
    pub memory: Option<String>,
    /// Garbage collector (e.g., "-XX:+UseParallelGC")
    pub garbage_collector: Option<String>,
    /// Additional JVM arguments
    pub additional_args: Vec<String>,
}

impl JavaOptions {
    /// Parse Java options from string
    pub fn parse(options_str: &str) -> GatkResult<Self> {
        let mut memory = None;
        let mut garbage_collector = None;
        let mut additional_args = Vec::new();

        let args: Vec<&str> = options_str.split_whitespace().collect();

        for arg in args {
            if let Some(stripped) = arg.strip_prefix("-Xmx") {
                memory = Some(stripped.to_string());
            } else if arg.starts_with("-XX:+Use") && arg.contains("GC") {
                garbage_collector = Some(arg.to_string());
            } else {
                additional_args.push(arg.to_string());
            }
        }

        Ok(Self {
            memory,
            garbage_collector,
            additional_args,
        })
    }

    /// Get memory in GB if specified
    pub fn memory_gb(&self) -> Option<f64> {
        self.memory.as_ref().and_then(|mem| {
            if mem.ends_with('G') || mem.ends_with('g') {
                mem[..mem.len() - 1].parse().ok()
            } else if mem.ends_with('M') || mem.ends_with('m') {
                mem[..mem.len() - 1]
                    .parse::<f64>()
                    .ok()
                    .map(|mb| mb / 1024.0)
            } else {
                None
            }
        })
    }
}

/// Spark options configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SparkOptions {
    /// Master URL (e.g., "local\[*\]", "yarn")
    pub master: Option<String>,
    /// Application name
    pub app_name: Option<String>,
    /// Number of executors
    pub num_executors: Option<u32>,
    /// Executor memory
    pub executor_memory: Option<String>,
    /// Driver memory
    pub driver_memory: Option<String>,
    /// Additional Spark arguments
    pub additional_args: Vec<String>,
}

impl SparkOptions {
    /// Parse Spark options from arguments
    pub fn parse(args: &[String]) -> GatkResult<Self> {
        let mut master = None;
        let mut app_name = None;
        let mut num_executors = None;
        let mut executor_memory = None;
        let mut driver_memory = None;
        let mut additional_args = Vec::new();

        let mut i = 0;
        while i < args.len() {
            match args[i].as_str() {
                "--master" if i + 1 < args.len() => {
                    master = Some(args[i + 1].clone());
                    i += 2;
                }
                "--name" if i + 1 < args.len() => {
                    app_name = Some(args[i + 1].clone());
                    i += 2;
                }
                "--num-executors" if i + 1 < args.len() => {
                    num_executors = args[i + 1].parse().ok();
                    i += 2;
                }
                "--executor-memory" if i + 1 < args.len() => {
                    executor_memory = Some(args[i + 1].clone());
                    i += 2;
                }
                "--driver-memory" if i + 1 < args.len() => {
                    driver_memory = Some(args[i + 1].clone());
                    i += 2;
                }
                _ => {
                    additional_args.push(args[i].clone());
                    i += 1;
                }
            }
        }

        Ok(Self {
            master,
            app_name,
            num_executors,
            executor_memory,
            driver_memory,
            additional_args,
        })
    }
}

/// Tool-specific configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ToolConfig {
    /// Tool name
    pub tool_name: String,

    /// Input files
    pub inputs: InputConfig,

    /// Output configuration
    pub outputs: OutputConfig,

    /// Tool-specific parameters
    pub parameters: HashMap<String, String>,
}

/// Input configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct InputConfig {
    /// Reference genome file
    pub reference: Option<String>,

    /// Input BAM/SAM files
    pub input_files: Vec<String>,

    /// Intervals to process
    pub intervals: Option<String>,

    /// Pedigree file (for family-based analysis)
    pub pedigree: Option<String>,

    /// Known sites for recalibration
    pub known_sites: Vec<String>,
}

/// Output configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OutputConfig {
    /// Output VCF file
    pub output_vcf: Option<String>,

    /// Output GVCF file
    pub output_gvcf: Option<String>,

    /// Output mode (VCF, GVCF, BP_RESOLUTION)
    pub output_mode: String,

    /// Emit reference confidence
    pub emit_ref_confidence: Option<String>,

    /// Create output variant index
    pub create_output_variant_index: bool,
}

/// Global configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GlobalConfig {
    /// Number of threads
    pub num_threads: Option<usize>,

    /// Memory limit in GB
    pub memory_limit: Option<f64>,

    /// Verbosity level
    pub verbosity: VerbosityLevel,

    /// Quiet mode
    pub quiet: bool,

    /// Validate inputs
    pub validate_inputs: bool,

    /// Create index files
    pub create_index: bool,

    /// Disable auto-index creation
    pub disable_auto_index: bool,
}

/// Verbosity levels
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum VerbosityLevel {
    Error,
    Warning,
    Info,
    Debug,
    Trace,
}

impl VerbosityLevel {
    /// Convert from string
    pub fn parse(s: &str) -> GatkResult<Self> {
        match s.to_lowercase().as_str() {
            "error" => Ok(VerbosityLevel::Error),
            "warning" => Ok(VerbosityLevel::Warning),
            "info" => Ok(VerbosityLevel::Info),
            "debug" => Ok(VerbosityLevel::Debug),
            "trace" => Ok(VerbosityLevel::Trace),
            _ => Err(GatkError::configuration(format!(
                "Invalid verbosity level: {}",
                s
            ))),
        }
    }

    /// Convert to string
    pub fn to_str(self) -> &'static str {
        match self {
            VerbosityLevel::Error => "ERROR",
            VerbosityLevel::Warning => "WARNING",
            VerbosityLevel::Info => "INFO",
            VerbosityLevel::Debug => "DEBUG",
            VerbosityLevel::Trace => "TRACE",
        }
    }

    /// Convert to tracing level
    pub fn to_tracing_level(self) -> tracing::Level {
        match self {
            VerbosityLevel::Error => tracing::Level::ERROR,
            VerbosityLevel::Warning => tracing::Level::WARN,
            VerbosityLevel::Info => tracing::Level::INFO,
            VerbosityLevel::Debug => tracing::Level::DEBUG,
            VerbosityLevel::Trace => tracing::Level::TRACE,
        }
    }
}

impl std::str::FromStr for VerbosityLevel {
    type Err = GatkError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::parse(s)
    }
}

impl Default for GlobalConfig {
    fn default() -> Self {
        Self {
            num_threads: None,
            memory_limit: None,
            verbosity: VerbosityLevel::Info,
            quiet: false,
            validate_inputs: true,
            create_index: true,
            disable_auto_index: false,
        }
    }
}

/// HaplotypeCaller specific configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HaplotypeCallerConfig {
    /// Minimum base quality score
    pub min_base_quality_score: u8,

    /// Minimum mapping quality
    pub min_mapping_quality: u32,

    /// Maximum number of alternate alleles
    pub max_alternate_alleles: u32,

    /// Stand call confidence threshold
    pub stand_call_confidence: f64,

    /// Stand emit confidence threshold
    pub stand_emit_confidence: f64,

    /// Use original base qualities
    pub original_base_qualities: bool,

    /// Don't use soft clipped bases
    pub dont_use_soft_clipped_bases: bool,

    /// PairHMM implementation
    pub pair_hmm: Option<String>,

    /// Enable paired-end assembly
    pub enable_paired_end_assembly: bool,

    /// Assembly region parameters
    pub assembly_region: AssemblyRegionConfig,
}

/// Assembly / active-region traversal parameters (GATK `AssemblyRegionArgumentCollection` + `BandPassActivityProfile` wiring).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AssemblyRegionConfig {
    /// Minimum assembly region size (`--min-assembly-region-size`, default 50 in GATK).
    pub min_assembly_region_size: u32,

    /// Maximum assembly region size (`--max-assembly-region-size`, default 300 in GATK).
    pub max_assembly_region_size: u32,

    /// Minimum probability for a locus to be considered active (`--active-probability-threshold`).
    #[serde(alias = "active_region_threshold")]
    pub active_prob_threshold: f64,

    /// Upper limit on probability-mass propagation for active/inactive boundaries (`--max-prob-propagation-distance`).
    #[serde(default = "serde_default_max_prob_propagation_distance")]
    pub max_prob_propagation_distance: u32,

    /// Cap on Gaussian band-pass half-width in bp (`BandPassActivityProfile.MAX_FILTER_SIZE`, fixed 50 in GATK HC).
    #[serde(default = "serde_default_activity_profile_max_filter_size")]
    pub activity_profile_max_filter_size: u32,

    /// Gaussian σ for the activity band-pass (`BandPassActivityProfile.DEFAULT_SIGMA`, fixed 17 in GATK HC).
    #[serde(default = "serde_default_activity_profile_sigma")]
    pub activity_profile_sigma: f64,

    /// Maximum trails per active region (extension hook; not part of GATK’s assembly-region argument bundle).
    pub max_trails_per_active_region: u32,

    /// Assembly region padding (`--assembly-region-padding`).
    pub assembly_region_padding: u32,
}

impl AssemblyRegionConfig {
    /// GATK `AssemblyRegionArgumentCollection.DEFAULT_ACTIVE_PROB_THRESHOLD`.
    pub const GATK_DEFAULT_ACTIVE_PROB_THRESHOLD: f64 = 0.002;
    /// GATK `AssemblyRegionArgumentCollection.DEFAULT_MAX_PROB_PROPAGATION_DISTANCE`.
    pub const GATK_DEFAULT_MAX_PROB_PROPAGATION_DISTANCE: u32 = 50;
    /// GATK `BandPassActivityProfile.MAX_FILTER_SIZE`.
    pub const GATK_ACTIVITY_PROFILE_MAX_FILTER_SIZE: u32 = 50;
    /// GATK `BandPassActivityProfile.DEFAULT_SIGMA`.
    pub const GATK_ACTIVITY_PROFILE_SIGMA: f64 = 17.0;
}

fn serde_default_max_prob_propagation_distance() -> u32 {
    AssemblyRegionConfig::GATK_DEFAULT_MAX_PROB_PROPAGATION_DISTANCE
}

fn serde_default_activity_profile_max_filter_size() -> u32 {
    AssemblyRegionConfig::GATK_ACTIVITY_PROFILE_MAX_FILTER_SIZE
}

fn serde_default_activity_profile_sigma() -> f64 {
    AssemblyRegionConfig::GATK_ACTIVITY_PROFILE_SIGMA
}

impl Default for HaplotypeCallerConfig {
    fn default() -> Self {
        Self {
            min_base_quality_score: 10,
            min_mapping_quality: 20,
            max_alternate_alleles: 6,
            stand_call_confidence: 30.0,
            stand_emit_confidence: 10.0,
            original_base_qualities: false,
            dont_use_soft_clipped_bases: false,
            pair_hmm: None,
            enable_paired_end_assembly: true,
            assembly_region: AssemblyRegionConfig::default(),
        }
    }
}

impl Default for AssemblyRegionConfig {
    fn default() -> Self {
        Self {
            min_assembly_region_size: 50,
            max_assembly_region_size: 300,
            active_prob_threshold: Self::GATK_DEFAULT_ACTIVE_PROB_THRESHOLD,
            max_prob_propagation_distance: Self::GATK_DEFAULT_MAX_PROB_PROPAGATION_DISTANCE,
            activity_profile_max_filter_size: Self::GATK_ACTIVITY_PROFILE_MAX_FILTER_SIZE,
            activity_profile_sigma: Self::GATK_ACTIVITY_PROFILE_SIGMA,
            max_trails_per_active_region: 3,
            assembly_region_padding: 100,
        }
    }
}

impl GatkConfig {
    /// Create a new configuration
    pub fn new(tool_name: String) -> Self {
        Self {
            java_options: None,
            spark_options: None,
            tool_config: ToolConfig {
                tool_name,
                inputs: InputConfig {
                    reference: None,
                    input_files: Vec::new(),
                    intervals: None,
                    pedigree: None,
                    known_sites: Vec::new(),
                },
                outputs: OutputConfig {
                    output_vcf: None,
                    output_gvcf: None,
                    output_mode: "VCF".to_string(),
                    emit_ref_confidence: None,
                    create_output_variant_index: true,
                },
                parameters: HashMap::new(),
            },
            global_config: GlobalConfig::default(),
        }
    }

    /// Create a new configuration with Java options
    pub fn with_java_options(tool_name: String, java_options_str: &str) -> GatkResult<Self> {
        let java_options = JavaOptions::parse(java_options_str)?;
        let mut config = Self::new(tool_name);
        config.java_options = Some(java_options);
        Ok(config)
    }

    /// Create a new configuration with Spark options
    pub fn with_spark_options(tool_name: String, spark_args: &[String]) -> GatkResult<Self> {
        let spark_options = SparkOptions::parse(spark_args)?;
        let mut config = Self::new(tool_name);
        config.spark_options = Some(spark_options);
        Ok(config)
    }

    /// Create a new configuration with both Java and Spark options
    pub fn with_all_options(
        tool_name: String,
        java_options_str: Option<&str>,
        spark_args: Option<&[String]>,
    ) -> GatkResult<Self> {
        let mut config = Self::new(tool_name);

        if let Some(java_str) = java_options_str {
            config.java_options = Some(JavaOptions::parse(java_str)?);
        }

        if let Some(spark) = spark_args {
            config.spark_options = Some(SparkOptions::parse(spark)?);
        }

        Ok(config)
    }

    /// Validate the configuration
    pub fn validate(&self) -> GatkResult<()> {
        // Validate tool name
        if self.tool_config.tool_name.is_empty() {
            return Err(GatkError::configuration("Tool name cannot be empty"));
        }

        // Validate inputs
        if self.tool_config.inputs.input_files.is_empty() {
            return Err(GatkError::configuration(
                "At least one input file is required",
            ));
        }

        // Validate outputs
        if self.tool_config.outputs.output_vcf.is_none()
            && self.tool_config.outputs.output_gvcf.is_none()
        {
            return Err(GatkError::configuration(
                "Output VCF or GVCF file is required",
            ));
        }

        // Validate output mode
        match self.tool_config.outputs.output_mode.as_str() {
            "VCF" | "GVCF" | "BP_RESOLUTION" => {}
            _ => {
                return Err(GatkError::configuration(format!(
                    "Invalid output mode: {}",
                    self.tool_config.outputs.output_mode
                )))
            }
        }

        // Validate Java options if present
        if let Some(ref java_opts) = self.java_options {
            if let Some(ref memory) = java_opts.memory {
                if !memory.starts_with(|c: char| c.is_ascii_digit())
                    || !memory.ends_with(['G', 'g', 'M', 'm'])
                {
                    return Err(GatkError::configuration(format!(
                        "Invalid Java memory format: {}",
                        memory
                    )));
                }
            }
        }

        // Validate Spark options if present
        if let Some(ref spark_opts) = self.spark_options {
            if let Some(ref num_executors) = spark_opts.num_executors {
                if *num_executors == 0 {
                    return Err(GatkError::configuration(
                        "Number of Spark executors must be greater than 0",
                    ));
                }
            }

            if let Some(ref master) = spark_opts.master {
                if master.is_empty() {
                    return Err(GatkError::configuration("Spark master URL cannot be empty"));
                }
            }
        }

        // Validate memory limit
        if let Some(memory) = self.global_config.memory_limit {
            if memory <= 0.0 {
                return Err(GatkError::configuration("Memory limit must be positive"));
            }
        }

        // Validate thread count
        if let Some(threads) = self.global_config.num_threads {
            if threads == 0 {
                return Err(GatkError::configuration(
                    "Thread count must be greater than 0",
                ));
            }
        }

        Ok(())
    }

    /// Get HaplotypeCaller configuration
    pub fn get_haplotypecaller_config(&self) -> GatkResult<HaplotypeCallerConfig> {
        if self.tool_config.tool_name != "HaplotypeCaller" {
            return Err(GatkError::configuration(
                "Not a HaplotypeCaller configuration",
            ));
        }

        let mut config = HaplotypeCallerConfig::default();

        // Parse parameters
        for (key, value) in &self.tool_config.parameters {
            match key.as_str() {
                "min_base_quality_score" => {
                    config.min_base_quality_score = value.parse().map_err(|_| {
                        GatkError::configuration(format!(
                            "Invalid min_base_quality_score: {}",
                            value
                        ))
                    })?;
                }
                "min_mapping_quality" => {
                    config.min_mapping_quality = value.parse().map_err(|_| {
                        GatkError::configuration(format!("Invalid min_mapping_quality: {}", value))
                    })?;
                }
                "max_alternate_alleles" => {
                    config.max_alternate_alleles = value.parse().map_err(|_| {
                        GatkError::configuration(format!(
                            "Invalid max_alternate_alleles: {}",
                            value
                        ))
                    })?;
                }
                "stand_call_confidence" => {
                    config.stand_call_confidence = value.parse().map_err(|_| {
                        GatkError::configuration(format!(
                            "Invalid stand_call_confidence: {}",
                            value
                        ))
                    })?;
                }
                "stand_emit_confidence" => {
                    config.stand_emit_confidence = value.parse().map_err(|_| {
                        GatkError::configuration(format!(
                            "Invalid stand_emit_confidence: {}",
                            value
                        ))
                    })?;
                }
                "original_base_qualities" => {
                    config.original_base_qualities = value.parse().map_err(|_| {
                        GatkError::configuration(format!(
                            "Invalid original_base_qualities: {}",
                            value
                        ))
                    })?;
                }
                "dont_use_soft_clipped_bases" => {
                    config.dont_use_soft_clipped_bases = value.parse().map_err(|_| {
                        GatkError::configuration(format!(
                            "Invalid dont_use_soft_clipped_bases: {}",
                            value
                        ))
                    })?;
                }
                "pair_hmm" => {
                    config.pair_hmm = Some(value.clone());
                }
                _ => {
                    // Unknown parameter - could log a warning
                }
            }
        }

        Ok(config)
    }

    /// Set a tool parameter
    pub fn set_parameter(&mut self, key: String, value: String) {
        self.tool_config.parameters.insert(key, value);
    }

    /// Get a tool parameter
    pub fn get_parameter(&self, key: &str) -> Option<&String> {
        self.tool_config.parameters.get(key)
    }

    /// Add input file
    pub fn add_input_file(&mut self, file: String) {
        self.tool_config.inputs.input_files.push(file);
    }

    /// Set reference file
    pub fn set_reference(&mut self, reference: String) {
        self.tool_config.inputs.reference = Some(reference);
    }

    /// Set output VCF file
    pub fn set_output_vcf(&mut self, output: String) {
        self.tool_config.outputs.output_vcf = Some(output);
    }

    /// Set output mode
    pub fn set_output_mode(&mut self, mode: String) {
        self.tool_config.outputs.output_mode = mode;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_config_creation() {
        let config = GatkConfig::new("HaplotypeCaller".to_string());
        assert_eq!(config.tool_config.tool_name, "HaplotypeCaller");
        assert_eq!(config.tool_config.outputs.output_mode, "VCF");
    }

    #[test]
    fn test_config_validation() {
        let mut config = GatkConfig::new("HaplotypeCaller".to_string());

        // Should fail - no input files
        assert!(config.validate().is_err());

        // Add input file
        config.add_input_file("test.bam".to_string());
        assert!(config.validate().is_err()); // Still fails - no output

        // Set output
        config.set_output_vcf("test.vcf".to_string());
        assert!(config.validate().is_ok());
    }

    #[test]
    fn test_verbosity_levels() {
        let level = VerbosityLevel::parse("info").unwrap();
        assert_eq!(level, VerbosityLevel::Info);
        assert_eq!(level.to_str(), "INFO");
        assert_eq!(level.to_tracing_level(), tracing::Level::INFO);
    }

    #[test]
    fn test_haplotypecaller_config() {
        let mut config = GatkConfig::new("HaplotypeCaller".to_string());
        config.add_input_file("test.bam".to_string());
        config.set_output_vcf("test.vcf".to_string());
        config.set_parameter("min_base_quality_score".to_string(), "20".to_string());

        let hc_config = config.get_haplotypecaller_config().unwrap();
        assert_eq!(hc_config.min_base_quality_score, 20);
    }

    #[test]
    fn assembly_region_config_deserializes_legacy_active_region_threshold_key() {
        let json = r#"{"min_assembly_region_size":50,"max_assembly_region_size":300,"active_region_threshold":0.01,"max_trails_per_active_region":3,"assembly_region_padding":100}"#;
        let c: AssemblyRegionConfig = serde_json::from_str(json).unwrap();
        assert_eq!(c.active_prob_threshold, 0.01);
        assert_eq!(
            c.max_prob_propagation_distance,
            AssemblyRegionConfig::GATK_DEFAULT_MAX_PROB_PROPAGATION_DISTANCE
        );
        assert_eq!(
            c.activity_profile_max_filter_size,
            AssemblyRegionConfig::GATK_ACTIVITY_PROFILE_MAX_FILTER_SIZE
        );
        assert_eq!(
            c.activity_profile_sigma,
            AssemblyRegionConfig::GATK_ACTIVITY_PROFILE_SIGMA
        );
    }

    #[test]
    fn assembly_region_config_serde_roundtrip_preserves_gatk_defaults() {
        let c = AssemblyRegionConfig::default();
        let v = serde_json::to_value(&c).unwrap();
        let c2: AssemblyRegionConfig = serde_json::from_value(v).unwrap();
        assert_eq!(c2.active_prob_threshold, 0.002);
        assert_eq!(c2.max_prob_propagation_distance, 50);
        assert_eq!(c2.activity_profile_max_filter_size, 50);
        assert_eq!(c2.activity_profile_sigma, 17.0);
    }
}
