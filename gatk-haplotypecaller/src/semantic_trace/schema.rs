//! Semantic-trace event schema (`gatk_rs.hc.semantic_trace/v1`).

use serde::{Deserialize, Serialize};
use serde_json::Value;

/// Schema identifier written on every NDJSON line.
pub const SCHEMA_ID: &str = "gatk_rs.hc.semantic_trace/v1";

/// Implementation label for Rust emitters.
pub const TRACE_IMPL_RUST: &str = "rust";

/// Implementation label for Java / projected emitters.
pub const TRACE_IMPL_JAVA: &str = "java";

/// Ordered production checkpoints (earlier stages sort before later for first-divergence).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SemanticStage {
    /// Band-pass activity profile cut (pre-region materialization).
    ActivityProfile,
    /// Materialized active/inactive assembly-region boundary.
    ActiveRegion,
    /// Assembly graph / result-set metrics.
    AssemblyGraph,
    /// Reference haplotype path.
    ReferencePath,
    /// Candidate haplotypes after assembly/trim.
    CandidateHaplotypes,
    /// Read×haplotype PairHMM likelihoods.
    ReadLikelihoods,
    /// Site genotype likelihoods / PL.
    GenotypeLikelihoods,
    /// Inactive reference-confidence model path.
    InactiveRcm,
    /// VCF records emitted for the region.
    VcfEmission,
}

impl SemanticStage {
    /// Wire name matching serde `snake_case`.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::ActivityProfile => "activity_profile",
            Self::ActiveRegion => "active_region",
            Self::AssemblyGraph => "assembly_graph",
            Self::ReferencePath => "reference_path",
            Self::CandidateHaplotypes => "candidate_haplotypes",
            Self::ReadLikelihoods => "read_likelihoods",
            Self::GenotypeLikelihoods => "genotype_likelihoods",
            Self::InactiveRcm => "inactive_rcm",
            Self::VcfEmission => "vcf_emission",
        }
    }

    /// Pipeline order index (lower = earlier).
    pub fn order_index(self) -> u8 {
        match self {
            Self::ActivityProfile => 0,
            Self::ActiveRegion => 1,
            Self::AssemblyGraph => 2,
            Self::ReferencePath => 3,
            Self::CandidateHaplotypes => 4,
            Self::ReadLikelihoods => 5,
            Self::GenotypeLikelihoods => 6,
            Self::InactiveRcm => 7,
            Self::VcfEmission => 8,
        }
    }
}

/// Genomic region key shared across stages.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct RegionKey {
    pub contig: String,
    /// 1-based inclusive start (unpadded active/inactive span).
    pub start: u64,
    /// 1-based inclusive end.
    pub end: u64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub is_active: Option<bool>,
}

/// One NDJSON semantic checkpoint.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SemanticTraceEvent {
    pub schema: String,
    /// Monotonic sequence number within the process (per sink).
    pub seq: u64,
    /// `"rust"` or `"java"`.
    #[serde(rename = "impl")]
    pub impl_name: String,
    pub stage: SemanticStage,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub region: Option<RegionKey>,
    pub payload: Value,
}
