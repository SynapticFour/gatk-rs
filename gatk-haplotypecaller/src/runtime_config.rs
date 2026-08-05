//! Central runtime configuration parsed from environment variables (Sprint I-4).
//! Production algorithm modules must not call `std::env::var` directly — use this
//! module or [`crate::parity_harness`] (harness-gated flags only).
//! Inventory: `docs/ARCHITECTURE.md`.

/// Execution / CLI surface flags (always readable in production builds).
/// # Invariants
/// `legacy_provisional` is rejected at runtime when set (Sprint B removal).
/// `scaffold_output` and `activate_output_opt_out` are independent emit toggles.
/// # Ownership
/// [`Clone`] sub-config owned by [`RuntimeConfig`]; parsed once at process edges.
/// # Mutation
/// Immutable after [`RuntimeConfig::from_env`]; must not be mutated mid-run for parity.
/// # Biological assumptions
/// None documented (process control flags only).
/// # Java equivalence
/// Rust-native env surface; no direct Java argument collection (see `ENV_INVENTORY.md`).
#[derive(Debug, Clone, Default)]
pub struct ExecutionConfig {
    /// `GATK_RS_HC_SCAFFOLD_OUTPUT=1` — header-only VCF.
    pub scaffold_output: bool,
    /// `GATK_RS_HC_ACTIVATE_OUTPUT=0` — explicit opt-out of variant emit.
    pub activate_output_opt_out: bool,
    /// `GATK_RS_HC_LEGACY_PROVISIONAL=1` — rejected (Sprint B removal).
    pub legacy_provisional: bool,
    /// `GATK_RS_HC_SEQUENTIAL=1` — force assembly-region apply batch size 1 (Peak-RSS).
    pub sequential_regions: bool,
    /// `GATK_RS_HC_LARGE_REGION_READS` — flush a region alone at/above this read count.
    /// `None` → caller default (`run.rs` `LARGE_REGION_READS_SEQUENTIAL_DEFAULT`).
    pub large_region_reads: Option<usize>,
}

/// Diagnostics that must never change emit/genotype results.
/// # Invariants
/// Debug flags affect logging/tracing only; documented as non-behavior-changing.
/// `debug_tile_overlaps` may be O(tiles×BAM) and is for scaffold diagnostics only.
/// # Ownership
/// [`Clone`] sub-config owned by [`RuntimeConfig`].
/// # Mutation
/// Immutable after env parse; helpers may re-read env on each call for dump-only paths.
/// # Biological assumptions
/// None documented.
/// # Java equivalence
/// Rust-native diagnostics env vars; not mirrored in GATK CLI.
#[derive(Debug, Clone, Default)]
pub struct DebugConfig {
    /// `GATK_RS_STRICT_CLUSTER_DEBUG=1` — stderr tracing for cluster materialize.
    pub strict_cluster_debug: bool,
    /// `GATK_RS_HC_DEBUG_TILE_OVERLAPS=1` — O(tiles×BAM) scaffold overlap count (R4-1).
    pub debug_tile_overlaps: bool,
    /// `GATK_RS_SEMANTIC_TRACE=<path>` — observe-only NDJSON semantic checkpoints.
    /// Never changes genotype/emit results; see `semantic_trace` module.
    pub semantic_trace_path: Option<String>,
}

/// Full runtime config. Prefer [`RuntimeConfig::from_env`] at process edges.
/// # Invariants
/// Production algorithm modules must not call `std::env::var` directly (use this module).
/// Harness-gated flags remain in [`crate::parity_harness`], not here.
/// # Ownership
/// Owns [`ExecutionConfig`] and [`DebugConfig`]; cheap to clone at startup.
/// # Mutation
/// Immutable after construction except dump helpers that re-parse env for debug booleans.
/// # Biological assumptions
/// None documented.
/// # Java equivalence
/// Rust-native centralized env config (Sprint I-4); GATK uses CLI args instead.
#[derive(Debug, Clone, Default)]
pub struct RuntimeConfig {
    pub execution: ExecutionConfig,
    pub debug: DebugConfig,
}

impl RuntimeConfig {
    /// Parse production-facing env vars once. Harness flags stay in [`crate::parity_harness`].
    pub fn from_env() -> Self {
        Self {
            execution: ExecutionConfig {
                scaffold_output: env_truthy("GATK_RS_HC_SCAFFOLD_OUTPUT"),
                activate_output_opt_out: env_is("GATK_RS_HC_ACTIVATE_OUTPUT", &["0", "false"]),
                legacy_provisional: env_truthy("GATK_RS_HC_LEGACY_PROVISIONAL"),
                sequential_regions: env_truthy("GATK_RS_HC_SEQUENTIAL"),
                large_region_reads: std::env::var("GATK_RS_HC_LARGE_REGION_READS")
                    .ok()
                    .and_then(|s| s.parse().ok()),
            },
            debug: DebugConfig {
                strict_cluster_debug: env_truthy("GATK_RS_STRICT_CLUSTER_DEBUG"),
                debug_tile_overlaps: env_truthy("GATK_RS_HC_DEBUG_TILE_OVERLAPS"),
                semantic_trace_path: std::env::var("GATK_RS_SEMANTIC_TRACE")
                    .ok()
                    .filter(|p| !p.is_empty() && p != "0" && !p.eq_ignore_ascii_case("off")),
            },
        }
    }
}

fn env_truthy(name: &str) -> bool {
    std::env::var(name)
        .map(|v| matches!(v.as_str(), "1" | "true" | "TRUE" | "yes" | "YES"))
        .unwrap_or(false)
}

fn env_is(name: &str, values: &[&str]) -> bool {
    std::env::var(name).ok().is_some_and(|v| {
        values
            .iter()
            .any(|want| v == *want || v.eq_ignore_ascii_case(want))
    })
}

/// Process-wide debug helpers (cheap; reads env each call — dump/debug only).
pub fn strict_cluster_debug_enabled() -> bool {
    RuntimeConfig::from_env().debug.strict_cluster_debug
}

/// `GATK_RS_HC_SEQUENTIAL=1` — force one assembly region at a time (16 GiB hosts).
pub fn hc_force_sequential_regions() -> bool {
    RuntimeConfig::from_env().execution.sequential_regions
}

/// Region read-count at which apply flushes alone. Override via `GATK_RS_HC_LARGE_REGION_READS`.
pub fn large_region_reads_sequential(default: usize) -> usize {
    RuntimeConfig::from_env()
        .execution
        .large_region_reads
        .unwrap_or(default)
}
