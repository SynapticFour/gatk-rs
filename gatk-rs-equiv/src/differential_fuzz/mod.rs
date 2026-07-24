//! Differential fuzzer: synthetic BAM → Java GATK4 vs gatk-rs → shrink → fixture + gh issue.
//! Seed decoding is shared with `fuzz/fuzz_targets/hc_differential.rs` via `#[path]`.

mod github;
mod runner;
mod scenario;
mod shrink;
mod synth;

pub use runner::{run_campaign, run_from_cli, DiffFuzzArgs, EvalBins};
pub use scenario::{scenario_from_bytes, scenario_to_seed_bytes, Scenario};
pub use synth::materialize_scenario;

/// Smoke entry: decode a scenario from raw bytes (never panics).
pub fn fuzz_scenario_smoke(data: &[u8]) {
    let scenario = scenario_from_bytes(data);
    let _ = scenario_to_seed_bytes(&scenario);
}
