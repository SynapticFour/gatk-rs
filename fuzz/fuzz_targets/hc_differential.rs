//! LibFuzzer entry for HC differential scenario decoding.
//! Shares the exact `Scenario` / `scenario_from_bytes` implementation with
//! `gatk-rs-equiv` via `#[path]` (no heavy bio deps in this binary).
//! Full Java↔Rust evaluation, shrinking, fixtures, and `gh` issues belong to:
//! `gatk-rs-equiv differential-fuzz` / `fuzz/run_hc_differential.sh`

#![no_main]

#[path = "../../gatk-rs-equiv/src/differential_fuzz/scenario.rs"]
mod scenario;

use libfuzzer_sys::fuzz_target;
use scenario::{scenario_from_bytes, scenario_to_seed_bytes};

fuzz_target!(|data: &[u8]| {
    let s = scenario_from_bytes(data);
    // Touch encode/decode so the optimizer cannot drop either path.
    let _ = scenario_to_seed_bytes(&s);
});
