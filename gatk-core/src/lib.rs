//! GATK-RS Core Library
//! This crate contains the core data structures and utilities for GATK-RS.
//! Prefer module-path imports (`gatk_core::io::…`, `gatk_core::types::…`) over
//! root barrel imports. Benchmarking/parallel scaffolding is available as modules
//! but is not re-exported at the crate root.
#![allow(ambiguous_glob_reexports)]
#![allow(clippy::result_large_err)]

// Allow in-crate tests/modules to reference this crate by name.
extern crate self as gatk_core;

pub mod io;
pub mod math;
pub mod memory;
pub mod reference;
pub mod types;
pub mod utils;
/// GATK-compatible VariantFiltration hard-filtering (INFO/QUAL JEXL subset).
pub mod variant_filtration;

/// Benchmark harness helpers (prefer depending on this module path explicitly).
pub mod benchmarking;
/// Incomplete / experimental parallel job scaffolding (not a stable product API).
pub mod parallel;

/// Test fixtures and helpers. Prefer `#[cfg(test)]` / dev-dependency usage over
/// embedding this module in production binaries.
#[doc(hidden)]
pub mod tests;
#[doc(hidden)]
pub use tests::integration_test_helpers as integration;

// Stable domain / utility re-exports (explicit, not benchmarking/parallel globs).
pub use math::*;
pub use memory::*;
pub use types::*;
pub use utils::*;

#[cfg(test)]
mod smoke_tests {
    #[test]
    fn test_basic_functionality() {
        let sum = 1 + 1;
        assert_eq!(sum, 2);
    }
}
