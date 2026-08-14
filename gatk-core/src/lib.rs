//! GATK-RS Core Library
//! This crate contains the core data structures and utilities for GATK-RS.
//! Prefer module-path imports (`gatk_core::io::…`, `gatk_core::types::…`) over
//! root barrel imports.
#![allow(clippy::result_large_err)]

// Allow in-crate tests/modules to reference this crate by name.
extern crate self as gatk_core;

pub mod io;
pub mod memory;
pub mod reference;
pub mod types;
pub mod utils;
/// GATK-compatible VariantFiltration hard-filtering (INFO/QUAL JEXL subset).
pub mod variant_filtration;

/// Test fixtures and helpers. Prefer `#[cfg(test)]` / dev-dependency usage over
/// embedding this module in production binaries.
#[doc(hidden)]
pub mod tests;
#[doc(hidden)]
pub use tests::integration_test_helpers as integration;

pub use memory::MemoryMappedFile;
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
