//! Core site INFO key names (AC, AN, AF, NS, DP).
//! Computation lives in [`crate::annotator::annotate_parity_v1_site`].

/// INFO keys for the HC core / ChromosomeCounts subset.
#[allow(dead_code)] // public key catalog; emission uses annotate_parity_v1_site today
pub const CORE_INFO_KEYS: &[&str] = &["AC", "AN", "AF", "NS", "DP"];
