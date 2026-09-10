//! Experimental SeqGraph k-best **resource** policy (6R.140).
//!
//! Separates requested **K** (completed paths / visit budget) from **frontier
//! resource safety**. This is instrumentation / architecture only.
//!
//! **No k-best memory contract has been selected** (6R.139). The production
//! default remains `legacy_1024`: `BinaryHeap.len() >= 1024` Peak-RSS proxy.
//! Do not treat env vars or diagnostic byte ceilings as a CLAIM_MATRIX claim.

use crate::runtime_config;
use gatk_common::{GatkError, GatkResult};

/// Historical production heap-object proxy. Do not raise or lower this value here.
pub const LEGACY_MAX_HEAP_PATHS: usize = 1_024;

/// Diagnostic-only generous byte ceiling for 6R.140 canonical experiments.
/// **Not** a product budget. Canonical copied frontier was ~839 KiB (6R.138).
pub const DIAGNOSTIC_GENEROUS_BYTE_BUDGET: usize = 8 * 1024 * 1024;

/// How much SeqGraph k-best frontier resource is permitted.
///
/// Independent of `requested_k`. Java has no analogue.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SeqKbestResourcePolicy {
    /// Exact current production refusal: `BinaryHeap.len() >= 1024`.
    Legacy1024,
    /// No frontier-resource refusal (path-edge / expansion Peak-RSS caps remain).
    /// Forensic/diagnostic only — must not become the production default.
    UnboundedDiagnostic,
    /// Experimental live-frontier byte ceiling. Accounting is experimental
    /// (see [`SeqKbestRunDiagnostics`] and `docs/parity/KBEST_RESOURCE_POLICY.md`).
    ByteBudget { max_frontier_bytes: usize },
}

impl SeqKbestResourcePolicy {
    pub fn name(self) -> &'static str {
        match self {
            Self::Legacy1024 => "legacy_1024",
            Self::UnboundedDiagnostic => "unbounded_diagnostic",
            Self::ByteBudget { .. } => "byte_budget",
        }
    }

    pub fn resource_limit_display(self) -> String {
        match self {
            Self::Legacy1024 => format!("heap_len<{LEGACY_MAX_HEAP_PATHS}"),
            Self::UnboundedDiagnostic => "none".into(),
            Self::ByteBudget { max_frontier_bytes } => {
                format!("frontier_bytes<{max_frontier_bytes}")
            }
        }
    }

    /// Parse experimental env. Unset → [`Self::Legacy1024`].
    /// Invalid names/budgets **fail closed** (no silent fallback).
    pub fn from_runtime() -> GatkResult<Self> {
        parse_policy(
            runtime_config::experimental_kbest_policy_raw().as_deref(),
            runtime_config::experimental_kbest_memory_budget_raw().as_deref(),
        )
    }
}

pub fn parse_policy(
    policy_raw: Option<&str>,
    budget_raw: Option<&str>,
) -> GatkResult<SeqKbestResourcePolicy> {
    match policy_raw {
        None => {
            if budget_raw.is_some() {
                return Err(GatkError::configuration(
                    "GATK_RS_EXPERIMENTAL_KBEST_MEMORY_BUDGET is set but GATK_RS_EXPERIMENTAL_KBEST_POLICY is unset; set policy=byte_budget or unset the budget (experimental; not a product contract)",
                ));
            }
            Ok(SeqKbestResourcePolicy::Legacy1024)
        }
        Some("legacy_1024") => {
            if budget_raw.is_some() {
                return Err(GatkError::configuration(
                    "GATK_RS_EXPERIMENTAL_KBEST_MEMORY_BUDGET is only valid with GATK_RS_EXPERIMENTAL_KBEST_POLICY=byte_budget",
                ));
            }
            Ok(SeqKbestResourcePolicy::Legacy1024)
        }
        Some("unbounded_diagnostic") => {
            if budget_raw.is_some() {
                return Err(GatkError::configuration(
                    "GATK_RS_EXPERIMENTAL_KBEST_MEMORY_BUDGET is only valid with GATK_RS_EXPERIMENTAL_KBEST_POLICY=byte_budget",
                ));
            }
            Ok(SeqKbestResourcePolicy::UnboundedDiagnostic)
        }
        Some("byte_budget") => {
            let Some(raw) = budget_raw else {
                return Err(GatkError::configuration(
                    "GATK_RS_EXPERIMENTAL_KBEST_POLICY=byte_budget requires GATK_RS_EXPERIMENTAL_KBEST_MEMORY_BUDGET=<positive integer bytes>",
                ));
            };
            let bytes: usize = raw.parse().map_err(|_| {
                GatkError::configuration(format!(
                    "GATK_RS_EXPERIMENTAL_KBEST_MEMORY_BUDGET={raw:?} is not a usize byte count"
                ))
            })?;
            if bytes == 0 {
                return Err(GatkError::configuration(
                    "GATK_RS_EXPERIMENTAL_KBEST_MEMORY_BUDGET must be > 0",
                ));
            }
            Ok(SeqKbestResourcePolicy::ByteBudget {
                max_frontier_bytes: bytes,
            })
        }
        Some(other) => Err(GatkError::configuration(format!(
            "GATK_RS_EXPERIMENTAL_KBEST_POLICY={other:?} is invalid; expected unset, legacy_1024, unbounded_diagnostic, or byte_budget (experimental; not a product contract)"
        ))),
    }
}

/// Structured diagnostics for one SeqGraph k-best run.
///
/// **Frontier allocation** fields are experimental layout accounting, not malloc RSS.
/// **Process RSS** is sampled around the operation and is **not** k-best-attributable
/// merely because it was recorded here (6R.138).
#[derive(Debug, Clone)]
pub struct SeqKbestRunDiagnostics {
    pub requested_k: usize,
    pub completed_paths: usize,
    pub policy: &'static str,
    pub resource_limit: String,
    pub resource_limit_hit: bool,
    /// Generated, not-yet-popped distinct partials (copied: equals heap entries;
    /// compact lazy: heap heads + deferred siblings).
    pub peak_logical_frontier: usize,
    pub peak_heap_entries: usize,
    /// Live heap items × `size_of::<HeapItem>()` + Σ `edges.capacity()×16`.
    /// Does **not** include malloc rounding, BinaryHeap buffer slack, result paths,
    /// or the SeqGraph. Experimental.
    pub peak_frontier_bytes: usize,
    /// Live `heap.len() × size_of::<PathState>()` (Vec header included; payload not).
    pub peak_pathstate_bytes: usize,
    /// Σ live `edges.capacity() × size_of::<(usize,usize)>()`.
    pub peak_edge_payload_bytes: usize,
    pub peak_edge_payload_len_bytes: usize,
    /// `BinaryHeap.capacity() × size_of::<HeapItem>()` (slack not in the byte-budget check).
    pub peak_heap_vec_capacity_bytes: usize,
    pub paths_refused: usize,
    pub expansions: usize,
    pub rss_before_mib: Option<f64>,
    pub rss_peak_mib: Option<f64>,
    pub rss_after_mib: Option<f64>,
    pub accounting: &'static str,
}

impl SeqKbestRunDiagnostics {
    pub fn to_log_line(&self) -> String {
        format!(
            "KBEST_RESOURCE_DIAG requested_k={} completed_paths={} policy={} resource_limit={} resource_limit_hit={} peak_logical_frontier={} peak_heap_entries={} peak_frontier_bytes={} peak_pathstate_bytes={} peak_edge_payload_bytes={} paths_refused={} expansions={} rss_before_MiB={} rss_peak_MiB={} rss_after_MiB={} accounting={} (frontier_allocation≠process_RSS; experimental; not a product contract)",
            self.requested_k,
            self.completed_paths,
            self.policy,
            self.resource_limit,
            self.resource_limit_hit,
            self.peak_logical_frontier,
            self.peak_heap_entries,
            self.peak_frontier_bytes,
            self.peak_pathstate_bytes,
            self.peak_edge_payload_bytes,
            self.paths_refused,
            self.expansions,
            fmt_rss(self.rss_before_mib),
            fmt_rss(self.rss_peak_mib),
            fmt_rss(self.rss_after_mib),
            self.accounting,
        )
    }
}

fn fmt_rss(v: Option<f64>) -> String {
    v.map(|x| format!("{x:.3}")).unwrap_or_else(|| "NA".into())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_unset_is_legacy_1024() {
        let p = parse_policy(None, None).expect("unset");
        assert_eq!(p, SeqKbestResourcePolicy::Legacy1024);
        assert_eq!(p.name(), "legacy_1024");
    }

    #[test]
    fn explicit_legacy_matches_default() {
        assert_eq!(
            parse_policy(Some("legacy_1024"), None).unwrap(),
            SeqKbestResourcePolicy::Legacy1024
        );
    }

    #[test]
    fn unbounded_parses() {
        assert_eq!(
            parse_policy(Some("unbounded_diagnostic"), None).unwrap(),
            SeqKbestResourcePolicy::UnboundedDiagnostic
        );
    }

    #[test]
    fn byte_budget_requires_positive_integer() {
        assert_eq!(
            parse_policy(Some("byte_budget"), Some("1048576")).unwrap(),
            SeqKbestResourcePolicy::ByteBudget {
                max_frontier_bytes: 1_048_576
            }
        );
        assert!(parse_policy(Some("byte_budget"), None).is_err());
        assert!(parse_policy(Some("byte_budget"), Some("0")).is_err());
        assert!(parse_policy(Some("byte_budget"), Some("nope")).is_err());
    }

    #[test]
    fn invalid_policy_fails_closed() {
        let err = parse_policy(Some("raise_1024"), None).unwrap_err();
        let s = err.to_string();
        assert!(s.contains("invalid"), "{s}");
        assert!(s.contains("legacy_1024"), "{s}");
    }

    #[test]
    fn budget_without_byte_policy_fails() {
        assert!(parse_policy(None, Some("100")).is_err());
        assert!(parse_policy(Some("legacy_1024"), Some("100")).is_err());
        assert!(parse_policy(Some("unbounded_diagnostic"), Some("100")).is_err());
    }
}
