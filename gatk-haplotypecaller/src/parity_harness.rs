//! Parity harness flags — opt-in env toggles for L3/L4 sign-off only.
//! Production `strict_java` builds must not set these. Integration tests that need
//! harness behavior must compile with `--features parity_harness`.
//! See [`docs/ARCHITECTURE.md`](../../docs/ARCHITECTURE.md) (harness env flags).

/// Env vars that change Rust behavior only under [`harness_env_allowed`].
pub const HARNESS_ENV_FLAGS: &[&str] = &[
    "P12_PHASE_E",
    "P12_BASELINE_EMIT_FILTER",
    "GATK_RS_P12_EVENT_REGISTRY",
    "GATK_RS_P12_ENSURE_BRIDGES",
    "P12_L4_JAVA_FORMAT",
    "GATK_RS_ENABLE_READ_SUPPLEMENT",
    "GATK_RS_ENABLE_REF_MOTIF",
    "GATK_RS_ENABLE_CLUSTER_INJECT",
    "GATK_RS_ASM8_ONLY",
    "GATK_RS_HC_GIVEN_VCF",
];

/// True when parity harness env vars may affect Rust behavior.
pub fn harness_env_allowed() -> bool {
    cfg!(any(test, feature = "parity_harness"))
}

pub fn env_flag_true(name: &str) -> bool {
    if !harness_env_allowed() {
        return false;
    }
    std::env::var(name)
        .ok()
        .is_some_and(|v| matches!(v.as_str(), "1" | "true" | "TRUE" | "yes" | "YES"))
}

pub fn env_flag_set(name: &str) -> bool {
    if !harness_env_allowed() {
        return false;
    }
    std::env::var(name)
        .ok()
        .is_some_and(|v| !v.is_empty() && v != "0")
}

/// Harness-only string env (e.g. `GATK_RS_HC_GIVEN_VCF` path).
pub fn env_string(name: &str) -> Option<String> {
    if !harness_env_allowed() {
        return None;
    }
    std::env::var(name)
        .ok()
        .filter(|v| !v.is_empty() && v != "0")
}

#[cfg(debug_assertions)]
pub fn warn_if_harness_flags_set() {
    for flag in HARNESS_ENV_FLAGS {
        if std::env::var(flag)
            .ok()
            .is_some_and(|v| !v.is_empty() && v != "0")
        {
            tracing::warn!(
                "Harness flag {flag} is set; ignored unless built with `--features parity_harness` \
                 (see docs/ARCHITECTURE.md)"
            );
        }
    }
}

#[cfg(not(debug_assertions))]
pub fn warn_if_harness_flags_set() {}
