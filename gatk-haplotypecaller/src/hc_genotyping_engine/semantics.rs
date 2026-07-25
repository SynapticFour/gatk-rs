//! Genotyping semantics — Sprint **K-2** single source of mode truth.
//! Maps onto [`crate::engine::CallRegionMode`] for the production pipeline. Prefer
//! [`GenotypingSemantics`] / [`HcGenotypingConfig::is_java_compatible`](crate::hc_genotyping_engine::HcGenotypingConfig::is_java_compatible)
//! over raw booleans.

/// How genotyping/emit should behave relative to GATK 4.4 Java HC.
/// # Invariants
/// [`Self::JavaCompatible`] is the only production mode on release builds.
/// Maps 1:1 from [`crate::engine::CallRegionMode`] for pipeline mode selection.
/// # Ownership
/// [`Copy`] enum stored on [`HcGenotypingConfig`](crate::hc_genotyping_engine::HcGenotypingConfig).
/// # Mutation
/// Selected at engine construction; not changed mid-region in production.
/// # Biological assumptions
/// Java-compatible mode disables read-bridge and parity-only genotype rescues.
/// # Java equivalence
/// Rust-native mode enum (Sprint K-2); `JavaCompatible` ↔ pinned GATK 4.4 HC behavior.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum GenotypingSemantics {
    /// Production: PairHMM + Java `calculateGenotypes` emit only (`CallRegionMode::StrictJava`).
    #[default]
    JavaCompatible,
    /// Parity/L4 experiments — not production.
    ParityExperimental,
    /// Legacy N1 read-bridge genotyping — not production.
    LegacyReadBridges,
}

impl GenotypingSemantics {
    pub const fn is_java_compatible(self) -> bool {
        matches!(self, Self::JavaCompatible)
    }

    /// Map from [`crate::engine::CallRegionMode`].
    pub fn from_call_region_mode(mode: crate::engine::CallRegionMode) -> Self {
        match mode {
            crate::engine::CallRegionMode::StrictJava => Self::JavaCompatible,
            #[cfg(any(test, feature = "parity_harness"))]
            #[allow(deprecated)]
            crate::engine::CallRegionMode::ParityAligned => Self::ParityExperimental,
            #[cfg(any(test, feature = "parity_harness"))]
            crate::engine::CallRegionMode::LegacyReadBridges => Self::LegacyReadBridges,
        }
    }
}
