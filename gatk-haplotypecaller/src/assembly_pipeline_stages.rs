//! Typed assembly → genotyping pipeline stages (Sprint **K-5**).
//! Documents the Java mutation sequence that `call_region` mirrors. The double
//! [`crate::read_event_discovery::sync_assembly_events_from_haplotype_cigars_with_harvest`]
//! call is intentional: EventMap must be current **before** `filterAlleles` and again
//! **after** haplotype set changes.

/// Stages in the post-assemble `call_region` path (Java `HaplotypeCallerEngine` order).
/// # Invariants
/// Stages are ordered markers for EventMap sync / filter / realign / evidence change; not a state machine enum stored on the engine.
/// # Ownership
/// [`Copy`] stage tag for traces and sync helpers.
/// # Mutation
/// Immutable discriminant passed into sync calls.
/// # Biological assumptions
/// None — pipeline orchestration marker.
/// # Java equivalence
/// Documents GATK `HaplotypeCallerEngine.callRegion` mutation order (Sprint K-5).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CallRegionAssemblyStage {
    /// Haplotypes assembled; EventMap materialized from CIGARs (pre-trim or post-trim).
    EventMapMaterialized,
    /// After `AlleleFilteringHC.filterAlleles` — haplotype set may have shrunk.
    AlleleFiltered,
    /// After `realignReadsToTheirBestHaplotype`.
    Realigned,
    /// After `changeEvidence` — likelihood evidence points at best haplotypes.
    EvidenceChanged,
}

impl CallRegionAssemblyStage {
    /// Human-readable Java analogue for traces/docs.
    pub const fn java_analogue(self) -> &'static str {
        match self {
            Self::EventMapMaterialized => {
                "AssemblyResultSet.regenerateVariationEvents / EventMap.buildEventMapsForHaplotypes"
            }
            Self::AlleleFiltered => "AlleleFilteringHC.filterAlleles",
            Self::Realigned => "AssemblyBasedCallerUtils.realignReadsToTheirBestHaplotype",
            Self::EvidenceChanged => "ReadLikelihoods.changeEvidence",
        }
    }
}

/// Why EventMap sync runs twice around allele filtering (do not “simplify” without Java proof).
pub const EVENT_MAP_SYNC_AROUND_FILTER_RATIONALE: &str = "\
Java: EventMap on haplotypes → filterAlleles → (haplotype set changes) → EventMap must be rebuilt \
before realign/changeEvidence. Rust mirrors this with sync(EventMapMaterialized) before filter and \
sync again after AlleleFiltered.";
