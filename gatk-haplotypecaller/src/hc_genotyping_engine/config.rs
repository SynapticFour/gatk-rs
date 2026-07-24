//! [`HcGenotypingConfig`] — production vs parity genotyping knobs.
//! Split out in Sprint I-2a for navigability. Sprint **K-2**: [`GenotypingSemantics`] is the
//! mode source of truth.

use crate::genotyping::BiallelicDiploidPriorModel;

pub use super::semantics::GenotypingSemantics;

/// Default GATK `informativeReadOverlapMargin` (bases).
pub const DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN: i32 = 2;
/// Default GATK `standard-min-confidence-threshold-for-calling` (GQ/QUAL gate).
pub const DEFAULT_STAND_EMIT_CONFIDENCE: f64 = 10.0;

/// HC genotyping configuration for [`super::genotype_active_region`].
/// # Invariants
/// [`Self::strict_java`] sets [`GenotypingSemantics::JavaCompatible`] with bridges/ rescue flags off.
/// `stand_emit_confidence` matches GATK standard min GQ/QUAL threshold (default 10).
/// # Ownership
/// Cloneable config snapshot threaded through `call_region` and genotype engines.
/// # Mutation
/// Callers clone/adjust before a region; engines read immutably per region.
/// # Biological assumptions
/// Diploid biallelic priors unless multiallelic paths extend GL vectors at emit time.
/// # Java equivalence
/// GATK `HaplotypeCallerArgumentCollection` genotyping slice + `assignGenotypeLikelihoods` / emit gates.
#[derive(Debug, Clone)]
pub struct HcGenotypingConfig {
    pub priors: BiallelicDiploidPriorModel,
    /// GATK `informativeReadOverlapMargin`.
    pub informative_read_overlap_margin: i32,
    /// GATK `disableSpanningEventGenotyping` (default false → spanning enabled).
    pub disable_spanning_event_genotyping: bool,
    /// Minimum GQ to keep a site (GATK `standard-min-confidence-threshold-for-calling`).
    pub stand_emit_confidence: f64,
    /// N1 bridge: genotype SNPs from read AD when no alt-hap support (off for Java parity).
    pub enable_sparse_read_genotype: bool,
    /// N1 bridge: relax GQ/emit from read depth (off for Java parity).
    pub enable_read_style_emit: bool,
    /// Parity/legacy: genotype `variation_events` list only (Java uses hap EventMap position walk).
    pub genotype_stored_events_only: bool,
    /// Sprint K-2: single source of genotyping mode truth.
    pub semantics: GenotypingSemantics,
    /// Parity/L4 experiments: P12 VCF-shaped GL rescue (never on [`Self::strict_java`]).
    pub enable_l4_emit_gl_rescue: bool,
}

impl Default for HcGenotypingConfig {
    fn default() -> Self {
        Self::strict_java()
    }
}

impl HcGenotypingConfig {
    /// Production Java-compatible genotyping (`GenotypingSemantics::JavaCompatible`).
    pub fn is_java_compatible(&self) -> bool {
        self.semantics.is_java_compatible()
    }

    /// Backward-compatible alias for [`Self::is_java_compatible`].
    #[inline]
    pub fn enable_java_strict(&self) -> bool {
        self.is_java_compatible()
    }

    /// GATK 4.4 — `assignGenotypeLikelihoods` EventMap walk + `calculateGenotypes` emit gate.
    pub fn strict_java() -> Self {
        Self {
            priors: BiallelicDiploidPriorModel::default(),
            informative_read_overlap_margin: DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN,
            disable_spanning_event_genotyping: false,
            stand_emit_confidence: DEFAULT_STAND_EMIT_CONFIDENCE,
            enable_sparse_read_genotype: false,
            enable_read_style_emit: false,
            genotype_stored_events_only: false,
            semantics: GenotypingSemantics::JavaCompatible,
            enable_l4_emit_gl_rescue: false,
        }
    }

    /// Deprecated alias — use [`Self::strict_java`] (no emit-GL rescue).
    #[deprecated(note = "use HcGenotypingConfig::strict_java()")]
    pub fn strict_java_emit() -> Self {
        Self::strict_java()
    }

    /// Parity/L4: optional GL rescue for FORMAT experiments (not Java).
    /// Sprint **L-4**: `cfg(test)` or `--features parity_harness` only.
    #[cfg(any(test, feature = "parity_harness"))]
    pub fn strict_java_l4() -> Self {
        Self {
            enable_l4_emit_gl_rescue: true,
            semantics: GenotypingSemantics::ParityExperimental,
            ..Self::parity_aligned()
        }
    }

    /// Parity/L4 experiments (not production).
    /// Sprint **L-4**: `cfg(test)` or `--features parity_harness` only.
    #[cfg(any(test, feature = "parity_harness"))]
    pub fn parity_aligned() -> Self {
        Self {
            semantics: GenotypingSemantics::ParityExperimental,
            ..Self::strict_java()
        }
    }

    /// Legacy N1 bridges (read pileup genotype + relaxed emit).
    /// Sprint **L-4**: `cfg(test)` or `--features parity_harness` only.
    #[cfg(any(test, feature = "parity_harness"))]
    pub fn legacy_read_bridges() -> Self {
        Self {
            enable_sparse_read_genotype: true,
            enable_read_style_emit: true,
            genotype_stored_events_only: false,
            semantics: GenotypingSemantics::LegacyReadBridges,
            ..Self::parity_aligned()
        }
    }

    /// Align genotyping semantics with [`crate::engine::CallRegionMode`].
    pub fn with_call_region_mode(mut self, mode: crate::engine::CallRegionMode) -> Self {
        self.semantics = GenotypingSemantics::from_call_region_mode(mode);
        self
    }
}
