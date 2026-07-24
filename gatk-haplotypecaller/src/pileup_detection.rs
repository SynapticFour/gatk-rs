//! GATK `PileupDetectionArgumentCollection` — HC default path + parity supplement.
//! Standard Illumina `HaplotypeCaller` runs with `usePileupDetection=false`. When false, Java does
//! not call `processPileupAlleles`; assembly variation comes from haplotype CIGARs only.
//! Rust parity bridge: [`supplement_assembly_pileup_events_from_reads`] (strict read pileup) after
//! trim when [`PileupDetectionConfig::enable_event_supplement`] is set.

use crate::alignment::SwParameters;
use crate::assembly_result_set::AssemblyResultSet;
use crate::read_event_discovery::supplement_assembly_pileup_events_from_reads;
use gatk_common::GatkResult;
use rust_htslib::bam::Record;

/// GATK `PileupDetectionArgumentCollection.usePileupDetection` default for HC.
pub const GATK_HC_USE_PILEUP_DETECTION_DEFAULT: bool = false;

/// Whether pileup-based forced alleles are enabled (parity: must match Java for the run).
/// # Invariants
/// Default HC has `use_pileup_detection == false` (GATK 4.4); supplement is off unless explicitly enabled.
/// `enable_event_supplement` is a Rust bridge only; must not be set when strict Java parity is required.
/// # Ownership
/// [`Copy`] config threaded through assembly supplement hooks.
/// # Mutation
/// Immutable per region; assembly result sets are mutated by supplement functions when enabled.
/// # Biological assumptions
/// Standard HC variation comes from haplotype CIGARs; pileup supplement fills ref-only assembly gaps.
/// # Java equivalence
/// GATK `PileupDetectionArgumentCollection.usePileupDetection`; `enable_event_supplement` is Rust-native.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PileupDetectionConfig {
    pub use_pileup_detection: bool,
    /// Rust bridge: strict read-pileup events + haps when assembly CIGARs miss variants (A1).
    pub enable_event_supplement: bool,
}

impl Default for PileupDetectionConfig {
    fn default() -> Self {
        Self::gatk_haplotype_caller_defaults()
    }
}

impl PileupDetectionConfig {
    /// GATK 4.4 default HC — pileup detection off; no Rust event supplement (Java parity).
    pub fn gatk_haplotype_caller_defaults() -> Self {
        Self {
            use_pileup_detection: GATK_HC_USE_PILEUP_DETECTION_DEFAULT,
            enable_event_supplement: false,
        }
    }
}

/// Run pileup supplement only when assembly CIGARs produced no variation events (P0 rust-only).
pub fn should_run_pileup_supplement(
    assembly: &AssemblyResultSet,
    apply_bases: &[u8],
    apply_pad: u64,
) -> bool {
    let events = crate::event_map::collect_variation_events(
        &assembly.haplotypes,
        apply_bases,
        apply_pad,
        &assembly.contig,
        assembly.max_mnp_distance(),
    );
    events.is_empty()
}

/// A1 pileup supplement (see [`supplement_assembly_pileup_events_from_reads`]).
pub fn supplement_pileup_events_into_assembly(
    assembly: &mut AssemblyResultSet,
    reads: &[Record],
    active_start_1based: u64,
    active_end_1based: u64,
    sw: &SwParameters,
) -> GatkResult<()> {
    supplement_assembly_pileup_events_from_reads(
        assembly,
        reads,
        active_start_1based,
        active_end_1based,
        sw,
    )
}
