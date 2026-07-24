//! Read×haplotype likelihood matrix cell (PairHMM / callRegion).
//! Leaf domain type shared by the engine and genotyping — kept out of [`crate::engine`]
//! so genotyping does not depend on the full callRegion engine module.

use crate::bio_ids::{HaplotypeIndex, ReadIndex};

/// One read×haplotype log10 likelihood from `callRegion`.
/// # Invariants
/// Indices refer to genotyping-read / haplotype lists on the region outcome.
/// `log10_likelihood` is a finite PairHMM (or engine) score in log10 space.
/// # Ownership
/// Owned scalar triple; matrices are flat vectors of these rows.
/// # Mutation
/// Immutable after PairHMM scoring.
/// # Biological assumptions
/// Likelihood of observing the read given the haplotype sequence.
/// # Java equivalence
/// GATK `ReadLikelihoods` matrix cell from `callRegion` PairHMM.
#[derive(Debug, Clone, PartialEq)]
pub struct RegionReadLikelihood {
    pub read_index: ReadIndex,
    pub haplotype_index: HaplotypeIndex,
    pub log10_likelihood: f64,
}
