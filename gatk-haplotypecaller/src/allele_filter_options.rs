//! Typed options for [`crate::allele_filtering::filter_assembly_and_likelihoods`] (Sprint K-1 / R2).

use crate::genome_loc::GenomePosition;

/// Active-region span for allele filtering (both ends or neither).
/// # Compiler-enforced
/// One-sided spans are impossible: either unrestricted or inclusive `[start, end]`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ActiveRegionSpan {
    Unrestricted,
    Inclusive {
        start: GenomePosition,
        end: GenomePosition,
    },
}

/// Options for haplotype allele filtering after PairHMM.
/// # Invariants
/// [`ActiveRegionSpan`] pairs start/end when filtering is span-scoped.
/// `strict_java_snp_rank_only` enables Java SNP-only ranking rules when true.
/// # Ownership
/// [`Copy`] options; coordinates use [`GenomePosition`] newtypes.
/// # Mutation
/// Immutable per filter invocation.
/// # Biological assumptions
/// Filters reduce haplotype alleles to those supported by read likelihoods within the active region.
/// # Java equivalence
/// GATK `AlleleFilteringHC.filterAlleles` options (Sprint K-1 / R2 rust-native wrapper).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AlleleFilterOptions {
    /// When true, rank/filter using Java SNP-only rules (`strict_java` path).
    pub strict_java_snp_rank_only: bool,
    pub span: ActiveRegionSpan,
}

impl AlleleFilterOptions {
    pub fn strict_java_span(active_start_1based: u64, active_end_1based: u64) -> Self {
        Self {
            strict_java_snp_rank_only: true,
            span: ActiveRegionSpan::Inclusive {
                start: GenomePosition::new_1based(active_start_1based),
                end: GenomePosition::new_1based(active_end_1based),
            },
        }
    }

    pub fn from_strict_java(
        is_strict_java: bool,
        active_start_1based: Option<u64>,
        active_end_1based: Option<u64>,
    ) -> Self {
        let span = match (active_start_1based, active_end_1based) {
            (Some(s), Some(e)) => ActiveRegionSpan::Inclusive {
                start: GenomePosition::new_1based(s),
                end: GenomePosition::new_1based(e),
            },
            // Half-spans previously ignored the filter (unrestricted); preserve that.
            _ => ActiveRegionSpan::Unrestricted,
        };
        Self {
            strict_java_snp_rank_only: is_strict_java,
            span,
        }
    }

    /// Unscoped / non-strict (tests and legacy).
    pub fn unrestricted() -> Self {
        Self {
            strict_java_snp_rank_only: false,
            span: ActiveRegionSpan::Unrestricted,
        }
    }

    #[inline]
    pub fn active_start_1based(self) -> Option<u64> {
        match self.span {
            ActiveRegionSpan::Inclusive { start, .. } => Some(start.get()),
            ActiveRegionSpan::Unrestricted => None,
        }
    }

    #[inline]
    pub fn active_end_1based(self) -> Option<u64> {
        match self.span {
            ActiveRegionSpan::Inclusive { end, .. } => Some(end.get()),
            ActiveRegionSpan::Unrestricted => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn half_span_options_become_unrestricted() {
        let o = AlleleFilterOptions::from_strict_java(true, Some(10), None);
        assert_eq!(o.span, ActiveRegionSpan::Unrestricted);
        assert_eq!(o.active_start_1based(), None);
        assert_eq!(o.active_end_1based(), None);
    }

    #[test]
    fn both_ends_form_inclusive_span() {
        let o = AlleleFilterOptions::strict_java_span(10, 20);
        assert_eq!(o.active_start_1based(), Some(10));
        assert_eq!(o.active_end_1based(), Some(20));
        assert!(matches!(o.span, ActiveRegionSpan::Inclusive { .. }));
    }
}
