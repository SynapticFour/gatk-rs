//! `AssemblyRegionWalker.apply` vs `HaplotypeCallerEngine.callRegion` disposition.
//! Java: [`AssemblyRegionWalker`](https://github.com/broadinstitute/gatk/blob/master/src/main/java/org/broadinstitute/hellbender/engine/AssemblyRegionWalker.java)
//! invokes `apply` once per region from [`AssemblyRegionIterator`](https://github.com/broadinstitute/gatk/blob/master/src/main/java/org/broadinstitute/hellbender/engine/AssemblyRegionIterator.java).
//! [`HaplotypeCallerEngine.callRegion`](https://github.com/broadinstitute/gatk/blob/master/src/main/java/org/broadinstitute/hellbender/tools/walkers/haplotypecaller/HaplotypeCallerEngine.java)
//! takes the **inactive fast path** (no graph / local assembly) when `!region.isActive`, and otherwise runs assembly + genotyping.

use crate::assembly_region_iterator::AssemblyRegion;

/// How HC-style `callRegion` would treat this region (assembly vs reference-only fast path).
/// # Invariants
/// Exactly one of active full vs inactive fast path per region (`region.is_active` bit).
/// Inactive path skips `assembleReads` in Java parity semantics.
/// # Ownership
/// [`Copy`] enum derived from [`AssemblyRegion`] by [`call_disposition`].
/// # Mutation
/// Immutable tag computed per region at apply/call time.
/// # Biological assumptions
/// Inactive regions lack sufficient activity evidence for local re-assembly.
/// # Java equivalence
/// GATK `HaplotypeCallerEngine.callRegion` inactive guard vs full assembly path.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AssemblyRegionCallDisposition {
    /// `region.isActive` — full assembly + genotyping path (GATK `callRegion` after the inactive guard).
    ActiveFull,
    /// `!region.isActive` — `referenceModelForNoVariation` path; **no** `assembleReads` (GATK fast path).
    InactiveReferenceFastPath,
}

#[inline]
pub fn call_disposition(region: &AssemblyRegion) -> AssemblyRegionCallDisposition {
    if region.is_active {
        AssemblyRegionCallDisposition::ActiveFull
    } else {
        AssemblyRegionCallDisposition::InactiveReferenceFastPath
    }
}

/// First region suitable for ASM-1 parity dumps: active if present, else inactive with reads.
#[cfg(any(feature = "dev-dumps", test))]
pub fn select_region_for_asm_dump(regions: &[AssemblyRegion]) -> Option<&AssemblyRegion> {
    regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            )
        })
        .or_else(|| {
            regions.iter().find(|r| {
                !r.reads.is_empty()
                    && matches!(
                        call_disposition(r),
                        AssemblyRegionCallDisposition::InactiveReferenceFastPath
                    )
            })
        })
}

/// Counts aligned with Java “one `apply` per iterator region” + inactive fast path in `callRegion`.
/// # Invariants
/// `total_apply == inactive_fast_path + active_full` (debug_assert in [`Self::from_regions`]).
/// One count increment per [`AssemblyRegion`] in the iterator stream.
/// # Ownership
/// [`Copy`] aggregate stats; regions are not retained.
/// # Mutation
/// Immutable snapshot after traversal; built by scanning region dispositions.
/// # Biological assumptions
/// None documented (walker bookkeeping for parity dumps).
/// # Java equivalence
/// GATK `AssemblyRegionWalker.apply` per-region invocations + `callRegion` inactive fast-path tally.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct WalkerApplyStats {
    /// Same as number of `AssemblyRegionWalker.apply` invocations for this shard/interval stream.
    pub total_apply: usize,
    pub inactive_fast_path: usize,
    pub active_full: usize,
}

impl WalkerApplyStats {
    pub fn from_regions(regions: &[AssemblyRegion]) -> Self {
        let mut s = WalkerApplyStats::default();
        s.total_apply = regions.len();
        for r in regions {
            match call_disposition(r) {
                AssemblyRegionCallDisposition::InactiveReferenceFastPath => {
                    s.inactive_fast_path += 1;
                }
                AssemblyRegionCallDisposition::ActiveFull => {
                    s.active_full += 1;
                }
            }
        }
        debug_assert_eq!(s.total_apply, s.inactive_fast_path + s.active_full);
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genome_loc::GenomePosition;

    fn region(active: bool) -> AssemblyRegion {
        AssemblyRegion {
            contig: "chr1".into(),
            start: GenomePosition::new_1based(1),
            end: GenomePosition::new_1based(10),
            is_active: active,
            extended_start: GenomePosition::new_1based(1),
            extended_end: GenomePosition::new_1based(10),
            extension: 0,
            reads: Vec::new(),
            read_qnames: Vec::new(),
            reference: crate::reference_context::ReferenceContext::empty(),
            features: crate::feature_context::FeatureContext::empty(),
            pileup_loci: Vec::new(),
        }
    }

    #[test]
    fn disposition_matches_is_active() {
        assert_eq!(
            call_disposition(&region(false)),
            AssemblyRegionCallDisposition::InactiveReferenceFastPath
        );
        assert_eq!(
            call_disposition(&region(true)),
            AssemblyRegionCallDisposition::ActiveFull
        );
    }

    #[test]
    fn stats_sum_invariant() {
        let v = vec![region(false), region(true), region(false)];
        let st = WalkerApplyStats::from_regions(&v);
        assert_eq!(st.total_apply, 3);
        assert_eq!(st.inactive_fast_path, 2);
        assert_eq!(st.active_full, 1);
    }
}
