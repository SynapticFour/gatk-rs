//! PairHMM region scoring used by `call_region` (split from `engine.rs` for N-3).
use super::*;
use crate::assembly_region_finalize::{
    clip_finalized_reads_in_place, finalize_region_reads_for_assembly,
    gatk_min_tail_quality_for_assembly,
};

pub(super) fn compute_region_read_likelihoods(
    region: &AssemblyRegion,
    haplotypes: &[Haplotype],
    config: &HcLikelihoodEngineConfig,
    apply_normalize: bool,
    pre_finalized: Option<Vec<rust_htslib::bam::Record>>,
) -> GatkResult<Vec<RegionReadLikelihood>> {
    if haplotypes.is_empty() {
        return Ok(Vec::new());
    }
    // A2: consume assemble finalize buffer when present (clip in place — no second owned copy).
    let finalized = if let Some(mut pre) = pre_finalized.filter(|p| !p.is_empty()) {
        clip_finalized_reads_in_place(&mut pre, region);
        pre
    } else {
        finalize_region_reads_for_assembly(
            &region.reads,
            region,
            true,
            gatk_min_tail_quality_for_assembly(10),
            false,
        )
    };
    let active_span = Some((region.start.get(), region.end.get()));
    // Trim/hard-clip can drop sparse-BAM reads that still overlap the active locus.
    if finalized.is_empty() && !region.reads.is_empty() {
        let out = score_pairhmm_from_records(region.reads.as_slice(), haplotypes, config)?;
        return Ok(post_process_pairhmm_likelihoods(
            out,
            region.reads.as_slice(),
            haplotypes,
            apply_normalize,
            active_span,
        ));
    }
    let out = score_pairhmm_from_records(&finalized, haplotypes, config)?;
    Ok(post_process_pairhmm_likelihoods(
        out,
        &finalized,
        haplotypes,
        apply_normalize,
        active_span,
    ))
}
