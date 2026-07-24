//! Sample ploidy vs activity-evaluation ploidy (L2 parity `c5-ploidy`).

use gatk_common::{GatkError, GatkResult};
use std::io::Write;

/// GATK `HaplotypeCallerEngine.MINIMUM_PUTATIVE_PLOIDY_FOR_ACTIVE_REGION_DISCOVERY`.
pub const GATK_MINIMUM_PUTATIVE_PLOIDY_FOR_ACTIVE_REGION_DISCOVERY: u32 = 2;

/// Ploidy used in `isActive` (constant engine ploidy, floored at minimum).
pub fn activity_evaluation_ploidy(sample_ploidy: u32) -> u32 {
    sample_ploidy.max(GATK_MINIMUM_PUTATIVE_PLOIDY_FOR_ACTIVE_REGION_DISCOVERY)
}

/// Ploidy for genotyping (`HomogeneousPloidyModel` / `getPloidyToUseAtThisSite` on uniform samples).
pub fn genotyping_ploidy(sample_ploidy: u32) -> u32 {
    sample_ploidy
}

pub fn dump_ploidy_resolution_tsv(sample_ploidy: u32, out: &mut impl Write) -> GatkResult<()> {
    writeln!(out, "sample_ploidy\t{sample_ploidy}")
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    writeln!(
        out,
        "activity_eval_ploidy\t{}",
        activity_evaluation_ploidy(sample_ploidy)
    )
    .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    writeln!(
        out,
        "genotyping_ploidy\t{}",
        genotyping_ploidy(sample_ploidy)
    )
    .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    Ok(())
}
