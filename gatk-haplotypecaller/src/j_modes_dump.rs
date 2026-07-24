//! HC output mode emission decisions.

use crate::genotyping::{decide_locus_emission, summarize_no_variation_region, EmitMode};
use gatk_common::GatkResult;
use std::io::Write;

pub fn dump_emit_mode_decision_tsv(
    mode: EmitMode,
    has_variant: bool,
    locus_count: usize,
    out: &mut impl Write,
) -> GatkResult<()> {
    let decision = decide_locus_emission(mode, has_variant);
    let summary = summarize_no_variation_region(mode, locus_count);
    writeln!(
        out,
        "emit_mode\t{}",
        match mode {
            EmitMode::Vcf => "VCF",
            EmitMode::Gvcf => "GVCF",
            EmitMode::BpResolution => "BP_RESOLUTION",
        }
    )?;
    writeln!(out, "has_variant\t{has_variant}")?;
    writeln!(out, "locus_decision\t{decision:?}")?;
    writeln!(
        out,
        "no_var_emit_blocks\t{}",
        summary.reference_blocks_emitted
    )?;
    writeln!(
        out,
        "no_var_emit_sites\t{}",
        summary.reference_sites_emitted
    )?;
    Ok(())
}

/// J-D05 — DRAGEN output-mode branch scaffold (off on default HC).
pub fn dump_dragen_mode_branch_tsv(out: &mut impl Write) -> GatkResult<()> {
    writeln!(out, "dragen_mode_active\tfalse")?;
    writeln!(out, "emit_mode_default\tVCF")?;
    writeln!(out, "read_shard_dragen_pipeline\tfalse")?;
    Ok(())
}
