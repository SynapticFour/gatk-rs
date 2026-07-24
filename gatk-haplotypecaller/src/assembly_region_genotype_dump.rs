//! Active region → assembly + genotyping slice (Java `HcFullParityGateDump` aligned).

use crate::allele_subsetting::subset_alleles_for_genotyping;
use crate::assembly_based_caller::{assemble_reads, AssembleReadsArgs};
use crate::assembly_region_assemble_dump::AssemblyRegionHaplotypeTarget;
use crate::parity_region_genotype::{
    parity_java_aligned_hap_log10_sums, parity_java_aligned_read_rows,
    parity_region_genotype_from_rows_with_gl_mode, write_parity_region_genotype_dump,
};
use crate::read_model::ReadFilterParams;
use crate::walker_apply::call_disposition;
use crate::walker_apply::AssemblyRegionCallDisposition;
use crate::walker_traversal::{
    flatten_assembly_regions, traverse_assembly_region_walker, WalkerTraversalConfig,
};
use gatk_common::{GatkError, GatkResult};
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use std::io::Write;
use std::path::Path;

fn write_region_header(
    out: &mut impl Write,
    region: &crate::assembly_region_iterator::AssemblyRegion,
) -> GatkResult<()> {
    writeln!(out, "region_contig\t{}", region.contig)?;
    writeln!(out, "region_start\t{}", region.start.get())?;
    writeln!(out, "region_end\t{}", region.end.get())?;
    writeln!(out, "is_active\t{}", region.is_active)?;
    Ok(())
}

fn write_java_genotype_dump(
    out: &mut impl Write,
    haplotypes: &[crate::haplotype::Haplotype],
    reads: &[rust_htslib::bam::Record],
    max_allele_count: Option<usize>,
) -> GatkResult<()> {
    let n_haps = haplotypes.len();
    if n_haps == 0 {
        writeln!(out, "haplotype_count\t0")?;
        writeln!(out, "read_count\t{}", reads.len())?;
        writeln!(out, "genotyped\tfalse")?;
        return Ok(());
    }
    let rows = parity_java_aligned_read_rows(reads, haplotypes)?;
    let is_ref: Vec<bool> = haplotypes.iter().map(|h| h.is_reference).collect();
    let legacy_gl = max_allele_count.is_none();
    let gt_dump = parity_region_genotype_from_rows_with_gl_mode(&rows, &is_ref, legacy_gl)?;
    write_parity_region_genotype_dump(out, &gt_dump)?;
    if let Some(max_allele_count) = max_allele_count {
        let agg = parity_java_aligned_hap_log10_sums(&rows)?;
        for (i, sum) in agg.haplotype_log10_sums.iter().enumerate() {
            writeln!(out, "haplotype_{i}_log10_sum\t{sum}")?;
        }
        for (i, h) in haplotypes.iter().enumerate() {
            writeln!(out, "haplotype_{i}_is_reference\t{}", h.is_reference)?;
            writeln!(out, "haplotype_{i}_bases\t{}", h.sequence_string())?;
        }
        writeln!(out, "max_allele_count\t{max_allele_count}")?;
        let kept = subset_alleles_for_genotyping(haplotypes, &agg, max_allele_count)?;
        writeln!(
            out,
            "kept_indices\t{}",
            kept.iter()
                .map(|i| i.to_string())
                .collect::<Vec<_>>()
                .join(",")
        )?;
        writeln!(
            out,
            "trim_triggered\t{}",
            haplotypes.len() > max_allele_count
        )?;
    }
    Ok(())
}

/// First active region: `HaplotypeCallerEngine::callRegion` genotype dump (G.2 / E2E genotyping).
pub fn dump_assembly_region_genotype_tsv(
    ref_fasta: &Path,
    bam_path: &Path,
    interval_cli: &str,
    padding: u64,
    target: AssemblyRegionHaplotypeTarget,
    out: &mut impl Write,
) -> GatkResult<()> {
    let dict = SequenceDictionary::from_fasta_path(ref_fasta)?;
    let specs = parse_intervals_cli_string(&dict, interval_cli)?;
    let filters = ReadFilterParams::gatk_standard_hc();
    let cfg = WalkerTraversalConfig::gatk_haplotype_caller_production(padding);
    let walk = traverse_assembly_region_walker(&dict, &specs, ref_fasta, bam_path, &filters, &cfg)?;
    let regions = flatten_assembly_regions(&walk);

    let region = match target {
        AssemblyRegionHaplotypeTarget::Active => regions
            .iter()
            .find(|r| {
                matches!(
                    call_disposition(r),
                    AssemblyRegionCallDisposition::ActiveFull
                )
            })
            .ok_or_else(|| GatkError::argument("no active assembly region in interval"))?,
        AssemblyRegionHaplotypeTarget::Inactive => regions
            .iter()
            .find(|r| !r.is_active)
            .ok_or_else(|| GatkError::argument("no inactive assembly region in interval"))?,
    };

    write_region_header(out, region)?;
    if target == AssemblyRegionHaplotypeTarget::Inactive || !region.is_active {
        writeln!(out, "genotyped\tfalse")?;
        writeln!(out, "haplotype_count\t0")?;
        writeln!(out, "read_count\t0")?;
        return Ok(());
    }

    // Java `assemblyRegionGenotype` uses `assembleReads` + `HcParityRegionGenotype`, not full `callRegion`.
    let assemble_args = AssembleReadsArgs::java_parity_gate_default("-");
    let mut ref_cache = gatk_core::reference::ReferenceWindowCache::new(ref_fasta.to_path_buf(), 4);
    let assembly = assemble_reads(region, &dict, &mut ref_cache, &assemble_args)?;
    write_java_genotype_dump(out, &assembly.haplotypes, &region.reads, None)
}

/// G-D05 live path: matches Java `assemblyRegionGenotypeSubset` (`assembleReads` + parity PairHMM + subset).
pub fn dump_assembly_region_genotype_subset_tsv(
    ref_fasta: &Path,
    bam_path: &Path,
    interval_cli: &str,
    padding: u64,
    target: AssemblyRegionHaplotypeTarget,
    max_allele_count: usize,
    assembly_profile: &str,
    out: &mut impl Write,
) -> GatkResult<()> {
    let dict = SequenceDictionary::from_fasta_path(ref_fasta)?;
    let specs = parse_intervals_cli_string(&dict, interval_cli)?;
    let filters = ReadFilterParams::gatk_standard_hc();
    let cfg = WalkerTraversalConfig::gatk_haplotype_caller_production(padding);
    let walk = traverse_assembly_region_walker(&dict, &specs, ref_fasta, bam_path, &filters, &cfg)?;
    let regions = flatten_assembly_regions(&walk);

    let region = match target {
        AssemblyRegionHaplotypeTarget::Active => regions
            .iter()
            .find(|r| {
                matches!(
                    call_disposition(r),
                    AssemblyRegionCallDisposition::ActiveFull
                )
            })
            .ok_or_else(|| GatkError::argument("no active assembly region in interval"))?,
        AssemblyRegionHaplotypeTarget::Inactive => regions
            .iter()
            .find(|r| !r.is_active)
            .ok_or_else(|| GatkError::argument("no inactive assembly region in interval"))?,
    };

    write_region_header(out, region)?;
    if target == AssemblyRegionHaplotypeTarget::Inactive || !region.is_active {
        writeln!(out, "genotyped\tfalse")?;
        writeln!(out, "haplotype_count\t0")?;
        writeln!(out, "read_count\t0")?;
        return Ok(());
    }

    let profile_for_call = if assembly_profile == "-" {
        "default"
    } else {
        assembly_profile
    };
    let profile_label = if assembly_profile == "default" {
        "-"
    } else {
        assembly_profile
    };
    writeln!(out, "assembly_profile\t{profile_label}")?;

    let assemble_args = AssembleReadsArgs::java_parity_gate_default(profile_for_call);
    let mut ref_cache = gatk_core::reference::ReferenceWindowCache::new(ref_fasta.to_path_buf(), 4);
    let assembly = assemble_reads(region, &dict, &mut ref_cache, &assemble_args)?;
    write_java_genotype_dump(
        out,
        &assembly.haplotypes,
        &region.reads,
        Some(max_allele_count),
    )
}
