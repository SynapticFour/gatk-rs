//! Active assembly region → `ReadThreadingAssembler` haplotypes (GATK `callRegion` assembly slice).

use crate::assembly_based_caller::AssembleReadsArgs;
use crate::assembly_graph_dump::write_haplotype_rows_with_ref_recovery;
use crate::assembly_region_finalize::{
    assembly_reads_for_java_materialize_dump, assembly_reference_read,
    reference_haplotype_for_assembly_region,
};
use crate::assembly_region_iterator::AssemblyRegion;
use crate::read_model::ReadFilterParams;
use crate::read_threading_assembler::probe_seq_graph_kmer_attempts;
use crate::walker_apply::{
    call_disposition, select_region_for_asm_dump, AssemblyRegionCallDisposition,
};
use crate::walker_traversal::{
    flatten_assembly_regions, traverse_assembly_region_walker, WalkerTraversalConfig,
};
use gatk_common::{GatkError, GatkResult};
use gatk_core::reference::ReferenceWindowCache;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use std::io::Write;
use std::path::Path;

/// Which assembly region to dump (first match in interval stream).
/// # Invariants
/// `Active` selects first `ActiveFull` disposition; `Inactive` dumps metadata only (no assembly).
/// # Ownership
/// [`Copy`] CLI target enum.
/// # Mutation
/// Immutable dump selector.
/// # Biological assumptions
/// None — parity dump routing.
/// # Java equivalence
/// Rust-native E2E dump selector over GATK `callRegion` active vs inactive branches.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AssemblyRegionHaplotypeTarget {
    /// First region with `ActiveFull` disposition (GATK `assembleReads` path).
    Active,
    /// First inactive region — metadata only, no assembly (`E2E.5`).
    Inactive,
}

impl AssemblyRegionHaplotypeTarget {
    pub fn parse_cli(s: &str) -> Option<Self> {
        match s {
            "active" => Some(Self::Active),
            "inactive" => Some(Self::Inactive),
            _ => None,
        }
    }
}

fn write_region_meta(
    out: &mut impl Write,
    region: &AssemblyRegion,
    status: &str,
    kmer_size: usize,
) -> GatkResult<()> {
    writeln!(out, "region_contig\t{}", region.contig)
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    writeln!(out, "region_start\t{}", region.start.get())
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    writeln!(out, "region_end\t{}", region.end.get())
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    writeln!(out, "is_active\t{}", region.is_active)
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    writeln!(out, "status\t{status}").map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    writeln!(out, "kmer_size\t{kmer_size}")
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    Ok(())
}

/// First matching region in interval stream → haplotype TSV (E2E / `callRegion` assembly slice).
pub fn dump_assembly_region_haplotypes_tsv(
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
    let mut ref_cache = ReferenceWindowCache::new(ref_fasta.to_path_buf(), 4);

    let region = match target {
        AssemblyRegionHaplotypeTarget::Active => select_region_for_asm_dump(&regions)
            .ok_or_else(|| GatkError::argument("no active assembly region in interval"))?,
        AssemblyRegionHaplotypeTarget::Inactive => regions
            .iter()
            .find(|r| !r.is_active)
            .ok_or_else(|| GatkError::argument("no inactive assembly region in interval"))?,
    };

    if target == AssemblyRegionHaplotypeTarget::Inactive || !region.is_active {
        write_region_meta(out, region, "inactive_skip", 0)?;
        return Ok(());
    }

    let reference = assembly_reference_read(&dict, &mut ref_cache, region)?;
    let assemble_args = AssembleReadsArgs::default();
    let assembly_set = crate::assembly_based_caller::assemble_reads(
        region,
        &dict,
        &mut ref_cache,
        &assemble_args,
    )?;
    // Java `emitAssemblyRegionHaplotypesFromMaterial`: `isVariationPresent` is still false
    // after `runLocalAssembly` (SeqGraph `findBestPaths` uses `add(h)` only — scores stay 0).
    // Kmer metadata is 0. Match frozen java_dumps/e2e/*.tsv (scores always 0).
    write_region_meta(out, region, "just_assembled_reference", 0)?;
    let ref_bases = reference.bases.as_slice();
    let mut haplotypes = if assembly_set.haplotypes.len() == 1
        && assembly_set.haplotypes[0].is_reference
        && assembly_set.haplotypes[0].bases.as_slice() != ref_bases
    {
        // Cyclic k=85 graphs can mark a stitched ref path `is_reference` with non-ref-length bases.
        vec![reference_haplotype_for_assembly_region(
            &reference, region, &dict,
        )]
    } else {
        assembly_set.haplotypes
    };
    for h in &mut haplotypes {
        h.score = 0.0;
    }
    write_haplotype_rows_with_ref_recovery(&haplotypes, ref_bases, &reference.bases, out)
}

/// Per-kmer SeqGraph probe for the first active region (assembly divergence diagnostics).
pub fn dump_assembly_region_kmer_probe_tsv(
    ref_fasta: &Path,
    bam_path: &Path,
    interval_cli: &str,
    padding: u64,
    out: &mut impl Write,
) -> GatkResult<()> {
    let dict = SequenceDictionary::from_fasta_path(ref_fasta)?;
    let specs = parse_intervals_cli_string(&dict, interval_cli)?;
    let filters = ReadFilterParams::gatk_standard_hc();
    let cfg = WalkerTraversalConfig::gatk_haplotype_caller_production(padding);
    let walk = traverse_assembly_region_walker(&dict, &specs, ref_fasta, bam_path, &filters, &cfg)?;
    let regions = flatten_assembly_regions(&walk);
    let mut ref_cache = ReferenceWindowCache::new(ref_fasta.to_path_buf(), 4);

    let region = select_region_for_asm_dump(&regions)
        .ok_or_else(|| GatkError::argument("no assembly region with reads in interval"))?;
    if !matches!(
        call_disposition(region),
        AssemblyRegionCallDisposition::ActiveFull
    ) {
        writeln!(out, "warn\tregion_inactive_kmer_probe_on_inactive_reads")
            .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    }

    let reference = assembly_reference_read(&dict, &mut ref_cache, region)?;
    let assemble_args = AssembleReadsArgs::default();
    let reads = assembly_reads_for_java_materialize_dump(&region.reads);

    writeln!(out, "region_contig\t{}", region.contig)
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    writeln!(out, "region_start\t{}", region.start.get())
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    writeln!(out, "region_end\t{}", region.end.get())
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    writeln!(out, "padded_ref_len\t{}", reference.bases.len())
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    writeln!(out, "read_count\t{}", reads.len())
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;

    let rows = probe_seq_graph_kmer_attempts(&reference, &reads, &assemble_args.assembler)?;
    writeln!(
        out,
        "phase\tkmer\tallow_low_complexity\tallow_non_unique\toutcome\tthread_nodes\tthread_edges\tcleanup_status\thas_ref_source\thas_ref_sink\tref_path_matches\tkbest_paths\textracted_haps\tnon_ref_haps\tpath_bases_len\tpath_eq_ref"
    )
    .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    for r in &rows {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.phase,
            r.kmer_size,
            r.allow_low_complexity,
            r.allow_non_unique_ref,
            r.outcome,
            r.thread_nodes,
            r.thread_edges,
            r.cleanup_status,
            r.has_ref_source,
            r.has_ref_sink,
            r.ref_path_matches,
            r.kbest_paths,
            r.extracted_haps,
            r.non_ref_haps,
            r.path_bases_len,
            r.path_eq_ref_bases,
        )
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    }
    Ok(())
}
