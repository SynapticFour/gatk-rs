//! E2E integration: active region → assembly → native `Log10PairHMM` read×haplotype matrix.

use crate::assembly_based_caller::{assemble_reads, AssembleReadsArgs};
use crate::assembly_region_assemble_dump::AssemblyRegionHaplotypeTarget;
use crate::assembly_region_finalize::{
    finalize_region_reads_for_assembly, gatk_min_tail_quality_for_assembly,
};
use crate::assembly_region_iterator::AssemblyRegion;
use crate::pairhmm_log10::log10_pairhmm_likelihood_parity_defaults;
use crate::pairhmm_qual::cap_read_base_qualities;
use crate::read_model::ReadFilterParams;
use crate::walker_apply::call_disposition;
use crate::walker_apply::AssemblyRegionCallDisposition;
use crate::walker_traversal::{
    flatten_assembly_regions, traverse_assembly_region_walker, WalkerTraversalConfig,
};
use gatk_common::{GatkError, GatkResult};
use gatk_core::reference::{parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary};
use std::io::Write;
use std::path::Path;

fn format_ll(v: f64) -> String {
    if v.is_infinite() && v.is_sign_negative() {
        "-inf".to_string()
    } else {
        format!("{v:.17}")
    }
}

fn write_region_header(out: &mut impl Write, region: &AssemblyRegion) -> GatkResult<()> {
    writeln!(out, "region_contig\t{}", region.contig)?;
    writeln!(out, "region_start\t{}", region.start.get())?;
    writeln!(out, "region_end\t{}", region.end.get())?;
    writeln!(out, "is_active\t{}", region.is_active)?;
    Ok(())
}

/// First active region: assemble haplotypes, score each read×haplotype with GATK `Log10PairHMM` defaults.
pub fn dump_assembly_region_pairhmm_likelihoods_tsv(
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
        writeln!(out, "read_count\t0")?;
        writeln!(out, "haplotype_count\t0")?;
        return Ok(());
    }

    let args = AssembleReadsArgs::default();
    let finalized = finalize_region_reads_for_assembly(
        &region.reads,
        region,
        args.correct_overlapping_base_qualities,
        gatk_min_tail_quality_for_assembly(args.assembler.min_base_quality),
        false,
    );
    let ars = assemble_reads(region, &dict, &mut ref_cache, &args)?;
    let haplotypes: Vec<&[u8]> = ars.haplotypes.iter().map(|h| h.bases.as_slice()).collect();

    writeln!(out, "read_count\t{}", finalized.len())?;
    writeln!(out, "haplotype_count\t{}", haplotypes.len())?;
    writeln!(out, "read_index\thaplotype_index\tlog10_likelihood")?;

    const BQ_THRESHOLD: u8 = 18;
    let mut read_records: Vec<_> = finalized.into_iter().collect();
    read_records.sort_by(|a, b| a.qname().cmp(b.qname()).then_with(|| a.pos().cmp(&b.pos())));

    for (ri, rec) in read_records.iter().enumerate() {
        let bases = rec.seq().as_bytes();
        let mut quals = rec.qual().to_vec();
        let mapq = rec.mapq();
        cap_read_base_qualities(&mut quals, mapq, BQ_THRESHOLD, false);
        for (hi, hap) in haplotypes.iter().enumerate() {
            let ll = log10_pairhmm_likelihood_parity_defaults(&bases, &quals, hap)?;
            writeln!(out, "{ri}\t{hi}\t{}", format_ll(ll))?;
        }
    }
    Ok(())
}
