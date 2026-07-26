//! ASM-1: per-read `finalizeRegion` parity (`assembly-region-finalize-reads`).

use crate::assembly_region_finalize::{
    assembly_reads_for_production, gatk_min_tail_quality_for_assembly,
};
use crate::assembly_region_iterator::AssemblyRegion;
use crate::read_downsample::compare_read_coordinates_java;
use crate::read_model::ReadFilterParams;
use crate::read_pre_len::unclipped_read_length;
use crate::read_unclip::alignment_end_1based;
use crate::walker::GATK_DEFAULT_ASSEMBLY_REGION_PADDING;
use crate::walker_apply::{
    call_disposition, select_region_for_asm_dump, AssemblyRegionCallDisposition,
};
use crate::walker_traversal::{
    flatten_assembly_regions, traverse_assembly_region_walker, WalkerTraversalConfig,
};
use gatk_common::{GatkError, GatkResult};
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use rust_htslib::bam;
use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::Read as _;
use std::io::Write;
use std::path::Path;

#[derive(Debug, Clone)]
struct FinalizeReadRow {
    qname: String,
    phase: &'static str,
    flags: u16,
    mapq: u8,
    start: i64,
    end: i64,
    cigar: String,
    seq_len: usize,
    unclipped_len: usize,
    fragment_length: i64,
    unmapped: bool,
}

fn format_cigar(rec: &bam::Record) -> String {
    rec.cigar()
        .iter()
        .map(|c| {
            let op = match c {
                Cigar::Match(_) => 'M',
                Cigar::Ins(_) => 'I',
                Cigar::Del(_) => 'D',
                Cigar::SoftClip(_) => 'S',
                Cigar::HardClip(_) => 'H',
                Cigar::Equal(_) => '=',
                Cigar::Diff(_) => 'X',
                Cigar::RefSkip(_) => 'N',
                Cigar::Pad(_) => 'P',
            };
            format!("{}{}", crate::read_unclip::cigar_len(c), op)
        })
        .collect()
}

fn finalize_read_row(rec: &bam::Record, phase: &'static str) -> FinalizeReadRow {
    let unmapped = rec.tid() < 0 || rec.is_unmapped() || rec.seq().is_empty();
    let start = if unmapped { 0 } else { rec.pos() + 1 };
    let end = if unmapped {
        0
    } else {
        i64::from(alignment_end_1based(rec))
    };
    FinalizeReadRow {
        qname: String::from_utf8_lossy(rec.qname()).into_owned(),
        phase,
        flags: rec.flags(),
        mapq: rec.mapq(),
        start,
        end,
        cigar: if unmapped {
            String::new()
        } else {
            format_cigar(rec)
        },
        seq_len: rec.seq().as_bytes().len(),
        unclipped_len: unclipped_read_length(rec),
        fragment_length: rec.insert_size(),
        unmapped,
    }
}

fn select_first_active_region(regions: &[AssemblyRegion]) -> GatkResult<&AssemblyRegion> {
    select_region_for_asm_dump(regions)
        .ok_or_else(|| GatkError::argument("no assembly region with reads in interval"))
}

/// First active region: raw iterator reads vs production `finalizeRegion` rows.
pub fn dump_assembly_region_finalize_reads_tsv(
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
    let region = select_first_active_region(&regions)?;

    let min_base_quality = 10u8;
    let min_tail = gatk_min_tail_quality_for_assembly(min_base_quality);
    let finalized = assembly_reads_for_production(&region.reads, region, min_tail, true, false);

    let header = bam::Reader::from_path(bam_path)
        .map_err(|e| GatkError::generic(format!("open bam: {e}")))?
        .header()
        .clone();
    let sort_java = |reads: &[bam::Record], phase: &'static str| -> Vec<FinalizeReadRow> {
        let mut sorted: Vec<&bam::Record> = reads.iter().collect();
        sorted.sort_by(|a, b| {
            compare_read_coordinates_java(a, b, &header)
                .cmp(&0)
                .then_with(|| a.qname().cmp(b.qname()))
        });
        sorted
            .into_iter()
            .map(|r| finalize_read_row(r, phase))
            .collect()
    };
    let raw_rows = sort_java(&region.reads, "raw");
    let fin_rows = sort_java(&finalized, "finalize");

    writeln!(out, "region_contig\t{}", region.contig)
        .map_err(|e| GatkError::generic(format!("write: {e}")))?;
    writeln!(out, "region_start\t{}", region.start.get())
        .map_err(|e| GatkError::generic(format!("write: {e}")))?;
    writeln!(out, "region_end\t{}", region.end.get())
        .map_err(|e| GatkError::generic(format!("write: {e}")))?;
    writeln!(out, "extended_start\t{}", region.extended_start.get())
        .map_err(|e| GatkError::generic(format!("write: {e}")))?;
    writeln!(out, "extended_end\t{}", region.extended_end.get())
        .map_err(|e| GatkError::generic(format!("write: {e}")))?;
    writeln!(out, "read_path\tfinalize").map_err(|e| GatkError::generic(format!("write: {e}")))?;
    writeln!(out, "raw_read_count\t{}", raw_rows.len())
        .map_err(|e| GatkError::generic(format!("write: {e}")))?;
    writeln!(out, "finalize_read_count\t{}", fin_rows.len())
        .map_err(|e| GatkError::generic(format!("write: {e}")))?;
    if !matches!(
        call_disposition(region),
        AssemblyRegionCallDisposition::ActiveFull
    ) {
        writeln!(out, "warn\tregion_inactive_using_reads_for_asm_diagnostic")
            .map_err(|e| GatkError::generic(format!("write: {e}")))?;
    }
    writeln!(
        out,
        "read\tqname\tphase\tflags\tmapq\tstart\tend\tcigar\tseq_len\tunclipped_len\tfragment_length\tunmapped"
    )
    .map_err(|e| GatkError::generic(format!("write: {e}")))?;
    for row in raw_rows.iter().chain(fin_rows.iter()) {
        writeln!(
            out,
            "read\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            row.qname,
            row.phase,
            row.flags,
            row.mapq,
            row.start,
            row.end,
            row.cigar,
            row.seq_len,
            row.unclipped_len,
            row.fragment_length,
            row.unmapped
        )
        .map_err(|e| GatkError::generic(format!("write: {e}")))?;
    }
    Ok(())
}

/// Explicit production finalize path label for assembly-stages (same graph as default finalize).
pub fn dump_assembly_region_assembly_stages_finalize_tsv(
    ref_fasta: &Path,
    bam_path: &Path,
    interval_cli: &str,
    padding: u64,
    out: &mut impl Write,
) -> GatkResult<()> {
    writeln!(out, "read_path\tfinalize").map_err(|e| GatkError::generic(format!("write: {e}")))?;
    crate::assembly_region_stages_dump::dump_assembly_region_assembly_stages_tsv(
        ref_fasta,
        bam_path,
        interval_cli,
        padding,
        out,
    )
}

pub fn default_finalize_reads_padding() -> u64 {
    GATK_DEFAULT_ASSEMBLY_REGION_PADDING
}
