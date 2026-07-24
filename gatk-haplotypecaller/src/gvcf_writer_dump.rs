//! GVCF writer / header parity dumps.

use crate::genotyping::{
    build_gvcf_blocks, gvcf_block_to_record_fields, validate_joint_compatibility_gvcf_records,
    ReferenceConfidenceLocus,
};
use crate::gvcf_writer::{
    dump_gvcf_writer_state_tsv, gatk_hc_gvcf_header_lines, GvcfWriter, GvcfWriterConfig,
    GATK_HC_DEFAULT_GQB,
};
use gatk_common::{GatkError, GatkResult};
use std::io::Write;
use std::path::Path;

/// H.2.1 — agreed `GATKVCFHeaderLines` subset.
pub fn dump_gvcf_header_tsv(
    contig: &str,
    contig_length: u64,
    out: &mut impl Write,
) -> GatkResult<()> {
    let lines = gatk_hc_gvcf_header_lines(contig, contig_length);
    writeln!(out, "contig\t{contig}")?;
    writeln!(out, "contig_length\t{contig_length}")?;
    writeln!(out, "header_line_count\t{}", lines.len())?;
    for (i, line) in lines.iter().enumerate() {
        writeln!(out, "header_{i}\t{line}")?;
    }
    Ok(())
}

/// H.2.1 — build writer block records from p8-style locus fixture TSV.
pub fn dump_gvcf_writer_from_loci_fixture_tsv(
    fixture_path: &Path,
    out: &mut impl Write,
) -> GatkResult<()> {
    let raw = std::fs::read_to_string(fixture_path)
        .map_err(|e| GatkError::io(format!("read {}", fixture_path.display()), e))?;
    let mut loci: Vec<ReferenceConfidenceLocus> = Vec::new();
    let mut case_id = String::new();
    for line in raw.lines() {
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = t.split('\t').collect();
        if cols.len() < 4 {
            return Err(GatkError::argument(format!("bad fixture row: {t}")));
        }
        case_id = cols[0].to_string();
        loci.push(ReferenceConfidenceLocus {
            position_1based: cols[1]
                .parse()
                .map_err(|_| GatkError::argument("bad pos"))?,
            gq: cols[2].parse().map_err(|_| GatkError::argument("bad gq"))?,
            dp: cols[3].parse().map_err(|_| GatkError::argument("bad dp"))?,
        });
    }
    let bands: Vec<i32> = GATK_HC_DEFAULT_GQB.to_vec();
    let blocks = build_gvcf_blocks(&loci, &bands)?;
    let mut writer = GvcfWriter::with_config(&GvcfWriterConfig {
        sample_name: case_id.clone(),
        ..GvcfWriterConfig::default()
    });
    for block in &blocks {
        writer.records.push(gvcf_block_to_record_fields(block)?);
    }
    writeln!(out, "case_id\t{case_id}")?;
    dump_gvcf_writer_state_tsv(&writer, out)
}

/// H-D01/D02 — merged gVCF block pseudo-VCF lines + joint-compat (L5 scaffold).
pub fn dump_gvcf_l5_merged_tsv(fixture_path: &Path, out: &mut impl Write) -> GatkResult<()> {
    let raw = std::fs::read_to_string(fixture_path)
        .map_err(|e| GatkError::io(format!("read {}", fixture_path.display()), e))?;
    let mut loci: Vec<ReferenceConfidenceLocus> = Vec::new();
    let contig = "chr1";
    for line in raw.lines() {
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = t.split('\t').collect();
        if cols.len() < 4 {
            return Err(GatkError::argument(format!("bad fixture row: {t}")));
        }
        loci.push(ReferenceConfidenceLocus {
            position_1based: cols[1]
                .parse()
                .map_err(|_| GatkError::argument("bad pos"))?,
            gq: cols[2].parse().map_err(|_| GatkError::argument("bad gq"))?,
            dp: cols[3].parse().map_err(|_| GatkError::argument("bad dp"))?,
        });
    }
    let bands: Vec<i32> = GATK_HC_DEFAULT_GQB.to_vec();
    let blocks = build_gvcf_blocks(&loci, &bands)?;
    let records: Vec<_> = blocks
        .iter()
        .map(gvcf_block_to_record_fields)
        .collect::<Result<_, _>>()?;
    let compat = validate_joint_compatibility_gvcf_records(&records)?;
    writeln!(out, "contig\t{contig}")?;
    writeln!(out, "record_count\t{}", records.len())?;
    writeln!(out, "joint_compatible\t{}", compat.compatible)?;
    for (i, rec) in records.iter().enumerate() {
        writeln!(
            out,
            "vcf_line_{i}\t{contig}\t{}\t.\t<NON_REF>\t.\t.\tEND={}\tMIN_DP={}\tMAX_DP={}\tGQ_BAND={}\tMIN_RGQ={}",
            rec.start_1based,
            rec.end_info,
            rec.min_dp,
            rec.max_dp,
            rec.gq_band_upper,
            rec.min_rgq,
        )?;
    }
    Ok(())
}
