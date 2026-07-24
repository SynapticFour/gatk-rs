//! GATK `GVCFWriter` + `GATKVCFHeaderLines` slice.

use crate::genotyping::{
    build_gvcf_blocks_hc_emit, gvcf_block_to_record_fields, EmitMode, GvcfBlockRecordFields,
    ReferenceConfidenceLocus,
};
#[cfg(any(feature = "dev-dumps", test))]
use gatk_common::GatkError;
use gatk_common::GatkResult;
#[cfg(any(feature = "dev-dumps", test))]
use std::io::Write;

/// Default GQB bands aligned with GATK HC gVCF (`--gqb` defaults, parity subset).
pub const GATK_HC_DEFAULT_GQB: &[i32] = &[
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,
    27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50,
    51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 70, 80, 90, 99,
];

/// Agreed HC gVCF header lines for parity v1 (subset of `GATKVCFHeaderLines`).
pub fn gatk_hc_gvcf_header_lines(contig: &str, length: u64) -> Vec<String> {
    let mut lines = vec![
        format!("##reference=file://{contig}.fa"),
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">".to_string(),
        "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">".to_string(),
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">".to_string(),
        "##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description=\"Minimum DP\">".to_string(),
        "##FORMAT=<ID=MAX_DP,Number=1,Type=Integer,Description=\"Maximum DP\">".to_string(),
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End coordinate\">".to_string(),
        format!("##contig=<ID={contig},length={length}>"),
    ];
    for band in GATK_HC_DEFAULT_GQB {
        lines.push(format!("##GQB={band}"));
    }
    lines
}

/// GVCF writer configuration.
/// # Invariants
/// `gq_bands` are ascending GQB thresholds (HC defaults from [`GATK_HC_DEFAULT_GQB`]).
/// Contig length used for header contig lines.
/// # Ownership
/// Owns sample/contig names and GQB band vector.
/// # Mutation
/// Snapshot for writer construction.
/// # Biological assumptions
/// Configures single-sample gVCF block emission for reference confidence.
/// # Java equivalence
/// GATK `GVCFWriter` / HC `--gqb` configuration slice.
#[derive(Debug, Clone)]
pub struct GvcfWriterConfig {
    pub sample_name: String,
    pub contig: String,
    pub contig_length: u64,
    pub gq_bands: Vec<i32>,
}

impl Default for GvcfWriterConfig {
    fn default() -> Self {
        Self {
            sample_name: "SAMPLE".to_string(),
            contig: "chr1".to_string(),
            contig_length: 1_000_000,
            gq_bands: GATK_HC_DEFAULT_GQB.to_vec(),
        }
    }
}

/// In-memory gVCF writer (block records only; no streaming I/O yet).
/// # Invariants
/// `records` grow monotonically via reference-confidence write helpers.
/// Header lines match [`gatk_hc_gvcf_header_lines`] for the configured contig.
/// # Ownership
/// Owns header lines and block record fields.
/// # Mutation
/// Append-only record buffer during write calls.
/// # Biological assumptions
/// Emits compressed hom-ref blocks (and related gVCF fields) for one sample.
/// # Java equivalence
/// GATK `GVCFWriter` + `GATKVCFHeaderLines` subset (in-memory parity slice).
#[derive(Debug, Clone, Default)]
pub struct GvcfWriter {
    pub header_lines: Vec<String>,
    pub records: Vec<GvcfBlockRecordFields>,
}

impl GvcfWriter {
    pub fn with_config(config: &GvcfWriterConfig) -> Self {
        Self {
            header_lines: gatk_hc_gvcf_header_lines(&config.contig, config.contig_length),
            records: Vec::new(),
        }
    }

    pub fn write_reference_confidence_loci(
        &mut self,
        loci: &[ReferenceConfidenceLocus],
        gq_bands: &[i32],
    ) -> GatkResult<()> {
        let blocks = build_gvcf_blocks_hc_emit(loci, gq_bands)?;
        for block in &blocks {
            self.records.push(gvcf_block_to_record_fields(block)?);
        }
        Ok(())
    }

    pub fn emit_mode_for_writer(&self) -> EmitMode {
        EmitMode::Gvcf
    }
}

/// Serialize header + block records to a parity TSV (stable gate format).
#[cfg(any(feature = "dev-dumps", test))]
pub fn dump_gvcf_writer_state_tsv(writer: &GvcfWriter, out: &mut impl Write) -> GatkResult<()> {
    writeln!(out, "header_line_count\t{}", writer.header_lines.len())
        .map_err(|e| GatkError::generic(format!("write: {e}")))?;
    for (i, line) in writer.header_lines.iter().enumerate() {
        writeln!(out, "header_{i}\t{line}")
            .map_err(|e| GatkError::generic(format!("write: {e}")))?;
    }
    writeln!(out, "record_count\t{}", writer.records.len())
        .map_err(|e| GatkError::generic(format!("write: {e}")))?;
    writeln!(out, "pos\tend\tmin_dp\tmax_dp\tgq_band_upper\tmin_rgq")
        .map_err(|e| GatkError::generic(format!("write: {e}")))?;
    for rec in &writer.records {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}",
            rec.start_1based, rec.end_info, rec.min_dp, rec.max_dp, rec.gq_band_upper, rec.min_rgq
        )
        .map_err(|e| GatkError::generic(format!("write: {e}")))?;
    }
    Ok(())
}
