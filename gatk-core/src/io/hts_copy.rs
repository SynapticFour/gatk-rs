//! SAM/BAM copy and validation using HTSlib (real binary I/O).
//! The in-crate `BamReader` / `BamWriter` paths are intentionally lightweight;
//! alignment round-trips and BAM validation for parity use this module.

use crate::io::bam::{BamHeader, ReferenceSequence};
use crate::reference::SequenceDictionary;
use gatk_common::{GatkError, GatkResult};
use rust_htslib::bam::Read as HtsRead;
use rust_htslib::bam::{self, Format};
use std::path::Path;

fn output_format(path: &Path) -> Format {
    match path
        .extension()
        .and_then(|e| e.to_str())
        .map(|s| s.to_ascii_lowercase())
        .as_deref()
    {
        Some("sam") => Format::Sam,
        Some("cram") => Format::Cram,
        _ => Format::Bam,
    }
}

fn is_cram_path(path: &Path) -> bool {
    path.extension()
        .and_then(|e| e.to_str())
        .map(|s| s.eq_ignore_ascii_case("cram"))
        .unwrap_or(false)
}

/// Copy all alignment records from `input` to `output`, preserving the header.
/// Output format is inferred from `output`'s extension (`.sam` → SAM, otherwise BAM).
pub fn copy_alignments_with_htslib(input: &Path, output: &Path) -> GatkResult<u64> {
    copy_alignments_with_htslib_with_reference(input, output, None)
}

/// Copy all alignment records from `input` to `output`, preserving the header.
/// When either side is CRAM, `reference_fasta` must be set and readable.
pub fn copy_alignments_with_htslib_with_reference(
    input: &Path,
    output: &Path,
    reference_fasta: Option<&Path>,
) -> GatkResult<u64> {
    let needs_reference = is_cram_path(input) || is_cram_path(output);
    if needs_reference && reference_fasta.is_none() {
        return Err(GatkError::argument(
            "CRAM input/output requires --reference-fasta (or equivalent reference path)",
        ));
    }

    let mut reader = bam::Reader::from_path(input).map_err(|e| {
        GatkError::generic(format!("Failed to open alignment input {:?}: {}", input, e))
    })?;
    if let Some(reference_fasta) = reference_fasta {
        if is_cram_path(input) {
            reader.set_reference(reference_fasta).map_err(|e| {
                GatkError::generic(format!(
                    "Failed to set CRAM reference {:?} for input {:?}: {}",
                    reference_fasta, input, e
                ))
            })?;
        }
    }

    let header = bam::Header::from_template(reader.header());
    let fmt = output_format(output);

    let mut writer = bam::Writer::from_path(output, &header, fmt).map_err(|e| {
        GatkError::generic(format!(
            "Failed to open alignment output {:?}: {}",
            output, e
        ))
    })?;
    if let Some(reference_fasta) = reference_fasta {
        if is_cram_path(output) {
            writer.set_reference(reference_fasta).map_err(|e| {
                GatkError::generic(format!(
                    "Failed to set CRAM reference {:?} for output {:?}: {}",
                    reference_fasta, output, e
                ))
            })?;
        }
    }

    let mut n = 0u64;
    for result in reader.records() {
        let record = result.map_err(|e| {
            GatkError::generic(format!(
                "Failed to read alignment record from {:?}: {}",
                input, e
            ))
        })?;
        writer.write(&record).map_err(|e| {
            GatkError::generic(format!(
                "Failed to write alignment record to {:?}: {}",
                output, e
            ))
        })?;
        n += 1;
    }
    Ok(n)
}

fn parse_sq_line(line: &str) -> GatkResult<ReferenceSequence> {
    let mut name = String::new();
    let mut length: u64 = 0;
    let mut md5 = None;

    for field in line.split('\t').skip(1) {
        if let Some((k, v)) = field.split_once(':') {
            match k {
                "SN" => name = v.to_string(),
                "LN" => {
                    length = v.parse().map_err(|_| {
                        GatkError::generic(format!(
                            "Invalid @SQ LN value in BAM header line: {line}"
                        ))
                    })?;
                }
                "M5" => md5 = Some(v.to_string()),
                _ => {}
            }
        }
    }

    if name.is_empty() {
        return Err(GatkError::generic(format!(
            "Missing SN in @SQ BAM header line: {line}"
        )));
    }

    Ok(ReferenceSequence {
        name,
        length,
        md5,
        assembly: None,
        uri: None,
        species: None,
    })
}

fn bam_header_from_hts_view(view: &bam::HeaderView) -> GatkResult<BamHeader> {
    let text = std::str::from_utf8(view.as_bytes())
        .map_err(|e| GatkError::generic(format!("BAM header is not valid UTF-8: {e}")))?;

    let mut header = BamHeader::default();
    for line in text.lines() {
        if line.starts_with("@SQ") {
            header.reference_sequences.push(parse_sq_line(line)?);
        }
    }
    Ok(header)
}

/// Normalize a SAM/BAM header line for stable comparisons: keep the record tag
/// (`@HD`, `@SQ`, …) and sort tab-separated optional fields lexicographically.
fn normalize_sam_header_line_fields(line: &str) -> String {
    let line = line.trim_end();
    let mut parts = line.splitn(2, '\t');
    let tag = parts.next().unwrap_or("");
    match parts.next() {
        None | Some("") => tag.to_string(),
        Some(rest) => {
            let mut fields: Vec<&str> = rest.split('\t').filter(|s| !s.is_empty()).collect();
            fields.sort_unstable();
            if fields.is_empty() {
                tag.to_string()
            } else {
                format!("{tag}\t{}", fields.join("\t"))
            }
        }
    }
}

/// Return canonical `@HD`, `@SQ`, and `@RG` header lines for an alignment file.
/// This mirrors the subset used in `scripts/parity/compare_bam_alignment_parity.py`
/// (`@HD` / `@SQ` / `@RG` only) and applies deterministic field ordering so
/// roundtrips through HTSlib can be compared without brittle raw-line equality.
pub fn alignment_header_canonical_hd_sq_rg(
    input: &Path,
) -> GatkResult<(Vec<String>, Vec<String>, Vec<String>)> {
    let reader = bam::Reader::from_path(input).map_err(|e| {
        GatkError::generic(format!("Failed to open alignment input {:?}: {}", input, e))
    })?;
    let text = std::str::from_utf8(reader.header().as_bytes())
        .map_err(|e| GatkError::generic(format!("BAM header is not valid UTF-8: {e}")))?;

    let mut hd = Vec::new();
    let mut sq = Vec::new();
    let mut rg = Vec::new();
    for line in text.lines() {
        let line = line.trim_end();
        if line.is_empty() {
            continue;
        }
        let tag = line.split('\t').next().unwrap_or("");
        let normalized = normalize_sam_header_line_fields(line);
        match tag {
            "@HD" => hd.push(normalized),
            "@SQ" => sq.push(normalized),
            "@RG" => rg.push(normalized),
            _ => {}
        }
    }
    hd.sort();
    sq.sort();
    rg.sort();
    Ok((hd, sq, rg))
}

/// Validate a BAM (or SAM) file by reading the header and all records via HTSlib.
/// When `dictionary` is set, `@SQ` lines are checked against the reference dictionary.
pub fn validate_bam_file(input: &Path, dictionary: Option<&SequenceDictionary>) -> GatkResult<u64> {
    let mut reader = bam::Reader::from_path(input).map_err(|e| {
        GatkError::generic(format!("Failed to open alignment file {:?}: {}", input, e))
    })?;

    let hdr = bam_header_from_hts_view(reader.header())?;
    if let Some(dict) = dictionary {
        dict.validate_bam_header(&hdr)?;
    }

    let mut count = 0u64;
    for record in reader.records() {
        let _ = record.map_err(|e| {
            GatkError::generic(format!("Invalid alignment record in {:?}: {}", input, e))
        })?;
        count += 1;
    }

    if count == 0 {
        return Err(GatkError::generic(format!(
            "Alignment file {:?} contains no records",
            input
        )));
    }

    Ok(count)
}

/// Count records in a 1-based inclusive region using BAM index random access.
pub fn count_records_in_region_indexed(
    input: &Path,
    contig: &str,
    start_1based: u64,
    end_1based_inclusive: u64,
) -> GatkResult<u64> {
    if start_1based == 0 {
        return Err(GatkError::argument("Region start must be >= 1"));
    }
    if start_1based > end_1based_inclusive {
        return Err(GatkError::argument("Region start must be <= end"));
    }

    let mut reader = bam::IndexedReader::from_path(input).map_err(|e| {
        GatkError::generic(format!(
            "Failed to open indexed alignment file {:?}: {}",
            input, e
        ))
    })?;
    let tid = reader
        .header()
        .tid(contig.as_bytes())
        .ok_or_else(|| GatkError::argument(format!("Contig not found in BAM header: {contig}")))?;
    let start0 = start_1based - 1;
    let end0 = end_1based_inclusive;
    reader.fetch((tid, start0, end0)).map_err(|e| {
        GatkError::generic(format!(
            "Failed indexed fetch on {:?} for {}:{}-{}: {}",
            input, contig, start_1based, end_1based_inclusive, e
        ))
    })?;

    let mut count = 0u64;
    for result in reader.records() {
        let _ = result.map_err(|e| {
            GatkError::generic(format!(
                "Failed while iterating indexed region on {:?}: {}",
                input, e
            ))
        })?;
        count += 1;
    }
    Ok(count)
}

/// Return QNAMEs for records in a 1-based inclusive indexed region.
pub fn qnames_in_region_indexed(
    input: &Path,
    contig: &str,
    start_1based: u64,
    end_1based_inclusive: u64,
) -> GatkResult<Vec<String>> {
    if start_1based == 0 {
        return Err(GatkError::argument("Region start must be >= 1"));
    }
    if start_1based > end_1based_inclusive {
        return Err(GatkError::argument("Region start must be <= end"));
    }

    let mut reader = bam::IndexedReader::from_path(input).map_err(|e| {
        GatkError::generic(format!(
            "Failed to open indexed alignment file {:?}: {}",
            input, e
        ))
    })?;
    let tid = reader
        .header()
        .tid(contig.as_bytes())
        .ok_or_else(|| GatkError::argument(format!("Contig not found in BAM header: {contig}")))?;
    let start0 = start_1based - 1;
    let end0 = end_1based_inclusive;
    reader.fetch((tid, start0, end0)).map_err(|e| {
        GatkError::generic(format!(
            "Failed indexed fetch on {:?} for {}:{}-{}: {}",
            input, contig, start_1based, end_1based_inclusive, e
        ))
    })?;

    let mut out = Vec::new();
    for result in reader.records() {
        let rec = result.map_err(|e| {
            GatkError::generic(format!(
                "Failed while iterating indexed region on {:?}: {}",
                input, e
            ))
        })?;
        out.push(String::from_utf8_lossy(rec.qname()).into_owned());
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn normalize_sam_header_line_fields_is_order_invariant() {
        let a = normalize_sam_header_line_fields("@HD\tSO:coordinate\tVN:1.6");
        let b = normalize_sam_header_line_fields("@HD\tVN:1.6\tSO:coordinate");
        assert_eq!(a, b);
        assert!(a.starts_with("@HD\t"));
    }

    #[test]
    fn copy_minimal_sam_roundtrip() {
        let dir = tempfile::tempdir().unwrap();
        let sam = dir.path().join("in.sam");
        let out_sam = dir.path().join("out.sam");
        let mut f = std::fs::File::create(&sam).unwrap();
        writeln!(
            f,
            "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:32\nr1\t0\tchr1\t1\t60\t32M\t*\t0\t0\tACGTACGTACGTACGTACGTACGTACGTACGT\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"
        )
        .unwrap();

        let n = copy_alignments_with_htslib(&sam, &out_sam).unwrap();
        assert_eq!(n, 1);
        let text = std::fs::read_to_string(&out_sam).unwrap();
        assert!(text.contains("r1\t0\tchr1"));
    }
}
