//! Bind alignments to traversal tiles.

use crate::active_region::TraversalTile;
use crate::read_model::{passes_hc_read_filters_with_header, ReadFilterParams};
use gatk_common::{GatkError, GatkResult};
use rust_htslib::bam;
use rust_htslib::bam::Read as _;
use std::path::Path;

/// 0-based half-open interval `[start, end)` on the reference for a mapped primary record.
fn reference_span0(
    record: &bam::Record,
    header: &bam::HeaderView,
    filters: &ReadFilterParams,
) -> Option<(i64, i64)> {
    if !passes_hc_read_filters_with_header(record, header, filters) {
        return None;
    }
    let start = record.pos();
    let end = record.cigar().end_pos();
    if end <= start {
        return None;
    }
    Some((start, end))
}

/// Tile as 0-based half-open from 1-based inclusive GATK-style `[start, end]`.
fn tile_span0(tile: &TraversalTile) -> (i64, i64) {
    let ts = tile.start.saturating_sub(1) as i64;
    let te = tile.end as i64; // exclusive 0-based end (1-based inclusive `end` → next base index)
    (ts, te)
}

pub fn half_open_overlaps(a0: i64, a1: i64, b0: i64, b1: i64) -> bool {
    a0 < b1 && b0 < a1
}

/// GATK-style 1-based inclusive `[start1, end1]` → 0-based half-open `[lo, hi)` on the reference.
pub fn closed_interval_1based_to_ref_span0(start1: u64, end1: u64) -> (i64, i64) {
    let lo = start1.saturating_sub(1) as i64;
    let hi = end1 as i64;
    (lo, hi)
}

/// GATK `IntervalUtils.isAfter(loc, span, dict)` for a 1-based locus vs a closed padded span.
pub fn locus_is_strictly_after_closed_interval_1based(pos1: u64, span_end1: u64) -> bool {
    pos1 > span_end1
}

/// GATK `IntervalUtils.isAfter(read, span, dict)` — alignment **start** (1-based) is past `span` end.
pub fn record_is_strictly_after_closed_interval_1based(
    record: &bam::Record,
    header: &bam::HeaderView,
    contig: &str,
    span_end1: u64,
    filters: &ReadFilterParams,
) -> bool {
    let rn = String::from_utf8_lossy(header.tid2name(record.tid() as u32));
    if rn != contig {
        return true;
    }
    let Some((r0, _)) = reference_span0(record, header, filters) else {
        return true;
    };
    let start1 = r0.saturating_add(1) as u64;
    start1 > span_end1
}

/// 1-based inclusive alignment end from CIGAR (0 if unmapped/empty).
pub fn record_alignment_end_1based(
    record: &bam::Record,
    header: &bam::HeaderView,
    filters: &ReadFilterParams,
) -> Option<u64> {
    let (r0, r1) = reference_span0(record, header, filters)?;
    let end1 = r1.max(r0 + 1) as u64;
    Some(end1)
}

/// 1-based inclusive alignment start.
pub fn record_alignment_start_1based(
    record: &bam::Record,
    header: &bam::HeaderView,
    filters: &ReadFilterParams,
) -> Option<u64> {
    reference_span0(record, header, filters).map(|(r0, _)| r0.saturating_add(1) as u64)
}

/// True if the mapped primary record overlaps the closed 1-based interval on `contig`.
pub fn record_overlaps_closed_interval_1based(
    record: &bam::Record,
    header: &bam::HeaderView,
    contig: &str,
    start1: u64,
    end1: u64,
    filters: &ReadFilterParams,
) -> bool {
    let rn = String::from_utf8_lossy(header.tid2name(record.tid() as u32));
    if rn != contig {
        return false;
    }
    let Some((r0, r1)) = reference_span0(record, header, filters) else {
        return false;
    };
    let (i0, i1) = closed_interval_1based_to_ref_span0(start1, end1);
    half_open_overlaps(r0, r1, i0, i1)
}

/// Deterministic filtered ingress order for reads with aligned reference span.
/// The returned order follows file iteration order from HTSlib and is used as a
/// contract for ingress reproducibility checks.
pub fn filtered_read_iteration_order(
    bam_path: &Path,
    filters: &ReadFilterParams,
) -> GatkResult<Vec<(i32, i64, String)>> {
    let mut reader =
        bam::Reader::from_path(bam_path).map_err(|e| GatkError::generic(format!("{e}")))?;
    let header = reader.header().clone();
    let mut out = Vec::new();
    for res in reader.records() {
        let rec = res.map_err(|e| GatkError::generic(format!("{e}")))?;
        if reference_span0(&rec, &header, filters).is_some() {
            out.push((
                rec.tid(),
                rec.pos(),
                String::from_utf8_lossy(rec.qname()).into_owned(),
            ));
        }
    }
    Ok(out)
}

/// Count primary alignments whose reference span overlaps `tile` (linear BAM scan).
pub fn count_reads_overlapping_tile(
    bam_path: &Path,
    tile: &TraversalTile,
    filters: &ReadFilterParams,
) -> GatkResult<usize> {
    let mut reader =
        bam::Reader::from_path(bam_path).map_err(|e| GatkError::generic(format!("{e}")))?;
    let tid = reader.header().tid(tile.contig.as_bytes()).ok_or_else(|| {
        GatkError::argument(format!(
            "Contig {} not found in BAM header dictionary",
            tile.contig
        ))
    })? as i32;
    let (t0, t1) = tile_span0(tile);
    let mut n = 0usize;
    let header = reader.header().clone();
    for res in reader.records() {
        let rec = res.map_err(|e| GatkError::generic(format!("{e}")))?;
        if rec.tid() != tid {
            continue;
        }
        if let Some((r0, r1)) = reference_span0(&rec, &header, filters) {
            if half_open_overlaps(r0, r1, t0, t1) {
                n += 1;
            }
        }
    }
    Ok(n)
}

/// Sum [`count_reads_overlapping_tile`] across all tiles (read may be counted in multiple tiles).
pub fn total_read_tile_overlaps(
    bam_path: &Path,
    tiles: &[TraversalTile],
    filters: &ReadFilterParams,
) -> GatkResult<usize> {
    let mut sum = 0usize;
    for tile in tiles {
        sum += count_reads_overlapping_tile(bam_path, tile, filters)?;
    }
    Ok(sum)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn sample_bam() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../parity/fixtures/sample.bam")
    }

    #[test]
    fn sample_read_overlaps_chr1_tile() {
        let tile = TraversalTile {
            contig: "chr1".to_string(),
            start: 1,
            end: 32,
        };
        let n = count_reads_overlapping_tile(&sample_bam(), &tile, &ReadFilterParams::default())
            .unwrap();
        assert!(n >= 1);
    }

    #[test]
    fn filtered_iteration_order_is_stable_across_calls() {
        let p = ReadFilterParams {
            min_mapping_quality: 0,
            exclude_duplicates: false,
            exclude_secondary: false,
            exclude_supplementary: false,
        };
        let first = filtered_read_iteration_order(&sample_bam(), &p).unwrap();
        let second = filtered_read_iteration_order(&sample_bam(), &p).unwrap();
        assert_eq!(first, second);
    }
}
