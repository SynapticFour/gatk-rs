//! HC read-filter decision dump.

use crate::read_model::{
    standard_hc_read_filter_failure_index, GATK_HC_DEFAULT_MIN_MAPPING_QUALITY,
    STANDARD_HC_READ_FILTER_JAVA_NAMES,
};
use gatk_common::{GatkError, GatkResult};
use rust_htslib::bam;
use rust_htslib::bam::Read as _;
use std::io::Write;
use std::path::Path;

/// L2 delimiter between per-read rows and the `CountingReadFilter`-style summary (matches Java gate dump).
pub const HC_READ_FILTER_COUNT_SECTION: &str = "---HC_READ_FILTER_COUNTS---";

/// TSV: `qname`, `flags`, `mapq`, `passes_hc_filter` (`true`/`false`).
/// Matches `HcFullParityGateDump read-filters`: `HaplotypeCallerEngine.makeStandardHCReadFilters`
/// at default MQ threshold ([`GATK_HC_DEFAULT_MIN_MAPPING_QUALITY`]).
/// Appends [`HC_READ_FILTER_COUNT_SECTION`] plus `filter\tfiltered_count` rows mirroring Java
/// `CountingReadFilter` leaf counts (AND short-circuit order).
pub fn dump_hc_read_filter_tsv(alignment_path: &Path, out: &mut impl Write) -> GatkResult<()> {
    writeln!(out, "qname\tflags\tmapq\tpasses_hc_filter")
        .map_err(|e| GatkError::generic(format!("write: {e}")))?;
    let mut reader = bam::Reader::from_path(alignment_path)
        .map_err(|e| GatkError::generic(format!("open: {e}")))?;
    let header = reader.header().clone();
    let mut counts = [0u64; STANDARD_HC_READ_FILTER_JAVA_NAMES.len()];
    for res in reader.records() {
        let rec = res.map_err(|e| GatkError::generic(format!("read: {e}")))?;
        let fail = standard_hc_read_filter_failure_index(
            &rec,
            &header,
            GATK_HC_DEFAULT_MIN_MAPPING_QUALITY,
        );
        let pass = fail.is_none();
        if let Some(i) = fail {
            counts[i] += 1;
        }
        writeln!(
            out,
            "{}\t{}\t{}\t{}",
            String::from_utf8_lossy(rec.qname()),
            rec.flags(),
            rec.mapq(),
            if pass { "true" } else { "false" }
        )
        .map_err(|e| GatkError::generic(format!("write: {e}")))?;
    }
    writeln!(out, "{HC_READ_FILTER_COUNT_SECTION}")
        .map_err(|e| GatkError::generic(format!("write: {e}")))?;
    writeln!(out, "filter\tfiltered_count")
        .map_err(|e| GatkError::generic(format!("write: {e}")))?;
    for (i, name) in STANDARD_HC_READ_FILTER_JAVA_NAMES.iter().enumerate() {
        writeln!(out, "{name}\t{}", counts[i])
            .map_err(|e| GatkError::generic(format!("write: {e}")))?;
    }
    Ok(())
}
