//! PRE.3 gate dump: `read-pre-mq`.

use crate::read_pre_mq::{passes_assembly_mq_filter, GATK_ASSEMBLY_MQ_FILTER_THRESHOLD};
use gatk_common::{GatkError, GatkResult};
use rust_htslib::bam;
use rust_htslib::bam::Read as _;
use std::io::Write;
use std::path::Path;

pub fn dump_read_pre_mq_tsv(
    alignment_path: &Path,
    mq_threshold: u8,
    out: &mut impl Write,
) -> GatkResult<()> {
    writeln!(out, "qname\tmapq\tmq_threshold\tpasses_mq_filter")
        .map_err(|e| GatkError::generic(format!("write: {e}")))?;

    let mut reader = bam::Reader::from_path(alignment_path)
        .map_err(|e| GatkError::generic(format!("open: {e}")))?;
    for res in reader.records() {
        let rec = res.map_err(|e| GatkError::generic(format!("read: {e}")))?;
        writeln!(
            out,
            "{}\t{}\t{}\t{}",
            String::from_utf8_lossy(rec.qname()),
            rec.mapq(),
            mq_threshold,
            passes_assembly_mq_filter(&rec, mq_threshold)
        )
        .map_err(|e| GatkError::generic(format!("write: {e}")))?;
    }
    Ok(())
}

pub fn dump_read_pre_mq_tsv_default(alignment_path: &Path, out: &mut impl Write) -> GatkResult<()> {
    dump_read_pre_mq_tsv(alignment_path, GATK_ASSEMBLY_MQ_FILTER_THRESHOLD, out)
}
