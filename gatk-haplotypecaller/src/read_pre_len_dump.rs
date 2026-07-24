//! PRE.2 gate dump: `read-pre-len`.

use crate::read_pre_len::{passes_read_length_filter, unclipped_read_length};
use gatk_common::{GatkError, GatkResult};
use rust_htslib::bam;
use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::Read as _;
use std::io::Write;
use std::path::Path;

pub fn dump_read_pre_len_tsv(alignment_path: &Path, out: &mut impl Write) -> GatkResult<()> {
    writeln!(
        out,
        "qname\tcigar\tread_length\tunclipped_length\tpasses_len_filter"
    )
    .map_err(|e| GatkError::generic(format!("write: {e}")))?;

    let mut reader = bam::Reader::from_path(alignment_path)
        .map_err(|e| GatkError::generic(format!("open: {e}")))?;
    for res in reader.records() {
        let rec = res.map_err(|e| GatkError::generic(format!("read: {e}")))?;
        let unclipped = unclipped_read_length(&rec);
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}",
            String::from_utf8_lossy(rec.qname()),
            format_cigar(&rec),
            rec.seq().as_bytes().len(),
            unclipped,
            passes_read_length_filter(&rec)
        )
        .map_err(|e| GatkError::generic(format!("write: {e}")))?;
    }
    Ok(())
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
