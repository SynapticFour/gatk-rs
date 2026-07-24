//! PRE.4 gate dump: `read-pre-overlap`.

use crate::fragment_overlap::{
    clean_overlapping_read_pairs, format_quals, overlapping_pairs_indices,
};
use gatk_common::{GatkError, GatkResult};
use rust_htslib::bam;
use rust_htslib::bam::Read as _;
use std::collections::HashMap;
use std::io::Write;
use std::path::Path;

pub fn dump_read_pre_overlap_tsv(alignment_path: &Path, out: &mut impl Write) -> GatkResult<()> {
    writeln!(out, "qname\tqual_in\tqual_out\toverlap_pair")
        .map_err(|e| GatkError::generic(format!("write: {e}")))?;

    let mut reader = bam::Reader::from_path(alignment_path)
        .map_err(|e| GatkError::generic(format!("open: {e}")))?;
    let mut records = Vec::new();
    let mut qual_in = HashMap::new();
    for res in reader.records() {
        let rec = res.map_err(|e| GatkError::generic(format!("read: {e}")))?;
        qual_in.insert(rec.qname().to_vec(), format_quals(&rec));
        records.push(rec);
    }
    records.sort_by_key(|r| (r.tid(), r.pos()));
    let pair_indices = overlapping_pairs_indices(&records)?;
    let mut in_pair = std::collections::HashSet::new();
    for (a, b) in &pair_indices {
        in_pair.insert(records[*a].qname().to_vec());
        in_pair.insert(records[*b].qname().to_vec());
    }
    clean_overlapping_read_pairs(&mut records, true)?;
    for rec in &records {
        let qname = rec.qname();
        let qual_in_s = qual_in.get(qname).map(String::as_str).unwrap_or("");
        writeln!(
            out,
            "{}\t{}\t{}\t{}",
            String::from_utf8_lossy(qname),
            qual_in_s,
            format_quals(rec),
            in_pair.contains(qname)
        )
        .map_err(|e| GatkError::generic(format!("write: {e}")))?;
    }
    Ok(())
}
