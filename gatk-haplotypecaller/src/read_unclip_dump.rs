//! PRE.1 gate dump: `read-pre-softclip`.

use crate::read_unclip::{
    apply_hc_softclip_pre_step, HcSoftclipPolicy, ORIGINAL_SOFTCLIP_END_TAG,
    ORIGINAL_SOFTCLIP_START_TAG,
};
use gatk_common::{GatkError, GatkResult};
use rust_htslib::bam;
use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::Read as _;
use std::io::Write;
use std::path::Path;

pub fn dump_read_pre_softclip_tsv(
    alignment_path: &Path,
    dont_use_soft_clipped_bases: bool,
    override_softclip_fragment_check: bool,
    out: &mut impl Write,
) -> GatkResult<()> {
    let policy = HcSoftclipPolicy {
        dont_use_soft_clipped_bases,
        override_softclip_fragment_check,
    };
    writeln!(
        out,
        "qname\tflags\tfragment_length\tcigar_in\tcigar_out\tseq_len_in\tseq_len_out\taction\tos\toe"
    )
    .map_err(|e| GatkError::generic(format!("write: {e}")))?;

    let mut reader = bam::Reader::from_path(alignment_path)
        .map_err(|e| GatkError::generic(format!("open: {e}")))?;
    for res in reader.records() {
        let rec = res.map_err(|e| GatkError::generic(format!("read: {e}")))?;
        let cigar_in = format_cigar(&rec);
        let len_in = rec.seq().as_bytes().len();
        let (mut out_rec, action, tags) = apply_hc_softclip_pre_step(&rec, &policy);
        if let Some((os, oe)) = tags {
            out_rec.remove_aux(ORIGINAL_SOFTCLIP_START_TAG).ok();
            out_rec.remove_aux(ORIGINAL_SOFTCLIP_END_TAG).ok();
            out_rec
                .push_aux(
                    ORIGINAL_SOFTCLIP_START_TAG,
                    rust_htslib::bam::record::Aux::I32(os),
                )
                .map_err(|e| GatkError::generic(e.to_string()))?;
            out_rec
                .push_aux(
                    ORIGINAL_SOFTCLIP_END_TAG,
                    rust_htslib::bam::record::Aux::I32(oe),
                )
                .map_err(|e| GatkError::generic(e.to_string()))?;
        }
        let (os_s, oe_s) = aux_tags(&out_rec);
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            String::from_utf8_lossy(rec.qname()),
            rec.flags(),
            rec.insert_size(),
            cigar_in,
            format_cigar(&out_rec),
            len_in,
            out_rec.seq().as_bytes().len(),
            action,
            os_s,
            oe_s,
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

fn aux_tags(rec: &bam::Record) -> (String, String) {
    let os = rec
        .aux(ORIGINAL_SOFTCLIP_START_TAG)
        .ok()
        .and_then(|a| i32_from_aux(a).ok())
        .map(|v| v.to_string())
        .unwrap_or_default();
    let oe = rec
        .aux(ORIGINAL_SOFTCLIP_END_TAG)
        .ok()
        .and_then(|a| i32_from_aux(a).ok())
        .map(|v| v.to_string())
        .unwrap_or_default();
    (os, oe)
}

fn i32_from_aux(aux: rust_htslib::bam::record::Aux<'_>) -> Result<i32, ()> {
    match aux {
        rust_htslib::bam::record::Aux::I32(v) => Ok(v),
        _ => Err(()),
    }
}
