//! HC shard read-pipeline dump: IUPAC pre-transform, filter, post-transform.
//! Downsample is **not** included (isolated in D.2). Matches `HcFullParityGateDump read-shard-pipeline`.

use crate::dragen_mq::apply_dragen_mapping_quality_transform;
use crate::read_model::{
    passes_hc_read_filters_with_header, ReadFilterParams, GATK_HC_DEFAULT_MIN_MAPPING_QUALITY,
};
use crate::read_transformer::apply_iupac_strict_transform;
use gatk_common::{GatkError, GatkResult};
use rust_htslib::bam;
use rust_htslib::bam::Read as _;
use std::io::Write;
use std::path::Path;

fn seq_to_string(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).into_owned()
}

/// TSV per read: `qname`, `flags`, `mapq`, `seq_raw`, `seq_after_pre`, `passes_hc_filter`, `seq_after_post`.
/// Pipeline stages mirror `MultiIntervalLocalReadShard` without downsampler:
/// `makeStandardHCReadTransformer` → `makeStandardHCReadFilters` → identity post-transform.
pub fn dump_read_shard_pipeline_tsv(
    alignment_path: &Path,
    apply_dragen_mq: bool,
    out: &mut impl Write,
) -> GatkResult<()> {
    writeln!(
        out,
        "qname\tflags\tmapq\tseq_raw\tseq_after_pre\tpasses_hc_filter\tseq_after_post"
    )
    .map_err(|e| GatkError::generic(format!("write: {e}")))?;

    let filters = ReadFilterParams::gatk_standard_hc();
    let mut reader = bam::Reader::from_path(alignment_path)
        .map_err(|e| GatkError::generic(format!("open: {e}")))?;
    let header = reader.header().clone();

    for res in reader.records() {
        let mut rec = res.map_err(|e| GatkError::generic(format!("read: {e}")))?;
        let seq_raw = seq_to_string(&rec.seq().as_bytes());
        apply_iupac_strict_transform(std::slice::from_mut(&mut rec));
        let seq_after_pre = seq_to_string(&rec.seq().as_bytes());
        let pass = passes_hc_read_filters_with_header(&rec, &header, &filters);
        let mut post_rec = rec.clone();
        if apply_dragen_mq {
            apply_dragen_mapping_quality_transform(std::slice::from_mut(&mut post_rec));
        }
        let seq_after_post = seq_to_string(&post_rec.seq().as_bytes());
        let mapq_after = post_rec.mapq();
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            String::from_utf8_lossy(rec.qname()),
            rec.flags(),
            if apply_dragen_mq {
                mapq_after
            } else {
                rec.mapq()
            },
            seq_raw,
            seq_after_pre,
            if pass { "true" } else { "false" },
            seq_after_post,
        )
        .map_err(|e| GatkError::generic(format!("write: {e}")))?;
    }

    let _ = GATK_HC_DEFAULT_MIN_MAPPING_QUALITY;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dump_header_present() {
        let repo = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../parity/fixtures/hc-full-parity/d4/d4_iupac_mixed.sam");
        if !repo.exists() {
            return;
        }
        let mut buf = Vec::new();
        dump_read_shard_pipeline_tsv(&repo, false, &mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.starts_with("qname\tflags\tmapq\tseq_raw"));
        assert!(s.contains("r_mixed"));
        assert!(s.contains("ANNN"));
    }
}
