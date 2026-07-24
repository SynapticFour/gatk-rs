//! Per-read pileup element dump for RCM / `isActive` parity (ASM act-2).

use crate::activity_scoring::is_alt_before_assembly;
use crate::locus_iterator::{IntervalLocusIterator, LocusPileupState};
use crate::pileup_element::pileup_element_flags_at_ref;
use crate::read_model::ReadFilterParams;
use crate::read_transformer::{apply_shard_read_pipeline, ShardReadPipelineConfig};
use crate::walker::make_read_shards;
use gatk_common::{GatkError, GatkResult};
use gatk_core::reference::{parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary};
use std::io::Write;
use std::path::Path;

/// One row per pileup observation: `contig`, `pos`, `read`, `read_base`, `ref_base`, `qual`, flags, `is_alt`.
pub fn dump_locus_pileup_detail_tsv(
    reference_fasta: &Path,
    bam_path: &Path,
    intervals_cli: &str,
    assembly_region_padding: u64,
    out: &mut impl Write,
    read_filters: &ReadFilterParams,
) -> GatkResult<()> {
    writeln!(
        out,
        "contig\tpos\tread\tread_base\tref_base\tqual\tis_del\tbefore_del\tafter_del\tbefore_ins\tafter_ins\tnext_to_sc\tis_alt"
    )
    .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;

    let dict = SequenceDictionary::from_fasta_path(reference_fasta)?;
    let specs = parse_intervals_cli_string(&dict, intervals_cli)?;
    let pipeline = ShardReadPipelineConfig::gatk_haplotype_caller_production();
    let mut ref_cache = ReferenceWindowCache::new(reference_fasta.to_path_buf(), 4);
    let mut rng = crate::read_downsample::GatkJavaRng::reset_gatk_default();

    for spec in &specs {
        let (c, s, e) = spec
            .resolve_closed_ends(&dict)
            .map_err(|e| GatkError::argument(e.to_string()))?;
        let shard = make_read_shards(&dict, std::slice::from_ref(spec), assembly_region_padding)?
            .into_iter()
            .find(|sh| sh.contig == c)
            .ok_or_else(|| GatkError::argument(format!("no shard for {c}")))?;
        let (header, mut records) =
            crate::assembly_region_iterator::load_records_for_shard_raw(bam_path, &shard)?;
        apply_shard_read_pipeline(
            &mut records,
            Some(&header),
            read_filters,
            &pipeline,
            &mut rng,
        )?;
        let ref_window = ref_cache
            .get_interval_bytes(&dict, &c, s, e)
            .map_err(|e| GatkError::generic(e.to_string()))?;
        let mut pileup_state = LocusPileupState::from_records(&records, &header, &c, read_filters);
        for pos1 in IntervalLocusIterator::from_closed_interval(s, e) {
            let ref_base = *ref_window
                .get((pos1 - s) as usize)
                .ok_or_else(|| GatkError::argument("reference window index out of range"))?;
            pileup_state.advance_to(&records, read_filters, pos1)?;
            let ref_pos0 = pos1.saturating_sub(1) as i64;
            for &idx in &pileup_state.active {
                let rec = &records[idx];
                let cigar: Vec<_> = rec.cigar().iter().copied().collect();
                let seq = rec.seq();
                let qual = rec.qual();
                let Some(flags) =
                    pileup_element_flags_at_ref(rec.pos(), &cigar, &seq.as_bytes(), qual, ref_pos0)
                else {
                    continue;
                };
                let is_alt = is_alt_before_assembly(
                    flags.read_base,
                    ref_base,
                    flags.is_deletion,
                    flags.is_before_deletion_start,
                    flags.is_after_deletion_end,
                    flags.is_before_insertion,
                    flags.is_after_insertion,
                    flags.is_next_to_soft_clip,
                );
                let read_base = if flags.is_deletion {
                    '-'
                } else {
                    flags.read_base as char
                };
                writeln!(
                    out,
                    "{c}\t{pos1}\t{}\t{read_base}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{is_alt}",
                    String::from_utf8_lossy(rec.qname()),
                    ref_base as char,
                    flags.qual,
                    flags.is_deletion,
                    flags.is_before_deletion_start,
                    flags.is_after_deletion_end,
                    flags.is_before_insertion,
                    flags.is_after_insertion,
                    flags.is_next_to_soft_clip,
                )
                .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
            }
        }
    }
    Ok(())
}
