//! GATK `ReadTransformer` hooks for shard read preparation.
//! Java `MultiIntervalLocalReadShard` order:
//! pre-transform → read filter → post-transform → downsampler.

use crate::read_downsample::{
    apply_positional_downsampler, GatkJavaRng, PositionalDownsamplerConfig,
};
use crate::read_model::{passes_hc_read_filters_with_header, ReadFilterParams};
use gatk_common::GatkError;
use rust_htslib::bam;

/// Shard-level read pipeline matching GATK HC defaults when enabled.
/// # Invariants
/// Pipeline order matches Java: pre-transform → filter → post-transform → downsampler.
/// [`Self::parity_l1_goldens`] disables IUPAC transform and downsampling for L1 gates.
/// # Ownership
/// [`Clone`] config owns nested [`PositionalDownsamplerConfig`]; BAM records are mutated in place.
/// # Mutation
/// Config immutable; [`apply_shard_read_pipeline`] mutates record bases/MAPQ and may downsample.
/// # Biological assumptions
/// Strict IUPAC→N avoids ambiguous bases before HC filters; DRAGEN transform maps `XQ` to MAPQ when enabled.
/// # Java equivalence
/// GATK `MultiIntervalLocalReadShard` + `HaplotypeCallerEngine.makeStandardHCReadTransformer`.
#[derive(Debug, Clone)]
pub struct ShardReadPipelineConfig {
    /// `HaplotypeCallerEngine.makeStandardHCReadTransformer` (`IUPACReadTransformer` strict).
    pub apply_iupac_pre_transform: bool,
    /// HC post-transform: identity, or DRAGEN `XQ` → MAPQ when true.
    pub apply_dragen_mapq_transform: bool,
    pub downsample: PositionalDownsamplerConfig,
}

impl ShardReadPipelineConfig {
    /// B.2/B.3 L1 goldens: no IUPAC transform, downsampling off.
    pub fn parity_l1_goldens() -> Self {
        Self {
            apply_iupac_pre_transform: false,
            apply_dragen_mapq_transform: false,
            downsample: PositionalDownsamplerConfig::disabled(),
        }
    }

    /// HC production shard pipeline (B.4.1–B.4.3).
    pub fn gatk_haplotype_caller_production() -> Self {
        Self {
            apply_iupac_pre_transform: true,
            apply_dragen_mapq_transform: false,
            downsample: PositionalDownsamplerConfig::gatk_haplotype_caller_defaults(),
        }
    }
}

impl Default for ShardReadPipelineConfig {
    fn default() -> Self {
        Self::parity_l1_goldens()
    }
}

/// In-place strict IUPAC → `N` on read bases (`BaseUtils.convertIUPACtoN(..., true, false)`).
pub fn apply_iupac_strict_transform(records: &mut [bam::Record]) {
    use rust_htslib::bam::record::CigarString;
    for rec in records.iter_mut() {
        let qname = rec.qname().to_vec();
        let cigar = CigarString(rec.cigar().iter().copied().collect());
        let qual = rec.qual().to_vec();
        let mut seq = rec.seq().as_bytes().to_vec();
        convert_iupac_to_n_strict(&mut seq);
        rec.set(&qname, Some(&cigar), &seq, qual.as_slice());
    }
}

/// GATK `BaseUtils.convertIUPACtoN` with `errorOnBadReferenceBase=true`, `ignoreConversionOfFirstByte=false`.
fn convert_iupac_to_n_strict(bases: &mut [u8]) {
    for base in bases.iter_mut() {
        if iupac_maps_to_n(*base) {
            *base = b'N';
        }
    }
}

/// Mirrors GATK `baseIndexWithIupacMap`: A/C/G/T (and `*`→A in simple map) stay; IUPAC and N map to N.
fn iupac_maps_to_n(base: u8) -> bool {
    match base {
        b'A' | b'a' | b'C' | b'c' | b'G' | b'g' | b'T' | b't' => false,
        b'N' | b'n' | b'R' | b'r' | b'Y' | b'y' | b'M' | b'm' | b'K' | b'k' | b'W' | b'w'
        | b'S' | b's' | b'B' | b'b' | b'D' | b'd' | b'H' | b'h' | b'V' | b'v' => true,
        _ => false,
    }
}

/// Load one contig with the HC production shard pipeline (parity / activity / pileup dumps).
pub fn load_contig_records_hc_production(
    bam_path: &std::path::Path,
    contig: &str,
    filters: &ReadFilterParams,
    rng: &mut GatkJavaRng,
) -> Result<(rust_htslib::bam::HeaderView, Vec<bam::Record>), gatk_common::GatkError> {
    use crate::assembly_region_iterator::load_all_records_for_contig_raw;
    use gatk_common::GatkError;
    let (header, mut records) = load_all_records_for_contig_raw(bam_path, contig)
        .map_err(|e| GatkError::generic(e.to_string()))?;
    apply_shard_read_pipeline(
        &mut records,
        Some(&header),
        filters,
        &ShardReadPipelineConfig::gatk_haplotype_caller_production(),
        rng,
    )?;
    Ok((header, records))
}

/// Apply GATK shard iterator pipeline to a contig-local record list.
/// `header`: pass `Some` for Java-identical `WellformedReadFilter` + full chains when
/// [`ReadFilterParams::resolved_hc_filter_set`] is `Some`; use `None` only in narrow unit tests
/// with [`ReadFilterParams`] that intentionally do **not** resolve (field-only fallback).
pub fn apply_shard_read_pipeline(
    records: &mut Vec<bam::Record>,
    header: Option<&bam::HeaderView>,
    filters: &ReadFilterParams,
    pipeline: &ShardReadPipelineConfig,
    rng: &mut GatkJavaRng,
) -> Result<(), GatkError> {
    if pipeline.apply_iupac_pre_transform {
        apply_iupac_strict_transform(records);
    }
    records.retain(|rec| match header {
        Some(h) => passes_hc_read_filters_with_header(rec, h, filters),
        None => {
            debug_assert!(
                filters.resolved_hc_filter_set().is_none(),
                "production pipeline must pass BAM header when filters resolve to a Java chain"
            );
            crate::read_model::passes_hc_read_filters(rec, filters)
        }
    });
    if pipeline.apply_dragen_mapq_transform {
        crate::dragen_mq::apply_dragen_mapping_quality_transform(records);
    }
    apply_positional_downsampler(records, header, &pipeline.downsample, rng)
        .map_err(GatkError::generic)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam;

    fn record_with_seq(seq: &[u8]) -> bam::Record {
        use rust_htslib::bam::record::{Cigar, CigarString};
        let qual = vec![30u8; seq.len()];
        let mut rec = bam::Record::new();
        rec.set(
            b"r1",
            Some(&CigarString::from(vec![Cigar::Match(seq.len() as u32)])),
            seq,
            &qual,
        );
        rec
    }

    #[test]
    fn iupac_strict_converts_ambiguity_codes() {
        let mut rec = record_with_seq(b"AMKW");
        apply_iupac_strict_transform(std::slice::from_mut(&mut rec));
        assert_eq!(rec.seq().as_bytes(), b"ANNN");
    }

    #[test]
    fn iupac_strict_leaves_atcg() {
        let mut rec = record_with_seq(b"ACGT");
        apply_iupac_strict_transform(std::slice::from_mut(&mut rec));
        assert_eq!(rec.seq().as_bytes(), b"ACGT");
    }

    #[test]
    fn pipeline_order_filter_before_downsample() {
        use crate::read_model::{ReadFilterParams, FLAG_NOT_PRIMARY};
        let mut records = vec![
            record_with_seq(b"AAAAAAAAAA"),
            record_with_seq(b"AAAAAAAAAA"),
        ];
        for rec in &mut records {
            rec.set_mapq(60);
        }
        records[1].set_flags(FLAG_NOT_PRIMARY);
        let pipeline = ShardReadPipelineConfig {
            apply_iupac_pre_transform: false,
            apply_dragen_mapq_transform: false,
            downsample: PositionalDownsamplerConfig {
                max_reads_per_alignment_start: 1,
                non_random_downsampling_mode: true,
                rng_seed: 0,
            },
        };
        let mut rng = GatkJavaRng::reset_gatk_default();
        apply_shard_read_pipeline(
            &mut records,
            None,
            &ReadFilterParams {
                min_mapping_quality: 20,
                exclude_duplicates: false,
                exclude_secondary: true,
                exclude_supplementary: true,
            },
            &pipeline,
            &mut rng,
        )
        .unwrap();
        assert_eq!(records.len(), 1);
    }
}
