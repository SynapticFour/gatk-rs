//! HQ soft-clip mean dump.
//! Emits the same **RCM** running mean as `ReferenceConfidenceModel#calcGenotypeLikelihoodsOfRefVsAny`,
//! not a read-level heuristic.

use crate::activity_scoring::{
    hq_soft_clip_running_mean_rcm_path, HaplotypeCallerActivityScoringParams,
};
use crate::assembly_regions_dump::load_contig_records_linear;
use crate::locus_iterator::LocusPileupWalker;
use crate::read_model::ReadFilterParams;
use gatk_common::{GatkError, GatkResult};
use gatk_core::reference::{parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary};
use std::io::Write;
use std::path::Path;

/// TSV: `contig`, `pos`, `hq_soft_clip_mean` (8 decimal places, `0` when none).
pub fn dump_hq_soft_clip_mean_tsv(
    reference_fasta: &Path,
    bam_path: &Path,
    intervals_cli: &str,
    out: &mut impl Write,
    read_filters: &ReadFilterParams,
) -> GatkResult<()> {
    writeln!(out, "contig\tpos\thq_soft_clip_mean")
        .map_err(|e| GatkError::generic(format!("write: {e}")))?;

    let dict = SequenceDictionary::from_fasta_path(reference_fasta)?;
    let specs = parse_intervals_cli_string(&dict, intervals_cli)?;
    let mut records_by_contig = std::collections::HashMap::new();
    let mut ref_cache = ReferenceWindowCache::new(reference_fasta.to_path_buf(), 4);
    let mut rng = crate::read_downsample::GatkJavaRng::reset_gatk_default();
    let scoring = HaplotypeCallerActivityScoringParams::default();

    for spec in &specs {
        let (c, s, e) = spec
            .resolve_closed_ends(&dict)
            .map_err(|e| GatkError::argument(e.to_string()))?;
        let ref_window = ref_cache
            .get_interval_bytes(&dict, &c, s, e)
            .map_err(|e| GatkError::generic(e.to_string()))?;
        if !records_by_contig.contains_key(&c) {
            records_by_contig.insert(
                c.clone(),
                load_contig_records_linear(bam_path, &c, read_filters, &mut rng)?,
            );
        }
        let (header, records) = records_by_contig.get(&c).expect("cache");
        let mut walker = LocusPileupWalker::new(records, header, &c, read_filters);
        for pos1 in s..=e {
            let ref_base = *ref_window
                .get((pos1 - s) as usize)
                .ok_or_else(|| GatkError::argument("reference window index out of range"))?;
            let pile = walker.pileup_at(pos1, ref_base)?;
            let mean = hq_soft_clip_running_mean_rcm_path(&pile, &scoring);
            let formatted = if mean == 0.0 {
                "0".to_string()
            } else {
                format!("{mean:.8}")
            };
            writeln!(out, "{c}\t{pos1}\t{formatted}")
                .map_err(|e| GatkError::generic(format!("write: {e}")))?;
        }
    }
    Ok(())
}
