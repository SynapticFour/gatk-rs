//! Dense smoothed activity profile export for parity with GATK `HaplotypeCaller --assembly-region-out`.
//! Walks every 1-based position in `-L`, builds per-locus pileups, feeds [`BandPassActivityProfile`],
//! then writes `position\\tactive_prob` (smoothed) for comparison to the Java IGV derivative.

use crate::activity_profile::{
    ActivityProfileState, ActivityProfileStateKind, BandPassActivityProfile,
    BandPassActivityProfileParams,
};
use crate::activity_scoring::{
    calc_ref_vs_any_log10_genotype_likelihoods, HaplotypeCallerActivityScoringParams,
};
use crate::minimal_genotyping::{
    calculate_single_sample_ref_vs_any_active_state_profile_value,
    cap_genotype_likelihoods_by_hom_ref,
};
use crate::read_header_semantics::ReadHeaderSemantics;
use crate::read_transformer::{apply_shard_read_pipeline, ShardReadPipelineConfig};

use crate::assembly_region_evaluator::{
    add_locus_for_smoothed_activity, hc_activity_after_locus_advance,
};
use crate::assembly_region_iterator::load_records_for_shard_raw;
use crate::feature_context::FeatureDataSources;
use crate::locus_iterator::{IntervalLocusIterator, LocusPileupState};
use crate::read_binding::record_overlaps_closed_interval_1based;
use crate::read_model::ReadFilterParams;
use crate::walker::{make_read_shards, ReadShard, GATK_DEFAULT_ASSEMBLY_REGION_PADDING};
use gatk_common::{GatkError, GatkResult};
use gatk_core::reference::{
    parse_intervals_cli_string, IntervalSpec, ReferenceWindowCache, SequenceDictionary,
};
use rust_htslib::bam;
use std::fs::File;
use std::io::Write;
use std::path::Path;

/// GATK `--alleles` feature-input name in [`FeatureDataSources::load_vcf_source`].
pub const GATK_HC_ALLELES_FEATURE_SOURCE: &str = "alleles";

/// Optional `--force-calling-alleles-present` input for raw-activity parity.
/// # Invariants
/// Borrows [`FeatureDataSources`] for the dump lifetime; does not own allele resources.
/// # Ownership
/// Lifetime-bound borrow of feature sources plus force-call-filtered flag.
/// # Mutation
/// Immutable dump input.
/// # Biological assumptions
/// Force-calling alleles can mark loci active even without strong pileup evidence.
/// # Java equivalence
/// GATK `HaplotypeCallerArgumentCollection.forceCallFiltered` + alleles feature source.
#[derive(Clone, Copy)]
pub struct ForceCallingAllelesDump<'a> {
    pub sources: &'a FeatureDataSources,
    /// `HaplotypeCallerArgumentCollection.forceCallFiltered`
    pub force_call_filtered: bool,
}

/// Padded-shard read set + LIBS pileup cursor (Java `HcFullParityGateDump#makeLocusIterator` / `AssemblyRegionIterator`).
struct ShardActivityWalk {
    header: bam::HeaderView,
    records: Vec<bam::Record>,
    semantics: ReadHeaderSemantics,
    pileup_state: LocusPileupState,
}

fn open_shard_activity_walk(
    bam_path: &Path,
    shard: &ReadShard,
    read_filters: &ReadFilterParams,
    pipeline: &ShardReadPipelineConfig,
    rng: &mut crate::read_downsample::GatkJavaRng,
) -> GatkResult<ShardActivityWalk> {
    let (header, mut records) = match load_records_for_shard_raw(bam_path, shard) {
        Ok(pair) => pair,
        Err(_) => {
            let (header, all_records) =
                load_contig_records_linear(bam_path, &shard.contig, read_filters, rng)?;
            let filtered: Vec<bam::Record> = all_records
                .into_iter()
                .filter(|r| {
                    shard.padded_spans.iter().any(|&(rs, re)| {
                        record_overlaps_closed_interval_1based(
                            r,
                            &header,
                            &shard.contig,
                            rs,
                            re,
                            read_filters,
                        )
                    })
                })
                .collect();
            (header, filtered)
        }
    };
    apply_shard_read_pipeline(&mut records, Some(&header), read_filters, pipeline, rng)?;
    let semantics = ReadHeaderSemantics::from_bam_header_view(&header)?;
    let pileup_state =
        LocusPileupState::from_records(&records, &header, &shard.contig, read_filters);
    Ok(ShardActivityWalk {
        header,
        records,
        semantics,
        pileup_state,
    })
}

fn shard_for_spec(
    dict: &SequenceDictionary,
    spec: &IntervalSpec,
    assembly_region_padding: u64,
) -> GatkResult<ReadShard> {
    make_read_shards(dict, std::slice::from_ref(spec), assembly_region_padding)?
        .into_iter()
        .next()
        .ok_or_else(|| GatkError::argument("no read shard for interval"))
}

fn force_active_at_locus(
    force_calling: Option<&ForceCallingAllelesDump<'_>>,
    contig: &str,
    pos1: u64,
) -> bool {
    force_calling.is_some_and(|fc| {
        fc.sources.force_calling_allele_overlaps_locus(
            GATK_HC_ALLELES_FEATURE_SOURCE,
            contig,
            pos1,
            fc.force_call_filtered,
        )
    })
}

/// Fixed-width probability formatting for L1 parity TSVs.
pub fn format_activity_prob(v: f64) -> String {
    if v == 0.0 {
        return "0".to_string();
    }
    format!("{:.8}", v)
}

fn format_activity_kind(kind: ActivityProfileStateKind) -> &'static str {
    match kind {
        ActivityProfileStateKind::Plain => "none",
        ActivityProfileStateKind::HighQualitySoftClips { .. } => "hq_soft_clips",
    }
}

fn write_raw_activity_row(out: &mut impl Write, st: &ActivityProfileState) -> GatkResult<()> {
    // Match pinned GATK `HcFullParityGateDump` raw-activity (both columns use `isActiveProb`).
    let prob = format_activity_prob(st.active_prob);
    writeln!(
        out,
        "{}\t{}\t{}\t{}\t{}",
        st.contig,
        st.pos,
        &prob,
        &prob,
        format_activity_kind(st.kind()),
    )
    .map_err(|e| GatkError::generic(format!("write tsv: {e}")))
}

/// Write one line per 1-based position in the given intervals: `contig\\tpos\\tsmoothed_active_prob`.
pub fn dump_smoothed_activity_tsv(
    reference_fasta: &Path,
    bam_path: &Path,
    intervals_cli: &str,
    out_tsv: &Path,
    read_filters: &ReadFilterParams,
) -> GatkResult<()> {
    let mut out = File::create(out_tsv)
        .map_err(|e| GatkError::generic(format!("create {}: {e}", out_tsv.display())))?;
    dump_smoothed_activity_profile_tsv(
        reference_fasta,
        bam_path,
        intervals_cli,
        GATK_DEFAULT_ASSEMBLY_REGION_PADDING,
        &mut out,
        read_filters,
    )
}

/// Band-pass **smoothed** activity per locus.
pub fn dump_smoothed_activity_profile_tsv(
    reference_fasta: &Path,
    bam_path: &Path,
    intervals_cli: &str,
    assembly_region_padding: u64,
    out: &mut impl Write,
    read_filters: &ReadFilterParams,
) -> GatkResult<()> {
    writeln!(out, "contig\tpos\tsmoothed_active_prob")
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    dump_smoothed_profile_rows(
        reference_fasta,
        bam_path,
        intervals_cli,
        assembly_region_padding,
        read_filters,
        |contig, pos, st, threshold| {
            writeln!(
                out,
                "{contig}\t{pos}\t{}",
                format_activity_prob(st.active_prob)
            )
            .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
            let _ = threshold;
            Ok(())
        },
    )
}

/// Per-locus binary `isActive` after smoothing (`active_prob > threshold`).
pub fn dump_active_locus_tsv(
    reference_fasta: &Path,
    bam_path: &Path,
    intervals_cli: &str,
    assembly_region_padding: u64,
    out: &mut impl Write,
    read_filters: &ReadFilterParams,
) -> GatkResult<()> {
    writeln!(out, "contig\tpos\tis_active")
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    dump_smoothed_profile_rows(
        reference_fasta,
        bam_path,
        intervals_cli,
        assembly_region_padding,
        read_filters,
        |contig, pos, st, threshold| {
            let fudge = std::env::var("PARITY_HC_ACTIVE_THRESHOLD_FUDGE")
                .ok()
                .and_then(|s| s.parse::<f64>().ok())
                .unwrap_or(0.0);
            let active = st.active_prob > threshold - fudge;
            writeln!(
                out,
                "{contig}\t{pos}\t{}",
                if active { "true" } else { "false" }
            )
            .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
            Ok(())
        },
    )
}

fn dump_smoothed_profile_rows(
    reference_fasta: &Path,
    bam_path: &Path,
    intervals_cli: &str,
    assembly_region_padding: u64,
    read_filters: &ReadFilterParams,
    mut write_row: impl FnMut(&str, u64, &ActivityProfileState, f64) -> GatkResult<()>,
) -> GatkResult<()> {
    let dict = SequenceDictionary::from_fasta_path(reference_fasta)?;
    let specs = parse_intervals_cli_string(&dict, intervals_cli)?;
    if specs.is_empty() {
        return Err(GatkError::argument("no intervals after parse"));
    }

    let scoring = HaplotypeCallerActivityScoringParams::default();
    let bp_params = BandPassActivityProfileParams::gatk_haplotype_caller_defaults();
    let threshold = bp_params.active_prob_threshold;
    let pipeline = ShardReadPipelineConfig::gatk_haplotype_caller_production();
    let mut ref_cache = ReferenceWindowCache::new(reference_fasta.to_path_buf(), 4);
    let mut rng = crate::read_downsample::GatkJavaRng::reset_gatk_default();

    for spec in &specs {
        let (c, s, e) = spec
            .resolve_closed_ends(&dict)
            .map_err(|e| GatkError::argument(e.to_string()))?;
        let ref_window = ref_cache
            .get_interval_bytes(&dict, &c, s, e)
            .map_err(|e| GatkError::generic(e.to_string()))?;
        let contig_len = dict
            .contig(&c)
            .ok_or_else(|| GatkError::argument(format!("unknown contig {c}")))?
            .length;
        let shard = shard_for_spec(&dict, spec, assembly_region_padding)?;
        let mut walk =
            open_shard_activity_walk(bam_path, &shard, read_filters, &pipeline, &mut rng)?;

        let mut prof = BandPassActivityProfile::new(c.as_str(), contig_len, bp_params.clone());
        for pos1 in IntervalLocusIterator::from_closed_interval(s, e) {
            let ref_base = *ref_window
                .get((pos1 - s) as usize)
                .ok_or_else(|| GatkError::argument("reference window index out of range"))?;
            add_locus_for_smoothed_activity(
                &mut prof,
                &walk.records,
                &walk.header,
                &walk.semantics,
                &c,
                pos1,
                read_filters,
                ref_base,
                &scoring,
                Some(&mut walk.pileup_state),
                false,
            )?;
        }

        let rs = prof
            .region_start()
            .ok_or_else(|| GatkError::generic("activity profile empty"))?;
        for (i, st) in prof.states().iter().enumerate() {
            let pos = rs + i as u64;
            if pos < s || pos > e {
                continue;
            }
            write_row(&c, pos, st, threshold)?;
        }
    }

    Ok(())
}

/// Per-locus **raw** `ActivityProfileState` (no band-pass), one row per 1-based position in `-L`.
/// Schema: `contig`, `pos`, `active_prob`, `original_active_prob`, `kind` (`none` | `hq_soft_clips`).
/// Parity gate (`hc_full_parity_gate_dump raw-activity`).
pub fn dump_raw_activity_profile_tsv(
    reference_fasta: &Path,
    bam_path: &Path,
    intervals_cli: &str,
    out: &mut impl Write,
    read_filters: &ReadFilterParams,
) -> GatkResult<()> {
    dump_raw_activity_profile_tsv_with_force_calling(
        reference_fasta,
        bam_path,
        intervals_cli,
        GATK_DEFAULT_ASSEMBLY_REGION_PADDING,
        out,
        read_filters,
        None,
    )
}

/// [`dump_raw_activity_profile_tsv`] with optional `--force-calling-alleles-present` / `--alleles` features.
pub fn dump_raw_activity_profile_tsv_with_contamination(
    reference_fasta: &Path,
    bam_path: &Path,
    intervals_cli: &str,
    contamination_fraction: f64,
    out: &mut impl Write,
    read_filters: &ReadFilterParams,
) -> GatkResult<()> {
    let mut scoring = HaplotypeCallerActivityScoringParams::default();
    scoring.contamination_fraction_to_filter = contamination_fraction;
    writeln!(out, "contig\tpos\tactive_prob\toriginal_active_prob\tkind")
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    dump_raw_activity_rows_with_scoring(
        reference_fasta,
        bam_path,
        intervals_cli,
        GATK_DEFAULT_ASSEMBLY_REGION_PADDING,
        read_filters,
        &scoring,
        None,
        out,
    )
}

pub fn dump_raw_activity_profile_tsv_with_force_calling(
    reference_fasta: &Path,
    bam_path: &Path,
    intervals_cli: &str,
    assembly_region_padding: u64,
    out: &mut impl Write,
    read_filters: &ReadFilterParams,
    force_calling: Option<ForceCallingAllelesDump<'_>>,
) -> GatkResult<()> {
    writeln!(out, "contig\tpos\tactive_prob\toriginal_active_prob\tkind")
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    dump_raw_activity_rows_with_scoring(
        reference_fasta,
        bam_path,
        intervals_cli,
        assembly_region_padding,
        read_filters,
        &HaplotypeCallerActivityScoringParams::default(),
        force_calling,
        out,
    )
}

fn dump_raw_activity_rows_with_scoring(
    reference_fasta: &Path,
    bam_path: &Path,
    intervals_cli: &str,
    assembly_region_padding: u64,
    read_filters: &ReadFilterParams,
    scoring: &HaplotypeCallerActivityScoringParams,
    force_calling: Option<ForceCallingAllelesDump<'_>>,
    out: &mut impl Write,
) -> GatkResult<()> {
    let dict = SequenceDictionary::from_fasta_path(reference_fasta)?;
    let specs = parse_intervals_cli_string(&dict, intervals_cli)?;
    if specs.is_empty() {
        return Err(GatkError::argument("no intervals after parse"));
    }

    let pipeline = ShardReadPipelineConfig::gatk_haplotype_caller_production();
    let mut ref_cache = ReferenceWindowCache::new(reference_fasta.to_path_buf(), 4);
    let mut rng = crate::read_downsample::GatkJavaRng::reset_gatk_default();

    for spec in &specs {
        let (c, s, e) = spec
            .resolve_closed_ends(&dict)
            .map_err(|e| GatkError::argument(e.to_string()))?;
        let ref_window = ref_cache
            .get_interval_bytes(&dict, &c, s, e)
            .map_err(|e| GatkError::generic(e.to_string()))?;
        let shard = shard_for_spec(&dict, spec, assembly_region_padding)?;
        let mut walk =
            open_shard_activity_walk(bam_path, &shard, read_filters, &pipeline, &mut rng)?;

        for pos1 in IntervalLocusIterator::from_closed_interval(s, e) {
            let ref_base = *ref_window
                .get((pos1 - s) as usize)
                .ok_or_else(|| GatkError::argument("reference window index out of range"))?;
            walk.pileup_state
                .advance_to(&walk.records, read_filters, pos1)?;
            let st = hc_activity_after_locus_advance(
                &c,
                pos1,
                &mut walk.pileup_state,
                &walk.records,
                &walk.semantics,
                scoring,
                ref_base,
                force_active_at_locus(force_calling.as_ref(), &c, pos1),
            )?;
            write_raw_activity_row(out, &st)?;
        }
    }

    Ok(())
}

/// Per-locus log10 genotype likelihoods + activity (/ `c4-gl`).
/// Schema: `contig`, `pos`, `gl0`…`gl{ploidy}`, `active_prob`, `original_active_prob`.
/// Matches Java `HcFullParityGateDump#walkGenotypeLikelihood`: one padded shard per user interval,
/// production read pipeline, capped GL vector for output and activity.
pub fn dump_genotype_likelihood_activity_tsv(
    reference_fasta: &Path,
    bam_path: &Path,
    intervals_cli: &str,
    assembly_region_padding: u64,
    out: &mut impl Write,
    read_filters: &ReadFilterParams,
) -> GatkResult<()> {
    dump_genotype_likelihood_activity_tsv_inner(
        reference_fasta,
        bam_path,
        intervals_cli,
        assembly_region_padding,
        out,
        read_filters,
    )
}

pub fn dump_genotype_likelihood_activity_tsv_default_padding(
    reference_fasta: &Path,
    bam_path: &Path,
    intervals_cli: &str,
    out: &mut impl Write,
    read_filters: &ReadFilterParams,
) -> GatkResult<()> {
    dump_genotype_likelihood_activity_tsv(
        reference_fasta,
        bam_path,
        intervals_cli,
        GATK_DEFAULT_ASSEMBLY_REGION_PADDING,
        out,
        read_filters,
    )
}

fn dump_genotype_likelihood_activity_tsv_inner(
    reference_fasta: &Path,
    bam_path: &Path,
    intervals_cli: &str,
    assembly_region_padding: u64,
    out: &mut impl Write,
    read_filters: &ReadFilterParams,
) -> GatkResult<()> {
    let scoring = HaplotypeCallerActivityScoringParams::default();
    let ploidy = scoring.sample_ploidy.as_u32();
    let mut header_line = String::from("contig\tpos");
    for i in 0..=ploidy {
        header_line.push_str(&format!("\tgl{i}"));
    }
    header_line.push_str("\tactive_prob\toriginal_active_prob");
    writeln!(out, "{header_line}").map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;

    let dict = SequenceDictionary::from_fasta_path(reference_fasta)?;
    let specs = parse_intervals_cli_string(&dict, intervals_cli)?;
    if specs.is_empty() {
        return Err(GatkError::argument("no intervals after parse"));
    }

    let mut ref_cache = ReferenceWindowCache::new(reference_fasta.to_path_buf(), 4);
    let mut rng = crate::read_downsample::GatkJavaRng::reset_gatk_default();
    let pipeline = ShardReadPipelineConfig::gatk_haplotype_caller_production();

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
            let pile = pileup_state.pileup_at(&records, read_filters, pos1, ref_base)?;
            let gl_raw = calc_ref_vs_any_log10_genotype_likelihoods(ploidy, &pile, &scoring);
            let gl_dump = cap_genotype_likelihoods_by_hom_ref(&gl_raw);
            // Java `genotypeLikelihoodActivity`: capped raw GL columns; activity on capped GL.
            let active_prob =
                calculate_single_sample_ref_vs_any_active_state_profile_value(&gl_dump, &scoring);
            let max_gl = gl_dump.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            let original = max_gl - gl_dump[0];
            let mut row = format!("{c}\t{pos1}");
            for g in &gl_dump {
                row.push('\t');
                row.push_str(&format_activity_prob(*g));
            }
            row.push('\t');
            row.push_str(&format_activity_prob(active_prob));
            row.push('\t');
            row.push_str(&format_activity_prob(original));
            writeln!(out, "{row}").map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
        }
    }

    Ok(())
}

/// Per-locus LIBS pileup depth (`AlignmentContext.size`) / `b2-locus`.
/// Schema: `contig`, `pos`, `pileup_depth`. Walks every 1-based position in user `-L` intervals
/// (grouped per contig like GATK `MultiIntervalLocalReadShard`), with reads bounded to padded spans.
pub fn dump_locus_pileup_tsv(
    reference_fasta: &Path,
    bam_path: &Path,
    intervals_cli: &str,
    assembly_region_padding: u64,
    out: &mut impl Write,
    read_filters: &ReadFilterParams,
) -> GatkResult<()> {
    writeln!(out, "contig\tpos\tpileup_depth")
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;

    let dict = SequenceDictionary::from_fasta_path(reference_fasta)?;
    let specs = parse_intervals_cli_string(&dict, intervals_cli)?;
    if specs.is_empty() {
        return Err(GatkError::argument("no intervals after parse"));
    }

    let shards = make_read_shards(&dict, &specs, assembly_region_padding)?;
    let mut ref_cache = ReferenceWindowCache::new(reference_fasta.to_path_buf(), 4);
    let mut rng = crate::read_downsample::GatkJavaRng::reset_gatk_default();

    for shard in &shards {
        let (header, all_records) = crate::read_transformer::load_contig_records_hc_production(
            bam_path,
            &shard.contig,
            read_filters,
            &mut rng,
        )?;
        let filtered: Vec<bam::Record> = all_records
            .into_iter()
            .filter(|r| {
                shard.padded_spans.iter().any(|&(rs, re)| {
                    record_overlaps_closed_interval_1based(
                        r,
                        &header,
                        &shard.contig,
                        rs,
                        re,
                        read_filters,
                    )
                })
            })
            .collect();

        let mut pileup_state =
            LocusPileupState::from_records(&filtered, &header, &shard.contig, read_filters);

        for &(s, e) in &shard.user_spans {
            if pileup_state.last_pos1.is_some_and(|prev| s < prev) {
                pileup_state.reset_cursor();
            }
            let ref_window = ref_cache
                .get_interval_bytes(&dict, &shard.contig, s, e)
                .map_err(|e| GatkError::generic(e.to_string()))?;
            for pos1 in IntervalLocusIterator::from_closed_interval(s, e) {
                let ref_base = *ref_window
                    .get((pos1 - s) as usize)
                    .ok_or_else(|| GatkError::argument("reference window index out of range"))?;
                let depth = pileup_state
                    .pileup_at(&filtered, read_filters, pos1, ref_base)?
                    .len();
                writeln!(out, "{}\t{}\t{}", shard.contig, pos1, depth)
                    .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
            }
        }
    }

    Ok(())
}

/// [`dump_locus_pileup_tsv`] with GATK default assembly-region padding.
pub fn dump_locus_pileup_tsv_default_padding(
    reference_fasta: &Path,
    bam_path: &Path,
    intervals_cli: &str,
    out: &mut impl Write,
    read_filters: &ReadFilterParams,
) -> GatkResult<()> {
    dump_locus_pileup_tsv(
        reference_fasta,
        bam_path,
        intervals_cli,
        GATK_DEFAULT_ASSEMBLY_REGION_PADDING,
        out,
        read_filters,
    )
}

pub(crate) fn load_contig_records_linear(
    bam_path: &Path,
    contig: &str,
    filters: &ReadFilterParams,
    rng: &mut crate::read_downsample::GatkJavaRng,
) -> GatkResult<(bam::HeaderView, Vec<bam::Record>)> {
    crate::read_transformer::load_contig_records_hc_production(bam_path, contig, filters, rng)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::read_model::ReadFilterParams;

    #[test]
    fn format_activity_prob_is_stable() {
        assert_eq!(format_activity_prob(0.0), "0");
        assert_eq!(format_activity_prob(0.123456789), "0.12345679");
    }

    /// P12 cluster: per-locus activity must match genotype-likelihood dump (not joint multisample 0.76).
    #[test]
    fn p12_boundary_activity_matches_genotype_likelihood_path() {
        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../parity/realworld");
        let ref_fa = root.join("assets/hs37d5.simple.fa");
        let bam = root.join("na12878_20k_b37/NA12878_20k.b37.bam");
        if !bam.exists() {
            return;
        }
        let dict = SequenceDictionary::from_fasta_path(&ref_fa).unwrap();
        let spec = parse_intervals_cli_string(&dict, "2:92307220-92307220")
            .unwrap()
            .remove(0);
        let scoring = HaplotypeCallerActivityScoringParams::default();
        let filters = ReadFilterParams::gatk_standard_hc();
        let shard = shard_for_spec(&dict, &spec, GATK_DEFAULT_ASSEMBLY_REGION_PADDING).unwrap();
        let mut rng = crate::read_downsample::GatkJavaRng::reset_gatk_default();
        let pipeline = ShardReadPipelineConfig::gatk_haplotype_caller_production();
        let mut walk =
            open_shard_activity_walk(&bam, &shard, &filters, &pipeline, &mut rng).unwrap();
        let (c, s, e) = spec.resolve_closed_ends(&dict).unwrap();
        let mut ref_cache = ReferenceWindowCache::new(ref_fa.clone(), 4);
        let ref_base = ref_cache.get_interval_bytes(&dict, &c, s, e).unwrap()[0];
        walk.pileup_state
            .advance_to(&walk.records, &filters, s)
            .unwrap();
        let st = hc_activity_after_locus_advance(
            &c,
            s,
            &mut walk.pileup_state,
            &walk.records,
            &walk.semantics,
            &scoring,
            ref_base,
            false,
        )
        .unwrap();
        let pile = walk
            .pileup_state
            .pileup_observations(&walk.records, ref_base)
            .unwrap();
        let gl_raw = calc_ref_vs_any_log10_genotype_likelihoods(
            scoring.sample_ploidy.as_u32(),
            &pile,
            &scoring,
        );
        let gl = cap_genotype_likelihoods_by_hom_ref(&gl_raw);
        let active = calculate_single_sample_ref_vs_any_active_state_profile_value(&gl, &scoring);
        assert_eq!(
            active,
            st.active_prob,
            "gl_raw={gl_raw:?} gl={gl:?} unique_samples={} single_header={}",
            walk.semantics.unique_sample_count(),
            walk.semantics.is_single_sample_header()
        );
        assert!(
            st.active_prob < 0.01,
            "inactive prefix: st.active_prob={} gl_raw={gl_raw:?} gl_capped={gl:?}",
            st.active_prob
        );
    }
}
