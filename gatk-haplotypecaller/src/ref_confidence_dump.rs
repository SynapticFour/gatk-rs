//! ReferenceConfidenceModel parity dumps.

use crate::assembly_regions_dump::format_activity_prob;
use crate::engine::{CallRegionArgs, HaplotypeCallerEngine};
use crate::genotyping::{build_gvcf_blocks_hc_emit, gvcf_block_to_record_fields, EmitMode};
use crate::gvcf_writer::GATK_HC_DEFAULT_GQB;
use crate::locus_iterator::{IntervalLocusIterator, LocusPileupState};
use crate::minimal_genotyping::calculate_single_sample_ref_vs_any_active_state_profile_value;
use crate::read_model::ReadFilterParams;
use crate::read_transformer::{apply_shard_read_pipeline, ShardReadPipelineConfig};
use crate::ref_confidence::capped_genotype_likelihoods_by_hom_ref;
use crate::ref_confidence::{
    reference_confidence_locus_from_pileup, reference_model_for_no_variation_region,
    ReferenceConfidenceConfig,
};
use crate::reference_vcf_emit::active_region_reference_confidence_loci;
use crate::walker::make_read_shards;
use crate::walker_apply::call_disposition;
use crate::walker_apply::AssemblyRegionCallDisposition;
use crate::walker_traversal::{
    flatten_assembly_regions, traverse_assembly_region_walker, WalkerTraversalConfig,
};
use gatk_common::{GatkError, GatkResult};
use gatk_core::reference::{parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary};
use rust_htslib::bam::Read;
use std::io::Write;
use std::path::Path;

/// H.1.1 — per-locus RCM GL + GQ + DP + activity (C.4 GL walk + reference confidence fields).
pub fn dump_reference_confidence_locus_tsv(
    reference_fasta: &Path,
    bam_path: &Path,
    interval_cli: &str,
    assembly_region_padding: u64,
    out: &mut impl Write,
    read_filters: &ReadFilterParams,
) -> GatkResult<()> {
    let config = ReferenceConfidenceConfig::default();
    let scoring = &config.scoring;
    let ploidy = scoring.sample_ploidy.as_u32();
    let mut header_line = String::from("contig\tpos");
    for i in 0..=ploidy {
        header_line.push_str(&format!("\tgl{i}"));
    }
    header_line.push_str("\tgq\tdp\tactive_prob");
    writeln!(out, "{header_line}").map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;

    let dict = SequenceDictionary::from_fasta_path(reference_fasta)?;
    let specs = parse_intervals_cli_string(&dict, interval_cli)?;
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
            let gl = capped_genotype_likelihoods_by_hom_ref(&pile, &config);
            let active_prob =
                calculate_single_sample_ref_vs_any_active_state_profile_value(&gl, scoring);
            let detail = reference_confidence_locus_from_pileup(pos1 as usize, &pile, &config);
            let mut row = format!("{c}\t{pos1}");
            for g in &gl {
                row.push('\t');
                row.push_str(&format_activity_prob(*g));
            }
            row.push_str(&format!("\t{}\t{}", detail.locus.gq, detail.locus.dp));
            row.push('\t');
            row.push_str(&format_activity_prob(active_prob));
            writeln!(out, "{row}").map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
        }
    }
    Ok(())
}

const DEFAULT_STAND_EMIT_CONFIDENCE: f64 = 10.0;

/// Active `call_region` hybrid RCM loci (Java `getPileupsOverReference` path).
pub fn dump_call_region_active_rcm_loci_tsv(
    ref_fasta: &Path,
    bam_path: &Path,
    interval_cli: &str,
    padding: u64,
    out: &mut impl Write,
) -> GatkResult<()> {
    writeln!(out, "contig\tpos\tgq\tdp").map_err(|e| GatkError::generic(format!("write: {e}")))?;
    let dict = SequenceDictionary::from_fasta_path(ref_fasta)?;
    let specs = parse_intervals_cli_string(&dict, interval_cli)?;
    let filters = ReadFilterParams::gatk_standard_hc();
    let cfg = WalkerTraversalConfig::gatk_haplotype_caller_production(padding);
    let walk = traverse_assembly_region_walker(&dict, &specs, ref_fasta, bam_path, &filters, &cfg)?;
    let regions = flatten_assembly_regions(&walk);
    let args = CallRegionArgs::default();
    let config = ReferenceConfidenceConfig::default();
    let header = rust_htslib::bam::Reader::from_path(bam_path)
        .map_err(|e| GatkError::generic(format!("open bam: {e}")))?
        .header()
        .clone();
    let mut ref_cache = ReferenceWindowCache::new(ref_fasta.to_path_buf(), 4);

    for region in regions.iter().filter(|r| {
        matches!(
            call_disposition(r),
            AssemblyRegionCallDisposition::ActiveFull
        )
    }) {
        let Some(outcome) = HaplotypeCallerEngine::call_region(region, &dict, ref_fasta, &args)?
        else {
            continue;
        };
        let loci = active_region_reference_confidence_loci(
            region,
            &outcome,
            DEFAULT_STAND_EMIT_CONFIDENCE,
            &header,
            &dict,
            &mut ref_cache,
            &filters,
            &config,
            crate::ref_confidence::ClusterRcmEvidenceMode::Production,
        )?;
        writeln!(
            out,
            "region\t{}\t{}\t{}\t{}",
            region.contig,
            region.start.get(),
            region.end.get(),
            loci.len()
        )
        .map_err(|e| GatkError::generic(format!("write: {e}")))?;
        for locus in &loci {
            writeln!(
                out,
                "{}\t{}\t{}\t{}",
                region.contig, locus.position_1based, locus.gq, locus.dp
            )
            .map_err(|e| GatkError::generic(format!("write: {e}")))?;
        }
        let blocks = build_gvcf_blocks_hc_emit(&loci, GATK_HC_DEFAULT_GQB)?;
        writeln!(
            out,
            "block_header\tstart\tend\tmin_dp\tmax_dp\tgq_band\tmin_rgq"
        )
        .map_err(|e| GatkError::generic(format!("write: {e}")))?;
        for block in &blocks {
            let fields = gvcf_block_to_record_fields(block)?;
            writeln!(
                out,
                "block\t{}\t{}\t{}\t{}\t{}\t{}",
                fields.start_1based,
                fields.end_info,
                fields.min_dp,
                fields.max_dp,
                fields.gq_band_upper,
                fields.min_rgq,
            )
            .map_err(|e| GatkError::generic(format!("write: {e}")))?;
        }
    }
    Ok(())
}

/// H.1.2 — first inactive region `referenceModelForNoVariation` summary.
pub fn dump_inactive_reference_model_tsv(
    ref_fasta: &Path,
    bam_path: &Path,
    interval_cli: &str,
    padding: u64,
    emit_mode: EmitMode,
    out: &mut impl Write,
) -> GatkResult<()> {
    let dict = SequenceDictionary::from_fasta_path(ref_fasta)?;
    let specs = parse_intervals_cli_string(&dict, interval_cli)?;
    let filters = ReadFilterParams::gatk_standard_hc();
    let cfg = WalkerTraversalConfig::gatk_haplotype_caller_production(padding);
    let walk = traverse_assembly_region_walker(&dict, &specs, ref_fasta, bam_path, &filters, &cfg)?;
    let regions = flatten_assembly_regions(&walk);
    let region = regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::InactiveReferenceFastPath
            )
        })
        .ok_or_else(|| GatkError::argument("no inactive assembly region in interval"))?;

    let shard = crate::walker::make_read_shards(&dict, &specs, padding)?
        .into_iter()
        .find(|s| s.contig == region.contig)
        .ok_or_else(|| GatkError::argument("no shard for inactive region contig"))?;
    let (header, _) =
        crate::assembly_region_iterator::load_records_for_shard_raw(bam_path, &shard)?;

    let mut ref_cache = ReferenceWindowCache::new(ref_fasta.to_path_buf(), 4);
    let config = ReferenceConfidenceConfig::default();
    let outcome = reference_model_for_no_variation_region(
        region,
        &header,
        &config,
        &filters,
        &mut ref_cache,
        &dict,
        emit_mode,
    )?;

    writeln!(out, "region_contig\t{}", outcome.region_contig)?;
    writeln!(out, "region_start\t{}", outcome.region_start)?;
    writeln!(out, "region_end\t{}", outcome.region_end)?;
    writeln!(out, "is_active\tfalse")?;
    writeln!(out, "path\treferenceModelForNoVariation")?;
    let emit_label = match outcome.emit_mode {
        EmitMode::Gvcf => "GVCF",
        EmitMode::Vcf => "VCF",
        EmitMode::BpResolution => "BP_RESOLUTION",
    };
    writeln!(out, "emit_mode\t{emit_label}")?;
    writeln!(out, "locus_count\t{}", outcome.loci.len())?;
    writeln!(
        out,
        "reference_blocks_emitted\t{}",
        outcome.summary.reference_blocks_emitted
    )?;
    writeln!(
        out,
        "reference_sites_emitted\t{}",
        outcome.summary.reference_sites_emitted
    )?;
    if let Some(first) = outcome.loci.first() {
        writeln!(out, "first_gq\t{}", first.gq)?;
        writeln!(out, "first_dp\t{}", first.dp)?;
    }
    if let Some(last) = outcome.loci.last() {
        writeln!(out, "last_gq\t{}", last.gq)?;
        writeln!(out, "last_dp\t{}", last.dp)?;
    }
    if outcome.emit_mode == EmitMode::Gvcf && !outcome.loci.is_empty() {
        let blocks = build_gvcf_blocks_hc_emit(&outcome.loci, GATK_HC_DEFAULT_GQB)?;
        writeln!(out, "gvcf_block_count\t{}", blocks.len())?;
        writeln!(out, "locus_header\tpos\tgq\tdp")?;
        for locus in &outcome.loci {
            writeln!(
                out,
                "locus\t{}\t{}\t{}",
                locus.position_1based, locus.gq, locus.dp
            )?;
        }
        writeln!(
            out,
            "block_header\tstart\tend\tmin_dp\tmax_dp\tgq_band\tmin_rgq"
        )?;
        for block in &blocks {
            let fields = gvcf_block_to_record_fields(block)?;
            writeln!(
                out,
                "block\t{}\t{}\t{}\t{}\t{}\t{}",
                fields.start_1based,
                fields.end_info,
                fields.min_dp,
                fields.max_dp,
                fields.gq_band_upper,
                fields.min_rgq,
            )?;
        }
    }
    Ok(())
}
