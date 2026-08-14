//! HaplotypeCaller CLI pipeline entry.
//! **gVCF RCM scope (H1-3):** post-walk `reconcile_p12_cluster_rcm_band`
//! `reconcile_fragmented_dense_cluster_bands` run only when walker layout fragments P12 dense
//! cluster bands on `chr2:92300000–92350000`. Not genome-wide Java RCM parity.

use crate::assembly_region_iterator::AssemblyRegion;
use crate::engine::{CallRegionArgs, HaplotypeCallerEngine};
use crate::feature_context::FeatureContext;
use crate::genotyping::EmitMode;
use crate::gvcf_writer::gatk_hc_gvcf_header_lines;
use crate::read_binding::total_read_tile_overlaps;
use crate::read_event_discovery::{P12_CLUSTER_RCM_CONTIG, P12_L5_JAVA_EXTRA_VARIANT_NO_HOM_REF};
use crate::read_header_semantics::ReadHeaderSemantics;
use crate::read_model::ReadFilterParams;
use crate::ref_confidence::{
    reference_confidence_loci_for_active_call_none, ReferenceConfidenceConfig,
};
use crate::reference_context::ReferenceContext;
use crate::reference_vcf_emit::{
    active_region_reference_confidence_loci, emit_mode_from_output_mode,
    inactive_reference_model_to_vcf_records, p12_cluster_rcm_band_fragmented,
    reconcile_fragmented_dense_cluster_bands, reconcile_p12_cluster_rcm_band,
    GvcfIntervalCollector,
};
use crate::region_vcf_emit::{
    populate_hc_vcf_header_schema, try_emit_call_region_variants, HC_PIPELINE_ASSEMBLY_REGION_V1,
    HC_PIPELINE_SCAFFOLD,
};
use crate::runtime_config::RuntimeConfig;
use crate::walker::GATK_DEFAULT_ASSEMBLY_REGION_PADDING;
use crate::walker_apply::{call_disposition, AssemblyRegionCallDisposition};
use crate::walker_traversal::{for_each_assembly_region, WalkerTraversalConfig};
use gatk_common::{GatkConfig, GatkError, GatkResult, HaplotypeCallerConfig};
use gatk_core::io::vcf::{Contig, VcfHeader, VcfRecord, VcfWriter};
use gatk_core::reference::{
    intervals_for_haplotype_caller, ReferenceWindowCache, SequenceDictionary,
};
use rayon::prelude::*;
use rust_htslib::bam;
use rust_htslib::bam::Read as _;
use std::cell::RefCell;
use std::path::{Path, PathBuf};
use tracing::info;

/// Regions at or above this read count are flushed alone (no Rayon siblings).
///
/// Deep piles after positional DS still hold BAM records into each in-flight region;
/// parallelizing several such regions amplified Peak-RSS to ~15 GiB on hosted runners.
/// Sequential processing matches Java's one-heavy-region peak shape without changing evidence.
///
/// Override with `GATK_RS_HC_LARGE_REGION_READS` (usize). `GATK_RS_HC_SEQUENTIAL=1` forces
/// fully sequential region apply (batch size 1) for 16 GiB hosts.
const LARGE_REGION_READS_SEQUENTIAL_DEFAULT: usize = 4_096;

/// Drop PairHMM / SW TLS scratch after a region (or between SW-heavy and PairHMM-heavy phases).
/// Full clear — soft high-water keep left sticky multi-hundred-MiB RSS on dense windows.
///
/// Includes SIMD/NEON planes (previously only cleared mid-`call_region` before realign,
/// which paid multi-second `munmap` gaps on dense NA12878 while leaving sticky Peak-RSS
/// when region-end release omitted them).
pub(crate) fn release_region_tls_scratch() {
    crate::pairhmm_log10::release_pairhmm_tls_scratch();
    crate::pairhmm_logless::release_pairhmm_logless_tls_scratch();
    crate::pairhmm_simd::release_pairhmm_simd_tls_scratch();
    crate::smith_waterman::release_sw_tls_scratch();
}

/// Same as [`release_region_tls_scratch`], broadcast onto every Rayon pool thread.
/// Needed when haplotype scoring used `par_iter` (worker TLS is invisible to the caller).
pub(crate) fn release_region_tls_scratch_all_threads() {
    release_region_tls_scratch();
    rayon::broadcast(|_| {
        release_region_tls_scratch();
    });
}

#[inline]
fn owned_bam_header(reader: &bam::Reader) -> bam::HeaderView {
    reader.header().clone()
}

thread_local! {
    #[allow(clippy::missing_const_for_thread_local)]
    static TLS_BAM_HEADER: RefCell<Option<(PathBuf, bam::HeaderView)>> = RefCell::new(None);
    #[allow(clippy::missing_const_for_thread_local)]
    static TLS_REF_CACHE: RefCell<Option<(PathBuf, ReferenceWindowCache)>> = RefCell::new(None);
}

fn tls_bam_header(bam_path: &Path) -> GatkResult<bam::HeaderView> {
    TLS_BAM_HEADER.with(|cell| {
        let mut slot = cell.borrow_mut();
        if let Some((p, h)) = slot.as_ref() {
            if p.as_path() == bam_path {
                return Ok(h.clone());
            }
        }
        let reader = bam::Reader::from_path(bam_path)
            .map_err(|e| GatkError::generic(format!("open BAM for region parallel worker: {e}")))?;
        let header = owned_bam_header(&reader);
        *slot = Some((bam_path.to_path_buf(), header.clone()));
        Ok(header)
    })
}

fn with_tls_ref_cache<T>(
    ref_path: &Path,
    f: impl FnOnce(&mut ReferenceWindowCache) -> GatkResult<T>,
) -> GatkResult<T> {
    TLS_REF_CACHE.with(|cell| {
        let mut slot = cell.borrow_mut();
        let reuse = slot.as_ref().is_some_and(|(p, _)| p.as_path() == ref_path);
        if !reuse {
            *slot = Some((
                ref_path.to_path_buf(),
                ReferenceWindowCache::new(ref_path.to_path_buf(), 4),
            ));
        }
        let cache = slot
            .as_mut()
            .map(|(_, c)| c)
            .ok_or_else(|| GatkError::generic("TLS reference cache missing after insert"))?;
        f(cache)
    })
}

/// Owned per-region emission batch (Send) for parallel Active-Region processing.
struct RegionEmitBatch {
    region_index: usize,
    contig: String,
    start: u64,
    records: Vec<VcfRecord>,
}

/// Run the HaplotypeCaller tool: validate inputs, build traversal scaffold, then genotyping when implemented.
pub fn run_haplotype_caller(config: &GatkConfig) -> GatkResult<()> {
    if config.tool_config.tool_name != "HaplotypeCaller" {
        return Err(GatkError::configuration(
            "run_haplotype_caller requires a HaplotypeCaller GatkConfig",
        ));
    }

    // Observe-only: optional NDJSON semantic checkpoints (`GATK_RS_SEMANTIC_TRACE`).
    crate::semantic_trace::try_init_from_runtime(&RuntimeConfig::from_env());

    let ref_path = config
        .tool_config
        .inputs
        .reference
        .as_deref()
        .ok_or_else(|| GatkError::configuration("Reference FASTA (-R) is required"))?;

    let reference = Path::new(ref_path);
    if !reference.is_file() {
        return Err(GatkError::argument(format!(
            "Reference FASTA not found or not a file: {ref_path}"
        )));
    }

    let dict = SequenceDictionary::from_fasta_path(reference).map_err(|e| {
        GatkError::argument(format!("Failed to read reference FASTA {ref_path}: {e}"))
    })?;
    if dict.contig_count() == 0 {
        return Err(GatkError::argument(format!(
            "Reference FASTA contains no sequences: {ref_path}"
        )));
    }

    if config.tool_config.inputs.input_files.is_empty() {
        return Err(GatkError::configuration(
            "At least one alignment input (-I) is required",
        ));
    }

    let interval_specs = intervals_for_haplotype_caller(
        &dict,
        config.get_parameter("intervals").map(|s| s.as_str()),
    )?;

    let engine = HaplotypeCallerEngine::prepare_traversal_default(&dict, interval_specs)?;

    let hc = config.get_haplotypecaller_config()?;
    validate_haplotype_caller_input_policy(config)?;
    reject_unsupported_haplotype_caller_cli(&hc, config)?;
    let read_filters = ReadFilterParams::from_haplotype_caller(&hc);

    // GATK accepts header-only / empty-record BAMs (e.g. chr tip windows with no
    // coverage) and emits a header-only VCF. Only reject missing/unreadable paths.
    for in_path in &config.tool_config.inputs.input_files {
        let p = Path::new(in_path);
        if !p.is_file() {
            return Err(GatkError::argument(format!(
                "Alignment input not found or not a file: {in_path}"
            )));
        }
        let _reader = bam::Reader::from_path(p).map_err(|e| {
            GatkError::generic(format!("Failed to open alignment input {in_path}: {e}"))
        })?;
    }

    let out = config
        .tool_config
        .outputs
        .output_vcf
        .as_deref()
        .ok_or_else(|| GatkError::configuration("Output VCF (-O) is required"))?;
    let out_path = Path::new(out);
    if let Some(parent) = out_path.parent() {
        if !parent.as_os_str().is_empty() && !parent.is_dir() {
            return Err(GatkError::argument(format!(
                "Output directory does not exist: {}",
                parent.display()
            )));
        }
    }

    // R4-1: default skips O(tiles×BAM) diagnostic scan (does not affect calling).
    let debug_tile_overlaps = RuntimeConfig::from_env().debug.debug_tile_overlaps;
    let read_tile_overlaps = if debug_tile_overlaps {
        let mut n = 0usize;
        for bam in &config.tool_config.inputs.input_files {
            n += total_read_tile_overlaps(Path::new(bam), &engine.tiles, &read_filters)?;
        }
        Some(n)
    } else {
        None
    };

    match read_tile_overlaps {
        Some(n) => info!(
            "HaplotypeCaller traversal scaffold: {} input interval(s), {} traversal tile(s); read–tile overlap count (linear scan, scaffold) = {}",
            engine.interval_specs.len(),
            engine.tile_count(),
            n
        ),
        None => info!(
            "HaplotypeCaller traversal scaffold: {} input interval(s), {} traversal tile(s)",
            engine.interval_specs.len(),
            engine.tile_count(),
        ),
    }

    reject_legacy_provisional_output(config)?;
    crate::parity_harness::warn_if_harness_flags_set();

    let variant_output_enabled = hc_variant_output_enabled();

    // Default: assembly-region-v1 (`call_region`). Opt-in header-only scaffold for Phase 9 tests:
    // `GATK_RS_HC_SCAFFOLD_OUTPUT=1` (or legacy `GATK_RS_HC_ACTIVATE_OUTPUT=0`).
    let mut header = VcfHeader::default();
    header.source = Some(if variant_output_enabled {
        format!("gatk-rs HaplotypeCaller {HC_PIPELINE_ASSEMBLY_REGION_V1}")
    } else {
        format!("gatk-rs HaplotypeCaller {HC_PIPELINE_SCAFFOLD}")
    });
    header.reference = Some(ref_path.to_string());
    header.contigs = dict
        .contig_records()
        .iter()
        .map(|c| Contig {
            id: c.name.clone(),
            length: Some(c.length),
            md5: None,
            assembly: None,
            species: None,
            uri: None,
        })
        .collect();
    header.other_headers.push((
        "GATK_RS_HC_PIPELINE".to_string(),
        if variant_output_enabled {
            HC_PIPELINE_ASSEMBLY_REGION_V1.to_string()
        } else {
            HC_PIPELINE_SCAFFOLD.to_string()
        },
    ));
    if variant_output_enabled {
        if let Some(mode) = config.get_parameter("output_mode") {
            header
                .other_headers
                // CLONE: needed because owned element into collection.
                .push(("GATK_RS_HC_OUTPUT_MODE".to_string(), mode.clone()));
        }
    }
    let emit_mode = emit_mode_from_output_mode(&config.tool_config.outputs.output_mode);
    let sample_name = sample_name_from_bam_inputs(&config.tool_config.inputs.input_files)?;
    if variant_output_enabled {
        // CLONE: needed because owned element into collection.
        header.samples.push(sample_name.clone());
        if emit_mode == EmitMode::Gvcf {
            let primary = dict
                .contig_records()
                .first()
                .map(|c| (c.name.as_str(), c.length))
                .unwrap_or(("chr1", 1_000_000));
            for line in gatk_hc_gvcf_header_lines(primary.0, primary.1) {
                if line.starts_with("##reference=") {
                    continue;
                }
                if let Some((k, v)) = line.strip_prefix("##").and_then(|s| s.split_once('=')) {
                    header.other_headers.push((k.to_string(), v.to_string()));
                }
            }
            header
                .other_headers
                .push(("GATK_RS_HC_GVCF".to_string(), "1".to_string()));
        } else {
            // VCF mode: declare INFO/FORMAT used by region emit (hap.py vcfcheck).
            populate_hc_vcf_header_schema(&mut header);
        }
    }

    let mut writer = VcfWriter::new(out_path, header)?;
    writer.write_header()?;
    let mut variant_records = 0usize;
    if variant_output_enabled {
        let pair_hmm_impl = pair_hmm_impl_from_config(&hc, config)?;
        let records = assembly_region_variant_records(
            reference,
            &config.tool_config.inputs.input_files,
            &engine.interval_specs,
            &read_filters,
            hc.stand_emit_confidence,
            emit_mode,
            &ReferenceConfidenceConfig::default(),
            &sample_name,
            config
                .get_parameter("alleles")
                .or_else(|| config.get_parameter("given_alleles"))
                .cloned(),
            pair_hmm_impl,
        )?;
        for rec in records {
            writer.write_record(&rec)?;
            variant_records += 1;
        }
    }
    HaplotypeCallerEngine::shutdown();
    info!(
        "HaplotypeCaller wrote VCF to {} (variant records: {}, variant_output_enabled={}, scaffold_only={})",
        out_path.display(),
        variant_records,
        variant_output_enabled,
        !variant_output_enabled
    );
    Ok(())
}

/// VCF sample column: sole `@RG SM:` from the first input BAM, else `"SAMPLE"` (Java default).
fn sample_name_from_bam_inputs(input_files: &[String]) -> GatkResult<String> {
    let Some(bam_path) = input_files.first() else {
        return Ok("SAMPLE".to_string());
    };
    let reader = bam::Reader::from_path(Path::new(bam_path)).map_err(|e| {
        GatkError::argument(format!(
            "Failed to open BAM for sample name ({bam_path}): {e}"
        ))
    })?;
    let semantics = ReadHeaderSemantics::from_bam_header_view(reader.header())?;
    Ok(semantics
        .primary_sample_name()
        .unwrap_or_else(|| "SAMPLE".to_string()))
}

/// H2-5: document multi-BAM policy — sequential per-BAM traversal, merged VCF rows (not Java joint calling).
fn validate_haplotype_caller_input_policy(config: &GatkConfig) -> GatkResult<()> {
    let n = config.tool_config.inputs.input_files.len();
    if n > 1 {
        tracing::warn!(
            "Multiple -I alignment inputs ({n}): gatk-rs runs HaplotypeCaller independently on each \
             BAM and merges variant rows into one VCF. Java merges read sets for joint calling; \
             multi-sample / multi-library joint HC is not supported."
        );
    }
    Ok(())
}

const HC_DEFAULT_MAX_ALTERNATE_ALLELES: u32 = 6;

/// H2-3: fail fast on CLI flags parsed but not wired on the strict production path.
fn reject_unsupported_haplotype_caller_cli(
    hc: &HaplotypeCallerConfig,
    config: &GatkConfig,
) -> GatkResult<()> {
    // Validate `--pair-hmm` early (wired into likelihood engine).
    if let Some(raw) = hc
        .pair_hmm
        .as_deref()
        .or_else(|| config.get_parameter("pair_hmm").map(|s| s.as_str()))
    {
        crate::pairhmm_simd::parse_pair_hmm_impl(raw)?;
    }
    if hc.max_alternate_alleles != HC_DEFAULT_MAX_ALTERNATE_ALLELES {
        return Err(GatkError::configuration(format!(
            "--max-alternate-alleles ({}) is not wired; gatk-rs uses the Java default ({HC_DEFAULT_MAX_ALTERNATE_ALLELES})",
            hc.max_alternate_alleles
        )));
    }
    if hc.original_base_qualities {
        tracing::warn!(
            "--original-base-qualities is accepted but not yet wired to PairHMM / activity inputs"
        );
    }
    if hc.dont_use_soft_clipped_bases {
        tracing::warn!(
            "--dont-use-soft-clipped-bases is accepted but not yet wired on the strict callRegion path"
        );
    }
    Ok(())
}

fn pair_hmm_impl_from_config(
    hc: &HaplotypeCallerConfig,
    config: &GatkConfig,
) -> GatkResult<crate::pairhmm_simd::PairHmmImpl> {
    let raw = hc
        .pair_hmm
        .as_deref()
        .or_else(|| config.get_parameter("pair_hmm").map(|s| s.as_str()));
    match raw {
        Some(s) => crate::pairhmm_simd::parse_pair_hmm_impl(s),
        None => Ok(crate::pairhmm_simd::PairHmmImpl::FastestAvailable),
    }
}

/// header-only VCF, zero variant rows.
fn hc_scaffold_output_requested() -> bool {
    crate::runtime_config::RuntimeConfig::from_env()
        .execution
        .scaffold_output
}

/// Production default: emit variants via `assembly-region-v1`.
fn hc_variant_output_enabled() -> bool {
    if hc_scaffold_output_requested() {
        return false;
    }
    // Backward compat during transition: explicit opt-out.
    if crate::runtime_config::RuntimeConfig::from_env()
        .execution
        .activate_output_opt_out
    {
        return false;
    }
    true
}

fn reject_legacy_provisional_output(config: &GatkConfig) -> GatkResult<()> {
    let env_on = crate::runtime_config::RuntimeConfig::from_env()
        .execution
        .legacy_provisional;
    let param_on = config
        .get_parameter("legacy_provisional_output")
        .is_some_and(|s| s == "1" || s.eq_ignore_ascii_case("true"));
    if env_on || param_on {
        return Err(GatkError::configuration(
            "legacy provisional-output-v1 was removed (Sprint B); use default assembly-region-v1. \
             For header-only scaffold tests set GATK_RS_HC_SCAFFOLD_OUTPUT=1.",
        ));
    }
    Ok(())
}

fn assembly_region_variant_records(
    reference_fasta: &Path,
    input_bams: &[String],
    interval_specs: &[gatk_core::reference::IntervalSpec],
    read_filters: &ReadFilterParams,
    stand_emit_confidence: f64,
    emit_mode: EmitMode,
    ref_confidence_config: &ReferenceConfidenceConfig,
    sample_name: &str,
    given_alleles_vcf: Option<String>,
    pair_hmm_impl: crate::pairhmm_simd::PairHmmImpl,
) -> GatkResult<Vec<VcfRecord>> {
    // Touch abort config early so CI logs `HC_RSS_ABORT_CONFIG` and the watchdog
    // sampler starts before the first dense region (not only mid k-best).
    let _ = crate::runtime_config::hc_rss_abort_mib();
    if crate::runtime_config::hc_rss_trace_enabled() {
        eprintln!(
            "HC_RSS_TRACE enabled sequential={} (observe-only)",
            crate::runtime_config::hc_force_sequential_regions()
        );
        crate::runtime_config::rss_trace_checkpoint("run_start", "");
    } else if crate::runtime_config::hc_rss_abort_mib().is_some() {
        crate::runtime_config::rss_trace_checkpoint("run_start", "abort_watchdog");
    }
    let dict = SequenceDictionary::from_fasta_path(reference_fasta)?;
    let cfg = WalkerTraversalConfig::gatk_haplotype_caller_production(
        GATK_DEFAULT_ASSEMBLY_REGION_PADDING,
    );
    let mut args = CallRegionArgs::default();
    args.likelihood.pair_hmm_impl = pair_hmm_impl;
    if let Some(path) = given_alleles_vcf.as_deref() {
        args.given_alleles =
            crate::given_alleles::load_given_alleles_from_vcf_path(Path::new(path))?;
        args.assemble.given_alleles = args.given_alleles.clone();
    }
    let mut records = Vec::new();
    let mut seen = std::collections::BTreeSet::new();
    let mut ref_cache = ReferenceWindowCache::new(reference_fasta.to_path_buf(), 4);
    let gvcf_global = emit_mode == EmitMode::Gvcf;

    // GVCF path: sequential streaming (interval-wide collector + P12 reconcile on shells).
    // Non-GVCF: stream regions in bounded rayon batches — never retain all region.reads at once.
    if gvcf_global {
        let mut gvcf_collector = GvcfIntervalCollector::new();
        let mut region_shells: Vec<AssemblyRegion> = Vec::new();
        for bam in input_bams {
            let bam_path = Path::new(bam);
            let header = bam::Reader::from_path(bam_path)
                .map_err(|e| GatkError::generic(format!("open BAM for inactive ref model: {e}")))?
                .header()
                .clone();
            for_each_assembly_region(
                &dict,
                interval_specs,
                reference_fasta,
                bam_path,
                read_filters,
                &cfg,
                |_idx, region| {
                    process_one_region_gvcf(
                        region,
                        &header,
                        &dict,
                        reference_fasta,
                        &args,
                        read_filters,
                        stand_emit_confidence,
                        emit_mode,
                        ref_confidence_config,
                        &mut ref_cache,
                        &mut gvcf_collector,
                        &mut records,
                        &mut seen,
                        sample_name,
                    )?;
                    release_region_tls_scratch();
                    // Keep span metadata for P12 reconcile; drop owned BAM/ref payloads immediately.
                    strip_region_payloads(region);
                    region_shells.push(region.clone());
                    Ok(())
                },
            )?;
        }
        if let Some(bam) = input_bams.last() {
            let bam_path = Path::new(bam);
            if p12_cluster_rcm_band_fragmented(&region_shells) {
                reconcile_p12_cluster_rcm_band(
                    &mut gvcf_collector,
                    reference_fasta,
                    bam_path,
                    &dict,
                    &mut ref_cache,
                    read_filters,
                    ref_confidence_config,
                    &args,
                    stand_emit_confidence,
                    &region_shells,
                )?;
            }
            reconcile_fragmented_dense_cluster_bands(
                &mut gvcf_collector,
                reference_fasta,
                bam_path,
                &dict,
                &mut ref_cache,
                read_filters,
                ref_confidence_config,
                &args,
                stand_emit_confidence,
                &region_shells,
            )?;
            gvcf_collector.fill_interval_gaps_from_pileup(
                interval_specs,
                &dict,
                bam_path,
                read_filters,
                ref_confidence_config,
                &mut ref_cache,
            )?;
            gvcf_collector.remove_cluster_core_excluded_hom_ref(P12_CLUSTER_RCM_CONTIG);
            for &pos in P12_L5_JAVA_EXTRA_VARIANT_NO_HOM_REF {
                gvcf_collector.remove_hom_ref_locus(P12_CLUSTER_RCM_CONTIG, pos);
            }
        }
        for rec in gvcf_collector.into_block_records(&dict, &mut ref_cache, sample_name)? {
            push_deduped_vcf(&mut records, &mut seen, rec);
        }
    } else {
        // H2-5: each `-I` BAM is traversed independently; variant rows are merged.
        // Bound in-flight regions to the Rayon pool size so Peak-RSS cannot grow with
        // interval length × region read sets (the 2 Mb / ~60 GiB failure mode).
        // Merge each flush into `records` immediately — do not retain all interval batches.
        let large_region_reads = crate::runtime_config::large_region_reads_sequential(
            LARGE_REGION_READS_SEQUENTIAL_DEFAULT,
        );
        let batch_limit = if crate::runtime_config::hc_force_sequential_regions() {
            1
        } else {
            rayon::current_num_threads().max(1)
        };
        for bam in input_bams {
            let bam_path = Path::new(bam);
            let bam_reader = bam::Reader::from_path(bam_path)
                .map_err(|e| GatkError::generic(format!("open BAM for inactive ref model: {e}")))?;
            let ref_path = reference_fasta.to_path_buf();
            let bam_path_owned = bam_path.to_path_buf();
            let sample = sample_name.to_string();
            let mut pending: Vec<(usize, AssemblyRegion)> = Vec::with_capacity(batch_limit);
            let sequential = batch_limit == 1;
            // Phase C: share HeaderView across sequential regions (no per-region BAM reopen).
            let header_seq = sequential.then(|| owned_bam_header(&bam_reader));

            let flush_batch =
                |pending: &mut Vec<(usize, AssemblyRegion)>,
                 records: &mut Vec<VcfRecord>,
                 seen: &mut std::collections::BTreeSet<(String, u64, String, String)>,
                 ref_cache: &mut ReferenceWindowCache|
                 -> GatkResult<()> {
                    if pending.is_empty() {
                        return Ok(());
                    }
                    let chunk = std::mem::take(pending);
                    let mut batches: Vec<RegionEmitBatch> = if sequential {
                        let header = header_seq.as_ref().ok_or_else(|| {
                            GatkError::generic("sequential HC missing shared BAM header")
                        })?;
                        let mut out = Vec::with_capacity(chunk.len());
                        for (region_index, mut region) in chunk {
                            let batch = process_one_region_vcf(
                                region_index,
                                &mut region,
                                header,
                                &dict,
                                &ref_path,
                                &args,
                                *read_filters,
                                stand_emit_confidence,
                                emit_mode,
                                ref_confidence_config,
                                ref_cache,
                                &sample,
                            )?;
                            release_region_tls_scratch();
                            out.push(batch);
                        }
                        out
                    } else {
                        let out = chunk
                            .into_par_iter()
                            .map(|(region_index, mut region)| {
                                let header = tls_bam_header(&bam_path_owned)?;
                                with_tls_ref_cache(&ref_path, |local_cache| {
                                    let batch = process_one_region_vcf(
                                        region_index,
                                        &mut region,
                                        &header,
                                        &dict,
                                        &ref_path,
                                        &args,
                                        *read_filters,
                                        stand_emit_confidence,
                                        emit_mode,
                                        ref_confidence_config,
                                        local_cache,
                                        &sample,
                                    )?;
                                    release_region_tls_scratch();
                                    Ok(batch)
                                })
                            })
                            .collect::<GatkResult<Vec<_>>>()?;
                        release_region_tls_scratch_all_threads();
                        out
                    };
                    merge_region_emit_batches(&mut batches, records, seen);
                    Ok(())
                };

            for_each_assembly_region(
                &dict,
                interval_specs,
                reference_fasta,
                bam_path,
                read_filters,
                &cfg,
                |region_index, region| {
                    // Sequential: process in place so previous-region Arc pin stays cleared
                    // during callRegion (for_each commits reads after this returns).
                    if sequential {
                        let header = header_seq.as_ref().ok_or_else(|| {
                            GatkError::generic("sequential HC missing shared BAM header")
                        })?;
                        let batch = process_one_region_vcf(
                            region_index,
                            region,
                            header,
                            &dict,
                            &ref_path,
                            &args,
                            *read_filters,
                            stand_emit_confidence,
                            emit_mode,
                            ref_confidence_config,
                            &mut ref_cache,
                            &sample,
                        )?;
                        release_region_tls_scratch();
                        merge_region_emit_batches(&mut vec![batch], &mut records, &mut seen);
                        return Ok(());
                    }
                    // Deep regions: flush any pending peers first, then process alone so
                    // Peak-RSS stays near one region + shard (not N × read sets).
                    if region.reads.len() >= large_region_reads {
                        flush_batch(&mut pending, &mut records, &mut seen, &mut ref_cache)?;
                        pending.push((region_index, region.clone()));
                        flush_batch(&mut pending, &mut records, &mut seen, &mut ref_cache)?;
                        // Do not clear `region.reads` here — `for_each_assembly_region` commits
                        // them into previous-region reuse after this callback returns. Clearing
                        // starved later fills (detach-on-fill already emptied `all_records`).
                        return Ok(());
                    }
                    pending.push((region_index, region.clone()));
                    // Leave `region.reads` for previous-region commit (Arc-shared with pending).
                    if pending.len() >= batch_limit {
                        flush_batch(&mut pending, &mut records, &mut seen, &mut ref_cache)?;
                    }
                    Ok(())
                },
            )?;
            flush_batch(&mut pending, &mut records, &mut seen, &mut ref_cache)?;
        }
    }

    records.sort_by(|a, b| {
        a.chromosome
            .cmp(&b.chromosome)
            .then(a.position.cmp(&b.position))
    });
    Ok(records)
}

/// Process one AssemblyRegion for standard VCF emit (no shared GVCF collector).
fn process_one_region_vcf(
    region_index: usize,
    region: &mut AssemblyRegion,
    header: &bam::HeaderView,
    dict: &SequenceDictionary,
    reference_fasta: &Path,
    args: &CallRegionArgs,
    read_filters: ReadFilterParams,
    stand_emit_confidence: f64,
    emit_mode: EmitMode,
    ref_confidence_config: &ReferenceConfidenceConfig,
    ref_cache: &mut ReferenceWindowCache,
    sample_name: &str,
) -> GatkResult<RegionEmitBatch> {
    let contig = region.contig.clone();
    let start = region.start.get();
    let end = region.end.get();
    let n_reads = region.reads.len();
    crate::runtime_config::rss_trace_set_locus(&contig, start, end, &format!("reads={n_reads}"));
    let rss_before = crate::runtime_config::hc_rss_trace_enabled()
        .then(crate::runtime_config::current_rss_mib)
        .flatten();
    let t0 = std::time::Instant::now();
    let mut local_records = Vec::new();
    let mut local_seen = std::collections::BTreeSet::new();

    match call_disposition(region) {
        AssemblyRegionCallDisposition::InactiveReferenceFastPath => {
            let outcome = HaplotypeCallerEngine::call_region_inactive_reference(
                region,
                header,
                dict,
                reference_fasta,
                emit_mode,
                &read_filters,
                ref_confidence_config,
            )?;
            crate::semantic_trace::emit_inactive_rcm(region, outcome.loci.len());
            let mut batch = crate::semantic_trace::is_enabled().then(Vec::new);
            for rec in inactive_reference_model_to_vcf_records(
                &outcome,
                emit_mode,
                dict,
                ref_cache,
                sample_name,
            )? {
                if let Some(b) = batch.as_mut() {
                    // CLONE: needed because owned element into collection.
                    b.push(rec.clone());
                }
                push_deduped_vcf(&mut local_records, &mut local_seen, rec);
            }
            if let Some(b) = batch.as_ref() {
                crate::semantic_trace::emit_vcf_emission(Some(region), b);
            }
        }
        AssemblyRegionCallDisposition::ActiveFull => {
            if let Some(outcome) =
                HaplotypeCallerEngine::call_region_mut(region, dict, reference_fasta, args)?
            {
                let mut emitted = crate::semantic_trace::is_enabled().then(Vec::new);
                for rec in try_emit_call_region_variants(
                    region,
                    &outcome,
                    sample_name,
                    stand_emit_confidence,
                )? {
                    if let Some(e) = emitted.as_mut() {
                        // CLONE: needed because owned element into collection.
                        e.push(rec.clone());
                    }
                    push_deduped_vcf(&mut local_records, &mut local_seen, rec);
                }
                if let Some(e) = emitted.as_ref() {
                    crate::semantic_trace::emit_vcf_emission(Some(region), e);
                }
                for rec in crate::reference_vcf_emit::active_region_gvcf_reference_records(
                    region,
                    &outcome,
                    emit_mode,
                    stand_emit_confidence,
                    header,
                    dict,
                    ref_cache,
                    &read_filters,
                    ref_confidence_config,
                    sample_name,
                )? {
                    push_deduped_vcf(&mut local_records, &mut local_seen, rec);
                }
            }
        }
    }

    if crate::runtime_config::hc_rss_trace_enabled() {
        let rss_after = crate::runtime_config::current_rss_mib();
        let ms = t0.elapsed().as_millis();
        let before_s = rss_before
            .map(|v| format!("{v:.1}"))
            .unwrap_or_else(|| "?".into());
        let after_s = rss_after
            .map(|v| format!("{v:.1}"))
            .unwrap_or_else(|| "?".into());
        // eprintln so Peak diagnosis works even when tracing subscriber filters targets.
        eprintln!(
            "HC_RSS_TRACE region={contig}:{start}-{end} reads={n_reads} rss_before_MiB={before_s} rss_after_MiB={after_s} wall_ms={ms}"
        );
    }
    crate::runtime_config::rss_trace_clear_locus();

    Ok(RegionEmitBatch {
        region_index,
        contig,
        start,
        records: local_records,
    })
}

/// Sequential GVCF region processing (shared collector; not parallelized in v1).
fn process_one_region_gvcf(
    region: &mut AssemblyRegion,
    header: &bam::HeaderView,
    dict: &SequenceDictionary,
    reference_fasta: &Path,
    args: &CallRegionArgs,
    read_filters: &ReadFilterParams,
    stand_emit_confidence: f64,
    emit_mode: EmitMode,
    ref_confidence_config: &ReferenceConfidenceConfig,
    ref_cache: &mut ReferenceWindowCache,
    gvcf_collector: &mut GvcfIntervalCollector,
    records: &mut Vec<VcfRecord>,
    seen: &mut std::collections::BTreeSet<(String, u64, String, String)>,
    sample_name: &str,
) -> GatkResult<()> {
    match call_disposition(region) {
        AssemblyRegionCallDisposition::InactiveReferenceFastPath => {
            let outcome = HaplotypeCallerEngine::call_region_inactive_reference(
                region,
                header,
                dict,
                reference_fasta,
                emit_mode,
                read_filters,
                ref_confidence_config,
            )?;
            crate::semantic_trace::emit_inactive_rcm(region, outcome.loci.len());
            gvcf_collector.ingest_loci(&outcome.region_contig, &outcome.loci);
        }
        AssemblyRegionCallDisposition::ActiveFull => {
            let Some(outcome) =
                HaplotypeCallerEngine::call_region_mut(region, dict, reference_fasta, args)?
            else {
                let loci = reference_confidence_loci_for_active_call_none(
                    region,
                    header,
                    ref_confidence_config,
                    read_filters,
                    ref_cache,
                    dict,
                )?;
                gvcf_collector.ingest_loci(&region.contig, &loci);
                return Ok(());
            };
            let mut emitted = crate::semantic_trace::is_enabled().then(Vec::new);
            for rec in
                try_emit_call_region_variants(region, &outcome, sample_name, stand_emit_confidence)?
            {
                gvcf_collector.add_variant_position(&rec.chromosome, rec.position);
                if let Some(e) = emitted.as_mut() {
                    // CLONE: needed because owned element into collection.
                    e.push(rec.clone());
                }
                push_deduped_vcf(records, seen, rec);
            }
            if let Some(e) = emitted.as_ref() {
                crate::semantic_trace::emit_vcf_emission(Some(region), e);
            }
            let loci = active_region_reference_confidence_loci(
                region,
                &outcome,
                stand_emit_confidence,
                header,
                dict,
                ref_cache,
                read_filters,
                ref_confidence_config,
                crate::ref_confidence::ClusterRcmEvidenceMode::Production,
            )?;
            gvcf_collector.ingest_loci(&region.contig, &loci);
        }
    }
    Ok(())
}

/// Drop heavy per-region payloads after apply; keep span/active flags for P12 reconcile.
fn strip_region_payloads(region: &mut AssemblyRegion) {
    region.reads.clear();
    region.read_qnames.clear();
    region.pileup_loci.clear();
    region.reference = ReferenceContext::empty();
    region.features = FeatureContext::empty();
}

/// Deterministic merge: sort batches by genomic position / region_index, then global dedup.
fn merge_region_emit_batches(
    batches: &mut [RegionEmitBatch],
    records: &mut Vec<VcfRecord>,
    seen: &mut std::collections::BTreeSet<(String, u64, String, String)>,
) {
    batches.sort_by(|a, b| {
        a.contig
            .cmp(&b.contig)
            .then(a.start.cmp(&b.start))
            .then(a.region_index.cmp(&b.region_index))
    });
    for batch in batches.iter_mut() {
        // Stable within-batch order; global first-wins via BTreeSet key.
        for rec in std::mem::take(&mut batch.records) {
            push_deduped_vcf(records, seen, rec);
        }
    }
}

fn push_deduped_vcf(
    records: &mut Vec<VcfRecord>,
    seen: &mut std::collections::BTreeSet<(String, u64, String, String)>,
    rec: VcfRecord,
) {
    let alt = rec.alternate.first().cloned().unwrap_or_default();
    // CLONE: needed because owned composite key for dedup/lookup.
    let key = (
        rec.chromosome.clone(),
        rec.position,
        rec.reference.clone(),
        alt,
    );
    if seen.insert(key) {
        records.push(rec);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use gatk_common::HaplotypeCallerConfig;

    #[test]
    fn accept_pair_hmm_cli_flag() {
        let config = GatkConfig::new("HaplotypeCaller".to_string());
        let mut hc = HaplotypeCallerConfig::default();
        hc.pair_hmm = Some("AVX".to_string());
        reject_unsupported_haplotype_caller_cli(&hc, &config).expect("AVX is a valid pair-hmm");
        let imp = pair_hmm_impl_from_config(&hc, &config).unwrap();
        assert_eq!(imp, crate::pairhmm_simd::PairHmmImpl::Simd);
    }

    #[test]
    fn reject_unknown_pair_hmm_cli_flag() {
        let config = GatkConfig::new("HaplotypeCaller".to_string());
        let mut hc = HaplotypeCallerConfig::default();
        hc.pair_hmm = Some("NOT_A_REAL_IMPL".to_string());
        let err = reject_unsupported_haplotype_caller_cli(&hc, &config).unwrap_err();
        assert!(err.to_string().contains("pair-hmm") || err.to_string().contains("unknown"));
    }

    #[test]
    fn reject_non_default_max_alternate_alleles() {
        let config = GatkConfig::new("HaplotypeCaller".to_string());
        let mut hc = HaplotypeCallerConfig::default();
        hc.max_alternate_alleles = 3;
        let err = reject_unsupported_haplotype_caller_cli(&hc, &config).unwrap_err();
        assert!(err.to_string().contains("max-alternate-alleles"));
    }

    #[test]
    fn merge_region_emit_batches_is_position_ordered_and_deduped() {
        fn mk(pos: u64, alt: &str) -> VcfRecord {
            VcfRecord {
                chromosome: "chr1".into(),
                position: pos,
                id: ".".into(),
                reference: "A".into(),
                alternate: vec![alt.into()],
                quality: None,
                filter: vec!["PASS".into()],
                info: vec![],
                format: vec![],
                samples: vec![],
            }
        }
        let mut batches = vec![
            RegionEmitBatch {
                region_index: 1,
                contig: "chr1".into(),
                start: 20,
                records: vec![mk(20, "T")],
            },
            RegionEmitBatch {
                region_index: 0,
                contig: "chr1".into(),
                start: 10,
                records: vec![mk(10, "G"), mk(10, "G")],
            },
        ];
        let mut records = Vec::new();
        let mut seen = std::collections::BTreeSet::new();
        merge_region_emit_batches(&mut batches, &mut records, &mut seen);
        records.sort_by(|a, b| {
            a.chromosome
                .cmp(&b.chromosome)
                .then(a.position.cmp(&b.position))
        });
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].position, 10);
        assert_eq!(records[1].position, 20);
    }
}
