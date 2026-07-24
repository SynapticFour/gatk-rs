//! Small dumps for `scripts/parity` gates (read shards, assembly regions, apply summary).
//! Usage:
//! `read-shards <ref.fa> <interval_cli> [padding]` — prints `contig\\tspan_start\\tspan_end` per padded span.
//! `assembly-regions <ref.fa> <bam> <interval_cli> [padding]` — region rows via full walker traversal (B.2).
//! `apply-summary` / `walker-traversal-summary` — [`WalkerApplyStats`] with HC production shard pipeline (B.3 / B.4).
//! `raw-activity` / `smoothed-activity` / `active-locus` — activity dumps.

use gatk_common::try_init_from_gatk_rs_hc_trace;
use gatk_core::reference::{parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary};
use gatk_haplotypecaller::genome_loc::GenomePosition;
use gatk_haplotypecaller::{
    dump_active_locus_tsv, dump_af_em_tsv, dump_allele_biased_evidence_locus_tsv,
    dump_allele_subsetting_tsv, dump_annotate_core_tsv, dump_annotation_manifest_tsv,
    dump_annotation_plugin_tsv, dump_as_annotations_tsv, dump_assembler_args_tsv,
    dump_assembly_assemble_tsv, dump_assembly_debug_stub_tsv,
    dump_assembly_graph_dangling_summary_tsv, dump_assembly_graph_edges_tsv,
    dump_assembly_graph_low_quality_tsv, dump_assembly_graph_multi_kmer_edges_tsv,
    dump_assembly_graph_non_unique_summary_tsv, dump_assembly_graph_pruned_summary_tsv,
    dump_assembly_haplotype_cigars_tsv, dump_assembly_haplotypes_cap_tsv,
    dump_assembly_haplotypes_production_tsv, dump_assembly_haplotypes_tsv,
    dump_assembly_junction_haplotypes_tsv, dump_assembly_kbest_paths_tsv,
    dump_assembly_region_assembly_stages_finalize_tsv, dump_assembly_region_assembly_stages_tsv,
    dump_assembly_region_features_tsv, dump_assembly_region_finalize_reads_tsv,
    dump_assembly_region_genotype_subset_tsv, dump_assembly_region_genotype_tsv,
    dump_assembly_region_haplotypes_tsv, dump_assembly_region_kbest_paths_tsv,
    dump_assembly_region_kmer_probe_tsv, dump_assembly_region_pairhmm_likelihoods_tsv,
    dump_assembly_region_pileup_track_tsv, dump_assembly_region_reads_tsv,
    dump_assembly_region_reference_tsv, dump_assembly_region_trim_tsv,
    dump_assembly_seqgraph_summary_tsv, dump_bamout_stub_tsv, dump_call_region_active_rcm_loci_tsv,
    dump_call_region_format_tsv, dump_call_region_vcf_tsv, dump_depth_per_sample_hc_tsv,
    dump_dragen_mode_branch_tsv, dump_dragstr_calibration_tsv, dump_emit_mode_decision_tsv,
    dump_excess_het_tsv, dump_force_calling_genotype_tsv, dump_genotype_format_tsv,
    dump_genotype_likelihood_activity_tsv, dump_genotype_limits_tsv, dump_genotype_phasing_tsv,
    dump_genotyping_aggregate_tsv, dump_gvcf_header_tsv, dump_gvcf_l5_merged_tsv,
    dump_gvcf_writer_from_loci_fixture_tsv, dump_hc_read_filter_tsv, dump_hq_soft_clip_mean_tsv,
    dump_inactive_reference_model_tsv, dump_likelihood_engine_config_tsv,
    dump_likelihood_pcr_read_tsv, dump_locus_pileup_detail_tsv, dump_locus_pileup_tsv,
    dump_pairhmm_bq_cap_tsv, dump_pairhmm_haplotype_filter_tsv, dump_pairhmm_likelihoods_tsv,
    dump_pairhmm_native_likelihoods_tsv, dump_pcr_error_model_tsv, dump_ploidy_resolution_tsv,
    dump_positional_downsample_summary_tsv, dump_raw_activity_profile_tsv_with_contamination,
    dump_raw_activity_profile_tsv_with_force_calling, dump_read_error_correction_tsv,
    dump_read_pre_len_tsv, dump_read_pre_mq_tsv, dump_read_pre_overlap_tsv,
    dump_read_pre_softclip_tsv, dump_read_shard_pipeline_tsv, dump_ref_confidence_merge_case,
    dump_reference_confidence_locus_tsv, dump_smoothed_activity_profile_tsv,
    dump_standard_annotations_tsv, dump_subset_alleles_integration_tsv, dump_subset_alleles_pl_tsv,
    dump_subset_alleles_vc_tsv, dump_target_allele_counts_tsv, dump_variant_format_from_gl_ad_tsv,
    dump_variant_vcf_from_gl_ad_tsv, flatten_assembly_regions, load_trim_variants_tsv,
    traverse_assembly_region_walker, AssemblyGraphPruningParams, AssemblyRegion,
    AssemblyRegionHaplotypeTarget, AssemblyRegionTrimmer, AssemblyRegionTrimmerConfig,
    FeatureContext, FeatureDataSources, ForceCallingAllelesDump, ReadFilterParams,
    ReferenceContext, WalkerApplyStats, WalkerTraversalConfig,
    GATK_DEFAULT_ASSEMBLY_REGION_PADDING, GATK_HC_ALLELES_FEATURE_SOURCE,
};
use rust_htslib::bam;
use rust_htslib::bam::Read as _;
use std::env;
use std::io;
use std::path::{Path, PathBuf};
use std::process::{Command, ExitCode};

fn usage_and_exit() -> ! {
    eprintln!(
        "Usage:\n  … Phase B–D: (see prior gates)\n  … Phase E: assembly-graph … | assembly-haplotype-cigars …\n  … Phase F: pairhmm-likelihoods <cases.tsv>"
    );
    std::process::exit(2);
}

fn main() -> ExitCode {
    let _ = try_init_from_gatk_rs_hc_trace();
    gatk_haplotypecaller::semantic_trace::try_init_from_runtime(
        &gatk_haplotypecaller::runtime_config::RuntimeConfig::from_env(),
    );
    if let Err(e) = run() {
        eprintln!("{e}");
        return ExitCode::from(1);
    }
    ExitCode::SUCCESS
}

fn next_arg(it: &mut impl Iterator<Item = String>, what: &str) -> String {
    it.next().unwrap_or_else(|| {
        eprintln!("missing {what}");
        usage_and_exit();
    })
}

fn parse_padding_token(token: Option<String>) -> Result<u64, String> {
    match token {
        None => Ok(GATK_DEFAULT_ASSEMBLY_REGION_PADDING),
        Some(s) if s == "-" => Ok(GATK_DEFAULT_ASSEMBLY_REGION_PADDING),
        Some(s) => s.parse().map_err(|_| "padding: invalid u64".to_string()),
    }
}

fn parse_padding(it: &mut impl Iterator<Item = String>) -> Result<u64, String> {
    parse_padding_token(it.next())
}

fn parse_padding_and_region_target(
    it: &mut impl Iterator<Item = String>,
) -> Result<(u64, AssemblyRegionHaplotypeTarget), String> {
    let first = it.next();
    match first.as_deref() {
        None => Ok((
            GATK_DEFAULT_ASSEMBLY_REGION_PADDING,
            AssemblyRegionHaplotypeTarget::Active,
        )),
        Some("active") => Ok((
            GATK_DEFAULT_ASSEMBLY_REGION_PADDING,
            AssemblyRegionHaplotypeTarget::Active,
        )),
        Some("inactive") => Ok((
            GATK_DEFAULT_ASSEMBLY_REGION_PADDING,
            AssemblyRegionHaplotypeTarget::Inactive,
        )),
        Some("-") => {
            let target = it
                .next()
                .and_then(|s| AssemblyRegionHaplotypeTarget::parse_cli(&s))
                .unwrap_or(AssemblyRegionHaplotypeTarget::Active);
            Ok((GATK_DEFAULT_ASSEMBLY_REGION_PADDING, target))
        }
        Some(s) => {
            let padding: u64 = s.parse().map_err(|_| "padding: invalid u64".to_string())?;
            let target = it
                .next()
                .and_then(|t| AssemblyRegionHaplotypeTarget::parse_cli(&t))
                .unwrap_or(AssemblyRegionHaplotypeTarget::Active);
            Ok((padding, target))
        }
    }
}

/// Match `run_hc_full_parity_java_refresh.sh`: SAM fixtures → coordinate-sorted indexed BAM.
/// Rebuilds cached BAM when the source SAM is newer (avoids stale gates after editing fixtures).
fn resolve_alignment_for_parity(bam: &Path) -> Result<PathBuf, String> {
    if bam.extension().and_then(|s| s.to_str()) != Some("sam") {
        return Ok(bam.to_path_buf());
    }
    let cache = Path::new(env!("CARGO_MANIFEST_DIR")).join("../parity/build/sam-indexed-bam");
    std::fs::create_dir_all(&cache).map_err(|e| format!("mkdir sam cache: {e}"))?;
    let stem = bam.file_stem().and_then(|s| s.to_str()).unwrap_or("input");
    let out = cache.join(format!("{stem}.bam"));
    let sam_meta = std::fs::metadata(bam).map_err(|e| format!("stat sam: {e}"))?;
    let need_rebuild = match std::fs::metadata(&out) {
        Ok(bam_meta) => {
            let sam_t = sam_meta.modified().map_err(|e| format!("sam mtime: {e}"))?;
            let bam_t = bam_meta.modified().map_err(|e| format!("bam mtime: {e}"))?;
            sam_t > bam_t
        }
        Err(_) => true,
    };
    if !need_rebuild && out.is_file() {
        return Ok(out);
    }
    if out.is_file() {
        let _ = std::fs::remove_file(&out);
        let _ = std::fs::remove_file(out.with_extension("bam.bai"));
    }
    eprintln!("[hc_full_parity_gate_dump] sam->bam {}", bam.display());
    let sort_status = Command::new("sh")
        .arg("-c")
        .arg(format!(
            "samtools view -bS {} | samtools sort -o {}",
            bam.display(),
            out.display()
        ))
        .status()
        .map_err(|e| format!("samtools: {e}"))?;
    if !sort_status.success() {
        return Err(format!("samtools sort failed for {}", bam.display()));
    }
    let idx_status = Command::new("samtools")
        .args(["index", out.to_str().unwrap()])
        .status()
        .map_err(|e| format!("samtools index: {e}"))?;
    if !idx_status.success() {
        return Err(format!("samtools index failed for {}", out.display()));
    }
    Ok(out)
}

fn load_context(
    ref_fa: &Path,
    interval_cli: &str,
) -> Result<(SequenceDictionary, Vec<gatk_core::reference::IntervalSpec>), String> {
    let dict = SequenceDictionary::from_fasta_path(ref_fa).map_err(|e| format!("dict: {e}"))?;
    let specs =
        parse_intervals_cli_string(&dict, interval_cli).map_err(|e| format!("intervals: {e}"))?;
    Ok((dict, specs))
}

fn run_walker(
    ref_fa: &Path,
    bam: &Path,
    interval_cli: &str,
    padding: u64,
    force_active: bool,
    features_vcf: Option<&Path>,
    track_pileups: bool,
) -> Result<(Vec<AssemblyRegion>, WalkerApplyStats), String> {
    let (dict, specs) = load_context(ref_fa, interval_cli)?;
    let alignment = resolve_alignment_for_parity(bam)?;
    let filters = ReadFilterParams::gatk_standard_hc();
    let mut cfg = WalkerTraversalConfig::gatk_haplotype_caller_production(padding);
    cfg.force_active = force_active;
    cfg.track_pileups = track_pileups;
    if let Some(vcf) = features_vcf {
        let mut sources = FeatureDataSources::default();
        sources
            .load_vcf_source("alleles", vcf)
            .map_err(|e| format!("feature VCF: {e}"))?;
        cfg.feature_sources = Some(sources);
    }
    let walk = traverse_assembly_region_walker(&dict, &specs, ref_fa, &alignment, &filters, &cfg)
        .map_err(|e| format!("traverse: {e}"))?;
    let regions = flatten_assembly_regions(&walk);
    Ok((regions, walk.apply_stats))
}

fn print_regions(regions: &[AssemblyRegion]) {
    println!("contig\tstart\tend\tis_active\textended_start\textended_end\textension");
    for r in regions {
        println!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.contig,
            r.start.get(),
            r.end.get(),
            r.is_active,
            r.extended_start.get(),
            r.extended_end.get(),
            r.extension
        );
    }
}

fn print_apply_summary(st: &WalkerApplyStats) {
    println!("total_apply\tinactive_fast_path\tactive_full");
    println!(
        "{}\t{}\t{}",
        st.total_apply, st.inactive_fast_path, st.active_full
    );
}

fn run() -> Result<(), String> {
    let mut it = env::args().skip(1);
    let cmd = it.next().unwrap_or_else(|| usage_and_exit());
    match cmd.as_str() {
        "read-shards" => {
            let ref_path = next_arg(&mut it, "ref.fa");
            let ref_fa = Path::new(&ref_path);
            let interval_cli = next_arg(&mut it, "interval_cli");
            let padding = parse_padding(&mut it)?;
            let (dict, specs) = load_context(ref_fa, interval_cli.as_str())?;
            let shards = gatk_haplotypecaller::make_read_shards(&dict, &specs, padding)
                .map_err(|e| format!("shards: {e}"))?;
            for s in &shards {
                for &(lo, hi) in &s.padded_spans {
                    println!("{}\t{}\t{}", s.contig, lo, hi);
                }
            }
        }
        "assembly-regions" | "assembly-regions-force-active" => {
            let ref_path = next_arg(&mut it, "ref.fa");
            let ref_fa = Path::new(&ref_path);
            let bam_path = next_arg(&mut it, "bam");
            let bam = resolve_alignment_for_parity(Path::new(&bam_path))?;
            let interval_cli = next_arg(&mut it, "interval_cli");
            let padding = parse_padding(&mut it)?;
            let force_active = cmd == "assembly-regions-force-active";
            let (regions, _) = run_walker(
                ref_fa,
                &bam,
                interval_cli.as_str(),
                padding,
                force_active,
                None,
                false,
            )?;
            print_regions(&regions);
        }
        "assembly-region-reads" => {
            let ref_path = next_arg(&mut it, "ref.fa");
            let ref_fa = Path::new(&ref_path);
            let bam_path = next_arg(&mut it, "bam");
            let bam = resolve_alignment_for_parity(Path::new(&bam_path))?;
            let interval_cli = next_arg(&mut it, "interval_cli");
            let padding = parse_padding(&mut it)?;
            let (regions, _) = run_walker(
                ref_fa,
                &bam,
                interval_cli.as_str(),
                padding,
                false,
                None,
                false,
            )?;
            let reader = bam::Reader::from_path(&bam).map_err(|e| format!("open bam: {e}"))?;
            let hdr = reader.header().clone();
            let filters = ReadFilterParams::gatk_standard_hc();
            let mut stdout = io::stdout().lock();
            dump_assembly_region_reads_tsv(&regions, &mut stdout, &hdr, &filters)
                .map_err(|e| format!("assembly-region-reads: {e}"))?;
        }
        "assembly-region-reference" => {
            let ref_path = next_arg(&mut it, "ref.fa");
            let ref_fa = Path::new(&ref_path);
            let bam_path = next_arg(&mut it, "bam");
            let bam = resolve_alignment_for_parity(Path::new(&bam_path))?;
            let interval_cli = next_arg(&mut it, "interval_cli");
            let padding = parse_padding(&mut it)?;
            let (regions, _) = run_walker(
                ref_fa,
                &bam,
                interval_cli.as_str(),
                padding,
                false,
                None,
                false,
            )?;
            let mut stdout = io::stdout().lock();
            dump_assembly_region_reference_tsv(&regions, &mut stdout)
                .map_err(|e| format!("assembly-region-reference: {e}"))?;
        }
        "assembly-region-features" => {
            let ref_path = next_arg(&mut it, "ref.fa");
            let ref_fa = Path::new(&ref_path);
            let bam_path = next_arg(&mut it, "bam");
            let bam = resolve_alignment_for_parity(Path::new(&bam_path))?;
            let interval_cli = next_arg(&mut it, "interval_cli");
            let padding = parse_padding(&mut it)?;
            let features_vcf_path = it.next().filter(|s| s != "-");
            let features_vcf = features_vcf_path.as_deref().map(Path::new);
            let (regions, _) = run_walker(
                ref_fa,
                &bam,
                interval_cli.as_str(),
                padding,
                false,
                features_vcf,
                false,
            )?;
            let mut stdout = io::stdout().lock();
            dump_assembly_region_features_tsv(&regions, &mut stdout)
                .map_err(|e| format!("assembly-region-features: {e}"))?;
        }
        "assembly-region-pileup-track" => {
            let ref_path = next_arg(&mut it, "ref.fa");
            let ref_fa = Path::new(&ref_path);
            let bam_path = next_arg(&mut it, "bam");
            let bam = resolve_alignment_for_parity(Path::new(&bam_path))?;
            let interval_cli = next_arg(&mut it, "interval_cli");
            let padding = parse_padding(&mut it)?;
            let track = matches!(
                it.next().as_deref(),
                Some("1") | Some("true") | Some("track")
            );
            let (regions, _) = run_walker(
                ref_fa,
                &bam,
                interval_cli.as_str(),
                padding,
                false,
                None,
                track,
            )?;
            let mut stdout = io::stdout().lock();
            dump_assembly_region_pileup_track_tsv(&regions, track, &mut stdout)
                .map_err(|e| format!("assembly-region-pileup-track: {e}"))?;
        }
        "assembly-region-trim" => {
            let ref_path = next_arg(&mut it, "ref.fa");
            let ref_fa = Path::new(&ref_path);
            let contig = next_arg(&mut it, "contig");
            let start = GenomePosition::new_1based(
                next_arg(&mut it, "start")
                    .parse()
                    .map_err(|_| "start: invalid u64".to_string())?,
            );
            let end = GenomePosition::new_1based(
                next_arg(&mut it, "end")
                    .parse()
                    .map_err(|_| "end: invalid u64".to_string())?,
            );
            let ext_start: u64 = next_arg(&mut it, "ext_start")
                .parse()
                .map_err(|_| "ext_start: invalid u64".to_string())?;
            let ext_end: u64 = next_arg(&mut it, "ext_end")
                .parse()
                .map_err(|_| "ext_end: invalid u64".to_string())?;
            let variants_arg = it.next().unwrap_or_else(|| "-".to_string());
            let legacy = matches!(
                it.next().as_deref(),
                Some("1") | Some("true") | Some("legacy")
            );
            let dict =
                SequenceDictionary::from_fasta_path(ref_fa).map_err(|e| format!("dict: {e}"))?;
            let mut ref_cache = ReferenceWindowCache::new(ref_fa.to_path_buf(), 4);
            let mut region = AssemblyRegion {
                contig: contig.clone(),
                start,
                end,
                is_active: true,
                extended_start: GenomePosition::new_1based(ext_start),
                extended_end: GenomePosition::new_1based(ext_end),
                extension: 100,
                reads: Vec::new(),
                read_qnames: Vec::new(),
                reference: ReferenceContext::empty(),
                features: FeatureContext::empty(),
                pileup_loci: Vec::new(),
            };
            region.reference =
                ReferenceContext::from_interval(&dict, &mut ref_cache, &contig, ext_start, ext_end)
                    .map_err(|e| format!("reference: {e}"))?;
            let variants = if variants_arg == "-" {
                Vec::new()
            } else {
                load_trim_variants_tsv(Path::new(&variants_arg))?
            };
            let mut trim_cfg = AssemblyRegionTrimmerConfig::gatk_defaults();
            trim_cfg.enable_legacy_assembly_region_trimming = legacy;
            let trimmer = AssemblyRegionTrimmer::new(trim_cfg, &dict, &contig);
            let trim_result = trimmer.trim(&region, &variants, Some(&region.reference));
            let trimmed = AssemblyRegionTrimmer::apply_trim(&region, &trim_result);
            let mut stdout = io::stdout().lock();
            dump_assembly_region_trim_tsv(&region, &trim_result, &trimmed, &mut stdout)?;
        }
        "apply-summary" => {
            let ref_path = next_arg(&mut it, "ref.fa");
            let ref_fa = Path::new(&ref_path);
            let bam_path = next_arg(&mut it, "bam");
            let bam = resolve_alignment_for_parity(Path::new(&bam_path))?;
            let interval_cli = next_arg(&mut it, "interval_cli");
            let padding = parse_padding(&mut it)?;
            let (_, st) = run_walker(
                ref_fa,
                &bam,
                interval_cli.as_str(),
                padding,
                false,
                None,
                false,
            )?;
            print_apply_summary(&st);
        }
        "walker-traversal-summary" => {
            let ref_path = next_arg(&mut it, "ref.fa");
            let ref_fa = Path::new(&ref_path);
            let bam_path = next_arg(&mut it, "bam");
            let bam = resolve_alignment_for_parity(Path::new(&bam_path))?;
            let interval_cli = next_arg(&mut it, "interval_cli");
            let padding = parse_padding(&mut it)?;
            let (_, st) = run_walker(
                ref_fa,
                &bam,
                interval_cli.as_str(),
                padding,
                false,
                None,
                false,
            )?;
            print_apply_summary(&st);
        }
        "locus-pileup" => {
            let ref_path = next_arg(&mut it, "ref.fa");
            let ref_fa = Path::new(&ref_path);
            let bam_path = next_arg(&mut it, "bam");
            let bam = resolve_alignment_for_parity(Path::new(&bam_path))?;
            let interval_cli = next_arg(&mut it, "interval_cli");
            let padding = parse_padding(&mut it)?;
            let filters = ReadFilterParams::gatk_standard_hc();
            let mut stdout = io::stdout().lock();
            dump_locus_pileup_tsv(
                ref_fa,
                &bam,
                interval_cli.as_str(),
                padding,
                &mut stdout,
                &filters,
            )
            .map_err(|e| format!("locus-pileup: {e}"))?;
        }
        "locus-pileup-detail" => {
            let ref_path = next_arg(&mut it, "ref.fa");
            let ref_fa = Path::new(&ref_path);
            let bam_path = next_arg(&mut it, "bam");
            let bam = resolve_alignment_for_parity(Path::new(&bam_path))?;
            let interval_cli = next_arg(&mut it, "interval_cli");
            let padding = parse_padding(&mut it)?;
            let filters = ReadFilterParams::gatk_standard_hc();
            let mut stdout = io::stdout().lock();
            dump_locus_pileup_detail_tsv(
                ref_fa,
                &bam,
                interval_cli.as_str(),
                padding,
                &mut stdout,
                &filters,
            )
            .map_err(|e| format!("locus-pileup-detail: {e}"))?;
        }
        "genotype-likelihoods" => {
            let ref_path = next_arg(&mut it, "ref.fa");
            let ref_fa = Path::new(&ref_path);
            let bam_path = next_arg(&mut it, "bam");
            let bam = resolve_alignment_for_parity(Path::new(&bam_path))?;
            let interval_cli = next_arg(&mut it, "interval_cli");
            let filters = ReadFilterParams::gatk_standard_hc();
            let mut stdout = io::stdout().lock();
            let padding = parse_padding(&mut it)?;
            dump_genotype_likelihood_activity_tsv(
                ref_fa,
                &bam,
                interval_cli.as_str(),
                padding,
                &mut stdout,
                &filters,
            )
            .map_err(|e| format!("genotype-likelihoods: {e}"))?;
        }
        "raw-activity" | "smoothed-activity" | "active-locus" => {
            let ref_path = next_arg(&mut it, "ref.fa");
            let ref_fa = Path::new(&ref_path);
            let bam_path = next_arg(&mut it, "bam");
            let bam = resolve_alignment_for_parity(Path::new(&bam_path))?;
            let interval_cli = next_arg(&mut it, "interval_cli");
            let padding = parse_padding(&mut it)?;
            let filters = ReadFilterParams::gatk_standard_hc();
            let mut stdout = io::stdout().lock();
            match cmd.as_str() {
                "raw-activity" => dump_raw_activity_profile_tsv_with_force_calling(
                    ref_fa,
                    &bam,
                    interval_cli.as_str(),
                    padding,
                    &mut stdout,
                    &filters,
                    None,
                ),
                "smoothed-activity" => dump_smoothed_activity_profile_tsv(
                    ref_fa,
                    &bam,
                    interval_cli.as_str(),
                    padding,
                    &mut stdout,
                    &filters,
                ),
                "active-locus" => dump_active_locus_tsv(
                    ref_fa,
                    &bam,
                    interval_cli.as_str(),
                    padding,
                    &mut stdout,
                    &filters,
                ),
                _ => unreachable!(),
            }
            .map_err(|e| format!("{cmd}: {e}"))?;
        }
        "raw-activity-force" => {
            let ref_path = next_arg(&mut it, "ref.fa");
            let ref_fa = Path::new(&ref_path);
            let bam_path = next_arg(&mut it, "bam");
            let bam = resolve_alignment_for_parity(Path::new(&bam_path))?;
            let interval_cli = next_arg(&mut it, "interval_cli");
            let vcf_path = next_arg(&mut it, "alleles.vcf");
            let filters = ReadFilterParams::gatk_standard_hc();
            let mut stdout = io::stdout().lock();
            let mut sources = FeatureDataSources::default();
            sources
                .load_vcf_source(GATK_HC_ALLELES_FEATURE_SOURCE, Path::new(&vcf_path))
                .map_err(|e| format!("raw-activity-force: load VCF: {e}"))?;
            let fc = ForceCallingAllelesDump {
                sources: &sources,
                force_call_filtered: false,
            };
            let padding = parse_padding(&mut it)?;
            dump_raw_activity_profile_tsv_with_force_calling(
                ref_fa,
                &bam,
                interval_cli.as_str(),
                padding,
                &mut stdout,
                &filters,
                Some(fc),
            )
            .map_err(|e| format!("raw-activity-force: {e}"))?;
        }
        "read-filters" => {
            let bam_path = next_arg(&mut it, "bam");
            let bam = resolve_alignment_for_parity(Path::new(&bam_path))?;
            let mut stdout = io::stdout().lock();
            dump_hc_read_filter_tsv(&bam, &mut stdout).map_err(|e| format!("read-filters: {e}"))?;
        }
        "read-shard-pipeline" => {
            let bam_path = next_arg(&mut it, "bam");
            let bam = resolve_alignment_for_parity(Path::new(&bam_path))?;
            let mut stdout = io::stdout().lock();
            dump_read_shard_pipeline_tsv(&bam, false, &mut stdout)
                .map_err(|e| format!("read-shard-pipeline: {e}"))?;
        }
        "read-shard-pipeline-dragen" => {
            let bam_path = next_arg(&mut it, "bam");
            let bam = resolve_alignment_for_parity(Path::new(&bam_path))?;
            let mut stdout = io::stdout().lock();
            dump_read_shard_pipeline_tsv(&bam, true, &mut stdout)
                .map_err(|e| format!("read-shard-pipeline-dragen: {e}"))?;
        }
        "read-pre-softclip" => {
            let bam_path = next_arg(&mut it, "bam");
            let bam = resolve_alignment_for_parity(Path::new(&bam_path))?;
            let dont_use = it.next().map(|s| s == "1").unwrap_or(false);
            let override_frag = it.next().map(|s| s == "1").unwrap_or(false);
            let mut stdout = io::stdout().lock();
            dump_read_pre_softclip_tsv(&bam, dont_use, override_frag, &mut stdout)
                .map_err(|e| format!("read-pre-softclip: {e}"))?;
        }
        "read-pre-len" => {
            let bam_path = next_arg(&mut it, "bam");
            let bam = resolve_alignment_for_parity(Path::new(&bam_path))?;
            let mut stdout = io::stdout().lock();
            dump_read_pre_len_tsv(&bam, &mut stdout).map_err(|e| format!("read-pre-len: {e}"))?;
        }
        "read-pre-mq" => {
            let bam_path = next_arg(&mut it, "bam");
            let bam = resolve_alignment_for_parity(Path::new(&bam_path))?;
            let mq_threshold: u8 = it
                .next()
                .unwrap_or_else(|| "20".to_string())
                .parse()
                .map_err(|_| "read-pre-mq: mq_threshold must be u8".to_string())?;
            let mut stdout = io::stdout().lock();
            dump_read_pre_mq_tsv(&bam, mq_threshold, &mut stdout)
                .map_err(|e| format!("read-pre-mq: {e}"))?;
        }
        "read-pre-overlap" => {
            let bam_path = next_arg(&mut it, "bam");
            let bam = resolve_alignment_for_parity(Path::new(&bam_path))?;
            let mut stdout = io::stdout().lock();
            dump_read_pre_overlap_tsv(&bam, &mut stdout)
                .map_err(|e| format!("read-pre-overlap: {e}"))?;
        }
        "downsample-positional" => {
            let bam_path = next_arg(&mut it, "bam");
            let bam = resolve_alignment_for_parity(Path::new(&bam_path))?;
            let cap: u32 = next_arg(&mut it, "cap")
                .parse()
                .map_err(|_| "cap: invalid u32".to_string())?;
            let mode = it.next().unwrap_or_else(|| "non-random".to_string());
            let non_random = match mode.as_str() {
                "non-random" | "nonrandom" => true,
                "random" => false,
                other => {
                    return Err(format!(
                        "downsample-positional: unknown mode {other:?} (use non-random or random)"
                    ));
                }
            };
            let mut stdout = io::stdout().lock();
            dump_positional_downsample_summary_tsv(&bam, cap, non_random, &mut stdout)
                .map_err(|e| format!("downsample-positional: {e}"))?;
        }
        "allele-biased-target-counts" => {
            let counts_csv = next_arg(&mut it, "allele_counts_csv");
            let num_remove: usize = next_arg(&mut it, "num_reads_to_remove")
                .parse()
                .map_err(|_| "num_reads_to_remove: invalid usize".to_string())?;
            let mut stdout = io::stdout().lock();
            dump_target_allele_counts_tsv(counts_csv.as_str(), num_remove, &mut stdout)
                .map_err(|e| format!("allele-biased-target-counts: {e}"))?;
        }
        "allele-biased-evidence" => {
            let ref_path = next_arg(&mut it, "ref.fa");
            let ref_fa = Path::new(&ref_path);
            let bam_path = next_arg(&mut it, "bam");
            let bam = resolve_alignment_for_parity(Path::new(&bam_path))?;
            let contig = next_arg(&mut it, "contig");
            let pos1: u64 = next_arg(&mut it, "pos1")
                .parse()
                .map_err(|_| "pos1: invalid u64".to_string())?;
            let contamination: f64 = next_arg(&mut it, "contamination_fraction")
                .parse()
                .map_err(|_| "contamination_fraction: invalid f64".to_string())?;
            let filters = ReadFilterParams::default();
            let mut stdout = io::stdout().lock();
            dump_allele_biased_evidence_locus_tsv(
                ref_fa,
                &bam,
                contig.as_str(),
                pos1,
                contamination,
                &mut stdout,
                &filters,
            )
            .map_err(|e| format!("allele-biased-evidence: {e}"))?;
        }
        "raw-activity-contam" => {
            let ref_path = next_arg(&mut it, "ref.fa");
            let ref_fa = Path::new(&ref_path);
            let bam_path = next_arg(&mut it, "bam");
            let bam = resolve_alignment_for_parity(Path::new(&bam_path))?;
            let interval_cli = next_arg(&mut it, "interval_cli");
            let contamination: f64 = next_arg(&mut it, "contamination_fraction")
                .parse()
                .map_err(|_| "contamination_fraction: invalid f64".to_string())?;
            let filters = ReadFilterParams::gatk_standard_hc();
            let mut stdout = io::stdout().lock();
            dump_raw_activity_profile_tsv_with_contamination(
                ref_fa,
                &bam,
                interval_cli.as_str(),
                contamination,
                &mut stdout,
                &filters,
            )
            .map_err(|e| format!("raw-activity-contam: {e}"))?;
        }
        "soft-clip-mean" => {
            let ref_path = next_arg(&mut it, "ref.fa");
            let ref_fa = Path::new(&ref_path);
            let bam_path = next_arg(&mut it, "bam");
            let bam = resolve_alignment_for_parity(Path::new(&bam_path))?;
            let interval_cli = next_arg(&mut it, "interval_cli");
            let filters = ReadFilterParams::default();
            let mut stdout = io::stdout().lock();
            dump_hq_soft_clip_mean_tsv(ref_fa, &bam, interval_cli.as_str(), &mut stdout, &filters)
                .map_err(|e| format!("soft-clip-mean: {e}"))?;
        }
        "assembly-graph" => {
            let reads_path = next_arg(&mut it, "reads.tsv");
            let reads = Path::new(&reads_path);
            let kmer_size: usize = it
                .next()
                .map(|s| {
                    s.parse()
                        .map_err(|_| "kmer_size: invalid usize".to_string())
                })
                .transpose()?
                .unwrap_or(3);
            let min_qual: u8 = it
                .next()
                .map(|s| {
                    s.parse()
                        .map_err(|_| "min_base_quality: invalid u8".to_string())
                })
                .transpose()?
                .unwrap_or(10);
            let mut stdout = io::stdout().lock();
            dump_assembly_graph_edges_tsv(reads, kmer_size, min_qual, &mut stdout)
                .map_err(|e| format!("assembly-graph: {e}"))?;
        }
        "assembly-graph-multi" => {
            let reads_path = next_arg(&mut it, "reads.tsv");
            let reads = Path::new(&reads_path);
            let kmer_list = next_arg(&mut it, "kmer_sizes_csv");
            let kmer_sizes: Vec<usize> = kmer_list
                .split(',')
                .filter(|s| !s.is_empty())
                .map(|s| {
                    s.parse()
                        .map_err(|_| format!("kmer_sizes: invalid usize in {s}"))
                })
                .collect::<Result<Vec<_>, _>>()?;
            if kmer_sizes.is_empty() {
                return Err("kmer_sizes: need at least one value".to_string());
            }
            let min_qual: u8 = it
                .next()
                .map(|s| {
                    s.parse()
                        .map_err(|_| "min_base_quality: invalid u8".to_string())
                })
                .transpose()?
                .unwrap_or(10);
            let mut stdout = io::stdout().lock();
            dump_assembly_graph_multi_kmer_edges_tsv(reads, &kmer_sizes, min_qual, &mut stdout)
                .map_err(|e| format!("assembly-graph-multi: {e}"))?;
        }
        "assembly-graph-summary" => {
            let reads_path = next_arg(&mut it, "reads.tsv");
            let reads = Path::new(&reads_path);
            let kmer_size: usize = next_arg(&mut it, "kmer_size")
                .parse()
                .map_err(|_| "kmer_size: invalid usize".to_string())?;
            let min_qual: u8 = next_arg(&mut it, "min_base_quality")
                .parse()
                .map_err(|_| "min_base_quality: invalid u8".to_string())?;
            let min_prune: u32 = next_arg(&mut it, "min_prune")
                .parse()
                .map_err(|_| "min_prune: invalid u32".to_string())?;
            let adaptive = match next_arg(&mut it, "adaptive").as_str() {
                "0" | "false" => false,
                "1" | "true" => true,
                s => return Err(format!("adaptive: expected 0|1|true|false, got {s}")),
            };
            let mut pruning = AssemblyGraphPruningParams::gatk_haplotype_caller_defaults();
            pruning.min_prune_factor = min_prune;
            pruning.use_adaptive_pruning = adaptive;
            let mut stdout = io::stdout().lock();
            dump_assembly_graph_pruned_summary_tsv(
                reads,
                kmer_size,
                min_qual,
                &pruning,
                &mut stdout,
            )
            .map_err(|e| format!("assembly-graph-summary: {e}"))?;
        }
        "assembly-graph-dangling-summary" => {
            let ref_tsv = next_arg(&mut it, "ref.tsv");
            let reads_tsv = next_arg(&mut it, "reads.tsv");
            let ref_path = Path::new(&ref_tsv);
            let reads_path = Path::new(&reads_tsv);
            let kmer_size: usize = next_arg(&mut it, "kmer_size")
                .parse()
                .map_err(|_| "kmer_size: invalid usize".to_string())?;
            let min_qual: u8 = next_arg(&mut it, "min_base_quality")
                .parse()
                .map_err(|_| "min_base_quality: invalid u8".to_string())?;
            let min_prune: u32 = next_arg(&mut it, "min_prune")
                .parse()
                .map_err(|_| "min_prune: invalid u32".to_string())?;
            let min_dangling: usize = next_arg(&mut it, "min_dangling")
                .parse()
                .map_err(|_| "min_dangling: invalid usize".to_string())?;
            let recover_heads = match next_arg(&mut it, "recover_heads").as_str() {
                "0" | "false" => false,
                "1" | "true" => true,
                s => return Err(format!("recover_heads: expected 0|1|true|false, got {s}")),
            };
            let mut stdout = io::stdout().lock();
            dump_assembly_graph_dangling_summary_tsv(
                ref_path,
                reads_path,
                kmer_size,
                min_qual,
                min_prune,
                min_dangling,
                recover_heads,
                &mut stdout,
            )
            .map_err(|e| format!("assembly-graph-dangling-summary: {e}"))?;
        }
        "assembly-graph-non-unique-summary" => {
            let ref_tsv = next_arg(&mut it, "ref.tsv|->");
            let reads_tsv = next_arg(&mut it, "reads.tsv");
            let reads_path = Path::new(&reads_tsv);
            let ref_path = if ref_tsv == "-" {
                None
            } else {
                Some(Path::new(&ref_tsv))
            };
            let kmer_size: usize = next_arg(&mut it, "kmer_size")
                .parse()
                .map_err(|_| "kmer_size: invalid usize".to_string())?;
            let min_qual: u8 = next_arg(&mut it, "min_base_quality")
                .parse()
                .map_err(|_| "min_base_quality: invalid u8".to_string())?;
            let mut stdout = io::stdout().lock();
            dump_assembly_graph_non_unique_summary_tsv(
                ref_path,
                reads_path,
                kmer_size,
                min_qual,
                &mut stdout,
            )
            .map_err(|e| format!("assembly-graph-non-unique-summary: {e}"))?;
        }
        "assembly-haplotype-cigars" => {
            let ref_tsv = next_arg(&mut it, "ref.tsv");
            let haps_tsv = next_arg(&mut it, "haplotypes.tsv");
            let mut stdout = io::stdout().lock();
            dump_assembly_haplotype_cigars_tsv(
                Path::new(&ref_tsv),
                Path::new(&haps_tsv),
                &mut stdout,
            )
            .map_err(|e| format!("assembly-haplotype-cigars: {e}"))?;
        }
        "assembly-haplotypes" => {
            let ref_tsv = next_arg(&mut it, "ref.tsv");
            let reads_tsv = next_arg(&mut it, "reads.tsv");
            let kmer_size: usize = next_arg(&mut it, "kmer_size")
                .parse()
                .map_err(|_| "kmer_size: invalid usize".to_string())?;
            let min_qual: u8 = next_arg(&mut it, "min_base_quality")
                .parse()
                .map_err(|_| "min_base_quality: invalid u8".to_string())?;
            let min_prune: u32 = next_arg(&mut it, "min_prune")
                .parse()
                .map_err(|_| "min_prune: invalid u32".to_string())?;
            let min_dangling: usize = next_arg(&mut it, "min_dangling")
                .parse()
                .map_err(|_| "min_dangling: invalid usize".to_string())?;
            let recover_heads = match next_arg(&mut it, "recover_heads").as_str() {
                "0" | "false" => false,
                "1" | "true" => true,
                s => return Err(format!("recover_heads: expected 0|1|true|false, got {s}")),
            };
            let mut stdout = io::stdout().lock();
            dump_assembly_haplotypes_tsv(
                Path::new(&ref_tsv),
                Path::new(&reads_tsv),
                kmer_size,
                min_qual,
                min_prune,
                min_dangling,
                recover_heads,
                &mut stdout,
            )
            .map_err(|e| format!("assembly-haplotypes: {e}"))?;
        }
        "assembly-kbest-paths" => {
            let ref_tsv = next_arg(&mut it, "ref.tsv");
            let reads_tsv = next_arg(&mut it, "reads.tsv");
            let kmer_size: usize = next_arg(&mut it, "kmer_size")
                .parse()
                .map_err(|_| "kmer_size: invalid usize".to_string())?;
            let min_qual: u8 = next_arg(&mut it, "min_base_quality")
                .parse()
                .map_err(|_| "min_base_quality: invalid u8".to_string())?;
            let min_prune: u32 = next_arg(&mut it, "min_prune")
                .parse()
                .map_err(|_| "min_prune: invalid u32".to_string())?;
            let min_dangling: usize = next_arg(&mut it, "min_dangling")
                .parse()
                .map_err(|_| "min_dangling: invalid usize".to_string())?;
            let recover_heads = match next_arg(&mut it, "recover_heads").as_str() {
                "0" | "false" => false,
                "1" | "true" => true,
                s => return Err(format!("recover_heads: expected 0|1|true|false, got {s}")),
            };
            let max_haps: usize = next_arg(&mut it, "max_haplotypes")
                .parse()
                .map_err(|_| "max_haplotypes: invalid usize".to_string())?;
            let mut stdout = io::stdout().lock();
            dump_assembly_kbest_paths_tsv(
                Path::new(&ref_tsv),
                Path::new(&reads_tsv),
                kmer_size,
                min_qual,
                min_prune,
                min_dangling,
                recover_heads,
                max_haps,
                &mut stdout,
            )
            .map_err(|e| format!("assembly-kbest-paths: {e}"))?;
        }
        "assembly-junction-haplotypes" => {
            let ref_tsv = next_arg(&mut it, "ref.tsv");
            let reads_tsv = next_arg(&mut it, "reads.tsv");
            let kmer_size: usize = next_arg(&mut it, "kmer_size")
                .parse()
                .map_err(|_| "kmer_size: invalid usize".to_string())?;
            let min_qual: u8 = next_arg(&mut it, "min_base_quality")
                .parse()
                .map_err(|_| "min_base_quality: invalid u8".to_string())?;
            let recover_edges = match next_arg(&mut it, "recover_edges").as_str() {
                "0" | "false" => false,
                "1" | "true" => true,
                s => return Err(format!("recover_edges: expected 0|1|true|false, got {s}")),
            };
            let max_haps: usize = next_arg(&mut it, "max_haplotypes")
                .parse()
                .map_err(|_| "max_haplotypes: invalid usize".to_string())?;
            let mut stdout = io::stdout().lock();
            dump_assembly_junction_haplotypes_tsv(
                Path::new(&ref_tsv),
                Path::new(&reads_tsv),
                kmer_size,
                min_qual,
                recover_edges,
                max_haps,
                &mut stdout,
            )
            .map_err(|e| format!("assembly-junction-haplotypes: {e}"))?;
        }
        "assembly-seqgraph-summary" => {
            let ref_tsv = next_arg(&mut it, "ref.tsv");
            let reads_tsv = next_arg(&mut it, "reads.tsv");
            let kmer_size: usize = next_arg(&mut it, "kmer_size")
                .parse()
                .map_err(|_| "kmer_size: invalid usize".to_string())?;
            let min_qual: u8 = next_arg(&mut it, "min_base_quality")
                .parse()
                .map_err(|_| "min_base_quality: invalid u8".to_string())?;
            let min_prune: u32 = next_arg(&mut it, "min_prune")
                .parse()
                .map_err(|_| "min_prune: invalid u32".to_string())?;
            let min_dangling: usize = next_arg(&mut it, "min_dangling")
                .parse()
                .map_err(|_| "min_dangling: invalid usize".to_string())?;
            let recover_heads = match next_arg(&mut it, "recover_heads").as_str() {
                "0" | "false" => false,
                "1" | "true" => true,
                s => return Err(format!("recover_heads: expected 0|1|true|false, got {s}")),
            };
            let mut stdout = io::stdout().lock();
            dump_assembly_seqgraph_summary_tsv(
                Path::new(&ref_tsv),
                Path::new(&reads_tsv),
                kmer_size,
                min_qual,
                min_prune,
                min_dangling,
                recover_heads,
                &mut stdout,
            )
            .map_err(|e| format!("assembly-seqgraph-summary: {e}"))?;
        }
        "read-error-correction" => {
            let reads_tsv = next_arg(&mut it, "reads.tsv");
            let log_odds: f64 = next_arg(&mut it, "log_odds_threshold")
                .parse()
                .map_err(|_| "log_odds_threshold: invalid f64".to_string())?;
            let mut stdout = io::stdout().lock();
            dump_read_error_correction_tsv(Path::new(&reads_tsv), log_odds, &mut stdout)
                .map_err(|e| format!("read-error-correction: {e}"))?;
        }
        "assembly-assemble" => {
            let ref_tsv = next_arg(&mut it, "ref.tsv");
            let reads_tsv = next_arg(&mut it, "reads.tsv");
            let mut stdout = io::stdout().lock();
            dump_assembly_assemble_tsv(Path::new(&ref_tsv), Path::new(&reads_tsv), &mut stdout)
                .map_err(|e| format!("assembly-assemble: {e}"))?;
        }
        "assembler-args" => {
            let mut stdout = io::stdout().lock();
            dump_assembler_args_tsv(&mut stdout).map_err(|e| format!("assembler-args: {e}"))?;
        }
        "assembly-graph-low-quality" => {
            let ref_tsv = next_arg(&mut it, "ref.tsv");
            let reads_tsv = next_arg(&mut it, "reads.tsv");
            let kmer_size: usize = next_arg(&mut it, "kmer_size")
                .parse()
                .map_err(|_| "kmer_size: invalid usize".to_string())?;
            let mut stdout = io::stdout().lock();
            dump_assembly_graph_low_quality_tsv(
                Path::new(&ref_tsv),
                Path::new(&reads_tsv),
                kmer_size,
                &mut stdout,
            )
            .map_err(|e| format!("assembly-graph-low-quality: {e}"))?;
        }
        "likelihood-engine-config" => {
            let mut stdout = io::stdout().lock();
            dump_likelihood_engine_config_tsv(&mut stdout)
                .map_err(|e| format!("likelihood-engine-config: {e}"))?;
        }
        "pcr-error-model" => {
            let mut stdout = io::stdout().lock();
            dump_pcr_error_model_tsv(&mut stdout).map_err(|e| format!("pcr-error-model: {e}"))?;
        }
        "likelihood-pcr-read" => {
            let cases_path = next_arg(&mut it, "cases.tsv");
            let mut stdout = io::stdout().lock();
            dump_likelihood_pcr_read_tsv(Path::new(&cases_path), &mut stdout)
                .map_err(|e| format!("likelihood-pcr-read: {e}"))?;
        }
        "assembly-region-haplotypes" => {
            let ref_path = next_arg(&mut it, "ref.fa");
            let bam_path = next_arg(&mut it, "bam");
            let interval_cli = next_arg(&mut it, "interval_cli");
            let (padding, target) = parse_padding_and_region_target(&mut it)?;
            let bam = resolve_alignment_for_parity(Path::new(&bam_path))?;
            let mut stdout = io::stdout().lock();
            dump_assembly_region_haplotypes_tsv(
                Path::new(&ref_path),
                &bam,
                interval_cli.as_str(),
                padding,
                target,
                &mut stdout,
            )
            .map_err(|e| format!("assembly-region-haplotypes: {e}"))?;
        }
        "assembly-region-kmer-probe" => {
            let ref_path = next_arg(&mut it, "ref.fa");
            let bam_path = next_arg(&mut it, "bam");
            let interval_cli = next_arg(&mut it, "interval_cli");
            let padding = parse_padding(&mut it)?;
            let bam = resolve_alignment_for_parity(Path::new(&bam_path))?;
            let mut stdout = io::stdout().lock();
            dump_assembly_region_kmer_probe_tsv(
                Path::new(&ref_path),
                &bam,
                interval_cli.as_str(),
                padding,
                &mut stdout,
            )
            .map_err(|e| format!("assembly-region-kmer-probe: {e}"))?;
        }
        "assembly-region-assembly-stages" => {
            let ref_path = next_arg(&mut it, "ref.fa");
            let bam_path = next_arg(&mut it, "bam");
            let interval_cli = next_arg(&mut it, "interval_cli");
            let padding = parse_padding(&mut it)?;
            let bam = resolve_alignment_for_parity(Path::new(&bam_path))?;
            let mut stdout = io::stdout().lock();
            dump_assembly_region_assembly_stages_tsv(
                Path::new(&ref_path),
                &bam,
                interval_cli.as_str(),
                padding,
                &mut stdout,
            )
            .map_err(|e| format!("assembly-region-assembly-stages: {e}"))?;
        }
        "assembly-region-assembly-stages-finalize" => {
            let ref_path = next_arg(&mut it, "ref.fa");
            let bam_path = next_arg(&mut it, "bam");
            let interval_cli = next_arg(&mut it, "interval_cli");
            let padding = parse_padding(&mut it)?;
            let bam = resolve_alignment_for_parity(Path::new(&bam_path))?;
            let mut stdout = io::stdout().lock();
            dump_assembly_region_assembly_stages_finalize_tsv(
                Path::new(&ref_path),
                &bam,
                interval_cli.as_str(),
                padding,
                &mut stdout,
            )
            .map_err(|e| format!("assembly-region-assembly-stages-finalize: {e}"))?;
        }
        "assembly-region-finalize-reads" => {
            let ref_path = next_arg(&mut it, "ref.fa");
            let bam_path = next_arg(&mut it, "bam");
            let interval_cli = next_arg(&mut it, "interval_cli");
            let padding = parse_padding(&mut it)?;
            let bam = resolve_alignment_for_parity(Path::new(&bam_path))?;
            let mut stdout = io::stdout().lock();
            dump_assembly_region_finalize_reads_tsv(
                Path::new(&ref_path),
                &bam,
                interval_cli.as_str(),
                padding,
                &mut stdout,
            )
            .map_err(|e| format!("assembly-region-finalize-reads: {e}"))?;
        }
        "assembly-region-kbest-paths" => {
            let ref_path = next_arg(&mut it, "ref.fa");
            let bam_path = next_arg(&mut it, "bam");
            let interval_cli = next_arg(&mut it, "interval_cli");
            let padding = parse_padding(&mut it)?;
            let bam = resolve_alignment_for_parity(Path::new(&bam_path))?;
            let mut stdout = io::stdout().lock();
            dump_assembly_region_kbest_paths_tsv(
                Path::new(&ref_path),
                &bam,
                interval_cli.as_str(),
                padding,
                &mut stdout,
            )
            .map_err(|e| format!("assembly-region-kbest-paths: {e}"))?;
        }
        "assembly-region-pairhmm-likelihoods" => {
            let ref_path = next_arg(&mut it, "ref.fa");
            let bam_path = next_arg(&mut it, "bam");
            let interval_cli = next_arg(&mut it, "interval_cli");
            let (padding, target) = parse_padding_and_region_target(&mut it)?;
            let bam = resolve_alignment_for_parity(Path::new(&bam_path))?;
            let mut stdout = io::stdout().lock();
            dump_assembly_region_pairhmm_likelihoods_tsv(
                Path::new(&ref_path),
                &bam,
                interval_cli.as_str(),
                padding,
                target,
                &mut stdout,
            )
            .map_err(|e| format!("assembly-region-pairhmm-likelihoods: {e}"))?;
        }
        "pairhmm-native-likelihoods" => {
            let cases_path = next_arg(&mut it, "cases.tsv");
            let mut stdout = io::stdout().lock();
            dump_pairhmm_native_likelihoods_tsv(Path::new(&cases_path), &mut stdout)
                .map_err(|e| format!("pairhmm-native-likelihoods: {e}"))?;
        }
        "assembly-haplotypes-cap" => {
            let ref_tsv = next_arg(&mut it, "ref.tsv");
            let reads_tsv = next_arg(&mut it, "reads.tsv");
            let kmer_size: usize = next_arg(&mut it, "kmer_size")
                .parse()
                .map_err(|_| "kmer_size: invalid usize".to_string())?;
            let min_qual: u8 = next_arg(&mut it, "min_base_quality")
                .parse()
                .map_err(|_| "min_base_quality: invalid u8".to_string())?;
            let min_prune: u32 = next_arg(&mut it, "min_prune")
                .parse()
                .map_err(|_| "min_prune: invalid u32".to_string())?;
            let min_dangling: usize = next_arg(&mut it, "min_dangling")
                .parse()
                .map_err(|_| "min_dangling: invalid usize".to_string())?;
            let recover_heads = match next_arg(&mut it, "recover_heads").as_str() {
                "0" | "false" => false,
                "1" | "true" => true,
                s => return Err(format!("recover_heads: expected 0|1|true|false, got {s}")),
            };
            let max_haps: usize = next_arg(&mut it, "max_haplotypes")
                .parse()
                .map_err(|_| "max_haplotypes: invalid usize".to_string())?;
            let mut stdout = io::stdout().lock();
            dump_assembly_haplotypes_cap_tsv(
                Path::new(&ref_tsv),
                Path::new(&reads_tsv),
                kmer_size,
                min_qual,
                min_prune,
                min_dangling,
                recover_heads,
                max_haps,
                &mut stdout,
            )
            .map_err(|e| format!("assembly-haplotypes-cap: {e}"))?;
        }
        "assembly-haplotypes-production" => {
            let ref_tsv = next_arg(&mut it, "ref.tsv");
            let reads_tsv = next_arg(&mut it, "reads.tsv");
            let kmer_size: usize = next_arg(&mut it, "kmer_size")
                .parse()
                .map_err(|_| "kmer_size: invalid usize".to_string())?;
            let min_qual: u8 = next_arg(&mut it, "min_base_quality")
                .parse()
                .map_err(|_| "min_base_quality: invalid u8".to_string())?;
            let min_prune: u32 = next_arg(&mut it, "min_prune")
                .parse()
                .map_err(|_| "min_prune: invalid u32".to_string())?;
            let min_dangling: usize = next_arg(&mut it, "min_dangling")
                .parse()
                .map_err(|_| "min_dangling: invalid usize".to_string())?;
            let recover_heads = match next_arg(&mut it, "recover_heads").as_str() {
                "0" | "false" => false,
                "1" | "true" => true,
                s => return Err(format!("recover_heads: expected 0|1|true|false, got {s}")),
            };
            let mut stdout = io::stdout().lock();
            dump_assembly_haplotypes_production_tsv(
                Path::new(&ref_tsv),
                Path::new(&reads_tsv),
                kmer_size,
                min_qual,
                min_prune,
                min_dangling,
                recover_heads,
                &mut stdout,
            )
            .map_err(|e| format!("assembly-haplotypes-production: {e}"))?;
        }
        "pairhmm-likelihoods" => {
            let cases_tsv = next_arg(&mut it, "cases.tsv");
            let mut stdout = io::stdout().lock();
            dump_pairhmm_likelihoods_tsv(Path::new(&cases_tsv), &mut stdout)
                .map_err(|e| format!("pairhmm-likelihoods: {e}"))?;
        }
        "pairhmm-bq-cap" => {
            let cases_tsv = next_arg(&mut it, "cases.tsv");
            let mut stdout = io::stdout().lock();
            dump_pairhmm_bq_cap_tsv(Path::new(&cases_tsv), &mut stdout)
                .map_err(|e| format!("pairhmm-bq-cap: {e}"))?;
        }
        "pairhmm-haplotype-filter" => {
            let cases_tsv = next_arg(&mut it, "cases.tsv");
            let mut stdout = io::stdout().lock();
            dump_pairhmm_haplotype_filter_tsv(Path::new(&cases_tsv), &mut stdout)
                .map_err(|e| format!("pairhmm-haplotype-filter: {e}"))?;
        }
        "genotyping-aggregate" => {
            let cases_tsv = next_arg(&mut it, "pairhmm_cases.tsv");
            let mut stdout = io::stdout().lock();
            dump_genotyping_aggregate_tsv(Path::new(&cases_tsv), &mut stdout)
                .map_err(|e| format!("genotyping-aggregate: {e}"))?;
        }
        "genotype-format" => {
            let fixture = next_arg(&mut it, "fixture.tsv");
            let mut stdout = io::stdout().lock();
            dump_genotype_format_tsv(Path::new(&fixture), &mut stdout)
                .map_err(|e| format!("genotype-format: {e}"))?;
        }
        "assembly-region-genotype" => {
            let ref_path = next_arg(&mut it, "ref.fa");
            let bam_path = next_arg(&mut it, "bam");
            let interval_cli = next_arg(&mut it, "interval_cli");
            let (padding, target) = parse_padding_and_region_target(&mut it)?;
            let bam = resolve_alignment_for_parity(Path::new(&bam_path))?;
            let mut stdout = io::stdout().lock();
            dump_assembly_region_genotype_tsv(
                Path::new(&ref_path),
                &bam,
                interval_cli.as_str(),
                padding,
                target,
                &mut stdout,
            )
            .map_err(|e| format!("assembly-region-genotype: {e}"))?;
        }
        "assembly-region-genotype-subset" => {
            let ref_path = next_arg(&mut it, "ref.fa");
            let bam_path = next_arg(&mut it, "bam");
            let interval_cli = next_arg(&mut it, "interval_cli");
            let (padding, target) = parse_padding_and_region_target(&mut it)?;
            let assembly_profile = next_arg(&mut it, "assembly_profile");
            let max_allele_count: usize = next_arg(&mut it, "max_allele_count")
                .parse()
                .map_err(|_| "max_allele_count: invalid usize".to_string())?;
            let bam = resolve_alignment_for_parity(Path::new(&bam_path))?;
            let mut stdout = io::stdout().lock();
            dump_assembly_region_genotype_subset_tsv(
                Path::new(&ref_path),
                &bam,
                interval_cli.as_str(),
                padding,
                target,
                max_allele_count,
                assembly_profile.as_str(),
                &mut stdout,
            )
            .map_err(|e| format!("assembly-region-genotype-subset: {e}"))?;
        }
        "gvcf-merge-ref-confidence" => {
            let case_id = next_arg(&mut it, "case_id");
            let mut stdout = io::stdout().lock();
            dump_ref_confidence_merge_case(case_id.as_str(), &mut stdout)
                .map_err(|e| format!("gvcf-merge-ref-confidence: {e}"))?;
        }
        "reference-confidence-locus" => {
            let ref_path = next_arg(&mut it, "ref.fa");
            let bam_path = next_arg(&mut it, "bam");
            let interval_cli = next_arg(&mut it, "interval_cli");
            let padding: u64 = it
                .next()
                .map(|s| s.parse().map_err(|_| "padding: invalid u64".to_string()))
                .transpose()?
                .unwrap_or(0);
            let bam = resolve_alignment_for_parity(Path::new(&bam_path))?;
            let filters = ReadFilterParams::gatk_standard_hc();
            let mut stdout = io::stdout().lock();
            dump_reference_confidence_locus_tsv(
                Path::new(&ref_path),
                &bam,
                interval_cli.as_str(),
                padding,
                &mut stdout,
                &filters,
            )
            .map_err(|e| format!("reference-confidence-locus: {e}"))?;
        }
        "call-region-active-rcm-loci" => {
            let ref_path = next_arg(&mut it, "ref.fa");
            let bam_path = next_arg(&mut it, "bam");
            let interval_cli = next_arg(&mut it, "interval_cli");
            let padding: u64 = it
                .next()
                .map(|s| s.parse().map_err(|_| "padding: invalid u64".to_string()))
                .transpose()?
                .unwrap_or(0);
            let bam = resolve_alignment_for_parity(Path::new(&bam_path))?;
            let mut stdout = io::stdout().lock();
            dump_call_region_active_rcm_loci_tsv(
                Path::new(&ref_path),
                &bam,
                interval_cli.as_str(),
                padding,
                &mut stdout,
            )
            .map_err(|e| format!("call-region-active-rcm-loci: {e}"))?;
        }
        "inactive-reference-model" => {
            let ref_path = next_arg(&mut it, "ref.fa");
            let bam_path = next_arg(&mut it, "bam");
            let interval_cli = next_arg(&mut it, "interval_cli");
            let padding: u64 = it
                .next()
                .map(|s| s.parse().map_err(|_| "padding: invalid u64".to_string()))
                .transpose()?
                .unwrap_or(0);
            let bam = resolve_alignment_for_parity(Path::new(&bam_path))?;
            let mut stdout = io::stdout().lock();
            dump_inactive_reference_model_tsv(
                Path::new(&ref_path),
                &bam,
                interval_cli.as_str(),
                padding,
                gatk_haplotypecaller::genotyping::EmitMode::Gvcf,
                &mut stdout,
            )
            .map_err(|e| format!("inactive-reference-model: {e}"))?;
        }
        "gvcf-header" => {
            let contig = next_arg(&mut it, "contig");
            let length: u64 = next_arg(&mut it, "contig_length")
                .parse()
                .map_err(|_| "contig_length: invalid u64".to_string())?;
            let mut stdout = io::stdout().lock();
            dump_gvcf_header_tsv(contig.as_str(), length, &mut stdout)
                .map_err(|e| format!("gvcf-header: {e}"))?;
        }
        "gvcf-writer-blocks" => {
            let fixture = next_arg(&mut it, "loci_fixture.tsv");
            let mut stdout = io::stdout().lock();
            dump_gvcf_writer_from_loci_fixture_tsv(Path::new(&fixture), &mut stdout)
                .map_err(|e| format!("gvcf-writer-blocks: {e}"))?;
        }
        "ploidy-resolution" => {
            let sample_ploidy: u32 = next_arg(&mut it, "sample_ploidy")
                .parse()
                .map_err(|_| "sample_ploidy: invalid u32".to_string())?;
            let mut stdout = io::stdout().lock();
            dump_ploidy_resolution_tsv(sample_ploidy, &mut stdout)
                .map_err(|e| format!("ploidy-resolution: {e}"))?;
        }
        "annotate-core" => {
            let alt_count: usize = next_arg(&mut it, "alt_allele_count")
                .parse()
                .map_err(|_| "alt_allele_count: invalid usize".to_string())?;
            let samples_path = next_arg(&mut it, "samples.tsv");
            let mut stdout = io::stdout().lock();
            dump_annotate_core_tsv(alt_count, Path::new(&samples_path), &mut stdout)
                .map_err(|e| format!("annotate-core: {e}"))?;
        }
        "annotation-manifest" => {
            let mut stdout = io::stdout().lock();
            dump_annotation_manifest_tsv(&mut stdout)
                .map_err(|e| format!("annotation-manifest: {e}"))?;
        }
        "call-region-vcf" => {
            let ref_path = next_arg(&mut it, "ref.fa");
            let bam_path = next_arg(&mut it, "bam");
            let interval_cli = next_arg(&mut it, "interval_cli");
            let padding: u64 = it
                .next()
                .map(|s| s.parse().map_err(|_| "padding: invalid u64".to_string()))
                .transpose()?
                .unwrap_or(0);
            let bam = resolve_alignment_for_parity(Path::new(&bam_path))?;
            let mut stdout = io::stdout().lock();
            dump_call_region_vcf_tsv(
                Path::new(&ref_path),
                &bam,
                interval_cli.as_str(),
                padding,
                &mut stdout,
            )
            .map_err(|e| format!("call-region-vcf: {e}"))?;
        }
        "call-region-format" => {
            let ref_path = next_arg(&mut it, "ref.fa");
            let bam_path = next_arg(&mut it, "bam");
            let interval_cli = next_arg(&mut it, "interval_cli");
            let padding: u64 = it
                .next()
                .map(|s| s.parse().map_err(|_| "padding: invalid u64".to_string()))
                .transpose()?
                .unwrap_or(0);
            let bam = resolve_alignment_for_parity(Path::new(&bam_path))?;
            let mut stdout = io::stdout().lock();
            dump_call_region_format_tsv(
                Path::new(&ref_path),
                &bam,
                interval_cli.as_str(),
                padding,
                &mut stdout,
            )
            .map_err(|e| format!("call-region-format: {e}"))?;
        }
        "variant-vcf-from-gl-ad" => {
            let contig = next_arg(&mut it, "contig");
            let pos: u64 = next_arg(&mut it, "pos")
                .parse()
                .map_err(|_| "pos: invalid u64".to_string())?;
            let ref_allele = next_arg(&mut it, "ref_allele");
            let alt_allele = next_arg(&mut it, "alt_allele");
            let gl_csv = next_arg(&mut it, "gl_csv");
            let ad_csv = next_arg(&mut it, "ad_csv");
            let mut stdout = io::stdout().lock();
            dump_variant_vcf_from_gl_ad_tsv(
                contig.as_str(),
                pos,
                ref_allele.as_str(),
                alt_allele.as_str(),
                gl_csv.as_str(),
                ad_csv.as_str(),
                &mut stdout,
            )
            .map_err(|e| format!("variant-vcf-from-gl-ad: {e}"))?;
        }
        "variant-format-from-gl-ad" => {
            let contig = next_arg(&mut it, "contig");
            let pos: u64 = next_arg(&mut it, "pos")
                .parse()
                .map_err(|_| "pos: invalid u64".to_string())?;
            let ref_allele = next_arg(&mut it, "ref_allele");
            let alt_allele = next_arg(&mut it, "alt_allele");
            let gl_csv = next_arg(&mut it, "gl_csv");
            let ad_csv = next_arg(&mut it, "ad_csv");
            let mut stdout = io::stdout().lock();
            dump_variant_format_from_gl_ad_tsv(
                contig.as_str(),
                pos,
                ref_allele.as_str(),
                alt_allele.as_str(),
                gl_csv.as_str(),
                ad_csv.as_str(),
                &mut stdout,
            )
            .map_err(|e| format!("variant-format-from-gl-ad: {e}"))?;
        }
        "af-em" => {
            let fixture = next_arg(&mut it, "fixture.tsv");
            let mut stdout = io::stdout().lock();
            dump_af_em_tsv(Path::new(&fixture), &mut stdout).map_err(|e| format!("af-em: {e}"))?;
        }
        "genotype-limits" => {
            let ploidy: u32 = next_arg(&mut it, "ploidy")
                .parse()
                .map_err(|_| "ploidy: invalid u32".to_string())?;
            let max_gt: u32 = next_arg(&mut it, "max_genotype_count")
                .parse()
                .map_err(|_| "max_genotype_count: invalid u32".to_string())?;
            let mut stdout = io::stdout().lock();
            dump_genotype_limits_tsv(ploidy, max_gt, &mut stdout)
                .map_err(|e| format!("genotype-limits: {e}"))?;
        }
        "genotype-phasing" => {
            let alleles = next_arg(&mut it, "alleles_csv");
            let enabled = next_arg(&mut it, "phasing_enabled") == "1";
            let phase_set = it
                .next()
                .and_then(|s| if s == "-" { None } else { s.parse().ok() });
            let mut stdout = io::stdout().lock();
            dump_genotype_phasing_tsv(alleles.as_str(), enabled, phase_set, &mut stdout)
                .map_err(|e| format!("genotype-phasing: {e}"))?;
        }
        "force-calling-genotype" => {
            let vcf = next_arg(&mut it, "alleles.vcf");
            let contig = next_arg(&mut it, "contig");
            let pos: u64 = next_arg(&mut it, "pos")
                .parse()
                .map_err(|_| "pos: invalid u64".to_string())?;
            let filtered = next_arg(&mut it, "force_call_filtered") == "1";
            let mut stdout = io::stdout().lock();
            dump_force_calling_genotype_tsv(
                Path::new(&vcf),
                contig.as_str(),
                pos,
                filtered,
                &mut stdout,
            )
            .map_err(|e| format!("force-calling-genotype: {e}"))?;
        }
        "allele-subsetting" => {
            let sums = next_arg(&mut it, "haplotype_log10_sums");
            let is_ref = next_arg(&mut it, "is_reference_flags");
            let max_alleles: usize = next_arg(&mut it, "max_allele_count")
                .parse()
                .map_err(|_| "max_allele_count: invalid usize".to_string())?;
            let mut stdout = io::stdout().lock();
            dump_allele_subsetting_tsv(sums.as_str(), is_ref.as_str(), max_alleles, &mut stdout)
                .map_err(|e| format!("allele-subsetting: {e}"))?;
        }
        "subset-alleles-pl" => {
            let fixture = next_arg(&mut it, "fixture.tsv");
            let mut stdout = io::stdout().lock();
            dump_subset_alleles_pl_tsv(Path::new(&fixture), &mut stdout)
                .map_err(|e| format!("subset-alleles-pl: {e}"))?;
        }
        "subset-alleles-vc" => {
            let fixture = next_arg(&mut it, "fixture.tsv");
            let mut stdout = io::stdout().lock();
            dump_subset_alleles_vc_tsv(Path::new(&fixture), &mut stdout)
                .map_err(|e| format!("subset-alleles-vc: {e}"))?;
        }
        "subset-alleles-integration" => {
            let sums = next_arg(&mut it, "haplotype_log10_sums");
            let is_ref = next_arg(&mut it, "is_reference_flags");
            let max_alleles: usize = next_arg(&mut it, "max_allele_count")
                .parse()
                .map_err(|_| "max_allele_count: invalid usize".to_string())?;
            let fixture = next_arg(&mut it, "vc_fixture.tsv");
            let mut stdout = io::stdout().lock();
            dump_subset_alleles_integration_tsv(
                sums.as_str(),
                is_ref.as_str(),
                max_alleles,
                Path::new(&fixture),
                &mut stdout,
            )
            .map_err(|e| format!("subset-alleles-integration: {e}"))?;
        }
        "standard-annotations" => {
            let ref_fw: u32 = next_arg(&mut it, "ref_fw")
                .parse()
                .map_err(|_| "ref_fw".to_string())?;
            let ref_rv: u32 = next_arg(&mut it, "ref_rv")
                .parse()
                .map_err(|_| "ref_rv".to_string())?;
            let alt_fw: u32 = next_arg(&mut it, "alt_fw")
                .parse()
                .map_err(|_| "alt_fw".to_string())?;
            let alt_rv: u32 = next_arg(&mut it, "alt_rv")
                .parse()
                .map_err(|_| "alt_rv".to_string())?;
            let qual: f64 = next_arg(&mut it, "qual")
                .parse()
                .map_err(|_| "qual".to_string())?;
            let dp: i32 = next_arg(&mut it, "dp")
                .parse()
                .map_err(|_| "dp".to_string())?;
            let ref_bqs = next_arg(&mut it, "ref_bqs");
            let alt_bqs = next_arg(&mut it, "alt_bqs");
            let ref_pos = it.next().unwrap_or_else(|| "-".to_string());
            let alt_pos = it.next().unwrap_or_else(|| "-".to_string());
            let ref_mq = it.next().unwrap_or_else(|| "-".to_string());
            let alt_mq = it.next().unwrap_or_else(|| "-".to_string());
            let mut stdout = io::stdout().lock();
            dump_standard_annotations_tsv(
                ref_fw,
                ref_rv,
                alt_fw,
                alt_rv,
                qual,
                dp,
                ref_bqs.as_str(),
                alt_bqs.as_str(),
                ref_pos.as_str(),
                alt_pos.as_str(),
                ref_mq.as_str(),
                alt_mq.as_str(),
                &mut stdout,
            )
            .map_err(|e| format!("standard-annotations: {e}"))?;
        }
        "emit-mode-decision" => {
            let mode = next_arg(&mut it, "mode");
            let has_variant = next_arg(&mut it, "has_variant") == "1";
            let locus_count: usize = next_arg(&mut it, "locus_count")
                .parse()
                .map_err(|_| "locus_count".to_string())?;
            let emit_mode = match mode.as_str() {
                "GVCF" => gatk_haplotypecaller::genotyping::EmitMode::Gvcf,
                "BP_RESOLUTION" => gatk_haplotypecaller::genotyping::EmitMode::BpResolution,
                _ => gatk_haplotypecaller::genotyping::EmitMode::Vcf,
            };
            let mut stdout = io::stdout().lock();
            dump_emit_mode_decision_tsv(emit_mode, has_variant, locus_count, &mut stdout)
                .map_err(|e| format!("emit-mode-decision: {e}"))?;
        }
        "as-annotations" => {
            let site_af: f64 = next_arg(&mut it, "site_af")
                .parse()
                .map_err(|_| "site_af".to_string())?;
            let site_qual: f64 = next_arg(&mut it, "site_qual")
                .parse()
                .map_err(|_| "site_qual".to_string())?;
            let mut stdout = io::stdout().lock();
            dump_as_annotations_tsv(site_af, site_qual, &mut stdout)
                .map_err(|e| format!("as-annotations: {e}"))?;
        }
        "excess-het" => {
            let ref_count: u32 = next_arg(&mut it, "ref_count")
                .parse()
                .map_err(|_| "ref_count".to_string())?;
            let het: u32 = next_arg(&mut it, "het_count")
                .parse()
                .map_err(|_| "het_count".to_string())?;
            let hom: u32 = next_arg(&mut it, "hom_count")
                .parse()
                .map_err(|_| "hom_count".to_string())?;
            let mut stdout = io::stdout().lock();
            dump_excess_het_tsv(ref_count, het, hom, &mut stdout)
                .map_err(|e| format!("excess-het: {e}"))?;
        }
        "depth-per-sample-hc" => {
            let ad_csv = next_arg(&mut it, "ad_csv");
            let mut stdout = io::stdout().lock();
            dump_depth_per_sample_hc_tsv(ad_csv.as_str(), &mut stdout)
                .map_err(|e| format!("depth-per-sample-hc: {e}"))?;
        }
        "annotation-plugin" => {
            let plugin = next_arg(&mut it, "plugin");
            let ref_fw: u32 = next_arg(&mut it, "ref_fw")
                .parse()
                .map_err(|_| "ref_fw".to_string())?;
            let ref_rv: u32 = next_arg(&mut it, "ref_rv")
                .parse()
                .map_err(|_| "ref_rv".to_string())?;
            let alt_fw: u32 = next_arg(&mut it, "alt_fw")
                .parse()
                .map_err(|_| "alt_fw".to_string())?;
            let alt_rv: u32 = next_arg(&mut it, "alt_rv")
                .parse()
                .map_err(|_| "alt_rv".to_string())?;
            let qual: f64 = next_arg(&mut it, "qual")
                .parse()
                .map_err(|_| "qual".to_string())?;
            let dp: i32 = next_arg(&mut it, "dp")
                .parse()
                .map_err(|_| "dp".to_string())?;
            let ref_bqs = next_arg(&mut it, "ref_bqs");
            let alt_bqs = next_arg(&mut it, "alt_bqs");
            let mut stdout = io::stdout().lock();
            dump_annotation_plugin_tsv(
                plugin.as_str(),
                ref_fw,
                ref_rv,
                alt_fw,
                alt_rv,
                qual,
                dp,
                ref_bqs.as_str(),
                alt_bqs.as_str(),
                &mut stdout,
            )
            .map_err(|e| format!("annotation-plugin: {e}"))?;
        }
        "bamout-stub" => {
            let enabled = next_arg(&mut it, "enabled") == "1";
            let count: usize = next_arg(&mut it, "write_count")
                .parse()
                .map_err(|_| "write_count".to_string())?;
            let mut stdout = io::stdout().lock();
            dump_bamout_stub_tsv(enabled, count, &mut stdout)
                .map_err(|e| format!("bamout-stub: {e}"))?;
        }
        "dragen-mode-branch" => {
            let mut stdout = io::stdout().lock();
            dump_dragen_mode_branch_tsv(&mut stdout)
                .map_err(|e| format!("dragen-mode-branch: {e}"))?;
        }
        "gvcf-l5-merged" => {
            let fixture = next_arg(&mut it, "fixture.tsv");
            let mut stdout = io::stdout().lock();
            dump_gvcf_l5_merged_tsv(Path::new(&fixture), &mut stdout)
                .map_err(|e| format!("gvcf-l5-merged: {e}"))?;
        }
        "dragstr-calibration" => {
            let loaded = next_arg(&mut it, "params_loaded") == "1";
            let mut stdout = io::stdout().lock();
            dump_dragstr_calibration_tsv(loaded, &mut stdout)
                .map_err(|e| format!("dragstr-calibration: {e}"))?;
        }
        "assembly-debug-stub" => {
            let failure_bam = next_arg(&mut it, "assembly_failure_bam") == "1";
            let graph_dot = next_arg(&mut it, "graph_dot") == "1";
            let mut stdout = io::stdout().lock();
            dump_assembly_debug_stub_tsv(failure_bam, graph_dot, &mut stdout)
                .map_err(|e| format!("assembly-debug-stub: {e}"))?;
        }
        _ => usage_and_exit(),
    }
    Ok(())
}
