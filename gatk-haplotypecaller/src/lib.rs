//! GATK-RS HaplotypeCaller Library
//! This crate contains the HaplotypeCaller implementation for GATK-RS.
//! # Public API surface
//! Default features expose the **HC engine / product** API used by `gatk-cli`
//! (run/engine/region/assembly/genotyping/PairHMM/VCF emit). Internal parity TSV
//! exporters require the `dev-dumps` feature (off by default; enabled via
//! `gatk-cli` and `parity_harness`). P12 fixture tables and discovery/compat
//! helpers are crate-private unless `parity_harness` is enabled (integration
//! tests + the `hc_full_parity_gate` example).
#![allow(clippy::result_large_err)]
// Without `parity_harness`, P12/compat/dump scaffolding stays compiled but unwired
// into the product path — suppress that dead-code noise. CI uses `--all-features`.
#![cfg_attr(not(feature = "parity_harness"), allow(dead_code, unused_imports))]

// --------------------------------------------------------------------------
// Product modules (stable / embedding surface)
// --------------------------------------------------------------------------
pub mod active_region;
pub mod activity_profile;
pub mod activity_scoring;
pub mod alignment;
pub mod allele_filter_options;
pub mod allele_filtering;
pub mod assembly;
pub mod assembly_based_caller;
pub mod assembly_dangling_recovery;
pub mod assembly_pipeline_stages;
pub mod assembly_pruning;
pub mod assembly_region_evaluator;
pub mod assembly_region_finalize;
pub mod assembly_region_iterator;
pub mod assembly_region_trimmer;
pub mod assembly_result_set;
pub mod bio_ids;
pub mod cigar;
pub mod cigar_builder;
/// Multi-sample gVCF merge (GATK CombineGVCFs algorithm slice).
pub mod combine_gvcfs;
pub mod emit_gates;
pub mod engine;
pub mod event_map;
pub mod feature_context;
pub mod genome_loc;
/// Joint genotyping of combined gVCFs (GATK GenotypeGVCFs algorithm slice).
pub mod genotype_gvcfs;
pub mod genotype_site;
pub mod genotyping;
pub mod given_alleles;
pub mod gvcf_writer;
pub mod haplotype;
pub mod haplotype_cigar;
pub mod hc_allele_mapping;
pub mod hc_emit_policy;
pub mod hc_genotyping_engine;
pub mod hq_soft_clip;
pub mod junction_kbest;
pub mod junction_tree_graph;
pub mod kbest_haplotype;
pub mod likelihood_engine;
pub mod locus_iterator;
pub mod minimal_genotyping;
pub mod multiallelic_emit;
pub mod pairhmm;
pub mod pairhmm_log10;
pub mod pairhmm_logless;
pub mod pairhmm_qual;
pub mod pairhmm_simd;
pub mod pcr_error_model;
pub mod pileup_detection;
pub mod pileup_element;
pub mod read_assembly_filter;
pub mod read_downsample;
pub mod read_header_semantics;
pub mod read_model;
pub mod read_projection;
pub mod read_realignment;
pub mod read_threading_assembler;
pub mod read_threading_graph;
pub mod read_transformer;
pub mod read_validation;
pub mod ref_confidence;
pub mod ref_confidence_merger;
pub mod reference_context;
pub mod reference_vcf_emit;
/// Leaf likelihood-row type (breaks engine ↔ genotyping module cycle).
pub mod region_read_likelihood;
pub mod region_vcf_emit;
pub mod run;
/// Process env → typed config (Sprint I-4). Prefer over scattered `std::env::var`.
pub mod runtime_config;
/// Optional observe-only semantic checkpoints (`GATK_RS_SEMANTIC_TRACE`).
pub mod semantic_trace;
pub mod seq_graph;
pub mod seq_kbest_haplotype;
/// Arc-backed BAM records shared between the shard read cache and assembly regions.
pub mod shared_bam;
pub mod smith_waterman;
pub mod variant_site_hc_annotations;
pub mod walker;
pub mod walker_apply;
pub mod walker_traversal;
pub mod worker_threads;

// --------------------------------------------------------------------------
// Implementation-detail modules (never part of the supported embedding API)
// --------------------------------------------------------------------------
pub(crate) mod af_calc;
pub(crate) mod allele_downsample;
pub(crate) mod allele_subsetting;
pub(crate) mod allele_subsetting_pl;
pub(crate) mod annotator;
pub(crate) mod bamout;
pub(crate) mod dragen_mq;
pub(crate) mod event_map_rebuild;
pub mod fragment_overlap;
pub(crate) mod gatk_well_rng;
pub(crate) mod genotype_limits;
pub(crate) mod hc_joint_is_active;
pub(crate) mod long_homopolymer_collapsing;
pub(crate) mod mann_whitney_u;
pub(crate) mod nearby_kmer_error_corrector;
pub(crate) mod read_binding;
pub(crate) mod read_error_correction;
pub(crate) mod read_optional_tags;
pub(crate) mod read_pre_len;
pub(crate) mod read_pre_mate;
pub(crate) mod read_pre_mq;
pub mod read_unclip;
pub(crate) mod region_pileup;
pub(crate) mod seq_graph_simplify;

// --------------------------------------------------------------------------
// Parity / dump / P12 quarantine — public only with `parity_harness`
// --------------------------------------------------------------------------
#[cfg(feature = "dev-dumps")]
pub mod annotator_dump;
#[cfg(all(test, not(feature = "dev-dumps")))]
pub(crate) mod annotator_dump;
#[cfg(feature = "dev-dumps")]
pub mod assembler_args_dump;
#[cfg(all(test, not(feature = "dev-dumps")))]
pub(crate) mod assembler_args_dump;
#[cfg(feature = "dev-dumps")]
pub mod assembly_debug_dump;
#[cfg(all(test, not(feature = "dev-dumps")))]
pub(crate) mod assembly_debug_dump;
#[cfg(feature = "dev-dumps")]
pub mod assembly_graph_dump;
#[cfg(all(test, not(feature = "dev-dumps")))]
pub(crate) mod assembly_graph_dump;
#[cfg(feature = "dev-dumps")]
pub mod assembly_region_assemble_dump;
#[cfg(all(test, not(feature = "dev-dumps")))]
pub(crate) mod assembly_region_assemble_dump;
#[cfg(feature = "dev-dumps")]
pub mod assembly_region_finalize_dump;
#[cfg(all(test, not(feature = "dev-dumps")))]
pub(crate) mod assembly_region_finalize_dump;
#[cfg(feature = "dev-dumps")]
pub mod assembly_region_genotype_dump;
#[cfg(all(test, not(feature = "dev-dumps")))]
pub(crate) mod assembly_region_genotype_dump;
#[cfg(feature = "dev-dumps")]
pub mod assembly_region_pairhmm_dump;
#[cfg(all(test, not(feature = "dev-dumps")))]
pub(crate) mod assembly_region_pairhmm_dump;
#[cfg(feature = "dev-dumps")]
pub mod assembly_region_stages_dump;
#[cfg(all(test, not(feature = "dev-dumps")))]
pub(crate) mod assembly_region_stages_dump;
#[cfg(feature = "dev-dumps")]
pub mod assembly_regions_dump;
#[cfg(all(test, not(feature = "dev-dumps")))]
pub(crate) mod assembly_regions_dump;
/// Java-compatibility / P12-interval semantics (Sprint I-1 quarantine).
#[cfg(feature = "parity_harness")]
pub mod compatibility;
#[cfg(not(feature = "parity_harness"))]
pub(crate) mod compatibility;
#[cfg(feature = "dev-dumps")]
pub mod genotyping_dump;
#[cfg(all(test, not(feature = "dev-dumps")))]
pub(crate) mod genotyping_dump;
#[cfg(feature = "dev-dumps")]
pub mod gvcf_writer_dump;
#[cfg(all(test, not(feature = "dev-dumps")))]
pub(crate) mod gvcf_writer_dump;
#[cfg(feature = "dev-dumps")]
pub mod hq_soft_clip_dump;
#[cfg(all(test, not(feature = "dev-dumps")))]
pub(crate) mod hq_soft_clip_dump;
#[cfg(feature = "dev-dumps")]
pub mod j_modes_dump;
#[cfg(all(test, not(feature = "dev-dumps")))]
pub(crate) mod j_modes_dump;
#[cfg(feature = "dev-dumps")]
pub mod likelihood_engine_dump;
#[cfg(all(test, not(feature = "dev-dumps")))]
pub(crate) mod likelihood_engine_dump;
#[cfg(feature = "dev-dumps")]
pub mod locus_pileup_detail_dump;
#[cfg(all(test, not(feature = "dev-dumps")))]
pub(crate) mod locus_pileup_detail_dump;
#[cfg(feature = "parity_harness")]
pub mod p12_java_format_fixup;
#[cfg(not(feature = "parity_harness"))]
pub(crate) mod p12_java_format_fixup;
#[cfg(feature = "dev-dumps")]
pub mod pairhmm_dump;
#[cfg(all(test, not(feature = "dev-dumps")))]
pub(crate) mod pairhmm_dump;
#[cfg(feature = "dev-dumps")]
pub mod pairhmm_f3_dump;
#[cfg(all(test, not(feature = "dev-dumps")))]
pub(crate) mod pairhmm_f3_dump;
#[cfg(feature = "dev-dumps")]
pub mod pairhmm_native_dump;
#[cfg(all(test, not(feature = "dev-dumps")))]
pub(crate) mod pairhmm_native_dump;
#[cfg(feature = "parity_harness")]
pub mod parity_harness;
#[cfg(not(feature = "parity_harness"))]
pub(crate) mod parity_harness;
#[cfg(feature = "parity_harness")]
pub mod parity_region_genotype;
#[cfg(not(feature = "parity_harness"))]
pub(crate) mod parity_region_genotype;
#[cfg(feature = "dev-dumps")]
pub mod ploidy_dump;
#[cfg(all(test, not(feature = "dev-dumps")))]
pub(crate) mod ploidy_dump;
#[cfg(feature = "dev-dumps")]
pub mod pre_dragstr_dump;
#[cfg(all(test, not(feature = "dev-dumps")))]
pub(crate) mod pre_dragstr_dump;
/// Read-event fallback / P12 discovery — not part of Java `HaplotypeCaller.callRegion` API.
#[doc(hidden)]
#[cfg(feature = "parity_harness")]
pub mod read_event_discovery;
#[doc(hidden)]
#[cfg(not(feature = "parity_harness"))]
pub(crate) mod read_event_discovery;
#[cfg(feature = "dev-dumps")]
pub mod read_filter_dump;
#[cfg(all(test, not(feature = "dev-dumps")))]
pub(crate) mod read_filter_dump;
#[cfg(feature = "dev-dumps")]
pub mod read_pre_len_dump;
#[cfg(all(test, not(feature = "dev-dumps")))]
pub(crate) mod read_pre_len_dump;
#[cfg(feature = "dev-dumps")]
pub mod read_pre_mq_dump;
#[cfg(all(test, not(feature = "dev-dumps")))]
pub(crate) mod read_pre_mq_dump;
#[cfg(feature = "dev-dumps")]
pub mod read_pre_overlap_dump;
#[cfg(all(test, not(feature = "dev-dumps")))]
pub(crate) mod read_pre_overlap_dump;
#[cfg(feature = "dev-dumps")]
pub mod read_transformer_dump;
#[cfg(all(test, not(feature = "dev-dumps")))]
pub(crate) mod read_transformer_dump;
#[cfg(feature = "dev-dumps")]
pub mod read_unclip_dump;
#[cfg(all(test, not(feature = "dev-dumps")))]
pub(crate) mod read_unclip_dump;
#[cfg(feature = "dev-dumps")]
pub mod ref_confidence_dump;
#[cfg(all(test, not(feature = "dev-dumps")))]
pub(crate) mod ref_confidence_dump;
#[cfg(feature = "dev-dumps")]
pub mod ref_confidence_merger_dump;
#[cfg(all(test, not(feature = "dev-dumps")))]
pub(crate) mod ref_confidence_merger_dump;
#[cfg(feature = "dev-dumps")]
pub mod region_vcf_dump;
#[cfg(all(test, not(feature = "dev-dumps")))]
pub(crate) mod region_vcf_dump;

/// Crate-internal path for P12 tables (call sites use `crate::java_hc_site_semantics::…`).
/// Not part of the default public API; with `parity_harness` the same path is `pub`.
#[cfg(feature = "parity_harness")]
pub use compatibility::java_hc_site_semantics;
#[cfg(not(feature = "parity_harness"))]
pub(crate) use compatibility::java_hc_site_semantics;

// --------------------------------------------------------------------------
// Product re-exports (default feature set)
// --------------------------------------------------------------------------
pub use active_region::{TraversalTile, DEFAULT_TRAVERSAL_TILE_BP};
pub use activity_profile::{
    activity_profile_base_process_state, adaptive_filter_size, band_pass_process_state,
    make_gaussian_kernel, normal_distribution, ActivityEvidence, ActivityProfileRegion,
    ActivityProfileState, ActivityProfileStateKind, BandPassActivityProfile,
    BandPassActivityProfileParams, PositiveSigma, GATK_BAND_PASS_DEFAULT_SIGMA,
    GATK_BAND_PASS_MAX_FILTER_SIZE, GATK_BAND_PASS_MIN_PROB_TO_KEEP_IN_FILTER,
    GATK_DEFAULT_ACTIVE_PROB_THRESHOLD, GATK_DEFAULT_MAX_PROB_PROPAGATION_DISTANCE,
};
pub use activity_scoring::{
    approximate_log10_sum_log10_pair, calc_ref_vs_any_log10_genotype_likelihoods,
    calculate_multi_sample_any_non_ref_posterior,
    calculate_single_sample_biallelic_non_ref_posterior,
    haplotype_caller_activity_profile_state_multi_sample,
    haplotype_caller_activity_profile_state_single_sample, is_alt_after_assembly,
    is_alt_before_assembly, log10_sum_log10, log_binomial_coefficient_natural,
    normalize_from_log10_to_linear_space, HaplotypeCallerActivityScoringParams, PileupObservation,
    AVERAGE_HQ_SOFTCLIPS_HQ_BASES_THRESHOLD, LOG10_ONE_THIRD, REF_MODEL_DELETION_QUAL,
};
pub use alignment::*;
pub use alignment::{calculate_haplotype_cigar, Cigar, CigarOperator, SwParameters};
pub use allele_downsample::{
    apply_contamination_to_pileup, select_allele_biased_evidence_indices, target_allele_counts,
};
pub use allele_filter_options::{ActiveRegionSpan, AlleleFilterOptions};
pub use allele_filtering::{
    filter_assembly_and_likelihoods, MAX_NON_REF_HAPLOTYPES_FOR_GENOTYPING,
};
pub use annotator::{
    annotate_parity_v1_site, AnnotatedSite, VariantAnnotationContext, PARITY_V1_FORMAT_KEYS,
    PARITY_V1_INFO_KEYS,
};
pub use assembly::*;
pub use assembly::{AssemblyGraphPruningParams, AssemblyGraphSummary};
pub use assembly_based_caller::{
    assemble_reads, assemble_reads_with_finalized, call_region_assemble, AssembleReadsArgs,
    AssembledRegion,
};
pub use assembly_pipeline_stages::{
    CallRegionAssemblyStage, EVENT_MAP_SYNC_AROUND_FILTER_RATIONALE,
};
pub use assembly_region_evaluator::{add_locus_for_smoothed_activity, evaluate_hc_activity_state};
pub use assembly_region_iterator::{
    load_all_records_for_contig, load_all_records_for_contig_raw, load_records_for_shard_raw,
    refuse_oversized_assembly_region_reads, sync_read_qnames, AssemblyRegion,
    AssemblyRegionIterator, AssemblyRegionIteratorConfig, MAX_READS_PER_ASSEMBLY_REGION,
};
pub use assembly_region_trimmer::{
    load_trim_variants_tsv, trim_assembly_region, AssemblyRegionTrimResult, AssemblyRegionTrimmer,
    AssemblyRegionTrimmerConfig, TrimVariant,
};
/// CLI diagnostic (`gatk-cli DumpSmoothedActivity`) — requires `dev-dumps`.
#[cfg(feature = "dev-dumps")]
pub use assembly_regions_dump::dump_smoothed_activity_tsv;
pub use assembly_result_set::AssemblyResultSet;
pub use bamout::{BamoutWriter, BamoutWriterConfig};
pub use bio_ids::{
    AlleleDepth, AlleleIndex, DiploidGenotypeIndex, GenotypeQuality, HaplotypeIndex, KmerSize,
    MappingQuality, PadOffset0, PhredLikelihood, Ploidy, ReadCoordinate, ReadDepth, ReadIndex,
    ReferenceCoordinate, SampleIndex,
};
pub use combine_gvcfs::{run_combine_gvcfs, CombineGvcfsArgs};
pub use engine::{CallRegionArgs, CallRegionMode, CallRegionOutcome, HaplotypeCallerEngine};
pub use event_map::{Event, EventMap, IndelSpan};
pub use feature_context::{FeatureContext, FeatureDataSources, FeatureLocatable};
pub use gatk_well_rng::{Well19937c, GATK_WELL19937C_SEED};
pub use genome_loc::{GenomeLoc, GenomePosition};
pub use genotype_gvcfs::{run_genotype_gvcfs, GenotypeGvcfsArgs, DEFAULT_STAND_CALL_CONF};
pub use genotyping::*;
pub use given_alleles::{merge_given_alleles_into_assembly, GatkGivenAllele};
pub use gvcf_writer::{
    gatk_hc_gvcf_header_lines, GvcfWriter, GvcfWriterConfig, GATK_HC_DEFAULT_GQB,
};
pub use haplotype::Haplotype;
pub use hc_genotyping_engine::{
    biallelic_genotype_log10_likelihoods_gatk, diagnose_genotype_variation_event,
    genotype_active_region, java_emit_af_decision, java_vcf_shaped_rescue_gl,
    marginalize_rows_to_biallelic_alleles, subset_biallelic_haplotype_indices,
    GenotypeRejectReason, GenotypingSemantics, HcGenotypingConfig, InformativeAd,
    JavaEmitAfDecision, RegionGenotypeResult, SparsePlShape,
    DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN, DEFAULT_STAND_EMIT_CONFIDENCE,
};
#[cfg(feature = "dev-dumps")]
pub use hc_genotyping_engine::{format_locus_genotype_pl_dump, pairhmm_locus_trace_dump};
pub use hq_soft_clip::{
    count_high_quality_soft_clip_bases_rcm, hq_soft_clip_running_mean_at_locus,
    max_hq_soft_clip_bases, RCM_HQ_SOFT_CLIP_QUAL_THRESHOLD,
};
pub use junction_kbest::{find_junction_best_haplotypes, JunctionKBestPath};
pub use junction_tree_graph::build_junction_tree_graph_from_ref_and_reads;
pub use kbest_haplotype::{
    find_best_haplotypes, find_best_haplotypes_for_assembly,
    find_best_haplotypes_preserving_cycles, KBestPath,
};
pub use likelihood_engine::{
    prepare_read_quals_for_pairhmm_inplace, score_read_against_haplotypes,
    HcLikelihoodEngineConfig, HcLikelihoodImplementation,
};
pub use locus_iterator::{IntervalLocusIterator, LocusPileupState, LocusPileupWalker};
pub use minimal_genotyping::{
    calculate_single_sample_ref_vs_any_active_state_profile_value,
    haplotype_caller_activity_profile_state_minimal_genotyping,
};
pub use pairhmm::{
    pairhmm_fp_eq, pairhmm_log10_likelihood, pairhmm_log10_likelihood_slices,
    pairhmm_log10_likelihoods_vectorized, pairhmm_log10_likelihoods_vectorized_slices,
    PairHmmFpPolicy, PairHmmInput, PairHmmParams,
};
pub use pairhmm_log10::{
    log10_pairhmm_likelihood as log10_pairhmm_likelihood_exact,
    log10_pairhmm_likelihood_parity_defaults, GATK_PARITY_DEFAULT_DEL_QUAL,
    GATK_PARITY_DEFAULT_GCP, GATK_PARITY_DEFAULT_INS_QUAL,
};
pub use pairhmm_logless::{
    logless_pairhmm_likelihood, logless_pairhmm_likelihood_parity_defaults, INITIAL_CONDITION,
    INITIAL_CONDITION_LOG10,
};
pub use pairhmm_simd::{
    best_simd_backend, parse_pair_hmm_impl, resolve_pair_hmm_impl, score_read_haps_logless,
    PairHmmBackend, PairHmmImpl,
};
pub use read_downsample::{
    apply_positional_downsampler, GatkJavaRng, PositionalDownsamplerConfig,
    GATK_DEFAULT_MAX_READS_PER_ALIGNMENT_START, GATK_JAVA_RANDOM_SEED,
};
/// Production ASM-8 / P12-bridge env flags (public so integration tests need no `parity_harness`).
pub use read_event_discovery::{
    strict_java_asm8_only_enabled, strict_java_p12_ensure_bridges_enabled,
};
pub use read_header_semantics::{ReadHeaderSemantics, ResolvedReadHeaderSemantics};
pub use read_model::{
    mapq_passes_minimum, passes_hc_read_filter_set, passes_hc_read_filters,
    passes_hc_read_filters_fields, passes_hc_read_filters_with_header,
    passes_printreads_parity_filters, standard_hc_read_filter_failure_index, HcReadFilterSet,
    ReadFilterParams, FLAG_DUPLICATE, FLAG_NOT_PRIMARY, FLAG_SEGMENT_UNMAPPED, FLAG_SUPPLEMENTARY,
    FLAG_VENDOR_QUALITY_FAILED, GATK_HC_DEFAULT_MIN_MAPPING_QUALITY, MAPPING_QUALITY_UNAVAILABLE,
    STANDARD_HC_READ_FILTER_JAVA_NAMES,
};
pub use read_projection::{
    cigar_hard_clip_length, cigar_soft_clip_ends, query_index_at_reference_position,
    reference_position_at_query_index,
};
pub use read_threading_assembler::{
    assemble_from_ref_and_reads, audit_threading_dangling_recovery,
    AssemblyResult as ThreadingAssemblyResult, AssemblyStatus, ReadThreadingAssemblerArgs,
    ThreadingDanglingAudit, DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH,
};
pub use read_threading_assembler::{
    audit_kbest_extract, path_is_too_divergent_from_reference, KbestExtractAuditRow,
    KbestExtractReject,
};
pub use read_threading_assembler::{probe_seq_graph_kmer_attempts, SeqGraphKmerProbeRow};
pub use read_threading_graph::{
    reference_has_non_unique_kmers, threading_non_unique_summary, ThreadingNonUniqueSummary,
};
pub use read_transformer::{
    apply_iupac_strict_transform, apply_shard_read_pipeline, load_contig_records_hc_production,
    ShardReadPipelineConfig,
};
pub use read_validation::validate_mapped_read_sanity;
pub use ref_confidence::{
    calc_genotype_likelihoods_of_ref_vs_any, reference_gq_from_log10_gl,
    reference_model_for_no_variation_region, InactiveReferenceModelOutcome,
    ReferenceConfidenceConfig, ReferenceConfidenceLocusDetail,
};
pub use reference_context::ReferenceContext;
pub use region_pileup::RegionPileupLocus;
pub use region_read_likelihood::RegionReadLikelihood;
pub use region_vcf_emit::{
    try_emit_call_region_variant, try_emit_call_region_variants, HC_PIPELINE_ASSEMBLY_REGION_V1,
    HC_PIPELINE_LEGACY_PROVISIONAL, HC_PIPELINE_SCAFFOLD,
};
pub use run::run_haplotype_caller;
pub use shared_bam::{
    empty_shared_record, empty_shared_record_ref, into_unique_records, is_empty_shared_record,
    record_make_mut, share_record, share_records, BamRecordSlot, SharedBamRecord,
};
pub use walker::{
    make_read_shards, make_read_shards_default_padding, ReadShard,
    GATK_DEFAULT_ASSEMBLY_REGION_PADDING,
};
pub use walker_apply::{call_disposition, AssemblyRegionCallDisposition, WalkerApplyStats};
pub use walker_traversal::{
    collect_assembly_regions, drain_assembly_regions_for_shard, flatten_assembly_regions,
    for_each_assembly_region, into_assembly_regions, traverse_assembly_region_walker,
    WalkerShardTraversal, WalkerTraversal, WalkerTraversalConfig,
};
pub use worker_threads::init_worker_threads;

// --------------------------------------------------------------------------
// Dev dump re-exports (`dev-dumps`; not part of the default product API)
// --------------------------------------------------------------------------
#[cfg(feature = "dev-dumps")]
pub use allele_downsample::{dump_allele_biased_evidence_locus_tsv, dump_target_allele_counts_tsv};
#[cfg(feature = "dev-dumps")]
pub use annotator_dump::{
    dump_annotate_core_tsv, dump_annotation_manifest_tsv, dump_annotation_plugin_tsv,
    dump_as_annotations_tsv, dump_depth_per_sample_hc_tsv, dump_excess_het_tsv,
    dump_standard_annotations_tsv,
};
#[cfg(feature = "dev-dumps")]
pub use assembler_args_dump::dump_assembler_args_tsv;
#[cfg(feature = "dev-dumps")]
pub use assembly_debug_dump::dump_assembly_debug_stub_tsv;
#[cfg(feature = "dev-dumps")]
pub use assembly_graph_dump::{
    dump_assembly_assemble_tsv, dump_assembly_graph_dangling_summary_tsv,
    dump_assembly_graph_edges_tsv, dump_assembly_graph_low_quality_tsv,
    dump_assembly_graph_multi_kmer_edges_tsv, dump_assembly_graph_non_unique_summary_tsv,
    dump_assembly_graph_pruned_summary_tsv, dump_assembly_haplotype_cigars_tsv,
    dump_assembly_haplotypes_cap_tsv, dump_assembly_haplotypes_production_tsv,
    dump_assembly_haplotypes_tsv, dump_assembly_junction_haplotypes_tsv,
    dump_assembly_kbest_paths_tsv, dump_assembly_seqgraph_summary_tsv,
    dump_read_error_correction_tsv, load_assembly_reads_tsv, load_assembly_ref_tsv,
    write_assembly_graph_summary_tsv, write_dangling_recovery_summary_tsv,
    write_threading_non_unique_summary_tsv,
};
#[cfg(feature = "dev-dumps")]
pub use assembly_region_assemble_dump::{
    dump_assembly_region_haplotypes_tsv, dump_assembly_region_kmer_probe_tsv,
    AssemblyRegionHaplotypeTarget,
};
#[cfg(feature = "dev-dumps")]
pub use assembly_region_finalize_dump::{
    dump_assembly_region_assembly_stages_finalize_tsv, dump_assembly_region_finalize_reads_tsv,
};
#[cfg(feature = "dev-dumps")]
pub use assembly_region_genotype_dump::{
    dump_assembly_region_genotype_subset_tsv, dump_assembly_region_genotype_tsv,
};
#[cfg(feature = "dev-dumps")]
pub use assembly_region_iterator::{
    dump_assembly_region_features_tsv, dump_assembly_region_pileup_track_tsv,
    dump_assembly_region_reads_tsv, dump_assembly_region_reference_tsv,
};
#[cfg(feature = "dev-dumps")]
pub use assembly_region_pairhmm_dump::dump_assembly_region_pairhmm_likelihoods_tsv;
#[cfg(feature = "dev-dumps")]
pub use assembly_region_stages_dump::{
    dump_assembly_region_assembly_stages_tsv, dump_assembly_region_kbest_paths_tsv,
};
#[cfg(feature = "dev-dumps")]
pub use assembly_region_trimmer::dump_assembly_region_trim_tsv;
#[cfg(feature = "dev-dumps")]
pub use assembly_regions_dump::{
    dump_active_locus_tsv, dump_genotype_likelihood_activity_tsv, dump_locus_pileup_tsv,
    dump_locus_pileup_tsv_default_padding, dump_raw_activity_profile_tsv,
    dump_raw_activity_profile_tsv_with_contamination,
    dump_raw_activity_profile_tsv_with_force_calling, dump_smoothed_activity_profile_tsv,
    format_activity_prob, ForceCallingAllelesDump, GATK_HC_ALLELES_FEATURE_SOURCE,
};
#[cfg(feature = "dev-dumps")]
pub use bamout::dump_bamout_stub_tsv;
#[cfg(feature = "dev-dumps")]
pub use genotyping_dump::{
    dump_af_em_tsv, dump_allele_subsetting_tsv, dump_force_calling_genotype_tsv,
    dump_genotype_format_tsv, dump_genotype_limits_tsv, dump_genotype_phasing_tsv,
    dump_genotyping_aggregate_tsv, dump_subset_alleles_integration_tsv, dump_subset_alleles_pl_tsv,
    dump_subset_alleles_vc_tsv,
};
#[cfg(feature = "dev-dumps")]
pub use gvcf_writer_dump::{
    dump_gvcf_header_tsv, dump_gvcf_l5_merged_tsv, dump_gvcf_writer_from_loci_fixture_tsv,
};
#[cfg(feature = "dev-dumps")]
pub use hq_soft_clip_dump::dump_hq_soft_clip_mean_tsv;
#[cfg(feature = "dev-dumps")]
pub use j_modes_dump::{dump_dragen_mode_branch_tsv, dump_emit_mode_decision_tsv};
#[cfg(feature = "dev-dumps")]
pub use likelihood_engine_dump::{
    dump_likelihood_engine_config_tsv, dump_likelihood_pcr_read_tsv, dump_pcr_error_model_tsv,
};
#[cfg(feature = "dev-dumps")]
pub use locus_pileup_detail_dump::dump_locus_pileup_detail_tsv;
#[cfg(feature = "dev-dumps")]
pub use pairhmm_dump::{dump_pairhmm_likelihoods_tsv, load_pairhmm_cases_tsv};
#[cfg(feature = "dev-dumps")]
pub use pairhmm_f3_dump::{dump_pairhmm_bq_cap_tsv, dump_pairhmm_haplotype_filter_tsv};
#[cfg(feature = "dev-dumps")]
pub use pairhmm_native_dump::dump_pairhmm_native_likelihoods_tsv;
#[cfg(feature = "dev-dumps")]
pub use ploidy_dump::dump_ploidy_resolution_tsv;
#[cfg(feature = "dev-dumps")]
pub use pre_dragstr_dump::dump_dragstr_calibration_tsv;
/// Deterministic ingress order helpers (foundation Step 35 + HC tiling).
pub use read_binding::{
    count_reads_overlapping_tile, filtered_read_iteration_order, total_read_tile_overlaps,
};
#[cfg(feature = "dev-dumps")]
pub use read_downsample::dump_positional_downsample_summary_tsv;
#[cfg(feature = "dev-dumps")]
pub use read_filter_dump::{dump_hc_read_filter_tsv, HC_READ_FILTER_COUNT_SECTION};
#[cfg(feature = "parity_harness")]
pub use read_optional_tags::{
    format_optional_tag_field, parse_optional_tag_field, OptionalTagValue,
};
#[cfg(feature = "dev-dumps")]
pub use read_pre_len_dump::dump_read_pre_len_tsv;
#[cfg(feature = "dev-dumps")]
pub use read_pre_mq_dump::dump_read_pre_mq_tsv;
#[cfg(feature = "dev-dumps")]
pub use read_pre_overlap_dump::dump_read_pre_overlap_tsv;
#[cfg(feature = "dev-dumps")]
pub use read_transformer_dump::dump_read_shard_pipeline_tsv;
#[cfg(feature = "dev-dumps")]
pub use read_unclip_dump::dump_read_pre_softclip_tsv;
#[cfg(feature = "dev-dumps")]
pub use ref_confidence_dump::{
    dump_call_region_active_rcm_loci_tsv, dump_inactive_reference_model_tsv,
    dump_reference_confidence_locus_tsv,
};
#[cfg(feature = "dev-dumps")]
pub use ref_confidence_merger_dump::dump_ref_confidence_merge_case;
#[cfg(feature = "dev-dumps")]
pub use region_vcf_dump::{
    dump_call_region_format_tsv, dump_call_region_vcf_tsv, dump_variant_format_from_gl_ad_tsv,
    dump_variant_vcf_from_gl_ad_tsv,
};

#[cfg(test)]
mod tests {
    #[test]
    fn test_basic_functionality() {
        let sum = 1 + 1;
        assert_eq!(sum, 2);
    }
}
