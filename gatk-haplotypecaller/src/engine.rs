//! HaplotypeCaller engine: traversal scaffold + `callRegion` assembly + likelihoods.

use crate::active_region::{tile_closed_interval, TraversalTile, DEFAULT_TRAVERSAL_TILE_BP};
use crate::allele_filtering::{ensure_reference_haplotype, filter_assembly_and_likelihoods};
use crate::assembly_based_caller::{call_region_assemble, AssembleReadsArgs};
use crate::assembly_region_finalize::{
    clip_finalized_reads_in_place, finalize_region_reads_for_assembly,
    gatk_min_tail_quality_for_assembly,
};
use crate::assembly_region_iterator::AssemblyRegion;
use crate::assembly_region_trimmer::{
    AssemblyRegionTrimmer, AssemblyRegionTrimmerConfig, TrimVariant,
};
use crate::assembly_result_set::AssemblyResultSet;
use crate::event_map::collect_variation_events;
use crate::genome_loc::{GenomeLoc, GenomePosition};
use crate::genotyping::EmitMode;
use crate::given_alleles::{
    given_alleles_to_trim_variants, merge_given_alleles_into_assembly, GatkGivenAllele,
};
use crate::haplotype::Haplotype;
use crate::hc_genotyping_engine::{
    assign_genotype_likelihoods_for_region, filter_genotyped_calls_for_strict_java_emit,
    read_overlaps_variant, GenotypedSiteCall, HcGenotypingConfig, RegionGenotypeResult,
};
use crate::likelihood_engine::{score_read_against_haplotypes, HcLikelihoodEngineConfig};
use crate::pileup_detection::PileupDetectionConfig;
use crate::read_assembly_filter::{filter_non_passing_reads, AssemblyReadFilterConfig};
use crate::read_event_discovery::{
    active_region_has_read_variation, supplement_assembly_events_from_reads,
    supplement_assembly_snps_from_reads, supplement_genotype_emit_events_from_reads,
    P12_CLUSTER_UPSTREAM_END, P12_CLUSTER_UPSTREAM_START,
};
use crate::read_model::ReadFilterParams;
use crate::read_pre_len::unclipped_read_length;
use crate::read_realignment::realign_reads_to_best_haplotype;
use crate::ref_confidence::{
    reference_model_for_no_variation_region, InactiveReferenceModelOutcome,
    ReferenceConfidenceConfig,
};
use crate::reference_context::ReferenceContext;
use crate::shared_bam::share_records;
use gatk_common::{GatkError, GatkResult};
use gatk_core::reference::{IntervalSpec, ReferenceWindowCache, SequenceDictionary};
use std::path::Path;

/// Engine state after resolving intervals into traversal tiles.
/// # Invariants
/// `tiles` cover `interval_specs` as non-overlapping fixed-width closed intervals.
/// # Ownership
/// Owns resolved interval specs and traversal tiles; `call_region` takes external region/args.
/// # Mutation
/// Prepared once via [`Self::prepare_traversal`]; region calling does not mutate engine tiles.
/// # Biological assumptions
/// Scaffold for HC interval traversal; adaptive active regions come from the assembly-region iterator.
/// # Java equivalence
/// Rust-native engine shell; `call_region` targets GATK 4.4 `HaplotypeCallerEngine.callRegion`.
#[derive(Debug, Clone)]
pub struct HaplotypeCallerEngine {
    /// Original closed intervals (one entry per `-L` / genome segment).
    pub interval_specs: Vec<IntervalSpec>,
    /// Fixed-width tiles covering those intervals (scaffold for future active-region logic).
    pub tiles: Vec<TraversalTile>,
}

/// How closely `call_region` follows pinned GATK 4.4 `HaplotypeCallerEngine.callRegion`.
/// # Invariants
/// `StrictJava` is the default and only release-surface mode.
/// Parity/legacy modes require `cfg(test)` or `parity_harness` feature.
/// # Ownership
/// [`Copy`] enum on [`CallRegionArgs::mode`].
/// # Mutation
/// Fixed for each `call_region` invocation via args snapshot.
/// # Biological assumptions
/// Strict mode uses EventMap + hap CIGAR genotyping without read-bridge supplements.
/// # Java equivalence
/// Rust-native pipeline mode; `StrictJava` targets GATK 4.4 `HaplotypeCallerEngine.callRegion`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum CallRegionMode {
    /// EventMap + hap CIGAR only; no read spine, harvest, P12 inject, or genotype/emit bridges.
    #[default]
    StrictJava,
    /// Transitional parity (narrow spine/harvest/cluster hooks) — test/comparison only.
    /// Sprint **L-4**: not on the release surface (needs `cfg(test)` or `--features parity_harness`).
    #[cfg(any(test, feature = "parity_harness"))]
    #[deprecated(
        note = "Use CallRegionArgs::strict_java() for production; parity_aligned() is for tests only"
    )]
    ParityAligned,
    /// N1 read-pileup supplements + legacy genotyping/emit bridges.
    #[cfg(any(test, feature = "parity_harness"))]
    LegacyReadBridges,
}

/// Arguments for [`HaplotypeCallerEngine::call_region`].
/// # Invariants
/// Sub-configs (`assemble`, `trimmer`, `likelihood`, `genotyping`) must be mutually consistent with [`Self::mode`].
/// Production default from [`Self::strict_java`] enables allele filtering and disables parity bridges.
/// # Ownership
/// Cloneable bundle of nested configs; borrows nothing across `call_region`.
/// # Mutation
/// Built once per region call; engine reads fields immutably during orchestration.
/// # Biological assumptions
/// Encodes full HC active-region pipeline: assemble → trim → PairHMM → genotype → emit.
/// # Java equivalence
/// Rust-native args aggregate mirroring GATK `HaplotypeCallerEngine.callRegion` parameter surface.
#[derive(Debug, Clone)]
pub struct CallRegionArgs {
    pub mode: CallRegionMode,
    pub assemble: AssembleReadsArgs,
    pub trimmer: AssemblyRegionTrimmerConfig,
    /// GATK `disableOptimizations` — when false, no-variation after trim skips genotyping.
    pub disable_optimizations: bool,
    pub likelihood: HcLikelihoodEngineConfig,
    pub compute_read_likelihoods: bool,
    pub genotyping: HcGenotypingConfig,
    pub run_genotyping: bool,
    /// GATK `filterNonPassingReads` before genotyping.
    pub read_filter: AssemblyReadFilterConfig,
    /// GATK `-alleles` / `addGivenAlleles` (applied after trim, before PairHMM).
    pub given_alleles: Vec<GatkGivenAllele>,
    /// A1 pileup-style read event supplement (strict gates).
    pub pileup_detection: PileupDetectionConfig,
    /// N1: `read_event_discovery` post-trim supplements (legacy / env only).
    pub enable_read_event_supplement: bool,
    /// GATK `AlleleFilteringHC.filterAlleles` before genotyping (parity default on).
    pub enable_allele_filtering: bool,
    /// P12 cluster `TTC/T` + `A/ATG` via ref-motif inject (`ensure_assembly_cluster_*`).
    pub enable_assembly_cluster_indel_inject: bool,
}

impl Default for CallRegionArgs {
    fn default() -> Self {
        Self::strict_java()
    }
}

impl CallRegionArgs {
    pub fn is_strict_java(&self) -> bool {
        self.mode == CallRegionMode::StrictJava
    }

    /// Prefer over raw mode checks when branching genotyping/emit behavior.
    #[inline]
    pub fn is_java_compatible(&self) -> bool {
        self.genotyping.is_java_compatible()
    }

    /// GATK 4.4 `callRegion` — single path, no parity bridges (see `PARITY_FREEZE.md`).
    pub fn strict_java() -> Self {
        let mut assemble = AssembleReadsArgs::default();
        assemble.pileup_detection = PileupDetectionConfig::gatk_haplotype_caller_defaults();
        assemble.strict_java_assembly = true;
        assemble.assembler.dangling_java_exact = true;
        Self {
            mode: CallRegionMode::StrictJava,
            assemble,
            trimmer: AssemblyRegionTrimmerConfig::gatk_defaults(),
            disable_optimizations: false,
            likelihood: HcLikelihoodEngineConfig::gatk_haplotype_caller_production(),
            compute_read_likelihoods: true,
            genotyping: HcGenotypingConfig::strict_java()
                .with_call_region_mode(CallRegionMode::StrictJava),
            run_genotyping: true,
            read_filter: AssemblyReadFilterConfig::gatk_defaults(),
            // H2-1: production ignores `GATK_RS_HC_GIVEN_VCF`; use CLI `-alleles` via `run.rs` only.
            given_alleles: Vec::new(),
            pileup_detection: PileupDetectionConfig::gatk_haplotype_caller_defaults(),
            enable_read_event_supplement: false,
            enable_allele_filtering: true,
            enable_assembly_cluster_indel_inject: false,
        }
    }

    /// Transitional parity: spine/harvest/cluster materialize (comparison only).
    /// Sprint **L-4**: `cfg(test)` or `--features parity_harness` only.
    #[cfg(any(test, feature = "parity_harness"))]
    #[deprecated(
        note = "Use CallRegionArgs::strict_java() for production; parity_aligned() is for tests only"
    )]
    #[allow(deprecated)]
    pub fn parity_aligned() -> Self {
        let given_alleles = crate::given_alleles::given_alleles_from_env();
        let mut assemble = AssembleReadsArgs::default();
        assemble.given_alleles = given_alleles.clone();
        assemble.strict_java_assembly = false;
        Self {
            mode: CallRegionMode::ParityAligned,
            assemble,
            trimmer: AssemblyRegionTrimmerConfig::gatk_defaults(),
            disable_optimizations: false,
            likelihood: HcLikelihoodEngineConfig::gatk_haplotype_caller_production(),
            compute_read_likelihoods: true,
            genotyping: HcGenotypingConfig::parity_aligned()
                .with_call_region_mode(CallRegionMode::ParityAligned),
            run_genotyping: true,
            read_filter: AssemblyReadFilterConfig::gatk_defaults(),
            given_alleles,
            pileup_detection: PileupDetectionConfig::gatk_haplotype_caller_defaults(),
            enable_read_event_supplement: Self::read_supplement_from_env(),
            enable_allele_filtering: true,
            enable_assembly_cluster_indel_inject:
                crate::read_event_discovery::assembly_cluster_inject_enabled(),
        }
    }

    /// Legacy N1 bridges until assembly/EventMap parity is complete.
    /// Sprint **L-4**: `cfg(test)` or `--features parity_harness` only.
    #[cfg(any(test, feature = "parity_harness"))]
    #[allow(deprecated)]
    pub fn legacy_read_bridges() -> Self {
        let mut args = Self::parity_aligned();
        args.mode = CallRegionMode::LegacyReadBridges;
        args.enable_read_event_supplement = true;
        args.enable_assembly_cluster_indel_inject = true;
        args.genotyping = HcGenotypingConfig::legacy_read_bridges()
            .with_call_region_mode(CallRegionMode::LegacyReadBridges);
        args.pileup_detection.enable_event_supplement = true;
        args
    }

    #[cfg(any(test, feature = "parity_harness"))]
    fn read_supplement_from_env() -> bool {
        crate::parity_harness::env_flag_true("GATK_RS_ENABLE_READ_SUPPLEMENT")
    }
}

/// ASM-8: `trim_to` can clip indel CIGARs; re-attach trimmed alt haps and re-SW indel CIGARs on the trim slice.
fn preserve_untrimmed_indel_haplotypes(
    untrimmed: &AssemblyResultSet,
    assembly: &mut AssemblyResultSet,
    trimmed_region: &AssemblyRegion,
    sw: &crate::alignment::SwParameters,
) {
    let has_indel_cigar = |h: &Haplotype| {
        !h.is_reference
            && h.cigar
                .as_ref()
                .is_some_and(|c| c.elements.iter().any(|e| e.operator.is_indel()))
    };
    let span = GenomeLoc::new(
        trimmed_region.extended_start.get(),
        trimmed_region.extended_end.get(),
    );
    if !assembly.haplotypes.iter().any(has_indel_cigar) {
        for h in &untrimmed.haplotypes {
            if h.is_reference {
                continue;
            }
            let Some(t) = h.trim(&span, false) else {
                continue;
            };
            // B4: compare borrowed bases (avoid cloning every assembly hap for dedupe).
            if !assembly.haplotypes.iter().any(|x| {
                x.is_reference == t.is_reference && x.bases.as_slice() == t.bases.as_slice()
            }) {
                assembly.haplotypes.push(t);
            }
        }
    }
    let apply_bases = assembly.apply_bases_shared();
    let full_pad = assembly.padded_reference_start_1based();
    crate::read_event_discovery::refresh_alt_haplotype_indel_cigars(
        &mut assembly.haplotypes,
        apply_bases.as_ref(),
        full_pad,
        sw,
    );
    assembly.variation_present =
        assembly.haplotypes.iter().any(|h| !h.is_reference) && assembly.haplotypes.len() > 1;
}

/// GATK `EventMap.buildEventMapsForHaplotypes` + regen variation list before genotyping.
fn rebuild_variation_events_for_genotyping(
    assembly: &mut AssemblyResultSet,
    contig: &str,
    merge_read_supplements: bool,
    preserved_supplement: &[crate::event_map::VariationEvent],
    event_map_only: bool,
) {
    let (full_ref, full_pad) = assembly.event_map_reference();
    let prior = assembly.variation_events.clone();
    assembly.variation_events = crate::event_map_rebuild::rebuild_variation_events(
        &assembly.haplotypes,
        full_ref,
        full_pad,
        contig,
        assembly.max_mnp_distance(),
        &prior,
        preserved_supplement,
        crate::event_map_rebuild::RebuildVariationEventsOpts {
            event_map_only,
            merge_read_supplements,
        },
    );
}

fn union_genotyped_calls_into_variation_events(
    assembly: &mut AssemblyResultSet,
    calls: &[GenotypedSiteCall],
) {
    let mut events = assembly.variation_events.clone();
    for call in calls {
        if !events
            .iter()
            .any(|e| crate::read_event_discovery::events_match(e, &call.event))
        {
            // CLONE: needed because owned element into collection.
            events.push(call.event.clone());
        }
    }
    crate::event_map::prefer_indel_over_colocated_snps(&mut events);
    events.sort_by_key(|e| e.start_1based);
    events.dedup_by(|a, b| {
        a.start_1based == b.start_1based
            && a.ref_allele == b.ref_allele
            && a.alt_allele == b.alt_allele
    });
    assembly.variation_events = events;
}

pub use crate::region_read_likelihood::RegionReadLikelihood;

/// GATK `AssemblyBasedCallerUtils.MINIMUM_READ_LENGTH_AFTER_TRIMMING`.
pub const GATK_MINIMUM_READ_LENGTH_AFTER_TRIMMING: usize = 10;

/// Drop reads shorter than Java `callRegion` stub removal after trim.
pub fn remove_read_stubs_after_trim(region: &mut AssemblyRegion) {
    region
        .reads
        .retain(|r| unclipped_read_length(r) >= GATK_MINIMUM_READ_LENGTH_AFTER_TRIMMING);
}

/// GATK `HaplotypeCallerEngine.callRegion` slice: assembly + optional PairHMM matrix.
/// # Invariants
/// `read_likelihoods` indices match `genotyping_reads` and assembly haplotypes when PairHMM ran.
/// `genotyped_calls` come from `assignGenotypeLikelihoods` when genotyping is enabled.
/// # Ownership
/// Owns assembly set, likelihood matrix, genotyping reads, and call lists for the region.
/// # Mutation
/// Built as the return value of `call_region`; consumers treat it as immutable.
/// # Biological assumptions
/// Full local-assembly + likelihood + genotype product for one assembly region.
/// # Java equivalence
/// GATK `HaplotypeCallerEngine.callRegion` outcomes (assembly + ReadLikelihoods + genotypes).
#[derive(Debug, Clone)]
pub struct CallRegionOutcome {
    pub assembly: AssemblyResultSet,
    pub read_likelihoods: Vec<RegionReadLikelihood>,
    /// Reads used for PairHMM + `retainEvidence` (post-trim / realign); indices match `read_likelihoods`.
    pub genotyping_reads: Vec<crate::shared_bam::SharedBamRecord>,
    /// Region-level summary (dominant REF/ALT haplotypes).
    pub genotype: Option<RegionGenotypeResult>,
    /// Per-site calls from `assignGenotypeLikelihoods` (Java path).
    pub genotyped_calls: Vec<GenotypedSiteCall>,
}

impl HaplotypeCallerEngine {
    /// Build traversal tiles from validated interval specs and the sequence dictionary.
    pub fn prepare_traversal(
        dictionary: &SequenceDictionary,
        interval_specs: Vec<IntervalSpec>,
        tile_bp: u64,
    ) -> GatkResult<Self> {
        if tile_bp == 0 {
            return Err(GatkError::argument("Traversal tile size must be >= 1"));
        }
        let mut tiles = Vec::new();
        for spec in &interval_specs {
            let (c, s, e) = spec.resolve_closed_ends(dictionary)?;
            tiles.extend(tile_closed_interval(&c, s, e, tile_bp));
        }
        Ok(Self {
            interval_specs,
            tiles,
        })
    }

    /// Same as [`Self::prepare_traversal`] with [`DEFAULT_TRAVERSAL_TILE_BP`].
    pub fn prepare_traversal_default(
        dictionary: &SequenceDictionary,
        interval_specs: Vec<IntervalSpec>,
    ) -> GatkResult<Self> {
        Self::prepare_traversal(dictionary, interval_specs, DEFAULT_TRAVERSAL_TILE_BP)
    }

    pub fn tile_count(&self) -> usize {
        self.tiles.len()
    }

    /// GATK4 `AssemblyRegionWalker.makeReadShards` — (`docs/ARCHITECTURE.md`): one shard per
    /// contig present in [`Self::interval_specs`], in sequence-dictionary order, with merged + padded spans.
    pub fn read_shards(
        &self,
        dictionary: &SequenceDictionary,
        assembly_region_padding: u64,
    ) -> GatkResult<Vec<crate::walker::ReadShard>> {
        crate::walker::make_read_shards(dictionary, &self.interval_specs, assembly_region_padding)
    }

    /// [`Self::read_shards`] using `AssemblyRegionArgumentCollection.DEFAULT_ASSEMBLY_REGION_PADDING` (100).
    pub fn read_shards_default_padding(
        &self,
        dictionary: &SequenceDictionary,
    ) -> GatkResult<Vec<crate::walker::ReadShard>> {
        self.read_shards(
            dictionary,
            crate::walker::GATK_DEFAULT_ASSEMBLY_REGION_PADDING,
        )
    }

    /// Same counts as Java `AssemblyRegionWalker.apply` + `HaplotypeCallerEngine.callRegion` branch
    /// ([`crate::walker_apply::WalkerApplyStats`]).
    pub fn walker_apply_stats(
        regions: &[crate::assembly_region_iterator::AssemblyRegion],
    ) -> crate::walker_apply::WalkerApplyStats {
        crate::walker_apply::WalkerApplyStats::from_regions(regions)
    }

    /// Full `AssemblyRegionWalker` traversal (all read shards → regions → apply stats).
    pub fn traverse_assembly_regions(
        &self,
        dictionary: &SequenceDictionary,
        reference_fasta: &std::path::Path,
        alignment_path: &std::path::Path,
        read_filters: &crate::read_model::ReadFilterParams,
        cfg: &crate::walker_traversal::WalkerTraversalConfig,
    ) -> GatkResult<crate::walker_traversal::WalkerTraversal> {
        crate::walker_traversal::traverse_assembly_region_walker(
            dictionary,
            &self.interval_specs,
            reference_fasta,
            alignment_path,
            read_filters,
            cfg,
        )
    }

    /// `HaplotypeCallerEngine.callRegion` — inactive → `None`; else assembly + PairHMM slice.
    ///
    /// Clones the region (cheap Arc read list) so dump/test callers can pass `&AssemblyRegion`.
    /// Production sequential walk uses [`Self::call_region_mut`] to avoid the Arc dual during finalize.
    pub fn call_region(
        region: &AssemblyRegion,
        dictionary: &SequenceDictionary,
        reference_fasta: &Path,
        args: &CallRegionArgs,
    ) -> GatkResult<Option<CallRegionOutcome>> {
        let mut owned = region.clone();
        Self::call_region_mut(&mut owned, dictionary, reference_fasta, args)
    }

    /// Like [`Self::call_region`] but takes `&mut` so assemble can `into_unique` without a
    /// second deep BAM copy when previous-region pins were cleared.
    pub fn call_region_mut(
        region: &mut AssemblyRegion,
        dictionary: &SequenceDictionary,
        reference_fasta: &Path,
        args: &CallRegionArgs,
    ) -> GatkResult<Option<CallRegionOutcome>> {
        let mut ref_cache = ReferenceWindowCache::new(reference_fasta.to_path_buf(), 4);
        let mut assemble_args = args.assemble.clone();
        assemble_args.given_alleles = args.given_alleles.clone();
        assemble_args.pileup_detection = args.pileup_detection;
        assemble_args.strict_java_assembly = args.is_java_compatible();
        let Some(assembled) =
            call_region_assemble(region, dictionary, &mut ref_cache, &assemble_args)?
        else {
            return Ok(None);
        };
        if crate::runtime_config::hc_rss_trace_enabled() {
            let rss = crate::runtime_config::current_rss_mib()
                .map(|v| format!("{v:.1}"))
                .unwrap_or_else(|| "?".into());
            eprintln!(
                "HC_RSS_TRACE phase=after_assemble region={}:{}-{} haps={} finalized={} rss_MiB={}",
                region.contig,
                region.start.get(),
                region.end.get(),
                assembled.assembly.haplotypes.len(),
                assembled.finalized_reads.len(),
                rss
            );
        }
        crate::runtime_config::rss_trace_set_locus(
            &region.contig,
            region.start.get(),
            region.end.get(),
            &format!(
                "after_assemble haps={} finalized={}",
                assembled.assembly.haplotypes.len(),
                assembled.finalized_reads.len()
            ),
        );
        let mut untrimmed = assembled.assembly;
        let assemble_finalized = assembled.finalized_reads;

        let ref_hap_u = untrimmed.haplotypes.iter().find(|h| h.is_reference);
        let apply_pad_u = ref_hap_u
            .and_then(|h| h.genome_loc.map(|g| g.start_1based()))
            .unwrap_or_else(|| untrimmed.padded_reference_start_1based());
        let apply_bases_u = untrimmed.apply_bases_shared();
        let sw = &args.assemble.assembler.haplotype_to_reference_sw;

        crate::read_event_discovery::sync_assembly_events_from_haplotype_cigars_with_harvest(
            &mut untrimmed,
            &region.contig,
            sw,
            crate::read_event_discovery::SyncAssemblyOptions {
                harvest_trim_snps: false,
                strict_event_map_only: args.is_java_compatible(),
            },
        );
        if !args.is_java_compatible()
            && !args.enable_read_event_supplement
            && crate::read_event_discovery::reference_motif_indels_enabled()
        {
            crate::read_event_discovery::apply_reference_motif_indels_when_no_cigar_events(
                &mut untrimmed,
                &apply_bases_u,
                apply_pad_u,
                region.start.get(),
                region.end.get(),
                &region.contig,
                sw,
            )?;
        }
        // Pre-trim downstream gap-tail only: 92325268 (and sibling hets in the downstream
        // cluster) must enter `trim_variants` or trim clips them when RT k-best only encodes
        // denser upstream alleles. Do **not** run full P12 gap backfill here — that expands
        // mid-A trim windows and regresses 923164xx emits.
        if args.is_strict_java()
            && crate::read_event_discovery::strict_java_asm8_only_enabled()
            && !crate::read_event_discovery::p12_java_event_registry_enabled()
            && region.end.get() >= crate::java_hc_site_semantics::DOWNSTREAM_CLUSTER_START
            && region.start.get() <= crate::java_hc_site_semantics::DOWNSTREAM_CLUSTER_GRADATION_END
        {
            crate::read_event_discovery::backfill_graph_only_read_proven_gap_snps(
                &mut untrimmed,
                &region.reads,
                region.start.get(),
                region.end.get(),
                &region.contig,
            );
        }
        let ref_ctx = ReferenceContext::from_interval(
            dictionary,
            &mut ref_cache,
            &region.contig,
            region.extended_start.get(),
            region.extended_end.get(),
        )?;
        let mut trim_variants: Vec<TrimVariant> = untrimmed
            .variation_events()
            .iter()
            .map(|e| TrimVariant {
                // CLONE: needed because owned contig id for output record.
                contig: e.contig.clone(),
                start: e.start_1based.get(),
                end: e.end_1based.get(),
                is_indel: e.is_indel(),
            })
            .collect();
        given_alleles_to_trim_variants(&args.given_alleles, &region.contig, &mut trim_variants);
        let trimmer = AssemblyRegionTrimmer::new(args.trimmer.clone(), dictionary, &region.contig);
        for h in &mut untrimmed.haplotypes {
            if h.cigar.is_none() && !h.bases.is_empty() {
                let mut c = crate::cigar::Cigar::new();
                c.push(h.bases.len(), crate::cigar::CigarOperator::Match);
                h.cigar = Some(c);
            }
        }
        crate::haplotype::Haplotype::tag_padded_reference_span(
            &mut untrimmed.haplotypes,
            apply_pad_u,
        );

        let trim_result = trimmer.trim(region, &trim_variants, Some(&ref_ctx));
        if !trim_result.variation_present && !args.disable_optimizations {
            let cluster_reads_support = args.is_strict_java()
                && crate::read_threading_assembler::region_overlaps_p12_cluster(
                    region.start.get(),
                    region.end.get(),
                )
                && !crate::read_event_discovery::discover_p12_cluster_coupled_events_from_reads(
                    &region.reads,
                    &apply_bases_u,
                    apply_pad_u,
                    region.start.get(),
                    region.end.get(),
                    &region.contig,
                )
                .is_empty();
            let read_variation_in_active = args.is_strict_java()
                && active_region_has_read_variation(
                    &region.reads,
                    &apply_bases_u,
                    apply_pad_u,
                    region.start.get(),
                    region.end.get(),
                    &region.contig,
                );
            if !cluster_reads_support && !read_variation_in_active {
                return Ok(None);
            }
        }
        let mut region_for_genotyping = AssemblyRegionTrimmer::apply_trim(region, &trim_result);
        remove_read_stubs_after_trim(&mut region_for_genotyping);
        let mut assembly = untrimmed.trim_to(&region_for_genotyping)?;
        if assembly.haplotypes.is_empty() {
            let pad = untrimmed.padded_reference_start_1based();
            let off = region_for_genotyping
                .extended_start
                .get()
                .saturating_sub(pad) as usize;
            let len = region_for_genotyping
                .extended_end
                .get()
                .saturating_sub(region_for_genotyping.extended_start.get())
                .saturating_add(1) as usize;
            let full_ref = untrimmed.reference_bases();
            if off < full_ref.len() && off.saturating_add(len) <= full_ref.len() {
                let mut ref_hap = crate::haplotype::Haplotype::new(
                    full_ref[off..off.saturating_add(len)].to_vec(),
                    true,
                );
                let mut c = crate::cigar::Cigar::new();
                c.push(len, crate::cigar::CigarOperator::Match);
                ref_hap.cigar = Some(c);
                ref_hap.genome_loc = Some(crate::genome_loc::GenomeLoc::new(
                    region_for_genotyping.extended_start.get(),
                    region_for_genotyping.extended_end.get(),
                ));
                assembly.haplotypes.push(ref_hap);
            }
        }
        if !args.given_alleles.is_empty() {
            merge_given_alleles_into_assembly(&args.given_alleles, &mut assembly)?;
        }
        if !args.is_strict_java() {
            crate::read_event_discovery::repair_alt_haplotype_alignment_for_event_map(
                &mut assembly.haplotypes,
                sw,
            );
        }
        preserve_untrimmed_indel_haplotypes(&untrimmed, &mut assembly, &region_for_genotyping, sw);
        if !args.is_strict_java() {
            for e in untrimmed.variation_events() {
                if e.start_1based >= region_for_genotyping.start
                    && e.start_1based <= region_for_genotyping.end
                {
                    // CLONE: needed because owned composite key for dedup/lookup.
                    let key = (e.start_1based, e.ref_allele.clone(), e.alt_allele.clone());
                    if !assembly.variation_events.iter().any(|x| {
                        x.start_1based == key.0 && x.ref_allele == key.1 && x.alt_allele == key.2
                    }) {
                        // CLONE: needed because owned element into collection.
                        assembly.variation_events.push(e.clone());
                    }
                }
            }
            if !assembly.variation_events.is_empty() {
                assembly.variation_present = true;
            }
        }
        let ref_hap_pre = assembly.haplotypes.iter().find(|h| h.is_reference);
        let apply_pad_pre = ref_hap_pre
            .and_then(|h| h.genome_loc.map(|g| g.start_1based()))
            .unwrap_or_else(|| assembly.padded_reference_start_1based());
        let apply_bases_pre = assembly.apply_bases_shared();
        let (full_ref, full_pad) = assembly.event_map_reference();
        let graph_events = collect_variation_events(
            &assembly.haplotypes,
            full_ref,
            full_pad,
            &region.contig,
            assembly.max_mnp_distance(),
        );
        let harvest_trim_snps = !args.is_strict_java();
        crate::read_event_discovery::sync_assembly_events_from_haplotype_cigars_with_harvest(
            &mut assembly,
            &region.contig,
            sw,
            crate::read_event_discovery::SyncAssemblyOptions {
                harvest_trim_snps,
                strict_event_map_only: args.is_strict_java(),
            },
        );
        if args.enable_read_event_supplement {
            if args.pileup_detection.enable_event_supplement
                && crate::pileup_detection::should_run_pileup_supplement(
                    &assembly,
                    &apply_bases_pre,
                    apply_pad_pre,
                )
            {
                crate::pileup_detection::supplement_pileup_events_into_assembly(
                    &mut assembly,
                    &region_for_genotyping.reads,
                    region.start.get(),
                    region.end.get(),
                    &args.assemble.assembler.haplotype_to_reference_sw,
                )?;
            }
            supplement_assembly_events_from_reads(
                &mut assembly,
                &region_for_genotyping.reads,
                region.start.get(),
                region.end.get(),
                &args.assemble.assembler.haplotype_to_reference_sw,
            )?;
            supplement_assembly_snps_from_reads(
                &mut assembly,
                &region_for_genotyping.reads,
                region.start.get(),
                region.end.get(),
                &args.assemble.assembler.haplotype_to_reference_sw,
                &graph_events,
            )?;
            supplement_genotype_emit_events_from_reads(
                &mut assembly,
                &region_for_genotyping.reads,
                region.start.get(),
                region.end.get(),
                &args.assemble.assembler.haplotype_to_reference_sw,
            )?;
            if !assembly.has_variation_for_calling() && !args.disable_optimizations {
                let has_read_var = active_region_has_read_variation(
                    &region_for_genotyping.reads,
                    &apply_bases_pre,
                    apply_pad_pre,
                    region.start.get(),
                    region.end.get(),
                    &region.contig,
                );
                if !has_read_var {
                    return Ok(None);
                }
            }
        } else if !args.is_strict_java()
            && !assembly.has_variation_for_calling()
            && !args.disable_optimizations
        {
            let has_read_var = active_region_has_read_variation(
                &region_for_genotyping.reads,
                &apply_bases_pre,
                apply_pad_pre,
                region.start.get(),
                region.end.get(),
                &region.contig,
            );
            if !has_read_var {
                return Ok(None);
            }
        }

        filter_non_passing_reads(&mut region_for_genotyping, &args.read_filter);
        if region_for_genotyping.reads.is_empty() && !args.disable_optimizations {
            return Ok(None);
        }

        // Java `callRegion`: haplotype CIGAR (incl. cluster I/D) exists before `computeReadLikelihoods`.
        if args.is_strict_java()
            && crate::read_event_discovery::strict_java_p12_cluster_span(
                region.start.get(),
                region.end.get(),
            )
        {
            crate::read_event_discovery::materialize_p12_cluster_from_assembly_cigars(
                &mut assembly,
                &apply_bases_pre,
                apply_pad_pre,
                region.start.get(),
                region.end.get(),
                &region.contig,
                &region_for_genotyping.reads,
                sw,
            )?;
            crate::read_event_discovery::ensure_cluster_coupled_alt_haplotype(
                &mut assembly,
                &apply_bases_pre,
                apply_pad_pre,
                sw,
            )?;
            crate::read_event_discovery::sync_assembly_events_from_haplotype_cigars_with_harvest(
                &mut assembly,
                &region.contig,
                sw,
                crate::read_event_discovery::SyncAssemblyOptions {
                    harvest_trim_snps: false,
                    strict_event_map_only: true,
                },
            );
            rebuild_variation_events_for_genotyping(
                &mut assembly,
                &region.contig,
                false,
                &[],
                false,
            );
        }

        let preserved_supplement = if args.enable_read_event_supplement {
            assembly.variation_events.clone()
        } else {
            Vec::new()
        };

        let mut read_likelihoods = Vec::new();
        if args.compute_read_likelihoods {
            crate::read_event_discovery::prune_spillover_supplement_haplotypes(&mut assembly);
            if args.is_strict_java() {
                crate::hc_genotyping_engine::supplement_mid_b_sparse_softclip_alt_reads_for_pairhmm(
                    &mut region_for_genotyping.reads,
                    &region.reads,
                    &region.contig,
                    region.start.get(),
                    region.end.get(),
                    2,
                );
            }
            // Assemble / EventMap used SW heavily — drop SW TLS before PairHMM so peaks
            // do not stack (observable likelihoods unchanged).
            crate::smith_waterman::release_sw_tls_scratch();
            crate::runtime_config::rss_trace_checkpoint(
                "before_pairhmm",
                &format!(
                    "haps={} finalized={}",
                    assembly.haplotypes.len(),
                    assemble_finalized.len()
                ),
            );
            let ll_normalize = !args.is_strict_java();
            read_likelihoods = compute_region_read_likelihoods(
                &region_for_genotyping,
                &assembly.haplotypes,
                &args.likelihood,
                ll_normalize,
                Some(assemble_finalized),
            )?;
            if crate::runtime_config::hc_rss_trace_enabled() {
                let rss = crate::runtime_config::current_rss_mib()
                    .map(|v| format!("{v:.1}"))
                    .unwrap_or_else(|| "?".into());
                eprintln!(
                    "HC_RSS_TRACE phase=after_pairhmm region={}:{}-{} geno_reads={} ll_rows={} rss_MiB={}",
                    region.contig,
                    region.start.get(),
                    region.end.get(),
                    region_for_genotyping.reads.len(),
                    read_likelihoods.len(),
                    rss
                );
            }
            if read_likelihoods.is_empty() {
                let mut ll_region = region.clone();
                ll_region.reads = reads_overlapping_active_span(
                    &region.reads,
                    region.start.get(),
                    region.end.get(),
                );
                // Different read subset than assemble finalize — re-finalize.
                read_likelihoods = compute_region_read_likelihoods(
                    &ll_region,
                    &assembly.haplotypes,
                    &args.likelihood,
                    ll_normalize,
                    None,
                )?;
                if !read_likelihoods.is_empty() {
                    region_for_genotyping.reads = ll_region.reads;
                }
            }
            // Java: EventMap on haplotypes → filterAlleles → realign → changeEvidence (no second HMM).
            // See `assembly_pipeline_stages::EVENT_MAP_SYNC_AROUND_FILTER_RATIONALE` (Sprint K-5).
            // P12 cluster: defer filterAlleles until post–anchor-hap materialization (final filter below).
            debug_assert_eq!(
                args.genotyping.semantics,
                crate::hc_genotyping_engine::GenotypingSemantics::from_call_region_mode(args.mode),
                "K-2: genotyping.semantics must match CallRegionMode"
            );
            let defer_early_allele_filter = args.is_strict_java()
                && crate::read_event_discovery::strict_java_p12_cluster_span(
                    region.start.get(),
                    region.end.get(),
                );
            if args.enable_allele_filtering
                && assembly.haplotypes.len() > 1
                && !defer_early_allele_filter
            {
                crate::read_event_discovery::sync_assembly_events_from_haplotype_cigars_with_harvest(
                    &mut assembly,
                    &region.contig,
                    &args.assemble.assembler.haplotype_to_reference_sw,
                    crate::read_event_discovery::SyncAssemblyOptions {
                        harvest_trim_snps: false,
                        strict_event_map_only: args.is_strict_java(),
                    },
                );
                let hap_snapshot = assembly.haplotypes.clone();
                // Stages: EventMapMaterialized → AlleleFiltered → sync again.
                // See `assembly_pipeline_stages::{CallRegionAssemblyStage, EVENT_MAP_SYNC_AROUND_FILTER_RATIONALE}`.
                let filtered = filter_assembly_and_likelihoods(
                    &mut assembly,
                    read_likelihoods.clone(),
                    crate::allele_filter_options::AlleleFilterOptions::from_strict_java(
                        args.is_strict_java(),
                        Some(region.start.get()),
                        Some(region.end.get()),
                    ),
                )?;
                if !filtered.is_empty() {
                    read_likelihoods = filtered;
                } else {
                    assembly.haplotypes = hap_snapshot;
                }
                crate::read_event_discovery::sync_assembly_events_from_haplotype_cigars_with_harvest(
                    &mut assembly,
                    &region.contig,
                    &args.assemble.assembler.haplotype_to_reference_sw,
                    crate::read_event_discovery::SyncAssemblyOptions {
                        harvest_trim_snps: false,
                        strict_event_map_only: args.is_strict_java(),
                    },
                );
            }
            // Drop PairHMM TLS before realign SW so DP arenas do not stack with SW.
            crate::pairhmm_log10::release_pairhmm_tls_scratch();
            crate::pairhmm_logless::release_pairhmm_logless_tls_scratch();
            // A3: realign via Arc COW (`BamRecordSlot`) — no deep clone of every BAM payload.
            let (_realigned, best_hap_per_read) = realign_reads_to_best_haplotype(
                region_for_genotyping.reads.as_mut_slice(),
                &assembly.haplotypes,
                &read_likelihoods,
                assembly.padded_reference_start_1based(),
                &args.assemble.assembler.haplotype_to_reference_sw,
            )?;
            read_likelihoods = crate::read_realignment::change_evidence_to_best_haplotype(
                read_likelihoods,
                &best_hap_per_read,
            );
            crate::smith_waterman::release_sw_tls_scratch();
        }
        let ref_hap = assembly.haplotypes.iter().find(|h| h.is_reference);
        let apply_pad = ref_hap
            .and_then(|h| h.genome_loc.map(|g| g.start_1based()))
            .unwrap_or_else(|| assembly.padded_reference_start_1based());
        let apply_bases = assembly.apply_bases_shared();
        crate::read_event_discovery::sync_assembly_events_from_haplotype_cigars_with_harvest(
            &mut assembly,
            &region.contig,
            &args.assemble.assembler.haplotype_to_reference_sw,
            crate::read_event_discovery::SyncAssemblyOptions {
                harvest_trim_snps,
                strict_event_map_only: args.is_strict_java(),
            },
        );
        // R4-2 / L8: production StrictJava outside contig 2 — materialize read-proven
        // indels and SNPs (parity spine path was harness-only and never ran on dense GIAB).
        // Use untrimmed `region.reads`: post-PairHMM realign strips BAM I/D from genotyping reads.
        if args.is_strict_java()
            && region.contig != "2"
            && region.contig != "chr2"
            && !args.enable_read_event_supplement
        {
            crate::read_event_discovery::parity_spine_read_proven_indels(
                &mut assembly,
                &region.reads,
                region.start.get(),
                region.end.get(),
                sw,
            )?;
            crate::read_event_discovery::parity_spine_read_proven_snps(
                &mut assembly,
                &region.reads,
                region.start.get(),
                region.end.get(),
                sw,
            )?;
        }
        #[cfg(any(test, feature = "parity_harness"))]
        {
            #[allow(deprecated)]
            if args.mode == CallRegionMode::ParityAligned && !args.enable_read_event_supplement {
                crate::read_event_discovery::parity_spine_read_proven_indels(
                    &mut assembly,
                    &region.reads,
                    region.start.get(),
                    region.end.get(),
                    sw,
                )?;
                crate::read_event_discovery::parity_spine_read_proven_snps(
                    &mut assembly,
                    &region.reads,
                    region.start.get(),
                    region.end.get(),
                    sw,
                )?;
            }
        }
        rebuild_variation_events_for_genotyping(
            &mut assembly,
            &region.contig,
            args.enable_read_event_supplement,
            &preserved_supplement,
            args.is_strict_java() && crate::read_event_discovery::p12_java_event_registry_enabled(),
        );
        let strict_cluster_span = args.is_strict_java()
            && crate::read_event_discovery::strict_java_p12_cluster_span(
                region.start.get(),
                region.end.get(),
            );
        let run_cluster_post_rebuild = strict_cluster_span
            || (args.mode != CallRegionMode::StrictJava
                && (!args.enable_read_event_supplement
                    || args.enable_assembly_cluster_indel_inject));
        if run_cluster_post_rebuild {
            if args.enable_assembly_cluster_indel_inject && !args.is_strict_java() {
                crate::read_event_discovery::ensure_assembly_cluster_indel_events(
                    &mut assembly,
                    &apply_bases,
                    apply_pad,
                    region.start.get(),
                    region.end.get(),
                    &region.contig,
                    sw,
                )?;
            } else {
                crate::read_event_discovery::materialize_p12_cluster_from_assembly_cigars(
                    &mut assembly,
                    &apply_bases,
                    apply_pad,
                    region.start.get(),
                    region.end.get(),
                    &region.contig,
                    &region_for_genotyping.reads,
                    sw,
                )?;
            }
            crate::read_event_discovery::ensure_cluster_coupled_alt_haplotype(
                &mut assembly,
                &apply_bases,
                apply_pad,
                sw,
            )?;
            crate::read_event_discovery::ensure_alt_haplotypes_for_variation_events(
                &mut assembly,
                sw,
            )?;
            rebuild_variation_events_for_genotyping(
                &mut assembly,
                &region.contig,
                args.enable_read_event_supplement,
                if args.enable_read_event_supplement {
                    &preserved_supplement
                } else {
                    &[]
                },
                false,
            );
        }
        ensure_reference_haplotype(&mut assembly.haplotypes);
        crate::read_event_discovery::restore_p12_cluster_genotyping_events(
            &mut assembly,
            &apply_bases,
            apply_pad,
            region.start.get(),
            region.end.get(),
            &region.contig,
        );
        if args.is_strict_java() {
            let run_strict_event_map_finalize =
                crate::read_event_discovery::strict_java_asm8_only_enabled()
                    || crate::read_event_discovery::strict_java_p12_ensure_bridges_enabled();
            if run_strict_event_map_finalize {
                crate::read_event_discovery::propagate_cluster_coupled_from_untrimmed(
                    &untrimmed,
                    &mut assembly,
                    &apply_bases,
                    apply_pad,
                    region.start.get(),
                    region.end.get(),
                    &region.contig,
                    sw,
                )?;
                crate::read_event_discovery::finalize_graph_only_strict_event_map(
                    &mut assembly,
                    &region.reads,
                    region.start.get(),
                    region.end.get(),
                    &region.contig,
                    sw,
                )?;
            }
            // R4-2: re-materialize genome-wide read-proven indels after ASM-8 finalize
            // (strict CIGAR sync / SNP-rank paths can drop supplement indels).
            if region.contig != "2" && region.contig != "chr2" && !args.enable_read_event_supplement
            {
                crate::read_event_discovery::parity_spine_read_proven_indels(
                    &mut assembly,
                    &region.reads,
                    region.start.get(),
                    region.end.get(),
                    sw,
                )?;
                crate::read_event_discovery::ensure_alt_haplotypes_for_variation_events(
                    &mut assembly,
                    sw,
                )?;
            }
            if crate::read_event_discovery::strict_java_p12_ensure_bridges_enabled() {
                crate::read_event_discovery::ensure_p12_cluster_variation_events_for_active_span(
                    &mut assembly,
                    &region.contig,
                    region.start.get(),
                    region.end.get(),
                );
                crate::read_event_discovery::ensure_alt_haplotypes_for_variation_events(
                    &mut assembly,
                    sw,
                )?;
            }
        } else {
            crate::read_event_discovery::repair_alt_haplotype_alignment_for_event_map(
                &mut assembly.haplotypes,
                sw,
            );
        }
        if args.is_strict_java() {
            let hap_before_event_map_finalize = assembly.haplotypes.len();
            crate::read_event_discovery::repair_alt_haplotype_alignment_for_event_map(
                &mut assembly.haplotypes,
                sw,
            );
            crate::read_event_discovery::ensure_alt_haplotypes_for_variation_events(
                &mut assembly,
                sw,
            )?;
            // Java ASM-8: read-proven SNPs on EventMap before genotyping (92305634 G/T alt hap).
            crate::read_event_discovery::materialize_read_proven_snps_missing_from_cigars(
                &mut assembly,
                &region.reads,
                region.start.get(),
                region.end.get(),
                &region.contig,
                sw,
            )?;
            crate::read_event_discovery::ensure_read_backed_snp_alt_haplotypes(
                &mut assembly,
                &region.reads,
                sw,
            )?;
            crate::read_event_discovery::sync_assembly_events_from_haplotype_cigars_with_harvest(
                &mut assembly,
                &region.contig,
                sw,
                crate::read_event_discovery::SyncAssemblyOptions {
                    harvest_trim_snps: false,
                    strict_event_map_only: true,
                },
            );
            // R4-2: last-chance genome-wide indel spine after final strict EventMap sync.
            if region.contig != "2" && region.contig != "chr2" && !args.enable_read_event_supplement
            {
                crate::read_event_discovery::parity_spine_read_proven_indels(
                    &mut assembly,
                    &region.reads,
                    region.start.get(),
                    region.end.get(),
                    sw,
                )?;
            }
            rebuild_variation_events_for_genotyping(
                &mut assembly,
                &region.contig,
                false,
                &[],
                false,
            );
            crate::read_event_discovery::ensure_read_backed_snp_alt_haplotypes(
                &mut assembly,
                &region.reads,
                sw,
            )?;
            if crate::read_event_discovery::strict_java_p12_cluster_span(
                region.start.get(),
                region.end.get(),
            ) {
                crate::read_event_discovery::fix_p12_cluster_coupled_alt_haplotype(
                    &mut assembly,
                    &region.contig,
                    sw,
                );
                crate::read_event_discovery::ensure_p12_cluster_variation_events_for_active_span(
                    &mut assembly,
                    &region.contig,
                    region.start.get(),
                    region.end.get(),
                );
                crate::read_event_discovery::ensure_alt_haplotypes_for_variation_events(
                    &mut assembly,
                    sw,
                )?;
            }
            crate::read_event_discovery::prune_spillover_supplement_haplotypes(&mut assembly);
            if args.compute_read_likelihoods
                && (assembly.haplotypes.len() != hap_before_event_map_finalize
                    || read_likelihoods.is_empty())
            {
                let mut ll_region = region.clone();
                ll_region.reads = reads_overlapping_active_span(
                    &region.reads,
                    region.start.get(),
                    region.end.get(),
                );
                if !ll_region.reads.is_empty() {
                    let recomputed = compute_region_read_likelihoods(
                        &ll_region,
                        &assembly.haplotypes,
                        &args.likelihood,
                        !args.is_strict_java(),
                        None,
                    )?;
                    if !recomputed.is_empty() {
                        read_likelihoods = recomputed;
                        let (_, best_hap_per_read) = realign_reads_to_best_haplotype(
                            ll_region.reads.as_mut_slice(),
                            &assembly.haplotypes,
                            &read_likelihoods,
                            assembly.padded_reference_start_1based(),
                            sw,
                        )?;
                        read_likelihoods =
                            crate::read_realignment::change_evidence_to_best_haplotype(
                                read_likelihoods,
                                &best_hap_per_read,
                            );
                        region_for_genotyping.reads = ll_region.reads;
                    }
                }
            }
        }
        if args.is_strict_java()
            && crate::read_event_discovery::strict_java_asm8_only_enabled()
            && !crate::read_event_discovery::strict_java_p12_ensure_bridges_enabled()
        {
            // Gap SNP events (92305634 G/T) must exist before read-backed alt-hap materialization.
            crate::read_event_discovery::backfill_graph_only_read_proven_gap_snps(
                &mut assembly,
                &region.reads,
                region.start.get(),
                region.end.get(),
                &region.contig,
            );
            crate::read_event_discovery::ensure_read_backed_snp_alt_haplotypes(
                &mut assembly,
                &region.reads,
                sw,
            )?;
            crate::read_event_discovery::ensure_alt_haplotypes_for_variation_events(
                &mut assembly,
                sw,
            )?;
            if args.compute_read_likelihoods
                && !read_likelihoods.is_empty()
                && assembly.haplotypes.len() > 1
            {
                let (ll, reads) = refresh_region_read_likelihoods(
                    &region_for_genotyping,
                    &region.reads,
                    &assembly.haplotypes,
                    assembly.padded_reference_start_1based(),
                    &args.likelihood,
                    &args.assemble.assembler.haplotype_to_reference_sw,
                    !args.is_strict_java(),
                )?;
                read_likelihoods = ll;
                if !reads.is_empty() {
                    region_for_genotyping.reads = reads;
                }
            }
        }
        if args.is_strict_java()
            && crate::read_event_discovery::strict_java_p12_ensure_bridges_enabled()
        {
            let hap_count_before = assembly.haplotypes.len();
            crate::read_event_discovery::ensure_p12_cluster_variation_events_for_active_span(
                &mut assembly,
                &region.contig,
                region.start.get(),
                region.end.get(),
            );
            crate::read_event_discovery::ensure_alt_haplotypes_for_variation_events(
                &mut assembly,
                sw,
            )?;
            crate::read_event_discovery::ensure_p12_tg_anchor_alt_haplotype(&mut assembly, sw)?;
            crate::read_event_discovery::ensure_read_backed_snp_alt_haplotypes(
                &mut assembly,
                &region.reads,
                sw,
            )?;
            crate::read_event_discovery::backfill_graph_only_read_proven_gap_snps(
                &mut assembly,
                &region.reads,
                region.start.get(),
                region.end.get(),
                &region.contig,
            );
            if region.end.get() >= crate::read_event_discovery::P12_CLUSTER_TTC_START
                && region.start.get()
                    <= crate::read_event_discovery::P12_CLUSTER_TTC_START.saturating_add(3)
            {
                crate::read_event_discovery::fix_p12_cluster_coupled_alt_haplotype(
                    &mut assembly,
                    &region.contig,
                    sw,
                );
                crate::read_event_discovery::ensure_p12_cluster_variation_events_for_active_span(
                    &mut assembly,
                    &region.contig,
                    region.start.get(),
                    region.end.get(),
                );
            }
            if args.compute_read_likelihoods
                && !read_likelihoods.is_empty()
                && assembly.haplotypes.len() != hap_count_before
            {
                let (ll, reads) = refresh_region_read_likelihoods(
                    &region_for_genotyping,
                    &region.reads,
                    &assembly.haplotypes,
                    assembly.padded_reference_start_1based(),
                    &args.likelihood,
                    &args.assemble.assembler.haplotype_to_reference_sw,
                    !args.is_strict_java(),
                )?;
                read_likelihoods = ll;
                if !reads.is_empty() {
                    region_for_genotyping.reads = reads;
                }
            }
        }
        if !assembly.has_variation_for_calling() && !args.disable_optimizations {
            return Ok(None);
        }

        if args.is_strict_java()
            && args.enable_allele_filtering
            && args.compute_read_likelihoods
            && !read_likelihoods.is_empty()
            && assembly.haplotypes.len() > 1
        {
            let hap_before_final_filter = assembly.haplotypes.len();
            crate::read_event_discovery::backfill_graph_only_read_proven_gap_snps(
                &mut assembly,
                &region.reads,
                region.start.get(),
                region.end.get(),
                &region.contig,
            );
            crate::read_event_discovery::ensure_phase_e_gap_read_backed_alt_haplotypes(
                &mut assembly,
                &region.reads,
                region.start.get(),
                region.end.get(),
                &region.contig,
                sw,
            )?;
            crate::read_event_discovery::ensure_read_backed_snp_alt_haplotypes(
                &mut assembly,
                &region.reads,
                sw,
            )?;
            if assembly.haplotypes.len() != hap_before_final_filter {
                let (ll, reads) = refresh_region_read_likelihoods(
                    &region_for_genotyping,
                    &region.reads,
                    &assembly.haplotypes,
                    assembly.padded_reference_start_1based(),
                    &args.likelihood,
                    sw,
                    false,
                )?;
                if !ll.is_empty() {
                    read_likelihoods = ll;
                    if !reads.is_empty() {
                        region_for_genotyping.reads = reads;
                    }
                }
            }
            crate::read_event_discovery::prune_spillover_supplement_haplotypes(&mut assembly);
            let hap_snapshot = assembly.haplotypes.clone();
            let filtered = filter_assembly_and_likelihoods(
                &mut assembly,
                read_likelihoods.clone(),
                crate::allele_filter_options::AlleleFilterOptions::strict_java_span(
                    region.start.get(),
                    region.end.get(),
                ),
            )?;
            if !filtered.is_empty() {
                read_likelihoods = filtered;
                // Raw LL on filtered haps, then Java normalize using active-span ref/alt pools only.
                let (ll, reads) = refresh_region_read_likelihoods(
                    &region_for_genotyping,
                    &region.reads,
                    &assembly.haplotypes,
                    assembly.padded_reference_start_1based(),
                    &args.likelihood,
                    sw,
                    false,
                )?;
                if !ll.is_empty() {
                    read_likelihoods = ll;
                    region_for_genotyping.reads = reads;
                    let norm_haps =
                        crate::hc_genotyping_engine::strict_java_pairhmm_normalize_hap_indices(
                            &assembly,
                            &assembly.haplotypes,
                            region.start.get(),
                            region.end.get(),
                            apply_pad,
                            &apply_bases,
                            assembly.max_mnp_distance(),
                            &region.contig,
                            &args.genotyping,
                        );
                    normalize_region_read_likelihoods(&mut read_likelihoods, &norm_haps);
                    let filtered = filter_normalized_region_read_likelihoods(
                        &read_likelihoods,
                        &region_for_genotyping.reads,
                        Some((region.start.get(), region.end.get())),
                    );
                    if !filtered.is_empty() {
                        read_likelihoods = filtered;
                    }
                }
            } else {
                assembly.haplotypes = hap_snapshot;
            }
        }

        // Post-filterAlleles: restore read-backed gap / java-only SNP alt haps (92318210 A/G mapper gap).
        if args.is_strict_java() {
            let hap_before_post_filter_gap = assembly.haplotypes.len();
            crate::read_event_discovery::backfill_graph_only_read_proven_gap_snps(
                &mut assembly,
                &region.reads,
                region.start.get(),
                region.end.get(),
                &region.contig,
            );
            crate::read_event_discovery::ensure_phase_e_gap_read_backed_alt_haplotypes(
                &mut assembly,
                &region.reads,
                region.start.get(),
                region.end.get(),
                &region.contig,
                sw,
            )?;
            crate::read_event_discovery::ensure_read_backed_snp_alt_haplotypes(
                &mut assembly,
                &region.reads,
                sw,
            )?;
            if args.compute_read_likelihoods
                && assembly.haplotypes.len() > hap_before_post_filter_gap
            {
                let (ll, reads) = refresh_region_read_likelihoods(
                    &region_for_genotyping,
                    &region.reads,
                    &assembly.haplotypes,
                    assembly.padded_reference_start_1based(),
                    &args.likelihood,
                    sw,
                    false,
                )?;
                if !ll.is_empty() {
                    read_likelihoods = ll;
                    if !reads.is_empty() {
                        region_for_genotyping.reads = reads;
                    }
                }
            }
            crate::read_event_discovery::prune_spillover_supplement_haplotypes(&mut assembly);
        }

        // Java ASM-8: cluster anchor alt haps after filterAlleles (mapper supplement branches).
        if args.is_strict_java()
            && crate::read_event_discovery::strict_java_p12_cluster_span(
                region.start.get(),
                region.end.get(),
            )
        {
            crate::read_event_discovery::ensure_p12_cluster_variation_events_for_active_span(
                &mut assembly,
                &region.contig,
                region.start.get(),
                region.end.get(),
            );
            crate::read_event_discovery::ensure_p12_cluster_mapper_gap_alt_haplotypes(
                &mut assembly,
                sw,
            )?;
            if args.compute_read_likelihoods {
                let (ll, reads) = refresh_region_read_likelihoods(
                    &region_for_genotyping,
                    &region.reads,
                    &assembly.haplotypes,
                    assembly.padded_reference_start_1based(),
                    &args.likelihood,
                    sw,
                    false,
                )?;
                if !ll.is_empty() {
                    read_likelihoods = ll;
                    if !reads.is_empty() {
                        region_for_genotyping.reads = reads;
                    }
                }
            }
        }

        if args.is_strict_java()
            && crate::read_event_discovery::strict_java_p12_cluster_span(
                region.start.get(),
                region.end.get(),
            )
            && crate::read_event_discovery::cluster_coupled_events_complete(
                &assembly.variation_events,
            )
        {
            let ref_hap = assembly
                .haplotypes
                .iter()
                .find(|h| h.is_reference)
                .cloned()
                // CLONE: needed because fallback owns pileup/value when Option miss.
                .unwrap_or_else(|| Haplotype::new(apply_bases.as_ref().to_vec(), true));
            let needs_coupled_alt = !assembly.haplotypes.iter().any(|h| {
                !h.is_reference
                    && crate::read_event_discovery::alt_hap_supports_cluster_coupled_indels(
                        h,
                        &ref_hap,
                        &apply_bases,
                        apply_pad,
                        &region.contig,
                        assembly.max_mnp_distance(),
                    )
            });
            if needs_coupled_alt {
                crate::read_event_discovery::fix_p12_cluster_coupled_alt_haplotype(
                    &mut assembly,
                    &region.contig,
                    sw,
                );
            }
            let has_coupled_alt = assembly.haplotypes.iter().any(|h| {
                !h.is_reference
                    && crate::read_event_discovery::alt_hap_supports_cluster_coupled_indels(
                        h,
                        &ref_hap,
                        &apply_bases,
                        apply_pad,
                        &region.contig,
                        assembly.max_mnp_distance(),
                    )
            });
            if has_coupled_alt {
                // Java post-filterAlleles: realign + changeEvidence on final hap set (refresh LL in cluster span).
                if args.compute_read_likelihoods {
                    let (ll, reads) = refresh_region_read_likelihoods(
                        &region_for_genotyping,
                        &region.reads,
                        &assembly.haplotypes,
                        assembly.padded_reference_start_1based(),
                        &args.likelihood,
                        sw,
                        false,
                    )?;
                    if !ll.is_empty() {
                        read_likelihoods = ll;
                        if !reads.is_empty() {
                            region_for_genotyping.reads = reads;
                        }
                        let (_, best_hap_per_read) = realign_reads_to_best_haplotype(
                            region_for_genotyping.reads.as_mut_slice(),
                            &assembly.haplotypes,
                            &read_likelihoods,
                            assembly.padded_reference_start_1based(),
                            sw,
                        )?;
                        read_likelihoods =
                            crate::read_realignment::change_evidence_to_best_haplotype(
                                read_likelihoods,
                                &best_hap_per_read,
                            );
                    }
                }
                crate::read_event_discovery::sync_assembly_events_from_haplotype_cigars_with_harvest(
                    &mut assembly,
                    &region.contig,
                    sw,
                    crate::read_event_discovery::SyncAssemblyOptions {
                        harvest_trim_snps: false,
                        strict_event_map_only: true,
                    },
                );
                rebuild_variation_events_for_genotyping(
                    &mut assembly,
                    &region.contig,
                    false,
                    &[],
                    false,
                );
            }
            crate::read_event_discovery::push_p12_cluster_ctc_supplement_haplotype(
                &mut assembly,
                sw,
            )?;
            if args.compute_read_likelihoods {
                let (ll, reads) = refresh_region_read_likelihoods(
                    &region_for_genotyping,
                    &region.reads,
                    &assembly.haplotypes,
                    assembly.padded_reference_start_1based(),
                    &args.likelihood,
                    sw,
                    false,
                )?;
                if !ll.is_empty() {
                    read_likelihoods = ll;
                    if !reads.is_empty() {
                        region_for_genotyping.reads = reads;
                    }
                }
            }
        }

        crate::read_event_discovery::restore_p12_cluster_genotyping_events(
            &mut assembly,
            &apply_bases,
            apply_pad,
            region.start.get(),
            region.end.get(),
            &region.contig,
        );
        crate::read_event_discovery::restore_p12_phase_e_genotyping_events(
            &mut assembly,
            &region.reads,
            region.start.get(),
            region.end.get(),
            &region.contig,
        );
        let genotyping_events = assembly.variation_events.clone();

        let (genotype, genotyped_calls) = if args.run_genotyping
            && (!read_likelihoods.is_empty() || !genotyping_events.is_empty())
        {
            // L4.2 read split (Java parity):
            // `region_for_genotyping.reads`: PairHMM + retainEvidence indices (post-trim/realign)
            // `region.reads`: untrimmed active-region pileup for AD / sparse rescue only
            let gt_region = crate::genotype_site::GenotypeSiteRegion {
                likelihoods: &read_likelihoods,
                likelihood_reads: &region_for_genotyping.reads,
                pileup_reads: &region.reads,
                supplemental_pileup_reads: if args.is_strict_java() {
                    Some(region.reads.as_slice())
                } else {
                    None
                },
                haplotypes: &assembly.haplotypes,
                ref_bytes: &apply_bases,
                pad_start_1based: GenomePosition::new_1based(apply_pad),
                full_reference_bases: assembly.reference_bases(),
                full_reference_pad_1based: GenomePosition::new_1based(
                    assembly.padded_reference_start_1based(),
                ),
                active_start_1based: region.start,
                active_end_1based: region.end,
                contig: &region.contig,
                max_mnp_distance: assembly.max_mnp_distance(),
                stored_events: &genotyping_events,
                graph_events: &graph_events,
            };
            let gt = assign_genotype_likelihoods_for_region(
                gt_region.likelihoods,
                gt_region.likelihood_reads,
                gt_region.pileup_reads,
                gt_region.supplemental_pileup_reads,
                gt_region.haplotypes,
                gt_region.ref_bytes,
                gt_region.pad_start_1based.get(),
                gt_region.full_reference_bases,
                gt_region.full_reference_pad_1based.get(),
                gt_region.active_start_1based.get(),
                gt_region.active_end_1based.get(),
                gt_region.contig,
                gt_region.max_mnp_distance,
                &args.genotyping,
                gt_region.stored_events,
                gt_region.graph_events,
            )?;
            (Some(gt.region_summary), gt.calls)
        } else {
            (None, Vec::new())
        };
        let mut genotyped_calls = genotyped_calls;
        if args.is_strict_java() {
            filter_genotyped_calls_for_strict_java_emit(
                &mut genotyped_calls,
                &region.reads,
                &assembly,
                &args.genotyping,
            )?;
        }
        if !genotyped_calls.is_empty() {
            union_genotyped_calls_into_variation_events(&mut assembly, &genotyped_calls);
        }
        // Observe-only semantic checkpoints (no algorithm effect).
        if crate::semantic_trace::is_enabled() {
            crate::semantic_trace::emit_post_assemble(region, &assembly);
            crate::semantic_trace::emit_read_likelihoods(
                region,
                &read_likelihoods,
                assembly.haplotypes.len(),
            );
            crate::semantic_trace::emit_genotype_likelihoods(region, &genotyped_calls);
        }
        Ok(Some(CallRegionOutcome {
            assembly,
            read_likelihoods,
            genotyping_reads: region_for_genotyping.reads.clone(),
            genotype,
            genotyped_calls,
        }))
    }

    /// Assembly-only hook (legacy parity name).
    /// Inactive `callRegion` → `referenceModelForNoVariation` (no assembly).
    pub fn call_region_inactive_reference(
        region: &AssemblyRegion,
        header: &rust_htslib::bam::HeaderView,
        dictionary: &SequenceDictionary,
        reference_fasta: &Path,
        emit_mode: EmitMode,
        read_filters: &ReadFilterParams,
        ref_confidence_config: &ReferenceConfidenceConfig,
    ) -> GatkResult<InactiveReferenceModelOutcome> {
        let mut ref_cache = ReferenceWindowCache::new(reference_fasta.to_path_buf(), 4);
        reference_model_for_no_variation_region(
            region,
            header,
            ref_confidence_config,
            read_filters,
            &mut ref_cache,
            dictionary,
            emit_mode,
        )
    }

    /// Assembly-only hook (legacy parity name).
    pub fn call_region_assemble(
        region: &AssemblyRegion,
        dictionary: &SequenceDictionary,
        reference_fasta: &Path,
        args: &AssembleReadsArgs,
    ) -> GatkResult<Option<AssemblyResultSet>> {
        let mut owned = region.clone();
        let mut ref_cache = ReferenceWindowCache::new(reference_fasta.to_path_buf(), 4);
        Ok(call_region_assemble(&mut owned, dictionary, &mut ref_cache, args)?.map(|a| a.assembly))
    }

    /// GATK `HaplotypeCallerEngine.shutdown` — native PairHMM teardown hook (J.1.2).
    pub fn shutdown() {
        // Native LIBS PairHMM has no persistent global state in gatk-rs yet.
    }
}

fn reads_overlapping_active_span(
    reads: &[crate::shared_bam::SharedBamRecord],
    active_start_1based: u64,
    active_end_1based: u64,
) -> Vec<crate::shared_bam::SharedBamRecord> {
    reads
        .iter()
        .filter(|r| read_overlaps_variant(r, active_start_1based, active_end_1based, 0))
        .cloned()
        .collect()
}

fn refresh_region_read_likelihoods(
    region: &AssemblyRegion,
    source_reads: &[crate::shared_bam::SharedBamRecord],
    haplotypes: &[Haplotype],
    padded_reference_start_1based: u64,
    config: &HcLikelihoodEngineConfig,
    sw: &crate::alignment::SwParameters,
    apply_normalize: bool,
) -> GatkResult<(
    Vec<RegionReadLikelihood>,
    Vec<crate::shared_bam::SharedBamRecord>,
)> {
    if haplotypes.is_empty() {
        return Ok((Vec::new(), Vec::new()));
    }
    let mut work = region.clone();
    work.reads = reads_overlapping_active_span(source_reads, region.start.get(), region.end.get());
    if work.reads.is_empty() {
        return Ok((Vec::new(), Vec::new()));
    }
    // Refresh uses a different overlapping subset — re-finalize (no assemble buffer).
    let ll = compute_region_read_likelihoods(&work, haplotypes, config, apply_normalize, None)?;
    if ll.is_empty() {
        return Ok((ll, work.reads));
    }
    let (_realigned, best_hap_per_read) = realign_reads_to_best_haplotype(
        work.reads.as_mut_slice(),
        haplotypes,
        &ll,
        padded_reference_start_1based,
        sw,
    )?;
    let ll = crate::read_realignment::change_evidence_to_best_haplotype(ll, &best_hap_per_read);
    Ok((ll, work.reads))
}

/// GATK `--phred-scaled-global-read-mismapping-rate` default 45 → log10 error prob.
const LOG10_GLOBAL_READ_MISMATCHING_RATE: f64 = -4.5;
const EXPECTED_ERROR_RATE_PER_BASE: f64 = 0.02;
const READ_DISQUALIFICATION_LOG10_PER_ERROR: f64 = -4.0;
/// P12 cluster upstream: Java retains long marginal reads (92305716 class, max LL ~−13).
const P12_CLUSTER_UPSTREAM_MARGINAL_READ_KEEP_LOG10: f64 = -13.5;

fn retain_marginal_p12_cluster_upstream_read(
    best: f64,
    qual_len: usize,
    rec: &rust_htslib::bam::Record,
) -> bool {
    qual_len >= 100
        && best >= P12_CLUSTER_UPSTREAM_MARGINAL_READ_KEEP_LOG10
        && read_overlaps_variant(rec, P12_CLUSTER_UPSTREAM_START, P12_CLUSTER_UPSTREAM_END, 0)
}

/// P12 sparse BAM: retain soft-unclipped reads in active span when LL is marginal (92318227).
fn retain_marginal_sparse_softclip_read(
    best: f64,
    qual_len: usize,
    rec: &rust_htslib::bam::Record,
    active_start: u64,
    active_end: u64,
) -> bool {
    use crate::hc_genotyping_engine::{
        java_alignment_read_overlaps_interval, soft_unclipped_read_overlaps_interval,
    };
    let soft_overlaps = soft_unclipped_read_overlaps_interval(rec, active_start, active_end, 2);
    let align_overlaps = java_alignment_read_overlaps_interval(rec, active_start, active_end, 2);
    if !soft_overlaps {
        return false;
    }
    // Java retains marginal soft-clip reads that fail the static threshold (92318227/92318325).
    if !align_overlaps && qual_len >= 20 && best >= P12_CLUSTER_UPSTREAM_MARGINAL_READ_KEEP_LOG10 {
        return true;
    }
    qual_len >= 50 && best >= P12_CLUSTER_UPSTREAM_MARGINAL_READ_KEEP_LOG10 && !align_overlaps
}

/// When static filter drops every read, retain all marginal soft-clip reads in the active span.
fn retain_marginal_sparse_softclip_likelihoods<R: std::borrow::Borrow<rust_htslib::bam::Record>>(
    ll: &[RegionReadLikelihood],
    reads: &[R],
    active_span: Option<(u64, u64)>,
) -> Vec<RegionReadLikelihood> {
    let Some((active_start, active_end)) = active_span else {
        return Vec::new();
    };
    if ll.is_empty() {
        return Vec::new();
    }
    let max_read = ll.iter().map(|e| e.read_index.get()).max().unwrap_or(0);
    let mut keep = std::collections::BTreeSet::new();
    for read_idx in 0..=max_read {
        let Some(rec) = reads.get(read_idx) else {
            continue;
        };
        let rec = rec.borrow();
        let best = ll
            .iter()
            .filter(|e| e.read_index.get() == read_idx)
            .map(|e| e.log10_likelihood)
            .fold(f64::NEG_INFINITY, f64::max);
        if !best.is_finite() {
            continue;
        }
        let qual_len = rec.qual().len().max(1);
        if retain_marginal_sparse_softclip_read(best, qual_len, rec, active_start, active_end) {
            keep.insert(read_idx);
        }
    }
    ll.iter()
        .filter(|e| keep.contains(&e.read_index.get()))
        .cloned()
        .collect()
}

/// Java `ReadLikelihoodCalculationEngine.log10MinTrueLikelihood` (capLikelihoods=true).
fn log10_min_true_likelihood_for_read(read_len: usize) -> f64 {
    let max_errors = (read_len as f64 * EXPECTED_ERROR_RATE_PER_BASE)
        .ceil()
        .min(2.0);
    max_errors * READ_DISQUALIFICATION_LOG10_PER_ERROR
}

/// Hap indices eligible for PairHMM normalize max (drop spillover alts far longer than ref).
fn pairhmm_eligible_haplotype_indices(haplotypes: &[Haplotype]) -> Vec<usize> {
    let ref_len = haplotypes
        .iter()
        .find(|h| h.is_reference)
        .map(|h| h.bases.len())
        .unwrap_or(0);
    if ref_len == 0 {
        return haplotypes.iter().enumerate().map(|(i, _)| i).collect();
    }
    let max_alt_len = ref_len.saturating_add(8);
    haplotypes
        .iter()
        .enumerate()
        .filter(|(_, h)| h.is_reference || h.bases.len() <= max_alt_len)
        .map(|(i, _)| i)
        .collect()
}

/// Java `AlleleLikelihoods.normalizeLikelihoods` per read (symmetric ref allowed).
fn normalize_region_read_likelihoods(
    ll: &mut [RegionReadLikelihood],
    eligible_hap_indices: &[usize],
) {
    if ll.is_empty() {
        return;
    }
    let eligible: std::collections::BTreeSet<usize> =
        eligible_hap_indices.iter().copied().collect();
    let max_read = ll.iter().map(|e| e.read_index.get()).max().unwrap_or(0);
    for read_idx in 0..=max_read {
        let mut best = f64::NEG_INFINITY;
        for entry in ll.iter().filter(|e| e.read_index.get() == read_idx) {
            if !eligible.contains(&entry.haplotype_index.get()) {
                continue;
            }
            if entry.log10_likelihood > best {
                best = entry.log10_likelihood;
            }
        }
        if !best.is_finite() {
            for entry in ll.iter().filter(|e| e.read_index.get() == read_idx) {
                if entry.log10_likelihood > best {
                    best = entry.log10_likelihood;
                }
            }
        }
        if !best.is_finite() {
            continue;
        }
        let floor = best + LOG10_GLOBAL_READ_MISMATCHING_RATE;
        for entry in ll.iter_mut().filter(|e| e.read_index.get() == read_idx) {
            if entry.log10_likelihood < floor {
                entry.log10_likelihood = floor;
            }
        }
    }
}

/// Java `AlleleLikelihoods.filterPoorlyModeledEvidence` (static threshold, dynamic off).
fn filter_poorly_modeled_region_read_likelihoods<
    R: std::borrow::Borrow<rust_htslib::bam::Record>,
>(
    ll: &[RegionReadLikelihood],
    reads: &[R],
    active_span: Option<(u64, u64)>,
) -> Vec<RegionReadLikelihood> {
    if ll.is_empty() {
        return Vec::new();
    }
    let max_read = ll.iter().map(|e| e.read_index.get()).max().unwrap_or(0);
    let mut keep = std::collections::BTreeSet::new();
    for read_idx in 0..=max_read {
        let Some(rec) = reads.get(read_idx) else {
            continue;
        };
        let rec = rec.borrow();
        let best = ll
            .iter()
            .filter(|e| e.read_index.get() == read_idx)
            .map(|e| e.log10_likelihood)
            .fold(f64::NEG_INFINITY, f64::max);
        if !best.is_finite() {
            continue;
        }
        // Java: HMM base-qual array length when present, else read length.
        let qual_len = rec.qual().len().max(1);
        let mut retain = best >= log10_min_true_likelihood_for_read(qual_len)
            || retain_marginal_p12_cluster_upstream_read(best, qual_len, rec);
        if !retain {
            if let Some((active_start, active_end)) = active_span {
                retain = retain_marginal_sparse_softclip_read(
                    best,
                    qual_len,
                    rec,
                    active_start,
                    active_end,
                );
            }
        }
        if retain {
            keep.insert(read_idx);
        }
    }
    ll.iter()
        .filter(|e| keep.contains(&e.read_index.get()))
        .cloned()
        .collect()
}

/// When every read fails the static threshold, Java may still keep P12 upstream marginal evidence.
fn retain_marginal_cluster_upstream_likelihoods<
    R: std::borrow::Borrow<rust_htslib::bam::Record>,
>(
    ll: &[RegionReadLikelihood],
    reads: &[R],
) -> Vec<RegionReadLikelihood> {
    if ll.is_empty() {
        return Vec::new();
    }
    let max_read = ll.iter().map(|e| e.read_index.get()).max().unwrap_or(0);
    let mut keep = std::collections::BTreeSet::new();
    for read_idx in 0..=max_read {
        let Some(rec) = reads.get(read_idx) else {
            continue;
        };
        let rec = rec.borrow();
        let best = ll
            .iter()
            .filter(|e| e.read_index.get() == read_idx)
            .map(|e| e.log10_likelihood)
            .fold(f64::NEG_INFINITY, f64::max);
        if !best.is_finite() {
            continue;
        }
        let qual_len = rec.qual().len().max(1);
        if retain_marginal_p12_cluster_upstream_read(best, qual_len, rec) {
            keep.insert(read_idx);
        }
    }
    ll.iter()
        .filter(|e| keep.contains(&e.read_index.get()))
        .cloned()
        .collect()
}

fn filter_normalized_region_read_likelihoods<R: std::borrow::Borrow<rust_htslib::bam::Record>>(
    ll: &[RegionReadLikelihood],
    reads: &[R],
    active_span: Option<(u64, u64)>,
) -> Vec<RegionReadLikelihood> {
    let filtered = filter_poorly_modeled_region_read_likelihoods(ll, reads, active_span);
    if !filtered.is_empty() {
        return filtered;
    }
    let sparse_softclip = retain_marginal_sparse_softclip_likelihoods(ll, reads, active_span);
    if !sparse_softclip.is_empty() {
        return sparse_softclip;
    }
    let marginal = retain_marginal_cluster_upstream_likelihoods(ll, reads);
    if !marginal.is_empty() {
        return marginal;
    }
    // Java `filterPoorlyModeledEvidence`: do not retain the full matrix when no read passes.
    Vec::new()
}

fn post_process_pairhmm_likelihoods<R: std::borrow::Borrow<rust_htslib::bam::Record>>(
    mut ll: Vec<RegionReadLikelihood>,
    reads: &[R],
    haplotypes: &[Haplotype],
    apply_normalize: bool,
    active_span: Option<(u64, u64)>,
) -> Vec<RegionReadLikelihood> {
    if !apply_normalize {
        // strict_java: normalize + filter after allele filtering (Java order).
        return ll;
    }
    let eligible = pairhmm_eligible_haplotype_indices(haplotypes);
    normalize_region_read_likelihoods(&mut ll, &eligible);
    filter_normalized_region_read_likelihoods(&ll, reads, active_span)
}

/// Score PairHMM from BAM records without `AssemblyRead` / UTF-8 `String` rematerialization.
///
/// # Observable contract
/// Same finalizeRegion evidence and PairHMM inputs as the prior `records_to_assembly_reads` path
/// (BAM seq/qual bytes are ASCII ACGTN — identical to `String::from_utf8_lossy` for valid records).
fn score_pairhmm_from_records<R: std::borrow::Borrow<rust_htslib::bam::Record>>(
    reads: &[R],
    haplotypes: &[Haplotype],
    config: &HcLikelihoodEngineConfig,
) -> GatkResult<Vec<RegionReadLikelihood>> {
    let eligible = pairhmm_eligible_haplotype_indices(haplotypes);
    // L12-A3: zero-copy hap membership for PairHMM (no post-prune `Vec<u8>` rematerialize).
    let hap_refs: Vec<&[u8]> = eligible
        .iter()
        .map(|&hi| haplotypes[hi].bases.as_slice())
        .collect();
    let mut out = Vec::with_capacity(reads.len() * eligible.len());
    for (ri, rec) in reads.iter().enumerate() {
        let rec = rec.borrow();
        // One BAM packed-seq decode; no intermediate UTF-8 String.
        let bases = rec.seq().as_bytes();
        let scores =
            score_read_against_haplotypes(config, &bases, rec.qual(), rec.mapq(), &hap_refs)?;
        for (score_i, &hi) in eligible.iter().enumerate() {
            out.push(RegionReadLikelihood {
                read_index: crate::bio_ids::ReadIndex::new(ri),
                haplotype_index: crate::bio_ids::HaplotypeIndex::new(hi),
                log10_likelihood: scores[score_i],
            });
        }
    }
    Ok(out)
}

fn compute_region_read_likelihoods(
    region: &AssemblyRegion,
    haplotypes: &[Haplotype],
    config: &HcLikelihoodEngineConfig,
    apply_normalize: bool,
    pre_finalized: Option<Vec<rust_htslib::bam::Record>>,
) -> GatkResult<Vec<RegionReadLikelihood>> {
    if haplotypes.is_empty() {
        return Ok(Vec::new());
    }
    // A2: consume assemble finalize buffer when present (clip in place — no second owned copy).
    let finalized = if let Some(mut pre) = pre_finalized.filter(|p| !p.is_empty()) {
        clip_finalized_reads_in_place(&mut pre, region);
        pre
    } else {
        finalize_region_reads_for_assembly(
            &region.reads,
            region,
            true,
            gatk_min_tail_quality_for_assembly(10),
            false,
        )
    };
    let active_span = Some((region.start.get(), region.end.get()));
    // Trim/hard-clip can drop sparse-BAM reads that still overlap the active locus (P12 92305634).
    if finalized.is_empty() && !region.reads.is_empty() {
        let out = score_pairhmm_from_records(region.reads.as_slice(), haplotypes, config)?;
        return Ok(post_process_pairhmm_likelihoods(
            out,
            region.reads.as_slice(),
            haplotypes,
            apply_normalize,
            active_span,
        ));
    }
    let out = score_pairhmm_from_records(&finalized, haplotypes, config)?;
    Ok(post_process_pairhmm_likelihoods(
        out,
        &finalized,
        haplotypes,
        apply_normalize,
        active_span,
    ))
}

#[cfg(test)]
mod pairhmm_post_process_tests {
    use super::*;

    #[test]
    fn filter_drops_poorly_modeled_read_keeps_informative() {
        let mut ll = vec![
            RegionReadLikelihood {
                read_index: crate::bio_ids::ReadIndex::new(0),
                haplotype_index: crate::bio_ids::HaplotypeIndex::new(0),
                log10_likelihood: -50.0,
            },
            RegionReadLikelihood {
                read_index: crate::bio_ids::ReadIndex::new(0),
                haplotype_index: crate::bio_ids::HaplotypeIndex::new(1),
                log10_likelihood: -7.0,
            },
            RegionReadLikelihood {
                read_index: crate::bio_ids::ReadIndex::new(1),
                haplotype_index: crate::bio_ids::HaplotypeIndex::new(0),
                log10_likelihood: -49.0,
            },
            RegionReadLikelihood {
                read_index: crate::bio_ids::ReadIndex::new(1),
                haplotype_index: crate::bio_ids::HaplotypeIndex::new(1),
                log10_likelihood: -44.0,
            },
        ];
        normalize_region_read_likelihoods(&mut ll, &[0, 1]);
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let qual = b"########################################################################################";
        let mut keep = rust_htslib::bam::Record::new();
        keep.set(b"ok", None, seq, qual);
        let mut drop = rust_htslib::bam::Record::new();
        drop.set(b"bad", None, seq, qual);
        let filtered = filter_poorly_modeled_region_read_likelihoods(
            &ll,
            &crate::shared_bam::share_records(vec![keep, drop]),
            None,
        );
        assert!(
            filtered.iter().any(|e| e.read_index.get() == 0),
            "informative read kept"
        );
        assert!(
            !filtered.iter().any(|e| e.read_index.get() == 1),
            "poorly modeled read dropped"
        );
    }

    #[test]
    fn filter_drops_all_reads_when_none_pass_threshold() {
        let mut ll = vec![
            RegionReadLikelihood {
                read_index: crate::bio_ids::ReadIndex::new(0),
                haplotype_index: crate::bio_ids::HaplotypeIndex::new(0),
                log10_likelihood: -16.5526,
            },
            RegionReadLikelihood {
                read_index: crate::bio_ids::ReadIndex::new(0),
                haplotype_index: crate::bio_ids::HaplotypeIndex::new(1),
                log10_likelihood: -13.0042,
            },
        ];
        normalize_region_read_likelihoods(&mut ll, &[0, 1]);
        let mut rec = rust_htslib::bam::Record::new();
        rec.set(
            b"r2",
            None,
            b"ACGTACGTACGTACGTACGT",
            b"####################",
        );
        let reads = [crate::shared_bam::share_record(rec)];
        let filtered = filter_poorly_modeled_region_read_likelihoods(&ll, &reads, None);
        assert!(
            filtered.is_empty(),
            "92305716 read-2 class max LL -13.0 below Java static threshold -8.0"
        );
        let post = post_process_pairhmm_likelihoods(ll.clone(), &reads, &[], true, None);
        assert!(
            post.is_empty(),
            "Java does not retain full matrix when no read passes filterPoorlyModeledEvidence"
        );
    }
}

/// G-6 audit: when trim drops indel CIGARs, untrimmed indel haps are re-attached (ASM-8 / CR-3).
#[cfg(test)]
mod asm_trim_audit_tests {
    use super::*;
    use crate::alignment::SwParameters;
    use crate::assembly_region_iterator::AssemblyRegion;
    use crate::assembly_result_set::AssemblyResultSet;
    use crate::cigar::{Cigar, CigarOperator};
    use crate::genome_loc::{GenomeLoc, GenomePosition};
    use crate::haplotype::Haplotype;
    use crate::read_threading_assembler::{AssemblyResult, AssemblyStatus};

    #[test]
    fn preserve_untrimmed_indel_haplotypes_reattaches_when_trim_loses_indel_cigar() {
        let sw = SwParameters::gatk_haplotype_to_reference();
        let ref_bases = b"ACGTACGTACGTACGT".to_vec();
        let span = GenomeLoc::new(100, 115);
        let mut ref_cigar = Cigar::new();
        ref_cigar.push(ref_bases.len(), CigarOperator::Match);
        // CLONE: needed because haplotype constructor takes owned bases.
        let mut untrimmed_ref = Haplotype::new(ref_bases.clone(), true);
        // CLONE: needed because haplotype owns CIGAR.
        untrimmed_ref.cigar = Some(ref_cigar.clone());
        untrimmed_ref.genome_loc = Some(span);
        let mut untrimmed_alt = Haplotype::new(b"ACGTACGTTACGTACGT".to_vec(), false);
        let mut indel_cigar = Cigar::new();
        indel_cigar.push(8, CigarOperator::Match);
        indel_cigar.push(1, CigarOperator::Insertion);
        indel_cigar.push(8, CigarOperator::Match);
        untrimmed_alt.cigar = Some(indel_cigar);
        untrimmed_alt.alignment_start_hap_wrt_ref = 0;
        let untrimmed = AssemblyResultSet::from_assembly_for_calling(
            &AssemblyResult {
                status: AssemblyStatus::AssembledSomeVariation,
                kmer_size: 10,
                haplotypes: vec![untrimmed_alt, untrimmed_ref],
                event_maps: Vec::new(),
            },
            ref_bases.as_slice(),
            100,
            "2",
            0,
        );
        let trimmed_region = AssemblyRegion {
            contig: "2".into(),
            start: GenomePosition::new_1based(102),
            end: GenomePosition::new_1based(110),
            is_active: true,
            extended_start: GenomePosition::new_1based(100),
            extended_end: GenomePosition::new_1based(115),
            extension: 0,
            reads: Vec::new(),
            read_qnames: Vec::new(),
            reference: crate::reference_context::ReferenceContext::empty(),
            features: crate::feature_context::FeatureContext::empty(),
            pileup_loci: Vec::new(),
        };
        let mut trimmed_ref = Haplotype::new(ref_bases[2..11].to_vec(), true);
        trimmed_ref.cigar = Some({
            let mut c = Cigar::new();
            c.push(9, CigarOperator::Match);
            c
        });
        trimmed_ref.genome_loc = Some(GenomeLoc::new(100, 115));
        let mut assembly = AssemblyResultSet::from_assembly_for_calling(
            &AssemblyResult {
                status: AssemblyStatus::AssembledSomeVariation,
                kmer_size: 10,
                haplotypes: vec![trimmed_ref],
                event_maps: Vec::new(),
            },
            &ref_bases[2..11],
            100,
            "2",
            0,
        );
        preserve_untrimmed_indel_haplotypes(&untrimmed, &mut assembly, &trimmed_region, &sw);
        assert!(
            assembly
                .haplotypes
                .iter()
                .filter(|h| !h.is_reference)
                .any(|h| {
                    h.cigar
                        .as_ref()
                        .is_some_and(|c| c.elements.iter().any(|e| e.operator.is_indel()))
                }),
            "G-6: indel alt hap must survive trim when assembly lost indel CIGAR"
        );
    }
}
