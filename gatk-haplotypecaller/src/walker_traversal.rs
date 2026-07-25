//! GATK4 `AssemblyRegionWalker#traverse` orchestration.
//! Chains B.1 read shards → B.2 [`AssemblyRegionIterator`] → B.3 [`WalkerApplyStats`].

use crate::assembly_region_iterator::{
    AssemblyRegion, AssemblyRegionIterator, AssemblyRegionIteratorConfig,
};
use crate::feature_context::FeatureDataSources;
use crate::read_downsample::{GatkJavaRng, PositionalDownsamplerConfig};
use crate::read_model::ReadFilterParams;
use crate::read_transformer::{apply_shard_read_pipeline, ShardReadPipelineConfig};
use crate::walker::{make_read_shards, ReadShard};
use crate::walker_apply::WalkerApplyStats;
use gatk_common::GatkResult;
use gatk_core::reference::{IntervalSpec, SequenceDictionary};
use std::path::Path;

/// Walker-level configuration (shard padding + iterator + optional downsampling).
/// # Invariants
/// `assembly_region_padding` equals iterator `assemblyRegionExtension` in GATK defaults helper.
/// Shard pipeline and downsample configs stay consistent when using production preset.
/// # Ownership
/// Owns nested iterator, pipeline, and optional feature source handles.
/// # Mutation
/// Immutable for a traversal pass; regions and stats are output separately.
/// # Biological assumptions
/// Padding/extension must cover reads influencing activity and assembly near interval edges.
/// # Java equivalence
/// GATK `AssemblyRegionWalker` + `AssemblyRegionArgumentCollection` + HC read-transform/downsample hooks.
#[derive(Debug, Clone)]
pub struct WalkerTraversalConfig {
    pub assembly_region_padding: u64,
    pub iterator: AssemblyRegionIteratorConfig,
    pub downsample: PositionalDownsamplerConfig,
    /// Pre/post transform + downsampling order for shard record load (B.4).
    pub shard_pipeline: ShardReadPipelineConfig,
    /// GATK `assemblyRegionArgs.forceActive` (applied when building iterator config).
    pub force_active: bool,
    pub feature_sources: Option<FeatureDataSources>,
    /// GATK `shouldTrackPileupsForAssemblyRegions`.
    pub track_pileups: bool,
}

impl WalkerTraversalConfig {
    /// GATK uses the same value for read-shard padding and `assemblyRegionExtension` on the iterator.
    pub fn gatk_haplotype_caller_defaults(padding: u64) -> Self {
        let mut iterator = AssemblyRegionIteratorConfig::gatk_haplotype_caller_defaults();
        iterator.assembly_region_extension = padding.min(u32::MAX as u64) as u32;
        Self {
            assembly_region_padding: padding,
            iterator,
            downsample: PositionalDownsamplerConfig::disabled(),
            shard_pipeline: ShardReadPipelineConfig::parity_l1_goldens(),
            force_active: false,
            feature_sources: None,
            track_pileups: false,
        }
    }

    /// HC production defaults: IUPAC pre-transform, filter, positional DS cap 50 (B.4.1–B.4.3).
    pub fn gatk_haplotype_caller_production(padding: u64) -> Self {
        let mut cfg = Self::gatk_haplotype_caller_defaults(padding);
        let pipeline = ShardReadPipelineConfig::gatk_haplotype_caller_production();
        cfg.downsample = pipeline.downsample;
        cfg.shard_pipeline = pipeline;
        cfg
    }
}

/// Regions emitted for one B.1 [`ReadShard`].
/// # Invariants
/// `regions` preserves iterator emission order for the shard's padded spans.
/// `shard` metadata matches the shard used to load reads and build the iterator.
/// # Ownership
/// Owns shard copy and full region vector including embedded reads.
/// # Mutation
/// Immutable traversal result per shard after drain completes.
/// # Biological assumptions
/// None documented (structural walker output).
/// # Java equivalence
/// GATK per-shard assembly region list from `AssemblyRegionIterator`.
#[derive(Debug, Clone)]
pub struct WalkerShardTraversal {
    pub shard: ReadShard,
    pub regions: Vec<AssemblyRegion>,
}

/// Full walker pass over all shards for the given intervals.
/// # Invariants
/// `apply_stats.total_apply` equals the sum of regions across all shards when built via
/// [`traverse_assembly_region_walker`].
/// Shard order follows reference dictionary among contigs in user intervals.
/// # Ownership
/// Owns shard traversals and aggregate apply stats.
/// # Mutation
/// Immutable result of a complete walker traversal.
/// # Biological assumptions
/// None documented (orchestration container).
/// # Java equivalence
/// GATK `AssemblyRegionWalker#traverse` over all read shards.
#[derive(Debug, Clone)]
pub struct WalkerTraversal {
    pub shards: Vec<WalkerShardTraversal>,
    pub apply_stats: WalkerApplyStats,
}

/// Drain all assembly regions for one shard (shared by parity dumps and traversal).
pub fn drain_assembly_regions_for_shard(
    shard: &ReadShard,
    dictionary: &SequenceDictionary,
    reference_fasta: &Path,
    alignment_path: &Path,
    read_filters: &ReadFilterParams,
    iterator_cfg: &AssemblyRegionIteratorConfig,
    shard_pipeline: &ShardReadPipelineConfig,
    rng: &mut GatkJavaRng,
) -> GatkResult<Vec<AssemblyRegion>> {
    let (header, mut records) =
        crate::assembly_region_iterator::load_records_for_shard_raw(alignment_path, shard)?;
    apply_shard_read_pipeline(
        &mut records,
        Some(&header),
        read_filters,
        shard_pipeline,
        rng,
    )?;
    let mut iter = AssemblyRegionIterator::try_new(
        shard,
        dictionary,
        reference_fasta,
        records,
        header,
        *read_filters,
        iterator_cfg.clone(),
    )?;
    let mut out = Vec::new();
    while let Some(r) = iter.next_region()? {
        out.push(r);
    }
    Ok(out)
}

/// `AssemblyRegionWalker` over `interval_specs`: one iterator pass per dictionary-order shard.
pub fn traverse_assembly_region_walker(
    dictionary: &SequenceDictionary,
    interval_specs: &[IntervalSpec],
    reference_fasta: &Path,
    alignment_path: &Path,
    read_filters: &ReadFilterParams,
    cfg: &WalkerTraversalConfig,
) -> GatkResult<WalkerTraversal> {
    let read_shards = make_read_shards(dictionary, interval_specs, cfg.assembly_region_padding)?;
    let mut iterator_cfg = cfg.iterator.clone();
    iterator_cfg.force_active = cfg.force_active;
    iterator_cfg.feature_sources = cfg.feature_sources.clone();
    iterator_cfg.track_pileups = cfg.track_pileups;
    let mut per_shard = Vec::with_capacity(read_shards.len());
    let mut apply_stats = WalkerApplyStats::default();
    let mut rng = GatkJavaRng::reset_gatk_default();
    for shard in &read_shards {
        let regions = drain_assembly_regions_for_shard(
            shard,
            dictionary,
            reference_fasta,
            alignment_path,
            read_filters,
            &iterator_cfg,
            &cfg.shard_pipeline,
            &mut rng,
        )?;
        // R4-1: accumulate stats without cloning region.reads (was O(regions) BAM clones).
        let shard_stats = WalkerApplyStats::from_regions(&regions);
        apply_stats.total_apply += shard_stats.total_apply;
        apply_stats.inactive_fast_path += shard_stats.inactive_fast_path;
        apply_stats.active_full += shard_stats.active_full;
        per_shard.push(WalkerShardTraversal {
            // CLONE: needed because shard owned by parallel/worker task.
            shard: shard.clone(),
            regions,
        });
    }
    Ok(WalkerTraversal {
        shards: per_shard,
        apply_stats,
    })
}

/// All [`AssemblyRegion`] values in shard order (clones; prefer [`into_assembly_regions`] in production).
pub fn flatten_assembly_regions(walk: &WalkerTraversal) -> Vec<AssemblyRegion> {
    walk.shards
        .iter()
        .flat_map(|s| s.regions.iter().cloned())
        .collect()
}

/// Move regions out of a traversal (no BAM record clones).
pub fn into_assembly_regions(walk: WalkerTraversal) -> Vec<AssemblyRegion> {
    walk.shards.into_iter().flat_map(|s| s.regions).collect()
}

/// Stream assembly regions shard-by-shard, invoking `on_region` as each region is finalized.
///
/// Unlike [`collect_assembly_regions`] / [`into_assembly_regions`], this does **not** retain
/// prior regions (and their owned BAM record clones) for the whole interval — peak memory stays
/// near one shard's `all_records` plus the in-flight region(s) the callback holds.
pub fn for_each_assembly_region<F>(
    dictionary: &SequenceDictionary,
    interval_specs: &[IntervalSpec],
    reference_fasta: &Path,
    alignment_path: &Path,
    read_filters: &ReadFilterParams,
    cfg: &WalkerTraversalConfig,
    mut on_region: F,
) -> GatkResult<()>
where
    F: FnMut(usize, AssemblyRegion) -> GatkResult<()>,
{
    let read_shards = make_read_shards(dictionary, interval_specs, cfg.assembly_region_padding)?;
    let mut iterator_cfg = cfg.iterator.clone();
    iterator_cfg.force_active = cfg.force_active;
    iterator_cfg.feature_sources = cfg.feature_sources.clone();
    iterator_cfg.track_pileups = cfg.track_pileups;
    let mut rng = GatkJavaRng::reset_gatk_default();
    let mut region_index = 0usize;
    for shard in &read_shards {
        let (header, mut records) =
            crate::assembly_region_iterator::load_records_for_shard_raw(alignment_path, shard)?;
        apply_shard_read_pipeline(
            &mut records,
            Some(&header),
            read_filters,
            &cfg.shard_pipeline,
            &mut rng,
        )?;
        let mut iter = AssemblyRegionIterator::try_new(
            shard,
            dictionary,
            reference_fasta,
            records,
            header,
            *read_filters,
            iterator_cfg.clone(),
        )?;
        while let Some(region) = iter.next_region()? {
            on_region(region_index, region)?;
            region_index += 1;
        }
    }
    Ok(())
}

/// Convenience: full walker traversal then flattened region list (B.2 dumps / multi-span fixtures).
pub fn collect_assembly_regions(
    dictionary: &SequenceDictionary,
    interval_specs: &[IntervalSpec],
    reference_fasta: &Path,
    alignment_path: &Path,
    read_filters: &ReadFilterParams,
    cfg: &WalkerTraversalConfig,
) -> GatkResult<Vec<AssemblyRegion>> {
    Ok(into_assembly_regions(traverse_assembly_region_walker(
        dictionary,
        interval_specs,
        reference_fasta,
        alignment_path,
        read_filters,
        cfg,
    )?))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::walker::GATK_DEFAULT_ASSEMBLY_REGION_PADDING;
    use gatk_core::reference::parse_intervals_cli_string;

    #[test]
    fn traversal_matches_direct_drain_on_chr1_fixture() {
        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../parity/fixtures");
        let ref_fa = root.join("reference.fa");
        let bam = root.join("sample.bam");
        let dict = SequenceDictionary::from_fasta_path(&ref_fa).unwrap();
        let specs = parse_intervals_cli_string(&dict, "chr1:5-15").unwrap();
        let filters = ReadFilterParams::default();
        let iter_cfg = AssemblyRegionIteratorConfig::gatk_haplotype_caller_defaults();
        let pipeline = ShardReadPipelineConfig::parity_l1_goldens();
        let shards = make_read_shards(&dict, &specs, GATK_DEFAULT_ASSEMBLY_REGION_PADDING).unwrap();
        assert_eq!(shards.len(), 1);
        let mut rng = GatkJavaRng::reset_gatk_default();
        let direct = drain_assembly_regions_for_shard(
            &shards[0], &dict, &ref_fa, &bam, &filters, &iter_cfg, &pipeline, &mut rng,
        )
        .unwrap();
        let walk = traverse_assembly_region_walker(
            &dict,
            &specs,
            &ref_fa,
            &bam,
            &filters,
            &WalkerTraversalConfig::gatk_haplotype_caller_defaults(
                GATK_DEFAULT_ASSEMBLY_REGION_PADDING,
            ),
        )
        .unwrap();
        assert_eq!(walk.shards.len(), 1);
        assert_eq!(walk.shards[0].regions, direct);
        assert_eq!(walk.apply_stats, WalkerApplyStats::from_regions(&direct));
    }
}
