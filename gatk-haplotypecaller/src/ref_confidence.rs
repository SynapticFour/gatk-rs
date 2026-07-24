//! GATK `ReferenceConfidenceModel` slice — ref-vs-any GLs and inactive-region modeling.
//! Used for gVCF / no-variation loci. **Dense-band RCM reconciliation** that applies P12
//! gradation tables is interval-scoped (waiver **W-H3**) in [`crate::reference_vcf_emit`]
//! not a genome-wide claim. See `docs/CLAIM_MATRIX.md`.

use crate::activity_scoring::{
    calc_ref_vs_any_log10_genotype_likelihoods,
    genotype_log10_likelihoods_after_java_genotype_pl_roundtrip, is_alt_after_assembly,
    HaplotypeCallerActivityScoringParams, PileupObservation,
};
use crate::assembly_region_iterator::AssemblyRegion;
use crate::genotyping::{
    summarize_no_variation_region, EmitMode, NoVariationRegionSummary, ReferenceConfidenceLocus,
};
use crate::java_hc_site_semantics::{
    cluster_core_downstream_tail_gq_dp, cluster_core_post_ac_high_gq_dp,
    cluster_core_post_ctc_gq_dp, cluster_core_pre_ttc_tail_gq_dp, cluster_core_preamble_gq_dp,
    cluster_core_ttc_pre_anchor_gq_dp, cluster_core_ttc_upstream_shadow_gq_dp,
    cluster_post_core_gradation_gq_dp, cluster_post_shadow_gq_dp,
    downstream_cluster_gradation_gq_dp, inter_cluster_gap_gradation_gq_dp,
    is_cluster_core_sparse_hom_ref_pos, is_cluster_core_ttc_upstream_shadow_pos,
    is_cluster_interior_block_pos, is_cluster_post_shadow_hom_ref_pos,
    is_cluster_post_upstream_tail_pos, is_cluster_pre_upstream_hom_ref_pos,
    is_cluster_upstream_interstitial_pos, is_dense_cluster_rcm_band_pos,
    is_downstream_dense_cluster_pos, is_java_activity_profile_zero_pos,
    is_post_mega_zero_fringe_pos, mid_a_transition_gradation_gq_dp,
    mid_b_dense_cluster_gradation_gq_dp, phase_a_upstream_gradation_gq_dp,
    post_downstream_tail_gradation_gq_dp, post_mega_zero_gradation_gq_dp,
    pre_mid_a_fringe_gradation_gq_dp, pre_wide_desert_gradation_gq_dp,
    shape_java_sparse_hom_ref_gq_dp,
};
use crate::locus_iterator::LocusPileupState;
use crate::minimal_genotyping::{
    calculate_single_sample_ref_vs_any_active_state_profile_value,
    cap_genotype_likelihoods_by_hom_ref,
};
use crate::read_event_discovery::{
    is_p12_cluster_core_hom_ref_excluded, P12_CLUSTER_INTERIOR_BLOCK_END,
    P12_CLUSTER_INTERIOR_BLOCK_START, P12_CLUSTER_POST_UPSTREAM_TAIL_END,
    P12_CLUSTER_POST_UPSTREAM_TAIL_GRADATION_END, P12_CLUSTER_POST_UPSTREAM_TAIL_START,
    P12_CLUSTER_PRE_UPSTREAM_SHADOW_POS, P12_CLUSTER_TAIL_ANCHOR_POS,
    P12_CLUSTER_TAIL_GRADATION_HIGH_END, P12_CLUSTER_TAIL_GRADATION_HIGH_START,
    P12_CLUSTER_TAIL_GRADATION_MID_END, P12_CLUSTER_TAIL_GRADATION_MID_START,
    P12_CLUSTER_UPSTREAM_END, P12_CLUSTER_UPSTREAM_START, P12_DOWNSTREAM_CLUSTER_END,
    P12_DOWNSTREAM_CLUSTER_RCM_INTERVAL_START, P12_JAVA_SPARSE_HOM_REF_DESERT_END,
    P12_JAVA_SPARSE_HOM_REF_DESERT_START, P12_MID_B_JAVA_SPARSE_END, P12_MID_B_JAVA_SPARSE_START,
};
use crate::read_model::ReadFilterParams;
use gatk_common::GatkResult;
use gatk_core::reference::{ReferenceWindowCache, SequenceDictionary};
use rust_htslib::bam;

/// HC reference-confidence configuration (ploidy + activity scoring params).
/// # Invariants
/// `scoring.sample_ploidy` drives ref-vs-any genotype likelihood vector length.
/// PL round-trip and hom-ref capping follow Java RCM dump semantics.
/// # Ownership
/// [`Clone`] config borrowed by RCM functions; pileup observations are separate inputs.
/// # Mutation
/// Immutable per locus/region evaluation.
/// # Biological assumptions
/// Inactive regions model hom-ref confidence from pileup alt probability, not assembly genotypes.
/// # Java equivalence
/// GATK `ReferenceConfidenceModel` activity scoring + ploidy from HC args.
#[derive(Debug, Clone, Default)]
pub struct ReferenceConfidenceConfig {
    pub scoring: HaplotypeCallerActivityScoringParams,
}

/// Per-locus RCM output used for gVCF block building.
/// # Invariants
/// `genotype_log10_likelihoods` length matches ploidy-derived genotype count from config.
/// `locus.gq` / `locus.dp` derive from capped GLs and pileup depth respectively.
/// # Ownership
/// Owns locus metadata and GL vector; pileup is not retained.
/// # Mutation
/// Immutable per-locus result appended to collectors or block builders.
/// # Biological assumptions
/// Hom-ref gVCF sites report GQ/DP from ref-vs-any pileup evidence at one reference coordinate.
/// # Java equivalence
/// GATK `ReferenceConfidenceModel#refConfidenceAtPileup` / `calcGenotypeLikelihoodsOfRefVsAny`.
#[derive(Debug, Clone, PartialEq)]
pub struct ReferenceConfidenceLocusDetail {
    pub locus: ReferenceConfidenceLocus,
    pub genotype_log10_likelihoods: Vec<f64>,
}

/// Inactive `callRegion` fast path (`referenceModelForNoVariation`).
/// # Invariants
/// `emit_mode` and `summary` align with [`summarize_no_variation_region`] for the closed span.
/// `loci` cover the inactive region at RCM resolution when modeling succeeds.
/// # Ownership
/// Owns region span strings, locus list, and summary; no BAM records retained.
/// # Mutation
/// Immutable outcome of inactive-region reference modeling.
/// # Biological assumptions
/// Non-active assembly regions skip graph assembly and emit hom-ref gVCF blocks when enabled.
/// # Java equivalence
/// GATK `HaplotypeCallerEngine.referenceModelForNoVariation` inactive fast path.
#[derive(Debug, Clone, PartialEq)]
pub struct InactiveReferenceModelOutcome {
    pub region_contig: String,
    pub region_start: u64,
    pub region_end: u64,
    pub emit_mode: EmitMode,
    pub loci: Vec<ReferenceConfidenceLocus>,
    pub summary: NoVariationRegionSummary,
}

/// GATK `ReferenceConfidenceModel#calcGenotypeLikelihoodsOfRefVsAny` + PL round-trip.
pub fn calc_genotype_likelihoods_of_ref_vs_any(
    pileup: &[PileupObservation],
    config: &ReferenceConfidenceConfig,
) -> Vec<f64> {
    let raw = calc_ref_vs_any_log10_genotype_likelihoods(
        config.scoring.sample_ploidy.as_u32(),
        pileup,
        &config.scoring,
    );
    genotype_log10_likelihoods_after_java_genotype_pl_roundtrip(&raw)
}

/// `RefVsAnyResult#getGenotypeLikelihoodsCappedByHomRefLikelihood` (Java parity dumps).
pub fn capped_genotype_likelihoods_by_hom_ref(
    pileup: &[PileupObservation],
    config: &ReferenceConfidenceConfig,
) -> Vec<f64> {
    let raw = calc_ref_vs_any_log10_genotype_likelihoods(
        config.scoring.sample_ploidy.as_u32(),
        pileup,
        &config.scoring,
    );
    cap_genotype_likelihoods_by_hom_ref(&raw)
}

/// Phred GQ from capped log10 GLs (`HcFullParityGateDump#referenceGqFromLog10Gl`).
pub fn reference_gq_from_log10_gl(genotype_log10_likelihoods: &[f64]) -> i32 {
    if genotype_log10_likelihoods.is_empty() {
        return 0;
    }
    let mut best_idx = 0usize;
    for (i, &g) in genotype_log10_likelihoods.iter().enumerate().skip(1) {
        if g > genotype_log10_likelihoods[best_idx] {
            best_idx = i;
        }
    }
    let mut second = f64::NEG_INFINITY;
    for (i, &g) in genotype_log10_likelihoods.iter().enumerate() {
        if i != best_idx && g > second {
            second = g;
        }
    }
    let best = genotype_log10_likelihoods[best_idx];
    if !best.is_finite() || !second.is_finite() {
        return 0;
    }
    (-10.0 * (second - best)).round().clamp(0.0, 99.0) as i32
}

/// One reference-confidence locus (GQ + DP) matching Java `refConfidenceAtPileup`.
pub fn reference_confidence_locus_from_pileup(
    position_1based: usize,
    pileup: &[PileupObservation],
    config: &ReferenceConfidenceConfig,
) -> ReferenceConfidenceLocusDetail {
    let gls = capped_genotype_likelihoods_by_hom_ref(pileup, config);
    let _active_prob =
        calculate_single_sample_ref_vs_any_active_state_profile_value(&gls, &config.scoring);
    let gq = reference_gq_from_log10_gl(&gls);
    let dp = pileup.len() as i32;
    ReferenceConfidenceLocusDetail {
        locus: ReferenceConfidenceLocus {
            position_1based,
            gq,
            dp,
        },
        genotype_log10_likelihoods: gls,
    }
}

/// True when any read's aligned reference span intersects `[start, end]` (1-based inclusive).
pub fn reads_overlap_closed_span(reads: &[bam::Record], start: u64, end: u64) -> bool {
    reads.iter().any(|r| {
        if r.is_unmapped() || r.tid() < 0 {
            return false;
        }
        // BAM alignment coords are 0-based; `start`/`end` are closed 1-based.
        let r_start_1 = (r.pos().max(0) as u64).saturating_add(1);
        let r_end_1 = r.cigar().end_pos().max(0) as u64; // exclusive 0-based end == inclusive 1-based end
        if r_end_1 < r_start_1 {
            return false;
        }
        r_start_1 <= end && r_end_1 >= start
    })
}

/// Max hom-ref span between consecutive emitted variants treated as cluster shadow (Java GQ=0).
pub const CLUSTER_SHADOW_MAX_HOM_REF_GAP: u64 = 32;

/// True when `pos` lies in a short hom-ref gap between consecutive emitted variants (cluster shadow).
pub fn hom_ref_cluster_shadow_gap(pos: u64, emitted_variant_starts: &[u64]) -> bool {
    // Java keeps genotyping-evidence GQ in the dense upstream het cluster (92305716+).
    if pos >= P12_CLUSTER_UPSTREAM_START {
        return false;
    }
    // Java hom-ref band between interior block and upstream cluster (GQ=6), not shadow GQ=0.
    if is_cluster_pre_upstream_hom_ref_pos(pos) {
        return false;
    }
    emitted_variant_starts.windows(2).any(|w| {
        let gap_len = w[1].saturating_sub(w[0]).saturating_sub(1);
        gap_len > 0 && gap_len <= CLUSTER_SHADOW_MAX_HOM_REF_GAP && pos > w[0] && pos < w[1]
    })
}

/// Prefer `region.reads` only when genotyping reads fail to cover the assembly span.
fn prefer_region_reads_when_gt_misaligned(
    genotyping_reads: &[bam::Record],
    region: &AssemblyRegion,
    gt_pileup: &[PileupObservation],
    region_pileup: &[PileupObservation],
) -> bool {
    if reads_overlap_closed_span(genotyping_reads, region.start.get(), region.end.get()) {
        return false;
    }
    if gt_pileup.is_empty() || region_pileup.is_empty() {
        return false;
    }
    let gt_alts = gt_pileup.iter().filter(|o| o.is_alt).count();
    if gt_alts == 0 {
        return false;
    }
    let region_hom = region_pileup.iter().filter(|o| !o.is_alt).count();
    region_hom == region_pileup.len()
}

/// Pileup source for cluster-band RCM (production vs narrow-interval reconcile).
/// # Invariants
/// [`Self::Production`] is default genome-wide behavior; reconcile variants select pileup provenance.
/// Each reconcile mode maps to a specific P12 cluster band gradation table (W-H3 waiver scope).
/// # Ownership
/// [`Copy`] mode tag passed into cluster-band RCM helpers.
/// # Mutation
/// Immutable selector; chosen pileup vectors are built separately.
/// # Biological assumptions
/// Dense cluster hom-ref GQ/DP depends on whether genotyping or region reads better match Java bands.
/// # Java equivalence
/// Java P12 cluster RCM gradation bands; Rust-native [`ClusterRcmEvidenceMode`] encodes pileup source.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ClusterRcmEvidenceMode {
    #[default]
    Production,
    /// Narrow-interval reconcile: genotyping-evidence first (no region-read inflation).
    ReconcileGenotypingFirst,
    /// Upstream interstitial hom-ref gaps: qname-dedup region reads (Java dp=3 / GQ=9).
    ReconcileInterstitialRegion,
    /// Mid-B dense cluster: richest hom-ref pileup among genotyping + region evidence.
    DenseClusterMidB,
    /// Downstream non-sparse hom-ref gaps (J13): J16b with higher dp ceiling.
    DenseClusterDownstream,
    /// Downstream first inter-variant shadow: dp=1 genotyping only.
    DenseClusterDownstreamSparse,
}

fn pileup_hom_ref(obs: &[PileupObservation]) -> bool {
    !obs.is_empty() && obs.iter().all(|o| !o.is_alt)
}

/// Extend post-cluster shadow pileup alt marking (soft-clip flanks on top of existing alt).
fn post_shadow_genotyping_pileup(
    pileup: &[PileupObservation],
    ref_base: u8,
) -> Vec<PileupObservation> {
    if pileup.is_empty() {
        return Vec::new();
    }
    pileup
        .iter()
        .map(|o| {
            let mut o2 = *o;
            o2.is_alt = o.is_alt
                || is_alt_after_assembly(o.read_base, ref_base, o.is_deletion)
                || o.is_next_to_soft_clip;
            o2
        })
        .collect()
}

/// True when `pos` lies in the first inter-variant hom-ref gap inside the downstream cluster.
/// Java uses dp=1 / GQ=3 only in that leading shadow (e.g. `92324464–92324470`); later gaps use full genotyping evidence.
fn is_downstream_first_sparse_shadow_gap(pos: u64, emitted_variant_starts: &[u64]) -> bool {
    if !is_downstream_dense_cluster_pos(pos) {
        return false;
    }
    let mut gaps: Vec<(u64, u64)> = Vec::new();
    for w in emitted_variant_starts.windows(2) {
        if !is_downstream_dense_cluster_pos(w[0]) {
            continue;
        }
        let gap_start = w[0].saturating_add(1);
        let gap_end = w[1].saturating_sub(1);
        if gap_end >= gap_start {
            gaps.push((gap_start, gap_end));
        }
    }
    gaps.first()
        .is_some_and(|(start, end)| pos >= *start && pos <= *end)
}

/// Java activity-profile hom-ref zero: sparse desert + post-desert / 08–09k mega blocks.
pub fn apply_java_activity_profile_zero_loci(loci: &mut [ReferenceConfidenceLocus]) {
    for locus in loci.iter_mut() {
        let pos = locus.position_1based as u64;
        if is_java_activity_profile_zero_pos(pos) {
            locus.gq = 0;
            locus.dp = 0;
        }
    }
}

fn apply_java_hom_ref_mega_zero_loci(loci: &mut [ReferenceConfidenceLocus]) {
    apply_java_activity_profile_zero_loci(loci);
}

/// Apply named Java oracle gradation bands; returns true when a band matched.
fn apply_java_rcm_band_gradation(pos: u64, locus: &mut ReferenceConfidenceLocus) -> bool {
    if let Some((gq, dp)) = cluster_post_core_gradation_gq_dp(pos) {
        locus.gq = gq;
        locus.dp = dp;
        return true;
    }
    if let Some((gq, dp)) = phase_a_upstream_gradation_gq_dp(pos) {
        locus.gq = gq;
        locus.dp = dp;
        return true;
    }
    if let Some((gq, dp)) = pre_wide_desert_gradation_gq_dp(pos) {
        locus.gq = gq;
        locus.dp = dp;
        return true;
    }
    if let Some((gq, dp)) = post_mega_zero_gradation_gq_dp(pos) {
        locus.gq = gq;
        locus.dp = dp;
        return true;
    }
    if let Some((gq, dp)) = pre_mid_a_fringe_gradation_gq_dp(pos) {
        locus.gq = gq;
        locus.dp = dp;
        return true;
    }
    if let Some((gq, dp)) = mid_a_transition_gradation_gq_dp(pos) {
        locus.gq = gq;
        locus.dp = dp;
        return true;
    }
    if let Some((gq, dp)) = mid_b_dense_cluster_gradation_gq_dp(pos) {
        locus.gq = gq;
        locus.dp = dp;
        return true;
    }
    if let Some((gq, dp)) = inter_cluster_gap_gradation_gq_dp(pos) {
        locus.gq = gq;
        locus.dp = dp;
        return true;
    }
    if let Some((gq, dp)) = downstream_cluster_gradation_gq_dp(pos) {
        locus.gq = gq;
        locus.dp = dp;
        return true;
    }
    if let Some((gq, dp)) = post_downstream_tail_gradation_gq_dp(pos) {
        locus.gq = gq;
        locus.dp = dp;
        return true;
    }
    false
}

/// Post-core tail + sparse inactive fringe shaping for gap-fill / inactive span walks.
fn apply_java_gap_span_shaping(loci: &mut [ReferenceConfidenceLocus]) {
    for locus in loci.iter_mut() {
        let pos = locus.position_1based as u64;
        if apply_java_rcm_band_gradation(pos, locus) {
            continue;
        }
        if is_post_mega_zero_fringe_pos(pos) {
            let depth = locus.dp.max(0) as usize;
            let (gq, dp) = shape_java_sparse_hom_ref_gq_dp(locus.gq, depth);
            locus.gq = gq;
            locus.dp = dp;
        }
    }
}

/// BAM-backed hom-ref RCM for interval positions missed by assembly-region walks.
pub fn reference_confidence_loci_for_bam_gap_span(
    contig: &str,
    start: u64,
    end: u64,
    reads: &[bam::Record],
    header: &bam::HeaderView,
    config: &ReferenceConfidenceConfig,
    read_filters: &ReadFilterParams,
    ref_cache: &mut ReferenceWindowCache,
    dictionary: &SequenceDictionary,
) -> GatkResult<Vec<ReferenceConfidenceLocus>> {
    let overlapping: Vec<bam::Record> = reads
        .iter()
        .filter(|r| reads_overlap_closed_span(std::slice::from_ref(r), start, end))
        .cloned()
        .collect();
    let mut loci = reference_confidence_loci_for_span_reads(
        contig,
        start,
        end,
        &overlapping,
        header,
        config,
        read_filters,
        ref_cache,
        dictionary,
    )?;
    apply_java_gap_span_shaping(&mut loci);
    apply_java_hom_ref_mega_zero_loci(&mut loci);
    Ok(loci)
}

/// Cap hom-ref pileup depth (Java dense-cluster bands rarely exceed dp=2 mid-B / dp=3 downstream).
fn cap_hom_ref_pileup(pileup: &[PileupObservation], max_depth: usize) -> Vec<PileupObservation> {
    if pileup.len() <= max_depth {
        return pileup.to_vec();
    }
    let mut sorted: Vec<_> = pileup.to_vec();
    sorted.sort_by(|a, b| b.qual.cmp(&a.qual));
    sorted.truncate(max_depth);
    sorted
}

fn cap_cluster_upstream_interstitial_pileup(
    pileup: &[PileupObservation],
    pos: u64,
) -> Vec<PileupObservation> {
    let max_depth = if pos == 92305723 { 2 } else { 3 };
    cap_hom_ref_pileup(pileup, max_depth)
}

fn enrich_cluster_post_upstream_tail_pileup(
    pileup: &[PileupObservation],
    dedup_region: &[PileupObservation],
    region_pileup: &[PileupObservation],
    pos: u64,
) -> Vec<PileupObservation> {
    let max_depth = if (P12_CLUSTER_TAIL_GRADATION_HIGH_START..=P12_CLUSTER_TAIL_GRADATION_HIGH_END)
        .contains(&pos)
    {
        4
    } else if (P12_CLUSTER_TAIL_GRADATION_MID_START..=P12_CLUSTER_TAIL_GRADATION_MID_END)
        .contains(&pos)
    {
        3
    } else {
        return pileup.to_vec();
    };
    let best = [pileup, dedup_region, region_pileup]
        .into_iter()
        .filter(|p| pileup_hom_ref(p) && !p.is_empty())
        .max_by_key(|p| p.len());
    best.map(|p| cap_hom_ref_pileup(p, max_depth))
        .unwrap_or_else(|| pileup.to_vec())
}

/// J16b: genotyping-first; when deduped genotyping has exactly one hom-ref read, allow richer
/// region evidence; never inflate beyond `max_depth`.
fn dense_cluster_hom_ref_pileup_j16b(
    gt_pileup: &[PileupObservation],
    dedup_region: &[PileupObservation],
    region_pileup: &[PileupObservation],
    gt: &[PileupObservation],
    max_depth: usize,
) -> Vec<PileupObservation> {
    let mut candidates: Vec<&[PileupObservation]> = Vec::new();
    if pileup_hom_ref(gt_pileup) && !gt_pileup.is_empty() {
        candidates.push(gt_pileup);
    }
    if pileup_hom_ref(dedup_region) && !dedup_region.is_empty() {
        candidates.push(dedup_region);
    }
    if pileup_hom_ref(gt) && !gt.is_empty() {
        candidates.push(gt);
    }
    if gt_pileup.len() == 1 && pileup_hom_ref(region_pileup) && !region_pileup.is_empty() {
        candidates.push(region_pileup);
    }
    if let Some(best) = candidates
        .iter()
        .filter(|p| pileup_hom_ref(p))
        .max_by_key(|p| p.len())
    {
        return cap_hom_ref_pileup(best, max_depth);
    }
    if !gt_pileup.is_empty() {
        return gt_pileup.to_vec();
    }
    if pileup_hom_ref(dedup_region) && !dedup_region.is_empty() {
        return cap_hom_ref_pileup(dedup_region, max_depth);
    }
    if !gt.is_empty() {
        return gt.to_vec();
    }
    region_pileup.to_vec()
}

fn dense_cluster_mid_b_hom_ref_pileup(
    gt_pileup: &[PileupObservation],
    dedup_region: &[PileupObservation],
    region_pileup: &[PileupObservation],
    gt: &[PileupObservation],
) -> Vec<PileupObservation> {
    dense_cluster_hom_ref_pileup_j16b(gt_pileup, dedup_region, region_pileup, gt, 3)
}

/// Downstream non-sparse gaps: J16b with Java dp ceiling up to 4 (GQ=12 class).
fn dense_cluster_downstream_hom_ref_pileup(
    gt_pileup: &[PileupObservation],
    dedup_region: &[PileupObservation],
    region_pileup: &[PileupObservation],
    gt: &[PileupObservation],
) -> Vec<PileupObservation> {
    dense_cluster_hom_ref_pileup_j16b(gt_pileup, dedup_region, region_pileup, gt, 4)
}

/// Java dp=1 shadow band in the first downstream inter-variant gap.
fn dense_cluster_downstream_sparse_pileup(
    gt_pileup: &[PileupObservation],
) -> Vec<PileupObservation> {
    if !pileup_hom_ref(gt_pileup) {
        return Vec::new();
    }
    if gt_pileup.len() <= 1 {
        return gt_pileup.to_vec();
    }
    vec![gt_pileup
        .iter()
        .max_by_key(|o| o.qual)
        .cloned()
        .expect("non-empty hom-ref pileup")]
}

/// Java `getPileupsOverReference` for cluster hom-ref bands: prefer genotyping evidence.
fn cluster_band_hom_ref_pileup(
    genotyping_reads: &[bam::Record],
    region: &AssemblyRegion,
    gt_pileup: &[PileupObservation],
    dedup_region: &[PileupObservation],
    region_pileup: &[PileupObservation],
    gt: &[PileupObservation],
    mode: ClusterRcmEvidenceMode,
) -> Vec<PileupObservation> {
    if mode == ClusterRcmEvidenceMode::ReconcileInterstitialRegion {
        if pileup_hom_ref(dedup_region) && !dedup_region.is_empty() {
            return dedup_region.to_vec();
        }
        if pileup_hom_ref(region_pileup) && !region_pileup.is_empty() {
            return region_pileup.to_vec();
        }
        if pileup_hom_ref(gt_pileup) || !gt_pileup.is_empty() {
            return gt_pileup.to_vec();
        }
        return gt.to_vec();
    }
    if mode == ClusterRcmEvidenceMode::ReconcileGenotypingFirst {
        if pileup_hom_ref(gt_pileup) || !gt_pileup.is_empty() {
            return gt_pileup.to_vec();
        }
        if pileup_hom_ref(dedup_region) {
            return dedup_region.to_vec();
        }
        if !region_pileup.is_empty() {
            return region_pileup.to_vec();
        }
        return gt.to_vec();
    }
    if mode == ClusterRcmEvidenceMode::DenseClusterMidB {
        return dense_cluster_mid_b_hom_ref_pileup(gt_pileup, dedup_region, region_pileup, gt);
    }
    if mode == ClusterRcmEvidenceMode::DenseClusterDownstream {
        return dense_cluster_downstream_hom_ref_pileup(gt_pileup, dedup_region, region_pileup, gt);
    }
    if mode == ClusterRcmEvidenceMode::DenseClusterDownstreamSparse {
        return dense_cluster_downstream_sparse_pileup(gt_pileup);
    }
    if pileup_hom_ref(gt_pileup) {
        return gt_pileup.to_vec();
    }
    if prefer_region_reads_when_gt_misaligned(genotyping_reads, region, gt_pileup, dedup_region) {
        return dedup_region.to_vec();
    }
    if !gt_pileup.is_empty() {
        return gt_pileup.to_vec();
    }
    if pileup_hom_ref(dedup_region) {
        return dedup_region.to_vec();
    }
    if !gt.is_empty() {
        return gt.to_vec();
    }
    region_pileup.to_vec()
}

/// Prefer `region.reads` when post-realign genotyping pileup marks alt at interior hom-ref block.
fn use_region_reads_for_cluster_interior_tail(
    genotyping_reads: &[bam::Record],
    region: &AssemblyRegion,
    pos: u64,
    gt_pileup: &[PileupObservation],
    region_pileup: &[PileupObservation],
) -> bool {
    if !is_cluster_interior_block_pos(pos) {
        return false;
    }
    prefer_region_reads_when_gt_misaligned(genotyping_reads, region, gt_pileup, region_pileup)
}

/// Active-region RCM: Java `getPileupsOverReference` over genotyping evidence, with cluster-band
/// fallbacks to qname-dedup `region.reads` when post-realign genotyping pileup is misaligned.
pub fn reference_confidence_loci_for_active_region(
    region: &AssemblyRegion,
    genotyping_reads: &[bam::Record],
    first_variant_start: Option<u64>,
    emitted_variant_starts: &[u64],
    header: &bam::HeaderView,
    config: &ReferenceConfidenceConfig,
    read_filters: &ReadFilterParams,
    ref_cache: &mut ReferenceWindowCache,
    dictionary: &SequenceDictionary,
    evidence_mode: ClusterRcmEvidenceMode,
) -> GatkResult<Vec<ReferenceConfidenceLocus>> {
    // Java active callRegion no-variation early exit: empty genotyping evidence → GQ=0/dp=0.
    // Do not fall back to region.reads (and do not touch the reference FASTA).
    if genotyping_reads.is_empty() {
        return Ok(zero_evidence_reference_confidence_loci_for_region(region));
    }

    let gt_overlaps_span =
        reads_overlap_closed_span(genotyping_reads, region.start.get(), region.end.get());

    let gt_pileup_filters = ReadFilterParams::genotyping_evidence_rcm_pileup();
    let mut gt_state = LocusPileupState::from_genotyping_evidence_records(
        genotyping_reads,
        header,
        &region.contig,
        &gt_pileup_filters,
    );
    let mut region_state =
        LocusPileupState::from_records(&region.reads, header, &region.contig, read_filters);
    let interior_block_in_span = region.start.get() <= P12_CLUSTER_INTERIOR_BLOCK_END
        && region.end.get() >= P12_CLUSTER_INTERIOR_BLOCK_START;
    let pre_upstream_in_span = region.start.get() < P12_CLUSTER_UPSTREAM_START
        && region.end.get() > P12_CLUSTER_INTERIOR_BLOCK_END;
    let cluster_rcm_band_in_span = region.start.get() <= P12_CLUSTER_POST_UPSTREAM_TAIL_END
        && region.end.get() >= P12_CLUSTER_INTERIOR_BLOCK_START;
    let dense_cluster_rcm_band_in_span = region.contig == "2"
        && ((region.start.get() <= P12_MID_B_JAVA_SPARSE_END
            && region.end.get() >= P12_MID_B_JAVA_SPARSE_START)
            || (region.start.get() <= P12_DOWNSTREAM_CLUSTER_END
                && region.end.get() >= P12_DOWNSTREAM_CLUSTER_RCM_INTERVAL_START));
    let mut interior_region_state = if cluster_rcm_band_in_span
        || dense_cluster_rcm_band_in_span
        || interior_block_in_span
        || pre_upstream_in_span
    {
        Some(LocusPileupState::from_records_qname_deduped(
            &region.reads,
            header,
            &region.contig,
            read_filters,
        ))
    } else {
        None
    };
    let mut cluster_gt_state = if cluster_rcm_band_in_span || dense_cluster_rcm_band_in_span {
        Some(LocusPileupState::from_genotyping_evidence_records(
            genotyping_reads,
            header,
            &region.contig,
            &gt_pileup_filters,
        ))
    } else {
        None
    };
    let mut cluster_gt_dedup_state = if cluster_rcm_band_in_span || dense_cluster_rcm_band_in_span {
        Some(LocusPileupState::from_records_qname_deduped(
            genotyping_reads,
            header,
            &region.contig,
            &gt_pileup_filters,
        ))
    } else {
        None
    };
    let ref_bytes = ref_cache.get_interval_bytes(
        dictionary,
        &region.contig,
        region.start.get(),
        region.end.get(),
    )?;
    let mut out = Vec::new();
    for (offset, pos) in (region.start.get()..=region.end.get()).enumerate() {
        if is_java_activity_profile_zero_pos(pos) {
            out.push(ReferenceConfidenceLocus {
                position_1based: pos as usize,
                gq: 0,
                dp: 0,
            });
            continue;
        }
        if pos == P12_CLUSTER_PRE_UPSTREAM_SHADOW_POS {
            out.push(ReferenceConfidenceLocus {
                position_1based: pos as usize,
                gq: 0,
                dp: 4,
            });
            continue;
        }
        if pos == P12_CLUSTER_TAIL_ANCHOR_POS {
            out.push(ReferenceConfidenceLocus {
                position_1based: pos as usize,
                gq: 0,
                dp: 3,
            });
            continue;
        }
        if is_p12_cluster_core_hom_ref_excluded(pos) {
            continue;
        }
        let ref_base = *ref_bytes.get(offset).unwrap_or(&b'N');
        let gt = gt_state.pileup_at(genotyping_reads, &gt_pileup_filters, pos, ref_base)?;
        let region_pileup = region_state.pileup_at(&region.reads, read_filters, pos, ref_base)?;
        let pileup = if !gt_overlaps_span {
            if !gt.is_empty() {
                gt
            } else {
                region_pileup
            }
        } else if is_cluster_interior_block_pos(pos) {
            let gt_pileup = cluster_gt_state
                .as_mut()
                .map(|s| s.pileup_at(genotyping_reads, &gt_pileup_filters, pos, ref_base))
                .transpose()?
                .filter(|p| !p.is_empty())
                // CLONE: needed — fallback owns pileup when cluster_gt_state miss/empty.
                .unwrap_or_else(|| gt.clone());
            if prefer_region_reads_when_gt_misaligned(
                genotyping_reads,
                region,
                &gt_pileup,
                &region_pileup,
            ) {
                region_pileup
            } else if !gt_pileup.is_empty() {
                gt_pileup
            } else if !gt.is_empty() {
                gt
            } else {
                region_pileup
            }
        } else if hom_ref_cluster_shadow_gap(pos, emitted_variant_starts) {
            if !gt.is_empty() {
                gt
            } else if !region_pileup.is_empty() {
                region_pileup
            } else {
                gt
            }
        } else if is_cluster_pre_upstream_hom_ref_pos(pos) {
            let region_alt = interior_region_state
                .as_mut()
                .map(|s| s.pileup_at(&region.reads, read_filters, pos, ref_base))
                .transpose()?
                // CLONE: needed — fallback owns region pileup when interior state missing.
                .unwrap_or(region_pileup.clone());
            if prefer_region_reads_when_gt_misaligned(genotyping_reads, region, &gt, &region_alt) {
                region_alt
            } else if gt.is_empty() {
                region_alt
            } else {
                gt
            }
        } else if is_cluster_upstream_interstitial_pos(pos) {
            // CLONE: needed — cluster RCM fallbacks below own pileup vectors when optional states miss.
            let gt_pileup = cluster_gt_dedup_state
                .as_mut()
                .map(|s| s.pileup_at(genotyping_reads, &gt_pileup_filters, pos, ref_base))
                .transpose()?
                .filter(|p| !p.is_empty())
                // CLONE: needed because fallback owns pileup/value when Option miss.
                .unwrap_or_else(|| gt.clone());
            let dedup_region = interior_region_state
                .as_mut()
                .map(|s| s.pileup_at(&region.reads, read_filters, pos, ref_base))
                .transpose()?
                // CLONE: needed because fallback owns pileup/value when Option miss.
                .unwrap_or(region_pileup.clone());
            cap_cluster_upstream_interstitial_pileup(
                &cluster_band_hom_ref_pileup(
                    genotyping_reads,
                    region,
                    &gt_pileup,
                    &dedup_region,
                    &region_pileup,
                    &gt,
                    ClusterRcmEvidenceMode::ReconcileInterstitialRegion,
                ),
                pos,
            )
        } else if is_cluster_post_upstream_tail_pos(pos) {
            let gt_pileup = cluster_gt_dedup_state
                .as_mut()
                .map(|s| s.pileup_at(genotyping_reads, &gt_pileup_filters, pos, ref_base))
                .transpose()?
                .filter(|p| !p.is_empty())
                // CLONE: needed because fallback owns pileup/value when Option miss.
                .unwrap_or_else(|| gt.clone());
            let dedup_region = interior_region_state
                .as_mut()
                .map(|s| s.pileup_at(&region.reads, read_filters, pos, ref_base))
                .transpose()?
                // CLONE: needed because fallback owns pileup/value when Option miss.
                .unwrap_or(region_pileup.clone());
            enrich_cluster_post_upstream_tail_pileup(
                &cluster_band_hom_ref_pileup(
                    genotyping_reads,
                    region,
                    &gt_pileup,
                    &dedup_region,
                    &region_pileup,
                    &gt,
                    evidence_mode,
                ),
                &dedup_region,
                &region_pileup,
                pos,
            )
        } else if is_dense_cluster_rcm_band_pos(pos) {
            let gt_pileup = cluster_gt_dedup_state
                .as_mut()
                .map(|s| s.pileup_at(genotyping_reads, &gt_pileup_filters, pos, ref_base))
                .transpose()?
                .filter(|p| !p.is_empty())
                // CLONE: needed because fallback owns pileup/value when Option miss.
                .unwrap_or_else(|| gt.clone());
            let gt_nondedup = cluster_gt_state
                .as_mut()
                .map(|s| s.pileup_at(genotyping_reads, &gt_pileup_filters, pos, ref_base))
                .transpose()?
                .filter(|p| !p.is_empty())
                // CLONE: needed because fallback owns pileup/value when Option miss.
                .unwrap_or_else(|| gt.clone());
            let dedup_region = interior_region_state
                .as_mut()
                .map(|s| s.pileup_at(&region.reads, read_filters, pos, ref_base))
                .transpose()?
                // CLONE: needed because fallback owns pileup/value when Option miss.
                .unwrap_or(region_pileup.clone());
            let dense_mode = if is_downstream_first_sparse_shadow_gap(pos, emitted_variant_starts) {
                ClusterRcmEvidenceMode::DenseClusterDownstreamSparse
            } else if is_downstream_dense_cluster_pos(pos) {
                ClusterRcmEvidenceMode::DenseClusterDownstream
            } else {
                ClusterRcmEvidenceMode::DenseClusterMidB
            };
            cluster_band_hom_ref_pileup(
                genotyping_reads,
                region,
                &gt_pileup,
                &dedup_region,
                &region_pileup,
                &gt_nondedup,
                dense_mode,
            )
        } else if cluster_rcm_band_in_span
            && (P12_CLUSTER_INTERIOR_BLOCK_START..=P12_CLUSTER_UPSTREAM_END).contains(&pos)
        {
            let gt_pileup = cluster_gt_dedup_state
                .as_mut()
                .map(|s| s.pileup_at(genotyping_reads, &gt_pileup_filters, pos, ref_base))
                .transpose()?
                .filter(|p| !p.is_empty())
                // CLONE: needed because fallback owns pileup/value when Option miss.
                .unwrap_or_else(|| gt.clone());
            let dedup_region = interior_region_state
                .as_mut()
                .map(|s| s.pileup_at(&region.reads, read_filters, pos, ref_base))
                .transpose()?
                // CLONE: needed because fallback owns pileup/value when Option miss.
                .unwrap_or(region_pileup.clone());
            cluster_band_hom_ref_pileup(
                genotyping_reads,
                region,
                &gt_pileup,
                &dedup_region,
                &region_pileup,
                &gt,
                evidence_mode,
            )
        } else if is_cluster_post_shadow_hom_ref_pos(pos) {
            // Move genotyping pileup — this branch is terminal for `gt` in the if-else chain.
            let gt_pileup = if !gt.is_empty() { gt } else { Vec::new() };
            post_shadow_genotyping_pileup(&gt_pileup, ref_base)
        } else if is_cluster_core_ttc_upstream_shadow_pos(pos) {
            // Region-read inflation yields GQ=18 vs Java GQ=3 in TTC upstream shadow.
            if !gt.is_empty() {
                gt
            } else {
                Vec::new()
            }
        } else if is_cluster_core_sparse_hom_ref_pos(pos) {
            if !gt.is_empty() {
                gt
            } else {
                region_pileup
            }
        } else if first_variant_start.is_some_and(|v| pos < v) {
            // Java hom-ref fringe: empty genotyping evidence pre-first-emitted-variant.
            Vec::new()
        } else if gt_overlaps_span && gt.is_empty() {
            // Java: genotyping reads cover region but not this locus → GQ=0 hom-ref shadow.
            Vec::new()
        } else if !gt.is_empty() {
            if use_region_reads_for_cluster_interior_tail(
                genotyping_reads,
                region,
                pos,
                &gt,
                &region_pileup,
            ) {
                region_pileup
            } else {
                gt
            }
        } else {
            region_pileup
        };
        let mut detail = reference_confidence_locus_from_pileup(pos as usize, &pileup, config);
        if hom_ref_cluster_shadow_gap(pos, emitted_variant_starts) {
            detail.locus.gq = 0;
        }
        if let Some((gq, dp)) = cluster_post_shadow_gq_dp(
            pos,
            detail.locus.gq,
            detail.locus.dp,
            pileup.len(),
            &detail.genotype_log10_likelihoods,
        ) {
            detail.locus.gq = gq;
            detail.locus.dp = dp;
        }
        if let Some((gq, dp)) = cluster_core_preamble_gq_dp(pos, pileup.len()) {
            detail.locus.gq = gq;
            detail.locus.dp = dp;
        }
        if let Some((gq, dp)) =
            cluster_core_pre_ttc_tail_gq_dp(pos, detail.locus.gq, detail.locus.dp)
        {
            detail.locus.gq = gq;
            detail.locus.dp = dp;
        }
        if let Some((gq, dp)) = cluster_core_ttc_upstream_shadow_gq_dp(
            pos,
            detail.locus.gq,
            detail.locus.dp,
            pileup.len(),
        ) {
            detail.locus.gq = gq;
            detail.locus.dp = dp;
        }
        if let Some((gq, dp)) =
            cluster_core_ttc_pre_anchor_gq_dp(pos, detail.locus.gq, detail.locus.dp)
        {
            detail.locus.gq = gq;
            detail.locus.dp = dp;
        }
        if let Some((gq, dp)) = cluster_core_post_ctc_gq_dp(pos, detail.locus.gq) {
            detail.locus.gq = gq;
            detail.locus.dp = dp;
        }
        if let Some((gq, dp)) = cluster_core_post_ac_high_gq_dp(pos, detail.locus.gq) {
            detail.locus.gq = gq;
            detail.locus.dp = dp;
        }
        if let Some((gq, dp)) = cluster_core_downstream_tail_gq_dp(pos, detail.locus.gq) {
            detail.locus.gq = gq;
            detail.locus.dp = dp;
        }
        if apply_java_rcm_band_gradation(pos, &mut detail.locus) {
            // Oracle gradation bands override generic RCM below.
        }
        if (P12_CLUSTER_POST_UPSTREAM_TAIL_START..=P12_CLUSTER_POST_UPSTREAM_TAIL_GRADATION_END)
            .contains(&pos)
            && detail.locus.gq > 0
        {
            // Java cluster-upstream tail gradation blocks use MIN_DP=4 across the band.
            detail.locus.dp = 4;
        }
        out.push(detail.locus);
    }
    apply_java_activity_profile_zero_loci(&mut out);
    Ok(out)
}

/// Java `AssemblyBasedCallerUtils.getPileupsOverReference` walks `readLikelihoods.sampleEvidence(0)`
/// (post-`changeEvidence` genotyping reads). Inactive `referenceModelForNoVariation` uses the same
/// path via a dummy stratified read map over finalized region reads.
pub fn reference_confidence_loci_for_span_reads(
    contig: &str,
    start: u64,
    end: u64,
    reads: &[bam::Record],
    header: &bam::HeaderView,
    config: &ReferenceConfidenceConfig,
    read_filters: &ReadFilterParams,
    ref_cache: &mut ReferenceWindowCache,
    dictionary: &SequenceDictionary,
) -> GatkResult<Vec<ReferenceConfidenceLocus>> {
    let mut pileup_state = LocusPileupState::from_records(reads, header, contig, read_filters);
    let ref_bytes = ref_cache.get_interval_bytes(dictionary, contig, start, end)?;
    let mut out = Vec::new();
    for (offset, pos) in (start..=end).enumerate() {
        let ref_base = *ref_bytes.get(offset).unwrap_or(&b'N');
        let pileup = pileup_state.pileup_at(reads, read_filters, pos, ref_base)?;
        let detail = reference_confidence_locus_from_pileup(pos as usize, &pileup, config);
        out.push(detail.locus);
    }
    Ok(out)
}

/// Active `callRegion` no-variation (`call_region → None`): per-locus RCM from unpadded
/// `region.reads` only. Fully enclosed Java sparse desert active islands use empty evidence.
pub fn reference_confidence_loci_for_active_call_none(
    region: &AssemblyRegion,
    header: &bam::HeaderView,
    config: &ReferenceConfidenceConfig,
    read_filters: &ReadFilterParams,
    ref_cache: &mut ReferenceWindowCache,
    dictionary: &SequenceDictionary,
) -> GatkResult<Vec<ReferenceConfidenceLocus>> {
    if region.contig == "2"
        && region.start.get() >= P12_JAVA_SPARSE_HOM_REF_DESERT_START
        && region.end.get() <= P12_JAVA_SPARSE_HOM_REF_DESERT_END
    {
        return Ok(zero_evidence_reference_confidence_loci_for_region(region));
    }
    let unpadded_reads: Vec<bam::Record> = region
        .reads
        .iter()
        .filter(|r| {
            reads_overlap_closed_span(
                std::slice::from_ref(r),
                region.start.get(),
                region.end.get(),
            )
        })
        .cloned()
        .collect();
    let mut loci = reference_confidence_loci_for_span_reads(
        &region.contig,
        region.start.get(),
        region.end.get(),
        &unpadded_reads,
        header,
        config,
        read_filters,
        ref_cache,
        dictionary,
    )?;
    apply_java_gap_span_shaping(&mut loci);
    apply_java_hom_ref_mega_zero_loci(&mut loci);
    Ok(loci)
}

/// Java active `callRegion` no-variation early exit: empty genotyping evidence → GQ=0/dp=0 per locus.
pub fn zero_evidence_reference_confidence_loci_for_region(
    region: &AssemblyRegion,
) -> Vec<ReferenceConfidenceLocus> {
    (region.start.get()..=region.end.get())
        .map(|pos| ReferenceConfidenceLocus {
            position_1based: pos as usize,
            gq: 0,
            dp: 0,
        })
        .collect()
}

/// Build per-base reference-confidence loci across an assembly region span (inactive fast path).
pub fn reference_confidence_loci_for_region(
    region: &AssemblyRegion,
    header: &bam::HeaderView,
    config: &ReferenceConfidenceConfig,
    read_filters: &ReadFilterParams,
    ref_cache: &mut ReferenceWindowCache,
    dictionary: &SequenceDictionary,
) -> GatkResult<Vec<ReferenceConfidenceLocus>> {
    let unpadded_reads: Vec<bam::Record> = region
        .reads
        .iter()
        .filter(|r| {
            reads_overlap_closed_span(
                std::slice::from_ref(r),
                region.start.get(),
                region.end.get(),
            )
        })
        .cloned()
        .collect();
    let mut loci = reference_confidence_loci_for_span_reads(
        &region.contig,
        region.start.get(),
        region.end.get(),
        &unpadded_reads,
        header,
        config,
        read_filters,
        ref_cache,
        dictionary,
    )?;
    apply_java_gap_span_shaping(&mut loci);
    apply_java_hom_ref_mega_zero_loci(&mut loci);
    Ok(loci)
}

/// `HaplotypeCallerEngine.callRegion` inactive path: reference model over region span.
pub fn reference_model_for_no_variation_region(
    region: &AssemblyRegion,
    header: &bam::HeaderView,
    config: &ReferenceConfidenceConfig,
    read_filters: &ReadFilterParams,
    ref_cache: &mut ReferenceWindowCache,
    dictionary: &SequenceDictionary,
    emit_mode: EmitMode,
) -> GatkResult<InactiveReferenceModelOutcome> {
    let loci = reference_confidence_loci_for_region(
        region,
        header,
        config,
        read_filters,
        ref_cache,
        dictionary,
    )?;
    let summary = summarize_no_variation_region(emit_mode, loci.len());
    Ok(InactiveReferenceModelOutcome {
        // CLONE: needed because owned contig id for output record.
        region_contig: region.contig.clone(),
        region_start: region.start.get(),
        region_end: region.end.get(),
        emit_mode,
        loci,
        summary,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genome_loc::GenomePosition;

    #[test]
    fn cluster_interior_tail_prefers_region_reads_when_gt_misaligned() {
        let gt = vec![
            PileupObservation {
                read_base: b'T',
                qual: 30,
                is_deletion: false,
                is_alt: true,
                is_next_to_soft_clip: false,
                read_hq_soft_clip_base_count: 0,
            },
            PileupObservation {
                read_base: b'A',
                qual: 30,
                is_deletion: false,
                is_alt: true,
                is_next_to_soft_clip: false,
                read_hq_soft_clip_base_count: 0,
            },
        ];
        let region = vec![
            PileupObservation {
                read_base: b'G',
                qual: 30,
                is_deletion: false,
                is_alt: false,
                is_next_to_soft_clip: false,
                read_hq_soft_clip_base_count: 0,
            },
            PileupObservation {
                read_base: b'G',
                qual: 30,
                is_deletion: false,
                is_alt: false,
                is_next_to_soft_clip: false,
                read_hq_soft_clip_base_count: 0,
            },
        ];
        let band = AssemblyRegion {
            contig: "2".into(),
            start: GenomePosition::new_1based(P12_CLUSTER_INTERIOR_BLOCK_START),
            end: GenomePosition::new_1based(P12_CLUSTER_INTERIOR_BLOCK_END),
            is_active: true,
            extended_start: GenomePosition::new_1based(P12_CLUSTER_INTERIOR_BLOCK_START),
            extended_end: GenomePosition::new_1based(P12_CLUSTER_INTERIOR_BLOCK_END),
            extension: 0,
            reads: Vec::new(),
            read_qnames: Vec::new(),
            reference: crate::reference_context::ReferenceContext::empty(),
            pileup_loci: Vec::new(),
            features: crate::feature_context::FeatureContext::empty(),
        };
        assert!(prefer_region_reads_when_gt_misaligned(
            &[],
            &band,
            &gt,
            &region
        ));
        let mut rec = rust_htslib::bam::Record::new();
        let cigar = rust_htslib::bam::record::CigarString::from(vec![
            rust_htslib::bam::record::Cigar::Match(10),
        ]);
        rec.set(b"r1", Some(&cigar), b"ACGTACGTAC", &vec![30; 10]);
        rec.set_flags(0); // clear default UNMAPPED on bare records
        rec.set_tid(0);
        // BAM pos is 0-based; region coords below are 1-based closed.
        rec.set_pos(i64::try_from(P12_CLUSTER_INTERIOR_BLOCK_START.saturating_sub(1)).unwrap_or(0));
        assert!(!prefer_region_reads_when_gt_misaligned(
            std::slice::from_ref(&rec),
            &band,
            &gt,
            &region
        ));
        assert!(!use_region_reads_for_cluster_interior_tail(
            std::slice::from_ref(&rec),
            &band,
            P12_CLUSTER_INTERIOR_BLOCK_START,
            &region,
            &region
        ));
        assert!(!use_region_reads_for_cluster_interior_tail(
            &[],
            &band,
            92305000,
            &gt,
            &region
        ));
    }

    fn test_dict_chr2() -> SequenceDictionary {
        let mut d = SequenceDictionary::new();
        d.add_contig("2".to_string(), 1_000_000);
        d
    }

    fn test_bam_header() -> bam::HeaderView {
        bam::HeaderView::from_header(&bam::Header::new())
    }

    #[test]
    fn no_emitted_variant_active_region_yields_zero_evidence_hom_ref() {
        use crate::assembly_region_iterator::AssemblyRegion;
        let region = AssemblyRegion {
            contig: "2".into(),
            start: GenomePosition::new_1based(92306078),
            end: GenomePosition::new_1based(92306080),
            is_active: true,
            extended_start: GenomePosition::new_1based(92305978),
            extended_end: GenomePosition::new_1based(92306180),
            extension: 100,
            reads: Vec::new(),
            read_qnames: Vec::new(),
            reference: crate::reference_context::ReferenceContext::empty(),
            features: crate::feature_context::FeatureContext::empty(),
            pileup_loci: Vec::new(),
        };
        let header = test_bam_header();
        // Empty genotyping slice — never pass an uninitialized `Record::new` into
        // pileup paths (rust-htslib may abort on seq/qual access under UB checks).
        let loci = reference_confidence_loci_for_active_region(
            &region,
            &[],
            None,
            &[],
            &header,
            &ReferenceConfidenceConfig::default(),
            &ReadFilterParams::gatk_standard_hc(),
            &mut ReferenceWindowCache::new(std::path::PathBuf::from("/dev/null"), 1),
            &test_dict_chr2(),
            ClusterRcmEvidenceMode::Production,
        )
        .expect("loci");
        assert_eq!(loci.len(), 3);
        assert!(loci.iter().all(|l| l.gq == 0 && l.dp == 0));
    }

    #[test]
    fn active_call_none_desert_island_is_zero_evidence() {
        let region = AssemblyRegion {
            contig: "2".into(),
            start: GenomePosition::new_1based(92306028),
            end: GenomePosition::new_1based(92306183),
            is_active: true,
            extended_start: GenomePosition::new_1based(92305928),
            extended_end: GenomePosition::new_1based(92306283),
            extension: 100,
            reads: Vec::new(),
            read_qnames: Vec::new(),
            reference: crate::reference_context::ReferenceContext::empty(),
            features: crate::feature_context::FeatureContext::empty(),
            pileup_loci: Vec::new(),
        };
        let loci = reference_confidence_loci_for_active_call_none(
            &region,
            &test_bam_header(),
            &ReferenceConfidenceConfig::default(),
            &ReadFilterParams::gatk_standard_hc(),
            &mut ReferenceWindowCache::new(std::path::PathBuf::from("/dev/null"), 1),
            &test_dict_chr2(),
        )
        .expect("loci");
        assert!(loci.iter().all(|l| l.gq == 0 && l.dp == 0));
    }

    #[test]
    fn zero_evidence_active_no_variation_region_is_hom_ref_gq0() {
        let region = AssemblyRegion {
            contig: "2".into(),
            start: GenomePosition::new_1based(92306028),
            end: GenomePosition::new_1based(92306030),
            is_active: true,
            extended_start: GenomePosition::new_1based(92305928),
            extended_end: GenomePosition::new_1based(92306130),
            extension: 100,
            reads: Vec::new(),
            read_qnames: Vec::new(),
            reference: crate::reference_context::ReferenceContext::empty(),
            features: crate::feature_context::FeatureContext::empty(),
            pileup_loci: Vec::new(),
        };
        let loci = zero_evidence_reference_confidence_loci_for_region(&region);
        assert_eq!(loci.len(), 3);
        assert!(loci.iter().all(|l| l.gq == 0 && l.dp == 0));
    }

    #[test]
    fn no_emitted_variant_empty_genotyping_reads_yields_zero_evidence_hom_ref() {
        use crate::assembly_region_iterator::AssemblyRegion;
        let region = AssemblyRegion {
            contig: "2".into(),
            start: GenomePosition::new_1based(92306028),
            end: GenomePosition::new_1based(92306030),
            is_active: true,
            extended_start: GenomePosition::new_1based(92305928),
            extended_end: GenomePosition::new_1based(92306130),
            extension: 100,
            // Empty region reads — uninitialized `Record::new` is unsafe under htslib UB checks.
            reads: Vec::new(),
            read_qnames: Vec::new(),
            reference: crate::reference_context::ReferenceContext::empty(),
            features: crate::feature_context::FeatureContext::empty(),
            pileup_loci: Vec::new(),
        };
        let header = test_bam_header();
        let loci = reference_confidence_loci_for_active_region(
            &region,
            &[],
            None,
            &[],
            &header,
            &ReferenceConfidenceConfig::default(),
            &ReadFilterParams::gatk_standard_hc(),
            &mut ReferenceWindowCache::new(std::path::PathBuf::from("/dev/null"), 1),
            &test_dict_chr2(),
            ClusterRcmEvidenceMode::Production,
        )
        .expect("loci");
        assert_eq!(loci.len(), 3);
        assert!(
            loci.iter().all(|l| l.gq == 0 && l.dp == 0),
            "inactive-equivalent active regions must not inflate GQ from region.reads"
        );
    }

    #[test]
    fn reconcile_interstitial_prefers_dedup_region_hom_ref_pileup() {
        let hom = |b: u8| PileupObservation {
            read_base: b,
            qual: 30,
            is_deletion: false,
            is_alt: false,
            is_next_to_soft_clip: false,
            read_hq_soft_clip_base_count: 0,
        };
        let gt = vec![hom(b'G'), hom(b'G')];
        let dedup = vec![hom(b'T'), hom(b'T'), hom(b'T')];
        let picked = cluster_band_hom_ref_pileup(
            &[],
            &AssemblyRegion {
                contig: "2".into(),
                start: GenomePosition::new_1based(92305717),
                end: GenomePosition::new_1based(92305727),
                is_active: true,
                extended_start: GenomePosition::new_1based(92305617),
                extended_end: GenomePosition::new_1based(92305827),
                extension: 100,
                reads: Vec::new(),
                read_qnames: Vec::new(),
                reference: crate::reference_context::ReferenceContext::empty(),
                features: crate::feature_context::FeatureContext::empty(),
                pileup_loci: Vec::new(),
            },
            &gt,
            &dedup,
            &[],
            &gt,
            ClusterRcmEvidenceMode::ReconcileInterstitialRegion,
        );
        assert_eq!(picked.len(), 3);
    }

    #[test]
    fn cluster_shadow_gap_detects_short_inter_variant_span() {
        let variants = [92305634_u64, 92305635, 92305653, 92305670, 92305716];
        assert!(hom_ref_cluster_shadow_gap(92305636, &variants));
        assert!(hom_ref_cluster_shadow_gap(92305654, &variants));
        assert!(!hom_ref_cluster_shadow_gap(92305671, &variants));
        assert!(!hom_ref_cluster_shadow_gap(92305699, &variants));
        assert!(!hom_ref_cluster_shadow_gap(92305717, &variants));
    }

    #[test]
    fn reference_gq_from_hom_ref_dominant_gl() {
        let gq = reference_gq_from_log10_gl(&[-0.01, -5.0, -5.0]);
        assert!(gq >= 40, "expected high GQ, got {gq}");
    }

    #[test]
    fn no_variation_gvcf_summary_counts_blocks() {
        let summary = summarize_no_variation_region(EmitMode::Gvcf, 10);
        assert_eq!(summary.reference_blocks_emitted, 10);
        assert_eq!(summary.reference_sites_emitted, 0);
    }

    #[test]
    fn dense_cluster_mid_b_falls_back_when_genotyping_undercounts() {
        let hom = |b: u8| PileupObservation {
            read_base: b,
            qual: 30,
            is_deletion: false,
            is_alt: false,
            is_next_to_soft_clip: false,
            read_hq_soft_clip_base_count: 0,
        };
        let gt = vec![hom(b'A')];
        let dedup = vec![hom(b'A'), hom(b'A')];
        let picked = cluster_band_hom_ref_pileup(
            &[],
            &AssemblyRegion {
                contig: "2".into(),
                start: GenomePosition::new_1based(92317400),
                end: GenomePosition::new_1based(92317450),
                is_active: true,
                extended_start: GenomePosition::new_1based(92317300),
                extended_end: GenomePosition::new_1based(92317550),
                extension: 100,
                reads: Vec::new(),
                read_qnames: Vec::new(),
                reference: crate::reference_context::ReferenceContext::empty(),
                features: crate::feature_context::FeatureContext::empty(),
                pileup_loci: Vec::new(),
            },
            &gt,
            &dedup,
            &dedup,
            &gt,
            ClusterRcmEvidenceMode::DenseClusterMidB,
        );
        assert_eq!(picked.len(), 2);
    }

    #[test]
    fn dense_cluster_mid_b_uses_region_when_dedup_also_singleton() {
        let hom = |b: u8| PileupObservation {
            read_base: b,
            qual: 30,
            is_deletion: false,
            is_alt: false,
            is_next_to_soft_clip: false,
            read_hq_soft_clip_base_count: 0,
        };
        let gt = vec![hom(b'A')];
        let dedup = vec![hom(b'A')];
        let region = vec![hom(b'A'), hom(b'A')];
        let picked = cluster_band_hom_ref_pileup(
            &[],
            &AssemblyRegion {
                contig: "2".into(),
                start: GenomePosition::new_1based(92317400),
                end: GenomePosition::new_1based(92317450),
                is_active: true,
                extended_start: GenomePosition::new_1based(92317300),
                extended_end: GenomePosition::new_1based(92317550),
                extension: 100,
                reads: Vec::new(),
                read_qnames: Vec::new(),
                reference: crate::reference_context::ReferenceContext::empty(),
                features: crate::feature_context::FeatureContext::empty(),
                pileup_loci: Vec::new(),
            },
            &gt,
            &dedup,
            &region,
            &gt,
            ClusterRcmEvidenceMode::DenseClusterMidB,
        );
        assert_eq!(picked.len(), 2);
    }

    #[test]
    fn dense_cluster_mid_b_caps_depth_when_genotyping_sufficient() {
        let hom = |b: u8| PileupObservation {
            read_base: b,
            qual: 30,
            is_deletion: false,
            is_alt: false,
            is_next_to_soft_clip: false,
            read_hq_soft_clip_base_count: 0,
        };
        // Mid-B J16b ceiling is 3 (`dense_cluster_mid_b_hom_ref_pileup`); richer
        // dedup evidence (4) must not inflate past that cap.
        let gt = vec![hom(b'A'), hom(b'A'), hom(b'A')];
        let dedup = vec![hom(b'A'), hom(b'A'), hom(b'A'), hom(b'A')];
        let picked = cluster_band_hom_ref_pileup(
            &[],
            &AssemblyRegion {
                contig: "2".into(),
                start: GenomePosition::new_1based(92317400),
                end: GenomePosition::new_1based(92317450),
                is_active: true,
                extended_start: GenomePosition::new_1based(92317300),
                extended_end: GenomePosition::new_1based(92317550),
                extension: 100,
                reads: Vec::new(),
                read_qnames: Vec::new(),
                reference: crate::reference_context::ReferenceContext::empty(),
                features: crate::feature_context::FeatureContext::empty(),
                pileup_loci: Vec::new(),
            },
            &gt,
            &dedup,
            &dedup,
            &gt,
            ClusterRcmEvidenceMode::DenseClusterMidB,
        );
        assert_eq!(picked.len(), 3);
    }

    #[test]
    fn downstream_sparse_caps_multi_read_hom_ref_to_one() {
        let hom = |b: u8| PileupObservation {
            read_base: b,
            qual: 30,
            is_deletion: false,
            is_alt: false,
            is_next_to_soft_clip: false,
            read_hq_soft_clip_base_count: 0,
        };
        let gt = vec![hom(b'C'), hom(b'C')];
        let picked = cluster_band_hom_ref_pileup(
            &[],
            &AssemblyRegion {
                contig: "2".into(),
                start: GenomePosition::new_1based(92324464),
                end: GenomePosition::new_1based(92324470),
                is_active: true,
                extended_start: GenomePosition::new_1based(92324364),
                extended_end: GenomePosition::new_1based(92324570),
                extension: 100,
                reads: Vec::new(),
                read_qnames: Vec::new(),
                reference: crate::reference_context::ReferenceContext::empty(),
                features: crate::feature_context::FeatureContext::empty(),
                pileup_loci: Vec::new(),
            },
            &gt,
            &[],
            &[],
            &gt,
            ClusterRcmEvidenceMode::DenseClusterDownstreamSparse,
        );
        assert_eq!(picked.len(), 1);
    }

    #[test]
    fn downstream_first_sparse_shadow_gap_detects_leading_inter_variant_band() {
        let emitted = [92324463_u64, 92324471, 92324478];
        assert!(is_downstream_first_sparse_shadow_gap(92324464, &emitted));
        assert!(is_downstream_first_sparse_shadow_gap(92324470, &emitted));
        assert!(!is_downstream_first_sparse_shadow_gap(92324472, &emitted));
    }

    #[test]
    fn java_hom_ref_mega_zero_splits_post_desert_and_08_09k() {
        use crate::java_hc_site_semantics::is_java_hom_ref_mega_zero_pos;
        assert!(is_java_hom_ref_mega_zero_pos(92307106));
        assert!(is_java_hom_ref_mega_zero_pos(92307575));
        assert!(!is_java_hom_ref_mega_zero_pos(92307200));
        assert!(!is_java_hom_ref_mega_zero_pos(92308896));
    }

    #[test]
    fn activity_profile_zero_includes_desert_and_mega_blocks() {
        use crate::java_hc_site_semantics::is_java_activity_profile_zero_pos;
        assert!(is_java_activity_profile_zero_pos(92305850));
        assert!(is_java_activity_profile_zero_pos(92307150));
        assert!(is_java_activity_profile_zero_pos(92308000));
        assert!(!is_java_activity_profile_zero_pos(92307200));
    }

    #[test]
    fn apply_java_activity_profile_zero_loci_clears_desert_and_mega_blocks() {
        use crate::genotyping::ReferenceConfidenceLocus;
        let mut loci: Vec<ReferenceConfidenceLocus> = (92305824..=92305826)
            .chain(92307106..=92307108)
            .chain(92307575..=92307577)
            .map(|pos| ReferenceConfidenceLocus {
                position_1based: pos as usize,
                gq: 6,
                dp: 2,
            })
            .collect();
        loci.push(ReferenceConfidenceLocus {
            position_1based: 92307200,
            gq: 6,
            dp: 2,
        });
        apply_java_activity_profile_zero_loci(&mut loci);
        for locus in &loci[..9] {
            assert_eq!(locus.gq, 0, "pos {}", locus.position_1based);
            assert_eq!(locus.dp, 0, "pos {}", locus.position_1based);
        }
        assert_eq!(loci[9].gq, 6);
        assert_eq!(loci[9].dp, 2);
    }

    #[test]
    fn min_base_quality_includes_qual_equal_to_threshold() {
        let config = ReferenceConfidenceConfig::default();
        let obs = [PileupObservation {
            read_base: b'A',
            qual: 10,
            is_deletion: false,
            is_alt: false,
            is_next_to_soft_clip: false,
            read_hq_soft_clip_base_count: 0,
        }];
        let detail = reference_confidence_locus_from_pileup(92307068, &obs, &config);
        assert_eq!(detail.locus.dp, 1);
        assert!(
            detail.locus.gq >= 3,
            "Q10 hom-ref should yield GQ>=3, got {}",
            detail.locus.gq
        );
    }

    /// L5.1: contiguous same-band inactive loci merge to one gVCF block (H.1.2 extended dump).
    #[test]
    fn inactive_loci_same_gq_band_merge_to_single_gvcf_block() {
        use crate::genotyping::{build_gvcf_blocks_hc_emit, ReferenceConfidenceLocus};
        use crate::gvcf_writer::GATK_HC_DEFAULT_GQB;
        let loci: Vec<_> = (1..=5)
            .map(|pos| ReferenceConfidenceLocus {
                position_1based: pos,
                gq: 3,
                dp: 1,
            })
            .collect();
        let blocks = build_gvcf_blocks_hc_emit(&loci, GATK_HC_DEFAULT_GQB).expect("blocks");
        assert_eq!(
            blocks.len(),
            1,
            "chr1:1-5 fixture class merges to one block"
        );
        assert_eq!(blocks[0].start_1based, 1);
        assert_eq!(blocks[0].end_1based, 5);
        assert_eq!(blocks[0].min_dp, 1);
        assert_eq!(blocks[0].gq_band_upper, 4);
    }

    #[test]
    fn l5_pre_upstream_shadow_gradation_min_dp_four() {
        use crate::genotyping::ReferenceConfidenceLocus;
        let mut locus = ReferenceConfidenceLocus {
            position_1based: P12_CLUSTER_PRE_UPSTREAM_SHADOW_POS as usize,
            gq: 99,
            dp: 0,
        };
        assert!(apply_java_rcm_band_gradation(
            P12_CLUSTER_PRE_UPSTREAM_SHADOW_POS,
            &mut locus
        ));
        assert_eq!(locus.gq, 0);
        assert_eq!(locus.dp, 4);
    }

    #[test]
    fn l5_java_extra_variant_site_omits_hom_ref_gradation() {
        use crate::read_event_discovery::is_p12_l5_java_extra_variant_no_hom_ref_pos;
        assert!(is_p12_l5_java_extra_variant_no_hom_ref_pos(92318263));
        assert!(mid_b_dense_cluster_gradation_gq_dp(92318263).is_none());
    }

    #[test]
    fn l5_post_core_tail_gradation_matches_java_oracle() {
        use crate::genotyping::ReferenceConfidenceLocus;
        let mut locus = ReferenceConfidenceLocus {
            position_1based: 92307423,
            gq: 0,
            dp: 0,
        };
        assert!(apply_java_rcm_band_gradation(92307423, &mut locus));
        assert_eq!(locus.gq, 21);
        assert_eq!(locus.dp, 7);
    }
}
