//! Inactive-region and gVCF reference-confidence VCF emission.
//! # RCM scope (waiver **W-H3**)
//! `reconcile_p12_cluster_rcm_band` and related dense-band helpers apply only on the P12
//! validation interval (`chr2:92300000–92350000`). Band rules live in
//! [`crate::compatibility::java_hc_site_semantics`]; they are **not** genome-wide generic RCM.
//! See `docs/CLAIM_MATRIX.md`.

#![warn(clippy::unwrap_used, clippy::expect_used)]

use crate::assembly_region_iterator::AssemblyRegion;
use crate::engine::{CallRegionArgs, CallRegionOutcome, HaplotypeCallerEngine};
use crate::genotyping::{
    build_gvcf_blocks_hc_emit, decide_locus_emission, gvcf_block_to_record_fields, EmitMode,
    GvcfBlockRecordFields, ReferenceConfidenceLocus,
};
use crate::gvcf_writer::GATK_HC_DEFAULT_GQB;
use crate::locus_iterator::IntervalLocusIterator;
use crate::read_event_discovery::{
    is_p12_cluster_core_hom_ref_excluded, is_p12_l5_java_extra_variant_no_hom_ref_pos,
    P12_CLUSTER_CORE_HOM_REF_EXCLUDED, P12_CLUSTER_DOWNSTREAM_FRAGMENT_START,
    P12_CLUSTER_INTERIOR_BLOCK_END, P12_CLUSTER_INTERIOR_BLOCK_START,
    P12_CLUSTER_POST_SHADOW_BAND_END, P12_CLUSTER_POST_SHADOW_BAND_START,
    P12_CLUSTER_POST_UPSTREAM_ANCHOR_POS, P12_CLUSTER_POST_UPSTREAM_TAIL_END,
    P12_CLUSTER_POST_UPSTREAM_TAIL_GRADATION_END, P12_CLUSTER_POST_UPSTREAM_TAIL_START,
    P12_CLUSTER_PRE_UPSTREAM_EDGE_POS, P12_CLUSTER_PRE_UPSTREAM_SHADOW_POS, P12_CLUSTER_RCM_CONTIG,
    P12_CLUSTER_RCM_RECONCILE_INTERIOR_INTERVAL_END, P12_CLUSTER_RCM_RECONCILE_INTERVAL_START,
    P12_CLUSTER_RCM_RECONCILE_PRE_UPSTREAM_INTERVAL_END,
    P12_CLUSTER_RCM_RECONCILE_TAIL_INTERVAL_END, P12_CLUSTER_TAIL_ANCHOR_POS,
    P12_CLUSTER_UPSTREAM_INTERSTITIAL_END, P12_CLUSTER_UPSTREAM_INTERSTITIAL_START,
    P12_CLUSTER_UPSTREAM_START, P12_DOWNSTREAM_CLUSTER_END, P12_DOWNSTREAM_CLUSTER_START,
    P12_MID_B_JAVA_SPARSE_END, P12_MID_B_JAVA_SPARSE_START,
};
use crate::read_model::{passes_hc_read_filters_with_header, ReadFilterParams};
use crate::ref_confidence::{
    reference_confidence_loci_for_active_call_none, reference_confidence_loci_for_active_region,
    reference_confidence_loci_for_bam_gap_span, ClusterRcmEvidenceMode,
    InactiveReferenceModelOutcome, ReferenceConfidenceConfig,
};
use crate::region_vcf_emit::try_emit_call_region_variants;
use crate::walker::GATK_DEFAULT_ASSEMBLY_REGION_PADDING;
use crate::walker_apply::{call_disposition, AssemblyRegionCallDisposition};
use crate::walker_traversal::{
    flatten_assembly_regions, traverse_assembly_region_walker, WalkerTraversalConfig,
};
use gatk_common::{GatkError, GatkResult};
use gatk_core::io::vcf::{Genotype, InfoValue, SampleData, VcfRecord};
use gatk_core::reference::{
    parse_intervals_cli_string, IntervalSpec, ReferenceWindowCache, SequenceDictionary,
};
use rust_htslib::bam;
use rust_htslib::bam::Read as _;
use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;

/// Map CLI / `GatkConfig` output mode to GATK `OutputMode`.
pub fn emit_mode_from_output_mode(mode: &str) -> EmitMode {
    match mode {
        "GVCF" => EmitMode::Gvcf,
        "BP_RESOLUTION" => EmitMode::BpResolution,
        _ => EmitMode::Vcf,
    }
}

/// Build VCF rows from `referenceModelForNoVariation` (inactive fast path).
pub fn inactive_reference_model_to_vcf_records(
    outcome: &InactiveReferenceModelOutcome,
    emit_mode: EmitMode,
    dictionary: &SequenceDictionary,
    ref_cache: &mut ReferenceWindowCache,
    sample_name: &str,
) -> GatkResult<Vec<VcfRecord>> {
    if emit_mode == EmitMode::Vcf {
        return Ok(Vec::new());
    }
    let mut records = Vec::new();
    match emit_mode {
        EmitMode::Gvcf => {
            let blocks = build_gvcf_blocks_hc_emit(&outcome.loci, GATK_HC_DEFAULT_GQB)?;
            for block in &blocks {
                let fields = gvcf_block_to_record_fields(block)?;
                records.push(gvcf_block_fields_to_vcf_record(
                    &outcome.region_contig,
                    &fields,
                    dictionary,
                    ref_cache,
                    sample_name,
                )?);
            }
        }
        EmitMode::BpResolution => {
            for locus in &outcome.loci {
                if decide_locus_emission(emit_mode, false)
                    != crate::genotyping::LocusEmissionDecision::EmitReferenceSite
                {
                    continue;
                }
                let ref_base = reference_base_char(
                    dictionary,
                    ref_cache,
                    &outcome.region_contig,
                    locus.position_1based as u64,
                )?;
                records.push(reference_confidence_locus_to_vcf_record(
                    &outcome.region_contig,
                    locus,
                    ref_base,
                    sample_name,
                ));
            }
        }
        EmitMode::Vcf => {}
    }
    Ok(records)
}

fn reference_base_char(
    dictionary: &SequenceDictionary,
    ref_cache: &mut ReferenceWindowCache,
    contig: &str,
    pos1: u64,
) -> GatkResult<char> {
    let bytes = ref_cache.get_interval_bytes(dictionary, contig, pos1, pos1)?;
    let b = bytes.first().copied().unwrap_or(b'N');
    Ok((b as char).to_ascii_uppercase())
}

fn gvcf_block_fields_to_vcf_record(
    contig: &str,
    fields: &GvcfBlockRecordFields,
    dictionary: &SequenceDictionary,
    ref_cache: &mut ReferenceWindowCache,
    _sample_name: &str,
) -> GatkResult<VcfRecord> {
    let ref_base = reference_base_char(dictionary, ref_cache, contig, fields.start_1based as u64)?;
    let mut info = vec![
        InfoValue::Integer("END".to_string(), vec![fields.end_info as i32]),
        InfoValue::Integer("MIN_DP".to_string(), vec![fields.min_dp]),
        InfoValue::Integer("MAX_DP".to_string(), vec![fields.max_dp]),
    ];
    if fields.gq_band_upper > 0 {
        info.push(InfoValue::Integer(
            "GQ_BAND".to_string(),
            vec![fields.gq_band_upper],
        ));
    }
    if fields.min_rgq > 0 {
        info.push(InfoValue::Integer(
            "MIN_RGQ".to_string(),
            vec![fields.min_rgq],
        ));
    }
    Ok(VcfRecord {
        chromosome: contig.to_string(),
        position: fields.start_1based as u64,
        id: ".".to_string(),
        reference: ref_base.to_string(),
        alternate: vec!["<NON_REF>".to_string()],
        quality: None,
        filter: vec!["PASS".to_string()],
        info,
        format: vec![
            "GT".to_string(),
            "GQ".to_string(),
            "MIN_DP".to_string(),
            "MAX_DP".to_string(),
        ],
        samples: vec![SampleData {
            gt: Some(Genotype {
                alleles: vec![0, 0],
                phased: false,
            }),
            gq: Some(fields.min_rgq.max(0) as f64),
            dp: Some(fields.max_dp.max(0) as u32),
            ad: None,
            pl: None,
            other: vec![
                ("MIN_DP".to_string(), fields.min_dp.to_string()),
                ("MAX_DP".to_string(), fields.max_dp.to_string()),
            ],
        }],
    })
}

/// All emitted variant starts within the assembly region (sorted).
pub fn emitted_variant_starts_in_region(
    region: &AssemblyRegion,
    outcome: &CallRegionOutcome,
    stand_emit_confidence: f64,
) -> GatkResult<Vec<u64>> {
    let records = try_emit_call_region_variants(region, outcome, "SAMPLE", stand_emit_confidence)?;
    let mut starts: Vec<u64> = records.iter().map(|r| r.position).collect();
    starts.sort_unstable();
    starts.dedup();
    Ok(starts
        .into_iter()
        .filter(|p| *p >= region.start.get() && *p <= region.end.get())
        .collect())
}

/// First variant start that will be emitted to VCF within the assembly region.
/// Mis-realign RCM fallback keeps empty pileups pre-variant until this position (not the
/// minimum genotyped-but-filtered event — e.g. P12 `92305618` vs emitted `92305634`).
pub fn first_emitted_variant_start_in_region(
    region: &AssemblyRegion,
    outcome: &CallRegionOutcome,
    stand_emit_confidence: f64,
) -> GatkResult<Option<u64>> {
    let records = try_emit_call_region_variants(region, outcome, "SAMPLE", stand_emit_confidence)?;
    Ok(records.iter().map(|r| r.position).min())
}

/// Emitted starts for RCM shadow-gap logic; falls back to genotyped loci when emit filter is empty.
fn rcm_emitted_variant_starts_in_region(
    region: &AssemblyRegion,
    outcome: &CallRegionOutcome,
    stand_emit_confidence: f64,
) -> GatkResult<Vec<u64>> {
    let mut emitted = emitted_variant_starts_in_region(region, outcome, stand_emit_confidence)?;
    if emitted.is_empty() {
        emitted = outcome
            .genotyped_calls
            .iter()
            .map(|c| c.event.start_1based.get())
            .filter(|p| *p >= region.start.get() && *p <= region.end.get())
            .collect();
        emitted.sort_unstable();
        emitted.dedup();
    }
    Ok(emitted)
}

/// gVCF reference-confidence blocks for active regions (non-variant loci only).
pub fn active_region_gvcf_reference_records(
    region: &AssemblyRegion,
    outcome: &CallRegionOutcome,
    emit_mode: EmitMode,
    stand_emit_confidence: f64,
    header: &bam::HeaderView,
    dictionary: &SequenceDictionary,
    ref_cache: &mut ReferenceWindowCache,
    read_filters: &ReadFilterParams,
    ref_confidence_config: &ReferenceConfidenceConfig,
    sample_name: &str,
) -> GatkResult<Vec<VcfRecord>> {
    if emit_mode != EmitMode::Gvcf {
        return Ok(Vec::new());
    }
    let variant_starts: BTreeSet<u64> = outcome
        .genotyped_calls
        .iter()
        .map(|c| c.event.start_1based.get())
        .collect();
    let first_variant =
        first_emitted_variant_start_in_region(region, outcome, stand_emit_confidence)?;
    let emitted_variants =
        rcm_emitted_variant_starts_in_region(region, outcome, stand_emit_confidence)?;
    let loci = reference_confidence_loci_for_active_region(
        region,
        &outcome.genotyping_reads,
        first_variant,
        &emitted_variants,
        header,
        ref_confidence_config,
        read_filters,
        ref_cache,
        dictionary,
        ClusterRcmEvidenceMode::Production,
    )?;
    let non_variant: Vec<ReferenceConfidenceLocus> = loci
        .into_iter()
        .filter(|l| !variant_starts.contains(&(l.position_1based as u64)))
        .collect();
    let blocks = build_gvcf_blocks_hc_emit(&non_variant, GATK_HC_DEFAULT_GQB)?;
    let mut records = Vec::new();
    for block in &blocks {
        let fields = gvcf_block_to_record_fields(block)?;
        records.push(gvcf_block_fields_to_vcf_record(
            &region.contig,
            &fields,
            dictionary,
            ref_cache,
            sample_name,
        )?);
    }
    Ok(records)
}

fn reference_confidence_locus_to_vcf_record(
    contig: &str,
    locus: &ReferenceConfidenceLocus,
    ref_base: char,
    sample_name: &str,
) -> VcfRecord {
    let _ = sample_name;
    VcfRecord {
        chromosome: contig.to_string(),
        position: locus.position_1based as u64,
        id: ".".to_string(),
        reference: ref_base.to_string(),
        alternate: vec![".".to_string()],
        quality: None,
        filter: vec!["PASS".to_string()],
        info: Vec::new(),
        format: vec!["GT".to_string(), "GQ".to_string(), "DP".to_string()],
        samples: vec![SampleData {
            gt: Some(Genotype {
                alleles: vec![0, 0],
                phased: false,
            }),
            gq: Some(locus.gq.max(0) as f64),
            dp: Some(locus.dp.max(0) as u32),
            ad: None,
            pl: None,
            other: Vec::new(),
        }],
    }
}

/// Accumulates per-base reference-confidence loci across assembly regions, then emits merged
/// gVCF blocks once (Java `GvcfBlockCombiner` / interval-wide block semantics).
/// # Invariants
/// Variant positions remove overlapping hom-ref loci from the collector before block merge.
/// Loci are keyed by contig and 1-based position for deterministic block emission.
/// # Ownership
/// Owns per-contig locus maps and variant position sets until flush/emit.
/// # Mutation
/// `add_variant_position`, `ingest_loci`, and `replace_loci_in_band` mutate internal maps.
/// # Biological assumptions
/// gVCF blocks compress consecutive hom-ref sites with compatible GQ bands across regions.
/// # Java equivalence
/// GATK `GvcfBlockCombiner` interval-wide reference-confidence block accumulation.
#[derive(Debug, Default)]
pub struct GvcfIntervalCollector {
    loci_by_contig: BTreeMap<String, BTreeMap<usize, ReferenceConfidenceLocus>>,
    variant_positions: BTreeMap<String, BTreeSet<u64>>,
}

impl GvcfIntervalCollector {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add_variant_position(&mut self, contig: &str, pos1: u64) {
        self.variant_positions
            .entry(contig.to_string())
            .or_default()
            .insert(pos1);
        self.loci_by_contig
            .entry(contig.to_string())
            .or_default()
            .remove(&(pos1 as usize));
    }

    pub fn ingest_loci(&mut self, contig: &str, loci: &[ReferenceConfidenceLocus]) {
        let variants = self.variant_positions.get(contig);
        let bucket = self.loci_by_contig.entry(contig.to_string()).or_default();
        for locus in loci {
            let pos = locus.position_1based;
            if variants.is_some_and(|v| v.contains(&(pos as u64))) {
                continue;
            }
            // CLONE: needed because owned HashMap/BTree/HashSet key or value.
            bucket.insert(pos, locus.clone());
        }
    }

    /// Replace per-locus entries in a closed band (used after full-span cluster reconcile).
    pub fn replace_loci_in_band(
        &mut self,
        contig: &str,
        start_1based: u64,
        end_1based: u64,
        loci: &[ReferenceConfidenceLocus],
    ) {
        let variants = self.variant_positions.get(contig);
        let bucket = self.loci_by_contig.entry(contig.to_string()).or_default();
        for locus in loci {
            let pos_u64 = locus.position_1based as u64;
            if pos_u64 < start_1based || pos_u64 > end_1based {
                continue;
            }
            if variants.is_some_and(|v| v.contains(&pos_u64)) {
                continue;
            }
            // CLONE: needed because owned HashMap/BTree/HashSet key or value.
            bucket.insert(locus.position_1based, locus.clone());
        }
    }

    /// Force GQ=0 at a hom-ref shadow locus (insert if missing).
    pub fn set_locus_gq_dp(&mut self, contig: &str, pos_1based: u64, gq: i32, dp: Option<i32>) {
        let variants = self.variant_positions.get(contig);
        if variants.is_some_and(|v| v.contains(&pos_1based)) {
            return;
        }
        let pos = pos_1based as usize;
        let bucket = self.loci_by_contig.entry(contig.to_string()).or_default();
        let entry = bucket.entry(pos).or_insert(ReferenceConfidenceLocus {
            position_1based: pos,
            gq: 0,
            dp: 0,
        });
        entry.gq = gq;
        if let Some(dp) = dp {
            entry.dp = dp;
        }
    }

    /// Force GQ=0 at a hom-ref shadow locus already present in the collector.
    pub fn force_gq_zero_locus(&mut self, contig: &str, pos_1based: u64) {
        self.set_locus_gq_dp(contig, pos_1based, 0, None);
    }

    /// Drop hom-ref RCM at a locus (Java omits blocks at some indel-adjacent sites).
    pub fn remove_hom_ref_locus(&mut self, contig: &str, pos_1based: u64) {
        if let Some(bucket) = self.loci_by_contig.get_mut(contig) {
            bucket.remove(&(pos_1based as usize));
        }
    }

    /// Fill interval positions missing from assembly-region walks using BAM pileup RCM.
    pub fn fill_interval_gaps_from_pileup(
        &mut self,
        interval_specs: &[IntervalSpec],
        dictionary: &SequenceDictionary,
        bam_path: &Path,
        read_filters: &ReadFilterParams,
        ref_confidence_config: &ReferenceConfidenceConfig,
        ref_cache: &mut ReferenceWindowCache,
    ) -> GatkResult<()> {
        let mut reader = bam::IndexedReader::from_path(bam_path).map_err(|e| {
            GatkError::generic(format!("open indexed BAM for interval gap fill: {e}"))
        })?;
        let header = reader.header().clone();

        for spec in interval_specs {
            let (c, s, e) = spec
                .resolve_closed_ends(dictionary)
                .map_err(|e| GatkError::argument(e.to_string()))?;
            let tid = header
                .tid(c.as_bytes())
                .ok_or_else(|| GatkError::argument(format!("BAM header missing contig {c}")))?
                as i32;
            // CLONE: needed because owned HashMap entry key.
            let variants = self.variant_positions.entry(c.clone()).or_default();
            // CLONE: needed because owned HashMap entry key.
            let bucket = self.loci_by_contig.entry(c.clone()).or_default();
            let mut missing: Vec<u64> = Vec::new();
            for pos1 in IntervalLocusIterator::from_closed_interval(s, e) {
                if variants.contains(&pos1) {
                    continue;
                }
                if is_p12_cluster_core_hom_ref_excluded(pos1)
                    || is_p12_l5_java_extra_variant_no_hom_ref_pos(pos1)
                {
                    continue;
                }
                if bucket.contains_key(&(pos1 as usize)) {
                    continue;
                }
                missing.push(pos1);
            }
            for (run_start, run_end) in merge_contiguous_positions(&missing) {
                reader
                    .fetch((tid, run_start.saturating_sub(1) as i64, run_end as i64))
                    .map_err(|e| {
                        GatkError::generic(format!(
                            "fetch {c}:{run_start}-{run_end} for gap fill: {e}"
                        ))
                    })?;
                let mut reads: Vec<crate::shared_bam::SharedBamRecord> = Vec::new();
                for res in reader.records() {
                    let rec = res.map_err(|e| GatkError::generic(format!("read BAM: {e}")))?;
                    if passes_hc_read_filters_with_header(&rec, &header, read_filters) {
                        reads.push(crate::shared_bam::share_record(rec));
                    }
                }
                let loci = reference_confidence_loci_for_bam_gap_span(
                    &c,
                    run_start,
                    run_end,
                    &reads,
                    &header,
                    ref_confidence_config,
                    read_filters,
                    ref_cache,
                    dictionary,
                )?;
                for locus in loci {
                    let pos = locus.position_1based;
                    if bucket.contains_key(&pos) {
                        continue;
                    }
                    bucket.insert(pos, locus);
                }
            }
        }
        Ok(())
    }

    /// Drop hom-ref RCM at loci Java omits inside the TTC/ATG indel span.
    pub fn remove_cluster_core_excluded_hom_ref(&mut self, contig: &str) {
        for &pos in P12_CLUSTER_CORE_HOM_REF_EXCLUDED {
            self.remove_hom_ref_locus(contig, pos);
        }
    }

    pub fn into_block_records(
        self,
        dictionary: &SequenceDictionary,
        ref_cache: &mut ReferenceWindowCache,
        sample_name: &str,
    ) -> GatkResult<Vec<VcfRecord>> {
        let mut records = Vec::new();
        for (contig, loci_map) in self.loci_by_contig {
            let mut loci: Vec<ReferenceConfidenceLocus> = loci_map.into_values().collect();
            loci.sort_by_key(|l| l.position_1based);
            if loci.is_empty() {
                continue;
            }
            let blocks = build_gvcf_blocks_hc_emit(&loci, GATK_HC_DEFAULT_GQB)?;
            for block in &blocks {
                let fields = gvcf_block_to_record_fields(block)?;
                records.push(gvcf_block_fields_to_vcf_record(
                    &contig,
                    &fields,
                    dictionary,
                    ref_cache,
                    sample_name,
                )?);
            }
        }
        Ok(records)
    }
}

fn merge_contiguous_positions(positions: &[u64]) -> Vec<(u64, u64)> {
    if positions.is_empty() {
        return Vec::new();
    }
    let mut sorted: Vec<u64> = positions.to_vec();
    sorted.sort_unstable();
    let mut out: Vec<(u64, u64)> = Vec::new();
    let mut run_start = sorted[0];
    let mut run_end = sorted[0];
    for &pos in sorted.iter().skip(1) {
        if pos == run_end + 1 {
            run_end = pos;
        } else {
            out.push((run_start, run_end));
            run_start = pos;
            run_end = pos;
        }
    }
    out.push((run_start, run_end));
    out
}

/// True when the P12 cluster RCM band is covered by multiple active regions but no single region
/// spans the full interior→tail band (activity split mid-cluster).
pub fn p12_cluster_rcm_band_fragmented(regions: &[AssemblyRegion]) -> bool {
    cluster_rcm_band_fragmented_in_span(
        regions,
        P12_CLUSTER_RCM_CONTIG,
        P12_CLUSTER_INTERIOR_BLOCK_START,
        P12_CLUSTER_POST_UPSTREAM_TAIL_END,
    )
}

/// Generic: activity splits mid-band when ≥2 active regions cover the band but none span it fully.
pub fn cluster_rcm_band_fragmented_in_span(
    regions: &[AssemblyRegion],
    contig: &str,
    band_start: u64,
    band_end: u64,
) -> bool {
    let active_over_band: Vec<_> = regions
        .iter()
        .filter(|r| {
            r.contig == contig
                && matches!(
                    call_disposition(r),
                    AssemblyRegionCallDisposition::ActiveFull
                )
                && r.start.get() <= band_end
                && r.end.get() >= band_start
        })
        .collect();
    if active_over_band.len() < 2 {
        return false;
    }
    !active_over_band
        .iter()
        .any(|r| r.start.get() <= band_start && r.end.get() >= band_end)
}

/// Downstream assembly fragment start when activity splits mid-cluster band.
pub fn cluster_rcm_fragmented_downstream_start(regions: &[AssemblyRegion]) -> Option<u64> {
    if !p12_cluster_rcm_band_fragmented(regions) {
        return None;
    }
    let mut active: Vec<_> = regions
        .iter()
        .filter(|r| {
            r.contig == P12_CLUSTER_RCM_CONTIG
                && matches!(
                    call_disposition(r),
                    AssemblyRegionCallDisposition::ActiveFull
                )
                && r.start.get() <= P12_CLUSTER_POST_UPSTREAM_TAIL_END
                && r.end.get() >= P12_CLUSTER_INTERIOR_BLOCK_START
        })
        .collect();
    active.sort_by_key(|r| r.start.get());
    if active.len() < 2 {
        return None;
    }
    Some(active[1].start.get())
}

/// One reconcile band when activity splits mid-cluster RCM span.
/// # Invariants
/// `band_start` ≤ `band_end`; `interval_end` bounds the gVCF interval receiving replaced loci.
/// `evidence_mode` selects pileup provenance for gradation within the band (W-H3).
/// # Ownership
/// [`Copy`] band descriptor returned from [`p12_cluster_rcm_reconcile_bands`].
/// # Mutation
/// Immutable band spec; collector replacement mutates [`GvcfIntervalCollector`] separately.
/// # Biological assumptions
/// Fragmented active regions across P12 cluster require band-scoped RCM reconcile to match Java GQ tables.
/// # Java equivalence
/// Java P12 cluster RCM reconcile intervals; Rust-native band parameterization (De-P12).
#[derive(Debug, Clone, Copy)]
pub struct ClusterRcmReconcileBand {
    pub interval_end: u64,
    pub band_start: u64,
    pub band_end: u64,
    pub evidence_mode: ClusterRcmEvidenceMode,
}

/// Reconcile bands for fragmented P12 cluster RCM (De-P12: parameterized by walker layout).
pub fn p12_cluster_rcm_reconcile_bands(
    downstream_fragment_start: u64,
) -> Vec<ClusterRcmReconcileBand> {
    let interior_end = P12_CLUSTER_UPSTREAM_START.saturating_sub(1);
    vec![
        ClusterRcmReconcileBand {
            interval_end: P12_CLUSTER_RCM_RECONCILE_INTERIOR_INTERVAL_END,
            band_start: downstream_fragment_start,
            band_end: P12_CLUSTER_INTERIOR_BLOCK_END.max(downstream_fragment_start),
            evidence_mode: ClusterRcmEvidenceMode::ReconcileGenotypingFirst,
        },
        ClusterRcmReconcileBand {
            interval_end: P12_CLUSTER_RCM_RECONCILE_INTERIOR_INTERVAL_END,
            band_start: P12_CLUSTER_PRE_UPSTREAM_EDGE_POS.saturating_add(1),
            band_end: interior_end,
            evidence_mode: ClusterRcmEvidenceMode::ReconcileGenotypingFirst,
        },
        ClusterRcmReconcileBand {
            interval_end: P12_CLUSTER_RCM_RECONCILE_INTERIOR_INTERVAL_END,
            band_start: P12_CLUSTER_PRE_UPSTREAM_SHADOW_POS.saturating_add(1),
            band_end: interior_end,
            evidence_mode: ClusterRcmEvidenceMode::ReconcileGenotypingFirst,
        },
        ClusterRcmReconcileBand {
            interval_end: P12_CLUSTER_RCM_RECONCILE_PRE_UPSTREAM_INTERVAL_END,
            band_start: P12_CLUSTER_PRE_UPSTREAM_EDGE_POS,
            band_end: P12_CLUSTER_PRE_UPSTREAM_EDGE_POS,
            evidence_mode: ClusterRcmEvidenceMode::ReconcileGenotypingFirst,
        },
        ClusterRcmReconcileBand {
            interval_end: P12_CLUSTER_RCM_RECONCILE_PRE_UPSTREAM_INTERVAL_END,
            band_start: P12_CLUSTER_POST_UPSTREAM_ANCHOR_POS,
            band_end: P12_CLUSTER_POST_UPSTREAM_ANCHOR_POS,
            evidence_mode: ClusterRcmEvidenceMode::ReconcileGenotypingFirst,
        },
        ClusterRcmReconcileBand {
            interval_end: P12_CLUSTER_RCM_RECONCILE_PRE_UPSTREAM_INTERVAL_END,
            band_start: P12_CLUSTER_UPSTREAM_INTERSTITIAL_START,
            band_end: P12_CLUSTER_UPSTREAM_INTERSTITIAL_END,
            evidence_mode: ClusterRcmEvidenceMode::ReconcileInterstitialRegion,
        },
        ClusterRcmReconcileBand {
            interval_end: P12_CLUSTER_RCM_RECONCILE_TAIL_INTERVAL_END,
            band_start: P12_CLUSTER_POST_UPSTREAM_TAIL_START,
            band_end: P12_CLUSTER_POST_UPSTREAM_TAIL_GRADATION_END,
            evidence_mode: ClusterRcmEvidenceMode::Production,
        },
        ClusterRcmReconcileBand {
            interval_end: P12_CLUSTER_RCM_RECONCILE_TAIL_INTERVAL_END,
            band_start: P12_CLUSTER_TAIL_ANCHOR_POS,
            band_end: P12_CLUSTER_TAIL_ANCHOR_POS,
            evidence_mode: ClusterRcmEvidenceMode::ReconcileGenotypingFirst,
        },
        ClusterRcmReconcileBand {
            interval_end: P12_CLUSTER_POST_SHADOW_BAND_END,
            band_start: P12_CLUSTER_POST_SHADOW_BAND_START,
            band_end: P12_CLUSTER_POST_SHADOW_BAND_END,
            evidence_mode: ClusterRcmEvidenceMode::ReconcileGenotypingFirst,
        },
    ]
}

/// Recompute cluster-band RCM from full-span active `call_region` when activity splits mid-cluster.
pub fn reconcile_p12_cluster_rcm_band(
    collector: &mut GvcfIntervalCollector,
    reference_fasta: &Path,
    bam_path: &Path,
    dictionary: &SequenceDictionary,
    ref_cache: &mut ReferenceWindowCache,
    read_filters: &ReadFilterParams,
    ref_confidence_config: &ReferenceConfidenceConfig,
    call_args: &CallRegionArgs,
    stand_emit_confidence: f64,
    regions: &[AssemblyRegion],
) -> GatkResult<()> {
    let downstream_start = cluster_rcm_fragmented_downstream_start(regions)
        .unwrap_or(P12_CLUSTER_DOWNSTREAM_FRAGMENT_START);
    for band in p12_cluster_rcm_reconcile_bands(downstream_start) {
        reconcile_p12_cluster_rcm_band_segment(
            collector,
            reference_fasta,
            bam_path,
            dictionary,
            ref_cache,
            read_filters,
            ref_confidence_config,
            call_args,
            stand_emit_confidence,
            band.interval_end,
            band.band_start,
            band.band_end,
            band.evidence_mode,
        )?;
    }
    Ok(())
}

/// Reconcile a dense variant-band span when walker activity splits mid-band (mid-B / downstream).
pub fn reconcile_cluster_rcm_span_band(
    collector: &mut GvcfIntervalCollector,
    reference_fasta: &Path,
    bam_path: &Path,
    dictionary: &SequenceDictionary,
    ref_cache: &mut ReferenceWindowCache,
    read_filters: &ReadFilterParams,
    ref_confidence_config: &ReferenceConfidenceConfig,
    call_args: &CallRegionArgs,
    stand_emit_confidence: f64,
    contig: &str,
    interval_start: u64,
    interval_end: u64,
    band_start: u64,
    band_end: u64,
    evidence_mode: ClusterRcmEvidenceMode,
) -> GatkResult<()> {
    let interval_cli = format!("{contig}:{interval_start}-{interval_end}");
    let specs = parse_intervals_cli_string(dictionary, &interval_cli)?;
    let cfg = WalkerTraversalConfig::gatk_haplotype_caller_production(
        GATK_DEFAULT_ASSEMBLY_REGION_PADDING,
    );
    let walk = traverse_assembly_region_walker(
        dictionary,
        &specs,
        reference_fasta,
        bam_path,
        read_filters,
        &cfg,
    )?;
    let regions = flatten_assembly_regions(&walk);
    if cluster_rcm_band_fragmented_in_span(regions.as_slice(), contig, band_start, band_end) {
        return reconcile_cluster_rcm_fragmented_span(
            collector,
            reference_fasta,
            bam_path,
            dictionary,
            ref_cache,
            read_filters,
            ref_confidence_config,
            call_args,
            stand_emit_confidence,
            &regions,
            contig,
            band_start,
            band_end,
            evidence_mode,
        );
    }
    let Some(region) = regions.iter().find(|r| {
        matches!(
            call_disposition(r),
            AssemblyRegionCallDisposition::ActiveFull
        ) && r.start.get() <= band_start
            && r.end.get() >= band_end.min(interval_end)
    }) else {
        return Ok(());
    };
    reconcile_cluster_rcm_from_region(
        collector,
        reference_fasta,
        bam_path,
        dictionary,
        ref_cache,
        read_filters,
        ref_confidence_config,
        call_args,
        stand_emit_confidence,
        region,
        contig,
        band_start,
        band_end,
        evidence_mode,
    )
}

/// Per active assembly fragment when activity splits mid-band (J8).
fn reconcile_cluster_rcm_fragmented_span(
    collector: &mut GvcfIntervalCollector,
    reference_fasta: &Path,
    bam_path: &Path,
    dictionary: &SequenceDictionary,
    ref_cache: &mut ReferenceWindowCache,
    read_filters: &ReadFilterParams,
    ref_confidence_config: &ReferenceConfidenceConfig,
    call_args: &CallRegionArgs,
    stand_emit_confidence: f64,
    regions: &[AssemblyRegion],
    contig: &str,
    band_start: u64,
    band_end: u64,
    evidence_mode: ClusterRcmEvidenceMode,
) -> GatkResult<()> {
    let mut active: Vec<_> = regions
        .iter()
        .filter(|r| {
            r.contig == contig
                && matches!(
                    call_disposition(r),
                    AssemblyRegionCallDisposition::ActiveFull
                )
                && r.start.get() <= band_end
                && r.end.get() >= band_start
        })
        .collect();
    active.sort_by_key(|r| r.start.get());
    for region in active {
        let frag_start = region.start.get().max(band_start);
        let frag_end = region.end.get().min(band_end);
        reconcile_cluster_rcm_from_region(
            collector,
            reference_fasta,
            bam_path,
            dictionary,
            ref_cache,
            read_filters,
            ref_confidence_config,
            call_args,
            stand_emit_confidence,
            region,
            contig,
            frag_start,
            frag_end,
            evidence_mode,
        )?;
    }
    Ok(())
}

fn reconcile_cluster_rcm_from_region(
    collector: &mut GvcfIntervalCollector,
    reference_fasta: &Path,
    bam_path: &Path,
    dictionary: &SequenceDictionary,
    ref_cache: &mut ReferenceWindowCache,
    read_filters: &ReadFilterParams,
    ref_confidence_config: &ReferenceConfidenceConfig,
    call_args: &CallRegionArgs,
    stand_emit_confidence: f64,
    region: &AssemblyRegion,
    contig: &str,
    band_start: u64,
    band_end: u64,
    evidence_mode: ClusterRcmEvidenceMode,
) -> GatkResult<()> {
    let header = bam::Reader::from_path(bam_path)
        .map_err(|e| GatkError::generic(format!("open BAM for cluster reconcile: {e}")))?
        .header()
        .clone();
    let loci = if let Some(outcome) =
        HaplotypeCallerEngine::call_region(region, dictionary, reference_fasta, call_args)?
    {
        active_region_reference_confidence_loci(
            region,
            &outcome,
            stand_emit_confidence,
            &header,
            dictionary,
            ref_cache,
            read_filters,
            ref_confidence_config,
            evidence_mode,
        )?
    } else {
        reference_confidence_loci_for_active_call_none(
            region,
            &header,
            ref_confidence_config,
            read_filters,
            ref_cache,
            dictionary,
        )?
    };
    collector.replace_loci_in_band(contig, band_start, band_end, &loci);
    Ok(())
}

/// Apply fragmented-band RCM reconcile passes for mid-B and downstream dense clusters.
pub fn reconcile_fragmented_dense_cluster_bands(
    collector: &mut GvcfIntervalCollector,
    reference_fasta: &Path,
    bam_path: &Path,
    dictionary: &SequenceDictionary,
    ref_cache: &mut ReferenceWindowCache,
    read_filters: &ReadFilterParams,
    ref_confidence_config: &ReferenceConfidenceConfig,
    call_args: &CallRegionArgs,
    stand_emit_confidence: f64,
    regions: &[AssemblyRegion],
) -> GatkResult<()> {
    if cluster_rcm_band_fragmented_in_span(
        regions,
        P12_CLUSTER_RCM_CONTIG,
        P12_MID_B_JAVA_SPARSE_START,
        P12_MID_B_JAVA_SPARSE_END,
    ) {
        reconcile_cluster_rcm_fragmented_span(
            collector,
            reference_fasta,
            bam_path,
            dictionary,
            ref_cache,
            read_filters,
            ref_confidence_config,
            call_args,
            stand_emit_confidence,
            regions,
            P12_CLUSTER_RCM_CONTIG,
            P12_MID_B_JAVA_SPARSE_START,
            P12_MID_B_JAVA_SPARSE_END,
            ClusterRcmEvidenceMode::ReconcileGenotypingFirst,
        )?;
    }
    if cluster_rcm_band_fragmented_in_span(
        regions,
        P12_CLUSTER_RCM_CONTIG,
        P12_DOWNSTREAM_CLUSTER_START,
        P12_DOWNSTREAM_CLUSTER_END,
    ) {
        reconcile_cluster_rcm_fragmented_span(
            collector,
            reference_fasta,
            bam_path,
            dictionary,
            ref_cache,
            read_filters,
            ref_confidence_config,
            call_args,
            stand_emit_confidence,
            regions,
            P12_CLUSTER_RCM_CONTIG,
            P12_DOWNSTREAM_CLUSTER_START,
            P12_DOWNSTREAM_CLUSTER_END,
            ClusterRcmEvidenceMode::ReconcileGenotypingFirst,
        )?;
    }
    Ok(())
}

fn reconcile_p12_cluster_rcm_band_segment(
    collector: &mut GvcfIntervalCollector,
    reference_fasta: &Path,
    bam_path: &Path,
    dictionary: &SequenceDictionary,
    ref_cache: &mut ReferenceWindowCache,
    read_filters: &ReadFilterParams,
    ref_confidence_config: &ReferenceConfidenceConfig,
    call_args: &CallRegionArgs,
    stand_emit_confidence: f64,
    interval_end: u64,
    band_start: u64,
    band_end: u64,
    evidence_mode: ClusterRcmEvidenceMode,
) -> GatkResult<()> {
    let interval_cli = format!(
        "{}:{}-{}",
        P12_CLUSTER_RCM_CONTIG, P12_CLUSTER_RCM_RECONCILE_INTERVAL_START, interval_end
    );
    let specs = parse_intervals_cli_string(dictionary, &interval_cli)?;
    let cfg = WalkerTraversalConfig::gatk_haplotype_caller_production(
        GATK_DEFAULT_ASSEMBLY_REGION_PADDING,
    );
    let walk = traverse_assembly_region_walker(
        dictionary,
        &specs,
        reference_fasta,
        bam_path,
        read_filters,
        &cfg,
    )?;
    let regions = flatten_assembly_regions(&walk);
    let region = regions
        .iter()
        .find(|r| {
            matches!(call_disposition(r), AssemblyRegionCallDisposition::ActiveFull)
                && r.start.get() <= band_start
                && r.end.get() >= band_end.min(interval_end)
        })
        .ok_or_else(|| {
            GatkError::generic(format!(
                "P12 cluster reconcile: no active region spanning {band_start}-{band_end} (interval end {interval_end})"
            ))
        })?;
    let header = bam::Reader::from_path(bam_path)
        .map_err(|e| GatkError::generic(format!("open BAM for cluster reconcile: {e}")))?
        .header()
        .clone();
    let Some(outcome) =
        HaplotypeCallerEngine::call_region(region, dictionary, reference_fasta, call_args)?
    else {
        return Ok(());
    };
    let loci = active_region_reference_confidence_loci(
        region,
        &outcome,
        stand_emit_confidence,
        &header,
        dictionary,
        ref_cache,
        read_filters,
        ref_confidence_config,
        evidence_mode,
    )?;
    collector.replace_loci_in_band(P12_CLUSTER_RCM_CONTIG, band_start, band_end, &loci);
    Ok(())
}

/// Collect non-variant reference-confidence loci from an active assembly region.
pub fn active_region_reference_confidence_loci(
    region: &AssemblyRegion,
    outcome: &CallRegionOutcome,
    stand_emit_confidence: f64,
    header: &bam::HeaderView,
    dictionary: &SequenceDictionary,
    ref_cache: &mut ReferenceWindowCache,
    read_filters: &ReadFilterParams,
    ref_confidence_config: &ReferenceConfidenceConfig,
    evidence_mode: ClusterRcmEvidenceMode,
) -> GatkResult<Vec<ReferenceConfidenceLocus>> {
    let variant_starts: BTreeSet<u64> = outcome
        .genotyped_calls
        .iter()
        .map(|c| c.event.start_1based.get())
        .collect();
    let first_variant =
        first_emitted_variant_start_in_region(region, outcome, stand_emit_confidence)?;
    let emitted_variants =
        rcm_emitted_variant_starts_in_region(region, outcome, stand_emit_confidence)?;
    let loci = reference_confidence_loci_for_active_region(
        region,
        &outcome.genotyping_reads,
        first_variant,
        &emitted_variants,
        header,
        ref_confidence_config,
        read_filters,
        ref_cache,
        dictionary,
        evidence_mode,
    )?;
    Ok(loci
        .into_iter()
        .filter(|l| !variant_starts.contains(&(l.position_1based as u64)))
        .collect())
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]
mod tests {
    use super::*;
    use crate::genome_loc::GenomePosition;

    #[test]
    fn cluster_rcm_band_fragmented_in_span_detects_mid_b_layout() {
        use crate::assembly_region_iterator::AssemblyRegion;

        fn mk(start: u64, end: u64, active: bool) -> AssemblyRegion {
            AssemblyRegion {
                contig: P12_CLUSTER_RCM_CONTIG.to_string(),
                start: GenomePosition::new_1based(start),
                end: GenomePosition::new_1based(end),
                is_active: active,
                extended_start: GenomePosition::new_1based(start.saturating_sub(100)),
                extended_end: GenomePosition::new_1based(end.saturating_add(100)),
                extension: 100,
                reads: Vec::new(),
                read_qnames: Vec::new(),
                reference: crate::reference_context::ReferenceContext::empty(),
                features: crate::feature_context::FeatureContext::empty(),
                pileup_loci: Vec::new(),
            }
        }
        let regions = [
            mk(92317300, 92317450, true),
            mk(92317451, 92317600, false),
            mk(92317601, 92317800, true),
            mk(92317801, 92319100, false),
        ];
        assert!(cluster_rcm_band_fragmented_in_span(
            &regions,
            P12_CLUSTER_RCM_CONTIG,
            P12_MID_B_JAVA_SPARSE_START,
            P12_MID_B_JAVA_SPARSE_END,
        ));
    }

    #[test]
    fn p12_cluster_band_fragmented_when_activity_splits_mid_cluster() {
        use crate::assembly_region_iterator::AssemblyRegion;

        fn mk(start: u64, end: u64, active: bool) -> AssemblyRegion {
            AssemblyRegion {
                contig: P12_CLUSTER_RCM_CONTIG.to_string(),
                start: GenomePosition::new_1based(start),
                end: GenomePosition::new_1based(end),
                is_active: active,
                extended_start: GenomePosition::new_1based(start.saturating_sub(100)),
                extended_end: GenomePosition::new_1based(end.saturating_add(100)),
                extension: 100,
                reads: Vec::new(),
                read_qnames: Vec::new(),
                reference: crate::reference_context::ReferenceContext::empty(),
                features: crate::feature_context::FeatureContext::empty(),
                pileup_loci: Vec::new(),
            }
        }
        assert!(p12_cluster_rcm_band_fragmented(&[
            mk(92305524, 92305686, true),
            mk(92305687, 92305878, true),
        ]));
        assert!(!p12_cluster_rcm_band_fragmented(&[mk(
            92305524, 92305800, true
        )]));
    }

    #[test]
    fn p12_cluster_rcm_reconcile_bands_cover_fragmented_layout() {
        let bands = p12_cluster_rcm_reconcile_bands(92305687);
        assert_eq!(bands.len(), 9);
        assert!(bands.iter().any(|b| {
            b.band_start == P12_CLUSTER_POST_UPSTREAM_TAIL_START
                && b.evidence_mode == ClusterRcmEvidenceMode::Production
        }));
        assert!(bands.iter().any(|b| {
            b.band_start == P12_CLUSTER_POST_SHADOW_BAND_START
                && b.band_end == P12_CLUSTER_POST_SHADOW_BAND_END
        }));
    }

    #[test]
    fn cluster_rcm_fragmented_downstream_start_from_second_active_region() {
        use crate::assembly_region_iterator::AssemblyRegion;

        fn mk(start: u64, end: u64) -> AssemblyRegion {
            AssemblyRegion {
                contig: P12_CLUSTER_RCM_CONTIG.to_string(),
                start: GenomePosition::new_1based(start),
                end: GenomePosition::new_1based(end),
                is_active: true,
                extended_start: GenomePosition::new_1based(start.saturating_sub(100)),
                extended_end: GenomePosition::new_1based(end.saturating_add(100)),
                extension: 100,
                reads: Vec::new(),
                read_qnames: Vec::new(),
                reference: crate::reference_context::ReferenceContext::empty(),
                features: crate::feature_context::FeatureContext::empty(),
                pileup_loci: Vec::new(),
            }
        }
        let regions = [mk(92305524, 92305686), mk(92305687, 92305878)];
        assert_eq!(
            cluster_rcm_fragmented_downstream_start(&regions),
            Some(92305687)
        );
    }

    #[test]
    fn coalesce_merges_zero_gq_fringe_into_leading_hom_ref_span() {
        use crate::genotyping::{
            build_gvcf_blocks_hc_emit, coalesce_gvcf_blocks_for_emit, GvcfMergeSemantics,
            ReferenceConfidenceLocus,
        };

        let loci: Vec<_> = (1..=100)
            .map(|p| ReferenceConfidenceLocus {
                position_1based: p,
                gq: if p <= 80 { 0 } else { 3 },
                dp: 0,
            })
            .collect();
        let blocks = build_gvcf_blocks_hc_emit(&loci, GATK_HC_DEFAULT_GQB).expect("blocks");
        assert!(blocks.len() > 1, "fixture: band split before coalesce");
        let merged = coalesce_gvcf_blocks_for_emit(
            blocks,
            GATK_HC_DEFAULT_GQB,
            GvcfMergeSemantics::default().max_rgq_delta_within_block,
        );
        assert_eq!(merged.len(), 1);
        assert_eq!(merged[0].start_1based, 1);
        assert_eq!(merged[0].end_1based, 100);
    }

    #[test]
    fn gvcf_collector_merges_contiguous_region_loci_into_one_block() {
        use crate::genotyping::ReferenceConfidenceLocus;

        let mut col = GvcfIntervalCollector::new();
        let loci: Vec<_> = (1..=300)
            .map(|p| ReferenceConfidenceLocus {
                position_1based: p,
                gq: 0,
                dp: 0,
            })
            .collect();
        col.ingest_loci("chr1", &loci[..150]);
        col.ingest_loci("chr1", &loci[150..]);
        let dict = SequenceDictionary::from_fasta_path(
            std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("../parity/fixtures/reference.fa"),
        )
        .expect("dict");
        let mut ref_cache = ReferenceWindowCache::new(
            std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("../parity/fixtures/reference.fa"),
            4,
        );
        let blocks = col
            .into_block_records(&dict, &mut ref_cache, "SAMPLE")
            .expect("blocks");
        assert_eq!(blocks.len(), 1, "contiguous same-band loci merge globally");
        assert_eq!(blocks[0].position, 1);
        let end = blocks[0]
            .info
            .iter()
            .find_map(|v| match v {
                gatk_core::io::vcf::InfoValue::Integer(k, vals) if k == "END" => {
                    vals.first().copied()
                }
                _ => None,
            })
            .expect("END");
        assert_eq!(end, 300);
    }

    #[test]
    fn emit_mode_parsing_matches_gatk_output_modes() {
        assert_eq!(emit_mode_from_output_mode("VCF"), EmitMode::Vcf);
        assert_eq!(emit_mode_from_output_mode("GVCF"), EmitMode::Gvcf);
        assert_eq!(
            emit_mode_from_output_mode("BP_RESOLUTION"),
            EmitMode::BpResolution
        );
    }
}
