//! Test-only `call_region` snapshot used by 6R.20+ audits.
use super::*;
use crate::assembly_region_trimmer::{AssemblyRegionTrimResult, TrimVariant};
use std::cell::RefCell;

#[derive(Clone, Debug)]
pub struct AuditEvent {
    pub start: u64,
    pub end: u64,
    pub ref_al: String,
    pub alt_al: String,
    pub is_indel: bool,
}

impl AuditEvent {
    pub fn from_variation(e: &crate::event_map::VariationEvent) -> Self {
        Self {
            start: e.start_1based.get(),
            end: e.end_1based.get(),
            ref_al: e.ref_allele.to_owned(),
            alt_al: e.alt_allele.to_owned(),
            is_indel: e.is_indel(),
        }
    }
}

#[derive(Clone, Debug)]
pub struct AuditTrimVar {
    pub start: u64,
    pub end: u64,
    pub is_indel: bool,
    pub overlaps_active: bool,
}

#[derive(Clone, Debug, Default)]
pub struct CallRegionAuditSnap {
    pub none_at: Option<String>,
    pub region_start: u64,
    pub region_end: u64,
    pub n_haps: usize,
    pub n_alt_haps: usize,
    pub hap_cigars: Vec<String>,
    pub events: Vec<AuditEvent>,
    pub trim_variants: Vec<AuditTrimVar>,
    pub n_trim_overlapping: usize,
    pub trim_variation_present: bool,
    pub trim_variant_start: Option<u64>,
    pub trim_variant_end: Option<u64>,
    pub trim_padded_start: Option<u64>,
    pub trim_padded_end: Option<u64>,
    pub cluster_reads_support: bool,
    pub read_variation_in_active: bool,
    pub disable_optimizations: bool,
    pub n_reads_after_filter: Option<usize>,
    pub has_variation_for_calling: Option<bool>,
    pub n_haps_after_trim: Option<usize>,
    pub n_events_after_trim: Option<usize>,
    pub events_after_trim: Vec<AuditEvent>,
    /// Test-only post-`trim_to` stage dumps (`after_trim_to` / `after_resync`).
    pub post_trim_stages: Vec<String>,
    pub needs_post_trim_resync: Option<bool>,
}

thread_local! {
    static SNAP: RefCell<CallRegionAuditSnap> = RefCell::new(CallRegionAuditSnap::default());
}

pub fn reset(region_start: u64, region_end: u64) {
    SNAP.with(|s| {
        *s.borrow_mut() = CallRegionAuditSnap {
            region_start,
            region_end,
            ..CallRegionAuditSnap::default()
        };
    });
}

pub fn take_call_region_audit() -> CallRegionAuditSnap {
    SNAP.with(|s| s.replace(CallRegionAuditSnap::default()))
}

pub fn note_none(at: &str) {
    SNAP.with(|s| s.borrow_mut().none_at = Some(at.to_string()));
}

pub fn with_mut(f: impl FnOnce(&mut CallRegionAuditSnap)) {
    SNAP.with(|s| f(&mut s.borrow_mut()));
}

pub fn record_after_trim(
    untrimmed: &AssemblyResultSet,
    trim_variants: &[TrimVariant],
    region: &AssemblyRegion,
    trim_result: &AssemblyRegionTrimResult,
    args: &CallRegionArgs,
    apply_bases: &[u8],
    apply_pad: u64,
) {
    let n_overlap = trim_variants
        .iter()
        .filter(|v| v.overlaps_active_region(region))
        .count();
    let cluster_reads_support = args.is_strict_java()
        && crate::read_threading_assembler::region_overlaps_p12_cluster(
            region.start.get(),
            region.end.get(),
        )
        && !crate::read_event_discovery::discover_p12_cluster_coupled_events_from_reads(
            &region.reads,
            apply_bases,
            apply_pad,
            region.start.get(),
            region.end.get(),
            &region.contig,
        )
        .is_empty();
    let read_variation_in_active = args.is_strict_java()
        && active_region_has_read_variation(
            &region.reads,
            apply_bases,
            apply_pad,
            region.start.get(),
            region.end.get(),
            &region.contig,
        );
    with_mut(|a| {
        a.n_haps = untrimmed.haplotypes.len();
        a.n_alt_haps = untrimmed
            .haplotypes
            .iter()
            .filter(|h| !h.is_reference)
            .count();
        a.hap_cigars = untrimmed
            .haplotypes
            .iter()
            .map(|h| {
                format!(
                    "ref={} len={} cigar={} loc={} align={} score={}",
                    h.is_reference,
                    h.bases.len(),
                    h.cigar
                        .as_ref()
                        .map(|c| c.to_gatk_string())
                        .unwrap_or_default(),
                    h.genome_loc
                        .map(|g| format!("{}-{}", g.start_1based(), g.end_1based()))
                        .unwrap_or_else(|| "None".to_string()),
                    h.alignment_start_hap_wrt_ref,
                    h.score
                )
            })
            .collect();
        a.events = untrimmed
            .variation_events()
            .iter()
            .map(AuditEvent::from_variation)
            .collect();
        a.trim_variants = trim_variants
            .iter()
            .map(|v| AuditTrimVar {
                start: v.start,
                end: v.end,
                is_indel: v.is_indel,
                overlaps_active: v.overlaps_active_region(region),
            })
            .collect();
        a.n_trim_overlapping = n_overlap;
        a.trim_variation_present = trim_result.variation_present;
        a.trim_variant_start = trim_result.variant_start;
        a.trim_variant_end = trim_result.variant_end;
        a.trim_padded_start = trim_result.padded_variant_start;
        a.trim_padded_end = trim_result.padded_variant_end;
        a.cluster_reads_support = cluster_reads_support;
        a.read_variation_in_active = read_variation_in_active;
        a.disable_optimizations = args.disable_optimizations;
    });
}

pub fn record_after_read_filter(assembly: &AssemblyResultSet, n_reads: usize) {
    with_mut(|a| {
        a.n_reads_after_filter = Some(n_reads);
        a.has_variation_for_calling = Some(assembly.has_variation_for_calling());
        a.n_haps_after_trim = Some(assembly.haplotypes.len());
        a.n_events_after_trim = Some(assembly.variation_events().len());
        a.events_after_trim = assembly
            .variation_events()
            .iter()
            .map(AuditEvent::from_variation)
            .collect();
    });
}

pub fn note_hap_stage(label: &str, assembly: &AssemblyResultSet) {
    let n_alt = assembly
        .haplotypes
        .iter()
        .filter(|h| !h.is_reference)
        .count();
    let line = format!(
        "{label} n_haps={} n_alt={} n_events={} has_variation_for_calling={} is_variation_present={}",
        assembly.haplotypes.len(),
        n_alt,
        assembly.variation_events().len(),
        assembly.has_variation_for_calling(),
        assembly.is_variation_present(),
    );
    with_mut(|a| a.post_trim_stages.push(line));
}

pub fn note_resync(assembly: &AssemblyResultSet, needs_post_trim_resync: bool) {
    with_mut(|a| a.needs_post_trim_resync = Some(needs_post_trim_resync));
    note_hap_stage("after_resync", assembly);
}

pub fn record_no_variation_for_calling(assembly: &AssemblyResultSet, none_at: &str) {
    with_mut(|a| {
        a.has_variation_for_calling = Some(false);
        a.n_haps_after_trim = Some(assembly.haplotypes.len());
        a.n_events_after_trim = Some(assembly.variation_events().len());
        a.events_after_trim = assembly
            .variation_events()
            .iter()
            .map(AuditEvent::from_variation)
            .collect();
    });
    note_none(none_at);
}
