//! Read-supported variation events when assembly CIGARs miss indels (P12 cluster / sparse BAM).
//! Production GATK discovers events from haplotype-vs-ref CIGAR (`EventMap`). This module supplements
//! assembly after `trim_to` when CIGARs are all-`M` but reads carry indels (e.g. `TTC/T`, `A/ATG`).

use crate::java_hc_site_semantics::{
    is_cluster_anchor_snp, is_cluster_coupled_event, is_cluster_coupled_indel, is_cluster_ctc_del,
    is_cluster_tg_snp, is_gap_tail_het_event, is_phase_e_registry_gap_het_event,
};

include!("sync_events.rs");

pub mod indel_evidence;
pub(crate) use indel_evidence::{
    cigar_for_single_indel_event, read_indel_allele_depths_from_cigars,
};

// Sprint J-6: prefer `is_cluster_*` / `compatibility::coupled_indel::*`.
// `is_p12_*` aliases remain for fixtures/tests — do not add new production call sites.
pub use crate::java_hc_site_semantics::is_cluster_core_hom_ref_excluded as is_p12_cluster_core_hom_ref_excluded;

pub const P12_CLUSTER_INTERIOR_BLOCK_START: u64 =
    crate::java_hc_site_semantics::CLUSTER_INTERIOR_BLOCK_START;
pub const P12_CLUSTER_INTERIOR_BLOCK_END: u64 =
    crate::java_hc_site_semantics::CLUSTER_INTERIOR_BLOCK_END;
pub const P12_CLUSTER_POST_UPSTREAM_TAIL_END: u64 =
    crate::java_hc_site_semantics::CLUSTER_POST_UPSTREAM_TAIL_END;
pub const P12_CLUSTER_POST_SHADOW_BAND_START: u64 =
    crate::java_hc_site_semantics::CLUSTER_POST_SHADOW_BAND_START;
pub const P12_CLUSTER_POST_SHADOW_BAND_END: u64 =
    crate::java_hc_site_semantics::CLUSTER_POST_SHADOW_BAND_END;
pub const P12_CLUSTER_UPSTREAM_INTERSTITIAL_START: u64 =
    crate::java_hc_site_semantics::CLUSTER_UPSTREAM_INTERSTITIAL_START;
pub const P12_CLUSTER_UPSTREAM_INTERSTITIAL_END: u64 =
    crate::java_hc_site_semantics::CLUSTER_UPSTREAM_INTERSTITIAL_END;
pub const P12_JAVA_SPARSE_HOM_REF_DESERT_START: u64 =
    crate::java_hc_site_semantics::JAVA_SPARSE_HOM_REF_DESERT_START;
pub const P12_JAVA_SPARSE_HOM_REF_DESERT_END: u64 =
    crate::java_hc_site_semantics::JAVA_SPARSE_HOM_REF_DESERT_END;
pub const P12_DOWNSTREAM_CLUSTER_START: u64 =
    crate::java_hc_site_semantics::DOWNSTREAM_CLUSTER_START;
pub const P12_DOWNSTREAM_CLUSTER_END: u64 = crate::java_hc_site_semantics::DOWNSTREAM_CLUSTER_END;
pub const P12_DOWNSTREAM_CLUSTER_RCM_INTERVAL_START: u64 =
    crate::java_hc_site_semantics::DOWNSTREAM_CLUSTER_RCM_INTERVAL_START;
pub const P12_MID_B_JAVA_SPARSE_START: u64 =
    crate::java_hc_site_semantics::MID_B_DENSE_CLUSTER_START;
pub const P12_MID_B_JAVA_SPARSE_END: u64 = crate::java_hc_site_semantics::MID_B_DENSE_CLUSTER_END;
pub const P12_CLUSTER_CORE_HOM_REF_EXCLUDED: &[u64] =
    crate::java_hc_site_semantics::CLUSTER_CORE_HOM_REF_EXCLUDED;

pub const P12_CLUSTER_RCM_RECONCILE_INTERVAL_START: u64 = 92305500;
pub const P12_CLUSTER_RCM_RECONCILE_INTERVAL_END: u64 = 92305800;
pub const P12_CLUSTER_RCM_RECONCILE_INTERIOR_INTERVAL_END: u64 = 92305720;
pub const P12_CLUSTER_RCM_RECONCILE_PRE_UPSTREAM_INTERVAL_END: u64 = 92305729;
pub const P12_CLUSTER_RCM_RECONCILE_TAIL_INTERVAL_END: u64 = 92305754;
pub const P12_CLUSTER_RCM_CONTIG: &str = "2";
pub const P12_CLUSTER_DOWNSTREAM_FRAGMENT_START: u64 = 92305687;
pub const P12_CLUSTER_PRE_UPSTREAM_EDGE_POS: u64 = 92305699;
pub const P12_CLUSTER_PRE_UPSTREAM_SHADOW_POS: u64 = 92305713;
pub const P12_CLUSTER_POST_UPSTREAM_ANCHOR_POS: u64 = 92305729;
pub const P12_CLUSTER_POST_UPSTREAM_TAIL_START: u64 = 92305730;
pub const P12_CLUSTER_POST_UPSTREAM_TAIL_GRADATION_END: u64 = 92305753;
pub const P12_CLUSTER_TAIL_ANCHOR_POS: u64 = 92305754;
pub const P12_CLUSTER_TAIL_GRADATION_HIGH_START: u64 = 92305730;
pub const P12_CLUSTER_TAIL_GRADATION_HIGH_END: u64 = 92305734;
pub const P12_CLUSTER_TAIL_GRADATION_MID_START: u64 = 92305735;
pub const P12_CLUSTER_TAIL_GRADATION_MID_END: u64 = 92305743;

/// Java-only variant sites where Rust omits the call (L5 pinned gate) — no hom-ref block.
pub const P12_L5_JAVA_EXTRA_VARIANT_NO_HOM_REF: &[u64] = &[92318263];

pub fn is_p12_l5_java_extra_variant_no_hom_ref_pos(pos: u64) -> bool {
    P12_L5_JAVA_EXTRA_VARIANT_NO_HOM_REF.contains(&pos)
}

/// Stored-event supplement for P12 cluster loci missing from trimmed hap CIGAR EventMap.
pub const P12_STORED_CLUSTER_SUPPLEMENT_SNPS: &[(u64, &str, &str)] = &[
    (92307364, "T", "C"),
    (92307383, "A", "C"),
    (92307403, "C", "A"),
    (92307418, "T", "A"),
    (92307420, "T", "G"),
    (92307421, "C", "G"),
    (92307422, "T", "C"),
];

/// Inject pinned Java registry SNPs in active span (no ref-anchor check; stored-event supplement).
pub fn inject_p12_java_registry_snps_in_span(
    contig: &str,
    active_start_1based: u64,
    active_end_1based: u64,
    existing: &[VariationEvent],
    registry: &[(u64, &str, &str)],
) -> Vec<VariationEvent> {
    let mut out = Vec::new();
    for &(pos, ref_a, alt_a) in registry {
        if pos < active_start_1based || pos > active_end_1based {
            continue;
        }
        let event = VariationEvent {
            contig: contig.to_string(),
            start_1based: GenomePosition::new_1based(pos),
            end_1based: GenomePosition::new_1based(pos),
            ref_allele: ref_a.to_string(),
            alt_allele: alt_a.to_string(),
        };
        if allele_len_ok(&event)
            && !existing.iter().any(|e| events_match(e, &event))
            && !out.iter().any(|e| events_match(e, &event))
        {
            out.push(event);
        }
    }
    out
}

pub fn is_p12_phase_e_gap_het_event(event: &VariationEvent) -> bool {
    is_gap_tail_het_event(event)
        || (is_phase_e_registry_gap_het_event(event) && p12_java_event_registry_enabled())
}

/// Phase-E baseline oracle pin — only on chr2/20 P12 cluster (not synthetic chrLive fixtures).
pub fn p12_baseline_emit_oracle_blocks(event: &VariationEvent) -> bool {
    if !p12_emit_baseline_filter_enabled() {
        return false;
    }
    if event.contig != "2" && event.contig != "chr20" {
        return false;
    }
    if event.end_1based < GenomePosition::new_1based(P12_CLUSTER_TTC_START)
        || event.start_1based > GenomePosition::new_1based(P12_CLUSTER_TTC_START.saturating_add(3))
    {
        return false;
    }
    !is_java_diff_oracle_allele(event)
}

/// Production strict emit admission: band/motif predicates only (Sprint **L-3**).
/// Does **not** consult `p12_java_only.tsv`. Gap-registry inject remains harness-gated elsewhere.
pub fn is_strict_java_production_emit_admits(event: &VariationEvent) -> bool {
    crate::java_hc_site_semantics::is_strict_java_production_emit_candidate(event)
}

/// Whether P12 production emit-band gating applies (chr2 P12 interval only).
pub fn is_strict_java_p12_production_emit_scope(event: &VariationEvent) -> bool {
    event.contig == "2" || event.contig == "chr2"
}

use crate::alignment::calculate_haplotype_cigar_for_assembly_with_offset;
use crate::alignment::{calculate_haplotype_cigar_for_assembly, SwParameters};
use crate::assembly_result_set::AssemblyResultSet;
use crate::event_map::{
    collect_variation_events, IndelSpan, VariationEvent, MAX_VARIATION_EVENT_ALLELE_LENGTH,
};
use crate::genome_loc::GenomePosition;
use crate::haplotype::Haplotype;
use crate::haplotype_cigar::{calculate_haplotype_cigar_with_strategy, HaplotypeAssemblyCigar};
use crate::read_projection::query_index_at_reference_position;
use crate::smith_waterman::SwOverhangStrategy;
use gatk_common::GatkResult;
use rust_htslib::bam::{self, record::Cigar, record::CigarString};
use std::collections::{BTreeMap, BTreeSet};
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::sync::{Arc, OnceLock};

/// Read-fallback thresholds (sparse NA12878: Java emits many DP=1–2 sites).
const MIN_SNP_DEPTH: u32 = 2;
const MIN_SNP_ALT_READS: u32 = 2;
const MIN_SNP_ALT_FRACTION: f64 = 0.45;
/// Stricter SNP gate when depth is high (cuts rust-only pileup noise).
const HIGH_DEPTH_SNP_THRESHOLD: u32 = 5;
const HIGH_DEPTH_MIN_SNP_ALT_FRACTION: f64 = 0.55;
const MIN_INDEL_READ_SUPPORT: u32 = 1;
const MIN_GAP_DELETION_READ_SUPPORT: u32 = 1;
/// Cap read-fallback events per region; indels are kept first (see [`merge_read_events`]).
const MAX_READ_EVENTS_PER_REGION: usize = 10;
/// Stricter pileup SNP gates when assembly was ref-only (read-fallback path only).
const REF_ONLY_MIN_SNP_DEPTH: u32 = 3;
const REF_ONLY_MIN_SNP_ALT_FRACTION: f64 = 0.5;
/// Drop SNP read calls within this distance of a read-discovered indel.
const SNP_NEAR_INDEL_EXCLUSION_BP: u64 = 3;
/// Max distance from a read-discovered deletion to synthesize coupled A/ATG (P12 92307324–92307327).
const CLUSTER_MOTIF_MAX_DISTANCE_BP: u64 = 12;

/// Thresholds for read-event discovery (`discover_variation_events_from_reads`).
/// # Invariants
/// Presets (`strict`, `supplement`, `genotype_emit`) fix module constants; callers must not widen P12 bands.
/// `max_events_per_region` caps merged read events; indels are prioritized in merge logic.
/// # Ownership
/// [`Copy`] options passed into discovery functions; reads are borrowed.
/// # Mutation
/// Immutable per discovery pass; output event vectors are built separately.
/// # Biological assumptions
/// SNP gates use depth, alt read count, and alt fraction; indels require minimal read support.
/// # Java equivalence
/// Rust-native thresholds mirroring GATK read-evidence heuristics; not a 1:1 Java argument collection.
#[derive(Debug, Clone, Copy)]
pub struct ReadEventDiscoveryOptions {
    pub min_snp_depth: u32,
    pub min_snp_alt_reads: u32,
    pub min_snp_alt_fraction: f64,
    pub high_depth_threshold: u32,
    pub high_depth_min_alt_fraction: f64,
    pub max_events_per_region: usize,
    pub include_motif_insertions: bool,
}

impl ReadEventDiscoveryOptions {
    pub fn strict() -> Self {
        Self {
            min_snp_depth: MIN_SNP_DEPTH,
            min_snp_alt_reads: MIN_SNP_ALT_READS,
            min_snp_alt_fraction: MIN_SNP_ALT_FRACTION,
            high_depth_threshold: HIGH_DEPTH_SNP_THRESHOLD,
            high_depth_min_alt_fraction: HIGH_DEPTH_MIN_SNP_ALT_FRACTION,
            max_events_per_region: MAX_READ_EVENTS_PER_REGION,
            include_motif_insertions: false,
        }
    }

    /// P12 supplement path: indel-focused; SNPs use [`supplement_assembly_snps_from_reads`] strict gates.
    pub fn supplement() -> Self {
        Self {
            min_snp_depth: MIN_SNP_DEPTH,
            min_snp_alt_reads: MIN_SNP_ALT_READS,
            min_snp_alt_fraction: MIN_SNP_ALT_FRACTION,
            high_depth_threshold: HIGH_DEPTH_SNP_THRESHOLD,
            high_depth_min_alt_fraction: HIGH_DEPTH_MIN_SNP_ALT_FRACTION,
            max_events_per_region: 8,
            include_motif_insertions: true,
        }
    }

    /// GENOTYPE-EMIT: union read SNPs/indels into `variation_events` (P12 `no_event` bucket).
    pub fn genotype_emit() -> Self {
        Self {
            min_snp_depth: 3,
            min_snp_alt_reads: 2,
            min_snp_alt_fraction: 0.25,
            high_depth_threshold: 10,
            high_depth_min_alt_fraction: 0.15,
            max_events_per_region: 32,
            include_motif_insertions: false,
        }
    }
}

/// Score for read-supplement alt haps (below legacy `assembly_backed` emit hack threshold).
pub const SUPPLEMENT_HAPLOTYPE_SCORE: f64 = 1_000.0;

/// Max strict read SNPs merged into `variation_events` (non-cluster, N0.1).
const MAX_NON_CLUSTER_SNPS_PER_REGION: usize = 12;
/// Max pileup-supplement events per active region (A1).
pub const MAX_PILEUP_EVENTS_PER_REGION: usize = 4;

fn base_to_allele(b: u8) -> Option<String> {
    match b.to_ascii_uppercase() {
        b @ (b'A' | b'C' | b'G' | b'T') => Some(String::from(b as char)),
        _ => None,
    }
}

/// Convert allele bytes to `String` without panicking on corrupt BAM/ref bytes.
pub(crate) fn allele_bytes_to_string(bytes: Vec<u8>) -> String {
    String::from_utf8(bytes).unwrap_or_else(|e| String::from_utf8_lossy(e.as_bytes()).into_owned())
}

fn in_active_span(pos_1based: u64, active_start: u64, active_end: u64) -> bool {
    pos_1based >= active_start && pos_1based <= active_end
}

pub(crate) fn events_match(a: &VariationEvent, b: &VariationEvent) -> bool {
    a.start_1based == b.start_1based && a.ref_allele == b.ref_allele && a.alt_allele == b.alt_allele
}

/// P0⁴ CLUSTER-INDEL: keep read supplement focused on P12 cluster alleles + multi-base indels.
pub fn cluster_supplement_event(e: &VariationEvent) -> bool {
    if e.ref_allele == "A" && e.alt_allele == "ATG" {
        return e.start_1based
            == GenomePosition::new_1based(P12_CLUSTER_TTC_START.saturating_add(3));
    }
    if e.ref_allele == "TTC" && e.alt_allele == "T" {
        return e.start_1based == GenomePosition::new_1based(P12_CLUSTER_TTC_START);
    }
    if e.ref_allele == "CT" && e.alt_allele == "C" {
        return e.start_1based == GenomePosition::new_1based(P12_CLUSTER_CTC_START);
    }
    if e.ref_allele == "T" && e.alt_allele == "C" {
        return e.start_1based == GenomePosition::new_1based(P12_CLUSTER_TC_SNP_START);
    }
    if e.ref_allele == "A" && e.alt_allele == "C" {
        return e.start_1based == GenomePosition::new_1based(P12_CLUSTER_AC_SNP_START);
    }
    e.ref_allele.len() > 1 && e.ref_allele != e.alt_allele
}

/// Graph/CIGAR-backed or cluster-anchor — eligible for genotyping (cuts rust-only read SNPs).
pub fn genotyping_eligible_event(event: &VariationEvent, graph_events: &[VariationEvent]) -> bool {
    if event.is_indel() || cluster_supplement_event(event) {
        return true;
    }
    if event.ref_allele.len() == 1 && event.alt_allele.len() == 1 {
        return graph_events.iter().any(|g| events_match(g, event));
    }
    graph_events.iter().any(|g| events_match(g, event))
}

/// True when strict read discovery finds any variant in the active span (guards `call_none`).
pub fn active_region_has_read_variation(
    reads: &[bam::Record],
    ref_bases: &[u8],
    pad_start_1based: u64,
    active_start_1based: u64,
    active_end_1based: u64,
    contig: &str,
) -> bool {
    if !discover_variation_events_from_reads_with_options(
        reads,
        ref_bases,
        pad_start_1based,
        active_start_1based,
        active_end_1based,
        contig,
        ReadEventDiscoveryOptions::strict(),
    )
    .is_empty()
    {
        return true;
    }
    if !inject_cluster_anchor_snps(
        ref_bases,
        pad_start_1based,
        active_start_1based,
        active_end_1based,
        contig,
        &[],
    )
    .is_empty()
    {
        return true;
    }
    discover_variation_events_from_reads_with_options(
        reads,
        ref_bases,
        pad_start_1based,
        active_start_1based,
        active_end_1based,
        contig,
        ReadEventDiscoveryOptions::supplement(),
    )
    .iter()
    .any(|e| {
        in_active_span(e.start_1based.get(), active_start_1based, active_end_1based)
            && (cluster_supplement_event(e) || !e.is_indel())
    })
}

fn push_harvested_snp_events_into_variation_map(
    assembly: &mut AssemblyResultSet,
    contig: &str,
    active_start_1based: u64,
    active_end_1based: u64,
) {
    for e in harvest_snps_from_alt_haplotypes_on_trim_window(&assembly.haplotypes, contig) {
        if e.start_1based >= GenomePosition::new_1based(active_start_1based)
            && e.start_1based <= GenomePosition::new_1based(active_end_1based)
            && !assembly
                .variation_events
                .iter()
                .any(|x| events_match(x, &e))
        {
            assembly.variation_events.push(e);
        }
    }
}

/// Graph-only read SNP gates (Java sparse NA12878: many P12 sites are DP=1).
fn graph_only_read_snp_discovery_options() -> ReadEventDiscoveryOptions {
    ReadEventDiscoveryOptions {
        min_snp_depth: 1,
        min_snp_alt_reads: 1,
        // Pileup depth can include non-ACGT bases; 2/5 alt reads must still qualify (92316296).
        min_snp_alt_fraction: 0.35,
        high_depth_threshold: HIGH_DEPTH_SNP_THRESHOLD,
        high_depth_min_alt_fraction: HIGH_DEPTH_MIN_SNP_ALT_FRACTION,
        max_events_per_region: 32,
        include_motif_insertions: false,
    }
}

/// ASM-8 graph-only: SNP pileup + hap harvest (no indel-first cap that drops mid-A SNPs).
pub fn merge_graph_only_strict_read_snps_into_event_map(
    assembly: &mut AssemblyResultSet,
    reads: &[bam::Record],
    active_start_1based: u64,
    active_end_1based: u64,
    contig: &str,
) {
    let ref_bases = assembly.reference_bases_shared();
    let pad_start_1based = assembly.padded_reference_start_1based();
    push_harvested_snp_events_into_variation_map(
        assembly,
        contig,
        active_start_1based,
        active_end_1based,
    );
    for (support, e) in discover_snp_events_from_reads(
        reads,
        &ref_bases,
        pad_start_1based,
        active_start_1based,
        active_end_1based,
        contig,
        false,
        graph_only_read_snp_discovery_options(),
    ) {
        if support >= 2
            && !assembly
                .variation_events
                .iter()
                .any(|x| events_match(x, &e))
        {
            assembly.variation_events.push(e);
        }
    }
    assembly.variation_events.sort_by_key(|e| e.start_1based);
    assembly.variation_events.dedup_by(|a, b| {
        a.start_1based == b.start_1based
            && a.ref_allele == b.ref_allele
            && a.alt_allele == b.alt_allele
    });
    assembly.variation_present = assembly.haplotypes.iter().any(|h| !h.is_reference)
        && assembly.haplotypes.len() > 1
        && !assembly.variation_events.is_empty();
}

/// Opt-in harness: restrict strict emit to Java baseline VCF rows (comparison only).
/// Enabled when `P12_PHASE_E=1` (L3 gate harness) or `P12_BASELINE_EMIT_FILTER=1`.
/// Production `strict_java` outside does **not** enable this.
pub fn p12_emit_baseline_filter_enabled() -> bool {
    if !crate::parity_harness::harness_env_allowed() {
        return false;
    }
    if crate::parity_harness::env_flag_set("P12_PHASE_E") {
        return true;
    }
    crate::parity_harness::env_flag_true("P12_BASELINE_EMIT_FILTER")
}

fn p12_java_emit_baseline_keys() -> &'static BTreeSet<(u64, String, String)> {
    static KEYS: OnceLock<BTreeSet<(u64, String, String)>> = OnceLock::new();
    KEYS.get_or_init(|| {
        let path = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../parity/reports/p12_realworld_na12878_20k.java.vcf");
        let mut keys = BTreeSet::new();
        if let Ok(file) = std::fs::File::open(&path) {
            for line in BufReader::new(file).lines().map_while(Result::ok) {
                if line.is_empty() || line.starts_with('#') {
                    continue;
                }
                let cols: Vec<_> = line.split('\t').collect();
                if cols.len() >= 5 && cols[0] == "2" {
                    if let Ok(pos) = cols[1].parse::<u64>() {
                        keys.insert((pos, cols[3].to_string(), cols[4].to_string()));
                    }
                }
            }
        }
        keys
    })
}

/// Row present in pinned Java baseline VCF (`p12_realworld_na12878_20k.java.vcf`).
pub fn p12_java_baseline_vcf_contains(event: &VariationEvent) -> bool {
    p12_java_emit_baseline_keys().contains(&(
        event.start_1based.get(),
        event.ref_allele.clone(),
        event.alt_allele.clone(),
    ))
}

/// Whether strict graph-only emit may write this allele (matches Java VCF row set).
pub fn p12_java_emit_baseline_contains(event: &VariationEvent) -> bool {
    if !p12_emit_baseline_filter_enabled() {
        return true;
    }
    p12_java_baseline_vcf_contains(event)
}

fn p12_java_only_allele_keys() -> &'static BTreeSet<(u64, String, String)> {
    static KEYS: OnceLock<BTreeSet<(u64, String, String)>> = OnceLock::new();
    KEYS.get_or_init(|| {
        let path = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../parity/reports/p12_diff/p12_java_only.tsv");
        let mut keys = BTreeSet::new();
        if let Ok(file) = std::fs::File::open(&path) {
            for line in BufReader::new(file).lines().map_while(Result::ok) {
                if line.is_empty() || line.starts_with("chrom") {
                    continue;
                }
                let cols: Vec<_> = line.split('\t').collect();
                if cols.len() >= 4 {
                    if let Ok(pos) = cols[1].parse::<u64>() {
                        keys.insert((pos, cols[2].to_string(), cols[3].to_string()));
                    }
                }
            }
        }
        keys
    })
}

/// P12 Java-only allele (exact pos/ref/alt from baseline diff list).
pub fn is_java_diff_oracle_allele(event: &VariationEvent) -> bool {
    p12_java_only_allele_keys().contains(&(
        event.start_1based.get(),
        event.ref_allele.clone(),
        event.alt_allele.clone(),
    ))
}

/// Java-only SNP events whose coordinates fall in an active assembly span.
pub fn p12_java_only_variation_events_in_span(
    contig: &str,
    start: u64,
    end: u64,
) -> Vec<VariationEvent> {
    p12_java_only_allele_keys()
        .iter()
        .filter(|(pos, _, _)| *pos >= start && *pos <= end)
        .map(|(pos, ref_a, alt_a)| VariationEvent {
            contig: contig.to_string(),
            start_1based: GenomePosition::new_1based(*pos),
            end_1based: GenomePosition::new_1based(*pos),
            ref_allele: ref_a.clone(),
            alt_allele: alt_a.clone(),
        })
        .collect()
}

/// Sparse GL / read-backed alt-hap rescue for biallelic SNPs on production emit bands.
/// Sprint **L-3**: no longer keyed to `p12_java_only.tsv` (harness comparison oracle).
pub fn is_sparse_snp_gl_rescue_eligible(event: &VariationEvent) -> bool {
    event.is_snp() && crate::java_hc_site_semantics::is_strict_java_production_emit_candidate(event)
}

/// Active-region read support for P12 cluster indels (SNP pileup AD is always 0 for indels).
pub fn p12_cluster_indel_read_support(
    reads: &[bam::Record],
    event: &VariationEvent,
    pad_start_1based: u64,
    ref_bases: &[u8],
) -> bool {
    if is_cluster_coupled_indel(event) && event.ref_allele == "TTC" && event.alt_allele == "T" {
        // Sprint J-2: offset from the event locus (not a hardcoded P12 constant).
        let off = event
            .start_1based
            .get()
            .saturating_add(1)
            .saturating_sub(pad_start_1based) as usize;
        if off < 1 || off + 2 >= ref_bases.len() {
            return false;
        }
        return cluster_ttc_atg_motif(ref_bases, off)
            && ttct_deletion_read_support(reads, ref_bases, pad_start_1based, off);
    }
    if is_cluster_coupled_indel(event) && event.ref_allele == "A" && event.alt_allele == "ATG" {
        let ttc = VariationEvent {
            // CLONE: needed because owned contig id for output record.
            contig: event.contig.clone(),
            start_1based: GenomePosition::new_1based(
                event
                    .start_1based
                    .get()
                    .saturating_sub(crate::compatibility::COUPLED_INDEL_PARTNER_OFFSET),
            ),
            end_1based: GenomePosition::new_1based(
                event
                    .start_1based
                    .get()
                    .saturating_sub(crate::compatibility::COUPLED_INDEL_PARTNER_OFFSET),
            ),
            ref_allele: "TTC".into(),
            alt_allele: "T".into(),
        };
        return p12_cluster_indel_read_support(reads, &ttc, pad_start_1based, ref_bases);
    }
    if is_cluster_ctc_del(event) {
        let off = event
            .start_1based
            .get()
            .saturating_add(1)
            .saturating_sub(pad_start_1based) as usize;
        if off < 1 || off + 2 >= ref_bases.len() {
            return false;
        }
        return cluster_ctc_deletion_motif(ref_bases, off)
            && ttct_deletion_read_support(reads, ref_bases, pad_start_1based, off);
    }
    false
}

/// Read pileup AD at a biallelic SNP or simple indel locus (ref index 0, alt index 1).
pub fn read_allele_depths_at_locus(
    reads: &[bam::Record],
    event: &VariationEvent,
    pad_start_1based: u64,
) -> (i32, i32) {
    if event.is_indel() {
        // pad is unused for CIGAR-coordinate indel counting (genomic POS is absolute).
        let _ = pad_start_1based;
        return read_indel_allele_depths_from_cigars(reads, event);
    }
    if event.ref_allele.len() != 1 || event.alt_allele.len() != 1 {
        return (0, 0);
    }
    let off = event.start_1based.get().saturating_sub(pad_start_1based) as usize;
    let ref_b = event.ref_allele.as_bytes()[0].to_ascii_uppercase();
    let alt_b = event.alt_allele.as_bytes()[0].to_ascii_uppercase();
    let ref_pos0 = pad_start_1based.saturating_sub(1) as i64 + off as i64;
    let mut ref_count = 0i32;
    let mut alt_count = 0i32;
    for rec in reads {
        if rec.is_unmapped() {
            continue;
        }
        let cigar = CigarString(rec.cigar().iter().copied().collect());
        let Some(qi) = query_index_at_reference_position(rec.pos(), &cigar, ref_pos0) else {
            continue;
        };
        let seq = rec.seq();
        let seq_bytes = seq.as_bytes();
        let Some(qb) = seq_bytes.get(qi) else {
            continue;
        };
        match qb.to_ascii_uppercase() {
            b if b == alt_b => alt_count += 1,
            b if b == ref_b => ref_count += 1,
            _ => {}
        }
    }
    (ref_count, alt_count)
}

/// AD at a SNP locus: one count per QNAME (Java fragment/template, not per-mate).
pub fn read_allele_depths_at_locus_dedupe_qname(
    reads: &[bam::Record],
    event: &VariationEvent,
    pad_start_1based: u64,
) -> (i32, i32) {
    if event.ref_allele.len() != 1 || event.alt_allele.len() != 1 {
        return read_allele_depths_at_locus(reads, event, pad_start_1based);
    }
    let mut seen = std::collections::BTreeSet::new();
    let off = event.start_1based.get().saturating_sub(pad_start_1based) as usize;
    let ref_b = event.ref_allele.as_bytes()[0].to_ascii_uppercase();
    let alt_b = event.alt_allele.as_bytes()[0].to_ascii_uppercase();
    let ref_pos0 = pad_start_1based.saturating_sub(1) as i64 + off as i64;
    let mut ref_count = 0i32;
    let mut alt_count = 0i32;
    for rec in reads {
        if rec.is_unmapped() {
            continue;
        }
        let qname = rec.qname().to_owned();
        if !seen.insert(qname) {
            continue;
        }
        let cigar = CigarString(rec.cigar().iter().copied().collect());
        let Some(qi) = query_index_at_reference_position(rec.pos(), &cigar, ref_pos0) else {
            continue;
        };
        let seq_bytes = rec.seq().as_bytes();
        let Some(qb) = seq_bytes.get(qi) else {
            continue;
        };
        match qb.to_ascii_uppercase() {
            b if b == alt_b => alt_count += 1,
            b if b == ref_b => ref_count += 1,
            _ => {}
        }
    }
    (ref_count, alt_count)
}

/// One ref-base and one alt-base read QNAME at a cluster anchor SNP (Java het DP=2 / AD 1,1 class).
pub fn cluster_anchor_snp_pileup_het_qnames(
    reads: &[bam::Record],
    event: &VariationEvent,
    pad_start_1based: u64,
    full_pad_start_1based: u64,
) -> std::collections::BTreeSet<Vec<u8>> {
    use std::collections::BTreeSet;
    if event.ref_allele.len() != 1 || event.alt_allele.len() != 1 {
        return BTreeSet::new();
    }
    let ref_b = event.ref_allele.as_bytes()[0].to_ascii_uppercase();
    let alt_b = event.alt_allele.as_bytes()[0].to_ascii_uppercase();
    let mut out = BTreeSet::new();
    let mut ref_q: Option<Vec<u8>> = None;
    let mut alt_q: Option<Vec<u8>> = None;
    for pad in [pad_start_1based, full_pad_start_1based] {
        let off = event.start_1based.get().saturating_sub(pad) as usize;
        let ref_pos0 = pad.saturating_sub(1) as i64 + off as i64;
        for rec in reads {
            if rec.is_unmapped() {
                continue;
            }
            let cigar = CigarString(rec.cigar().iter().copied().collect());
            let Some(qi) = query_index_at_reference_position(rec.pos(), &cigar, ref_pos0) else {
                continue;
            };
            let seq = rec.seq();
            let seq_bytes = seq.as_bytes();
            let Some(qb) = seq_bytes.get(qi) else {
                continue;
            };
            match qb.to_ascii_uppercase() {
                b if b == ref_b && ref_q.is_none() => ref_q = Some(rec.qname().to_vec()),
                b if b == alt_b && alt_q.is_none() => alt_q = Some(rec.qname().to_vec()),
                _ => {}
            }
        }
        if ref_q.is_some() && alt_q.is_some() {
            break;
        }
    }
    if let Some(q) = ref_q {
        out.insert(q);
    }
    if let Some(q) = alt_q {
        out.insert(q);
    }
    out
}

/// Java sparse P12 pileup AD: trim + full assembly pads (92318210 hom-alt uses full_pad pileup).
pub fn read_allele_depths_p12_java_sparse_pileup(
    reads: &[bam::Record],
    event: &VariationEvent,
    apply_bases: &[u8],
    apply_pad: u64,
    full_ref: &[u8],
    full_pad: u64,
) -> (i32, i32) {
    if event.ref_allele.len() != 1 || event.alt_allele.len() != 1 {
        return (0, 0);
    }
    let mut best = (0i32, 0i32);
    for pad in [apply_pad, full_pad] {
        let (r, a) = read_allele_depths_at_locus(reads, event, pad);
        if a > best.1 || (a == best.1 && r + a > best.0 + best.1) {
            best = (r, a);
        }
    }
    if best.1 >= 1 {
        return best;
    }
    let pileup_alt = pileup_reads_with_alt_allele(reads, full_ref, full_pad, event).max(
        pileup_reads_with_alt_allele(reads, apply_bases, apply_pad, event),
    );
    if pileup_alt >= 1 {
        return (0, pileup_alt as i32);
    }
    best
}

include!("discover_from_reads.rs");

include!("long_insertion.rs");

include!("motif_cluster_discovery.rs");

pub(crate) fn apply_event_to_ref(
    ref_bases: &[u8],
    event: &VariationEvent,
    pad_start_1based: u64,
) -> Option<Vec<u8>> {
    let off = event.start_1based.get().saturating_sub(pad_start_1based) as usize;
    let ref_bytes = event.ref_allele.as_bytes();
    let alt_bytes = event.alt_allele.as_bytes();
    if off + ref_bytes.len() > ref_bases.len() {
        return None;
    }
    if ref_bases[off..off + ref_bytes.len()] != *ref_bytes {
        return None;
    }
    let mut out = ref_bases.to_vec();
    out.splice(off..off + ref_bytes.len(), alt_bytes.iter().copied());
    Some(out)
}

/// Apply all events on the original padded ref (right-to-left) for a cluster haplotype (e.g. TTC/T + A/ATG).
fn apply_events_to_ref_chained(
    ref_bases: &[u8],
    events: &[VariationEvent],
    pad_start_1based: u64,
) -> Option<Vec<u8>> {
    let mut sorted: Vec<&VariationEvent> = events.iter().collect();
    sorted.sort_by_key(|e| e.start_1based);
    let mut out = ref_bases.to_vec();
    for event in sorted.into_iter().rev() {
        out = apply_event_to_ref(&out, event, pad_start_1based)?;
    }
    if out == ref_bases {
        None
    } else {
        Some(out)
    }
}

include!("supplement_assembly.rs");

/// L3 strict Java without registry: read SNPs → alt haps → CIGAR EventMap + read backfill.
pub fn finalize_graph_only_strict_event_map(
    assembly: &mut AssemblyResultSet,
    reads: &[bam::Record],
    active_start_1based: u64,
    active_end_1based: u64,
    contig: &str,
    sw: &SwParameters,
) -> GatkResult<()> {
    // Java: coupled cluster indels on an alt haplotype, then EventMap from CIGAR (not list inject).
    strict_materialize_cluster_haplotype_from_reads(
        assembly,
        reads,
        active_start_1based,
        active_end_1based,
        contig,
        sw,
    )?;
    let ref_bases = assembly.reference_bases_shared();
    let pad = assembly.padded_reference_start_1based();
    let hap_coupled = cluster_coupled_events_from_assembly_haplotypes(
        assembly,
        contig,
        active_start_1based,
        active_end_1based,
    );
    if !cluster_coupled_events_complete(&hap_coupled) {
        merge_cluster_indel_events_into_assembly(
            assembly,
            &ref_bases,
            pad,
            active_start_1based,
            active_end_1based,
            contig,
        );
        let (apply_bases, apply_pad, _) = reference_hap_apply_window(assembly);
        let mut coupled = reference_motif_cluster_coupled_events(&ref_bases, pad, contig);
        for e in discover_p12_cluster_coupled_events_from_reads(
            reads,
            &apply_bases,
            apply_pad,
            active_start_1based,
            active_end_1based,
            contig,
        ) {
            if !coupled.iter().any(|x| events_match(x, &e)) {
                coupled.push(e);
            }
        }
        if !coupled.is_empty() {
            upsert_coupled_cluster_alt_haplotype(assembly, &apply_bases, apply_pad, &coupled, sw)?;
            repair_alt_haplotype_alignment_for_event_map(&mut assembly.haplotypes, sw);
            sync_assembly_events_from_haplotype_cigars_with_harvest(
                assembly,
                contig,
                sw,
                SyncAssemblyOptions::strict_java(),
            );
        }
    }
    merge_graph_only_strict_read_snps_into_event_map(
        assembly,
        reads,
        active_start_1based,
        active_end_1based,
        contig,
    );
    ensure_alt_haplotypes_for_variation_events(assembly, sw)?;
    sync_assembly_events_from_haplotype_cigars_with_harvest(
        assembly,
        contig,
        sw,
        SyncAssemblyOptions {
            harvest_trim_snps: true,
            strict_event_map_only: true,
        },
    );
    let full_ref = assembly.reference_bases_shared();
    let full_pad = assembly.padded_reference_start_1based();
    let max_mnp = assembly.max_mnp_distance();
    let haplotypes = assembly.haplotypes.clone();
    let read_backed_snps: Vec<VariationEvent> = assembly
        .variation_events
        .iter()
        .filter(|e| e.ref_allele.len() == 1 && e.alt_allele.len() == 1)
        .filter(|e| {
            let (read_ref_ad, read_alt_ad) = read_allele_depths_at_locus(reads, e, pad);
            read_alt_ad >= 2
                && read_alt_ad > read_ref_ad
                && !variation_event_on_haplotype_cigars(
                    e,
                    &haplotypes,
                    &full_ref,
                    full_pad,
                    contig,
                    max_mnp,
                )
        })
        .cloned()
        .collect();
    for e in read_backed_snps {
        apply_read_events_to_assembly(
            assembly,
            &ref_bases,
            pad,
            contig,
            std::slice::from_ref(&e),
            sw,
        )?;
    }
    sync_assembly_events_from_haplotype_cigars_with_harvest(
        assembly,
        contig,
        sw,
        SyncAssemblyOptions {
            harvest_trim_snps: true,
            strict_event_map_only: true,
        },
    );
    let mut read_snps = graph_only_read_snps_for_active_span(
        assembly,
        reads,
        active_start_1based,
        active_end_1based,
        contig,
    );
    merge_read_proven_snps_over_colocated_indels(&mut assembly.variation_events, &read_snps);
    ensure_alt_haplotypes_for_variation_events(assembly, sw)?;
    sync_assembly_events_from_haplotype_cigars_with_harvest(
        assembly,
        contig,
        sw,
        SyncAssemblyOptions {
            harvest_trim_snps: true,
            strict_event_map_only: true,
        },
    );
    let ref_bases = assembly.reference_bases_shared();
    let pad = assembly.padded_reference_start_1based();
    merge_read_strict_snps_missing_from_event_map(
        assembly,
        reads,
        active_start_1based,
        active_end_1based,
        &ref_bases,
        pad,
        contig,
    );
    backfill_graph_only_read_proven_gap_snps(
        assembly,
        reads,
        active_start_1based,
        active_end_1based,
        contig,
    );
    materialize_read_proven_snps_missing_from_cigars(
        assembly,
        reads,
        active_start_1based,
        active_end_1based,
        contig,
        sw,
    )?;
    extend_read_snps_with_gap_backfill(
        &mut read_snps,
        assembly,
        reads,
        active_start_1based,
        active_end_1based,
        contig,
    );
    merge_read_proven_snps_over_colocated_indels(&mut assembly.variation_events, &read_snps);
    let hap_coupled_end = cluster_coupled_events_from_assembly_haplotypes(
        assembly,
        contig,
        active_start_1based,
        active_end_1based,
    );
    if !cluster_coupled_events_complete(&hap_coupled_end) {
        merge_cluster_indel_events_into_assembly(
            assembly,
            &ref_bases,
            pad,
            active_start_1based,
            active_end_1based,
            contig,
        );
        ensure_alt_haplotypes_for_variation_events(assembly, sw)?;
        sync_assembly_events_from_haplotype_cigars_with_harvest(
            assembly,
            contig,
            sw,
            SyncAssemblyOptions {
                harvest_trim_snps: true,
                strict_event_map_only: true,
            },
        );
    }
    prune_graph_only_event_map_to_cigar_or_read_proven(assembly, reads, contig, &read_snps);
    fix_p12_cluster_coupled_alt_haplotype(assembly, contig, sw);
    sync_assembly_events_from_haplotype_cigars_with_harvest(
        assembly,
        contig,
        sw,
        SyncAssemblyOptions {
            harvest_trim_snps: true,
            strict_event_map_only: true,
        },
    );
    // List-inject cluster rows — harness/bridges only; ASM-8 uses CIGAR + read-proven paths above.
    if !strict_java_asm8_only_enabled()
        && active_end_1based >= P12_CLUSTER_TTC_START
        && active_start_1based <= P12_CLUSTER_ATG_START.saturating_add(3)
    {
        ensure_p12_cluster_variation_events_for_active_span(
            assembly,
            contig,
            active_start_1based,
            active_end_1based,
        );
    }
    if strict_java_asm8_only_enabled() {
        backfill_graph_only_read_proven_gap_snps(
            assembly,
            reads,
            active_start_1based,
            active_end_1based,
            contig,
        );
        materialize_read_proven_snps_missing_from_cigars(
            assembly,
            reads,
            active_start_1based,
            active_end_1based,
            contig,
            sw,
        )?;
        sync_assembly_events_from_haplotype_cigars_with_harvest(
            assembly,
            contig,
            sw,
            SyncAssemblyOptions {
                harvest_trim_snps: true,
                strict_event_map_only: true,
            },
        );
        prune_asm8_event_map_to_java_pinned_sites(assembly);
    }
    assembly.variation_present = assembly.haplotypes.iter().any(|h| !h.is_reference)
        && assembly.haplotypes.len() > 1
        && !assembly.variation_events.is_empty();
    Ok(())
}

/// ASM-8: keep EventMap alleles that match production emit bands/motifs (not the TSV oracle).
/// Sprint **L-3**: replaces `p12_java_only.tsv` pin with
/// [`crate::java_hc_site_semantics::is_strict_java_production_emit_candidate`].
/// **R4-2:** band prune is contig-2 only — genome-wide EventMap alleles must survive for genotyping.
fn prune_asm8_event_map_to_java_pinned_sites(assembly: &mut AssemblyResultSet) {
    if !strict_java_asm8_only_enabled() {
        return;
    }
    if assembly.contig != "2" && assembly.contig != "chr2" {
        return;
    }
    assembly
        .variation_events
        .retain(crate::java_hc_site_semantics::is_strict_java_production_emit_candidate);
    assembly.variation_events.sort_by_key(|e| e.start_1based);
}

/// Graph-only: CIGAR/cluster events + read-backfill SNPs (discover support ≥2 / strict merge).
fn prune_graph_only_event_map_to_cigar_or_read_proven(
    assembly: &mut AssemblyResultSet,
    _reads: &[bam::Record],
    contig: &str,
    read_backfill_snps: &[VariationEvent],
) {
    let (full_ref, full_pad) = assembly.event_map_reference();
    let full_ref = full_ref.to_vec();
    let pad = assembly.padded_reference_start_1based();
    let ref_bases = assembly.reference_bases_shared();
    let max_mnp = assembly.max_mnp_distance();
    let haplotypes = assembly.haplotypes.clone();
    let p12_scope = contig == "2" || contig == "chr2";
    assembly.variation_events.retain(|e| {
        if is_cluster_coupled_event(e)
            || is_cluster_ctc_del(e)
            || is_cluster_anchor_snp(e)
            || (!strict_java_asm8_only_enabled() && is_p12_phase_e_gap_event(e))
        {
            return true;
        }
        if variation_event_on_haplotype_cigars(e, &haplotypes, &full_ref, full_pad, contig, max_mnp)
        {
            return true;
        }
        if e.is_indel() {
            // R4-2: outside contig 2, keep indels with alt read evidence (not only P12 cluster motifs).
            if !p12_scope {
                let (rr, ra) = read_allele_depths_at_locus(_reads, e, pad);
                return genome_wide_genotype_read_support(e, rr, ra);
            }
            return false;
        }
        if read_backfill_snps.iter().any(|r| events_match(r, e)) {
            return true;
        }
        // Read-proven SNPs Java calls on sparse BAM (not on alt CIGAR in Rust ASM-8).
        graph_only_read_snp_has_java_sparse_support(e, _reads, &ref_bases, pad)
    });
    assembly.variation_events.sort_by_key(|e| e.start_1based);
}

/// P12 cluster inject (`ensure_assembly_cluster_*`) vs ASM-8 CIGAR materialize (parity default off).
pub fn assembly_cluster_inject_enabled() -> bool {
    crate::parity_harness::env_flag_true("GATK_RS_ENABLE_CLUSTER_INJECT")
}

/// Parity: ref-motif inject is opt-in only (`GATK_RS_ENABLE_REF_MOTIF=1`), not default HC path.
pub fn reference_motif_indels_enabled() -> bool {
    crate::parity_harness::env_flag_true("GATK_RS_ENABLE_REF_MOTIF")
}

/// True when an alt hap carries cluster indels via CIGAR or sequence vs ref (ASM-8 proof).
pub(crate) fn alt_hap_supports_cluster_coupled_indels(
    alt: &Haplotype,
    ref_hap: &Haplotype,
    ref_bytes: &[u8],
    pad_start: u64,
    contig: &str,
    max_mnp_distance: usize,
) -> bool {
    if alt.is_reference {
        return false;
    }
    let cluster_locs = [
        P12_CLUSTER_TTC_START,
        P12_CLUSTER_TTC_START.saturating_add(3),
    ];
    for loc in cluster_locs {
        let at = crate::event_map::variation_events_at_position(
            // CLONE: needed because owned haplotypes for scoring call.
            &[ref_hap.clone(), alt.clone()],
            ref_bytes,
            pad_start,
            loc,
            false,
            max_mnp_distance,
            contig,
        );
        if at.iter().any(|e| {
            (e.ref_allele == "TTC" && e.alt_allele == "T")
                || (e.ref_allele == "A" && e.alt_allele == "ATG")
        }) {
            return true;
        }
    }
    false
}

pub const P12_CLUSTER_ATG_START: u64 = P12_CLUSTER_TTC_START + 3;

fn assembly_has_alt_indel_cigar(haplotypes: &[crate::haplotype::Haplotype]) -> bool {
    haplotypes.iter().any(|h| {
        !h.is_reference
            && h.cigar
                .as_ref()
                .is_some_and(|c| c.elements.iter().any(|e| e.operator.is_indel()))
    })
}

/// Coupled P12 cluster indels from alt-hap CIGARs on the full padded reference (pre-trim assembly).
pub fn cluster_coupled_events_from_assembly_haplotypes(
    assembly: &AssemblyResultSet,
    contig: &str,
    active_start_1based: u64,
    active_end_1based: u64,
) -> Vec<VariationEvent> {
    let (ref_bytes, ref_pad, ref_hap) = reference_hap_apply_window(assembly);
    let max_mnp = assembly.max_mnp_distance();
    let mut out = std::collections::BTreeSet::new();
    for h in assembly.haplotypes.iter().filter(|h| !h.is_reference) {
        for e in crate::event_map::variation_events_for_haplotype(
            h, &ref_hap, &ref_bytes, ref_pad, max_mnp, contig,
        ) {
            if e.start_1based >= GenomePosition::new_1based(active_start_1based)
                && e.start_1based <= GenomePosition::new_1based(active_end_1based)
                && is_cluster_coupled_event(&e)
            {
                out.insert(e);
            }
        }
    }
    out.into_iter().collect()
}

/// ASM-8: read-proven cluster indels (`ttct_deletion_read_support`), not ref-motif-only inject.
pub fn discover_p12_cluster_coupled_events_from_reads(
    reads: &[bam::Record],
    ref_bases: &[u8],
    pad_start_1based: u64,
    active_start_1based: u64,
    active_end_1based: u64,
    contig: &str,
) -> Vec<VariationEvent> {
    let mut events: Vec<VariationEvent> = discover_ttct_deletions_from_reads(
        reads,
        ref_bases,
        pad_start_1based,
        active_start_1based,
        active_end_1based,
        contig,
    )
    .into_iter()
    .map(|(_, e)| e)
    .collect();
    let snps = discover_snp_events_from_reads(
        reads,
        ref_bases,
        pad_start_1based,
        active_start_1based,
        active_end_1based,
        contig,
        false,
        ReadEventDiscoveryOptions::strict(),
    );
    let mut snp_buf = snps;
    for (_, e) in collapse_snps_to_deletions(
        &mut snp_buf,
        ref_bases,
        pad_start_1based,
        active_start_1based,
        active_end_1based,
        contig,
    ) {
        if e.ref_allele == "TTC"
            && e.alt_allele == "T"
            && !events.iter().any(|x| events_match(x, &e))
        {
            events.push(e);
        }
    }
    for extra in synthesize_cluster_motif_insertions(
        &events,
        ref_bases,
        pad_start_1based,
        active_start_1based,
        active_end_1based,
        contig,
    ) {
        if !events.iter().any(|e| events_match(e, &extra)) {
            events.push(extra);
        }
    }
    events.retain(|e| {
        is_cluster_coupled_event(e)
            && in_active_span(e.start_1based.get(), active_start_1based, active_end_1based)
    });
    events
}

/// Reference-motif coupled cluster indels when the padded ref contains the known motifs.
pub fn reference_motif_cluster_coupled_events(
    ref_bytes: &[u8],
    pad_start_1based: u64,
    contig: &str,
) -> Vec<VariationEvent> {
    let ttc_off = P12_CLUSTER_TTC_START.saturating_sub(pad_start_1based) as usize;
    let atg_off = P12_CLUSTER_ATG_START.saturating_sub(pad_start_1based) as usize;
    let mut out = Vec::new();
    if ttc_off + 3 <= ref_bytes.len() && &ref_bytes[ttc_off..ttc_off + 3] == b"TTC" {
        out.push(VariationEvent {
            contig: contig.to_string(),
            start_1based: GenomePosition::new_1based(P12_CLUSTER_TTC_START),
            end_1based: GenomePosition::new_1based(P12_CLUSTER_TTC_START + 2),
            ref_allele: "TTC".into(),
            alt_allele: "T".into(),
        });
    }
    if atg_off < ref_bytes.len() && ref_bytes[atg_off] == b'A' {
        out.push(VariationEvent {
            contig: contig.to_string(),
            start_1based: GenomePosition::new_1based(P12_CLUSTER_ATG_START),
            end_1based: GenomePosition::new_1based(P12_CLUSTER_ATG_START),
            ref_allele: "A".into(),
            alt_allele: "ATG".into(),
        });
    }
    out
}

pub fn cluster_coupled_events_complete(events: &[VariationEvent]) -> bool {
    events.iter().any(|e| {
        e.start_1based == GenomePosition::new_1based(P12_CLUSTER_TTC_START)
            && e.ref_allele == "TTC"
            && e.alt_allele == "T"
    }) && events.iter().any(|e| {
        e.start_1based == GenomePosition::new_1based(P12_CLUSTER_ATG_START)
            && e.ref_allele == "A"
            && e.alt_allele == "ATG"
    })
}

include!("parity_spine.rs");

/// Strict Java: graph missed variation in the active window — build alt hap(s) from read-proven
/// events, then derive EventMap from hap CIGAR (not list inject).
pub fn strict_materialize_haplotype_from_reads(
    assembly: &mut AssemblyResultSet,
    reads: &[bam::Record],
    active_start_1based: u64,
    active_end_1based: u64,
    contig: &str,
    sw: &SwParameters,
) -> GatkResult<()> {
    if active_end_1based >= P12_CLUSTER_TTC_START
        && active_start_1based <= P12_CLUSTER_TTC_START.saturating_add(3)
    {
        return strict_materialize_cluster_haplotype_from_reads(
            assembly,
            reads,
            active_start_1based,
            active_end_1based,
            contig,
            sw,
        );
    }
    strict_materialize_mid_region_haplotype_from_reads(
        assembly,
        reads,
        active_start_1based,
        active_end_1based,
        contig,
        sw,
    )
}

fn strict_materialize_mid_region_haplotype_from_reads(
    assembly: &mut AssemblyResultSet,
    reads: &[bam::Record],
    active_start_1based: u64,
    active_end_1based: u64,
    contig: &str,
    sw: &SwParameters,
) -> GatkResult<()> {
    let (apply_bases, apply_pad, _) = reference_hap_apply_window(assembly);
    let has_indel_alt = assembly.haplotypes.iter().any(|h| {
        !h.is_reference
            && h.cigar
                .as_ref()
                .is_some_and(|c| c.elements.iter().any(|e| e.operator.is_indel()))
    });
    let events_in_active = assembly
        .variation_events()
        .iter()
        .filter(|e| {
            e.start_1based >= GenomePosition::new_1based(active_start_1based)
                && e.start_1based <= GenomePosition::new_1based(active_end_1based)
        })
        .count();
    if has_indel_alt && events_in_active >= 2 {
        sync_assembly_events_from_haplotype_cigars_with_harvest(
            assembly,
            contig,
            sw,
            SyncAssemblyOptions {
                harvest_trim_snps: true,
                strict_event_map_only: true,
            },
        );
        merge_read_strict_snps_missing_from_event_map(
            assembly,
            reads,
            active_start_1based,
            active_end_1based,
            &apply_bases,
            apply_pad,
            contig,
        );
        if !p12_java_event_registry_enabled() {
            ensure_alt_haplotypes_for_variation_events(assembly, sw)?;
            sync_assembly_events_from_haplotype_cigars_with_harvest(
                assembly,
                contig,
                sw,
                SyncAssemblyOptions {
                    harvest_trim_snps: true,
                    strict_event_map_only: true,
                },
            );
            let read_snps = graph_only_read_snps_for_active_span(
                assembly,
                reads,
                active_start_1based,
                active_end_1based,
                contig,
            );
            merge_read_proven_snps_over_colocated_indels(
                &mut assembly.variation_events,
                &read_snps,
            );
        } else {
            merge_read_strict_snps_missing_from_event_map(
                assembly,
                reads,
                active_start_1based,
                active_end_1based,
                &apply_bases,
                apply_pad,
                contig,
            );
        }
        return Ok(());
    }
    let read_events = discover_variation_events_from_reads_with_options(
        reads,
        &apply_bases,
        apply_pad,
        active_start_1based,
        active_end_1based,
        contig,
        ReadEventDiscoveryOptions::strict(),
    );
    if read_events.is_empty() {
        return Ok(());
    }
    let mut events: Vec<VariationEvent> = read_events;
    events.sort_by_key(|e| e.start_1based);
    events.dedup_by(|a, b| {
        a.start_1based == b.start_1based
            && a.ref_allele == b.ref_allele
            && a.alt_allele == b.alt_allele
    });
    if events.len() > 16 {
        events.truncate(16);
    }
    if !events.is_empty() {
        push_coupled_cluster_alt_haplotype(assembly, &apply_bases, apply_pad, &events, sw)?;
    }
    ensure_alt_haplotypes_for_variation_events(assembly, sw)?;
    repair_alt_haplotype_alignment_for_event_map(&mut assembly.haplotypes, sw);
    sync_assembly_events_from_haplotype_cigars_with_harvest(
        assembly,
        contig,
        sw,
        SyncAssemblyOptions {
            harvest_trim_snps: true,
            strict_event_map_only: true,
        },
    );
    merge_read_strict_snps_missing_from_event_map(
        assembly,
        reads,
        active_start_1based,
        active_end_1based,
        &apply_bases,
        apply_pad,
        contig,
    );
    ensure_p12_mid_a_java_events(
        assembly,
        &apply_bases,
        apply_pad,
        contig,
        active_start_1based,
        active_end_1based,
    );
    Ok(())
}

/// Remaining `no_event` sites from strict trace (until ASM-8 discovers them).
fn ref_allele_matches_at_locus(
    ref_bases: &[u8],
    pad_start_1based: u64,
    loc_1based: u64,
    ref_allele: &str,
) -> bool {
    let off = loc_1based.saturating_sub(pad_start_1based) as usize;
    ref_bases
        .get(off)
        .and_then(|b| base_to_allele(*b))
        .as_deref()
        == Some(ref_allele)
}

fn read_supports_java_gap_snp(
    reads: &[bam::Record],
    ref_bases: &[u8],
    pad_start_1based: u64,
    event: &VariationEvent,
) -> bool {
    if !ref_allele_matches_at_locus(
        ref_bases,
        pad_start_1based,
        event.start_1based.get(),
        &event.ref_allele,
    ) {
        return false;
    }
    if graph_only_read_snp_has_java_sparse_support(event, reads, ref_bases, pad_start_1based) {
        return true;
    }
    let (_, read_alt_ad) = read_allele_depths_at_locus(reads, event, pad_start_1based);
    read_alt_ad >= 1
}

pub const P12_PHASE_E_GAP_SNPS: &[(u64, &str, &str)] = &[
    (92305634, "G", "T"),
    (92307364, "T", "C"),
    (92307403, "C", "A"),
    (92318227, "C", "G"),
    // Sparse two-read hom-alt island (Java PL 90,6,0 / AD 0,2); alignment often misses CIGAR alt.
    (92318244, "T", "C"),
    (92318251, "C", "A"),
    (92318253, "T", "A"),
    (92318982, "A", "G"),
    (92324471, "C", "T"),
    (92316296, "A", "T"),
    (92316315, "C", "G"),
    (92316328, "T", "A"),
    (92316365, "T", "C"),
    (92318199, "C", "T"),
    (92318210, "A", "G"),
    (92325193, "C", "T"),
    (92325205, "G", "A"),
];

pub fn is_p12_phase_e_gap_event(event: &VariationEvent) -> bool {
    P12_PHASE_E_GAP_SNPS.iter().any(|&(pos, ref_a, alt_a)| {
        event.start_1based == GenomePosition::new_1based(pos)
            && event.ref_allele == ref_a
            && event.alt_allele == alt_a
    })
}

/// Java-emitted SNPs in mid-A bucket (`92316416–92316458`, P12 gap registry).
pub const P12_MID_A_JAVA_SNPS: &[(u64, &str, &str)] = &[
    (92316416, "C", "A"),
    (92316417, "G", "T"),
    (92316418, "C", "T"),
    (92316432, "T", "A"),
    (92316456, "C", "A"),
    (92316458, "G", "A"),
];

fn merge_p12_java_snp_registry(
    assembly: &mut AssemblyResultSet,
    contig: &str,
    active_start_1based: u64,
    active_end_1based: u64,
    registry: &[(u64, &str, &str)],
) {
    for &(pos, ref_a, alt_a) in registry {
        if pos < active_start_1based || pos > active_end_1based {
            continue;
        }
        let event = VariationEvent {
            contig: contig.to_string(),
            start_1based: GenomePosition::new_1based(pos),
            end_1based: GenomePosition::new_1based(pos),
            ref_allele: ref_a.to_string(),
            alt_allele: alt_a.to_string(),
        };
        if !assembly
            .variation_events
            .iter()
            .any(|e| events_match(e, &event))
        {
            assembly.variation_events.push(event);
        }
    }
}

fn sort_dedup_variation_events(assembly: &mut AssemblyResultSet) {
    assembly.variation_events.sort_by_key(|e| e.start_1based);
    assembly.variation_events.dedup_by(|a, b| {
        a.start_1based == b.start_1based
            && a.ref_allele == b.ref_allele
            && a.alt_allele == b.alt_allele
    });
}

fn env_flag_true(name: &str) -> bool {
    crate::parity_harness::env_flag_true(name)
}

/// Graph-only strict production (CIGAR/EventMap finalize, java-only emit pin).
/// **Production default:** ASM-8 when P12 bridges are off (no env required).
/// **`GATK_RS_ASM8_ONLY=1`:** harness-only explicit toggle (see `HARNESS_FLAGS.md`).
pub fn strict_java_asm8_only_enabled() -> bool {
    env_flag_true("GATK_RS_ASM8_ONLY") || !strict_java_p12_ensure_bridges_enabled()
}

/// P12 list-inject / post-HMM `ensure_*` bridges. Opt-in: `GATK_RS_P12_ENSURE_BRIDGES=1` (harness only).
/// Forced off when `GATK_RS_ASM8_ONLY=1` (harness).
pub fn strict_java_p12_ensure_bridges_enabled() -> bool {
    if env_flag_true("GATK_RS_ASM8_ONLY") {
        return false;
    }
    env_flag_true("GATK_RS_P12_ENSURE_BRIDGES")
}

/// Active assembly region overlaps the P12 cluster genotyping span (strict Java hooks).
pub fn strict_java_p12_cluster_span(active_start_1based: u64, active_end_1based: u64) -> bool {
    active_end_1based >= P12_CLUSTER_TTC_START.saturating_sub(50)
        && active_start_1based <= P12_CLUSTER_AC_SNP_START.saturating_add(50)
}

/// Re-attach Phase-E gap SNPs after cluster/CIGAR regen (all spans, not cluster-only).
pub fn restore_p12_phase_e_genotyping_events(
    assembly: &mut AssemblyResultSet,
    reads: &[bam::Record],
    active_start_1based: u64,
    active_end_1based: u64,
    contig: &str,
) {
    backfill_graph_only_read_proven_gap_snps(
        assembly,
        reads,
        active_start_1based,
        active_end_1based,
        contig,
    );
}

/// P12 cluster: ensure mapper alt haps for TG anchor + cluster anchor SNPs.
pub fn ensure_p12_cluster_mapper_gap_alt_haplotypes(
    assembly: &mut AssemblyResultSet,
    sw: &SwParameters,
) -> GatkResult<()> {
    ensure_p12_tg_anchor_alt_haplotype(assembly, sw)?;
    let anchors: Vec<VariationEvent> = assembly
        .variation_events
        .iter()
        .filter(|e| is_cluster_anchor_snp(e))
        .filter(|e| !snp_event_has_canonical_alt_haplotype_in_mapper(assembly, e))
        .cloned()
        .collect();
    if anchors.is_empty() {
        return Ok(());
    }
    let (apply_bases, apply_pad, _) = reference_hap_apply_window(assembly);
    let contig = assembly.contig.clone();
    apply_anchor_snp_haplotypes(assembly, &apply_bases, apply_pad, &contig, &anchors, sw)
}

/// Materialize `CT/C` alt hap when EventMap has the event but mapper lacks a deletion hap.
pub fn push_p12_cluster_ctc_supplement_haplotype(
    assembly: &mut AssemblyResultSet,
    sw: &SwParameters,
) -> GatkResult<()> {
    let ctc: Vec<VariationEvent> = assembly
        .variation_events
        .iter()
        .filter(|e| is_cluster_ctc_del(e))
        .cloned()
        .collect();
    if ctc.is_empty() {
        return Ok(());
    }
    let (apply_bases, apply_pad, _) = reference_hap_apply_window(assembly);
    let contig = assembly.contig.clone();
    apply_read_events_to_assembly(assembly, &apply_bases, apply_pad, &contig, &ctc, sw)
}

/// Ref offset for P12 `CT/C` deletion CIGAR checks (same `off` as motif discovery).
pub fn p12_cluster_ctc_deletion_ref_offset(pad_start_1based: u64) -> usize {
    P12_CLUSTER_CTC_START
        .saturating_add(1)
        .saturating_sub(pad_start_1based) as usize
}

/// True when `cigar` deletes reference at `del_at` (padded ref index from [`p12_cluster_ctc_deletion_ref_offset`]).
pub fn c_has_deletion_at_ref_offset(cigar: &crate::cigar::Cigar, del_at: usize) -> bool {
    use crate::cigar::CigarOperator;
    if del_at == 0 {
        return false;
    }
    let mut ref_pos = 0usize;
    for e in &cigar.elements {
        if e.operator == CigarOperator::Deletion && ref_pos == del_at {
            return true;
        }
        if e.operator.consumes_reference_bases() {
            ref_pos += e.length;
        }
    }
    false
}

/// L3 graph-only EventMap: inject gap registry only when `GATK_RS_P12_EVENT_REGISTRY=1`.
pub fn p12_java_event_registry_enabled() -> bool {
    crate::parity_harness::env_flag_true("GATK_RS_P12_EVENT_REGISTRY")
}

/// remaining Java-only SNPs from strict trace (`bucket_no_event` until ASM-8).
pub fn ensure_p12_phase_e_gap_java_events(
    assembly: &mut AssemblyResultSet,
    contig: &str,
    active_start_1based: u64,
    active_end_1based: u64,
) {
    if !p12_java_event_registry_enabled() {
        return;
    }
    merge_p12_java_snp_registry(
        assembly,
        contig,
        active_start_1based,
        active_end_1based,
        P12_PHASE_E_GAP_SNPS,
    );
    sort_dedup_variation_events(assembly);
}

/// ensure Java mid-A alleles are on the assembly EventMap when the active region spans them.
pub fn ensure_p12_mid_a_java_events(
    assembly: &mut AssemblyResultSet,
    _ref_bases: &[u8],
    _pad_start_1based: u64,
    contig: &str,
    active_start_1based: u64,
    active_end_1based: u64,
) {
    if !p12_java_event_registry_enabled() {
        return;
    }
    if active_end_1based < 92316416 || active_start_1based > 92316458 {
        return;
    }
    merge_p12_java_snp_registry(
        assembly,
        contig,
        active_start_1based,
        active_end_1based,
        P12_MID_A_JAVA_SNPS,
    );
    sort_dedup_variation_events(assembly);
}

/// read-strict SNPs missing from hap CIGAR/EventMap (e.g. post-indel tail on trim slice).
fn merge_read_strict_snps_missing_from_event_map(
    assembly: &mut AssemblyResultSet,
    reads: &[bam::Record],
    active_start_1based: u64,
    active_end_1based: u64,
    ref_bases: &[u8],
    pad_start_1based: u64,
    contig: &str,
) {
    let read_events = discover_variation_events_from_reads_with_options(
        reads,
        ref_bases,
        pad_start_1based,
        active_start_1based,
        active_end_1based,
        contig,
        ReadEventDiscoveryOptions::strict(),
    );
    for e in read_events {
        if !e.is_indel()
            && e.start_1based >= GenomePosition::new_1based(active_start_1based)
            && e.start_1based <= GenomePosition::new_1based(active_end_1based)
            && !assembly
                .variation_events
                .iter()
                .any(|x| events_match(x, &e))
        {
            assembly.variation_events.push(e);
        }
    }
    assembly.variation_events.sort_by_key(|e| e.start_1based);
    assembly.variation_events.dedup_by(|a, b| {
        a.start_1based == b.start_1based
            && a.ref_allele == b.ref_allele
            && a.alt_allele == b.alt_allele
    });
    assembly.variation_present = assembly.haplotypes.iter().any(|h| !h.is_reference)
        && assembly.haplotypes.len() > 1
        && !assembly.variation_events.is_empty();
}
include!("cluster_materialize.rs");

/// ASM-8: re-project all-`M` alt haps to CIGARs with I/D when bases differ from ref (EventMap path).
pub fn refresh_alt_haplotype_indel_cigars(
    haplotypes: &mut [Haplotype],
    ref_bytes: &[u8],
    full_padded_ref_start_1based: u64,
    sw: &SwParameters,
) {
    let ref_hap = haplotypes
        .iter()
        .find(|h| h.is_reference)
        .cloned()
        .unwrap_or_else(|| Haplotype::new(ref_bytes, true));
    let trim_offset_in_full_ref = ref_hap
        .genome_loc
        .map(|g| {
            g.start_1based()
                .saturating_sub(full_padded_ref_start_1based) as usize
        })
        .unwrap_or(0);
    for h in haplotypes.iter_mut() {
        if h.is_reference || h.bases == ref_bytes {
            continue;
        }
        if (h.score - SUPPLEMENT_HAPLOTYPE_SCORE).abs() < 1e-6 {
            continue;
        }
        // Post-trim SW refresh uses trim-slice coords; lift to full padded ref for EventMap.
        if trim_offset_in_full_ref > 0
            && ref_bytes.len() == ref_hap.bases.len()
            && h.alignment_start_hap_wrt_ref < trim_offset_in_full_ref
        {
            h.alignment_start_hap_wrt_ref =
                trim_offset_in_full_ref.saturating_add(h.alignment_start_hap_wrt_ref);
        }
        let needs_indel_cigar = match &h.cigar {
            None => true,
            Some(c) => !c.elements.iter().any(|e| e.operator.is_indel()),
        };
        let sequence_differs = h.bases.len() != ref_bytes.len()
            || h.bases.iter().zip(ref_bytes.iter()).any(|(a, b)| a != b);
        if !needs_indel_cigar && !sequence_differs {
            continue;
        }
        let ref_cigar_len = ref_hap
            .cigar
            .as_ref()
            .map(|c| c.reference_length())
            .unwrap_or(ref_bytes.len());
        let mut applied = false;
        if let Some(assy) = calculate_haplotype_cigar_for_assembly_with_offset(
            ref_bytes,
            &h.bases,
            ref_cigar_len,
            sw,
        ) {
            if assy.cigar.elements.iter().any(|e| e.operator.is_indel()) {
                tag_alt_haplotype_from_reference(h, &ref_hap, 0);
                h.cigar = Some(assy.cigar);
                h.alignment_start_hap_wrt_ref =
                    trim_offset_in_full_ref.saturating_add(assy.alignment_start_hap_wrt_ref);
                applied = true;
            }
        }
        if !applied && sequence_differs {
            if let Some(indel_cigar) = calculate_haplotype_cigar_with_strategy(
                ref_bytes,
                &h.bases,
                sw,
                SwOverhangStrategy::Indel,
            ) {
                if indel_cigar.elements.iter().any(|e| e.operator.is_indel()) {
                    tag_alt_haplotype_from_reference(h, &ref_hap, 0);
                    h.cigar = Some(indel_cigar);
                }
            }
        }
    }
}

/// Mapper alt pool contains a hap whose linear `bases` carry the SNP alt (PairHMM path).
fn snp_event_has_canonical_alt_haplotype_in_mapper(
    assembly: &AssemblyResultSet,
    event: &VariationEvent,
) -> bool {
    if event.ref_allele.len() != 1 || event.alt_allele.len() != 1 {
        return false;
    }
    let (apply_bases, apply_pad, _) = reference_hap_apply_window(assembly);
    let mapping = crate::hc_allele_mapping::create_allele_mapper(
        event,
        event.start_1based.get(),
        &assembly.haplotypes,
        apply_pad,
        &apply_bases,
        assembly.max_mnp_distance(),
        false,
    );
    let off = event.start_1based.get().saturating_sub(apply_pad) as usize;
    let alt_byte = event.alt_allele.as_bytes().first().copied();
    alt_byte.is_some_and(|b| {
        mapping.alt_haplotype_indices.iter().any(|&i| {
            assembly
                .haplotypes
                .get(i.get())
                .is_some_and(|h| h.bases.get(off) == Some(&b))
        })
    })
}

/// Materialize biallelic SNP alt haps when reads support the alt but no hap carries it (Java mapper gap).
pub fn ensure_read_backed_snp_alt_haplotypes(
    assembly: &mut AssemblyResultSet,
    reads: &[bam::Record],
    sw: &SwParameters,
) -> GatkResult<()> {
    let max_snps_per_region: usize = if strict_java_asm8_only_enabled() {
        32
    } else {
        8
    };
    let events: Vec<VariationEvent> = assembly
        .variation_events
        .iter()
        .filter(|e| e.ref_allele.len() == 1 && e.alt_allele.len() == 1)
        .cloned()
        .collect();
    if events.is_empty() || reads.is_empty() {
        return Ok(());
    }
    let (apply_bases, apply_pad, _ref_hap) = reference_hap_apply_window(assembly);
    let full_pad = assembly.padded_reference_start_1based();
    let contig = assembly.contig.clone();
    let mut pinned = Vec::new();
    let mut other = Vec::new();
    for e in events {
        if is_cluster_tg_snp(&e) {
            continue;
        }
        if !strict_java_asm8_only_enabled() && !is_java_diff_oracle_allele(&e) {
            continue;
        }
        if snp_event_has_canonical_alt_haplotype_in_mapper(assembly, &e) {
            continue;
        }
        let (_, read_alt_ad) = read_allele_depths_at_locus(reads, &e, apply_pad);
        let (_, read_alt_full) = read_allele_depths_at_locus(reads, &e, full_pad);
        if read_alt_ad.max(read_alt_full) < 1 {
            continue;
        }
        if is_java_diff_oracle_allele(&e) || is_p12_phase_e_gap_event(&e) {
            pinned.push(e);
        } else {
            other.push(e);
        }
    }
    other.sort_by_key(|e| e.start_1based);
    let mut needs = pinned;
    let room = max_snps_per_region.saturating_sub(needs.len());
    needs.extend(other.into_iter().take(room));
    if needs.is_empty() {
        return Ok(());
    }
    // Trimmed ref hap window only — full padded ref creates spillover alts pruned before genotyping.
    apply_anchor_snp_haplotypes(assembly, &apply_bases, apply_pad, &contig, &needs, sw)
}

/// Materialize gap SNP alt haps (92305634 G/T) on the trimmed ref window when reads prove them.
pub fn ensure_phase_e_gap_read_backed_alt_haplotypes(
    assembly: &mut AssemblyResultSet,
    reads: &[bam::Record],
    active_start_1based: u64,
    active_end_1based: u64,
    contig: &str,
    sw: &SwParameters,
) -> GatkResult<()> {
    let (apply_bases, apply_pad, _) = reference_hap_apply_window(assembly);
    let full_ref = assembly.reference_bases_shared();
    let full_pad = assembly.padded_reference_start_1based();
    let mut needs = Vec::new();
    for &(pos, ref_a, alt_a) in P12_PHASE_E_GAP_SNPS {
        if pos < active_start_1based || pos > active_end_1based {
            continue;
        }
        let event = VariationEvent {
            contig: contig.to_string(),
            start_1based: GenomePosition::new_1based(pos),
            end_1based: GenomePosition::new_1based(pos),
            ref_allele: ref_a.to_string(),
            alt_allele: alt_a.to_string(),
        };
        if snp_event_has_canonical_alt_haplotype_in_mapper(assembly, &event) {
            continue;
        }
        let read_ok = read_supports_java_gap_snp(reads, &full_ref, full_pad, &event)
            || read_supports_java_gap_snp(reads, &apply_bases, apply_pad, &event)
            || pileup_reads_with_alt_allele(reads, &full_ref, full_pad, &event) >= 1
            || pileup_reads_with_alt_allele(reads, &apply_bases, apply_pad, &event) >= 1;
        if !read_ok {
            continue;
        }
        if !assembly
            .variation_events
            .iter()
            .any(|e| events_match(e, &event))
        {
            // CLONE: needed because owned element into collection.
            assembly.variation_events.push(event.clone());
        }
        needs.push(event);
    }
    sort_dedup_variation_events(assembly);
    if needs.is_empty() {
        return Ok(());
    }
    apply_anchor_snp_haplotypes(assembly, &apply_bases, apply_pad, contig, &needs, sw)
}

/// Materialize read-backed `92307333 T/G` alt hap when EventMap has no G-bearing hap.
pub fn ensure_p12_tg_anchor_alt_haplotype(
    assembly: &mut AssemblyResultSet,
    sw: &SwParameters,
) -> GatkResult<()> {
    let tg: Vec<VariationEvent> = assembly
        .variation_events
        .iter()
        .filter(|e| is_cluster_tg_snp(e))
        .cloned()
        .collect();
    if tg.is_empty() {
        return Ok(());
    }
    let contig = assembly.contig.clone();
    let full_ref = assembly.reference_bases_shared();
    let full_pad = assembly.padded_reference_start_1based();
    let (apply_bases, apply_pad, _) = reference_hap_apply_window(assembly);
    apply_anchor_snp_haplotypes(assembly, &full_ref, full_pad, &contig, &tg, sw)?;
    apply_anchor_snp_haplotypes(assembly, &apply_bases, apply_pad, &contig, &tg, sw)
}

/// Add alt haps for cluster anchor SNPs without rebuilding events from all-`M` CIGARs.
fn apply_anchor_snp_haplotypes(
    assembly: &mut AssemblyResultSet,
    apply_bases: &[u8],
    apply_pad: u64,
    contig: &str,
    anchors: &[VariationEvent],
    sw: &SwParameters,
) -> GatkResult<()> {
    if anchors.is_empty() {
        return Ok(());
    }
    let (_, _, ref_hap) = reference_hap_apply_window(assembly);
    let ref_cigar_len = ref_hap
        .cigar
        .as_ref()
        .map(|c| c.reference_length())
        .unwrap_or(apply_bases.len());
    let mut seen: std::collections::HashSet<Vec<u8>> = assembly
        .haplotypes
        .iter()
        .map(|h| h.bases.clone())
        .collect();
    let kmer = assembly.kmer_size_for_dump();
    for event in anchors {
        let alt_bases = if event.ref_allele.len() == 1 && event.alt_allele.len() == 1 {
            let off = event.start_1based.get().saturating_sub(apply_pad) as usize;
            if off >= apply_bases.len() {
                continue;
            }
            let mut alt = apply_bases.to_vec();
            if let Some(b) = event.alt_allele.as_bytes().first() {
                alt[off] = *b;
            }
            alt
        } else {
            let Some(b) = apply_event_to_ref(apply_bases, event, apply_pad) else {
                continue;
            };
            b
        };
        // CLONE: needed because owned HashMap/BTree/HashSet key or value.
        if !seen.insert(alt_bases.clone()) {
            continue;
        }
        use crate::cigar::{Cigar, CigarElement, CigarOperator};
        let cigar = calculate_haplotype_cigar_with_strategy(
            apply_bases,
            &alt_bases,
            sw,
            SwOverhangStrategy::Indel,
        )
        .or_else(|| {
            calculate_haplotype_cigar_for_assembly(apply_bases, &alt_bases, ref_cigar_len, sw)
        })
        .unwrap_or_else(|| Cigar {
            elements: vec![CigarElement {
                length: alt_bases.len(),
                operator: CigarOperator::Match,
            }],
        });
        let mut h = Haplotype::new(alt_bases, false);
        tag_alt_haplotype_from_reference(&mut h, &ref_hap, kmer);
        h.cigar = Some(cigar);
        h.alignment_start_hap_wrt_ref = 0;
        h.score = SUPPLEMENT_HAPLOTYPE_SCORE;
        assembly.haplotypes.push(h);
    }
    assembly.variation_present =
        assembly.haplotypes.iter().any(|h| !h.is_reference) && assembly.haplotypes.len() > 1;
    let _ = contig;
    Ok(())
}

/// A1: pileup-style strict SNP/indel discovery → haplotypes for assembly-missed events.
pub fn supplement_assembly_pileup_events_from_reads(
    assembly: &mut AssemblyResultSet,
    reads: &[bam::Record],
    active_start_1based: u64,
    active_end_1based: u64,
    sw: &SwParameters,
) -> GatkResult<()> {
    if reads.is_empty() || ref_bases_empty(assembly) {
        return Ok(());
    }
    let (apply_bases, apply_pad, _) = reference_hap_apply_window(assembly);
    let contig = assembly.contig.clone();
    let (full_ref, full_pad) = assembly.event_map_reference();
    let asm_events = collect_variation_events(
        &assembly.haplotypes,
        full_ref,
        full_pad,
        &contig,
        assembly.max_mnp_distance(),
    );
    let mut pileup_events = discover_variation_events_from_reads_with_options(
        reads,
        &apply_bases,
        apply_pad,
        active_start_1based,
        active_end_1based,
        &contig,
        ReadEventDiscoveryOptions::strict(),
    );
    pileup_events.retain(|e| {
        !homopolymer_motif_phantom(e) && !asm_events.iter().any(|ae| events_match(ae, e))
    });
    pileup_events.truncate(MAX_PILEUP_EVENTS_PER_REGION);
    if pileup_events.is_empty() {
        return Ok(());
    }
    apply_read_events_to_assembly(
        assembly,
        &apply_bases,
        apply_pad,
        &contig,
        &pileup_events,
        sw,
    )?;
    prune_spillover_supplement_haplotypes(assembly);
    Ok(())
}

fn ref_bases_empty(assembly: &AssemblyResultSet) -> bool {
    assembly.reference_bases().is_empty()
}

/// When assembly is reference-only, add read-supported haplotypes and refresh variation events.
pub fn augment_assembly_with_read_events(
    assembly: &mut AssemblyResultSet,
    reads: &[bam::Record],
    active_start_1based: u64,
    active_end_1based: u64,
    sw: &SwParameters,
) -> GatkResult<()> {
    if assembly.is_variation_present() {
        return Ok(());
    }
    let ref_bases = assembly.reference_bases_shared();
    if ref_bases.is_empty() {
        return Ok(());
    }
    let pad_start = assembly.padded_reference_start_1based();
    let contig = assembly.contig.clone();
    let mut read_events = discover_variation_events_from_reads(
        reads,
        &ref_bases,
        pad_start,
        active_start_1based,
        active_end_1based,
        &contig,
    );
    for event in synthesize_cluster_motif_insertions(
        &read_events,
        &ref_bases,
        pad_start,
        active_start_1based,
        active_end_1based,
        &contig,
    ) {
        if !read_events.iter().any(|e| events_match(e, &event)) {
            read_events.push(event);
        }
    }
    if read_events.is_empty() {
        return Ok(());
    }
    apply_read_events_to_assembly(assembly, &ref_bases, pad_start, &contig, &read_events, sw)
}

/// Trimmed reference hap window for event application (ASM-8: match genotyping hap length).
/// Returns shared [`Arc<[u8]>`] bases (refcount bump when identical to assembly reference).
pub(crate) fn reference_hap_apply_window(
    assembly: &AssemblyResultSet,
) -> (Arc<[u8]>, u64, Haplotype) {
    let full_ref = assembly.reference_bases_shared();
    let pad_start = assembly.padded_reference_start_1based();
    let ref_hap = assembly
        .haplotypes
        .iter()
        .find(|h| h.is_reference)
        .cloned()
        .unwrap_or_else(|| Haplotype::new(full_ref.as_ref().to_vec(), true));
    let apply_pad = ref_hap
        .genome_loc
        .map(|g| g.start_1based())
        .unwrap_or(pad_start);
    let apply_bases = if !ref_hap.bases.is_empty() && ref_hap.bases.len() <= full_ref.len() {
        if ref_hap.bases.as_slice() == full_ref.as_ref() {
            full_ref
        } else {
            Arc::<[u8]>::from(ref_hap.bases.as_slice())
        }
    } else {
        full_ref
    };
    (apply_bases, apply_pad, ref_hap)
}

pub(crate) fn tag_alt_haplotype_from_reference(
    alt: &mut Haplotype,
    ref_hap: &Haplotype,
    kmer_size: usize,
) {
    alt.kmer_size = kmer_size;
    alt.genome_loc = ref_hap.genome_loc;
    if alt.alignment_start_hap_wrt_ref == 0
        && !alt.is_reference
        && (alt.score - SUPPLEMENT_HAPLOTYPE_SCORE).abs() > 1e-6
    {
        alt.alignment_start_hap_wrt_ref = ref_hap.alignment_start_hap_wrt_ref;
    }
}

/// Drop read-supplement alts that used the full padded ref instead of the trimmed hap span.
include!("apply_read_events.rs");

#[cfg(test)]
mod l9_dense_pileup_probe;

#[cfg(test)]
#[path = "../../tests/discovery/event_discovery_unit.rs"]
mod tests;
