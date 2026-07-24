//! Read model semantics aligned with GATK `HaplotypeCallerEngine` read filters.
//! Java reference: `HaplotypeCallerEngine.makeStandardHCReadFilters` and `WellformedReadFilter`
//! (`CigarUtils.isGood`, `ReadFilterLibrary`, `ReadUtils.alignmentAgreesWithHeader`).

#![warn(clippy::unwrap_used, clippy::expect_used)]

use gatk_common::HaplotypeCallerConfig;
use rust_htslib::bam;
use rust_htslib::bam::record::Cigar;

/// SAM FLAG bits used for HC-style alignment decisions (SAMv1).
pub const FLAG_SEGMENT_UNMAPPED: u16 = 0x0004;
pub const FLAG_NOT_PRIMARY: u16 = 0x0100;
pub const FLAG_VENDOR_QUALITY_FAILED: u16 = 0x0200;
pub const FLAG_DUPLICATE: u16 = 0x0400;
pub const FLAG_SUPPLEMENTARY: u16 = 0x0800;

/// GATK `QualityUtils.MAPPING_QUALITY_UNAVAILABLE` — **rejected** by standard HC (after min-MQ check).
pub const MAPPING_QUALITY_UNAVAILABLE: u8 = 255;

/// GATK `HaplotypeCallerEngine.DEFAULT_READ_QUALITY_FILTER_THRESHOLD`.
pub const GATK_HC_DEFAULT_MIN_MAPPING_QUALITY: u8 = 20;

/// Which fully-defined Java read-filter chain to mirror (requires BAM [`bam::HeaderView`]).
/// # Invariants
/// `Standard` retains supplementary alignments; `SoftClipMinimal` filters them (parity shard path).
/// # Ownership
/// [`Copy`] enum selected via [`ReadFilterParams::resolved_hc_filter_set`].
/// # Mutation
/// Immutable filter-set discriminant.
/// # Biological assumptions
/// Filters remove low-quality / malformed reads before assembly and genotyping evidence.
/// # Java equivalence
/// GATK `HaplotypeCallerEngine.makeStandardHCReadFilters` vs soft-clip shard configuration.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HcReadFilterSet {
    /// `HaplotypeCallerEngine.makeStandardHCReadFilters` — supplementary alignments **allowed**.
    Standard,
    /// `HcFullParityGateDump.configureSoftClipReadShard` — supplementary **filtered**.
    SoftClipMinimal,
}

/// Parameters for overlap counting / legacy field-only checks (no `WellformedReadFilter`).
/// # Invariants
/// `min_mapping_quality` applies before MQ=255 unavailable rejection in full header chains.
/// `exclude_supplementary == false` matches Java standard HC (supplementary allowed).
/// # Ownership
/// [`Copy`] snapshot from CLI / [`HaplotypeCallerConfig`] conversion.
/// # Mutation
/// Immutable for the duration of a filter pass.
/// # Biological assumptions
/// Primary-alignment read set for HC with MAPQ and duplicate/secondary flags as configured.
/// # Java equivalence
/// GATK `HaplotypeCallerEngine.makeStandardHCReadFilters` parameterization slice.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ReadFilterParams {
    pub min_mapping_quality: u8,
    pub exclude_duplicates: bool,
    pub exclude_secondary: bool,
    pub exclude_supplementary: bool,
}

impl Default for ReadFilterParams {
    fn default() -> Self {
        Self {
            min_mapping_quality: 20,
            exclude_duplicates: true,
            exclude_secondary: true,
            exclude_supplementary: true,
        }
    }
}

impl ReadFilterParams {
    pub fn from_haplotype_caller(hc: &HaplotypeCallerConfig) -> Self {
        Self {
            min_mapping_quality: hc.min_mapping_quality.min(u8::MAX as u32) as u8,
            exclude_duplicates: true,
            exclude_secondary: true,
            // Java `makeStandardHCReadFilters` does **not** apply `NOT_SUPPLEMENTARY_ALIGNMENT`.
            exclude_supplementary: false,
        }
    }

    /// Structural flags matching Java standard HC filter list (supplementary retained).
    pub fn gatk_standard_hc() -> Self {
        Self {
            min_mapping_quality: GATK_HC_DEFAULT_MIN_MAPPING_QUALITY,
            exclude_duplicates: true,
            exclude_secondary: true,
            exclude_supplementary: false,
        }
    }

    /// Java `getPileupsOverReference` — `readLikelihoods.sampleEvidence` reads are already
    /// shard-filtered and realigned; RCM pileup must not re-apply GoodCigar/Wellformed.
    pub fn genotyping_evidence_rcm_pileup() -> Self {
        Self {
            min_mapping_quality: 0,
            exclude_duplicates: false,
            exclude_secondary: false,
            exclude_supplementary: false,
        }
    }

    /// Resolve to a [`HcReadFilterSet`] when `exclude_*` flags match a known Java bundle; else `None`.
    pub fn resolved_hc_filter_set(&self) -> Option<HcReadFilterSet> {
        if !(self.exclude_duplicates && self.exclude_secondary) {
            return None;
        }
        Some(if self.exclude_supplementary {
            HcReadFilterSet::SoftClipMinimal
        } else {
            HcReadFilterSet::Standard
        })
    }
}

/// `true` when `flags` indicate unmapped, secondary, or supplementary per [`ReadFilterParams`].
#[inline]
pub fn flags_exclude_from_primary_alignment(flags: u16, p: &ReadFilterParams) -> bool {
    if flags & FLAG_SEGMENT_UNMAPPED != 0 {
        return true;
    }
    if p.exclude_secondary && flags & FLAG_NOT_PRIMARY != 0 {
        return true;
    }
    if p.exclude_supplementary && flags & FLAG_SUPPLEMENTARY != 0 {
        return true;
    }
    false
}

/// Legacy helper: MQ threshold only (does **not** model `MAPPING_QUALITY_AVAILABLE`).
#[inline]
pub fn mapq_passes_minimum(mapq: u8, min_mapping_quality: u8) -> bool {
    mapq >= min_mapping_quality
}

/// Field-level version **without** `WellformedReadFilter` / full Java chain — for unit tests only.
pub fn passes_hc_read_filters_fields(flags: u16, mapq: u8, p: &ReadFilterParams) -> bool {
    if flags_exclude_from_primary_alignment(flags, p) {
        return false;
    }
    if p.exclude_duplicates && (flags & FLAG_DUPLICATE) != 0 {
        return false;
    }
    mapq_passes_minimum(mapq, p.min_mapping_quality)
}

/// Legacy ingress without header — **not** Java-identical when well-formedness matters.
pub fn passes_hc_read_filters(rec: &bam::Record, p: &ReadFilterParams) -> bool {
    passes_hc_read_filters_fields(rec.flags(), rec.mapq(), p)
}

/// Full Java filter chain when [`ReadFilterParams::resolved_hc_filter_set`] is `Some`, else [`passes_hc_read_filters`].
pub fn passes_hc_read_filters_with_header(
    rec: &bam::Record,
    header: &bam::HeaderView,
    p: &ReadFilterParams,
) -> bool {
    match p.resolved_hc_filter_set() {
        Some(set) => passes_hc_read_filter_set(rec, header, set, p.min_mapping_quality),
        None => passes_hc_read_filters(rec, p),
    }
}

/// Dispatch Java-identical chains (requires header for `WellformedReadFilter`).
pub fn passes_hc_read_filter_set(
    rec: &bam::Record,
    header: &bam::HeaderView,
    set: HcReadFilterSet,
    min_mapping_quality: u8,
) -> bool {
    match set {
        HcReadFilterSet::Standard => {
            passes_standard_hc_read_filters_inner(rec, header, min_mapping_quality)
        }
        HcReadFilterSet::SoftClipMinimal => {
            passes_soft_clip_shard_read_filters(rec, header, min_mapping_quality)
        }
    }
}

/// Java `Class.getSimpleName` for each delegate in `HaplotypeCallerEngine.makeStandardHCReadFilters`,
/// in the same order as `ReadFilter.fromList(...)` / `CountingReadFilter.fromList(...)` short-circuit evaluation.
pub const STANDARD_HC_READ_FILTER_JAVA_NAMES: [&str; 9] = [
    "MappingQualityReadFilter",
    "MappingQualityAvailableReadFilter",
    "MappedReadFilter",
    "NotSecondaryAlignmentReadFilter",
    "NotDuplicateReadFilter",
    "PassesVendorQualityCheckReadFilter",
    "NonZeroReferenceLengthAlignmentReadFilter",
    "GoodCigarReadFilter",
    "WellformedReadFilter",
];

/// Index into [`STANDARD_HC_READ_FILTER_JAVA_NAMES`] for the first failing predicate (Java `ReadFilterAnd`
/// short-circuit order). `None` if all pass.
pub fn standard_hc_read_filter_failure_index(
    rec: &bam::Record,
    header: &bam::HeaderView,
    min_mq: u8,
) -> Option<usize> {
    let mq = rec.mapq();
    if mq < min_mq {
        return Some(0);
    }
    if mq == MAPPING_QUALITY_UNAVAILABLE {
        return Some(1);
    }
    if rec.is_unmapped() {
        return Some(2);
    }
    let flags = rec.flags();
    if flags & FLAG_NOT_PRIMARY != 0 {
        return Some(3);
    }
    if flags & FLAG_DUPLICATE != 0 {
        return Some(4);
    }
    if flags & FLAG_VENDOR_QUALITY_FAILED != 0 {
        return Some(5);
    }
    let cigar: Vec<Cigar> = rec.cigar().to_vec();
    if !cigar_non_zero_reference_consumption(&cigar) {
        return Some(6);
    }
    if !cigar_is_good(&cigar) {
        return Some(7);
    }
    if !wellformed_read_filter(rec, header, &cigar) {
        return Some(8);
    }
    None
}

/// `HaplotypeCallerEngine.makeStandardHCReadFilters` in order (exact list; MQ params from caller).
fn passes_standard_hc_read_filters_inner(
    rec: &bam::Record,
    header: &bam::HeaderView,
    min_mq: u8,
) -> bool {
    standard_hc_read_filter_failure_index(rec, header, min_mq).is_none()
}

/// `HcFullParityGateDump.configureSoftClipReadShard` — **Wellformed first**, then MQ / duplicate / secondary / supplementary / mapped.
fn passes_soft_clip_shard_read_filters(
    rec: &bam::Record,
    header: &bam::HeaderView,
    min_mq: u8,
) -> bool {
    let cigar: Vec<Cigar> = rec.cigar().to_vec();
    if !wellformed_read_filter(rec, header, &cigar) {
        return false;
    }
    let mq = rec.mapq();
    if mq == MAPPING_QUALITY_UNAVAILABLE {
        return false;
    }
    if mq < min_mq {
        return false;
    }
    let flags = rec.flags();
    if flags & FLAG_DUPLICATE != 0 {
        return false;
    }
    if flags & FLAG_NOT_PRIMARY != 0 {
        return false;
    }
    if flags & FLAG_SUPPLEMENTARY != 0 {
        return false;
    }
    if rec.is_unmapped() {
        return false;
    }
    true
}

fn cigar_op_is_clip(op: &Cigar) -> bool {
    matches!(op, Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_))
}

fn cigar_has_consecutive_indels(cigar: &[Cigar]) -> bool {
    let mut prev_indel = false;
    for op in cigar {
        let is_indel = matches!(op, Cigar::Ins(_) | Cigar::Del(_));
        if prev_indel && is_indel {
            return true;
        }
        prev_indel = is_indel;
    }
    false
}

fn cigar_starts_or_ends_with_deletion_ignoring_clips(cigar: &[Cigar]) -> bool {
    let first_non_clip = cigar.iter().find(|op| !cigar_op_is_clip(op));
    if matches!(first_non_clip, Some(Cigar::Del(_))) {
        return true;
    }
    let last_non_clip = cigar.iter().rev().find(|op| !cigar_op_is_clip(op));
    matches!(last_non_clip, Some(Cigar::Del(_)))
}

/// GATK `CigarUtils.isGood` (BAM record has already-parseable CIGAR; skips htsjdk `Cigar.isValid` deep check).
fn cigar_is_good(cigar: &[Cigar]) -> bool {
    if cigar_has_consecutive_indels(cigar)
        || cigar_starts_or_ends_with_deletion_ignoring_clips(cigar)
    {
        return false;
    }
    true
}

fn cigar_non_zero_reference_consumption(cigar: &[Cigar]) -> bool {
    cigar.iter().any(|op| {
        let (consumes, len) = match op {
            Cigar::Match(n)
            | Cigar::Del(n)
            | Cigar::RefSkip(n)
            | Cigar::Equal(n)
            | Cigar::Diff(n) => (true, *n),
            Cigar::Ins(n) | Cigar::SoftClip(n) | Cigar::HardClip(n) => (false, *n),
            Cigar::Pad(n) => (true, *n),
        };
        consumes && len > 0
    })
}

fn cigar_contains_ref_skip_n(cigar: &[Cigar]) -> bool {
    cigar.iter().any(|op| matches!(op, Cigar::RefSkip(_)))
}

fn cigar_read_length_htsjdk(cigar: &[Cigar]) -> usize {
    cigar
        .iter()
        .map(|op| match op {
            Cigar::Match(n)
            | Cigar::Ins(n)
            | Cigar::SoftClip(n)
            | Cigar::Equal(n)
            | Cigar::Diff(n) => *n as usize,
            _ => 0usize,
        })
        .sum()
}

fn alignment_agrees_with_header(rec: &bam::Record, header: &bam::HeaderView) -> bool {
    if rec.is_unmapped() {
        return true;
    }
    let tid = rec.tid();
    if tid < 0 {
        return false;
    }
    let Some(ctg_len) = header.target_len(tid as u32) else {
        return false;
    };
    let start1 = rec.pos() as u64 + 1;
    start1 <= ctg_len
}

fn has_read_group(rec: &bam::Record) -> bool {
    rec.aux(b"RG").is_ok()
}

/// GATK `WellformedReadFilter` AND-chain.
fn wellformed_read_filter(rec: &bam::Record, header: &bam::HeaderView, cigar: &[Cigar]) -> bool {
    if !(rec.is_unmapped() || rec.pos() + 1 > 0) {
        return false;
    }
    if !(rec.is_unmapped() || rec.cigar().end_pos() >= rec.pos()) {
        return false;
    }
    if !alignment_agrees_with_header(rec, header) {
        return false;
    }
    if !has_read_group(rec) {
        return false;
    }
    let seq_len = rec.seq().len();
    let qual_len = rec.qual().len();
    if seq_len != qual_len {
        return false;
    }
    if !(rec.is_unmapped() || seq_len == cigar_read_length_htsjdk(cigar)) {
        return false;
    }
    if seq_len == 0 {
        return false;
    }
    if cigar_contains_ref_skip_n(cigar) {
        return false;
    }
    true
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]
mod tests {
    use super::*;
    use rust_htslib::bam::Read as _;

    #[test]
    fn mapq_below_min_fails() {
        assert!(!mapq_passes_minimum(19, 20));
        assert!(mapq_passes_minimum(20, 20));
    }

    #[test]
    fn flags_supplementary_secondary_unmapped_excluded() {
        let p = ReadFilterParams::default();
        assert!(flags_exclude_from_primary_alignment(FLAG_SUPPLEMENTARY, &p));
        assert!(flags_exclude_from_primary_alignment(FLAG_NOT_PRIMARY, &p));
        assert!(flags_exclude_from_primary_alignment(
            FLAG_SEGMENT_UNMAPPED,
            &p
        ));
        assert!(!flags_exclude_from_primary_alignment(0, &p));
        assert!(!flags_exclude_from_primary_alignment(FLAG_DUPLICATE, &p));
        let gatk = ReadFilterParams::gatk_standard_hc();
        assert!(!flags_exclude_from_primary_alignment(
            FLAG_SUPPLEMENTARY,
            &gatk
        ));
    }

    #[test]
    fn read_filter_params_from_hc_clamps_mapq() {
        let hc = gatk_common::HaplotypeCallerConfig {
            min_mapping_quality: 500,
            ..Default::default()
        };
        let p = ReadFilterParams::from_haplotype_caller(&hc);
        assert_eq!(p.min_mapping_quality, 255);
    }

    #[test]
    fn field_level_filter_matrix() {
        let p = ReadFilterParams {
            min_mapping_quality: 20,
            exclude_duplicates: true,
            exclude_secondary: true,
            exclude_supplementary: true,
        };
        assert!(passes_hc_read_filters_fields(0, 20, &p));
        assert!(passes_hc_read_filters_fields(
            0,
            MAPPING_QUALITY_UNAVAILABLE,
            &p
        ));
        assert!(!passes_hc_read_filters_fields(0, 19, &p));
        assert!(!passes_hc_read_filters_fields(FLAG_NOT_PRIMARY, 60, &p));
        assert!(!passes_hc_read_filters_fields(FLAG_SUPPLEMENTARY, 60, &p));
        assert!(!passes_hc_read_filters_fields(
            FLAG_SEGMENT_UNMAPPED,
            60,
            &p
        ));
        assert!(!passes_hc_read_filters_fields(FLAG_DUPLICATE, 60, &p));
    }

    #[test]
    fn duplicate_can_be_allowed() {
        let p = ReadFilterParams {
            min_mapping_quality: 20,
            exclude_duplicates: false,
            exclude_secondary: true,
            exclude_supplementary: true,
        };
        assert!(passes_hc_read_filters_fields(FLAG_DUPLICATE, 60, &p));
    }

    #[test]
    fn standard_hc_failure_index_matches_read_filter_set_on_fixture() {
        use std::path::Path;
        let sam =
            Path::new(env!("CARGO_MANIFEST_DIR")).join("../parity/fixtures/read_filter_slice.sam");
        let mut reader = bam::Reader::from_path(&sam).unwrap();
        let header = reader.header().clone();
        for res in reader.records() {
            let rec = res.unwrap();
            let pass = passes_hc_read_filter_set(
                &rec,
                &header,
                HcReadFilterSet::Standard,
                GATK_HC_DEFAULT_MIN_MAPPING_QUALITY,
            );
            let fail = standard_hc_read_filter_failure_index(
                &rec,
                &header,
                GATK_HC_DEFAULT_MIN_MAPPING_QUALITY,
            );
            assert_eq!(
                pass,
                fail.is_none(),
                "qname={}",
                String::from_utf8_lossy(rec.qname())
            );
        }
    }
}
