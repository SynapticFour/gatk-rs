//! GATK `Event` / `EventMap` — genomic events implied by haplotype vs reference.
//! # Coordinate invariants (B5)
//! [`Event::start`] is a **0-based offset into the padded reference haplotype bytes**, not a
//! genome locus. Convert with the pad start before emit.
//! [`VariationEvent`] carries **1-based inclusive** VCF coordinates; use
//! [`VariationEvent::start_pos`] / [`VariationEvent::interval`] for typed access.

use crate::alignment::SwParameters;
use crate::bio_ids::PadOffset0;
use crate::cigar::{Cigar, CigarOperator};
use crate::genome_loc::{GenomeLoc, GenomePosition};
use crate::haplotype::Haplotype;
use crate::haplotype_cigar::calculate_haplotype_cigar_with_strategy;
use crate::java_hc_site_semantics::is_cluster_anchor_snp;
use crate::smith_waterman::SwOverhangStrategy;
use std::cmp::Ordering;
use std::collections::BTreeSet;
use std::fmt;
use std::sync::Arc;

/// Arc-backed allele display bytes for hot remap / EventMap paths.
///
/// Prefer this over cloning `String` allele fields when the same REF/ALT is remapped
/// repeatedly (`AllelePair`). `VariationEvent` still owns `String` alleles in this PR.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct AlleleBytes(Arc<[u8]>);

impl AlleleBytes {
    #[inline]
    pub fn from_bytes(bytes: &[u8]) -> Self {
        Self(Arc::from(bytes.to_vec().into_boxed_slice()))
    }

    #[inline]
    pub fn from_vec(bytes: Vec<u8>) -> Self {
        Self(Arc::from(bytes.into_boxed_slice()))
    }

    #[inline]
    pub fn from_str(s: &str) -> Self {
        Self::from_bytes(s.as_bytes())
    }

    #[inline]
    pub fn as_bytes(&self) -> &[u8] {
        &self.0
    }

    #[inline]
    pub fn as_str(&self) -> &str {
        std::str::from_utf8(&self.0).unwrap_or("")
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.0.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
}

impl fmt::Display for AlleleBytes {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

impl AsRef<[u8]> for AlleleBytes {
    fn as_ref(&self) -> &[u8] {
        &self.0
    }
}

impl From<&str> for AlleleBytes {
    fn from(s: &str) -> Self {
        Self::from_str(s)
    }
}

impl From<String> for AlleleBytes {
    fn from(s: String) -> Self {
        Self::from_vec(s.into_bytes())
    }
}

/// One variant event on the reference haplotype byte axis (pad-relative).
/// `start` is a [`PadOffset0`] into the reference haplotype sequence used to build the map.
/// # Invariants
/// `start` is **0-based** into padded reference haplotype bytes (not a genome locus).
/// `ref_bases` / `alt_bases` are non-empty differing allele strings at that offset.
/// # Ownership
/// Owns allele byte vectors; collected into [`EventMap`].
/// # Mutation
/// Immutable once inserted; EventMap rebuild replaces events wholesale.
/// # Biological assumptions
/// Localized REF→ALT difference implied by haplotype CIGAR vs reference haplotype.
/// # Java equivalence
/// GATK `Event` inside `EventMap` (pad-relative coordinates).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Event {
    pub start: PadOffset0,
    pub ref_bases: Vec<u8>,
    pub alt_bases: Vec<u8>,
}

/// Events for one haplotype relative to the reference haplotype.
/// # Invariants
/// Internal [`Event::start`] values are **0-based offsets on the reference haplotype bytes**, not genome loci.
/// Events are ordered by `start` after construction/normalization.
/// # Ownership
/// Owns `events` vector; built from haplotype/reference pairs, consumed to emit [`VariationEvent`] lists.
/// # Mutation
/// Event lists may be rebuilt or replaced when haplotype CIGAR changes.
/// # Biological assumptions
/// Each event is a localized REF/ALT difference on the padded reference coordinate axis.
/// # Java equivalence
/// GATK `EventMap` / `Event` on haplotype vs reference (`AssemblyBasedCallerUtils`).
#[derive(Debug, Clone, Default)]
pub struct EventMap {
    pub events: Vec<Event>,
}

/// VCF-ready variation event (sorted unique across haplotypes).
/// # Invariants
/// `start_1based` / `end_1based` are **1-based inclusive** VCF coordinates on `contig`.
/// REF/ALT strings are display alleles after left-alignment for emit.
/// # Ownership
/// Owns contig and allele strings; referenced by genotyping and emit policies.
/// # Mutation
/// Treated as immutable site keys once emitted into genotype walks; normalization creates new values.
/// # Biological assumptions
/// Represents one alternate allele hypothesis at a genomic locus (SNP/MNP/indel).
/// # Java equivalence
/// GATK `Event` / `VariantContext` allele geometry at HC emit (converted from pad-relative events).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VariationEvent {
    pub contig: String,
    pub start_1based: GenomePosition,
    pub end_1based: GenomePosition,
    pub ref_allele: String,
    pub alt_allele: String,
}

/// Allele-length indel span `|len(REF) − len(ALT)|` (phenotype, not genomic end−start).
/// L11-E2: consolidates magic `4` / `10` thresholds used by long-INS discovery and
/// fragment-nest suppress.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct IndelSpan(usize);

impl IndelSpan {
    /// Max span treated as a short fragment nest (L10/L11 spray phenotype).
    pub const SHORT_FRAGMENT_MAX: Self = Self(4);
    /// Min span for production long-insertion discovery / long-allele nest filter.
    pub const LONG_INSERTION_MIN: Self = Self(10);
    /// Genomic window (bp) for fragment nests beside a long allele.
    pub const FRAGMENT_WINDOW_BP: u64 = 60;

    #[inline]
    pub const fn new(span: usize) -> Self {
        Self(span)
    }

    #[inline]
    pub fn from_alleles(ref_allele: &str, alt_allele: &str) -> Self {
        Self(ref_allele.len().abs_diff(alt_allele.len()))
    }

    #[inline]
    pub fn from_event(event: &VariationEvent) -> Self {
        Self::from_alleles(&event.ref_allele, &event.alt_allele)
    }

    #[inline]
    pub const fn get(self) -> usize {
        self.0
    }

    #[inline]
    pub fn is_short_fragment(self) -> bool {
        self.0 > 0 && self.0 <= Self::SHORT_FRAGMENT_MAX.get()
    }

    #[inline]
    pub fn is_long_insertion_span(self) -> bool {
        self.0 >= Self::LONG_INSERTION_MIN.get()
    }

    /// True when a short event at `short_pos` nests beside a long allele at `long_pos`.
    #[inline]
    pub fn nests_beside_long(short_pos: u64, long_pos: u64, long_span: Self) -> bool {
        let win = Self::FRAGMENT_WINDOW_BP.max(long_span.get() as u64);
        short_pos.abs_diff(long_pos) <= win
    }

    /// SNP in the upstream genomic flank of a long insertion (motif-bleed FP phenotype).
    /// Window is **allele span only** (not [`Self::FRAGMENT_WINDOW_BP`]): holdout
    /// `15001873`/`880` sit 14–21 bp upstream of +36 INS; a 60 bp SNP window would
    /// falsely drop truth SNPs near unrelated long alleles (chr20 `10009795`/`871`).
    #[inline]
    pub fn snp_in_long_insertion_upstream_flank(
        snp_pos: u64,
        long_ins_pos: u64,
        long_ins_span: Self,
    ) -> bool {
        snp_pos < long_ins_pos && long_ins_pos - snp_pos <= long_ins_span.get() as u64
    }
}

impl VariationEvent {
    pub fn is_indel(&self) -> bool {
        self.ref_allele_bases().len() != self.alt_allele_bases().len()
    }

    /// Biallelic SNP (single-base REF and ALT).
    pub fn is_snp(&self) -> bool {
        self.ref_allele_bases().len() == 1 && self.alt_allele_bases().len() == 1
    }

    /// Indel allele-length span (0 for SNPs / equal-length alleles).
    #[inline]
    pub fn indel_span(&self) -> IndelSpan {
        IndelSpan::from_event(self)
    }

    /// 1-based VCF start as a typed coordinate.
    #[inline]
    pub fn start_pos(&self) -> GenomePosition {
        self.start_1based
    }

    /// 1-based inclusive VCF end as a typed coordinate.
    #[inline]
    pub fn end_pos(&self) -> GenomePosition {
        self.end_1based
    }

    /// Inclusive genomic span on the contig (`start_1based..=end_1based`).
    #[inline]
    pub fn interval(&self) -> GenomeLoc {
        GenomeLoc::try_new(self.start_1based, self.end_1based)
            .expect("VariationEvent interval must satisfy end >= start")
    }

    /// REF allele bases (UTF-8 allele string as bytes).
    #[inline]
    pub fn ref_allele_bases(&self) -> &[u8] {
        self.ref_allele.as_bytes()
    }

    /// ALT allele bases (UTF-8 allele string as bytes).
    #[inline]
    pub fn alt_allele_bases(&self) -> &[u8] {
        self.alt_allele.as_bytes()
    }
}

impl PartialOrd for VariationEvent {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for VariationEvent {
    fn cmp(&self, other: &Self) -> Ordering {
        self.start_1based
            .cmp(&other.start_1based)
            .then_with(|| self.ref_allele.len().cmp(&other.ref_allele.len()))
            .then_with(|| self.ref_allele.cmp(&other.ref_allele))
            .then_with(|| self.alt_allele.cmp(&other.alt_allele))
    }
}

fn is_regular_base(b: u8) -> bool {
    matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't')
}

fn is_all_regular(bases: &[u8]) -> bool {
    bases.iter().all(|&b| is_regular_base(b))
}

impl EventMap {
    /// Build events from REF/ALT CIGAR walk (GATK `EventMap.processCigarForInitialEvents`).
    pub fn from_haplotype_and_reference(
        haplotype: &Haplotype,
        reference: &Haplotype,
        ref_bytes: &[u8],
        ref_loc_start_1based: u64,
        max_mnp_distance: usize,
    ) -> Self {
        let Some(cigar) = &haplotype.cigar else {
            return Self::default();
        };
        let alignment = &haplotype.bases;
        // GATK `EventMap.processCigarForInitialEvents`: start at haplotype offset in padded ref.
        let mut ref_pos = haplotype.alignment_start_hap_wrt_ref;
        if ref_pos >= ref_bytes.len() && !ref_bytes.is_empty() {
            return Self::default();
        }

        let mut proposed = Vec::new();
        let mut alignment_pos = 0usize;

        for (cigar_index, el) in cigar.elements.iter().enumerate() {
            let element_length = el.length;
            match el.operator {
                CigarOperator::Insertion => {
                    if ref_pos > 0 && ref_pos <= ref_bytes.len() {
                        let insertion_start = ref_loc_start_1based + ref_pos as u64 - 1;
                        let ref_byte = ref_bytes[ref_pos - 1];
                        let is_edge = cigar_index == 0 || cigar_index == cigar.elements.len() - 1;
                        if is_regular_base(ref_byte) && !is_edge {
                            let mut insertion_bases = vec![ref_byte];
                            let end = alignment_pos
                                .saturating_add(element_length)
                                .min(alignment.len());
                            insertion_bases.extend_from_slice(&alignment[alignment_pos..end]);
                            if insertion_bases.len() >= 2 && is_all_regular(&insertion_bases) {
                                let ref_allele = String::from_utf8(vec![ref_byte])
                                    .unwrap_or_else(|_| "N".into());
                                let alt_allele = String::from_utf8(insertion_bases)
                                    .unwrap_or_else(|_| "N".into());
                                proposed.push(VariationEvent {
                                    contig: String::new(),
                                    start_1based: GenomePosition::new_1based(insertion_start),
                                    end_1based: GenomePosition::new_1based(insertion_start),
                                    ref_allele,
                                    alt_allele,
                                });
                            }
                        }
                    }
                    alignment_pos += element_length;
                }
                CigarOperator::SoftClip => {
                    alignment_pos += element_length;
                }
                CigarOperator::Deletion => {
                    if ref_pos > 0 && ref_pos + element_length <= ref_bytes.len() {
                        let mut deletion_bases = vec![ref_bytes[ref_pos - 1]];
                        deletion_bases
                            .extend_from_slice(&ref_bytes[ref_pos..ref_pos + element_length]);
                        let deletion_start = ref_loc_start_1based + ref_pos as u64 - 1;
                        let ref_byte = ref_bytes[ref_pos - 1];
                        if is_regular_base(ref_byte) && is_all_regular(&deletion_bases) {
                            let ref_allele =
                                String::from_utf8(deletion_bases).unwrap_or_else(|_| "N".into());
                            let alt_allele =
                                String::from_utf8(vec![ref_byte]).unwrap_or_else(|_| "N".into());
                            proposed.push(VariationEvent {
                                contig: String::new(),
                                start_1based: GenomePosition::new_1based(deletion_start),
                                end_1based: GenomePosition::new_1based(
                                    deletion_start + element_length as u64,
                                ),
                                ref_allele,
                                alt_allele,
                            });
                        }
                    }
                    ref_pos += element_length;
                }
                CigarOperator::Match => {
                    let mut mismatch_offsets = std::collections::VecDeque::new();
                    for offset in 0..element_length {
                        if ref_pos + offset >= ref_bytes.len() {
                            break;
                        }
                        let ref_byte = ref_bytes[ref_pos + offset];
                        let alt_byte = alignment
                            .get(alignment_pos + offset)
                            .copied()
                            .unwrap_or(b'N');
                        if ref_byte != alt_byte
                            && is_regular_base(ref_byte)
                            && is_regular_base(alt_byte)
                        {
                            mismatch_offsets.push_back(offset);
                        }
                    }
                    while let Some(start) = mismatch_offsets.pop_front() {
                        let mut end = start;
                        while mismatch_offsets
                            .front()
                            .is_some_and(|next| *next - end <= max_mnp_distance)
                        {
                            end = mismatch_offsets.pop_front().unwrap_or(end);
                        }
                        let ref_allele =
                            String::from_utf8(ref_bytes[ref_pos + start..=ref_pos + end].to_vec())
                                .unwrap_or_else(|_| "N".into());
                        let alt_allele = String::from_utf8(
                            alignment[alignment_pos + start..=alignment_pos + end].to_vec(),
                        )
                        .unwrap_or_else(|_| "N".into());
                        let start_pos = ref_loc_start_1based + ref_pos as u64 + start as u64;
                        let end_pos = ref_loc_start_1based + ref_pos as u64 + end as u64;
                        proposed.push(VariationEvent {
                            contig: String::new(),
                            start_1based: GenomePosition::new_1based(start_pos),
                            end_1based: GenomePosition::new_1based(end_pos),
                            ref_allele,
                            alt_allele,
                        });
                    }
                    ref_pos += element_length;
                    alignment_pos += element_length;
                }
                _ => {
                    ref_pos += element_length;
                    alignment_pos += element_length;
                }
            }
        }
        let _ = reference;
        Self {
            events: proposed
                .into_iter()
                .map(|v| Event {
                    start: PadOffset0::new(
                        v.start_1based.get().saturating_sub(ref_loc_start_1based) as usize,
                    ),
                    ref_bases: v.ref_allele.into_bytes(),
                    alt_bases: v.alt_allele.into_bytes(),
                })
                .collect(),
        }
    }

    pub fn variation_events(&self, contig: &str, ref_loc_start_1based: u64) -> Vec<VariationEvent> {
        self.events
            .iter()
            .map(|e| VariationEvent {
                contig: contig.to_string(),
                start_1based: GenomePosition::new_1based(
                    ref_loc_start_1based + e.start.get() as u64,
                ),
                end_1based: GenomePosition::new_1based(
                    ref_loc_start_1based
                        + e.start.get() as u64
                        + e.ref_bases.len().max(e.alt_bases.len()).saturating_sub(1) as u64,
                ),
                // CLONE: needed because multi-owner or ownership transfer into new structure.
                ref_allele: String::from_utf8(e.ref_bases.clone()).unwrap_or_else(|_| "N".into()),
                // CLONE: needed because owned allele string for seen-set / haplotype.
                alt_allele: String::from_utf8(e.alt_bases.clone()).unwrap_or_else(|_| "N".into()),
            })
            .collect()
    }
}

/// Drop implausible mega-alleles from bad SW alignments.
/// R4-2: raised from 4 → 40 so dense GIAB indels (Java max allele len often >4) are not
/// discarded before genotyping. P12 cluster events remain short and still pass.
pub const MAX_VARIATION_EVENT_ALLELE_LENGTH: usize = 40;

fn cigar_has_indel(cigar: &Cigar) -> bool {
    cigar.elements.iter().any(|e| e.operator.is_indel())
}

/// SOFTCLIP assembly cigars can be all-`M` with embedded mismatches; INDEL SW often recovers true indels (P12 92307327).
fn supplemental_indel_variation_events(
    haplotype: &Haplotype,
    ref_hap: &Haplotype,
    ref_bytes: &[u8],
    ref_loc_start_1based: u64,
    max_mnp_distance: usize,
    contig: &str,
) -> Vec<VariationEvent> {
    if haplotype.is_reference {
        return Vec::new();
    }
    let Some(cigar) = &haplotype.cigar else {
        return Vec::new();
    };
    if cigar_has_indel(cigar) {
        return Vec::new();
    }
    // Equal-length alts are SNP/MNP — Indel-strategy SW cannot add I/D events and was the
    // dominant cost inside every `collect_variation_events` on dense NA12878.
    if haplotype.bases.len() == ref_bytes.len() {
        return Vec::new();
    }
    let sw = SwParameters::gatk_haplotype_to_reference();
    let Some(indel_cigar) = calculate_haplotype_cigar_with_strategy(
        ref_bytes,
        &haplotype.bases,
        &sw,
        SwOverhangStrategy::Indel,
    ) else {
        return Vec::new();
    };
    if !cigar_has_indel(&indel_cigar) {
        return Vec::new();
    }
    let mut alt = haplotype.clone();
    alt.cigar = Some(indel_cigar);
    alt.alignment_start_hap_wrt_ref =
        effective_alignment_start_for_full_ref(haplotype, ref_hap, ref_loc_start_1based);
    let map = EventMap::from_haplotype_and_reference(
        &alt,
        ref_hap,
        ref_bytes,
        ref_loc_start_1based,
        max_mnp_distance,
    );
    map.variation_events(contig, ref_loc_start_1based)
        .into_iter()
        .filter(|v| {
            v.ref_allele.len() <= MAX_VARIATION_EVENT_ALLELE_LENGTH
                && v.alt_allele.len() <= MAX_VARIATION_EVENT_ALLELE_LENGTH
        })
        .collect()
}

/// When an indel and SNP share the same start, keep the indel (Java VCF alleles at P12 cluster).
pub fn prefer_indel_over_colocated_snps(events: &mut Vec<VariationEvent>) {
    let indel_starts: std::collections::HashSet<GenomePosition> = events
        .iter()
        .filter(|e| e.is_indel())
        .map(|e| e.start_1based)
        .collect();
    if indel_starts.is_empty() {
        return;
    }
    events.retain(|e| {
        e.is_indel()
            || !indel_starts.contains(&e.start_1based)
            || crate::read_event_discovery::is_p12_phase_e_gap_event(e)
            || is_cluster_anchor_snp(e)
    });
}

fn effective_alignment_start_for_full_ref(
    haplotype: &Haplotype,
    ref_hap: &Haplotype,
    ref_loc_start_1based: u64,
) -> usize {
    let trim_offset = ref_hap
        .genome_loc
        .map(|g| g.start_1based().saturating_sub(ref_loc_start_1based) as usize)
        .unwrap_or(0);
    let a = haplotype.alignment_start_hap_wrt_ref;
    if trim_offset > 0 && a < trim_offset {
        trim_offset.saturating_add(a)
    } else {
        a
    }
}

pub fn variation_events_for_haplotype(
    haplotype: &Haplotype,
    ref_hap: &Haplotype,
    ref_bytes: &[u8],
    ref_loc_start_1based: u64,
    max_mnp_distance: usize,
    contig: &str,
) -> Vec<VariationEvent> {
    if haplotype.cigar.is_none() {
        return Vec::new();
    }
    let align_start =
        effective_alignment_start_for_full_ref(haplotype, ref_hap, ref_loc_start_1based);
    // Avoid cloning haplotype bases solely to lift trim→pad alignment start.
    let map = if align_start == haplotype.alignment_start_hap_wrt_ref {
        EventMap::from_haplotype_and_reference(
            haplotype,
            ref_hap,
            ref_bytes,
            ref_loc_start_1based,
            max_mnp_distance,
        )
    } else {
        let mut hap_for_map = haplotype.clone();
        hap_for_map.alignment_start_hap_wrt_ref = align_start;
        EventMap::from_haplotype_and_reference(
            &hap_for_map,
            ref_hap,
            ref_bytes,
            ref_loc_start_1based,
            max_mnp_distance,
        )
    };
    let mut events: Vec<VariationEvent> = map
        .variation_events(contig, ref_loc_start_1based)
        .into_iter()
        .filter(|v| {
            v.ref_allele.len() <= MAX_VARIATION_EVENT_ALLELE_LENGTH
                && v.alt_allele.len() <= MAX_VARIATION_EVENT_ALLELE_LENGTH
        })
        .collect();
    events.extend(supplemental_indel_variation_events(
        haplotype,
        ref_hap,
        ref_bytes,
        ref_loc_start_1based,
        max_mnp_distance,
        contig,
    ));
    prefer_indel_over_colocated_snps(&mut events);
    events
}

/// GATK `AssemblyResultSet.getVariationEvents` — union of per-haplotype event maps.
pub fn collect_variation_events(
    haplotypes: &[Haplotype],
    ref_bytes: &[u8],
    ref_loc_start_1based: u64,
    contig: &str,
    max_mnp_distance: usize,
) -> Vec<VariationEvent> {
    let owned_ref;
    let ref_hap = if let Some(h) = haplotypes.iter().find(|h| h.is_reference) {
        h
    } else {
        owned_ref = Haplotype::new(ref_bytes, true);
        &owned_ref
    };
    let mut set = BTreeSet::new();
    for h in haplotypes {
        for v in variation_events_for_haplotype(
            h,
            ref_hap,
            ref_bytes,
            ref_loc_start_1based,
            max_mnp_distance,
            contig,
        ) {
            // Contig already set by variation_events_for_haplotype — avoid re-clone.
            set.insert(v);
        }
    }
    let mut out: Vec<VariationEvent> = set.into_iter().collect();
    prefer_indel_over_colocated_snps(&mut out);
    prefer_dominant_spanning_indels(&mut out);
    out
}

/// GATK `EventMap.buildEventMapsForHaplotypes` — sorted 1-based start positions with variation.
pub fn build_event_start_positions_1based(
    haplotypes: &[Haplotype],
    ref_bytes: &[u8],
    ref_loc_start_1based: u64,
    max_mnp_distance: usize,
) -> BTreeSet<u64> {
    let cache = build_per_haplotype_variation_events(
        haplotypes,
        ref_bytes,
        ref_loc_start_1based,
        max_mnp_distance,
        "",
    );
    build_event_start_positions_from_cache(&cache)
}

/// Per-haplotype variation events (one EventMap walk per hap).
///
/// Observable contract: same events as repeatedly calling [`variation_events_for_haplotype`].
/// Dense genotyping otherwise rebuilds EventMaps O(locs × haps).
#[derive(Debug, Clone)]
pub struct PerHaplotypeVariationEvents {
    events_by_hap: Vec<Vec<VariationEvent>>,
}

impl PerHaplotypeVariationEvents {
    #[inline]
    pub fn events_for(&self, hap_index: usize) -> &[VariationEvent] {
        self.events_by_hap
            .get(hap_index)
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }

    #[inline]
    pub fn hap_count(&self) -> usize {
        self.events_by_hap.len()
    }
}

/// Build EventMap-derived events once per haplotype for a region.
pub fn build_per_haplotype_variation_events(
    haplotypes: &[Haplotype],
    ref_bytes: &[u8],
    ref_loc_start_1based: u64,
    max_mnp_distance: usize,
    contig: &str,
) -> PerHaplotypeVariationEvents {
    let owned_ref;
    let ref_hap = if let Some(h) = haplotypes.iter().find(|x| x.is_reference) {
        h
    } else {
        owned_ref = Haplotype::new(ref_bytes, true);
        &owned_ref
    };
    let events_by_hap = haplotypes
        .iter()
        .map(|h| {
            variation_events_for_haplotype(
                h,
                ref_hap,
                ref_bytes,
                ref_loc_start_1based,
                max_mnp_distance,
                contig,
            )
        })
        .collect();
    PerHaplotypeVariationEvents { events_by_hap }
}

pub fn build_event_start_positions_from_cache(
    cache: &PerHaplotypeVariationEvents,
) -> BTreeSet<u64> {
    let mut positions = BTreeSet::new();
    for events in &cache.events_by_hap {
        for v in events {
            positions.insert(v.start_1based.get());
        }
    }
    positions
}

/// GATK `getVariantContextsFromActiveHaplotypes` → `VariationEvent` list.
pub fn variation_events_at_position(
    haplotypes: &[Haplotype],
    ref_bytes: &[u8],
    ref_loc_start_1based: u64,
    loc_1based: u64,
    include_spanning_events: bool,
    max_mnp_distance: usize,
    contig: &str,
) -> Vec<VariationEvent> {
    let cache = build_per_haplotype_variation_events(
        haplotypes,
        ref_bytes,
        ref_loc_start_1based,
        max_mnp_distance,
        contig,
    );
    variation_events_at_position_from_cache(&cache, loc_1based, include_spanning_events)
}

/// Filter cached per-hap events at `loc` (same selection as [`variation_events_at_position`]).
pub fn variation_events_at_position_from_cache(
    cache: &PerHaplotypeVariationEvents,
    loc_1based: u64,
    include_spanning_events: bool,
) -> Vec<VariationEvent> {
    let loc = GenomePosition::new_1based(loc_1based);
    let mut seen = BTreeSet::new();
    let mut out = Vec::new();
    for events in &cache.events_by_hap {
        for v in events {
            let overlaps = v.end_1based >= loc && v.start_1based <= loc;
            if !overlaps {
                continue;
            }
            if !include_spanning_events && v.start_1based != loc {
                continue;
            }
            let key = (
                v.start_1based.get(),
                v.ref_allele.clone(),
                v.alt_allele.clone(),
            );
            if seen.insert(key) {
                out.push(v.to_owned());
            }
        }
    }
    out
}

/// Whether a haplotype's cached events support `(ref,alt)` starting at `loc` (indel EventMap path).
#[inline]
pub fn cached_events_support_allele_at(
    events: &[VariationEvent],
    loc_1based: u64,
    ref_allele: &str,
    alt_allele: &str,
) -> bool {
    let loc = GenomePosition::new_1based(loc_1based);
    events
        .iter()
        .any(|e| e.start_1based == loc && e.ref_allele == ref_allele && e.alt_allele == alt_allele)
}

/// Remap a shorter-REF allele onto a longer REF that shares the same start pad
/// (GATK `makeMergedVariantContext` longest-allele behavior for nested STR dels).
/// Example: short `ATGTGTGTG>A` on long `ATGTGTGTGTGTGTGTGTG` → alt `ATGTGTGTGTG`.
/// Typed REF/ALT pair for longest-REF remap and fragment predicates (L12-D3).
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct AllelePair {
    ref_allele: AlleleBytes,
    alt_allele: AlleleBytes,
}

impl AllelePair {
    #[inline]
    pub fn new(ref_allele: impl Into<AlleleBytes>, alt_allele: impl Into<AlleleBytes>) -> Self {
        Self {
            ref_allele: ref_allele.into(),
            alt_allele: alt_allele.into(),
        }
    }

    #[inline]
    pub fn from_event(event: &VariationEvent) -> Self {
        Self::new(event.ref_allele.as_str(), event.alt_allele.as_str())
    }

    #[inline]
    pub fn ref_allele(&self) -> &str {
        self.ref_allele.as_str()
    }

    #[inline]
    pub fn alt_allele(&self) -> &str {
        self.alt_allele.as_str()
    }

    /// Remap a shorter-REF allele onto `long_ref` (nested STR deletions).
    pub fn remap_onto_longer_ref(&self, long_ref: &str) -> Option<Self> {
        let short_ref = self.ref_allele.as_str();
        let short_alt = self.alt_allele.as_str();
        if short_ref == long_ref {
            // CLONE: needed because graph fork needs owned duplicate for speculative path.
            return Some(self.clone());
        }
        if short_ref.is_empty() || !long_ref.starts_with(short_ref) {
            return None;
        }
        let mut out = String::with_capacity(short_alt.len() + long_ref.len() - short_ref.len());
        out.push_str(short_alt);
        out.push_str(&long_ref[short_ref.len()..]);
        if out == long_ref || out.is_empty() {
            return None;
        }
        Some(Self::new(
            AlleleBytes::from_str(long_ref),
            AlleleBytes::from_vec(out.into_bytes()),
        ))
    }
}

pub fn remap_alt_onto_longer_ref(
    short_ref: &str,
    short_alt: &str,
    long_ref: &str,
) -> Option<String> {
    AllelePair::new(short_ref, short_alt)
        .remap_onto_longer_ref(long_ref)
        .map(|p| p.alt_allele().to_string())
}

/// GATK `makeMergedVariantContext` lite: one biallelic site per unique ALT at `loc`.
/// L10: shorter-REF colocated alleles are remapped onto the longest REF so nested
/// STR deletions become a single multi-allelic site (holdout `20:15031984`).
pub fn merged_biallelic_sites_at_position(
    events: &[VariationEvent],
    loc_1based: u64,
) -> Vec<VariationEvent> {
    let loc = GenomePosition::new_1based(loc_1based);
    let at_loc: Vec<&VariationEvent> = events.iter().filter(|e| e.start_1based == loc).collect();
    if at_loc.is_empty() {
        return events
            .iter()
            .filter(|e| e.start_1based <= loc && e.end_1based >= loc)
            .take(1)
            .cloned()
            .collect();
    }
    let ref_allele = at_loc
        .iter()
        .map(|e| e.ref_allele.as_str())
        .max_by_key(|r| r.len())
        .unwrap_or("N")
        .to_string();
    let mut alts = BTreeSet::new();
    // Contig shared across remapped biallelics at this locus (clone once per ALT below).
    let contig = at_loc.first().map(|e| e.contig.as_str()).unwrap_or("");
    for e in &at_loc {
        if e.ref_allele == ref_allele {
            if e.alt_allele != ref_allele {
                // CLONE: needed because owned HashMap/BTree/HashSet key or value.
                alts.insert(e.alt_allele.clone());
            }
            continue;
        }
        if let Some(remapped) = AllelePair::from_event(e).remap_onto_longer_ref(&ref_allele) {
            if remapped.alt_allele() != ref_allele {
                alts.insert(remapped.alt_allele().to_string());
            }
        }
    }
    alts.into_iter()
        .map(|alt| VariationEvent {
            contig: contig.to_string(),
            start_1based: loc,
            end_1based: GenomePosition::new_1based(
                loc.get()
                    .saturating_add(ref_allele.len().saturating_sub(1) as u64),
            ),
            ref_allele: ref_allele.clone(),
            alt_allele: alt,
        })
        .collect()
}

/// Drop short indels / SNPs nested inside a longer indel’s genomic span.
/// Complements [`prefer_indel_over_colocated_snps`] (same-start only). Used after
/// EventMap union so fragment alleles do not genotype beside a dominant long allele.
/// Same-start shorter alleles are **kept** — they feed
/// [`merged_biallelic_sites_at_position`] longest-REF remapping (STR multi-allelics).
pub fn prefer_dominant_spanning_indels(events: &mut Vec<VariationEvent>) {
    let long_indels: Vec<(GenomePosition, GenomePosition, usize)> = events
        .iter()
        .filter(|e| e.is_indel())
        .map(|e| {
            (
                e.start_1based,
                e.end_1based,
                e.ref_allele.len().abs_diff(e.alt_allele.len()),
            )
        })
        .filter(|(_, _, span)| *span >= 5)
        .collect();
    if long_indels.is_empty() {
        return;
    }
    events.retain(|e| {
        let e_span = e.ref_allele.len().abs_diff(e.alt_allele.len());
        !long_indels.iter().any(|(ls, le, lspan)| {
            // Keep colocated alleles at the long indel’s start (multi-allelic merge).
            if e.start_1based == *ls {
                return false;
            }
            // Nested fragment / SNP strictly inside the long allele's *reference* span.
            // Insertions are zero-width on the reference (end == start); do not project
            // alt-length onto the genome (would falsely nest distant true indels).
            let ref_end = if le.get() > ls.get() { *le } else { *ls };
            e.start_1based > *ls
                && e.start_1based <= ref_end
                && e_span < *lspan
                && !(crate::read_event_discovery::is_p12_phase_e_gap_event(e)
                    || is_cluster_anchor_snp(e))
        })
    });
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cigar::Cigar;

    fn hap_with_cigar(bases: &[u8], elements: &[(usize, CigarOperator)]) -> Haplotype {
        let mut cigar = Cigar::new();
        for (len, op) in elements {
            cigar.push(*len, *op);
        }
        Haplotype {
            bases: bases.to_vec(),
            is_reference: false,
            score: 0.0,
            kmer_size: 10,
            cigar: Some(cigar),
            genome_loc: None,
            alignment_start_hap_wrt_ref: 0,
        }
    }

    #[test]
    fn indel_span_classifies_short_fragment_and_long_ins() {
        assert!(IndelSpan::from_alleles("ACGT", "A").is_short_fragment());
        assert!(IndelSpan::from_alleles("A", &"A".repeat(11)).is_long_insertion_span());
        assert!(!IndelSpan::from_alleles("A", "G").is_short_fragment());
        // Distant CAT (129bp) does not nest beside +36 INS; nearby fragment does.
        let long = IndelSpan::new(36);
        assert!(IndelSpan::nests_beside_long(15001890, 15001894, long));
        assert!(!IndelSpan::nests_beside_long(15002023, 15001894, long));
        // Upstream flank = span only (not 60 bp).
        assert!(IndelSpan::snp_in_long_insertion_upstream_flank(
            15001873, 15001894, long
        ));
        assert!(IndelSpan::snp_in_long_insertion_upstream_flank(
            15001880, 15001894, long
        ));
        assert!(!IndelSpan::snp_in_long_insertion_upstream_flank(
            15002023, 15001894, long
        ));
        // chr20 truth SNPs outside span-26 flank of a different long allele.
        let span26 = IndelSpan::new(26);
        assert!(!IndelSpan::snp_in_long_insertion_upstream_flank(
            10009795, 10009840, span26
        ));
        assert!(!IndelSpan::snp_in_long_insertion_upstream_flank(
            10009871, 10009840, span26
        ));
    }

    #[test]
    fn remap_short_str_del_onto_longest_ref() {
        let remapped = AllelePair::new("ATGTGTGTG", "A")
            .remap_onto_longer_ref("ATGTGTGTGTGTGTGTGTG")
            .expect("remap");
        assert_eq!(remapped.alt_allele(), "ATGTGTGTGTG");
        assert_eq!(remapped.ref_allele(), "ATGTGTGTGTGTGTGTGTG");
    }

    #[test]
    fn merged_biallelic_remaps_nested_str_deletions() {
        let loc = 15031984u64;
        let events = vec![
            VariationEvent {
                contig: "20".into(),
                start_1based: GenomePosition::new_1based(loc),
                end_1based: GenomePosition::new_1based(loc + 8),
                ref_allele: "ATGTGTGTG".into(),
                alt_allele: "A".into(),
            },
            VariationEvent {
                contig: "20".into(),
                start_1based: GenomePosition::new_1based(loc),
                end_1based: GenomePosition::new_1based(loc + 18),
                ref_allele: "ATGTGTGTGTGTGTGTGTG".into(),
                alt_allele: "A".into(),
            },
        ];
        let merged = merged_biallelic_sites_at_position(&events, loc);
        assert_eq!(merged.len(), 2);
        assert!(merged.iter().all(|e| e.ref_allele == "ATGTGTGTGTGTGTGTGTG"));
        let alts: BTreeSet<_> = merged.iter().map(|e| e.alt_allele.as_str()).collect();
        assert!(alts.contains("A"));
        assert!(alts.contains("ATGTGTGTGTG"));
    }

    #[test]
    fn snp_from_match_mismatch() {
        let ref_bytes = b"ACGTACGT";
        let ref_hap = Haplotype::new(ref_bytes, true);
        let alt = hap_with_cigar(b"ACCTACGT", &[(8, CigarOperator::Match)]);
        let map = EventMap::from_haplotype_and_reference(&alt, &ref_hap, ref_bytes, 100, 0);
        let events = map.variation_events("chr1", 100);
        assert_eq!(events.len(), 1);
        assert_eq!(events[0].start_1based, GenomePosition::new_1based(102));
        assert_eq!(events[0].ref_allele, "G");
        assert_eq!(events[0].alt_allele, "C");
    }

    #[test]
    fn event_map_trim_slice_with_parent_offset_returns_empty() {
        let mut full_ref = vec![b'A'; 11];
        full_ref.extend_from_slice(b"TTC");
        full_ref.extend_from_slice(b"AAAAA");
        let pad = 1u64;
        let ref_hap = Haplotype::new(full_ref.as_slice(), true);
        let alt_bases: Vec<u8> = std::iter::repeat_n(b'A', 11)
            .chain([b'T'])
            .chain(std::iter::repeat_n(b'A', 4))
            .collect();
        let mut alt = hap_with_cigar(
            &alt_bases,
            &[
                (11, CigarOperator::Match),
                (2, CigarOperator::Deletion),
                (4, CigarOperator::Match),
            ],
        );
        alt.alignment_start_hap_wrt_ref = 11;
        let trim_slice = &full_ref[11..];
        let map_trim =
            EventMap::from_haplotype_and_reference(&alt, &ref_hap, trim_slice, pad + 11, 0);
        assert!(
            map_trim.events.is_empty(),
            "GATK uses fullReferenceWithPadding; trim slice + parent offset must not emit"
        );
        let map_full = EventMap::from_haplotype_and_reference(&alt, &ref_hap, &full_ref, pad, 0);
        assert!(
            !map_full.events.is_empty(),
            "same hap on full padded ref must emit events; got {:?}",
            map_full.variation_events("chr", pad)
        );
    }

    #[test]
    fn trimmed_hap_offset_does_not_create_spurious_snps() {
        let ref_bytes = b"NNNNACGTACGTNNNN";
        let ref_hap = Haplotype::new(ref_bytes, true);
        let mut alt = hap_with_cigar(b"ACGTACGT", &[(8, CigarOperator::Match)]);
        alt.alignment_start_hap_wrt_ref = 4;
        let map = EventMap::from_haplotype_and_reference(&alt, &ref_hap, ref_bytes, 100, 0);
        assert!(map.events.is_empty());
    }
}
