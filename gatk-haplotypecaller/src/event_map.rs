//! GATK `Event` / `EventMap` — genomic events implied by haplotype vs reference.
//! # Coordinate invariants (B5)
//! [`Event::start`] is a **0-based offset into the padded reference haplotype bytes**, not a
//! genome locus. Convert with the pad start before emit.
//! [`VariationEvent`] carries **1-based inclusive** VCF coordinates
//! (`end = start + REF.len() − 1`, GATK 4.4); use
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
use std::collections::{BTreeMap, BTreeSet, HashSet};
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
/// `end_1based = start_1based + REF.len() − 1` (insertions have `start == end`).
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

    /// GATK 4.4 / HTSJDK VCF inclusive end: `start + REF.len() − 1` (insertion ⇒ `start == end`).
    #[inline]
    pub fn vcf_end_1based(start_1based: u64, ref_allele: &str) -> u64 {
        start_1based.saturating_add(ref_allele.len().saturating_sub(1) as u64)
    }

    /// Biallelic event with Java 4.4 `VariantContext` start/end from alleles.
    pub fn from_alleles(
        contig: impl Into<String>,
        start_1based: u64,
        ref_allele: impl Into<String>,
        alt_allele: impl Into<String>,
    ) -> Self {
        let ref_allele = ref_allele.into();
        let end = Self::vcf_end_1based(start_1based, &ref_allele);
        Self {
            contig: contig.into(),
            start_1based: GenomePosition::new_1based(start_1based),
            end_1based: GenomePosition::new_1based(end),
            ref_allele,
            alt_allele: alt_allele.into(),
        }
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

/// Failure matching Java `Utils.validateArg` in `EventMap.makeBlock`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MakeBlockError(pub String);

impl fmt::Display for MakeBlockError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

fn is_simple_snp(e: &VariationEvent) -> bool {
    e.ref_allele.len() == 1 && e.alt_allele.len() == 1
}

/// HTSJDK `VariantContext.isSimpleInsertion`: biallelic indel with REF length 1.
fn is_simple_insertion(e: &VariationEvent) -> bool {
    e.is_indel() && e.ref_allele.len() == 1
}

/// HTSJDK `VariantContext.isSimpleDeletion`: biallelic indel with ALT length 1.
fn is_simple_deletion(e: &VariationEvent) -> bool {
    e.is_indel() && e.alt_allele.len() == 1
}

/// GATK 4.4 `EventMap.makeBlock(vc1, vc2)`.
///
/// `vc1` is already stored at this start; `vc2` is the newly added event
/// (`addVC` always calls `makeBlock(prev, vc)`). Encounter order matters.
pub fn make_block(
    vc1: &VariationEvent,
    vc2: &VariationEvent,
) -> Result<VariationEvent, MakeBlockError> {
    if vc1.start_1based != vc2.start_1based {
        return Err(MakeBlockError(format!(
            "vc1 and 2 must have the same start but got {} and {}",
            vc1.start_1based.get(),
            vc2.start_1based.get()
        )));
    }
    if !is_simple_snp(vc1) {
        let ok = (is_simple_deletion(vc1) && is_simple_insertion(vc2))
            || (is_simple_insertion(vc1) && is_simple_deletion(vc2));
        if !ok {
            return Err(MakeBlockError(format!(
                "Can only merge single insertion with deletion (or vice versa) but got {}→{} merging with {}→{}",
                vc1.ref_allele, vc1.alt_allele, vc2.ref_allele, vc2.alt_allele
            )));
        }
    } else if is_simple_snp(vc2) {
        return Err(MakeBlockError(format!(
            "vc1 is {}→{} but vc2 is a SNP, which implies there's been some terrible bug in the cigar {}→{}",
            vc1.ref_allele, vc1.alt_allele, vc2.ref_allele, vc2.alt_allele
        )));
    }

    if is_simple_snp(vc1) {
        if vc1.ref_allele == vc2.ref_allele {
            if vc2.alt_allele.len() < 2 {
                return Err(MakeBlockError(
                    "insertion alt must include padding base".into(),
                ));
            }
            let alt = format!("{}{}", vc1.alt_allele, &vc2.alt_allele[1..]);
            let mut out = VariationEvent::from_alleles(
                vc1.contig.as_str(),
                vc1.start_1based.get(),
                vc1.ref_allele.as_str(),
                alt,
            );
            // Java: VariantContextBuilder(vc1) keeps the SNP stop (start == end).
            out.end_1based = vc1.end_1based;
            Ok(out)
        } else {
            let mut out = VariationEvent::from_alleles(
                vc1.contig.as_str(),
                vc1.start_1based.get(),
                vc2.ref_allele.as_str(),
                vc1.alt_allele.as_str(),
            );
            out.end_1based = vc2.end_1based;
            Ok(out)
        }
    } else {
        let (insertion, deletion) = if is_simple_insertion(vc1) {
            (vc1, vc2)
        } else {
            (vc2, vc1)
        };
        let mut out = VariationEvent::from_alleles(
            vc1.contig.as_str(),
            vc1.start_1based.get(),
            deletion.ref_allele.as_str(),
            insertion.alt_allele.as_str(),
        );
        out.end_1based = deletion.end_1based;
        Ok(out)
    }
}

/// GATK 4.4 `EventMap.addVC(vc, merge=true)` over a proposed-event sequence.
///
/// Same-start events fold with [`make_block`] in **encounter order**. Output is start-sorted.
/// A pair Java would reject leaves the first event (malformed CIGAR; HC does not throw).
pub fn add_vc_merge(proposed: Vec<VariationEvent>) -> Vec<VariationEvent> {
    let mut by_start: BTreeMap<u64, VariationEvent> = BTreeMap::new();
    for vc in proposed {
        let start = vc.start_1based.get();
        match by_start.remove(&start) {
            Some(prev) => match make_block(&prev, &vc) {
                Ok(merged) => {
                    by_start.insert(start, merged);
                }
                Err(_) => {
                    by_start.insert(start, prev);
                }
            },
            None => {
                by_start.insert(start, vc);
            }
        }
    }
    by_start.into_values().collect()
}

/// GATK 4.4 `EventMap.getOverlappingEvents(loc)`.
///
/// `start <= loc <= end`. If a simple deletion ends at `loc` and a simple insertion also
/// overlaps, the deletion is dropped (insertion kept).
pub fn overlapping_events(events: &[VariationEvent], loc_1based: u64) -> Vec<VariationEvent> {
    let loc = GenomePosition::new_1based(loc_1based);
    let mut overlapping: Vec<VariationEvent> = events
        .iter()
        .filter(|e| e.start_1based <= loc && e.end_1based >= loc)
        .cloned()
        .collect();
    let del_ending: Vec<usize> = overlapping
        .iter()
        .enumerate()
        .filter(|(_, e)| is_simple_deletion(e) && e.end_1based == loc)
        .map(|(i, _)| i)
        .collect();
    let contains_insertion = overlapping.iter().any(is_simple_insertion);
    if !del_ending.is_empty() && contains_insertion {
        overlapping.remove(del_ending[0]);
    }
    overlapping
}

impl EventMap {
    /// Build events from REF/ALT CIGAR walk, then GATK 4.4 `addVC(merge=true)` / `makeBlock`
    /// (`EventMap.processCigarForInitialEvents`).
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
        let merged = add_vc_merge(proposed);
        Self {
            events: merged
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
                    (ref_loc_start_1based + e.start.get() as u64)
                        .saturating_add(e.ref_bases.len().saturating_sub(1) as u64),
                ),
                // CLONE: needed because multi-owner or ownership transfer into new structure.
                ref_allele: String::from_utf8(e.ref_bases.clone()).unwrap_or_else(|_| "N".into()),
                // CLONE: needed because owned allele string for seen-set / haplotype.
                alt_allele: String::from_utf8(e.alt_bases.clone()).unwrap_or_else(|_| "N".into()),
            })
            .collect()
    }
}

/// Search bound for **read-event discovery** (pileup / motif scans), not EventMap.
///
/// GATK 4.4.0.0 `EventMap.processCigarForInitialEvents` (SHA `2dbc0258`) has no
/// allele-length cap. A previous Rust-only filter (`len <= 40`) dropped CIGAR
/// deletions such as `171D` (REF length 172) and is not applied to CIGAR-derived
/// events (6R.47).
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
    // Java `EventMap` uses the haplotype CIGAR only (`processCigarForInitialEvents`).
    // After trim the alt is often the same length as the *trimmed* reference haplotype
    // while `ref_bytes` is still the untrimmed padded window. Re-SW vs that window
    // invents a spanning deletion (6R.47 then keeps it — no 40 bp cap) which
    // `prefer_dominant_spanning_indels` uses to drop SNPs Java still emits.
    if haplotype.bases.len() == ref_hap.bases.len() {
        return Vec::new();
    }
    // Equal-length vs the padded window: SNP/MNP — Indel-strategy SW cannot add I/D.
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
    let mut events: Vec<VariationEvent> = map.variation_events(contig, ref_loc_start_1based);
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
    let out: Vec<VariationEvent> = set.into_iter().collect();
    // Java `getAllVariantContexts` is a TreeSet union of per-haplotype EventMaps
    // (`AssemblyResultSet.regenerateVariationEvents`). It does not drop SNPs nested
    // inside another haplotype's spanning indel (6R.50), and it does not drop a SNP
    // when a *different* haplotype has an indel at the same start (6R.57).
    // `prefer_indel_over_colocated_snps` remains on per-haplotype EventMap construction.
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
        // GATK `createAlleleMapping`: SPAN_DEL is not extendable (`Allele.extend` skipped).
        if short_alt == "*" {
            return Some(Self::new(
                AlleleBytes::from_str(long_ref),
                AlleleBytes::from_str("*"),
            ));
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

/// GATK `makeMergedVariantContext` allele list at `loc` (6R.61).
///
/// Longest REF, then native longest-REF alts in encounter order, then remapped
/// shorter-REF alts (`createAlleleMapping`). Returns `None` unless a shorter REF
/// was remapped (same-REF multi-alts stay on the biallelic walk).
pub fn merged_alleles_for_genotyping(
    events: &[VariationEvent],
    loc_1based: u64,
) -> Option<(String, Vec<String>)> {
    let loc = GenomePosition::new_1based(loc_1based);
    let at_loc: Vec<&VariationEvent> = events.iter().filter(|e| e.start_1based == loc).collect();
    if at_loc.len() < 2 {
        return None;
    }
    let long_ref = at_loc
        .iter()
        .map(|e| e.ref_allele.as_str())
        .max_by_key(|r| r.len())
        .unwrap_or("")
        .to_string();
    if long_ref.is_empty() {
        return None;
    }
    if !at_loc.iter().any(|e| e.ref_allele.len() < long_ref.len()) {
        return None;
    }
    let mut alts = Vec::new();
    let mut seen: HashSet<String> = HashSet::new();
    for e in &at_loc {
        if e.alt_allele == "*" {
            continue;
        }
        if e.ref_allele == long_ref && e.alt_allele != long_ref && seen.insert(e.alt_allele.clone())
        {
            alts.push(e.alt_allele.clone());
        }
    }
    for e in &at_loc {
        if e.ref_allele == long_ref || e.alt_allele == "*" {
            continue;
        }
        if let Some(remapped) = remap_alt_onto_longer_ref(&e.ref_allele, &e.alt_allele, &long_ref) {
            if remapped != long_ref && remapped != "*" && seen.insert(remapped.clone()) {
                alts.push(remapped);
            }
        }
    }
    // Java `simpleMerge` / `createAlleleMapping`: SPAN_DEL is not extendable and remains
    // an explicit genotyping allele after `replaceSpanDels` (6R.85). Unused-ALT subset
    // may drop it before VCF emission; it is not remapped onto the longest REF.
    if at_loc.iter().any(|e| e.alt_allele == "*") && seen.insert("*".to_string()) {
        alts.push("*".to_string());
    }
    if alts.len() < 2 {
        return None;
    }
    Some((long_ref, alts))
}

/// Colocated SNP + indel after longest-REF remap (not L10 nested-STR dels, not same-REF multi-alts).
pub fn is_colocated_snp_indel_merged_site(long_ref: &str, alts: &[String]) -> bool {
    if alts.len() < 2 {
        return false;
    }
    let snp_like = alts.iter().any(|a| a.len() == long_ref.len());
    let indel = alts.iter().any(|a| a.len() != long_ref.len());
    snp_like && indel
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

    /// Pre-6R.47 Rust-only EventMap gate (must not be used on CIGAR events).
    fn legacy_allele_length_40_keeps(allele_len: usize) -> bool {
        allele_len <= MAX_VARIATION_EVENT_ALLELE_LENGTH
    }

    #[test]
    fn eventmap_short_deletion_is_emitted() {
        let ref_bytes = b"ACGTAAAA";
        let ref_hap = Haplotype::new(ref_bytes, true);
        let alt = hap_with_cigar(
            b"ACAAAA",
            &[
                (2, CigarOperator::Match),
                (2, CigarOperator::Deletion),
                (4, CigarOperator::Match),
            ],
        );
        let events = variation_events_for_haplotype(&alt, &ref_hap, ref_bytes, 100, 0, "chr");
        let del = events
            .iter()
            .find(|e| e.is_indel())
            .expect("Java EventMap emits a simple deletion");
        assert_eq!(del.ref_allele, "CGT");
        assert_eq!(del.alt_allele, "C");
        assert_eq!(del.start_1based.get(), 101);
        assert_eq!(del.end_1based.get(), 103);
    }

    #[test]
    fn eventmap_deletion_longer_than_forty_is_kept() {
        const D: usize = 50;
        let mut ref_bytes = vec![b'A'; 10];
        ref_bytes.extend(std::iter::repeat_n(b'C', D));
        ref_bytes.extend(std::iter::repeat_n(b'T', 10));
        let ref_hap = Haplotype::new(ref_bytes.clone(), true);
        let mut alt_bases = vec![b'A'; 10];
        alt_bases.extend(std::iter::repeat_n(b'T', 10));
        let alt = hap_with_cigar(
            &alt_bases,
            &[
                (10, CigarOperator::Match),
                (D, CigarOperator::Deletion),
                (10, CigarOperator::Match),
            ],
        );
        let events = variation_events_for_haplotype(&alt, &ref_hap, &ref_bytes, 1, 0, "chr");
        let del = events
            .iter()
            .find(|e| e.is_indel())
            .expect("Java EventMap has no 40 bp allele cap");
        assert_eq!(del.ref_allele.len(), D + 1);
        assert_eq!(del.alt_allele.len(), 1);
        assert!(
            !legacy_allele_length_40_keeps(del.ref_allele.len()),
            "old 40-base cap would have dropped this CIGAR deletion"
        );
    }

    /// Structure of `28M171D160M`: long D plus SNPs on the downstream match.
    #[test]
    fn eventmap_long_deletion_then_snps_emits_both() {
        const M1: usize = 28;
        const D: usize = 171;
        const M2: usize = 160;
        let mut ref_bytes = vec![b'A'; M1 + D + M2];
        let m2 = M1 + D;
        ref_bytes[m2 + 12] = b'A';
        ref_bytes[m2 + 30] = b'G';
        ref_bytes[m2 + 47] = b'T';
        let mut alt_bases = vec![b'A'; M1];
        let mut tail = vec![b'A'; M2];
        tail[12] = b'G';
        tail[30] = b'C';
        tail[47] = b'A';
        alt_bases.extend_from_slice(&tail);
        let alt = hap_with_cigar(
            &alt_bases,
            &[
                (M1, CigarOperator::Match),
                (D, CigarOperator::Deletion),
                (M2, CigarOperator::Match),
            ],
        );
        let ref_hap = Haplotype::new(ref_bytes.clone(), true);
        let events = variation_events_for_haplotype(&alt, &ref_hap, &ref_bytes, 1, 0, "chr");
        let del = events
            .iter()
            .find(|e| e.ref_allele.len() == D + 1)
            .expect("171D REF length 172 must be present");
        assert_eq!(del.alt_allele.len(), 1);
        assert!(!legacy_allele_length_40_keeps(del.ref_allele.len()));
        let snps: Vec<_> = events
            .iter()
            .filter(|e| e.is_snp())
            .map(|e| (e.ref_allele.as_str(), e.alt_allele.as_str()))
            .collect();
        assert!(snps.contains(&("A", "G")));
        assert!(snps.contains(&("G", "C")));
        assert!(snps.contains(&("T", "A")));
    }

    /// 6R.49: Java EventMap does not re-SW a trimmed SNP haplotype against the untrimmed
    /// padded reference. A Match-only CIGAR whose bases match the *reference haplotype*
    /// length must emit the SNP and must not invent a spanning deletion.
    #[test]
    fn eventmap_trimmed_snp_hap_does_not_invent_indel_against_padded_ref() {
        let mut padded = vec![b'N'; 10];
        padded.extend_from_slice(b"ACGTACGT");
        padded.extend(std::iter::repeat_n(b'N', 22));
        let mut ref_hap = hap_with_cigar(b"ACGTACGT", &[(8, CigarOperator::Match)]);
        ref_hap.is_reference = true;
        ref_hap.genome_loc = Some(GenomeLoc::new(110, 117));
        ref_hap.alignment_start_hap_wrt_ref = 10;
        let mut alt = hap_with_cigar(b"ACCTACGT", &[(8, CigarOperator::Match)]);
        alt.genome_loc = Some(GenomeLoc::new(110, 117));
        alt.alignment_start_hap_wrt_ref = 10;
        let events = collect_variation_events(&[ref_hap, alt], &padded, 100, "chr", 0);
        let snp = events
            .iter()
            .find(|e| e.is_snp())
            .unwrap_or_else(|| panic!("expected SNP from 8M CIGAR, got {events:?}"));
        assert_eq!(snp.start_1based.get(), 112);
        assert_eq!(snp.ref_allele, "G");
        assert_eq!(snp.alt_allele, "C");
        assert!(
            events
                .iter()
                .all(|e| e.ref_allele.len() <= 8 && e.alt_allele.len() <= 8),
            "must not invent a spanning deletion vs padding: {events:?}"
        );
    }

    /// 6R.59: Java `createAlleleMapping` / `simpleMerge` analogue (coordinate-free).
    #[test]
    fn colocated_snp_tc_and_deletion_tgt_remap_keeps_both_on_longest_ref() {
        let snp = VariationEvent::from_alleles("20", 1000, "T", "C");
        let del = VariationEvent::from_alleles("20", 1000, "TG", "T");
        assert_eq!(
            remap_alt_onto_longer_ref("T", "C", "TG").as_deref(),
            Some("CG")
        );
        let merged = merged_biallelic_sites_at_position(&[snp, del], 1000);
        let keys: Vec<(&str, &str)> = merged
            .iter()
            .map(|e| (e.ref_allele.as_str(), e.alt_allele.as_str()))
            .collect();
        assert!(keys.contains(&("TG", "CG")));
        assert!(keys.contains(&("TG", "T")));
    }

    #[test]
    fn merge_colocated_snp_indel_genotype_input_is_longest_ref_then_alts() {
        let snp = VariationEvent::from_alleles("20", 1000, "T", "C");
        let del = VariationEvent::from_alleles("20", 1000, "TG", "T");
        let (long_ref, alts) =
            merged_alleles_for_genotyping(&[snp, del], 1000).expect("merged site");
        assert_eq!(long_ref, "TG");
        assert_eq!(alts, vec!["T".to_string(), "CG".to_string()]);
        assert!(is_colocated_snp_indel_merged_site(&long_ref, &alts));
    }

    /// 6R.67: Java `buildEventMapsForHaplotypes(haps, fullRef, paddedRefLoc)`.
    /// Alignment start is an offset into that padded array. Rebuilding EventMaps
    /// against the trimmed apply window drops the colocated SNP+indel from the
    /// `getVariantContextsFromActiveHaplotypes` analogue, so pre-genotype merge
    /// never sees `[TG, T, CG]`.
    #[test]
    fn colocated_snp_indel_eventmaps_require_full_padded_ref_not_apply_window() {
        let mut full_ref = vec![b'A'; 10];
        full_ref.extend_from_slice(b"TG");
        full_ref.extend(std::iter::repeat_n(b'A', 12));
        let full_pad = 100u64;
        let loc = 110u64;
        let apply_off = 5usize;
        let apply_pad = full_pad + apply_off as u64;
        let apply_bases = &full_ref[apply_off..];

        let mut ref_hap = hap_with_cigar(&full_ref, &[(full_ref.len(), CigarOperator::Match)]);
        ref_hap.is_reference = true;

        let mut snp_bases = full_ref.clone();
        snp_bases[10] = b'C';
        let snp = hap_with_cigar(&snp_bases, &[(snp_bases.len(), CigarOperator::Match)]);

        let mut del_bases = full_ref.clone();
        del_bases.remove(11);
        let del = hap_with_cigar(
            &del_bases,
            &[
                (11, CigarOperator::Match),
                (1, CigarOperator::Deletion),
                (12, CigarOperator::Match),
            ],
        );

        let haps = [ref_hap, snp, del];
        let full_events = collect_variation_events(&haps, &full_ref, full_pad, "20", 0);
        let apply_events = collect_variation_events(&haps, apply_bases, apply_pad, "20", 0);

        let has_snp = |evs: &[VariationEvent]| {
            evs.iter()
                .any(|e| e.start_1based.get() == loc && e.ref_allele == "T" && e.alt_allele == "C")
        };
        let has_del = |evs: &[VariationEvent]| {
            evs.iter()
                .any(|e| e.start_1based.get() == loc && e.ref_allele == "TG" && e.alt_allele == "T")
        };
        assert!(
            has_snp(&full_events) && has_del(&full_events),
            "full padded ref EventMaps must carry T/C and TG/T at loc: {full_events:?}"
        );

        let (long_ref, alts) =
            merged_alleles_for_genotyping(&full_events, loc).expect("pre-genotype merge");
        assert_eq!(long_ref, "TG");
        assert_eq!(alts, vec!["T".to_string(), "CG".to_string()]);
        assert!(is_colocated_snp_indel_merged_site(&long_ref, &alts));

        assert!(
            !(has_snp(&apply_events) && has_del(&apply_events)),
            "trim/apply EventMap must not be Java's merge input: {apply_events:?}"
        );
        assert!(
            merged_alleles_for_genotyping(&apply_events, loc).is_none(),
            "apply-window union must not form the colocated merge: {apply_events:?}"
        );
    }

    #[test]
    fn merge_colocated_is_not_hardcoded_to_t_g() {
        let snp = VariationEvent::from_alleles("20", 2000, "A", "C");
        let del = VariationEvent::from_alleles("20", 2000, "AC", "A");
        let (long_ref, alts) =
            merged_alleles_for_genotyping(&[snp, del], 2000).expect("merged site");
        assert_eq!(long_ref, "AC");
        assert_eq!(alts, vec!["A".to_string(), "CC".to_string()]);
        assert!(is_colocated_snp_indel_merged_site(&long_ref, &alts));
    }

    #[test]
    fn nested_str_dels_are_not_snp_indel_joint_sites() {
        let short = VariationEvent::from_alleles("20", 100, "ATGTGTGTG", "A");
        let long = VariationEvent::from_alleles("20", 100, "ATGTGTGTGTGTGTGTGTG", "A");
        let (long_ref, alts) = merged_alleles_for_genotyping(&[short, long], 100).expect("remap");
        assert!(!is_colocated_snp_indel_merged_site(&long_ref, &alts));
    }

    #[test]
    fn same_ref_multi_alt_is_not_pre_genotype_merge() {
        let a = VariationEvent::from_alleles("20", 50, "G", "GTT");
        let b = VariationEvent::from_alleles("20", 50, "G", "GTTT");
        assert!(merged_alleles_for_genotyping(&[a, b], 50).is_none());
    }

    #[test]
    fn span_del_star_is_not_extended_and_is_genotyping_alt() {
        let snp = VariationEvent::from_alleles("20", 1000, "T", "C");
        let del = VariationEvent::from_alleles("20", 1000, "TG", "T");
        let star = VariationEvent::from_alleles("20", 1000, "T", "*");
        assert_eq!(
            remap_alt_onto_longer_ref("T", "*", "TG").as_deref(),
            Some("*"),
            "SPAN_DEL is not padded onto the longest REF"
        );
        let (long_ref, alts) =
            merged_alleles_for_genotyping(&[snp, del, star], 1000).expect("merged site");
        assert_eq!(long_ref, "TG");
        assert_eq!(
            alts,
            vec!["T".to_string(), "CG".to_string(), "*".to_string()],
            "Java simpleMerge keeps * as a genotyping allele after replaceSpanDels"
        );
    }

    /// Java `getAllVariantContexts` keeps a SNP from a Match haplotype even when another
    /// haplotype carries a spanning deletion over the same bases.
    #[test]
    fn eventmap_union_keeps_snp_nested_in_another_haplotype_deletion() {
        let ref_bytes = b"ACGTACGTACGTACGTACGT";
        let mut ref_hap = hap_with_cigar(ref_bytes, &[(20, CigarOperator::Match)]);
        ref_hap.is_reference = true;
        let del = hap_with_cigar(
            b"ACGTACGTACGT",
            &[
                (6, CigarOperator::Match),
                (8, CigarOperator::Deletion),
                (6, CigarOperator::Match),
            ],
        );
        let mut snp_bases = ref_bytes.to_vec();
        snp_bases[10] = b'T'; // genomic start 100+10 = 110, REF G ALT T
        let snp_hap = hap_with_cigar(&snp_bases, &[(20, CigarOperator::Match)]);
        let events = collect_variation_events(&[ref_hap, del, snp_hap], ref_bytes, 100, "chr", 0);
        let snp = events
            .iter()
            .find(|e| e.is_snp() && e.start_1based.get() == 110);
        assert!(
            snp.is_some(),
            "Java EventMap union keeps the SNP; got {events:?}"
        );
        assert!(
            events
                .iter()
                .any(|e| e.is_indel() && e.ref_allele.len() > 1),
            "deletion from the other haplotype must remain: {events:?}"
        );
    }

    /// Java `getAllVariantContexts` keeps a SNP from a Match haplotype even when another
    /// haplotype carries an insertion that starts at the same VCF coordinate.
    #[test]
    fn eventmap_union_keeps_snp_when_another_haplotype_has_insertion_at_same_start() {
        let ref_bytes = b"ACGTACGT";
        let mut ref_hap = hap_with_cigar(ref_bytes, &[(8, CigarOperator::Match)]);
        ref_hap.is_reference = true;
        let mut snp_bases = ref_bytes.to_vec();
        snp_bases[2] = b'A'; // genomic 102, REF G ALT A
        let snp_hap = hap_with_cigar(&snp_bases, &[(8, CigarOperator::Match)]);
        let ins_hap = hap_with_cigar(
            b"ACGTTTACGT",
            &[
                (3, CigarOperator::Match),
                (2, CigarOperator::Insertion),
                (5, CigarOperator::Match),
            ],
        );
        let events =
            collect_variation_events(&[ref_hap, snp_hap, ins_hap], ref_bytes, 100, "chr", 0);
        let snp = events
            .iter()
            .find(|e| e.is_snp() && e.start_1based.get() == 102);
        assert!(
            snp.is_some(),
            "Java EventMap union keeps the SNP beside a same-start insertion; got {events:?}"
        );
        assert!(
            events
                .iter()
                .any(|e| e.is_indel() && e.start_1based.get() == 102),
            "insertion from the other haplotype must remain: {events:?}"
        );
    }

    /// Length-changing Match-only CIGARs still take Indel-strategy SW (P12 A>ATG class).
    #[test]
    fn eventmap_length_changing_match_cigar_still_gets_supplemental_indel() {
        let ref_bytes = b"ACGTACGT";
        let mut ref_hap = hap_with_cigar(ref_bytes, &[(8, CigarOperator::Match)]);
        ref_hap.is_reference = true;
        let alt = hap_with_cigar(b"ACGTTACGT", &[(9, CigarOperator::Match)]);
        let events = variation_events_for_haplotype(&alt, &ref_hap, ref_bytes, 100, 0, "chr");
        assert!(
            events.iter().any(|e| e.is_indel()),
            "M-only CIGAR with extra bases must still recover an indel: {events:?}"
        );
    }

    #[test]
    fn eventmap_internal_insertion_is_emitted() {
        let ref_bytes = b"ACGTACGT";
        let ref_hap = Haplotype::new(ref_bytes, true);
        let alt = hap_with_cigar(
            b"ACGTTTACGT",
            &[
                (4, CigarOperator::Match),
                (2, CigarOperator::Insertion),
                (4, CigarOperator::Match),
            ],
        );
        let events = variation_events_for_haplotype(&alt, &ref_hap, ref_bytes, 100, 0, "chr");
        let ins = events
            .iter()
            .find(|e| e.is_indel())
            .expect("Java keeps resolved internal I");
        assert_eq!(ins.ref_allele, "T");
        assert_eq!(ins.alt_allele, "TTT");
        assert_eq!(ins.start_1based.get(), ins.end_1based.get());
    }

    #[test]
    fn eventmap_leading_insertion_is_skipped() {
        let ref_bytes = b"ACGTACGT";
        let ref_hap = Haplotype::new(ref_bytes, true);
        let alt = hap_with_cigar(
            b"TTACGTACGT",
            &[(2, CigarOperator::Insertion), (8, CigarOperator::Match)],
        );
        let events = variation_events_for_haplotype(&alt, &ref_hap, ref_bytes, 100, 0, "chr");
        assert!(
            events.iter().all(|e| !e.is_indel()),
            "Java skips I at cigarIndex==0; got {events:?}"
        );
    }
}

#[cfg(test)]
#[path = "../tests/event_map/event_map_java44.rs"]
mod event_map_java44;

#[cfg(test)]
#[path = "../tests/event_map/event_map_java44_parity_test.rs"]
mod event_map_java44_parity_test;

#[cfg(test)]
#[path = "../tests/event_map/event_map_p12_6r14_test.rs"]
mod p12_6r14_eventmap_tests;

#[cfg(test)]
#[path = "../tests/event_map/event_map_p12_6r15_test.rs"]
mod p12_6r15_audit_tests;

#[cfg(test)]
#[path = "../tests/event_map/event_map_p12_6r16_test.rs"]
mod p12_6r16_trimmed_audit_tests;

#[cfg(test)]
#[path = "../tests/event_map/event_map_p12_6r17_test.rs"]
mod p12_6r17_mapper_audit_tests;

#[cfg(test)]
#[path = "../tests/event_map/event_map_p12_6r20_test.rs"]
mod p12_6r20_mid_b_audit_tests;

#[cfg(test)]
#[path = "../tests/event_map/event_map_p12_6r21_test.rs"]
mod p12_6r21_assembly_audit_tests;

#[cfg(test)]
#[path = "../tests/event_map/event_map_p12_6r22_test.rs"]
mod p12_6r22_dangling_audit_tests;

#[cfg(test)]
#[path = "../tests/event_map/event_map_p12_6r23_test.rs"]
mod p12_6r23_threading_audit_tests;

#[cfg(test)]
#[path = "../tests/event_map/event_map_p12_6r24_test.rs"]
mod p12_6r24_oracle_provenance_tests;

#[cfg(test)]
#[path = "../tests/event_map/event_map_p12_6r25_test.rs"]
mod p12_6r25_java_ref_k25_tests;

#[cfg(test)]
#[path = "../tests/event_map/event_map_p12_6r26_test.rs"]
mod p12_6r26_java_ref_gate_production_tests;

#[cfg(test)]
#[path = "../tests/event_map/event_map_p12_6r27_test.rs"]
mod p12_6r27_call_region_none_audit_tests;

#[cfg(test)]
#[path = "../tests/event_map/event_map_p12_6r28_test.rs"]
mod p12_6r28_allele_filtering_tests;

#[cfg(test)]
#[path = "../tests/event_map/event_map_p12_6r29_test.rs"]
mod p12_6r29_extra_snp_emission_tests;

#[cfg(test)]
#[path = "../tests/event_map/event_map_p12_6r30_test.rs"]
mod p12_6r30_seqgraph_kbest_topology_tests;

#[cfg(test)]
#[path = "../tests/event_map/event_map_p12_6r31_test.rs"]
mod p12_6r31_rt_cleanup_target_tracking_tests;

#[cfg(test)]
#[path = "../tests/event_map/event_map_p12_6r32_test.rs"]
mod p12_6r32_dangling_head_parity_tests;

#[cfg(test)]
#[path = "../tests/event_map/event_map_p12_6r33_test.rs"]
mod p12_6r33_prefix_match_legacy_parity_tests;

#[cfg(test)]
#[path = "../tests/event_map/event_map_p12_6r34_test.rs"]
mod p12_6r34_path_bases_parity_tests;

#[cfg(test)]
#[path = "../tests/event_map/event_map_p12_6r35_test.rs"]
mod p12_6r35_path_bases_holdout_tests;

#[cfg(test)]
#[path = "../tests/event_map/event_map_p12_6r36_test.rs"]
mod p12_6r36_path_bases_java_parity_tests;

#[cfg(test)]
#[path = "../tests/event_map/event_map_p12_6r37_test.rs"]
mod p12_6r37_cleaned_graph_divergence_tests;

#[cfg(test)]
#[path = "../tests/event_map/event_map_p12_6r38_test.rs"]
mod p12_6r38_eventmap_vcf_parity_tests;

#[cfg(test)]
#[path = "../tests/event_map/event_map_p12_6r39_test.rs"]
mod p12_6r39_trim_max_end_parity_tests;

#[cfg(test)]
#[path = "../tests/event_map/event_map_p12_6r40_test.rs"]
mod p12_6r40_af_qual_mleac_tests;

#[cfg(test)]
#[path = "../tests/event_map/event_map_p12_6r41_test.rs"]
mod p12_6r41_qual_by_depth_jitter_tests;

#[cfg(test)]
#[path = "../tests/event_map/event_map_p12_6r48_test.rs"]
mod p12_6r48_qual_by_depth_rng_stream_tests;
