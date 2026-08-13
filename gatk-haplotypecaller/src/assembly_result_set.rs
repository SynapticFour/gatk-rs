//! GATK `AssemblyResultSet` — multi-kmer assembly outcomes for one region (E2E.4).

use crate::assembly_region_iterator::AssemblyRegion;
use crate::event_map::{collect_variation_events, VariationEvent};
use crate::genome_loc::GenomeLoc;
use crate::haplotype::{haplotype_size_and_base_order, Haplotype};
use crate::read_threading_assembler::{AssemblyResult, AssemblyStatus};
use gatk_common::GatkResult;
use std::collections::{BTreeSet, HashMap};
use std::sync::Arc;

/// Default GATK HC `--max-mnp-distance` (HaplotypeCallerArgumentCollection).
pub const DEFAULT_MAX_MNP_DISTANCE: usize = 0;

/// Rust mirror of GATK `AssemblyResultSet` (haplotype list + variation / kmer metadata).
/// # Invariants
/// `variation_present` stays false until calling-path event regeneration ([`Self::from_assembly_for_calling`]).
/// Haplotypes share the padded reference coordinate system for EventMap/CIGAR.
/// # Ownership
/// Owns haplotypes, contig, and variation events. Padded reference bases are
/// [`Arc<[u8]>`] — immutable and cheap to share across trim / apply / genotyping stages.
/// # Mutation
/// Calling path regenerates variation events and may rewrite haplotype metadata in place.
/// # Biological assumptions
/// Multi-kmer assembly outcomes for one active region before genotyping.
/// # Java equivalence
/// GATK `AssemblyResultSet` (`regenerateVariationEvents`, haplotype list).
#[derive(Debug, Clone)]
pub struct AssemblyResultSet {
    pub haplotypes: Vec<Haplotype>,
    pub(crate) variation_present: bool,
    kmer_sizes: BTreeSet<usize>,
    padded_reference_start_1based: u64,
    /// Immutable padded reference window (shared; not mutated after construction).
    reference_bases: Arc<[u8]>,
    pub(crate) variation_events: Vec<VariationEvent>,
    pub(crate) contig: String,
    max_mnp_distance: usize,
}

impl AssemblyResultSet {
    /// Build from a single `assemble_from_ref_and_reads` outcome (parity gate metadata only).
    /// Java `HcFullParityGateDump` reads `isVariationPresent` before event-map regeneration;
    /// keep `variation_present` false until [`Self::from_assembly_for_calling`].
    pub fn from_assembly_result(result: &AssemblyResult) -> Self {
        let _ = result.status;
        let _ = AssemblyStatus::AssembledSomeVariation;
        Self {
            haplotypes: result.haplotypes.clone(),
            variation_present: false,
            kmer_sizes: BTreeSet::new(),
            padded_reference_start_1based: 0,
            reference_bases: Arc::from([]),
            variation_events: Vec::new(),
            contig: String::new(),
            max_mnp_distance: DEFAULT_MAX_MNP_DISTANCE,
        }
    }

    /// Production path: padded reference + variation events (GATK `regenerateVariationEvents`).
    pub fn from_assembly_for_calling(
        result: &AssemblyResult,
        reference_bases: impl Into<Arc<[u8]>>,
        padded_reference_start_1based: u64,
        contig: &str,
        max_mnp_distance: usize,
    ) -> Self {
        // CLONE: shared API keeps `AssemblyResult` for dump callers; production prefers
        // [`Self::from_assembly_for_calling_owned`].
        Self::from_assembly_for_calling_owned(
            result.status,
            result.kmer_size,
            result.haplotypes.clone(),
            reference_bases,
            padded_reference_start_1based,
            contig,
            max_mnp_distance,
        )
    }

    /// Like [`Self::from_assembly_for_calling`] but takes haplotypes by move (no deep clone).
    pub fn from_assembly_for_calling_owned(
        status: AssemblyStatus,
        kmer_size: usize,
        mut haplotypes: Vec<Haplotype>,
        reference_bases: impl Into<Arc<[u8]>>,
        padded_reference_start_1based: u64,
        contig: &str,
        max_mnp_distance: usize,
    ) -> Self {
        let reference_bases: Arc<[u8]> = reference_bases.into();
        let variation_present = matches!(status, AssemblyStatus::AssembledSomeVariation)
            || (haplotypes.iter().any(|h| !h.is_reference) && haplotypes.len() > 1);
        let mut kmer_sizes = BTreeSet::new();
        if kmer_size > 0 {
            kmer_sizes.insert(kmer_size);
        }
        Haplotype::tag_padded_reference_span(&mut haplotypes, padded_reference_start_1based);
        let variation_events = if variation_present {
            collect_variation_events(
                &haplotypes,
                reference_bases.as_ref(),
                padded_reference_start_1based,
                contig,
                max_mnp_distance,
            )
        } else {
            Vec::new()
        };
        Self {
            haplotypes,
            variation_present,
            kmer_sizes,
            padded_reference_start_1based,
            reference_bases,
            variation_events,
            contig: contig.to_string(),
            max_mnp_distance,
        }
    }

    pub fn is_variation_present(&self) -> bool {
        self.variation_present && self.haplotypes.len() > 1
    }

    /// Whether `call_region` should proceed (Java trimmer `isVariationPresent` or events on assembly).
    pub fn has_variation_for_calling(&self) -> bool {
        !self.variation_events.is_empty() || self.is_variation_present()
    }

    pub fn variation_events(&self) -> &[VariationEvent] {
        &self.variation_events
    }

    pub fn padded_reference_start_1based(&self) -> u64 {
        self.padded_reference_start_1based
    }

    pub fn reference_bases(&self) -> &[u8] {
        self.reference_bases.as_ref()
    }

    /// Shared handle to the padded reference (refcount bump only — no byte copy).
    pub fn reference_bases_shared(&self) -> Arc<[u8]> {
        Arc::clone(&self.reference_bases)
    }

    /// Reference-haplotype apply window as shared bytes.
    /// Prefer the assembly [`Arc`] when the ref haplotype bases match; otherwise allocate a
    /// fresh [`Arc`] from the haplotype (still shared for later stages).
    pub fn apply_bases_shared(&self) -> Arc<[u8]> {
        if let Some(h) = self.haplotypes.iter().find(|h| h.is_reference) {
            if !h.bases.is_empty() {
                if h.bases.as_slice() == self.reference_bases.as_ref() {
                    return self.reference_bases_shared();
                }
                return Arc::<[u8]>::from(h.bases.as_slice());
            }
        }
        self.reference_bases_shared()
    }

    /// GATK `EventMap.buildEventMapsForHaplotypes`: full padded reference + start (not trim slice).
    pub fn event_map_reference(&self) -> (&[u8], u64) {
        (
            self.reference_bases.as_ref(),
            self.padded_reference_start_1based,
        )
    }

    /// GATK `getMinimumKmerSize` — `None` when no kmer was registered (dump uses 0).
    pub fn minimum_kmer_size(&self) -> Option<usize> {
        self.kmer_sizes.iter().copied().next()
    }

    #[doc(hidden)]
    pub fn status_for_dump(&self) -> &'static str {
        if self.is_variation_present() {
            "assembled_some_variation"
        } else {
            "just_assembled_reference"
        }
    }

    #[doc(hidden)]
    pub fn kmer_size_for_dump(&self) -> usize {
        if !self.is_variation_present() {
            0
        } else {
            self.minimum_kmer_size().unwrap_or(0)
        }
    }

    /// GATK `--max-mnp-distance` used when regenerating variation events after read augmentation.
    pub fn max_mnp_distance(&self) -> usize {
        self.max_mnp_distance
    }

    /// GATK `AssemblyResultSet.trimTo(trimmedAssemblyRegion)` — haplotypes clipped to padded span.
    pub fn trim_to(&self, trimmed_region: &AssemblyRegion) -> GatkResult<Self> {
        let span = GenomeLoc::new(
            trimmed_region.extended_start.get(),
            trimmed_region.extended_end.get(),
        );
        // GATK `trimDownHaplotypes`: dedupe by (bases, is_reference), not sequence alone.
        let mut trimmed_list: Vec<Haplotype> = Vec::new();
        let mut index_by_hap_key: HashMap<(Vec<u8>, bool), usize> = HashMap::new();
        for h in &self.haplotypes {
            let Some(t) = h.trim(&span, false) else {
                continue;
            };
            // CLONE: needed because owned composite key for dedup/lookup.
            let key = (t.bases.clone(), t.is_reference);
            if let Some(&idx) = index_by_hap_key.get(&key) {
                if h.is_reference {
                    trimmed_list[idx] = t;
                }
            } else {
                index_by_hap_key.insert(key, trimmed_list.len());
                trimmed_list.push(t);
            }
        }
        if trimmed_list.is_empty() {
            let pad = self.padded_reference_start_1based();
            let off = trimmed_region.extended_start.get().saturating_sub(pad) as usize;
            let len = span.reference_span_length() as usize;
            let full_ref = self.reference_bases();
            if off < full_ref.len() && off.saturating_add(len) <= full_ref.len() {
                let mut ref_hap =
                    Haplotype::new(full_ref[off..off.saturating_add(len)].to_vec(), true);
                let mut c = crate::cigar::Cigar::new();
                c.push(len, crate::cigar::CigarOperator::Match);
                ref_hap.cigar = Some(c);
                ref_hap.genome_loc = Some(span);
                trimmed_list.push(ref_hap);
            } else if let Some(ref_h) = self.haplotypes.iter().find(|h| h.is_reference) {
                if let Some(t) = ref_h.trim(&span, false) {
                    trimmed_list.push(t);
                }
            }
        }
        trimmed_list.sort_by(haplotype_size_and_base_order);
        let variation_present =
            trimmed_list.iter().any(|h| !h.is_reference) && trimmed_list.len() > 1;
        let (full_ref, full_pad) = self.event_map_reference();
        let mut variation_events = if variation_present {
            collect_variation_events(
                &trimmed_list,
                full_ref,
                full_pad,
                &self.contig,
                self.max_mnp_distance,
            )
        } else {
            Vec::new()
        };
        for e in &self.variation_events {
            if e.start_1based >= trimmed_region.start && e.start_1based <= trimmed_region.end {
                // CLONE: needed because owned composite key for dedup/lookup.
                let key = (e.start_1based, e.ref_allele.clone(), e.alt_allele.clone());
                if !variation_events.iter().any(|x| {
                    x.start_1based == key.0 && x.ref_allele == key.1 && x.alt_allele == key.2
                }) {
                    // CLONE: needed because owned element into collection.
                    variation_events.push(e.clone());
                }
            }
        }
        crate::event_map::prefer_indel_over_colocated_snps(&mut variation_events);
        variation_events.sort();
        variation_events.dedup();
        let variation_present = variation_present && !variation_events.is_empty();
        Ok(Self {
            haplotypes: trimmed_list,
            variation_present,
            kmer_sizes: self.kmer_sizes.clone(),
            padded_reference_start_1based: self.padded_reference_start_1based,
            // Arc clone — share padded reference bytes with the untrimmed assembly.
            reference_bases: Arc::clone(&self.reference_bases),
            variation_events,
            // CLONE: needed because owned contig id for output record.
            contig: self.contig.clone(),
            max_mnp_distance: self.max_mnp_distance,
        })
    }
}
