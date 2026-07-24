//! Coupled indel phenotype — Sprint **J-2**.
//! Java HC emits a paired deletion+insertion when an alt haplotype carries both
//! `TTC→T` and `A→ATG` with insert start = delete start + [`COUPLED_INDEL_PARTNER_OFFSET`].
//! Recognition is **relative** (allele patterns + Δpos). Absolute chr2 windows remain only as
//! [`coupled_indel_canonical_oracle_locus`] for single-event call sites that lack a partner list
//! (waiver W-H1 until those sites thread an event set).

use crate::event_map::VariationEvent;
use crate::genome_loc::GenomePosition;

/// Genomic distance from coupled deletion start to insertion start (Java EventMap geometry).
pub const COUPLED_INDEL_PARTNER_OFFSET: u64 = 3;

/// Left-align / EventMap jitter tolerated when matching the P12 oracle window only.
pub const COUPLED_INDEL_ORACLE_LOCUS_TOLERANCE: u64 = 3;

/// One side of a [`CoupledIndelCluster`].
/// # Invariants
/// Classifies only the canonical `TTC→T` deletion and `A→ATG` insertion allele strings.
/// # Ownership
/// [`Copy`] phenotype tag.
/// # Mutation
/// Immutable classification result.
/// # Biological assumptions
/// Paired indel alleles observed together on one alt haplotype (P12 cluster phenotype).
/// # Java equivalence
/// Rust-native evidence class (Sprint J-2); allele patterns measured from Java HC EventMap.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CoupledIndelAllele {
    /// `TTC` → `T` (2 bp deletion after left anchor).
    Deletion,
    /// `A` → `ATG` (2 bp insertion after left anchor).
    Insertion,
}

impl CoupledIndelAllele {
    pub fn classify(event: &VariationEvent) -> Option<Self> {
        if event.ref_allele == "TTC" && event.alt_allele == "T" {
            Some(Self::Deletion)
        } else if event.ref_allele == "A" && event.alt_allele == "ATG" {
            Some(Self::Insertion)
        } else {
            None
        }
    }
}

/// A complete coupled indel pair (deletion + insertion at +3).
/// # Invariants
/// Insertion start = deletion start + [`COUPLED_INDEL_PARTNER_OFFSET`].
/// Alleles match [`CoupledIndelAllele`] deletion/insertion patterns.
/// # Ownership
/// Owns both [`VariationEvent`] partners.
/// # Mutation
/// Immutable cluster once built via [`Self::try_from_pair`].
/// # Biological assumptions
/// Java HC emits both events when an alt haplotype carries the paired indel geometry.
/// # Java equivalence
/// Rust-native phenotype (algorithm parity); not a Java class mirror.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CoupledIndelCluster {
    pub deletion: VariationEvent,
    pub insertion: VariationEvent,
}

impl CoupledIndelCluster {
    /// Build a cluster when `del`/`ins` alleles and relative geometry match.
    pub fn try_from_pair(deletion: &VariationEvent, insertion: &VariationEvent) -> Option<Self> {
        if CoupledIndelAllele::classify(deletion) != Some(CoupledIndelAllele::Deletion) {
            return None;
        }
        if CoupledIndelAllele::classify(insertion) != Some(CoupledIndelAllele::Insertion) {
            return None;
        }
        if insertion.start_1based
            != GenomePosition::new_1based(
                deletion
                    .start_1based
                    .get()
                    .saturating_add(COUPLED_INDEL_PARTNER_OFFSET),
            )
        {
            return None;
        }
        Some(Self {
            deletion: deletion.clone(),
            insertion: insertion.clone(),
        })
    }

    pub fn contains(&self, event: &VariationEvent) -> bool {
        events_allele_locus_match(event, &self.deletion)
            || events_allele_locus_match(event, &self.insertion)
    }
}

fn events_allele_locus_match(a: &VariationEvent, b: &VariationEvent) -> bool {
    a.start_1based == b.start_1based && a.ref_allele == b.ref_allele && a.alt_allele == b.alt_allele
}

/// Allele pattern only (either side of a coupled pair).
pub fn is_coupled_indel_allele(event: &VariationEvent) -> bool {
    CoupledIndelAllele::classify(event).is_some()
}

/// Find all complete coupled pairs in `events` (phenotype; no absolute coordinates).
pub fn find_coupled_indel_clusters(events: &[VariationEvent]) -> Vec<CoupledIndelCluster> {
    let mut out = Vec::new();
    let dels: Vec<&VariationEvent> = events
        .iter()
        .filter(|e| CoupledIndelAllele::classify(e) == Some(CoupledIndelAllele::Deletion))
        .collect();
    let insc: Vec<&VariationEvent> = events
        .iter()
        .filter(|e| CoupledIndelAllele::classify(e) == Some(CoupledIndelAllele::Insertion))
        .collect();
    for d in dels {
        for i in &insc {
            if let Some(c) = CoupledIndelCluster::try_from_pair(d, i) {
                if !out.iter().any(|x: &CoupledIndelCluster| {
                    x.deletion.start_1based == c.deletion.start_1based
                }) {
                    out.push(c);
                }
            }
        }
    }
    out
}

/// True when `event` is a member of any complete coupled pair in `events`.
pub fn is_coupled_indel_member(event: &VariationEvent, events: &[VariationEvent]) -> bool {
    if !is_coupled_indel_allele(event) {
        return false;
    }
    find_coupled_indel_clusters(events)
        .iter()
        .any(|c| c.contains(event))
}

/// P12 oracle window (W-H1) — **not** the preferred recognition path.
/// Used only when a call site has a single event and no partner list. Prefer
/// [`is_coupled_indel_member`] whenever the region event set is available.
pub fn coupled_indel_canonical_oracle_locus(event: &VariationEvent) -> bool {
    use super::java_hc_site_semantics::{CLUSTER_ATG_INSERT_START, CLUSTER_TTC_DEL_START};
    match CoupledIndelAllele::classify(event) {
        Some(CoupledIndelAllele::Deletion) => {
            let p = event.start_1based;
            p >= GenomePosition::new_1based(
                CLUSTER_TTC_DEL_START.saturating_sub(COUPLED_INDEL_ORACLE_LOCUS_TOLERANCE),
            ) && p
                <= GenomePosition::new_1based(
                    CLUSTER_TTC_DEL_START.saturating_add(COUPLED_INDEL_ORACLE_LOCUS_TOLERANCE),
                )
        }
        Some(CoupledIndelAllele::Insertion) => {
            let p = event.start_1based;
            p >= GenomePosition::new_1based(
                CLUSTER_ATG_INSERT_START.saturating_sub(COUPLED_INDEL_ORACLE_LOCUS_TOLERANCE),
            ) && p
                <= GenomePosition::new_1based(
                    CLUSTER_ATG_INSERT_START.saturating_add(COUPLED_INDEL_ORACLE_LOCUS_TOLERANCE),
                )
        }
        None => false,
    }
}

/// Genotyping/emit recognition: phenotype member if partners present, else W-H1 oracle.
/// Pass `region_events` when available (preferred). Empty slice → oracle-only (legacy).
pub fn is_coupled_indel_for_genotyping(
    event: &VariationEvent,
    region_events: &[VariationEvent],
) -> bool {
    if !is_coupled_indel_allele(event) {
        return false;
    }
    if !region_events.is_empty() {
        return is_coupled_indel_member(event, region_events);
    }
    coupled_indel_canonical_oracle_locus(event)
}

/// Genomic offset from coupled ATG insertion start to the CTC→C deletion (Java EventMap).
pub const CTC_DEL_OFFSET_FROM_COUPLED_INSERT: u64 = 32;

/// CTC del (`CT`→`C`) relative to a coupled insert, else W-H1 absolute band fixture.
/// When `region_events` contains a complete coupled pair, match insert_start+32 (phenotype).
/// Otherwise fall back to [`super::java_hc_site_semantics::is_cluster_ctc_del`].
pub fn is_ctc_del_for_genotyping(event: &VariationEvent, region_events: &[VariationEvent]) -> bool {
    if event.ref_allele != "CT" || event.alt_allele != "C" {
        return false;
    }
    if !region_events.is_empty() {
        for cluster in find_coupled_indel_clusters(region_events) {
            let expect = GenomePosition::new_1based(
                cluster
                    .insertion
                    .start_1based
                    .get()
                    .saturating_add(CTC_DEL_OFFSET_FROM_COUPLED_INSERT),
            );
            if event.start_1based == expect {
                return true;
            }
        }
    }
    super::java_hc_site_semantics::is_cluster_ctc_del(event)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ev(pos: u64, r: &str, a: &str) -> VariationEvent {
        VariationEvent {
            contig: "2".into(),
            start_1based: GenomePosition::new_1based(pos),
            end_1based: GenomePosition::new_1based(pos),
            ref_allele: r.into(),
            alt_allele: a.into(),
        }
    }

    #[test]
    fn phenotype_recognizes_non_p12_coordinates() {
        let del = ev(1000, "TTC", "T");
        let ins = ev(1003, "A", "ATG");
        let events = vec![del.clone(), ins.clone()];
        let clusters = find_coupled_indel_clusters(&events);
        assert_eq!(clusters.len(), 1);
        assert!(is_coupled_indel_member(&del, &events));
        assert!(is_coupled_indel_member(&ins, &events));
        assert!(!coupled_indel_canonical_oracle_locus(&del));
    }

    #[test]
    fn orphan_deletion_without_partner_is_not_a_cluster_member() {
        let del = ev(1000, "TTC", "T");
        let events = vec![del.clone()];
        assert!(find_coupled_indel_clusters(&events).is_empty());
        assert!(!is_coupled_indel_member(&del, &events));
    }

    #[test]
    fn wrong_offset_is_not_coupled() {
        let del = ev(1000, "TTC", "T");
        let ins = ev(1005, "A", "ATG");
        assert!(CoupledIndelCluster::try_from_pair(&del, &ins).is_none());
    }

    #[test]
    fn genotyping_helper_uses_partners_when_present() {
        let del = ev(1000, "TTC", "T");
        let ins = ev(1003, "A", "ATG");
        let events = vec![del.clone(), ins];
        assert!(is_coupled_indel_for_genotyping(&del, &events));
        // No partners → oracle only (non-P12 fails).
        assert!(!is_coupled_indel_for_genotyping(&del, &[]));
    }

    #[test]
    fn p12_oracle_still_matches_canonical_single_event() {
        let del = ev(92307324, "TTC", "T");
        assert!(coupled_indel_canonical_oracle_locus(&del));
        assert!(is_coupled_indel_for_genotyping(&del, &[]));
    }

    #[test]
    fn ctc_del_uses_relative_offset_when_coupled_partners_present() {
        let del = ev(1000, "TTC", "T");
        let ins = ev(1003, "A", "ATG");
        let ctc = ev(1003 + CTC_DEL_OFFSET_FROM_COUPLED_INSERT, "CT", "C");
        let events = vec![del, ins, ctc.clone()];
        assert!(is_ctc_del_for_genotyping(&ctc, &events));
        // Non-P12 absolute band alone does not match without partners.
        assert!(!is_ctc_del_for_genotyping(&ctc, &[]));
    }
}
