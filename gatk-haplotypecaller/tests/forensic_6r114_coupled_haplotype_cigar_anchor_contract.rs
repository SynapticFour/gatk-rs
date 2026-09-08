//! 6R.114 coordinate-free: coupled deletion/insertion CIGAR anchors on the
//! preceding reference base so EventMap replay reconstructs the indel alleles.
//!
//! GATK 4.4 `EventMap.processCigarForInitialEvents` records I/D against the
//! preceding ref base. A coupled `TTC/T` + `A/ATG` encoding (`TTCA` → `TATG`)
//! therefore needs `Match` through the T of TTC (`lead = ttc_off + 1`), not
//! `Match(ttc_off)` which starts `2D` on the anchor and leftover-Match SNPs `A/G`.
//!
//! Does not pin genomic `92307327` as the contract. Sequence relationship only.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r114_coupled_haplotype_cigar_anchor_contract
//! ```

use gatk_haplotypecaller::cigar::{Cigar, CigarOperator};
use gatk_haplotypecaller::event_map::{
    build_per_haplotype_variation_events, variation_events_for_haplotype, VariationEvent,
};
use gatk_haplotypecaller::genome_loc::GenomeLoc;
use gatk_haplotypecaller::haplotype::Haplotype;
use gatk_haplotypecaller::hc_allele_mapping::create_allele_mapper_with_events;

const REF_MOTIF: &[u8] = b"TTCATGA";
const HAP_MOTIF: &[u8] = b"TATGTGA";
const PREFIX: usize = 79;
const REF_LEN: usize = 191;

fn fill_window(prefix: usize, motif: &[u8], len: usize) -> Vec<u8> {
    let mut v = vec![b'C'; len];
    assert!(prefix + motif.len() <= len);
    v[prefix..prefix + motif.len()].copy_from_slice(motif);
    v
}

fn cigar_anchor(lead: usize, ref_len: usize) -> Cigar {
    let tail = ref_len.saturating_sub(lead.saturating_add(3));
    let mut c = Cigar::new();
    if lead > 0 {
        c.push(lead, CigarOperator::Match);
    }
    c.push(2, CigarOperator::Deletion);
    c.push(1, CigarOperator::Match);
    c.push(2, CigarOperator::Insertion);
    if tail > 0 {
        c.push(tail, CigarOperator::Match);
    }
    c
}

fn hap_with_cigar(bases: Vec<u8>, cigar: Cigar, is_ref: bool, pad: u64) -> Haplotype {
    let len = bases.len() as u64;
    let mut h = Haplotype::new(bases, is_ref);
    h.cigar = Some(cigar);
    h.alignment_start_hap_wrt_ref = 0;
    h.genome_loc = Some(GenomeLoc::new(
        pad,
        pad.saturating_add(len).saturating_sub(1),
    ));
    h
}

fn replay_keys(events: &[VariationEvent], start: u64) -> Vec<String> {
    events
        .iter()
        .filter(|e| e.start_1based.get() == start)
        .map(|e| format!("{}/{}", e.ref_allele, e.alt_allele))
        .collect()
}

fn has_allele(events: &[VariationEvent], start: u64, r: &str, a: &str) -> bool {
    events
        .iter()
        .any(|e| e.start_1based.get() == start && e.ref_allele == r && e.alt_allele == a)
}

/// Case A — proven geometry: prefix 79, window 191, `TTCATGA` vs `TATGTGA`.
#[test]
fn forensic_6r114_anchor_match_replays_coupled_indels_not_snp() {
    let pad: u64 = 1;
    let ref_bases = fill_window(PREFIX, REF_MOTIF, REF_LEN);
    let hap_bases = fill_window(PREFIX, HAP_MOTIF, REF_LEN);
    let mut ref_m = Cigar::new();
    ref_m.push(REF_LEN, CigarOperator::Match);
    let ref_hap = hap_with_cigar(ref_bases.clone(), ref_m, true, pad);

    let old = cigar_anchor(PREFIX, REF_LEN);
    let new = cigar_anchor(PREFIX + 1, REF_LEN);
    assert_eq!(old.to_gatk_string(), "79M2D1M2I109M");
    assert_eq!(new.to_gatk_string(), "80M2D1M2I108M");
    assert_eq!(old.reference_length(), REF_LEN);
    assert_eq!(new.reference_length(), REF_LEN);
    assert_eq!(old.read_length(), REF_LEN);
    assert_eq!(new.read_length(), REF_LEN);

    let ttc_start = pad + PREFIX as u64;
    let atg_start = ttc_start + 3;

    let old_hap = hap_with_cigar(hap_bases.clone(), old, false, pad);
    let old_ev = variation_events_for_haplotype(&old_hap, &ref_hap, &ref_bases, pad, 0, "synth");
    assert_eq!(
        replay_keys(&old_ev, atg_start),
        vec!["A/G".to_string()],
        "lead=ttc_off leftover Match is A/G (6R.113 historical)"
    );
    assert!(
        !has_allele(&old_ev, atg_start, "A", "ATG"),
        "79M2D1M2I must not reconstruct A/ATG"
    );

    let new_hap = hap_with_cigar(hap_bases, new, false, pad);
    let new_ev = variation_events_for_haplotype(&new_hap, &ref_hap, &ref_bases, pad, 0, "synth");
    assert!(
        has_allele(&new_ev, ttc_start, "TTC", "T"),
        "lead=ttc_off+1 must replay TTC/T"
    );
    assert!(
        has_allele(&new_ev, atg_start, "A", "ATG"),
        "lead=ttc_off+1 must replay A/ATG"
    );
    assert!(
        !has_allele(&new_ev, atg_start, "A", "G"),
        "EventMap-compatible CIGAR must not leave A/G at the insertion loc"
    );
}

/// Case B — conservation: motif at window start (`ttc_off = 0`) and away from start.
#[test]
fn forensic_6r114_cigar_consumption_preserved_at_window_edge_and_interior() {
    for prefix in [0usize, 5, 79] {
        let ref_len = prefix.saturating_add(20);
        let lead = prefix + 1;
        let c = cigar_anchor(lead, ref_len);
        assert_eq!(
            c.reference_length(),
            ref_len,
            "prefix={prefix} ref consumption"
        );
        assert_eq!(
            c.read_length(),
            ref_len,
            "prefix={prefix} hap consumption (equal-length coupled alt)"
        );
        let s = c.to_gatk_string();
        assert!(
            s.contains("2D1M2I"),
            "prefix={prefix} keeps 2D1M2I core, got {s}"
        );
        assert!(
            s.starts_with(&format!("{lead}M")) || lead == 0,
            "prefix={prefix} lead Match includes anchor, got {s}"
        );
    }
}

/// Mapper: EventMap-compatible CIGAR puts the insertion hap in the ATG alt pool.
#[test]
fn forensic_6r114_compatible_cigar_fills_insertion_mapper_pool() {
    let pad: u64 = 1;
    let ref_bases = fill_window(PREFIX, REF_MOTIF, REF_LEN);
    let hap_bases = fill_window(PREFIX, HAP_MOTIF, REF_LEN);
    let mut ref_m = Cigar::new();
    ref_m.push(REF_LEN, CigarOperator::Match);
    let ref_hap = hap_with_cigar(ref_bases.clone(), ref_m, true, pad);
    let alt = hap_with_cigar(hap_bases, cigar_anchor(PREFIX + 1, REF_LEN), false, pad);
    let haps = vec![ref_hap, alt];
    let atg_start = pad + PREFIX as u64 + 3;
    let event = VariationEvent::from_alleles("synth", atg_start, "A", "ATG");
    let cache = build_per_haplotype_variation_events(&haps, &ref_bases, pad, 0, "synth");
    let before = cigar_anchor(PREFIX, REF_LEN);
    let old_alt = hap_with_cigar(fill_window(PREFIX, HAP_MOTIF, REF_LEN), before, false, pad);
    let old_haps = vec![
        hap_with_cigar(
            ref_bases.clone(),
            {
                let mut c = Cigar::new();
                c.push(REF_LEN, CigarOperator::Match);
                c
            },
            true,
            pad,
        ),
        old_alt,
    ];
    let old_cache = build_per_haplotype_variation_events(&old_haps, &ref_bases, pad, 0, "synth");
    let old_map = create_allele_mapper_with_events(
        &event,
        atg_start,
        &old_haps,
        pad,
        &ref_bases,
        0,
        false,
        Some(&old_cache),
    );
    assert!(
        old_map.alt_haplotype_indices.is_empty(),
        "79M leftover A/G must not fill the ATG mapper pool"
    );

    let mapping = create_allele_mapper_with_events(
        &event,
        atg_start,
        &haps,
        pad,
        &ref_bases,
        0,
        false,
        Some(&cache),
    );
    assert!(
        !mapping.alt_haplotype_indices.is_empty(),
        "ATG mapper alt pool must be > 0 after EventMap-compatible CIGAR"
    );
}
