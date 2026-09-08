//! 6R.119 coordinate-free: after `trim_to`, SNP-anchor supplements must not
//! materialize a haplotype whose length equals the full padded reference when
//! the live reference haplotype is already the shorter apply window.
//!
//! GATK 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`) has no
//! `ensure_p12_tg_anchor_alt_haplotype`. Other Rust SNP-anchor paths already
//! apply to the trimmed apply window only. Applying a SNP to `full_ref` after
//! trim produced a score-1000 object tagged with apply-window `genome_loc`.
//!
//! The helper still responds only to the cluster T/G event class; this test
//! uses that class as the trigger, not as an emit pin. Sequence relationship:
//! padded ref longer than the trim span; SNP inside the trim span.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --lib forensic_6r119
//! ```

use super::{ensure_p12_tg_anchor_alt_haplotype, P12_CLUSTER_TG_SNP_START};
use crate::alignment::SwParameters;
use crate::assembly_region_iterator::AssemblyRegion;
use crate::assembly_result_set::AssemblyResultSet;
use crate::cigar::{Cigar, CigarOperator};
use crate::feature_context::FeatureContext;
use crate::genome_loc::{GenomeLoc, GenomePosition};
use crate::haplotype::Haplotype;
use crate::read_threading_assembler::AssemblyStatus;
use crate::reference_context::ReferenceContext;

const PAD_BEFORE: u64 = 40;
const FULL_LEN: usize = 80;
const TRIM_BEFORE: u64 = 20;
const TRIM_LEN: usize = 40;

fn dummy_region(start: u64, end: u64) -> AssemblyRegion {
    AssemblyRegion {
        contig: "synth".into(),
        start: GenomePosition::new_1based(start),
        end: GenomePosition::new_1based(end),
        is_active: true,
        extended_start: GenomePosition::new_1based(start),
        extended_end: GenomePosition::new_1based(end),
        extension: 0,
        reads: Vec::new(),
        read_qnames: Vec::new(),
        reference: ReferenceContext::empty(),
        features: FeatureContext::empty(),
        pileup_loci: Vec::new(),
    }
}

fn hap_all_m(bases: Vec<u8>, is_ref: bool, loc: GenomeLoc) -> Haplotype {
    let mut h = Haplotype::new(bases, is_ref);
    let mut c = Cigar::new();
    c.push(h.bases.len(), CigarOperator::Match);
    h.cigar = Some(c);
    h.genome_loc = Some(loc);
    h.alignment_start_hap_wrt_ref = 0;
    h
}

/// After trim, the TG-anchor helper must not push a full-padded-ref SNP haplotype.
#[test]
fn forensic_6r119_tg_anchor_does_not_materialize_full_padded_ref_snp() {
    let snp = P12_CLUSTER_TG_SNP_START;
    let pad = snp.saturating_sub(PAD_BEFORE);
    let full_end = pad + FULL_LEN as u64 - 1;
    let trim_start = snp.saturating_sub(TRIM_BEFORE);
    let trim_end = trim_start + TRIM_LEN as u64 - 1;
    assert!(trim_start >= pad && trim_end <= full_end);
    assert!(snp >= trim_start && snp <= trim_end);

    let snp_off = (snp - pad) as usize;
    let mut full_ref = vec![b'C'; FULL_LEN];
    full_ref[snp_off] = b'T';
    let mut full_alt = full_ref.clone();
    full_alt[snp_off] = b'G';

    let loc = GenomeLoc::new(pad, full_end);
    let ref_hap = hap_all_m(full_ref.clone(), true, loc);
    let alt_hap = hap_all_m(full_alt, false, loc);
    let assembly = AssemblyResultSet::from_assembly_for_calling_owned(
        AssemblyStatus::AssembledSomeVariation,
        25,
        vec![ref_hap, alt_hap],
        full_ref,
        pad,
        "synth",
        0,
    );
    assert!(
        assembly
            .variation_events()
            .iter()
            .any(|e| e.start_1based.get() == snp && e.ref_allele == "T" && e.alt_allele == "G"),
        "untrimmed EventMap must carry the T/G event class: {:?}",
        assembly
            .variation_events()
            .iter()
            .map(|e| format!("{}/{}", e.ref_allele, e.alt_allele))
            .collect::<Vec<_>>()
    );

    let trimmed = assembly
        .trim_to(&dummy_region(trim_start, trim_end))
        .expect("trim");
    let apply_len = trimmed
        .haplotypes
        .iter()
        .find(|h| h.is_reference)
        .map(|h| h.bases.len())
        .expect("trimmed ref");
    assert_eq!(apply_len, TRIM_LEN);
    assert_eq!(trimmed.reference_bases().len(), FULL_LEN);
    assert!(trimmed.haplotypes.iter().all(|h| h.bases.len() == TRIM_LEN));

    let mut live = trimmed;
    let sw = SwParameters::gatk_haplotype_to_reference();
    ensure_p12_tg_anchor_alt_haplotype(&mut live, &sw).expect("tg anchor");

    let n_full_len = live
        .haplotypes
        .iter()
        .filter(|h| !h.is_reference && h.bases.len() == FULL_LEN)
        .count();
    let n_full_len_supplement = live
        .haplotypes
        .iter()
        .filter(|h| !h.is_reference && h.bases.len() == FULL_LEN && (h.score - 1000.0).abs() < 1e-6)
        .count();
    assert_eq!(
        n_full_len, 0,
        "post-trim TG-anchor must not re-insert a full-padded-ref SNP haplotype"
    );
    assert_eq!(n_full_len_supplement, 0);
    assert!(
        live.haplotypes.iter().all(|h| h.bases.len() == apply_len),
        "live haplotypes must stay on the apply-window length"
    );
}
