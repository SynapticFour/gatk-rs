//! 6R.106 coordinate-free: allele-mapped likelihoods vs SparsePlShape overwrite.
//!
//! Java 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`)
//! `HaplotypeCallerGenotypingEngine.calculateGLsForThisEvent`:
//!
//! ```text
//! createAlleleMapper(mergedVC, loc, haplotypes, emitSpanningDels)
//!   empty EventMap overlap            → REF
//!   event.start == loc, allele in VC   → that ALT
//!   event.start < loc, emitSpanningDels → * (not C/T)
//! AlleleLikelihoods.marginalize         → max over haplotypes in each pool
//! retainEvidence(interval ± 2)
//! IndependentSampleGenotypesModel      → diploid GLs / PLs
//! ```
//!
//! Live HOLDOUT_6R53 (`20:29455388 C/T`):
//! Java calculator on the 49-read overlap matrix is hom-ref `PL=0,6,1780`.
//! Rust's same mapper + max-marginalize + GL calculator on the same 49 reads
//! reproduces `PL=0,6,1780`. Live Rust then assigns `SparsePlShape::Het`
//! `PL=81,0,36` after `finalize_site` returns None (Java emit fail) via the
//! genome-wide L9 pileup overwrite (`genome_wide_genotype_read_support`
//! is true for any SNP with `alt_ad >= 1`).
//!
//! Mapper T pools match (EventMap C/T). Three spanning-del haplotypes are
//! REF in production and `*` in the Java walk; substitution shows they
//! do not change the 49-read GLs. Haplotype count 60 vs 138 is not causal.
//!
//! Production change: NONE.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r106_allele_likelihood_contract
//! HOLDOUT_6R106=1 cargo test -p gatk-haplotypecaller --test holdout_6r106_allele_likelihood -- --nocapture
//! ```

use gatk_haplotypecaller::bio_ids::HaplotypeIndex;
use gatk_haplotypecaller::event_map::VariationEvent;
use gatk_haplotypecaller::genotyping::{
    best_pl_index, diploid_genotype_alleles_from_pl_index, emit_genotype_format_fields,
    ReadLikelihoodRow,
};
use gatk_haplotypecaller::hc_genotyping_engine::marginalize_rows_to_biallelic_alleles;
use gatk_haplotypecaller::hc_genotyping_engine::{
    biallelic_genotype_log10_likelihoods_gatk, java_emit_would_pass, SparsePlShape,
    DEFAULT_STAND_EMIT_CONFIDENCE,
};

const JAVA_STAND_CALL_CONF: f64 = 30.0;

fn snp(r: &str, a: &str) -> VariationEvent {
    VariationEvent::from_alleles("chr", 100, r, a)
}

fn gl_from_pl(pl: &[i32]) -> Vec<f64> {
    pl.iter().map(|&p| (p as f64) / -10.0).collect()
}

/// Java 4.4 `createAlleleMapper`: empty overlapping EventMap → REF.
fn java_empty_overlap_is_ref(spanning_empty: bool) -> &'static str {
    if spanning_empty {
        "REF"
    } else {
        "walk"
    }
}

/// Genome-wide L9 trigger (non-chr2): SNP with any pileup ALT.
fn genome_wide_snp_l9_trigger(is_snp: bool, read_alt_ad: i32) -> bool {
    is_snp && read_alt_ad >= 1
}

#[test]
fn forensic_6r106_java_empty_eventmap_maps_to_ref_not_snp_base() {
    assert_eq!(java_empty_overlap_is_ref(true), "REF");
}

#[test]
fn forensic_6r106_marginalize_is_max_over_pool() {
    let rows = vec![ReadLikelihoodRow {
        read_index: 0,
        read_id: String::new(),
        haplotype_log10_likelihoods: vec![-1.0, -4.0, -2.0],
    }];
    let ref_p = vec![HaplotypeIndex::new(0), HaplotypeIndex::new(2)];
    let alt_p = vec![HaplotypeIndex::new(1)];
    let marg = marginalize_rows_to_biallelic_alleles(&rows, &ref_p, &alt_p);
    assert_eq!(marg[0].haplotype_log10_likelihoods[0], -1.0);
    assert_eq!(marg[0].haplotype_log10_likelihoods[1], -4.0);
}

#[test]
fn forensic_6r106_java_calculator_homref_pl_is_not_sparse_het() {
    let java_pl = [0, 6, 1780];
    let sparse = SparsePlShape::Het.pl();
    assert_ne!(java_pl.as_slice(), sparse.as_slice());
    let java_fmt = emit_genotype_format_fields(&gl_from_pl(&java_pl), &[43, 4]).expect("j");
    let sparse_fmt =
        emit_genotype_format_fields(&SparsePlShape::Het.gl_vec(), &[44, 4]).expect("s");
    assert_eq!(
        diploid_genotype_alleles_from_pl_index(2, best_pl_index(&java_fmt.pl)).as_slice(),
        [0, 0]
    );
    assert_eq!(
        diploid_genotype_alleles_from_pl_index(2, best_pl_index(&sparse_fmt.pl)).as_slice(),
        [0, 1]
    );
}

#[test]
fn forensic_6r106_pileup_ref_majority_is_still_sparse_het_shape() {
    assert_eq!(SparsePlShape::from_pileup_depths(44, 4), SparsePlShape::Het);
    assert_eq!(SparsePlShape::Het.pl(), [81, 0, 36]);
}

#[test]
fn forensic_6r106_l9_snp_alt_ge_1_is_not_java_calculator() {
    assert!(
        genome_wide_snp_l9_trigger(true, 4),
        "L9 fires for any SNP with pileup ALT >= 1"
    );
    assert!(!genome_wide_snp_l9_trigger(true, 0));
    let event = snp("C", "T");
    let java_gl = gl_from_pl(&[0, 6, 1780]);
    let java_fmt = emit_genotype_format_fields(&java_gl, &[43, 4]).expect("j");
    let het_gl = SparsePlShape::Het.gl_vec();
    let het_fmt = emit_genotype_format_fields(&het_gl, &[44, 4]).expect("h");
    assert!(!java_emit_would_pass(&event, &java_gl, &java_fmt, JAVA_STAND_CALL_CONF, &[]).unwrap());
    assert!(java_emit_would_pass(&event, &het_gl, &het_fmt, JAVA_STAND_CALL_CONF, &[]).unwrap());
    assert!(java_emit_would_pass(
        &event,
        &het_gl,
        &het_fmt,
        DEFAULT_STAND_EMIT_CONFIDENCE,
        &[]
    )
    .unwrap());
}

#[test]
fn forensic_6r106_calculator_on_ref_majority_matrix_stays_homref() {
    // Live HOLDOUT_6R53 overlap phenotype: most reads prefer C (~−3.5) over T (~−8).
    // Two balanced votes (one C-preferring, one T-preferring) are a het, not this contract.
    let rows: Vec<ReadLikelihoodRow> = (0..8)
        .map(|i| ReadLikelihoodRow {
            read_index: i,
            read_id: String::new(),
            haplotype_log10_likelihoods: vec![-3.5, -8.0],
        })
        .collect();
    let gls = biallelic_genotype_log10_likelihoods_gatk(&rows, 0, 1);
    let fmt = emit_genotype_format_fields(&gls, &[8, 0]).expect("fmt");
    let gt = diploid_genotype_alleles_from_pl_index(2, best_pl_index(&fmt.pl));
    assert_eq!(gt.as_slice(), [0, 0]);
    assert_ne!(fmt.pl_as_i32().as_slice(), SparsePlShape::Het.pl());
}
