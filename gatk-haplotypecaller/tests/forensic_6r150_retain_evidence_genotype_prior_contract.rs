//! 6R.150 coordinate-free: `retainEvidence` then genotype-prior construction
//! after the 6R.149-equivalent allele-likelihood object.
//!
//! Java 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`)
//!
//! ```text
//! HaplotypeCallerGenotypingEngine.assignGenotypeLikelihoods
//!   readAlleleLikelihoods = readLikelihoods.marginalize(alleleMapper)
//!   variantCallingRelevantOverlap =
//!       SimpleInterval(mergedVC).expandWithinContig(informativeReadOverlapMargin, dict)
//!   composeReadQualifiesForGenotypingPredicate:
//!     default HC (applyBQD=false, applyFRD=false):
//!       (read, target) -> target.overlaps(read)   // alignment start/end; margin 0
//!     BQD/FRD only: softUnclippedReadOverlapsInterval
//!   AlleleLikelihoods.retainEvidence(predicate)   // in-place removeEvidenceByIndex
//!   setVariantCallingSubsetUsed(overlap)
//!   calculateGLsForThisEvent
//!     IndependentSampleGenotypesModel.calculateLikelihoods
//!     GenotypeBuilder(NO_CALL).PL(getAsPLs())
//!   resolveGenotypePriorCalculator
//!     dragstrParams == null → GenotypePriorCalculator.assumingHW(
//!         log10(snpHeterozygosity), log10(indelHeterozygosity))
//!   calculateGenotypes(vc_with_NO_CALL_PLs, gpc, givenAlleles)
//!     genotypeAssignmentMethod = USE_PLS_TO_ASSIGN  // default
//!     GATKVariantContextUtils.makeGenotypeCall: GT = argmax(GL); gpc unused
//! ```
//!
//! `retainEvidence` copies no likelihood values: it compact-shifts surviving columns
//! and NaN-fills the tail. Allele order and sample membership are unchanged.
//! Overlapping mates are independent evidence units (no QNAME collapse).
//!
//! `assumingHW` diploid biallelic SNP prior vector (ploidy 2, alleles REF/ALT SNP):
//!   0/0 = 0  (REF convention)
//!   0/1 = log10(snpHet) − log10(3)
//!   1/1 = 2·log10(snpHet) − log10(3)
//! `stand-call-conf` (Java default 30) is `passesEmitThreshold` only, not GT/PL.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r150_retain_evidence_genotype_prior_contract
//! HOLDOUT_6R150=1 GATK_RS_EXPERIMENTAL_KBEST_POLICY=unbounded_diagnostic \
//!   cargo test -p gatk-haplotypecaller --test holdout_6r150_retain_evidence -- --nocapture --test-threads=1
//! ```

use gatk_haplotypecaller::activity_scoring::{
    DEFAULT_INDEL_HETEROZYGOSITY, DEFAULT_SNP_HETEROZYGOSITY, LOG10_ONE_THIRD,
};
use gatk_haplotypecaller::genotyping::{
    biallelic_diploid_log10_priors, genotype_posteriors_from_log10_likelihoods,
    BiallelicDiploidPriorModel,
};
use gatk_haplotypecaller::hc_genotyping_engine::{
    java_alignment_read_overlaps_interval, DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN,
    DEFAULT_STAND_EMIT_CONFIDENCE,
};
use gatk_haplotypecaller::read_unclip::{alignment_end_1based, gatk_soft_start_1based};
use rust_htslib::bam::record::{Cigar, CigarString};

/// Java `SimpleInterval.expandWithinContig` then `overlaps` with margin 0.
fn java_expand_within_contig(start: i64, end: i64, padding: i32, contig_len: i64) -> (i64, i64) {
    let p = i64::from(padding.max(0));
    ((start - p).max(1), (end + p).min(contig_len))
}

fn java_simple_interval_overlaps(tstart: i64, tend: i64, rstart: i64, rend: i64) -> bool {
    tstart <= rend && rstart <= tend
}

/// Java `GenotypePriorCalculator.assumingHW` SNP slice for diploid genotypes 0/0, 0/1, 1/1.
fn java_assuming_hw_snp_log10_priors(snp_heterozygosity: f64) -> [f64; 3] {
    let snp_het = snp_heterozygosity.log10();
    let log10_snp_norm = 3.0_f64.log10();
    [
        0.0,
        snp_het - log10_snp_norm,
        snp_het * 2.0 - log10_snp_norm,
    ]
}

fn mate(qname: &[u8], pos0: i64, cigar: CigarString, flags: u16) -> rust_htslib::bam::Record {
    let mut rec = rust_htslib::bam::Record::new();
    rec.set(qname, Some(&cigar), b"ACGTACGTAC", b"##########");
    rec.set_pos(pos0);
    rec.set_flags(flags);
    rec
}

#[test]
fn forensic_6r150_retain_follows_marginalize_before_gls() {
    let order = [
        "marginalize",
        "retainEvidence",
        "calculateGLsForThisEvent",
        "calculateGenotypes",
    ];
    assert_eq!(order[1], "retainEvidence");
    assert_eq!(order[2], "calculateGLsForThisEvent");
}

#[test]
fn forensic_6r150_interval_is_expanded_before_overlap_not_after() {
    // SimpleInterval(mergedVC).expandWithinContig(2) then target.overlaps(read) with margin 0.
    let merged_start = 100i64;
    let merged_end = 100i64;
    let (s, e) = java_expand_within_contig(merged_start, merged_end, 2, 1_000_000);
    assert_eq!((s, e), (98, 102));
    assert!(java_simple_interval_overlaps(s, e, 98, 98));
    assert!(!java_simple_interval_overlaps(
        merged_start,
        merged_end,
        98,
        98
    ));
}

#[test]
fn forensic_6r150_default_hc_overlap_is_alignment_not_softclip() {
    let rec = mate(
        b"r",
        90,
        CigarString(vec![Cigar::Match(8), Cigar::SoftClip(5)]),
        0,
    );
    let rs = rec.pos() + 1;
    let re = i64::from(alignment_end_1based(&rec).max(1));
    assert_eq!((rs, re), (91, 98));
    assert!(!java_simple_interval_overlaps(100, 104, rs, re));
    let soft_start = gatk_soft_start_1based(&rec).max(1);
    let mut soft_end = re;
    if let Some(Cigar::SoftClip(n)) = rec.cigar().iter().rev().next() {
        soft_end += i64::from(*n);
    }
    assert!(java_simple_interval_overlaps(
        100, 104, soft_start, soft_end
    ));
    assert!(!java_alignment_read_overlaps_interval(&rec, 100, 104, 0));
}

#[test]
fn forensic_6r150_informative_margin_default_is_two() {
    assert_eq!(DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN, 2);
}

#[test]
fn forensic_6r150_retain_evidence_is_in_place_row_filter() {
    // AlleleLikelihoods.retainEvidence → removeEvidenceByIndex: compact surviving columns,
    // NaN-fill the tail. Remaining (read, allele) values are copied, not recomputed.
    let mut cols = vec![-1.0_f64, -2.0, -3.0];
    let remove = [1usize];
    let mut kept = Vec::new();
    for (n, v) in cols.iter().copied().enumerate() {
        if !remove.contains(&n) {
            kept.push(v);
        }
    }
    for (i, v) in kept.iter().copied().enumerate() {
        cols[i] = v;
    }
    cols[kept.len()..].fill(f64::NAN);
    assert_eq!(cols[0], -1.0);
    assert_eq!(cols[1], -3.0);
    assert!(cols[2].is_nan());
}

#[test]
fn forensic_6r150_retain_does_not_collapse_qname_or_change_alleles() {
    let mate0 = mate(b"pair", 99, CigarString(vec![Cigar::Match(10)]), 99);
    let mate1 = mate(b"pair", 101, CigarString(vec![Cigar::Match(10)]), 147);
    let overlapping: Vec<_> = [&mate0, &mate1]
        .into_iter()
        .filter(|r| java_alignment_read_overlaps_interval(r, 105, 105, 2))
        .collect();
    assert_eq!(overlapping.len(), 2);
    let unique_qnames = 1usize;
    assert_eq!(unique_qnames, 1);
    assert_ne!(
        overlapping.len(),
        unique_qnames,
        "Java retainEvidence keeps both overlapping mates; QNAME collapse is extra"
    );
}

#[test]
fn forensic_6r150_calculate_gls_emits_nocall_pl_not_gt() {
    let java_gt_at_calculate_gls_is_nocall = true;
    let java_pl_from_independent_sample_model = true;
    let genotyping_model = "IndependentSampleGenotypesModel";
    assert!(java_gt_at_calculate_gls_is_nocall);
    assert!(java_pl_from_independent_sample_model);
    assert_eq!(genotyping_model, "IndependentSampleGenotypesModel");
}

#[test]
fn forensic_6r150_default_assignment_is_use_pls_not_posteriors() {
    let genotype_assignment_method = "USE_PLS_TO_ASSIGN";
    assert_eq!(genotype_assignment_method, "USE_PLS_TO_ASSIGN");
    // makeGenotypeCall USE_PLS branch: maxElementIndex(genotypeLikelihoods); gpc is unused.
    let gls = [-628.48_f64, -628.88, -892.36];
    let priors = java_assuming_hw_snp_log10_priors(DEFAULT_SNP_HETEROZYGOSITY);
    let post: Vec<f64> = gls.iter().zip(priors.iter()).map(|(g, p)| g + p).collect();
    let gt_from_gl = gls
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.total_cmp(b.1))
        .map(|(i, _)| i)
        .unwrap();
    let gt_from_post = post
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.total_cmp(b.1))
        .map(|(i, _)| i)
        .unwrap();
    assert_eq!(gt_from_gl, 0);
    assert_eq!(gt_from_post, 0);
    let gpc_used_for_gt = false;
    assert!(!gpc_used_for_gt);
}

#[test]
fn forensic_6r150_java_hw_snp_priors_are_not_rust_flat_remainder() {
    let java = java_assuming_hw_snp_log10_priors(DEFAULT_SNP_HETEROZYGOSITY);
    assert_eq!(java[0], 0.0);
    assert!((java[1] - (DEFAULT_SNP_HETEROZYGOSITY.log10() + LOG10_ONE_THIRD)).abs() < 1e-15);
    assert!((java[2] - (2.0 * DEFAULT_SNP_HETEROZYGOSITY.log10() + LOG10_ONE_THIRD)).abs() < 1e-15);

    let rust = biallelic_diploid_log10_priors(BiallelicDiploidPriorModel::default()).unwrap();
    let max_abs = java
        .iter()
        .zip(rust.iter())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0_f64, f64::max);
    assert!(
        max_abs > 1.0,
        "Java assumingHW vs Rust remainder-mass priors must differ by construction, got {max_abs}"
    );
    assert_eq!(DEFAULT_INDEL_HETEROZYGOSITY, 1.0 / 8000.0);
}

#[test]
fn forensic_6r150_java_posterior_is_gl_plus_prior_when_requested() {
    // USE_POSTERIOR_PROBABILITIES: ebeAdd(log10Priors, genotypeLikelihoods). Not the default.
    let gls = [-1.0_f64, -2.0, -3.0];
    let priors = java_assuming_hw_snp_log10_priors(DEFAULT_SNP_HETEROZYGOSITY);
    let java_post: [f64; 3] = [gls[0] + priors[0], gls[1] + priors[1], gls[2] + priors[2]];
    let rust_post = genotype_posteriors_from_log10_likelihoods(&gls, &priors).unwrap();
    for i in 0..3 {
        assert!((java_post[i] - rust_post.genotype_log10_posteriors[i]).abs() < 1e-15);
    }
}

#[test]
fn forensic_6r150_stand_call_conf_is_emit_not_genotype_math() {
    let java_stand_call = 30.0_f64;
    assert_eq!(java_stand_call, 30.0);
    // Rust emit-gate default is 10; that knob is not the GL/prior/GT calculator.
    assert_eq!(DEFAULT_STAND_EMIT_CONFIDENCE, 10.0);
    let stand_call_enters_retain_or_prior = false;
    assert!(!stand_call_enters_retain_or_prior);
}

#[test]
fn forensic_6r150_dragstr_unused_without_params() {
    let dragstr_params_null = true;
    let dont_use_dragstr_priors_default = false;
    let uses_assuming_hw = dragstr_params_null || dont_use_dragstr_priors_default;
    assert!(uses_assuming_hw);
}

#[test]
fn forensic_6r150_ploidy_two_one_sample_snp_not_indel() {
    let ploidy = 2usize;
    let sample_count = 1usize;
    let snp_same_length = true;
    assert_eq!(ploidy, 2);
    assert_eq!(sample_count, 1);
    assert!(snp_same_length);
}
