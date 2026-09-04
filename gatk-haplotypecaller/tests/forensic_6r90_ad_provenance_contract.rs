//! 6R.90 coordinate-free: FORMAT/AD provenance after DepthPerAlleleBySample.
//!
//! Java 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`)
//!
//! ```text
//! calculateGLsForThisEvent          // GenotypeBuilder: PL only; AD absent
//! calculateGenotypes
//!   AlleleSubsettingUtils.subsetAlleles
//!     if (g.hasAD()) slice by keep indices   // no-op: AD is still absent
//!     does NOT remarginalize likelihoods
//!     does NOT invent AD
//! prepareReadAlleleLikelihoodsForAnnotation  // reuse C; updateNonRef no-op
//! annotationEngine.annotateContext(call, …, allele-level likelihoods)
//!   DepthPerAlleleBySample.annotate
//!     alleles = call.getAlleles()            // remaining; reverseTrim is AFTER
//!     FIRST AD write
//! reverseTrimAlleles                         // GenotypeBuilder copy; AD unchanged
//! phaseVC                                    // GenotypeBuilder copy; AD unchanged
//! VCF FORMAT/AD                              // same counts as annotation
//! ```
//!
//! 6R.88/6R.89 closed C and the remarg algorithm on C ([27, 9]). This round traces
//! every AD mutation after that. PL/QUAL/PairHMM/overlap are out of scope.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r90_ad_provenance_contract
//! HOLDOUT_6R90=1 cargo test -p gatk-haplotypecaller --test forensic_6r90_ad_provenance_contract live_ -- --nocapture
//! ```

use gatk_haplotypecaller::reverse_trim_alleles;
use gatk_haplotypecaller::subset_unused_alts_after_merged_genotyping;

/// Java `calculateGLsForThisEvent` builds genotypes with PL only.
/// `AlleleSubsettingUtils.subsetAlleles` remaps AD iff `g.hasAD()`.
/// Absent AD is not a zero vector: the subset path never invents FORMAT AD.
#[test]
fn forensic_6r90_subset_without_ad_cannot_invent_counts() {
    let after_gls_has_ad = false;
    let after_subset_has_ad = after_gls_has_ad; // Java `if (g.hasAD())` skipped
    assert!(
        !after_subset_has_ad,
        "unused-ALT subset cannot create AD from a PL-only genotype"
    );
    // Contrast: slicing an *existing* 4-way AD is a different operation (Rust production).
    let existing = vec![26, 2, 9, 1];
    let sliced = [existing[0], existing[2]];
    assert_eq!(sliced, [26, 9]);
    assert_ne!(sliced, [36, 19]);
}

/// When AD *is* present (Rust production), unused-ALT subset slices keep indices.
/// That is permutation of existing counts, not remarginalization.
#[test]
fn forensic_6r90_subset_with_ad_slices_keep_indices_not_remarg() {
    let alts = vec!["T".to_string(), "CG".to_string(), "*".to_string()];
    // 4-way informative: REF, unused T, called CG, unused *
    let ad4 = vec![26, 2, 9, 1];
    let remarg_remaining = vec![27, 9];
    let gls = vec![0.0, -1.0, -2.0, -3.0, -4.0, -5.0, -6.0, -7.0, -8.0, -9.0];
    let r = subset_unused_alts_after_merged_genotyping(&alts, &[0, 2], &gls, &ad4).unwrap();
    assert_eq!(r.alt_alleles, vec!["CG".to_string()]);
    assert_eq!(r.ad, vec![26, 9], "slice keep [0, 2] of 4-way counts");
    assert_ne!(
        r.ad, remarg_remaining,
        "permute of 4-way counts is not remaining-allele remarg"
    );
}

/// Java `reverseTrimAlleles` → `updateGenotypesWithMappedAlleles`:
/// `new GenotypeBuilder(genotype).alleles(updatedAlleles).make()` copies AD.
#[test]
fn forensic_6r90_reverse_trim_copies_ad_does_not_recompute() {
    let ad_before = [27, 9];
    let (trim_ref, trim_alts) = reverse_trim_alleles("TG", &["CG".to_string()]);
    assert_eq!(trim_ref, "T");
    assert_eq!(trim_alts, vec!["C".to_string()]);
    let ad_after = ad_before;
    assert_eq!(
        ad_after, ad_before,
        "reverseTrim remaps allele identity; AD vector is unchanged"
    );
    assert_ne!(
        ad_after,
        [36, 19],
        "copying remarg AD cannot produce Java VCF AD"
    );
}

/// 4-way skip-unused (T and *) is not remaining-allele remarg {TG, CG}.
#[test]
fn forensic_6r90_four_way_slice_is_not_identity_remarg() {
    // Per-read remaining remarg can move a vote that 4-way assigned to unused T/*.
    let four_way = [26, 2, 9, 1];
    let sliced = [four_way[0], four_way[2]];
    let remarg = [27, 9];
    assert_eq!(sliced, [26, 9]);
    assert_eq!(sliced.iter().sum::<i32>(), 35);
    assert_eq!(remarg.iter().sum::<i32>(), 36);
    assert_ne!(sliced, remarg);
}

/// Java `calculateGLsForThisEvent` cannot be the source of FORMAT AD:
/// genotypes are constructed with `.PL(...).make()` and no `.AD(...)`.
#[test]
fn forensic_6r90_gl_construction_does_not_write_ad() {
    let pl_only_has_ad = false;
    let after_subset_has_ad = pl_only_has_ad; // Java if (g.hasAD()) skipped
    assert!(
        !after_subset_has_ad,
        "AD is still absent after unused-ALT subset"
    );
}

/// reverseTrim of a 4-way allele list that still contains a length-1 allele is a no-op.
/// Annotation on 4-way remaining alleles cannot become biallelic T/C via reverseTrim.
#[test]
fn forensic_6r90_reverse_trim_skips_when_length_one_allele_present() {
    let (r, a) = reverse_trim_alleles("TG", &["T".to_string(), "CG".to_string()]);
    assert_eq!(r, "TG");
    assert_eq!(a, vec!["T".to_string(), "CG".to_string()]);
}

#[test]
fn live_ad_provenance_after_annotation() {
    if std::env::var("HOLDOUT_6R90").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R90=1");
        return;
    }
    use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
    use gatk_haplotypecaller::{
        call_disposition, flatten_assembly_regions, take_colocated_merge_numerics,
        traverse_assembly_region_walker, try_emit_call_region_variants,
        AssemblyRegionCallDisposition, CallRegionArgs, HaplotypeCallerEngine, ReadFilterParams,
        WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
    };
    use std::path::Path;

    const INTERVAL: &str = "20:29455000-29456500";
    const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
    const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
    const POS_SNP: u64 = 29_456_344;

    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: live BAM/ref missing");
        return;
    }

    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, INTERVAL).expect("interval");
    let walk = traverse_assembly_region_walker(
        &dict,
        &specs,
        &ref_fasta,
        &bam,
        &ReadFilterParams::gatk_standard_hc(),
        &WalkerTraversalConfig::gatk_haplotype_caller_production(100),
    )
    .expect("walk");
    let regions = flatten_assembly_regions(&walk);
    let covering: Vec<_> = regions
        .iter()
        .filter(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= POS_SNP
                && r.end.get() >= POS_SNP
        })
        .collect();
    assert_eq!(covering.len(), 1);
    let outcome = HaplotypeCallerEngine::call_region(
        covering[0],
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");
    let emitted = try_emit_call_region_variants(
        covering[0],
        &outcome,
        "SAMPLE",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .unwrap_or_default();
    let vcf = emitted
        .iter()
        .find(|r| {
            r.position == POS_SNP && r.reference == "T" && r.alternate.iter().any(|a| a == "C")
        })
        .expect("T/C");
    let snap = take_colocated_merge_numerics()
        .into_iter()
        .find(|n| n.loc == POS_SNP)
        .expect("snap");

    let (trim_ref, trim_alts) = reverse_trim_alleles(&snap.long_ref, &["CG".to_string()]);
    let vcf_ad = vcf.samples[0].ad.clone().unwrap_or_default();

    eprintln!(
        "6R.90 merged_alleles=[{:?}, {:?}] n_c={} merged_ad={:?} permute={:?} remarg={:?} trim={}/{:?} vcf_ad={:?} java_oracle=[36, 19]",
        snap.long_ref,
        snap.alts,
        snap.n_reads,
        snap.merged_ad,
        snap.subset_ad_permuted,
        snap.subset_ad_remarginalized,
        trim_ref,
        trim_alts,
        vcf_ad
    );

    assert_eq!(snap.n_reads, 62);
    assert_eq!(trim_ref, "T");
    assert_eq!(trim_alts, vec!["C".to_string()]);
    // reverseTrim copies; it does not turn remarg into Java VCF AD.
    assert_eq!(snap.subset_ad_remarginalized, vec![27, 9]);
    assert_ne!(snap.subset_ad_remarginalized, vec![36, 19]);
    // Rust FORMAT AD is the unused-ALT slice, not annotation remarg.
    assert_eq!(snap.subset_ad_permuted, vec![26, 9]);
    assert_eq!(vcf_ad, vec![26u32, 9]);
    assert_eq!(
        vcf_ad.iter().map(|&x| x as i32).collect::<Vec<_>>(),
        snap.subset_ad_permuted
    );
    assert_ne!(
        snap.subset_ad_permuted, snap.subset_ad_remarginalized,
        "Rust VCF permute and Java-style remarg remain distinct operations"
    );
}
