//! 6R.310 holdout: a two-base repeat, not the canonical A>AT event.
//! Skipped unless `HOLDOUT_6R310=1`. Does not run HOLDOUT_6R243.
//!
//! ```text
//! HOLDOUT_6R310=1 cargo test -p gatk-haplotypecaller --test holdout_6r310_str_production_patch -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::SequenceDictionary;
use gatk_haplotypecaller::{
    tandem_repeat_at_event, AssemblyRegion, AssemblyRegionTrimmer, AssemblyRegionTrimmerConfig,
    GenomePosition, ReferenceContext, TrimVariant,
};

#[test]
fn holdout_6r310_str_production_patch() {
    if std::env::var("HOLDOUT_6R310").ok().as_deref() != Some("1") {
        eprintln!("6R310 holdout skipped; set HOLDOUT_6R310=1");
        return;
    }
    // Anchor G, unit AT. Four reference copies, five once the inserted AT is counted.
    let bases = b"GATATATAT";
    let rep = tandem_repeat_at_event(50, bases, 50, b"G", b"GAT").expect("repeat");
    assert_eq!(rep.unit, b"AT");
    assert_eq!(rep.counts, vec![4, 5]);
    assert_eq!(rep.padding_bases(), 10);
    assert_ne!(
        (rep.unit.as_slice(), rep.counts.as_slice()),
        (b"T".as_slice(), [8, 9].as_slice())
    );

    let mut dict = SequenceDictionary::new();
    dict.add_contig("chr1".into(), 500);
    let trimmer =
        AssemblyRegionTrimmer::new(AssemblyRegionTrimmerConfig::gatk_defaults(), &dict, "chr1");
    let region = AssemblyRegion {
        contig: "chr1".into(),
        start: GenomePosition::new_1based(40),
        end: GenomePosition::new_1based(80),
        is_active: true,
        extended_start: GenomePosition::new_1based(1),
        extended_end: GenomePosition::new_1based(500),
        extension: 100,
        reads: Vec::new(),
        read_qnames: Vec::new(),
        reference: ReferenceContext::empty(),
        features: gatk_haplotypecaller::FeatureContext::empty(),
        pileup_loci: Vec::new(),
    };
    let vars = vec![TrimVariant {
        contig: "chr1".into(),
        start: 50,
        end: 50,
        is_indel: true,
        ref_allele: "G".into(),
        alt_allele: "GAT".into(),
    }];
    let ctx = ReferenceContext {
        contig: "chr1".into(),
        start: 50,
        end: 58,
        window_start: 50,
        window_end: 58,
        bases: gatk_haplotypecaller::reference_context::SharedBases::from_slice(bases),
    };
    let res = trimmer.trim(&region, &vars, Some(&ctx));
    assert_eq!(res.padded_variant_start, Some(1));
    assert_eq!(res.padded_variant_end, Some(50 + 75 + 10));
    println!("6R310\tclassification\tSTR_PRODUCTION_PATCH_REPRODUCES_JAVA_GEOMETRY");
}
