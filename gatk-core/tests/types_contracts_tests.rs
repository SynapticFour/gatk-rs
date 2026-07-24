use gatk_core::types::{
    Allele, Base, BaseQuality, GenomicInterval, GenomicPosition, Genotype, ReadQuality,
    SequenceRead, VariantContext, VariantType,
};

#[test]
fn base_roundtrip_and_complement_contract() {
    let cases = [
        ('A', Base::A, Base::T),
        ('C', Base::C, Base::G),
        ('G', Base::G, Base::C),
        ('T', Base::T, Base::A),
        ('N', Base::N, Base::N),
    ];

    for (ch, base, complement) in cases {
        assert_eq!(Base::from_char(ch), Some(base));
        assert_eq!(base.to_char(), ch);
        assert_eq!(base.complement(), complement);
    }
    assert_eq!(Base::from_char('x'), None);
}

#[test]
fn allele_normalization_contract() {
    let allele = Allele::from_string("aTn").expect("valid base symbols should parse");
    assert_eq!(allele.to_string(), "ATN");
    assert_eq!(allele.length(), 3);
    assert!(!allele.is_empty());
    assert!(Allele::from_string("ATX").is_none());
}

#[test]
fn genomic_interval_contains_and_length_contract() {
    let interval = GenomicInterval::new(1, 100, 110);
    assert_eq!(interval.length(), 11);
    assert!(interval.contains(GenomicPosition {
        contig: 1,
        position: 100
    }));
    assert!(interval.contains(GenomicPosition {
        contig: 1,
        position: 110
    }));
    assert!(!interval.contains(GenomicPosition {
        contig: 1,
        position: 111
    }));
    assert!(!interval.contains(GenomicPosition {
        contig: 2,
        position: 105
    }));
}

#[test]
fn variant_type_contract() {
    let pos = GenomicPosition {
        contig: 1,
        position: 1000,
    };
    let ref_a = Allele::from_string("A").unwrap();

    let snp = VariantContext::new(pos, ref_a.clone(), vec![Allele::from_string("T").unwrap()]);
    assert_eq!(snp.variant_type(), VariantType::SNP);

    let ins = VariantContext::new(pos, ref_a.clone(), vec![Allele::from_string("AT").unwrap()]);
    assert_eq!(ins.variant_type(), VariantType::Insertion);

    let del = VariantContext::new(
        pos,
        Allele::from_string("AT").unwrap(),
        vec![Allele::from_string("A").unwrap()],
    );
    assert_eq!(del.variant_type(), VariantType::Deletion);
}

#[test]
fn sequence_read_valid_dna_contract() {
    let position = GenomicPosition {
        contig: 1,
        position: 50,
    };
    let qualities = ReadQuality::new(
        gatk_core::MappingQuality::Score(60),
        vec![BaseQuality::new(30); 4],
    );

    let valid = SequenceRead::new(
        "read-valid".to_string(),
        vec![Base::A, Base::C, Base::G, Base::T],
        qualities.clone(),
        position,
        false,
        false,
    );
    assert!(valid.is_valid_dna());

    let ambiguous = SequenceRead::new(
        "read-ambiguous".to_string(),
        vec![Base::A, Base::N, Base::G, Base::T],
        qualities,
        position,
        false,
        false,
    );
    assert!(!ambiguous.is_valid_dna());
}

#[test]
fn genotype_classification_contract() {
    let hom_ref = Genotype::new(vec![0, 0], 2);
    assert!(hom_ref.is_hom_ref());
    assert!(!hom_ref.is_hom_var());
    assert!(!hom_ref.is_het());

    let hom_var = Genotype::new(vec![1, 1], 2);
    assert!(!hom_var.is_hom_ref());
    assert!(hom_var.is_hom_var());
    assert!(!hom_var.is_het());

    let het = Genotype::new(vec![0, 1], 2);
    assert!(!het.is_hom_ref());
    assert!(!het.is_hom_var());
    assert!(het.is_het());
}
