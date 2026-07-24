//! Java `callRegion`: given alleles inform trim before `addGivenAlleles` on trimmed assembly.

use gatk_haplotypecaller::assembly_region_trimmer::TrimVariant;
use gatk_haplotypecaller::given_alleles::{given_alleles_to_trim_variants, GatkGivenAllele};

#[test]
fn given_alleles_extend_trim_variants_without_duplicate() {
    let given = vec![GatkGivenAllele {
        contig: "2".into(),
        start_1based: 100,
        end_1based: 100,
        ref_allele: "A".into(),
        alt_alleles: vec!["G".into()],
    }];
    let mut trim = vec![TrimVariant {
        contig: "2".into(),
        start: 100,
        end: 100,
        is_indel: false,
    }];
    given_alleles_to_trim_variants(&given, "2", &mut trim);
    assert_eq!(trim.len(), 1);
    trim.clear();
    given_alleles_to_trim_variants(&given, "2", &mut trim);
    assert_eq!(trim.len(), 1);
    assert!(!trim[0].is_indel);
}

#[test]
fn given_indel_marked_for_trimmer() {
    let given = vec![GatkGivenAllele {
        contig: "2".into(),
        start_1based: 200,
        end_1based: 202,
        ref_allele: "TTC".into(),
        alt_alleles: vec!["T".into()],
    }];
    let mut trim = Vec::new();
    given_alleles_to_trim_variants(&given, "2", &mut trim);
    assert_eq!(trim.len(), 1);
    assert!(trim[0].is_indel);
}
