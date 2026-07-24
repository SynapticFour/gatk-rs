//! N3: reference haplotype preserved through allele filtering.

use gatk_haplotypecaller::allele_filtering::ensure_reference_haplotype;
use gatk_haplotypecaller::haplotype::Haplotype;

#[test]
fn ensure_reference_restores_flag() {
    let mut haps = vec![
        Haplotype::new(b"ACGT", false),
        Haplotype::new(b"ACGT", false),
    ];
    ensure_reference_haplotype(&mut haps);
    assert!(haps[0].is_reference);
}
