use gatk_core::io::fasta::FastaSequence;
use gatk_core::tests::{test_allele, test_position, test_read, test_sequence, test_variant};
use gatk_core::{Base, SequenceRead, VariantType};

#[test]
fn test_helper_allele_construction() {
    let allele = test_allele("AtX");
    assert_eq!(allele.to_string(), "ATN");
    assert_eq!(allele.length(), 3);
}

#[test]
fn test_helper_sequence_construction() {
    let seq = test_sequence("aTcgnx");
    assert_eq!(
        seq,
        vec![Base::A, Base::T, Base::C, Base::G, Base::N, Base::N]
    );
}

#[test]
fn test_variant_helper_semantics() {
    let snp = test_variant("chr2", 42, "A", "T");
    assert_eq!(snp.position, test_position("chr2", 42));
    assert_eq!(snp.variant_type(), VariantType::SNP);
    assert_eq!(snp.reference.to_string(), "A");
    assert_eq!(snp.alternate_alleles[0].to_string(), "T");
}

#[test]
fn test_sequence_read_valid_dna_behavior() {
    let read: SequenceRead = test_read("r1", "ATCG", &[30, 30, 30, 30]);
    assert!(read.is_valid_dna());

    let read_with_ambiguous = test_read("r2", "ATNG", &[30, 30, 30, 30]);
    assert!(!read_with_ambiguous.is_valid_dna());
}

#[test]
fn test_fasta_reverse_complement_and_validation() {
    let seq = FastaSequence::new("chr1".to_string(), b"ATCGNatcgn".to_vec());
    assert!(seq.is_valid_dna());

    let rc = seq.reverse_complement();
    assert_eq!(rc.as_bytes(), b"ncgatNCGAT");
}
