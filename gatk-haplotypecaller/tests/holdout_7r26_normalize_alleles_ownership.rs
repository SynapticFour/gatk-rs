//! 7R.26 holdout: `left_align_indels` consumes the `IndexRange` values
//! `normalize_alleles` already wrote. Coordinate-free. Skipped unless
//! `HOLDOUT_7R26=1`.
//!
//! The composed CIGAR is a short homopolymer analogue of the 7R.24 `1I1D`
//! block that is fully shifted through a match run. The old temporary slice
//! kept the pre-call ranges and emitted `3D1I19M2I2M`. Preserving the
//! normalized ranges emits `5M2D13M2I2M`.
//!
//! ```text
//! HOLDOUT_7R26=1 cargo test -p gatk-haplotypecaller --test holdout_7r26_normalize_alleles_ownership -- --nocapture --test-threads=1
//! ```

use gatk_haplotypecaller::cigar::{Cigar, CigarOperator};
use gatk_haplotypecaller::haplotype_cigar::left_align_indels_for_read;

#[test]
fn normalized_ranges_survive_into_the_cigar() {
    if std::env::var_os("HOLDOUT_7R26").is_none() {
        eprintln!("HOLDOUT_7R26 unset; not run");
        return;
    }
    let reference = b"ATAGTGAGAGAGAGAGAGCTAGCAGA";
    let read = b"ATAGTGAGAGAGAGAGCTAGCAGA";
    let mut composed = Cigar::new();
    for (length, operator) in [
        (5, CigarOperator::Match),
        (2, CigarOperator::Deletion),
        (8, CigarOperator::Match),
        (1, CigarOperator::Insertion),
        (1, CigarOperator::Deletion),
        (4, CigarOperator::Match),
        (2, CigarOperator::Insertion),
        (2, CigarOperator::Match),
    ] {
        composed.push(length, operator);
    }
    let left = left_align_indels_for_read(&composed, reference, read, 0);
    assert_eq!(left.cigar.to_gatk_string(), "5M2D13M2I2M");
    assert_ne!(left.cigar.to_gatk_string(), "3D1I19M2I2M");
}
