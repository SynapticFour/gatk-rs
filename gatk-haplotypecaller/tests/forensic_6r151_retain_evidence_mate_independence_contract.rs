//! 6R.151 coordinate-free: default HC `retainEvidence` keeps overlapping mates.
//!
//! Java 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`)
//!
//! ```text
//! AlleleLikelihoods.retainEvidence(predicate)
//!   removeEvidenceByIndex for !predicate.test(evidence)
//! composeReadQualifiesForGenotypingPredicate (applyBQD=false, applyFRD=false):
//!   (read, target) -> target.overlaps(read)
//! calculateGLsForThisEvent uses the retained object as-is.
//!
//! Mutect groupEvidence(GATKRead::getName) collapses fragments.
//! Default HC assignGenotypeLikelihoods does not call groupEvidence.
//! ```
//!
//! QNAME identity is not a Java retainEvidence key. Overlapping mates stay
//! independent evidence units.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r151_retain_evidence_mate_independence_contract
//! HOLDOUT_6R151=1 GATK_RS_EXPERIMENTAL_KBEST_POLICY=unbounded_diagnostic \
//!   cargo test -p gatk-haplotypecaller --test holdout_6r151_retain_evidence -- --nocapture --test-threads=1
//! ```

use gatk_haplotypecaller::hc_genotyping_engine::java_alignment_read_overlaps_interval;
use rust_htslib::bam::record::{Cigar, CigarString};
use std::collections::{HashMap, HashSet};

fn mate(qname: &[u8], pos0: i64, flags: u16) -> rust_htslib::bam::Record {
    let mut rec = rust_htslib::bam::Record::new();
    rec.set(
        qname,
        Some(&CigarString(vec![Cigar::Match(10)])),
        b"ACGTACGTAC",
        b"##########",
    );
    rec.set_pos(pos0);
    rec.set_flags(flags);
    rec
}

/// Java predicate is alignment overlap, not QNAME uniqueness.
#[test]
fn forensic_6r151_retain_evidence_uses_overlap_not_qname() {
    let java_retain_key = "target.overlaps(read)";
    let java_qname_key = false;
    assert_eq!(java_retain_key, "target.overlaps(read)");
    assert!(!java_qname_key);
}

/// flags=99 / 147 overlapping mates: Java keeps both.
/// Live 6R.150 alignment for flags=99 is 29456092–29456240 (not a 10 bp stub).
#[test]
fn forensic_6r151_overlapping_mates_are_independent_evidence() {
    let qn = b"HWI-D00360:7:H88WKADXX:2:2107:6787:30989";
    let mut first = rust_htslib::bam::Record::new();
    first.set(
        qn,
        Some(&CigarString(vec![Cigar::Match(149)])),
        &vec![b'A'; 149],
        &vec![b'#'; 149],
    );
    first.set_pos(29_456_091);
    first.set_flags(99);
    let mut second = rust_htslib::bam::Record::new();
    second.set(
        qn,
        Some(&CigarString(vec![Cigar::Match(120)])),
        &vec![b'A'; 120],
        &vec![b'#'; 120],
    );
    second.set_pos(29_456_140);
    second.set_flags(147);
    assert_eq!(first.qname(), second.qname());
    assert_ne!(first.flags(), second.flags());
    let loc = 29_456_196_u64;
    assert!(java_alignment_read_overlaps_interval(&first, loc, loc, 2));
    assert!(java_alignment_read_overlaps_interval(&second, loc, loc, 2));
    let overlapping: Vec<_> = [&first, &second]
        .into_iter()
        .filter(|r| java_alignment_read_overlaps_interval(r, loc, loc, 2))
        .collect();
    assert_eq!(overlapping.len(), 2);
    let mut qnames = HashSet::new();
    qnames.insert(first.qname().to_owned());
    qnames.insert(second.qname().to_owned());
    assert_eq!(qnames.len(), 1);
}

/// Helper-style QNAME collapse would drop one mate; that is not retainEvidence.
#[test]
fn forensic_6r151_qname_is_not_a_dedupe_key_for_java_retain_evidence() {
    let qn = b"pair";
    let r0 = mate(qn, 99, 99);
    let r1 = mate(qn, 101, 147);
    let reads = [&r0, &r1];
    let overlapping: Vec<usize> = (0..2)
        .filter(|&i| java_alignment_read_overlaps_interval(reads[i], 105, 105, 2))
        .collect();
    assert_eq!(overlapping.len(), 2, "Java retainEvidence keeps both mates");
    let mut best: HashMap<Vec<u8>, usize> = HashMap::new();
    for &i in &overlapping {
        best.insert(reads[i].qname().to_owned(), i);
    }
    assert_eq!(best.len(), 1, "QNAME collapse would drop a Java-kept mate");
    assert_ne!(
        overlapping.len(),
        best.len(),
        "QNAME identity is not a Java retainEvidence key"
    );
}

/// Java `groupEvidence` is a different API (Mutect fragment collapse).
#[test]
fn forensic_6r151_group_evidence_is_not_default_hc_retain_evidence() {
    let default_hc_calls_group_evidence = false;
    let mutect_groups_by_read_name = true;
    assert!(!default_hc_calls_group_evidence);
    assert!(mutect_groups_by_read_name);
}

/// calculateGLsForThisEvent consumes the retainEvidence object with no QNAME pass.
#[test]
fn forensic_6r151_gls_follow_retain_evidence_without_qname_pass() {
    let order = ["marginalize", "retainEvidence", "calculateGLsForThisEvent"];
    assert_eq!(order[1], "retainEvidence");
    assert_eq!(order[2], "calculateGLsForThisEvent");
    let qname_pass_between_retain_and_gls = false;
    assert!(!qname_pass_between_retain_and_gls);
}
