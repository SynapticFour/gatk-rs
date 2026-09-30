//! 6R.311 holdout: a generic cross-contig mapped mate, not the canonical read.
//! Skipped unless `HOLDOUT_6R311=1`. Does not run HOLDOUT_6R243.
//!
//! ```text
//! HOLDOUT_6R311=1 cargo test -p gatk-haplotypecaller --test holdout_6r311_failed_mate_production_patch -- --nocapture --test-threads=1
//! ```

use gatk_haplotypecaller::{
    failed_mate_evidence_keys, passes_mate_on_same_contig_or_no_mapped_mate,
};
use rust_htslib::bam::record::{Cigar, CigarString};
use rust_htslib::bam::{HeaderView, Record};
use std::sync::Arc;

fn rec(qname: &str, paired: bool, unmapped: bool, mate_unmapped: bool, mtid: i32) -> Record {
    let mut r = Record::new();
    r.set_header(Arc::new(HeaderView::from_bytes(
        b"@HD\tVN:1.0\n@SQ\tSN:alpha\tLN:1000\n@SQ\tSN:beta\tLN:1000\n",
    )));
    r.set(
        qname.as_bytes(),
        Some(&CigarString::from(vec![Cigar::Match(8)])),
        b"ACGTACGT",
        &vec![20u8; 8],
    );
    r.set_tid(0);
    r.set_mtid(mtid);
    if paired {
        r.set_paired();
    }
    if unmapped {
        r.set_unmapped();
    } else {
        r.unset_unmapped();
    }
    if mate_unmapped {
        r.set_mate_unmapped();
    } else {
        r.unset_mate_unmapped();
    }
    r
}

#[test]
fn holdout_6r311_failed_mate_production_patch() {
    if std::env::var("HOLDOUT_6R311").ok().as_deref() != Some("1") {
        eprintln!("6R311 holdout skipped; set HOLDOUT_6R311=1");
        return;
    }
    let cross = rec("holdout-cross-contig", true, false, false, 1);
    let same = rec("holdout-same-contig", true, false, false, 0);
    let unmapped_mate = rec("holdout-unmapped-mate", true, false, true, 1);
    let unpaired = rec("holdout-unpaired", false, false, false, 1);
    assert!(!passes_mate_on_same_contig_or_no_mapped_mate(&cross));
    assert!(passes_mate_on_same_contig_or_no_mapped_mate(&same));
    assert!(passes_mate_on_same_contig_or_no_mapped_mate(&unmapped_mate));
    assert!(passes_mate_on_same_contig_or_no_mapped_mate(&unpaired));
    assert_ne!(cross.qname(), b"HISEQ1:11:H8GV6ADXX:2:1103:14252:55237");

    let originals = [cross, same, unmapped_mate, unpaired];
    let keys = failed_mate_evidence_keys(&originals);
    assert_eq!(keys.len(), 1);
    let mut evidence = originals.to_vec();
    evidence.retain(|r| !keys.contains(&(r.qname().to_vec(), r.flags())));
    assert_eq!(evidence.len(), 3);
    assert!(evidence
        .iter()
        .all(|r| r.qname() != b"holdout-cross-contig"));
    assert!(evidence.iter().any(|r| r.qname() == b"holdout-same-contig"));
    assert!(evidence
        .iter()
        .any(|r| r.qname() == b"holdout-unmapped-mate"));
    assert!(evidence.iter().any(|r| r.qname() == b"holdout-unpaired"));
    println!("6R311\tclassification\tFAILED_MATE_PRODUCTION_PATCH_REPRODUCES_JAVA_MEMBERSHIP");
}
