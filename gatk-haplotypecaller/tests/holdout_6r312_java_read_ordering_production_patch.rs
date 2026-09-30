//! 6R.312 holdout: forward-before-reverse at one clipped start, not the
//! canonical 29455560 pair. Skipped unless `HOLDOUT_6R312=1`.
//! Does not run HOLDOUT_6R243.
//!
//! ```text
//! HOLDOUT_6R312=1 cargo test -p gatk-haplotypecaller --test holdout_6r312_java_read_ordering_production_patch -- --nocapture --test-threads=1
//! ```

use gatk_haplotypecaller::assembly_region_finalize::{
    clip_finalized_reads_in_place, java_read_coordinate_compare,
    set_forensic_6r308_java_read_coordinate_order,
};
use gatk_haplotypecaller::{AssemblyRegion, GenomePosition};
use rust_htslib::bam::record::{Cigar, CigarString};
use rust_htslib::bam::{HeaderView, Record};
use std::sync::Arc;

fn rec(name: &str, flags: u16) -> Record {
    let mut r = Record::new();
    r.set_header(Arc::new(HeaderView::from_bytes(
        b"@HD\tVN:1.0\n@SQ\tSN:chrZ\tLN:3000\n",
    )));
    r.set(
        name.as_bytes(),
        Some(&CigarString::from(vec![Cigar::Match(24)])),
        b"GGGGGGGGGGGGGGGGGGGGGGGG",
        &vec![25u8; 24],
    );
    r.set_tid(0);
    r.set_pos(799);
    r.set_flags(flags);
    r.set_mapq(15);
    r
}

#[test]
fn holdout_6r312_java_read_ordering_production_patch() {
    if std::env::var("HOLDOUT_6R312").ok().as_deref() != Some("1") {
        eprintln!("6R312 holdout skipped; set HOLDOUT_6R312=1");
        return;
    }
    set_forensic_6r308_java_read_coordinate_order(true);
    let forward = rec("zeta-forward", 99);
    let reverse = rec("alpha-reverse", 83);
    assert_ne!(forward.qname(), b"HISEQ1:11:H8GV6ADXX:1:1116:4033:65919");
    assert_ne!(reverse.qname(), b"HISEQ1:11:H8GV6ADXX:1:1109:9994:7054");
    assert_eq!(forward.pos(), reverse.pos());
    assert!(java_read_coordinate_compare(&forward, &reverse) < 0);
    let region = AssemblyRegion {
        contig: "chrZ".into(),
        start: GenomePosition::new_1based(700),
        end: GenomePosition::new_1based(900),
        is_active: true,
        extended_start: GenomePosition::new_1based(1),
        extended_end: GenomePosition::new_1based(3000),
        extension: 50,
        reads: Vec::new(),
        read_qnames: Vec::new(),
        reference: gatk_haplotypecaller::ReferenceContext::empty(),
        features: gatk_haplotypecaller::FeatureContext::empty(),
        pileup_loci: Vec::new(),
    };
    let mut name_order = vec![reverse.clone(), forward.clone()];
    name_order.sort_by(|a, b| a.qname().cmp(b.qname()));
    assert_eq!(name_order[0].qname(), b"alpha-reverse");
    let mut production = vec![reverse, forward];
    clip_finalized_reads_in_place(&mut production, &region);
    assert_eq!(production[0].qname(), b"zeta-forward");
    assert_eq!(production[0].flags(), 99);
    assert_eq!(production[1].qname(), b"alpha-reverse");
    assert_eq!(production[1].flags(), 83);
    println!("6R312\tclassification\tJAVA_READ_ORDERING_PRODUCTION_PATCH_REPRODUCES_JAVA");
}
