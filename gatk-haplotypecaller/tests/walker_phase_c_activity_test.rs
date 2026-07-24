//! smoothed activity + binary is_active gates align on fixtures.

use gatk_haplotypecaller::{
    activity_profile::GATK_DEFAULT_ACTIVE_PROB_THRESHOLD, dump_active_locus_tsv,
    dump_smoothed_activity_profile_tsv, format_activity_prob, ReadFilterParams,
};
use std::io::Cursor;

#[test]
fn active_locus_matches_smoothed_threshold_on_p5_snp_fixture() {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../parity/fixtures");
    let ref_fa = root.join("p5_live_reference.fa");
    let bam = root.join("p5_live_case_snp.sam");
    let interval = "chrLive:1-24";
    let filters = ReadFilterParams::default();

    let mut smoothed = Cursor::new(Vec::new());
    use gatk_haplotypecaller::GATK_DEFAULT_ASSEMBLY_REGION_PADDING;
    dump_smoothed_activity_profile_tsv(
        &ref_fa,
        &bam,
        interval,
        GATK_DEFAULT_ASSEMBLY_REGION_PADDING,
        &mut smoothed,
        &filters,
    )
    .unwrap();
    let mut active = Cursor::new(Vec::new());
    dump_active_locus_tsv(
        &ref_fa,
        &bam,
        interval,
        GATK_DEFAULT_ASSEMBLY_REGION_PADDING,
        &mut active,
        &filters,
    )
    .unwrap();

    let smooth_s = String::from_utf8(smoothed.into_inner()).unwrap();
    let active_s = String::from_utf8(active.into_inner()).unwrap();
    let mut smooth_rows = smooth_s.lines().skip(1);
    for line in active_s.lines().skip(1) {
        let Some(sm) = smooth_rows.next() else {
            panic!("fewer smoothed rows than active rows");
        };
        let mut sa = sm.split('\t');
        let contig_s = sa.next().unwrap();
        let pos_s: u64 = sa.next().unwrap().parse().unwrap();
        let prob: f64 = sa.next().unwrap().parse().unwrap();
        let mut aa = line.split('\t');
        assert_eq!(contig_s, aa.next().unwrap());
        assert_eq!(pos_s, aa.next().unwrap().parse::<u64>().unwrap());
        let is_active = aa.next().unwrap() == "true";
        assert_eq!(is_active, prob > GATK_DEFAULT_ACTIVE_PROB_THRESHOLD);
        let _ = format_activity_prob(prob);
    }
    assert!(smooth_rows.next().is_none());
}

#[test]
fn c2_golden_file_matches_dump_chr1_5_15() {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../parity/fixtures");
    let ref_fa = root.join("reference.fa");
    let bam = root.join("sample.bam");
    let golden =
        std::fs::read_to_string(root.join("hc-full-parity/c2/expected/chr1_5_15.sample_bam.tsv"))
            .unwrap();
    let mut actual = Cursor::new(Vec::new());
    use gatk_haplotypecaller::GATK_DEFAULT_ASSEMBLY_REGION_PADDING;
    dump_smoothed_activity_profile_tsv(
        &ref_fa,
        &bam,
        "chr1:5-15",
        GATK_DEFAULT_ASSEMBLY_REGION_PADDING,
        &mut actual,
        &ReadFilterParams::default(),
    )
    .unwrap();
    assert_eq!(String::from_utf8(actual.into_inner()).unwrap(), golden);
}
