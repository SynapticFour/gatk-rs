use gatk_haplotypecaller::{
    passes_hc_read_filters_fields, query_index_at_reference_position,
    reference_position_at_query_index, ReadFilterParams, FLAG_DUPLICATE, FLAG_NOT_PRIMARY,
    FLAG_SEGMENT_UNMAPPED, FLAG_SUPPLEMENTARY,
};
use rust_htslib::bam::record::{Cigar, CigarString};

fn cig(v: Vec<Cigar>) -> CigarString {
    CigarString(v)
}

#[test]
fn hc_filter_contract_primary_and_mapq() {
    let p = ReadFilterParams {
        min_mapping_quality: 20,
        exclude_duplicates: true,
        exclude_secondary: true,
        exclude_supplementary: true,
    };
    assert!(passes_hc_read_filters_fields(0, 30, &p));
    assert!(passes_hc_read_filters_fields(0, 255, &p));
    assert!(!passes_hc_read_filters_fields(0, 19, &p));
}

#[test]
fn hc_filter_contract_flag_exclusions() {
    let p = ReadFilterParams::default();
    for flags in [
        FLAG_DUPLICATE,
        FLAG_NOT_PRIMARY,
        FLAG_SUPPLEMENTARY,
        FLAG_SEGMENT_UNMAPPED,
    ] {
        assert!(!passes_hc_read_filters_fields(flags, 60, &p));
    }
}

#[test]
fn projection_contract_match_blocks() {
    let c = cig(vec![Cigar::Match(5), Cigar::Ins(1), Cigar::Match(5)]);
    let start = 100i64;
    assert_eq!(reference_position_at_query_index(start, &c, 0), Some(100));
    assert_eq!(reference_position_at_query_index(start, &c, 5), None);
    assert_eq!(reference_position_at_query_index(start, &c, 6), Some(105));
    assert_eq!(query_index_at_reference_position(start, &c, 104), Some(4));
    assert_eq!(query_index_at_reference_position(start, &c, 105), Some(6));
}
