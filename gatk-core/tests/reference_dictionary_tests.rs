use gatk_core::io::{BamHeader, Contig, ReferenceSequence, SamHeader, VcfHeader};
use gatk_core::reference::{parse_intervals_cli_string, IntervalSpec, SequenceDictionary};
use std::fs;
use std::path::PathBuf;

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("src")
        .join("tests")
        .join("test_data")
        .join(name)
}

#[test]
fn dictionary_loads_reference_and_validates_interval() {
    let dict = SequenceDictionary::from_fasta_path(fixture("reference.fa")).unwrap();
    assert!(dict.contig_count() > 0);

    let iv = IntervalSpec::parse("chr1:1-10").unwrap();
    assert!(dict.validate_interval(&iv).is_ok());
}

#[test]
fn dictionary_rejects_unknown_contig() {
    let dict = SequenceDictionary::from_fasta_path(fixture("reference.fa")).unwrap();
    let iv = IntervalSpec::parse("chr999:1-10").unwrap();
    assert!(dict.validate_interval(&iv).is_err());
}

#[test]
fn parses_interval_list_file() {
    let temp = std::env::temp_dir().join("gatk_rs_intervals.list");
    fs::write(&temp, "chr1:1-10\n# comment\nchr1:20-30\n").unwrap();
    let intervals = IntervalSpec::parse_list_file(&temp).unwrap();
    assert_eq!(intervals.len(), 2);
    let _ = fs::remove_file(temp);
}

#[test]
fn validates_vcf_and_bam_headers_against_dictionary() {
    let mut dict = SequenceDictionary::new();
    dict.add_contig("chr1".to_string(), 1000);

    let mut vcf_header = VcfHeader::default();
    vcf_header.contigs.push(Contig {
        id: "chr1".to_string(),
        length: Some(1000),
        md5: None,
        assembly: None,
        species: None,
        uri: None,
    });
    assert!(dict.validate_vcf_header(&vcf_header).is_ok());

    let mut bam_header = BamHeader::default();
    bam_header.reference_sequences.push(ReferenceSequence {
        name: "chr1".to_string(),
        length: 1000,
        md5: None,
        assembly: None,
        uri: None,
        species: None,
    });
    assert!(dict.validate_bam_header(&bam_header).is_ok());

    let mut sam_header = SamHeader::default();
    sam_header.reference_sequences.push(ReferenceSequence {
        name: "chr1".to_string(),
        length: 1000,
        md5: None,
        assembly: None,
        uri: None,
        species: None,
    });
    assert!(dict.validate_sam_header(&sam_header).is_ok());
}

#[test]
fn parses_picard_style_interval_list_file() {
    let temp = std::env::temp_dir().join("gatk_rs_picard.interval_list");
    fs::write(
        &temp,
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\nchr1\t10\t20\t+\tname1\n",
    )
    .unwrap();
    let intervals = IntervalSpec::parse_list_file(&temp).unwrap();
    assert_eq!(intervals.len(), 1);
    assert_eq!(intervals[0].contig, "chr1");
    assert_eq!(intervals[0].start, Some(10));
    assert_eq!(intervals[0].end, Some(20));
    let _ = fs::remove_file(temp);
}

#[test]
fn rejects_interval_list_sq_length_mismatch_against_dictionary() {
    let temp = std::env::temp_dir().join("gatk_rs_picard_bad_sq.interval_list");
    fs::write(
        &temp,
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:999\nchr1\t10\t20\t+\tname1\n",
    )
    .unwrap();

    let mut dict = SequenceDictionary::new();
    dict.add_contig("chr1".to_string(), 1000);
    let err = IntervalSpec::parse_list_file_with_dictionary(&temp, &dict).unwrap_err();
    assert!(format!("{err}").contains("length mismatch"));

    let _ = fs::remove_file(temp);
}

#[test]
fn rejects_interval_list_sq_missing_sn() {
    let temp = std::env::temp_dir().join("gatk_rs_picard_missing_sn.interval_list");
    fs::write(
        &temp,
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tLN:1000\nchr1\t10\t20\t+\tname1\n",
    )
    .unwrap();

    let mut dict = SequenceDictionary::new();
    dict.add_contig("chr1".to_string(), 1000);
    let err = IntervalSpec::parse_list_file_with_dictionary(&temp, &dict).unwrap_err();
    assert!(format!("{err}").contains("missing SN"));

    let _ = fs::remove_file(temp);
}

#[test]
fn rejects_interval_list_sq_missing_ln() {
    let temp = std::env::temp_dir().join("gatk_rs_picard_missing_ln.interval_list");
    fs::write(
        &temp,
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\nchr1\t10\t20\t+\tname1\n",
    )
    .unwrap();

    let mut dict = SequenceDictionary::new();
    dict.add_contig("chr1".to_string(), 1000);
    let err = IntervalSpec::parse_list_file_with_dictionary(&temp, &dict).unwrap_err();
    assert!(format!("{err}").contains("missing LN"));

    let _ = fs::remove_file(temp);
}

#[test]
fn rejects_interval_list_sq_duplicate_sn() {
    let temp = std::env::temp_dir().join("gatk_rs_picard_dup_sn.interval_list");
    fs::write(
        &temp,
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tSN:chr1\tLN:1000\nchr1\t10\t20\t+\tname1\n",
    )
    .unwrap();

    let mut dict = SequenceDictionary::new();
    dict.add_contig("chr1".to_string(), 1000);
    let err = IntervalSpec::parse_list_file_with_dictionary(&temp, &dict).unwrap_err();
    assert!(format!("{err}").contains("duplicate SN"));

    let _ = fs::remove_file(temp);
}

#[test]
fn rejects_interval_list_sq_non_key_value_field() {
    let temp = std::env::temp_dir().join("gatk_rs_picard_bad_sq_field.interval_list");
    fs::write(
        &temp,
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\tBADFIELD\nchr1\t10\t20\t+\tname1\n",
    )
    .unwrap();

    let mut dict = SequenceDictionary::new();
    dict.add_contig("chr1".to_string(), 1000);
    let err = IntervalSpec::parse_list_file_with_dictionary(&temp, &dict).unwrap_err();
    assert!(format!("{err}").contains("KEY:VALUE"));

    let _ = fs::remove_file(temp);
}

#[test]
fn rejects_interval_list_invalid_strand_value() {
    let temp = std::env::temp_dir().join("gatk_rs_picard_bad_strand.interval_list");
    fs::write(
        &temp,
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\nchr1\t10\t20\tX\tname1\n",
    )
    .unwrap();

    let mut dict = SequenceDictionary::new();
    dict.add_contig("chr1".to_string(), 1000);
    let err = IntervalSpec::parse_list_file_with_dictionary(&temp, &dict).unwrap_err();
    assert!(format!("{err}").contains("strand"));

    let _ = fs::remove_file(temp);
}

#[test]
fn rejects_interval_list_row_outside_contig_bounds() {
    let temp = std::env::temp_dir().join("gatk_rs_picard_out_of_bounds.interval_list");
    fs::write(
        &temp,
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\nchr1\t995\t1005\t+\tname1\n",
    )
    .unwrap();

    let mut dict = SequenceDictionary::new();
    dict.add_contig("chr1".to_string(), 1000);
    let err = IntervalSpec::parse_list_file_with_dictionary(&temp, &dict).unwrap_err();
    assert!(format!("{err}").contains("exceeds contig"));

    let _ = fs::remove_file(temp);
}

#[test]
fn accepts_interval_list_with_chr_alias_against_dictionary() {
    let temp = std::env::temp_dir().join("gatk_rs_picard_alias.interval_list");
    fs::write(
        &temp,
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\nchr1\t10\t20\t+\tname1\n",
    )
    .unwrap();

    let mut dict = SequenceDictionary::new();
    dict.add_contig("1".to_string(), 1000);
    let intervals = IntervalSpec::parse_list_file_with_dictionary(&temp, &dict).unwrap();
    assert_eq!(intervals.len(), 1);

    let _ = fs::remove_file(temp);
}

#[test]
fn interval_boundary_edge_cases_match_1_based_contract() {
    let mut dict = SequenceDictionary::new();
    dict.add_contig("chr1".to_string(), 10);

    let at_start = IntervalSpec::parse("chr1:1-1").unwrap();
    let at_end = IntervalSpec::parse("chr1:10-10").unwrap();
    assert!(dict.validate_interval(&at_start).is_ok());
    assert!(dict.validate_interval(&at_end).is_ok());

    assert!(IntervalSpec::parse("chr1:0-1").is_err());
    assert!(IntervalSpec::parse("chr1:5-4").is_err());
}

#[test]
fn interval_compatibility_rejects_exclusion_only_expression() {
    let mut dict = SequenceDictionary::new();
    dict.add_contig("chr1".to_string(), 100);
    let err = parse_intervals_cli_string(&dict, "^chr1:1-10").unwrap_err();
    assert!(format!("{err}").contains("only -L token"));
}
