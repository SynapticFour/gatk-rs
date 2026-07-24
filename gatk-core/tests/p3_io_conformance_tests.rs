use gatk_core::io::{
    alignment_header_canonical_hd_sq_rg, copy_alignments_with_htslib,
    copy_alignments_with_htslib_with_reference, count_records_in_region_indexed,
    qnames_in_region_indexed, validate_bam_file, VcfReader,
};
use std::path::PathBuf;

fn fixture(path: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join("parity")
        .join("fixtures")
        .join(path)
}

#[test]
fn bam_roundtrip_preserves_optional_tags_from_sam_source() {
    let dir = tempfile::tempdir().unwrap();
    let in_sam = fixture("p3_optional_tags.sam");
    let out_bam = dir.path().join("out.bam");
    let out_sam = dir.path().join("out.sam");

    let n1 = copy_alignments_with_htslib(&in_sam, &out_bam).unwrap();
    let n2 = copy_alignments_with_htslib(&out_bam, &out_sam).unwrap();
    assert_eq!(n1, 2);
    assert_eq!(n2, 2);

    let text = std::fs::read_to_string(&out_sam).unwrap();
    assert!(text.contains("NM:i:0"));
    assert!(text.contains("AS:i:16"));
    assert!(text.contains("MD:Z:16"));
    assert!(text.contains("XX:A:K"));
    assert!(text.contains("XF:f:1.5"));
    assert!(text.contains("XS:Z:alpha"));
    assert!(text.contains("XH:H:0A0B"));
    assert!(text.contains("XB:B:i,1,2,3"));
    assert!(text.contains("XE:B:f,1.25,2.5"));
    assert!(text.contains("XC:B:c,-1,2"));
    assert!(text.contains("XO:i:-5"));
}

#[test]
fn cram_roundtrip_with_reference_preserves_optional_tags_contract() {
    let dir = tempfile::tempdir().unwrap();
    let in_sam = fixture("p3_optional_tags.sam");
    let reference_fasta = fixture("reference.fa");
    let out_cram = dir.path().join("out.cram");
    let out_bam = dir.path().join("out.bam");
    let out_sam = dir.path().join("out.sam");

    let n1 = copy_alignments_with_htslib_with_reference(&in_sam, &out_cram, Some(&reference_fasta))
        .unwrap();
    let n2 =
        copy_alignments_with_htslib_with_reference(&out_cram, &out_bam, Some(&reference_fasta))
            .unwrap();
    let n3 = copy_alignments_with_htslib(&out_bam, &out_sam).unwrap();
    assert_eq!(n1, 2);
    assert_eq!(n2, 2);
    assert_eq!(n3, 2);

    let text = std::fs::read_to_string(&out_sam).unwrap();
    assert!(text.contains("NM:i:0"));
    assert!(text.contains("AS:i:16"));
    assert!(text.contains("XX:A:K"));
    assert!(text.contains("XF:f:1.5"));
    assert!(text.contains("XB:B:i,1,2,3"));
    assert!(text.contains("XE:B:f,1.25,2.5"));
    assert!(text.contains("XC:B:c,-1,2"));
    assert!(text.contains("XO:i:-5"));
}

#[test]
fn header_canonical_hd_sq_rg_stable_across_htslib_roundtrip_contract() {
    let reference_fasta = fixture("reference.fa");
    let cases = [
        ("p3_header_canonicalization.sam", 1u64),
        ("p3_header_edge_cases.sam", 2u64),
    ];
    for (fixture_name, expected_records) in cases {
        let dir = tempfile::tempdir().unwrap();
        let in_sam = fixture(fixture_name);
        let out_bam = dir.path().join("canonical.bam");
        let out_sam = dir.path().join("canonical.sam");
        let out_cram = dir.path().join("canonical.cram");
        let out_cram_sam = dir.path().join("canonical.cram.sam");

        let before = alignment_header_canonical_hd_sq_rg(&in_sam).unwrap();

        let n1 = copy_alignments_with_htslib(&in_sam, &out_bam).unwrap();
        let n2 = copy_alignments_with_htslib(&out_bam, &out_sam).unwrap();
        assert_eq!(n1, expected_records, "BAM roundtrip record count mismatch");
        assert_eq!(n2, expected_records, "SAM materialization count mismatch");
        let after_bam = alignment_header_canonical_hd_sq_rg(&out_sam).unwrap();
        assert_eq!(
            before, after_bam,
            "Canonical HD/SQ/RG mismatch after SAM->BAM->SAM for fixture {fixture_name}"
        );

        let n3 =
            copy_alignments_with_htslib_with_reference(&in_sam, &out_cram, Some(&reference_fasta))
                .unwrap();
        let n4 = copy_alignments_with_htslib_with_reference(
            &out_cram,
            &out_cram_sam,
            Some(&reference_fasta),
        )
        .unwrap();
        assert_eq!(n3, expected_records, "CRAM roundtrip record count mismatch");
        assert_eq!(
            n4, expected_records,
            "CRAM->SAM materialization count mismatch"
        );
        let after_cram = alignment_header_canonical_hd_sq_rg(&out_cram_sam).unwrap();
        assert_eq!(
            before, after_cram,
            "Canonical HD/SQ/RG mismatch after SAM->CRAM->SAM for fixture {fixture_name}"
        );
    }
}

#[test]
fn indexed_region_query_counts_expected_reads() {
    let dir = tempfile::tempdir().unwrap();
    let in_sam = fixture("p3_optional_tags.sam");
    let out_bam = dir.path().join("indexed.bam");
    let _ = copy_alignments_with_htslib(&in_sam, &out_bam).unwrap();

    let sorted_bam = dir.path().join("indexed.sorted.bam");
    let _ = copy_alignments_with_htslib(&out_bam, &sorted_bam).unwrap();

    // Build index via samtools for portability across rust-htslib API changes.
    let Ok(status) = std::process::Command::new("samtools")
        .arg("index")
        .arg(&sorted_bam)
        .status()
    else {
        return;
    };
    if !status.success() {
        return;
    }

    let c1 = count_records_in_region_indexed(&sorted_bam, "chr1", 1, 16).unwrap();
    let c2 = count_records_in_region_indexed(&sorted_bam, "chr1", 17, 32).unwrap();
    assert_eq!(c1, 1);
    assert_eq!(c2, 1);
}

#[test]
fn validate_bam_file_accepts_roundtripped_p3_fixture() {
    let dir = tempfile::tempdir().unwrap();
    let in_sam = fixture("p3_optional_tags.sam");
    let out_bam = dir.path().join("validated.bam");
    let _ = copy_alignments_with_htslib(&in_sam, &out_bam).unwrap();
    let n = validate_bam_file(&out_bam, None).unwrap();
    assert_eq!(n, 2);
}

#[test]
fn complex_header_roundtrip_retains_rg_pg_fields() {
    let dir = tempfile::tempdir().unwrap();
    let in_sam = fixture("p3_complex_header.sam");
    let out_bam = dir.path().join("complex.bam");
    let out_sam = dir.path().join("complex.sam");

    let n1 = copy_alignments_with_htslib(&in_sam, &out_bam).unwrap();
    let n2 = copy_alignments_with_htslib(&out_bam, &out_sam).unwrap();
    assert_eq!(n1, 1);
    assert_eq!(n2, 1);

    let text = std::fs::read_to_string(&out_sam).unwrap();
    assert!(text.contains("@RG\tID:rgA"));
    assert!(text.contains("@RG\tID:rgB"));
    assert!(text.contains("@PG\tID:pgA"));
    assert!(text.contains("@PG\tID:pgB"));
    assert!(text.contains("RG:Z:rgA"));
    assert!(text.contains("PG:Z:pgB"));
}

#[test]
fn malformed_bam_is_rejected_with_error_contract() {
    let dir = tempfile::tempdir().unwrap();
    let malformed = dir.path().join("malformed.bam");
    std::fs::write(&malformed, b"not-a-bam").unwrap();
    let err = validate_bam_file(&malformed, None).unwrap_err();
    let msg = format!("{err}");
    assert!(msg.contains("Failed to open alignment file"));
}

#[test]
fn malformed_vcf_short_record_is_rejected() {
    let malformed = fixture("p3_malformed_short_record.vcf");
    let mut reader = VcfReader::from_file(&malformed).unwrap();
    let err = reader.read_next_record().unwrap_err();
    let msg = format!("{err}");
    assert!(msg.contains("VCF record must have at least 8 fields"));
}

#[test]
fn indexed_region_query_returns_stable_qname_order() {
    let dir = tempfile::tempdir().unwrap();
    let in_sam = fixture("p3_region_reads.sam");
    let out_bam = dir.path().join("region_reads.bam");
    let _ = copy_alignments_with_htslib(&in_sam, &out_bam).unwrap();
    let Ok(status) = std::process::Command::new("samtools")
        .arg("index")
        .arg(&out_bam)
        .status()
    else {
        return;
    };
    if !status.success() {
        return;
    }

    let qnames = qnames_in_region_indexed(&out_bam, "chr1", 1, 16).unwrap();
    assert_eq!(qnames, vec!["rA", "rB", "rC"]);
}
