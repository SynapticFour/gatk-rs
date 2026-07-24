//! Malformed BAM/VCF/FASTA/interval inputs must yield `Err(GatkError)` — never panic.
//! Complements `docs/ARCHITECTURE.md` category (b).

use gatk_common::{GatkConfig, GatkError};
use gatk_core::io::bam::BamReader;
use gatk_core::io::vcf::VcfReader;
use gatk_core::reference::{IntervalSpec, SequenceDictionary};
use gatk_haplotypecaller::fragment_overlap::overlapping_pairs_indices;
use gatk_haplotypecaller::read_header_semantics::ReadHeaderSemantics;
use gatk_haplotypecaller::run_haplotype_caller;
use gatk_haplotypecaller::smith_waterman::{align, SwOverhangStrategy, SwParameters};
use gatk_haplotypecaller::validate_mapped_read_sanity;
use rust_htslib::bam;
use rust_htslib::bam::record::{Cigar, CigarString};
use std::fs;
use std::path::PathBuf;
use tempfile::tempdir;

fn assert_err_no_panic<T, E: std::fmt::Debug>(r: Result<T, E>, label: &str) {
    assert!(r.is_err(), "{label}: expected Err, got Ok");
}

#[test]
fn missing_bam_path_returns_err() {
    let missing = PathBuf::from("/tmp/gatk-rs-definitely-missing-alignment.bam");
    assert_err_no_panic(BamReader::from_file(&missing), "missing bam");
}

#[test]
fn empty_file_as_bam_returns_err() {
    let dir = tempdir().unwrap();
    let bam = dir.path().join("empty.bam");
    fs::write(&bam, b"").unwrap();
    assert_err_no_panic(BamReader::from_file(&bam), "empty bam");
}

#[test]
fn truncated_bam_magic_returns_err() {
    let dir = tempdir().unwrap();
    let bam = dir.path().join("bad_magic.bam");
    fs::write(&bam, b"NOTBAM\x01\x00").unwrap();
    assert_err_no_panic(BamReader::from_file(&bam), "bad bam magic");
}

#[test]
fn missing_reference_path_returns_err() {
    let mut cfg = GatkConfig::new("HaplotypeCaller".to_string());
    cfg.set_reference("/tmp/gatk-rs-missing-ref.fa".to_string());
    cfg.add_input_file("/tmp/gatk-rs-missing.bam".to_string());
    cfg.set_output_vcf(
        tempdir()
            .unwrap()
            .path()
            .join("out.vcf")
            .display()
            .to_string(),
    );
    assert_err_no_panic(run_haplotype_caller(&cfg), "missing reference");
}

#[test]
fn empty_fasta_returns_err_from_hc_run() {
    let dir = tempdir().unwrap();
    let fa = dir.path().join("empty.fa");
    fs::write(&fa, b"").unwrap();
    // Dictionary open may succeed with 0 contigs; HC run must reject.
    let dict = SequenceDictionary::from_fasta_path(&fa);
    if let Ok(d) = &dict {
        assert_eq!(d.contig_count(), 0);
    }
    let bam = dir.path().join("tiny.sam");
    fs::write(
        &bam,
        b"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:10\nread1\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n",
    )
    .unwrap();
    let mut cfg = GatkConfig::new("HaplotypeCaller".to_string());
    cfg.set_reference(fa.display().to_string());
    cfg.add_input_file(bam.display().to_string());
    cfg.set_output_vcf(dir.path().join("out.vcf").display().to_string());
    let err = run_haplotype_caller(&cfg).expect_err("empty fasta must fail cleanly");
    let msg = format!("{err}");
    assert!(
        msg.contains("no sequences")
            || msg.contains("Failed to read reference")
            || matches!(
                err,
                GatkError::Argument { .. } | GatkError::Configuration { .. } | GatkError::Io { .. }
            ),
        "unexpected error shape: {err:?}"
    );
}

#[test]
fn fasta_header_only_no_bases_rejected_or_empty_dict() {
    let dir = tempdir().unwrap();
    let fa = dir.path().join("header_only.fa");
    fs::write(&fa, b">chr1\n").unwrap();
    // Either parse error or empty/zero-length contig — never panic.
    let _ = SequenceDictionary::from_fasta_path(&fa);
}

#[test]
fn missing_vcf_path_returns_err() {
    let missing = PathBuf::from("/tmp/gatk-rs-missing.vcf");
    assert_err_no_panic(VcfReader::from_file(&missing), "missing vcf");
}

#[test]
fn garbage_vcf_body_returns_err_on_read() {
    let dir = tempdir().unwrap();
    let vcf = dir.path().join("garbage.vcf");
    // Minimal header so open succeeds; body row is truncated / invalid.
    fs::write(
        &vcf,
        b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nNOT_A_RECORD\n",
    )
    .unwrap();
    let mut r = VcfReader::from_file(&vcf).expect("header-only open");
    assert_err_no_panic(r.read_all_records(), "garbage vcf body");
}

#[test]
fn empty_vcf_file_returns_err() {
    let dir = tempdir().unwrap();
    let vcf = dir.path().join("empty.vcf");
    fs::write(&vcf, b"").unwrap();
    assert_err_no_panic(VcfReader::from_file(&vcf), "empty vcf");
}

#[test]
fn invalid_interval_syntax_returns_err() {
    assert_err_no_panic(IntervalSpec::parse("chr1:20-10"), "inverted interval");
    assert_err_no_panic(IntervalSpec::parse("chr1:not-a-number"), "bad coords");
}

#[test]
fn unknown_contig_interval_vs_dict_returns_err() {
    let mut d = SequenceDictionary::new();
    d.add_contig("chr1".to_string(), 100);
    let iv = IntervalSpec::parse("chr9:1-10").unwrap();
    assert_err_no_panic(d.validate_interval(&iv), "unknown contig");
}

#[test]
fn malformed_read_lengths_return_err() {
    assert_err_no_panic(validate_mapped_read_sanity(0, 0, 1, 10), "empty bases");
    assert_err_no_panic(validate_mapped_read_sanity(10, 3, 1, 10), "qual mismatch");
    assert_err_no_panic(validate_mapped_read_sanity(10, 10, 20, 10), "negative span");
}

#[test]
fn unsorted_reads_overlap_returns_err_not_panic() {
    let mut r1 = bam::Record::new();
    let mut r2 = bam::Record::new();
    let cigar = CigarString::from(vec![Cigar::Match(5)]);
    r1.set(b"a", Some(&cigar), b"ACGTA", &[30; 5]);
    r1.set_tid(0);
    r1.set_pos(10);
    r2.set(b"b", Some(&cigar), b"ACGTA", &[30; 5]);
    r2.set_tid(0);
    r2.set_pos(1); // earlier than r1 → unsorted
    let err = overlapping_pairs_indices(&[r1, r2]).expect_err("unsorted must Err");
    assert!(matches!(err, GatkError::Read { .. }));
}

#[test]
fn smith_waterman_empty_input_returns_err_not_panic() {
    assert_err_no_panic(
        align(
            b"",
            b"ACGT",
            &SwParameters::gatk_haplotype_to_reference(),
            SwOverhangStrategy::Indel,
        ),
        "empty ref",
    );
}

#[test]
fn broken_sam_header_rg_validation_returns_err() {
    // @RG without ID is invalid for ReadHeaderSemantics.
    let header = "@HD\tVN:1.6\n@RG\tSM:only_sample\n";
    let res = ReadHeaderSemantics::from_sam_header_text(header);
    assert_err_no_panic(res, "rg missing ID");
}

#[test]
fn hc_run_missing_alignment_input_returns_err() {
    let dir = tempdir().unwrap();
    let fa = dir.path().join("ref.fa");
    fs::write(&fa, b">chr1\nACGTACGTAC\n").unwrap();
    let mut cfg = GatkConfig::new("HaplotypeCaller".to_string());
    cfg.set_reference(fa.display().to_string());
    // no -I
    cfg.set_output_vcf(dir.path().join("out.vcf").display().to_string());
    assert_err_no_panic(run_haplotype_caller(&cfg), "missing -I");
}

#[test]
fn hc_run_empty_alignment_file_returns_err() {
    let dir = tempdir().unwrap();
    let fa = dir.path().join("ref.fa");
    fs::write(&fa, b">chr1\nACGTACGTAC\n").unwrap();
    let bam = dir.path().join("empty_records.sam");
    // Header-only SAM (no alignments).
    fs::write(&bam, b"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:10\n").unwrap();
    let mut cfg = GatkConfig::new("HaplotypeCaller".to_string());
    cfg.set_reference(fa.display().to_string());
    cfg.add_input_file(bam.display().to_string());
    cfg.set_output_vcf(dir.path().join("out.vcf").display().to_string());
    let err = run_haplotype_caller(&cfg).expect_err("empty alignment");
    let msg = format!("{err}");
    assert!(
        msg.contains("no records") || msg.contains("Alignment"),
        "unexpected: {err:?}"
    );
}
