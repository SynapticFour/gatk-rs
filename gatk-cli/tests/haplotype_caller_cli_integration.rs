//! Phase 9 (112): CLI integration tests for HaplotypeCaller argument surface and output modes.

use assert_cmd::Command;
use std::fs;
use std::path::Path;

fn repo_root() -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

#[test]
fn haplotype_caller_help_exits_zero() {
    Command::cargo_bin("gatk-rs")
        .unwrap()
        .args(["HaplotypeCaller", "--help"])
        .assert()
        .success()
        .stdout(predicates::str::contains(
            "Call germline SNPs and indels via local re-assembly of haplotypes",
        ));
}

#[test]
fn haplotype_caller_emit_ref_confidence_long_form_accepted() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.vcf");
    Command::cargo_bin("gatk-rs")
        .unwrap()
        .current_dir(repo_root())
        .env("GATK_RS_HC_SCAFFOLD_OUTPUT", "1")
        .args([
            "HaplotypeCaller",
            "-R",
            "parity/fixtures/reference.fa",
            "-I",
            "parity/fixtures/sample.bam",
            "-O",
            out.to_str().unwrap(),
            "-L",
            "chr1:1-32",
            "--emit-ref-confidence",
            "NONE",
        ])
        .assert()
        .success();
    let body = fs::read_to_string(&out).unwrap();
    assert!(body.contains("##GATK_RS_HC_PIPELINE=scaffold-v1"));
}

#[test]
fn haplotype_caller_accepts_repeated_interval_flags() {
    // GIAB ci-subset / Java GATK pass many `-L` tokens; clap must Append, not reject.
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.vcf");
    Command::cargo_bin("gatk-rs")
        .unwrap()
        .current_dir(repo_root())
        .env("GATK_RS_HC_SCAFFOLD_OUTPUT", "1")
        .args([
            "HaplotypeCaller",
            "-R",
            "parity/fixtures/reference.fa",
            "-I",
            "parity/fixtures/sample.bam",
            "-O",
            out.to_str().unwrap(),
            "-L",
            "chr1:1-16",
            "-L",
            "chr1:17-32",
        ])
        .assert()
        .success();
    assert!(out.exists());
}

#[test]
fn haplotype_caller_invalid_interval_user_error_exit_2() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("bad.vcf");
    Command::cargo_bin("gatk-rs")
        .unwrap()
        .current_dir(repo_root())
        .args([
            "HaplotypeCaller",
            "-R",
            "parity/fixtures/reference.fa",
            "-I",
            "parity/fixtures/sample.bam",
            "-O",
            out.to_str().unwrap(),
            "-L",
            "chr999:1-10",
        ])
        .assert()
        .code(2);
}

#[test]
fn haplotype_caller_scaffold_output_matches_golden() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("hc.vcf");
    Command::cargo_bin("gatk-rs")
        .unwrap()
        .current_dir(repo_root())
        .env("GATK_RS_HC_SCAFFOLD_OUTPUT", "1")
        .args([
            "HaplotypeCaller",
            "-R",
            "parity/fixtures/reference.fa",
            "-I",
            "parity/fixtures/sample.bam",
            "-O",
            out.to_str().unwrap(),
            "-L",
            "chr1:1-32",
        ])
        .assert()
        .success();

    let got = fs::read_to_string(&out).unwrap();
    let golden =
        fs::read_to_string(repo_root().join("parity/expected/p9_hc_scaffold_golden.vcf")).unwrap();
    assert_eq!(got, golden, "scaffold VCF must match frozen golden header");
}

#[test]
fn haplotype_caller_default_pipeline_is_assembly_region_v1() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("hc.vcf");
    Command::cargo_bin("gatk-rs")
        .unwrap()
        .current_dir(repo_root())
        .args([
            "HaplotypeCaller",
            "-R",
            "parity/fixtures/reference.fa",
            "-I",
            "parity/fixtures/sample.bam",
            "-O",
            out.to_str().unwrap(),
            "-L",
            "chr1:1-32",
        ])
        .assert()
        .success();
    let body = fs::read_to_string(&out).unwrap();
    assert!(body.contains("##GATK_RS_HC_PIPELINE=assembly-region-v1"));
}
