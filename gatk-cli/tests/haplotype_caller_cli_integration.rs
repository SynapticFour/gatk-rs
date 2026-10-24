//! Phase 9 (112): CLI integration tests for HaplotypeCaller argument surface and output modes.

use assert_cmd::Command;
use std::fs;
use std::path::Path;

fn repo_root() -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

#[test]
fn haplotype_caller_help_pair_hmm_omitted_default_is_fastest_available() {
    let asserted = Command::cargo_bin("gatk-rs")
        .unwrap()
        .args(["HaplotypeCaller", "--help"])
        .assert()
        .success();
    let help = String::from_utf8(asserted.get_output().stdout.clone()).unwrap();
    let pair_hmm = help
        .lines()
        .skip_while(|line| !line.contains("--pair-hmm"))
        .take(14)
        .collect::<Vec<_>>()
        .join("\n");
    assert!(
        pair_hmm.contains("FASTEST_AVAILABLE"),
        "omitted --pair-hmm must be described as FASTEST_AVAILABLE:\n{pair_hmm}"
    );
    assert!(
        pair_hmm.contains("LOG10_PAIRHMM"),
        "LOG10_PAIRHMM must remain an explicit option:\n{pair_hmm}"
    );
    assert!(
        !pair_hmm.contains("LOG10_PAIRHMM` (default)")
            && !pair_hmm.contains("LOG10_PAIRHMM (default)"),
        "LOG10_PAIRHMM must not be described as the default:\n{pair_hmm}"
    );
}

#[test]
fn haplotype_caller_plain_vcf_is_text_and_compressed_suffixes_are_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let plain = dir.path().join("sample.vcf");
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
            plain.to_str().unwrap(),
            "-L",
            "chr1:1-32",
        ])
        .assert()
        .success();
    let bytes = fs::read(&plain).unwrap();
    assert!(bytes.starts_with(b"##fileformat=VCF"));
    assert!(!bytes.starts_with(&[0x1f, 0x8b]));

    for name in ["sample.vcf.gz", "sample.bgz", "sample.bgzip", "sample.bcf"] {
        let seeded = dir.path().join(name);
        fs::write(&seeded, b"KEEP").unwrap();
        let asserted = Command::cargo_bin("gatk-rs")
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
                seeded.to_str().unwrap(),
                "-L",
                "chr1:1-32",
            ])
            .assert()
            .code(2);
        let stderr = String::from_utf8(asserted.get_output().stderr.clone()).unwrap();
        assert!(
            stderr.contains("unsupported output format"),
            "{name}: {stderr}"
        );
        assert!(stderr.contains("plain-text VCF"), "{name}: {stderr}");
        assert!(stderr.contains("not implemented"), "{name}: {stderr}");
        assert_eq!(fs::read(&seeded).unwrap(), b"KEEP", "{name} was modified");

        let absent = dir.path().join(format!("absent-{name}"));
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
                absent.to_str().unwrap(),
                "-L",
                "chr1:1-32",
            ])
            .assert()
            .code(2);
        assert!(!absent.exists(), "{name} was created");
    }
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

fn hc_fixture_args(out: &Path) -> Vec<String> {
    vec![
        "HaplotypeCaller".to_string(),
        "-R".to_string(),
        "parity/fixtures/reference.fa".to_string(),
        "-I".to_string(),
        "parity/fixtures/sample.bam".to_string(),
        "-O".to_string(),
        out.to_str().unwrap().to_string(),
        "-L".to_string(),
        "chr1:1-32".to_string(),
    ]
}

fn assert_hc_flag_rejected(extra: &[&str], flag: &str) {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.vcf");
    let mut args = hc_fixture_args(&out);
    args.extend(extra.iter().map(|s| (*s).to_string()));
    let asserted = Command::cargo_bin("gatk-rs")
        .unwrap()
        .current_dir(repo_root())
        .args(&args)
        .assert()
        .code(2);
    let stderr = String::from_utf8(asserted.get_output().stderr.clone()).unwrap();
    assert!(
        stderr.contains(flag) && stderr.contains("not implemented"),
        "stderr for {flag}:\n{stderr}"
    );
    assert!(!out.exists(), "{flag} must not create the output VCF");
}

#[test]
fn haplotype_caller_rejects_non_default_stand_call_conf_without_creating_vcf() {
    assert_hc_flag_rejected(&["--stand-call-conf", "20"], "--stand-call-conf");
}

#[test]
fn haplotype_caller_rejects_non_default_min_base_quality_without_creating_vcf() {
    assert_hc_flag_rejected(
        &["--min-base-quality-score", "20"],
        "--min-base-quality-score",
    );
}

#[test]
fn haplotype_caller_rejects_original_base_qualities_without_creating_vcf() {
    assert_hc_flag_rejected(&["--original-base-qualities"], "--original-base-qualities");
}

#[test]
fn haplotype_caller_rejects_dont_use_soft_clipped_bases_without_creating_vcf() {
    assert_hc_flag_rejected(
        &["--dont-use-soft-clipped-bases"],
        "--dont-use-soft-clipped-bases",
    );
}

#[test]
fn haplotype_caller_explicit_production_defaults_still_call() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("hc.vcf");
    let mut args = hc_fixture_args(&out);
    args.extend([
        "--stand-call-conf".to_string(),
        "30".to_string(),
        "--min-base-quality-score".to_string(),
        "10".to_string(),
    ]);
    Command::cargo_bin("gatk-rs")
        .unwrap()
        .current_dir(repo_root())
        .env("GATK_RS_HC_SCAFFOLD_OUTPUT", "1")
        .args(&args)
        .assert()
        .success();
    let body = fs::read_to_string(&out).unwrap();
    assert!(body.contains("##GATK_RS_HC_PIPELINE=scaffold-v1"));
}

#[test]
fn haplotype_caller_rejects_non_default_max_alternate_alleles_without_creating_vcf() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.vcf");
    let mut args = hc_fixture_args(&out);
    args.extend(["--max-alternate-alleles".to_string(), "3".to_string()]);
    let asserted = Command::cargo_bin("gatk-rs")
        .unwrap()
        .current_dir(repo_root())
        .args(&args)
        .assert()
        .code(2);
    let stderr = String::from_utf8(asserted.get_output().stderr.clone()).unwrap();
    assert!(
        stderr.contains("--max-alternate-alleles"),
        "stderr:\n{stderr}"
    );
    assert!(!out.exists());
}
