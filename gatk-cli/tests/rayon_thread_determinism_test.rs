//! Required CI gate: HaplotypeCaller VCF must be byte-identical across Rayon thread counts.
//! Rayon's global pool is initialized once per process from `RAYON_NUM_THREADS` / `--threads`,
//! so each trial runs as a **subprocess** of `gatk-rs` with a fresh env.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use tempfile::tempdir;

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../parity/fixtures")
        .join(name)
}

fn normalize_vcf_for_compare(bytes: &[u8]) -> Vec<u8> {
    let text = String::from_utf8_lossy(bytes);
    let mut out = String::new();
    for line in text.lines() {
        if line.starts_with("##fileDate=")
            || line.starts_with("##source=")
            || line.starts_with("##GATKCommandLine")
        {
            continue;
        }
        out.push_str(line);
        out.push('\n');
    }
    out.into_bytes()
}

fn run_hc(args: &[&str], env_threads: Option<u32>, out_vcf: &Path) -> Vec<u8> {
    let bin = env!("CARGO_BIN_EXE_gatk-rs");
    let ref_fa = fixture("reference.fa");
    let bam = fixture("sample.bam");
    let mut cmd = Command::new(bin);
    cmd.env_remove("JAVA_TOOL_OPTIONS");
    if let Some(t) = env_threads {
        cmd.env("RAYON_NUM_THREADS", t.to_string());
    } else {
        cmd.env_remove("RAYON_NUM_THREADS");
    }
    let status = cmd
        .args(["HaplotypeCaller", "-R"])
        .arg(&ref_fa)
        .arg("-I")
        .arg(&bam)
        .arg("-O")
        .arg(out_vcf)
        .arg("-L")
        .arg("chr1:1-32")
        .args(args)
        .status()
        .unwrap_or_else(|e| panic!("spawn gatk-rs: {e}"));
    assert!(
        status.success(),
        "HaplotypeCaller failed: {status:?} args={args:?}"
    );
    fs::read(out_vcf).unwrap_or_else(|e| panic!("read {}: {e}", out_vcf.display()))
}

#[test]
fn haplotype_caller_vcf_byte_identical_across_rayon_thread_counts() {
    let dir = tempdir().expect("tempdir");
    let thread_counts = [1u32, 2, 4, 8, 16];
    let mut normalized: Vec<(u32, Vec<u8>)> = Vec::with_capacity(thread_counts.len());

    for (i, &threads) in thread_counts.iter().enumerate() {
        let out = dir.path().join(format!("hc_env_t{threads}_{i}.vcf"));
        let raw = run_hc(&[], Some(threads), &out);
        assert!(!raw.is_empty(), "empty VCF for RAYON_NUM_THREADS={threads}");
        let norm = normalize_vcf_for_compare(&raw);
        assert!(
            String::from_utf8_lossy(&norm).contains("#CHROM"),
            "missing VCF header for threads={threads}"
        );
        normalized.push((threads, norm));
    }

    let (t0, ref_bytes) = &normalized[0];
    for (t, bytes) in &normalized[1..] {
        assert_eq!(
            bytes, ref_bytes,
            "VCF not byte-identical: RAYON_NUM_THREADS={t} vs {t0} (after normalizing volatile headers)"
        );
    }
}

#[test]
fn haplotype_caller_vcf_byte_identical_across_threads_cli_flag() {
    let dir = tempdir().expect("tempdir");
    let thread_counts = [1u32, 2, 4, 8, 16];
    let mut normalized: Vec<(u32, Vec<u8>)> = Vec::with_capacity(thread_counts.len());

    for (i, &threads) in thread_counts.iter().enumerate() {
        let out = dir.path().join(format!("hc_cli_t{threads}_{i}.vcf"));
        let t = threads.to_string();
        let raw = run_hc(&["--threads", &t], None, &out);
        assert!(!raw.is_empty(), "empty VCF for --threads={threads}");
        let norm = normalize_vcf_for_compare(&raw);
        assert!(
            String::from_utf8_lossy(&norm).contains("#CHROM"),
            "missing VCF header for --threads={threads}"
        );
        normalized.push((threads, norm));
    }

    let (t0, ref_bytes) = &normalized[0];
    for (t, bytes) in &normalized[1..] {
        assert_eq!(
            bytes, ref_bytes,
            "VCF not byte-identical: --threads={t} vs {t0}"
        );
    }

    // CLI --threads and env RAYON_NUM_THREADS must agree on the same fixture.
    let env_out = dir.path().join("hc_env_ref.vcf");
    let env_raw = run_hc(&[], Some(4), &env_out);
    assert_eq!(
        normalize_vcf_for_compare(&env_raw),
        *ref_bytes,
        "--threads path must match RAYON_NUM_THREADS path"
    );
}
