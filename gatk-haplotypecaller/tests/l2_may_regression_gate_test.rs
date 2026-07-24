//! Permanent L2 regression gates for the two cases that failed in the
//! 2026-05-17 `last_run.log` (`passed=143 failed=2`).
//! Historical failures (strict L2 vs frozen Java dumps):
//! `e0-assemble/p5_case1_assemble`: Rust emitted `status=just_assembled_reference`
//! while Java oracle is `status=failed` (short-ref multi-kmer gate dump).
//! `e2e/p5_indel_chrindel`: Rust emitted `status=assembled_some_variation`
//! `kmer_size=10` while Java materialize dump is `just_assembled_reference`
//! `kmer_size=0`.
//! These tests call the same dump entry points as `run_hc_full_parity_l2.sh` and
//! must remain even when the full L2 battery is green. Run with:
//! `cargo test -p gatk-haplotypecaller --features parity_harness --test l2_may_regression_gate_test`

use gatk_haplotypecaller::{
    dump_assembly_assemble_tsv, dump_assembly_region_haplotypes_tsv, AssemblyRegionHaplotypeTarget,
};
use std::collections::HashMap;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .canonicalize()
        .expect("repo root")
}

fn kv_meta(tsv: &str) -> HashMap<String, String> {
    let mut out = HashMap::new();
    for line in tsv.lines() {
        let mut parts = line.splitn(2, '\t');
        let Some(k) = parts.next() else { continue };
        let Some(v) = parts.next() else { continue };
        // Stop at haplotype table header (L2 compare also treats leading alpha keys as meta).
        if k == "rank" {
            break;
        }
        out.insert(k.to_string(), v.to_string());
    }
    out
}

/// Exact May L2 failure #1: e0-assemble / p5_case1_assemble.
#[test]
fn l2_e0_assemble_p5_case1_matches_java_failed_status() {
    let repo = repo_root();
    let ref_tsv = repo.join("parity/fixtures/hc-full-parity/e4/p5_case1_ref.tsv");
    let reads_tsv = repo.join("parity/fixtures/hc-full-parity/e4/p5_case1_reads.tsv");
    let java = repo.join(
        "parity/fixtures/hc-full-parity/java_dumps/e0-assemble/p5_case1_assemble_2dbc0258.tsv",
    );
    assert!(ref_tsv.is_file(), "missing {}", ref_tsv.display());
    assert!(reads_tsv.is_file(), "missing {}", reads_tsv.display());
    assert!(java.is_file(), "missing {}", java.display());

    let mut rust_buf = Vec::new();
    dump_assembly_assemble_tsv(&ref_tsv, &reads_tsv, &mut rust_buf).expect("rust dump");
    let rust_tsv = String::from_utf8(rust_buf).expect("utf8");
    let java_tsv = fs::read_to_string(&java).expect("java dump");

    let rust_kv = kv_meta(&rust_tsv);
    let java_kv = kv_meta(&java_tsv);
    assert_eq!(
        rust_kv.get("status"),
        java_kv.get("status"),
        "e0-assemble/p5_case1_assemble status drift\nrust:\n{rust_tsv}\njava:\n{java_tsv}"
    );
    assert_eq!(
        rust_kv.get("kmer_size"),
        java_kv.get("kmer_size"),
        "e0-assemble/p5_case1_assemble kmer_size drift\nrust:\n{rust_tsv}\njava:\n{java_tsv}"
    );
    assert_eq!(rust_kv.get("status").map(String::as_str), Some("failed"));
    assert_eq!(rust_kv.get("kmer_size").map(String::as_str), Some("25"));
}

/// Exact May L2 failure #2: e2e / p5_indel_chrindel.
/// Fixture bytes are embedded so the gate does not depend on gitignored `*.fa`/`*.sam`.
#[test]
fn l2_e2e_p5_indel_chrindel_matches_java_just_assembled_reference() {
    let repo = repo_root();
    let java =
        repo.join("parity/fixtures/hc-full-parity/java_dumps/e2e/p5_indel_chrindel_2dbc0258.tsv");
    assert!(java.is_file(), "missing {}", java.display());

    let dir = tempfile::tempdir().expect("tempdir");
    let fa = dir.path().join("p5_live_reference_indel.fa");
    let sam = dir.path().join("p5_live_case_indel.sam");
    let bam = dir.path().join("p5_live_case_indel.bam");

    fs::write(
        &fa,
        ">chrIndel\nTGCATGACTGATCGTACGATTCGAGCTAGTCGATCGATGCTAGCTAGGCTAACGTTAGCTAGTAACTG\n",
    )
    .expect("write fa");
    // Matching tracked `.fa.fai` layout (name, length, offset, linebases, linewidth).
    fs::write(
        dir.path().join("p5_live_reference_indel.fa.fai"),
        "chrIndel\t68\t10\t68\t69\n",
    )
    .expect("write fai");
    fs::write(
        dir.path().join("p5_live_reference_indel.dict"),
        "@HD\tVN:1.6\n@SQ\tSN:chrIndel\tLN:68\n",
    )
    .expect("write dict");

    // Exact local fixture content from May/July L2 runs (gitignored `*.sam` on disk).
    let mut sam_file = fs::File::create(&sam).expect("sam");
    writeln!(sam_file, "@HD\tVN:1.6\tSO:coordinate").unwrap();
    writeln!(sam_file, "@SQ\tSN:chrIndel\tLN:68").unwrap();
    writeln!(sam_file, "@RG\tID:rg1\tSM:s1\tPL:ILLUMINA").unwrap();
    writeln!(
        sam_file,
        "ri1\t0\tchrIndel\t5\t60\t20M\t*\t0\t0\tTGACTGATCGTACGATTCGA\tFFFFFFFFFFFFFFFFFFFF\tRG:Z:rg1"
    )
    .unwrap();
    writeln!(
        sam_file,
        "ri2\t0\tchrIndel\t5\t60\t10M1I9M\t*\t0\t0\tTGACTGATCGATACGATTCG\tFFFFFFFFFFFFFFFFFFFF\tRG:Z:rg1"
    )
    .unwrap();
    writeln!(
        sam_file,
        "ri3\t0\tchrIndel\t5\t60\t10M1I9M\t*\t0\t0\tTGACTGATCGATACGATTCG\tFFFFFFFFFFFFFFFFFFFF\tRG:Z:rg1"
    )
    .unwrap();
    writeln!(
        sam_file,
        "ri4\t0\tchrIndel\t5\t60\t10M1D10M\t*\t0\t0\tTGACTGATCGCGATTCGAGC\tFFFFFFFFFFFFFFFFFFFF\tRG:Z:rg1"
    )
    .unwrap();
    writeln!(
        sam_file,
        "ri5\t0\tchrIndel\t5\t60\t20M\t*\t0\t0\tTGACTGATCGTACGATTCGA\tFFFFFFFFFFFFFFFFFFFF\tRG:Z:rg1"
    )
    .unwrap();
    drop(sam_file);

    let status = Command::new("sh")
        .arg("-c")
        .arg(format!(
            "samtools view -bS {} | samtools sort -o {}",
            sam.display(),
            bam.display()
        ))
        .status()
        .expect("samtools");
    assert!(
        status.success(),
        "sam->bam failed for embedded p5_indel fixture"
    );
    let idx = Command::new("samtools")
        .args(["index", bam.to_str().unwrap()])
        .status()
        .expect("samtools index");
    assert!(idx.success(), "samtools index failed");

    let mut rust_buf = Vec::new();
    dump_assembly_region_haplotypes_tsv(
        &fa,
        &bam,
        "chrIndel:1-40",
        0,
        AssemblyRegionHaplotypeTarget::Active,
        &mut rust_buf,
    )
    .expect("rust e2e dump");
    let rust_tsv = String::from_utf8(rust_buf).expect("utf8");
    let java_tsv = fs::read_to_string(&java).expect("java dump");

    let rust_kv = kv_meta(&rust_tsv);
    let java_kv = kv_meta(&java_tsv);
    assert_eq!(
        rust_kv.get("status"),
        java_kv.get("status"),
        "e2e/p5_indel_chrindel status drift\nrust:\n{rust_tsv}\njava:\n{java_tsv}"
    );
    assert_eq!(
        rust_kv.get("kmer_size"),
        java_kv.get("kmer_size"),
        "e2e/p5_indel_chrindel kmer_size drift\nrust:\n{rust_tsv}\njava:\n{java_tsv}"
    );
    assert_eq!(
        rust_kv.get("status").map(String::as_str),
        Some("just_assembled_reference")
    );
    assert_eq!(rust_kv.get("kmer_size").map(String::as_str), Some("0"));
}
