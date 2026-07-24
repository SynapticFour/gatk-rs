//! J1: P12 sparse desert assembly-region tiling matches Java (activity profile parity).
//! Run: `P12_REFERENCE=…/hs37d5.simple.fa cargo test -p gatk-haplotypecaller p12_desert_assembly_regions_java_parity --release`

use std::path::Path;
use std::process::Command;

fn p12_paths() -> Option<(std::path::PathBuf, std::path::PathBuf)> {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
    let ref_path = std::env::var("P12_REFERENCE")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|_| root.join("parity/realworld/assets/hs37d5.simple.fa"));
    let ref_path = if ref_path.is_absolute() {
        ref_path
    } else {
        root.join(ref_path)
    };
    let bam = root.join("parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam");
    if !ref_path.is_file() || !bam.is_file() {
        return None;
    }
    Some((ref_path, bam))
}

fn dump_assembly_regions(
    bin: &Path,
    ref_fasta: &Path,
    bam: &Path,
    interval: &str,
) -> Option<String> {
    let out = Command::new(bin)
        .arg("assembly-regions")
        .arg(ref_fasta)
        .arg(bam)
        .arg(interval)
        .output()
        .ok()?;
    if !out.status.success() {
        eprintln!("dump failed: {}", String::from_utf8_lossy(&out.stderr));
        return None;
    }
    Some(String::from_utf8_lossy(&out.stdout).into_owned())
}

#[test]
fn p12_desert_assembly_regions_java_parity() {
    let Some((ref_fasta, bam)) = p12_paths() else {
        eprintln!("skip: set P12_REFERENCE");
        return;
    };
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
    let target = std::env::var("CARGO_TARGET_DIR")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|_| root.join("target"));
    let rust_bin = target.join("release/examples/hc_full_parity_gate_dump");
    if !rust_bin.is_file() {
        eprintln!("skip: build hc_full_parity_gate_dump release example first");
        return;
    }
    let java_dump = root.join("scripts/parity/run_hc_full_parity_java_dump.sh");
    if !java_dump.is_file() {
        eprintln!("skip: java dump script missing");
        return;
    }

    let interval = "2:92305800-92307000";
    let rust_tsv = dump_assembly_regions(&rust_bin, &ref_fasta, &bam, interval)
        .expect("rust assembly-regions dump");
    let java_out = Command::new(&java_dump)
        .arg("assembly-regions")
        .arg(&ref_fasta)
        .arg(&bam)
        .arg(interval)
        .output()
        .expect("java dump");
    assert!(
        java_out.status.success(),
        "java dump: {}",
        String::from_utf8_lossy(&java_out.stderr)
    );
    let java_tsv = String::from_utf8_lossy(&java_out.stdout);

    assert_eq!(
        rust_tsv.trim(),
        java_tsv.trim(),
        "desert assembly-region tiling must match Java"
    );
}
