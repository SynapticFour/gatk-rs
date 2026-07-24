//! ASM-8-only production experiment: CIGAR/EventMap path without P12 `ensure_*` bridges.
//! ```bash
//! export P12_PHASE_E=1 GATK_RS_ASM8_ONLY=1 P12_REFERENCE=parity/realworld/assets/hs37d5.simple.fa
//! cargo test -p gatk-haplotypecaller --test p12_asm8_only_gate_test --release -- --ignored --nocapture
//! ```
//! Target (production): `shared=66`, `rust_only=0` **without** `P12_PHASE_E` whitelist.
//! Harness (`P12_PHASE_E=1`): records regression vs bridge-on path until CIGAR parity is complete.

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions,
    read_event_discovery::{strict_java_asm8_only_enabled, strict_java_p12_ensure_bridges_enabled},
    traverse_assembly_region_walker, try_emit_call_region_variants, AssemblyRegionCallDisposition,
    CallRegionArgs, HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
};
use std::collections::BTreeSet;
use std::io::{BufRead, BufReader};
use std::path::Path;

const STAND_EMIT: f64 = 10.0;
/// Target 66/66 once ASM-8 CIGAR finalize matches bridge-on path.
const MIN_SHARED_ASM8_HARNESS: usize = 66;
const MAX_RUST_ONLY_HARNESS: usize = 0;

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
        eprintln!("skip: P12_REFERENCE / BAM missing");
        return None;
    }
    Some((ref_path, bam))
}

fn load_java_positions() -> BTreeSet<u64> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../parity/fixtures/p12-java-production-emit/p12_production_emit_sites.tsv");
    let mut out = BTreeSet::new();
    if let Ok(file) = std::fs::File::open(path) {
        for line in BufReader::new(file).lines().skip(1).flatten() {
            let cols: Vec<_> = line.split('\t').collect();
            if cols.len() >= 2 {
                if let Ok(pos) = cols[1].parse::<u64>() {
                    out.insert(pos);
                }
            }
        }
    }
    out
}

#[test]
fn p12_asm8_only_bridge_flags() {
    let asm8_env = std::env::var("GATK_RS_ASM8_ONLY").ok();
    let bridges_env = std::env::var("GATK_RS_P12_ENSURE_BRIDGES").ok();
    if asm8_env
        .as_deref()
        .is_some_and(|v| matches!(v, "1" | "true" | "TRUE" | "yes" | "YES"))
    {
        assert!(strict_java_asm8_only_enabled());
        assert!(!strict_java_p12_ensure_bridges_enabled());
    } else if bridges_env
        .as_deref()
        .is_some_and(|v| matches!(v, "1" | "true" | "TRUE" | "yes" | "YES"))
    {
        assert!(!strict_java_asm8_only_enabled());
        assert!(strict_java_p12_ensure_bridges_enabled());
    } else {
        assert!(
            strict_java_asm8_only_enabled(),
            "default: graph-only production"
        );
        assert!(
            !strict_java_p12_ensure_bridges_enabled(),
            "default: bridges off"
        );
    }
}

#[test]
#[ignore = "ASM-8 harness: long; GATK_RS_ASM8_ONLY=1; P12_PHASE_E optional (emit whitelist if set)"]
fn p12_asm8_only_gate() {
    if !strict_java_asm8_only_enabled() {
        eprintln!("skip: export GATK_RS_ASM8_ONLY=1");
        return;
    }
    let Some((ref_fasta, bam)) = p12_paths() else {
        return;
    };
    let java_positions = load_java_positions();
    assert_eq!(java_positions.len(), 66);
    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, "2:92300000-92350000").expect("interval");
    let walk = traverse_assembly_region_walker(
        &dict,
        &specs,
        &ref_fasta,
        &bam,
        &ReadFilterParams::gatk_standard_hc(),
        &WalkerTraversalConfig::gatk_haplotype_caller_production(100),
    )
    .expect("walk");
    let regions = flatten_assembly_regions(&walk);
    let args = CallRegionArgs::strict_java();

    let mut emitted = BTreeSet::new();
    let mut shared = 0usize;
    for region in &regions {
        if !matches!(
            call_disposition(region),
            AssemblyRegionCallDisposition::ActiveFull
        ) {
            continue;
        }
        let Some(outcome) =
            HaplotypeCallerEngine::call_region(region, &dict, &ref_fasta, &args).expect("call")
        else {
            continue;
        };
        for rec in
            try_emit_call_region_variants(region, &outcome, "SAMPLE", STAND_EMIT).expect("emit")
        {
            if java_positions.contains(&rec.position) {
                shared += 1;
            }
            emitted.insert(rec.position);
        }
    }
    let rust_only = emitted.len().saturating_sub(shared);
    eprintln!("=== P12 ASM-8-only gate (harness) ===");
    eprintln!(
        "p12_ensure_bridges\t{}",
        strict_java_p12_ensure_bridges_enabled()
    );
    eprintln!("shared_with_java\t{shared}");
    eprintln!("rust_emitted\t{}", emitted.len());
    eprintln!("rust_only_delta\t{rust_only}");
    eprintln!("min_shared_floor\t{MIN_SHARED_ASM8_HARNESS}");
    assert!(
        shared >= MIN_SHARED_ASM8_HARNESS,
        "ASM-8-only shared {shared} < floor {MIN_SHARED_ASM8_HARNESS}"
    );
    assert!(
        rust_only <= MAX_RUST_ONLY_HARNESS,
        "rust_only {rust_only} > {MAX_RUST_ONLY_HARNESS} (harness whitelist)"
    );
}
