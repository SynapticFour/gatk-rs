//! P12 parity regression gate (N5): shared / rust-only bounds after N-priority work.
//! Registry on (E.1): `P12_PHASE_E=1 GATK_RS_P12_EVENT_REGISTRY=1 P12_REFERENCE=… cargo test … p12_parity_gate --release -- --ignored --nocapture`
//! Graph-only: omit `GATK_RS_P12_EVENT_REGISTRY` (floor ~53 shared).

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions,
    read_event_discovery::p12_java_event_registry_enabled, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
};
use std::collections::BTreeSet;
use std::io::{BufRead, BufReader};
use std::path::Path;

const STAND_EMIT: f64 = 10.0;
/// M5 target when `GATK_RS_P12_EVENT_REGISTRY=1` (signed at 66/66 site trace).
const MIN_SHARED_REGISTRY_ON: usize = 66;
/// Graph-only ASM-8 floor (no list inject; target 66/66).
const MIN_SHARED_GRAPH_ONLY: usize = 66;
/// Rust-only position cap (M5 target 0 — graph-only emit matches Java baseline VCF).
const MAX_RUST_ONLY_DELTA: usize = 0;

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
    if !ref_path.is_file() {
        eprintln!("skip: P12_REFERENCE not found: {}", ref_path.display());
        eprintln!("  use the real FASTA path (not /path/to/... placeholder)");
        return None;
    }
    if !bam.is_file() {
        eprintln!("skip: P12 BAM not found: {}", bam.display());
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
#[ignore = "Phase E: long L3 gate; run with P12_PHASE_E=1 and --ignored"]
fn p12_parity_gate() {
    if std::env::var("P12_PHASE_E").is_err() {
        eprintln!("skip p12_parity_gate: export P12_PHASE_E=1 (any non-empty value)");
        return;
    }
    let Some((ref_fasta, bam)) = p12_paths() else {
        return;
    };
    let java_positions = load_java_positions();
    let interval = "2:92300000-92350000";
    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, interval).expect("interval");
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
    let registry_on = p12_java_event_registry_enabled();
    let min_shared = if registry_on {
        MIN_SHARED_REGISTRY_ON
    } else {
        MIN_SHARED_GRAPH_ONLY
    };
    eprintln!("=== P12 parity gate ===");
    eprintln!("shared_with_java\t{shared}");
    eprintln!("rust_emitted\t{}", emitted.len());
    eprintln!("rust_only_delta\t{rust_only}");
    eprintln!(
        "p12_event_registry\t{}",
        if registry_on { "on" } else { "off" }
    );
    eprintln!("min_shared\t{min_shared}");
    if shared < min_shared {
        let missing: Vec<_> = java_positions
            .iter()
            .filter(|pos| !emitted.contains(pos))
            .copied()
            .collect();
        eprintln!("missing_java_positions\t{}", missing.len());
        for pos in &missing {
            eprintln!("missing_java_pos\t{pos}");
        }
    }
    assert_eq!(
        shared,
        min_shared,
        "shared {shared} != {min_shared} (registry {})",
        if registry_on { "on" } else { "off" }
    );
    assert_eq!(
        emitted.len(),
        min_shared,
        "rust_emitted {} != {min_shared}",
        emitted.len()
    );
    assert_eq!(
        rust_only, MAX_RUST_ONLY_DELTA,
        "rust_only delta {rust_only}"
    );
}
