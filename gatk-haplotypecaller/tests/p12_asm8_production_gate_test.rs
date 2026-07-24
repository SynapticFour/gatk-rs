//! ASM-8 production gate: CIGAR/EventMap finalize without `P12_PHASE_E` or post-HMM `ensure_*` bridges.
//! ```bash
//! env -u P12_PHASE_E GATK_RS_ASM8_ONLY=1 P12_REFERENCE="$PWD/parity/realworld/assets/hs37d5.simple.fa" \
//! env -u P12_PHASE_E GATK_RS_ASM8_ONLY=1 \
//! cargo test -p gatk-haplotypecaller --test p12_asm8_production_gate_test \
//! p12_asm8_production_parity_gate --release -- --ignored --nocapture
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions,
    read_event_discovery::{
        p12_emit_baseline_filter_enabled, strict_java_asm8_only_enabled,
        strict_java_p12_ensure_bridges_enabled,
    },
    traverse_assembly_region_walker, try_emit_call_region_variants, AssemblyRegionCallDisposition,
    CallRegionArgs, HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
};
use std::collections::BTreeSet;
use std::io::{BufRead, BufReader};
use std::path::Path;

const STAND_EMIT: f64 = 10.0;
const MIN_SHARED: usize = 66;
const MAX_RUST_ONLY: usize = 0;

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
    if ref_path.is_file() && bam.is_file() {
        Some((ref_path, bam))
    } else {
        None
    }
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
#[ignore = "ASM-8 production: long; GATK_RS_ASM8_ONLY=1; unset P12_PHASE_E"]
fn p12_asm8_production_parity_gate() {
    if std::env::var("P12_PHASE_E").is_ok() {
        eprintln!("skip: unset P12_PHASE_E");
        return;
    }
    if !strict_java_asm8_only_enabled() {
        eprintln!("skip: graph-only production off (set GATK_RS_ASM8_ONLY=1 or unset GATK_RS_P12_ENSURE_BRIDGES)");
        return;
    }
    assert!(!strict_java_p12_ensure_bridges_enabled());
    assert!(!p12_emit_baseline_filter_enabled());

    let Some((ref_fasta, bam)) = p12_paths() else {
        eprintln!("skip: set P12_REFERENCE");
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
    eprintln!("=== P12 ASM-8 production gate ===");
    eprintln!("shared_with_java\t{shared}");
    eprintln!("rust_emitted\t{}", emitted.len());
    eprintln!("rust_only_delta\t{rust_only}");
    eprintln!(
        "p12_ensure_bridges\t{}",
        if strict_java_p12_ensure_bridges_enabled() {
            "on"
        } else {
            "off"
        }
    );
    if shared < MIN_SHARED {
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
    assert_eq!(shared, MIN_SHARED, "ASM-8 shared {shared} != {MIN_SHARED}");
    assert_eq!(
        emitted.len(),
        MIN_SHARED,
        "ASM-8 rust_emitted {} != {MIN_SHARED}",
        emitted.len()
    );
    assert_eq!(rust_only, MAX_RUST_ONLY, "ASM-8 rust_only {rust_only}");
}
