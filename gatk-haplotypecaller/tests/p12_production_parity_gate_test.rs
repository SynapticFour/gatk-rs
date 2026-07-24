//! L3b production parity: strict_java graph-only emit vs Java **without** `P12_PHASE_E`.
//! Default: bridges off (ASM-8 graph-only). Legacy bridges: `GATK_RS_P12_ENSURE_BRIDGES=1`.
//! ```bash
//! env -u P12_PHASE_E P12_REFERENCE=parity/realworld/assets/hs37d5.simple.fa \
//! cargo test -p gatk-haplotypecaller --test p12_production_parity_gate_test \
//! p12_production_parity_gate --release -- --ignored --nocapture
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions,
    read_event_discovery::{
        p12_emit_baseline_filter_enabled, p12_java_event_registry_enabled,
        strict_java_asm8_only_enabled, strict_java_p12_ensure_bridges_enabled,
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

fn load_java_variant_keys() -> BTreeSet<(u64, String, String)> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../parity/fixtures/p12-java-production-emit/p12_production_emit_sites.tsv");
    let mut out = BTreeSet::new();
    if let Ok(file) = std::fs::File::open(path) {
        for line in BufReader::new(file).lines().skip(1).flatten() {
            let cols: Vec<_> = line.split('\t').collect();
            if cols.len() >= 4 {
                if let Ok(pos) = cols[1].parse::<u64>() {
                    out.insert((pos, cols[2].into(), cols[3].into()));
                }
            }
        }
    }
    out
}

#[test]
#[ignore = "L3b production: long gate (~6 min); do not set P12_PHASE_E"]
fn p12_production_parity_gate() {
    if std::env::var("P12_PHASE_E").is_ok() {
        eprintln!("skip production gate: P12_PHASE_E is set (use L3a p12_parity_gate instead)");
        return;
    }
    let Some((ref_fasta, bam)) = p12_paths() else {
        eprintln!("skip: set P12_REFERENCE");
        return;
    };
    assert!(
        !p12_emit_baseline_filter_enabled(),
        "production gate must not use P12_PHASE_E emit whitelist"
    );
    if std::env::var("GATK_RS_P12_ENSURE_BRIDGES").is_err() {
        assert!(
            !strict_java_p12_ensure_bridges_enabled(),
            "L3b sign-off: unset GATK_RS_P12_ENSURE_BRIDGES (graph-only default)"
        );
        assert!(
            strict_java_asm8_only_enabled(),
            "L3b sign-off: graph-only production must be active"
        );
    }

    let java_keys = load_java_variant_keys();
    assert_eq!(
        java_keys.len(),
        MIN_SHARED,
        "p12_java_only.tsv must list {MIN_SHARED} variants"
    );
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
            let alt = rec.alternate.first().cloned().unwrap_or_default();
            let key = (rec.position, rec.reference.clone(), alt);
            emitted.insert(key);
        }
    }
    let shared = emitted.iter().filter(|k| java_keys.contains(*k)).count();
    let rust_only = emitted.iter().filter(|k| !java_keys.contains(*k)).count();
    eprintln!("=== P12 production parity gate (L3b) ===");
    eprintln!("shared_with_java\t{shared}");
    eprintln!("rust_emitted\t{}", emitted.len());
    eprintln!("rust_only_delta\t{rust_only}");
    eprintln!(
        "p12_emit_whitelist\t{}",
        if p12_emit_baseline_filter_enabled() {
            "on"
        } else {
            "off"
        }
    );
    eprintln!(
        "p12_event_registry\t{}",
        if p12_java_event_registry_enabled() {
            "on"
        } else {
            "off"
        }
    );
    eprintln!(
        "p12_ensure_bridges\t{}",
        if strict_java_p12_ensure_bridges_enabled() {
            "on"
        } else {
            "off"
        }
    );
    eprintln!(
        "asm8_only\t{}",
        if strict_java_asm8_only_enabled() {
            "on"
        } else {
            "off"
        }
    );
    if shared < MIN_SHARED {
        let missing: Vec<_> = java_keys
            .iter()
            .filter(|k| !emitted.contains(*k))
            .map(|(pos, _, _)| *pos)
            .collect();
        eprintln!("missing_java_variants\t{}", missing.len());
        for pos in &missing {
            eprintln!("missing_java_pos\t{pos}");
        }
    }
    if rust_only > MAX_RUST_ONLY {
        for (pos, ref_a, alt_a) in emitted.iter().filter(|k| !java_keys.contains(*k)) {
            eprintln!("rust_only_variant\t{pos}\t{ref_a}\t{alt_a}");
        }
    }
    assert_eq!(shared, MIN_SHARED, "shared {shared} != {MIN_SHARED}");
    assert_eq!(
        emitted.len(),
        MIN_SHARED,
        "rust_emitted {} != {MIN_SHARED}",
        emitted.len()
    );
    assert_eq!(rust_only, MAX_RUST_ONLY, "rust_only delta {rust_only}");
}
