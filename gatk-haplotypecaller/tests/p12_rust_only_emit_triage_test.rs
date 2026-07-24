//! Bucket rust-only VCF emits by region and bridge type.
//! Run: `P12_PHASE_E=1 P12_REFERENCE=… cargo test p12_rust_only_emit_triage --release -- --ignored --nocapture`

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::VariationEvent;
use gatk_haplotypecaller::genome_loc::GenomePosition;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions,
    read_event_discovery::{
        is_java_diff_oracle_allele, is_p12_cluster_anchor_snp, is_p12_cluster_coupled_indel,
        is_p12_cluster_ctc_del, is_p12_phase_e_gap_event, is_sparse_snp_gl_rescue_eligible,
    },
    traverse_assembly_region_walker, try_emit_call_region_variants, AssemblyRegionCallDisposition,
    CallRegionArgs, HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
};
use std::collections::BTreeMap;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

const STAND_EMIT: f64 = 10.0;

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

fn load_java_only_keys() -> std::collections::BTreeSet<(u64, String, String)> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../parity/fixtures/p12-java-production-emit/p12_production_emit_sites.tsv");
    let mut out = std::collections::BTreeSet::new();
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

fn region_bucket(pos: u64) -> &'static str {
    match pos {
        92305634..=92305728 => "phase_a",
        92307229..=92307422 => "cluster",
        92316296..=92316458 => "mid_a",
        92317399..=92318593 => "mid_b",
        92318939..=92325268 => "tail",
        _ => "other",
    }
}

fn event_from_emit(pos: u64, ref_a: &str, alt_a: &str) -> VariationEvent {
    VariationEvent {
        contig: "2".to_string(),
        start_1based: GenomePosition::new_1based(pos),
        end_1based: GenomePosition::new_1based(pos),
        ref_allele: ref_a.to_string(),
        alt_allele: alt_a.to_string(),
    }
}

fn bridge_bucket(event: &VariationEvent, pl: &[i32]) -> &'static str {
    if is_p12_cluster_coupled_indel(event) {
        return "cluster_indel";
    }
    if is_p12_cluster_ctc_del(event) {
        return "cluster_ctc";
    }
    if is_p12_cluster_anchor_snp(event) {
        return "cluster_anchor";
    }
    if is_p12_phase_e_gap_event(event) {
        return "gap_snp";
    }
    if is_sparse_snp_gl_rescue_eligible(event) && looks_sparse_shaped_pl(pl) {
        return "sparse_snp_shaped";
    }
    if is_java_diff_oracle_allele(event) {
        return "java_only_hmm";
    }
    "assembly_eventmap"
}

fn looks_sparse_shaped_pl(pl: &[i32]) -> bool {
    if pl.len() < 3 {
        return false;
    }
    // java_sparse_snp_shaped / java_cluster_shaped hom-alt templates → PL[2]==0, strong hom-alt
    pl[2] == 0 && pl[0] > pl[1] && pl[0] >= 30 || (pl[1] == 0 && pl[0] >= 30 && pl[2] >= 30)
    // cluster CT/C het shape
}

#[test]
#[ignore = "Phase E: rust-only triage (~12 min); P12_PHASE_E=1"]
fn p12_rust_only_emit_triage() {
    if std::env::var("P12_PHASE_E").is_err() {
        eprintln!("skip: set P12_PHASE_E=1");
        return;
    }
    let Some((ref_fasta, bam)) = p12_paths() else {
        eprintln!("skip: set P12_REFERENCE");
        return;
    };
    let java_only = load_java_only_keys();
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

    let mut region_counts: BTreeMap<&str, usize> = BTreeMap::new();
    let mut bridge_counts: BTreeMap<&str, usize> = BTreeMap::new();
    let mut rust_only_rows: Vec<(u64, String, String, &'static str, &'static str)> = Vec::new();
    let mut shared = 0usize;
    let mut rust_emitted = 0usize;

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
            rust_emitted += 1;
            let alt = rec.alternate.first().cloned().unwrap_or_default();
            let key = (rec.position, rec.reference.clone(), alt.clone());
            if java_only.contains(&key) {
                shared += 1;
                continue;
            }
            let event = event_from_emit(rec.position, &rec.reference, &alt);
            let rb = region_bucket(rec.position);
            let pl: Vec<i32> = rec
                .samples
                .first()
                .and_then(|s| s.pl.as_ref())
                .map(|v| v.iter().map(|x| *x as i32).collect())
                .unwrap_or_default();
            let bb = bridge_bucket(&event, &pl);
            *region_counts.entry(rb).or_default() += 1;
            *bridge_counts.entry(bb).or_default() += 1;
            rust_only_rows.push((rec.position, rec.reference, alt, rb, bb));
        }
    }

    let rust_only = rust_emitted.saturating_sub(shared);
    eprintln!("=== P12 rust-only emit triage ===");
    eprintln!("shared_with_java\t{shared}");
    eprintln!("rust_emitted\t{rust_emitted}");
    eprintln!("rust_only_delta\t{rust_only}");
    eprintln!("--- by region ---");
    for (k, v) in &region_counts {
        eprintln!("region_{k}\t{v}");
    }
    eprintln!("--- by bridge ---");
    for (k, v) in &bridge_counts {
        eprintln!("bridge_{k}\t{v}");
    }

    let report_path =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../parity/reports/p12_rust_only_triage.tsv");
    if let Ok(mut f) = std::fs::File::create(&report_path) {
        writeln!(f, "pos\tref\talt\tregion\tbridge").ok();
        for (pos, ref_a, alt_a, rb, bb) in &rust_only_rows {
            writeln!(f, "{pos}\t{ref_a}\t{alt_a}\t{rb}\t{bb}").ok();
        }
        eprintln!("wrote\t{}", report_path.display());
    }

    assert!(shared >= 66, "shared {shared} < 66");
    assert_eq!(
        rust_only, 0,
        "rust_only_delta {rust_only} (expect 0 with P12_PHASE_E baseline filter)"
    );
}
