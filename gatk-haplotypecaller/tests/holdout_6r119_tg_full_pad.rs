//! 6R.119 live dump: post-trim full-padded-ref TG haplotype must be absent.
//! Skipped unless `HOLDOUT_6R119=1`. Does not reopen 6R.113–6R.118.
//!
//! ```text
//! HOLDOUT_6R119=1 P12_REFERENCE=$PWD/parity/realworld/assets/hs37d5.simple.fa \
//!   cargo test -p gatk-haplotypecaller --test holdout_6r119_tg_full_pad -- --nocapture
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::{collect_variation_events, VariationEvent};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    AssemblyRegionCallDisposition, CallRegionArgs, HaplotypeCallerEngine, ReadFilterParams,
    WalkerTraversalConfig,
};
use serde_json::json;
use std::path::{Path, PathBuf};

const INTERVAL: &str = "2:92300000-92350000";
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 92_307_327;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn fnv1a64_hex(bases: &[u8]) -> String {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in bases {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    format!("{h:016x}")
}

fn events_at(events: &[VariationEvent], loc: u64) -> Vec<String> {
    events
        .iter()
        .filter(|e| e.start_1based.get() == loc)
        .map(|e| format!("{}/{}", e.ref_allele, e.alt_allele))
        .collect()
}

#[test]
fn holdout_6r119_live_no_full_pad_tg_supplement() {
    if std::env::var("HOLDOUT_6R119").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R119=1");
        return;
    }
    let root = repo_root();
    let ref_fasta = std::env::var("P12_REFERENCE")
        .map(PathBuf::from)
        .unwrap_or_else(|_| root.join(REF_REL));
    let bam = root.join(BAM_REL);
    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, INTERVAL).expect("interval");
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
    let covering = regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= TARGET
                && r.end.get() >= TARGET
        })
        .expect("ActiveFull covering target");
    let args = CallRegionArgs::strict_java();
    let outcome = HaplotypeCallerEngine::call_region(covering, &dict, &ref_fasta, &args)
        .expect("call")
        .expect("outcome");
    let full_pad = outcome.assembly.padded_reference_start_1based();
    let full_ref = outcome.assembly.reference_bases();
    let haps = &outcome.assembly.haplotypes;
    let max_mnp = outcome.assembly.max_mnp_distance();
    let union = collect_variation_events(haps, full_ref, full_pad, "2", max_mnp);
    let n_full_len_supplement = haps
        .iter()
        .filter(|h| {
            !h.is_reference && h.bases.len() == full_ref.len() && (h.score - 1000.0).abs() < 1e-6
        })
        .count();
    let per: Vec<_> = haps
        .iter()
        .enumerate()
        .map(|(i, h)| {
            json!({
                "idx": i,
                "hash": fnv1a64_hex(&h.bases),
                "len": h.bases.len(),
                "is_reference": h.is_reference,
                "score": h.score,
                "cigar": h.cigar.as_ref().map(|c| c.to_gatk_string()).unwrap_or_else(|| ".".into()),
            })
        })
        .collect();
    let dump = json!({
        "mission": "6R.119",
        "hap_count": haps.len(),
        "full_ref_len": full_ref.len(),
        "n_full_len_supplement_1000": n_full_len_supplement,
        "eventmap_union_at_target": events_at(&union, TARGET),
        "per_haplotype": per,
    });
    println!("{}", serde_json::to_string_pretty(&dump).expect("json"));
    assert_eq!(
        n_full_len_supplement, 0,
        "post-trim full-padded-ref score-1000 SNP haplotype must be absent"
    );
}
