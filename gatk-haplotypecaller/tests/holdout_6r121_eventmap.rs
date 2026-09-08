//! 6R.121 forensic dump: EventMap at HOLDOUT_6R53 remainder `20:29456196 A/T`.
//! Inventory A only. Skipped unless `HOLDOUT_6R121=1`. Production change: NONE.
//!
//! ```text
//! HOLDOUT_6R121=1 cargo test -p gatk-haplotypecaller --test holdout_6r121_eventmap -- --nocapture
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::{
    build_per_haplotype_variation_events, merged_alleles_for_genotyping, overlapping_events,
    variation_events_at_position_from_cache, variation_events_for_haplotype, EventMap,
    VariationEvent,
};
use gatk_haplotypecaller::hc_allele_mapping::{hap_base_at_ref_locus, replace_span_del_events};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    AssemblyRegionCallDisposition, CallRegionArgs, Haplotype, HaplotypeCallerEngine,
    ReadFilterParams, WalkerTraversalConfig,
};
use serde_json::{json, Value};
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const JAVA_VCF_REL: &str = "parity/reports/6r43/chr20_tiny/java.vcf";
const RUST_VCF_REL: &str = "parity/reports/6r43/chr20_tiny/rust.vcf";
const TARGET: u64 = 29_456_196;

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

fn cigar_str(h: &Haplotype) -> String {
    h.cigar
        .as_ref()
        .map(|c| c.to_gatk_string())
        .unwrap_or_else(|| ".".to_string())
}

fn cigar_has_indel(h: &Haplotype) -> bool {
    h.cigar
        .as_ref()
        .is_some_and(|c| c.elements.iter().any(|e| e.operator.is_indel()))
}

fn event_json(e: &VariationEvent) -> Value {
    json!({
        "start": e.start_1based.get(),
        "end": e.end_1based.get(),
        "ref": e.ref_allele,
        "alt": e.alt_allele,
        "indel": e.is_indel(),
        "kind": if e.is_indel() { "indel" } else { "snp" },
    })
}

fn is_at(e: &VariationEvent) -> bool {
    e.start_1based.get() == TARGET && e.ref_allele == "A" && e.alt_allele == "T"
}

fn vcf_record(path: &Path, pos: u64) -> Option<Value> {
    for line in std::fs::read_to_string(path).unwrap_or_default().lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<_> = line.split('\t').collect();
        if f.len() < 10 || f[0] != "20" {
            continue;
        }
        if f[1].parse::<u64>().ok()? != pos {
            continue;
        }
        let fmt: Vec<_> = f[8].split(':').collect();
        let samp: Vec<_> = f[9].split(':').collect();
        let get = |k: &str| {
            fmt.iter()
                .position(|x| *x == k)
                .and_then(|i| samp.get(i))
                .unwrap_or(&".")
                .to_string()
        };
        return Some(json!({
            "pos": pos,
            "alleles": format!("{}/{}", f[3], f[4]),
            "qual": f[5],
            "gt": get("GT"),
            "ad": get("AD"),
            "pl": get("PL"),
        }));
    }
    None
}

fn events_at(events: &[VariationEvent], loc: u64) -> Vec<Value> {
    events
        .iter()
        .filter(|e| {
            e.start_1based.get() == loc
                || (e.start_1based.get() <= loc && e.end_1based.get() >= loc)
        })
        .map(event_json)
        .collect()
}

fn cigar_only_events(
    h: &Haplotype,
    ref_hap: &Haplotype,
    ref_bytes: &[u8],
    pad: u64,
) -> Vec<VariationEvent> {
    if h.cigar.is_none() {
        return Vec::new();
    }
    EventMap::from_haplotype_and_reference(h, ref_hap, ref_bytes, pad, 0)
        .variation_events("20", pad)
}

fn hap_dump(
    i: usize,
    h: &Haplotype,
    ref_hap: &Haplotype,
    ref_bytes: &[u8],
    pad: u64,
    hap_events: &[VariationEvent],
) -> Value {
    let cigar_only = cigar_only_events(h, ref_hap, ref_bytes, pad);
    let hap_base = hap_base_at_ref_locus(h, pad, TARGET);
    json!({
        "idx": i,
        "hash": fnv1a64_hex(&h.bases),
        "len": h.bases.len(),
        "is_reference": h.is_reference,
        "score": h.score,
        "cigar": cigar_str(h),
        "cigar_has_indel": cigar_has_indel(h),
        "align_start": h.alignment_start_hap_wrt_ref,
        "genome_loc": h.genome_loc.map(|g| format!("{}-{}", g.start_1based(), g.end_1based())),
        "hap_base_at_target": hap_base.map(|b| (b as char).to_string()),
        "eventmap_at_or_overlapping_target": events_at(hap_events, TARGET),
        "cigar_replay_at_or_overlapping_target": events_at(&cigar_only, TARGET),
        "has_at": hap_events.iter().any(is_at),
        "cigar_replay_has_at": cigar_only.iter().any(is_at),
    })
}

#[test]
fn holdout_6r121_eventmap_at_29456196() {
    if std::env::var("HOLDOUT_6R121").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R121=1");
        return;
    }
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    let java_vcf = root.join(JAVA_VCF_REL);
    let rust_vcf = root.join(RUST_VCF_REL);
    assert!(ref_fasta.is_file(), "missing {}", ref_fasta.display());
    assert!(bam.is_file(), "missing {}", bam.display());

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
        .expect("ActiveFull containing 20:29456196");

    let args = CallRegionArgs::strict_java();
    let outcome = HaplotypeCallerEngine::call_region(covering, &dict, &ref_fasta, &args)
        .expect("call")
        .expect("ActiveFull outcome");

    let pad = outcome.assembly.padded_reference_start_1based();
    let ref_bytes = outcome.assembly.reference_bases();
    let haps = &outcome.assembly.haplotypes;
    let ref_hap = haps
        .iter()
        .find(|h| h.is_reference)
        .expect("reference haplotype");
    let hap_cache = build_per_haplotype_variation_events(
        haps,
        ref_bytes,
        pad,
        outcome.assembly.max_mnp_distance(),
        "20",
    );

    let union_at: Vec<Value> = outcome
        .assembly
        .variation_events()
        .iter()
        .filter(|e| e.start_1based.get() == TARGET)
        .map(event_json)
        .collect();
    let overlapping: Vec<Value> = outcome
        .assembly
        .variation_events()
        .iter()
        .filter(|e| e.start_1based.get() <= TARGET && e.end_1based.get() >= TARGET)
        .map(event_json)
        .collect();

    let mut per_hap = Vec::new();
    let mut n_hap_at = 0usize;
    let mut n_cigar_replay_at = 0usize;
    for (i, h) in haps.iter().enumerate() {
        let evs = hap_cache.events_for(i);
        if evs.iter().any(is_at) {
            n_hap_at += 1;
        }
        let cigar_only = cigar_only_events(h, ref_hap, ref_bytes, pad);
        if cigar_only.iter().any(is_at) {
            n_cigar_replay_at += 1;
        }
        per_hap.push(hap_dump(i, h, ref_hap, ref_bytes, pad, evs));
    }

    let ref_base = {
        let off = TARGET.saturating_sub(pad) as usize;
        ref_bytes.get(off).copied().map(|b| (b as char).to_string())
    };

    let collect_replay_union: Vec<VariationEvent> = haps
        .iter()
        .flat_map(|h| variation_events_for_haplotype(h, ref_hap, ref_bytes, pad, 0, "20"))
        .filter(is_at)
        .collect();

    let at_start = variation_events_at_position_from_cache(&hap_cache, TARGET, false);
    let at_spanning = variation_events_at_position_from_cache(&hap_cache, TARGET, true);
    let overlap = overlapping_events(outcome.assembly.variation_events(), TARGET);
    let replaced = replace_span_del_events(&at_spanning, TARGET, pad, ref_bytes);
    let merge_b = merged_alleles_for_genotyping(&replaced, TARGET);
    let merge_single = if replaced.len() == 1 {
        Some(format!(
            "{}/{}",
            replaced[0].ref_allele, replaced[0].alt_allele
        ))
    } else {
        None
    };

    let doc = json!({
        "target": TARGET,
        "vcf": {
            "java": vcf_record(&java_vcf, TARGET),
            "rust": vcf_record(&rust_vcf, TARGET),
        },
        "active": [covering.start.get(), covering.end.get()],
        "extended": [covering.extended_start.get(), covering.extended_end.get()],
        "n_reads": covering.reads.len(),
        "pad": pad,
        "ref_len": ref_bytes.len(),
        "ref_base_at_target": ref_base,
        "hap_count": haps.len(),
        "union_events_starting_at_target": union_at,
        "union_events_overlapping_target": overlapping,
        "union_has_at": outcome.assembly.variation_events().iter().any(is_at),
        "n_hap_with_at": n_hap_at,
        "n_hap_cigar_replay_has_at": n_cigar_replay_at,
        "variation_events_for_haplotype_at_n": collect_replay_union.len(),
        "next_object_after_eventmap": {
            "events_starting_at_loc": at_start.iter().map(event_json).collect::<Vec<_>>(),
            "events_spanning_loc": at_spanning.iter().map(event_json).collect::<Vec<_>>(),
            "overlapping_union": overlap.iter().map(event_json).collect::<Vec<_>>(),
            "after_replace_span_dels": replaced.iter().map(event_json).collect::<Vec<_>>(),
            "merged_alleles_for_genotyping": merge_b.map(|(r, a)| json!({"ref": r, "alts": a})),
            "single_event_merged": merge_single,
        },
        "haps_with_at_or_t_base": per_hap.iter().filter(|h| {
            h["has_at"].as_bool() == Some(true)
                || h["hap_base_at_target"].as_str() == Some("T")
        }).cloned().collect::<Vec<_>>(),
        "all_haps": per_hap,
    });
    println!("{}", serde_json::to_string_pretty(&doc).unwrap());
}
