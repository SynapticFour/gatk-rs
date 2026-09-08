//! 6R.122 forensic dump: `create_allele_mapper` at `20:29456196 A/T`.
//! Mapper only. Skipped unless `HOLDOUT_6R122=1`. Production change: NONE.
//!
//! ```text
//! HOLDOUT_6R122=1 cargo test -p gatk-haplotypecaller --test holdout_6r122_mapper -- --nocapture
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::{
    build_per_haplotype_variation_events, overlapping_events, VariationEvent,
};
use gatk_haplotypecaller::hc_allele_mapping::{
    create_allele_mapper_with_events, hap_base_at_ref_locus,
};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    AssemblyRegionCallDisposition, CallRegionArgs, Haplotype, HaplotypeCallerEngine,
    ReadFilterParams, WalkerTraversalConfig,
};
use serde_json::{json, Value};
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const JAVA_VCF_REL: &str = "parity/reports/6r43/chr20_tiny/java.vcf";
const RUST_VCF_REL: &str = "parity/reports/6r43/chr20_tiny/rust.vcf";
const TARGET: u64 = 29_456_196;
const EMIT_SPANNING_DELS: bool = true;

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
        return Some(json!({
            "alleles": format!("{}/{}", f[3], f[4]),
            "qual": f[5],
        }));
    }
    None
}

fn mapper_reason(spanning: &[VariationEvent], hap_base: Option<char>, assigned: &str) -> String {
    if spanning.is_empty() {
        return match (assigned, hap_base) {
            ("ALT", Some('T')) => "empty_span_snp_base_eq_alt".to_string(),
            ("REF", Some('A')) => "empty_span_snp_base_eq_ref_or_java_empty_ref".to_string(),
            ("REF", _) => "empty_span_default_ref".to_string(),
            _ => format!("empty_span_assigned_{assigned}_base_{hap_base:?}"),
        };
    }
    let loc = TARGET;
    let mut has_exact = false;
    let mut has_prior = false;
    let mut other = Vec::new();
    for ev in spanning {
        let s = ev.start_1based.get();
        if s == loc && ev.ref_allele == "A" && ev.alt_allele == "T" {
            has_exact = true;
        } else if s < loc {
            has_prior = true;
            other.push(format!("{}:{}/{}", s, ev.ref_allele, ev.alt_allele));
        } else {
            other.push(format!("{}:{}/{}", s, ev.ref_allele, ev.alt_allele));
        }
    }
    if has_exact && assigned == "ALT" {
        return "eventmap_start_eq_loc_exact_A/T".to_string();
    }
    if has_prior {
        return format!("overlapping_prior_start other={other:?} assigned={assigned}");
    }
    format!("spanning_unclassified other={other:?} assigned={assigned}")
}

fn java_predicate_from_spanning(spanning: &[VariationEvent]) -> &'static str {
    if spanning.is_empty() {
        return "REF (empty overlapping EventMap)";
    }
    for ev in spanning {
        if ev.start_1based.get() == TARGET {
            if ev.ref_allele == "A" && ev.alt_allele == "T" {
                return "ALT T (start==loc and alt in mergedVC)";
            }
            return "unmapped or skip (start==loc allele not in merged keys / longer REF)";
        }
        if EMIT_SPANNING_DELS {
            return "SPAN_DEL * (start<loc, emitSpanningDels)";
        }
        return "REF (start<loc, !emitSpanningDels)";
    }
    "unmapped (no matching spanning event)"
}

#[test]
fn holdout_6r122_create_allele_mapper() {
    if std::env::var("HOLDOUT_6R122").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R122=1");
        return;
    }
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
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
        .expect("ActiveFull containing 20:29456196");
    let args = CallRegionArgs::strict_java();
    let outcome = HaplotypeCallerEngine::call_region(covering, &dict, &ref_fasta, &args)
        .expect("call")
        .expect("ActiveFull outcome");
    let pad = outcome.assembly.padded_reference_start_1based();
    let ref_bytes = outcome.assembly.reference_bases();
    let haps = &outcome.assembly.haplotypes;
    let hap_cache = build_per_haplotype_variation_events(
        haps,
        ref_bytes,
        pad,
        outcome.assembly.max_mnp_distance(),
        "20",
    );
    let merged = VariationEvent::from_alleles("20", TARGET, "A", "T");
    let mapping = create_allele_mapper_with_events(
        &merged,
        TARGET,
        haps,
        pad,
        ref_bytes,
        outcome.assembly.max_mnp_distance(),
        EMIT_SPANNING_DELS,
        Some(&hap_cache),
    );
    let ref_set: BTreeSet<usize> = mapping
        .ref_haplotype_indices
        .iter()
        .map(|i| i.get())
        .collect();
    let alt_set: BTreeSet<usize> = mapping
        .alt_haplotype_indices
        .iter()
        .map(|i| i.get())
        .collect();

    let mut rows = Vec::new();
    let mut n_unmapped = 0usize;
    let mut n_empty_span_alt = 0usize;
    let mut n_prior_span = 0usize;
    for (i, h) in haps.iter().enumerate() {
        let spanning = overlapping_events(hap_cache.events_for(i), TARGET);
        let hap_base = hap_base_at_ref_locus(h, pad, TARGET).map(|b| b as char);
        let assigned = if alt_set.contains(&i) {
            "ALT"
        } else if ref_set.contains(&i) {
            "REF"
        } else {
            n_unmapped += 1;
            "unmapped"
        };
        if spanning.is_empty() && assigned == "ALT" {
            n_empty_span_alt += 1;
        }
        if spanning.iter().any(|e| e.start_1based.get() < TARGET) {
            n_prior_span += 1;
        }
        rows.push(json!({
            "idx": i,
            "hash": fnv1a64_hex(&h.bases),
            "len": h.bases.len(),
            "is_reference": h.is_reference,
            "cigar": cigar_str(h),
            "hap_base": hap_base.map(|c| c.to_string()),
            "spanning": spanning.iter().map(|e| format!(
                "{}:{}/{}",
                e.start_1based.get(), e.ref_allele, e.alt_allele
            )).collect::<Vec<_>>(),
            "assigned": assigned,
            "rust_reason": mapper_reason(&spanning, hap_base, assigned),
            "java_predicate_if_same_spanning": java_predicate_from_spanning(&spanning),
        }));
    }

    let doc = json!({
        "target": TARGET,
        "active": [covering.start.get(), covering.end.get()],
        "vcf": {
            "java": vcf_record(&root.join(JAVA_VCF_REL), TARGET),
            "rust": vcf_record(&root.join(RUST_VCF_REL), TARGET),
        },
        "mapper_input": {
            "merged_ref": mapping.ref_allele,
            "merged_alt": mapping.alt_allele,
            "allele_order": ["A", "T"],
            "emit_spanning_dels": EMIT_SPANNING_DELS,
        },
        "pool_counts": {
            "hap_count": haps.len(),
            "ref": mapping.ref_haplotype_indices.len(),
            "alt": mapping.alt_haplotype_indices.len(),
            "unmapped": n_unmapped,
            "empty_span_assigned_alt": n_empty_span_alt,
            "haps_with_prior_start_overlap": n_prior_span,
        },
        "haps": rows,
    });
    println!("{}", serde_json::to_string_pretty(&doc).unwrap());
}
