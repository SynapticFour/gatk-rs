//! 6R.104 forensic dump: allele inventories at HOLDOUT_6R53 `20:29455388 C/T`.
//!
//! Distinguishes inventories A (EventMap union), B (merged genotyping alleles),
//! C (emitted VCF). Skipped unless `HOLDOUT_6R104=1`.
//!
//! Original 6R.104 finding (pre-6R.107): A and B already contain `C/T` on both
//! sides; inventory C was the first unequal object (Java VCF absent, Rust extra
//! `C/T`). 6R.106/6R.107 proved that extra was the L9 SparsePlShape overwrite of
//! calculator `PL=0,6,1780`. After HomAltStrong gating, inventory C matches Java:
//! neither VCF emits the site. EventMap presence is not VCF emission.
//!
//! ```text
//! HOLDOUT_6R104=1 cargo test -p gatk-haplotypecaller --test holdout_6r104_holdout_allele_set -- --nocapture
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::{
    build_per_haplotype_variation_events, merged_alleles_for_genotyping, overlapping_events,
    variation_events_at_position_from_cache, VariationEvent,
};
use gatk_haplotypecaller::hc_allele_mapping::{
    hap_base_at_ref_locus, replace_span_del_events, SPAN_DEL_ALLELE,
};
use gatk_haplotypecaller::hc_genotyping_engine::DEFAULT_STAND_EMIT_CONFIDENCE;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs, Haplotype,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
};
use serde_json::{json, Value};
use std::collections::hash_map::DefaultHasher;
use std::collections::BTreeSet;
use std::hash::{Hash, Hasher};
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const JAVA_VCF_REL: &str = "parity/reports/6r43/chr20_tiny/java.vcf";
const RUST_VCF_REL: &str = "parity/reports/6r43/chr20_tiny/rust.vcf";
const COVERING: (u64, u64) = (29_455_300, 29_455_559);
const TARGET: u64 = 29_455_388;
const WINDOW: (u64, u64) = (29_455_370, 29_455_400);

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn parse_vcf_keys(path: &Path) -> BTreeSet<(u64, String, String)> {
    let mut out = BTreeSet::new();
    if !path.is_file() {
        return out;
    }
    for line in std::fs::read_to_string(path).unwrap_or_default().lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<_> = line.split('\t').collect();
        if f.len() < 5 || f[0] != "20" {
            continue;
        }
        out.insert((f[1].parse().unwrap(), f[3].to_string(), f[4].to_string()));
    }
    out
}

fn vcf_record_at(path: &Path, pos: u64) -> Option<Value> {
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

fn hap_hash(h: &Haplotype) -> String {
    let mut s = DefaultHasher::new();
    h.bases.hash(&mut s);
    format!("{:016x}", s.finish())
}

fn cigar_str(h: &Haplotype) -> String {
    h.cigar
        .as_ref()
        .map(|c| {
            c.elements
                .iter()
                .map(|e| format!("{}{}", e.length, e.operator.as_char()))
                .collect::<String>()
        })
        .unwrap_or_else(|| ".".to_string())
}

fn event_json(e: &VariationEvent) -> Value {
    json!({
        "start": e.start_1based.get(),
        "end": e.end_1based.get(),
        "ref": e.ref_allele,
        "alt": e.alt_allele,
        "indel": e.is_indel(),
        "kind": if e.alt_allele == SPAN_DEL_ALLELE {
            "span_del"
        } else if e.is_indel() {
            "indel"
        } else {
            "snp"
        },
    })
}

fn in_window(e: &VariationEvent) -> bool {
    let s = e.start_1based.get();
    s >= WINDOW.0 && s <= WINDOW.1
}

#[test]
fn holdout_6r104_first_allele_set_divergence() {
    if std::env::var("HOLDOUT_6R104").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R104=1");
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
    let covering = regions_covering(&regions);
    let args = CallRegionArgs::strict_java();
    let outcome = HaplotypeCallerEngine::call_region(covering, &dict, &ref_fasta, &args)
        .expect("call")
        .expect("ActiveFull outcome");

    let pad = outcome.assembly.padded_reference_start_1based();
    let ref_bytes = outcome.assembly.reference_bases();
    let haps = &outcome.assembly.haplotypes;
    let hap_events = build_per_haplotype_variation_events(
        haps,
        ref_bytes,
        pad,
        outcome.assembly.max_mnp_distance(),
        "20",
    );

    let union_a: Vec<Value> = outcome
        .assembly
        .variation_events()
        .iter()
        .filter(|e| in_window(e) || e.start_1based.get() <= TARGET && e.end_1based.get() >= TARGET)
        .map(event_json)
        .collect();
    let union_has_ct = outcome
        .assembly
        .variation_events()
        .iter()
        .any(|e| e.start_1based.get() == TARGET && e.ref_allele == "C" && e.alt_allele == "T");

    let mut per_hap = Vec::new();
    for (i, h) in haps.iter().enumerate() {
        let base = hap_base_at_ref_locus(h, pad, TARGET)
            .map(|b| (b as char).to_string())
            .unwrap_or_else(|| ".".to_string());
        let evs: Vec<Value> = hap_events
            .events_for(i)
            .iter()
            .filter(|e| {
                in_window(e) || e.start_1based.get() <= TARGET && e.end_1based.get() >= TARGET
            })
            .map(event_json)
            .collect();
        let has_ct = hap_events
            .events_for(i)
            .iter()
            .any(|e| e.start_1based.get() == TARGET && e.ref_allele == "C" && e.alt_allele == "T");
        per_hap.push(json!({
            "idx": i,
            "hash": hap_hash(h),
            "len": h.bases.len(),
            "is_reference": h.is_reference,
            "score": h.score,
            "cigar": cigar_str(h),
            "align_start": h.alignment_start_hap_wrt_ref,
            "base_at_target": base,
            "eventmap_ct": has_ct,
            "events_window_or_span": evs,
        }));
    }

    let at_start_only = variation_events_at_position_from_cache(&hap_events, TARGET, false);
    let at_spanning = variation_events_at_position_from_cache(&hap_events, TARGET, true);
    let overlap = overlapping_events(outcome.assembly.variation_events(), TARGET);
    let replaced = replace_span_del_events(&at_spanning, TARGET, pad, ref_bytes);
    let merge_b = merged_alleles_for_genotyping(&replaced, TARGET);
    let merge_from_overlap = merged_alleles_for_genotyping(
        &replace_span_del_events(&overlap, TARGET, pad, ref_bytes),
        TARGET,
    );

    let calls_at: Vec<Value> = outcome
        .genotyped_calls
        .iter()
        .filter(|c| {
            c.event.start_1based.get() == TARGET
                || (c.event.start_1based.get() <= TARGET && c.event.end_1based.get() >= TARGET)
        })
        .map(|c| {
            json!({
                "start": c.event.start_1based.get(),
                "end": c.event.end_1based.get(),
                "ref": c.event.ref_allele,
                "alt": c.event.alt_allele,
                "extra_alts": c.extra_alt_alleles,
                "post_merge_unused_alt_subset": c.post_merge_unused_alt_subset,
                "ad": c.genotype.format.ad_as_i32(),
                "pl": c.genotype.format.pl_as_i32(),
                "gq": c.genotype.format.gq.as_i32(),
            })
        })
        .collect();

    let emitted =
        try_emit_call_region_variants(covering, &outcome, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .unwrap_or_default();
    let emit_c: Vec<Value> = emitted
        .iter()
        .filter(|r| r.position >= WINDOW.0 && r.position <= WINDOW.1)
        .map(|r| {
            json!({
                "pos": r.position,
                "ref": r.reference,
                "alts": r.alternate,
            })
        })
        .collect();

    let n_t = haps
        .iter()
        .filter(|h| hap_base_at_ref_locus(h, pad, TARGET) == Some(b'T'))
        .count();
    let n_c = haps
        .iter()
        .filter(|h| hap_base_at_ref_locus(h, pad, TARGET) == Some(b'C'))
        .count();

    let jk = parse_vcf_keys(&java_vcf);
    let rk = parse_vcf_keys(&rust_vcf);
    let in_span = |k: &(u64, String, String)| k.0 >= COVERING.0 && k.0 <= COVERING.1;
    let cov_j: BTreeSet<_> = jk.iter().filter(|k| in_span(k)).cloned().collect();
    let cov_r: BTreeSet<_> = rk.iter().filter(|k| in_span(k)).cloned().collect();

    let doc = json!({
        "holdout": "20:29455388 C/T",
        "covering": COVERING,
        "n_haps": haps.len(),
        "pad": pad,
        "final_vcf": {
            "java": vcf_record_at(&java_vcf, TARGET),
            "rust": vcf_record_at(&rust_vcf, TARGET),
            "covering_java_only": cov_j.difference(&cov_r).cloned().collect::<Vec<_>>(),
            "covering_rust_only": cov_r.difference(&cov_j).cloned().collect::<Vec<_>>(),
        },
        "hap_bases_at_target": {"C": n_c, "T": n_t},
        "inventory_a_eventmap_union_window_or_span": union_a,
        "inventory_a_has_ct": union_has_ct,
        "per_haplotype": per_hap,
        "inventory_b": {
            "events_at_start_only": at_start_only.iter().map(event_json).collect::<Vec<_>>(),
            "events_with_spanning": at_spanning.iter().map(event_json).collect::<Vec<_>>(),
            "overlapping_union": overlap.iter().map(event_json).collect::<Vec<_>>(),
            "after_replace_span_dels": replaced.iter().map(event_json).collect::<Vec<_>>(),
            "merged_alleles_for_genotyping": merge_b.as_ref().map(|(r, alts)| json!({"ref": r, "alts": alts})),
            "merged_from_overlap": merge_from_overlap.as_ref().map(|(r, alts)| json!({"ref": r, "alts": alts})),
            "colocated_merge_applicable": merge_b.is_some(),
            "star_present_after_replace": replaced.iter().any(|e| e.alt_allele == SPAN_DEL_ALLELE),
            "snp_present_after_replace": replaced.iter().any(|e| e.ref_allele == "C" && e.alt_allele == "T"),
        },
        "inventory_b_genotyped_calls_at_or_spanning": calls_at,
        "inventory_c_emit_window": emit_c,
    });
    println!("{}", serde_json::to_string_pretty(&doc).unwrap());

    assert!(
        vcf_record_at(&java_vcf, TARGET).is_none(),
        "Java final VCF has no record at 20:29455388"
    );
    assert!(
        vcf_record_at(&rust_vcf, TARGET).is_none(),
        "6R.107: inventory C must match Java — no frozen Rust VCF C/T at 20:29455388"
    );
    assert!(
        !emit_c.iter().any(|r| r["pos"] == TARGET),
        "6R.107: live emit window must not contain 20:29455388"
    );
    assert!(
        !cov_r
            .iter()
            .any(|(p, r, a)| *p == TARGET && r == "C" && a == "T"),
        "6R.107: covering rust.vcf must not be rust-only at C/T"
    );
    assert!(
        union_has_ct,
        "Rust EventMap union (inventory A) contains C/T"
    );
    assert!(
        replaced.iter().any(|e| e.alt_allele == SPAN_DEL_ALLELE)
            && replaced
                .iter()
                .any(|e| e.ref_allele == "C" && e.alt_allele == "T"),
        "Rust replaceSpanDels at loc is C/* beside C/T"
    );
    assert!(
        merge_b.is_none(),
        "same-REF SNP+* is not 6R.61 colocated merge"
    );
    assert!(
        calls_at.is_empty(),
        "6R.107: calculator hom-ref must not become a genotyped C/T call: {calls_at:?}"
    );
}

fn regions_covering(
    regions: &[gatk_haplotypecaller::AssemblyRegion],
) -> &gatk_haplotypecaller::AssemblyRegion {
    regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() == COVERING.0
                && r.end.get() == COVERING.1
        })
        .expect("covering ActiveFull")
}
