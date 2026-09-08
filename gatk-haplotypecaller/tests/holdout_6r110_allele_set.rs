//! 6R.110 forensic dump: HOLDOUT_6R53 remainder `20:29455649 T/TGTTTG` vs `T/TGTTTG,TTTG`.
//! Inventory A–D before GLs. Skipped unless `HOLDOUT_6R110=1`.
//! Does not reopen 6R.107 / 6R.109.
//!
//! ```text
//! HOLDOUT_6R110=1 cargo test -p gatk-haplotypecaller --test holdout_6r110_allele_set -- --nocapture
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::{
    build_per_haplotype_variation_events, merged_alleles_for_genotyping,
    merged_biallelic_sites_at_position, overlapping_events,
    variation_events_at_position_from_cache, VariationEvent,
};
use gatk_haplotypecaller::hc_allele_mapping::{
    create_allele_mapper_with_events, replace_span_del_events, SPAN_DEL_ALLELE,
};
use gatk_haplotypecaller::hc_genotyping_engine::DEFAULT_STAND_EMIT_CONFIDENCE;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs, Haplotype,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
};
use serde_json::{json, Value};
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const JAVA_VCF_REL: &str = "parity/reports/6r43/chr20_tiny/java.vcf";
const RUST_VCF_REL: &str = "parity/reports/6r43/chr20_tiny/rust.vcf";
const REGION: (u64, u64) = (29_455_560, 29_455_744);
const TARGET: u64 = 29_455_649;
const NEARBY_INS: u64 = 29_455_644;

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
        "len_ref": e.ref_allele.len(),
        "len_alt": e.alt_allele.len(),
        "kind": if e.alt_allele == SPAN_DEL_ALLELE {
            "span_del"
        } else if e.is_indel() {
            "indel"
        } else {
            "snp"
        },
    })
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
            "pos": pos,
            "alleles": format!("{}/{}", f[3], f[4]),
            "qual": f[5],
        }));
    }
    None
}

fn is_tgt(e: &VariationEvent) -> bool {
    e.start_1based.get() == TARGET && e.ref_allele == "T" && e.alt_allele == "TGTTTG"
}

fn is_tttg(e: &VariationEvent) -> bool {
    e.start_1based.get() == TARGET && e.ref_allele == "T" && e.alt_allele == "TTTG"
}

fn is_aat(e: &VariationEvent) -> bool {
    e.start_1based.get() == NEARBY_INS && e.ref_allele == "A" && e.alt_allele == "AT"
}

#[test]
fn holdout_6r110_allele_set_dump() {
    if std::env::var("HOLDOUT_6R110").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R110=1");
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
            ) && r.start.get() == REGION.0
                && r.end.get() == REGION.1
        })
        .expect("ActiveFull");
    let args = CallRegionArgs::strict_java();
    let outcome = HaplotypeCallerEngine::call_region(covering, &dict, &ref_fasta, &args)
        .expect("call")
        .expect("outcome");
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

    let union_at: Vec<Value> = outcome
        .assembly
        .variation_events()
        .iter()
        .filter(|e| e.start_1based.get() == TARGET)
        .map(event_json)
        .collect();
    let union_has_tgt = outcome.assembly.variation_events().iter().any(is_tgt);
    let union_has_tttg = outcome.assembly.variation_events().iter().any(is_tttg);
    let union_has_aat = outcome.assembly.variation_events().iter().any(is_aat);

    let mut per_hap = Vec::new();
    let mut n_tgt = 0usize;
    let mut n_tttg = 0usize;
    let mut n_both = 0usize;
    let mut n_aat_and_tttg = 0usize;
    for (i, h) in haps.iter().enumerate() {
        let evs: Vec<Value> = hap_events
            .events_for(i)
            .iter()
            .filter(|e| e.start_1based.get() == TARGET)
            .map(event_json)
            .collect();
        let has_tgt = hap_events.events_for(i).iter().any(is_tgt);
        let has_tttg = hap_events.events_for(i).iter().any(is_tttg);
        let has_aat = hap_events.events_for(i).iter().any(is_aat);
        if has_tgt {
            n_tgt += 1;
        }
        if has_tttg {
            n_tttg += 1;
        }
        if has_tgt && has_tttg {
            n_both += 1;
        }
        if has_aat && has_tttg {
            n_aat_and_tttg += 1;
        }
        per_hap.push(json!({
            "idx": i,
            "hash": fnv1a64_hex(&h.bases),
            "len": h.bases.len(),
            "is_reference": h.is_reference,
            "score": h.score,
            "cigar": cigar_str(h),
            "align_start": h.alignment_start_hap_wrt_ref,
            "has_TGTTTG": has_tgt,
            "has_TTTG": has_tttg,
            "has_A_AT": has_aat,
            "events_at_649": evs,
        }));
    }

    let at_start = variation_events_at_position_from_cache(&hap_events, TARGET, false);
    let at_spanning = variation_events_at_position_from_cache(&hap_events, TARGET, true);
    let overlap = overlapping_events(outcome.assembly.variation_events(), TARGET);
    let replaced = replace_span_del_events(&at_spanning, TARGET, pad, ref_bytes);
    let merge_b = merged_alleles_for_genotyping(&replaced, TARGET);
    let biallelics = merged_biallelic_sites_at_position(&replaced, TARGET);

    let event_tgt = VariationEvent::from_alleles("20", TARGET, "T", "TGTTTG");
    let event_tttg = VariationEvent::from_alleles("20", TARGET, "T", "TTTG");
    let map_tgt = create_allele_mapper_with_events(
        &event_tgt,
        TARGET,
        haps,
        pad,
        ref_bytes,
        outcome.assembly.max_mnp_distance(),
        true,
        Some(&hap_events),
    );
    let map_tttg = create_allele_mapper_with_events(
        &event_tttg,
        TARGET,
        haps,
        pad,
        ref_bytes,
        outcome.assembly.max_mnp_distance(),
        true,
        Some(&hap_events),
    );

    let live: Vec<Value> = outcome
        .genotyped_calls
        .iter()
        .filter(|c| c.event.start_1based.get() == TARGET)
        .map(|c| {
            json!({
                "ref": c.event.ref_allele,
                "alt": c.event.alt_allele,
                "extra_alts": c.extra_alt_alleles,
                "ad": c.genotype.format.ad_as_i32(),
                "pl": c.genotype.format.pl_as_i32(),
                "post_merge_unused_alt_subset": c.post_merge_unused_alt_subset,
            })
        })
        .collect();
    let emitted =
        try_emit_call_region_variants(covering, &outcome, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .unwrap_or_default();
    let emit_c: Vec<Value> = emitted
        .iter()
        .filter(|r| r.position == TARGET)
        .map(|r| {
            json!({
                "pos": r.position,
                "ref": r.reference,
                "alts": r.alternate,
            })
        })
        .collect();

    let unique_alts: BTreeSet<String> = at_start
        .iter()
        .map(|e| format!("{}/{}", e.ref_allele, e.alt_allele))
        .collect();

    let doc = json!({
        "holdout": "20:29455649 T/TGTTTG vs T/TGTTTG,TTTG",
        "active_full": REGION,
        "hap_count": haps.len(),
        "frozen_vcf": {
            "java_649": vcf_record(&java_vcf, TARGET),
            "rust_649": vcf_record(&rust_vcf, TARGET),
        },
        "inventory_a_eventmap": {
            "union_at_649": union_at,
            "union_has_T_TGTTTG": union_has_tgt,
            "union_has_T_TTTG": union_has_tttg,
            "union_has_A_AT_644": union_has_aat,
            "events_at_start_only": at_start.iter().map(event_json).collect::<Vec<_>>(),
            "events_with_spanning": at_spanning.iter().map(event_json).collect::<Vec<_>>(),
            "overlap_union": overlap.iter().filter(|e| e.start_1based.get() == TARGET).map(event_json).collect::<Vec<_>>(),
            "unique_start_alleles": unique_alts,
        },
        "inventory_b_haplotypes": {
            "n_with_TGTTTG": n_tgt,
            "n_with_TTTG": n_tttg,
            "n_with_both": n_both,
            "n_A_AT_and_TTTG": n_aat_and_tttg,
        },
        "inventory_c_d_merge": {
            "after_replace_span_dels": replaced.iter().map(event_json).collect::<Vec<_>>(),
            "merged_alleles_for_genotyping": merge_b,
            "merged_biallelic_sites": biallelics.iter().map(event_json).collect::<Vec<_>>(),
        },
        "inventory_e_mapper": {
            "TGTTTG_n_ref": map_tgt.ref_haplotype_indices.len(),
            "TGTTTG_n_alt": map_tgt.alt_haplotype_indices.len(),
            "TGTTTG_alt_idx": map_tgt.alt_haplotype_indices.iter().map(|i| i.get()).collect::<Vec<_>>(),
            "TTTG_n_ref": map_tttg.ref_haplotype_indices.len(),
            "TTTG_n_alt": map_tttg.alt_haplotype_indices.len(),
            "TTTG_alt_idx": map_tttg.alt_haplotype_indices.iter().map(|i| i.get()).collect::<Vec<_>>(),
        },
        "live_genotype_at_649": live,
        "live_emit_649": emit_c,
        "inventory_c_d_merge_uses_joint_gls": merge_b.as_ref().map(|(r, a)| {
            gatk_haplotypecaller::event_map::merged_site_uses_joint_gls(&replaced, TARGET, r, a)
        }),
        "per_hap": per_hap,
    });
    println!("{}", serde_json::to_string_pretty(&doc).unwrap());

    assert!(
        vcf_record(&java_vcf, TARGET)
            .as_ref()
            .is_some_and(|v| v["alleles"] == "T/TGTTTG"),
        "frozen Java VCF is T/TGTTTG"
    );
    assert_eq!(
        live.len(),
        1,
        "same-REF multi-indel is one joint genotyping call, not two biallelics: {live:?}"
    );
    assert!(
        !live
            .iter()
            .any(|c| c["alt"] == "TTTG" && c["extra_alts"] == json!([])),
        "unused insertion must not remain as a standalone genotyped site"
    );
    assert_eq!(emit_c.len(), 1, "one VCF record at the loc: {emit_c:?}");
    assert_eq!(
        emit_c[0]["alts"],
        json!(["TGTTTG"]),
        "Java calculateOutputAlleleSubset keeps TGTTTG only: {emit_c:?}"
    );
}
