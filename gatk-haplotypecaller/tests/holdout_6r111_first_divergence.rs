//! 6R.111 forensic dump: HOLDOUT 6R.53 remainder at `20:29455745`.
//! Inventory A (EventMap) first. Skipped unless `HOLDOUT_6R111=1`.
//! Does not reopen 6R.107–6R.110.
//!
//! ```text
//! HOLDOUT_6R111=1 cargo test -p gatk-haplotypecaller --test holdout_6r111_first_divergence -- --nocapture
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
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const JAVA_VCF_REL: &str = "parity/reports/6r43/chr20_tiny/java.vcf";
const RUST_VCF_REL: &str = "parity/reports/6r43/chr20_tiny/rust.vcf";
const PREV_REGION: (u64, u64) = (29_455_560, 29_455_744);
const REGION: (u64, u64) = (29_455_745, 29_455_993);
const TARGET: u64 = 29_455_745;
const WINDOW: (u64, u64) = (29_455_725, 29_455_765);
const CLOSED_644: u64 = 29_455_644;
const CLOSED_649: u64 = 29_455_649;

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

fn vcf_records_in(path: &Path, lo: u64, hi: u64) -> Vec<Value> {
    let mut out = Vec::new();
    if !path.is_file() {
        return out;
    }
    for line in std::fs::read_to_string(path).unwrap_or_default().lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<_> = line.split('\t').collect();
        if f.len() < 10 || f[0] != "20" {
            continue;
        }
        let Ok(pos) = f[1].parse::<u64>() else {
            continue;
        };
        if pos < lo || pos > hi {
            continue;
        }
        out.push(json!({
            "pos": pos,
            "ref": f[3],
            "alts": f[4],
            "qual": f[5],
            "format": f[8],
            "sample": f[9],
        }));
    }
    out
}

#[test]
fn holdout_6r111_first_divergence_dump() {
    if std::env::var("HOLDOUT_6R111").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R111=1");
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
    let active: Vec<Value> = regions
        .iter()
        .filter(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            )
        })
        .map(|r| {
            json!({
                "start": r.start.get(),
                "end": r.end.get(),
                "covers_target": r.start.get() <= TARGET && r.end.get() >= TARGET,
            })
        })
        .collect();
    let covering = regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() == REGION.0
                && r.end.get() == REGION.1
        })
        .expect("ActiveFull 20:29455745-29455993");
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
    let ref_at_target = {
        let off = TARGET.saturating_sub(pad) as usize;
        ref_bytes
            .get(off)
            .copied()
            .map(|b| (b as char).to_string())
            .unwrap_or_else(|| ".".to_string())
    };

    let union_window: Vec<Value> = outcome
        .assembly
        .variation_events()
        .iter()
        .filter(|e| {
            let s = e.start_1based.get();
            (s >= WINDOW.0 && s <= WINDOW.1)
                || (e.start_1based.get() <= TARGET && e.end_1based.get() >= TARGET)
        })
        .map(event_json)
        .collect();
    let union_at_loc: Vec<Value> = outcome
        .assembly
        .variation_events()
        .iter()
        .filter(|e| e.start_1based.get() == TARGET)
        .map(event_json)
        .collect();
    let overlap = overlapping_events(outcome.assembly.variation_events(), TARGET);
    let at_start = variation_events_at_position_from_cache(&hap_events, TARGET, false);
    let at_spanning = variation_events_at_position_from_cache(&hap_events, TARGET, true);
    let replaced = replace_span_del_events(&at_spanning, TARGET, pad, ref_bytes);
    let merge_b = merged_alleles_for_genotyping(&replaced, TARGET);
    let biallelics = merged_biallelic_sites_at_position(&replaced, TARGET);

    let mut per_hap = Vec::new();
    let mut n_with_event_at_loc = 0usize;
    let mut n_with_spanning = 0usize;
    for (i, h) in haps.iter().enumerate() {
        let evs: Vec<&VariationEvent> = hap_events
            .events_for(i)
            .iter()
            .filter(|e| {
                e.start_1based.get() == TARGET
                    || (e.start_1based.get() <= TARGET && e.end_1based.get() >= TARGET)
            })
            .collect();
        let at = evs.iter().any(|e| e.start_1based.get() == TARGET);
        let span = !evs.is_empty();
        if at {
            n_with_event_at_loc += 1;
        }
        if span {
            n_with_spanning += 1;
        }
        if evs.is_empty() && !h.is_reference {
            continue;
        }
        per_hap.push(json!({
            "idx": i,
            "hash": fnv1a64_hex(&h.bases),
            "len": h.bases.len(),
            "is_reference": h.is_reference,
            "score": h.score,
            "cigar": cigar_str(h),
            "align_start": h.alignment_start_hap_wrt_ref,
            "events_at_or_span_target": evs.iter().copied().map(event_json).collect::<Vec<_>>(),
        }));
    }

    let mapper_rows: Vec<Value> = at_start
        .iter()
        .map(|e| {
            let mapping = create_allele_mapper_with_events(
                e,
                TARGET,
                haps,
                pad,
                ref_bytes,
                outcome.assembly.max_mnp_distance(),
                true,
                Some(&hap_events),
            );
            json!({
                "event": event_json(e),
                "n_ref": mapping.ref_haplotype_indices.len(),
                "n_alt": mapping.alt_haplotype_indices.len(),
                "alt_idx": mapping.alt_haplotype_indices.iter().map(|i| i.get()).collect::<Vec<_>>(),
            })
        })
        .collect();

    let live: Vec<Value> = outcome
        .genotyped_calls
        .iter()
        .filter(|c| {
            let p = c.event.start_1based.get();
            p == TARGET || (p >= WINDOW.0 && p <= WINDOW.1)
        })
        .map(|c| {
            json!({
                "pos": c.event.start_1based.get(),
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
    let emit_window: Vec<Value> = emitted
        .iter()
        .filter(|r| r.position >= WINDOW.0 && r.position <= WINDOW.1)
        .map(|r| {
            json!({
                "pos": r.position,
                "ref": r.reference,
                "alts": r.alternate,
                "qual": r.quality,
                "gt": r.samples.first().and_then(|s| s.gt.as_ref().map(|g| g.alleles.clone())),
                "ad": r.samples.first().and_then(|s| s.ad.clone()),
                "pl": r.samples.first().and_then(|s| s.pl.clone()),
            })
        })
        .collect();
    let emit_at_loc: Vec<Value> = emitted
        .iter()
        .filter(|r| r.position == TARGET)
        .map(|r| {
            json!({
                "pos": r.position,
                "ref": r.reference,
                "alts": r.alternate,
                "qual": r.quality,
            })
        })
        .collect();

    let prev = regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() == PREV_REGION.0
                && r.end.get() == PREV_REGION.1
        })
        .expect("ActiveFull 20:29455560-29455744");
    let prev_outcome = HaplotypeCallerEngine::call_region(prev, &dict, &ref_fasta, &args)
        .expect("prev call")
        .expect("prev outcome");
    let prev_emitted =
        try_emit_call_region_variants(prev, &prev_outcome, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .unwrap_or_default();
    let live_644 = prev_emitted.iter().any(|r| {
        r.position == CLOSED_644 && r.reference == "A" && r.alternate.iter().any(|a| a == "AT")
    });
    let live_649: Vec<String> = prev_emitted
        .iter()
        .filter(|r| r.position == CLOSED_649)
        .map(|r| format!("{}/{}", r.reference, r.alternate.join(",")))
        .collect();

    assert!(
        vcf_records_in(&java_vcf, TARGET, TARGET).is_empty(),
        "frozen Java VCF must not contain 20:29455745"
    );
    assert!(
        union_at_loc.is_empty() && at_spanning.is_empty() && overlap.is_empty(),
        "EventMap at 20:29455745 must be empty (Java union/spanning/merged_vc also empty)"
    );
    assert!(
        emit_at_loc.is_empty() && live.is_empty(),
        "Rust must not genotype or emit 20:29455745: emit={emit_at_loc:?} gt={live:?}"
    );
    assert!(
        !live_644,
        "6R.109 closed: 20:29455644 A/AT must stay absent"
    );
    assert_eq!(
        live_649,
        vec!["T/TGTTTG".to_string()],
        "6R.110 closed: 20:29455649 must stay T/TGTTTG: {live_649:?}"
    );

    let doc = json!({
        "holdout": "20:29455745",
        "prev_active_full": PREV_REGION,
        "active_full": REGION,
        "active_fulls": active,
        "ref_at_target": ref_at_target,
        "hap_count": haps.len(),
        "union_n": outcome.assembly.variation_events().len(),
        "pad": pad,
        "step0_frozen": {
            "java_window": vcf_records_in(&java_vcf, WINDOW.0, WINDOW.1),
            "rust_window": vcf_records_in(&rust_vcf, WINDOW.0, WINDOW.1),
            "java_at_loc": vcf_records_in(&java_vcf, TARGET, TARGET),
            "rust_at_loc": vcf_records_in(&rust_vcf, TARGET, TARGET),
            "java_closed_644": vcf_records_in(&java_vcf, CLOSED_644, CLOSED_644),
            "java_closed_649": vcf_records_in(&java_vcf, CLOSED_649, CLOSED_649),
        },
        "inventory_a_eventmap": {
            "union_window_or_span": union_window,
            "union_at_loc": union_at_loc,
            "overlap_union": overlap.iter().map(event_json).collect::<Vec<_>>(),
            "events_at_start_only": at_start.iter().map(event_json).collect::<Vec<_>>(),
            "events_with_spanning": at_spanning.iter().map(event_json).collect::<Vec<_>>(),
        },
        "inventory_b_haplotypes": {
            "n_with_event_at_loc": n_with_event_at_loc,
            "n_with_spanning_or_at": n_with_spanning,
            "relevant_haps": per_hap,
        },
        "inventory_c_d_merge": {
            "after_replace_span_dels": replaced.iter().map(event_json).collect::<Vec<_>>(),
            "merged_alleles_for_genotyping": merge_b,
            "merged_biallelic_sites": biallelics.iter().map(event_json).collect::<Vec<_>>(),
        },
        "inventory_e_mapper": mapper_rows,
        "live_genotype_window": live,
        "live_emit_window": emit_window,
        "live_emit_at_loc": emit_at_loc,
        "closed_sites": {
            "live_644_emitted": live_644,
            "live_649": live_649,
        },
    });
    println!("{}", serde_json::to_string_pretty(&doc).unwrap());
}
