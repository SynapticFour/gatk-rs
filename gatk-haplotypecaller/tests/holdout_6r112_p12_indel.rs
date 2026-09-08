//! 6R.112 forensic dump: P12 L3a miss `2:92307327 A/ATG`.
//! EventMap first; L9 only if upstream matches. Skipped unless `HOLDOUT_6R112=1`.
//! Does not reopen 6R.107 / 6R.109 / 6R.110.
//!
//! ```text
//! HOLDOUT_6R112=1 P12_REFERENCE=$PWD/parity/realworld/assets/hs37d5.simple.fa \
//!   cargo test -p gatk-haplotypecaller --test holdout_6r112_p12_indel -- --nocapture
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::{
    build_per_haplotype_variation_events, merged_alleles_for_genotyping, overlapping_events,
    variation_events_at_position_from_cache, VariationEvent,
};
use gatk_haplotypecaller::genotyping::best_pl_index;
use gatk_haplotypecaller::hc_allele_mapping::{
    create_allele_mapper_with_events, replace_span_del_events, SPAN_DEL_ALLELE,
};
use gatk_haplotypecaller::hc_genotyping_engine::{
    java_emit_would_pass, l9_may_overwrite_pairhmm_gls_after_emit_fail, SparsePlShape,
    DEFAULT_STAND_EMIT_CONFIDENCE,
};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs, Haplotype,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
};
use serde_json::{json, Value};
use std::path::{Path, PathBuf};

const INTERVAL: &str = "2:92300000-92350000";
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 92_307_327;
const P12_STAND_EMIT: f64 = 10.0;
const JAVA_GT: &str = "1/1";
const JAVA_PL: &str = "45,3,0";
const JAVA_AD: &str = "0,1";
const JAVA_QUAL: &str = "35.44";

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
        "contig": e.contig,
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
        "p12_scope": e.contig == "2" || e.contig == "chr2",
    })
}

fn is_atg(e: &VariationEvent) -> bool {
    e.start_1based.get() == TARGET && e.ref_allele == "A" && e.alt_allele == "ATG"
}

/// Pre-6R.109 L9: P12 contig-2 scope false; SNP HomAltStrong; indel genome-wide support → true.
fn l9_pre_6r109(event: &VariationEvent, read_ref_ad: i32, read_alt_ad: i32) -> bool {
    if event.contig == "2" || event.contig == "chr2" {
        return false;
    }
    if event.is_snp() {
        SparsePlShape::pileup_is_hom_alt_strong(read_ref_ad, read_alt_ad)
    } else {
        read_alt_ad >= 2
    }
}

#[test]
fn holdout_6r112_p12_indel_dump() {
    if std::env::var("HOLDOUT_6R112").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R112=1");
        return;
    }
    let root = repo_root();
    let ref_fasta = std::env::var("P12_REFERENCE")
        .map(PathBuf::from)
        .unwrap_or_else(|_| root.join(REF_REL));
    let bam = root.join(BAM_REL);
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
        .expect("ActiveFull covering 2:92307327");
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
        "2",
    );

    let union_atg = outcome.assembly.variation_events().iter().any(is_atg);
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

    let mut n_hap_atg = 0usize;
    let mut hap_atg = Vec::new();
    for (i, h) in haps.iter().enumerate() {
        let has = hap_events.events_for(i).iter().any(is_atg);
        if has {
            n_hap_atg += 1;
            hap_atg.push(json!({
                "idx": i,
                "hash": fnv1a64_hex(&h.bases),
                "len": h.bases.len(),
                "is_reference": h.is_reference,
                "cigar": cigar_str(h),
            }));
        }
    }

    let event = VariationEvent::from_alleles("2", TARGET, "A", "ATG");
    let mapping = create_allele_mapper_with_events(
        &event,
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
            })
        })
        .collect();
    let atg_call = outcome.genotyped_calls.iter().find(|c| is_atg(&c.event));
    let calc_pl = atg_call.map(|c| c.genotype.format.pl_as_i32());
    let calc_ad = atg_call.map(|c| c.genotype.format.ad_as_i32());
    let calc_gls = atg_call
        .map(|c| c.genotype.genotype_log10_likelihoods.clone())
        .unwrap_or_default();
    let calculator_is_hom_ref = calc_pl.as_ref().is_some_and(|pl| {
        !pl.is_empty()
            && best_pl_index(
                &pl.iter()
                    .map(|&p| {
                        gatk_haplotypecaller::bio_ids::PhredLikelihood::from_i32_saturating(p)
                    })
                    .collect::<Vec<_>>(),
            ) == 0
    });
    let region_evs = outcome.assembly.variation_events();
    let java_emit_pass = atg_call.and_then(|c| {
        java_emit_would_pass(
            &c.event,
            &c.genotype.genotype_log10_likelihoods,
            &c.genotype.format,
            P12_STAND_EMIT,
            region_evs,
        )
        .ok()
    });
    let java_emit_pass_30 = atg_call.and_then(|c| {
        java_emit_would_pass(
            &c.event,
            &c.genotype.genotype_log10_likelihoods,
            &c.genotype.format,
            DEFAULT_STAND_EMIT_CONFIDENCE,
            region_evs,
        )
        .ok()
    });

    let (read_ref_ad, read_alt_ad) = calc_ad
        .as_ref()
        .map(|ad| {
            (
                ad.first().copied().unwrap_or(0),
                ad.get(1).copied().unwrap_or(0),
            )
        })
        .unwrap_or((0, 0));

    let l9_now = l9_may_overwrite_pairhmm_gls_after_emit_fail(
        &event,
        read_ref_ad,
        read_alt_ad,
        calculator_is_hom_ref,
    );
    let l9_pre = l9_pre_6r109(&event, read_ref_ad, read_alt_ad);

    let emitted = try_emit_call_region_variants(covering, &outcome, "SAMPLE", P12_STAND_EMIT)
        .unwrap_or_default();
    let emit_atg: Vec<Value> = emitted
        .iter()
        .filter(|r| r.position == TARGET)
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

    let doc = json!({
        "holdout": "2:92307327 A/ATG",
        "note_contig": "P12 L3a oracle is contig 2 (not 20)",
        "active_full": [covering.start.get(), covering.end.get()],
        "hap_count": haps.len(),
        "java_oracle": {
            "emits": true,
            "alleles": "A/ATG",
            "gt": JAVA_GT,
            "pl": JAVA_PL,
            "ad": JAVA_AD,
            "qual": JAVA_QUAL,
        },
        "inventory_a_eventmap": {
            "union_has_A_ATG": union_atg,
            "union_at_loc": union_at_loc,
            "overlap": overlap.iter().map(event_json).collect::<Vec<_>>(),
            "events_at_start_only": at_start.iter().map(event_json).collect::<Vec<_>>(),
            "events_with_spanning": at_spanning.iter().map(event_json).collect::<Vec<_>>(),
        },
        "inventory_b_haplotypes": {
            "n_with_A_ATG": n_hap_atg,
            "haps": hap_atg,
        },
        "inventory_c_d_merge": {
            "after_replace_span_dels": replaced.iter().map(event_json).collect::<Vec<_>>(),
            "merged_alleles_for_genotyping": merge_b,
        },
        "inventory_e_mapper": {
            "n_ref": mapping.ref_haplotype_indices.len(),
            "n_alt": mapping.alt_haplotype_indices.len(),
            "alt_idx": mapping.alt_haplotype_indices.iter().map(|i| i.get()).collect::<Vec<_>>(),
        },
        "read_ad": {
            "from_live_genotype": [read_ref_ad, read_alt_ad],
            "p12_scope": event.contig == "2" || event.contig == "chr2",
        },
        "live_genotype": live,
        "calculator": {
            "called": atg_call.is_some(),
            "pl": calc_pl,
            "ad": calc_ad,
            "gls": calc_gls,
            "calculator_is_hom_ref": calculator_is_hom_ref,
            "java_emit_stand10": java_emit_pass,
            "java_emit_stand30": java_emit_pass_30,
        },
        "l9_counterfactual": {
            "current_6r109": l9_now,
            "pre_6r109": l9_pre,
        },
        "live_emit": emit_atg,
        "rust_emits_A_ATG": emit_atg.iter().any(|r| r["alts"] == json!(["ATG"])),
    });
    println!("{}", serde_json::to_string_pretty(&doc).unwrap());

    assert!(
        union_atg,
        "Rust stored EventMap union still lists A/ATG (Java getVariationEvents also has it)"
    );
    assert!(
        !l9_now && !l9_pre,
        "P12 contig-2 L9 is off before and after 6R.109: now={l9_now} pre={l9_pre}"
    );
    // 6R.114 may restore A/ATG emit via CIGAR replay. 6R.112's L9-not-causal
    // conclusion is the holdout contract, not the pre-fix missing emit.
}
