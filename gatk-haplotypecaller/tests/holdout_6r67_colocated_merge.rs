//! 6R.67 forensic: why colocated SNP/indel merge does not fire at `20:29456344`
//! after the 6R.66 empty-EventMap pad-slice removal.
//!
//! Skipped unless `HOLDOUT_6R67=1`.
//!
//! ```text
//! HOLDOUT_6R67=1 cargo test -p gatk-haplotypecaller --test holdout_6r67_colocated_merge -- --nocapture
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::{
    build_per_haplotype_variation_events, is_colocated_snp_indel_merged_site,
    merged_alleles_for_genotyping, variation_events_at_position_from_cache, VariationEvent,
};
use gatk_haplotypecaller::hc_allele_mapping::{
    create_allele_mapper_with_events, hap_base_at_ref_locus, replace_span_del_events,
};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, take_colocated_merge_numerics,
    traverse_assembly_region_walker, try_emit_call_region_variants, AssemblyRegionCallDisposition,
    CallRegionArgs, GenomePosition, HaplotypeCallerEngine, HcGenotypingConfig, ReadFilterParams,
    WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use serde_json::json;
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const POS: u64 = 29_456_344;
const WINDOW_LO: u64 = 29_456_340;
const WINDOW_HI: u64 = 29_456_350;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn ev_key(e: &VariationEvent) -> String {
    format!("{} {}/{}", e.start_1based.get(), e.ref_allele, e.alt_allele)
}

fn is_snp_tc(e: &VariationEvent) -> bool {
    e.start_1based.get() == POS && e.ref_allele == "T" && e.alt_allele == "C"
}

fn is_del_tgt(e: &VariationEvent) -> bool {
    e.start_1based.get() == POS && e.ref_allele == "TG" && e.alt_allele == "T"
}

#[test]
fn holdout_6r67_colocated_merge_after_empty_eventmap() {
    if std::env::var("HOLDOUT_6R67").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R67=1");
        return;
    }
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    assert!(ref_fasta.is_file() && bam.is_file());

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
    let covering: Vec<_> = regions
        .iter()
        .filter(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= POS
                && r.end.get() >= POS
        })
        .collect();
    assert_eq!(covering.len(), 1);
    let region = covering[0];
    let args = CallRegionArgs::strict_java();
    let outcome = HaplotypeCallerEngine::call_region(region, &dict, &ref_fasta, &args)
        .expect("call")
        .expect("Some");

    let emitted =
        try_emit_call_region_variants(region, &outcome, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .unwrap_or_default();
    let nearby: Vec<_> = emitted
        .iter()
        .filter(|r| r.position >= WINDOW_LO && r.position <= WINDOW_HI)
        .map(|r| {
            json!({
                "pos": r.position,
                "ref": r.reference,
                "alt": r.alternate,
                "gt": r.samples.first().and_then(|s| s.gt.as_ref()).map(|g| g.alleles.clone()),
                "pl": r.samples.first().and_then(|s| s.pl.clone()),
            })
        })
        .collect();

    let apply_pad = outcome
        .assembly
        .haplotypes
        .iter()
        .find(|h| h.is_reference)
        .and_then(|h| h.genome_loc.map(|g| g.start_1based()))
        .unwrap_or_else(|| outcome.assembly.padded_reference_start_1based());
    let full_pad = outcome.assembly.padded_reference_start_1based();
    let apply_bytes = outcome.assembly.apply_bases_shared();
    let full_ref = outcome.assembly.reference_bases();
    let haps = &outcome.assembly.haplotypes;
    let contig = region.contig.as_str();
    let max_mnp = outcome.assembly.max_mnp_distance();

    let hap_events_apply = build_per_haplotype_variation_events(
        haps,
        apply_bytes.as_ref(),
        apply_pad,
        max_mnp,
        contig,
    );
    let hap_events_full =
        build_per_haplotype_variation_events(haps, full_ref, full_pad, max_mnp, contig);

    let emit_spanning = !HcGenotypingConfig::strict_java().disable_spanning_event_genotyping;
    let cache_at_loc =
        variation_events_at_position_from_cache(&hap_events_apply, POS, emit_spanning);
    let cache_full_at_loc =
        variation_events_at_position_from_cache(&hap_events_full, POS, emit_spanning);

    let mut hap_rows = Vec::new();
    let mut n_cache_events = 0usize;
    let mut nearest: Option<i64> = None;
    for (i, h) in haps.iter().enumerate() {
        let evs_apply = hap_events_apply.events_for(i);
        let evs = hap_events_full.events_for(i);
        n_cache_events += evs_apply.len();
        let overlapping: Vec<_> = evs
            .iter()
            .filter(|e| e.end_1based.get() >= POS && e.start_1based.get() <= POS)
            .map(ev_key)
            .collect();
        let window: Vec<_> = evs
            .iter()
            .filter(|e| e.start_1based.get() >= WINDOW_LO && e.start_1based.get() <= WINDOW_HI)
            .map(ev_key)
            .collect();
        for e in evs {
            let d = e.start_1based.get() as i64 - POS as i64;
            nearest = Some(match nearest {
                None => d,
                Some(cur) if d.abs() < cur.abs() => d,
                Some(cur) => cur,
            });
        }
        hap_rows.push(json!({
            "i": i,
            "is_ref": h.is_reference,
            "base_344": hap_base_at_ref_locus(h, apply_pad, POS).map(|b| (b as char).to_string()),
            "base_345": hap_base_at_ref_locus(h, apply_pad, POS + 1).map(|b| (b as char).to_string()),
            "n_events": evs.len(),
            "overlap_at_loc": overlapping,
            "window": window,
            "has_snp_tc": evs.iter().any(is_snp_tc),
            "has_del_tgt": evs.iter().any(is_del_tgt),
        }));
    }

    let stored: Vec<_> = outcome
        .assembly
        .variation_events()
        .iter()
        .filter(|e| e.end_1based.get() >= POS && e.start_1based.get() <= POS)
        .cloned()
        .collect();
    let stored_window: Vec<_> = outcome
        .assembly
        .variation_events()
        .iter()
        .filter(|e| e.start_1based.get() >= WINDOW_LO && e.start_1based.get() <= WINDOW_HI)
        .map(ev_key)
        .collect();

    let loc_pos = GenomePosition::new_1based(POS);
    // Production 6R.67 walk: genotyping EventMap cache is the full padded ref.
    let mut raw = cache_full_at_loc.clone();
    let mut seen: std::collections::HashSet<(u64, String, String)> = raw
        .iter()
        .map(|e| {
            (
                e.start_1based.get(),
                e.ref_allele.clone(),
                e.alt_allele.clone(),
            )
        })
        .collect();
    for e in outcome.assembly.variation_events() {
        let overlaps = e.end_1based >= loc_pos && e.start_1based <= loc_pos;
        if !overlaps {
            continue;
        }
        if !emit_spanning && e.start_1based != loc_pos {
            continue;
        }
        let key = (
            e.start_1based.get(),
            e.ref_allele.clone(),
            e.alt_allele.clone(),
        );
        if seen.insert(key) {
            raw.push(e.clone());
        }
    }
    let at_loc = replace_span_del_events(&raw, POS, apply_pad, apply_bytes.as_ref());
    let merged = merged_alleles_for_genotyping(&at_loc, POS);
    let colocated = merged
        .as_ref()
        .map(|(r, a)| is_colocated_snp_indel_merged_site(r, a))
        .unwrap_or(false);

    let mut pool_sizes: Option<Vec<usize>> = None;
    if let Some((long_ref, alts)) = merged.as_ref() {
        let mut sizes = Vec::new();
        let end_1based = POS.saturating_add(long_ref.len().saturating_sub(1) as u64);
        for alt in alts {
            let ev = VariationEvent {
                contig: contig.to_string(),
                start_1based: GenomePosition::new_1based(POS),
                end_1based: GenomePosition::new_1based(end_1based),
                ref_allele: long_ref.clone(),
                alt_allele: alt.clone(),
            };
            let mapping = create_allele_mapper_with_events(
                &ev,
                POS,
                haps,
                apply_pad,
                apply_bytes.as_ref(),
                max_mnp,
                emit_spanning,
                Some(&hap_events_full),
            );
            sizes.push(mapping.alt_haplotype_indices.len());
        }
        pool_sizes = Some(sizes);
    }

    let live = take_colocated_merge_numerics();
    let merge_fired = live.iter().any(|n| n.loc == POS);

    let skip = if merged.is_none() {
        "merged_alleles_for_genotyping_none"
    } else if !colocated {
        "not_colocated_snp_indel_merged_site"
    } else if pool_sizes
        .as_ref()
        .is_some_and(|s| s.iter().all(|&n| n == 0))
    {
        "empty_alt_pools_not_applicable"
    } else if merge_fired {
        "merge_fired"
    } else {
        "merge_eligible_but_no_live_numerics"
    };

    let doc = json!({
        "locus": "20:29456344",
        "pads": {
            "apply_pad": apply_pad,
            "full_pad": full_pad,
            "apply_eq_full": apply_pad == full_pad,
        },
        "n_haps": haps.len(),
        "n_cache_events_apply": n_cache_events,
        "nearest_event_start_delta": nearest,
        "hap_cache_at_loc": cache_at_loc.iter().map(ev_key).collect::<Vec<_>>(),
        "hap_cache_fullpad_at_loc": cache_full_at_loc.iter().map(ev_key).collect::<Vec<_>>(),
        "stored_at_loc": stored.iter().map(ev_key).collect::<Vec<_>>(),
        "stored_window": stored_window,
        "union_at_loc": at_loc.iter().map(ev_key).collect::<Vec<_>>(),
        "union_has_snp_tc": at_loc.iter().any(is_snp_tc),
        "union_has_del_tgt": at_loc.iter().any(is_del_tgt),
        "any_hap_has_snp_tc": hap_rows.iter().any(|r| r["has_snp_tc"] == true),
        "any_hap_has_del_tgt": hap_rows.iter().any(|r| r["has_del_tgt"] == true),
        "merged_alleles": merged.as_ref().map(|(r, a)| json!({"long_ref": r, "alts": a})),
        "colocated": colocated,
        "mapper_alt_pool_sizes": pool_sizes,
        "skip": skip,
        "merge_fired": merge_fired,
        "vcf_nearby": nearby,
        "haps": hap_rows,
    });
    eprintln!("{}", serde_json::to_string_pretty(&doc).expect("json"));

    let vcf_tc = emitted
        .iter()
        .any(|r| r.position == POS && r.reference == "T" && r.alternate.iter().any(|a| a == "C"));
    assert!(
        cache_full_at_loc.iter().any(is_snp_tc) && cache_full_at_loc.iter().any(is_del_tgt),
        "full-pad EventMaps must carry T/C and TG/T: {cache_full_at_loc:?}"
    );
    assert!(
        merge_fired,
        "6R.67: colocated merge must run before genotype; skip={skip} merged={merged:?} union={:?}",
        at_loc.iter().map(ev_key).collect::<Vec<_>>()
    );
    assert!(
        vcf_tc,
        "canonical representation after merge+unused-ALT+reverseTrim must be T/C; nearby={nearby:?}"
    );
}
