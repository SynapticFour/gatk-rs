//! 6R.115 forensic dump: first post-CIGAR divergence at `2:92307327`.
//! EventMap union → merged alleles → mapper pools → live hap population.
//! Likelihoods/GLs only if A–D match (printed, not investigated here).
//! Skipped unless `HOLDOUT_6R115=1`. Does not reopen 6R.107–6R.114.
//!
//! ```text
//! HOLDOUT_6R115=1 P12_REFERENCE=$PWD/parity/realworld/assets/hs37d5.simple.fa \
//!   cargo test -p gatk-haplotypecaller --test holdout_6r115_post_cigar -- --nocapture
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::{
    build_per_haplotype_variation_events, collect_variation_events,
    is_colocated_snp_indel_merged_site, is_same_ref_multi_indel_merged_site,
    merged_alleles_for_genotyping, merged_biallelic_sites_at_position, merged_site_uses_joint_gls,
    overlapping_events, variation_events_at_position_from_cache, VariationEvent,
};
use gatk_haplotypecaller::hc_allele_mapping::{
    create_allele_mapper_with_events, hap_base_at_ref_locus, replace_span_del_events,
    SPAN_DEL_ALLELE,
};
use gatk_haplotypecaller::hc_genotyping_engine::DEFAULT_STAND_EMIT_CONFIDENCE;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs, Haplotype,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
};
use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

const INTERVAL: &str = "2:92300000-92350000";
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 92_307_327;
const P12_STAND_EMIT: f64 = 10.0;
const MOTIF_LO: u64 = 92_307_324;
const MOTIF_HI: u64 = 92_307_330;

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

fn is_atg(e: &VariationEvent) -> bool {
    e.start_1based.get() == TARGET && e.ref_allele == "A" && e.alt_allele == "ATG"
}

fn is_ag(e: &VariationEvent) -> bool {
    e.start_1based.get() == TARGET && e.ref_allele == "A" && e.alt_allele == "G"
}

fn linear_motif(h: &Haplotype, pad: u64) -> String {
    let mut s = String::new();
    for loc in MOTIF_LO..=MOTIF_HI {
        match hap_base_at_ref_locus(h, pad, loc) {
            Some(b) => s.push(b as char),
            None => s.push('.'),
        }
    }
    s
}

fn mapper_pool(mapping: &gatk_haplotypecaller::hc_allele_mapping::AlleleHaplotypeMapping) -> Value {
    json!({
        "ref_allele": mapping.ref_allele,
        "alt_allele": mapping.alt_allele,
        "ref_idx": mapping
            .ref_haplotype_indices
            .iter()
            .map(|i| i.get())
            .collect::<Vec<_>>(),
        "alt_idx": mapping
            .alt_haplotype_indices
            .iter()
            .map(|i| i.get())
            .collect::<Vec<_>>(),
        "n_ref": mapping.ref_haplotype_indices.len(),
        "n_alt": mapping.alt_haplotype_indices.len(),
    })
}

#[test]
fn holdout_6r115_post_cigar_dump() {
    if std::env::var("HOLDOUT_6R115").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R115=1");
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
    let full_pad = outcome.assembly.padded_reference_start_1based();
    let full_ref = outcome.assembly.reference_bases();
    let apply_bases = outcome.assembly.apply_bases_shared();
    let apply_pad = outcome
        .assembly
        .haplotypes
        .iter()
        .find(|h| h.is_reference)
        .and_then(|h| h.genome_loc)
        .map(|g| g.start_1based())
        .unwrap_or(full_pad);
    let haps = &outcome.assembly.haplotypes;
    let max_mnp = outcome.assembly.max_mnp_distance();
    let hap_events = build_per_haplotype_variation_events(haps, full_ref, full_pad, max_mnp, "2");
    let cigar_union = collect_variation_events(haps, full_ref, full_pad, "2", max_mnp);

    let stored_at: Vec<Value> = outcome
        .assembly
        .variation_events()
        .iter()
        .filter(|e| e.start_1based.get() == TARGET)
        .map(event_json)
        .collect();
    let cigar_at: Vec<Value> = cigar_union
        .iter()
        .filter(|e| e.start_1based.get() == TARGET)
        .map(event_json)
        .collect();

    let mut n_atg = 0usize;
    let mut n_ag = 0usize;
    let mut n_both = 0usize;
    let mut n_neither = 0usize;
    let mut hashes = BTreeSet::new();
    let mut per_hap = Vec::new();
    let mut atg_hashes = BTreeSet::new();
    let mut ag_hashes = BTreeSet::new();
    for (i, h) in haps.iter().enumerate() {
        let evs: Vec<_> = hap_events
            .events_for(i)
            .iter()
            .filter(|e| e.start_1based.get() == TARGET)
            .cloned()
            .collect();
        let has_atg = evs.iter().any(is_atg);
        let has_ag = evs.iter().any(is_ag);
        if has_atg {
            n_atg += 1;
            atg_hashes.insert(fnv1a64_hex(&h.bases));
        }
        if has_ag {
            n_ag += 1;
            ag_hashes.insert(fnv1a64_hex(&h.bases));
        }
        if has_atg && has_ag {
            n_both += 1;
        }
        if !has_atg && !has_ag {
            n_neither += 1;
        }
        hashes.insert(fnv1a64_hex(&h.bases));
        per_hap.push(json!({
            "idx": i,
            "hash": fnv1a64_hex(&h.bases),
            "len": h.bases.len(),
            "is_reference": h.is_reference,
            "score": h.score,
            "cigar": cigar_str(h),
            "align_start": h.alignment_start_hap_wrt_ref,
            "genome_loc": h.genome_loc.map(|g| format!("{}-{}", g.start_1based(), g.end_1based())),
            "linear_motif_92307324_330": linear_motif(h, full_pad),
            "has_A_ATG": has_atg,
            "has_A_G": has_ag,
            "events_at_327": evs.iter().map(event_json).collect::<Vec<_>>(),
        }));
    }

    let at_start = variation_events_at_position_from_cache(&hap_events, TARGET, false);
    let at_spanning = variation_events_at_position_from_cache(&hap_events, TARGET, true);
    let overlap_stored = overlapping_events(outcome.assembly.variation_events(), TARGET);
    let mut raw = at_spanning.clone();
    {
        let mut seen: BTreeSet<(u64, String, String)> = raw
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
            let overlaps = e.end_1based.get() >= TARGET && e.start_1based.get() <= TARGET;
            if !overlaps {
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
    }
    let replaced = replace_span_del_events(&raw, TARGET, apply_pad, apply_bases.as_ref());
    let merge_b = merged_alleles_for_genotyping(&replaced, TARGET);
    let biallelics = merged_biallelic_sites_at_position(&replaced, TARGET);
    let (joint_gls, colocated, same_ref_multi) = match &merge_b {
        Some((long_ref, alts)) => (
            merged_site_uses_joint_gls(&replaced, TARGET, long_ref, alts),
            is_colocated_snp_indel_merged_site(long_ref, alts),
            is_same_ref_multi_indel_merged_site(&replaced, TARGET, long_ref, alts),
        ),
        None => (false, false, false),
    };

    let map_atg = create_allele_mapper_with_events(
        &VariationEvent::from_alleles("2", TARGET, "A", "ATG"),
        TARGET,
        haps,
        apply_pad,
        apply_bases.as_ref(),
        max_mnp,
        true,
        Some(&hap_events),
    );
    let map_g = create_allele_mapper_with_events(
        &VariationEvent::from_alleles("2", TARGET, "A", "G"),
        TARGET,
        haps,
        apply_pad,
        apply_bases.as_ref(),
        max_mnp,
        true,
        Some(&hap_events),
    );
    let map_a_ref = create_allele_mapper_with_events(
        &VariationEvent::from_alleles("2", TARGET, "A", "A"),
        TARGET,
        haps,
        apply_pad,
        apply_bases.as_ref(),
        max_mnp,
        true,
        Some(&hap_events),
    );

    let mut assignment: BTreeMap<usize, Vec<String>> = BTreeMap::new();
    for i in map_atg.alt_haplotype_indices.iter().map(|x| x.get()) {
        assignment.entry(i).or_default().push("ATG".to_string());
    }
    for i in map_g.alt_haplotype_indices.iter().map(|x| x.get()) {
        assignment.entry(i).or_default().push("G".to_string());
    }
    for i in map_atg.ref_haplotype_indices.iter().map(|x| x.get()) {
        assignment
            .entry(i)
            .or_default()
            .push("A_via_ATG_mapper".to_string());
    }
    for i in map_g.ref_haplotype_indices.iter().map(|x| x.get()) {
        assignment
            .entry(i)
            .or_default()
            .push("A_via_G_mapper".to_string());
    }
    let dual: Vec<usize> = assignment
        .iter()
        .filter(|(_, v)| v.iter().any(|s| s == "ATG") && v.iter().any(|s| s == "G"))
        .map(|(k, _)| *k)
        .collect();

    let live: Vec<Value> = outcome
        .genotyped_calls
        .iter()
        .filter(|c| c.event.start_1based.get() == TARGET)
        .map(|c| {
            json!({
                "ref": c.event.ref_allele,
                "alt": c.event.alt_allele,
                "extra_alts": c.extra_alt_alleles,
                "post_merge_unused_alt_subset": c.post_merge_unused_alt_subset,
                "ad": c.genotype.format.ad_as_i32(),
                "pl": c.genotype.format.pl_as_i32(),
                "gls": c.genotype.genotype_log10_likelihoods,
            })
        })
        .collect();
    let emitted = try_emit_call_region_variants(covering, &outcome, "SAMPLE", P12_STAND_EMIT)
        .unwrap_or_default();
    let emit_c: Vec<Value> = emitted
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
    let emit_30 =
        try_emit_call_region_variants(covering, &outcome, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .unwrap_or_default();
    let emit_30_c: Vec<Value> = emit_30
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

    let unique_alts: BTreeSet<String> = replaced
        .iter()
        .filter(|e| e.start_1based.get() == TARGET)
        .map(|e| format!("{}/{}", e.ref_allele, e.alt_allele))
        .collect();

    let dump = json!({
        "mission": "6R.115",
        "target": "2:92307327",
        "region": format!("{}:{}-{}", covering.contig, covering.start.get(), covering.end.get()),
        "full_pad": full_pad,
        "apply_pad": apply_pad,
        "full_ref_len": full_ref.len(),
        "apply_ref_len": apply_bases.len(),
        "hap_count": haps.len(),
        "max_mnp_distance": max_mnp,
        "A_EVENTMAP": {
            "stored_union_at_loc": stored_at,
            "cigar_replay_union_at_loc": cigar_at,
            "hap_eventmap_at_start": at_start.iter().map(event_json).collect::<Vec<_>>(),
            "hap_eventmap_spanning": at_spanning.iter().map(event_json).collect::<Vec<_>>(),
            "stored_overlapping": overlap_stored.iter().map(event_json).collect::<Vec<_>>(),
            "raw_hap_plus_stored": raw.iter().map(event_json).collect::<Vec<_>>(),
            "replaced_span_dels": replaced.iter().map(event_json).collect::<Vec<_>>(),
            "unique_alleles": unique_alts.iter().cloned().collect::<Vec<_>>(),
            "n_haps_A_ATG": n_atg,
            "n_haps_A_G": n_ag,
            "n_haps_both": n_both,
            "n_haps_neither": n_neither,
            "atg_hashes": atg_hashes.iter().cloned().collect::<Vec<_>>(),
            "ag_hashes": ag_hashes.iter().cloned().collect::<Vec<_>>(),
        },
        "B_MERGED": {
            "merged_alleles_for_genotyping": merge_b.as_ref().map(|(r, a)| json!({
                "long_ref": r,
                "alts": a,
            })),
            "merged_site_uses_joint_gls": joint_gls,
            "is_colocated_snp_indel": colocated,
            "is_same_ref_multi_indel": same_ref_multi,
            "biallelic_sites": biallelics.iter().map(event_json).collect::<Vec<_>>(),
            "A_G_independently_retained": unique_alts.contains("A/G"),
            "A_ATG_retained": unique_alts.contains("A/ATG"),
        },
        "C_MAPPER": {
            "A_ATG": mapper_pool(&map_atg),
            "A_G": mapper_pool(&map_g),
            "dummy_A_A": mapper_pool(&map_a_ref),
            "dual_ATG_and_G_hap_idx": dual,
        },
        "D_LIVE_HAPS": {
            "unique_seq_hashes": hashes.iter().cloned().collect::<Vec<_>>(),
            "per_haplotype": per_hap,
            "assignment_by_idx": assignment.iter().map(|(k, v)| json!({
                "idx": k,
                "alleles": v,
            })).collect::<Vec<_>>(),
        },
        "LIVE_CALL": {
            "genotyped": live,
            "emit_stand10": emit_c,
            "emit_stand30": emit_30_c,
        },
    });
    println!("{}", serde_json::to_string_pretty(&dump).expect("json"));
}
