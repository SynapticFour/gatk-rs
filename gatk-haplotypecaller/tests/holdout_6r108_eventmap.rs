//! 6R.108 forensic dump: HOLDOUT_6R53 remainder `20:29455644 A/AT`.
//!
//! Inventory A (EventMap union) before GLs / L9. Skipped unless `HOLDOUT_6R108=1`.
//! Production change: NONE. Does not reopen 6R.107.
//!
//! ```text
//! HOLDOUT_6R108=1 cargo test -p gatk-haplotypecaller --test holdout_6r108_eventmap -- --nocapture
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::bio_ids::HaplotypeIndex;
use gatk_haplotypecaller::event_map::{
    build_per_haplotype_variation_events, merged_alleles_for_genotyping, overlapping_events,
    variation_events_at_position_from_cache, VariationEvent,
};
use gatk_haplotypecaller::genotyping::{
    best_pl_index, diploid_genotype_alleles_from_pl_index, emit_genotype_format_fields,
};
use gatk_haplotypecaller::hc_allele_mapping::{
    create_allele_mapper_with_events, replace_span_del_events, SPAN_DEL_ALLELE,
};
use gatk_haplotypecaller::hc_genotyping_engine::{
    biallelic_genotype_log10_likelihoods_gatk, java_alignment_read_overlaps_interval,
    l9_may_overwrite_pairhmm_gls_after_emit_fail, marginalize_rows_to_biallelic_alleles,
    region_likelihoods_to_rows, SparsePlShape, DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN,
    DEFAULT_STAND_EMIT_CONFIDENCE,
};
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
const TARGET: u64 = 29_455_644;
const ADJACENT: u64 = 29_455_649;
const WINDOW: (u64, u64) = (29_455_620, 29_455_680);
const LOG10_GLOBAL_READ_MISMATCHING_RATE: f64 = -4.5;

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

fn is_aat(e: &VariationEvent) -> bool {
    e.start_1based.get() == TARGET && e.ref_allele == "A" && e.alt_allele == "AT"
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

fn apply_mismapping_floor(marg: &mut [gatk_haplotypecaller::genotyping::ReadLikelihoodRow]) {
    for row in marg {
        let lr = row.haplotype_log10_likelihoods[0];
        let la = row.haplotype_log10_likelihoods[1];
        if !lr.is_finite() || !la.is_finite() {
            continue;
        }
        let best = lr.max(la);
        let floor = best + LOG10_GLOBAL_READ_MISMATCHING_RATE;
        if lr < floor {
            row.haplotype_log10_likelihoods[0] = floor;
        }
        if la < floor {
            row.haplotype_log10_likelihoods[1] = floor;
        }
    }
}

fn gl_pl_from_pools(
    likelihoods: &[gatk_haplotypecaller::RegionReadLikelihood],
    n_haps: usize,
    ref_p: &[HaplotypeIndex],
    alt_p: &[HaplotypeIndex],
) -> (Vec<f64>, Vec<i32>, Vec<i32>) {
    let rows = region_likelihoods_to_rows(likelihoods, n_haps);
    let mut marg = marginalize_rows_to_biallelic_alleles(&rows, ref_p, alt_p);
    apply_mismapping_floor(&mut marg);
    let gls = biallelic_genotype_log10_likelihoods_gatk(&marg, 0, 1);
    let fmt = emit_genotype_format_fields(&gls, &[0, 0]).expect("fmt");
    let gt = diploid_genotype_alleles_from_pl_index(2, best_pl_index(&fmt.pl));
    (gls, fmt.pl_as_i32(), gt)
}

fn filter_overlapping_likelihoods(
    all: &[gatk_haplotypecaller::RegionReadLikelihood],
    reads: &[gatk_haplotypecaller::shared_bam::SharedBamRecord],
) -> Vec<gatk_haplotypecaller::RegionReadLikelihood> {
    all.iter()
        .filter(|row| {
            reads.get(row.read_index.get()).is_some_and(|r| {
                java_alignment_read_overlaps_interval(
                    r.as_ref(),
                    TARGET,
                    TARGET,
                    DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN,
                )
            })
        })
        .cloned()
        .collect()
}

#[test]
fn holdout_6r108_eventmap_aat_dump() {
    if std::env::var("HOLDOUT_6R108").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R108=1");
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
    let region = regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() == REGION.0
                && r.end.get() == REGION.1
        })
        .expect("ActiveFull 20:29455560-29455744");
    let args = CallRegionArgs::strict_java();
    let outcome = HaplotypeCallerEngine::call_region(region, &dict, &ref_fasta, &args)
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

    let union: Vec<Value> = outcome
        .assembly
        .variation_events()
        .iter()
        .filter(|e| {
            let s = e.start_1based.get();
            s >= WINDOW.0 && s <= WINDOW.1
                || e.start_1based.get() <= TARGET && e.end_1based.get() >= TARGET
        })
        .map(event_json)
        .collect();
    let union_has_aat = outcome.assembly.variation_events().iter().any(is_aat);
    let union_at_644: Vec<Value> = outcome
        .assembly
        .variation_events()
        .iter()
        .filter(|e| e.start_1based.get() == TARGET)
        .map(event_json)
        .collect();
    let union_at_649: Vec<Value> = outcome
        .assembly
        .variation_events()
        .iter()
        .filter(|e| e.start_1based.get() == ADJACENT)
        .map(event_json)
        .collect();

    let mut per_hap = Vec::new();
    let mut n_hap_aat = 0usize;
    for (i, h) in haps.iter().enumerate() {
        let evs: Vec<Value> = hap_events
            .events_for(i)
            .iter()
            .filter(|e| {
                e.start_1based.get() == TARGET
                    || e.start_1based.get() <= TARGET && e.end_1based.get() >= TARGET
            })
            .map(event_json)
            .collect();
        let has_aat = hap_events.events_for(i).iter().any(is_aat);
        if has_aat {
            n_hap_aat += 1;
        }
        per_hap.push(json!({
            "idx": i,
            "hash": fnv1a64_hex(&h.bases),
            "len": h.bases.len(),
            "is_reference": h.is_reference,
            "score": h.score,
            "cigar": cigar_str(h),
            "align_start": h.alignment_start_hap_wrt_ref,
            "eventmap_aat": has_aat,
            "events_at_or_span_644": evs,
        }));
    }

    let at_start = variation_events_at_position_from_cache(&hap_events, TARGET, false);
    let at_spanning = variation_events_at_position_from_cache(&hap_events, TARGET, true);
    let overlap = overlapping_events(outcome.assembly.variation_events(), TARGET);
    let replaced = replace_span_del_events(&at_spanning, TARGET, pad, ref_bytes);
    let merge_b = merged_alleles_for_genotyping(&replaced, TARGET);
    let merge_from_overlap = merged_alleles_for_genotyping(
        &replace_span_del_events(&overlap, TARGET, pad, ref_bytes),
        TARGET,
    );

    let event = VariationEvent::from_alleles("20", TARGET, "A", "AT");
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

    let ll_all = &outcome.read_likelihoods;
    let ll_overlap = filter_overlapping_likelihoods(ll_all, &outcome.genotyping_reads);
    let (calc_gl, calc_pl, calc_gt) =
        if mapping.ref_haplotype_indices.is_empty() && mapping.alt_haplotype_indices.is_empty() {
            (vec![], vec![], vec![])
        } else {
            gl_pl_from_pools(
                &ll_overlap,
                haps.len(),
                &mapping.ref_haplotype_indices,
                &mapping.alt_haplotype_indices,
            )
        };

    let live_ad = outcome
        .genotyped_calls
        .iter()
        .find(|c| is_aat(&c.event))
        .map(|c| c.genotype.format.ad_as_i32());
    let (pileup_ref, pileup_alt) = match live_ad.as_deref() {
        Some([r, a, ..]) => (*r, *a),
        _ => (i32::MIN, i32::MIN),
    };
    let l9 = if pileup_ref == i32::MIN {
        None
    } else {
        Some(l9_may_overwrite_pairhmm_gls_after_emit_fail(
            &event,
            pileup_ref,
            pileup_alt,
            calc_gt.as_slice() == [0, 0],
        ))
    };
    let shape = if pileup_ref == i32::MIN {
        None
    } else {
        Some(format!(
            "{:?}",
            SparsePlShape::from_pileup_depths(pileup_ref, pileup_alt)
        ))
    };

    let calls: Vec<Value> = outcome
        .genotyped_calls
        .iter()
        .filter(|c| c.event.start_1based.get() == TARGET)
        .map(|c| {
            json!({
                "start": c.event.start_1based.get(),
                "end": c.event.end_1based.get(),
                "ref": c.event.ref_allele,
                "alt": c.event.alt_allele,
                "extra_alts": c.extra_alt_alleles,
                "ad": c.genotype.format.ad_as_i32(),
                "pl": c.genotype.format.pl_as_i32(),
                "gq": c.genotype.format.gq.as_i32(),
                "gl": c.genotype.genotype_log10_likelihoods,
            })
        })
        .collect();

    let emitted =
        try_emit_call_region_variants(region, &outcome, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .unwrap_or_default();
    let emit_c: Vec<Value> = emitted
        .iter()
        .filter(|r| r.position == TARGET || r.position == ADJACENT)
        .map(|r| {
            json!({
                "pos": r.position,
                "ref": r.reference,
                "alts": r.alternate,
            })
        })
        .collect();

    let per_hap_event_set: BTreeSet<(u64, String, String)> = (0..haps.len())
        .flat_map(|i| hap_events.events_for(i).iter())
        .filter(|e| e.start_1based.get() == TARGET)
        .map(|e| {
            (
                e.start_1based.get(),
                e.ref_allele.clone(),
                e.alt_allele.clone(),
            )
        })
        .collect();
    let union_event_set: BTreeSet<(u64, String, String)> = outcome
        .assembly
        .variation_events()
        .iter()
        .filter(|e| e.start_1based.get() == TARGET)
        .map(|e| {
            (
                e.start_1based.get(),
                e.ref_allele.clone(),
                e.alt_allele.clone(),
            )
        })
        .collect();
    let supplemental_only: Vec<_> = union_event_set
        .difference(&per_hap_event_set)
        .cloned()
        .collect();

    let doc = json!({
        "holdout": "20:29455644 A/AT",
        "active_full": REGION,
        "hap_count": haps.len(),
        "pad": pad,
        "frozen_vcf": {
            "java_644": vcf_record(&java_vcf, TARGET),
            "rust_644": vcf_record(&rust_vcf, TARGET),
            "java_649": vcf_record(&java_vcf, ADJACENT),
            "rust_649": vcf_record(&rust_vcf, ADJACENT),
        },
        "inventory_a_eventmap": {
            "union_window_or_span": union,
            "union_at_644": union_at_644,
            "union_has_A_AT": union_has_aat,
            "per_hap_n_with_A_AT": n_hap_aat,
            "per_hap_unique_at_644": per_hap_event_set.iter().map(|(s,r,a)| format!("{s}:{r}/{a}")).collect::<Vec<_>>(),
            "union_not_in_any_hap_eventmap": supplemental_only.iter().map(|(s,r,a)| format!("{s}:{r}/{a}")).collect::<Vec<_>>(),
            "events_at_start_only_n": at_start.len(),
            "events_at_start_only": at_start.iter().map(event_json).collect::<Vec<_>>(),
            "events_with_spanning_n": at_spanning.len(),
            "events_with_spanning": at_spanning.iter().map(event_json).collect::<Vec<_>>(),
            "overlap_union_n": overlap.len(),
            "overlap_union": overlap.iter().map(event_json).collect::<Vec<_>>(),
        },
        "inventory_a_context_649": {
            "union_at_649": union_at_649,
        },
        "inventory_b_merged": {
            "after_replace_span_dels": replaced.iter().map(event_json).collect::<Vec<_>>(),
            "merged_alleles_for_genotyping": merge_b,
            "merged_from_union_overlap": merge_from_overlap,
        },
        "allele_mapper": {
            "ref": mapping.ref_allele,
            "alt": mapping.alt_allele,
            "n_ref": mapping.ref_haplotype_indices.len(),
            "n_alt": mapping.alt_haplotype_indices.len(),
            "ref_idx": mapping.ref_haplotype_indices.iter().map(|i| i.get()).collect::<Vec<_>>(),
            "alt_idx": mapping.alt_haplotype_indices.iter().map(|i| i.get()).collect::<Vec<_>>(),
        },
        "calculate_gls_from_mapper": {
            "log10_gl": calc_gl,
            "pl": calc_pl,
            "gt_from_pl": calc_gt,
        },
        "l9_read_only": {
            "note": "AD is live genotyped FORMAT AD, not a new pileup walk",
            "live_ad": live_ad,
            "sparse_from_live_ad": shape,
            "l9_may_overwrite_pairhmm_after_emit_fail": l9,
            "event_is_indel": event.is_indel(),
            "event_is_snp": event.is_snp(),
        },
        "live_genotype_at_644": calls,
        "live_emit_644_649": emit_c,
        "per_hap": per_hap,
    });
    println!("{}", serde_json::to_string_pretty(&doc).unwrap());

    assert!(
        vcf_record(&java_vcf, TARGET).is_none(),
        "frozen Java VCF must not contain 20:29455644"
    );
    let rust_frozen = vcf_record(&rust_vcf, TARGET).expect("frozen rust VCF has 29455644");
    assert_eq!(rust_frozen["alleles"], "A/AT");
}
