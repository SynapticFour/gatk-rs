//! 6R.64 forensic: numerical AD/PL/QUAL at `20:29456344 T/C` after representation match.
//!
//! No production algorithm change. Skipped unless `HOLDOUT_6R64=1`.
//!
//! ```text
//! HOLDOUT_6R64=1 cargo test -p gatk-haplotypecaller --test holdout_6r64_qual_ad_pl_forensic -- --nocapture
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::{
    build_per_haplotype_variation_events, merged_alleles_for_genotyping,
    variation_events_at_position_from_cache,
};
use gatk_haplotypecaller::hc_allele_mapping::replace_span_del_events;
use gatk_haplotypecaller::variant_site_hc_annotations::qual_from_af_calculation;
use gatk_haplotypecaller::{
    audit_colocated_snp_indel_merge_numerics, call_disposition, flatten_assembly_regions,
    take_colocated_merge_numerics, traverse_assembly_region_walker, try_emit_call_region_variants,
    AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition, HaplotypeCallerEngine,
    HcGenotypingConfig, ReadFilterParams, WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use serde_json::json;
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const POS_SNP: u64 = 29_456_344;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

#[test]
fn holdout_6r64_qual_ad_pl_forensic_29456344() {
    if std::env::var("HOLDOUT_6R64").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R64=1");
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
            ) && r.start.get() <= POS_SNP
                && r.end.get() >= POS_SNP
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
    let vcf = emitted.iter().find(|r| {
        r.position == POS_SNP && r.reference == "T" && r.alternate.iter().any(|a| a == "C")
    });
    let vcf = vcf.expect("6R.63 representation T/C must remain");
    let vcf_pl = vcf.samples[0].pl.clone().unwrap_or_default();
    let vcf_ad = vcf.samples[0].ad.clone().unwrap_or_default();
    let vcf_qual = vcf.quality;

    let pad = outcome
        .assembly
        .haplotypes
        .iter()
        .find(|h| h.is_reference)
        .and_then(|h| h.genome_loc.map(|g| g.start_1based()))
        .unwrap_or_else(|| outcome.assembly.padded_reference_start_1based());
    let ref_bytes = outcome.assembly.apply_bases_shared();
    let hap_events = build_per_haplotype_variation_events(
        &outcome.assembly.haplotypes,
        ref_bytes.as_ref(),
        pad,
        outcome.assembly.max_mnp_distance(),
        &region.contig,
    );
    let emit_spanning = !HcGenotypingConfig::strict_java().disable_spanning_event_genotyping;
    let mut raw = variation_events_at_position_from_cache(&hap_events, POS_SNP, emit_spanning);
    let hap_cache_events: Vec<String> = raw
        .iter()
        .map(|e| format!("{}/{}", e.ref_allele, e.alt_allele))
        .collect();
    // Same union as production `merge_stored_variation_events_at_position`.
    let loc_pos = GenomePosition::new_1based(POS_SNP);
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
    let at_loc = replace_span_del_events(&raw, POS_SNP, pad, ref_bytes.as_ref());
    let merged = merged_alleles_for_genotyping(&at_loc, POS_SNP);
    let site_call = outcome.genotyped_calls.iter().find(|c| {
        c.event.start_1based.get() == POS_SNP
            && (c.event.ref_allele == "T" || c.event.ref_allele == "TG")
    });
    let live = take_colocated_merge_numerics();
    let numerics = live
        .iter()
        .find(|n| n.loc == POS_SNP)
        .cloned()
        .or_else(|| {
            audit_colocated_snp_indel_merge_numerics(
                &at_loc,
                POS_SNP,
                &outcome.read_likelihoods,
                &outcome.genotyping_reads,
                &outcome.assembly.haplotypes,
                ref_bytes.as_ref(),
                pad,
                region.start.get(),
                region.end.get(),
                outcome.assembly.max_mnp_distance(),
                &HcGenotypingConfig::strict_java(),
                Some(&hap_events),
            )
            .ok()
            .flatten()
        })
        .unwrap_or_else(|| {
            let events: Vec<_> = at_loc
                .iter()
                .map(|e| format!("{} {}/{}", e.start_1based.get(), e.ref_allele, e.alt_allele))
                .collect();
            panic!("colocated merge numerics missing: live={live:?} merged={merged:?} events={events:?} site={site_call:?}");
        });

    let rust_gl = [-29.8, 0.0, -110.3];
    let java_gl = [-54.2, 0.0, -135.3];
    let qual_from_rust_pl = qual_from_af_calculation(&rust_gl).ok();
    let qual_from_java_pl = qual_from_af_calculation(&java_gl).ok();

    let doc = json!({
        "locus": "20:29456344 T/C",
        "vcf": {
            "ref": vcf.reference,
            "alt": vcf.alternate,
            "gt": vcf.samples[0].gt.as_ref().map(|g| g.alleles.clone()),
            "ad": vcf_ad,
            "pl": vcf_pl,
            "qual": vcf_qual,
        },
        "java_oracle": { "ad": [36, 19], "pl": [542, 0, 1353], "qual": 510.06 },
        "hap_cache_events": hap_cache_events,
        "merged_alleles": merged,
        "numerics": {
            "long_ref": numerics.long_ref,
            "alts": numerics.alts,
            "n_reads": numerics.n_reads,
            "pool_sizes": numerics.pool_sizes,
            "merged_pl": numerics.merged_pl,
            "merged_gls": numerics.merged_gls,
            "merged_ad": numerics.merged_ad,
            "assigned_gt": numerics.assigned_gt,
            "subset_pl": numerics.subset_pl,
            "subset_ad_permuted": numerics.subset_ad_permuted,
            "subset_ad_remarginalized": numerics.subset_ad_remarginalized,
            "n_uninformative_3way": numerics.n_uninformative_3way,
        },
        "qual_counterfactual": {
            "from_rust_pl_298_0_1103": qual_from_rust_pl,
            "from_java_pl_542_0_1353": qual_from_java_pl,
        },
        "site_call": site_call.map(|c| json!({
            "ref": c.event.ref_allele,
            "alt": c.event.alt_allele,
            "extra_alts": c.extra_alt_alleles,
            "post_merge_unused_alt_subset": c.post_merge_unused_alt_subset,
            "gls": c.genotype.genotype_log10_likelihoods,
            "pl": c.genotype.format.pl_as_i32(),
            "ad": c.genotype.format.ad_as_i32(),
            "read_count": c.genotype.aggregation.read_count,
            "n_hap_sums": c.genotype.aggregation.haplotype_log10_sums.len(),
        })),
        "events_at_loc": at_loc.iter().map(|e| format!("{}/{}", e.ref_allele, e.alt_allele)).collect::<Vec<_>>(),
        "n_haps": outcome.assembly.haplotypes.len(),
        "n_likelihood_rows": outcome.read_likelihoods.len(),
    });
    println!("{}", serde_json::to_string_pretty(&doc).unwrap());

    assert_eq!(vcf.reference, "T");
    assert!(vcf.alternate.iter().any(|a| a == "C"));
    assert_eq!(
        vcf_pl,
        vec![542u32, 0, 1353],
        "6R.100 scored-evidence bind matches Java PL"
    );
    assert_eq!(vcf_ad, vec![36u32, 19], "6R.101 remaining remarg FORMAT AD");
    assert_eq!(numerics.subset_pl, vec![542, 0, 1353]);
    assert_eq!(numerics.subset_ad_permuted, vec![34, 17]);
    // Remarginalized AD is diagnostic; equality vs permute is a finding, not a failure.
}
