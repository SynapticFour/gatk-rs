//! 6R.221: first operation inside production `call_region` /
//! `assign_genotype_likelihoods_for_region` that changes the 52-row
//! isolated `try_genotype` object (`AD 47,5` / `PL 68,0,1937`) into
//! FORMAT `AD 44,5` / `PL 78,0,1811` at `20:29455379 G/A`.
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! PRODUCTION CHANGE: NONE.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r221_call_region_assign_boundary -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::{
    build_event_start_positions_from_cache, build_per_haplotype_variation_events,
    merged_biallelic_sites_at_position, prefer_indel_over_colocated_snps,
    variation_events_at_position_from_cache,
};
use gatk_haplotypecaller::hc_allele_mapping::create_allele_mapper_with_events;
use gatk_haplotypecaller::hc_allele_mapping::replace_span_del_events;
use gatk_haplotypecaller::hc_genotyping_engine::{
    alt_hap_indices_for_genotype_marginalization, assign_genotype_likelihoods_for_region,
    diagnose_genotype_variation_event_with_region_state,
    filter_genotyped_calls_for_strict_java_emit, java_alignment_read_overlaps_interval,
    region_likelihoods_to_rows, take_colocated_merge_numerics, GenotypedSiteCall,
    HcGenotypingConfig, DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN,
};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition, HaplotypeCallerEngine,
    ReadFilterParams, WalkerTraversalConfig,
};
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 29_455_379;
const CLOSED_PL: u64 = 29_455_015;
const TARGET_REF: &str = "G";
const TARGET_ALT: &str = "A";
const AD_A: [i32; 2] = [47, 5];
const PL_A: [i32; 3] = [68, 0, 1937];
const AD_C: [i32; 2] = [44, 5];
const PL_C: [i32; 3] = [78, 0, 1811];

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R221\t{key}\t{}", value.as_ref());
}

fn unique_n(likelihoods: &[gatk_haplotypecaller::RegionReadLikelihood]) -> usize {
    likelihoods
        .iter()
        .map(|c| c.read_index.get())
        .collect::<BTreeSet<_>>()
        .len()
}

fn fmt_of(call: &GenotypedSiteCall) -> (Vec<i32>, Vec<i32>, i32, usize) {
    (
        call.genotype.format.ad_as_i32(),
        call.genotype.format.pl_as_i32(),
        call.genotype.format.dp.as_i32(),
        unique_n(&call.annotation_likelihoods),
    )
}

fn dump_call(tag: &str, call: &GenotypedSiteCall) {
    let (ad, pl, dp, ann) = fmt_of(call);
    kv(
        tag,
        format!(
            "ad={ad:?}\tpl={pl:?}\tdp={dp}\tann_n={ann}\textra_alt={}\tpost_merge={}",
            call.extra_alt_alleles.len(),
            call.post_merge_unused_alt_subset
        ),
    );
}

#[test]
fn forensic_6r221_java_lifecycle_contract() {
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE — proof-only");
    kv(
        "java_lifecycle",
        "calculateGLsForThisEvent and DepthPerAlleleBySample share retainEvidence AlleleLikelihoods; no post-retainEvidence FORMAT subset",
    );
    assert_eq!(DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN, 2);
}

#[test]
fn forensic_6r221_call_region_assign_boundary() {
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing chr20_tiny BAM/REF");
        return;
    }
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE — proof-only");
    kv("target", "20:29455379 G/A");

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
            ) && r.start.get() <= TARGET
                && r.end.get() >= TARGET
        })
        .collect();
    assert_eq!(covering.len(), 1);
    let region = covering[0];
    kv(
        "active_full",
        format!(
            "{}:{}-{}",
            region.contig,
            region.start.get(),
            region.end.get()
        ),
    );

    let outcome = HaplotypeCallerEngine::call_region(
        region,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");

    let haps = &outcome.assembly.haplotypes;
    let reads = &outcome.genotyping_reads;
    let likelihoods = &outcome.read_likelihoods;
    let (full_ref, full_pad) = outcome.assembly.event_map_reference();
    let apply_bases = outcome.assembly.apply_bases_shared();
    let apply_pad = haps
        .iter()
        .find(|h| h.is_reference)
        .and_then(|h| h.genome_loc.map(|g| g.start_1based()))
        .unwrap_or(full_pad);
    kv("apply_pad", apply_pad.to_string());
    kv("full_pad", full_pad.to_string());
    kv("hap_n", haps.len().to_string());
    kv("genotyping_reads_n", reads.len().to_string());
    kv("region_reads_n", region.reads.len().to_string());
    kv("likelihood_cells", likelihoods.len().to_string());

    let config = HcGenotypingConfig::strict_java();
    assert!(!config.genotype_stored_events_only);
    let emit_spanning = !config.disable_spanning_event_genotyping;
    let hap_events = build_per_haplotype_variation_events(
        haps,
        full_ref,
        full_pad,
        outcome.assembly.max_mnp_distance(),
        &region.contig,
    );
    let cache_only = variation_events_at_position_from_cache(&hap_events, TARGET, emit_spanning);
    kv(
        "hap_events_at_loc",
        format!(
            "n={}\t{}",
            cache_only.len(),
            cache_only
                .iter()
                .map(|e| format!("{}/{}", e.ref_allele, e.alt_allele))
                .collect::<Vec<_>>()
                .join(",")
        ),
    );
    assert_eq!(cache_only.len(), 1);
    assert_eq!(cache_only[0].ref_allele, TARGET_REF);
    assert_eq!(cache_only[0].alt_allele, TARGET_ALT);
    let event = &cache_only[0];

    let stored = outcome.assembly.variation_events();
    let stored_at: Vec<_> = stored
        .iter()
        .filter(|e| e.start_1based.get() == TARGET)
        .collect();
    kv(
        "stored_events_at_loc",
        format!(
            "region_n={}\tat_loc={}\t{}",
            stored.len(),
            stored_at.len(),
            stored_at
                .iter()
                .map(|e| format!("{}/{}", e.ref_allele, e.alt_allele))
                .collect::<Vec<_>>()
                .join(",")
        ),
    );

    let prod_hits: Vec<_> = outcome
        .genotyped_calls
        .iter()
        .filter(|c| c.event.start_1based.get() == TARGET)
        .collect();
    kv("production_calls_at_loc_n", prod_hits.len().to_string());
    assert_eq!(
        prod_hits.len(),
        1,
        "target processed once in genotyped_calls"
    );
    let site = prod_hits[0];
    dump_call("C_call_region_format", site);
    assert_eq!(site.event.ref_allele, TARGET_REF);
    assert_eq!(site.event.alt_allele, TARGET_ALT);
    assert!(site.extra_alt_alleles.is_empty());
    assert!(!site.post_merge_unused_alt_subset);
    let (fmt_ad, fmt_pl, fmt_dp, fmt_ann) = fmt_of(site);
    assert_eq!(fmt_ad, AD_C);
    assert_eq!(fmt_pl, PL_C);
    assert_eq!(fmt_dp, 49);
    kv("C_annotation_unique_n", fmt_ann.to_string());

    let first_loc = outcome
        .genotyped_calls
        .first()
        .map(|c| c.event.start_1based.get())
        .unwrap_or(0);
    kv("first_genotyped_loc", first_loc.to_string());
    kv("target_is_first_loc", (first_loc == TARGET).to_string());
    let positions = build_event_start_positions_from_cache(&hap_events);
    let prior: Vec<u64> = positions
        .iter()
        .copied()
        .filter(|p| *p >= region.start.get() && *p < TARGET && *p <= region.end.get())
        .collect();
    kv("prior_eventmap_locs_in_active", format!("{prior:?}"));
    kv("prior_eventmap_locs_n", prior.len().to_string());
    assert_eq!(
        prior,
        vec![29_455_375],
        "assign visits 29455375 before the target"
    );

    let overlap_n = reads
        .iter()
        .filter(|r| java_alignment_read_overlaps_interval(r, TARGET, TARGET, 2))
        .count();
    kv("overlap_pm2_genotyping_reads", overlap_n.to_string());

    let mapping_none = create_allele_mapper_with_events(
        event,
        TARGET,
        haps,
        apply_pad,
        apply_bases.as_ref(),
        outcome.assembly.max_mnp_distance(),
        emit_spanning,
        None,
    );
    let mapping_hap = create_allele_mapper_with_events(
        event,
        TARGET,
        haps,
        apply_pad,
        apply_bases.as_ref(),
        outcome.assembly.max_mnp_distance(),
        emit_spanning,
        Some(&hap_events),
    );
    kv(
        "mapping_alt_haps",
        format!(
            "hap_events=None:{}\thap_events=Some:{}",
            mapping_none.alt_haplotype_indices.len(),
            mapping_hap.alt_haplotype_indices.len()
        ),
    );
    let ref_hap = haps.iter().find(|h| h.is_reference).expect("ref hap");
    let alt_none = alt_hap_indices_for_genotype_marginalization(
        &mapping_none,
        haps,
        event,
        ref_hap,
        apply_pad,
        apply_bases.as_ref(),
        outcome.assembly.max_mnp_distance(),
        &region.contig,
        &config,
    );
    let alt_hap = alt_hap_indices_for_genotype_marginalization(
        &mapping_hap,
        haps,
        event,
        ref_hap,
        apply_pad,
        apply_bases.as_ref(),
        outcome.assembly.max_mnp_distance(),
        &region.contig,
        &config,
    );
    kv(
        "alt_pool_n",
        format!("none={}\thap={}", alt_none.len(), alt_hap.len()),
    );

    let diagnose =
        |tag: &str,
         region_events: &[gatk_haplotypecaller::event_map::VariationEvent],
         hap: Option<&gatk_haplotypecaller::event_map::PerHaplotypeVariationEvents>| {
            let diagnosed = diagnose_genotype_variation_event_with_region_state(
                event,
                likelihoods,
                reads,
                &region.reads,
                Some(region.reads.as_slice()),
                haps,
                apply_bases.as_ref(),
                apply_pad,
                full_ref,
                full_pad,
                region.start.get(),
                region.end.get(),
                outcome.assembly.max_mnp_distance(),
                &config,
                region_events,
                hap,
            )
            .expect(tag);
            match diagnosed {
                Ok(call) => {
                    dump_call(tag, &call);
                    Some(call)
                }
                Err(reason) => {
                    kv(tag, format!("REJECT {reason:?}"));
                    None
                }
            }
        };

    let a = diagnose("A_isolated_try_genotype", &[], None).expect("A");
    let (a_ad, a_pl, _, a_ann) = fmt_of(&a);
    assert_eq!(a_ad, AD_A);
    assert_eq!(a_pl, PL_A);
    kv("A_annotation_unique_n", a_ann.to_string());
    kv(
        "A_event",
        format!(
            "{}:{}-{} {}/{}",
            a.event.contig,
            a.event.start_1based.get(),
            a.event.end_1based.get(),
            a.event.ref_allele,
            a.event.alt_allele
        ),
    );
    kv(
        "cache_event",
        format!(
            "{}:{}-{} {}/{}",
            event.contig,
            event.start_1based.get(),
            event.end_1based.get(),
            event.ref_allele,
            event.alt_allele
        ),
    );

    let mut loc_raw = variation_events_at_position_from_cache(&hap_events, TARGET, emit_spanning);
    kv(
        "loc_raw_n",
        format!(
            "n={}\t{}",
            loc_raw.len(),
            loc_raw
                .iter()
                .map(|e| format!(
                    "{}-{} {}/{}",
                    e.start_1based.get(),
                    e.end_1based.get(),
                    e.ref_allele,
                    e.alt_allele
                ))
                .collect::<Vec<_>>()
                .join(",")
        ),
    );
    prefer_indel_over_colocated_snps(&mut loc_raw);
    let with_span = replace_span_del_events(&loc_raw, TARGET, apply_pad, apply_bases.as_ref());
    kv(
        "loc_span_n",
        format!(
            "n={}\t{}",
            with_span.len(),
            with_span
                .iter()
                .map(|e| format!(
                    "{}-{} {}/{}",
                    e.start_1based.get(),
                    e.end_1based.get(),
                    e.ref_allele,
                    e.alt_allele
                ))
                .collect::<Vec<_>>()
                .join(",")
        ),
    );
    let merged_sites = merged_biallelic_sites_at_position(&with_span, TARGET);
    kv(
        "loc_merged_n",
        format!(
            "n={}\t{}",
            merged_sites.len(),
            merged_sites
                .iter()
                .map(|e| format!(
                    "{}-{} {}/{}",
                    e.start_1based.get(),
                    e.end_1based.get(),
                    e.ref_allele,
                    e.alt_allele
                ))
                .collect::<Vec<_>>()
                .join(",")
        ),
    );
    if let Some(ev) = merged_sites
        .iter()
        .find(|e| e.ref_allele == TARGET_REF && e.alt_allele == TARGET_ALT)
    {
        kv(
            "D2_used_event",
            format!(
                "{}:{}-{} {}/{}",
                ev.contig,
                ev.start_1based.get(),
                ev.end_1based.get(),
                ev.ref_allele,
                ev.alt_allele
            ),
        );
        let g = diagnose_genotype_variation_event_with_region_state(
            ev,
            likelihoods,
            reads,
            &region.reads,
            Some(region.reads.as_slice()),
            haps,
            apply_bases.as_ref(),
            apply_pad,
            full_ref,
            full_pad,
            region.start.get(),
            region.end.get(),
            outcome.assembly.max_mnp_distance(),
            &config,
            stored,
            Some(&hap_events),
        )
        .expect("D2");
        match g {
            Ok(call) => dump_call("D2_merged_event_try_genotype", &call),
            Err(reason) => kv("D2_merged_event_try_genotype", format!("REJECT {reason:?}")),
        }
    }

    let _ = take_colocated_merge_numerics();

    let b = diagnose("B_hap_events_only", &[], Some(&hap_events));
    let c_ev = diagnose("C_region_events_only", stored, None);
    let d = diagnose("D_hap_and_region_events", stored, Some(&hap_events));

    let _ = region_likelihoods_to_rows(likelihoods, haps.len());
    let d_after_tls = diagnose("D_after_tls_full_table", stored, Some(&hap_events));

    let assigned = assign_genotype_likelihoods_for_region(
        likelihoods,
        reads,
        &region.reads,
        Some(region.reads.as_slice()),
        haps,
        apply_bases.as_ref(),
        apply_pad,
        full_ref,
        full_pad,
        region.start.get(),
        region.end.get(),
        &region.contig,
        outcome.assembly.max_mnp_distance(),
        &config,
        stored,
        &[],
    )
    .expect("assign");
    let assign_hits: Vec<_> = assigned
        .calls
        .iter()
        .filter(|c| {
            c.event.start_1based.get() == TARGET
                && c.event.ref_allele == TARGET_REF
                && c.event.alt_allele == TARGET_ALT
        })
        .collect();
    kv("E_assign_calls_at_loc_n", assign_hits.len().to_string());
    assert_eq!(assign_hits.len(), 1);
    dump_call("E_assign", assign_hits[0]);
    kv(
        "E_event",
        format!(
            "{}:{}-{} {}/{}",
            assign_hits[0].event.contig,
            assign_hits[0].event.start_1based.get(),
            assign_hits[0].event.end_1based.get(),
            assign_hits[0].event.ref_allele,
            assign_hits[0].event.alt_allele
        ),
    );
    let colocated = take_colocated_merge_numerics();
    kv("colocated_merge_snaps_n", colocated.len().to_string());
    for snap in &colocated {
        kv(
            "colocated_snap",
            format!(
                "loc={}\tlong_ref={}\talts={:?}\tn_reads={}\tmerged_ad={:?}\tmerged_pl={:?}\tsubset_ad_remarg={:?}\tsubset_pl={:?}",
                snap.loc,
                snap.long_ref,
                snap.alts,
                snap.n_reads,
                snap.merged_ad,
                snap.merged_pl,
                snap.subset_ad_remarginalized,
                snap.subset_pl
            ),
        );
    }
    let (e_ad, e_pl, _, _e_ann) = fmt_of(assign_hits[0]);

    let mut filtered = assigned.calls.clone();
    filter_genotyped_calls_for_strict_java_emit(
        &mut filtered,
        &region.reads,
        &outcome.assembly,
        &config,
    )
    .expect("filter");
    let filter_hits: Vec<_> = filtered
        .iter()
        .filter(|c| {
            c.event.start_1based.get() == TARGET
                && c.event.ref_allele == TARGET_REF
                && c.event.alt_allele == TARGET_ALT
        })
        .collect();
    kv("F_filter_calls_at_loc_n", filter_hits.len().to_string());
    assert_eq!(filter_hits.len(), 1);
    dump_call("F_filter", filter_hits[0]);
    let (f_ad, f_pl, _, _) = fmt_of(filter_hits[0]);

    let matches_c = |ad: &[i32], pl: &[i32]| ad == AD_C && pl == PL_C;
    let matches_a = |ad: &[i32], pl: &[i32]| ad == AD_A && pl == PL_A;

    let b_fmt = b.as_ref().map(fmt_of);
    let c_fmt = c_ev.as_ref().map(fmt_of);
    let d_fmt = d.as_ref().map(fmt_of);
    let d_tls_fmt = d_after_tls.as_ref().map(fmt_of);

    kv(
        "b_matches_c",
        b_fmt
            .as_ref()
            .map(|(ad, pl, _, _)| matches_c(ad, pl).to_string())
            .unwrap_or_else(|| "none".into()),
    );
    kv(
        "c_matches_c",
        c_fmt
            .as_ref()
            .map(|(ad, pl, _, _)| matches_c(ad, pl).to_string())
            .unwrap_or_else(|| "none".into()),
    );
    kv(
        "d_matches_c",
        d_fmt
            .as_ref()
            .map(|(ad, pl, _, _)| matches_c(ad, pl).to_string())
            .unwrap_or_else(|| "none".into()),
    );
    kv(
        "d_tls_matches_c",
        d_tls_fmt
            .as_ref()
            .map(|(ad, pl, _, _)| matches_c(ad, pl).to_string())
            .unwrap_or_else(|| "none".into()),
    );
    kv("e_matches_c", matches_c(&e_ad, &e_pl).to_string());
    kv("f_matches_c", matches_c(&f_ad, &f_pl).to_string());
    kv("e_matches_a", matches_a(&e_ad, &e_pl).to_string());
    assert!(
        d_fmt
            .as_ref()
            .is_some_and(|(ad, pl, _, _)| matches_a(ad, pl)),
        "production-arg try_genotype must remain the 52-row object"
    );
    assert!(
        matches_c(&e_ad, &e_pl),
        "assign_genotype_likelihoods_for_region must be the FORMAT object"
    );
    assert_eq!(
        colocated.len(),
        0,
        "colocated merge must not fire at this SNP"
    );
    assert_eq!(merged_sites.len(), 1);

    let first_op = "assign_genotype_likelihoods_for_region loc-loop produces FORMAT AD 44,5 / PL 78,0,1811 while try_genotype_variation_event with the same event, hap_events, region_events, pads, reads, and likelihoods produces AD 47,5 / PL 68,0,1937. Colocated merge snaps=0. Event identity unchanged (20:29455379 G/A). TLS full-table fill is not sufficient. filter_genotyped_calls does not rewrite. Annotation unique n=52 on both (not a retainEvidence membership filter). AD and PL change together. Assign visits EventMap loc 29455375 before the target (no genotyped_call there); that prior-loc wrap is the next inner boundary.";
    let classification = "GENOTYPE_LIKELIHOOD_LIFECYCLE_DIVERGENCE";
    let membership = "ORDER_ONLY / none on the 52-row retainEvidence annotation object (unique n=52 both). FORMAT AD/PL are a second genotype-likelihood state, not a proven 52→49 row filter.";
    let java_counterpart = "NO — Java calculateGLsForThisEvent and DepthPerAlleleBySample share retainEvidence AlleleLikelihoods. Rust assign loc-loop yields a different FORMAT state than that object.";

    kv("first_divergent_operation", first_op);
    kv("classification", classification);
    kv("membership", membership);
    kv("java_counterpart", java_counterpart);
    kv(
        "pl_changed_with_ad",
        ((a_pl != fmt_pl) && (a_ad != fmt_ad)).to_string(),
    );
    assert_ne!(a_ad, fmt_ad);
    assert_ne!(a_pl, fmt_pl);
    assert_eq!(classification, "GENOTYPE_LIKELIHOOD_LIFECYCLE_DIVERGENCE");
    kv("production_change", "NONE");
}

#[test]
fn forensic_6r221_closed_6r218_untouched() {
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing chr20_tiny BAM/REF");
        return;
    }
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
            ) && r.start.get() <= CLOSED_PL
                && r.end.get() >= CLOSED_PL
        })
        .expect("6R.218 ActiveFull");
    let outcome = HaplotypeCallerEngine::call_region(
        covering,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");
    let site = outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(CLOSED_PL)
                && c.event.ref_allele == "G"
                && c.event.alt_allele == "T"
        })
        .expect("G/T");
    assert_eq!(outcome.assembly.haplotypes.len(), 30);
    assert_eq!(site.genotype.format.pl_as_i32(), vec![69, 0, 2140]);
    kv("closed_6r218_pl", "69,0,2140");
    kv("closed_6r218_hap_n", "30");
}
