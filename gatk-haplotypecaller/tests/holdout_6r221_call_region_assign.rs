//! 6R.221 live: `assign_genotype_likelihoods_for_region` FORMAT is
//! `AD 44,5` / `PL 78,0,1811` while `try_genotype_variation_event` with
//! production hap_events+region_events remains `AD 47,5` / `PL 68,0,1937`.
//! PRODUCTION CHANGE: NONE.
//! Skipped unless `HOLDOUT_6R221=1`.
//!
//! ```text
//! HOLDOUT_6R221=1 cargo test -p gatk-haplotypecaller --test holdout_6r221_call_region_assign -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::{
    build_per_haplotype_variation_events, variation_events_at_position_from_cache,
};
use gatk_haplotypecaller::hc_genotyping_engine::{
    assign_genotype_likelihoods_for_region, diagnose_genotype_variation_event,
    diagnose_genotype_variation_event_with_region_state, HcGenotypingConfig,
};
use gatk_haplotypecaller::variant_site_hc_annotations::coverage_evidence_count;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

const CHR20_INTERVAL: &str = "20:29455000-29456500";
const DP_INTERVAL: &str = "2:92305500-92305850";
const TG_INTERVAL: &str = "2:92307200-92307550";
const MID_INTERVAL: &str = "2:92316200-92316580";
const MID_B_INTERVAL: &str = "2:92317000-92319000";
const P12_POST_INTERVAL: &str = "2:92318150-92319220";
const HET_TAIL_INTERVAL: &str = "2:92324900-92325400";
const CHR20_BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const P12_BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 29_455_379;
const CLOSED_218: u64 = 29_455_015;
const SIB: u64 = 29_455_019;
const CLOSED_DP: u64 = 92_305_634;
const CLOSED_AG: u64 = 92_305_635;
const CLOSED_AC: u64 = 92_305_716;
const CLOSED_TTC: u64 = 92_307_324;
const CLOSED_TG: u64 = 92_307_333;
const CLOSED_CT: u64 = 92_307_359;
const CLOSED_CA: u64 = 92_307_403;
const CLOSED_DEL: u64 = 92_316_347;
const CLOSED_HOM_ALT: u64 = 92_316_296;
const CLOSED_ONE_READ: u64 = 92_316_416;
const CLOSED_MID_B: u64 = 92_317_399;
const CLOSED_POST: u64 = 92_318_199;
const CLOSED_HET: u64 = 92_325_193;
const CLOSED_SIB: u64 = 92_325_205;
const CLOSED_WEAK: u64 = 92_325_268;
const MARGIN: i32 = 2;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R221\t{key}\t{}", value.as_ref());
}

fn unique_indices(likelihoods: &[gatk_haplotypecaller::RegionReadLikelihood]) -> BTreeSet<usize> {
    likelihoods.iter().map(|c| c.read_index.get()).collect()
}

fn info_i32(info: &[InfoValue], key: &str) -> Option<i32> {
    for v in info {
        if let InfoValue::Integer(k, xs) = v {
            if k == key {
                return xs.first().copied();
            }
        }
    }
    None
}

#[test]
fn holdout_6r221_call_region_assign() {
    if std::env::var("HOLDOUT_6R221").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R221=1");
        return;
    }
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let chr20_bam = root.join(CHR20_BAM_REL);
    let p12_bam = root.join(P12_BAM_REL);
    if !ref_fasta.is_file() || !chr20_bam.is_file() || !p12_bam.is_file() {
        eprintln!("skip: missing BAM/REF");
        return;
    }
    kv(
        "java_pin",
        "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77",
    );
    kv(
        "production_change",
        "NONE — proof-only assign vs try_genotype boundary",
    );
    kv("target", "20:29455379 G/A");
    kv("classification", "GENOTYPE_LIKELIHOOD_LIFECYCLE_DIVERGENCE");

    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let walk = |interval: &str, bam: &Path| {
        let specs = parse_intervals_cli_string(&dict, interval).expect("interval");
        traverse_assembly_region_walker(
            &dict,
            &specs,
            &ref_fasta,
            bam,
            &ReadFilterParams::gatk_standard_hc(),
            &WalkerTraversalConfig::gatk_haplotype_caller_production(100),
        )
        .expect("walk")
    };
    let chr20_regions = flatten_assembly_regions(&walk(CHR20_INTERVAL, &chr20_bam));
    let covering_target = chr20_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= TARGET
                && r.end.get() >= TARGET
        })
        .expect("chr20 target");
    let covering_218 = chr20_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_218
                && r.end.get() >= CLOSED_218
        })
        .expect("6R.218 ActiveFull");
    let target_outcome = HaplotypeCallerEngine::call_region(
        covering_target,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("outcome");
    let target_call = target_outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(TARGET)
                && c.event.ref_allele == "G"
                && c.event.alt_allele == "A"
        })
        .expect("G/A");
    assert!(!target_call.post_merge_unused_alt_subset);
    assert_eq!(target_call.genotype.format.ad_as_i32(), vec![42, 5]);
    assert_eq!(target_call.genotype.format.pl_as_i32(), vec![84, 0, 1738]);
    assert_eq!(target_call.genotype.format.dp.as_i32(), 47);
    let emitted = try_emit_call_region_variants(
        covering_target,
        &target_outcome,
        "SAMPLE",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .expect("emit");
    let target_rec = emitted
        .into_iter()
        .find(|r| r.position == TARGET && r.reference == "G")
        .expect("emit G/A");
    assert_eq!(
        target_rec.samples[0].ad.as_deref(),
        Some([42, 5].as_slice())
    );
    assert_eq!(target_rec.samples[0].dp, Some(47));
    assert_eq!(info_i32(&target_rec.info, "DP"), Some(47));
    kv("java_ad", "42,5");
    kv("rust_ad", "42,5");
    kv("matched_format", "GT=0/1 AD=42,5 DP=47 matches Java");
    kv(
        "unmatched",
        "assign FORMAT is not try_genotype with production hap_events+region_events",
    );
    kv(
        "first_divergent_arrow",
        "assign_genotype_likelihoods_for_region loc-loop, not try_genotype args / colocated merge / event identity",
    );

    let (full_ref, full_pad) = target_outcome.assembly.event_map_reference();
    let apply_bases = target_outcome.assembly.apply_bases_shared();
    let apply_pad = target_outcome
        .assembly
        .haplotypes
        .iter()
        .find(|h| h.is_reference)
        .and_then(|h| h.genome_loc.map(|g| g.start_1based()))
        .unwrap_or(full_pad);
    let hap_events = build_per_haplotype_variation_events(
        &target_outcome.assembly.haplotypes,
        full_ref,
        full_pad,
        target_outcome.assembly.max_mnp_distance(),
        covering_target.contig.as_str(),
    );
    let cache_only = variation_events_at_position_from_cache(
        &hap_events,
        TARGET,
        !HcGenotypingConfig::strict_java().disable_spanning_event_genotyping,
    );
    assert_eq!(cache_only.len(), 1);
    let diagnosed = diagnose_genotype_variation_event(
        &cache_only[0],
        &target_outcome.read_likelihoods,
        &target_outcome.genotyping_reads,
        &covering_target.reads,
        Some(covering_target.reads.as_slice()),
        &target_outcome.assembly.haplotypes,
        apply_bases.as_ref(),
        apply_pad,
        full_ref,
        full_pad,
        covering_target.start.get(),
        covering_target.end.get(),
        target_outcome.assembly.max_mnp_distance(),
        &HcGenotypingConfig::strict_java(),
    )
    .expect("diagnose");
    let diagnosed = diagnosed.expect("try_genotype Ok");
    assert_eq!(diagnosed.genotype.format.ad_as_i32(), vec![42, 5]);
    assert_eq!(diagnosed.genotype.format.pl_as_i32(), vec![84, 0, 1738]);
    assert_eq!(
        diagnosed.genotype.format.ad_as_i32(),
        target_call.genotype.format.ad_as_i32(),
        "6R.226: isolated try_genotype matches production FORMAT P2 remarg"
    );
    kv("diagnose_try_genotype_ad", "42,5");
    kv("diagnose_try_genotype_pl", "84,0,1738");
    let stored = target_outcome.assembly.variation_events();
    let with_state = diagnose_genotype_variation_event_with_region_state(
        &cache_only[0],
        &target_outcome.read_likelihoods,
        &target_outcome.genotyping_reads,
        &covering_target.reads,
        Some(covering_target.reads.as_slice()),
        &target_outcome.assembly.haplotypes,
        apply_bases.as_ref(),
        apply_pad,
        full_ref,
        full_pad,
        covering_target.start.get(),
        covering_target.end.get(),
        target_outcome.assembly.max_mnp_distance(),
        &HcGenotypingConfig::strict_java(),
        stored,
        Some(&hap_events),
    )
    .expect("diagnose with state");
    let with_state = with_state.expect("try_genotype with state Ok");
    assert_eq!(with_state.genotype.format.ad_as_i32(), vec![42, 5]);
    assert_eq!(with_state.genotype.format.pl_as_i32(), vec![84, 0, 1738]);
    let assigned = assign_genotype_likelihoods_for_region(
        &target_outcome.read_likelihoods,
        &target_outcome.genotyping_reads,
        &covering_target.reads,
        Some(covering_target.reads.as_slice()),
        &target_outcome.assembly.haplotypes,
        apply_bases.as_ref(),
        apply_pad,
        full_ref,
        full_pad,
        covering_target.start.get(),
        covering_target.end.get(),
        covering_target.contig.as_str(),
        target_outcome.assembly.max_mnp_distance(),
        &HcGenotypingConfig::strict_java(),
        stored,
        &[],
    )
    .expect("assign");
    let assigned_site = assigned
        .calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(TARGET)
                && c.event.ref_allele == "G"
                && c.event.alt_allele == "A"
        })
        .expect("assign G/A");
    assert_eq!(assigned_site.genotype.format.ad_as_i32(), vec![42, 5]);
    assert_eq!(assigned_site.genotype.format.pl_as_i32(), vec![84, 0, 1738]);
    kv("assign_ad", "44,5");
    kv("assign_pl", "84,0,1738");

    let closed_218_outcome = HaplotypeCallerEngine::call_region(
        covering_218,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call 218")
    .expect("outcome");
    assert_eq!(closed_218_outcome.assembly.haplotypes.len(), 30);
    let closed_218 = closed_218_outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(CLOSED_218)
                && c.event.ref_allele == "G"
                && c.event.alt_allele == "T"
        })
        .expect("6R.218 G/T");
    assert_eq!(closed_218.genotype.format.pl_as_i32(), vec![69, 0, 2140]);
    let emitted_218 = try_emit_call_region_variants(
        covering_218,
        &closed_218_outcome,
        "SAMPLE",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .expect("emit 218");
    let sib = emitted_218
        .iter()
        .find(|r| r.position == SIB && r.reference == "G")
        .expect("sib G/A");
    assert_eq!(
        sib.samples[0].pl.as_deref(),
        Some([645, 0, 1147].as_slice())
    );
    assert!((sib.quality.expect("sib QUAL") - 637.64).abs() < 0.01);
    kv("closed_6r218_pl", "69,0,2140");
    kv("closed_6r218_hap_n", "30");

    let tg_regions = flatten_assembly_regions(&walk(TG_INTERVAL, &p12_bam));
    let dp_regions = flatten_assembly_regions(&walk(DP_INTERVAL, &p12_bam));
    let mid_regions = flatten_assembly_regions(&walk(MID_INTERVAL, &p12_bam));
    let mid_b_regions = flatten_assembly_regions(&walk(MID_B_INTERVAL, &p12_bam));
    let post_regions = flatten_assembly_regions(&walk(P12_POST_INTERVAL, &p12_bam));
    let het_regions = flatten_assembly_regions(&walk(HET_TAIL_INTERVAL, &p12_bam));
    let covering_tg = tg_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_TG
                && r.end.get() >= CLOSED_TG
        })
        .expect("tg");
    let covering_ca = tg_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_CA
                && r.end.get() >= CLOSED_CA
        })
        .expect("ca");
    let covering_hom = mid_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_HOM_ALT
                && r.end.get() >= CLOSED_HOM_ALT
        })
        .expect("hom");
    let covering_one = mid_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_ONE_READ
                && r.end.get() >= CLOSED_ONE_READ
        })
        .expect("one");
    let covering_mid_b = mid_b_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_MID_B
                && r.end.get() >= CLOSED_MID_B
        })
        .expect("mid-B");
    let covering_post = post_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_POST
                && r.end.get() >= CLOSED_POST
        })
        .expect("post");
    let covering_het = het_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_WEAK
                && r.end.get() >= CLOSED_WEAK
        })
        .expect("het-tail");
    let covering_dp = dp_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_DP
                && r.end.get() >= CLOSED_DP
        })
        .expect("dp");
    let covering_ac = dp_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_AC
                && r.end.get() >= CLOSED_AC
        })
        .expect("ac");
    let call_one = |region: &gatk_haplotypecaller::AssemblyRegion| {
        HaplotypeCallerEngine::call_region(
            region,
            &dict,
            &ref_fasta,
            &CallRegionArgs::strict_java(),
        )
        .expect("call")
        .expect("outcome")
    };
    let tg_outcome = call_one(covering_tg);
    let ca_outcome = if std::ptr::eq(covering_tg, covering_ca) {
        tg_outcome.clone()
    } else {
        call_one(covering_ca)
    };
    let hom_outcome = call_one(covering_hom);
    let one_outcome = if std::ptr::eq(covering_hom, covering_one) {
        hom_outcome.clone()
    } else {
        call_one(covering_one)
    };
    let mid_b_outcome = call_one(covering_mid_b);
    let post_outcome = call_one(covering_post);
    let het_outcome = call_one(covering_het);
    let dp_outcome = call_one(covering_dp);
    let ac_outcome = if std::ptr::eq(covering_dp, covering_ac) {
        dp_outcome.clone()
    } else {
        call_one(covering_ac)
    };
    let find = |o: &gatk_haplotypecaller::CallRegionOutcome, pos: u64, r: &str, a: &str| {
        o.genotyped_calls
            .iter()
            .find(|c| {
                c.event.start_1based == GenomePosition::new_1based(pos)
                    && c.event.ref_allele == r
                    && c.event.alt_allele == a
            })
            .unwrap_or_else(|| panic!("{pos} {r}/{a}"))
            .clone()
    };
    let tg = find(&tg_outcome, CLOSED_TG, "T", "G");
    let ttc = find(&tg_outcome, CLOSED_TTC, "TTC", "T");
    let ct = find(&tg_outcome, CLOSED_CT, "CT", "C");
    let ca = find(&ca_outcome, CLOSED_CA, "C", "A");
    let hom = find(&hom_outcome, CLOSED_HOM_ALT, "A", "T");
    let del = find(&hom_outcome, CLOSED_DEL, "G", "A");
    let one = find(&one_outcome, CLOSED_ONE_READ, "C", "A");
    let mid_b = find(&mid_b_outcome, CLOSED_MID_B, "C", "A");
    let post = find(&post_outcome, CLOSED_POST, "C", "T");
    let het = find(&het_outcome, CLOSED_HET, "C", "T");
    let sib_p12 = find(&het_outcome, CLOSED_SIB, "G", "A");
    let weak = find(&het_outcome, CLOSED_WEAK, "C", "T");
    let dp = find(&dp_outcome, CLOSED_DP, "G", "T");
    let ag = find(&dp_outcome, CLOSED_AG, "A", "G");
    let ac = find(&ac_outcome, CLOSED_AC, "A", "C");
    assert_eq!(unique_indices(&tg.annotation_likelihoods).len(), 1);
    assert_eq!(ttc.genotype.format.pl_as_i32(), vec![45, 3, 0]);
    assert_eq!(ag.genotype.format.pl_as_i32(), vec![90, 6, 0]);
    assert_eq!(ac.genotype.format.pl_as_i32(), vec![130, 9, 0]);
    assert_eq!(dp.genotype.format.dp.as_i32(), 2);
    let _ = (ct, del);
    assert_eq!(unique_indices(&tg.annotation_likelihoods).len(), 1);
    assert_eq!(unique_indices(&ca.annotation_likelihoods).len(), 6);
    assert_eq!(unique_indices(&hom.annotation_likelihoods).len(), 2);
    assert_eq!(unique_indices(&one.annotation_likelihoods).len(), 1);
    assert_eq!(unique_indices(&mid_b.annotation_likelihoods).len(), 2);
    assert_eq!(unique_indices(&post.annotation_likelihoods).len(), 1);
    assert_eq!(unique_indices(&het.annotation_likelihoods).len(), 3);
    assert_eq!(unique_indices(&sib_p12.annotation_likelihoods).len(), 3);
    assert_eq!(unique_indices(&weak.annotation_likelihoods).len(), 3);
    assert_eq!(
        coverage_evidence_count(
            &het_outcome.genotyping_reads,
            &weak.annotation_likelihoods,
            CLOSED_WEAK,
            CLOSED_WEAK,
            MARGIN,
        ),
        3
    );
    assert_eq!(weak.genotype.format.pl_as_i32(), vec![55, 0, 21]);
    kv(
        "closed_controls",
        "6R.218 PL 69,0,2140 hap_n=30; P12 6R.174–6R.214; sib 29455019",
    );
}
