//! 6R.193: emit BaseQRankSum from the same fillQualsFromLikelihood object as
//! ReadPosRankSum. Finite zero is `Some(0.0)` and must be inserted.
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! Target `2:92307359 CT/C`. MQRankSum closed by 6R.194.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r193_baseq_rank_sum_emission -- --nocapture --test-threads=1
//! HOLDOUT_6R193=1 cargo test -p gatk-haplotypecaller --test holdout_6r193_baseq_rank_sum -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::variant_site_hc_annotations::{
    baseq_rank_sum_quals_from_likelihoods, read_pos_rank_sum_quals_from_likelihoods,
    HcStrandBiasLikelihoods,
};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
const GAP_INTERVAL: &str = "2:92305500-92305850";
const TG_INTERVAL: &str = "2:92307200-92307550";
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 92_307_359;
const CLOSED_GT: u64 = 92_305_634;
const CLOSED_AG: u64 = 92_305_635;
const CLOSED_AC: u64 = 92_305_716;
const CLOSED_INDEL: u64 = 92_307_324;
const CLOSED_TG: u64 = 92_307_333;
const CLOSED_MID: u64 = 92_316_347;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R193\t{key}\t{}", value.as_ref());
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

fn info_f64(info: &[InfoValue], key: &str) -> Option<f64> {
    for v in info {
        if let InfoValue::Float(k, xs) = v {
            if k == key {
                return xs.first().copied();
            }
        }
    }
    None
}

fn info_has(info: &[InfoValue], key: &str) -> bool {
    info.iter().any(|v| match v {
        InfoValue::Flag(k)
        | InfoValue::Integer(k, _)
        | InfoValue::Float(k, _)
        | InfoValue::String(k, _)
        | InfoValue::Character(k, _) => k == key,
    })
}

#[test]
fn forensic_6r193_source_contracts() {
    kv("java_pin", JAVA_PIN);
    kv(
        "production_change",
        "BaseQRankSum from fillQualsFromLikelihood on the same annotation likelihoods as ReadPos",
    );
    kv("classification", "F — WRONG EMISSION PREDICATE (closed)");
    kv("target", "2:92307359 CT/C");

    let ann = include_str!("../src/variant_site_hc_annotations.rs");
    let plugin = include_str!("../src/annotator/plugins/rank_sum_baseq.rs");
    let rp = include_str!("../src/annotator/plugins/read_pos_rank_sum.rs");
    let emit = include_str!("../src/region_vcf_emit.rs");
    let pipe = include_str!("../src/hc_genotyping_engine/genotype_site_pipeline.rs");
    let early = include_str!("../src/hc_genotyping_engine/genotype_site_early_template.rs");
    let production_plugin = plugin.split("#[cfg(test)]").next().expect("plugin");

    let fill = ann
        .split("fn fill_ranksum_quals_from_likelihood")
        .nth(1)
        .expect("shared fillQuals")
        .split("pub fn rms_mapping_quality_sample_reads")
        .next()
        .expect("fill body");
    assert!(
        fill.contains("evidence.likelihoods")
            && fill.contains("marginalize_rows_to_biallelic_alleles")
            && fill.contains("LOG_10_INFORMATIVE_THRESHOLD")
            && fill.contains("java_read_pos_rank_element")
            && fill.contains("java_baseq_rank_element"),
        "ReadPos and BaseQ share fillQualsFromLikelihood membership"
    );
    assert!(
        fill.contains("read_base_quality_at_ref_coord_1based")
            || fill.contains("java_baseq_rank_element"),
        "BaseQ element is getReadBaseQualityAtReferenceCoordinate, not ReadPos"
    );
    assert!(
        !fill.contains("region.reads") && !fill.contains("read_offset_evidence_at_site"),
        "must not reconstruct pileup"
    );
    assert!(
        !fill.contains("92307359") && !fill.contains("92305759"),
        "no locus special case"
    );

    let annotate = ann
        .split("pub fn annotate_hc_variant_site")
        .nth(1)
        .expect("annotate");
    let annotate_body = annotate
        .split("fn read_offset_evidence_at_site")
        .next()
        .expect("body");
    assert!(
        annotate_body.contains("baseq_rank_sum_quals_from_likelihoods")
            && annotate_body.contains("read_pos_rank_sum_quals_from_likelihoods")
            && annotate_body.contains("rank_sum_baseq::base_quality_rank_sum"),
        "annotate binds BaseQ lists from the same strand_likelihoods as ReadPos"
    );
    assert!(
        annotate_body.contains("mapping_quality_rank_sum"),
        "6R.194 wires MQRankSum from the same fillQuals membership"
    );

    assert!(
        production_plugin.contains(
            "pub fn base_quality_rank_sum(ref_quals: &[f64], alt_quals: &[f64]) -> Option<f64>"
        ) && production_plugin.contains("if alt.is_empty() || reference.is_empty()")
            && production_plugin.contains("return None")
            && production_plugin.contains("result.z.is_nan()"),
        "6R.172 Option semantics on BaseQ"
    );
    assert!(
        rp.contains(
            "pub fn read_pos_rank_sum(ref_positions: &[f64], alt_positions: &[f64]) -> Option<f64>"
        ),
        "ReadPosRankSum formula/Option stays"
    );

    let info_fn = emit
        .split("fn hc_info_values")
        .nth(1)
        .expect("hc_info_values");
    let info_body = info_fn
        .split("fn try_emit_call_region_variants")
        .next()
        .expect("body");
    assert!(
        info_body.contains("if let Some(z) = ann.base_q_rank_sum")
            && info_body.contains("InfoValue::Float(\"BaseQRankSum\""),
        "Some(z) including 0.0 inserts BaseQRankSum"
    );
    assert!(
        !info_body.contains("if z != 0") && !info_body.contains("if z == 0"),
        "must not omit legitimate zero"
    );
    assert!(
        info_body.contains("InfoValue::Float(\"MQRankSum\""),
        "6R.194 emits MQRankSum from the same Some(z) contract"
    );
    assert!(
        pipe.contains("let annotation = annotation_likelihoods_from_stored_haplotypes")
            && early.contains("if is_cluster_tg_snp(&event)"),
        "6R.189 SiteScore attach and 6R.186 cluster-TG stay"
    );
    assert!(
        !ann.contains("92307359") && !emit.contains("92307359") && !pipe.contains("92307359"),
        "no locus-specific 6R.193 production patch"
    );
}

#[test]
fn forensic_6r193_baseq_rank_sum_emission() {
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing P12 ref/BAM");
        return;
    }
    kv("java_pin", JAVA_PIN);
    kv(
        "production_change",
        "BaseQRankSum from fillQualsFromLikelihood on the same annotation likelihoods as ReadPos",
    );
    kv("target", "2:92307359 CT/C");
    kv(
        "java_path",
        "StandardAnnotation → BaseQualityRankSumTest → RankSumTest.fillQualsFromLikelihood → getReadBaseQualityAtReferenceCoordinate → finite z=0.0 → INFO insert",
    );
    kv(
        "rust_path",
        "annotate_hc_variant_site → baseq_rank_sum_quals_from_likelihoods(strand_likelihoods) → Option Some(0.0) → hc_info_values insert",
    );

    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let gap_specs = parse_intervals_cli_string(&dict, GAP_INTERVAL).expect("gap");
    let tg_specs = parse_intervals_cli_string(&dict, TG_INTERVAL).expect("tg");
    let gap_walk = traverse_assembly_region_walker(
        &dict,
        &gap_specs,
        &ref_fasta,
        &bam,
        &ReadFilterParams::gatk_standard_hc(),
        &WalkerTraversalConfig::gatk_haplotype_caller_production(100),
    )
    .expect("gap walk");
    let tg_walk = traverse_assembly_region_walker(
        &dict,
        &tg_specs,
        &ref_fasta,
        &bam,
        &ReadFilterParams::gatk_standard_hc(),
        &WalkerTraversalConfig::gatk_haplotype_caller_production(100),
    )
    .expect("tg walk");
    let gap_regions = flatten_assembly_regions(&gap_walk);
    let tg_regions = flatten_assembly_regions(&tg_walk);
    let covering_ag = gap_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_AG
                && r.end.get() >= CLOSED_AG
        })
        .expect("covering A/G");
    let covering_ac = gap_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_AC
                && r.end.get() >= CLOSED_AC
        })
        .expect("covering A/C");
    let covering_tg = tg_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= TARGET
                && r.end.get() >= TARGET
        })
        .expect("tg covering");

    let ag_outcome = HaplotypeCallerEngine::call_region(
        covering_ag,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("ag call")
    .expect("ag outcome");
    let ac_outcome = HaplotypeCallerEngine::call_region(
        covering_ac,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("ac call")
    .expect("ac outcome");
    let tg_outcome = HaplotypeCallerEngine::call_region(
        covering_tg,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("tg call")
    .expect("tg outcome");

    let closed_ag = ag_outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(CLOSED_AG)
                && c.event.ref_allele == "A"
                && c.event.alt_allele == "G"
        })
        .expect("A/G");
    let closed_ac = ac_outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(CLOSED_AC)
                && c.event.ref_allele == "A"
                && c.event.alt_allele == "C"
        })
        .expect("A/C");
    let closed_tg = tg_outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(CLOSED_TG)
                && c.event.ref_allele == "T"
                && c.event.alt_allele == "G"
        })
        .expect("T/G");
    let target = tg_outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(TARGET)
                && c.event.ref_allele == "CT"
                && c.event.alt_allele == "C"
        })
        .expect("CT/C");
    let closed_indel = tg_outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(CLOSED_INDEL)
                && c.event.ref_allele == "TTC"
                && c.event.alt_allele == "T"
        })
        .expect("TTC/T");

    let ag_n = unique_indices(&closed_ag.annotation_likelihoods).len();
    let ac_n = unique_indices(&closed_ac.annotation_likelihoods).len();
    let tg_n = unique_indices(&closed_tg.annotation_likelihoods).len();
    kv("closed_92305635_annotation_n", ag_n.to_string());
    kv("closed_92305716_annotation_n", ac_n.to_string());
    kv("closed_92307333_annotation_n", tg_n.to_string());
    assert_eq!(ag_n, 3, "6R.189 n=3 must not change");
    assert_eq!(ac_n, 4, "6R.191 annotation object must not change");
    assert_eq!(tg_n, 1, "6R.186 n=1 must not change");
    assert_eq!(closed_ag.genotype.format.ad_as_i32(), vec![0, 2]);
    assert_eq!(closed_ac.genotype.format.ad_as_i32(), vec![0, 3]);
    assert_eq!(closed_indel.genotype.format.ad_as_i32(), vec![0, 1]);
    assert_eq!(target.genotype.format.ad_as_i32(), vec![1, 1]);
    assert_eq!(target.genotype.format.dp.as_i32(), 2);

    let (full_ref, full_pad) = tg_outcome.assembly.event_map_reference();
    let apply_pad = tg_outcome
        .assembly
        .haplotypes
        .iter()
        .find(|h| h.is_reference)
        .and_then(|h| h.genome_loc.as_ref().map(|g| g.start.get()))
        .unwrap_or(full_pad);
    let likelihoods = if target.annotation_likelihoods.is_empty() {
        tg_outcome.read_likelihoods.as_slice()
    } else {
        target.annotation_likelihoods.as_slice()
    };
    let evidence = HcStrandBiasLikelihoods {
        reads: &tg_outcome.genotyping_reads,
        likelihoods,
        haplotypes: &tg_outcome.assembly.haplotypes,
        contig: covering_tg.contig.as_str(),
        ref_bytes: tg_outcome.assembly.reference_bases(),
        pad_start_1based: apply_pad,
        full_ref_bytes: full_ref,
        full_pad_1based: full_pad,
        max_mnp_distance: tg_outcome.assembly.max_mnp_distance(),
        emit_spanning_dels: true,
    };
    let (ref_bq, alt_bq) = baseq_rank_sum_quals_from_likelihoods(&evidence, TARGET, "CT", "C");
    let (ref_rp, alt_rp) = read_pos_rank_sum_quals_from_likelihoods(&evidence, TARGET, "CT", "C");
    kv(
        "ranksum_lists",
        format!(
            "BaseQ REF={} ALT={} ReadPos REF={} ALT={}",
            ref_bq.len(),
            alt_bq.len(),
            ref_rp.len(),
            alt_rp.len()
        ),
    );
    assert_eq!(ref_bq.len(), 1, "Java BaseQ REF evidence count is 1");
    assert_eq!(alt_bq.len(), 1, "Java BaseQ ALT evidence count is 1");
    assert_eq!(ref_rp.len(), 1);
    assert_eq!(alt_rp.len(), 1);
    assert_eq!(
        ref_bq, alt_bq,
        "equal singleton base qualities are the Java z=0.000 input"
    );
    kv("baseq_option", "Some(0.0) from equal REF/ALT lists");

    let ag_emitted = try_emit_call_region_variants(
        covering_ag,
        &ag_outcome,
        "NA12878",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .expect("ag emit");
    let ac_emitted = try_emit_call_region_variants(
        covering_ac,
        &ac_outcome,
        "NA12878",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .expect("ac emit");
    let tg_emitted = try_emit_call_region_variants(
        covering_tg,
        &tg_outcome,
        "NA12878",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .expect("tg emit");
    let gt_rec = ag_emitted
        .iter()
        .find(|r| r.position == CLOSED_GT && r.reference == "G")
        .expect("G/T emit");
    let ag_rec = ag_emitted
        .iter()
        .find(|r| r.position == CLOSED_AG && r.reference == "A")
        .expect("A/G emit");
    let ac_rec = ac_emitted
        .iter()
        .find(|r| r.position == CLOSED_AC && r.reference == "A")
        .expect("A/C emit");
    let indel_rec = tg_emitted
        .iter()
        .find(|r| r.position == CLOSED_INDEL && r.reference == "TTC")
        .expect("TTC/T emit");
    let tg_rec = tg_emitted
        .iter()
        .find(|r| r.position == CLOSED_TG && r.reference == "T")
        .expect("T/G emit");
    let target_rec = tg_emitted
        .iter()
        .find(|r| r.position == TARGET && r.reference == "CT")
        .expect("CT/C emit");

    assert_eq!(info_i32(&gt_rec.info, "DP"), Some(3));
    assert_eq!(info_i32(&ag_rec.info, "DP"), Some(3));
    assert!((info_f64(&ag_rec.info, "MQ").unwrap_or(-1.0) - 41.96).abs() < 0.005);
    assert!((info_f64(&ag_rec.info, "SOR").unwrap_or(-1.0) - 0.693).abs() < 0.002);
    assert!(!info_has(&ac_rec.info, "ReadPosRankSum"));
    assert!(!info_has(&ac_rec.info, "BaseQRankSum"));
    assert_eq!(info_i32(&indel_rec.info, "DP"), Some(1));
    assert_eq!(info_i32(&tg_rec.info, "DP"), Some(1));
    assert!((info_f64(&tg_rec.info, "MQ").unwrap_or(-1.0) - 44.0).abs() < 0.005);

    let sample = target_rec.samples.first().expect("sample");
    assert_eq!(
        sample.gt.as_ref().map(|g| g.to_string()).as_deref(),
        Some("0/1")
    );
    assert_eq!(sample.ad.as_deref(), Some(&[1u32, 1][..]));
    assert_eq!(sample.dp.map(|d| d as i32), Some(2));
    assert_eq!(sample.gq.map(|g| g as i32), Some(39));
    assert!((target_rec.quality.unwrap_or(0.0) - 31.64).abs() < 0.05);
    assert!(info_has(&target_rec.info, "ReadPosRankSum"));
    let bq = info_f64(&target_rec.info, "BaseQRankSum").expect("BaseQRankSum inserted");
    assert!(
        bq.abs() < 0.0005,
        "finite zero must emit as BaseQRankSum, got {bq}"
    );
    let mqrs = info_f64(&target_rec.info, "MQRankSum").expect("6R.194 MQRankSum");
    assert_eq!(
        (mqrs * 1000.0).round(),
        -674.0,
        "6R.194 emits MQRankSum=-0.674, got {mqrs}"
    );
    kv(
        "after",
        "BaseQRankSum=0.000 emitted; MQRankSum=-0.674 closed by 6R.194; FORMAT/QUAL/ReadPos unchanged",
    );
    kv("classification", "F — WRONG EMISSION PREDICATE (closed)");

    let mid_specs = parse_intervals_cli_string(&dict, "2:92316200-92316580").expect("mid");
    let mid_walk = traverse_assembly_region_walker(
        &dict,
        &mid_specs,
        &ref_fasta,
        &bam,
        &ReadFilterParams::gatk_standard_hc(),
        &WalkerTraversalConfig::gatk_haplotype_caller_production(100),
    )
    .expect("mid walk");
    let mid_regions = flatten_assembly_regions(&mid_walk);
    let covering_mid = mid_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_MID
                && r.end.get() >= CLOSED_MID
        })
        .expect("covering mid");
    let mid_outcome = HaplotypeCallerEngine::call_region(
        covering_mid,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("mid call")
    .expect("mid outcome");
    let mid_emitted = try_emit_call_region_variants(
        covering_mid,
        &mid_outcome,
        "NA12878",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .expect("mid emit");
    let mid_rec = mid_emitted
        .iter()
        .find(|r| r.position == CLOSED_MID && r.reference == "G")
        .expect("G/A emit");
    assert_eq!(info_i32(&mid_rec.info, "DP"), Some(3));
    assert!((info_f64(&mid_rec.info, "MQ").unwrap_or(-1.0) - 40.25).abs() < 0.005);
    assert!(!info_has(&mid_rec.info, "BaseQRankSum"));
}
