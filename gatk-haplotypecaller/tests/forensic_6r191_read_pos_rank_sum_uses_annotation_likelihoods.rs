//! 6R.191: ReadPosRankSum consumes the per-variant annotation AlleleLikelihoods
//! (`fillQualsFromLikelihood`), not `region.reads` pileup.
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! Target `2:92305716 A/C`. Formula/undefined semantics stay 6R.172.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r191_read_pos_rank_sum_uses_annotation_likelihoods -- --nocapture --test-threads=1
//! HOLDOUT_6R191=1 cargo test -p gatk-haplotypecaller --test holdout_6r191_read_pos_rank_sum -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::variant_site_hc_annotations::{
    read_pos_rank_sum_quals_from_likelihoods, HcStrandBiasLikelihoods,
};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    HaplotypeCallerEngine, HcGenotypingConfig, ReadFilterParams, RegionReadLikelihood,
    WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
const INTERVAL: &str = "2:92305500-92305850";
const TG_INTERVAL: &str = "2:92307200-92307550";
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 92_305_716;
const MERGED_REF: &str = "A";
const MERGED_ALT: &str = "C";
const CLOSED_SNP: u64 = 92_305_634;
const CLOSED_AG: u64 = 92_305_635;
const CLOSED_INDEL: u64 = 92_307_324;
const CLOSED_TG: u64 = 92_307_333;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R191\t{key}\t{}", value.as_ref());
}

fn unique_indices(likelihoods: &[RegionReadLikelihood]) -> BTreeSet<usize> {
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
fn forensic_6r191_source_contracts() {
    kv("java_pin", JAVA_PIN);
    kv(
        "production_change",
        "ReadPosRankSum fillQualsFromLikelihood on annotation AlleleLikelihoods",
    );
    kv("classification", "A — WRONG SOURCE OBJECT (closed)");
    kv("target", "2:92305716 A/C");

    let ann = include_str!("../src/variant_site_hc_annotations.rs");
    let plugin = include_str!("../src/annotator/plugins/read_pos_rank_sum.rs");
    let emit = include_str!("../src/region_vcf_emit.rs");
    let pipe = include_str!("../src/hc_genotyping_engine/genotype_site_pipeline.rs");
    let early = include_str!("../src/hc_genotyping_engine/genotype_site_early_template.rs");
    let production_plugin = plugin.split("#[cfg(test)]").next().expect("plugin");

    let helper = ann
        .split("fn fill_ranksum_quals_from_likelihood")
        .nth(1)
        .expect("6R.191 helper")
        .split("pub fn rms_mapping_quality_sample_reads")
        .next()
        .expect("helper body");
    assert!(
        helper.contains("evidence.likelihoods")
            && helper.contains("marginalize_rows_to_biallelic_alleles")
            && helper.contains("LOG_10_INFORMATIVE_THRESHOLD")
            && helper.contains("java_read_pos_rank_element"),
        "helper is fillQualsFromLikelihood on evidence.likelihoods"
    );
    assert!(
        !helper.contains("region.reads") && !helper.contains("read_offset_evidence_at_site"),
        "helper must not reconstruct pileup"
    );
    assert!(
        !helper.contains("92305716") && !helper.contains("92305635"),
        "helper has no locus special case"
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
        annotate_body.contains("read_pos_rank_sum_quals_from_likelihoods(")
            && annotate_body.contains(
                "let rp = read_pos_rank_sum::read_pos_rank_sum(&ref_positions, &alt_positions)"
            ),
        "annotate binds RankSum lists from strand_likelihoods, then existing rank_sum_z"
    );
    assert!(
        !annotate_body.contains("read_offset_evidence_at_site(region,"),
        "annotate must not still pileup RankSum from region.reads"
    );

    assert!(
        production_plugin.contains(
            "pub fn read_pos_rank_sum(ref_positions: &[f64], alt_positions: &[f64]) -> Option<f64>"
        ) && production_plugin.contains("if alt.is_empty() || reference.is_empty()")
            && production_plugin.contains("return None")
            && production_plugin.contains("result.z.is_nan()"),
        "6R.172 formula/undefined semantics unchanged"
    );
    assert!(
        emit.contains("if let Some(z) = ann.read_pos_rank_sum")
            && emit.contains("likelihoods: &call.annotation_likelihoods"),
        "emit still inserts Some(z) and still binds the 6R.189 attached object"
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
        info_body.contains("BaseQRankSum") && info_body.contains("InfoValue::Float(\"MQRankSum\""),
        "6R.193 emits BaseQRankSum; 6R.194 emits MQRankSum"
    );
    assert!(
        pipe.contains("let annotation = annotation_likelihoods_from_stored_haplotypes")
            && early.contains("if is_cluster_tg_snp(&event)"),
        "6R.189 SiteScore attach and 6R.186 cluster-TG stay"
    );
    assert!(
        ann.contains("6R.174: INFO DP is Java `Coverage.evidenceCount`")
            && ann.contains("mq_rms_of_sample_evidence")
            && ann.contains("strand_bias_sample_counts"),
        "DP/MQ/SOR formulas must stay closed"
    );
    assert!(
        !ann.contains("92305716") && !emit.contains("92305716") && !pipe.contains("92305716"),
        "no locus-specific 6R.191 production patch"
    );
}

#[test]
fn forensic_6r191_read_pos_rank_sum_uses_annotation_likelihoods() {
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing P12 ref/BAM");
        return;
    }
    kv("java_pin", JAVA_PIN);
    kv("target", "2:92305716 A/C");
    kv(
        "before",
        "annotation n=4; pileup REF=2 ALT=3; ReadPosRankSum=1.645",
    );

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
    let covering_closed = regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_AG
                && r.end.get() >= CLOSED_SNP
        })
        .expect("covering closed");
    let covering = regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= TARGET
                && r.end.get() >= TARGET
        })
        .expect("covering");
    let closed_outcome = HaplotypeCallerEngine::call_region(
        covering_closed,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("closed call")
    .expect("closed outcome");
    let outcome = HaplotypeCallerEngine::call_region(
        covering,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("outcome");

    let closed_gt = closed_outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(CLOSED_SNP)
                && c.event.ref_allele == "G"
                && c.event.alt_allele == "T"
        })
        .expect("G/T");
    let closed_ag = closed_outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(CLOSED_AG)
                && c.event.ref_allele == "A"
                && c.event.alt_allele == "G"
        })
        .expect("A/G");
    let target = outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(TARGET)
                && c.event.ref_allele == MERGED_REF
                && c.event.alt_allele == MERGED_ALT
        })
        .expect("A/C");

    let closed_ag_ann = unique_indices(&closed_ag.annotation_likelihoods);
    let target_ann = unique_indices(&target.annotation_likelihoods);
    kv(
        "closed_92305635_annotation_n",
        closed_ag_ann.len().to_string(),
    );
    kv("target_annotation_n", target_ann.len().to_string());
    assert_eq!(closed_ag_ann.len(), 3, "6R.189 n=3 must not change");
    assert_eq!(
        target_ann.len(),
        4,
        "annotation object itself must not change"
    );
    assert_eq!(closed_ag.genotype.format.ad_as_i32(), vec![0, 2]);
    assert_eq!(closed_ag.genotype.format.dp.as_i32(), 2);
    assert_eq!(closed_ag.genotype.format.gq.as_i32(), 6);
    assert_eq!(closed_ag.genotype.format.pl_as_i32(), vec![90, 6, 0]);
    assert_eq!(target.genotype.format.ad_as_i32(), vec![0, 3]);
    assert_eq!(target.genotype.format.dp.as_i32(), 3);
    assert_eq!(target.genotype.format.pl_as_i32(), vec![130, 9, 0]);
    assert_eq!(closed_gt.genotype.format.dp.as_i32(), 2);

    let (full_ref, full_pad) = outcome.assembly.event_map_reference();
    let apply_pad = outcome
        .assembly
        .haplotypes
        .iter()
        .find(|h| h.is_reference)
        .and_then(|h| h.genome_loc.as_ref().map(|g| g.start.get()))
        .unwrap_or(full_pad);
    let cfg = HcGenotypingConfig::default();
    let evidence = HcStrandBiasLikelihoods {
        reads: &outcome.genotyping_reads,
        likelihoods: &target.annotation_likelihoods,
        haplotypes: &outcome.assembly.haplotypes,
        contig: &covering.contig,
        ref_bytes: outcome.assembly.reference_bases(),
        pad_start_1based: apply_pad,
        full_ref_bytes: full_ref,
        full_pad_1based: full_pad,
        max_mnp_distance: outcome.assembly.max_mnp_distance(),
        emit_spanning_dels: !cfg.disable_spanning_event_genotyping,
    };
    kv(
        "ranksum_source",
        "GenotypedSiteCall.annotation_likelihoods via read_pos_rank_sum_quals_from_likelihoods",
    );
    let (java_ref, java_alt) =
        read_pos_rank_sum_quals_from_likelihoods(&evidence, TARGET, MERGED_REF, MERGED_ALT);
    kv(
        "likelihood_counts",
        format!("REF={} ALT={}", java_ref.len(), java_alt.len()),
    );
    kv("likelihood_ref_pos", format!("{java_ref:?}"));
    kv("likelihood_alt_pos", format!("{java_alt:?}"));
    assert_eq!(java_ref.len(), 0, "Java-equivalent REF list is empty");
    assert_eq!(java_alt.len(), 3, "Java-equivalent ALT list is 3");

    let closed_emitted = try_emit_call_region_variants(
        covering_closed,
        &closed_outcome,
        "NA12878",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .expect("closed emit");
    let emitted =
        try_emit_call_region_variants(covering, &outcome, "NA12878", DEFAULT_STAND_EMIT_CONFIDENCE)
            .expect("emit");
    let closed_gt_rec = closed_emitted
        .iter()
        .find(|r| r.position == CLOSED_SNP && r.reference == "G")
        .expect("G/T emit");
    let closed_ag_rec = closed_emitted
        .iter()
        .find(|r| r.position == CLOSED_AG && r.reference == "A")
        .expect("A/G emit");
    let target_rec = emitted
        .iter()
        .find(|r| r.position == TARGET && r.reference == MERGED_REF)
        .expect("A/C emit");
    let ag_sample = closed_ag_rec.samples.first().expect("sample");
    let target_sample = target_rec.samples.first().expect("sample");
    assert_eq!(
        ag_sample.gt.as_ref().map(|g| g.to_string()).as_deref(),
        Some("1/1")
    );
    assert_eq!(ag_sample.ad.as_deref(), Some(&[0u32, 2][..]));
    assert_eq!(ag_sample.dp.map(|d| d as i32), Some(2));
    assert!((closed_ag_rec.quality.unwrap_or(0.0) - 78.32).abs() < 0.05);
    assert_eq!(info_i32(&closed_gt_rec.info, "DP"), Some(3));
    assert_eq!(info_i32(&closed_ag_rec.info, "DP"), Some(3));
    let ag_mq = info_f64(&closed_ag_rec.info, "MQ").unwrap_or(-1.0);
    let ag_sor = info_f64(&closed_ag_rec.info, "SOR").unwrap_or(-1.0);
    assert!((ag_mq - 41.96).abs() < 0.005);
    assert!((ag_sor - 0.693).abs() < 0.002);
    assert_eq!(
        target_sample.gt.as_ref().map(|g| g.to_string()).as_deref(),
        Some("1/1")
    );
    assert_eq!(target_sample.ad.as_deref(), Some(&[0u32, 3][..]));
    assert_eq!(target_sample.dp.map(|d| d as i32), Some(3));
    assert!((target_rec.quality.unwrap_or(0.0) - 116.84).abs() < 0.05);
    assert_eq!(info_i32(&target_rec.info, "DP"), Some(4));
    let target_mq = info_f64(&target_rec.info, "MQ").unwrap_or(-1.0);
    assert!((target_mq - 35.15).abs() < 0.005);
    assert!(
        !info_has(&target_rec.info, "ReadPosRankSum"),
        "empty REF → undefined → omit"
    );
    assert!(!info_has(&target_rec.info, "BaseQRankSum"));
    assert!(!info_has(&target_rec.info, "MQRankSum"));
    kv(
        "after",
        "annotation n=4; likelihoods REF=0 ALT=3; ReadPosRankSum absent",
    );

    let tg_specs = parse_intervals_cli_string(&dict, TG_INTERVAL).expect("tg");
    let tg_walk = traverse_assembly_region_walker(
        &dict,
        &tg_specs,
        &ref_fasta,
        &bam,
        &ReadFilterParams::gatk_standard_hc(),
        &WalkerTraversalConfig::gatk_haplotype_caller_production(100),
    )
    .expect("tg walk");
    let tg_regions = flatten_assembly_regions(&tg_walk);
    let tg_covering = tg_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_TG
                && r.end.get() >= CLOSED_INDEL
        })
        .expect("tg covering");
    let tg_outcome = HaplotypeCallerEngine::call_region(
        tg_covering,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("tg call")
    .expect("tg outcome");
    let tg_call = tg_outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(CLOSED_TG)
                && c.event.ref_allele == "T"
                && c.event.alt_allele == "G"
        })
        .expect("T/G");
    let tg_ann = unique_indices(&tg_call.annotation_likelihoods);
    assert_eq!(tg_ann.len(), 1, "6R.186 n=1 stays");
    let tg_emitted = try_emit_call_region_variants(
        tg_covering,
        &tg_outcome,
        "NA12878",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .expect("tg emit");
    let indel = tg_emitted
        .iter()
        .find(|r| r.position == CLOSED_INDEL && r.reference == "TTC")
        .expect("TTC/T");
    let tg_rec = tg_emitted
        .iter()
        .find(|r| r.position == CLOSED_TG && r.reference == "T")
        .expect("T/G emit");
    assert_eq!(info_i32(&indel.info, "DP"), Some(1));
    assert_eq!(info_i32(&tg_rec.info, "DP"), Some(1));
    let tg_mq = info_f64(&tg_rec.info, "MQ").unwrap_or(-1.0);
    let tg_sor = info_f64(&tg_rec.info, "SOR").unwrap_or(-1.0);
    assert!((tg_mq - 44.0).abs() < 0.005);
    assert!((tg_sor - 1.609).abs() < 0.002);
}
