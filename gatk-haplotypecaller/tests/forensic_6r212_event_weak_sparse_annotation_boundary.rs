//! 6R.212: attach Java loc-loop annotation on the weak-sparse het
//! (`event_weak_sparse_het_pl`) FORMAT path at `2:92325268 C/T`.
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r212_event_weak_sparse_annotation_boundary -- --nocapture --test-threads=1
//! HOLDOUT_6R212=1 cargo test -p gatk-haplotypecaller --test holdout_6r212_event_weak_sparse_annotation_boundary -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::hc_genotyping_engine::java_alignment_read_overlaps_interval;
use gatk_haplotypecaller::variant_site_hc_annotations::{
    coverage_evidence_count, mq_rank_sum_quals_from_likelihoods, rms_mapping_quality_raw,
    rms_mapping_quality_sample_mapqs, strand_bias_contingency_table, HcStrandBiasLikelihoods,
    STRAND_ODDS_RATIO_MIN_COUNT,
};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    HaplotypeCallerEngine, ReadFilterParams, RegionReadLikelihood, WalkerTraversalConfig,
    DEFAULT_STAND_EMIT_CONFIDENCE,
};
use rust_htslib::bam::Record;
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
const JAVA_COVERAGE_ANNOTATE: &str = r#"
if (likelihoods == null || likelihoods.evidenceCount() == 0) {
    return Collections.emptyMap();
}
final int depth = likelihoods.evidenceCount();
return Collections.singletonMap(getKeyNames().get(0), String.format("%d", depth));
"#;
const JAVA_PREPARE_ANN: &str = r#"
if (hcArgs.useFilteredReadMapForAnnotations || !configuration.isSampleContaminationPresent()) {
    readAlleleLikelihoodsForAnnotations = readAlleleLikelihoodsForGenotyping;
} else {
    readAlleleLikelihoodsForAnnotations = readHaplotypeLikelihoods.marginalize(alleleMapper);
    readAlleleLikelihoodsForAnnotations.retainEvidence(relevantReadsOverlap::overlaps);
}
readAlleleLikelihoodsForAnnotations.addEvidence(overlappingFilteredReads,0);
"#;
const TG_INTERVAL: &str = "2:92307200-92307550";
const MID_INTERVAL: &str = "2:92316200-92316580";
const MID_B_INTERVAL: &str = "2:92317000-92319000";
const P12_POST_INTERVAL: &str = "2:92318150-92319220";
const HET_TAIL_INTERVAL: &str = "2:92324900-92325400";
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 92_325_268;
const TARGET_REF: &str = "C";
const TARGET_ALT: &str = "T";
const KEEP_A: (&str, u16, u8) = ("H06HDADXX130110:2:1101:10045:25649", 81, 24);
const KEEP_B: (&str, u16, u8) = ("H06JUADXX130110:1:1101:10030:91220", 83, 25);
const KEEP_C: (&str, u16, u8) = ("H06HDADXX130110:1:1101:10049:91610", 83, 34);
const DROP_FWD: (&str, u16, u8) = ("H06HDADXX130110:1:1101:10049:91610", 163, 34);
const CLOSED_TG: u64 = 92_307_333;
const CLOSED_CA: u64 = 92_307_403;
const CLOSED_HOM_ALT: u64 = 92_316_296;
const CLOSED_ONE_READ: u64 = 92_316_416;
const CLOSED_MID_B: u64 = 92_317_399;
const CLOSED_POST: u64 = 92_318_199;
const CLOSED_HET: u64 = 92_325_193;
const CLOSED_SIB: u64 = 92_325_205;
const MARGIN: i32 = 2;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R212\t{key}\t{}", value.as_ref());
}

fn qname(rec: &Record) -> String {
    String::from_utf8_lossy(rec.qname()).into_owned()
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

fn rec_at<'a>(
    emitted: &'a [gatk_core::io::vcf::VcfRecord],
    pos: u64,
    r: &str,
) -> &'a gatk_core::io::vcf::VcfRecord {
    emitted
        .iter()
        .find(|rec| rec.position == pos && rec.reference == r)
        .unwrap_or_else(|| panic!("emit {pos} {r}"))
}

fn weak_arm(early: &str) -> &str {
    early
        .split("if event_weak_sparse_het_pl(&event)")
        .nth(1)
        .expect("weak arm")
        .split("if is_cluster_tc_snp(&event)")
        .next()
        .expect("weak body")
}

fn rec_meta(reads: &[gatk_haplotypecaller::SharedBamRecord], idx: usize) -> (String, u16, u8) {
    let rec: &Record = reads[idx].as_ref();
    (qname(rec), rec.flags(), rec.mapq())
}

#[test]
fn forensic_6r212_event_weak_sparse_annotation_boundary_contract() {
    kv("java_pin", JAVA_PIN);
    kv(
        "production_change",
        "one annotation-object construction arrow on event_weak_sparse_het_pl",
    );
    kv("classification", "ANNOTATION_SOURCE_OBJECT_DIVERGENCE");
    assert!(
        JAVA_COVERAGE_ANNOTATE.contains("likelihoods.evidenceCount()")
            && !JAVA_COVERAGE_ANNOTATE.contains("overlaps"),
        "Java Coverage.annotate is evidenceCount of the passed AlleleLikelihoods"
    );
    assert!(
        JAVA_PREPARE_ANN.contains("readAlleleLikelihoodsForGenotyping")
            && JAVA_PREPARE_ANN.contains("addEvidence(overlappingFilteredReads,0)"),
        "default HC annotation reuses genotyping likelihoods (contamination off)"
    );

    let early = include_str!("../src/hc_genotyping_engine/genotype_site_early_template.rs");
    let emit = include_str!("../src/region_vcf_emit.rs");
    let ann = include_str!("../src/variant_site_hc_annotations.rs");
    let pipe = include_str!("../src/hc_genotyping_engine/genotype_site_pipeline.rs");
    let semantics = include_str!("../src/compatibility/java_hc_site_semantics.rs");

    let arm = weak_arm(early);
    assert!(
        arm.contains("vec![-5.5, 0.0, -2.1]") && arm.contains("genotype_from_java_shaped_gls"),
        "weak-sparse het still uses the shaped FORMAT path PL 55,0,21"
    );
    assert!(
        arm.contains("annotation_likelihoods_from_stored_haplotypes")
            && arm.contains("with_annotation_likelihoods"),
        "6R.212 attaches the stored-haplotype loc-loop object on this lifecycle"
    );
    assert!(
        !arm.contains("annotation_likelihoods_from_stored_unique_evidence"),
        "6R.199 stored-unique must not be used on this branch"
    );
    assert!(
        !arm.contains("92325268")
            && !arm.contains(KEEP_A.0)
            && !arm.contains("FLAG=163")
            && !arm.contains("FLAG=81"),
        "production arm must not QNAME/FLAG/locus-filter rows"
    );

    let het_arm = early
        .split("if is_p12_phase_e_gap_het_event(&event)")
        .nth(1)
        .expect("het arm")
        .split("if is_cluster_downstream_snp(&event)")
        .next()
        .expect("het body");
    assert!(
        het_arm.contains("annotation_likelihoods_from_stored_haplotypes"),
        "6R.211 gap-tail het attach stays"
    );
    let gap_arm = early
        .split("if gap_sparse_shaped_early")
        .nth(1)
        .expect("gap_sparse_shaped_early")
        .split("let (_trim_pileup_ref, trim_pileup_alt)")
        .next()
        .expect("gap body");
    assert!(
        gap_arm.contains("annotation_likelihoods_from_stored_haplotypes"),
        "6R.210 gap-sparse attach stays"
    );
    let retain_fn = early
        .split("fn annotation_likelihoods_from_stored_haplotypes")
        .nth(1)
        .expect("186 helper")
        .split("fn annotation_likelihoods_from_stored_unique_evidence")
        .next()
        .expect("186 body");
    assert!(
        retain_fn.contains("marginalize_rows_to_biallelic_alleles")
            && retain_fn.contains("likelihood_subset_for_event")
            && retain_fn.contains("per_variant_annotation_likelihoods")
            && !retain_fn.contains("per_variant_annotation_likelihoods(&subset, reads)")
            && !retain_fn.contains("92325268"),
        "helper stays general Java marginalize + retainEvidence(±2); 6R.209 skip of QNAME collapse stays"
    );
    assert!(
        pipe.contains("let annotation = annotation_likelihoods_from_stored_haplotypes"),
        "6R.189 SiteScore attach stays"
    );
    assert!(
        semantics.contains("pub fn event_weak_sparse_het_pl")
            && semantics
                .split("pub fn event_weak_sparse_het_pl")
                .nth(1)
                .expect("fn")
                .split("pub fn is_gap_tail_het_event")
                .next()
                .expect("body")
                .contains("DOWNSTREAM_CLUSTER_GRADATION_END"),
        "target is the existing weak-sparse het class, not a new coordinate predicate"
    );
    assert!(
        !early.contains("92325268") && !ann.contains("92325268") && !emit.contains("92325268"),
        "no locus-specific 6R.212 production patch"
    );
    let cov = ann
        .split("pub fn coverage_evidence_count")
        .nth(1)
        .expect("coverage fn")
        .split("fn mq_rms_of_sample_evidence")
        .next()
        .expect("coverage body");
    assert!(
        !cov.contains("java_alignment_read_overlaps_interval"),
        "6R.201 Coverage cardinality stays"
    );
    assert!(
        ann.contains("fn mq_rms_of_sample_evidence")
            && ann.contains("strand_bias_sample_counts")
            && !ann.contains("DP = 3"),
        "MQ/SOR/Coverage formulas are not patched in this round"
    );
    let ge = include_str!("../src/hc_genotyping_engine/mod.rs");
    assert!(
        ge.contains("if qnames.len() == 1 && best_ll.len() > 1"),
        "6R.180 same-QNAME collapse remains on FORMAT-subset / stored-unique callers"
    );
}

#[test]
fn forensic_6r212_event_weak_sparse_annotation_boundary() {
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
        "one annotation-object construction arrow on event_weak_sparse_het_pl",
    );
    kv("target", "2:92325268 C/T");
    kv("classification", "ANNOTATION_SOURCE_OBJECT_DIVERGENCE");

    let java_vcf = root.join("parity/reports/6r43/p12_het_tail/java.vcf");
    if java_vcf.is_file() {
        let text = std::fs::read_to_string(&java_vcf).expect("java.vcf");
        let line = text
            .lines()
            .find(|l| l.contains("\t92325268\t") && l.contains("\tC\tT\t"))
            .expect("frozen Java 4.4.0.0 p12_het_tail 2:92325268 C/T");
        let fields: Vec<&str> = line.split('\t').collect();
        kv("java_vcf_qual", fields[5]);
        kv("java_vcf_info", fields[7]);
        kv("java_vcf_format", format!("{} {}", fields[8], fields[9]));
        assert_eq!(fields[5], "47.64");
        assert!(
            fields[7].contains("DP=3")
                && fields[7].contains("MQ=28.03")
                && fields[7].contains("SOR=1.179")
                && fields[7].contains("MQRankSum=0.000")
        );
        assert!(
            fields[9].contains("0/1") && fields[9].contains("1,2") && fields[9].contains("55,0,21")
        );
    }

    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let walk = |interval: &str| {
        let specs = parse_intervals_cli_string(&dict, interval).expect("interval");
        traverse_assembly_region_walker(
            &dict,
            &specs,
            &ref_fasta,
            &bam,
            &ReadFilterParams::gatk_standard_hc(),
            &WalkerTraversalConfig::gatk_haplotype_caller_production(100),
        )
        .expect("walk")
    };
    let tg_regions = flatten_assembly_regions(&walk(TG_INTERVAL));
    let mid_regions = flatten_assembly_regions(&walk(MID_INTERVAL));
    let mid_b_regions = flatten_assembly_regions(&walk(MID_B_INTERVAL));
    let post_regions = flatten_assembly_regions(&walk(P12_POST_INTERVAL));
    let het_regions = flatten_assembly_regions(&walk(HET_TAIL_INTERVAL));
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
    let covering_target = het_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= TARGET
                && r.end.get() >= TARGET
        })
        .expect("het-tail");

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
    let outcome = call_one(covering_target);

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
    let ca = find(&ca_outcome, CLOSED_CA, "C", "A");
    let hom = find(&hom_outcome, CLOSED_HOM_ALT, "A", "T");
    let one = find(&one_outcome, CLOSED_ONE_READ, "C", "A");
    let mid_b = find(&mid_b_outcome, CLOSED_MID_B, "C", "A");
    let post = find(&post_outcome, CLOSED_POST, "C", "T");
    let het = find(&outcome, CLOSED_HET, "C", "T");
    let sib = find(&outcome, CLOSED_SIB, "G", "A");
    let target = find(&outcome, TARGET, TARGET_REF, TARGET_ALT);

    kv(
        "branch",
        "event_weak_sparse_het_pl attaches annotation_likelihoods_from_stored_haplotypes",
    );

    let stored = unique_indices(&outcome.read_likelihoods);
    let attached = unique_indices(&target.annotation_likelihoods);
    let mut retain = BTreeSet::new();
    for idx in &stored {
        if java_alignment_read_overlaps_interval(
            &outcome.genotyping_reads[*idx],
            TARGET,
            TARGET,
            MARGIN,
        ) {
            retain.insert(*idx);
        }
    }
    kv("stored_hap_n", stored.len().to_string());
    kv("java_loc_loop_retainEvidence_n", retain.len().to_string());
    kv("attached_annotation_n", attached.len().to_string());

    assert_eq!(target.genotype.format.ad_as_i32(), vec![1, 2]);
    assert_eq!(target.genotype.format.dp.as_i32(), 3);
    assert_eq!(target.genotype.format.gq.as_i32(), 21);
    assert_eq!(target.genotype.format.pl_as_i32(), vec![55, 0, 21]);

    assert_eq!(unique_indices(&tg.annotation_likelihoods).len(), 1);
    assert_eq!(unique_indices(&ca.annotation_likelihoods).len(), 6);
    assert_eq!(unique_indices(&hom.annotation_likelihoods).len(), 2);
    assert_eq!(unique_indices(&one.annotation_likelihoods).len(), 1);
    assert_eq!(unique_indices(&mid_b.annotation_likelihoods).len(), 2);
    assert_eq!(unique_indices(&post.annotation_likelihoods).len(), 1);
    assert_eq!(unique_indices(&het.annotation_likelihoods).len(), 3);
    assert_eq!(unique_indices(&sib.annotation_likelihoods).len(), 3);

    assert_eq!(stored.len(), 4, "het-tail stored unique n=4");
    assert_eq!(retain.len(), 3, "Java loc-loop retainEvidence(±2) n=3");
    assert_eq!(
        attached.len(),
        3,
        "6R.212: loc-loop annotation is retainEvidence n=3"
    );
    assert_eq!(attached, retain);
    assert!(!target.annotation_likelihoods.is_empty());

    let attached_meta: Vec<_> = attached
        .iter()
        .copied()
        .map(|i| rec_meta(&outcome.genotyping_reads, i))
        .collect();
    let has = |want: (&str, u16, u8)| {
        attached_meta
            .iter()
            .any(|m| m.0 == want.0 && m.1 == want.1 && m.2 == want.2)
    };
    assert!(has(KEEP_A) && has(KEEP_B) && has(KEEP_C));
    assert!(
        !attached_meta
            .iter()
            .any(|m| m.0 == DROP_FWD.0 && m.1 == DROP_FWD.1 && m.2 == DROP_FWD.2),
        "FLAG=163 that misses overlap must not survive loc-loop retainEvidence"
    );
    kv(
        "dropped_fwd",
        format!(
            "{} FLAG={} MAPQ={} start=92325173 end=92325230 overlap=false",
            DROP_FWD.0, DROP_FWD.1, DROP_FWD.2
        ),
    );

    let attached_dp = coverage_evidence_count(
        &outcome.genotyping_reads,
        &target.annotation_likelihoods,
        TARGET,
        TARGET,
        MARGIN,
    );
    let stored_dp = coverage_evidence_count(
        &outcome.genotyping_reads,
        &outcome.read_likelihoods,
        TARGET,
        TARGET,
        MARGIN,
    );
    kv("coverage_attached_n", attached_dp.to_string());
    kv("coverage_stored_unique_n", stored_dp.to_string());
    assert_eq!(attached_dp, 3);
    assert_eq!(stored_dp, 4);
    assert_eq!(attached_dp, attached.len() as i32);

    let attached_mqs =
        rms_mapping_quality_sample_mapqs(&outcome.genotyping_reads, &target.annotation_likelihoods);
    kv("attached_mapq_list", format!("{attached_mqs:?}"));
    assert_eq!(attached_mqs, vec![24, 25, 34]);
    let mq_rms = rms_mapping_quality_raw(&attached_mqs).map(|t| (t.2 * 100.0).round() / 100.0);
    assert_eq!(mq_rms, Some(28.03));

    let (full_ref, full_pad) = outcome.assembly.event_map_reference();
    let apply_pad = outcome
        .assembly
        .haplotypes
        .iter()
        .find(|h| h.is_reference)
        .and_then(|h| h.genome_loc.as_ref().map(|g| g.start.get()))
        .unwrap_or(full_pad);
    let attached_evidence = HcStrandBiasLikelihoods {
        reads: &outcome.genotyping_reads,
        likelihoods: &target.annotation_likelihoods,
        haplotypes: &outcome.assembly.haplotypes,
        contig: &covering_target.contig,
        ref_bytes: outcome.assembly.reference_bases(),
        pad_start_1based: apply_pad,
        full_ref_bytes: full_ref,
        full_pad_1based: full_pad,
        max_mnp_distance: outcome.assembly.max_mnp_distance(),
        emit_spanning_dels: true,
    };
    let attached_table = strand_bias_contingency_table(
        &attached_evidence,
        TARGET,
        TARGET_REF,
        TARGET_ALT,
        STRAND_ODDS_RATIO_MIN_COUNT,
    );
    kv(
        "attached_sor_table",
        format!(
            "[{},{};{},{}]",
            attached_table.0, attached_table.1, attached_table.2, attached_table.3
        ),
    );
    assert_eq!(attached_table, (0, 1, 0, 2));
    let (ref_mq, alt_mq) =
        mq_rank_sum_quals_from_likelihoods(&attached_evidence, TARGET, TARGET_REF, TARGET_ALT);
    kv("attached_mqrs_ref", format!("{ref_mq:?}"));
    kv("attached_mqrs_alt", format!("{alt_mq:?}"));
    assert_eq!(ref_mq, vec![25.0]);
    assert_eq!(alt_mq, vec![24.0, 34.0]);

    let emit = |region: &gatk_haplotypecaller::AssemblyRegion,
                o: &gatk_haplotypecaller::CallRegionOutcome| {
        try_emit_call_region_variants(region, o, "NA12878", DEFAULT_STAND_EMIT_CONFIDENCE)
            .expect("emit")
    };
    let tg_emitted = emit(covering_tg, &tg_outcome);
    let ca_emitted = emit(covering_ca, &ca_outcome);
    let hom_emitted = emit(covering_hom, &hom_outcome);
    let one_emitted = emit(covering_one, &one_outcome);
    let mid_b_emitted = emit(covering_mid_b, &mid_b_outcome);
    let post_emitted = emit(covering_post, &post_outcome);
    let target_emitted = emit(covering_target, &outcome);

    let tg_rec = rec_at(&tg_emitted, CLOSED_TG, "T");
    assert_eq!(info_i32(&tg_rec.info, "DP"), Some(1));
    assert!((info_f64(&tg_rec.info, "MQ").expect("mq") - 44.00).abs() < 0.005);
    let ca_rec = rec_at(&ca_emitted, CLOSED_CA, "C");
    assert_eq!(info_i32(&ca_rec.info, "DP"), Some(6));
    assert!((info_f64(&ca_rec.info, "MQ").expect("mq") - 40.58).abs() < 0.005);
    let hom_rec = rec_at(&hom_emitted, CLOSED_HOM_ALT, "A");
    assert_eq!(info_i32(&hom_rec.info, "DP"), Some(2));
    assert!((info_f64(&hom_rec.info, "MQ").expect("mq") - 47.00).abs() < 0.005);
    let one_rec = rec_at(&one_emitted, CLOSED_ONE_READ, "C");
    assert_eq!(info_i32(&one_rec.info, "DP"), Some(1));
    assert!((info_f64(&one_rec.info, "MQ").expect("mq") - 21.00).abs() < 0.005);
    let mid_b_rec = rec_at(&mid_b_emitted, CLOSED_MID_B, "C");
    assert_eq!(info_i32(&mid_b_rec.info, "DP"), Some(2));
    assert!((info_f64(&mid_b_rec.info, "MQ").expect("mq") - 27.00).abs() < 0.005);
    let post_rec = rec_at(&post_emitted, CLOSED_POST, "C");
    assert_eq!(info_i32(&post_rec.info, "DP"), Some(1));
    assert!((info_f64(&post_rec.info, "MQ").expect("mq") - 24.00).abs() < 0.005);
    let het_rec = rec_at(&target_emitted, CLOSED_HET, "C");
    assert_eq!(info_i32(&het_rec.info, "DP"), Some(3));
    assert!((info_f64(&het_rec.info, "MQ").expect("mq") - 28.03).abs() < 0.005);
    let sib_rec = rec_at(&target_emitted, CLOSED_SIB, "G");
    assert_eq!(info_i32(&sib_rec.info, "DP"), Some(3));
    assert!((info_f64(&sib_rec.info, "MQ").expect("mq") - 28.03).abs() < 0.005);

    let rec = rec_at(&target_emitted, TARGET, TARGET_REF);
    assert!((rec.quality.expect("QUAL") - 47.64).abs() < 0.005);
    let dp = info_i32(&rec.info, "DP").expect("DP");
    let mq = info_f64(&rec.info, "MQ").expect("MQ");
    let sor = info_f64(&rec.info, "SOR").expect("SOR");
    kv("rust_info_dp", dp.to_string());
    kv("java_info_dp", "3");
    kv("rust_mq", format!("{mq:.2}"));
    kv("java_mq", "28.03");
    kv("rust_sor", format!("{sor:.3}"));
    kv("java_sor", "1.179");
    let sample = rec.samples.first().expect("sample");
    assert_eq!(
        sample.gt.as_ref().map(|g| g.to_string()).as_deref(),
        Some("0/1")
    );
    assert_eq!(sample.ad.as_deref(), Some(&[1u32, 2][..]));
    assert_eq!(sample.dp, Some(3));
    assert_eq!(sample.gq, Some(21.0));
    assert_eq!(sample.pl.as_deref(), Some(&[55u32, 0, 21][..]));
    assert_eq!(dp, 3);
    assert!((mq - 28.03).abs() < 0.005);
    assert!((sor - 1.179).abs() < 0.002);
    let mqrs = info_f64(&rec.info, "MQRankSum").expect("MQRankSum still emitted");
    assert_eq!((mqrs * 1000.0).round(), 0.0);

    kv(
        "first_arrow",
        "event_weak_sparse_het_pl FORMAT-only return left annotation empty; emit used region-wide n=4 including non-overlapping FLAG=163",
    );
}
