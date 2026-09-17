//! 6R.199: cluster-downstream early-template attaches the Java-equivalent
//! loc-loop annotation object from stored unique haplotype evidence.
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! Target `2:92307403 C/A`. Does **not** reuse the 6R.186 retainEvidence subset.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r199_cluster_downstream_annotation_object -- --nocapture --test-threads=1
//! HOLDOUT_6R199=1 cargo test -p gatk-haplotypecaller --test holdout_6r199_cluster_downstream_annotation -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::hc_genotyping_engine::java_alignment_read_overlaps_interval;
use gatk_haplotypecaller::variant_site_hc_annotations::{
    coverage_evidence_count, rms_mapping_quality_raw, rms_mapping_quality_sample_mapqs,
};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegion, AssemblyRegionCallDisposition, CallRegionArgs,
    GenomePosition, HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
    DEFAULT_STAND_EMIT_CONFIDENCE,
};
use rust_htslib::bam::Record;
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
const GAP_INTERVAL: &str = "2:92305500-92305850";
const TG_INTERVAL: &str = "2:92307200-92307550";
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 92_307_403;
const TARGET_REF: &str = "C";
const TARGET_ALT: &str = "A";
const CLOSED_GT: u64 = 92_305_634;
const CLOSED_AG: u64 = 92_305_635;
const CLOSED_AC: u64 = 92_305_716;
const CLOSED_INDEL: u64 = 92_307_324;
const CLOSED_TG: u64 = 92_307_333;
const CLOSED_QUAL: u64 = 92_307_359;
const MARGIN: i32 = 2;
const EXTRA_QNAME_A: &str = "H06HDADXX130110:1:1101:10061:17286";
const EXTRA_FLAG_A: u16 = 83;
const EXTRA_QNAME_B: &str = "H06JUADXX130110:1:1101:10011:51168";
const EXTRA_FLAG_B: u16 = 81;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R199\t{key}\t{}", value.as_ref());
}

fn qname(rec: &Record) -> String {
    String::from_utf8_lossy(rec.qname()).into_owned()
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

fn format_info(info: &[InfoValue]) -> String {
    let mut parts = Vec::new();
    for v in info {
        parts.push(match v {
            InfoValue::Flag(k) => k.clone(),
            InfoValue::Integer(k, xs) => format!(
                "{k}={}",
                xs.iter()
                    .map(|x| x.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            ),
            InfoValue::Float(k, xs) => format!(
                "{k}={}",
                xs.iter()
                    .map(|x| format!("{x}"))
                    .collect::<Vec<_>>()
                    .join(",")
            ),
            InfoValue::String(k, xs) => format!("{k}={}", xs.join(",")),
            InfoValue::Character(k, xs) => format!(
                "{k}={}",
                xs.iter()
                    .map(|c| c.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            ),
        });
    }
    parts.join(";")
}

fn has_member(
    reads: &[gatk_haplotypecaller::SharedBamRecord],
    idxs: &BTreeSet<usize>,
    name: &str,
    flag: u16,
) -> bool {
    idxs.iter().any(|idx| {
        let rec: &Record = reads[*idx].as_ref();
        rec.flags() == flag && qname(rec) == name
    })
}

#[test]
fn forensic_6r199_source_contracts() {
    kv("java_pin", JAVA_PIN);
    kv(
        "production_change",
        "cluster-downstream early-template attaches stored unique hap evidence",
    );
    let early = include_str!("../src/hc_genotyping_engine/genotype_site_early_template.rs");
    let fin = include_str!("../src/hc_genotyping_engine/genotype_finalize.rs");
    let pipe = include_str!("../src/hc_genotyping_engine/genotype_site_pipeline.rs");
    let emit = include_str!("../src/region_vcf_emit.rs");
    let ann = include_str!("../src/variant_site_hc_annotations.rs");
    let ge = include_str!("../src/hc_genotyping_engine/mod.rs");

    let ds_arm = early
        .split("if is_cluster_downstream_snp(&event)")
        .nth(1)
        .expect("downstream arm");
    let ds_body = ds_arm
        .split("if is_cluster_tg_snp(&event)")
        .next()
        .expect("downstream body");
    assert!(
        ds_body.contains("java_cluster_downstream_shaped_genotype")
            && ds_body.contains("finish_strict_java_shaped_site_call")
            && ds_body.contains("annotation_likelihoods_from_stored_unique_evidence")
            && ds_body.contains("with_annotation_likelihoods")
            && !ds_body.contains("annotation_likelihoods_from_stored_haplotypes")
            && !ds_body.contains("likelihood_subset_for_event"),
        "cluster-downstream attaches stored unique evidence, not 6R.186 retain"
    );

    let unique_fn = early
        .split("fn annotation_likelihoods_from_stored_unique_evidence")
        .nth(1)
        .expect("unique helper");
    assert!(
        unique_fn.contains("per_variant_annotation_likelihoods(likelihoods, reads)")
            && !unique_fn.contains("likelihood_subset_for_event")
            && !unique_fn.contains("92307403")
            && !unique_fn.contains(EXTRA_QNAME_A)
            && !unique_fn.contains(EXTRA_QNAME_B),
        "unique helper reuses stored hap cells; no retainEvidence / locus / QNAME pin"
    );

    let retain_fn = early
        .split("fn annotation_likelihoods_from_stored_haplotypes")
        .nth(1)
        .expect("186 helper")
        .split("fn annotation_likelihoods_from_stored_unique_evidence")
        .next()
        .expect("186 body");
    assert!(
        retain_fn.contains("likelihood_subset_for_event"),
        "6R.186 retain helper stays unchanged"
    );

    let tg_arm = early
        .split("if is_cluster_tg_snp(&event)")
        .nth(1)
        .expect("tg arm");
    let tg_body = tg_arm
        .split("if is_cluster_tc_snp(&event)")
        .next()
        .expect("tg body");
    assert!(
        tg_body.contains("annotation_likelihoods_from_stored_haplotypes")
            && tg_body.contains("with_annotation_likelihoods"),
        "6R.186 cluster-TG attach stays"
    );
    assert!(
        pipe.contains("let annotation = annotation_likelihoods_from_stored_haplotypes")
            && pipe.contains("with_annotation_likelihoods(annotation)"),
        "6R.189 SiteScore attach stays"
    );

    let fin_fn = fin
        .split("fn finish_strict_java_shaped_site_call")
        .nth(1)
        .expect("finish fn");
    let fin_body = fin_fn
        .split("fn finalize_strict_java_variation_genotype")
        .next()
        .expect("finish body");
    assert!(
        !fin_body.contains("with_annotation_likelihoods"),
        "shared finish helper does not globally attach annotation_likelihoods"
    );
    assert!(
        emit.contains("if call.annotation_likelihoods.is_empty()")
            && emit.contains("likelihoods: region_wide.likelihoods"),
        "6R.180 empty-object fallback is unchanged"
    );
    assert!(
        ann.contains("coverage_evidence_count")
            && ge.contains("fn per_variant_annotation_likelihoods"),
        "Coverage formula and unique-evidence helper stay"
    );
    assert!(
        !early.contains("92307403") && !pipe.contains("92307403") && !ann.contains("92307403"),
        "no locus-specific 6R.199 production patch"
    );
}

#[test]
fn forensic_6r199_cluster_downstream_annotation_object() {
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
        "cluster-downstream early-template attaches stored unique hap evidence",
    );
    kv("target", "2:92307403 C/A");
    kv(
        "classification",
        "E — DIFFERENT BRANCH / SPECIAL-PATH SEMANTICS",
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
    let covering_gap = |pos: u64| {
        gap_regions
            .iter()
            .find(|r| {
                matches!(
                    call_disposition(r),
                    AssemblyRegionCallDisposition::ActiveFull
                ) && r.start.get() <= pos
                    && r.end.get() >= pos
            })
            .unwrap_or_else(|| panic!("gap covering {pos}"))
    };
    let covering_tg_int = |pos: u64| {
        tg_regions
            .iter()
            .find(|r| {
                matches!(
                    call_disposition(r),
                    AssemblyRegionCallDisposition::ActiveFull
                ) && r.start.get() <= pos
                    && r.end.get() >= pos
            })
            .unwrap_or_else(|| panic!("tg covering {pos}"))
    };
    let covering_gt = covering_gap(CLOSED_GT);
    let covering_ag = covering_gap(CLOSED_AG);
    let covering_ac = covering_gap(CLOSED_AC);
    let covering_indel = covering_tg_int(CLOSED_INDEL);
    let covering_tg = covering_tg_int(CLOSED_TG);
    let covering_qual = covering_tg_int(CLOSED_QUAL);
    let covering = covering_tg_int(TARGET);

    let call = |region: &AssemblyRegion| {
        HaplotypeCallerEngine::call_region(
            region,
            &dict,
            &ref_fasta,
            &CallRegionArgs::strict_java(),
        )
        .expect("call")
        .expect("outcome")
    };
    let reuse =
        |a: &AssemblyRegion, b: &AssemblyRegion, oa: &gatk_haplotypecaller::CallRegionOutcome| {
            if std::ptr::eq(a, b) {
                oa.clone()
            } else {
                call(b)
            }
        };
    let gap_outcome = call(covering_gt);
    let ag_outcome = reuse(covering_gt, covering_ag, &gap_outcome);
    let ac_outcome = reuse(covering_gt, covering_ac, &gap_outcome);
    let indel_outcome = call(covering_indel);
    let tg_snp_outcome = reuse(covering_indel, covering_tg, &indel_outcome);
    let qual_outcome = reuse(covering_indel, covering_qual, &indel_outcome);
    let outcome = reuse(covering_indel, covering, &indel_outcome);

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
    let gt_snp = find(&gap_outcome, CLOSED_GT, "G", "T");
    let ag = find(&ag_outcome, CLOSED_AG, "A", "G");
    let ac = find(&ac_outcome, CLOSED_AC, "A", "C");
    let indel = find(&indel_outcome, CLOSED_INDEL, "TTC", "T");
    let tg = find(&tg_snp_outcome, CLOSED_TG, "T", "G");
    let qual = find(&qual_outcome, CLOSED_QUAL, "CT", "C");
    let target = find(&outcome, TARGET, TARGET_REF, TARGET_ALT);

    let sem = include_str!("../src/compatibility/java_hc_site_semantics.rs");
    assert!(
        sem.contains("(92307403, \"C\", \"A\")"),
        "target is CLUSTER_DOWNSTREAM_SNPS"
    );
    assert!(
        sem.contains("CLUSTER_TG_SNP_START: u64 = 92307333"),
        "control T/G remains cluster-TG"
    );

    let stored = unique_indices(&outcome.read_likelihoods);
    let attached = unique_indices(&target.annotation_likelihoods);
    kv("stored_hap_n", stored.len().to_string());
    kv("attached_n", attached.len().to_string());
    assert_eq!(stored.len(), 6, "stored haplotype unique evidence is n=6");
    assert!(
        !target.annotation_likelihoods.is_empty(),
        "annotation object is non-empty"
    );
    assert_eq!(
        attached.len(),
        6,
        "attached annotation object is stored unique n=6"
    );
    assert_eq!(
        attached, stored,
        "attached membership is the stored unique hap evidence"
    );

    let retain_n = stored
        .iter()
        .filter(|idx| {
            java_alignment_read_overlaps_interval(
                &outcome.genotyping_reads[**idx],
                TARGET,
                TARGET,
                MARGIN,
            )
        })
        .count();
    kv("retainEvidence_n", retain_n.to_string());
    assert_eq!(
        retain_n, 4,
        "6R.186 retain recipe remains n=4; not used here"
    );
    assert!(
        has_member(
            &outcome.genotyping_reads,
            &attached,
            EXTRA_QNAME_A,
            EXTRA_FLAG_A
        ) && has_member(
            &outcome.genotyping_reads,
            &attached,
            EXTRA_QNAME_B,
            EXTRA_FLAG_B
        ),
        "the two extra stored rows remain in the annotation object"
    );

    let java_coverage_n = attached.len();
    let rust_overlap_dp = coverage_evidence_count(
        &outcome.genotyping_reads,
        &target.annotation_likelihoods,
        TARGET,
        TARGET,
        MARGIN,
    );
    kv("coverage_evidenceCount_n", java_coverage_n.to_string());
    kv("coverage_overlap_filter_n", rust_overlap_dp.to_string());
    assert_eq!(
        java_coverage_n, 6,
        "Coverage sees six evidence rows on the attached object"
    );

    let mqs =
        rms_mapping_quality_sample_mapqs(&outcome.genotyping_reads, &target.annotation_likelihoods);
    let rms = rms_mapping_quality_raw(&mqs).map(|t| t.2).expect("rms");
    kv("mq_n", mqs.len().to_string());
    kv("mq_rms", format!("{rms:.2}"));
    assert_eq!(mqs.len(), 6, "MQ sees the six-read stored population");
    assert!(
        (rms - 40.58).abs() < 0.01,
        "MQ RMS of the attached object is Java 40.58, not retain n=4 ~44.15"
    );

    assert_eq!(unique_indices(&gt_snp.annotation_likelihoods).len(), 3);
    assert_eq!(unique_indices(&ag.annotation_likelihoods).len(), 3);
    assert_eq!(unique_indices(&ac.annotation_likelihoods).len(), 4);
    assert_eq!(unique_indices(&indel.annotation_likelihoods).len(), 1);
    assert_eq!(unique_indices(&tg.annotation_likelihoods).len(), 1);
    assert_eq!(unique_indices(&qual.annotation_likelihoods).len(), 0);

    assert_eq!(target.genotype.format.ad_as_i32(), vec![2, 4]);
    assert_eq!(target.genotype.format.dp.as_i32(), 6);
    assert_eq!(target.genotype.format.gq.as_i32(), 72);
    assert_eq!(target.genotype.format.pl_as_i32(), vec![162, 0, 72]);

    let emitted =
        try_emit_call_region_variants(covering, &outcome, "NA12878", DEFAULT_STAND_EMIT_CONFIDENCE)
            .expect("emit");
    let rec = emitted
        .iter()
        .find(|r| r.position == TARGET && r.reference == TARGET_REF)
        .expect("C/A emit");
    let qual_v = rec.quality.expect("QUAL");
    kv("qual", format!("{qual_v:.2}"));
    kv("info", format_info(&rec.info));
    assert!((qual_v - 154.64).abs() < 0.005, "QUAL remains 154.64");
    let sample = rec.samples.first().expect("sample");
    assert_eq!(
        sample.gt.as_ref().map(|g| g.to_string()).as_deref(),
        Some("0/1")
    );
    assert_eq!(sample.ad.as_deref(), Some(&[2u32, 4][..]));
    assert_eq!(sample.dp, Some(6));
    assert_eq!(sample.gq, Some(72.0));
    assert_eq!(sample.pl.as_deref(), Some(&[162u32, 0, 72][..]));

    kv("info_dp", format!("{:?}", info_i32(&rec.info, "DP")));
    kv("info_mq", format!("{:?}", info_f64(&rec.info, "MQ")));
    kv(
        "info_baseq",
        format!("{:?}", info_f64(&rec.info, "BaseQRankSum")),
    );
    kv(
        "info_mqrs",
        format!("{:?}", info_f64(&rec.info, "MQRankSum")),
    );
    kv(
        "info_readpos",
        format!("{:?}", info_f64(&rec.info, "ReadPosRankSum")),
    );
    kv(
        "info_has_baseq",
        info_has(&rec.info, "BaseQRankSum").to_string(),
    );
    kv(
        "info_has_readpos",
        info_has(&rec.info, "ReadPosRankSum").to_string(),
    );

    let indel_emitted = try_emit_call_region_variants(
        covering_indel,
        &indel_outcome,
        "NA12878",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .expect("indel emit");
    let tg_emitted = try_emit_call_region_variants(
        covering_tg,
        &tg_snp_outcome,
        "NA12878",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .expect("tg emit");
    let qual_emitted = try_emit_call_region_variants(
        covering_qual,
        &qual_outcome,
        "NA12878",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .expect("qual emit");
    let indel_rec = indel_emitted
        .iter()
        .find(|r| r.position == CLOSED_INDEL && r.reference == "TTC")
        .expect("TTC/T");
    let tg_rec = tg_emitted
        .iter()
        .find(|r| r.position == CLOSED_TG && r.reference == "T")
        .expect("T/G");
    let pin = qual_emitted
        .iter()
        .find(|r| r.position == CLOSED_QUAL && r.reference == "CT")
        .expect("CT/C");
    assert_eq!(info_i32(&indel_rec.info, "DP"), Some(1));
    assert_eq!(info_i32(&tg_rec.info, "DP"), Some(1));
    assert_eq!(info_i32(&pin.info, "DP"), Some(2));
    assert!((pin.quality.expect("pin QUAL") - 31.60).abs() < 0.005);
    assert!((info_f64(&pin.info, "QD").expect("QD") - 15.80).abs() < 0.005);
    assert_eq!(
        (info_f64(&pin.info, "BaseQRankSum").expect("BaseQ") * 1000.0).round(),
        0.0
    );
    assert_eq!(
        (info_f64(&pin.info, "MQRankSum").expect("MQRS") * 1000.0).round(),
        -674.0
    );
    assert!(info_has(&pin.info, "ReadPosRankSum"));
}
