//! 6R.210: attach Java loc-loop annotation on the gap-sparse shaped-early
//! FORMAT path at `2:92318199 C/T`.
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r210_annotation_source_boundary -- --nocapture --test-threads=1
//! HOLDOUT_6R210=1 cargo test -p gatk-haplotypecaller --test holdout_6r210_annotation_source_boundary -- --nocapture --test-threads=1
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
const JAVA_MQRS: &str = r#"
if (alt.isEmpty() || reference.isEmpty()) {
    return Collections.emptyMap();
}
"#;
const GAP_INTERVAL: &str = "2:92305500-92305850";
const TG_INTERVAL: &str = "2:92307200-92307550";
const MID_INTERVAL: &str = "2:92316200-92316580";
const MID_B_INTERVAL: &str = "2:92317000-92319000";
const P12_POST_INTERVAL: &str = "2:92318150-92319220";
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 92_318_199;
const TARGET_REF: &str = "C";
const TARGET_ALT: &str = "T";
const JAVA_QNAME: &str = "H06JUADXX130110:1:1101:10002:20630";
const JAVA_FLAG: u16 = 99;
const JAVA_MAPQ: u8 = 24;
const CLOSED_AG: u64 = 92_305_635;
const CLOSED_TG: u64 = 92_307_333;
const CLOSED_QUAL: u64 = 92_307_359;
const CLOSED_CA: u64 = 92_307_403;
const CLOSED_HOM_ALT: u64 = 92_316_296;
const CLOSED_ONE_READ: u64 = 92_316_416;
const CLOSED_MID_B: u64 = 92_317_399;
const MARGIN: i32 = 2;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R210\t{key}\t{}", value.as_ref());
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

fn info_has(info: &[InfoValue], key: &str) -> bool {
    info.iter().any(|v| match v {
        InfoValue::Flag(k)
        | InfoValue::Integer(k, _)
        | InfoValue::Float(k, _)
        | InfoValue::String(k, _)
        | InfoValue::Character(k, _) => k == key,
    })
}

fn java_calculate_sor(ref_fw: u32, ref_rv: u32, alt_fw: u32, alt_rv: u32) -> f64 {
    let t00 = f64::from(ref_fw) + 1.0;
    let t01 = f64::from(ref_rv) + 1.0;
    let t10 = f64::from(alt_fw) + 1.0;
    let t11 = f64::from(alt_rv) + 1.0;
    let ratio = (t00 / t01) * (t11 / t10) + (t01 / t00) * (t10 / t11);
    let ref_ratio = t00.min(t01) / t00.max(t01);
    let alt_ratio = t10.min(t11) / t10.max(t11);
    ratio.ln() + ref_ratio.ln() - alt_ratio.ln()
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

fn gap_sparse_arm(early: &str) -> &str {
    early
        .split("if gap_sparse_shaped_early")
        .nth(1)
        .expect("gap_sparse_shaped_early")
        .split("let (_trim_pileup_ref, trim_pileup_alt)")
        .next()
        .expect("gap body")
}

#[test]
fn forensic_6r210_annotation_source_boundary_contract() {
    kv("java_pin", JAVA_PIN);
    kv(
        "production_change",
        "one annotation-object construction arrow on gap_sparse_shaped_early",
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
        "default HC annotation reuses genotyping likelihoods (contamination off) then addEvidence filtered overlap"
    );
    assert!(
        JAVA_MQRS.contains("alt.isEmpty()") && JAVA_MQRS.contains("reference.isEmpty()"),
        "Java MappingQualityRankSumTest omits when REF or ALT fillQuals list is empty"
    );

    let early = include_str!("../src/hc_genotyping_engine/genotype_site_early_template.rs");
    let emit = include_str!("../src/region_vcf_emit.rs");
    let ann = include_str!("../src/variant_site_hc_annotations.rs");
    let pipe = include_str!("../src/hc_genotyping_engine/genotype_site_pipeline.rs");
    let semantics = include_str!("../src/compatibility/java_hc_site_semantics.rs");
    let discovery = include_str!("../src/read_event_discovery/mod.rs");

    let gap_arm = gap_sparse_arm(early);
    assert!(
        gap_arm.contains("apply_sparse_shaped_hom_alt_rescue(0, fmt_alt, config)")
            && gap_arm.contains("Some((0, fmt_alt))"),
        "gap-sparse shaped-early still uses the shaped FORMAT path"
    );
    assert!(
        gap_arm.contains("annotation_likelihoods_from_stored_haplotypes")
            && gap_arm.contains("with_annotation_likelihoods"),
        "6R.210 attaches the stored-haplotype loc-loop object on this lifecycle"
    );
    assert!(
        !gap_arm.contains("annotation_likelihoods_from_stored_unique_evidence"),
        "6R.199 stored-unique must not be used on this branch"
    );
    assert!(
        !gap_arm.contains("92318199")
            && !gap_arm.contains(JAVA_QNAME)
            && !gap_arm.contains("FLAG=99")
            && !gap_arm.contains("10088:6763"),
        "production arm must not QNAME/FLAG/locus-filter rows"
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
            && !retain_fn.contains("92318199"),
        "helper stays general Java marginalize + retainEvidence(±2); 6R.209 skip of QNAME collapse stays"
    );

    let one_read_arm = early
        .split("if is_mid_a_one_read_hom_alt_site(&event)")
        .nth(1)
        .expect("one-read arm")
        .split("if is_p12_phase_e_two_read_hom_alt_site(&event)")
        .next()
        .expect("arm body");
    assert!(
        one_read_arm.contains("annotation_likelihoods_from_stored_haplotypes"),
        "6R.207 one-read attach stays; this site is not that arm"
    );
    let two_read_arm = early
        .split("if is_p12_phase_e_two_read_hom_alt_site(&event)")
        .nth(1)
        .expect("two-read arm")
        .split("let outside_trim")
        .next()
        .expect("two-read body");
    assert!(
        two_read_arm.contains("annotation_likelihoods_from_stored_haplotypes"),
        "6R.205 two-read attach stays; this site is not that arm"
    );
    let ds_body = early
        .split("if is_cluster_downstream_snp(&event)")
        .nth(1)
        .expect("downstream arm")
        .split("if is_cluster_tg_snp(&event)")
        .next()
        .expect("downstream body");
    assert!(
        ds_body.contains("annotation_likelihoods_from_stored_unique_evidence")
            && !ds_body.contains("annotation_likelihoods_from_stored_haplotypes"),
        "6R.199 cluster-downstream stays stored-unique"
    );
    let tg_body = early
        .split("if is_cluster_tg_snp(&event)")
        .nth(1)
        .expect("tg arm")
        .split("if is_cluster_tc_snp(&event)")
        .next()
        .expect("tg body");
    assert!(
        tg_body.contains("annotation_likelihoods_from_stored_haplotypes"),
        "6R.186 cluster-TG attach stays"
    );
    assert!(
        pipe.contains("let annotation = annotation_likelihoods_from_stored_haplotypes")
            && pipe.contains("with_annotation_likelihoods(annotation)"),
        "6R.189 SiteScore attach stays"
    );
    assert!(
        !semantics
            .split("pub fn is_mid_a_one_read_hom_alt_site")
            .nth(1)
            .expect("one-read fn")
            .split("pub fn is_mid_a_two_read_hom_alt_site")
            .next()
            .expect("body")
            .contains("92318199"),
        "6R.207 one-read predicate does not include this locus"
    );
    assert!(
        !semantics
            .split("pub fn is_java_sparse_two_read_hom_alt_site")
            .nth(1)
            .expect("two-read fn")
            .split("#[cfg(test)]")
            .next()
            .expect("body")
            .contains("92318199"),
        "6R.205 two-read predicate does not include this locus"
    );
    assert!(
        discovery.contains("(92318199, \"C\", \"T\")"),
        "target is a P12 gap-registry SNP, not mid-A / two-read / cluster-TG"
    );
    assert!(
        !early.contains("92318199") && !ann.contains("92318199") && !emit.contains("92318199"),
        "no locus-specific 6R.210 production patch"
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
        "6R.201 Coverage cardinality stays: no second ±2 overlap filter"
    );
    assert!(
        ann.contains("fn mq_rms_of_sample_evidence")
            && ann.contains("strand_bias_sample_counts")
            && !ann.contains("0.693")
            && !ann.contains("DP = 1"),
        "MQ/SOR/Coverage formulas are not patched in this round"
    );
    assert!(
        emit.contains("if call.annotation_likelihoods.is_empty()")
            && emit.contains("likelihoods: &call.annotation_likelihoods"),
        "emit still consumes the attached object when present"
    );
    let ge = include_str!("../src/hc_genotyping_engine/mod.rs");
    assert!(
        ge.contains("if qnames.len() == 1 && best_ll.len() > 1"),
        "6R.180 same-QNAME collapse remains on FORMAT-subset / stored-unique callers"
    );
}

#[test]
fn forensic_6r210_annotation_source_boundary() {
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
        "one annotation-object construction arrow on gap_sparse_shaped_early",
    );
    kv("target", "2:92318199 C/T");
    kv("classification", "ANNOTATION_SOURCE_OBJECT_DIVERGENCE");

    let java_vcf = root.join("parity/reports/6r43/p12_post/java.vcf");
    if java_vcf.is_file() {
        let text = std::fs::read_to_string(&java_vcf).expect("java.vcf");
        let line = text
            .lines()
            .find(|l| l.contains("\t92318199\t") && l.contains("\tC\tT\t"))
            .expect("frozen Java 4.4.0.0 p12_post 2:92318199 C/T");
        let fields: Vec<&str> = line.split('\t').collect();
        kv("java_vcf_qual", fields[5]);
        kv("java_vcf_info", fields[7]);
        kv("java_vcf_format", format!("{} {}", fields[8], fields[9]));
        assert_eq!(fields[5], "35.48");
        assert!(
            fields[7].contains("DP=1")
                && fields[7].contains("MQ=24.00")
                && fields[7].contains("SOR=1.609")
                && !fields[7].contains("MQRankSum")
        );
        assert!(
            fields[9].contains("1/1") && fields[9].contains("0,1") && fields[9].contains("45,3,0")
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
    let gap_regions = flatten_assembly_regions(&walk(GAP_INTERVAL));
    let tg_regions = flatten_assembly_regions(&walk(TG_INTERVAL));
    let mid_regions = flatten_assembly_regions(&walk(MID_INTERVAL));
    let mid_b_regions = flatten_assembly_regions(&walk(MID_B_INTERVAL));
    let post_regions = flatten_assembly_regions(&walk(P12_POST_INTERVAL));
    let covering_ag = gap_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_AG
                && r.end.get() >= CLOSED_AG
        })
        .expect("gap covering A/G");
    let covering_tg = tg_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_TG
                && r.end.get() >= CLOSED_TG
        })
        .expect("tg covering T/G");
    let covering_qual = tg_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_QUAL
                && r.end.get() >= CLOSED_QUAL
        })
        .expect("tg covering CT/C");
    let covering_ca = tg_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_CA
                && r.end.get() >= CLOSED_CA
        })
        .expect("tg covering C/A");
    let covering_hom = mid_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_HOM_ALT
                && r.end.get() >= CLOSED_HOM_ALT
        })
        .expect("mid covering A/T");
    let covering_one = mid_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_ONE_READ
                && r.end.get() >= CLOSED_ONE_READ
        })
        .expect("mid covering C/A");
    let covering_mid_b = mid_b_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_MID_B
                && r.end.get() >= CLOSED_MID_B
        })
        .expect("mid-B covering C/A");
    let covering_target = post_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= TARGET
                && r.end.get() >= TARGET
        })
        .expect("p12_post covering 92318199");

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
    let ag_outcome = call_one(covering_ag);
    let tg_outcome = call_one(covering_tg);
    let qual_outcome = if std::ptr::eq(covering_tg, covering_qual) {
        tg_outcome.clone()
    } else {
        call_one(covering_qual)
    };
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
    let ag = find(&ag_outcome, CLOSED_AG, "A", "G");
    let tg = find(&tg_outcome, CLOSED_TG, "T", "G");
    let qual = find(&qual_outcome, CLOSED_QUAL, "CT", "C");
    let ca = find(&ca_outcome, CLOSED_CA, "C", "A");
    let hom = find(&hom_outcome, CLOSED_HOM_ALT, "A", "T");
    let one = find(&one_outcome, CLOSED_ONE_READ, "C", "A");
    let mid_b = find(&mid_b_outcome, CLOSED_MID_B, "C", "A");
    let target = find(&outcome, TARGET, TARGET_REF, TARGET_ALT);

    kv(
        "branch",
        "gap_sparse_shaped_early attaches annotation_likelihoods_from_stored_haplotypes",
    );
    kv(
        "active_region",
        format!(
            "{}-{}",
            covering_target.start.get(),
            covering_target.end.get()
        ),
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

    assert_eq!(target.genotype.format.ad_as_i32(), vec![0, 1]);
    assert_eq!(target.genotype.format.dp.as_i32(), 1);
    assert_eq!(target.genotype.format.gq.as_i32(), 3);
    assert_eq!(target.genotype.format.pl_as_i32(), vec![45, 3, 0]);
    kv("format", "GT=1/1 AD=0,1 DP=1 GQ=3 PL=45,3,0");

    assert_eq!(unique_indices(&ag.annotation_likelihoods).len(), 3);
    assert_eq!(unique_indices(&tg.annotation_likelihoods).len(), 1);
    assert_eq!(unique_indices(&qual.annotation_likelihoods).len(), 0);
    assert_eq!(unique_indices(&ca.annotation_likelihoods).len(), 6);
    assert_eq!(unique_indices(&hom.annotation_likelihoods).len(), 2);
    assert_eq!(unique_indices(&one.annotation_likelihoods).len(), 1);
    assert_eq!(unique_indices(&mid_b.annotation_likelihoods).len(), 2);

    assert_eq!(stored.len(), 4, "p12_post stored unique n=4");
    assert_eq!(retain.len(), 1, "Java loc-loop retainEvidence(±2) n=1");
    assert_eq!(
        attached.len(),
        1,
        "6R.210: loc-loop annotation is retainEvidence n=1"
    );
    assert_eq!(attached, retain);
    assert!(
        !target.annotation_likelihoods.is_empty(),
        "region-wide fallback is not taken; annotation is attached"
    );

    let rec_meta = |idx: usize| {
        let rec: &Record = outcome.genotyping_reads[idx].as_ref();
        (qname(rec), rec.flags(), rec.mapq(), rec.is_reverse())
    };
    let attached_meta: Vec<_> = attached.iter().copied().map(rec_meta).collect();
    assert_eq!(attached_meta.len(), 1);
    assert_eq!(attached_meta[0].0, JAVA_QNAME);
    assert_eq!(attached_meta[0].1, JAVA_FLAG);
    assert_eq!(attached_meta[0].2, JAVA_MAPQ);
    assert!(!attached_meta[0].3);
    for idx in &stored {
        let meta = rec_meta(*idx);
        if meta.0 != JAVA_QNAME || meta.1 != JAVA_FLAG {
            assert!(
                !attached.contains(idx),
                "non-overlapping stored unique row must not survive loc-loop retainEvidence"
            );
        }
    }

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
    assert_eq!(attached_dp, 1, "Coverage consumes attached n=1");
    assert_eq!(stored_dp, 4);
    assert_eq!(
        attached_dp,
        attached.len() as i32,
        "Coverage does not apply a second ±2 filter"
    );

    let attached_mqs =
        rms_mapping_quality_sample_mapqs(&outcome.genotyping_reads, &target.annotation_likelihoods);
    kv("attached_mapq_list", format!("{attached_mqs:?}"));
    assert_eq!(attached_mqs, vec![JAVA_MAPQ]);
    let mq_rms = rms_mapping_quality_raw(&attached_mqs).map(|t| (t.2 * 100.0).round() / 100.0);
    assert_eq!(mq_rms, Some(24.0));

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
    assert_eq!(attached_table, (0, 0, 1, 0));
    let java_sor = java_calculate_sor(
        attached_table.0,
        attached_table.1,
        attached_table.2,
        attached_table.3,
    );
    assert!((java_sor - 1.609).abs() < 0.002);

    let (ref_mq, alt_mq) =
        mq_rank_sum_quals_from_likelihoods(&attached_evidence, TARGET, TARGET_REF, TARGET_ALT);
    kv("attached_mqrs_ref", format!("{ref_mq:?}"));
    kv("attached_mqrs_alt", format!("{alt_mq:?}"));
    assert!(
        ref_mq.is_empty() && alt_mq == [24.0],
        "Java MQRankSum omits: empty REF fillQuals on the loc-loop object"
    );

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
    let target_emitted = emit(covering_target, &outcome);

    let tg_rec = rec_at(&tg_emitted, CLOSED_TG, "T");
    assert_eq!(info_i32(&tg_rec.info, "DP"), Some(1));
    assert!((info_f64(&tg_rec.info, "MQ").expect("mq") - 44.00).abs() < 0.005);
    assert!((info_f64(&tg_rec.info, "SOR").expect("sor") - 1.609).abs() < 0.002);
    let ca_rec = rec_at(&ca_emitted, CLOSED_CA, "C");
    assert_eq!(info_i32(&ca_rec.info, "DP"), Some(6));
    assert!((info_f64(&ca_rec.info, "MQ").expect("mq") - 40.58).abs() < 0.005);
    assert!((info_f64(&ca_rec.info, "SOR").expect("sor") - 0.105).abs() < 0.002);
    let hom_rec = rec_at(&hom_emitted, CLOSED_HOM_ALT, "A");
    assert_eq!(info_i32(&hom_rec.info, "DP"), Some(2));
    assert!((info_f64(&hom_rec.info, "MQ").expect("mq") - 47.00).abs() < 0.005);
    assert!((info_f64(&hom_rec.info, "SOR").expect("sor") - 2.303).abs() < 0.002);
    let one_rec = rec_at(&one_emitted, CLOSED_ONE_READ, "C");
    assert_eq!(info_i32(&one_rec.info, "DP"), Some(1));
    assert!((info_f64(&one_rec.info, "MQ").expect("mq") - 21.00).abs() < 0.005);
    assert!((info_f64(&one_rec.info, "SOR").expect("sor") - 1.609).abs() < 0.002);
    let mid_b_rec = rec_at(&mid_b_emitted, CLOSED_MID_B, "C");
    assert_eq!(info_i32(&mid_b_rec.info, "DP"), Some(2));
    assert!((info_f64(&mid_b_rec.info, "MQ").expect("mq") - 27.00).abs() < 0.005);
    assert!((info_f64(&mid_b_rec.info, "SOR").expect("sor") - 0.693).abs() < 0.002);

    let rec = rec_at(&target_emitted, TARGET, TARGET_REF);
    assert!((rec.quality.expect("QUAL") - 35.48).abs() < 0.005);
    let dp = info_i32(&rec.info, "DP").expect("DP");
    let mq = info_f64(&rec.info, "MQ").expect("MQ");
    let sor = info_f64(&rec.info, "SOR").expect("SOR");
    kv("rust_info_dp", dp.to_string());
    kv("java_info_dp", "1");
    kv("rust_mq", format!("{mq:.2}"));
    kv("java_mq", "24.00");
    kv("rust_sor", format!("{sor:.3}"));
    kv("java_sor", "1.609");
    let sample = rec.samples.first().expect("sample");
    assert_eq!(
        sample.gt.as_ref().map(|g| g.to_string()).as_deref(),
        Some("1/1")
    );
    assert_eq!(sample.ad.as_deref(), Some(&[0u32, 1][..]));
    assert_eq!(sample.dp, Some(1));
    assert_eq!(sample.gq, Some(3.0));
    assert_eq!(sample.pl.as_deref(), Some(&[45u32, 3, 0][..]));
    assert_eq!(dp, 1, "INFO DP from attached n=1");
    assert!((mq - 24.00).abs() < 0.005, "MQ RMS of [24] is 24.00");
    assert!((sor - 1.609).abs() < 0.002, "INFO SOR from [0,0;1,0]");
    assert!(
        !info_has(&rec.info, "MQRankSum"),
        "Java omits MQRankSum: empty REF fillQuals on loc-loop n=1"
    );

    kv(
        "first_arrow",
        "gap_sparse_shaped_early FORMAT-only return left annotation empty; emit used region-wide n=4",
    );
    kv(
        "java_survivor",
        format!("{JAVA_QNAME} FLAG={JAVA_FLAG} FWD MAPQ={JAVA_MAPQ}"),
    );
}
