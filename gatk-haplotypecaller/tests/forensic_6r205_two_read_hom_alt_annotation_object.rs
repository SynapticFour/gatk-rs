//! 6R.205: two-read hom-alt early-template attaches the Java-equivalent
//! loc-loop annotation object (`annotation_likelihoods_from_stored_haplotypes`).
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! Target `2:92316296 A/T`. One object-construction arrow. Coverage is not
//! retuned (6R.201). 6R.199 stored-unique is not used on this branch.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r205_two_read_hom_alt_annotation_object -- --nocapture --test-threads=1
//! HOLDOUT_6R205=1 cargo test -p gatk-haplotypecaller --test holdout_6r205_two_read_hom_alt_annotation -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::hc_genotyping_engine::java_alignment_read_overlaps_interval;
use gatk_haplotypecaller::variant_site_hc_annotations::coverage_evidence_count;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    HaplotypeCallerEngine, ReadFilterParams, RegionReadLikelihood, SharedBamRecord,
    WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
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
const GAP_INTERVAL: &str = "2:92305500-92305850";
const TG_INTERVAL: &str = "2:92307200-92307550";
const MID_INTERVAL: &str = "2:92316200-92316580";
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 92_316_296;
const TARGET_REF: &str = "A";
const TARGET_ALT: &str = "T";
const CLOSED_GT: u64 = 92_305_634;
const CLOSED_AG: u64 = 92_305_635;
const CLOSED_AC: u64 = 92_305_716;
const CLOSED_INDEL: u64 = 92_307_324;
const CLOSED_TG: u64 = 92_307_333;
const CLOSED_QUAL: u64 = 92_307_359;
const CLOSED_CA: u64 = 92_307_403;
const CLOSED_MID: u64 = 92_316_347;
const MARGIN: i32 = 2;
const JAVA_QNAME_A: &str = "H06HDADXX130110:1:1101:10034:45116";
const JAVA_FLAG_A: u16 = 99;
const JAVA_QNAME_B: &str = "H06HDADXX130110:2:1101:10025:49248";
const JAVA_FLAG_B: u16 = 99;
const EXCLUDED_QNAME: &str = "H06HDADXX130110:2:1101:10046:78083";
const EXCLUDED_FLAG: u16 = 147;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R205\t{key}\t{}", value.as_ref());
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

fn alignment_end_1based(rec: &Record) -> i64 {
    let start = rec.pos() + 1;
    let mut ref_len = 0i64;
    for c in rec.cigar().iter() {
        match c {
            rust_htslib::bam::record::Cigar::Match(n)
            | rust_htslib::bam::record::Cigar::Equal(n)
            | rust_htslib::bam::record::Cigar::Diff(n)
            | rust_htslib::bam::record::Cigar::Del(n)
            | rust_htslib::bam::record::Cigar::RefSkip(n) => ref_len += i64::from(*n),
            _ => {}
        }
    }
    if ref_len > 0 {
        start + ref_len - 1
    } else {
        start
    }
}

fn dump_reads(label: &str, reads: &[SharedBamRecord], idxs: &BTreeSet<usize>) {
    for idx in idxs {
        let rec: &Record = reads[*idx].as_ref();
        let overlap = java_alignment_read_overlaps_interval(rec, TARGET, TARGET, MARGIN);
        kv(
            label,
            format!(
                "idx={idx} qname={} flag={} mapq={} contig={} start={} end={} cigar={} overlap_pm2={}",
                qname(rec),
                rec.flags(),
                rec.mapq(),
                rec.tid(),
                rec.pos() + 1,
                alignment_end_1based(rec),
                rec.cigar(),
                overlap
            ),
        );
    }
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

fn two_read_arm_body(early: &str) -> &str {
    early
        .split("if is_p12_phase_e_two_read_hom_alt_site(&event)")
        .nth(1)
        .expect("two-read arm")
        .split("let outside_trim")
        .next()
        .expect("arm body")
}

#[test]
fn forensic_6r205_two_read_hom_alt_annotation_object_contract() {
    kv("java_pin", JAVA_PIN);
    kv(
        "production_change",
        "one annotation-object construction arrow on is_p12_phase_e_two_read_hom_alt_site",
    );
    kv("classification", "ANNOTATION_SOURCE_OBJECT_DIVERGENCE");
    assert!(
        JAVA_COVERAGE_ANNOTATE.contains("likelihoods.evidenceCount()")
            && !JAVA_COVERAGE_ANNOTATE.contains("overlaps"),
        "Java Coverage.annotate is evidenceCount of the passed AlleleLikelihoods"
    );

    let early = include_str!("../src/hc_genotyping_engine/genotype_site_early_template.rs");
    let emit = include_str!("../src/region_vcf_emit.rs");
    let ann = include_str!("../src/variant_site_hc_annotations.rs");
    let pipe = include_str!("../src/hc_genotyping_engine/genotype_site_pipeline.rs");
    let fin = include_str!("../src/hc_genotyping_engine/genotype_finalize.rs");

    let two_read_arm = two_read_arm_body(early);
    assert!(
        two_read_arm.contains("finish_strict_java_shaped_site_call"),
        "two-read hom-alt still uses the shaped FORMAT path"
    );
    assert!(
        two_read_arm.contains("annotation_likelihoods_from_stored_haplotypes")
            && two_read_arm.contains("with_annotation_likelihoods"),
        "6R.205 attaches the stored-haplotype loc-loop object"
    );
    assert!(
        !two_read_arm.contains("annotation_likelihoods_from_stored_unique_evidence"),
        "6R.199 stored-unique helper must not be used on this branch"
    );
    assert!(
        !two_read_arm.contains(JAVA_QNAME_A)
            && !two_read_arm.contains(JAVA_QNAME_B)
            && !two_read_arm.contains(EXCLUDED_QNAME)
            && !two_read_arm.contains("FLAG=147")
            && !two_read_arm.contains("92316296"),
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
            && !retain_fn.contains("92316296")
            && !retain_fn.contains(EXCLUDED_QNAME),
        "helper stays general Java marginalize + retainEvidence(±2)"
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
        "6R.199 cluster-downstream still uses stored unique, not retain"
    );

    let tg_body = early
        .split("if is_cluster_tg_snp(&event)")
        .nth(1)
        .expect("tg arm")
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
        cov.contains("seen.insert(cell.read_index.get())")
            && !cov.contains("92316296")
            && !cov.contains("DP = 2")
            && !cov.contains("info_dp"),
        "Coverage is unique evidenceCount; INFO DP is not patched"
    );
    assert!(
        ann.contains("fn mq_rms_of_sample_evidence")
            && ann.contains("strand_bias_sample_counts")
            && !ann.contains("47.00")
            && !ann.contains("2.303"),
        "MQ/SOR formulas are not directly patched"
    );
    assert!(
        emit.contains("if call.annotation_likelihoods.is_empty()")
            && emit.contains("likelihoods: &call.annotation_likelihoods"),
        "emit still consumes the attached object when present"
    );

    let fin_fn = fin
        .split("fn finish_strict_java_shaped_site_call")
        .nth(1)
        .expect("finish fn");
    let fin_body = fin_fn.split("\nfn ").next().expect("finish body");
    assert!(
        !fin_body.contains("with_annotation_likelihoods"),
        "shared finish helper does not globally attach annotation_likelihoods"
    );
    assert!(
        !ann.contains("92316296") && !early.contains("92316296") && !emit.contains("92316296"),
        "no locus-specific 6R.205 production patch"
    );
}

#[test]
fn forensic_6r205_two_read_hom_alt_annotation_object() {
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
        "one annotation-object construction arrow",
    );
    kv("target", "2:92316296 A/T");
    kv(
        "branch",
        "is_p12_phase_e_two_read_hom_alt_site / is_java_sparse_two_read_hom_alt_site",
    );

    let java_vcf = root.join("parity/reports/6r43/p12_mid_a/java.vcf");
    if java_vcf.is_file() {
        let text = std::fs::read_to_string(&java_vcf).expect("java.vcf");
        let line = text
            .lines()
            .find(|l| l.contains("\t92316296\t") && l.contains("\tA\tT\t"))
            .expect("frozen Java 4.4.0.0 p12_mid_a 2:92316296 A/T");
        let fields: Vec<&str> = line.split('\t').collect();
        kv("java_vcf_qual", fields[5]);
        kv("java_vcf_info", fields[7]);
        kv("java_vcf_format", format!("{} {}", fields[8], fields[9]));
        assert!(
            fields[7].contains("DP=2"),
            "pinned Java INFO DP=2, got {}",
            fields[7]
        );
        assert!(
            fields[7].contains("MQ=47.00") && fields[7].contains("SOR=2.303"),
            "pinned Java MQ/SOR recorded, got {}",
            fields[7]
        );
        assert_eq!(fields[5], "78.32");
        assert!(
            fields[9].contains("1/1") && fields[9].contains("0,2") && fields[9].contains("90,6,0")
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
    let covering_tg = |pos: u64| {
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
    let covering_mid = |pos: u64| {
        mid_regions
            .iter()
            .find(|r| {
                matches!(
                    call_disposition(r),
                    AssemblyRegionCallDisposition::ActiveFull
                ) && r.start.get() <= pos
                    && r.end.get() >= pos
            })
            .unwrap_or_else(|| panic!("mid covering {pos}"))
    };

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
    let covering_gt = covering_gap(CLOSED_GT);
    let covering_ag = covering_gap(CLOSED_AG);
    let covering_ac = covering_gap(CLOSED_AC);
    let covering_indel = covering_tg(CLOSED_INDEL);
    let covering_tg_snp = covering_tg(CLOSED_TG);
    let covering_qual = covering_tg(CLOSED_QUAL);
    let covering_ca = covering_tg(CLOSED_CA);
    let covering_target = covering_mid(TARGET);
    let covering_closed_mid = covering_mid(CLOSED_MID);

    let gap_outcome = call_one(covering_gt);
    let ag_outcome = if std::ptr::eq(covering_gt, covering_ag) {
        gap_outcome.clone()
    } else {
        call_one(covering_ag)
    };
    let ac_outcome = if std::ptr::eq(covering_gt, covering_ac) {
        gap_outcome.clone()
    } else {
        call_one(covering_ac)
    };
    let indel_outcome = call_one(covering_indel);
    let tg_outcome = if std::ptr::eq(covering_indel, covering_tg_snp) {
        indel_outcome.clone()
    } else {
        call_one(covering_tg_snp)
    };
    let qual_outcome = if std::ptr::eq(covering_indel, covering_qual) {
        indel_outcome.clone()
    } else {
        call_one(covering_qual)
    };
    let ca_outcome = if std::ptr::eq(covering_indel, covering_ca) {
        indel_outcome.clone()
    } else {
        call_one(covering_ca)
    };
    let outcome = call_one(covering_target);
    let mid_outcome = if std::ptr::eq(covering_target, covering_closed_mid) {
        outcome.clone()
    } else {
        call_one(covering_closed_mid)
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
    let gt_snp = find(&gap_outcome, CLOSED_GT, "G", "T");
    let ag = find(&ag_outcome, CLOSED_AG, "A", "G");
    let ac = find(&ac_outcome, CLOSED_AC, "A", "C");
    let indel = find(&indel_outcome, CLOSED_INDEL, "TTC", "T");
    let tg = find(&tg_outcome, CLOSED_TG, "T", "G");
    let qual = find(&qual_outcome, CLOSED_QUAL, "CT", "C");
    let ca = find(&ca_outcome, CLOSED_CA, "C", "A");
    let mid = find(&mid_outcome, CLOSED_MID, "G", "A");
    let target = find(&outcome, TARGET, TARGET_REF, TARGET_ALT);

    let semantics = include_str!("../src/compatibility/java_hc_site_semantics.rs");
    assert!(
        semantics.contains("s == GenomePosition::new_1based(92316296) && r == \"A\" && a == \"T\""),
        "target is pinned in is_java_sparse_two_read_hom_alt_site"
    );

    let stored = unique_indices(&outcome.read_likelihoods);
    let attached = unique_indices(&target.annotation_likelihoods);
    kv("stored_hap_n", stored.len().to_string());
    kv("attached_annotation_n", attached.len().to_string());
    kv(
        "attached_cells",
        target.annotation_likelihoods.len().to_string(),
    );
    assert!(
        !target.annotation_likelihoods.is_empty(),
        "after 6R.205 the two-read hom-alt arm attaches annotation_likelihoods"
    );
    assert_eq!(stored.len(), 3, "region-wide stored unique remains n=3");
    assert_eq!(
        attached.len(),
        2,
        "attached Java-equivalent annotation object is n=2"
    );

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
    kv("java_loc_loop_retainEvidence_n", retain.len().to_string());
    dump_reads("stored_unique", &outcome.genotyping_reads, &stored);
    dump_reads("attached_annotation", &outcome.genotyping_reads, &attached);
    dump_reads("java_retainEvidence", &outcome.genotyping_reads, &retain);
    assert_eq!(
        retain.len(),
        2,
        "Java loc-loop retainEvidence(±2) on stored unique is n=2"
    );
    assert_eq!(
        attached, retain,
        "attached object membership equals retainEvidence(±2) on stored unique"
    );
    let extra: BTreeSet<usize> = stored.difference(&attached).copied().collect();
    assert_eq!(extra.len(), 1, "exactly one stored unique row is excluded");
    let extra_idx = *extra.iter().next().expect("extra idx");
    let extra_rec: &Record = outcome.genotyping_reads[extra_idx].as_ref();
    assert_eq!(qname(extra_rec), EXCLUDED_QNAME);
    assert_eq!(extra_rec.flags(), EXCLUDED_FLAG);
    assert_eq!(extra_rec.mapq(), 21);
    assert!(
        !java_alignment_read_overlaps_interval(extra_rec, TARGET, TARGET, MARGIN),
        "FLAG147 MAPQ21 row misses retainEvidence overlap at 92316296"
    );
    assert!(
        !attached.contains(&extra_idx),
        "FLAG147 row is not in the annotation object"
    );
    assert!(
        attached.iter().any(|&idx| {
            let rec: &Record = outcome.genotyping_reads[idx].as_ref();
            rec.flags() == JAVA_FLAG_A && qname(rec) == JAVA_QNAME_A
        }) && attached.iter().any(|&idx| {
            let rec: &Record = outcome.genotyping_reads[idx].as_ref();
            rec.flags() == JAVA_FLAG_B && qname(rec) == JAVA_QNAME_B
        }),
        "the two Java-counted MAPQ=47 overlapping rows remain"
    );

    let attached_dp = coverage_evidence_count(
        &outcome.genotyping_reads,
        &target.annotation_likelihoods,
        TARGET,
        TARGET,
        MARGIN,
    );
    let fallback_dp = coverage_evidence_count(
        &outcome.genotyping_reads,
        &outcome.read_likelihoods,
        TARGET,
        TARGET,
        MARGIN,
    );
    kv("coverage_attached_n", attached_dp.to_string());
    kv("coverage_stored_unique_n", fallback_dp.to_string());
    assert_eq!(attached_dp, 2, "Coverage receives the n=2 attached object");
    assert_eq!(
        fallback_dp, 3,
        "stored unique Coverage cardinality remains 3 (not used for this site)"
    );
    assert_eq!(
        attached_dp,
        attached.len() as i32,
        "Coverage output equals annotation input cardinality; no second overlap filter"
    );

    assert_eq!(target.genotype.format.ad_as_i32(), vec![0, 2]);
    assert_eq!(target.genotype.format.dp.as_i32(), 2);
    assert_eq!(target.genotype.format.gq.as_i32(), 6);
    assert_eq!(target.genotype.format.pl_as_i32(), vec![90, 6, 0]);
    kv("format", "GT=1/1 AD=0,2 DP=2 GQ=6 PL=90,6,0");

    assert_eq!(
        unique_indices(&gt_snp.annotation_likelihoods).len(),
        3,
        "neighbor G/T two-read arm now attaches retainEvidence n=3"
    );
    assert_eq!(unique_indices(&ag.annotation_likelihoods).len(), 3);
    assert_eq!(unique_indices(&ac.annotation_likelihoods).len(), 4);
    assert_eq!(unique_indices(&indel.annotation_likelihoods).len(), 1);
    assert_eq!(unique_indices(&tg.annotation_likelihoods).len(), 1);
    assert_eq!(unique_indices(&qual.annotation_likelihoods).len(), 0);
    assert_eq!(unique_indices(&ca.annotation_likelihoods).len(), 6);
    assert_eq!(unique_indices(&mid.annotation_likelihoods).len(), 3);

    let emit = |region: &gatk_haplotypecaller::AssemblyRegion,
                o: &gatk_haplotypecaller::CallRegionOutcome| {
        try_emit_call_region_variants(region, o, "NA12878", DEFAULT_STAND_EMIT_CONFIDENCE)
            .expect("emit")
    };
    let gap_emitted = emit(covering_ag, &ag_outcome);
    let ac_emitted = emit(covering_ac, &ac_outcome);
    let indel_emitted = emit(covering_indel, &indel_outcome);
    let tg_emitted = emit(covering_tg_snp, &tg_outcome);
    let qual_emitted = emit(covering_qual, &qual_outcome);
    let ca_emitted = emit(covering_ca, &ca_outcome);
    let target_emitted = emit(covering_target, &outcome);
    let mid_emitted = emit(covering_closed_mid, &mid_outcome);

    assert_eq!(
        info_i32(&rec_at(&gap_emitted, CLOSED_GT, "G").info, "DP"),
        Some(3)
    );
    assert_eq!(
        info_i32(&rec_at(&gap_emitted, CLOSED_AG, "A").info, "DP"),
        Some(3)
    );
    assert_eq!(
        info_i32(&rec_at(&ac_emitted, CLOSED_AC, "A").info, "DP"),
        Some(4)
    );
    assert_eq!(
        info_i32(&rec_at(&indel_emitted, CLOSED_INDEL, "TTC").info, "DP"),
        Some(1)
    );
    assert_eq!(
        info_i32(&rec_at(&tg_emitted, CLOSED_TG, "T").info, "DP"),
        Some(1)
    );
    let pin = rec_at(&qual_emitted, CLOSED_QUAL, "CT");
    assert_eq!(info_i32(&pin.info, "DP"), Some(2));
    assert!((pin.quality.expect("qual") - 31.60).abs() < 0.005);
    assert!((info_f64(&pin.info, "QD").expect("QD") - 15.80).abs() < 0.005);
    assert!(info_has(&pin.info, "BaseQRankSum"));
    assert!(info_has(&pin.info, "MQRankSum"));
    assert!(info_has(&pin.info, "ReadPosRankSum"));
    let ca_rec = rec_at(&ca_emitted, CLOSED_CA, "C");
    assert_eq!(info_i32(&ca_rec.info, "DP"), Some(6));
    assert_eq!(
        (info_f64(&ca_rec.info, "BaseQRankSum").expect("bq") * 1000.0).round(),
        -1834.0
    );
    assert_eq!(
        (info_f64(&ca_rec.info, "ReadPosRankSum").expect("rp") * 1000.0).round(),
        1282.0
    );
    assert_eq!(
        (info_f64(&ca_rec.info, "MQRankSum").expect("mqrs") * 1000.0).round(),
        1834.0
    );
    assert!((info_f64(&ca_rec.info, "MQ").expect("mq") - 40.58).abs() < 0.005);
    assert_eq!(
        info_i32(&rec_at(&mid_emitted, CLOSED_MID, "G").info, "DP"),
        Some(3)
    );

    let rec = rec_at(&target_emitted, TARGET, TARGET_REF);
    assert!((rec.quality.expect("QUAL") - 78.32).abs() < 0.005);
    assert_eq!(
        info_i32(&rec.info, "DP"),
        Some(2),
        "INFO DP is Coverage of the attached n=2 object"
    );
    let mq = info_f64(&rec.info, "MQ").expect("MQ");
    let sor = info_f64(&rec.info, "SOR").expect("SOR");
    kv("rust_info_dp", "2");
    kv("java_info_dp", "2");
    kv("rust_mq", format!("{mq:.2}"));
    kv("rust_sor", format!("{sor:.3}"));
    assert!(
        (mq - 47.00).abs() < 0.005,
        "MQ follows the n=2 object, got {mq}"
    );
    assert!(
        (sor - 2.303).abs() < 0.002,
        "SOR follows the n=2 object, got {sor}"
    );
    let sample = rec.samples.first().expect("sample");
    assert_eq!(
        sample.gt.as_ref().map(|g| g.to_string()).as_deref(),
        Some("1/1")
    );
    assert_eq!(sample.ad.as_deref(), Some(&[0u32, 2][..]));
    assert_eq!(sample.dp, Some(2));
    assert_eq!(sample.gq, Some(6.0));
    assert_eq!(sample.pl.as_deref(), Some(&[90u32, 6, 0][..]));
    kv(
        "first_arrow",
        "two-read hom-alt early-template attaches stored-haplotype retainEvidence n=2; Coverage consumes that object",
    );
}
