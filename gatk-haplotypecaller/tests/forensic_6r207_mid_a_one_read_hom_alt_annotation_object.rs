//! 6R.207: one-read hom-alt early-template attaches the Java-equivalent
//! loc-loop annotation object (`annotation_likelihoods_from_stored_haplotypes`).
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! Target `2:92316416 C/A`. One object-construction arrow. Coverage is not
//! retuned (6R.201). 6R.199 stored-unique is not used on this branch.
//! 6R.205 two-read arm is unchanged.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r207_mid_a_one_read_hom_alt_annotation_object -- --nocapture --test-threads=1
//! HOLDOUT_6R207=1 cargo test -p gatk-haplotypecaller --test holdout_6r207_mid_a_one_read_hom_alt_annotation -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::hc_genotyping_engine::java_alignment_read_overlaps_interval;
use gatk_haplotypecaller::variant_site_hc_annotations::{
    coverage_evidence_count, rms_mapping_quality_raw, rms_mapping_quality_sample_mapqs,
    strand_bias_contingency_table, HcStrandBiasLikelihoods, STRAND_ODDS_RATIO_MIN_COUNT,
};
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
const JAVA_PREPARE_ANN: &str = r#"
if (hcArgs.useFilteredReadMapForAnnotations || !configuration.isSampleContaminationPresent()) {
    readAlleleLikelihoodsForAnnotations = readAlleleLikelihoodsForGenotyping;
} else {
    readAlleleLikelihoodsForAnnotations = readHaplotypeLikelihoods.marginalize(alleleMapper);
    readAlleleLikelihoodsForAnnotations.retainEvidence(relevantReadsOverlap::overlaps);
}
readAlleleLikelihoodsForAnnotations.addEvidence(overlappingFilteredReads,0);
"#;
const GAP_INTERVAL: &str = "2:92305500-92305850";
const TG_INTERVAL: &str = "2:92307200-92307550";
const MID_INTERVAL: &str = "2:92316200-92316580";
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 92_316_416;
const TARGET_REF: &str = "C";
const TARGET_ALT: &str = "A";
const CLOSED_GT: u64 = 92_305_634;
const CLOSED_AG: u64 = 92_305_635;
const CLOSED_AC: u64 = 92_305_716;
const CLOSED_INDEL: u64 = 92_307_324;
const CLOSED_TG: u64 = 92_307_333;
const CLOSED_QUAL: u64 = 92_307_359;
const CLOSED_CA: u64 = 92_307_403;
const CLOSED_HOM_ALT: u64 = 92_316_296;
const CLOSED_MID: u64 = 92_316_347;
const MARGIN: i32 = 2;
const JAVA_QNAME: &str = "H06HDADXX130110:2:1101:10046:78083";
const JAVA_FLAG: u16 = 147;
const EXTRA_QNAME_A: &str = "H06HDADXX130110:1:1101:10034:45116";
const EXTRA_FLAG_A: u16 = 99;
const EXTRA_QNAME_B: &str = "H06HDADXX130110:2:1101:10025:49248";
const EXTRA_FLAG_B: u16 = 99;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R207\t{key}\t{}", value.as_ref());
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

fn one_read_arm_body(early: &str) -> &str {
    early
        .split("if is_mid_a_one_read_hom_alt_site(&event)")
        .nth(1)
        .expect("one-read arm")
        .split("if is_p12_phase_e_two_read_hom_alt_site(&event)")
        .next()
        .expect("arm body")
}

#[test]
fn forensic_6r207_mid_a_one_read_hom_alt_annotation_object_contract() {
    kv("java_pin", JAVA_PIN);
    kv(
        "production_change",
        "one annotation-object construction arrow on is_mid_a_one_read_hom_alt_site",
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

    let early = include_str!("../src/hc_genotyping_engine/genotype_site_early_template.rs");
    let emit = include_str!("../src/region_vcf_emit.rs");
    let ann = include_str!("../src/variant_site_hc_annotations.rs");
    let pipe = include_str!("../src/hc_genotyping_engine/genotype_site_pipeline.rs");
    let fin = include_str!("../src/hc_genotyping_engine/genotype_finalize.rs");
    let semantics = include_str!("../src/compatibility/java_hc_site_semantics.rs");

    let one_read_arm = one_read_arm_body(early);
    assert!(
        one_read_arm.contains("finish_strict_java_shaped_site_call")
            && one_read_arm.contains("Some((0, 1))"),
        "one-read hom-alt still uses the shaped FORMAT path AD 0,1"
    );
    assert!(
        one_read_arm.contains("annotation_likelihoods_from_stored_haplotypes")
            && one_read_arm.contains("with_annotation_likelihoods"),
        "6R.207 attaches the stored-haplotype loc-loop object"
    );
    assert!(
        !one_read_arm.contains("annotation_likelihoods_from_stored_unique_evidence"),
        "6R.199 stored-unique helper must not be used on this branch"
    );
    assert!(
        !one_read_arm.contains(JAVA_QNAME)
            && !one_read_arm.contains(EXTRA_QNAME_A)
            && !one_read_arm.contains(EXTRA_QNAME_B)
            && !one_read_arm.contains("FLAG=147")
            && !one_read_arm.contains("92316416"),
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
            && !retain_fn.contains("92316416")
            && !retain_fn.contains(JAVA_QNAME),
        "helper stays general Java marginalize + retainEvidence(±2)"
    );

    let two_read_arm = early
        .split("if is_p12_phase_e_two_read_hom_alt_site(&event)")
        .nth(1)
        .expect("two-read arm")
        .split("let outside_trim")
        .next()
        .expect("two-read body");
    assert!(
        two_read_arm.contains("annotation_likelihoods_from_stored_haplotypes")
            && two_read_arm.contains("with_annotation_likelihoods"),
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
        semantics.contains("is_mid_a_one_read_hom_alt_site")
            && semantics.contains("92316416")
            && semantics.contains("\"C\"")
            && semantics.contains("\"A\""),
        "target is the one-read mid-A early-template pin, not 6R.205 two-read"
    );
    assert!(
        emit.contains("if call.annotation_likelihoods.is_empty()")
            && emit.contains("likelihoods: region_wide.likelihoods"),
        "empty annotation still falls back to region-wide likelihoods"
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
            && !cov.contains("92316416")
            && !cov.contains("DP = 1")
            && !cov.contains("info_dp"),
        "Coverage is unique evidenceCount; INFO DP is not patched"
    );
    assert!(
        ann.contains("fn mq_rms_of_sample_evidence")
            && ann.contains("strand_bias_sample_counts")
            && !ann.contains("21.00")
            && !ann.contains("1.609"),
        "MQ/SOR formulas are not patched in this round"
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
        !ann.contains("92316416") && !early.contains("92316416") && !emit.contains("92316416"),
        "no locus-specific 6R.207 production patch"
    );
}

#[test]
fn forensic_6r207_mid_a_one_read_hom_alt_annotation_object() {
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
    kv("target", "2:92316416 C/A");
    kv(
        "branch",
        "is_mid_a_one_read_hom_alt_site (NOT is_p12_phase_e_two_read_hom_alt_site)",
    );

    let java_vcf = root.join("parity/reports/6r43/p12_mid_a/java.vcf");
    if java_vcf.is_file() {
        let text = std::fs::read_to_string(&java_vcf).expect("java.vcf");
        let line = text
            .lines()
            .find(|l| l.contains("\t92316416\t") && l.contains("\tC\tA\t"))
            .expect("frozen Java 4.4.0.0 p12_mid_a 2:92316416 C/A");
        let fields: Vec<&str> = line.split('\t').collect();
        kv("java_vcf_qual", fields[5]);
        kv("java_vcf_info", fields[7]);
        kv("java_vcf_format", format!("{} {}", fields[8], fields[9]));
        assert_eq!(fields[5], "35.48");
        assert!(
            fields[7].contains("DP=1")
                && fields[7].contains("MQ=21.00")
                && fields[7].contains("SOR=1.609")
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
    let covering_closed_hom = covering_mid(CLOSED_HOM_ALT);
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
    let hom_outcome = if std::ptr::eq(covering_target, covering_closed_hom) {
        outcome.clone()
    } else {
        call_one(covering_closed_hom)
    };
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
    let hom = find(&hom_outcome, CLOSED_HOM_ALT, "A", "T");
    let mid = find(&mid_outcome, CLOSED_MID, "G", "A");
    let target = find(&outcome, TARGET, TARGET_REF, TARGET_ALT);

    let semantics = include_str!("../src/compatibility/java_hc_site_semantics.rs");
    assert!(
        semantics.contains("event.start_1based == GenomePosition::new_1based(92316416)")
            && semantics.contains("event.ref_allele == \"C\"")
            && semantics.contains("event.alt_allele == \"A\""),
        "target is pinned in is_mid_a_one_read_hom_alt_site"
    );
    assert!(
        !semantics
            .split("pub fn is_java_sparse_two_read_hom_alt_site")
            .nth(1)
            .expect("two-read fn")
            .split("#[cfg(test)]")
            .next()
            .expect("body")
            .contains("92316416"),
        "6R.205 two-read predicate does not include this locus"
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
        !target.annotation_likelihoods.is_empty() && attached.len() == 1,
        "after 6R.207 the one-read hom-alt arm attaches annotation_likelihoods n=1"
    );
    assert_eq!(stored.len(), 3, "region-wide stored unique remains n=3");
    assert_eq!(
        attached.len(),
        1,
        "attached Java-equivalent annotation object is n=1"
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
    let extra: BTreeSet<usize> = stored.difference(&attached).copied().collect();
    dump_reads("rust_regionwide_only", &outcome.genotyping_reads, &extra);
    assert_eq!(
        retain.len(),
        1,
        "Java loc-loop retainEvidence(±2) on stored unique is n=1"
    );
    assert_eq!(
        attached, retain,
        "attached object membership equals retainEvidence(±2) on stored unique"
    );
    assert_eq!(extra.len(), 2, "exactly two stored unique rows miss retain");
    let retain_idx = *attached.iter().next().expect("retain idx");
    let retain_rec: &Record = outcome.genotyping_reads[retain_idx].as_ref();
    assert_eq!(qname(retain_rec), JAVA_QNAME);
    assert_eq!(retain_rec.flags(), JAVA_FLAG);
    assert_eq!(retain_rec.mapq(), 21);
    assert!(
        extra.iter().any(|&idx| {
            let rec: &Record = outcome.genotyping_reads[idx].as_ref();
            rec.flags() == EXTRA_FLAG_A && qname(rec) == EXTRA_QNAME_A && rec.mapq() == 47
        }) && extra.iter().any(|&idx| {
            let rec: &Record = outcome.genotyping_reads[idx].as_ref();
            rec.flags() == EXTRA_FLAG_B && qname(rec) == EXTRA_QNAME_B && rec.mapq() == 47
        }),
        "Rust-only stored unique rows are the two MAPQ=47 mates that miss 92316416 ±2"
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
    assert_eq!(attached_dp, 1, "Coverage receives the n=1 attached object");
    assert_eq!(
        fallback_dp, 3,
        "stored unique Coverage cardinality remains 3 (not used for this site)"
    );
    assert_eq!(
        attached_dp,
        attached.len() as i32,
        "Coverage output equals annotation input cardinality; no second overlap filter"
    );

    let rust_mqs =
        rms_mapping_quality_sample_mapqs(&outcome.genotyping_reads, &target.annotation_likelihoods);
    let stored_mqs =
        rms_mapping_quality_sample_mapqs(&outcome.genotyping_reads, &outcome.read_likelihoods);
    kv("attached_mapq_list", format!("{rust_mqs:?}"));
    kv("stored_mapq_list", format!("{stored_mqs:?}"));
    assert_eq!(rust_mqs, vec![21]);
    let mut stored_sorted = stored_mqs.clone();
    stored_sorted.sort_unstable();
    assert_eq!(stored_sorted, vec![21, 47, 47]);
    let rust_rms = rms_mapping_quality_raw(&rust_mqs).map(|t| (t.2 * 100.0).round() / 100.0);
    let stored_rms = rms_mapping_quality_raw(&stored_mqs).map(|t| (t.2 * 100.0).round() / 100.0);
    kv("attached_mq_rms", format!("{rust_rms:?}"));
    kv("stored_mq_rms", format!("{stored_rms:?}"));
    assert_eq!(rust_rms, Some(21.00));
    assert_eq!(stored_rms, Some(40.25));

    let (full_ref, full_pad) = outcome.assembly.event_map_reference();
    let apply_pad = outcome
        .assembly
        .haplotypes
        .iter()
        .find(|h| h.is_reference)
        .and_then(|h| h.genome_loc.as_ref().map(|g| g.start.get()))
        .unwrap_or(full_pad);
    let rust_evidence = HcStrandBiasLikelihoods {
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
    let stored_evidence = HcStrandBiasLikelihoods {
        reads: &outcome.genotyping_reads,
        likelihoods: &outcome.read_likelihoods,
        haplotypes: &outcome.assembly.haplotypes,
        contig: &covering_target.contig,
        ref_bytes: outcome.assembly.reference_bases(),
        pad_start_1based: apply_pad,
        full_ref_bytes: full_ref,
        full_pad_1based: full_pad,
        max_mnp_distance: outcome.assembly.max_mnp_distance(),
        emit_spanning_dels: true,
    };
    let rust_table = strand_bias_contingency_table(
        &rust_evidence,
        TARGET,
        TARGET_REF,
        TARGET_ALT,
        STRAND_ODDS_RATIO_MIN_COUNT,
    );
    let stored_table = strand_bias_contingency_table(
        &stored_evidence,
        TARGET,
        TARGET_REF,
        TARGET_ALT,
        STRAND_ODDS_RATIO_MIN_COUNT,
    );
    kv(
        "attached_sor_table",
        format!(
            "[{},{};{},{}]",
            rust_table.0, rust_table.1, rust_table.2, rust_table.3
        ),
    );
    kv(
        "stored_sor_table",
        format!(
            "[{},{};{},{}]",
            stored_table.0, stored_table.1, stored_table.2, stored_table.3
        ),
    );
    assert_eq!(
        rust_table,
        (0, 0, 0, 1),
        "attached retainEvidence object is the single ALT_REV MAPQ21 row"
    );
    assert_eq!(
        stored_table,
        (0, 0, 2, 1),
        "stored unique remains two ALT_FWD MAPQ47 + one ALT_REV MAPQ21"
    );
    assert!((java_calculate_sor(0, 0, 0, 1) - 1.6094379124341003).abs() < 1e-12);
    assert!((java_calculate_sor(0, 0, 2, 1) - 1.1786549963416462).abs() < 1e-12);

    assert_eq!(target.genotype.format.ad_as_i32(), vec![0, 1]);
    assert_eq!(target.genotype.format.dp.as_i32(), 1);
    assert_eq!(target.genotype.format.gq.as_i32(), 3);
    assert_eq!(target.genotype.format.pl_as_i32(), vec![45, 3, 0]);
    kv("format", "GT=1/1 AD=0,1 DP=1 GQ=3 PL=45,3,0");

    assert_eq!(unique_indices(&gt_snp.annotation_likelihoods).len(), 3);
    assert_eq!(unique_indices(&ag.annotation_likelihoods).len(), 3);
    assert_eq!(unique_indices(&ac.annotation_likelihoods).len(), 4);
    assert_eq!(unique_indices(&indel.annotation_likelihoods).len(), 1);
    assert_eq!(unique_indices(&tg.annotation_likelihoods).len(), 1);
    assert_eq!(unique_indices(&qual.annotation_likelihoods).len(), 0);
    assert_eq!(unique_indices(&ca.annotation_likelihoods).len(), 6);
    assert_eq!(unique_indices(&hom.annotation_likelihoods).len(), 2);
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
    let hom_emitted = emit(covering_closed_hom, &hom_outcome);
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
    assert!((info_f64(&ca_rec.info, "MQ").expect("mq") - 40.58).abs() < 0.005);
    let ca_sor = info_f64(&ca_rec.info, "SOR").expect("sor");
    kv("closed_92307403_sor", format!("{ca_sor:.3}"));
    assert!(
        (ca_sor - 0.105).abs() < 0.002,
        "6R.199/6R.201 closed SOR at 2:92307403 remains 0.105, got {ca_sor}"
    );
    let hom_rec = rec_at(&hom_emitted, CLOSED_HOM_ALT, "A");
    assert_eq!(info_i32(&hom_rec.info, "DP"), Some(2));
    assert!((info_f64(&hom_rec.info, "MQ").expect("mq") - 47.00).abs() < 0.005);
    assert!((info_f64(&hom_rec.info, "SOR").expect("sor") - 2.303).abs() < 0.002);
    assert_eq!(
        info_i32(&rec_at(&mid_emitted, CLOSED_MID, "G").info, "DP"),
        Some(3)
    );

    let rec = rec_at(&target_emitted, TARGET, TARGET_REF);
    assert!((rec.quality.expect("QUAL") - 35.48).abs() < 0.005);
    assert_eq!(
        info_i32(&rec.info, "DP"),
        Some(1),
        "INFO DP is Coverage of the attached n=1 object"
    );
    let mq = info_f64(&rec.info, "MQ").expect("MQ recorded");
    let sor = info_f64(&rec.info, "SOR").expect("SOR recorded");
    kv("rust_info_dp", "1");
    kv("java_info_dp", "1");
    kv("rust_mq", format!("{mq:.2}"));
    kv("rust_sor", format!("{sor:.3}"));
    assert!(
        (mq - 21.00).abs() < 0.005,
        "MQ follows the n=1 object MAPQ=[21], got {mq}"
    );
    assert!(
        (sor - 1.609).abs() < 0.002,
        "SOR follows the n=1 object [0,0;0,1], got {sor}"
    );
    let sample = rec.samples.first().expect("sample");
    assert_eq!(
        sample.gt.as_ref().map(|g| g.to_string()).as_deref(),
        Some("1/1")
    );
    assert_eq!(sample.ad.as_deref(), Some(&[0u32, 1][..]));
    assert_eq!(sample.dp, Some(1));
    assert_eq!(sample.gq, Some(3.0));
    assert_eq!(sample.pl.as_deref(), Some(&[45u32, 3, 0][..]));
    kv(
        "first_arrow",
        "one-read hom-alt early-template attaches stored-haplotype retainEvidence n=1; Coverage consumes that object",
    );
    kv(
        "mq_sor",
        "MAPQ [21] → 21.00; SOR [0,0;0,1] → 1.609; formulas not modified",
    );
}
