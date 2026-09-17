//! 6R.209: loc-loop retainEvidence keeps both same-QNAME mates at `2:92317399 C/A`.
//! One production arrow: do not apply the 6R.180 same-QNAME collapse on the
//! Java loc-loop annotation object (`annotation_likelihoods_from_stored_haplotypes`).
//! The 6R.180 collapse stays on FORMAT-subset / stored-unique callers.
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r209_preserve_loc_loop_mate_evidence -- --nocapture --test-threads=1
//! HOLDOUT_6R209=1 cargo test -p gatk-haplotypecaller --test holdout_6r209_preserve_loc_loop_mate_evidence -- --nocapture --test-threads=1
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
const MID_B_INTERVAL: &str = "2:92317000-92319000";
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 92_317_399;
const TARGET_REF: &str = "C";
const TARGET_ALT: &str = "A";
const JAVA_QNAME: &str = "H06HDADXX130110:1:1101:10040:45938";
const JAVA_FWD_FLAG: u16 = 99;
const JAVA_REV_FLAG: u16 = 147;
const CLOSED_GT: u64 = 92_305_634;
const CLOSED_AG: u64 = 92_305_635;
const CLOSED_AC: u64 = 92_305_716;
const CLOSED_INDEL: u64 = 92_307_324;
const CLOSED_TG: u64 = 92_307_333;
const CLOSED_QUAL: u64 = 92_307_359;
const CLOSED_CA: u64 = 92_307_403;
const CLOSED_HOM_ALT: u64 = 92_316_296;
const CLOSED_ONE_READ: u64 = 92_316_416;
const CLOSED_MID: u64 = 92_316_347;
const MARGIN: i32 = 2;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R209\t{key}\t{}", value.as_ref());
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
                "idx={idx} qname={} flag={} mapq={} contig={} start={} end={} cigar={} strand={} overlap_pm2={}",
                qname(rec),
                rec.flags(),
                rec.mapq(),
                rec.tid(),
                rec.pos() + 1,
                alignment_end_1based(rec),
                rec.cigar(),
                if rec.is_reverse() { "REV" } else { "FWD" },
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

#[test]
fn forensic_6r209_preserve_loc_loop_mate_evidence_contract() {
    kv("java_pin", JAVA_PIN);
    kv(
        "production_change",
        "one lifecycle-scoped removal of same-QNAME collapse from the Java loc-loop annotation object",
    );
    kv(
        "classification",
        "ANNOTATION_SOURCE_OBJECT_DIVERGENCE / loc-loop retainEvidence preserves same-QNAME mates",
    );
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
    let semantics = include_str!("../src/compatibility/java_hc_site_semantics.rs");

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
    assert!(
        pipe.contains("let annotation = annotation_likelihoods_from_stored_haplotypes")
            && pipe.contains("with_annotation_likelihoods(annotation)"),
        "6R.189 SiteScore attach stays (default loc-loop object)"
    );
    assert!(
        !semantics
            .split("pub fn is_java_sparse_two_read_hom_alt_site")
            .nth(1)
            .expect("two-read fn")
            .split("#[cfg(test)]")
            .next()
            .expect("body")
            .contains("92317399"),
        "6R.205 two-read predicate does not include this locus"
    );
    assert!(
        !semantics
            .split("pub fn is_mid_a_one_read_hom_alt_site")
            .nth(1)
            .expect("one-read fn")
            .split("pub fn is_mid_a_two_read_hom_alt_site")
            .next()
            .expect("body")
            .contains("92317399"),
        "6R.207 one-read predicate does not include this locus"
    );
    assert!(
        semantics.contains("MID_B_DENSE_CLUSTER_START: u64 = 92317399")
            && semantics.contains("is_mid_b_java_sparse_snp"),
        "target is the canonical mid-B sparse SNP start, not mid-A"
    );
    assert!(
        !early.contains("92317399") && !ann.contains("92317399") && !emit.contains("92317399"),
        "no locus-specific 6R.209 production patch"
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
            && !ann.contains("DP = 2"),
        "MQ/SOR/Coverage formulas are not patched in this round"
    );
    assert!(
        emit.contains("if call.annotation_likelihoods.is_empty()")
            && emit.contains("likelihoods: &call.annotation_likelihoods"),
        "emit still consumes the attached object when present"
    );
    let helper = early
        .split("fn annotation_likelihoods_from_stored_haplotypes(")
        .nth(1)
        .expect("stored-hap helper")
        .split("fn annotation_likelihoods_from_stored_unique_evidence")
        .next()
        .expect("helper body");
    assert!(
        helper.contains("likelihood_subset_for_event")
            && !helper.contains("per_variant_annotation_likelihoods(&subset, reads)"),
        "loc-loop helper retainEvidence membership is not QNAME-collapsed"
    );
    let unique_fn = early
        .split("fn annotation_likelihoods_from_stored_unique_evidence")
        .nth(1)
        .expect("unique helper");
    assert!(
        unique_fn.contains("per_variant_annotation_likelihoods(likelihoods, reads)"),
        "6R.199 stored-unique still uses the 6R.180 collapse helper"
    );
    assert!(
        pipe.contains("per_variant_annotation_likelihoods(&subset, likelihood_reads)"),
        "6R.180 FORMAT-subset collapse stays on the specialized pipeline caller"
    );
    let ge = include_str!("../src/hc_genotyping_engine/mod.rs");
    assert!(
        ge.contains("if qnames.len() == 1 && best_ll.len() > 1"),
        "6R.180 same-QNAME collapse is not deleted globally"
    );
    assert!(
        !early
            .split("if gap_sparse_read_genotype")
            .nth(1)
            .unwrap_or("")
            .contains("92317399"),
        "no mid-B coordinate special-case in later try_shaped arms"
    );
}

#[test]
fn forensic_6r209_preserve_loc_loop_mate_evidence() {
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
        "one lifecycle-scoped removal of same-QNAME collapse from the Java loc-loop annotation object",
    );
    kv("target", "2:92317399 C/A");

    let java_vcf = root.join("parity/reports/6r43/ctrl_mid_b/java.vcf");
    if java_vcf.is_file() {
        let text = std::fs::read_to_string(&java_vcf).expect("java.vcf");
        let line = text
            .lines()
            .find(|l| l.contains("\t92317399\t") && l.contains("\tC\tA\t"))
            .expect("frozen Java 4.4.0.0 ctrl_mid_b 2:92317399 C/A");
        let fields: Vec<&str> = line.split('\t').collect();
        kv("java_vcf_qual", fields[5]);
        kv("java_vcf_info", fields[7]);
        kv("java_vcf_format", format!("{} {}", fields[8], fields[9]));
        assert_eq!(fields[5], "78.32");
        assert!(
            fields[7].contains("DP=2")
                && fields[7].contains("MQ=27.00")
                && fields[7].contains("SOR=0.693")
        );
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
    let mid_b_regions = flatten_assembly_regions(&walk(MID_B_INTERVAL));
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
    let covering_mid_b = mid_b_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= TARGET
                && r.end.get() >= TARGET
        })
        .expect("mid-B covering 92317399");

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
    let covering_hom = covering_mid(CLOSED_HOM_ALT);
    let covering_one = covering_mid(CLOSED_ONE_READ);
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
    let hom_outcome = call_one(covering_hom);
    let one_outcome = if std::ptr::eq(covering_hom, covering_one) {
        hom_outcome.clone()
    } else {
        call_one(covering_one)
    };
    let mid_outcome = if std::ptr::eq(covering_hom, covering_closed_mid) {
        hom_outcome.clone()
    } else {
        call_one(covering_closed_mid)
    };
    let outcome = call_one(covering_mid_b);

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
    let one = find(&one_outcome, CLOSED_ONE_READ, "C", "A");
    let mid = find(&mid_outcome, CLOSED_MID, "G", "A");
    let target = find(&outcome, TARGET, TARGET_REF, TARGET_ALT);

    kv(
        "branch",
        "try_shaped=None; SiteScore 6R.189 annotation_likelihoods_from_stored_haplotypes",
    );
    kv(
        "active_region",
        format!(
            "{}-{}",
            covering_mid_b.start.get(),
            covering_mid_b.end.get()
        ),
    );

    let stored = unique_indices(&outcome.read_likelihoods);
    let attached = unique_indices(&target.annotation_likelihoods);
    kv("stored_hap_n", stored.len().to_string());
    kv("attached_annotation_n", attached.len().to_string());
    kv(
        "attached_cells",
        target.annotation_likelihoods.len().to_string(),
    );
    kv(
        "format_ad",
        format!("{:?}", target.genotype.format.ad_as_i32()),
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
    dump_reads("stored_not_attached", &outcome.genotyping_reads, &extra);

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

    let attached_mqs =
        rms_mapping_quality_sample_mapqs(&outcome.genotyping_reads, &target.annotation_likelihoods);
    let stored_mqs =
        rms_mapping_quality_sample_mapqs(&outcome.genotyping_reads, &outcome.read_likelihoods);
    kv("attached_mapq_list", format!("{attached_mqs:?}"));
    kv("stored_mapq_list", format!("{stored_mqs:?}"));
    kv(
        "attached_mq_rms",
        format!(
            "{:?}",
            rms_mapping_quality_raw(&attached_mqs).map(|t| (t.2 * 100.0).round() / 100.0)
        ),
    );

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
        contig: &covering_mid_b.contig,
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
        contig: &covering_mid_b.contig,
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
            attached_table.0, attached_table.1, attached_table.2, attached_table.3
        ),
    );
    kv(
        "stored_sor_table",
        format!(
            "[{},{};{},{}]",
            stored_table.0, stored_table.1, stored_table.2, stored_table.3
        ),
    );
    kv(
        "attached_sor_formula",
        format!(
            "{:.3}",
            java_calculate_sor(
                attached_table.0,
                attached_table.1,
                attached_table.2,
                attached_table.3
            )
        ),
    );

    assert_eq!(target.genotype.format.ad_as_i32(), vec![0, 2]);
    assert_eq!(target.genotype.format.dp.as_i32(), 2);
    assert_eq!(target.genotype.format.gq.as_i32(), 6);
    assert_eq!(target.genotype.format.pl_as_i32(), vec![90, 6, 0]);
    kv("format", "GT=1/1 AD=0,2 DP=2 GQ=6 PL=90,6,0");

    assert_eq!(unique_indices(&gt_snp.annotation_likelihoods).len(), 3);
    assert_eq!(unique_indices(&ag.annotation_likelihoods).len(), 3);
    assert_eq!(unique_indices(&ac.annotation_likelihoods).len(), 4);
    assert_eq!(unique_indices(&indel.annotation_likelihoods).len(), 1);
    assert_eq!(unique_indices(&tg.annotation_likelihoods).len(), 1);
    assert_eq!(unique_indices(&qual.annotation_likelihoods).len(), 0);
    assert_eq!(unique_indices(&ca.annotation_likelihoods).len(), 6);
    assert_eq!(unique_indices(&hom.annotation_likelihoods).len(), 2);
    assert_eq!(unique_indices(&one.annotation_likelihoods).len(), 1);
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
    let hom_emitted = emit(covering_hom, &hom_outcome);
    let one_emitted = emit(covering_one, &one_outcome);
    let mid_emitted = emit(covering_closed_mid, &mid_outcome);
    let target_emitted = emit(covering_mid_b, &outcome);

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
    let tg_rec = rec_at(&tg_emitted, CLOSED_TG, "T");
    assert!((info_f64(&tg_rec.info, "MQ").expect("tg mq") - 44.00).abs() < 0.005);
    assert!((info_f64(&tg_rec.info, "SOR").expect("tg sor") - 1.609).abs() < 0.002);
    let pin = rec_at(&qual_emitted, CLOSED_QUAL, "CT");
    assert_eq!(info_i32(&pin.info, "DP"), Some(2));
    assert!((pin.quality.expect("qual") - 31.60).abs() < 0.005);
    assert!(info_has(&pin.info, "BaseQRankSum"));
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
    assert_eq!(
        info_i32(&rec_at(&mid_emitted, CLOSED_MID, "G").info, "DP"),
        Some(3)
    );

    let rec = rec_at(&target_emitted, TARGET, TARGET_REF);
    assert!((rec.quality.expect("QUAL") - 78.32).abs() < 0.005);
    let dp = info_i32(&rec.info, "DP").expect("DP");
    let mq = info_f64(&rec.info, "MQ").expect("MQ");
    let sor = info_f64(&rec.info, "SOR").expect("SOR");
    kv("rust_info_dp", dp.to_string());
    kv("java_info_dp", "2");
    kv("rust_mq", format!("{mq:.2}"));
    kv("java_mq", "27.00");
    kv("rust_sor", format!("{sor:.3}"));
    kv("java_sor", "0.693");
    let sample = rec.samples.first().expect("sample");
    assert_eq!(
        sample.gt.as_ref().map(|g| g.to_string()).as_deref(),
        Some("1/1")
    );
    assert_eq!(sample.ad.as_deref(), Some(&[0u32, 2][..]));
    assert_eq!(sample.dp, Some(2));
    assert_eq!(sample.gq, Some(6.0));
    assert_eq!(sample.pl.as_deref(), Some(&[90u32, 6, 0][..]));

    assert_eq!(stored.len(), 2, "A stored unique n=2");
    assert_eq!(retain.len(), 2, "Java loc-loop retainEvidence(±2) n=2");
    assert_eq!(
        attached.len(),
        2,
        "C loc-loop annotation preserves retainEvidence n=2"
    );
    assert_eq!(attached, retain, "attached membership is retainEvidence");
    assert!(
        !target.annotation_likelihoods.is_empty(),
        "D region-wide fallback is not taken; annotation is attached"
    );
    assert_eq!(attached_dp, 2, "E Coverage consumes attached n=2");
    assert_eq!(stored_dp, 2, "stored unique Coverage n=2");
    assert_eq!(
        attached_dp,
        attached.len() as i32,
        "Coverage does not apply a second ±2 filter"
    );
    assert_eq!(dp, 2, "INFO DP from attached n=2");
    assert!((mq - 27.00).abs() < 0.005, "MQ RMS of [27,27] is 27.00");
    assert!((sor - 0.693).abs() < 0.002, "INFO SOR from [0,0;1,1]");
    assert_eq!(
        attached_table,
        (0, 0, 1, 1),
        "F SOR from both FWD and REV mates"
    );
    assert_eq!(stored_table, (0, 0, 1, 1));
    let java_sor = java_calculate_sor(
        attached_table.0,
        attached_table.1,
        attached_table.2,
        attached_table.3,
    );
    assert!(
        (java_sor - 0.693).abs() < 0.002,
        "SOR formula on the attached table"
    );
    assert!(
        (java_sor - sor).abs() < 0.002,
        "INFO SOR is calculateSOR of the attached table, not a formula write"
    );

    let rec_meta = |idx: usize| {
        let rec: &Record = outcome.genotyping_reads[idx].as_ref();
        (qname(rec), rec.flags(), rec.mapq(), rec.is_reverse())
    };
    let retain_meta: Vec<_> = retain.iter().copied().map(rec_meta).collect();
    let attached_meta: Vec<_> = attached.iter().copied().map(rec_meta).collect();
    let qname_n = attached_meta
        .iter()
        .map(|m| m.0.clone())
        .collect::<BTreeSet<_>>()
        .len();
    assert_eq!(
        qname_n, 1,
        "same-QNAME cardinality is 1 and must not collapse"
    );
    assert!(
        retain_meta.iter().all(|m| m.0 == JAVA_QNAME && m.2 == 27),
        "retainEvidence is the mate pair MAPQ=27"
    );
    assert!(
        retain_meta.iter().any(|m| m.1 == JAVA_FWD_FLAG && !m.3)
            && retain_meta.iter().any(|m| m.1 == JAVA_REV_FLAG && m.3),
        "retainEvidence keeps both FLAG=99 FWD and FLAG=147 REV"
    );
    assert_eq!(attached_meta.len(), 2);
    assert!(
        attached_meta.iter().any(|m| m.1 == JAVA_FWD_FLAG && !m.3)
            && attached_meta.iter().any(|m| m.1 == JAVA_REV_FLAG && m.3),
        "both GATKRead mates survive on the loc-loop annotation object"
    );
    assert!(
        extra.is_empty(),
        "no retained row is dropped after collapse"
    );
    assert_eq!(target.genotype.format.ad_as_i32(), vec![0, 2]);
    assert_eq!(
        unique_indices(&tg.annotation_likelihoods).len(),
        1,
        "6R.186 cluster-TG stays n=1"
    );

    kv(
        "first_arrow",
        "loc-loop retainEvidence n=2 same-QNAME mates preserved; 6R.180 collapse not applied",
    );
    kv(
        "classification",
        "ANNOTATION_SOURCE_OBJECT_DIVERGENCE / loc-loop unique GATKRead membership",
    );
}
