//! 6R.198: proof-only. Per-variant annotation object lifecycle at `2:92307403 C/A`.
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! PRODUCTION CHANGE: NONE.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r198_per_variant_annotation_object -- --nocapture --test-threads=1
//! HOLDOUT_6R198=1 cargo test -p gatk-haplotypecaller --test holdout_6r198_per_variant_object -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::{build_per_haplotype_variation_events, VariationEvent};
use gatk_haplotypecaller::hc_allele_mapping::create_allele_mapper_with_events;
use gatk_haplotypecaller::hc_genotyping_engine::java_alignment_read_overlaps_interval;
use gatk_haplotypecaller::variant_site_hc_annotations::{
    baseq_rank_sum_quals_from_likelihoods, coverage_evidence_count,
    mq_rank_sum_quals_from_likelihoods, read_pos_rank_sum_quals_from_likelihoods,
    rms_mapping_quality_raw, rms_mapping_quality_sample_mapqs, HcStrandBiasLikelihoods,
};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, marginalize_rows_to_biallelic_alleles,
    region_likelihoods_to_rows, traverse_assembly_region_walker, try_emit_call_region_variants,
    AssemblyRegion, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    HaplotypeCallerEngine, ReadFilterParams, RegionReadLikelihood, SharedBamRecord,
    WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
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

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R198\t{key}\t{}", value.as_ref());
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
        kv(
            label,
            format!(
                "idx={idx} qname={} flag={} mapq={} start={} end={} cigar={}",
                qname(rec),
                rec.flags(),
                rec.mapq(),
                rec.pos() + 1,
                alignment_end_1based(rec),
                rec.cigar()
            ),
        );
    }
}

#[test]
fn forensic_6r198_no_production_change() {
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");
    kv(
        "classification",
        "E — DIFFERENT BRANCH / SPECIAL-PATH SEMANTICS (cluster-downstream early-template never constructs the loc-loop annotation object)",
    );

    let early = include_str!("../src/hc_genotyping_engine/genotype_site_early_template.rs");
    let pipe = include_str!("../src/hc_genotyping_engine/genotype_site_pipeline.rs");
    let emit = include_str!("../src/region_vcf_emit.rs");
    let ann = include_str!("../src/variant_site_hc_annotations.rs");

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
        "6R.199 cluster-downstream attaches stored unique evidence, not 6R.186 retain"
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
    assert!(
        emit.contains("if call.annotation_likelihoods.is_empty()")
            && emit.contains("likelihoods: region_wide.likelihoods"),
        "6R.180 empty-object fallback stays; 6R.198 does not retarget it"
    );
    assert!(
        ann.contains("coverage_evidence_count")
            && ann.contains("baseq_rank_sum_quals_from_likelihoods")
            && ann.contains("mq_rank_sum_quals_from_likelihoods"),
        "DP/RankSum formulas stay closed"
    );
    assert!(
        !early.contains("92307403") && !pipe.contains("92307403") && !ann.contains("92307403"),
        "no locus-specific 6R.198 production patch"
    );
}

#[test]
fn forensic_6r198_object_lifecycle() {
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing P12 ref/BAM");
        return;
    }
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");
    kv("target", "2:92307403 C/A");

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
    let covering_for = |pos: u64| {
        tg_regions
            .iter()
            .find(|r| {
                matches!(
                    call_disposition(r),
                    AssemblyRegionCallDisposition::ActiveFull
                ) && r.start.get() <= pos
                    && r.end.get() >= pos
            })
            .unwrap_or_else(|| panic!("covering {pos}"))
    };
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
    let covering_ag = covering_gap(CLOSED_AG);
    let covering_ac = covering_gap(CLOSED_AC);
    let covering_indel = covering_for(CLOSED_INDEL);
    let covering_tg = covering_for(CLOSED_TG);
    let covering_qual = covering_for(CLOSED_QUAL);
    let covering = covering_for(TARGET);

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
    let ag_outcome = call(covering_ag);
    let ac_outcome = if std::ptr::eq(covering_ag, covering_ac) {
        None
    } else {
        Some(call(covering_ac))
    };
    let ac_src = ac_outcome.as_ref().unwrap_or(&ag_outcome);
    let indel_outcome = call(covering_indel);
    let tg_snp_outcome = call(covering_tg);
    let qual_outcome = call(covering_qual);
    let outcome = call(covering);

    let gt_snp = ag_outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(CLOSED_GT)
                && c.event.ref_allele == "G"
                && c.event.alt_allele == "T"
        })
        .expect("G/T");
    let ag = ag_outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(CLOSED_AG)
                && c.event.ref_allele == "A"
                && c.event.alt_allele == "G"
        })
        .expect("A/G");
    let ac = ac_src
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(CLOSED_AC)
                && c.event.ref_allele == "A"
                && c.event.alt_allele == "C"
        })
        .expect("A/C");
    let indel = indel_outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(CLOSED_INDEL)
                && c.event.ref_allele == "TTC"
                && c.event.alt_allele == "T"
        })
        .expect("TTC/T");
    let tg = tg_snp_outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(CLOSED_TG)
                && c.event.ref_allele == "T"
                && c.event.alt_allele == "G"
        })
        .expect("T/G");
    let qual = qual_outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(CLOSED_QUAL)
                && c.event.ref_allele == "CT"
                && c.event.alt_allele == "C"
        })
        .expect("CT/C");
    let target = outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(TARGET)
                && c.event.ref_allele == TARGET_REF
                && c.event.alt_allele == TARGET_ALT
        })
        .expect("C/A");

    assert_eq!(
        target.event.start_1based,
        GenomePosition::new_1based(TARGET)
    );
    assert_eq!(target.event.ref_allele, TARGET_REF);
    assert_eq!(target.event.alt_allele, TARGET_ALT);
    assert_eq!(tg.event.start_1based, GenomePosition::new_1based(CLOSED_TG));
    assert_eq!(tg.event.ref_allele, "T");
    assert_eq!(tg.event.alt_allele, "G");
    let sem = include_str!("../src/compatibility/java_hc_site_semantics.rs");
    assert!(
        sem.contains("(92307403, \"C\", \"A\")"),
        "target is CLUSTER_DOWNSTREAM_SNPS"
    );
    assert!(
        sem.contains("CLUSTER_TG_SNP_START: u64 = 92307333"),
        "control T/G remains cluster-TG"
    );

    let gt_n = unique_indices(&gt_snp.annotation_likelihoods);
    let ag_n = unique_indices(&ag.annotation_likelihoods);
    let ac_n = unique_indices(&ac.annotation_likelihoods);
    let indel_n = unique_indices(&indel.annotation_likelihoods);
    let tg_n = unique_indices(&tg.annotation_likelihoods);
    let qual_n = unique_indices(&qual.annotation_likelihoods);
    let attached = unique_indices(&target.annotation_likelihoods);
    kv("control_92305634_ann_n", gt_n.len().to_string());
    kv("control_92305635_ann_n", ag_n.len().to_string());
    kv("control_92305716_ann_n", ac_n.len().to_string());
    kv("control_92307324_ann_n", indel_n.len().to_string());
    kv("control_92307333_ann_n", tg_n.len().to_string());
    kv("control_92307359_ann_n", qual_n.len().to_string());
    kv("target_attached_n", attached.len().to_string());
    assert_eq!(
        gt_n.len(),
        3,
        "neighbor G/T two-read arm attaches retainEvidence n=3 (6R.205)"
    );
    assert_eq!(ag_n.len(), 3, "6R.189 SiteScore attach stays");
    assert_eq!(ac_n.len(), 4, "6R.191 A/C annotation object stays");
    assert_eq!(indel_n.len(), 1, "6R.180 TTC/T attach stays");
    assert_eq!(tg_n.len(), 1, "6R.186 cluster-TG attach stays");
    assert_eq!(qual_n.len(), 0, "CT/C early-template still empty");
    assert_eq!(
        attached.len(),
        6,
        "6R.199 cluster-downstream attaches stored unique n=6"
    );

    let stored = unique_indices(&outcome.read_likelihoods);
    kv("stored_hap_n", stored.len().to_string());
    kv(
        "stored_hap_columns",
        outcome.assembly.haplotypes.len().to_string(),
    );
    kv(
        "genotyping_reads_n",
        outcome.genotyping_reads.len().to_string(),
    );
    kv("covering_reads_n", covering.reads.len().to_string());

    let mut stored_retain = BTreeSet::new();
    for idx in &stored {
        if java_alignment_read_overlaps_interval(
            &outcome.genotyping_reads[*idx],
            TARGET,
            TARGET,
            MARGIN,
        ) {
            stored_retain.insert(*idx);
        }
    }
    let mut geno_overlap = BTreeSet::new();
    for (idx, rec) in outcome.genotyping_reads.iter().enumerate() {
        if java_alignment_read_overlaps_interval(rec, TARGET, TARGET, MARGIN) {
            geno_overlap.insert(idx);
        }
    }
    let mut covering_overlap = BTreeSet::new();
    for (idx, rec) in covering.reads.iter().enumerate() {
        if java_alignment_read_overlaps_interval(rec, TARGET, TARGET, MARGIN) {
            covering_overlap.insert(idx);
        }
    }
    kv("stored_retainEvidence_n", stored_retain.len().to_string());
    kv("genotyping_reads_overlap_n", geno_overlap.len().to_string());
    kv(
        "covering_reads_overlap_n",
        covering_overlap.len().to_string(),
    );
    dump_reads("stored_all", &outcome.genotyping_reads, &stored);
    dump_reads("stored_retain", &outcome.genotyping_reads, &stored_retain);
    dump_reads(
        "genotyping_overlap",
        &outcome.genotyping_reads,
        &geno_overlap,
    );
    dump_reads("covering_overlap", &covering.reads, &covering_overlap);
    let add_evidence_candidates: BTreeSet<usize> =
        geno_overlap.difference(&stored).copied().collect();
    kv(
        "addEvidence_candidate_n",
        add_evidence_candidates.len().to_string(),
    );
    dump_reads(
        "addEvidence_candidate",
        &outcome.genotyping_reads,
        &add_evidence_candidates,
    );

    let (full_ref, full_pad) = outcome.assembly.event_map_reference();
    let apply_pad = outcome
        .assembly
        .haplotypes
        .iter()
        .find(|h| h.is_reference)
        .and_then(|h| h.genome_loc.as_ref().map(|g| g.start.get()))
        .unwrap_or(full_pad);
    let hap_cache = build_per_haplotype_variation_events(
        &outcome.assembly.haplotypes,
        full_ref,
        full_pad,
        outcome.assembly.max_mnp_distance(),
        &covering.contig,
    );
    let mapping = create_allele_mapper_with_events(
        &VariationEvent::from_alleles(&covering.contig, TARGET, TARGET_REF, TARGET_ALT),
        TARGET,
        &outcome.assembly.haplotypes,
        apply_pad,
        outcome.assembly.reference_bases(),
        outcome.assembly.max_mnp_distance(),
        true,
        Some(&hap_cache),
    );
    kv(
        "allele_mapper_ref_haps",
        mapping.ref_haplotype_indices.len().to_string(),
    );
    kv(
        "allele_mapper_alt_haps",
        mapping.alt_haplotype_indices.len().to_string(),
    );
    let rows =
        region_likelihoods_to_rows(&outcome.read_likelihoods, outcome.assembly.haplotypes.len());
    let marg = marginalize_rows_to_biallelic_alleles(
        &rows,
        &mapping.ref_haplotype_indices,
        &mapping.alt_haplotype_indices,
    );
    let marg_idxs: BTreeSet<usize> = marg.iter().map(|r| r.read_index).collect();
    let retain: BTreeSet<usize> = marg_idxs
        .iter()
        .copied()
        .filter(|&idx| {
            outcome
                .genotyping_reads
                .get(idx)
                .is_some_and(|r| java_alignment_read_overlaps_interval(r, TARGET, TARGET, MARGIN))
        })
        .collect();
    kv("after_marginalize_n", marg_idxs.len().to_string());
    kv("after_retain_evidence_n", retain.len().to_string());
    assert_eq!(
        retain, stored_retain,
        "retainEvidence on the marginalized stored object is overlap of stored unique reads"
    );
    let diagnostic: Vec<RegionReadLikelihood> = outcome
        .read_likelihoods
        .iter()
        .filter(|c| retain.contains(&c.read_index.get()))
        .cloned()
        .collect();
    let diagnostic_n = unique_indices(&diagnostic);
    kv("diagnostic_java_reuse_n", diagnostic_n.len().to_string());
    assert_eq!(
        diagnostic_n, retain,
        "Java contamination-off reuse keeps retained stored-hap rows"
    );
    for idx in &retain {
        let rec = &outcome.genotyping_reads[*idx];
        let row = marg.iter().find(|r| r.read_index == *idx).expect("marg");
        let ll_ref = row
            .haplotype_log10_likelihoods
            .first()
            .copied()
            .unwrap_or(f64::NEG_INFINITY);
        let ll_alt = row
            .haplotype_log10_likelihoods
            .get(1)
            .copied()
            .unwrap_or(f64::NEG_INFINITY);
        let best = if ll_alt > ll_ref { "A" } else { "C" };
        kv(
            "retained_best_allele",
            format!(
                "idx={idx} qname={} best={best} C*={ll_ref:.6} A={ll_alt:.6}",
                qname(rec)
            ),
        );
    }

    let fallback = unique_indices(&outcome.read_likelihoods);
    let fallback_dp = coverage_evidence_count(
        &outcome.genotyping_reads,
        &outcome.read_likelihoods,
        TARGET,
        TARGET,
        MARGIN,
    );
    let diagnostic_dp = coverage_evidence_count(
        &outcome.genotyping_reads,
        &diagnostic,
        TARGET,
        TARGET,
        MARGIN,
    );
    kv("fallback_unique_n", fallback.len().to_string());
    kv("fallback_coverage_dp", fallback_dp.to_string());
    kv("diagnostic_coverage_dp", diagnostic_dp.to_string());
    kv("stored_unique_as_evidenceCount", stored.len().to_string());
    let stored_dropped: BTreeSet<usize> = stored.difference(&retain).copied().collect();
    kv(
        "stored_dropped_by_retain_n",
        stored_dropped.len().to_string(),
    );
    dump_reads(
        "stored_dropped_by_retain",
        &outcome.genotyping_reads,
        &stored_dropped,
    );
    let genotyping_not_in_stored: BTreeSet<usize> = (0..outcome.genotyping_reads.len())
        .filter(|i| !stored.contains(i))
        .collect();
    kv(
        "genotyping_not_in_stored_n",
        genotyping_not_in_stored.len().to_string(),
    );
    dump_reads(
        "genotyping_not_in_stored",
        &outcome.genotyping_reads,
        &genotyping_not_in_stored,
    );
    let stored_keys: BTreeSet<(String, u16)> = stored
        .iter()
        .filter_map(|i| {
            outcome
                .genotyping_reads
                .get(*i)
                .map(|r| (qname(r), r.flags()))
        })
        .collect();
    for rec in &covering.reads {
        let key = (qname(rec), rec.flags());
        if !stored_keys.contains(&key) {
            kv(
                "covering_not_in_stored",
                format!(
                    "qname={} flag={} mapq={} start={} end={} cigar={}",
                    key.0,
                    rec.flags(),
                    rec.mapq(),
                    rec.pos() + 1,
                    alignment_end_1based(rec),
                    rec.cigar()
                ),
            );
        }
    }
    assert_eq!(
        fallback_dp, 6,
        "6R.201 Coverage of stored hap unique evidence is n=6 (no second overlap filter)"
    );
    assert_eq!(
        diagnostic_dp, 4,
        "6R.186 helper (stored hap → marginalize → retainEvidence ±2) is n=4 at this SNP, not Java INFO DP=6"
    );
    assert_eq!(stored.len(), 6, "stored hap unique evidence is n=6");
    assert_eq!(retain.len(), 4);

    let (full_ref_b, full_pad_b) = outcome.assembly.event_map_reference();
    let apply_pad_b = outcome
        .assembly
        .haplotypes
        .iter()
        .find(|h| h.is_reference)
        .and_then(|h| h.genome_loc.as_ref().map(|g| g.start.get()))
        .unwrap_or(full_pad_b);
    let evidence_diag = HcStrandBiasLikelihoods {
        reads: &outcome.genotyping_reads,
        likelihoods: &diagnostic,
        haplotypes: &outcome.assembly.haplotypes,
        contig: &covering.contig,
        ref_bytes: outcome.assembly.reference_bases(),
        pad_start_1based: apply_pad_b,
        full_ref_bytes: full_ref_b,
        full_pad_1based: full_pad_b,
        max_mnp_distance: outcome.assembly.max_mnp_distance(),
        emit_spanning_dels: true,
    };
    let evidence_fb = HcStrandBiasLikelihoods {
        reads: &outcome.genotyping_reads,
        likelihoods: &outcome.read_likelihoods,
        haplotypes: &outcome.assembly.haplotypes,
        contig: &covering.contig,
        ref_bytes: outcome.assembly.reference_bases(),
        pad_start_1based: apply_pad_b,
        full_ref_bytes: full_ref_b,
        full_pad_1based: full_pad_b,
        max_mnp_distance: outcome.assembly.max_mnp_distance(),
        emit_spanning_dels: true,
    };
    let (rp_ref, rp_alt) =
        read_pos_rank_sum_quals_from_likelihoods(&evidence_diag, TARGET, TARGET_REF, TARGET_ALT);
    let (bq_ref, bq_alt) =
        baseq_rank_sum_quals_from_likelihoods(&evidence_diag, TARGET, TARGET_REF, TARGET_ALT);
    let (mq_ref, mq_alt) =
        mq_rank_sum_quals_from_likelihoods(&evidence_diag, TARGET, TARGET_REF, TARGET_ALT);
    let (f_rp_ref, f_rp_alt) =
        read_pos_rank_sum_quals_from_likelihoods(&evidence_fb, TARGET, TARGET_REF, TARGET_ALT);
    let (f_bq_ref, f_bq_alt) =
        baseq_rank_sum_quals_from_likelihoods(&evidence_fb, TARGET, TARGET_REF, TARGET_ALT);
    let (f_mq_ref, f_mq_alt) =
        mq_rank_sum_quals_from_likelihoods(&evidence_fb, TARGET, TARGET_REF, TARGET_ALT);
    kv(
        "diag_ranksum_lists",
        format!(
            "ReadPos REF={} ALT={} BaseQ REF={} ALT={} MQRS REF={} ALT={}",
            rp_ref.len(),
            rp_alt.len(),
            bq_ref.len(),
            bq_alt.len(),
            mq_ref.len(),
            mq_alt.len()
        ),
    );
    kv(
        "fallback_ranksum_lists",
        format!(
            "ReadPos REF={} ALT={} BaseQ REF={} ALT={} MQRS REF={} ALT={}",
            f_rp_ref.len(),
            f_rp_alt.len(),
            f_bq_ref.len(),
            f_bq_alt.len(),
            f_mq_ref.len(),
            f_mq_alt.len()
        ),
    );
    kv(
        "diag_both_ranksum_sides",
        format!(
            "ReadPos={} BaseQ={} MQRS={}",
            !rp_ref.is_empty() && !rp_alt.is_empty(),
            !bq_ref.is_empty() && !bq_alt.is_empty(),
            !mq_ref.is_empty() && !mq_alt.is_empty()
        ),
    );
    kv(
        "fallback_omits_ranksum_side",
        format!(
            "ReadPos_empty_side={} BaseQ_empty_side={}",
            f_rp_ref.is_empty() || f_rp_alt.is_empty(),
            f_bq_ref.is_empty() || f_bq_alt.is_empty()
        ),
    );

    let diag_mqs = rms_mapping_quality_sample_mapqs(&outcome.genotyping_reads, &diagnostic);
    let fb_mqs =
        rms_mapping_quality_sample_mapqs(&outcome.genotyping_reads, &outcome.read_likelihoods);
    let diag_rms = rms_mapping_quality_raw(&diag_mqs).map(|t| t.2);
    let fb_rms = rms_mapping_quality_raw(&fb_mqs).map(|t| t.2);
    kv("diag_mq_rms", format!("{diag_rms:?}"));
    kv("fallback_mq_rms", format!("{fb_rms:?}"));
    let fb_rms = fb_rms.expect("stored n=6 MQ RMS");
    let diag_rms = diag_rms.expect("retain n=4 MQ RMS");
    assert!(
        (fb_rms - 40.58).abs() < 0.01,
        "Java MQ 40.58 is RMS of stored unique n=6 MAPQs, not retain n=4"
    );
    assert!(
        (diag_rms - 44.15).abs() < 0.01,
        "6R.186 retain n=4 would change MQ to ~44.15"
    );

    let emitted =
        try_emit_call_region_variants(covering, &outcome, "NA12878", DEFAULT_STAND_EMIT_CONFIDENCE)
            .expect("emit");
    let rec = emitted
        .iter()
        .find(|r| r.position == TARGET && r.reference == TARGET_REF)
        .expect("C/A emit");
    assert_eq!(info_i32(&rec.info, "DP"), Some(6));
    let bq = info_f64(&rec.info, "BaseQRankSum").expect("6R.203 BaseQ");
    assert_eq!((bq * 1000.0).round(), -1834.0);
    let rprs = info_f64(&rec.info, "ReadPosRankSum").expect("6R.203 ReadPos");
    assert_eq!((rprs * 1000.0).round(), 1282.0);
    let mqrs = info_f64(&rec.info, "MQRankSum").expect("MQRankSum fallback still emits");
    assert_eq!((mqrs * 1000.0).round(), 1834.0);
    let sample = rec.samples.first().expect("sample");
    assert_eq!(
        sample.gt.as_ref().map(|g| g.to_string()).as_deref(),
        Some("0/1")
    );
    assert_eq!(sample.ad.as_deref(), Some(&[2u32, 4][..]));
    assert_eq!(sample.pl.as_deref(), Some(&[162u32, 0, 72][..]));
    assert_eq!(target.genotype.format.pl_as_i32(), vec![162, 0, 72]);
    kv(
        "first_arrow",
        "SiteEarlyTemplate is_cluster_downstream_snp returns shaped FORMAT without annotation_likelihoods_from_stored_haplotypes; Java loc-loop always constructs marg→retain→reuse→addEvidence. Helper retain at this SNP is n=4, not Java evidenceCount 6.",
    );
    kv(
        "recommendation_6r199",
        "construct the loc-loop annotation object on the cluster-downstream early-template. Do not copy 6R.186 retain as-is (n=4, MQ would become 44.15). Java evidenceCount is stored unique n=6. Do not retarget the empty fallback.",
    );
}
