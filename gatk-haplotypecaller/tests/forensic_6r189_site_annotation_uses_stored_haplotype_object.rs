//! 6R.189: SiteScore annotation binds to the Java loc-loop stored-hap object
//! (`marginalize` then `retainEvidence(±2)`), not the FORMAT genotyping subset.
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! Target `2:92305635 A/G`. FORMAT subset is unchanged.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r189_site_annotation_uses_stored_haplotype_object -- --nocapture --test-threads=1
//! HOLDOUT_6R189=1 cargo test -p gatk-haplotypecaller --test holdout_6r189_site_annotation_stored_hap -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::{build_per_haplotype_variation_events, VariationEvent};
use gatk_haplotypecaller::hc_allele_mapping::create_allele_mapper_with_events;
use gatk_haplotypecaller::hc_genotyping_engine::java_alignment_read_overlaps_interval;
use gatk_haplotypecaller::read_realignment::LOG_10_INFORMATIVE_THRESHOLD;
use gatk_haplotypecaller::variant_site_hc_annotations::coverage_evidence_count;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    HaplotypeCallerEngine, ReadFilterParams, RegionReadLikelihood, WalkerTraversalConfig,
    DEFAULT_STAND_EMIT_CONFIDENCE,
};
use rust_htslib::bam::Record;
use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::path::{Path, PathBuf};

const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
const INTERVAL: &str = "2:92305500-92305850";
const TG_INTERVAL: &str = "2:92307200-92307550";
const MID_INTERVAL: &str = "2:92316200-92316580";
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 92_305_635;
const MERGED_REF: &str = "A";
const MERGED_ALT: &str = "G";
const CLOSED_SNP: u64 = 92_305_634;
const CLOSED_INDEL: u64 = 92_307_324;
const CLOSED_TG: u64 = 92_307_333;
const CLOSED_MID: u64 = 92_316_347;
const MARGIN: i32 = 2;
const MAPQ25_QNAME: &str = "H06JUADXX130110:1:1101:10018:4569";
const MATE_QNAME: &str = "H06HDADXX130110:2:1101:10029:53752";

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R189\t{key}\t{}", value.as_ref());
}

fn unique_indices(likelihoods: &[RegionReadLikelihood]) -> BTreeSet<usize> {
    likelihoods.iter().map(|c| c.read_index.get()).collect()
}

fn qname(rec: &Record) -> String {
    String::from_utf8_lossy(rec.qname()).into_owned()
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

/// Copy of production `per_variant_annotation_likelihoods` (FORMAT-subset reconstruction).
fn hypothetical_annotation_from_subset(
    subset: &[RegionReadLikelihood],
    reads: &[gatk_haplotypecaller::SharedBamRecord],
) -> Vec<RegionReadLikelihood> {
    if subset.is_empty() {
        return Vec::new();
    }
    let mut best_ll: BTreeMap<usize, f64> = BTreeMap::new();
    for cell in subset {
        let idx = cell.read_index.get();
        let ll = cell.log10_likelihood;
        best_ll
            .entry(idx)
            .and_modify(|m| {
                if ll > *m {
                    *m = ll;
                }
            })
            .or_insert(ll);
    }
    let mut qnames: BTreeSet<Vec<u8>> = BTreeSet::new();
    for &idx in best_ll.keys() {
        if let Some(rec) = reads.get(idx) {
            qnames.insert(rec.qname().to_owned());
        }
    }
    if qnames.len() == 1 && best_ll.len() > 1 {
        let Some((&keep_idx, _)) = best_ll.iter().max_by(|a, b| a.1.total_cmp(b.1)) else {
            return subset.to_vec();
        };
        return subset
            .iter()
            .filter(|cell| cell.read_index.get() == keep_idx)
            .cloned()
            .collect();
    }
    subset.to_vec()
}

#[test]
fn forensic_6r189_source_contracts() {
    kv("java_pin", JAVA_PIN);
    kv(
        "production_change",
        "SiteScore annotation_likelihoods from stored haplotypes (marginalize then retainEvidence)",
    );
    kv("classification", "A2 — WRONG SUBSET (closed)");
    kv("target", "2:92305635 A/G");

    let early = include_str!("../src/hc_genotyping_engine/genotype_site_early_template.rs");
    let pipe = include_str!("../src/hc_genotyping_engine/genotype_site_pipeline.rs");
    let ge = include_str!("../src/hc_genotyping_engine/mod.rs");
    let emit = include_str!("../src/region_vcf_emit.rs");
    let ann = include_str!("../src/variant_site_hc_annotations.rs");
    let engine = include_str!("../src/engine.rs");

    let helper = early
        .split("fn annotation_likelihoods_from_stored_haplotypes")
        .nth(1)
        .expect("helper");
    assert!(
        helper.contains("marginalize_rows_to_biallelic_alleles")
            && helper.contains("likelihood_subset_for_event"),
        "helper is stored hap → marginalize → retainEvidence"
    );
    assert!(
        !helper.contains("92305635")
            && !helper.contains("92305634")
            && !helper.contains(MAPQ25_QNAME)
            && !helper.contains(MATE_QNAME),
        "helper has no locus/QNAME/FLAG special case"
    );

    let sitescore = pipe
        .split("// L9: PairHMM genotype failed Java emit")
        .next()
        .expect("L9");
    let sitescore_attach = sitescore
        .split("let annotation = annotation_likelihoods_from_stored_haplotypes")
        .nth(1)
        .expect("SiteScore stored-hap attach");
    assert!(
        sitescore_attach.contains("likelihoods,")
            && sitescore.contains("with_annotation_likelihoods(annotation)"),
        "SiteScore first argument is the stored hap matrix, not FORMAT subset"
    );
    assert!(
        !sitescore_attach.contains("per_variant_annotation_likelihoods(&subset"),
        "SiteScore does not attach the FORMAT genotyping subset"
    );
    assert!(
        pipe.contains("let keep_qnames:") && pipe.contains("subset = narrowed"),
        "FORMAT keep_qnames subset construction is unchanged"
    );
    assert!(
        ge.contains("fn per_variant_annotation_likelihoods"),
        "6R.180 helper remains for FORMAT-subset collapse where still used"
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
        "6R.186 cluster-TG attachment stays"
    );
    assert!(
        early.contains("if is_p12_phase_e_two_read_hom_alt_site(&event)")
            && early.contains("finish_strict_java_shaped_site_call"),
        "neighbor G/T early-template path stays"
    );
    assert!(
        emit.contains("if call.annotation_likelihoods.is_empty()")
            && emit.contains("likelihoods: &call.annotation_likelihoods"),
        "emit still consumes the attached object when present"
    );
    assert!(
        ann.contains("6R.174: INFO DP is Java `Coverage.evidenceCount`")
            && ann.contains("mq_rms_of_sample_evidence")
            && ann.contains("strand_bias_sample_counts"),
        "DP/MQ/SOR formulas must stay closed"
    );
    assert!(
        engine.contains("fn apply_java_order_normalize_and_filter"),
        "6R.185 stored-matrix lifecycle stays closed"
    );
    assert!(
        !ann.contains("92305635")
            && !emit.contains("92305635")
            && !pipe.contains("92305635")
            && !engine.contains("92305635"),
        "no locus-specific 6R.189 production patch"
    );
}

#[test]
fn forensic_6r189_site_annotation_uses_stored_haplotype_object() {
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing P12 ref/BAM");
        return;
    }
    kv("java_pin", JAVA_PIN);
    kv("target", "2:92305635 A/G");

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
    let outcome = HaplotypeCallerEngine::call_region(
        covering,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("outcome");

    let stored = unique_indices(&outcome.read_likelihoods);
    kv("stored_hap_n", stored.len().to_string());
    assert_eq!(stored.len(), 3, "stored hap matrix is Java n=3");

    let closed = outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(CLOSED_SNP)
                && c.event.ref_allele == "G"
                && c.event.alt_allele == "T"
        })
        .expect("closed G/T");
    let target = outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(TARGET)
                && c.event.ref_allele == MERGED_REF
                && c.event.alt_allele == MERGED_ALT
        })
        .expect("target A/G");
    assert_eq!(
        unique_indices(&closed.annotation_likelihoods).len(),
        3,
        "neighbor 92305634 two-read arm attaches retainEvidence n=3 (6R.205)"
    );

    let attached = unique_indices(&target.annotation_likelihoods);
    kv("attached_n", attached.len().to_string());
    kv(
        "attached_cells",
        target.annotation_likelihoods.len().to_string(),
    );
    let hap_cols: BTreeSet<usize> = target
        .annotation_likelihoods
        .iter()
        .map(|c| c.haplotype_index.get())
        .collect();
    kv(
        "attached_hap_columns",
        format!("{:?}", hap_cols.iter().collect::<Vec<_>>()),
    );
    assert_eq!(attached.len(), 3, "annotation object unique n=3");
    assert_eq!(
        target.genotype.format.ad_as_i32(),
        vec![0, 2],
        "FORMAT AD unchanged"
    );
    assert_eq!(target.genotype.format.dp.as_i32(), 2, "FORMAT DP unchanged");
    assert_eq!(
        target.genotype.format.pl_as_i32(),
        vec![90, 6, 0],
        "FORMAT PL unchanged"
    );
    assert_eq!(target.genotype.format.gq.as_i32(), 6, "FORMAT GQ unchanged");

    let mut overlap: BTreeSet<usize> = BTreeSet::new();
    for idx in &stored {
        let rec = &outcome.genotyping_reads[*idx];
        if java_alignment_read_overlaps_interval(rec, TARGET, TARGET, MARGIN) {
            overlap.insert(*idx);
        }
    }
    kv("java_retainEvidence_n", overlap.len().to_string());
    assert_eq!(overlap.len(), 3);
    assert_eq!(
        attached, overlap,
        "annotation membership is stored-hap retainEvidence(±2), not FORMAT subset"
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
        &VariationEvent::from_alleles(&covering.contig, TARGET, MERGED_REF, MERGED_ALT),
        TARGET,
        &outcome.assembly.haplotypes,
        apply_pad,
        outcome.assembly.reference_bases(),
        outcome.assembly.max_mnp_distance(),
        true,
        Some(&hap_cache),
    );
    let rows = gatk_haplotypecaller::region_likelihoods_to_rows(
        &outcome.read_likelihoods,
        outcome.assembly.haplotypes.len(),
    );
    let marg = gatk_haplotypecaller::marginalize_rows_to_biallelic_alleles(
        &rows,
        &mapping.ref_haplotype_indices,
        &mapping.alt_haplotype_indices,
    );
    let mut alt_favoring = Vec::new();
    for row in &marg {
        if !overlap.contains(&row.read_index) {
            continue;
        }
        let lr = row
            .haplotype_log10_likelihoods
            .first()
            .copied()
            .unwrap_or(f64::NEG_INFINITY);
        let la = row
            .haplotype_log10_likelihoods
            .get(1)
            .copied()
            .unwrap_or(f64::NEG_INFINITY);
        if la > lr && (la - lr) > LOG_10_INFORMATIVE_THRESHOLD {
            alt_favoring.push((row.read_index, lr));
        }
    }
    alt_favoring.sort_by(|a, b| a.1.total_cmp(&b.1));
    let keep_n = 2.min(alt_favoring.len()).max(1);
    let keep_qnames: HashSet<Vec<u8>> = alt_favoring
        .iter()
        .take(keep_n)
        .filter_map(|(idx, _)| {
            outcome
                .genotyping_reads
                .get(*idx)
                .map(|r| r.qname().to_owned())
        })
        .collect();
    let format_subset: Vec<RegionReadLikelihood> = outcome
        .read_likelihoods
        .iter()
        .filter(|c| {
            outcome
                .genotyping_reads
                .get(c.read_index.get())
                .is_some_and(|r| keep_qnames.contains(r.qname()))
        })
        .cloned()
        .collect();
    let format_unique = unique_indices(&format_subset);
    let after_collapse = unique_indices(&hypothetical_annotation_from_subset(
        &format_subset,
        &outcome.genotyping_reads,
    ));
    kv("format_subset_n", format_unique.len().to_string());
    kv(
        "per_variant_of_format_subset_n",
        after_collapse.len().to_string(),
    );
    assert_eq!(format_unique.len(), 2, "FORMAT keep_qnames still n=2 mates");
    assert_eq!(
        after_collapse.len(),
        1,
        "FORMAT/per-variant collapse still n=1"
    );
    assert_ne!(
        after_collapse, attached,
        "annotation object is not the FORMAT genotyping subset"
    );

    let mapq25 = stored.iter().copied().find(|&idx| {
        qname(&outcome.genotyping_reads[idx]) == MAPQ25_QNAME
            && outcome.genotyping_reads[idx].mapq() == 25
    });
    let mate_rev = stored.iter().copied().find(|&idx| {
        qname(&outcome.genotyping_reads[idx]) == MATE_QNAME
            && outcome.genotyping_reads[idx].flags() == 147
    });
    let mate_fwd = stored.iter().copied().find(|&idx| {
        qname(&outcome.genotyping_reads[idx]) == MATE_QNAME
            && outcome.genotyping_reads[idx].flags() == 99
    });
    let mapq25 = mapq25.expect("MAPQ=25");
    let mate_rev = mate_rev.expect("FLAG=147");
    let mate_fwd = mate_fwd.expect("FLAG=99");
    for idx in [mapq25, mate_fwd, mate_rev] {
        let rec = &outcome.genotyping_reads[idx];
        kv(
            "READ",
            format!(
                "QNAME={} FLAG={} MAPQ={} start={} end={} stored=true java_ann=true rust_ann={}",
                qname(rec),
                rec.flags(),
                rec.mapq(),
                rec.pos() + 1,
                gatk_haplotypecaller::read_unclip::alignment_end_1based(rec),
                attached.contains(&idx)
            ),
        );
        assert!(attached.contains(&idx));
        assert!(!after_collapse.contains(&idx) || idx == mate_fwd);
    }
    assert!(
        after_collapse.contains(&mate_fwd) && !after_collapse.contains(&mapq25),
        "FORMAT subset still drops MAPQ=25"
    );

    let rust_dp = coverage_evidence_count(
        &outcome.genotyping_reads,
        &target.annotation_likelihoods,
        TARGET,
        TARGET,
        MARGIN,
    );
    kv("coverage_on_attached", rust_dp.to_string());
    assert_eq!(rust_dp, 3);

    let emitted =
        try_emit_call_region_variants(covering, &outcome, "NA12878", DEFAULT_STAND_EMIT_CONFIDENCE)
            .expect("emit");
    let closed_rec = emitted
        .iter()
        .find(|r| r.position == CLOSED_SNP && r.reference == "G")
        .expect("closed emit");
    let rec = emitted
        .iter()
        .find(|r| r.position == TARGET && r.reference == MERGED_REF)
        .expect("target emit");
    let sample = rec.samples.first().expect("sample");
    assert_eq!(
        sample.gt.as_ref().map(|g| g.to_string()).as_deref(),
        Some("1/1")
    );
    assert_eq!(sample.ad.as_deref(), Some(&[0u32, 2][..]));
    assert_eq!(sample.dp.map(|d| d as i32), Some(2));
    assert_eq!(sample.gq.map(|g| g as i32), Some(6));
    assert_eq!(sample.pl.as_deref(), Some(&[90u32, 6, 0][..]));
    assert!((rec.quality.unwrap_or(0.0) - 78.32).abs() < 0.05);
    assert_eq!(info_i32(&closed_rec.info, "DP"), Some(3));
    assert_eq!(info_i32(&rec.info, "DP"), Some(3));
    let mq = info_f64(&rec.info, "MQ").unwrap_or(-1.0);
    let sor = info_f64(&rec.info, "SOR").unwrap_or(-1.0);
    let closed_mq = info_f64(&closed_rec.info, "MQ").unwrap_or(-1.0);
    let closed_sor = info_f64(&closed_rec.info, "SOR").unwrap_or(-1.0);
    assert!((mq - 41.96).abs() < 0.005, "MQ={mq}");
    assert!((sor - 0.693).abs() < 0.002, "SOR={sor}");
    assert!((closed_mq - 41.96).abs() < 0.005);
    assert!((closed_sor - 0.693).abs() < 0.002);
    kv("format", "GT=1/1 AD=0,2 DP=2 GQ=6 PL=90,6,0 QUAL=78.32");
    kv("info", format!("DP=3 MQ={mq:.2} SOR={sor:.3}"));

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
    let tg = tg_emitted
        .iter()
        .find(|r| r.position == CLOSED_TG && r.reference == "T")
        .expect("T/G");
    assert_eq!(info_i32(&indel.info, "DP"), Some(1));
    assert_eq!(info_i32(&tg.info, "DP"), Some(1));
    let indel_mq = info_f64(&indel.info, "MQ").unwrap_or(-1.0);
    let tg_mq = info_f64(&tg.info, "MQ").unwrap_or(-1.0);
    assert!((indel_mq - 44.0).abs() < 0.005);
    assert!((tg_mq - 44.0).abs() < 0.005);
    let tg_call = tg_outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(CLOSED_TG)
                && c.event.ref_allele == "T"
                && c.event.alt_allele == "G"
        })
        .expect("T/G call");
    let tg_ann = unique_indices(&tg_call.annotation_likelihoods);
    assert_eq!(tg_ann.len(), 1, "6R.186 n=1 stays");

    let mid_specs = parse_intervals_cli_string(&dict, MID_INTERVAL).expect("mid");
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
    let mid_covering = mid_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_MID
                && r.end.get() >= CLOSED_MID
        })
        .expect("mid covering");
    let mid_outcome = HaplotypeCallerEngine::call_region(
        mid_covering,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("mid call")
    .expect("mid outcome");
    let mid_emitted = try_emit_call_region_variants(
        mid_covering,
        &mid_outcome,
        "NA12878",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .expect("mid emit");
    let mid_rec = mid_emitted
        .iter()
        .find(|r| r.position == CLOSED_MID && r.reference == "G")
        .expect("2:92316347 G/A");
    let mid_mq = info_f64(&mid_rec.info, "MQ").unwrap_or(-1.0);
    let mid_sor = info_f64(&mid_rec.info, "SOR").unwrap_or(-1.0);
    assert!((mid_mq - 40.25).abs() < 0.005, "MQ={mid_mq}");
    assert!((mid_sor - 1.179).abs() < 0.002, "SOR={mid_sor}");
}
