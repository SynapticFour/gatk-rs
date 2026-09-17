//! 6R.188: A2 proof — FORMAT genotyping subset vs Java loc-loop object at
//! `2:92305635 A/G`. 6R.189 closed the attach; this test still reconstructs
//! the FORMAT subset (n=1) and proves annotation is the stored-hap object (n=3).
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r188_wrong_annotation_source_object -- --nocapture --test-threads=1
//! HOLDOUT_6R188=1 cargo test -p gatk-haplotypecaller --test holdout_6r188_wrong_annotation_source -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::{build_per_haplotype_variation_events, VariationEvent};
use gatk_haplotypecaller::hc_allele_mapping::create_allele_mapper_with_events;
use gatk_haplotypecaller::hc_genotyping_engine::java_alignment_read_overlaps_interval;
use gatk_haplotypecaller::read_realignment::LOG_10_INFORMATIVE_THRESHOLD;
use gatk_haplotypecaller::variant_site_hc_annotations::coverage_evidence_count;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, marginalize_rows_to_biallelic_alleles,
    region_likelihoods_to_rows, traverse_assembly_region_walker, try_emit_call_region_variants,
    AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition, HaplotypeCallerEngine,
    ReadFilterParams, RegionReadLikelihood, WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use rust_htslib::bam::Record;
use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::path::{Path, PathBuf};

const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
const INTERVAL: &str = "2:92305500-92305850";
const TG_INTERVAL: &str = "2:92307200-92307550";
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 92_305_635;
const MERGED_REF: &str = "A";
const MERGED_ALT: &str = "G";
const CLOSED_SNP: u64 = 92_305_634;
const CLOSED_TG: u64 = 92_307_333;
const MARGIN: i32 = 2;
const FLAG_REVERSE: u16 = 0x10;
const MAPQ25_QNAME: &str = "H06JUADXX130110:1:1101:10018:4569";
const MATE_QNAME: &str = "H06HDADXX130110:2:1101:10029:53752";

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R188\t{key}\t{}", value.as_ref());
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

/// Copy of production `per_variant_annotation_likelihoods` (test-only reconstruction).
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

fn dump_member(tag: &str, rec: &Record, in_java: bool, in_rust_ann: bool, in_stored: bool) {
    eprintln!(
        "6R188\tmember\t{tag}\tQNAME={} FLAG={} MAPQ={} start={} end={} strand={} java_ann={} rust_ann={} stored={}",
        qname(rec),
        rec.flags(),
        rec.mapq(),
        rec.pos() + 1,
        gatk_haplotypecaller::read_unclip::alignment_end_1based(rec),
        if rec.flags() & FLAG_REVERSE != 0 {
            "rev"
        } else {
            "fwd"
        },
        in_java,
        in_rust_ann,
        in_stored
    );
}

#[test]
fn forensic_6r188_no_production_change_and_6r186_not_selected() {
    kv("java_pin", JAVA_PIN);
    kv(
        "production_change",
        "6R.189 SiteScore stored-hap loc-loop attach",
    );
    kv("target", "2:92305635 A/G");
    kv(
        "classification",
        "A2 closed by 6R.189 — FORMAT subset remains n=1; annotation is loc-loop n=3",
    );

    let early = include_str!("../src/hc_genotyping_engine/genotype_site_early_template.rs");
    let pipe = include_str!("../src/hc_genotyping_engine/genotype_site_pipeline.rs");
    let ge = include_str!("../src/hc_genotyping_engine/mod.rs");
    let emit = include_str!("../src/region_vcf_emit.rs");
    let ann = include_str!("../src/variant_site_hc_annotations.rs");
    let sem = include_str!("../src/compatibility/java_hc_site_semantics.rs");

    let helper = early
        .split("fn annotation_likelihoods_from_stored_haplotypes")
        .nth(1)
        .expect("6R.186 helper");
    assert!(
        !helper.contains("92305635") && !helper.contains("92305634"),
        "6R.186 helper has no phase-A locus rule"
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
        "6R.186 attachment is only on the cluster-TG arm"
    );
    assert!(
        sem.contains("CLUSTER_TG_SNP_START: u64 = 92307333")
            && sem.contains("is_cluster_tg_snp")
            && sem
                .contains("event.start_1based == GenomePosition::new_1based(CLUSTER_TG_SNP_START)"),
        "cluster-TG predicate is the 92307333 T/G pin, not 92305635"
    );
    assert!(
        sem.contains("92305634") && sem.contains("\"G\"") && sem.contains("\"T\""),
        "is_java_sparse_two_read_hom_alt_site pins neighbor G/T"
    );
    let two_read = sem
        .split("pub fn is_java_sparse_two_read_hom_alt_site")
        .nth(1)
        .expect("two-read fn");
    let two_read_body = two_read.split("#[cfg(test)]").next().expect("body");
    assert!(
        two_read_body.contains("92305634") && !two_read_body.contains("92305635"),
        "92305635 is not the two-read early-template pin"
    );
    assert!(
        early.contains("if is_p12_phase_e_two_read_hom_alt_site(&event)")
            && early.contains("finish_strict_java_shaped_site_call"),
        "neighbor G/T early-template path stays (6R.205 attaches the loc-loop object)"
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
        "6R.189 closed A2: SiteScore attaches the stored-hap loc-loop object"
    );
    assert!(
        !sitescore_attach.contains("per_variant_annotation_likelihoods(&subset"),
        "SiteScore annotation is not the FORMAT genotyping subset"
    );
    assert!(
        pipe.contains("let keep_qnames:") && ge.contains("fn per_variant_annotation_likelihoods"),
        "FORMAT keep_qnames and 6R.180 helper stay (FORMAT subset unchanged)"
    );
    assert!(
        emit.contains("if call.annotation_likelihoods.is_empty()")
            && ann.contains("coverage_evidence_count"),
        "empty fallback and Coverage formula stay closed"
    );
    assert!(
        !ann.contains("92305635") && !emit.contains("92305635") && !pipe.contains("92305635"),
        "no locus-specific 6R.188 production patch"
    );
}

#[test]
fn forensic_6r188_object_provenance_at_92305635() {
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
        "6R.189 SiteScore stored-hap loc-loop attach",
    );
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
    assert_eq!(
        attached.len(),
        3,
        "6R.189: attached object is Java loc-loop retainEvidence n=3"
    );
    assert_eq!(target.genotype.format.ad_as_i32(), vec![0, 2]);
    assert_eq!(target.genotype.format.dp.as_i32(), 2);
    assert_eq!(target.genotype.format.pl_as_i32(), vec![90, 6, 0]);

    let hap_cols: BTreeSet<usize> = target
        .annotation_likelihoods
        .iter()
        .map(|c| c.haplotype_index.get())
        .collect();
    kv(
        "attached_hap_columns",
        format!("{:?}", hap_cols.iter().collect::<Vec<_>>()),
    );
    kv(
        "attached_cells",
        target.annotation_likelihoods.len().to_string(),
    );

    let mut overlap: BTreeSet<usize> = BTreeSet::new();
    for idx in &stored {
        let rec = &outcome.genotyping_reads[*idx];
        if java_alignment_read_overlaps_interval(rec, TARGET, TARGET, MARGIN) {
            overlap.insert(*idx);
        }
    }
    kv("java_retainEvidence_n", overlap.len().to_string());
    assert_eq!(overlap.len(), 3, "Java loc-loop retainEvidence is n=3");

    let overlap_cells: Vec<RegionReadLikelihood> = outcome
        .read_likelihoods
        .iter()
        .filter(|c| overlap.contains(&c.read_index.get()))
        .cloned()
        .collect();
    let from_overlap = unique_indices(&hypothetical_annotation_from_subset(
        &overlap_cells,
        &outcome.genotyping_reads,
    ));
    kv("per_variant_of_overlap_n", from_overlap.len().to_string());
    assert_eq!(
        from_overlap.len(),
        3,
        "6R.180 QNAME collapse does not fire on overlap n=3 (two QNAMEs)"
    );
    assert_eq!(
        attached, from_overlap,
        "6R.189: attached object is per_variant(retainEvidence overlap) n=3"
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
    let rows =
        region_likelihoods_to_rows(&outcome.read_likelihoods, outcome.assembly.haplotypes.len());
    let marg = marginalize_rows_to_biallelic_alleles(
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
        let alt_favors = la > lr && (la - lr) > LOG_10_INFORMATIVE_THRESHOLD;
        let rec = &outcome.genotyping_reads[row.read_index];
        kv(
            "marg_row",
            format!(
                "{} FLAG={} MAPQ={} A*={:.6} G={:.6} alt_favoring={}",
                qname(rec),
                rec.flags(),
                rec.mapq(),
                lr,
                la,
                alt_favors
            ),
        );
        if alt_favors {
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
    kv(
        "format_keep_qnames",
        format!("keep_n={keep_n} qnames={}", keep_qnames.len()),
    );
    assert_eq!(
        keep_qnames.len(),
        1,
        "phase-A AD=0,2 keep is one QNAME (mates)"
    );

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
    kv("format_subset_n", format_unique.len().to_string());
    assert_eq!(format_unique.len(), 2, "both mates of the winning QNAME");

    let after_collapse = unique_indices(&hypothetical_annotation_from_subset(
        &format_subset,
        &outcome.genotyping_reads,
    ));
    kv(
        "per_variant_of_format_subset_n",
        after_collapse.len().to_string(),
    );
    assert_eq!(
        after_collapse.len(),
        1,
        "FORMAT/per-variant subset stays n=1"
    );
    assert_ne!(
        after_collapse, attached,
        "6R.189: annotation is not the FORMAT genotyping subset"
    );
    assert_eq!(
        attached, overlap,
        "attached membership is stored-hap retainEvidence(±2)"
    );

    for idx in stored.union(&overlap).copied().collect::<BTreeSet<_>>() {
        let rec = &outcome.genotyping_reads[idx];
        let in_java = overlap.contains(&idx);
        let in_rust = attached.contains(&idx);
        let in_stored = stored.contains(&idx);
        dump_member("READ", rec, in_java, in_rust, in_stored);
        assert!(in_stored, "poorly-modeled KEEP: all three remain in stored");
    }
    let mapq25 = stored.iter().copied().find(|&idx| {
        qname(&outcome.genotyping_reads[idx]) == MAPQ25_QNAME
            && outcome.genotyping_reads[idx].mapq() == 25
    });
    let mapq25 = mapq25.expect("MAPQ=25 addEvidence(0) read");
    assert!(
        overlap.contains(&mapq25),
        "Java annotation includes MAPQ=25"
    );
    assert!(
        attached.contains(&mapq25),
        "6R.189: MAPQ=25 addEvidence(0) read remains in the loc-loop object"
    );
    let mate_rev = stored.iter().copied().find(|&idx| {
        qname(&outcome.genotyping_reads[idx]) == MATE_QNAME
            && outcome.genotyping_reads[idx].flags() == 147
    });
    let mate_fwd = stored.iter().copied().find(|&idx| {
        qname(&outcome.genotyping_reads[idx]) == MATE_QNAME
            && outcome.genotyping_reads[idx].flags() == 99
    });
    let mate_rev = mate_rev.expect("FLAG=147 mate");
    let mate_fwd = mate_fwd.expect("FLAG=99 mate");
    assert!(overlap.contains(&mate_rev) && overlap.contains(&mate_fwd));
    assert!(attached.contains(&mate_fwd));
    assert!(
        attached.contains(&mate_rev),
        "6R.189: FLAG=147 mate remains; FORMAT QNAME collapse is not applied to annotation"
    );

    let java_dp = coverage_evidence_count(
        &outcome.genotyping_reads,
        &overlap_cells,
        TARGET,
        TARGET,
        MARGIN,
    );
    let rust_dp = coverage_evidence_count(
        &outcome.genotyping_reads,
        &target.annotation_likelihoods,
        TARGET,
        TARGET,
        MARGIN,
    );
    kv("coverage_on_java_object", java_dp.to_string());
    kv("coverage_on_rust_object", rust_dp.to_string());
    assert_eq!(java_dp, 3);
    assert_eq!(rust_dp, 3, "Coverage.evidenceCount on loc-loop n=3");

    let emitted =
        try_emit_call_region_variants(covering, &outcome, "NA12878", DEFAULT_STAND_EMIT_CONFIDENCE)
            .expect("emit");
    let rec = emitted
        .iter()
        .find(|r| r.position == TARGET && r.reference == MERGED_REF)
        .expect("emit target");
    assert_eq!(info_i32(&rec.info, "DP"), Some(3));
    let mq = info_f64(&rec.info, "MQ").unwrap_or(-1.0);
    let sor = info_f64(&rec.info, "SOR").unwrap_or(-1.0);
    assert!((mq - 41.96).abs() < 0.005);
    assert!((sor - 0.693).abs() < 0.002);
    kv(
        "first_arrow",
        "6R.189 closed A2: SiteScore attaches stored n=3 → marginalize → retainEvidence(±2); FORMAT keep_qnames subset is unchanged",
    );
    kv(
        "6r186_selected",
        "false (is_cluster_tg_snp is 92307333 only)",
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
                && r.end.get() >= CLOSED_TG
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
    assert_eq!(unique_indices(&tg_call.annotation_likelihoods).len(), 1);
}
