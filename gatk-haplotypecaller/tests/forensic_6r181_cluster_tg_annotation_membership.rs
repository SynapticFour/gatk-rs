//! 6R.181: why `2:92307333 T/G` has no Java-equivalent per-variant annotation
//! `AlleleLikelihoods` (proof-only).
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! PRODUCTION CHANGE: NONE.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r181_cluster_tg_annotation_membership -- --nocapture --test-threads=1
//! HOLDOUT_6R181=1 cargo test -p gatk-haplotypecaller --test holdout_6r181_cluster_tg -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::{build_per_haplotype_variation_events, VariationEvent};
use gatk_haplotypecaller::hc_allele_mapping::create_allele_mapper_with_events;
use gatk_haplotypecaller::hc_genotyping_engine::java_alignment_read_overlaps_interval;
use gatk_haplotypecaller::read_assembly_filter::{passes_assembly_read, AssemblyReadFilterConfig};
use gatk_haplotypecaller::read_model::MAPPING_QUALITY_UNAVAILABLE;
use gatk_haplotypecaller::read_realignment::LOG_10_INFORMATIVE_THRESHOLD;
use gatk_haplotypecaller::variant_site_hc_annotations::{
    coverage_evidence_count, rms_mapping_quality_raw, rms_mapping_quality_sample_mapqs,
    strand_bias_contingency_table, HcStrandBiasLikelihoods, STRAND_ODDS_RATIO_MIN_COUNT,
};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, marginalize_rows_to_biallelic_alleles,
    region_likelihoods_to_rows, traverse_assembly_region_walker, try_emit_call_region_variants,
    AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition, HaplotypeCallerEngine,
    ReadFilterParams, RegionReadLikelihood, SharedBamRecord, WalkerTraversalConfig,
    DEFAULT_STAND_EMIT_CONFIDENCE,
};
use rust_htslib::bam::Record;
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

const INTERVAL: &str = "2:92307200-92307550";
const CLOSED_SNP: u64 = 92_305_634;
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const JAVA_VCF_REL: &str = "parity/reports/6r43/p12_indel_mix/java.vcf";
const TARGET: u64 = 92_307_333;
const MERGED_REF: &str = "T";
const MERGED_ALT: &str = "G";
const CLOSED_INDEL: u64 = 92_307_324;
const CLOSED_INDEL_REF: &str = "TTC";
const CLOSED_INDEL_ALT: &str = "T";
const FLAG_REVERSE: u16 = 0x10;
const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
const MARGIN: i32 = 2;
const JAVA_MQ44_QNAME: &str = "H06JUADXX130110:1:1101:10052:88682";

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R181\t{key}\t{}", value.as_ref());
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

fn qname(rec: &Record) -> String {
    String::from_utf8_lossy(rec.qname()).into_owned()
}

fn strand_label(rec: &Record) -> &'static str {
    if rec.flags() & FLAG_REVERSE != 0 {
        "rev"
    } else {
        "fwd"
    }
}

fn parse_info_map(info: &str) -> BTreeMap<String, String> {
    let mut out = BTreeMap::new();
    for piece in info.split(';') {
        if piece.is_empty() {
            continue;
        }
        match piece.split_once('=') {
            Some((k, v)) => {
                out.insert(k.to_string(), v.to_string());
            }
            None => {
                out.insert(piece.to_string(), String::new());
            }
        }
    }
    out
}

fn java_target_record(path: &Path) -> (BTreeMap<String, String>, String, String, String) {
    let text = std::fs::read_to_string(path).expect("java.vcf");
    for line in text.lines() {
        if line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 10 {
            continue;
        }
        if f[0] == "2" && f[1] == "92307333" && f[3] == MERGED_REF && f[4] == MERGED_ALT {
            return (
                parse_info_map(f[7]),
                f[5].to_string(),
                f[8].to_string(),
                f[9].to_string(),
            );
        }
    }
    panic!("java.vcf missing 2:92307333 T/G");
}

fn mate_ok_java(rec: &Record) -> bool {
    if !rec.is_paired() || rec.is_mate_unmapped() || rec.is_unmapped() {
        return true;
    }
    rec.tid() == rec.mtid()
}

fn dump_read(prefix: &str, rec: &Record) {
    let cfg = AssemblyReadFilterConfig::gatk_defaults();
    eprintln!(
        "6R181\t{prefix}\tQNAME={} FLAG={} MAPQ={} strand={} start={} end={} CIGAR={} tid={} mtid={} mate_ok={} asm_fail={}",
        qname(rec),
        rec.flags(),
        rec.mapq(),
        strand_label(rec),
        rec.pos() + 1,
        gatk_haplotypecaller::read_unclip::alignment_end_1based(rec),
        rec.cigar(),
        rec.tid(),
        rec.mtid(),
        mate_ok_java(rec),
        !passes_assembly_read(rec, &cfg),
    );
}

fn unique_likelihood_indices(likelihoods: &[RegionReadLikelihood]) -> BTreeSet<usize> {
    likelihoods.iter().map(|c| c.read_index.get()).collect()
}

fn real_pairhmm(likelihoods: &[RegionReadLikelihood], idx: usize) -> bool {
    likelihoods.iter().any(|c| c.read_index.get() == idx)
}

/// Copy of production `per_variant_annotation_likelihoods` (test-only reconstruction).
fn hypothetical_annotation_from_subset(
    subset: &[RegionReadLikelihood],
    reads: &[SharedBamRecord],
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
fn forensic_6r181_source_contracts_no_production_change() {
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");
    let early = include_str!("../src/hc_genotyping_engine/genotype_site_early_template.rs");
    let fin = include_str!("../src/hc_genotyping_engine/genotype_finalize.rs");
    let pipe = include_str!("../src/hc_genotyping_engine/genotype_site_pipeline.rs");
    let ge = include_str!("../src/hc_genotyping_engine/mod.rs");
    let emit = include_str!("../src/region_vcf_emit.rs");
    let ann = include_str!("../src/variant_site_hc_annotations.rs");
    let ll = include_str!("../src/engine_likelihoods.rs");
    assert!(
        early.contains("if is_cluster_tg_snp(&event)")
            && early.contains("finish_strict_java_shaped_site_call"),
        "cluster-TG still selected in SiteEarlyTemplate::try_shaped"
    );
    assert!(
        fin.contains("fn finish_strict_java_shaped_site_call")
            && fin.contains("GenotypedSiteCall::new(event, genotype)"),
        "shaped path still constructs GenotypedSiteCall::new"
    );
    let fin_fn = fin
        .split("fn finish_strict_java_shaped_site_call")
        .nth(1)
        .expect("finish fn");
    let fin_body = fin_fn
        .split("fn finalize_strict_java_variation_genotype")
        .next()
        .expect("body");
    assert!(
        !fin_body.contains("with_annotation_likelihoods"),
        "early-template finish must not attach annotation_likelihoods"
    );
    assert!(
        pipe.contains("SiteEarlyTemplate::try_shaped") && pipe.contains("return Ok(Some(call));"),
        "early template still returns before subset construction"
    );
    let after_early = pipe
        .split("SiteEarlyTemplate::try_shaped")
        .nth(1)
        .expect("try_shaped");
    let subset_after = after_early.find("likelihood_subset_for_event");
    let return_after = after_early.find("return Ok(Some(call));");
    assert!(
        return_after.is_some() && subset_after.is_some_and(|s| return_after.unwrap() < s),
        "cluster-TG return precedes likelihood_subset_for_event"
    );
    assert!(
        emit.contains("if call.annotation_likelihoods.is_empty()")
            && emit.contains("likelihoods: region_wide.likelihoods"),
        "empty annotation_likelihoods still falls back to region-wide"
    );
    assert!(
        ge.contains("fn per_variant_annotation_likelihoods")
            && pipe.contains("with_annotation_likelihoods"),
        "6R.180 helper still exists on the later SiteScore path only"
    );
    assert!(
        ann.contains("6R.174: INFO DP is Java `Coverage.evidenceCount`")
            && ann.contains("mq_rms_of_sample_evidence")
            && ann.contains("strand_bias_sample_counts"),
        "DP/MQ/SOR formulas must stay closed"
    );
    assert!(
        ll.contains("score_pairhmm_from_records_java_mate_contig"),
        "6R.176 mate-contig gate must stay closed"
    );
    assert!(
        emit.contains("INBREEDING_COEFF_MIN_SAMPLES"),
        "6R.178 InbreedingCoeff gate must stay closed"
    );
    assert!(
        !ann.contains("92307333") && !emit.contains("92307333") && !pipe.contains("92307333"),
        "no locus-specific 6R.181 production patch"
    );
    assert!((java_calculate_sor(0, 0, 0, 1) - 1.6094379124341003).abs() < 1e-12);
}

#[test]
fn forensic_6r181_cluster_tg_annotation_membership() {
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    let java_vcf = root.join(JAVA_VCF_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing P12 ref/BAM");
        return;
    }

    kv("java_pin", JAVA_PIN);
    kv("target", "2:92307333 T/G");
    kv("production_change", "NONE");
    kv(
        "classification_target",
        "ANNOTATION_INPUT_MEMBERSHIP_DIVERGENCE",
    );

    if java_vcf.is_file() {
        let (jinfo, jqual, jfmt, jsample) = java_target_record(&java_vcf);
        kv("java_qual", jqual.clone());
        kv("java_format_keys", &jfmt);
        kv("java_sample", jsample.clone());
        kv("java_info_dp", jinfo.get("DP").cloned().unwrap_or_default());
        kv("java_info_mq", jinfo.get("MQ").cloned().unwrap_or_default());
        kv(
            "java_info_sor",
            jinfo.get("SOR").cloned().unwrap_or_default(),
        );
        kv(
            "java_rp",
            if jinfo.contains_key("ReadPosRankSum") {
                jinfo["ReadPosRankSum"].clone()
            } else {
                "OMITTED".into()
            },
        );
        kv(
            "java_ic",
            if jinfo.contains_key("InbreedingCoeff") {
                jinfo["InbreedingCoeff"].clone()
            } else {
                "OMITTED".into()
            },
        );
        assert_eq!(jinfo.get("DP").map(String::as_str), Some("1"));
        assert_eq!(jinfo.get("MQ").map(String::as_str), Some("44.00"));
        assert_eq!(jinfo.get("SOR").map(String::as_str), Some("1.609"));
        assert!(!jinfo.contains_key("ReadPosRankSum"));
        assert!(!jinfo.contains_key("InbreedingCoeff"));
        let fields: Vec<&str> = jfmt.split(':').collect();
        let vals: Vec<&str> = jsample.split(':').collect();
        let get = |name: &str| {
            fields
                .iter()
                .zip(vals.iter())
                .find(|(k, _)| **k == name)
                .map(|(_, v)| *v)
                .unwrap_or("")
        };
        assert_eq!(get("GT"), "1/1");
        assert_eq!(get("AD"), "0,1");
        assert_eq!(get("DP"), "1");
        assert_eq!(get("GQ"), "3");
        assert_eq!(get("PL"), "45,3,0");
        assert!((jqual.parse::<f64>().unwrap_or(0.0) - 35.48).abs() < 0.02);
    }

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
    kv(
        "active_region",
        format!(
            "{}:{}-{}",
            covering.contig,
            covering.start.get(),
            covering.end.get()
        ),
    );
    let outcome = HaplotypeCallerEngine::call_region(
        covering,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("outcome");

    let call = outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(TARGET)
                && c.event.ref_allele == MERGED_REF
                && c.event.alt_allele == MERGED_ALT
        })
        .expect("T/G call");
    kv(
        "tg_format",
        format!(
            "AD={:?} DP={} GQ={} PL={:?}",
            call.genotype.format.ad_as_i32(),
            call.genotype.format.dp.as_i32(),
            call.genotype.format.gq.as_i32(),
            call.genotype.format.pl_as_i32()
        ),
    );
    assert_eq!(call.genotype.format.ad_as_i32(), vec![0, 1]);
    assert_eq!(call.genotype.format.dp.as_i32(), 1);
    assert_eq!(call.genotype.format.gq.as_i32(), 3);
    assert_eq!(call.genotype.format.pl_as_i32(), vec![45, 3, 0]);
    kv(
        "tg_annotation_likelihoods_n_cells",
        call.annotation_likelihoods.len().to_string(),
    );
    let tg_unique: BTreeSet<usize> = unique_likelihood_indices(&call.annotation_likelihoods);
    kv("tg_annotation_unique_n", tg_unique.len().to_string());
    assert!(
        call.annotation_likelihoods.is_empty() && tg_unique.is_empty(),
        "cluster-TG path must leave annotation_likelihoods empty"
    );

    let indel = outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(CLOSED_INDEL)
                && c.event.ref_allele == CLOSED_INDEL_REF
                && c.event.alt_allele == CLOSED_INDEL_ALT
        })
        .expect("TTC/T call");
    let indel_unique: BTreeSet<usize> = unique_likelihood_indices(&indel.annotation_likelihoods);
    kv("ttc_annotation_unique_n", indel_unique.len().to_string());
    assert_eq!(indel_unique.len(), 1, "6R.180 TTC/T attachment stays");

    let region_unique = unique_likelihood_indices(&outcome.read_likelihoods);
    kv("region_wide_unique_n", region_unique.len().to_string());
    kv(
        "genotyping_reads_n",
        outcome.genotyping_reads.len().to_string(),
    );

    let var_end = TARGET;
    let overlap: BTreeSet<usize> = outcome
        .genotyping_reads
        .iter()
        .enumerate()
        .filter(|(_, r)| java_alignment_read_overlaps_interval(r, TARGET, var_end, MARGIN))
        .map(|(i, _)| i)
        .filter(|i| region_unique.contains(i))
        .collect();
    kv("overlap_m2_unique_n", overlap.len().to_string());

    let overlap_cells: Vec<RegionReadLikelihood> = outcome
        .read_likelihoods
        .iter()
        .filter(|c| overlap.contains(&c.read_index.get()))
        .cloned()
        .collect();
    let hypo = hypothetical_annotation_from_subset(&overlap_cells, &outcome.genotyping_reads);
    let hypo_unique = unique_likelihood_indices(&hypo);
    kv(
        "hypothetical_overlap_then_6r180_helper_n",
        hypo_unique.len().to_string(),
    );

    let orig_by_key: BTreeMap<(Vec<u8>, u16), &Record> = covering
        .reads
        .iter()
        .map(|r| ((r.qname().to_vec(), r.flags()), r.as_ref()))
        .collect();
    for rec in covering.reads.iter() {
        dump_read("REGION", rec);
    }
    for rec in outcome.genotyping_reads.iter() {
        dump_read("GENOTYPING", rec);
    }

    let (full_ref, full_pad) = outcome.assembly.event_map_reference();
    let apply_pad = outcome
        .assembly
        .haplotypes
        .iter()
        .find(|h| h.is_reference)
        .and_then(|h| h.genome_loc.as_ref().map(|g| g.start.get()))
        .unwrap_or(full_pad);
    let evidence = HcStrandBiasLikelihoods {
        reads: &outcome.genotyping_reads,
        likelihoods: &outcome.read_likelihoods,
        haplotypes: &outcome.assembly.haplotypes,
        contig: &covering.contig,
        ref_bytes: outcome.assembly.reference_bases(),
        pad_start_1based: apply_pad,
        full_ref_bytes: full_ref,
        full_pad_1based: full_pad,
        max_mnp_distance: outcome.assembly.max_mnp_distance(),
        emit_spanning_dels: true,
    };
    let hap_cache = build_per_haplotype_variation_events(
        evidence.haplotypes,
        evidence.full_ref_bytes,
        evidence.full_pad_1based,
        evidence.max_mnp_distance,
        evidence.contig,
    );
    let mapping = create_allele_mapper_with_events(
        &VariationEvent::from_alleles(evidence.contig, TARGET, MERGED_REF, MERGED_ALT),
        TARGET,
        evidence.haplotypes,
        evidence.pad_start_1based,
        evidence.ref_bytes,
        evidence.max_mnp_distance,
        evidence.emit_spanning_dels,
        Some(&hap_cache),
    );
    kv(
        "haplotypes",
        format!(
            "n={} ref_pool={:?} alt_pool={:?}",
            evidence.haplotypes.len(),
            mapping
                .ref_haplotype_indices
                .iter()
                .map(|i| i.get())
                .collect::<Vec<_>>(),
            mapping
                .alt_haplotype_indices
                .iter()
                .map(|i| i.get())
                .collect::<Vec<_>>()
        ),
    );
    let rows = region_likelihoods_to_rows(evidence.likelihoods, evidence.haplotypes.len());
    let marg = marginalize_rows_to_biallelic_alleles(
        &rows,
        &mapping.ref_haplotype_indices,
        &mapping.alt_haplotype_indices,
    );
    let allele_by_idx: BTreeMap<usize, String> = marg
        .iter()
        .map(|row| {
            let lls = &row.haplotype_log10_likelihoods;
            let ll_ref = lls.first().copied().unwrap_or(f64::NEG_INFINITY);
            let ll_alt = lls.get(1).copied().unwrap_or(f64::NEG_INFINITY);
            let (best_is_ref, best, second) = if ll_ref > ll_alt {
                (true, ll_ref, ll_alt)
            } else if ll_alt > ll_ref {
                (false, ll_alt, ll_ref)
            } else {
                (true, ll_ref, ll_alt)
            };
            let gap = if second.is_finite() {
                best - second
            } else {
                f64::INFINITY
            };
            let informative = gap > LOG_10_INFORMATIVE_THRESHOLD;
            let allele = if !informative {
                "NON_INF"
            } else if best_is_ref {
                "REF"
            } else {
                "ALT"
            };
            (row.read_index, allele.to_string())
        })
        .collect();

    let mut mapq44: Vec<String> = Vec::new();
    for &idx in &region_unique {
        let Some(rec) = outcome.genotyping_reads.get(idx) else {
            continue;
        };
        let orig = orig_by_key
            .get(&(rec.qname().to_vec(), rec.flags()))
            .copied();
        let mate_ok = orig.map(mate_ok_java).unwrap_or_else(|| mate_ok_java(rec));
        let ov = overlap.contains(&idx);
        let in_tg_ann = tg_unique.contains(&idx);
        let in_ttc_ann = indel_unique.contains(&idx);
        let in_hypo = hypo_unique.contains(&idx);
        let allele = allele_by_idx
            .get(&idx)
            .cloned()
            .unwrap_or_else(|| "?".into());
        if rec.mapq() == 44 {
            mapq44.push(format!("QNAME={} FLAG={}", qname(rec), rec.flags()));
        }
        kv(
            "LL_READ",
            format!(
                "idx={idx} QNAME={} FLAG={} MAPQ={} strand={} start={} end={} CIGAR={} mate_ok={} pairhmm_real={} overlap_tg={} allele={} tg_annotation={} ttc_annotation={} hypo_overlap_helper={} contributes_fallback_dp={} contributes_fallback_mq={}",
                qname(rec),
                rec.flags(),
                rec.mapq(),
                strand_label(rec),
                rec.pos() + 1,
                gatk_haplotypecaller::read_unclip::alignment_end_1based(rec),
                rec.cigar(),
                mate_ok,
                real_pairhmm(&outcome.read_likelihoods, idx),
                ov,
                allele,
                in_tg_ann,
                in_ttc_ann,
                in_hypo,
                ov,
                rec.mapq() != MAPPING_QUALITY_UNAVAILABLE,
            ),
        );
    }
    kv("mapq44_reads", mapq44.join(" | "));
    assert_eq!(
        mapq44.len(),
        1,
        "Java MQ=44.00 uniquely identifies one read"
    );
    assert!(mapq44[0].contains(JAVA_MQ44_QNAME) && mapq44[0].contains("FLAG=83"));

    let fallback_dp = coverage_evidence_count(
        &outcome.genotyping_reads,
        &outcome.read_likelihoods,
        TARGET,
        var_end,
        MARGIN,
    );
    let fallback_mqs =
        rms_mapping_quality_sample_mapqs(&outcome.genotyping_reads, &outcome.read_likelihoods);
    let fallback_mq = rms_mapping_quality_raw(&fallback_mqs);
    kv("fallback_dp", fallback_dp.to_string());
    kv("fallback_mapqs", format!("{fallback_mqs:?}"));
    kv("fallback_mq", format!("{fallback_mq:?}"));
    assert_eq!(fallback_dp, 1);
    assert_eq!(region_unique.len(), 2);
    assert_eq!(
        overlap.len(),
        1,
        "stored n=2 retainEvidence overlap is the MAPQ=44 read; annotation still not attached"
    );
    assert_eq!(
        hypo_unique.len(),
        1,
        "overlap on stored n=2 is n=1; 6R.181 B still does not attach it"
    );
    let fallback_mq_printed = fallback_mq.map(|t| (t.2 * 100.0).round() / 100.0);
    assert_eq!(fallback_mq_printed, Some(42.05));

    let sor_table = strand_bias_contingency_table(
        &evidence,
        TARGET,
        MERGED_REF,
        MERGED_ALT,
        STRAND_ODDS_RATIO_MIN_COUNT,
    );
    kv(
        "fallback_sor_table",
        format!(
            "[{},{};{},{}]",
            sor_table.0, sor_table.1, sor_table.2, sor_table.3
        ),
    );
    let fallback_sor = java_calculate_sor(sor_table.0, sor_table.1, sor_table.2, sor_table.3);
    kv("fallback_sor", format!("{fallback_sor}"));
    assert!((fallback_sor - 0.6931471805599453).abs() < 1e-9);

    kv(
        "java_annotation_n",
        "1 — MAPQ=44 FLAG=83 ALT_REV (Coverage/MQ/SOR observables; independent of TTC/T object)",
    );
    kv(
        "rust_variant_local",
        "None — GenotypedSiteCall.annotation_likelihoods empty; likelihood_subset_for_event never ran",
    );
    kv(
        "first_arrow",
        "SiteEarlyTemplate::try_shaped is_cluster_tg_snp → finish_strict_java_shaped_site_call → GenotypedSiteCall::new without subset construction",
    );
    kv("classification", "B — MISSING OBJECT CONSTRUCTION");

    let emitted =
        try_emit_call_region_variants(covering, &outcome, "NA12878", DEFAULT_STAND_EMIT_CONFIDENCE)
            .expect("emit");
    let rec = emitted
        .iter()
        .find(|r| {
            r.position == TARGET
                && r.reference == MERGED_REF
                && r.alternate.first().map(String::as_str) == Some(MERGED_ALT)
        })
        .expect("T/G record");
    let sample = rec.samples.first().expect("sample");
    assert_eq!(
        sample.gt.as_ref().map(|g| g.to_string()).as_deref(),
        Some("1/1")
    );
    assert_eq!(sample.ad.as_deref(), Some(&[0u32, 1][..]));
    assert_eq!(sample.dp.map(|d| d as i32), Some(1));
    assert_eq!(sample.gq.map(|g| g as i32), Some(3));
    assert_eq!(sample.pl.as_deref(), Some(&[45u32, 3, 0][..]));
    assert!((rec.quality.unwrap_or(0.0) - 35.48).abs() < 0.05);
    assert_eq!(info_i32(&rec.info, "DP"), Some(1));
    let mq = info_f64(&rec.info, "MQ").unwrap_or(-1.0);
    let sor = info_f64(&rec.info, "SOR").unwrap_or(-1.0);
    kv("emit_mq", format!("{mq}"));
    kv("emit_sor", format!("{sor}"));
    assert!((mq - 42.05).abs() < 0.005, "MQ={mq}");
    assert!((sor - 0.6931471805599453).abs() < 1e-9, "SOR={sor}");

    let closed_indel = emitted
        .iter()
        .find(|r| r.position == CLOSED_INDEL && r.reference == CLOSED_INDEL_REF)
        .expect("TTC/T");
    assert_eq!(info_i32(&closed_indel.info, "DP"), Some(1));
    let closed_mq = info_f64(&closed_indel.info, "MQ").unwrap_or(-1.0);
    let closed_sor = info_f64(&closed_indel.info, "SOR").unwrap_or(-1.0);
    assert!((closed_mq - 44.0).abs() < 0.005);
    assert!((closed_sor - 1.6094379124341003).abs() < 1e-3 || (closed_sor - 1.609).abs() < 0.002);

    let closed_specs = parse_intervals_cli_string(&dict, "2:92305500-92305850").expect("closed");
    let closed_walk = traverse_assembly_region_walker(
        &dict,
        &closed_specs,
        &ref_fasta,
        &bam,
        &ReadFilterParams::gatk_standard_hc(),
        &WalkerTraversalConfig::gatk_haplotype_caller_production(100),
    )
    .expect("closed walk");
    let closed_regions = flatten_assembly_regions(&closed_walk);
    let closed_covering = closed_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_SNP
                && r.end.get() >= CLOSED_SNP
        })
        .expect("closed covering");
    let closed_outcome = HaplotypeCallerEngine::call_region(
        closed_covering,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("closed call")
    .expect("closed outcome");
    let closed_emitted = try_emit_call_region_variants(
        closed_covering,
        &closed_outcome,
        "NA12878",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .expect("closed emit");
    let closed_rec = closed_emitted
        .iter()
        .find(|r| r.position == CLOSED_SNP)
        .expect("closed G/T");
    assert_eq!(info_i32(&closed_rec.info, "DP"), Some(3), "6R.174 stays");
    assert!(
        !info_has(&closed_rec.info, "InbreedingCoeff"),
        "6R.178 stays"
    );
    let closed_snp_sor = info_f64(&closed_rec.info, "SOR").unwrap_or(-1.0);
    assert!(
        (closed_snp_sor - 0.6931471805599453).abs() < 1e-6
            || (closed_snp_sor - 0.693).abs() < 0.002
    );
    kv(
        "closed_snp",
        format!(
            "DP={} SOR={closed_snp_sor} IC_omitted",
            info_i32(&closed_rec.info, "DP").unwrap_or(-1)
        ),
    );
}
