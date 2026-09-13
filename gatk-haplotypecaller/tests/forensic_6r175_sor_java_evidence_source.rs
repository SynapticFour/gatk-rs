//! 6R.175: first causal SOR arrow at `2:92305634 G/T` (proof-only).
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! 6R.174 INFO DP stays closed. FORMAT/QUAL/FS/MQ/ReadPosRankSum stay closed.
//! PRODUCTION CHANGE: NONE.
//!
//! Java SOR 0.693 from informative table `[0,0;1,1]` (the two 53752 mates).
//! 6R.175 proved Rust `[0,0;1,2]` because FLAG=145 (mate on contig 9) was
//! still PairHMM-informative ALT_REV. 6R.176 closed that membership gate.
//! Formula is the same.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r175_sor_java_evidence_source -- --nocapture --test-threads=1
//! HOLDOUT_6R175=1 cargo test -p gatk-haplotypecaller --test holdout_6r175_sor -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::{build_per_haplotype_variation_events, VariationEvent};
use gatk_haplotypecaller::fragment_overlap::read_base_at_ref_coord_1based;
use gatk_haplotypecaller::hc_allele_mapping::create_allele_mapper_with_events;
use gatk_haplotypecaller::hc_genotyping_engine::java_alignment_read_overlaps_interval;
use gatk_haplotypecaller::read_assembly_filter::{passes_assembly_read, AssemblyReadFilterConfig};
use gatk_haplotypecaller::read_realignment::LOG_10_INFORMATIVE_THRESHOLD;
use gatk_haplotypecaller::variant_site_hc_annotations::{
    strand_bias_contingency_table, HcStrandBiasLikelihoods, STRAND_ODDS_RATIO_MIN_COUNT,
};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, marginalize_rows_to_biallelic_alleles,
    region_likelihoods_to_rows, traverse_assembly_region_walker, try_emit_call_region_variants,
    AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition, HaplotypeCallerEngine,
    HcGenotypingConfig, ReadFilterParams, WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use rust_htslib::bam::Record;
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

const INTERVAL: &str = "2:92305500-92305850";
const CLOSED_INTERVAL: &str = "2:92316200-92316580";
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 92_305_634;
const CLOSED: u64 = 92_316_347;
const MERGED_REF: &str = "G";
const MERGED_ALT: &str = "T";
const FLAG_REVERSE: u16 = 0x10;
const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
const THIRD_READ_QNAME: &str = "H06JUADXX130110:1:1101:10018:4569";

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R175\t{key}\t{}", value.as_ref());
}

/// Java `StrandOddsRatio.calculateSOR` (pseudocount 1, natural log).
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

fn strand_label(rec: &Record) -> &'static str {
    if rec.flags() & FLAG_REVERSE != 0 {
        "rev"
    } else {
        "fwd"
    }
}

fn qname(rec: &Record) -> String {
    String::from_utf8_lossy(rec.qname()).into_owned()
}

fn mate_on_same_contig_or_unmapped(rec: &Record) -> bool {
    if !rec.is_paired() || rec.is_mate_unmapped() || rec.is_unmapped() {
        return true;
    }
    rec.tid() == rec.mtid()
}

fn dump_read(prefix: &str, rec: &Record) {
    let cfg = AssemblyReadFilterConfig::gatk_defaults();
    eprintln!(
        "6R175\t{prefix}\tQNAME={} FLAG={} MAPQ={} strand={} start={} end={} CIGAR={} tid={} mtid={} mate_unmapped={} mate_ok={} asm_fail={}",
        qname(rec),
        rec.flags(),
        rec.mapq(),
        strand_label(rec),
        rec.pos() + 1,
        gatk_haplotypecaller::read_unclip::alignment_end_1based(rec),
        rec.cigar(),
        rec.tid(),
        rec.mtid(),
        rec.is_mate_unmapped(),
        mate_on_same_contig_or_unmapped(rec),
        !passes_assembly_read(rec, &cfg),
    );
}

#[test]
fn forensic_6r175_java_sor_formula_contract() {
    kv("java_pin", JAVA_PIN);
    kv(
        "java_path",
        "StrandOddsRatio.calculateAnnotationFromLikelihoods → StrandBiasTest.getContingencyTable(MIN_COUNT=0) → calculateSOR",
    );
    assert_eq!(STRAND_ODDS_RATIO_MIN_COUNT, 0);
    assert_eq!(LOG_10_INFORMATIVE_THRESHOLD, 0.2);
    let balanced = java_calculate_sor(0, 0, 1, 1);
    assert!(
        (balanced - 0.6931471805599453).abs() < 1e-12,
        "Java [0,0;1,1] must be ln(2)≈0.693, got {balanced}"
    );
    let rust_like = java_calculate_sor(0, 0, 2, 1);
    assert!(
        (rust_like - 1.1786549963416462).abs() < 1e-9,
        "Java [0,0;2,1] must print SOR=1.179, got {rust_like}"
    );
    let rust_like_swap = java_calculate_sor(0, 0, 1, 2);
    assert!(
        (rust_like_swap - rust_like).abs() < 1e-12,
        "[0,0;1,2] is formula-symmetric with [0,0;2,1]"
    );
    let src = include_str!("../src/variant_site_hc_annotations.rs");
    assert!(
        src.contains("6R.166: FS/SOR use Java `getContingencyTable`")
            && src.contains("strand_bias_sample_counts"),
        "production SOR still uses the 6R.166 helper"
    );
    assert!(
        !src.contains("92305634"),
        "no locus-specific SOR / DP patch"
    );
    let sor_src = include_str!("../src/annotator/plugins/strand_odds_ratio.rs");
    assert!(
        sor_src.contains("PSEUDOCOUNT: f64 = 1.0") && sor_src.contains(".ln()"),
        "production calculateSOR is Java pseudocount-1 + ln"
    );
}

#[test]
fn forensic_6r175_sor_java_evidence_source() {
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing P12 ref/BAM");
        return;
    }
    kv("java_pin", JAVA_PIN);
    kv("target", "2:92305634 G/T");
    kv("java_sor", "0.693");
    kv("rust_sor_before_6r176", "1.179");
    kv("production_change", "NONE");

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
        .expect("ActiveFull covering target");

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
        .expect("genotyped G/T");
    let fmt = &call.genotype.format;
    assert_eq!(fmt.ad_as_i32(), vec![0, 2], "FORMAT AD closed");
    assert_eq!(fmt.pl_as_i32(), vec![90, 6, 0], "FORMAT PL closed");
    assert_eq!(fmt.dp.as_i32(), 2, "FORMAT DP closed");
    assert_eq!(fmt.gq.as_i32(), 6, "FORMAT GQ closed");

    let cfg = HcGenotypingConfig::default();
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
        emit_spanning_dels: !cfg.disable_spanning_event_genotyping,
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
    for (i, h) in evidence.haplotypes.iter().enumerate() {
        let events: Vec<String> = hap_cache
            .events_for(i)
            .iter()
            .map(|e| {
                format!(
                    "{}:{} {}/{}",
                    e.contig,
                    e.start_1based.get(),
                    e.ref_allele,
                    e.alt_allele
                )
            })
            .collect();
        kv(
            "hap",
            format!(
                "i={i} is_ref={} events={}",
                h.is_reference,
                events.join(",")
            ),
        );
    }

    let rows = region_likelihoods_to_rows(evidence.likelihoods, evidence.haplotypes.len());
    let marg = marginalize_rows_to_biallelic_alleles(
        &rows,
        &mapping.ref_haplotype_indices,
        &mapping.alt_haplotype_indices,
    );

    let likelihood_reads: BTreeSet<usize> = rows.iter().map(|r| r.read_index).collect();
    let mut retain = BTreeSet::new();
    for cell in evidence.likelihoods {
        let idx = cell.read_index.get();
        let Some(rec) = evidence.reads.get(idx) else {
            continue;
        };
        if java_alignment_read_overlaps_interval(
            rec,
            TARGET,
            TARGET,
            cfg.informative_read_overlap_margin,
        ) {
            retain.insert(idx);
        }
    }

    kv("region_reads", covering.reads.len().to_string());
    kv("genotyping_reads", evidence.reads.len().to_string());
    kv(
        "likelihood_unique_reads",
        likelihood_reads.len().to_string(),
    );
    kv("retainEvidence_unique", retain.len().to_string());
    for rec in evidence.reads {
        dump_read("GENOTYPING", rec);
    }

    let mut pileup = (0u32, 0u32, 0u32, 0u32);
    let mut full_table = (0u32, 0u32, 0u32, 0u32);
    let mut retain_table = (0u32, 0u32, 0u32, 0u32);
    let mut informative_members = Vec::new();
    let mut retain_informative = Vec::new();
    let mut seen_qnames: BTreeMap<String, usize> = BTreeMap::new();

    for rec in &covering.reads {
        dump_read("REGION", rec);
        *seen_qnames.entry(qname(rec)).or_insert(0) += 1;
        let Some(base) = read_base_at_ref_coord_1based(rec, TARGET as i32) else {
            continue;
        };
        let reverse = rec.flags() & FLAG_REVERSE != 0;
        let allele = if base.eq_ignore_ascii_case(&b'T') {
            "ALT"
        } else if base.eq_ignore_ascii_case(&b'G') {
            "REF"
        } else {
            "OTHER"
        };
        match (allele, reverse) {
            ("REF", false) => pileup.0 += 1,
            ("REF", true) => pileup.1 += 1,
            ("ALT", false) => pileup.2 += 1,
            ("ALT", true) => pileup.3 += 1,
            _ => {}
        }
        kv(
            "pileup_overlap",
            format!(
                "QNAME={} FLAG={} strand={} pileup={} cell={}",
                qname(rec),
                rec.flags(),
                strand_label(rec),
                base as char,
                allele
            ),
        );
    }
    for (name, n) in &seen_qnames {
        if *n > 1 {
            kv("multi_mate_qname", format!("{name} n={n}"));
        }
    }

    for row in &marg {
        let Some(rec) = evidence.reads.get(row.read_index) else {
            continue;
        };
        let lls = &row.haplotype_log10_likelihoods;
        let ll_ref = lls.first().copied().unwrap_or(f64::NEG_INFINITY);
        let ll_alt = lls.get(1).copied().unwrap_or(f64::NEG_INFINITY);
        let overlaps = java_alignment_read_overlaps_interval(
            rec,
            TARGET,
            TARGET,
            cfg.informative_read_overlap_margin,
        );
        let in_likelihoods = likelihood_reads.contains(&row.read_index);
        let pileup_base = read_base_at_ref_coord_1based(rec, TARGET as i32)
            .map(|b| (b as char).to_string())
            .unwrap_or_else(|| "-".to_string());
        if !ll_ref.is_finite() && !ll_alt.is_finite() {
            dump_read("LIKELIHOOD_NONFINITE", rec);
            continue;
        }
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
        let reverse = rec.flags() & FLAG_REVERSE != 0;
        let allele = if best_is_ref { "REF" } else { "ALT" };
        let cell = if !informative {
            "none"
        } else if best_is_ref {
            if reverse {
                "REF_REV"
            } else {
                "REF_FWD"
            }
        } else if reverse {
            "ALT_REV"
        } else {
            "ALT_FWD"
        };
        kv(
            "likelihood_read",
            format!(
                "QNAME={} FLAG={} MAPQ={} strand={} start={} end={} CIGAR={} pileup={} ll_ref={ll_ref:.4} ll_alt={ll_alt:.4} gap={gap:.4} informative={informative} allele={allele} cell={cell} retain={overlaps} in_matrix={in_likelihoods}",
                qname(rec),
                rec.flags(),
                rec.mapq(),
                strand_label(rec),
                rec.pos() + 1,
                gatk_haplotypecaller::read_unclip::alignment_end_1based(rec),
                rec.cigar(),
                pileup_base
            ),
        );
        if informative {
            match (best_is_ref, reverse) {
                (true, false) => full_table.0 += 1,
                (true, true) => full_table.1 += 1,
                (false, false) => full_table.2 += 1,
                (false, true) => full_table.3 += 1,
            }
            informative_members.push(format!(
                "{} FLAG={} strand={} allele={allele} cell={cell}",
                qname(rec),
                rec.flags(),
                strand_label(rec)
            ));
            if overlaps {
                match (best_is_ref, reverse) {
                    (true, false) => retain_table.0 += 1,
                    (true, true) => retain_table.1 += 1,
                    (false, false) => retain_table.2 += 1,
                    (false, true) => retain_table.3 += 1,
                }
                retain_informative.push(format!(
                    "{} FLAG={} strand={} allele={allele} cell={cell}",
                    qname(rec),
                    rec.flags(),
                    strand_label(rec)
                ));
            }
        } else if overlaps {
            kv(
                "retain_uninformative",
                format!(
                    "QNAME={} FLAG={} strand={} gap={gap:.4} pileup={}",
                    qname(rec),
                    rec.flags(),
                    strand_label(rec),
                    pileup_base
                ),
            );
        }
        if qname(rec) == THIRD_READ_QNAME {
            kv(
                "third_read",
                format!(
                    "FLAG={} strand={} pileup={} informative={informative} cell={cell} retain={overlaps} gap={gap:.4}",
                    rec.flags(),
                    strand_label(rec),
                    pileup_base
                ),
            );
        }
    }

    let production = strand_bias_contingency_table(
        &evidence,
        TARGET,
        MERGED_REF,
        MERGED_ALT,
        STRAND_ODDS_RATIO_MIN_COUNT,
    );
    kv(
        "rust_production_table",
        format!(
            "[{},{};{},{}]",
            production.0, production.1, production.2, production.3
        ),
    );
    kv(
        "full_matrix_informative",
        format!(
            "[{},{};{},{}] members={informative_members:?}",
            full_table.0, full_table.1, full_table.2, full_table.3
        ),
    );
    kv(
        "retainEvidence_informative",
        format!(
            "[{},{};{},{}] members={retain_informative:?}",
            retain_table.0, retain_table.1, retain_table.2, retain_table.3
        ),
    );
    kv(
        "pileup_table",
        format!("[{},{};{},{}]", pileup.0, pileup.1, pileup.2, pileup.3),
    );

    let sor_prod = java_calculate_sor(production.0, production.1, production.2, production.3);
    let sor_retain = java_calculate_sor(
        retain_table.0,
        retain_table.1,
        retain_table.2,
        retain_table.3,
    );
    let sor_pileup = java_calculate_sor(pileup.0, pileup.1, pileup.2, pileup.3);
    kv("formula_on_production_table", format!("{sor_prod:.9}"));
    kv("formula_on_retain_table", format!("{sor_retain:.9}"));
    kv("formula_on_pileup_table", format!("{sor_pileup:.9}"));
    kv(
        "formula_on_java_expected_[0,0;1,1]",
        format!("{:.9}", java_calculate_sor(0, 0, 1, 1)),
    );

    let emitted =
        try_emit_call_region_variants(covering, &outcome, "NA12878", DEFAULT_STAND_EMIT_CONFIDENCE)
            .expect("emit");
    let rec = emitted
        .iter()
        .find(|r| r.position == TARGET)
        .expect("emitted");
    let sample = rec.samples.first().expect("sample");
    let gt = sample
        .gt
        .as_ref()
        .map(|g| {
            g.alleles
                .iter()
                .map(|a| a.to_string())
                .collect::<Vec<_>>()
                .join("/")
        })
        .unwrap_or_default();
    let ad = sample
        .ad
        .as_ref()
        .map(|v| {
            v.iter()
                .map(|x| x.to_string())
                .collect::<Vec<_>>()
                .join(",")
        })
        .unwrap_or_default();
    let pl = sample
        .pl
        .as_ref()
        .map(|v| {
            v.iter()
                .map(|x| x.to_string())
                .collect::<Vec<_>>()
                .join(",")
        })
        .unwrap_or_default();
    let fs = info_f64(&rec.info, "FS").unwrap_or(-1.0);
    let sor = info_f64(&rec.info, "SOR").unwrap_or(-1.0);
    let mq = info_f64(&rec.info, "MQ").unwrap_or(-1.0);
    kv(
        "emit",
        format!(
            "GT={gt} AD={ad} FORMAT_DP={:?} GQ={:?} PL={pl} QUAL={:?} INFO_DP={:?} FS={fs} SOR={sor} MQ={mq} RP={:?}",
            sample.dp,
            sample.gq,
            rec.quality,
            info_i32(&rec.info, "DP"),
            info_f64(&rec.info, "ReadPosRankSum")
        ),
    );
    assert_eq!(gt, "1/1");
    assert_eq!(ad, "0,2");
    assert_eq!(sample.dp.map(|d| d as i32), Some(2));
    assert_eq!(info_i32(&rec.info, "DP"), Some(3), "6R.174 INFO DP stays 3");
    assert!(
        (rec.quality.unwrap_or(0.0) - 78.32).abs() < 0.02,
        "QUAL stays 78.32, got {:?}",
        rec.quality
    );
    assert!(
        (sor - 0.6931471805599453).abs() < 1e-6 || (sor - 0.693).abs() < 0.002,
        "after 6R.176, production SOR is Java 0.693, got {sor}"
    );
    assert_eq!(
        production, full_table,
        "production helper must equal the dumped full-matrix informative table"
    );
    assert_eq!(
        production,
        (0, 0, 1, 1),
        "after 6R.176, production SOR table is Java [0,0;1,1]"
    );
    assert_eq!(
        retain_table,
        (0, 0, 1, 1),
        "retainEvidence informative table matches Java after excluding FLAG=145"
    );
    assert!(
        !informative_members
            .iter()
            .any(|m| m.contains(THIRD_READ_QNAME)),
        "FLAG=145 must not be an informative SOR cell after 6R.176, got {informative_members:?}"
    );
    assert!(
        (java_calculate_sor(production.0, production.1, production.2, production.3) - 0.693).abs()
            < 0.002,
        "Java formula on the production table must print 0.693"
    );
    assert!(
        (java_calculate_sor(0, 0, 1, 1) - 0.693).abs() < 0.002,
        "Java formula on [0,0;1,1] must reproduce 0.693"
    );

    // Closed 6R.166 site must not regress.
    let specs2 = parse_intervals_cli_string(&dict, CLOSED_INTERVAL).expect("closed interval");
    let walk2 = traverse_assembly_region_walker(
        &dict,
        &specs2,
        &ref_fasta,
        &bam,
        &ReadFilterParams::gatk_standard_hc(),
        &WalkerTraversalConfig::gatk_haplotype_caller_production(100),
    )
    .expect("walk2");
    let regions2 = flatten_assembly_regions(&walk2);
    let covering2 = regions2
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED
                && r.end.get() >= CLOSED
        })
        .expect("closed covering");
    let outcome2 = HaplotypeCallerEngine::call_region(
        covering2,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call2")
    .expect("outcome2");
    let emitted2 = try_emit_call_region_variants(
        covering2,
        &outcome2,
        "NA12878",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .expect("emit2");
    let rec2 = emitted2
        .iter()
        .find(|r| r.position == CLOSED)
        .expect("closed emit");
    let fs2 = info_f64(&rec2.info, "FS").unwrap_or(-1.0);
    let sor2 = info_f64(&rec2.info, "SOR").unwrap_or(-1.0);
    let mq2 = info_f64(&rec2.info, "MQ").unwrap_or(-1.0);
    assert!((fs2 - 0.0).abs() < 1e-6, "6R.166 FS stays 0");
    assert!((sor2 - 1.179).abs() < 0.002, "6R.166 SOR stays 1.179");
    assert!((mq2 - 40.25).abs() < 1e-9, "6R.168 MQ stays 40.25");
    kv("closed_92316347", "FORMAT/FS/SOR/MQ unchanged");
}
