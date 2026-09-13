//! 6R.165: first causal INFO annotation arrows at `2:92316347 G/A`.
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`, GKL 0.8.8.
//! Production change in 6R.165: NONE. 6R.164 haplotype materialization stays closed.
//! 6R.166 later switched production FS/SOR onto the Java likelihoods table;
//! this file still documents the pileup `[0,2;3,2]` vs likelihoods `[0,0;2,1]`
//! split. MQ production evidence moved in 6R.167; the pileup mean 36.4 remains
//! the 6R.165 reconstructed ALT-pileup value.
//!
//! Java INFO: FS=0 MQ=40.25 SOR=1.179
//! 6R.165 pileup reconstruction: FS≈3.68 MQ=36.4 SOR≈0.636
//! FORMAT/QUAL already match and are not reopened.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r165_info_annotation_divergence -- --test-threads=1 --nocapture
//! HOLDOUT_6R165=1 cargo test -p gatk-haplotypecaller --test holdout_6r165_info_annotation -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::{build_per_haplotype_variation_events, VariationEvent};
use gatk_haplotypecaller::fragment_overlap::read_base_at_ref_coord_1based;
use gatk_haplotypecaller::hc_allele_mapping::create_allele_mapper_with_events;
use gatk_haplotypecaller::read_realignment::LOG_10_INFORMATIVE_THRESHOLD;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, region_likelihoods_to_rows,
    traverse_assembly_region_walker, try_emit_call_region_variants, AssemblyRegionCallDisposition,
    CallRegionArgs, GenomePosition, HaplotypeCallerEngine, ReadFilterParams, ReadLikelihoodRow,
    WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use rust_htslib::bam::Record;
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

const INTERVAL: &str = "2:92316200-92316580";
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 92_316_347;
const MERGED_REF: &str = "G";
const MERGED_ALT: &str = "A";
const FLAG_REVERSE: u16 = 0x10;
const FLAG_SECONDARY: u16 = 0x100;
const FLAG_DUP: u16 = 0x400;
const FLAG_SUPPLEMENTARY: u16 = 0x800;
const MAPPING_QUALITY_UNAVAILABLE: u8 = 255;
const JAVA_FS_MIN_COUNT: u32 = 2;
const JAVA_SOR_MIN_COUNT: u32 = 0;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
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

/// Java `RMSMappingQuality.makeFinalizedAnnotationString`: `sqrt(sum(mq²)/n)`.
fn java_rms_mq(mapqs: &[u8]) -> (usize, u64, f64) {
    let usable: Vec<u8> = mapqs
        .iter()
        .copied()
        .filter(|&mq| mq != MAPPING_QUALITY_UNAVAILABLE)
        .collect();
    let n = usable.len();
    let sum_sq: u64 = usable.iter().map(|&mq| u64::from(mq) * u64::from(mq)).sum();
    let rms = if n == 0 {
        0.0
    } else {
        (sum_sq as f64 / n as f64).sqrt()
    };
    (n, sum_sq, rms)
}

fn read_id(rec: &Record) -> String {
    format!(
        "QNAME={} FLAG={} pos={} CIGAR={} MAPQ={}",
        String::from_utf8_lossy(rec.qname()),
        rec.flags(),
        rec.pos() + 1,
        rec.cigar(),
        rec.mapq()
    )
}

fn strand_label(rec: &Record) -> &'static str {
    if rec.flags() & FLAG_REVERSE != 0 {
        "rev"
    } else {
        "fwd"
    }
}

struct StrandTable {
    ref_fw: u32,
    ref_rv: u32,
    alt_fw: u32,
    alt_rv: u32,
    members: Vec<String>,
}

impl StrandTable {
    fn cells(&self) -> (u32, u32, u32, u32) {
        (self.ref_fw, self.ref_rv, self.alt_fw, self.alt_rv)
    }
    fn total(&self) -> u32 {
        self.ref_fw + self.ref_rv + self.alt_fw + self.alt_rv
    }
}

fn info_key(v: &InfoValue) -> &str {
    match v {
        InfoValue::Flag(k)
        | InfoValue::Integer(k, _)
        | InfoValue::Float(k, _)
        | InfoValue::String(k, _)
        | InfoValue::Character(k, _) => k.as_str(),
    }
}

/// Rust production `read_strand_evidence_at_site`: pileup G/A on `region.reads`.
fn rust_pileup_strand_table<'a>(
    reads: impl IntoIterator<Item = &'a Record>,
) -> (StrandTable, Vec<u8>) {
    let mut table = StrandTable {
        ref_fw: 0,
        ref_rv: 0,
        alt_fw: 0,
        alt_rv: 0,
        members: Vec::new(),
    };
    let mut alt_mapqs = Vec::new();
    let ref_b = MERGED_REF.as_bytes()[0];
    let alt_b = MERGED_ALT.as_bytes()[0];
    let pos = TARGET as i32;
    for rec in reads {
        let Some(base) = read_base_at_ref_coord_1based(rec, pos) else {
            continue;
        };
        let reverse = rec.flags() & FLAG_REVERSE != 0;
        let supports_alt = base.eq_ignore_ascii_case(&alt_b);
        let supports_ref = base.eq_ignore_ascii_case(&ref_b);
        if supports_alt {
            if reverse {
                table.alt_rv += 1;
            } else {
                table.alt_fw += 1;
            }
            alt_mapqs.push(rec.mapq());
            table.members.push(format!(
                "{} strand={} allele=ALT pileup={}",
                read_id(rec),
                strand_label(rec),
                base as char
            ));
        } else if supports_ref {
            if reverse {
                table.ref_rv += 1;
            } else {
                table.ref_fw += 1;
            }
            table.members.push(format!(
                "{} strand={} allele=REF pileup={}",
                read_id(rec),
                strand_label(rec),
                base as char
            ));
        }
    }
    (table, alt_mapqs)
}

/// Java `StrandBiasTest.getContingencyTable`: informative best allele + strand of evidence.
fn java_informative_strand_table(
    reads: &[Record],
    rows: &[ReadLikelihoodRow],
    ref_haps: &[usize],
    alt_haps: &[usize],
    min_count: u32,
) -> StrandTable {
    let mut sample = StrandTable {
        ref_fw: 0,
        ref_rv: 0,
        alt_fw: 0,
        alt_rv: 0,
        members: Vec::new(),
    };
    for row in rows {
        let Some(rec) = reads.get(row.read_index) else {
            continue;
        };
        let lls = &row.haplotype_log10_likelihoods;
        let ll_ref = ref_haps
            .iter()
            .filter_map(|&i| lls.get(i).copied())
            .fold(f64::NEG_INFINITY, f64::max);
        let ll_alt = alt_haps
            .iter()
            .filter_map(|&i| lls.get(i).copied())
            .fold(f64::NEG_INFINITY, f64::max);
        if !ll_ref.is_finite() && !ll_alt.is_finite() {
            continue;
        }
        let (best_is_ref, best, second) = if ll_ref > ll_alt {
            (true, ll_ref, ll_alt)
        } else if ll_alt > ll_ref {
            (false, ll_alt, ll_ref)
        } else {
            // Java `bestAllelesBreakingTies`: REF priority 1.0 vs ALT 0 when tied.
            (true, ll_ref, ll_alt)
        };
        let gap = if second.is_finite() {
            best - second
        } else {
            f64::INFINITY
        };
        // Near-tie: Java `searchBestAllele` prefers REF when gap < 0.2.
        let allele_is_ref = if gap < LOG_10_INFORMATIVE_THRESHOLD {
            true
        } else {
            best_is_ref
        };
        let informative = gap > LOG_10_INFORMATIVE_THRESHOLD;
        if !informative {
            continue;
        }
        let reverse = rec.flags() & FLAG_REVERSE != 0;
        if allele_is_ref {
            if reverse {
                sample.ref_rv += 1;
            } else {
                sample.ref_fw += 1;
            }
        } else if reverse {
            sample.alt_rv += 1;
        } else {
            sample.alt_fw += 1;
        }
        sample.members.push(format!(
            "{} strand={} allele={} gap={gap:.4} informative=true",
            read_id(rec),
            strand_label(rec),
            if allele_is_ref { "REF" } else { "ALT" }
        ));
    }
    if sample.total() > min_count {
        sample
    } else {
        StrandTable {
            ref_fw: 0,
            ref_rv: 0,
            alt_fw: 0,
            alt_rv: 0,
            members: sample.members,
        }
    }
}

#[test]
fn forensic_6r165_java_fs_class_and_formula_contract() {
    // Pinned SHA 2dbc0258: FisherStrand extends StrandBiasTest.
    // annotate → calculateAnnotationFromLikelihoods → getContingencyTable(MIN_COUNT=2)
    // → pValueForContingencyTable → QualityUtils.phredScaleErrorRate → "%.3f"
    assert_eq!(JAVA_FS_MIN_COUNT, 2);
    assert_eq!(LOG_10_INFORMATIVE_THRESHOLD, 0.2);
    // Zero REF row: FisherExactTest lo>=hi → p=1 → FS=0. Not a rounding-to-zero artifact.
    let java_table = (0u32, 0u32, 2u32, 1u32);
    assert_eq!(java_table.0 + java_table.1, 0);
}

#[test]
fn forensic_6r165_java_sor_shares_table_not_min_count() {
    // StrandOddsRatio.calculateAnnotationFromLikelihoods → getContingencyTable(MIN_COUNT=0)
    // then calculateSOR with pseudocount 1, ln, "%.3f".
    assert_eq!(JAVA_SOR_MIN_COUNT, 0);
    let sor = java_calculate_sor(0, 0, 2, 1);
    assert!(
        (sor - 1.179).abs() < 0.001,
        "Java table [0,0;2,1] must print SOR=1.179, got {sor}"
    );
    assert!(
        (java_calculate_sor(0, 0, 3, 2) - 1.179).abs() > 0.05,
        "a different ALT split must not accidentally match Java SOR"
    );
}

#[test]
fn forensic_6r165_java_mq_is_rms_of_all_likelihood_evidence() {
    // RMSMappingQuality.annotate → calculateRawData:
    // for each likelihoods.sampleEvidence read with MQ != 255: squareSum += mq*mq
    // MQ = sqrt(squareSum / n), "%.2f". Not informative-only. Not ALT-only. Not a mean.
    assert_eq!(MAPPING_QUALITY_UNAVAILABLE, 255);
    let (n, sum_sq, rms) = java_rms_mq(&[37, 37, 37, 60]);
    assert_eq!(n, 4);
    assert_eq!(sum_sq, 37 * 37 * 3 + 60 * 60);
    assert!((rms - (sum_sq as f64 / 4.0).sqrt()).abs() < 1e-12);
}

#[test]
fn forensic_6r165_mq_pileup_vs_likelihoods_documented() {
    let src = include_str!("../src/variant_site_hc_annotations.rs");
    assert!(
        src.contains("fn read_strand_evidence_at_site(")
            && src.contains("for rec in &region.reads")
            && src.contains("mq_sum += u64::from(rec.mapq())"),
        "6R.165 ALT-pileup walk remains as the documented pre-6R.167 MQ path"
    );
    assert!(
        src.contains("strand_bias_contingency_table")
            && src.contains("6R.166: FS/SOR use Java `getContingencyTable`"),
        "6R.166 moved FS/SOR off the pileup 2x2; 6R.165 still documents that split"
    );
    assert!(
        src.contains("rms_mapping_quality_sample_mapqs")
            && src.contains("6R.167: MQ membership is Java")
            && src.contains("6R.168: MQ aggregation is Java"),
        "6R.167 membership stays; 6R.168 owns RMS"
    );
    assert!(
        !src.contains("92316347"),
        "no locus-specific annotation exception"
    );
}

#[test]
fn forensic_6r165_info_annotation_divergence() {
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing P12 ref/BAM");
        return;
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
        .expect("ActiveFull covering target");

    let call_args = CallRegionArgs::strict_java();
    let outcome = HaplotypeCallerEngine::call_region(covering, &dict, &ref_fasta, &call_args)
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
        .expect("genotyped G/A");
    let fmt = &call.genotype.format;
    assert_eq!(
        fmt.ad_as_i32(),
        vec![0, 3],
        "6R.164 FORMAT AD must stay closed"
    );
    assert_eq!(
        fmt.pl_as_i32(),
        vec![135, 9, 0],
        "6R.164 FORMAT PL must stay closed"
    );

    let emitted =
        try_emit_call_region_variants(covering, &outcome, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .expect("emit");
    let rec = emitted
        .iter()
        .find(|r| {
            r.position == TARGET
                && r.reference == MERGED_REF
                && r.alternate.first().map(String::as_str) == Some(MERGED_ALT)
        })
        .expect("emitted G/A");
    let info_keys: BTreeSet<String> = rec.info.iter().map(|v| info_key(v).to_string()).collect();
    eprintln!("6R165\temitted_info\t{:?}", rec.info);
    eprintln!(
        "6R165\trust_only_info_tags\t{:?}",
        info_keys
            .iter()
            .filter(|k| !matches!(
                k.as_str(),
                "AC" | "AF"
                    | "AN"
                    | "DP"
                    | "ExcessHet"
                    | "FS"
                    | "MLEAC"
                    | "MLEAF"
                    | "MQ"
                    | "QD"
                    | "SOR"
            ))
            .collect::<Vec<_>>()
    );

    let (ann_fs, ann_mq, ann_sor, ann_qual) = {
        let mut fs = 0.0;
        let mut mq = 0.0;
        let mut sor = 0.0;
        for v in &rec.info {
            match v {
                InfoValue::Float(k, xs) if k == "FS" => fs = xs[0],
                InfoValue::Float(k, xs) if k == "MQ" => mq = xs[0],
                InfoValue::Float(k, xs) if k == "SOR" => sor = xs[0],
                _ => {}
            }
        }
        (fs, mq, sor, rec.quality.unwrap_or(0.0))
    };
    eprintln!(
        "6R165\trust_final\tFS={:.5} MQ={} SOR={:.5} QUAL={:.5}",
        ann_fs, ann_mq, ann_sor, ann_qual
    );

    let (rust_table, rust_alt_mapqs) =
        rust_pileup_strand_table(covering.reads.iter().map(|r| r.as_ref()));
    let rust_mq_mean = if rust_alt_mapqs.is_empty() {
        0.0
    } else {
        f64::from(rust_alt_mapqs.iter().map(|&m| u32::from(m)).sum::<u32>())
            / rust_alt_mapqs.len() as f64
    };
    eprintln!(
        "6R165\trust_pileup_table\t[{},{};{},{}] members={}",
        rust_table.ref_fw,
        rust_table.ref_rv,
        rust_table.alt_fw,
        rust_table.alt_rv,
        rust_table.members.len()
    );
    for m in &rust_table.members {
        eprintln!("6R165\trust_pileup_read\t{m}");
    }
    eprintln!(
        "6R165\trust_alt_mapqs\t{:?} mean={}",
        rust_alt_mapqs, rust_mq_mean
    );

    let haps = &outcome.assembly.haplotypes;
    let (full_ref, full_pad) = outcome.assembly.event_map_reference();
    let apply_pad = haps
        .iter()
        .find(|h| h.is_reference)
        .and_then(|h| h.genome_loc.as_ref().map(|g| g.start.get()))
        .unwrap_or(full_pad);
    let hap_cache = build_per_haplotype_variation_events(
        haps,
        full_ref,
        full_pad,
        outcome.assembly.max_mnp_distance(),
        "2",
    );
    let mapping = create_allele_mapper_with_events(
        &VariationEvent::from_alleles("2", TARGET, MERGED_REF, MERGED_ALT),
        TARGET,
        haps,
        apply_pad,
        outcome.assembly.reference_bases(),
        outcome.assembly.max_mnp_distance(),
        true,
        Some(&hap_cache),
    );
    let ref_haps: Vec<usize> = mapping
        .ref_haplotype_indices
        .iter()
        .map(|i| i.get())
        .collect();
    let alt_haps: Vec<usize> = mapping
        .alt_haplotype_indices
        .iter()
        .map(|i| i.get())
        .collect();
    let geno: Vec<Record> = outcome
        .genotyping_reads
        .iter()
        .map(|r| (**r).clone())
        .collect();
    let rows = region_likelihoods_to_rows(&outcome.read_likelihoods, haps.len());
    let java_fs_table =
        java_informative_strand_table(&geno, &rows, &ref_haps, &alt_haps, JAVA_FS_MIN_COUNT);
    let java_sor_table =
        java_informative_strand_table(&geno, &rows, &ref_haps, &alt_haps, JAVA_SOR_MIN_COUNT);
    eprintln!(
        "6R165\tjava_fs_table\t[{},{};{},{}] members={}",
        java_fs_table.ref_fw,
        java_fs_table.ref_rv,
        java_fs_table.alt_fw,
        java_fs_table.alt_rv,
        java_fs_table.members.len()
    );
    for m in &java_fs_table.members {
        eprintln!("6R165\tjava_informative_read\t{m}");
    }
    eprintln!(
        "6R165\tjava_sor_table\t[{},{};{},{}]",
        java_sor_table.ref_fw, java_sor_table.ref_rv, java_sor_table.alt_fw, java_sor_table.alt_rv
    );

    let likelihood_mapqs: Vec<u8> = {
        let mut seen = BTreeSet::new();
        let mut mqs = Vec::new();
        for row in &rows {
            if !seen.insert(row.read_index) {
                continue;
            }
            if let Some(rec) = geno.get(row.read_index) {
                mqs.push(rec.mapq());
                eprintln!(
                    "6R165\tjava_mq_evidence\t{} strand={} dup={} sec={} sup={}",
                    read_id(rec),
                    strand_label(rec),
                    rec.flags() & FLAG_DUP != 0,
                    rec.flags() & FLAG_SECONDARY != 0,
                    rec.flags() & FLAG_SUPPLEMENTARY != 0
                );
            }
        }
        mqs
    };
    let (mq_n, mq_sum_sq, mq_rms) = java_rms_mq(&likelihood_mapqs);
    eprintln!(
        "6R165\tjava_mq_agg\tn={mq_n} sum_sq={mq_sum_sq} rms={mq_rms:.5} printed={:.2}",
        mq_rms
    );

    // FORMAT/QUAL closed. 6R.165 pileup vs likelihoods tables remain different.
    // 6R.166 production FS/SOR consume the likelihoods table.
    // 6R.168 production MQ is Java RMS 40.25 on 6R.167 sampleEvidence.
    assert!(
        (ann_mq - 40.25).abs() < 1e-12,
        "after 6R.168, production MQ is Java RMS 40.25, got {ann_mq}"
    );
    assert!(
        ann_fs < 0.02,
        "after 6R.166, production FS must be Java 0 from [0,0;2,1], got {ann_fs}"
    );
    assert!(
        (ann_sor - 1.179).abs() < 0.002,
        "after 6R.166, production SOR must be Java 1.179 from [0,0;2,1], got {ann_sor}"
    );

    assert_ne!(
        rust_table.cells(),
        java_fs_table.cells(),
        "FS/SOR first arrow is the 2x2 table, not the Fisher/SOR formula"
    );
    assert_eq!(
        java_fs_table.cells(),
        java_sor_table.cells(),
        "at this site FS MIN_COUNT=2 and SOR MIN_COUNT=0 admit the same 3 informative reads"
    );
    assert_eq!(
        rust_table.cells(),
        (0, 2, 3, 2),
        "Rust region.reads pileup table must be [0,2;3,2]"
    );
    assert_eq!(
        java_fs_table.cells(),
        (0, 0, 2, 1),
        "Java informative best-allele strand table must be [0,0;2,1] (AD=0,3; FS=0; SOR=1.179)"
    );
    assert_eq!(rust_alt_mapqs, vec![47, 47, 46, 21, 21]);
    assert_eq!(likelihood_mapqs, vec![47, 47, 21]);
    let java_sor = java_calculate_sor(
        java_fs_table.ref_fw,
        java_fs_table.ref_rv,
        java_fs_table.alt_fw,
        java_fs_table.alt_rv,
    );
    assert!(
        (java_sor - 1.179).abs() < 0.001,
        "Java reconstructed SOR from [0,0;2,1] must be 1.179, got {java_sor}"
    );
    assert!(
        (java_calculate_sor(
            rust_table.ref_fw,
            rust_table.ref_rv,
            rust_table.alt_fw,
            rust_table.alt_rv
        ) - 0.636)
            .abs()
            < 0.002,
        "6R.165 pileup table [0,2;3,2] still yields SOR≈0.636"
    );
    assert!(
        (ann_sor - java_sor).abs() < 1e-9,
        "production SOR must match Java formula on the likelihoods table"
    );

    assert!(
        (rust_mq_mean - 36.4).abs() < 1e-9,
        "6R.165 ALT-pileup reconstruction must remain mean 36.4, got {rust_mq_mean}",
    );
    assert!(
        (ann_mq - rust_mq_mean).abs() > 1.0,
        "production MQ must no longer equal the ALT-pileup mean ({ann_mq} vs {rust_mq_mean})",
    );
    let rust_mq_ids: BTreeSet<String> = rust_table
        .members
        .iter()
        .filter(|m| m.contains("allele=ALT"))
        .cloned()
        .collect();
    let java_mq_ids: BTreeSet<String> = {
        let mut seen = BTreeSet::new();
        let mut ids = BTreeSet::new();
        for row in &rows {
            if !seen.insert(row.read_index) {
                continue;
            }
            if let Some(rec) = geno.get(row.read_index) {
                if rec.mapq() != MAPPING_QUALITY_UNAVAILABLE {
                    ids.insert(read_id(rec));
                }
            }
        }
        ids
    };
    assert_ne!(
        rust_alt_mapqs.len(),
        mq_n,
        "MQ first arrow is evidence membership (ALT pileup vs all likelihood evidence)"
    );
    assert!(
        rust_mq_ids.len() != java_mq_ids.len()
            || rust_alt_mapqs.iter().copied().map(u64::from).sum::<u64>() * mq_n as u64
                != mq_sum_sq,
        "MQ membership and/or aggregation must differ"
    );
    assert!(
        (mq_rms - 40.25).abs() < 0.005,
        "Java RMS over likelihood evidence must print MQ=40.25, got {mq_rms}"
    );

    // Extra INFO tags are recorded, not investigated. 6R.178 omits InbreedingCoeff at n=1.
    assert!(
        info_keys.contains("ReadPosRankSum"),
        "ReadPosRankSum remains out of 6R.165 scope: {info_keys:?}"
    );
    assert!(
        !info_keys.contains("InbreedingCoeff"),
        "6R.178: one-sample records omit InbreedingCoeff, got {info_keys:?}"
    );
}
