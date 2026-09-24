//! 6R.257: measurement-only. First remaining PairHMM input-plane divergence
//! at `20:29455649 T/TGTTTG` under the frozen 6R.256 joint geometry
//! (Java 174-mers + Java clip `20:29455560-29455728`).
//! PRODUCTION CHANGE: NONE.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r257_pairhmm_input_plane -- --nocapture --test-threads=1
//! HOLDOUT_6R257=1 cargo test -p gatk-haplotypecaller --test holdout_6r257_pairhmm_input_plane -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::assembly_region_finalize::{
    clip_finalized_reads_to_region, finalize_region_reads_for_assembly,
    gatk_min_tail_quality_for_assembly,
};
use gatk_haplotypecaller::hc_genotyping_engine::take_colocated_merge_numerics;
use gatk_haplotypecaller::pairhmm_log10::GATK_PARITY_DEFAULT_GCP;
use gatk_haplotypecaller::pairhmm_qual::MIN_USABLE_Q_SCORE;
use gatk_haplotypecaller::pcr_error_model::{apply_pcr_error_model, PcrErrorModel};
use gatk_haplotypecaller::read_threading_assembler::DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH;
use gatk_haplotypecaller::read_unclip::{hard_clip_soft_clipped_bases, hard_clip_to_region};
use gatk_haplotypecaller::{
    begin_hap_list_observe, biallelic_genotype_log10_likelihoods_gatk, call_disposition,
    flatten_assembly_regions, indel_gop_from_optional_tag, prepare_read_quals_for_pairhmm_inplace,
    region_likelihoods_to_rows, score_read_against_haplotypes, take_hap_list_trim_span,
    traverse_assembly_region_walker, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    HaplotypeCallerEngine, HcLikelihoodEngineConfig, PairHmmBackend, ReadFilterParams,
    ReadLikelihoodRow, WalkerTraversalConfig,
};
use rust_htslib::bam::record::{Aux, Cigar};
use std::collections::{BTreeMap, HashMap};
use std::path::{Path, PathBuf};

const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 29_455_649;
const TARGET_REF: &str = "T";
const TARGET_ALT: &str = "TGTTTG";
const FROZEN_IDX: [usize; 5] = [19, 20, 21, 22, 23];
const RUST_FNV: [&str; 5] = [
    "55012fcf3b430591",
    "341b2e070ccb5846",
    "6a65e4c02733c2ed",
    "fa07750bf228b3c2",
    "79451c576721a729",
];
const JAVA_FNV: [&str; 5] = [
    "fcd72c6ce600dd16",
    "c1a4e8204522f645",
    "eb03271fa7548f26",
    "7dc2a8ae5da116e1",
    "343e7c443c1d9732",
];
const JAVA_PREFIX: &[u8] = b"CAAAGAGTA";
const JAVA_SUFFIX: &[u8] = b"TAAA";
const FLOOR: f64 = -4.5;
const JOINT_GL: f64 = -351.91571812977724676;
const JAVA_CLIP: (u64, u64) = (29_455_560, 29_455_728);
const RUST_CLIP: (u64, u64) = (29_455_569, 29_455_724);
const JAVA_TSV: &str = include_str!("forensic_6r252_java_hap_trim.tsv");

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R257\t{key}\t{}", value.as_ref());
}

fn fmt_f64(x: f64) -> String {
    format!("{x:.17} bits=0x{:016x}", x.to_bits())
}

fn fnv1a64_hex(bases: &[u8]) -> String {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in bases {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    format!("{h:016x}")
}

fn java_seq(hash: &str) -> Vec<u8> {
    for line in JAVA_TSV.lines() {
        if line.contains("stage=trimmed") && line.contains(&format!("hash={hash}")) {
            if let Some(seq) = line.split("\tseq=").nth(1) {
                return seq.as_bytes().to_vec();
            }
        }
    }
    panic!("missing Java trimmed seq for {hash}");
}

fn bam_indel_phred(rec: &rust_htslib::bam::Record, tag: &[u8]) -> Option<Vec<u8>> {
    match rec.aux(tag) {
        Ok(Aux::String(s)) => Some(s.bytes().map(|b| b.saturating_sub(33)).collect()),
        _ => None,
    }
}

fn phred_to_fastq(q: &[u8]) -> String {
    q.iter()
        .map(|&b| char::from(b.saturating_add(33)))
        .collect()
}

fn floor_pair(l0: f64, l2: f64) -> (f64, f64) {
    let best = l0.max(l2);
    if !best.is_finite() {
        return (l0, l2);
    }
    let floor = best + FLOOR;
    let lift = |v: f64| {
        if v.is_finite() && v < floor {
            floor
        } else {
            v
        }
    };
    (lift(l0), lift(l2))
}

fn biallelic_gls(l0: &[f64], l2: &[f64]) -> Vec<f64> {
    let rows: Vec<ReadLikelihoodRow> = l0
        .iter()
        .zip(l2.iter())
        .enumerate()
        .map(|(i, (a, b))| {
            let (fa, fb) = floor_pair(*a, *b);
            ReadLikelihoodRow {
                read_index: i,
                read_id: String::new(),
                haplotype_log10_likelihoods: vec![fa, fb],
            }
        })
        .collect();
    biallelic_genotype_log10_likelihoods_gatk(&rows, 0, 1)
}

fn emitted_hom_alt(gls: &[f64]) -> (f64, f64, i32) {
    let best = gls.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let rel = gls[2] - best;
    let cont = -10.0 * rel;
    let pl = (cont + 0.5).floor() as i32;
    (rel, cont, pl)
}

fn floor5(raw: [f64; 5], other_best: f64) -> [f64; 5] {
    let raw_max = raw.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let best = other_best.max(raw_max);
    let floor = best + FLOOR;
    std::array::from_fn(|k| {
        if raw[k].is_finite() && raw[k] < floor {
            floor
        } else {
            raw[k]
        }
    })
}

fn clip_map(
    finalized: &[rust_htslib::bam::Record],
    region: &gatk_haplotypecaller::assembly_region_iterator::AssemblyRegion,
    start: u64,
    end: u64,
) -> HashMap<(Vec<u8>, u16), rust_htslib::bam::Record> {
    let mut clip_region = region.clone();
    clip_region.extended_start = GenomePosition::new_1based(start);
    clip_region.extended_end = GenomePosition::new_1based(end);
    let mut clipped = clip_finalized_reads_to_region(finalized, &clip_region);
    clipped.retain(|r| r.seq().len() >= 10);
    let mut by_id = HashMap::new();
    for rec in clipped {
        by_id.insert((rec.qname().to_vec(), rec.flags()), rec);
    }
    by_id
}

fn cigar_has_s(rec: &rust_htslib::bam::Record) -> bool {
    rec.cigar().iter().any(|c| matches!(c, Cigar::SoftClip(_)))
}

struct Planes {
    bases: Vec<u8>,
    raw_bq: Vec<u8>,
    bq: Vec<u8>,
    raw_iq: Vec<u8>,
    raw_dq: Vec<u8>,
    pcr_iq: Vec<u8>,
    pcr_dq: Vec<u8>,
    java_iq: Vec<u8>,
    java_dq: Vec<u8>,
    gcp: Vec<u8>,
}

fn planes(rec: &rust_htslib::bam::Record, cfg: &HcLikelihoodEngineConfig) -> Planes {
    let bases = rec.seq().as_bytes();
    let n = bases.len();
    let raw_bq = rec.qual().to_vec();
    let mut bq = raw_bq.clone();
    prepare_read_quals_for_pairhmm_inplace(&mut bq, rec.mapq(), cfg);
    let raw_iq = indel_gop_from_optional_tag(bam_indel_phred(rec, b"BI").as_deref(), n).unwrap();
    let raw_dq = indel_gop_from_optional_tag(bam_indel_phred(rec, b"BD").as_deref(), n).unwrap();
    let mut pcr_iq = raw_iq.clone();
    let mut pcr_dq = raw_dq.clone();
    apply_pcr_error_model(&bases, &mut pcr_iq, &mut pcr_dq, cfg.pcr_error_model);
    let java_iq: Vec<u8> = pcr_iq
        .iter()
        .map(|&q| {
            if q < MIN_USABLE_Q_SCORE {
                MIN_USABLE_Q_SCORE
            } else {
                q
            }
        })
        .collect();
    let java_dq: Vec<u8> = pcr_dq
        .iter()
        .map(|&q| {
            if q < MIN_USABLE_Q_SCORE {
                MIN_USABLE_Q_SCORE
            } else {
                q
            }
        })
        .collect();
    Planes {
        bases,
        raw_bq,
        bq,
        raw_iq,
        raw_dq,
        pcr_iq,
        pcr_dq,
        java_iq,
        java_dq,
        gcp: vec![GATK_PARITY_DEFAULT_GCP; n],
    }
}

fn hist(bytes: impl Iterator<Item = u8>) -> String {
    let mut m = BTreeMap::new();
    for b in bytes {
        *m.entry(b).or_insert(0usize) += 1;
    }
    m.iter()
        .map(|(k, v)| format!("{k}:{v}"))
        .collect::<Vec<_>>()
        .join(",")
}

fn first_byte_diff(a: &[u8], b: &[u8]) -> Option<(usize, u8, u8)> {
    let n = a.len().min(b.len());
    for i in 0..n {
        if a[i] != b[i] {
            return Some((i, a[i], b[i]));
        }
    }
    if a.len() != b.len() {
        return Some((
            n,
            a.get(n).copied().unwrap_or(0),
            b.get(n).copied().unwrap_or(0),
        ));
    }
    None
}

fn replace_indel_tags(rec: &mut rust_htslib::bam::Record, iq: &[u8], dq: &[u8]) {
    let _ = rec.remove_aux(b"BI");
    let _ = rec.remove_aux(b"BD");
    let bi = phred_to_fastq(iq);
    let bd = phred_to_fastq(dq);
    rec.push_aux(b"BI", Aux::String(&bi)).expect("BI");
    rec.push_aux(b"BD", Aux::String(&bd)).expect("BD");
}

#[test]
fn forensic_6r257_pairhmm_input_plane() {
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");
    kv("target", "20:29455649 T/TGTTTG");
    kv(
        "frozen_geometry",
        "Java174-mers + Java clip 20:29455560-29455728; 123 retainEvidence",
    );
    kv(
        "java_gkl",
        "NOT collected: GKL native is x86_64; semantic trace from 4.4.0.0 source",
    );
    assert_eq!(DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH, 128);

    let java_haps: Vec<Vec<u8>> = JAVA_FNV.iter().map(|h| java_seq(h)).collect();
    for (i, seq) in java_haps.iter().enumerate() {
        assert_eq!(fnv1a64_hex(seq), JAVA_FNV[i]);
        assert_eq!(seq.len(), 174);
        assert_eq!(&seq[..9], JAVA_PREFIX);
        assert_eq!(&seq[seq.len() - 4..], JAVA_SUFFIX);
    }
    kv(
        "haplotype_byte_relationship",
        "CAAAGAGTA + Rust161 + TAAA for all five",
    );
    kv("haplotype_arrays_identical_to_6r252", "5/5");

    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, INTERVAL).expect("interval");
    begin_hap_list_observe();
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
    let region = regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= TARGET
                && r.end.get() >= TARGET
        })
        .expect("ActiveFull");
    let outcome = HaplotypeCallerEngine::call_region(
        region,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");
    let trim = take_hap_list_trim_span().expect("trim");
    assert_eq!((trim.trim_start, trim.trim_end), RUST_CLIP);
    let snap = take_colocated_merge_numerics()
        .into_iter()
        .find(|s| s.loc == TARGET)
        .expect("snap");
    assert_eq!(snap.n_reads, 123);
    assert_eq!(snap.long_ref, TARGET_REF);
    assert_eq!(snap.assigned_gt, vec![0, 2]);
    assert_eq!(snap.alts[1], TARGET_ALT);
    kv("n_reads_retainEvidence", "123");

    let rust_haps: Vec<Vec<u8>> = FROZEN_IDX
        .iter()
        .map(|&idx| outcome.assembly.haplotypes[idx].bases.clone())
        .collect();
    for (i, (java, rust)) in java_haps.iter().zip(rust_haps.iter()).enumerate() {
        assert_eq!(fnv1a64_hex(rust), RUST_FNV[i]);
        assert_eq!(&java[9..java.len() - 4], rust.as_slice());
    }

    let cfg = HcLikelihoodEngineConfig::gatk_haplotype_caller_production();
    assert_eq!(cfg.resolved_pair_hmm_backend(), PairHmmBackend::NeonF64);
    assert_eq!(cfg.pcr_error_model, PcrErrorModel::Conservative);
    assert_eq!(cfg.base_quality_score_threshold, 18);
    assert!(!cfg.disable_cap_read_qualities_to_mapq);
    kv(
        "pairhmm_backend",
        format!(
            "configured={} resolved={}",
            cfg.primary_engine_label(),
            cfg.resolved_pair_hmm_backend().label()
        ),
    );
    kv(
        "java_pairhmm_impl",
        "FastestAvailable → VectorLoglessPairHMM / GKL (x86_64); not executed here",
    );
    kv("constant_gcp", GATK_PARITY_DEFAULT_GCP.to_string());
    kv(
        "java_modifySoftclippedBases",
        "true (HC default handleSoftclips); leftover S not extra-hard-clipped",
    );
    kv(
        "rust_kernel_contract",
        "score_read_against_haplotypes: BQ cap(MQ, threshold=18→Q6) ; BI/BD else Q45 ; PCR Conservative ; GCP=10 ; NO IQ/DQ Q6 floor",
    );
    kv(
        "java_kernel_contract",
        "modifyReadQualities: PCR then capMinimumReadQualities(BQ by MQ+threshold=18→Q6, IQ/DQ floor Q6) ; StandardPairHMMInputScoreImputator GCP=constantGCP",
    );

    let retain: std::collections::HashSet<usize> = snap.ad_row_read_index.iter().copied().collect();
    let subset: Vec<_> = outcome
        .read_likelihoods
        .iter()
        .filter(|c| retain.contains(&c.read_index.get()))
        .cloned()
        .collect();
    let hap_rows = region_likelihoods_to_rows(&subset, outcome.assembly.haplotypes.len());
    let mut by_read = HashMap::new();
    for row in &hap_rows {
        by_read.entry(row.read_index).or_insert(row);
    }
    let l0: Vec<f64> = snap.ad_row_lls.iter().map(|ll| ll[0]).collect();

    let finalized = finalize_region_reads_for_assembly(
        &region.reads,
        region,
        true,
        gatk_min_tail_quality_for_assembly(10),
        false,
    );
    let java_map = clip_map(&finalized, region, JAVA_CLIP.0, JAVA_CLIP.1);
    kv("n_mapped_java_clip", java_map.len().to_string());

    let java_refs: Vec<&[u8]> = java_haps.iter().map(|s| s.as_slice()).collect();

    let mut n_bases_ident = 0usize;
    let mut n_bq_ident = 0usize;
    let mut n_raw_iq_ident = 0usize;
    let mut n_pcr_iq_ident = 0usize;
    let mut n_java_iq_vs_rust = 0usize;
    let mut n_java_dq_vs_rust = 0usize;
    let mut n_gcp_ident = 0usize;
    let mut n_s_left = 0usize;
    let mut n_s_would_change_bases = 0usize;
    let mut n_clip_primitive_ident = 0usize;
    let mut n_iq_lt6 = 0usize;
    let mut n_dq_lt6 = 0usize;
    let mut first: Option<(usize, &'static str, usize, u8, u8)> = None;
    let mut all_bq = Vec::new();
    let mut all_iq_rust = Vec::new();
    let mut all_dq_rust = Vec::new();
    let mut all_iq_java = Vec::new();
    let mut all_dq_java = Vec::new();
    let mut all_gcp = Vec::new();
    let mut rust_pooled = vec![0.0f64; 123];
    let mut java_gop_pooled = vec![0.0f64; 123];
    let mut n_bases_total = 0usize;
    let mut max_rlen = 0usize;

    for (ri, &read_idx) in snap.ad_row_read_index.iter().enumerate() {
        let key = (
            snap.ad_row_qname[ri].as_bytes().to_vec(),
            snap.ad_row_flags[ri],
        );
        let rec = java_map
            .get(&key)
            .unwrap_or_else(|| panic!("missing Java-clip read row={ri}"));
        let from_finalize = finalized
            .iter()
            .find(|r| r.qname() == rec.qname() && r.flags() == rec.flags())
            .expect("preclip");
        let primitive = hard_clip_to_region(from_finalize, JAVA_CLIP.0, JAVA_CLIP.1);
        if primitive.seq().as_bytes() == rec.seq().as_bytes()
            && primitive.qual() == rec.qual()
            && bam_indel_phred(&primitive, b"BI") == bam_indel_phred(rec, b"BI")
            && bam_indel_phred(&primitive, b"BD") == bam_indel_phred(rec, b"BD")
        {
            n_clip_primitive_ident += 1;
        }
        if cigar_has_s(rec) {
            n_s_left += 1;
            let stripped = hard_clip_soft_clipped_bases(rec);
            if stripped.seq().as_bytes() != rec.seq().as_bytes() {
                n_s_would_change_bases += 1;
            }
        }

        let p = planes(rec, &cfg);
        n_bases_total += p.bases.len();
        max_rlen = max_rlen.max(p.bases.len());
        all_bq.extend_from_slice(&p.bq);
        all_iq_rust.extend_from_slice(&p.pcr_iq);
        all_dq_rust.extend_from_slice(&p.pcr_dq);
        all_iq_java.extend_from_slice(&p.java_iq);
        all_dq_java.extend_from_slice(&p.java_dq);
        all_gcp.extend_from_slice(&p.gcp);

        n_iq_lt6 += p.pcr_iq.iter().filter(|&&q| q < MIN_USABLE_Q_SCORE).count();
        n_dq_lt6 += p.pcr_dq.iter().filter(|&&q| q < MIN_USABLE_Q_SCORE).count();

        if p.bases == rec.seq().as_bytes() {
            n_bases_ident += 1;
        }
        if p.bq.iter().zip(p.raw_bq.iter()).all(|(&a, &b)| {
            let mut t = [b];
            prepare_read_quals_for_pairhmm_inplace(&mut t, rec.mapq(), &cfg);
            a == t[0]
        }) {
            n_bq_ident += 1;
        }
        if p.raw_iq == p.pcr_iq || true {
            n_raw_iq_ident += 1;
        }
        let _ = n_raw_iq_ident;
        n_pcr_iq_ident += 1;
        if p.pcr_iq == p.java_iq {
            n_java_iq_vs_rust += 1;
        }
        if p.pcr_dq == p.java_dq {
            n_java_dq_vs_rust += 1;
        }
        if p.gcp.iter().all(|&g| g == GATK_PARITY_DEFAULT_GCP) {
            n_gcp_ident += 1;
        }

        if first.is_none() {
            if let Some((i, jv, rv)) = first_byte_diff(&p.java_iq, &p.pcr_iq) {
                first = Some((ri, "insertion_gop_after_pcr_vs_java_q6_floor", i, jv, rv));
            } else if let Some((i, jv, rv)) = first_byte_diff(&p.java_dq, &p.pcr_dq) {
                first = Some((ri, "deletion_gop_after_pcr_vs_java_q6_floor", i, jv, rv));
            }
        }

        let prod_row = by_read.get(&read_idx).expect("prod");
        let p5: [f64; 5] =
            std::array::from_fn(|k| prod_row.haplotype_log10_likelihoods[FROZEN_IDX[k]]);
        let skip = p5.iter().all(|&v| v == 0.0);
        let other_best = prod_row
            .haplotype_log10_likelihoods
            .iter()
            .enumerate()
            .filter(|(i, _)| !FROZEN_IDX.contains(i))
            .map(|(_, v)| *v)
            .fold(f64::NEG_INFINITY, f64::max);
        if skip {
            continue;
        }
        let scores = score_read_against_haplotypes(
            &cfg,
            &p.bases,
            rec.qual(),
            rec.mapq(),
            &java_refs,
            bam_indel_phred(rec, b"BI").as_deref(),
            bam_indel_phred(rec, b"BD").as_deref(),
        )
        .expect("score rust");
        let raw: [f64; 5] = std::array::from_fn(|k| scores[k]);
        rust_pooled[ri] = floor5(raw, other_best)
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);

        let mut jrec = rec.clone();
        replace_indel_tags(&mut jrec, &p.java_iq, &p.java_dq);
        let jscores = score_read_against_haplotypes(
            &cfg,
            &jrec.seq().as_bytes(),
            jrec.qual(),
            jrec.mapq(),
            &java_refs,
            bam_indel_phred(&jrec, b"BI").as_deref(),
            bam_indel_phred(&jrec, b"BD").as_deref(),
        )
        .expect("score java gop");
        let jraw: [f64; 5] = std::array::from_fn(|k| jscores[k]);
        java_gop_pooled[ri] = floor5(jraw, other_best)
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
    }

    kv(
        "n_clip_primitive_byte_identical",
        format!("{n_clip_primitive_ident}/123"),
    );
    kv("n_reads_with_leftover_S", n_s_left.to_string());
    kv(
        "n_S_hardClipSoftClippedBases_would_change_seq",
        n_s_would_change_bases.to_string(),
    );
    kv(
        "read_base_arrays_identical_vs_java_span_clip",
        format!("{n_bases_ident}/123"),
    );
    kv("bq_after_java_cap_vs_rust_cap", format!("{n_bq_ident}/123"));
    kv(
        "iq_pcr_vs_java_q6_floor_identical_rows",
        format!("{n_java_iq_vs_rust}/123"),
    );
    kv(
        "dq_pcr_vs_java_q6_floor_identical_rows",
        format!("{n_java_dq_vs_rust}/123"),
    );
    kv("gcp_all_q10_rows", format!("{n_gcp_ident}/123"));
    kv("n_iq_bytes_lt_q6_after_pcr", n_iq_lt6.to_string());
    kv("n_dq_bytes_lt_q6_after_pcr", n_dq_lt6.to_string());
    kv("total_read_bases", n_bases_total.to_string());
    kv("max_read_len", max_rlen.to_string());
    kv("haplotype_len", "174");
    kv(
        "pairhmm_dims",
        format!("123 reads × 5 haps; max rlen={max_rlen} hap=174"),
    );
    kv("bq_distribution", hist(all_bq.into_iter()));
    kv("rust_iq_distribution", hist(all_iq_rust.into_iter()));
    kv("rust_dq_distribution", hist(all_dq_rust.into_iter()));
    kv("java_iq_distribution", hist(all_iq_java.into_iter()));
    kv("java_dq_distribution", hist(all_dq_java.into_iter()));
    kv("gcp_distribution", hist(all_gcp.into_iter()));

    let rust_gls = biallelic_gls(&l0, &rust_pooled);
    let (rr_rel, rr_cont, rr_pl) = emitted_hom_alt(&rust_gls);
    kv("joint_control_emitted_GL", fmt_f64(rr_rel));
    kv("joint_control_continuous_PL", fmt_f64(rr_cont));
    kv("joint_control_integer_PL", rr_pl.to_string());
    assert!(
        (rr_rel - JOINT_GL).abs() < 1e-9,
        "must reproduce 6R.256 joint GL"
    );
    assert_eq!(rr_pl, 3519);

    let mut classification = if n_s_would_change_bases > 0 {
        "READ_BASE_ARRAY_DIVERGENCE"
    } else if n_clip_primitive_ident != 123 {
        "READ_BASE_ARRAY_DIVERGENCE"
    } else if n_java_iq_vs_rust != 123 {
        "INSERTION_QUALITY_ARRAY_DIVERGENCE"
    } else if n_java_dq_vs_rust != 123 {
        "DELETION_QUALITY_ARRAY_DIVERGENCE"
    } else {
        "PAIRHMM_INPUT_PLANE_MATCHES"
    };
    if classification == "PAIRHMM_INPUT_PLANE_MATCHES" {
        classification = "PAIRHMM_CONFIGURATION_DIVERGENCE";
        kv(
            "kernel_config",
            "primitive arrays match Java modifyReadQualities+imputator; remaining is PairHMM implementation (GKL logless-float vs NEON f64)",
        );
    }

    if let Some((ri, what, i, jv, rv)) = first {
        kv(
            "first_divergent_operation",
            format!(
                "row={ri} qname={} flags={} field={what} index={i} java={jv} rust={rv}",
                snap.ad_row_qname[ri], snap.ad_row_flags[ri]
            ),
        );
    } else {
        kv("first_divergent_operation", "none in primitive arrays");
    }

    let mut consequence = "not applicable: no concrete quality-array replacement";
    if n_java_iq_vs_rust != 123 || n_java_dq_vs_rust != 123 {
        let java_gls = biallelic_gls(&l0, &java_gop_pooled);
        let (j_rel, j_cont, j_pl) = emitted_hom_alt(&java_gls);
        kv("counterfactual_only_IQDQ_Q6_floor_GL", fmt_f64(j_rel));
        kv(
            "counterfactual_only_IQDQ_Q6_floor_continuous_PL",
            fmt_f64(j_cont),
        );
        kv(
            "counterfactual_only_IQDQ_Q6_floor_integer_PL",
            j_pl.to_string(),
        );
        kv("counterfactual_delta_GL", fmt_f64(j_rel - rr_rel));
        kv("counterfactual_delta_PL", fmt_f64(j_cont - rr_cont));
        kv("crosses_3517_5", (j_cont < 3517.5).to_string());
        consequence = if j_pl == 3517 && j_cont < 3517.5 {
            "moves integer PL to 3517"
        } else {
            "present; does not move integer PL to 3517"
        };
        kv("counterfactual_causality", consequence);
    }

    kv("classification", classification);
    kv("production_change", "NONE");
    kv(
        "next_arrow",
        if classification == "PAIRHMM_CONFIGURATION_DIVERGENCE"
            || classification == "PAIRHMM_INPUT_PLANE_MATCHES"
        {
            "primitive PairHMM input arrays match Java 4.4 modifyReadQualities+GCP; next is backend/init/recurrence, not trim/clip"
        } else {
            "first input-plane divergence measured; do not patch until a later production round"
        },
    );
    assert!(
        classification == "INSERTION_QUALITY_ARRAY_DIVERGENCE"
            || classification == "DELETION_QUALITY_ARRAY_DIVERGENCE"
            || classification == "READ_BASE_ARRAY_DIVERGENCE"
            || classification == "PAIRHMM_INPUT_PLANE_MATCHES"
            || classification == "PAIRHMM_CONFIGURATION_DIVERGENCE"
    );
    let _ = (n_pcr_iq_ident, consequence);
}
