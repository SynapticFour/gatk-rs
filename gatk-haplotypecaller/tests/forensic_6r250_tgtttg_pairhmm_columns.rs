//! 6R.250: proof-only. Trace the five TGTTTG PairHMM columns at
//! `20:29455649 T/TGTTTG`.
//!
//! 6R.249 closed pooling: the allele-row is `max` of columns 19–23.
//! This round audits PairHMM inputs/outputs for those five haplotypes.
//! PRODUCTION CHANGE: NONE.
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r250_tgtttg_pairhmm_columns -- --nocapture --test-threads=1
//! HOLDOUT_6R250=1 cargo test -p gatk-haplotypecaller --test holdout_6r250_tgtttg_pairhmm_columns -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::{
    build_per_haplotype_variation_events, overlapping_events, remap_alt_onto_longer_ref,
};
use gatk_haplotypecaller::hc_genotyping_engine::take_colocated_merge_numerics;
use gatk_haplotypecaller::pcr_error_model::apply_pcr_error_model;
use gatk_haplotypecaller::read_threading_assembler::DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH;
use gatk_haplotypecaller::{
    begin_hap_list_observe, begin_likelihood_pipeline_observe, call_disposition,
    flatten_assembly_regions, indel_gop_from_optional_tag, prepare_read_quals_for_pairhmm_inplace,
    region_likelihoods_to_rows, take_hap_list_snaps, take_hap_list_trim_span,
    take_likelihood_pipeline_cells, traverse_assembly_region_walker, AssemblyRegionCallDisposition,
    CallRegionArgs, GenomePosition, Haplotype, HaplotypeCallerEngine, HcLikelihoodEngineConfig,
    PairHmmBackend, ReadFilterParams, WalkerTraversalConfig, GATK_PARITY_DEFAULT_GCP,
};
use rust_htslib::bam::record::Aux;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs;
use std::path::{Path, PathBuf};

const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 29_455_649;
const TARGET_REF: &str = "T";
const TARGET_ALT: &str = "TGTTTG";
const UNUSED_ALT: &str = "TTTG";
const FROZEN_IDX: [usize; 5] = [19, 20, 21, 22, 23];
const FROZEN_FNV: [&str; 5] = [
    "55012fcf3b430591",
    "341b2e070ccb5846",
    "6a65e4c02733c2ed",
    "fa07750bf228b3c2",
    "79451c576721a729",
];
const REP_QNAME: &str = "HISEQ1:12:H8GVUADXX:1:2107:7301:87830";
const REP_FLAGS: u16 = 83;
const REP_WINNER_IDX: usize = 22;
const REP_WINNER_LL: f64 = -5.02439935627938894;
const RUST_HOM_ALT_GL: f64 = -351.75162674083281900;
const KNOWN_ALIGNED_RESIDUAL: f64 = 7.41e-6;
const EMPTY_POOL: f64 = -50.0;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R250\t{key}\t{}", value.as_ref());
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

fn cigar_str(h: &Haplotype) -> String {
    h.cigar
        .as_ref()
        .map(|c| {
            c.elements
                .iter()
                .map(|e| format!("{}{}", e.length, e.operator.as_char()))
                .collect::<String>()
        })
        .unwrap_or_else(|| ".".to_string())
}

fn fn_body<'a>(src: &'a str, sig: &str) -> &'a str {
    let start = src.find(sig).unwrap_or(0);
    let rest = &src[start..];
    let rel = rest[sig.len()..]
        .find("\nfn ")
        .map(|i| sig.len() + i)
        .unwrap_or(rest.len());
    &rest[..rel]
}

fn java_mapper_pools(
    n_haps: usize,
    hap_events: &gatk_haplotypecaller::event_map::PerHaplotypeVariationEvents,
    loc: u64,
    long_ref: &str,
    alts: &[String],
) -> Vec<Vec<usize>> {
    let loc_pos = GenomePosition::new_1based(loc);
    let mut pools: Vec<Vec<usize>> = vec![Vec::new(); 1 + alts.len()];
    for i in 0..n_haps {
        let spanning = overlapping_events(hap_events.events_for(i), loc);
        if spanning.is_empty() {
            pools[0].push(i);
            continue;
        }
        let mut in_alt = vec![false; alts.len()];
        for ev in &spanning {
            if ev.start_1based == loc_pos {
                if ev.ref_allele.len() == long_ref.len() {
                    if let Some(ai) = alts.iter().position(|a| a == &ev.alt_allele) {
                        in_alt[ai] = true;
                    }
                } else if ev.ref_allele.len() < long_ref.len() {
                    if let Some(remapped) =
                        remap_alt_onto_longer_ref(&ev.ref_allele, &ev.alt_allele, long_ref)
                    {
                        if let Some(ai) = alts.iter().position(|a| a == &remapped) {
                            in_alt[ai] = true;
                        }
                    }
                }
            } else {
                break;
            }
        }
        if !in_alt.iter().any(|&b| b) {
            continue;
        }
        for (ai, hit) in in_alt.iter().enumerate() {
            if *hit {
                pools[ai + 1].push(i);
            }
        }
    }
    pools
}

fn pool_max(indices: &[usize], lls: &[f64]) -> f64 {
    let ll = indices
        .iter()
        .filter_map(|&i| lls.get(i).copied())
        .fold(f64::NEG_INFINITY, f64::max);
    if ll.is_finite() {
        ll
    } else {
        EMPTY_POOL
    }
}

fn bam_indel_phred(rec: &rust_htslib::bam::Record, tag: &[u8]) -> Option<Vec<u8>> {
    match rec.aux(tag) {
        Ok(Aux::String(s)) => Some(s.bytes().map(|b| b.saturating_sub(33)).collect()),
        _ => None,
    }
}

fn pairhmm_planes(
    rec: &rust_htslib::bam::Record,
    cfg: &HcLikelihoodEngineConfig,
) -> (Vec<u8>, Vec<u8>, Vec<u8>, Vec<u8>, Vec<u8>) {
    let bases = rec.seq().as_bytes();
    let mapq = rec.mapq();
    let mut bq = rec.qual().to_vec();
    prepare_read_quals_for_pairhmm_inplace(&mut bq, mapq, cfg);
    let mut ins = indel_gop_from_optional_tag(bam_indel_phred(rec, b"BI").as_deref(), bases.len())
        .expect("BI");
    let mut del = indel_gop_from_optional_tag(bam_indel_phred(rec, b"BD").as_deref(), bases.len())
        .expect("BD");
    apply_pcr_error_model(&bases, &mut ins, &mut del, cfg.pcr_error_model);
    let gcp = vec![GATK_PARITY_DEFAULT_GCP; bases.len()];
    (bases, bq, ins, del, gcp)
}

#[test]
fn forensic_6r250_source_pairhmm_lifecycle() {
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");
    assert_eq!(DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH, 128);

    let ll_eng =
        fs::read_to_string(repo_root().join("gatk-haplotypecaller/src/likelihood_engine.rs"))
            .expect("likelihood_engine.rs");
    let score = fn_body(&ll_eng, "pub fn score_read_against_haplotypes");
    assert!(score.contains("prepare_read_quals_for_pairhmm_inplace"));
    assert!(score.contains("apply_pcr_error_model"));
    assert!(score.contains("fill_indel_gop_from_optional_tag"));
    assert!(score.contains("scratch.gcp[..n].fill(GATK_PARITY_DEFAULT_GCP)"));
    assert!(
        !score.contains("cigar"),
        "PairHMM kernel takes haplotype bases, not CIGAR"
    );
    kv(
        "rust_pairhmm_prep",
        "per-read: BQ cap(threshold=18, MAPQ) → GOP from BI/BD or Q45 → PCR conservative on IQ/DQ → GCP=10; then score all hap bases with that single plane",
    );

    let eng_ll =
        fs::read_to_string(repo_root().join("gatk-haplotypecaller/src/engine_likelihoods.rs"))
            .expect("engine_likelihoods.rs");
    assert!(eng_ll.contains("clip_finalized_reads_in_place"));
    assert!(eng_ll.contains("score_pairhmm_from_records_java_mate_contig"));
    kv(
        "rust_read_lifecycle",
        "assemble finalize → hardClipToRegion on trim span → unclipped-length filter → mate-contig skip (0 cells) → PairHMM",
    );
    kv(
        "java_lifecycle",
        "HaplotypeCallerEngine.callRegion: finalizeRegion + hardClipToRegion → PairHMMLikelihoodCalculationEngine.computeReadLikelihoods → pairHMM.computeLog10Likelihoods → normalizeLikelihoods(-4.5) → filterPoorlyModeledEvidence → createAlleleMapper → marginalize",
    );
    let cfg = HcLikelihoodEngineConfig::gatk_haplotype_caller_production();
    kv(
        "pairhmm_backend",
        format!(
            "configured={} resolved={}",
            cfg.primary_engine_label(),
            cfg.resolved_pair_hmm_backend().label()
        ),
    );
    assert_eq!(
        cfg.resolved_pair_hmm_backend(),
        PairHmmBackend::NeonF64,
        "this host is aarch64 NEON f64"
    );
}

#[test]
fn forensic_6r250_live_five_pairhmm_columns() {
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");
    kv("target", "20:29455649 T/TGTTTG");
    kv(
        "direct_java_executable",
        "NO: gatk_image_absent; not pulled",
    );

    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, INTERVAL).expect("interval");
    begin_hap_list_observe();
    begin_likelihood_pipeline_observe();
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
    let hap_snaps = take_hap_list_snaps();
    let trim = take_hap_list_trim_span();
    let pipe_cells = take_likelihood_pipeline_cells();
    let snap = take_colocated_merge_numerics()
        .into_iter()
        .find(|s| s.loc == TARGET)
        .expect("merge snapshot");

    assert_eq!(snap.n_reads, 123);
    assert_eq!(snap.long_ref, TARGET_REF);
    assert_eq!(snap.assigned_gt, vec![0, 2]);
    assert_eq!(snap.alts, [UNUSED_ALT, TARGET_ALT]);
    kv("n_reads_retainEvidence", "123");
    kv("n_pairhmm_region", snap.n_pairhmm_reads.to_string());
    kv("n_haps_genotyping", snap.n_haps.to_string());

    let haps = &outcome.assembly.haplotypes;
    for (i, &idx) in FROZEN_IDX.iter().enumerate() {
        let h = &haps[idx];
        let hash = fnv1a64_hex(&h.bases);
        assert_eq!(hash, FROZEN_FNV[i], "frozen FNV H{i}");
        assert_eq!(h.bases.len(), 161);
        assert!(!h.is_reference);
        kv(
            &format!("hap_H{i}"),
            format!(
                "idx={idx} fnv={hash} is_ref=false len=161 kmer={} score={:.6} cigar={} align_start={} loc={}-{} pairhmm_input=bases_only",
                h.kmer_size,
                h.score,
                cigar_str(h),
                h.alignment_start_hap_wrt_ref,
                h.genome_loc.map(|g| g.start_1based()).unwrap_or(0),
                h.genome_loc.map(|g| g.end_1based()).unwrap_or(0),
            ),
        );
        kv(
            &format!("hap_H{i}_bases"),
            String::from_utf8_lossy(&h.bases).into_owned(),
        );
    }

    let stages: Vec<String> = hap_snaps
        .iter()
        .map(|s| format!("{}:n={}", s.stage, s.n))
        .collect();
    kv("hap_list_stages", stages.join(" | "));
    if let Some(t) = &trim {
        kv(
            "trim_span",
            format!(
                "active={}-{} trim={}-{}",
                t.active_start, t.active_end, t.trim_start, t.trim_end
            ),
        );
    }
    let pairhmm_input = hap_snaps
        .iter()
        .find(|s| s.stage == "pairhmm_input")
        .expect("pairhmm_input stage");
    for (i, &idx) in FROZEN_IDX.iter().enumerate() {
        let col = pairhmm_input
            .columns
            .iter()
            .find(|c| c.index == idx)
            .expect("idx at pairhmm_input");
        assert_eq!(format!("{:016x}", col.fnv1a), FROZEN_FNV[i]);
        assert_eq!(col.len, 161);
        assert_eq!(col.cigar, "81M5I75M");
    }
    kv(
        "hap_identity_pairhmm_vs_genotype",
        "YES: idx 19-23 FNV/len/CIGAR identical at pairhmm_input and genotyping",
    );

    let hap_events = build_per_haplotype_variation_events(
        haps,
        outcome.assembly.reference_bases(),
        outcome.assembly.padded_reference_start_1based(),
        outcome.assembly.max_mnp_distance(),
        "20",
    );
    let alts = vec![UNUSED_ALT.to_string(), TARGET_ALT.to_string()];
    let pools = java_mapper_pools(haps.len(), &hap_events, TARGET, "T", &alts);
    assert_eq!(pools[2], FROZEN_IDX);

    let cfg = HcLikelihoodEngineConfig::gatk_haplotype_caller_production();
    kv(
        "pairhmm_mode",
        format!(
            "Java GKL float / Rust {}",
            cfg.resolved_pair_hmm_backend().label()
        ),
    );
    kv("pcr_error_model", format!("{:?}", cfg.pcr_error_model));
    kv("bq_threshold", cfg.base_quality_score_threshold.to_string());
    kv("gcp", GATK_PARITY_DEFAULT_GCP.to_string());

    let retain: HashSet<usize> = snap.ad_row_read_index.iter().copied().collect();
    assert_eq!(retain.len(), 123);
    let subset: Vec<_> = outcome
        .read_likelihoods
        .iter()
        .filter(|c| retain.contains(&c.read_index.get()))
        .cloned()
        .collect();
    let hap_rows = region_likelihoods_to_rows(&subset, haps.len());
    let mut by_read = HashMap::new();
    for row in &hap_rows {
        by_read.entry(row.read_index).or_insert(row);
    }

    let mut col_min = [f64::INFINITY; 5];
    let mut col_max = [f64::NEG_INFINITY; 5];
    let mut col_sum = [0.0f64; 5];
    let mut n_unique = 0usize;
    let mut n_ties = 0usize;
    let mut tie_sets: BTreeMap<u8, usize> = BTreeMap::new();
    let mut winner_only = [0usize; 5];
    let mut max_spread = 0.0f64;
    let mut n_identical_five = 0usize;
    let mut first_rep: Option<usize> = None;

    for (ri, &read_idx) in snap.ad_row_read_index.iter().enumerate() {
        let row = by_read.get(&read_idx).expect("row");
        let five: [f64; 5] =
            std::array::from_fn(|k| row.haplotype_log10_likelihoods[FROZEN_IDX[k]]);
        let mx = five.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let mn = five.iter().copied().fold(f64::INFINITY, f64::min);
        max_spread = max_spread.max(mx - mn);
        let mut mask = 0u8;
        for (k, &v) in five.iter().enumerate() {
            col_min[k] = col_min[k].min(v);
            col_max[k] = col_max[k].max(v);
            col_sum[k] += v;
            if v == mx {
                mask |= 1 << k;
            }
        }
        *tie_sets.entry(mask).or_insert(0) += 1;
        let n_at = mask.count_ones() as usize;
        if n_at == 1 {
            n_unique += 1;
            winner_only[mask.trailing_zeros() as usize] += 1;
        } else {
            n_ties += 1;
        }
        if five.windows(2).all(|w| w[0] == w[1]) {
            n_identical_five += 1;
        }
        let qname = snap.ad_row_qname.get(ri).map(String::as_str).unwrap_or("");
        let flags = snap.ad_row_flags.get(ri).copied().unwrap_or(0);
        if qname == REP_QNAME && flags == REP_FLAGS {
            first_rep = Some(ri);
            kv(
                "rep_row_five_ll",
                format!(
                    "row={ri} read_index={read_idx} five=[{}] max={} winner_idx={}",
                    five.iter()
                        .map(|v| fmt_f64(*v))
                        .collect::<Vec<_>>()
                        .join(", "),
                    fmt_f64(mx),
                    FROZEN_IDX[five
                        .iter()
                        .enumerate()
                        .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
                        .map(|x| x.0)
                        .unwrap_or(0)],
                ),
            );
            assert!((mx - REP_WINNER_LL).abs() < 1e-12);
            assert_eq!(
                FROZEN_IDX[five
                    .iter()
                    .enumerate()
                    .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
                    .map(|x| x.0)
                    .unwrap()],
                REP_WINNER_IDX
            );
        }
        let pooled = pool_max(&pools[2], &row.haplotype_log10_likelihoods);
        assert_eq!(pooled, mx);
        assert!((pooled - snap.ad_row_lls[ri][2]).abs() <= 1e-15);
    }
    assert!(first_rep.is_some(), "representative read in 123-row object");
    kv("n_unique_max", n_unique.to_string());
    kv("n_ties", n_ties.to_string());
    kv("n_identical_five", n_identical_five.to_string());
    kv("max_five_column_spread", fmt_f64(max_spread));
    kv("winner_only_H0_H4", format!("{winner_only:?}"));
    for (mask, n) in &tie_sets {
        let members: Vec<String> = (0..5)
            .filter(|k| mask & (1 << k) != 0)
            .map(|k| format!("H{k}"))
            .collect();
        kv(
            &format!("tie_set_{mask:05b}"),
            format!("{} n={n}", members.join("/")),
        );
    }
    for k in 0..5 {
        kv(
            &format!("col_hap{}_stats", FROZEN_IDX[k]),
            format!(
                "rows=123 min={} max={} mean={}",
                fmt_f64(col_min[k]),
                fmt_f64(col_max[k]),
                fmt_f64(col_sum[k] / 123.0)
            ),
        );
    }

    let post: Vec<_> = pipe_cells
        .iter()
        .filter(|c| c.stage == "post_kernel")
        .collect();
    kv("post_kernel_cells", post.len().to_string());
    kv(
        "post_kernel_n_haps",
        post.first()
            .map(|c| c.n_haps.to_string())
            .unwrap_or_default(),
    );
    kv(
        "post_kernel_n_reads",
        post.first()
            .map(|c| c.n_reads.to_string())
            .unwrap_or_default(),
    );
    let frozen_fnv_u: HashSet<u64> = FROZEN_FNV
        .iter()
        .map(|s| u64::from_str_radix(s, 16).unwrap())
        .collect();
    let post_five: Vec<_> = post
        .iter()
        .filter(|c| frozen_fnv_u.contains(&c.hap_fnv))
        .collect();
    kv("post_kernel_five_hap_cells", post_five.len().to_string());

    let qnames: HashSet<(String, u16)> = snap
        .ad_row_qname
        .iter()
        .cloned()
        .zip(snap.ad_row_flags.iter().copied())
        .collect();
    assert_eq!(qnames.len(), 123);
    let post_five_retain: Vec<_> = post_five
        .iter()
        .filter(|c| qnames.contains(&(c.qname.clone(), c.flags)))
        .collect();
    kv(
        "post_kernel_five_hap_on_123",
        post_five_retain.len().to_string(),
    );

    let mut n_len = 0usize;
    let mut n_planes_ok = 0usize;
    let geno = &outcome.genotyping_reads;
    let mut seen = HashSet::new();
    for &idx in &snap.ad_row_read_index {
        if !seen.insert(idx) {
            continue;
        }
        let rec = geno.get(idx).expect("geno read");
        n_len += rec.seq_len();
        let (bases, bq, ins, del, gcp) = pairhmm_planes(rec, &cfg);
        assert_eq!(bases.len(), bq.len());
        assert_eq!(bq.len(), ins.len());
        assert_eq!(ins.len(), del.len());
        assert_eq!(del.len(), gcp.len());
        if gcp.iter().all(|&g| g == GATK_PARITY_DEFAULT_GCP) {
            n_planes_ok += 1;
        }
        let q = String::from_utf8_lossy(rec.qname()).into_owned();
        if q == REP_QNAME && rec.flags() == REP_FLAGS {
            kv(
                "rep_pairhmm_read_input",
                format!(
                    "len={} mapq={} bq_fnv={} ins_fnv={} del_fnv={} gcp_all_10={} has_BI={} has_BD={}",
                    bases.len(),
                    rec.mapq(),
                    fnv1a64_hex(&bq),
                    fnv1a64_hex(&ins),
                    fnv1a64_hex(&del),
                    gcp.iter().all(|&g| g == 10),
                    bam_indel_phred(rec, b"BI").is_some(),
                    bam_indel_phred(rec, b"BD").is_some(),
                ),
            );
            kv(
                "rep_five_haps_share_read_plane",
                "YES: one BQ/IQ/DQ/GCP plane per read, five haplotype base arrays",
            );
        }
    }
    kv("n_123_quality_planes", n_planes_ok.to_string());
    kv(
        "mean_pairhmm_read_len",
        format!("{:.3}", n_len as f64 / 123.0),
    );
    kv(
        "five_column_shared_inputs",
        "identical read order, sequences, BQ, IQ, DQ, PCR, GCP, mode, padding; haplotype bases differ",
    );

    let stacked = KNOWN_ALIGNED_RESIDUAL * 123.0;
    let rust_pl_cont = -10.0 * RUST_HOM_ALT_GL;
    let min_gl_to_cross = (rust_pl_cont - 3517.5) / 10.0;
    kv("rust_hom_alt_gl", fmt_f64(RUST_HOM_ALT_GL));
    kv("rust_continuous_pl", format!("{rust_pl_cont:.17}"));
    kv(
        "min_gl_displacement_to_cross_3517.5",
        format!("{min_gl_to_cross:.17}"),
    );
    kv(
        "known_backend_residual_per_cell",
        format!("{KNOWN_ALIGNED_RESIDUAL:.3e} (6R.82 aligned max |delta|)"),
    );
    kv("stacked_residual_123", format!("{stacked:.17}"));
    kv(
        "known_backend_residual_sufficient",
        if stacked < min_gl_to_cross {
            "NO: 7.41e-6 × 123 same-sign stack is smaller than the 3517.5 crossing"
        } else {
            "YES"
        },
    );
    assert!(
        stacked < min_gl_to_cross,
        "stacked GKL residual cannot explain the PL rounding crossing"
    );

    kv(
        "haplotype_inputs_java_equivalent",
        "UNVERIFIED vs executable; source path is assembly k-best bases after Java-equivalent trim; FNV identity preserved into PairHMM",
    );
    kv(
        "read_inputs_java_equivalent",
        "UNVERIFIED vs executable sequences; membership 123 = retainEvidence (6R.248); clip/hardClipToRegion is Java-equivalent",
    );
    kv(
        "quality_inputs_java_equivalent",
        "YES at source: BQ threshold 18, MAPQ cap, PCR Conservative, GOP Q45, GCP 10 (PairHMMLikelihoodCalculationEngine)",
    );
    kv(
        "transforms",
        "k-best hap → trim_to (Java-eq) → PairHMM hap.bases (Java-eq; CIGAR unused) | reads → finalize+hardClipToRegion (Java-eq) → BQ/PCR (Java-eq) → NEON f64 kernel (not GKL float)",
    );
    kv(
        "first_divergent_operation",
        "not pooling; not allele floor; not GKL residual stack; remaining unverified Java hap sequences / GKL-float cells for these five FNVs",
    );
    kv(
        "classification",
        "HAPLOTYPE_INPUT_DIVERGENCE (PAIRHMM_NUMERICAL_DIVERGENCE insufficient to cross 3517.5; Java hap dump absent)",
    );
}
