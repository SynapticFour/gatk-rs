//! 6R.249: proof-only. Trace the TGTTTG allele-row at `20:29455649 T/TGTTTG`.
//!
//! 6R.248 closed the homozygous calculator. The remaining split is already in
//! per-read `L(read|TGTTTG)` after pooling the five-haplotype TGTTTG EventMap
//! mapper pool. This round reconstructs that pool from the live mapper, not
//! from the allele string.
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! PRODUCTION CHANGE: NONE. Do not patch PairHMM, pooling, mapper, floor, or GLs.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r249_tgtttg_allele_row -- --nocapture --test-threads=1
//! HOLDOUT_6R249=1 cargo test -p gatk-haplotypecaller --test holdout_6r249_tgtttg_allele_row -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::{
    build_per_haplotype_variation_events, overlapping_events, remap_alt_onto_longer_ref,
};
use gatk_haplotypecaller::hc_genotyping_engine::take_colocated_merge_numerics;
use gatk_haplotypecaller::read_threading_assembler::DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, region_likelihoods_to_rows,
    traverse_assembly_region_walker, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    Haplotype, HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
};
use std::collections::{HashMap, HashSet};
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
const EMPTY_POOL: f64 = -50.0;
const FLOOR_CAP: f64 = -4.5;
const SAMPLE_ROWS: usize = 8;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R249\t{key}\t{}", value.as_ref());
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

fn src_contains(rel: &str, needle: &str) -> bool {
    fs::read_to_string(repo_root().join(rel))
        .unwrap_or_default()
        .contains(needle)
}

/// GATK 4.4 `createAlleleMapper` EventMap walk (production
/// `java_create_allele_mapper_pools`). Dual alt membership is allowed.
fn java_create_allele_mapper_pools(
    n_haps: usize,
    hap_events: &gatk_haplotypecaller::event_map::PerHaplotypeVariationEvents,
    loc: u64,
    long_ref: &str,
    alts: &[String],
    emit_spanning: bool,
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
        let mut in_ref = false;
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
            } else if emit_spanning {
                if let Some(ai) = alts.iter().position(|a| a == "*") {
                    in_alt[ai] = true;
                }
                break;
            } else {
                in_ref = true;
                break;
            }
        }
        let n_alts_hit = in_alt.iter().filter(|&&b| b).count();
        if !in_ref && n_alts_hit == 0 {
            continue;
        }
        if in_ref {
            pools[0].push(i);
        }
        for (ai, hit) in in_alt.iter().enumerate() {
            if *hit {
                pools[ai + 1].push(i);
            }
        }
    }
    pools
}

fn pool_max_log10(indices: &[usize], lls: &[f64]) -> f64 {
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

#[test]
fn forensic_6r249_source_pool_max_is_java_marginalize_max() {
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");
    assert_eq!(DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH, 128);

    let mod_rs = fs::read_to_string(
        repo_root().join("gatk-haplotypecaller/src/hc_genotyping_engine/mod.rs"),
    )
    .expect("mod.rs");
    let pool_body = fn_body(&mod_rs, "fn pool_max_log10");
    assert!(
        pool_body.contains("fold(f64::NEG_INFINITY, f64::max)"),
        "pool_max_log10 is per-read max over haplotype LLs in the allele pool"
    );
    assert!(
        !pool_body.contains("log10_sum_log10"),
        "pool_max_log10 is not log-sum-exp"
    );
    kv(
        "rust_pool_max_log10",
        "L(read|allele)=max_i L(read|hap_i) over mapper pool; empty → -50",
    );

    let assign = fs::read_to_string(
        repo_root().join("gatk-haplotypecaller/src/hc_genotyping_engine/genotype_assign.rs"),
    )
    .expect("genotype_assign.rs");
    let merge = fn_body(&assign, "fn try_genotype_colocated_snp_indel_merge");
    assert!(merge.contains("pool_max_log10(pool, row)"));
    assert!(merge.contains("apply_java_marginal_normalize_n(&mut marg)"));
    let coloc = fn_body(&assign, "fn colocated_merge_allele_pools");
    assert!(
        coloc.contains("java_create_allele_mapper_pools"),
        "production merge uses the Java EventMap mapper when the cache is present"
    );
    let java_mapper = fn_body(&assign, "fn java_create_allele_mapper_pools");
    assert!(java_mapper.contains("overlapping_events(hap_events.events_for(i), loc)"));
    assert!(java_mapper.contains("remap_alt_onto_longer_ref"));
    kv(
        "rust_mapper",
        "colocated_merge_allele_pools → java_create_allele_mapper_pools (GATK 4.4 createAlleleMapper)",
    );
    kv(
        "java_pool",
        "AssemblyBasedCallerUtils.createAlleleMapper + AlleleLikelihoods.marginalize: max of old-allele LLs (not log-sum-exp); empty → -Inf then empty-pool path",
    );
    assert!(src_contains(
        "gatk-haplotypecaller/src/hc_genotyping_engine/mod.rs",
        "const LOG10_GLOBAL_READ_MISMATCHING_RATE: f64 = -4.5"
    ));
}

#[test]
fn forensic_6r249_live_tgtttg_five_hap_pool() {
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");
    kv("target", "20:29455649 T/TGTTTG");

    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing chr20_tiny BAM/REF");
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
    kv(
        "active_full",
        format!(
            "{}:{}-{}",
            region.contig,
            region.start.get(),
            region.end.get()
        ),
    );
    let outcome = HaplotypeCallerEngine::call_region(
        region,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");
    let snap = take_colocated_merge_numerics()
        .into_iter()
        .find(|s| s.loc == TARGET)
        .expect("merge snapshot");

    assert_eq!(snap.long_ref, TARGET_REF);
    assert_eq!(snap.alts, [UNUSED_ALT, TARGET_ALT]);
    assert_eq!(snap.assigned_gt, vec![0, 2]);
    assert_eq!(snap.remaining_keep_indices, vec![0, 2]);
    kv("merged_allele_0", "T (REF)");
    kv("merged_allele_1", "TTTG (unused ALT)");
    kv("merged_allele_2", "TGTTTG (kept ALT → emitted ALT1)");
    kv(
        "mapping",
        "merged 2/2 = TGTTTG/TGTTTG → unused-ALT subset [0,3,5] → emitted 1/1",
    );
    kv("assigned_gt", format!("{:?}", snap.assigned_gt));
    kv("keep_indices", format!("{:?}", snap.remaining_keep_indices));
    kv("pool_sizes", format!("{:?}", snap.pool_sizes));
    kv(
        "java_style_pool_sizes",
        format!("{:?}", snap.java_style_pool_sizes),
    );
    kv("n_reads", snap.n_reads.to_string());
    kv("n_haps", snap.n_haps.to_string());
    kv(
        "n_allele_floor_clips",
        snap.n_allele_floor_clips.to_string(),
    );
    kv(
        "hap_event_signatures",
        format!("{:?}", snap.hap_event_signatures_at_loc),
    );
    kv(
        "n_haps_unassigned_java",
        snap.n_haps_unassigned_java.to_string(),
    );
    kv(
        "n_haps_in_multiple_java_alts",
        snap.n_haps_in_multiple_java_alts.to_string(),
    );
    kv("subset_pl", format!("{:?}", snap.subset_pl));
    kv(
        "subset_pl_java_style_pools",
        format!("{:?}", snap.subset_pl_java_style_pools),
    );
    kv(
        "subset_pl_no_allele_floor",
        format!("{:?}", snap.subset_pl_no_allele_floor),
    );

    assert_eq!(snap.n_reads, 123);
    assert_eq!(snap.n_haps, outcome.assembly.haplotypes.len());
    assert_eq!(snap.pool_sizes.len(), 3);
    assert_eq!(
        snap.pool_sizes[2], 5,
        "production TGTTTG pool is five haplotypes"
    );
    assert_eq!(
        snap.java_style_pool_sizes, snap.pool_sizes,
        "production pools are the Java EventMap mapper"
    );
    assert_eq!(snap.n_allele_floor_clips, 0);
    assert_eq!(snap.subset_pl_no_allele_floor, snap.subset_pl);

    let haps = &outcome.assembly.haplotypes;
    let pad = outcome.assembly.padded_reference_start_1based();
    let ref_bytes = outcome.assembly.reference_bases();
    let hap_events = build_per_haplotype_variation_events(
        haps,
        ref_bytes,
        pad,
        outcome.assembly.max_mnp_distance(),
        "20",
    );
    let alts = vec![UNUSED_ALT.to_string(), TARGET_ALT.to_string()];
    let pools = java_create_allele_mapper_pools(haps.len(), &hap_events, TARGET, "T", &alts, true);
    kv(
        "reconstructed_pool_sizes",
        format!(
            "REF={} TTTG={} TGTTTG={}",
            pools[0].len(),
            pools[1].len(),
            pools[2].len()
        ),
    );
    assert_eq!(
        pools.iter().map(|p| p.len()).collect::<Vec<_>>(),
        snap.pool_sizes
    );

    let loc_pos = GenomePosition::new_1based(TARGET);
    kv("rust_REF_hap_indices", format!("{:?}", pools[0]));
    kv("rust_TTTG_hap_indices", format!("{:?}", pools[1]));
    kv("rust_TGTTTG_hap_indices", format!("{:?}", pools[2]));

    let tg = &pools[2];
    assert_eq!(tg.len(), 5);
    let mut dual_tttg = 0usize;
    let mut dual_ref = 0usize;
    for (rank, &hi) in tg.iter().enumerate() {
        let h = &haps[hi];
        let spanning = overlapping_events(hap_events.events_for(hi), TARGET);
        let at_loc: Vec<String> = spanning
            .iter()
            .filter(|e| e.start_1based == loc_pos)
            .map(|e| format!("{}/{}", e.ref_allele, e.alt_allele))
            .collect();
        let remapped: Vec<String> = spanning
            .iter()
            .filter(|e| e.start_1based == loc_pos)
            .filter_map(|e| {
                if e.ref_allele.len() == 1 {
                    Some(e.alt_allele.clone())
                } else {
                    remap_alt_onto_longer_ref(&e.ref_allele, &e.alt_allele, "T")
                }
            })
            .collect();
        let role = if h.is_reference {
            "reference haplotype"
        } else {
            "assembly k-best haplotype"
        };
        kv(
            &format!("tgtttg_H{rank}"),
            format!(
                "idx={hi} fnv={hash} is_ref={is_ref} len={len} kmer={kmer} score={score:.6} cigar={cigar} align_start={align} loc={loc} role={role} events_at_loc={events} remapped={remap} spanning_n={span_n} ll_col={hi}",
                hash = fnv1a64_hex(&h.bases),
                is_ref = h.is_reference,
                len = h.bases.len(),
                kmer = h.kmer_size,
                score = h.score,
                cigar = cigar_str(h),
                align = h.alignment_start_hap_wrt_ref,
                loc = h.genome_loc.map(|g| format!("{}-{}", g.start_1based(), g.end_1based())).unwrap_or_else(|| ".".to_string()),
                events = at_loc.join("+"),
                remap = remapped.join("+"),
                span_n = spanning.len(),
            ),
        );
        kv(
            &format!("tgtttg_H{rank}_bases"),
            String::from_utf8_lossy(&h.bases).into_owned(),
        );
        if pools[1].contains(&hi) {
            dual_tttg += 1;
        }
        if pools[0].contains(&hi) {
            dual_ref += 1;
        }
        assert!(
            remapped.iter().any(|a| a == TARGET_ALT),
            "TGTTTG pool member must carry a remapped/native EventMap alt TGTTTG at loc"
        );
        assert!(
            !h.is_reference,
            "TGTTTG pool must not contain the reference haplotype"
        );
    }
    kv("tgtttg_also_in_TTTG", dual_tttg.to_string());
    kv("tgtttg_also_in_REF", dual_ref.to_string());
    kv(
        "subsetting_before_pool",
        "none: hap_rows keep all region hap columns; mapper selects five indices; pool_max_log10 is max over those columns",
    );

    let mut rust_ref_events = Vec::new();
    for &hi in &pools[0] {
        let at_loc: Vec<String> = overlapping_events(hap_events.events_for(hi), TARGET)
            .iter()
            .filter(|e| e.start_1based == loc_pos)
            .map(|e| format!("{}/{}", e.ref_allele, e.alt_allele))
            .collect();
        rust_ref_events.push(format!(
            "{hi}:{}",
            if at_loc.is_empty() {
                "none".to_string()
            } else {
                at_loc.join("+")
            }
        ));
    }
    kv("rust_REF_membership", rust_ref_events.join("; "));
    let mut rust_tttg_events = Vec::new();
    for &hi in &pools[1] {
        let at_loc: Vec<String> = overlapping_events(hap_events.events_for(hi), TARGET)
            .iter()
            .filter(|e| e.start_1based == loc_pos)
            .map(|e| format!("{}/{}", e.ref_allele, e.alt_allele))
            .collect();
        rust_tttg_events.push(format!(
            "{hi}:{} fnv={}",
            at_loc.join("+"),
            fnv1a64_hex(&haps[hi].bases)
        ));
    }
    kv("rust_TTTG_membership", rust_tttg_events.join("; "));
    kv(
        "java_membership_dump",
        "gatk_image_absent: no executable Java haplotype list for ActiveFull 20:29455560-29455744",
    );

    let retain: HashSet<usize> = snap.ad_row_read_index.iter().copied().collect();
    assert_eq!(retain.len(), 123);
    let subset: Vec<_> = outcome
        .read_likelihoods
        .iter()
        .filter(|c| retain.contains(&c.read_index.get()))
        .cloned()
        .collect();
    let hap_rows = region_likelihoods_to_rows(&subset, haps.len());
    let mut by_read: HashMap<usize, &gatk_haplotypecaller::genotyping::ReadLikelihoodRow> =
        HashMap::new();
    for row in &hap_rows {
        by_read.entry(row.read_index).or_insert(row);
    }
    assert_eq!(hap_rows.len(), 123);

    let mut n_unique_max = 0usize;
    let mut n_ties = 0usize;
    let mut winner_counts = [0usize; 5];
    let mut n_empty = 0usize;
    let mut spreads: Vec<f64> = Vec::new();
    let mut n_raw_ne_floor = 0usize;
    let mut n_floor_clip_allele2 = 0usize;
    let mut first_floor_row: Option<usize> = None;
    let mut first_tie_row: Option<usize> = None;
    let mut first_unique_row: Option<usize> = None;
    let mut n_reconstruct_mismatch = 0usize;
    let mut first_mismatch: Option<usize> = None;
    let mut max_abs_floor_delta = 0.0f64;
    let mut sum_abs_floor_delta = 0.0f64;

    for (ri, &read_idx) in snap.ad_row_read_index.iter().enumerate() {
        let row = by_read.get(&read_idx).expect("hap row");
        let five: Vec<f64> = tg
            .iter()
            .map(|&hi| {
                row.haplotype_log10_likelihoods
                    .get(hi)
                    .copied()
                    .unwrap_or(f64::NEG_INFINITY)
            })
            .collect();
        let raw = pool_max_log10(tg, &row.haplotype_log10_likelihoods);
        if !raw.is_finite() || raw == EMPTY_POOL {
            n_empty += 1;
        }
        let n_at_max = five.iter().filter(|&&v| v == raw).count();
        if n_at_max == 1 {
            n_unique_max += 1;
            let mut argmax = 0usize;
            let mut best5 = f64::NEG_INFINITY;
            for (k, &v) in five.iter().enumerate() {
                if v > best5 {
                    best5 = v;
                    argmax = k;
                }
            }
            winner_counts[argmax] += 1;
            if first_unique_row.is_none() {
                first_unique_row = Some(ri);
                let qname = snap.ad_row_qname.get(ri).cloned().unwrap_or_default();
                let flags = snap.ad_row_flags.get(ri).copied().unwrap_or(0);
                kv(
                    "first_unique_max_row",
                    format!(
                        "row={ri} read_index={read_idx} qname={qname} flags={flags} five=[{}] max={} winner=idx{} H{argmax} ad_row2={}",
                        five.iter().map(|v| fmt_f64(*v)).collect::<Vec<_>>().join(", "),
                        fmt_f64(raw),
                        tg[argmax],
                        fmt_f64(snap.ad_row_lls[ri][2]),
                    ),
                );
            }
        } else if n_at_max > 1 {
            n_ties += 1;
            if first_tie_row.is_none() {
                first_tie_row = Some(ri);
            }
        }
        let finite: Vec<f64> = five.iter().copied().filter(|v| v.is_finite()).collect();
        if finite.len() >= 2 {
            let mn = finite.iter().copied().fold(f64::INFINITY, f64::min);
            let mx = finite.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            spreads.push(mx - mn);
        }
        let pooled3: Vec<f64> = pools
            .iter()
            .map(|p| pool_max_log10(p, &row.haplotype_log10_likelihoods))
            .collect();
        let best = pooled3.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let floor = best + FLOOR_CAP;
        let after = if pooled3[2].is_finite() && pooled3[2] < floor {
            n_floor_clip_allele2 += 1;
            if first_floor_row.is_none() {
                first_floor_row = Some(ri);
            }
            floor
        } else {
            pooled3[2]
        };
        let ad = snap.ad_row_lls[ri][2];
        let d = (after - ad).abs();
        if d > 1e-15 {
            n_reconstruct_mismatch += 1;
            if first_mismatch.is_none() {
                first_mismatch = Some(ri);
            }
        }
        let fd = (raw - ad).abs();
        max_abs_floor_delta = max_abs_floor_delta.max(fd);
        sum_abs_floor_delta += fd;
        if (raw - ad).abs() > 1e-12 {
            n_raw_ne_floor += 1;
        }
        if ri < SAMPLE_ROWS {
            let qname = snap.ad_row_qname.get(ri).cloned().unwrap_or_default();
            let flags = snap.ad_row_flags.get(ri).copied().unwrap_or(0);
            let mut argmax = 0usize;
            let mut best5 = f64::NEG_INFINITY;
            for (k, &v) in five.iter().enumerate() {
                if v > best5 {
                    best5 = v;
                    argmax = k;
                }
            }
            kv(
                &format!("row{ri}_pool"),
                format!(
                    "read_index={read_idx} qname={qname} flags={flags} five=[{}] max={} argmax=H{argmax} n_at_max={n_at_max} ad_row2={} raw_vs_ad={fd:.3e}",
                    five.iter().map(|v| fmt_f64(*v)).collect::<Vec<_>>().join(", "),
                    fmt_f64(raw),
                    fmt_f64(ad),
                ),
            );
        }
    }

    let spread_min = spreads.iter().copied().fold(f64::INFINITY, f64::min);
    let spread_max = spreads.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let spread_mean = if spreads.is_empty() {
        0.0
    } else {
        spreads.iter().sum::<f64>() / spreads.len() as f64
    };
    kv("pool_rows", "123");
    kv("n_unique_max", n_unique_max.to_string());
    kv("n_ties", n_ties.to_string());
    kv(
        "unique_max_winner_counts_H0_H4",
        format!("{winner_counts:?}"),
    );
    kv("n_empty_pool", n_empty.to_string());
    kv("first_unique_row", format!("{first_unique_row:?}"));
    kv("first_tie_row", format!("{first_tie_row:?}"));
    kv(
        "five_hap_spread_min_max_mean",
        format!("{spread_min:.17} / {spread_max:.17} / {spread_mean:.17}"),
    );
    kv("n_raw_ne_ad_row", n_raw_ne_floor.to_string());
    kv(
        "n_allele2_floor_clips_reconstructed",
        n_floor_clip_allele2.to_string(),
    );
    kv("first_floor_row", format!("{first_floor_row:?}"));
    kv(
        "n_reconstruct_mismatch_vs_ad_row",
        n_reconstruct_mismatch.to_string(),
    );
    kv("first_reconstruct_mismatch", format!("{first_mismatch:?}"));
    kv(
        "raw_vs_ad_row_max_mean_abs",
        format!(
            "{:.3e} / {:.3e}",
            max_abs_floor_delta,
            sum_abs_floor_delta / 123.0
        ),
    );
    kv(
        "allele_floor_control",
        "n_allele_floor_clips=0; subset_pl_no_allele_floor == subset_pl; reconstructed allele2 clips=0; floor is a no-op on TGTTTG",
    );
    kv(
        "calculator_residual",
        "4.55e-13 non-causal (6R.248); not re-opened",
    );

    assert_eq!(n_empty, 0);
    assert_eq!(
        n_reconstruct_mismatch, 0,
        "pool_max_log10 + floor reconstructs ad_row allele 2"
    );
    assert_eq!(n_floor_clip_allele2, 0);
    assert_eq!(n_unique_max + n_ties, 123);

    kv(
        "pooling_operation_identical_to_java_source",
        "YES: both are per-read max over mapper-pool haplotype LLs; not log-sum-exp",
    );
    kv(
        "haplotype_membership_algorithm_identical",
        "YES on this 24-hap Rust EventMap: production pools == java_create_allele_mapper_pools == snap.java_style_pool_sizes",
    );
    kv(
        "haplotype_population_vs_java_executable",
        "UNTESTED: GATK 4.4 image absent; no Java hap list / hap-LL dump for this ActiveFull",
    );
    kv(
        "pooled_column_vs_java",
        "UNTESTED: no Java TGTTTG column dump; Rust column equals reconstructed max of the five Rust hap LLs",
    );
    kv(
        "first_unequal_operation",
        "per-read L(read|TGTTTG)=max of five EventMap-mapped hap PairHMM LLs; Java cell comparison absent",
    );
    kv("case_c_pooling_divergence", "NO");
    kv(
        "classification",
        "ALLELE_LIKELIHOOD_INPUT_DIVERGENCE (Case A haplotype-population identity vs Java executable untested; Case C pooling operator closed as identical max)",
    );
}
