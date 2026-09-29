//! 6R.251: proof-only. Exact five TGTTTG haplotype *base sequences* at
//! `20:29455649 T/TGTTTG`.
//!
//! 6R.250 closed the shared read/quality plane and showed the known
//! GKL-float residual cannot cross 3517.5. This round asks whether the
//! five Rust `hap.bases` arrays are Java-equivalent sequences.
//! PRODUCTION CHANGE: NONE. Direct Java executable: NO.
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r251_haplotype_sequence_identity -- --nocapture --test-threads=1
//! HOLDOUT_6R251=1 cargo test -p gatk-haplotypecaller --test holdout_6r251_haplotype_sequence_identity -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary};
use gatk_haplotypecaller::event_map::{
    build_per_haplotype_variation_events, overlapping_events, VariationEvent,
};
use gatk_haplotypecaller::hc_genotyping_engine::take_colocated_merge_numerics;
use gatk_haplotypecaller::read_threading_assembler::DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH;
use gatk_haplotypecaller::{
    assemble_reads_with_finalized, begin_hap_list_observe, begin_trim_io_observe, call_disposition,
    flatten_assembly_regions, traverse_assembly_region_walker, AssembleReadsArgs,
    AssemblyRegionCallDisposition, CallRegionArgs, CigarOperator, GenomeLoc, Haplotype,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
};
use std::collections::BTreeMap;
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

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R251\t{key}\t{}", value.as_ref());
}

fn fnv1a64_hex(bases: &[u8]) -> String {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in bases {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    format!("{h:016x}")
}

fn hex_bytes(bases: &[u8]) -> String {
    bases.iter().map(|b| format!("{b:02x}")).collect()
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

fn apply_events_to_window(
    window: &[u8],
    window_start_1based: u64,
    events: &[VariationEvent],
) -> Result<Vec<u8>, String> {
    let mut sorted: Vec<&VariationEvent> = events.iter().collect();
    sorted.sort_by_key(|e| e.start_1based.get());
    let mut out = window.to_vec();
    for ev in sorted.into_iter().rev() {
        let start = ev.start_1based.get();
        if start < window_start_1based {
            return Err(format!(
                "event start {start} precedes window {window_start_1based}"
            ));
        }
        let off = (start - window_start_1based) as usize;
        let rb = ev.ref_allele.as_bytes();
        let ab = ev.alt_allele.as_bytes();
        if off + rb.len() > out.len() {
            return Err(format!(
                "event {}/{} at {start} overruns window len={}",
                ev.ref_allele,
                ev.alt_allele,
                out.len()
            ));
        }
        if &out[off..off + rb.len()] != rb {
            return Err(format!(
                "event {}/{} at {start} ref mismatch at off={off}",
                ev.ref_allele, ev.alt_allele
            ));
        }
        out.splice(off..off + rb.len(), ab.iter().copied());
    }
    Ok(out)
}

fn first_diff(a: &[u8], b: &[u8]) -> Option<usize> {
    let n = a.len().min(b.len());
    for i in 0..n {
        if a[i] != b[i] {
            return Some(i);
        }
    }
    if a.len() != b.len() {
        Some(n)
    } else {
        None
    }
}

fn common_prefix(a: &[u8], b: &[u8]) -> usize {
    a.iter().zip(b.iter()).take_while(|(x, y)| x == y).count()
}

fn common_suffix(a: &[u8], b: &[u8]) -> usize {
    a.iter()
        .rev()
        .zip(b.iter().rev())
        .take_while(|(x, y)| x == y)
        .count()
}

fn hamming(a: &[u8], b: &[u8]) -> usize {
    let n = a.len().min(b.len());
    let mut d = a.len().abs_diff(b.len());
    for i in 0..n {
        if a[i] != b[i] {
            d += 1;
        }
    }
    d
}

fn mismatch_runs(a: &[u8], b: &[u8]) -> Vec<(usize, usize)> {
    let n = a.len().min(b.len());
    let mut runs = Vec::new();
    let mut i = 0usize;
    while i < n {
        if a[i] == b[i] {
            i += 1;
            continue;
        }
        let start = i;
        while i < n && a[i] != b[i] {
            i += 1;
        }
        runs.push((start, i));
    }
    if a.len() != b.len() {
        let lo = n;
        let hi = a.len().max(b.len());
        runs.push((lo, hi));
    }
    runs
}

fn cigar_insertions(h: &Haplotype) -> Vec<(usize, String)> {
    let mut out = Vec::new();
    let Some(cigar) = h.cigar.as_ref() else {
        return out;
    };
    let mut hap_pos = 0usize;
    for el in &cigar.elements {
        match el.operator {
            CigarOperator::Insertion | CigarOperator::SoftClip => {
                if el.operator == CigarOperator::Insertion {
                    let end = (hap_pos + el.length).min(h.bases.len());
                    out.push((
                        hap_pos,
                        String::from_utf8_lossy(&h.bases[hap_pos..end]).into_owned(),
                    ));
                }
                hap_pos += el.length;
            }
            CigarOperator::Match => {
                hap_pos += el.length;
            }
            CigarOperator::Deletion => {}
            _ => {
                if el.operator.consumes_read_bases() {
                    hap_pos += el.length;
                }
            }
        }
    }
    out
}

fn cigar_match_mismatches(h: &Haplotype, window_ref: &[u8]) -> Vec<(usize, u8, u8)> {
    let mut out = Vec::new();
    let Some(cigar) = h.cigar.as_ref() else {
        return out;
    };
    let mut hap_pos = 0usize;
    let mut ref_pos = 0usize;
    for el in &cigar.elements {
        match el.operator {
            CigarOperator::Match => {
                for _ in 0..el.length {
                    if hap_pos < h.bases.len() && ref_pos < window_ref.len() {
                        let hb = h.bases[hap_pos];
                        let rb = window_ref[ref_pos];
                        if hb != rb {
                            out.push((hap_pos, rb, hb));
                        }
                    }
                    hap_pos += 1;
                    ref_pos += 1;
                }
            }
            CigarOperator::Insertion | CigarOperator::SoftClip => {
                hap_pos += el.length;
            }
            CigarOperator::Deletion => {
                ref_pos += el.length;
            }
            _ => {
                if el.operator.consumes_read_bases() {
                    hap_pos += el.length;
                }
                if el.operator.consumes_reference_bases() {
                    ref_pos += el.length;
                }
            }
        }
    }
    out
}

#[test]
fn forensic_6r251_source_haplotype_lifecycle() {
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");
    kv(
        "direct_java_executable",
        "NO: gatk_image_absent; not pulled",
    );
    assert_eq!(DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH, 128);

    let rta = fs::read_to_string(
        repo_root().join("gatk-haplotypecaller/src/read_threading_assembler.rs"),
    )
    .expect("read_threading_assembler.rs");
    let extract = fn_body(&rta, "pub fn extract_haplotypes_from_seq_kbest_paths");
    assert!(extract.contains("graph.path_bases_bytes(path.start, &path.edges)"));
    assert!(extract.contains("Haplotype::new(bases, path.is_reference)"));
    assert!(extract.contains("calculate_haplotype_cigar_for_assembly_with_offset"));
    kv(
        "java_kbest_haplotype",
        "GATK 4.4 KBestHaplotype.haplotype(): new Haplotype(getBases(), isReference()); Path.getBases() concatenates vertex sequences; SW CIGAR is annotation, not a base rewrite",
    );
    kv(
        "rust_kbest_materialize",
        "extract_haplotypes_from_seq_kbest_paths: path_bases_bytes → Haplotype::new(bases) → SW CIGAR/offset; bases are the graph walk",
    );

    let seqg = fs::read_to_string(repo_root().join("gatk-haplotypecaller/src/seq_graph.rs"))
        .expect("seq_graph.rs");
    let path_bases = fn_body(&seqg, "pub fn path_bases_bytes");
    assert!(path_bases.contains("self.vertices[first].sequence.to_vec()"));
    assert!(path_bases.contains("bases.extend_from_slice(&self.vertices[to].sequence)"));
    kv(
        "path_bases",
        "SeqGraph.path_bases_bytes concatenates vertex payloads; no inserted-allele rewrite",
    );

    let hap_src = fs::read_to_string(repo_root().join("gatk-haplotypecaller/src/haplotype.rs"))
        .expect("haplotype.rs");
    let trim = fn_body(&hap_src, "pub fn trim(");
    assert!(trim.contains("get_bases_covering_ref_interval"));
    assert!(trim.contains("trim_cigar_by_reference"));
    kv(
        "java_haplotype_trim",
        "GATK 4.4 Haplotype.trim(loc, ignoreRefState): AlignmentUtils.getBasesCoveringRefInterval then drop leading/trailing insertion CIGAR ops; subsequence of getBases(), not an allele rewrite",
    );
    kv(
        "rust_haplotype_trim",
        "Haplotype::trim: get_bases_covering_ref_interval slice + CIGAR trim; can drop flanking insertion *bases* that sit outside the ref interval; cannot invent new interior bases",
    );

    let ars =
        fs::read_to_string(repo_root().join("gatk-haplotypecaller/src/assembly_result_set.rs"))
            .expect("assembly_result_set.rs");
    let trim_to = fn_body(&ars, "pub fn trim_to(");
    assert!(trim_to.contains("h.trim(&span, true)"));
    assert!(trim_to.contains("index_by_bases"));
    kv(
        "java_trim_to",
        "GATK 4.4 AssemblyResultSet.trimDownHaplotypes: h.trim(span, true) then HashMap keyed by Haplotype.equals (bases+ref flag; ignoreRefState forces non-ref so uniqueness is sequence-only)",
    );
    kv(
        "hardclip_target",
        "hardClipToRegion / clip_finalized_reads_in_place mutates READS, not haplotype bases",
    );

    let ll = fs::read_to_string(repo_root().join("gatk-haplotypecaller/src/likelihood_engine.rs"))
        .expect("likelihood_engine.rs");
    let score = fn_body(&ll, "pub fn score_read_against_haplotypes");
    assert!(score.contains("haplotype_bases: &[&[u8]]"));
    kv(
        "java_pairhmm_hap_input",
        "GATK 4.4 PairHMMLikelihoodCalculationEngine.computeLog10Likelihoods uses Haplotype.getBases() after trimTo; EventMap is not the PairHMM input",
    );
}

#[test]
fn forensic_6r251_haplotype_sequence_identity() {
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
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing chr20_tiny BAM/REF");
        return;
    }
    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, INTERVAL).expect("interval");
    begin_hap_list_observe();
    begin_trim_io_observe();
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
    let hap_snaps = gatk_haplotypecaller::take_hap_list_snaps();
    let trim_span = gatk_haplotypecaller::take_hap_list_trim_span();
    let trim_io = gatk_haplotypecaller::take_trim_io();
    let snap = take_colocated_merge_numerics()
        .into_iter()
        .find(|s| s.loc == TARGET)
        .expect("merge snapshot");
    assert_eq!(snap.n_reads, 123);
    assert_eq!(snap.long_ref, TARGET_REF);
    assert_eq!(snap.alts, [UNUSED_ALT, TARGET_ALT]);
    assert_eq!(snap.pool_sizes[2], 5);

    let mut region_for_assemble = region.clone();
    let mut ref_cache = ReferenceWindowCache::new(ref_fasta.clone(), 4);
    let mut assemble_args = AssembleReadsArgs::default();
    assemble_args.strict_java_assembly = true;
    assemble_args.assembler.dangling_java_exact = true;
    let assembled = assemble_reads_with_finalized(
        &mut region_for_assemble,
        &dict,
        &mut ref_cache,
        &assemble_args,
    )
    .expect("assemble");
    let untrimmed = &assembled.assembly.haplotypes;

    let haps = &outcome.assembly.haplotypes;
    let pad = outcome.assembly.padded_reference_start_1based();
    let full_ref = outcome.assembly.reference_bases();
    let trim = trim_span.expect("trim span");
    kv(
        "trim_span",
        format!(
            "active={}-{} trim={}-{} pad_start={}",
            trim.active_start, trim.active_end, trim.trim_start, trim.trim_end, pad
        ),
    );
    let ref_hap = haps
        .iter()
        .find(|h| h.is_reference)
        .expect("trimmed reference haplotype");
    let trim_loc = ref_hap.genome_loc.expect("ref genome_loc");
    assert_eq!(trim_loc.start_1based(), trim.trim_start);
    assert_eq!(trim_loc.end_1based(), trim.trim_end);
    kv(
        "trimmed_ref",
        format!(
            "len={} fnv={} cigar={} align_start={}",
            ref_hap.bases.len(),
            fnv1a64_hex(&ref_hap.bases),
            cigar_str(ref_hap),
            ref_hap.alignment_start_hap_wrt_ref
        ),
    );

    let hap_events = build_per_haplotype_variation_events(
        haps,
        full_ref,
        pad,
        outcome.assembly.max_mnp_distance(),
        "20",
    );

    let stages: Vec<String> = hap_snaps
        .iter()
        .map(|s| format!("{}:n={}", s.stage, s.n))
        .collect();
    kv("hap_list_stages", stages.join(" | "));

    let mut all_eventmap_exact = true;
    let mut all_insertion_tgtttg = true;
    let mut all_trim_subsequence = true;
    let mut inserted: Vec<String> = Vec::new();
    let mut seqs: Vec<Vec<u8>> = Vec::new();

    for (i, &idx) in FROZEN_IDX.iter().enumerate() {
        let h = &haps[idx];
        let hash = fnv1a64_hex(&h.bases);
        assert_eq!(hash, FROZEN_FNV[i], "frozen FNV H{i}");
        assert_eq!(h.bases.len(), 161);
        assert_eq!(cigar_str(h), "81M5I75M");
        assert!(!h.is_reference);
        seqs.push(h.bases.clone());

        let bases_s = String::from_utf8_lossy(&h.bases).into_owned();
        kv(
            &format!("H{i}_identity"),
            format!(
                "idx={idx} fnv={hash} len=161 cigar=81M5I75M is_ref=false kmer={} score={:.6} align_start={} loc={}-{}",
                h.kmer_size,
                h.score,
                h.alignment_start_hap_wrt_ref,
                h.genome_loc.map(|g| g.start_1based()).unwrap_or(0),
                h.genome_loc.map(|g| g.end_1based()).unwrap_or(0),
            ),
        );
        kv(&format!("H{i}_bases"), &bases_s);
        kv(&format!("H{i}_bases_hex"), hex_bytes(&h.bases));

        let events = hap_events.events_for(idx);
        let ev_s: Vec<String> = events
            .iter()
            .map(|e| {
                format!(
                    "{}:{} {}/{}",
                    e.start_1based.get(),
                    e.end_1based.get(),
                    e.ref_allele,
                    e.alt_allele
                )
            })
            .collect();
        kv(&format!("H{i}_eventmap"), ev_s.join(" | "));
        let spanning = overlapping_events(events, TARGET);
        let has_tgtttg = spanning
            .iter()
            .any(|e| e.ref_allele == TARGET_REF && e.alt_allele == TARGET_ALT);
        assert!(has_tgtttg, "H{i} must carry EventMap T/TGTTTG");

        let ins = cigar_insertions(h);
        kv(
            &format!("H{i}_cigar_insertions"),
            ins.iter()
                .map(|(p, s)| format!("hap_pos={p} bases={s}"))
                .collect::<Vec<_>>()
                .join(" ; "),
        );
        let ins5 = ins
            .iter()
            .find(|(p, s)| *p == 81 && s.len() == 5)
            .map(|(_, s)| s.as_str())
            .unwrap_or("");
        inserted.push(ins5.to_string());
        kv(
            &format!("H{i}_inserted_bases"),
            format!(
                "ref_allele=T alt_allele=TGTTTG cigar_inserted={ins5} hap_pos=81 genome={}",
                TARGET
            ),
        );
        if ins5 != "GTTTG" {
            all_insertion_tgtttg = false;
        }

        let mm = cigar_match_mismatches(h, &ref_hap.bases);
        kv(
            &format!("H{i}_match_mismatches"),
            if mm.is_empty() {
                "none".to_string()
            } else {
                mm.iter()
                    .map(|(p, r, a)| format!("{p}:{}->{}", *r as char, *a as char))
                    .collect::<Vec<_>>()
                    .join(",")
            },
        );

        let window_events: Vec<VariationEvent> = events
            .iter()
            .filter(|e| {
                e.start_1based.get() >= trim.trim_start && e.start_1based.get() <= trim.trim_end
            })
            .cloned()
            .collect();
        match apply_events_to_window(&ref_hap.bases, trim.trim_start, &window_events) {
            Ok(recon) => {
                let diff = first_diff(&recon, &h.bases);
                kv(
                    &format!("H{i}_eventmap_recon"),
                    format!(
                        "recon_len={} hap_len=161 first_diff={} exact={}",
                        recon.len(),
                        diff.map(|d| d.to_string()).unwrap_or_else(|| "none".into()),
                        diff.is_none() && recon.len() == h.bases.len()
                    ),
                );
                if diff.is_some() || recon.len() != h.bases.len() {
                    all_eventmap_exact = false;
                } else {
                    assert_eq!(recon, h.bases, "H{i} EventMap reconstruction");
                }
            }
            Err(e) => {
                all_eventmap_exact = false;
                kv(&format!("H{i}_eventmap_recon"), format!("FAILED: {e}"));
            }
        }

        let only_tgtttg = vec![VariationEvent {
            contig: "20".into(),
            start_1based: gatk_haplotypecaller::GenomePosition::new_1based(TARGET),
            end_1based: gatk_haplotypecaller::GenomePosition::new_1based(TARGET),
            ref_allele: TARGET_REF.into(),
            alt_allele: TARGET_ALT.into(),
        }];
        let only = apply_events_to_window(&ref_hap.bases, trim.trim_start, &only_tgtttg)
            .expect("T/TGTTTG apply");
        kv(
            &format!("H{i}_only_tgtttg_recon_exact"),
            (only == h.bases).to_string(),
        );

        let out_fnv = u64::from_str_radix(FROZEN_FNV[i], 16).unwrap();
        let io_hits: Vec<_> = trim_io.iter().filter(|r| r.output_fnv == out_fnv).collect();
        kv(
            &format!("H{i}_trim_io"),
            io_hits
                .iter()
                .map(|r| {
                    format!(
                        "input_idx={} input_fnv={:016x} input_len={} outcome={} output_len={}",
                        r.input_idx, r.input_fnv, r.input_len, r.outcome, r.output_len
                    )
                })
                .collect::<Vec<_>>()
                .join(" | "),
        );
        assert!(
            io_hits.iter().any(|r| r.outcome == "kept"),
            "H{i} must be a kept trim_to output"
        );

        let parent = io_hits
            .iter()
            .find(|r| r.outcome == "kept")
            .and_then(|r| untrimmed.get(r.input_idx));
        if let Some(u) = parent {
            kv(
                &format!("H{i}_untrimmed"),
                format!(
                    "len={} fnv={} cigar={} kmer={} score={:.6} align_start={}",
                    u.bases.len(),
                    fnv1a64_hex(&u.bases),
                    cigar_str(u),
                    u.kmer_size,
                    u.score,
                    u.alignment_start_hap_wrt_ref
                ),
            );
            let span = GenomeLoc::new(trim.trim_start, trim.trim_end);
            let retrim = u.trim(&span, true).expect("retrim");
            kv(
                &format!("H{i}_trim_arrow"),
                format!(
                    "untrimmed_len={} trimmed_len=161 can_alter_interior=false retrim_fnv={} matches_pairhmm={}",
                    u.bases.len(),
                    fnv1a64_hex(&retrim.bases),
                    retrim.bases == h.bases
                ),
            );
            if retrim.bases != h.bases {
                all_trim_subsequence = false;
            } else {
                assert_eq!(retrim.bases, h.bases);
            }
        } else {
            all_trim_subsequence = false;
            kv(&format!("H{i}_untrimmed"), "NOT_FOUND in assemble pass");
        }

        let assemble_rank = hap_snaps
            .iter()
            .find(|s| s.stage == "after_assemble")
            .and_then(|s| {
                io_hits.iter().find(|r| r.outcome == "kept").and_then(|r| {
                    s.columns.get(r.input_idx).map(|c| {
                        format!(
                            "assemble_idx={} fnv={:016x} len={}",
                            c.index, c.fnv1a, c.len
                        )
                    })
                })
            })
            .unwrap_or_else(|| "absent".into());
        kv(&format!("H{i}_kbest_rank"), assemble_rank);

        for stage in [
            "after_trim_to",
            "after_preserve_untrimmed",
            "after_prune_spillover",
            "pairhmm_input",
        ] {
            let col = hap_snaps.iter().find(|s| s.stage == stage).and_then(|s| {
                s.columns
                    .iter()
                    .find(|c| format!("{:016x}", c.fnv1a) == hash)
            });
            kv(
                &format!("H{i}_{stage}"),
                col.map(|c| {
                    format!(
                        "idx={} len={} cigar={} loc={}-{}",
                        c.index, c.len, c.cigar, c.loc_start, c.loc_end
                    )
                })
                .unwrap_or_else(|| "absent".into()),
            );
        }
    }

    kv("all_five_inserted_bases", inserted.join(","));
    kv(
        "identical_tgtttg_insertion",
        (inserted.iter().all(|s| s == "GTTTG")).to_string(),
    );
    assert!(
        all_insertion_tgtttg,
        "all five CIGAR insertions must be GTTTG (T→TGTTTG)"
    );
    assert!(all_eventmap_exact, "EventMap reconstruction must be exact");
    assert!(
        all_trim_subsequence,
        "re-trim of untrimmed parent must equal PairHMM hap.bases"
    );

    kv("pairwise_header", "i j hamming prefix suffix runs");
    for a in 0..5 {
        for b in (a + 1)..5 {
            let ha = &seqs[a];
            let hb = &seqs[b];
            let runs = mismatch_runs(ha, hb);
            let run_s: Vec<String> = runs.iter().map(|(lo, hi)| format!("{lo}..{hi}")).collect();
            kv(
                &format!("pair_H{a}_H{b}"),
                format!(
                    "hamming={} prefix={} suffix={} runs={}",
                    hamming(ha, hb),
                    common_prefix(ha, hb),
                    common_suffix(ha, hb),
                    if run_s.is_empty() {
                        "none".into()
                    } else {
                        run_s.join(",")
                    }
                ),
            );
        }
    }

    let unique_seqs: std::collections::BTreeSet<_> = seqs.iter().cloned().collect();
    kv("n_unique_sequences", unique_seqs.len().to_string());
    assert_eq!(unique_seqs.len(), 5, "five distinct 161-mers");

    let mut vs_ref: BTreeMap<usize, usize> = BTreeMap::new();
    for (i, s) in seqs.iter().enumerate() {
        vs_ref.insert(
            i,
            cigar_match_mismatches(&haps[FROZEN_IDX[i]], &ref_hap.bases).len(),
        );
        kv(
            &format!("H{i}_vs_trimmed_ref_hamming_in_M"),
            vs_ref[&i].to_string(),
        );
    }

    kv(
        "hardclip_on_haps",
        "NO: hardClipToRegion is a read transform; haplotype bases after trim_to are PairHMM input",
    );
    kv(
        "proven_from_java_source",
        "KBestHaplotype.haplotype copies Path.getBases(); Haplotype.trim is getBasesCoveringRefInterval subsequence; trimDownHaplotypes does not rewrite interior bases; PairHMM uses getBases() after trim",
    );
    kv(
        "proven_in_rust",
        "five 161-mers dumped; EventMap invert equals hap.bases; CIGAR insertion is GTTTG on all five; trim re-application equals PairHMM bases; stages keep FNV",
    );
    kv(
        "not_executable_verified",
        "Java did not materialize these five concrete byte strings in this environment; GATK image absent",
    );
    kv(
        "first_divergent_operation",
        "none inside Rust construction/trim: sequences are EventMap-consistent subsequence of k-best path bases; Java 161-mers unverified",
    );
    kv(
        "classification",
        "HAPLOTYPE_SEQUENCE_INTERNALLY_CONSISTENT (Java byte sequences UNVERIFIED; parent remaining split still HAPLOTYPE_INPUT_DIVERGENCE)",
    );
    kv(
        "next_smallest_arrow",
        "Java 4.4 SeqGraph Path.getBases() / assembleReads haplotype list for ActiveFull 20:29455560-29455744 — whether Java emits these five 161-mers (requires Java hap dump; executable unavailable)",
    );
}
