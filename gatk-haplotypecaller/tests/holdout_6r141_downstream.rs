//! 6R.141: first haplotype-population boundary after unbounded SeqGraph k-best.
//!
//! Diagnostic-only policy. Does not change the production default.
//! Skipped unless `HOLDOUT_6R141=1`.
//!
//! ```text
//! HOLDOUT_6R141=1 GATK_RS_EXPERIMENTAL_KBEST_POLICY=unbounded_diagnostic \
//!   cargo test -p gatk-haplotypecaller --test holdout_6r141_downstream -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary};
use gatk_haplotypecaller::assembly_region_finalize::{
    create_graph_reference_read, records_to_assembly_reads,
};
use gatk_haplotypecaller::cigar::{Cigar, CigarOperator};
use gatk_haplotypecaller::kbest_haplotype::KBestPath;
use gatk_haplotypecaller::read_threading_assembler::{
    build_threading_graph_for_seq_assembly, extract_haplotypes_from_seq_kbest_paths,
};
use gatk_haplotypecaller::seq_graph::{SeqGraph, SeqGraphCleanupStatus};
use gatk_haplotypecaller::seq_kbest_haplotype::find_best_haplotypes_seq_graph_with_policy;
use gatk_haplotypecaller::seq_kbest_resource_policy::SeqKbestResourcePolicy;
use gatk_haplotypecaller::{
    assemble_reads_with_finalized, call_disposition, flatten_assembly_regions,
    traverse_assembly_region_walker, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    Haplotype, ReadFilterParams, WalkerTraversalConfig,
};
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const JAVA_KBEST_DUMP_REL: &str = "parity/giab/runs/local-pairhmm-diff/6r133_java.txt";
const JAVA_TRIM_DUMP_REL: &str = "parity/giab/runs/local-pairhmm-diff/6r131_java.txt";
const ORACLE_TSV_REL: &str = "docs/parity/6R.141_KBEST_ORACLE.tsv";
const TARGET: u64 = 29_456_196;
const JAVA_ACTIVE_START: u64 = 29_455_995;
const JAVA_ACTIVE_END: u64 = 29_456_293;
const JAVA_PAD_START: u64 = 29_455_895;
const JAVA_PAD_END: u64 = 29_456_393;
const JAVA_TRIM_START: u64 = 29_455_977;
const JAVA_TRIM_END: u64 = 29_456_301;
const EXPECTED_REF_HASH: &str = "659741f99b7f78a5";
const EXPECTED_FINALIZED: usize = 379;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R141\t{key}\t{}", value.as_ref());
}

fn fnv1a64(bases: &[u8]) -> u64 {
    let mut h = 0xcbf29ce484222325u64;
    for &b in bases {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x100000001b3);
    }
    h
}

fn hex(h: u64) -> String {
    format!("{h:016x}")
}

fn hap_hash(h: &Haplotype) -> String {
    hex(fnv1a64(&h.bases))
}

fn path_hash(graph: &SeqGraph, p: &KBestPath) -> String {
    hex(fnv1a64(&graph.path_bases_bytes(p.start, &p.edges)))
}

fn unique_of(haps: &[Haplotype]) -> BTreeSet<String> {
    haps.iter().map(hap_hash).collect()
}

fn unique_paths(graph: &SeqGraph, paths: &[KBestPath]) -> BTreeSet<String> {
    paths.iter().map(|p| path_hash(graph, p)).collect()
}

struct EnvGuard {
    key: &'static str,
    prior: Option<String>,
}

impl EnvGuard {
    fn set(key: &'static str, value: &str) -> Self {
        let prior = std::env::var(key).ok();
        std::env::set_var(key, value);
        Self { key, prior }
    }
}

impl Drop for EnvGuard {
    fn drop(&mut self) {
        match &self.prior {
            Some(v) => std::env::set_var(self.key, v),
            None => std::env::remove_var(self.key),
        }
    }
}

struct NamedRow {
    hash: String,
    score: f64,
    is_ref: bool,
    len: usize,
    cigar: String,
    align: usize,
    loc: String,
}

fn hap_row(h: &Haplotype) -> NamedRow {
    NamedRow {
        hash: hap_hash(h),
        score: h.score,
        is_ref: h.is_reference,
        len: h.bases.len(),
        cigar: h
            .cigar
            .as_ref()
            .map(|c| c.to_gatk_string())
            .unwrap_or_else(|| ".".into()),
        align: h.alignment_start_hap_wrt_ref,
        loc: h
            .genome_loc
            .map(|g| format!("{}-{}", g.start_1based(), g.end_1based()))
            .unwrap_or_else(|| ".".into()),
    }
}

fn dump_haps(stage: &str, haps: &[Haplotype]) {
    let uniq = unique_of(haps);
    let n_ref = haps.iter().filter(|h| h.is_reference).count();
    let mut counts = BTreeMap::<String, usize>::new();
    for h in haps {
        *counts.entry(hap_hash(h)).or_insert(0) += 1;
    }
    let n_dup = counts.values().filter(|n| **n > 1).count();
    kv(
        "set",
        format!(
            "stage={stage}\tn_cols={}\tunique={}\tn_ref={n_ref}\tdup_seqs={n_dup}",
            haps.len(),
            uniq.len()
        ),
    );
    for (i, h) in haps.iter().enumerate() {
        let r = hap_row(h);
        kv(
            "hap",
            format!(
                "stage={stage}\tidx={i}\thash={}\tisRef={}\tlen={}\tscore={}\talign={}\tloc={}\tcigar={}",
                r.hash, r.is_ref, r.len, r.score, r.align, r.loc, r.cigar
            ),
        );
    }
}

fn compare(label: &str, rust: &BTreeSet<String>, java: &BTreeSet<String>) -> (usize, usize, usize) {
    let common = rust.intersection(java).count();
    let java_only = java.difference(rust).count();
    let rust_only = rust.difference(java).count();
    kv(
        "compare",
        format!(
            "label={label}\trust={}\tjava={}\tCOMMON={common}\tJAVA_ONLY={java_only}\tRUST_ONLY={rust_only}",
            rust.len(),
            java.len()
        ),
    );
    for h in java.difference(rust).take(8) {
        kv("java_only", format!("label={label}\thash={h}"));
    }
    for h in rust.difference(java).take(8) {
        kv("rust_only", format!("label={label}\thash={h}"));
    }
    (common, java_only, rust_only)
}

fn order_mismatch(a: &[String], b: &[String]) -> usize {
    let n = a.len().min(b.len());
    let mut m = a.len().abs_diff(b.len());
    for i in 0..n {
        if a[i] != b[i] {
            m += 1;
        }
    }
    m
}

fn parse_kv_line(line: &str) -> BTreeMap<String, String> {
    let mut m = BTreeMap::new();
    for tok in line.split('\t') {
        if let Some((k, v)) = tok.split_once('=') {
            m.insert(k.to_string(), v.to_string());
        }
    }
    m
}

struct JavaKbest {
    ordered: Vec<String>,
    scores: BTreeMap<String, f64>,
    is_ref: BTreeMap<String, bool>,
    set: BTreeSet<String>,
}

fn load_java_kbest(path: &Path) -> JavaKbest {
    let mut by_rank: BTreeMap<usize, (String, f64, bool)> = BTreeMap::new();
    let text = std::fs::read_to_string(path).expect("java kbest dump");
    for line in text.lines() {
        if !line.starts_with("6R133\tkbest\t") {
            continue;
        }
        let m = parse_kv_line(line);
        if m.get("kmer").map(|s| s.as_str()) != Some("25") {
            continue;
        }
        if m.get("K").map(|s| s.as_str()) != Some("128") {
            continue;
        }
        let rank: usize = m.get("rank").and_then(|s| s.parse().ok()).expect("rank");
        let hash = m.get("hash").cloned().expect("hash");
        let score: f64 = m
            .get("score")
            .and_then(|s| s.parse().ok())
            .unwrap_or(f64::NAN);
        let is_ref = m.get("isRef").map(|s| s == "true").unwrap_or(false);
        by_rank.insert(rank, (hash, score, is_ref));
    }
    let mut ordered = Vec::new();
    let mut scores = BTreeMap::new();
    let mut is_ref = BTreeMap::new();
    for (_, (h, sc, ir)) in by_rank {
        scores.insert(h.clone(), sc);
        is_ref.insert(h.clone(), ir);
        ordered.push(h);
    }
    let set = ordered.iter().cloned().collect();
    JavaKbest {
        ordered,
        scores,
        is_ref,
        set,
    }
}

fn load_oracle_tsv(path: &Path) -> Vec<String> {
    let text = std::fs::read_to_string(path).expect("oracle tsv");
    let mut out = Vec::new();
    for line in text.lines() {
        if line.starts_with('#') || line.starts_with("rank\t") || line.is_empty() {
            continue;
        }
        let hash = line.split('\t').nth(1).expect("oracle hash").to_string();
        out.push(hash);
    }
    out
}

struct JavaAssemble {
    ordered: Vec<String>,
    is_ref: BTreeMap<String, bool>,
    scores: BTreeMap<String, String>,
    set: BTreeSet<String>,
}

fn load_java_assemble(path: &Path) -> JavaAssemble {
    let text = std::fs::read_to_string(path).expect("java assemble dump");
    let mut by_idx: BTreeMap<usize, (String, bool, String)> = BTreeMap::new();
    for line in text.lines() {
        if !line.starts_with("6R133\tassemble\t") {
            continue;
        }
        let m = parse_kv_line(line);
        let idx: usize = m.get("idx").and_then(|s| s.parse().ok()).expect("idx");
        let hash = m.get("hash").cloned().expect("hash");
        let is_ref = m.get("isRef").map(|s| s == "true").unwrap_or(false);
        let score = m.get("score").cloned().unwrap_or_else(|| ".".into());
        by_idx.insert(idx, (hash, is_ref, score));
    }
    let mut ordered = Vec::new();
    let mut is_ref = BTreeMap::new();
    let mut scores = BTreeMap::new();
    for (_, (h, ir, sc)) in by_idx {
        is_ref.insert(h.clone(), ir);
        scores.insert(h.clone(), sc);
        ordered.push(h);
    }
    let set = ordered.iter().cloned().collect();
    JavaAssemble {
        ordered,
        is_ref,
        scores,
        set,
    }
}

struct JavaUntrimmed {
    cigars: BTreeMap<String, String>,
    align: BTreeMap<String, usize>,
    is_ref: BTreeMap<String, bool>,
    ordered: Vec<String>,
    set: BTreeSet<String>,
}

fn load_java_untrimmed(path: &Path) -> JavaUntrimmed {
    let text = std::fs::read_to_string(path).expect("java trim dump");
    let mut by_idx: BTreeMap<usize, String> = BTreeMap::new();
    let mut cigars = BTreeMap::new();
    let mut align = BTreeMap::new();
    let mut is_ref = BTreeMap::new();
    for line in text.lines() {
        if !line.starts_with("6R131\thap\t") {
            continue;
        }
        let m = parse_kv_line(line);
        if m.get("stage").map(|s| s.as_str()) != Some("untrimmed") {
            continue;
        }
        let idx: usize = m.get("idx").and_then(|s| s.parse().ok()).unwrap_or(0);
        let hash = m.get("hash").cloned().expect("hash");
        if let Some(c) = m.get("cigar") {
            cigars.insert(hash.clone(), c.clone());
        }
        if let Some(a) = m.get("alignStart").and_then(|s| s.parse().ok()) {
            align.insert(hash.clone(), a);
        }
        is_ref.insert(
            hash.clone(),
            m.get("isRef").map(|s| s == "true").unwrap_or(false),
        );
        by_idx.insert(idx, hash);
    }
    let ordered: Vec<String> = by_idx.into_values().collect();
    let set = ordered.iter().cloned().collect();
    JavaUntrimmed {
        cigars,
        align,
        is_ref,
        ordered,
        set,
    }
}

fn load_java_trimmed(path: &Path) -> (Vec<String>, BTreeSet<String>) {
    let text = std::fs::read_to_string(path).expect("java trim dump");
    let mut by_idx: BTreeMap<usize, String> = BTreeMap::new();
    for line in text.lines() {
        if !line.starts_with("6R131\thap\t") {
            continue;
        }
        let m = parse_kv_line(line);
        if m.get("stage").map(|s| s.as_str()) != Some("trimmed") {
            continue;
        }
        if let (Some(idx), Some(h)) = (m.get("idx").and_then(|s| s.parse().ok()), m.get("hash")) {
            by_idx.insert(idx, h.clone());
        }
    }
    let ordered: Vec<String> = by_idx.into_values().collect();
    let set = ordered.iter().cloned().collect();
    (ordered, set)
}

fn score_close(a: f64, b: f64) -> bool {
    if a.is_nan() && b.is_nan() {
        return true;
    }
    if a.is_nan() || b.is_nan() {
        return false;
    }
    (a - b).abs() <= 1e-8
}

#[test]
fn holdout_6r141_downstream() {
    if std::env::var("HOLDOUT_6R141").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R141=1");
        return;
    }
    let _policy = EnvGuard::set("GATK_RS_EXPERIMENTAL_KBEST_POLICY", "unbounded_diagnostic");
    kv(
        "parity_mode",
        "unbounded_diagnostic\tnote=diagnostic_not_production_default",
    );

    let root = repo_root();
    let java_kbest = load_java_kbest(&root.join(JAVA_KBEST_DUMP_REL));
    let oracle = load_oracle_tsv(&root.join(ORACLE_TSV_REL));
    let java_assemble = load_java_assemble(&root.join(JAVA_KBEST_DUMP_REL));
    let java_untrimmed = load_java_untrimmed(&root.join(JAVA_TRIM_DUMP_REL));
    let (java_trim_order, java_trimmed) = load_java_trimmed(&root.join(JAVA_TRIM_DUMP_REL));
    assert_eq!(java_kbest.set.len(), 128);
    assert_eq!(oracle.len(), 128);
    assert_eq!(java_kbest.ordered, oracle);
    kv("java_kbest_loaded", java_kbest.set.len().to_string());
    kv("oracle_tsv_loaded", oracle.len().to_string());
    kv(
        "java_assemble_vs_kbest_unique",
        format!(
            "assemble={}\tkbest={}\tequal={}",
            java_assemble.set.len(),
            java_kbest.set.len(),
            java_assemble.set == java_kbest.set
        ),
    );
    kv(
        "java_untrimmed_vs_kbest_unique",
        format!(
            "untrimmed={}\tkbest={}\tequal={}",
            java_untrimmed.set.len(),
            java_kbest.set.len(),
            java_untrimmed.set == java_kbest.set
        ),
    );

    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
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
        .expect("ActiveFull")
        .clone();
    let mut java_bounds = covering;
    java_bounds.start = GenomePosition::new_1based(JAVA_ACTIVE_START);
    java_bounds.end = GenomePosition::new_1based(JAVA_ACTIVE_END);
    java_bounds.extended_start = GenomePosition::new_1based(JAVA_PAD_START);
    java_bounds.extended_end = GenomePosition::new_1based(JAVA_PAD_END);

    let args = CallRegionArgs::strict_java();
    let mut ref_cache = ReferenceWindowCache::new(ref_fasta.clone(), 4);
    let mut owned = java_bounds.clone();
    let assembled =
        assemble_reads_with_finalized(&mut owned, &dict, &mut ref_cache, &args.assemble)
            .expect("assemble");
    let graph_ref_hash = hex(fnv1a64(assembled.assembly.reference_bases()));
    assert_eq!(graph_ref_hash, EXPECTED_REF_HASH);
    assert_eq!(assembled.finalized_reads.len(), EXPECTED_FINALIZED);

    let padded_ref = gatk_haplotypecaller::assembly_region_finalize::assembly_reference_read(
        &dict,
        &mut ref_cache,
        &java_bounds,
    )
    .expect("pad ref");
    let graph_ref = create_graph_reference_read(&padded_ref, &java_bounds, &dict);
    let graph_reads = records_to_assembly_reads(&assembled.finalized_reads);
    let mut assembler = args.assemble.assembler.clone();
    assembler.dangling_java_exact = true;
    let rt = build_threading_graph_for_seq_assembly(
        &graph_ref,
        &graph_reads,
        25,
        &assembler,
        false,
        false,
    )
    .expect("rt")
    .expect("k=25 graph");
    let mut seq = SeqGraph::from_assembly_graph(&rt);
    seq.clean_non_ref_paths();
    assert_eq!(
        seq.cleanup_seq_graph(),
        SeqGraphCleanupStatus::AssembledSomeVariation
    );
    kv(
        "seqgraph",
        format!("nodes={}\tedges={}", seq.node_count(), seq.edge_count()),
    );

    let report = find_best_haplotypes_seq_graph_with_policy(
        &seq,
        128,
        SeqKbestResourcePolicy::UnboundedDiagnostic,
    )
    .expect("kbest");
    let rust_kbest_set = unique_paths(&seq, &report.paths);
    let rust_kbest_order: Vec<String> = report.paths.iter().map(|p| path_hash(&seq, p)).collect();
    let (k_common, k_java_only, k_rust_only) =
        compare("kbest_unique", &rust_kbest_set, &java_kbest.set);
    kv(
        "kbest_order",
        format!(
            "rust_n={}\tjava_n={}\tpositional_mismatches={}\tidentical_order={}",
            rust_kbest_order.len(),
            java_kbest.ordered.len(),
            order_mismatch(&rust_kbest_order, &java_kbest.ordered),
            rust_kbest_order == java_kbest.ordered
        ),
    );
    if rust_kbest_order != java_kbest.ordered {
        for i in 0..rust_kbest_order.len().min(java_kbest.ordered.len()) {
            if rust_kbest_order[i] != java_kbest.ordered[i] {
                kv(
                    "kbest_order_first_mismatch",
                    format!(
                        "idx={i}\trust={}\tjava={}",
                        rust_kbest_order[i], java_kbest.ordered[i]
                    ),
                );
                break;
            }
        }
    }
    for (i, p) in report.paths.iter().enumerate() {
        let h = path_hash(&seq, p);
        kv(
            "kbest",
            format!(
                "rank={i}\thash={h}\tscore={:.8}\tisRef={}\tjava_rank={}\tjava_score={}\tjava_isRef={}",
                p.score,
                p.is_reference,
                java_kbest
                    .ordered
                    .iter()
                    .position(|x| x == &h)
                    .map(|r| r.to_string())
                    .unwrap_or_else(|| ".".into()),
                java_kbest
                    .scores
                    .get(&h)
                    .map(|s| format!("{s:.8}"))
                    .unwrap_or_else(|| ".".into()),
                java_kbest
                    .is_ref
                    .get(&h)
                    .map(|b| b.to_string())
                    .unwrap_or_else(|| ".".into())
            ),
        );
    }
    let mut score_mismatch = 0usize;
    let mut isref_mismatch = 0usize;
    for p in &report.paths {
        let h = path_hash(&seq, p);
        if let Some(&js) = java_kbest.scores.get(&h) {
            if !score_close(p.score, js) {
                score_mismatch += 1;
            }
        }
        if let Some(&jr) = java_kbest.is_ref.get(&h) {
            if p.is_reference != jr {
                isref_mismatch += 1;
            }
        }
    }
    kv(
        "kbest_field_mismatch",
        format!("score={score_mismatch}\tisRef={isref_mismatch}"),
    );
    assert_eq!(k_common, 128, "unbounded k-best unique set must match Java");
    assert_eq!(k_java_only, 0);
    assert_eq!(k_rust_only, 0);
    assert!(!report.diagnostics.resource_limit_hit);
    kv(
        "kbest_oracle",
        format!(
            "COMMON={k_common}\tJAVA_ONLY={k_java_only}\tRUST_ONLY={k_rust_only}\tpeak_logical={}\tlimit_hit={}",
            report.diagnostics.peak_logical_frontier, report.diagnostics.resource_limit_hit
        ),
    );

    let mut ref_hap = Haplotype::new(graph_ref.bases.as_slice(), true);
    let mut ref_cigar = Cigar::new();
    ref_cigar.push(ref_hap.bases.len(), CigarOperator::Match);
    ref_hap.cigar = Some(ref_cigar);
    let extracted = extract_haplotypes_from_seq_kbest_paths(
        &report.paths,
        &seq,
        25,
        &ref_hap,
        ref_hap.bases.len(),
        &args.assemble.assembler.haplotype_to_reference_sw,
    )
    .expect("extract");
    dump_haps("extract", &extracted);
    let rust_extract_set = unique_of(&extracted);
    let extract_vs_kbest = compare("extract_vs_kbest", &rust_extract_set, &rust_kbest_set);
    let extract_vs_java = compare("extract_vs_java_kbest", &rust_extract_set, &java_kbest.set);
    kv(
        "extract_drops",
        format!(
            "paths={}\textracted_cols={}\textracted_unique={}\tdropped_cols={}",
            report.paths.len(),
            extracted.len(),
            rust_extract_set.len(),
            report.paths.len().saturating_sub(extracted.len())
        ),
    );

    let mut cigar_mismatch = 0usize;
    let mut cigar_missing = 0usize;
    let mut align_mismatch = 0usize;
    for h in &extracted {
        let hash = hap_hash(h);
        let cigar = h
            .cigar
            .as_ref()
            .map(|c| c.to_gatk_string())
            .unwrap_or_else(|| ".".into());
        match java_untrimmed.cigars.get(&hash) {
            Some(jc) if jc == &cigar => {}
            Some(jc) => {
                cigar_mismatch += 1;
                if cigar_mismatch <= 8 {
                    kv(
                        "cigar_mismatch",
                        format!(
                            "stage=extract\thash={hash}\trust={cigar}\tjava={jc}\talign_rust={}\talign_java={}",
                            h.alignment_start_hap_wrt_ref,
                            java_untrimmed
                                .align
                                .get(&hash)
                                .map(|a| a.to_string())
                                .unwrap_or_else(|| ".".into())
                        ),
                    );
                }
            }
            None => cigar_missing += 1,
        }
        if let Some(&ja) = java_untrimmed.align.get(&hash) {
            if h.alignment_start_hap_wrt_ref != ja {
                align_mismatch += 1;
            }
        }
    }
    kv(
        "extract_cigar_vs_java_untrimmed",
        format!(
            "mismatch={cigar_mismatch}\tmissing_java={cigar_missing}\talign_mismatch={align_mismatch}"
        ),
    );

    dump_haps("assemble", &assembled.assembly.haplotypes);
    let rust_assemble_set = unique_of(&assembled.assembly.haplotypes);
    let assemble_vs_extract = compare("assemble_vs_extract", &rust_assemble_set, &rust_extract_set);
    let assemble_vs_java = compare(
        "assemble_vs_java_assemble",
        &rust_assemble_set,
        &java_assemble.set,
    );
    let rust_assemble_order: Vec<String> =
        assembled.assembly.haplotypes.iter().map(hap_hash).collect();
    kv(
        "assemble_order",
        format!(
            "rust_n={}\tjava_n={}\tpositional_mismatches={}\tidentical_order={}",
            rust_assemble_order.len(),
            java_assemble.ordered.len(),
            order_mismatch(&rust_assemble_order, &java_assemble.ordered),
            rust_assemble_order == java_assemble.ordered
        ),
    );
    if rust_assemble_order != java_assemble.ordered {
        for i in 0..rust_assemble_order.len().min(java_assemble.ordered.len()) {
            if rust_assemble_order[i] != java_assemble.ordered[i] {
                kv(
                    "assemble_order_first_mismatch",
                    format!(
                        "idx={i}\trust={}\tjava={}",
                        rust_assemble_order[i], java_assemble.ordered[i]
                    ),
                );
                break;
            }
        }
    }
    let mut assemble_isref_mismatch = 0usize;
    let mut assemble_cigar_mismatch = 0usize;
    for h in &assembled.assembly.haplotypes {
        let hash = hap_hash(h);
        if let Some(&jr) = java_assemble.is_ref.get(&hash) {
            if h.is_reference != jr {
                assemble_isref_mismatch += 1;
                kv(
                    "assemble_isref_mismatch",
                    format!("hash={hash}\trust={}\tjava={jr}", h.is_reference),
                );
            }
        }
        let cigar = h
            .cigar
            .as_ref()
            .map(|c| c.to_gatk_string())
            .unwrap_or_else(|| ".".into());
        if let Some(jc) = java_untrimmed.cigars.get(&hash) {
            if jc != &cigar {
                assemble_cigar_mismatch += 1;
                if assemble_cigar_mismatch <= 8 {
                    kv(
                        "cigar_mismatch",
                        format!("stage=assemble\thash={hash}\trust={cigar}\tjava={jc}"),
                    );
                }
            }
        }
    }
    kv(
        "assemble_fields",
        format!(
            "isref_mismatch={assemble_isref_mismatch}\tcigar_mismatch={assemble_cigar_mismatch}\tevents={}",
            assembled.assembly.variation_events().len()
        ),
    );
    for (i, ev) in assembled.assembly.variation_events().iter().enumerate() {
        kv(
            "event",
            format!(
                "idx={i}\tstart={}\tend={}\tref={}\talt={}",
                ev.start_1based.get(),
                ev.end_1based.get(),
                ev.ref_allele,
                ev.alt_allele
            ),
        );
    }

    let mut trim_region = java_bounds.clone();
    trim_region.extended_start = GenomePosition::new_1based(JAVA_TRIM_START);
    trim_region.extended_end = GenomePosition::new_1based(JAVA_TRIM_END);
    let trimmed = assembled.assembly.trim_to(&trim_region).expect("trim_to");
    dump_haps("trim", &trimmed.haplotypes);
    let rust_trim_set = unique_of(&trimmed.haplotypes);
    let rust_trim_order: Vec<String> = trimmed.haplotypes.iter().map(hap_hash).collect();
    let trim_vs_java = compare("trim_vs_java_trimmed", &rust_trim_set, &java_trimmed);
    kv(
        "trim_order",
        format!(
            "rust_n={}\tjava_n={}\tpositional_mismatches={}\tidentical_order={}",
            rust_trim_order.len(),
            java_trim_order.len(),
            order_mismatch(&rust_trim_order, &java_trim_order),
            rust_trim_order == java_trim_order
        ),
    );
    let mut java_assemble_from_kbest = java_kbest.ordered[1..].to_vec();
    java_assemble_from_kbest.push(java_kbest.ordered[0].clone());
    kv(
        "assemble_permutation",
        format!(
            "rust_assemble_eq_java_kbest={}\tjava_assemble_eq_kbest_tail_plus_ref={}",
            rust_assemble_order == java_kbest.ordered,
            java_assemble.ordered == java_assemble_from_kbest
        ),
    );

    let extract_pop_changed = extract_vs_kbest.1 != 0 || extract_vs_kbest.2 != 0;
    let assemble_pop_changed = assemble_vs_extract.1 != 0 || assemble_vs_extract.2 != 0;
    let assemble_vs_java_split = assemble_vs_java.1 != 0 || assemble_vs_java.2 != 0;
    let cigar_div = cigar_mismatch > 0 || assemble_cigar_mismatch > 0;
    let trim_split = trim_vs_java.1 != 0 || trim_vs_java.2 != 0;
    let kbest_order_only = rust_kbest_order != java_kbest.ordered && k_common == 128;
    let assemble_order_only =
        rust_assemble_order != java_assemble.ordered && assemble_vs_java.0 == 128;
    let trim_order_match = rust_trim_order == java_trim_order;

    let (classification, first_divergence) = if k_common != 128 {
        (
            "KBEST_UNIQUE_REGRESSION",
            "unbounded k-best unique set vs Java top-128",
        )
    } else if extract_pop_changed {
        (
            "EXTRACT_POPULATION_DIVERGENCE",
            "extract_haplotypes_from_seq_kbest_paths",
        )
    } else if assemble_pop_changed || assemble_vs_java_split {
        (
            "ASSEMBLE_POPULATION_DIVERGENCE",
            "assemble_reads after k-best (merge_rt / normalize / AssemblyResultSet)",
        )
    } else if cigar_div && trim_split {
        (
            "CIGAR_CAUSAL_FOR_TRIM",
            "haplotype CIGAR / SW alignment after extract (causal for trim_to unique set)",
        )
    } else if trim_split {
        ("TRIM_POPULATION_DIVERGENCE", "AssemblyResultSet::trim_to")
    } else {
        (
            "NO_DIVERGENCE",
            "none: unique haplotype population identical through trim_to",
        )
    };
    kv(
        "classification",
        format!(
            "class={classification}\tfirst_divergence={first_divergence}\textract_pop_changed={extract_pop_changed}\tassemble_pop_changed={assemble_pop_changed}\tassemble_vs_java_split={assemble_vs_java_split}\tcigar_div={cigar_div}\ttrim_split={trim_split}\tkbest_order_only={kbest_order_only}\tassemble_order_only={assemble_order_only}\ttrim_order_match={trim_order_match}\talign_mismatch={align_mismatch}\talign_note=java_untrimmed_500_vs_rust_0_not_causal_for_unique_set"
        ),
    );
    kv(
        "summary",
        format!(
            "kbest_java_common={k_common}/128\tkbest_java_only={k_java_only}\tkbest_rust_only={k_rust_only}\textract_unique={}\tassemble_unique={}\ttrim_unique={}\ttrim_COMMON={}\ttrim_JAVA_ONLY={}\ttrim_RUST_ONLY={}",
            rust_extract_set.len(),
            rust_assemble_set.len(),
            rust_trim_set.len(),
            trim_vs_java.0,
            trim_vs_java.1,
            trim_vs_java.2
        ),
    );

    assert_eq!(
        rust_kbest_order, java_kbest.ordered,
        "unbounded k-best completed-path order must match Java rank order"
    );
    assert_eq!(
        rust_assemble_order, java_kbest.ordered,
        "Rust assemble list is k-best order with in-place REF flag"
    );
    assert_eq!(
        java_assemble.ordered, java_assemble_from_kbest,
        "Java assemble list is k-best[1..] + canonical REF last"
    );
    assert_eq!(extract_vs_kbest, (128, 0, 0));
    assert_eq!(assemble_vs_extract, (128, 0, 0));
    assert_eq!(assemble_vs_java, (128, 0, 0));
    assert_eq!(cigar_mismatch, 0);
    assert_eq!(assemble_cigar_mismatch, 0);
    assert_eq!(assemble_isref_mismatch, 0);
    assert_eq!(trim_vs_java, (25, 0, 0));
    assert_eq!(rust_trim_order, java_trim_order);
    let _ = (
        extract_vs_java,
        java_untrimmed.ordered,
        java_assemble.scores,
        java_untrimmed.is_ref,
        trimmed,
        kbest_order_only,
        assemble_order_only,
    );
}
