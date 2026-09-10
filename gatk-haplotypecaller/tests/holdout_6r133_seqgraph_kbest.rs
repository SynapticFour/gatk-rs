//! 6R.133 forensic: SeqGraph path existence vs K=128 selection at `20:29456196`.
//! Skipped unless `HOLDOUT_6R133=1`. Production change: NONE.
//!
//! ```text
//! HOLDOUT_6R133=1 cargo test -p gatk-haplotypecaller --test holdout_6r133_seqgraph_kbest -- --nocapture
//! ```

use gatk_core::reference::{parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary};
use gatk_haplotypecaller::assembly_region_finalize::{
    create_graph_reference_read, records_to_assembly_reads,
};
use gatk_haplotypecaller::read_threading_assembler::build_threading_graph_for_seq_assembly;
use gatk_haplotypecaller::seq_graph::{SeqGraph, SeqGraphCleanupStatus};
use gatk_haplotypecaller::seq_kbest_haplotype::{
    find_best_haplotypes_seq_graph, find_best_haplotypes_seq_graph_forensic, SeqKbestCapPolicy,
};
use gatk_haplotypecaller::{
    assemble_reads_with_finalized, call_disposition, flatten_assembly_regions,
    traverse_assembly_region_walker, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    Haplotype, ReadFilterParams, WalkerTraversalConfig,
};
use std::collections::{BTreeMap, BTreeSet, HashSet, VecDeque};
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const JAVA_ASSEMBLE_DUMP_REL: &str = "parity/giab/runs/local-pairhmm-diff/6r131_java.txt";
const JAVA_KBEST_DUMP_REL: &str = "parity/giab/runs/local-pairhmm-diff/6r133_java.txt";
const TARGET: u64 = 29_456_196;
const JAVA_ACTIVE_START: u64 = 29_455_995;
const JAVA_ACTIVE_END: u64 = 29_456_293;
const JAVA_PAD_START: u64 = 29_455_895;
const JAVA_PAD_END: u64 = 29_456_393;
const EXPECTED_REF_HASH: &str = "659741f99b7f78a5";
const EXPECTED_FINALIZED: usize = 379;
const WORST_PARENT_SLICES: &[&str] = &[
    "249e5a1d0f4250be",
    "d6c499e8eb382a52",
    "c0ea90412035f7c1",
    "534008c250ea7844",
];
const KS: &[usize] = &[128, 256, 512, 1024];

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R133\t{key}\t{}", value.as_ref());
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

fn unique_of(haps: &[Haplotype]) -> BTreeSet<String> {
    haps.iter().map(|h| hex(fnv1a64(&h.bases))).collect()
}

fn load_java_untrimmed_hashes(path: &Path) -> BTreeSet<String> {
    let mut out = BTreeSet::new();
    let Ok(text) = std::fs::read_to_string(path) else {
        return out;
    };
    for line in text.lines() {
        if !line.starts_with("6R131\thap\t") {
            continue;
        }
        let mut stage = None;
        let mut hash = None;
        for tok in line.split('\t').skip(2) {
            if let Some(v) = tok.strip_prefix("stage=") {
                stage = Some(v);
            }
            if let Some(v) = tok.strip_prefix("hash=") {
                hash = Some(v.to_string());
            }
        }
        if stage == Some("untrimmed") {
            if let Some(h) = hash {
                out.insert(h);
            }
        }
    }
    out
}

struct JavaDump {
    assemble: BTreeMap<String, Vec<u8>>,
    kbest: BTreeMap<(usize, usize), Vec<(String, f64, usize)>>, // (kmer, K) -> (hash, score, rank)
    graphs: BTreeMap<usize, DumpGraph>,
}

struct DumpGraph {
    nodes: usize,
    edges: usize,
    verts: Vec<Vec<u8>>,
    outgoing: Vec<Vec<usize>>,
    source: Option<usize>,
    sink: Option<usize>,
}

fn load_java_kbest_dump(path: &Path) -> Option<JavaDump> {
    let text = std::fs::read_to_string(path).ok()?;
    if !text.contains("6R133\t") {
        return None;
    }
    let mut assemble = BTreeMap::new();
    let mut kbest: BTreeMap<(usize, usize), Vec<(String, f64, usize)>> = BTreeMap::new();
    let mut graphs: BTreeMap<usize, DumpGraph> = BTreeMap::new();
    for line in text.lines() {
        let line = if let Some(i) = line.find("6R133\t") {
            &line[i..]
        } else {
            continue;
        };
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 2 {
            continue;
        }
        match parts[1] {
            "assemble" => {
                let mut hash = None;
                let mut seq = None;
                for tok in parts.iter().skip(2) {
                    if let Some(v) = tok.strip_prefix("hash=") {
                        hash = Some(v.to_string());
                    }
                    if let Some(v) = tok.strip_prefix("seq=") {
                        seq = Some(v.as_bytes().to_vec());
                    }
                }
                if let (Some(h), Some(s)) = (hash, seq) {
                    assemble.insert(h, s);
                }
            }
            "seq" => {
                let mut kmer = None;
                let mut nodes = None;
                let mut edges = None;
                for tok in parts.iter().skip(2) {
                    if let Some(v) = tok.strip_prefix("kmer=") {
                        kmer = v.parse().ok();
                    }
                    if let Some(v) = tok.strip_prefix("nodes=") {
                        nodes = v.parse().ok();
                    }
                    if let Some(v) = tok.strip_prefix("edges=") {
                        edges = v.parse().ok();
                    }
                }
                if let (Some(k), Some(n), Some(e)) = (kmer, nodes, edges) {
                    graphs.entry(k).or_insert_with(|| DumpGraph {
                        nodes: n,
                        edges: e,
                        verts: Vec::new(),
                        outgoing: Vec::new(),
                        source: None,
                        sink: None,
                    });
                    if let Some(g) = graphs.get_mut(&k) {
                        g.nodes = n;
                        g.edges = e;
                    }
                }
            }
            "vtx" => {
                let mut kmer = None;
                let mut id = None;
                let mut seq = None;
                let mut is_source = false;
                let mut is_sink = false;
                for tok in parts.iter().skip(2) {
                    if let Some(v) = tok.strip_prefix("kmer=") {
                        kmer = v.parse().ok();
                    }
                    if let Some(v) = tok.strip_prefix("id=") {
                        id = v.parse().ok();
                    }
                    if let Some(v) = tok.strip_prefix("seq=") {
                        seq = Some(v.as_bytes().to_vec());
                    }
                    if *tok == "is_source=true" {
                        is_source = true;
                    }
                    if *tok == "is_sink=true" {
                        is_sink = true;
                    }
                }
                if let (Some(k), Some(i), Some(s)) = (kmer, id, seq) {
                    let g = graphs.entry(k).or_insert_with(|| DumpGraph {
                        nodes: 0,
                        edges: 0,
                        verts: Vec::new(),
                        outgoing: Vec::new(),
                        source: None,
                        sink: None,
                    });
                    if g.verts.len() <= i {
                        g.verts.resize(i + 1, Vec::new());
                        g.outgoing.resize(i + 1, Vec::new());
                    }
                    g.verts[i] = s;
                    if is_source {
                        g.source = Some(i);
                    }
                    if is_sink {
                        g.sink = Some(i);
                    }
                }
            }
            "edge" => {
                let mut kmer: Option<usize> = None;
                let mut from: Option<usize> = None;
                let mut to: Option<usize> = None;
                for tok in parts.iter().skip(2) {
                    if let Some(v) = tok.strip_prefix("kmer=") {
                        kmer = v.parse().ok();
                    }
                    if let Some(v) = tok.strip_prefix("from=") {
                        from = v.parse().ok();
                    }
                    if let Some(v) = tok.strip_prefix("to=") {
                        to = v.parse().ok();
                    }
                }
                if let (Some(k), Some(f), Some(t)) = (kmer, from, to) {
                    let g = graphs.entry(k).or_insert_with(|| DumpGraph {
                        nodes: 0,
                        edges: 0,
                        verts: Vec::new(),
                        outgoing: Vec::new(),
                        source: None,
                        sink: None,
                    });
                    let n = g.verts.len().max(f + 1).max(t + 1);
                    g.verts.resize(n, Vec::new());
                    g.outgoing.resize(n, Vec::<usize>::new());
                    g.outgoing[f].push(t);
                }
            }
            "kbest" => {
                let mut kmer = None;
                let mut kcap = None;
                let mut rank = None;
                let mut hash = None;
                let mut score = None;
                for tok in parts.iter().skip(2) {
                    if let Some(v) = tok.strip_prefix("kmer=") {
                        kmer = v.parse().ok();
                    }
                    if let Some(v) = tok.strip_prefix("K=") {
                        kcap = v.parse().ok();
                    }
                    if let Some(v) = tok.strip_prefix("rank=") {
                        rank = v.parse().ok();
                    }
                    if let Some(v) = tok.strip_prefix("hash=") {
                        hash = Some(v.to_string());
                    }
                    if let Some(v) = tok.strip_prefix("score=") {
                        score = v.parse().ok();
                    }
                }
                if let (Some(kmer), Some(kcap), Some(rank), Some(hash), Some(score)) =
                    (kmer, kcap, rank, hash, score)
                {
                    kbest
                        .entry((kmer, kcap))
                        .or_default()
                        .push((hash, score, rank));
                }
            }
            _ => {}
        }
    }
    Some(JavaDump {
        assemble,
        kbest,
        graphs,
    })
}

#[derive(Debug)]
struct WalkEvidence {
    present: bool,
    first_missing: String,
}

fn exact_walk_generic(
    n: usize,
    verts: &[Vec<u8>],
    outgoing: &[Vec<usize>],
    source: Option<usize>,
    sink: Option<usize>,
    target: &[u8],
) -> WalkEvidence {
    let Some(source) = source else {
        return WalkEvidence {
            present: false,
            first_missing: "no_source".into(),
        };
    };
    let Some(sink) = sink else {
        return WalkEvidence {
            present: false,
            first_missing: "no_sink".into(),
        };
    };
    if source >= n || sink >= n {
        return WalkEvidence {
            present: false,
            first_missing: "endpoint_oob".into(),
        };
    }
    let src = &verts[source];
    if !target.starts_with(src) {
        return WalkEvidence {
            present: false,
            first_missing: format!(
                "source_prefix_mismatch src_len={} tgt_len={}",
                src.len(),
                target.len()
            ),
        };
    }
    let mut q = VecDeque::new();
    let mut seen: HashSet<(usize, usize)> = HashSet::new();
    q.push_back((source, src.len()));
    seen.insert((source, src.len()));
    let mut stuck = String::from("exhausted");
    while let Some((v, off)) = q.pop_front() {
        if off == target.len() {
            if v == sink {
                return WalkEvidence {
                    present: true,
                    first_missing: String::new(),
                };
            }
            for &to in outgoing.get(v).map(|s| s.as_slice()).unwrap_or(&[]) {
                if verts.get(to).is_some_and(|s| s.is_empty()) && seen.insert((to, off)) {
                    q.push_back((to, off));
                }
            }
            continue;
        }
        let outs = outgoing.get(v).map(|s| s.as_slice()).unwrap_or(&[]);
        let mut any = false;
        for &to in outs {
            let seq = verts.get(to).map(|s| s.as_slice()).unwrap_or(&[]);
            if off + seq.len() <= target.len() && &target[off..off + seq.len()] == seq {
                any = true;
                if seen.insert((to, off + seq.len())) {
                    q.push_back((to, off + seq.len()));
                }
            }
        }
        if !any && stuck == "exhausted" {
            stuck = format!("no_outgoing_match v={v} off={off} n_out={}", outs.len());
        }
    }
    WalkEvidence {
        present: false,
        first_missing: stuck,
    }
}

fn exact_walk_seqgraph(graph: &SeqGraph, target: &[u8]) -> WalkEvidence {
    let verts: Vec<Vec<u8>> = graph
        .vertices()
        .iter()
        .map(|v| v.sequence.clone())
        .collect();
    let mut outgoing = vec![Vec::new(); verts.len()];
    for e in graph.edges() {
        if e.from < outgoing.len() {
            outgoing[e.from].push(e.to);
        }
    }
    exact_walk_generic(
        verts.len(),
        &verts,
        &outgoing,
        graph.reference_source_vertex(),
        graph.reference_sink_vertex(),
        target,
    )
}

fn compare(label: &str, rust: &BTreeSet<String>, java: &BTreeSet<String>) {
    kv(
        "compare",
        format!(
            "label={label}\trust={}\tjava={}\tCOMMON={}\tJAVA_ONLY={}\tRUST_ONLY={}",
            rust.len(),
            java.len(),
            rust.intersection(java).count(),
            java.difference(rust).count(),
            rust.difference(java).count()
        ),
    );
}

fn pick_spread<'a, T>(items: &'a [T], n: usize) -> Vec<&'a T> {
    if items.is_empty() || n == 0 {
        return Vec::new();
    }
    if items.len() <= n {
        return items.iter().collect();
    }
    let mut out = Vec::new();
    for i in 0..n {
        let idx = i * (items.len() - 1) / (n - 1);
        out.push(&items[idx]);
    }
    out
}

fn build_cleaned_seqgraph(
    graph_ref: &gatk_haplotypecaller::assembly::AssemblyRead,
    graph_reads: &[gatk_haplotypecaller::assembly::AssemblyRead],
    kmer: usize,
    assembler: &gatk_haplotypecaller::ReadThreadingAssemblerArgs,
) -> Option<(SeqGraph, SeqGraphCleanupStatus)> {
    let rt = build_threading_graph_for_seq_assembly(
        graph_ref,
        graph_reads,
        kmer,
        assembler,
        false,
        false,
    )
    .ok()??;
    let mut seq = SeqGraph::from_assembly_graph(&rt);
    seq.clean_non_ref_paths();
    let status = seq.cleanup_seq_graph();
    Some((seq, status))
}

fn first_k_for(hash: &str, by_k: &BTreeMap<usize, BTreeMap<String, (usize, f64)>>) -> String {
    for &k in KS {
        if let Some(m) = by_k.get(&k) {
            if let Some((rank, score)) = m.get(hash) {
                return format!("K={k}\trank={rank}\tscore={score:.9}");
            }
        }
    }
    "never".into()
}

#[test]
fn holdout_6r133_seqgraph_kbest() {
    if std::env::var("HOLDOUT_6R133").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R133=1");
        return;
    }
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    let java_hashes = load_java_untrimmed_hashes(&root.join(JAVA_ASSEMBLE_DUMP_REL));
    let java_dump = load_java_kbest_dump(&root.join(JAVA_KBEST_DUMP_REL));
    kv("java_untrimmed_loaded", java_hashes.len().to_string());
    kv(
        "java_kbest_dump",
        if java_dump.is_some() {
            "present"
        } else {
            "absent"
        },
    );

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
    kv(
        "matched_input",
        format!(
            "graph_ref={graph_ref_hash}\tlen={}\tfinalized={}",
            assembled.assembly.reference_bases().len(),
            assembled.finalized_reads.len()
        ),
    );
    assert_eq!(graph_ref_hash, EXPECTED_REF_HASH);
    assert_eq!(assembled.finalized_reads.len(), EXPECTED_FINALIZED);

    let rust_assemble = unique_of(&assembled.assembly.haplotypes);
    compare("assemble_vs_java_untrimmed", &rust_assemble, &java_hashes);
    let rust_only: Vec<_> = rust_assemble.difference(&java_hashes).cloned().collect();
    let java_only: Vec<_> = java_hashes.difference(&rust_assemble).cloned().collect();
    kv(
        "private",
        format!(
            "JAVA_ONLY={}\tRUST_ONLY={}",
            java_only.len(),
            rust_only.len()
        ),
    );

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

    let mut rust_graphs: BTreeMap<usize, SeqGraph> = BTreeMap::new();
    for &kmer in &[10usize, 25] {
        match build_cleaned_seqgraph(&graph_ref, &graph_reads, kmer, &assembler) {
            None => kv("seqgraph", format!("kmer={kmer}\tstatus=no_graph")),
            Some((seq, status)) => {
                kv(
                    "seqgraph",
                    format!(
                        "kmer={kmer}\tstatus={status:?}\tnodes={}\tedges={}\tsrc={:?}\tsink={:?}",
                        seq.node_count(),
                        seq.edge_count(),
                        seq.reference_source_vertex(),
                        seq.reference_sink_vertex()
                    ),
                );
                if status == SeqGraphCleanupStatus::AssembledSomeVariation {
                    rust_graphs.insert(kmer, seq);
                }
            }
        }
    }

    let mut rust_kbest: BTreeMap<usize, BTreeMap<usize, BTreeMap<String, (usize, f64)>>> =
        BTreeMap::new();
    let mut rust_k128_paths: BTreeMap<usize, Vec<(String, f64, usize, Vec<u8>)>> = BTreeMap::new();
    for (&kmer, seq) in &rust_graphs {
        for &k in KS {
            let report = find_best_haplotypes_seq_graph_forensic(
                seq,
                k,
                k,
                SeqKbestCapPolicy::unbounded(),
                &[],
            )
            .unwrap_or_else(|_| panic!("kbest kmer={kmer} K={k}"));
            kv(
                "kbest_run",
                format!(
                    "kmer={kmer}\tK={k}\tn={}\texpansions={}\tmax_heap={}\tvisit_refused={}\tskip_heap_pop={}\tskip_heap_exp={}\tskip_exp_cap={}",
                    report.paths.len(),
                    report.expansions,
                    report.max_heap,
                    report.vertex_visit_refused,
                    report.skip_heap_full_at_pop,
                    report.skip_heap_full_at_expand,
                    report.skip_expansion_cap
                ),
            );
            if k == 128 && !report.paths.is_empty() {
                let cut = &report.paths[report.paths.len().min(128).saturating_sub(1)];
                kv(
                    "k128_cutoff",
                    format!(
                        "kmer={kmer}\tscore={:.9}\tn_edges={}",
                        cut.score,
                        cut.edges.len()
                    ),
                );
            }
            let mut map = BTreeMap::new();
            let mut listed = Vec::new();
            for (rank, p) in report.paths.iter().enumerate() {
                let bases = seq.path_bases_bytes(p.start, &p.edges);
                let h = hex(fnv1a64(&bases));
                map.entry(h.clone()).or_insert((rank, p.score));
                if k == 128 {
                    listed.push((h, p.score, rank, bases));
                }
            }
            if k == 128 {
                rust_k128_paths.insert(kmer, listed);
            }
            rust_kbest.entry(kmer).or_default().insert(k, map);
        }
        let prod = find_best_haplotypes_seq_graph_forensic(
            seq,
            128,
            128,
            SeqKbestCapPolicy::production(),
            &[],
        )
        .expect("prod k128");
        kv(
            "kbest_prod",
            format!(
                "kmer={kmer}\tn={}\texpansions={}\tmax_heap={}\tskip_heap_pop={}\tskip_heap_exp={}\tskip_exp_cap={}",
                prod.paths.len(),
                prod.expansions,
                prod.max_heap,
                prod.skip_heap_full_at_pop,
                prod.skip_heap_full_at_expand,
                prod.skip_expansion_cap
            ),
        );
    }

    let mut rust_union_k128 = BTreeSet::new();
    for listed in rust_k128_paths.values() {
        for (h, _, _, _) in listed {
            rust_union_k128.insert(h.clone());
        }
    }
    compare("rust_k128_union_vs_java", &rust_union_k128, &java_hashes);
    compare(
        "rust_k128_union_vs_assemble",
        &rust_union_k128,
        &rust_assemble,
    );

    if let Some(seq) = rust_graphs.get(&25) {
        let prod_paths =
            find_best_haplotypes_seq_graph(seq, 128).expect("production find_best_haplotypes");
        let prod_set: BTreeSet<String> = prod_paths
            .iter()
            .map(|p| hex(fnv1a64(&seq.path_bases_bytes(p.start, &p.edges))))
            .collect();
        compare("prod_k128_vs_java", &prod_set, &java_hashes);
        compare("prod_k128_vs_unbounded_k128", &prod_set, &rust_union_k128);
        compare("prod_k128_vs_assemble", &prod_set, &rust_assemble);

        let mut ref_hap = Haplotype::new(graph_ref.bases.clone(), true);
        let mut ref_cigar = gatk_haplotypecaller::Cigar::new();
        ref_cigar.push(
            ref_hap.bases.len(),
            gatk_haplotypecaller::CigarOperator::Match,
        );
        ref_hap.cigar = Some(ref_cigar);
        let ref_cigar_len = ref_hap.cigar.as_ref().unwrap().reference_length();
        let extracted = gatk_haplotypecaller::read_threading_assembler::extract_haplotypes_from_seq_kbest_paths(
            &prod_paths,
            seq,
            25,
            &ref_hap,
            ref_cigar_len,
            &assembler.haplotype_to_reference_sw,
        )
        .expect("extract");
        let ext_set = unique_of(&extracted);
        kv(
            "extract",
            format!(
                "n_paths={}\tn_extracted={}\tunique={}",
                prod_paths.len(),
                extracted.len(),
                ext_set.len()
            ),
        );
        compare("extract_prod_vs_java", &ext_set, &java_hashes);
        compare("extract_prod_vs_assemble", &ext_set, &rust_assemble);

        let seq_only =
            gatk_haplotypecaller::assemble_from_ref_and_reads(&graph_ref, &graph_reads, &assembler)
                .expect("assemble_from_ref_and_reads");
        let seq_only_set = unique_of(&seq_only.haplotypes);
        kv(
            "assemble_from_ref",
            format!(
                "n={}\tunique={}\tstatus={:?}",
                seq_only.haplotypes.len(),
                seq_only_set.len(),
                seq_only.status
            ),
        );
        compare("assemble_from_ref_vs_java", &seq_only_set, &java_hashes);
        compare(
            "assemble_from_ref_vs_assemble_reads",
            &seq_only_set,
            &rust_assemble,
        );
        compare("assemble_from_ref_vs_extract", &seq_only_set, &ext_set);
    }

    for &kmer in rust_graphs.keys() {
        if let Some(m) = rust_kbest.get(&kmer).and_then(|m| m.get(&128)) {
            let set: BTreeSet<_> = m.keys().cloned().collect();
            compare(&format!("rust_k{kmer}_K128_vs_java"), &set, &java_hashes);
        }
    }

    let java_seqs: BTreeMap<String, Vec<u8>> = java_dump
        .as_ref()
        .map(|d| d.assemble.clone())
        .unwrap_or_default();
    let rust_seq_by_hash: BTreeMap<String, Vec<u8>> = assembled
        .assembly
        .haplotypes
        .iter()
        .map(|h| (hex(fnv1a64(&h.bases)), h.bases.clone()))
        .collect();

    let mut rust_only_scored: Vec<(String, f64, usize, usize)> = Vec::new();
    for (&kmer, listed) in &rust_k128_paths {
        for (h, score, rank, _) in listed {
            if rust_only.iter().any(|x| x == h) {
                rust_only_scored.push((h.clone(), *score, *rank, kmer));
            }
        }
    }
    rust_only_scored.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
    let rust_reps = pick_spread(&rust_only_scored, 5);
    kv("rust_only_scored_n", rust_only_scored.len().to_string());

    let mut needles: Vec<(String, String, Vec<u8>)> = Vec::new();
    for r in rust_reps {
        if let Some(seq) = rust_seq_by_hash.get(&r.0) {
            needles.push((format!("RUST_ONLY:{}", r.0), r.0.clone(), seq.clone()));
        }
    }
    for h in WORST_PARENT_SLICES {
        if let Some(seq) = rust_seq_by_hash.get(*h) {
            needles.push((format!("PARENT_SLICE:{h}"), (*h).to_string(), seq.clone()));
        } else {
            kv("parent_slice_missing", *h);
        }
    }

    let mut java_only_scored: Vec<(String, f64)> = Vec::new();
    if let Some(dump) = java_dump.as_ref() {
        for h in &java_only {
            let mut best: Option<f64> = None;
            for rows in dump.kbest.values() {
                for (hh, score, _) in rows {
                    if hh == h {
                        best = Some(best.map(|b: f64| b.max(*score)).unwrap_or(*score));
                    }
                }
            }
            java_only_scored.push((h.clone(), best.unwrap_or(f64::NEG_INFINITY)));
        }
        java_only_scored.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        for r in pick_spread(&java_only_scored, 5) {
            if let Some(seq) = java_seqs.get(&r.0) {
                needles.push((format!("JAVA_ONLY:{}", r.0), r.0.clone(), seq.clone()));
            }
        }
        for (h, seq) in &java_seqs {
            if java_only.iter().any(|x| x == h)
                && !needles.iter().any(|n| n.1 == *h)
                && needles
                    .iter()
                    .filter(|n| n.0.starts_with("JAVA_ONLY:"))
                    .count()
                    < 5
            {
                needles.push((format!("JAVA_ONLY:{h}"), h.clone(), seq.clone()));
            }
        }
    }

    kv("needles_n", needles.len().to_string());
    for (label, hash, seq) in &needles {
        kv(
            "needle",
            format!("label={label}\thash={hash}\tlen={}", seq.len()),
        );
        for (&kmer, seqg) in &rust_graphs {
            let walk = exact_walk_seqgraph(seqg, seq);
            let first = rust_kbest
                .get(&kmer)
                .map(|m| first_k_for(hash, m))
                .unwrap_or_else(|| "no_graph".into());
            let k128 = rust_kbest
                .get(&kmer)
                .and_then(|m| m.get(&128))
                .and_then(|m| m.get(hash));
            kv(
                "rust_walk",
                format!(
                    "label={label}\tkmer={kmer}\tpresent={}\tmissing={}\tfirst_K={first}\tk128={}",
                    walk.present,
                    walk.first_missing,
                    k128.map(|(r, s)| format!("rank={r} score={s:.9}"))
                        .unwrap_or_else(|| "absent".into())
                ),
            );
        }
        if let Some(dump) = java_dump.as_ref() {
            for (&kmer, g) in &dump.graphs {
                let walk =
                    exact_walk_generic(g.verts.len(), &g.verts, &g.outgoing, g.source, g.sink, seq);
                let mut first = "never".to_string();
                for &k in KS {
                    if let Some(rows) = dump.kbest.get(&(kmer, k)) {
                        if let Some((_, score, rank)) = rows.iter().find(|(h, _, _)| h == hash) {
                            first = format!("K={k}\trank={rank}\tscore={score:.8}");
                            break;
                        }
                    }
                }
                kv(
                    "java_walk",
                    format!(
                        "label={label}\tkmer={kmer}\tnodes={}\tedges={}\tpresent={}\tmissing={}\tfirst_K={first}",
                        g.nodes, g.edges, walk.present, walk.first_missing
                    ),
                );
            }
        }
    }

    let mut rust_present_not_sel = 0usize;
    let mut rust_absent = 0usize;
    let mut java_present_not_sel = 0usize;
    let mut java_absent = 0usize;
    if let Some(dump) = java_dump.as_ref() {
        for h in &java_only {
            let Some(seq) = java_seqs.get(h) else {
                continue;
            };
            let mut any_walk = false;
            let mut selected_k128 = false;
            for seqg in rust_graphs.values() {
                if exact_walk_seqgraph(seqg, seq).present {
                    any_walk = true;
                }
            }
            for m in rust_kbest.values() {
                if let Some(map) = m.get(&128) {
                    if map.contains_key(h) {
                        selected_k128 = true;
                    }
                }
            }
            if !any_walk {
                rust_absent += 1;
            } else if !selected_k128 {
                rust_present_not_sel += 1;
            }
            let _ = dump;
        }
        for h in &rust_only {
            let Some(seq) = rust_seq_by_hash.get(h) else {
                continue;
            };
            let mut any_walk = false;
            let mut selected_k128 = false;
            for g in dump.graphs.values() {
                if exact_walk_generic(g.verts.len(), &g.verts, &g.outgoing, g.source, g.sink, seq)
                    .present
                {
                    any_walk = true;
                }
            }
            for ((_, kcap), rows) in &dump.kbest {
                if *kcap == 128 && rows.iter().any(|(hh, _, _)| hh == h) {
                    selected_k128 = true;
                }
            }
            if !any_walk {
                java_absent += 1;
            } else if !selected_k128 {
                java_present_not_sel += 1;
            }
        }
    }
    kv(
        "accounting_private",
        format!(
            "java_only_on_rust_walk_not_kbest={rust_present_not_sel}\tjava_only_absent_from_rust_graph={rust_absent}\trust_only_on_java_walk_not_kbest={java_present_not_sel}\trust_only_absent_from_java_graph={java_absent}\tjava_only_n={}\trust_only_n={}",
            java_only.len(),
            rust_only.len()
        ),
    );

    for (&kmer, seq) in &rust_graphs {
        kv(
            "graph_size",
            format!(
                "side=rust\tkmer={kmer}\tvertices={}\tedges={}",
                seq.node_count(),
                seq.edge_count()
            ),
        );
        let k128_n = rust_kbest
            .get(&kmer)
            .and_then(|m| m.get(&128))
            .map(|m| m.len())
            .unwrap_or(0);
        kv(
            "completed_vs_k128",
            format!(
                "side=rust\tkmer={kmer}\tscored_candidates_K1024={}\tK128_unique={}",
                rust_kbest
                    .get(&kmer)
                    .and_then(|m| m.get(&1024))
                    .map(|m| m.len())
                    .unwrap_or(0),
                k128_n
            ),
        );
    }
    if let Some(dump) = java_dump.as_ref() {
        for (&kmer, g) in &dump.graphs {
            kv(
                "graph_size",
                format!(
                    "side=java\tkmer={kmer}\tvertices={}\tedges={}",
                    g.nodes, g.edges
                ),
            );
            for &k in KS {
                if let Some(rows) = dump.kbest.get(&(kmer, k)) {
                    let uniq: BTreeSet<_> = rows.iter().map(|(h, _, _)| h.clone()).collect();
                    kv(
                        "java_kbest_n",
                        format!(
                            "kmer={kmer}\tK={k}\tn={}\tunique={}",
                            rows.len(),
                            uniq.len()
                        ),
                    );
                }
            }
        }
    }
}
