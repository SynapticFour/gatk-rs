//! 6R.143: PairHMM haplotype/read inputs vs VariantContext construction.
//! Diagnostic-only `unbounded_diagnostic`. Production k-best policy unchanged.
//! Skipped unless `HOLDOUT_6R143=1`.
//!
//! ```text
//! HOLDOUT_6R143=1 GATK_RS_EXPERIMENTAL_KBEST_POLICY=unbounded_diagnostic \
//!   cargo test -p gatk-haplotypecaller --test holdout_6r143_likelihood_input -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary};
use gatk_haplotypecaller::event_map::{
    build_per_haplotype_variation_events, variation_events_at_position_from_cache, VariationEvent,
};
use gatk_haplotypecaller::{
    assemble_reads_with_finalized, begin_hap_list_observe, begin_likelihood_pipeline_observe,
    call_disposition, flatten_assembly_regions, take_hap_list_snaps,
    take_likelihood_pipeline_cells, take_likelihood_pipeline_snaps,
    traverse_assembly_region_walker, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    Haplotype, HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
};
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const JAVA_LL_DUMP_REL: &str = "parity/giab/runs/local-pairhmm-diff/6r143_java.txt";
const JAVA_EM_DUMP_REL: &str = "parity/giab/runs/local-pairhmm-diff/6r142_java.txt";
const TARGET: u64 = 29_456_196;
const JAVA_ACTIVE_START: u64 = 29_455_995;
const JAVA_ACTIVE_END: u64 = 29_456_293;
const JAVA_PAD_START: u64 = 29_455_895;
const JAVA_PAD_END: u64 = 29_456_393;
const JAVA_TRIM_START: u64 = 29_455_977;
const JAVA_TRIM_END: u64 = 29_456_301;
const EXPECTED_REF_HASH: &str = "659741f99b7f78a5";
const EXPECTED_TRIM_REF_HASH: &str = "9af99c7ab7cfc3ee";
const EXPECTED_FINALIZED: usize = 379;
const MAX_MNP: usize = 0;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R143\t{key}\t{}", value.as_ref());
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

#[derive(Clone, Debug)]
struct JavaHap {
    idx: usize,
    hash: String,
    is_ref: bool,
    align_start: usize,
    cigar: String,
}

#[derive(Clone, Debug)]
struct JavaRead {
    qname: String,
    flags: String,
    bq: String,
    iq: String,
    dq: String,
    gcp: String,
    bases: String,
}

fn parse_kv_fields(parts: &[&str]) -> BTreeMap<String, String> {
    let mut m = BTreeMap::new();
    for tok in parts {
        if let Some((k, v)) = tok.split_once('=') {
            m.insert(k.to_string(), v.to_string());
        }
    }
    m
}

fn load_java_haps(path: &Path) -> Vec<JavaHap> {
    let text = std::fs::read_to_string(path).expect("6r143 java dump");
    let mut out = Vec::new();
    for line in text.lines() {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 4 || parts[0] != "6R143" || parts[1] != "hap" {
            continue;
        }
        if !parts[2].chars().all(|c| c.is_ascii_digit()) {
            continue;
        }
        let kv = parse_kv_fields(&parts[4..]);
        out.push(JavaHap {
            idx: parts[2].parse().expect("hap idx"),
            hash: parts[3].to_string(),
            is_ref: kv.get("isRef").map(|s| s == "true").unwrap_or(false),
            align_start: kv
                .get("alignStart")
                .and_then(|s| s.parse().ok())
                .unwrap_or(0),
            cigar: kv.get("cigar").cloned().unwrap_or_else(|| ".".into()),
        });
    }
    out
}

fn load_java_reads(path: &Path, kind: &str) -> Vec<JavaRead> {
    let text = std::fs::read_to_string(path).expect("6r143 java dump");
    let mut out = Vec::new();
    for line in text.lines() {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 4 || parts[0] != "6R143" || parts[1] != kind {
            continue;
        }
        if !parts[2].chars().all(|c| c.is_ascii_digit()) {
            continue;
        }
        let kv = parse_kv_fields(&parts[3..]);
        out.push(JavaRead {
            qname: kv.get("qname").cloned().unwrap_or_default(),
            flags: kv.get("flags").cloned().unwrap_or_default(),
            bq: kv.get("bq").cloned().unwrap_or_default(),
            iq: kv.get("iq").cloned().unwrap_or_default(),
            dq: kv.get("dq").cloned().unwrap_or_default(),
            gcp: kv.get("gcp").cloned().unwrap_or_default(),
            bases: kv.get("bases").cloned().unwrap_or_default(),
        });
    }
    out
}

fn load_java_em_union_at_loc(path: &Path, loc: u64) -> BTreeSet<(u64, String, String)> {
    let text = std::fs::read_to_string(path).expect("6r142 java dump");
    let mut out = BTreeSet::new();
    for line in text.lines() {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 3 || parts[0] != "6R142" || parts[1] != "event" {
            continue;
        }
        let kv = parse_kv_fields(&parts[2..]);
        if kv.get("stage").map(|s| s.as_str()) != Some("trimmed") {
            continue;
        }
        let start: u64 = kv.get("start").and_then(|s| s.parse().ok()).unwrap_or(0);
        let end: u64 = kv.get("end").and_then(|s| s.parse().ok()).unwrap_or(start);
        if start <= loc && loc <= end {
            out.insert((
                start,
                kv.get("ref").cloned().unwrap_or_default(),
                kv.get("alt").cloned().unwrap_or_default(),
            ));
        }
    }
    out
}

fn load_java_hap_event_contrib(path: &Path, loc: u64) -> (BTreeSet<String>, usize, usize) {
    let text = std::fs::read_to_string(path).expect("6r142 java dump");
    let mut hashes = BTreeSet::new();
    let mut n_ref_events = 0usize;
    let mut n_haps = 0usize;
    let mut ref_hash = String::new();
    for line in text.lines() {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 3 || parts[0] != "6R142" {
            continue;
        }
        if parts[1] == "hap" {
            let kv = parse_kv_fields(&parts[2..]);
            if kv.get("stage").map(|s| s.as_str()) != Some("trimmed") {
                continue;
            }
            n_haps += 1;
            if kv.get("isRef").map(|s| s.as_str()) == Some("true") {
                ref_hash = kv.get("hash").cloned().unwrap_or_default();
            }
        }
        if parts[1] == "event" {
            let kv = parse_kv_fields(&parts[2..]);
            if kv.get("stage").map(|s| s.as_str()) != Some("trimmed") {
                continue;
            }
            let start: u64 = kv.get("start").and_then(|s| s.parse().ok()).unwrap_or(0);
            if start == loc {
                if let Some(h) = kv.get("hash") {
                    hashes.insert(h.clone());
                    if *h == ref_hash {
                        n_ref_events += 1;
                    }
                }
            }
        }
    }
    (hashes, n_ref_events, n_haps)
}

fn ve_key(e: &VariationEvent) -> (u64, String, String) {
    (
        e.start_1based.get(),
        e.ref_allele.clone(),
        e.alt_allele.clone(),
    )
}

fn compare_seq(label: &str, rust: &[String], java: &[String]) -> (usize, usize, usize, bool) {
    let rs: BTreeSet<_> = rust.iter().cloned().collect();
    let js: BTreeSet<_> = java.iter().cloned().collect();
    let common = rs.intersection(&js).count();
    let java_only = js.difference(&rs).count();
    let rust_only = rs.difference(&js).count();
    let order = rust == java;
    kv(
        "compare",
        format!(
            "label={label}\trust={}\tjava={}\tCOMMON={common}\tJAVA_ONLY={java_only}\tRUST_ONLY={rust_only}\torder_identical={order}",
            rust.len(),
            java.len()
        ),
    );
    for h in js.difference(&rs).take(8) {
        kv("java_only", format!("label={label}\t{h}"));
    }
    for h in rs.difference(&js).take(8) {
        kv("rust_only", format!("label={label}\t{h}"));
    }
    (common, java_only, rust_only, order)
}

fn parse_kernel_dump(path: &Path, n_haps: usize) -> Vec<(String, String, String, String, String)> {
    let text = std::fs::read_to_string(path).unwrap_or_default();
    let mut rows = Vec::new();
    let mut i = 0usize;
    for line in text.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let mut it = line.split(' ');
        let _hap = it.next();
        let read = it.next().unwrap_or("");
        let bq = it.next().unwrap_or("");
        let iq = it.next().unwrap_or("");
        let dq = it.next().unwrap_or("");
        let gcp = it.next().unwrap_or("");
        if n_haps == 0 || i % n_haps == 0 {
            rows.push((
                read.to_string(),
                bq.to_string(),
                iq.to_string(),
                dq.to_string(),
                gcp.to_string(),
            ));
        }
        i += 1;
    }
    rows
}

#[test]
fn holdout_6r143_likelihood_input() {
    if std::env::var("HOLDOUT_6R143").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R143=1");
        return;
    }
    let _policy = EnvGuard::set("GATK_RS_EXPERIMENTAL_KBEST_POLICY", "unbounded_diagnostic");
    kv(
        "parity_mode",
        "unbounded_diagnostic\tnote=diagnostic_not_production_default",
    );

    let root = repo_root();
    let java_haps = load_java_haps(&root.join(JAVA_LL_DUMP_REL));
    let java_orig = load_java_reads(&root.join(JAVA_LL_DUMP_REL), "orig");
    let java_proc = load_java_reads(&root.join(JAVA_LL_DUMP_REL), "proc");
    assert_eq!(java_haps.len(), 25, "Java PairHMM hap_count");
    kv("java_hap_count", java_haps.len().to_string());
    kv(
        "java_unique",
        java_haps
            .iter()
            .map(|h| h.hash.as_str())
            .collect::<BTreeSet<_>>()
            .len()
            .to_string(),
    );
    kv("java_orig_n", java_orig.len().to_string());
    kv("java_proc_n", java_proc.len().to_string());
    kv(
        "java_eventmap_at_pairhmm",
        "null\tgetVariantContextsFromActiveHaplotypes=NPE_not_invoked_before_kernel",
    );

    let java_order: Vec<String> = java_haps.iter().map(|h| h.hash.clone()).collect();
    let java_ref = java_haps.iter().find(|h| h.is_ref).expect("java ref");
    kv(
        "java_ref",
        format!(
            "idx={}\thash={}\talignStart={}\tcigar={}",
            java_ref.idx, java_ref.hash, java_ref.align_start, java_ref.cigar
        ),
    );
    assert_eq!(java_ref.hash, EXPECTED_TRIM_REF_HASH);

    let dump_path = std::env::temp_dir().join("6r143_rust_pairhmm_inputs.txt");
    let _dump_env = EnvGuard::set(
        "GATK_RS_PAIRHMM_INPUT_DUMP",
        dump_path.to_str().expect("utf8 dump path"),
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
    begin_hap_list_observe();
    begin_likelihood_pipeline_observe();
    let _ = HaplotypeCallerEngine::call_region(&java_bounds, &dict, &ref_fasta, &args)
        .expect("call")
        .expect("outcome");
    let snaps = take_hap_list_snaps();
    let pipe = take_likelihood_pipeline_snaps();
    let cells = take_likelihood_pipeline_cells();

    for s in &snaps {
        let uniq: BTreeSet<_> = s.columns.iter().map(|c| hex(c.fnv1a)).collect();
        kv(
            "snap",
            format!("stage={}\tn_cols={}\tunique={}", s.stage, s.n, uniq.len()),
        );
    }
    let pairhmm = snaps
        .iter()
        .find(|s| s.stage == "pairhmm_input")
        .expect("pairhmm_input snap");
    let after_trim = snaps
        .iter()
        .find(|s| s.stage == "after_trim_to")
        .expect("after_trim_to");
    let after_preserve = snaps
        .iter()
        .find(|s| s.stage == "after_preserve_untrimmed")
        .expect("after_preserve");
    kv(
        "preserve_delta",
        format!(
            "after_trim_to={}\tafter_preserve={}\tpairhmm_input={}",
            after_trim.n, after_preserve.n, pairhmm.n
        ),
    );
    assert_eq!(
        after_trim.n, after_preserve.n,
        "preserve_untrimmed must not re-insert REF bases as a second PairHMM column"
    );
    assert_eq!(after_preserve.n, pairhmm.n);

    let rust_order: Vec<String> = pairhmm.columns.iter().map(|c| hex(c.fnv1a)).collect();
    let (common, java_only, rust_only, order) =
        compare_seq("pairhmm_hap_seq", &rust_order, &java_order);
    kv(
        "hap_meta",
        format!(
            "rust_n_ref={}\tjava_n_ref={}\trust_cigars={}\tjava_cigars={}",
            pairhmm.columns.iter().filter(|c| c.is_reference).count(),
            java_haps.iter().filter(|h| h.is_ref).count(),
            pairhmm
                .columns
                .iter()
                .map(|c| c.cigar.as_str())
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect::<Vec<_>>()
                .join(","),
            java_haps
                .iter()
                .map(|h| h.cigar.as_str())
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect::<Vec<_>>()
                .join(",")
        ),
    );
    let ref_flag_match = pairhmm
        .columns
        .iter()
        .zip(java_haps.iter())
        .all(|(c, j)| c.is_reference == j.is_ref);
    let cigar_match = pairhmm
        .columns
        .iter()
        .zip(java_haps.iter())
        .all(|(c, j)| c.cigar == j.cigar);
    let loc_match = pairhmm
        .columns
        .iter()
        .all(|c| c.loc_start == JAVA_TRIM_START && c.loc_end == JAVA_TRIM_END);
    kv(
        "hap_fields",
        format!("isRef_identical={ref_flag_match}\tcigar_identical={cigar_match}\tloc_identical={loc_match}\talignStart_java=582\talignStart_rust=graph_pad_offset_representation_only")
    );
    assert_eq!(common, 25);
    assert_eq!(java_only, 0);
    assert_eq!(rust_only, 0);
    assert!(order, "PairHMM haplotype order must match Java");
    assert!(ref_flag_match);
    assert!(cigar_match);

    let mut rust_reads: Vec<(String, u16)> = Vec::new();
    if let Some(first) = cells
        .iter()
        .min_by_key(|c| (c.seq, c.read_index, c.hap_index))
    {
        let seq0 = first.seq;
        let mut by_idx: BTreeMap<usize, (String, u16)> = BTreeMap::new();
        for c in &cells {
            if c.seq == seq0 {
                by_idx
                    .entry(c.read_index)
                    .or_insert((c.qname.clone(), c.flags));
            }
        }
        rust_reads = by_idx.into_values().collect();
    }
    kv("rust_pairhmm_reads", rust_reads.len().to_string());
    for s in &pipe {
        kv(
            "pipe_snap",
            format!(
                "seq={}\tstage={}\tn_reads={}\tn_haps={}",
                s.seq, s.stage, s.n_reads, s.n_haps
            ),
        );
    }
    let java_qnames: Vec<String> = java_orig.iter().map(|r| r.qname.clone()).collect();
    let rust_qnames: Vec<String> = rust_reads.iter().map(|(q, _)| q.clone()).collect();
    let (rq_common, rq_j, rq_r, rq_order) =
        compare_seq("pairhmm_orig_qname", &rust_qnames, &java_qnames);
    let java_keys: Vec<String> = java_orig
        .iter()
        .map(|r| format!("{}:{}", r.qname, r.flags))
        .collect();
    let rust_keys: Vec<String> = rust_reads.iter().map(|(q, f)| format!("{q}:{f}")).collect();
    let (rk_common, rk_j, rk_r, rk_order) =
        compare_seq("pairhmm_orig_qname_flags", &rust_keys, &java_keys);
    if !rk_order {
        for i in 0..rust_keys.len().min(java_keys.len()) {
            if rust_keys[i] != java_keys[i] {
                kv(
                    "orig_order_first_mismatch",
                    format!("idx={i}\trust={}\tjava={}", rust_keys[i], java_keys[i]),
                );
                break;
            }
        }
    }
    kv(
        "orig_order",
        format!(
            "qname_set_COMMON={rq_common}\tJAVA_ONLY={rq_j}\tRUST_ONLY={rq_r}\tqname_order={rq_order}\tkey_COMMON={rk_common}\tkey_JAVA_ONLY={rk_j}\tkey_RUST_ONLY={rk_r}\tkey_order={rk_order}"
        ),
    );

    let kernel = parse_kernel_dump(&dump_path, pairhmm.n);
    kv("rust_kernel_rows", kernel.len().to_string());
    let mut java_queue: BTreeMap<(String, String), Vec<usize>> = BTreeMap::new();
    for (i, r) in java_orig.iter().enumerate() {
        java_queue
            .entry((r.qname.clone(), r.flags.clone()))
            .or_default()
            .push(i);
    }
    let mut kernel_eq = 0usize;
    let mut kernel_ne = 0usize;
    let mut kernel_unmatched = 0usize;
    for (i, (q, f)) in rust_reads.iter().enumerate() {
        let key = (q.clone(), f.to_string());
        let Some(idx) = java_queue.get_mut(&key).and_then(|v| {
            if v.is_empty() {
                None
            } else {
                Some(v.remove(0))
            }
        }) else {
            kernel_unmatched += 1;
            continue;
        };
        let jr = &java_proc[idx];
        let Some((bases, bq, iq, dq, gcp)) = kernel.get(i) else {
            kernel_ne += 1;
            continue;
        };
        if bases == &jr.bases && bq == &jr.bq && iq == &jr.iq && dq == &jr.dq && gcp == &jr.gcp {
            kernel_eq += 1;
        } else {
            kernel_ne += 1;
            if kernel_ne <= 4 {
                kv(
                    "kernel_mismatch",
                    format!(
                        "rust_idx={i}\tjava_idx={idx}\tqname={q}\tbases_eq={}\tbq_eq={}\tiq_eq={}\tdq_eq={}\tgcp_eq={}",
                        bases == &jr.bases,
                        bq == &jr.bq,
                        iq == &jr.iq,
                        dq == &jr.dq,
                        gcp == &jr.gcp
                    ),
                );
            }
        }
    }
    kv(
        "kernel_planes_aligned",
        format!(
            "equal={kernel_eq}\tne={kernel_ne}\tunmatched={kernel_unmatched}\tjava_proc={}\trust_kernel={}",
            java_proc.len(),
            kernel.len()
        ),
    );

    let mut ref_cache = ReferenceWindowCache::new(ref_fasta.clone(), 4);
    let mut owned = java_bounds.clone();
    let assembled =
        assemble_reads_with_finalized(&mut owned, &dict, &mut ref_cache, &args.assemble)
            .expect("assemble");
    assert_eq!(
        hex(fnv1a64(assembled.assembly.reference_bases())),
        EXPECTED_REF_HASH
    );
    assert_eq!(assembled.finalized_reads.len(), EXPECTED_FINALIZED);
    let mut trim_region = java_bounds.clone();
    trim_region.extended_start = GenomePosition::new_1based(JAVA_TRIM_START);
    trim_region.extended_end = GenomePosition::new_1based(JAVA_TRIM_END);
    let trimmed = assembled.assembly.trim_to(&trim_region).expect("trim_to");
    let rust_trim_order: Vec<String> = trimmed.haplotypes.iter().map(hap_hash).collect();
    assert_eq!(rust_trim_order, java_order);

    let (graph_ref, graph_pad) = trimmed.event_map_reference();
    let cache = build_per_haplotype_variation_events(
        &trimmed.haplotypes,
        graph_ref,
        graph_pad,
        MAX_MNP,
        "20",
    );
    let vc_span = variation_events_at_position_from_cache(&cache, TARGET, true);
    let vc_no_span = variation_events_at_position_from_cache(&cache, TARGET, false);
    kv(
        "vc_at_loc",
        format!(
            "span_n={}\tno_span_n={}\tkeys={}",
            vc_span.len(),
            vc_no_span.len(),
            vc_span
                .iter()
                .map(|e| format!("{}:{}/{}", e.start_1based.get(), e.ref_allele, e.alt_allele))
                .collect::<Vec<_>>()
                .join(",")
        ),
    );
    let java_vc = load_java_em_union_at_loc(&root.join(JAVA_EM_DUMP_REL), TARGET);
    let rust_vc: BTreeSet<_> = vc_span.iter().map(ve_key).collect();
    kv(
        "vc_compare",
        format!(
            "rust={}\tjava={}\tCOMMON={}\tJAVA_ONLY={}\tRUST_ONLY={}",
            rust_vc.len(),
            java_vc.len(),
            rust_vc.intersection(&java_vc).count(),
            java_vc.difference(&rust_vc).count(),
            rust_vc.difference(&java_vc).count()
        ),
    );
    assert_eq!(rust_vc, java_vc);
    assert_eq!(vc_span.len(), vc_no_span.len());
    let (contrib, n_ref_events, n_haps) =
        load_java_hap_event_contrib(&root.join(JAVA_EM_DUMP_REL), TARGET);
    kv(
        "vc_contrib",
        format!(
            "haps_with_event_at_loc={}\tref_hap_events={n_ref_events}\ttrimmed_haps={n_haps}",
            contrib.len()
        ),
    );
    assert_eq!(n_haps, 25);
    assert_eq!(n_ref_events, 0);

    kv(
        "summary",
        format!(
            "haplotypes_java={}/25\thaplotypes_rust={}/25\tvariant_contexts_equal={}\tpairhmm_haplotype_inputs_equal={}\tpairhmm_read_identity_equal={}\tpairhmm_read_order_identical={}\tkernel_planes_aligned_equal={}\tfirst_divergence={}",
            java_haps.len(),
            pairhmm.n,
            rust_vc == java_vc,
            common == 25 && java_only == 0 && rust_only == 0 && order && ref_flag_match,
            rk_j == 0 && rk_r == 0 && rust_reads.len() == java_orig.len(),
            rk_order,
            kernel_eq == java_proc.len() && kernel_ne == 0 && kernel_unmatched == 0,
            if common != 25 || java_only != 0 || rust_only != 0 || !order {
                "pairhmm_haplotype_list"
            } else if rust_vc != java_vc {
                "variant_context_construction"
            } else if rk_j != 0 || rk_r != 0 {
                "pairhmm_read_population"
            } else if kernel_ne != 0 || kernel_unmatched != 0 {
                "pairhmm_kernel_planes"
            } else if !rk_order {
                "pairhmm_read_order_noncausal_for_kernel_planes"
            } else {
                "NONE"
            }
        ),
    );
    assert_eq!(rk_j, 0, "JAVA_ONLY orig QNAME+flags");
    assert_eq!(rk_r, 0, "RUST_ONLY orig QNAME+flags");
    assert_eq!(rust_reads.len(), java_orig.len());
    assert_eq!(kernel_unmatched, 0);
    assert_eq!(
        kernel_ne, 0,
        "PairHMM kernel BQ/IQ/DQ/GCP/bases aligned by QNAME+flags"
    );
    assert_eq!(kernel_eq, java_proc.len());
}
