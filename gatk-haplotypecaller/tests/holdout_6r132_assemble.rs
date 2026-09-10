//! 6R.132 forensic dump: untrimmed assemble population vs Java window at `20:29456196`.
//! Skipped unless `HOLDOUT_6R132=1`. Production change: NONE.
//!
//! ```text
//! HOLDOUT_6R132=1 cargo test -p gatk-haplotypecaller --test holdout_6r132_assemble -- --nocapture
//! ```

use gatk_core::reference::{parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary};
use gatk_haplotypecaller::{
    assemble_reads, call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    AssemblyRegionCallDisposition, CallRegionArgs, GenomeLoc, GenomePosition, Haplotype,
    ReadFilterParams, WalkerTraversalConfig,
};
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const JAVA_DUMP_REL: &str = "parity/giab/runs/local-pairhmm-diff/6r131_java.txt";
const TARGET: u64 = 29_456_196;
const JAVA_ACTIVE_START: u64 = 29_455_995;
const JAVA_ACTIVE_END: u64 = 29_456_293;
const JAVA_PAD_START: u64 = 29_455_895;
const JAVA_PAD_END: u64 = 29_456_393;
const WORST_PARENTS: &[&str] = &[
    "a24d349c1a423e3b",
    "275eab049a12c68f",
    "39c53d4bf0725810",
    "c7953f2f7ade1e75",
];

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R132\t{key}\t{}", value.as_ref());
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

fn load_java_untrimmed(path: &Path) -> BTreeSet<String> {
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

fn dump_haps(stage: &str, haps: &[Haplotype]) {
    let uniq = unique_of(haps);
    let n_ref = haps.iter().filter(|h| h.is_reference).count();
    let mut counts = std::collections::BTreeMap::<String, usize>::new();
    for h in haps {
        *counts.entry(hex(fnv1a64(&h.bases))).or_insert(0) += 1;
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
        let loc = h
            .genome_loc
            .map(|g| format!("{}-{}", g.start_1based(), g.end_1based()))
            .unwrap_or_else(|| ".".into());
        let cigar = h
            .cigar
            .as_ref()
            .map(|c| c.to_gatk_string())
            .unwrap_or_else(|| ".".into());
        kv(
            "hap",
            format!(
                "stage={stage}\tidx={i}\thash={}\tisRef={}\tlen={}\tscore={}\tloc={loc}\tcigar={cigar}",
                hex(fnv1a64(&h.bases)),
                h.is_reference,
                h.bases.len(),
                h.score
            ),
        );
    }
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

fn slice_to_window(haps: &[Haplotype], start: u64, end: u64) -> Vec<Haplotype> {
    let span = GenomeLoc::new(start, end);
    haps.iter().filter_map(|h| h.trim(&span, true)).collect()
}

#[test]
fn holdout_6r132_assemble() {
    if std::env::var("HOLDOUT_6R132").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R132=1");
        return;
    }
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    let java_dump = root.join(JAVA_DUMP_REL);
    let java_untrimmed = load_java_untrimmed(&java_dump);
    kv("java_untrimmed_loaded", java_untrimmed.len().to_string());

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
    for r in &regions {
        kv(
            "region",
            format!(
                "{}:{}-{}\tactive={}\tdisp={:?}\text={}-{}",
                r.contig,
                r.start.get(),
                r.end.get(),
                r.is_active,
                call_disposition(r),
                r.extended_start.get(),
                r.extended_end.get()
            ),
        );
    }
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
    kv(
        "covering",
        format!(
            "active={}-{}\text={}-{}\tn_reads={}",
            covering.start.get(),
            covering.end.get(),
            covering.extended_start.get(),
            covering.extended_end.get(),
            covering.reads.len()
        ),
    );
    let mut qnames: Vec<String> = covering
        .reads
        .iter()
        .map(|r| String::from_utf8_lossy(r.as_ref().qname()).into_owned())
        .collect();
    qnames.sort();
    kv("covering_qnames_n", qnames.len().to_string());
    for q in &qnames {
        kv("qname", q.as_str());
    }

    let args = CallRegionArgs::strict_java();
    let mut ref_cache = ReferenceWindowCache::new(ref_fasta.clone(), 4);
    let assembly =
        assemble_reads(&covering, &dict, &mut ref_cache, &args.assemble).expect("assemble");
    let graph_ref = assembly.reference_bases();
    kv(
        "graph_ref",
        format!(
            "hash={}\tlen={}\tpad_start={}",
            hex(fnv1a64(graph_ref)),
            graph_ref.len(),
            assembly.padded_reference_start_1based()
        ),
    );
    dump_haps("assemble", &assembly.haplotypes);
    let rust_untrimmed = unique_of(&assembly.haplotypes);
    compare(
        "assemble_vs_java_untrimmed",
        &rust_untrimmed,
        &java_untrimmed,
    );

    let sliced = slice_to_window(&assembly.haplotypes, JAVA_PAD_START, JAVA_PAD_END);
    kv(
        "slice_to_java_pad",
        format!(
            "n_ok={}\tn_drop={}",
            sliced.len(),
            assembly.haplotypes.len().saturating_sub(sliced.len())
        ),
    );
    dump_haps("slice_java_pad", &sliced);
    let rust_sliced = unique_of(&sliced);
    compare(
        "slice_java_pad_vs_java_untrimmed",
        &rust_sliced,
        &java_untrimmed,
    );

    for p in WORST_PARENTS {
        let orig = assembly
            .haplotypes
            .iter()
            .find(|h| hex(fnv1a64(&h.bases)) == *p);
        match orig {
            None => kv("worst_parent", format!("hash={p}\tpresent=false")),
            Some(h) => {
                let sliced_h = h.trim(&GenomeLoc::new(JAVA_PAD_START, JAVA_PAD_END), true);
                let sh = sliced_h.as_ref().map(|t| hex(fnv1a64(&t.bases)));
                kv(
                    "worst_parent",
                    format!(
                        "hash={p}\tpresent=true\tlen={}\tcigar={}\tslice={}\tin_java_untrimmed={}",
                        h.bases.len(),
                        h.cigar
                            .as_ref()
                            .map(|c| c.to_gatk_string())
                            .unwrap_or_else(|| ".".into()),
                        sh.as_deref().unwrap_or("dropped"),
                        sh.as_ref()
                            .map(|s| java_untrimmed.contains(s))
                            .unwrap_or(false)
                    ),
                );
            }
        }
    }

    let mut java_bounds = covering.clone();
    java_bounds.start = GenomePosition::new_1based(JAVA_ACTIVE_START);
    java_bounds.end = GenomePosition::new_1based(JAVA_ACTIVE_END);
    java_bounds.extended_start = GenomePosition::new_1based(JAVA_PAD_START);
    java_bounds.extended_end = GenomePosition::new_1based(JAVA_PAD_END);
    kv(
        "counterfactual_region",
        format!(
            "active={}-{}\text={}-{}",
            java_bounds.start.get(),
            java_bounds.end.get(),
            java_bounds.extended_start.get(),
            java_bounds.extended_end.get()
        ),
    );
    let mut java_owned = java_bounds.clone();
    let cf_fin = gatk_haplotypecaller::assemble_reads_with_finalized(
        &mut java_owned,
        &dict,
        &mut ref_cache,
        &args.assemble,
    )
    .expect("counterfactual assemble_with_finalized");
    kv(
        "counterfactual_finalized_n",
        cf_fin.finalized_reads.len().to_string(),
    );
    let cf = cf_fin.assembly;
    kv(
        "counterfactual_graph_ref",
        format!(
            "hash={}\tlen={}\tpad_start={}",
            hex(fnv1a64(cf.reference_bases())),
            cf.reference_bases().len(),
            cf.padded_reference_start_1based()
        ),
    );
    dump_haps("counterfactual_java_window", &cf.haplotypes);
    let cf_set = unique_of(&cf.haplotypes);
    compare(
        "counterfactual_java_window_vs_java_untrimmed",
        &cf_set,
        &java_untrimmed,
    );
    compare(
        "counterfactual_vs_native_assemble",
        &cf_set,
        &rust_untrimmed,
    );
}
