//! 6R.130 forensic dump: haplotype lists immediately before PairHMM at `20:29456196`.
//! Skipped unless `HOLDOUT_6R130=1`. Production change: NONE.
//!
//! ```text
//! HOLDOUT_6R130=1 cargo test -p gatk-haplotypecaller --test holdout_6r130_hap_list -- --nocapture
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::{
    begin_hap_list_observe, begin_likelihood_pipeline_observe, call_disposition,
    flatten_assembly_regions, take_hap_list_snaps, take_likelihood_pipeline_snaps,
    traverse_assembly_region_walker, AssemblyRegionCallDisposition, CallRegionArgs,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
};
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 29_456_196;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R130\t{key}\t{}", value.as_ref());
}

#[test]
fn holdout_6r130_hap_list() {
    if std::env::var("HOLDOUT_6R130").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R130=1");
        return;
    }
    let root = repo_root();
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
        .expect("ActiveFull");
    kv(
        "active",
        format!(
            "{}:{}-{}",
            covering.contig,
            covering.start.get(),
            covering.end.get()
        ),
    );
    let args = CallRegionArgs::strict_java();
    begin_hap_list_observe();
    begin_likelihood_pipeline_observe();
    let _ = HaplotypeCallerEngine::call_region(covering, &dict, &ref_fasta, &args)
        .expect("call")
        .expect("outcome");
    let snaps = take_hap_list_snaps();
    let pipe = take_likelihood_pipeline_snaps();
    kv("n_snaps", snaps.len().to_string());
    for s in &snaps {
        let mut uniq = std::collections::BTreeSet::new();
        for c in &s.columns {
            uniq.insert(format!("{:016x}", c.fnv1a));
        }
        kv(
            "snap",
            format!("stage={}\tn_cols={}\tunique={}", s.stage, s.n, uniq.len()),
        );
        for c in &s.columns {
            kv(
                "hap",
                format!(
                    "stage={}\tidx={}\thash={:016x}\tisRef={}\tlen={}",
                    s.stage, c.index, c.fnv1a, c.is_reference, c.len
                ),
            );
        }
    }
    for s in &pipe {
        kv(
            "pipe_snap",
            format!(
                "seq={}\tstage={}\tn_reads={}\tn_haps={}\tn_ll={}",
                s.seq, s.stage, s.n_reads, s.n_haps, s.n_ll_entries
            ),
        );
    }
}
