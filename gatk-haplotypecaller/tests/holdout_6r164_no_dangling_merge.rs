//! 6R.164: live target proof that dangling-tail splice remains and extra hap is gone.
//! Skipped unless `HOLDOUT_6R164=1`.
//!
//! ```text
//! HOLDOUT_6R164=1 cargo test -p gatk-haplotypecaller --test holdout_6r164_no_dangling_merge -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary};
use gatk_haplotypecaller::assembly_based_caller::{
    assemble_reads_with_finalized, AssembleReadsArgs,
};
use gatk_haplotypecaller::assembly_region_finalize::{
    create_graph_reference_read, records_to_assembly_reads,
};
use gatk_haplotypecaller::read_threading_assembler::audit_threading_dangling_recovery;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use std::path::{Path, PathBuf};

const INTERVAL: &str = "2:92316200-92316580";
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 92_316_347;
const MERGED_REF: &str = "G";
const MERGED_ALT: &str = "A";
const JAVA_REF: &str = "9b4092f3b20a5de7";
const JAVA_ALT: &str = "8c9a6fc8f302f4a3";
const EXTRA: &str = "3a53b2a941bbdc43";
const JAVA_GT: &str = "1/1";
const JAVA_AD: &str = "0,3";
const JAVA_PL: &str = "135,9,0";

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R164\t{key}\t{}", value.as_ref());
}

fn fnv1a64_hex(bases: &[u8]) -> String {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in bases {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x100000001b3);
    }
    format!("{h:016x}")
}

#[test]
fn holdout_6r164_no_dangling_merge() {
    if std::env::var("HOLDOUT_6R164").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R164=1");
        return;
    }
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing P12 ref/BAM");
        return;
    }

    kv(
        "java_pin",
        "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77",
    );
    kv("variant", format!("2:{TARGET} G/A"));

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
        .expect("ActiveFull covering 2:92316347");

    let mut region = covering.clone();
    let mut ref_cache = ReferenceWindowCache::new(ref_fasta.clone(), 4);
    let mut args = AssembleReadsArgs::default();
    args.strict_java_assembly = true;
    let assembled =
        assemble_reads_with_finalized(&mut region, &dict, &mut ref_cache, &args).expect("assemble");
    let hashes: Vec<String> = assembled
        .assembly
        .haplotypes
        .iter()
        .map(|h| fnv1a64_hex(&h.bases))
        .collect();
    kv("assemble_n", hashes.len().to_string());
    for (i, (h, hash)) in assembled
        .assembly
        .haplotypes
        .iter()
        .zip(hashes.iter())
        .enumerate()
    {
        kv(
            "assemble_hap",
            format!(
                "idx={i}\thash={hash}\tisRef={}\tlen={}\tscore={}\tkmer={}",
                h.is_reference,
                h.bases.len(),
                h.score,
                h.kmer_size
            ),
        );
    }
    assert_eq!(hashes.len(), 2);
    assert!(hashes.iter().any(|h| h == JAVA_REF));
    assert!(hashes.iter().any(|h| h == JAVA_ALT));
    assert!(!hashes.iter().any(|h| h == EXTRA));

    let graph_ref = create_graph_reference_read(
        &{
            gatk_haplotypecaller::assembly_region_finalize::assembly_reference_read(
                &dict,
                &mut ReferenceWindowCache::new(ref_fasta.clone(), 4),
                covering,
            )
            .expect("ref")
        },
        covering,
        &dict,
    );
    let reads = records_to_assembly_reads(&assembled.finalized_reads);
    let mut assembler = args.assembler.clone();
    assembler.dangling_java_exact = true;
    let audit = audit_threading_dangling_recovery(&graph_ref, &reads, 35, &assembler, false, false)
        .expect("audit")
        .expect("k=35");
    kv(
        "dangling_audit_k35",
        format!(
            "edges_before={}\tedges_after={}\ttails={}/{}\theads={}/{}",
            audit.edges_before,
            audit.edges_after,
            audit.tails_recovered,
            audit.tails_attempted,
            audit.heads_recovered,
            audit.heads_attempted
        ),
    );
    assert_eq!(audit.tails_recovered, 1);
    assert!(audit.edges_after > audit.edges_before);

    let call_args = CallRegionArgs::strict_java();
    let outcome = HaplotypeCallerEngine::call_region(covering, &dict, &ref_fasta, &call_args)
        .expect("call")
        .expect("outcome");
    let in_eventmap = outcome.assembly.variation_events().iter().any(|e| {
        e.start_1based == GenomePosition::new_1based(TARGET)
            && e.ref_allele == MERGED_REF
            && e.alt_allele == MERGED_ALT
    });
    kv("eventmap_has_G_A", in_eventmap.to_string());
    kv(
        "haplotype_population",
        format!(
            "n_haps={} n_pairhmm_rows={} geno_reads={}",
            outcome.assembly.haplotypes.len(),
            outcome.read_likelihoods.len(),
            outcome.genotyping_reads.len()
        ),
    );

    let call = outcome.genotyped_calls.iter().find(|c| {
        c.event.start_1based == GenomePosition::new_1based(TARGET)
            && c.event.ref_allele == MERGED_REF
            && c.event.alt_allele == MERGED_ALT
    });
    match call {
        Some(c) => {
            let fmt = &c.genotype.format;
            kv(
                "genotyped_calls_FORMAT",
                format!(
                    "PL={:?} AD={:?} GQ={} DP={}",
                    fmt.pl_as_i32(),
                    fmt.ad_as_i32(),
                    fmt.gq.as_i32(),
                    fmt.dp.as_i32()
                ),
            );
        }
        None => kv("genotyped_calls_FORMAT", "ABSENT"),
    }

    let emitted =
        try_emit_call_region_variants(covering, &outcome, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .expect("emit");
    let rec = emitted.iter().find(|r| {
        r.position == TARGET
            && r.reference == MERGED_REF
            && r.alternate.first().map(String::as_str) == Some(MERGED_ALT)
    });
    match rec {
        Some(r) => {
            let sample = r.samples.first();
            let gt = sample
                .and_then(|s| s.gt.as_ref())
                .map(|g| {
                    g.alleles
                        .iter()
                        .map(|a| a.to_string())
                        .collect::<Vec<_>>()
                        .join("/")
                })
                .unwrap_or_default();
            let ad = sample
                .and_then(|s| s.ad.as_ref())
                .map(|v| {
                    v.iter()
                        .map(|x| x.to_string())
                        .collect::<Vec<_>>()
                        .join(",")
                })
                .unwrap_or_default();
            let pl = sample
                .and_then(|s| s.pl.as_ref())
                .map(|v| {
                    v.iter()
                        .map(|x| x.to_string())
                        .collect::<Vec<_>>()
                        .join(",")
                })
                .unwrap_or_default();
            let dp = sample
                .and_then(|s| s.dp.map(|d| d.to_string()))
                .unwrap_or_default();
            let gq = sample
                .and_then(|s| s.gq.map(|g| format!("{g:.0}")))
                .unwrap_or_default();
            kv(
                "vcf",
                format!(
                    "QUAL={:?} FILTER={:?} GT={gt} AD={ad} DP={dp} GQ={gq} PL={pl}",
                    r.quality, r.filter
                ),
            );
            kv("vcf_info", format!("{:?}", r.info));
            kv(
                "vcf_vs_java",
                format!(
                    "java=GT={JAVA_GT} AD={JAVA_AD} PL={JAVA_PL} QUAL=121.84 rust=GT={gt} AD={ad} PL={pl} QUAL={:?}",
                    r.quality
                ),
            );
        }
        None => kv("vcf", "ABSENT"),
    }

    let java_vcf = root.join("parity/reports/6r43/p12_mid_a/java.vcf");
    if java_vcf.is_file() {
        let java_text = std::fs::read_to_string(&java_vcf).unwrap_or_default();
        let mut java_keys = Vec::new();
        for line in java_text.lines() {
            if line.starts_with('#') || line.is_empty() {
                continue;
            }
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 5 {
                continue;
            }
            java_keys.push((
                f[0].to_string(),
                f[1].parse::<u64>().unwrap_or(0),
                f[3].to_string(),
                f[4].split(',').next().unwrap_or(f[4]).to_string(),
            ));
        }
        let rust_keys: Vec<_> = emitted
            .iter()
            .map(|r| {
                (
                    r.chromosome.clone(),
                    r.position,
                    r.reference.clone(),
                    r.alternate.first().cloned().unwrap_or_default(),
                )
            })
            .collect();
        let first_java_only = java_keys.iter().find(|k| !rust_keys.contains(k));
        let first_rust_only = rust_keys.iter().find(|k| !java_keys.contains(k));
        kv("p12_mid_a_java_n", java_keys.len().to_string());
        kv("p12_mid_a_rust_n", rust_keys.len().to_string());
        kv(
            "p12_mid_a_first_java_only",
            first_java_only
                .map(|(c, p, r, a)| format!("{c}:{p} {r}/{a}"))
                .unwrap_or_else(|| "none".into()),
        );
        kv(
            "p12_mid_a_first_rust_only",
            first_rust_only
                .map(|(c, p, r, a)| format!("{c}:{p} {r}/{a}"))
                .unwrap_or_else(|| "none".into()),
        );

        let mut java_fmt = std::collections::BTreeMap::new();
        for line in java_text.lines() {
            if line.starts_with('#') || line.is_empty() {
                continue;
            }
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 10 {
                continue;
            }
            let mut sample = std::collections::BTreeMap::new();
            for (k, v) in f[8].split(':').zip(f[9].split(':')) {
                sample.insert(k.to_string(), v.to_string());
            }
            let key = (
                f[0].to_string(),
                f[1].parse::<u64>().unwrap_or(0),
                f[3].to_string(),
                f[4].split(',').next().unwrap_or(f[4]).to_string(),
            );
            java_fmt.insert(
                key,
                (
                    f[5].to_string(),
                    sample.get("GT").cloned().unwrap_or_default(),
                    sample.get("AD").cloned().unwrap_or_default(),
                    sample.get("PL").cloned().unwrap_or_default(),
                ),
            );
        }
        let mut first_fmt_diff = "none".to_string();
        for rec in &emitted {
            let alt = rec.alternate.first().cloned().unwrap_or_default();
            let key = (
                rec.chromosome.clone(),
                rec.position,
                rec.reference.clone(),
                alt,
            );
            let Some((jqual, jgt, jad, jpl)) = java_fmt.get(&key) else {
                continue;
            };
            let sample = rec.samples.first();
            let gt = sample
                .and_then(|s| s.gt.as_ref())
                .map(|g| {
                    g.alleles
                        .iter()
                        .map(|a| a.to_string())
                        .collect::<Vec<_>>()
                        .join("/")
                })
                .unwrap_or_default();
            let ad = sample
                .and_then(|s| s.ad.as_ref())
                .map(|v| {
                    v.iter()
                        .map(|x| x.to_string())
                        .collect::<Vec<_>>()
                        .join(",")
                })
                .unwrap_or_default();
            let pl = sample
                .and_then(|s| s.pl.as_ref())
                .map(|v| {
                    v.iter()
                        .map(|x| x.to_string())
                        .collect::<Vec<_>>()
                        .join(",")
                })
                .unwrap_or_default();
            let rqual = rec.quality.map(|q| format!("{q:.2}")).unwrap_or_default();
            if gt != *jgt || ad != *jad || pl != *jpl {
                first_fmt_diff = format!(
                    "{}:{} {}/{} java=GT={} AD={} PL={} rust=GT={} AD={} PL={}",
                    key.0, key.1, key.2, key.3, jgt, jad, jpl, gt, ad, pl
                );
                break;
            }
            let _ = (jqual, rqual);
        }
        kv("p12_mid_a_first_format_diff", first_fmt_diff);
    }
}
