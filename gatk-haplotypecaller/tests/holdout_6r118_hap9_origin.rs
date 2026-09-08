//! 6R.118 forensic dump: origin of live hap 9 (`2ebf0f5626c63365`, 358M).
//! Skipped unless `HOLDOUT_6R118=1`. Does not reopen 6R.113–6R.117.
//!
//! ```text
//! HOLDOUT_6R118=1 P12_REFERENCE=$PWD/parity/realworld/assets/hs37d5.simple.fa \
//!   cargo test -p gatk-haplotypecaller --test holdout_6r118_hap9_origin -- --nocapture
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::{
    collect_variation_events, variation_events_for_haplotype, VariationEvent,
};
use gatk_haplotypecaller::genome_loc::GenomeLoc;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    AssemblyRegionCallDisposition, CallRegionArgs, Haplotype, HaplotypeCallerEngine,
    ReadFilterParams, WalkerTraversalConfig,
};
use serde_json::json;
use std::path::{Path, PathBuf};

const INTERVAL: &str = "2:92300000-92350000";
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 92_307_327;
const HAP9_HASH: &str = "2ebf0f5626c63365";

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
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
        .map(|c| c.to_gatk_string())
        .unwrap_or_else(|| ".".to_string())
}

fn hap_row(i: usize, h: &Haplotype) -> serde_json::Value {
    json!({
        "idx": i,
        "hash": fnv1a64_hex(&h.bases),
        "len": h.bases.len(),
        "is_reference": h.is_reference,
        "score": h.score,
        "score_is_supplement_1000": (h.score - 1000.0).abs() < 1e-6,
        "cigar": cigar_str(h),
        "align_start": h.alignment_start_hap_wrt_ref,
        "kmer_size": h.kmer_size,
        "genome_loc": h.genome_loc.map(|g| format!("{}-{}", g.start_1based(), g.end_1based())),
    })
}

fn mismatches(a: &[u8], b: &[u8], pad: u64) -> Vec<String> {
    let n = a.len().min(b.len());
    let mut out = Vec::new();
    for i in 0..n {
        if a[i].to_ascii_uppercase() != b[i].to_ascii_uppercase() {
            out.push(format!(
                "{}:{}>{}",
                pad + i as u64,
                b[i] as char,
                a[i] as char
            ));
            if out.len() >= 40 {
                break;
            }
        }
    }
    out
}

fn events_at(events: &[VariationEvent], loc: u64) -> Vec<String> {
    events
        .iter()
        .filter(|e| e.start_1based.get() == loc)
        .map(|e| format!("{}/{}", e.ref_allele, e.alt_allele))
        .collect()
}

#[test]
fn holdout_6r118_hap9_origin_dump() {
    if std::env::var("HOLDOUT_6R118").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R118=1");
        return;
    }
    let root = repo_root();
    let ref_fasta = std::env::var("P12_REFERENCE")
        .map(PathBuf::from)
        .unwrap_or_else(|_| root.join(REF_REL));
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
        .expect("ActiveFull covering 2:92307327");
    let args = CallRegionArgs::strict_java();
    let assembled =
        HaplotypeCallerEngine::call_region_assemble(covering, &dict, &ref_fasta, &args.assemble)
            .expect("assemble")
            .expect("untrimmed assembly");
    let untrimmed_rows: Vec<_> = assembled
        .haplotypes
        .iter()
        .enumerate()
        .map(|(i, h)| hap_row(i, h))
        .collect();
    let untrimmed_has_hap9 = assembled
        .haplotypes
        .iter()
        .any(|h| fnv1a64_hex(&h.bases) == HAP9_HASH);
    let untrimmed_n_358 = assembled
        .haplotypes
        .iter()
        .filter(|h| h.bases.len() == assembled.reference_bases().len())
        .count();

    let outcome = HaplotypeCallerEngine::call_region(covering, &dict, &ref_fasta, &args)
        .expect("call")
        .expect("outcome");
    let full_pad = outcome.assembly.padded_reference_start_1based();
    let full_ref = outcome.assembly.reference_bases();
    let apply_bases = outcome.assembly.apply_bases_shared();
    let apply_pad = outcome
        .assembly
        .haplotypes
        .iter()
        .find(|h| h.is_reference)
        .and_then(|h| h.genome_loc)
        .map(|g| g.start_1based())
        .unwrap_or(full_pad);
    let apply_end = outcome
        .assembly
        .haplotypes
        .iter()
        .find(|h| h.is_reference)
        .and_then(|h| h.genome_loc)
        .map(|g| g.end_1based())
        .unwrap_or(apply_pad);
    let haps = &outcome.assembly.haplotypes;
    let max_mnp = outcome.assembly.max_mnp_distance();
    let union = collect_variation_events(haps, full_ref, full_pad, "2", max_mnp);
    let hap9 = haps.iter().find(|h| fnv1a64_hex(&h.bases) == HAP9_HASH);
    let ref_hap = haps.iter().find(|h| h.is_reference);

    let mut hap9_dump = json!(null);
    if let Some(h) = hap9 {
        let mm_full = mismatches(&h.bases, full_ref, full_pad);
        let same_len_apply = h.bases.len() == apply_bases.len();
        let mm_apply = if same_len_apply {
            mismatches(&h.bases, apply_bases.as_ref(), apply_pad)
        } else {
            Vec::new()
        };
        let trim_span = GenomeLoc::new(apply_pad, apply_end);
        let trimmed = h.trim(&trim_span, true);
        let trim_hash = trimmed.as_ref().map(|t| fnv1a64_hex(&t.bases));
        let trim_len = trimmed.as_ref().map(|t| t.bases.len());
        let trim_cigar = trimmed
            .as_ref()
            .and_then(|t| t.cigar.as_ref().map(|c| c.to_gatk_string()));
        let trim_matches_live = trimmed.as_ref().is_some_and(|t| {
            haps.iter()
                .any(|x| x.bases == t.bases && x.is_reference == t.is_reference)
        });
        let untrimmed_match_idx: Vec<usize> = assembled
            .haplotypes
            .iter()
            .enumerate()
            .filter(|(_, u)| u.bases == h.bases)
            .map(|(i, _)| i)
            .collect();
        let replay = if let Some(rh) = ref_hap {
            events_at(
                &variation_events_for_haplotype(h, rh, full_ref, full_pad, max_mnp, "2"),
                TARGET,
            )
        } else {
            Vec::new()
        };
        hap9_dump = json!({
            "found": true,
            "live": hap_row(haps.iter().position(|x| fnv1a64_hex(&x.bases) == HAP9_HASH).unwrap_or(9), h),
            "hamming_vs_full_ref": if h.bases.len() == full_ref.len() {
                Some(
                    h.bases
                        .iter()
                        .zip(full_ref.iter())
                        .filter(|(a, b)| a.to_ascii_uppercase() != b.to_ascii_uppercase())
                        .count(),
                )
            } else {
                None
            },
            "n_mismatch_listed": mm_full.len(),
            "len_eq_full_ref": h.bases.len() == full_ref.len(),
            "mismatches_vs_full_ref": mm_full,
            "len_eq_apply": same_len_apply,
            "mismatches_vs_apply": mm_apply,
            "java_trim_span": format!("{apply_pad}-{apply_end}"),
            "trim_contains": h.genome_loc.is_some_and(|g| g.contains(&trim_span)),
            "after_java_trim_hash": trim_hash,
            "after_java_trim_len": trim_len,
            "after_java_trim_cigar": trim_cigar,
            "trimmed_matches_existing_live_hap": trim_matches_live,
            "present_in_untrimmed_idx": untrimmed_match_idx,
            "eventmap_at_327": replay,
        });
    }

    let dump = json!({
        "mission": "6R.118",
        "region": format!("{}:{}-{}", covering.contig, covering.start.get(), covering.end.get()),
        "full_pad": full_pad,
        "full_ref_len": full_ref.len(),
        "apply_pad": apply_pad,
        "apply_ref_len": apply_bases.len(),
        "untrimmed": {
            "hap_count": assembled.haplotypes.len(),
            "full_ref_len": assembled.reference_bases().len(),
            "padded_start": assembled.padded_reference_start_1based(),
            "has_hap9_hash": untrimmed_has_hap9,
            "n_len_eq_full_ref": untrimmed_n_358,
            "n_supplement_1000": assembled.haplotypes.iter().filter(|h| (h.score - 1000.0).abs() < 1e-6).count(),
            "per_haplotype": untrimmed_rows,
        },
        "live": {
            "hap_count": haps.len(),
            "n_supplement_1000": haps.iter().filter(|h| (h.score - 1000.0).abs() < 1e-6).count(),
            "n_len_eq_full_ref": haps.iter().filter(|h| h.bases.len() == full_ref.len()).count(),
            "eventmap_union_at_327": events_at(&union, TARGET),
            "per_haplotype": haps.iter().enumerate().map(|(i, h)| hap_row(i, h)).collect::<Vec<_>>(),
        },
        "hap9": hap9_dump,
    });
    println!("{}", serde_json::to_string_pretty(&dump).expect("json"));
}
