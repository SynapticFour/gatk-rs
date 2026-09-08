//! 6R.116 forensic dump: live haplotype population / CIGAR representation at
//! `2:92307327`. Skipped unless `HOLDOUT_6R116=1`. Does not reopen 6R.114 / 6R.115.
//!
//! ```text
//! HOLDOUT_6R116=1 P12_REFERENCE=$PWD/parity/realworld/assets/hs37d5.simple.fa \
//!   cargo test -p gatk-haplotypecaller --test holdout_6r116_hap_population -- --nocapture
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::{
    build_per_haplotype_variation_events, variation_events_for_haplotype, VariationEvent,
};
use gatk_haplotypecaller::haplotype_cigar::calculate_haplotype_cigar_java_padded_sw;
use gatk_haplotypecaller::hc_allele_mapping::hap_base_at_ref_locus;
use gatk_haplotypecaller::smith_waterman::{SwOverhangStrategy, SwParameters};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    AssemblyRegionCallDisposition, CallRegionArgs, Haplotype, HaplotypeCallerEngine,
    ReadFilterParams, WalkerTraversalConfig,
};
use serde_json::{json, Value};
use std::path::{Path, PathBuf};

const INTERVAL: &str = "2:92300000-92350000";
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 92_307_327;
const MOTIF_LO: u64 = 92_307_324;
const MOTIF_HI: u64 = 92_307_330;

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

fn cigar_has_indel(h: &Haplotype) -> bool {
    h.cigar
        .as_ref()
        .is_some_and(|c| c.elements.iter().any(|e| e.operator.is_indel()))
}

fn cigar_is_all_m(h: &Haplotype) -> bool {
    h.cigar.as_ref().is_some_and(|c| {
        c.elements.len() == 1
            && matches!(
                c.elements[0].operator,
                gatk_haplotypecaller::cigar::CigarOperator::Match
            )
    })
}

fn linear_motif(h: &Haplotype, pad: u64) -> String {
    let mut s = String::new();
    for loc in MOTIF_LO..=MOTIF_HI {
        match hap_base_at_ref_locus(h, pad, loc) {
            Some(b) => s.push(b as char),
            None => s.push('.'),
        }
    }
    s
}

fn hamming(a: &[u8], b: &[u8]) -> Option<usize> {
    if a.len() != b.len() {
        return None;
    }
    Some(
        a.iter()
            .zip(b)
            .filter(|(x, y)| x.to_ascii_uppercase() != y.to_ascii_uppercase())
            .count(),
    )
}

fn events_at(events: &[VariationEvent], loc: u64) -> Vec<String> {
    events
        .iter()
        .filter(|e| e.start_1based.get() == loc)
        .map(|e| format!("{}/{}", e.ref_allele, e.alt_allele))
        .collect()
}

#[test]
fn holdout_6r116_hap_population_dump() {
    if std::env::var("HOLDOUT_6R116").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R116=1");
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
    let haps = &outcome.assembly.haplotypes;
    let max_mnp = outcome.assembly.max_mnp_distance();
    let hap_events = build_per_haplotype_variation_events(haps, full_ref, full_pad, max_mnp, "2");
    let sw = SwParameters::gatk_haplotype_to_reference();
    let coupled = haps
        .iter()
        .find(|h| !h.is_reference && cigar_has_indel(h) && (h.score - 1000.0).abs() < 1e-6);
    let ref_hap = haps.iter().find(|h| h.is_reference).expect("ref hap");

    let mut rows = Vec::new();
    for (i, h) in haps.iter().enumerate() {
        let motif = linear_motif(h, full_pad);
        let motif_is_tatg = motif.starts_with("TATG") || motif == "T..ATGA";
        let vs_coupled = coupled.and_then(|c| hamming(&h.bases, &c.bases));
        let vs_apply = hamming(&h.bases, apply_bases.as_ref());
        let vs_full = hamming(&h.bases, full_ref);
        let java_sw_apply = calculate_haplotype_cigar_java_padded_sw(
            apply_bases.as_ref(),
            &h.bases,
            &sw,
            SwOverhangStrategy::Indel,
        )
        .map(|r| r.cigar.to_gatk_string());
        let java_sw_soft = calculate_haplotype_cigar_java_padded_sw(
            apply_bases.as_ref(),
            &h.bases,
            &sw,
            SwOverhangStrategy::SoftClip,
        )
        .map(|r| r.cigar.to_gatk_string());
        let replay = events_at(hap_events.events_for(i), TARGET);
        let per_hap_replay =
            variation_events_for_haplotype(h, ref_hap, full_ref, full_pad, max_mnp, "2");
        rows.push(json!({
            "idx": i,
            "hash": fnv1a64_hex(&h.bases),
            "len": h.bases.len(),
            "is_reference": h.is_reference,
            "score": h.score,
            "score_is_supplement_1000": (h.score - 1000.0).abs() < 1e-6,
            "cigar": cigar_str(h),
            "cigar_has_indel": cigar_has_indel(h),
            "cigar_is_all_m": cigar_is_all_m(h),
            "align_start": h.alignment_start_hap_wrt_ref,
            "genome_loc": h.genome_loc.map(|g| format!("{}-{}", g.start_1based(), g.end_1based())),
            "linear_motif": motif,
            "motif_is_tatgtga_class": motif_is_tatg,
            "hamming_vs_coupled": vs_coupled,
            "hamming_vs_apply_ref": vs_apply,
            "hamming_vs_full_ref": vs_full,
            "java_padded_sw_indel": java_sw_apply,
            "java_padded_sw_softclip": java_sw_soft,
            "eventmap_at_327": replay,
            "per_hap_replay_at_327": events_at(&per_hap_replay, TARGET),
            "produces_A_ATG": replay.iter().any(|s| s == "A/ATG"),
            "produces_A_G": replay.iter().any(|s| s == "A/G"),
        }));
    }

    let unique: std::collections::BTreeSet<String> =
        haps.iter().map(|h| fnv1a64_hex(&h.bases)).collect();
    let dump = json!({
        "mission": "6R.116",
        "region": format!("{}:{}-{}", covering.contig, covering.start.get(), covering.end.get()),
        "full_pad": full_pad,
        "apply_pad": apply_pad,
        "full_ref_len": full_ref.len(),
        "apply_ref_len": apply_bases.len(),
        "hap_count": haps.len(),
        "unique_seq_count": unique.len(),
        "n_all_m": haps.iter().filter(|h| cigar_is_all_m(h)).count(),
        "n_indel_cigar": haps.iter().filter(|h| cigar_has_indel(h)).count(),
        "n_supplement_score": haps.iter().filter(|h| (h.score - 1000.0).abs() < 1e-6).count(),
        "per_haplotype": rows,
    });
    println!("{}", serde_json::to_string_pretty(&dump).expect("json"));
}
