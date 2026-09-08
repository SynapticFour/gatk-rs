//! 6R.117 live dump: equal-length CIGAR SW fix at `2:92307327`.
//! Skipped unless `HOLDOUT_6R117=1`. Does not reopen 6R.113–6R.116.
//!
//! ```text
//! HOLDOUT_6R117=1 P12_REFERENCE=$PWD/parity/realworld/assets/hs37d5.simple.fa \
//!   cargo test -p gatk-haplotypecaller --test holdout_6r117_equal_length_cigar -- --nocapture
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::{
    build_per_haplotype_variation_events, collect_variation_events, VariationEvent,
};
use gatk_haplotypecaller::hc_allele_mapping::{
    create_allele_mapper_with_events, hap_base_at_ref_locus,
};
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

fn events_at(events: &[VariationEvent], loc: u64) -> Vec<String> {
    events
        .iter()
        .filter(|e| e.start_1based.get() == loc)
        .map(|e| format!("{}/{}", e.ref_allele, e.alt_allele))
        .collect()
}

#[test]
fn holdout_6r117_live_equal_length_cigar_dump() {
    if std::env::var("HOLDOUT_6R117").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R117=1");
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
    let union = collect_variation_events(haps, full_ref, full_pad, "2", max_mnp);
    let union_at = events_at(&union, TARGET);

    let mut rows = Vec::new();
    for (i, h) in haps.iter().enumerate() {
        let motif = linear_motif(h, full_pad);
        let replay = events_at(hap_events.events_for(i), TARGET);
        rows.push(json!({
            "idx": i,
            "hash": fnv1a64_hex(&h.bases),
            "len": h.bases.len(),
            "is_reference": h.is_reference,
            "score": h.score,
            "cigar": cigar_str(h),
            "cigar_is_all_m": cigar_is_all_m(h),
            "linear_motif": motif,
            "eventmap_at_327": replay,
            "produces_A_ATG": replay.iter().any(|s| s == "A/ATG"),
            "produces_A_G": replay.iter().any(|s| s == "A/G"),
        }));
    }

    let map_atg = create_allele_mapper_with_events(
        &VariationEvent::from_alleles("2", TARGET, "A", "ATG"),
        TARGET,
        haps,
        apply_pad,
        apply_bases.as_ref(),
        max_mnp,
        true,
        Some(&hap_events),
    );
    let map_g = create_allele_mapper_with_events(
        &VariationEvent::from_alleles("2", TARGET, "A", "G"),
        TARGET,
        haps,
        apply_pad,
        apply_bases.as_ref(),
        max_mnp,
        true,
        Some(&hap_events),
    );

    let dump = json!({
        "mission": "6R.117",
        "region": format!("{}:{}-{}", covering.contig, covering.start.get(), covering.end.get()),
        "hap_count": haps.len(),
        "eventmap_union_at_327": union_at,
        "union_has_A_ATG": union_at.iter().any(|s| s == "A/ATG"),
        "union_has_A_G": union_at.iter().any(|s| s == "A/G"),
        "mapper_ATG_alt_idx": map_atg.alt_haplotype_indices.iter().map(|i| i.get()).collect::<Vec<_>>(),
        "mapper_G_alt_idx": map_g.alt_haplotype_indices.iter().map(|i| i.get()).collect::<Vec<_>>(),
        "per_haplotype": rows,
    });
    println!("{}", serde_json::to_string_pretty(&dump).expect("json"));
}
