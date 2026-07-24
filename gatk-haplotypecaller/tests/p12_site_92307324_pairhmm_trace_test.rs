//! L4: Java vs Rust marginalize → retainEvidence trace for P12 cluster 92307324/92307327.
//! Run: `P12_PHASE_E=1 P12_REFERENCE=…/hs37d5.simple.fa cargo test -p gatk-haplotypecaller p12_site_92307324_pairhmm --release -- --ignored --nocapture`

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::genome_loc::GenomePosition;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, pairhmm_locus_trace_dump,
    read_event_discovery::p12_cluster_coupled_indel_supporting_read_qnames,
    traverse_assembly_region_walker, try_emit_call_region_variants, AssemblyRegionCallDisposition,
    CallRegionArgs, HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
};
use std::path::Path;

const TTC_POS: u64 = 92307324;
const ATG_POS: u64 = 92307327;
const REGION_LO: u64 = 92307229;
const REGION_HI: u64 = 92307386;

fn p12_paths() -> Option<(std::path::PathBuf, std::path::PathBuf)> {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
    let ref_path = std::env::var("P12_REFERENCE")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|_| root.join("parity/realworld/assets/hs37d5.simple.fa"));
    let ref_path = if ref_path.is_absolute() {
        ref_path
    } else {
        root.join(ref_path)
    };
    let bam = root.join("parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam");
    if ref_path.is_file() && bam.is_file() {
        Some((ref_path, bam))
    } else {
        None
    }
}

fn trace_site(
    label: &str,
    event: gatk_haplotypecaller::event_map::VariationEvent,
    outcome: &gatk_haplotypecaller::CallRegionOutcome,
    region: &gatk_haplotypecaller::AssemblyRegion,
    pad: u64,
    ref_bases: &[u8],
    gt_cfg: &gatk_haplotypecaller::HcGenotypingConfig,
) {
    let full_ref = outcome.assembly.reference_bases();
    let full_pad = outcome.assembly.padded_reference_start_1based();
    let support =
        p12_cluster_coupled_indel_supporting_read_qnames(&region.reads, &event, full_ref, full_pad);
    eprintln!(
        "\n=== {label} pileup_indel_support_qnames (full ref pad) count={} ===",
        support.len()
    );
    for q in &support {
        eprintln!("  qname={}", String::from_utf8_lossy(q));
    }

    let dump = pairhmm_locus_trace_dump(
        &event,
        &outcome.read_likelihoods,
        &outcome.genotyping_reads,
        &outcome.assembly.haplotypes,
        ref_bases,
        pad,
        region.start.get(),
        region.end.get(),
        outcome.assembly.max_mnp_distance(),
        gt_cfg,
    )
    .expect("trace");
    eprintln!("\n=== {label} PairHMM trace ===\n{dump}");
}

#[test]
#[ignore = "P12 BAM: 92307324 cluster PairHMM trace (~3+ min)"]
fn p12_site_92307324_pairhmm_trace() {
    if std::env::var("P12_PHASE_E").is_err() {
        eprintln!("skip: set P12_PHASE_E=1");
        return;
    }
    let Some((ref_fasta, bam)) = p12_paths() else {
        eprintln!("skip: P12_REFERENCE / BAM");
        return;
    };

    eprintln!("=== Java target (p12_realworld_na12878_20k.java.vcf) ===");
    eprintln!("92307324 TTC/T GT=1/1 PL=45,3,0 GQ=3 AD=0,1 DP=1 QUAL=35.44");
    eprintln!("92307327 A/ATG GT=1/1 PL=45,3,0 GQ=3 AD=0,1 DP=1 QUAL=35.44");
    eprintln!("Java: marginalize(alleleMapper) → retainEvidence(target.overlaps(read)) → calculateGLsForThisEvent");

    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, "2:92300000-92350000").expect("interval");
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
            ) && r.start.get() <= REGION_LO
                && r.end.get() >= REGION_HI
        })
        .expect("active region");
    eprintln!(
        "region={}-{} pileup_reads={} genotyping_reads pending",
        region.start.get(),
        region.end.get(),
        region.reads.len()
    );

    let args = CallRegionArgs::strict_java();
    let gt_cfg = args.genotyping.clone();
    let outcome = HaplotypeCallerEngine::call_region(region, &dict, &ref_fasta, &args)
        .expect("call")
        .expect("outcome");

    eprintln!(
        "genotyping_reads={} read_ll_rows={} haplotypes={}",
        outcome.genotyping_reads.len(),
        outcome.read_likelihoods.len(),
        outcome.assembly.haplotypes.len()
    );
    for (i, h) in outcome.assembly.haplotypes.iter().enumerate() {
        eprintln!(
            "  hap[{i}] ref={} len={} cigar={}",
            h.is_reference,
            h.bases.len(),
            h.cigar
                .as_ref()
                .map(|c| c.to_gatk_string())
                .unwrap_or_default()
        );
    }

    let ref_hap = outcome
        .assembly
        .haplotypes
        .iter()
        .find(|h| h.is_reference)
        .expect("ref hap");
    let pad = ref_hap
        .genome_loc
        .map(|g| g.start_1based())
        .unwrap_or_else(|| outcome.assembly.padded_reference_start_1based());

    for (label, pos, ref_a, alt_a) in [
        ("TTC/T", TTC_POS, "TTC", "T"),
        ("A/ATG", ATG_POS, "A", "ATG"),
    ] {
        let event = outcome
            .assembly
            .variation_events()
            .iter()
            .find(|e| {
                e.start_1based == GenomePosition::new_1based(pos)
                    && e.ref_allele == ref_a
                    && e.alt_allele == alt_a
            })
            .cloned()
            .unwrap_or_else(|| gatk_haplotypecaller::event_map::VariationEvent {
                contig: "2".into(),
                start_1based: GenomePosition::new_1based(pos),
                end_1based: GenomePosition::new_1based(pos),
                ref_allele: ref_a.into(),
                alt_allele: alt_a.into(),
            });
        trace_site(label, event, &outcome, region, pad, &ref_hap.bases, &gt_cfg);
    }

    let emitted = try_emit_call_region_variants(region, &outcome, "SAMPLE", 10.0).expect("emit");
    eprintln!("\n=== Rust VCF emit ===");
    for pos in [TTC_POS, ATG_POS] {
        if let Some(rec) = emitted.iter().find(|r| r.position == pos) {
            let s = rec.samples.first().expect("sample");
            eprintln!(
                "POS={pos} GT={:?} PL={:?} GQ={:?} AD={:?} DP={:?} QUAL={:?}",
                s.gt, s.pl, s.gq, s.ad, s.dp, rec.quality
            );
        } else {
            eprintln!("POS={pos} NOT EMITTED");
        }
    }
}
