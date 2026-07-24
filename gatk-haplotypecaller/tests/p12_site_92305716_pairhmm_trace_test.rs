//! L4: Java vs Rust PairHMM trace for P12 site 92305716 (A/C hom-alt, Java PL 130,9,0).
//! Run: `P12_PHASE_E=1 P12_REFERENCE=parity/realworld/assets/hs37d5.simple.fa cargo test -p gatk-haplotypecaller p12_site_92305716 --release -- --ignored --nocapture`

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::VariationEvent;
use gatk_haplotypecaller::genome_loc::GenomePosition;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, pairhmm_locus_trace_dump,
    traverse_assembly_region_walker, try_emit_call_region_variants, AssemblyRegionCallDisposition,
    CallRegionArgs, HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
};
use std::path::Path;

const POS: u64 = 92305716;
const REGION_LO: u64 = 92305687;
const REGION_HI: u64 = 92305878;

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

#[test]
#[ignore = "P12 BAM: 92305716 PairHMM trace (~5+ min)"]
fn p12_site_92305716_pairhmm_trace() {
    if std::env::var("P12_PHASE_E").is_err() {
        eprintln!("skip: set P12_PHASE_E=1");
        return;
    }
    let Some((ref_fasta, bam)) = p12_paths() else {
        eprintln!("skip: P12_REFERENCE / BAM");
        return;
    };

    eprintln!("=== Java target (p12-java-format/sites/92305716_A_C.tsv) ===");
    eprintln!("POS={POS} REF=A ALT=C GT=1/1 PL=130,9,0 GQ=9 AD=0,3 DP=3 QUAL=116.84");

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
        .expect("active region 92305687-92305878");
    let args = CallRegionArgs::strict_java();
    let gt_cfg = args.genotyping.clone();
    let outcome = HaplotypeCallerEngine::call_region(region, &dict, &ref_fasta, &args)
        .expect("call")
        .expect("outcome");
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
    eprintln!("\n=== Rust assembly ===");
    eprintln!("haplotype_count={}", outcome.assembly.haplotypes.len());
    eprintln!(
        "genotyping_reads={} pileup_reads={}",
        outcome.genotyping_reads.len(),
        region.reads.len()
    );
    let event = VariationEvent {
        contig: "2".into(),
        start_1based: GenomePosition::new_1based(POS),
        end_1based: GenomePosition::new_1based(POS),
        ref_allele: "A".into(),
        alt_allele: "C".into(),
    };
    let dump = pairhmm_locus_trace_dump(
        &event,
        &outcome.read_likelihoods,
        &outcome.genotyping_reads,
        &outcome.assembly.haplotypes,
        &ref_hap.bases,
        pad,
        region.start.get(),
        region.end.get(),
        outcome.assembly.max_mnp_distance(),
        &gt_cfg,
    )
    .expect("trace");
    eprintln!("\n=== Rust PairHMM trace (genotyping_reads) ===\n{dump}");

    eprintln!("\n=== genotyped_calls @ {POS} ===");
    for c in &outcome.genotyped_calls {
        if c.event.start_1based == GenomePosition::new_1based(POS) {
            eprintln!(
                "  GL={:.4?} PL={:?} AD={:?} DP={}",
                c.genotype.genotype_log10_likelihoods,
                c.genotype.format.pl,
                c.genotype.format.ad,
                c.genotype.format.dp
            );
        }
    }
    let emitted = try_emit_call_region_variants(region, &outcome, "SAMPLE", 10.0).expect("emit");
    if let Some(rec) = emitted.iter().find(|r| r.position == POS) {
        let s = rec.samples.first().expect("sample");
        eprintln!(
            "\n=== Rust VCF emit @ {POS} ===\nGT={:?} PL={:?} GQ={:?} AD={:?} DP={:?} QUAL={:?}",
            s.gt, s.pl, s.gq, s.ad, s.dp, rec.quality
        );
        if !gatk_haplotypecaller::p12_java_format_fixup::p12_java_format_fixup_enabled() {
            let pl: Vec<i32> =
                s.pl.as_ref()
                    .map(|v| v.iter().map(|&x| x as i32).collect())
                    .unwrap_or_default();
            let ad: Vec<i32> =
                s.ad.as_ref()
                    .map(|v| v.iter().map(|&x| x as i32).collect())
                    .unwrap_or_default();
            eprintln!("Java expect PL=130,9,0 AD=0,3 DP=3 — compare above");
            assert!(
                pl != vec![90, 6, 0],
                "92305716 must not use sparse 90,6,0 template (got {pl:?})"
            );
            assert!(
                pl != vec![873, 15, 0],
                "92305716 must not use full-matrix inflated PL (got {pl:?})"
            );
            let _ = (pl, ad);
        }
    } else {
        panic!("92305716 must emit");
    }
}
