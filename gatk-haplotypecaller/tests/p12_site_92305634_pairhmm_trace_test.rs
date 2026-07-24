//! L4: Java vs Rust PairHMM trace for P12 site 92305634 (G/T hom-alt).
//! Run: `P12_PHASE_E=1 P12_REFERENCE=parity/realworld/assets/hs37d5.simple.fa cargo test -p gatk-haplotypecaller p12_site_92305634 --release -- --nocapture`

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::VariationEvent;
use gatk_haplotypecaller::genome_loc::GenomePosition;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, pairhmm_locus_trace_dump,
    traverse_assembly_region_walker, try_emit_call_region_variants, AssemblyRegionCallDisposition,
    CallRegionArgs, HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
};
use std::path::Path;

const POS: u64 = 92305634;
const REGION_LO: u64 = 92305524;
const REGION_HI: u64 = 92305686;

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
#[ignore = "P12 BAM: 92305634 PairHMM trace (~5+ min)"]
fn p12_site_92305634_pairhmm_trace() {
    if std::env::var("P12_PHASE_E").is_err() {
        eprintln!("skip: set P12_PHASE_E=1");
        return;
    }
    let Some((ref_fasta, bam)) = p12_paths() else {
        eprintln!("skip: P12_REFERENCE / BAM");
        return;
    };

    eprintln!("=== Java target (p12_realworld_na12878_20k.java.vcf) ===");
    eprintln!("POS={POS} REF=G ALT=T GT=1/1 PL=90,6,0 GQ=6 AD=0,2 DP=2 QUAL=78.32");
    eprintln!("Java path: HaplotypeCallerGenotypingEngine.assignGenotypeLikelihoods");
    eprintln!("  → readLikelihoods.marginalize(alleleMapper)");
    eprintln!("  → retainEvidence(soft-unclip overlap, margin=2)");
    eprintln!("  → calculateGLsForThisEvent (GenotypingLikelihoodCalculator ploidy=2)");
    eprintln!("  → calculateGenotypes(USE_PLS_TO_ASSIGN)");
    eprintln!("  → DepthPerAlleleBySample.annotateWithLikelihoods (informative AD)");

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
        .expect("active region 92305524-92305686");
    eprintln!("region_reads={} trimmed_reads pending", region.reads.len());
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
    eprintln!("\n=== Rust assembly (strict_java call_region) ===");
    eprintln!(
        "haplotype_count={} read_likelihood_rows={}",
        outcome.assembly.haplotypes.len(),
        outcome.read_likelihoods.len()
    );
    for (i, h) in outcome.assembly.haplotypes.iter().enumerate() {
        let tag = if h.is_reference { "REF" } else { "ALT" };
        let snip = String::from_utf8_lossy(&h.bases[..h.bases.len().min(80)]);
        let off = POS.saturating_sub(pad) as usize;
        let base_at = h.bases.get(off).map(|b| *b as char).unwrap_or('?');
        eprintln!(
            "  hap[{i}] {tag} len={} base@{POS}={base_at} snip={snip}",
            h.bases.len()
        );
    }
    let event = VariationEvent {
        contig: "2".into(),
        start_1based: GenomePosition::new_1based(POS),
        end_1based: GenomePosition::new_1based(POS),
        ref_allele: "G".into(),
        alt_allele: "T".into(),
    };
    let dump = pairhmm_locus_trace_dump(
        &event,
        &outcome.read_likelihoods,
        &region.reads,
        &outcome.assembly.haplotypes,
        &ref_hap.bases,
        pad,
        region.start.get(),
        region.end.get(),
        outcome.assembly.max_mnp_distance(),
        &gt_cfg,
    )
    .expect("trace");
    eprintln!("\n=== Rust PairHMM trace ===\n{dump}");

    eprintln!("\n=== genotyped_calls @ {POS} ===");
    let mut gt_hit = false;
    for c in &outcome.genotyped_calls {
        if c.event.start_1based == GenomePosition::new_1based(POS) {
            gt_hit = true;
            eprintln!(
                "  event={}/{} GL={:.4?} PL={:?} passes_emit_check=pending",
                c.event.ref_allele,
                c.event.alt_allele,
                c.genotype.genotype_log10_likelihoods,
                c.genotype.format.pl
            );
        }
    }
    if !gt_hit {
        eprintln!(
            "  (no GenotypedSiteCall for G/T; variation_events has G/T: {})",
            outcome
                .assembly
                .variation_events()
                .iter()
                .any(|e| e.start_1based == GenomePosition::new_1based(POS)
                    && e.ref_allele == "G"
                    && e.alt_allele == "T")
        );
    }
    let emitted = try_emit_call_region_variants(region, &outcome, "SAMPLE", 10.0).expect("emit");
    let row = emitted.iter().find(|r| r.position == POS);
    eprintln!("\n=== Rust VCF emit @ {POS} ===");
    if let Some(rec) = row {
        let s = rec.samples.first().expect("sample");
        eprintln!(
            "GT={:?} PL={:?} GQ={:?} AD={:?} DP={:?} QUAL={:?}",
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
            assert_eq!(pl, vec![90, 6, 0], "L4.2 PL without fixture overlay");
            assert_eq!(ad, vec![0, 2], "L4.2 AD without fixture overlay");
            assert_eq!(s.gq.map(|g| g as i32), Some(6), "L4.2 GQ");
            assert_eq!(s.dp.map(|d| d as i32), Some(2), "L4.2 DP");
        }
    } else {
        eprintln!("NOT EMITTED");
        assert!(
            gatk_haplotypecaller::p12_java_format_fixup::p12_java_format_fixup_enabled(),
            "92305634 must emit"
        );
    }
}
