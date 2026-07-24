//! Cluster emit probe: TTC/T, T/G, CT/C on cluster vs full-walker regions.
//! Run: `P12_PHASE_E=1 P12_REFERENCE=…/hs37d5.simple.fa cargo test -p gatk-haplotypecaller p12_site_92307324_emit_probe --release -- --nocapture`

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::genome_loc::GenomePosition;
use gatk_haplotypecaller::{
    call_disposition, diagnose_genotype_variation_event, flatten_assembly_regions,
    hc_allele_mapping::create_allele_mapper, hc_emit_policy::explain_strict_java_emit_gates,
    read_event_discovery::read_allele_depths_at_locus, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
};

fn p12_paths() -> Option<(std::path::PathBuf, std::path::PathBuf)> {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
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

fn probe_site(pos: u64, ref_a: &str, alt_a: &str, label: &str, interval: &str) {
    let Some((ref_fasta, bam)) = p12_paths() else {
        eprintln!("skip: P12_REFERENCE / BAM");
        return;
    };
    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, interval).expect("interval");
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
            ) && r.start.get() <= pos
                && r.end.get() >= pos
        })
        .expect("active region");
    let args = CallRegionArgs::strict_java();
    let gt_cfg = args.genotyping.clone();

    eprintln!(
        "=== {label} {pos} {ref_a}/{alt_a} interval={interval} region={}-{} reads={} ===",
        region.start.get(),
        region.end.get(),
        region.reads.len()
    );

    let outcome = HaplotypeCallerEngine::call_region(region, &dict, &ref_fasta, &args)
        .expect("call")
        .expect("outcome");

    let event_hit = outcome.assembly.variation_events().iter().find(|e| {
        e.start_1based == GenomePosition::new_1based(pos)
            && e.ref_allele == ref_a
            && e.alt_allele == alt_a
    });
    eprintln!(
        "{label}\tevent_map\thit={} events={}",
        event_hit.is_some(),
        outcome.assembly.variation_events().len()
    );

    let genotyped = outcome.genotyped_calls.iter().find(|c| {
        c.event.start_1based == GenomePosition::new_1based(pos)
            && c.event.ref_allele == ref_a
            && c.event.alt_allele == alt_a
    });
    eprintln!(
        "{label}\tgenotyped\thit={} total_calls={}",
        genotyped.is_some(),
        outcome.genotyped_calls.len()
    );
    if let Some(c) = genotyped {
        eprintln!(
            "{label}\tcall\tPL={:?} AD={:?}",
            c.genotype.format.pl, c.genotype.format.ad
        );
    }

    let emitted = try_emit_call_region_variants(region, &outcome, "SAMPLE", 10.0).expect("emit");
    let emit_hit = emitted.iter().any(|r| {
        r.position == pos
            && r.reference == ref_a
            && r.alternate.first().map(String::as_str) == Some(alt_a)
    });
    eprintln!("{label}\temit\thit={} rows={}", emit_hit, emitted.len());

    if let Some(ev) = event_hit {
        let ref_hap = outcome
            .assembly
            .haplotypes
            .iter()
            .find(|h| h.is_reference)
            .expect("ref hap");
        let pad = ref_hap
            .genome_loc
            .as_ref()
            .map(|g| g.start_1based())
            .unwrap_or_else(|| outcome.assembly.padded_reference_start_1based());
        let mapping = create_allele_mapper(
            ev,
            ev.start_1based.get(),
            &outcome.assembly.haplotypes,
            pad,
            &ref_hap.bases,
            outcome.assembly.max_mnp_distance(),
            true,
        );
        eprintln!(
            "{label}\tallele_map\tref_haps={:?} alt_haps={:?} haps={}",
            mapping.ref_haplotype_indices,
            mapping.alt_haplotype_indices,
            outcome.assembly.haplotypes.len()
        );
        let off = pos.saturating_sub(pad) as usize;
        for (i, h) in outcome.assembly.haplotypes.iter().enumerate() {
            eprintln!(
                "{label}\thap{i}\tlen={} base@off={:?} cigar={}",
                h.bases.len(),
                h.bases.get(off),
                h.cigar
                    .as_ref()
                    .map(|c| c.to_gatk_string())
                    .unwrap_or_default()
            );
        }
        match diagnose_genotype_variation_event(
            ev,
            &outcome.read_likelihoods,
            &outcome.genotyping_reads,
            &outcome.genotyping_reads,
            Some(region.reads.as_slice()),
            &outcome.assembly.haplotypes,
            &ref_hap.bases,
            pad,
            outcome.assembly.reference_bases(),
            outcome.assembly.padded_reference_start_1based(),
            region.start.get(),
            region.end.get(),
            outcome.assembly.max_mnp_distance(),
            &gt_cfg,
        )
        .expect("diagnose")
        {
            Ok(call) => {
                let (rr, ra) = read_allele_depths_at_locus(&region.reads, ev, pad);
                let gates = explain_strict_java_emit_gates(
                    &call.event,
                    &call.genotype.genotype_log10_likelihoods,
                    &call.genotype.format,
                    gt_cfg.stand_emit_confidence,
                    gt_cfg.genotype_stored_events_only,
                    rr,
                    ra,
                    &[],
                )
                .expect("gates");
                eprintln!(
                    "{label}\tdiagnose\tok PL={:?} {gates}",
                    call.genotype.format.pl
                );
            }
            Err(reason) => eprintln!("{label}\tdiagnose\treject={reason:?}"),
        }
    }
}

#[test]
fn p12_site_92307324_emit_probe() {
    if std::env::var("P12_PHASE_E").is_err() {
        eprintln!("skip: set P12_PHASE_E=1");
        return;
    }
    probe_site(
        92307364,
        "T",
        "C",
        "cluster_interval",
        "2:92307228-92307400",
    );
    probe_site(92307364, "T", "C", "full_walker", "2:92300000-92350000");
    probe_site(
        92307324,
        "TTC",
        "T",
        "cluster_interval",
        "2:92307228-92307400",
    );
    probe_site(92307324, "TTC", "T", "full_walker", "2:92300000-92350000");
    probe_site(92307333, "T", "G", "full_walker", "2:92300000-92350000");
    probe_site(92307359, "CT", "C", "full_walker", "2:92300000-92350000");
}
