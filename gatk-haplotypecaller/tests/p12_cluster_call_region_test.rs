//! P12 cluster active region: `call_region` variation events + genotyped calls.
//! Run: `P12_REFERENCE=…/hs37d5.simple.fa cargo test -p gatk-haplotypecaller p12_cluster_call_region --release -- --nocapture`

use gatk_core::reference::ReferenceWindowCache;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::assembly_region_finalize::{
    assembly_reference_read, finalize_region_reads_for_assembly,
    gatk_min_tail_quality_for_assembly, records_to_assembly_reads,
};
use gatk_haplotypecaller::cigar::{Cigar, CigarOperator};
use gatk_haplotypecaller::genome_loc::GenomePosition;
use gatk_haplotypecaller::haplotype::Haplotype;
use gatk_haplotypecaller::read_threading_assembler::{
    audit_threading_dangling_recovery, build_threading_graph_for_haplotype_dump,
};
use gatk_haplotypecaller::{
    assemble_reads, audit_kbest_extract, call_disposition, diagnose_genotype_variation_event,
    find_best_haplotypes_for_assembly, flatten_assembly_regions, format_locus_genotype_pl_dump,
    traverse_assembly_region_walker, AssembleReadsArgs, AssemblyRegionCallDisposition,
    CallRegionArgs, HaplotypeCallerEngine, HcGenotypingConfig, ReadFilterParams,
    ReadThreadingAssemblerArgs, WalkerTraversalConfig,
};
use std::path::Path;

fn p12_cluster_debug_enabled() -> bool {
    std::env::var("GATK_RS_P12_CLUSTER_DEBUG")
        .ok()
        .is_some_and(|v| matches!(v.as_str(), "1" | "true" | "TRUE"))
}

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
    if !ref_path.is_file() || !bam.is_file() {
        return None;
    }
    Some((ref_path, bam))
}

#[test]
fn p12_cluster_call_region() {
    let Some((ref_fasta, bam)) = p12_paths() else {
        eprintln!("skip: set P12_REFERENCE");
        return;
    };
    let interval = "2:92307228-92307400";
    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, interval).expect("interval");
    let filters = ReadFilterParams::gatk_standard_hc();
    let cfg = WalkerTraversalConfig::gatk_haplotype_caller_production(100);
    let walk = traverse_assembly_region_walker(&dict, &specs, &ref_fasta, &bam, &filters, &cfg)
        .expect("walk");
    let regions = flatten_assembly_regions(&walk);
    let region = regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= 92307324
                && r.end.get() >= 92307327
        })
        .expect("cluster active region");

    // Optional slow path: `GATK_RS_CIGAR_EX=1` runs extra assemble + KBest audit (~30s).
    if std::env::var("GATK_RS_CIGAR_EX")
        .ok()
        .is_some_and(|v| matches!(v.as_str(), "1" | "true" | "TRUE"))
    {
        let mut ref_cache = ReferenceWindowCache::new(ref_fasta.clone(), 4);
        let untrimmed =
            assemble_reads(region, &dict, &mut ref_cache, &AssembleReadsArgs::default())
                .expect("assemble");
        eprintln!(
            "full_asm\thaps={} events={}",
            untrimmed.haplotypes.len(),
            untrimmed.variation_events().len()
        );
        let reference = assembly_reference_read(&dict, &mut ref_cache, region).expect("ref");
        let finalized = finalize_region_reads_for_assembly(
            &region.reads,
            region,
            true,
            gatk_min_tail_quality_for_assembly(10),
            false,
        );
        let reads = records_to_assembly_reads(&finalized);
        let mut asm_args = ReadThreadingAssemblerArgs::default();
        asm_args.kmer_sizes = vec![85];
        if let Some(graph) =
            build_threading_graph_for_haplotype_dump(&reference, &reads, 85, &asm_args, true, true)
                .expect("graph")
        {
            let mut ref_hap = Haplotype::new(reference.bases.as_slice(), true);
            let mut ref_cigar = Cigar::new();
            ref_cigar.push(ref_hap.bases.len(), CigarOperator::Match);
            ref_hap.cigar = Some(ref_cigar);
            let ref_cigar_len = ref_hap.cigar.as_ref().unwrap().reference_length();
            let (paths, graph) = find_best_haplotypes_for_assembly(graph, 128).expect("kbest");
            eprintln!("cigar_ex\tkbest_paths={}", paths.len());
            for row in audit_kbest_extract(
                &paths,
                &graph,
                &ref_hap,
                ref_cigar_len,
                &asm_args.haplotype_to_reference_sw,
            ) {
                eprintln!(
                    "cigar_ex_path\tidx={} eq_ref={} outcome={:?}",
                    row.path_index, row.eq_ref_bases, row.outcome
                );
            }
        }
    }

    // ASM-1 L2: dangling recovery on cluster graph (optional diagnostics).
    if p12_cluster_debug_enabled() {
        let mut ref_cache = ReferenceWindowCache::new(ref_fasta.clone(), 4);
        let reference = assembly_reference_read(&dict, &mut ref_cache, region).expect("ref");
        let finalized = finalize_region_reads_for_assembly(
            &region.reads,
            region,
            true,
            gatk_min_tail_quality_for_assembly(10),
            false,
        );
        let reads = records_to_assembly_reads(&finalized);
        let asm_args = ReadThreadingAssemblerArgs::default();
        for &kmer in &[85usize, 25, 10] {
            if let Some(audit) =
                audit_threading_dangling_recovery(&reference, &reads, kmer, &asm_args, true, true)
                    .expect("dangling audit")
            {
                eprintln!(
                    "cluster_dangling\tk={kmer}\ttails {}/{} heads {}/{} edges {} -> {}",
                    audit.tails_recovered,
                    audit.tails_attempted,
                    audit.heads_recovered,
                    audit.heads_attempted,
                    audit.edges_before,
                    audit.edges_after
                );
            } else {
                eprintln!("cluster_dangling\tk={kmer}\tno_graph");
            }
        }
    }

    let strict_args = CallRegionArgs::strict_java();
    let outcome = HaplotypeCallerEngine::call_region(region, &dict, &ref_fasta, &strict_args)
        .expect("call")
        .expect("some outcome");

    if p12_cluster_debug_enabled() {
        eprintln!(
            "variation_present={} events={} genotyped_calls={}",
            outcome.assembly.is_variation_present(),
            outcome.assembly.variation_events().len(),
            outcome.genotyped_calls.len()
        );
        for c in &outcome.genotyped_calls {
            eprintln!(
                "genotyped\t{} {}/{} gl={:?} pl={:?}",
                c.event.start_1based.get(),
                c.event.ref_allele,
                c.event.alt_allele,
                c.genotype.genotype_log10_likelihoods,
                c.genotype.format.pl
            );
        }
        for (i, h) in outcome.assembly.haplotypes.iter().enumerate() {
            eprintln!(
                "hap{i}\tref={}\tlen={}\talign={}\tpad={}\tcigar={}",
                h.is_reference,
                h.bases.len(),
                h.alignment_start_hap_wrt_ref,
                outcome.assembly.padded_reference_start_1based(),
                h.cigar
                    .as_ref()
                    .map(|c| c.to_gatk_string())
                    .unwrap_or_default()
            );
            if std::env::var("GATK_RS_ASM8_DEBUG")
                .ok()
                .is_some_and(|v| matches!(v.as_str(), "1" | "true" | "TRUE"))
                && !h.is_reference
                && h.cigar.as_ref().is_some_and(|c| {
                    c.to_gatk_string().contains('D') || c.to_gatk_string().contains('I')
                })
            {
                use gatk_haplotypecaller::event_map::EventMap;
                let full_ref = outcome.assembly.reference_bases();
                let full_pad = outcome.assembly.padded_reference_start_1based();
                let ref_hap = outcome
                    .assembly
                    .haplotypes
                    .iter()
                    .find(|x| x.is_reference)
                    .expect("ref");
                let evs = EventMap::from_haplotype_and_reference(
                    h,
                    ref_hap,
                    full_ref,
                    full_pad,
                    outcome.assembly.max_mnp_distance(),
                )
                .variation_events("2", full_pad);
                for e in evs.iter().filter(|e| {
                    e.start_1based.get() >= 92307320 && e.start_1based.get() <= 92307335
                }) {
                    eprintln!(
                        "asm8_eventmap\thap={i}\t{} {}/{}",
                        e.start_1based.get(),
                        e.ref_allele,
                        e.alt_allele
                    );
                }
            }
        }
        let active_lo = region.start.get();
        let active_hi = region.end.get();
        let in_active = outcome
            .assembly
            .variation_events()
            .iter()
            .filter(|e| e.start_1based.get() >= active_lo && e.start_1based.get() <= active_hi)
            .count();
        eprintln!("events_in_active_span\t{in_active}");
        if outcome.assembly.haplotypes.len() >= 2 {
            let alt = &outcome.assembly.haplotypes[0];
            let rf = &outcome.assembly.haplotypes[1];
            let mut diffs = 0usize;
            for (i, (a, b)) in alt.bases.iter().zip(rf.bases.iter()).enumerate() {
                if a != b {
                    diffs += 1;
                    if diffs <= 5 {
                        let g = outcome.assembly.padded_reference_start_1based()
                            + alt.alignment_start_hap_wrt_ref as u64
                            + i as u64;
                        eprintln!("base_diff\tgenomic={g} i={i} alt={a} ref={b}");
                    }
                }
            }
            eprintln!("base_diff_count\t{diffs} len={}", alt.bases.len());
        }
        for e in outcome.assembly.variation_events() {
            eprintln!(
                "event\t{} {}/{}",
                e.start_1based.get(),
                e.ref_allele,
                e.alt_allele
            );
        }
        for c in &outcome.genotyped_calls {
            eprintln!(
                "call\t{} {}/{}",
                c.event.start_1based.get(),
                c.event.ref_allele,
                c.event.alt_allele
            );
        }
    }

    if std::env::var("P12_CLUSTER_PL")
        .ok()
        .is_some_and(|v| matches!(v.as_str(), "1" | "true" | "TRUE"))
    {
        let ref_hap = outcome
            .assembly
            .haplotypes
            .iter()
            .find(|h| h.is_reference)
            .expect("ref");
        let pad = ref_hap
            .genome_loc
            .as_ref()
            .map(|g| g.start_1based())
            .unwrap_or_else(|| outcome.assembly.padded_reference_start_1based());
        let gt_cfg = HcGenotypingConfig::strict_java();
        for (pos, ra, aa) in [
            ("92307324", "TTC", "T"),
            ("92307327", "A", "ATG"),
            ("92307333", "T", "G"),
        ] {
            let pos: u64 = pos.parse().unwrap();
            let Some(ev) = outcome.assembly.variation_events().iter().find(|e| {
                e.start_1based == GenomePosition::new_1based(pos)
                    && e.ref_allele == ra
                    && e.alt_allele == aa
            }) else {
                eprintln!("pl_dump\t{pos}\tmissing_event");
                continue;
            };
            let dump = format_locus_genotype_pl_dump(
                ev,
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
            .expect("pl dump");
            eprintln!("pl_dump\t{pos}\n{dump}");
            if pos == 92307324 {
                let ttc_off = pos.saturating_sub(pad) as usize;
                eprintln!(
                    "pl_dump\tref_slice_at_ttc={:?}",
                    ref_hap.bases.get(ttc_off..ttc_off.saturating_add(3))
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
            ) {
                Ok(Ok(_)) => eprintln!("pl_dump\t{pos}\tdiagnose=ok"),
                Ok(Err(r)) => eprintln!("pl_dump\t{pos}\tdiagnose={r:?}"),
                Err(e) => eprintln!("pl_dump\t{pos}\terr={e}"),
            }
        }
    }

    assert!(
        outcome.assembly.haplotypes.len() > 1,
        "strict_java: expected alt haplotypes"
    );

    let strict_phase_a = std::env::var("GATK_RS_PHASE_A_STRICT")
        .ok()
        .is_some_and(|v| matches!(v.as_str(), "1" | "true" | "TRUE"));
    let strict_genotype = std::env::var("GATK_RS_PHASE_A_STRICT_GENOTYPE")
        .ok()
        .is_some_and(|v| matches!(v.as_str(), "1" | "true" | "TRUE"));
    if strict_phase_a {
        let indel_alt_haps = outcome
            .assembly
            .haplotypes
            .iter()
            .filter(|h| {
                !h.is_reference
                    && h.cigar
                        .as_ref()
                        .is_some_and(|c| c.elements.iter().any(|e| e.operator.is_indel()))
            })
            .count();
        if p12_cluster_debug_enabled() {
            eprintln!("strict_phase_b\tindel_alt_haps={indel_alt_haps}");
        }
        assert!(
            indel_alt_haps >= 1 || outcome.assembly.haplotypes.len() >= 3,
            "Phase B.2: expect alt hap(s) with I/D CIGAR or multiple alts; haps={}",
            outcome.assembly.haplotypes.len()
        );
        assert!(
            outcome.assembly.variation_events().iter().any(|e| {
                e.start_1based == GenomePosition::new_1based(92307324)
                    && e.ref_allele == "TTC"
                    && e.alt_allele == "T"
            }),
            "Phase B.3: TTC/T from EventMap (no inject); events={:?}",
            outcome.assembly.variation_events()
        );
        if strict_genotype {
            let emitted: std::collections::BTreeSet<_> = outcome
                .genotyped_calls
                .iter()
                .map(|c| {
                    (
                        c.event.start_1based.get(),
                        c.event.ref_allele.clone(),
                        c.event.alt_allele.clone(),
                    )
                })
                .collect();
            assert!(
                emitted.contains(&(92307324, "TTC".into(), "T".into())),
                "TTC/T genotyped"
            );
            assert!(
                emitted.contains(&(92307327, "A".into(), "ATG".into())),
                "A/ATG genotyped"
            );
            assert!(
                emitted.contains(&(92307333, "T".into(), "G".into())),
                "T/G genotyped"
            );
        } else if p12_cluster_debug_enabled() {
            eprintln!(
                "parity_phase_a: events ok; genotyped_calls={} (Phase E: GATK_RS_PHASE_A_STRICT_GENOTYPE=1)",
                outcome.genotyped_calls.len()
            );
        }
    } else if p12_cluster_debug_enabled() {
        eprintln!(
            "parity_phase_a: events={} genotyped={} (set GATK_RS_PHASE_A_STRICT=1 when ASM-1 lands)",
            outcome.assembly.variation_events().len(),
            outcome.genotyped_calls.len()
        );
    }
}

#[cfg(feature = "parity_harness")]
#[test]
fn p12_cluster_call_region_legacy_bridges() {
    let Some((ref_fasta, bam)) = p12_paths() else {
        eprintln!("skip: set P12_REFERENCE");
        return;
    };
    let interval = "2:92307228-92307422";
    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, interval).expect("interval");
    let filters = ReadFilterParams::gatk_standard_hc();
    let cfg = WalkerTraversalConfig::gatk_haplotype_caller_production(100);
    let walk = traverse_assembly_region_walker(&dict, &specs, &ref_fasta, &bam, &filters, &cfg)
        .expect("walk");
    let regions = flatten_assembly_regions(&walk);
    let region = regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= 92307324
                && r.end.get() >= 92307327
        })
        .expect("cluster active region");

    let outcome = HaplotypeCallerEngine::call_region(
        region,
        &dict,
        &ref_fasta,
        &CallRegionArgs::legacy_read_bridges(),
    )
    .expect("call")
    .expect("some outcome");

    assert!(
        outcome.genotyped_calls.iter().any(|c| c.event.start_1based
            == GenomePosition::new_1based(92307324)
            && c.event.ref_allele == "TTC"),
        "legacy bridges TTC/T"
    );
}
