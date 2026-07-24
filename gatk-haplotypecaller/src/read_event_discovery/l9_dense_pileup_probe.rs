use super::*;
use rust_htslib::bam::{self, Read as _};
use std::path::PathBuf;

fn dense_bam() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../parity/realworld/na12878_giab_window_b37/NA12878_giab_window.b37.bam")
}

fn dense_ref() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../parity/realworld/assets/hs37d5.simple.fa")
}

fn load_reads_at(pos: u64) -> Vec<bam::Record> {
    let mut reader = bam::IndexedReader::from_path(dense_bam()).expect("bam");
    let tid = reader.header().tid(b"20").expect("tid") as u32;
    reader
        .fetch((tid, (pos - 1) as i64, pos as i64))
        .expect("fetch");
    let mut out = Vec::new();
    for r in reader.records() {
        out.push(r.expect("rec"));
    }
    out
}

#[test]
fn pileup_matches_samtools_at_gt_flip_sites() {
    if !dense_bam().is_file() {
        eprintln!("skip: dense BAM missing");
        return;
    }
    for (pos, ref_a, alt_a, min_ref, min_alt) in [
        (10009227u64, "A", "G", 12i32, 12i32),
        (10012384, "T", "C", 10, 10),
        (10012636, "G", "C", 8, 8),
    ] {
        let reads = load_reads_at(pos);
        let event = VariationEvent {
            contig: "20".into(),
            start_1based: GenomePosition::new_1based(pos),
            end_1based: GenomePosition::new_1based(pos),
            ref_allele: ref_a.into(),
            alt_allele: alt_a.into(),
        };
        // pad cancels for SNPs; use genomic pad=1.
        let (rr, ra) = read_allele_depths_at_locus(&reads, &event, 1);
        let (dr, da) = read_allele_depths_at_locus_dedupe_qname(&reads, &event, 1);
        eprintln!(
            "L9 pileup {pos} {ref_a}>{alt_a}: n={} raw={rr},{ra} dedupe={dr},{da}",
            reads.len()
        );
        assert!(
            rr >= min_ref && ra >= min_alt,
            "{pos}: raw pileup {rr},{ra} expected roughly balanced (min {min_ref},{min_alt})"
        );
        assert!(
            dr >= min_ref.saturating_sub(2) && da >= min_alt.saturating_sub(2),
            "{pos}: dedupe pileup {dr},{da} unexpectedly alt-skewed"
        );
    }
}

/// Stage dump for residual SNP FNs (discovery vs EventMap vs genotype).
#[test]
fn snp_fn_sites_discover_and_call_region_stage() {
    use crate::engine::{CallRegionArgs, HaplotypeCallerEngine};
    use crate::read_model::ReadFilterParams;
    use crate::walker_traversal::{flatten_assembly_regions, traverse_assembly_region_walker};
    use crate::{call_disposition, AssemblyRegionCallDisposition, WalkerTraversalConfig};
    use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};

    if !dense_bam().is_file() || !dense_ref().is_file() {
        eprintln!("skip: dense assets missing");
        return;
    }
    let ref_fa = dense_ref();
    let bam = dense_bam();
    let dict = SequenceDictionary::from_fasta_path(&ref_fa).expect("dict");
    let filters = ReadFilterParams::gatk_standard_hc();
    let cfg = WalkerTraversalConfig::gatk_haplotype_caller_production(100);
    let args = CallRegionArgs::strict_java();

    for (interval, pos, ref_a, alt_a) in [
        ("20:10001400-10001700", 10001474u64, "C", "T"),
        ("20:10008100-10008300", 10008221, "T", "C"),
        ("20:10036900-10037200", 10037037, "C", "T"),
    ] {
        let specs = parse_intervals_cli_string(&dict, interval).expect("interval");
        let walk = traverse_assembly_region_walker(&dict, &specs, &ref_fa, &bam, &filters, &cfg)
            .expect("walk");
        let regions = flatten_assembly_regions(&walk);
        let Some(region) = regions.iter().find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= pos
                && r.end.get() >= pos
        }) else {
            eprintln!("L9-FN {pos}: NO active region covering site in {interval}");
            continue;
        };

        let event = VariationEvent {
            contig: "20".into(),
            start_1based: GenomePosition::new_1based(pos),
            end_1based: GenomePosition::new_1based(pos),
            ref_allele: ref_a.into(),
            alt_allele: alt_a.into(),
        };
        let (rr, ra) = read_allele_depths_at_locus(&region.reads, &event, 1);

        let outcome = HaplotypeCallerEngine::call_region(region, &dict, &ref_fa, &args)
            .expect("call_region")
            .expect("outcome");
        let target = GenomePosition::new_1based(pos);
        let in_events =
            outcome.assembly.variation_events().iter().any(|e| {
                e.start_1based == target && e.ref_allele == ref_a && e.alt_allele == alt_a
            });
        let in_calls = outcome.genotyped_calls.iter().any(|c| {
            c.event.start_1based == target
                && c.event.ref_allele == ref_a
                && c.event.alt_allele == alt_a
        });
        let on_cigar = variation_event_on_haplotype_cigars(
            &event,
            &outcome.assembly.haplotypes,
            outcome.assembly.reference_bases(),
            outcome.assembly.padded_reference_start_1based(),
            "20",
            outcome.assembly.max_mnp_distance(),
        );
        let ref_hap = outcome.assembly.haplotypes.iter().find(|h| h.is_reference);
        let apply_bases = ref_hap
            .map(|h| h.bases.as_slice())
            .unwrap_or_else(|| outcome.assembly.reference_bases());
        let apply_pad = ref_hap
            .and_then(|h| h.genome_loc.map(|g| g.start_1based()))
            .unwrap_or_else(|| outcome.assembly.padded_reference_start_1based());
        let mapping = crate::hc_allele_mapping::create_allele_mapper(
            &event,
            pos,
            &outcome.assembly.haplotypes,
            apply_pad,
            apply_bases,
            outcome.assembly.max_mnp_distance(),
            false,
        );
        let mapping_full = crate::hc_allele_mapping::create_allele_mapper(
            &event,
            pos,
            &outcome.assembly.haplotypes,
            outcome.assembly.padded_reference_start_1based(),
            outcome.assembly.reference_bases(),
            outcome.assembly.max_mnp_distance(),
            false,
        );
        let nearby: Vec<_> = outcome
            .assembly
            .variation_events()
            .iter()
            .filter(|e| e.start_1based.get().abs_diff(pos) <= 200)
            .map(|e| format!("{}:{}>{}", e.start_1based.get(), e.ref_allele, e.alt_allele))
            .collect();
        let call_dbg = outcome
            .genotyped_calls
            .iter()
            .find(|c| c.event.start_1based == target)
            .map(|c| {
                format!(
                    "AD={:?} PL={:?} GQ={}",
                    c.genotype.format.ad_as_i32(),
                    c.genotype.format.pl_as_i32(),
                    c.genotype.format.gq.as_i32()
                )
            })
            .unwrap_or_else(|| "NO_CALL".into());
        eprintln!(
            "L9-FN {pos} {ref_a}>{alt_a}: region={}:{}-{} reads={} pileup_AD={rr},{ra} \
                 on_cigar={on_cigar} in_event_map={in_events} in_genotyped_calls={in_calls} \
                 mapper_trim_alts={} mapper_full_alts={} call={call_dbg} nearby_events={nearby:?}",
            region.contig,
            region.start.get(),
            region.end.get(),
            region.reads.len(),
            mapping.alt_haplotype_indices.len(),
            mapping_full.alt_haplotype_indices.len(),
        );
        assert!(ra >= 2, "{pos}: expected pileup alt>=2, got {rr},{ra}");
        assert!(
            on_cigar && in_events,
            "L9-FN {pos}: expected EventMap/CIGAR discovery"
        );
        assert!(
                in_calls,
                "L9-FN {pos}: expected genotyped call after L9 pileup rescue / Class-A3; got {call_dbg}"
            );
    }
}

/// Diagnose Class-A2 GT-flip sites: why PairHMM informative AD collapses to alt
/// (pre-Class-A2 PL ~1/1) despite balanced pileup. Prefer 10009227.
#[test]
fn class_a2_hap_pool_collapse_probe() {
    use crate::engine::{CallRegionArgs, HaplotypeCallerEngine};
    use crate::event_map::variation_events_for_haplotype;
    use crate::hc_allele_mapping::{create_allele_mapper, hap_base_at_ref_locus};
    use crate::hc_genotyping_engine::{
        alt_hap_indices_for_genotype_marginalization, biallelic_allele_depths_from_rows,
        marginalize_rows_to_biallelic_alleles, ref_hap_indices_for_genotype_marginalization,
        region_likelihoods_to_rows,
    };
    use crate::read_model::ReadFilterParams;
    use crate::walker_traversal::{flatten_assembly_regions, traverse_assembly_region_walker};
    use crate::{call_disposition, AssemblyRegionCallDisposition, WalkerTraversalConfig};
    use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};

    if !dense_bam().is_file() || !dense_ref().is_file() {
        eprintln!("skip: dense assets missing");
        return;
    }
    let ref_fa = dense_ref();
    let bam = dense_bam();
    let dict = SequenceDictionary::from_fasta_path(&ref_fa).expect("dict");
    let filters = ReadFilterParams::gatk_standard_hc();
    let cfg = WalkerTraversalConfig::gatk_haplotype_caller_production(100);
    let args = CallRegionArgs::strict_java();
    let gt_cfg = args.genotyping.clone();

    for (interval, pos, ref_a, alt_a) in [
        ("20:10009000-10009400", 10009227u64, "A", "G"),
        ("20:10012200-10012500", 10012384, "T", "C"),
        ("20:10012450-10012800", 10012636, "G", "C"),
    ] {
        let specs = parse_intervals_cli_string(&dict, interval).expect("interval");
        let walk = traverse_assembly_region_walker(&dict, &specs, &ref_fa, &bam, &filters, &cfg)
            .expect("walk");
        let regions = flatten_assembly_regions(&walk);
        let Some(region) = regions.iter().find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= pos
                && r.end.get() >= pos
        }) else {
            panic!("L9-A2 {pos}: NO active region covering site in {interval}");
        };

        let outcome = HaplotypeCallerEngine::call_region(region, &dict, &ref_fa, &args)
            .expect("call_region")
            .expect("outcome");
        let event = VariationEvent {
            contig: "20".into(),
            start_1based: GenomePosition::new_1based(pos),
            end_1based: GenomePosition::new_1based(pos),
            ref_allele: ref_a.into(),
            alt_allele: alt_a.into(),
        };
        let ref_hap = outcome
            .assembly
            .haplotypes
            .iter()
            .find(|h| h.is_reference)
            .expect("ref hap");
        let apply_bases = ref_hap.bases.as_slice();
        let apply_pad = ref_hap
            .genome_loc
            .map(|g| g.start_1based())
            .unwrap_or_else(|| outcome.assembly.padded_reference_start_1based());
        let full_pad = outcome.assembly.padded_reference_start_1based();
        let max_mnp = outcome.assembly.max_mnp_distance();

        let mapping = create_allele_mapper(
            &event,
            pos,
            &outcome.assembly.haplotypes,
            apply_pad,
            apply_bases,
            max_mnp,
            false,
        );
        let ref_pool = ref_hap_indices_for_genotype_marginalization(
            &mapping,
            &outcome.assembly.haplotypes,
            &gt_cfg,
            Some(&event),
        );
        let alt_pool = alt_hap_indices_for_genotype_marginalization(
            &mapping,
            &outcome.assembly.haplotypes,
            &event,
            ref_hap,
            apply_pad,
            apply_bases,
            max_mnp,
            "20",
            &gt_cfg,
        );

        for (i, h) in outcome.assembly.haplotypes.iter().enumerate() {
            let base = hap_base_at_ref_locus(h, apply_pad, pos)
                .map(|b| (b as char).to_ascii_uppercase())
                .unwrap_or('?');
            let nearby_full: Vec<_> = variation_events_for_haplotype(
                h,
                ref_hap,
                outcome.assembly.reference_bases(),
                full_pad,
                max_mnp,
                "20",
            )
            .into_iter()
            .filter(|e| e.start_1based.get().abs_diff(pos) <= 50)
            .map(|e| format!("{}:{}>{}", e.start_1based.get(), e.ref_allele, e.alt_allele))
            .collect();
            eprintln!(
                    "L9-A2 {pos} hap[{i}] ref={} base={base} align0={} pool_ref={} pool_alt={} nearby={nearby_full:?}",
                    h.is_reference,
                    h.alignment_start_hap_wrt_ref,
                    mapping.ref_haplotype_indices.iter().any(|x| x.get() == i),
                    mapping.alt_haplotype_indices.iter().any(|x| x.get() == i),
                );
        }

        let pileup_region = read_allele_depths_at_locus(&region.reads, &event, 1);
        let rows = region_likelihoods_to_rows(
            &outcome.read_likelihoods,
            outcome.assembly.haplotypes.len(),
        );
        let marg = marginalize_rows_to_biallelic_alleles(&rows, &ref_pool, &alt_pool);
        let info_ad = biallelic_allele_depths_from_rows(&marg, 0, 1);
        let call = outcome
            .genotyped_calls
            .iter()
            .find(|c| {
                c.event.start_1based == GenomePosition::new_1based(pos)
                    && c.event.ref_allele == ref_a
                    && c.event.alt_allele == alt_a
            })
            .expect("expected genotyped call");
        let pl = call.genotype.format.pl_as_i32();
        let ad = call.genotype.format.ad_as_i32();
        eprintln!(
            "L9-A2 {pos} {ref_a}>{alt_a}: mapper_ref={:?} mapper_alt={:?} \
                 pileup={:?} informative_AD={info_ad:?} emitted AD={ad:?} PL={pl:?}",
            mapping.ref_haplotype_indices, mapping.alt_haplotype_indices, pileup_region,
        );

        assert!(
            pileup_region.0 >= 8 && pileup_region.1 >= 8,
            "expected balanced pileup at {pos}, got {:?}",
            pileup_region
        );
        assert!(
            info_ad[0] >= 5 && info_ad[1] >= 5,
            "L9-A2 {pos}: expected balanced informative AD after pad fix, got {info_ad:?}"
        );
        assert_eq!(pl.len(), 3, "biallelic PL");
        assert_eq!(pl[1], 0, "{pos}: het PL must be best; got {pl:?}");
        assert!(
            pl != [81, 0, 36],
            "{pos}: must not fall back to SparsePlShape Het; got {pl:?}"
        );
        assert!(
            pl[0] >= 100 && pl[2] >= 100,
            "{pos}: expected Java-scale PairHMM het PLs, got {pl:?}"
        );
    }
}
