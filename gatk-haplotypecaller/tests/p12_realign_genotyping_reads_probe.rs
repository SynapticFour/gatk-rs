//! P3.5 / P1.3: post-realign genotyping read positions vs assembly span + cluster tail RCM.
//! Run: `P12_REFERENCE=parity/realworld/assets/hs37d5.simple.fa cargo test -p gatk-haplotypecaller p12_realign_genotyping --release -- --ignored --nocapture`

use gatk_core::reference::{parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary};
use gatk_haplotypecaller::locus_iterator::LocusPileupState;
use gatk_haplotypecaller::ref_confidence::{
    reads_overlap_closed_span, reference_confidence_loci_for_active_region, ClusterRcmEvidenceMode,
    ReferenceConfidenceConfig,
};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    AssemblyRegionCallDisposition, CallRegionArgs, HaplotypeCallerEngine, ReadFilterParams,
    WalkerTraversalConfig,
};
use rust_htslib::bam::Read;
use std::path::Path;

const PROBE_POS: u64 = 92305671;
const TAIL_START: u64 = 92305699;
const TAIL_END: u64 = 92305754;
const TAIL_ANCHOR: u64 = 92305754;

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
#[ignore = "P12 BAM: realign genotyping read span probe"]
fn p12_realign_genotyping_reads_overlap_assembly_span() {
    let Some((ref_fasta, bam)) = p12_paths() else {
        eprintln!("skip: P12_REFERENCE / BAM");
        return;
    };

    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, "2:92305500-92305720").expect("interval");
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
            ) && r.start.get() <= PROBE_POS
                && r.end.get() >= PROBE_POS
        })
        .expect("active region");

    let outcome = HaplotypeCallerEngine::call_region(
        region,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("outcome");

    let pad = outcome.assembly.padded_reference_start_1based();
    eprintln!(
        "region {}:{}-{} pad_start={pad}",
        region.contig,
        region.start.get(),
        region.end.get()
    );
    for (i, h) in outcome.assembly.haplotypes.iter().enumerate() {
        eprintln!(
            "hap{i}\tref={}\talign_start={}\tlen={}\tcigar={}",
            h.is_reference,
            h.alignment_start_hap_wrt_ref,
            h.bases.len(),
            h.cigar
                .as_ref()
                .map(|c| c.to_gatk_string())
                .unwrap_or_default()
        );
    }
    for (i, r) in outcome.genotyping_reads.iter().enumerate() {
        let pos_1 = r.pos().max(0) as u64 + 1;
        let end_1 = r.cigar().end_pos().max(0) as u64;
        eprintln!(
            "gt_read{i}\tpos_1based={pos_1}\tend_1based={end_1}\tcigar={:?}",
            r.cigar().iter().collect::<Vec<_>>()
        );
    }
    let overlaps = reads_overlap_closed_span(
        &outcome.genotyping_reads,
        region.start.get(),
        region.end.get(),
    );
    eprintln!("gt_overlaps={overlaps}");
    assert!(
        overlaps,
        "genotyping reads should overlap assembly span after Java-aligned realign"
    );
}

#[test]
#[ignore = "P12 BAM: cluster tail RCM after realign (92305699-92305754)"]
fn p12_realign_cluster_tail_rcm_probe() {
    let Some((ref_fasta, bam)) = p12_paths() else {
        eprintln!("skip: P12_REFERENCE / BAM");
        return;
    };

    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, "2:92305500-92305800").expect("interval");
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
            ) && r.start.get() <= TAIL_ANCHOR
                && r.end.get() >= TAIL_START
        })
        .expect("active region covering cluster tail");

    let header = rust_htslib::bam::Reader::from_path(&bam)
        .expect("bam")
        .header()
        .clone();

    let outcome = HaplotypeCallerEngine::call_region(
        region,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("outcome");

    assert!(
        reads_overlap_closed_span(&outcome.genotyping_reads, TAIL_START, TAIL_END),
        "genotyping reads must cover cluster tail band"
    );

    let emitted: Vec<u64> = outcome
        .genotyped_calls
        .iter()
        .map(|c| c.event.start_1based.get())
        .collect();
    let first_variant = emitted.iter().copied().min();

    let filters = ReadFilterParams::genotyping_evidence_rcm_pileup();
    let mut ref_cache = ReferenceWindowCache::new(ref_fasta.clone(), 4);
    let loci = reference_confidence_loci_for_active_region(
        region,
        &outcome.genotyping_reads,
        first_variant,
        &emitted,
        &header,
        &ReferenceConfidenceConfig::default(),
        &filters,
        &mut ref_cache,
        &dict,
        ClusterRcmEvidenceMode::Production,
    )
    .expect("loci");

    let locus_at = |pos: u64| {
        loci.iter()
            .find(|l| l.position_1based as u64 == pos)
            .unwrap_or_else(|| panic!("missing RCM locus at {pos}"))
    };

    let gt_filters = ReadFilterParams::genotyping_evidence_rcm_pileup();
    let mut gt_state = LocusPileupState::from_genotyping_evidence_records(
        &outcome.genotyping_reads,
        &header,
        &region.contig,
        &gt_filters,
    );
    let ref_cache_bytes = ref_cache
        .get_interval_bytes(&dict, &region.contig, TAIL_START, TAIL_END)
        .expect("ref bytes");

    for (offset, pos) in (TAIL_START..=TAIL_END).enumerate() {
        let ref_base = *ref_cache_bytes.get(offset).unwrap_or(&b'N');
        let pileup = gt_state
            .pileup_at(&outcome.genotyping_reads, &gt_filters, pos, ref_base)
            .expect("pileup");
        if pos == 92305729 || pos == TAIL_ANCHOR {
            eprintln!(
                "tail pos={pos} gt_dp={} hom_ref={} locus_gq={} locus_dp={}",
                pileup.len(),
                pileup.iter().filter(|o| !o.is_alt).count(),
                locus_at(pos).gq,
                locus_at(pos).dp
            );
        }
    }

    let l_92305729 = locus_at(92305729);
    assert!(
        l_92305729.gq >= 5 && l_92305729.gq <= 7,
        "92305729 Java GQ=6; got {}",
        l_92305729.gq
    );
    assert!(
        l_92305729.dp >= 2,
        "92305729 Java MIN_DP=2; got dp={}",
        l_92305729.dp
    );

    let l_tail = locus_at(TAIL_ANCHOR);
    assert_eq!(l_tail.gq, 0, "92305754 Java GQ=0 shadow");
    assert!(
        l_tail.dp >= 2,
        "92305754 Java MIN_DP=3; got dp={}",
        l_tail.dp
    );
}
