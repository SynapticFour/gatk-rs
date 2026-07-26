//! P3.6: interior RCM uses post-realign genotyping evidence pileup (Java getPileupsOverReference).
//! Run: `P12_REFERENCE=parity/realworld/assets/hs37d5.simple.fa cargo test -p gatk-haplotypecaller p12_rcm_interior --release -- --ignored --nocapture`

use gatk_core::reference::{parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary};
use gatk_haplotypecaller::locus_iterator::LocusPileupState;
use gatk_haplotypecaller::pileup_element::pileup_element_flags_at_ref;
use gatk_haplotypecaller::read_model::ReadFilterParams;
use gatk_haplotypecaller::read_projection::query_index_at_reference_position;
use gatk_haplotypecaller::ref_confidence::{
    reads_overlap_closed_span, reference_confidence_loci_for_active_region, ClusterRcmEvidenceMode,
    ReferenceConfidenceConfig,
};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions,
    reference_vcf_emit::{emitted_variant_starts_in_region, first_emitted_variant_start_in_region},
    traverse_assembly_region_walker, AssemblyRegionCallDisposition, CallRegionArgs,
    HaplotypeCallerEngine, WalkerTraversalConfig,
};
use rust_htslib::bam::Read;
use std::path::Path;

const INTERIOR_POS: u64 = 92305671;

fn query_index_for_ref_offset(
    alignment_start: i64,
    cigar: &rust_htslib::bam::record::CigarString,
    ref_pos0: i64,
) -> Option<usize> {
    query_index_at_reference_position(alignment_start, cigar, ref_pos0)
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
    if ref_path.is_file() && bam.is_file() {
        Some((ref_path, bam))
    } else {
        None
    }
}

#[test]
#[ignore = "P12 BAM: interior RCM after Java-aligned realign"]
fn p12_rcm_interior_genotyping_vs_region_pileup() {
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
            ) && r.start.get() <= INTERIOR_POS
                && r.end.get() >= INTERIOR_POS
        })
        .expect("active region");

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

    let gt_overlaps = reads_overlap_closed_span(
        &outcome.genotyping_reads,
        region.start.get(),
        region.end.get(),
    );
    let first_variant = first_emitted_variant_start_in_region(region, &outcome, 10.0)
        .expect("first emitted variant");
    let emitted_variants =
        emitted_variant_starts_in_region(region, &outcome, 10.0).expect("emitted variants");

    let gt_filters = ReadFilterParams::genotyping_evidence_rcm_pileup();
    let mut ref_cache_diag = ReferenceWindowCache::new(ref_fasta.clone(), 4);
    let ref_bytes = ref_cache_diag
        .get_interval_bytes(&dict, &region.contig, region.start.get(), region.end.get())
        .expect("ref bytes");
    let ref_base = *ref_bytes
        .get((INTERIOR_POS - region.start.get()) as usize)
        .unwrap_or(&b'N');
    let mut gt_state = LocusPileupState::from_genotyping_evidence_records(
        &outcome.genotyping_reads,
        &header,
        &region.contig,
        &gt_filters,
    );
    let gt_only = gt_state
        .pileup_at(
            &outcome.genotyping_reads,
            &gt_filters,
            INTERIOR_POS,
            ref_base,
        )
        .expect("gt pileup");

    for (i, h) in outcome.assembly.haplotypes.iter().enumerate() {
        eprintln!(
            "hap{i} ref={} align_start={} len={} cigar={}",
            h.is_reference,
            h.alignment_start_hap_wrt_ref,
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
    let hap2 = &outcome.assembly.haplotypes[2];
    for off in [0usize, 40, 56, 60] {
        let rb = ref_hap.bases.get(off).copied().unwrap_or(b'?') as char;
        let h2 = hap2.bases.get(off).copied().unwrap_or(b'?') as char;
        eprintln!("hap_base_offset {off}: ref_hap={rb} hap2={h2}");
    }

    let hap_priorities: Vec<f64> = outcome
        .assembly
        .haplotypes
        .iter()
        .map(|h| {
            let reference_term = if h.is_reference { 1.0 } else { 0.0 };
            let cigar_term = h
                .cigar
                .as_ref()
                .map(|c| 1.0 - c.elements.len() as f64)
                .unwrap_or(0.0);
            reference_term + cigar_term
        })
        .collect();
    let best_hap_for_read = |ri: usize| -> usize {
        const INFORMATIVE: f64 = 0.2;
        let ll = |hi: usize| -> f64 {
            outcome
                .read_likelihoods
                .iter()
                .find(|r| r.read_index.get() == ri && r.haplotype_index.get() == hi)
                .map(|r| r.log10_likelihood)
                .unwrap_or(f64::NEG_INFINITY)
        };
        let hap_count = outcome.assembly.haplotypes.len();
        let mut best_a = 0usize;
        let mut second_a = 0usize;
        let mut best_ll = ll(0);
        let mut second_ll = f64::NEG_INFINITY;
        for a in 1..hap_count {
            let candidate = ll(a);
            if candidate > best_ll {
                second_a = best_a;
                second_ll = best_ll;
                best_a = a;
                best_ll = candidate;
            } else if candidate > second_ll {
                second_a = a;
                second_ll = candidate;
            }
        }
        if best_ll - second_ll < INFORMATIVE {
            let mut tie_best = best_a;
            let mut tie_best_pri = hap_priorities.get(best_a).copied().unwrap_or(0.0);
            let mut tie_second_pri = hap_priorities.get(second_a).copied().unwrap_or(0.0);
            for a in 0..hap_count {
                let candidate = ll(a);
                if a == best_a || best_ll - candidate > INFORMATIVE {
                    continue;
                }
                let pri = hap_priorities.get(a).copied().unwrap_or(0.0);
                if pri > tie_best_pri {
                    tie_second_pri = tie_best_pri;
                    tie_best = a;
                    tie_best_pri = pri;
                } else if pri > tie_second_pri {
                    tie_second_pri = pri;
                }
            }
            let _ = tie_second_pri;
            best_a = tie_best;
        }
        best_a
    };

    for (i, r) in outcome.genotyping_reads.iter().enumerate() {
        let pos_1 = r.pos().max(0) as u64 + 1;
        let cigar: Vec<_> = r.cigar().iter().copied().collect();
        let flags = pileup_element_flags_at_ref(
            r.pos(),
            &cigar,
            &r.seq().as_bytes(),
            r.qual(),
            INTERIOR_POS.saturating_sub(1) as i64,
        );
        let gt_base = flags.as_ref().map(|f| f.read_base as char).unwrap_or('?');
        let region_r = region.reads.iter().find(|rr| rr.qname() == r.qname());
        let region_base = region_r
            .and_then(|rr| {
                let c: Vec<_> = rr.cigar().iter().copied().collect();
                pileup_element_flags_at_ref(
                    rr.pos(),
                    &c,
                    &rr.seq().as_bytes(),
                    rr.qual(),
                    INTERIOR_POS.saturating_sub(1) as i64,
                )
            })
            .map(|f| f.read_base as char)
            .unwrap_or('?');
        let best_hi = best_hap_for_read(i);
        let region_cigar = region_r
            .map(|rr| format!("{:?}", rr.cigar().iter().collect::<Vec<_>>()))
            .unwrap_or_else(|| "?".into());
        eprintln!(
            "read{i} {} best_hap={best_hi} ref={} pos={pos_1} region_base={region_base} gt_base={gt_base} region_cigar={region_cigar} gt_cigar={cigar:?}",
            String::from_utf8_lossy(r.qname()),
            outcome.assembly.haplotypes[best_hi].is_reference,
        );
        if let Some(rr) = region_r {
            let orig_pos = rr.pos().max(0) as u64 + 1;
            let orig_off = INTERIOR_POS.saturating_sub(orig_pos) as usize;
            let gt_off = INTERIOR_POS.saturating_sub(pos_1) as usize;
            let seq = r.seq().as_bytes();
            let (orig_qi, gt_qi) = (
                query_index_for_ref_offset(
                    rr.pos(),
                    &rr.cigar(),
                    INTERIOR_POS.saturating_sub(1) as i64,
                ),
                query_index_for_ref_offset(
                    r.pos(),
                    &r.cigar(),
                    INTERIOR_POS.saturating_sub(1) as i64,
                ),
            );
            eprintln!(
                "  read{i} orig_off={orig_off} gt_off={gt_off} orig_qi={orig_qi:?} gt_qi={gt_qi:?} base@orig_qi={} base@gt_qi={}",
                orig_qi.and_then(|qi| seq.get(qi)).map(|b| *b as char).unwrap_or('?'),
                gt_qi.and_then(|qi| seq.get(qi)).map(|b| *b as char).unwrap_or('?'),
            );
            for hi in 0..outcome.assembly.haplotypes.len().min(4) {
                let ll = outcome
                    .read_likelihoods
                    .iter()
                    .find(|x| x.read_index.get() == i && x.haplotype_index.get() == hi)
                    .map(|x| x.log10_likelihood)
                    .unwrap_or(f64::NEG_INFINITY);
                eprintln!("  read{i} hap{hi} log10_ll={ll}");
            }
        }
    }

    let filters = ReadFilterParams::gatk_standard_hc();
    let mut region_state =
        LocusPileupState::from_records(&region.reads, &header, &region.contig, &filters);
    let region_only = region_state
        .pileup_at(&region.reads, &filters, INTERIOR_POS, ref_base)
        .expect("region pileup");
    eprintln!(
        "region.reads pileup dp={} at {INTERIOR_POS}",
        region_only.len()
    );
    for (i, obs) in region_only.iter().enumerate() {
        let base = if obs.is_deletion {
            '-'
        } else {
            obs.read_base as char
        };
        eprintln!(
            "region_pileup[{i}] base={base} ref={} is_alt={}",
            ref_base as char, obs.is_alt
        );
    }

    let mut ref_cache = ReferenceWindowCache::new(ref_fasta.clone(), 4);
    let config = ReferenceConfidenceConfig::default();
    let loci = reference_confidence_loci_for_active_region(
        region,
        &outcome.genotyping_reads,
        first_variant,
        &emitted_variants,
        &header,
        &config,
        &filters,
        &mut ref_cache,
        &dict,
        ClusterRcmEvidenceMode::Production,
    )
    .expect("loci");
    let interior = loci
        .iter()
        .find(|l| l.position_1based as u64 == INTERIOR_POS)
        .expect("interior locus");
    let ref_char = ref_base as char;
    for (i, obs) in gt_only.iter().enumerate() {
        let base = if obs.is_deletion {
            '-'
        } else {
            obs.read_base as char
        };
        eprintln!(
            "gt_pileup[{i}] base={base} ref={ref_char} qual={} is_alt={}",
            obs.qual, obs.is_alt
        );
    }
    let hom_ref = gt_only.iter().filter(|o| !o.is_alt).count();
    eprintln!(
        "region {}:{}-{} gt_overlaps={gt_overlaps} gt_only_dp={} hom_ref={hom_ref} interior gq={} dp={} (Java gq=12 dp=4)",
        region.contig, region.start.get(), region.end.get(), gt_only.len(), interior.gq, interior.dp
    );

    assert!(gt_overlaps, "P3.5: genotyping reads overlap assembly span");
    assert!(
        gt_only.len() >= 4,
        "P3.6: genotyping pileup depth (Java MIN_DP~4); got {}",
        gt_only.len()
    );
    assert_eq!(
        interior.dp,
        gt_only.len() as i32,
        "interior RCM must use genotyping evidence pileup, not region.reads fallback"
    );
    // P3.7: Java hom-ref interior block GQ=12 MIN_DP=4 at 92305671.
    assert!(
        hom_ref >= 3,
        "P3.7: post-realign pileup should be mostly hom-ref at interior; got hom_ref={hom_ref}/{}",
        gt_only.len()
    );
    assert!(
        interior.gq >= 10,
        "P3.7: interior GQ should match Java ~12; got {}",
        interior.gq
    );
}
