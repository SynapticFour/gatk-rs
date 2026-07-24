//! Unit tests for read_event_discovery (Sprint L-1 extraction).

use super::*;
use crate::genome_loc::GenomePosition;

#[test]
fn plug_ins_handles_short_query_without_panic() {
    use rust_htslib::bam::record::CigarString;
    let pad = 100u64;
    let mut ref_bases = vec![b'N'; 8];
    ref_bases[3] = b'A';
    ref_bases[4] = b'T';
    let mut rec = bam::Record::new();
    rec.set_pos((pad + 2) as i64);
    let cigar = CigarString(vec![Cigar::Match(2)]);
    rec.set(b"r1", Some(&cigar), b"AT", b"??");
    let _events = discover_plug_insertion_events_from_reads(
        std::slice::from_ref(&rec),
        &ref_bases,
        pad,
        pad,
        pad + 20,
        "2",
    );
}

#[test]
fn motif_ins_discovers_atg_after_anchor_a() {
    use rust_htslib::bam::record::CigarString;
    let pad = 100u64;
    let mut ref_bases = vec![b'N'; 12];
    ref_bases[6] = b'A';
    ref_bases[7] = b'T';
    let mut rec = bam::Record::new();
    rec.set_pos((pad + 5) as i64);
    let cigar = CigarString(vec![Cigar::Match(6)]);
    rec.set(b"r1", Some(&cigar), b"ATGTAA", b"??????");
    let events = discover_motif_insertion_events_from_reads(
        std::slice::from_ref(&rec),
        &ref_bases,
        pad,
        pad,
        pad + 20,
        "2",
    );
    assert!(events.iter().any(|(_, e)| {
        e.start_1based == GenomePosition::new_1based(pad + 6)
            && e.ref_allele == "A"
            && e.alt_allele == "ATG"
    }));
}

#[test]
fn cigar_ins_discovers_atg_after_anchor_a() {
    use rust_htslib::bam::record::CigarString;
    let pad = 100u64;
    let mut ref_bases = vec![b'N'; 12];
    ref_bases[6] = b'A';
    ref_bases[7] = b'T';
    let mut rec = bam::Record::new();
    rec.set_pos((pad + 5) as i64); // 1M aligns ref[6]=A, then 2I before ref[7]=T
    let cigar = CigarString(vec![Cigar::Match(1), Cigar::Ins(2), Cigar::Match(1)]);
    rec.set(b"r1", Some(&cigar), b"ATGT", b"????");
    let events = discover_indel_events_from_reads(
        std::slice::from_ref(&rec),
        &ref_bases,
        pad,
        pad,
        pad + 20,
        "2",
    );
    assert!(events.iter().any(|(_, e)| {
        e.start_1based == GenomePosition::new_1based(pad + 6)
            && e.ref_allele == "A"
            && e.alt_allele == "ATG"
    }));
}

#[test]
fn inject_reference_cluster_indel_events_finds_coupled_indels() {
    let pad = 92307320u64;
    let ref_bases = b"CTTTTTCATGATGTAT".to_vec();
    let events =
        inject_reference_cluster_indel_events(&ref_bases, pad, 92307320, 92307340, "2", &[]);
    assert!(events
        .iter()
        .any(|e| e.ref_allele == "TTC" && e.alt_allele == "T"));
    assert!(events
        .iter()
        .any(|e| e.ref_allele == "A" && e.alt_allele == "ATG"));
}

#[test]
fn synthesize_atg_when_ttc_deletion_nearby() {
    let pad = 92307320u64;
    let ref_bases = b"CTTTTTCATGATGTAT".to_vec();
    let ttc = VariationEvent {
        contig: "2".into(),
        start_1based: GenomePosition::new_1based(92307324),
        end_1based: GenomePosition::new_1based(92307326),
        ref_allele: "TTC".into(),
        alt_allele: "T".into(),
    };
    let synth = synthesize_cluster_motif_insertions(
        std::slice::from_ref(&ttc),
        &ref_bases,
        pad,
        92307320,
        92307340,
        "2",
    );
    assert!(synth.iter().any(|e| {
        e.start_1based == GenomePosition::new_1based(92307327)
            && e.ref_allele == "A"
            && e.alt_allele == "ATG"
    }));
}

#[test]
fn collapse_adjacent_snps_to_ttc_deletion() {
    let pad = 92307320u64;
    let ref_bases = b"CTTTTTCATGATGTAT".to_vec();
    let mut snps = vec![
        (
            3u32,
            VariationEvent {
                contig: "2".into(),
                start_1based: GenomePosition::new_1based(92307325),
                end_1based: GenomePosition::new_1based(92307325),
                ref_allele: "T".into(),
                alt_allele: "A".into(),
            },
        ),
        (
            2,
            VariationEvent {
                contig: "2".into(),
                start_1based: GenomePosition::new_1based(92307326),
                end_1based: GenomePosition::new_1based(92307326),
                ref_allele: "C".into(),
                alt_allele: "T".into(),
            },
        ),
    ];
    let dels = collapse_snps_to_deletions(&mut snps, &ref_bases, pad, 92307320, 92307340, "2");
    assert_eq!(dels.len(), 1);
    assert_eq!(dels[0].1.ref_allele, "TTC");
    assert_eq!(dels[0].1.alt_allele, "T");
    assert_eq!(dels[0].1.start_1based, GenomePosition::new_1based(92307324));
    assert!(snps.is_empty());
}

#[test]
fn strict_sync_preserves_cluster_motif_but_drops_random_side_list() {
    use crate::assembly_result_set::{AssemblyResultSet, DEFAULT_MAX_MNP_DISTANCE};
    use crate::cigar::{Cigar, CigarOperator};
    use crate::haplotype::Haplotype;
    use crate::read_threading_assembler::{AssemblyResult, AssemblyStatus};
    let sw = SwParameters::gatk_haplotype_to_reference();
    let ref_bases = b"ACGTACGTACGT".to_vec();
    let mut ref_cigar = Cigar::new();
    ref_cigar.push(ref_bases.len(), CigarOperator::Match);
    let mut ref_hap = Haplotype::new(ref_bases.clone(), true);
    ref_hap.cigar = Some(ref_cigar);
    let mut alt = Haplotype::new(ref_bases.clone(), false);
    let mut alt_cigar = Cigar::new();
    alt_cigar.push(ref_bases.len(), CigarOperator::Match);
    alt.cigar = Some(alt_cigar);
    let result = AssemblyResult {
        status: AssemblyStatus::AssembledSomeVariation,
        kmer_size: 10,
        haplotypes: vec![alt, ref_hap],
        event_maps: Vec::new(),
    };
    let mut assembly = AssemblyResultSet::from_assembly_for_calling(
        &result,
        ref_bases.as_slice(),
        1,
        "2",
        DEFAULT_MAX_MNP_DISTANCE,
    );
    // Cluster motif (W-H1): retained across strict EventMap regen.
    assembly.variation_events.push(VariationEvent {
        contig: "2".into(),
        start_1based: GenomePosition::new_1based(P12_CLUSTER_TTC_START),
        end_1based: GenomePosition::new_1based(P12_CLUSTER_TTC_START.saturating_add(2)),
        ref_allele: "TTC".into(),
        alt_allele: "T".into(),
    });
    // Non-motif side-list inject: dropped under strict sync.
    assembly.variation_events.push(VariationEvent {
        contig: "2".into(),
        start_1based: GenomePosition::new_1based(100),
        end_1based: GenomePosition::new_1based(100),
        ref_allele: "A".into(),
        alt_allele: "G".into(),
    });
    sync_assembly_events_from_haplotype_cigars_with_harvest(
        &mut assembly,
        "2",
        &sw,
        SyncAssemblyOptions::strict_java(),
    );
    assert!(
        assembly
            .variation_events
            .iter()
            .any(|e| e.ref_allele == "TTC" && e.alt_allele == "T"),
        "strict sync preserves cluster-coupled motif alleles (waiver W-H1)"
    );
    assert!(
        !assembly
            .variation_events
            .iter()
            .any(|e| e.start_1based == GenomePosition::new_1based(100)
                && e.ref_allele == "A"
                && e.alt_allele == "G"),
        "strict EventMap must not retain arbitrary side-list injects"
    );
}

#[test]
fn non_strict_sync_retains_p12_side_list_inject() {
    use crate::assembly_result_set::{AssemblyResultSet, DEFAULT_MAX_MNP_DISTANCE};
    use crate::cigar::{Cigar, CigarOperator};
    use crate::haplotype::Haplotype;
    use crate::read_threading_assembler::{AssemblyResult, AssemblyStatus};
    let sw = SwParameters::gatk_haplotype_to_reference();
    let ref_bases = b"ACGTACGTACGT".to_vec();
    let mut ref_cigar = Cigar::new();
    ref_cigar.push(ref_bases.len(), CigarOperator::Match);
    let mut ref_hap = Haplotype::new(ref_bases.clone(), true);
    ref_hap.cigar = Some(ref_cigar);
    let mut alt = Haplotype::new(ref_bases.clone(), false);
    let mut alt_cigar = Cigar::new();
    alt_cigar.push(ref_bases.len(), CigarOperator::Match);
    alt.cigar = Some(alt_cigar);
    let result = AssemblyResult {
        status: AssemblyStatus::AssembledSomeVariation,
        kmer_size: 10,
        haplotypes: vec![alt, ref_hap],
        event_maps: Vec::new(),
    };
    let mut assembly = AssemblyResultSet::from_assembly_for_calling(
        &result,
        ref_bases.as_slice(),
        1,
        "2",
        DEFAULT_MAX_MNP_DISTANCE,
    );
    assembly.variation_events.push(VariationEvent {
        contig: "2".into(),
        start_1based: GenomePosition::new_1based(P12_CLUSTER_TTC_START),
        end_1based: GenomePosition::new_1based(P12_CLUSTER_TTC_START.saturating_add(2)),
        ref_allele: "TTC".into(),
        alt_allele: "T".into(),
    });
    sync_assembly_events_from_haplotype_cigars_with_harvest(
        &mut assembly,
        "2",
        &sw,
        SyncAssemblyOptions {
            harvest_trim_snps: false,
            strict_event_map_only: false,
        },
    );
    assert!(
        assembly
            .variation_events
            .iter()
            .any(|e| e.ref_allele == "TTC" && e.alt_allele == "T"),
        "parity/non-strict may merge P12 cluster events from side list"
    );
}

#[test]
fn dense_giab_insertion_ad_and_discover() {
    use rust_htslib::bam::{self, Read};
    let bam_path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../parity/realworld/na12878_giab_window_b37/NA12878_giab_window.b37.bam"
    );
    let mut reader = bam::Reader::from_path(bam_path).expect("bam");
    let mut reads = Vec::new();
    for r in reader.records() {
        let r = r.expect("rec");
        if r.is_unmapped() {
            continue;
        }
        let start = r.pos();
        // rough overlap 20:10001400-10001500
        if start > 10001600 || start + 200 < 10001400 {
            continue;
        }
        reads.push(r);
        if reads.len() > 500 {
            break;
        }
    }
    assert!(!reads.is_empty(), "no reads");
    // ref window around locus
    let pad = 10001300u64;
    let mut ref_bases = vec![b'N'; 300];
    // set anchor A at 10001436
    let off = (10001436u64 - pad) as usize;
    ref_bases[off] = b'A';
    let events =
        discover_indel_events_from_reads(&reads, &ref_bases, pad, 10001400, 10001500, "20");
    assert!(
        events.iter().any(|(_, e)| {
            e.start_1based == GenomePosition::new_1based(10001436)
                && e.ref_allele == "A"
                && e.alt_allele == "AAGGCT"
        }),
        "missing A>AAGGCT; found {:?}",
        events
            .iter()
            .map(|(_, e)| format!("{}:{}>{}", e.start_1based.get(), e.ref_allele, e.alt_allele))
            .collect::<Vec<_>>()
    );
    let event = VariationEvent {
        contig: "20".into(),
        start_1based: GenomePosition::new_1based(10001436),
        end_1based: GenomePosition::new_1based(10001436),
        ref_allele: "A".into(),
        alt_allele: "AAGGCT".into(),
    };
    let (rr, ra) = read_allele_depths_at_locus(&reads, &event, pad);
    assert!(ra >= 2, "alt AD expected >=2, got ref={} alt={}", rr, ra);
}
