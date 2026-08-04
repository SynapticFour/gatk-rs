//! Gold: Java finalizeRegion row for H06HDADXX130110:1:1101:10073:74407 (flag 83) on P12 cluster.

use gatk_haplotypecaller::assembly_region_finalize::{
    finalize_region_reads_for_assembly, gatk_min_tail_quality_for_assembly,
};
use gatk_haplotypecaller::assembly_region_iterator::AssemblyRegion;
use gatk_haplotypecaller::feature_context::FeatureContext;
use gatk_haplotypecaller::genome_loc::GenomePosition;
use gatk_haplotypecaller::read_unclip::{
    apply_hc_softclip_pre_step, hard_clip_adaptor_sequence, hard_clip_low_qual_ends,
    HcSoftclipPolicy,
};
use gatk_haplotypecaller::reference_context::ReferenceContext;
use rust_htslib::bam::{self, Read};
use std::path::PathBuf;

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn load_reverse_74407() -> Option<bam::Record> {
    let bam = repo_root().join("parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam");
    // Realworld BAM is gitignored / fetched locally — skip on bare CI checkouts.
    if !bam.is_file() {
        return None;
    }
    let mut r = bam::Reader::from_path(&bam).expect("bam");
    let mut rec = bam::Record::new();
    while let Some(res) = r.read(&mut rec) {
        res.expect("read record");
        if std::str::from_utf8(rec.qname()).unwrap().contains("74407") && rec.flags() & 16 != 0 {
            return Some(rec);
        }
    }
    panic!("reverse 74407 not found");
}

fn format_cigar(rec: &bam::Record) -> String {
    rec.cigar()
        .iter()
        .map(|c| {
            let op = match c {
                rust_htslib::bam::record::Cigar::Match(_) => 'M',
                rust_htslib::bam::record::Cigar::Ins(_) => 'I',
                rust_htslib::bam::record::Cigar::Del(_) => 'D',
                rust_htslib::bam::record::Cigar::SoftClip(_) => 'S',
                rust_htslib::bam::record::Cigar::HardClip(_) => 'H',
                rust_htslib::bam::record::Cigar::Equal(_) => '=',
                rust_htslib::bam::record::Cigar::Diff(_) => 'X',
                rust_htslib::bam::record::Cigar::RefSkip(_) => 'N',
                rust_htslib::bam::record::Cigar::Pad(_) => 'P',
            };
            format!(
                "{}{}",
                match c {
                    rust_htslib::bam::record::Cigar::Match(n)
                    | rust_htslib::bam::record::Cigar::Ins(n)
                    | rust_htslib::bam::record::Cigar::Del(n)
                    | rust_htslib::bam::record::Cigar::SoftClip(n)
                    | rust_htslib::bam::record::Cigar::HardClip(n)
                    | rust_htslib::bam::record::Cigar::Equal(n)
                    | rust_htslib::bam::record::Cigar::Diff(n)
                    | rust_htslib::bam::record::Cigar::RefSkip(n)
                    | rust_htslib::bam::record::Cigar::Pad(n) => *n,
                },
                op
            )
        })
        .collect()
}

#[test]
fn p12_74407_reverse_finalize_matches_java_cigar() {
    let Some(original) = load_reverse_74407() else {
        return;
    };
    let policy = HcSoftclipPolicy::haplotype_caller_defaults();
    let min_tail = gatk_min_tail_quality_for_assembly(10);

    let (reverted, _, _) = apply_hc_softclip_pre_step(&original, &policy);
    assert_eq!(format_cigar(&reverted), "250M");
    assert_eq!(reverted.pos(), 92307240);

    let low_qual = hard_clip_low_qual_ends(&reverted, min_tail);
    let adaptor = hard_clip_adaptor_sequence(&low_qual);
    assert_eq!(format_cigar(&adaptor), "31H219M");
    assert_eq!(adaptor.pos(), 92307271);
    assert_eq!(adaptor.seq().as_bytes().len(), 219);

    let region = AssemblyRegion {
        contig: "2".into(),
        start: GenomePosition::new_1based(92307228),
        end: GenomePosition::new_1based(92307400),
        extended_start: GenomePosition::new_1based(92307128),
        extended_end: GenomePosition::new_1based(92307500),
        extension: 100,
        reads: vec![gatk_haplotypecaller::share_record(original.clone())],
        read_qnames: vec![],
        reference: ReferenceContext::empty(),
        features: FeatureContext::empty(),
        pileup_loci: vec![],
        is_active: true,
    };
    let finalized = finalize_region_reads_for_assembly(
        &[gatk_haplotypecaller::share_record(original)],
        &region,
        true,
        min_tail,
        false,
    );
    let rev = finalized
        .iter()
        .find(|r| r.flags() & 16 != 0)
        .expect("reverse finalized");
    assert_eq!(rev.pos() + 1, 92307272, "1-based start");
    assert_eq!(format_cigar(rev), "31H219M");
}
