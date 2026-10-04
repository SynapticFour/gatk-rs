//! 7R.26 production closure at `20:29455951 T>C`.
//! Policy env stays unset. No injected CIGAR, likelihood, or PL.

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::hc_genotyping_engine::java_alignment_read_overlaps_interval;
use gatk_haplotypecaller::read_unclip::alignment_end_1based;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    AssemblyRegionCallDisposition, CallRegionArgs, HaplotypeCallerEngine, ReadFilterParams,
    WalkerTraversalConfig,
};
use rust_htslib::bam::Record;
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

const QNAME: &str = "HISEQ1:9:H8962ADXX:1:1107:11064:51936";
const FLAG: u16 = 147;
const SITE: u64 = 29_455_951;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn cigar_of(rec: &Record) -> String {
    rec.cigar().to_string()
}

#[test]
fn production_realignment_excludes_diagnostic_read_and_matches_java_gl() {
    assert!(std::env::var_os("GATK_RS_EXPERIMENTAL_KBEST_POLICY").is_none());
    let root = repo_root();
    let ref_fasta = root.join("parity/realworld/assets/hs37d5.simple.fa");
    let bam = root.join("parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam");
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: canonical fixtures absent");
        return;
    }
    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).unwrap();
    let specs = parse_intervals_cli_string(&dict, "20:29455000-29456500").unwrap();
    let walk = traverse_assembly_region_walker(
        &dict,
        &specs,
        &ref_fasta,
        &bam,
        &ReadFilterParams::gatk_standard_hc(),
        &WalkerTraversalConfig::gatk_haplotype_caller_production(100),
    )
    .unwrap();
    let region = flatten_assembly_regions(&walk)
        .into_iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= SITE
                && r.end.get() >= SITE
        })
        .expect("region covering 29455951");
    let outcome = HaplotypeCallerEngine::call_region(
        &region,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .unwrap()
    .unwrap();

    let read = outcome
        .genotyping_reads
        .iter()
        .find(|r| {
            let rec = r.as_ref();
            String::from_utf8_lossy(rec.qname()) == QNAME && rec.flags() == FLAG
        })
        .expect("diagnostic read")
        .as_ref();
    let start = read.pos() + 1;
    let end = i64::from(alignment_end_1based(read));
    let cigar = cigar_of(read);
    eprintln!("7R26\tdiagnostic\tstart={start}\tend={end}\tcigar={cigar}");
    assert_eq!(cigar, "5M2D136M2I5M");
    assert_eq!(start, 29_455_800);
    assert_eq!(end, 29_455_947);

    let matrix_reads: BTreeSet<usize> = outcome
        .read_likelihoods
        .iter()
        .map(|cell| cell.read_index.get())
        .collect();
    let overlap = outcome
        .genotyping_reads
        .iter()
        .enumerate()
        .filter(|(index, rec)| {
            matrix_reads.contains(index)
                && java_alignment_read_overlaps_interval(rec.as_ref(), SITE, SITE, 2)
        })
        .count();
    let read_overlaps = java_alignment_read_overlaps_interval(read, SITE, SITE, 2);
    eprintln!(
        "7R26\toverlap_count\t{overlap}\tmatrix_reads\t{}\tdiagnostic_overlaps\t{read_overlaps}",
        matrix_reads.len()
    );
    assert!(!read_overlaps);
    assert_eq!(overlap, 111);

    let call = outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based.get() == SITE
                && c.event.ref_allele == "T"
                && c.event.alt_allele == "C"
        })
        .expect("20:29455951 T/C");
    let gl = &call.genotype.genotype_log10_likelihoods;
    let pl = call.genotype.format.pl_as_i32();
    eprintln!(
        "7R26\tgl\t{}\tbits={:#x},{:#x},{:#x}",
        gl.iter()
            .map(|v| format!("{v:.17}"))
            .collect::<Vec<_>>()
            .join(","),
        gl.first().copied().unwrap_or(0.0).to_bits(),
        gl.get(1).copied().unwrap_or(0.0).to_bits(),
        gl.get(2).copied().unwrap_or(0.0).to_bits()
    );
    eprintln!("7R26\tpl\t{pl:?}");
    eprintln!(
        "7R26\tannotation_likelihoods\t{}",
        call.annotation_likelihoods.len()
    );
    assert_eq!(gl.len(), 3);
    assert_eq!(gl[0].to_bits(), 0xc0836ff5bb98d511);
    assert_eq!(gl[1].to_bits(), 0xc077ccd6ec7a68b8);
    assert_eq!(gl[2].to_bits(), 0xc0818577a4d086f0);
    assert_eq!(pl, vec![2412, 0, 1799]);
}
