//! 6R.305: localize why one read is absent from Java's genotyping matrix.
//! The original measurement found Rust reinserting it as a 0,0,0 row.
//! 6R.311 production omits that row. The BAM identity below is unchanged.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r305_extra_read_membership -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::read_assembly_filter::{passes_assembly_read, AssemblyReadFilterConfig};
use gatk_haplotypecaller::{
    begin_likelihood_pipeline_observe, call_disposition, flatten_assembly_regions,
    set_forensic_6r294_padded_span, take_colocated_merge_numerics, take_likelihood_pipeline_cells,
    traverse_assembly_region_walker, AssemblyRegionCallDisposition, CallRegionArgs,
    HaplotypeCallerEngine, LikelihoodPipelineCell, ReadFilterParams, WalkerTraversalConfig,
};
use rust_htslib::bam::{self, Read};
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const JAVA_TSV: &str = "gatk-haplotypecaller/tests/6r297_java_hap_ll.tsv";
const LOCUS: u64 = 29_455_649;
const ALT: &str = "TGTTTG";
const JAVA_SPAN: (u64, u64) = (29_455_560, 29_455_728);
const QNAME: &str = "HISEQ1:11:H8GV6ADXX:2:1103:14252:55237";
const FLAGS: u16 = 97;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R305\t{key}\t{}", value.as_ref());
}

fn clear_other_forensics() {
    gatk_haplotypecaller::hc_genotyping_engine::set_forensic_6r287_global_hap_floor(None, None);
    gatk_haplotypecaller::hc_genotyping_engine::set_forensic_6r288_haplotype_substitute(
        None, None, 0, None,
    );
    gatk_haplotypecaller::likelihood_engine::set_forensic_6r289_target(None, 0);
    gatk_haplotypecaller::likelihood_engine::set_forensic_6r289_substitute(None, None);
}

fn java_excludes_this_read(text: &str) -> bool {
    !text.lines().any(|line| {
        let p: Vec<_> = line.split('\t').collect();
        p.len() >= 4 && p[0] == "6R297" && p[1] == "allele_ll" && p[2] == QNAME && p[3] == "97"
    })
}

fn target_cells<'a>(cells: &'a [LikelihoodPipelineCell]) -> Vec<&'a LikelihoodPipelineCell> {
    cells
        .iter()
        .filter(|c| c.stage == "post_kernel" && c.qname == QNAME && c.flags == FLAGS)
        .collect()
}

#[test]
fn forensic_6r305_extra_read_membership() {
    let root = repo_root();
    let bam_path = root.join(BAM_REL);
    let mut bam = bam::Reader::from_path(&bam_path).expect("bam");
    let header = bam::Header::from_template(bam.header());
    let text = header.to_hashmap();
    let sq = text.get("SQ").expect("SQ");
    let mut rec = bam::Record::new();
    let mut found = None;
    while let Some(r) = bam.read(&mut rec) {
        r.expect("read");
        if std::str::from_utf8(rec.qname()).ok() == Some(QNAME) && rec.flags() == FLAGS {
            found = Some(rec.clone());
            break;
        }
    }
    let rec = found.expect("frozen BAM record");
    let tid_name = |tid: i32| -> String {
        if tid < 0 {
            return "*".into();
        }
        sq.get(tid as usize)
            .and_then(|row| row.get("SN"))
            .map(|s| s.to_string())
            .unwrap_or_else(|| format!("tid={tid}"))
    };
    let start_1 = rec.pos() + 1;
    let end_1 = rec.cigar().end_pos();
    let qual = rec.qual().to_vec();
    let mate_same = rec.tid() == rec.mtid();
    let mate_predicate =
        !rec.is_paired() || rec.is_mate_unmapped() || rec.is_unmapped() || mate_same;
    let passes_filter = passes_assembly_read(&rec, &AssemblyReadFilterConfig::gatk_defaults());

    let java_text = std::fs::read_to_string(root.join(JAVA_TSV)).unwrap();
    assert!(java_excludes_this_read(&java_text));

    clear_other_forensics();
    set_forensic_6r294_padded_span(Some(JAVA_SPAN.0), Some(JAVA_SPAN.1));
    begin_likelihood_pipeline_observe();
    let ref_fasta = root.join(REF_REL);
    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, INTERVAL).expect("interval");
    let walk = traverse_assembly_region_walker(
        &dict,
        &specs,
        &ref_fasta,
        &bam_path,
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
            ) && r.start.get() <= LOCUS
                && r.end.get() >= LOCUS
        })
        .expect("ActiveFull");
    let in_region = region
        .reads
        .iter()
        .any(|r| std::str::from_utf8(r.qname()).ok() == Some(QNAME) && r.flags() == FLAGS);
    let outcome = HaplotypeCallerEngine::call_region(
        region,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");
    set_forensic_6r294_padded_span(None, None);
    let _ = outcome;
    let cells = take_likelihood_pipeline_cells();
    let numerics = take_colocated_merge_numerics();
    let site = numerics
        .iter()
        .find(|n| n.loc == LOCUS && n.alts.iter().any(|a| a == ALT))
        .expect("site");
    let mut row = None;
    for i in 0..site.ad_row_qname.len() {
        if site.ad_row_qname[i] == QNAME && site.ad_row_flags[i] == FLAGS {
            row = Some(site.ad_row_lls[i].clone());
        }
    }
    // 6R.305 originally observed `row == [0, 0, 0]` and a kernel of zeros.
    let kernel = target_cells(&cells);

    kv("qname", QNAME);
    kv("flags", FLAGS.to_string());
    kv("contig", tid_name(rec.tid()));
    kv("start", start_1.to_string());
    kv("end", end_1.to_string());
    kv("cigar", rec.cigar().to_string());
    kv("mapq", rec.mapq().to_string());
    kv("base_qual_len", qual.len().to_string());
    kv(
        "base_quals",
        qual.iter()
            .map(|q| q.to_string())
            .collect::<Vec<_>>()
            .join(","),
    );
    kv("mate_contig", tid_name(rec.mtid()));
    kv("mate_start", rec.mpos().saturating_add(1).to_string());
    kv("insert_size", rec.insert_size().to_string());
    kv("paired", rec.is_paired().to_string());
    kv("proper_pair", rec.is_proper_pair().to_string());
    kv("unmapped", rec.is_unmapped().to_string());
    kv("mate_unmapped", rec.is_mate_unmapped().to_string());
    kv("reverse", rec.is_reverse().to_string());
    kv("mate_reverse", rec.is_mate_reverse().to_string());
    kv("secondary", rec.is_secondary().to_string());
    kv("supplementary", rec.is_supplementary().to_string());
    kv("duplicate", rec.is_duplicate().to_string());
    kv("vendor_fail", rec.is_quality_check_failed().to_string());
    kv("read1", rec.is_first_in_template().to_string());
    kv("region_reads", region.reads.len().to_string());
    kv("in_assembly_region", in_region.to_string());
    kv("passes_assembly_read", passes_filter.to_string());
    kv("mate_predicate", mate_predicate.to_string());
    kv("java_allele_ll_present", "false");
    kv(
        "java_filter",
        "HaplotypeCallerEngine.filterNonPassingReads:948 MATE_ON_SAME_CONTIG_OR_NO_MAPPED_MATE.test == false; removeAll at 956. Active-region trace 450 -> 449, this read removed.",
    );
    kv("genotyping_reads", site.n_reads.to_string());
    kv(
        "genotyping_alleles",
        if row.is_none() {
            "ABSENT".to_string()
        } else {
            row.as_ref()
                .unwrap()
                .iter()
                .map(|v| format!("{v:.17}"))
                .collect::<Vec<_>>()
                .join(",")
        },
    );
    kv("post_kernel_cells", kernel.len().to_string());
    kv(
        "historical_rust_reinsertion",
        "6R.305 observed score_pairhmm_from_records_java_mate_contig inserting log10_likelihood 0.0",
    );
    kv(
        "historical_classification",
        "READ_MEMBERSHIP_JAVA_EXCLUDES_RUST_RETAINS",
    );
    kv("production_membership", "ABSENT");
    kv("production_change", "6R.311");

    assert_eq!(tid_name(rec.tid()), "20");
    assert_eq!(start_1, 29_455_546);
    assert_eq!(end_1, 29_455_693);
    assert_eq!(rec.cigar().to_string(), "148M");
    assert_eq!(rec.mapq(), 20);
    assert_eq!(tid_name(rec.mtid()), "7");
    assert_eq!(rec.mpos() + 1, 57_771_705);
    assert_eq!(rec.insert_size(), 0);
    assert!(rec.is_paired());
    assert!(!rec.is_proper_pair());
    assert!(!rec.is_unmapped());
    assert!(!rec.is_mate_unmapped());
    assert!(!rec.is_reverse());
    assert!(rec.is_mate_reverse());
    assert!(!rec.is_secondary());
    assert!(!rec.is_supplementary());
    assert!(!rec.is_duplicate());
    assert!(!rec.is_quality_check_failed());
    assert!(rec.is_first_in_template());
    assert!(!mate_predicate);
    assert!(!passes_filter);
    assert!(in_region);
    assert!(row.is_none(), "failed-mate genotyping row");
    assert!(kernel.is_empty(), "failed-mate PairHMM cells");
    assert!(start_1 <= JAVA_SPAN.1 as i64 && end_1 >= JAVA_SPAN.0 as i64);
}
