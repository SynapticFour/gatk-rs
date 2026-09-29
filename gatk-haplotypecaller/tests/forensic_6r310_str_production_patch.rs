//! 6R.310: production STR trim padding uses Java's anchored tandem-repeat count.
//! Failed-mate exclusion and the Java read-coordinate sort stay off.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r310_str_production_patch -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary};
use gatk_haplotypecaller::assembly_region_finalize::set_forensic_6r308_java_read_coordinate_order;
use gatk_haplotypecaller::{
    begin_hap_list_observe, begin_trim_calc_observe, call_disposition, flatten_assembly_regions,
    set_forensic_6r294_padded_span, set_forensic_6r306_exclude_failed_mate,
    take_colocated_merge_numerics, take_hap_list_snaps, take_hap_list_trim_span,
    take_trim_calc_snap, tandem_repeat_at_event, traverse_assembly_region_walker,
    AssemblyRegionCallDisposition, CallRegionArgs, HapListSnap, HaplotypeCallerEngine,
    ReadFilterParams, ReferenceContext, WalkerTraversalConfig,
};
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const LOCUS: u64 = 29_455_649;
const ALT: &str = "TGTTTG";
const EVENT: u64 = 29_455_644;
const WINDOW_START: u64 = 29_455_460;
const WINDOW_END: u64 = 29_455_844;
const JAVA_SPAN: (u64, u64) = (29_455_560, 29_455_728);
const DIAG_QNAME: &str = "HISEQ1:11:H8GV6ADXX:1:2116:18670:99941";
const MATE_QNAME: &str = "HISEQ1:11:H8GV6ADXX:2:1103:14252:55237";

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R310\t{key}\t{}", value.as_ref());
}

fn pairhmm_snap(snaps: &[HapListSnap]) -> &HapListSnap {
    snaps
        .iter()
        .find(|s| s.stage == "pairhmm_input")
        .expect("pairhmm_input haplotypes")
}

#[test]
fn generic_anchor_is_excluded() {
    let rep = tandem_repeat_at_event(100, b"CATTTTT", 101, b"A", b"AT").expect("repeat");
    assert_eq!(rep.unit, b"T");
    assert_eq!(rep.counts, vec![5, 6]);
    assert_ne!(rep.unit, b"A");
}

#[test]
fn generic_alternate_allele_extends_the_repeat() {
    let rep = tandem_repeat_at_event(100, b"CATTT", 101, b"A", b"ATTT").expect("repeat");
    assert_eq!(rep.counts, vec![3, 6]);
    assert_eq!(rep.padding_bases(), 6);
}

#[test]
fn generic_non_str_indel_has_no_repeat() {
    assert!(tandem_repeat_at_event(100, b"CATGC", 101, b"A", b"AG").is_none());
}

#[test]
fn generic_snp_is_not_a_tandem_repeat() {
    assert!(tandem_repeat_at_event(100, b"CATTTTTTTTT", 101, b"A", b"G").is_none());
}

#[test]
fn forensic_6r310_str_production_patch() {
    set_forensic_6r294_padded_span(None, None);
    set_forensic_6r306_exclude_failed_mate(false);
    set_forensic_6r308_java_read_coordinate_order(false);
    gatk_haplotypecaller::hc_genotyping_engine::set_forensic_6r287_global_hap_floor(None, None);
    gatk_haplotypecaller::hc_genotyping_engine::set_forensic_6r288_haplotype_substitute(
        None, None, 0, None,
    );
    gatk_haplotypecaller::likelihood_engine::set_forensic_6r289_substitute(None, None);
    gatk_haplotypecaller::likelihood_engine::set_forensic_6r289_target(Some(DIAG_QNAME), 99);

    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let mut cache = ReferenceWindowCache::new(&ref_fasta, 4);
    let ctx = ReferenceContext::from_interval(&dict, &mut cache, "20", WINDOW_START, WINDOW_END)
        .expect("ref");
    let rep = tandem_repeat_at_event(ctx.window_start, ctx.bases.as_slice(), EVENT, b"A", b"AT")
        .expect("canonical repeat");
    assert_eq!(rep.unit, b"T");
    assert_eq!(rep.counts, vec![8, 9]);
    assert_eq!(rep.padding_bases(), 9);
    kv("rust_repeat_unit", "T");
    kv("rust_repeat_counts", format!("{:?}", rep.counts));
    kv("rust_repeat_count", "9");

    let bam_path = root.join(BAM_REL);
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

    begin_trim_calc_observe();
    begin_hap_list_observe();
    let outcome = HaplotypeCallerEngine::call_region(
        region,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");
    let _ = outcome;
    let trim = take_trim_calc_snap().expect("trim snap");
    let applied = take_hap_list_trim_span().expect("applied span");
    let snaps = take_hap_list_snaps();
    let read_snaps = gatk_haplotypecaller::likelihood_engine::take_forensic_6r289_snaps();
    gatk_haplotypecaller::likelihood_engine::set_forensic_6r289_target(None, 0);

    let event = trim
        .events
        .iter()
        .find(|e| e.start == EVENT && e.is_indel)
        .expect("A>AT event");
    assert_eq!(event.ref_al, "A");
    assert_eq!(event.alt_al, "AT");
    kv(
        "trim_variant_alleles_reached_event_list",
        format!("{}>{}", event.ref_al, event.alt_al),
    );

    assert_eq!((trim.padded_start, trim.padded_end), JAVA_SPAN);
    assert_eq!((applied.trim_start, applied.trim_end), JAVA_SPAN);
    kv(
        "production_padded_span",
        format!("{}-{}", trim.padded_start, trim.padded_end),
    );

    let read = read_snaps.first().expect("diagnostic read");
    let read_repr = format!(
        "{} @ {} len={}",
        read.cigar,
        read.pos_1based,
        read.read_bases.len()
    );
    assert_eq!(read_repr, "51H97M @ 29455560 len=97");
    kv("diagnostic_read", &read_repr);

    let hap = pairhmm_snap(&snaps);
    let ref_h = hap
        .columns
        .iter()
        .find(|c| c.is_reference)
        .expect("ref hap");
    let lens: BTreeSet<_> = hap.columns.iter().map(|c| c.len).collect();
    let tgtttg_lens: Vec<usize> = (19..24).map(|i| hap.columns[i].len).collect();
    assert_eq!(ref_h.len, 169);
    assert_eq!(ref_h.cigar, "169M");
    assert_eq!(lens, BTreeSet::from([169, 170, 172, 174]));
    assert_eq!(tgtttg_lens, vec![174, 174, 174, 174, 174]);
    assert_eq!(hap.n, 24);
    kv(
        "reference_haplotype",
        format!("{} {}", ref_h.cigar, ref_h.len),
    );
    kv(
        "haplotype_lengths",
        hap.columns
            .iter()
            .map(|c| c.len.to_string())
            .collect::<Vec<_>>()
            .join(","),
    );
    kv("tgtttg_indices", "19-23");

    let numerics = take_colocated_merge_numerics();
    let site = numerics
        .iter()
        .find(|n| n.loc == LOCUS && n.alts.iter().any(|a| a == ALT))
        .expect("merge numerics");
    let mate_present = site.ad_row_qname.iter().any(|q| q == MATE_QNAME);
    kv("genotyping_reads", site.n_reads.to_string());
    kv("pairhmm_reads", site.n_pairhmm_reads.to_string());
    kv("failed_mate_still_present", mate_present.to_string());
    kv("subset_pl", format!("{:?}", site.subset_pl));
    kv("merged_pl", format!("{:?}", site.merged_pl));
    kv(
        "classification",
        "STR_PRODUCTION_PATCH_REPRODUCES_JAVA_GEOMETRY",
    );
    let _ = site;
}
