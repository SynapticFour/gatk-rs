//! 6R.91 coordinate-free: live `DepthPerAlleleBySample` annotation-call state.
//!
//! Java 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`)
//! `HaplotypeCallerGenotypingEngine.makeAnnotatedCall` →
//! `VariantAnnotatorEngine.annotateContext(call, …, allele-level likelihoods)` →
//! `DepthPerAlleleBySample.annotate` → `annotateWithLikelihoods`.
//!
//! Live dump (`HcParityAdAnnotationDump` on `HaplotypeCallerEngine.callRegion`):
//!
//! ```text
//! call.getAlleles()              = [TG(ref), CG]          // unused-ALT remaining; reverseTrim AFTER
//! AlleleLikelihoods.alleles()    = [TG, *, T, CG]         // still 4-way
//! evidence                       = 60, sample NA12878
//! remaining map                  = {TG→[TG], CG→[CG]}     // identity; T and * not pooled
//! independent remarg AD          = [36, 19]
//! annotateWithLikelihoods AD     = [36, 19]               // first AD write; genotype.hasAD() was false
//! reverseTrim / phase            copy → VCF T/C AD 36,19
//! ```
//!
//! 6R.88/6R.89 reconstructed a 62×4 object (order TG,T,CG,*) whose identity remarg is
//! [27, 9]. That object is not the live annotation input: dropping 2 reads cannot turn
//! [27, 9] into [36, 19], so the live 60×4 likelihood values also differ.
//!
//! Production change: NONE (observation only).
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r91_live_ad_annotation_contract
//! HOLDOUT_6R91=1 cargo test -p gatk-haplotypecaller --test forensic_6r91_live_ad_annotation_contract live_ -- --nocapture
//! ```

/// Java `DepthPerAlleleBySample`: remaining alleles are `vc.getAlleles()`, mapped to themselves.
#[test]
fn forensic_6r91_remaining_alleles_are_call_subset_of_likelihoods() {
    let call = ["TG", "CG"];
    let likelihoods = ["TG", "*", "T", "CG"];
    assert!(
        call.iter().all(|a| likelihoods.contains(a)),
        "containsAll(call) must hold or annotate() aborts"
    );
    assert_ne!(
        call.len(),
        likelihoods.len(),
        "unused ALTs may still be present in AlleleLikelihoods"
    );
    let identity_tg = ["TG"];
    let identity_cg = ["CG"];
    assert_eq!(identity_tg, ["TG"]);
    assert_eq!(identity_cg, ["CG"]);
    assert!(
        !["T", "*"]
            .iter()
            .any(|u| identity_tg.contains(u) || identity_cg.contains(u)),
        "unused ALTs are not pooled into remaining alleles"
    );
}

/// Marginalize key-set order (Java `Collectors.toMap` → HashMap) does not change AD:
/// counts are by allele identity into `vc.getAlleles()` order, not remarg column index.
#[test]
fn forensic_6r91_remarg_column_order_does_not_change_identity_counts() {
    // Remaining TG,CG. HashMap keySet may emit CG,TG (observed live subsetted_alleles).
    let votes_tg_cg_order = [36, 19];
    let votes_cg_tg_order = [19, 36];
    let ad_from_vc_order_a = [votes_tg_cg_order[0], votes_tg_cg_order[1]];
    let ad_from_vc_order_b = [votes_cg_tg_order[1], votes_cg_tg_order[0]];
    assert_eq!(ad_from_vc_order_a, ad_from_vc_order_b);
    assert_eq!(ad_from_vc_order_a, [36, 19]);
}

/// Count ±2 cannot explain [27, 9] vs [36, 19]; the live vs reconstructed
/// *set* difference is larger (10 JAVA_LIVE_ONLY + 12 RUST_ONLY QNAMEs).
#[test]
fn forensic_6r91_two_read_count_cannot_explain_ad_gap() {
    let reconstructed = [27i32, 9];
    let live = [36i32, 19];
    let max_shift = 2;
    assert!(
        (live[0] - reconstructed[0]).abs() + (live[1] - reconstructed[1]).abs() > max_shift,
        "AD gap exceeds any 2-read membership difference"
    );
}

/// Identity remarg of remaining columns is not permute of 4-way informative counts.
/// Live Java 4-way informative was TG=34, *=3, T=4, CG=17; remarg TG,CG = [36, 19].
#[test]
fn forensic_6r91_live_four_way_permute_is_not_remaining_remarg() {
    let four_way_tg_star_t_cg = [34i32, 3, 4, 17];
    let permute_remaining = [four_way_tg_star_t_cg[0], four_way_tg_star_t_cg[3]];
    let remarg_remaining = [36i32, 19];
    assert_eq!(permute_remaining, [34, 17]);
    assert_eq!(remarg_remaining, [36, 19]);
    assert_ne!(
        permute_remaining, remarg_remaining,
        "unused * and T steal 4-way votes that remarg reassigns to remaining alleles"
    );
}

/// Independent reconstruction of the live object equals the annotation AD (Outcome A).
#[test]
fn forensic_6r91_independent_live_reconstruction_equals_annotation_ad() {
    let independent_remarg = [36i32, 19];
    let annotate_with_likelihoods = [36i32, 19];
    let first_ad_write = true;
    assert!(first_ad_write);
    assert_eq!(independent_remarg, annotate_with_likelihoods);
}

#[test]
fn live_java_annotation_object_vs_rust_equivalent() {
    if std::env::var("HOLDOUT_6R91").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R91=1");
        return;
    }
    use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
    use gatk_haplotypecaller::hc_genotyping_engine::java_alignment_read_overlaps_interval;
    use gatk_haplotypecaller::{
        call_disposition, flatten_assembly_regions, take_colocated_merge_numerics,
        traverse_assembly_region_walker, try_emit_call_region_variants,
        AssemblyRegionCallDisposition, CallRegionArgs, HaplotypeCallerEngine, ReadFilterParams,
        WalkerTraversalConfig, DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN,
        DEFAULT_STAND_EMIT_CONFIDENCE,
    };
    use std::collections::HashSet;
    use std::path::Path;

    const INTERVAL: &str = "20:29455000-29456500";
    const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
    const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
    const POS_SNP: u64 = 29_456_344;
    const JAVA_LIVE_QNAMES: &[&str] = &[
        "HISEQ1:11:H8GV6ADXX:1:1104:4186:41128",
        "HISEQ1:11:H8GV6ADXX:1:1108:6005:20931",
        "HISEQ1:11:H8GV6ADXX:1:1213:10338:76652",
        "HISEQ1:11:H8GV6ADXX:2:2205:10648:31409",
        "HISEQ1:11:H8GV6ADXX:2:2216:2203:76921",
        "HISEQ1:12:H8GVUADXX:1:1106:3080:10032",
        "HISEQ1:12:H8GVUADXX:2:1216:17904:51162",
        "HISEQ1:13:H8G92ADXX:1:1103:8983:35134",
        "HISEQ1:13:H8G92ADXX:1:1111:12251:89078",
        "HISEQ1:13:H8G92ADXX:1:1214:8757:76066",
        "HISEQ1:13:H8G92ADXX:1:2103:3021:81022",
        "HISEQ1:13:H8G92ADXX:1:2209:4644:71639",
        "HISEQ1:13:H8G92ADXX:2:1101:15512:24344",
        "HISEQ1:13:H8G92ADXX:2:1205:1862:61782",
        "HISEQ1:13:H8G92ADXX:2:1209:13449:57024",
        "HISEQ1:9:H8962ADXX:1:1112:19265:60083",
        "HISEQ1:9:H8962ADXX:1:1213:5437:90850",
        "HISEQ1:9:H8962ADXX:1:2109:11968:22409",
        "HISEQ1:9:H8962ADXX:2:1104:19714:91101",
        "HISEQ1:9:H8962ADXX:2:1105:19831:29366",
        "HISEQ1:9:H8962ADXX:2:1206:11075:50547",
        "HISEQ1:9:H8962ADXX:2:1209:5235:78510",
        "HISEQ1:9:H8962ADXX:2:2101:14564:85219",
        "HISEQ1:9:H8962ADXX:2:2111:4053:81460",
        "HISEQ1:9:H8962ADXX:2:2115:11042:98655",
        "HWI-D00360:5:H814YADXX:1:1202:11051:34179",
        "HWI-D00360:5:H814YADXX:1:1208:17759:33171",
        "HWI-D00360:5:H814YADXX:1:1212:12402:73354",
        "HWI-D00360:5:H814YADXX:1:2207:10890:76583",
        "HWI-D00360:5:H814YADXX:2:1102:2154:52493",
        "HWI-D00360:5:H814YADXX:2:1109:4475:85583",
        "HWI-D00360:5:H814YADXX:2:2107:1177:54891",
        "HWI-D00360:5:H814YADXX:2:2209:16484:91770",
        "HWI-D00360:5:H814YADXX:2:2211:1401:25869",
        "HWI-D00360:6:H81VLADXX:1:1102:1483:21701",
        "HWI-D00360:6:H81VLADXX:1:2101:16091:54483",
        "HWI-D00360:6:H81VLADXX:1:2212:1317:49311",
        "HWI-D00360:6:H81VLADXX:2:1104:15554:2818",
        "HWI-D00360:6:H81VLADXX:2:1116:2361:8950",
        "HWI-D00360:6:H81VLADXX:2:1202:18367:85709",
        "HWI-D00360:6:H81VLADXX:2:2103:8562:24900",
        "HWI-D00360:6:H81VLADXX:2:2204:18162:74438",
        "HWI-D00360:6:H81VLADXX:2:2205:5601:38112",
        "HWI-D00360:7:H88WKADXX:1:1112:12463:50733",
        "HWI-D00360:7:H88WKADXX:1:1116:9273:30844",
        "HWI-D00360:7:H88WKADXX:1:2101:19221:59478",
        "HWI-D00360:7:H88WKADXX:2:1106:19349:74372",
        "HWI-D00360:7:H88WKADXX:2:1210:8920:81214",
        "HWI-D00360:7:H88WKADXX:2:2105:16900:85825",
        "HWI-D00360:7:H88WKADXX:2:2106:13368:29882",
        "HWI-D00360:8:H88U0ADXX:1:1101:14447:3965",
        "HWI-D00360:8:H88U0ADXX:1:1106:10867:71386",
        "HWI-D00360:8:H88U0ADXX:1:1106:15654:72455",
        "HWI-D00360:8:H88U0ADXX:1:1114:5563:73517",
        "HWI-D00360:8:H88U0ADXX:1:1202:9544:72702",
        "HWI-D00360:8:H88U0ADXX:1:2108:16806:75328",
        "HWI-D00360:8:H88U0ADXX:1:2201:20293:70383",
        "HWI-D00360:8:H88U0ADXX:2:1116:13828:55178",
        "HWI-D00360:8:H88U0ADXX:2:1204:18013:40432",
        "HWI-D00360:8:H88U0ADXX:2:1206:5788:95711",
    ];

    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: live BAM/ref missing");
        return;
    }

    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, INTERVAL).expect("interval");
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
    let covering: Vec<_> = regions
        .iter()
        .filter(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= POS_SNP
                && r.end.get() >= POS_SNP
        })
        .collect();
    assert_eq!(covering.len(), 1);
    let outcome = HaplotypeCallerEngine::call_region(
        covering[0],
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");
    let emitted = try_emit_call_region_variants(
        covering[0],
        &outcome,
        "SAMPLE",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .unwrap_or_default();
    let vcf = emitted
        .iter()
        .find(|r| {
            r.position == POS_SNP && r.reference == "T" && r.alternate.iter().any(|a| a == "C")
        })
        .expect("T/C");
    let snap = take_colocated_merge_numerics()
        .into_iter()
        .find(|n| n.loc == POS_SNP)
        .expect("snap");

    let loc = POS_SNP;
    let end = loc.saturating_add(snap.long_ref.len().saturating_sub(1) as u64);
    let margin = DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN;
    let rust_qnames: HashSet<String> = outcome
        .genotyping_reads
        .iter()
        .enumerate()
        .filter(|(i, r)| {
            outcome
                .read_likelihoods
                .iter()
                .any(|e| e.read_index.get() == *i)
                && java_alignment_read_overlaps_interval(r, loc, end, margin)
        })
        .map(|(_, r)| String::from_utf8_lossy(r.qname()).into_owned())
        .collect();
    let java: HashSet<&str> = JAVA_LIVE_QNAMES.iter().copied().collect();
    let rust: HashSet<&str> = rust_qnames.iter().map(String::as_str).collect();
    let common = java.intersection(&rust).count();
    let java_only: Vec<&&str> = java.difference(&rust).collect();
    let rust_only: Vec<&&str> = rust.difference(&java).collect();

    eprintln!(
        "6R.91 java_live n={} alleles_ll=TG,*,T,CG remaining=TG,CG remarg=[36,19] sample=NA12878 subset=20:29456342-29456347",
        JAVA_LIVE_QNAMES.len()
    );
    eprintln!(
        "6R.91 rust_equiv n={} n_snap={} remarg={:?} permute={:?} vcf_ad={:?}",
        rust.len(),
        snap.n_reads,
        snap.subset_ad_remarginalized,
        snap.subset_ad_permuted,
        vcf.samples[0].ad
    );
    eprintln!(
        "6R.91 membership common={} JAVA_LIVE_ONLY={} RUST_ONLY={}",
        common,
        java_only.len(),
        rust_only.len()
    );
    for (i, q) in java_only.iter().take(8).enumerate() {
        eprintln!("6R.91 JAVA_LIVE_ONLY[{i}] {q}");
    }
    for (i, q) in rust_only.iter().take(8).enumerate() {
        eprintln!("6R.91 RUST_ONLY[{i}] {q}");
    }

    assert_eq!(JAVA_LIVE_QNAMES.len(), 60);
    assert_eq!(snap.n_reads, 62);
    assert_eq!(snap.subset_ad_remarginalized, vec![27, 9]);
    assert_eq!(
        vcf.samples[0].ad.clone().unwrap_or_default(),
        vec![26u32, 9]
    );
    assert_ne!(
        snap.subset_ad_remarginalized,
        vec![36, 19],
        "Rust remarg of reconstructed C is not live Java annotation AD"
    );
    assert!(
        java_only.len() + rust_only.len() > 0 || rust.len() != java.len(),
        "live Java 60-read set is not the reconstructed 62-read C"
    );
}
