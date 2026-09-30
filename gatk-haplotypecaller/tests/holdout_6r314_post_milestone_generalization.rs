//! 6R.314: independent holdouts for the four closed production semantics.
//!
//! Not the canonical witness `20:29455649 T/TGTTTG` and not `20:29455644 A>AT`.
//! Skipped unless `HOLDOUT_6R314=1`. Does not run `HOLDOUT_6R243`.
//!
//! ```text
//! HOLDOUT_6R314=1 cargo test -p gatk-haplotypecaller --test holdout_6r314_post_milestone_generalization -- --nocapture --test-threads=1
//! cargo test -p gatk-haplotypecaller --lib holdout_6r314_failed_mate_zero_row -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::SequenceDictionary;
use gatk_haplotypecaller::assembly_region_finalize::{
    clip_finalized_reads_in_place, java_read_coordinate_compare,
    set_forensic_6r308_java_read_coordinate_order,
};
use gatk_haplotypecaller::{
    approximate_log10_sum_log10_pair, biallelic_genotype_log10_likelihoods_gatk,
    failed_mate_evidence_keys, log10_sum_log10, passes_mate_on_same_contig_or_no_mapped_mate,
    tandem_repeat_at_event, AssemblyRegion, AssemblyRegionTrimmer, AssemblyRegionTrimmerConfig,
    GenomePosition, ReadLikelihoodRow, ReferenceContext, TrimVariant,
};
use rust_htslib::bam::record::{Cigar, CigarString};
use rust_htslib::bam::{HeaderView, Record};
use std::sync::Arc;

fn armed() -> bool {
    if std::env::var("HOLDOUT_6R314").ok().as_deref() != Some("1") {
        eprintln!("6R314 holdout skipped; set HOLDOUT_6R314=1");
        return false;
    }
    true
}

fn fail(class: &str, checkpoint: &str, detail: &str) -> ! {
    panic!("6R314\tclassification\t{class}\tcheckpoint\t{checkpoint}\tdetail\t{detail}");
}

fn region() -> AssemblyRegion {
    AssemblyRegion {
        contig: "chrHold".into(),
        start: GenomePosition::new_1based(3000),
        end: GenomePosition::new_1based(5000),
        is_active: true,
        extended_start: GenomePosition::new_1based(1),
        extended_end: GenomePosition::new_1based(20000),
        extension: 100,
        reads: Vec::new(),
        read_qnames: Vec::new(),
        reference: ReferenceContext::empty(),
        features: gatk_haplotypecaller::FeatureContext::empty(),
        pileup_loci: Vec::new(),
    }
}

fn trimmer() -> AssemblyRegionTrimmer {
    let mut dict = SequenceDictionary::new();
    dict.add_contig("chrHold".into(), 20000);
    AssemblyRegionTrimmer::new(
        AssemblyRegionTrimmerConfig::gatk_defaults(),
        &dict,
        "chrHold",
    )
}

fn ctx(start: u64, bases: &[u8]) -> ReferenceContext {
    ReferenceContext {
        contig: "chrHold".into(),
        start,
        end: start + bases.len() as u64 - 1,
        window_start: start,
        window_end: start + bases.len() as u64 - 1,
        bases: gatk_haplotypecaller::reference_context::SharedBases::from_slice(bases),
    }
}

struct StrCase {
    name: &'static str,
    bases: &'static [u8],
    event: u64,
    end: u64,
    indel: bool,
    refer: &'static [u8],
    alt: &'static [u8],
    unit: Option<&'static [u8]>,
    counts: Option<[usize; 2]>,
    pad: u32,
}

fn str_cases() -> Vec<StrCase> {
    vec![
        StrCase {
            name: "unit-len-1-insertion",
            bases: b"CGGGGG",
            event: 4000,
            end: 4000,
            indel: true,
            refer: b"C",
            alt: b"CG",
            unit: Some(b"G"),
            counts: Some([5, 6]),
            pad: 75 + 6,
        },
        StrCase {
            name: "unit-len-gt1-insertion",
            bases: b"GCAGCAGCAG",
            event: 4000,
            end: 4000,
            indel: true,
            refer: b"G",
            alt: b"GCAG",
            unit: Some(b"CAG"),
            counts: Some([3, 4]),
            pad: 75 + 12,
        },
        StrCase {
            name: "unit-len-1-deletion",
            bases: b"ATTTTT",
            event: 4000,
            end: 4001,
            indel: true,
            refer: b"AT",
            alt: b"A",
            unit: Some(b"T"),
            counts: Some([5, 4]),
            pad: 75 + 5,
        },
        StrCase {
            name: "unit-len-gt1-deletion",
            bases: b"GATATAT",
            event: 4000,
            end: 4002,
            indel: true,
            refer: b"GAT",
            alt: b"G",
            unit: Some(b"AT"),
            counts: Some([3, 2]),
            pad: 75 + 6,
        },
        StrCase {
            name: "alt-extends-repeat",
            bases: b"ATTTT",
            event: 4000,
            end: 4000,
            indel: true,
            refer: b"A",
            alt: b"ATT",
            unit: Some(b"T"),
            counts: Some([4, 6]),
            pad: 75 + 6,
        },
        StrCase {
            name: "anchor-excluded",
            bases: b"XCCCCC",
            event: 4000,
            end: 4000,
            indel: true,
            refer: b"X",
            alt: b"XC",
            unit: Some(b"C"),
            counts: Some([5, 6]),
            pad: 75 + 6,
        },
        StrCase {
            name: "non-str-indel",
            bases: b"CATGC",
            event: 4001,
            end: 4001,
            indel: true,
            refer: b"A",
            alt: b"AG",
            unit: None,
            counts: None,
            pad: 75,
        },
        StrCase {
            name: "snp-equal-length",
            bases: b"ACGT",
            event: 4000,
            end: 4000,
            indel: false,
            refer: b"A",
            alt: b"G",
            unit: None,
            counts: None,
            pad: 20,
        },
        StrCase {
            name: "insertion-without-tandem",
            bases: b"CACGT",
            event: 4000,
            end: 4000,
            indel: true,
            refer: b"C",
            alt: b"CT",
            unit: None,
            counts: None,
            pad: 75,
        },
    ]
}

#[test]
fn holdout_6r314_str() {
    if !armed() {
        return;
    }
    let trimmer = trimmer();
    let region = region();
    for case in str_cases() {
        assert_ne!(case.event, 29_455_644, "canonical A>AT locus");
        assert_ne!(case.counts, Some([8, 9]), "canonical T run counts");
        let found =
            tandem_repeat_at_event(case.event, case.bases, case.event, case.refer, case.alt);
        match (case.unit, case.counts) {
            (Some(unit), Some(counts)) => {
                let found = found.unwrap_or_else(|| {
                    fail(
                        "STR_HOLDOUT_DIVERGENCE",
                        case.name,
                        "tandem_repeat_at_event returned None",
                    )
                });
                if found.unit.as_slice() != unit || found.counts.as_slice() != counts {
                    fail(
                        "STR_HOLDOUT_DIVERGENCE",
                        case.name,
                        &format!("unit {:?} counts {:?}", found.unit, found.counts),
                    );
                }
                let expect_add = counts.into_iter().max().unwrap() * unit.len();
                if found.padding_bases() != expect_add {
                    fail(
                        "STR_HOLDOUT_DIVERGENCE",
                        case.name,
                        &format!("padding_bases {}", found.padding_bases()),
                    );
                }
            }
            (None, None) => {
                if found.is_some() {
                    fail(
                        "STR_HOLDOUT_DIVERGENCE",
                        case.name,
                        "expected no tandem repeat",
                    );
                }
            }
            _ => fail("STR_HOLDOUT_DIVERGENCE", case.name, "case table"),
        }
        let reference = ctx(case.event, case.bases);
        let vars = vec![TrimVariant {
            contig: "chrHold".into(),
            start: case.event,
            end: case.end,
            is_indel: case.indel,
            ref_allele: String::from_utf8(case.refer.to_vec()).unwrap(),
            alt_allele: String::from_utf8(case.alt.to_vec()).unwrap(),
        }];
        let res = trimmer.trim(&region, &vars, Some(&reference));
        let start = case.event - u64::from(case.pad);
        let end = case.end + u64::from(case.pad);
        if res.padded_variant_start != Some(start) || res.padded_variant_end != Some(end) {
            fail(
                "STR_HOLDOUT_DIVERGENCE",
                case.name,
                &format!(
                    "span {:?}..={:?} expected {start}..={end}",
                    res.padded_variant_start, res.padded_variant_end
                ),
            );
        }
        println!("6R314\tcheckpoint\t{}\tPASS", case.name);
    }
    println!("6R314\tclassification\tSTR_HOLDOUT_PASS");
}

fn mate_rec(
    qname: &str,
    paired: bool,
    unmapped: bool,
    mate_unmapped: bool,
    mtid: i32,
    mapq: u8,
    len: usize,
) -> Record {
    let mut r = Record::new();
    r.set_header(Arc::new(HeaderView::from_bytes(
        b"@HD\tVN:1.0\n@SQ\tSN:hold-a\tLN:1000\n@SQ\tSN:hold-b\tLN:1000\n",
    )));
    let bases = vec![b'A'; len];
    let qual = vec![20u8; len];
    r.set(
        qname.as_bytes(),
        Some(&CigarString::from(vec![Cigar::Match(len as u32)])),
        &bases,
        &qual,
    );
    r.set_tid(0);
    r.set_pos(10);
    r.set_mtid(mtid);
    r.set_mapq(mapq);
    if paired {
        r.set_paired();
    }
    if unmapped {
        r.set_unmapped();
    } else {
        r.unset_unmapped();
    }
    if mate_unmapped {
        r.set_mate_unmapped();
    } else {
        r.unset_mate_unmapped();
    }
    r
}

#[test]
fn holdout_6r314_failed_mate() {
    if !armed() {
        return;
    }
    let unpaired = mate_rec("hold314-unpaired", false, false, false, 1, 20, 16);
    let unmapped_mate = mate_rec("hold314-unmapped-mate", true, false, true, 1, 20, 16);
    let unmapped_read = mate_rec("hold314-unmapped-read", true, true, false, 1, 20, 16);
    let same = mate_rec("hold314-same", true, false, false, 0, 20, 16);
    let cross = mate_rec("hold314-cross", true, false, false, 1, 20, 16);
    let low_mapq = mate_rec("hold314-low-mapq", true, false, false, 0, 1, 16);
    let short = mate_rec("hold314-short", true, false, false, 0, 20, 8);
    let cases = [
        ("unpaired", &unpaired, true),
        ("unmapped-mate", &unmapped_mate, true),
        ("unmapped-read", &unmapped_read, true),
        ("same-contig", &same, true),
        ("different-contig", &cross, false),
        ("low-mapq", &low_mapq, true),
        ("short-read", &short, true),
    ];
    for (name, rec, keep) in cases {
        if passes_mate_on_same_contig_or_no_mapped_mate(rec) != keep {
            fail("FAILED_MATE_HOLDOUT_DIVERGENCE", name, "predicate");
        }
        println!("6R314\tcheckpoint\tpredicate-{name}\tPASS");
    }
    assert_ne!(cross.qname(), b"HISEQ1:11:H8GV6ADXX:2:1103:14252:55237");
    let originals = [
        unpaired,
        unmapped_mate,
        unmapped_read,
        same,
        cross,
        low_mapq,
        short,
    ];
    let keys = failed_mate_evidence_keys(&originals);
    if keys.len() != 1 || !keys.contains(&(b"hold314-cross".to_vec(), originals[4].flags())) {
        fail(
            "FAILED_MATE_HOLDOUT_DIVERGENCE",
            "evidence-keys",
            &format!("keys {keys:?}"),
        );
    }
    let mut evidence = originals.to_vec();
    evidence.retain(|r| !keys.contains(&(r.qname().to_vec(), r.flags())));
    if evidence.iter().any(|r| r.qname() == b"hold314-cross") {
        fail(
            "FAILED_MATE_HOLDOUT_DIVERGENCE",
            "pairhmm-membership",
            "cross-contig read retained",
        );
    }
    for name in [
        "hold314-unpaired",
        "hold314-unmapped-mate",
        "hold314-unmapped-read",
        "hold314-same",
        "hold314-low-mapq",
        "hold314-short",
    ] {
        if !evidence.iter().any(|r| r.qname() == name.as_bytes()) {
            fail("FAILED_MATE_HOLDOUT_DIVERGENCE", "pairhmm-membership", name);
        }
    }
    println!("6R314\tcheckpoint\tmembership\tPASS");
    println!("6R314\tclassification\tFAILED_MATE_HOLDOUT_PASS");
}

fn order_rec(
    name: &str,
    flags: u16,
    tid: i32,
    pos: i64,
    mapq: u8,
    mtid: i32,
    mpos: i64,
    tlen: i64,
) -> Record {
    let mut r = Record::new();
    r.set_header(Arc::new(HeaderView::from_bytes(
        b"@HD\tVN:1.0\n@SQ\tSN:chrHold\tLN:20000\n@SQ\tSN:chrMate\tLN:20000\n",
    )));
    r.set(
        name.as_bytes(),
        Some(&CigarString::from(vec![Cigar::Match(24)])),
        b"ACGTACGTACGTACGTACGTACGT",
        &vec![25u8; 24],
    );
    r.set_tid(tid);
    r.set_pos(pos);
    r.set_flags(flags);
    r.set_mapq(mapq);
    r.set_mtid(mtid);
    r.set_mpos(mpos);
    r.set_insert_size(tlen);
    r
}

fn order_region() -> AssemblyRegion {
    AssemblyRegion {
        contig: "chrHold".into(),
        start: GenomePosition::new_1based(1500),
        end: GenomePosition::new_1based(2500),
        is_active: true,
        extended_start: GenomePosition::new_1based(1),
        extended_end: GenomePosition::new_1based(20000),
        extension: 50,
        reads: Vec::new(),
        read_qnames: Vec::new(),
        reference: ReferenceContext::empty(),
        features: gatk_haplotypecaller::FeatureContext::empty(),
        pileup_loci: Vec::new(),
    }
}

fn assert_order(first: &Record, second: &Record, field: &str) {
    let cmp = java_read_coordinate_compare(first, second);
    if cmp >= 0 {
        fail(
            "READ_ORDER_HOLDOUT_DIVERGENCE",
            field,
            &format!("compare {cmp}"),
        );
    }
    println!("6R314\tcheckpoint\tfield-{field}\tPASS");
}

#[test]
fn holdout_6r314_read_order() {
    if !armed() {
        return;
    }
    set_forensic_6r308_java_read_coordinate_order(true);
    let base = order_rec("read", 99, 0, 2000, 20, 0, 10, 100);
    let mut other_tid = base.clone();
    other_tid.set_tid(1);
    assert_order(&base, &other_tid, "contig");
    let mut later = base.clone();
    later.set_pos(2001);
    assert_order(&base, &later, "start");
    let forward = order_rec("zeta-forward", 99, 0, 2000, 20, 0, 10, 100);
    let reverse = order_rec("aaa-reverse", 83, 0, 2000, 60, 0, 10, 100);
    assert_ne!(forward.qname(), b"HISEQ1:11:H8GV6ADXX:1:1116:4033:65919");
    assert_order(&forward, &reverse, "forward-before-reverse");
    let name_a = order_rec("alpha", 99, 0, 2000, 20, 0, 10, 100);
    let name_b = order_rec("beta", 99, 0, 2000, 20, 0, 10, 100);
    assert_order(&name_a, &name_b, "name");
    let flags_lo = order_rec("alpha", 97, 0, 2000, 20, 0, 10, 100);
    let flags_hi = order_rec("alpha", 99, 0, 2000, 20, 0, 10, 100);
    assert_order(&flags_lo, &flags_hi, "flags");
    let mapq_lo = order_rec("alpha", 99, 0, 2000, 10, 0, 10, 100);
    let mapq_hi = order_rec("alpha", 99, 0, 2000, 30, 0, 10, 100);
    assert_order(&mapq_lo, &mapq_hi, "mapq");
    let mate0 = order_rec("alpha", 99, 0, 2000, 20, 0, 10, 100);
    let mate1 = order_rec("alpha", 99, 0, 2000, 20, 1, 10, 100);
    assert_order(&mate0, &mate1, "mate-ref");
    let mate_near = order_rec("alpha", 99, 0, 2000, 20, 0, 4, 100);
    let mate_far = order_rec("alpha", 99, 0, 2000, 20, 0, 40, 100);
    assert_order(&mate_near, &mate_far, "mate-start");
    let short_frag = order_rec("alpha", 99, 0, 2000, 20, 0, 10, 40);
    let long_frag = order_rec("alpha", 99, 0, 2000, 20, 0, 10, 90);
    assert_order(&short_frag, &long_frag, "fragment-length");

    let expected = [
        order_rec("alpha", 97, 0, 2000, 30, 0, 10, 100),
        order_rec("alpha", 99, 0, 2000, 10, 0, 10, 100),
        order_rec("alpha", 99, 0, 2000, 30, 0, 10, 100),
        order_rec("mate-pos", 99, 0, 2000, 20, 0, 4, 100),
        order_rec("mate-pos", 99, 0, 2000, 20, 0, 40, 100),
        order_rec("mate-ref", 99, 0, 2000, 20, 0, 10, 100),
        order_rec("mate-ref", 99, 0, 2000, 20, 1, 10, 100),
        order_rec("zeta-forward", 99, 0, 2000, 20, 0, 10, 100),
        order_rec("aaa-reverse", 83, 0, 2000, 60, 0, 10, 100),
    ];
    let mut by_comparator = expected.to_vec();
    by_comparator.sort_by(|a, b| java_read_coordinate_compare(a, b).cmp(&0));
    for (i, (hand, cmp)) in expected.iter().zip(by_comparator.iter()).enumerate() {
        if hand.qname() != cmp.qname() || hand.flags() != cmp.flags() || hand.mapq() != cmp.mapq() {
            fail(
                "READ_ORDER_HOLDOUT_DIVERGENCE",
                "hand-vs-comparator",
                &format!("index {i}"),
            );
        }
    }
    let mut production = expected.to_vec();
    production.reverse();
    clip_finalized_reads_in_place(&mut production, &order_region());
    if production.len() != expected.len() {
        fail(
            "READ_ORDER_HOLDOUT_DIVERGENCE",
            "clip-membership",
            &format!("len {}", production.len()),
        );
    }
    let mut positional = 0usize;
    for (i, (got, want)) in production.iter().zip(expected.iter()).enumerate() {
        if got.qname() != want.qname()
            || got.flags() != want.flags()
            || got.mapq() != want.mapq()
            || got.mtid() != want.mtid()
            || got.mpos() != want.mpos()
        {
            positional += 1;
            if positional == 1 {
                fail(
                    "READ_ORDER_HOLDOUT_DIVERGENCE",
                    "identity-permutation",
                    &format!(
                        "index {i} got {} flags {} want {} flags {}",
                        String::from_utf8_lossy(got.qname()),
                        got.flags(),
                        String::from_utf8_lossy(want.qname()),
                        want.flags()
                    ),
                );
            }
        }
    }
    println!("6R314\tcheckpoint\tidentity-permutation\tYES");
    println!("6R314\tcheckpoint\tpositional-differences\t0");

    set_forensic_6r308_java_read_coordinate_order(false);
    let mut old = expected.to_vec();
    clip_finalized_reads_in_place(&mut old, &order_region());
    set_forensic_6r308_java_read_coordinate_order(true);
    if old.first().map(|r| r.qname()) != Some(b"aaa-reverse".as_slice()) {
        fail(
            "READ_ORDER_HOLDOUT_DIVERGENCE",
            "name-sort-counterfactual",
            "pre-6R.312 sort did not put aaa-reverse first",
        );
    }
    if production.first().map(|r| r.qname()) == Some(b"aaa-reverse".as_slice()) {
        fail(
            "READ_ORDER_HOLDOUT_DIVERGENCE",
            "production-order",
            "reverse read sorted first",
        );
    }
    println!("6R314\tclassification\tREAD_ORDER_HOLDOUT_PASS");
}

/// Independent GATK 4.4 `approximateLog10SumLog10` + `JacobianLogTable`.
/// Out-of-range table slots return `0.0`, matching Java `get`.
fn java_fast_round(d: f64) -> usize {
    if d > 0.0 {
        (d + 0.5) as usize
    } else {
        (d - 0.5) as usize
    }
}

fn java_jacobian(diff: f64) -> f64 {
    const STEP: f64 = 0.0001;
    const MAX_TOLERANCE: f64 = 8.0;
    let n = (MAX_TOLERANCE / STEP) as usize + 1;
    let idx = java_fast_round(diff * (1.0 / STEP));
    if idx >= n {
        return 0.0;
    }
    (1.0 + 10_f64.powf(-(idx as f64) * STEP)).log10()
}

fn java_pair(a: f64, b: f64) -> f64 {
    let (smaller, larger) = if a > b { (b, a) } else { (a, b) };
    if smaller == f64::NEG_INFINITY {
        return larger;
    }
    let diff = larger - smaller;
    if diff < 8.0 {
        larger + java_jacobian(diff)
    } else {
        larger
    }
}

fn bits(x: f64) -> u64 {
    x.to_bits()
}

#[test]
fn holdout_6r314_jacobian() {
    if !armed() {
        return;
    }
    let larger = -1.25_f64;
    let probes = [
        ("diff-0", larger, larger),
        ("diff-lt-8", larger, larger - 0.5),
        ("diff-eq-8", larger, larger - 8.0),
        ("diff-gt-8", larger, larger - 8.5),
        ("table-step", larger, larger - 0.0001),
        ("table-half-step", larger, larger - 0.00035),
    ];
    for (name, a, b) in probes {
        let prod = approximate_log10_sum_log10_pair(a, b);
        let java = java_pair(a, b);
        if bits(prod) != bits(java) {
            fail(
                "JACOBIAN_HOLDOUT_DIVERGENCE",
                name,
                &format!("prod {prod:.20} java {java:.20}"),
            );
        }
        println!("6R314\tcheckpoint\t{name}\tPASS");
    }
    let at_8 = approximate_log10_sum_log10_pair(larger, larger - 8.0);
    if bits(at_8) != bits(larger) {
        fail(
            "JACOBIAN_HOLDOUT_DIVERGENCE",
            "cutoff-8",
            "diff 8 changed the larger value",
        );
    }
    let above = approximate_log10_sum_log10_pair(larger - 9.0, larger);
    if bits(above) != bits(larger) {
        fail(
            "JACOBIAN_HOLDOUT_DIVERGENCE",
            "cutoff-above-8",
            "diff > 8 changed the larger value",
        );
    }
    let near = approximate_log10_sum_log10_pair(larger, larger);
    if bits(near) != bits(larger + 2.0_f64.log10()) {
        fail(
            "JACOBIAN_HOLDOUT_DIVERGENCE",
            "diff-0-closed-form",
            "equal arguments did not add log10(2)",
        );
    }
    let inf = approximate_log10_sum_log10_pair(f64::NEG_INFINITY, larger);
    if bits(inf) != bits(larger) {
        fail(
            "JACOBIAN_HOLDOUT_DIVERGENCE",
            "neg-inf",
            "negative infinity was not dropped",
        );
    }
    println!("6R314\tcheckpoint\tcutoff\tPASS");

    let row = ReadLikelihoodRow {
        read_index: 0,
        read_id: String::new(),
        haplotype_log10_likelihoods: vec![-1.0, -1.00015],
    };
    let gl = biallelic_genotype_log10_likelihoods_gatk(std::slice::from_ref(&row), 0, 1);
    if bits(gl[0]) != bits(-1.0) || bits(gl[2]) != bits(-1.00015) {
        fail(
            "JACOBIAN_HOLDOUT_DIVERGENCE",
            "homozygote",
            &format!("{gl:?}"),
        );
    }
    let het = java_pair(-1.0, -1.00015) - 2.0_f64.log10();
    if bits(gl[1]) != bits(het) {
        fail(
            "JACOBIAN_HOLDOUT_DIVERGENCE",
            "heterozygote",
            &format!("gl1 {} java {}", gl[1], het),
        );
    }
    let analytic = log10_sum_log10(&[-1.0, -1.00015]) - 2.0_f64.log10();
    if bits(gl[1]) == bits(analytic) {
        fail(
            "JACOBIAN_HOLDOUT_DIVERGENCE",
            "heterozygote-not-analytic",
            "Jacobian het matched analytic log10_sum_log10",
        );
    }
    println!("6R314\tcheckpoint\tbiallelic-hom-het\tPASS");

    let tri = ReadLikelihoodRow {
        read_index: 0,
        read_id: String::new(),
        haplotype_log10_likelihoods: vec![-1.0, -1.5, -4.0],
    };
    let dips = gatk_haplotypecaller::hc_genotyping_engine::diploid_genotype_log10_likelihoods_from_allele_rows(
        std::slice::from_ref(&tri),
        3,
    );
    let logp = 2.0_f64.log10();
    let expect = [
        -1.0,
        java_pair(-1.0, -1.5) - logp,
        -1.5,
        java_pair(-1.0, -4.0) - logp,
        java_pair(-1.5, -4.0) - logp,
        -4.0,
    ];
    for (i, (got, want)) in dips.iter().zip(expect.iter()).enumerate() {
        if bits(*got) != bits(*want) {
            fail(
                "JACOBIAN_HOLDOUT_DIVERGENCE",
                "diploid-3",
                &format!("index {i} got {got} want {want}"),
            );
        }
    }
    println!("6R314\tcheckpoint\tdiploid-3\tPASS");
    println!("6R314\tclassification\tJACOBIAN_HOLDOUT_PASS");
}

#[test]
fn holdout_6r314_composite() {
    if !armed() {
        return;
    }
    // chrHold:4500 G>GCAG. Not 20:29455649 and not the canonical A>AT pad.
    let event = 4500_u64;
    let bases = b"GCAGCAGCAG";
    let rep = tandem_repeat_at_event(event, bases, event, b"G", b"GCAG").unwrap_or_else(|| {
        fail(
            "COMPOSITE_HOLDOUT_DIVERGENCE",
            "trim-geometry",
            "no tandem repeat",
        )
    });
    if rep.unit.as_slice() != b"CAG" || rep.counts.as_slice() != [3, 4] {
        fail(
            "COMPOSITE_HOLDOUT_DIVERGENCE",
            "trim-geometry",
            &format!("{:?} {:?}", rep.unit, rep.counts),
        );
    }
    let pad = 75 + rep.padding_bases() as u32;
    let res = trimmer().trim(
        &region(),
        &[TrimVariant {
            contig: "chrHold".into(),
            start: event,
            end: event,
            is_indel: true,
            ref_allele: "G".into(),
            alt_allele: "GCAG".into(),
        }],
        Some(&ctx(event, bases)),
    );
    if res.padded_variant_start != Some(event - u64::from(pad))
        || res.padded_variant_end != Some(event + u64::from(pad))
    {
        fail(
            "COMPOSITE_HOLDOUT_DIVERGENCE",
            "trim-geometry",
            &format!(
                "{:?} {:?}",
                res.padded_variant_start, res.padded_variant_end
            ),
        );
    }
    println!("6R314\tcheckpoint\ttrim-geometry\tPASS");

    let forward = order_rec("zeta-forward", 99, 0, 2000, 20, 0, 10, 100);
    let reverse = order_rec("aaa-reverse", 83, 0, 2000, 20, 0, 10, 100);
    let mut cross = order_rec("hold314-cross", 99, 0, 2000, 20, 1, 10, 100);
    cross.set_mtid(1);
    let originals = [forward.clone(), reverse.clone(), cross.clone()];
    let keys = failed_mate_evidence_keys(&originals);
    if keys.len() != 1 || !keys.contains(&(b"hold314-cross".to_vec(), cross.flags())) {
        fail(
            "COMPOSITE_HOLDOUT_DIVERGENCE",
            "evidence-membership",
            "mate keys",
        );
    }
    let mut evidence = originals.to_vec();
    evidence.retain(|r| !keys.contains(&(r.qname().to_vec(), r.flags())));
    if evidence.len() != 2 || evidence.iter().any(|r| r.qname() == b"hold314-cross") {
        fail(
            "COMPOSITE_HOLDOUT_DIVERGENCE",
            "evidence-membership",
            "cross read still present",
        );
    }
    println!("6R314\tcheckpoint\tevidence-membership\tPASS");

    set_forensic_6r308_java_read_coordinate_order(true);
    let mut clipped = evidence;
    clip_finalized_reads_in_place(&mut clipped, &order_region());
    if clipped.len() != 2
        || clipped[0].qname() != b"zeta-forward"
        || clipped[0].flags() != 99
        || clipped[1].qname() != b"aaa-reverse"
        || clipped[1].flags() != 83
    {
        fail(
            "COMPOSITE_HOLDOUT_DIVERGENCE",
            "clipped-read-order",
            "forward did not precede reverse",
        );
    }
    println!("6R314\tcheckpoint\tclipped-read-order\tPASS");

    let rows = [
        ReadLikelihoodRow {
            read_index: 0,
            read_id: String::new(),
            haplotype_log10_likelihoods: vec![-3.0, -3.4],
        },
        ReadLikelihoodRow {
            read_index: 1,
            read_id: String::new(),
            haplotype_log10_likelihoods: vec![-1.0, -5.0],
        },
    ];
    let gl = biallelic_genotype_log10_likelihoods_gatk(&rows, 0, 1);
    let denom = 2.0 * 2.0_f64.log10();
    let want = [
        (-3.0 - 1.0 + 2.0 * 2.0_f64.log10()) - denom,
        (java_pair(-3.0, -3.4) + java_pair(-1.0, -5.0)) - denom,
        (-3.4 - 5.0 + 2.0 * 2.0_f64.log10()) - denom,
    ];
    if gl.len() != 3
        || gl
            .iter()
            .zip(want.iter())
            .any(|(g, w)| bits(*g) != bits(*w))
    {
        fail(
            "COMPOSITE_HOLDOUT_DIVERGENCE",
            "genotype-likelihoods",
            &format!("{gl:?} vs {want:?}"),
        );
    }
    println!("6R314\tcheckpoint\tgenotype-likelihoods\tPASS");
    println!("6R314\tclassification\tCOMPOSITE_HOLDOUT_PASS");
}
