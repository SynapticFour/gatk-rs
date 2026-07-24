use gatk_haplotypecaller::{pairhmm_log10_likelihoods_vectorized, PairHmmParams};
use std::collections::BTreeMap;
use std::fs;
use std::path::PathBuf;

#[derive(Debug, Clone)]
struct CaseRow {
    case_id: String,
    read_bases: String,
    read_base_quals: Vec<u8>,
    read_mapq: u8,
    haplotype: String,
}

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("workspace root")
        .to_path_buf()
}

fn load_cases(path: &PathBuf) -> Vec<CaseRow> {
    fs::read_to_string(path)
        .expect("cases fixture")
        .lines()
        .filter(|l| !l.trim().is_empty() && !l.starts_with('#'))
        .map(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            let quals = c[2]
                .split(',')
                .map(|q| q.parse::<u8>().expect("qual"))
                .collect::<Vec<_>>();
            CaseRow {
                case_id: c[0].to_string(),
                read_bases: c[1].to_string(),
                read_base_quals: quals,
                read_mapq: c[3].parse::<u8>().expect("mapq"),
                haplotype: c[4].to_string(),
            }
        })
        .collect()
}

fn load_expected(path: &PathBuf) -> BTreeMap<String, f64> {
    fs::read_to_string(path)
        .expect("expected fixture")
        .lines()
        .filter(|l| {
            let t = l.trim();
            !t.is_empty() && !t.starts_with('#') && !t.starts_with("case_id\t")
        })
        .map(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            (
                c[0].to_string(),
                c[1].parse::<f64>().expect("expected likelihood"),
            )
        })
        .collect()
}

#[test]
fn pairhmm_likelihood_vector_matches_frozen_java_dump_fixture() {
    let root = repo_root();
    let cases = load_cases(&root.join("parity/fixtures/p6_pairhmm_case1_reads.tsv"));
    let expected =
        load_expected(&root.join("parity/expected/p6_pairhmm_case1.java_likelihoods.tsv"));

    assert!(!cases.is_empty());
    let read_bases = &cases[0].read_bases;
    let read_base_quals = &cases[0].read_base_quals;
    let read_mapq = cases[0].read_mapq;
    let haplotypes = cases
        .iter()
        .map(|c| c.haplotype.clone())
        .collect::<Vec<_>>();

    let got = pairhmm_log10_likelihoods_vectorized(
        read_bases,
        read_base_quals,
        read_mapq,
        &haplotypes,
        &PairHmmParams::default(),
    )
    .expect("vectorized likelihoods");

    for (idx, case) in cases.iter().enumerate() {
        let exp = expected
            .get(&case.case_id)
            .unwrap_or_else(|| panic!("missing expected for {}", case.case_id));
        let delta = (got[idx] - exp).abs();
        assert!(
            delta <= 1e-9,
            "case={} got={} expected={} delta={}",
            case.case_id,
            got[idx],
            exp,
            delta
        );
    }
}
