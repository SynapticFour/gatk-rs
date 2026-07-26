use gatk_haplotypecaller::{
    build_gvcf_blocks, gvcf_block_to_record_fields, GvcfBlockRecordFields, ReferenceConfidenceLocus,
};
use std::collections::BTreeMap;
use std::fs;
use std::path::PathBuf;

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("workspace root")
        .to_path_buf()
}

fn load_loci_by_case(path: &PathBuf) -> BTreeMap<String, Vec<ReferenceConfidenceLocus>> {
    let mut out = BTreeMap::<String, Vec<ReferenceConfidenceLocus>>::new();
    for line in fs::read_to_string(path)
        .expect("fixture")
        .lines()
        .filter(|l| !l.trim().is_empty() && !l.starts_with('#'))
    {
        let c: Vec<&str> = line.split('\t').collect();
        out.entry(c[0].to_string())
            .or_default()
            .push(ReferenceConfidenceLocus {
                position_1based: c[1].parse::<usize>().expect("position"),
                gq: c[2].parse::<i32>().expect("gq"),
                dp: c[3].parse::<i32>().expect("dp"),
            });
    }
    out
}

fn load_expected(path: &PathBuf) -> BTreeMap<String, Vec<GvcfBlockRecordFields>> {
    let mut out = BTreeMap::<String, Vec<GvcfBlockRecordFields>>::new();
    for line in fs::read_to_string(path)
        .expect("expected")
        .lines()
        .filter(|l| !l.trim().is_empty() && !l.starts_with('#'))
    {
        let c: Vec<&str> = line.split('\t').collect();
        out.entry(c[0].to_string())
            .or_default()
            .push(GvcfBlockRecordFields {
                start_1based: c[1].parse::<usize>().expect("start"),
                end_info: c[2].parse::<usize>().expect("end"),
                gq_band_upper: c[3].parse::<i32>().expect("band"),
                min_rgq: c[4].parse::<i32>().expect("min_rgq"),
                min_dp: c[5].parse::<i32>().expect("min_dp"),
                max_dp: c[6].parse::<i32>().expect("max_dp"),
            });
    }
    out
}

#[test]
fn gvcf_block_emission_matches_frozen_java_smoke_fixture() {
    let root = repo_root();
    let by_case = load_loci_by_case(&root.join("parity/fixtures/p8_gvcf_blocks_smoke.tsv"));
    let expected = load_expected(&root.join("parity/expected/p8_gvcf_blocks_smoke.java.tsv"));
    let gq_bands = [9, 19, 29, 99];

    assert_eq!(by_case.len(), expected.len());
    for (case_id, loci) in &by_case {
        let blocks = build_gvcf_blocks(loci, &gq_bands).expect("build blocks");
        let got = blocks
            .iter()
            .map(|b| gvcf_block_to_record_fields(b).expect("record fields"))
            .collect::<Vec<_>>();
        let exp = expected
            .get(case_id)
            .unwrap_or_else(|| panic!("missing expected for case {case_id}"));
        assert_eq!(got, *exp, "case={case_id}");
    }
}
