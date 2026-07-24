use gatk_haplotypecaller::{dump_hc_read_filter_tsv, HC_READ_FILTER_COUNT_SECTION};
use std::path::Path;

#[test]
fn read_filter_dump_appends_counting_read_filter_section() {
    let sam =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../parity/fixtures/read_filter_slice.sam");
    let mut buf = Vec::new();
    dump_hc_read_filter_tsv(&sam, &mut buf).unwrap();
    let s = String::from_utf8(buf).unwrap();
    assert!(
        s.contains(HC_READ_FILTER_COUNT_SECTION),
        "expected delimiter {HC_READ_FILTER_COUNT_SECTION:?}"
    );
    assert!(s.contains("MappingQualityReadFilter\t2"));
    assert!(s.contains("WellformedReadFilter\t2"));
    assert!(s.lines().filter(|l| l == &"MappedReadFilter\t0").count() == 1);
}
