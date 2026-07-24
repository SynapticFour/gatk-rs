use gatk_haplotypecaller::{filtered_read_iteration_order, ReadFilterParams};
use std::path::PathBuf;

fn read_order_sam() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../parity/fixtures/read_order.sam")
}

#[test]
fn ingress_order_contract_is_stable_across_thread_count_env() {
    let path = read_order_sam();
    let p = ReadFilterParams {
        min_mapping_quality: 0,
        exclude_duplicates: false,
        exclude_secondary: false,
        exclude_supplementary: false,
    };

    // Current ingress path is single-threaded; this guards against accidental
    // future non-determinism introduced via parallelized ingestion.
    std::env::set_var("RAYON_NUM_THREADS", "1");
    let one_thread = filtered_read_iteration_order(&path, &p).unwrap();
    std::env::set_var("RAYON_NUM_THREADS", "8");
    let eight_threads = filtered_read_iteration_order(&path, &p).unwrap();
    assert_eq!(one_thread, eight_threads);
    assert_eq!(
        one_thread
            .iter()
            .map(|(_, _, q)| q.as_str())
            .collect::<Vec<_>>(),
        vec!["ordA", "ordB", "ordC", "ordD"]
    );
}
