use gatk_haplotypecaller::{passes_hc_read_filters_fields, ReadFilterParams};

#[derive(Clone, Copy)]
struct Fixture {
    flags: u16,
    mapq: u8,
    java_expected_accept: bool,
}

#[test]
fn read_filter_synthetic_diff_contract_matches_java_fixtures() {
    // Expected outcomes are frozen from Java GATK HaplotypeCaller read-filter behavior
    // for a synthetic fixture matrix (Phase 2, Step 33 target harness).
    let fixtures = [
        Fixture {
            flags: 0,
            mapq: 60,
            java_expected_accept: true,
        },
        Fixture {
            flags: 0x0400, // duplicate
            mapq: 60,
            java_expected_accept: false,
        },
        Fixture {
            flags: 0x0100, // secondary
            mapq: 60,
            java_expected_accept: false,
        },
        Fixture {
            flags: 0x0800, // supplementary
            mapq: 60,
            java_expected_accept: false,
        },
        Fixture {
            flags: 0x0004, // unmapped
            mapq: 60,
            java_expected_accept: false,
        },
        Fixture {
            flags: 0,
            mapq: 19,
            java_expected_accept: false,
        },
        Fixture {
            flags: 0,
            mapq: 20,
            java_expected_accept: true,
        },
        Fixture {
            flags: 0,
            mapq: 255,
            java_expected_accept: true,
        },
    ];

    let params = ReadFilterParams {
        min_mapping_quality: 20,
        exclude_duplicates: true,
        exclude_secondary: true,
        exclude_supplementary: true,
    };

    for f in fixtures {
        let rust_accept = passes_hc_read_filters_fields(f.flags, f.mapq, &params);
        assert_eq!(
            rust_accept, f.java_expected_accept,
            "read-filter parity mismatch for flags=0x{:04x}, mapq={}",
            f.flags, f.mapq
        );
    }
}
