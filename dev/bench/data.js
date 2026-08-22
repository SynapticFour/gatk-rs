window.BENCHMARK_DATA = {
  "lastUpdate": 1787369253529,
  "repoUrl": "https://github.com/SynapticFour/gatk-rs",
  "entries": {
    "Benchmark": [
      {
        "commit": {
          "author": {
            "name": "SynapticFour",
            "username": "SynapticFour",
            "email": "contact@synapticfour.com"
          },
          "committer": {
            "name": "GitHub",
            "username": "web-flow",
            "email": "noreply@github.com"
          },
          "id": "0a450d3059bc81a543d6b32e66676a79c4b8766c",
          "message": "deps: bundle Dependabot (itertools + GHA bumps) (#80)\n\nBundles green Dependabot updates (#73, #75–#79). Leaves #74 rust-htslib 1.x for a separate fix.",
          "timestamp": "2026-08-03T15:48:05Z",
          "url": "https://github.com/SynapticFour/gatk-rs/commit/0a450d3059bc81a543d6b32e66676a79c4b8766c"
        },
        "date": 1785816598181,
        "tool": "cargo",
        "benches": [
          {
            "name": "sam_parsing/parsing/100",
            "value": 88624,
            "range": "± 5257",
            "unit": "ns/iter"
          },
          {
            "name": "sam_parsing/parsing/1000",
            "value": 965421,
            "range": "± 42109",
            "unit": "ns/iter"
          },
          {
            "name": "sam_parsing/parsing/10000",
            "value": 9760822,
            "range": "± 445514",
            "unit": "ns/iter"
          },
          {
            "name": "sam_iterator/iterator",
            "value": 681938,
            "range": "± 37896",
            "unit": "ns/iter"
          },
          {
            "name": "sam_writing/writing",
            "value": 1877849,
            "range": "± 5347382",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/0",
            "value": 116,
            "range": "± 5",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/1",
            "value": 160,
            "range": "± 8",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/2",
            "value": 125,
            "range": "± 7",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/3",
            "value": 137,
            "range": "± 7",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/4",
            "value": 160,
            "range": "± 7",
            "unit": "ns/iter"
          },
          {
            "name": "optional_fields/parsing",
            "value": 3483,
            "range": "± 187",
            "unit": "ns/iter"
          },
          {
            "name": "framework_overhead/create_result",
            "value": 22,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "framework_overhead/create_comparison",
            "value": 0,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "metrics_collection/metrics_collector",
            "value": 1058241,
            "range": "± 480",
            "unit": "ns/iter"
          },
          {
            "name": "suite_operations/create_suite",
            "value": 11066095,
            "range": "± 602878",
            "unit": "ns/iter"
          },
          {
            "name": "suite_operations/suite_with_results",
            "value": 10461056,
            "range": "± 494057",
            "unit": "ns/iter"
          },
          {
            "name": "suite_operations/suite_serialization",
            "value": 56825,
            "range": "± 2622",
            "unit": "ns/iter"
          },
          {
            "name": "performance_analysis/analyze_small_dataset",
            "value": 63,
            "range": "± 3",
            "unit": "ns/iter"
          },
          {
            "name": "performance_analysis/analyze_large_dataset",
            "value": 8344,
            "range": "± 430",
            "unit": "ns/iter"
          },
          {
            "name": "regression_detection/detect_no_regression",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "regression_detection/detect_regression",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "report_generation/generate_json_report",
            "value": 164207,
            "range": "± 72228",
            "unit": "ns/iter"
          },
          {
            "name": "report_generation/generate_markdown_report",
            "value": 176491,
            "range": "± 401051",
            "unit": "ns/iter"
          },
          {
            "name": "report_generation/generate_csv_report",
            "value": 164188,
            "range": "± 90294",
            "unit": "ns/iter"
          },
          {
            "name": "dataset_operations/create_dataset_info",
            "value": 67,
            "range": "± 4",
            "unit": "ns/iter"
          },
          {
            "name": "dataset_operations/validate_dataset_info",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "comparison_operations/create_comparison",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "comparison_operations/check_targets",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_gc_content",
            "value": 747883,
            "range": "± 38842",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_revcomp",
            "value": 310077,
            "range": "± 20056",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_subsequence",
            "value": 104,
            "range": "± 6",
            "unit": "ns/iter"
          },
          {
            "name": "memory_pool_allocation",
            "value": 54,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "cache_put_get",
            "value": 76298,
            "range": "± 4128",
            "unit": "ns/iter"
          },
          {
            "name": "interval_query",
            "value": 29472173,
            "range": "± 1312919",
            "unit": "ns/iter"
          },
          {
            "name": "variant_type",
            "value": 13580,
            "range": "± 713",
            "unit": "ns/iter"
          },
          {
            "name": "quality_error_prob",
            "value": 1049,
            "range": "± 53",
            "unit": "ns/iter"
          },
          {
            "name": "log_addition",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "stream_processing",
            "value": 356,
            "range": "± 18",
            "unit": "ns/iter"
          },
          {
            "name": "memory_mapped_access",
            "value": 3242,
            "range": "± 195",
            "unit": "ns/iter"
          },
          {
            "name": "hamming_distance",
            "value": 7,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/100",
            "value": 23815,
            "range": "± 1307",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/100",
            "value": 41377,
            "range": "± 2593",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/100",
            "value": 21191,
            "range": "± 990",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/1000",
            "value": 226598,
            "range": "± 11933",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/1000",
            "value": 318068,
            "range": "± 17355",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/1000",
            "value": 191547,
            "range": "± 9288",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/10000",
            "value": 2299977,
            "range": "± 102926",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/10000",
            "value": 3121078,
            "range": "± 137776",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/10000",
            "value": 1909377,
            "range": "± 101468",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/100",
            "value": 26187,
            "range": "± 1283",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/100",
            "value": 35814,
            "range": "± 1535",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/100",
            "value": 19919,
            "range": "± 1017",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/1000",
            "value": 235229,
            "range": "± 9044",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/1000",
            "value": 256886,
            "range": "± 10775",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/1000",
            "value": 175147,
            "range": "± 8408",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/10000",
            "value": 2383035,
            "range": "± 84595",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/10000",
            "value": 2447227,
            "range": "± 120806",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/10000",
            "value": 1664363,
            "range": "± 79818",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/gc_content",
            "value": 50,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/reverse_complement",
            "value": 49,
            "range": "± 3",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/subsequence",
            "value": 17,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/count_pattern",
            "value": 150,
            "range": "± 9",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/average_quality",
            "value": 27,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/min_quality",
            "value": 53,
            "range": "± 3",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/max_quality",
            "value": 63,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/median_quality",
            "value": 54,
            "range": "± 3",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/trim_quality",
            "value": 46,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/reverse_complement",
            "value": 86,
            "range": "± 4",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/is_valid",
            "value": 59,
            "range": "± 3",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/build_index",
            "value": 298596,
            "range": "± 12425",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/random_access",
            "value": 1936,
            "range": "± 87",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/save_index",
            "value": 605350,
            "range": "± 527220",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/load_index",
            "value": 227169,
            "range": "± 8325",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/filter_by_quality",
            "value": 280081,
            "range": "± 11463",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/filter_by_length",
            "value": 170232,
            "range": "± 8696",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/sample_reads",
            "value": 168810,
            "range": "± 10713",
            "unit": "ns/iter"
          },
          {
            "name": "io_writing/write_fasta",
            "value": 369002,
            "range": "± 670390",
            "unit": "ns/iter"
          },
          {
            "name": "io_writing/write_fastq",
            "value": 4105006,
            "range": "± 4591275",
            "unit": "ns/iter"
          },
          {
            "name": "memory_usage/parse_memory_usage",
            "value": 2344523,
            "range": "± 148736",
            "unit": "ns/iter"
          },
          {
            "name": "macro_sum_bench",
            "value": 1,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/100",
            "value": 76712,
            "range": "± 4455",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/1000",
            "value": 862497,
            "range": "± 41383",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/10000",
            "value": 8961561,
            "range": "± 453400",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_iterator/iterator",
            "value": 540542,
            "range": "± 27911",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_writing/writing",
            "value": 2265831,
            "range": "± 834390",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/0",
            "value": 75,
            "range": "± 3",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/1",
            "value": 76,
            "range": "± 4",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/2",
            "value": 72,
            "range": "± 4",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/3",
            "value": 100,
            "range": "± 5",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/4",
            "value": 84,
            "range": "± 4",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/5",
            "value": 43,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/6",
            "value": 97,
            "range": "± 5",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/7",
            "value": 495,
            "range": "± 22",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/0",
            "value": 68,
            "range": "± 3",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/1",
            "value": 68,
            "range": "± 3",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/2",
            "value": 31,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/3",
            "value": 68,
            "range": "± 3",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/4",
            "value": 68,
            "range": "± 3",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/5",
            "value": 71,
            "range": "± 4",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/6",
            "value": 70,
            "range": "± 4",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/7",
            "value": 67,
            "range": "± 3",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/variant_type_detection",
            "value": 2049,
            "range": "± 110",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/allele_frequency_access",
            "value": 2611,
            "range": "± 147",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/sample_data_access",
            "value": 1086,
            "range": "± 66",
            "unit": "ns/iter"
          },
          {
            "name": "band_pass_normalized_kernel_gatk_defaults",
            "value": 812,
            "range": "± 25",
            "unit": "ns/iter"
          },
          {
            "name": "band_pass_add_then_pop_ready_256_loci",
            "value": 263605,
            "range": "± 7615",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_medium_depth_k10",
            "value": 172765,
            "range": "± 4335",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_high_depth_k10",
            "value": 1322106,
            "range": "± 47331",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/4",
            "value": 11444408,
            "range": "± 622572",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/16",
            "value": 41655596,
            "range": "± 2132820",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/64",
            "value": 125915231,
            "range": "± 12330128",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_haps/8",
            "value": 1562175,
            "range": "± 66644",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_haps/8",
            "value": 3083931,
            "range": "± 136372",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_haps/8",
            "value": 3193145,
            "range": "± 154380",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_haps/8",
            "value": 30393088,
            "range": "± 1710988",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_haps/32",
            "value": 6317718,
            "range": "± 265432",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_haps/32",
            "value": 12529715,
            "range": "± 642168",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_haps/32",
            "value": 12841326,
            "range": "± 576229",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_haps/32",
            "value": 120880907,
            "range": "± 4536904",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_haps/64",
            "value": 12616598,
            "range": "± 537182",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_haps/64",
            "value": 24867699,
            "range": "± 1041363",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_haps/64",
            "value": 25270944,
            "range": "± 1438423",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_haps/64",
            "value": 243974503,
            "range": "± 9799195",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/value_copy_stages/3072x12",
            "value": 4110,
            "range": "± 201",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/arc_share_stages/3072x12",
            "value": 3141,
            "range": "± 150",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/assembly_result_set_arc_fanout/3072",
            "value": 160,
            "range": "± 9",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/assembly_result_set_value_fanout/3072",
            "value": 1060,
            "range": "± 57",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/64x48",
            "value": 11420,
            "range": "± 600",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/64x48",
            "value": 12387,
            "range": "± 2582",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/128x96",
            "value": 50833,
            "range": "± 2743",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/128x96",
            "value": 44955,
            "range": "± 2084",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/256x192",
            "value": 180851,
            "range": "± 9252",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/256x192",
            "value": 169449,
            "range": "± 7954",
            "unit": "ns/iter"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "name": "SynapticFour",
            "username": "SynapticFour",
            "email": "contact@synapticfour.com"
          },
          "committer": {
            "name": "GitHub",
            "username": "web-flow",
            "email": "noreply@github.com"
          },
          "id": "fa35f207af1f91b24bf55368d280cb42d5e73307",
          "message": "perf(hc): Arc-share BAM reads and shrink Peak-RSS for GIAB (#81)\n\nArc-share BAM reads, sequential regions, TLS shrink, GIAB disk hygiene for Peak-RSS on 16 GiB hosts.",
          "timestamp": "2026-08-04T16:49:06Z",
          "url": "https://github.com/SynapticFour/gatk-rs/commit/fa35f207af1f91b24bf55368d280cb42d5e73307"
        },
        "date": 1785903281339,
        "tool": "cargo",
        "benches": [
          {
            "name": "sam_parsing/parsing/100",
            "value": 108465,
            "range": "± 4842",
            "unit": "ns/iter"
          },
          {
            "name": "sam_parsing/parsing/1000",
            "value": 1122619,
            "range": "± 13123",
            "unit": "ns/iter"
          },
          {
            "name": "sam_parsing/parsing/10000",
            "value": 10747984,
            "range": "± 169942",
            "unit": "ns/iter"
          },
          {
            "name": "sam_iterator/iterator",
            "value": 879243,
            "range": "± 15081",
            "unit": "ns/iter"
          },
          {
            "name": "sam_writing/writing",
            "value": 3290592,
            "range": "± 244691",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/0",
            "value": 162,
            "range": "± 12",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/1",
            "value": 199,
            "range": "± 7",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/2",
            "value": 164,
            "range": "± 5",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/3",
            "value": 169,
            "range": "± 5",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/4",
            "value": 202,
            "range": "± 12",
            "unit": "ns/iter"
          },
          {
            "name": "optional_fields/parsing",
            "value": 7298,
            "range": "± 197",
            "unit": "ns/iter"
          },
          {
            "name": "framework_overhead/create_result",
            "value": 27,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "framework_overhead/create_comparison",
            "value": 1,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "metrics_collection/metrics_collector",
            "value": 1062076,
            "range": "± 1091",
            "unit": "ns/iter"
          },
          {
            "name": "suite_operations/create_suite",
            "value": 17422390,
            "range": "± 438278",
            "unit": "ns/iter"
          },
          {
            "name": "suite_operations/suite_with_results",
            "value": 17137714,
            "range": "± 430624",
            "unit": "ns/iter"
          },
          {
            "name": "suite_operations/suite_serialization",
            "value": 74376,
            "range": "± 2283",
            "unit": "ns/iter"
          },
          {
            "name": "performance_analysis/analyze_small_dataset",
            "value": 83,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "performance_analysis/analyze_large_dataset",
            "value": 8855,
            "range": "± 83",
            "unit": "ns/iter"
          },
          {
            "name": "regression_detection/detect_no_regression",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "regression_detection/detect_regression",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "report_generation/generate_json_report",
            "value": 180259,
            "range": "± 27172",
            "unit": "ns/iter"
          },
          {
            "name": "report_generation/generate_markdown_report",
            "value": 98115,
            "range": "± 11757",
            "unit": "ns/iter"
          },
          {
            "name": "report_generation/generate_csv_report",
            "value": 95789,
            "range": "± 7277",
            "unit": "ns/iter"
          },
          {
            "name": "dataset_operations/create_dataset_info",
            "value": 81,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "dataset_operations/validate_dataset_info",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "comparison_operations/create_comparison",
            "value": 1,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "comparison_operations/check_targets",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_gc_content",
            "value": 814284,
            "range": "± 4592",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_revcomp",
            "value": 297683,
            "range": "± 7921",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_subsequence",
            "value": 164,
            "range": "± 7",
            "unit": "ns/iter"
          },
          {
            "name": "memory_pool_allocation",
            "value": 71,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "cache_put_get",
            "value": 111472,
            "range": "± 908",
            "unit": "ns/iter"
          },
          {
            "name": "interval_query",
            "value": 41749499,
            "range": "± 1312549",
            "unit": "ns/iter"
          },
          {
            "name": "variant_type",
            "value": 16203,
            "range": "± 178",
            "unit": "ns/iter"
          },
          {
            "name": "quality_error_prob",
            "value": 1214,
            "range": "± 26",
            "unit": "ns/iter"
          },
          {
            "name": "log_addition",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "stream_processing",
            "value": 560,
            "range": "± 8",
            "unit": "ns/iter"
          },
          {
            "name": "memory_mapped_access",
            "value": 9715,
            "range": "± 170",
            "unit": "ns/iter"
          },
          {
            "name": "hamming_distance",
            "value": 9,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/100",
            "value": 29446,
            "range": "± 1007",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/100",
            "value": 51414,
            "range": "± 692",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/100",
            "value": 26459,
            "range": "± 353",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/1000",
            "value": 242941,
            "range": "± 3331",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/1000",
            "value": 379974,
            "range": "± 11817",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/1000",
            "value": 211149,
            "range": "± 3897",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/10000",
            "value": 2459955,
            "range": "± 21456",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/10000",
            "value": 3472925,
            "range": "± 38333",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/10000",
            "value": 2094964,
            "range": "± 45538",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/100",
            "value": 33250,
            "range": "± 87",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/100",
            "value": 44330,
            "range": "± 200",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/100",
            "value": 27251,
            "range": "± 961",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/1000",
            "value": 273014,
            "range": "± 5788",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/1000",
            "value": 299960,
            "range": "± 818",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/1000",
            "value": 221201,
            "range": "± 12898",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/10000",
            "value": 2675385,
            "range": "± 35384",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/10000",
            "value": 2954235,
            "range": "± 77015",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/10000",
            "value": 2185436,
            "range": "± 93090",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/gc_content",
            "value": 63,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/reverse_complement",
            "value": 59,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/subsequence",
            "value": 23,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/count_pattern",
            "value": 225,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/average_quality",
            "value": 62,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/min_quality",
            "value": 45,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/max_quality",
            "value": 49,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/median_quality",
            "value": 45,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/trim_quality",
            "value": 53,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/reverse_complement",
            "value": 94,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/is_valid",
            "value": 92,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/build_index",
            "value": 342222,
            "range": "± 15123",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/random_access",
            "value": 4901,
            "range": "± 25",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/save_index",
            "value": 223781,
            "range": "± 34946",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/load_index",
            "value": 303795,
            "range": "± 833",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/filter_by_quality",
            "value": 332712,
            "range": "± 3739",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/filter_by_length",
            "value": 187213,
            "range": "± 2107",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/sample_reads",
            "value": 185736,
            "range": "± 657",
            "unit": "ns/iter"
          },
          {
            "name": "io_writing/write_fasta",
            "value": 216169,
            "range": "± 24164",
            "unit": "ns/iter"
          },
          {
            "name": "io_writing/write_fastq",
            "value": 8770881,
            "range": "± 780146",
            "unit": "ns/iter"
          },
          {
            "name": "memory_usage/parse_memory_usage",
            "value": 2462150,
            "range": "± 11118",
            "unit": "ns/iter"
          },
          {
            "name": "macro_sum_bench",
            "value": 1,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/100",
            "value": 90948,
            "range": "± 2161",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/1000",
            "value": 876021,
            "range": "± 7682",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/10000",
            "value": 9448102,
            "range": "± 166301",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_iterator/iterator",
            "value": 626120,
            "range": "± 6275",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_writing/writing",
            "value": 3875459,
            "range": "± 314304",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/0",
            "value": 106,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/1",
            "value": 96,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/2",
            "value": 93,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/3",
            "value": 134,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/4",
            "value": 113,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/5",
            "value": 59,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/6",
            "value": 114,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/7",
            "value": 476,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/0",
            "value": 91,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/1",
            "value": 91,
            "range": "± 3",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/2",
            "value": 47,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/3",
            "value": 91,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/4",
            "value": 91,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/5",
            "value": 91,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/6",
            "value": 91,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/7",
            "value": 91,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/variant_type_detection",
            "value": 2472,
            "range": "± 3",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/allele_frequency_access",
            "value": 4235,
            "range": "± 19",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/sample_data_access",
            "value": 1527,
            "range": "± 24",
            "unit": "ns/iter"
          },
          {
            "name": "band_pass_normalized_kernel_gatk_defaults",
            "value": 991,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "band_pass_add_then_pop_ready_256_loci",
            "value": 327220,
            "range": "± 1220",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_medium_depth_k10",
            "value": 252407,
            "range": "± 3151",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_high_depth_k10",
            "value": 1753947,
            "range": "± 74335",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/4",
            "value": 13896274,
            "range": "± 872750",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/16",
            "value": 47231120,
            "range": "± 794282",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/64",
            "value": 139272704,
            "range": "± 19748313",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_haps/8",
            "value": 2129199,
            "range": "± 27967",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_haps/8",
            "value": 4940475,
            "range": "± 63790",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_haps/8",
            "value": 4341779,
            "range": "± 77074",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_haps/8",
            "value": 36472596,
            "range": "± 730728",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_haps/32",
            "value": 8502107,
            "range": "± 32314",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_haps/32",
            "value": 19636393,
            "range": "± 462695",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_haps/32",
            "value": 17367917,
            "range": "± 22683",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_haps/32",
            "value": 145812527,
            "range": "± 1053846",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_haps/64",
            "value": 17002365,
            "range": "± 23133",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_haps/64",
            "value": 38744777,
            "range": "± 204087",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_haps/64",
            "value": 34736903,
            "range": "± 33329",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_haps/64",
            "value": 291617254,
            "range": "± 1740461",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/value_copy_stages/3072x12",
            "value": 5089,
            "range": "± 121",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/arc_share_stages/3072x12",
            "value": 3936,
            "range": "± 32",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/assembly_result_set_arc_fanout/3072",
            "value": 56,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/assembly_result_set_value_fanout/3072",
            "value": 1245,
            "range": "± 15",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/64x48",
            "value": 17719,
            "range": "± 429",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/64x48",
            "value": 17543,
            "range": "± 715",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/128x96",
            "value": 70039,
            "range": "± 358",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/128x96",
            "value": 69475,
            "range": "± 2870",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/256x192",
            "value": 268732,
            "range": "± 1020",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/256x192",
            "value": 266068,
            "range": "± 649",
            "unit": "ns/iter"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "name": "SynapticFour",
            "username": "SynapticFour",
            "email": "contact@synapticfour.com"
          },
          "committer": {
            "name": "GitHub",
            "username": "web-flow",
            "email": "noreply@github.com"
          },
          "id": "79100ddcf72b345136fead0c4edbe5be18306c03",
          "message": "perf(hc): Phase A ownership showcase for Peak-RSS (#83)\n\n* perf(hc): score PairHMM without AssemblyRead Strings and reuse empty BAM sentinel\n\nCut Peak-RSS rematerialization on the likelihood path and progressive-release\nfragmentation; route GATK_RS_HC_SEQUENTIAL via runtime_config for N-2.\n\nCo-authored-by: Cursor <cursoragent@cursor.com>\n\n* perf(hc): Phase A ownership showcase for Peak-RSS\n\nByte-native assembly, single finalize buffer, COW unique realign ownership,\nsequential hap scoring under GATK_RS_HC_SEQUENTIAL, and optional jemalloc —\nwithout flipping PairHMM defaults or widening P12 bands.\n\nCo-authored-by: Cursor <cursoragent@cursor.com>\n\n* perf(hc): Phase B SIMD pack reuse + measured Peak-RSS\n\nReuse Logless DP scratch across haplotypes, borrow k-best graphs when acyclic,\nextend phenotype benches/tests, and record bomb/50kb/holdout RSS (100kb still\naborts — leave GIAB ci-subset unsigned).\n\nCo-authored-by: Cursor <cursoragent@cursor.com>\n\n---------\n\nCo-authored-by: Cursor <cursoragent@cursor.com>",
          "timestamp": "2026-08-05T17:24:41Z",
          "url": "https://github.com/SynapticFour/gatk-rs/commit/79100ddcf72b345136fead0c4edbe5be18306c03"
        },
        "date": 1785990315477,
        "tool": "cargo",
        "benches": [
          {
            "name": "sam_parsing/parsing/100",
            "value": 111814,
            "range": "± 1170",
            "unit": "ns/iter"
          },
          {
            "name": "sam_parsing/parsing/1000",
            "value": 1125004,
            "range": "± 68144",
            "unit": "ns/iter"
          },
          {
            "name": "sam_parsing/parsing/10000",
            "value": 11046527,
            "range": "± 59315",
            "unit": "ns/iter"
          },
          {
            "name": "sam_iterator/iterator",
            "value": 889326,
            "range": "± 7463",
            "unit": "ns/iter"
          },
          {
            "name": "sam_writing/writing",
            "value": 3224459,
            "range": "± 114201",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/0",
            "value": 124,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/1",
            "value": 172,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/2",
            "value": 134,
            "range": "± 4",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/3",
            "value": 141,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/4",
            "value": 171,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "optional_fields/parsing",
            "value": 6491,
            "range": "± 46",
            "unit": "ns/iter"
          },
          {
            "name": "framework_overhead/create_result",
            "value": 27,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "framework_overhead/create_comparison",
            "value": 1,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "metrics_collection/metrics_collector",
            "value": 1065408,
            "range": "± 1572",
            "unit": "ns/iter"
          },
          {
            "name": "suite_operations/create_suite",
            "value": 15262687,
            "range": "± 379062",
            "unit": "ns/iter"
          },
          {
            "name": "suite_operations/suite_with_results",
            "value": 15412669,
            "range": "± 322313",
            "unit": "ns/iter"
          },
          {
            "name": "suite_operations/suite_serialization",
            "value": 75894,
            "range": "± 1907",
            "unit": "ns/iter"
          },
          {
            "name": "performance_analysis/analyze_small_dataset",
            "value": 82,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "performance_analysis/analyze_large_dataset",
            "value": 7604,
            "range": "± 50",
            "unit": "ns/iter"
          },
          {
            "name": "regression_detection/detect_no_regression",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "regression_detection/detect_regression",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "report_generation/generate_json_report",
            "value": 254478,
            "range": "± 16323",
            "unit": "ns/iter"
          },
          {
            "name": "report_generation/generate_markdown_report",
            "value": 131187,
            "range": "± 13304",
            "unit": "ns/iter"
          },
          {
            "name": "report_generation/generate_csv_report",
            "value": 135695,
            "range": "± 10726",
            "unit": "ns/iter"
          },
          {
            "name": "dataset_operations/create_dataset_info",
            "value": 75,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "dataset_operations/validate_dataset_info",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "comparison_operations/create_comparison",
            "value": 1,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "comparison_operations/check_targets",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_gc_content",
            "value": 719324,
            "range": "± 9555",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_revcomp",
            "value": 258383,
            "range": "± 3201",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_subsequence",
            "value": 147,
            "range": "± 6",
            "unit": "ns/iter"
          },
          {
            "name": "memory_pool_allocation",
            "value": 70,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "cache_put_get",
            "value": 107434,
            "range": "± 1821",
            "unit": "ns/iter"
          },
          {
            "name": "interval_query",
            "value": 41532845,
            "range": "± 935041",
            "unit": "ns/iter"
          },
          {
            "name": "variant_type",
            "value": 16813,
            "range": "± 259",
            "unit": "ns/iter"
          },
          {
            "name": "quality_error_prob",
            "value": 1281,
            "range": "± 7",
            "unit": "ns/iter"
          },
          {
            "name": "log_addition",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "stream_processing",
            "value": 564,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "memory_mapped_access",
            "value": 8741,
            "range": "± 163",
            "unit": "ns/iter"
          },
          {
            "name": "hamming_distance",
            "value": 8,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/100",
            "value": 31938,
            "range": "± 217",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/100",
            "value": 53764,
            "range": "± 1357",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/100",
            "value": 27048,
            "range": "± 543",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/1000",
            "value": 276227,
            "range": "± 1398",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/1000",
            "value": 397855,
            "range": "± 2242",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/1000",
            "value": 228249,
            "range": "± 2781",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/10000",
            "value": 2753125,
            "range": "± 38339",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/10000",
            "value": 3732882,
            "range": "± 32971",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/10000",
            "value": 2290576,
            "range": "± 60950",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/100",
            "value": 33629,
            "range": "± 205",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/100",
            "value": 46396,
            "range": "± 231",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/100",
            "value": 24647,
            "range": "± 137",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/1000",
            "value": 331783,
            "range": "± 3163",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/1000",
            "value": 312935,
            "range": "± 1570",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/1000",
            "value": 262298,
            "range": "± 1297",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/10000",
            "value": 3182159,
            "range": "± 36176",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/10000",
            "value": 2865212,
            "range": "± 38698",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/10000",
            "value": 2615989,
            "range": "± 40028",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/gc_content",
            "value": 57,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/reverse_complement",
            "value": 55,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/subsequence",
            "value": 20,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/count_pattern",
            "value": 209,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/average_quality",
            "value": 54,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/min_quality",
            "value": 40,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/max_quality",
            "value": 45,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/median_quality",
            "value": 42,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/trim_quality",
            "value": 51,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/reverse_complement",
            "value": 89,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/is_valid",
            "value": 83,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/build_index",
            "value": 364151,
            "range": "± 6130",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/random_access",
            "value": 4391,
            "range": "± 13",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/save_index",
            "value": 305065,
            "range": "± 18981",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/load_index",
            "value": 292809,
            "range": "± 4887",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/filter_by_quality",
            "value": 387952,
            "range": "± 3289",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/filter_by_length",
            "value": 245848,
            "range": "± 1825",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/sample_reads",
            "value": 241576,
            "range": "± 19814",
            "unit": "ns/iter"
          },
          {
            "name": "io_writing/write_fasta",
            "value": 378897,
            "range": "± 22336",
            "unit": "ns/iter"
          },
          {
            "name": "io_writing/write_fastq",
            "value": 8356388,
            "range": "± 187131",
            "unit": "ns/iter"
          },
          {
            "name": "memory_usage/parse_memory_usage",
            "value": 2789690,
            "range": "± 324758",
            "unit": "ns/iter"
          },
          {
            "name": "macro_sum_bench",
            "value": 1,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/100",
            "value": 92952,
            "range": "± 699",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/1000",
            "value": 927058,
            "range": "± 6073",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/10000",
            "value": 9977760,
            "range": "± 272689",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_iterator/iterator",
            "value": 621373,
            "range": "± 6475",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_writing/writing",
            "value": 3595365,
            "range": "± 76307",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/0",
            "value": 99,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/1",
            "value": 97,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/2",
            "value": 94,
            "range": "± 4",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/3",
            "value": 120,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/4",
            "value": 112,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/5",
            "value": 58,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/6",
            "value": 115,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/7",
            "value": 491,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/0",
            "value": 85,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/1",
            "value": 85,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/2",
            "value": 42,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/3",
            "value": 85,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/4",
            "value": 85,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/5",
            "value": 85,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/6",
            "value": 85,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/7",
            "value": 85,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/variant_type_detection",
            "value": 2506,
            "range": "± 12",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/allele_frequency_access",
            "value": 3855,
            "range": "± 10",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/sample_data_access",
            "value": 1377,
            "range": "± 8",
            "unit": "ns/iter"
          },
          {
            "name": "band_pass_normalized_kernel_gatk_defaults",
            "value": 987,
            "range": "± 4",
            "unit": "ns/iter"
          },
          {
            "name": "band_pass_add_then_pop_ready_256_loci",
            "value": 299779,
            "range": "± 3016",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_medium_depth_k10",
            "value": 235805,
            "range": "± 1380",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_high_depth_k10",
            "value": 1632288,
            "range": "± 5679",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/4",
            "value": 13666849,
            "range": "± 132610",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/16",
            "value": 48065252,
            "range": "± 1035004",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/64",
            "value": 185468398,
            "range": "± 16091481",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/8",
            "value": 498767,
            "range": "± 1922",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/8",
            "value": 214640,
            "range": "± 1247",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/8",
            "value": 1079962,
            "range": "± 24726",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/8",
            "value": 9476518,
            "range": "± 27802",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/32",
            "value": 1994208,
            "range": "± 2808",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/32",
            "value": 853172,
            "range": "± 15460",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/32",
            "value": 4292349,
            "range": "± 24924",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/32",
            "value": 37916722,
            "range": "± 1596029",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/64",
            "value": 3985586,
            "range": "± 5205",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/64",
            "value": 1702511,
            "range": "± 10873",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/64",
            "value": 8566904,
            "range": "± 323397",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/64",
            "value": 75820711,
            "range": "± 152809",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/128",
            "value": 7966264,
            "range": "± 19120",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/128",
            "value": 3172309,
            "range": "± 19627",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/128",
            "value": 17141065,
            "range": "± 20815",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/128",
            "value": 151750133,
            "range": "± 342795",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/8",
            "value": 1936636,
            "range": "± 184728",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/8",
            "value": 4642436,
            "range": "± 11602",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/8",
            "value": 4309691,
            "range": "± 8476",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/8",
            "value": 38360863,
            "range": "± 173171",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/32",
            "value": 7768530,
            "range": "± 347437",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/32",
            "value": 18337471,
            "range": "± 185228",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/32",
            "value": 17139398,
            "range": "± 46967",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/32",
            "value": 153439797,
            "range": "± 1863896",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/64",
            "value": 15768443,
            "range": "± 34313",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/64",
            "value": 37350096,
            "range": "± 982120",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/64",
            "value": 34153848,
            "range": "± 114530",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/64",
            "value": 306733098,
            "range": "± 827356",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/128",
            "value": 31486233,
            "range": "± 77094",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/128",
            "value": 75020413,
            "range": "± 354200",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/128",
            "value": 68645896,
            "range": "± 112733",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/128",
            "value": 613473979,
            "range": "± 1861426",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/8",
            "value": 4384966,
            "range": "± 8063",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/8",
            "value": 7538326,
            "range": "± 36066",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/8",
            "value": 9523531,
            "range": "± 19776",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/8",
            "value": 86507279,
            "range": "± 223197",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/32",
            "value": 17543978,
            "range": "± 25704",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/32",
            "value": 36616966,
            "range": "± 291491",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/32",
            "value": 37835669,
            "range": "± 79703",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/32",
            "value": 346347147,
            "range": "± 2900184",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/64",
            "value": 34851053,
            "range": "± 51463",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/64",
            "value": 71850090,
            "range": "± 2326237",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/64",
            "value": 76096131,
            "range": "± 2531458",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/64",
            "value": 693256268,
            "range": "± 7737122",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/128",
            "value": 69734692,
            "range": "± 344619",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/128",
            "value": 147748029,
            "range": "± 1924137",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/128",
            "value": 152402372,
            "range": "± 206937",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/128",
            "value": 1388536056,
            "range": "± 5518454",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/value_copy_stages/3072x12",
            "value": 4889,
            "range": "± 16",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/arc_share_stages/3072x12",
            "value": 3469,
            "range": "± 9",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/assembly_result_set_arc_fanout/3072",
            "value": 53,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/assembly_result_set_value_fanout/3072",
            "value": 1224,
            "range": "± 9",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/64x48",
            "value": 19302,
            "range": "± 484",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/64x48",
            "value": 19131,
            "range": "± 574",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/128x96",
            "value": 69685,
            "range": "± 1704",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/128x96",
            "value": 68886,
            "range": "± 1882",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/256x192",
            "value": 270323,
            "range": "± 2741",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/256x192",
            "value": 268874,
            "range": "± 613",
            "unit": "ns/iter"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "name": "SynapticFour",
            "username": "SynapticFour",
            "email": "contact@synapticfour.com"
          },
          "committer": {
            "name": "GitHub",
            "username": "web-flow",
            "email": "noreply@github.com"
          },
          "id": "4dce5e14b1273f502d811c2c84274295575201b2",
          "message": "perf(hc): fail-closed DP caps and ownership fixes for 100kb RSS (#84)\n\nClip finalize in place, realign via SharedBam COW, stream SeqGraphs, share\nBAM header when sequential, and refuse PairHMM/SW above 8M cells so dense\nwindows no longer climb toward multi-GiB Peak-RSS on 16 GiB hosts.\n\nCo-authored-by: Cursor <cursoragent@cursor.com>",
          "timestamp": "2026-08-06T04:26:46Z",
          "url": "https://github.com/SynapticFour/gatk-rs/commit/4dce5e14b1273f502d811c2c84274295575201b2"
        },
        "date": 1786076662080,
        "tool": "cargo",
        "benches": [
          {
            "name": "sam_parsing/parsing/100",
            "value": 107990,
            "range": "± 1363",
            "unit": "ns/iter"
          },
          {
            "name": "sam_parsing/parsing/1000",
            "value": 1080902,
            "range": "± 20283",
            "unit": "ns/iter"
          },
          {
            "name": "sam_parsing/parsing/10000",
            "value": 10600350,
            "range": "± 177600",
            "unit": "ns/iter"
          },
          {
            "name": "sam_iterator/iterator",
            "value": 887710,
            "range": "± 6326",
            "unit": "ns/iter"
          },
          {
            "name": "sam_writing/writing",
            "value": 3259660,
            "range": "± 42606",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/0",
            "value": 125,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/1",
            "value": 176,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/2",
            "value": 139,
            "range": "± 3",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/3",
            "value": 143,
            "range": "± 8",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/4",
            "value": 168,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "optional_fields/parsing",
            "value": 7412,
            "range": "± 15",
            "unit": "ns/iter"
          },
          {
            "name": "framework_overhead/create_result",
            "value": 26,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "framework_overhead/create_comparison",
            "value": 1,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "metrics_collection/metrics_collector",
            "value": 1064373,
            "range": "± 1409",
            "unit": "ns/iter"
          },
          {
            "name": "suite_operations/create_suite",
            "value": 16409360,
            "range": "± 167798",
            "unit": "ns/iter"
          },
          {
            "name": "suite_operations/suite_with_results",
            "value": 16338046,
            "range": "± 193012",
            "unit": "ns/iter"
          },
          {
            "name": "suite_operations/suite_serialization",
            "value": 74969,
            "range": "± 1575",
            "unit": "ns/iter"
          },
          {
            "name": "performance_analysis/analyze_small_dataset",
            "value": 84,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "performance_analysis/analyze_large_dataset",
            "value": 8860,
            "range": "± 18",
            "unit": "ns/iter"
          },
          {
            "name": "regression_detection/detect_no_regression",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "regression_detection/detect_regression",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "report_generation/generate_json_report",
            "value": 137813,
            "range": "± 11939",
            "unit": "ns/iter"
          },
          {
            "name": "report_generation/generate_markdown_report",
            "value": 118695,
            "range": "± 6976",
            "unit": "ns/iter"
          },
          {
            "name": "report_generation/generate_csv_report",
            "value": 124547,
            "range": "± 15483",
            "unit": "ns/iter"
          },
          {
            "name": "dataset_operations/create_dataset_info",
            "value": 81,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "dataset_operations/validate_dataset_info",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "comparison_operations/create_comparison",
            "value": 1,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "comparison_operations/check_targets",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_gc_content",
            "value": 814458,
            "range": "± 10460",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_revcomp",
            "value": 297325,
            "range": "± 6458",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_subsequence",
            "value": 171,
            "range": "± 6",
            "unit": "ns/iter"
          },
          {
            "name": "memory_pool_allocation",
            "value": 72,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "cache_put_get",
            "value": 110401,
            "range": "± 573",
            "unit": "ns/iter"
          },
          {
            "name": "interval_query",
            "value": 40968183,
            "range": "± 1247927",
            "unit": "ns/iter"
          },
          {
            "name": "variant_type",
            "value": 16203,
            "range": "± 52",
            "unit": "ns/iter"
          },
          {
            "name": "quality_error_prob",
            "value": 1214,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "log_addition",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "stream_processing",
            "value": 666,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "memory_mapped_access",
            "value": 9699,
            "range": "± 190",
            "unit": "ns/iter"
          },
          {
            "name": "hamming_distance",
            "value": 9,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/100",
            "value": 29453,
            "range": "± 206",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/100",
            "value": 52014,
            "range": "± 1582",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/100",
            "value": 26105,
            "range": "± 262",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/1000",
            "value": 243227,
            "range": "± 2965",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/1000",
            "value": 382526,
            "range": "± 5766",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/1000",
            "value": 212403,
            "range": "± 924",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/10000",
            "value": 2442204,
            "range": "± 27076",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/10000",
            "value": 3586732,
            "range": "± 13513",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/10000",
            "value": 2096024,
            "range": "± 10628",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/100",
            "value": 33737,
            "range": "± 65",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/100",
            "value": 43648,
            "range": "± 216",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/100",
            "value": 26313,
            "range": "± 110",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/1000",
            "value": 273680,
            "range": "± 1668",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/1000",
            "value": 303015,
            "range": "± 2244",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/1000",
            "value": 209854,
            "range": "± 1090",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/10000",
            "value": 2634247,
            "range": "± 15602",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/10000",
            "value": 2824528,
            "range": "± 26561",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/10000",
            "value": 2049188,
            "range": "± 18433",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/gc_content",
            "value": 63,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/reverse_complement",
            "value": 60,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/subsequence",
            "value": 24,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/count_pattern",
            "value": 227,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/average_quality",
            "value": 62,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/min_quality",
            "value": 45,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/max_quality",
            "value": 49,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/median_quality",
            "value": 45,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/trim_quality",
            "value": 54,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/reverse_complement",
            "value": 93,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/is_valid",
            "value": 92,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/build_index",
            "value": 349080,
            "range": "± 6232",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/random_access",
            "value": 4864,
            "range": "± 71",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/save_index",
            "value": 220132,
            "range": "± 15821",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/load_index",
            "value": 299497,
            "range": "± 935",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/filter_by_quality",
            "value": 334567,
            "range": "± 788",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/filter_by_length",
            "value": 190278,
            "range": "± 630",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/sample_reads",
            "value": 191922,
            "range": "± 1638",
            "unit": "ns/iter"
          },
          {
            "name": "io_writing/write_fasta",
            "value": 266041,
            "range": "± 7030",
            "unit": "ns/iter"
          },
          {
            "name": "io_writing/write_fastq",
            "value": 8675175,
            "range": "± 572967",
            "unit": "ns/iter"
          },
          {
            "name": "memory_usage/parse_memory_usage",
            "value": 2470277,
            "range": "± 62419",
            "unit": "ns/iter"
          },
          {
            "name": "macro_sum_bench",
            "value": 1,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/100",
            "value": 89556,
            "range": "± 7094",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/1000",
            "value": 866491,
            "range": "± 12204",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/10000",
            "value": 8788398,
            "range": "± 165090",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_iterator/iterator",
            "value": 611800,
            "range": "± 4846",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_writing/writing",
            "value": 3810734,
            "range": "± 74606",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/0",
            "value": 107,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/1",
            "value": 96,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/2",
            "value": 93,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/3",
            "value": 131,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/4",
            "value": 111,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/5",
            "value": 59,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/6",
            "value": 114,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/7",
            "value": 478,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/0",
            "value": 91,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/1",
            "value": 91,
            "range": "± 3",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/2",
            "value": 47,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/3",
            "value": 91,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/4",
            "value": 91,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/5",
            "value": 91,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/6",
            "value": 91,
            "range": "± 5",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/7",
            "value": 92,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/variant_type_detection",
            "value": 2475,
            "range": "± 4",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/allele_frequency_access",
            "value": 4249,
            "range": "± 9",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/sample_data_access",
            "value": 1522,
            "range": "± 15",
            "unit": "ns/iter"
          },
          {
            "name": "band_pass_normalized_kernel_gatk_defaults",
            "value": 995,
            "range": "± 51",
            "unit": "ns/iter"
          },
          {
            "name": "band_pass_add_then_pop_ready_256_loci",
            "value": 324030,
            "range": "± 1699",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_medium_depth_k10",
            "value": 244194,
            "range": "± 213",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_high_depth_k10",
            "value": 1675903,
            "range": "± 21033",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/4",
            "value": 13550207,
            "range": "± 259275",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/16",
            "value": 47736456,
            "range": "± 1118233",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/64",
            "value": 139348298,
            "range": "± 19803445",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/8",
            "value": 543523,
            "range": "± 9085",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/8",
            "value": 223977,
            "range": "± 3380",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/8",
            "value": 1092779,
            "range": "± 10487",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/8",
            "value": 9037067,
            "range": "± 28556",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/32",
            "value": 2167750,
            "range": "± 4405",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/32",
            "value": 837420,
            "range": "± 12685",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/32",
            "value": 4338720,
            "range": "± 7169",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/32",
            "value": 36028958,
            "range": "± 42330",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/64",
            "value": 4334765,
            "range": "± 99406",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/64",
            "value": 1666884,
            "range": "± 6974",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/64",
            "value": 8675003,
            "range": "± 10130",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/64",
            "value": 72100190,
            "range": "± 1545854",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/128",
            "value": 8661666,
            "range": "± 159012",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/128",
            "value": 3324362,
            "range": "± 10634",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/128",
            "value": 17334439,
            "range": "± 14932",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/128",
            "value": 144146995,
            "range": "± 1210959",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/8",
            "value": 2124949,
            "range": "± 3101",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/8",
            "value": 4873465,
            "range": "± 12895",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/8",
            "value": 4333737,
            "range": "± 3799",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/8",
            "value": 36490782,
            "range": "± 76823",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/32",
            "value": 8498449,
            "range": "± 177122",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/32",
            "value": 19586115,
            "range": "± 34697",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/32",
            "value": 17256772,
            "range": "± 30827",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/32",
            "value": 145722328,
            "range": "± 1587476",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/64",
            "value": 16990952,
            "range": "± 13908",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/64",
            "value": 38976416,
            "range": "± 61421",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/64",
            "value": 34493769,
            "range": "± 67398",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/64",
            "value": 291534375,
            "range": "± 2498207",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/128",
            "value": 34007958,
            "range": "± 72475",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/128",
            "value": 77964928,
            "range": "± 1717008",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/128",
            "value": 68933068,
            "range": "± 941415",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/128",
            "value": 582889287,
            "range": "± 2200414",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/8",
            "value": 4716248,
            "range": "± 6191",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/8",
            "value": 7818577,
            "range": "± 17658",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/8",
            "value": 9658938,
            "range": "± 46863",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/8",
            "value": 82451237,
            "range": "± 1364254",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/32",
            "value": 18859399,
            "range": "± 31979",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/32",
            "value": 25013754,
            "range": "± 89468",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/32",
            "value": 38503698,
            "range": "± 33787",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/32",
            "value": 328966930,
            "range": "± 2590248",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/64",
            "value": 37928250,
            "range": "± 230072",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/64",
            "value": 50189570,
            "range": "± 136016",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/64",
            "value": 77332826,
            "range": "± 410485",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/64",
            "value": 657798612,
            "range": "± 3296773",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/128",
            "value": 75822784,
            "range": "± 97503",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/128",
            "value": 100127966,
            "range": "± 1255298",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/128",
            "value": 154597279,
            "range": "± 122163",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/128",
            "value": 1316396274,
            "range": "± 3441061",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/value_copy_stages/3072x12",
            "value": 5146,
            "range": "± 110",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/arc_share_stages/3072x12",
            "value": 3930,
            "range": "± 6",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/assembly_result_set_arc_fanout/3072",
            "value": 56,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/assembly_result_set_value_fanout/3072",
            "value": 1163,
            "range": "± 3",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/64x48",
            "value": 17881,
            "range": "± 149",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/64x48",
            "value": 17934,
            "range": "± 616",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/128x96",
            "value": 69702,
            "range": "± 181",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/128x96",
            "value": 69057,
            "range": "± 200",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/256x192",
            "value": 268252,
            "range": "± 2076",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/256x192",
            "value": 265226,
            "range": "± 581",
            "unit": "ns/iter"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "name": "SynapticFour",
            "username": "SynapticFour",
            "email": "contact@synapticfour.com"
          },
          "committer": {
            "name": "GitHub",
            "username": "web-flow",
            "email": "noreply@github.com"
          },
          "id": "4dce5e14b1273f502d811c2c84274295575201b2",
          "message": "perf(hc): fail-closed DP caps and ownership fixes for 100kb RSS (#84)\n\nClip finalize in place, realign via SharedBam COW, stream SeqGraphs, share\nBAM header when sequential, and refuse PairHMM/SW above 8M cells so dense\nwindows no longer climb toward multi-GiB Peak-RSS on 16 GiB hosts.\n\nCo-authored-by: Cursor <cursoragent@cursor.com>",
          "timestamp": "2026-08-06T04:26:46Z",
          "url": "https://github.com/SynapticFour/gatk-rs/commit/4dce5e14b1273f502d811c2c84274295575201b2"
        },
        "date": 1786159507672,
        "tool": "cargo",
        "benches": [
          {
            "name": "sam_parsing/parsing/100",
            "value": 84901,
            "range": "± 1280",
            "unit": "ns/iter"
          },
          {
            "name": "sam_parsing/parsing/1000",
            "value": 879522,
            "range": "± 19523",
            "unit": "ns/iter"
          },
          {
            "name": "sam_parsing/parsing/10000",
            "value": 8351888,
            "range": "± 54084",
            "unit": "ns/iter"
          },
          {
            "name": "sam_iterator/iterator",
            "value": 637280,
            "range": "± 13103",
            "unit": "ns/iter"
          },
          {
            "name": "sam_writing/writing",
            "value": 2603196,
            "range": "± 89448",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/0",
            "value": 102,
            "range": "± 6",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/1",
            "value": 137,
            "range": "± 7",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/2",
            "value": 111,
            "range": "± 6",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/3",
            "value": 114,
            "range": "± 5",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/4",
            "value": 136,
            "range": "± 5",
            "unit": "ns/iter"
          },
          {
            "name": "optional_fields/parsing",
            "value": 5596,
            "range": "± 21",
            "unit": "ns/iter"
          },
          {
            "name": "framework_overhead/create_result",
            "value": 21,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "framework_overhead/create_comparison",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "metrics_collection/metrics_collector",
            "value": 1059826,
            "range": "± 920",
            "unit": "ns/iter"
          },
          {
            "name": "suite_operations/create_suite",
            "value": 13415303,
            "range": "± 212715",
            "unit": "ns/iter"
          },
          {
            "name": "suite_operations/suite_with_results",
            "value": 13693208,
            "range": "± 327502",
            "unit": "ns/iter"
          },
          {
            "name": "suite_operations/suite_serialization",
            "value": 59203,
            "range": "± 345",
            "unit": "ns/iter"
          },
          {
            "name": "performance_analysis/analyze_small_dataset",
            "value": 65,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "performance_analysis/analyze_large_dataset",
            "value": 7456,
            "range": "± 11",
            "unit": "ns/iter"
          },
          {
            "name": "regression_detection/detect_no_regression",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "regression_detection/detect_regression",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "report_generation/generate_json_report",
            "value": 165214,
            "range": "± 290883",
            "unit": "ns/iter"
          },
          {
            "name": "report_generation/generate_markdown_report",
            "value": 99144,
            "range": "± 43465",
            "unit": "ns/iter"
          },
          {
            "name": "report_generation/generate_csv_report",
            "value": 98011,
            "range": "± 76124",
            "unit": "ns/iter"
          },
          {
            "name": "dataset_operations/create_dataset_info",
            "value": 63,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "dataset_operations/validate_dataset_info",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "comparison_operations/create_comparison",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "comparison_operations/check_targets",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_gc_content",
            "value": 638157,
            "range": "± 20587",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_revcomp",
            "value": 231448,
            "range": "± 7471",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_subsequence",
            "value": 128,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "memory_pool_allocation",
            "value": 58,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "cache_put_get",
            "value": 84490,
            "range": "± 2086",
            "unit": "ns/iter"
          },
          {
            "name": "interval_query",
            "value": 32223598,
            "range": "± 1198568",
            "unit": "ns/iter"
          },
          {
            "name": "variant_type",
            "value": 12564,
            "range": "± 13",
            "unit": "ns/iter"
          },
          {
            "name": "quality_error_prob",
            "value": 954,
            "range": "± 23",
            "unit": "ns/iter"
          },
          {
            "name": "log_addition",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "stream_processing",
            "value": 442,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "memory_mapped_access",
            "value": 7151,
            "range": "± 282",
            "unit": "ns/iter"
          },
          {
            "name": "hamming_distance",
            "value": 7,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/100",
            "value": 22882,
            "range": "± 684",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/100",
            "value": 40614,
            "range": "± 919",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/100",
            "value": 20210,
            "range": "± 124",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/1000",
            "value": 189332,
            "range": "± 3663",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/1000",
            "value": 297943,
            "range": "± 839",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/1000",
            "value": 162572,
            "range": "± 458",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/10000",
            "value": 1914455,
            "range": "± 12074",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/10000",
            "value": 2764835,
            "range": "± 8219",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/10000",
            "value": 1630976,
            "range": "± 5407",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/100",
            "value": 25150,
            "range": "± 169",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/100",
            "value": 35041,
            "range": "± 644",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/100",
            "value": 20015,
            "range": "± 81",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/1000",
            "value": 207684,
            "range": "± 1448",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/1000",
            "value": 241110,
            "range": "± 1040",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/1000",
            "value": 160641,
            "range": "± 650",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/10000",
            "value": 2037896,
            "range": "± 11963",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/10000",
            "value": 2259374,
            "range": "± 8005",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/10000",
            "value": 1594069,
            "range": "± 8130",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/gc_content",
            "value": 49,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/reverse_complement",
            "value": 48,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/subsequence",
            "value": 18,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/count_pattern",
            "value": 161,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/average_quality",
            "value": 48,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/min_quality",
            "value": 34,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/max_quality",
            "value": 38,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/median_quality",
            "value": 34,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/trim_quality",
            "value": 41,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/reverse_complement",
            "value": 71,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/is_valid",
            "value": 71,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/build_index",
            "value": 263414,
            "range": "± 3146",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/random_access",
            "value": 3764,
            "range": "± 10",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/save_index",
            "value": 517791,
            "range": "± 121720",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/load_index",
            "value": 233457,
            "range": "± 2978",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/filter_by_quality",
            "value": 256932,
            "range": "± 4654",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/filter_by_length",
            "value": 146245,
            "range": "± 1779",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/sample_reads",
            "value": 146538,
            "range": "± 377",
            "unit": "ns/iter"
          },
          {
            "name": "io_writing/write_fasta",
            "value": 502302,
            "range": "± 236599",
            "unit": "ns/iter"
          },
          {
            "name": "io_writing/write_fastq",
            "value": 6852971,
            "range": "± 2007298",
            "unit": "ns/iter"
          },
          {
            "name": "memory_usage/parse_memory_usage",
            "value": 1922623,
            "range": "± 29561",
            "unit": "ns/iter"
          },
          {
            "name": "macro_sum_bench",
            "value": 1,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/100",
            "value": 69992,
            "range": "± 515",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/1000",
            "value": 679133,
            "range": "± 31791",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/10000",
            "value": 7214293,
            "range": "± 67753",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_iterator/iterator",
            "value": 474902,
            "range": "± 4305",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_writing/writing",
            "value": 3098282,
            "range": "± 123972",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/0",
            "value": 82,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/1",
            "value": 75,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/2",
            "value": 72,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/3",
            "value": 101,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/4",
            "value": 86,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/5",
            "value": 46,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/6",
            "value": 88,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/7",
            "value": 358,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/0",
            "value": 73,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/1",
            "value": 73,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/2",
            "value": 36,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/3",
            "value": 73,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/4",
            "value": 73,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/5",
            "value": 73,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/6",
            "value": 73,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/7",
            "value": 73,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/variant_type_detection",
            "value": 1919,
            "range": "± 32",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/allele_frequency_access",
            "value": 3281,
            "range": "± 7",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/sample_data_access",
            "value": 1181,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "band_pass_normalized_kernel_gatk_defaults",
            "value": 769,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "band_pass_add_then_pop_ready_256_loci",
            "value": 245998,
            "range": "± 1630",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_medium_depth_k10",
            "value": 187407,
            "range": "± 129",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_high_depth_k10",
            "value": 1293371,
            "range": "± 15164",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/4",
            "value": 10477472,
            "range": "± 187661",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/16",
            "value": 37146768,
            "range": "± 695523",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/64",
            "value": 108622081,
            "range": "± 13572651",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/8",
            "value": 421754,
            "range": "± 413",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/8",
            "value": 166716,
            "range": "± 1496",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/8",
            "value": 854191,
            "range": "± 1550",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/8",
            "value": 6996901,
            "range": "± 6317",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/32",
            "value": 1686280,
            "range": "± 4897",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/32",
            "value": 675170,
            "range": "± 8470",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/32",
            "value": 3386748,
            "range": "± 5255",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/32",
            "value": 27955266,
            "range": "± 21598",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/64",
            "value": 3371346,
            "range": "± 2165",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/64",
            "value": 1342760,
            "range": "± 7660",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/64",
            "value": 6830974,
            "range": "± 48270",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/64",
            "value": 55977756,
            "range": "± 73710",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/128",
            "value": 6749873,
            "range": "± 13645",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/128",
            "value": 2697352,
            "range": "± 7479",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/128",
            "value": 13649445,
            "range": "± 97332",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/128",
            "value": 111879182,
            "range": "± 88920",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/8",
            "value": 1653415,
            "range": "± 29932",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/8",
            "value": 3828405,
            "range": "± 17134",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/8",
            "value": 3366764,
            "range": "± 4025",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/8",
            "value": 28293197,
            "range": "± 372954",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/32",
            "value": 6612704,
            "range": "± 6535",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/32",
            "value": 15760988,
            "range": "± 217998",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/32",
            "value": 13483911,
            "range": "± 43029",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/32",
            "value": 113220802,
            "range": "± 479300",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/64",
            "value": 13223346,
            "range": "± 145695",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/64",
            "value": 30740705,
            "range": "± 1284670",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/64",
            "value": 26919031,
            "range": "± 102564",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/64",
            "value": 226223176,
            "range": "± 127422",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/128",
            "value": 26464876,
            "range": "± 174559",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/128",
            "value": 60745395,
            "range": "± 355884",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/128",
            "value": 53866091,
            "range": "± 208577",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/128",
            "value": 452782697,
            "range": "± 2576455",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/8",
            "value": 3667301,
            "range": "± 3326",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/8",
            "value": 6496773,
            "range": "± 158411",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/8",
            "value": 7534558,
            "range": "± 372298",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/8",
            "value": 63853199,
            "range": "± 82886",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/32",
            "value": 14667246,
            "range": "± 21301",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/32",
            "value": 25849008,
            "range": "± 2249551",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/32",
            "value": 30017862,
            "range": "± 72199",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/32",
            "value": 255550587,
            "range": "± 472509",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/64",
            "value": 29385043,
            "range": "± 31823",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/64",
            "value": 51821771,
            "range": "± 883662",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/64",
            "value": 60241681,
            "range": "± 131722",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/64",
            "value": 510996564,
            "range": "± 2421826",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/128",
            "value": 59082150,
            "range": "± 51357",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/128",
            "value": 103594836,
            "range": "± 1171450",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/128",
            "value": 120433230,
            "range": "± 363072",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/128",
            "value": 1022088511,
            "range": "± 958035",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/value_copy_stages/3072x12",
            "value": 3947,
            "range": "± 5",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/arc_share_stages/3072x12",
            "value": 3045,
            "range": "± 21",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/assembly_result_set_arc_fanout/3072",
            "value": 43,
            "range": "± 3",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/assembly_result_set_value_fanout/3072",
            "value": 910,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/64x48",
            "value": 13732,
            "range": "± 41",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/64x48",
            "value": 13574,
            "range": "± 360",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/128x96",
            "value": 54173,
            "range": "± 1068",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/128x96",
            "value": 53583,
            "range": "± 103",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/256x192",
            "value": 212479,
            "range": "± 634",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/256x192",
            "value": 208931,
            "range": "± 1061",
            "unit": "ns/iter"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "name": "SynapticFour",
            "username": "SynapticFour",
            "email": "contact@synapticfour.com"
          },
          "committer": {
            "name": "GitHub",
            "username": "web-flow",
            "email": "noreply@github.com"
          },
          "id": "4dce5e14b1273f502d811c2c84274295575201b2",
          "message": "perf(hc): fail-closed DP caps and ownership fixes for 100kb RSS (#84)\n\nClip finalize in place, realign via SharedBam COW, stream SeqGraphs, share\nBAM header when sequential, and refuse PairHMM/SW above 8M cells so dense\nwindows no longer climb toward multi-GiB Peak-RSS on 16 GiB hosts.\n\nCo-authored-by: Cursor <cursoragent@cursor.com>",
          "timestamp": "2026-08-06T04:26:46Z",
          "url": "https://github.com/SynapticFour/gatk-rs/commit/4dce5e14b1273f502d811c2c84274295575201b2"
        },
        "date": 1786246156469,
        "tool": "cargo",
        "benches": [
          {
            "name": "sam_parsing/parsing/100",
            "value": 83407,
            "range": "± 630",
            "unit": "ns/iter"
          },
          {
            "name": "sam_parsing/parsing/1000",
            "value": 912721,
            "range": "± 5974",
            "unit": "ns/iter"
          },
          {
            "name": "sam_parsing/parsing/10000",
            "value": 9046846,
            "range": "± 44439",
            "unit": "ns/iter"
          },
          {
            "name": "sam_iterator/iterator",
            "value": 638574,
            "range": "± 1473",
            "unit": "ns/iter"
          },
          {
            "name": "sam_writing/writing",
            "value": 1619138,
            "range": "± 25758",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/0",
            "value": 106,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/1",
            "value": 154,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/2",
            "value": 118,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/3",
            "value": 123,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/4",
            "value": 155,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "optional_fields/parsing",
            "value": 3264,
            "range": "± 11",
            "unit": "ns/iter"
          },
          {
            "name": "framework_overhead/create_result",
            "value": 20,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "framework_overhead/create_comparison",
            "value": 0,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "metrics_collection/metrics_collector",
            "value": 1057269,
            "range": "± 504",
            "unit": "ns/iter"
          },
          {
            "name": "suite_operations/create_suite",
            "value": 9692428,
            "range": "± 96751",
            "unit": "ns/iter"
          },
          {
            "name": "suite_operations/suite_with_results",
            "value": 9705037,
            "range": "± 74325",
            "unit": "ns/iter"
          },
          {
            "name": "suite_operations/suite_serialization",
            "value": 52807,
            "range": "± 409",
            "unit": "ns/iter"
          },
          {
            "name": "performance_analysis/analyze_small_dataset",
            "value": 59,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "performance_analysis/analyze_large_dataset",
            "value": 7532,
            "range": "± 15",
            "unit": "ns/iter"
          },
          {
            "name": "regression_detection/detect_no_regression",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "regression_detection/detect_regression",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "report_generation/generate_json_report",
            "value": 161188,
            "range": "± 8930",
            "unit": "ns/iter"
          },
          {
            "name": "report_generation/generate_markdown_report",
            "value": 114775,
            "range": "± 15314",
            "unit": "ns/iter"
          },
          {
            "name": "report_generation/generate_csv_report",
            "value": 110395,
            "range": "± 5844",
            "unit": "ns/iter"
          },
          {
            "name": "dataset_operations/create_dataset_info",
            "value": 59,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "dataset_operations/validate_dataset_info",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "comparison_operations/create_comparison",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "comparison_operations/check_targets",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_gc_content",
            "value": 688372,
            "range": "± 5525",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_revcomp",
            "value": 294731,
            "range": "± 950",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_subsequence",
            "value": 96,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "memory_pool_allocation",
            "value": 55,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "cache_put_get",
            "value": 72220,
            "range": "± 587",
            "unit": "ns/iter"
          },
          {
            "name": "interval_query",
            "value": 28299973,
            "range": "± 1468012",
            "unit": "ns/iter"
          },
          {
            "name": "variant_type",
            "value": 13282,
            "range": "± 45",
            "unit": "ns/iter"
          },
          {
            "name": "quality_error_prob",
            "value": 993,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "log_addition",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "stream_processing",
            "value": 347,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "memory_mapped_access",
            "value": 3140,
            "range": "± 54",
            "unit": "ns/iter"
          },
          {
            "name": "hamming_distance",
            "value": 6,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/100",
            "value": 22820,
            "range": "± 265",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/100",
            "value": 38710,
            "range": "± 1838",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/100",
            "value": 19749,
            "range": "± 69",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/1000",
            "value": 213811,
            "range": "± 783",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/1000",
            "value": 294043,
            "range": "± 871",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/1000",
            "value": 182380,
            "range": "± 838",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/10000",
            "value": 2173567,
            "range": "± 8226",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/10000",
            "value": 2952566,
            "range": "± 13535",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/10000",
            "value": 1818642,
            "range": "± 20337",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/100",
            "value": 24670,
            "range": "± 59",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/100",
            "value": 31741,
            "range": "± 250",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/100",
            "value": 17941,
            "range": "± 84",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/1000",
            "value": 218687,
            "range": "± 427",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/1000",
            "value": 233915,
            "range": "± 2117",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/1000",
            "value": 161665,
            "range": "± 272",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/10000",
            "value": 2242395,
            "range": "± 10767",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/10000",
            "value": 2337409,
            "range": "± 10089",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/10000",
            "value": 1591304,
            "range": "± 4992",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/gc_content",
            "value": 48,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/reverse_complement",
            "value": 42,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/subsequence",
            "value": 16,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/count_pattern",
            "value": 139,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/average_quality",
            "value": 24,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/min_quality",
            "value": 48,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/max_quality",
            "value": 59,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/median_quality",
            "value": 41,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/trim_quality",
            "value": 45,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/reverse_complement",
            "value": 92,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/is_valid",
            "value": 70,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/build_index",
            "value": 286889,
            "range": "± 1175",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/random_access",
            "value": 1766,
            "range": "± 5",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/save_index",
            "value": 514306,
            "range": "± 72629",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/load_index",
            "value": 206126,
            "range": "± 559",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/filter_by_quality",
            "value": 256853,
            "range": "± 2866",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/filter_by_length",
            "value": 156139,
            "range": "± 798",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/sample_reads",
            "value": 154993,
            "range": "± 811",
            "unit": "ns/iter"
          },
          {
            "name": "io_writing/write_fasta",
            "value": 502706,
            "range": "± 37634",
            "unit": "ns/iter"
          },
          {
            "name": "io_writing/write_fastq",
            "value": 3585920,
            "range": "± 30924",
            "unit": "ns/iter"
          },
          {
            "name": "memory_usage/parse_memory_usage",
            "value": 2176665,
            "range": "± 5270",
            "unit": "ns/iter"
          },
          {
            "name": "macro_sum_bench",
            "value": 1,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/100",
            "value": 69067,
            "range": "± 566",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/1000",
            "value": 771347,
            "range": "± 2418",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/10000",
            "value": 8204321,
            "range": "± 33812",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_iterator/iterator",
            "value": 505363,
            "range": "± 5084",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_writing/writing",
            "value": 1985472,
            "range": "± 24969",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/0",
            "value": 71,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/1",
            "value": 69,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/2",
            "value": 67,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/3",
            "value": 94,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/4",
            "value": 82,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/5",
            "value": 41,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/6",
            "value": 90,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/7",
            "value": 425,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/0",
            "value": 68,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/1",
            "value": 68,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/2",
            "value": 29,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/3",
            "value": 69,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/4",
            "value": 69,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/5",
            "value": 68,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/6",
            "value": 68,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/7",
            "value": 68,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/variant_type_detection",
            "value": 1971,
            "range": "± 8",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/allele_frequency_access",
            "value": 2545,
            "range": "± 19",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/sample_data_access",
            "value": 1004,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "band_pass_normalized_kernel_gatk_defaults",
            "value": 749,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "band_pass_add_then_pop_ready_256_loci",
            "value": 255556,
            "range": "± 1287",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_medium_depth_k10",
            "value": 160003,
            "range": "± 337",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_high_depth_k10",
            "value": 1129296,
            "range": "± 5589",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/4",
            "value": 10391601,
            "range": "± 359349",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/16",
            "value": 29075578,
            "range": "± 4598274",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/64",
            "value": 150606303,
            "range": "± 10706565",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/8",
            "value": 378434,
            "range": "± 1376",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/8",
            "value": 176321,
            "range": "± 855",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/8",
            "value": 765146,
            "range": "± 1128",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/8",
            "value": 7098056,
            "range": "± 21813",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/32",
            "value": 1509900,
            "range": "± 12613",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/32",
            "value": 680169,
            "range": "± 1827",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/32",
            "value": 3036815,
            "range": "± 6758",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/32",
            "value": 28403461,
            "range": "± 82415",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/64",
            "value": 3016937,
            "range": "± 6255",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/64",
            "value": 1354908,
            "range": "± 2225",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/64",
            "value": 6062145,
            "range": "± 13741",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/64",
            "value": 56835295,
            "range": "± 95978",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/128",
            "value": 6035980,
            "range": "± 6420",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/128",
            "value": 2697664,
            "range": "± 50845",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/128",
            "value": 12109659,
            "range": "± 19176",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/128",
            "value": 113665465,
            "range": "± 281215",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/8",
            "value": 1471765,
            "range": "± 18324",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/8",
            "value": 2946337,
            "range": "± 11003",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/8",
            "value": 3043653,
            "range": "± 15342",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/8",
            "value": 28621727,
            "range": "± 73310",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/32",
            "value": 5881195,
            "range": "± 22029",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/32",
            "value": 11865079,
            "range": "± 39333",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/32",
            "value": 12100777,
            "range": "± 22705",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/32",
            "value": 114476634,
            "range": "± 335357",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/64",
            "value": 11769132,
            "range": "± 21470",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/64",
            "value": 23585355,
            "range": "± 62438",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/64",
            "value": 24175823,
            "range": "± 30284",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/64",
            "value": 229034184,
            "range": "± 182135",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/128",
            "value": 23546679,
            "range": "± 173979",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/128",
            "value": 47476308,
            "range": "± 127451",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/128",
            "value": 48230904,
            "range": "± 71079",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/128",
            "value": 458219294,
            "range": "± 626981",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/8",
            "value": 3463534,
            "range": "± 4948",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/8",
            "value": 5040493,
            "range": "± 13279",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/8",
            "value": 7298283,
            "range": "± 8441",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/8",
            "value": 64749623,
            "range": "± 261259",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/32",
            "value": 13886248,
            "range": "± 37632",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/32",
            "value": 20266262,
            "range": "± 64167",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/32",
            "value": 29052371,
            "range": "± 53098",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/32",
            "value": 258680934,
            "range": "± 243941",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/64",
            "value": 27754865,
            "range": "± 35283",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/64",
            "value": 40760135,
            "range": "± 429215",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/64",
            "value": 58040012,
            "range": "± 62562",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/64",
            "value": 518183234,
            "range": "± 564190",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/128",
            "value": 55634946,
            "range": "± 483529",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/128",
            "value": 81289135,
            "range": "± 203494",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/128",
            "value": 116008108,
            "range": "± 124924",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/128",
            "value": 1036334432,
            "range": "± 519793",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/value_copy_stages/3072x12",
            "value": 3948,
            "range": "± 11",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/arc_share_stages/3072x12",
            "value": 2954,
            "range": "± 6",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/assembly_result_set_arc_fanout/3072",
            "value": 151,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/assembly_result_set_value_fanout/3072",
            "value": 1113,
            "range": "± 22",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/64x48",
            "value": 11620,
            "range": "± 35",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/64x48",
            "value": 11584,
            "range": "± 286",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/128x96",
            "value": 49003,
            "range": "± 1659",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/128x96",
            "value": 48928,
            "range": "± 146",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/256x192",
            "value": 188824,
            "range": "± 507",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/256x192",
            "value": 189524,
            "range": "± 476",
            "unit": "ns/iter"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "name": "SynapticFour",
            "username": "SynapticFour",
            "email": "contact@synapticfour.com"
          },
          "committer": {
            "name": "GitHub",
            "username": "web-flow",
            "email": "noreply@github.com"
          },
          "id": "87bea2944d2f5ae1f8e6a63ed33b5a27c5b7a034",
          "message": "fix(hc): declare INFO/FORMAT in non-gVCF HaplotypeCaller headers (#87)\n\nGIAB smoke hap.py vcfcheck rejected rust.vcf because body sites carried\nAC/AF/… and GT:GQ:DP:AD:PL without matching ##INFO/##FORMAT lines\n(contig-only header). Populate the schema used by region emit.\n\nCo-authored-by: Cursor <cursoragent@cursor.com>",
          "timestamp": "2026-08-10T03:35:30Z",
          "url": "https://github.com/SynapticFour/gatk-rs/commit/87bea2944d2f5ae1f8e6a63ed33b5a27c5b7a034"
        },
        "date": 1786334985196,
        "tool": "cargo",
        "benches": [
          {
            "name": "sam_parsing/parsing/100",
            "value": 108313,
            "range": "± 1759",
            "unit": "ns/iter"
          },
          {
            "name": "sam_parsing/parsing/1000",
            "value": 1088500,
            "range": "± 16486",
            "unit": "ns/iter"
          },
          {
            "name": "sam_parsing/parsing/10000",
            "value": 10632781,
            "range": "± 36250",
            "unit": "ns/iter"
          },
          {
            "name": "sam_iterator/iterator",
            "value": 847757,
            "range": "± 6104",
            "unit": "ns/iter"
          },
          {
            "name": "sam_writing/writing",
            "value": 3186375,
            "range": "± 225222",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/0",
            "value": 126,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/1",
            "value": 171,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/2",
            "value": 136,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/3",
            "value": 144,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/4",
            "value": 170,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "optional_fields/parsing",
            "value": 7188,
            "range": "± 25",
            "unit": "ns/iter"
          },
          {
            "name": "framework_overhead/create_result",
            "value": 27,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "framework_overhead/create_comparison",
            "value": 1,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "metrics_collection/metrics_collector",
            "value": 1062062,
            "range": "± 1278",
            "unit": "ns/iter"
          },
          {
            "name": "suite_operations/create_suite",
            "value": 17116974,
            "range": "± 223407",
            "unit": "ns/iter"
          },
          {
            "name": "suite_operations/suite_with_results",
            "value": 16784560,
            "range": "± 239888",
            "unit": "ns/iter"
          },
          {
            "name": "suite_operations/suite_serialization",
            "value": 73938,
            "range": "± 539",
            "unit": "ns/iter"
          },
          {
            "name": "performance_analysis/analyze_small_dataset",
            "value": 83,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "performance_analysis/analyze_large_dataset",
            "value": 8854,
            "range": "± 8",
            "unit": "ns/iter"
          },
          {
            "name": "regression_detection/detect_no_regression",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "regression_detection/detect_regression",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "report_generation/generate_json_report",
            "value": 168086,
            "range": "± 8115",
            "unit": "ns/iter"
          },
          {
            "name": "report_generation/generate_markdown_report",
            "value": 118663,
            "range": "± 5867",
            "unit": "ns/iter"
          },
          {
            "name": "report_generation/generate_csv_report",
            "value": 85740,
            "range": "± 4215",
            "unit": "ns/iter"
          },
          {
            "name": "dataset_operations/create_dataset_info",
            "value": 80,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "dataset_operations/validate_dataset_info",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "comparison_operations/create_comparison",
            "value": 1,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "comparison_operations/check_targets",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_gc_content",
            "value": 815214,
            "range": "± 12151",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_revcomp",
            "value": 297751,
            "range": "± 4442",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_subsequence",
            "value": 167,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "memory_pool_allocation",
            "value": 73,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "cache_put_get",
            "value": 109742,
            "range": "± 588",
            "unit": "ns/iter"
          },
          {
            "name": "interval_query",
            "value": 40648369,
            "range": "± 1246386",
            "unit": "ns/iter"
          },
          {
            "name": "variant_type",
            "value": 15890,
            "range": "± 179",
            "unit": "ns/iter"
          },
          {
            "name": "quality_error_prob",
            "value": 1216,
            "range": "± 34",
            "unit": "ns/iter"
          },
          {
            "name": "log_addition",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "stream_processing",
            "value": 595,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "memory_mapped_access",
            "value": 9848,
            "range": "± 157",
            "unit": "ns/iter"
          },
          {
            "name": "hamming_distance",
            "value": 9,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/100",
            "value": 29811,
            "range": "± 145",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/100",
            "value": 51959,
            "range": "± 2685",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/100",
            "value": 26242,
            "range": "± 216",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/1000",
            "value": 246733,
            "range": "± 668",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/1000",
            "value": 376938,
            "range": "± 8073",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/1000",
            "value": 207606,
            "range": "± 1401",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/10000",
            "value": 2459549,
            "range": "± 10780",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/10000",
            "value": 3566809,
            "range": "± 72260",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/10000",
            "value": 2084115,
            "range": "± 22441",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/100",
            "value": 33404,
            "range": "± 114",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/100",
            "value": 42508,
            "range": "± 172",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/100",
            "value": 27164,
            "range": "± 247",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/1000",
            "value": 279207,
            "range": "± 2445",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/1000",
            "value": 286571,
            "range": "± 1076",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/1000",
            "value": 216980,
            "range": "± 2234",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/10000",
            "value": 2730452,
            "range": "± 26501",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/10000",
            "value": 2721268,
            "range": "± 24608",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/10000",
            "value": 2152796,
            "range": "± 14568",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/gc_content",
            "value": 63,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/reverse_complement",
            "value": 60,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/subsequence",
            "value": 23,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/count_pattern",
            "value": 225,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/average_quality",
            "value": 62,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/min_quality",
            "value": 45,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/max_quality",
            "value": 49,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/median_quality",
            "value": 67,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/trim_quality",
            "value": 53,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/reverse_complement",
            "value": 95,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/is_valid",
            "value": 92,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/build_index",
            "value": 344922,
            "range": "± 4852",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/random_access",
            "value": 4940,
            "range": "± 33",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/save_index",
            "value": 216252,
            "range": "± 19132",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/load_index",
            "value": 302736,
            "range": "± 913",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/filter_by_quality",
            "value": 332454,
            "range": "± 3019",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/filter_by_length",
            "value": 186749,
            "range": "± 949",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/sample_reads",
            "value": 189155,
            "range": "± 798",
            "unit": "ns/iter"
          },
          {
            "name": "io_writing/write_fasta",
            "value": 275137,
            "range": "± 10417",
            "unit": "ns/iter"
          },
          {
            "name": "io_writing/write_fastq",
            "value": 8747506,
            "range": "± 700470",
            "unit": "ns/iter"
          },
          {
            "name": "memory_usage/parse_memory_usage",
            "value": 2482033,
            "range": "± 10424",
            "unit": "ns/iter"
          },
          {
            "name": "macro_sum_bench",
            "value": 1,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/100",
            "value": 90695,
            "range": "± 1094",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/1000",
            "value": 890617,
            "range": "± 7636",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/10000",
            "value": 9470252,
            "range": "± 199285",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_iterator/iterator",
            "value": 610618,
            "range": "± 5420",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_writing/writing",
            "value": 3763497,
            "range": "± 22043",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/0",
            "value": 107,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/1",
            "value": 97,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/2",
            "value": 94,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/3",
            "value": 132,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/4",
            "value": 116,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/5",
            "value": 60,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/6",
            "value": 114,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/7",
            "value": 494,
            "range": "± 9",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/0",
            "value": 91,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/1",
            "value": 91,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/2",
            "value": 47,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/3",
            "value": 91,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/4",
            "value": 91,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/5",
            "value": 92,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/6",
            "value": 91,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/7",
            "value": 91,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/variant_type_detection",
            "value": 2203,
            "range": "± 12",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/allele_frequency_access",
            "value": 4249,
            "range": "± 6",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/sample_data_access",
            "value": 1517,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "band_pass_normalized_kernel_gatk_defaults",
            "value": 1007,
            "range": "± 9",
            "unit": "ns/iter"
          },
          {
            "name": "band_pass_add_then_pop_ready_256_loci",
            "value": 324533,
            "range": "± 851",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_medium_depth_k10",
            "value": 259352,
            "range": "± 17062",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_high_depth_k10",
            "value": 1672010,
            "range": "± 1575",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/4",
            "value": 13590684,
            "range": "± 1373967",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/16",
            "value": 47510414,
            "range": "± 763449",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/64",
            "value": 164750035,
            "range": "± 21507286",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/8",
            "value": 541263,
            "range": "± 835",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/8",
            "value": 220399,
            "range": "± 1329",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/8",
            "value": 1090724,
            "range": "± 4060",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/8",
            "value": 9007818,
            "range": "± 16282",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/32",
            "value": 2163343,
            "range": "± 4603",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/32",
            "value": 829929,
            "range": "± 1773",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/32",
            "value": 4332444,
            "range": "± 6344",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/32",
            "value": 36025357,
            "range": "± 90805",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/64",
            "value": 4326693,
            "range": "± 6052",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/64",
            "value": 1651901,
            "range": "± 7017",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/64",
            "value": 8656957,
            "range": "± 9107",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/64",
            "value": 72066513,
            "range": "± 163810",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/128",
            "value": 8655357,
            "range": "± 9908",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/128",
            "value": 3298459,
            "range": "± 6124",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/128",
            "value": 17305002,
            "range": "± 19480",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/128",
            "value": 144102479,
            "range": "± 1125338",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/8",
            "value": 2123392,
            "range": "± 2237",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/8",
            "value": 4858989,
            "range": "± 12980",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/8",
            "value": 4328975,
            "range": "± 7445",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/8",
            "value": 36420164,
            "range": "± 114702",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/32",
            "value": 8470472,
            "range": "± 6801",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/32",
            "value": 19475577,
            "range": "± 128365",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/32",
            "value": 17244517,
            "range": "± 18556",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/32",
            "value": 145640418,
            "range": "± 457818",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/64",
            "value": 16941017,
            "range": "± 17692",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/64",
            "value": 39007233,
            "range": "± 67799",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/64",
            "value": 34478001,
            "range": "± 58706",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/64",
            "value": 291403662,
            "range": "± 353713",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/128",
            "value": 33879397,
            "range": "± 53955",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/128",
            "value": 77461255,
            "range": "± 202632",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/128",
            "value": 69158415,
            "range": "± 124435",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/128",
            "value": 582881492,
            "range": "± 1228544",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/8",
            "value": 4718255,
            "range": "± 8068",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/8",
            "value": 9410909,
            "range": "± 29992",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/8",
            "value": 9671459,
            "range": "± 8643",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/8",
            "value": 82186313,
            "range": "± 187818",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/32",
            "value": 18865898,
            "range": "± 20966",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/32",
            "value": 38261498,
            "range": "± 363058",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/32",
            "value": 38574368,
            "range": "± 42405",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/32",
            "value": 328873130,
            "range": "± 570236",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/64",
            "value": 37958323,
            "range": "± 36266",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/64",
            "value": 76249000,
            "range": "± 403194",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/64",
            "value": 77395366,
            "range": "± 228751",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/64",
            "value": 657775339,
            "range": "± 487348",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/128",
            "value": 75958535,
            "range": "± 1177892",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/128",
            "value": 152849356,
            "range": "± 857585",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/128",
            "value": 154746295,
            "range": "± 123114",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/128",
            "value": 1315704005,
            "range": "± 694288",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/value_copy_stages/3072x12",
            "value": 5213,
            "range": "± 31",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/arc_share_stages/3072x12",
            "value": 3921,
            "range": "± 13",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/assembly_result_set_arc_fanout/3072",
            "value": 74,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/assembly_result_set_value_fanout/3072",
            "value": 1155,
            "range": "± 7",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/64x48",
            "value": 17656,
            "range": "± 438",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/64x48",
            "value": 17560,
            "range": "± 393",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/128x96",
            "value": 70208,
            "range": "± 3488",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/128x96",
            "value": 69525,
            "range": "± 213",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/256x192",
            "value": 267481,
            "range": "± 821",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/256x192",
            "value": 265219,
            "range": "± 1145",
            "unit": "ns/iter"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "name": "SynapticFour",
            "username": "SynapticFour",
            "email": "contact@synapticfour.com"
          },
          "committer": {
            "name": "GitHub",
            "username": "web-flow",
            "email": "noreply@github.com"
          },
          "id": "4b008904b4cdbb678e62cca1f84560fa96ac7060",
          "message": "perf(hc): apply 8k graph node cap on primary assemble paths (#97)\n\nRT-supplement already skipped k-best above 8k nodes; primary SeqGraph\nand RT assemble did not, so bushy dense 1 Mb GIAB shards could still\nexpand into multi-GiB Peak. Share MAX_ASSEMBLY_GRAPH_NODES across\nprimary SeqGraph/RT and supplement extract.\n\nCo-authored-by: Cursor <cursoragent@cursor.com>",
          "timestamp": "2026-08-10T15:21:08Z",
          "url": "https://github.com/SynapticFour/gatk-rs/commit/4b008904b4cdbb678e62cca1f84560fa96ac7060"
        },
        "date": 1786420139819,
        "tool": "cargo",
        "benches": [
          {
            "name": "sam_parsing/parsing/100",
            "value": 61579,
            "range": "± 1789",
            "unit": "ns/iter"
          },
          {
            "name": "sam_parsing/parsing/1000",
            "value": 712403,
            "range": "± 11478",
            "unit": "ns/iter"
          },
          {
            "name": "sam_parsing/parsing/10000",
            "value": 6760583,
            "range": "± 189799",
            "unit": "ns/iter"
          },
          {
            "name": "sam_iterator/iterator",
            "value": 482779,
            "range": "± 28131",
            "unit": "ns/iter"
          },
          {
            "name": "sam_writing/writing",
            "value": 2639351,
            "range": "± 4351911",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/0",
            "value": 63,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/1",
            "value": 94,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/2",
            "value": 71,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/3",
            "value": 73,
            "range": "± 3",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/4",
            "value": 97,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "optional_fields/parsing",
            "value": 4902,
            "range": "± 123",
            "unit": "ns/iter"
          },
          {
            "name": "framework_overhead/create_result",
            "value": 13,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "framework_overhead/create_comparison",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "metrics_collection/metrics_collector",
            "value": 1059935,
            "range": "± 8804",
            "unit": "ns/iter"
          },
          {
            "name": "suite_operations/create_suite",
            "value": 13249835,
            "range": "± 675016",
            "unit": "ns/iter"
          },
          {
            "name": "suite_operations/suite_with_results",
            "value": 12374156,
            "range": "± 1063476",
            "unit": "ns/iter"
          },
          {
            "name": "suite_operations/suite_serialization",
            "value": 41698,
            "range": "± 1479",
            "unit": "ns/iter"
          },
          {
            "name": "performance_analysis/analyze_small_dataset",
            "value": 48,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "performance_analysis/analyze_large_dataset",
            "value": 5127,
            "range": "± 200",
            "unit": "ns/iter"
          },
          {
            "name": "regression_detection/detect_no_regression",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "regression_detection/detect_regression",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "report_generation/generate_json_report",
            "value": 249499,
            "range": "± 917650",
            "unit": "ns/iter"
          },
          {
            "name": "report_generation/generate_markdown_report",
            "value": 237237,
            "range": "± 2458487",
            "unit": "ns/iter"
          },
          {
            "name": "report_generation/generate_csv_report",
            "value": 190223,
            "range": "± 1376965",
            "unit": "ns/iter"
          },
          {
            "name": "dataset_operations/create_dataset_info",
            "value": 53,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "dataset_operations/validate_dataset_info",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "comparison_operations/create_comparison",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "comparison_operations/check_targets",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_gc_content",
            "value": 489167,
            "range": "± 18391",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_revcomp",
            "value": 209529,
            "range": "± 5464",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_subsequence",
            "value": 104,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "memory_pool_allocation",
            "value": 36,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "cache_put_get",
            "value": 61878,
            "range": "± 1823",
            "unit": "ns/iter"
          },
          {
            "name": "interval_query",
            "value": 24290328,
            "range": "± 587889",
            "unit": "ns/iter"
          },
          {
            "name": "variant_type",
            "value": 7299,
            "range": "± 240",
            "unit": "ns/iter"
          },
          {
            "name": "quality_error_prob",
            "value": 709,
            "range": "± 14",
            "unit": "ns/iter"
          },
          {
            "name": "log_addition",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "stream_processing",
            "value": 359,
            "range": "± 11",
            "unit": "ns/iter"
          },
          {
            "name": "memory_mapped_access",
            "value": 6608,
            "range": "± 167",
            "unit": "ns/iter"
          },
          {
            "name": "hamming_distance",
            "value": 5,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/100",
            "value": 18051,
            "range": "± 647",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/100",
            "value": 33506,
            "range": "± 1484",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/100",
            "value": 14581,
            "range": "± 289",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/1000",
            "value": 146667,
            "range": "± 5621",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/1000",
            "value": 249426,
            "range": "± 9100",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/1000",
            "value": 118459,
            "range": "± 3965",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/10000",
            "value": 1497018,
            "range": "± 23070",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/10000",
            "value": 2319119,
            "range": "± 112724",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/10000",
            "value": 1115108,
            "range": "± 28681",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/100",
            "value": 19676,
            "range": "± 488",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/100",
            "value": 26371,
            "range": "± 907",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/100",
            "value": 14980,
            "range": "± 647",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/1000",
            "value": 159415,
            "range": "± 7626",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/1000",
            "value": 166291,
            "range": "± 3007",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/1000",
            "value": 117300,
            "range": "± 6937",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/10000",
            "value": 1554697,
            "range": "± 59157",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/10000",
            "value": 1535160,
            "range": "± 46370",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/10000",
            "value": 1122068,
            "range": "± 60949",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/gc_content",
            "value": 33,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/reverse_complement",
            "value": 33,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/subsequence",
            "value": 13,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/count_pattern",
            "value": 117,
            "range": "± 4",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/average_quality",
            "value": 26,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/min_quality",
            "value": 30,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/max_quality",
            "value": 42,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/median_quality",
            "value": 22,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/trim_quality",
            "value": 26,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/reverse_complement",
            "value": 62,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/is_valid",
            "value": 54,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/build_index",
            "value": 213013,
            "range": "± 5989",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/random_access",
            "value": 3488,
            "range": "± 65",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/save_index",
            "value": 729071,
            "range": "± 1031660",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/load_index",
            "value": 188589,
            "range": "± 4255",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/filter_by_quality",
            "value": 190527,
            "range": "± 2218",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/filter_by_length",
            "value": 131695,
            "range": "± 15268",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/sample_reads",
            "value": 100836,
            "range": "± 2340",
            "unit": "ns/iter"
          },
          {
            "name": "io_writing/write_fasta",
            "value": 605464,
            "range": "± 1188369",
            "unit": "ns/iter"
          },
          {
            "name": "io_writing/write_fastq",
            "value": 6676460,
            "range": "± 4314144",
            "unit": "ns/iter"
          },
          {
            "name": "memory_usage/parse_memory_usage",
            "value": 1552627,
            "range": "± 31367",
            "unit": "ns/iter"
          },
          {
            "name": "macro_sum_bench",
            "value": 1,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/100",
            "value": 51554,
            "range": "± 1681",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/1000",
            "value": 587597,
            "range": "± 13576",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/10000",
            "value": 6286479,
            "range": "± 118188",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_iterator/iterator",
            "value": 342826,
            "range": "± 8436",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_writing/writing",
            "value": 3102471,
            "range": "± 2187647",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/0",
            "value": 62,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/1",
            "value": 56,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/2",
            "value": 55,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/3",
            "value": 74,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/4",
            "value": 63,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/5",
            "value": 32,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/6",
            "value": 60,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/7",
            "value": 311,
            "range": "± 9",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/0",
            "value": 52,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/1",
            "value": 53,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/2",
            "value": 24,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/3",
            "value": 53,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/4",
            "value": 52,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/5",
            "value": 53,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/6",
            "value": 53,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/7",
            "value": 53,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/variant_type_detection",
            "value": 1581,
            "range": "± 48",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/allele_frequency_access",
            "value": 2797,
            "range": "± 54",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/sample_data_access",
            "value": 951,
            "range": "± 39",
            "unit": "ns/iter"
          },
          {
            "name": "band_pass_normalized_kernel_gatk_defaults",
            "value": 592,
            "range": "± 33",
            "unit": "ns/iter"
          },
          {
            "name": "band_pass_add_then_pop_ready_256_loci",
            "value": 190916,
            "range": "± 2019",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_medium_depth_k10",
            "value": 137593,
            "range": "± 7302",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_high_depth_k10",
            "value": 1006090,
            "range": "± 7423",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/4",
            "value": 8702568,
            "range": "± 384045",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/16",
            "value": 22528808,
            "range": "± 472855",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/64",
            "value": 116374041,
            "range": "± 4658309",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/8",
            "value": 335892,
            "range": "± 16542",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/8",
            "value": 121035,
            "range": "± 3301",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/8",
            "value": 703010,
            "range": "± 16004",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/8",
            "value": 5151179,
            "range": "± 130187",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/32",
            "value": 1318614,
            "range": "± 24246",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/32",
            "value": 481392,
            "range": "± 6883",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/32",
            "value": 2859771,
            "range": "± 140240",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/32",
            "value": 20344099,
            "range": "± 758119",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/64",
            "value": 2651542,
            "range": "± 58442",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/64",
            "value": 1015215,
            "range": "± 57365",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/64",
            "value": 5734886,
            "range": "± 204284",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/64",
            "value": 42481509,
            "range": "± 1756001",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/128",
            "value": 5386399,
            "range": "± 143157",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/128",
            "value": 2023445,
            "range": "± 24428",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/128",
            "value": 11663926,
            "range": "± 179037",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/128",
            "value": 85300183,
            "range": "± 1167033",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/8",
            "value": 1335925,
            "range": "± 20959",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/8",
            "value": 3601350,
            "range": "± 85444",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/8",
            "value": 2957909,
            "range": "± 165315",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/8",
            "value": 21851767,
            "range": "± 331128",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/32",
            "value": 5282660,
            "range": "± 115889",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/32",
            "value": 14159125,
            "range": "± 1206766",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/32",
            "value": 11309271,
            "range": "± 387366",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/32",
            "value": 85592144,
            "range": "± 2671541",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/64",
            "value": 10704085,
            "range": "± 134197",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/64",
            "value": 28737978,
            "range": "± 327611",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/64",
            "value": 23355880,
            "range": "± 548993",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/64",
            "value": 172893124,
            "range": "± 5417296",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/128",
            "value": 21634804,
            "range": "± 736358",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/128",
            "value": 58148328,
            "range": "± 657389",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/128",
            "value": 48060181,
            "range": "± 963115",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/128",
            "value": 349192043,
            "range": "± 6620675",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/8",
            "value": 3037943,
            "range": "± 33423",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/8",
            "value": 5844020,
            "range": "± 132222",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/8",
            "value": 6665391,
            "range": "± 56163",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/8",
            "value": 49108341,
            "range": "± 1573610",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/32",
            "value": 12216378,
            "range": "± 130955",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/32",
            "value": 23008259,
            "range": "± 446164",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/32",
            "value": 26778223,
            "range": "± 265400",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/32",
            "value": 195197912,
            "range": "± 1901295",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/64",
            "value": 24183446,
            "range": "± 321469",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/64",
            "value": 46342306,
            "range": "± 873729",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/64",
            "value": 53810994,
            "range": "± 483821",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/64",
            "value": 392324969,
            "range": "± 6947980",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/128",
            "value": 48297910,
            "range": "± 1616528",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/128",
            "value": 92852257,
            "range": "± 6070273",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/128",
            "value": 107216722,
            "range": "± 1761662",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/128",
            "value": 790280425,
            "range": "± 8960652",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/value_copy_stages/3072x12",
            "value": 3323,
            "range": "± 43",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/arc_share_stages/3072x12",
            "value": 2617,
            "range": "± 50",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/assembly_result_set_arc_fanout/3072",
            "value": 127,
            "range": "± 3",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/assembly_result_set_value_fanout/3072",
            "value": 782,
            "range": "± 9",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/64x48",
            "value": 10716,
            "range": "± 237",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/64x48",
            "value": 10610,
            "range": "± 284",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/128x96",
            "value": 39287,
            "range": "± 784",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/128x96",
            "value": 36599,
            "range": "± 437",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/256x192",
            "value": 144088,
            "range": "± 1563",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/256x192",
            "value": 141573,
            "range": "± 3652",
            "unit": "ns/iter"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "name": "SynapticFour",
            "username": "SynapticFour",
            "email": "contact@synapticfour.com"
          },
          "committer": {
            "name": "GitHub",
            "username": "web-flow",
            "email": "noreply@github.com"
          },
          "id": "830f241ecac80a8f4004cc320ca3b4dee88fa08f",
          "message": "perf(hc): cut RT graph-build Peak via ownership, not abort (#100)\n\nShare kmer bytes (Arc), build the threading graph once, and stop\nunbounded dangling ref walks so dense chr20 windows stay tens of MiB.\n\nCo-authored-by: Cursor <cursoragent@cursor.com>",
          "timestamp": "2026-08-11T18:36:02Z",
          "url": "https://github.com/SynapticFour/gatk-rs/commit/830f241ecac80a8f4004cc320ca3b4dee88fa08f"
        },
        "date": 1786508175399,
        "tool": "cargo",
        "benches": [
          {
            "name": "sam_parsing/parsing/100",
            "value": 113944,
            "range": "± 969",
            "unit": "ns/iter"
          },
          {
            "name": "sam_parsing/parsing/1000",
            "value": 1167227,
            "range": "± 42625",
            "unit": "ns/iter"
          },
          {
            "name": "sam_parsing/parsing/10000",
            "value": 11510575,
            "range": "± 42797",
            "unit": "ns/iter"
          },
          {
            "name": "sam_iterator/iterator",
            "value": 858683,
            "range": "± 8546",
            "unit": "ns/iter"
          },
          {
            "name": "sam_writing/writing",
            "value": 3232628,
            "range": "± 62113",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/0",
            "value": 134,
            "range": "± 8",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/1",
            "value": 179,
            "range": "± 12",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/2",
            "value": 141,
            "range": "± 7",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/3",
            "value": 145,
            "range": "± 4",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/4",
            "value": 182,
            "range": "± 6",
            "unit": "ns/iter"
          },
          {
            "name": "optional_fields/parsing",
            "value": 6602,
            "range": "± 78",
            "unit": "ns/iter"
          },
          {
            "name": "framework_overhead/create_result",
            "value": 25,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "framework_overhead/create_comparison",
            "value": 1,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "metrics_collection/metrics_collector",
            "value": 1066959,
            "range": "± 4207",
            "unit": "ns/iter"
          },
          {
            "name": "suite_operations/create_suite",
            "value": 16021864,
            "range": "± 414371",
            "unit": "ns/iter"
          },
          {
            "name": "suite_operations/suite_with_results",
            "value": 16237378,
            "range": "± 421884",
            "unit": "ns/iter"
          },
          {
            "name": "suite_operations/suite_serialization",
            "value": 75056,
            "range": "± 991",
            "unit": "ns/iter"
          },
          {
            "name": "performance_analysis/analyze_small_dataset",
            "value": 82,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "performance_analysis/analyze_large_dataset",
            "value": 7681,
            "range": "± 79",
            "unit": "ns/iter"
          },
          {
            "name": "regression_detection/detect_no_regression",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "regression_detection/detect_regression",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "report_generation/generate_json_report",
            "value": 198847,
            "range": "± 41136",
            "unit": "ns/iter"
          },
          {
            "name": "report_generation/generate_markdown_report",
            "value": 183722,
            "range": "± 11833",
            "unit": "ns/iter"
          },
          {
            "name": "report_generation/generate_csv_report",
            "value": 182481,
            "range": "± 38981",
            "unit": "ns/iter"
          },
          {
            "name": "dataset_operations/create_dataset_info",
            "value": 74,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "dataset_operations/validate_dataset_info",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "comparison_operations/create_comparison",
            "value": 1,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "comparison_operations/check_targets",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_gc_content",
            "value": 719342,
            "range": "± 8625",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_revcomp",
            "value": 258991,
            "range": "± 1473",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_subsequence",
            "value": 199,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "memory_pool_allocation",
            "value": 70,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "cache_put_get",
            "value": 109197,
            "range": "± 710",
            "unit": "ns/iter"
          },
          {
            "name": "interval_query",
            "value": 43287852,
            "range": "± 2351769",
            "unit": "ns/iter"
          },
          {
            "name": "variant_type",
            "value": 16226,
            "range": "± 26",
            "unit": "ns/iter"
          },
          {
            "name": "quality_error_prob",
            "value": 1274,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "log_addition",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "stream_processing",
            "value": 533,
            "range": "± 3",
            "unit": "ns/iter"
          },
          {
            "name": "memory_mapped_access",
            "value": 8887,
            "range": "± 352",
            "unit": "ns/iter"
          },
          {
            "name": "hamming_distance",
            "value": 8,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/100",
            "value": 31072,
            "range": "± 263",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/100",
            "value": 55082,
            "range": "± 892",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/100",
            "value": 26147,
            "range": "± 265",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/1000",
            "value": 271255,
            "range": "± 2000",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/1000",
            "value": 422315,
            "range": "± 2210",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/1000",
            "value": 220708,
            "range": "± 625",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/10000",
            "value": 2692811,
            "range": "± 15449",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/10000",
            "value": 4057654,
            "range": "± 23575",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/10000",
            "value": 2215730,
            "range": "± 20431",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/100",
            "value": 32733,
            "range": "± 193",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/100",
            "value": 43867,
            "range": "± 287",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/100",
            "value": 24179,
            "range": "± 107",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/1000",
            "value": 271586,
            "range": "± 1580",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/1000",
            "value": 299246,
            "range": "± 1049",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/1000",
            "value": 201931,
            "range": "± 4886",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/10000",
            "value": 2555362,
            "range": "± 31023",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/10000",
            "value": 2746052,
            "range": "± 22057",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/10000",
            "value": 1979858,
            "range": "± 27878",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/gc_content",
            "value": 57,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/reverse_complement",
            "value": 56,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/subsequence",
            "value": 20,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/count_pattern",
            "value": 210,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/average_quality",
            "value": 53,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/min_quality",
            "value": 41,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/max_quality",
            "value": 44,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/median_quality",
            "value": 42,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/trim_quality",
            "value": 50,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/reverse_complement",
            "value": 88,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/is_valid",
            "value": 82,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/build_index",
            "value": 359925,
            "range": "± 2440",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/random_access",
            "value": 4597,
            "range": "± 91",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/save_index",
            "value": 305047,
            "range": "± 11110",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/load_index",
            "value": 288459,
            "range": "± 28985",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/filter_by_quality",
            "value": 328436,
            "range": "± 2081",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/filter_by_length",
            "value": 184437,
            "range": "± 446",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/sample_reads",
            "value": 183816,
            "range": "± 720",
            "unit": "ns/iter"
          },
          {
            "name": "io_writing/write_fasta",
            "value": 266282,
            "range": "± 32790",
            "unit": "ns/iter"
          },
          {
            "name": "io_writing/write_fastq",
            "value": 8475412,
            "range": "± 176506",
            "unit": "ns/iter"
          },
          {
            "name": "memory_usage/parse_memory_usage",
            "value": 2687124,
            "range": "± 117746",
            "unit": "ns/iter"
          },
          {
            "name": "macro_sum_bench",
            "value": 1,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/100",
            "value": 91634,
            "range": "± 1213",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/1000",
            "value": 900851,
            "range": "± 10507",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/10000",
            "value": 9811934,
            "range": "± 259178",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_iterator/iterator",
            "value": 623672,
            "range": "± 3734",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_writing/writing",
            "value": 3780096,
            "range": "± 73241",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/0",
            "value": 98,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/1",
            "value": 93,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/2",
            "value": 91,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/3",
            "value": 122,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/4",
            "value": 105,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/5",
            "value": 70,
            "range": "± 5",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/6",
            "value": 108,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/7",
            "value": 473,
            "range": "± 7",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/0",
            "value": 86,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/1",
            "value": 85,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/2",
            "value": 41,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/3",
            "value": 85,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/4",
            "value": 85,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/5",
            "value": 85,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/6",
            "value": 85,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/7",
            "value": 85,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/variant_type_detection",
            "value": 2507,
            "range": "± 17",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/allele_frequency_access",
            "value": 4770,
            "range": "± 37",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/sample_data_access",
            "value": 1384,
            "range": "± 8",
            "unit": "ns/iter"
          },
          {
            "name": "band_pass_normalized_kernel_gatk_defaults",
            "value": 994,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "band_pass_add_then_pop_ready_256_loci",
            "value": 296084,
            "range": "± 2130",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_medium_depth_k10",
            "value": 234545,
            "range": "± 308",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_high_depth_k10",
            "value": 1664987,
            "range": "± 9388",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/4",
            "value": 13725555,
            "range": "± 85118",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/16",
            "value": 48873366,
            "range": "± 836257",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/64",
            "value": 177122417,
            "range": "± 22773994",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/8",
            "value": 498861,
            "range": "± 1366",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/8",
            "value": 213212,
            "range": "± 1494",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/8",
            "value": 1091593,
            "range": "± 15458",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/8",
            "value": 9459422,
            "range": "± 31895",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/32",
            "value": 1990588,
            "range": "± 5011",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/32",
            "value": 816706,
            "range": "± 7180",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/32",
            "value": 4305808,
            "range": "± 17369",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/32",
            "value": 37869081,
            "range": "± 77701",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/64",
            "value": 3981638,
            "range": "± 15958",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/64",
            "value": 1627370,
            "range": "± 7795",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/64",
            "value": 8608413,
            "range": "± 25499",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/64",
            "value": 75762068,
            "range": "± 205619",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/128",
            "value": 7972064,
            "range": "± 37114",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/128",
            "value": 3250505,
            "range": "± 14048",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/128",
            "value": 17182281,
            "range": "± 239103",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/128",
            "value": 151412801,
            "range": "± 391999",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/8",
            "value": 1958360,
            "range": "± 23740",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/8",
            "value": 4909516,
            "range": "± 104825",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/8",
            "value": 4299872,
            "range": "± 13055",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/8",
            "value": 38267091,
            "range": "± 137297",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/32",
            "value": 7826574,
            "range": "± 25801",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/32",
            "value": 18773153,
            "range": "± 193048",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/32",
            "value": 17151133,
            "range": "± 55915",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/32",
            "value": 153281370,
            "range": "± 401166",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/64",
            "value": 15659145,
            "range": "± 48637",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/64",
            "value": 39405301,
            "range": "± 893747",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/64",
            "value": 34263612,
            "range": "± 165712",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/64",
            "value": 306866056,
            "range": "± 1163858",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/128",
            "value": 31132363,
            "range": "± 96040",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/128",
            "value": 75803732,
            "range": "± 1302412",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/128",
            "value": 68400892,
            "range": "± 239557",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/128",
            "value": 613863172,
            "range": "± 1433692",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/8",
            "value": 4324492,
            "range": "± 20117",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/8",
            "value": 8082597,
            "range": "± 245026",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/8",
            "value": 9572227,
            "range": "± 55961",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/8",
            "value": 86665908,
            "range": "± 284618",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/32",
            "value": 17289841,
            "range": "± 61099",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/32",
            "value": 31323492,
            "range": "± 990806",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/32",
            "value": 38179877,
            "range": "± 183069",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/32",
            "value": 346620646,
            "range": "± 1407786",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/64",
            "value": 34223906,
            "range": "± 113934",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/64",
            "value": 61689901,
            "range": "± 1470425",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/64",
            "value": 75958405,
            "range": "± 261913",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/64",
            "value": 692746135,
            "range": "± 1490080",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/128",
            "value": 68522296,
            "range": "± 153814",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/128",
            "value": 124421322,
            "range": "± 3613465",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/128",
            "value": 151847704,
            "range": "± 435592",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/128",
            "value": 1384324513,
            "range": "± 3163132",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/value_copy_stages/3072x12",
            "value": 4904,
            "range": "± 27",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/arc_share_stages/3072x12",
            "value": 3476,
            "range": "± 15",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/assembly_result_set_arc_fanout/3072",
            "value": 62,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/assembly_result_set_value_fanout/3072",
            "value": 1217,
            "range": "± 10",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/64x48",
            "value": 18667,
            "range": "± 77",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/64x48",
            "value": 18429,
            "range": "± 1146",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/128x96",
            "value": 70245,
            "range": "± 791",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/128x96",
            "value": 70314,
            "range": "± 371",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/256x192",
            "value": 271540,
            "range": "± 870",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/256x192",
            "value": 270807,
            "range": "± 604",
            "unit": "ns/iter"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "name": "SynapticFour",
            "username": "SynapticFour",
            "email": "contact@synapticfour.com"
          },
          "committer": {
            "name": "GitHub",
            "username": "web-flow",
            "email": "noreply@github.com"
          },
          "id": "007ae487a8bbd2ac3e7512b17d48ea1ab63dfdc3",
          "message": "fix(hc): attach match CIGAR on just-reference haplotypes (#105)\n\nRT-first seeds from just_reference_result; a cigarless seed ref won\ndedup and emptied the E2E ref cigar column (p5_indel_chrindel after #104).\n\nCo-authored-by: Cursor <cursoragent@cursor.com>",
          "timestamp": "2026-08-12T16:11:07Z",
          "url": "https://github.com/SynapticFour/gatk-rs/commit/007ae487a8bbd2ac3e7512b17d48ea1ab63dfdc3"
        },
        "date": 1786593931873,
        "tool": "cargo",
        "benches": [
          {
            "name": "sam_parsing/parsing/100",
            "value": 68774,
            "range": "± 338",
            "unit": "ns/iter"
          },
          {
            "name": "sam_parsing/parsing/1000",
            "value": 751400,
            "range": "± 4638",
            "unit": "ns/iter"
          },
          {
            "name": "sam_parsing/parsing/10000",
            "value": 7579358,
            "range": "± 42274",
            "unit": "ns/iter"
          },
          {
            "name": "sam_iterator/iterator",
            "value": 527474,
            "range": "± 7156",
            "unit": "ns/iter"
          },
          {
            "name": "sam_writing/writing",
            "value": 1456882,
            "range": "± 4347932",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/0",
            "value": 82,
            "range": "± 9",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/1",
            "value": 116,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/2",
            "value": 92,
            "range": "± 5",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/3",
            "value": 94,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/4",
            "value": 119,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "optional_fields/parsing",
            "value": 2694,
            "range": "± 9",
            "unit": "ns/iter"
          },
          {
            "name": "framework_overhead/create_result",
            "value": 16,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "framework_overhead/create_comparison",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "metrics_collection/metrics_collector",
            "value": 1060148,
            "range": "± 5629",
            "unit": "ns/iter"
          },
          {
            "name": "suite_operations/create_suite",
            "value": 8690332,
            "range": "± 157976",
            "unit": "ns/iter"
          },
          {
            "name": "suite_operations/suite_with_results",
            "value": 8915063,
            "range": "± 179439",
            "unit": "ns/iter"
          },
          {
            "name": "suite_operations/suite_serialization",
            "value": 42606,
            "range": "± 262",
            "unit": "ns/iter"
          },
          {
            "name": "performance_analysis/analyze_small_dataset",
            "value": 49,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "performance_analysis/analyze_large_dataset",
            "value": 6472,
            "range": "± 6",
            "unit": "ns/iter"
          },
          {
            "name": "regression_detection/detect_no_regression",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "regression_detection/detect_regression",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "report_generation/generate_json_report",
            "value": 647487,
            "range": "± 1649607",
            "unit": "ns/iter"
          },
          {
            "name": "report_generation/generate_markdown_report",
            "value": 228299,
            "range": "± 1475695",
            "unit": "ns/iter"
          },
          {
            "name": "report_generation/generate_csv_report",
            "value": 257400,
            "range": "± 1454538",
            "unit": "ns/iter"
          },
          {
            "name": "dataset_operations/create_dataset_info",
            "value": 52,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "dataset_operations/validate_dataset_info",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "comparison_operations/create_comparison",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "comparison_operations/check_targets",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_gc_content",
            "value": 600242,
            "range": "± 20571",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_revcomp",
            "value": 245661,
            "range": "± 6506",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_subsequence",
            "value": 86,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "memory_pool_allocation",
            "value": 36,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "cache_put_get",
            "value": 60698,
            "range": "± 2321",
            "unit": "ns/iter"
          },
          {
            "name": "interval_query",
            "value": 25533884,
            "range": "± 1134040",
            "unit": "ns/iter"
          },
          {
            "name": "variant_type",
            "value": 9806,
            "range": "± 74",
            "unit": "ns/iter"
          },
          {
            "name": "quality_error_prob",
            "value": 833,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "log_addition",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "stream_processing",
            "value": 305,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "memory_mapped_access",
            "value": 2476,
            "range": "± 51",
            "unit": "ns/iter"
          },
          {
            "name": "hamming_distance",
            "value": 6,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/100",
            "value": 18646,
            "range": "± 73",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/100",
            "value": 34890,
            "range": "± 2458",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/100",
            "value": 15768,
            "range": "± 584",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/1000",
            "value": 173081,
            "range": "± 639",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/1000",
            "value": 275437,
            "range": "± 712",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/1000",
            "value": 142741,
            "range": "± 1644",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/10000",
            "value": 1789692,
            "range": "± 5836",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/10000",
            "value": 2584778,
            "range": "± 61884",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/10000",
            "value": 1437443,
            "range": "± 65425",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/100",
            "value": 19815,
            "range": "± 138",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/100",
            "value": 28139,
            "range": "± 1579",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/100",
            "value": 13891,
            "range": "± 70",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/1000",
            "value": 174443,
            "range": "± 4978",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/1000",
            "value": 203313,
            "range": "± 5723",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/1000",
            "value": 127129,
            "range": "± 8411",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/10000",
            "value": 1821512,
            "range": "± 45922",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/10000",
            "value": 2084663,
            "range": "± 15003",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/10000",
            "value": 1258511,
            "range": "± 82270",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/gc_content",
            "value": 41,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/reverse_complement",
            "value": 35,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/subsequence",
            "value": 13,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/count_pattern",
            "value": 89,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/average_quality",
            "value": 21,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/min_quality",
            "value": 41,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/max_quality",
            "value": 51,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/median_quality",
            "value": 38,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/trim_quality",
            "value": 30,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/reverse_complement",
            "value": 67,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/is_valid",
            "value": 70,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/build_index",
            "value": 235513,
            "range": "± 454",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/random_access",
            "value": 1424,
            "range": "± 9",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/save_index",
            "value": 6013588,
            "range": "± 23397004",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/load_index",
            "value": 177961,
            "range": "± 4733",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/filter_by_quality",
            "value": 210958,
            "range": "± 1130",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/filter_by_length",
            "value": 120597,
            "range": "± 634",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/sample_reads",
            "value": 120444,
            "range": "± 552",
            "unit": "ns/iter"
          },
          {
            "name": "io_writing/write_fasta",
            "value": 734521,
            "range": "± 2267368",
            "unit": "ns/iter"
          },
          {
            "name": "io_writing/write_fastq",
            "value": 3059829,
            "range": "± 1108516",
            "unit": "ns/iter"
          },
          {
            "name": "memory_usage/parse_memory_usage",
            "value": 1782861,
            "range": "± 6494",
            "unit": "ns/iter"
          },
          {
            "name": "macro_sum_bench",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/100",
            "value": 58180,
            "range": "± 1240",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/1000",
            "value": 620150,
            "range": "± 18113",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/10000",
            "value": 6902250,
            "range": "± 57735",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_iterator/iterator",
            "value": 397229,
            "range": "± 2018",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_writing/writing",
            "value": 1891962,
            "range": "± 11988651",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/0",
            "value": 59,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/1",
            "value": 56,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/2",
            "value": 55,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/3",
            "value": 71,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/4",
            "value": 63,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/5",
            "value": 35,
            "range": "± 3",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/6",
            "value": 70,
            "range": "± 4",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/7",
            "value": 371,
            "range": "± 3",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/0",
            "value": 54,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/1",
            "value": 54,
            "range": "± 3",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/2",
            "value": 25,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/3",
            "value": 55,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/4",
            "value": 55,
            "range": "± 3",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/5",
            "value": 54,
            "range": "± 3",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/6",
            "value": 53,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/7",
            "value": 53,
            "range": "± 3",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/variant_type_detection",
            "value": 1690,
            "range": "± 84",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/allele_frequency_access",
            "value": 2248,
            "range": "± 66",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/sample_data_access",
            "value": 942,
            "range": "± 24",
            "unit": "ns/iter"
          },
          {
            "name": "band_pass_normalized_kernel_gatk_defaults",
            "value": 654,
            "range": "± 8",
            "unit": "ns/iter"
          },
          {
            "name": "band_pass_add_then_pop_ready_256_loci",
            "value": 222934,
            "range": "± 6263",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_medium_depth_k10",
            "value": 150971,
            "range": "± 856",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_high_depth_k10",
            "value": 1121858,
            "range": "± 53042",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/4",
            "value": 6197194,
            "range": "± 294411",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/16",
            "value": 32667511,
            "range": "± 1186302",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/64",
            "value": 127726209,
            "range": "± 4427744",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/8",
            "value": 317645,
            "range": "± 12053",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/8",
            "value": 139724,
            "range": "± 6198",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/8",
            "value": 625838,
            "range": "± 18792",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/8",
            "value": 6105060,
            "range": "± 43294",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/32",
            "value": 1240580,
            "range": "± 8236",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/32",
            "value": 558917,
            "range": "± 26966",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/32",
            "value": 2541426,
            "range": "± 115108",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/32",
            "value": 24916703,
            "range": "± 188297",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/64",
            "value": 2533085,
            "range": "± 107590",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/64",
            "value": 1101626,
            "range": "± 14207",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/64",
            "value": 4965557,
            "range": "± 191091",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/64",
            "value": 49966676,
            "range": "± 1808012",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/128",
            "value": 5072892,
            "range": "± 140728",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/128",
            "value": 2148143,
            "range": "± 90526",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/128",
            "value": 9915873,
            "range": "± 265607",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/128",
            "value": 100121546,
            "range": "± 5574597",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/8",
            "value": 1221102,
            "range": "± 52299",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/8",
            "value": 2531537,
            "range": "± 20425",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/8",
            "value": 2481927,
            "range": "± 12958",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/8",
            "value": 24713651,
            "range": "± 1254817",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/32",
            "value": 4771358,
            "range": "± 150585",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/32",
            "value": 10079532,
            "range": "± 303024",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/32",
            "value": 10076311,
            "range": "± 503879",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/32",
            "value": 100349003,
            "range": "± 4061260",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/64",
            "value": 9767832,
            "range": "± 58478",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/64",
            "value": 20692736,
            "range": "± 896399",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/64",
            "value": 20153526,
            "range": "± 72415",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/64",
            "value": 201950393,
            "range": "± 5982389",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/128",
            "value": 21713753,
            "range": "± 1062094",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/128",
            "value": 41124442,
            "range": "± 873215",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/128",
            "value": 40310411,
            "range": "± 1768540",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/128",
            "value": 403863365,
            "range": "± 9611260",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/8",
            "value": 2935123,
            "range": "± 9772",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/8",
            "value": 5028394,
            "range": "± 49000",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/8",
            "value": 6181599,
            "range": "± 12028",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/8",
            "value": 55706811,
            "range": "± 130751",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/32",
            "value": 11476928,
            "range": "± 135952",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/32",
            "value": 19799949,
            "range": "± 28337",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/32",
            "value": 26712190,
            "range": "± 809642",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/32",
            "value": 228999375,
            "range": "± 7263058",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/64",
            "value": 23391009,
            "range": "± 744508",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/64",
            "value": 41009299,
            "range": "± 1351167",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/64",
            "value": 49620806,
            "range": "± 880167",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/64",
            "value": 446216562,
            "range": "± 6123284",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/128",
            "value": 45807733,
            "range": "± 925897",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/128",
            "value": 80450622,
            "range": "± 6265173",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/128",
            "value": 100746895,
            "range": "± 4522927",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/128",
            "value": 908642789,
            "range": "± 10708423",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/value_copy_stages/3072x12",
            "value": 3402,
            "range": "± 34",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/arc_share_stages/3072x12",
            "value": 2631,
            "range": "± 16",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/assembly_result_set_arc_fanout/3072",
            "value": 132,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/assembly_result_set_value_fanout/3072",
            "value": 860,
            "range": "± 13",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/64x48",
            "value": 11731,
            "range": "± 125",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/64x48",
            "value": 11917,
            "range": "± 1162",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/128x96",
            "value": 49017,
            "range": "± 3583",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/128x96",
            "value": 52892,
            "range": "± 4176",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/256x192",
            "value": 193280,
            "range": "± 1954",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/256x192",
            "value": 193149,
            "range": "± 1866",
            "unit": "ns/iter"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "name": "SynapticFour",
            "username": "SynapticFour",
            "email": "contact@synapticfour.com"
          },
          "committer": {
            "name": "GitHub",
            "username": "web-flow",
            "email": "noreply@github.com"
          },
          "id": "03445555bcea4bd70a4c75f31226b47ed90833ec",
          "message": "perf(hc): beat Java wall phase4 — indel CIGAR + loser path cuts (#109)\n\n* perf(hc): Indel-first assembly CIGAR for length-changing alts\n\nSkip wasted SoftClip SW when ref/alt lengths differ (dense extract\nphenotype); keep SoftClip fallback if Indel yields no I/D. Skip dangling\nmerge re-SW when the merge CIGAR already carries indels.\n\nCo-authored-by: Cursor <cursoragent@cursor.com>\n\n* perf(hc): cut loser kbest-strip, genotype EventMap, PairHMM pack, SW TLS\n\nDense CI losers spent wall in cycle-strip DFS, EventMap rebuilds, singleton\nPairHMM fallback, and mid-realign SW arena drops. Cache per-hap events,\nborrow graph outs + one remap, look-ahead equal-length packs with TLS\nscratch, and keep SW gap planes warm across realign.\n\nCo-authored-by: Cursor <cursoragent@cursor.com>\n\n---------\n\nCo-authored-by: Cursor <cursoragent@cursor.com>",
          "timestamp": "2026-08-14T03:31:21Z",
          "url": "https://github.com/SynapticFour/gatk-rs/commit/03445555bcea4bd70a4c75f31226b47ed90833ec"
        },
        "date": 1786680972835,
        "tool": "cargo",
        "benches": [
          {
            "name": "sam_parsing/parsing/100",
            "value": 98005,
            "range": "± 1677",
            "unit": "ns/iter"
          },
          {
            "name": "sam_parsing/parsing/1000",
            "value": 1129020,
            "range": "± 13514",
            "unit": "ns/iter"
          },
          {
            "name": "sam_parsing/parsing/10000",
            "value": 11198093,
            "range": "± 50573",
            "unit": "ns/iter"
          },
          {
            "name": "sam_iterator/iterator",
            "value": 821289,
            "range": "± 2729",
            "unit": "ns/iter"
          },
          {
            "name": "sam_writing/writing",
            "value": 1932239,
            "range": "± 23290",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/0",
            "value": 113,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/1",
            "value": 168,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/2",
            "value": 122,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/3",
            "value": 129,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/4",
            "value": 169,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "optional_fields/parsing",
            "value": 3981,
            "range": "± 17",
            "unit": "ns/iter"
          },
          {
            "name": "framework_overhead/create_result",
            "value": 22,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "framework_overhead/create_comparison",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "metrics_collection/metrics_collector",
            "value": 1061915,
            "range": "± 929",
            "unit": "ns/iter"
          },
          {
            "name": "suite_operations/create_suite",
            "value": 12979772,
            "range": "± 399924",
            "unit": "ns/iter"
          },
          {
            "name": "suite_operations/suite_with_results",
            "value": 13240993,
            "range": "± 360039",
            "unit": "ns/iter"
          },
          {
            "name": "suite_operations/suite_serialization",
            "value": 61670,
            "range": "± 1876",
            "unit": "ns/iter"
          },
          {
            "name": "performance_analysis/analyze_small_dataset",
            "value": 86,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "performance_analysis/analyze_large_dataset",
            "value": 10618,
            "range": "± 11",
            "unit": "ns/iter"
          },
          {
            "name": "regression_detection/detect_no_regression",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "regression_detection/detect_regression",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "report_generation/generate_json_report",
            "value": 185426,
            "range": "± 29606",
            "unit": "ns/iter"
          },
          {
            "name": "report_generation/generate_markdown_report",
            "value": 123500,
            "range": "± 7689",
            "unit": "ns/iter"
          },
          {
            "name": "report_generation/generate_csv_report",
            "value": 86120,
            "range": "± 6345",
            "unit": "ns/iter"
          },
          {
            "name": "dataset_operations/create_dataset_info",
            "value": 62,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "dataset_operations/validate_dataset_info",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "comparison_operations/create_comparison",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "comparison_operations/check_targets",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_gc_content",
            "value": 763274,
            "range": "± 14301",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_revcomp",
            "value": 348565,
            "range": "± 2421",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_subsequence",
            "value": 104,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "memory_pool_allocation",
            "value": 60,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "cache_put_get",
            "value": 91776,
            "range": "± 748",
            "unit": "ns/iter"
          },
          {
            "name": "interval_query",
            "value": 33939331,
            "range": "± 586836",
            "unit": "ns/iter"
          },
          {
            "name": "variant_type",
            "value": 13294,
            "range": "± 123",
            "unit": "ns/iter"
          },
          {
            "name": "quality_error_prob",
            "value": 1339,
            "range": "± 3",
            "unit": "ns/iter"
          },
          {
            "name": "log_addition",
            "value": 0,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "stream_processing",
            "value": 416,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "memory_mapped_access",
            "value": 3930,
            "range": "± 175",
            "unit": "ns/iter"
          },
          {
            "name": "hamming_distance",
            "value": 7,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/100",
            "value": 25905,
            "range": "± 299",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/100",
            "value": 48589,
            "range": "± 713",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/100",
            "value": 21609,
            "range": "± 95",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/1000",
            "value": 248105,
            "range": "± 1063",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/1000",
            "value": 387189,
            "range": "± 1217",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/1000",
            "value": 198459,
            "range": "± 810",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/10000",
            "value": 2560961,
            "range": "± 17880",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/10000",
            "value": 3863258,
            "range": "± 16297",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/10000",
            "value": 1996786,
            "range": "± 39098",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/100",
            "value": 27907,
            "range": "± 321",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/100",
            "value": 35072,
            "range": "± 310",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/100",
            "value": 19356,
            "range": "± 220",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/1000",
            "value": 246853,
            "range": "± 7026",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/1000",
            "value": 260773,
            "range": "± 743",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/1000",
            "value": 171753,
            "range": "± 2707",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/10000",
            "value": 2557029,
            "range": "± 37274",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/10000",
            "value": 2620084,
            "range": "± 9746",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/10000",
            "value": 1747848,
            "range": "± 11712",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/gc_content",
            "value": 55,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/reverse_complement",
            "value": 55,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/subsequence",
            "value": 17,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/count_pattern",
            "value": 120,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/average_quality",
            "value": 49,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/min_quality",
            "value": 50,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/max_quality",
            "value": 63,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/median_quality",
            "value": 38,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/trim_quality",
            "value": 43,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/reverse_complement",
            "value": 82,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/is_valid",
            "value": 80,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/build_index",
            "value": 332114,
            "range": "± 1228",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/random_access",
            "value": 2087,
            "range": "± 4",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/save_index",
            "value": 229841,
            "range": "± 54680",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/load_index",
            "value": 263910,
            "range": "± 2524",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/filter_by_quality",
            "value": 308927,
            "range": "± 806",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/filter_by_length",
            "value": 167375,
            "range": "± 1135",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/sample_reads",
            "value": 167079,
            "range": "± 1381",
            "unit": "ns/iter"
          },
          {
            "name": "io_writing/write_fasta",
            "value": 266163,
            "range": "± 16355",
            "unit": "ns/iter"
          },
          {
            "name": "io_writing/write_fastq",
            "value": 3789792,
            "range": "± 107450",
            "unit": "ns/iter"
          },
          {
            "name": "memory_usage/parse_memory_usage",
            "value": 2579984,
            "range": "± 14456",
            "unit": "ns/iter"
          },
          {
            "name": "macro_sum_bench",
            "value": 1,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/100",
            "value": 82611,
            "range": "± 2074",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/1000",
            "value": 914941,
            "range": "± 24548",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/10000",
            "value": 9936689,
            "range": "± 256406",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_iterator/iterator",
            "value": 603402,
            "range": "± 1614",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_writing/writing",
            "value": 2260142,
            "range": "± 16589",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/0",
            "value": 83,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/1",
            "value": 79,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/2",
            "value": 78,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/3",
            "value": 104,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/4",
            "value": 91,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/5",
            "value": 47,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/6",
            "value": 95,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/7",
            "value": 541,
            "range": "± 20",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/0",
            "value": 75,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/1",
            "value": 76,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/2",
            "value": 35,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/3",
            "value": 75,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/4",
            "value": 76,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/5",
            "value": 75,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/6",
            "value": 76,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/7",
            "value": 77,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/variant_type_detection",
            "value": 2260,
            "range": "± 54",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/allele_frequency_access",
            "value": 3763,
            "range": "± 5",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/sample_data_access",
            "value": 1568,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "band_pass_normalized_kernel_gatk_defaults",
            "value": 1022,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "band_pass_add_then_pop_ready_256_loci",
            "value": 281061,
            "range": "± 1958",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_medium_depth_k10",
            "value": 227998,
            "range": "± 2051",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_high_depth_k10",
            "value": 1702539,
            "range": "± 10202",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/4",
            "value": 14518437,
            "range": "± 926703",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/16",
            "value": 52626562,
            "range": "± 1559582",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/64",
            "value": 198727041,
            "range": "± 15336325",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/8",
            "value": 432927,
            "range": "± 573",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/8",
            "value": 194796,
            "range": "± 1626",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/8",
            "value": 895173,
            "range": "± 1761",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/8",
            "value": 9812685,
            "range": "± 15428",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/32",
            "value": 1731337,
            "range": "± 5626",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/32",
            "value": 768034,
            "range": "± 27324",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/32",
            "value": 3549282,
            "range": "± 5604",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/32",
            "value": 39273954,
            "range": "± 81902",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/64",
            "value": 3462547,
            "range": "± 6422",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/64",
            "value": 1528747,
            "range": "± 6177",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/64",
            "value": 7088401,
            "range": "± 18927",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/64",
            "value": 78565370,
            "range": "± 125205",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/128",
            "value": 6924115,
            "range": "± 7382",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/128",
            "value": 3043069,
            "range": "± 8287",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/128",
            "value": 14162944,
            "range": "± 24188",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/128",
            "value": 157111525,
            "range": "± 242457",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/8",
            "value": 1660058,
            "range": "± 6382",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/8",
            "value": 3429425,
            "range": "± 29164",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/8",
            "value": 3719998,
            "range": "± 6247",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/8",
            "value": 39721037,
            "range": "± 138348",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/32",
            "value": 6637184,
            "range": "± 26232",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/32",
            "value": 13206680,
            "range": "± 129678",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/32",
            "value": 14793468,
            "range": "± 22244",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/32",
            "value": 158865503,
            "range": "± 356971",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/64",
            "value": 13268823,
            "range": "± 52158",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/64",
            "value": 26343541,
            "range": "± 197852",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/64",
            "value": 29573678,
            "range": "± 38751",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/64",
            "value": 317948157,
            "range": "± 541311",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/128",
            "value": 26538455,
            "range": "± 58551",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/128",
            "value": 53281894,
            "range": "± 620870",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/128",
            "value": 59136887,
            "range": "± 760421",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/128",
            "value": 635626573,
            "range": "± 920549",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/8",
            "value": 4074745,
            "range": "± 140253",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/8",
            "value": 6358715,
            "range": "± 284244",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/8",
            "value": 8480523,
            "range": "± 32200",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/8",
            "value": 90159625,
            "range": "± 193129",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/32",
            "value": 16293078,
            "range": "± 28795",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/32",
            "value": 24043527,
            "range": "± 1401354",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/32",
            "value": 33723611,
            "range": "± 29565",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/32",
            "value": 359819777,
            "range": "± 960571",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/64",
            "value": 32593267,
            "range": "± 111992",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/64",
            "value": 47701002,
            "range": "± 1362825",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/64",
            "value": 67378684,
            "range": "± 128402",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/64",
            "value": 720515125,
            "range": "± 2182259",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/128",
            "value": 65191535,
            "range": "± 85769",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/128",
            "value": 97334444,
            "range": "± 3166975",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/128",
            "value": 134007947,
            "range": "± 305312",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/128",
            "value": 1441007871,
            "range": "± 2888716",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/value_copy_stages/3072x12",
            "value": 4094,
            "range": "± 22",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/arc_share_stages/3072x12",
            "value": 3207,
            "range": "± 6",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/assembly_result_set_arc_fanout/3072",
            "value": 177,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/assembly_result_set_value_fanout/3072",
            "value": 1162,
            "range": "± 11",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/64x48",
            "value": 16663,
            "range": "± 920",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/64x48",
            "value": 16415,
            "range": "± 240",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/128x96",
            "value": 67133,
            "range": "± 257",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/128x96",
            "value": 67519,
            "range": "± 702",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/256x192",
            "value": 260128,
            "range": "± 551",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/256x192",
            "value": 263268,
            "range": "± 801",
            "unit": "ns/iter"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "name": "SynapticFour",
            "username": "SynapticFour",
            "email": "contact@synapticfour.com"
          },
          "committer": {
            "name": "GitHub",
            "username": "web-flow",
            "email": "noreply@github.com"
          },
          "id": "5db42deab92ab38bf379022520dbabdfb11283c3",
          "message": "perf(giab): honest Peak probe + phase6 wall rematch baseline (#112)\n\nFix hc-mem-probe to track max HC RSS (skip docker wrappers) and sample on\nexit; prefer probe Peak in the ci-subset summarizer. Document cancelled\nphase5 CI median (~1.57×), TRACE share rematch, holdout F1 gap, and the\nnext beat-Java wall bets without speculative P12-adjacent cuts.\n\nCo-authored-by: Cursor <cursoragent@cursor.com>",
          "timestamp": "2026-08-14T21:58:55Z",
          "url": "https://github.com/SynapticFour/gatk-rs/commit/5db42deab92ab38bf379022520dbabdfb11283c3"
        },
        "date": 1786763948388,
        "tool": "cargo",
        "benches": [
          {
            "name": "sam_parsing/parsing/100",
            "value": 110485,
            "range": "± 902",
            "unit": "ns/iter"
          },
          {
            "name": "sam_parsing/parsing/1000",
            "value": 1119426,
            "range": "± 6646",
            "unit": "ns/iter"
          },
          {
            "name": "sam_parsing/parsing/10000",
            "value": 10998360,
            "range": "± 149183",
            "unit": "ns/iter"
          },
          {
            "name": "sam_iterator/iterator",
            "value": 851781,
            "range": "± 27446",
            "unit": "ns/iter"
          },
          {
            "name": "sam_writing/writing",
            "value": 3317122,
            "range": "± 446462",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/0",
            "value": 131,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/1",
            "value": 174,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/2",
            "value": 137,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/3",
            "value": 142,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/4",
            "value": 174,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "optional_fields/parsing",
            "value": 7062,
            "range": "± 169",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_gc_content",
            "value": 814844,
            "range": "± 6460",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_revcomp",
            "value": 297886,
            "range": "± 565",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_subsequence",
            "value": 165,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "variant_type",
            "value": 15850,
            "range": "± 34",
            "unit": "ns/iter"
          },
          {
            "name": "quality_error_prob",
            "value": 1216,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "memory_mapped_access",
            "value": 9233,
            "range": "± 129",
            "unit": "ns/iter"
          },
          {
            "name": "hamming_distance",
            "value": 9,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/100",
            "value": 29321,
            "range": "± 686",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/100",
            "value": 54085,
            "range": "± 258",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/100",
            "value": 26050,
            "range": "± 346",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/1000",
            "value": 242141,
            "range": "± 1202",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/1000",
            "value": 409351,
            "range": "± 1958",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/1000",
            "value": 206759,
            "range": "± 2158",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/10000",
            "value": 2431072,
            "range": "± 19327",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/10000",
            "value": 3793224,
            "range": "± 23043",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/10000",
            "value": 2066588,
            "range": "± 33539",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/100",
            "value": 32013,
            "range": "± 73",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/100",
            "value": 41009,
            "range": "± 188",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/100",
            "value": 25165,
            "range": "± 257",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/1000",
            "value": 267572,
            "range": "± 1807",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/1000",
            "value": 274003,
            "range": "± 1341",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/1000",
            "value": 197618,
            "range": "± 960",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/10000",
            "value": 2600039,
            "range": "± 37920",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/10000",
            "value": 2572450,
            "range": "± 36064",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/10000",
            "value": 1967146,
            "range": "± 15802",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/gc_content",
            "value": 63,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/reverse_complement",
            "value": 60,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/subsequence",
            "value": 23,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/count_pattern",
            "value": 209,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/average_quality",
            "value": 62,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/min_quality",
            "value": 44,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/max_quality",
            "value": 49,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/median_quality",
            "value": 67,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/trim_quality",
            "value": 54,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/reverse_complement",
            "value": 92,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/is_valid",
            "value": 92,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/build_index",
            "value": 347688,
            "range": "± 1178",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/random_access",
            "value": 4934,
            "range": "± 66",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/save_index",
            "value": 222874,
            "range": "± 15235",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/load_index",
            "value": 312808,
            "range": "± 2566",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/filter_by_quality",
            "value": 329426,
            "range": "± 1156",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/filter_by_length",
            "value": 183689,
            "range": "± 623",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/sample_reads",
            "value": 184158,
            "range": "± 768",
            "unit": "ns/iter"
          },
          {
            "name": "io_writing/write_fasta",
            "value": 208943,
            "range": "± 12682",
            "unit": "ns/iter"
          },
          {
            "name": "io_writing/write_fastq",
            "value": 8977418,
            "range": "± 311134",
            "unit": "ns/iter"
          },
          {
            "name": "memory_usage/parse_memory_usage",
            "value": 2435998,
            "range": "± 196183",
            "unit": "ns/iter"
          },
          {
            "name": "macro_sum_bench",
            "value": 1,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/100",
            "value": 92436,
            "range": "± 941",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/1000",
            "value": 876465,
            "range": "± 6801",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/10000",
            "value": 9429022,
            "range": "± 217066",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_iterator/iterator",
            "value": 632740,
            "range": "± 7445",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_writing/writing",
            "value": 3956253,
            "range": "± 33027",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/0",
            "value": 113,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/1",
            "value": 104,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/2",
            "value": 101,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/3",
            "value": 131,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/4",
            "value": 117,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/5",
            "value": 63,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/6",
            "value": 120,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/7",
            "value": 481,
            "range": "± 3",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/0",
            "value": 91,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/1",
            "value": 91,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/2",
            "value": 47,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/3",
            "value": 91,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/4",
            "value": 91,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/5",
            "value": 91,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/6",
            "value": 91,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/7",
            "value": 91,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/variant_type_detection",
            "value": 2472,
            "range": "± 8",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/allele_frequency_access",
            "value": 4234,
            "range": "± 8",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/sample_data_access",
            "value": 1518,
            "range": "± 5",
            "unit": "ns/iter"
          },
          {
            "name": "band_pass_normalized_kernel_gatk_defaults",
            "value": 995,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "band_pass_add_then_pop_ready_256_loci",
            "value": 329438,
            "range": "± 722",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_medium_depth_k10",
            "value": 255961,
            "range": "± 1029",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_high_depth_k10",
            "value": 1836814,
            "range": "± 3503",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/4",
            "value": 24877059,
            "range": "± 227694",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/16",
            "value": 99512842,
            "range": "± 176499",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/64",
            "value": 398154691,
            "range": "± 1552830",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/8",
            "value": 548024,
            "range": "± 885",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/8",
            "value": 272055,
            "range": "± 6884",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/8",
            "value": 1037678,
            "range": "± 11173",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/8",
            "value": 9018792,
            "range": "± 11557",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/32",
            "value": 2191686,
            "range": "± 5226",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/32",
            "value": 1071073,
            "range": "± 14233",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/32",
            "value": 4131250,
            "range": "± 3911",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/32",
            "value": 36098684,
            "range": "± 649201",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/64",
            "value": 4384063,
            "range": "± 6209",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/64",
            "value": 2131597,
            "range": "± 7480",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/64",
            "value": 8255303,
            "range": "± 8571",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/64",
            "value": 72157847,
            "range": "± 153708",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/128",
            "value": 8768724,
            "range": "± 23881",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/128",
            "value": 4258787,
            "range": "± 15151",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/128",
            "value": 16494432,
            "range": "± 17164",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/128",
            "value": 144278967,
            "range": "± 134984",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/8",
            "value": 2138637,
            "range": "± 3830",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/8",
            "value": 1035097,
            "range": "± 17888",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/8",
            "value": 4106361,
            "range": "± 81994",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/8",
            "value": 36433800,
            "range": "± 52680",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/32",
            "value": 8544771,
            "range": "± 13026",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/32",
            "value": 4141423,
            "range": "± 12442",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/32",
            "value": 16390317,
            "range": "± 295732",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/32",
            "value": 145740968,
            "range": "± 234618",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/64",
            "value": 17087935,
            "range": "± 24032",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/64",
            "value": 8263703,
            "range": "± 19493",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/64",
            "value": 32777391,
            "range": "± 594915",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/64",
            "value": 291563745,
            "range": "± 482914",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/128",
            "value": 34166012,
            "range": "± 62439",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/128",
            "value": 16525734,
            "range": "± 42679",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/128",
            "value": 65521498,
            "range": "± 85061",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/128",
            "value": 583229450,
            "range": "± 544002",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/8",
            "value": 4731952,
            "range": "± 6282",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/8",
            "value": 2508490,
            "range": "± 9351",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/8",
            "value": 9115263,
            "range": "± 10925",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/8",
            "value": 82283273,
            "range": "± 193314",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/32",
            "value": 18931434,
            "range": "± 34340",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/32",
            "value": 10000724,
            "range": "± 78129",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/32",
            "value": 36411201,
            "range": "± 39688",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/32",
            "value": 329234696,
            "range": "± 333581",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/64",
            "value": 37866929,
            "range": "± 36973",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/64",
            "value": 19987071,
            "range": "± 95938",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/64",
            "value": 72826346,
            "range": "± 113991",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/64",
            "value": 658436953,
            "range": "± 2567261",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/128",
            "value": 76027747,
            "range": "± 127832",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/128",
            "value": 39893466,
            "range": "± 136289",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/128",
            "value": 145921342,
            "range": "± 92822",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/128",
            "value": 1316750497,
            "range": "± 1889450",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/value_copy_stages/3072x12",
            "value": 5113,
            "range": "± 37",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/arc_share_stages/3072x12",
            "value": 3923,
            "range": "± 3",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/assembly_result_set_arc_fanout/3072",
            "value": 56,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/assembly_result_set_value_fanout/3072",
            "value": 1204,
            "range": "± 6",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/64x48",
            "value": 17342,
            "range": "± 219",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/64x48",
            "value": 17218,
            "range": "± 91",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/128x96",
            "value": 67968,
            "range": "± 291",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/128x96",
            "value": 67639,
            "range": "± 154",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/256x192",
            "value": 260865,
            "range": "± 20975",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/256x192",
            "value": 260019,
            "range": "± 602",
            "unit": "ns/iter"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "name": "SynapticFour",
            "username": "SynapticFour",
            "email": "contact@synapticfour.com"
          },
          "committer": {
            "name": "GitHub",
            "username": "web-flow",
            "email": "noreply@github.com"
          },
          "id": "24a30c1c967ef26261e417b1db2485c755ea2c72",
          "message": "perf(hc): phase9 NEON by_len PairHMM packs + TRACE occupancy (#118)\n\nPhase9: NEON by_len hap packs + TRACE neon_pack2/leftover occupancy. Same Logless numerics; hap leftover already ~6.7% so next leaf is read-axis packs.",
          "timestamp": "2026-08-15T21:23:47Z",
          "url": "https://github.com/SynapticFour/gatk-rs/commit/24a30c1c967ef26261e417b1db2485c755ea2c72"
        },
        "date": 1786850610073,
        "tool": "cargo",
        "benches": [
          {
            "name": "sam_parsing/parsing/100",
            "value": 119087,
            "range": "± 4983",
            "unit": "ns/iter"
          },
          {
            "name": "sam_parsing/parsing/1000",
            "value": 1206884,
            "range": "± 19973",
            "unit": "ns/iter"
          },
          {
            "name": "sam_parsing/parsing/10000",
            "value": 11787513,
            "range": "± 47218",
            "unit": "ns/iter"
          },
          {
            "name": "sam_iterator/iterator",
            "value": 915582,
            "range": "± 3771",
            "unit": "ns/iter"
          },
          {
            "name": "sam_writing/writing",
            "value": 3215942,
            "range": "± 52613",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/0",
            "value": 124,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/1",
            "value": 178,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/2",
            "value": 133,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/3",
            "value": 140,
            "range": "± 6",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/4",
            "value": 186,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "optional_fields/parsing",
            "value": 6634,
            "range": "± 23",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_gc_content",
            "value": 715073,
            "range": "± 4253",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_revcomp",
            "value": 258732,
            "range": "± 540",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_subsequence",
            "value": 147,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "variant_type",
            "value": 16536,
            "range": "± 442",
            "unit": "ns/iter"
          },
          {
            "name": "quality_error_prob",
            "value": 1279,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "memory_mapped_access",
            "value": 8473,
            "range": "± 120",
            "unit": "ns/iter"
          },
          {
            "name": "hamming_distance",
            "value": 8,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/100",
            "value": 33831,
            "range": "± 166",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/100",
            "value": 56553,
            "range": "± 487",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/100",
            "value": 25353,
            "range": "± 76",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/1000",
            "value": 298660,
            "range": "± 1029",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/1000",
            "value": 450302,
            "range": "± 1664",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/1000",
            "value": 216717,
            "range": "± 2007",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/10000",
            "value": 2975460,
            "range": "± 34288",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/10000",
            "value": 4089956,
            "range": "± 10147",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/10000",
            "value": 2183527,
            "range": "± 26862",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/100",
            "value": 38196,
            "range": "± 147",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/100",
            "value": 49585,
            "range": "± 264",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/100",
            "value": 23696,
            "range": "± 62",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/1000",
            "value": 328648,
            "range": "± 1362",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/1000",
            "value": 345677,
            "range": "± 1328",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/1000",
            "value": 193304,
            "range": "± 533",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/10000",
            "value": 3111052,
            "range": "± 21723",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/10000",
            "value": 3205330,
            "range": "± 27448",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/10000",
            "value": 1888781,
            "range": "± 28393",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/gc_content",
            "value": 57,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/reverse_complement",
            "value": 56,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/subsequence",
            "value": 20,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/count_pattern",
            "value": 204,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/average_quality",
            "value": 54,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/min_quality",
            "value": 40,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/max_quality",
            "value": 41,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/median_quality",
            "value": 66,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/trim_quality",
            "value": 50,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/reverse_complement",
            "value": 89,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/is_valid",
            "value": 83,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/build_index",
            "value": 369073,
            "range": "± 1719",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/random_access",
            "value": 4394,
            "range": "± 35",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/save_index",
            "value": 286896,
            "range": "± 16144",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/load_index",
            "value": 302984,
            "range": "± 675",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/filter_by_quality",
            "value": 366211,
            "range": "± 12109",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/filter_by_length",
            "value": 181334,
            "range": "± 2527",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/sample_reads",
            "value": 180621,
            "range": "± 1069",
            "unit": "ns/iter"
          },
          {
            "name": "io_writing/write_fasta",
            "value": 357287,
            "range": "± 17242",
            "unit": "ns/iter"
          },
          {
            "name": "io_writing/write_fastq",
            "value": 8366257,
            "range": "± 67131",
            "unit": "ns/iter"
          },
          {
            "name": "memory_usage/parse_memory_usage",
            "value": 3001825,
            "range": "± 112074",
            "unit": "ns/iter"
          },
          {
            "name": "macro_sum_bench",
            "value": 2,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/100",
            "value": 92298,
            "range": "± 1342",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/1000",
            "value": 880342,
            "range": "± 5787",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/10000",
            "value": 9574351,
            "range": "± 88926",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_iterator/iterator",
            "value": 649104,
            "range": "± 1778",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_writing/writing",
            "value": 3603485,
            "range": "± 57029",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/0",
            "value": 98,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/1",
            "value": 93,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/2",
            "value": 90,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/3",
            "value": 124,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/4",
            "value": 105,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/5",
            "value": 56,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/6",
            "value": 108,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/7",
            "value": 479,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/0",
            "value": 85,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/1",
            "value": 85,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/2",
            "value": 41,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/3",
            "value": 86,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/4",
            "value": 85,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/5",
            "value": 85,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/6",
            "value": 85,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/7",
            "value": 85,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/variant_type_detection",
            "value": 2571,
            "range": "± 13",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/allele_frequency_access",
            "value": 4084,
            "range": "± 33",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/sample_data_access",
            "value": 1380,
            "range": "± 3",
            "unit": "ns/iter"
          },
          {
            "name": "band_pass_normalized_kernel_gatk_defaults",
            "value": 1000,
            "range": "± 11",
            "unit": "ns/iter"
          },
          {
            "name": "band_pass_add_then_pop_ready_256_loci",
            "value": 284890,
            "range": "± 434",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_medium_depth_k10",
            "value": 242136,
            "range": "± 408",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_high_depth_k10",
            "value": 1752091,
            "range": "± 5316",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/4",
            "value": 23856619,
            "range": "± 692580",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/16",
            "value": 95117322,
            "range": "± 2074925",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/64",
            "value": 380436919,
            "range": "± 2270684",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/8",
            "value": 504394,
            "range": "± 4171",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/8",
            "value": 267275,
            "range": "± 3453",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/8",
            "value": 961453,
            "range": "± 2604",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/8",
            "value": 9448520,
            "range": "± 12572",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/32",
            "value": 2012927,
            "range": "± 3975",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/32",
            "value": 1017032,
            "range": "± 15729",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/32",
            "value": 3817492,
            "range": "± 6517",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/32",
            "value": 37805448,
            "range": "± 93354",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/64",
            "value": 4024675,
            "range": "± 11444",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/64",
            "value": 2018509,
            "range": "± 14320",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/64",
            "value": 7630571,
            "range": "± 12715",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/64",
            "value": 75582119,
            "range": "± 91603",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/128",
            "value": 8053446,
            "range": "± 14780",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/128",
            "value": 4186343,
            "range": "± 11151",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/128",
            "value": 15285128,
            "range": "± 19400",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/128",
            "value": 151225676,
            "range": "± 446672",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/8",
            "value": 1977146,
            "range": "± 8536",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/8",
            "value": 1050646,
            "range": "± 18888",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/8",
            "value": 3780567,
            "range": "± 10276",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/8",
            "value": 38167059,
            "range": "± 375281",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/32",
            "value": 7755745,
            "range": "± 28575",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/32",
            "value": 4104346,
            "range": "± 15717",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/32",
            "value": 15017454,
            "range": "± 30188",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/32",
            "value": 152776783,
            "range": "± 212778",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/64",
            "value": 15512629,
            "range": "± 42841",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/64",
            "value": 8190579,
            "range": "± 18569",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/64",
            "value": 30021622,
            "range": "± 96040",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/64",
            "value": 305741876,
            "range": "± 2005510",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/128",
            "value": 31449259,
            "range": "± 64758",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/128",
            "value": 15758548,
            "range": "± 59386",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/128",
            "value": 59720383,
            "range": "± 179996",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/128",
            "value": 611515867,
            "range": "± 1512241",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/8",
            "value": 4313464,
            "range": "± 11942",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/8",
            "value": 2311102,
            "range": "± 16948",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/8",
            "value": 8284080,
            "range": "± 22347",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/8",
            "value": 86168120,
            "range": "± 134045",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/32",
            "value": 17337441,
            "range": "± 66569",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/32",
            "value": 9488378,
            "range": "± 79167",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/32",
            "value": 33458672,
            "range": "± 66179",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/32",
            "value": 345014477,
            "range": "± 340380",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/64",
            "value": 34973002,
            "range": "± 196600",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/64",
            "value": 18366363,
            "range": "± 63849",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/64",
            "value": 66784899,
            "range": "± 146077",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/64",
            "value": 690032368,
            "range": "± 603066",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/128",
            "value": 69617478,
            "range": "± 172892",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/128",
            "value": 37606655,
            "range": "± 205186",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/128",
            "value": 132973449,
            "range": "± 202678",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/128",
            "value": 1380436040,
            "range": "± 2748054",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/value_copy_stages/3072x12",
            "value": 4896,
            "range": "± 90",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/arc_share_stages/3072x12",
            "value": 3471,
            "range": "± 9",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/assembly_result_set_arc_fanout/3072",
            "value": 53,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/assembly_result_set_value_fanout/3072",
            "value": 1221,
            "range": "± 5",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/64x48",
            "value": 16695,
            "range": "± 32",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/64x48",
            "value": 16195,
            "range": "± 247",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/128x96",
            "value": 61844,
            "range": "± 189",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/128x96",
            "value": 60072,
            "range": "± 315",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/256x192",
            "value": 237877,
            "range": "± 477",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/256x192",
            "value": 231829,
            "range": "± 515",
            "unit": "ns/iter"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "name": "SynapticFour",
            "username": "SynapticFour",
            "email": "contact@synapticfour.com"
          },
          "committer": {
            "name": "GitHub",
            "username": "web-flow",
            "email": "noreply@github.com"
          },
          "id": "6a478ca1845f39ee15550390a1fe1994e9ebbfac",
          "message": "perf(hc): collapse multi-pass pileup AD and cut wall leaves (#120)\n\n* perf(hc): collapse multi-pass pileup AD and cut wall leaves\n\nReuse CIGAR/seq decode and equal-pad AD scans in try_genotype, fix NEON\nprefix-reuse TRACE, prefer longer hap prefixes / AVX2 reuse, skip SNP\nspine with no alt haps, and tighten SoftClip last_index_of.\n\nCo-authored-by: Cursor <cursoragent@cursor.com>\n\n* perf(hc): further AD softclip/pad reuse and NEON PairHMM TLS\n\nCollapse remaining softclip/trim AD double-scans and equivalent-pad\nrescans; reuse AD decode cache on alignment/anchor paths. NEON: TLS\nby_len map, in-place hap sort, leftover singles via score_one_hap.\n\nCo-authored-by: Cursor <cursoragent@cursor.com>\n\n* docs(parity-site): GIAB wall-losers history (7991f2817e78)\n\n* perf(hc): softclip AD TLS + AVX2 PairHMM hygiene; wall/callrate docs\n\nExtend AdDecodeCache with SoftClip-as-ref base lookup for softclip/alignment\nscanners; mirror NEON TLS by_len and leftover score_one_hap on AVX2. Record\nphase10 wall-losers baseline and stage-classify remaining Java-only undercall.\n\nCo-authored-by: Cursor <cursoragent@cursor.com>\n\n* fix(hc): restore SNP parity spine when assembly is ref-only\n\nThe no-alt early return blocked materialize_alt_haps from creating SNP\nalt haplotypes, so call_region returned None on the p5 j2 emit fixture.\n\nCo-authored-by: Cursor <cursoragent@cursor.com>\n\n---------\n\nCo-authored-by: Cursor <cursoragent@cursor.com>\nCo-authored-by: github-actions[bot] <41898282+github-actions[bot]@users.noreply.github.com>",
          "timestamp": "2026-08-16T15:37:53Z",
          "url": "https://github.com/SynapticFour/gatk-rs/commit/6a478ca1845f39ee15550390a1fe1994e9ebbfac"
        },
        "date": 1786937015214,
        "tool": "cargo",
        "benches": [
          {
            "name": "sam_parsing/parsing/100",
            "value": 110761,
            "range": "± 2234",
            "unit": "ns/iter"
          },
          {
            "name": "sam_parsing/parsing/1000",
            "value": 1109917,
            "range": "± 6528",
            "unit": "ns/iter"
          },
          {
            "name": "sam_parsing/parsing/10000",
            "value": 11014365,
            "range": "± 78232",
            "unit": "ns/iter"
          },
          {
            "name": "sam_iterator/iterator",
            "value": 901898,
            "range": "± 16142",
            "unit": "ns/iter"
          },
          {
            "name": "sam_writing/writing",
            "value": 3352853,
            "range": "± 33430",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/0",
            "value": 129,
            "range": "± 3",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/1",
            "value": 172,
            "range": "± 5",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/2",
            "value": 135,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/3",
            "value": 142,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/4",
            "value": 174,
            "range": "± 3",
            "unit": "ns/iter"
          },
          {
            "name": "optional_fields/parsing",
            "value": 7220,
            "range": "± 73",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_gc_content",
            "value": 815977,
            "range": "± 19986",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_revcomp",
            "value": 297112,
            "range": "± 1821",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_subsequence",
            "value": 164,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "variant_type",
            "value": 15852,
            "range": "± 40",
            "unit": "ns/iter"
          },
          {
            "name": "quality_error_prob",
            "value": 1215,
            "range": "± 3",
            "unit": "ns/iter"
          },
          {
            "name": "memory_mapped_access",
            "value": 9319,
            "range": "± 433",
            "unit": "ns/iter"
          },
          {
            "name": "hamming_distance",
            "value": 9,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/100",
            "value": 28927,
            "range": "± 119",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/100",
            "value": 54079,
            "range": "± 843",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/100",
            "value": 25565,
            "range": "± 91",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/1000",
            "value": 242567,
            "range": "± 1311",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/1000",
            "value": 409006,
            "range": "± 2465",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/1000",
            "value": 203288,
            "range": "± 5359",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/10000",
            "value": 2423222,
            "range": "± 39137",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/10000",
            "value": 3772619,
            "range": "± 31186",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/10000",
            "value": 2070283,
            "range": "± 18027",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/100",
            "value": 33086,
            "range": "± 301",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/100",
            "value": 41281,
            "range": "± 604",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/100",
            "value": 25079,
            "range": "± 346",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/1000",
            "value": 271076,
            "range": "± 4220",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/1000",
            "value": 276271,
            "range": "± 1325",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/1000",
            "value": 201203,
            "range": "± 3161",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/10000",
            "value": 2716457,
            "range": "± 97090",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/10000",
            "value": 2569935,
            "range": "± 23552",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/10000",
            "value": 1968790,
            "range": "± 28304",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/gc_content",
            "value": 63,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/reverse_complement",
            "value": 61,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/subsequence",
            "value": 23,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/count_pattern",
            "value": 209,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/average_quality",
            "value": 62,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/min_quality",
            "value": 45,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/max_quality",
            "value": 49,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/median_quality",
            "value": 67,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/trim_quality",
            "value": 54,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/reverse_complement",
            "value": 93,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/is_valid",
            "value": 92,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/build_index",
            "value": 351787,
            "range": "± 4087",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/random_access",
            "value": 4918,
            "range": "± 23",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/save_index",
            "value": 239598,
            "range": "± 19822",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/load_index",
            "value": 314825,
            "range": "± 1317",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/filter_by_quality",
            "value": 336275,
            "range": "± 4227",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/filter_by_length",
            "value": 187159,
            "range": "± 2097",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/sample_reads",
            "value": 186796,
            "range": "± 6168",
            "unit": "ns/iter"
          },
          {
            "name": "io_writing/write_fasta",
            "value": 228614,
            "range": "± 26698",
            "unit": "ns/iter"
          },
          {
            "name": "io_writing/write_fastq",
            "value": 9007175,
            "range": "± 632024",
            "unit": "ns/iter"
          },
          {
            "name": "memory_usage/parse_memory_usage",
            "value": 2461042,
            "range": "± 25679",
            "unit": "ns/iter"
          },
          {
            "name": "macro_sum_bench",
            "value": 1,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/100",
            "value": 90860,
            "range": "± 587",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/1000",
            "value": 872946,
            "range": "± 4865",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/10000",
            "value": 9500357,
            "range": "± 154013",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_iterator/iterator",
            "value": 666850,
            "range": "± 6107",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_writing/writing",
            "value": 3963244,
            "range": "± 92021",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/0",
            "value": 105,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/1",
            "value": 96,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/2",
            "value": 93,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/3",
            "value": 130,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/4",
            "value": 110,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/5",
            "value": 60,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/6",
            "value": 114,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/7",
            "value": 479,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/0",
            "value": 92,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/1",
            "value": 92,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/2",
            "value": 47,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/3",
            "value": 93,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/4",
            "value": 93,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/5",
            "value": 92,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/6",
            "value": 93,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/7",
            "value": 92,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/variant_type_detection",
            "value": 2473,
            "range": "± 12",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/allele_frequency_access",
            "value": 4250,
            "range": "± 47",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/sample_data_access",
            "value": 1523,
            "range": "± 45",
            "unit": "ns/iter"
          },
          {
            "name": "band_pass_normalized_kernel_gatk_defaults",
            "value": 993,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "band_pass_add_then_pop_ready_256_loci",
            "value": 316394,
            "range": "± 851",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_medium_depth_k10",
            "value": 252478,
            "range": "± 1115",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_high_depth_k10",
            "value": 1836507,
            "range": "± 6355",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/4",
            "value": 24938622,
            "range": "± 336228",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/16",
            "value": 99470310,
            "range": "± 1094543",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/64",
            "value": 397780623,
            "range": "± 1411326",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/8",
            "value": 548653,
            "range": "± 5787",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/8",
            "value": 299599,
            "range": "± 295",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/8",
            "value": 1038227,
            "range": "± 8898",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/8",
            "value": 9127749,
            "range": "± 155005",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/32",
            "value": 2193513,
            "range": "± 45279",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/32",
            "value": 881152,
            "range": "± 6698",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/32",
            "value": 4133828,
            "range": "± 5046",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/32",
            "value": 36554708,
            "range": "± 686268",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/64",
            "value": 4387341,
            "range": "± 39872",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/64",
            "value": 1395889,
            "range": "± 8091",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/64",
            "value": 8260252,
            "range": "± 6830",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/64",
            "value": 72495350,
            "range": "± 1442390",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/128",
            "value": 8774376,
            "range": "± 114489",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/128",
            "value": 1983182,
            "range": "± 6390",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/128",
            "value": 16509266,
            "range": "± 23805",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/128",
            "value": 145925601,
            "range": "± 2969558",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/8",
            "value": 2136715,
            "range": "± 9706",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/8",
            "value": 1205530,
            "range": "± 6015",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/8",
            "value": 4086774,
            "range": "± 5480",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/8",
            "value": 36901442,
            "range": "± 779668",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/32",
            "value": 8540612,
            "range": "± 9153",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/32",
            "value": 3792441,
            "range": "± 14973",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/32",
            "value": 16313899,
            "range": "± 79429",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/32",
            "value": 145974821,
            "range": "± 2776405",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/64",
            "value": 17083062,
            "range": "± 121739",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/64",
            "value": 6712517,
            "range": "± 79088",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/64",
            "value": 32628289,
            "range": "± 123387",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/64",
            "value": 297734643,
            "range": "± 5940762",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/128",
            "value": 34205512,
            "range": "± 224392",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/128",
            "value": 10763382,
            "range": "± 40076",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/128",
            "value": 65510356,
            "range": "± 797677",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/128",
            "value": 585146058,
            "range": "± 11892184",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/8",
            "value": 4734245,
            "range": "± 5972",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/8",
            "value": 2709103,
            "range": "± 32634",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/8",
            "value": 9097116,
            "range": "± 9123",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/8",
            "value": 84153340,
            "range": "± 1408070",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/32",
            "value": 18932838,
            "range": "± 32761",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/32",
            "value": 8722355,
            "range": "± 31693",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/32",
            "value": 36462131,
            "range": "± 287209",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/32",
            "value": 334027338,
            "range": "± 6374963",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/64",
            "value": 38019337,
            "range": "± 48103",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/64",
            "value": 15953026,
            "range": "± 195702",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/64",
            "value": 72911874,
            "range": "± 301232",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/64",
            "value": 669447577,
            "range": "± 12553537",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/128",
            "value": 76009217,
            "range": "± 87463",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/128",
            "value": 27723029,
            "range": "± 372679",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/128",
            "value": 145815827,
            "range": "± 1677684",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/128",
            "value": 1341519014,
            "range": "± 25361678",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/value_copy_stages/3072x12",
            "value": 5096,
            "range": "± 47",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/arc_share_stages/3072x12",
            "value": 3923,
            "range": "± 43",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/assembly_result_set_arc_fanout/3072",
            "value": 56,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/assembly_result_set_value_fanout/3072",
            "value": 1203,
            "range": "± 12",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/64x48",
            "value": 16515,
            "range": "± 206",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/64x48",
            "value": 16525,
            "range": "± 35",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/128x96",
            "value": 67982,
            "range": "± 282",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/128x96",
            "value": 67875,
            "range": "± 146",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/256x192",
            "value": 252812,
            "range": "± 847",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/256x192",
            "value": 252381,
            "range": "± 875",
            "unit": "ns/iter"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "name": "SynapticFour",
            "username": "SynapticFour",
            "email": "contact@synapticfour.com"
          },
          "committer": {
            "name": "GitHub",
            "username": "web-flow",
            "email": "noreply@github.com"
          },
          "id": "67df80140b9a5c5a2c2ff42238adc54a5c5e8b2e",
          "message": "ci: stop parity-site jobs pushing protected main (#124)\n\nPublish history.json via Pages instead of git push to protected main. Skip paid genomewide VM when CLOUD_PROVIDER is unset.",
          "timestamp": "2026-08-17T17:58:17Z",
          "url": "https://github.com/SynapticFour/gatk-rs/commit/67df80140b9a5c5a2c2ff42238adc54a5c5e8b2e"
        },
        "date": 1787023215247,
        "tool": "cargo",
        "benches": [
          {
            "name": "sam_parsing/parsing/100",
            "value": 111882,
            "range": "± 1395",
            "unit": "ns/iter"
          },
          {
            "name": "sam_parsing/parsing/1000",
            "value": 1138393,
            "range": "± 23542",
            "unit": "ns/iter"
          },
          {
            "name": "sam_parsing/parsing/10000",
            "value": 11146019,
            "range": "± 374322",
            "unit": "ns/iter"
          },
          {
            "name": "sam_iterator/iterator",
            "value": 865850,
            "range": "± 14337",
            "unit": "ns/iter"
          },
          {
            "name": "sam_writing/writing",
            "value": 3362116,
            "range": "± 164446",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/0",
            "value": 128,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/1",
            "value": 174,
            "range": "± 13",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/2",
            "value": 136,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/3",
            "value": 145,
            "range": "± 3",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/4",
            "value": 173,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "optional_fields/parsing",
            "value": 7149,
            "range": "± 34",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_gc_content",
            "value": 814868,
            "range": "± 21465",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_revcomp",
            "value": 297753,
            "range": "± 3543",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_subsequence",
            "value": 163,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "variant_type",
            "value": 15838,
            "range": "± 28",
            "unit": "ns/iter"
          },
          {
            "name": "quality_error_prob",
            "value": 1215,
            "range": "± 3",
            "unit": "ns/iter"
          },
          {
            "name": "memory_mapped_access",
            "value": 9295,
            "range": "± 168",
            "unit": "ns/iter"
          },
          {
            "name": "hamming_distance",
            "value": 9,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/100",
            "value": 29286,
            "range": "± 181",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/100",
            "value": 53882,
            "range": "± 186",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/100",
            "value": 26035,
            "range": "± 243",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/1000",
            "value": 242299,
            "range": "± 1155",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/1000",
            "value": 406596,
            "range": "± 3609",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/1000",
            "value": 207729,
            "range": "± 1166",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/10000",
            "value": 2453166,
            "range": "± 27697",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/10000",
            "value": 3746541,
            "range": "± 44394",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/10000",
            "value": 2093921,
            "range": "± 15920",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/100",
            "value": 32439,
            "range": "± 209",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/100",
            "value": 41417,
            "range": "± 523",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/100",
            "value": 25037,
            "range": "± 158",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/1000",
            "value": 264191,
            "range": "± 2809",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/1000",
            "value": 274729,
            "range": "± 1403",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/1000",
            "value": 200968,
            "range": "± 1302",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/10000",
            "value": 2588520,
            "range": "± 30844",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/10000",
            "value": 2578404,
            "range": "± 18053",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/10000",
            "value": 1963681,
            "range": "± 17831",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/gc_content",
            "value": 63,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/reverse_complement",
            "value": 61,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/subsequence",
            "value": 23,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/count_pattern",
            "value": 226,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/average_quality",
            "value": 62,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/min_quality",
            "value": 45,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/max_quality",
            "value": 49,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/median_quality",
            "value": 67,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/trim_quality",
            "value": 52,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/reverse_complement",
            "value": 93,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/is_valid",
            "value": 92,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/build_index",
            "value": 348301,
            "range": "± 1404",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/random_access",
            "value": 4931,
            "range": "± 20",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/save_index",
            "value": 221601,
            "range": "± 10950",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/load_index",
            "value": 314148,
            "range": "± 8017",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/filter_by_quality",
            "value": 327577,
            "range": "± 1379",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/filter_by_length",
            "value": 182546,
            "range": "± 1054",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/sample_reads",
            "value": 182583,
            "range": "± 2090",
            "unit": "ns/iter"
          },
          {
            "name": "io_writing/write_fasta",
            "value": 277221,
            "range": "± 12061",
            "unit": "ns/iter"
          },
          {
            "name": "io_writing/write_fastq",
            "value": 8896149,
            "range": "± 44456",
            "unit": "ns/iter"
          },
          {
            "name": "memory_usage/parse_memory_usage",
            "value": 2399277,
            "range": "± 160247",
            "unit": "ns/iter"
          },
          {
            "name": "macro_sum_bench",
            "value": 1,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/100",
            "value": 90819,
            "range": "± 1295",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/1000",
            "value": 868527,
            "range": "± 12174",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/10000",
            "value": 9185318,
            "range": "± 87199",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_iterator/iterator",
            "value": 625497,
            "range": "± 9230",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_writing/writing",
            "value": 3855281,
            "range": "± 32003",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/0",
            "value": 105,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/1",
            "value": 97,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/2",
            "value": 93,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/3",
            "value": 129,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/4",
            "value": 110,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/5",
            "value": 59,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/6",
            "value": 114,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/7",
            "value": 476,
            "range": "± 4",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/0",
            "value": 94,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/1",
            "value": 94,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/2",
            "value": 66,
            "range": "± 3",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/3",
            "value": 94,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/4",
            "value": 94,
            "range": "± 4",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/5",
            "value": 95,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/6",
            "value": 94,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/7",
            "value": 94,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/variant_type_detection",
            "value": 2474,
            "range": "± 6",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/allele_frequency_access",
            "value": 4241,
            "range": "± 41",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/sample_data_access",
            "value": 1523,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "band_pass_normalized_kernel_gatk_defaults",
            "value": 992,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "band_pass_add_then_pop_ready_256_loci",
            "value": 330532,
            "range": "± 1018",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_medium_depth_k10",
            "value": 253502,
            "range": "± 1418",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_high_depth_k10",
            "value": 1822128,
            "range": "± 3124",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/4",
            "value": 24927144,
            "range": "± 467312",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/16",
            "value": 99384513,
            "range": "± 258987",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/64",
            "value": 397732611,
            "range": "± 426397",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/8",
            "value": 548239,
            "range": "± 1369",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/8",
            "value": 299654,
            "range": "± 1126",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/8",
            "value": 1036671,
            "range": "± 15665",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/8",
            "value": 9016230,
            "range": "± 30774",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/32",
            "value": 2193204,
            "range": "± 6194",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/32",
            "value": 881152,
            "range": "± 1822",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/32",
            "value": 4128030,
            "range": "± 3196",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/32",
            "value": 36057819,
            "range": "± 131467",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/64",
            "value": 4386534,
            "range": "± 60142",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/64",
            "value": 1395207,
            "range": "± 7529",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/64",
            "value": 8251892,
            "range": "± 10491",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/64",
            "value": 72109488,
            "range": "± 99372",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/128",
            "value": 8759811,
            "range": "± 14155",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/128",
            "value": 1982864,
            "range": "± 3243",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/128",
            "value": 16499093,
            "range": "± 231630",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/128",
            "value": 144228174,
            "range": "± 219147",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/8",
            "value": 2138709,
            "range": "± 2645",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/8",
            "value": 1207987,
            "range": "± 5950",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/8",
            "value": 4095440,
            "range": "± 6581",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/8",
            "value": 36435188,
            "range": "± 41019",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/32",
            "value": 8543532,
            "range": "± 8647",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/32",
            "value": 3795488,
            "range": "± 3452",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/32",
            "value": 16387886,
            "range": "± 14988",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/32",
            "value": 145751789,
            "range": "± 573720",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/64",
            "value": 17081260,
            "range": "± 24878",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/64",
            "value": 6711021,
            "range": "± 6678",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/64",
            "value": 32764079,
            "range": "± 419019",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/64",
            "value": 291546406,
            "range": "± 1077594",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/128",
            "value": 34166099,
            "range": "± 555471",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/128",
            "value": 10770564,
            "range": "± 14201",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/128",
            "value": 65517047,
            "range": "± 55741",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/128",
            "value": 583559726,
            "range": "± 3199296",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/8",
            "value": 4732328,
            "range": "± 11719",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/8",
            "value": 2709835,
            "range": "± 34164",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/8",
            "value": 9115566,
            "range": "± 11504",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/8",
            "value": 82266954,
            "range": "± 768331",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/32",
            "value": 18926119,
            "range": "± 24190",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/32",
            "value": 8725043,
            "range": "± 12278",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/32",
            "value": 36378231,
            "range": "± 50782",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/32",
            "value": 329020036,
            "range": "± 679421",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/64",
            "value": 37835353,
            "range": "± 42607",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/64",
            "value": 15950231,
            "range": "± 17321",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/64",
            "value": 72727753,
            "range": "± 237429",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/64",
            "value": 658251711,
            "range": "± 1912989",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/128",
            "value": 76112528,
            "range": "± 148423",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/128",
            "value": 27715067,
            "range": "± 35298",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/128",
            "value": 145995733,
            "range": "± 115180",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/128",
            "value": 1316631668,
            "range": "± 2321639",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/value_copy_stages/3072x12",
            "value": 5099,
            "range": "± 20",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/arc_share_stages/3072x12",
            "value": 3923,
            "range": "± 5",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/assembly_result_set_arc_fanout/3072",
            "value": 56,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/assembly_result_set_value_fanout/3072",
            "value": 1209,
            "range": "± 4",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/64x48",
            "value": 16564,
            "range": "± 181",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/64x48",
            "value": 16596,
            "range": "± 191",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/128x96",
            "value": 66150,
            "range": "± 901",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/128x96",
            "value": 66082,
            "range": "± 308",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/256x192",
            "value": 292387,
            "range": "± 970",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/256x192",
            "value": 290026,
            "range": "± 1134",
            "unit": "ns/iter"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "name": "SynapticFour",
            "username": "SynapticFour",
            "email": "contact@synapticfour.com"
          },
          "committer": {
            "name": "GitHub",
            "username": "web-flow",
            "email": "noreply@github.com"
          },
          "id": "41caa6dbe0f5cd78967cd1d9d110c05fd292f9a5",
          "message": "perf(hc): SIMD hapStartIndex PairHMM and cut redundant EventMap/SW wall (#127)\n\nDense losers were scalar in the prefix-reuse PairHMM inner loop and paid\na post-spine EventMap rebuild plus a second SNP pileup that Java does not.\nKeep PairHMM/SW TLS high-water (still well under Java Peak RSS).\n\nCo-authored-by: Cursor <cursoragent@cursor.com>",
          "timestamp": "2026-08-18T15:51:54Z",
          "url": "https://github.com/SynapticFour/gatk-rs/commit/41caa6dbe0f5cd78967cd1d9d110c05fd292f9a5"
        },
        "date": 1787109371980,
        "tool": "cargo",
        "benches": [
          {
            "name": "sam_parsing/parsing/100",
            "value": 85966,
            "range": "± 1664",
            "unit": "ns/iter"
          },
          {
            "name": "sam_parsing/parsing/1000",
            "value": 892165,
            "range": "± 3408",
            "unit": "ns/iter"
          },
          {
            "name": "sam_parsing/parsing/10000",
            "value": 8518348,
            "range": "± 65098",
            "unit": "ns/iter"
          },
          {
            "name": "sam_iterator/iterator",
            "value": 628068,
            "range": "± 4553",
            "unit": "ns/iter"
          },
          {
            "name": "sam_writing/writing",
            "value": 2686388,
            "range": "± 2654966",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/0",
            "value": 99,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/1",
            "value": 137,
            "range": "± 4",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/2",
            "value": 109,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/3",
            "value": 115,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/4",
            "value": 138,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "optional_fields/parsing",
            "value": 5538,
            "range": "± 12",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_gc_content",
            "value": 637232,
            "range": "± 1234",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_revcomp",
            "value": 230546,
            "range": "± 1475",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_subsequence",
            "value": 126,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "variant_type",
            "value": 12833,
            "range": "± 26",
            "unit": "ns/iter"
          },
          {
            "name": "quality_error_prob",
            "value": 950,
            "range": "± 21",
            "unit": "ns/iter"
          },
          {
            "name": "memory_mapped_access",
            "value": 7572,
            "range": "± 163",
            "unit": "ns/iter"
          },
          {
            "name": "hamming_distance",
            "value": 7,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/100",
            "value": 22552,
            "range": "± 157",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/100",
            "value": 40588,
            "range": "± 149",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/100",
            "value": 19446,
            "range": "± 88",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/1000",
            "value": 187732,
            "range": "± 1059",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/1000",
            "value": 307056,
            "range": "± 3064",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/1000",
            "value": 159267,
            "range": "± 569",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/10000",
            "value": 1872187,
            "range": "± 41562",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/10000",
            "value": 2795546,
            "range": "± 28979",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/10000",
            "value": 1561381,
            "range": "± 20504",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/100",
            "value": 25061,
            "range": "± 67",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/100",
            "value": 34806,
            "range": "± 152",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/100",
            "value": 19489,
            "range": "± 248",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/1000",
            "value": 198807,
            "range": "± 816",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/1000",
            "value": 234170,
            "range": "± 767",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/1000",
            "value": 159115,
            "range": "± 1089",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/10000",
            "value": 1949692,
            "range": "± 10337",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/10000",
            "value": 2219349,
            "range": "± 19210",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/10000",
            "value": 1549579,
            "range": "± 72713",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/gc_content",
            "value": 49,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/reverse_complement",
            "value": 48,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/subsequence",
            "value": 18,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/count_pattern",
            "value": 161,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/average_quality",
            "value": 48,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/min_quality",
            "value": 34,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/max_quality",
            "value": 38,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/median_quality",
            "value": 34,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/trim_quality",
            "value": 40,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/reverse_complement",
            "value": 71,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/is_valid",
            "value": 71,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/build_index",
            "value": 271648,
            "range": "± 1142",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/random_access",
            "value": 3846,
            "range": "± 9",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/save_index",
            "value": 524023,
            "range": "± 412222",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/load_index",
            "value": 240161,
            "range": "± 7904",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/filter_by_quality",
            "value": 252017,
            "range": "± 845",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/filter_by_length",
            "value": 140566,
            "range": "± 1008",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/sample_reads",
            "value": 140913,
            "range": "± 1343",
            "unit": "ns/iter"
          },
          {
            "name": "io_writing/write_fasta",
            "value": 510169,
            "range": "± 334813",
            "unit": "ns/iter"
          },
          {
            "name": "io_writing/write_fastq",
            "value": 7161274,
            "range": "± 5428858",
            "unit": "ns/iter"
          },
          {
            "name": "memory_usage/parse_memory_usage",
            "value": 1871868,
            "range": "± 16147",
            "unit": "ns/iter"
          },
          {
            "name": "macro_sum_bench",
            "value": 1,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/100",
            "value": 69410,
            "range": "± 696",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/1000",
            "value": 674372,
            "range": "± 9958",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/10000",
            "value": 7170707,
            "range": "± 79728",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_iterator/iterator",
            "value": 455990,
            "range": "± 2404",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_writing/writing",
            "value": 3122246,
            "range": "± 243841",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/0",
            "value": 83,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/1",
            "value": 75,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/2",
            "value": 73,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/3",
            "value": 101,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/4",
            "value": 87,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/5",
            "value": 47,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/6",
            "value": 88,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/7",
            "value": 372,
            "range": "± 3",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/0",
            "value": 74,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/1",
            "value": 74,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/2",
            "value": 36,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/3",
            "value": 74,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/4",
            "value": 74,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/5",
            "value": 74,
            "range": "± 5",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/6",
            "value": 74,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/7",
            "value": 75,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/variant_type_detection",
            "value": 1645,
            "range": "± 6",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/allele_frequency_access",
            "value": 3283,
            "range": "± 53",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/sample_data_access",
            "value": 1180,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "band_pass_normalized_kernel_gatk_defaults",
            "value": 771,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "band_pass_add_then_pop_ready_256_loci",
            "value": 250368,
            "range": "± 833",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_medium_depth_k10",
            "value": 199054,
            "range": "± 7359",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_high_depth_k10",
            "value": 1427797,
            "range": "± 3270",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/4",
            "value": 19335706,
            "range": "± 223904",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/16",
            "value": 77223273,
            "range": "± 178255",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/64",
            "value": 308980731,
            "range": "± 2219519",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/8",
            "value": 333072,
            "range": "± 2084",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/8",
            "value": 110577,
            "range": "± 184",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/8",
            "value": 741344,
            "range": "± 1559",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/8",
            "value": 7010339,
            "range": "± 17777",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/32",
            "value": 1329875,
            "range": "± 34079",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/32",
            "value": 330537,
            "range": "± 438",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/32",
            "value": 2953550,
            "range": "± 3325",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/32",
            "value": 28035556,
            "range": "± 495664",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/64",
            "value": 2655880,
            "range": "± 7440",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/64",
            "value": 522672,
            "range": "± 3121",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/64",
            "value": 5906767,
            "range": "± 133190",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/64",
            "value": 56066344,
            "range": "± 150101",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/128",
            "value": 5319941,
            "range": "± 12896",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/128",
            "value": 751210,
            "range": "± 13116",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/128",
            "value": 11806787,
            "range": "± 10676",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/128",
            "value": 112138345,
            "range": "± 167956",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/8",
            "value": 1328428,
            "range": "± 10832",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/8",
            "value": 450476,
            "range": "± 297",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/8",
            "value": 2974650,
            "range": "± 2787",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/8",
            "value": 28348690,
            "range": "± 24262",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/32",
            "value": 5315573,
            "range": "± 14303",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/32",
            "value": 1433614,
            "range": "± 6230",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/32",
            "value": 11880301,
            "range": "± 21497",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/32",
            "value": 113315217,
            "range": "± 153134",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/64",
            "value": 10615654,
            "range": "± 37926",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/64",
            "value": 2543028,
            "range": "± 1844",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/64",
            "value": 23762035,
            "range": "± 39232",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/64",
            "value": 226756858,
            "range": "± 188769",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/128",
            "value": 21220445,
            "range": "± 120658",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/128",
            "value": 4068959,
            "range": "± 3619",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/128",
            "value": 47493337,
            "range": "± 57749",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/128",
            "value": 453496344,
            "range": "± 1754681",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/8",
            "value": 2984749,
            "range": "± 7312",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/8",
            "value": 1037231,
            "range": "± 1298",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/8",
            "value": 6705569,
            "range": "± 6600",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/8",
            "value": 63994561,
            "range": "± 225560",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/32",
            "value": 11921991,
            "range": "± 14925",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/32",
            "value": 3390580,
            "range": "± 5141",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/32",
            "value": 26802879,
            "range": "± 65504",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/32",
            "value": 255990868,
            "range": "± 2391211",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/64",
            "value": 23859889,
            "range": "± 172685",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/64",
            "value": 6207118,
            "range": "± 4580",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/64",
            "value": 53603784,
            "range": "± 78432",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/64",
            "value": 511870271,
            "range": "± 3842196",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/128",
            "value": 47715061,
            "range": "± 75182",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/128",
            "value": 10777155,
            "range": "± 226292",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/128",
            "value": 107248428,
            "range": "± 225375",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/128",
            "value": 1025665624,
            "range": "± 3884804",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/value_copy_stages/3072x12",
            "value": 3947,
            "range": "± 18",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/arc_share_stages/3072x12",
            "value": 3041,
            "range": "± 58",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/assembly_result_set_arc_fanout/3072",
            "value": 43,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/assembly_result_set_value_fanout/3072",
            "value": 943,
            "range": "± 25",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/64x48",
            "value": 8440,
            "range": "± 108",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/64x48",
            "value": 8268,
            "range": "± 94",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/128x96",
            "value": 32727,
            "range": "± 124",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/128x96",
            "value": 32672,
            "range": "± 158",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/256x192",
            "value": 123784,
            "range": "± 563",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/256x192",
            "value": 123622,
            "range": "± 602",
            "unit": "ns/iter"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "name": "SynapticFour",
            "username": "SynapticFour",
            "email": "contact@synapticfour.com"
          },
          "committer": {
            "name": "GitHub",
            "username": "web-flow",
            "email": "noreply@github.com"
          },
          "id": "7611aec4df1ae21e300af62dc5e2ebc801421465",
          "message": "perf(hc): wall ledger tip — AD memo, reshape cache, RT/SW/PairHMM TLS (#128)\n\n* perf(hc): wall ledger tip — AD memo, reshape cache, RT keys, SW/PairHMM TLS\n\nCut genotype multi-pass AD rescans and per-allele likelihood reshape; pack\nk-mer keys; drop dead PairHMM prior planes and TLS-reuse transitions; add\nproduction profiler + performance ledger for wall-losers rematch.\n\nLocal mega rematch (21:9825–9828k): assign Σ ~1.2s (was ~130–172s TRACE class).\n\nbaseline-bump: lower clone ratchet 439→414 after tip reduced ownership churn\nCo-authored-by: Cursor <cursoragent@cursor.com>\n\n* docs(perf): drop gitignored runs/ markdown links for CI doc-link gate\n\nCo-authored-by: Cursor <cursoragent@cursor.com>\n\n* docs(perf): record wall-losers 1.15× rematch and PairHMM A/B revert\n\nSign the tip rematch (median 1.15× / Σ 1.27×), loser-head top-3 profiles, and\ndocument prefix-vs-pack knobs as REVERT (occupancy saturated). Name PairHMM\nmin-haps constants without behavior change; gate striped SW behind oracle.\n\nCo-authored-by: Cursor <cursoragent@cursor.com>\n\n* chore(pairhmm): cfg-gate prefix min-haps consts per SIMD arch\n\nSilence unused-const warnings on aarch64 (AVX2) and x86 (NEON).\n\nCo-authored-by: Cursor <cursoragent@cursor.com>\n\n---------\n\nCo-authored-by: Cursor <cursoragent@cursor.com>",
          "timestamp": "2026-08-19T15:47:07Z",
          "url": "https://github.com/SynapticFour/gatk-rs/commit/7611aec4df1ae21e300af62dc5e2ebc801421465"
        },
        "date": 1787196335083,
        "tool": "cargo",
        "benches": [
          {
            "name": "sam_parsing/parsing/100",
            "value": 112267,
            "range": "± 354",
            "unit": "ns/iter"
          },
          {
            "name": "sam_parsing/parsing/1000",
            "value": 1125905,
            "range": "± 8159",
            "unit": "ns/iter"
          },
          {
            "name": "sam_parsing/parsing/10000",
            "value": 11073265,
            "range": "± 36514",
            "unit": "ns/iter"
          },
          {
            "name": "sam_iterator/iterator",
            "value": 867700,
            "range": "± 4989",
            "unit": "ns/iter"
          },
          {
            "name": "sam_writing/writing",
            "value": 3331024,
            "range": "± 27430",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/0",
            "value": 138,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/1",
            "value": 180,
            "range": "± 4",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/2",
            "value": 144,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/3",
            "value": 153,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/4",
            "value": 177,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "optional_fields/parsing",
            "value": 7171,
            "range": "± 69",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_gc_content",
            "value": 814786,
            "range": "± 3661",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_revcomp",
            "value": 297291,
            "range": "± 408",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_subsequence",
            "value": 164,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "variant_type",
            "value": 15865,
            "range": "± 29",
            "unit": "ns/iter"
          },
          {
            "name": "quality_error_prob",
            "value": 1215,
            "range": "± 3",
            "unit": "ns/iter"
          },
          {
            "name": "memory_mapped_access",
            "value": 9281,
            "range": "± 213",
            "unit": "ns/iter"
          },
          {
            "name": "hamming_distance",
            "value": 9,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/100",
            "value": 29281,
            "range": "± 308",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/100",
            "value": 51160,
            "range": "± 193",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/100",
            "value": 25950,
            "range": "± 319",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/1000",
            "value": 242843,
            "range": "± 445",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/1000",
            "value": 389082,
            "range": "± 2534",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/1000",
            "value": 204053,
            "range": "± 2333",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/10000",
            "value": 2429852,
            "range": "± 21126",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/10000",
            "value": 3613090,
            "range": "± 27000",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/10000",
            "value": 2036283,
            "range": "± 28608",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/100",
            "value": 33519,
            "range": "± 109",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/100",
            "value": 42983,
            "range": "± 160",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/100",
            "value": 25433,
            "range": "± 151",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/1000",
            "value": 267792,
            "range": "± 1650",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/1000",
            "value": 291663,
            "range": "± 4411",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/1000",
            "value": 205813,
            "range": "± 1713",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/10000",
            "value": 2687322,
            "range": "± 27121",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/10000",
            "value": 2662698,
            "range": "± 25469",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/10000",
            "value": 2007593,
            "range": "± 35224",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/gc_content",
            "value": 63,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/reverse_complement",
            "value": 60,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/subsequence",
            "value": 23,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/count_pattern",
            "value": 209,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/average_quality",
            "value": 62,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/min_quality",
            "value": 45,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/max_quality",
            "value": 49,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/median_quality",
            "value": 45,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/trim_quality",
            "value": 54,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/reverse_complement",
            "value": 91,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/is_valid",
            "value": 92,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/build_index",
            "value": 344302,
            "range": "± 6011",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/random_access",
            "value": 4978,
            "range": "± 11",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/save_index",
            "value": 219579,
            "range": "± 23399",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/load_index",
            "value": 303998,
            "range": "± 19993",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/filter_by_quality",
            "value": 333527,
            "range": "± 1013",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/filter_by_length",
            "value": 186903,
            "range": "± 650",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/sample_reads",
            "value": 190125,
            "range": "± 1224",
            "unit": "ns/iter"
          },
          {
            "name": "io_writing/write_fasta",
            "value": 270297,
            "range": "± 12322",
            "unit": "ns/iter"
          },
          {
            "name": "io_writing/write_fastq",
            "value": 9045193,
            "range": "± 380835",
            "unit": "ns/iter"
          },
          {
            "name": "memory_usage/parse_memory_usage",
            "value": 2438010,
            "range": "± 16832",
            "unit": "ns/iter"
          },
          {
            "name": "macro_sum_bench",
            "value": 1,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/100",
            "value": 91871,
            "range": "± 579",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/1000",
            "value": 890098,
            "range": "± 4992",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/10000",
            "value": 9539930,
            "range": "± 100790",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_iterator/iterator",
            "value": 621703,
            "range": "± 7868",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_writing/writing",
            "value": 3890986,
            "range": "± 20514",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/0",
            "value": 106,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/1",
            "value": 96,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/2",
            "value": 93,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/3",
            "value": 130,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/4",
            "value": 111,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/5",
            "value": 60,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/6",
            "value": 113,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/7",
            "value": 489,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/0",
            "value": 91,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/1",
            "value": 91,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/2",
            "value": 47,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/3",
            "value": 92,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/4",
            "value": 91,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/5",
            "value": 91,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/6",
            "value": 91,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/7",
            "value": 91,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/variant_type_detection",
            "value": 2475,
            "range": "± 13",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/allele_frequency_access",
            "value": 4243,
            "range": "± 10",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/sample_data_access",
            "value": 1524,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "band_pass_normalized_kernel_gatk_defaults",
            "value": 1000,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "band_pass_add_then_pop_ready_256_loci",
            "value": 323933,
            "range": "± 1135",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_medium_depth_k10",
            "value": 1619087,
            "range": "± 3524",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_high_depth_k10",
            "value": 9777019,
            "range": "± 41166",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_depth/full_pipeline_low_k10/16",
            "value": 635619,
            "range": "± 1004",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_depth/threading_build_low_k10/16",
            "value": 432858,
            "range": "± 900",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_depth/full_pipeline_low_k25/16",
            "value": 657475,
            "range": "± 1416",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_depth/threading_build_low_k25/16",
            "value": 449879,
            "range": "± 710",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_depth/full_pipeline_medium_k10/64",
            "value": 1620580,
            "range": "± 4317",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_depth/threading_build_medium_k10/64",
            "value": 1298332,
            "range": "± 1480",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_depth/full_pipeline_medium_k25/64",
            "value": 1782609,
            "range": "± 9217",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_depth/threading_build_medium_k25/64",
            "value": 1388710,
            "range": "± 4774",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_depth/full_pipeline_high_k10/512",
            "value": 9754589,
            "range": "± 37984",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_depth/threading_build_high_k10/512",
            "value": 9343641,
            "range": "± 7866",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_depth/full_pipeline_high_k25/512",
            "value": 9408271,
            "range": "± 16122",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_depth/threading_build_high_k25/512",
            "value": 8943625,
            "range": "± 40436",
            "unit": "ns/iter"
          },
          {
            "name": "kmer_key_representation/A_arc_bytes/10",
            "value": 395882,
            "range": "± 498",
            "unit": "ns/iter"
          },
          {
            "name": "kmer_key_representation/B_packed_key/10",
            "value": 356116,
            "range": "± 638",
            "unit": "ns/iter"
          },
          {
            "name": "kmer_key_representation/A_arc_bytes/25",
            "value": 433106,
            "range": "± 409",
            "unit": "ns/iter"
          },
          {
            "name": "kmer_key_representation/B_packed_key/25",
            "value": 361947,
            "range": "± 564",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/4",
            "value": 24866847,
            "range": "± 46236",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/16",
            "value": 99414565,
            "range": "± 149521",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/64",
            "value": 397938293,
            "range": "± 325002",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/8",
            "value": 423445,
            "range": "± 614",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/8",
            "value": 141560,
            "range": "± 310",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/8",
            "value": 881457,
            "range": "± 780",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/8",
            "value": 9020301,
            "range": "± 32182",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/32",
            "value": 1692117,
            "range": "± 3939",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/32",
            "value": 422990,
            "range": "± 797",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/32",
            "value": 3519716,
            "range": "± 2783",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/32",
            "value": 36078564,
            "range": "± 70562",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/64",
            "value": 3384389,
            "range": "± 4574",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/64",
            "value": 667840,
            "range": "± 1591",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/64",
            "value": 7034451,
            "range": "± 13296",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/64",
            "value": 72158700,
            "range": "± 226003",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/128",
            "value": 6765307,
            "range": "± 5569",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/128",
            "value": 959896,
            "range": "± 1136",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/128",
            "value": 14066985,
            "range": "± 101207",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/128",
            "value": 144341700,
            "range": "± 260693",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/8",
            "value": 1712378,
            "range": "± 2772",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/8",
            "value": 578037,
            "range": "± 534",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/8",
            "value": 3606905,
            "range": "± 5912",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/8",
            "value": 36441941,
            "range": "± 105537",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/32",
            "value": 6840654,
            "range": "± 8034",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/32",
            "value": 1845979,
            "range": "± 3977",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/32",
            "value": 14412450,
            "range": "± 11866",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/32",
            "value": 145767109,
            "range": "± 220198",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/64",
            "value": 13672013,
            "range": "± 16384",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/64",
            "value": 3272467,
            "range": "± 3942",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/64",
            "value": 28813971,
            "range": "± 40356",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/64",
            "value": 291649271,
            "range": "± 805829",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/128",
            "value": 27340723,
            "range": "± 55693",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/128",
            "value": 5230008,
            "range": "± 6226",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/128",
            "value": 57610700,
            "range": "± 61630",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/128",
            "value": 583603579,
            "range": "± 583689",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/8",
            "value": 3853198,
            "range": "± 12529",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/8",
            "value": 1310174,
            "range": "± 3651",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/8",
            "value": 8121322,
            "range": "± 9448",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/8",
            "value": 82241473,
            "range": "± 129540",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/32",
            "value": 15399148,
            "range": "± 8384",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/32",
            "value": 4276663,
            "range": "± 4252",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/32",
            "value": 32463634,
            "range": "± 29731",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/32",
            "value": 329017195,
            "range": "± 266463",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/64",
            "value": 30770327,
            "range": "± 33988",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/64",
            "value": 7842389,
            "range": "± 10017",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/64",
            "value": 65020606,
            "range": "± 102584",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/64",
            "value": 657875712,
            "range": "± 389076",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/128",
            "value": 61551730,
            "range": "± 82607",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/128",
            "value": 13629728,
            "range": "± 18912",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/128",
            "value": 130040630,
            "range": "± 106038",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/128",
            "value": 1316095842,
            "range": "± 891206",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/logless_scalar_1x1/1",
            "value": 214270,
            "range": "± 256",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/simd_1x1/1",
            "value": 209303,
            "range": "± 899",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/wavefront_1x1/1",
            "value": 102598,
            "range": "± 154",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/logless_scalar_1x8/8",
            "value": 1711110,
            "range": "± 3765",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/simd_1x8/8",
            "value": 581053,
            "range": "± 614",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/wavefront_1x8/8",
            "value": 790642,
            "range": "± 658",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/logless_scalar_1x16/16",
            "value": 3419839,
            "range": "± 11756",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/simd_1x16/16",
            "value": 1024781,
            "range": "± 3261",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/wavefront_1x16/16",
            "value": 1586108,
            "range": "± 2319",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/logless_scalar_1x32/32",
            "value": 6835862,
            "range": "± 8523",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/simd_1x32/32",
            "value": 1857854,
            "range": "± 3568",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/wavefront_1x32/32",
            "value": 3168835,
            "range": "± 7172",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/logless_scalar_1x64/64",
            "value": 13667556,
            "range": "± 10787",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/simd_1x64/64",
            "value": 3296472,
            "range": "± 4390",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/wavefront_1x64/64",
            "value": 6334876,
            "range": "± 5396",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/logless_scalar_tile_8x32/8",
            "value": 54644031,
            "range": "± 58891",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/simd_tile_8x32/8",
            "value": 14857450,
            "range": "± 14869",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/wavefront_tile_8x32/8",
            "value": 25325016,
            "range": "± 115003",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/logless_scalar_tile_16x32/16",
            "value": 109332543,
            "range": "± 101238",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/simd_tile_16x32/16",
            "value": 29720806,
            "range": "± 38971",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/wavefront_tile_16x32/16",
            "value": 50620536,
            "range": "± 135020",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/logless_scalar_realistic_h/8",
            "value": 1880003,
            "range": "± 4701",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/simd_realistic_h/8",
            "value": 154820,
            "range": "± 243",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/wavefront_realistic_h/8",
            "value": 887635,
            "range": "± 782",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/logless_scalar_realistic_h/32",
            "value": 7514312,
            "range": "± 9017",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/simd_realistic_h/32",
            "value": 186368,
            "range": "± 430",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/wavefront_realistic_h/32",
            "value": 3517807,
            "range": "± 2276",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/logless_scalar_realistic_h/64",
            "value": 15032492,
            "range": "± 14248",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/simd_realistic_h/64",
            "value": 228642,
            "range": "± 675",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/wavefront_realistic_h/64",
            "value": 7032959,
            "range": "± 17847",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/value_copy_stages/3072x12",
            "value": 5104,
            "range": "± 28",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/arc_share_stages/3072x12",
            "value": 3921,
            "range": "± 5",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/assembly_result_set_arc_fanout/3072",
            "value": 56,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/assembly_result_set_value_fanout/3072",
            "value": 1182,
            "range": "± 8",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/64x48",
            "value": 10359,
            "range": "± 70",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/64x48",
            "value": 10292,
            "range": "± 24",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/128x96",
            "value": 43039,
            "range": "± 119",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/128x96",
            "value": 43125,
            "range": "± 93",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/256x192",
            "value": 173919,
            "range": "± 1118",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/256x192",
            "value": 173733,
            "range": "± 660",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/200x151",
            "value": 108556,
            "range": "± 207",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/200x151",
            "value": 108808,
            "range": "± 400",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/100x100",
            "value": 35536,
            "range": "± 361",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/100x100",
            "value": 35704,
            "range": "± 92",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_hc/padded_indel/80",
            "value": 36213,
            "range": "± 110",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_hc/padded_indel/150",
            "value": 108238,
            "range": "± 860",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_hc/padded_indel/250",
            "value": 255927,
            "range": "± 919",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_hc/read_to_hap_soft/hap120_read100",
            "value": 40836,
            "range": "± 116",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_hc/read_to_hap_soft/hap200_read151",
            "value": 102783,
            "range": "± 677",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_hc/read_to_hap_soft/hap280_read151",
            "value": 141825,
            "range": "± 380",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_hc/exact_substring_fast_path",
            "value": 224,
            "range": "± 1",
            "unit": "ns/iter"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "name": "SynapticFour",
            "username": "SynapticFour",
            "email": "contact@synapticfour.com"
          },
          "committer": {
            "name": "GitHub",
            "username": "web-flow",
            "email": "noreply@github.com"
          },
          "id": "7611aec4df1ae21e300af62dc5e2ebc801421465",
          "message": "perf(hc): wall ledger tip — AD memo, reshape cache, RT/SW/PairHMM TLS (#128)\n\n* perf(hc): wall ledger tip — AD memo, reshape cache, RT keys, SW/PairHMM TLS\n\nCut genotype multi-pass AD rescans and per-allele likelihood reshape; pack\nk-mer keys; drop dead PairHMM prior planes and TLS-reuse transitions; add\nproduction profiler + performance ledger for wall-losers rematch.\n\nLocal mega rematch (21:9825–9828k): assign Σ ~1.2s (was ~130–172s TRACE class).\n\nbaseline-bump: lower clone ratchet 439→414 after tip reduced ownership churn\nCo-authored-by: Cursor <cursoragent@cursor.com>\n\n* docs(perf): drop gitignored runs/ markdown links for CI doc-link gate\n\nCo-authored-by: Cursor <cursoragent@cursor.com>\n\n* docs(perf): record wall-losers 1.15× rematch and PairHMM A/B revert\n\nSign the tip rematch (median 1.15× / Σ 1.27×), loser-head top-3 profiles, and\ndocument prefix-vs-pack knobs as REVERT (occupancy saturated). Name PairHMM\nmin-haps constants without behavior change; gate striped SW behind oracle.\n\nCo-authored-by: Cursor <cursoragent@cursor.com>\n\n* chore(pairhmm): cfg-gate prefix min-haps consts per SIMD arch\n\nSilence unused-const warnings on aarch64 (AVX2) and x86 (NEON).\n\nCo-authored-by: Cursor <cursoragent@cursor.com>\n\n---------\n\nCo-authored-by: Cursor <cursoragent@cursor.com>",
          "timestamp": "2026-08-19T15:47:07Z",
          "url": "https://github.com/SynapticFour/gatk-rs/commit/7611aec4df1ae21e300af62dc5e2ebc801421465"
        },
        "date": 1787283198208,
        "tool": "cargo",
        "benches": [
          {
            "name": "sam_parsing/parsing/100",
            "value": 109115,
            "range": "± 280",
            "unit": "ns/iter"
          },
          {
            "name": "sam_parsing/parsing/1000",
            "value": 1132746,
            "range": "± 2798",
            "unit": "ns/iter"
          },
          {
            "name": "sam_parsing/parsing/10000",
            "value": 10728891,
            "range": "± 29427",
            "unit": "ns/iter"
          },
          {
            "name": "sam_iterator/iterator",
            "value": 789414,
            "range": "± 8935",
            "unit": "ns/iter"
          },
          {
            "name": "sam_writing/writing",
            "value": 3290390,
            "range": "± 21505",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/0",
            "value": 130,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/1",
            "value": 172,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/2",
            "value": 141,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/3",
            "value": 148,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/4",
            "value": 174,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "optional_fields/parsing",
            "value": 7140,
            "range": "± 68",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_gc_content",
            "value": 814724,
            "range": "± 2967",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_revcomp",
            "value": 297627,
            "range": "± 1029",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_subsequence",
            "value": 161,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "variant_type",
            "value": 16223,
            "range": "± 158",
            "unit": "ns/iter"
          },
          {
            "name": "quality_error_prob",
            "value": 1215,
            "range": "± 3",
            "unit": "ns/iter"
          },
          {
            "name": "memory_mapped_access",
            "value": 9682,
            "range": "± 234",
            "unit": "ns/iter"
          },
          {
            "name": "hamming_distance",
            "value": 9,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/100",
            "value": 29935,
            "range": "± 132",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/100",
            "value": 50405,
            "range": "± 173",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/100",
            "value": 27063,
            "range": "± 315",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/1000",
            "value": 243732,
            "range": "± 912",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/1000",
            "value": 383836,
            "range": "± 2601",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/1000",
            "value": 215966,
            "range": "± 3705",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/10000",
            "value": 2453377,
            "range": "± 31488",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/10000",
            "value": 3565913,
            "range": "± 15233",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/10000",
            "value": 2142980,
            "range": "± 31687",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/100",
            "value": 32280,
            "range": "± 231",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/100",
            "value": 43248,
            "range": "± 429",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/100",
            "value": 25243,
            "range": "± 213",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/1000",
            "value": 263649,
            "range": "± 752",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/1000",
            "value": 285477,
            "range": "± 955",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/1000",
            "value": 199452,
            "range": "± 1527",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/10000",
            "value": 2600021,
            "range": "± 18039",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/10000",
            "value": 2664445,
            "range": "± 18193",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/10000",
            "value": 1983101,
            "range": "± 18284",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/gc_content",
            "value": 63,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/reverse_complement",
            "value": 60,
            "range": "± 4",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/subsequence",
            "value": 23,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/count_pattern",
            "value": 208,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/average_quality",
            "value": 62,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/min_quality",
            "value": 45,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/max_quality",
            "value": 49,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/median_quality",
            "value": 44,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/trim_quality",
            "value": 54,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/reverse_complement",
            "value": 91,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/is_valid",
            "value": 92,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/build_index",
            "value": 347124,
            "range": "± 3370",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/random_access",
            "value": 4939,
            "range": "± 18",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/save_index",
            "value": 224737,
            "range": "± 14938",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/load_index",
            "value": 310043,
            "range": "± 2327",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/filter_by_quality",
            "value": 330907,
            "range": "± 12708",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/filter_by_length",
            "value": 186366,
            "range": "± 918",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/sample_reads",
            "value": 189420,
            "range": "± 651",
            "unit": "ns/iter"
          },
          {
            "name": "io_writing/write_fasta",
            "value": 210430,
            "range": "± 12659",
            "unit": "ns/iter"
          },
          {
            "name": "io_writing/write_fastq",
            "value": 9129248,
            "range": "± 108947",
            "unit": "ns/iter"
          },
          {
            "name": "memory_usage/parse_memory_usage",
            "value": 2444663,
            "range": "± 15393",
            "unit": "ns/iter"
          },
          {
            "name": "macro_sum_bench",
            "value": 1,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/100",
            "value": 89938,
            "range": "± 1565",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/1000",
            "value": 857790,
            "range": "± 2466",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/10000",
            "value": 9208798,
            "range": "± 89412",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_iterator/iterator",
            "value": 595222,
            "range": "± 4435",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_writing/writing",
            "value": 3834084,
            "range": "± 19224",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/0",
            "value": 107,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/1",
            "value": 99,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/2",
            "value": 97,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/3",
            "value": 131,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/4",
            "value": 111,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/5",
            "value": 63,
            "range": "± 4",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/6",
            "value": 116,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/7",
            "value": 455,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/0",
            "value": 96,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/1",
            "value": 95,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/2",
            "value": 50,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/3",
            "value": 95,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/4",
            "value": 95,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/5",
            "value": 95,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/6",
            "value": 95,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/7",
            "value": 95,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/variant_type_detection",
            "value": 2472,
            "range": "± 6",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/allele_frequency_access",
            "value": 4234,
            "range": "± 9",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/sample_data_access",
            "value": 1523,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "band_pass_normalized_kernel_gatk_defaults",
            "value": 994,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "band_pass_add_then_pop_ready_256_loci",
            "value": 318980,
            "range": "± 1775",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_medium_depth_k10",
            "value": 1601763,
            "range": "± 4398",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_high_depth_k10",
            "value": 9572451,
            "range": "± 107425",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_depth/full_pipeline_low_k10/16",
            "value": 625971,
            "range": "± 627",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_depth/threading_build_low_k10/16",
            "value": 426492,
            "range": "± 609",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_depth/full_pipeline_low_k25/16",
            "value": 651412,
            "range": "± 1298",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_depth/threading_build_low_k25/16",
            "value": 444791,
            "range": "± 2101",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_depth/full_pipeline_medium_k10/64",
            "value": 1597424,
            "range": "± 2105",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_depth/threading_build_medium_k10/64",
            "value": 1279803,
            "range": "± 2010",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_depth/full_pipeline_medium_k25/64",
            "value": 1758301,
            "range": "± 2308",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_depth/threading_build_medium_k25/64",
            "value": 1377237,
            "range": "± 2148",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_depth/full_pipeline_high_k10/512",
            "value": 9591624,
            "range": "± 44251",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_depth/threading_build_high_k10/512",
            "value": 9189127,
            "range": "± 34063",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_depth/full_pipeline_high_k25/512",
            "value": 9232240,
            "range": "± 16501",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_depth/threading_build_high_k25/512",
            "value": 8770429,
            "range": "± 60386",
            "unit": "ns/iter"
          },
          {
            "name": "kmer_key_representation/A_arc_bytes/10",
            "value": 396859,
            "range": "± 1555",
            "unit": "ns/iter"
          },
          {
            "name": "kmer_key_representation/B_packed_key/10",
            "value": 410556,
            "range": "± 956",
            "unit": "ns/iter"
          },
          {
            "name": "kmer_key_representation/A_arc_bytes/25",
            "value": 434089,
            "range": "± 736",
            "unit": "ns/iter"
          },
          {
            "name": "kmer_key_representation/B_packed_key/25",
            "value": 413097,
            "range": "± 1063",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/4",
            "value": 24841656,
            "range": "± 27698",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/16",
            "value": 99436002,
            "range": "± 123073",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/64",
            "value": 397854694,
            "range": "± 530810",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/8",
            "value": 423111,
            "range": "± 481",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/8",
            "value": 141493,
            "range": "± 485",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/8",
            "value": 942716,
            "range": "± 1349",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/8",
            "value": 9011431,
            "range": "± 18242",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/32",
            "value": 1691190,
            "range": "± 2840",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/32",
            "value": 423396,
            "range": "± 572",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/32",
            "value": 3766354,
            "range": "± 3059",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/32",
            "value": 36039454,
            "range": "± 89537",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/64",
            "value": 3380881,
            "range": "± 2868",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/64",
            "value": 668022,
            "range": "± 834",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/64",
            "value": 7538021,
            "range": "± 15135",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/64",
            "value": 72091966,
            "range": "± 83026",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/128",
            "value": 6761244,
            "range": "± 8829",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/128",
            "value": 961589,
            "range": "± 766",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/128",
            "value": 15079801,
            "range": "± 17545",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/128",
            "value": 144131132,
            "range": "± 175601",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/8",
            "value": 1708353,
            "range": "± 7422",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/8",
            "value": 577070,
            "range": "± 1019",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/8",
            "value": 3809076,
            "range": "± 4350",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/8",
            "value": 36403302,
            "range": "± 88274",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/32",
            "value": 6824451,
            "range": "± 9393",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/32",
            "value": 1843155,
            "range": "± 4005",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/32",
            "value": 15225211,
            "range": "± 22493",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/32",
            "value": 145597865,
            "range": "± 422937",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/64",
            "value": 13645602,
            "range": "± 18374",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/64",
            "value": 3265778,
            "range": "± 4258",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/64",
            "value": 30447496,
            "range": "± 27057",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/64",
            "value": 291362040,
            "range": "± 915822",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/128",
            "value": 27280556,
            "range": "± 35141",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/128",
            "value": 5219088,
            "range": "± 9070",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/128",
            "value": 60889247,
            "range": "± 145738",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/128",
            "value": 582704442,
            "range": "± 335917",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/8",
            "value": 3845450,
            "range": "± 2962",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/8",
            "value": 1313550,
            "range": "± 1800",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/8",
            "value": 8590274,
            "range": "± 11679",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/8",
            "value": 82226746,
            "range": "± 111048",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/32",
            "value": 15378482,
            "range": "± 17313",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/32",
            "value": 4288210,
            "range": "± 7090",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/32",
            "value": 34375478,
            "range": "± 37327",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/32",
            "value": 328872168,
            "range": "± 391589",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/64",
            "value": 30764035,
            "range": "± 58235",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/64",
            "value": 7849638,
            "range": "± 6306",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/64",
            "value": 68710540,
            "range": "± 63233",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/64",
            "value": 657626099,
            "range": "± 962980",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/128",
            "value": 61515515,
            "range": "± 48685",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/128",
            "value": 13644444,
            "range": "± 20194",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/128",
            "value": 137371146,
            "range": "± 123354",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/128",
            "value": 1315773211,
            "range": "± 714716",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/logless_scalar_1x1/1",
            "value": 213638,
            "range": "± 480",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/simd_1x1/1",
            "value": 210904,
            "range": "± 226",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/wavefront_1x1/1",
            "value": 102474,
            "range": "± 121",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/logless_scalar_1x8/8",
            "value": 1706031,
            "range": "± 6089",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/simd_1x8/8",
            "value": 581112,
            "range": "± 663",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/wavefront_1x8/8",
            "value": 790062,
            "range": "± 823",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/logless_scalar_1x16/16",
            "value": 3411750,
            "range": "± 3309",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/simd_1x16/16",
            "value": 1024810,
            "range": "± 4342",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/wavefront_1x16/16",
            "value": 1583216,
            "range": "± 2082",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/logless_scalar_1x32/32",
            "value": 6820828,
            "range": "± 9193",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/simd_1x32/32",
            "value": 1856840,
            "range": "± 7610",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/wavefront_1x32/32",
            "value": 3164176,
            "range": "± 2173",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/logless_scalar_1x64/64",
            "value": 13636931,
            "range": "± 18195",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/simd_1x64/64",
            "value": 3294588,
            "range": "± 4563",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/wavefront_1x64/64",
            "value": 6296795,
            "range": "± 17347",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/logless_scalar_tile_8x32/8",
            "value": 54576686,
            "range": "± 57734",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/simd_tile_8x32/8",
            "value": 14850642,
            "range": "± 14145",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/wavefront_tile_8x32/8",
            "value": 25184204,
            "range": "± 14139",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/logless_scalar_tile_16x32/16",
            "value": 109160311,
            "range": "± 352780",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/simd_tile_16x32/16",
            "value": 29697916,
            "range": "± 23492",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/wavefront_tile_16x32/16",
            "value": 50355839,
            "range": "± 45254",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/logless_scalar_realistic_h/8",
            "value": 1880993,
            "range": "± 2997",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/simd_realistic_h/8",
            "value": 155043,
            "range": "± 376",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/wavefront_realistic_h/8",
            "value": 887267,
            "range": "± 555",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/logless_scalar_realistic_h/32",
            "value": 7520567,
            "range": "± 7914",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/simd_realistic_h/32",
            "value": 186299,
            "range": "± 474",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/wavefront_realistic_h/32",
            "value": 3515804,
            "range": "± 2938",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/logless_scalar_realistic_h/64",
            "value": 15038232,
            "range": "± 10074",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/simd_realistic_h/64",
            "value": 227809,
            "range": "± 442",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/wavefront_realistic_h/64",
            "value": 7027942,
            "range": "± 7175",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/value_copy_stages/3072x12",
            "value": 5097,
            "range": "± 9",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/arc_share_stages/3072x12",
            "value": 3921,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/assembly_result_set_arc_fanout/3072",
            "value": 56,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/assembly_result_set_value_fanout/3072",
            "value": 1166,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/64x48",
            "value": 10327,
            "range": "± 22",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/64x48",
            "value": 10272,
            "range": "± 21",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/128x96",
            "value": 42864,
            "range": "± 159",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/128x96",
            "value": 43057,
            "range": "± 202",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/256x192",
            "value": 168279,
            "range": "± 241",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/256x192",
            "value": 167881,
            "range": "± 488",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/200x151",
            "value": 104282,
            "range": "± 208",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/200x151",
            "value": 104370,
            "range": "± 439",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/100x100",
            "value": 36489,
            "range": "± 366",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/100x100",
            "value": 35834,
            "range": "± 107",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_hc/padded_indel/80",
            "value": 36324,
            "range": "± 230",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_hc/padded_indel/150",
            "value": 103770,
            "range": "± 475",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_hc/padded_indel/250",
            "value": 255514,
            "range": "± 637",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_hc/read_to_hap_soft/hap120_read100",
            "value": 40899,
            "range": "± 199",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_hc/read_to_hap_soft/hap200_read151",
            "value": 102189,
            "range": "± 663",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_hc/read_to_hap_soft/hap280_read151",
            "value": 140983,
            "range": "± 459",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_hc/exact_substring_fast_path",
            "value": 225,
            "range": "± 1",
            "unit": "ns/iter"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "name": "SynapticFour",
            "username": "SynapticFour",
            "email": "contact@synapticfour.com"
          },
          "committer": {
            "name": "GitHub",
            "username": "web-flow",
            "email": "noreply@github.com"
          },
          "id": "7611aec4df1ae21e300af62dc5e2ebc801421465",
          "message": "perf(hc): wall ledger tip — AD memo, reshape cache, RT/SW/PairHMM TLS (#128)\n\n* perf(hc): wall ledger tip — AD memo, reshape cache, RT keys, SW/PairHMM TLS\n\nCut genotype multi-pass AD rescans and per-allele likelihood reshape; pack\nk-mer keys; drop dead PairHMM prior planes and TLS-reuse transitions; add\nproduction profiler + performance ledger for wall-losers rematch.\n\nLocal mega rematch (21:9825–9828k): assign Σ ~1.2s (was ~130–172s TRACE class).\n\nbaseline-bump: lower clone ratchet 439→414 after tip reduced ownership churn\nCo-authored-by: Cursor <cursoragent@cursor.com>\n\n* docs(perf): drop gitignored runs/ markdown links for CI doc-link gate\n\nCo-authored-by: Cursor <cursoragent@cursor.com>\n\n* docs(perf): record wall-losers 1.15× rematch and PairHMM A/B revert\n\nSign the tip rematch (median 1.15× / Σ 1.27×), loser-head top-3 profiles, and\ndocument prefix-vs-pack knobs as REVERT (occupancy saturated). Name PairHMM\nmin-haps constants without behavior change; gate striped SW behind oracle.\n\nCo-authored-by: Cursor <cursoragent@cursor.com>\n\n* chore(pairhmm): cfg-gate prefix min-haps consts per SIMD arch\n\nSilence unused-const warnings on aarch64 (AVX2) and x86 (NEON).\n\nCo-authored-by: Cursor <cursoragent@cursor.com>\n\n---------\n\nCo-authored-by: Cursor <cursoragent@cursor.com>",
          "timestamp": "2026-08-19T15:47:07Z",
          "url": "https://github.com/SynapticFour/gatk-rs/commit/7611aec4df1ae21e300af62dc5e2ebc801421465"
        },
        "date": 1787369249665,
        "tool": "cargo",
        "benches": [
          {
            "name": "sam_parsing/parsing/100",
            "value": 117721,
            "range": "± 2406",
            "unit": "ns/iter"
          },
          {
            "name": "sam_parsing/parsing/1000",
            "value": 1158604,
            "range": "± 13857",
            "unit": "ns/iter"
          },
          {
            "name": "sam_parsing/parsing/10000",
            "value": 11537683,
            "range": "± 479711",
            "unit": "ns/iter"
          },
          {
            "name": "sam_iterator/iterator",
            "value": 875757,
            "range": "± 14580",
            "unit": "ns/iter"
          },
          {
            "name": "sam_writing/writing",
            "value": 3205490,
            "range": "± 64410",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/0",
            "value": 125,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/1",
            "value": 170,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/2",
            "value": 136,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/3",
            "value": 142,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "cigar_parsing/parse/4",
            "value": 169,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "optional_fields/parsing",
            "value": 6462,
            "range": "± 42",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_gc_content",
            "value": 717712,
            "range": "± 13796",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_revcomp",
            "value": 259345,
            "range": "± 585",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_subsequence",
            "value": 195,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "variant_type",
            "value": 16840,
            "range": "± 81",
            "unit": "ns/iter"
          },
          {
            "name": "quality_error_prob",
            "value": 1279,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "memory_mapped_access",
            "value": 8583,
            "range": "± 149",
            "unit": "ns/iter"
          },
          {
            "name": "hamming_distance",
            "value": 8,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/100",
            "value": 29821,
            "range": "± 312",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/100",
            "value": 52505,
            "range": "± 173",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/100",
            "value": 25540,
            "range": "± 163",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/1000",
            "value": 258143,
            "range": "± 1153",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/1000",
            "value": 410490,
            "range": "± 1567",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/1000",
            "value": 217531,
            "range": "± 1439",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/buffered/10000",
            "value": 2596539,
            "range": "± 23167",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/memory_mapped/10000",
            "value": 3903663,
            "range": "± 21966",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_parsing/iterator/10000",
            "value": 2190380,
            "range": "± 31700",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/100",
            "value": 32003,
            "range": "± 123",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/100",
            "value": 44862,
            "range": "± 238",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/100",
            "value": 23337,
            "range": "± 103",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/1000",
            "value": 262643,
            "range": "± 3528",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/1000",
            "value": 305637,
            "range": "± 825",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/1000",
            "value": 190994,
            "range": "± 1170",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/buffered/10000",
            "value": 2518552,
            "range": "± 16702",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/memory_mapped/10000",
            "value": 2824912,
            "range": "± 13107",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_parsing/iterator/10000",
            "value": 1869455,
            "range": "± 20683",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/gc_content",
            "value": 57,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/reverse_complement",
            "value": 55,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/subsequence",
            "value": 20,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "sequence_operations/count_pattern",
            "value": 205,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/average_quality",
            "value": 54,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/min_quality",
            "value": 40,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/max_quality",
            "value": 42,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/median_quality",
            "value": 48,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/trim_quality",
            "value": 50,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/reverse_complement",
            "value": 89,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_operations/is_valid",
            "value": 85,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/build_index",
            "value": 351512,
            "range": "± 1541",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/random_access",
            "value": 4344,
            "range": "± 9",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/save_index",
            "value": 299686,
            "range": "± 13866",
            "unit": "ns/iter"
          },
          {
            "name": "fasta_indexing/load_index",
            "value": 295436,
            "range": "± 3365",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/filter_by_quality",
            "value": 328929,
            "range": "± 8643",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/filter_by_length",
            "value": 179189,
            "range": "± 2632",
            "unit": "ns/iter"
          },
          {
            "name": "fastq_filtering/sample_reads",
            "value": 180435,
            "range": "± 902",
            "unit": "ns/iter"
          },
          {
            "name": "io_writing/write_fasta",
            "value": 300754,
            "range": "± 57303",
            "unit": "ns/iter"
          },
          {
            "name": "io_writing/write_fastq",
            "value": 8505987,
            "range": "± 564147",
            "unit": "ns/iter"
          },
          {
            "name": "memory_usage/parse_memory_usage",
            "value": 2592808,
            "range": "± 26200",
            "unit": "ns/iter"
          },
          {
            "name": "macro_sum_bench",
            "value": 1,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/100",
            "value": 90079,
            "range": "± 565",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/1000",
            "value": 879358,
            "range": "± 2113",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_parsing/parsing/10000",
            "value": 9523768,
            "range": "± 134244",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_iterator/iterator",
            "value": 616494,
            "range": "± 2464",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_writing/writing",
            "value": 3594898,
            "range": "± 53981",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/0",
            "value": 101,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/1",
            "value": 92,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/2",
            "value": 90,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/3",
            "value": 123,
            "range": "± 2",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/4",
            "value": 105,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/5",
            "value": 57,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/6",
            "value": 107,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_info_parsing/parsing/7",
            "value": 484,
            "range": "± 9",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/0",
            "value": 85,
            "range": "± 1",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/1",
            "value": 85,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/2",
            "value": 41,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/3",
            "value": 85,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/4",
            "value": 85,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/5",
            "value": 85,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/6",
            "value": 85,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_genotype_parsing/parsing/7",
            "value": 85,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/variant_type_detection",
            "value": 2534,
            "range": "± 66",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/allele_frequency_access",
            "value": 4545,
            "range": "± 30",
            "unit": "ns/iter"
          },
          {
            "name": "vcf_operations/sample_data_access",
            "value": 1375,
            "range": "± 10",
            "unit": "ns/iter"
          },
          {
            "name": "band_pass_normalized_kernel_gatk_defaults",
            "value": 987,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "band_pass_add_then_pop_ready_256_loci",
            "value": 292058,
            "range": "± 624",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_medium_depth_k10",
            "value": 1482115,
            "range": "± 4961",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_high_depth_k10",
            "value": 8926876,
            "range": "± 80847",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_depth/full_pipeline_low_k10/16",
            "value": 577266,
            "range": "± 2165",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_depth/threading_build_low_k10/16",
            "value": 396178,
            "range": "± 1120",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_depth/full_pipeline_low_k25/16",
            "value": 604911,
            "range": "± 4810",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_depth/threading_build_low_k25/16",
            "value": 419026,
            "range": "± 14180",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_depth/full_pipeline_medium_k10/64",
            "value": 1479106,
            "range": "± 4714",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_depth/threading_build_medium_k10/64",
            "value": 1191526,
            "range": "± 3683",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_depth/full_pipeline_medium_k25/64",
            "value": 1647994,
            "range": "± 3218",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_depth/threading_build_medium_k25/64",
            "value": 1291222,
            "range": "± 2572",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_depth/full_pipeline_high_k10/512",
            "value": 8938144,
            "range": "± 25546",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_depth/threading_build_high_k10/512",
            "value": 8568296,
            "range": "± 25598",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_depth/full_pipeline_high_k25/512",
            "value": 8634366,
            "range": "± 28070",
            "unit": "ns/iter"
          },
          {
            "name": "assembly_graph_depth/threading_build_high_k25/512",
            "value": 8200750,
            "range": "± 15599",
            "unit": "ns/iter"
          },
          {
            "name": "kmer_key_representation/A_arc_bytes/10",
            "value": 378791,
            "range": "± 1610",
            "unit": "ns/iter"
          },
          {
            "name": "kmer_key_representation/B_packed_key/10",
            "value": 332478,
            "range": "± 882",
            "unit": "ns/iter"
          },
          {
            "name": "kmer_key_representation/A_arc_bytes/25",
            "value": 420054,
            "range": "± 733",
            "unit": "ns/iter"
          },
          {
            "name": "kmer_key_representation/B_packed_key/25",
            "value": 338123,
            "range": "± 710",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/4",
            "value": 23820763,
            "range": "± 239234",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/16",
            "value": 95378287,
            "range": "± 1544248",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_scaffold_vectorized/batch/64",
            "value": 381968400,
            "range": "± 2034143",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/8",
            "value": 382012,
            "range": "± 945",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/8",
            "value": 129307,
            "range": "± 310",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/8",
            "value": 861515,
            "range": "± 2176",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/8",
            "value": 9457948,
            "range": "± 9143",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/32",
            "value": 1526873,
            "range": "± 9963",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/32",
            "value": 388357,
            "range": "± 1172",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/32",
            "value": 3427017,
            "range": "± 11291",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/32",
            "value": 37861195,
            "range": "± 136232",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/64",
            "value": 3053961,
            "range": "± 8531",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/64",
            "value": 617433,
            "range": "± 1709",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/64",
            "value": 6847184,
            "range": "± 19812",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/64",
            "value": 75694363,
            "range": "± 315290",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r100_h/128",
            "value": 6106752,
            "range": "± 22653",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r100_h/128",
            "value": 889853,
            "range": "± 2284",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r100_h/128",
            "value": 13694497,
            "range": "± 15125",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r100_h/128",
            "value": 151326890,
            "range": "± 241980",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/8",
            "value": 1521740,
            "range": "± 8522",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/8",
            "value": 523162,
            "range": "± 895",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/8",
            "value": 3451601,
            "range": "± 9916",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/8",
            "value": 38215394,
            "range": "± 87313",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/32",
            "value": 6086443,
            "range": "± 21052",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/32",
            "value": 1675474,
            "range": "± 4327",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/32",
            "value": 13798788,
            "range": "± 37369",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/32",
            "value": 152897384,
            "range": "± 466287",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/64",
            "value": 12171472,
            "range": "± 65986",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/64",
            "value": 2978533,
            "range": "± 5100",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/64",
            "value": 27591391,
            "range": "± 76321",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/64",
            "value": 306008662,
            "range": "± 1589749",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r200_h/128",
            "value": 24345346,
            "range": "± 64709",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r200_h/128",
            "value": 4776769,
            "range": "± 7885",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r200_h/128",
            "value": 55320505,
            "range": "± 159602",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r200_h/128",
            "value": 612294895,
            "range": "± 1685025",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/8",
            "value": 3407629,
            "range": "± 11900",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/8",
            "value": 1178866,
            "range": "± 7594",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/8",
            "value": 7741211,
            "range": "± 194645",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/8",
            "value": 86267255,
            "range": "± 207348",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/32",
            "value": 13626680,
            "range": "± 22208",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/32",
            "value": 3854177,
            "range": "± 5411",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/32",
            "value": 30938959,
            "range": "± 759643",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/32",
            "value": 345513224,
            "range": "± 745568",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/64",
            "value": 27262246,
            "range": "± 127874",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/64",
            "value": 7069327,
            "range": "± 170928",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/64",
            "value": 62025682,
            "range": "± 114284",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/64",
            "value": 690857385,
            "range": "± 2599553",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/logless_scalar_r300_h/128",
            "value": 54556849,
            "range": "± 253734",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_r300_h/128",
            "value": 12340170,
            "range": "± 55813",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/simd_f32_r300_h/128",
            "value": 123799518,
            "range": "± 253296",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_logless_simd/log10_scalar_r300_h/128",
            "value": 1381612229,
            "range": "± 4198741",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/logless_scalar_1x1/1",
            "value": 190372,
            "range": "± 494",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/simd_1x1/1",
            "value": 187876,
            "range": "± 246",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/wavefront_1x1/1",
            "value": 92242,
            "range": "± 171",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/logless_scalar_1x8/8",
            "value": 1521973,
            "range": "± 9440",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/simd_1x8/8",
            "value": 522482,
            "range": "± 507",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/wavefront_1x8/8",
            "value": 718409,
            "range": "± 1071",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/logless_scalar_1x16/16",
            "value": 3041234,
            "range": "± 3216",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/simd_1x16/16",
            "value": 925609,
            "range": "± 7714",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/wavefront_1x16/16",
            "value": 1433525,
            "range": "± 2588",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/logless_scalar_1x32/32",
            "value": 6082362,
            "range": "± 14724",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/simd_1x32/32",
            "value": 1680060,
            "range": "± 2501",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/wavefront_1x32/32",
            "value": 2849985,
            "range": "± 13418",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/logless_scalar_1x64/64",
            "value": 12166806,
            "range": "± 23776",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/simd_1x64/64",
            "value": 2984559,
            "range": "± 4648",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/wavefront_1x64/64",
            "value": 5697115,
            "range": "± 10177",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/logless_scalar_tile_8x32/8",
            "value": 48666164,
            "range": "± 96737",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/simd_tile_8x32/8",
            "value": 13429792,
            "range": "± 16934",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/wavefront_tile_8x32/8",
            "value": 22794384,
            "range": "± 25493",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/logless_scalar_tile_16x32/16",
            "value": 97325792,
            "range": "± 152333",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/simd_tile_16x32/16",
            "value": 26876282,
            "range": "± 28618",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/wavefront_tile_16x32/16",
            "value": 45598563,
            "range": "± 143649",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/logless_scalar_realistic_h/8",
            "value": 1669645,
            "range": "± 3822",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/simd_realistic_h/8",
            "value": 141422,
            "range": "± 266",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/wavefront_realistic_h/8",
            "value": 809815,
            "range": "± 8814",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/logless_scalar_realistic_h/32",
            "value": 6677812,
            "range": "± 8520",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/simd_realistic_h/32",
            "value": 171898,
            "range": "± 473",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/wavefront_realistic_h/32",
            "value": 3228534,
            "range": "± 7280",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/logless_scalar_realistic_h/64",
            "value": 13354887,
            "range": "± 59666",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/simd_realistic_h/64",
            "value": 212762,
            "range": "± 1489",
            "unit": "ns/iter"
          },
          {
            "name": "pairhmm_wavefront/wavefront_realistic_h/64",
            "value": 6419684,
            "range": "± 7819",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/value_copy_stages/3072x12",
            "value": 4888,
            "range": "± 7",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/arc_share_stages/3072x12",
            "value": 3458,
            "range": "± 9",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/assembly_result_set_arc_fanout/3072",
            "value": 53,
            "range": "± 0",
            "unit": "ns/iter"
          },
          {
            "name": "shared_reference_arc/assembly_result_set_value_fanout/3072",
            "value": 1194,
            "range": "± 14",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/64x48",
            "value": 11306,
            "range": "± 23",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/64x48",
            "value": 10994,
            "range": "± 1466",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/128x96",
            "value": 41729,
            "range": "± 94",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/128x96",
            "value": 40340,
            "range": "± 249",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/256x192",
            "value": 160338,
            "range": "± 483",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/256x192",
            "value": 155281,
            "range": "± 1410",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/200x151",
            "value": 99550,
            "range": "± 208",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/200x151",
            "value": 96319,
            "range": "± 574",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/soft_clip/100x100",
            "value": 34724,
            "range": "± 113",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_align/indel/100x100",
            "value": 33652,
            "range": "± 102",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_hc/padded_indel/80",
            "value": 35221,
            "range": "± 260",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_hc/padded_indel/150",
            "value": 96212,
            "range": "± 282",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_hc/padded_indel/250",
            "value": 236917,
            "range": "± 957",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_hc/read_to_hap_soft/hap120_read100",
            "value": 42987,
            "range": "± 104",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_hc/read_to_hap_soft/hap200_read151",
            "value": 103320,
            "range": "± 221",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_hc/read_to_hap_soft/hap280_read151",
            "value": 148821,
            "range": "± 242",
            "unit": "ns/iter"
          },
          {
            "name": "smith_waterman_hc/exact_substring_fast_path",
            "value": 212,
            "range": "± 0",
            "unit": "ns/iter"
          }
        ]
      }
    ]
  }
}