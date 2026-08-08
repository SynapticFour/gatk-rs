window.BENCHMARK_DATA = {
  "lastUpdate": 1786159510271,
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
      }
    ]
  }
}