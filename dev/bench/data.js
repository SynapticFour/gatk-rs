window.BENCHMARK_DATA = {
  "lastUpdate": 1785816600384,
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
      }
    ]
  }
}