//! gatk-rs CLI entry point (independent community project — not Broad GATK).
//! GATK-style flags are accepted for pipeline familiarity; this binary is native
//! Rust and does not launch the Broad GATK JVM.

use anyhow::Result;
use clap::{error::ErrorKind, Parser, Subcommand};
use gatk_common::{gatk_cli_exit_code, GatkConfig};
use gatk_core::io::{
    copy_alignments_with_htslib, count_records_in_region_indexed, qnames_in_region_indexed,
    validate_bam_file, FastaReader, VcfReader,
};
use gatk_core::reference::{
    count_acgtn_histogram_for_interval_specs, parse_intervals_cli_string, SequenceDictionary,
};
use gatk_core::variant_filtration::{
    gatk_indel_hard_filters, gatk_snp_hard_filters, run_variant_filtration, zip_filter_pairs,
    FilterSpec, VariantFiltrationArgs,
};
use gatk_haplotypecaller::{
    dump_smoothed_activity_tsv, init_worker_threads, passes_hc_read_filters_with_header,
    run_combine_gvcfs, run_genotype_gvcfs, run_haplotype_caller, CombineGvcfsArgs,
    GenotypeGvcfsArgs, ReadFilterParams, ReadHeaderSemantics, DEFAULT_STAND_CALL_CONF,
};
use rust_htslib::bam;
use rust_htslib::bam::Read as _;
use std::path::{Path, PathBuf};
use std::process;
use tracing::info;
use tracing_subscriber::{layer::SubscriberExt, util::SubscriberInitExt};
mod benchmarking;
use benchmarking::*;

const DISCLAIMER_HELP: &str = "\
DISCLAIMER: gatk-rs is an independent, community-driven reimplementation and is \
not affiliated with, endorsed by, or supported by the Broad Institute. \
\"GATK\" is a trademark of the Broad Institute; this project's name and branding \
will be revisited if requested. GATK-style flags (including --java-options) are \
for pipeline familiarity only — this binary is native Rust and does not start a JVM. \
See NOTICE.md and docs/CLAIM_MATRIX.md.";

#[derive(Parser)]
#[command(
    name = "gatk-rs",
    version = env!("CARGO_PKG_VERSION"),
    author = "GATK-RS Contributors",
    about = "Independent community Rust reimplementation of GATK-style tools (Alpha; scoped HC parity — not Broad GATK)",
    long_about = "gatk-rs is an independent, community-driven reimplementation and is not \
affiliated with, endorsed by, or supported by the Broad Institute. \
\"GATK\" is a trademark of the Broad Institute; this project's name and branding \
will be revisited if requested.\n\n\
HaplotypeCaller parity is validated on limited genomic regions and fixtures — \
see docs/CLAIM_MATRIX.md. Not a genome-wide drop-in or \
bitwise-identical-everywhere claim.\n\n\
GATK-style CLI flags exist for interoperability with existing pipelines; \
--java-options is accepted for familiarity and does not launch a JVM.",
    after_help = DISCLAIMER_HELP
)]
struct Cli {
    /// Accepted for familiarity with Broad GATK launchers only; native Rust — does not start a JVM (not Broad-affiliated)
    #[arg(long = "java-options", value_name = "OPTIONS")]
    java_options: Option<String>,

    /// Tool to execute
    #[command(subcommand)]
    tool: Tool,
}

#[derive(Subcommand)]
enum Tool {
    /// Call germline SNPs and indels via local re-assembly of haplotypes
    /// Default: `assembly-region-v1` emits variants via full `call_region`. Set
    /// `GATK_RS_HC_SCAFFOLD_OUTPUT=1` for header-only scaffold.
    #[command(
        visible_alias = "HaplotypeCaller",
        after_long_help = "Notes:\n\
- Independent community tool — not Broad GATK (see `gatk-rs --help` disclaimer / NOTICE.md).\n\
- Default pipeline is assembly-region-v1 (variant emission via call_region).\n\
- Set GATK_RS_HC_SCAFFOLD_OUTPUT=1 for header-only scaffold-v1 (Phase 9 golden tests).\n\
- GATK-style flags such as -ERC/--emit-ref-confidence are accepted for pipeline familiarity only.\n\
- `-alleles` / `--alleles` forces given alleles from a VCF (see docs/CLAIM_MATRIX.md T3-5).\n\
- Other deferred Java flags (-bamout, DRAGSTR, DRAGEN, --assembly-region-out): see docs/CLAIM_MATRIX.md\n"
    )]
    HaplotypeCaller {
        /// Reference genome FASTA (GATK: `-R`)
        #[arg(short = 'R', long)]
        reference: String,

        /// Input BAM/SAM (GATK: `-I`; repeat for multiple inputs)
        #[arg(short = 'I', long, num_args = 1..)]
        input: Vec<String>,

        /// Output VCF (GATK: `-O`)
        #[arg(short = 'O', long)]
        output: String,

        /// Genomic intervals (GATK: `-L`)
        #[arg(short = 'L', long)]
        intervals: Option<String>,

        /// Variant calling output mode: VCF | GVCF | BP_RESOLUTION (GATK: `--output-mode` / `-mode`)
        #[arg(long = "output-mode", alias = "mode", default_value = "VCF")]
        output_mode: String,

        /// Minimum base quality (GATK: `--min-base-quality-score`)
        #[arg(long = "min-base-quality-score", default_value = "10")]
        min_base_quality_score: u8,

        /// Emit reference confidence: NONE | GVCF | BP_RESOLUTION (GATK: `--emit-ref-confidence`; native GATK also documents `-ERC`, which is not expressible as a single clap short flag — use the long form here).
        #[arg(long = "emit-ref-confidence")]
        emit_ref_confidence: Option<String>,

        /// PairHMM implementation (GATK: `--pair-hmm-implementation`):
        /// `LOG10_PAIRHMM` (default) | `LOGLESS_HMM` | `AVX`/`SIMD` | `SIMD_F32` | `FASTEST_AVAILABLE`.
        /// Alias: `--pairhmm-impl scalar|logless|simd|fastest`.
        #[arg(
            long = "pair-hmm",
            visible_aliases = ["pair-hmm-implementation", "pairhmm-impl"]
        )]
        pair_hmm: Option<String>,

        /// Use original base qualities
        #[arg(long)]
        original_base_qualities: bool,

        /// Don't use soft clipped bases
        #[arg(long)]
        dont_use_soft_clipped_bases: bool,

        /// Maximum number of alternative haplotypes (not wired — Java default 6 always used)
        #[arg(long, default_value = "6", hide = true)]
        max_alternate_alleles: u32,

        /// Minimum mapping quality (GATK: `--min-mapping-quality`)
        #[arg(long = "min-mapping-quality", default_value = "20")]
        min_mapping_quality: u32,

        /// Minimum phred-scaled confidence to call variants (GATK: `--standard-min-confidence-threshold-for-calling`)
        #[arg(
            long = "standard-min-confidence-threshold-for-calling",
            visible_alias = "stand-call-conf",
            default_value = "30.0"
        )]
        stand_call_confidence: f64,

        /// Minimum phred-scaled confidence to emit variants (GATK: `--standard-min-confidence-threshold-for-emitting`)
        #[arg(
            long = "standard-min-confidence-threshold-for-emitting",
            visible_alias = "stand-emit-conf",
            default_value = "10.0"
        )]
        stand_emit_confidence: f64,

        /// VCF of alleles to force into assembly (GATK: `-alleles` / `--alleles`)
        #[arg(short = 'A', long = "alleles", value_name = "VCF")]
        alleles: Option<String>,

        /// Worker threads for Active-Region / PairHMM parallelism (GATK-shaped `--threads` / `-nt`).
        /// Default: number of logical CPUs. `0` also means auto-detect.
        /// When set, wins over ambient `RAYON_NUM_THREADS`.
        #[arg(short = 't', long = "threads", visible_alias = "nt", value_name = "N")]
        threads: Option<usize>,
    },

    /// PrintReads - Print reads from SAM/BAM file to output
    #[command(visible_alias = "PrintReads")]
    PrintReads {
        /// Input BAM/SAM file
        #[arg(short = 'I', long)]
        input: String,

        /// Output BAM/SAM file
        #[arg(short = 'O', long)]
        output: String,
    },

    /// Combine per-sample gVCFs into a multi-sample gVCF (no joint genotyping).
    /// Merges reference-confidence blocks and remaps PL/AD onto a unified allele
    /// set (GATK CombineGVCFs algorithm). GenotypeGVCFs is a separate step.
    #[command(visible_alias = "CombineGVCFs")]
    CombineGvcfs {
        /// Reference genome FASTA (GATK: `-R`)
        #[arg(short = 'R', long)]
        reference: String,

        /// Input gVCF(s); repeat `-V` once per sample (GATK: `-V`)
        #[arg(short = 'V', long = "variant", required = true)]
        variant: Vec<String>,

        /// Output combined gVCF (GATK: `-O`)
        #[arg(short = 'O', long)]
        output: String,

        /// Optional intervals (GATK: `-L`)
        #[arg(short = 'L', long)]
        intervals: Option<String>,
    },

    /// Hard-filter variants with GATK-style filter expressions (not VQSR).
    /// Soft-filters by writing FILTER tags. Does **not** replace VQSR algorithmically
    /// it is the pragmatic fallback for smaller cohorts where VQSR cannot be trained
    /// cleanly (aligned with official GATK Best Practices guidance).
    /// Prefer one annotation per `--filter-expression` (GATK recommendation). Use
    /// `--preset snp` or `--preset indel` for the official hard-filter tables, or pass
    /// explicit `--filter-expression` / `--filter-name` pairs (Java GATK CLI compatible).
    #[command(visible_alias = "VariantFiltration")]
    VariantFiltration {
        /// Input VCF (GATK: `-V`)
        #[arg(short = 'V', long = "variant")]
        variant: String,

        /// Output VCF (GATK: `-O`)
        #[arg(short = 'O', long)]
        output: String,

        /// Optional reference (accepted for GATK familiarity; unused for hard filters)
        #[arg(short = 'R', long)]
        reference: Option<String>,

        /// Filter expression on INFO/QUAL (repeatable; GATK: `--filter-expression` / `-filter`)
        #[arg(long = "filter-expression", visible_alias = "filter")]
        filter_expression: Vec<String>,

        /// Name written to FILTER when the paired expression fails (repeatable)
        #[arg(long = "filter-name")]
        filter_name: Vec<String>,

        /// Apply official GATK hard-filter table: `snp` or `indel`
        /// (merged with any explicit --filter-expression pairs)
        #[arg(long = "preset")]
        preset: Option<String>,
    },

    /// Joint-genotype a multi-sample gVCF into a final cohort VCF.
    /// Uses cohort allele-frequency EM over sample PLs (GATK GenotypeGVCFs).
    #[command(visible_alias = "GenotypeGVCFs")]
    GenotypeGvcfs {
        /// Reference genome FASTA (GATK: `-R`)
        #[arg(short = 'R', long)]
        reference: String,

        /// Combined multi-sample gVCF (GATK: `-V`)
        #[arg(short = 'V', long = "variant")]
        variant: String,

        /// Output cohort VCF (GATK: `-O`)
        #[arg(short = 'O', long)]
        output: String,

        /// Optional intervals (GATK: `-L`)
        #[arg(short = 'L', long)]
        intervals: Option<String>,

        /// Minimum QUAL to emit (GATK: `--standard-min-confidence-threshold-for-calling`)
        #[arg(long = "stand-call-conf", default_value_t = DEFAULT_STAND_CALL_CONF)]
        stand_call_conf: f64,
    },

    /// Dump band-pass smoothed per-base activity (for parity with GATK `--assembly-region-out` analysis).
    #[command(visible_alias = "DumpSmoothedActivity")]
    DumpSmoothedActivity {
        #[arg(short = 'R', long)]
        reference: String,
        #[arg(short = 'I', long)]
        input: String,
        #[arg(short = 'L', long)]
        intervals: String,
        #[arg(short = 'O', long)]
        output: String,
    },

    /// FilterReads - Apply HC-style ingress read filters and emit filtered SAM/BAM.
    #[command(visible_alias = "FilterReads")]
    FilterReads {
        /// Input BAM/SAM file
        #[arg(short = 'I', long)]
        input: String,

        /// Output BAM/SAM file
        #[arg(short = 'O', long)]
        output: String,

        /// Optional 1-based inclusive region `contig:start-end` (requires BAI; matches GATK `-L` subsetting)
        #[arg(short = 'L', long)]
        interval: Option<String>,

        /// Minimum mapping quality (MQ 255 passes as unknown, matching HC semantics)
        #[arg(long, default_value = "20")]
        min_mapping_quality: u8,

        /// Include duplicate reads
        #[arg(long, default_value_t = false)]
        include_duplicates: bool,
    },

    /// List all available tools
    #[command(name = "--list")]
    ListTools,

    /// Print tool version information
    #[command(name = "--version")]
    Version,

    /// Count A/C/G/T/N bases in the reference over intervals (GATK CountBasesInReference subset).
    #[command(visible_alias = "CountBasesInReference")]
    CountBasesInReference {
        #[arg(short = 'R', long)]
        reference: String,
        #[arg(short = 'L', long)]
        intervals: Option<String>,
    },

    /// Count reads in one indexed region (BAM+BAI required).
    #[command(visible_alias = "CountReadsInRegion")]
    CountReadsInRegion {
        #[arg(short = 'I', long)]
        input: String,
        #[arg(short = 'L', long)]
        region: String,
    },

    /// List read QNAMEs in one indexed region (BAM+BAI required).
    #[command(visible_alias = "ListReadsInRegion")]
    ListReadsInRegion {
        #[arg(short = 'I', long)]
        input: String,
        #[arg(short = 'L', long)]
        region: String,
    },

    /// Validate input files
    #[command(visible_alias = "Validate")]
    Validate {
        /// Input file to validate
        input: String,

        /// File type (BAM, SAM, VCF, FASTA)
        #[arg(short = 't', long)]
        file_type: String,

        /// Optional reference FASTA for dictionary compatibility checks
        #[arg(short = 'R', long)]
        reference: Option<String>,
    },

    /// Benchmarking commands
    Benchmark {
        #[command(flatten)]
        args: BenchmarkingArgs,
    },
}

/// Resolve HC worker threads and initialize Rayon's global pool (via haplotypecaller).
/// Explicit `--threads N` / `--nt N` (N>0) wins over ambient `RAYON_NUM_THREADS`
/// `--threads 0` → auto-detect (`num_cpus`)
/// Omitted → honor `RAYON_NUM_THREADS` if set, else `num_cpus`
fn init_haplotype_caller_threads(threads: Option<usize>) -> usize {
    let n = match threads {
        Some(0) => num_cpus::get().max(1),
        Some(n) => n.max(1),
        None => std::env::var("RAYON_NUM_THREADS")
            .ok()
            .and_then(|s| s.parse::<usize>().ok())
            .filter(|&v| v > 0)
            .unwrap_or_else(|| num_cpus::get().max(1)),
    };
    init_worker_threads(n, threads.is_some());
    info!("HaplotypeCaller worker threads: {n}");
    n
}

fn main() {
    // Initialize logging
    init_logging();

    let cli = match Cli::try_parse() {
        Ok(cli) => cli,
        Err(err) => {
            let _ = err.print();
            let code = match err.kind() {
                ErrorKind::DisplayHelp | ErrorKind::DisplayVersion => 0,
                _ => 2,
            };
            process::exit(code);
        }
    };

    // Compatibility shim: do not imply a JVM / Broad GATK launch.
    if cli.java_options.is_some() {
        eprintln!(
            "note: --java-options is accepted for Broad GATK launcher familiarity only; \
gatk-rs is a native Rust binary and does not start a JVM \
(not affiliated with, endorsed by, or supported by the Broad Institute)."
        );
    }

    let result: Result<()> = (|| {
        match cli.tool {
            Tool::HaplotypeCaller {
                reference,
                input,
                output,
                intervals,
                output_mode,
                min_base_quality_score,
                emit_ref_confidence,
                pair_hmm,
                original_base_qualities,
                dont_use_soft_clipped_bases,
                max_alternate_alleles,
                min_mapping_quality,
                stand_call_confidence,
                stand_emit_confidence,
                alleles,
                threads,
            } => {
                let hc_run: gatk_common::GatkResult<()> = (|| {
                    let n_threads = init_haplotype_caller_threads(threads);
                    let mut config = match &cli.java_options {
                        Some(java_opts) => {
                            GatkConfig::with_java_options("HaplotypeCaller".to_string(), java_opts)?
                        }
                        None => GatkConfig::new("HaplotypeCaller".to_string()),
                    };
                    config.global_config.num_threads = Some(n_threads);

                    for input_file in &input {
                        config.add_input_file(input_file.clone());
                    }
                    config.set_reference(reference.clone());
                    config.set_output_vcf(output.clone());

                    if let Some(intervals_val) = &intervals {
                        let dict = SequenceDictionary::from_fasta_path(&reference)?;
                        let _specs = parse_intervals_cli_string(&dict, intervals_val)?;
                        config.set_parameter("intervals".to_string(), intervals_val.clone());
                    }
                    config.set_output_mode(output_mode.clone());
                    config.set_parameter("output_mode".to_string(), output_mode);
                    config.set_parameter(
                        "min_base_quality_score".to_string(),
                        min_base_quality_score.to_string(),
                    );
                    config.set_parameter(
                        "min_mapping_quality".to_string(),
                        min_mapping_quality.to_string(),
                    );
                    config.set_parameter(
                        "max_alternate_alleles".to_string(),
                        max_alternate_alleles.to_string(),
                    );
                    config.set_parameter(
                        "stand_call_confidence".to_string(),
                        stand_call_confidence.to_string(),
                    );
                    config.set_parameter(
                        "stand_emit_confidence".to_string(),
                        stand_emit_confidence.to_string(),
                    );

                    if let Some(alleles_vcf) = alleles {
                        config.set_parameter("alleles".to_string(), alleles_vcf);
                    }

                    if let Some(emit_ref_conf) = emit_ref_confidence {
                        config.set_parameter(
                            "emit_ref_confidence".to_string(),
                            emit_ref_conf.clone(),
                        );
                        // GATK -ERC GVCF / BP_RESOLUTION implies matching output mode when not overridden.
                        match emit_ref_conf.as_str() {
                            "GVCF" => config.set_output_mode("GVCF".to_string()),
                            "BP_RESOLUTION" => config.set_output_mode("BP_RESOLUTION".to_string()),
                            _ => {}
                        }
                    }
                    if let Some(pair_hmm_val) = pair_hmm {
                        config.set_parameter("pair_hmm".to_string(), pair_hmm_val);
                    }
                    if original_base_qualities {
                        config.set_parameter(
                            "original_base_qualities".to_string(),
                            "true".to_string(),
                        );
                    }
                    if dont_use_soft_clipped_bases {
                        config.set_parameter(
                            "dont_use_soft_clipped_bases".to_string(),
                            "true".to_string(),
                        );
                    }

                    config.validate()?;

                    info!("Starting HaplotypeCaller");
                    info!("Reference: {}", reference);
                    info!("Input files: {:?}", input);
                    info!("Output: {}", output);

                    if let Some(ref java_opts) = config.java_options {
                        info!("Java memory: {:?}", java_opts.memory);
                        info!("Java GC: {:?}", java_opts.garbage_collector);
                        if !java_opts.additional_args.is_empty() {
                            info!("Additional Java args: {:?}", java_opts.additional_args);
                        }
                    }

                    run_haplotype_caller(&config)
                })();

                match hc_run {
                    Ok(()) => Ok(()),
                    Err(e) => {
                        eprintln!("{}", e.display_with_context());
                        process::exit(gatk_cli_exit_code(&e));
                    }
                }
            }

            Tool::CountBasesInReference {
                reference,
                intervals,
            } => {
                let dict = SequenceDictionary::from_fasta_path(&reference).map_err(|e| {
                    anyhow::anyhow!("Failed to derive sequence dictionary from reference: {e}")
                })?;
                let specs = match &intervals {
                    None => dict.whole_genome_interval_specs(),
                    Some(s) => parse_intervals_cli_string(&dict, s)
                        .map_err(|e| anyhow::anyhow!("Invalid intervals: {e}"))?,
                };
                let counts = count_acgtn_histogram_for_interval_specs(&reference, &dict, &specs)
                    .map_err(|e| anyhow::anyhow!("Count failed: {e}"))?;
                let labels = ['A', 'C', 'G', 'T', 'N'];
                for (i, &n) in counts.iter().enumerate() {
                    if i == 4 && n == 0 {
                        continue;
                    }
                    println!("{} : {}", labels[i], n);
                }
                Ok(())
            }

            Tool::CountReadsInRegion { input, region } => {
                let (contig, start, end) = parse_single_region(&region)?;
                let count = count_records_in_region_indexed(
                    std::path::Path::new(&input),
                    &contig,
                    start,
                    end,
                )
                .map_err(|e| anyhow::anyhow!("Indexed region count failed: {e}"))?;
                println!("COUNT : {}", count);
                Ok(())
            }

            Tool::ListReadsInRegion { input, region } => {
                let (contig, start, end) = parse_single_region(&region)?;
                let qnames =
                    qnames_in_region_indexed(std::path::Path::new(&input), &contig, start, end)
                        .map_err(|e| anyhow::anyhow!("Indexed region list failed: {e}"))?;
                for q in qnames {
                    println!("{q}");
                }
                Ok(())
            }

            Tool::PrintReads { input, output } => {
                info!("Starting PrintReads");
                info!("Input: {}", input);
                info!("Output: {}", output);

                let n = copy_alignments_with_htslib(
                    std::path::Path::new(&input),
                    std::path::Path::new(&output),
                )?;
                info!("PrintReads copied {} records", n);
                Ok(())
            }

            Tool::CombineGvcfs {
                reference,
                variant,
                output,
                intervals,
            } => {
                info!("Starting CombineGVCFs");
                info!("Reference: {}", reference);
                info!("Inputs: {:?}", variant);
                info!("Output: {}", output);
                run_combine_gvcfs(&CombineGvcfsArgs {
                    reference: PathBuf::from(reference),
                    variant_paths: variant.into_iter().map(PathBuf::from).collect(),
                    output: PathBuf::from(output),
                    intervals,
                })
                .map_err(|e| anyhow::anyhow!("{e}"))?;
                Ok(())
            }

            Tool::VariantFiltration {
                variant,
                output,
                reference,
                filter_expression,
                filter_name,
                preset,
            } => {
                info!("Starting VariantFiltration");
                let mut filters: Vec<FilterSpec> = Vec::new();
                if let Some(p) = preset.as_deref() {
                    match p.to_ascii_lowercase().as_str() {
                        "snp" | "snps" => filters.extend(gatk_snp_hard_filters()),
                        "indel" | "indels" => filters.extend(gatk_indel_hard_filters()),
                        other => {
                            return Err(anyhow::anyhow!(
                                "unknown --preset '{other}' (use snp or indel)"
                            ));
                        }
                    }
                }
                if !filter_expression.is_empty() || !filter_name.is_empty() {
                    filters.extend(
                        zip_filter_pairs(&filter_expression, &filter_name)
                            .map_err(|e| anyhow::anyhow!("{e}"))?,
                    );
                }
                if filters.is_empty() {
                    return Err(anyhow::anyhow!(
                        "VariantFiltration requires --preset snp|indel and/or \
                         matching --filter-expression / --filter-name pairs"
                    ));
                }
                run_variant_filtration(&VariantFiltrationArgs {
                    variant: PathBuf::from(variant),
                    output: PathBuf::from(output),
                    filters,
                    reference: reference.map(PathBuf::from),
                })
                .map_err(|e| anyhow::anyhow!("{e}"))?;
                Ok(())
            }

            Tool::GenotypeGvcfs {
                reference,
                variant,
                output,
                intervals,
                stand_call_conf,
            } => {
                info!("Starting GenotypeGVCFs");
                info!("Reference: {}", reference);
                info!("Input: {}", variant);
                info!("Output: {}", output);
                run_genotype_gvcfs(&GenotypeGvcfsArgs {
                    reference: PathBuf::from(reference),
                    variant: PathBuf::from(variant),
                    output: PathBuf::from(output),
                    intervals,
                    stand_call_conf,
                    include_non_variant_sites: false,
                })
                .map_err(|e| anyhow::anyhow!("{e}"))?;
                Ok(())
            }

            Tool::DumpSmoothedActivity {
                reference,
                input,
                intervals,
                output,
            } => {
                info!("Starting DumpSmoothedActivity");
                let p = ReadFilterParams::default();
                dump_smoothed_activity_tsv(
                    Path::new(&reference),
                    Path::new(&input),
                    &intervals,
                    Path::new(&output),
                    &p,
                )
                .map_err(|e| anyhow::anyhow!("{e}"))?;
                info!("Wrote smoothed activity TSV to {}", output);
                Ok(())
            }

            Tool::FilterReads {
                input,
                output,
                interval,
                min_mapping_quality,
                include_duplicates,
            } => {
                info!("Starting FilterReads");
                info!("Input: {}", input);
                info!("Output: {}", output);
                let params = ReadFilterParams {
                    min_mapping_quality,
                    exclude_duplicates: !include_duplicates,
                    exclude_secondary: true,
                    exclude_supplementary: true,
                };
                let out_format = if output.to_ascii_lowercase().ends_with(".sam") {
                    bam::Format::Sam
                } else {
                    bam::Format::Bam
                };
                let mut kept = 0usize;
                if let Some(reg) = interval {
                    let (contig, start_1, end_1) = parse_single_region(&reg)?;
                    let mut reader = bam::IndexedReader::from_path(&input)
                        .map_err(|e| anyhow::anyhow!("Indexed BAM required for -L: {e}"))?;
                    let header = bam::Header::from_template(reader.header());
                    let tid = reader.header().tid(contig.as_bytes()).ok_or_else(|| {
                        anyhow::anyhow!("Contig {contig} not found in BAM header")
                    })?;
                    let start0 = start_1 - 1;
                    let end0 = end_1;
                    reader
                        .fetch((tid, start0 as i64, end0 as i64))
                        .map_err(|e| anyhow::anyhow!("fetch failed: {e}"))?;
                    let hdr_view = reader.header().clone();
                    let mut writer =
                        bam::Writer::from_path(&output, &header, out_format).map_err(|e| {
                            anyhow::anyhow!("Failed to create output alignment file: {e}")
                        })?;
                    for rec in reader.records() {
                        let rec =
                            rec.map_err(|e| anyhow::anyhow!("Failed reading input record: {e}"))?;
                        if passes_hc_read_filters_with_header(&rec, &hdr_view, &params) {
                            writer.write(&rec).map_err(|e| {
                                anyhow::anyhow!("Failed writing output record: {e}")
                            })?;
                            kept += 1;
                        }
                    }
                } else {
                    let mut reader = bam::Reader::from_path(&input)
                        .map_err(|e| anyhow::anyhow!("Failed to open input alignment file: {e}"))?;
                    let header = bam::Header::from_template(reader.header());
                    let hdr_view = reader.header().clone();
                    let mut writer =
                        bam::Writer::from_path(&output, &header, out_format).map_err(|e| {
                            anyhow::anyhow!("Failed to create output alignment file: {e}")
                        })?;
                    for rec in reader.records() {
                        let rec =
                            rec.map_err(|e| anyhow::anyhow!("Failed reading input record: {e}"))?;
                        if passes_hc_read_filters_with_header(&rec, &hdr_view, &params) {
                            writer.write(&rec).map_err(|e| {
                                anyhow::anyhow!("Failed writing output record: {e}")
                            })?;
                            kept += 1;
                        }
                    }
                }
                println!("FilterReads kept {} records", kept);
                Ok(())
            }

            Tool::ListTools => {
                println!("Available tools:");
                println!("  HaplotypeCaller        - Call germline SNPs and indels via local re-assembly of haplotypes");
                println!(
                    "  CombineGVCFs           - Merge per-sample gVCFs into a multi-sample gVCF"
                );
                println!(
                    "  GenotypeGVCFs          - Joint-genotype a multi-sample gVCF into a cohort VCF"
                );
                println!(
                    "  VariantFiltration     - Hard-filter variants (GATK expressions; not VQSR)"
                );
                println!(
                    "  CountBasesInReference  - Count A/C/G/T/N in the reference over intervals"
                );
                println!("  CountReadsInRegion     - Count reads in one indexed region");
                println!("  ListReadsInRegion      - List read names in one indexed region");
                println!("  PrintReads            - Print reads from SAM/BAM file to output");
                println!(
                    "  DumpSmoothedActivity  - Per-base smoothed activity TSV (assembly-region parity)"
                );
                println!("  FilterReads           - Apply HC-style ingress read filters");
                println!("  Validate              - Validate input file formats");
                println!();
                println!("More tools will be added as implementation progresses.");

                Ok(())
            }

            Tool::Version => {
                // Do not impersonate Broad GATK. Expose our version and the pinned
                // Java oracle version used for differential parity testing only.
                println!(
                    "gatk-rs {} (independent community project — not Broad Institute GATK)",
                    env!("CARGO_PKG_VERSION")
                );
                println!(
                    "Pinned Java GATK oracle for parity tests: 4.4.0.0 (see GATK_PINNED_SHA / NOTICE.md)"
                );
                println!("{DISCLAIMER_HELP}");

                Ok(())
            }

            Tool::Validate {
                input,
                file_type,
                reference,
            } => {
                info!("Validating {} file: {}", file_type, input);
                let dictionary = if let Some(reference_path) = &reference {
                    Some(
                        SequenceDictionary::from_fasta_path(reference_path).map_err(|e| {
                            anyhow::anyhow!("Reference dictionary creation failed: {e}")
                        })?,
                    )
                } else {
                    None
                };

                match file_type.to_uppercase().as_str() {
                    "BAM" => {
                        let count =
                            validate_bam_file(std::path::Path::new(&input), dictionary.as_ref())
                                .map_err(|e| anyhow::anyhow!("BAM validation failed: {e}"))?;
                        println!(
                            "BAM validation passed for file: {} ({} records)",
                            input, count
                        );
                    }
                    "SAM" => {
                        let mut reader =
                            gatk_core::io::SamReader::from_file(&input).map_err(|e| {
                                anyhow::anyhow!("SAM validation failed to open file: {e}")
                            })?;
                        if let Some(dict) = &dictionary {
                            dict.validate_sam_header(reader.header()).map_err(|e| {
                                anyhow::anyhow!("SAM/reference dictionary mismatch: {e}")
                            })?;
                        }
                        let semantics = ReadHeaderSemantics::from_sam_header_text(
                            &sam_header_to_text(reader.header()),
                        )
                        .map_err(|e| anyhow::anyhow!("SAM header semantics invalid: {e}"))?;

                        let mut records_seen = 0usize;
                        while let Some(record) = reader.read_next_record().map_err(|e| {
                            anyhow::anyhow!("SAM validation failed while reading: {e}")
                        })? {
                            records_seen += 1;
                            if record.cigar != "*" && record.seq != "*" {
                                let expected_seq_bases =
                                    query_bases_consumed_by_cigar(&record.cigar).map_err(|e| {
                                        anyhow::anyhow!(
                                            "SAM validation failed for record {}: {e}",
                                            record.qname
                                        )
                                    })?;
                                if expected_seq_bases != record.seq.len() {
                                    return Err(anyhow::anyhow!(
                                        "SAM validation failed for record {}: CIGAR/SEQ length mismatch (expected {} query bases from CIGAR '{}', observed SEQ length {})",
                                        record.qname,
                                        expected_seq_bases,
                                        record.cigar,
                                        record.seq.len()
                                    ));
                                }
                            }
                            let rg_tag = record.get_optional_value("RG");
                            let pg_tag = record.get_optional_value("PG");
                            semantics
                                .validate_record_links(rg_tag.as_deref(), pg_tag.as_deref())
                                .map_err(|e| {
                                    anyhow::anyhow!(
                                        "SAM read-header semantics failed for record {}: {e}",
                                        record.qname
                                    )
                                })?;
                        }
                        if records_seen == 0 {
                            return Err(anyhow::anyhow!(
                                "SAM validation failed while reading: no alignment records found"
                            ));
                        }
                        println!("SAM validation passed for file: {}", input);
                    }
                    "VCF" => {
                        let mut reader = VcfReader::from_file(&input).map_err(|e| {
                            anyhow::anyhow!("VCF validation failed to open file: {e}")
                        })?;
                        if let Some(dict) = &dictionary {
                            dict.validate_vcf_header(reader.header()).map_err(|e| {
                                anyhow::anyhow!("VCF/reference dictionary mismatch: {e}")
                            })?;
                        }
                        let _ = reader.read_next_record().map_err(|e| {
                            anyhow::anyhow!("VCF validation failed while reading: {e}")
                        })?;
                        println!("VCF validation passed for file: {}", input);
                    }
                    "FASTA" => {
                        let mut reader = FastaReader::from_file(&input).map_err(|e| {
                            anyhow::anyhow!("FASTA validation failed to open file: {e}")
                        })?;
                        let _ = reader.read_next_sequence().map_err(|e| {
                            anyhow::anyhow!("FASTA validation failed while reading: {e}")
                        })?;
                        println!("FASTA validation passed for file: {}", input);
                    }
                    _ => {
                        return Err(anyhow::anyhow!("Unsupported file type: {}", file_type));
                    }
                }

                Ok(())
            }

            Tool::Benchmark { args } => {
                info!("Running benchmarking command");
                run_benchmarking_command(args)?;
                Ok(())
            }
        }
    })();

    if let Err(e) = result {
        eprintln!("{e}");
        // Match GATK user-facing failure convention.
        process::exit(2);
    }
}

fn init_logging() {
    tracing_subscriber::registry()
        .with(
            tracing_subscriber::EnvFilter::try_from_default_env().unwrap_or_else(|_| "info".into()),
        )
        .with(tracing_subscriber::fmt::layer())
        .init();
}

fn sam_header_to_text(header: &gatk_core::io::SamHeader) -> String {
    let mut lines = Vec::new();
    lines.push("@HD\tVN:1.6".to_string());
    for rg in &header.read_groups {
        let mut line = format!("@RG\tID:{}", rg.id);
        if let Some(sm) = &rg.sample {
            line.push_str(&format!("\tSM:{sm}"));
        }
        lines.push(line);
    }
    for pg in &header.programs {
        lines.push(format!("@PG\tID:{}", pg.id));
    }
    lines.join("\n")
}

fn parse_single_region(region: &str) -> Result<(String, u64, u64)> {
    let (contig, span) = region
        .split_once(':')
        .ok_or_else(|| anyhow::anyhow!("Invalid region '{}': expected contig:start-end", region))?;
    let span = span.replace(',', "");
    let (start_s, end_s) = span
        .split_once('-')
        .ok_or_else(|| anyhow::anyhow!("Invalid region '{}': expected start-end", region))?;
    let start = start_s
        .parse::<u64>()
        .map_err(|_| anyhow::anyhow!("Invalid region start '{}'", start_s))?;
    let end = end_s
        .parse::<u64>()
        .map_err(|_| anyhow::anyhow!("Invalid region end '{}'", end_s))?;
    if start == 0 {
        return Err(anyhow::anyhow!(
            "Invalid region '{}': start must be >= 1",
            region
        ));
    }
    if start > end {
        return Err(anyhow::anyhow!(
            "Invalid region '{}': start must be <= end",
            region
        ));
    }
    Ok((contig.to_string(), start, end))
}

fn query_bases_consumed_by_cigar(cigar: &str) -> Result<usize> {
    let mut total = 0usize;
    let mut num = String::new();
    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            num.push(ch);
            continue;
        }
        if num.is_empty() {
            return Err(anyhow::anyhow!(
                "Invalid CIGAR '{}': missing length before op '{}'",
                cigar,
                ch
            ));
        }
        let len: usize = num
            .parse()
            .map_err(|_| anyhow::anyhow!("Invalid CIGAR '{}': bad op length '{}'", cigar, num))?;
        match ch {
            'M' | 'I' | 'S' | '=' | 'X' => total += len,
            'D' | 'N' | 'H' | 'P' => {}
            _ => {
                return Err(anyhow::anyhow!(
                    "Invalid CIGAR '{}': unknown op '{}'",
                    cigar,
                    ch
                ))
            }
        }
        num.clear();
    }
    if !num.is_empty() {
        return Err(anyhow::anyhow!(
            "Invalid CIGAR '{}': trailing length without op '{}'",
            cigar,
            num
        ));
    }
    Ok(total)
}
