//! 6R.43 holdout: production `run_haplotype_caller` + assembly checkpoints vs frozen Java VCFs.
//!
//! Skipped unless `HOLDOUT_6R43=1`. Does not change production algorithms.
//!
//! ```text
//! HOLDOUT_6R43=1 cargo test -p gatk-haplotypecaller --test holdout_6r43_test -- --nocapture
//! ```

use gatk_common::config::GatkConfig;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, run_haplotype_caller,
    traverse_assembly_region_walker, AssemblyRegionCallDisposition, CallRegionArgs,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
};
use serde_json::{json, Value};
use std::fs;
use std::path::{Path, PathBuf};

const PANEL: &str = include_str!("../../scripts/parity/6r43_holdout_panel.json");

struct RegionSpec {
    id: String,
    interval: String,
    bam: String,
}

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn load_panel() -> Vec<RegionSpec> {
    let v: Value = serde_json::from_str(PANEL).expect("panel json");
    v["regions"]
        .as_array()
        .expect("regions")
        .iter()
        .map(|r| RegionSpec {
            id: r["id"].as_str().unwrap().to_string(),
            interval: r["interval"].as_str().unwrap().to_string(),
            bam: r["bam"].as_str().unwrap().to_string(),
        })
        .collect()
}

fn dump_checkpoints(id: &str, interval: &str, ref_fasta: &Path, bam: &Path, out_json: &Path) {
    let dict = SequenceDictionary::from_fasta_path(ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, interval).expect("interval");
    let walk = traverse_assembly_region_walker(
        &dict,
        &specs,
        ref_fasta,
        bam,
        &ReadFilterParams::gatk_standard_hc(),
        &WalkerTraversalConfig::gatk_haplotype_caller_production(100),
    )
    .expect("walk");
    let regions = flatten_assembly_regions(&walk);
    let args = CallRegionArgs::strict_java();
    let mut regions_json = Vec::new();
    for region in &regions {
        let disp = format!("{:?}", call_disposition(region));
        let mut row = json!({
            "contig": region.contig,
            "active": [region.start.get(), region.end.get()],
            "extended": [region.extended_start.get(), region.extended_end.get()],
            "n_reads": region.reads.len(),
            "disposition": disp,
        });
        if matches!(
            call_disposition(region),
            AssemblyRegionCallDisposition::ActiveFull
        ) {
            match HaplotypeCallerEngine::call_region(region, &dict, ref_fasta, &args) {
                Ok(Some(outcome)) => {
                    let k = outcome.assembly.kmer_size_for_dump();
                    let events: Vec<Value> = outcome
                        .assembly
                        .variation_events()
                        .iter()
                        .map(|e| {
                            json!({
                                "pos": e.start_1based.get(),
                                "ref": e.ref_allele,
                                "alt": e.alt_allele,
                            })
                        })
                        .collect();
                    let haps: Vec<Value> = outcome
                        .assembly
                        .haplotypes
                        .iter()
                        .map(|h| {
                            json!({
                                "is_ref": h.is_reference,
                                "len": h.bases.len(),
                                "k": h.kmer_size,
                                "loc": h.genome_loc.map(|g| [g.start.get(), g.end.get()]),
                            })
                        })
                        .collect();
                    row["kmer"] = json!(k);
                    row["n_haplotypes"] = json!(outcome.assembly.haplotypes.len());
                    row["n_events"] = json!(events.len());
                    row["events"] = json!(events);
                    row["haplotypes"] = json!(haps);
                    row["n_pairhmm_rows"] = json!(outcome.read_likelihoods.len());
                }
                Ok(None) => {
                    row["call_region"] = json!("None");
                }
                Err(e) => {
                    row["call_region_err"] = json!(e.to_string());
                }
            }
        }
        regions_json.push(row);
    }
    let doc = json!({
        "id": id,
        "interval": interval,
        "n_regions": regions.len(),
        "n_active_full": regions.iter().filter(|r| {
            matches!(call_disposition(r), AssemblyRegionCallDisposition::ActiveFull)
        }).count(),
        "regions": regions_json,
    });
    fs::write(out_json, serde_json::to_string_pretty(&doc).unwrap()).expect("write checkpoints");
}

#[test]
fn holdout_6r43_production_vcf_and_checkpoints() {
    if std::env::var("HOLDOUT_6R43").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R43=1 to run 6R.43 holdouts");
        return;
    }
    let root = repo_root();
    let ref_fasta = root.join("parity/realworld/assets/hs37d5.simple.fa");
    if !ref_fasta.is_file() {
        eprintln!("skip: missing {}", ref_fasta.display());
        return;
    }
    let out_root = root.join("parity/reports/6r43");
    fs::create_dir_all(&out_root).expect("out root");
    let only = std::env::var("HOLDOUT_ONLY").ok();
    for spec in load_panel() {
        if only.as_ref().is_some_and(|o| o != &spec.id) {
            continue;
        }
        let bam = root.join(&spec.bam);
        if !bam.is_file() {
            eprintln!("skip {}: missing BAM", spec.id);
            continue;
        }
        let dest_dir = out_root.join(&spec.id);
        fs::create_dir_all(&dest_dir).expect("region dir");
        let rust_vcf = dest_dir.join("rust.vcf");
        let mut cfg = GatkConfig::new("HaplotypeCaller".to_string());
        cfg.set_reference(ref_fasta.to_string_lossy().into_owned());
        cfg.add_input_file(bam.to_string_lossy().into_owned());
        cfg.set_output_vcf(rust_vcf.to_string_lossy().into_owned());
        cfg.set_parameter("intervals".to_string(), spec.interval.clone());
        eprintln!("[6R.43] RUST {} {}", spec.id, spec.interval);
        run_haplotype_caller(&cfg).unwrap_or_else(|e| panic!("{}: {e}", spec.id));
        dump_checkpoints(
            &spec.id,
            &spec.interval,
            &ref_fasta,
            &bam,
            &dest_dir.join("rust_checkpoints.json"),
        );
        assert!(
            rust_vcf.is_file(),
            "{} rust.vcf missing after run_haplotype_caller",
            spec.id
        );
    }
}

#[test]
fn holdout_6r43_panel_is_well_formed() {
    let specs = load_panel();
    assert!(specs.len() >= 8, "panel too small");
    assert!(specs.iter().any(|s| s.id == "ctrl_mid_b"));
    assert!(specs.iter().any(|s| s.id == "chr20_tiny"));
    assert!(specs.iter().any(|s| s.id == "chr21_w10"));
    assert!(specs.iter().any(|s| s.interval.starts_with("20:")));
    assert!(specs.iter().any(|s| s.interval.starts_with("21:")));
}

#[test]
fn holdout_6r43_java_fixture_has_rng_and_het_classes() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../parity/fixtures/p12-java-format/all_sites.tsv");
    let text = fs::read_to_string(path).expect("java format fixture");
    assert!(text.contains("0/1"), "need heterozygous Java sites");
    assert!(text.contains("TTC\tT"), "need deletion");
    assert!(text.contains("A\tATG"), "need insertion");
    assert!(text.contains("116.84"), "need high-QUAL jitter-class site");
}
