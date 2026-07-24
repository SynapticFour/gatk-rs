//! Warm-up path: real fixtures, scaffold + default assembly-region emission.

use gatk_common::GatkConfig;
use gatk_haplotypecaller::run_haplotype_caller;
use std::fs;
use std::path::PathBuf;
use std::sync::{Mutex, OnceLock};

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../parity/fixtures")
        .join(name)
}

fn hc_env_lock() -> &'static Mutex<()> {
    static LOCK: OnceLock<Mutex<()>> = OnceLock::new();
    LOCK.get_or_init(|| Mutex::new(()))
}

fn clear_hc_output_env() {
    std::env::remove_var("GATK_RS_HC_SCAFFOLD_OUTPUT");
    std::env::remove_var("GATK_RS_HC_ACTIVATE_OUTPUT");
}

#[test]
fn haplotype_caller_warmup_writes_scaffold_vcf() {
    let _guard = hc_env_lock().lock().expect("lock");
    clear_hc_output_env();
    std::env::set_var("GATK_RS_HC_SCAFFOLD_OUTPUT", "1");
    let ref_fa = fixture("reference.fa");
    let bam = fixture("sample.bam");
    let out = std::env::temp_dir().join("gatk_rs_hc_warmup_out.vcf");

    let mut config = GatkConfig::new("HaplotypeCaller".to_string());
    config.set_reference(ref_fa.to_string_lossy().to_string());
    config.add_input_file(bam.to_string_lossy().to_string());
    config.set_output_vcf(out.to_string_lossy().to_string());
    config.set_parameter("intervals".to_string(), "chr1:1-32".to_string());
    config.validate().expect("valid config");

    run_haplotype_caller(&config).expect("scaffold run");
    clear_hc_output_env();
    let text = fs::read_to_string(&out).expect("read out vcf");
    assert!(text.contains("##fileformat=VCFv4.2"));
    assert!(text.contains("##GATK_RS_HC_PIPELINE=scaffold-v1"));
    assert!(text.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"));
}

#[test]
fn haplotype_caller_default_emits_assembly_region_pipeline() {
    let _guard = hc_env_lock().lock().expect("lock");
    clear_hc_output_env();
    let ref_fa = fixture("reference.fa");
    let bam = fixture("sample.bam");
    let out = std::env::temp_dir().join("gatk_rs_hc_default_out.vcf");

    let mut config = GatkConfig::new("HaplotypeCaller".to_string());
    config.set_reference(ref_fa.to_string_lossy().to_string());
    config.add_input_file(bam.to_string_lossy().to_string());
    config.set_output_vcf(out.to_string_lossy().to_string());
    config.set_parameter("intervals".to_string(), "chr1:1-32".to_string());
    config.validate().expect("valid config");

    run_haplotype_caller(&config).expect("default run");
    let text = fs::read_to_string(&out).expect("read out vcf");
    assert!(text.contains("##GATK_RS_HC_PIPELINE=assembly-region-v1"));
    assert!(text.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"));
}
