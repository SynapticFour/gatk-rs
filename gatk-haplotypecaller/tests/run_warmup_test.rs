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

#[test]
fn haplotype_caller_chrlive_qd_follows_shared_reservoir_stream() {
    let _guard = hc_env_lock().lock().expect("lock");
    clear_hc_output_env();
    let prev_kbest = std::env::var("GATK_RS_EXPERIMENTAL_KBEST_POLICY").ok();
    let prev_rss = std::env::var("GATK_RS_HC_RSS_ABORT_MIB").ok();
    std::env::remove_var("GATK_RS_EXPERIMENTAL_KBEST_POLICY");
    std::env::remove_var("GATK_RS_HC_RSS_ABORT_MIB");

    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("..");
    let ref_fa = root.join("parity/fixtures/p5_live_reference.fa");
    let bam = root.join("parity/build/sam-indexed-bam/p11_java_positive.bam");
    let out = std::env::temp_dir().join("gatk_rs_hc_chrlive_qd_stream.vcf");

    let mut config = GatkConfig::new("HaplotypeCaller".to_string());
    config.set_reference(ref_fa.to_string_lossy().to_string());
    config.add_input_file(bam.to_string_lossy().to_string());
    config.set_output_vcf(out.to_string_lossy().to_string());
    config.set_parameter("intervals".to_string(), "chrLive:1-63".to_string());
    config.validate().expect("valid config");

    let result = run_haplotype_caller(&config);

    match prev_kbest {
        Some(v) => std::env::set_var("GATK_RS_EXPERIMENTAL_KBEST_POLICY", v),
        None => std::env::remove_var("GATK_RS_EXPERIMENTAL_KBEST_POLICY"),
    }
    match prev_rss {
        Some(v) => std::env::set_var("GATK_RS_HC_RSS_ABORT_MIB", v),
        None => std::env::remove_var("GATK_RS_HC_RSS_ABORT_MIB"),
    }
    result.expect("chrLive HaplotypeCaller");

    let text = fs::read_to_string(&out).expect("read chrLive vcf");
    let line = text
        .lines()
        .find(|ln| ln.starts_with("chrLive\t15\t"))
        .expect("chrLive:15 record");
    let cols: Vec<&str> = line.split('\t').collect();
    assert_eq!(cols[3], "T");
    assert_eq!(cols[4], "A");
    assert_eq!(cols[5], "2224.06");
    let qd: f64 = cols[7]
        .split(';')
        .find_map(|kv| kv.strip_prefix("QD="))
        .expect("QD")
        .parse()
        .expect("qd f64");
    assert_eq!(qd, 26.670352933736748);
    assert_eq!(format!("{qd:.2}"), "26.67");
    let sample = cols[9];
    assert!(sample.contains("2238,151,0"), "{sample}");
    assert!(sample.contains("0,50"), "{sample}");
    let format_keys: Vec<&str> = cols[8].split(':').collect();
    let format_vals: Vec<&str> = sample.split(':').collect();
    let dp = format_keys
        .iter()
        .position(|k| *k == "DP")
        .map(|i| format_vals[i])
        .expect("FORMAT DP");
    let ad = format_keys
        .iter()
        .position(|k| *k == "AD")
        .map(|i| format_vals[i])
        .expect("AD");
    let pl = format_keys
        .iter()
        .position(|k| *k == "PL")
        .map(|i| format_vals[i])
        .expect("PL");
    assert_eq!(dp, "50");
    assert_eq!(ad, "0,50");
    assert_eq!(pl, "2238,151,0");
}
