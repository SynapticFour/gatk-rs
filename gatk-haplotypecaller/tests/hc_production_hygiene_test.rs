//! production surface hygiene (given alleles, sample name).

use gatk_common::GatkConfig;
use gatk_haplotypecaller::engine::CallRegionArgs;
use gatk_haplotypecaller::run_haplotype_caller;
use std::path::Path;

#[test]
fn strict_java_ignores_given_alleles_env() {
    let _guard = EnvGuard::set("GATK_RS_HC_GIVEN_VCF", "parity/fixtures/none.vcf");
    let args = CallRegionArgs::strict_java();
    assert!(
        args.given_alleles.is_empty(),
        "production strict_java must not load GATK_RS_HC_GIVEN_VCF"
    );
}

#[test]
fn sample_name_from_first_bam_rg_sm() {
    let repo = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
    let bam = repo.join("parity/build/sam-indexed-bam/p11_java_positive.bam");
    if !bam.is_file() {
        eprintln!("skip: build p11 bam first (sam-indexed-bam staging)");
        return;
    }
    let mut config = GatkConfig::new("HaplotypeCaller".to_string());
    config.add_input_file(bam.to_string_lossy().into_owned());
    config.set_reference(
        repo.join("parity/fixtures/p5_live_reference.fa")
            .to_string_lossy()
            .into_owned(),
    );
    config.set_output_vcf("/tmp/hc_hygiene_sample.vcf".to_string());
    config.set_parameter("intervals".to_string(), "chrLive:1-63".to_string());
    std::env::set_var("GATK_RS_HC_SCAFFOLD_OUTPUT", "1");
    let result = run_haplotype_caller(&config);
    std::env::remove_var("GATK_RS_HC_SCAFFOLD_OUTPUT");
    result.expect("scaffold HC run");
}

struct EnvGuard {
    key: &'static str,
    prior: Option<String>,
}

impl EnvGuard {
    fn set(key: &'static str, value: &str) -> Self {
        let prior = std::env::var(key).ok();
        std::env::set_var(key, value);
        Self { key, prior }
    }
}

impl Drop for EnvGuard {
    fn drop(&mut self) {
        match &self.prior {
            Some(v) => std::env::set_var(self.key, v),
            None => std::env::remove_var(self.key),
        }
    }
}
