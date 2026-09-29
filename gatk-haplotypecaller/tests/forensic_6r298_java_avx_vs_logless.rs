//! 6R.298: provenance of the stored AVX allele-T value versus the
//! LOGLESS 24-column matrix. No Rust or Java production change.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r298_java_avx_vs_logless -- --nocapture --test-threads=1
//! ```

use std::path::{Path, PathBuf};
use std::process::Command;

const AVX_TSV: &str = "gatk-haplotypecaller/tests/6r285_java_genotyping.tsv";
const LOGLESS_TSV: &str = "gatk-haplotypecaller/tests/6r297_java_hap_ll.tsv";
const QNAME: &str = "HWI-D00360:5:H814YADXX:2:2207:19511:63503";
const FLAGS: &str = "163";
const AVX_ALLELE_T_BITS: u64 = 0xc0029ba800000000;
const LOGLESS_T_BITS: u64 = 0xc0029ba6bfb37500;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R298\t{key}\t{}", value.as_ref());
}

fn fmt_f(x: f64) -> String {
    format!("{x:.17}")
}

fn cell_bits(cell: &str) -> f64 {
    let hex = cell.split("bits=").nth(1).unwrap().trim();
    f64::from_bits(u64::from_str_radix(hex.trim_start_matches("0x"), 16).unwrap())
}

fn is_exact_f32(x: f64) -> bool {
    (x as f32 as f64).to_bits() == x.to_bits()
}

fn find_row<'a>(text: &'a str, tag: &str, kind: &str) -> Vec<&'a str> {
    for line in text.lines() {
        let p: Vec<_> = line.split('\t').collect();
        if p.len() > 4 && p[0] == tag && p[1] == kind && p[2] == QNAME && p[3] == FLAGS {
            return p;
        }
    }
    panic!("missing {tag} {kind}");
}

fn avx_probe() -> String {
    let java = "/opt/homebrew/opt/openjdk/bin/java";
    let javac = "/opt/homebrew/opt/openjdk/bin/javac";
    let jar = "/tmp/gatk440/gatk-package-4.4.0.0-local.jar";
    let dir = std::env::temp_dir().join("avx-probe-298");
    std::fs::create_dir_all(&dir).unwrap();
    let src = dir.join("AvxProbe298.java");
    std::fs::write(
        &src,
        r#"
import org.broadinstitute.gatk.nativebindings.pairhmm.PairHMMNativeArguments;
import org.broadinstitute.hellbender.utils.pairhmm.PairHMM;
public class AvxProbe298 {
  public static void main(String[] args) {
    PairHMMNativeArguments nativeArgs = new PairHMMNativeArguments();
    nativeArgs.maxNumberOfThreads = 1;
    nativeArgs.useDoublePrecision = false;
    try {
      PairHMM hmm = PairHMM.Implementation.AVX_LOGLESS_CACHING.makeNewHMM(nativeArgs);
      System.out.println("AVX_CONSTRUCTED " + hmm.getClass().getName());
    } catch (Throwable t) {
      System.out.println("AVX_UNAVAILABLE " + t.getClass().getName());
      System.out.println(String.valueOf(t.getMessage()));
    }
  }
}
"#,
    )
    .unwrap();
    let compile = Command::new(javac)
        .args(["-proc:none", "-cp", jar, "-d"])
        .arg(&dir)
        .arg(&src)
        .output()
        .expect("javac");
    assert!(compile.status.success(), "javac failed");
    let run = Command::new(java)
        .args([
            "--add-opens",
            "java.base/java.lang=ALL-UNNAMED",
            "--enable-native-access=ALL-UNNAMED",
            "-cp",
        ])
        .arg(format!("{}:{jar}", dir.display()))
        .arg("AvxProbe298")
        .output()
        .expect("java");
    String::from_utf8_lossy(&run.stdout).into_owned() + &String::from_utf8_lossy(&run.stderr)
}

#[test]
fn forensic_6r298_java_avx_vs_logless() {
    let root = repo_root();
    let avx = std::fs::read_to_string(root.join(AVX_TSV)).unwrap();
    let logless = std::fs::read_to_string(root.join(LOGLESS_TSV)).unwrap();
    let avx_hap = find_row(&avx, "6R285", "hap_ll");
    let avx_al = find_row(&avx, "6R285", "allele_ll");
    let ll_hap = find_row(&logless, "6R297", "hap_ll");
    let ll_al = find_row(&logless, "6R297", "allele_ll");
    assert_eq!(avx_hap.len(), 9, "AVX hap_ll stores five columns, not 0-13");
    assert_eq!(ll_hap[4], "57");
    let avx_t = cell_bits(avx_al[4]);
    let ll_t = cell_bits(ll_al[4]);
    assert_eq!(avx_t.to_bits(), AVX_ALLELE_T_BITS);
    assert_eq!(ll_t.to_bits(), LOGLESS_T_BITS);
    let ll_haps_2_5: Vec<f64> = (2..6).map(|h| cell_bits(ll_hap[5 + h])).collect();
    assert!(ll_haps_2_5.iter().all(|v| v.to_bits() == LOGLESS_T_BITS));
    let avx_floor = cell_bits(avx_hap[4]);
    assert_eq!((avx_floor + 4.5).to_bits(), avx_t.to_bits());
    assert!(avx_hap[4..]
        .iter()
        .map(|c| cell_bits(c))
        .all(|v| v.to_bits() == avx_floor.to_bits()));

    let mut avx_cells = 0usize;
    let mut avx_f32 = 0usize;
    for line in avx.lines() {
        let p: Vec<_> = line.split('\t').collect();
        if p.first() != Some(&"6R285")
            || (p.get(1) != Some(&"hap_ll") && p.get(1) != Some(&"allele_ll"))
        {
            continue;
        }
        for cell in p.iter().skip(4) {
            if !cell.contains("bits=") {
                continue;
            }
            avx_cells += 1;
            if is_exact_f32(cell_bits(cell)) {
                avx_f32 += 1;
            }
        }
    }
    assert!(avx_cells > 0);
    assert_eq!(avx_f32, avx_cells);

    let probe = avx_probe();
    assert!(probe.contains("AVX_UNAVAILABLE"), "probe output: {probe}");
    assert!(probe.contains("Machine does not support AVX PairHMM."));

    kv("producer", "Downstream285.Dump.dumpEvent allele_ll");
    kv("producer_class", "Downstream285");
    kv("producer_method", "dumpEvent");
    kv(
        "producer_lines",
        "Downstream285.java:311-320 via matrix.get; bits() at 163-166",
    );
    kv("hap_ll_method", "dumpFrozenColumns at Downstream285.java:255, called from assignGenotypeLikelihoods line 199 before super");
    kv(
        "backend",
        "PairHMM.Implementation.AVX_LOGLESS_CACHING at Downstream285.java:66",
    );
    kv("fastest_available", "NOT_SELECTED");
    kv("use_double_precision", "false at Downstream285.java:74");
    kv("native_threads", "1 at Downstream285.java:70");
    kv(
        "input_matrix",
        "AlleleLikelihoods filled by PairHMMLikelihoodCalculationEngine.computeReadLikelihoods",
    );
    kv("normalization", "AlleleLikelihoods.normalizeLikelihoods line 417, called from PairHMMLikelihoodCalculationEngine.java:200");
    kv("marginalization", "AlleleLikelihoods.marginalize line 719, called from HaplotypeCallerGenotypingEngine.java:191; per-allele max at AlleleLikelihoods.java:772-774");
    kv(
        "matrix_get_type",
        "LikelihoodMatrix.get(int,int) returns double",
    );
    kv("serialization", "Double.doubleToRawLongBits");
    kv("logless_haps_2_5", fmt_f(ll_t));
    kv("logless_allele_T", fmt_f(ll_t));
    kv("avx_haps_2_5", "NOT_CAPTURED");
    kv("avx_haps_0_13", "NOT_CAPTURED");
    kv("avx_post_norm_19_23", fmt_f(avx_floor));
    kv("avx_allele_T", fmt_f(avx_t));
    kv("stored_allele_delta", fmt_f(ll_t - avx_t));
    kv(
        "first_stored_divergence",
        "post-normalization hap_ll columns 19-23; raw PairHMM output is not in either file",
    );
    kv("avx_cells_exact_f32", format!("{avx_f32}/{avx_cells}"));
    kv("same_input_matrix_comparison", "NOT_RUN");
    kv("max_avx_vs_logless_delta", "NOT_AVAILABLE");
    kv("avx_on_this_host", "UNAVAILABLE x86_64 libgkl on arm64");
    kv("classification", "JAVA_AVX_CAPTURE_NOT_REPRODUCIBLE");
    kv("production_change", "NONE");
}
