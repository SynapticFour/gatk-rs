//! 7R.22: unset `GATK_RS_EXPERIMENTAL_KBEST_POLICY` is the Java-equivalent unbounded
//! frontier policy. Canonical production region is not end-equalized.
//! `legacy_1024` stays selectable. HOLDOUT_6R243 is not this test.

use gatk_core::io::vcf::VcfRecord;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::seq_kbest_resource_policy::SeqKbestResourcePolicy;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegion, AssemblyRegionCallDisposition, CallRegionArgs,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const JAVA_VCF_REL: &str = "parity/reports/6r43/chr20_tiny/java.vcf";
const OLD_RUST_VCF_REL: &str = "parity/reports/6r43/chr20_tiny/rust.vcf";
const JAVA_LL_REL: &str = "parity/reports/7r13/java_hap_ll.tsv";
const TARGET_SCORE: f64 = -3.0279336364708747;
const TARGET_BASES: &str = "AATTATTATAATAGCTGTAAAATAGTTATATGTGTTTCAAATTTCATACTCCAAAAAGTTTGGATATAGATGGAAGATACATAACTATATGTATAGTGAGATAGAGAGCTAGCAGAATAGATAAAAAATAATTTTATAATAAAAGAATACTATATCTGCTTTTTTAATTGGAAGTATTAAATTTATCTATGTGACTGCATTATTTCCTAAAAGGTCACATTTTAAAAAAATTTTAAAAAACATTAGAAGATTGAGTCATTATAATTTATTTCTTAAATTTTAATTATATGTATATCTATGAATAATTTATG";

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("7R22\t{key}\t{}", value.as_ref());
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct SiteKey {
    chrom: String,
    pos: u64,
    ref_a: String,
    alt: String,
}

impl SiteKey {
    fn label(&self) -> String {
        format!("{}:{} {}/{}", self.chrom, self.pos, self.ref_a, self.alt)
    }
}

fn parse_vcf_keys(path: &Path) -> BTreeMap<SiteKey, String> {
    let mut out = BTreeMap::new();
    let text = fs::read_to_string(path).unwrap_or_default();
    for line in text.lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let c: Vec<_> = line.split('\t').collect();
        if c.len() < 5 {
            continue;
        }
        let key = SiteKey {
            chrom: c[0].to_string(),
            pos: c[1].parse().unwrap_or(0),
            ref_a: c[3].to_string(),
            alt: c[4].to_string(),
        };
        out.insert(key, line.to_string());
    }
    out
}

fn record_key(rec: &VcfRecord) -> SiteKey {
    SiteKey {
        chrom: rec.chromosome.clone(),
        pos: rec.position,
        ref_a: rec.reference.clone(),
        alt: rec.alternate.join(","),
    }
}

fn covers_site(region: &AssemblyRegion) -> bool {
    region.start.get() <= 29_455_902 && region.end.get() >= 29_455_902
}

#[test]
fn default_unbounded_canonical_hc_matches_established_java_event() {
    let policy = SeqKbestResourcePolicy::from_runtime().expect("policy");
    assert_eq!(policy, SeqKbestResourcePolicy::UnboundedDiagnostic);
    assert_eq!(policy.name(), "unbounded_diagnostic");
    assert!(
        std::env::var_os("GATK_RS_EXPERIMENTAL_KBEST_POLICY").is_none(),
        "this proof requires the production default, not an explicit policy env"
    );
    assert!(
        gatk_haplotypecaller::runtime_config::hc_rss_abort_mib().is_none(),
        "RSS abort must stay unset"
    );

    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    let java_vcf = root.join(JAVA_VCF_REL);
    if !ref_fasta.is_file() || !bam.is_file() || !java_vcf.is_file() {
        eprintln!("skip: canonical BAM, reference, or Java VCF absent");
        return;
    }

    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, INTERVAL).expect("interval");
    let walk = traverse_assembly_region_walker(
        &dict,
        &specs,
        &ref_fasta,
        &bam,
        &ReadFilterParams::gatk_standard_hc(),
        &WalkerTraversalConfig::gatk_haplotype_caller_production(100),
    )
    .expect("walk");
    let regions = flatten_assembly_regions(&walk);
    let args = CallRegionArgs::strict_java();
    let java = parse_vcf_keys(&java_vcf);
    let old_rust = parse_vcf_keys(&root.join(OLD_RUST_VCF_REL));
    let mut live: BTreeMap<SiteKey, VcfRecord> = BTreeMap::new();
    let mut saw_site = false;

    for region in regions.iter().filter(|r| {
        matches!(
            call_disposition(r),
            AssemblyRegionCallDisposition::ActiveFull
        )
    }) {
        kv(
            "active_region",
            format!(
                "{}:{}-{} reads={}",
                region.contig,
                region.start.get(),
                region.end.get(),
                region.reads.len()
            ),
        );
        let outcome = HaplotypeCallerEngine::call_region(region, &dict, &ref_fasta, &args)
            .expect("call")
            .expect("active region returns Some");
        if covers_site(region) {
            saw_site = true;
            assert_eq!(region.end.get(), 29_455_993);
            let haps = &outcome.assembly.haplotypes;
            assert_eq!(haps.len(), 76, "trimmed haplotype count");
            let target = &haps[5];
            assert!(!target.is_reference);
            assert_eq!(target.bases.len(), 311);
            assert_eq!(target.score.to_bits(), TARGET_SCORE.to_bits());
            assert_eq!(std::str::from_utf8(&target.bases).unwrap(), TARGET_BASES);
            let ref_i = haps.iter().position(|h| h.is_reference).expect("reference");
            assert_eq!(ref_i, 59);
            assert_eq!(haps[ref_i].bases.len(), 314);
            let mut matrix_reads = BTreeSet::new();
            for cell in &outcome.read_likelihoods {
                matrix_reads.insert(cell.read_index.get());
            }
            kv(
                "matrix",
                format!(
                    "genotyping_reads={} matrix_reads={} cells={}",
                    outcome.genotyping_reads.len(),
                    matrix_reads.len(),
                    outcome.read_likelihoods.len()
                ),
            );
            assert_eq!(matrix_reads.len(), 332);
            assert_eq!(outcome.read_likelihoods.len(), 76 * 332);
            compare_java_matrix(&root.join(JAVA_LL_REL), &outcome);
            compare_java_haps(&root.join(JAVA_LL_REL), haps);
        }
        let emitted = try_emit_call_region_variants(
            region,
            &outcome,
            "SAMPLE",
            DEFAULT_STAND_EMIT_CONFIDENCE,
        )
        .expect("emit");
        for rec in emitted {
            live.insert(record_key(&rec), rec);
        }
    }
    assert!(saw_site, "walker region covering 29455902");

    let target = SiteKey {
        chrom: "20".into(),
        pos: 29_455_902,
        ref_a: "G".into(),
        alt: "A".into(),
    };
    let rec = live.get(&target).expect("20:29455902 G/A emitted");
    let sample = rec.samples.first().expect("sample");
    assert_eq!(
        sample.gt.as_ref().map(|g| g.alleles.as_slice()),
        Some(&[0, 1][..])
    );
    assert_eq!(sample.pl.as_deref(), Some(&[203, 0, 3081][..]));
    assert_eq!(sample.ad.as_deref(), Some(&[95, 15][..]));
    assert_eq!(sample.dp, Some(110));
    let qual = rec.quality.expect("qual");
    kv("canonical_qual", format!("{qual:.12}"));
    assert_eq!(format!("{qual:.12}"), "195.639586842884");
    let java_line = java.get(&target).expect("Java has 20:29455902 G/A");
    assert!(
        java_line.split('\t').nth(5) == Some("195.64"),
        "Java QUAL print: {java_line}"
    );

    let live_keys: BTreeSet<_> = live.keys().cloned().collect();
    let java_keys: BTreeSet<_> = java.keys().cloned().collect();
    let old_keys: BTreeSet<_> = old_rust.keys().cloned().collect();
    let java_only: Vec<_> = java_keys.difference(&live_keys).cloned().collect();
    let rust_only: Vec<_> = live_keys.difference(&java_keys).cloned().collect();
    let old_rust_only: BTreeSet<_> = old_keys.difference(&java_keys).cloned().collect();
    let previously_common: BTreeSet<_> = java_keys.intersection(&old_keys).cloned().collect();
    let mut format_diff = 0usize;
    let mut first_format: Option<String> = None;
    for key in java_keys.intersection(&live_keys) {
        let Some(line) = java.get(key) else {
            continue;
        };
        let Some(rec) = live.get(key) else {
            continue;
        };
        let cols: Vec<_> = line.split('\t').collect();
        let sample = rec.samples.first();
        let gt = sample
            .and_then(|s| s.gt.as_ref())
            .map(|g| {
                if g.alleles.len() == 2 {
                    format!("{}/{}", g.alleles[0], g.alleles[1])
                } else {
                    format!("{:?}", g.alleles)
                }
            })
            .unwrap_or_default();
        let pl = sample
            .and_then(|s| s.pl.as_ref())
            .map(|v| {
                v.iter()
                    .map(|x| x.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            })
            .unwrap_or_default();
        let ad = sample
            .and_then(|s| s.ad.as_ref())
            .map(|v| {
                v.iter()
                    .map(|x| x.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            })
            .unwrap_or_default();
        let dp = sample
            .and_then(|s| s.dp)
            .map(|v| v.to_string())
            .unwrap_or_default();
        let qual = rec
            .quality
            .map(|q| format!("{q:.2}"))
            .unwrap_or_else(|| ".".into());
        let j_qual = cols.get(5).copied().unwrap_or("");
        let fmt = cols.get(8).copied().unwrap_or("");
        let mut jgt = "";
        let mut jad = "";
        let mut jdp = "";
        let mut jpl = "";
        if let Some(sample_col) = cols.get(9) {
            let keys: Vec<_> = fmt.split(':').collect();
            let vals: Vec<_> = sample_col.split(':').collect();
            for (k, v) in keys.iter().zip(vals.iter()) {
                match *k {
                    "GT" => jgt = v,
                    "AD" => jad = v,
                    "DP" => jdp = v,
                    "PL" => jpl = v,
                    _ => {}
                }
            }
        }
        if gt != jgt || pl != jpl || ad != jad || dp != jdp || qual != j_qual {
            format_diff += 1;
            if first_format.is_none() {
                first_format = Some(format!(
                    "{} GT {jgt}/{gt} PL {jpl}/{pl} AD {jad}/{ad} DP {jdp}/{dp} QUAL {j_qual}/{qual}",
                    key.label()
                ));
            }
        }
    }
    kv("common_format_or_qual_diff", format_diff.to_string());
    if let Some(first) = &first_format {
        kv("first_format_or_qual_diff", first);
    }
    kv("java_records", java_keys.len().to_string());
    kv("live_records", live_keys.len().to_string());
    kv("java_only", java_only.len().to_string());
    for k in &java_only {
        kv("java_only_site", k.label());
    }
    kv("rust_only", rust_only.len().to_string());
    for k in &rust_only {
        kv("rust_only_site", k.label());
    }
    let dropped: Vec<_> = previously_common.difference(&live_keys).cloned().collect();
    for k in &dropped {
        kv("dropped_previously_common", k.label());
    }
    assert!(
        dropped.is_empty(),
        "previously common Java sites missing after the default switch: {dropped:?}"
    );
    let new_rust_only: Vec<_> = rust_only
        .iter()
        .filter(|k| !old_rust_only.contains(*k))
        .cloned()
        .collect();
    for k in &new_rust_only {
        kv("new_rust_only", k.label());
    }
    assert!(
        new_rust_only.is_empty(),
        "new Rust-only calls after unbounded default: {new_rust_only:?}"
    );
    assert!(!java_only
        .iter()
        .any(|k| k.pos == 29_455_902 && k.alt == "A"));
}

fn compare_java_haps(path: &Path, haps: &[gatk_haplotypecaller::Haplotype]) {
    if !path.is_file() {
        return;
    }
    let text = fs::read_to_string(path).expect("java hap tsv");
    let mut rows = Vec::new();
    for line in text.lines() {
        let c: Vec<_> = line.split('\t').collect();
        if c.len() >= 7 && c[1] == "hap" {
            rows.push(c[6].to_string());
        }
    }
    assert_eq!(rows.len(), haps.len());
    for (i, (java, rust)) in rows.iter().zip(haps.iter()).enumerate() {
        let bases = std::str::from_utf8(&rust.bases).unwrap();
        assert_eq!(bases, java, "haplotype index {i}");
    }
    kv("haplotype_inventory", "76/76 identical index order");
}

fn compare_java_matrix(path: &Path, outcome: &gatk_haplotypecaller::CallRegionOutcome) {
    if !path.is_file() {
        return;
    }
    let text = fs::read_to_string(path).expect("java ll");
    let mut java: BTreeMap<(String, u16, usize), f64> = BTreeMap::new();
    let mut n_rows = 0usize;
    for line in text.lines() {
        let c: Vec<_> = line.split('\t').collect();
        if c.len() < 16 || c[1] != "hap_ll" {
            continue;
        }
        n_rows += 1;
        let qname = c[2].to_string();
        let flag: u16 = c[3].parse().expect("flag");
        let mut hap = 0usize;
        let mut i = 15usize;
        while i < c.len() {
            let v: f64 = c[i].parse().expect("ll");
            java.insert((qname.clone(), flag, hap), v);
            hap += 1;
            i += 2;
        }
        assert_eq!(hap, 76);
    }
    assert_eq!(n_rows, 332);
    let mut mismatches = 0usize;
    for cell in &outcome.read_likelihoods {
        let read = &outcome.genotyping_reads[cell.read_index.get()];
        let qname = String::from_utf8_lossy(read.qname()).into_owned();
        let flag = read.flags();
        let key = (qname, flag, cell.haplotype_index.get());
        let expected = java
            .get(&key)
            .unwrap_or_else(|| panic!("missing java cell {key:?}"));
        if expected.to_bits() != cell.log10_likelihood.to_bits()
            && *expected != cell.log10_likelihood
        {
            mismatches += 1;
            if mismatches == 1 {
                kv(
                    "first_ll_mismatch",
                    format!(
                        "{key:?} java={expected:.16} rust={:.16}",
                        cell.log10_likelihood
                    ),
                );
            }
        }
    }
    assert_eq!(mismatches, 0);
    assert_eq!(java.len(), 76 * 332);
    kv("pairhmm", "76x332 all cells match parsed f64");
}
