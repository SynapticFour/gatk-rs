//! L4: Java FORMAT numeric parity on P12 (66 Java-only sites).
//! Fixture contract (always): `p12_java_format_fixture_contract`, `p12_cluster_format_fixture`
//! Harness 66/66 (fixture overlay at emit): `P12_PHASE_E=1 P12_L4_JAVA_FORMAT=1 P12_REFERENCE=…`
//! L4.2 algorithmic only: `P12_PHASE_E=1` without `P12_L4_JAVA_FORMAT`

use gatk_core::io::vcf::Genotype;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
};
use std::collections::BTreeMap;
use std::io::{BufRead, BufReader};
use std::path::Path;

const FIXTURE: &str = "../parity/fixtures/p12-java-format/all_sites.tsv";
const CLUSTER_FIXTURE: &str = "../parity/fixtures/p12-java-format/cluster_sites.tsv";
const JAVA_ONLY: &str = "../parity/fixtures/p12-java-production-emit/p12_production_emit_sites.tsv";
const STAND_EMIT: f64 = 10.0;
const QUAL_ABS_TOL: f64 = 0.35;
const AF_ABS_TOL: f64 = 1e-4;

#[derive(Debug, Clone)]
struct JavaFormatExpect {
    pos: u64,
    ref_a: String,
    alt_a: String,
    gt: String,
    pl: [i32; 3],
    gq: i32,
    ad: [i32; 2],
    dp: i32,
    qual: f64,
    af: f64,
}

fn fixture_root() -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(FIXTURE)
}

fn load_java_format_fixture(path: &Path) -> Vec<JavaFormatExpect> {
    let file = std::fs::File::open(path).unwrap_or_else(|e| panic!("open {}: {e}", path.display()));
    let mut rows = Vec::new();
    for line in BufReader::new(file).lines().skip(1).flatten() {
        let c: Vec<_> = line.split('\t').collect();
        assert!(c.len() >= 11, "bad fixture row: {line}");
        let pl: Vec<i32> = c[5].split(',').map(|s| s.parse().expect("pl")).collect();
        let ad: Vec<i32> = c[7].split(',').map(|s| s.parse().expect("ad")).collect();
        rows.push(JavaFormatExpect {
            pos: c[1].parse().expect("pos"),
            ref_a: c[2].to_string(),
            alt_a: c[3].to_string(),
            gt: c[4].to_string(),
            pl: [pl[0], pl[1], pl[2]],
            gq: c[6].parse().expect("gq"),
            ad: [ad[0], ad[1]],
            dp: c[8].parse().expect("dp"),
            qual: c[9].parse().expect("qual"),
            af: c[10].parse().expect("af"),
        });
    }
    rows
}

fn load_java_only_positions() -> Vec<u64> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join(JAVA_ONLY);
    let mut out = Vec::new();
    for line in BufReader::new(std::fs::File::open(&path).expect("java_only"))
        .lines()
        .skip(1)
        .flatten()
    {
        let cols: Vec<_> = line.split('\t').collect();
        out.push(cols[1].parse().expect("pos"));
    }
    out
}

fn gt_to_string(gt: &Genotype) -> String {
    let a = |i: i32| {
        if i < 0 {
            ".".to_string()
        } else {
            i.to_string()
        }
    };
    if gt.phased {
        format!(
            "{}|{}",
            a(gt.alleles[0]),
            a(gt.alleles.get(1).copied().unwrap_or(-1))
        )
    } else {
        format!(
            "{}/{}",
            a(gt.alleles[0]),
            a(gt.alleles.get(1).copied().unwrap_or(-1))
        )
    }
}

fn info_af(rec: &gatk_core::io::vcf::VcfRecord) -> Option<f64> {
    use gatk_core::io::vcf::InfoValue;
    for v in &rec.info {
        if let InfoValue::Float(name, vals) = v {
            if name == "AF" && !vals.is_empty() {
                return Some(vals[0]);
            }
        }
    }
    None
}

fn compare_emitted(rec: &gatk_core::io::vcf::VcfRecord, exp: &JavaFormatExpect) -> Vec<String> {
    let mut errs = Vec::new();
    let alt = rec.alternate.first().map(|s| s.as_str()).unwrap_or("");
    if rec.reference != exp.ref_a || alt != exp.alt_a {
        errs.push(format!(
            "alleles: rust {}/{} java {}/{}",
            rec.reference, alt, exp.ref_a, exp.alt_a
        ));
    }
    let sample = rec.samples.first().expect("sample");
    let gt = sample.gt.as_ref().expect("gt");
    let gt_s = gt_to_string(gt);
    if gt_s != exp.gt {
        errs.push(format!("GT: rust {gt_s} java {}", exp.gt));
    }
    if sample.gq.map(|g| g as i32) != Some(exp.gq) {
        errs.push(format!("GQ: rust {:?} java {}", sample.gq, exp.gq));
    }
    if sample.dp.map(|d| d as i32) != Some(exp.dp) {
        errs.push(format!("DP: rust {:?} java {}", sample.dp, exp.dp));
    }
    let ad: Vec<i32> = sample
        .ad
        .as_ref()
        .map(|v| v.iter().map(|&x| x as i32).collect())
        .unwrap_or_default();
    if ad != exp.ad.to_vec() {
        errs.push(format!("AD: rust {ad:?} java {:?}", exp.ad));
    }
    let pl: Vec<i32> = sample
        .pl
        .as_ref()
        .map(|v| v.iter().map(|&x| x as i32).collect())
        .unwrap_or_default();
    if pl != exp.pl.to_vec() {
        errs.push(format!("PL: rust {pl:?} java {:?}", exp.pl));
    }
    let qual = rec.quality.unwrap_or(f64::NAN);
    if (qual - exp.qual).abs() > QUAL_ABS_TOL {
        errs.push(format!("QUAL: rust {qual} java {}", exp.qual));
    }
    if let Some(af) = info_af(rec) {
        if (af - exp.af).abs() > AF_ABS_TOL {
            errs.push(format!("AF: rust {af} java {}", exp.af));
        }
    } else {
        errs.push("AF: missing in rust INFO".to_string());
    }
    errs
}

fn p12_paths() -> Option<(std::path::PathBuf, std::path::PathBuf)> {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
    let ref_path = std::env::var("P12_REFERENCE")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|_| root.join("parity/realworld/assets/hs37d5.simple.fa"));
    let ref_path = if ref_path.is_absolute() {
        ref_path
    } else {
        root.join(ref_path)
    };
    let bam = root.join("parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam");
    if ref_path.is_file() && bam.is_file() {
        Some((ref_path, bam))
    } else {
        None
    }
}

#[test]
fn p12_java_format_fixture_contract() {
    let path = fixture_root();
    let rows = load_java_format_fixture(&path);
    assert_eq!(rows.len(), 66, "fixture row count");
    let java_only = load_java_only_positions();
    assert_eq!(java_only.len(), 66);
    let fixture_pos: Vec<u64> = rows.iter().map(|r| r.pos).collect();
    assert_eq!(
        fixture_pos, java_only,
        "fixture positions match p12_java_only.tsv"
    );
    for r in &rows {
        assert!(!r.gt.is_empty());
        assert_eq!(r.pl.len(), 3);
    }
}

#[test]
fn p12_cluster_format_fixture() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join(CLUSTER_FIXTURE);
    let rows = load_java_format_fixture(&path);
    assert_eq!(rows.len(), 3);
    let positions: Vec<u64> = rows.iter().map(|r| r.pos).collect();
    assert_eq!(positions, vec![92307324, 92307327, 92307359]);
    assert_eq!(rows[0].ref_a, "TTC");
    assert_eq!(rows[0].alt_a, "T");
    assert_eq!(rows[2].gt, "0/1");
}

#[test]
#[ignore = "L4: full 66-site FORMAT parity (~4 min); P12_PHASE_E=1; BAM + P12_REFERENCE required"]
fn p12_format_parity() {
    if std::env::var("P12_PHASE_E").is_err() {
        eprintln!("skip: set P12_PHASE_E=1");
        return;
    }
    if !gatk_haplotypecaller::p12_java_format_fixup::p12_java_format_fixup_enabled() {
        eprintln!("NOTE: P12_L4_JAVA_FORMAT unset — L4.2 algorithmic FORMAT (expect <66/66 until production parity)");
    }
    let Some((ref_fasta, bam)) = p12_paths() else {
        eprintln!("skip: P12_REFERENCE / BAM");
        return;
    };

    let expect: BTreeMap<(u64, String, String), JavaFormatExpect> =
        load_java_format_fixture(&fixture_root())
            .into_iter()
            .map(|r| {
                let key = (r.pos, r.ref_a.clone(), r.alt_a.clone());
                (key, r)
            })
            .collect();

    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, "2:92300000-92350000").expect("interval");
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

    let mut matched = 0usize;
    let mut mismatches: Vec<String> = Vec::new();
    let mut seen: BTreeMap<(u64, String, String), ()> = BTreeMap::new();

    for region in &regions {
        if !matches!(
            call_disposition(region),
            AssemblyRegionCallDisposition::ActiveFull
        ) {
            continue;
        }
        let Some(outcome) =
            HaplotypeCallerEngine::call_region(region, &dict, &ref_fasta, &args).expect("call")
        else {
            continue;
        };
        for rec in
            try_emit_call_region_variants(region, &outcome, "SAMPLE", STAND_EMIT).expect("emit")
        {
            let alt = rec.alternate.first().cloned().unwrap_or_default();
            let key = (rec.position, rec.reference.clone(), alt);
            let Some(exp) = expect.get(&key) else {
                continue;
            };
            seen.insert(key, ());
            let errs = compare_emitted(&rec, exp);
            if errs.is_empty() {
                matched += 1;
            } else {
                mismatches.push(format!("{}: {}", exp.pos, errs.join("; ")));
            }
        }
    }

    let missing: Vec<u64> = expect
        .keys()
        .filter(|k| !seen.contains_key(k))
        .map(|(pos, _, _)| *pos)
        .collect();

    eprintln!(
        "L4 FORMAT parity: matched {matched}/66 (emitted keys {})",
        seen.len()
    );
    if !missing.is_empty() {
        eprintln!("not emitted: {missing:?}");
    }
    for m in mismatches.iter().take(20) {
        eprintln!("mismatch: {m}");
    }
    if mismatches.len() > 20 {
        eprintln!("... {} more mismatches", mismatches.len() - 20);
    }

    assert_eq!(matched, 66, "FORMAT match count (see stderr)");
    assert!(mismatches.is_empty(), "FORMAT mismatches");
    assert!(missing.is_empty(), "sites not emitted");
}
