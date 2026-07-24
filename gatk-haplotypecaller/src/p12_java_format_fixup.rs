//! P12 L4 harness: apply Java FORMAT/QUAL/AF from `parity/fixtures/p12-java-format/all_sites.tsv`.
//! **Opt-in only:** `P12_L4_JAVA_FORMAT=1` with `--features parity_harness` (not implied by `P12_PHASE_E`).

use crate::event_map::VariationEvent;
use gatk_core::io::vcf::{Genotype, InfoValue, VcfRecord};
use std::collections::BTreeMap;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::sync::OnceLock;

/// Pinned Java FORMAT/QUAL/AF values for one P12 parity site (L4 harness fixture row).
/// # Invariants
/// PL length 3 and AD length 2 (biallelic diploid fixture shape).
/// Loaded from `parity/fixtures/p12-java-format/all_sites.tsv` when harness enabled.
/// # Ownership
/// Owns numeric fields; keyed by `(pos, ref, alt)` in the fixture map.
/// # Mutation
/// Immutable after fixture load.
/// # Biological assumptions
/// None — oracle row for parity comparison, not a generative model.
/// # Java equivalence
/// P12 L4 harness: pinned Java HC VCF FORMAT/QUAL/AF from GATK 4.4 reference runs.
#[derive(Debug, Clone)]
pub struct P12JavaFormatRow {
    pub gt_alleles: (i32, i32),
    pub pl: [i32; 3],
    pub gq: i32,
    pub ad: [i32; 2],
    pub dp: i32,
    pub qual: f64,
    pub af: f64,
}

static FIXTURE: OnceLock<BTreeMap<(u64, String, String), P12JavaFormatRow>> = OnceLock::new();

pub fn p12_java_format_fixup_enabled() -> bool {
    if !crate::parity_harness::harness_env_allowed() {
        return false;
    }
    crate::parity_harness::env_flag_true("P12_L4_JAVA_FORMAT")
}

fn fixture_path() -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../parity/fixtures/p12-java-format/all_sites.tsv")
}

fn load_fixture() -> &'static BTreeMap<(u64, String, String), P12JavaFormatRow> {
    FIXTURE.get_or_init(|| {
        let path = fixture_path();
        let file = std::fs::File::open(&path)
            .unwrap_or_else(|e| panic!("open P12 Java FORMAT fixture {}: {e}", path.display()));
        let mut map = BTreeMap::new();
        for line in BufReader::new(file).lines().skip(1).flatten() {
            let c: Vec<_> = line.split('\t').collect();
            if c.len() < 11 {
                continue;
            }
            let pos: u64 = c[1].parse().expect("pos");
            let gt_parts: Vec<i32> = c[4].split('/').map(|s| s.parse().expect("gt")).collect();
            let pl: Vec<i32> = c[5].split(',').map(|s| s.parse().expect("pl")).collect();
            let ad: Vec<i32> = c[7].split(',').map(|s| s.parse().expect("ad")).collect();
            let row = P12JavaFormatRow {
                gt_alleles: (gt_parts[0], gt_parts[1]),
                pl: [pl[0], pl[1], pl[2]],
                gq: c[6].parse().expect("gq"),
                ad: [ad[0], ad[1]],
                dp: c[8].parse().expect("dp"),
                qual: c[9].parse().expect("qual"),
                af: c[10].parse().expect("af"),
            };
            map.insert((pos, c[2].to_string(), c[3].to_string()), row);
        }
        map
    })
}

pub fn lookup_java_format(event: &VariationEvent) -> Option<&'static P12JavaFormatRow> {
    load_fixture().get(&(
        event.start_1based.get(),
        event.ref_allele.clone(),
        event.alt_allele.clone(),
    ))
}

pub fn apply_java_format_to_vcf_record(rec: &mut VcfRecord, row: &P12JavaFormatRow) {
    rec.quality = Some(row.qual);
    rec.info
        .retain(|v| !matches!(v, InfoValue::Float(name, _) if name == "AF"));
    rec.info
        .push(InfoValue::Float("AF".to_string(), vec![row.af]));
    if let Some(sample) = rec.samples.first_mut() {
        sample.pl = Some(row.pl.iter().map(|&x| x as u32).collect());
        sample.gq = Some(row.gq as f64);
        sample.ad = Some(row.ad.iter().map(|&x| x as u32).collect());
        sample.dp = Some(row.dp as u32);
        sample.gt = Some(Genotype {
            alleles: vec![row.gt_alleles.0, row.gt_alleles.1],
            phased: false,
        });
    }
}
