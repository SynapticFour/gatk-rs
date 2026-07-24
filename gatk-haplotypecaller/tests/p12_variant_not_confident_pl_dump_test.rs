//! PL dumps at `variant_not_confident` loci (L2 diagnosis).
//! Run: `P12_REFERENCE=… cargo test -p gatk-haplotypecaller p12_variant_not_confident_pl --release -- --nocapture`

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::VariationEvent;
use gatk_haplotypecaller::genome_loc::GenomePosition;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, format_locus_genotype_pl_dump,
    traverse_assembly_region_walker, AssemblyRegion, AssemblyRegionCallDisposition, CallRegionArgs,
    HaplotypeCallerEngine, HcGenotypingConfig, ReadFilterParams, WalkerTraversalConfig,
};
use std::path::Path;

struct PlDumpCase {
    interval: &'static str,
    region_lo: u64,
    region_hi: u64,
    pos: u64,
    ref_allele: &'static str,
    alt_allele: &'static str,
}

const PL_DUMP_CASES: &[PlDumpCase] = &[
    PlDumpCase {
        interval: "2:92307228-92307400",
        region_lo: 92307324,
        region_hi: 92307364,
        pos: 92307364,
        ref_allele: "T",
        alt_allele: "C",
    },
    PlDumpCase {
        interval: "2:92300000-92350000",
        region_lo: 92316227,
        region_hi: 92316475,
        pos: 92316315,
        ref_allele: "C",
        alt_allele: "G",
    },
    PlDumpCase {
        interval: "2:92300000-92350000",
        region_lo: 92316227,
        region_hi: 92316475,
        pos: 92316328,
        ref_allele: "T",
        alt_allele: "A",
    },
    PlDumpCase {
        interval: "2:92325000-92325350",
        region_lo: 92325071,
        region_hi: 92325332,
        pos: 92325193,
        ref_allele: "C",
        alt_allele: "T",
    },
    PlDumpCase {
        interval: "2:92325000-92325350",
        region_lo: 92325071,
        region_hi: 92325332,
        pos: 92325205,
        ref_allele: "G",
        alt_allele: "A",
    },
];

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
    if !ref_path.is_file() || !bam.is_file() {
        return None;
    }
    Some((ref_path, bam))
}

fn find_region<'a>(regions: &'a [AssemblyRegion], lo: u64, hi: u64) -> &'a AssemblyRegion {
    regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= lo
                && r.end.get() >= hi
        })
        .unwrap_or_else(|| panic!("active region covering {lo}-{hi}"))
}

#[test]
fn p12_variant_not_confident_pl_dump() {
    let Some((ref_fasta, bam)) = p12_paths() else {
        eprintln!("skip: set P12_REFERENCE");
        return;
    };
    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let args = CallRegionArgs::strict_java();
    let config = HcGenotypingConfig::strict_java();
    for case in PL_DUMP_CASES {
        let specs = parse_intervals_cli_string(&dict, case.interval).expect("interval");
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
        let region = find_region(&regions, case.region_lo, case.region_hi);
        let outcome = HaplotypeCallerEngine::call_region(region, &dict, &ref_fasta, &args)
            .expect("call")
            .expect("call_region Some");
        let ref_hap = outcome
            .assembly
            .haplotypes
            .iter()
            .find(|h| h.is_reference)
            .expect("ref hap");
        let pad = ref_hap
            .genome_loc
            .map(|g| g.start_1based())
            .unwrap_or_else(|| outcome.assembly.padded_reference_start_1based());
        let ref_bytes = ref_hap.bases.clone();
        let event = VariationEvent {
            contig: "2".into(),
            start_1based: GenomePosition::new_1based(case.pos),
            end_1based: GenomePosition::new_1based(case.pos),
            ref_allele: case.ref_allele.into(),
            alt_allele: case.alt_allele.into(),
        };
        let dump = format_locus_genotype_pl_dump(
            &event,
            &outcome.read_likelihoods,
            &region.reads,
            &outcome.assembly.haplotypes,
            &ref_bytes,
            pad,
            region.start.get(),
            region.end.get(),
            outcome.assembly.max_mnp_distance(),
            &config,
        )
        .expect("pl dump");
        eprintln!("--- pl_dump {} ---\n{dump}", case.pos);
        assert!(
            dump.contains("GL\t") || dump.contains("reject\t"),
            "PL dump must include GL or reject for {}:\n{dump}",
            case.pos
        );
        if case.pos == 92307364 {
            assert!(
                dump.contains("passes_emit\tfalse") || dump.contains("reject\t"),
                "92307364 PL dump should show not-confident emit (variant_not_confident bucket)"
            );
        }
    }
}
