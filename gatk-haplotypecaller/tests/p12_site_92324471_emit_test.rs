//! L4: emit probe for P12 site 92324471 (graph-only gap C/T, Java PL 45,3,0 AD 0,1).

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::genome_loc::GenomePosition;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
};
use std::path::Path;

const POS: u64 = 92324471;

#[test]
#[ignore = "P12 BAM: 92324471 emit"]
fn p12_site_92324471_emit() {
    if std::env::var("P12_PHASE_E").is_err() {
        return;
    }
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
    let dict = SequenceDictionary::from_fasta_path(&ref_path).expect("dict");
    let specs = parse_intervals_cli_string(&dict, "2:92300000-92350000").expect("interval");
    let walk = traverse_assembly_region_walker(
        &dict,
        &specs,
        &ref_path,
        &bam,
        &ReadFilterParams::gatk_standard_hc(),
        &WalkerTraversalConfig::gatk_haplotype_caller_production(100),
    )
    .expect("walk");
    let args = CallRegionArgs::strict_java();
    let mut genotyped = false;
    let mut emitted = false;
    for region in flatten_assembly_regions(&walk) {
        if !matches!(
            call_disposition(&region),
            AssemblyRegionCallDisposition::ActiveFull
        ) {
            continue;
        }
        if region.start.get() > POS || region.end.get() < POS {
            continue;
        }
        let Some(outcome) =
            HaplotypeCallerEngine::call_region(&region, &dict, &ref_path, &args).expect("call")
        else {
            continue;
        };
        if outcome
            .genotyped_calls
            .iter()
            .any(|c| c.event.start_1based == GenomePosition::new_1based(POS))
        {
            genotyped = true;
            for c in &outcome.genotyped_calls {
                if c.event.start_1based == GenomePosition::new_1based(POS) {
                    eprintln!(
                        "genotyped PL={:?} AD={:?}",
                        c.genotype.format.pl, c.genotype.format.ad
                    );
                }
            }
            let recs =
                try_emit_call_region_variants(&region, &outcome, "SAMPLE", 10.0).expect("emit");
            if let Some(rec) = recs.iter().find(|r| r.position == POS) {
                emitted = true;
                let s = rec.samples.first().unwrap();
                eprintln!("VCF PL={:?} AD={:?} DP={:?}", s.pl, s.ad, s.dp);
            }
        }
    }
    eprintln!("genotyped={genotyped} emitted={emitted}");
    assert!(genotyped, "92324471 must be genotyped");
    assert!(emitted, "92324471 must be emitted");
}
