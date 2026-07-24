//! Compare exact Java vs Rust AF/emit decision values at the five remaining P12 divergent sites.
//! Run: `P12_PHASE_E=1 P12_REFERENCE=… cargo test -p gatk-haplotypecaller p12_java_emit_values --release -- --ignored --nocapture`

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::genome_loc::GenomePosition;
use gatk_haplotypecaller::hc_emit_policy::explain_strict_java_emit_gates;
use gatk_haplotypecaller::hc_genotyping_engine::{java_emit_af_decision, JavaEmitAfDecision};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
};
use std::path::Path;

const STAND: f64 = 10.0;

struct JavaVcfExpect {
    pos: u64,
    ref_a: &'static str,
    alt_a: &'static str,
    qual: f64,
    pl: [i32; 3],
    gq: i32,
    ad: [i32; 2],
}

const FIVE: &[JavaVcfExpect] = &[
    JavaVcfExpect {
        pos: 92316296,
        ref_a: "A",
        alt_a: "T",
        qual: 78.32,
        pl: [90, 6, 0],
        gq: 6,
        ad: [0, 2],
    },
    JavaVcfExpect {
        pos: 92316315,
        ref_a: "C",
        alt_a: "G",
        qual: 78.32,
        pl: [90, 6, 0],
        gq: 6,
        ad: [0, 2],
    },
    JavaVcfExpect {
        pos: 92316328,
        ref_a: "T",
        alt_a: "A",
        qual: 78.32,
        pl: [90, 6, 0],
        gq: 6,
        ad: [0, 2],
    },
    JavaVcfExpect {
        pos: 92325193,
        ref_a: "C",
        alt_a: "T",
        qual: 73.64,
        pl: [81, 0, 36],
        gq: 36,
        ad: [1, 2],
    },
    JavaVcfExpect {
        pos: 92325205,
        ref_a: "G",
        alt_a: "A",
        qual: 73.64,
        pl: [81, 0, 36],
        gq: 36,
        ad: [1, 2],
    },
];

fn pl_to_gl(pl: [i32; 3]) -> [f64; 3] {
    let min = pl[0].min(pl[1]).min(pl[2]);
    [
        (pl[0] - min) as f64 / -10.0,
        (pl[1] - min) as f64 / -10.0,
        (pl[2] - min) as f64 / -10.0,
    ]
}

fn print_decision(label: &str, d: &JavaEmitAfDecision, java: &JavaVcfExpect) {
    eprintln!("--- {label} {} {}/{} ---", java.pos, java.ref_a, java.alt_a);
    eprintln!(
        "java_VCF\tQUAL={}\tPL={:?}\tGQ={}\tAD={:?}",
        java.qual, java.pl, java.gq, java.ad
    );
    eprintln!("rust_raw_GL\t{:?}", d.gl_raw);
    eprintln!("rust_java_GL\t{:?}", d.gl_java_pl_roundtrip);
    eprintln!(
        "AF\tlog10P_no_variant={:.6}\tcall_conf_log10={:.6}\talt_plausible={}\tsite_monomorphic={}",
        d.log10_posterior_no_variant, d.call_conf_log10, d.alt_plausible, d.site_is_monomorphic
    );
    eprintln!(
        "emit\tlog10_vc_conf={:.6}\tphred_scaled={:.2}\tpasses_emit={}",
        d.log10_vc_confidence, d.phred_scaled, d.passes_emit
    );
    eprintln!(
        "delta_QUAL\t{:.2} (rust phred - java QUAL)",
        d.phred_scaled - java.qual
    );
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
#[ignore = "Phase E: five-site Java vs Rust AF value dump"]
fn p12_java_emit_values() {
    if std::env::var("P12_PHASE_E").is_err() {
        eprintln!("skip: set P12_PHASE_E=1");
        return;
    }
    let Some((ref_fasta, bam)) = p12_paths() else {
        return;
    };

    eprintln!("=== Java PL-derived GL reference (invert PL min=0) ===");
    for j in FIVE {
        let gl = pl_to_gl(j.pl);
        let d = java_emit_af_decision(&gl, STAND).expect("java ref gl");
        print_decision("java_PL_inverted", &d, j);
        assert!(
            d.passes_emit,
            "{}: Java VCF PL must pass Rust java_emit_af_decision",
            j.pos
        );
        // HTSJDK PL round-trip: within ~0.35 phred of Java VCF QUAL (see unit tests).
        assert!(
            (d.phred_scaled - j.qual).abs() < 0.35,
            "{}: phred {:.2} vs java QUAL {:.2}",
            j.pos,
            d.phred_scaled,
            j.qual
        );
    }

    eprintln!("\n=== Rust call_region genotypes (actual pipeline) ===");
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
        for java in FIVE {
            let Some(call) = outcome.genotyped_calls.iter().find(|c| {
                c.event.start_1based == GenomePosition::new_1based(java.pos)
                    && c.event.ref_allele == java.ref_a
                    && c.event.alt_allele == java.alt_a
            }) else {
                continue;
            };
            let gl = &call.genotype.genotype_log10_likelihoods;
            let d = java_emit_af_decision(gl, STAND).expect("af");
            print_decision("rust_pipeline", &d, java);
            eprintln!(
                "rust_FMT\tPL={:?}\tGQ={}\tAD={:?}",
                call.genotype.format.pl, call.genotype.format.gq, call.genotype.format.ad
            );
            let gates = explain_strict_java_emit_gates(
                &call.event,
                gl,
                &call.genotype.format,
                STAND,
                true,
                0,
                0,
                &[],
            )
            .expect("gates");
            eprintln!("strict_gates\t{gates}");
            let recs =
                try_emit_call_region_variants(region, &outcome, "SAMPLE", STAND).expect("emit");
            let emitted = recs.iter().any(|r| {
                r.position == java.pos
                    && r.reference == java.ref_a
                    && r.alternate.first().map(|s| s.as_str()) == Some(java.alt_a)
            });
            eprintln!("vcf_emitted\t{emitted}\n");
        }
    }
}
