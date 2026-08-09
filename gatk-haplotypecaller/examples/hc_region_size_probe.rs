//! Memory-safe probe: stream assembly regions on a small GIAB window and report
//! max reference / haplotype / read lengths **without** running PairHMM.
//!
//! ```text
//! cargo run -p gatk-haplotypecaller --release --example hc_region_size_probe -- \
//!   parity/realworld/assets/hs37d5.simple.fa \
//!   parity/realworld/na12878_giab_window_mem_2mb_b37/NA12878_giab_window.b37.bam \
//!   20:10000000-10050000
//! ```

use gatk_core::reference::{parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary};
use gatk_haplotypecaller::{
    assemble_reads, for_each_assembly_region, AssembleReadsArgs, ReadFilterParams,
    WalkerTraversalConfig, GATK_DEFAULT_ASSEMBLY_REGION_PADDING,
};
use std::env;
use std::path::Path;
use std::process;

const ABORT_LEN: usize = 50_000;

fn main() {
    let mut args = env::args().skip(1);
    let Some(ref_fa) = args.next() else {
        eprintln!("usage: hc_region_size_probe <ref.fa> <bam> <interval>");
        process::exit(2);
    };
    let Some(bam) = args.next() else {
        eprintln!("usage: hc_region_size_probe <ref.fa> <bam> <interval>");
        process::exit(2);
    };
    let Some(interval) = args.next() else {
        eprintln!("usage: hc_region_size_probe <ref.fa> <bam> <interval>");
        process::exit(2);
    };

    let dict = SequenceDictionary::from_fasta_path(&ref_fa).expect("dict");
    let specs = parse_intervals_cli_string(&dict, &interval).expect("interval");
    let cfg = WalkerTraversalConfig::gatk_haplotype_caller_production(
        GATK_DEFAULT_ASSEMBLY_REGION_PADDING,
    );
    let filters = ReadFilterParams::default();
    let mut ref_cache = ReferenceWindowCache::new(Path::new(&ref_fa).to_path_buf(), 4);
    let assemble_args = AssembleReadsArgs::default();

    let mut n_regions = 0usize;
    let mut n_active = 0usize;
    let mut max_ref = 0usize;
    let mut max_ext = 0usize;
    let mut max_reads = 0usize;
    let mut max_read_seq = 0usize;
    let mut max_hap = 0usize;

    for_each_assembly_region(
        &dict,
        &specs,
        Path::new(&ref_fa),
        Path::new(&bam),
        &filters,
        &cfg,
        |_idx, region| {
            n_regions += 1;
            let ext_len = region
                .extended_end
                .get()
                .saturating_sub(region.extended_start.get())
                .saturating_add(1) as usize;
            max_ext = max_ext.max(ext_len);
            max_ref = max_ref.max(region.reference.bases.len());
            max_reads = max_reads.max(region.reads.len());
            for r in &region.reads {
                max_read_seq = max_read_seq.max(r.seq().len() as usize);
            }
            if region.reference.bases.len() >= ABORT_LEN || ext_len >= ABORT_LEN {
                eprintln!(
                    "ABORT contig-scale region {} {}:{}-{} ref_bases={} ext_len={}",
                    region.contig,
                    region.contig,
                    region.start.get(),
                    region.end.get(),
                    region.reference.bases.len(),
                    ext_len
                );
                process::exit(3);
            }
            if region.is_active {
                n_active += 1;
                eprintln!(
                    "assemble_begin active={n_active} {}:{}-{} ext={}-{} reads={} region_ref={}",
                    region.contig,
                    region.start.get(),
                    region.end.get(),
                    region.extended_start.get(),
                    region.extended_end.get(),
                    region.reads.len(),
                    region.reference.bases.len()
                );
                let set = assemble_reads(region, &dict, &mut ref_cache, &assemble_args)?;
                let mut region_max_hap = 0usize;
                for h in &set.haplotypes {
                    region_max_hap = region_max_hap.max(h.bases.len());
                    max_hap = max_hap.max(h.bases.len());
                    if h.bases.len() >= ABORT_LEN {
                        eprintln!(
                            "ABORT contig-scale haplotype {}:{}-{} hap_len={} is_ref={}",
                            region.contig,
                            region.start.get(),
                            region.end.get(),
                            h.bases.len(),
                            h.is_reference
                        );
                        process::exit(4);
                    }
                }
                max_ref = max_ref.max(set.reference_bases().len());
                if set.reference_bases().len() >= ABORT_LEN {
                    eprintln!(
                        "ABORT contig-scale assembly ref {}:{}-{} len={}",
                        region.contig,
                        region.start.get(),
                        region.end.get(),
                        set.reference_bases().len()
                    );
                    process::exit(5);
                }
                eprintln!(
                    "assemble_end active={n_active} {}:{}-{} haps={} max_hap={region_max_hap} asm_ref={}",
                    region.contig,
                    region.start.get(),
                    region.end.get(),
                    set.haplotypes.len(),
                    set.reference_bases().len()
                );
                drop(set);
            }
            Ok(())
        },
    )
    .expect("traversal");

    println!("interval={interval}");
    println!("regions={n_regions}");
    println!("active={n_active}");
    println!("max_extended_span_bp={max_ext}");
    println!("max_region_reference_bases={max_ref}");
    println!("max_reads_per_region={max_reads}");
    println!("max_read_seq_len={max_read_seq}");
    println!("max_haplotype_len={max_hap}");
}
