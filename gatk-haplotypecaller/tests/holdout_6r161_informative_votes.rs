//! 6R.161: per-read informative allele votes at `2:92316347 G/A`.
//! Production `call_region` (no diagnostic k-best). Skipped unless `HOLDOUT_6R161=1`.
//!
//! Does not patch production. Does not investigate the later `45,3,0` overwrite.
//!
//! ```text
//! HOLDOUT_6R161=1 cargo test -p gatk-haplotypecaller --test holdout_6r161_informative_votes -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::{build_per_haplotype_variation_events, VariationEvent};
use gatk_haplotypecaller::hc_allele_mapping::create_allele_mapper_with_events;
use gatk_haplotypecaller::hc_genotyping_engine::java_alignment_read_overlaps_interval;
use gatk_haplotypecaller::query_index_at_reference_position;
use gatk_haplotypecaller::read_realignment::LOG_10_INFORMATIVE_THRESHOLD;
use gatk_haplotypecaller::{
    begin_hap_list_observe, begin_likelihood_pipeline_observe, begin_poorly_modeled_observe,
    call_disposition, flatten_assembly_regions, marginalize_rows_to_biallelic_alleles,
    region_likelihoods_to_rows, take_hap_list_snaps, take_likelihood_pipeline_snaps,
    take_poorly_modeled_observe, traverse_assembly_region_walker, AssemblyRegionCallDisposition,
    CallRegionArgs, GenomePosition, HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
    DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN,
};
use rust_htslib::bam::record::CigarString;
use rust_htslib::bam::Record;
use std::path::{Path, PathBuf};

const INTERVAL: &str = "2:92316200-92316580";
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 92_316_347;
const MERGED_REF: &str = "G";
const MERGED_ALT: &str = "A";
const JAVA_INFORMATIVE: f64 = 0.2;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R161\t{key}\t{}", value.as_ref());
}

fn fnv1a64_hex(bases: &[u8]) -> String {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in bases {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x100000001b3);
    }
    format!("{h:016x}")
}

fn snp_base_at(rec: &Record) -> Option<u8> {
    let loc0 = (TARGET as i64).saturating_sub(1);
    if rec.is_unmapped() {
        return None;
    }
    let cigar = CigarString(rec.cigar().iter().copied().collect());
    let qi = query_index_at_reference_position(rec.pos(), &cigar, loc0)?;
    let seq = rec.seq();
    if qi >= seq.len() {
        return None;
    }
    Some(seq.as_bytes()[qi].to_ascii_uppercase())
}

/// Java `searchBestAllele` + REF priority 1.0 + `isInformative` (`confidence > 0.2`).
fn java_best_breaking_ties(lls: &[f64], ref_i: usize) -> (usize, usize, f64, bool) {
    let n = lls.len();
    if n == 0 {
        return (0, 0, 0.0, false);
    }
    let mut best_i = 0usize;
    let mut second_i = 0usize;
    let mut best = lls[0];
    let mut second = f64::NEG_INFINITY;
    for a in 1..n {
        let cand = lls[a];
        if cand > best {
            second_i = best_i;
            second = best;
            best_i = a;
            best = cand;
        } else if cand > second {
            second_i = a;
            second = cand;
        }
    }
    let priorities: Vec<f64> = (0..n).map(|i| if i == ref_i { 1.0 } else { 0.0 }).collect();
    if best - second < JAVA_INFORMATIVE {
        let mut best_pri = priorities[best_i];
        let mut second_pri = priorities[second_i];
        for a in 0..n {
            let cand = lls[a];
            if a == best_i || best - cand > JAVA_INFORMATIVE {
                continue;
            }
            let pri = priorities[a];
            if pri > best_pri {
                second_i = best_i;
                best_i = a;
                second_pri = best_pri;
                best_pri = pri;
            } else if pri > second_pri {
                second_i = a;
                second_pri = pri;
            }
        }
    }
    let best_ll = lls[best_i];
    let second_ll = if second_i != best_i {
        lls[second_i]
    } else {
        f64::NEG_INFINITY
    };
    let conf = if best_ll == second_ll {
        0.0
    } else {
        best_ll - second_ll
    };
    (best_i, second_i, conf, conf > JAVA_INFORMATIVE)
}

fn allele_name(i: usize) -> &'static str {
    match i {
        0 => "REF",
        1 => "ALT",
        _ => "OTHER",
    }
}

#[test]
fn holdout_6r161_informative_per_read_votes() {
    if std::env::var("HOLDOUT_6R161").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R161=1");
        return;
    }
    assert_ne!(
        std::env::var("GATK_RS_EXPERIMENTAL_KBEST_POLICY")
            .ok()
            .as_deref(),
        Some("unbounded_diagnostic"),
        "6R.161 uses production k-best (legacy_1024), matching the 6R.160 VCF corpus"
    );

    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing P12 ref/BAM");
        return;
    }

    kv(
        "java_pin",
        "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77",
    );
    kv("variant", format!("2:{TARGET} {MERGED_REF}/{MERGED_ALT}"));
    kv("java_vcf_ad", "0,3");
    kv("rust_pre_overwrite_ad_6r160", "2,5");
    kv("downstream_sparse_overwrite", "AD=0,1 PL=45,3,0 NOT_CHASED");

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
    let covering = regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= TARGET
                && r.end.get() >= TARGET
        })
        .expect("ActiveFull covering 2:92316347");
    kv(
        "active_full",
        format!(
            "{}:{}-{}\textended={}-{}\tn_reads={}",
            covering.contig,
            covering.start.get(),
            covering.end.get(),
            covering.extended_start.get(),
            covering.extended_end.get(),
            covering.reads.len()
        ),
    );

    begin_hap_list_observe();
    begin_likelihood_pipeline_observe();
    begin_poorly_modeled_observe();
    let args = CallRegionArgs::strict_java();
    let outcome = HaplotypeCallerEngine::call_region(covering, &dict, &ref_fasta, &args)
        .expect("call")
        .expect("outcome");
    let hap_snaps = take_hap_list_snaps();
    let pipe_snaps = take_likelihood_pipeline_snaps();
    let poorly = take_poorly_modeled_observe();

    for snap in &hap_snaps {
        kv("hap_list", format!("stage={} n={}", snap.stage, snap.n));
        for c in &snap.columns {
            kv(
                "hap_col",
                format!(
                    "stage={}\tidx={}\tisRef={}\tlen={}\tfnv={:016x}\tloc={}-{}\tcigar={}",
                    snap.stage,
                    c.index,
                    c.is_reference,
                    c.len,
                    c.fnv1a,
                    c.loc_start,
                    c.loc_end,
                    c.cigar
                ),
            );
        }
    }
    for snap in &pipe_snaps {
        kv(
            "ll_pipeline",
            format!(
                "stage={}\tn_reads={}\tn_haps={}\tn_ll={}",
                snap.stage, snap.n_reads, snap.n_haps, snap.n_ll_entries
            ),
        );
    }
    let mut poorly_drop = 0usize;
    for row in &poorly {
        kv(
            "poorly_modeled",
            format!(
                "qname={}\tflags={}\tmax_ll={:.6}\tthresh={:.6}\tjava_keep={}\trust_keep={}\textra_retain={}",
                row.qname,
                row.flags,
                row.max_ll,
                row.threshold,
                row.java_equiv_keep,
                row.rust_keep,
                row.extra_retain
            ),
        );
        if !row.rust_keep {
            poorly_drop += 1;
        }
    }
    kv("poorly_modeled_dropped", format!("{poorly_drop}"));

    let haps = &outcome.assembly.haplotypes;
    let (full_ref, full_pad) = outcome.assembly.event_map_reference();
    let apply_pad = outcome
        .assembly
        .haplotypes
        .iter()
        .find(|h| h.is_reference)
        .and_then(|h| h.genome_loc.as_ref().map(|g| g.start.get()))
        .unwrap_or(full_pad);
    let loc_off = TARGET.saturating_sub(apply_pad) as usize;
    kv(
        "haplotype_population",
        format!(
            "n_haps={} n_pairhmm_cells={} pad={apply_pad} geno_reads={}",
            haps.len(),
            outcome.read_likelihoods.len(),
            outcome.genotyping_reads.len()
        ),
    );
    for (i, h) in haps.iter().enumerate() {
        let base = h
            .bases
            .get(loc_off)
            .map(|b| (*b as char).to_ascii_uppercase())
            .unwrap_or('?');
        kv(
            "haplotype",
            format!(
                "idx={i}\tisRef={}\tlen={}\tfnv={}\tbase_at_target={base}",
                h.is_reference,
                h.bases.len(),
                fnv1a64_hex(&h.bases)
            ),
        );
    }

    let at_pos: Vec<_> = outcome
        .assembly
        .variation_events()
        .iter()
        .filter(|e| e.start_1based == GenomePosition::new_1based(TARGET))
        .collect();
    kv(
        "eventmap_at_target",
        at_pos
            .iter()
            .map(|e| format!("{}/{}", e.ref_allele, e.alt_allele))
            .collect::<Vec<_>>()
            .join(","),
    );

    let event = VariationEvent::from_alleles("2", TARGET, MERGED_REF, MERGED_ALT);
    let hap_cache = build_per_haplotype_variation_events(
        haps,
        full_ref,
        full_pad,
        outcome.assembly.max_mnp_distance(),
        "2",
    );
    let mapping = create_allele_mapper_with_events(
        &event,
        TARGET,
        haps,
        apply_pad,
        outcome.assembly.reference_bases(),
        outcome.assembly.max_mnp_distance(),
        true,
        Some(&hap_cache),
    );
    kv(
        "allele_mapper",
        format!(
            "n_ref_haps={} n_alt_haps={}",
            mapping.ref_haplotype_indices.len(),
            mapping.alt_haplotype_indices.len()
        ),
    );

    let rust_rows = region_likelihoods_to_rows(&outcome.read_likelihoods, haps.len());
    let rust_marg = marginalize_rows_to_biallelic_alleles(
        &rust_rows,
        &mapping.ref_haplotype_indices,
        &mapping.alt_haplotype_indices,
    );

    let mut ad_ref = 0i32;
    let mut ad_alt = 0i32;
    let mut n_overlap = 0usize;
    let mut n_informative = 0usize;
    println!("6R161\tvote_header\tqname\tflags\tpileup\tretain\tlr\tla\tbest\tsecond\tconf\tinformative\tad_contrib");
    for row in &rust_marg {
        let Some(rec) = outcome.genotyping_reads.get(row.read_index) else {
            continue;
        };
        let rec = rec.as_ref();
        let retain = java_alignment_read_overlaps_interval(
            rec,
            TARGET,
            TARGET,
            DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN,
        );
        if !retain {
            continue;
        }
        n_overlap += 1;
        let lr = row.haplotype_log10_likelihoods[0];
        let la = row.haplotype_log10_likelihoods[1];
        let (best_i, second_i, conf, inf) = java_best_breaking_ties(&[lr, la], 0);
        let pileup = snp_base_at(rec)
            .map(|b| (b as char).to_string())
            .unwrap_or_else(|| ".".into());
        let contrib = if inf {
            n_informative += 1;
            if best_i == 0 {
                ad_ref += 1;
                "REF"
            } else {
                ad_alt += 1;
                "ALT"
            }
        } else {
            "."
        };
        let qname = String::from_utf8_lossy(rec.qname());
        println!(
            "6R161\tvote\t{qname}\t{}\t{pileup}\t{retain}\t{lr:.8}\t{la:.8}\t{}\t{}\t{conf:.8}\t{inf}\t{contrib}",
            rec.flags(),
            allele_name(best_i),
            allele_name(second_i),
        );
        let rust_gap = (lr - la).abs();
        let rust_vote = if rust_gap > LOG_10_INFORMATIVE_THRESHOLD {
            if lr > la {
                "REF"
            } else {
                "ALT"
            }
        } else {
            "."
        };
        if rust_vote != contrib {
            kv(
                "vote_semantics_mismatch",
                format!(
                    "{qname}\tflags={}\tjava_tie={contrib}\trust_abs={rust_vote}",
                    rec.flags()
                ),
            );
        }
    }
    kv(
        "informative_AD_retainEvidence",
        format!("n_overlap={n_overlap} n_informative={n_informative} AD={ad_ref},{ad_alt}"),
    );
    kv("production_change", "NONE");
}
