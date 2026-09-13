//! 6R.165 live dump: FS / MQ / SOR annotation evidence at `2:92316347 G/A`.
//! 6R.165 documented pileup vs likelihoods INFO membership.
//! Production FS/SOR moved in 6R.166; MQ still pileup. Skipped unless `HOLDOUT_6R165=1`.
//!
//! ```text
//! HOLDOUT_6R165=1 cargo test -p gatk-haplotypecaller --test holdout_6r165_info_annotation -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::{build_per_haplotype_variation_events, VariationEvent};
use gatk_haplotypecaller::fragment_overlap::read_base_at_ref_coord_1based;
use gatk_haplotypecaller::hc_allele_mapping::create_allele_mapper_with_events;
use gatk_haplotypecaller::read_realignment::LOG_10_INFORMATIVE_THRESHOLD;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, region_likelihoods_to_rows,
    traverse_assembly_region_walker, try_emit_call_region_variants, AssemblyRegionCallDisposition,
    CallRegionArgs, GenomePosition, HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
    DEFAULT_STAND_EMIT_CONFIDENCE,
};
use rust_htslib::bam::Record;
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

const INTERVAL: &str = "2:92316200-92316580";
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 92_316_347;
const MERGED_REF: &str = "G";
const MERGED_ALT: &str = "A";
const FLAG_REVERSE: u16 = 0x10;
const FLAG_SECONDARY: u16 = 0x100;
const FLAG_DUP: u16 = 0x400;
const FLAG_SUPPLEMENTARY: u16 = 0x800;
const MQ_UNAVAILABLE: u8 = 255;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R165\t{key}\t{}", value.as_ref());
}

fn info_key(v: &InfoValue) -> &str {
    match v {
        InfoValue::Flag(k)
        | InfoValue::Integer(k, _)
        | InfoValue::Float(k, _)
        | InfoValue::String(k, _)
        | InfoValue::Character(k, _) => k.as_str(),
    }
}

fn read_id(rec: &Record) -> String {
    format!(
        "QNAME={} FLAG={} pos={} CIGAR={} MAPQ={}",
        String::from_utf8_lossy(rec.qname()),
        rec.flags(),
        rec.pos() + 1,
        rec.cigar(),
        rec.mapq()
    )
}

fn java_calculate_sor(ref_fw: u32, ref_rv: u32, alt_fw: u32, alt_rv: u32) -> f64 {
    let t00 = f64::from(ref_fw) + 1.0;
    let t01 = f64::from(ref_rv) + 1.0;
    let t10 = f64::from(alt_fw) + 1.0;
    let t11 = f64::from(alt_rv) + 1.0;
    let ratio = (t00 / t01) * (t11 / t10) + (t01 / t00) * (t10 / t11);
    let ref_ratio = t00.min(t01) / t00.max(t01);
    let alt_ratio = t10.min(t11) / t10.max(t11);
    ratio.ln() + ref_ratio.ln() - alt_ratio.ln()
}

#[test]
fn holdout_6r165_info_annotation() {
    if std::env::var("HOLDOUT_6R165").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R165=1");
        return;
    }
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
    kv("variant", format!("2:{TARGET} G/A"));
    kv("production_change", "NONE");

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
        .expect("ActiveFull covering target");

    let call_args = CallRegionArgs::strict_java();
    let outcome = HaplotypeCallerEngine::call_region(covering, &dict, &ref_fasta, &call_args)
        .expect("call")
        .expect("outcome");
    kv(
        "counts",
        format!(
            "region_reads={} geno_reads={} pairhmm_cells={} haps={}",
            covering.reads.len(),
            outcome.genotyping_reads.len(),
            outcome.read_likelihoods.len(),
            outcome.assembly.haplotypes.len()
        ),
    );

    for rec in covering.reads.iter().map(|r| r.as_ref()) {
        let pileup = read_base_at_ref_coord_1based(rec, TARGET as i32)
            .map(|b| (b as char).to_string())
            .unwrap_or_else(|| "-".into());
        kv(
            "region_read",
            format!(
                "{} strand={} pileup={pileup} dup={} sec={} sup={}",
                read_id(rec),
                if rec.flags() & FLAG_REVERSE != 0 {
                    "rev"
                } else {
                    "fwd"
                },
                rec.flags() & FLAG_DUP != 0,
                rec.flags() & FLAG_SECONDARY != 0,
                rec.flags() & FLAG_SUPPLEMENTARY != 0
            ),
        );
    }

    let call = outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(TARGET)
                && c.event.ref_allele == MERGED_REF
                && c.event.alt_allele == MERGED_ALT
        })
        .expect("genotyped G/A");
    let fmt = &call.genotype.format;
    kv(
        "format",
        format!(
            "AD={:?} PL={:?} GT_closed=1/1",
            fmt.ad_as_i32(),
            fmt.pl_as_i32()
        ),
    );

    let emitted =
        try_emit_call_region_variants(covering, &outcome, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .expect("emit");
    let rec = emitted
        .iter()
        .find(|r| {
            r.position == TARGET
                && r.reference == MERGED_REF
                && r.alternate.first().map(String::as_str) == Some(MERGED_ALT)
        })
        .expect("emitted G/A");
    kv("emitted_info", format!("{:?}", rec.info));
    let extra: Vec<_> = rec
        .info
        .iter()
        .map(info_key)
        .filter(|k| {
            !matches!(
                *k,
                "AC" | "AF"
                    | "AN"
                    | "DP"
                    | "ExcessHet"
                    | "FS"
                    | "MLEAC"
                    | "MLEAF"
                    | "MQ"
                    | "QD"
                    | "SOR"
            )
        })
        .collect();
    kv("rust_only_info_tags", extra.join(","));

    let (ann_fs, ann_mq, ann_sor, ann_qual) = {
        let mut fs = 0.0;
        let mut mq = 0.0;
        let mut sor = 0.0;
        for v in &rec.info {
            match v {
                InfoValue::Float(k, xs) if k == "FS" => fs = xs[0],
                InfoValue::Float(k, xs) if k == "MQ" => mq = xs[0],
                InfoValue::Float(k, xs) if k == "SOR" => sor = xs[0],
                _ => {}
            }
        }
        (fs, mq, sor, rec.quality.unwrap_or(0.0))
    };
    kv(
        "rust_final",
        format!(
            "FS={:.5} MQ={} SOR={:.5} QUAL={:.5}",
            ann_fs, ann_mq, ann_sor, ann_qual
        ),
    );
    kv("java_final", "FS=0 MQ=40.25 SOR=1.179 QUAL=121.84");

    let ref_b = MERGED_REF.as_bytes()[0];
    let alt_b = MERGED_ALT.as_bytes()[0];
    let mut rf = 0u32;
    let mut rr = 0u32;
    let mut af = 0u32;
    let mut ar = 0u32;
    let mut alt_mq = Vec::new();
    for rec in covering.reads.iter().map(|r| r.as_ref()) {
        let Some(base) = read_base_at_ref_coord_1based(rec, TARGET as i32) else {
            continue;
        };
        let reverse = rec.flags() & FLAG_REVERSE != 0;
        if base.eq_ignore_ascii_case(&alt_b) {
            if reverse {
                ar += 1;
            } else {
                af += 1;
            }
            alt_mq.push(rec.mapq());
            kv("rust_fs_member", format!("{} allele=ALT", read_id(rec)));
        } else if base.eq_ignore_ascii_case(&ref_b) {
            if reverse {
                rr += 1;
            } else {
                rf += 1;
            }
            kv("rust_fs_member", format!("{} allele=REF", read_id(rec)));
        }
    }
    kv(
        "rust_strand_table",
        format!(
            "[{rf},{rr};{af},{ar}] SOR={}",
            java_calculate_sor(rf, rr, af, ar)
        ),
    );
    let rust_mq = if alt_mq.is_empty() {
        0.0
    } else {
        f64::from(alt_mq.iter().map(|&m| u32::from(m)).sum::<u32>()) / alt_mq.len() as f64
    };
    kv(
        "rust_mq_agg",
        format!("alt_n={} mapqs={:?} mean={rust_mq}", alt_mq.len(), alt_mq),
    );

    let haps = &outcome.assembly.haplotypes;
    let (full_ref, full_pad) = outcome.assembly.event_map_reference();
    let apply_pad = haps
        .iter()
        .find(|h| h.is_reference)
        .and_then(|h| h.genome_loc.as_ref().map(|g| g.start.get()))
        .unwrap_or(full_pad);
    let hap_cache = build_per_haplotype_variation_events(
        haps,
        full_ref,
        full_pad,
        outcome.assembly.max_mnp_distance(),
        "2",
    );
    let mapping = create_allele_mapper_with_events(
        &VariationEvent::from_alleles("2", TARGET, MERGED_REF, MERGED_ALT),
        TARGET,
        haps,
        apply_pad,
        outcome.assembly.reference_bases(),
        outcome.assembly.max_mnp_distance(),
        true,
        Some(&hap_cache),
    );
    let ref_haps: Vec<usize> = mapping
        .ref_haplotype_indices
        .iter()
        .map(|i| i.get())
        .collect();
    let alt_haps: Vec<usize> = mapping
        .alt_haplotype_indices
        .iter()
        .map(|i| i.get())
        .collect();
    kv(
        "allele_mapper",
        format!("ref_haps={ref_haps:?} alt_haps={alt_haps:?}"),
    );

    let geno: Vec<Record> = outcome
        .genotyping_reads
        .iter()
        .map(|r| (**r).clone())
        .collect();
    let rows = region_likelihoods_to_rows(&outcome.read_likelihoods, haps.len());
    let mut jrf = 0u32;
    let mut jrr = 0u32;
    let mut jaf = 0u32;
    let mut jar = 0u32;
    let mut seen = BTreeSet::new();
    let mut mq_list = Vec::new();
    for row in &rows {
        if seen.insert(row.read_index) {
            if let Some(rec) = geno.get(row.read_index) {
                if rec.mapq() != MQ_UNAVAILABLE {
                    mq_list.push(rec.mapq());
                }
                kv(
                    "java_mq_member",
                    format!(
                        "{} strand={}",
                        read_id(rec),
                        if rec.flags() & FLAG_REVERSE != 0 {
                            "rev"
                        } else {
                            "fwd"
                        }
                    ),
                );
            }
        }
        let Some(rec) = geno.get(row.read_index) else {
            continue;
        };
        let lls = &row.haplotype_log10_likelihoods;
        let ll_ref = ref_haps
            .iter()
            .filter_map(|&i| lls.get(i).copied())
            .fold(f64::NEG_INFINITY, f64::max);
        let ll_alt = alt_haps
            .iter()
            .filter_map(|&i| lls.get(i).copied())
            .fold(f64::NEG_INFINITY, f64::max);
        let (best_is_ref, best, second) = if ll_ref >= ll_alt {
            (true, ll_ref, ll_alt)
        } else {
            (false, ll_alt, ll_ref)
        };
        let gap = if second.is_finite() {
            best - second
        } else {
            f64::INFINITY
        };
        let informative = gap > LOG_10_INFORMATIVE_THRESHOLD;
        let allele = if gap < LOG_10_INFORMATIVE_THRESHOLD {
            "REF"
        } else if best_is_ref {
            "REF"
        } else {
            "ALT"
        };
        kv(
            "java_best_allele",
            format!(
                "{} allele={allele} informative={informative} gap={gap:.4} ll_ref={ll_ref:.4} ll_alt={ll_alt:.4}",
                read_id(rec)
            ),
        );
        if !informative {
            continue;
        }
        let reverse = rec.flags() & FLAG_REVERSE != 0;
        if allele == "REF" {
            if reverse {
                jrr += 1;
            } else {
                jrf += 1;
            }
        } else if reverse {
            jar += 1;
        } else {
            jaf += 1;
        }
    }
    kv(
        "java_strand_table",
        format!(
            "[{jrf},{jrr};{jaf},{jar}] FS_MIN_COUNT=2 SOR={}",
            java_calculate_sor(jrf, jrr, jaf, jar)
        ),
    );
    let n = mq_list.len();
    let sum_sq: u64 = mq_list.iter().map(|&m| u64::from(m) * u64::from(m)).sum();
    let rms = if n == 0 {
        0.0
    } else {
        (sum_sq as f64 / n as f64).sqrt()
    };
    kv(
        "java_mq_agg",
        format!("n={n} mapqs={mq_list:?} sum_sq={sum_sq} rms={rms:.5} printed={rms:.2}"),
    );
    kv(
        "first_fs_unequal",
        format!("rust=[{rf},{rr};{af},{ar}] java=[{jrf},{jrr};{jaf},{jar}]"),
    );
    kv(
        "first_mq_unequal",
        format!(
            "rust_alt_n={} rust_mean={rust_mq} java_n={n} java_rms={rms:.5}",
            alt_mq.len()
        ),
    );
}
