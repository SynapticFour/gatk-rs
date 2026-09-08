//! 6R.124 forensic dump: haplotype-likelihood matrix provenance at `20:29456196`.
//! Dumps the objects entering the likelihood engine (not PairHMM DP).
//! Skipped unless `HOLDOUT_6R124=1`. Production change: NONE.
//!
//! ```text
//! HOLDOUT_6R124=1 cargo test -p gatk-haplotypecaller --test holdout_6r124_ll_inputs -- --nocapture
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::hc_genotyping_engine::strict_java_pairhmm_normalize_hap_indices;
use gatk_haplotypecaller::pcr_error_model::apply_pcr_error_model;
use gatk_haplotypecaller::{
    begin_likelihood_pipeline_observe, call_disposition, flatten_assembly_regions,
    indel_gop_from_optional_tag, prepare_read_quals_for_pairhmm_inplace,
    region_likelihoods_to_rows, take_likelihood_pipeline_cells, take_likelihood_pipeline_snaps,
    traverse_assembly_region_walker, AssemblyRegionCallDisposition, CallRegionArgs,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig, GATK_PARITY_DEFAULT_GCP,
};
use rust_htslib::bam::record::{Aux, Cigar};
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 29_456_196;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn fnv1a64_hex(data: &[u8]) -> String {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in data {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    format!("{h:016x}")
}

fn fastq(q: &[u8]) -> String {
    q.iter()
        .map(|&b| char::from(b.saturating_add(33)))
        .collect()
}

fn bam_bqsr_indel_quals_phred(rec: &rust_htslib::bam::Record, tag: &[u8]) -> Option<Vec<u8>> {
    match rec.aux(tag) {
        Ok(Aux::String(s)) => Some(s.bytes().map(|b| b.saturating_sub(33)).collect()),
        _ => None,
    }
}

fn cigar_str(rec: &rust_htslib::bam::Record) -> String {
    format!("{}", rec.cigar())
}

fn alignment_end_1based(rec: &rust_htslib::bam::Record) -> i64 {
    let mut ref_len: i64 = 0;
    for c in rec.cigar().iter() {
        match c {
            Cigar::Match(n)
            | Cigar::Equal(n)
            | Cigar::Diff(n)
            | Cigar::Del(n)
            | Cigar::RefSkip(n) => ref_len += i64::from(*n),
            _ => {}
        }
    }
    rec.pos() + ref_len
}

fn leading_clip(rec: &rust_htslib::bam::Record) -> i64 {
    let mut n: i64 = 0;
    for c in rec.cigar().iter() {
        match c {
            Cigar::SoftClip(k) | Cigar::HardClip(k) => n += i64::from(*k),
            _ => break,
        }
    }
    n
}

fn trailing_clip(rec: &rust_htslib::bam::Record) -> i64 {
    let mut n: i64 = 0;
    for c in rec.cigar().iter().rev() {
        match c {
            Cigar::SoftClip(k) | Cigar::HardClip(k) => n += i64::from(*k),
            _ => break,
        }
    }
    n
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R124\t{key}\t{}", value.as_ref());
}

#[test]
fn holdout_6r124_ll_inputs() {
    if std::env::var("HOLDOUT_6R124").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R124=1");
        return;
    }
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
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
        .expect("ActiveFull");
    let args = CallRegionArgs::strict_java();
    begin_likelihood_pipeline_observe();
    let outcome = HaplotypeCallerEngine::call_region(covering, &dict, &ref_fasta, &args)
        .expect("call")
        .expect("outcome");
    let pipe_cells = take_likelihood_pipeline_cells();
    let pipe_snaps = take_likelihood_pipeline_snaps();
    let haps = &outcome.assembly.haplotypes;
    let cfg = &args.likelihood;
    kv(
        "active",
        format!(
            "{}:{}-{}",
            covering.contig,
            covering.start.get(),
            covering.end.get()
        ),
    );
    kv("hap_count", haps.len().to_string());
    kv(
        "orig_evidence_count",
        outcome.genotyping_reads.len().to_string(),
    );
    kv("pcr_error_model", format!("{:?}", cfg.pcr_error_model));
    kv("bq_threshold", cfg.base_quality_score_threshold.to_string());
    kv(
        "disable_cap_mapq",
        cfg.disable_cap_read_qualities_to_mapq.to_string(),
    );
    kv("constant_gcp", GATK_PARITY_DEFAULT_GCP.to_string());
    kv(
        "pairhmm_backend",
        format!("{:?}", cfg.resolved_pair_hmm_backend()),
    );

    for (i, h) in haps.iter().enumerate() {
        let loc = h
            .genome_loc
            .map(|g| format!("{}-{}", g.start_1based(), g.end_1based()))
            .unwrap_or_else(|| ".".to_string());
        let cigar = h
            .cigar
            .as_ref()
            .map(|c| c.to_gatk_string())
            .unwrap_or_else(|| ".".to_string());
        kv(
            "hap",
            format!(
                "{i}\t{}\tlen={}\tisRef={}\tloc={loc}\talignStart={}\tcigar={cigar}\tbases={}",
                fnv1a64_hex(&h.bases),
                h.bases.len(),
                h.is_reference,
                h.alignment_start_hap_wrt_ref,
                String::from_utf8_lossy(&h.bases)
            ),
        );
    }

    let pad = outcome.assembly.padded_reference_start_1based();
    kv("pad_start", pad.to_string());

    // 6R.127: live normalize-mask object (same function as engine.rs; TEST dump only).
    let ref_hap = haps.iter().find(|h| h.is_reference);
    let apply_pad = ref_hap
        .and_then(|h| h.genome_loc.map(|g| g.start_1based()))
        .unwrap_or(pad);
    let apply_bases = outcome.assembly.apply_bases_shared();
    let mut mask = strict_java_pairhmm_normalize_hap_indices(
        &outcome.assembly,
        haps,
        covering.start.get(),
        covering.end.get(),
        apply_pad,
        apply_bases.as_ref(),
        outcome.assembly.max_mnp_distance(),
        &covering.contig,
        &args.genotyping,
    );
    mask.sort_unstable();
    kv("norm_mask_n", mask.len().to_string());
    kv("likelihood_hap_n", haps.len().to_string());
    kv("production_norm_p_n", haps.len().to_string());
    for &i in &mask {
        let h = &haps[i];
        kv(
            "norm_mask_hap",
            format!(
                "idx={i}\thash={}\tisRef={}\tlen={}",
                fnv1a64_hex(&h.bases),
                h.is_reference,
                h.bases.len()
            ),
        );
    }
    let mut n_snp = 0u32;
    let mut n_indel = 0u32;
    let mut n_span = 0u32;
    for e in outcome.assembly.variation_events() {
        let in_span = e.start_1based.get() >= covering.start.get()
            && e.start_1based.get() <= covering.end.get();
        if in_span {
            n_span += 1;
        }
        let kind = if e.ref_allele.len() == 1 && e.alt_allele.len() == 1 {
            n_snp += 1;
            "snp"
        } else if e.is_indel() {
            n_indel += 1;
            "indel"
        } else {
            "other"
        };
        if in_span {
            kv(
                "norm_mask_event",
                format!(
                    "pos={}\tkind={kind}\tref={}\talt={}\tin_span=true",
                    e.start_1based.get(),
                    e.ref_allele,
                    e.alt_allele
                ),
            );
        }
    }
    kv(
        "norm_mask_events",
        format!("span={n_span}\tsnp_all={n_snp}\tindel_all={n_indel}"),
    );

    for (i, rec) in outcome.genotyping_reads.iter().enumerate() {
        let qname = String::from_utf8_lossy(rec.qname()).into_owned();
        let bases = rec.seq().as_bytes();
        let raw_bq = rec.qual();
        let ins_tag = bam_bqsr_indel_quals_phred(rec, b"BI");
        let del_tag = bam_bqsr_indel_quals_phred(rec, b"BD");
        let raw_iq = ins_tag.clone().unwrap_or_else(|| vec![45u8; bases.len()]);
        let raw_dq = del_tag.clone().unwrap_or_else(|| vec![45u8; bases.len()]);
        let start = rec.pos() + 1;
        let end = alignment_end_1based(rec);
        let u_start = start - leading_clip(rec);
        let u_end = end + trailing_clip(rec);
        kv(
            "orig",
            format!(
                "{i}\tqname={qname}\tflags={}\tstart={start}\tend={end}\tuStart={u_start}\tuEnd={u_end}\tmapq={}\tlen={}\tcigar={}\tbasesHash={}\tbqHash={}\tiqHash={}\tdqHash={}\tgcpHash=.\thmmHash=.\tbases={}\tbq={}\tiq={}\tdq={}\tgcp=.\thmmBq=.",
                rec.flags(),
                rec.mapq(),
                bases.len(),
                cigar_str(rec),
                fnv1a64_hex(&bases),
                fnv1a64_hex(raw_bq),
                fnv1a64_hex(&raw_iq),
                fnv1a64_hex(&raw_dq),
                String::from_utf8_lossy(&bases),
                fastq(raw_bq),
                fastq(&raw_iq),
                fastq(&raw_dq),
            ),
        );

        let n = bases.len();
        let mut bq = raw_bq.to_vec();
        prepare_read_quals_for_pairhmm_inplace(&mut bq, rec.mapq(), cfg);
        let mut ins = indel_gop_from_optional_tag(ins_tag.as_deref(), n).expect("ins");
        let mut del = indel_gop_from_optional_tag(del_tag.as_deref(), n).expect("del");
        if !cfg.uses_dragstr_pair_hmm() {
            apply_pcr_error_model(&bases, &mut ins, &mut del, cfg.pcr_error_model);
        }
        let gcp = vec![GATK_PARITY_DEFAULT_GCP; n];
        kv(
            "proc",
            format!(
                "{i}\tqname={qname}\tflags={}\tstart={start}\tend={end}\tuStart={u_start}\tuEnd={u_end}\tmapq={}\tlen={n}\tcigar={}\tbasesHash={}\tbqHash={}\tiqHash={}\tdqHash={}\tgcpHash={}\thmmHash={}\tbases={}\tbq={}\tiq={}\tdq={}\tgcp={}\thmmBq={}",
                rec.flags(),
                rec.mapq(),
                cigar_str(rec),
                fnv1a64_hex(&bases),
                fnv1a64_hex(&bq),
                fnv1a64_hex(&ins),
                fnv1a64_hex(&del),
                fnv1a64_hex(&gcp),
                fnv1a64_hex(&bq),
                String::from_utf8_lossy(&bases),
                fastq(&bq),
                fastq(&ins),
                fastq(&del),
                fastq(&gcp),
                fastq(&bq),
            ),
        );
    }

    let hap_hashes: Vec<String> = haps.iter().map(|h| fnv1a64_hex(&h.bases)).collect();
    let mut rows = region_likelihoods_to_rows(&outcome.read_likelihoods, haps.len());
    kv("stored_n_hap", haps.len().to_string());
    kv("stored_n_ev", rows.len().to_string());
    kv(
        "stored_evidence_count",
        outcome.genotyping_reads.len().to_string(),
    );
    for row in &mut rows {
        if let Some(rec) = outcome.genotyping_reads.get(row.read_index) {
            let qname = String::from_utf8_lossy(rec.qname()).into_owned();
            let start = rec.pos() + 1;
            for (a, &ll) in row.haplotype_log10_likelihoods.iter().enumerate() {
                kv(
                    "stored",
                    format!(
                        "{qname}\tflags={}\tstart={start}\thap={}\tll={ll:.12}",
                        rec.flags(),
                        hap_hashes[a]
                    ),
                );
            }
        }
    }
    kv("pipe_n", pipe_cells.len().to_string());
    for s in &pipe_snaps {
        kv(
            "snap",
            format!(
                "seq={}\tstage={}\tn_reads={}\tn_haps={}\tn_ll={}",
                s.seq, s.stage, s.n_reads, s.n_haps, s.n_ll_entries
            ),
        );
    }
    for c in &pipe_cells {
        kv(
            "pipe",
            format!(
                "stage={}\tqname={}\tflags={}\thap={:016x}\tll={:.12}\tbits={:016x}\tf32wide={}\tfinite={}\tread_index={}\thap_index={}",
                c.stage,
                c.qname,
                c.flags,
                c.hap_fnv,
                c.log10_likelihood,
                c.log10_likelihood.to_bits(),
                f64::from(c.log10_likelihood as f32).to_bits() == c.log10_likelihood.to_bits(),
                c.log10_likelihood.is_finite(),
                c.read_index,
                c.hap_index
            ),
        );
    }
}
