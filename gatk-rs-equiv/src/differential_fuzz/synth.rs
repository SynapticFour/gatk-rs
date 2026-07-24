//! Materialize a [`Scenario`] into FASTA + BAM (+ index) under a work directory.

use super::scenario::Scenario;
use anyhow::{bail, Context, Result};
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

const CONTIG: &str = "synth1";

#[derive(Debug, Clone)]
pub struct Materialized {
    pub dir: PathBuf,
    pub reference: PathBuf,
    pub bam: PathBuf,
    pub interval: String,
}

struct Rng(u64);
impl Rng {
    fn next_u32(&mut self) -> u32 {
        // xorshift64*
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        ((x.wrapping_mul(0x2545F4914F6CDD1D)) >> 32) as u32
    }
    fn chance(&mut self, p255: u8) -> bool {
        (self.next_u32() % 256) < u32::from(p255)
    }
    fn gen_range(&mut self, lo: u32, hi: u32) -> u32 {
        if hi <= lo {
            return lo;
        }
        lo + self.next_u32() % (hi - lo + 1)
    }
}

fn base_at(rng: &mut Rng) -> u8 {
    b"ACGT"[rng.gen_range(0, 3) as usize]
}

fn mutate_base(b: u8, rng: &mut Rng) -> u8 {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut out = bases[rng.gen_range(0, 3) as usize];
    if out == b {
        out = bases[((rng.next_u32() as usize) + 1) % 4];
    }
    out
}

fn phred_char(mean: u8, jitter: u8, rng: &mut Rng) -> u8 {
    let j = if jitter == 0 {
        0
    } else {
        rng.gen_range(0, u32::from(jitter)) as i32 - (jitter as i32 / 2)
    };
    let q = (i32::from(mean) + j).clamp(2, 40) as u8;
    q + 33
}

/// Write FASTA + BAM for `scenario` into `dir`.
pub fn materialize_scenario(scenario: &Scenario, dir: &Path) -> Result<Materialized> {
    fs::create_dir_all(dir)?;
    let mut rng = Rng(scenario.seed ^ 0x9E37_79B9_7F4A_7C15);

    // Reference with mild complexity (avoid pure homopolymer).
    let mut reference = Vec::with_capacity(scenario.ref_len as usize);
    for i in 0..scenario.ref_len {
        if i > 0 && i % 17 == 0 {
            reference.push(reference[i as usize - 1]); // short homopolymer run
        } else {
            reference.push(base_at(&mut rng));
        }
    }
    if scenario.plant_snp && scenario.ref_len > 40 {
        let mid = (scenario.ref_len / 2) as usize;
        reference[mid] = mutate_base(reference[mid], &mut rng);
    }
    if scenario.plant_indel && scenario.ref_len > 60 {
        let mid = (scenario.ref_len / 2 + 7) as usize;
        // Delete one base in haplotype space by rewriting flanking — keep ref as-is;
        // plant indel via reads instead (below).
        let _ = mid;
    }

    let fa = dir.join("reference.fa");
    {
        let mut f = fs::File::create(&fa)?;
        writeln!(f, ">{CONTIG}")?;
        for chunk in reference.chunks(80) {
            f.write_all(chunk)?;
            f.write_all(b"\n")?;
        }
    }
    // faidx via samtools when available; else write trivial.fai
    if Command::new("samtools")
        .args(["faidx", fa.to_str().unwrap()])
        .status()
        .map(|s| s.success())
        .unwrap_or(false)
    {
        // ok
    } else {
        let fai = dir.join("reference.fa.fai");
        fs::write(
            fai,
            format!(
                "{CONTIG}\t{}\t{}\t80\t81\n",
                scenario.ref_len,
                CONTIG.len() + 2
            ),
        )?;
    }

    let sam = dir.join("reads.sam");
    let mut sam_f = fs::File::create(&sam)?;
    writeln!(sam_f, "@HD\tVN:1.6\tSO:coordinate")?;
    writeln!(sam_f, "@SQ\tSN:{CONTIG}\tLN:{}", scenario.ref_len)?;
    writeln!(sam_f, "@RG\tID:synth\tSM:SYNTH\tPL:ILLUMINA")?;

    let max_start = scenario.ref_len.saturating_sub(scenario.read_len).max(1);
    let n = scenario.n_reads.max(scenario.coverage);

    for i in 0..n {
        let start0 = rng.gen_range(0, max_start); // 0-based
        let pos1 = start0 + 1; // SAM 1-based
        let mapq =
            scenario.mapq_min + (rng.gen_range(0, u32::from(scenario.mapq_span.max(1))) as u8);
        let mapq = mapq.min(60);

        let mut seq = reference[start0 as usize..(start0 + scenario.read_len) as usize].to_vec();
        let mut cigar = format!("{}M", scenario.read_len);

        // Soft-clip
        if rng.chance(scenario.softclip_p) && scenario.softclip_max > 0 && scenario.read_len > 20 {
            let sc = rng
                .gen_range(1, u32::from(scenario.softclip_max))
                .min(scenario.read_len / 4);
            for b in seq.iter_mut().take(sc as usize) {
                *b = base_at(&mut rng);
            }
            cigar = format!("{sc}S{}M", scenario.read_len - sc);
        }

        // Small indel in read vs reference
        if rng.chance(scenario.indel_p) && scenario.indel_max > 0 && seq.len() > 20 {
            let k = rng.gen_range(1, u32::from(scenario.indel_max)) as usize;
            let at = seq.len() / 2;
            if rng.chance(128) {
                // insertion
                for _ in 0..k {
                    seq.insert(at, base_at(&mut rng));
                }
                cigar = format!("{}M{}I{}M", at, k, seq.len() - at - k);
            } else if at + k < seq.len() {
                seq.drain(at..at + k);
                cigar = format!("{}M{}D{}M", at, k, seq.len() - at);
            }
        }

        // Plant indel haplotype near mid for some reads
        if scenario.plant_indel && i % 3 == 0 && seq.len() > 30 {
            let at = seq.len() / 2;
            seq.insert(at, b'A');
            // crude CIGAR rebuild for planted insertion
            cigar = format!("{}M1I{}M", at, seq.len() - at - 1);
        }

        let qual: String = (0..seq.len())
            .map(|_| phred_char(scenario.bq_mean, scenario.bq_jitter, &mut rng) as char)
            .collect();
        let seq_s =
            String::from_utf8(seq).unwrap_or_else(|_| "N".repeat(scenario.read_len as usize));

        let flag = if rng.chance(scenario.mate_overlap_p) {
            // paired, proper pair, mate reverse (simplified)
            99u16
        } else {
            0u16
        };
        let qname = format!("r{i:04}");
        let mate_pos = if flag != 0 {
            (pos1 as i64 + i32::from(scenario.read_len as i16 / 2) as i64).max(1)
        } else {
            0
        };
        let tlen = if flag != 0 {
            i32::from(scenario.read_len as i16)
        } else {
            0
        };

        writeln!(
            sam_f,
            "{qname}\t{flag}\t{CONTIG}\t{pos1}\t{mapq}\t{cigar}\t{}\t{mate_pos}\t{tlen}\t{seq_s}\t{qual}\tRG:Z:synth",
            if flag != 0 { CONTIG } else { "*" }
        )?;

        // overlapping mate
        if flag != 0 {
            let mate_start = (start0 + scenario.read_len / 3).min(max_start);
            let mpos1 = mate_start + 1;
            let mut mseq = reference[mate_start as usize
                ..(mate_start + scenario.read_len).min(scenario.ref_len) as usize]
                .to_vec();
            while mseq.len() < scenario.read_len as usize {
                mseq.push(b'N');
            }
            let mqual: String = (0..mseq.len())
                .map(|_| phred_char(scenario.bq_mean, scenario.bq_jitter, &mut rng) as char)
                .collect();
            let mseq_s = String::from_utf8_lossy(&mseq).into_owned();
            let mcigar = format!("{}M", mseq.len());
            writeln!(
                sam_f,
                "{qname}\t147\t{CONTIG}\t{mpos1}\t{mapq}\t{mcigar}\t{CONTIG}\t{pos1}\t{}\t{mseq_s}\t{mqual}\tRG:Z:synth",
                -(tlen)
            )?;
        }
    }
    drop(sam_f);

    let bam_unsorted = dir.join("reads.unsorted.bam");
    let bam = dir.join("reads.bam");
    let status = Command::new("samtools")
        .args([
            "view",
            "-bS",
            sam.to_str().unwrap(),
            "-o",
            bam_unsorted.to_str().unwrap(),
        ])
        .status()
        .context("samtools view (SAM→BAM) — samtools required on PATH")?;
    if !status.success() {
        bail!("samtools view failed: {status}");
    }
    let status = Command::new("samtools")
        .args([
            "sort",
            "-o",
            bam.to_str().unwrap(),
            bam_unsorted.to_str().unwrap(),
        ])
        .status()
        .context("samtools sort")?;
    if !status.success() {
        bail!("samtools sort failed: {status}");
    }
    let _ = fs::remove_file(&bam_unsorted);
    let status = Command::new("samtools")
        .args(["index", bam.to_str().unwrap()])
        .status()
        .context("samtools index")?;
    if !status.success() {
        bail!("samtools index failed: {status}");
    }

    let interval = format!("{CONTIG}:1-{}", scenario.ref_len);
    Ok(Materialized {
        dir: dir.to_path_buf(),
        reference: fa,
        bam,
        interval,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::differential_fuzz::scenario::scenario_from_bytes;
    use tempfile::tempdir;

    #[test]
    fn materialize_smoke() {
        if Command::new("samtools").arg("--version").output().is_err() {
            eprintln!("skip: samtools missing");
            return;
        }
        let dir = tempdir().unwrap();
        let s = scenario_from_bytes(&[9, 8, 7, 6, 5, 4, 3, 2, 1, 0]);
        let m = materialize_scenario(&s, dir.path()).unwrap();
        assert!(m.bam.is_file());
        assert!(m.reference.is_file());
    }
}
