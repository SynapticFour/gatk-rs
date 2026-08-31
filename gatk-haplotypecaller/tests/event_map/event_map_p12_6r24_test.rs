//! 6R.24 TEST-ONLY: Java oracle provenance for mid-B C/A vs Rust 2-read island.
//! Does not change production k, unique-kmer gates, dangling, EventMap, or W-H1.

#[cfg(test)]
mod traces {
    use crate::assembly::{AssemblyGraphParams, AssemblyRead};
    use crate::assembly_region_finalize::{
        assembly_reference_read, finalize_region_reads_for_assembly,
        gatk_min_tail_quality_for_assembly, padded_reference_loc, records_to_assembly_reads,
        GATK_REFERENCE_PADDING_FOR_ASSEMBLY,
    };
    use crate::bio_ids::KmerSize;
    use crate::event_map::collect_variation_events;
    use crate::haplotype::Haplotype;
    use crate::read_model::{passes_hc_read_filters, ReadFilterParams};
    use crate::read_projection::query_index_at_reference_position;
    use crate::read_threading_graph::{
        assembly_graph_from_ref_and_reads_threading_with_summary, reference_has_non_unique_kmers,
    };
    use crate::walker_apply::{call_disposition, AssemblyRegionCallDisposition};
    use gatk_core::reference::{
        parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary,
    };
    use rust_htslib::bam::record::CigarString;
    use rust_htslib::bam::{IndexedReader, Read};
    use std::collections::HashSet;
    use std::path::Path;

    const SITE_CA: u64 = 92_317_399;
    const SITE_TC: u64 = 92_317_407;
    const SITE_GC: u64 = 92_317_412;
    const JAVA_ACTIVE: (u64, u64) = (92_317_262, 92_317_491);
    const JAVA_EXTENDED: (u64, u64) = (92_317_162, 92_317_591);

    fn fixture_paths() -> Option<(std::path::PathBuf, std::path::PathBuf)> {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
        let ref_path = root.join("parity/realworld/assets/hs37d5.simple.fa");
        let bam = root.join("parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam");
        if !ref_path.is_file() || !bam.is_file() {
            return None;
        }
        Some((ref_path, bam))
    }

    fn unique_kmers(bases: &[u8], k: usize) -> bool {
        !reference_has_non_unique_kmers(
            &AssemblyRead {
                bases: bases.to_vec(),
                base_quals: vec![30; bases.len()],
            },
            k,
        )
    }

    fn n_identical_kmers(seq: &[u8], ref_bases: &[u8], k: usize) -> (usize, usize) {
        if seq.len() < k || ref_bases.len() < k {
            return (0, 0);
        }
        let ref_set: HashSet<&[u8]> = (0..=ref_bases.len() - k)
            .map(|i| &ref_bases[i..i + k])
            .collect();
        let n = seq.len() - k + 1;
        let ident = (0..n).filter(|&i| ref_set.contains(&seq[i..i + k])).count();
        (n, ident)
    }

    fn graph_params(k: usize) -> AssemblyGraphParams {
        AssemblyGraphParams {
            kmer_size: KmerSize::try_from_usize(k).expect("k"),
            min_base_quality: 10,
            min_edge_weight: 1,
            dangling_path_max_nodes: 0,
            max_haplotypes: 128,
            max_haplotype_bases: 4096,
            start_threading_only_at_existing_vertex: false,
        }
    }

    #[test]
    fn six_r24_mid_b_oracle_provenance_audit() {
        let Some((ref_fasta, bam_path)) = fixture_paths() else {
            eprintln!("Real-data mid-B comparison unavailable");
            return;
        };

        eprintln!("ORACLE_FIXTURE source=parity/fixtures/p12-java-format/all_sites.tsv");
        eprintln!(
            "ORACLE_GENERATOR scripts/parity/run_p12_realworld_na12878_20k.sh docker gatk:4.4.0.0"
        );
        eprintln!("ORACLE_INTERVAL default P12_INTERVAL=2:92300000-92350000");
        eprintln!(
            "ORACLE_BAM parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam (same as Rust tests)"
        );
        eprintln!("LIVE_JAVA_RERUN docker broadinstitute/gatk:4.4.0.0 -L 2:92317000-92319000");
        eprintln!(
            "LIVE_JAVA_SITES 92317399 C/A QUAL=78.32 GT=1/1 AD=0,2 DP=2 GQ=6 PL=90,6,0 MQ=27.00"
        );
        eprintln!(
            "LIVE_JAVA_ASSEMBLY active=2:{}-{} nReads=2 extended=2:{}-{} k=25 (k=10 skipped non-unique REF haplotype)",
            JAVA_ACTIVE.0, JAVA_ACTIVE.1, JAVA_EXTENDED.0, JAVA_EXTENDED.1
        );

        let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
        let walk_iv = parse_intervals_cli_string(&dict, "2:92317000-92319000").expect("iv");
        let p12_iv = parse_intervals_cli_string(&dict, "2:92300000-92350000").expect("p12");
        let filters = ReadFilterParams::gatk_standard_hc();
        let cfg =
            crate::walker_traversal::WalkerTraversalConfig::gatk_haplotype_caller_production(100);

        let walk_small = crate::walker_traversal::traverse_assembly_region_walker(
            &dict, &walk_iv, &ref_fasta, &bam_path, &filters, &cfg,
        )
        .expect("walk small");
        let walk_p12 = crate::walker_traversal::traverse_assembly_region_walker(
            &dict, &p12_iv, &ref_fasta, &bam_path, &filters, &cfg,
        )
        .expect("walk p12");

        let region = crate::walker_traversal::flatten_assembly_regions(&walk_small)
            .into_iter()
            .find(|r| {
                matches!(
                    call_disposition(r),
                    AssemblyRegionCallDisposition::ActiveFull
                ) && r.start.get() <= SITE_CA
                    && r.end.get() >= SITE_CA
            })
            .expect("ActiveFull mid-B");
        let region_p12 = crate::walker_traversal::flatten_assembly_regions(&walk_p12)
            .into_iter()
            .find(|r| {
                matches!(
                    call_disposition(r),
                    AssemblyRegionCallDisposition::ActiveFull
                ) && r.start.get() <= SITE_CA
                    && r.end.get() >= SITE_CA
            })
            .expect("ActiveFull mid-B on P12 interval");

        eprintln!(
            "RUST_REGION walker=2:92317000-92319000 active={}-{} extended={}-{} n_raw_reads={} contig={}",
            region.start.get(),
            region.end.get(),
            region.extended_start.get(),
            region.extended_end.get(),
            region.reads.len(),
            region.contig
        );
        eprintln!(
            "RUST_REGION_P12_INTERVAL active={}-{} extended={}-{} n_raw_reads={}",
            region_p12.start.get(),
            region_p12.end.get(),
            region_p12.extended_start.get(),
            region_p12.extended_end.get(),
            region_p12.reads.len()
        );
        assert_eq!(region.start.get(), JAVA_ACTIVE.0);
        assert_eq!(region.end.get(), JAVA_ACTIVE.1);
        assert_eq!(region.extended_start.get(), JAVA_EXTENDED.0);
        assert_eq!(region.extended_end.get(), JAVA_EXTENDED.1);
        assert_eq!(region.reads.len(), 2);
        assert_eq!(
            region_p12.reads.len(),
            2,
            "full P12 -L does not add reads to this ActiveFull region"
        );

        let mut ref_cache = ReferenceWindowCache::new(ref_fasta.clone(), 4);
        let padded = assembly_reference_read(&dict, &mut ref_cache, &region).expect("pad");
        let (pad_start, _) = padded_reference_loc(&region, &dict);
        let ext_off = (region.extended_start.get() - pad_start) as usize;
        let ext_len = (region.extended_end.get() - region.extended_start.get() + 1) as usize;
        let ext_ref = &padded.bases[ext_off..ext_off + ext_len];
        eprintln!(
            "REF_WINDOWS pad_start={pad_start} pad_len={} ext_len={ext_len} java_hap_len=430 extra_pad={}",
            padded.bases.len(),
            GATK_REFERENCE_PADDING_FOR_ASSEMBLY
        );
        for k in [10usize, 25, 85] {
            eprintln!(
                "UNIQUE_KMER k={k} padded_1430_unique={} extended_430_unique={}",
                unique_kmers(&padded.bases, k),
                unique_kmers(ext_ref, k)
            );
        }
        assert!(
            !unique_kmers(&padded.bases, 25),
            "production unique-kmer gate sees non-unique k=25 on ±500 padded REF"
        );
        assert!(
            unique_kmers(ext_ref, 25),
            "Java createGraph REF haplotype (extended 430bp) has unique k=25"
        );
        assert!(
            !unique_kmers(ext_ref, 10),
            "Java skipped k=10 on the same 430bp haplotype"
        );

        let finalized = finalize_region_reads_for_assembly(
            &region.reads,
            &region,
            true,
            gatk_min_tail_quality_for_assembly(10),
            false,
        );
        let assembly_reads = records_to_assembly_reads(&finalized);
        assert_eq!(assembly_reads.len(), 2);

        let mut bam = IndexedReader::from_path(&bam_path).expect("bam");
        let _ = bam.set_reference(&ref_fasta);
        let tid = bam.header().tid(b"2").expect("contig 2") as u32;
        bam.fetch((tid, (SITE_CA - 1) as i64, SITE_GC as i64))
            .expect("fetch");
        let bam_hits: Vec<rust_htslib::bam::Record> =
            bam.records().filter_map(|r| r.ok()).collect();
        eprintln!("BAM_OVERLAPPING_SITES n={}", bam_hits.len());
        assert_eq!(bam_hits.len(), 2);

        for (i, rec) in bam_hits.iter().enumerate() {
            let qname = String::from_utf8_lossy(rec.qname()).into_owned();
            let raw_cigar = CigarString(rec.cigar().iter().copied().collect());
            let fin = &finalized[i];
            let fin_cigar = CigarString(fin.cigar().iter().copied().collect());
            let (n85, id85) = n_identical_kmers(&rec.seq().as_bytes(), &padded.bases, 85);
            let (n25, id25) = n_identical_kmers(&rec.seq().as_bytes(), ext_ref, 25);
            let class = if id25 > 0 && id25 < n25 {
                "BOTH"
            } else if id25 == n25 && n25 > 0 {
                "REF_CONNECTING"
            } else if n25 > 0 && id25 == 0 {
                "ALT_SUPPORTING"
            } else {
                "NEITHER"
            };
            eprintln!(
                "READ_UNIVERSE[{i}] qname={qname} flags={} mapq={} start1={} raw_cigar={raw_cigar} seq_len={} \
                 hc_filter={} fin_cigar={fin_cigar} fin_len={} raw_k85_ident={id85}/{n85} ext_k25_ident={id25}/{n25} class={class}",
                rec.flags(),
                rec.mapq(),
                rec.pos() + 1,
                rec.seq().len(),
                passes_hc_read_filters(rec, &filters),
                fin.seq().len()
            );
            for site in [SITE_CA, SITE_TC, SITE_GC] {
                let qi =
                    query_index_at_reference_position(rec.pos(), &raw_cigar, (site - 1) as i64);
                eprintln!("  overlap_site {site} qi={qi:?}");
            }
        }

        bam.fetch((tid, 92_316_661, 92_318_091)).expect("fetch pad");
        let mut extra = 0usize;
        let mut extra_pass = 0usize;
        for rec in bam.records().filter_map(|r| r.ok()) {
            extra += 1;
            if passes_hc_read_filters(&rec, &filters) {
                extra_pass += 1;
            }
        }
        eprintln!(
            "BAM_PADDED_SPAN 2:92316662-92318091 n_records={extra} n_pass_hc_filters={extra_pass} (includes the 2 site-overlapping reads)"
        );

        let ext_assembly_ref = AssemblyRead {
            bases: ext_ref.to_vec(),
            base_quals: vec![30; ext_ref.len()],
        };
        let (g25, _) = assembly_graph_from_ref_and_reads_threading_with_summary(
            &ext_assembly_ref,
            &assembly_reads,
            &graph_params(25),
        )
        .expect("k=25 on Java-length REF haplotype — test-only");
        let ref_set: HashSet<Vec<u8>> = (0..=ext_ref.len().saturating_sub(25))
            .map(|i| ext_ref[i..i + 25].to_vec())
            .collect();
        let mut n_alt = 0usize;
        let mut n_alt_on_ref_node = 0usize;
        for ar in &assembly_reads {
            if ar.bases.len() < 25 {
                continue;
            }
            for i in 0..=ar.bases.len() - 25 {
                let km = &ar.bases[i..i + 25];
                if ref_set.contains(km) {
                    continue;
                }
                n_alt += 1;
                if let Some(id) = g25.vertex_id_for_kmer(km) {
                    if g25.ref_nodes.contains(&id) {
                        n_alt_on_ref_node += 1;
                    }
                }
            }
        }
        eprintln!(
            "TEST_ONLY k=25 graph on 430bp REF haplotype nodes={} edges={} read_alt_25mers={n_alt} alt_on_ref_node={n_alt_on_ref_node}",
            g25.node_count(),
            g25.edge_count()
        );

        let mut ref_hap = Haplotype::new(ext_ref.to_vec(), true);
        let mut cig = crate::cigar::Cigar::new();
        cig.push(ext_ref.len(), crate::cigar::CigarOperator::Match);
        ref_hap.cigar = Some(cig);
        let events = collect_variation_events(
            std::slice::from_ref(&ref_hap),
            ext_ref,
            region.extended_start.get(),
            &region.contig,
            0,
        );
        eprintln!(
            "NOTE REF-only EventMap on 430bp hap n={} (expect 0); Java EventMap had C/A+T/C+G/C from k=25 ALT hap",
            events.len()
        );

        eprintln!(
            "PREPROCESS Java finalizeRegion: well-defined fragment => revertSoftClippedBases then assembler hardClipSoftClippedBases; \
             Rust finalize hard-clips S. Live Java still assembled k=25 430M haplotypes — soft-clip revert is not the C/A source."
        );
        eprintln!(
            "PILEUP usePileupDetection default false; Java log shows EventMap from k=25 haplotypes — NOT A CANDIDATE"
        );
        eprintln!(
            "CLASSIFICATION C: Rust unique-kmer gate uses ±500 padded REF (1430bp, k=25 non-unique) while Java createGraph uses refHaplotype (430bp, k=25 unique). Same BAM/interval/2 reads."
        );
    }
}
