//! 6R.21 TEST-ONLY: mid-B assembly / read-threading audit.
//! Does not change production assembler, EventMap, k, waivers, or W-H1.

#[cfg(test)]
mod traces {
    use crate::assembly::{
        AssemblyGraph, AssemblyGraphParams, AssemblyGraphPruningParams, AssemblyRead,
    };
    use crate::assembly_based_caller::{assemble_reads_with_finalized, AssembleReadsArgs};
    use crate::assembly_region_finalize::{
        assembly_reference_read, finalize_region_reads_for_assembly,
        gatk_min_tail_quality_for_assembly, padded_reference_loc, records_to_assembly_reads,
    };
    use crate::bio_ids::KmerSize;
    use crate::cigar::CigarOperator;
    use crate::event_map::collect_variation_events;
    use crate::haplotype::Haplotype;
    use crate::kbest_haplotype::find_best_haplotypes_for_assembly;
    use crate::read_projection::query_index_at_reference_position;
    use crate::read_threading_assembler::{
        assemble_from_ref_and_reads, build_threading_graph_for_haplotype_dump,
        build_threading_graph_for_seq_assembly, extract_haplotypes_from_seq_kbest_paths,
        extract_rt_haplotypes_before_remove_paths, AssemblyScoringContext,
        ReadThreadingAssemblerArgs,
    };
    use crate::read_threading_graph::{
        assembly_graph_from_ref_and_reads_threading_with_summary, reference_has_non_unique_kmers,
    };
    use crate::seq_graph::SeqGraph;
    use crate::seq_kbest_haplotype::find_best_haplotypes_seq_graph;
    use crate::walker_apply::{call_disposition, AssemblyRegionCallDisposition};
    use gatk_core::reference::{
        parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary,
    };
    use rust_htslib::bam::record::CigarString;
    use std::collections::HashSet;
    use std::path::Path;

    const SITE: u64 = 92_317_399;
    const CTX: usize = 16;

    struct MidB {
        dict: SequenceDictionary,
        ref_fasta: std::path::PathBuf,
        region_start: u64,
        region_end: u64,
        ext_start: u64,
        ext_end: u64,
        pad_start: u64,
        contig: String,
        n_raw_reads: usize,
        reference: AssemblyRead,
        raw_reads: Vec<crate::shared_bam::SharedBamRecord>,
        finalized: Vec<rust_htslib::bam::Record>,
        assembly_reads: Vec<AssemblyRead>,
        region: crate::assembly_region_iterator::AssemblyRegion,
    }

    fn fixture_paths() -> Option<(std::path::PathBuf, std::path::PathBuf)> {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
        let ref_path = root.join("parity/realworld/assets/hs37d5.simple.fa");
        let bam = root.join("parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam");
        if !ref_path.is_file() || !bam.is_file() {
            return None;
        }
        Some((ref_path, bam))
    }

    fn load_mid_b() -> Option<MidB> {
        let (ref_fasta, bam) = fixture_paths()?;
        let dict = SequenceDictionary::from_fasta_path(&ref_fasta).ok()?;
        let specs = parse_intervals_cli_string(&dict, "2:92317000-92319000").ok()?;
        let walk = crate::walker_traversal::traverse_assembly_region_walker(
            &dict,
            &specs,
            &ref_fasta,
            &bam,
            &crate::read_model::ReadFilterParams::gatk_standard_hc(),
            &crate::walker_traversal::WalkerTraversalConfig::gatk_haplotype_caller_production(100),
        )
        .ok()?;
        let regions = crate::walker_traversal::flatten_assembly_regions(&walk);
        let region = regions
            .iter()
            .find(|r| {
                matches!(
                    call_disposition(r),
                    AssemblyRegionCallDisposition::ActiveFull
                ) && r.start.get() <= SITE
                    && r.end.get() >= SITE
            })?
            .clone();
        let mut ref_cache = ReferenceWindowCache::new(ref_fasta.clone(), 4);
        let reference = assembly_reference_read(&dict, &mut ref_cache, &region).ok()?;
        let (pad_start, _) = padded_reference_loc(&region, &dict);
        let finalized = finalize_region_reads_for_assembly(
            &region.reads,
            &region,
            true,
            gatk_min_tail_quality_for_assembly(10),
            false,
        );
        let assembly_reads = records_to_assembly_reads(&finalized);
        Some(MidB {
            dict,
            ref_fasta,
            region_start: region.start.get(),
            region_end: region.end.get(),
            ext_start: region.extended_start.get(),
            ext_end: region.extended_end.get(),
            pad_start,
            contig: region.contig.clone(),
            n_raw_reads: region.reads.len(),
            reference,
            raw_reads: region.reads.clone(),
            finalized,
            assembly_reads,
            region,
        })
    }

    fn ascii(bytes: &[u8]) -> String {
        String::from_utf8_lossy(bytes).into_owned()
    }

    fn site_offset(pad_start: u64) -> usize {
        (SITE.saturating_sub(pad_start)) as usize
    }

    fn overlapping_kmer_starts(seq_len: usize, off: usize, k: usize) -> Vec<usize> {
        if seq_len < k || off >= seq_len {
            return Vec::new();
        }
        let first = off.saturating_sub(k - 1);
        let last = off.min(seq_len - k);
        if first > last {
            Vec::new()
        } else {
            (first..=last).collect()
        }
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

    fn prod_args(mid: &MidB) -> ReadThreadingAssemblerArgs {
        let mut a = ReadThreadingAssemblerArgs::default();
        a.dangling_java_exact = true;
        a.scoring = Some(AssemblyScoringContext {
            padded_reference_start_1based: mid.pad_start,
            active_start_1based: mid.region_start,
            active_end_1based: mid.region_end,
            contig: mid.contig.clone(),
        });
        a
    }

    fn dump_bam_read(
        label: &str,
        rec: &rust_htslib::bam::Record,
        pad_start: u64,
        ref_bases: &[u8],
    ) {
        let cigar = CigarString(rec.cigar().iter().copied().collect());
        let qname = String::from_utf8_lossy(rec.qname()).into_owned();
        let start1 = rec.pos() as u64 + 1;
        let end1 = rec.cigar().end_pos() as u64;
        let seq = rec.seq().as_bytes();
        let quals = rec.qual();
        let site0 = (SITE - 1) as i64;
        let qi = query_index_at_reference_position(rec.pos(), &cigar, site0);
        eprintln!(
            "READ[{label}] qname={qname} start1={start1} end1={end1} cigar={cigar} mapq={} flags={} reverse={} seq_len={}",
            rec.mapq(),
            rec.flags(),
            rec.is_reverse(),
            seq.len()
        );
        match qi {
            Some(i) if i < seq.len() => {
                let q = if i < quals.len() { quals[i] } else { 255 };
                let lo = i.saturating_sub(CTX);
                let hi = (i + CTX + 1).min(seq.len());
                let qlo = lo.min(quals.len());
                let qhi = hi.min(quals.len());
                eprintln!(
                    "  LOCUS qi={i} base={} qual={q} usable_q10={} ctx={}-{} seq={} quals={:?}",
                    seq[i] as char,
                    q >= 10,
                    lo,
                    hi,
                    ascii(&seq[lo..hi]),
                    &quals[qlo..qhi]
                );
            }
            Some(i) => eprintln!("  LOCUS qi={i} OUT_OF_SEQ len={}", seq.len()),
            None => eprintln!("  LOCUS no query base (deletion / unaligned / clipped off)"),
        }
        let off = site_offset(pad_start);
        if off < ref_bases.len() {
            let rlo = off.saturating_sub(CTX);
            let rhi = (off + CTX + 1).min(ref_bases.len());
            eprintln!(
                "  REF_CTX off={off} base={} window={}",
                ref_bases[off] as char,
                ascii(&ref_bases[rlo..rhi])
            );
        }
    }

    fn kmer_presence(
        graph: &AssemblyGraph,
        kmers: &[Vec<u8>],
    ) -> (usize, Vec<(String, u32, Vec<u32>)>) {
        let mut present = 0usize;
        let mut rows = Vec::new();
        for kmer in kmers {
            if let Some(id) = graph.vertex_id_for_kmer(kmer) {
                present += 1;
                let support = graph.nodes()[id].support;
                let out_w: Vec<u32> = graph
                    .outgoing_nodes(id)
                    .into_iter()
                    .filter_map(|t| {
                        graph
                            .edges_sorted()
                            .into_iter()
                            .find(|e| e.from == id && e.to == t)
                            .map(|e| e.support)
                    })
                    .collect();
                rows.push((ascii(kmer), support, out_w));
            }
        }
        (present, rows)
    }

    fn snp_kmers(ref_bases: &[u8], off: usize, k: usize, alt: u8) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {
        let mut refs = Vec::new();
        let mut alts = Vec::new();
        for start in overlapping_kmer_starts(ref_bases.len(), off, k) {
            let mut rk = ref_bases[start..start + k].to_vec();
            refs.push(rk.clone());
            let pos = off - start;
            rk[pos] = alt;
            alts.push(rk);
        }
        (refs, alts)
    }

    fn hap_summary(label: &str, haps: &[Haplotype], ref_bases: &[u8], pad_start: u64) {
        let n_alt = haps.iter().filter(|h| !h.is_reference).count();
        let has_alt_bases = haps.iter().any(|h| h.bases.as_slice() != ref_bases);
        eprintln!(
            "HAPS[{label}] n={} n_alt_flag={n_alt} any_diff_bases={has_alt_bases}",
            haps.len()
        );
        for (i, h) in haps.iter().enumerate() {
            let cig = h.cigar.as_ref().map(|c| {
                c.elements
                    .iter()
                    .map(|e| {
                        let op = match e.operator {
                            CigarOperator::Match => 'M',
                            CigarOperator::Insertion => 'I',
                            CigarOperator::Deletion => 'D',
                            CigarOperator::SoftClip => 'S',
                            _ => '?',
                        };
                        format!("{}{op}", e.length)
                    })
                    .collect::<String>()
            });
            eprintln!(
                "  HAP[{i}] ref={} k={} len={} cigar={:?} align_off={}",
                h.is_reference,
                h.kmer_size,
                h.bases.len(),
                cig,
                h.alignment_start_hap_wrt_ref
            );
        }
        let events = collect_variation_events(haps, ref_bases, pad_start, "2", 0);
        eprintln!("  EVENTS n={}", events.len());
        for e in &events {
            eprintln!(
                "    {}-{} {}→{}",
                e.start_1based.get(),
                e.end_1based.get(),
                e.ref_allele,
                e.alt_allele
            );
        }
        let has_ca = events
            .iter()
            .any(|e| e.start_1based.get() == SITE && e.ref_allele == "C" && e.alt_allele == "A");
        eprintln!("  HAS_C_A={has_ca}");
    }

    fn unique_kmer_dup(seq: &[u8], k: usize) -> (usize, usize) {
        if seq.len() < k {
            return (0, 0);
        }
        let total = seq.len() - k + 1;
        let mut set = HashSet::new();
        for i in 0..total {
            set.insert(&seq[i..i + k]);
        }
        (set.len(), total)
    }

    fn audit_k(mid: &MidB, k: usize, args: &ReadThreadingAssemblerArgs) {
        let off = site_offset(mid.pad_start);
        let (uniq, total) = unique_kmer_dup(&mid.reference.bases, k);
        eprintln!(
            "K={k} REF_KMERS unique={uniq}/{total} non_unique={}",
            uniq < total
        );
        let (ref_ks, alt_ks) = snp_kmers(&mid.reference.bases, off, k, b'A');
        eprintln!(
            "  SNP_WINDOWS n={} first_ref={} first_alt={}",
            ref_ks.len(),
            ref_ks.first().map(|s| ascii(s)).unwrap_or_default(),
            alt_ks.first().map(|s| ascii(s)).unwrap_or_default()
        );

        let params = graph_params(k);
        let (raw, summary) = assembly_graph_from_ref_and_reads_threading_with_summary(
            &mid.reference,
            &mid.assembly_reads,
            &params,
        )
        .expect("raw graph");
        eprintln!(
            "  RAW nodes={} edges={} low_complexity={} cycles={}",
            raw.node_count(),
            raw.edge_count(),
            summary.is_low_complexity,
            raw.has_cycle()
        );
        let (n_ref, ref_rows) = kmer_presence(&raw, &ref_ks);
        let (n_alt, alt_rows) = kmer_presence(&raw, &alt_ks);
        eprintln!(
            "  RAW REF_KMER_NODES={n_ref}/{} ALT_KMER_NODES={n_alt}/{}",
            ref_ks.len(),
            alt_ks.len()
        );
        for (s, sup, out) in &alt_rows {
            eprintln!("    ALT_NODE {s} support={sup} out_w={out:?}");
        }
        if n_alt == 0 {
            for (s, sup, _) in ref_rows.iter().take(3) {
                eprintln!("    REF_NODE {s} support={sup}");
            }
        }

        for (ri, ar) in mid.assembly_reads.iter().enumerate() {
            let rec = &mid.finalized[ri];
            let cigar = CigarString(rec.cigar().iter().copied().collect());
            let Some(qi) = query_index_at_reference_position(rec.pos(), &cigar, (SITE - 1) as i64)
            else {
                eprintln!("  READ_KMER[{ri}] no locus qi");
                continue;
            };
            let starts = overlapping_kmer_starts(ar.bases.len(), qi, k);
            let mut present = 0usize;
            let mut sample = None;
            for st in &starts {
                let km = ar.bases[*st..*st + k].to_vec();
                if raw.vertex_id_for_kmer(&km).is_some() {
                    present += 1;
                    if sample.is_none() {
                        sample = Some(ascii(&km));
                    }
                }
            }
            eprintln!(
                "  READ_KMER[{ri}] qi={qi} windows={} in_raw_graph={present} sample={}",
                starts.len(),
                sample.unwrap_or_default()
            );
        }

        let mut pruned = raw.clone();
        let mut pruning = AssemblyGraphPruningParams::gatk_haplotype_caller_defaults();
        pruning.min_prune_factor = args.min_prune_factor;
        let n_pruned = pruned.apply_pruning(&pruning);
        let (n_ref_p, _) = kmer_presence(&pruned, &ref_ks);
        let (n_alt_p, alt_p_rows) = kmer_presence(&pruned, &alt_ks);
        eprintln!(
            "  PRUNE factor={} edges_removed={n_pruned} nodes={} edges={} REF_K={n_ref_p} ALT_K={n_alt_p} cycles={}",
            args.min_prune_factor,
            pruned.node_count(),
            pruned.edge_count(),
            pruned.has_cycle()
        );
        for (s, sup, out) in &alt_p_rows {
            eprintln!("    PRUNED_ALT {s} support={sup} out_w={out:?}");
        }

        match build_threading_graph_for_haplotype_dump(
            &mid.reference,
            &mid.assembly_reads,
            k,
            args,
            args.allow_low_complexity_graphs,
            args.allow_non_unique_kmers_in_ref,
        ) {
            Ok(Some(g)) => {
                let (n_ref_d, _) = kmer_presence(&g, &ref_ks);
                let (n_alt_d, alt_d) = kmer_presence(&g, &alt_ks);
                eprintln!(
                    "  DUMP_RT nodes={} edges={} REF_K={n_ref_d} ALT_K={n_alt_d} cycles={} merge_haps={}",
                    g.node_count(),
                    g.edge_count(),
                    g.has_cycle(),
                    g.dangling_merge_haps.len()
                );
                for (s, sup, out) in &alt_d {
                    eprintln!("    DUMP_ALT {s} support={sup} out_w={out:?}");
                }
                for (ri, ar) in mid.assembly_reads.iter().enumerate() {
                    let rec = &mid.finalized[ri];
                    let cigar = CigarString(rec.cigar().iter().copied().collect());
                    if let Some(qi) =
                        query_index_at_reference_position(rec.pos(), &cigar, (SITE - 1) as i64)
                    {
                        let starts = overlapping_kmer_starts(ar.bases.len(), qi, k);
                        let present = starts
                            .iter()
                            .filter(|&&st| g.vertex_id_for_kmer(&ar.bases[st..st + k]).is_some())
                            .count();
                        eprintln!(
                            "  DUMP_READ_KMER[{ri}] windows={} in_dump_graph={present}",
                            starts.len()
                        );
                    }
                }
                let n_nodes = g.node_count();
                let n_edges = g.edge_count();
                match find_best_haplotypes_for_assembly(g, args.num_best_haplotypes_per_graph) {
                    Ok((paths, _)) => {
                        let n_alt_paths = paths.iter().filter(|p| !p.is_reference).count();
                        eprintln!(
                            "  KBEST_RT n_paths={} n_alt_flag={} (from nodes={n_nodes} edges={n_edges})",
                            paths.len(),
                            n_alt_paths
                        );
                    }
                    Err(e) => eprintln!("  KBEST_RT err={e}"),
                }
            }
            Ok(None) => eprintln!("  DUMP_RT None (failed gate)"),
            Err(e) => eprintln!("  DUMP_RT err={e}"),
        }

        match build_threading_graph_for_seq_assembly(
            &mid.reference,
            &mid.assembly_reads,
            k,
            args,
            args.allow_low_complexity_graphs,
            args.allow_non_unique_kmers_in_ref,
        ) {
            Ok(Some(g)) => {
                let (n_alt_s, _) = kmer_presence(&g, &alt_ks);
                eprintln!(
                    "  SEQ_RT nodes={} edges={} ALT_K={n_alt_s} cycles={}",
                    g.node_count(),
                    g.edge_count(),
                    g.has_cycle()
                );
                let mut seq = SeqGraph::from_assembly_graph(&g);
                eprintln!(
                    "  SEQ_FROM nodes={} edges={}",
                    seq.node_count(),
                    seq.edge_count()
                );
                seq.clean_non_ref_paths();
                let status = seq.cleanup_seq_graph();
                eprintln!(
                    "  SEQ_CLEAN status={status:?} nodes={} edges={}",
                    seq.node_count(),
                    seq.edge_count()
                );
                if let Ok(paths) =
                    find_best_haplotypes_seq_graph(&seq, args.num_best_haplotypes_per_graph)
                {
                    let n_alt_paths = paths.iter().filter(|p| !p.is_reference).count();
                    eprintln!(
                        "  KBEST_SEQ n_paths={} n_alt_flag={}",
                        paths.len(),
                        n_alt_paths
                    );
                    let mut ref_hap = Haplotype::new(mid.reference.bases.as_slice(), true);
                    let mut cig = crate::cigar::Cigar::new();
                    cig.push(ref_hap.bases.len(), CigarOperator::Match);
                    ref_hap.cigar = Some(cig);
                    let rlen = ref_hap.bases.len();
                    if let Ok(haps) = extract_haplotypes_from_seq_kbest_paths(
                        &paths,
                        &seq,
                        k,
                        &ref_hap,
                        rlen,
                        &args.haplotype_to_reference_sw,
                    ) {
                        hap_summary(
                            &format!("seq_kbest_k{k}"),
                            &haps,
                            &mid.reference.bases,
                            mid.pad_start,
                        );
                    }
                }
            }
            Ok(None) => eprintln!("  SEQ_RT None (cycle / complexity / non-unique abort)"),
            Err(e) => eprintln!("  SEQ_RT err={e}"),
        }

        match extract_rt_haplotypes_before_remove_paths(
            &mid.reference,
            &mid.assembly_reads,
            args,
            k,
            args.allow_low_complexity_graphs,
            args.allow_non_unique_kmers_in_ref,
        ) {
            Ok(haps) => hap_summary(
                &format!("rt_before_remove_k{k}"),
                &haps,
                &mid.reference.bases,
                mid.pad_start,
            ),
            Err(e) => eprintln!("  RT_BEFORE_REMOVE k={k} err={e}"),
        }
    }

    fn run_assemble(label: &str, mid: &MidB, args: &ReadThreadingAssemblerArgs) {
        match assemble_from_ref_and_reads(&mid.reference, &mid.assembly_reads, args) {
            Ok(res) => {
                eprintln!(
                    "ASSEMBLE[{label}] status={:?} k={} n_haps={}",
                    res.status,
                    res.kmer_size,
                    res.haplotypes.len()
                );
                hap_summary(label, &res.haplotypes, &mid.reference.bases, mid.pad_start);
            }
            Err(e) => eprintln!("ASSEMBLE[{label}] err={e}"),
        }
    }

    #[test]
    fn six_r21_mid_b_assembly_read_threading_audit() {
        let Some(mid) = load_mid_b() else {
            eprintln!("Real-data mid-B comparison unavailable");
            return;
        };
        let off = site_offset(mid.pad_start);
        eprintln!("=== 6R.21 INPUT ===");
        eprintln!(
            "REGION contig={} active={}..{} extended={}..{} pad_start={} n_raw={} n_finalized={} n_assembly={}",
            mid.contig,
            mid.region_start,
            mid.region_end,
            mid.ext_start,
            mid.ext_end,
            mid.pad_start,
            mid.n_raw_reads,
            mid.finalized.len(),
            mid.assembly_reads.len()
        );
        let ctx = AssemblyScoringContext {
            padded_reference_start_1based: mid.pad_start,
            active_start_1based: mid.region_start,
            active_end_1based: mid.region_end,
            contig: mid.contig.clone(),
        };
        eprintln!(
            "SCORING p12_cluster={} l_gate={} rt_first_skipped={}",
            ctx.overlaps_p12_cluster(),
            ctx.overlaps_p12_l_gate_interval(),
            ctx.overlaps_p12_cluster() || ctx.overlaps_p12_l_gate_interval()
        );
        eprintln!(
            "REF_LEN={} SITE_OFF={off} REF_BASE={}",
            mid.reference.bases.len(),
            if off < mid.reference.bases.len() {
                mid.reference.bases[off] as char
            } else {
                '?'
            }
        );

        eprintln!("=== RAW BAM (pre-finalize) ===");
        for (i, rec) in mid.raw_reads.iter().enumerate() {
            dump_bam_read(&format!("raw{i}"), rec, mid.pad_start, &mid.reference.bases);
        }
        eprintln!("=== FINALIZED BAM ===");
        for (i, rec) in mid.finalized.iter().enumerate() {
            dump_bam_read(&format!("fin{i}"), rec, mid.pad_start, &mid.reference.bases);
        }
        eprintln!("=== ASSEMBLY READS (seq/qual only) ===");
        for (i, r) in mid.assembly_reads.iter().enumerate() {
            let n_lt10 = r.base_quals.iter().filter(|&&q| q < 10).count();
            eprintln!(
                "AREAD[{i}] len={} q<10={n_lt10} head={} tail={}",
                r.bases.len(),
                ascii(&r.bases[..r.bases.len().min(24)]),
                ascii(&r.bases[r.bases.len().saturating_sub(24)..])
            );
        }

        let args = prod_args(&mid);
        eprintln!("=== CONTROL A — production k=10,25 prune=2 SeqGraph ===");
        audit_k(&mid, 10, &args);
        audit_k(&mid, 25, &args);
        run_assemble("A_prod", &mid, &args);

        eprintln!("=== assemble_reads_with_finalized production ===");
        let mut owned = mid.region.clone();
        let mut ref_cache = ReferenceWindowCache::new(mid.ref_fasta.clone(), 4);
        match assemble_reads_with_finalized(
            &mut owned,
            &mid.dict,
            &mut ref_cache,
            &AssembleReadsArgs::default(),
        ) {
            Ok(assembled) => {
                let set = assembled.assembly;
                eprintln!(
                    "FULL_ASSEMBLE n_haps={} n_events={} has_var={}",
                    set.haplotypes.len(),
                    set.variation_events.len(),
                    set.has_variation_for_calling()
                );
                hap_summary(
                    "full_assemble",
                    &set.haplotypes,
                    &mid.reference.bases,
                    mid.pad_start,
                );
            }
            Err(e) => eprintln!("FULL_ASSEMBLE err={e}"),
        }

        eprintln!("=== CONTROL B — min_prune_factor=1 ===");
        let mut b = args.clone();
        b.min_prune_factor = 1;
        audit_k(&mid, 10, &b);
        audit_k(&mid, 25, &b);
        run_assemble("B_prune1", &mid, &b);

        eprintln!("=== CONTROL B2 — use_seq_graph=false ===");
        let mut b2 = args.clone();
        b2.use_seq_graph = false;
        run_assemble("B2_rt_only", &mid, &b2);

        eprintln!("=== CONTROL C — diagnostic k (harness only) ===");
        for k in [8usize, 11, 15, 21, 35] {
            audit_k(&mid, k, &args);
            let mut ck = args.clone();
            ck.kmer_sizes = vec![k];
            ck.dont_increase_kmer_sizes_for_cycles = true;
            run_assemble(&format!("C_k{k}"), &mid, &ck);
        }

        eprintln!("=== CONTROL D — isolated two assembly reads, k=10 prune=1 ===");
        let mut d = args.clone();
        d.min_prune_factor = 1;
        d.use_seq_graph = false;
        d.kmer_sizes = vec![10];
        d.dont_increase_kmer_sizes_for_cycles = true;
        run_assemble("D_isolate_rt_k10_p1", &mid, &d);

        eprintln!("=== REF non-unique k-mer gate (Java determineNonUniqueKmers) ===");
        for k in [10usize, 25, 35, 45, 55, 65, 75, 85] {
            eprintln!(
                "UNIQ_GATE k={k} has_non_unique={}",
                reference_has_non_unique_kmers(&mid.reference, k)
            );
        }

        eprintln!("=== CONTROL E — last-attempt flags (allow non-unique + low-complexity) ===");
        let mut e = args.clone();
        e.allow_non_unique_kmers_in_ref = true;
        e.allow_low_complexity_graphs = true;
        for k in [10usize, 25, 35, 85] {
            audit_k(&mid, k, &e);
        }
        e.kmer_sizes = vec![85];
        e.dont_increase_kmer_sizes_for_cycles = true;
        run_assemble("E_k85_last_attempt_flags", &mid, &e);

        eprintln!("=== CONTROL G — k=85 last-attempt, dangling_java_exact=false ===");
        let mut gctl = e.clone();
        gctl.kmer_sizes = vec![85];
        gctl.dont_increase_kmer_sizes_for_cycles = true;
        gctl.dangling_java_exact = false;
        audit_k(&mid, 85, &gctl);
        run_assemble("G_k85_dangling_not_exact", &mid, &gctl);

        eprintln!("=== CONTROL F — k=10 SeqGraph/RT with last-attempt flags ===");
        let mut f = e.clone();
        f.kmer_sizes = vec![10];
        f.use_seq_graph = true;
        run_assemble("F_k10_gates_off_seq", &mid, &f);
        f.use_seq_graph = false;
        run_assemble("F_k10_gates_off_rt", &mid, &f);

        if let Ok(res) = assemble_from_ref_and_reads(&mid.reference, &mid.assembly_reads, &args) {
            eprintln!("=== SNP BASE ON A_prod HAPS off={off} ===");
            for (i, h) in res.haplotypes.iter().enumerate() {
                let b = if off < h.bases.len() {
                    h.bases[off] as char
                } else {
                    '?'
                };
                let ctx_lo = off.saturating_sub(8);
                let ctx_hi = (off + 9).min(h.bases.len());
                eprintln!(
                    "  HAP[{i}] ref={} snp_base={b} ctx={}",
                    h.is_reference,
                    ascii(&h.bases[ctx_lo..ctx_hi])
                );
            }
        }

        assert_eq!(mid.n_raw_reads, 2, "mid-B fixture is the 2-read region");
        assert_eq!(
            mid.reference.bases[off], b'C',
            "reference at 92317399 is C (oracle REF)"
        );
    }
}
