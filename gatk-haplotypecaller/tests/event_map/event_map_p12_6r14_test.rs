//! 6R.14 TEST-ONLY: control haplotype EventMap / trim / allele-mapper vs GATK 4.4 source walk.
//! Does not call W-H1 constructors. Does not change production algorithms.

#[cfg(test)]
mod traces {
    use crate::assembly_region_iterator::AssemblyRegion;
    use crate::assembly_result_set::{AssemblyResultSet, DEFAULT_MAX_MNP_DISTANCE};
    use crate::cigar::{Cigar, CigarOperator};
    use crate::event_map::{
        prefer_indel_over_colocated_snps, variation_events_for_haplotype, EventMap, VariationEvent,
    };
    use crate::feature_context::FeatureContext;
    use crate::genome_loc::{GenomeLoc, GenomePosition};
    use crate::haplotype::Haplotype;
    use crate::hc_allele_mapping::create_allele_mapper;
    use crate::read_event_discovery::{
        reference_motif_cluster_coupled_events, P12_CLUSTER_ATG_START, P12_CLUSTER_TTC_START,
        SUPPLEMENT_HAPLOTYPE_SCORE,
    };
    use crate::read_threading_assembler::AssemblyStatus;
    use crate::reference_context::ReferenceContext;
    use std::collections::BTreeMap;
    use std::path::Path;

    const ALT_WIN: &[u8] = b"CTTTTATGTGATGTAT";
    const REAL_P12_ACTIVE_START: u64 = P12_CLUSTER_TTC_START - 96;
    const REAL_P12_ACTIVE_END: u64 = P12_CLUSTER_TTC_START + 76;

    #[derive(Clone, Debug, PartialEq, Eq)]
    struct JavaVc {
        start: u64,
        end: u64,
        ref_al: String,
        alt_al: String,
    }

    fn load_real_p12_ref() -> Option<(Vec<u8>, u64)> {
        use crate::assembly_region_finalize::{
            assembly_reference_read, finalize_region_reads_for_assembly,
            gatk_min_tail_quality_for_assembly, records_to_assembly_reads,
        };
        use crate::read_model::ReadFilterParams;
        use crate::walker_apply::{call_disposition, AssemblyRegionCallDisposition};
        use crate::walker_traversal::{
            flatten_assembly_regions, traverse_assembly_region_walker, WalkerTraversalConfig,
        };
        use gatk_core::reference::{
            parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary,
        };

        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
        let ref_path = root.join("parity/realworld/assets/hs37d5.simple.fa");
        let bam = root.join("parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam");
        if !ref_path.is_file() || !bam.is_file() {
            return None;
        }
        let dict = SequenceDictionary::from_fasta_path(&ref_path).ok()?;
        let interval = format!("2:{REAL_P12_ACTIVE_START}-{REAL_P12_ACTIVE_END}");
        let specs = parse_intervals_cli_string(&dict, &interval).ok()?;
        let walk = traverse_assembly_region_walker(
            &dict,
            &specs,
            &ref_path,
            &bam,
            &ReadFilterParams::gatk_standard_hc(),
            &WalkerTraversalConfig::gatk_haplotype_caller_production(100),
        )
        .ok()?;
        let regions = flatten_assembly_regions(&walk);
        let region = regions.iter().find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= P12_CLUSTER_TTC_START
                && r.end.get() >= P12_CLUSTER_ATG_START
        })?;
        let mut ref_cache = ReferenceWindowCache::new(ref_path.clone(), 4);
        let reference = assembly_reference_read(&dict, &mut ref_cache, region).ok()?;
        let _finalized = finalize_region_reads_for_assembly(
            &region.reads,
            region,
            true,
            gatk_min_tail_quality_for_assembly(10),
            false,
        );
        let _reads = records_to_assembly_reads(&_finalized);
        let pad = region
            .extended_start
            .get()
            .saturating_sub(crate::assembly_region_finalize::GATK_REFERENCE_PADDING_FOR_ASSEMBLY)
            .max(1);
        Some((reference.bases, pad))
    }

    /// Control bytes matching the W-H1 product (REF splice), not a production call.
    fn control_alt_bases(ref_bytes: &[u8], pad: u64) -> Option<Vec<u8>> {
        let ttc_off = P12_CLUSTER_TTC_START.saturating_sub(pad) as usize;
        let atg_off = P12_CLUSTER_ATG_START.saturating_sub(pad) as usize;
        if ttc_off + 3 > ref_bytes.len() || atg_off >= ref_bytes.len() {
            return None;
        }
        if &ref_bytes[ttc_off..ttc_off + 3] != b"TTC" {
            return None;
        }
        if ref_bytes[atg_off] != b'A' && ref_bytes[atg_off] != b'a' {
            return None;
        }
        let mut out = ref_bytes.to_vec();
        out.remove(ttc_off + 1);
        out.remove(ttc_off + 1);
        let atg_adj = atg_off.saturating_sub(2);
        if !out
            .get(atg_adj)
            .copied()
            .unwrap_or(0)
            .eq_ignore_ascii_case(&b'A')
        {
            return None;
        }
        out.insert(atg_adj + 1, b'T');
        out.insert(atg_adj + 2, b'G');
        Some(out)
    }

    fn forced_cigar(pad: u64, ref_len: usize) -> Cigar {
        let ttc_off = P12_CLUSTER_TTC_START.saturating_sub(pad) as usize;
        let tail = ref_len.saturating_sub(ttc_off + 3);
        let mut c = Cigar::new();
        if ttc_off > 0 {
            c.push(ttc_off, CigarOperator::Match);
        }
        c.push(2, CigarOperator::Deletion);
        c.push(1, CigarOperator::Match);
        c.push(2, CigarOperator::Insertion);
        if tail > 0 {
            c.push(tail, CigarOperator::Match);
        }
        c
    }

    fn control_haplotype(ref_bytes: &[u8], pad: u64) -> Option<Haplotype> {
        let bases = control_alt_bases(ref_bytes, pad)?;
        let cigar = forced_cigar(pad, ref_bytes.len());
        let mut h = Haplotype::new(bases, false);
        h.cigar = Some(cigar);
        h.score = SUPPLEMENT_HAPLOTYPE_SCORE;
        h.alignment_start_hap_wrt_ref = 0;
        h.genome_loc = Some(GenomeLoc::new(
            pad,
            pad.saturating_add(ref_bytes.len() as u64).saturating_sub(1),
        ));
        Some(h)
    }

    fn is_regular(b: u8) -> bool {
        matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't')
    }

    /// GATK 4.4 `EventMap.processCigarForInitialEvents` before `addVC` merge.
    fn java_source_eventmap_proposed(hap: &Haplotype, ref_bytes: &[u8], pad: u64) -> Vec<JavaVc> {
        let cigar = hap.cigar.as_ref().expect("control cigar");
        let alignment = hap.bases.as_slice();
        let mut ref_pos = hap.alignment_start_hap_wrt_ref;
        let mut alignment_pos = 0usize;
        let mut proposed = Vec::new();
        for (cigar_index, el) in cigar.elements.iter().enumerate() {
            let n = el.length;
            match el.operator {
                CigarOperator::Insertion => {
                    if ref_pos > 0 && ref_pos <= ref_bytes.len() {
                        let insertion_start = pad + ref_pos as u64 - 1;
                        let ref_byte = ref_bytes[ref_pos - 1];
                        let is_edge = cigar_index == 0 || cigar_index == cigar.elements.len() - 1;
                        if is_regular(ref_byte) && !is_edge {
                            let mut ins = vec![ref_byte];
                            let end = alignment_pos.saturating_add(n).min(alignment.len());
                            ins.extend_from_slice(&alignment[alignment_pos..end]);
                            if ins.len() >= 2 && ins.iter().copied().all(is_regular) {
                                proposed.push(JavaVc {
                                    start: insertion_start,
                                    end: insertion_start,
                                    ref_al: String::from_utf8(vec![ref_byte]).unwrap(),
                                    alt_al: String::from_utf8(ins).unwrap(),
                                });
                            }
                        }
                    }
                    alignment_pos += n;
                }
                CigarOperator::Deletion => {
                    if ref_pos > 0 && ref_pos + n <= ref_bytes.len() {
                        let mut del = vec![ref_bytes[ref_pos - 1]];
                        del.extend_from_slice(&ref_bytes[ref_pos..ref_pos + n]);
                        let deletion_start = pad + ref_pos as u64 - 1;
                        if is_regular(ref_bytes[ref_pos - 1]) && del.iter().copied().all(is_regular)
                        {
                            proposed.push(JavaVc {
                                start: deletion_start,
                                end: deletion_start + n as u64,
                                ref_al: String::from_utf8(del).unwrap(),
                                alt_al: String::from_utf8(vec![ref_bytes[ref_pos - 1]]).unwrap(),
                            });
                        }
                    }
                    ref_pos += n;
                }
                CigarOperator::Match => {
                    let mut mismatches = Vec::new();
                    for offset in 0..n {
                        if ref_pos + offset >= ref_bytes.len() {
                            break;
                        }
                        let rb = ref_bytes[ref_pos + offset];
                        let ab = alignment
                            .get(alignment_pos + offset)
                            .copied()
                            .unwrap_or(b'N');
                        if rb != ab && is_regular(rb) && is_regular(ab) {
                            mismatches.push(offset);
                        }
                    }
                    let mut i = 0;
                    while i < mismatches.len() {
                        let start = mismatches[i];
                        let mut end = start;
                        i += 1;
                        while i < mismatches.len() && mismatches[i] - end <= 0 {
                            end = mismatches[i];
                            i += 1;
                        }
                        proposed.push(JavaVc {
                            start: pad + ref_pos as u64 + start as u64,
                            end: pad + ref_pos as u64 + end as u64,
                            ref_al: String::from_utf8(
                                ref_bytes[ref_pos + start..=ref_pos + end].to_vec(),
                            )
                            .unwrap(),
                            alt_al: String::from_utf8(
                                alignment[alignment_pos + start..=alignment_pos + end].to_vec(),
                            )
                            .unwrap(),
                        });
                    }
                    ref_pos += n;
                    alignment_pos += n;
                }
                _ => {
                    ref_pos += n;
                    alignment_pos += n;
                }
            }
        }
        proposed
    }

    /// GATK 4.4 `addVC(merge=true)` / `makeBlock` over proposed events.
    fn java_source_eventmap_merged(proposed: &[JavaVc]) -> Vec<JavaVc> {
        let mut by_start: BTreeMap<u64, JavaVc> = BTreeMap::new();
        for vc in proposed {
            if let Some(prev) = by_start.remove(&vc.start) {
                by_start.insert(vc.start, java_make_block(prev, vc.clone()));
            } else {
                by_start.insert(vc.start, vc.clone());
            }
        }
        by_start.into_values().collect()
    }

    /// GATK `EventMap.makeBlock` (SNP+indel or ins+del at the same start).
    fn java_make_block(vc1: JavaVc, vc2: JavaVc) -> JavaVc {
        let vc1_snp = vc1.ref_al.len() == 1 && vc1.alt_al.len() == 1;
        if vc1_snp {
            if vc1.ref_al == vc2.ref_al {
                JavaVc {
                    start: vc1.start,
                    end: vc1.end,
                    ref_al: vc1.ref_al,
                    alt_al: format!("{}{}", vc1.alt_al, &vc2.alt_al[1..]),
                }
            } else {
                JavaVc {
                    start: vc1.start,
                    end: vc2.end,
                    ref_al: vc2.ref_al,
                    alt_al: vc1.alt_al,
                }
            }
        } else {
            let (ins, del) = if vc1.alt_al.len() > vc1.ref_al.len() {
                (vc1, vc2)
            } else {
                (vc2, vc1)
            };
            JavaVc {
                start: del.start,
                end: del.end,
                ref_al: del.ref_al,
                alt_al: ins.alt_al,
            }
        }
    }

    fn dump_events(label: &str, events: &[VariationEvent]) {
        eprintln!("EVENTS {label} n={}", events.len());
        for e in events {
            eprintln!(
                "  {} {} REF={} ALT={} indel={}",
                e.start_1based.get(),
                e.end_1based.get(),
                e.ref_allele,
                e.alt_allele,
                e.is_indel()
            );
        }
    }

    fn dummy_region(start: u64, end: u64, ext_start: u64, ext_end: u64) -> AssemblyRegion {
        AssemblyRegion {
            contig: "2".into(),
            start: GenomePosition::new_1based(start),
            end: GenomePosition::new_1based(end),
            is_active: true,
            extended_start: GenomePosition::new_1based(ext_start),
            extended_end: GenomePosition::new_1based(ext_end),
            extension: 100,
            reads: Vec::new(),
            read_qnames: Vec::new(),
            reference: ReferenceContext::empty(),
            features: FeatureContext::empty(),
            pileup_loci: Vec::new(),
        }
    }

    #[test]
    fn six_r14_control_haplotype_eventmap_trim_mapper() {
        let Some((ref_bytes, pad)) = load_real_p12_ref() else {
            eprintln!("Real-data P12 comparison unavailable");
            return;
        };
        let Some(alt) = control_haplotype(&ref_bytes, pad) else {
            panic!("control haplotype construction failed");
        };
        let mut ref_hap = Haplotype::new(ref_bytes.as_slice(), true);
        let mut rc = Cigar::new();
        rc.push(ref_bytes.len(), CigarOperator::Match);
        ref_hap.cigar = Some(rc);
        ref_hap.genome_loc = alt.genome_loc;

        eprintln!("=== 6R.14 CONTROL HAPLOTYPE ===");
        eprintln!(
            "pad={pad} ref_len={} hap_len={} cigar={:?} alt_win={} score={}",
            ref_bytes.len(),
            alt.bases.len(),
            alt.cigar.as_ref().map(|c| c.to_gatk_string()),
            alt.bases.windows(ALT_WIN.len()).any(|w| w == ALT_WIN),
            alt.score
        );
        let ttc_off = (P12_CLUSTER_TTC_START - pad) as usize;
        eprintln!(
            "REF[ttc-1..ttc+4]={:?} HAP[ttc-1..ttc+4]={:?}",
            ref_bytes.get(ttc_off.saturating_sub(1)..ttc_off.saturating_add(5)),
            alt.bases
                .get(ttc_off.saturating_sub(1)..ttc_off.saturating_add(5)),
        );
        eprintln!(
            "1M compare hap[{ttc_off}]={:?} vs ref[{}]={:?}",
            alt.bases.get(ttc_off).map(|b| *b as char),
            ttc_off + 2,
            ref_bytes.get(ttc_off + 2).map(|b| *b as char),
        );

        let java_proposed = java_source_eventmap_proposed(&alt, &ref_bytes, pad);
        eprintln!("JAVA_SOURCE_WALK_UNMERGED n={}", java_proposed.len());
        for v in &java_proposed {
            eprintln!(
                "  java_unmerged {}-{} REF={} ALT={}",
                v.start, v.end, v.ref_al, v.alt_al
            );
        }
        let java_merged = java_source_eventmap_merged(&java_proposed);
        eprintln!("JAVA_SOURCE_WALK_MERGED n={}", java_merged.len());
        for v in &java_merged {
            eprintln!(
                "  java_merged {}-{} REF={} ALT={}",
                v.start, v.end, v.ref_al, v.alt_al
            );
        }
        let motif = reference_motif_cluster_coupled_events(&ref_bytes, pad, "2");
        dump_events("rust_reference_motif_cluster_coupled_events", &motif);

        let raw_map = EventMap::from_haplotype_and_reference(&alt, &ref_hap, &ref_bytes, pad, 0);
        eprintln!("RUST_EVENTMAP_KEYS_PAD_OFFSET n={}", raw_map.events.len());
        for e in &raw_map.events {
            eprintln!(
                "  pad_off={} REF={} ALT={}",
                e.start.get(),
                String::from_utf8_lossy(&e.ref_bases),
                String::from_utf8_lossy(&e.alt_bases)
            );
        }
        let raw_events = raw_map.variation_events("2", pad);
        dump_events("rust_from_haplotype_raw", &raw_events);

        let hap_events = variation_events_for_haplotype(&alt, &ref_hap, &ref_bytes, pad, 0, "2");
        dump_events("rust_variation_events_for_haplotype", &hap_events);

        let mut after_prefer = raw_events.clone();
        prefer_indel_over_colocated_snps(&mut after_prefer);
        dump_events("rust_raw_after_prefer_indel", &after_prefer);

        let rust_starts: Vec<u64> = hap_events.iter().map(|e| e.start_1based.get()).collect();
        let java_starts: Vec<u64> = java_merged.iter().map(|v| v.start).collect();
        eprintln!("STARTS rust={rust_starts:?} java_source={java_starts:?}");

        let same_start_merge = java_merged.iter().any(|j| {
            raw_events
                .iter()
                .filter(|e| e.start_1based.get() == j.start)
                .count()
                > 1
        });
        eprintln!("JAVA_MAKEBLOCK_NEEDED_SAME_START={same_start_merge}");

        let trim_span = GenomeLoc::new(REAL_P12_ACTIVE_START, REAL_P12_ACTIVE_END);
        let trimmed = alt
            .trim(&trim_span, true)
            .expect("control hap must trim into active span");
        eprintln!(
            "TRIM hap_len {}→{} cigar={:?} align_start={}",
            alt.bases.len(),
            trimmed.bases.len(),
            trimmed.cigar.as_ref().map(|c| c.to_gatk_string()),
            trimmed.alignment_start_hap_wrt_ref
        );
        let trim_events =
            variation_events_for_haplotype(&trimmed, &ref_hap, &ref_bytes, pad, 0, "2");
        dump_events("rust_trimmed_variation_events", &trim_events);

        let assembly = AssemblyResultSet::from_assembly_for_calling_owned(
            AssemblyStatus::AssembledSomeVariation,
            85,
            vec![ref_hap.clone(), alt.clone()],
            ref_bytes.clone(),
            pad,
            "2",
            DEFAULT_MAX_MNP_DISTANCE,
        );
        let region = dummy_region(
            REAL_P12_ACTIVE_START,
            REAL_P12_ACTIVE_END,
            REAL_P12_ACTIVE_START.saturating_sub(100),
            REAL_P12_ACTIVE_END.saturating_add(100),
        );
        let trimmed_set = assembly.trim_to(&region).expect("trim_to");
        eprintln!(
            "TRIM_TO n_haps={} events={}",
            trimmed_set.haplotypes.len(),
            trimmed_set.variation_events().len()
        );
        dump_events("trim_to_variation_events", trimmed_set.variation_events());

        let haps = [ref_hap.clone(), alt.clone()];
        eprintln!("=== HAPLOTYPE-TO-ALLELE MAPPER ===");
        for (label, events) in [
            ("eventmap_hap_events", hap_events.as_slice()),
            ("motif_coupled", motif.as_slice()),
        ] {
            eprintln!("MAPPER_GROUP {label}");
            for e in events {
                let mapper =
                    create_allele_mapper(e, e.start_1based.get(), &haps, pad, &ref_bytes, 0, true);
                eprintln!(
                    "  loc={} {}:{}  ref_haps={:?} alt_haps={:?}",
                    e.start_1based.get(),
                    e.ref_allele,
                    e.alt_allele,
                    mapper.ref_haplotype_indices,
                    mapper.alt_haplotype_indices
                );
            }
        }

        let rust_raw_matches_java_merged = java_merged.len() == raw_events.len()
            && java_merged.iter().zip(raw_events.iter()).all(|(j, r)| {
                j.start == r.start_1based.get()
                    && j.end == r.end_1based.get()
                    && j.ref_al == r.ref_allele
                    && j.alt_al == r.alt_allele
            });
        eprintln!(
            "JAVA_UNMERGED_N={} (CIGAR walk replica; production EventMap is post-addVC)",
            java_proposed.len()
        );
        eprintln!("RUST_RAW_MATCHES_JAVA_MERGED={rust_raw_matches_java_merged}");
        eprintln!(
            "RUST_HAP_EVENTS_MATCH_JAVA_MERGED={}",
            hap_events.len() == java_merged.len()
                && hap_events.iter().zip(java_merged.iter()).all(|(r, j)| {
                    r.start_1based.get() == j.start
                        && r.ref_allele == j.ref_al
                        && r.alt_allele == j.alt_al
                        && r.end_1based.get() == j.end
                })
        );

        assert!(
            alt.bases.windows(ALT_WIN.len()).any(|w| w == ALT_WIN),
            "control hap must contain canonical ALT_WIN"
        );
        assert_eq!(
            alt.cigar.as_ref().map(|c| c.to_gatk_string()).as_deref(),
            Some("696M2D1M2I674M")
        );
        assert!(
            rust_raw_matches_java_merged,
            "production EventMap must match Java 4.4 addVC(merge=true) / makeBlock"
        );
        assert_eq!(
            java_merged
                .iter()
                .map(|v| (v.start, v.ref_al.as_str(), v.alt_al.as_str()))
                .collect::<Vec<_>>(),
            vec![
                (P12_CLUSTER_TTC_START - 1, "TTT", "T"),
                (P12_CLUSTER_TTC_START + 2, "C", "TAT"),
                (P12_CLUSTER_ATG_START, "A", "G"),
            ]
        );
        assert_eq!(
            hap_events
                .iter()
                .map(|e| (
                    e.start_1based.get(),
                    e.end_1based.get(),
                    e.ref_allele.as_str(),
                    e.alt_allele.as_str()
                ))
                .collect::<Vec<_>>(),
            vec![
                (
                    P12_CLUSTER_TTC_START - 1,
                    P12_CLUSTER_TTC_START + 1,
                    "TTT",
                    "T"
                ),
                (
                    P12_CLUSTER_TTC_START + 2,
                    P12_CLUSTER_TTC_START + 2,
                    "C",
                    "TAT"
                ),
                (P12_CLUSTER_ATG_START, P12_CLUSTER_ATG_START, "A", "G"),
            ]
        );
        assert!(
            !hap_events
                .iter()
                .any(|e| e.ref_allele == "TTC" && e.alt_allele == "T"),
            "EventMap must not emit biological TTC→T; that allele is motif-injected"
        );
        assert!(
            !hap_events
                .iter()
                .any(|e| e.ref_allele == "A" && e.alt_allele == "ATG"),
            "EventMap must not emit biological A→ATG; that allele is motif-injected"
        );
        assert_eq!(motif.len(), 2);
        assert_eq!(motif[0].ref_allele, "TTC");
        assert_eq!(motif[1].alt_allele, "ATG");
        let _ = rust_starts;
        let _ = java_starts;
        let _ = same_start_merge;
    }
}
