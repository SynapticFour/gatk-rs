//! GATK `DRAGENMappingQualityReadTransformer` (/ HC `--transform-dragen-mapping-quality`).

use rust_htslib::bam;

/// SAM optional tag for DRAGEN extended mapping quality (`XQ`).
pub const EXTENDED_MAPPING_QUALITY_TAG: &[u8] = b"XQ";

const MQ_TABLE_X: [i32; 6] = [0, 30, 60, 100, 200, 256];
const MQ_TABLE_Y: [i32; 6] = [0, 30, 40, 45, 50, 50];

/// GATK `DRAGENMappingQualityReadTransformer.mapMappingQualityToPhredLikelihoodScore`.
pub fn map_dragen_mq_to_phred(val: i32) -> u8 {
    for i in 1..MQ_TABLE_X.len() {
        if val <= MQ_TABLE_X[i] {
            let xfactor =
                (val - MQ_TABLE_X[i - 1]) as f64 / (MQ_TABLE_X[i] - MQ_TABLE_X[i - 1]) as f64;
            let score =
                MQ_TABLE_Y[i - 1] as f64 + xfactor * (MQ_TABLE_Y[i] - MQ_TABLE_Y[i - 1]) as f64;
            return score.round() as u8;
        }
    }
    50
}

/// If `XQ` is present, copy extended MQ into the record's MAPQ field (Java transformer).
pub fn apply_dragen_mapping_quality_transform(records: &mut [bam::Record]) {
    for rec in records.iter_mut() {
        let xq = match rec.aux(EXTENDED_MAPPING_QUALITY_TAG) {
            Ok(rust_htslib::bam::record::Aux::I32(v)) => v,
            Ok(rust_htslib::bam::record::Aux::U8(v)) => i32::from(v),
            Ok(rust_htslib::bam::record::Aux::U16(v)) => i32::from(v),
            Ok(rust_htslib::bam::record::Aux::U32(v)) => v as i32,
            _ => continue,
        };
        rec.set_mapq(map_dragen_mq_to_phred(xq));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn table_interpolation_matches_java_knots() {
        assert_eq!(map_dragen_mq_to_phred(0), 0);
        assert_eq!(map_dragen_mq_to_phred(30), 30);
        assert_eq!(map_dragen_mq_to_phred(60), 40);
        assert_eq!(map_dragen_mq_to_phred(100), 45);
        assert_eq!(map_dragen_mq_to_phred(200), 50);
        assert_eq!(map_dragen_mq_to_phred(250), 50);
    }
}
