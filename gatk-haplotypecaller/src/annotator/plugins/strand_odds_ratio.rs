//! GATK `StrandOddsRatio` (production).

const PSEUDOCOUNT: f64 = 1.0;

/// SOR from ref/alt × forward/reverse counts (GATK `StrandOddsRatio.calculateSOR`).
pub fn strand_odds_ratio(ref_fw: u32, ref_rv: u32, alt_fw: u32, alt_rv: u32) -> f64 {
    let table = [[ref_fw, ref_rv], [alt_fw, alt_rv]];
    calculate_sor(&table)
}

fn calculate_sor(table: &[[u32; 2]; 2]) -> f64 {
    let t00 = table[0][0] as f64 + PSEUDOCOUNT;
    let t01 = table[0][1] as f64 + PSEUDOCOUNT;
    let t10 = table[1][0] as f64 + PSEUDOCOUNT;
    let t11 = table[1][1] as f64 + PSEUDOCOUNT;
    let ratio = (t00 / t01) * (t11 / t10) + (t01 / t00) * (t10 / t11);
    let ref_ratio = t00.min(t01) / t00.max(t01);
    let alt_ratio = t10.min(t11) / t10.max(t11);
    (ratio).ln() + ref_ratio.ln() - alt_ratio.ln()
}
