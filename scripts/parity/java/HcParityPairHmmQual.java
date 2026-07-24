import org.broadinstitute.hellbender.utils.QualityUtils;

/** GATK `PairHMMLikelihoodCalculationEngine.capMinimumReadQualities` (base qualities only). */
public final class HcParityPairHmmQual {

    private HcParityPairHmmQual() {}

    public static byte[] capBaseQualities(
            final byte[] quals,
            final int mapq,
            final byte baseQualityScoreThreshold,
            final boolean disableCapReadQualitiesToMapq) {
        final byte[] out = quals.clone();
        for (int i = 0; i < out.length; i++) {
            if (!disableCapReadQualitiesToMapq) {
                out[i] = (byte) Math.min(0xff & out[i], mapq);
            }
            out[i] =
                    setToFixedValueIfTooLow(
                            out[i], baseQualityScoreThreshold, QualityUtils.MIN_USABLE_Q_SCORE);
        }
        return out;
    }

    private static byte setToFixedValueIfTooLow(
            final byte currentVal, final byte minQual, final byte fixedQual) {
        return currentVal < minQual ? fixedQual : currentVal;
    }
}
