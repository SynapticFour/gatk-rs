package org.broadinstitute.hellbender.tools.walkers.genotyper;

import org.broadinstitute.hellbender.utils.MathUtils;

import java.util.List;

/** EM iteration counter aligned with Rust {@code calculate_biallelic_af_em} (g2-af parity). */
final class HcParityAfEmIterations {

    private HcParityAfEmIterations() {}

    static final class EmMetrics {
        final int altAc;
        final double af;
        final double log10PNoVariant;
        final int iterations;

        EmMetrics(int altAc, double af, double log10PNoVariant, int iterations) {
            this.altAc = altAc;
            this.af = af;
            this.log10PNoVariant = log10PNoVariant;
            this.iterations = iterations;
        }
    }

    static EmMetrics computeBiallelicMetrics(
            final double[] log10Gl, final GenotypeCalculationArgumentCollection genotypeArgs) {
        final double threshold = 0.1;
        final double refPseudo =
                genotypeArgs.snpHeterozygosity
                        / Math.pow(genotypeArgs.heterozygosityStandardDeviation, 2);
        final double snpPseudo = genotypeArgs.snpHeterozygosity * refPseudo;
        final double flat = -Math.log10(2.0);
        double[] log10Af = new double[] {flat, flat};
        double[] alleleCounts = new double[] {0.0, 0.0};
        int iterations = 0;
        while (true) {
            iterations++;
            final double[] newCounts = effectiveAlleleCounts(List.of(log10Gl), log10Af);
            double maxDiff = 0.0;
            for (int i = 0; i < 2; i++) {
                maxDiff = Math.max(maxDiff, Math.abs(newCounts[i] - alleleCounts[i]));
            }
            alleleCounts = newCounts;
            if (maxDiff <= threshold || iterations > 100) {
                break;
            }
            final double[] posteriorPseudo =
                    new double[] {
                        refPseudo + alleleCounts[0], snpPseudo + alleleCounts[1]
                    };
            log10Af = log10DirichletMeanWeights(posteriorPseudo);
        }
        double log10PNoVariant = 0.0;
        for (final double[] gl : List.of(log10Gl)) {
            if (gl.length >= 3) {
                log10PNoVariant += log10NormalizedGenotypePosteriors(gl, log10Af)[0];
            }
        }
        log10PNoVariant = Math.min(log10PNoVariant, 0.0);
        final double total = alleleCounts[0] + alleleCounts[1];
        final double af = total > 0.0 ? alleleCounts[1] / total : 0.0;
        return new EmMetrics(
                (int) Math.round(alleleCounts[1]), af, log10PNoVariant, iterations);
    }

    static int countBiallelicEmIterations(
            final double[] log10Gl, final GenotypeCalculationArgumentCollection genotypeArgs) {
        return computeBiallelicMetrics(log10Gl, genotypeArgs).iterations;
    }

    private static double[] log10DirichletMeanWeights(final double[] pseudocounts) {
        double sum = 0.0;
        for (final double c : pseudocounts) {
            sum += c;
        }
        final double[] out = new double[pseudocounts.length];
        for (int i = 0; i < pseudocounts.length; i++) {
            out[i] = Math.log10(Math.max(pseudocounts[i] / sum, 1e-300));
        }
        return out;
    }

    private static double[] effectiveAlleleCounts(
            final List<double[]> samplesGl, final double[] log10Af) {
        double log10Ref = Double.NEGATIVE_INFINITY;
        double log10Alt = Double.NEGATIVE_INFINITY;
        for (final double[] gl : samplesGl) {
            if (gl.length < 3) {
                continue;
            }
            final double[] post = log10NormalizedGenotypePosteriors(gl, log10Af);
            for (int gi = 0; gi < 3; gi++) {
                final int[] ac = diploidAlleleCounts(gi);
                if (ac[0] > 0) {
                    log10Ref = log10SumPair(log10Ref, post[gi] + Math.log10(ac[0]));
                }
                if (ac[1] > 0) {
                    log10Alt = log10SumPair(log10Alt, post[gi] + Math.log10(ac[1]));
                }
            }
        }
        return new double[] {Math.pow(10.0, log10Ref), Math.pow(10.0, log10Alt)};
    }

    private static int[] diploidAlleleCounts(final int genotypeIndex) {
        switch (genotypeIndex) {
            case 0:
                return new int[] {2, 0};
            case 1:
                return new int[] {1, 1};
            default:
                return new int[] {0, 2};
        }
    }

    private static double[] log10NormalizedGenotypePosteriors(
            final double[] gl, final double[] log10Af) {
        final double[] log10Post =
                new double[] {
                    Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY
                };
        for (int gi = 0; gi < Math.min(3, gl.length); gi++) {
            final int[] ac = diploidAlleleCounts(gi);
            final double logCombo = MathUtils.log10BinomialCoefficient(2, ac[1]);
            log10Post[gi] = gl[gi] + logCombo + ac[0] * log10Af[0] + ac[1] * log10Af[1];
        }
        return MathUtils.normalizeLog10(log10Post);
    }

    private static double log10SumPair(final double a, final double b) {
        return MathUtils.log10SumLog10(new double[] {a, b});
    }
}
