//! Mathematical utilities for genomic analysis

use crate::types::*;
use std::f64::consts::PI;

/// Natural logarithm of 10
const LN_10: f64 = std::f64::consts::LN_10;

/// Log-likelihood calculations
pub mod likelihood {
    use super::*;

    /// Calculate log-likelihood ratio
    pub fn log_likelihood_ratio(log_likelihood1: f64, log_likelihood2: f64) -> f64 {
        log_likelihood1 - log_likelihood2
    }

    /// Convert log probability to probability
    pub fn log_to_prob(log_p: f64) -> f64 {
        log_p.exp()
    }

    /// Convert probability to log probability (natural log)
    pub fn prob_to_log(p: f64) -> f64 {
        p.ln()
    }

    /// Calculate log10 probability to natural log probability
    pub fn log10_to_ln(log10_p: f64) -> f64 {
        log10_p * LN_10
    }

    /// Calculate natural log probability to log10 probability
    pub fn ln_to_log10(ln_p: f64) -> f64 {
        ln_p / LN_10
    }

    /// Safe log addition (log(a + b) from log(a) and log(b))
    pub fn log_add(log_a: f64, log_b: f64) -> f64 {
        if log_a.is_infinite() && log_a < 0.0 {
            return log_b;
        }
        if log_b.is_infinite() && log_b < 0.0 {
            return log_a;
        }

        let max_log = log_a.max(log_b);
        let min_log = log_a.min(log_b);

        max_log + (1.0 + (min_log - max_log).exp()).ln()
    }

    /// Safe log subtraction (log(a - b) from log(a) and log(b), assuming a > b)
    pub fn log_subtract(log_a: f64, log_b: f64) -> f64 {
        if log_a <= log_b {
            return f64::NEG_INFINITY;
        }

        log_a + (1.0 - (log_b - log_a).exp()).ln()
    }
}

/// Phred score calculations
pub mod phred {
    /// Convert Phred quality score to error probability
    pub fn phred_to_error_probability(phred: u8) -> f64 {
        10.0f64.powf(-(phred as f64) / 10.0)
    }

    /// Convert error probability to Phred quality score
    pub fn error_probability_to_phred(error_prob: f64) -> u8 {
        if error_prob <= 0.0 {
            return 93; // Maximum Phred score
        }
        if error_prob >= 1.0 {
            return 0;
        }

        let phred = (-10.0 * error_prob.log10()).clamp(0.0, 93.0);
        phred as u8
    }

    /// Convert Phred score to log10 probability
    pub fn phred_to_log10_prob(phred: u8) -> f64 {
        -(phred as f64) / 10.0
    }

    /// Convert log10 probability to Phred score
    pub fn log10_prob_to_phred(log10_p: f64) -> u8 {
        (-log10_p * 10.0).clamp(0.0, 93.0) as u8
    }
}

/// Statistical calculations
pub mod stats {
    use super::*;

    /// Calculate binomial probability
    pub fn binomial_probability(n: u32, k: u32, p: f64) -> f64 {
        if k > n || !(0.0..=1.0).contains(&p) {
            return 0.0;
        }

        let ln_p = p.ln();
        let ln_1_minus_p = (1.0 - p).ln();

        let mut ln_prob = ln_factorial(n) - ln_factorial(k) - ln_factorial(n - k);
        ln_prob += (k as f64) * ln_p + ((n - k) as f64) * ln_1_minus_p;

        ln_prob.exp()
    }

    /// Calculate binomial log probability
    pub fn binomial_log_probability(n: u32, k: u32, p: f64) -> f64 {
        if k > n || p <= 0.0 || p >= 1.0 {
            return f64::NEG_INFINITY;
        }

        let ln_p = p.ln();
        let ln_1_minus_p = (1.0 - p).ln();

        ln_factorial(n) - ln_factorial(k) - ln_factorial(n - k)
            + (k as f64) * ln_p
            + ((n - k) as f64) * ln_1_minus_p
    }

    /// Calculate Poisson probability
    pub fn poisson_probability(k: u32, lambda: f64) -> f64 {
        if lambda < 0.0 {
            return 0.0;
        }

        let ln_prob = (k as f64) * lambda.ln() - lambda - ln_factorial(k);
        ln_prob.exp()
    }

    /// Calculate Poisson log probability
    pub fn poisson_log_probability(k: u32, lambda: f64) -> f64 {
        if lambda < 0.0 {
            return f64::NEG_INFINITY;
        }

        (k as f64) * lambda.ln() - lambda - ln_factorial(k)
    }

    /// Calculate normal probability density function
    pub fn normal_pdf(x: f64, mean: f64, std_dev: f64) -> f64 {
        if std_dev <= 0.0 {
            return 0.0;
        }

        let variance = std_dev * std_dev;
        let exponent = -((x - mean).powi(2)) / (2.0 * variance);
        (1.0 / (std_dev * (2.0 * PI).sqrt())) * exponent.exp()
    }

    /// Calculate normal log probability density function
    pub fn normal_log_pdf(x: f64, mean: f64, std_dev: f64) -> f64 {
        if std_dev <= 0.0 {
            return f64::NEG_INFINITY;
        }

        let variance = std_dev * std_dev;
        let exponent = -((x - mean).powi(2)) / (2.0 * variance);
        (-std_dev.ln() - 0.5 * (2.0 * PI).ln()) + exponent
    }

    /// Calculate log factorial using Stirling's approximation for large numbers
    pub fn ln_factorial(n: u32) -> f64 {
        if n <= 20 {
            // Use exact values for small n
            let mut result = 0.0;
            for i in 2..=n {
                result += (i as f64).ln();
            }
            result
        } else {
            // Use Stirling's approximation for large n
            let n_f = n as f64;
            n_f * n_f.ln() - n_f + 0.5 * (2.0 * PI * n_f).ln() + 1.0 / (12.0 * n_f)
                - 1.0 / (360.0 * n_f.powi(3))
        }
    }

    /// Calculate beta function
    pub fn beta_function(a: f64, b: f64) -> f64 {
        if a <= 0.0 || b <= 0.0 {
            return f64::NEG_INFINITY;
        }

        (ln_gamma(a) + ln_gamma(b) - ln_gamma(a + b)).exp()
    }

    /// Calculate log gamma function
    pub fn ln_gamma(x: f64) -> f64 {
        // Lanczos approximation
        const COEFFICIENTS: [f64; 9] = [
            0.999_999_999_999_809_9,
            676.5203681218851,
            -1259.1392167224028,
            771.323_428_777_653_1,
            -176.615_029_162_140_6,
            12.507343278686905,
            -0.13857109526572012,
            9.984_369_578_019_572e-6,
            1.5056327351493116e-7,
        ];

        if x < 0.5 {
            return PI / (PI * x).sin() - ln_gamma(1.0 - x);
        }

        let x_minus_1 = x - 1.0;
        let mut sum = COEFFICIENTS[0];

        for (i, coefficient) in COEFFICIENTS.iter().enumerate().skip(1) {
            sum += coefficient / (x_minus_1 + i as f64);
        }

        let t = x_minus_1 + 7.5;
        let sqrt_2pi = (2.0 * PI).sqrt();

        (sqrt_2pi * t.powf(x_minus_1 + 0.5) * (-t).exp() * sum).ln()
    }
}

/// Matrix operations for genotype likelihoods
pub mod matrix {
    use ndarray::Array2;

    /// Create a 2x2 genotype likelihood matrix
    pub fn genotype_likelihood_matrix_2x2() -> Array2<f64> {
        Array2::zeros((2, 2))
    }

    /// Calculate genotype likelihoods from read data
    pub fn calculate_genotype_likelihoods(
        read_likelihoods: &[f64],
        genotype_combinations: &[(usize, usize)],
    ) -> Vec<f64> {
        genotype_combinations
            .iter()
            .map(|&(i, j)| {
                if i == j {
                    read_likelihoods[i]
                } else {
                    // For heterozygous genotypes, average the likelihoods
                    (read_likelihoods[i] + read_likelihoods[j]) / 2.0
                }
            })
            .collect()
    }

    /// Normalize likelihoods to sum to 1
    pub fn normalize_likelihoods(likelihoods: &mut [f64]) {
        let sum: f64 = likelihoods.iter().sum();
        if sum > 0.0 {
            for likelihood in likelihoods.iter_mut() {
                *likelihood /= sum;
            }
        }
    }

    /// Calculate posterior probabilities from likelihoods and priors
    pub fn calculate_posteriors(likelihoods: &[f64], priors: &[f64]) -> Vec<f64> {
        let mut posteriors = Vec::with_capacity(likelihoods.len());

        for (&likelihood, &prior) in likelihoods.iter().zip(priors.iter()) {
            posteriors.push(likelihood * prior);
        }

        normalize_likelihoods(&mut posteriors);
        posteriors
    }
}

/// Quality score recalibration utilities
pub mod recalibration {
    use super::*;

    /// Recalibrate base quality scores
    pub fn recalibrate_base_quality(
        original_quality: BaseQuality,
        _reported_quality: BaseQuality,
        empirical_quality: BaseQuality,
        num_observations: u32,
    ) -> BaseQuality {
        if num_observations == 0 {
            return original_quality;
        }

        // Simple empirical Bayesian recalibration
        let weight = (num_observations as f64) / (num_observations as f64 + 10.0);
        let recalibrated_q = weight * empirical_quality.value() as f64
            + (1.0 - weight) * original_quality.value() as f64;

        BaseQuality::new(recalibrated_q as u8)
    }

    /// Calculate empirical quality from observed errors
    pub fn empirical_quality_from_errors(num_errors: u32, total_observations: u32) -> BaseQuality {
        if total_observations == 0 {
            return BaseQuality::new(0);
        }

        let error_rate = num_errors as f64 / total_observations as f64;
        BaseQuality::new(phred::error_probability_to_phred(error_rate))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_phred_conversions() {
        let _q = BaseQuality::new(20);
        let error_prob = phred::phred_to_error_probability(20);
        assert!((error_prob - 0.01).abs() < 0.001);

        let q_back = phred::error_probability_to_phred(error_prob);
        assert_eq!(q_back, 20);
    }

    #[test]
    fn test_log_add() {
        let log_a = 0.3679_f64.ln();
        let log_b = 0.1353_f64.ln();

        let log_sum = likelihood::log_add(log_a, log_b);
        let expected_sum = (0.3679_f64 + 0.1353_f64).ln();

        assert!((log_sum - expected_sum).abs() < 0.001);
    }

    #[test]
    fn test_binomial_probability() {
        let prob = stats::binomial_probability(10, 5, 0.5);
        let expected = 252.0 * 0.5_f64.powi(10); // C(10,5) * 0.5^10

        assert!((prob - expected).abs() < 0.001);
    }

    #[test]
    fn test_normal_pdf() {
        let pdf = stats::normal_pdf(0.0, 0.0, 1.0);
        let expected = 1.0 / (2.0 * PI).sqrt();

        assert!((pdf - expected).abs() < 0.001);
    }

    #[test]
    fn test_genotype_likelihoods() {
        let read_likelihoods = vec![0.8, 0.2];
        let genotype_combinations = vec![(0, 0), (0, 1), (1, 1)];

        let genotype_likelihoods =
            matrix::calculate_genotype_likelihoods(&read_likelihoods, &genotype_combinations);

        assert_eq!(genotype_likelihoods, vec![0.8, 0.5, 0.2]);
    }
}
