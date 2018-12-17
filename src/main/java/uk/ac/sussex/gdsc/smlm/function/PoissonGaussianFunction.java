/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */

package uk.ac.sussex.gdsc.smlm.function;

import org.apache.commons.math3.util.FastMath;

/**
 * Implements the probability density function for a Poisson-Gaussian Mixture. The Gaussian is
 * assumed to have mean of zero. If no mean (zero or below) is provided for the Poisson distribution
 * then the probability density function matches that of the Gaussian.
 *
 * <p>The implementation uses the saddle-point approximation described in Snyder, et al (1995)
 * Compensation for readout noise in CCD images. J.Opt. Soc. Am. 12, 272-283. The method is adapted
 * from the C source code provided in the appendix.
 *
 * <p>The likelihood function is designed to model on-chip amplification of a EMCCD/CCD/sCMOS camera
 * which captures a Poisson process of emitted light, converted to electrons on the camera chip,
 * amplified by a gain and then read with Gaussian noise.
 */
public class PoissonGaussianFunction implements LikelihoodFunction, LogLikelihoodFunction {
  /**
   * The inverse of the on-chip gain multiplication factor.
   */
  final double alpha;

  /**
   * The log of {@link #alpha}.
   */
  final double logAlpha;

  private static final double EPSILON = 1e-4; // 1e-6

  /** The Gaussian normalisation: <code>1 / Math.sqrt(2 * Math.PI)</code>. */
  static final double NORMALISATION = 1 / Math.sqrt(2 * Math.PI);
  /** The log Gaussian normalisation: <code>Math.log(1 / Math.sqrt(2 * Math.PI))</code>. */
  static final double LOG_NORMALISATION = Math.log(NORMALISATION);

  /**
   * Number of Picard iterations to use.
   */
  private static final int NUM_PICARD = 3;

  private boolean usePicardApproximation;
  private final double mu;
  private final double sigmaSquared;
  private final boolean noPoisson;

  private final double probabilityNormalisation;
  private final double logNormalisation;

  /**
   * Instantiates a new poisson gaussian function.
   *
   * @param alpha The inverse of the on-chip gain multiplication factor
   * @param mu The mean of the Poisson distribution before gain
   * @param sigmaSquared The variance of the Gaussian distribution at readout (must be positive)
   */
  private PoissonGaussianFunction(double alpha, double mu, double sigmaSquared) {
    if (sigmaSquared <= 0) {
      throw new IllegalArgumentException("Gaussian variance must be strictly positive");
    }
    alpha = Math.abs(alpha);
    noPoisson = (mu <= 0);

    // Apply gain to the mean and readout standard deviation.
    // This compresses the probability distribution by alpha. Thus we can compute the
    // probability using a Poisson or Poisson-Gaussian mixture and then compress the
    // output probability so the cumulative probability is 1 over the uncompressed range.
    sigmaSquared *= (alpha * alpha);

    // System.out.printf("mu=%f s^2=%f\n", mu, sigmaSquared);

    this.alpha = alpha;
    this.mu = mu;
    this.sigmaSquared = sigmaSquared;

    probabilityNormalisation =
        ((noPoisson) ? getProbabilityNormalisation(sigmaSquared) : 1) * alpha;

    logAlpha = Math.log(alpha);
    logNormalisation = ((noPoisson) ? getLogNormalisation(sigmaSquared)
        // Note that this uses the LOG_NORMALISATION (not zero) as the logProbability function
        // uses the getPseudoLikelihood function. It would be zero if the static logProbability
        // was called.
        : LOG_NORMALISATION) + logAlpha;
  }

  /**
   * Creates the with standard deviation.
   *
   * @param alpha The inverse of the on-chip gain multiplication factor
   * @param mu The mean of the Poisson distribution before gain
   * @param sd The standard deviation of the Gaussian distribution at readout
   * @return the poisson gaussian function
   * @throws IllegalArgumentException if the mean or variance is zero or below
   */
  public static PoissonGaussianFunction createWithStandardDeviation(final double alpha,
      final double mu, final double sd) {
    return new PoissonGaussianFunction(alpha, mu, sd * sd);
  }

  /**
   * Creates the with variance.
   *
   * @param alpha The inverse of the on-chip gain multiplication factor
   * @param mu The mean of the Poisson distribution before gain
   * @param var The variance of the Gaussian distribution at readout (must be positive)
   * @return the poisson gaussian function
   * @throws IllegalArgumentException if the mean or variance is zero or below
   */
  public static PoissonGaussianFunction createWithVariance(final double alpha, final double mu,
      final double var) {
    return new PoissonGaussianFunction(alpha, mu, var);
  }

  /**
   * Get the probability of observation x.
   *
   * @param x The observation value (after gain)
   * @return The probability
   */
  public double probability(double x) {
    return probability(x, mu);
  }

  /**
   * Get the probability of observation x.
   *
   * <p>Note the normalisation will be based on the mu used in the constructor of the function. So
   * it may be wrong if the mu was below zero on construction and is above zero now, or vice-versa.
   *
   * @param x The observation value (after gain)
   * @param mu The mean of the Poisson distribution (before gain)
   * @return The probability
   */
  public double probability(double x, double mu) {
    // convert to photons
    x *= alpha;

    return (noPoisson) ? FastMath.exp(-0.5 * x * x / sigmaSquared) * probabilityNormalisation
        : getProbability(x, mu, sigmaSquared, usePicardApproximation) * probabilityNormalisation;
  }

  /**
   * Get the log(p) of observation x.
   *
   * @param x The observation value
   * @return The log of the probability
   */
  public double logProbability(double x) {
    return logProbability(x, mu);
  }

  /**
   * Get the log(p) of observation x.
   *
   * <p>Note the normalisation will be based on the mu used in the constructor of the function. So
   * it may be wrong if the mu was below zero on construction and is above zero now, or vice-versa.
   *
   * @param x The observation value
   * @param mu The mean of the Poisson distribution
   * @return The log of the probability
   */
  public double logProbability(double x, double mu) {
    // convert to photons
    x *= alpha;

    return (noPoisson) ? (-0.5 * x * x / sigmaSquared) + logNormalisation
        : getPseudoLikelihood(x, mu, sigmaSquared, usePicardApproximation) + logNormalisation;
  }

  /**
   * Get the probability of observation x.
   *
   * @param x The observation value
   * @param mu The mean of the Poisson distribution (if zero or below the probability is that of the
   *        Gaussian)
   * @param sigmaSquared The variance of the Gaussian distribution (must be positive)
   * @param usePicardApproximation Use the Picard approximation for the initial saddle point. The
   *        default is Pade.
   * @return The probability
   */
  public static double probability(final double x, final double mu, final double sigmaSquared,
      final boolean usePicardApproximation) {
    // If no Poisson mean then just use the Gaussian
    if (mu <= 0) {
      return FastMath.exp(-0.5 * x * x / sigmaSquared) * getProbabilityNormalisation(sigmaSquared);
    }
    return getProbability(x, mu, sigmaSquared, usePicardApproximation);
  }

  /**
   * Gets the probability normalisation for the Gaussian distribution so the cumulative probability
   * is 1.
   *
   * @param sigmaSquared the sigma squared
   * @return the probability normalisation
   */
  static double getProbabilityNormalisation(double sigmaSquared) {
    return NORMALISATION / Math.sqrt(sigmaSquared);
  }

  /**
   * Get the probability of observation x.
   *
   * @param x The observation value
   * @param mu The mean of the Poisson distribution (must be positive)
   * @param sigmaSquared The variance of the Gaussian distribution (must be positive)
   * @param usePicardApproximation Use the Picard approximation for the initial saddle point. The
   *        default is Pade.
   * @return The probability
   */
  private static double getProbability(final double x, final double mu, final double sigmaSquared,
      final boolean usePicardApproximation) {
    double saddlepoint =
        (usePicardApproximation) ? picard(x, mu, sigmaSquared) : pade(x, mu, sigmaSquared);
    saddlepoint = newton_iteration(x, mu, sigmaSquared, saddlepoint);
    final double logP = sp_approx(x, mu, sigmaSquared, saddlepoint);
    return FastMath.exp(logP) * NORMALISATION;
  }

  /**
   * Get the log(p) of observation x.
   *
   * @param x The observation value
   * @param mu The mean of the Poisson distribution (if zero or below the probability is that of the
   *        Gaussian)
   * @param sigmaSquared The variance of the Gaussian distribution (must be positive)
   * @param usePicardApproximation Use the Picard approximation for the initial saddle point. The
   *        default is Pade.
   * @return The log of the probability
   */
  public static double logProbability(final double x, final double mu, final double sigmaSquared,
      final boolean usePicardApproximation) {
    // If no Poisson mean then just use the Gaussian
    if (mu <= 0) {
      return (-0.5 * x * x / sigmaSquared) + getLogNormalisation(sigmaSquared);
    }

    return getPseudoLikelihood(x, mu, sigmaSquared, usePicardApproximation) + LOG_NORMALISATION;
  }

  /**
   * Gets the log probability normalisation for the Gaussian distribution.
   *
   * @param sigmaSquared the sigma squared
   * @return the log probability normalisation
   */
  static double getLogNormalisation(double sigmaSquared) {
    return LOG_NORMALISATION - Math.log(sigmaSquared) * 0.5;
  }

  /**
   * Get a pseudo-likelihood of observation x.
   *
   * <p>This is equivalent to the {@link #logProbability(double)} without the normalisation of the
   * probability density function to 1. It differs by a constant value of -log(1 / sqrt(2 * PI)).
   * This function is suitable for use as the likelihood function in maximum likelihood estimation
   * since all values will differ by the same constant but will evaluate faster.
   *
   * @param x The observation value
   * @param mu The mean of the Poisson distribution (if zero or below the probability is that of the
   *        Gaussian)
   * @param sigmaSquared The variance of the Gaussian distribution (must be positive)
   * @param usePicardApproximation Use the Picard approximation for the initial saddle point. The
   *        default is Pade.
   * @return The probability
   */
  public static double pseudoLikelihood(final double x, final double mu, final double sigmaSquared,
      final boolean usePicardApproximation) {
    // If no Poisson mean then just use the Gaussian
    if (mu <= 0) {
      return (-0.5 * x * x / sigmaSquared);
    }

    return getPseudoLikelihood(x, mu, sigmaSquared, usePicardApproximation);
  }

  /**
   * Get a pseudo-likelihood of observation x.
   *
   * <p>This is equivalent to the {@link #logProbability(double)} without the normalisation of the
   * probability density function to 1. It differs by a constant value of -log(1 / sqrt(2 * PI)).
   * This function is suitable for use as the likelihood function in maximum likelihood estimation
   * since all values will differ by the same constant but will evaluate faster.
   *
   * @param x The observation value
   * @param mu The mean of the Poisson distribution (must be positive)
   * @param sigmaSquared The variance of the Gaussian distribution (must be positive)
   * @param usePicardApproximation Use the Picard approximation for the initial saddle point. The
   *        default is Pade.
   * @return The probability
   */
  private static double getPseudoLikelihood(final double x, final double mu,
      final double sigmaSquared, final boolean usePicardApproximation) {
    double saddlepoint =
        (usePicardApproximation) ? picard(x, mu, sigmaSquared) : pade(x, mu, sigmaSquared);
    saddlepoint = newton_iteration(x, mu, sigmaSquared, saddlepoint);
    return sp_approx(x, mu, sigmaSquared, saddlepoint);
  }

  /**
   * Return the initial saddle point estimated by the Pade approximation.
   *
   * @param x the x
   * @param mu the mu
   * @param sigmaSquared the sigma squared
   * @return The saddle point
   */
  static double pade(final double x, final double mu, final double sigmaSquared) {
    final double bterm = x - 2 * sigmaSquared - mu;

    // Original code
    // return -Math.log(0.5 * (bterm + Math.sqrt(bterm * bterm + 4 * mu * (2 * sigmaSquared + x))) /
    // mu);

    // Check for negative sqrt
    final double argument_to_sqrt = bterm * bterm + 4 * mu * (2 * sigmaSquared + x);
    // This check is needed if x is very negative
    if (argument_to_sqrt < 0) {
      // Revert to Taylor approximation
      return (mu - x) / (mu + sigmaSquared);
    }

    // Check for negative log
    final double argument_to_log = 0.5 * (bterm + Math.sqrt(argument_to_sqrt)) / mu;
    if (argument_to_log <= 0) {
      // Revert to Taylor approximation
      return (mu - x) / (mu + sigmaSquared);
    }
    return -Math.log(argument_to_log);
  }

  /**
   * Return the initial saddle point estimated by the Picard approximation.
   *
   * @param x the x
   * @param mu the mu
   * @param sigmaSquared the sigma squared
   * @return The saddle point
   */
  static double picard(final double x, final double mu, final double sigmaSquared) {
    // Use Taylor approximation to obtain the starting point for Picard iteration
    final double taylor = (mu - x) / (mu + sigmaSquared);
    double saddlepoint = taylor;
    for (int i = 0; i < NUM_PICARD; i++) {
      final double argument_to_log = mu / (x + sigmaSquared * saddlepoint);
      if (argument_to_log <= 0) {
        // Break out of loop if argument to log goes negative
        return taylor;
      }
      saddlepoint = Math.log(argument_to_log);
    }
    return saddlepoint;
  }

  /**
   * Returns the saddlepoint found by Newton iteration for a given x, mu, sigmaSquared and an
   * initial estimate of the saddle point (found with either the Pade or Picard approach).
   *
   * @param x the x
   * @param mu the mu
   * @param sigmaSquared the sigma squared
   * @param initial_saddlepoint the initial saddlepoint
   * @return The saddle point
   */
  static double newton_iteration(final double x, final double mu, final double sigmaSquared,
      final double initial_saddlepoint) {
    double change = 0;
    double saddlepoint = initial_saddlepoint;

    // Original code can infinite loop
    // do
    // {
    // final double mu_exp_minus_s = mu * FastMath.exp(-saddlepoint);
    // change = (x + sigmaSquared * saddlepoint - mu_exp_minus_s) / (sigmaSquared + mu_exp_minus_s);
    // saddlepoint -= change;
    // } while (FastMath.abs(change) > EPSILON * FastMath.abs(saddlepoint));
    // return saddlepoint;

    // Limit the number of iterations
    for (int i = 0; i < 20; i++) {
      final double mu_exp_minus_s = mu * FastMath.exp(-saddlepoint);
      change = (x + sigmaSquared * saddlepoint - mu_exp_minus_s) / (sigmaSquared + mu_exp_minus_s);
      saddlepoint -= change;
      if (FastMath.abs(change) <= EPSILON * FastMath.abs(saddlepoint)) {
        return saddlepoint;
      }
    }

    // This happens when we cannot converge
    // System.out.printf("No Newton convergence: x=%f, mu=%f, s2=%f, %s -> %s : logP=%f, p=%f\n", x,
    // mu, sigmaSquared,
    // initial_saddlepoint, saddlepoint, sp_approx(x, mu, sigmaSquared, initial_saddlepoint),
    // FastMath.exp(sp_approx(x, mu, sigmaSquared, initial_saddlepoint)) * NORMALISATION);

    return saddlepoint; // initial_saddlepoint;
  }

  /**
   * Return the saddlepoint approximation to the log of p(x,mu,sigmaSquared) given the saddle point
   * found by the Newton iteration. Remember the sqrt(2*PI) factor has been left out.
   *
   * @param x the x
   * @param mu the mu
   * @param sigmaSquared the sigma squared
   * @param saddlepoint the saddlepoint
   * @return The saddlepoint approximation
   */
  static double sp_approx(final double x, final double mu, final double sigmaSquared,
      final double saddlepoint) {
    final double mu_exp_minus_s = mu * FastMath.exp(-saddlepoint);
    final double phi2 = sigmaSquared + mu_exp_minus_s;
    return -mu + (saddlepoint * (x + 0.5 * sigmaSquared * saddlepoint)) + mu_exp_minus_s
        - 0.5 * Math.log(phi2);
  }

  /**
   * Return if using the Picard approximation for the initial saddle point.
   *
   * @return True if using the Picard approximation
   */
  public boolean isUsePicardApproximation() {
    return usePicardApproximation;
  }

  /**
   * Specify whether to use the Picard approximation for the initial saddle point. The alternative
   * is Pade.
   *
   * @param usePicardApproximation True to use the Picard approximation
   */
  public void setUsePicardApproximation(boolean usePicardApproximation) {
    this.usePicardApproximation = usePicardApproximation;
  }

  /** {@inheritDoc} */
  @Override
  public double likelihood(double x, double mu) {
    return probability(x, mu);
  }

  /** {@inheritDoc} */
  @Override
  public double logLikelihood(double x, double mu) {
    return logProbability(x, mu);
  }
}
