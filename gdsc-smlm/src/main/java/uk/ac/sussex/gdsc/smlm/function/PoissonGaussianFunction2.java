/*-
 * #%L
 * Genome Damage and Stability Centre SMLM Package
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2023 Alex Herbert
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

import uk.ac.sussex.gdsc.smlm.utils.StdMath;

/**
 * Implements the probability density function for a Poisson-Gaussian Mixture. The Gaussian is
 * assumed to have mean of zero. If no mean (zero or below) is provided for the Poisson distribution
 * then the probability density function matches that of the Gaussian.
 *
 * <p>The implementation uses the saddle-point approximation described in Snyder, et al (1995)
 * Compensation for readout noise in CCD images. J.Opt. Soc. Am. 12, 272-283. The method is adapted
 * from the C source code provided in the appendix.
 *
 * <p>This is just a wrapper for the PoissonGaussianFunction that handles smart switching between a
 * Gaussian likelihood function and a Poisson-Gaussian likelihood function.
 *
 * <p>The likelihood function is designed to model on-chip amplification of a EMCCD/CCD/sCMOS camera
 * which captures a Poisson process of emitted light, converted to electrons on the camera chip,
 * amplified by a gain and then read with Gaussian noise.
 */
public final class PoissonGaussianFunction2 implements LikelihoodFunction, LogLikelihoodFunction {
  /**
   * The inverse of the on-chip gain multiplication factor.
   */
  final double alpha;

  private boolean usePicardApproximation;
  private final double sigmasquared;

  private final double probabilityNormalisation;
  private final double logNormalisation;
  private final double probabilityNormalisationNoPoisson;
  private final double logNormalisationNoPoisson;

  /**
   * Instantiates a new poisson gaussian function 2.
   *
   * @param alpha The inverse of the on-chip gain multiplication factor
   * @param sigmasquared The variance of the Gaussian distribution at readout (must be positive)
   */
  private PoissonGaussianFunction2(double alpha, double sigmasquared) {
    if (sigmasquared <= 0) {
      throw new IllegalArgumentException("Gaussian variance must be strictly positive");
    }
    alpha = Math.abs(alpha);

    // Apply gain to the readout standard deviation.
    // This compresses the probability distribution by alpha. Thus we can compute the
    // probability using a Poisson or Poisson-Gaussian mixture and then compress the
    // output probability so the cumulative probability is 1 over the uncompressed range.
    sigmasquared *= (alpha * alpha);

    this.alpha = alpha;
    this.sigmasquared = sigmasquared;

    // As per PoissonGaussianFunction
    probabilityNormalisation = PoissonGaussianFunction.NORMALISATION * alpha;
    logNormalisation = PoissonGaussianFunction.LOG_NORMALISATION + Math.log(alpha);
    probabilityNormalisationNoPoisson =
        PoissonGaussianFunction.getProbabilityNormalisation(sigmasquared) * alpha;
    logNormalisationNoPoisson =
        PoissonGaussianFunction.getLogNormalisation(sigmasquared) + Math.log(alpha);
  }

  /**
   * Creates the with standard deviation.
   *
   * @param alpha The inverse of the on-chip gain multiplication factor
   * @param sd The standard deviation of the Gaussian distribution at readout
   * @return the poisson gaussian function 2
   * @throws IllegalArgumentException if the variance is zero or below
   */
  public static PoissonGaussianFunction2 createWithStandardDeviation(final double alpha,
      final double sd) {
    return new PoissonGaussianFunction2(alpha, sd * sd);
  }

  /**
   * Creates the with variance.
   *
   * @param alpha The inverse of the on-chip gain multiplication factor
   * @param var The variance of the Gaussian distribution at readout (must be positive)
   * @return the poisson gaussian function 2
   * @throws IllegalArgumentException if the variance is zero or below
   */
  public static PoissonGaussianFunction2 createWithVariance(final double alpha, final double var) {
    return new PoissonGaussianFunction2(alpha, var);
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

  @Override
  public double likelihood(double x, double mu) {
    // convert to photons
    x *= alpha;
    if (mu <= 0) {
      // If no Poisson mean then just use the Gaussian
      return StdMath.exp(-0.5 * x * x / sigmasquared) * probabilityNormalisationNoPoisson;
    }

    // e *= alpha;
    double saddlepoint =
        (usePicardApproximation) ? PoissonGaussianFunction.picard(x, mu, sigmasquared)
            : PoissonGaussianFunction.pade(x, mu, sigmasquared);
    saddlepoint = PoissonGaussianFunction.newtonIteration(x, mu, sigmasquared, saddlepoint);
    final double logP = PoissonGaussianFunction.spApprox(x, mu, sigmasquared, saddlepoint);
    return StdMath.exp(logP) * probabilityNormalisation;
  }

  @Override
  public double logLikelihood(double x, double mu) {
    // convert to photons
    x *= alpha;
    if (mu <= 0) {
      // If no Poisson mean then just use the Gaussian
      return (-0.5 * x * x / sigmasquared) + logNormalisationNoPoisson;
    }

    // e *= alpha;
    double saddlepoint =
        (usePicardApproximation) ? PoissonGaussianFunction.picard(x, mu, sigmasquared)
            : PoissonGaussianFunction.pade(x, mu, sigmasquared);
    saddlepoint = PoissonGaussianFunction.newtonIteration(x, mu, sigmasquared, saddlepoint);
    final double logP = PoissonGaussianFunction.spApprox(x, mu, sigmasquared, saddlepoint);
    return logP + logNormalisation;
  }
}
