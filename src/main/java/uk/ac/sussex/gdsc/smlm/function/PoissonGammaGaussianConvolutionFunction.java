/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2022 Alex Herbert
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
import uk.ac.sussex.gdsc.core.utils.MathUtils;

/**
 * Implements the probability density function for a Poisson-Gamma-Gaussian Mixture. The Gaussian is
 * assumed to have mean of zero. If no mean (zero or below) is provided for the Poisson distribution
 * then the probability density function matches that of the Gaussian.
 *
 * <p>The implementation uses full convolution described from Ulbrich &amp; Isacoff (2007). Nature
 * Methods 4, 319-321, SI equation 3:<br> P(D) = A Sum_q e^-u * u^q / q! * 1/sqrt(2pi var) * e ^
 * -((D-q*g)^2 / 2*var)<br> Where:<br> A = normalisation constant var = the variance of the pixel
 * <br> g = the gain of the pixel <br> u = the function value (expected number of photons) <br> D =
 * the observed value at the pixel
 *
 * <p>The likelihood function is designed to model on-chip amplification of a EMCCD/CCD/sCMOS camera
 * which captures a Poisson process of emitted light, converted to electrons on the camera chip,
 * amplified by a gain and then read with Gaussian noise.
 */
public class PoissonGammaGaussianConvolutionFunction
    implements LikelihoodFunction, LogLikelihoodFunction {
  /**
   * The on-chip gain multiplication factor.
   */
  final double gain;

  private final double var;
  private final double sd;
  private final double range;
  private final double twoVar;

  private final double logNormalisationGaussian;

  /**
   * Instantiates a new poisson gaussian convolution function.
   *
   * @param alpha The inverse of the on-chip gain multiplication factor
   * @param variance The variance of the Gaussian distribution at readout (must be positive)
   * @param isVariance Set to true if the input parameter is variance; otherwise is is standard
   *        deviation
   */
  private PoissonGammaGaussianConvolutionFunction(double alpha, double variance,
      boolean isVariance) {
    if (variance <= 0) {
      throw new IllegalArgumentException("Gaussian variance must be strictly positive");
    }
    alpha = Math.abs(alpha);

    this.gain = 1.0 / alpha;
    if (isVariance) {
      sd = Math.sqrt(variance);
      this.var = variance;
    } else {
      sd = variance;
      this.var = sd * sd;
    }
    twoVar = 2 * var;

    // Use a range to cover the Gaussian convolution
    range = 5 * this.sd;

    // Determine the normalisation factor A in the event that the probability
    // distribution is being used as a discrete distribution.
    logNormalisationGaussian = PoissonGaussianFunction.getLogNormalisation(var);
  }

  /**
   * Creates the with standard deviation.
   *
   * @param alpha The inverse of the on-chip gain multiplication factor
   * @param sd The standard deviation of the Gaussian distribution at readout
   * @return the poisson gaussian function 2
   * @throws IllegalArgumentException if the variance is zero or below
   */
  public static PoissonGammaGaussianConvolutionFunction
      createWithStandardDeviation(final double alpha, final double sd) {
    return new PoissonGammaGaussianConvolutionFunction(alpha, sd, false);
  }

  /**
   * Creates the with variance.
   *
   * @param alpha The inverse of the on-chip gain multiplication factor
   * @param var The variance of the Gaussian distribution at readout (must be positive)
   * @return the poisson gaussian function 2
   * @throws IllegalArgumentException if the variance is zero or below
   */
  public static PoissonGammaGaussianConvolutionFunction createWithVariance(final double alpha,
      final double var) {
    return new PoissonGammaGaussianConvolutionFunction(alpha, var, true);
  }

  @Override
  public double likelihood(final double x, final double mu) {
    if (mu <= 0) {
      // If no Poisson mean then just use the Gaussian
      return FastMath.exp((-x * x / twoVar) + logNormalisationGaussian);
    }

    // Note:
    // This is a convolution of two continuous probability distributions.
    // It does not compute a good estimate when the variance is small since the
    // single point approximation to the gaussian is not valid. This is not
    // relevant for a EM-CCD since the variance is likely to be above 10 counts.
    // It also underestimates the cumulative distribution (sum < 1) when the Poisson
    // mean is close to 1 or the gain is small (<4) due to underestimation in the
    // Poisson-Gamma distribution.

    // Use a range to cover the Gaussian convolution
    final double max = x + range;
    if (max < 0) {
      return 0;
    }
    double min = x - range;
    if (min < 0) {
      min = 0;
    }

    return computeP(x, mu, max, min);
  }

  private double computeP(final double x, final double mu, double max, double min) {
    final int cmax = (int) Math.ceil(max);
    final int cmin = (int) Math.floor(min);

    if (cmin == cmax) {
      // Edge case with no range
      return FastMath.exp(
          // Poisson-Gamma
          PoissonGammaFunction.logPoissonGamma(cmin, mu, gain)
              // Gaussian
              - (MathUtils.pow2(cmin - x) / twoVar) + logNormalisationGaussian);
    }

    double pvalue = 0;

    // Overcome the problem with small variance using a set number of steps to
    // cover the range. This effectively makes the Poisson-Gamma a continuous
    // probability distribution.
    // Note:
    // This does not seem to be valid. The Poisson-Gamma is a discrete PMF.
    // The CameraModelAnalysis plugin works with full integration if this function
    // is computed using integer steps.
    //
    // This is computing:
    // Poisson-Gamma PMF(c) x Gaussian PDF(c-o)
    //
    // The solution is to compute:
    // Poisson-Gamma PMF(c) x Gaussian CDF(c-o-0.5,c-o+0.5)
    //
    // This can be done in the PoissonGammaGaussianFunction.

    // if (s < 0)
    // {
    // double step = (max - min) / 10;
    //
    // for (int i = 0; i <= 10; i++)
    // {
    // double c = min + i * step;
    // p += FastMath.exp(
    // // Poisson-Gamma
    // logPoissonGamma(c, e, g)
    // // Gaussian
    // - (Maths.pow2(c - o) / var_by_2) + logNormalisationGaussian);
    // }
    // p *= step;
    // }
    // else
    // {

    for (int c = cmin; c <= cmax; c++) {
      pvalue += FastMath.exp(
          // Poisson-Gamma
          PoissonGammaFunction.logPoissonGamma(c, mu, gain)
              // Gaussian
              - (MathUtils.pow2(c - x) / twoVar) + logNormalisationGaussian);
    }

    // }

    return pvalue;
  }

  @Override
  public double logLikelihood(double x, double mu) {
    if (mu <= 0) {
      // If no Poisson mean then just use the Gaussian
      return (-x * x / twoVar) + logNormalisationGaussian;
    }

    final double max = x + range;
    if (max < 0) {
      return Double.NEGATIVE_INFINITY;
    }
    double min = x - range;
    if (min < 0) {
      min = 0;
    }

    return Math.log(computeP(x, mu, max, min));
  }
}
