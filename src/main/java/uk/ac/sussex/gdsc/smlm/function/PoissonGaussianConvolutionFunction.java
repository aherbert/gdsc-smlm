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

import uk.ac.sussex.gdsc.core.utils.MathUtils;

/**
 * Implements the probability density function for a Poisson-Gaussian Mixture. The Gaussian is
 * assumed to have mean of zero. If no mean (zero or below) is provided for the Poisson distribution
 * then the probability density function matches that of the Gaussian. <p> The implementation uses
 * full convolution described from Huang, et al (2013), Supplementary Notes Eq 1.1:<br> P(D) = A
 * Sum_q e^-u * u^q / q! * 1/sqrt(2pi var) * e ^ -((D-q*g)^2 / 2*var)<br> Where:<br> A =
 * normalisation constant var = the variance of the pixel <br> g = the gain of the pixel <br> u =
 * the function value (expected number of photons) <br> D = the observed value at the pixel <p> The
 * likelihood function is designed to model on-chip amplification of a EMCCD/CCD/sCMOS camera which
 * captures a Poisson process of emitted light, converted to electrons on the camera chip, amplified
 * by a gain and then read with Gaussian noise.
 */
public class PoissonGaussianConvolutionFunction
    implements LikelihoodFunction, LogLikelihoodFunction {
  private static final LogFactorial logFactorial = new LogFactorial();

  /**
   * The on-chip gain multiplication factor.
   */
  final double g;

  private final double var;
  private final double s;
  private final double var_by_2;
  private final double sqrt_var_by_2;

  private final double logNormalisationGaussian;

  private boolean computePMF = false;

  /**
   * Instantiates a new poisson gaussian convolution function.
   *
   * @param alpha The inverse of the on-chip gain multiplication factor
   * @param variance The variance of the Gaussian distribution at readout (must be positive)
   * @param isVariance Set to true if the input parameter is variance; otherwise is is standard
   *        deviation
   */
  private PoissonGaussianConvolutionFunction(double alpha, double variance, boolean isVariance) {
    if (variance <= 0) {
      throw new IllegalArgumentException("Gaussian variance must be strictly positive");
    }
    alpha = Math.abs(alpha);

    this.g = 1.0 / alpha;
    if (isVariance) {
      s = Math.sqrt(variance);
      this.var = variance;
    } else {
      s = variance;
      this.var = s * s;
    }
    var_by_2 = var * 2;
    sqrt_var_by_2 = Math.sqrt(var_by_2);

    // Determine the normalisation factor A in the event that the probability
    // distribution is being used as a discrete distribution.
    logNormalisationGaussian = PoissonGaussianFunction.getLogNormalisation(var);
  }

  /**
   * Creates the with standard deviation.
   *
   * @param alpha The inverse of the on-chip gain multiplication factor
   * @param s The standard deviation of the Gaussian distribution at readout
   * @return the poisson gaussian function 2
   * @throws IllegalArgumentException if the variance is zero or below
   */
  public static PoissonGaussianConvolutionFunction createWithStandardDeviation(final double alpha,
      final double s) {
    return new PoissonGaussianConvolutionFunction(alpha, s, false);
  }

  /**
   * Creates the with variance.
   *
   * @param alpha The inverse of the on-chip gain multiplication factor
   * @param var The variance of the Gaussian distribution at readout (must be positive)
   * @return the poisson gaussian function 2
   * @throws IllegalArgumentException if the variance is zero or below
   */
  public static PoissonGaussianConvolutionFunction createWithVariance(final double alpha,
      final double var) {
    return new PoissonGaussianConvolutionFunction(alpha, var, true);
  }

  /**
   * {@inheritDoc} <p> The output is a PDF or PMF depending on the value of {@link #isComputePMF()}.
   * If set to true the function does not error if the input x is non-integer.
   *
   * @see uk.ac.sussex.gdsc.smlm.function.LikelihoodFunction#likelihood(double, double)
   * @see #isComputePMF()
   */
  @Override
  public double likelihood(double o, double e) {
    if (e <= 0) {
      // If no Poisson mean then just use the Gaussian
      if (computePMF) {
        final double x = Math.round(o);
        return (gaussianCDF(x + 0.5) - gaussianCDF(x - 0.5)) * 0.5;
      }
      return FastMath.exp((-0.5 * o * o / var) + logNormalisationGaussian);
    }
    // Use same nomenclature as Huang et al

    final double u = e; // expected photoelectrons
    final double D = o; // Camera counts
    // g == gain
    // var = readout variance

    // This is the probability of a Poisson convolved with a Gaussian.
    // Evaluate the Poisson only in the range where the Gaussian is significant.
    // I.e. when the Gaussian probability is zero then the Poisson is not relevant.
    // Use +/- 5 SD
    // x = D - q*g => q = (D-x) / g
    int qmax = (int) Math.ceil((D + 5 * s) / g);
    if (qmax < 0) {
      return 0;
    }
    int qmin = (int) Math.floor((D - 5 * s) / g);
    if (qmin < 0) {
      qmin = 0;
      // Collision check to avoid double computing
      if (qmax == 0) {
        qmax++;
      }
    }

    // Note: If D is camera counts then it will likely be limited to a 16-bit range
    // Assuming the gain is at least 1 then the max q is:
    // 65536 + 5 * s => This is an acceptable table size to pre-compute the log
    // factorial if s is reasonable.

    logFactorial.ensureRange(qmin, qmax);

    final double logu = Math.log(u);
    double p = 0;

    // Optionally use the error function for a full convolution between
    // the Poisson PMF and Gaussian PDF
    if (computePMF) {
      for (int q = qmin; q <= qmax; q++) {
        final double x = getX(D, q);
        p +=
            // Poisson PMF
            FastMath.exp(q * logu - u - logFactorial.getLogFUnsafe(q)) *
            // Gaussian CDF
                (gaussianCDF(x + 0.5) - gaussianCDF(x - 0.5)) * 0.5;
      }
    } else {
      for (int q = qmin; q <= qmax; q++) {
        p += FastMath.exp(
            // Poisson
            q * logu - u - logFactorial.getLogFUnsafe(q)
            // Gaussian
                - (MathUtils.pow2(D - q * g) / var_by_2) + logNormalisationGaussian);
      }
    }

    // Determine normalisation
    // Note: This is needed when using this as a discrete probability distribution,
    // e.g. input observed count is integer

    return p;
  }

  private double getX(final double D, int q) {
    // Do not round to compute the convolution point x
    // return Math.round(D - q * g);
    return D - q * g;
  }

  /**
   * Gaussian CDF.
   *
   * @param x the x
   * @return the cumulative density
   */
  double gaussianCDF(final double x) {
    // return org.apache.commons.math3.special.CDF.erf(x / sqrt_var_by_2);
    // This may not be precise enough.
    // Absolute error is <3e-7. Not sure what relative error is.
    // The standard CDF is much slower.
    return Erf.erf(x / sqrt_var_by_2);
  }

  /**
   * {@inheritDoc} <p> The output is a PDF or PMF depending on the value of {@link #isComputePMF()}.
   * If set to true the function does not error if the input x is non-integer.
   *
   * @see uk.ac.sussex.gdsc.smlm.function.LogLikelihoodFunction#logLikelihood(double, double)
   * @see #isComputePMF()
   */
  @Override
  public double logLikelihood(double o, double e) {
    // As above but return the log

    if (e <= 0) {
      // If no Poisson mean then just use the Gaussian
      if (computePMF) {
        final double x = Math.round(o);
        return Math.log((gaussianCDF(x + 0.5) - gaussianCDF(x - 0.5)) * 0.5);
      }
      return (-0.5 * o * o / var) + logNormalisationGaussian;
    }
    final double u = e; // expected photoelectrons
    final double D = o; // Camera counts
    int qmax = (int) Math.ceil((D + 5 * s) / g);
    if (qmax < 0) {
      return Double.NEGATIVE_INFINITY;
    }
    int qmin = (int) Math.floor((D - 5 * s) / g);
    if (qmin < 0) {
      qmin = 0;
      // Collision check to avoid double computing
      if (qmax == 0) {
        qmax++;
      }
    }
    logFactorial.ensureRange(qmin, qmax);
    final double logu = Math.log(u);
    double p = 0;
    if (computePMF) {
      for (int q = qmin; q <= qmax; q++) {
        final double x = getX(D, q);
        p +=
            // Poisson PMF
            FastMath.exp(q * logu - u - logFactorial.getLogFUnsafe(q)) *
            // Gaussian CDF
                (gaussianCDF(x + 0.5) - gaussianCDF(x - 0.5)) * 0.5;
      }
    } else {
      for (int q = qmin; q <= qmax; q++) {
        p += FastMath.exp(
            // Poisson
            q * logu - u - logFactorial.getLogFUnsafe(q)
            // Gaussian
                - (MathUtils.pow2(D - q * g) / var_by_2) + logNormalisationGaussian);
      }
    }
    return Math.log(p);
  }

  /**
   * Checks if computing a PMF(X=x) (for integer x) using the Gaussian CDF to convolve with the
   * Poisson PMF. <p> The default is a PDF(X=x). If set to true the function
   * {@link #likelihood(double, double)} does not error if the input x is non-integer.
   *
   * @return true, if computing a PMF(X=x).
   */
  public boolean isComputePMF() {
    return computePMF;
  }

  /**
   * Set to true if computing a PMF(X=x) (for integer x) using the Gaussian CDF to convolve with the
   * Poisson PMF. <p> The default is a PDF(X=x). If set to true the function
   * {@link #likelihood(double, double)} does not error if the input x is non-integer.
   *
   * @param computePMF the new use CDF flag
   */
  public void setComputePMF(boolean computePMF) {
    this.computePMF = computePMF;
  }
}
