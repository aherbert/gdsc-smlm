/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2019 Alex Herbert
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

import uk.ac.sussex.gdsc.smlm.math3.distribution.FastPoissonDistribution;

import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.FastMath;

/**
 * Implements the probability density function for a Poisson distribution.
 *
 * <p>This is a simple implementation of the LikelihoodFunction interface.
 *
 * <p>The likelihood function is designed to model on-chip amplification of a EMCCD/CCD/sCMOS camera
 * which captures a Poisson process of emitted light, converted to electrons on the camera chip,
 * amplified by a gain and then read.
 *
 * <p>This function allows any input value X to generate a probability, even if the value is between
 * steps of the scaled PMF of the poisson distribution, i.e. using a gain of 2 the integers 1, 3, 5,
 * etc should have no probability from the scaled Poisson PMF but this function portions the
 * probability accordingly.
 */
public class InterpolatedPoissonFunction
    implements GradientLikelihoodFunction, LogLikelihoodFunction {
  private final FastPoissonDistribution pd;

  /**
   * The inverse of the on-chip gain multiplication factor.
   */
  final double alpha;

  /**
   * The log of the inverse on-chip gain multiplication factor.
   */
  final double logAlpha;

  /**
   * Allow non-integer observed values.
   */
  final boolean nonInteger;

  /**
   * Instantiates a new interpolated poisson function.
   *
   * @param alpha The inverse of the on-chip gain multiplication factor
   * @param nonInteger Allow non-integer observed values
   */
  public InterpolatedPoissonFunction(double alpha, boolean nonInteger) {
    this.alpha = Math.abs(alpha);
    logAlpha = Math.log(alpha);
    this.nonInteger = nonInteger;
    pd = (nonInteger) ? null : new FastPoissonDistribution(1);
  }

  @Override
  public double likelihood(double x, double mu) {
    if (x < 0 || mu <= 0) {
      return 0;
    }

    // convert to photons
    x *= alpha;

    // Allow non-integer observed value using the gamma function to provide a factorial for
    // non-integer values
    // PMF(l,k) = e^-l * l^k / gamma(k+1)
    // log(PMF) = -l + k * log(l) - logGamma(k+1)
    if (nonInteger) {
      // return (FastMath.exp(-e) * Math.pow(e, o) / factorial(o)) * alpha;

      final double ll = -mu + x * Math.log(mu) - logFactorial(x);
      return FastMath.exp(ll) * alpha;
    }

    pd.setMeanUnsafe(mu);
    return pd.probability((int) x) * alpha;
  }

  /**
   * {@inheritDoc}
   *
   * <p>When evaluating using non-integer values the gradient for {@code 0 < o/gain < 1} is
   * incorrect. In this region there is no definition for the factorial required for the derivative
   * {@code (o/gain - 1)!}.
   */
  @Override
  public double likelihood(double x, double mu, double[] dpDmu) {
    if (x < 0 || mu <= 0) {
      dpDmu[0] = 0;
      return 0;
    }

    // convert to photons
    x *= alpha;

    // PMF(l,k) = e^-l * l^k / k!
    // PMF'(l,k) = e^-l * k*l^(k-1) / k! + -e^-l * l^k / k!
    // = e^-l * l^(k-1) / (k-1)! - e^-l * l^k / k!
    // = PMF(l,k-1) - PMF(l,k)

    double lk;
    double lkm1;

    if (nonInteger) {
      final double loge = Math.log(mu);
      double ll = -mu + x * loge - logFactorial(x);
      lk = FastMath.exp(ll);
      if (x == mu) {
        // Special case
        dpDmu[0] = 0;
        return lk * alpha;
      }

      if (x >= 1) {
        // In contrast to the logFactorial(o-1)
        // this continues to use the logGamma function even when o-1 < 1.
        // It creates the correct gradient down to o==1.
        ll = -mu + (x - 1) * loge - Gamma.logGamma(x);
        lkm1 = FastMath.exp(ll);
      } else if (x > 0) {
        // x is between 0 and 1.
        // There is no definition for the factorial (x-1)!

        // ll = -e - Gamma.logGamma(x);
        // ll = -e - Math.abs(x - 1) * loge - Gamma.logGamma(x);

        // This continues to be the best match to the numerical gradient
        // even though it is impossible. It works because the gamma function is still
        // defined when x>0.
        ll = -mu + (x - 1) * loge - Gamma.logGamma(x);
        lkm1 = FastMath.exp(ll);
      } else {
        lkm1 = 0;
      }
    } else {
      final int k = (int) x;
      pd.setMeanUnsafe(mu);
      lk = pd.probability(k);
      if (k == mu) {
        // Special case
        dpDmu[0] = 0;
        return lk * alpha;
      }

      if (k != 0) {
        lkm1 = pd.probability(k - 1);
      } else {
        lkm1 = 0;
      }
    }

    dpDmu[0] = (lkm1 - lk) * alpha;
    return lk * alpha;
  }

  /**
   * Return the log of the factorial for the given real number, using the gamma function.
   *
   * @param value the number
   * @return the log factorial
   */
  public static double logFactorial(double value) {
    if (value <= 1) {
      return 0;
    }
    return Gamma.logGamma(value + 1);
  }

  /**
   * Return the factorial for the given real number, using the gamma function.
   *
   * @param value the number
   * @return the factorial
   */
  public static double factorial(double value) {
    if (value <= 1) {
      return 1;
    }
    return Gamma.gamma(value + 1);
  }

  @Override
  public double logLikelihood(double x, double mu) {
    if (x < 0 || mu <= 0) {
      return Double.NEGATIVE_INFINITY;
    }

    // convert to photons
    x *= alpha;

    // Allow non-integer observed value using the gamma function to provide a factorial for
    // non-integer values
    // PMF(l,k) = e^-l * l^k / gamma(k+1)
    // log(PMF) = -l + k * log(l) - logGamma(k+1)
    if (nonInteger) {
      // return (FastMath.exp(-e) * Math.pow(e, o) / factorial(o)) * alpha;

      final double ll = -mu + x * Math.log(mu) - logFactorial(x);
      return ll + logAlpha;
    }

    pd.setMeanUnsafe(mu);
    return pd.logProbability((int) x) + logAlpha;
  }
}
