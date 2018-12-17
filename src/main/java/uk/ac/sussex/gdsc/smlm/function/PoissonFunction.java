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

import uk.ac.sussex.gdsc.smlm.math3.distribution.CustomPoissonDistribution;

import org.apache.commons.math3.special.Gamma;

/**
 * Implements the probability density function for a Poisson distribution.
 *
 * <p>This is a simple implementation of the LikelihoodFunction interface.
 *
 * <p>The likelihood function is designed to model on-chip amplification of a EMCCD/CCD/sCMOS camera
 * which captures a Poisson process of emitted light, converted to electrons on the camera chip,
 * amplified by a gain and then read.
 *
 * <p>This function is a scaled Poisson PMF, i.e. using a gain of 2 the integers 1, 3, 5, etc should
 * have no probability from the scaled Poisson PMF.
 */
public class PoissonFunction implements LikelihoodFunction, LogLikelihoodFunction {
  private final CustomPoissonDistribution pd;

  /**
   * The inverse of the on-chip gain multiplication factor.
   */
  final double alpha;

  /**
   * Set to true if the gain was above 1 (inverse gain below 1) and the Poisson distribution should
   * be expanded.
   */
  final boolean expand;

  /**
   * Instantiates a new poisson function.
   *
   * @param alpha The inverse of the on-chip gain multiplication factor
   */
  public PoissonFunction(double alpha) {
    this.alpha = Math.abs(alpha);
    expand = (alpha < 1);
    pd = new CustomPoissonDistribution(null, 1);
  }

  /**
   * {@inheritDoc}
   *
   * <p>This is a PMF.
   */
  @Override
  public double likelihood(double observed, double mu) {
    if (mu <= 0) {
      return 0;
    }

    final int x = getX(observed);
    if (x < 0) {
      return 0;
    }

    // Compute the integer intervals of the Poisson to sample
    // The entire P(x) of the scaled Poisson is assigned to the nearest integer
    // Find the limits of the integer range.
    final double min = (x - 0.5) * alpha;
    final double max = (x + 0.5) * alpha;

    // The first integer that would be rounded to x
    int imin = (int) Math.ceil(min);
    if (imin < 0) {
      imin = 0;
    }

    if (expand) {
      // The PMF was expanded so either 1 or 0 values fall in this range
      // When rounding an equality at the upper edge should be assigned
      // to the next interval.
      if (imin >= max) {
        return 0;
      }

      pd.setMeanUnsafe(mu);
      return pd.probability(imin);
    }
    // The PMF was contracted so 1 or more values fall in this range

    int imax = (int) Math.floor(max);
    // When rounding an equality at the upper edge should be assigned
    // to the next interval.
    if (imax == max) {
      imax--;
    }

    pd.setMeanUnsafe(mu);
    if (imin == imax) {
      return pd.probability(imin);
    }

    double pvalue = 0;
    for (int i = imin; i <= imax; i++) {
      pvalue += pd.probability(i);
    }
    return pvalue;
  }

  private static int getX(double observed) {
    // XXX: Change to throw an exception if o is not an integer?
    return (int) Math.round(observed);
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

  /** {@inheritDoc} */
  @Override
  public double logLikelihood(double observed, double mu) {
    // As above but with log output
    if (mu <= 0) {
      return Double.NEGATIVE_INFINITY;
    }

    final int x = getX(observed);
    if (x < 0) {
      return Double.NEGATIVE_INFINITY;
    }

    final double min = (x - 0.5) * alpha;
    final double max = (x + 0.5) * alpha;

    int imin = (int) Math.ceil(min);
    if (imin < 0) {
      imin = 0;
    }

    if (expand) {
      if (imin >= max) {
        return Double.NEGATIVE_INFINITY;
      }

      pd.setMeanUnsafe(mu);
      return pd.logProbability(imin);
    }

    int imax = (int) Math.floor(max);
    if (imax == max) {
      imax--;
    }

    pd.setMeanUnsafe(mu);
    if (imin == imax) {
      return pd.logProbability(imin);
    }

    double pvalue = 0;
    for (int i = imin; i <= imax; i++) {
      pvalue += pd.probability(i);
    }
    return Math.log(pvalue);
  }
}
