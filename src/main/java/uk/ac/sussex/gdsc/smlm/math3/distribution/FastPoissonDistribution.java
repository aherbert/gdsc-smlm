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

package uk.ac.sussex.gdsc.smlm.math3.distribution;

import org.apache.commons.math3.exception.NotStrictlyPositiveException;
import org.apache.commons.math3.exception.util.LocalizedFormats;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.MathUtils;

/**
 * Implementation of the Poisson distribution for fast computation of the probability mass function.
 *
 * <p>Adapted from org.apache.commons.math3.distribution.PoissonDistribution.
 *
 * @see <a href="http://en.wikipedia.org/wiki/Poisson_distribution">Poisson distribution
 *      (Wikipedia)</a>
 * @see <a href="http://mathworld.wolfram.com/PoissonDistribution.html">Poisson distribution
 *      (MathWorld)</a>
 */
public final class FastPoissonDistribution {
  /** Mean of the distribution. */
  private double mean;

  /**
   * Creates a new Poisson distribution with specified mean.
   *
   * @param mean the Poisson mean
   * @throws NotStrictlyPositiveException if {@code mean <= 0}.
   */
  public FastPoissonDistribution(double mean) {
    setMean(mean);
  }

  /**
   * Get the mean for the distribution.
   *
   * @return the mean for the distribution.
   */
  public double getMean() {
    return mean;
  }

  /**
   * Sets the mean.
   *
   * @param mean Poisson mean.
   * @throws NotStrictlyPositiveException if {@code mean <= 0}.
   */
  public void setMean(double mean) {
    if (mean <= 0) {
      throw new NotStrictlyPositiveException(LocalizedFormats.MEAN, mean);
    }
    this.mean = mean;
  }

  /**
   * Sets the mean.
   *
   * <p>Does not throw an exception if mean is not strictly positive.
   *
   * @param mean Poisson mean.
   */
  public void setMeanUnsafe(double mean) {
    this.mean = mean;
  }

  /**
   * For a random variable {@code X} whose values are distributed according to this distribution,
   * this method returns {@code P(X = x)}. In other words, this method represents the probability
   * mass function (PMF) for the distribution.
   *
   * @param x the point at which the PMF is evaluated
   * @return the value of the probability mass function at {@code x}
   */
  public double probability(int x) {
    final double logProbability = logProbability(x);
    return logProbability == Double.NEGATIVE_INFINITY ? 0 : FastMath.exp(logProbability);
  }

  /**
   * For a random variable {@code X} whose values are distributed according to this distribution,
   * this method returns {@code log(P(X = x))}, where {@code log} is the natural logarithm. In other
   * words, this method represents the logarithm of the probability mass function (PMF) for the
   * distribution.
   *
   * @param x the point at which the PMF is evaluated
   * @return the logarithm of the value of the probability mass function at {@code x}
   */
  public double logProbability(int x) {
    double ret;
    if (x < 0 || x == Integer.MAX_VALUE) {
      ret = Double.NEGATIVE_INFINITY;
    } else if (x == 0) {
      ret = -mean;
    } else {
      ret = -SaddlePointExpansionCopy.getStirlingError(x)
          - SaddlePointExpansionCopy.getDeviancePart(x, mean) - 0.5 * FastMath.log(MathUtils.TWO_PI)
          - 0.5 * FastMath.log(x);
    }
    return ret;
  }
}
