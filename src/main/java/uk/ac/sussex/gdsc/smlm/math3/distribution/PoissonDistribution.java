/*
 * Licensed to the Apache Software Foundation (ASF) under one or more contributor license
 * agreements. See the NOTICE file distributed with this work for additional information regarding
 * copyright ownership. The ASF licenses this file to You under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance with the License. You may obtain a
 * copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software distributed under the License
 * is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express
 * or implied. See the License for the specific language governing permissions and limitations under
 * the License.
 */

package uk.ac.sussex.gdsc.smlm.math3.distribution;

import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.MathUtils;
import uk.ac.sussex.gdsc.core.utils.ValidationUtils;

/**
 * Implementation of the Poisson distribution for computation of the probability mass function.
 *
 * <p>Adapted from org.apache.commons.math3.distribution.PoissonDistribution. The code has been
 * updated to: remove the sampling functionality and requirement for the random generator; allow the
 * mean to be altered using properties; use {@link java.lang.Math} and remove the {@code FastMath}
 * dependency.
 *
 * @see <a href="http://en.wikipedia.org/wiki/Poisson_distribution">Poisson distribution
 *      (Wikipedia)</a>
 * @see <a href="http://mathworld.wolfram.com/PoissonDistribution.html">Poisson distribution
 *      (MathWorld)</a>
 */
public final class PoissonDistribution {
  /** Mean of the distribution. */
  private double mean;

  /**
   * Creates a new Poisson distribution with specified mean.
   *
   * @param mean the Poisson mean
   * @throws IllegalArgumentException if {@code mean <= 0}.
   */
  public PoissonDistribution(double mean) {
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
   * @throws IllegalArgumentException if {@code mean <= 0}.
   */
  public void setMean(double mean) {
    ValidationUtils.checkStrictlyPositive(mean, "Mean");
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
    return logProbability == Double.NEGATIVE_INFINITY ? 0 : Math.exp(logProbability);
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
    if (x < 0 || x == Integer.MAX_VALUE) {
      return Double.NEGATIVE_INFINITY;
    }
    if (x == 0) {
      return -mean;
    }
    return -SaddlePointExpansionCopy.getStirlingError(x)
        - SaddlePointExpansionCopy.getDeviancePart(x, mean) - 0.5 * Math.log(MathUtils.TWO_PI)
        - 0.5 * Math.log(x);
  }

  /**
   * For a random variable {@code X} whose values are distributed according to this distribution,
   * this method returns {@code P(X <= x)}. In other words, this method represents the (cumulative)
   * distribution function (CDF) for this distribution.
   *
   * @param x the point at which the CDF is evaluated
   * @return the probability that a random variable with this distribution takes a value less than
   *         or equal to {@code x}
   */
  public double cumulativeProbability(int x) {
    if (x < 0) {
      return 0;
    }
    if (x == Integer.MAX_VALUE) {
      return 1;
    }
    return Gamma.regularizedGammaQ((double) x + 1, mean, 1e-12, Integer.MAX_VALUE);
  }

  /**
   * Computes the quantile function of this distribution. For a random variable {@code X}
   * distributed according to this distribution, the returned value is:
   *
   * <ul> <li><code>inf{x in Z | P(X<=x) >= p}</code> for {@code 0 < p <= 1},</li> <li><code>inf{x
   * in Z | P(X<=x) > 0}</code> for {@code p = 0}.</li> </ul>
   *
   * <p>If the result exceeds the range of the data type {@code int}, then {@code 0} or
   * {@code Integer.MAX_VALUE} is returned.
   *
   * @param probability the cumulative probability
   * @return the smallest {@code p}-quantile of this distribution
   * @throws IllegalArgumentException if {@code p < 0} or {@code p > 1}
   */
  public int inverseCumulativeProbability(final double probability) {
    ValidationUtils.checkArgument(probability >= 0 && probability <= 1, "Invalid probability: %f",
        probability);
    if (probability == 0.0) {
      return 0;
    }
    if (probability == 1.0) {
      return Integer.MAX_VALUE;
    }

    // Use the one-sided Chebyshev inequality to narrow the bracket.
    final double mu = mean;
    final double sigma = Math.sqrt(mean);

    double range = Math.sqrt((1.0 - probability) / probability);
    final double tmp = mu - range * sigma;

    // Using -1 ensures cumulativeProbability(lower) < p, which
    // is required for the solving step.
    final int lower = (tmp > -1) ? ((int) Math.ceil(tmp)) - 1 : -1;

    range = 1.0 / range;
    final int upper = (int) Math.floor(mu + range * sigma);

    return solveInverseCumulativeProbability(probability, lower, upper);
  }

  /**
   * This is a utility function used by {@link #inverseCumulativeProbability(double)}. It assumes
   * {@code 0 < p < 1} and that the inverse cumulative probability lies in the bracket {@code
   * (lower, upper]}. The implementation does simple bisection to find the smallest
   * {@code p}-quantile <code>inf{x in Z | P(X<=x) >= p}</code>.
   *
   * @param probability the cumulative probability
   * @param lower a value satisfying {@code cumulativeProbability(lower) < p}
   * @param upper a value satisfying {@code p <= cumulativeProbability(upper)}
   * @return the smallest {@code p}-quantile of this distribution
   */
  private int solveInverseCumulativeProbability(final double probability, int lower, int upper) {
    while (lower + 1 < upper) {
      final int mid = (lower + upper) >> 1;
      final double pm = cumulativeProbability(mid);
      if (pm >= probability) {
        upper = mid;
      } else {
        lower = mid;
      }
    }
    return upper;
  }
}
