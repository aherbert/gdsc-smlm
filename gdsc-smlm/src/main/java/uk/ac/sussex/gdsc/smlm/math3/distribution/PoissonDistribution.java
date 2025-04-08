//@formatter:off
/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
//@formatter:on

package uk.ac.sussex.gdsc.smlm.math3.distribution;

import uk.ac.sussex.gdsc.core.utils.ValidationUtils;

/**
 * Implementation of the Poisson distribution for computation of the probability mass function.
 *
 * <p>This implementation is an adaption of the Apache Commons Math 3 class
 * {@link org.apache.commons.math3.distribution.PoissonDistribution}. The code has been updated to:
 * remove the sampling functionality and requirement for the random generator; allow the mean to be
 * altered using properties; use {@link java.lang.Math} and remove the {@code FastMath} dependency.
 *
 * <p>Note: This class has been superseded by the Apache Commons Statistics implementation.
 * Cumulative probability computations are delegated. The probability and log probability
 * computation are maintained to provide zero memory allocation computation of values for a
 * single-use mean. The mean is specified using the class setters.
 *
 * @see <a href="http://en.wikipedia.org/wiki/Poisson_distribution">Poisson distribution
 *      (Wikipedia)</a>
 * @see <a href="http://mathworld.wolfram.com/PoissonDistribution.html">Poisson distribution
 *      (MathWorld)</a>
 * @see org.apache.commons.statistics.distribution.PoissonDistribution
 */
public final class PoissonDistribution {
  /** 0.5 * ln(2 * pi). Computed to 25-digits precision. */
  private static final double HALF_LOG_TWO_PI = 0.9189385332046727417803297;
  /** Mean of the distribution. */
  private double mean;

  /**
   * Creates a new Poisson distribution with the specified mean.
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
    return Math.exp(logProbability(x));
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
    if (x < 0) {
      return Double.NEGATIVE_INFINITY;
    } else if (x == 0) {
      return -mean;
    }
    return -SaddlePointExpansionCopy.getStirlingError(x)
        - SaddlePointExpansionCopy.getDeviancePart(x, mean) - HALF_LOG_TWO_PI - 0.5 * Math.log(x);
  }
}
