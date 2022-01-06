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

/**
 * Calculate the Fisher information for a Poisson-Gaussian distribution using an approximation of
 * the Poisson (mean=theta) as a Gaussian (u=theta, var=theta). The combined variance of two
 * convolved Gaussians is the sum of the variance.
 *
 * <p>Note that this is equivalent to the reverse by modelling the Gaussian noise s as a Poisson
 * with variance s*s for a Poisson-Poisson Fisher information since the combined variance of two
 * convolved Poissons is the sum of the variance.
 */
public class PoissonGaussianApproximationFisherInformation extends BasePoissonFisherInformation {
  /** The variance of the Gaussian. */
  public final double variance;

  /**
   * Instantiates a new poisson gaussian fisher information.
   *
   * @param sd the standard deviation of the Gaussian
   * @throws IllegalArgumentException If the standard deviation is not strictly positive
   */
  public PoissonGaussianApproximationFisherInformation(double sd) {
    if (!(sd > 0 && sd <= Double.MAX_VALUE)) {
      throw new IllegalArgumentException("Gaussian variance must be strictly positive");
    }
    this.variance = sd * sd;
  }

  /**
   * Copy constructor.
   *
   * @param source the source
   */
  protected PoissonGaussianApproximationFisherInformation(
      PoissonGaussianApproximationFisherInformation source) {
    this.variance = source.variance;
  }

  /**
   * {@inheritDoc}
   *
   * <p>Gets the approximate Poisson-Gaussian Fisher information. Approximate the Poisson as a
   * Gaussian (u=theta, var=theta) and convolve with a Gaussian (u=0,var=s*s). Gaussian-Gaussian
   * convolution: var1 * var2 =&gt; var = var1+var2. The Fisher information of Gaussian mean is
   * 1/variance. The Poisson-Gaussian Fisher information is therefore 1 / (theta + s*s).
   */
  @Override
  public double getFisherInformation(double theta) {
    if (theta <= 0) {
      throw new IllegalArgumentException("Poisson mean must be positive");
    }
    return 1.0 / (theta + variance);
  }

  @Override
  public double getAlpha(double theta) {
    if (theta <= 0) {
      throw new IllegalArgumentException("Poisson mean must be positive");
    }
    return theta / (theta + variance);
  }

  @Override
  public PoissonGaussianApproximationFisherInformation copy() {
    return new PoissonGaussianApproximationFisherInformation(this);
  }
}
