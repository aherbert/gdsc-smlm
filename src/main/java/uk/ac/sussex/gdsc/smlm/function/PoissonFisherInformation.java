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

/**
 * Calculate the Fisher information for a Poisson distribution.
 *
 * <p><a
 * href="https://en.wikipedia.org/wiki/Poisson_distribution">https://en.wikipedia.org/wiki/Poisson_distribution</a>
 */
public class PoissonFisherInformation extends BasePoissonFisherInformation {
  /**
   * {@inheritDoc}
   *
   * <p>The input parameter refers to the mean of the Poisson distribution. The Fisher information
   * is 1/mean.
   */
  @Override
  public double getFisherInformation(double theta) {
    if (theta <= 0) {
      throw new IllegalArgumentException("Poisson mean must be positive");
    }
    return 1.0 / theta;
  }

  /**
   * Gets the Poisson Fisher information.
   *
   * @param theta the poisson mean
   * @return the poisson Fisher information
   */
  public static double getPoissonI(double theta) {
    if (theta <= 0) {
      throw new IllegalArgumentException("Poisson mean must be positive");
    }
    return 1.0 / theta;
  }

  @Override
  public double getAlpha(double theta) {
    return 1;
  }

  @Override
  public PoissonFisherInformation copy() {
    // No state so no need to copy
    return this;
  }
}
