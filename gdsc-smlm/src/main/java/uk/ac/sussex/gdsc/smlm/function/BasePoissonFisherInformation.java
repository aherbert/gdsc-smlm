/*-
 * #%L
 * Genome Damage and Stability Centre SMLM Package
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2025 Alex Herbert
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
 * Base class for any probability distribution that can calculate the Fisher information for a
 * Poisson distributed mean.
 *
 * <p><a
 * href="https://en.wikipedia.org/wiki/Poisson_distribution">https://en.wikipedia.org/wiki/Poisson_distribution</a>
 */
public abstract class BasePoissonFisherInformation implements FisherInformation {
  /**
   * The lowest value for the mean that can be computed. This is the lowest value where the
   * reciprocal is not infinity
   */
  public static final double MIN_MEAN = Double.longBitsToDouble(0x4000000000001L);

  @Override
  public boolean isValid(double theta) {
    return theta >= MIN_MEAN;
  }

  /**
   * Gets the alpha scale of the Poisson Fisher information. This is the Fisher information relative
   * to the Fisher information of a pure Poisson distribution:
   *
   * <pre>
   * alpha = FI / (Poisson FI) = FI * theta
   * </pre>
   *
   * @param theta parameter θ of a distribution that models X
   * @return the alpha scale
   * @throws IllegalArgumentException if the parameter is not in the valid range
   */
  public abstract double getAlpha(double theta);

  /**
   * Create a copy.
   *
   * @return the copy
   */
  public abstract BasePoissonFisherInformation copy();
}
