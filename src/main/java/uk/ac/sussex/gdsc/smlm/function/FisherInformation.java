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

/**
 * Calculate the Fisher information: the amount of information that an observable random variable X
 * carries about an unknown parameter θ of a distribution that models X.
 *
 * <pre>
 * I = E [ (d log f(X;θ) / dθi) (d log f(X;θ) / dθj) | θ ]
 * E = Expected value
 * </pre>
 */
public interface FisherInformation {
  /**
   * Gets the fisher information: the amount of information that an observable random variable X
   * carries about an unknown parameter θ of a distribution that models X.
   *
   * @param t parameter θ of a distribution that models X
   * @return the fisher information
   * @throws IllegalArgumentException if the parameter is not in the valid range
   */
  public double getFisherInformation(double t) throws IllegalArgumentException;

  /**
   * Checks if the parameter θ is in a valid range to compute a representable value. <p> If not true
   * then it would be expected that {@link #getFisherInformation(double)} will: throw an exception;
   * compute zero; or compute infinity.
   *
   * @param t parameter θ of a distribution that models X
   * @return true, if a representable value can be computed
   */
  public boolean isValid(double t);
}
