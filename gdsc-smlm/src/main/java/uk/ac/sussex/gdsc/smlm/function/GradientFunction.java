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
 * Defines a function that can compute gradients.
 */
public interface GradientFunction {
  /**
   * Set the predictor coefficients (a) that will be used to predict each value. Allows the function
   * to perform initialisation.
   *
   * @param a An array of coefficients
   */
  void initialise(double[] a);

  /**
   * The function will evaluate the gradient for up to {@code n} parameters where
   * {@code n <= a.length}. This method returns the indices that are evaluated.
   *
   * @return The gradient indices
   */
  int[] gradientIndices();

  /**
   * Gets the number of gradients. The function will evaluate this many partial derivatives.
   *
   * @return the number of gradients
   */
  int getNumberOfGradients();
}
