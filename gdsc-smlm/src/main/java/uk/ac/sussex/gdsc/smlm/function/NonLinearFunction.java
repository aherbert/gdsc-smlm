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
 * Defines the non-linear fitting function.
 */
public interface NonLinearFunction extends GradientFunction {
  /**
   * The non-linear fitting function. Produce an output predicted value for a given input predictor
   * (x) and partial gradient for each of the coefficients (a).
   *
   * @param x Predictor
   * @param dyda Partial gradient of function with respect to each coefficient identified by
   *        {@link #gradientIndices()} . Note: dyda.length must be &gt;= to gradientIndices().length
   * @return The predicted value y
   */
  double eval(int x, double[] dyda);

  /**
   * The non-linear fitting function. Produce an output predicted value for a given input predictor
   * (x).
   *
   * @param x Predictor
   * @return The predicted value y
   */
  double eval(int x);

  /**
   * The non-linear fitting function. Produce an output predicted value for a given input predictor
   * (x) and partial gradient for each of the coefficients (a).
   *
   * @param x Predictor
   * @param dyda Partial gradient of function with respect to each coefficient identified by
   *        {@link #gradientIndices()} . Note: dyda.length must be &gt;= to gradientIndices().length
   * @param weight The output weight. Equivalent to the expected variance of the predicted value.
   *        This should not be zero to avoid divide by zero error.
   * @return The predicted value y
   * @throws NullPointerException If the output weight argument is null
   * @throws ArrayIndexOutOfBoundsException If the output weight argument is length 0
   */
  double evalw(int x, double[] dyda, double[] weight);

  /**
   * The non-linear fitting function. Produce an output predicted value for a given input predictor
   * (x).
   *
   * @param x Predictor
   * @param weight The output weight. Equivalent to the expected variance of the predicted value.
   *        This should not be zero to avoid divide by zero error.
   * @return The predicted value y
   */
  double evalw(int x, double[] weight);

  /**
   * Check if the function can compute weights.
   *
   * @return True if the {@link #evalw(int, double[], double[])} can compute weights other than 1
   */
  boolean canComputeWeights();
}
