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

import org.apache.commons.math3.analysis.MultivariateVectorFunction;

/**
 * Wrap the NonLinearFunction to allow use with the Apache Commons Math library.
 */
public class MultivariateVectorFunctionWrapper extends NonLinearFunctionWrapper
    implements MultivariateVectorFunction {
  /**
   * Instantiates a new multivariate vector function wrapper.
   *
   * @param fun The function
   * @param a The parameters
   * @param n The number of data points to evaluate
   */
  public MultivariateVectorFunctionWrapper(NonLinearFunction fun, double[] a, int n) {
    super(fun, a, n);
  }

  /** {@inheritDoc} */
  @Override
  public double[] value(double[] point) {
    return computeValues(point);
  }
}
