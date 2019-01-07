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

import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;

/**
 * Wraps a set of function values to implement the forEach procedure.
 */
public class PrecomputedGradient1Function extends PrecomputedValueFunction
    implements Gradient1Function {
  /** The gradient indices. */
  protected final int[] gradientIndices;

  /** The first order gradient. */
  protected final double[][] g1;

  /**
   * Instantiates a new pre-computed value function.
   *
   * @param values the pre-computed values
   * @param g1 the first order gradient
   * @throws IllegalArgumentException if the values length does not match the function size
   */
  public PrecomputedGradient1Function(double[] values, double[][] g1) {
    super(values);
    final int numberOfGradients = checkGradient(g1);
    gradientIndices = SimpleArrayUtils.natural(numberOfGradients);
    this.g1 = g1;
  }

  /**
   * Check the gradient has the correct length for the function values.
   *
   * @param g the gradient
   * @return the number of gradients
   */
  protected int checkGradient(double[][] g) {
    if (g == null) {
      throw new IllegalArgumentException("Gradient is null");
    }
    if (g.length != values.length) {
      throw new IllegalArgumentException("Gradient is not same size as values");
    }
    if (g.length == 0) {
      return 0;
    }
    if (g[0] == null) {
      throw new IllegalArgumentException("Gradient[0][] is null");
    }
    final int n = g[0].length;
    for (int i = 1; i < g.length; i++) {
      if (g[i] == null || g[i].length != n) {
        throw new IllegalArgumentException("Gradient[" + i + "][] is incorrect size");
      }
    }
    return n;
  }

  /**
   * Gets a reference to the first order gradients.
   *
   * @return the first order gradients
   */
  public double[][] getGradient1Ref() {
    return g1;
  }

  @Override
  public void initialise(double[] a) {
    // Ignore
  }

  @Override
  public int[] gradientIndices() {
    return gradientIndices;
  }

  @Override
  public int getNumberOfGradients() {
    return gradientIndices.length;
  }

  @Override
  public void initialise1(double[] a) {
    // Ignore
  }

  @Override
  public void forEach(Gradient1Procedure procedure) {
    for (int i = 0; i < values.length; i++) {
      procedure.execute(values[i], g1[i]);
    }
  }
}
