/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2020 Alex Herbert
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
 * Wraps a set of function values to implement the forEach procedure.
 */
public class PrecomputedGradient2Function extends PrecomputedGradient1Function
    implements Gradient2Function {
  /** The second order gradient. */
  protected final double[][] g2;

  /**
   * Instantiates a new pre-computed value function.
   *
   * @param values the pre-computed values
   * @param g1 the first order gradient
   * @param g2 the second order gradient
   * @throws IllegalArgumentException if the values length does not match the function size
   */
  public PrecomputedGradient2Function(double[] values, double[][] g1, double[][] g2) {
    super(values, g1);
    checkGradient(g2);
    this.g2 = g2;
  }

  /**
   * Gets a reference to the first order gradients.
   *
   * @return the first order gradients
   */
  public double[][] getGradient2Ref() {
    return g2;
  }

  @Override
  public void initialise2(double[] a) {
    // Ignore
  }

  @Override
  public void forEach(Gradient2Procedure procedure) {
    for (int i = 0; i < values.length; i++) {
      procedure.execute(values[i], g1[i], g2[i]);
    }
  }
}
