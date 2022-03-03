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

import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import java.util.Arrays;

/**
 * Allow optimisation using Apache Commons Math 3 Optimiser.
 */
public abstract class OptimiserFunction {
  // TODO - Make these private and optimise data access using DoubleArrayList.elements()

  /** The x. */
  protected DoubleArrayList x;
  /** The y. */
  protected DoubleArrayList y;

  /**
   * Adds the point.
   *
   * @param x the x
   * @param y the y
   */
  public void addPoint(double x, double y) {
    if (this.x == null) {
      this.x = new DoubleArrayList();
      this.y = new DoubleArrayList();
    }
    this.x.add(x);
    this.y.add(y);
  }

  /**
   * Adds the data.
   *
   * @param x the x
   * @param y the y
   */
  public void addData(double[] x, double[] y) {
    this.x = DoubleArrayList.wrap(x);
    this.y = DoubleArrayList.wrap(y);
  }

  /**
   * Gets the x data.
   *
   * @return the x
   */
  public double[] getX() {
    return x.toDoubleArray();
  }

  /**
   * Gets the y data.
   *
   * @return the y
   */
  public double[] getY() {
    return y.toDoubleArray();
  }

  /**
   * Gets the weights array. This is an array filled with ones.
   *
   * @return the weights
   */
  public double[] getWeights() {
    final double[] w = new double[y.size()];
    Arrays.fill(w, 1);
    return w;
  }

  /**
   * Get the size.
   *
   * @return the size
   */
  public int size() {
    return x.size();
  }
}
