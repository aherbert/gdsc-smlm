/*-
 * #%L
 * Genome Damage and Stability Centre SMLM Package
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2023 Alex Herbert
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
 * Base class providing support methods to allow univariate function optimisation.
 */
public class OptimiserFunction {
  // TODO - Make these private

  /** The x.
   * @deprecated use {@link #getX(int)} to access the data values. */
  @Deprecated
  protected DoubleArrayList x;
  /** The y. 
   * @deprecated use {@link #getY(int)} to access the data values. */
  @Deprecated
  protected DoubleArrayList y;

  /**
   * Create an instance.
   */
  protected OptimiserFunction() {
    // Intentionally empty
  }

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
   * Adds the data. All previous data is replaced.
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
   * Gets the x data at the specified index.
   * <p>Warning: This performs no bounds checking and may return invalid values outside of {@code [0, size)}.
   *
   * @param i the index
   * @return the x value
   * @see #size()
   */
  public double getX(int i) {
    return x.elements()[i];
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
   * Gets the y data at the specified index.
   * <p>Warning: This performs no bounds checking and may return invalid values outside of {@code [0, size)}.
   *
   * @param i the index
   * @return the y value
   * @see #size()
   */
  public double getY(int i) {
    return y.elements()[i];
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
