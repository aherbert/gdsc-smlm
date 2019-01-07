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

package uk.ac.sussex.gdsc.smlm.function.cspline;

/**
 * Create a cubic spline function.
 */
public class CubicSplineFunctionFactory {
  /**
   * Instantiates a new cubic spline function to model n points.
   *
   * @param splineData the spline data
   * @param maxx The maximum x value of the 2-dimensional data
   * @param maxy The maximum y value of the 2-dimensional data
   * @param cx the x centre of the spline data
   * @param cy the y centre of the spline data
   * @param cz the z centre of the spline data
   * @param scale the scale of the spline data
   * @param n the number of points
   * @return the cubic spline function
   * @throws IllegalArgumentException If the function does not have an integer grid spacing from the
   *         origin
   */
  public static CubicSplineFunction createCubicSplineFunction(CubicSplineData splineData, int maxx,
      int maxy, double cx, double cy, double cz, int scale, int n) {
    if (n == 1) {
      return new SingleCubicSplineFunction(splineData, maxx, maxy, cx, cy, cz, scale);
    }
    final MultiCubicSplineFunction f =
        new MultiCubicSplineFunction(splineData, maxx, maxy, cx, cy, cz, scale);
    f.setN(n);
    return f;
  }
}
