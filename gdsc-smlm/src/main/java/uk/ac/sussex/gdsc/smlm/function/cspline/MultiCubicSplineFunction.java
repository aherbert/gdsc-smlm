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

package uk.ac.sussex.gdsc.smlm.function.cspline;

import java.util.Arrays;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.smlm.function.Gradient1Procedure;
import uk.ac.sussex.gdsc.smlm.function.Gradient2Procedure;
import uk.ac.sussex.gdsc.smlm.function.ValueProcedure;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;

/**
 * Represent a cubic spline function for multiple points.
 */
public class MultiCubicSplineFunction extends CubicSplineFunction {
  /** The number of splines to draw. */
  private int numberOfSplines = 1;

  /** The gradient indices, This can be cached for the same n. */
  private int[] gradientIndices;

  /** The n target splines. This is cached to re-use memory */
  private TargetSpline[] targetSplines = new TargetSpline[0];

  /**
   * The working splines for the current evaluation.
   */
  private TargetSpline[] working;

  /**
   * The number of working splines. This is <=numberOfSplines depending on whether the spline is
   * within the target region for each point.
   */
  private int workingCount;

  /**
   * The working splines for the current evaluation of the current y-index.
   */
  private TargetSpline[] workingY;

  /**
   * The number of working splines for the current y-index. This is <=w depending on whether the
   * working spline is within the target region for Y..
   */
  private int workingCountY;

  /**
   * Instantiates a new cubic spline function.
   *
   * @param splineData the spline data
   * @param maxx The maximum x value of the 2-dimensional data
   * @param maxy The maximum y value of the 2-dimensional data
   * @throws IllegalArgumentException If the function does not have an integer grid spacing from the
   *         origin
   */
  public MultiCubicSplineFunction(CubicSplineData splineData, int maxx, int maxy) {
    super(splineData, maxx, maxy);
  }

  /**
   * Instantiates a new cubic spline function.
   *
   * @param splineData the spline data
   * @param maxx The maximum x value of the 2-dimensional data
   * @param maxy The maximum y value of the 2-dimensional data
   * @param cx the x centre of the spline data
   * @param cy the y centre of the spline data
   * @param cz the z centre of the spline data
   * @param scale the scale of the spline data
   * @throws IllegalArgumentException If the function does not have an integer grid spacing from the
   *         origin
   */
  public MultiCubicSplineFunction(CubicSplineData splineData, int maxx, int maxy, double cx,
      double cy, double cz, int scale) {
    super(splineData, maxx, maxy, cx, cy, cz, scale);
  }

  @Override
  public int getN() {
    return numberOfSplines;
  }

  /**
   * Sets the number of splines to draw.
   *
   * @param n the new number of splines to draw
   * @throws IllegalArgumentException If the number is not strictly positive
   */
  public void setN(int n) {
    if (n < 1) {
      throw new IllegalArgumentException();
    }
    if (n != this.numberOfSplines) {
      gradientIndices = null;
    }
    this.numberOfSplines = n;
  }

  @Override
  public int[] gradientIndices() {
    if (gradientIndices == null) {
      gradientIndices = SimpleArrayUtils.natural(getNumberOfGradients());
    }
    return gradientIndices;
  }

  @Override
  public int getNumberOfGradients() {
    return 1 + 4 * numberOfSplines;
  }

  @Override
  protected void initialise(double[] parameters, int order) {
    tb = parameters[PeakResult.BACKGROUND];
    // Ensure we have enough room
    if (targetSplines.length < numberOfSplines) {
      int index = targetSplines.length;
      targetSplines = Arrays.copyOf(targetSplines, numberOfSplines); // Preserve memory space
      final boolean sp = splines[0][0].isSinglePrecision();
      while (index < numberOfSplines) {
        targetSplines[index++] = (sp) ? new FloatTargetSpline() : new DoubleTargetSpline();
      }
      working = new TargetSpline[numberOfSplines];
      workingY = new TargetSpline[numberOfSplines];
    }
    // Convert the target parameters to spline offset tables
    workingCount = 0;
    for (int i = 0, j = PeakResult.INTENSITY; i < numberOfSplines; i++) {
      final double tI = parameters[j++];
      final double tX = parameters[j++];
      final double tY = parameters[j++];
      final double tZ = parameters[j++];
      if (targetSplines[i].initialise(i, tI, tX, tY, tZ, order)) {
        working[workingCount++] = targetSplines[i];
      }
    }
  }

  @Override
  public void forEach(ValueProcedure procedure) {
    for (int n = 0; n < workingCount; n++) {
      working[n].reset();
    }

    for (int y = 0; y < maxy; y++) {
      // Get the working targets for this Y
      workingCountY = 0;
      for (int n = 0; n < workingCount; n++) {
        if (working[n].isNextYActive()) {
          workingY[workingCountY++] = working[n];
        }
      }

      if (workingCountY == 0) {
        for (int x = 0; x < maxx; x++) {
          procedure.execute(tb);
        }
      } else {
        for (int x = 0; x < maxx; x++) {
          double intensity = tb;
          for (int n = 0; n < workingCountY; n++) {
            intensity += workingY[n].value(x);
          }
          procedure.execute(intensity);
        }
      }
    }
  }

  @Override
  public void forEach(Gradient1Procedure procedure) {
    final double[] duda = new double[getNumberOfGradients()];
    duda[0] = 1.0;

    for (int n = 0; n < workingCount; n++) {
      working[n].reset();
    }

    for (int y = 0; y < maxy; y++) {
      // Get the working targets for this Y
      workingCountY = 0;
      for (int n = 0; n < workingCount; n++) {
        if (working[n].isNextYActive(duda)) {
          workingY[workingCountY++] = working[n];
        }
      }

      if (workingCountY == 0) {
        for (int x = 0; x < maxx; x++) {
          procedure.execute(tb, duda);
        }
      } else {
        for (int x = 0; x < maxx; x++) {
          double intensity = tb;
          for (int n = 0; n < workingCountY; n++) {
            intensity += workingY[n].value(x, duda);
          }
          procedure.execute(intensity, duda);
        }
      }
    }
  }

  @Override
  public void forEach(Gradient2Procedure procedure) {
    final double[] duda = new double[getNumberOfGradients()];
    final double[] d2uda2 = new double[getNumberOfGradients()];
    duda[0] = 1.0;

    for (int n = 0; n < workingCount; n++) {
      working[n].reset();
    }

    for (int y = 0; y < maxy; y++) {
      // Get the working targets for this Y
      workingCountY = 0;
      for (int n = 0; n < workingCount; n++) {
        if (working[n].isNextYActive(duda, d2uda2)) {
          workingY[workingCountY++] = working[n];
        }
      }

      if (workingCountY == 0) {
        for (int x = 0; x < maxx; x++) {
          procedure.execute(tb, duda, d2uda2);
        }
      } else {
        for (int x = 0; x < maxx; x++) {
          double intensity = tb;
          for (int n = 0; n < workingCountY; n++) {
            intensity += workingY[n].value(x, duda, d2uda2);
          }
          procedure.execute(intensity, duda, d2uda2);
        }
      }
    }
  }

  @Override
  public boolean isNodeBoundary(int gradientIndex) {
    final int parameterIndex = gradientIndices()[gradientIndex];
    if (parameterIndex == BACKGROUND) {
      return false;
    }

    final int dimension = (parameterIndex - 1) % PARAMETERS_PER_PEAK;
    if (dimension == 0) {
      return false; // Signal
    }

    final int peak = getPeak(parameterIndex);
    for (int n = 0; n < workingCount; n++) {
      if (working[n].id == peak) {
        return working[n].isNodeBoundary(dimension - 1);
      }
    }
    return false;
  }
}
