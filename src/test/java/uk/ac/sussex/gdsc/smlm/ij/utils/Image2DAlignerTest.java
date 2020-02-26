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

package uk.ac.sussex.gdsc.smlm.ij.utils;

import java.util.logging.Level;
import java.util.logging.Logger;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import uk.ac.sussex.gdsc.smlm.function.StandardFloatValueProcedure;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import uk.ac.sussex.gdsc.test.utils.functions.FunctionUtils;

@SuppressWarnings({"javadoc"})
public class Image2DAlignerTest {
  private static Logger logger;

  @BeforeAll
  public static void beforeAll() {
    logger = Logger.getLogger(Image2DAlignerTest.class.getName());
  }

  @AfterAll
  public static void afterAll() {
    logger = null;
  }

  private static FloatImage2D createData(int x, int y, double cx, double cy) {
    final Gaussian2DFunction f = GaussianFunctionFactory.create2D(1, x, y,
        GaussianFunctionFactory.FIT_ERF_FREE_CIRCLE, null);
    final double[] a = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
    a[Gaussian2DFunction.SIGNAL] = 1;
    a[Gaussian2DFunction.X_POSITION] = cx;
    a[Gaussian2DFunction.Y_POSITION] = cy;
    a[Gaussian2DFunction.X_SD] = 1.2;
    a[Gaussian2DFunction.Y_SD] = 1.1;
    final StandardFloatValueProcedure p = new StandardFloatValueProcedure();
    return new FloatImage2D(x, y, p.getValues(f, a));
  }

  @Test
  public void canCorrelateSquareImage() {
    canCorrelate(16, 16, false);
  }

  @Test
  public void canCorrelateNonSquareImage() {
    canCorrelate(16, 32, false);
  }

  @Test
  public void canCorrelateNonPow2Image() {
    canCorrelate(17, 29, false);
  }

  @Test
  public void canCorrelateSquareImageUsingImageJImage() {
    canCorrelate(16, 16, true);
  }

  @Test
  public void canCorrelateNonSquareImageUsingImageJImage() {
    canCorrelate(16, 32, true);
  }

  @Test
  public void canCorrelateNonPow2ImageUsingImageJImage() {
    canCorrelate(17, 29, true);
  }

  private static void canCorrelate(int maxx, int maxy, boolean ijMode) {
    final double cx = (maxx - 1) / 2.0;
    final double cy = (maxy - 1) / 2.0;

    final double[] shift = new double[] {-1, -0.5, 0, 1.5, 2};

    final Image2D reference = createData(maxx, maxy, cx, cy);
    final Image2DAligner a = new Image2DAligner();
    if (ijMode) {
      a.setReference(reference.getImageProcessor());
    } else {
      a.setReference(reference);
    }

    for (final double sy : shift) {
      for (final double sx : shift) {
        canCorrelate(a, maxx, maxy, cx, cy, cx + sx, cy + sy, 0.25, 0, 1, ijMode);
        canCorrelate(a, maxx, maxy, cx, cy, cx + sx, cy + sy, 0.25, 5, 0.25, ijMode);
      }
    }
  }

  private static void canCorrelate(Image2DAligner a, int maxx, int maxy, double cx1, double cy1,
      double cx2, double cy2, double window, int refinements, double tolerance, boolean ijMode) {
    double[] result;
    final Image2D target = createData(maxx, maxy, cx2, cy2);

    // Utils.display("Ref", reference.getImageProcessor());
    // Utils.display("Tar", target.getImageProcessor());

    a.setEdgeWindow(window);
    // a.setSearchMode(uk.ac.sussex.gdsc.smlm.ij.utils.StackAligner.SearchMode.BINARY);
    // a.setRelativeThreshold(1e-6);

    final double[] e = new double[] {cx1 - cx2, cy1 - cy2};
    final int cx = maxx / 2 - (int) Math.round(e[0]);
    final int cy = maxy / 2 - (int) Math.round(e[1]);
    final int index = target.getIndex(cx, cy);

    // Debug the convergence
    // for (int i = 0; i <= refinements; i++)
    // {
    // result = a.align(target.copy(), i, error);
    // c = a.getCorrelation();
    //
    // logger.fine(FunctionUtils.getSupplier("e %s %g, o %s", java.util.Arrays.toString(e),
    // c.get(index),
    // java.util.Arrays.toString(result));
    // }

    // Test
    if (ijMode) {
      result = a.align(target.getImageProcessor(), refinements);
    } else {
      result = a.align(target, refinements);
    }
    Image2D correlation = a.getCorrelation();
    if (logger.isLoggable(Level.FINE)) {
      logger.fine(FunctionUtils.getSupplier("e %s %g, o %s", java.util.Arrays.toString(e),
          correlation.get(index), java.util.Arrays.toString(result)));
    }

    for (int i = 0; i < 2; i++) {
      Assertions.assertEquals(e[i], result[i], tolerance);
    }
  }
}
