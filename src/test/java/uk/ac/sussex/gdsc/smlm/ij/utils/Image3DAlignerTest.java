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
import uk.ac.sussex.gdsc.smlm.function.gaussian.QuadraticAstigmatismZModel;
import uk.ac.sussex.gdsc.test.utils.functions.FunctionUtils;

@SuppressWarnings({"javadoc"})
public class Image3DAlignerTest {
  private static Logger logger;

  @BeforeAll
  public static void beforeAll() {
    logger = Logger.getLogger(Image3DAlignerTest.class.getName());
  }

  @AfterAll
  public static void afterAll() {
    logger = null;
  }

  // TODO - Make this test the StackAligner with sub-pixel accuracy and non power of 2 images

  static final double gamma = 2.5;
  static final int zDepth = 5;
  protected QuadraticAstigmatismZModel zModel = new QuadraticAstigmatismZModel(gamma, zDepth);

  private FloatImage3D createData(int x, int y, int z, double cx, double cy, double cz) {
    final Gaussian2DFunction f =
        GaussianFunctionFactory.create2D(1, x, y, GaussianFunctionFactory.FIT_ASTIGMATISM, zModel);
    final int length = x * y;
    final float[] data = new float[z * length];
    final double[] a = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
    a[Gaussian2DFunction.SIGNAL] = 1;
    a[Gaussian2DFunction.X_POSITION] = cx;
    a[Gaussian2DFunction.Y_POSITION] = cy;
    a[Gaussian2DFunction.X_SD] = 1;
    a[Gaussian2DFunction.Y_SD] = 1;
    final StandardFloatValueProcedure p = new StandardFloatValueProcedure();
    for (int zz = 0; zz < z; zz++) {
      final double dz = zz - cz;
      // if (zz == 0 || zz == z - 1)
      // logger.fine(FunctionUtils.getSupplier("%f %f %f", dz, zModel.getSx(dz), zModel.getSy(dz));
      a[Gaussian2DFunction.Z_POSITION] = dz;
      p.getValues(f, a, data, zz * length);
    }
    return new FloatImage3D(x, y, z, data);
  }

  @Test
  public void canCorrelatePow2Image() {
    canCorrelate(16, 16, 32, false);
  }

  @Test
  public void canCorrelateNonPow2Image() {
    canCorrelate(15, 17, 29, false);
  }

  @Test
  public void canCorrelatePow2ImageUsingImageJImage() {
    canCorrelate(16, 16, 32, true);
  }

  @Test
  public void canCorrelateNonPow2ImageUsingImageJImage() {
    canCorrelate(15, 17, 29, true);
  }

  private void canCorrelate(int maxx, int maxy, int maxz, boolean ijMode) {
    // Not as much information in the z dimension due to axila imaging.
    // Q. How to address in real PSF alignment? Perhaps the z-dimension should be
    // sampled N-times more than the XY dimension.

    final double cx = (maxx - 1) / 2.0;
    final double cy = (maxy - 1) / 2.0;
    final double cz = (maxz - 1) / 2.0;

    final double[] shift = new double[] {0, 1, 1.5, 2};

    final Image3D reference = createData(maxx, maxy, maxz, cx, cy, cz);
    final Image3DAligner a = new Image3DAligner();
    if (ijMode) {
      a.setReference(reference.getImageStack());
    } else {
      a.setReference(reference);
    }

    for (final double sz : shift) {
      for (final double sy : shift) {
        for (final double sx : shift) {
          canCorrelate(a, maxx, maxy, maxz, cx, cy, cz, cx + sx, cy + sy, cz + sz, 0.25, 0, 1e-2, 1,
              ijMode);
          canCorrelate(a, maxx, maxy, maxz, cx, cy, cz, cx + sx, cy + sy, cz + sz, 0.25, 5, 1e-2,
              0.25, ijMode);
        }
      }
    }
  }

  private void canCorrelate(Image3DAligner a, int maxx, int maxy, int maxz, double cx1, double cy1,
      double cz1, double cx2, double cy2, double cz2, double window, int refinements, double error,
      double tolerance, boolean ijMode) {
    double[] result;
    final Image3D target = createData(maxx, maxy, maxz, cx2, cy2, cz2);

    // Utils.display("Ref", reference.getImageStack());
    // Utils.display("Tar", target.getImageStack());

    a.setEdgeWindow(window);
    // a.setSearchMode(uk.ac.sussex.gdsc.smlm.ij.utils.StackAligner.SearchMode.BINARY);
    // a.setRelativeThreshold(1e-6);

    final double[] e = new double[] {cx1 - cx2, cy1 - cy2, cz1 - cz2};
    final int cx = maxx / 2 - (int) Math.round(e[0]);
    final int cy = maxy / 2 - (int) Math.round(e[1]);
    final int cz = maxz / 2 - (int) Math.round(e[2]);
    final int index = target.getIndex(cx, cy, cz);

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
      result = a.align(target.getImageStack(), refinements, error);
    } else {
      result = a.align(target, refinements, error);
    }
    Image3D correlation = a.getCorrelation();
    if (logger.isLoggable(Level.FINE)) {
      logger.fine(FunctionUtils.getSupplier("e %s %g, o %s", java.util.Arrays.toString(e),
          correlation.get(index), java.util.Arrays.toString(result)));
    }

    for (int i = 0; i < 3; i++) {
      Assertions.assertEquals(e[i], result[i], tolerance);
    }
  }
}
