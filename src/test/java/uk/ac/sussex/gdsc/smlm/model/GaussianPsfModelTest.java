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

package uk.ac.sussex.gdsc.smlm.model;

import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.smlm.function.gaussian.AstigmatismZModel;
import uk.ac.sussex.gdsc.smlm.function.gaussian.HoltzerAstigmatismZModel;

import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import java.util.Arrays;
import java.util.logging.Logger;

@SuppressWarnings({"javadoc"})
public class GaussianPsfModelTest {
  @SuppressWarnings("unused")
  private static Logger logger;

  @BeforeAll
  public static void beforeAll() {
    logger = Logger.getLogger(GaussianPsfModelTest.class.getName());
  }

  @AfterAll
  public static void afterAll() {
    logger = null;
  }

  @Test
  public void canComputeValue() {
    // Use a reasonable z-depth function from the Smith, et al (2010) paper (page 377)
    final double sx = 1.08;
    final double sy = 1.01;
    final double gamma = 0.389;
    final double d = 0.531;
    final double Ax = -0.0708;
    final double Bx = -0.073;
    final double Ay = 0.164;
    final double By = 0.0417;
    final AstigmatismZModel zModel =
        HoltzerAstigmatismZModel.create(sx, sy, gamma, d, Ax, Bx, Ay, By);

    final int maxx = 21;
    final int maxy = 21;
    final double[] e = new double[maxx * maxy];
    final double[] o = new double[maxx * maxy];
    final double[] o2 = new double[maxx * maxy];
    final double[][] g = new double[maxx * maxy][];

    final PsfModel psf = new GaussianPsfModel(zModel);

    final double c = maxx * 0.5;
    for (int i = -1; i <= 1; i++) {
      final double x0 = c + i * 0.33;
      for (int j = -1; j <= 1; j++) {
        final double x1 = c + j * 0.33;
        for (int k = -1; k <= 1; k++) {
          final double x2 = k * 0.33;
          Arrays.fill(e, 0);
          psf.create3D(e, maxx, maxy, 1, x0, x1, x2, null);
          psf.getValue(maxx, maxy, x0, x1, x2, o);
          Assertions.assertEquals(1, MathUtils.sum(o), 1e-3);
          psf.getValueAndGradient(maxx, maxy, x0, x1, x2, o2, g);
          for (int ii = 0; ii < e.length; ii++) {
            if (e[ii] == 0) {
              // Same insertion into the blank data region
              Assertions.assertTrue(o[ii] == 0);

              // Only check where there is a reasonable amount of signal
            } else if (e[ii] > 1e-8) {
              final double error = DoubleEquality.relativeError(e[ii], o[ii]);
              // logger.fine(FunctionUtils.getSupplier("[%d,%d] %g == %g %g", ii/maxx, ii%maxx,
              // e[ii], o[ii], error);
              // We expect a small error since the ErfGaussian2DFunction uses a
              // fast approximation of the Erf(..) (the error function). The PSFModel
              // uses the Apache commons implementation.
              if (error > 1e-8) {
                Assertions.fail(String.format("[%d] %s != %s  error = %f", ii,
                    Double.toString(e[ii]), Double.toString(o[ii]), error));
              }
            }
          }
        }
      }
    }
  }

  @Test
  public void canComputeValueAndGradient() {
    // Use a reasonable z-depth function from the Smith, et al (2010) paper (page 377)
    final double sx = 1.08;
    final double sy = 1.01;
    final double gamma = 0.389;
    final double d = 0.531;
    final double Ax = -0.0708;
    final double Bx = -0.073;
    final double Ay = 0.164;
    final double By = 0.0417;
    final AstigmatismZModel zModel =
        HoltzerAstigmatismZModel.create(sx, sy, gamma, d, Ax, Bx, Ay, By);

    // zModel = new NullAstigmatismZModel(1.2, 1.2);

    final int maxx = 21;
    final int maxy = 21;
    final double[] o = new double[maxx * maxy];
    final double[] o2 = new double[maxx * maxy];
    final double[] o3 = new double[maxx * maxy];
    final double[][] g = new double[maxx * maxy][];
    final double[][] g2 = new double[maxx * maxy][3];
    final double[] dx = new double[] {1e-4, 1e-4, 1e-4};

    final PsfModel psf = new GaussianPsfModel(zModel);

    final double c = maxx * 0.5;
    for (int i = -1; i <= 1; i++) {
      final double x0 = c + i * 0.33;
      for (int j = -1; j <= 1; j++) {
        final double x1 = c + j * 0.33;
        for (int k = -1; k <= 1; k++) {
          final double x2 = k * 0.33;
          psf.getValue(maxx, maxy, x0, x1, x2, o);
          psf.getValueAndGradient(maxx, maxy, x0, x1, x2, o2, g);
          Assertions.assertArrayEquals(o, o2);

          // Check verses default manual gradient
          psf.computeValueAndGradient(maxx, maxy, x0, x1, x2, o3, g2, dx);
          Assertions.assertArrayEquals(o, o3);

          for (int ii = 0; ii < o.length; ii++) {
            if (o[ii] > 1e-4) {
              for (int l = 0; l < 3; l++) {
                final double error = DoubleEquality.relativeError(g[ii][l], g2[ii][l]);
                // logger.fine(FunctionUtils.getSupplier("[%d,%d] %g == %g %g", ii, l, g[ii][l],
                // g2[ii][l], error);
                if (error > 5e-3) {
                  Assertions.fail(String.format("[%d] %s != %s  error = %f", ii,
                      Double.toString(g[ii][l]), Double.toString(g2[ii][l]), error));
                }
              }
            }
          }
        }
      }
    }
  }

  @Test
  public void canComputeValueWithDifferentRange() {
    // Use a reasonable z-depth function from the Smith, et al (2010) paper (page 377)
    final double sx = 1.08;
    final double sy = 1.01;
    final double gamma = 0.389;
    final double d = 0.531;
    final double Ax = -0.0708;
    final double Bx = -0.073;
    final double Ay = 0.164;
    final double By = 0.0417;
    final AstigmatismZModel zModel =
        HoltzerAstigmatismZModel.create(sx, sy, gamma, d, Ax, Bx, Ay, By);

    final int maxx = 21;
    final int maxy = 21;
    final double[] o = new double[maxx * maxy];
    final double[] o2 = new double[maxx * maxy];

    final PsfModel psf = new GaussianPsfModel(zModel);
    final GaussianPsfModel psf2 = new GaussianPsfModel(zModel);
    psf2.setRange(40);

    final double c = maxx * 0.5;
    for (int i = -1; i <= 1; i++) {
      final double x0 = c + i * 0.33;
      for (int j = -1; j <= 1; j++) {
        final double x1 = c + j * 0.33;
        for (int k = -1; k <= 1; k++) {
          final double x2 = k * 0.33;
          psf.getValue(maxx, maxy, x0, x1, x2, o);
          psf2.getValue(maxx, maxy, x0, x1, x2, o2);
          int extra = 0;
          int zero = 0;
          for (int ii = 0; ii < o.length; ii++) {
            if (o[ii] == 0) {
              // PSF has a larger range
              zero++;
              if (o2[ii] != 0) {
                extra++;
              }
            } else {
              Assertions.assertEquals(o[ii], o2[ii]);
            }
          }
          Assertions.assertFalse(0 == extra, "No extra");
          final double f = (double) extra / zero;
          // logger.fine(FunctionUtils.getSupplier("Extra %d/%d (%g)", extra, zero, f);
          Assertions.assertTrue(f > 0.4, "Not much extra");
        }
      }
    }
  }
}
