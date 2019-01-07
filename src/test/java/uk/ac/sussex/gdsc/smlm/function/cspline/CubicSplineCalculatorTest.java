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

package uk.ac.sussex.gdsc.smlm.function.cspline;

import uk.ac.sussex.gdsc.core.math.interpolation.CubicSplinePosition;
import uk.ac.sussex.gdsc.core.math.interpolation.CustomTricubicFunction;
import uk.ac.sussex.gdsc.core.math.interpolation.CustomTricubicInterpolatingFunction;
import uk.ac.sussex.gdsc.core.math.interpolation.CustomTricubicInterpolator;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;

import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

@SuppressWarnings({"javadoc"})
public class CubicSplineCalculatorTest {
  @Test
  public void canComputeCoefficientsForDistanceFunction() {
    final double[] e = new double[64];
    int c = 0;
    for (int k = 0; k < 4; k++) {
      for (int j = 0; j < 4; j++) {
        for (int i = 0; i < 4; i++) {
          e[c++] = Math.sqrt(i * i + j * j + k * k);
        }
      }
    }
    final CustomTricubicFunction f = CustomTricubicFunction.create(e);
    final CubicSplinePosition[] s = new CubicSplinePosition[4];
    for (int i = 0; i < 4; i++) {
      s[i] = new CubicSplinePosition((double) i / 3);
    }
    final double[][][] value = new double[4][4][4];
    final double[] b = new double[64];
    c = 0;
    for (int k = 0; k < 4; k++) {
      for (int j = 0; j < 4; j++) {
        for (int i = 0; i < 4; i++) {
          value[i][j][k] = f.value(s[i], s[j], s[k]);
          b[c++] = value[i][j][k];
        }
      }
    }

    final CubicSplineCalculator calc = new CubicSplineCalculator();
    double[] o = calc.compute(value);
    Assertions.assertArrayEquals(e, o, 1e-6);

    o = calc.compute(b);
    Assertions.assertArrayEquals(e, o, 1e-6);
  }

  @Test
  public void canComputeCoefficientsForGaussianFunction() {
    final int x = 4;
    final int y = 4;
    final int z = 4;
    final double xscale = 1;
    final double yscale = 0.5;
    final double zscale = 2.0;
    final double[] xval = SimpleArrayUtils.newArray(x, 0, xscale);
    final double[] yval = SimpleArrayUtils.newArray(y, 0, yscale);
    final double[] zval = SimpleArrayUtils.newArray(z, 0, zscale);
    final double[][][] fval = createData(x, y, z, null);
    final CustomTricubicInterpolatingFunction f1 =
        new CustomTricubicInterpolator().interpolate(xval, yval, zval, fval);

    final double[] e = f1.getCoefficients(1, 1, 1);

    final CustomTricubicFunction f = CustomTricubicFunction.create(e);
    final CubicSplinePosition[] s = new CubicSplinePosition[4];
    for (int i = 0; i < 4; i++) {
      s[i] = new CubicSplinePosition((double) i / 3);
    }
    final double[][][] value = new double[4][4][4];
    final double[] b = new double[64];
    int c = 0;
    for (int k = 0; k < 4; k++) {
      for (int j = 0; j < 4; j++) {
        for (int i = 0; i < 4; i++) {
          value[i][j][k] = f.value(s[i], s[j], s[k]);
          b[c++] = value[i][j][k];
        }
      }
    }

    final CubicSplineCalculator calc = new CubicSplineCalculator();
    double[] o = calc.compute(value);
    Assertions.assertArrayEquals(e, o, 1e-6);

    o = calc.compute(b);
    Assertions.assertArrayEquals(e, o, 1e-6);
  }

  double[][][] createData(int x, int y, int z, UniformRandomProvider r) {
    final double[][][] fval = new double[x][y][z];
    // Create a 2D Gaussian
    double s = 1.0;
    double cx = x / 2.0;
    double cy = y / 2.0;
    if (r != null) {
      s += r.nextDouble() - 0.5;
      cx += r.nextDouble() - 0.5;
      cy += r.nextDouble() - 0.5;
    }
    final double[] otherx = new double[x];
    for (int zz = 0; zz < z; zz++) {
      final double s2 = 2 * s * s;
      for (int xx = 0; xx < x; xx++) {
        otherx[xx] = MathUtils.pow2(xx - cx) / s2;
      }
      for (int yy = 0; yy < y; yy++) {
        final double othery = MathUtils.pow2(yy - cy) / s2;
        for (int xx = 0; xx < x; xx++) {
          fval[xx][yy][zz] = Math.exp(otherx[xx] + othery);
        }
      }
      // Move Gaussian
      s += 0.1;
      cx += 0.1;
      cy -= 0.05;
    }
    return fval;
  }
}
