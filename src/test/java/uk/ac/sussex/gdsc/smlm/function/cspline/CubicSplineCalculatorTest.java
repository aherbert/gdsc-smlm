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

package uk.ac.sussex.gdsc.smlm.function.cspline;

import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import uk.ac.sussex.gdsc.core.math.interpolation.CubicSplinePosition;
import uk.ac.sussex.gdsc.core.math.interpolation.CustomTricubicFunction;
import uk.ac.sussex.gdsc.core.math.interpolation.CustomTricubicFunctionUtils;
import uk.ac.sussex.gdsc.core.math.interpolation.CustomTricubicInterpolatingFunction;
import uk.ac.sussex.gdsc.core.math.interpolation.CustomTricubicInterpolator;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;

@SuppressWarnings({"javadoc"})
public class CubicSplineCalculatorTest {
  @Test
  public void canComputeCoefficientsForDistanceFunction() {
    final double[] exp = new double[64];
    int count = 0;
    for (int k = 0; k < 4; k++) {
      for (int j = 0; j < 4; j++) {
        for (int i = 0; i < 4; i++) {
          exp[count++] = Math.sqrt(i * i + j * j + k * k);
        }
      }
    }
    final CustomTricubicFunction f = CustomTricubicFunctionUtils.create(exp);
    final CubicSplinePosition[] s = new CubicSplinePosition[4];
    for (int i = 0; i < 4; i++) {
      s[i] = new CubicSplinePosition((double) i / 3);
    }
    final double[][][] value = new double[4][4][4];
    final double[] b = new double[64];
    count = 0;
    for (int k = 0; k < 4; k++) {
      for (int j = 0; j < 4; j++) {
        for (int i = 0; i < 4; i++) {
          value[i][j][k] = f.value(s[i], s[j], s[k]);
          b[count++] = value[i][j][k];
        }
      }
    }

    final CubicSplineCalculator calc = new CubicSplineCalculator();
    double[] obs = calc.compute(value);
    Assertions.assertArrayEquals(exp, obs, 1e-6);

    obs = calc.compute(b);
    Assertions.assertArrayEquals(exp, obs, 1e-6);
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

    final double[] exp = new double[64];
    f1.getSplineNode(1, 1,1).getCoefficients(exp);

    final CustomTricubicFunction f = CustomTricubicFunctionUtils.create(exp);
    final CubicSplinePosition[] s = new CubicSplinePosition[4];
    for (int i = 0; i < 4; i++) {
      s[i] = new CubicSplinePosition((double) i / 3);
    }
    final double[][][] value = new double[4][4][4];
    final double[] b = new double[64];
    int count = 0;
    for (int k = 0; k < 4; k++) {
      for (int j = 0; j < 4; j++) {
        for (int i = 0; i < 4; i++) {
          value[i][j][k] = f.value(s[i], s[j], s[k]);
          b[count++] = value[i][j][k];
        }
      }
    }

    final CubicSplineCalculator calc = new CubicSplineCalculator();
    double[] obs = calc.compute(value);
    Assertions.assertArrayEquals(exp, obs, 1e-6);

    obs = calc.compute(b);
    Assertions.assertArrayEquals(exp, obs, 1e-6);
  }

  double[][][] createData(int x, int y, int z, UniformRandomProvider rng) {
    final double[][][] fval = new double[x][y][z];
    // Create a 2D Gaussian
    double sd = 1.0;
    double cx = x / 2.0;
    double cy = y / 2.0;
    if (rng != null) {
      sd += rng.nextDouble() - 0.5;
      cx += rng.nextDouble() - 0.5;
      cy += rng.nextDouble() - 0.5;
    }
    final double[] otherx = new double[x];
    for (int zz = 0; zz < z; zz++) {
      final double s2 = 2 * sd * sd;
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
      sd += 0.1;
      cx += 0.1;
      cy -= 0.05;
    }
    return fval;
  }
}
