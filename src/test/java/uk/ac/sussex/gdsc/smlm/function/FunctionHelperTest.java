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

import java.util.Arrays;
import java.util.logging.Logger;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultHelper;
import uk.ac.sussex.gdsc.test.utils.functions.FunctionUtils;

@SuppressWarnings({"javadoc"})
class FunctionHelperTest {
  @Test
  void canGetMeanValue() {
    final int n = 10;
    final double[] values = SimpleArrayUtils.newArray(n, 1.0, 1.0);
    Assertions.assertEquals(10, FunctionHelper.getMeanValue(values.clone(), 0));
    final double total = sum(values, n);
    Assertions.assertEquals(total / n, FunctionHelper.getMeanValue(values.clone(), 1));
    for (int i = 1; i < n; i++) {
      final double sum = sum(values, i);
      Assertions.assertEquals(sum / i, FunctionHelper.getMeanValue(values.clone(), sum / total));
    }
  }

  private static double sum(double[] values, int top) {
    double sum = 0;
    int count = 0;
    for (int i = values.length; i-- > 0 && count < top; count++) {
      sum += values[i];
    }
    return sum;
  }

  @Test
  void canGetFractionalMeanValue() {
    final int n = 10;
    final double[] values = SimpleArrayUtils.newDoubleArray(n, 1.0);
    Assertions.assertEquals(1, FunctionHelper.getMeanValue(values.clone(), 0));
    for (int i = 1; i < n; i++) {
      final double f = (double) i / n;
      Assertions.assertEquals(1, FunctionHelper.getMeanValue(values.clone(), f));
      Assertions.assertEquals(1, FunctionHelper.getMeanValue(values.clone(), f - 0.5));
    }

    Arrays.fill(values, 5, n, 2);
    // sum = 5*1 + 5*2 = 15
    Assertions.assertEquals(2, FunctionHelper.getMeanValue(values.clone(), 5.0 / 15));
    Assertions.assertEquals(2, FunctionHelper.getMeanValue(values.clone(), 10.0 / 15));
    Assertions.assertEquals(11.0 / 6, FunctionHelper.getMeanValue(values.clone(), 11.0 / 15));
    Assertions.assertEquals(11.5 / 6.5, FunctionHelper.getMeanValue(values.clone(), 11.5 / 15));
  }

  @Test
  void canGetXValue() {
    final int n = 10;
    final double[] values = SimpleArrayUtils.newArray(n, 1.0, 1.0);
    Assertions.assertEquals(0, FunctionHelper.getXValue(values.clone(), 0));
    Assertions.assertEquals(n, FunctionHelper.getXValue(values.clone(), 1));
    final double total = sum(values, n);
    for (int i = 1; i < n; i++) {
      final double sum = sum(values, i);
      Assertions.assertEquals(i, FunctionHelper.getXValue(values.clone(), sum / total));
    }
  }

  @Test
  void canGetFractionalXValue() {
    final int n = 10;
    final double[] values = SimpleArrayUtils.newDoubleArray(n, 1.0);
    Assertions.assertEquals(0, FunctionHelper.getXValue(values.clone(), 0));
    Assertions.assertEquals(n, FunctionHelper.getXValue(values.clone(), 1));
    for (int i = 1; i < n; i++) {
      final double f = (double) i / n;
      Assertions.assertEquals(i, FunctionHelper.getXValue(values.clone(), f));
      Assertions.assertEquals(i - 0.5, FunctionHelper.getXValue(values.clone(), f - 0.05), 1e-8);
    }

    Arrays.fill(values, 5, n, 2);
    // sum = 5*1 + 5*2 = 15
    Assertions.assertEquals(2.5, FunctionHelper.getXValue(values.clone(), 5.0 / 15));
    Assertions.assertEquals(5, FunctionHelper.getXValue(values.clone(), 10.0 / 15));
    Assertions.assertEquals(6, FunctionHelper.getXValue(values.clone(), 11.0 / 15));
    Assertions.assertEquals(6.5, FunctionHelper.getXValue(values.clone(), 11.5 / 15));
  }

  @Test
  void canGetMeanValueForGaussian() {
    final float intensity = 100;
    // Realistic standard deviations.
    // Only test the highest
    final float[] s = {1, 1.2f, 1.5f, 2, 2.5f};
    final int n_1 = s.length - 1;
    // Flag to indicate that all levels should be run and the difference reported
    final boolean debug = false;
    final Logger logger = (debug) ? Logger.getLogger(FunctionHelperTest.class.getName()) : null;
    for (int i = (debug) ? 0 : n_1; i < s.length; i++) {
      final float sx = s[i];
      for (int j = i; j < s.length; j++) {
        final float sy = s[j];
        final int size = 1 + 2 * (int) Math.ceil(Math.max(sx, sy) * 4);
        final float[] a = Gaussian2DPeakResultHelper.createParams(0.f, intensity, size / 2f,
            size / 2f, 0.f, sx, sy, 0);
        final Gaussian2DFunction f = GaussianFunctionFactory.create2D(1, size, size,
            GaussianFunctionFactory.FIT_FREE_CIRCLE, null);
        final double[] values = f.computeValues(SimpleArrayUtils.toDouble(a));
        // ImagePlus imp = new ImagePlus("gauss", new FloatProcessor(size, size, values));
        // double cx = size / 2.;
        // Shape shape = new Ellipse2D.Double(cx - sx, cx - sy, 2 * sx, 2 * sy);
        // imp.setRoi(new ShapeRoi(shape));
        // IJ.save(imp, "/Users/ah403/1.tif");
        final double scale = MathUtils.sum(values) / intensity;
        for (int range = 1; range <= 3; range++) {
          final double e = Gaussian2DPeakResultHelper.getMeanSignalUsingR(intensity, sx, sy, range);
          final double o = FunctionHelper.getMeanValue(values.clone(),
              scale * Gaussian2DPeakResultHelper.cumulative2D(range));
          if (debug) {
            logger.fine(FunctionUtils.getSupplier("%g,%g   %d  %g %g  %g", sx, sy, range, e, o,
                DoubleEquality.relativeError(e, o)));
          }
          // Only test the highest
          if (i == n_1) {
            Assertions.assertEquals(e, o, e * 0.025);
          }
        }
      }
    }
  }
}
