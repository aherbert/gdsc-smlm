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

package uk.ac.sussex.gdsc.smlm.utils;

import org.apache.commons.math3.util.FastMath;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.test.api.TestAssertions;
import uk.ac.sussex.gdsc.test.api.TestHelper;
import uk.ac.sussex.gdsc.test.api.function.DoubleDoubleBiPredicate;

@SuppressWarnings({"javadoc"})
class GaussianKernelTest {
  @Test
  void canGetConversionFactor() {
    final DoubleDoubleBiPredicate predicate = TestHelper.doublesAreClose(1e-10, 0);
    for (int i = 0; i < 5; i++) {
      final double s = 0.33 * (1 << i);

      final double norm = 1.0 / (Math.sqrt(2 * Math.PI) * s);
      final double var2 = 2 * s * s;

      final GaussianKernel k = new GaussianKernel(s);
      final double[] o = k.getGaussianKernel(1, 5, false);
      final double f = k.getConversionFactor(o);
      for (int u = o.length / 2, x = u; x >= 0; x--) {
        final double e = norm * FastMath.exp(-(x - u) * (x - u) / var2);
        TestAssertions.assertTest(e, f * o[x], predicate);
      }
    }
  }

  @Test
  void canComputeGaussianKernelIncScaleIncRange() {
    for (int i = 0; i < 5; i++) {
      final double s = 0.33 * (1 << i);

      final GaussianKernel k = new GaussianKernel(s);

      for (int scale = 1; scale < 16; scale *= 2) {
        for (int range = 3; range < 5; range++) {
          for (final boolean edgeCorrection : new boolean[] {true, false}) {
            final double[] e = GaussianKernel.makeGaussianKernel(s * scale, range, edgeCorrection);
            final double[] o = k.getGaussianKernel(scale, range, edgeCorrection);

            Assertions.assertArrayEquals(e, o);
          }
        }
      }
    }
  }

  @Test
  void canComputeGaussianKernelDecScaleIncRange() {
    for (int i = 0; i < 5; i++) {
      final double s = 0.33 * (1 << i);

      final GaussianKernel k = new GaussianKernel(s);

      for (int scale = 16; scale > 1; scale /= 2) {
        for (int range = 3; range < 5; range++) {
          for (final boolean edgeCorrection : new boolean[] {true, false}) {
            final double[] e = GaussianKernel.makeGaussianKernel(s * scale, range, edgeCorrection);
            final double[] o = k.getGaussianKernel(scale, range, edgeCorrection);

            Assertions.assertArrayEquals(e, o);
          }
        }
      }
    }
  }

  @Test
  void canComputeGaussianKernelIncRangeIncScale() {
    for (int i = 0; i < 5; i++) {
      final double s = 0.33 * (1 << i);

      final GaussianKernel k = new GaussianKernel(s);

      for (int range = 3; range < 5; range++) {
        for (int scale = 1; scale < 16; scale *= 2) {
          for (final boolean edgeCorrection : new boolean[] {true, false}) {
            final double[] e = GaussianKernel.makeGaussianKernel(s * scale, range, edgeCorrection);
            final double[] o = k.getGaussianKernel(scale, range, edgeCorrection);

            Assertions.assertArrayEquals(e, o);
          }
        }
      }
    }
  }

  @Test
  void canComputeGaussianKernelIncRangeDecScale() {
    for (int i = 0; i < 5; i++) {
      final double s = 0.33 * (1 << i);

      final GaussianKernel k = new GaussianKernel(s);

      for (int range = 3; range < 5; range++) {
        for (int scale = 16; scale > 1; scale /= 2) {
          for (final boolean edgeCorrection : new boolean[] {true, false}) {
            final double[] e = GaussianKernel.makeGaussianKernel(s * scale, range, edgeCorrection);
            final double[] o = k.getGaussianKernel(scale, range, edgeCorrection);

            Assertions.assertArrayEquals(e, o);
          }
        }
      }
    }
  }

  @Test
  void canComputeDownscaleGaussianKernelIncScaleIncRange() {
    final DoubleDoubleBiPredicate predicate = TestHelper.doublesAreClose(1e-10, 0);
    for (int i = 0; i < 5; i++) {
      final double s = 0.33 * (1 << i);

      final GaussianKernel k = new GaussianKernel(s);

      for (int scale = 1; scale <= 5; scale++) {
        for (int range = 3; range < 5; range++) {
          for (final boolean edgeCorrection : new boolean[] {true, false}) {
            final double[] e = GaussianKernel.makeGaussianKernel(s / scale, range, edgeCorrection);
            final double[] o = k.getDownscaleGaussianKernel(scale, range, edgeCorrection);

            if (MathUtils.isPow2(scale)) {
              Assertions.assertArrayEquals(e, o);
            } else {
              TestAssertions.assertArrayTest(e, o, predicate);
            }
          }
        }
      }
    }
  }

  @Test
  void canComputeDownscaleGaussianKernelDecScaleIncRange() {
    final DoubleDoubleBiPredicate predicate = TestHelper.doublesAreClose(1e-10, 0);
    for (int i = 0; i < 5; i++) {
      final double s = 0.33 * (1 << i);

      final GaussianKernel k = new GaussianKernel(s);

      for (int scale = 5; scale >= 1; scale--) {
        for (int range = 3; range < 5; range++) {
          for (final boolean edgeCorrection : new boolean[] {true, false}) {
            final double[] e = GaussianKernel.makeGaussianKernel(s / scale, range, edgeCorrection);
            final double[] o = k.getDownscaleGaussianKernel(scale, range, edgeCorrection);

            if (MathUtils.isPow2(scale)) {
              Assertions.assertArrayEquals(e, o);
            } else {
              TestAssertions.assertArrayTest(e, o, predicate);
            }
          }
        }
      }
    }
  }

  @Test
  void canComputeDownscaleGaussianKernelIncRangeIncScale() {
    final DoubleDoubleBiPredicate predicate = TestHelper.doublesAreClose(1e-10, 0);
    for (int i = 0; i < 5; i++) {
      final double s = 0.33 * (1 << i);

      final GaussianKernel k = new GaussianKernel(s);

      for (int range = 3; range < 5; range++) {
        for (int scale = 1; scale <= 5; scale++) {
          for (final boolean edgeCorrection : new boolean[] {true, false}) {
            final double[] e = GaussianKernel.makeGaussianKernel(s / scale, range, edgeCorrection);
            final double[] o = k.getDownscaleGaussianKernel(scale, range, edgeCorrection);

            if (MathUtils.isPow2(scale)) {
              Assertions.assertArrayEquals(e, o);
            } else {
              TestAssertions.assertArrayTest(e, o, predicate);
            }
          }
        }
      }
    }
  }

  @Test
  void canComputeDownscaleGaussianKernelIncRangeDecScale() {
    final DoubleDoubleBiPredicate predicate = TestHelper.doublesAreClose(1e-10, 0);
    for (int i = 0; i < 5; i++) {
      final double s = 0.33 * (1 << i);

      final GaussianKernel k = new GaussianKernel(s);

      for (int range = 3; range < 5; range++) {
        for (int scale = 5; scale >= 1; scale--) {
          for (final boolean edgeCorrection : new boolean[] {true, false}) {
            final double[] e = GaussianKernel.makeGaussianKernel(s / scale, range, edgeCorrection);
            final double[] o = k.getDownscaleGaussianKernel(scale, range, edgeCorrection);

            if (MathUtils.isPow2(scale)) {
              Assertions.assertArrayEquals(e, o);
            } else {
              TestAssertions.assertArrayTest(e, o, predicate);
            }
          }
        }
      }
    }
  }
}
