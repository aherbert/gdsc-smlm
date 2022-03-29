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

package uk.ac.sussex.gdsc.smlm.filters;

import ij.plugin.filter.Convolver;
import ij.process.FloatProcessor;
import java.awt.Rectangle;
import java.util.logging.Logger;
import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;
import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.rng.RandomUtils;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngFactory;
import uk.ac.sussex.gdsc.test.utils.RandomSeed;
import uk.ac.sussex.gdsc.test.utils.TestLogging.TestLevel;
import uk.ac.sussex.gdsc.test.utils.functions.FormatSupplier;

@SuppressWarnings({"javadoc"})
class KernelFilterTest {
  private static Logger logger;

  @BeforeAll
  public static void beforeAll() {
    logger = Logger.getLogger(KernelFilterTest.class.getName());
  }

  @AfterAll
  public static void afterAll() {
    logger = null;
  }

  static final int size = 256;
  int[] borders = {0, 1, 2, 3, 5, 10};

  static float[] createKernel(int kw, int kh) {
    // Simple linear ramp
    final float[] k = new float[kw * kh];
    final int cx = kw / 2;
    final int cy = kh / 2;
    for (int y = 0, i = 0; y < kh; y++) {
      final int dy2 = MathUtils.pow2(cy - y);
      for (int x = 0; x < kw; x++) {
        final int dx2 = MathUtils.pow2(cx - x);
        k[i++] = (float) Math.sqrt(dx2 + dy2);
      }
    }
    // Invert
    final float max = k[0];
    for (int i = 0; i < k.length; i++) {
      k[i] = max - k[i];
    }
    return k;
  }

  private abstract static class FilterWrapper {
    final float[] kernel;
    final int kw;
    final int kh;
    String name;

    FilterWrapper(String name, float[] kernel, int kw, int kh) {
      this.name = name + " " + kw + "x" + kh;
      this.kernel = kernel;
      this.kw = kw;
      this.kh = kh;
    }

    String getName() {
      return name;
    }

    abstract float[] filter(float[] data, int border);

    abstract void setWeights(float[] weights);
  }

  private static class ConvolverWrapper extends FilterWrapper {
    Convolver kf = new Convolver();

    ConvolverWrapper(float[] kernel, int kw, int kh) {
      super(Convolver.class.getSimpleName(), kernel, kw, kh);
    }

    @Override
    float[] filter(float[] data, int border) {
      final FloatProcessor fp = new FloatProcessor(size, size, data);
      if (border > 0) {
        final Rectangle roi = new Rectangle(border, border, size - 2 * border, size - 2 * border);
        fp.setRoi(roi);
      }
      kf.convolveFloat(fp, kernel, kw, kh);
      return data;
    }

    @Override
    void setWeights(float[] weights) {
      // Ignored
    }
  }

  private static class KernelFilterWrapper extends FilterWrapper {
    KernelFilter kf = new KernelFilter(kernel, kw, kh);

    KernelFilterWrapper(float[] kernel, int kw, int kh) {
      super(KernelFilterTest.class.getSimpleName(), kernel, kw, kh);
    }

    @Override
    float[] filter(float[] data, int border) {
      kf.convolve(data, size, size, border);
      return data;
    }

    @Override
    void setWeights(float[] weights) {
      kf.setWeights(weights, size, size);
    }
  }

  private static class ZeroKernelFilterWrapper extends FilterWrapper {
    ZeroKernelFilter kf = new ZeroKernelFilter(kernel, kw, kh);

    ZeroKernelFilterWrapper(float[] kernel, int kw, int kh) {
      super(ZeroKernelFilterWrapper.class.getSimpleName(), kernel, kw, kh);
    }

    @Override
    float[] filter(float[] data, int border) {
      kf.convolve(data, size, size, border);
      return data;
    }

    @Override
    void setWeights(float[] weights) {
      kf.setWeights(weights, size, size);
    }
  }

  @SeededTest
  void canRotate180() {
    for (int kw = 1; kw < 3; kw++) {
      for (int kh = 1; kh < 3; kh++) {
        final float[] kernel = createKernel(kw, kh);
        final FloatProcessor fp = new FloatProcessor(kw, kh, kernel.clone());
        fp.flipHorizontal();
        fp.flipVertical();
        KernelFilter.rotate180(kernel);
        Assertions.assertArrayEquals((float[]) fp.getPixels(), kernel);
      }
    }
  }

  @SeededTest
  void kernelFilterIsSameAsImageJFilter(RandomSeed seed) {
    final int kw = 5;
    final int kh = 5;
    final float[] kernel = createKernel(kw, kh);
    filter1IsSameAsFilter2(seed, new KernelFilterWrapper(kernel, kw, kh),
        new ConvolverWrapper(kernel, kw, kh), false, 1e-2);
  }

  @SeededTest
  void zeroKernelFilterIsSameAsImageJFilter(RandomSeed seed) {
    final int kw = 5;
    final int kh = 5;
    final float[] kernel = createKernel(kw, kh);
    filter1IsSameAsFilter2(seed, new ZeroKernelFilterWrapper(kernel, kw, kh),
        new ConvolverWrapper(kernel, kw, kh), true, 1e-2);
  }

  private void filter1IsSameAsFilter2(RandomSeed seed, FilterWrapper f1, FilterWrapper f2,
      boolean internal, double tolerance) {
    final UniformRandomProvider rand = RngFactory.create(seed.get());
    final float[] data = createData(rand, size, size);

    final int testBorder = (internal) ? f1.kw / 2 : 0;
    for (final int border : borders) {
      filter1IsSameAsFilter2(f1, f2, data, border, testBorder, tolerance);
    }
  }

  private static void filter1IsSameAsFilter2(FilterWrapper f1, FilterWrapper f2, float[] data,
      int border, int testBorder, double tolerance) {
    final float[] e = data.clone();
    f2.filter(e, border);
    final float[] o = data.clone();
    f1.filter(o, border);

    double max = 0;
    if (testBorder == 0) {
      for (int i = 0; i < e.length; i++) {
        final double d = DoubleEquality.relativeError(e[i], o[i]);
        if (max < d) {
          max = d;
        }
      }
    } else {
      final int limit = size - testBorder;
      for (int y = testBorder; y < limit; y++) {
        for (int x = testBorder, i = y * size + x; x < limit; x++, i++) {
          final double d = DoubleEquality.relativeError(e[i], o[i]);
          if (max < d) {
            max = d;
          }
        }
      }
    }

    logger.log(TestLevel.TEST_DEBUG,
        FormatSupplier.getSupplier("%s vs %s @ %d = %g", f1.getName(), f2.getName(), border, max));
    Assertions.assertTrue(max < tolerance);
  }

  private static float[] createData(UniformRandomProvider rg, int width, int height) {
    final float[] data = new float[width * height];
    for (int i = data.length; i-- > 0;) {
      data[i] = i;
    }

    RandomUtils.shuffle(data, rg);

    return data;
  }
}
