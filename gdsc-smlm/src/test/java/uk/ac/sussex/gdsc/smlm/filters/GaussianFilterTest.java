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

import ij.plugin.filter.GaussianBlur;
import ij.process.FloatProcessor;
import java.awt.Rectangle;
import java.util.logging.Logger;
import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.sampling.distribution.ContinuousSampler;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;
import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.rng.RandomUtils;
import uk.ac.sussex.gdsc.core.utils.rng.SamplerUtils;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngFactory;
import uk.ac.sussex.gdsc.test.utils.RandomSeed;
import uk.ac.sussex.gdsc.test.utils.TestLogging.TestLevel;
import uk.ac.sussex.gdsc.test.utils.functions.FunctionUtils;

@SuppressWarnings({"javadoc"})
class GaussianFilterTest {
  private static Logger logger;

  @BeforeAll
  public static void beforeAll() {
    logger = Logger.getLogger(GaussianFilterTest.class.getName());
  }

  @AfterAll
  public static void afterAll() {
    logger = null;
  }

  double[] sigmas = new double[] {12.4, 9.3, 5, 3.2, 2.1, 0.5};
  static final int size = 256;

  private abstract static class GFilter {
    boolean internal;
    String name;

    GFilter(String name, boolean internal) {
      this.name = name;
      this.internal = internal;
    }

    String getName() {
      if (internal) {
        return name + " internal";
      }
      return name;
    }

    float[] run(float[] data, double sigma) {
      if (internal) {
        return filterInternal(data, sigma);
      }
      return filter(data, sigma);
    }

    abstract float[] filter(float[] data, double sigma);

    abstract float[] filterInternal(float[] data, double sigma);

    abstract void setWeights(float[] weights);
  }

  private static class ImageJFilter extends GFilter {
    GaussianBlur gf = new GaussianBlur();

    ImageJFilter(boolean internal) {
      super(GaussianBlur.class.getSimpleName(), internal);
    }

    @Override
    float[] filter(float[] data, double sigma) {
      final FloatProcessor fp = new FloatProcessor(size, size, data);
      gf.blurGaussian(fp, sigma, sigma, GaussianFilter.DEFAULT_ACCURACY);
      return data;
    }

    @Override
    float[] filterInternal(float[] data, double sigma) {
      final FloatProcessor fp = new FloatProcessor(size, size, data);
      final int border = GaussianFilter.getBorder(sigma);
      final Rectangle roi = new Rectangle(border, border, size - 2 * border, size - 2 * border);
      fp.setRoi(roi);
      gf.blurGaussian(fp, sigma, sigma, GaussianFilter.DEFAULT_ACCURACY);
      return data;
    }

    @Override
    void setWeights(float[] weights) {
      // Ignored
    }
  }

  private static class FloatFilter extends GFilter {
    GaussianFilter gf = new GaussianFilter();

    FloatFilter(boolean internal) {
      super(GaussianFilter.class.getSimpleName(), internal);
    }

    @Override
    float[] filter(float[] data, double sigma) {
      gf.convolve(data, size, size, sigma);
      return data;
    }

    @Override
    float[] filterInternal(float[] data, double sigma) {
      gf.convolveInternal(data, size, size, sigma);
      return data;
    }

    @Override
    void setWeights(float[] weights) {
      gf.setWeights(weights, size, size);
    }
  }

  private static class DoubleFilter extends GFilter {
    DoubleGaussianFilter gf = new DoubleGaussianFilter();

    DoubleFilter(boolean internal) {
      super(DoubleGaussianFilter.class.getSimpleName(), internal);
    }

    @Override
    float[] filter(float[] data, double sigma) {
      gf.convolve(data, size, size, sigma);
      return data;
    }

    @Override
    float[] filterInternal(float[] data, double sigma) {
      gf.convolveInternal(data, size, size, sigma);
      return data;
    }

    @Override
    void setWeights(float[] weights) {
      gf.setWeights(weights, size, size);
    }
  }

  private static class DpFilter extends GFilter {
    DpGaussianFilter gf = new DpGaussianFilter();

    DpFilter(boolean internal) {
      super(DpGaussianFilter.class.getSimpleName(), internal);
    }

    @Override
    float[] filter(float[] data, double sigma) {
      gf.convolve(data, size, size, sigma);
      return data;
    }

    @Override
    float[] filterInternal(float[] data, double sigma) {
      gf.convolveInternal(data, size, size, sigma);
      return data;
    }

    @Override
    void setWeights(float[] weights) {
      gf.setWeights(weights, size, size);
    }
  }

  @SeededTest
  void floatFilterIsSameAsImageJFilter(RandomSeed seed) {
    filter1IsSameAsFilter2(seed, new FloatFilter(false), new ImageJFilter(false), false, 1e-2);
  }

  @SeededTest
  void floatFilterInternalIsSameAsImageJFilter(RandomSeed seed) {
    filter1IsSameAsFilter2(seed, new FloatFilter(true), new ImageJFilter(true), false, 1e-2);
  }

  @SeededTest
  void floatFilterIsSameAsDoubleFilter(RandomSeed seed) {
    filter1IsSameAsFilter2(seed, new FloatFilter(false), new DoubleFilter(false), false, 1e-2);
  }

  @SeededTest
  void floatFilterIsSameAsDoubleFilterWeighted(RandomSeed seed) {
    filter1IsSameAsFilter2(seed, new FloatFilter(false), new DoubleFilter(false), true, 1e-2);
  }

  @SeededTest
  void dpFloatFilterIsSameAsDoubleFilter(RandomSeed seed) {
    filter1IsSameAsFilter2(seed, new DpFilter(false), new DoubleFilter(false), false, 1e-2);
  }

  @SeededTest
  void dpFloatFilterIsSameAsDoubleFilterWeighted(RandomSeed seed) {
    filter1IsSameAsFilter2(seed, new DpFilter(false), new DoubleFilter(false), true, 1e-2);
  }

  private void filter1IsSameAsFilter2(RandomSeed seed, GFilter f1, GFilter f2, boolean weighted,
      double tolerance) {
    final UniformRandomProvider rand = RngFactory.create(seed.get());
    final float[] data = createData(rand, size, size);
    float[] weights = null;
    if (weighted) {
      final ContinuousSampler ed = SamplerUtils.createExponentialSampler(rand, 57);

      weights = new float[data.length];
      for (int i = 0; i < weights.length; i++) {
        weights[i] = (float) (1.0 / Math.max(0.01, ed.sample()));
      }
      // w[i] = (float) (1.0 / Math.max(0.01, rand.nextGaussian() * 0.2 + 2));
      // w[i] = 0.5f;
      f1.setWeights(weights);
      f2.setWeights(weights);
    }

    for (final double sigma : sigmas) {
      final float[] e = data.clone();
      f2.run(e, sigma);
      final float[] o = data.clone();
      f1.run(o, sigma);

      double max = 0;
      for (int i = 0; i < e.length; i++) {
        final double d = DoubleEquality.relativeError(e[i], o[i]);
        if (max < d) {
          max = d;
        }
      }

      logger.log(TestLevel.TEST_DEBUG, FunctionUtils.getSupplier("%s vs %s w=%b @ %.1f = %g",
          f1.getName(), f2.getName(), weighted, sigma, max));
      Assertions.assertTrue(max < tolerance);
    }
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
