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

package uk.ac.sussex.gdsc.smlm.filters;

import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.rng.RandomUtils;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.junit5.SpeedTag;
import uk.ac.sussex.gdsc.test.rng.RngUtils;
import uk.ac.sussex.gdsc.test.utils.BaseTimingTask;
import uk.ac.sussex.gdsc.test.utils.TestComplexity;
import uk.ac.sussex.gdsc.test.utils.TestLogUtils;
import uk.ac.sussex.gdsc.test.utils.TestSettings;
import uk.ac.sussex.gdsc.test.utils.TimingResult;
import uk.ac.sussex.gdsc.test.utils.TimingService;
import uk.ac.sussex.gdsc.test.utils.functions.FunctionUtils;

import ij.plugin.filter.GaussianBlur;
import ij.process.FloatProcessor;

import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.sampling.distribution.AhrensDieterExponentialSampler;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.BeforeAll;

import java.awt.Rectangle;
import java.util.logging.Level;
import java.util.logging.Logger;

@SuppressWarnings({"javadoc"})
public class GaussianFilterTest {
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
  public void floatFilterIsSameAsImageJFilter(RandomSeed seed) {
    filter1IsSameAsFilter2(seed, new FloatFilter(false), new ImageJFilter(false), false, 1e-2);
  }

  @SeededTest
  public void floatFilterInternalIsSameAsImageJFilter(RandomSeed seed) {
    filter1IsSameAsFilter2(seed, new FloatFilter(true), new ImageJFilter(true), false, 1e-2);
  }

  @SeededTest
  public void floatFilterIsSameAsDoubleFilter(RandomSeed seed) {
    filter1IsSameAsFilter2(seed, new FloatFilter(false), new DoubleFilter(false), false, 1e-2);
  }

  @SeededTest
  public void floatFilterIsSameAsDoubleFilterWeighted(RandomSeed seed) {
    filter1IsSameAsFilter2(seed, new FloatFilter(false), new DoubleFilter(false), true, 1e-2);
  }

  @SeededTest
  public void dpFloatFilterIsSameAsDoubleFilter(RandomSeed seed) {
    filter1IsSameAsFilter2(seed, new DpFilter(false), new DoubleFilter(false), false, 1e-2);
  }

  @SeededTest
  public void dpFloatFilterIsSameAsDoubleFilterWeighted(RandomSeed seed) {
    filter1IsSameAsFilter2(seed, new DpFilter(false), new DoubleFilter(false), true, 1e-2);
  }

  private void filter1IsSameAsFilter2(RandomSeed seed, GFilter f1, GFilter f2, boolean weighted,
      double tolerance) {
    final UniformRandomProvider rand = RngUtils.create(seed.getSeedAsLong());
    final float[] data = createData(rand, size, size);
    float[] weights = null;
    if (weighted) {
      final AhrensDieterExponentialSampler ed = new AhrensDieterExponentialSampler(rand, 57);

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

      logger.fine(FunctionUtils.getSupplier("%s vs %s w=%b @ %.1f = %g", f1.getName(), f2.getName(),
          weighted, sigma, max));
      Assertions.assertTrue(max < tolerance);
    }
  }

  private class MyTimingTask extends BaseTimingTask {
    GFilter filter;
    float[][] data;
    double sigma;

    public MyTimingTask(GFilter filter, float[][] data, double sigma) {
      super(filter.getName() + " " + sigma);
      this.filter = filter;
      this.data = data;
      this.sigma = sigma;
    }

    @Override
    public int getSize() {
      return data.length;
    }

    @Override
    public Object getData(int index) {
      return data[index].clone();
    }

    @Override
    public Object run(Object data) {
      final float[] d = (float[]) data;
      return filter.run(d, sigma);
    }
  }

  @SpeedTag
  @SeededTest
  public void floatFilterIsFasterThanDoubleFilter(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final UniformRandomProvider rg = RngUtils.create(seed.getSeedAsLong());

    final float[][] data = new float[10][];
    for (int i = 0; i < data.length; i++) {
      data[i] = createData(rg, size, size);
    }

    final TimingService ts = new TimingService();
    for (final double sigma : sigmas) {
      ts.execute(new MyTimingTask(new FloatFilter(false), data, sigma));
      ts.execute(new MyTimingTask(new DpFilter(false), data, sigma));
      ts.execute(new MyTimingTask(new DoubleFilter(false), data, sigma));
    }
    final int size = ts.getSize();
    ts.repeat();
    if (logger.isLoggable(Level.INFO)) {
      logger.info(ts.getReport(size));
    }
    final int n = size / sigmas.length;
    for (int i = 0, j = size; i < sigmas.length; i++, j += n) {
      for (int k = 1; k < n; k++) {
        final TimingResult slow = ts.get(j + k);
        final TimingResult fast = ts.get(j);
        logger.log(TestLogUtils.getTimingRecord(slow, fast));
      }
    }
  }

  @SpeedTag
  @SeededTest
  public void floatFilterInternalIsFasterThanDoubleFilterInternal(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.HIGH));

    final UniformRandomProvider rg = RngUtils.create(seed.getSeedAsLong());

    final float[][] data = new float[10][];
    for (int i = 0; i < data.length; i++) {
      data[i] = createData(rg, size, size);
    }

    final TimingService ts = new TimingService();
    for (final double sigma : sigmas) {
      ts.execute(new MyTimingTask(new FloatFilter(true), data, sigma));
      ts.execute(new MyTimingTask(new DpFilter(false), data, sigma));
      ts.execute(new MyTimingTask(new DoubleFilter(true), data, sigma));
    }
    final int size = ts.getSize();
    ts.repeat();
    if (logger.isLoggable(Level.INFO)) {
      logger.info(ts.getReport(size));
    }
    final int n = size / sigmas.length;
    for (int i = 0, j = size; i < sigmas.length; i++, j += n) {
      for (int k = 1; k < n; k++) {
        final TimingResult slow = ts.get(j + k);
        final TimingResult fast = ts.get(j);
        logger.log(TestLogUtils.getTimingRecord(slow, fast));
      }
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
