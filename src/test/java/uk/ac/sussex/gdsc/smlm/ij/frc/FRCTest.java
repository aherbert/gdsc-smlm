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

package uk.ac.sussex.gdsc.smlm.ij.frc;

import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.FloatEquality;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.rng.GaussianSamplerUtils;
import uk.ac.sussex.gdsc.smlm.ij.results.IJImagePeakResults;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngUtils;
import uk.ac.sussex.gdsc.test.utils.BaseTimingTask;
import uk.ac.sussex.gdsc.test.utils.TestComplexity;
import uk.ac.sussex.gdsc.test.utils.TestLogUtils;
import uk.ac.sussex.gdsc.test.utils.TestSettings;
import uk.ac.sussex.gdsc.test.utils.TimingResult;
import uk.ac.sussex.gdsc.test.utils.TimingService;

import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import org.apache.commons.math3.util.FastMath;
import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.sampling.PermutationSampler;
import org.apache.commons.rng.sampling.distribution.GaussianSampler;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.BeforeAll;

import java.awt.Rectangle;
import java.util.logging.Level;
import java.util.logging.Logger;

@SuppressWarnings({"javadoc"})
public class FRCTest {
  private static Logger logger;

  @BeforeAll
  public static void beforeAll() {
    logger = Logger.getLogger(FRCTest.class.getName());
  }

  @AfterAll
  public static void afterAll() {
    logger = null;
  }

  @SeededTest
  public void canComputeSine() {
    final int steps = 1000;
    final double delta = 2 * Math.PI / steps;
    for (int i = 0; i <= steps; i++) {
      final double a = i * delta;
      final double cosA = Math.cos(a);
      final double e = Math.sin(a);
      final double o = FRC.getSine(a, cosA);
      // logger.fine(FunctionUtils.getSupplier("%f %f ?= %f", a, e, o);
      Assertions.assertTrue(DoubleEquality.almostEqualRelativeOrAbsolute(o, e, 1e-6, 1e-10));
    }
  }

  @SeededTest
  public void canComputeMirrored(RandomSeed seed) {
    // Sample lines through an image to create a structure.
    final int size = 1024;
    final double[][] data = new double[size * 2][];
    final UniformRandomProvider r = RngUtils.create(seed.getSeedAsLong());
    final GaussianSampler gs = GaussianSamplerUtils.createGaussianSampler(r, 0, 5);
    for (int x = 0, y = 0, y2 = size, i = 0; x < size; x++, y++, y2--) {
      data[i++] = new double[] {x + gs.sample(), y + gs.sample()};
      data[i++] = new double[] {x + gs.sample(), y2 + gs.sample()};
    }
    // Create 2 images
    final Rectangle bounds = new Rectangle(0, 0, size, size);
    IJImagePeakResults i1 = createImage(bounds);
    IJImagePeakResults i2 = createImage(bounds);
    final int[] indices = SimpleArrayUtils.natural(data.length);
    PermutationSampler.shuffle(r, indices);
    for (final int i : indices) {
      final IJImagePeakResults image = i1;
      i1 = i2;
      i2 = image;
      image.add((float) data[i][0], (float) data[i][1], 1);
    }
    i1.end();
    i2.end();
    final ImageProcessor ip1 = i1.getImagePlus().getProcessor();
    final ImageProcessor ip2 = i2.getImagePlus().getProcessor();
    // Test
    final FRC frc = new FRC();
    FloatProcessor[] fft1;
    FloatProcessor[] fft2;
    fft1 = frc.getComplexFFT(ip1);
    fft2 = frc.getComplexFFT(ip2);

    final float[] dataA1 = (float[]) fft1[0].getPixels();
    final float[] dataB1 = (float[]) fft1[1].getPixels();
    final float[] dataA2 = (float[]) fft2[0].getPixels();
    final float[] dataB2 = (float[]) fft2[1].getPixels();

    final float[] numeratorE = new float[dataA1.length];
    final float[] absFFT1E = new float[dataA1.length];
    final float[] absFFT2E = new float[dataA1.length];

    FRC.compute(numeratorE, absFFT1E, absFFT2E, dataA1, dataB1, dataA2, dataB2);

    Assertions.assertTrue(checkSymmetry(numeratorE, size), "numeratorE");
    Assertions.assertTrue(checkSymmetry(absFFT1E, size), "absFFT1E");
    Assertions.assertTrue(checkSymmetry(absFFT2E, size), "absFFT2E");

    final float[] numeratorA = new float[dataA1.length];
    final float[] absFFT1A = new float[dataA1.length];
    final float[] absFFT2A = new float[dataA1.length];
    FRC.computeMirrored(size, numeratorA, absFFT1A, absFFT2A, dataA1, dataB1, dataA2, dataB2);

    // for (int y=0, i=0; y<size; y++)
    // for (int x=0; x<size; x++, i++)
    // {
    // logger.fine(FunctionUtils.getSupplier("[%d,%d = %d] %f ?= %f", x, y, i, numeratorE[i],
    // numeratorA[i]);
    // }

    Assertions.assertArrayEquals(numeratorE, numeratorA, "numerator");
    Assertions.assertArrayEquals(absFFT1E, absFFT1A, "absFFT1");
    Assertions.assertArrayEquals(absFFT2E, absFFT2A, "absFFT2");

    FRC.computeMirroredFast(size, numeratorA, absFFT1A, absFFT2A, dataA1, dataB1, dataA2, dataB2);

    // Check this.
    for (int y = 1; y < size; y++) {
      for (int x = 1, i = y * size + 1; x < size; x++, i++) {
        Assertions.assertEquals(numeratorE[i], numeratorA[i], "numerator");
        Assertions.assertEquals(absFFT1E[i], absFFT1A[i], "absFFT1");
        Assertions.assertEquals(absFFT2E[i], absFFT2A[i], "absFFT2");
      }
    }
  }

  private static IJImagePeakResults createImage(Rectangle bounds) {
    final IJImagePeakResults i1 = new IJImagePeakResults("1", bounds, 1);
    i1.setDisplayImage(false);
    i1.begin();
    return i1;
  }

  /**
   * Check symmetry.
   *
   * @param data the data
   * @param size the size
   * @return true, if successful
   */
  private static boolean checkSymmetry(float[] data, int size) {
    // Symmetry is around the centre
    final int centre = size / 2;

    final float maxRelativeError = 1e-10f;
    final float maxAbsoluteError = 1e-16f;

    for (int y = centre, y2 = centre; y >= 0 && y2 < size; y--, y2++) {
      for (int x = centre, x2 = centre, i = size * y + x, j = size * y2 + x2; x >= 0 && x2 < size;
          x--, x2++, i--, j++) {
        if (data[i] != data[j] || !FloatEquality.almostEqualRelativeOrAbsolute(data[i], data[j],
            maxRelativeError, maxAbsoluteError)) {
          return false;
        }
      }
    }
    return true;
  }

  private abstract class MyTimingTask extends BaseTimingTask {
    public MyTimingTask(String name) {
      super(name);
    }

    @Override
    public int getSize() {
      return 1;
    }

    @Override
    public Object getData(int index) {
      return null;
    }
  }

  @SeededTest
  public void computeSineIsFaster() {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.HIGH));

    final int steps = 100000;
    final double delta = 2 * Math.PI / steps;
    final double[] angle = new double[steps + 1];
    final double[] cosAngle = new double[steps + 1];
    for (int i = 0; i <= steps; i++) {
      angle[i] = i * delta;
      cosAngle[i] = Math.cos(angle[i]);
    }

    final TimingService ts = new TimingService(100);
    ts.execute(new MyTimingTask("sin") {
      @Override
      public Object run(Object data) {
        double value = 0;
        for (int i = 0; i < angle.length; i++) {
          value += Math.sin(angle[i]);
        }
        return value;
      }
    });
    ts.execute(new MyTimingTask("FastMath.sin") {
      @Override
      public Object run(Object data) {
        double value = 0;
        for (int i = 0; i < angle.length; i++) {
          value += FastMath.sin(angle[i]);
        }
        return value;
      }
    });
    ts.execute(new MyTimingTask("getSine") {
      @Override
      public Object run(Object data) {
        double value = 0;
        for (int i = 0; i < angle.length; i++) {
          value += FRC.getSine(angle[i], cosAngle[i]);
        }
        return value;
      }
    });

    final int size = ts.getSize();
    ts.repeat(size);
    if (logger.isLoggable(Level.INFO)) {
      logger.info(ts.getReport(size));
    }

    Assertions.assertTrue(ts.get(-1).getMean() < ts.get(-2).getMean());
    Assertions.assertTrue(ts.get(-1).getMean() < ts.get(-3).getMean());
  }

  @SeededTest
  public void computeMirroredIsFaster(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    // Sample lines through an image to create a structure.
    final int N = 2048;
    final double[][] data = new double[N * 2][];
    final UniformRandomProvider r = RngUtils.create(seed.getSeedAsLong());
    final GaussianSampler gs = GaussianSamplerUtils.createGaussianSampler(r, 0, 5);
    for (int x = 0, y = 0, y2 = N, i = 0; x < N; x++, y++, y2--) {
      data[i++] = new double[] {x + gs.sample(), y + gs.sample()};
      data[i++] = new double[] {x + gs.sample(), y2 + gs.sample()};
    }
    // Create 2 images
    final Rectangle bounds = new Rectangle(0, 0, N, N);
    IJImagePeakResults i1 = createImage(bounds);
    IJImagePeakResults i2 = createImage(bounds);
    final int[] indices = SimpleArrayUtils.natural(data.length);
    PermutationSampler.shuffle(r, indices);
    for (final int i : indices) {
      final IJImagePeakResults image = i1;
      i1 = i2;
      i2 = image;
      image.add((float) data[i][0], (float) data[i][1], 1);
    }
    i1.end();
    i2.end();
    final ImageProcessor ip1 = i1.getImagePlus().getProcessor();
    final ImageProcessor ip2 = i2.getImagePlus().getProcessor();
    // Test
    final FRC frc = new FRC();
    FloatProcessor[] fft1;
    FloatProcessor[] fft2;
    fft1 = frc.getComplexFFT(ip1);
    fft2 = frc.getComplexFFT(ip2);

    final float[] dataA1 = (float[]) fft1[0].getPixels();
    final float[] dataB1 = (float[]) fft1[1].getPixels();
    final float[] dataA2 = (float[]) fft2[0].getPixels();
    final float[] dataB2 = (float[]) fft2[1].getPixels();

    final float[] numerator = new float[dataA1.length];
    final float[] absFFT1 = new float[dataA1.length];
    final float[] absFFT2 = new float[dataA1.length];

    final TimingService ts = new TimingService(10);
    ts.execute(new MyTimingTask("compute") {
      @Override
      public Object run(Object data) {
        FRC.compute(numerator, absFFT1, absFFT2, dataA1, dataB1, dataA2, dataB2);
        return null;
      }
    });
    ts.execute(new MyTimingTask("computeMirrored") {
      @Override
      public Object run(Object data) {
        FRC.computeMirrored(N, numerator, absFFT1, absFFT2, dataA1, dataB1, dataA2, dataB2);
        return null;
      }
    });
    ts.execute(new MyTimingTask("computeMirroredFast") {
      @Override
      public Object run(Object data) {
        FRC.computeMirroredFast(N, numerator, absFFT1, absFFT2, dataA1, dataB1, dataA2, dataB2);
        return null;
      }
    });

    final int size = ts.getSize();
    ts.repeat(size);
    if (logger.isLoggable(Level.INFO)) {
      logger.info(ts.getReport(size));
    }

    // The 'Fast' method may not always be faster so just log results
    final TimingResult slow = ts.get(-3);
    final TimingResult fast = ts.get(-2);
    final TimingResult fastest = ts.get(-1);
    logger.log(TestLogUtils.getTimingRecord(slow, fastest));
    logger.log(TestLogUtils.getTimingRecord(fast, fastest));
    // It should be faster than the non mirrored version
    Assertions.assertTrue(ts.get(-1).getMean() <= ts.get(-3).getMean());
  }
}
