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

import uk.ac.sussex.gdsc.core.utils.ImageWindow;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.smlm.filters.FhtFilter.Operation;
import uk.ac.sussex.gdsc.test.api.TestAssertions;
import uk.ac.sussex.gdsc.test.api.TestHelper;
import uk.ac.sussex.gdsc.test.api.function.FloatFloatBiPredicate;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngUtils;
import uk.ac.sussex.gdsc.test.utils.TestCounter;
import uk.ac.sussex.gdsc.test.utils.functions.IndexSupplier;

import ij.plugin.filter.EDM;
import ij.process.ByteProcessor;
import ij.process.FHT;
import ij.process.FloatProcessor;

import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

@SuppressWarnings({"javadoc"})
public class FhtFilterTest {
  @SeededTest
  public void canCorrelate(RandomSeed seed) {
    canFilter(seed, Operation.CORRELATION);
  }

  @SeededTest
  public void canConvolve(RandomSeed seed) {
    canFilter(seed, Operation.CONVOLUTION);
  }

  @SeededTest
  public void canDeconvolve(RandomSeed seed) {
    canFilter(seed, Operation.DECONVOLUTION);
  }

  private static void canFilter(RandomSeed seed, Operation operation) {
    final int size = 16;
    final int ex = 5;
    final int ey = 7;
    final int ox = 1;
    final int oy = 2;
    final UniformRandomProvider r = RngUtils.create(seed.getSeed());
    final FloatProcessor fp1 = createProcessor(size, ex, ey, 4, 4, r);
    // This is offset from the centre
    final FloatProcessor fp2 = createProcessor(size, size / 2 + ox, size / 2 + oy, 4, 4, r);

    final float[] input1 = ((float[]) fp1.getPixels()).clone();
    final float[] input2 = ((float[]) fp2.getPixels()).clone();

    final FHT fht1 = new FHT(fp1);
    fht1.transform();
    final FHT fht2 = new FHT(fp2);
    fht2.transform();

    FHT fhtE;
    switch (operation) {
      case CONVOLUTION:
        fhtE = fht1.multiply(fht2);
        break;
      case CORRELATION:
        fhtE = fht1.conjugateMultiply(fht2);
        break;
      case DECONVOLUTION:
        fhtE = fht1.divide(fht2);
        break;
      default:
        throw new RuntimeException();
    }
    fhtE.inverseTransform();
    fhtE.swapQuadrants();

    final float[] e = (float[]) fhtE.getPixels();
    if (operation == Operation.CORRELATION) {
      // Test the max correlation position
      final int max = SimpleArrayUtils.findMaxIndex(e);
      final int x = max % 16;
      final int y = max / 16;

      Assertions.assertEquals(ex, x + ox);
      Assertions.assertEquals(ey, y + oy);
    }

    // Test verses a spatial domain filter in the middle of the image
    if (operation != Operation.DECONVOLUTION) {
      double sum = 0;
      float[] i2 = input2;
      if (operation == Operation.CONVOLUTION) {
        i2 = i2.clone();
        KernelFilter.rotate180(i2);
      }
      for (int i = 0; i < input1.length; i++) {
        sum += input1[i] * i2[i];
      }
      // double exp = e[size / 2 * size + size / 2];
      // logger.fine(() -> String.format("Sum = %f vs [%d] %f", sum, size / 2 * size + size / 2,
      // exp);
      Assertions.assertEquals(sum, sum, 1e-3);
    }

    // Test the FHT filter
    final FhtFilter ff = new FhtFilter(input2, size, size);
    ff.setOperation(operation);
    ff.filter(input1, size, size);

    // There may be differences due to the use of the JTransforms library
    final double error = (operation == Operation.DECONVOLUTION) ? 5e-2 : 1e-4;
    final FloatFloatBiPredicate predicate = TestHelper.floatsAreClose(error, 0);

    // This tests everything and can fail easily depending on the random generator
    // due to edge artifacts.
    // TestAssertions.assertArrayTest(e, input1, TestHelper.almostEqualFloats(error, 0));

    // This tests the centre to ignore edge differences
    final int min = size / 4;
    final int max = size - min;
    int repeats = 0;
    for (int y = min; y < max; y++) {
      for (int x = min; x < max; x++) {
        repeats++;
      }
    }

    // Use a fail counter for a 'soft' test that detects major problems
    final int failureLimit = TestCounter.computeFailureLimit(repeats, 0.1);
    final TestCounter failCounter = new TestCounter(failureLimit);

    final IndexSupplier msg = new IndexSupplier(2);
    for (int y = min; y < max; y++) {
      msg.set(1, y);
      for (int x = min; x < max; x++) {
        final int xx = x;
        final int i = y * size + x;
        failCounter.run(() -> {
          TestAssertions.assertTest(e[i], input1[i], predicate, msg.set(0, xx));
        });
      }
    }
  }

  private static FloatProcessor createProcessor(int size, int x, int y, int width, int height,
      UniformRandomProvider rng) {
    final ByteProcessor bp = new ByteProcessor(size, size);
    bp.setColor(255);
    bp.fillOval(x, y, width, height);
    final EDM e = new EDM();
    final FloatProcessor fp = e.makeFloatEDM(bp, 0, true);
    if (rng != null) {
      final float[] d = (float[]) fp.getPixels();
      for (int i = 0; i < d.length; i++) {
        d[i] += rng.nextFloat() * 0.01;
      }
    }
    return fp;
  }

  @Test
  public void canWindow() {
    final int size = 16;
    final float[] in = SimpleArrayUtils.newFloatArray(size * size, 1);
    final FhtFilter f = new FhtFilter(new float[1], 1, 1);
    for (int i = 1; i < 5; i++) {
      final double[] wx = ImageWindow.tukeyEdge(size, i);
      final float[] e = ImageWindow.applyWindowSeparable(in, size, size, wx, wx);
      final float[] o = in.clone();
      f.applyBorder(o, size, size, i);
      Assertions.assertArrayEquals(e, o);
    }
  }
}
