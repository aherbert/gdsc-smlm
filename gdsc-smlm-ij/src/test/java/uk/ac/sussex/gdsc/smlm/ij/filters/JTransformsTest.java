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

package uk.ac.sussex.gdsc.smlm.ij.filters;

import ij.plugin.filter.EDM;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import org.apache.commons.rng.UniformRandomProvider;
import org.jtransforms.dht.FloatDHT_2D;
import org.jtransforms.fft.FloatFFT_2D;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import uk.ac.sussex.gdsc.core.ij.process.Fht;
import uk.ac.sussex.gdsc.smlm.ij.filters.FhtFilter.Operation;

@SuppressWarnings({"javadoc"})
class JTransformsTest {

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
  void canCorrelateUsingFft() {
    canComputeUsingFft(false);
  }

  @Test
  void canConvolveUsingFft() {
    canComputeUsingFft(true);
  }

  private static void canComputeUsingFft(boolean convolution) {
    final int size = 16;
    final int ex = 5;
    final int ey = 7;
    final int ox = 1;
    final int oy = 2;
    final FloatProcessor fp1 = createProcessor(size, ex, ey, 4, 4, null);
    final FloatProcessor fp2 = createProcessor(size, size / 2 + ox, size / 2 + oy, 4, 4, null);

    final float[] input1 = (float[]) fp1.getPixels();
    final float[] input2 = (float[]) fp2.getPixels();

    final FhtFilter ff = new FhtFilter(input2.clone(), size, size);
    ff.setOperation((convolution) ? Operation.CONVOLUTION : Operation.CORRELATION);
    final float[] e = input1.clone();
    ff.filter(e, size, size);

    // Do the same with JTransforms
    final float[] data1 = new float[input1.length * 2];
    final FloatFFT_2D fft = new FloatFFT_2D(size, size);
    System.arraycopy(input1, 0, data1, 0, input1.length);
    final float[] data2 = new float[data1.length];
    System.arraycopy(input2, 0, data2, 0, input2.length);

    fft.realForwardFull(data1);
    fft.realForwardFull(data2);

    // Multiply
    // https://en.wikipedia.org/wiki/Complex_number#Multiplication_and_division
    for (int i = 0; i < data2.length; i += 2) {
      final float a = data1[i];
      final float b = data1[i + 1];
      final float c = data2[i];
      // Get the conjugate for correlation
      final float d = (convolution) ? data2[i + 1] : -data2[i + 1];
      data1[i] = a * c - b * d;
      data1[i + 1] = b * c + a * d;
    }

    fft.complexInverse(data1, true);

    final float[] o = new float[e.length];
    for (int i = 0, j = 0; i < o.length; i++, j += 2) {
      o[i] = data1[j];
    }
    Fht.swapQuadrants(new FloatProcessor(size, size, o));

    Assertions.assertArrayEquals(e, o, 1e-3f);
  }

  @Test
  void canComputeFhtUsingJTransforms() {
    // Note: no need to test the correlation as the transformed data
    // is the same format as FHT so we just test that.

    final int size = 16;
    final int ex = 5;
    final int ey = 7;
    final int ox = 1;
    final int oy = 2;
    final FloatProcessor fp1 = createProcessor(size, ex, ey, 4, 4, null);
    final FloatProcessor fp2 = createProcessor(size, size / 2 + ox, size / 2 + oy, 4, 4, null);

    final float[] input1 = (float[]) fp1.getPixels();
    final float[] input2 = (float[]) fp2.getPixels();

    final Fht fht1 = new Fht(input1.clone(), size, false);
    final Fht fht2 = new Fht(input2.clone(), size, false);

    fht1.transform();
    fht2.transform();

    // Do the same with JTransforms
    final FloatDHT_2D dht = new FloatDHT_2D(size, size);

    dht.forward(input1);
    dht.forward(input2);

    Assertions.assertArrayEquals(fht1.getData(), input1, 1e-5f);
    Assertions.assertArrayEquals(fht2.getData(), input2, 1e-5f);
  }
}
