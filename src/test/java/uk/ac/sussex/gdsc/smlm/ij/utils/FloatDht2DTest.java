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

package uk.ac.sussex.gdsc.smlm.ij.utils;

import ij.process.FloatProcessor;
import java.util.Arrays;
import org.jtransforms.fft.FloatFFT_2D;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import uk.ac.sussex.gdsc.core.ij.process.Fht;
import uk.ac.sussex.gdsc.core.utils.FloatEquality;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.smlm.function.StandardFloatValueProcedure;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.GaussianFunctionFactory;

@SuppressWarnings({"javadoc"})
public class FloatDht2DTest {
  static int size = 16;
  static double centre = (size - 1) / 2.0;

  private static FloatDht2D createData(double cx, double cy) {
    final Gaussian2DFunction f = GaussianFunctionFactory.create2D(1, size, size,
        GaussianFunctionFactory.FIT_ERF_FREE_CIRCLE, null);
    final double[] a = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
    a[Gaussian2DFunction.SIGNAL] = 1;
    a[Gaussian2DFunction.X_POSITION] = cx;
    a[Gaussian2DFunction.Y_POSITION] = cy;
    a[Gaussian2DFunction.X_SD] = 1.2;
    a[Gaussian2DFunction.Y_SD] = 1.1;
    final StandardFloatValueProcedure p = new StandardFloatValueProcedure();
    p.getValues(f, a);
    return new FloatDht2D(size, size, p.values, false);
  }

  private static FloatDht2D createData() {
    return createData(centre, centre);
  }

  private static FloatDht2D createQuadrants(int width, int height) {
    return new FloatDht2D(createQuadrantsProcessor(width, height));
  }

  static FloatProcessor createQuadrantsProcessor(int width, int height) {
    final int w_2 = width / 2;
    final int stepH2 = height / 2;
    final FloatProcessor fp = new FloatProcessor(width, height);
    fill(fp, w_2, 0, w_2, stepH2, 1);
    fill(fp, 0, 0, w_2, stepH2, 2);
    fill(fp, 0, stepH2, w_2, stepH2, 3);
    fill(fp, w_2, stepH2, w_2, stepH2, 4);
    return fp;
  }

  static void fill(FloatProcessor fp, int x, int y, int width, int height, double value) {
    fp.setRoi(x, y, width, height);
    fp.setValue(value);
    fp.fill();
  }

  @Test
  public void canSwapQuadrants() {
    FloatDht2D dht;

    // Simple test
    final float[] data = new float[] {2, 1, 3, 4};
    dht = new FloatDht2D(2, 2, data.clone(), false);
    dht.swapQuadrants();
    checkQuadrants(data, dht.getData());

    final int[] test = new int[] {2, 4, 6};
    for (final int w : test) {
      for (final int h : test) {
        dht = createQuadrants(w, h);

        final float[] in = dht.getData().clone();

        dht.swapQuadrants();

        final float[] o = dht.getData();

        checkQuadrants(in, o);
      }
    }
  }

  private static void checkQuadrants(float[] in, float[] out) {
    final int[] swap = new int[9];
    swap[1] = 3;
    swap[2] = 4;
    swap[3] = 1;
    swap[4] = 2;
    for (int i = 0; i < in.length; i++) {
      Assertions.assertEquals(in[i], swap[(int) out[i]]);
    }
  }

  @Test
  public void canConvolveAndDeconvolve() {
    final FloatDht2D dht = createData();

    final float[] pixels = dht.getData().clone();
    dht.transform();

    final FloatDht2D copy = dht.copy();
    copy.initialiseFastMultiply();

    final FloatDht2D convolved = dht.multiply(dht);
    final FloatDht2D deconvolved = convolved.divide(dht);

    final FloatDht2D convolved2 = dht.multiply(copy);
    final FloatDht2D deconvolved2 = convolved.divide(copy);

    Assertions.assertArrayEquals(convolved.getData(), convolved2.getData());
    Assertions.assertArrayEquals(deconvolved.getData(), deconvolved2.getData());

    float[] exp = dht.getData();
    float[] obs = deconvolved.getData();
    for (int i = 0; i < exp.length; i++) {
      Assertions
          .assertTrue(FloatEquality.almostEqualRelativeOrAbsolute(exp[i], obs[i], 1e-6f, 1e-6f));
    }

    deconvolved.inverseTransform();

    // Test after reverse transform
    exp = pixels;
    obs = deconvolved.getData();

    for (int i = 0; i < exp.length; i++) {
      Assertions
          .assertTrue(FloatEquality.almostEqualRelativeOrAbsolute(exp[i], obs[i], 1e-7f, 1e-7f));
    }
  }

  @Test
  public void canCorrelate() {
    final FloatDht2D dht = createData();
    dht.transform();

    final FloatDht2D copy = dht.copy();
    copy.initialiseFastMultiply();

    // Centre of power spectrum
    final int icentre = size / 2;

    for (int y = -1; y <= 1; y++) {
      for (int x = -1; x <= 1; x++) {
        final FloatDht2D dht2 = createData(centre + x, centre + y);
        dht2.transform();

        final FloatDht2D correlation = dht2.conjugateMultiply(dht);
        final FloatDht2D correlation2 = dht2.conjugateMultiply(copy);
        Assertions.assertArrayEquals(correlation.getData(), correlation2.getData());

        correlation.inverseTransform();
        correlation.swapQuadrants();

        final float[] pixels = correlation.getData();

        final int i = SimpleArrayUtils.findMaxIndex(pixels);
        final int[] xy = correlation.getXy(i);

        // This is how far dht has to move to align with dht2.
        // To align dht2 with dht would be the opposite sign.
        final int ox = xy[0] - icentre;
        final int oy = xy[1] - icentre;
        // logger.fine(FunctionUtils.getSupplier("Shift [%d,%d], centre [%d,%d]", x, y, xy[0],
        // xy[1]);
        Assertions.assertEquals(x, ox);
        Assertions.assertEquals(y, oy);
      }
    }
  }

  @Test
  public void canConvertToDft() {
    final FloatDht2D dht = createData();
    final float[] input = dht.getData().clone();
    dht.transform();

    final FloatImage2D[] result = dht.toDft(null, null);

    final float rel = 1e-6f;
    final float abs = 1e-6f;

    // Test reverse transform
    final FloatDht2D dht2 = FloatDht2D.fromDft(result[0], result[1], null);

    final float[] e = dht.getData();
    final float[] o = dht2.getData();
    for (int i = 0; i < e.length; i++) {
      Assertions.assertTrue(FloatEquality.almostEqualRelativeOrAbsolute(e[i], o[i], rel, abs));
    }

    // Test verses full forward transform
    final FloatFFT_2D fft = new FloatFFT_2D(dht.nr, dht.nc);
    final float[] dft = Arrays.copyOf(input, 2 * e.length);
    fft.realForwardFull(dft);

    final float[] or = result[0].getData();
    final float[] oi = result[1].getData();
    for (int i = 0, j = 0; i < dft.length; i += 2, j++) {
      Assertions.assertTrue(FloatEquality.almostEqualRelativeOrAbsolute(dft[i], or[j], rel, abs));
      Assertions
          .assertTrue(FloatEquality.almostEqualRelativeOrAbsolute(dft[i + 1], oi[j], rel, abs));
    }
  }

  @Test
  public void canMatchFht2() {
    final FloatDht2D dht = createData();
    final FloatDht2D dht2 = createData(centre + 1, centre + 1);

    final float[] pixels = dht.getData().clone();
    final float[] pixels2 = dht2.getData().clone();

    final Fht fht = new Fht(pixels, size, false);
    final Fht fht2 = new Fht(pixels2, size, false);

    dht.transform();
    dht2.transform();
    fht.transform();
    fht2.transform();

    Assertions.assertArrayEquals(fht.getData(), dht.getData(), 1e-6f);
    Assertions.assertArrayEquals(fht2.getData(), dht2.getData(), 1e-6f);

    check("multiply", dht.multiply(dht2), fht.multiply(fht2), 1e-6f, 1e-6f);
    check("conjugateMultiply", dht.conjugateMultiply(dht2), fht.conjugateMultiply(fht2), 1e-6f,
        1e-6f);
    check("divide", dht.divide(dht2), fht.divide(fht2), 1e-3f, 1e-3f);
  }

  private static void check(String operation, FloatDht2D dht, Fht fht, float rel, float abs) {
    final float[] e = fht.getData();
    final float[] o = dht.getData();
    for (int i = 0; i < e.length; i++) {
      if (!FloatEquality.almostEqualRelativeOrAbsolute(e[i], o[i], rel, abs)) {
        Assertions.fail(String.format("%s [%d] %g vs %g = %g", operation, i, e[i], o[i],
            FloatEquality.relativeError(e[i], o[i])));
      }
    }
  }
}
