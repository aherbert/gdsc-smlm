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

package uk.ac.sussex.gdsc.smlm.ij.utils;

import ij.ImageStack;
import java.util.Arrays;
import org.jtransforms.fft.DoubleFFT_3D;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.smlm.function.StandardValueProcedure;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import uk.ac.sussex.gdsc.smlm.function.gaussian.QuadraticAstigmatismZModel;

@SuppressWarnings({"javadoc"})
class DoubleDht3DTest {
  static int size = 16;
  static double centre = (size - 1) / 2.0;

  static final double gamma = 2;
  static final int zDepth = 5;
  private static QuadraticAstigmatismZModel zModel = new QuadraticAstigmatismZModel(gamma, zDepth);

  private static DoubleDht3D createData(double cx, double cy, double cz) {
    final Gaussian2DFunction f = GaussianFunctionFactory.create2D(1, size, size,
        GaussianFunctionFactory.FIT_ASTIGMATISM, zModel);
    final int length = size * size;
    final double[] data = new double[size * length];
    final double[] a = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
    a[Gaussian2DFunction.SIGNAL] = 1;
    a[Gaussian2DFunction.X_POSITION] = cx;
    a[Gaussian2DFunction.Y_POSITION] = cy;
    a[Gaussian2DFunction.X_SD] = 1;
    a[Gaussian2DFunction.Y_SD] = 1;
    final StandardValueProcedure p = new StandardValueProcedure();
    for (int z = 0; z < size; z++) {
      a[Gaussian2DFunction.Z_POSITION] = z - cz;
      p.getValues(f, a, data, z * length);
    }
    return new DoubleDht3D(size, size, size, data, false);
  }

  private static DoubleDht3D createData() {
    return createData(centre, centre, centre);
  }

  private static DoubleDht3D createOctants(int width, int height, int depth) {
    return new DoubleDht3D(FloatDht3DTest.createOctantsStack(width, height, depth));
  }

  @Test
  void canSwapOctants() {
    DoubleDht3D dht;

    // Simple test
    final double[] data = new double[] {2, 1, 3, 4, 6, 5, 7, 8};
    dht = new DoubleDht3D(2, 2, 2, data.clone(), false);
    dht.swapOctants();
    checkOctants(data, dht.getData());

    final int[] test = new int[] {2, 4, 6};
    for (final int w : test) {
      for (final int h : test) {
        for (final int d : test) {
          dht = createOctants(w, h, d);

          final double[] in = dht.getData().clone();

          // This just tests that the swap of the DHT and the stack matches
          final ImageStack stack = dht.getImageStack();
          // uk.ac.sussex.gdsc.core.ij.Utils.display("Test", stack);
          dht.swapOctants();
          FloatDht3D.swapOctants(stack);

          final double[] e = new DoubleDht3D(stack).getData();
          final double[] o = dht.getData();

          checkOctants(in, o);

          Assertions.assertArrayEquals(e, o);
        }
      }
    }
  }

  private static void checkOctants(double[] in, double[] out) {
    final int[] swap = new int[9];
    swap[1] = 7;
    swap[2] = 8;
    swap[3] = 5;
    swap[4] = 6;
    swap[5] = 3;
    swap[6] = 4;
    swap[7] = 1;
    swap[8] = 2;
    for (int i = 0; i < in.length; i++) {
      Assertions.assertEquals(in[i], swap[(int) out[i]]);
    }
  }

  @Test
  void canConvolveAndDeconvolve() {
    final DoubleDht3D dht = createData();
    final double[] pixels = dht.getData().clone();
    dht.transform();

    final DoubleDht3D copy = dht.copy();
    copy.initialiseFastMultiply();

    final DoubleDht3D convolved = dht.multiply(dht);
    final DoubleDht3D deconvolved = convolved.divide(dht);

    final DoubleDht3D convolved2 = dht.multiply(copy);
    final DoubleDht3D deconvolved2 = convolved.divide(copy);

    Assertions.assertArrayEquals(convolved.getData(), convolved2.getData());
    Assertions.assertArrayEquals(deconvolved.getData(), deconvolved2.getData());

    double[] exp = dht.getData();
    double[] obs = deconvolved.getData();
    for (int i = 0; i < exp.length; i++) {
      Assertions
          .assertTrue(DoubleEquality.almostEqualRelativeOrAbsolute(exp[i], obs[i], 1e-6, 1e-6));
    }

    deconvolved.inverseTransform();

    // Test after reverse transform
    exp = pixels;
    obs = deconvolved.getData();

    for (int i = 0; i < exp.length; i++) {
      Assertions
          .assertTrue(DoubleEquality.almostEqualRelativeOrAbsolute(exp[i], obs[i], 1e-8, 1e-8));
    }
  }

  @Test
  void canCorrelate() {
    final DoubleDht3D dht = createData();
    dht.transform();

    final DoubleDht3D copy = dht.copy();
    copy.initialiseFastMultiply();

    // Centre of power spectrum
    final int icentre = size / 2;

    for (int z = -1; z <= 1; z++) {
      for (int y = -1; y <= 1; y++) {
        for (int x = -1; x <= 1; x++) {
          final DoubleDht3D dht2 = createData(centre + x, centre + y, centre + z);
          dht2.transform();

          final DoubleDht3D correlation = dht2.conjugateMultiply(dht);
          final DoubleDht3D correlation2 = dht2.conjugateMultiply(copy);
          Assertions.assertArrayEquals(correlation.getData(), correlation2.getData());

          correlation.inverseTransform();
          correlation.swapOctants();

          final double[] pixels = correlation.getData();

          final int i = SimpleArrayUtils.findMaxIndex(pixels);
          final int[] xyz = correlation.getXyz(i);

          // This is how far dht has to move to align with dht2.
          // To align dht2 with dht would be the opposite sign.
          final int ox = xyz[0] - icentre;
          final int oy = xyz[1] - icentre;
          final int oz = xyz[2] - icentre;
          // logger.fine(FunctionUtils.getSupplier("Shift [%d,%d,%d], centre [%d,%d,%d]", x, y, z,
          // xyz[0], xyz[1], xyz[2]);
          Assertions.assertEquals(x, ox);
          Assertions.assertEquals(y, oy);
          Assertions.assertEquals(z, oz);
        }
      }
    }
  }

  @Test
  void canConvertToDft() {
    final DoubleDht3D dht = createData();
    final double[] input = dht.getData().clone();
    dht.transform();

    final DoubleImage3D[] result = dht.toDft(null, null);

    final double rel = 1e-14;
    final double abs = 1e-14;

    // Test reverse transform
    final DoubleDht3D dht2 = DoubleDht3D.fromDft(result[0], result[1], null);

    final double[] e = dht.getData();
    final double[] o = dht2.getData();
    for (int i = 0; i < e.length; i++) {
      Assertions.assertTrue(DoubleEquality.almostEqualRelativeOrAbsolute(e[i], o[i], rel, abs));
    }

    // Test verses full forward transform
    final DoubleFFT_3D fft = new DoubleFFT_3D(dht.ns, dht.nr, dht.nc);
    final double[] dft = Arrays.copyOf(input, 2 * e.length);
    fft.realForwardFull(dft);

    final double[] or = result[0].getData();
    final double[] oi = result[1].getData();
    for (int i = 0, j = 0; i < dft.length; i += 2, j++) {
      Assertions.assertTrue(DoubleEquality.almostEqualRelativeOrAbsolute(dft[i], or[j], rel, abs));
      Assertions
          .assertTrue(DoubleEquality.almostEqualRelativeOrAbsolute(dft[i + 1], oi[j], rel, abs));
    }
  }
}
