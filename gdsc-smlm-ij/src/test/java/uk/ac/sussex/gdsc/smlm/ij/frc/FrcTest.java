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

package uk.ac.sussex.gdsc.smlm.ij.frc;

import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import java.awt.Rectangle;
import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.sampling.PermutationSampler;
import org.apache.commons.rng.sampling.distribution.SharedStateContinuousSampler;
import org.junit.jupiter.api.Assertions;
import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.FloatEquality;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.rng.SamplerUtils;
import uk.ac.sussex.gdsc.smlm.ij.results.ImageJImagePeakResults;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngFactory;
import uk.ac.sussex.gdsc.test.utils.RandomSeed;

@SuppressWarnings({"javadoc"})
class FrcTest {

  @SeededTest
  void canComputeSine() {
    final int steps = 1000;
    final double delta = 2 * Math.PI / steps;
    for (int i = 0; i <= steps; i++) {
      final double a = i * delta;
      final double cosA = Math.cos(a);
      final double e = Math.sin(a);
      final double o = Frc.getSine(a, cosA);
      // logger.fine(FunctionUtils.getSupplier("%f %f ?= %f", a, e, o);
      Assertions.assertTrue(DoubleEquality.almostEqualRelativeOrAbsolute(o, e, 1e-6, 1e-10));
    }
  }

  @SeededTest
  void canComputeMirrored(RandomSeed seed) {
    // Sample lines through an image to create a structure.
    final int size = 1024;
    final double[][] data = new double[size * 2][];
    final UniformRandomProvider r = RngFactory.create(seed.get());
    final SharedStateContinuousSampler gs = SamplerUtils.createGaussianSampler(r, 0, 5);
    for (int x = 0, y = 0, y2 = size, i = 0; x < size; x++, y++, y2--) {
      data[i++] = new double[] {x + gs.sample(), y + gs.sample()};
      data[i++] = new double[] {x + gs.sample(), y2 + gs.sample()};
    }
    // Create 2 images
    final Rectangle bounds = new Rectangle(0, 0, size, size);
    ImageJImagePeakResults i1 = createImage(bounds);
    ImageJImagePeakResults i2 = createImage(bounds);
    final int[] indices = SimpleArrayUtils.natural(data.length);
    PermutationSampler.shuffle(r, indices);
    for (final int i : indices) {
      final ImageJImagePeakResults image = i1;
      i1 = i2;
      i2 = image;
      image.add((float) data[i][0], (float) data[i][1], 1);
    }
    i1.end();
    i2.end();
    final ImageProcessor ip1 = i1.getImagePlus().getProcessor();
    final ImageProcessor ip2 = i2.getImagePlus().getProcessor();
    // Test
    final Frc frc = new Frc();
    FloatProcessor[] fft1;
    FloatProcessor[] fft2;
    fft1 = frc.getComplexFft(ip1);
    fft2 = frc.getComplexFft(ip2);

    final float[] dataA1 = (float[]) fft1[0].getPixels();
    final float[] dataB1 = (float[]) fft1[1].getPixels();
    final float[] dataA2 = (float[]) fft2[0].getPixels();
    final float[] dataB2 = (float[]) fft2[1].getPixels();

    final float[] numeratorE = new float[dataA1.length];
    final float[] absFft1E = new float[dataA1.length];
    final float[] absFft2E = new float[dataA1.length];

    Frc.compute(numeratorE, absFft1E, absFft2E, dataA1, dataB1, dataA2, dataB2);

    Assertions.assertTrue(checkSymmetry(numeratorE, size), "numeratorE");
    Assertions.assertTrue(checkSymmetry(absFft1E, size), "absFft1E");
    Assertions.assertTrue(checkSymmetry(absFft2E, size), "absFft2E");

    final float[] numeratorA = new float[dataA1.length];
    final float[] absFft1A = new float[dataA1.length];
    final float[] absFft2A = new float[dataA1.length];
    Frc.computeMirrored(size, numeratorA, absFft1A, absFft2A, dataA1, dataB1, dataA2, dataB2);

    // for (int y=0, i=0; y<size; y++)
    // for (int x=0; x<size; x++, i++)
    // {
    // logger.fine(FunctionUtils.getSupplier("[%d,%d = %d] %f ?= %f", x, y, i, numeratorE[i],
    // numeratorA[i]);
    // }

    Assertions.assertArrayEquals(numeratorE, numeratorA, "numerator");
    Assertions.assertArrayEquals(absFft1E, absFft1A, "absFft1");
    Assertions.assertArrayEquals(absFft2E, absFft2A, "absFft2");

    Frc.computeMirroredFast(size, numeratorA, absFft1A, absFft2A, dataA1, dataB1, dataA2, dataB2);

    // Check this.
    for (int y = 1; y < size; y++) {
      for (int x = 1, i = y * size + 1; x < size; x++, i++) {
        Assertions.assertEquals(numeratorE[i], numeratorA[i], "numerator");
        Assertions.assertEquals(absFft1E[i], absFft1A[i], "absFft1");
        Assertions.assertEquals(absFft2E[i], absFft2A[i], "absFft2");
      }
    }
  }

  private static ImageJImagePeakResults createImage(Rectangle bounds) {
    final ImageJImagePeakResults i1 = new ImageJImagePeakResults("1", bounds, 1);
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
}
