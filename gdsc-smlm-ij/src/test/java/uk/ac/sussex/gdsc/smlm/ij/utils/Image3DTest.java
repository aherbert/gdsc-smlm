/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Package
 *
 * Software for single molecule localisation microscopy (SMLM) in ImageJ
 * %%
 * Copyright (C) 2011 - 2025 Alex Herbert
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
import ij.process.ImageProcessor;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

@SuppressWarnings({"javadoc"})
public abstract class Image3DTest {
  protected abstract Image3D createData(int width, int height, int depth);

  protected abstract Image3D createEmptyData(int width, int height, int depth);

  @Test
  void canCrop() {
    canCrop(3, 4, 5, 6, 7, 8);
    canCrop(3, 4, 5, 1, 1, 1);
    canCrop(0, 0, 0, 1, 2, 3);
  }

  private void canCrop(int x, int y, int z, int width, int height, int depth) {
    final Image3D image = createData(x + width + 1, y + height + 1, z + depth + 1);
    final ImageStack stack = image.getImageStack();

    final Image3D croppedData = image.crop(x, y, z, width, height, depth);
    Assertions.assertEquals(croppedData.getWidth(), width);
    Assertions.assertEquals(croppedData.getHeight(), height);
    Assertions.assertEquals(croppedData.getSize(), depth);

    final Image3D croppedData2 = FloatImage3D.crop(stack, x, y, z, width, height, depth, null);
    assertEquals(croppedData, croppedData2);

    final ImageStack croppedStack = image.cropToStack(x, y, z, width, height, depth);
    Assertions.assertEquals(croppedStack.getWidth(), width);
    Assertions.assertEquals(croppedStack.getHeight(), height);
    Assertions.assertEquals(croppedStack.getSize(), depth);

    assertEquals(croppedData, new FloatImage3D(croppedStack));

    final ImageStack croppedStack2 = Image3D.cropToStack(stack, x, y, z, width, height, depth);
    assertEquals(croppedData, new FloatImage3D(croppedStack2));

    // Test it is the correct region
    final ImageStack originalStack = image.getImageStack();
    for (int zz = 0; zz < depth; zz++) {
      // Crop from the original data
      ImageProcessor fp = originalStack.getProcessor(1 + zz + z);
      fp.setRoi(x, y, width, height);
      fp = fp.crop();
      final float[] e = (float[]) fp.getPixels();

      // Compare to the cropped stack
      final float[] o = (float[]) croppedStack.getPixels(1 + zz);
      Assertions.assertArrayEquals(e, o);

      // Compare to the cropped data
      for (int i = 0, j = zz * e.length; i < e.length; i++, j++) {
        Assertions.assertEquals(e[i], croppedData.get(j));
      }
    }
  }

  @Test
  void canInsert() {
    canInsert(3, 4, 5, 6, 7, 8);
    canInsert(3, 4, 5, 1, 1, 1);
    canInsert(0, 0, 0, 1, 2, 3);
  }

  private void canInsert(int x, int y, int z, int width, int height, int depth) {
    // This test assumes that copy and crop work!
    for (final int pad : new int[] {0, 1}) {
      final Image3D image = createData(x + width + pad, y + height + pad, z + depth + pad);
      final Image3D image2 = image.copy();
      final Image3D image3 = image.copy();

      final Image3D blank = createEmptyData(width, height, depth);

      image.insert(x, y, z, blank);

      Image3D croppedData = image.crop(x, y, z, width, height, depth);

      assertEquals(croppedData, blank);

      final ImageStack blankStack = blank.getImageStack();
      image2.insert(x, y, z, blankStack);

      croppedData = image2.crop(x, y, z, width, height, depth);

      assertEquals(croppedData, blank);

      for (int s = 0; s < blankStack.getSize(); s++) {
        image3.insert(x, y, z + s, blankStack.getProcessor(1 + s));
      }

      croppedData = image3.crop(x, y, z, width, height, depth);

      assertEquals(croppedData, blank);
    }
  }

  @Test
  void canFindMin() {
    final Image3D image = createData(2, 2, 2);
    Assertions.assertEquals(0, image.findMinIndex(0, 0, 0, 2, 2, 2));
    Assertions.assertEquals(1, image.findMinIndex(1, 0, 0, 2, 2, 2));
    Assertions.assertEquals(2, image.findMinIndex(0, 1, 0, 2, 2, 2));
    Assertions.assertEquals(3, image.findMinIndex(1, 1, 0, 2, 2, 2));
    Assertions.assertEquals(4, image.findMinIndex(0, 0, 1, 2, 2, 2));
    Assertions.assertEquals(5, image.findMinIndex(1, 0, 1, 2, 2, 2));
    Assertions.assertEquals(6, image.findMinIndex(0, 1, 1, 2, 2, 2));
    Assertions.assertEquals(7, image.findMinIndex(1, 1, 1, 2, 2, 2));

    // Larger slices
    canFindMin(3, 4, 5, 6, 7, 8);
    canFindMin(3, 4, 5, 1, 1, 1);
    canFindMin(0, 0, 0, 1, 2, 3);
  }

  private void canFindMin(int x, int y, int z, int width, int height, int depth) {
    // This test assumes that crop works!
    for (final int pad : new int[] {0, 1}) {
      final Image3D image = createData(x + width + pad, y + height + pad, z + depth + pad);

      final Image3D croppedData = image.crop(x, y, z, width, height, depth);
      final int i = findMinIndex(croppedData);
      final int[] xyz = croppedData.getXyz(i);

      final int j = image.findMinIndex(x, y, z, width, height, depth);
      final int[] xyz2 = image.getXyz(j);

      Assertions.assertEquals(xyz[0] + x, xyz2[0]);
      Assertions.assertEquals(xyz[1] + y, xyz2[1]);
      Assertions.assertEquals(xyz[2] + z, xyz2[2]);
    }
  }

  @Test
  void canFindMax() {
    final Image3D image = createData(2, 2, 2);
    Assertions.assertEquals(7, image.findMaxIndex(0, 0, 0, 2, 2, 2));
    Assertions.assertEquals(6, image.findMaxIndex(0, 0, 0, 1, 2, 2));
    Assertions.assertEquals(5, image.findMaxIndex(0, 0, 0, 2, 1, 2));
    Assertions.assertEquals(4, image.findMaxIndex(0, 0, 0, 1, 1, 2));
    Assertions.assertEquals(3, image.findMaxIndex(0, 0, 0, 2, 2, 1));
    Assertions.assertEquals(2, image.findMaxIndex(0, 0, 0, 1, 2, 1));
    Assertions.assertEquals(1, image.findMaxIndex(0, 0, 0, 2, 1, 1));
    Assertions.assertEquals(0, image.findMaxIndex(0, 0, 0, 1, 1, 1));

    // Larger slices
    canFindMax(3, 4, 5, 6, 7, 8);
    canFindMax(3, 4, 5, 1, 1, 1);
    canFindMax(0, 0, 0, 1, 2, 3);
  }

  private void canFindMax(int x, int y, int z, int width, int height, int depth) {
    // This test assumes that crop works!
    for (final int pad : new int[] {0, 1}) {
      final Image3D image = createData(x + width + pad, y + height + pad, z + depth + pad);

      final Image3D croppedData = image.crop(x, y, z, width, height, depth);
      final int i = findMaxIndex(croppedData);
      final int[] xyz = croppedData.getXyz(i);

      final int j = image.findMaxIndex(x, y, z, width, height, depth);
      final int[] xyz2 = image.getXyz(j);

      Assertions.assertEquals(xyz[0] + x, xyz2[0]);
      Assertions.assertEquals(xyz[1] + y, xyz2[1]);
      Assertions.assertEquals(xyz[2] + z, xyz2[2]);
    }
  }

  @Test
  void canComputeSum() {
    // Bounds checks
    final Image3D image = createData(2, 2, 2);
    Assertions.assertEquals(36, image.computeSum(0, 0, 0, 2, 2, 2));
    Assertions.assertEquals(0, image.computeSum(0, 0, 0, 0, 0, 0));
    Assertions.assertEquals(1, image.computeSum(0, 0, 0, 1, 1, 1));
    Assertions.assertEquals(2, image.computeSum(1, 0, 0, 1, 1, 1));
    Assertions.assertEquals(0, image.computeSum(-10, 0, 0, 1, 1, 1));
    Assertions.assertEquals(0, image.computeSum(10, 0, 0, 1, 1, 1));
    Assertions.assertEquals(0, image.computeSum(0, 10, 0, 1, 1, 1));
    Assertions.assertEquals(0, image.computeSum(0, -10, 0, 1, 1, 1));
    Assertions.assertEquals(0, image.computeSum(0, 0, 10, 1, 1, 1));
    Assertions.assertEquals(0, image.computeSum(0, 0, -10, 1, 1, 1));

    // Larger slices
    canComputeSum(3, 4, 5, 6, 7, 8);
    canComputeSum(3, 4, 5, 1, 1, 1);
    canComputeSum(0, 0, 0, 1, 2, 3);
  }

  private void canComputeSum(int x, int y, int z, int width, int height, int depth) {
    // This test assumes that crop works!
    for (final int pad : new int[] {0, 1}) {
      final Image3D image = createData(x + width + pad, y + height + pad, z + depth + pad);

      final Image3D croppedData = image.crop(x, y, z, width, height, depth);
      final double e = sum(croppedData);
      final double o = image.computeSum(x, y, z, width, height, depth);

      Assertions.assertEquals(o, e);
    }
  }

  @Test
  void canComputeRollingSumTable() {
    final int width = 2;
    final int height = 3;
    final int depth = 4;
    final Image3D image = createData(width, height, depth);
    final double[] table = image.computeRollingSumTable(null);
    for (int dd = 1, i = 0; dd <= depth; dd++) {
      for (int hh = 1; hh <= height; hh++) {
        for (int ww = 1; ww <= width; ww++) {
          final double e = image.computeSum(0, 0, 0, ww, hh, dd);
          final double o = table[i++];
          Assertions.assertEquals(e, o, 1e-3);
        }
      }
    }
  }

  @Test
  void canComputeSumUsingTable() {
    // Bounds checks
    final Image3D image = createData(2, 2, 2);
    final double[] table = image.computeRollingSumTable(null);
    // testComputeSum(36, image,table, 0, 0, 0, 2, 2, 2);
    testComputeSum(36, image, table, 0, 0, 0, 5, 7, 9);
    testComputeSum(0, image, table, 0, 0, 0, 0, 0, 0);
    testComputeSum(1, image, table, 0, 0, 0, 1, 1, 1);
    testComputeSum(2, image, table, 1, 0, 0, 1, 1, 1);
    testComputeSum(0, image, table, -10, 0, 0, 1, 1, 1);
    testComputeSum(0, image, table, 10, 0, 0, 1, 1, 1);
    testComputeSum(0, image, table, 0, 10, 0, 1, 1, 1);
    testComputeSum(0, image, table, 0, -10, 0, 1, 1, 1);
    testComputeSum(0, image, table, 0, 0, 10, 1, 1, 1);
    testComputeSum(0, image, table, 0, 0, -10, 1, 1, 1);

    // Larger slices
    canComputeSumUsingTable(3, 4, 5, 6, 7, 8);
    canComputeSumUsingTable(3, 4, 5, 1, 1, 1);
    canComputeSumUsingTable(0, 0, 0, 1, 2, 3);
  }

  private void canComputeSumUsingTable(int x, int y, int z, int width, int height, int depth) {
    // This test assumes that crop works!
    for (final int pad : new int[] {0, 1}) {
      final Image3D image = createData(x + width + pad, y + height + pad, z + depth + pad);
      final double[] table = image.computeRollingSumTable(null);

      final Image3D croppedData = image.crop(x, y, z, width, height, depth);
      final double e = sum(croppedData);

      testComputeSum(e, image, table, x, y, z, width, height, depth);
    }
  }

  private static void testComputeSum(double exp, Image3D image, double[] table, int x, int y, int z,
      int width, int height, int depth) {
    final double o1 = image.computeSum(table, x, y, z, width, height, depth);
    final double o2 = image.computeSumFast(table, x, y, z, width, height, depth);

    // This may be different due to floating point error
    // but we are adding integers so it should be OK
    Assertions.assertEquals(exp, o1);
    Assertions.assertEquals(exp, o2);
  }

  @Test
  void canFill() {
    canFill(3, 4, 5, 6, 7, 8);
    canFill(3, 4, 5, 1, 1, 1);
    canFill(0, 0, 0, 1, 2, 3);
  }

  private void canFill(int x, int y, int z, int width, int height, int depth) {
    // This test assumes that copy, crop and insert work!
    for (final int pad : new int[] {0, 1}) {
      final Image3D image = createEmptyData(x + width + pad, y + height + pad, z + depth + pad);
      image.fill(1);
      image.fill(x, y, z, width, height, depth, 2);

      final Image3D image2 = createEmptyData(x + width + pad, y + height + pad, z + depth + pad);
      image2.fill(1);
      final Image3D blank = createEmptyData(width, height, depth);
      blank.fill(2);
      image2.insert(x, y, z, blank);

      final Image3D croppedData = image.crop(x, y, z, width, height, depth);

      assertEquals(croppedData, blank);

      assertEquals(image, image2);
    }
  }

  @Test
  void canFillOutside() {
    canFillOutside(3, 4, 5, 6, 7, 8);
    canFillOutside(3, 4, 5, 1, 1, 1);
    canFillOutside(0, 0, 0, 1, 2, 3);
  }

  private void canFillOutside(int x, int y, int z, int width, int height, int depth) {
    // This test assumes that copy, crop and insert work!
    for (final int pad : new int[] {0, 1}) {
      final Image3D image = createEmptyData(x + width + pad, y + height + pad, z + depth + pad);
      image.fill(1);
      image.fillOutside(x, y, z, width, height, depth, 2);

      final Image3D image2 = createEmptyData(x + width + pad, y + height + pad, z + depth + pad);
      image2.fill(2);
      final Image3D blank = createEmptyData(width, height, depth);
      blank.fill(1);
      image2.insert(x, y, z, blank);

      final Image3D croppedData = image.crop(x, y, z, width, height, depth);

      assertEquals(croppedData, blank);

      assertEquals(image, image2);
    }
  }

  private static void assertEquals(Image3D i1, Image3D i2) {
    for (int i = i1.getDataLength(); i-- > 0;) {
      if (i1.get(i) != i2.get(i)) {
        Assertions.assertEquals(i1.get(i), i2.get(i), "Not equal @ " + i);
      }
    }
  }

  private static int findMinIndex(Image3D image) {
    int min = 0;
    for (int i = image.getDataLength(); i-- > 0;) {
      if (image.get(i) < image.get(min)) {
        min = i;
      }
    }
    return min;
  }

  private static int findMaxIndex(Image3D image) {
    int max = 0;
    for (int i = image.getDataLength(); i-- > 0;) {
      if (image.get(i) > image.get(max)) {
        max = i;
      }
    }
    return max;
  }

  private static double sum(Image3D image) {
    double sum = 0;
    for (int i = image.getDataLength(); i-- > 0;) {
      sum += image.get(i);
    }
    return sum;
  }
}
