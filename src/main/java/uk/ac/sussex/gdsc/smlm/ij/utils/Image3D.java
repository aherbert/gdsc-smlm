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

import uk.ac.sussex.gdsc.core.utils.MathUtils;

import ij.ImageStack;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

/**
 * Store a 3D image in a single array.
 */
public abstract class Image3D {
  /**
   * The largest array size for which a regular 1D Java array is used to store the data (2^30).
   */
  public static final int MAX_SIZE_OF_32_BIT_ARRAY = 1073741824;

  /** The number of slices (max z). */
  public final int ns;
  /** The number of rows (max y). */
  public final int nr;
  /** The number of columns (max x). */
  public final int nc;

  /** The number of rows multiplied by the number of columns. */
  public final int nrByNc;

  /**
   * Instantiates a new 3D image.
   *
   * @param nc the number of columns
   * @param nr the number of rows
   * @param ns the number of slices
   * @throws IllegalArgumentException If the combined dimensions is too large for an array
   */
  public Image3D(int nc, int nr, int ns) {
    createData(checkSize(nc, nr, ns, true));
    this.nc = nc;
    this.nr = nr;
    this.ns = ns;
    nrByNc = nr * nc;
  }

  /**
   * Instantiates a new 3D image.
   *
   * @param stack the stack
   * @throws IllegalArgumentException If the combined dimensions is too large for an array
   */
  public Image3D(ImageStack stack) {
    nc = stack.getWidth();
    nr = stack.getHeight();
    ns = stack.getSize();
    createData(checkSize(nc, nr, ns, true));
    nrByNc = nr * nc;
    if (stack.getBitDepth() == 32) {
      for (int s = 0; s < ns; s++) {
        copyFrom((float[]) stack.getPixels(s + 1), 0, nrByNc, s * nrByNc);
      }
    } else {
      for (int s = 1, i = 0; s <= ns; s++) {
        final ImageProcessor ip = stack.getProcessor(s);
        for (int j = 0; i < nrByNc; j++) {
          setf(i++, ip.getf(j));
        }
      }
    }
  }

  /**
   * Instantiates a new 3D image.
   *
   * <p>It is assumed that the sub-class will correctly create the data storage.
   *
   * @param nc the number of columns
   * @param nr the number of rows
   * @param ns the number of slices
   * @param nrByNc the number of rows multiplied by the number of columns
   */
  protected Image3D(int nc, int nr, int ns, int nrByNc) {
    // No checks as this is used internally
    this.nc = nc;
    this.nr = nr;
    this.ns = ns;
    this.nrByNc = nrByNc;
  }

  /**
   * Copy constructor.
   *
   * <p>It is assumed that the sub-class will correctly create the data storage.
   *
   * @param source the source
   */
  protected Image3D(Image3D source) {
    this.nc = source.nc;
    this.nr = source.nr;
    this.ns = source.ns;
    this.nrByNc = source.nrByNc;
  }

  /**
   * Check the size can fit in a 1D array.
   *
   * @param nc the number of columns
   * @param nr the number of rows
   * @param ns the number of slices
   * @param raiseException Set to true to raise an exception if too large
   * @return the size, or -1 if too large
   * @throws IllegalArgumentException if too large (optional)
   */
  public static int checkSize(int nc, int nr, int ns, boolean raiseException) {
    if (nc < 0 || nr < 0 || ns < 0) {
      if (raiseException) {
        throw new IllegalArgumentException("Negative dimensions");
      }
      return -1;
    }
    final long size = (long) ns * nr * nc;
    if (size > MAX_SIZE_OF_32_BIT_ARRAY) {
      if (raiseException) {
        throw new IllegalArgumentException("3D data too large");
      }
      return -1;
    }
    return (int) size;
  }

  /**
   * Creates the array to store the data.
   *
   * @param size the size of the array
   */
  protected abstract void createData(int size);

  /**
   * Copy the data from the given index to the buffer.
   *
   * <p>Utility method to handle conversion with ImageJ ImageProcessor objects.
   *
   * @param index the index
   * @param buffer the buffer
   * @param bufferIndex the buffer index
   * @param size the size to copy
   */
  protected abstract void copyTo(int index, float[] buffer, int bufferIndex, int size);

  /**
   * Copy the data from the given buffer to the given index.
   *
   * <p>Utility method to handle conversion with ImageJ ImageProcessor objects.
   *
   * @param buffer the buffer
   * @param bufferIndex the buffer index
   * @param size the size to copy
   * @param index the index
   */
  protected abstract void copyFrom(float[] buffer, int bufferIndex, int size, int index);

  /**
   * Gets the value at the given index.
   *
   * @param index the index
   * @return the value
   */
  public abstract double get(int index);

  /**
   * Sets the value at the given index.
   *
   * @param index the index
   * @param value the value
   */
  public abstract void set(int index, double value);

  /**
   * Gets the value at the given index.
   *
   * @param index the index
   * @return the value
   */
  public abstract float getf(int index);

  /**
   * Sets the value at the given index.
   *
   * @param index the index
   * @param value the value
   */
  public abstract void setf(int index, float value);

  /**
   * Return a copy of the 3D image.
   *
   * @return the copy
   */
  public abstract Image3D copy();

  /**
   * Gets the width (the number of columns).
   *
   * @return the width
   */
  public int getWidth() {
    return nc;
  }

  /**
   * Gets the height (the number of rows).
   *
   * @return the height
   */
  public int getHeight() {
    return nr;
  }

  /**
   * Gets the size (the number of slices).
   *
   * @return the size
   */
  public int getSize() {
    return ns;
  }

  /**
   * Gets the data length.
   *
   * @return the data length
   */
  public int getDataLength() {
    return ns * nrByNc;
  }

  /**
   * Convert to an image stack.
   *
   * @return the image stack
   */
  public ImageStack getImageStack() {
    final ImageStack stack = new ImageStack(nc, nr);
    for (int s = 0; s < ns; s++) {
      final float[] pixels = new float[nrByNc];
      copyTo(s * nrByNc, pixels, 0, nrByNc);
      stack.addSlice(null, pixels);
    }
    return stack;
  }

  /**
   * Gets the xyz components of the index.
   *
   * @param index the index
   * @return the xyz components
   * @throws IllegalArgumentException if the index is not within the data
   */
  public int[] getXyz(int index) {
    if (index < 0 || index >= getDataLength()) {
      throw new IllegalArgumentException(
          "Index in not in the correct range: 0 <= index < " + getDataLength());
    }
    final int[] xyz = new int[3];
    xyz[2] = index / nrByNc;
    final int xy = index % nrByNc;
    xyz[1] = xy / nc;
    xyz[0] = xy % nc;
    return xyz;
  }

  /**
   * Gets the xyz components of the index.
   *
   * @param index the index
   * @param xyz the xyz components (must be an array of at least length 3)
   * @throws IllegalArgumentException if the index is not within the data
   */
  public void getXyz(int index, int[] xyz) {
    if (index < 0 || index >= getDataLength()) {
      throw new IllegalArgumentException(
          "Index in not in the correct range: 0 <= index < " + getDataLength());
    }
    xyz[2] = index / nrByNc;
    final int xy = index % nrByNc;
    xyz[1] = xy / nc;
    xyz[0] = xy % nc;
  }

  /**
   * Gets the index using the xyz components.
   *
   * @param x the x
   * @param y the y
   * @param z the z
   * @return the index
   * @throws IllegalArgumentException if the index is not within the data
   */
  public int getIndex(int x, int y, int z) {
    if (x < 0 || x >= nc || y < 0 || y >= nr || z < 0 || z >= ns) {
      throw new IllegalArgumentException("Index in not inside the image");
    }
    return index(x, y, z);
  }

  /**
   * Gets the index using the xyz components.
   *
   * @param x the x
   * @param y the y
   * @param z the z
   * @return the index
   */
  protected int index(int x, int y, int z) {
    return z * nrByNc + y * nc + x;
  }

  /**
   * Crop a sub-region of the data. The target dimensions must be positive.
   *
   * @param x the x index
   * @param y the y index
   * @param z the z index
   * @param width the width
   * @param height the height
   * @param depth the depth
   * @return the cropped data
   * @throws IllegalArgumentException if the region is not within the data
   */
  public abstract Image3D crop(int x, int y, int z, int width, int height, int depth);

  /**
   * Crop a sub-region of the data into the given image. The target dimensions must be positive.
   *
   * @param x the x index
   * @param y the y index
   * @param z the z index
   * @param image the image
   * @return the cropped data
   * @throws IllegalArgumentException if the region is not within the data
   */
  public Image3D crop(int x, int y, int z, Image3D image) {
    // Check the region range
    final int width = image.getWidth();
    final int height = image.getHeight();
    final int depth = image.getSize();
    if (x < 0 || width < 1 || (long) x + width > nc || y < 0 || height < 1 || (long) y + height > nr
        || z < 0 || depth < 1 || (long) z + depth > ns) {
      throw new IllegalArgumentException("Region not within the data");
    }
    for (int s = 0, i = 0; s < depth; s++, z++) {
      int base = z * nrByNc + y * nc + x;
      for (int r = 0; r < height; r++) {
        for (int c = 0; c < width; c++) {
          image.set(i++, get(base + c));
        }
        base += nc;
      }
    }
    return image;
  }

  /**
   * Crop a sub-region of the data. The target dimensions must be positive.
   *
   * @param x the x index
   * @param y the y index
   * @param z the z index
   * @param width the width
   * @param height the height
   * @param depth the depth
   * @return the cropped data
   * @throws IllegalArgumentException if the region is not within the data
   */
  public ImageStack cropToStack(int x, int y, int z, int width, int height, int depth) {
    // Check the region range
    if (x < 0 || width < 1 || (long) x + width > nc || y < 0 || height < 1 || (long) y + height > nr
        || z < 0 || depth < 1 || (long) z + depth > ns) {
      throw new IllegalArgumentException("Region not within the data");
    }
    final int size = width * height;
    final ImageStack stack = new ImageStack(width, height, depth);
    for (int s = 0; s < depth; s++, z++) {
      int base = z * nrByNc + y * nc + x;
      final float[] region = new float[size];
      for (int r = 0, i = 0; r < height; r++) {
        copyTo(base, region, i, width);
        base += nc;
        i += width;
      }
      stack.setPixels(region, 1 + s);
    }
    return stack;
  }

  /**
   * Crop a sub-region of the data. The target dimensions must be positive.
   *
   * @param stack the stack
   * @param x the x index
   * @param y the y index
   * @param z the z index
   * @param width the width
   * @param height the height
   * @param depth the depth
   * @return the cropped data
   * @throws IllegalArgumentException if the region is not within the data
   */
  public static ImageStack cropToStack(ImageStack stack, int x, int y, int z, int width, int height,
      int depth) {
    final int nc = stack.getWidth();
    final int nr = stack.getHeight();
    final int ns = stack.getSize();

    // Check the region range
    if (x < 0 || width < 1 || (long) x + width > nc || y < 0 || height < 1 || (long) y + height > nr
        || z < 0 || depth < 1 || (long) z + depth > ns) {
      throw new IllegalArgumentException("Region not within the data");
    }
    final ImageStack stack2 = new ImageStack(width, height, depth);
    for (int s = 0; s < depth; s++, z++) {
      final ImageProcessor ip = stack.getProcessor(1 + z);
      ip.setRoi(x, y, width, height);
      stack2.setPixels(ip.crop().getPixels(), 1 + s);
    }
    return stack2;
  }

  /**
   * Insert a sub-region.
   *
   * @param x the x position
   * @param y the y position
   * @param z the z position
   * @param image the image
   * @throws IllegalArgumentException if the region is not within the data
   */
  public void insert(int x, int y, int z, Image3D image) {
    // Check the region range
    final int width = image.getWidth();
    final int height = image.getHeight();
    final int depth = image.getSize();
    if (width < 1 || height < 1 || depth < 1) {
      return;
    }
    if (x < 0 || (long) x + width > nc || y < 0 || (long) y + height > nr || z < 0
        || (long) z + depth > ns) {
      throw new IllegalArgumentException("Region not within the data");
    }
    for (int s = 0, i = 0; s < depth; s++, z++) {
      int base = z * nrByNc + y * nc + x;
      for (int r = 0; r < height; r++) {
        for (int c = 0; c < width; c++) {
          set(base + c, image.get(i++));
        }
        base += nc;
      }
    }
  }

  /**
   * Insert a sub-region.
   *
   * @param x the x position
   * @param y the y position
   * @param z the z position
   * @param stack the image stack
   * @throws IllegalArgumentException if the region is not within the data
   */
  public void insert(int x, int y, int z, ImageStack stack) {
    // Check the region range
    final int width = stack.getWidth();
    final int height = stack.getHeight();
    final int depth = stack.getSize();
    if (width < 1 || height < 1 || depth < 1) {
      return;
    }
    if (x < 0 || (long) x + width > nc || y < 0 || (long) y + height > nr || z < 0
        || (long) z + depth > ns) {
      throw new IllegalArgumentException("Region not within the data");
    }
    final boolean isFloat = stack.getBitDepth() == 32;
    final FloatProcessor fp = (isFloat) ? new FloatProcessor(width, height) : null;
    for (int s = 0; s < depth; s++, z++) {
      int base = z * nrByNc + y * nc + x;
      final float[] region = (float[]) ((isFloat) ? stack.getPixels(1 + s)
          : stack.getProcessor(1 + s).toFloat(0, fp).getPixels());
      for (int r = 0, i = 0; r < height; r++) {
        copyFrom(region, i, width, base);
        base += nc;
        i += width;
      }
    }
  }

  /**
   * Insert a sub-region.
   *
   * @param x the x position
   * @param y the y position
   * @param z the z position
   * @param image the image
   * @throws IllegalArgumentException if the region is not within the data
   */
  public void insert(int x, int y, int z, ImageProcessor image) {
    // Check the region range
    final int width = image.getWidth();
    final int height = image.getHeight();
    if (width < 1 || height < 1) {
      return;
    }
    if (x < 0 || (long) x + width > nc || y < 0 || (long) y + height > nr || z < 0 || z >= ns) {
      throw new IllegalArgumentException("Region not within the data");
    }
    final boolean isFloat = image.getBitDepth() == 32;
    int base = z * nrByNc + y * nc + x;
    final float[] region =
        (float[]) ((isFloat) ? image.getPixels() : image.toFloat(0, null).getPixels());
    for (int r = 0, i = 0; r < height; r++) {
      copyFrom(region, i, width, base);
      base += nc;
      i += width;
    }
  }

  /**
   * Compute 3D intersect with this object.
   *
   * <p>If any of width,height,depth are negative then the corresponding x,y,z is updated and the
   * width,height,depth is inverted. The maximum bounds of the given dimensions are then computed by
   * adding the width,height,depth to the x,y,z. The bounds are then clipped to the image dimensions
   * and the intersect returned.
   *
   * @param x the x index
   * @param y the y index
   * @param z the z index
   * @param width the width
   * @param height the height
   * @param depth the depth
   * @return [x,y,z,width,height,depth]
   */
  public int[] computeIntersect(int x, int y, int z, int width, int height, int depth) {
    if (width < 0) {
      width = -width;
      x = subtract(x, width);
    }
    if (height < 0) {
      height = -height;
      y = subtract(y, height);
    }
    if (depth < 0) {
      depth = -depth;
      z = subtract(z, depth);
    }
    // Compute 3D intersect with this object
    final int x2 = clip(nc, x, width);
    final int y2 = clip(nr, y, height);
    final int z2 = clip(ns, z, depth);
    x = clip(nc, x);
    y = clip(nr, y);
    z = clip(ns, z);
    return new int[] {x, y, z, x2 - x, y2 - y, z2 - z};
  }

  /**
   * Compute 3D intersect with this object or throw an exception if the intersect has no volume.
   *
   * <p>If any of width,height,depth are negative then the corresponding x,y,z is updated and the
   * width,height,depth is inverted. The maximum bounds of the given dimensions are then computed by
   * adding the width,height,depth to the x,y,z. The bounds are then clipped to the image dimensions
   * and the intersect returned.
   *
   * @param x the x index
   * @param y the y index
   * @param z the z index
   * @param width the width
   * @param height the height
   * @param depth the depth
   * @return [x,y,z,width,height,depth]
   * @throws IllegalArgumentException if the intersect has no volume
   */
  public int[] computeIntersectOrThrow(int x, int y, int z, int width, int height, int depth) {
    if (width < 0) {
      width = -width;
      x = subtract(x, width);
    }
    if (height < 0) {
      height = -height;
      y = subtract(y, height);
    }
    if (depth < 0) {
      depth = -depth;
      z = subtract(z, depth);
    }
    // Compute 3D intersect with this object
    final int x2 = clip(nc, x, width);
    final int y2 = clip(nr, y, height);
    final int z2 = clip(ns, z, depth);
    x = clip(nc, x);
    y = clip(nr, y);
    z = clip(ns, z);
    width = checkIntersectSize(x2 - x);
    height = checkIntersectSize(y2 - y);
    depth = checkIntersectSize(z2 - z);
    return new int[] {x, y, z, width, height, depth};
  }

  /**
   * Check size is above zero or throw an exception.
   *
   * @param size the size
   * @return the size
   * @throws IllegalArgumentException If the size if zero
   */
  private static int checkIntersectSize(int size) {
    if (size == 0) {
      throw new IllegalArgumentException("No intersect");
    }
    return size;
  }

  /**
   * Subtract the value avoiding underflow.
   *
   * @param value the value
   * @param subtraction the subtraction
   * @return the value
   */
  private static int subtract(int value, int subtraction) {
    // Avoid underflow
    final long v = (long) value - subtraction;
    return (v < Integer.MIN_VALUE) ? Integer.MIN_VALUE : (int) v;
  }

  /**
   * Return value clipped to within the given bounds (0 - upper).
   *
   * @param upper the upper limit
   * @param value the value
   * @param addition the addition
   * @return the clipped value
   */
  private static int clip(int upper, int value, int addition) {
    // Avoid overflow
    final long v = (long) value + addition;
    if (v < 0) {
      return 0;
    }
    if (v > upper) {
      return upper;
    }
    return (int) v;
  }

  /**
   * Return value clipped to within the given bounds (0 - upper).
   *
   * @param upper the upper limit
   * @param value the value
   * @return the clipped value
   */
  private static int clip(int upper, int value) {
    return MathUtils.clip(0, upper, value);
  }

  /**
   * Find the index of the minimum value in the region.
   *
   * @param x the x index
   * @param y the y index
   * @param z the z index
   * @param width the width
   * @param height the height
   * @param depth the depth
   * @return the min index
   * @throws IllegalArgumentException if there is no intersect
   * @throws IllegalStateException if the region is just NaN values
   */
  public int findMinIndex(int x, int y, int z, int width, int height, int depth) {
    final int[] intersect = computeIntersectOrThrow(x, y, z, width, height, depth);
    x = intersect[0];
    y = intersect[1];
    z = intersect[2];
    width = intersect[3];
    height = intersect[4];
    depth = intersect[5];
    int index = findValueIndex(x, y, z, width, height, depth);
    double min = get(index);
    for (int s = 0; s < depth; s++, z++) {
      int base = z * nrByNc + y * nc + x;
      for (int r = 0; r < height; r++) {
        for (int j = 0; j < width; j++) {
          if (get(base + j) < min) {
            index = base + j;
            min = get(index);
          }
        }
        base += nc;
      }
    }
    return index;
  }

  /**
   * Find the index of the minimum value in the region.
   *
   * @param x the x index
   * @param y the y index
   * @param z the z index
   * @param width the width
   * @param height the height
   * @param depth the depth
   * @return the min index
   * @throws IllegalArgumentException if there is no intersect
   * @throws IllegalStateException if the region is just NaN values
   */
  public int findMaxIndex(int x, int y, int z, int width, int height, int depth) {
    final int[] intersect = computeIntersectOrThrow(x, y, z, width, height, depth);
    x = intersect[0];
    y = intersect[1];
    z = intersect[2];
    width = intersect[3];
    height = intersect[4];
    depth = intersect[5];
    int index = findValueIndex(x, y, z, width, height, depth);
    double max = get(index);
    for (int s = 0; s < depth; s++, z++) {
      int base = z * nrByNc + y * nc + x;
      for (int r = 0; r < height; r++) {
        for (int j = 0; j < width; j++) {
          if (get(base + j) > max) {
            index = base + j;
            max = get(index);
          }
        }
        base += nc;
      }
    }
    return index;
  }

  /**
   * Find the index of the first non-NaN value in the region.
   *
   * @param x the x index
   * @param y the y index
   * @param z the z index
   * @param width the width
   * @param height the height
   * @param depth the depth
   * @return the index
   * @throws IllegalStateException if the region is just NaN values
   */
  private int findValueIndex(int x, int y, int z, int width, int height, int depth) {
    // Quick check without loops
    if (!Double.isNaN(get(z * nrByNc + y * nc + x))) {
      return z * nrByNc + y * nc + x;
    }

    for (int s = 0; s < depth; s++, z++) {
      int base = z * nrByNc + y * nc + x;
      for (int r = 0; r < height; r++) {
        for (int j = 0; j < width; j++) {
          if (!Double.isNaN(get(base + j))) {
            return base + j;
          }
        }
        base += nc;
      }
    }
    throw new IllegalStateException("Region is NaN");
  }

  /**
   * Compute the rolling sum table for use in
   * {@link #computeSum(double[], int, int, int, int, int, int)}.
   *
   * <p>This is a table of the sum of the volume from 0,0,0 to x,y,z inclusive.
   *
   * @param table the table (will be reused if the correct size)
   * @return the rolling sum table
   */
  public double[] computeRollingSumTable(double[] table) {
    if (table == null || table.length != getDataLength()) {
      table = new double[getDataLength()];
    }

    // First build a table for each XY slice
    for (int s = 0; s < ns; s++) {
      double sum = 0;
      int index = s * nrByNc;
      // Initialise first row sum
      // sum = rolling sum of (0 - column)
      for (int c = 0; c < nc; c++, index++) {
        sum += get(index);
        table[index] = sum;
      }
      // Remaining rows
      // sum = rolling sum of (0 - column) + sum of same position above
      for (int r = 1, ii = index - nc; r < nr; r++) {
        sum = 0;
        for (int c = 0; c < nc; c++, index++) {
          sum += get(index);
          // Add the sum from the previous row
          table[index] = sum + table[ii++];
        }
      }
    }

    // Now sum across slices
    // sum = rolling sum of (0,0 to column,row) + sum of same position above
    // => rolling sum of (0,0,0 to column,row,slice)
    for (int s = 1; s < ns; s++) {
      int i1 = s * nrByNc;
      int i2 = i1 - nrByNc;
      for (int j = 0; j < nrByNc; j++) {
        table[i1++] += table[i2++];
      }
    }

    return table;
  }

  /**
   * Compute the sum of the region.
   *
   * @param x the x index
   * @param y the y index
   * @param z the z index
   * @param width the width
   * @param height the height
   * @param depth the depth
   * @return the sum
   */
  public double computeSum(int x, int y, int z, int width, int height, int depth) {
    final int[] intersect = computeIntersect(x, y, z, width, height, depth);
    width = intersect[3];
    height = intersect[4];
    depth = intersect[5];
    // Recheck bounds
    if (width == 0 || height == 0 || depth == 0) {
      return 0;
    }
    x = intersect[0];
    y = intersect[1];
    z = intersect[2];
    double sum = 0;
    for (int s = 0; s < depth; s++, z++) {
      int base = z * nrByNc + y * nc + x;
      for (int r = 0; r < height; r++) {
        for (int j = 0; j < width; j++) {
          sum += get(base + j);
        }
        base += nc;
      }
    }
    return sum;
  }

  /**
   * Compute the sum of the region using the precomputed rolling sum table.
   *
   * @param table the rolling sum table
   * @param x the x index
   * @param y the y index
   * @param z the z index
   * @param width the width
   * @param height the height
   * @param depth the depth
   * @return the sum
   */
  public double computeSum(double[] table, int x, int y, int z, int width, int height, int depth) {
    final int[] intersect = computeIntersect(x, y, z, width, height, depth);
    width = intersect[3];
    height = intersect[4];
    depth = intersect[5];
    // Recheck bounds
    if (width == 0 || height == 0 || depth == 0) {
      return 0;
      // x = intersect[0];
      // y = intersect[1];
      // z = intersect[2];
    }

    // Compute sum from rolling sum using:
    // sum(x,y,z,width,height,depth) =
    // + s(x+width-1,y+height-1,z+depth-1)
    // - s(x-1,y+height-1,z+depth-1)
    // - s(x+width-1,y-1,z+depth-1)
    // + s(x-1,y-1,z+depth-1)
    // /* Stack above must be subtracted so reverse sign*/
    // - s(x+width-1,y+height-1,z-1)
    // + s(x-1,y+height-1,z-1)
    // + s(x+width-1,y-1,z-1)
    // - s(x-1,y-1,z-1)
    // Note:
    // s(i,j,k) = 0 when either i,j,k < 0
    // i = imax when i>imax
    // j = jmax when j>jmax
    // k = kmax when k>kmax

    final int x_1 = intersect[0] - 1;
    final int y_1 = intersect[1] - 1;
    final int z_1 = intersect[2] - 1;
    // The intersect has already checked the bounds
    // int x_w_1 = Math.min(x_1 + width, nc);
    // int y_h_1 = Math.min(y_1 + height, nr);
    // int z_d_1 = Math.min(z_1 + depth, ns);
    final int x_w_1 = x_1 + width;
    final int y_h_1 = y_1 + height;
    final int z_d_1 = z_1 + depth;

    // double sum = table[index(x_w_1, y_h_1, z_d_1)];
    // if (y_1 >= 0)
    // {
    // sum -= table[index(x_w_1, y_1, z_d_1)];
    // if (x_1 >= 0)
    // sum = sum + table[index(x_1, y_1, z_d_1)] - table[index(x_1, y_h_1, z_d_1)];
    // }
    // else if (x_1 >= 0)
    // {
    // sum -= table[index(x_1, y_h_1, z_d_1)];
    // }
    // if (z_1 >= 0)
    // {
    // sum -= table[index(x_w_1, y_h_1, z_1)];
    // if (y_1 >= 0)
    // {
    // sum += table[index(x_w_1, y_1, z_1)];
    // if (x_1 >= 0)
    // sum = sum - table[index(x_1, y_1, z_1)] + table[index(x_1, y_h_1, z_1)];
    // }
    // else if (x_1 >= 0)
    // {
    // sum += table[index(x_1, y_h_1, z_1)];
    // }
    // }
    // return sum;

    // This has been ordered to use the smallest sums first (i.e. closer to x,y,z than
    // x+width,y+height,z+depth)
    final int xw_yh_zd = index(x_w_1, y_h_1, z_d_1);
    if (z_1 >= 0) {
      final int xw_yh_z = xw_yh_zd - depth * nrByNc;
      double sum = 0;
      if (y_1 >= 0) {
        final int h_ = height * nc;
        if (x_1 >= 0) {
          sum = table[xw_yh_zd - width - h_] - table[xw_yh_z - width - h_] - table[xw_yh_zd - width]
              + table[xw_yh_z - width];
        }
        sum = sum + table[xw_yh_z - h_] - table[xw_yh_zd - h_];
      } else if (x_1 >= 0) {
        sum = table[xw_yh_z - width] - table[xw_yh_zd - width];
      }
      return sum + table[xw_yh_zd] - table[xw_yh_z];
    }
    double sum = 0;
    if (y_1 >= 0) {
      final int h_ = height * nc;
      if (x_1 >= 0) {
        sum = table[xw_yh_zd - width - h_] - table[xw_yh_zd - width];
      }
      sum -= table[xw_yh_zd - h_];
    } else if (x_1 >= 0) {
      sum = -table[xw_yh_zd - width];
    }
    return sum + table[xw_yh_zd];
  }

  /**
   * Compute the sum of the region using the precomputed rolling sum table. Assumes
   * x+width,y+height,z+depth will not overflow!
   *
   * @param table the rolling sum table
   * @param x the x index
   * @param y the y index
   * @param z the z index
   * @param width the width
   * @param height the height
   * @param depth the depth
   * @return the sum
   */
  public double computeSumFast(double[] table, int x, int y, int z, int width, int height,
      int depth) {
    if (width <= 0 || height <= 0 || depth <= 0 || x >= nc || y >= nr || z >= ns) {
      return 0;
    }

    // Compute sum from rolling sum using:
    // sum(x,y,z,width,height,depth) =
    // + s(x+width-1,y+height-1,z+depth-1)
    // - s(x-1,y+height-1,z+depth-1)
    // - s(x+width-1,y-1,z+depth-1)
    // + s(x-1,y-1,z+depth-1)
    // /* Stack above must be subtracted so reverse sign*/
    // - s(x+width-1,y+height-1,z-1)
    // + s(x-1,y+height-1,z-1)
    // + s(x+width-1,y-1,z-1)
    // - s(x-1,y-1,z-1)
    // Note:
    // s(i,j,k) = 0 when either i,j,k < 0
    // i = imax when i>imax
    // j = jmax when j>jmax
    // k = kmax when k>kmax

    // Compute bounds assuming width,height,depth is small and positive.
    int x1;
    int y1;
    int z1;
    int xw1;
    int yh1;
    int zd1;
    if (x < 0) {
      x1 = 0;
      xw1 = MathUtils.clip(0, nc, x + width);
    } else {
      x1 = x;
      xw1 = Math.min(nc, x + width);
    }
    width = xw1 - x1;
    if (width == 0) {
      return 0;
    }
    if (y < 0) {
      y1 = 0;
      yh1 = MathUtils.clip(0, nr, y + height);
    } else {
      y1 = y;
      yh1 = Math.min(nr, y + height);
    }
    height = yh1 - y1;
    if (height == 0) {
      return 0;
    }
    if (z < 0) {
      z1 = 0;
      zd1 = MathUtils.clip(0, ns, z + depth);
    } else {
      z1 = z;
      zd1 = Math.min(ns, z + depth);
    }
    depth = zd1 - z1;
    if (depth == 0) {
      return 0;
    }

    // Adjust for the -1
    x1--;
    y1--;
    z1--;
    xw1--;
    yh1--;
    zd1--;

    // This has been ordered to use the smallest sums first (i.e. closer to x,y,z than
    // x+width,y+height,z+depth)
    final int xw_yh_zd = index(xw1, yh1, zd1);
    if (z1 >= 0) {
      final int xw_yh_z = xw_yh_zd - depth * nrByNc;
      double sum = 0;
      if (y1 >= 0) {
        final int h_ = height * nc;
        if (x1 >= 0) {
          sum = table[xw_yh_zd - width - h_] - table[xw_yh_z - width - h_] - table[xw_yh_zd - width]
              + table[xw_yh_z - width];
        }
        sum = sum + table[xw_yh_z - h_] - table[xw_yh_zd - h_];
      } else if (x1 >= 0) {
        sum = table[xw_yh_z - width] - table[xw_yh_zd - width];
      }
      return sum + table[xw_yh_zd] - table[xw_yh_z];
    }
    double sum = 0;
    if (y1 >= 0) {
      final int h_ = height * nc;
      if (x1 >= 0) {
        sum = table[xw_yh_zd - width - h_] - table[xw_yh_zd - width];
      }
      sum -= table[xw_yh_zd - h_];
    } else if (x1 >= 0) {
      sum = -table[xw_yh_zd - width];
    }
    return sum + table[xw_yh_zd];
  }

  /**
   * Fill the image.
   *
   * @param value the value
   */
  public void fill(double value) {
    fill(0, getDataLength(), value);
  }

  /**
   * Fill the region.
   *
   * @param x the x index
   * @param y the y index
   * @param z the z index
   * @param width the width
   * @param height the height
   * @param depth the depth
   * @param value the value
   */
  public void fill(int x, int y, int z, int width, int height, int depth, double value) {
    final int[] intersect = computeIntersect(x, y, z, width, height, depth);
    width = intersect[3];
    height = intersect[4];
    depth = intersect[5];
    // Recheck bounds
    if (width == 0 || height == 0 || depth == 0) {
      return;
    }
    x = intersect[0];
    y = intersect[1];
    z = intersect[2];
    for (int s = 0; s < depth; s++, z++) {
      int base = z * nrByNc + y * nc + x;
      for (int r = 0; r < height; r++) {
        fill(base, width, value);
        base += nc;
      }
    }
  }

  /**
   * Fill with the given value from the given index.
   *
   * @param index the index
   * @param size the size to fill
   * @param value the value
   */
  protected abstract void fill(int index, int size, double value);

  /**
   * Fill outside the region. If the region is not within the image then the entire image is filled.
   *
   * @param x the x index
   * @param y the y index
   * @param z the z index
   * @param width the width
   * @param height the height
   * @param depth the depth
   * @param value the value
   */
  public void fillOutside(int x, int y, int z, int width, int height, int depth, double value) {
    final int[] intersect = computeIntersect(x, y, z, width, height, depth);
    width = intersect[3];
    height = intersect[4];
    depth = intersect[5];
    // Recheck bounds
    if (width == 0 || height == 0 || depth == 0) {
      fill(value);
      return;
    }
    x = intersect[0];
    y = intersect[1];
    z = intersect[2];

    // Before
    if (z > 0) {
      fill(0, z * nrByNc, value);
    }
    // After
    if (z + depth < ns) {
      fill((z + depth) * nrByNc, (ns - z - depth) * nrByNc, value);
    }

    final int y_p_h = y + height;
    final int fillYBefore = y * nc;
    final int yAfter = y_p_h * nc;
    final int fillYAfter = (nr - y_p_h) * nc;

    final int x_p_w = x + width;
    final int fillXBefore = x;
    final int fillXAfter = (nc - x_p_w);

    for (int s = 0; s < depth; s++, z++) {
      int base = z * nrByNc;

      if (fillYBefore != 0) {
        fill(base, fillYBefore, value);
      }
      if (fillYAfter != 0) {
        fill(base + yAfter, fillYAfter, value);
      }

      base += fillYBefore;

      for (int r = 0; r < height; r++) {
        if (fillXBefore != 0) {
          fill(base, fillXBefore, value);
        }
        if (fillXAfter != 0) {
          fill(base + x_p_w, fillXAfter, value);
        }
        base += nc;
      }
    }
  }
}
