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
import ij.process.ImageProcessor;
import uk.ac.sussex.gdsc.core.utils.MathUtils;

/**
 * Store a 2D image in a single array.
 */
public abstract class Image2D {
  /**
   * The largest array size for which a regular 1D Java array is used to store the data (2^30).
   */
  public static final int MAX_SIZE_OF_32_BIT_ARRAY = 1073741824;

  /** The number of rows (max y). */
  public final int nr;
  /** The number of columns (max x). */
  public final int nc;

  /**
   * Instantiates a new 2D image.
   *
   * @param nc the number of columns
   * @param nr the number of rows
   * @throws IllegalArgumentException If the combined dimensions is too large for an array
   */
  public Image2D(int nc, int nr) {
    createData(checkSize(nc, nr, true));
    this.nc = nc;
    this.nr = nr;
  }

  /**
   * Instantiates a new 2D image.
   *
   * @param image the image
   * @throws IllegalArgumentException If the combined dimensions is too large for an array
   */
  public Image2D(ImageProcessor image) {
    nc = image.getWidth();
    nr = image.getHeight();
    createData(checkSize(nc, nr, true));
    if (image.getBitDepth() == 32) {
      copyFrom((float[]) image.getPixels(), 0, nr * nc, 0);
    } else {
      for (int i = 0, size = nr * nc; i < size; i++) {
        setf(i, image.getf(i));
      }
    }
  }

  /**
   * Instantiates a new 2D image.
   *
   * <p>It is assumed that the sub-class will correctly create the data storage.
   *
   * @param nc the number of columns
   * @param nr the number of rows
   * @param dummy the dummy flag
   */
  protected Image2D(int nc, int nr, boolean dummy) {
    // No checks as this is used internally
    this.nc = nc;
    this.nr = nr;
  }

  /**
   * Copy constructor.
   *
   * <p>It is assumed that the sub-class will correctly create the data storage.
   *
   * @param source the source
   */
  protected Image2D(Image2D source) {
    nc = source.nc;
    nr = source.nr;
  }

  /**
   * Check the size can fit in a 1D array.
   *
   * @param nc the number of columns
   * @param nr the number of rows
   * @param raiseException Set to true to raise an exception if too large
   * @return the size, or -1 if too large
   * @throws IllegalArgumentException if too large (optional)
   */
  public static int checkSize(int nc, int nr, boolean raiseException) {
    if (nc < 0 || nr < 0) {
      if (raiseException) {
        throw new IllegalArgumentException("Negative dimensions");
      }
      return -1;
    }
    final long size = (long) nr * nc;
    if (size > MAX_SIZE_OF_32_BIT_ARRAY) {
      if (raiseException) {
        throw new IllegalArgumentException("2D data too large");
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
   * Return a copy of the 2D image.
   *
   * @return the copy
   */
  public abstract Image2D copy();

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
   * Gets the data length.
   *
   * @return the data length
   */
  public int getDataLength() {
    return nr * nc;
  }

  /**
   * Convert to an image processor.
   *
   * @return the image processor
   */
  public ImageProcessor getImageProcessor() {
    final float[] pixels = new float[getDataLength()];
    copyTo(0, pixels, 0, pixels.length);
    return new FloatProcessor(nc, nr, pixels);
  }

  /**
   * Gets the xy components of the index.
   *
   * @param index the index
   * @return the xy components
   * @throws IllegalArgumentException if the index is not within the data
   */
  public int[] getXy(int index) {
    if (index < 0 || index >= getDataLength()) {
      throw new IllegalArgumentException(
          "Index in not in the correct range: 0 <= i < " + getDataLength());
    }
    final int[] xyz = new int[2];
    xyz[1] = index / nc;
    xyz[0] = index % nc;
    return xyz;
  }

  /**
   * Gets the xy components of the index.
   *
   * @param index the index
   * @param xy the xy components (must be an array of at least length 2)
   * @throws IllegalArgumentException if the index is not within the data
   */
  public void getXy(int index, int[] xy) {
    if (index < 0 || index >= getDataLength()) {
      throw new IllegalArgumentException(
          "Index in not in the correct range: 0 <= i < " + getDataLength());
    }
    xy[1] = index / nc;
    xy[0] = index % nc;
  }

  /**
   * Gets the index using the xy components.
   *
   * @param x the x component
   * @param y the y component
   * @return the index
   * @throws IllegalArgumentException if the index is not within the data
   */
  public int getIndex(int x, int y) {
    if (x < 0 || x >= nc || y < 0 || y >= nr) {
      throw new IllegalArgumentException("Index in not inside the image");
    }
    return y * nc + x;
  }

  /**
   * Gets the index using the xy components.
   *
   * @param x the x component
   * @param y the y component
   * @return the index
   */
  protected int index(int x, int y) {
    return y * nc + x;
  }

  /**
   * Crop a sub-region of the data. The target dimensions must be positive.
   *
   * @param x the x index
   * @param y the y index
   * @param width the width
   * @param height the height
   * @return the cropped data
   * @throws IllegalArgumentException if the region is not within the data
   */
  public abstract Image2D crop(int x, int y, int width, int height);

  /**
   * Crop a sub-region of the data into the given image. The target dimensions must be positive.
   *
   * @param x the x index
   * @param y the y index
   * @param image the image
   * @return the cropped data
   * @throws IllegalArgumentException if the region is not within the data
   */
  public Image2D crop(int x, int y, Image2D image) {
    // Check the region range
    final int width = image.getWidth();
    final int height = image.getHeight();
    if (x < 0 || width < 1 || (long) x + width > nc || y < 0 || height < 1
        || (long) y + height > nr) {
      throw new IllegalArgumentException("Region not within the data");
    }
    int base = y * nc + x;
    for (int r = 0, i = 0; r < height; r++) {
      for (int c = 0; c < width; c++) {
        image.set(i++, get(base + c));
      }
      base += nc;
    }
    return image;
  }

  /**
   * Crop a sub-region of the data. The target dimensions must be positive.
   *
   * @param x the x index
   * @param y the y index
   * @param width the width
   * @param height the height
   * @return the cropped data
   * @throws IllegalArgumentException if the region is not within the data
   */
  public ImageProcessor cropToProcessor(int x, int y, int width, int height) {
    // Check the region range
    if (x < 0 || width < 1 || (long) x + width > nc || y < 0 || height < 1
        || (long) y + height > nr) {
      throw new IllegalArgumentException("Region not within the data");
    }
    final int size = width * height;
    int base = y * nc + x;
    final float[] region = new float[size];
    for (int r = 0, i = 0; r < height; r++) {
      copyTo(base, region, i, width);
      base += nc;
      i += width;
    }
    return new FloatProcessor(width, height, region);
  }

  /**
   * Insert a sub-region.
   *
   * @param x the x position
   * @param y the y position
   * @param image the image
   * @throws IllegalArgumentException if the region is not within the data
   */
  public void insert(int x, int y, Image2D image) {
    // Check the region range
    final int width = image.getWidth();
    final int height = image.getHeight();
    if (width < 1 || height < 1) {
      return;
    }
    if (x < 0 || (long) x + width > nc || y < 0 || (long) y + height > nr) {
      throw new IllegalArgumentException("Region not within the data");
    }
    int base = y * nc + x;
    for (int r = 0, i = 0; r < height; r++) {
      for (int c = 0; c < width; c++) {
        set(base + c, image.get(i++));
      }
      base += nc;
    }
  }

  /**
   * Insert a sub-region.
   *
   * @param x the x position
   * @param y the y position
   * @param image the image
   * @throws IllegalArgumentException if the region is not within the data
   */
  public void insert(int x, int y, ImageProcessor image) {
    // Check the region range
    final int width = image.getWidth();
    final int height = image.getHeight();
    if (width < 1 || height < 1) {
      return;
    }
    if (x < 0 || (long) x + width > nc || y < 0 || (long) y + height > nr) {
      throw new IllegalArgumentException("Region not within the data");
    }
    final boolean isFloat = image.getBitDepth() == 32;
    int base = y * nc + x;
    final float[] region =
        (float[]) ((isFloat) ? image.getPixels() : image.toFloat(0, null).getPixels());
    for (int r = 0, i = 0; r < height; r++) {
      copyFrom(region, i, width, base);
      base += nc;
      i += width;
    }
  }

  /**
   * Compute 2D intersect with this object.
   *
   * <p>If any of width,height are negative then the corresponding x,y is updated and the
   * width,height is inverted. The maximum bounds of the given dimensions are then computed by
   * adding the width,height to the x,y. The bounds are then clipped to the image dimensions and the
   * intersect returned.
   *
   * @param x the x index
   * @param y the y index
   * @param width the width
   * @param height the height
   * @return [x,y,width,height]
   */
  public int[] computeIntersect(int x, int y, int width, int height) {
    if (width < 0) {
      width = -width;
      x = subtract(x, width);
    }
    if (height < 0) {
      height = -height;
      y = subtract(y, height);
    }
    // Compute 2D intersect with this object
    final int x2 = clip(nc, x, width);
    final int y2 = clip(nr, y, height);
    x = clip(nc, x);
    y = clip(nr, y);
    return new int[] {x, y, x2 - x, y2 - y};
  }

  /**
   * Compute 2D intersect with this object or throw an exception if the intersect has no volume.
   *
   * <p>If any of width,height are negative then the corresponding x,y is updated and the
   * width,height is inverted. The maximum bounds of the given dimensions are then computed by
   * adding the width,height to the x,y. The bounds are then clipped to the image dimensions and the
   * intersect returned.
   *
   * @param x the x index
   * @param y the y index
   * @param width the width
   * @param height the height
   * @return [x,y,width,height]
   * @throws IllegalArgumentException if the intersect has no volume
   */
  public int[] computeIntersectOrThrow(int x, int y, int width, int height) {
    if (width < 0) {
      width = -width;
      x = subtract(x, width);
    }
    if (height < 0) {
      height = -height;
      y = subtract(y, height);
    }
    // Compute 2D intersect with this object
    final int x2 = clip(nc, x, width);
    final int y2 = clip(nr, y, height);
    x = clip(nc, x);
    y = clip(nr, y);
    width = checkIntersectSize(x2 - x);
    height = checkIntersectSize(y2 - y);
    return new int[] {x, y, width, height};
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
   * @param width the width
   * @param height the height
   * @return the min index
   * @throws IllegalArgumentException if there is no intersect
   * @throws IllegalStateException if the region is just NaN values
   */
  public int findMinIndex(int x, int y, int width, int height) {
    final int[] intersect = computeIntersectOrThrow(x, y, width, height);
    x = intersect[0];
    y = intersect[1];
    width = intersect[2];
    height = intersect[3];
    int index = findValueIndex(x, y, width, height);
    double min = get(index);
    int base = y * nc + x;
    for (int r = 0; r < height; r++) {
      for (int j = 0; j < width; j++) {
        if (get(base + j) < min) {
          index = base + j;
          min = get(index);
        }
      }
      base += nc;
    }
    return index;
  }

  /**
   * Find the index of the minimum value in the region.
   *
   * @param x the x index
   * @param y the y index
   * @param width the width
   * @param height the height
   * @return the min index
   * @throws IllegalArgumentException if there is no intersect
   * @throws IllegalStateException if the region is just NaN values
   */
  public int findMaxIndex(int x, int y, int width, int height) {
    final int[] intersect = computeIntersectOrThrow(x, y, width, height);
    x = intersect[0];
    y = intersect[1];
    width = intersect[2];
    height = intersect[3];
    int index = findValueIndex(x, y, width, height);
    double min = get(index);
    int base = y * nc + x;
    for (int r = 0; r < height; r++) {
      for (int j = 0; j < width; j++) {
        if (get(base + j) > min) {
          index = base + j;
          min = get(index);
        }
      }
      base += nc;
    }
    return index;
  }

  /**
   * Find the index of the first non-NaN value in the region.
   *
   * @param x the x index
   * @param y the y index
   * @param width the width
   * @param height the height
   * @return the index
   * @throws IllegalStateException if the region is just NaN values
   */
  private int findValueIndex(int x, int y, int width, int height) {
    // Quick check without loops
    if (!Double.isNaN(get(y * nc + x))) {
      return y * nc + x;
    }

    int base = y * nc + x;
    for (int r = 0; r < height; r++) {
      for (int j = 0; j < width; j++) {
        if (!Double.isNaN(get(base + j))) {
          return base + j;
        }
      }
      base += nc;
    }
    throw new IllegalStateException("Region is NaN");
  }

  /**
   * Compute the rolling sum table for use in {@link #computeSum(double[], int, int, int, int)}.
   *
   * <p>This is a table of the sum of the volume from 0,0 to x,y inclusive.
   *
   * @param table the table (will be reused if the correct size)
   * @return the rolling sum table
   */
  public double[] computeRollingSumTable(double[] table) {
    if (table == null || table.length != getDataLength()) {
      table = new double[getDataLength()];
    }

    double sum = 0;
    int index = 0;
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

    return table;
  }

  /**
   * Compute the sum of the region.
   *
   * @param x the x index
   * @param y the y index
   * @param width the width
   * @param height the height
   * @return the sum
   */
  public double computeSum(int x, int y, int width, int height) {
    final int[] intersect = computeIntersect(x, y, width, height);
    width = intersect[2];
    height = intersect[3];
    // Recheck bounds
    if (width == 0 || height == 0) {
      return 0;
    }
    x = intersect[0];
    y = intersect[1];
    double sum = 0;
    int base = y * nc + x;
    for (int r = 0; r < height; r++) {
      for (int j = 0; j < width; j++) {
        sum += get(base + j);
      }
      base += nc;
    }
    return sum;
  }

  /**
   * Compute the sum of the region using the precomputed rolling sum table.
   *
   * @param table the rolling sum table
   * @param x the x index
   * @param y the y index
   * @param width the width
   * @param height the height
   * @return the sum
   */
  public double computeSum(double[] table, int x, int y, int width, int height) {
    final int[] intersect = computeIntersect(x, y, width, height);
    width = intersect[2];
    height = intersect[3];
    // Recheck bounds
    if (width == 0 || height == 0) {
      return 0;
      // x = intersect[0];
      // y = intersect[1];
    }

    // Compute sum from rolling sum using:
    // sum(x,y,width,d) =
    // + s(x+width-1,y+height-1)
    // - s(x-1,y+height-1)
    // - s(x+width-1,y-1)
    // + s(x-1,y-1)
    // Note:
    // s(i,j) = 0 when either i,j < 0
    // i = imax when i>imax
    // j = jmax when j>jmax

    final int x_1 = intersect[0] - 1;
    final int y_1 = intersect[1] - 1;
    // The intersect has already checked the bounds
    // int x_w_1 = Math.min(x_1 + width, nc);
    // int y_h_1 = Math.min(y_1 + height, nr);
    final int x_w_1 = x_1 + width;
    final int y_h_1 = y_1 + height;

    // double sum = table[index(x_w_1, y_h_1)];
    // if (y_1 >= 0)
    // {
    // sum -= table[index(x_w_1, y_1)];
    // if (x_1 >= 0)
    // sum = sum + table[index(x_1, y_1)] - table[index(x_1, y_h_1)];
    // }
    // else if (x_1 >= 0)
    // {
    // sum -= table[index(x_1, y_h_1)];
    // }
    // return sum;

    // This has been ordered to use the smallest sums first (i.e. closer to x,y than
    // x+width,y+height)
    final int xw_yh = index(x_w_1, y_h_1);
    double sum = 0;
    if (y_1 >= 0) {
      final int h_ = height * nc;
      if (x_1 >= 0) {
        sum = table[xw_yh - width - h_] - table[xw_yh - width];
      }
      sum -= table[xw_yh - h_];
    } else if (x_1 >= 0) {
      sum = -table[xw_yh - width];
    }
    return sum + table[xw_yh];
  }

  /**
   * Compute the sum of the region using the precomputed rolling sum table. Assumes
   * x+width,y+height,z+d will not overflow!
   *
   * @param table the rolling sum table
   * @param x the x index
   * @param y the y index
   * @param width the width
   * @param height the height
   * @return the sum
   */
  public double computeSumFast(double[] table, int x, int y, int width, int height) {
    if (width <= 0 || height <= 0 || x >= nc || y >= nr) {
      return 0;
    }

    // Compute sum from rolling sum using:
    // sum(x,y,width,d) =
    // + s(x+width-1,y+height-1)
    // - s(x-1,y+height-1)
    // - s(x+width-1,y-1)
    // + s(x-1,y-1)
    // Note:
    // s(i,j) = 0 when either i,j < 0
    // i = imax when i>imax
    // j = jmax when j>jmax

    // Compute bounds assuming width,d is small and positive.
    int x1;
    int y1;
    int xw1;
    int yh1;
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

    // Adjust for the -1
    x1--;
    y1--;
    xw1--;
    yh1--;

    // This has been ordered to use the smallest sums first (i.e. closer to x,y than
    // x+width,y+height)
    final int xw1yh1 = index(xw1, yh1);
    double sum = 0;
    if (y1 >= 0) {
      final int heightOffset = height * nc;
      if (x1 >= 0) {
        sum = table[xw1yh1 - width - heightOffset] - table[xw1yh1 - width];
      }
      sum -= table[xw1yh1 - heightOffset];
    } else if (x1 >= 0) {
      sum = -table[xw1yh1 - width];
    }
    return sum + table[xw1yh1];
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
   * @param width the width
   * @param height the height
   * @param value the value
   */
  public void fill(int x, int y, int width, int height, double value) {
    final int[] intersect = computeIntersect(x, y, width, height);
    width = intersect[2];
    height = intersect[3];
    // Recheck bounds
    if (width == 0 || height == 0) {
      return;
    }
    x = intersect[0];
    y = intersect[1];
    int base = y * nc + x;
    for (int r = 0; r < height; r++) {
      fill(base, width, value);
      base += nc;
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
   * @param width the width
   * @param height the height
   * @param value the value
   */
  public void fillOutside(int x, int y, int width, int height, double value) {
    final int[] intersect = computeIntersect(x, y, width, height);
    width = intersect[2];
    height = intersect[3];
    // Recheck bounds
    if (width == 0 || height == 0) {
      fill(value);
      return;
    }
    x = intersect[0];
    y = intersect[1];

    final int y_p_h = y + height;
    final int fillYBefore = y * nc;
    final int fillYAfter = (nr - y_p_h) * nc;

    final int x_p_w = x + width;
    final int fillXBefore = x;
    final int fillXAfter = (nc - x_p_w);

    if (fillYBefore != 0) {
      fill(0, fillYBefore, value);
    }
    if (fillYAfter != 0) {
      fill(y_p_h * nc, fillYAfter, value);
    }

    int base = fillYBefore;

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
