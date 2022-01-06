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

import ij.process.ImageProcessor;

/**
 * Store a 2D image in a single double array. Forms a base for 2D DHT transform using the
 * JTransforms library.
 */
public class DoubleImage2D extends Image2D {
  /**
   * The data packed in yx order: y * nc + x .
   *
   * @see #index(int, int)
   */
  protected double[] data;

  /**
   * Instantiates a new 2D image.
   *
   * @param nc the number of columns
   * @param nr the number of rows
   * @throws IllegalArgumentException If the combined dimensions is too large for an array
   */
  public DoubleImage2D(int nc, int nr) {
    super(nc, nr);
  }

  /**
   * Instantiates a new 2D image.
   *
   * @param image the image
   * @throws IllegalArgumentException If the combined dimensions is too large for an array
   */
  public DoubleImage2D(ImageProcessor image) {
    super(image);
  }

  /**
   * Instantiates a new 2D image.
   *
   * @param nc the number of columns
   * @param nr the number of rows
   * @param data the data
   * @throws IllegalArgumentException If the data is not the correct length
   */
  public DoubleImage2D(int nc, int nr, double[] data) {
    // Avoid constructor that calls createData(int)
    super(nc, nr, false);
    if (data == null || data.length != checkSize(nc, nr, true)) {
      throw new IllegalArgumentException("Data is not correct length");
    }
    this.data = data;
  }

  /**
   * Instantiates a new 2D image.
   *
   * @param nc the number of columns
   * @param nr the number of rows
   * @param data the data
   * @param dummy the dummy flag
   */
  protected DoubleImage2D(int nc, int nr, double[] data, boolean dummy) {
    // No checks as this is used internally
    super(nc, nr, false);
    this.data = data;
  }

  /**
   * Copy constructor.
   *
   * @param source the source
   */
  protected DoubleImage2D(DoubleImage2D source) {
    super(source);
    this.data = source.data.clone();
  }

  @Override
  public DoubleImage2D copy() {
    return new DoubleImage2D(this);
  }

  @Override
  protected void createData(int size) {
    data = new double[size];
  }

  /**
   * Gets the data.
   *
   * @return the data
   */
  public double[] getData() {
    return data;
  }

  @Override
  public int getDataLength() {
    return data.length;
  }

  @Override
  public DoubleImage2D crop(int x, int y, int width, int height) {
    return crop(x, y, width, height, null);
  }

  /**
   * Crop a sub-region of the data. The target dimensions must be positive.
   *
   * @param x the x index
   * @param y the y index
   * @param width the width
   * @param height the height
   * @param region the cropped data (will be reused if the correct size)
   * @return the cropped data
   * @throws IllegalArgumentException if the region is not within the data
   */
  public DoubleImage2D crop(int x, int y, int width, int height, double[] region) {
    // Check the region range
    if (x < 0 || width < 1 || (long) x + width > nc || y < 0 || height < 1
        || (long) y + height > nr) {
      throw new IllegalArgumentException("Region not within the data");
    }
    final int size = height * width;
    if (region == null || region.length != size) {
      region = new double[size];
    }
    int base = y * nc + x;
    for (int r = 0, i = 0; r < height; r++) {
      System.arraycopy(data, base, region, i, width);
      base += nc;
      i += width;
    }
    return new DoubleImage2D(width, height, region, false);
  }

  /**
   * Crop a sub-region of the data. The target dimensions must be positive.
   *
   * @param image the image
   * @param x the x index
   * @param y the y index
   * @param width the width
   * @param height the height
   * @param region the cropped data (will be reused if the correct size)
   * @return the cropped data
   * @throws IllegalArgumentException if the region is not within the data
   */
  public static DoubleImage2D crop(ImageProcessor image, int x, int y, int width, int height,
      double[] region) {
    final int nc = image.getWidth();
    final int nr = image.getHeight();

    // Check the region range
    if (x < 0 || width < 1 || (long) x + width > nc || y < 0 || height < 1
        || (long) y + height > nr) {
      throw new IllegalArgumentException("Region not within the data");
    }
    final int size = checkSize(width, height, true);
    if (region == null || region.length != size) {
      region = new double[size];
    }
    int base = y * nc + x;
    for (int r = 0, i = 0; r < height; r++) {
      for (int c = 0; c < width; c++) {
        region[i++] = image.getf(base + c);
      }
      base += nc;
    }
    return new DoubleImage2D(width, height, region, false);
  }

  @Override
  public void insert(int x, int y, Image2D image) {
    if (image instanceof DoubleImage2D) {
      insert(x, y, (DoubleImage2D) image);
    } else {
      super.insert(x, y, image);
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
  public void insert(int x, int y, DoubleImage2D image) {
    // Check the region range
    final int width = image.getWidth();
    final int height = image.getHeight();
    if (width < 1 || height < 1) {
      return;
    }
    if (x < 0 || (long) x + width > nc || y < 0 || (long) y + height > nr) {
      throw new IllegalArgumentException("Region not within the data");
    }
    final double[] region = image.data;
    int base = y * nc + x;
    for (int r = 0, i = 0; r < height; r++) {
      System.arraycopy(region, i, data, base, width);
      base += nc;
      i += width;
    }
  }

  @Override
  protected void copyTo(int index, float[] buffer, int bufferIndex, int size) {
    while (size-- > 0) {
      buffer[bufferIndex++] = (float) data[index++];
    }
  }

  @Override
  protected void copyFrom(float[] buffer, int bufferIndex, int size, int index) {
    while (size-- > 0) {
      data[index++] = buffer[bufferIndex++];
    }
  }

  @Override
  public double get(int index) {
    return data[index];
  }

  @Override
  public void set(int index, double value) {
    data[index] = value;
  }

  @Override
  public float getf(int index) {
    return (float) data[index];
  }

  @Override
  public void setf(int index, float value) {
    data[index] = value;
  }

  @Override
  protected void fill(int index, int size, double value) {
    while (size-- > 0) {
      data[index++] = value;
    }
  }
}
