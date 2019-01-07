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

import ij.ImageStack;
import ij.process.ImageProcessor;

/**
 * Store a 3D image in a single double array. Forms a base for 3D DHT transform using the
 * JTransforms library.
 */
public class DoubleImage3D extends Image3D {
  /**
   * The data packed in zyx order: z * (nr * nc) + y * nc + x .
   *
   * @see #index(int, int, int)
   */
  protected double[] data;

  /**
   * Instantiates a new 3D image.
   *
   * @param nc the number of columns
   * @param nr the number of rows
   * @param ns the number of slices
   * @throws IllegalArgumentException If the combined dimensions is too large for an array
   */
  public DoubleImage3D(int nc, int nr, int ns) {
    super(nc, nr, ns);
  }

  /**
   * Instantiates a new 3D image.
   *
   * @param stack the stack
   * @throws IllegalArgumentException If the combined dimensions is too large for an array
   */
  public DoubleImage3D(ImageStack stack) {
    super(stack);
  }

  /**
   * Instantiates a new 3D image.
   *
   * @param nc the number of columns
   * @param nr the number of rows
   * @param ns the number of slices
   * @param data the data
   * @throws IllegalArgumentException If the data is not the correct length
   */
  public DoubleImage3D(int nc, int nr, int ns, double[] data) {
    // Avoid constructor that calls createData(int)
    super(nc, nr, ns, nr * nc);
    if (data == null || data.length != checkSize(nc, nr, ns, true)) {
      throw new IllegalArgumentException("Data is not correct length");
    }
    this.data = data;
  }

  /**
   * Instantiates a new 3D image.
   *
   * @param nc the number of columns
   * @param nr the number of rows
   * @param ns the number of slices
   * @param nrByNc the number of rows multiplied by the number of columns
   * @param data the data
   */
  protected DoubleImage3D(int nc, int nr, int ns, int nrByNc, double[] data) {
    // No checks as this is used internally
    super(nc, nr, ns, nrByNc);
    this.data = data;
  }

  /**
   * Copy constructor.
   *
   * @param source the source
   */
  protected DoubleImage3D(DoubleImage3D source) {
    super(source);
    this.data = source.data.clone();
  }

  /** {@inheritDoc} */
  @Override
  public DoubleImage3D copy() {
    return new DoubleImage3D(this);
  }

  /** {@inheritDoc} */
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

  /** {@inheritDoc} */
  @Override
  public int getDataLength() {
    return data.length;
  }

  /** {@inheritDoc} */
  @Override
  public DoubleImage3D crop(int x, int y, int z, int width, int height, int depth) {
    return crop(x, y, z, width, height, depth, null);
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
   * @param region the cropped data (will be reused if the correct size)
   * @return the cropped data
   * @throws IllegalArgumentException if the region is not within the data
   */
  public DoubleImage3D crop(int x, int y, int z, int width, int height, int depth,
      double[] region) {
    // Check the region range
    if (x < 0 || width < 1 || (long) x + width > nc || y < 0 || height < 1 || (long) y + height > nr
        || z < 0 || depth < 1 || (long) z + depth > ns) {
      throw new IllegalArgumentException("Region not within the data");
    }
    final int size = depth * height * width;
    if (region == null || region.length != size) {
      region = new double[size];
    }
    for (int s = 0, i = 0; s < depth; s++, z++) {
      int base = z * nrByNc + y * nc + x;
      for (int r = 0; r < height; r++) {
        System.arraycopy(data, base, region, i, width);
        base += nc;
        i += width;
      }
    }
    return new DoubleImage3D(width, height, depth, width * height, region);
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
   * @param region the cropped data (will be reused if the correct size)
   * @return the cropped data
   * @throws IllegalArgumentException if the region is not within the data
   */
  public static DoubleImage3D crop(ImageStack stack, int x, int y, int z, int width, int height,
      int depth, double[] region) {
    final int nc = stack.getWidth();
    final int nr = stack.getHeight();
    final int ns = stack.getSize();

    // Check the region range
    if (x < 0 || width < 1 || (long) x + width > nc || y < 0 || height < 1 || (long) y + height > nr
        || z < 0 || depth < 1 || (long) z + depth > ns) {
      throw new IllegalArgumentException("Region not within the data");
    }
    final int size = checkSize(width, height, depth, true);
    if (region == null || region.length != size) {
      region = new double[size];
    }
    for (int s = 0, i = 0; s < depth; s++, z++) {
      final ImageProcessor ip = stack.getProcessor(1 + z);
      int base = y * nc + x;
      for (int r = 0; r < height; r++) {
        for (int c = 0; c < width; c++) {
          region[i++] = ip.getf(base + c);
        }
        base += nc;
      }
    }
    return new DoubleImage3D(width, height, depth, width * height, region);
  }

  @Override
  public void insert(int x, int y, int z, Image3D image) {
    if (image instanceof DoubleImage3D) {
      insert(x, y, z, (DoubleImage3D) image);
    } else {
      super.insert(x, y, z, image);
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
  public void insert(int x, int y, int z, DoubleImage3D image) {
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
    final double[] region = image.data;
    for (int s = 0, i = 0; s < depth; s++, z++) {
      int base = z * nrByNc + y * nc + x;
      for (int r = 0; r < height; r++) {
        System.arraycopy(region, i, data, base, width);
        base += nc;
        i += width;
      }
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
