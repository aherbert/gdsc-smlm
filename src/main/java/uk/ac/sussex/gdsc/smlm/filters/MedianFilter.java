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

import org.apache.commons.math3.util.FastMath;
import uk.ac.sussex.gdsc.core.utils.FloatLinkedMedianWindow;

/**
 * Computes the block median for each point within the array.
 *
 * <p>block algorithm sweeps the entire (2n+1)*(2n+1) region around each pixel.
 *
 * <p>Note: Due to lack of small dimension checking the routines will fail if maxx or maxy are less
 * than 2. All routines are OK for 3x3 images and larger.
 */
public class MedianFilter {
  private float[] floatDataBuffer;

  /** The number of values above the guess. */
  private int above;
  /** The number of values below the guess. */
  private int below;
  /** Half of the total number of values that will be added to the data. */
  private int half;
  /** The guess for the median. */
  private float guess;
  /** Buffer to store all values above the guess. */
  private float[] aboveBuf;
  /** Buffer to store all values below the guess. */
  private float[] belowBuf;

  /**
   * Create a copy.
   *
   * @return the copy
   */
  public MedianFilter copy() {
    return new MedianFilter();
  }

  /**
   * Compute the block median within a 2n+1 size block around each point. Only pixels with a full
   * block are processed. Pixels within border regions are unchanged.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   * @param n The block size
   */
  public void blockMedianInternal(float[] data, final int maxx, final int maxy, final int n) {
    if (n == 1) {
      blockMedian3x3Internal(data, maxx, maxy);
    } else {
      blockMedianNxNInternal(data, maxx, maxy, n);
    }
  }

  /**
   * Compute the block median within a 2n+1 size block around each point. Only pixels with a full
   * block are processed. Pixels within border regions are unchanged.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   * @param n The block size
   */
  public void blockMedianNxNInternal(float[] data, final int maxx, final int maxy, final int n) {
    final int blockSize = 2 * n + 1;
    if (maxx < blockSize || maxy < blockSize) {
      return;
    }

    floatDataBuffer = floatBuffer(floatDataBuffer, data.length);

    final int[] offsets = new int[blockSize * blockSize - 1];
    for (int y = -n, d = 0; y <= n; y++) {
      for (int x = -n; x <= n; x++) {
        if (x != 0 || y != 0) {
          offsets[d] = maxx * y + x;
          d++;
        }
      }
    }

    init(offsets.length + 1, data[n * maxx + n]);

    for (int y = n; y < maxy - n; y++) {
      int index = y * maxx + n;
      for (int x = n; x < maxx - n; x++, index++) {
        reset();
        add(data[index]);

        // Sweep neighbourhood -
        // No check for boundaries as this should be an internal sweep.
        for (final int offset : offsets) {
          add(data[index + offset]);
        }

        floatDataBuffer[index] = getMedian();
      }
    }

    // Copy back
    for (int y = n; y < maxy - n; y++) {
      int index = y * maxx + n;
      for (int x = n; x < maxx - n; x++, index++) {
        data[index] = floatDataBuffer[index];
      }
    }
  }

  /**
   * Compute the block median within a 3x3 size block around each point. Only pixels with a full
   * block are processed. Pixels within border regions are unchanged.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   */
  public void blockMedian3x3Internal(float[] data, final int maxx, final int maxy) {
    final float[] newData = floatBuffer(floatDataBuffer, data.length);
    init(9, data[maxx + 1]);
    // Boundary control
    final int xlimit = maxx - 1;
    final int ylimit = maxy - 1;
    for (int y = 1; y < ylimit; y++) {
      int index1 = y * maxx + 1;
      int index0 = index1 - maxx;
      int index2 = index1 + maxx;
      for (int x = 1; x < xlimit; x++) {
        reset();
        add(data[index0 - 1]);
        add(data[index0]);
        add(data[index0 + 1]);
        add(data[index1 - 1]);
        add(data[index1]);
        add(data[index1 + 1]);
        add(data[index2 - 1]);
        add(data[index2]);
        add(data[index2 + 1]);
        newData[index1] = getMedian();
        index0++;
        index1++;
        index2++;
      }
    }

    // Copy back
    for (int y = 1; y < ylimit; y++) {
      int index = y * maxx + 1;
      for (int x = 1; x < xlimit; x++, index++) {
        data[index] = newData[index];
      }
    }
  }

  /**
   * Initialise the median buffers.
   *
   * @param size The total number of values to add (should be odd)
   * @param guess The guess for the median
   */
  private void init(int size, float guess) {
    if (aboveBuf == null || aboveBuf.length < size) {
      aboveBuf = new float[size];
      belowBuf = new float[size];
    }
    half = size / 2;
  }

  /**
   * Reset the median buffer counters.
   */
  private void reset() {
    above = below = 0;
  }

  /**
   * Add a value.
   *
   * @param value the vakue
   */
  private void add(float value) {
    if (value > guess) {
      aboveBuf[above++] = value;
    } else if (value < guess) {
      belowBuf[below++] = value;
    }
  }

  /**
   * Get median of values.
   *
   * @return The median
   */
  private float getMedian() {
    if (above > half) {
      guess = findNthLowestNumber(aboveBuf, above, above - half - 1);
    } else if (below > half) {
      guess = findNthLowestNumber(belowBuf, below, half);
    }

    // Debug. This is just a check and can be removed for production code
    // if (nAbove + nBelow == 2 * half + 1)
    // {
    // float[] values = new float[nAbove + nBelow];
    // for (int i = 0; i < nAbove; i++)
    // values[i] = aboveBuf[i];
    // for (int i = 0, j = nAbove; i < nBelow; i++, j++)
    // values[j] = belowBuf[i];
    // java.util.Arrays.sort(values);
    // if (values[half] != guess)
    // System.out.printf("Mistake: %f != %f\n", values[half], guess);
    // }

    return guess;
  }

  /**
   * Find the n-th lowest number in part of an array.
   *
   * <p>Copied from ij.plugins.filters.RankFilters.
   *
   * @param buf The input array. Only values 0 ... bufLength are read. <code>buf</code> will be
   *        modified.
   * @param bufLength Number of values in <code>buf</code> that should be read
   * @param n which value should be found; n=0 for the lowest, n=bufLength-1 for the highest
   * @return the value
   */
  private static final float findNthLowestNumber(float[] buf, int bufLength, int n) {
    // Hoare's find, algorithm, based on http://www.geocities.com/zabrodskyvlada/3alg.html
    // Contributed by Heinz Klar.
    // CHECKSTYLE.OFF: LocalVariableName
    int i;
    int j;
    int l = 0;
    int m = bufLength - 1;
    // CHECKSTYLE.ON: LocalVariableName
    float med = buf[n];

    while (l < m) {
      i = l;
      j = m;
      do {
        while (buf[i] < med) {
          i++;
        }
        while (med < buf[j]) {
          j--;
        }
        final float dum = buf[j];
        buf[j] = buf[i];
        buf[i] = dum;
        i++;
        j--;
      } while ((j >= n) && (i <= n));
      if (j < n) {
        l = i;
      }
      if (n < i) {
        m = j;
      }
      med = buf[n];
    }
    return med;
  }

  private static float[] floatBuffer(float[] buffer, int size) {
    if (buffer == null || buffer.length < size) {
      buffer = new float[size];
    }
    return buffer;
  }

  /**
   * Compute the block median within a 2n+1 size block around each point.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   * @param n The block size
   */
  public void blockMedian(float[] data, final int maxx, final int maxy, final int n) {
    if (n == 1) {
      blockMedian3x3(data, maxx, maxy);
    } else {
      blockMedianNxN(data, maxx, maxy, n);
    }
  }

  /**
   * Compute the block median within a 2n+1 size block around each point.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   * @param n The block size
   */
  public void blockMedianNxN(float[] data, final int maxx, final int maxy, final int n) {
    final int length = maxx * maxy;
    final float[] newData = floatBuffer(floatDataBuffer, length);

    // Boundary control
    final int xwidth = FastMath.min(n, maxx - 1);
    final int ywidth = FastMath.min(n, maxy - 1);
    final int xlimit = maxx - xwidth;
    final int ylimit = maxy - ywidth;

    final int[] offsets = new int[(2 * xwidth + 1) * (2 * ywidth + 1) - 1];
    final int[] xoffset = new int[offsets.length];
    final int[] yoffset = new int[offsets.length];
    for (int y = -ywidth, d = 0; y <= ywidth; y++) {
      for (int x = -xwidth; x <= xwidth; x++) {
        if (x != 0 || y != 0) {
          offsets[d] = maxx * y + x;
          xoffset[d] = x;
          yoffset[d] = y;
          d++;
        }
      }
    }

    init(offsets.length + 1, data[0]);

    int index = 0;
    for (int y = 0; y < maxy; y++) {
      for (int x = 0; x < maxx; x++, index++) {
        reset();
        add(data[index]);

        // Flag to indicate this pixels has a complete (2n+1) neighbourhood
        final boolean isInnerXy = (y >= ywidth && y < ylimit) && (x >= xwidth && x < xlimit);

        // Sweep neighbourhood
        if (isInnerXy) {
          for (final int offset : offsets) {
            add(data[index + offset]);
          }
        } else {
          for (int d = offsets.length; d-- > 0;) {
            // Get the pixel with boundary checking
            int yy = y + yoffset[d];
            int xx = x + xoffset[d];
            if (xx <= 0) {
              xx = 0;
            } else if (xx >= maxx) {
              xx = maxx - 1;
            }
            if (yy <= 0) {
              yy = 0;
            } else if (yy >= maxy) {
              yy = maxy - 1;
            }
            add(data[xx + yy * maxx]);
          }
        }

        newData[index] = getMedian();
      }
    }

    // Copy back
    System.arraycopy(newData, 0, data, 0, length);
  }

  /**
   * Compute the block median within a 3x3 size block around each point.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   */
  public void blockMedian3x3(float[] data, final int maxx, final int maxy) {
    final int length = maxx * maxy;
    final float[] newData = floatBuffer(floatDataBuffer, length);
    init(9, data[maxx + 1]);

    // Boundary control
    final int xlimit = maxx - 1;
    final int ylimit = maxy - 1;

    final int[] xoffset = {-1, 0, 1, -1, 1, -1, 0, 1};
    final int[] yoffset = {-1, -1, -1, 0, 0, 1, 1, 1};

    for (int y = 0; y < maxy; y++) {
      int index1 = y * maxx;
      int index0 = index1 - maxx;
      int index2 = index1 + maxx;
      final boolean isInnerY = (y > 0 && y < ylimit);
      for (int x = 0; x < maxx; x++) {
        reset();
        add(data[index1]);

        // Sweep neighbourhood
        if (isInnerY && (x > 0 && x < xlimit)) {
          add(data[index0 - 1]);
          add(data[index0]);
          add(data[index0 + 1]);
          add(data[index1 - 1]);
          add(data[index1 + 1]);
          add(data[index2 - 1]);
          add(data[index2]);
          add(data[index2 + 1]);
        } else {
          for (int d = xoffset.length; d-- > 0;) {
            // Get the pixel with boundary checking
            int yy = y + yoffset[d];
            final int xx = x + xoffset[d];
            if (yy < 0) {
              yy = 0;
            } else if (yy == maxy) {
              yy = ylimit;
            }
            if (xx < 0) {
              add(data[yy * maxx]);
            } else if (xx == maxx) {
              add(data[xlimit + yy * maxx]);
            } else {
              add(data[xx + yy * maxx]);
            }
          }
        }
        newData[index1] = getMedian();
        index0++;
        index1++;
        index2++;
      }
    }

    // Copy back
    System.arraycopy(newData, 0, data, 0, length);
  }

  /**
   * Compute the rolling median within a 2n+1 size rolling around each point. Only pixels with a
   * full block are processed. Pixels within border regions are unchanged.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   * @param n The rolling size
   */
  public void rollingMedianInternal(float[] data, final int maxx, final int maxy, final int n) {
    if (n == 1) {
      rollingMedian3x3Internal(data, maxx, maxy);
    } else {
      rollingMedianNxNInternal(data, maxx, maxy, n);
    }
  }

  /**
   * Compute the rolling median within a 2n+1 size rolling around each point. Only pixels with a
   * full block are processed. Pixels within border regions are unchanged.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   * @param n The rolling size
   */
  public void rollingMedianNxNInternal(float[] data, final int maxx, final int maxy, final int n) {
    final int blockSize = 2 * n + 1;
    if (maxx < blockSize || maxy < blockSize) {
      return;
    }

    final int length = maxx * maxy;
    final float[] newData = floatBuffer(floatDataBuffer, length);

    // Hold the pointers to the image data for nY rows
    final int[] p = new int[blockSize];
    // Buffer to hold the initial region
    final float[] values = new float[blockSize * blockSize];
    for (int y = n; y < maxy - n; y++) {
      // Set up the pointers to the image data at x=0, y=?
      int vi = 0;
      for (int yy = y - n; yy <= y + n; yy++) {
        p[vi] = maxx * yy;
        // zero the first column of the region
        values[vi++] = 0;
      }

      // Fill the initial region

      for (int x = -n; x < n; x++) {
        for (int d = 0; d < p.length; d++) {
          values[vi++] = data[p[d]++];
        }
      }

      // Initialise the rolling window
      final FloatLinkedMedianWindow window = new FloatLinkedMedianWindow(values);

      // For each position up to the limit, add the next column and increment
      int index = y * maxx + n;
      for (int x = n; x < maxx - n; x++) {
        for (int j = 0; j < p.length; j++) {
          window.add(data[p[j]++]);
        }
        newData[index++] = window.getMedian();
      }
    }

    // Copy back
    for (int y = n; y < maxy - n; y++) {
      int index = y * maxx + n;
      for (int x = n; x < maxx - n; x++, index++) {
        data[index] = newData[index];
      }
    }
  }

  /**
   * Compute the rolling median within a 3x3 size rolling around each point. Only pixels with a full
   * block are processed. Pixels within border regions are unchanged.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   */
  public void rollingMedian3x3Internal(float[] data, final int maxx, final int maxy) {
    final int length = maxx * maxy;
    final float[] newData = floatBuffer(floatDataBuffer, length);

    // Boundary control
    final int xlimit = maxx - 1;
    final int ylimit = maxy - 1;

    // Buffer to hold the initial region
    final float[] values = new float[9];

    for (int y = 1; y < ylimit; y++) {

      values[0] = 0;
      values[1] = 0;
      values[2] = 0;

      int p1 = y * maxx;
      int p0 = p1 - maxx;
      int p2 = p1 + maxx;
      values[3] = data[p0++];
      values[4] = data[p1++];
      values[5] = data[p2++];
      values[6] = data[p0++];
      values[7] = data[p1++];
      values[8] = data[p2++];

      // Initialise the rolling window
      final FloatLinkedMedianWindow window = new FloatLinkedMedianWindow(values);

      // For each position up to the limit, add the next column and increment
      int index = p1 - 1;
      for (int x = 1; x < xlimit; x++) {
        window.add(data[p0++]);
        window.add(data[p1++]);
        window.add(data[p2++]);
        newData[index++] = window.getMedian();
      }
    }

    // Copy back
    for (int y = 1; y < ylimit; y++) {
      int index = y * maxx + 1;
      for (int x = 1; x < xlimit; x++, index++) {
        data[index] = newData[index];
      }
    }
  }

  /**
   * Compute the rolling median within a 2n+1 size rolling around each point.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   * @param n The rolling size
   */
  public void rollingMedian(float[] data, final int maxx, final int maxy, final int n) {
    if (n == 1) {
      rollingMedian3x3(data, maxx, maxy);
    } else {
      rollingMedianNxN(data, maxx, maxy, n);
    }
  }

  /**
   * Compute the rolling median within a 2n+1 size rolling around each point.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   * @param n The rolling size
   */
  public void rollingMedianNxN(float[] data, final int maxx, final int maxy, final int n) {
    final int length = maxx * maxy;
    final float[] newData = floatBuffer(floatDataBuffer, length);

    // Boundary control
    final int xwidth = FastMath.min(n, maxx - 1);
    final int ywidth = FastMath.min(n, maxy - 1);
    final int xlimit = maxx - xwidth - 1;

    // The size of the region
    final int nX = (2 * xwidth + 1);
    final int nY = (2 * ywidth + 1);

    // Hold the pointers to the image data for nY rows
    final int[] p = new int[nY];
    // Buffer to hold the initial region
    final float[] values = new float[nX * nY];
    int index = 0;
    for (int y = 0; y < maxy; y++) {
      // Set up the pointers to the image data at x=0, y=?
      int vi = 0;
      for (int yy = y - ywidth; yy <= y + ywidth; yy++) {
        if (yy < 0) {
          p[vi] = 0;
        } else if (yy >= maxy) {
          p[vi] = length - maxx;
        } else {
          p[vi] = maxx * yy;
        }
        // zero the first column of the region
        values[vi++] = 0;
      }

      // Fill the initial region

      // The columns below x==0 use x=0
      for (int x = -xwidth; x < 0; x++) {
        for (final int pos : p) {
          values[vi++] = data[pos];
        }
      }
      // The remaining columns increment. Do not include x==xwidth
      for (int x = 0; x < xwidth; x++) {
        for (int d = 0; d < p.length; d++) {
          values[vi++] = data[p[d]++];
        }
      }

      // Initialise the rolling window
      final FloatLinkedMedianWindow window = new FloatLinkedMedianWindow(values);

      // For each position up to the limit, add the next column and increment
      for (int x = 0; x < xlimit; x++) {
        for (int j = 0; j < p.length; j++) {
          window.add(data[p[j]++]);
        }
        newData[index++] = window.getMedian();
      }

      // Add the last column but do not increment
      for (int x = xlimit; x < maxx; x++) {
        for (int j = 0; j < p.length; j++) {
          window.add(data[p[j]]);
        }
        newData[index++] = window.getMedian();
      }
    }

    // Copy back
    System.arraycopy(newData, 0, data, 0, length);
  }

  /**
   * Compute the rolling median within a 3x3 size rolling around each point.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   */
  public void rollingMedian3x3(float[] data, final int maxx, final int maxy) {
    final int length = maxx * maxy;
    final float[] newData = floatBuffer(floatDataBuffer, length);

    // Boundary control
    final int xlimit = maxx - 2;
    final int ylimit = maxy - 1;

    // Buffer to hold the initial region
    final float[] values = new float[9];

    int index = 0;
    for (int y = 0; y < maxy; y++) {
      values[0] = 0;
      values[1] = 0;
      values[2] = 0;

      // Fill the initial region

      // Set up the pointers to the image data at x=0, y=?
      int p1 = maxx * y;
      int p0 = (y > 0) ? p1 - maxx : p1;
      int p2 = (y < ylimit) ? p1 + maxx : p1;

      // The columns below x==0 use x=0
      values[3] = data[p0];
      values[4] = data[p1];
      values[5] = data[p2];

      // The remaining columns increment. Do not include x==xwidth
      values[6] = data[p0++];
      values[7] = data[p1++];
      values[8] = data[p2++];

      // Initialise the rolling window
      final FloatLinkedMedianWindow window = new FloatLinkedMedianWindow(values);

      // For each position up to the limit, add the next column and increment
      for (int x = 0; x < xlimit; x++) {
        window.add(data[p0++]);
        window.add(data[p1++]);
        window.add(data[p2++]);
        newData[index++] = window.getMedian();
      }

      // Add the last column but do not increment
      for (int x = xlimit; x < maxx; x++) {
        window.add(data[p0]);
        window.add(data[p1]);
        window.add(data[p2]);
        newData[index++] = window.getMedian();
      }
    }

    // Copy back
    System.arraycopy(newData, 0, data, 0, length);
  }
}
