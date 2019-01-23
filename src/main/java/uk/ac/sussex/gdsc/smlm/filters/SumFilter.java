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

/**
 * Computes the block sum for each point within the array.
 *
 * @deprecated replaced by BlockSumFilter
 */
@Deprecated
public class SumFilter {
  private float[] floatDataBuffer;
  private float[] floatRowBuffer;

  /**
   * Create a copy.
   *
   * @return the copy
   */
  public SumFilter copy() {
    return new SumFilter();
  }

  /**
   * Compute the block sum within a 2n+1 size block around each point. Only pixels with a full block
   * are processed. Pixels within border regions are unchanged.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   * @param n The block size
   */
  public void rollingBlockSumInternal(float[] data, final int maxx, final int maxy, final int n) {
    // Note: Speed tests show that this method is only marginally faster than
    // rollingBlockSumNxNInternal.
    // Sometimes it is slower. The intricacies of the java optimiser escape me.
    if (n == 1) {
      rollingBlockSum3x3Internal(data, maxx, maxy);
    } else {
      rollingBlockSumNxNInternal(data, maxx, maxy, n);
    }
  }

  /**
   * Compute the block sum within a 2n+1 size block around each point. Only pixels with a full block
   * are processed. Pixels within border regions are unchanged.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   * @param n The block size
   */
  public void rollingBlockSumNxNInternal(float[] data, final int maxx, final int maxy,
      final int n) {
    final int blockSize = 2 * n + 1;
    if (maxx < blockSize || maxy < blockSize) {
      return;
    }

    final float[] newData = floatBuffer(data.length);

    // X-direction
    for (int y = 0; y < maxy; y++) {
      // Initialise the rolling sum
      float sum = 0;

      int endIndex = y * maxx;
      int x = 0;
      while (x < blockSize) {
        sum += data[endIndex];
        endIndex++;
        x++;
      }

      // Rolling sum over the X-direction
      int startIndex = y * maxx;
      int centreIndex = startIndex + n;

      newData[centreIndex] = sum;

      while (x < maxx) {
        centreIndex++;

        sum += data[endIndex] - data[startIndex];

        newData[centreIndex] = sum;

        x++;
        startIndex++;
        endIndex++;
      }
    }

    // Y-direction.
    // Only sweep over the interior
    for (int x = n; x < maxx - n; x++) {
      // Initialise the rolling sum
      float sum = 0;

      int endIndex = x;
      int y = 0;
      while (y < blockSize) {
        sum += newData[endIndex];
        endIndex += maxx;
        y++;
      }

      // Rolling sum over the Y-direction
      int startIndex = x;
      int centreIndex = startIndex + n * maxx;

      data[centreIndex] = sum;

      while (y < maxy) {
        centreIndex += maxx;

        sum += newData[endIndex] - newData[startIndex];

        data[centreIndex] = sum;

        y++;
        startIndex += maxx;
        endIndex += maxx;
      }
    }
  }

  /**
   * Compute the block sum within a 3x3 size block around each point. Only pixels with a full block
   * are processed. Pixels within border regions are unchanged.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   */
  public void rollingBlockSum3x3Internal(float[] data, final int maxx, final int maxy) {
    if (maxx < 3 || maxy < 3) {
      return;
    }

    final float[] newData = floatBuffer(data.length);

    // X-direction
    for (int y = 0; y < maxy; y++) {
      // Initialise the rolling sum
      int startIndex = y * maxx;
      int centreIndex = startIndex + 1;
      int endIndex = centreIndex + 1;
      float sum = data[startIndex] + data[centreIndex] + data[endIndex];

      // Rolling sum over the X-direction
      newData[centreIndex++] = sum;

      for (int x = 0; x < maxx - 3; x++) {
        sum += data[++endIndex] - data[startIndex++];
        newData[centreIndex++] = sum;
      }
    }

    // Y-direction.
    // Only sweep over the interior
    for (int x = 1; x < maxx - 1; x++) {
      // Initialise the rolling sum
      int startIndex = x;
      int centreIndex = startIndex + maxx;
      int endIndex = centreIndex + maxx;
      float sum = newData[startIndex] + newData[centreIndex] + newData[endIndex];

      // Rolling sum over the Y-direction
      data[centreIndex] = sum;

      for (int y = 0; y < maxy - 3; y++) {
        centreIndex += maxx;
        endIndex += maxx;
        sum += newData[endIndex] - newData[startIndex];
        data[centreIndex] = sum;
        startIndex += maxx;
      }
    }
  }

  /**
   * Compute the block sum within a 2n+1 size block around each point. Only pixels with a full block
   * are processed. Pixels within border regions are unchanged.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   * @param n The block size
   */
  @Deprecated
  public void rollingBlockSumNxNInternalTransposed(float[] data, final int maxx, final int maxy,
      final int n) {
    final int blockSize = 2 * n + 1;
    if (maxx < blockSize || maxy < blockSize) {
      return;
    }

    final float[] newData = floatBuffer(data.length);

    // X-direction
    for (int y = 0; y < maxy; y++) {
      // Initialise the rolling sum
      float sum = 0;

      int endIndex = y * maxx;
      int x = 0;
      while (x < blockSize) {
        sum += data[endIndex];
        endIndex++;
        x++;
      }

      // Rolling sum over the X-direction
      int startIndex = y * maxx;
      int newCentreIndex = y + n * maxy;

      newData[newCentreIndex] = sum;

      while (x < maxx) {
        newCentreIndex += maxy;

        sum += data[endIndex] - data[startIndex];

        newData[newCentreIndex] = sum;

        x++;
        startIndex++;
        endIndex++;
      }
    }

    // Y-direction.
    // Only sweep over the interior
    for (int x = n; x < maxx - n; x++) {
      // Initialise the rolling sum
      float sum = 0;

      int endIndex = x * maxy;
      int y = 0;
      while (y < blockSize) {
        sum += newData[endIndex];
        endIndex++;
        y++;
      }

      // Rolling sum over the Y-direction
      int startIndex = x * maxy;
      int centreIndex = x + n * maxx;

      data[centreIndex] = sum;

      while (y < maxy) {
        centreIndex += maxx;

        sum += newData[endIndex] - newData[startIndex];

        data[centreIndex] = sum;

        y++;
        startIndex++;
        endIndex++;
      }
    }
  }

  /**
   * Compute the block sum within a 2n+1 size block around each point. Only pixels with a full block
   * are processed. Pixels within border regions are unchanged.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   * @param n The block size
   */
  public void stripedBlockSumInternal(float[] data, final int maxx, final int maxy, final int n) {
    if (n == 1) {
      stripedBlockSum3x3Internal(data, maxx, maxy);
    } else {
      stripedBlockSumNxNInternal(data, maxx, maxy, n);
    }
  }

  /**
   * Compute the block sum within a 2n+1 size block around each point. Only pixels with a full block
   * are processed. Pixels within border regions are unchanged.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   * @param n The block size
   */
  public void stripedBlockSumNxNInternal(float[] data, final int maxx, final int maxy,
      final int n) {
    final int blockSize = 2 * n + 1;
    if (maxx < blockSize || maxy < blockSize) {
      return;
    }

    final float[] newData = floatBuffer(data.length);

    // X-direction
    for (int y = 0; y < maxy; y++) {
      int index = y * maxx;
      for (int x = 0; x <= maxx - blockSize; x++, index++) {
        float sum = 0;
        for (int x2 = 0; x2 < blockSize; x2++) {
          sum += data[index + x2];
        }
        newData[(x + n) * maxy + y] = sum;
      }
    }

    // Y-direction.
    // Only sweep over the interior
    for (int x = n; x < maxx - n; x++) {
      int index = x * maxy;
      for (int y = 0; y <= maxy - blockSize; y++, index++) {
        float sum = 0;
        for (int y2 = 0; y2 < blockSize; y2++) {
          sum += newData[index + y2];
        }
        data[x + (y + n) * maxx] = sum;
      }
    }
  }

  /**
   * Compute the block sum within a 3x3 size block around each point. Only pixels with a full block
   * are processed. Pixels within border regions are unchanged.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   */
  public void stripedBlockSum3x3Internal(float[] data, final int maxx, final int maxy) {
    if (maxx < 3 || maxy < 3) {
      return;
    }

    final float[] newData = floatBuffer(data.length);

    // X-direction
    for (int y = 0; y < maxy; y++) {
      int index = y * maxx;
      int index2 = maxy + y;
      for (int x = 0; x <= maxx - 3; x++, index++) {
        final float sum = data[index] + data[index + 1] + data[index + 2];
        newData[index2] = sum;
        index2 += maxy;
      }
    }

    // Y-direction.
    // Only sweep over the interior
    for (int x = 1; x < maxx - 1; x++) {
      int index = x * maxy;
      int index2 = x + maxx;
      for (int y = 0; y <= maxy - 3; y++, index++) {
        final float sum = newData[index] + newData[index + 1] + newData[index + 2];
        data[index2] = sum;
        index2 += maxx;
      }
    }
  }

  /**
   * Compute the block sum within a 2n+1 size block around each point. Only pixels with a full block
   * are processed. Pixels within border regions are unchanged.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   * @param n The block size
   */
  public void blockSumInternal(float[] data, final int maxx, final int maxy, final int n) {
    if (n == 1) {
      blockSum3x3Internal(data, maxx, maxy);
    } else {
      blockSumNxNInternal(data, maxx, maxy, n);
    }
  }

  /**
   * Compute the block sum within a 2n+1 size block around each point. Only pixels with a full block
   * are processed. Pixels within border regions are unchanged.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   * @param n The block size
   */
  public void blockSumNxNInternal(float[] data, final int maxx, final int maxy, final int n) {
    final int blockSize = 2 * n + 1;
    if (maxx < blockSize || maxy < blockSize) {
      return;
    }

    final float[] newData = floatBuffer(data.length);

    final int[] offset = new int[blockSize * blockSize - 1];
    for (int y = -n, d = 0; y <= n; y++) {
      for (int x = -n; x <= n; x++) {
        if (x != 0 || y != 0) {
          offset[d] = maxx * y + x;
          d++;
        }
      }
    }

    for (int y = n; y < maxy - n; y++) {
      int index = y * maxx + n;
      for (int x = n; x < maxx - n; x++, index++) {
        float sum = data[index];

        // Sweep neighbourhood -
        // No check for boundaries as this should be an internal sweep.
        for (final int shift : offset) {
          sum += data[index + shift];
        }

        newData[index] = sum;
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
   * Compute the block sum within a 3x3 size block around each point. Only pixels with a full block
   * are processed. Pixels within border regions are unchanged.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   */
  public void blockSum3x3Internal(float[] data, final int maxx, final int maxy) {
    final float[] newData = floatBuffer(data.length);

    for (int y = 1; y < maxy - 1; y++) {
      int index0 = (y - 1) * maxx + 1;
      int index1 = y * maxx + 1;
      int index2 = (y + 1) * maxx + 1;
      for (int x = 1; x < maxx - 1; x++) {
        final float sum = data[index0 - 1] + data[index0] + data[index0 + 1] + data[index1 - 1]
            + data[index1] + data[index1 + 1] + data[index2 - 1] + data[index2] + data[index2 + 1];
        newData[index1] = sum;
        index0++;
        index1++;
        index2++;
      }
    }

    // Copy back
    for (int y = 1; y < maxy - 1; y++) {
      int index = y * maxx + 1;
      for (int x = 1; x < maxx - 1; x++, index++) {
        data[index] = newData[index];
      }
    }
  }

  private float[] floatBuffer(int size) {
    if (floatDataBuffer == null || floatDataBuffer.length < size) {
      floatDataBuffer = new float[size];
    }
    return floatDataBuffer;
  }

  /**
   * Compute the block sum within a 2n+1 size block around each point.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   * @param n The block size
   */
  public void rollingBlockSum(float[] data, final int maxx, final int maxy, final int n) {
    if (n == 1) {
      rollingBlockSum3x3(data, maxx, maxy);
    } else {
      rollingBlockSumNxN(data, maxx, maxy, n);
    }
  }

  /**
   * Compute the block sum within a 2n+1 size block around each point.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   * @param n The block size
   */
  public void rollingBlockSumNxN(float[] data, final int maxx, final int maxy, final int n) {
    final float[] newData = floatBuffer(data.length);

    // X-direction
    int width = maxx;
    int height = maxy;
    float[] inData = data;
    float[] outData = newData;

    float[] row = floatRowBuffer(width + 2 * n);
    for (int y = 0; y < height; y++) {
      extractRow(inData, y, width, n, row);

      // Initialise rolling sum
      float sum = (n + 1) * row[0];
      int endIndex = n + 1;
      for (int i = 0; i < n; i++) {
        sum += row[endIndex++];
      }

      int centreIndex = y;
      outData[centreIndex] = sum;

      // Rolling sum over the X-direction
      for (int x = 0; x < width - 1; x++) {
        sum += row[endIndex++] - row[x];
        centreIndex += height;
        outData[centreIndex] = sum;
      }
    }

    // Y-direction.
    width = maxy;
    height = maxx;
    inData = newData;
    outData = data;

    row = floatRowBuffer(width + 2 * n);
    for (int y = 0; y < height; y++) {
      extractRow(inData, y, width, n, row);

      // Initialise rolling sum
      float sum = (n + 1) * row[0];
      int endIndex = n + 1;
      for (int i = 0; i < n; i++) {
        sum += row[endIndex++];
      }

      int centreIndex = y;
      outData[centreIndex] = sum;

      // Rolling sum over the X-direction
      for (int x = 0; x < width - 1; x++) {
        sum += row[endIndex++] - row[x];
        centreIndex += height;
        outData[centreIndex] = sum;
      }
    }
  }

  /**
   * Compute the block sum within a 3x3 size block around each point.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   */
  public void rollingBlockSum3x3(float[] data, final int maxx, final int maxy) {
    final float[] newData = floatBuffer(data.length);

    // X-direction
    int width = maxx;
    int height = maxy;
    float[] inData = data;
    float[] outData = newData;

    float[] row = floatRowBuffer(width + 2);
    for (int y = 0; y < height; y++) {
      extractRow1(inData, y, width, row);

      // Initialise rolling sum
      float sum = 2 * row[0] + row[2];
      int endIndex = 3;

      int centreIndex = y;
      outData[centreIndex] = sum;

      // Rolling sum over the X-direction
      for (int x = 0; x < width - 1; x++) {
        sum += row[endIndex++] - row[x];
        centreIndex += height;
        outData[centreIndex] = sum;
      }
    }

    // Y-direction.
    width = maxy;
    height = maxx;
    inData = newData;
    outData = data;

    row = floatRowBuffer(width + 2);
    for (int y = 0; y < height; y++) {
      extractRow1(inData, y, width, row);

      // Initialise rolling sum
      float sum = 2 * row[0] + row[2];
      int endIndex = 3;

      int centreIndex = y;
      outData[centreIndex] = sum;

      // Rolling sum over the X-direction
      for (int x = 0; x < width - 1; x++) {
        sum += row[endIndex++] - row[x];
        centreIndex += height;
        outData[centreIndex] = sum;
      }
    }
  }

  /**
   * Compute the block sum within a 2n+1 size block around each point.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   * @param n The block size
   */
  public void stripedBlockSum(float[] data, final int maxx, final int maxy, final int n) {
    if (n == 1) {
      stripedBlockSum3x3(data, maxx, maxy);
    } else {
      stripedBlockSumNxN(data, maxx, maxy, n);
    }
  }

  /**
   * Compute the block sum within a 2w+1 size block around each point.
   *
   * <p>Uses a [[w*w, w, ..., w, w*w], [w, 1, ..., 1, w], [w*w, w, ..., w, w*w]] convolution kernel.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   * @param weight The block size
   */
  public void stripedBlockSum(float[] data, final int maxx, final int maxy, final float weight) {
    if (weight <= 1) {
      stripedBlockSum3x3(data, maxx, maxy, weight);
    } else {
      stripedBlockSumNxN(data, maxx, maxy, weight);
    }
  }

  /**
   * Compute the block sum within a 2n+1 size block around each point.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   * @param n The block size
   */
  public void stripedBlockSumNxN(float[] data, final int maxx, final int maxy, final int n) {
    final int blockSize = 2 * n + 1;

    final float[] newData = floatBuffer(data.length);

    // X-direction
    int width = maxx;
    int height = maxy;
    float[] inData = data;
    float[] outData = newData;

    float[] row = floatRowBuffer(width + 2 * n);
    for (int y = 0; y < height; y++) {
      extractRow(inData, y, width, n, row);

      int centreIndex = y;
      for (int x = 0; x < width; x++) {
        // Sum strips
        float sum = 0;

        for (int j = 0; j < blockSize; j++) {
          sum += row[x + j];
        }

        // Store result in transpose
        outData[centreIndex] = sum;
        centreIndex += height;
      }
    }

    // Y-direction.
    width = maxy;
    height = maxx;
    inData = newData;
    outData = data;

    row = floatRowBuffer(width + 2 * n);
    for (int y = 0; y < height; y++) {
      // Extract row (pad ends)
      extractRow(inData, y, width, n, row);

      int centreIndex = y;
      for (int x = 0; x < width; x++) {
        // Sum strips
        float sum = 0;

        for (int j = 0; j < blockSize; j++) {
          sum += row[x + j];
        }

        // Store result in transpose
        outData[centreIndex] = sum;
        centreIndex += height;
      }
    }
  }

  /**
   * Compute the block sum within a 2w+1 size block around each point.
   *
   * <p>Uses a [[w*w, w, ..., w, w*w], [w, 1, ..., 1, w], [w*w, w, ..., w, w*w]] convolution kernel.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   * @param weight The block size
   */
  public void stripedBlockSumNxN(float[] data, final int maxx, final int maxy, final float weight) {
    final int n = (int) weight;
    final int n1 = (n == weight) ? n : n + 1;

    if (n == n1) {
      // There is no edge
      stripedBlockSum(data, maxx, maxy, n);
      return;
    }

    final int blockSize = 2 * n1;

    final float[] newData = floatBuffer(data.length);

    final float w1 = weight - n;

    // NOTE:
    // To increase speed when sweeping the arrays and allow for reusing code:
    // newData is XY ordinal => x * maxy + y
    // data is YX ordinal => y * maxx + x

    // X-direction
    int width = maxx;
    int height = maxy;
    float[] inData = data;
    float[] outData = newData;

    float[] row = floatRowBuffer(width + 2 * n1);
    for (int y = 0; y < height; y++) {
      extractRow(inData, y, width, n1, row);

      int centreIndex = y;
      for (int x = 0; x < width; x++) {
        // Sum strips
        float sum = row[x] * w1;
        for (int j = 1; j < blockSize; j++) {
          sum += row[x + j];
        }
        sum += row[x + blockSize] * w1;

        // Store result in transpose
        outData[centreIndex] = sum;
        centreIndex += height;
      }
    }

    // Y-direction.
    width = maxy;
    height = maxx;
    inData = newData;
    outData = data;

    row = floatRowBuffer(width + 2 * n1);
    for (int y = 0; y < height; y++) {
      // Extract row (pad ends)
      extractRow(inData, y, width, n1, row);

      int centreIndex = y;
      for (int x = 0; x < width; x++) {
        // Sum strips
        float sum = row[x] * w1;
        for (int j = 1; j < blockSize; j++) {
          sum += row[x + j];
        }
        sum += row[x + blockSize] * w1;

        // Store result in transpose
        outData[centreIndex] = sum;
        centreIndex += height;
      }
    }
  }

  /**
   * Compute the block sum within a 3x3 size block around each point.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   */
  public void stripedBlockSum3x3(float[] data, final int maxx, final int maxy) {
    final float[] newData = floatBuffer(data.length);

    // X-direction
    int width = maxx;
    int height = maxy;
    float[] inData = data;
    float[] outData = newData;

    float[] row = floatRowBuffer(width + 2);
    for (int y = 0; y < height; y++) {
      extractRow1(inData, y, width, row);

      int centreIndex = y;
      for (int x = 0; x < width; x++) {
        // Sum strips
        final float sum = row[x] + row[x + 1] + row[x + 2];

        // Store result in transpose
        outData[centreIndex] = sum;
        centreIndex += height;
      }
    }

    // Y-direction.
    width = maxy;
    height = maxx;
    inData = newData;
    outData = data;

    row = floatRowBuffer(width + 2);
    for (int y = 0; y < height; y++) {
      // Extract row (pad ends)
      extractRow1(inData, y, width, row);

      int centreIndex = y;
      for (int x = 0; x < width; x++) {
        // Sum strips
        final float sum = row[x] + row[x + 1] + row[x + 2];

        // Store result in transpose
        outData[centreIndex] = sum;
        centreIndex += height;
      }
    }
  }

  /**
   * Compute the block sum within a 3x3 size block around each point.
   *
   * <p>Uses a [[w*w, w, w*w], [w, 1, w], [w*w, w, w*w]] convolution kernel.
   *
   * <p>Note: the input data is destructively modified.
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   * @param weight The weight
   */
  public void stripedBlockSum3x3(float[] data, final int maxx, final int maxy, final float weight) {
    final float[] newData = floatBuffer(data.length);

    // NOTE:
    // To increase speed when sweeping the arrays and allow for reusing code:
    // newData is XY ordinal => x * maxy + y
    // data is YX ordinal => y * maxx + x

    // X-direction
    int width = maxx;
    int height = maxy;
    float[] inData = data;
    float[] outData = newData;

    float[] row = floatRowBuffer(width + 2);
    for (int y = 0; y < height; y++) {
      extractRow1(inData, y, width, row);

      int centreIndex = y;
      for (int x = 0; x < width; x++) {
        // Sum strips
        // Store result in transpose
        outData[centreIndex] = weight * (row[x] + row[x + 2]) + row[x + 1];
        centreIndex += height;
      }
    }

    // Y-direction.
    width = maxy;
    height = maxx;
    inData = newData;
    outData = data;

    row = floatRowBuffer(width + 2);
    for (int y = 0; y < height; y++) {
      // Extract row (pad ends)
      extractRow1(inData, y, width, row);

      int centreIndex = y;
      for (int x = 0; x < width; x++) {
        // Sum strips
        // Store result in transpose
        outData[centreIndex] = weight * (row[x] + row[x + 2]) + row[x + 1];
        centreIndex += height;
      }
    }
  }

  private float[] floatRowBuffer(int size) {
    if (floatRowBuffer == null || floatRowBuffer.length < size) {
      floatRowBuffer = new float[size];
    }
    return floatRowBuffer;
  }

  private static void extractRow(float[] inData, int y, int width, final int n, float[] row) {
    int index = y * width;

    // Pad ends
    for (int i = 0; i < n; i++) {
      row[i] = inData[index];
      row[i + n + width] = inData[index + width - 1];
    }

    // Fill in data
    for (int x = 0, i = n; x < width; x++) {
      row[i++] = inData[index++];
    }
  }

  private static void extractRow1(float[] inData, int y, int width, float[] row) {
    int index = y * width;

    // Pad ends
    row[0] = inData[index];
    row[1 + width] = inData[index + width - 1];

    // Fill in data
    for (int x = 0, i = 1; x < width; x++) {
      row[i++] = inData[index++];
    }
  }

  /**
   * Compute the block sum within a 2n+1 size block around each point.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   * @param n The block size
   */
  public void blockSum(float[] data, final int maxx, final int maxy, final int n) {
    if (n == 1) {
      blockSum3x3(data, maxx, maxy);
    } else {
      blockSumNxN(data, maxx, maxy, n);
    }
  }

  /**
   * Compute the block sum within a 2n+1 size block around each point.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   * @param n The block size
   */
  public void blockSumNxN(float[] data, final int maxx, final int maxy, final int n) {
    final float[] newData = floatBuffer(data.length);

    // Boundary control
    final int xwidth = FastMath.min(n, maxx - 1);
    final int ywidth = FastMath.min(n, maxy - 1);
    final int xlimit = maxx - xwidth;
    final int ylimit = maxy - ywidth;

    final int[] offset = new int[(2 * xwidth + 1) * (2 * ywidth + 1) - 1];
    final int[] xoffset = new int[offset.length];
    final int[] yoffset = new int[offset.length];
    int offsetIndex = 0;
    for (int y = -ywidth; y <= ywidth; y++) {
      for (int x = -xwidth; x <= xwidth; x++) {
        if (x != 0 || y != 0) {
          offset[offsetIndex] = maxx * y + x;
          xoffset[offsetIndex] = x;
          yoffset[offsetIndex] = y;
          offsetIndex++;
        }
      }
    }

    int index = 0;
    for (int y = 0; y < maxy; y++) {
      for (int x = 0; x < maxx; x++, index++) {
        float sum = data[index];

        // Flag to indicate this pixels has a complete (2n+1) neighbourhood
        final boolean isInnerXy = (y >= ywidth && y < ylimit) && (x >= xwidth && x < xlimit);

        // Sweep neighbourhood
        if (isInnerXy) {
          for (offsetIndex = offset.length; offsetIndex-- > 0;) {
            sum += data[index + offset[offsetIndex]];
          }
        } else {
          for (offsetIndex = offset.length; offsetIndex-- > 0;) {
            // Get the pixel with boundary checking
            int yy = y + yoffset[offsetIndex];
            int xx = x + xoffset[offsetIndex];
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
            sum += data[xx + yy * maxx];
          }
        }
        newData[index] = sum;
      }
    }

    // Copy back
    System.arraycopy(newData, 0, data, 0, data.length);
  }

  /**
   * Compute the block sum within a 3x3 size block around each point.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   */
  public void blockSum3x3(float[] data, final int maxx, final int maxy) {
    final int n = 1;
    final float[] newData = floatBuffer(data.length);

    // Boundary control
    final int xwidth = FastMath.min(n, maxx - 1);
    final int ywidth = FastMath.min(n, maxy - 1);
    final int xlimit = maxx - xwidth;
    final int ylimit = maxy - ywidth;

    final int[] offset = new int[(2 * xwidth + 1) * (2 * ywidth + 1) - 1];
    final int[] xoffset = new int[offset.length];
    final int[] yoffset = new int[offset.length];
    int offsetIndex = 0;
    for (int y = -ywidth; y <= ywidth; y++) {
      for (int x = -xwidth; x <= xwidth; x++) {
        if (x != 0 || y != 0) {
          offset[offsetIndex] = maxx * y + x;
          xoffset[offsetIndex] = x;
          yoffset[offsetIndex] = y;
          offsetIndex++;
        }
      }
    }

    for (int y = 0; y < maxy; y++) {
      int index0 = (y - 1) * maxx;
      int index1 = y * maxx;
      int index2 = (y + 1) * maxx;
      for (int x = 0; x < maxx; x++) {
        // Flag to indicate this pixels has a complete (2n+1) neighbourhood
        final boolean isInnerXy = (y >= ywidth && y < ylimit) && (x >= xwidth && x < xlimit);

        // Sweep neighbourhood
        if (isInnerXy) {
          final float sum =
              data[index0 - 1] + data[index0] + data[index0 + 1] + data[index1 - 1] + data[index1]
                  + data[index1 + 1] + data[index2 - 1] + data[index2] + data[index2 + 1];
          newData[index1] = sum;
        } else {
          float sum = data[index1];

          for (offsetIndex = offset.length; offsetIndex-- > 0;) {
            // Get the pixel with boundary checking
            int yy = y + yoffset[offsetIndex];
            int xx = x + xoffset[offsetIndex];
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
            sum += data[xx + yy * maxx];
          }
          newData[index1] = sum;
        }
        index0++;
        index1++;
        index2++;
      }
    }

    // Copy back
    System.arraycopy(newData, 0, data, 0, data.length);
  }

  // ----------------------------------------------------
  // XXX
  // NOTE:
  // The following code is copied directly from above.
  // All 'float' have been replaced with 'int'.
  // ----------------------------------------------------
  // @CHECKSTYLE.OFF: OverloadMethodsDeclarationOrder
  private int[] intDataBuffer;
  private int[] intRowBuffer;

  /**
   * Compute the block sum within a 2n+1 size block around each point. Only pixels with a full block
   * are processed. Pixels within border regions are unchanged.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   * @param n The block size
   */
  public void rollingBlockSumInternal(int[] data, final int maxx, final int maxy, final int n) {
    // Note: Speed tests show that this method is only marginally faster than
    // rollingBlockSumNxNInternal.
    // Sometimes it is slower. The intricacies of the java optimiser escape me.
    if (n == 1) {
      rollingBlockSum3x3Internal(data, maxx, maxy);
    } else {
      rollingBlockSumNxNInternal(data, maxx, maxy, n);
    }
  }

  /**
   * Compute the block sum within a 2n+1 size block around each point. Only pixels with a full block
   * are processed. Pixels within border regions are unchanged.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   * @param n The block size
   */
  public void rollingBlockSumNxNInternal(int[] data, final int maxx, final int maxy, final int n) {
    final int blockSize = 2 * n + 1;
    if (maxx < blockSize || maxy < blockSize) {
      return;
    }

    final int[] newData = intBuffer(data.length);

    // X-direction
    for (int y = 0; y < maxy; y++) {
      // Initialise the rolling sum
      int sum = 0;

      int endIndex = y * maxx;
      int x = 0;
      while (x < blockSize) {
        sum += data[endIndex];
        endIndex++;
        x++;
      }

      // Rolling sum over the X-direction
      int startIndex = y * maxx;
      int centreIndex = startIndex + n;

      newData[centreIndex] = sum;

      while (x < maxx) {
        centreIndex++;

        sum += data[endIndex] - data[startIndex];

        newData[centreIndex] = sum;

        x++;
        startIndex++;
        endIndex++;
      }
    }

    // Y-direction.
    // Only sweep over the interior
    for (int x = n; x < maxx - n; x++) {
      // Initialise the rolling sum
      int sum = 0;

      int endIndex = x;
      int y = 0;
      while (y < blockSize) {
        sum += newData[endIndex];
        endIndex += maxx;
        y++;
      }

      // Rolling sum over the Y-direction
      int startIndex = x;
      int centreIndex = startIndex + n * maxx;

      data[centreIndex] = sum;

      while (y < maxy) {
        centreIndex += maxx;

        sum += newData[endIndex] - newData[startIndex];

        data[centreIndex] = sum;

        y++;
        startIndex += maxx;
        endIndex += maxx;
      }
    }
  }

  /**
   * Compute the block sum within a 3x3 size block around each point. Only pixels with a full block
   * are processed. Pixels within border regions are unchanged.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   */
  public void rollingBlockSum3x3Internal(int[] data, final int maxx, final int maxy) {
    if (maxx < 3 || maxy < 3) {
      return;
    }

    final int[] newData = intBuffer(data.length);

    // X-direction
    for (int y = 0; y < maxy; y++) {
      // Initialise the rolling sum
      int startIndex = y * maxx;
      int centreIndex = startIndex + 1;
      int endIndex = centreIndex + 1;
      int sum = data[startIndex] + data[centreIndex] + data[endIndex];

      // Rolling sum over the X-direction
      newData[centreIndex++] = sum;

      for (int x = 0; x < maxx - 3; x++) {
        sum += data[++endIndex] - data[startIndex++];
        newData[centreIndex++] = sum;
      }
    }

    // Y-direction.
    // Only sweep over the interior
    for (int x = 1; x < maxx - 1; x++) {
      // Initialise the rolling sum
      int startIndex = x;
      int centreIndex = startIndex + maxx;
      int endIndex = centreIndex + maxx;
      int sum = newData[startIndex] + newData[centreIndex] + newData[endIndex];

      // Rolling sum over the Y-direction
      data[centreIndex] = sum;

      for (int y = 0; y < maxy - 3; y++) {
        centreIndex += maxx;
        endIndex += maxx;
        sum += newData[endIndex] - newData[startIndex];
        data[centreIndex] = sum;
        startIndex += maxx;
      }
    }
  }

  /**
   * Compute the block sum within a 2n+1 size block around each point. Only pixels with a full block
   * are processed. Pixels within border regions are unchanged.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   * @param n The block size
   */
  @Deprecated
  public void rollingBlockSumNxNInternalTransposed(int[] data, final int maxx, final int maxy,
      final int n) {
    final int blockSize = 2 * n + 1;
    if (maxx < blockSize || maxy < blockSize) {
      return;
    }

    final int[] newData = intBuffer(data.length);

    // X-direction
    for (int y = 0; y < maxy; y++) {
      // Initialise the rolling sum
      int sum = 0;

      int endIndex = y * maxx;
      int x = 0;
      while (x < blockSize) {
        sum += data[endIndex];
        endIndex++;
        x++;
      }

      // Rolling sum over the X-direction
      int startIndex = y * maxx;
      int newCentreIndex = y + n * maxy;

      newData[newCentreIndex] = sum;

      while (x < maxx) {
        newCentreIndex += maxy;

        sum += data[endIndex] - data[startIndex];

        newData[newCentreIndex] = sum;

        x++;
        startIndex++;
        endIndex++;
      }
    }

    // Y-direction.
    // Only sweep over the interior
    for (int x = n; x < maxx - n; x++) {
      // Initialise the rolling sum
      int sum = 0;

      int endIndex = x * maxy;
      int y = 0;
      while (y < blockSize) {
        sum += newData[endIndex];
        endIndex++;
        y++;
      }

      // Rolling sum over the Y-direction
      int startIndex = x * maxy;
      int centreIndex = x + n * maxx;

      data[centreIndex] = sum;

      while (y < maxy) {
        centreIndex += maxx;

        sum += newData[endIndex] - newData[startIndex];

        data[centreIndex] = sum;

        y++;
        startIndex++;
        endIndex++;
      }
    }
  }

  /**
   * Compute the block sum within a 2n+1 size block around each point. Only pixels with a full block
   * are processed. Pixels within border regions are unchanged.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   * @param n The block size
   */
  public void stripedBlockSumInternal(int[] data, final int maxx, final int maxy, final int n) {
    if (n == 1) {
      stripedBlockSum3x3Internal(data, maxx, maxy);
    } else {
      stripedBlockSumNxNInternal(data, maxx, maxy, n);
    }
  }

  /**
   * Compute the block sum within a 2n+1 size block around each point. Only pixels with a full block
   * are processed. Pixels within border regions are unchanged.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   * @param n The block size
   */
  public void stripedBlockSumNxNInternal(int[] data, final int maxx, final int maxy, final int n) {
    final int blockSize = 2 * n + 1;
    if (maxx < blockSize || maxy < blockSize) {
      return;
    }

    final int[] newData = intBuffer(data.length);

    // X-direction
    for (int y = 0; y < maxy; y++) {
      int index = y * maxx;
      for (int x = 0; x <= maxx - blockSize; x++, index++) {
        int sum = 0;
        for (int x2 = 0; x2 < blockSize; x2++) {
          sum += data[index + x2];
        }
        newData[(x + n) * maxy + y] = sum;
      }
    }

    // Y-direction.
    // Only sweep over the interior
    for (int x = n; x < maxx - n; x++) {
      int index = x * maxy;
      for (int y = 0; y <= maxy - blockSize; y++, index++) {
        int sum = 0;
        for (int y2 = 0; y2 < blockSize; y2++) {
          sum += newData[index + y2];
        }
        data[x + (y + n) * maxx] = sum;
      }
    }
  }

  /**
   * Compute the block sum within a 3x3 size block around each point. Only pixels with a full block
   * are processed. Pixels within border regions are unchanged.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   */
  public void stripedBlockSum3x3Internal(int[] data, final int maxx, final int maxy) {
    if (maxx < 3 || maxy < 3) {
      return;
    }

    final int[] newData = intBuffer(data.length);

    // X-direction
    for (int y = 0; y < maxy; y++) {
      int index = y * maxx;
      int index2 = maxy + y;
      for (int x = 0; x <= maxx - 3; x++, index++) {
        final int sum = data[index] + data[index + 1] + data[index + 2];
        newData[index2] = sum;
        index2 += maxy;
      }
    }

    // Y-direction.
    // Only sweep over the interior
    for (int x = 1; x < maxx - 1; x++) {
      int index = x * maxy;
      int index2 = x + maxx;
      for (int y = 0; y <= maxy - 3; y++, index++) {
        final int sum = newData[index] + newData[index + 1] + newData[index + 2];
        data[index2] = sum;
        index2 += maxx;
      }
    }
  }

  /**
   * Compute the block sum within a 2n+1 size block around each point. Only pixels with a full block
   * are processed. Pixels within border regions are unchanged.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   * @param n The block size
   */
  public void blockSumInternal(int[] data, final int maxx, final int maxy, final int n) {
    if (n == 1) {
      blockSum3x3Internal(data, maxx, maxy);
    } else {
      blockSumNxNInternal(data, maxx, maxy, n);
    }
  }

  /**
   * Compute the block sum within a 2n+1 size block around each point. Only pixels with a full block
   * are processed. Pixels within border regions are unchanged.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   * @param n The block size
   */
  public void blockSumNxNInternal(int[] data, final int maxx, final int maxy, final int n) {
    final int blockSize = 2 * n + 1;
    if (maxx < blockSize || maxy < blockSize) {
      return;
    }

    final int[] newData = intBuffer(data.length);

    final int[] offset = new int[blockSize * blockSize - 1];
    for (int y = -n, d = 0; y <= n; y++) {
      for (int x = -n; x <= n; x++) {
        if (x != 0 || y != 0) {
          offset[d] = maxx * y + x;
          d++;
        }
      }
    }

    for (int y = n; y < maxy - n; y++) {
      int index = y * maxx + n;
      for (int x = n; x < maxx - n; x++, index++) {
        int sum = data[index];

        // Sweep neighbourhood -
        // No check for boundaries as this should be an internal sweep.
        for (final int shift : offset) {
          sum += data[index + shift];
        }

        newData[index] = sum;
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
   * Compute the block sum within a 3x3 size block around each point. Only pixels with a full block
   * are processed. Pixels within border regions are unchanged.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   */
  public void blockSum3x3Internal(int[] data, final int maxx, final int maxy) {
    final int[] newData = intBuffer(data.length);

    for (int y = 1; y < maxy - 1; y++) {
      int index0 = (y - 1) * maxx + 1;
      int index1 = y * maxx + 1;
      int index2 = (y + 1) * maxx + 1;
      for (int x = 1; x < maxx - 1; x++) {
        final int sum = data[index0 - 1] + data[index0] + data[index0 + 1] + data[index1 - 1]
            + data[index1] + data[index1 + 1] + data[index2 - 1] + data[index2] + data[index2 + 1];
        newData[index1] = sum;
        index0++;
        index1++;
        index2++;
      }
    }

    // Copy back
    for (int y = 1; y < maxy - 1; y++) {
      int index = y * maxx + 1;
      for (int x = 1; x < maxx - 1; x++, index++) {
        data[index] = newData[index];
      }
    }
  }

  private int[] intBuffer(int size) {
    if (intDataBuffer == null || intDataBuffer.length < size) {
      intDataBuffer = new int[size];
    }
    return intDataBuffer;
  }

  /**
   * Compute the block sum within a 2n+1 size block around each point.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   * @param n The block size
   */
  public void rollingBlockSum(int[] data, final int maxx, final int maxy, final int n) {
    if (n == 1) {
      rollingBlockSum3x3(data, maxx, maxy);
    } else {
      rollingBlockSumNxN(data, maxx, maxy, n);
    }
  }

  /**
   * Compute the block sum within a 2n+1 size block around each point.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   * @param n The block size
   */
  public void rollingBlockSumNxN(int[] data, final int maxx, final int maxy, final int n) {
    final int[] newData = intBuffer(data.length);

    // X-direction
    int width = maxx;
    int height = maxy;
    int[] inData = data;
    int[] outData = newData;

    int[] row = intRowBuffer(width + 2 * n);
    for (int y = 0; y < height; y++) {
      extractRow(inData, y, width, n, row);

      // Initialise rolling sum
      int sum = (n + 1) * row[0];
      int endIndex = n + 1;
      for (int i = 0; i < n; i++) {
        sum += row[endIndex++];
      }

      int centreIndex = y;
      outData[centreIndex] = sum;

      // Rolling sum over the X-direction
      for (int x = 0; x < width - 1; x++) {
        sum += row[endIndex++] - row[x];
        centreIndex += height;
        outData[centreIndex] = sum;
      }
    }

    // Y-direction.
    width = maxy;
    height = maxx;
    inData = newData;
    outData = data;

    row = intRowBuffer(width + 2 * n);
    for (int y = 0; y < height; y++) {
      extractRow(inData, y, width, n, row);

      // Initialise rolling sum
      int sum = (n + 1) * row[0];
      int endIndex = n + 1;
      for (int i = 0; i < n; i++) {
        sum += row[endIndex++];
      }

      int centreIndex = y;
      outData[centreIndex] = sum;

      // Rolling sum over the X-direction
      for (int x = 0; x < width - 1; x++) {
        sum += row[endIndex++] - row[x];
        centreIndex += height;
        outData[centreIndex] = sum;
      }
    }
  }

  /**
   * Compute the block sum within a 3x3 size block around each point.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   */
  public void rollingBlockSum3x3(int[] data, final int maxx, final int maxy) {
    final int[] newData = intBuffer(data.length);

    // X-direction
    int width = maxx;
    int height = maxy;
    int[] inData = data;
    int[] outData = newData;

    int[] row = intRowBuffer(width + 2);
    for (int y = 0; y < height; y++) {
      extractRow1(inData, y, width, row);

      // Initialise rolling sum
      int sum = 2 * row[0] + row[2];
      int endIndex = 3;

      int centreIndex = y;
      outData[centreIndex] = sum;

      // Rolling sum over the X-direction
      for (int x = 0; x < width - 1; x++) {
        sum += row[endIndex++] - row[x];
        centreIndex += height;
        outData[centreIndex] = sum;
      }
    }

    // Y-direction.
    width = maxy;
    height = maxx;
    inData = newData;
    outData = data;

    row = intRowBuffer(width + 2);
    for (int y = 0; y < height; y++) {
      extractRow1(inData, y, width, row);

      // Initialise rolling sum
      int sum = 2 * row[0] + row[2];
      int endIndex = 3;

      int centreIndex = y;
      outData[centreIndex] = sum;

      // Rolling sum over the X-direction
      for (int x = 0; x < width - 1; x++) {
        sum += row[endIndex++] - row[x];
        centreIndex += height;
        outData[centreIndex] = sum;
      }
    }
  }

  /**
   * Compute the block sum within a 2n+1 size block around each point.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   * @param n The block size
   */
  public void stripedBlockSum(int[] data, final int maxx, final int maxy, final int n) {
    if (n == 1) {
      stripedBlockSum3x3(data, maxx, maxy);
    } else {
      stripedBlockSumNxN(data, maxx, maxy, n);
    }
  }

  /**
   * Compute the block sum within a 2n+1 size block around each point.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   * @param n The block size
   */
  public void stripedBlockSumNxN(int[] data, final int maxx, final int maxy, final int n) {
    final int blockSize = 2 * n + 1;

    final int[] newData = intBuffer(data.length);

    // X-direction
    int width = maxx;
    int height = maxy;
    int[] inData = data;
    int[] outData = newData;

    int[] row = intRowBuffer(width + 2 * n);
    for (int y = 0; y < height; y++) {
      extractRow(inData, y, width, n, row);

      int centreIndex = y;
      for (int x = 0; x < width; x++) {
        // Sum strips
        int sum = 0;

        for (int j = 0; j < blockSize; j++) {
          sum += row[x + j];
        }

        // Store result in transpose
        outData[centreIndex] = sum;
        centreIndex += height;
      }
    }

    // Y-direction.
    width = maxy;
    height = maxx;
    inData = newData;
    outData = data;

    row = intRowBuffer(width + 2 * n);
    for (int y = 0; y < height; y++) {
      // Extract row (pad ends)
      extractRow(inData, y, width, n, row);

      int centreIndex = y;
      for (int x = 0; x < width; x++) {
        // Sum strips
        int sum = 0;

        for (int j = 0; j < blockSize; j++) {
          sum += row[x + j];
        }

        // Store result in transpose
        outData[centreIndex] = sum;
        centreIndex += height;
      }
    }
  }

  /**
   * Compute the block sum within a 3x3 size block around each point.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   */
  public void stripedBlockSum3x3(int[] data, final int maxx, final int maxy) {
    final int[] newData = intBuffer(data.length);

    // X-direction
    int width = maxx;
    int height = maxy;
    int[] inData = data;
    int[] outData = newData;

    int[] row = intRowBuffer(width + 2);
    for (int y = 0; y < height; y++) {
      extractRow1(inData, y, width, row);

      int centreIndex = y;
      for (int x = 0; x < width; x++) {
        // Sum strips
        final int sum = row[x] + row[x + 1] + row[x + 2];

        // Store result in transpose
        outData[centreIndex] = sum;
        centreIndex += height;
      }
    }

    // Y-direction.
    width = maxy;
    height = maxx;
    inData = newData;
    outData = data;

    row = intRowBuffer(width + 2);
    for (int y = 0; y < height; y++) {
      // Extract row (pad ends)
      extractRow1(inData, y, width, row);

      int centreIndex = y;
      for (int x = 0; x < width; x++) {
        // Sum strips
        final int sum = row[x] + row[x + 1] + row[x + 2];

        // Store result in transpose
        outData[centreIndex] = sum;
        centreIndex += height;
      }
    }
  }

  private int[] intRowBuffer(int size) {
    if (intRowBuffer == null || intRowBuffer.length < size) {
      intRowBuffer = new int[size];
    }
    return intRowBuffer;
  }

  private static void extractRow(int[] inData, int y, int width, final int n, int[] row) {
    int index = y * width;

    // Pad ends
    for (int i = 0; i < n; i++) {
      row[i] = inData[index];
      row[i + n + width] = inData[index + width - 1];
    }

    // Fill in data
    for (int x = 0, i = n; x < width; x++) {
      row[i++] = inData[index++];
    }
  }

  private static void extractRow1(int[] inData, int y, int width, int[] row) {
    int index = y * width;

    // Pad ends
    row[0] = inData[index];
    row[1 + width] = inData[index + width - 1];

    // Fill in data
    for (int x = 0, i = 1; x < width; x++) {
      row[i++] = inData[index++];
    }
  }

  /**
   * Compute the block sum within a 2n+1 size block around each point.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   * @param n The block size
   */
  public void blockSum(int[] data, final int maxx, final int maxy, final int n) {
    if (n == 1) {
      blockSum3x3(data, maxx, maxy);
    } else {
      blockSumNxN(data, maxx, maxy, n);
    }
  }

  /**
   * Compute the block sum within a 2n+1 size block around each point.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   * @param n The block size
   */
  public void blockSumNxN(int[] data, final int maxx, final int maxy, final int n) {
    final int[] newData = intBuffer(data.length);

    // Boundary control
    final int xwidth = FastMath.min(n, maxx - 1);
    final int ywidth = FastMath.min(n, maxy - 1);
    final int xlimit = maxx - xwidth;
    final int ylimit = maxy - ywidth;

    final int[] offset = new int[(2 * xwidth + 1) * (2 * ywidth + 1) - 1];
    final int[] xoffset = new int[offset.length];
    final int[] yoffset = new int[offset.length];
    int d = 0;
    for (int y = -ywidth; y <= ywidth; y++) {
      for (int x = -xwidth; x <= xwidth; x++) {
        if (x != 0 || y != 0) {
          offset[d] = maxx * y + x;
          xoffset[d] = x;
          yoffset[d] = y;
          d++;
        }
      }
    }

    int index = 0;
    for (int y = 0; y < maxy; y++) {
      for (int x = 0; x < maxx; x++, index++) {
        int sum = data[index];

        // Flag to indicate this pixels has a complete (2n+1) neighbourhood
        final boolean isInnerXy = (y >= ywidth && y < ylimit) && (x >= xwidth && x < xlimit);

        // Sweep neighbourhood
        if (isInnerXy) {
          for (d = offset.length; d-- > 0;) {
            sum += data[index + offset[d]];
          }
        } else {
          for (d = offset.length; d-- > 0;) {
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
            sum += data[xx + yy * maxx];
          }
        }
        newData[index] = sum;
      }
    }

    // Copy back
    System.arraycopy(newData, 0, data, 0, data.length);
  }

  /**
   * Compute the block sum within a 3x3 size block around each point.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   */
  public void blockSum3x3(int[] data, final int maxx, final int maxy) {
    final int n = 1;
    final int[] newData = intBuffer(data.length);

    // Boundary control
    final int xwidth = FastMath.min(n, maxx - 1);
    final int ywidth = FastMath.min(n, maxy - 1);
    final int xlimit = maxx - xwidth;
    final int ylimit = maxy - ywidth;

    final int[] offset = new int[(2 * xwidth + 1) * (2 * ywidth + 1) - 1];
    final int[] xoffset = new int[offset.length];
    final int[] yoffset = new int[offset.length];
    int d = 0;
    for (int y = -ywidth; y <= ywidth; y++) {
      for (int x = -xwidth; x <= xwidth; x++) {
        if (x != 0 || y != 0) {
          offset[d] = maxx * y + x;
          xoffset[d] = x;
          yoffset[d] = y;
          d++;
        }
      }
    }

    for (int y = 0; y < maxy; y++) {
      int index0 = (y - 1) * maxx;
      int index1 = y * maxx;
      int index2 = (y + 1) * maxx;
      for (int x = 0; x < maxx; x++) {
        // Flag to indicate this pixels has a complete (2n+1) neighbourhood
        final boolean isInnerXy = (y >= ywidth && y < ylimit) && (x >= xwidth && x < xlimit);

        // Sweep neighbourhood
        if (isInnerXy) {
          final int sum =
              data[index0 - 1] + data[index0] + data[index0 + 1] + data[index1 - 1] + data[index1]
                  + data[index1 + 1] + data[index2 - 1] + data[index2] + data[index2 + 1];
          newData[index1] = sum;
        } else {
          int sum = data[index1];

          for (d = offset.length; d-- > 0;) {
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
            sum += data[xx + yy * maxx];
          }
          newData[index1] = sum;
        }
        index0++;
        index1++;
        index2++;
      }
    }

    // Copy back
    System.arraycopy(newData, 0, data, 0, data.length);
  }
}
