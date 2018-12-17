/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
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

/**
 * Computes the area average for each point within the array.
 *
 * <p>The algorithm computes the average in the specified area. If the width is an integer the area
 * is defined using an entire block. If the width is non-integer it is rounded down and up to the
 * adjacent integers. The sum of each block is computed. The larger block is subtracted from the
 * smaller block to give the edge sum. The average is then computed using the sum of the smaller
 * block and a proportion of the edge sum.
 *
 *
 * <p>Note: The algorithm allocates a disproportionate weighting to the corner pixels. Each edge
 * pixel should get a weighting of w and corner pixels w*w. Using simple weighting of the inner and
 * outer blocks results in all the outer pixels receiving the same weight.
 *
 * <p>Note: Due to lack of small dimension checking the routines will fail if maxx or maxy are less
 * than 2. All routines are OK for 3x3 images and larger.
 */
public class AreaAverageFilter extends BaseWeightedFilter {
  // Use duplicate filters to support efficient caching of the weights
  private final BlockSumFilter sumFilter1;
  private final BlockSumFilter sumFilter2;
  private final BlockMeanFilter blockMeanFilter1;
  private final BlockMeanFilter blockMeanFilter2;

  private boolean simpleInterpolation;

  /**
   * Instantiates a new area average filter.
   */
  public AreaAverageFilter() {
    sumFilter1 = new BlockSumFilter();
    sumFilter2 = new BlockSumFilter();
    blockMeanFilter1 = new BlockMeanFilter();
    blockMeanFilter2 = new BlockMeanFilter();
  }

  /**
   * Copy constructor.
   *
   * @param source the source
   */
  protected AreaAverageFilter(AreaAverageFilter source) {
    super(source);
    sumFilter1 = source.sumFilter1.copy();
    sumFilter2 = source.sumFilter2.copy();
    blockMeanFilter1 = source.blockMeanFilter1.copy();
    blockMeanFilter2 = source.blockMeanFilter2.copy();
    simpleInterpolation = source.simpleInterpolation;
  }

  /**
   * Create a copy.
   *
   * @return the copy
   */
  public AreaAverageFilter copy() {
    return new AreaAverageFilter(this);
  }

  /**
   * Compute the block average within a 2w+1 size block around each point. Pixels within border
   * regions (width = ceil(w)) are unchanged.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   * @param width The width
   */
  public void areaAverageUsingSumsInternal(float[] data, final int maxx, final int maxy,
      final double width) {
    if (width <= 0) {
      return;
    }

    final int n = (int) width;
    final int n1 = (n == width) ? n : n + 1;
    final int blockSize = 2 * n1 + 1;
    if (maxx < blockSize || maxy < blockSize) {
      return;
    }

    if (n == n1) {
      // There is no edge
      blockMeanFilter1.rollingBlockFilterInternal(data, maxx, maxy, n);
      return;
    }

    // Calculate the sum in the n and n+1 regions
    final float[] sum1;
    if (n == 0) {
      sum1 = data;
    } else {
      sum1 = data.clone();
      sumFilter1.rollingBlockFilter(sum1, maxx, maxy, n);
    }
    final float[] sum2 = data.clone();
    sumFilter2.rollingBlockFilterInternal(sum2, maxx, maxy, n1);

    // Get the average by adding the inner sum to the weighted edge pixels.
    final double area = (2 * width + 1) * (2 * width + 1);

    final float edgeWeight;

    if (simpleInterpolation) {
      edgeWeight = (float) (width - n);
    } else {
      // Use the area to produce the weighting
      final int inner = (2 * n + 1) * (2 * n + 1);
      final int outer = (2 * n1 + 1) * (2 * n1 + 1);
      edgeWeight = (float) ((area - inner) / (outer - inner));
    }

    final float norm = (float) (1.0 / area);

    for (int y = n1; y < maxy - n1; y++) {
      int index = y * maxx + n1;
      for (int x = n1; x < maxx - n1; x++, index++) {
        data[index] = norm * (sum1[index] + edgeWeight * (sum2[index] - sum1[index]));
      }
    }
  }

  /**
   * Compute the block average within a 2w+1 size block around each point.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   * @param width The width
   */
  public void areaAverageUsingSums(float[] data, final int maxx, final int maxy,
      final double width) {
    if (width <= 0) {
      return;
    }

    final int n = (int) width;
    final int n1 = (n == width) ? n : n + 1;

    if (n == n1 || (maxx < n1 && maxy < n1)) {
      // There is no edge
      blockMeanFilter1.rollingBlockFilter(data, maxx, maxy, n);
      return;
    }

    // Calculate the sum in the n and n+1 regions
    final float[] sum1;
    if (n == 0) {
      sum1 = data;
    } else {
      sum1 = data.clone();
      sumFilter1.rollingBlockFilter(sum1, maxx, maxy, n);
    }
    final float[] sum2 = data.clone();
    sumFilter2.rollingBlockFilter(sum2, maxx, maxy, n1);

    // Get the average by adding the inner sum to the weighted edge pixels.
    final double area = (2 * width + 1) * (2 * width + 1);

    final float edgeWeight;

    if (simpleInterpolation) {
      edgeWeight = (float) (width - n);
    } else {
      // Use the area to produce the weighting
      final int inner = (2 * n + 1) * (2 * n + 1);
      final int outer = (2 * n1 + 1) * (2 * n1 + 1);
      edgeWeight = (float) ((area - inner) / (outer - inner));
    }

    final float norm = (float) (1.0 / area);

    for (int index = 0; index < sum1.length; index++) {
      data[index] = norm * (sum1[index] + edgeWeight * (sum2[index] - sum1[index]));
    }
  }

  /**
   * Compute the block average within a 2w+1 size block around each point. Pixels within border
   * regions (width = ceil(w)) are unchanged.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   * @param width The width
   */
  public void areaAverageUsingAveragesInternal(float[] data, final int maxx, final int maxy,
      final double width) {
    if (width <= 0) {
      return;
    }

    final int n = (int) width;
    final int n1 = (n == width) ? n : n + 1;
    final int blockSize = 2 * n1 + 1;
    if (maxx < blockSize || maxy < blockSize) {
      return;
    }

    if (n == n1) {
      // There is no edge
      blockMeanFilter1.rollingBlockFilterInternal(data, maxx, maxy, n);
      return;
    }

    // Calculate the sum in the n and n+1 regions
    final float[] av1;
    if (n == 0) {
      av1 = data;
    } else {
      av1 = data.clone();
      blockMeanFilter1.rollingBlockFilter(av1, maxx, maxy, n);
    }
    final float[] av2 = data.clone();
    blockMeanFilter2.rollingBlockFilterInternal(av2, maxx, maxy, n1);

    // Get the average by weighting the two
    final float outerWeight;

    if (simpleInterpolation) {
      outerWeight = (float) (width - n);
    } else {
      // Use the area to produce the weighting
      final int inner = (2 * n + 1) * (2 * n + 1);
      final int outer = (2 * n1 + 1) * (2 * n1 + 1);
      final double area = (2 * width + 1) * (2 * width + 1);
      outerWeight = (float) ((area - inner) / (outer - inner));
    }

    final float innerWeight = 1 - outerWeight;

    for (int y = n1; y < maxy - n1; y++) {
      int index = y * maxx + n1;
      for (int x = n1; x < maxx - n1; x++, index++) {
        data[index] = av1[index] * innerWeight + av2[index] * outerWeight;
      }
    }
  }

  /**
   * Compute the block average within a 2w+1 size block around each point.
   *
   * <p>Note: the input data is destructively modified
   *
   * @param data The input/output data (packed in YX order)
   * @param maxx The width of the data
   * @param maxy The height of the data
   * @param width The width
   */
  public void areaAverageUsingAverages(float[] data, final int maxx, final int maxy,
      final double width) {
    if (width <= 0) {
      return;
    }

    final int n = (int) width;
    final int n1 = (n == width) ? n : n + 1;

    if (n == n1 || (maxx < n1 && maxy < n1)) {
      // There is no edge
      blockMeanFilter1.rollingBlockFilter(data, maxx, maxy, n);
      return;
    }

    // Calculate the sum in the n and n+1 regions
    final float[] av1;
    if (n == 0) {
      av1 = data;
    } else {
      av1 = data.clone();
      blockMeanFilter1.rollingBlockFilter(av1, maxx, maxy, n);
    }
    final float[] av2 = data.clone();
    blockMeanFilter2.rollingBlockFilter(av2, maxx, maxy, n1);

    // Get the average by weighting the two
    final float outerWeight;

    if (simpleInterpolation) {
      outerWeight = (float) (width - n);
    } else {
      // Use the area to produce the weighting
      final int inner = (2 * n + 1) * (2 * n + 1);
      final int outer = (2 * n1 + 1) * (2 * n1 + 1);
      final double area = (2 * width + 1) * (2 * width + 1);
      outerWeight = (float) ((area - inner) / (outer - inner));
    }

    final float innerWeight = 1 - outerWeight;

    for (int index = 0; index < av1.length; index++) {
      data[index] = av1[index] * innerWeight + av2[index] * outerWeight;
    }
  }

  /**
   * Checks if using simple interpolation.
   *
   * @return true for simple interpolation.
   */
  public boolean isSimpleInterpolation() {
    return simpleInterpolation;
  }

  /**
   * The average for block size n and n+1 is linearly interpolated. Set this to true to use a weight
   * of (w-n) for the outer average. Set to false to use a weight based on the area of the edge
   * pixels in the (2w+1) region.
   *
   * @param simpleInterpolation true for simple interpolation
   */
  public void setSimpleInterpolation(boolean simpleInterpolation) {
    this.simpleInterpolation = simpleInterpolation;
  }

  /** {@inheritDoc} */
  @Override
  protected void newWeights() {
    sumFilter1.setWeights(weights, weightWidth, weightHeight);
    sumFilter2.setWeights(weights, weightWidth, weightHeight);
    blockMeanFilter1.setWeights(weights, weightWidth, weightHeight);
    blockMeanFilter2.setWeights(weights, weightWidth, weightHeight);
  }
}
