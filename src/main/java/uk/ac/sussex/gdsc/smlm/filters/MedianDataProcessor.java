/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2020 Alex Herbert
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

import java.util.Arrays;
import java.util.List;
import org.apache.commons.math3.util.FastMath;

/**
 * Identifies candidate spots (local maxima) in an image. The image is smoothed with an median box
 * filter.
 */
public class MedianDataProcessor extends DataProcessor {
  private final int smooth;
  private final MedianFilter filter;

  /**
   * Constructor.
   *
   * @param border The border to ignore for maxima
   * @param smooth The smoothing width to apply to the data
   */
  public MedianDataProcessor(int border, double smooth) {
    super(border);
    this.smooth = convert(smooth);
    filter = new MedianFilter();
  }

  /**
   * Copy constructor.
   *
   * @param source the source
   */
  protected MedianDataProcessor(MedianDataProcessor source) {
    super(source);
    smooth = source.smooth;
    filter = source.filter.copy();
  }

  @Override
  public MedianDataProcessor copy() {
    return new MedianDataProcessor(this);
  }

  /**
   * Convert the smoothing parameter to the value which is used for the MedianFilter. We only use
   * int smoothing. Values below zero are set to zero.
   *
   * @param smooth the smoothing parameter
   * @return The adjusted value
   * @see MedianFilter
   */
  public static int convert(double smooth) {
    if (smooth < 0) {
      return 0;
    }
    return (int) smooth;
  }

  @Override
  public boolean isWeighted() {
    return false;
  }

  @Override
  public void setWeights(float[] weights, int width, int height) {
    // Do nothing
  }

  @Override
  public boolean hasWeights() {
    return false;
  }

  @Override
  public float[] process(float[] data, int width, int height) {
    float[] smoothData = data;
    if (smooth > 0) {
      // Smoothing destructively modifies the data so create a copy
      smoothData = Arrays.copyOf(data, width * height);

      // Check upper limits are safe
      final int tmpSmooth = FastMath.min(smooth, FastMath.min(width, height) / 2);

      // JUnit speed tests show that the rolling median is not faster.
      // It used to be faster on windows less than 3.

      // if (tmpSmooth <= 3)
      // {
      // if (tmpSmooth <= getBorder())
      // {
      // filter.rollingMedianInternal(smoothData, width, height, tmpSmooth);
      // }
      // else
      // {
      // filter.rollingMedian(smoothData, width, height, tmpSmooth);
      // }
      // }
      // else
      // {
      if (tmpSmooth <= getBorder()) {
        filter.blockMedianInternal(smoothData, width, height, tmpSmooth);
      } else {
        filter.blockMedian(smoothData, width, height, tmpSmooth);
      }
    }
    return smoothData;
  }

  /**
   * Gets the smoothing width.
   *
   * @return the smoothing width.
   */
  public int getSmooth() {
    return smooth;
  }

  @Override
  public String getName() {
    return "Median";
  }

  @Override
  public List<String> getParameters() {
    final List<String> list = super.getParameters();
    list.add("smooth = " + smooth);
    return list;
  }

  @Override
  public double getSpread() {
    return 2 * smooth + 1;
  }
}
