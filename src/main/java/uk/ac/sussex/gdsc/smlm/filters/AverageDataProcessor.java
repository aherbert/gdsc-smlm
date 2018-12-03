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

import java.util.Arrays;
import java.util.List;

import uk.ac.sussex.gdsc.core.utils.MathUtils;

/**
 * Identifies candidate spots (local maxima) in an image. The image is smoothed with an average box
 * filter. When the box size is large then the smoothing switches to using an interpolation between
 * two block sizes. This is an approximation due to incorrect weighting of the corners. Note that at
 * large sizes the corners are a fraction of the total edge pixels so the difference is minor.
 *
 * @see uk.ac.sussex.gdsc.smlm.filters.AreaAverageFilter
 */
public class AverageDataProcessor extends DataProcessor {
  /**
   * Define the default smoothing size at which the smoothing switches to using an interpolation
   * between two block sizes. Below this limit the smoothing uses an exact mean filter.
   */
  public static final int AREA_FILTER_LIMIT = 3;

  private final float smooth;
  private final int iSmooth;
  private BlockMeanFilter blockMeanFilter = null;
  private AreaAverageFilter areaAverageFilter = null;

  /**
   * Constructor.
   *
   * @param border The border to ignore for maxima
   * @param smooth The smoothing width to apply to the data
   */
  public AverageDataProcessor(int border, double smooth) {
    this(border, smooth, AREA_FILTER_LIMIT);
  }

  /**
   * Constructor.
   *
   * @param border The border to ignore for maxima
   * @param smooth The smoothing width to apply to the data
   * @param areaFilterLimit The limit to switch to simple area interpolation
   */
  public AverageDataProcessor(int border, double smooth, int areaFilterLimit) {
    super(border);
    this.smooth = (float) convert(smooth);
    // Store the smoothing value as an integer
    iSmooth = ((int) smooth == smooth) ? (int) smooth : 0;

    // Only create the area filter if we need it
    if (iSmooth > 0) {
      blockMeanFilter = new BlockMeanFilter();
    } else if (smooth > areaFilterLimit) {
      areaAverageFilter = new AreaAverageFilter();
    } else {
      blockMeanFilter = new BlockMeanFilter();
    }
  }

  /**
   * Convert the smoothing parameter to the value which is used for the AreaAverageFilter. Values
   * below zero are set to zero.
   *
   * @param smooth the smoothing parameter
   * @return The adjusted value
   */
  public static double convert(double smooth) {
    if (smooth < 0) {
      return 0;
    }
    return smooth;
  }

  /** {@inheritDoc} */
  @Override
  public boolean isWeighted() {
    return true;
  }

  /** {@inheritDoc} */
  @Override
  public void setWeights(float[] weights, int width, int height) {
    final BaseWeightedFilter f = getFilter();
    if (f != null) {
      f.setWeights(weights, width, height);
    }
  }

  /** {@inheritDoc} */
  @Override
  public boolean hasWeights() {
    final BaseWeightedFilter f = getFilter();
    return (f != null) ? f.hasWeights() : false;
  }

  /** {@inheritDoc} */
  @Override
  public float[] process(float[] data, int width, int height) {
    float[] smoothData = data;
    if (smooth > 0) {
      // Smoothing destructively modifies the data so create a copy
      smoothData = Arrays.copyOf(data, width * height);

      // ADH 05-Jan-2017:
      // This was changed from 1 to 0. Previously if the iSmooth was 1 then
      // it would fall through to the striped block filter using a weight of 1.
      // This can be done using the rolling block algorithm instead.
      if (iSmooth > 0) {
        // Integer smoothing is faster using a rolling block algorithm
        if (smooth <= getBorder()) {
          blockMeanFilter.rollingBlockFilterInternal(smoothData, width, height, iSmooth);
        } else {
          blockMeanFilter.rollingBlockFilter(smoothData, width, height, iSmooth);
        }
      } else if (areaAverageFilter != null) {
        if (smooth <= getBorder()) {
          areaAverageFilter.areaAverageUsingSumsInternal(smoothData, width, height, smooth);
        } else {
          areaAverageFilter.areaAverageUsingSums(smoothData, width, height, smooth);
        }
      } else if (smooth <= getBorder()) {
        blockMeanFilter.stripedBlockFilterInternal(smoothData, width, height, smooth);
      } else {
        blockMeanFilter.stripedBlockFilter(smoothData, width, height, smooth);
      }
    }
    return smoothData;
  }

  private BaseWeightedFilter getFilter() {
    if (smooth > 0) {
      if (iSmooth > 1) {
        return blockMeanFilter;
      } else if (areaAverageFilter != null) {
        return areaAverageFilter;
      } else {
        return blockMeanFilter;
      }
    }
    return null;
  }

  /**
   * Gets the smooth.
   *
   * @return the smoothing width
   */
  public double getSmooth() {
    return smooth;
  }

  /** {@inheritDoc} */
  @Override
  public AverageDataProcessor clone() {
    final AverageDataProcessor f = (AverageDataProcessor) super.clone();
    // Ensure the object is duplicated and not passed by reference.
    if (blockMeanFilter != null) {
      f.blockMeanFilter = blockMeanFilter.clone();
    } else {
      f.areaAverageFilter = areaAverageFilter.clone();
    }
    return f;
  }

  /** {@inheritDoc} */
  @Override
  public String getName() {
    return "Average";
  }

  /** {@inheritDoc} */
  @Override
  public List<String> getParameters() {
    final List<String> list = super.getParameters();
    list.add("smooth = " + MathUtils.rounded(smooth));
    return list;
  }

  /** {@inheritDoc} */
  @Override
  public double getSpread() {
    return 2 * smooth + 1;
  }
}
