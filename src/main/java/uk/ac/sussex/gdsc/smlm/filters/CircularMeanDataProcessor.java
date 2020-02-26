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

import java.util.Arrays;
import java.util.List;
import uk.ac.sussex.gdsc.core.utils.MathUtils;

/**
 * Identifies candidate spots (local maxima) in an image. The image is smoothed with a circular mean
 * filter.
 */
public class CircularMeanDataProcessor extends DataProcessor {
  private final double radius;
  private final CircularMeanFilter filter;

  /**
   * Constructor.
   *
   * @param border The border to ignore for maxima
   * @param smooth The distance into neighbouring pixels to extend
   */
  public CircularMeanDataProcessor(int border, double smooth) {
    super(border);
    this.radius = getSigma(smooth);
    filter = new CircularMeanFilter();
  }

  /**
   * Copy constructor.
   *
   * @param source the source
   */
  protected CircularMeanDataProcessor(CircularMeanDataProcessor source) {
    super(source);
    radius = source.radius;
    filter = source.filter.copy();
  }

  @Override
  public CircularMeanDataProcessor copy() {
    return new CircularMeanDataProcessor(this);
  }

  /**
   * Get the radius for the desired smoothing distance.
   *
   * @param smooth the smoothing distance
   * @return the radius for the desired smoothing distance.
   */
  public static double getSigma(double smooth) {
    if (smooth < 0) {
      return 0;
    }
    return smooth;
  }

  @Override
  public boolean isWeighted() {
    return true;
  }

  @Override
  public void setWeights(float[] weights, int width, int height) {
    filter.setWeights(weights, width, height);
  }

  @Override
  public boolean hasWeights() {
    return filter.hasWeights();
  }

  @Override
  public float[] process(float[] data, int width, int height) {
    float[] smoothData = data;
    if (radius > 0) {
      // Smoothing destructively modifies the data so create a copy
      smoothData = Arrays.copyOf(data, width * height);
      if (CircularFilter.getBorder(radius) <= getBorder()) {
        filter.convolveInternal(smoothData, width, height, radius);
      } else {
        filter.convolve(smoothData, width, height, radius);
      }
    }
    return smoothData;
  }

  /**
   * Gets the radius.
   *
   * @return the smoothing radius
   */
  public double getRadius() {
    return radius;
  }

  @Override
  public String getName() {
    return "Circular Mean";
  }

  @Override
  public List<String> getParameters() {
    final List<String> list = super.getParameters();
    list.add("radius = " + MathUtils.rounded(radius));
    return list;
  }

  @Override
  public double getSpread() {
    return CircularFilter.getPixelRadius(radius) * 2;
  }
}
