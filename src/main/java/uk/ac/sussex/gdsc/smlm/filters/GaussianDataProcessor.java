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
 * Identifies candidate spots (local maxima) in an image. The image is smoothed with a Gaussian
 * filter.
 */
public class GaussianDataProcessor extends DataProcessor {
  private final double sigma;
  private final GaussianFilter filter;

  /**
   * Constructor.
   *
   * @param border The border to ignore for maxima
   * @param smooth The distance into neighbouring pixels to extend. The resulting standard deviation
   *        can be found using {@link #getSigma()}
   */
  public GaussianDataProcessor(int border, double smooth) {
    super(border);
    this.sigma = getSigma(smooth);
    filter = new GaussianFilter(0.02);
  }

  /**
   * Copy constructor.
   *
   * @param source the source
   */
  protected GaussianDataProcessor(GaussianDataProcessor source) {
    super(source);
    sigma = source.sigma;
    filter = source.filter.copy();
  }

  @Override
  public GaussianDataProcessor copy() {
    return new GaussianDataProcessor(this);
  }

  /**
   * Get the Gaussian standard deviation for the desired smoothing distance.
   *
   * @param smooth the smoothing distance
   * @return the Gaussian standard deviation for the desired smoothing distance.
   */
  public static double getSigma(double smooth) {
    if (smooth < 0) {
      return 0;
    }
    return smooth;
  }

  /**
   * Gets the Gaussian standard deviation (sigma).
   *
   * @return the Gaussian standard deviation.
   */
  public double getSigma() {
    return sigma;
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
    if (sigma > 0) {
      // Smoothing destructively modifies the data so create a copy
      smoothData = Arrays.copyOf(data, width * height);
      if (GaussianFilter.getBorder(sigma) <= getBorder()) {
        filter.convolveInternal(smoothData, width, height, sigma);
      } else {
        filter.convolve(smoothData, width, height, sigma);
      }
    }
    return smoothData;
  }

  @Override
  public String getName() {
    return "Gaussian";
  }

  @Override
  public List<String> getParameters() {
    final List<String> list = super.getParameters();
    list.add("sigma = " + MathUtils.rounded(sigma));
    list.add("width = " + MathUtils.rounded(filter.getHalfWidth(sigma)));
    return list;
  }

  @Override
  public double getSpread() {
    return 6 * sigma;
  }
}
