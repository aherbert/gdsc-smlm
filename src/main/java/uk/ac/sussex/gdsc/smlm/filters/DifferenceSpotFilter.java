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

import uk.ac.sussex.gdsc.core.utils.ValidationUtils;

import java.util.List;

/**
 * Identifies candidate spots (local maxima) in an image. The image is pre-processed with two
 * filters and the second subtracted from the first.
 */
public class DifferenceSpotFilter extends MaximaSpotFilter {
  private final DataProcessor processor1;
  private final DataProcessor processor2;

  /**
   * Constructor.
   *
   * @param search The search width for non-maximum suppression
   * @param border The border to ignore for maxima
   * @param processor1 The first data processor
   * @param processor2 The second data processor
   * @throws IllegalArgumentException if the spread of the second processor is smaller than the
   *         first
   */
  public DifferenceSpotFilter(int search, int border, DataProcessor processor1,
      DataProcessor processor2) {
    super(search, border);
    this.processor1 = ValidationUtils.checkNotNull(processor1, "Processor 1 is null");
    this.processor2 = ValidationUtils.checkNotNull(processor2, "Processor 2 is null");
    // This is a simple protection from invalid difference-of-smoothing.
    ValidationUtils.checkArgument(processor2.getSpread() > processor1.getSpread(),
        "Processor 2 acts on a smaller spread of data than processor 1");
  }

  /**
   * Copy constructor.
   *
   * @param source the source
   */
  protected DifferenceSpotFilter(DifferenceSpotFilter source) {
    super(source);
    processor1 = source.processor1.copy();
    processor2 = source.processor2.copy();
  }

  @Override
  public DifferenceSpotFilter copy() {
    return new DifferenceSpotFilter(this);
  }

  @Override
  public boolean isAbsoluteIntensity() {
    return false;
  }

  @Override
  public boolean isWeighted() {
    return processor1.isWeighted() || processor2.isWeighted();
  }

  @Override
  public void setWeights(float[] weights, int width, int height) {
    processor1.setWeights(weights, width, height);
    processor2.setWeights(weights, width, height);
  }

  @Override
  public boolean hasWeights() {
    return processor1.hasWeights() || processor2.hasWeights();
  }

  @Override
  public float[] preprocessData(float[] data, int width, int height) {
    final float[] data1 = processor1.process(data, width, height);
    final float[] data2 = processor2.process(data, width, height);
    for (int i = 0; i < data1.length; i++) {
      data1[i] -= data2[i];
    }
    return data1;
  }

  @Override
  public String getName() {
    return "Difference";
  }

  @Override
  public List<String> getParameters() {
    final List<String> list = super.getParameters();
    list.add("Filter 1 = " + processor1.getDescription());
    list.add("Filter 2 = " + processor2.getDescription());
    return list;
  }

  @Override
  public double getSpread() {
    return Math.max(processor1.getSpread(), processor2.getSpread());
  }
}
