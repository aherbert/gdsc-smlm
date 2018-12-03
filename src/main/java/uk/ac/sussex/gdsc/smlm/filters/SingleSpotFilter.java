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

import java.util.List;

/**
 * Identifies candidate spots (local maxima) in an image. The image is pre-processed with a single
 * filter.
 */
public class SingleSpotFilter extends MaximaSpotFilter {
  private DataProcessor processor;

  /**
   * Constructor.
   *
   * @param search The search width for non-maximum suppression
   * @param border The border to ignore for maxima
   * @param processor The data processor
   * @throws IllegalArgumentException if processor is null
   */
  public SingleSpotFilter(int search, int border, DataProcessor processor) {
    super(search, border);
    if (processor == null) {
      throw new IllegalArgumentException("Processor is null");
    }
    this.processor = processor;
  }

  /** {@inheritDoc} */
  @Override
  public boolean isAbsoluteIntensity() {
    return true;
  }

  /** {@inheritDoc} */
  @Override
  public boolean isWeighted() {
    return processor.isWeighted();
  }

  /** {@inheritDoc} */
  @Override
  public void setWeights(float[] weights, int width, int height) {
    processor.setWeights(weights, width, height);
  }

  /** {@inheritDoc} */
  @Override
  public boolean hasWeights() {
    return processor.hasWeights();
  }

  /** {@inheritDoc} */
  @Override
  public float[] preprocessData(float[] data, int width, int height) {
    return processor.process(data, width, height);
  }

  /** {@inheritDoc} */
  @Override
  public SingleSpotFilter clone() {
    final SingleSpotFilter f = (SingleSpotFilter) super.clone();
    // Ensure the object is duplicated and not passed by reference.
    f.processor = processor.clone();
    return f;
  }

  /** {@inheritDoc} */
  @Override
  public String getName() {
    return "Single";
  }

  /** {@inheritDoc} */
  @Override
  public List<String> getParameters() {
    final List<String> list = super.getParameters();
    list.add("Filter = " + processor.getDescription());
    return list;
  }

  /** {@inheritDoc} */
  @Override
  public double getSpread() {
    return processor.getSpread();
  }
}
