/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2022 Alex Herbert
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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Define a class to pre-process the data, ignoring a defined border.
 */
public abstract class DataProcessor {
  private final int border;

  /**
   * Instantiates a new data processor.
   *
   * @param border The border that can be ignored
   */
  public DataProcessor(int border) {
    this.border = border;
  }

  /**
   * Copy constructor.
   *
   * @param source the source
   */
  protected DataProcessor(DataProcessor source) {
    border = source.border;
  }

  /**
   * Checks if the data processor is weighted, i.e. supports {@link #setWeights(float[], int, int)}.
   *
   * @return true, if is weighted
   */
  public abstract boolean isWeighted();

  /**
   * Sets the weights of the data. This should be called before {@link #process (float[], int, int)}
   * is called with data samples.
   *
   * <p>Calling this in advance allows efficient caching of pre-computed weightings.
   *
   * @param weights the weights of the data (can be null)
   * @param width The width of the data
   * @param height The height of the data
   */
  public abstract void setWeights(float[] weights, int width, int height);

  /**
   * Checks for weights. Weights are set using {@link #setWeights(float[], int, int)}.
   *
   * @return true, if successful
   */
  public abstract boolean hasWeights();

  /**
   * Process the data.
   *
   * @param data The data
   * @param width The width of the data
   * @param height The height of the data
   * @return The new data
   */
  public abstract float[] process(float[] data, int width, int height);

  /**
   * Gets the border.
   *
   * @return the border.
   */
  public int getBorder() {
    return border;
  }

  /**
   * Create a copy.
   *
   * @return the copy
   */
  public abstract DataProcessor copy();

  /**
   * Gets the description of the processor and parameters.
   *
   * @return A description of the processor and parameters.
   */
  public String getDescription() {
    return getName() + ": " + Arrays.toString(getParameters().toArray());
  }

  /**
   * Gets the name of the filter.
   *
   * @return The name of the filter.
   */
  public abstract String getName();

  /**
   * Gets the parameters of the filter.
   *
   * @return The parameters of the filter.
   */
  public List<String> getParameters() {
    final ArrayList<String> list = new ArrayList<>();
    list.add("border = " + border);
    if (hasWeights()) {
      list.add("weighted");
    }
    return list;
  }

  /**
   * Get the width spread of data used to process each position.
   *
   * @return The spread
   */
  public abstract double getSpread();
}
