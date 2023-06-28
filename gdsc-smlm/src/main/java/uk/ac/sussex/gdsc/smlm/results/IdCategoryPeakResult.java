/*-
 * #%L
 * Genome Damage and Stability Centre SMLM Package
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2023 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.results;

/**
 * Specifies a peak fitting result with an Id and category.
 */
public class IdCategoryPeakResult extends IdPeakResult {
  private static final long serialVersionUID = 20210510L;

  private final int category;

  /**
   * Instantiates a new peak result.
   *
   * @param frame the frame
   * @param origX the original X position
   * @param origY the original Y position
   * @param origValue the original value
   * @param error the error
   * @param noise the noise
   * @param meanIntensity the mean intensity
   * @param params the params (must not be null and must have at least
   *        {@value uk.ac.sussex.gdsc.smlm.results.PeakResult#STANDARD_PARAMETERS} parameters)
   * @param paramsStdDev the params standard deviations (if not null must match the length of the
   *        params array)
   * @param id the id
   * @param category the category
   * @throws IllegalArgumentException the illegal argument exception if the parameters are invalid
   */
  public IdCategoryPeakResult(int frame, int origX, int origY, float origValue, double error,
      float noise, float meanIntensity, float[] params, float[] paramsStdDev, int id,
      int category) {
    super(frame, origX, origY, origValue, error, noise, meanIntensity, params, paramsStdDev, id);
    this.category = category;
  }

  /**
   * Copy constructor.
   *
   * @param source the source
   */
  protected IdCategoryPeakResult(IdCategoryPeakResult source) {
    super(source);
    this.category = source.category;
  }

  @Override
  public IdCategoryPeakResult copy() {
    return new IdCategoryPeakResult(this);
  }

  @Override
  public boolean hasCategory() {
    return true;
  }

  @Override
  public int getCategory() {
    return category;
  }
}
