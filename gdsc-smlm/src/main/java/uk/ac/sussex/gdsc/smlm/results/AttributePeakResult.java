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
 * Specifies a peak fitting result with assignable attributes. If the attributes are not set then
 * the default values are returned.
 */
public class AttributePeakResult extends PeakResult {
  private static final long serialVersionUID = 20190319L;

  // Provide assignable attributes for all the attribute methods in the super class

  // State flags
  private static final int FIELD_ID = 0x01;
  private static final int FIELD_END_FRAME = 0x02;
  private static final int FIELD_PRECISION = 0x04;
  private static final int FIELD_CATEGORY = 0x08;
  private int fields;

  private int id;
  private int category;
  private int endFrame;
  private float precision;

  @Override
  public boolean hasId() {
    return ((fields & FIELD_ID) == FIELD_ID);
  }

  @Override
  public boolean hasEndFrame() {
    return ((fields & FIELD_END_FRAME) == FIELD_END_FRAME);
  }

  @Override
  public boolean hasPrecision() {
    return ((fields & FIELD_PRECISION) == FIELD_PRECISION);
  }

  @Override
  public boolean hasCategory() {
    return ((fields & FIELD_CATEGORY) == FIELD_CATEGORY);
  }

  private void setHasId() {
    fields |= FIELD_ID;
  }

  private void setHasCategory() {
    fields |= FIELD_CATEGORY;
  }

  private void setHasEndFrame() {
    fields |= FIELD_END_FRAME;
  }

  private void setHasPrecision() {
    fields |= FIELD_PRECISION;
  }

  /**
   * Clear has id.
   */
  public void clearHasId() {
    fields = fields & ~FIELD_ID;
  }

  /**
   * Clear has category.
   */
  public void clearHasCategory() {
    fields = fields & ~FIELD_CATEGORY;
  }

  /**
   * Clear has end frame.
   */
  public void clearHasEndFrame() {
    fields = fields & ~FIELD_END_FRAME;
  }

  /**
   * Clear has precision.
   */
  public void clearHasPrecision() {
    fields = fields & ~FIELD_PRECISION;
  }

  @Override
  public int getId() {
    if (hasId()) {
      return id;
    }
    return super.getId();
  }

  /**
   * Sets the id.
   *
   * @param id the new id
   */
  public void setId(int id) {
    // Allow ID to be anything, including zero
    setHasId();
    this.id = id;
  }

  @Override
  public int getCategory() {
    if (hasCategory()) {
      return category;
    }
    return super.getCategory();
  }

  /**
   * Sets the category.
   *
   * @param category the new category
   */
  public void setCategory(int category) {
    // Allow category to be anything, including zero
    setHasCategory();
    this.category = category;
  }

  @Override
  public int getEndFrame() {
    if (hasEndFrame()) {
      return endFrame;
    }
    return super.getFrame();
  }

  /**
   * Sets the end frame.
   *
   * @param endFrame the new end frame
   */
  public void setEndFrame(int endFrame) {
    // End frame must be after the start frame.
    if (endFrame > super.getFrame()) {
      setHasEndFrame();
      this.endFrame = endFrame;
    } else {
      clearHasEndFrame();
    }
  }

  @Override
  public void setFrame(int frame) {
    // Set the new start frame
    super.setFrame(frame);
    // Validate the current end frame
    setEndFrame(endFrame);
  }

  @Override
  public double getPrecision() {
    if (hasPrecision()) {
      return precision;
    }
    return super.getPrecision();
  }

  /**
   * Sets the precision.
   *
   * @param precision the new precision
   */
  public void setPrecision(double precision) {
    if (precision >= 0) {
      setHasPrecision();
      this.precision = (float) precision;
    } else {
      clearHasPrecision();
    }
  }

  /**
   * Instantiates a new attribute peak result.
   *
   * @param startFrame the start frame
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
   * @throws IllegalArgumentException the illegal argument exception if the parameters are invalid
   */
  public AttributePeakResult(int startFrame, int origX, int origY, float origValue, double error,
      float noise, float meanIntensity, float[] params, float[] paramsStdDev) {
    super(startFrame, origX, origY, origValue, error, noise, meanIntensity, params, paramsStdDev);
  }

  /**
   * Instantiates a new attribute peak result.
   *
   * @param frame the frame
   * @param x the x position
   * @param y the y position
   * @param intensity the intensity
   */
  public AttributePeakResult(int frame, float x, float y, float intensity) {
    super(frame, x, y, intensity);
  }

  /**
   * Instantiates a new attribute peak result.
   *
   * @param x the x position
   * @param y the y position
   * @param intensity the intensity
   */
  public AttributePeakResult(float x, float y, float intensity) {
    super(x, y, intensity);
  }

  /**
   * Instantiates a new attribute peak result copying all the attributes from the result. This is a
   * deep copy of all the result data.
   *
   * @param peakResult the peak result
   */
  public AttributePeakResult(PeakResult peakResult) {
    super(peakResult);
    if (peakResult.hasId()) {
      setId(peakResult.getId());
    }
    if (peakResult.hasCategory()) {
      setCategory(peakResult.getCategory());
    }
    if (peakResult.hasEndFrame()) {
      setEndFrame(peakResult.getEndFrame());
    }
    if (peakResult.hasPrecision()) {
      setPrecision(peakResult.getPrecision());
    }
  }

  @Override
  public AttributePeakResult copy() {
    return new AttributePeakResult(this);
  }
}
