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

package uk.ac.sussex.gdsc.smlm.results;

/**
 * Specifies a peak fitting result with assignable attributes. If the attributes are not set then
 * the default values are returned.
 */
public class AttributePeakResult extends PeakResult {
  // Provide assignable attributes for all the attribute methods in the super class

  // State flags
  private static final int FIELD_ID = 0x00000001;
  private static final int FIELD_END_FRAME = 0x00000002;
  private static final int FIELD_PRECISION = 0x00000004;
  private int fields;

  private int id;
  private int endFrame;
  private float precision;

  /** {@inheritDoc} */
  @Override
  public boolean hasId() {
    return ((fields & FIELD_ID) == FIELD_ID);
  }

  /** {@inheritDoc} */
  @Override
  public boolean hasEndFrame() {
    return ((fields & FIELD_END_FRAME) == FIELD_END_FRAME);
  }

  /** {@inheritDoc} */
  @Override
  public boolean hasPrecision() {
    return ((fields & FIELD_PRECISION) == FIELD_PRECISION);
  }

  private void setHasId() {
    fields |= FIELD_ID;
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

  /** {@inheritDoc} */
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
    // if (id != super.getId())
    // {
    setHasId();
    this.id = id;
    // }
    // else
    // clearHasId();
  }

  /** {@inheritDoc} */
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

  /** {@inheritDoc} */
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
    super(peakResult.getFrame(), peakResult.getOrigX(), peakResult.getOrigY(),
        peakResult.getOrigValue(), peakResult.getError(), peakResult.getNoise(),
        peakResult.getMeanIntensity(), peakResult.getParameters().clone(),
        (peakResult.getParameterDeviations() == null) ? null
            : peakResult.getParameterDeviations().clone());
    if (peakResult.hasId()) {
      setId(peakResult.getId());
    }
    if (peakResult.hasEndFrame()) {
      setEndFrame(peakResult.getEndFrame());
    }
    if (peakResult.hasPrecision()) {
      setPrecision(peakResult.getPrecision());
    }
  }
}
