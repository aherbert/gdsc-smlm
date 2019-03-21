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

package uk.ac.sussex.gdsc.smlm.results;

/**
 * Specifies a peak fitting result that cannot be modified.
 *
 * <p>Any method that modifies the result will throw a data exception.
 */
public class ImmutablePeakResult extends AttributePeakResult {
  private static final long serialVersionUID = 20190319L;

  private static final String IMMUTABLE_MESSAGE = "This result is immutable";

  private final boolean built;

  /**
   * Instantiates a new peak result.
   *
   * @param peakResult the peak result
   * @throws IllegalArgumentException if the parameters are invalid
   */
  public ImmutablePeakResult(PeakResult peakResult) {
    super(peakResult);
    built = true;
  }

  @Override
  public void setBackground(float background) {
    throw new UnsupportedOperationException(IMMUTABLE_MESSAGE);
  }

  @Override
  public void setIntensity(float intensity) {
    throw new UnsupportedOperationException(IMMUTABLE_MESSAGE);
  }

  @Override
  public void setXPosition(float x) {
    throw new UnsupportedOperationException(IMMUTABLE_MESSAGE);
  }

  @Override
  public void setYPosition(float y) {
    throw new UnsupportedOperationException(IMMUTABLE_MESSAGE);
  }

  @Override
  public void setZPosition(float z) {
    throw new UnsupportedOperationException(IMMUTABLE_MESSAGE);
  }

  @Override
  public void setFrame(int frame) {
    throw new UnsupportedOperationException(IMMUTABLE_MESSAGE);
  }

  @Override
  public void setOrigX(int origX) {
    throw new UnsupportedOperationException(IMMUTABLE_MESSAGE);
  }

  @Override
  public void setOrigY(int origY) {
    throw new UnsupportedOperationException(IMMUTABLE_MESSAGE);
  }

  @Override
  public void setOrigValue(float origValue) {
    throw new UnsupportedOperationException(IMMUTABLE_MESSAGE);
  }

  @Override
  public void setError(double error) {
    throw new UnsupportedOperationException(IMMUTABLE_MESSAGE);
  }

  @Override
  public void setNoise(float noise) {
    throw new UnsupportedOperationException(IMMUTABLE_MESSAGE);
  }

  @Override
  public void setId(int id) {
    if (built) {
      throw new UnsupportedOperationException(IMMUTABLE_MESSAGE);
    }
    super.setId(id);
  }

  @Override
  public void setEndFrame(int endFrame) {
    if (built) {
      throw new UnsupportedOperationException(IMMUTABLE_MESSAGE);
    }
    super.setEndFrame(endFrame);
  }

  @Override
  public void setPrecision(double precision) {
    if (built) {
      throw new UnsupportedOperationException(IMMUTABLE_MESSAGE);
    }
    super.setPrecision(precision);
  }

  @Override
  public void setParameter(int index, float value) {
    throw new UnsupportedOperationException(IMMUTABLE_MESSAGE);
  }

  @Override
  public void setParameterDeviation(int index, float value) {
    throw new UnsupportedOperationException(IMMUTABLE_MESSAGE);
  }

  @Override
  public void clearHasEndFrame() {
    throw new UnsupportedOperationException(IMMUTABLE_MESSAGE);
  }

  @Override
  public void clearHasId() {
    throw new UnsupportedOperationException(IMMUTABLE_MESSAGE);
  }

  @Override
  public void clearHasPrecision() {
    throw new UnsupportedOperationException(IMMUTABLE_MESSAGE);
  }

  /**
   * Gets a copy of the parameters.
   *
   * @return the parameters
   */
  @Override
  public float[] getParameters() {
    return super.getParameters().clone();
  }

  /**
   * Gets a copy of the parameter deviations.
   *
   * @return the parameter deviations
   */
  @Override
  public float[] getParameterDeviations() {
    return (super.getParameterDeviations() == null) ? null : super.getParameterDeviations().clone();
  }

  @Override
  void resizeParameters(int length) {
    throw new UnsupportedOperationException(IMMUTABLE_MESSAGE);
  }

  @Override
  void resizeParameterDeviations(int length) {
    throw new UnsupportedOperationException(IMMUTABLE_MESSAGE);
  }
}
