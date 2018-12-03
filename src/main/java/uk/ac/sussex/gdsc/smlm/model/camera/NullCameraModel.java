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
package uk.ac.sussex.gdsc.smlm.model.camera;

import java.awt.Rectangle;

/**
 * A camera model with all pixels treated equally.
 *
 * @author Alex Herbert
 */
public class NullCameraModel extends BaseCameraModel {
  /** {@inheritDoc} */
  @Override
  public Rectangle getBounds() {
    return null;
  }

  /** {@inheritDoc} */
  @Override
  public void setOrigin(int x, int y) {
    // Ignore
  }

  /** {@inheritDoc} */
  @Override
  public CameraModel crop(Rectangle bounds, boolean resetOrigin) {
    return this;
  }

  /** {@inheritDoc} */
  @Override
  public boolean isPerPixelModel() {
    return false;
  }

  /** {@inheritDoc} */
  @Override
  public float[] getBias(Rectangle bounds) {
    return newArray(bounds, 0);
  }

  /** {@inheritDoc} */
  @Override
  public float[] getGain(Rectangle bounds) {
    return newArray(bounds, 1f);
  }

  /** {@inheritDoc} */
  @Override
  public float[] getVariance(Rectangle bounds) {
    return newArray(bounds, 0);
  }

  /** {@inheritDoc} */
  @Override
  public float[] getNormalisedVariance(Rectangle bounds) {
    return newArray(bounds, 0);
  }

  @Override
  public float getBias(int x, int y) {
    return 0f;
  }

  @Override
  public float getGain(int x, int y) {
    return 1f;
  }

  @Override
  public float getVariance(int x, int y) {
    return 0f;
  }

  @Override
  public float getNormalisedVariance(int x, int y) {
    return 0f;
  }

  /** {@inheritDoc} */
  @Override
  public double getMeanVariance(Rectangle bounds) {
    return 0d;
  }

  /** {@inheritDoc} */
  @Override
  public double getMeanNormalisedVariance(Rectangle bounds) {
    return 0d;
  }

  /** {@inheritDoc} */
  @Override
  public float[] getWeights(Rectangle bounds) {
    return newArray(bounds, 1f);
  }

  /** {@inheritDoc} */
  @Override
  public float[] getNormalisedWeights(Rectangle bounds) {
    return newArray(bounds, 1f);
  }

  /** {@inheritDoc} */
  @Override
  public void removeBias(Rectangle bounds, float[] data) {
    // Ignore
  }

  /** {@inheritDoc} */
  @Override
  public void removeGain(Rectangle bounds, float[] data) {
    // Ignore
  }

  /** {@inheritDoc} */
  @Override
  public void removeBiasAndGain(Rectangle bounds, float[] data) {
    // Ignore
  }

  /** {@inheritDoc} */
  @Override
  public void applyBias(Rectangle bounds, float[] data) {
    // Ignore
  }

  /** {@inheritDoc} */
  @Override
  public void applyGain(Rectangle bounds, float[] data) {
    // Ignore
  }

  /** {@inheritDoc} */
  @Override
  public void applyGainAndBias(Rectangle bounds, float[] data) {
    // Ignore
  }

  /** {@inheritDoc} */
  @Override
  public void removeBias(float[] data) {
    // Ignore
  }

  /** {@inheritDoc} */
  @Override
  public void removeGain(float[] data) {
    // Ignore
  }

  /** {@inheritDoc} */
  @Override
  public void removeBiasAndGain(float[] data) {
    // Ignore
  }

  /** {@inheritDoc} */
  @Override
  public void applyBias(float[] data) {
    // Ignore
  }

  /** {@inheritDoc} */
  @Override
  public void applyGain(float[] data) {
    // Ignore
  }

  /** {@inheritDoc} */
  @Override
  public void applyGainAndBias(float[] data) {
    // Ignore
  }

  /** {@inheritDoc} */
  @Override
  public NullCameraModel copy() {
    return this; // no state so no need to clone()
  }

  /** {@inheritDoc} */
  @Override
  protected NullCameraModel clone() {
    try {
      return (NullCameraModel) super.clone();
    } catch (final CloneNotSupportedException ex) {
      return null;
    }
  }
}
