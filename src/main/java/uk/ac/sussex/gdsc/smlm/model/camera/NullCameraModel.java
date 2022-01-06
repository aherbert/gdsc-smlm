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

package uk.ac.sussex.gdsc.smlm.model.camera;

import java.awt.Rectangle;

/**
 * A camera model with all pixels treated equally.
 */
public class NullCameraModel implements CameraModel {
  @Override
  public Rectangle getBounds() {
    return null;
  }

  @Override
  public void setOrigin(int x, int y) {
    // Ignore
  }

  @Override
  public CameraModel crop(Rectangle bounds, boolean resetOrigin) {
    return this;
  }

  @Override
  public boolean isPerPixelModel() {
    return false;
  }

  @Override
  public float[] getBias(Rectangle bounds) {
    return CameraModelUtils.newArray(bounds, 0);
  }

  @Override
  public float getBias(int x, int y) {
    return 0f;
  }

  @Override
  public float[] getGain(Rectangle bounds) {
    return CameraModelUtils.newArray(bounds, 1f);
  }

  @Override
  public float getGain(int x, int y) {
    return 1f;
  }

  @Override
  public float[] getVariance(Rectangle bounds) {
    return CameraModelUtils.newArray(bounds, 0);
  }

  @Override
  public float getVariance(int x, int y) {
    return 0f;
  }

  @Override
  public float[] getNormalisedVariance(Rectangle bounds) {
    return CameraModelUtils.newArray(bounds, 0);
  }

  @Override
  public float getNormalisedVariance(int x, int y) {
    return 0f;
  }

  @Override
  public double getMeanVariance(Rectangle bounds) {
    return 0d;
  }

  @Override
  public double getMeanNormalisedVariance(Rectangle bounds) {
    return 0d;
  }

  @Override
  public float[] getWeights(Rectangle bounds) {
    return CameraModelUtils.newArray(bounds, 1f);
  }

  @Override
  public float[] getNormalisedWeights(Rectangle bounds) {
    return CameraModelUtils.newArray(bounds, 1f);
  }

  @Override
  public void removeBias(Rectangle bounds, float[] data) {
    // Ignore
  }

  @Override
  public void removeBias(float[] data) {
    // Ignore
  }

  @Override
  public void removeGain(Rectangle bounds, float[] data) {
    // Ignore
  }

  @Override
  public void removeGain(float[] data) {
    // Ignore
  }

  @Override
  public void removeBiasAndGain(Rectangle bounds, float[] data) {
    // Ignore
  }

  @Override
  public void removeBiasAndGain(float[] data) {
    // Ignore
  }

  @Override
  public void applyBias(Rectangle bounds, float[] data) {
    // Ignore
  }

  @Override
  public void applyBias(float[] data) {
    // Ignore
  }

  @Override
  public void applyGain(Rectangle bounds, float[] data) {
    // Ignore
  }

  @Override
  public void applyGain(float[] data) {
    // Ignore
  }

  @Override
  public void applyGainAndBias(Rectangle bounds, float[] data) {
    // Ignore
  }

  @Override
  public void applyGainAndBias(float[] data) {
    // Ignore
  }

  @Override
  public NullCameraModel copy() {
    return this; // no state so no need to copy
  }
}
