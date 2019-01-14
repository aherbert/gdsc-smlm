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

package uk.ac.sussex.gdsc.smlm.model.camera;

import java.awt.Rectangle;

/**
 * Define the methods for manipulating camera pixel data.
 */
public interface CameraModel {
  /**
   * Gets the bounds of the camera pixel data. This could be null if the camera has infinite bounds,
   * e.g. all pixels are treated equally.
   *
   * @return the bounds
   */
  Rectangle getBounds();

  /**
   * Sets the origin. This updates the origin of the bounds.
   *
   * @param x the x
   * @param y the y
   */
  void setOrigin(int x, int y);

  /**
   * Crop the camera to the given bounds. The bounds are expected to fit within the camera bounds.
   *
   * <p>This can be used to create a more efficient representation if no data outside the bounds are
   * required.
   *
   * <p>The origin of the new model can optionally be reset to 0,0.
   *
   * <p>Note: If the bounds match the current bounds then the returned model may not be a copy.
   *
   * @param bounds the bounds
   * @param resetOrigin the reset origin flag
   * @return the camera model
   */
  CameraModel crop(Rectangle bounds, boolean resetOrigin);

  /**
   * Checks if is per pixel model. If false then all pixels are treated equally.
   *
   * @return true, if is per pixel model
   */
  boolean isPerPixelModel();

  /**
   * Gets the per-pixel camera bias (offset). The bounds are expected to fit within the camera
   * bounds.
   *
   * @param bounds the bounds
   * @return the bias
   */
  float[] getBias(Rectangle bounds);


  /**
   * Gets the per-pixel camera bias (offset). The coordinates are expected to fit within the camera
   * bounds.
   *
   * @param x the x
   * @param y the y
   * @return the bias
   */
  float getBias(int x, int y);

  /**
   * Gets the per-pixel camera gain (in count/photon). The bounds are expected to fit within the
   * camera bounds.
   *
   * @param bounds the bounds
   * @return the gain
   */
  float[] getGain(Rectangle bounds);

  /**
   * Gets the per-pixel camera gain (in count/photon). The coordinates are expected to fit within
   * the camera bounds.
   *
   * @param x the x
   * @param y the y
   * @return the gain
   */
  float getGain(int x, int y);

  /**
   * Gets the per-pixel variance. This is the variance of the pixel in camera counts. The bounds are
   * expected to fit within the camera bounds.
   *
   * @param bounds the bounds
   * @return the variance
   */
  float[] getVariance(Rectangle bounds);

  /**
   * Gets the per-pixel variance. This is the variance of the pixel in camera counts. The
   * coordinates are expected to fit within the camera bounds.
   *
   * @param x the x
   * @param y the y
   * @return the variance
   */
  float getVariance(int x, int y);

  /**
   * Gets the per-pixel normalised variance. This is the variance of the pixel in camera counts
   * divided by the squared gain, i.e. the variance in photon units. The bounds are expected to fit
   * within the camera bounds.
   *
   * @param bounds the bounds
   * @return the normalised variance
   */
  float[] getNormalisedVariance(Rectangle bounds);

  /**
   * Gets the per-pixel normalised variance. This is the variance of the pixel in camera counts
   * divided by the squared gain, i.e. the variance in photon units. The coordinates are expected to
   * fit within the camera bounds.
   *
   * @param x the x
   * @param y the y
   * @return the normalised variance
   */
  float getNormalisedVariance(int x, int y);

  /**
   * Gets the mean of the per-pixel variance. This is the variance of the pixel in camera counts.
   * The bounds are expected to fit within the camera bounds.
   *
   * @param bounds the bounds
   * @return the variance
   */
  double getMeanVariance(Rectangle bounds);

  /**
   * Gets the mean of the per-pixel normalised variance. This is the variance of the pixel in camera
   * counts divided by the squared gain, i.e. the variance in photon units. The bounds are expected
   * to fit within the camera bounds.
   *
   * @param bounds the bounds
   * @return the normalised variance
   */
  double getMeanNormalisedVariance(Rectangle bounds);

  /**
   * Gets the per-pixel weights, for example 1/variance.
   *
   * @param bounds the bounds
   * @return the weights
   */
  float[] getWeights(Rectangle bounds);

  /**
   * Gets the per-pixel normalised weights, for example 1/[normalised variance].
   *
   * @param bounds the bounds
   * @return the weights
   */
  float[] getNormalisedWeights(Rectangle bounds);

  /**
   * Remove the per-pixel camera bias (offset) from the crop of the camera data. The bounds are
   * expected to fit within the camera bounds.
   *
   * @param bounds the bounds of the data.
   * @param data the data
   */
  void removeBias(Rectangle bounds, float[] data);

  /**
   * Remove the per-pixel camera bias (offset) from the crop of the camera data. The data length is
   * expected to match the camera bounds.
   *
   * @param data the data
   */
  void removeBias(float[] data);

  /**
   * Remove the per-pixel gain from the crop of the camera data. The bounds are expected to fit
   * within the camera bounds.
   *
   * @param bounds the bounds of the data.
   * @param data the data
   */
  void removeGain(Rectangle bounds, float[] data);

  /**
   * Remove the per-pixel gain from the crop of the camera data. The data length is expected to
   * match the camera bounds.
   *
   * @param data the data
   */
  void removeGain(float[] data);

  /**
   * Remove the per-pixel camera bias (offset) and gain from the crop of the camera data. The bounds
   * are expected to fit within the camera bounds.
   *
   * @param bounds the bounds of the data.
   * @param data the data
   */
  void removeBiasAndGain(Rectangle bounds, float[] data);

  /**
   * Remove the per-pixel camera bias (offset) and gain from the crop of the camera data. The data
   * length is expected to match the camera bounds.
   *
   * @param data the data
   */
  void removeBiasAndGain(float[] data);

  /**
   * Apply the per-pixel camera bias (offset) to the crop of the camera data. The bounds are
   * expected to fit within the camera bounds.
   *
   * @param bounds the bounds of the data.
   * @param data the data
   */
  void applyBias(Rectangle bounds, float[] data);

  /**
   * Apply the per-pixel camera bias (offset) to the crop of the camera data. The data length is
   * expected to match the camera bounds.
   *
   * @param data the data
   */
  void applyBias(float[] data);

  /**
   * Apply the per-pixel gain to the crop of the camera data. The bounds are expected to fit within
   * the camera bounds.
   *
   * @param bounds the bounds of the data.
   * @param data the data
   */
  void applyGain(Rectangle bounds, float[] data);

  /**
   * Apply the per-pixel gain to the crop of the camera data. The data length is expected to match
   * the camera bounds.
   *
   * @param data the data
   */
  void applyGain(float[] data);

  /**
   * Apply the per-pixel gain and camera bias (offset) to the crop of the camera data. The bounds
   * are expected to fit within the camera bounds.
   *
   * @param bounds the bounds of the data.
   * @param data the data
   */
  void applyGainAndBias(Rectangle bounds, float[] data);

  /**
   * Apply the per-pixel gain and camera bias (offset) to the crop of the camera data. The data
   * length is expected to match the camera bounds.
   *
   * @param data the data
   */
  void applyGainAndBias(float[] data);

  /**
   * Copy this camera model. This is a deep copy of any structures.
   *
   * @return the copy
   */
  CameraModel copy();
}
