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
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;

/**
 * Define the methods for manipulating camera pixel data.
 */
public class PerPixelCameraModel implements CameraModel {
  private static final String BOUNDS_MUST_MATCH_FRAME_SIZE =
      "Bounds (width x height) must match the camera data frame size";
  private static final String FRAME_MUST_MATCH_MODEL_SIZE =
      "Camera data frame must match the camera model bounds";

  private final Rectangle cameraBounds;

  private final float[] bias;
  private final float[] gain;
  private final float[] variance;
  // This is computed when required
  private float[] varG2;

  /**
   * Instantiates a new per pixel camera model.
   *
   * @param width the width
   * @param height the height
   * @param bias the bias (in counts)
   * @param gain the gain (count/photon)
   * @param variance the variance (in counts)
   * @throws IllegalArgumentException If the data is not valid
   */
  public PerPixelCameraModel(int width, int height, float[] bias, float[] gain, float[] variance) {
    this(0, 0, width, height, bias, gain, variance);
  }

  /**
   * Instantiates a new per pixel camera model.
   *
   * @param xorigin the xorigin
   * @param yorigin the yorigin
   * @param width the width
   * @param height the height
   * @param bias the bias (in counts)
   * @param gain the gain (count/photon)
   * @param variance the variance (in counts)
   * @throws IllegalArgumentException If the data is not valid
   */
  public PerPixelCameraModel(int xorigin, int yorigin, int width, int height, float[] bias,
      float[] gain, float[] variance) {
    this(new Rectangle(xorigin, yorigin, width, height), bias, gain, variance, false, true);
  }

  /**
   * Instantiates a new per pixel camera model.
   *
   * @param bounds the bounds
   * @param bias the bias (in counts)
   * @param gain the gain (count/photon)
   * @param variance the variance (in counts)
   * @throws IllegalArgumentException If the data is not valid
   */
  public PerPixelCameraModel(Rectangle bounds, float[] bias, float[] gain, float[] variance) {
    this(bounds, bias, gain, variance, true, true);
  }

  /**
   * Instantiates a new per pixel camera model.
   *
   * @param bounds the bounds
   * @param bias the bias (in counts)
   * @param gain the gain (count/photon)
   * @param variance the variance (in counts)
   * @param cloneBounds Set to true to clone the bounds
   * @param cloneData Set to true to clone the data
   * @throws IllegalArgumentException If the data is not valid
   */
  private PerPixelCameraModel(Rectangle bounds, float[] bias, float[] gain, float[] variance,
      boolean cloneBounds, boolean cloneData) {
    if (bounds == null) {
      throw new IllegalArgumentException("Bounds must not be null");
    }
    checkBounds(bounds);
    cameraBounds = (cloneBounds) ? new Rectangle(bounds) : bounds;
    final int size = SimpleArrayUtils.check2DSize(bounds.width, bounds.height);
    checkArray(bias, size);
    checkArray(gain, size);
    checkArray(variance, size);
    if (cloneData) {
      this.bias = bias.clone();
      this.gain = gain.clone();
      this.variance = variance.clone();
    } else {
      this.bias = bias;
      this.gain = gain;
      this.variance = variance;
    }
    for (int i = 0; i < size; i++) {
      CameraModelUtils.checkBias(bias[i]);
      CameraModelUtils.checkGain(gain[i]);
      CameraModelUtils.checkVariance(variance[i]);
    }
  }

  /**
   * Check the bounds have positive coordinates and widths.
   *
   * @param bounds the bounds
   * @throws IllegalArgumentException If the data is not valid
   */
  private static void checkBounds(Rectangle bounds) {
    if (bounds.x < 0) {
      throw new IllegalArgumentException("Bounds must have positive x origin");
    }
    if (bounds.y < 0) {
      throw new IllegalArgumentException("Bounds must have positive y origin");
    }
    if (bounds.width < 0) {
      throw new IllegalArgumentException("Bounds must have positive width");
    }
    if (bounds.height < 0) {
      throw new IllegalArgumentException("Bounds must have positive height");
    }
  }

  /**
   * Instantiates a new per pixel camera model, copying all input fields.
   *
   * <p>This is an internally used copy constructor.
   *
   * @param duplicate a flag to indicate the data should be duplicated
   * @param bounds the bounds
   * @param bias the bias
   * @param gain the gain
   * @param variance the variance array
   * @param varG2 the normalised variance array
   */
  private PerPixelCameraModel(boolean duplicate, Rectangle bounds, float[] bias, float[] gain,
      float[] variance, float[] varG2) {
    if (duplicate) {
      cameraBounds = new Rectangle(bounds);
      this.bias = bias.clone();
      this.gain = gain.clone();
      this.variance = variance.clone();
      this.varG2 = (varG2 == null) ? null : varG2.clone();
    } else {
      cameraBounds = bounds;
      this.bias = bias;
      this.gain = gain;
      this.variance = variance;
      this.varG2 = varG2;
    }
  }

  /**
   * Creates a new per pixel camera model. The input arguments are not cloned.
   *
   * @param bounds the bounds
   * @param bias the bias (in counts)
   * @param gain the gain (count/photon)
   * @param variance the variance (in counts)
   * @return the per pixel camera model
   * @throws IllegalArgumentException If the data is not valid
   */
  public static PerPixelCameraModel create(Rectangle bounds, float[] bias, float[] gain,
      float[] variance) {
    return new PerPixelCameraModel(bounds, bias, gain, variance, false, false);
  }

  /**
   * Check array.
   *
   * @param array the array
   * @param size the size
   */
  private static void checkArray(float[] array, int size) {
    if (array == null || array.length != size) {
      throw new IllegalArgumentException("Input array must match the size of the input bounds");
    }
  }

  @Override
  public Rectangle getBounds() {
    return new Rectangle(cameraBounds);
  }

  @Override
  public void setOrigin(int x, int y) {
    cameraBounds.x = x;
    cameraBounds.y = y;
  }

  /**
   * Gets the x origin.
   *
   * @return the x origin
   */
  public int getXOrigin() {
    return cameraBounds.x;
  }

  /**
   * Gets the y origin.
   *
   * @return the y origin
   */
  public int getYOrigin() {
    return cameraBounds.y;
  }

  /**
   * Gets the width.
   *
   * @return the width
   */
  public int getWidth() {
    return cameraBounds.width;
  }

  /**
   * Gets the height.
   *
   * @return the height
   */
  public int getHeight() {
    return cameraBounds.height;
  }

  /**
   * Gets a copy of the bias for the current bounds.
   *
   * @return the bias
   */
  public float[] getBias() {
    return bias.clone();
  }

  @Override
  public float[] getBias(Rectangle bounds) {
    return getData(bounds, bias);
  }

  @Override
  public float getBias(int x, int y) {
    return getData(x, y, bias);
  }

  /**
   * Gets a copy of the gain for the current bounds.
   *
   * @return the gain
   */
  public float[] getGain() {
    return gain.clone();
  }

  @Override
  public float[] getGain(Rectangle bounds) {
    return getData(bounds, gain);
  }

  @Override
  public float getGain(int x, int y) {
    return getData(x, y, gain);
  }

  /**
   * Gets a copy of the variance for the current bounds.
   *
   * @return the variance
   */
  public float[] getVariance() {
    return variance.clone();
  }

  @Override
  public float[] getVariance(Rectangle bounds) {
    return getData(bounds, variance);
  }

  @Override
  public float getVariance(int x, int y) {
    return getData(x, y, variance);
  }

  /**
   * Gets a copy of the normalised variance for the current bounds.
   *
   * @return the normalised variance
   */
  public float[] getNormalisedVariance() {
    return getNormalisedVarianceInternal().clone();
  }

  @Override
  public float[] getNormalisedVariance(Rectangle bounds) {
    return getData(bounds, getNormalisedVarianceInternal());
  }

  @Override
  public float getNormalisedVariance(int x, int y) {
    return getData(x, y, getNormalisedVarianceInternal());
  }

  @Override
  public boolean isPerPixelModel() {
    return true;
  }

  /**
   * Initialise the model. This allows caching of precomputed values but it is not required to call
   * this method before using the model.
   */
  public void initialise() {
    getNormalisedVarianceInternal();
  }

  private float[] getNormalisedVarianceInternal() {
    if (varG2 == null) {
      createNormalisedVariance();
    }
    return varG2;
  }

  private synchronized void createNormalisedVariance() {
    if (varG2 == null) {
      final int size = variance.length;
      varG2 = new float[size];
      for (int i = 0; i < size; i++) {
        varG2[i] = variance[i] / (gain[i] * gain[i]);
      }
    }
  }

  @Override
  public double getMeanVariance(Rectangle bounds) {
    return getMean(bounds, variance);
  }

  @Override
  public double getMeanNormalisedVariance(Rectangle bounds) {
    return getMean(bounds, getNormalisedVarianceInternal());
  }

  /**
   * Return the weights as 1/variance.
   */
  @Override
  public float[] getWeights(Rectangle bounds) {
    return CameraModelUtils.toWeights(getVariance(bounds));
  }

  /**
   * Return the weights as 1/[normalised variance].
   */
  @Override
  public float[] getNormalisedWeights(Rectangle bounds) {
    return CameraModelUtils.toWeights(getNormalisedVariance(bounds));
  }

  /**
   * Gets the data from the values using the intersection of the bounds.
   *
   * @param bounds the bounds
   * @param values the values
   * @return the data
   */
  private float[] getData(Rectangle bounds, float[] values) {
    final Rectangle intersection = getIntersection(bounds);
    return getData(values, intersection, true);
  }

  /**
   * Crop the data from the per-pixel data using the given bounds.
   *
   * @param pixels the pixels
   * @param bounds the bounds
   * @param copy Set to true to copy the values (if the bounds cover all the pixel data)
   * @return the data
   */
  private float[] getData(float[] pixels, final Rectangle bounds, boolean copy) {
    if (bounds.x != 0 || bounds.y != 0 || bounds.width != cameraBounds.width
        || bounds.height != cameraBounds.height) {
      final float[] pixels2 = new float[bounds.width * bounds.height];
      final int width = cameraBounds.width;
      for (int ys = 0, offset1 = 0; ys < bounds.height; ys++) {
        for (int xs = 0, offset2 = (ys + bounds.y) * width + bounds.x; xs < bounds.width; xs++) {
          pixels2[offset1++] = pixels[offset2++];
        }
      }
      return pixels2;
    }
    return (copy) ? pixels.clone() : pixels;
  }

  /**
   * Gets the data value.
   *
   * @param x the x
   * @param y the y
   * @param data the data
   * @return the data value
   * @throws IllegalArgumentException If the coordinates are not inside the bounds
   */
  private float getData(int x, int y, float[] data) {
    x -= cameraBounds.x;
    y -= cameraBounds.y;
    if (x < 0 || y < 0 || x >= cameraBounds.width || y >= cameraBounds.height) {
      throw new IllegalArgumentException("Coordinates must be within the camera bounds");
    }
    return data[y * cameraBounds.width + x];
  }

  /**
   * Gets the mean of the data from the values using the intersection of the bounds.
   *
   * @param bounds the bounds
   * @param values the values
   * @return the mean of the data
   */
  private double getMean(Rectangle bounds, float[] values) {
    final Rectangle intersection = getIntersection(bounds);
    return getMean(values, intersection);
  }

  /**
   * Get the mean of the data from the per-pixel data using the given bounds.
   *
   * @param pixels the pixels
   * @param bounds the bounds
   * @return the mean of the data
   */
  private double getMean(float[] pixels, final Rectangle bounds) {
    double sum = 0;
    if (bounds.x != 0 || bounds.y != 0 || bounds.width != cameraBounds.width
        || bounds.height != cameraBounds.height) {
      final int width = cameraBounds.width;
      for (int ys = 0; ys < bounds.height; ys++) {
        for (int xs = 0, offset2 = (ys + bounds.y) * width + bounds.x; xs < bounds.width; xs++) {
          sum += pixels[offset2++];
        }
      }
      return sum / (bounds.height * bounds.width);
    }
    for (int i = pixels.length; i-- > 0;) {
      sum += pixels[i];
    }
    return sum / pixels.length;
  }

  /**
   * Gets the intersection between the target bounds and the camera bounds. The result is offset by
   * the camera bounds origin so can be used to crop the camera per-pixel data. I.e. if the camera
   * bounds are [3,4,10,10] and the bounds are [3,4,5,5] the result will be [0,0,5,5].
   *
   * @param bounds the bounds
   * @return the intersection
   */
  private Rectangle getIntersection(Rectangle bounds) {
    if (bounds == null) {
      throw new IllegalArgumentException("Bounds are null");
    }
    if (equalBounds(bounds)) {
      return new Rectangle(getWidth(), getHeight()); // Offset to the origin
    }
    checkBounds(bounds);
    // We avoid underflow since we have checked the bounds are positive integers
    final int minx = bounds.x - cameraBounds.x;
    final int miny = bounds.y - cameraBounds.y;
    if (minx < 0 || miny < 0) {
      throw new IllegalArgumentException("Bounds must be within the camera bounds");
    }
    // Avoid overflow using a long result
    final long maxx = (long) minx + bounds.width;
    final long maxy = (long) miny + bounds.height;
    if (maxx > cameraBounds.width || maxy > cameraBounds.height) {
      throw new IllegalArgumentException("Bounds must be within the camera bounds");
    }
    return new Rectangle(minx, miny, bounds.width, bounds.height);
  }

  /**
   * Check if the bounds are equal to the camera bounds.
   *
   * @param bounds the bounds
   * @return true, if successful
   */
  private boolean equalBounds(Rectangle bounds) {
    //@formatter:off
    return cameraBounds.x == bounds.x &&
         cameraBounds.y == bounds.y &&
         cameraBounds.width == bounds.width &&
         cameraBounds.height == bounds.height;
    //@formatter:on
  }

  @Override
  public void removeBias(Rectangle bounds, float[] data) {
    if (data == null) {
      return;
    }
    final Rectangle intersection = getIntersection(bounds);
    final float[] biasData = getData(this.bias, intersection, false);
    if (data.length != biasData.length) {
      throw new IllegalArgumentException(BOUNDS_MUST_MATCH_FRAME_SIZE);
    }
    for (int i = 0; i < data.length; i++) {
      data[i] -= biasData[i];
    }
  }

  @Override
  public void removeBias(float[] data) {
    if (data == null) {
      return;
    }
    if (data.length != bias.length) {
      throw new IllegalArgumentException(FRAME_MUST_MATCH_MODEL_SIZE);
    }
    for (int i = 0; i < data.length; i++) {
      data[i] -= bias[i];
    }
  }

  @Override
  public void removeGain(Rectangle bounds, float[] data) {
    if (data == null) {
      return;
    }
    final Rectangle intersection = getIntersection(bounds);
    final float[] gainData = getData(this.gain, intersection, false);
    if (data.length != gainData.length) {
      throw new IllegalArgumentException(BOUNDS_MUST_MATCH_FRAME_SIZE);
    }
    for (int i = 0; i < data.length; i++) {
      data[i] /= gainData[i];
    }
  }

  @Override
  public void removeGain(float[] data) {
    if (data == null) {
      return;
    }
    if (data.length != gain.length) {
      throw new IllegalArgumentException(FRAME_MUST_MATCH_MODEL_SIZE);
    }
    for (int i = 0; i < data.length; i++) {
      data[i] /= gain[i];
    }
  }

  @Override
  public void removeBiasAndGain(Rectangle bounds, float[] data) {
    if (data == null) {
      return;
    }
    final Rectangle intersection = getIntersection(bounds);
    final float[] biasData = getData(this.bias, intersection, false);
    if (data.length != biasData.length) {
      throw new IllegalArgumentException(BOUNDS_MUST_MATCH_FRAME_SIZE);
    }
    final float[] gainData = getData(this.gain, intersection, false);
    for (int i = 0; i < data.length; i++) {
      data[i] = (data[i] - biasData[i]) / gainData[i];
    }
  }

  @Override
  public void removeBiasAndGain(float[] data) {
    if (data == null) {
      return;
    }
    if (data.length != bias.length) {
      throw new IllegalArgumentException(FRAME_MUST_MATCH_MODEL_SIZE);
    }
    for (int i = 0; i < data.length; i++) {
      data[i] = (data[i] - bias[i]) / gain[i];
    }
  }

  @Override
  public void applyBias(Rectangle bounds, float[] data) {
    if (data == null) {
      return;
    }
    final Rectangle intersection = getIntersection(bounds);
    final float[] biasData = getData(this.bias, intersection, false);
    if (data.length != biasData.length) {
      throw new IllegalArgumentException(BOUNDS_MUST_MATCH_FRAME_SIZE);
    }
    for (int i = 0; i < data.length; i++) {
      data[i] += biasData[i];
    }
  }

  @Override
  public void applyBias(float[] data) {
    if (data == null) {
      return;
    }
    if (data.length != bias.length) {
      throw new IllegalArgumentException(FRAME_MUST_MATCH_MODEL_SIZE);
    }
    for (int i = 0; i < data.length; i++) {
      data[i] += bias[i];
    }
  }

  @Override
  public void applyGain(Rectangle bounds, float[] data) {
    if (data == null) {
      return;
    }
    final Rectangle intersection = getIntersection(bounds);
    final float[] gainData = getData(this.gain, intersection, false);
    if (data.length != gainData.length) {
      throw new IllegalArgumentException(BOUNDS_MUST_MATCH_FRAME_SIZE);
    }
    for (int i = 0; i < data.length; i++) {
      data[i] *= gainData[i];
    }
  }

  @Override
  public void applyGain(float[] data) {
    if (data == null) {
      return;
    }
    if (data.length != gain.length) {
      throw new IllegalArgumentException(FRAME_MUST_MATCH_MODEL_SIZE);
    }
    for (int i = 0; i < data.length; i++) {
      data[i] *= gain[i];
    }
  }

  @Override
  public void applyGainAndBias(Rectangle bounds, float[] data) {
    if (data == null) {
      return;
    }
    final Rectangle intersection = getIntersection(bounds);
    final float[] biasData = getData(this.bias, intersection, false);
    if (data.length != biasData.length) {
      throw new IllegalArgumentException(BOUNDS_MUST_MATCH_FRAME_SIZE);
    }
    final float[] gainData = getData(this.gain, intersection, false);
    for (int i = 0; i < data.length; i++) {
      data[i] = data[i] * gainData[i] + biasData[i];
    }
  }

  @Override
  public void applyGainAndBias(float[] data) {
    if (data == null) {
      return;
    }
    if (data.length != bias.length) {
      throw new IllegalArgumentException(FRAME_MUST_MATCH_MODEL_SIZE);
    }
    for (int i = 0; i < data.length; i++) {
      data[i] = data[i] * gain[i] + bias[i];
    }
  }

  @Override
  public CameraModel crop(Rectangle bounds, boolean resetOrigin) {
    if (bounds == null) {
      throw new IllegalArgumentException("Bounds are null");
    }
    if (equalBounds(bounds)) {
      if (resetOrigin) {
        final PerPixelCameraModel model = copy();
        model.setOrigin(0, 0);
        return model;
      }
      return this;
    }
    final Rectangle intersection = getIntersection(bounds);
    final float[] biasC = getData(this.bias, intersection, true);
    final float[] gainC = getData(this.gain, intersection, true);
    final float[] varianceC = getData(this.variance, intersection, true);
    final float[] varG2C = (this.varG2 == null) ? null : getData(this.varG2, intersection, true);
    return new PerPixelCameraModel(false, bounds, biasC, gainC, varianceC, varG2C);
  }

  @Override
  public PerPixelCameraModel copy() {
    // Deep copy
    return new PerPixelCameraModel(true, cameraBounds, bias, gain, variance, varG2);
  }
}
