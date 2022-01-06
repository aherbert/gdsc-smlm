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
public abstract class FixedPixelCameraModel implements CameraModel {
  /** The bias. */
  protected final float bias;

  /** The gain. */
  protected final float gain;

  /** The variance. */
  protected final float variance;

  /** The variance divided by the squared gain (variance./gain^2). */
  protected final float varG2;

  /**
   * Instantiates a new fixed pixel camera model.
   *
   * @param bias the bias (in counts)
   * @param gain the gain (count/photon)
   */
  public FixedPixelCameraModel(float bias, float gain) {
    this(bias, gain, 0F);
  }

  /**
   * Instantiates a new fixed pixel camera model.
   *
   * @param bias the bias (in counts)
   * @param gain the gain (count/photon)
   */
  public FixedPixelCameraModel(double bias, double gain) {
    this(bias, gain, 0D);
  }

  /**
   * Instantiates a new fixed pixel camera model.
   *
   * @param bias the bias (in counts)
   * @param gain the gain (count/photon)
   * @param variance the variance (in counts)
   */
  public FixedPixelCameraModel(float bias, float gain, float variance) {
    CameraModelUtils.checkBias(bias);
    CameraModelUtils.checkGain(gain);
    CameraModelUtils.checkVariance(variance);
    this.bias = bias;
    this.gain = gain;
    this.variance = variance;
    this.varG2 = variance / (gain * gain);
  }

  /**
   * Instantiates a new fixed pixel camera model.
   *
   * @param bias the bias (in counts)
   * @param gain the gain (count/photon)
   * @param variance the variance (in counts)
   */
  public FixedPixelCameraModel(double bias, double gain, double variance) {
    // Re-implement the constructor (rather than chaining)
    // to take advantage of double precision computation of varG2.

    // Cast to float then check
    this.bias = (float) bias;
    CameraModelUtils.checkBias(this.bias);
    this.gain = (float) gain;
    CameraModelUtils.checkGain(this.gain);
    this.variance = (float) variance;
    CameraModelUtils.checkVariance(this.variance);
    this.varG2 = (float) (variance / (gain * gain));
  }

  /**
   * Copy constructor.
   *
   * @param source the source
   */
  protected FixedPixelCameraModel(FixedPixelCameraModel source) {
    this.bias = source.bias;
    this.gain = source.gain;
    this.variance = source.variance;
    this.varG2 = source.varG2;
  }

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
    return CameraModelUtils.newArray(bounds, bias);
  }

  @Override
  public float getBias(int x, int y) {
    return bias;
  }

  @Override
  public float[] getGain(Rectangle bounds) {
    return CameraModelUtils.newArray(bounds, gain);
  }

  @Override
  public float getGain(int x, int y) {
    return gain;
  }

  @Override
  public float[] getVariance(Rectangle bounds) {
    return CameraModelUtils.newArray(bounds, variance);
  }

  @Override
  public float getVariance(int x, int y) {
    return variance;
  }

  @Override
  public float[] getNormalisedVariance(Rectangle bounds) {
    return CameraModelUtils.newArray(bounds, varG2);
  }

  @Override
  public float getNormalisedVariance(int x, int y) {
    return varG2;
  }

  @Override
  public double getMeanVariance(Rectangle bounds) {
    return variance;
  }

  @Override
  public double getMeanNormalisedVariance(Rectangle bounds) {
    return varG2;
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
    removeBias(data);
  }

  @Override
  public void removeBias(float[] data) {
    if (data == null) {
      return;
    }
    for (int i = 0; i < data.length; i++) {
      data[i] -= bias;
    }
  }

  @Override
  public void removeGain(Rectangle bounds, float[] data) {
    removeGain(data);
  }

  @Override
  public void removeGain(float[] data) {
    if (data == null) {
      return;
    }
    for (int i = 0; i < data.length; i++) {
      data[i] /= gain;
    }
  }

  @Override
  public void removeBiasAndGain(Rectangle bounds, float[] data) {
    removeBiasAndGain(data);
  }


  @Override
  public void removeBiasAndGain(float[] data) {
    if (data == null) {
      return;
    }
    for (int i = 0; i < data.length; i++) {
      data[i] = (data[i] - bias) / gain;
    }
  }

  @Override
  public void applyBias(Rectangle bounds, float[] data) {
    applyBias(data);
  }

  @Override
  public void applyBias(float[] data) {
    if (data == null) {
      return;
    }
    for (int i = 0; i < data.length; i++) {
      data[i] += bias;
    }
  }

  @Override
  public void applyGain(Rectangle bounds, float[] data) {
    applyGain(data);
  }

  @Override
  public void applyGain(float[] data) {
    if (data == null) {
      return;
    }
    for (int i = 0; i < data.length; i++) {
      data[i] *= gain;
    }
  }

  @Override
  public void applyGainAndBias(Rectangle bounds, float[] data) {
    applyGainAndBias(data);
  }

  @Override
  public void applyGainAndBias(float[] data) {
    if (data == null) {
      return;
    }
    for (int i = 0; i < data.length; i++) {
      data[i] = data[i] * gain + bias;
    }
  }
}
