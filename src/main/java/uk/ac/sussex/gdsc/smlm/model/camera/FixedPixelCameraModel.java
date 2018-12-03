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
public abstract class FixedPixelCameraModel extends BaseCameraModel {
  /** The bias. */
  protected final float bias;

  /** The gain. */
  protected final float gain;

  /** The variance. */
  protected final float variance;

  /** The variance divided by the squared gain (variance./gain^2). */
  protected final float var_g2;

  /**
   * Instantiates a new fixed pixel camera model.
   *
   * @param bias the bias (in counts)
   * @param gain the gain (count/photon)
   */
  public FixedPixelCameraModel(float bias, float gain) {
    this(bias, gain, 0f);
  }

  /**
   * Instantiates a new fixed pixel camera model.
   *
   * @param bias the bias (in counts)
   * @param gain the gain (count/photon)
   */
  public FixedPixelCameraModel(double bias, double gain) {
    this(bias, gain, 0d);
  }

  /**
   * Instantiates a new fixed pixel camera model.
   *
   * @param bias the bias (in counts)
   * @param gain the gain (count/photon)
   * @param variance the variance (in counts)
   */
  public FixedPixelCameraModel(float bias, float gain, float variance) {
    checkBias(bias);
    checkGain(gain);
    checkVariance(variance);
    this.bias = bias;
    this.gain = gain;
    this.variance = variance;
    this.var_g2 = variance / (gain * gain);
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
    // to take advantage of double precision computation of var_g2.

    // Cast to float then check
    this.bias = (float) bias;
    checkBias(this.bias);
    this.gain = (float) gain;
    checkGain(this.gain);
    this.variance = (float) variance;
    checkVariance(this.variance);
    this.var_g2 = (float) (variance / (gain * gain));
  }

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
    return newArray(bounds, bias);
  }

  /** {@inheritDoc} */
  @Override
  public float[] getGain(Rectangle bounds) {
    return newArray(bounds, gain);
  }

  /** {@inheritDoc} */
  @Override
  public float[] getVariance(Rectangle bounds) {
    return newArray(bounds, variance);
  }

  /** {@inheritDoc} */
  @Override
  public float[] getNormalisedVariance(Rectangle bounds) {
    return newArray(bounds, var_g2);
  }

  @Override
  public float getBias(int x, int y) {
    return bias;
  }

  @Override
  public float getGain(int x, int y) {
    return gain;
  }

  @Override
  public float getVariance(int x, int y) {
    return variance;
  }

  @Override
  public float getNormalisedVariance(int x, int y) {
    return var_g2;
  }

  /** {@inheritDoc} */
  @Override
  public double getMeanVariance(Rectangle bounds) {
    return variance;
  }

  /** {@inheritDoc} */
  @Override
  public double getMeanNormalisedVariance(Rectangle bounds) {
    return var_g2;
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
    removeBias(data);
  }

  /** {@inheritDoc} */
  @Override
  public void removeGain(Rectangle bounds, float[] data) {
    removeGain(data);
  }

  /** {@inheritDoc} */
  @Override
  public void removeBiasAndGain(Rectangle bounds, float[] data) {
    removeBiasAndGain(data);
  }

  /** {@inheritDoc} */
  @Override
  public void applyBias(Rectangle bounds, float[] data) {
    applyBias(data);
  }

  /** {@inheritDoc} */
  @Override
  public void applyGain(Rectangle bounds, float[] data) {
    applyGain(data);
  }

  /** {@inheritDoc} */
  @Override
  public void applyGainAndBias(Rectangle bounds, float[] data) {
    applyGainAndBias(data);
  }

  /** {@inheritDoc} */
  @Override
  public void removeBias(float[] data) {
    if (data == null) {
      return;
    }
    for (int i = 0; i < data.length; i++) {
      data[i] -= bias;
    }
  }

  /** {@inheritDoc} */
  @Override
  public void removeGain(float[] data) {
    if (data == null) {
      return;
    }
    for (int i = 0; i < data.length; i++) {
      data[i] /= gain;
    }
  }

  /** {@inheritDoc} */
  @Override
  public void removeBiasAndGain(float[] data) {
    if (data == null) {
      return;
    }
    for (int i = 0; i < data.length; i++) {
      data[i] = (data[i] - bias) / gain;
    }
  }

  /** {@inheritDoc} */
  @Override
  public void applyBias(float[] data) {
    if (data == null) {
      return;
    }
    for (int i = 0; i < data.length; i++) {
      data[i] += bias;
    }
  }

  /** {@inheritDoc} */
  @Override
  public void applyGain(float[] data) {
    if (data == null) {
      return;
    }
    for (int i = 0; i < data.length; i++) {
      data[i] *= gain;
    }
  }

  /** {@inheritDoc} */
  @Override
  public void applyGainAndBias(float[] data) {
    if (data == null) {
      return;
    }
    for (int i = 0; i < data.length; i++) {
      data[i] = data[i] * gain + bias;
    }
  }

  /** {@inheritDoc} */
  @Override
  public FixedPixelCameraModel copy() {
    return clone();
  }

  /** {@inheritDoc} */
  @Override
  protected FixedPixelCameraModel clone() {
    try {
      return (FixedPixelCameraModel) super.clone();
    } catch (final CloneNotSupportedException e) {
      return null;
    }
  }
}
