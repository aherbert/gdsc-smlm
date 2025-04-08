/*-
 * #%L
 * Genome Damage and Stability Centre SMLM Package
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2025 Alex Herbert
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
 * An EM-CCD camera model with all pixels treated equally.
 */
public class EmCcdCameraModel extends FixedPixelCameraModel {
  /**
   * Instantiates a new EM-CCD camera model.
   *
   * @param bias the bias (in counts)
   * @param gain the total gain (count/photon)
   */
  public EmCcdCameraModel(float bias, float gain) {
    this(bias, gain, 0f);
  }

  /**
   * Instantiates a new EM-CCD camera model.
   *
   * @param bias the bias (in counts)
   * @param gain the total gain (count/photon)
   */
  public EmCcdCameraModel(double bias, double gain) {
    this(bias, gain, 0d);
  }

  /**
   * Instantiates a new EM-CCD camera model.
   *
   * @param bias the bias (in counts)
   * @param gain the total gain (count/photon)
   * @param variance the variance (in counts)
   */
  public EmCcdCameraModel(float bias, float gain, float variance) {
    super(bias, gain, variance);
  }

  /**
   * Instantiates a new EM-CCD camera model.
   *
   * @param bias the bias (in counts)
   * @param gain the total gain (count/photon)
   * @param variance the variance (in counts)
   */
  public EmCcdCameraModel(double bias, double gain, double variance) {
    super(bias, gain, variance);
  }

  /**
   * Copy constructor.
   *
   * @param source the source
   */
  protected EmCcdCameraModel(EmCcdCameraModel source) {
    super(source);
  }

  /**
   * {@inheritDoc}
   *
   * <p>Note: This is an EM-CCD camera model. The normalised variance represents the effective read
   * noise in incident photons (i.e. before EM-gain). This can be combined with the expected shot
   * variance of a Poisson distribution (mean) scaled by the EM-amplification noise factor (2) to
   * obtain the total variance in photon units:
   *
   * <pre>
   * Total variance (photons) = [Poisson mean] * 2 + [normalised variance]
   * </pre>
   *
   * <p>This value multiplied by the [gain]^2 is the variance in counts.
   */
  @Override
  public float[] getNormalisedVariance(Rectangle bounds) {
    return super.getNormalisedVariance(bounds);
  }

  /**
   * {@inheritDoc}
   *
   * <p>Note: This is an EM-CCD camera model. The normalised variance represents the effective read
   * noise in incident photons (i.e. before EM-gain). This can be combined with the expected shot
   * variance of a Poisson distribution (mean) scaled by the EM-amplification noise factor (2) to
   * obtain the total variance in photon units:
   *
   * <pre>
   * Total variance (photons) = [Poisson mean] * 2 + [normalised variance]
   * </pre>
   *
   * <p>This value multiplied by the [gain]^2 is the variance in counts.
   */
  @Override
  public double getMeanNormalisedVariance(Rectangle bounds) {
    return super.getMeanNormalisedVariance(bounds);
  }

  @Override
  public EmCcdCameraModel copy() {
    return new EmCcdCameraModel(this);
  }
}
