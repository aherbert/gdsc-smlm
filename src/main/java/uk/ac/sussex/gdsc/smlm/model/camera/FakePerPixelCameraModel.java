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

/**
 * A per-pixel camera model with all pixels treated equally. Note that this concept is invalid since
 * this model reports itself as a per-pixel model even though all pixels are treated equally. This
 * allows testing algorithms that require a per-pixel model with a fixed-pixel model, e.g. for
 * performance comparison.
 *
 * @author Alex Herbert
 */
public class FakePerPixelCameraModel extends FixedPixelCameraModel {
  /**
   * Instantiates a new fake per pixel camera model.
   *
   * @param bias the bias (in counts)
   * @param gain the gain (count/photon)
   */
  public FakePerPixelCameraModel(float bias, float gain) {
    super(bias, gain);
  }

  /**
   * Instantiates a new fake per pixel camera model.
   *
   * @param bias the bias (in counts)
   * @param gain the gain (count/photon)
   * @param variance the variance (in counts)
   */
  public FakePerPixelCameraModel(float bias, float gain, float variance) {
    super(bias, gain, variance);
  }

  /**
   * Instantiates a new fake per pixel camera model.
   *
   * @param bias the bias (in counts)
   * @param gain the gain (count/photon)
   */
  public FakePerPixelCameraModel(double bias, double gain) {
    super(bias, gain);
  }

  /**
   * Instantiates a new fake per pixel camera model.
   *
   * @param bias the bias (in counts)
   * @param gain the gain (count/photon)
   * @param variance the variance (in counts)
   */
  public FakePerPixelCameraModel(double bias, double gain, double variance) {
    super(bias, gain, variance);
  }

  /** {@inheritDoc} */
  @Override
  public boolean isPerPixelModel() {
    return true;
  }
}
