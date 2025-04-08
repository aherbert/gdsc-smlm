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

package uk.ac.sussex.gdsc.smlm.function;

/**
 * Defines the expected variance of a signal recorded on a CCD or EM-CCD Camera. The model assumes a
 * Gaussian read noise, photon shot noise and an EM-gain noise factor.
 */
public class EmCcdCameraNoiseModel extends CameraNoiseModel {
  /**
   * Instantiates a new EM-CCD camera noise model.
   *
   * @param readNoise the read noise
   */
  public EmCcdCameraNoiseModel(final double readNoise) {
    super(readNoise);
  }

  /**
   * Instantiates a new EM-CCD camera noise model.
   *
   * @param readNoise the read noise
   * @param bias the bias
   */
  public EmCcdCameraNoiseModel(final double readNoise, final double bias) {
    super(readNoise, bias);
  }

  @Override
  public double variance(final double value) {
    return readNoise2 + Math.max(value - bias, 0.0) * 2.0;
  }

  @Override
  public boolean isEmCcd() {
    return true;
  }
}
