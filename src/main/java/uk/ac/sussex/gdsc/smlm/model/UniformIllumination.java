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

package uk.ac.sussex.gdsc.smlm.model;

/**
 * Specifies the same illumination for any position.
 */
public class UniformIllumination implements SpatialIllumination {
  private final double photons;
  private final double pulsePhotons;
  private final int pulseInterval;

  /**
   * Instantiates a new uniform illumination.
   *
   * @param photons The number of photons in a time frame
   */
  public UniformIllumination(double photons) {
    this(photons, 0, 0);
  }

  /**
   * Instantiates a new uniform illumination.
   *
   * @param photons The number of photons in a time frame
   * @param pulsePhotons The number of photons in a pulse
   * @param pulseInterval The interval between pulses (t=1 is the first pulse). Must be above 1.
   */
  public UniformIllumination(double photons, double pulsePhotons, int pulseInterval) {
    this.photons = photons;
    this.pulsePhotons = pulsePhotons;
    this.pulseInterval = pulseInterval;
  }

  /** {@inheritDoc} */
  @Override
  public double getPhotons(double[] xyz) {
    return photons;
  }

  /** {@inheritDoc} */
  @Override
  public double[] getPulsedPhotons(double[] xyz, int time) {

    if (pulseInterval > 1) {
      return new double[] {(time % pulseInterval == 1) ? pulsePhotons : 0, photons};
    }
    return new double[] {0, photons};
  }

  /** {@inheritDoc} */
  @Override
  public double getAveragePhotons() {
    if (pulseInterval > 1) {
      return photons + pulsePhotons / pulseInterval;
    }
    return photons;
  }
}
