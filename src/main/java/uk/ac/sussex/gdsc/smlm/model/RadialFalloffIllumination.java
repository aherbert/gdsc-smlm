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

package uk.ac.sussex.gdsc.smlm.model;

/**
 * Specifies the illumination with radial fall-off to simulate a wide field confocal microscope.
 *
 * <p>The intensity falls-off from the centre to the edge proportional to the square of the distance
 * to the centre. There is uniform intensity in the z-axis.
 */
public class RadialFalloffIllumination implements SpatialIllumination {
  private final double photons;
  private final double radius2;
  private final double pulsePhotons;
  private final int pulseInterval;

  /**
   * Assume the default refractive index.
   *
   * @param photons The photons in the focal plane
   * @param radius The radius where intensity is half that in the centre
   */
  public RadialFalloffIllumination(double photons, double radius) {
    this(photons, radius, 0, 0);
  }

  /**
   * Assume the default refractive index.
   *
   * @param photons The photons in the focal plane
   * @param radius The radius where intensity is half that in the centre
   * @param pulsePhotons The number of photons in a pulse
   * @param pulseInterval The interval between pulses (t=1 is the first pulse). Must be above 1.
   */
  public RadialFalloffIllumination(double photons, double radius, double pulsePhotons,
      int pulseInterval) {
    this.photons = photons;
    this.radius2 = radius * radius;
    this.pulsePhotons = pulsePhotons;
    this.pulseInterval = pulseInterval;
  }

  /** {@inheritDoc} */
  @Override
  public double getPhotons(double[] xyz) {
    return photons * getIntensity(xyz);
  }

  /**
   * Get the intensity of the illumination given the distance from the centre.
   *
   * @param xyz the xyz
   * @return the intensity
   */
  private double getIntensity(double[] xyz) {
    final double d2 = xyz[0] * xyz[0] + xyz[1] * xyz[1];
    return 1.0 / (1.0 + (d2 / radius2));
  }

  /** {@inheritDoc} */
  @Override
  public double[] getPulsedPhotons(double[] xyz, int t) {
    final double intensity = getIntensity(xyz);
    if (pulseInterval > 1) {
      return new double[] {(t % pulseInterval == 1) ? pulsePhotons * intensity : 0,
          photons * intensity};
    }
    return new double[] {0, photons * intensity};
  }

  /**
   * {@inheritDoc}
   *
   * @return The average intensity from the centre to the radius (specified in the constructor)
   */
  @Override
  public double getAveragePhotons() {
    // This should be the integral of the getIntensity() score from r = 0 ..
    // R divided by R
    // http://en.wikipedia.org/wiki/List_of_integrals_of_rational_functions#Integrands_of_the_form_xm_.2F_.28a_x2_.2B_b_x_.2B_c.29n
    //
    // f(x) = 1 / (ax^2 + bx + c)
    // F(x) = 2/sqrt(4ac - b^2) arctan((2ax + b) / sqrt(4ac - b^2)) + C (for
    // 4ac - b^2 > 0)
    //
    // For a = 1 / R^2, b = 0, c = 1
    // Simplifies:
    // F(x) = 2/sqrt(4/R^2) arctan((2x/R^2) / sqrt(4/R^2))
    // = 2/(2/R) arctan((2x/R^2) / (2/R))
    // = R arctan(x/R)
    // Solving for r=R - r=0:
    // Sum = R arctan(1) - R arctan(0)
    // = R arctan(1)
    // => Average = R arctan(1) / R = arctan(1) = Math.Pi / 4;
    if (pulseInterval > 1) {
      return (photons + pulsePhotons / pulseInterval) * (Math.PI / 4);
    }
    return photons * (Math.PI / 4);
  }
}
