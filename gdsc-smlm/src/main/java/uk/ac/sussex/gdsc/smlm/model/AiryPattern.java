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

package uk.ac.sussex.gdsc.smlm.model;

import uk.ac.sussex.gdsc.smlm.function.Bessel;
import uk.ac.sussex.gdsc.smlm.utils.StdMath;

/**
 * Simple class to calculate AiryPattern.
 */
public final class AiryPattern {
  /**
   * Convert Airy radius to Gaussian approximation.
   *
   * <p>See Abraham, et al (2009) Opt Express 17: 23352.
   */
  public static final double FACTOR = 1.323;

  /** No instances. */
  private AiryPattern() {}

  /**
   * Calculate the intensity of the AiryPattern at distance x from the centre.
   *
   * @param x the x
   * @return The intensity
   */
  public static double intensity(final double x) {
    if (x == 0) {
      return 1;
    }
    // I(x) = I0 * ( 2*J1(x) / x )^2
    // I0 is the maximum intensity of the pattern at the centre
    // Assuming I0 = 1:
    final double y = Bessel.j1(x) / x;
    return 4.0 * y * y;
  }

  /**
   * Calculate the intensity of the AiryPattern at distance x from the centre using a Gaussian
   * approximation.
   *
   * @param x the x
   * @return The intensity
   */
  public static double intensityGaussian(double x) {
    if (x == 0) {
      return 1;
    }
    x /= FACTOR;
    return StdMath.exp(-0.5 * (x * x));
  }

  /**
   * Calculate the total power of the AiryPattern at distance x from the centre.
   *
   * <p>Appears to be numerically unstable at {@code x<<1}.
   *
   * @param x the x
   * @return The total power
   */
  public static double power(final double x) {
    if (x == 0) {
      return 0;
    }
    final double j0 = Bessel.j0(x);
    final double j1 = Bessel.j1(x);
    return 1.0 - j0 * j0 - j1 * j1;
  }
}
