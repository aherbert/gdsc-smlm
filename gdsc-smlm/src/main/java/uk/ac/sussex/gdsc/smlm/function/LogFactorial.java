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

package uk.ac.sussex.gdsc.smlm.function;

import org.apache.commons.math3.special.Gamma;
import uk.ac.sussex.gdsc.core.utils.ValidationUtils;

/**
 * Compute the natural logarithm of the factorial function ({@code ln n!}).
 */
public final class LogFactorial {
  /** The maximum representable factorial. */
  private static final int MAX_FACTORIAL = 170;
  /** ln(sqrt(2 pi)). Computed to 25-digits precision. */
  private static final double LN_SQRT_TWO_PI = 0.9189385332046727417803297;
  /** Precomputed factor LANCZOS_G + 1.5. This has a value of 799 / 128. */
  private static final double LANCZOS_G_PLUS_1_5 = Gamma.LANCZOS_G + 1.5;

  /** No instances. */
  private LogFactorial() {}

  /**
   * Get the natural logarithm of the factorial of n.
   *
   * @param n argument (must be positive)
   * @return ln n!
   * @throws ArrayIndexOutOfBoundsException if n is negative
   */
  public static double value(int n) {
    if (n <= MAX_FACTORIAL) {
      return Math.log(Factorial.uncheckedValue(n));
    }
    return logGamma1p(n);
  }

  /**
   * Get the natural logarithm of the factorial of n as a real number.
   *
   * <p>This is computed using {@code ln(gamma(1+n))}.
   *
   * @param n argument (must be positive)
   * @return ln n!
   * @throws IllegalArgumentException if n is negative
   */
  public static double value(double n) {
    ValidationUtils.checkPositive(n);
    if (n <= MAX_FACTORIAL && Math.rint(n) == n) {
      return Math.log(Factorial.uncheckedValue((int) n));
    }
    return logGamma1p(n);
  }

  /**
   * Compute {@code log(Gamma(1+x))}.
   *
   * <p>Adapted from {@code org.apache.commons.math3.special.Gamma} logGamma. Removed support for
   * NaN and negative x. Updated the code to use x or (x+1) where appropriate.
   *
   * <p>For {@code x <= 7}, the implementation is based on the double precision implementation in
   * the <em>NSWC Library of Mathematics Subroutines</em>, {@code DGAMLN}; otherwise the
   * implementation is based on </p>
   *
   * <ul>
   *
   * <li><a href="http://mathworld.wolfram.com/GammaFunction.html">Gamma Function</a>, equation
   * (28).</li>
   *
   * <li><a href="https://mathworld.wolfram.com/LanczosApproximation.html"> Lanczos
   * Approximation</a>, equations (1) through (5).</li>
   *
   * <li><a href="http://www.numericana.com/answer/info/godfrey.htm">Paul Godfrey, A note on the
   * computation of the convergent Lanczos complex Gamma approximation</a></li>
   *
   * </ul>
   *
   * @param x Argument.
   * @return the value of {@code log(Gamma(1+x))}
   */
  private static double logGamma1p(double x) {
    if (x <= 1.5) {
      // No check for x == 0 as this is handled as a representable int
      return -Math.log1p(Gamma.invGamma1pm1(x));
    }
    if (x <= 7) {
      final int n = (int) Math.floor(x - 0.5);
      double prod = 1;
      for (int i = 0; i < n; i++) {
        prod *= x - i;
      }
      return Math.log(prod) - Math.log1p(Gamma.invGamma1pm1(x - n));
    }
    final double xp1 = x + 1;
    final double sum = Gamma.lanczos(xp1);
    final double tmp = x + LANCZOS_G_PLUS_1_5;
    return (x + 1.5) * Math.log(tmp) - tmp + LN_SQRT_TWO_PI + Math.log(sum / xp1);
  }
}
