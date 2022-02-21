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

import java.math.BigInteger;
import java.util.ArrayList;
import org.apache.commons.math3.special.Gamma;
import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import uk.ac.sussex.gdsc.test.api.TestAssertions;
import uk.ac.sussex.gdsc.test.api.TestHelper;
import uk.ac.sussex.gdsc.test.api.function.DoubleDoubleBiPredicate;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngUtils;
import uk.ac.sussex.gdsc.test.utils.RandomSeed;

@SuppressWarnings({"javadoc"})
class LogFactorialTest {
  /** The largest representable factorial using a double (n=170). */
  private static final int MAX_N = 170;

  @Test
  void testZero() {
    Assertions.assertEquals(0, LogFactorial.value(0), "0!");
    Assertions.assertEquals(0, LogFactorial.value(0.0), "0!");
  }

  @Test
  void testIllegalArgumentThrows() {
    // For performance this uses an IOOB exception
    Assertions.assertThrows(IndexOutOfBoundsException.class, () -> LogFactorial.value(-1));
    Assertions.assertThrows(IllegalArgumentException.class, () -> LogFactorial.value(-0.5));
    Assertions.assertThrows(IllegalArgumentException.class, () -> LogFactorial.value(Double.NaN));
  }

  @Test
  void testLogFactorialInteger() {
    // Start at 0!
    BigInteger value = BigInteger.ONE;
    final ArrayList<BigInteger> values = new ArrayList<>(MAX_N + 1);
    values.add(value);
    for (int n = 1; n <= MAX_N; n++) {
      // n! = (n-1)! * n
      value = value.multiply(BigInteger.valueOf(n));
      values.add(value);
      final double expected = Math.log(value.doubleValue());
      Assertions.assertEquals(expected, LogFactorial.value(n));
      Assertions.assertEquals(expected, LogFactorial.value((double) n));
    }

    // Compute some more by log summation
    double base = Math.log(value.doubleValue());
    final DoubleDoubleBiPredicate tol = TestHelper.doublesAreClose(1e-15);
    for (int n = MAX_N + 1; n < MAX_N * 3; n++) {
      base += Math.log(n);
      final double v1 = LogFactorial.value(n);
      final double v2 = LogFactorial.value((double) n);
      Assertions.assertEquals(v1, v2);
      TestAssertions.assertTest(base, v1, tol);
      // Prevent drift in the summation
      base = v1;
    }
  }

  /**
   * Test the factorial of a fractional number against Commons Math logGamma(1+n).
   */
  @SeededTest
  void testLogFactorialDouble(RandomSeed seed) {
    final UniformRandomProvider rng = RngUtils.create(seed.get());
    final DoubleDoubleBiPredicate tol = TestHelper.doublesAreClose(1e-14);
    for (int i = 0; i < 200; i++) {
      final double n = rng.nextDouble() * 200;
      final double expected = n <= 1.5 ? Gamma.logGamma1p(n) : Gamma.logGamma(1 + n);
      TestAssertions.assertTest(expected, LogFactorial.value(n), tol, () -> Double.toString(n));
    }
  }

  /**
   * Test the factorial of a fractional number against reference values.
   */
  @ParameterizedTest
  @CsvSource({
      // Computed using Matlab gammaln(1+n) with long double precision
      // @formatter:off
      "0.25, -0.098271836421813169",
      "0.75, -0.084401121020485553",
      "1.125, 0.057759851530343874",
      "1.75, 0.47521466691493708",
      "2.75, 1.4868155785934172",
      "3.25, 2.1144569274503713",
      "4.75, 4.3667160366222859",
      "6.25, 7.0521854507385395",
      "12.75, 21.903762491828793",
      "17.125, 33.86331080630174",
      "27.75, 67.053353891702798",
      "57.75, 179.43956662288721",
      "97.75, 353.39188783368269",
      "137.75, 544.11168543338499",
      "154.25, 626.38890784701778",
      "164.25, 677.12339194008723",
      "170.25, 707.85792962271012",
      "170.624, 709.78077443669906",
      // Limit
      "170.6243769563027, 709.78271289338397",
      // Infinite for n!
      "170.62437695630274, 709.78271289338409",
      "170.625, 709.78591682948354",
      "171.5, 714.28774629753423",
      "201.25, 869.86189472479612",
      "271.25, 1252.2956116469838",
      "411.25, 2068.075277298678",
      // @formatter:on
  })
  void testLogFactorial(double n, double expected) {
    final DoubleDoubleBiPredicate tol = TestHelper.doublesAreClose(1e-15);
    TestAssertions.assertTest(expected, LogFactorial.value(n), tol, () -> Double.toString(n));
  }
}
