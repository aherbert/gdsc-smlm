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
import org.apache.commons.math3.special.Gamma;
import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import uk.ac.sussex.gdsc.test.api.TestAssertions;
import uk.ac.sussex.gdsc.test.api.TestHelper;
import uk.ac.sussex.gdsc.test.api.function.DoubleDoubleBiPredicate;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngUtils;

@SuppressWarnings({"javadoc"})
class FactorialTest {
  /** The largest representable factorial using a double (n=170). */
  private static final int MAX_N = 170;

  @Test
  void testZero() {
    Assertions.assertEquals(1, Factorial.value(0), "0!");
    Assertions.assertEquals(1, Factorial.value(0.0), "0!");
  }

  @Test
  void testIllegalArgumentThrows() {
    // For performance this uses an IOOB exception
    Assertions.assertThrows(IndexOutOfBoundsException.class, () -> Factorial.value(-1));
    Assertions.assertThrows(IllegalArgumentException.class, () -> Factorial.value(-0.5));
    Assertions.assertThrows(IllegalArgumentException.class, () -> Factorial.value(Double.NaN));
  }

  @Test
  void testArgumentTooLarge() {
    // Long avoids overflow to negative
    for (long n = MAX_N + 1; n < Integer.MAX_VALUE; n *= 2) {
      Assertions.assertEquals(Double.POSITIVE_INFINITY, Factorial.value((int) n));
      Assertions.assertEquals(Double.POSITIVE_INFINITY, Factorial.value(n));
    }
  }

  @Test
  void testFactorialInteger() {
    // Start at 0!
    BigInteger value = BigInteger.ONE;
    for (int n = 1; n <= MAX_N; n++) {
      // n! = (n-1)! * n
      value = value.multiply(BigInteger.valueOf(n));
      final double expected = value.doubleValue();
      Assertions.assertEquals(expected, Factorial.value(n));
      Assertions.assertEquals(expected, Factorial.value((double) n));
    }
  }

  /**
   * Test the factorial of a fractional number against Commons Math gamma(1+n).
   */
  @SeededTest
  void testFactorialDouble(RandomSeed seed) {
    final UniformRandomProvider rng = RngUtils.create(seed.getSeed());
    final DoubleDoubleBiPredicate tol =
        TestHelper.doublesAreClose(1e-14).or(TestHelper.doublesEqual());
    for (int i = 0; i < 100; i++) {
      final double n = rng.nextDouble() * 180;
      final double expected = Gamma.gamma(1 + n);
      TestAssertions.assertTest(expected, Factorial.value(n), tol, () -> Double.toString(n));
    }
  }

  /**
   * Test the factorial of a fractional number against reference values.
   */
  @ParameterizedTest
  @CsvSource({
      // Computed using Boost gamma(1+n) with long double precision
      // @formatter:off
      "0.25, 0.906402477055477077939",
      "0.75, 0.919062526848883233914",
      "1.125, 1.05946053733091417382",
      "1.75, 1.60835942198554565977",
      "2.75, 4.42298841046025056293",
      "3.25, 8.28508514183522016498",
      "4.75, 78.7844810613232131788",
      "6.25, 1155.3810139199896867",
      "12.75, 3255990905.16494153091",
      "17.125, 508919418354999.991547",
      "27.75, 1.32099626069721824877e+29",
      "57.75, 8.50381139447823243123e+77",
      "97.75, 2.99327649630298942593e+153",
      "137.75, 2.01698432357024766181e+236",
      "154.25, 1.08954758738660925068e+272",
      "164.25, 1.17747769003558155385e+294",
      "170.25, 2.62296687867878940345e+307",
      "170.624, 1.79421175992481035917e+308",
      // Limit
      "170.6243769563027, 1.7976931348622301e+308",
      // Infinite
      "170.62437695630274, 1.79769313486249261288e+308",
      "170.625, 1.80346206550076047216e+308",
      "171.5, 1.62639753771045308624e+310",
      // @formatter:on
  })
  void testFactorial(double n, double expected) {
    final DoubleDoubleBiPredicate tol =
        TestHelper.doublesAreClose(1e-15).or(TestHelper.doublesEqual());
    TestAssertions.assertTest(expected, Factorial.value(n), tol, () -> Double.toString(n));
  }
}
