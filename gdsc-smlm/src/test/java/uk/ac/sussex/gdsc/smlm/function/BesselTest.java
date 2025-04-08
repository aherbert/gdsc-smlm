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

import java.util.function.DoubleUnaryOperator;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvFileSource;
import uk.ac.sussex.gdsc.test.api.Predicates;
import uk.ac.sussex.gdsc.test.api.TestAssertions;
import uk.ac.sussex.gdsc.test.api.function.DoubleDoubleBiPredicate;

@SuppressWarnings({"javadoc"})
class BesselTest {
  // Require doublesAreEqual to compare infinity == infinity
  private static final DoubleDoubleBiPredicate tolerance5e16 =
      Predicates.doublesAreEqual().or(Predicates.doublesAreRelativelyClose(5e-16));
  private static final DoubleDoubleBiPredicate tolerance1e15 =
      Predicates.doublesAreRelativelyClose(1e-15);
  private static final DoubleDoubleBiPredicate tolerance5e15 =
      Predicates.doublesAreRelativelyClose(5e-15);

  @ParameterizedTest
  @CsvFileSource(resources = {"bessel_i0.csv"})
  void testI0(double x, double expected) {
    assertBessel("i0", Bessel::i0, false, x, expected, tolerance5e16);
  }

  @ParameterizedTest
  @CsvFileSource(resources = {"bessel_i1.csv"})
  void testI1(double x, double expected) {
    assertBessel("i1", Bessel::i1, true, x, expected, tolerance5e16);
  }

  @ParameterizedTest
  @CsvFileSource(resources = {"bessel_j0.csv"})
  void testJ0(double x, double expected) {
    assertBessel("j0", Bessel::j0, false, x, expected, tolerance1e15);
  }

  @ParameterizedTest
  @CsvFileSource(resources = {"bessel_j1.csv"})
  void testJ1(double x, double expected) {
    assertBessel("j1", Bessel::j1, true, x, expected, tolerance5e15);
  }

  private static void assertBessel(String name, DoubleUnaryOperator bessel, boolean odd, double x,
      double expected, DoubleDoubleBiPredicate test) {
    final double i = bessel.applyAsDouble(x);
    TestAssertions.assertTest(expected, i, test, name);
    Assertions.assertEquals(odd ? -i : i, bessel.applyAsDouble(-x),
        () -> name + " is not " + (odd ? "odd" : "even"));
  }
}
