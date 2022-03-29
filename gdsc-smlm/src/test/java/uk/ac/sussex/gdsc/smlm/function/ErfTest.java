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

import java.util.logging.Logger;
import org.apache.commons.math3.util.Precision;
import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.smlm.utils.StdMath;
import uk.ac.sussex.gdsc.test.api.Predicates;
import uk.ac.sussex.gdsc.test.api.TestAssertions;
import uk.ac.sussex.gdsc.test.api.function.DoubleDoubleBiPredicate;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngFactory;
import uk.ac.sussex.gdsc.test.utils.RandomSeed;
import uk.ac.sussex.gdsc.test.utils.TestLogging;
import uk.ac.sussex.gdsc.test.utils.TestLogging.TestLevel;
import uk.ac.sussex.gdsc.test.utils.functions.FormatSupplier;

// TODO: Update this using the new Erf implementation in Commons Numbers 1.1

@SuppressWarnings({"javadoc"})
class ErfTest {
  private static Logger logger;

  @BeforeAll
  public static void beforeAll() {
    logger = Logger.getLogger(ErfTest.class.getName());
  }

  @AfterAll
  public static void afterAll() {
    logger = null;
  }

  //@formatter:off
  private abstract static class BaseErf {
    String name;
    BaseErf(String name) { this.name = name; }
    abstract double erf(double x);
    abstract double erf(double x1, double x2);
  }
  private static class ApacheErf extends BaseErf {
    ApacheErf() {  super("apache erf"); }
    @Override
    double erf(double x) { return org.apache.commons.math3.special.Erf.erf(x); }
    @Override
    double erf(double x1, double x2) { return org.apache.commons.math3.special.Erf.erf(x1, x2); }
  }
  private static class Erf extends BaseErf {
    Erf() {  super("erf"); }
    @Override
    double erf(double x) { return uk.ac.sussex.gdsc.smlm.function.Erf.erf(x); }
    @Override
    double erf(double x1, double x2) { return uk.ac.sussex.gdsc.smlm.function.Erf.erf(x1, x2); }
  }
  private static class Erf0 extends BaseErf {
    Erf0() { super("erf0"); }
    @Override
    double erf(double x) { return uk.ac.sussex.gdsc.smlm.function.Erf.erf0(x); }
    @Override
    double erf(double x1, double x2) { return uk.ac.sussex.gdsc.smlm.function.Erf.erf0(x1, x2); }
  }
  private static class Erf2 extends BaseErf {
    Erf2() { super("erf2"); }
    @Override
    double erf(double x) { return uk.ac.sussex.gdsc.smlm.function.Erf.erf2(x); }
    @Override
    double erf(double x1, double x2) { return uk.ac.sussex.gdsc.smlm.function.Erf.erf2(x1, x2); }
  }
  //@formatter:on

  @SeededTest
  void erf0xHasLowError(RandomSeed seed) {
    checkErfxHasLowError(seed, new Erf0(), 5e-4);
  }

  @SeededTest
  void erfxHasLowError(RandomSeed seed) {
    checkErfxHasLowError(seed, new Erf(), 3e-7);
  }

  @SeededTest
  void erf2xHasLowError(RandomSeed seed) {
    checkErfxHasLowError(seed, new Erf2(), 1.3e-4);
  }

  private static void checkErfxHasLowError(RandomSeed seed, BaseErf erf, double expected) {
    final UniformRandomProvider rg = RngFactory.create(seed.get());
    final int range = 8;
    double max = 0;

    for (int xi = -range; xi <= range; xi++) {
      for (int i = 0; i < 5; i++) {
        final double x = xi + rg.nextDouble();
        final double o = erf.erf(x);
        final double e = org.apache.commons.math3.special.Erf.erf(x);
        final double error = Math.abs(o - e);
        if (max < error) {
          max = error;
        }
        // logger.fine(FormatSupplier.getSupplier("x=%f, e=%f, o=%f, error=%f", x, e, o, error);
        Assertions.assertTrue(error < expected);
      }
    }
    logger
        .log(TestLogging.getRecord(TestLevel.TEST_INFO, "erfx %s max error = %g", erf.name, max));
  }

  @Test
  void erfApachexIndistinguishableFrom1() {
    checkErfxIndistinguishableFrom1(new ApacheErf());
  }

  @Test
  void erf0xIndistinguishableFrom1() {
    checkErfxIndistinguishableFrom1(new Erf0());
  }

  @Test
  void erfxIndistinguishableFrom1() {
    checkErfxIndistinguishableFrom1(new Erf());
  }

  @Test
  void erf2xIndistinguishableFrom1() {
    checkErfxIndistinguishableFrom1(new Erf2());
  }

  private static void checkErfxIndistinguishableFrom1(BaseErf erf) {
    Assumptions.assumeTrue(logger.isLoggable(TestLevel.TEST_INFO));

    // Find switch using a binary search
    double lower = 1;
    double upper = 40;
    while (DoubleEquality.complement(lower, upper) > 1) {
      final double mid = (upper + lower) * 0.5;
      final double o = erf.erf(mid);
      if (o == 1) {
        upper = mid;
      } else {
        lower = mid;
      }
    }

    logger.log(TestLevel.TEST_INFO,
        FormatSupplier.getSupplier("erfx %s indistinguishable from 1: x > %s, x >= %s", erf.name,
            Double.toString(lower), Double.toString(upper)));
  }

  @SeededTest
  void erf0xxHasLowError(RandomSeed seed) {
    checkErfxxHasLowError(seed, new Erf0(), 4e-2);
  }

  @SeededTest
  void erfxxHasLowError(RandomSeed seed) {
    checkErfxxHasLowError(seed, new Erf(), 7e-4);
  }

  @SeededTest
  void erf2xxHasLowError(RandomSeed seed) {
    checkErfxxHasLowError(seed, new Erf2(), 1.1e-2);
  }

  private static void checkErfxxHasLowError(RandomSeed seed, BaseErf erf, double expected) {
    final UniformRandomProvider rg = RngFactory.create(seed.get());

    final int range = 3;
    double max = 0;

    for (int xi = -range; xi <= range; xi++) {
      for (int xi2 = -range; xi2 <= range; xi2++) {
        for (int i = 0; i < 5; i++) {
          final double x = xi + rg.nextDouble();
          for (int j = 0; j < 5; j++) {
            final double x2 = xi2 + rg.nextDouble();

            final double o = erf.erf(x, x2);
            final double e = org.apache.commons.math3.special.Erf.erf(x, x2);
            final double error = Math.abs(o - e);
            if (max < error) {
              max = error;
            }
            // logger.fine(FormatSupplier.getSupplier("x=%f, x2=%f, e=%f, o=%f, error=%f", x, x2, e,
            // o, error);
            Assertions.assertTrue(error < expected);
          }
        }
      }
    }

    logger
        .log(TestLogging.getRecord(TestLevel.TEST_INFO, "erfxx %s max error = %g", erf.name, max));
  }

  @Test
  void erf0xxHasLowErrorForUnitBlocks() {
    checkErfxxHasLowErrorForUnitBlocks(new Erf0(), 5e-4);
  }

  @Test
  void erfxxHasLowErrorForUnitBlocks() {
    checkErfxxHasLowErrorForUnitBlocks(new Erf(), 5e-7);
  }

  @Test
  void erf2xxHasLowErrorForUnitBlocks() {
    checkErfxxHasLowErrorForUnitBlocks(new Erf2(), 1e-4);
  }

  private static void checkErfxxHasLowErrorForUnitBlocks(BaseErf erf, double expected) {
    final int range = 8;
    double max = 0;

    for (int xi = -range; xi <= range; xi++) {
      final double x = xi;
      final double x2 = xi + 1;
      final double o = erf.erf(x, x2);
      final double e = org.apache.commons.math3.special.Erf.erf(x, x2);
      final double error = Math.abs(o - e);
      if (max < error) {
        max = error;
      }
      // logger.fine(FormatSupplier.getSupplier("x=%f, x2=%f, e=%f, o=%f, error=%f", x, x2, e, o,
      // error);
      Assertions.assertTrue(error < expected);
    }

    logger.log(
        TestLogging.getRecord(TestLevel.TEST_INFO, "erfxx %s unit max error = %g", erf.name, max));
  }

  @Test
  void erf0xxHasLowerErrorThanGaussianApproximationForUnitBlocks() {
    checkErfxxHasLowerErrorThanGaussianApproximationForUnitBlocks(new Erf0());
  }

  @Test
  void erfxxHasLowerErrorThanGaussianApproximationForUnitBlocks() {
    checkErfxxHasLowerErrorThanGaussianApproximationForUnitBlocks(new Erf());
  }

  @Test
  void erf2xxHasLowerErrorThanGaussianApproximationForUnitBlocks() {
    checkErfxxHasLowerErrorThanGaussianApproximationForUnitBlocks(new Erf2());
  }

  private static void checkErfxxHasLowerErrorThanGaussianApproximationForUnitBlocks(BaseErf erf) {
    final int range = 5;
    double max = 0;
    double max2 = 0;

    // Standard deviation
    final double s = 1.3;
    final double twos2 = 2 * s * s;
    final double norm = 1 / (Math.PI * twos2);
    final double denom = 1.0 / (Math.sqrt(2.0) * s);

    double sum1 = 0;
    double sum2 = 0;
    double sum3 = 0;

    for (int x = -range; x <= range; x++) {
      final double o1 = 0.5 * erf.erf((x - 0.5) * denom, (x + 0.5) * denom);
      final double e1 =
          0.5 * org.apache.commons.math3.special.Erf.erf((x - 0.5) * denom, (x + 0.5) * denom);
      for (int y = -range; y <= range; y++) {
        final double o2 = 0.5 * erf.erf((y - 0.5) * denom, (y + 0.5) * denom);
        final double e2 =
            0.5 * org.apache.commons.math3.special.Erf.erf((y - 0.5) * denom, (y + 0.5) * denom);

        final double o = o1 * o2;
        final double e = e1 * e2;
        final double oo = norm * Math.exp(-(x * x + y * y) / twos2);

        sum1 += e;
        sum2 += o;
        sum3 += oo;

        final double absError = Math.abs(o - e);
        if (e < 1e-4 || absError < 1e-10) {
          continue;
        }
        final double error = DoubleEquality.relativeError(o, e);
        final double error2 = DoubleEquality.relativeError(oo, e);
        if (max < error) {
          max = error;
        }
        if (max2 < error2) {
          max2 = error2;
        }
        // logger.fine(FormatSupplier.getSupplier("x=%d, y=%d, e=%g, o=%g, o2=%g, error=%f,
        // error2=%f", x, y, e, o, oo, error, error2);
        Assertions.assertTrue(error < error2);
      }
    }

    Assertions.assertTrue(sum1 > 0.999, () -> erf.name + " Gaussian 2D integral is not 1");
    Assertions.assertTrue(DoubleEquality.relativeError(sum1, sum2) < 1e-3,
        () -> erf.name + " Erf approx integral is incorrect");
    Assertions.assertTrue(DoubleEquality.relativeError(sum1, sum3) < 1e-3,
        () -> erf.name + " Gaussian approx integral is incorrect");

    logger.log(TestLogging.getRecord(TestLevel.TEST_INFO,
        "%s Erf approx pixel unit max error = %f", erf.name, max));
    logger.log(TestLogging.getRecord(TestLevel.TEST_INFO,
        "%s Gaussian approx pixel unit max error = %f", erf.name, max2));
  }

  @Test
  void gaussianIntegralApproximatesErf() {
    final double x = 1.3;
    final double y = 2.2;
    final double s = 1.14;
    final int minx = (int) x;
    final int miny = (int) y;
    final int maxx = minx + 1;
    final int maxy = miny + 1;

    // Full integration using the Erf
    // Note: The PSF of a 2D Gaussian is described in Smith et all using a denominator
    // of (2.0 * s * s) for both x and Y directions. This is wrong. We need the
    // integral of the single Guassian in each dimension so the denomiator is (sqrt(2.0) * s).
    // See: Smith et al, (2010). Fast, single-molecule localisation that achieves
    // theoretically minimum uncertainty. Nature Methods 7, 373-375
    // (supplementary note).
    // final double denom = 1.0 / (2.0 * s * s); // As per Smith, etal (2010),

    final double denom = 1.0 / (Math.sqrt(2.0) * s);
    final double e1 = 0.5 * org.apache.commons.math3.special.Erf.erf(minx * denom, maxx * denom);
    final double e2 = 0.5 * org.apache.commons.math3.special.Erf.erf(miny * denom, maxy * denom);
    final double expected = e1 * e2;

    double observed = 0;
    // Numeric integration
    final double twos2 = 2 * s * s;
    final double norm = 1 / (Math.PI * twos2);
    for (int i = 0, steps = 1; i < 4; i++, steps = (int) Math.pow(10, i)) {
      // Gaussian is: exp(-(x * x + y * y) / twos2) over all x and y
      // But we can do this by separating x and y:
      // exp(-(x * x) / twos2) * exp(-(y * y) / twos2)

      // pre-compute
      final double[] ex = new double[steps];
      double sumey = 0;
      if (steps == 1) {
        // Use the actual values for x and y
        ex[0] = StdMath.exp(-(x * x) / twos2);
        sumey = StdMath.exp(-(y * y) / twos2);
      } else {
        for (int j = 0; j < steps; j++) {
          final double xx = minx + (double) j / steps;
          final double yy = miny + (double) j / steps;
          ex[j] = StdMath.exp(-(xx * xx) / twos2);
          sumey += StdMath.exp(-(yy * yy) / twos2);
        }
      }

      double sum = 0;
      for (int j = 0; j < steps; j++) {
        sum += ex[j] * sumey;
      }

      //// Check
      // double sum2 = 0;
      // for (int j = 0; j <= steps; j++)
      // {
      // double xx = minx + (double) j / steps;
      // for (int k = 0; k <= steps; k++)
      // {
      // double yy = miny + (double) k / steps;
      // sum2 += StdMath.exp(-(xx * xx + yy * yy) / twos2);
      // }
      // }
      // logger.fine(FormatSupplier.getSupplier("sum=%f, sum2=%f", sum, sum2);

      final int n = steps * steps;
      observed = norm * sum / n;
      logger.log(TestLogging.getRecord(TestLevel.TEST_INFO, "n=%d, e=%f, o=%f, error=%f", n,
          expected, observed, DoubleEquality.relativeError(expected, observed)));
    }

    TestAssertions.assertTest(expected, observed, Predicates.doublesAreClose(1e-2, 0));
  }

  @Test
  void analyticErfGradientCorrectForErfApproximation() {
    final BaseErf erf = new Erf();
    final int range = 7;
    final int steps = 10000;
    final double step = (double) range / steps;
    final double delta = 1e-3;
    final DoubleEquality eq = new DoubleEquality(5e-4, 1e-6);
    for (int i = 0; i < steps; i++) {
      final double x = i * step;
      final double x1 = x + Precision.representableDelta(x, delta);
      final double x2 = x - Precision.representableDelta(x, delta);
      final double o1 = erf.erf(x1);
      final double o2 = erf.erf(x2);
      final double delta2 = x1 - x2;
      final double g = (o1 - o2) / delta2;
      final double e = uk.ac.sussex.gdsc.smlm.function.Erf.erfDerivative(x);
      if (!eq.almostEqualRelativeOrAbsolute(e, g)) {
        Assertions.fail(x + " : " + e + " != " + g);
      }
    }
  }

  @Test
  void canComputePower4() {
    final DoubleDoubleBiPredicate equality = Predicates.doublesAreClose(1e-10, 0);
    for (int i = -10; i <= 10; i++) {
      for (final double d : new double[] {0, 0.1, 0.01, 0.001}) {
        final double f = i + d;
        final double e = Math.pow(f, 4);
        final double o = uk.ac.sussex.gdsc.smlm.function.Erf.pow4(f);
        TestAssertions.assertTest(e, o, equality, () -> "x=" + f);
      }
    }
  }

  @Test
  void canComputePower16() {
    final DoubleDoubleBiPredicate equality = Predicates.doublesAreClose(1e-10, 0);
    for (int i = -10; i <= 10; i++) {
      for (final double d : new double[] {0, 0.1, 0.01, 0.001}) {
        final double f = i + d;
        final double e = Math.pow(f, 16);
        final double o = uk.ac.sussex.gdsc.smlm.function.Erf.pow16(f);
        TestAssertions.assertTest(e, o, equality, () -> "x=" + f);
      }
    }
  }
}
