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

import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import java.util.Arrays;
import java.util.function.Supplier;
import java.util.logging.Logger;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.apache.commons.math3.analysis.integration.UnivariateIntegrator;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.test.api.TestAssertions;
import uk.ac.sussex.gdsc.test.api.TestHelper;
import uk.ac.sussex.gdsc.test.api.function.DoubleDoubleBiPredicate;
import uk.ac.sussex.gdsc.test.utils.TestComplexity;
import uk.ac.sussex.gdsc.test.utils.TestLogUtils;
import uk.ac.sussex.gdsc.test.utils.TestLogUtils.TestLevel;
import uk.ac.sussex.gdsc.test.utils.TestSettings;
import uk.ac.sussex.gdsc.test.utils.functions.FunctionUtils;

@SuppressWarnings({"javadoc"})
class PoissonGammaFunctionTest {
  private static Logger logger;

  @BeforeAll
  public static void beforeAll() {
    logger = Logger.getLogger(PoissonGammaFunctionTest.class.getName());
  }

  @AfterAll
  public static void afterAll() {
    logger = null;
  }

  static double[] gain = {6, 16, 30}; // ADU/electron above 1
  static double[] photons = {0.001, 0.1, 0.25, 0.5, 1, 2, 4, 10, 100, 1000};

  @Test
  void cumulativeProbabilityIsOneWithPmf() {
    for (final double g : gain) {
      for (final double p : photons) {
        cumulativeProbabilityIsOne(g, p, false);
      }
    }
  }

  @Test
  void cumulativeProbabilityIsOneWithPdf() {
    for (final double g : gain) {
      for (final double p : photons) {
        cumulativeProbabilityIsOne(g, p, true);
      }
    }
  }

  private static void cumulativeProbabilityIsOne(final double gain, final double mu, boolean pdf) {
    final double p2 = cumulativeProbability(gain, mu, pdf);

    if (pdf) {
      Assertions.assertEquals(1, p2, 0.02,
          () -> String.format("g=%f, mu=%f, pdf=%b", gain, mu, pdf));

      // This is not actually a PMF but is a PDF so requires integration.
      // This only works when the mean is above 2 if the gain is low
    } else if (mu > 2 || gain > 20) {
      Assertions.assertEquals(1, p2, 0.02,
          () -> String.format("g=%f, mu=%f, pdf=%b", gain, mu, pdf));
    }
  }

  private static double cumulativeProbability(final double gain, final double mu, boolean pdf) {
    final PoissonGammaFunction f = PoissonGammaFunction.createWithAlpha(1.0 / gain);

    double pvalue = 0;
    int min = 1;
    int max = 0;

    // Note: The input mu parameter is pre-gain.
    final double e = mu;

    final boolean debug = false;

    // Evaluate an initial range.
    // Gaussian should have >99% within +/- s
    // Poisson will have mean mu with a variance mu.
    // At large mu it is approximately normal so use 3 sqrt(mu) for the range added to the mean
    if (mu > 0) {
      // Note: The input s parameter is after-gain so adjust.
      final int[] range = PoissonGaussianFunctionTest.getRange(gain, mu, 0);
      min = range[0];
      max = range[1];
      for (int x = min; x <= max; x++) {
        final double pp = f.likelihood(x, e);
        // logger.fine(FunctionUtils.getSupplier("x=%d, p=%g", x, pp);
        if (debug) {
          logger.fine(FunctionUtils.getSupplier("x=%d, p=%f", x, pp));
        }
        pvalue += pp;
      }
      // if (p > 1.01)
      // Assertions.fail("P > 1: " + p);
    }

    // We have most of the likelihood density.
    // Now keep evaluating up and down until no difference
    final double changeTolerance = 1e-6;
    for (int x = min - 1;; x--) {
      min = x;
      final double pp = f.likelihood(x, e);
      // logger.fine(FunctionUtils.getSupplier("x=%d, p=%g", x, pp);
      if (debug) {
        logger.fine(FunctionUtils.getSupplier("x=%d, p=%f", x, pp));
      }
      pvalue += pp;
      if (pp == 0 || pp / pvalue < changeTolerance) {
        break;
      }
    }
    for (int x = max + 1;; x++) {
      max = x;
      final double pp = f.likelihood(x, e);
      // logger.fine(FunctionUtils.getSupplier("x=%d, p=%g", x, pp);
      if (debug) {
        logger.fine(FunctionUtils.getSupplier("x=%d, p=%f", x, pp));
      }
      pvalue += pp;
      if (pp == 0 || pp / pvalue < changeTolerance) {
        break;
      }
    }

    double p2 = pvalue;
    if (pdf) {
      // Do a formal integration
      if (debug) {
        if (pvalue < 0.98 || pvalue > 1.02) {
          logger.fine(FunctionUtils.getSupplier("g=%f, mu=%f, p=%f", gain, mu, pvalue));
        }
      }
      final UnivariateIntegrator in =
          new SimpsonIntegrator(1e-4, 1e-6, 4, SimpsonIntegrator.SIMPSON_MAX_ITERATIONS_COUNT);
      p2 = in.integrate(Integer.MAX_VALUE, new UnivariateFunction() {
        @Override
        public double value(double x) {
          // return f.likelihood(x, e);
          return PoissonGammaFunction.poissonGammaN(x, mu, gain);
        }
      }, min, max);

      p2 += PoissonGammaFunction.dirac(mu);
    }

    // if (p2 < 0.98 || p2 > 1.02)
    logger.log(
        TestLogUtils.getRecord(TestLevel.TEST_INFO, "g=%f, mu=%f, p=%f  %f", gain, mu, pvalue, p2));

    return p2;
  }

  @Test
  void probabilityMatchesLogProbability() {
    for (final double g : gain) {
      for (final double p : photons) {
        probabilityMatchesLogProbability(g, p);
      }
    }
  }

  private static void probabilityMatchesLogProbability(final double gain, double mu) {
    final PoissonGammaFunction f = PoissonGammaFunction.createWithAlpha(1.0 / gain);

    // Evaluate an initial range.
    // Gaussian should have >99% within +/- s
    // Poisson will have mean mu with a variance mu.
    // At large mu it is approximately normal so use 3 sqrt(mu) for the range added to the mean
    // Note: The input s parameter is after-gain so adjust.
    final int[] range = PoissonGaussianFunctionTest.getRange(gain, mu, 0);
    final int min = range[0];
    final int max = range[1];
    // Note: The input mu parameter is pre-gain.
    final double e = mu;
    final Supplier<String> msg = () -> String.format("g=%f, mu=%f", gain, mu);
    final DoubleDoubleBiPredicate predicate = TestHelper.doublesAreClose(1e-6, 0);
    for (int x = min; x <= max; x++) {
      final double p = f.likelihood(x, e);
      if (p == 0) {
        continue;
      }
      final double logP = f.logLikelihood(x, e);
      TestAssertions.assertTest(Math.log(p), logP, predicate, msg);
    }
  }

  @Test
  void canComputePoissonGammaGradientWithInteger() {
    for (int j = 0; j < gain.length; j++) {
      for (int i = 0; i < photons.length; i++) {
        canComputePoissonGammaGradient(gain[j], photons[i], false);
      }
    }
  }

  @Test
  void canComputePoissonGammaGradientWithReal() {
    for (int j = 0; j < gain.length; j++) {
      for (int i = 0; i < photons.length; i++) {
        canComputePoissonGammaGradient(gain[j], photons[i], true);
      }
    }
  }

  @SuppressWarnings("unused")
  private static void canComputePoissonGammaGradient(final double gain, final double mu,
      boolean nonInteger) {
    final double o = mu;
    final double delta = 1e-3; // * o;
    final double uo = o + delta;
    final double lo = o - delta;
    final double diff = uo - lo;

    // The numerical gradient is poor around the switch between the use of the
    // Bessel function and the approximation. So just count the errors.
    int fail = 0;
    double sum = 0;

    final int[] range = PoissonGaussianFunctionTest.getRange(gain, mu, 0);
    final int min = Math.max(0, range[0]);
    final int max = range[1];
    final double[] dp_dt = new double[1];
    final double[] dp_dt2 = new double[1];
    final double step = (nonInteger) ? 0.5 : 1;

    // When using the approximation the gradients are not as accurate
    final boolean approx = (2 * Math.sqrt(max * o / gain) > 709);
    final double tol = approx ? 0.05 : 1e-3;

    final DoubleArrayList list = new DoubleArrayList();
    if (min != 0) {
      list.add(0);
    }
    for (double x = min; x <= max; x += step) {
      list.add(x);
    }

    final DoubleDoubleBiPredicate predicate = TestHelper.doublesAreClose(1e-8, 0);
    for (final double x : list.toDoubleArray()) {
      final double p1 = PoissonGammaFunction.poissonGamma(x, o, gain);
      final double p2 = PoissonGammaFunction.poissonGamma(x, o, gain, dp_dt);
      Assertions.assertEquals(p1, p2);

      // Check partial gradient matches
      double p3 = PoissonGammaFunction.poissonGammaPartial(x, o, gain, dp_dt2);
      Assertions.assertEquals(p1, p3);
      TestAssertions.assertTest(dp_dt[0] + p1, dp_dt2[0], predicate);

      // Check no dirac gradient matches
      p3 = PoissonGammaFunction.poissonGammaN(x, o, gain, dp_dt2);
      if (x == 0) {
        final double dirac = PoissonGammaFunction.dirac(o);
        // Add the dirac contribution
        p3 += dirac;
        dp_dt2[0] -= dirac;
        TestAssertions.assertTest(p1, p3, predicate);
        TestAssertions.assertTest(dp_dt[0], dp_dt2[0], predicate);
      } else {
        Assertions.assertEquals(p1, p3);
        TestAssertions.assertTest(dp_dt[0], dp_dt2[0], predicate);
      }

      final double up = PoissonGammaFunction.poissonGamma(x, uo, gain);
      final double lp = PoissonGammaFunction.poissonGamma(x, lo, gain);

      final double eg = dp_dt[0];
      final double g = (up - lp) / diff;
      final double error = DoubleEquality.relativeError(g, eg);
      final double ox = x / gain;
      // logger.fine(FunctionUtils.getSupplier("g=%g, mu=%g, x=%g (ox=%g), p=%g g=%g %g error=%g",
      // gain, mu, x, ox, p1, g, eg,
      // error);

      if (error > tol) {
        fail++;
        sum += error;
      }
    }

    final double f = (double) fail / list.size();
    logger.log(TestLogUtils.getRecord(TestLevel.TEST_INFO, "g=%g, mu=%g, failures=%g, mean=%f",
        gain, mu, f, MathUtils.div0(sum, fail)));
    if (approx) {
      Assertions.assertTrue(f < 0.2);
    } else {
      Assertions.assertTrue(f < 0.01);
    }
  }

  @Test
  void canComputeSeparatelyAtC0() {
    final DoubleArrayList list = new DoubleArrayList();
    for (int exp = -12; exp < 6; exp++) {
      list.add(Math.pow(10, exp * 0.5));
    }
    for (int x = 2; x < 10; x++) {
      list.add(x / 10.0);
    }
    for (double x = 11; x <= 20; x++) {
      list.add(x / 10.0);
    }
    final double[] p = list.toDoubleArray();

    final boolean report =
        logger.isLoggable(TestLevel.TEST_INFO) && TestSettings.allow(TestComplexity.MEDIUM);
    if (report) {
      Arrays.sort(p);
    }

    final double m = 5;

    final DoubleDoubleBiPredicate predicate = TestHelper.doublesAreClose(1e-8, 0);
    for (final double x : p) {
      final double e = PoissonGammaFunction.poissonGamma(0, x, m);
      // Test the function can be separated into the dirac and the rest
      final double dirac = PoissonGammaFunction.dirac(x);
      final double p0 = PoissonGammaFunction.poissonGammaN(0, x, m);
      TestAssertions.assertTest(e, dirac + p0, predicate);

      // For reporting
      if (report) {
        final double p01 = PoissonGammaFunction.poissonGammaN(1e-10, x, m);

        logger.log(TestLogUtils.getRecord(TestLevel.TEST_INFO,
            "p=%g  Dirac=%s   p0=%s (dirac:p0=%s)   p01=%s  (p0:p01 = %s)", x, dirac, p0,
            dirac / p0,
            // uk.ac.sussex.gdsc.core.utils.DoubleEquality.relativeError(p0, dirac),
            p01, p0 / p01
        // uk.ac.sussex.gdsc.core.utils.DoubleEquality.relativeError(p0, p01)
        ));
      }
    }
  }
}
