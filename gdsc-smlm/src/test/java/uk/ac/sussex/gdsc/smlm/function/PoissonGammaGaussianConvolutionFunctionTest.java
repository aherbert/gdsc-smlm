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

import java.util.function.Supplier;
import java.util.logging.Logger;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.apache.commons.math3.analysis.integration.UnivariateIntegrator;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import uk.ac.sussex.gdsc.smlm.function.PoissonGammaGaussianFunction.ConvolutionMode;
import uk.ac.sussex.gdsc.test.api.Predicates;
import uk.ac.sussex.gdsc.test.api.TestAssertions;
import uk.ac.sussex.gdsc.test.api.function.DoubleDoubleBiPredicate;
import uk.ac.sussex.gdsc.test.utils.TestLogging;
import uk.ac.sussex.gdsc.test.utils.TestLogging.TestLevel;
import uk.ac.sussex.gdsc.test.utils.functions.FormatSupplier;

@SuppressWarnings({"javadoc"})
class PoissonGammaGaussianConvolutionFunctionTest {
  private static Logger logger;

  @BeforeAll
  public static void beforeAll() {
    logger = Logger.getLogger(PoissonGammaGaussianConvolutionFunctionTest.class.getName());
  }

  @AfterAll
  public static void afterAll() {
    logger = null;
  }

  static double[] gain = {6, 30}; // ADU/electron above 1
  static double[] photons = PoissonGaussianFunctionTest.photons;
  static double[] noise = {1, 10}; // ADUs

  @Test
  void cumulativeProbabilityIsOne() {
    for (final double g : gain) {
      for (final double p : photons) {
        for (final double s : noise) {
          cumulativeProbabilityIsOne(g, p, s);
        }
      }
    }
  }

  private static void cumulativeProbabilityIsOne(final double gain, final double mu,
      final double sd) {
    final double p2 = cumulativeProbability(gain, mu, sd);
    // This only works when the mean is above 2 if the gain is low
    if (mu > 2 || gain > 20) {
      Assertions.assertEquals(1, p2, 0.02, () -> String.format("g=%f, mu=%f, s=%f", gain, mu, sd));
    }
  }

  private static double cumulativeProbability(final double gain, final double mu, double sd) {
    final PoissonGammaGaussianConvolutionFunction f =
        PoissonGammaGaussianConvolutionFunction.createWithStandardDeviation(1.0 / gain, sd);

    final PoissonGammaGaussianFunction f2 = new PoissonGammaGaussianFunction(1.0 / gain, sd);
    f2.setConvolutionMode(ConvolutionMode.DISCRETE_PDF);

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
      final int[] range = PoissonGaussianFunctionTest.getRange(gain, mu, sd / gain);
      min = range[0];
      max = range[1];
      for (int x = min; x <= max; x++) {
        final double pp = f.likelihood(x, e);
        // logger.fine(FormatSupplier.getSupplier("x=%d, p=%g", x, pp);
        if (debug) {
          logger.fine(FormatSupplier.getSupplier("x=%d, p=%f   %f", x, pp, f2.likelihood(x, e)));
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
      // logger.fine(FormatSupplier.getSupplier("x=%d, p=%g", x, pp);
      if (debug) {
        logger.fine(FormatSupplier.getSupplier("x=%d, p=%f   %f", x, pp, f2.likelihood(x, e)));
      }
      pvalue += pp;
      if (pp == 0 || pp / pvalue < changeTolerance) {
        break;
      }
    }
    for (int x = max + 1;; x++) {
      max = x;
      final double pp = f.likelihood(x, e);
      // logger.fine(FormatSupplier.getSupplier("x=%d, p=%g", x, pp);
      if (debug) {
        logger.fine(FormatSupplier.getSupplier("x=%d, p=%f   %f", x, pp, f2.likelihood(x, e)));
      }
      pvalue += pp;
      if (pp == 0 || pp / pvalue < changeTolerance) {
        break;
      }
    }

    if (pvalue < 0.98 || pvalue > 1.02) {
      logger.log(TestLogging.getRecord(TestLevel.TEST_DEBUG, "g=%f, mu=%f, s=%f p=%f", gain, mu,
          sd, pvalue));
    }

    // Do a formal integration
    double p2 = 0;
    final UnivariateIntegrator in =
        new SimpsonIntegrator(1e-4, 1e-6, 4, SimpsonIntegrator.SIMPSON_MAX_ITERATIONS_COUNT);
    p2 = in.integrate(Integer.MAX_VALUE, new UnivariateFunction() {
      @Override
      public double value(double x) {
        return f.likelihood(x, e);
      }
    }, min, max);

    if (p2 < 0.98 || p2 > 1.02) {
      logger.log(TestLogging.getRecord(TestLevel.TEST_INFO, "g=%f, mu=%f, s=%f p=%f  %f", gain, mu,
          sd, pvalue, p2));
    }

    return p2;
  }

  @Test
  void probabilityMatchesLogProbability() {
    for (final double g : gain) {
      for (final double p : photons) {
        for (final double s : noise) {
          probabilityMatchesLogProbability(g, p, s);
        }
      }
    }
  }

  private static void probabilityMatchesLogProbability(final double gain, double mu, double sd) {
    final PoissonGammaGaussianConvolutionFunction f =
        PoissonGammaGaussianConvolutionFunction.createWithStandardDeviation(1.0 / gain, sd);

    // Evaluate an initial range.
    // Gaussian should have >99% within +/- s
    // Poisson will have mean mu with a variance mu.
    // At large mu it is approximately normal so use 3 sqrt(mu) for the range added to the mean
    // Note: The input s parameter is after-gain so adjust.
    final int[] range = PoissonGaussianFunctionTest.getRange(gain, mu, sd / gain);
    final int min = range[0];
    final int max = range[1];
    // Note: The input mu parameter is pre-gain.
    final double e = mu;
    final Supplier<String> msg = () -> String.format("g=%f, mu=%f, s=%f", gain, mu, sd);
    final DoubleDoubleBiPredicate predicate = Predicates.doublesAreClose(1e-3, 0);
    for (int x = min; x <= max; x++) {
      final double p = f.likelihood(x, e);
      if (p == 0) {
        continue;
      }
      final double logP = f.logLikelihood(x, e);
      TestAssertions.assertTest(Math.log(p), logP, predicate, msg);
    }
  }
}
