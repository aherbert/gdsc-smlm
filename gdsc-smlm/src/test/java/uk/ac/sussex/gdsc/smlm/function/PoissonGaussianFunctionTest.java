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

import java.util.function.Supplier;
import java.util.logging.Logger;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.apache.commons.math3.analysis.integration.UnivariateIntegrator;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.test.api.Predicates;
import uk.ac.sussex.gdsc.test.api.TestAssertions;
import uk.ac.sussex.gdsc.test.api.function.DoubleDoubleBiPredicate;
import uk.ac.sussex.gdsc.test.junit5.SpeedTag;
import uk.ac.sussex.gdsc.test.utils.TestComplexity;
import uk.ac.sussex.gdsc.test.utils.TestLogging;
import uk.ac.sussex.gdsc.test.utils.TestLogging.TestLevel;
import uk.ac.sussex.gdsc.test.utils.TestSettings;

@SuppressWarnings({"javadoc"})
class PoissonGaussianFunctionTest {
  private static Logger logger;

  @BeforeAll
  public static void beforeAll() {
    logger = Logger.getLogger(PoissonGaussianFunctionTest.class.getName());
  }

  @AfterAll
  public static void afterAll() {
    logger = null;
  }

  // Note: Realistic gain values:
  // Photometrics EM-CCD in CCD mode:
  // Gain = 3.0, 1.5, 0.5 e-/ADU => 0.33, 0.66, 2 ADU/e-
  // Noise = 14.3 to 5.8 electrons
  // Note that this has a large well capacity and may not be indicative of
  // a standard CCD
  // sCMOS camera:
  // Gain has a range of 1.44 to 5.8, mean = 1.737 ADU/e-
  // Noise variance 4 to 141 ADUs => 1.1 to 7 electrons

  static double[] gain = {0.25, 0.5, 1, 2, 4}; // ADU/electron
  static double[] photons = {-1, 0, 0.1, 0.25, 0.5, 1, 2, 4, 10, 100, 1000};
  static double[] noise = {1, 2, 4, 8}; // electrons

  @Test
  void cumulativeProbabilityIsOneWithPicard() {
    for (final double g : gain) {
      for (final double p : photons) {
        for (final double s : noise) {
          cumulativeProbabilityIsOne(g, p, s, true);
        }
      }
    }
  }

  @Test
  void cumulativeProbabilityIsOneWithPade() {
    for (final double g : gain) {
      for (final double p : photons) {
        for (final double s : noise) {
          cumulativeProbabilityIsOne(g, p, s, false);
        }
      }
    }
  }

  @Test
  void cumulativeProbabilityIsNotOneWhenMeanIsLowAndNoiseIsLow() {
    // The cumulative probability is poor for low mean and low noise.
    // It can over-predict or under predict. The pattern of over/under is unknown.
    // For example in the following:

    // OVER
    Assertions.assertTrue(1.02 < cumulativeProbability(1.7, 0.25, 0.01, true));
    Assertions.assertTrue(1.02 < cumulativeProbability(1.7, 0.25, 0.1, true));
    // OK
    Assertions.assertEquals(1, cumulativeProbability(1.7, 0.25, 0.3, true), 0.02);
    // UNDER
    Assertions.assertTrue(0.98 > cumulativeProbability(1.7, 0.25, 0.5, true));
    // OK
    Assertions.assertEquals(1, cumulativeProbability(1.7, 0.25, 0.75, true), 0.02);

    // Fine with higher mean
    Assertions.assertEquals(1, cumulativeProbability(1.7, 10, 0.01, true), 0.02);
    Assertions.assertEquals(1, cumulativeProbability(1.7, 10, 0.1, true), 0.02);
    Assertions.assertEquals(1, cumulativeProbability(1.7, 10, 0.3, true), 0.02);
    Assertions.assertEquals(1, cumulativeProbability(1.7, 10, 0.5, true), 0.02);
    Assertions.assertEquals(1, cumulativeProbability(1.7, 10, 0.75, true), 0.02);
  }


  private static void cumulativeProbabilityIsOne(final double gain, final double mu,
      final double sd, final boolean usePicard) {
    final double p2 = cumulativeProbability(gain, mu, sd, usePicard);
    Assertions.assertEquals(1, p2, 0.02, () -> String.format("g=%f, mu=%f, s=%f", gain, mu, sd));
  }

  private static double cumulativeProbability(final double gain, final double mu, final double sd,
      final boolean usePicard) {
    // Note: The input mu & s parameters are pre-gain.
    final PoissonGaussianFunction f =
        PoissonGaussianFunction.createWithStandardDeviation(1.0 / gain, mu, sd * gain);
    f.setUsePicardApproximation(usePicard);
    double pvalue = 0;
    int min = 1;
    int max = 0;

    // Evaluate an initial range.
    // Gaussian should have >99% within +/- s
    // Poisson will have mean mu with a variance mu.
    // At large mu it is approximately normal so use 3 sqrt(mu) for the range added to the mean
    if (mu > 0) {
      final int[] range = getRange(gain, mu, sd);
      min = range[0];
      max = range[1];
      for (int x = min; x <= max; x++) {
        final double pp = f.probability(x);
        // logger.fine(FunctionUtils.getSupplier("x=%d, p=%f", x, pp);
        pvalue += pp;
      }
      // if (p > 1.01)
      // Assertions.fail("P > 1: " + p);
    }

    // We have most of the probability density.
    // Now keep evaluating up and down until no difference
    final double changeTolerance = 1e-6;
    for (int x = min - 1;; x--) {
      min = x;
      final double pp = f.probability(x);
      // logger.fine(FunctionUtils.getSupplier("x=%d, p=%f", x, pp);
      pvalue += pp;
      if (pp == 0 || pp / pvalue < changeTolerance) {
        break;
      }
    }
    for (int x = max + 1;; x++) {
      max = x;
      final double pp = f.probability(x);
      // logger.fine(FunctionUtils.getSupplier("x=%d, p=%f", x, pp);
      pvalue += pp;
      if (pp == 0 || pp / pvalue < changeTolerance) {
        break;
      }
    }

    // Do a formal integration
    double p2 = 0;
    final UnivariateIntegrator in =
        new SimpsonIntegrator(1e-6, 1e-6, 4, SimpsonIntegrator.SIMPSON_MAX_ITERATIONS_COUNT);
    p2 = in.integrate(Integer.MAX_VALUE, new UnivariateFunction() {
      @Override
      public double value(double x) {
        return f.probability(x);
      }
    }, min, max);

    if (p2 < 0.98 || p2 > 1.02) {
      logger.log(TestLogging.getRecord(TestLevel.TEST_INFO, "g=%f, mu=%f, s=%f p=%f  %f", gain, mu,
          sd, pvalue, p2));
    }

    return p2;
  }

  static int[] getRange(final double gain, final double mu, final double sd) {
    // Evaluate an initial range.
    // Gaussian should have >99% within +/- s
    // Poisson will have mean mu with a variance mu.
    // At large mu it is approximately normal so use 3 sqrt(mu) for the range added to the mean
    final double range = Math.max(sd, Math.sqrt(mu));
    final int min = (int) Math.floor(gain * (mu - 3 * range));
    final int max = (int) Math.ceil(gain * (mu + 3 * range));
    return new int[] {min, max};
  }

  @Test
  void probabilityMatchesLogProbability() {
    for (final double g : gain) {
      for (final double p : photons) {
        for (final double s : noise) {
          probabilityMatchesLogProbability(g, p, s, true);
          probabilityMatchesLogProbability(g, p, s, false);
        }
      }
    }
  }

  private static void probabilityMatchesLogProbability(final double gain, final double mu,
      final double sd, final boolean usePicard) {
    // Note: The input mu & s parameters are pre-gain.
    final PoissonGaussianFunction f =
        PoissonGaussianFunction.createWithStandardDeviation(1.0 / gain, mu, sd * gain);
    f.setUsePicardApproximation(usePicard);

    // Evaluate an initial range.
    // Gaussian should have >99% within +/- s
    // Poisson will have mean mu with a variance mu.
    // At large mu it is approximately normal so use 3 sqrt(mu) for the range added to the mean
    final int[] range = getRange(gain, mu, sd);
    final int min = range[0];
    final int max = range[1];
    final Supplier<String> msg = () -> String.format("g=%f, mu=%f, s=%f", gain, mu, sd);
    final DoubleDoubleBiPredicate predicate = Predicates.doublesAreClose(1e-3, 0);
    for (int x = min; x <= max; x++) {
      final double p = f.probability(x);
      if (p == 0) {
        continue;
      }
      final double logP = f.logProbability(x);
      TestAssertions.assertTest(Math.log(p), logP, predicate, msg);
    }
  }

  @SpeedTag
  @Test
  void padeIsFaster() {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final double[] noise2 = new double[noise.length];
    for (int i = 0; i < noise.length; i++) {
      noise2[i] = noise[i] * noise[i];
    }

    final int size = 100;
    final double[][] x = new double[photons.length][size];
    for (int i = 0; i < photons.length; i++) {
      final double p = photons[i] * 2 / size;
      for (int j = 0; j < size; j++) {
        x[i][j] = p * j;
      }
    }

    final long t1 = getTime(noise2, x, true);
    final long t2 = getTime(noise2, x, false);

    logger.log(TestLogging.getRecord(TestLevel.TEST_INFO, "Picard %d : Pade %d (%fx)", t1, t2,
        t1 / (double) t2));
    Assertions.assertTrue(t2 < t1, () -> String.format("Picard %d < Pade %d", t1, t2));
  }

  private static long getTime(double[] noise2, double[][] x, final boolean usePicard) {
    final int size = x[0].length;

    // Warm up
    for (final double s2 : noise2) {
      for (int i = 0; i < photons.length; i++) {
        final double p = photons[i];
        for (int j = 0; j < size; j++) {
          PoissonGaussianFunction.probability(x[i][j], p, s2, usePicard);
        }
      }
    }

    // Time
    long t1 = System.nanoTime();
    for (final double s2 : noise2) {
      for (int i = 0; i < photons.length; i++) {
        final double p = photons[i];
        for (int j = 0; j < size; j++) {
          PoissonGaussianFunction.probability(x[i][j], p, s2, usePicard);
        }
      }
    }
    t1 = System.nanoTime() - t1;
    return t1;
  }

  @Test
  void staticMethodsMatchInstanceMethods() {
    for (final double g : gain) {
      for (final double p : photons) {
        for (final double s : noise) {
          staticMethodsMatchInstanceMethods(g, p, s, true);
          staticMethodsMatchInstanceMethods(g, p, s, false);
        }
      }
    }
  }

  private static void staticMethodsMatchInstanceMethods(final double gain, final double mu,
      final double sd, final boolean usePicard) {
    // Note: The input mu & s parameters are pre-gain.
    final PoissonGaussianFunction f =
        PoissonGaussianFunction.createWithStandardDeviation(1.0 / gain, mu, sd * gain);
    f.setUsePicardApproximation(usePicard);

    // Evaluate an initial range.
    // Gaussian should have >99% within +/- s
    // Poisson will have mean mu with a variance mu.
    // At large mu it is approximately normal so use 3 sqrt(mu) for the range added to the mean
    final int[] range = getRange(gain, mu, sd);
    final int min = range[0];
    final int max = range[1];
    final double logGain = Math.log(gain);
    final double s2 = MathUtils.pow2(sd);
    final Supplier<String> msg1 =
        () -> String.format("probability g=%f, mu=%f, s=%f", gain, mu, sd);
    final Supplier<String> msg2 =
        () -> String.format("logProbability g=%f, mu=%f, s=%f", gain, mu, sd);
    for (int x = min; x <= max; x++) {
      double p1 = f.probability(x);
      double p2 = PoissonGaussianFunction.probability(x / gain, mu, s2, usePicard) / gain;
      Assertions.assertEquals(p1, p2, 1e-10, msg1);

      p1 = f.logProbability(x);
      p2 = PoissonGaussianFunction.logProbability(x / gain, mu, s2, usePicard) - logGain;
      Assertions.assertEquals(p1, p2, 1e-10, msg2);
    }
  }
}
