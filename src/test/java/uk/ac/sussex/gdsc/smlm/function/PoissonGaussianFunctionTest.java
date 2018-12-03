package uk.ac.sussex.gdsc.smlm.function;

import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.test.api.TestAssertions;
import uk.ac.sussex.gdsc.test.api.TestHelper;
import uk.ac.sussex.gdsc.test.api.function.DoubleDoubleBiPredicate;
import uk.ac.sussex.gdsc.test.junit5.SpeedTag;
import uk.ac.sussex.gdsc.test.utils.TestComplexity;
import uk.ac.sussex.gdsc.test.utils.TestLogUtils;
import uk.ac.sussex.gdsc.test.utils.TestSettings;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.apache.commons.math3.analysis.integration.UnivariateIntegrator;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import java.util.function.Supplier;
import java.util.logging.Level;
import java.util.logging.Logger;

@SuppressWarnings({"javadoc"})
public class PoissonGaussianFunctionTest {
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
  public void cumulativeProbabilityIsOneWithPicard() {
    for (final double g : gain) {
      for (final double p : photons) {
        for (final double s : noise) {
          cumulativeProbabilityIsOne(g, p, s, true);
        }
      }
    }
  }

  @Test
  public void cumulativeProbabilityIsOneWithPade() {
    for (final double g : gain) {
      for (final double p : photons) {
        for (final double s : noise) {
          cumulativeProbabilityIsOne(g, p, s, false);
        }
      }
    }
  }

  @Test
  public void cumulativeProbabilityIsNotOneWhenMeanIsLowAndNoiseIsLow() {
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

  @Test
  public void probabilityMatchesLogProbability() {
    for (final double g : gain) {
      for (final double p : photons) {
        for (final double s : noise) {
          probabilityMatchesLogProbability(g, p, s, true);
          probabilityMatchesLogProbability(g, p, s, false);
        }
      }
    }
  }

  @SpeedTag
  @Test
  public void padeIsFaster() {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final double[] noise2 = new double[noise.length];
    for (int i = 0; i < noise.length; i++) {
      noise2[i] = noise[i] * noise[i];
    }

    final int N = 100;
    final double[][] x = new double[photons.length][N];
    for (int i = 0; i < photons.length; i++) {
      final double p = photons[i] * 2 / N;
      for (int j = 0; j < N; j++) {
        x[i][j] = p * j;
      }
    }

    final long t1 = getTime(noise2, N, x, true);
    final long t2 = getTime(noise2, N, x, false);

    logger.log(
        TestLogUtils.getRecord(Level.INFO, "Picard %d : Pade %d (%fx)", t1, t2, t1 / (double) t2));
    Assertions.assertTrue(t2 < t1, () -> String.format("Picard %d < Pade %d", t1, t2));
  }

  @Test
  public void staticMethodsMatchInstanceMethods() {
    for (final double g : gain) {
      for (final double p : photons) {
        for (final double s : noise) {
          staticMethodsMatchInstanceMethods(g, p, s, true);
          staticMethodsMatchInstanceMethods(g, p, s, false);
        }
      }
    }
  }

  private static long getTime(double[] noise2, int N, double[][] x, final boolean usePicard) {
    // Warm up
    for (final double s2 : noise2) {
      for (int i = 0; i < photons.length; i++) {
        final double p = photons[i];
        for (int j = 0; j < N; j++) {
          PoissonGaussianFunction.probability(x[i][j], p, s2, usePicard);
        }
      }
    }

    // Time
    long t1 = System.nanoTime();
    for (final double s2 : noise2) {
      for (int i = 0; i < photons.length; i++) {
        final double p = photons[i];
        for (int j = 0; j < N; j++) {
          PoissonGaussianFunction.probability(x[i][j], p, s2, usePicard);
        }
      }
    }
    t1 = System.nanoTime() - t1;
    return t1;
  }

  private static void cumulativeProbabilityIsOne(final double gain, final double mu, final double s,
      final boolean usePicard) {
    final double p2 = cumulativeProbability(gain, mu, s, usePicard);
    Assertions.assertEquals(1, p2, 0.02, () -> String.format("g=%f, mu=%f, s=%f", gain, mu, s));
  }

  private static double cumulativeProbability(final double gain, final double mu, final double s,
      final boolean usePicard) {
    // Note: The input mu & s parameters are pre-gain.
    final PoissonGaussianFunction f =
        PoissonGaussianFunction.createWithStandardDeviation(1.0 / gain, mu, s * gain);
    f.setUsePicardApproximation(usePicard);
    double p = 0;
    int min = 1;
    int max = 0;

    // Evaluate an initial range.
    // Gaussian should have >99% within +/- s
    // Poisson will have mean mu with a variance mu.
    // At large mu it is approximately normal so use 3 sqrt(mu) for the range added to the mean
    if (mu > 0) {
      final int[] range = getRange(gain, mu, s);
      min = range[0];
      max = range[1];
      for (int x = min; x <= max; x++) {
        final double pp = f.probability(x);
        // logger.fine(FunctionUtils.getSupplier("x=%d, p=%f", x, pp);
        p += pp;
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
      p += pp;
      if (pp == 0 || pp / p < changeTolerance) {
        break;
      }
    }
    for (int x = max + 1;; x++) {
      max = x;
      final double pp = f.probability(x);
      // logger.fine(FunctionUtils.getSupplier("x=%d, p=%f", x, pp);
      p += pp;
      if (pp == 0 || pp / p < changeTolerance) {
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
      logger.log(
          TestLogUtils.getRecord(Level.INFO, "g=%f, mu=%f, s=%f p=%f  %f", gain, mu, s, p, p2));
    }

    return p2;
  }

  static int[] getRange(final double gain, final double mu, final double s) {
    // Evaluate an initial range.
    // Gaussian should have >99% within +/- s
    // Poisson will have mean mu with a variance mu.
    // At large mu it is approximately normal so use 3 sqrt(mu) for the range added to the mean
    final double range = Math.max(s, Math.sqrt(mu));
    final int min = (int) Math.floor(gain * (mu - 3 * range));
    final int max = (int) Math.ceil(gain * (mu + 3 * range));
    return new int[] {min, max};
  }

  private static void probabilityMatchesLogProbability(final double gain, final double mu,
      final double s, final boolean usePicard) {
    // Note: The input mu & s parameters are pre-gain.
    final PoissonGaussianFunction f =
        PoissonGaussianFunction.createWithStandardDeviation(1.0 / gain, mu, s * gain);
    f.setUsePicardApproximation(usePicard);

    // Evaluate an initial range.
    // Gaussian should have >99% within +/- s
    // Poisson will have mean mu with a variance mu.
    // At large mu it is approximately normal so use 3 sqrt(mu) for the range added to the mean
    final int[] range = getRange(gain, mu, s);
    final int min = range[0];
    final int max = range[1];
    final Supplier<String> msg = () -> String.format("g=%f, mu=%f, s=%f", gain, mu, s);
    final DoubleDoubleBiPredicate predicate = TestHelper.doublesAreClose(1e-3, 0);
    for (int x = min; x <= max; x++) {
      final double p = f.probability(x);
      if (p == 0) {
        continue;
      }
      final double logP = f.logProbability(x);
      TestAssertions.assertTest(Math.log(p), logP, predicate, msg);
    }
  }

  private static void staticMethodsMatchInstanceMethods(final double gain, final double mu,
      final double s, final boolean usePicard) {
    // Note: The input mu & s parameters are pre-gain.
    final PoissonGaussianFunction f =
        PoissonGaussianFunction.createWithStandardDeviation(1.0 / gain, mu, s * gain);
    f.setUsePicardApproximation(usePicard);

    // Evaluate an initial range.
    // Gaussian should have >99% within +/- s
    // Poisson will have mean mu with a variance mu.
    // At large mu it is approximately normal so use 3 sqrt(mu) for the range added to the mean
    final int[] range = getRange(gain, mu, s);
    final int min = range[0];
    final int max = range[1];
    final double logGain = Math.log(gain);
    final double s2 = MathUtils.pow2(s);
    final Supplier<String> msg1 = () -> String.format("probability g=%f, mu=%f, s=%f", gain, mu, s);
    final Supplier<String> msg2 =
        () -> String.format("logProbability g=%f, mu=%f, s=%f", gain, mu, s);
    for (int x = min; x <= max; x++) {
      double p = f.probability(x);
      double pp = PoissonGaussianFunction.probability(x / gain, mu, s2, usePicard) / gain;
      Assertions.assertEquals(p, pp, 1e-10, msg1);

      p = f.logProbability(x);
      pp = PoissonGaussianFunction.logProbability(x / gain, mu, s2, usePicard) - logGain;
      Assertions.assertEquals(p, pp, 1e-10, msg2);
    }
  }
}
