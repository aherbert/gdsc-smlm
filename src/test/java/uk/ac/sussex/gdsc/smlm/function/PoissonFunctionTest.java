package uk.ac.sussex.gdsc.smlm.function;

import uk.ac.sussex.gdsc.smlm.math3.distribution.CustomPoissonDistribution;
import uk.ac.sussex.gdsc.test.api.TestAssertions;
import uk.ac.sussex.gdsc.test.api.TestHelper;
import uk.ac.sussex.gdsc.test.api.function.DoubleDoubleBiPredicate;
import uk.ac.sussex.gdsc.test.api.function.DoublePredicate;
import uk.ac.sussex.gdsc.test.utils.TestLogUtils;
import uk.ac.sussex.gdsc.test.utils.functions.FunctionUtils;

import gnu.trove.list.array.TDoubleArrayList;

import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import java.util.logging.Level;
import java.util.logging.Logger;

@SuppressWarnings({"unused", "javadoc"})
public class PoissonFunctionTest {
  private static Logger logger;

  @BeforeAll
  public static void beforeAll() {
    logger = Logger.getLogger(PoissonFunctionTest.class.getName());
  }

  @AfterAll
  public static void afterAll() {
    logger = null;
  }

  static double[] gain = {0.25, 0.5, 0.7, 1, 1.5, 1.7, 2, 2.2, 4, 8, 16};
  static double[] photons = {0.25, 0.5, 1, 2, 4, 10, 100, 1000};

  @Test
  public void cumulativeProbabilityIsOne() {
    for (int j = 0; j < gain.length; j++) {
      for (int i = 0; i < photons.length; i++) {
        final int[] result = cumulativeProbabilityIsOne(gain[j], photons[i]);
        logger.log(TestLogUtils.getRecord(Level.FINE, "minRange[%d][%d] = %d;", j, i, result[0]));
        logger.log(TestLogUtils.getRecord(Level.FINE, "maxRange[%d][%d] = %d;", j, i, result[1]));
      }
    }
  }

  private static int[] cumulativeProbabilityIsOne(final double gain, final double mu) {
    final double o = mu;

    final PoissonFunction f = new PoissonFunction(1.0 / gain);
    double p = 0;

    final TDoubleArrayList values = new TDoubleArrayList();

    double maxp = 0;
    int maxc = 0;

    // Evaluate an initial range.
    // Poisson will have mean mu with a variance mu.
    // At large mu it is approximately normal so use 3 sqrt(mu) for the range added to the mean

    final int[] range = getRange(gain, mu);
    int min = range[0];
    int max = range[1];
    for (int x = min; x <= max; x++) {
      final double pp = f.likelihood(x, o);
      // logger.fine(FunctionUtils.getSupplier("x=%d, p=%f", x, pp);
      p += pp;
      values.add(pp);
      if (maxp < pp) {
        maxp = pp;
        maxc = x;
      }
    }
    if (p > 1.01) {
      Assertions.fail("P > 1: " + p);
    }

    // We have most of the probability density.
    // Now keep evaluating up and down until no difference
    final double changeTolerance = 1e-6;
    if (min > 0) {
      values.reverse();
      for (int x = min - 1; x >= 0; x--) {
        min = x;
        final double pp = f.likelihood(x, o);
        // logger.fine(FunctionUtils.getSupplier("x=%d, p=%f", x, pp);
        p += pp;
        values.add(pp);
        if (maxp < pp) {
          maxp = pp;
          maxc = x;
        }
        if (pp == 0 || pp / p < changeTolerance) {
          break;
        }
      }
      values.reverse();
    }
    for (int x = max + 1;; x++) {
      max = x;
      final double pp = f.likelihood(x, o);
      // logger.fine(FunctionUtils.getSupplier("x=%d, p=%f", x, pp);
      p += pp;
      values.add(pp);
      if (maxp < pp) {
        maxp = pp;
        maxc = x;
      }
      if (pp == 0 || pp / p < changeTolerance) {
        break;
      }
    }

    // Find the range for 99.5% of the sum
    final double[] h = values.toArray();
    // Find cumulative
    for (int i = 1; i < h.length; i++) {
      h[i] += h[i - 1];
    }
    int minx = 0, maxx = h.length - 1;
    while (h[minx + 1] < 0.0025) {
      minx++;
    }
    while (h[maxx - 1] > 0.9975) {
      maxx--;
    }

    minx += min;
    maxx += min;

    logger
        .log(TestLogUtils.getRecord(Level.INFO, "g=%f, mu=%f, o=%f, p=%f, min=%d, %f @ %d, max=%d",
            gain, mu, o, p, minx, maxp, maxc, maxx));
    return new int[] {minx, maxx};
  }

  static int[] getRange(final double gain, final double mu) {
    // Evaluate an initial range.
    // Poisson will have mean mu with a variance mu.
    // At large mu it is approximately normal so use 3 sqrt(mu) for the range added to the mean
    final double range = Math.max(1, Math.sqrt(mu));
    final int min = Math.max(0, (int) Math.floor(gain * (mu - 3 * range)));
    final int max = (int) Math.ceil(gain * (mu + 3 * range));
    return new int[] {min, max};
  }

  @Test
  public void probabilityMatchesLogProbabilty() {
    for (int j = 0; j < gain.length; j++) {
      for (int i = 0; i < photons.length; i++) {
        probabilityMatchesLogProbabilty(gain[j], photons[i]);
      }
    }
  }

  private static void probabilityMatchesLogProbabilty(final double gain, final double mu) {
    final double o = mu;

    final PoissonFunction f = new PoissonFunction(1.0 / gain);
    final double p = 0;

    final int[] range = getRange(gain, mu);
    final int min = range[0];
    final int max = range[1];
    // Allow a relative difference or an exact match (for the case of -Infinity)
    final DoublePredicate isNegativeInfinity = TestHelper.doubleEquals(Double.NEGATIVE_INFINITY);
    final DoubleDoubleBiPredicate predicate = TestHelper.doublesAreClose(1e-8, 0)
        .or(TestHelper.and(isNegativeInfinity, isNegativeInfinity));
    for (int x = min; x <= max; x++) {
      final double v1 = Math.log(f.likelihood(x, o));
      final double v2 = f.logLikelihood(x, o);
      // System.out.printf("x=%d, v1=%s, v2=%s%n", x, v1, v2);

      TestAssertions.assertTest(v1, v2, predicate,
          FunctionUtils.getSupplier("g=%f, mu=%f, x=%d", gain, mu, x));
    }
  }

  @Test
  public void probabilityMatchesPoissonWithNoGain() {
    for (int i = 0; i < photons.length; i++) {
      probabilityMatchesPoissonWithNoGain(photons[i]);
    }
  }

  private static void probabilityMatchesPoissonWithNoGain(final double mu) {
    final double o = mu;

    final PoissonFunction f = new PoissonFunction(1.0);
    final CustomPoissonDistribution pd = new CustomPoissonDistribution(null, mu);

    final double p = 0;

    final int[] range = getRange(1, mu);
    final int min = range[0];
    final int max = range[1];
    final DoubleDoubleBiPredicate predicate = TestHelper.doublesAreClose(1e-8, 0);
    for (int x = min; x <= max; x++) {
      final double v1 = f.likelihood(x, o);
      final double v2 = pd.probability(x);

      TestAssertions.assertTest(v1, v2, predicate,
          FunctionUtils.getSupplier("g=%f, mu=%f, x=%d", gain, mu, x));
    }
  }
}
