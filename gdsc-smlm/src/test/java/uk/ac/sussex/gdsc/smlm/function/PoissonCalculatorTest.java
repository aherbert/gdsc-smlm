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

import gnu.trove.list.array.TDoubleArrayList;
import java.math.BigDecimal;
import java.math.MathContext;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.IterativeLegendreGaussIntegrator;
import org.apache.commons.math3.analysis.integration.UnivariateIntegrator;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.opentest4j.AssertionFailedError;
import uk.ac.sussex.gdsc.core.data.DataException;
import uk.ac.sussex.gdsc.core.math.QuadraticUtils;
import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.smlm.GdscSmlmTestUtils;
import uk.ac.sussex.gdsc.smlm.utils.StdMath;
import uk.ac.sussex.gdsc.test.api.TestAssertions;
import uk.ac.sussex.gdsc.test.api.TestHelper;
import uk.ac.sussex.gdsc.test.api.function.DoubleDoubleBiPredicate;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.junit5.SpeedTag;
import uk.ac.sussex.gdsc.test.rng.RngUtils;
import uk.ac.sussex.gdsc.test.utils.BaseTimingTask;
import uk.ac.sussex.gdsc.test.utils.RandomSeed;
import uk.ac.sussex.gdsc.test.utils.TestComplexity;
import uk.ac.sussex.gdsc.test.utils.TestLogUtils;
import uk.ac.sussex.gdsc.test.utils.TestSettings;
import uk.ac.sussex.gdsc.test.utils.TimingService;

@SuppressWarnings({"javadoc"})
class PoissonCalculatorTest {
  /** The logging level for verbose test messages. */
  private static final Level LOG_LEVEL = Level.FINE;
  private static Logger logger;

  @BeforeAll
  public static void beforeAll() {
    logger = Logger.getLogger(PoissonCalculatorTest.class.getName());
  }

  @AfterAll
  public static void afterAll() {
    logger = null;
  }

  static double[] photons = {1, 1.5, 2, 2.5, 3, 4, 5, 7.5, 10, 100, 1000};
  private static int maxx = 10;

  static double P_LIMIT = 0.999999;

  @Test
  void canComputeLikelihoodForIntegerData() {
    final DoubleDoubleBiPredicate predicate = TestHelper.doublesAreClose(1e-10, 0);
    for (final double u : photons) {
      final PoissonDistribution pd = new PoissonDistribution(u);
      for (int x = 0; x < 100; x++) {
        double expected = pd.probability(x);
        double observed = PoissonCalculator.likelihood(u, x);
        if (expected > 1e-100) {
          TestAssertions.assertTest(expected, observed, predicate);
        }
        expected = pd.logProbability(x);
        observed = PoissonCalculator.logLikelihood(u, x);
        TestAssertions.assertTest(expected, observed, predicate);
      }
    }
  }

  @Test
  void canComputeFastLikelihoodForIntegerData() {
    final DoubleDoubleBiPredicate predicate = TestHelper.doublesAreClose(1e-4, 0);
    for (final double u : photons) {
      final PoissonDistribution pd = new PoissonDistribution(u);
      for (int x = 0; x < 100; x++) {
        double expected = pd.probability(x);
        double observed = PoissonCalculator.fastLikelihood(u, x);
        if (expected > 1e-100) {
          TestAssertions.assertTest(expected, observed, predicate);
        }
        expected = pd.logProbability(x);
        observed = PoissonCalculator.fastLogLikelihood(u, x);
        TestAssertions.assertTest(expected, observed, predicate);
      }
    }
  }

  @Test
  void canComputeFastLog_FastLikelihoodForIntegerData() {
    final DoubleDoubleBiPredicate predicate = TestHelper.doublesAreClose(1e-4, 0);
    final FastLog fastLog = FastLogFactory.getFastLog();
    for (final double u : photons) {
      final PoissonDistribution pd = new PoissonDistribution(u);
      for (int x = 0; x < 100; x++) {
        double expected = pd.probability(x);
        double observed = PoissonCalculator.fastLikelihood(u, x, fastLog);
        if (expected > 1e-100) {
          TestAssertions.assertTest(expected, observed, predicate);
        }
        expected = pd.logProbability(x);
        observed = PoissonCalculator.fastLogLikelihood(u, x, fastLog);
        TestAssertions.assertTest(expected, observed, predicate);
      }
    }
  }

  private abstract static class PoissonFunction implements UnivariateFunction {
    double mu;

    PoissonFunction(double mu) {
      this.mu = mu;
    }

    @Override
    public double value(double x) {
      double value = likelihood(mu, x);
      // v = pgf.probability(x, mu);
      // logger.fine(FunctionUtils.getSupplier("x=%f, v=%f", x, v);
      return value;
    }

    abstract double likelihood(double mu, double x);
  }

  @Test
  void likelihoodCumulativeProbabilityIsOneWithRealDataForCountAbove4() {
    cumulativeProbabilityIsOneWithRealDataForCountAbove4(0);
  }

  @Test
  void fastLikelihoodCumulativeProbabilityIsOneWithRealDataForCountAbove4() {
    cumulativeProbabilityIsOneWithRealDataForCountAbove4(1);
  }

  @Test
  void fastLog_fastLikelihoodCumulativeProbabilityIsNotOneWithRealDataForCountAbove4() {
    Assertions.assertThrows(AssertionFailedError.class, () -> {
      cumulativeProbabilityIsOneWithRealDataForCountAbove4(2);
    });
  }

  private static void cumulativeProbabilityIsOneWithRealDataForCountAbove4(int function) {
    for (final double mu : photons) {
      // Determine upper limit for a Poisson
      final double max = new PoissonDistribution(mu).inverseCumulativeProbability(P_LIMIT);

      // Determine lower limit
      final double sd = Math.sqrt(mu);
      final double min = (int) Math.max(0, mu - 4 * sd);

      PoissonFunction fun;
      switch (function) {
        case 2:
          fun = new PoissonFunction(mu) {
            @Override
            double likelihood(double mu, double x) {
              return PoissonCalculator.fastLikelihood(mu, x, FastLogFactory.getFastLog());
            }
          };
          break;
        case 1:
          fun = new PoissonFunction(mu) {
            @Override
            double likelihood(double mu, double x) {
              return PoissonCalculator.fastLikelihood(mu, x);
            }
          };
          break;
        case 0:
          fun = new PoissonFunction(mu) {
            @Override
            double likelihood(double mu, double x) {
              return PoissonCalculator.likelihood(mu, x);
            }
          };
          break;
        default:
          throw new IllegalStateException();
      }

      cumulativeProbabilityIsOneWithRealData(min, max, mu >= 4, fun);
    }
  }

  private static void cumulativeProbabilityIsOneWithRealData(double min, double max, boolean test,
      PoissonFunction function) {
    final int integrationpoints = 10;
    final double relativeAccuracy = 1e-4;
    final double absoluteAccuracy = 1e-8;
    final int minimalIterationCount = 3;
    final int maximalIterationCount = 32;

    final UnivariateIntegrator in = new IterativeLegendreGaussIntegrator(integrationpoints,
        relativeAccuracy, absoluteAccuracy, minimalIterationCount, maximalIterationCount);
    // new SimpsonIntegrator();

    final double p = in.integrate(20000, function, min, max);

    logger.log(TestLogUtils.getRecord(LOG_LEVEL, "mu=%f, p=%f", function.mu, p));
    if (test) {
      Assertions.assertEquals(P_LIMIT, p, 0.02, () -> "mu=" + function.mu);
    }
  }

  private abstract static class BaseNonLinearFunction implements NonLinearFunction {
    double[] params;
    String name;

    BaseNonLinearFunction(String name) {
      this.name = name;
    }

    @Override
    public void initialise(double[] params) {
      this.params = params;
    }

    @Override
    public int[] gradientIndices() {
      return new int[1];
    }

    @Override
    public double evalw(int x, double[] dyda, double[] weights) {
      return 0;
    }

    @Override
    public double evalw(int x, double[] weights) {
      return 0;
    }

    @Override
    public double eval(int x, double[] dyda) {
      return 0;
    }

    @Override
    public boolean canComputeWeights() {
      return false;
    }

    @Override
    public int getNumberOfGradients() {
      return 1;
    }
  }

  @SeededTest
  void canComputeLogLikelihoodRatio(RandomSeed seed) {
    final double n2 = maxx * maxx * 0.5;
    // Functions must produce a strictly positive output so add background
    canComputeLogLikelihoodRatio(seed, new BaseNonLinearFunction("Quadratic") {
      @Override
      public double eval(int x) {
        return 0.1 + params[0] * (x - n2) * (x - n2);
      }
    });
    canComputeLogLikelihoodRatio(seed, new BaseNonLinearFunction("Gaussian") {
      @Override
      public double eval(int x) {
        return 0.1 + 100 * StdMath.exp(-0.5 * MathUtils.pow2(x - n2) / (params[0] * params[0]));
      }
    });
  }

  private static void canComputeLogLikelihoodRatio(RandomSeed seed, BaseNonLinearFunction nlf) {
    logger.log(TestLogUtils.getRecord(LOG_LEVEL, nlf.name));

    final int n = maxx * maxx;

    final double[] a = new double[] {1};

    // Simulate Poisson process
    nlf.initialise(a);
    final UniformRandomProvider rng = RngUtils.create(seed.get());
    final double[] x = new double[n];
    final double[] u = new double[n];
    for (int i = 0; i < n; i++) {
      u[i] = nlf.eval(i);
      if (u[i] > 0) {
        x[i] = GdscSmlmTestUtils.createPoissonSampler(rng, u[i]).sample();
      }
    }

    final DoubleDoubleBiPredicate predicate = TestHelper.doublesAreClose(1e-10, 0);

    double ll = PoissonCalculator.logLikelihood(u, x);
    final double mll = PoissonCalculator.maximumLogLikelihood(x);
    double llr = -2 * (ll - mll);
    double llr2 = PoissonCalculator.logLikelihoodRatio(u, x);
    logger.log(TestLogUtils.getRecord(LOG_LEVEL, "llr=%f, llr2=%f", llr, llr2));
    TestAssertions.assertTest(llr, llr2, predicate, "Log-likelihood ratio");

    final double[] op = new double[x.length];
    for (int i = 0; i < n; i++) {
      op[i] = PoissonCalculator.maximumLikelihood(x[i]);
    }

    // TestSettings.setLogLevel(uk.ac.sussex.gdsc.smlm.TestLog.Level.FINE);

    final int df = n - 1;
    final ChiSquaredDistributionTable table =
        ChiSquaredDistributionTable.createUpperTailed(0.05, df);
    final ChiSquaredDistributionTable table2 =
        ChiSquaredDistributionTable.createUpperTailed(0.001, df);
    if (logger.isLoggable(LOG_LEVEL)) {
      logger.log(TestLogUtils.getRecord(LOG_LEVEL, "Chi2 = %f (q=%.3f), %f (q=%.3f)  %f %b  %f",
          table.getCrititalValue(df), table.getSignificanceValue(), table2.getCrititalValue(df),
          table2.getSignificanceValue(), ChiSquaredDistributionTable.computeQValue(24, 2),
          ChiSquaredDistributionTable.createUpperTailed(0.05, 2).reject(24, 2),
          ChiSquaredDistributionTable.getChiSquared(1e-6, 2)

      ));
    }
    final TDoubleArrayList list = new TDoubleArrayList();
    final int imin = 5;
    final int imax = 15;
    for (int i = imin; i <= imax; i++) {
      a[0] = (double) i / 10;
      nlf.initialise(a);
      for (int j = 0; j < n; j++) {
        u[j] = nlf.eval(j);
      }

      ll = PoissonCalculator.logLikelihood(u, x);
      list.add(ll);
      llr = PoissonCalculator.logLikelihoodRatio(u, x);
      BigDecimal product = new BigDecimal(1);
      double ll2 = 0;
      for (int j = 0; j < n; j++) {
        final double p1 = PoissonCalculator.likelihood(u[j], x[j]);
        ll2 += Math.log(p1);
        final double ratio = p1 / op[j];
        product = product.multiply(new BigDecimal(ratio));
      }
      llr2 = -2 * Math.log(product.doubleValue());
      final double p = ChiSquaredDistributionTable.computePValue(llr, df);
      final double q = ChiSquaredDistributionTable.computeQValue(llr, df);
      logger.log(TestLogUtils.getRecord(LOG_LEVEL,
          "a=%f, ll=%f, ll2=%f, llr=%f, llr2=%f, product=%s, p=%f, q=%f "
              + "(reject=%b @ %.3f, reject=%b @ %.3f)",
          a[0], ll, ll2, llr, llr2, product.round(new MathContext(4)).toString(), p, q,
          table.reject(llr, df), table.getSignificanceValue(), table2.reject(llr, df),
          table2.getSignificanceValue()));

      // Only test if the product could be computed. Low ratios cause it to become
      // too small to store in a double.
      if (product.doubleValue() > 0) {
        TestAssertions.assertTest(ll, ll2, predicate, "Log-likelihood");
        TestAssertions.assertTest(llr, llr2, predicate, "Log-likelihood ratio");
      }
    }

    // Find max using quadratic fit
    final double[] data = list.toArray();
    int index = SimpleArrayUtils.findMaxIndex(data);
    final double maxa = (double) (imin + index) / 10;
    double fita = maxa;
    try {
      if (index == 0) {
        index++;
      }
      if (index == data.length - 1) {
        index--;
      }
      final int i1 = index - 1;
      final int i2 = index;
      final int i3 = index + 1;

      fita = QuadraticUtils.findMinMax((double) (imin + i1) / 10, data[i1],
          (double) (imin + i2) / 10, data[i2], (double) (imin + i3) / 10, data[i3]);
    } catch (final DataException ex) {
      // Ignore
    }

    // Allow a tolerance as the random data may alter the p-value computation.
    // Should allow it to be less than 2 increment either side of the answer.
    logger.log(TestLogUtils.getRecord(LOG_LEVEL, "max fit = %g => %g", maxa, fita));
    Assertions.assertEquals(1, fita, 0.199, "max");
  }

  @SeededTest
  void canComputeFastLog_LogLikelihoodRatio(RandomSeed seed) {
    final double n2 = maxx * maxx * 0.5;
    // Functions must produce a strictly positive output so add background
    canComputeFastLog_LogLikelihoodRatio(seed, new BaseNonLinearFunction("Quadratic") {
      @Override
      public double eval(int x) {
        return 0.1 + params[0] * (x - n2) * (x - n2);
      }
    });
    canComputeFastLog_LogLikelihoodRatio(seed, new BaseNonLinearFunction("Gaussian") {
      @Override
      public double eval(int x) {
        return 0.1 + 100 * StdMath.exp(-0.5 * MathUtils.pow2(x - n2) / (params[0] * params[0]));
      }
    });
  }

  private static void canComputeFastLog_LogLikelihoodRatio(RandomSeed seed,
      BaseNonLinearFunction nlf) {
    logger.log(TestLogUtils.getRecord(LOG_LEVEL, nlf.name));

    final int n = maxx * maxx;

    final double[] a = new double[] {1};

    // Simulate Poisson process
    nlf.initialise(a);
    final UniformRandomProvider rng = RngUtils.create(seed.get());
    final double[] x = new double[n];
    final double[] u = new double[n];
    for (int i = 0; i < n; i++) {
      u[i] = nlf.eval(i);
      if (u[i] > 0) {
        x[i] = GdscSmlmTestUtils.createPoissonSampler(rng, u[i]).sample();
      }
    }

    // Only test the LLR
    final double llr = PoissonCalculator.logLikelihoodRatio(u, x);
    final double llr2 = PoissonCalculator.logLikelihoodRatio(u, x, FastLogFactory.getFastLog());
    logger.log(TestLogUtils.getRecord(LOG_LEVEL, "llr=%f, llr2=%f", llr, llr2));
    // Approximately equal
    TestAssertions.assertTest(llr, llr2, TestHelper.doublesAreClose(5e-3, 0),
        "Log-likelihood ratio");
  }

  @SeededTest
  void cannotSubtractConstantBackgroundAndComputeLogLikelihoodRatio(RandomSeed seed) {
    final int n = maxx * maxx;
    final double n2 = n * 0.5;
    final double n3 = n * 0.33;
    final double n4 = n * 0.66;
    // Functions must produce a strictly positive output so add background
    cannotSubtractConstantBackgroundAndComputeLogLikelihoodRatio(seed,
        new BaseNonLinearFunction("Quadratic") {
          @Override
          public double eval(int x) {
            return 0.1 + params[0] * (x - n2) * (x - n2);
          }
        }, new BaseNonLinearFunction("Quadratic") {
          @Override
          public double eval(int x) {
            return 0.2 + 0.5 * params[0] * (x - n3) * (x - n3);
          }
        }, new BaseNonLinearFunction("Quadratic") {
          @Override
          public double eval(int x) {
            return 0.3 + 0.75 * params[0] * (x - n4) * (x - n4);
          }
        });
    cannotSubtractConstantBackgroundAndComputeLogLikelihoodRatio(seed,
        new BaseNonLinearFunction("Gaussian") {
          @Override
          public double eval(int x) {
            return 0.1 + 100 * StdMath.exp(-0.5 * MathUtils.pow2(x - n2) / (params[0] * params[0]));
          }
        }, new BaseNonLinearFunction("Gaussian") {
          @Override
          public double eval(int x) {
            return 0.2 + 50 * StdMath.exp(-0.5 * MathUtils.pow2(x - n3) / (params[0] * params[0]));
          }
        }, new BaseNonLinearFunction("Gaussian") {
          @Override
          public double eval(int x) {
            return 0.3 + 75 * StdMath.exp(-0.5 * MathUtils.pow2(x - n4) / (params[0] * params[0]));
          }
        });
  }

  private static void cannotSubtractConstantBackgroundAndComputeLogLikelihoodRatio(RandomSeed seed,
      BaseNonLinearFunction nlf1, BaseNonLinearFunction nlf2, BaseNonLinearFunction nlf3) {
    // System.out.println(nlf1.name);

    final int n = maxx * maxx;

    final double[] a = new double[] {1};

    // Simulate Poisson process of 3 combined functions
    nlf1.initialise(a);
    nlf2.initialise(a);
    nlf3.initialise(a);
    final UniformRandomProvider rng = RngUtils.create(seed.get());
    double[] x = SimpleArrayUtils.newArray(n, 0, 1.0);
    final double[] u = new double[x.length];
    final double[] b1 = new double[x.length];
    final double[] b2 = new double[x.length];
    final double[] b3 = new double[x.length];
    for (int i = 0; i < n; i++) {
      b1[i] = nlf1.eval(i);
      b2[i] = nlf2.eval(i);
      b3[i] = nlf3.eval(i);
      u[i] = b1[i] + b2[i] + b3[i];
      if (u[i] > 0) {
        x[i] = GdscSmlmTestUtils.createPoissonSampler(rng, u[i]).sample();
      }
    }

    // x is the target data
    // b1 is the already computed background
    // b2 is the first function to add to the model
    // b3 is the second function to add to the model

    // Adjust x to ensure that x - background is positive.
    double[] xb = subtract(x, b1);
    SimpleArrayUtils.apply(xb, z -> z < 0 ? 0 : z);
    x = add(xb, b1);

    // Compute the LLR of adding b3 to b2 given we already have b1 modelling data x
    final double[] b12 = add(b1, b2);
    final double ll1a = PoissonCalculator.logLikelihood(b12, x);
    final double ll2a = PoissonCalculator.logLikelihood(add(b12, b3), x);
    final double llra = -2 * (ll1a - ll2a);
    // logger.fine(FunctionUtils.getSupplier("x|(a+b+c) ll1=%f, ll2=%f, llra=%f", ll1a, ll2a, llra);

    // Compute the LLR of adding b3 to b2 given we already have x minus b1
    final double ll1b = PoissonCalculator.logLikelihood(b2, xb);
    final double ll2b = PoissonCalculator.logLikelihood(add(b2, b3), xb);
    final double llrb = -2 * (ll1b - ll2b);
    // logger.fine(FunctionUtils.getSupplier("x-a|(b+c) : ll1=%f, ll2=%f, llrb=%f", ll1b, ll2b,
    // llrb);

    // logger.fine(FunctionUtils.getSupplier("llr=%f (%g), llr2=%f (%g)", llra,
    // PoissonCalculator.computePValue(llra, 1), llrb,
    // PoissonCalculator.computePValue(llrb, 1));
    Assertions.assertThrows(AssertionError.class, () -> {
      TestAssertions.assertTest(llra, llrb, TestHelper.doublesAreClose(1e-10, 0),
          "Log-likelihood ratio");
    });
  }

  @Test
  void showRelativeErrorOfLogFactorialApproximation() {
    Assumptions.assumeTrue(logger.isLoggable(LOG_LEVEL));
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.HIGH));

    // The approximation is very bad close to 1
    double value = 1.0;
    for (int i = 1; i <= 100; i++) {
      value = Math.nextUp(value);
      showRelativeErrorOfLogFactorialApproximation(value);
    }

    // Error for the first LogFactorial computation switch point at x=1.5.
    // Below this level the error is very bad.
    // INFO: 1.5! = 0.2847 : [1.3192219378603063, 0.1925445003447508, 0.002597657209671918,
    // 2.866725480047198E-4, 8.04441744746754E-5, 4.192375182790675E-5]
    for (int i = 1; i <= 300; i++) {
      showRelativeErrorOfLogFactorialApproximation(1 + i / 100.0);
    }
    // INFO: 4.0! = 3.178 : [0.5137975858922471, 0.006541950896246824, 1.3423516145952628E-5,
    // 2.3333406200452226E-7, 1.0541480759747081E-8, 8.90185218881136E-10]

    // Error approaches machine epsilon using 3 terms:
    // INFO: 100.0! = 363.7 : [0.008858972036867302, 2.2910100241207718E-6, 7.63669686628571E-12,
    // 0.0, 1.5627513181378083E-16, 1.5627513181378083E-16]
    for (int i = 4; i <= 100; i++) {
      showRelativeErrorOfLogFactorialApproximation(i);
    }
  }

  private static void showRelativeErrorOfLogFactorialApproximation(double x) {
    final double e = LogFactorial.value(x);
    final double[] o = new double[6];
    final double[] error = new double[o.length];
    for (int i = 0; i < o.length; i++) {
      o[i] = PoissonCalculator.logFactorialApproximation(x, i);
      error[i] = DoubleEquality.relativeError(e, o[i]);
    }
    logger.log(TestLogUtils.getRecord(LOG_LEVEL, "%s! = %s : %s", Double.toString(x),
        MathUtils.rounded(e), Arrays.toString(error)));
  }

  @Test
  void showRelativeErrorOfFastLogLikelihood() {
    Assumptions.assumeTrue(logger.isLoggable(LOG_LEVEL));
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.HIGH));

    double value = 1.0;
    for (int i = 1; i <= 100; i++) {
      value = Math.nextUp(value);
      showRelativeErrorOfFastLogLikelihood(value);
    }
    for (int i = 1; i <= 300; i++) {
      showRelativeErrorOfFastLogLikelihood(1 + i / 100.0);
    }

    for (int i = 4; i <= 100; i++) {
      showRelativeErrorOfFastLogLikelihood(i);
    }
  }

  private static void showRelativeErrorOfFastLogLikelihood(double x) {
    for (final double factor : new double[] {0.5, 1, 2}) {
      final double u = x * factor;
      final double e = PoissonCalculator.logLikelihood(u, x);
      final double o = PoissonCalculator.fastLogLikelihood(u, x);
      final double error = DoubleEquality.relativeError(e, o);
      logger.log(TestLogUtils.getRecord(LOG_LEVEL, "ll(%s|%s) = %s : %s", Double.toString(x),
          Double.toString(u), MathUtils.rounded(e), Double.toString(error)));
    }
  }

  @Test
  void showRelativeErrorOfFastLog_FastLogLikelihood() {
    Assumptions.assumeTrue(logger.isLoggable(LOG_LEVEL));
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.HIGH));

    double value = 1.0;
    for (int i = 1; i <= 100; i++) {
      value = Math.nextUp(value);
      showRelativeErrorOfFastLog_FastLogLikelihood(value);
    }
    for (int i = 1; i <= 300; i++) {
      showRelativeErrorOfFastLog_FastLogLikelihood(1 + i / 100.0);
    }

    for (int i = 4; i <= 100; i++) {
      showRelativeErrorOfFastLog_FastLogLikelihood(i);
    }
  }

  private static void showRelativeErrorOfFastLog_FastLogLikelihood(double x) {
    final FastLog fastLog = FastLogFactory.getFastLog();
    for (final double factor : new double[] {0.5, 1, 2}) {
      final double u = x * factor;
      final double e = PoissonCalculator.logLikelihood(u, x);
      final double o = PoissonCalculator.fastLogLikelihood(u, x, fastLog);
      final double error = DoubleEquality.relativeError(e, o);
      logger.log(TestLogUtils.getRecord(LOG_LEVEL, "ll(%s|%s) = %s : %s", Double.toString(x),
          Double.toString(u), MathUtils.rounded(e), Double.toString(error)));
    }
  }

  @Test
  void showRelativeErrorOfFastLog_LogLikelihoodRatio() {
    Assumptions.assumeTrue(logger.isLoggable(LOG_LEVEL));
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.HIGH));

    double value = 1.0;
    for (int i = 1; i <= 100; i++) {
      value = Math.nextUp(value);
      showRelativeErrorOfFastLog_LogLikelihoodRatio(value);
    }
    for (int i = 1; i <= 300; i++) {
      showRelativeErrorOfFastLog_LogLikelihoodRatio(1 + i / 100.0);
    }

    for (int i = 4; i <= 100; i++) {
      showRelativeErrorOfFastLog_LogLikelihoodRatio(i);
    }
  }

  private static void showRelativeErrorOfFastLog_LogLikelihoodRatio(double x) {
    final FastLog fastLog = FastLogFactory.getFastLog();
    for (final double factor : new double[] {0.5, 1, 2}) {
      final double u = x * factor;
      final double e = PoissonCalculator.logLikelihoodRatio(u, x);
      final double o = PoissonCalculator.logLikelihoodRatio(u, x, fastLog);
      final double error = DoubleEquality.relativeError(e, o);
      logger.log(TestLogUtils.getRecord(LOG_LEVEL, "llr(%s|%s) = %s : %s", Double.toString(x),
          Double.toString(u), MathUtils.rounded(e), Double.toString(error)));
    }
  }

  @SeededTest
  void instanceAndFastMethodIsApproximatelyEqualToStaticMethod(RandomSeed seed) {
    final DoubleEquality eq = new DoubleEquality(3e-4, 0);
    final UniformRandomProvider rg = RngUtils.create(seed.get());
    // Test for different x. The calculator approximation begins
    final int n = 100;
    final double[] u = new double[n];
    final double[] x = new double[n];
    double expected;
    double observed;
    for (final double testx : new double[] {Math.nextDown(PoissonCalculator.APPROXIMATION_X),
        PoissonCalculator.APPROXIMATION_X, Math.nextUp(PoissonCalculator.APPROXIMATION_X),
        PoissonCalculator.APPROXIMATION_X * 1.1, PoissonCalculator.APPROXIMATION_X * 2,
        PoissonCalculator.APPROXIMATION_X * 10}) {
      final String X = Double.toString(testx);
      Arrays.fill(x, testx);
      final PoissonCalculator pc = new PoissonCalculator(x);
      expected = PoissonCalculator.maximumLogLikelihood(x);
      observed = pc.getMaximumLogLikelihood();
      logger.log(TestLogUtils.getRecord(LOG_LEVEL, "[%s] Instance MaxLL = %g vs %g (error = %g)", X,
          expected, observed, DoubleEquality.relativeError(expected, observed)));
      Assertions.assertTrue(eq.almostEqualRelativeOrAbsolute(expected, observed),
          () -> "Instance Max LL not equal: x=" + X);

      observed = PoissonCalculator.fastMaximumLogLikelihood(x);
      logger.log(TestLogUtils.getRecord(LOG_LEVEL, "[%s] Fast MaxLL = %g vs %g (error = %g)", X,
          expected, observed, DoubleEquality.relativeError(expected, observed)));
      Assertions.assertTrue(eq.almostEqualRelativeOrAbsolute(expected, observed),
          () -> "Fast Max LL not equal: x=" + X);

      // Generate data around the value
      for (int i = 0; i < n; i++) {
        u[i] = x[i] + rg.nextDouble() - 0.5;
      }

      expected = PoissonCalculator.logLikelihood(u, x);
      observed = pc.logLikelihood(u);
      logger.log(TestLogUtils.getRecord(LOG_LEVEL, "[%s] Instance LL = %g vs %g (error = %g)", X,
          expected, observed, DoubleEquality.relativeError(expected, observed)));
      Assertions.assertTrue(eq.almostEqualRelativeOrAbsolute(expected, observed),
          () -> "Instance LL not equal: x=" + X);

      observed = PoissonCalculator.fastLogLikelihood(u, x);
      logger.log(TestLogUtils.getRecord(LOG_LEVEL, "[%s] Fast LL = %g vs %g (error = %g)", X,
          expected, observed, DoubleEquality.relativeError(expected, observed)));
      Assertions.assertTrue(eq.almostEqualRelativeOrAbsolute(expected, observed),
          () -> "Fast LL not equal: x=" + X);

      expected = PoissonCalculator.logLikelihoodRatio(u, x);
      observed = pc.getLogLikelihoodRatio(observed);

      logger.log(TestLogUtils.getRecord(LOG_LEVEL, "[%s] Instance LLR = %g vs %g (error = %g)", X,
          expected, observed, DoubleEquality.relativeError(expected, observed)));
      Assertions.assertTrue(eq.almostEqualRelativeOrAbsolute(expected, observed),
          () -> "Instance LLR not equal: x=" + X);
    }
  }

  private abstract static class PcTimingTask extends BaseTimingTask {
    double[] x;
    double[] mean;
    int ll;
    int llr;

    public PcTimingTask(String name, double[] x, double[] mean, int ll, int llr) {
      super(String.format("%s ll=%d llr=%d", name, ll, llr));
      this.x = x;
      this.mean = mean;
      this.ll = ll;
      this.llr = llr;
    }

    @Override
    public int getSize() {
      return 1;
    }

    @Override
    public Object getData(int index) {
      return null;
    }
  }

  private static class StaticPcTimingTask extends PcTimingTask {
    public StaticPcTimingTask(double[] x, double[] mean, int ll, int llr) {
      super("static", x, mean, ll, llr);
    }

    @Override
    public Object run(Object data) {
      double value = 0;
      for (int i = 0; i < llr; i++) {
        value += PoissonCalculator.logLikelihoodRatio(mean, x);
      }
      for (int i = 0; i < ll; i++) {
        value += PoissonCalculator.logLikelihood(mean, x);
      }
      return value;
    }
  }

  private static class FastPcTimingTask extends PcTimingTask {
    public FastPcTimingTask(double[] x, double[] mean, int ll, int llr) {
      super("fast", x, mean, ll, llr);
    }

    @Override
    public Object run(Object data) {
      double value = 0;
      for (int i = 0; i < llr; i++) {
        value += PoissonCalculator.logLikelihoodRatio(mean, x);
      }
      for (int i = 0; i < ll; i++) {
        value += PoissonCalculator.fastLogLikelihood(mean, x);
      }
      return value;
    }
  }

  private static class FastLogPcTimingTask extends PcTimingTask {
    FastLog fastLog = FastLogFactory.getFastLog();

    public FastLogPcTimingTask(double[] x, double[] mean, int ll, int llr) {
      super("fastLog", x, mean, ll, llr);
    }

    @Override
    public Object run(Object data) {
      double value = 0;
      for (int i = 0; i < llr; i++) {
        value += PoissonCalculator.logLikelihoodRatio(mean, x, fastLog);
      }
      for (int i = 0; i < ll; i++) {
        value += PoissonCalculator.fastLogLikelihood(mean, x, fastLog);
      }
      return value;
    }
  }

  private static class InstancePcTimingTask extends PcTimingTask {
    int max;

    public InstancePcTimingTask(double[] x, double[] mean, int ll, int llr) {
      super("instance", x, mean, ll, llr);
      max = Math.max(llr, ll);
    }

    @Override
    public Object run(Object data) {
      final PoissonCalculator pc = new PoissonCalculator(x);
      double value = 0;
      // Use the fastest execution possible
      for (int i = 0; i < max; i++) {
        value += pc.pseudoLogLikelihood(mean);
      }
      if (llr > 0) {
        value += pc.getMaximumLogLikelihood();
      }
      return value;
    }
  }

  @SpeedTag
  @Test
  void instanceMethodIsFaster() {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final int n = 1000;
    final int m = 10;
    final double[] x = new double[n * m];
    final double[] u = new double[x.length];
    for (int i = 1, k = 0; i <= n; i++) {
      final double testx = 0.1 * i;
      // +/- 3SD of the expected
      final double sd = 3 * Math.sqrt(testx);
      final double min = Math.max(0.1, testx - sd);
      final double max = testx + sd;
      final double inc = (max - min) / (m - 1);
      for (int j = 0; j < m; j++, k++) {
        x[k] = testx;
        u[k] = min + j * inc;
      }
    }
    final double[] limits = MathUtils.limits(x);
    logger.log(
        TestLogUtils.getRecord(LOG_LEVEL, "Speed test x-range: %f - %f", limits[0], limits[1]));

    final TimingService ts = new TimingService(5);
    final int[] loops = new int[] {0, 1, 10};
    for (final int ll : loops) {
      for (final int llr : loops) {
        if (ll + llr == 0) {
          continue;
        }
        ts.execute(new StaticPcTimingTask(x, u, ll, llr));
        ts.execute(new FastPcTimingTask(x, u, ll, llr));
        ts.execute(new FastLogPcTimingTask(x, u, ll, llr));
        ts.execute(new InstancePcTimingTask(x, u, ll, llr));
      }
    }

    final int size = ts.getSize();
    ts.repeat(size);
    if (logger.isLoggable(LOG_LEVEL)) {
      logger.info(ts.getReport(size));
    }

    final int index = ts.getSize() - 1;
    Assertions.assertTrue(ts.get(index).getMean() < ts.get(index - 1).getMean());
  }

  private static double[] add(double[] v1, double[] v2) {
    final double[] result = new double[v1.length];
    for (int i = 0; i < v1.length; i++) {
      result[i] = v1[i] + v2[i];
    }
    return result;
  }

  private static double[] subtract(double[] v1, double[] v2) {
    final double[] result = new double[v1.length];
    for (int i = 0; i < v1.length; i++) {
      result[i] = v1[i] - v2[i];
    }
    return result;
  }
}
