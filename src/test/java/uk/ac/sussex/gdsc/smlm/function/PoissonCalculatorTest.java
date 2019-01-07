/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2019 Alex Herbert
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

import uk.ac.sussex.gdsc.core.data.DataException;
import uk.ac.sussex.gdsc.core.math.QuadraticUtils;
import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.RandomGeneratorAdapter;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.smlm.math3.distribution.CustomPoissonDistribution;
import uk.ac.sussex.gdsc.test.api.TestAssertions;
import uk.ac.sussex.gdsc.test.api.TestHelper;
import uk.ac.sussex.gdsc.test.api.function.DoubleDoubleBiPredicate;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.junit5.SpeedTag;
import uk.ac.sussex.gdsc.test.rng.RngUtils;
import uk.ac.sussex.gdsc.test.utils.BaseTimingTask;
import uk.ac.sussex.gdsc.test.utils.TestComplexity;
import uk.ac.sussex.gdsc.test.utils.TestLogUtils;
import uk.ac.sussex.gdsc.test.utils.TestSettings;
import uk.ac.sussex.gdsc.test.utils.TimingService;

import gnu.trove.list.array.TDoubleArrayList;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.IterativeLegendreGaussIntegrator;
import org.apache.commons.math3.analysis.integration.UnivariateIntegrator;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.opentest4j.AssertionFailedError;

import java.math.BigDecimal;
import java.math.MathContext;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;

@SuppressWarnings({"javadoc"})
public class PoissonCalculatorTest {
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
  public void canComputeLikelihoodForIntegerData() {
    final DoubleDoubleBiPredicate predicate = TestHelper.doublesAreClose(1e-10, 0);
    for (final double u : photons) {
      final PoissonDistribution pd = new PoissonDistribution(u);
      for (int x = 0; x < 100; x++) {
        double e = pd.probability(x);
        double o = PoissonCalculator.likelihood(u, x);
        if (e > 1e-100) {
          TestAssertions.assertTest(e, o, predicate);
        }
        e = pd.logProbability(x);
        o = PoissonCalculator.logLikelihood(u, x);
        TestAssertions.assertTest(e, o, predicate);
      }
    }
  }

  @Test
  public void canComputeFastLikelihoodForIntegerData() {
    final DoubleDoubleBiPredicate predicate = TestHelper.doublesAreClose(1e-4, 0);
    for (final double u : photons) {
      final PoissonDistribution pd = new PoissonDistribution(u);
      for (int x = 0; x < 100; x++) {
        double e = pd.probability(x);
        double o = PoissonCalculator.fastLikelihood(u, x);
        if (e > 1e-100) {
          TestAssertions.assertTest(e, o, predicate);
        }
        e = pd.logProbability(x);
        o = PoissonCalculator.fastLogLikelihood(u, x);
        TestAssertions.assertTest(e, o, predicate);
      }
    }
  }

  @Test
  public void canComputeFastLog_FastLikelihoodForIntegerData() {
    final DoubleDoubleBiPredicate predicate = TestHelper.doublesAreClose(1e-4, 0);
    final FastLog fastLog = FastLogFactory.getFastLog();
    for (final double u : photons) {
      final PoissonDistribution pd = new PoissonDistribution(u);
      for (int x = 0; x < 100; x++) {
        double e = pd.probability(x);
        double o = PoissonCalculator.fastLikelihood(u, x, fastLog);
        if (e > 1e-100) {
          TestAssertions.assertTest(e, o, predicate);
        }
        e = pd.logProbability(x);
        o = PoissonCalculator.fastLogLikelihood(u, x, fastLog);
        TestAssertions.assertTest(e, o, predicate);
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
      double v;
      v = likelihood(mu, x);
      // v = pgf.probability(x, mu);
      // logger.fine(FunctionUtils.getSupplier("x=%f, v=%f", x, v);
      return v;
    }

    abstract double likelihood(double mu, double x);
  }

  @Test
  public void cumulativeProbabilityIsOneWithRealDataForCountAbove4() {
    cumulativeProbabilityIsOneWithRealDataForCountAbove4(0);
  }

  @Test
  public void fastLikelihoodCumulativeProbabilityIsOneWithRealDataForCountAbove4() {
    cumulativeProbabilityIsOneWithRealDataForCountAbove4(1);
  }

  @Test
  public void fastLog_fastLikelihoodCumulativeProbabilityIsNotOneWithRealDataForCountAbove4() {
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

      PoissonFunction f;
      switch (function) {
        //@formatter:off
        case 2:
          f = new PoissonFunction(mu) { @Override
          double likelihood(double mu, double x) {
              return PoissonCalculator.fastLikelihood(mu, x, FastLogFactory.getFastLog());  } };
          break;
        case 1:
          f = new PoissonFunction(mu) { @Override
          double likelihood(double mu, double x) {
              return PoissonCalculator.fastLikelihood(mu, x);  } };
          break;
        case 0:
          f = new PoissonFunction(mu) { @Override
          double likelihood(double mu, double x) {
              return PoissonCalculator.likelihood(mu, x);  } };
          break;
        default:
          throw new IllegalStateException();
        //@formatter:on
      }

      cumulativeProbabilityIsOneWithRealData(min, max, mu >= 4, f);
    }
  }

  private static void cumulativeProbabilityIsOneWithRealData(double min, double max, boolean test,
      PoissonFunction f) {
    double p = 0;

    final int integrationpoints = 10;
    final double relativeAccuracy = 1e-4;
    final double absoluteAccuracy = 1e-8;
    final int minimalIterationCount = 3;
    final int maximalIterationCount = 32;

    final UnivariateIntegrator in = new IterativeLegendreGaussIntegrator(integrationpoints,
        relativeAccuracy, absoluteAccuracy, minimalIterationCount, maximalIterationCount);

    // new SimpsonIntegrator();

    p = in.integrate(20000, f, min, max);

    logger.log(TestLogUtils.getRecord(Level.INFO, "mu=%f, p=%f", f.mu, p));
    if (test) {
      Assertions.assertEquals(P_LIMIT, p, 0.02, () -> "mu=" + f.mu);
    }
  }

  private abstract static class BaseNonLinearFunction implements NonLinearFunction {
    double[] a;
    String name;

    BaseNonLinearFunction(String name) {
      this.name = name;
    }

    @Override
    public void initialise(double[] a) {
      this.a = a;
    }

    @Override
    public int[] gradientIndices() {
      return new int[1];
    }

    @Override
    public double eval(int x, double[] dyda, double[] w) {
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
    public double evalw(int x, double[] w) {
      return 0;
    }

    @Override
    public int getNumberOfGradients() {
      return 1;
    }
  }

  @SeededTest
  public void canComputeLogLikelihoodRatio(RandomSeed seed) {
    final double n2 = maxx * maxx * 0.5;
    // Functions must produce a strictly positive output so add background
    //@formatter:off
    canComputeLogLikelihoodRatio(seed, new BaseNonLinearFunction("Quadratic")
    {
      @Override
      public double eval(int x) {  return 0.1 + a[0] * (x-n2) * (x-n2); }
    });
    canComputeLogLikelihoodRatio(seed, new BaseNonLinearFunction("Gaussian")
    {
      @Override
      public double eval(int x) {  return 0.1 + 100 * FastMath.exp(-0.5 * MathUtils.pow2(x - n2) / (a[0] * a[0])); }
    });
    //@formatter:on
  }

  private static void canComputeLogLikelihoodRatio(RandomSeed seed, BaseNonLinearFunction nlf) {
    logger.log(TestLogUtils.getRecord(Level.INFO, nlf.name));

    final int n = maxx * maxx;

    final double[] a = new double[] {1};

    // Simulate Poisson process
    nlf.initialise(a);
    final CustomPoissonDistribution pd = new CustomPoissonDistribution(
        new RandomGeneratorAdapter(RngUtils.create(seed.getSeedAsLong())), 1);
    final double[] x = new double[n];
    final double[] u = new double[n];
    for (int i = 0; i < n; i++) {
      u[i] = nlf.eval(i);
      if (u[i] > 0) {
        pd.setMeanUnsafe(u[i]);
        x[i] = pd.sample();
      }
    }

    final DoubleDoubleBiPredicate predicate = TestHelper.doublesAreClose(1e-10, 0);

    double ll = PoissonCalculator.logLikelihood(u, x);
    final double mll = PoissonCalculator.maximumLogLikelihood(x);
    double llr = -2 * (ll - mll);
    double llr2 = PoissonCalculator.logLikelihoodRatio(u, x);
    logger.log(TestLogUtils.getRecord(Level.INFO, "llr=%f, llr2=%f", llr, llr2));
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
    if (logger.isLoggable(Level.INFO)) {
      logger.log(TestLogUtils.getRecord(Level.INFO, "Chi2 = %f (q=%.3f), %f (q=%.3f)  %f %b  %f",
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
      logger.log(TestLogUtils.getRecord(Level.INFO,
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
    int i = SimpleArrayUtils.findMaxIndex(data);
    final double maxa = (double) (imin + i) / 10;
    double fita = maxa;
    try {
      if (i == 0) {
        i++;
      }
      if (i == data.length - 1) {
        i--;
      }
      final int i1 = i - 1;
      final int i2 = i;
      final int i3 = i + 1;

      fita = QuadraticUtils.findMinMax((double) (imin + i1) / 10, data[i1],
          (double) (imin + i2) / 10, data[i2], (double) (imin + i3) / 10, data[i3]);
    } catch (final DataException ex) {
      // Ignore
    }

    // Allow a tolerance as the random data may alter the p-value computation.
    // Should allow it to be less than 2 increment either side of the answer.
    logger.log(TestLogUtils.getRecord(Level.INFO, "max fit = %g => %g", maxa, fita));
    Assertions.assertEquals(1, fita, 0.199, "max");
  }

  @SeededTest
  public void canComputeFastLog_LogLikelihoodRatio(RandomSeed seed) {
    final double n2 = maxx * maxx * 0.5;
    // Functions must produce a strictly positive output so add background
    //@formatter:off
    canComputeFastLog_LogLikelihoodRatio(seed, new BaseNonLinearFunction("Quadratic")
    {
      @Override
      public double eval(int x) {  return 0.1 + a[0] * (x-n2) * (x-n2); }
    });
    canComputeFastLog_LogLikelihoodRatio(seed,new BaseNonLinearFunction("Gaussian")
    {
      @Override
      public double eval(int x) {  return 0.1 + 100 * FastMath.exp(-0.5 * MathUtils.pow2(x - n2) / (a[0] * a[0])); }
    });
    //@formatter:on
  }

  private static void canComputeFastLog_LogLikelihoodRatio(RandomSeed seed,
      BaseNonLinearFunction nlf) {
    logger.log(TestLogUtils.getRecord(Level.INFO, nlf.name));

    final int n = maxx * maxx;

    final double[] a = new double[] {1};

    // Simulate Poisson process
    nlf.initialise(a);
    final CustomPoissonDistribution pd = new CustomPoissonDistribution(
        new RandomGeneratorAdapter(RngUtils.create(seed.getSeedAsLong())), 1);
    final double[] x = new double[n];
    final double[] u = new double[n];
    for (int i = 0; i < n; i++) {
      u[i] = nlf.eval(i);
      if (u[i] > 0) {
        pd.setMeanUnsafe(u[i]);
        x[i] = pd.sample();
      }
    }

    // Only test the LLR
    final double llr = PoissonCalculator.logLikelihoodRatio(u, x);
    final double llr2 = PoissonCalculator.logLikelihoodRatio(u, x, FastLogFactory.getFastLog());
    logger.log(TestLogUtils.getRecord(Level.INFO, "llr=%f, llr2=%f", llr, llr2));
    // Approximately equal
    TestAssertions.assertTest(llr, llr2, TestHelper.doublesAreClose(5e-3, 0),
        "Log-likelihood ratio");
  }

  @SeededTest
  public void cannotSubtractConstantBackgroundAndComputeLogLikelihoodRatio(RandomSeed seed) {
    final int n = maxx * maxx;
    final double n2 = n * 0.5;
    final double n3 = n * 0.33;
    final double n4 = n * 0.66;
    // Functions must produce a strictly positive output so add background
    //@formatter:off
    cannotSubtractConstantBackgroundAndComputeLogLikelihoodRatio(seed,
    new BaseNonLinearFunction("Quadratic")
    {
      @Override
      public double eval(int x) {  return 0.1 + a[0] * (x-n2) * (x-n2); }
    },
    new BaseNonLinearFunction("Quadratic")
    {
      @Override
      public double eval(int x) {  return 0.2 + 0.5 * a[0] * (x-n3) * (x-n3); }
    },
    new BaseNonLinearFunction("Quadratic")
    {
      @Override
      public double eval(int x) {  return 0.3 + 0.75 * a[0] * (x-n4) * (x-n4); }
    });
    cannotSubtractConstantBackgroundAndComputeLogLikelihoodRatio(seed,
    new BaseNonLinearFunction("Gaussian")
    {
      @Override
      public double eval(int x) {  return 0.1 + 100 * FastMath.exp(-0.5 * MathUtils.pow2(x - n2) / (a[0] * a[0])); }
    },
    new BaseNonLinearFunction("Gaussian")
    {
      @Override
      public double eval(int x) {  return 0.2 + 50 * FastMath.exp(-0.5 * MathUtils.pow2(x - n3) / (a[0] * a[0])); }
    },
    new BaseNonLinearFunction("Gaussian")
    {
      @Override
      public double eval(int x) {  return 0.3 + 75 * FastMath.exp(-0.5 * MathUtils.pow2(x - n4) / (a[0] * a[0])); }
    });
    //@formatter:on
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
    final CustomPoissonDistribution pd = new CustomPoissonDistribution(
        new RandomGeneratorAdapter(RngUtils.create(seed.getSeedAsLong())), 1);
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
        pd.setMeanUnsafe(u[i]);
        x[i] = pd.sample();
      }
    }

    // x is the target data
    // b1 is the already computed background
    // b2 is the first function to add to the model
    // b3 is the second function to add to the model

    // Compute the LLR of adding b3 to b2 given we already have b1 modelling data x
    final double[] b12 = add(b1, b2);
    final double ll1a = PoissonCalculator.logLikelihood(b12, x);
    final double ll2a = PoissonCalculator.logLikelihood(add(b12, b3), x);
    final double llra = -2 * (ll1a - ll2a);
    // logger.fine(FunctionUtils.getSupplier("x|(a+b+c) ll1=%f, ll2=%f, llra=%f", ll1a, ll2a, llra);

    // Compute the LLR of adding b3 to b2 given we already have x minus b1
    x = subtract(x, b1);
    final double ll1b = PoissonCalculator.logLikelihood(b2, x);
    final double ll2b = PoissonCalculator.logLikelihood(add(b2, b3), x);
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
  public void showRelativeErrorOfLogFactorialApproximation() {
    Assumptions.assumeTrue(logger.isLoggable(Level.INFO));
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.HIGH));

    double d = 1.0;
    for (int i = 1; i <= 100; i++) {
      d = Math.nextUp(d);
      showRelativeErrorOfLogFactorialApproximation(d);
    }
    for (int i = 1; i <= 300; i++) {
      showRelativeErrorOfLogFactorialApproximation(1 + i / 100.0);
    }

    for (int i = 4; i <= 100; i++) {
      showRelativeErrorOfLogFactorialApproximation(i);
    }
  }

  private static void showRelativeErrorOfLogFactorialApproximation(double x) {
    final double e = PoissonCalculator.logFactorial(x);
    final double[] o = new double[6];
    final double[] error = new double[o.length];
    for (int i = 0; i < o.length; i++) {
      o[i] = PoissonCalculator.logFactorialApproximation(x, i);
      error[i] = DoubleEquality.relativeError(e, o[i]);
    }
    logger.log(TestLogUtils.getRecord(Level.INFO, "%s! = %s : %s", Double.toString(x),
        MathUtils.rounded(e), Arrays.toString(error)));
  }

  @Test
  public void showRelativeErrorOfFastLogLikelihood() {
    Assumptions.assumeTrue(logger.isLoggable(Level.INFO));
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.HIGH));

    double d = 1.0;
    for (int i = 1; i <= 100; i++) {
      d = Math.nextUp(d);
      showRelativeErrorOfFastLogLikelihood(d);
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
      logger.log(TestLogUtils.getRecord(Level.INFO, "ll(%s|%s) = %s : %s", Double.toString(x),
          Double.toString(u), MathUtils.rounded(e), Double.toString(error)));
    }
  }

  @Test
  public void showRelativeErrorOfFastLog_FastLogLikelihood() {
    Assumptions.assumeTrue(logger.isLoggable(Level.INFO));
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.HIGH));

    double d = 1.0;
    for (int i = 1; i <= 100; i++) {
      d = Math.nextUp(d);
      showRelativeErrorOfFastLog_FastLogLikelihood(d);
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
      logger.log(TestLogUtils.getRecord(Level.INFO, "ll(%s|%s) = %s : %s", Double.toString(x),
          Double.toString(u), MathUtils.rounded(e), Double.toString(error)));
    }
  }

  @Test
  public void showRelativeErrorOfFastLog_LogLikelihoodRatio() {
    Assumptions.assumeTrue(logger.isLoggable(Level.INFO));
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.HIGH));

    double d = 1.0;
    for (int i = 1; i <= 100; i++) {
      d = Math.nextUp(d);
      showRelativeErrorOfFastLog_LogLikelihoodRatio(d);
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
      logger.log(TestLogUtils.getRecord(Level.INFO, "llr(%s|%s) = %s : %s", Double.toString(x),
          Double.toString(u), MathUtils.rounded(e), Double.toString(error)));
    }
  }

  @SeededTest
  public void instanceAndFastMethodIsApproximatelyEqualToStaticMethod(RandomSeed seed) {
    final DoubleEquality eq = new DoubleEquality(3e-4, 0);
    final UniformRandomProvider rg = RngUtils.create(seed.getSeedAsLong());
    // Test for different x. The calculator approximation begins
    final int n = 100;
    final double[] u = new double[n];
    final double[] x = new double[n];
    double e;
    double o;
    for (final double testx : new double[] {Math.nextDown(PoissonCalculator.APPROXIMATION_X),
        PoissonCalculator.APPROXIMATION_X, Math.nextUp(PoissonCalculator.APPROXIMATION_X),
        PoissonCalculator.APPROXIMATION_X * 1.1, PoissonCalculator.APPROXIMATION_X * 2,
        PoissonCalculator.APPROXIMATION_X * 10}) {
      final String X = Double.toString(testx);
      Arrays.fill(x, testx);
      final PoissonCalculator pc = new PoissonCalculator(x);
      e = PoissonCalculator.maximumLogLikelihood(x);
      o = pc.getMaximumLogLikelihood();
      logger.log(TestLogUtils.getRecord(Level.INFO, "[%s] Instance MaxLL = %g vs %g (error = %g)",
          X, e, o, DoubleEquality.relativeError(e, o)));
      Assertions.assertTrue(eq.almostEqualRelativeOrAbsolute(e, o),
          () -> "Instance Max LL not equal: x=" + X);

      o = PoissonCalculator.fastMaximumLogLikelihood(x);
      logger.log(TestLogUtils.getRecord(Level.INFO, "[%s] Fast MaxLL = %g vs %g (error = %g)", X, e,
          o, DoubleEquality.relativeError(e, o)));
      Assertions.assertTrue(eq.almostEqualRelativeOrAbsolute(e, o),
          () -> "Fast Max LL not equal: x=" + X);

      // Generate data around the value
      for (int i = 0; i < n; i++) {
        u[i] = x[i] + rg.nextDouble() - 0.5;
      }

      e = PoissonCalculator.logLikelihood(u, x);
      o = pc.logLikelihood(u);
      logger.log(TestLogUtils.getRecord(Level.INFO, "[%s] Instance LL = %g vs %g (error = %g)", X,
          e, o, DoubleEquality.relativeError(e, o)));
      Assertions.assertTrue(eq.almostEqualRelativeOrAbsolute(e, o),
          () -> "Instance LL not equal: x=" + X);

      o = PoissonCalculator.fastLogLikelihood(u, x);
      logger.log(TestLogUtils.getRecord(Level.INFO, "[%s] Fast LL = %g vs %g (error = %g)", X, e, o,
          DoubleEquality.relativeError(e, o)));
      Assertions.assertTrue(eq.almostEqualRelativeOrAbsolute(e, o),
          () -> "Fast LL not equal: x=" + X);

      e = PoissonCalculator.logLikelihoodRatio(u, x);
      o = pc.getLogLikelihoodRatio(o);

      logger.log(TestLogUtils.getRecord(Level.INFO, "[%s] Instance LLR = %g vs %g (error = %g)", X,
          e, o, DoubleEquality.relativeError(e, o)));
      Assertions.assertTrue(eq.almostEqualRelativeOrAbsolute(e, o),
          () -> "Instance LLR not equal: x=" + X);
    }
  }

  private abstract static class PCTimingTask extends BaseTimingTask {
    double[] x;
    double[] u;
    int ll;
    int llr;

    public PCTimingTask(String name, double[] x, double[] u, int ll, int llr) {
      super(String.format("%s ll=%d llr=%d", name, ll, llr));
      this.x = x;
      this.u = u;
      this.ll = ll;
      this.llr = llr;
    }

    @Override
    public int getSize() {
      return 1;
    }

    @Override
    public Object getData(int i) {
      return null;
    }
  }

  private static class StaticPCTimingTask extends PCTimingTask {
    public StaticPCTimingTask(double[] x, double[] u, int ll, int llr) {
      super("static", x, u, ll, llr);
    }

    @Override
    public Object run(Object data) {
      double value = 0;
      for (int i = 0; i < llr; i++) {
        value += PoissonCalculator.logLikelihoodRatio(u, x);
      }
      for (int i = 0; i < ll; i++) {
        value += PoissonCalculator.logLikelihood(u, x);
      }
      return value;
    }
  }

  private static class FastPCTimingTask extends PCTimingTask {
    public FastPCTimingTask(double[] x, double[] u, int ll, int llr) {
      super("fast", x, u, ll, llr);
    }

    @Override
    public Object run(Object data) {
      double value = 0;
      for (int i = 0; i < llr; i++) {
        value += PoissonCalculator.logLikelihoodRatio(u, x);
      }
      for (int i = 0; i < ll; i++) {
        value += PoissonCalculator.fastLogLikelihood(u, x);
      }
      return value;
    }
  }

  private static class FastLogPCTimingTask extends PCTimingTask {
    FastLog fastLog = FastLogFactory.getFastLog();

    public FastLogPCTimingTask(double[] x, double[] u, int ll, int llr) {
      super("fastLog", x, u, ll, llr);
    }

    @Override
    public Object run(Object data) {
      double value = 0;
      for (int i = 0; i < llr; i++) {
        value += PoissonCalculator.logLikelihoodRatio(u, x, fastLog);
      }
      for (int i = 0; i < ll; i++) {
        value += PoissonCalculator.fastLogLikelihood(u, x, fastLog);
      }
      return value;
    }
  }

  private static class InstancePCTimingTask extends PCTimingTask {
    int max;

    public InstancePCTimingTask(double[] x, double[] u, int ll, int llr) {
      super("instance", x, u, ll, llr);
      max = Math.max(llr, ll);
    }

    @Override
    public Object run(Object data) {
      final PoissonCalculator pc = new PoissonCalculator(x);
      double value = 0;
      // Use the fastest execution possible
      for (int i = 0; i < max; i++) {
        value += pc.pseudoLogLikelihood(u);
      }
      if (llr > 0) {
        value += pc.getMaximumLogLikelihood();
      }
      return value;
    }
  }

  @SpeedTag
  @Test
  public void instanceMethodIsFaster() {
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
        TestLogUtils.getRecord(Level.INFO, "Speed test x-range: %f - %f", limits[0], limits[1]));

    final TimingService ts = new TimingService(5);
    final int[] loops = new int[] {0, 1, 10};
    for (final int ll : loops) {
      for (final int llr : loops) {
        if (ll + llr == 0) {
          continue;
        }
        ts.execute(new StaticPCTimingTask(x, u, ll, llr));
        ts.execute(new FastPCTimingTask(x, u, ll, llr));
        ts.execute(new FastLogPCTimingTask(x, u, ll, llr));
        ts.execute(new InstancePCTimingTask(x, u, ll, llr));
      }
    }

    final int size = ts.getSize();
    ts.repeat(size);
    if (logger.isLoggable(Level.INFO)) {
      logger.info(ts.getReport(size));
    }

    final int index = ts.getSize() - 1;
    Assertions.assertTrue(ts.get(index).getMean() < ts.get(index - 1).getMean());
  }

  private static double[] add(double[] a, double[] b) {
    final double[] c = new double[a.length];
    for (int i = 0; i < a.length; i++) {
      c[i] = a[i] + b[i];
    }
    return c;
  }

  private static double[] subtract(double[] a, double[] b) {
    final double[] c = new double[a.length];
    for (int i = 0; i < a.length; i++) {
      c[i] = a[i] - b[i];
    }
    return c;
  }
}
