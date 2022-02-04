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
import java.util.concurrent.ConcurrentHashMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.apache.commons.math3.analysis.integration.UnivariateIntegrator;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.util.Precision;
import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.sampling.distribution.ContinuousSampler;
import org.apache.commons.rng.sampling.distribution.DiscreteSampler;
import org.apache.commons.rng.sampling.distribution.SharedStateContinuousSampler;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import uk.ac.sussex.gdsc.core.data.DataException;
import uk.ac.sussex.gdsc.core.math.QuadraticUtils;
import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.rng.SamplerUtils;
import uk.ac.sussex.gdsc.smlm.GdscSmlmTestUtils;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import uk.ac.sussex.gdsc.smlm.utils.StdMath;
import uk.ac.sussex.gdsc.test.api.TestAssertions;
import uk.ac.sussex.gdsc.test.api.TestHelper;
import uk.ac.sussex.gdsc.test.api.function.DoubleDoubleBiPredicate;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngUtils;
import uk.ac.sussex.gdsc.test.utils.TestLogUtils;
import uk.ac.sussex.gdsc.test.utils.functions.FunctionUtils;
import uk.ac.sussex.gdsc.test.utils.functions.IntArrayFormatSupplier;

@SuppressWarnings({"javadoc"})
class ScmosLikelihoodWrapperTest {
  private static Logger logger;
  private static ConcurrentHashMap<RandomSeed, Object> dataCache;

  @BeforeAll
  public static void beforeAll() {
    logger = Logger.getLogger(ScmosLikelihoodWrapperTest.class.getName());
    dataCache = new ConcurrentHashMap<>();
  }

  /**
   * Clear the data cache after all tests.
   */
  @AfterAll
  public static void afterAll() {
    dataCache.clear();
    dataCache = null;
    logger = null;
  }

  static final double P_LIMIT = 0.999999;

  private final double[] photons = {1, 1.5, 2, 2.5, 3, 4, 5, 7.5, 10, 100, 1000};

  DoubleEquality eqPerDatum = new DoubleEquality(5e12, 0.01);
  DoubleEquality eq = new DoubleEquality(5e-3, 0.001);

  static String[] NAME;

  static {
    NAME = new String[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
    for (int i = 0; i < NAME.length; i++) {
      NAME[i] = Gaussian2DFunction.getName(i);
    }
  }

  // Compute as per Numerical Recipes 5.7.
  // Approximate error accuracy in single precision: Ef
  // Step size for derivatives:
  // h ~ (Ef)^(1/3) * xc
  // xc is the characteristic scale over which x changes, assumed to be 1 (not x as per NR since x
  // is close to zero)
  private final double stepH = 0.01; // (double) (Math.pow(1e-3f, 1.0 / 3));

  private final int[] testx = {4, 5, 6};
  private final int[] testy = {4, 5, 6};
  // Do not test zero background since this is an edge case for the likelihood function
  private final double[] testbackgroundOptions = {0.1, 1, 10};
  private final double[] testsignal1Options = {15, 55, 105};
  private final double[] testangle1Options = {Math.PI / 5, Math.PI / 3};
  private final double[] testcx1Options = {4.9, 5.3};
  private final double[] testcy1Options = {4.8, 5.2};
  private final double[] testcz1Options = {-1.5, 1.0};
  private final double[][] testw1Options = {{1.1, 1.4}, {1.1, 1.7}, {1.5, 1.2}, {1.3, 1.7},};

  private double[] testbackground;
  private double[] testsignal1;
  private double[] testangle1;
  private double[] testcx1;
  private double[] testcy1;
  private double[] testcz1;
  private double[][] testw1;

  private static int maxx = 10;

  // Simulate per pixel noise
  private static float VAR = 57.9f;
  private static float G = 2.2f;
  private static float G_SD = 0.2f;
  private static float O = 100f;

  private static class SCcmosLikelihoodWrapperTestData {
    float[] var;
    float[] gain;
    float[] offset;
    float[] sd;
  }

  private static Object createData(RandomSeed source) {
    final int n = maxx * maxx;
    final SCcmosLikelihoodWrapperTestData data = new SCcmosLikelihoodWrapperTestData();
    data.var = new float[n];
    data.gain = new float[n];
    data.offset = new float[n];
    data.sd = new float[n];
    final UniformRandomProvider rg = RngUtils.create(source.getSeed());
    final DiscreteSampler pd = GdscSmlmTestUtils.createPoissonSampler(rg, O);
    final SharedStateContinuousSampler gs = SamplerUtils.createGaussianSampler(rg, G, G_SD);
    final ContinuousSampler ed = SamplerUtils.createExponentialSampler(rg, VAR);
    for (int i = 0; i < n; i++) {
      data.offset[i] = pd.sample();
      data.var[i] = (float) ed.sample();
      data.sd[i] = (float) Math.sqrt(data.var[i]);
      data.gain[i] = (float) gs.sample();
    }
    return data;
  }

  @SeededTest
  void fitFixedComputesGradientPerDatum(RandomSeed seed) {
    functionComputesGradientPerDatum(seed, GaussianFunctionFactory.FIT_FIXED);
  }

  @SeededTest
  void fitCircleComputesGradientPerDatum(RandomSeed seed) {
    functionComputesGradientPerDatum(seed, GaussianFunctionFactory.FIT_CIRCLE);
  }

  @SeededTest
  void fitFreeCircleComputesGradientPerDatum(RandomSeed seed) {
    functionComputesGradientPerDatum(seed, GaussianFunctionFactory.FIT_FREE_CIRCLE);
  }

  @SeededTest
  void fitEllipticalComputesGradientPerDatum(RandomSeed seed) {
    functionComputesGradientPerDatum(seed, GaussianFunctionFactory.FIT_ELLIPTICAL);
  }

  @SeededTest
  void fitNbFixedComputesGradientPerDatum(RandomSeed seed) {
    functionComputesGradientPerDatum(seed, GaussianFunctionFactory.FIT_SIMPLE_NB_FIXED);
  }

  @SeededTest
  void fitNbCircleComputesGradientPerDatum(RandomSeed seed) {
    functionComputesGradientPerDatum(seed, GaussianFunctionFactory.FIT_SIMPLE_NB_CIRCLE);
  }

  @SeededTest
  void fitNbFreeCircleComputesGradientPerDatum(RandomSeed seed) {
    functionComputesGradientPerDatum(seed, GaussianFunctionFactory.FIT_SIMPLE_NB_FREE_CIRCLE);
  }

  @SeededTest
  void fitNbEllipticalComputesGradientPerDatum(RandomSeed seed) {
    functionComputesGradientPerDatum(seed, GaussianFunctionFactory.FIT_SIMPLE_NB_ELLIPTICAL);
  }

  private void functionComputesGradientPerDatum(RandomSeed seed, int flags) {
    final Gaussian2DFunction f1 = GaussianFunctionFactory.create2D(1, maxx, maxx, flags, null);
    // Setup
    testbackground = testbackgroundOptions;
    testsignal1 = testsignal1Options;
    testcx1 = testcx1Options;
    testcy1 = testcy1Options;
    testcz1 = testcz1Options;
    testw1 = testw1Options;
    testangle1 = testangle1Options;
    if (!f1.evaluatesBackground()) {
      testbackground = new double[] {testbackground[0]};
    }
    if (!f1.evaluatesSignal()) {
      testsignal1 = new double[] {testsignal1[0]};
    }

    if (!f1.evaluatesZ()) {
      testcz1 = new double[] {0};
    }
    boolean noSecondWidth = false;
    if (!f1.evaluatesSD0()) {
      // Just use 1 width
      testw1 = new double[][] {testw1[0]};
      // If no width 0 then assume we have no width 1 as well
      noSecondWidth = true;
    } else if (!f1.evaluatesSD1()) {
      // No evaluation of second width needs only variation in width 0 so truncate
      testw1 = Arrays.copyOf(testw1, 2);
      noSecondWidth = true;
    }
    if (noSecondWidth) {
      for (int i = 0; i < testw1.length; i++) {
        testw1[i][1] = testw1[i][0];
      }
    }
    if (!f1.evaluatesAngle()) {
      testangle1 = new double[] {0};
    }

    if (f1.evaluatesBackground()) {
      functionComputesTargetGradientPerDatum(seed, f1, Gaussian2DFunction.BACKGROUND);
    }
    if (f1.evaluatesSignal()) {
      functionComputesTargetGradientPerDatum(seed, f1, Gaussian2DFunction.SIGNAL);
    }
    functionComputesTargetGradientPerDatum(seed, f1, Gaussian2DFunction.X_POSITION);
    functionComputesTargetGradientPerDatum(seed, f1, Gaussian2DFunction.Y_POSITION);
    if (f1.evaluatesZ()) {
      functionComputesTargetGradientPerDatum(seed, f1, Gaussian2DFunction.Z_POSITION);
    }
    if (f1.evaluatesSD0()) {
      functionComputesTargetGradientPerDatum(seed, f1, Gaussian2DFunction.X_SD);
    }
    if (f1.evaluatesSD1()) {
      functionComputesTargetGradientPerDatum(seed, f1, Gaussian2DFunction.Y_SD);
    }
    if (f1.evaluatesAngle()) {
      functionComputesTargetGradientPerDatum(seed, f1, Gaussian2DFunction.ANGLE);
    }
  }

  private void functionComputesTargetGradientPerDatum(RandomSeed seed, Gaussian2DFunction f1,
      int targetParameter) {
    final int[] indices = f1.gradientIndices();
    final int gradientIndex = findGradientIndex(f1, targetParameter);
    final double[] dyda = new double[indices.length];
    double[] params;

    ScmosLikelihoodWrapper ff1;

    final int n = maxx * maxx;
    int count = 0;
    int total = 0;

    final SCcmosLikelihoodWrapperTestData testData = (SCcmosLikelihoodWrapperTestData) dataCache
        .computeIfAbsent(seed, ScmosLikelihoodWrapperTest::createData);
    final float[] var = testData.var;
    final float[] g = testData.gain;
    final float[] o = testData.offset;
    final float[] sd = testData.sd;
    final UniformRandomProvider r = RngUtils.create(seed.getSeed());
    final SharedStateContinuousSampler gs = SamplerUtils.createGaussianSampler(r, 0, 1);

    for (final double background : testbackground) {
      for (final double signal1 : testsignal1) {
        for (final double cx1 : testcx1) {
          for (final double cy1 : testcy1) {
            for (final double cz1 : testcz1) {
              for (final double[] w1 : testw1) {
                for (final double angle1 : testangle1) {
                  params =
                      createParameters(background, signal1, cx1, cy1, cz1, w1[0], w1[1], angle1);

                  // Create y as a function we would want to move towards
                  final double[] a2 = params.clone();
                  a2[targetParameter] *= 1.1;
                  f1.initialise(a2);
                  final double[] data = new double[n];
                  for (int i = 0; i < n; i++) {
                    // Simulate sCMOS camera
                    final double u = f1.eval(i);
                    data[i] = GdscSmlmTestUtils.createPoissonSampler(r, u).sample() * g[i] + o[i]
                        + gs.sample() * sd[i];
                  }

                  ff1 = new ScmosLikelihoodWrapper(f1, params, data, n, var, g, o);

                  // Numerically solve gradient.
                  // Calculate the step size h to be an exact numerical representation
                  final double xx = params[targetParameter];

                  // Get h to minimise roundoff error
                  final double h = Precision.representableDelta(xx, stepH);

                  for (final int x : testx) {
                    for (final int y : testy) {
                      final int i = y * maxx + x;
                      params[targetParameter] = xx;
                      ff1.likelihood(getVariables(indices, params), dyda, i);

                      // Evaluate at (x+h) and (x-h)
                      params[targetParameter] = xx + h;
                      final double value2 = ff1.likelihood(getVariables(indices, params), i);

                      params[targetParameter] = xx - h;
                      final double value3 = ff1.likelihood(getVariables(indices, params), i);

                      final double gradient = (value2 - value3) / (2 * h);
                      boolean ok = Math.signum(gradient) == Math.signum(dyda[gradientIndex])
                          || Math.abs(gradient - dyda[gradientIndex]) < 0.1;
                      // logger.fine(FunctionUtils.getSupplier("[%s-%s]/2*%g : %g == %g", "" +
                      // value2, "" + value3, h, gradient,
                      // dyda[gradientIndex]);
                      if (!ok) {
                        Assertions.fail(
                            NAME[targetParameter] + ": " + gradient + " != " + dyda[gradientIndex]);
                      }
                      ok = eqPerDatum.almostEqualRelativeOrAbsolute(gradient, dyda[gradientIndex]);
                      if (ok) {
                        count++;
                      }
                      total++;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    final double p = (100.0 * count) / total;
    logger.log(TestLogUtils.getRecord(Level.INFO, "Per Datum %s : %s = %d / %d (%.2f)",
        f1.getClass().getSimpleName(), NAME[targetParameter], count, total, p));
    Assertions.assertTrue(p > 90,
        () -> NAME[targetParameter] + " fraction too low per datum: " + p);
  }

  @SeededTest
  void fitFixedComputesGradient(RandomSeed seed) {
    functionComputesGradient(seed, GaussianFunctionFactory.FIT_FIXED);
  }

  @SeededTest
  void fitCircleComputesGradient(RandomSeed seed) {
    functionComputesGradient(seed, GaussianFunctionFactory.FIT_CIRCLE);
  }

  @SeededTest
  void fitFreeCircleComputesGradient(RandomSeed seed) {
    functionComputesGradient(seed, GaussianFunctionFactory.FIT_FREE_CIRCLE);
  }

  @SeededTest
  void fitEllipticalComputesGradient(RandomSeed seed) {
    // The elliptical function gradient evaluation is worse
    final DoubleEquality tmp = eq;
    eq = eqPerDatum;
    functionComputesGradient(seed, GaussianFunctionFactory.FIT_ELLIPTICAL);
    eq = tmp;
  }

  @SeededTest
  void fitNbFixedComputesGradient(RandomSeed seed) {
    functionComputesGradient(seed, GaussianFunctionFactory.FIT_SIMPLE_NB_FIXED);
  }

  @SeededTest
  void fitNbCircleComputesGradient(RandomSeed seed) {
    functionComputesGradient(seed, GaussianFunctionFactory.FIT_SIMPLE_NB_CIRCLE);
  }

  @SeededTest
  void fitNbFreeCircleComputesGradient(RandomSeed seed) {
    functionComputesGradient(seed, GaussianFunctionFactory.FIT_SIMPLE_NB_FREE_CIRCLE);
  }

  @SeededTest
  void fitNbEllipticalComputesGradient(RandomSeed seed) {
    // The elliptical function gradient evaluation is worse
    final DoubleEquality tmp = eq;
    eq = eqPerDatum;
    functionComputesGradient(seed, GaussianFunctionFactory.FIT_SIMPLE_NB_ELLIPTICAL);
    eq = tmp;
  }

  private void functionComputesGradient(RandomSeed seed, int flags) {
    final Gaussian2DFunction f1 = GaussianFunctionFactory.create2D(1, maxx, maxx, flags, null);
    // Setup
    testbackground = testbackgroundOptions;
    testsignal1 = testsignal1Options;
    testcx1 = testcx1Options;
    testcy1 = testcy1Options;
    testcz1 = testcz1Options;
    testw1 = testw1Options;
    testangle1 = testangle1Options;
    if (!f1.evaluatesBackground()) {
      testbackground = new double[] {testbackground[0]};
    }
    if (!f1.evaluatesSignal()) {
      testsignal1 = new double[] {testsignal1[0]};
    }

    if (!f1.evaluatesZ()) {
      testcz1 = new double[] {0};
    }
    boolean noSecondWidth = false;
    if (!f1.evaluatesSD0()) {
      // Just use 1 width
      testw1 = new double[][] {testw1[0]};
      // If no width 0 then assume we have no width 1 as well
      noSecondWidth = true;
    } else if (!f1.evaluatesSD1()) {
      // No evaluation of second width needs only variation in width 0 so truncate
      testw1 = Arrays.copyOf(testw1, 2);
      noSecondWidth = true;
    }
    if (noSecondWidth) {
      for (int i = 0; i < testw1.length; i++) {
        testw1[i][1] = testw1[i][0];
      }
    }
    if (!f1.evaluatesAngle()) {
      testangle1 = new double[] {0};
    }

    final double fraction = 85;
    if (f1.evaluatesBackground()) {
      functionComputesTargetGradient(seed, f1, Gaussian2DFunction.BACKGROUND, fraction);
    }
    if (f1.evaluatesSignal()) {
      functionComputesTargetGradient(seed, f1, Gaussian2DFunction.SIGNAL, fraction);
    }
    functionComputesTargetGradient(seed, f1, Gaussian2DFunction.X_POSITION, fraction);
    functionComputesTargetGradient(seed, f1, Gaussian2DFunction.Y_POSITION, fraction);
    if (f1.evaluatesZ()) {
      functionComputesTargetGradient(seed, f1, Gaussian2DFunction.Z_POSITION, fraction);
    }
    if (f1.evaluatesSD0()) {
      functionComputesTargetGradient(seed, f1, Gaussian2DFunction.X_SD, fraction);
    }
    if (f1.evaluatesSD1()) {
      functionComputesTargetGradient(seed, f1, Gaussian2DFunction.Y_SD, fraction);
    }
    if (f1.evaluatesAngle()) {
      functionComputesTargetGradient(seed, f1, Gaussian2DFunction.ANGLE, fraction);
    }
  }

  private void functionComputesTargetGradient(RandomSeed seed, Gaussian2DFunction f1,
      int targetParameter, double threshold) {
    final int[] indices = f1.gradientIndices();
    final int gradientIndex = findGradientIndex(f1, targetParameter);
    final double[] dyda = new double[indices.length];
    double[] params;

    ScmosLikelihoodWrapper ff1;

    final int n = maxx * maxx;
    int count = 0;
    int total = 0;

    final SCcmosLikelihoodWrapperTestData testData = (SCcmosLikelihoodWrapperTestData) dataCache
        .computeIfAbsent(seed, ScmosLikelihoodWrapperTest::createData);
    final float[] var = testData.var;
    final float[] g = testData.gain;
    final float[] o = testData.offset;
    final float[] sd = testData.sd;
    final UniformRandomProvider r = RngUtils.create(seed.getSeed());
    final SharedStateContinuousSampler gs = SamplerUtils.createGaussianSampler(r, 0, 1);

    for (final double background : testbackground) {
      for (final double signal1 : testsignal1) {
        for (final double cx1 : testcx1) {
          for (final double cy1 : testcy1) {
            for (final double cz1 : testcz1) {
              for (final double[] w1 : testw1) {
                for (final double angle1 : testangle1) {
                  params =
                      createParameters(background, signal1, cx1, cy1, cz1, w1[0], w1[1], angle1);

                  // Create y as a function we would want to move towards
                  final double[] a2 = params.clone();
                  a2[targetParameter] *= 1.3;
                  f1.initialise(a2);
                  final double[] data = new double[n];
                  for (int i = 0; i < n; i++) {
                    // Simulate sCMOS camera
                    final double u = f1.eval(i);
                    data[i] = GdscSmlmTestUtils.createPoissonSampler(r, u).sample() * g[i] + o[i]
                        + gs.sample() * sd[i];
                  }

                  ff1 = new ScmosLikelihoodWrapper(f1, params, data, n, var, g, o);

                  // Numerically solve gradient.
                  // Calculate the step size h to be an exact numerical representation
                  final double xx = params[targetParameter];

                  // Get h to minimise roundoff error
                  final double h = Precision.representableDelta(xx, stepH);

                  ff1.likelihood(getVariables(indices, params), dyda);

                  // Evaluate at (x+h) and (x-h)
                  params[targetParameter] = xx + h;
                  final double value2 = ff1.likelihood(getVariables(indices, params));

                  params[targetParameter] = xx - h;
                  final double value3 = ff1.likelihood(getVariables(indices, params));

                  final double gradient = (value2 - value3) / (2 * h);
                  boolean ok = Math.signum(gradient) == Math.signum(dyda[gradientIndex])
                      || Math.abs(gradient - dyda[gradientIndex]) < 0.1;
                  // logger.fine(FunctionUtils.getSupplier("[%s-%s]/2*%g : %g == %g", "" + value2,
                  // "" + value3, h, gradient,
                  // dyda[gradientIndex]));
                  if (!ok) {
                    Assertions.fail(
                        NAME[targetParameter] + ": " + gradient + " != " + dyda[gradientIndex]);
                  }
                  ok = eq.almostEqualRelativeOrAbsolute(gradient, dyda[gradientIndex]);
                  if (ok) {
                    count++;
                  }
                  total++;

                }
              }
            }
          }
        }
      }
    }
    final double p = (100.0 * count) / total;
    logger.log(TestLogUtils.getRecord(Level.INFO, "%s : %s = %d / %d (%.2f)",
        f1.getClass().getSimpleName(), NAME[targetParameter], count, total, p));
    Assertions.assertTrue(p > threshold,
        FunctionUtils.getSupplier("%s fraction too low: %s", NAME[targetParameter], p));
  }

  private static double[] getVariables(int[] indices, double[] a) {
    final double[] variables = new double[indices.length];
    for (int i = 0; i < indices.length; i++) {
      variables[i] = a[indices[i]];
    }
    return variables;
  }

  private static int findGradientIndex(Gaussian2DFunction func, int targetParameter) {
    final int index = func.findGradientIndex(targetParameter);
    Assertions.assertTrue(index >= 0, "Cannot find gradient index");
    return index;
  }

  double[] createParameters(double... args) {
    return args;
  }

  @Test
  void cumulativeProbabilityIsOneWithRealDataForCountAbove8() {
    for (final double mu : photons) {
      // Determine upper limit for a Poisson
      double max = new PoissonDistribution(mu).inverseCumulativeProbability(P_LIMIT);

      // Determine lower limit
      final double sd = Math.sqrt(mu);
      double min = (int) Math.max(0, mu - 4 * sd);

      // Map to observed values using the gain and offset
      max = max * G + O;
      min = min * G + O;

      cumulativeProbabilityIsOneWithRealData(mu, min, max, mu > 8);
    }
  }

  private static void cumulativeProbabilityIsOneWithRealData(final double mu, double min,
      double max, boolean test) {
    // Test using a standard Poisson-Gaussian convolution
    // min = -max;
    // final PoissonGaussianFunction pgf = PoissonGaussianFunction.createWithVariance(1, 1, VAR);

    final UnivariateIntegrator in = new SimpsonIntegrator();

    final double pvalue =
        in.integrate(20000, x -> ScmosLikelihoodWrapper.likelihood(mu, VAR, G, O, x), min, max);

    // TestLog.fine(logger,"mu=%f, p=%f", mu, p);
    if (test) {
      Assertions.assertEquals(P_LIMIT, pvalue, 0.02, () -> "mu=" + mu);
    }
  }

  @Test
  void instanceLikelihoodMatches() {
    for (final double mu : photons) {
      instanceLikelihoodMatches(mu, mu > 8);
    }
  }

  private static void instanceLikelihoodMatches(final double mu, boolean test) {
    // Determine upper limit for a Poisson
    final int limit = new PoissonDistribution(mu).inverseCumulativeProbability(P_LIMIT);

    // Map to observed values using the gain and offset
    final double max = limit * G;

    final double step = 0.1;

    final int n = (int) Math.ceil(max / step);

    // Evaluate all values from (zero+offset) to large n
    final double[] k = SimpleArrayUtils.newArray(n, O, step);
    final double[] a = new double[0];
    final double[] gradient = new double[0];

    final float[] var = newArray(n, VAR);
    final float[] g = newArray(n, G);
    final float[] o = newArray(n, O);

    final NonLinearFunction nlf = new NonLinearFunction() {
      @Override
      public void initialise(double[] a) {
        // Ignore
      }

      @Override
      public int[] gradientIndices() {
        return new int[0];
      }

      @Override
      public double evalw(int x, double[] dyda, double[] weight) {
        return 0;
      }

      @Override
      public double evalw(int x, double[] weight) {
        return 0;
      }

      @Override
      public double eval(int x) {
        return mu;
      }

      @Override
      public double eval(int x, double[] dyda) {
        return mu;
      }

      @Override
      public boolean canComputeWeights() {
        return false;
      }

      @Override
      public int getNumberOfGradients() {
        return 0;
      }
    };
    ScmosLikelihoodWrapper func = new ScmosLikelihoodWrapper(nlf, a, k, n, var, g, o);

    final IntArrayFormatSupplier msg1 = new IntArrayFormatSupplier("computeLikelihood @ %d", 1);
    final IntArrayFormatSupplier msg2 =
        new IntArrayFormatSupplier("computeLikelihood+gradient @ %d", 1);
    double total = 0;
    double pvalue = 0;
    double maxp = 0;
    int maxi = 0;
    final DoubleDoubleBiPredicate predicate = TestHelper.doublesAreClose(1e-10, 0);
    for (int i = 0; i < n; i++) {
      final double nll = func.computeLikelihood(i);
      final double nll2 = func.computeLikelihood(gradient, i);
      final double nll3 =
          ScmosLikelihoodWrapper.negativeLogLikelihood(mu, var[i], g[i], o[i], k[i]);
      total += nll;
      TestAssertions.assertTest(nll3, nll, predicate, msg1.set(0, i));
      TestAssertions.assertTest(nll3, nll2, predicate, msg2.set(0, i));
      final double pp = StdMath.exp(-nll);
      if (maxp < pp) {
        maxp = pp;
        maxi = i;
        // TestLog.fine(logger,"mu=%f, e=%f, k=%f, pp=%f", mu, mu * G + O, k[i], pp);
      }
      pvalue += pp * step;
    }

    // Expected max of the distribution is the mode of the Poisson distribution.
    // This has two modes for integer input counts. We take the mean of those.
    // https://en.wikipedia.org/wiki/Poisson_distribution
    // Note that the shift of VAR/(G*G) is a constant applied to both the expected and
    // observed values and consequently cancels when predicting the max, i.e. we add
    // a constant count to the observed values and shift the distribution by the same
    // constant. We can thus compute the mode for the unshifted distribution.
    final double lambda = mu;
    final double mode1 = Math.floor(lambda);
    final double mode2 = Math.ceil(lambda) - 1;
    final double kmax = ((mode1 + mode2) * 0.5) * G + O; // Scale to observed values
    // TestLog.fine(logger,"mu=%f, p=%f, maxp=%f @ %f (expected=%f %f)", mu, p, maxp, k[maxi], kmax,
    // kmax - k[maxi]);
    TestAssertions.assertTest(kmax, k[maxi], TestHelper.doublesAreClose(1e-3, 0), "k-max");

    if (test) {
      Assertions.assertEquals(P_LIMIT, pvalue, 0.02, () -> "mu=" + mu);
    }

    // Check the function can compute the same total
    double sum;
    double sum2;
    sum = func.computeLikelihood();
    sum2 = func.computeLikelihood(gradient);
    TestAssertions.assertTest(total, sum, predicate, "computeLikelihood");
    TestAssertions.assertTest(total, sum2, predicate, "computeLikelihood with gradient");

    // Check the function can compute the same total after duplication
    func = func.build(nlf, a);
    sum = func.computeLikelihood();
    sum2 = func.computeLikelihood(gradient);
    TestAssertions.assertTest(total, sum, predicate, "computeLikelihood");
    TestAssertions.assertTest(total, sum2, predicate, "computeLikelihood with gradient");
  }

  private static float[] newArray(int n, float val) {
    final float[] a = new float[n];
    Arrays.fill(a, val);
    return a;
  }

  private abstract class BaseNonLinearFunction implements NonLinearFunction {
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
    public double evalw(int x, double[] dyda, double[] weight) {
      return 0;
    }

    @Override
    public double evalw(int x, double[] weight) {
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
  void canComputePValue(RandomSeed seed) {
    final double n2 = maxx * maxx * 0.5;
    //@formatter:off
    canComputePValue(seed,new BaseNonLinearFunction("Linear")
    {
      @Override
      public double eval(int x) {  return params[0] * (x-n2); }
    });
    canComputePValue(seed,new BaseNonLinearFunction("Quadratic")
    {
      @Override
      public double eval(int x) {  return params[0] * (x-n2) * (x-n2); }
    });
    canComputePValue(seed,new BaseNonLinearFunction("Linear+C")
    {
      @Override
      public double eval(int x) {  return 10 * params[0] + (x-n2); }
    });
    canComputePValue(seed,new BaseNonLinearFunction("Gaussian")
    {
      @Override
      public double eval(int x) {  return 100 * StdMath.exp(
          -0.5 * Math.pow(x - n2, 2) / (params[0] * params[0])); }
    });
    //@formatter:on
  }

  private static void canComputePValue(RandomSeed seed, BaseNonLinearFunction nlf) {
    logger.log(TestLogUtils.getRecord(Level.INFO, nlf.name));

    final int n = maxx * maxx;

    final double[] a = new double[] {1};

    // Simulate sCMOS camera
    nlf.initialise(a);

    final SCcmosLikelihoodWrapperTestData testData = (SCcmosLikelihoodWrapperTestData) dataCache
        .computeIfAbsent(seed, ScmosLikelihoodWrapperTest::createData);
    final float[] var = testData.var;
    final float[] g = testData.gain;
    final float[] o = testData.offset;
    final float[] sd = testData.sd;
    final UniformRandomProvider r = RngUtils.create(seed.getSeed());
    final SharedStateContinuousSampler gs = SamplerUtils.createGaussianSampler(r, 0, 1);

    final double[] k = SimpleArrayUtils.newArray(n, 0, 1.0);
    for (int i = 0; i < n; i++) {
      double mean = nlf.eval(i);
      if (mean > 0) {
        mean = GdscSmlmTestUtils.createPoissonSampler(r, mean).sample();
      }
      k[i] = mean * g[i] + o[i] + gs.sample() * sd[i];
    }

    final ScmosLikelihoodWrapper f = new ScmosLikelihoodWrapper(nlf, a, k, n, var, g, o);

    final double oll = f.computeObservedLikelihood();
    double oll2 = 0;
    final double[] op = new double[n];
    for (int j = 0; j < n; j++) {
      op[j] = ScmosLikelihoodWrapper.likelihood((k[j] - o[j]) / g[j], var[j], g[j], o[j], k[j]);
      oll2 -= Math.log(op[j]);
    }
    logger.log(TestLogUtils.getRecord(Level.INFO, "oll=%f, oll2=%f", oll, oll2));
    final DoubleDoubleBiPredicate predicate = TestHelper.doublesAreClose(1e-10, 0);
    TestAssertions.assertTest(oll2, oll, predicate, "Observed Log-likelihood");

    final TDoubleArrayList list = new TDoubleArrayList();
    final int imin = 5;
    final int imax = 15;
    for (int i = imin; i <= imax; i++) {
      a[0] = (double) i / 10;
      final double ll = f.likelihood(a);
      list.add(ll);
      final double llr = f.computeLogLikelihoodRatio(ll);
      BigDecimal product = new BigDecimal(1);
      double ll2 = 0;
      for (int j = 0; j < n; j++) {
        final double p1 = ScmosLikelihoodWrapper.likelihood(nlf.eval(j), var[j], g[j], o[j], k[j]);
        ll2 -= Math.log(p1);
        final double ratio = p1 / op[j];
        product = product.multiply(new BigDecimal(ratio));
      }
      final double llr2 = -2 * Math.log(product.doubleValue());
      final double q = f.computeQValue(ll);
      logger.log(TestLogUtils.getRecord(Level.INFO,
          "a=%f, ll=%f, ll2=%f, llr=%f, llr2=%f, product=%s, p=%f", a[0], ll, ll2, llr, llr2,
          product.round(new MathContext(4)).toString(), q));

      // Only value if the product could be computed. Low ratios cause it to becomes
      // too small to store in a double.
      if (product.doubleValue() > 0) {
        TestAssertions.assertTest(llr, llr2, predicate, "Log-likelihood");
      }
    }

    // Find min using quadratic fit
    final double[] data = list.toArray();
    int index = SimpleArrayUtils.findMinIndex(data);
    final double mina = (double) (imin + index) / 10;
    double fita = mina;
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
    logger.log(TestLogUtils.getRecord(Level.INFO, "min fit = %g => %g", mina, fita));
    Assertions.assertEquals(1, fita, 0.199, "min");
  }
}
