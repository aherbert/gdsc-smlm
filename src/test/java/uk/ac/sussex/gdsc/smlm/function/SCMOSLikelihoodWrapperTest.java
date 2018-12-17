package uk.ac.sussex.gdsc.smlm.function;

import uk.ac.sussex.gdsc.core.data.DataException;
import uk.ac.sussex.gdsc.core.math.QuadraticUtils;
import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.RandomGeneratorAdapter;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.rng.GaussianSamplerUtils;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import uk.ac.sussex.gdsc.smlm.math3.distribution.CustomPoissonDistribution;
import uk.ac.sussex.gdsc.test.api.TestAssertions;
import uk.ac.sussex.gdsc.test.api.TestHelper;
import uk.ac.sussex.gdsc.test.api.function.DoubleDoubleBiPredicate;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngUtils;
import uk.ac.sussex.gdsc.test.utils.TestLogUtils;
import uk.ac.sussex.gdsc.test.utils.functions.FunctionUtils;
import uk.ac.sussex.gdsc.test.utils.functions.IntArrayFormatSupplier;

import gnu.trove.list.array.TDoubleArrayList;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.apache.commons.math3.analysis.integration.UnivariateIntegrator;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.Precision;
import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.sampling.distribution.AhrensDieterExponentialSampler;
import org.apache.commons.rng.sampling.distribution.GaussianSampler;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import java.math.BigDecimal;
import java.math.MathContext;
import java.util.Arrays;
import java.util.concurrent.ConcurrentHashMap;
import java.util.function.Function;
import java.util.logging.Level;
import java.util.logging.Logger;

@SuppressWarnings({"javadoc"})
public class SCMOSLikelihoodWrapperTest implements Function<RandomSeed, Object> {
  private static Logger logger;
  private static ConcurrentHashMap<RandomSeed, Object> ConcurrentHashMap;

  @BeforeAll
  public static void beforeAll() {
    logger = Logger.getLogger(SCMOSLikelihoodWrapperTest.class.getName());
    ConcurrentHashMap = new ConcurrentHashMap<>();
  }

  @AfterAll
  public static void afterAll() {
    ConcurrentHashMap.clear();
    ConcurrentHashMap = null;
    logger = null;
  }

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
  private final double h_ = 0.01; // (double) (Math.pow(1e-3f, 1.0 / 3));

  private final int[] testx = new int[] {4, 5, 6};
  private final int[] testy = new int[] {4, 5, 6};
  // Do not test zero background since this is an edge case for the likelihood function
  private final double[] testbackground_ = new double[] {0.1, 1, 10};
  private final double[] testsignal1_ = new double[] {15, 55, 105};
  private final double[] testangle1_ = new double[] {Math.PI / 5, Math.PI / 3};
  private final double[] testcx1_ = new double[] {4.9, 5.3};
  private final double[] testcy1_ = new double[] {4.8, 5.2};
  private final double[] testcz1_ = new double[] {-1.5, 1.0};
  private final double[][] testw1_ =
      new double[][] {{1.1, 1.4}, {1.1, 1.7}, {1.5, 1.2}, {1.3, 1.7},};

  private double[] testbackground, testsignal1, testangle1, testcx1, testcy1, testcz1;
  private double[][] testw1;

  private static int maxx = 10;

  // Simulate per pixel noise
  private static float VAR = 57.9f;
  private static float G = 2.2f;
  private static float G_SD = 0.2f;
  private static float O = 100f;

  private class SCMOSLikelihoodWrapperTestData {
    float[] var;
    float[] g;
    float[] o;
    float[] sd;
  }

  @Override
  public Object apply(RandomSeed source) {
    final int n = maxx * maxx;
    final SCMOSLikelihoodWrapperTestData data = new SCMOSLikelihoodWrapperTestData();
    data.var = new float[n];
    data.g = new float[n];
    data.o = new float[n];
    data.sd = new float[n];
    final UniformRandomProvider rg = RngUtils.create(source.getSeed());
    final CustomPoissonDistribution pd =
        new CustomPoissonDistribution(new RandomGeneratorAdapter(rg), O);
    final GaussianSampler gs = GaussianSamplerUtils.createGaussianSampler(rg, G, G_SD);
    final AhrensDieterExponentialSampler ed = new AhrensDieterExponentialSampler(rg, VAR);
    for (int i = 0; i < n; i++) {
      data.o[i] = pd.sample();
      data.var[i] = (float) ed.sample();
      data.sd[i] = (float) Math.sqrt(data.var[i]);
      data.g[i] = (float) gs.sample();
    }
    return data;
  }

  @SeededTest
  public void fitFixedComputesGradientPerDatum(RandomSeed seed) {
    functionComputesGradientPerDatum(seed, GaussianFunctionFactory.FIT_FIXED);
  }

  @SeededTest
  public void fitCircleComputesGradientPerDatum(RandomSeed seed) {
    functionComputesGradientPerDatum(seed, GaussianFunctionFactory.FIT_CIRCLE);
  }

  @SeededTest
  public void fitFreeCircleComputesGradientPerDatum(RandomSeed seed) {
    functionComputesGradientPerDatum(seed, GaussianFunctionFactory.FIT_FREE_CIRCLE);
  }

  @SeededTest
  public void fitEllipticalComputesGradientPerDatum(RandomSeed seed) {
    functionComputesGradientPerDatum(seed, GaussianFunctionFactory.FIT_ELLIPTICAL);
  }

  @SeededTest
  public void fitNBFixedComputesGradientPerDatum(RandomSeed seed) {
    functionComputesGradientPerDatum(seed, GaussianFunctionFactory.FIT_SIMPLE_NB_FIXED);
  }

  @SeededTest
  public void fitNBCircleComputesGradientPerDatum(RandomSeed seed) {
    functionComputesGradientPerDatum(seed, GaussianFunctionFactory.FIT_SIMPLE_NB_CIRCLE);
  }

  @SeededTest
  public void fitNBFreeCircleComputesGradientPerDatum(RandomSeed seed) {
    functionComputesGradientPerDatum(seed, GaussianFunctionFactory.FIT_SIMPLE_NB_FREE_CIRCLE);
  }

  @SeededTest
  public void fitNBEllipticalComputesGradientPerDatum(RandomSeed seed) {
    functionComputesGradientPerDatum(seed, GaussianFunctionFactory.FIT_SIMPLE_NB_ELLIPTICAL);
  }

  private void functionComputesGradientPerDatum(RandomSeed seed, int flags) {
    final Gaussian2DFunction f1 = GaussianFunctionFactory.create2D(1, maxx, maxx, flags, null);
    // Setup
    // Setup
    testbackground = testbackground_;
    testsignal1 = testsignal1_;
    testcx1 = testcx1_;
    testcy1 = testcy1_;
    testcz1 = testcz1_;
    testw1 = testw1_;
    testangle1 = testangle1_;
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
    double[] a;

    SCMOSLikelihoodWrapper ff1;

    final int n = maxx * maxx;
    int count = 0, total = 0;

    final SCMOSLikelihoodWrapperTestData testData =
        (SCMOSLikelihoodWrapperTestData) ConcurrentHashMap.computeIfAbsent(seed, this);
    final float[] var = testData.var;
    final float[] g = testData.g;
    final float[] o = testData.o;
    final float[] sd = testData.sd;
    final UniformRandomProvider r = RngUtils.create(seed.getSeedAsLong());
    final CustomPoissonDistribution pd =
        new CustomPoissonDistribution(new RandomGeneratorAdapter(r), 1);
    final GaussianSampler gs = GaussianSamplerUtils.createGaussianSampler(r, 0, 1);

    for (final double background : testbackground) {
      for (final double signal1 : testsignal1) {
        for (final double cx1 : testcx1) {
          for (final double cy1 : testcy1) {
            for (final double cz1 : testcz1) {
              for (final double[] w1 : testw1) {
                for (final double angle1 : testangle1) {
                  a = createParameters(background, signal1, cx1, cy1, cz1, w1[0], w1[1], angle1);

                  // Create y as a function we would want to move towards
                  final double[] a2 = a.clone();
                  a2[targetParameter] *= 1.1;
                  f1.initialise(a2);
                  final double[] data = new double[n];
                  for (int i = 0; i < n; i++) {
                    // Simulate sCMOS camera
                    final double u = f1.eval(i);
                    pd.setMeanUnsafe(u);
                    data[i] = pd.sample() * g[i] + o[i] + gs.sample() * sd[i];
                  }

                  ff1 = new SCMOSLikelihoodWrapper(f1, a, data, n, var, g, o);

                  // Numerically solve gradient.
                  // Calculate the step size h to be an exact numerical representation
                  final double xx = a[targetParameter];

                  // Get h to minimise roundoff error
                  final double h = Precision.representableDelta(xx, h_);

                  for (final int x : testx) {
                    for (final int y : testy) {
                      final int i = y * maxx + x;
                      a[targetParameter] = xx;
                      ff1.likelihood(getVariables(indices, a), dyda, i);

                      // Evaluate at (x+h) and (x-h)
                      a[targetParameter] = xx + h;
                      final double value2 = ff1.likelihood(getVariables(indices, a), i);

                      a[targetParameter] = xx - h;
                      final double value3 = ff1.likelihood(getVariables(indices, a), i);

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
  public void fitFixedComputesGradient(RandomSeed seed) {
    functionComputesGradient(seed, GaussianFunctionFactory.FIT_FIXED);
  }

  @SeededTest
  public void fitCircleComputesGradient(RandomSeed seed) {
    functionComputesGradient(seed, GaussianFunctionFactory.FIT_CIRCLE);
  }

  @SeededTest
  public void fitFreeCircleComputesGradient(RandomSeed seed) {
    functionComputesGradient(seed, GaussianFunctionFactory.FIT_FREE_CIRCLE);
  }

  @SeededTest
  public void fitEllipticalComputesGradient(RandomSeed seed) {
    // The elliptical function gradient evaluation is worse
    final DoubleEquality tmp = eq;
    eq = eqPerDatum;
    functionComputesGradient(seed, GaussianFunctionFactory.FIT_ELLIPTICAL);
    eq = tmp;
  }

  @SeededTest
  public void fitNBFixedComputesGradient(RandomSeed seed) {
    functionComputesGradient(seed, GaussianFunctionFactory.FIT_SIMPLE_NB_FIXED);
  }

  @SeededTest
  public void fitNBCircleComputesGradient(RandomSeed seed) {
    functionComputesGradient(seed, GaussianFunctionFactory.FIT_SIMPLE_NB_CIRCLE);
  }

  @SeededTest
  public void fitNBFreeCircleComputesGradient(RandomSeed seed) {
    functionComputesGradient(seed, GaussianFunctionFactory.FIT_SIMPLE_NB_FREE_CIRCLE);
  }

  @SeededTest
  public void fitNBEllipticalComputesGradient(RandomSeed seed) {
    // The elliptical function gradient evaluation is worse
    final DoubleEquality tmp = eq;
    eq = eqPerDatum;
    functionComputesGradient(seed, GaussianFunctionFactory.FIT_SIMPLE_NB_ELLIPTICAL);
    eq = tmp;
  }

  private void functionComputesGradient(RandomSeed seed, int flags) {
    final Gaussian2DFunction f1 = GaussianFunctionFactory.create2D(1, maxx, maxx, flags, null);
    // Setup
    testbackground = testbackground_;
    testsignal1 = testsignal1_;
    testcx1 = testcx1_;
    testcy1 = testcy1_;
    testcz1 = testcz1_;
    testw1 = testw1_;
    testangle1 = testangle1_;
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
    double[] a;

    SCMOSLikelihoodWrapper ff1;

    final int n = maxx * maxx;
    int count = 0, total = 0;

    final SCMOSLikelihoodWrapperTestData testData =
        (SCMOSLikelihoodWrapperTestData) ConcurrentHashMap.computeIfAbsent(seed, this);
    final float[] var = testData.var;
    final float[] g = testData.g;
    final float[] o = testData.o;
    final float[] sd = testData.sd;
    final UniformRandomProvider r = RngUtils.create(seed.getSeedAsLong());
    final CustomPoissonDistribution pd =
        new CustomPoissonDistribution(new RandomGeneratorAdapter(r), 1);
    final GaussianSampler gs = GaussianSamplerUtils.createGaussianSampler(r, 0, 1);

    for (final double background : testbackground) {
      for (final double signal1 : testsignal1) {
        for (final double cx1 : testcx1) {
          for (final double cy1 : testcy1) {
            for (final double cz1 : testcz1) {
              for (final double[] w1 : testw1) {
                for (final double angle1 : testangle1) {
                  a = createParameters(background, signal1, cx1, cy1, cz1, w1[0], w1[1], angle1);

                  // Create y as a function we would want to move towards
                  final double[] a2 = a.clone();
                  a2[targetParameter] *= 1.3;
                  f1.initialise(a2);
                  final double[] data = new double[n];
                  for (int i = 0; i < n; i++) {
                    // Simulate sCMOS camera
                    final double u = f1.eval(i);
                    pd.setMeanUnsafe(u);
                    data[i] = pd.sample() * g[i] + o[i] + gs.sample() * sd[i];
                  }

                  ff1 = new SCMOSLikelihoodWrapper(f1, a, data, n, var, g, o);

                  // Numerically solve gradient.
                  // Calculate the step size h to be an exact numerical representation
                  final double xx = a[targetParameter];

                  // Get h to minimise roundoff error
                  final double h = Precision.representableDelta(xx, h_);

                  ff1.likelihood(getVariables(indices, a), dyda);

                  // Evaluate at (x+h) and (x-h)
                  a[targetParameter] = xx + h;
                  final double value2 = ff1.likelihood(getVariables(indices, a));

                  a[targetParameter] = xx - h;
                  final double value3 = ff1.likelihood(getVariables(indices, a));

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

  private static int findGradientIndex(Gaussian2DFunction f, int targetParameter) {
    final int i = f.findGradientIndex(targetParameter);
    Assertions.assertTrue(i >= 0, "Cannot find gradient index");
    return i;
  }

  double[] createParameters(double... args) {
    return args;
  }

  double P_LIMIT = 0.999999;

  @Test
  public void cumulativeProbabilityIsOneWithRealDataForCountAbove8() {
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

  private void cumulativeProbabilityIsOneWithRealData(final double mu, double min, double max,
      boolean test) {
    double p = 0;

    // Test using a standard Poisson-Gaussian convolution
    // min = -max;
    // final PoissonGaussianFunction pgf = PoissonGaussianFunction.createWithVariance(1, 1, VAR);

    final UnivariateIntegrator in = new SimpsonIntegrator();

    p = in.integrate(20000, new UnivariateFunction() {
      @Override
      public double value(double x) {
        double v;
        v = SCMOSLikelihoodWrapper.likelihood(mu, VAR, G, O, x);
        // v = pgf.probability(x, mu);
        // TestLog.fine(logger,"x=%f, v=%f", x, v);
        return v;
      }
    }, min, max);

    // TestLog.fine(logger,"mu=%f, p=%f", mu, p);
    if (test) {
      Assertions.assertEquals(P_LIMIT, p, 0.02, () -> "mu=" + mu);
    }
  }

  @Test
  public void instanceLikelihoodMatches() {
    for (final double mu : photons) {
      instanceLikelihoodMatches(mu, mu > 8);
    }
  }

  private void instanceLikelihoodMatches(final double mu, boolean test) {
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
      public double eval(int x, double[] dyda, double[] w) {
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
      public double evalw(int x, double[] w) {
        return 0;
      }

      @Override
      public int getNumberOfGradients() {
        return 0;
      }
    };
    SCMOSLikelihoodWrapper f = new SCMOSLikelihoodWrapper(nlf, a, k, n, var, g, o);

    final IntArrayFormatSupplier msg1 = new IntArrayFormatSupplier("computeLikelihood @ %d", 1);
    final IntArrayFormatSupplier msg2 =
        new IntArrayFormatSupplier("computeLikelihood+gradient @ %d", 1);
    double total = 0, p = 0;
    double maxp = 0;
    int maxi = 0;
    final DoubleDoubleBiPredicate predicate = TestHelper.doublesAreClose(1e-10, 0);
    for (int i = 0; i < n; i++) {
      final double nll = f.computeLikelihood(i);
      final double nll2 = f.computeLikelihood(gradient, i);
      final double nll3 =
          SCMOSLikelihoodWrapper.negativeLogLikelihood(mu, var[i], g[i], o[i], k[i]);
      total += nll;
      TestAssertions.assertTest(nll3, nll, predicate, msg1.set(0, i));
      TestAssertions.assertTest(nll3, nll2, predicate, msg2.set(0, i));
      final double pp = FastMath.exp(-nll);
      if (maxp < pp) {
        maxp = pp;
        maxi = i;
        // TestLog.fine(logger,"mu=%f, e=%f, k=%f, pp=%f", mu, mu * G + O, k[i], pp);
      }
      p += pp * step;
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
      Assertions.assertEquals(P_LIMIT, p, 0.02, () -> "mu=" + mu);
    }

    // Check the function can compute the same total
    double sum, sum2;
    sum = f.computeLikelihood();
    sum2 = f.computeLikelihood(gradient);
    TestAssertions.assertTest(total, sum, predicate, "computeLikelihood");
    TestAssertions.assertTest(total, sum2, predicate, "computeLikelihood with gradient");

    // Check the function can compute the same total after duplication
    f = f.build(nlf, a);
    sum = f.computeLikelihood();
    sum2 = f.computeLikelihood(gradient);
    TestAssertions.assertTest(total, sum, predicate, "computeLikelihood");
    TestAssertions.assertTest(total, sum2, predicate, "computeLikelihood with gradient");
  }

  private static float[] newArray(int n, float val) {
    final float[] a = new float[n];
    Arrays.fill(a, val);
    return a;
  }

  private abstract class BaseNonLinearFunction implements NonLinearFunction {
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
  public void canComputePValue(RandomSeed seed) {
    final double n2 = maxx * maxx * 0.5;
    //@formatter:off
    canComputePValue(seed,new BaseNonLinearFunction("Linear")
    {
      @Override
      public double eval(int x) {  return a[0] * (x-n2); }
    });
    canComputePValue(seed,new BaseNonLinearFunction("Quadratic")
    {
      @Override
      public double eval(int x) {  return a[0] * (x-n2) * (x-n2); }
    });
    canComputePValue(seed,new BaseNonLinearFunction("Linear+C")
    {
      @Override
      public double eval(int x) {  return 10 * a[0] + (x-n2); }
    });
    canComputePValue(seed,new BaseNonLinearFunction("Gaussian")
    {
      @Override
      public double eval(int x) {  return 100 * FastMath.exp(-0.5 * Math.pow(x - n2, 2) / (a[0] * a[0])); }
    });
    //@formatter:on
  }

  private void canComputePValue(RandomSeed seed, BaseNonLinearFunction nlf) {
    logger.log(TestLogUtils.getRecord(Level.INFO, nlf.name));

    final int n = maxx * maxx;

    final double[] a = new double[] {1};

    // Simulate sCMOS camera
    nlf.initialise(a);

    final SCMOSLikelihoodWrapperTestData testData =
        (SCMOSLikelihoodWrapperTestData) ConcurrentHashMap.computeIfAbsent(seed, this);
    final float[] var = testData.var;
    final float[] g = testData.g;
    final float[] o = testData.o;
    final float[] sd = testData.sd;
    final UniformRandomProvider r = RngUtils.create(seed.getSeedAsLong());
    final CustomPoissonDistribution pd =
        new CustomPoissonDistribution(new RandomGeneratorAdapter(r), 1);
    final GaussianSampler gs = GaussianSamplerUtils.createGaussianSampler(r, 0, 1);

    final double[] k = SimpleArrayUtils.newArray(n, 0, 1.0);
    for (int i = 0; i < n; i++) {
      double u = nlf.eval(i);
      if (u > 0) {
        pd.setMeanUnsafe(u);
        u = pd.sample();
      }
      k[i] = u * g[i] + o[i] + gs.sample() * sd[i];
    }

    final SCMOSLikelihoodWrapper f = new SCMOSLikelihoodWrapper(nlf, a, k, n, var, g, o);

    final double oll = f.computeObservedLikelihood();
    double oll2 = 0;
    final double[] op = new double[n];
    for (int j = 0; j < n; j++) {
      op[j] = SCMOSLikelihoodWrapper.likelihood((k[j] - o[j]) / g[j], var[j], g[j], o[j], k[j]);
      oll2 -= Math.log(op[j]);
    }
    logger.log(TestLogUtils.getRecord(Level.INFO, "oll=%f, oll2=%f", oll, oll2));
    final DoubleDoubleBiPredicate predicate = TestHelper.doublesAreClose(1e-10, 0);
    TestAssertions.assertTest(oll2, oll, predicate, "Observed Log-likelihood");

    final TDoubleArrayList list = new TDoubleArrayList();
    final int imin = 5, imax = 15;
    for (int i = imin; i <= imax; i++) {
      a[0] = (double) i / 10;
      final double ll = f.likelihood(a);
      list.add(ll);
      final double llr = f.computeLogLikelihoodRatio(ll);
      BigDecimal product = new BigDecimal(1);
      double ll2 = 0;
      for (int j = 0; j < n; j++) {
        final double p1 = SCMOSLikelihoodWrapper.likelihood(nlf.eval(j), var[j], g[j], o[j], k[j]);
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
    int i = SimpleArrayUtils.findMinIndex(data);
    final double mina = (double) (imin + i) / 10;
    double fita = mina;
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
    logger.log(TestLogUtils.getRecord(Level.INFO, "min fit = %g => %g", mina, fita));
    Assertions.assertEquals(1, fita, 0.199, "min");
  }
}
