/*-
 * #%L
 * Genome Damage and Stability Centre SMLM Package
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2023 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.function.gaussian;

import java.util.logging.Logger;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.test.api.Predicates;
import uk.ac.sussex.gdsc.test.api.TestAssertions;
import uk.ac.sussex.gdsc.test.api.function.DoubleDoubleBiPredicate;
import uk.ac.sussex.gdsc.test.junit5.SpeedTag;
import uk.ac.sussex.gdsc.test.utils.TestComplexity;
import uk.ac.sussex.gdsc.test.utils.TestLogging;
import uk.ac.sussex.gdsc.test.utils.TestLogging.TestLevel;
import uk.ac.sussex.gdsc.test.utils.TestSettings;

/**
 * Contains tests for the Gaussian functions in single or double precision
 *
 * <p>The tests show that there is very little (if any) time penalty when using double precision for
 * the calculations. However the precision of the single-precision functions is 1e-4 when using
 * reasonable Gaussian parameters. This could effect the convergence of optimisers/fitters if using
 * single precision math.
 */
@SuppressWarnings({"javadoc"})
class PrecisionTest {
  private static Logger logger;

  @BeforeAll
  public static void beforeAll() {
    logger = Logger.getLogger(PrecisionTest.class.getName());
  }

  @AfterAll
  public static void afterAll() {
    logger = null;
  }

  static final int SINGLE = 1;
  static final int DOUBLE = 2;

  private static final int MAX_ITER = 200000;

  int maxx = 10;
  // Use realistic values for a camera with a bias of 500
  static double[] params2 = new double[] {500.23, 300.12, 0, 5.12, 5.23, 1.11, 1.11};
  static float[] params1 = toFloat(params2);

  // Stripped down Gaussian functions copied from the
  // uk.ac.sussex.gdsc.smlm.fitting.function.gaussian package
  public abstract class Gaussian {
    public static final int BACKGROUND = 0;
    public static final int AMPLITUDE = 1;
    public static final int ANGLE = 2;
    public static final int X_POSITION = 3;
    public static final int Y_POSITION = 4;
    public static final int X_SD = 5;
    public static final int Y_SD = 6;

    int maxx;

    public Gaussian(int maxx) {
      this.maxx = maxx;
    }

    public void setMaxX(int maxx) {
      this.maxx = maxx;
    }
  }

  public interface DoublePrecision {
    void setMaxX(int maxx);

    void initialise(double[] a);

    double eval(final int x, final double[] dyda);

    double eval(final int x);
  }

  public interface SinglePrecision {
    void setMaxX(int maxx);

    void initialise(float[] a);

    float eval(final int x, final float[] dyda);

    float eval(final int x);
  }

  public class DoubleCircularGaussian extends Gaussian implements DoublePrecision {
    double background;
    double amplitude;
    double x0pos;
    double x1pos;

    double aa;
    double aa2;
    double ax;

    public DoubleCircularGaussian(int maxx) {
      super(maxx);
    }

    @Override
    public void initialise(double[] a) {
      background = a[BACKGROUND];
      amplitude = a[AMPLITUDE];
      x0pos = a[X_POSITION];
      x1pos = a[Y_POSITION];

      final double sx = a[X_SD];
      final double sx2 = sx * sx;
      final double sx3 = sx2 * sx;

      aa = -0.5 / sx2;
      aa2 = -2.0 * aa;

      // For the x-width gradient
      ax = 1.0 / sx3;
    }

    @Override
    public double eval(final int x, final double[] dyda) {
      dyda[0] = 1.0;

      final int x1 = x / maxx;
      final int x0 = x % maxx;

      return background + gaussian(x0, x1, dyda);
    }

    @Override
    public double eval(final int x) {
      final int x1 = x / maxx;
      final int x0 = x % maxx;

      final double dx = x0 - x0pos;
      final double dy = x1 - x1pos;

      return background + amplitude * Math.exp(aa * (dx * dx + dy * dy));
    }

    private double gaussian(final int x0, final int x1, final double[] dyda) {
      final double h = amplitude;

      final double dx = x0 - x0pos;
      final double dy = x1 - x1pos;
      final double dx2dy2 = dx * dx + dy * dy;

      dyda[1] = Math.exp(aa * (dx2dy2));
      final double y = h * dyda[1];
      final double yaa2 = y * aa2;
      dyda[2] = yaa2 * dx;
      dyda[3] = yaa2 * dy;

      dyda[4] = y * (ax * (dx2dy2));

      return y;
    }
  }

  public class SingleCircularGaussian extends Gaussian implements SinglePrecision {
    float background;
    float amplitude;
    float x0pos;
    float x1pos;

    float aa;
    float aa2;
    float ax;

    public SingleCircularGaussian(int maxx) {
      super(maxx);
    }

    @Override
    public void initialise(float[] a) {
      background = a[BACKGROUND];
      amplitude = a[AMPLITUDE];
      x0pos = a[X_POSITION];
      x1pos = a[Y_POSITION];

      final float sx = a[X_SD];
      final float sx2 = sx * sx;
      final float sx3 = sx2 * sx;

      aa = -0.5f / sx2;
      aa2 = -2.0f * aa;

      ax = 1.0f / sx3;
    }

    @Override
    public float eval(final int x, final float[] dyda) {
      dyda[0] = 1.0f;

      final int x1 = x / maxx;
      final int x0 = x % maxx;

      return background + gaussian(x0, x1, dyda);
    }

    @Override
    public float eval(final int x) {
      final int x1 = x / maxx;
      final int x0 = x % maxx;

      final float dx = x0 - x0pos;
      final float dy = x1 - x1pos;

      return background + amplitude * (float) (Math.exp(aa * (dx * dx + dy * dy)));
    }

    private float gaussian(final int x0, final int x1, final float[] dyda) {
      final float h = amplitude;

      final float dx = x0 - x0pos;
      final float dy = x1 - x1pos;
      final float dx2dy2 = dx * dx + dy * dy;

      dyda[1] = (float) Math.exp(aa * (dx2dy2));
      final float y = h * dyda[1];
      final float yaa2 = y * aa2;
      dyda[2] = yaa2 * dx;
      dyda[3] = yaa2 * dy;

      dyda[4] = y * (ax * (dx2dy2));

      return y;
    }
  }

  public class DoubleFixedGaussian extends Gaussian implements DoublePrecision {
    double width;

    double background;
    double amplitude;
    double x0pos;
    double x1pos;

    double aa;
    double aa2;

    public DoubleFixedGaussian(int maxx) {
      super(maxx);
    }

    @Override
    public void initialise(double[] a) {
      background = a[BACKGROUND];
      amplitude = a[AMPLITUDE];
      x0pos = a[X_POSITION];
      x1pos = a[Y_POSITION];
      width = a[X_SD];

      final double sx = a[X_SD];
      final double sx2 = sx * sx;

      aa = -0.5 / sx2;
      aa2 = -2.0 * aa;
    }

    @Override
    public double eval(final int x, final double[] dyda) {
      dyda[0] = 1.0;

      final int x1 = x / maxx;
      final int x0 = x % maxx;

      return background + gaussian(x0, x1, dyda);
    }

    @Override
    public double eval(final int x) {
      final int x1 = x / maxx;
      final int x0 = x % maxx;

      final double dx = x0 - x0pos;
      final double dy = x1 - x1pos;

      return background + amplitude * Math.exp(aa * (dx * dx + dy * dy));
    }

    private double gaussian(final int x0, final int x1, final double[] dyda) {
      final double h = amplitude;

      final double dx = x0 - x0pos;
      final double dy = x1 - x1pos;

      dyda[1] = Math.exp(aa * (dx * dx + dy * dy));
      final double y = h * dyda[1];
      final double yaa2 = y * aa2;
      dyda[2] = yaa2 * dx;
      dyda[3] = yaa2 * dy;

      return y;
    }
  }

  public class SingleFixedGaussian extends Gaussian implements SinglePrecision {
    float width;

    float background;
    float amplitude;
    float x0pos;
    float x1pos;

    float aa;
    float aa2;

    public SingleFixedGaussian(int maxx) {
      super(maxx);
    }

    @Override
    public void initialise(float[] a) {
      background = a[BACKGROUND];
      amplitude = a[AMPLITUDE];
      x0pos = a[X_POSITION];
      x1pos = a[Y_POSITION];
      width = a[X_SD];

      final float sx = a[X_SD];
      final float sx2 = sx * sx;

      aa = -0.5f / sx2;
      aa2 = -2.0f * aa;
    }

    @Override
    public float eval(final int x, final float[] dyda) {
      dyda[0] = 1.0f;

      final int x1 = x / maxx;
      final int x0 = x % maxx;

      return background + gaussian(x0, x1, dyda);
    }

    @Override
    public float eval(final int x) {
      final int x1 = x / maxx;
      final int x0 = x % maxx;

      final float dx = x0 - x0pos;
      final float dy = x1 - x1pos;

      return background + amplitude * (float) (Math.exp(aa * (dx * dx + dy * dy)));
    }

    private float gaussian(final int x0, final int x1, final float[] dyda) {
      final float h = amplitude;

      final float dx = x0 - x0pos;
      final float dy = x1 - x1pos;

      dyda[1] = (float) (Math.exp(aa * (dx * dx + dy * dy)));
      final float y = h * dyda[1];
      final float yaa2 = y * aa2;
      dyda[2] = yaa2 * dx;
      dyda[3] = yaa2 * dy;

      return y;
    }
  }

  @Test
  void circularFunctionPrecisionIs3sf() {
    functionsComputeSameValue(maxx, new SingleCircularGaussian(maxx),
        new DoubleCircularGaussian(maxx), 1e-3);
  }

  @Test
  void circularFunctionPrecisionIs4sf() {
    functionsComputeSameValue(maxx, new SingleCircularGaussian(maxx),
        new DoubleCircularGaussian(maxx), 1e-4);
  }

  @Test
  void circularFunctionPrecisionIsNot5sf() {
    Assertions.assertThrows(AssertionError.class, () -> {
      functionsComputeSameValue(maxx, new SingleCircularGaussian(maxx),
          new DoubleCircularGaussian(maxx), 1e-5);
    });
  }

  @Test
  void circularFunctionsPrecisionIsNot3sfAtLargeXy() {
    int maxx = this.maxx;
    try {
      maxx *= 2;
      while (maxx * maxx < Integer.MAX_VALUE) {
        logger.log(TestLogging.getRecord(TestLevel.TEST_INFO, "maxx = %d", maxx));
        functionsComputeSameValue(maxx, new SingleCircularGaussian(maxx),
            new DoubleCircularGaussian(maxx), 1e-3);
        maxx *= 2;
      }
    } catch (final AssertionError ex) {
      logger.log(TestLogging.getRecord(TestLevel.TEST_INFO, ex.getMessage()));
      // ex.printStackTrace();
      return;
    }
    Assertions.fail("Expected different value");
  }

  @SpeedTag
  @Test
  void circularDoublePrecisionIsFasterWithGradients() {
    isFasterWithGradients(maxx, new SingleCircularGaussian(maxx), new DoubleCircularGaussian(maxx),
        false, true);
  }

  @SpeedTag
  @Test
  void circularDoublePrecisionIsFaster() {
    isFaster(maxx, new SingleCircularGaussian(maxx), new DoubleCircularGaussian(maxx), false, true);
  }

  @SpeedTag
  @Test
  void circularDoublePrecisionIsFasterWithGradientsNoSum() {
    isFasterWithGradients(maxx, new SingleCircularGaussian(maxx), new DoubleCircularGaussian(maxx),
        true, true);
  }

  @SpeedTag
  @Test
  void circularDoublePrecisionIsFasterNoSum() {
    isFaster(maxx, new SingleCircularGaussian(maxx), new DoubleCircularGaussian(maxx), true, true);
  }

  @Test
  void fixedFunctionPrecisionIs3sf() {
    functionsComputeSameValue(maxx, new SingleFixedGaussian(maxx), new DoubleFixedGaussian(maxx),
        1e-3);
  }

  @Test
  void fixedFunctionPrecisionIs4sf() {
    functionsComputeSameValue(maxx, new SingleFixedGaussian(maxx), new DoubleFixedGaussian(maxx),
        1e-4);
  }

  @Test
  void fixedFunctionPrecisionIsNot6sf() {
    Assertions.assertThrows(AssertionError.class, () -> {
      functionsComputeSameValue(maxx, new SingleFixedGaussian(maxx), new DoubleFixedGaussian(maxx),
          1e-6);
    });
  }

  @Test
  void fixedFunctionsPrecisionIsNot3sfAtLargeXy() {
    int maxx = this.maxx;
    try {
      maxx *= 2;
      while (maxx * maxx < Integer.MAX_VALUE) {
        logger.log(TestLogging.getRecord(TestLevel.TEST_INFO, "maxx = %d", maxx));
        functionsComputeSameValue(maxx, new SingleFixedGaussian(maxx),
            new DoubleFixedGaussian(maxx), 1e-3);
        maxx *= 2;
      }
    } catch (final AssertionError ex) {
      logger.log(TestLogging.getRecord(TestLevel.TEST_INFO, ex.getMessage()));
      // ex.printStackTrace();
      return;
    }
    Assertions.fail("Expected different value");
  }

  @SpeedTag
  @Test
  void fixedDoublePrecisionIsFasterWithGradients() {
    isFasterWithGradients(maxx, new SingleFixedGaussian(maxx), new DoubleFixedGaussian(maxx), false,
        true);
  }

  @SpeedTag
  @Test
  void fixedDoublePrecisionIsFaster() {
    isFaster(maxx, new SingleFixedGaussian(maxx), new DoubleFixedGaussian(maxx), false, true);
  }

  @SpeedTag
  @Test
  void fixedDoublePrecisionIsFasterWithGradientsNoSum() {
    isFasterWithGradients(maxx, new SingleFixedGaussian(maxx), new DoubleFixedGaussian(maxx), true,
        true);
  }

  @SpeedTag
  @Test
  void fixedDoublePrecisionIsFasterNoSum() {
    isFaster(maxx, new SingleFixedGaussian(maxx), new DoubleFixedGaussian(maxx), true, true);
  }

  private static void functionsComputeSameValue(int maxx, SinglePrecision f1, DoublePrecision f2,
      final double precision) {
    f1.setMaxX(maxx);
    f2.setMaxX(maxx);
    final float[] p1 = params1.clone();
    final double[] p2 = params2.clone();
    p1[Gaussian.X_POSITION] = (float) (p2[Gaussian.X_POSITION] = (float) (0.123 + maxx / 2));
    p1[Gaussian.Y_POSITION] = (float) (p2[Gaussian.Y_POSITION] = (float) (0.789 + maxx / 2));
    f1.initialise(p1);
    f2.initialise(p2);
    final int n = p1.length;
    final float[] g1 = new float[n];
    final double[] g2 = new double[n];

    double t1 = 0;
    double t2 = 0;
    final double[] tg1 = new double[n];
    final double[] tg2 = new double[n];

    final DoubleDoubleBiPredicate predicate = Predicates.doublesAreClose(precision, 0);

    for (int i = 0, limit = maxx * maxx; i < limit; i++) {
      final float v1 = f1.eval(i);
      t1 += v1;
      final double v2 = f2.eval(i);
      t2 += v2;
      TestAssertions.assertTest(v2, v1, predicate, "Different values");
      final float vv1 = f1.eval(i, g1);
      final double vv2 = f2.eval(i, g2);
      Assertions.assertEquals(v1, vv1, "Different f1 values");
      Assertions.assertEquals(v2, vv2, "Different f2 values");
      for (int j = 0; j < n; j++) {
        tg1[j] += g1[j];
        tg2[j] += g2[j];
      }
      TestAssertions.assertArrayTest(g2, toDouble(g1), predicate, "Different gradients");
    }
    TestAssertions.assertArrayTest(tg2, tg1, predicate, "Different total gradients");
    TestAssertions.assertTest(t2, t1, predicate, "Different totals");
  }

  private static void isFasterWithGradients(int maxx, SinglePrecision f1, DoublePrecision f2,
      boolean noSum, boolean doubleFaster) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    f1.setMaxX(maxx);
    f2.setMaxX(maxx);
    final float[] p1 = params1.clone();
    final double[] p2 = params2.clone();
    p1[Gaussian.X_POSITION] = (float) (p2[Gaussian.X_POSITION] = (float) (0.123 + maxx / 2));
    p1[Gaussian.Y_POSITION] = (float) (p2[Gaussian.Y_POSITION] = (float) (0.789 + maxx / 2));

    long time1;
    long time2;

    if (noSum) {
      time1 = runSingleWithGradientsNoSum(maxx, f1, p1);
      time1 = runSingleWithGradientsNoSum(maxx, f1, p1);
      time1 += runSingleWithGradientsNoSum(maxx, f1, p1);
      time2 = runDoubleWithGradientsNoSum(maxx, f2, p2);
      time2 = runDoubleWithGradientsNoSum(maxx, f2, p2);
      time2 += runDoubleWithGradientsNoSum(maxx, f2, p2);
    } else {
      time1 = runSingleWithGradients(maxx, f1, p1);
      time1 = runSingleWithGradients(maxx, f1, p1);
      time1 += runSingleWithGradients(maxx, f1, p1);
      time2 = runDoubleWithGradients(maxx, f2, p2);
      time2 = runDoubleWithGradients(maxx, f2, p2);
      time2 += runDoubleWithGradients(maxx, f2, p2);
    }

    Class<?> c1;
    Class<?> c2;
    if (doubleFaster) {
      final long time = time1;
      time1 = time2;
      time2 = time;
      c1 = f2.getClass();
      c2 = f1.getClass();
    } else {
      c1 = f1.getClass();
      c2 = f2.getClass();
    }

    logger.log(
        TestLogging.getTimingRecord(((noSum) ? "No sum " : "") + "Gradient " + c1.getSimpleName(),
            time1, c2.getSimpleName(), time2));
  }

  @SuppressWarnings("unused")
  private static long runSingleWithGradients(int maxx, SinglePrecision function, float[] params) {
    function.initialise(params);
    final int n = params1.length;
    final float[] g = new float[n];
    final double[] tg = new double[n];

    final int limit = maxx * maxx;

    // Warm up
    for (int j = 0; j < 10; j++) {
      for (int i = 0; i < limit; i++) {
        function.eval(i, g);
      }
    }

    final long time = System.nanoTime();
    double sum = 0;
    for (int j = 0; j < MAX_ITER; j++) {
      sum = 0;
      for (int i = 0; i < limit; i++) {
        sum += function.eval(i, g);
        for (int k = 0; k < n; k++) {
          tg[k] += g[k];
        }
      }
    }
    return System.nanoTime() - time;
  }

  private static long runSingleWithGradientsNoSum(int maxx, SinglePrecision function,
      float[] params) {
    function.initialise(params);
    final float[] g = new float[params1.length];

    final int limit = maxx * maxx;

    // Warm up
    for (int j = 0; j < 10; j++) {
      for (int i = 0; i < limit; i++) {
        function.eval(i, g);
      }
    }

    final long time = System.nanoTime();
    for (int j = 0; j < MAX_ITER; j++) {
      for (int i = 0; i < limit; i++) {
        function.eval(i, g);
      }
    }
    return System.nanoTime() - time;
  }

  @SuppressWarnings("unused")
  private static long runDoubleWithGradients(int maxx, DoublePrecision function, double[] params) {
    function.initialise(params);
    final int n = params1.length;
    final double[] g = new double[n];
    final double[] tg = new double[n];

    final int limit = maxx * maxx;

    // Warm up
    for (int j = 0; j < 10; j++) {
      for (int i = 0; i < limit; i++) {
        function.eval(i, g);
      }
    }

    final long time = System.nanoTime();
    double sum = 0;
    for (int j = 0; j < MAX_ITER; j++) {
      sum = 0;
      for (int i = 0; i < limit; i++) {
        sum += function.eval(i, g);
        for (int k = 0; k < n; k++) {
          tg[k] += g[k];
        }
      }
    }
    return System.nanoTime() - time;
  }

  private static long runDoubleWithGradientsNoSum(int maxx, DoublePrecision function,
      double[] params) {
    function.initialise(params);
    final double[] g = new double[params1.length];

    final int limit = maxx * maxx;

    // Warm up
    for (int j = 0; j < 10; j++) {
      for (int i = 0; i < limit; i++) {
        function.eval(i, g);
      }
    }

    final long time = System.nanoTime();
    for (int j = 0; j < MAX_ITER; j++) {
      for (int i = 0; i < limit; i++) {
        function.eval(i, g);
      }
    }
    return System.nanoTime() - time;
  }

  private static void isFaster(int maxx, SinglePrecision f1, DoublePrecision f2, boolean noSum,
      boolean doubleFaster) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    f1.setMaxX(maxx);
    f2.setMaxX(maxx);
    final float[] p1 = params1.clone();
    final double[] p2 = params2.clone();
    p1[Gaussian.X_POSITION] = (float) (p2[Gaussian.X_POSITION] = (float) (0.123 + maxx / 2));
    p1[Gaussian.Y_POSITION] = (float) (p2[Gaussian.Y_POSITION] = (float) (0.789 + maxx / 2));

    long time1;
    long time2;
    if (noSum) {
      time1 = runSingleNoSum(maxx, f1, p1);
      time1 = runSingleNoSum(maxx, f1, p1);
      time1 += runSingleNoSum(maxx, f1, p1);
      time2 = runDoubleNoSum(maxx, f2, p2);
      time2 = runDoubleNoSum(maxx, f2, p2);
      time2 += runDoubleNoSum(maxx, f2, p2);
    } else {
      time1 = runSingle(maxx, f1, p1);
      time1 = runSingle(maxx, f1, p1);
      time1 += runSingle(maxx, f1, p1);
      time2 = runDouble(maxx, f2, p2);
      time2 = runDouble(maxx, f2, p2);
      time2 += runDouble(maxx, f2, p2);
    }

    Class<?> c1;
    Class<?> c2;
    if (doubleFaster) {
      final long time = time1;
      time1 = time2;
      time2 = time;
      c1 = f2.getClass();
      c2 = f1.getClass();
    } else {
      c1 = f1.getClass();
      c2 = f2.getClass();
    }

    logger.log(TestLogging.getTimingRecord(((noSum) ? "No sum " : "") + c1.getSimpleName(), time1,
        c2.getSimpleName(), time2));
  }

  @SuppressWarnings("unused")
  private static long runSingle(int maxx, SinglePrecision function, float[] params) {
    final int limit = maxx * maxx;

    // Warm up
    for (int j = 0; j < 10; j++) {
      function.initialise(params);
      for (int i = 0; i < limit; i++) {
        function.eval(i);
      }
    }

    final long time = System.nanoTime();
    double sum = 0;
    for (int j = 0; j < MAX_ITER; j++) {
      sum = 0;
      function.initialise(params);
      for (int i = 0; i < limit; i++) {
        sum += function.eval(i);
      }
    }
    return System.nanoTime() - time;
  }

  private static long runSingleNoSum(int maxx, SinglePrecision function, float[] params) {
    final int limit = maxx * maxx;

    // Warm up
    for (int j = 0; j < 10; j++) {
      function.initialise(params);
      for (int i = 0; i < limit; i++) {
        function.eval(i);
      }
    }

    final long time = System.nanoTime();
    for (int j = 0; j < MAX_ITER; j++) {
      function.initialise(params);
      for (int i = 0; i < limit; i++) {
        function.eval(i);
      }
    }
    return System.nanoTime() - time;
  }

  @SuppressWarnings("unused")
  private static long runDouble(int maxx, DoublePrecision function, double[] params) {
    final int limit = maxx * maxx;

    // Warm up
    for (int j = 0; j < 10; j++) {
      function.initialise(params);
      for (int i = 0; i < limit; i++) {
        function.eval(i);
      }
    }

    final long time = System.nanoTime();
    double sum = 0;
    for (int j = 0; j < MAX_ITER; j++) {
      sum = 0;
      function.initialise(params);
      for (int i = 0; i < limit; i++) {
        sum += function.eval(i);
      }
    }
    return System.nanoTime() - time;
  }

  private static long runDoubleNoSum(int maxx, DoublePrecision function, double[] params) {
    final int limit = maxx * maxx;

    // Warm up
    for (int j = 0; j < 10; j++) {
      function.initialise(params);
      for (int i = 0; i < limit; i++) {
        function.eval(i);
      }
    }

    final long time = System.nanoTime();
    for (int j = 0; j < MAX_ITER; j++) {
      function.initialise(params);
      for (int i = 0; i < limit; i++) {
        function.eval(i);
      }
    }
    return System.nanoTime() - time;
  }

  private static float[] toFloat(double[] data) {
    return SimpleArrayUtils.toFloat(data);
  }

  private static double[] toDouble(float[] data) {
    return SimpleArrayUtils.toDouble(data);
  }
}
