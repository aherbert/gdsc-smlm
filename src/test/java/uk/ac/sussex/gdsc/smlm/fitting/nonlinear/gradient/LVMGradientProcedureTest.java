package uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient;

import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.RandomGeneratorAdapter;
import uk.ac.sussex.gdsc.core.utils.RandomUtils;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.rng.GaussianSamplerUtils;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.LVMGradientProcedureFactory.Type;
import uk.ac.sussex.gdsc.smlm.function.DummyGradientFunction;
import uk.ac.sussex.gdsc.smlm.function.FakeGradientFunction;
import uk.ac.sussex.gdsc.smlm.function.FastLog;
import uk.ac.sussex.gdsc.smlm.function.FastLogFactory;
import uk.ac.sussex.gdsc.smlm.function.Gradient1Function;
import uk.ac.sussex.gdsc.smlm.function.OffsetGradient1Function;
import uk.ac.sussex.gdsc.smlm.function.TurboLog2;
import uk.ac.sussex.gdsc.smlm.function.ValueProcedure;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import uk.ac.sussex.gdsc.smlm.function.gaussian.erf.ErfGaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.erf.SingleFreeCircularErfGaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.math3.distribution.CustomPoissonDistribution;
import uk.ac.sussex.gdsc.test.api.TestAssertions;
import uk.ac.sussex.gdsc.test.api.TestHelper;
import uk.ac.sussex.gdsc.test.api.function.DoubleDoubleBiPredicate;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.junit5.SpeedTag;
import uk.ac.sussex.gdsc.test.rng.RngUtils;
import uk.ac.sussex.gdsc.test.utils.TestComplexity;
import uk.ac.sussex.gdsc.test.utils.TestCounter;
import uk.ac.sussex.gdsc.test.utils.TestLogUtils;
import uk.ac.sussex.gdsc.test.utils.TestSettings;
import uk.ac.sussex.gdsc.test.utils.TimingResult;
import uk.ac.sussex.gdsc.test.utils.functions.FunctionUtils;
import uk.ac.sussex.gdsc.test.utils.functions.IndexSupplier;
import uk.ac.sussex.gdsc.test.utils.functions.IntArrayFormatSupplier;

import org.apache.commons.math3.util.Precision;
import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.sampling.distribution.GaussianSampler;
import org.ejml.data.DenseMatrix64F;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Disabled;
import org.opentest4j.AssertionFailedError;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.logging.Logger;

/**
 * Contains speed tests for the methods for calculating the Hessian and gradient vector for use in
 * the LVM algorithm.
 */
@SuppressWarnings({"javadoc"})
public class LVMGradientProcedureTest {
  private static Logger logger;

  @BeforeAll
  public static void beforeAll() {
    logger = Logger.getLogger(LVMGradientProcedureTest.class.getName());
  }

  @AfterAll
  public static void afterAll() {
    logger = null;
  }

  static FastLog fastLog = null;

  static FastLog getFastLog() {
    if (fastLog == null) {
      // Default
      fastLog = FastLogFactory.getFastLog();
    }
    return fastLog;
  }

  int MAX_ITER = 20000;
  int blockWidth = 10;
  double noise = 0.3;
  double background = 0.5;
  // High signal required to enabled the SupportsPrecomputed test to distinguish
  // LSQ and MLE/WLSQ. With a value of 100 it is possible to get the same gradient
  // with 3 peaks as with 2 peaks minus the third peak from the data.
  double signal = 1000;
  double angle = Math.PI;
  double xpos = 5;
  double ypos = 5;
  double xwidth = 1.2;
  double ywidth = 1.2;

  private static double random(UniformRandomProvider r, double d) {
    return d - d * 0.1 + r.nextDouble() * 0.2;
  }

  @SeededTest
  public void gradientProcedureFactoryCreatesOptimisedProcedures() {
    final DummyGradientFunction[] f = new DummyGradientFunction[7];
    for (int i = 1; i < f.length; i++) {
      f[i] = new DummyGradientFunction(i);
    }

    final LVMGradientProcedureFactory.Type MLE = LVMGradientProcedureFactory.Type.MLE;
    final LVMGradientProcedureFactory.Type WLSQ = LVMGradientProcedureFactory.Type.WLSQ;
    final LVMGradientProcedureFactory.Type LSQ = LVMGradientProcedureFactory.Type.LSQ;
    final LVMGradientProcedureFactory.Type FMLE = LVMGradientProcedureFactory.Type.FastLogMLE;

    final FastLog fl = getFastLog();

    //@formatter:off

    // Generic factory
    final double[] y0 = new double[1];
    final double[] y1 = new double[]{ 1 };
    Assertions.assertEquals(LVMGradientProcedureFactory.create(y0, f[6], LSQ, fl).getClass(), LSQLVMGradientProcedure6.class);
    Assertions.assertEquals(LVMGradientProcedureFactory.create(y0, f[5], LSQ, fl).getClass(), LSQLVMGradientProcedure5.class);
    Assertions.assertEquals(LVMGradientProcedureFactory.create(y0, f[4], LSQ, fl).getClass(), LSQLVMGradientProcedure4.class);
    Assertions.assertEquals(LVMGradientProcedureFactory.create(y0, f[1], LSQ, fl).getClass(), LSQLVMGradientProcedure.class);
    Assertions.assertEquals(LVMGradientProcedureFactory.create(y0, f[6], MLE, fl).getClass(), MLELVMGradientProcedure6.class);
    Assertions.assertEquals(LVMGradientProcedureFactory.create(y0, f[5], MLE, fl).getClass(), MLELVMGradientProcedure5.class);
    Assertions.assertEquals(LVMGradientProcedureFactory.create(y0, f[4], MLE, fl).getClass(), MLELVMGradientProcedure4.class);
    Assertions.assertEquals(LVMGradientProcedureFactory.create(y0, f[1], MLE, fl).getClass(), MLELVMGradientProcedure.class);
    Assertions.assertEquals(LVMGradientProcedureFactory.create(y1, f[6], MLE, fl).getClass(), MLELVMGradientProcedureX6.class);
    Assertions.assertEquals(LVMGradientProcedureFactory.create(y1, f[5], MLE, fl).getClass(), MLELVMGradientProcedureX5.class);
    Assertions.assertEquals(LVMGradientProcedureFactory.create(y1, f[4], MLE, fl).getClass(), MLELVMGradientProcedureX4.class);
    Assertions.assertEquals(LVMGradientProcedureFactory.create(y1, f[1], MLE, fl).getClass(), MLELVMGradientProcedureX.class);
    Assertions.assertEquals(LVMGradientProcedureFactory.create(y0, f[6], WLSQ, fl).getClass(), WLSQLVMGradientProcedure6.class);
    Assertions.assertEquals(LVMGradientProcedureFactory.create(y0, f[5], WLSQ, fl).getClass(), WLSQLVMGradientProcedure5.class);
    Assertions.assertEquals(LVMGradientProcedureFactory.create(y0, f[4], WLSQ, fl).getClass(), WLSQLVMGradientProcedure4.class);
    Assertions.assertEquals(LVMGradientProcedureFactory.create(y0, f[1], WLSQ, fl).getClass(), WLSQLVMGradientProcedure.class);
    Assertions.assertEquals(LVMGradientProcedureFactory.create(y0, f[6], FMLE, fl).getClass(), FastLogMLELVMGradientProcedure6.class);
    Assertions.assertEquals(LVMGradientProcedureFactory.create(y0, f[5], FMLE, fl).getClass(), FastLogMLELVMGradientProcedure5.class);
    Assertions.assertEquals(LVMGradientProcedureFactory.create(y0, f[4], FMLE, fl).getClass(), FastLogMLELVMGradientProcedure4.class);
    Assertions.assertEquals(LVMGradientProcedureFactory.create(y0, f[1], FMLE, fl).getClass(), FastLogMLELVMGradientProcedure.class);
    Assertions.assertEquals(LVMGradientProcedureFactory.create(y1, f[6], FMLE, fl).getClass(), FastLogMLELVMGradientProcedureX6.class);
    Assertions.assertEquals(LVMGradientProcedureFactory.create(y1, f[5], FMLE, fl).getClass(), FastLogMLELVMGradientProcedureX5.class);
    Assertions.assertEquals(LVMGradientProcedureFactory.create(y1, f[4], FMLE, fl).getClass(), FastLogMLELVMGradientProcedureX4.class);
    Assertions.assertEquals(LVMGradientProcedureFactory.create(y1, f[1], FMLE, fl).getClass(), FastLogMLELVMGradientProcedureX.class);

    // Dedicated factories
    Assertions.assertEquals(LSQLVMGradientProcedureFactory.create(y0, f[6]).getClass(), LSQLVMGradientProcedure6.class);
    Assertions.assertEquals(LSQLVMGradientProcedureFactory.create(y0, f[5]).getClass(), LSQLVMGradientProcedure5.class);
    Assertions.assertEquals(LSQLVMGradientProcedureFactory.create(y0, f[4]).getClass(), LSQLVMGradientProcedure4.class);
    Assertions.assertEquals(LSQLVMGradientProcedureFactory.create(y0, f[1]).getClass(), LSQLVMGradientProcedure.class);
    Assertions.assertEquals(MLELVMGradientProcedureFactory.create(y0, f[6]).getClass(), MLELVMGradientProcedure6.class);
    Assertions.assertEquals(MLELVMGradientProcedureFactory.create(y0, f[5]).getClass(), MLELVMGradientProcedure5.class);
    Assertions.assertEquals(MLELVMGradientProcedureFactory.create(y0, f[4]).getClass(), MLELVMGradientProcedure4.class);
    Assertions.assertEquals(MLELVMGradientProcedureFactory.create(y0, f[1]).getClass(), MLELVMGradientProcedure.class);
    Assertions.assertEquals(MLELVMGradientProcedureFactory.create(y1, f[6]).getClass(), MLELVMGradientProcedureX6.class);
    Assertions.assertEquals(MLELVMGradientProcedureFactory.create(y1, f[5]).getClass(), MLELVMGradientProcedureX5.class);
    Assertions.assertEquals(MLELVMGradientProcedureFactory.create(y1, f[4]).getClass(), MLELVMGradientProcedureX4.class);
    Assertions.assertEquals(MLELVMGradientProcedureFactory.create(y1, f[1]).getClass(), MLELVMGradientProcedureX.class);
    Assertions.assertEquals(WLSQLVMGradientProcedureFactory.create(y0, null, f[6]).getClass(), WLSQLVMGradientProcedure6.class);
    Assertions.assertEquals(WLSQLVMGradientProcedureFactory.create(y0, null, f[5]).getClass(), WLSQLVMGradientProcedure5.class);
    Assertions.assertEquals(WLSQLVMGradientProcedureFactory.create(y0, null, f[4]).getClass(), WLSQLVMGradientProcedure4.class);
    Assertions.assertEquals(WLSQLVMGradientProcedureFactory.create(y0, null, f[1]).getClass(), WLSQLVMGradientProcedure.class);
    Assertions.assertEquals(MLELVMGradientProcedureFactory.create(y0, f[6], fl).getClass(), FastLogMLELVMGradientProcedure6.class);
    Assertions.assertEquals(MLELVMGradientProcedureFactory.create(y0, f[5], fl).getClass(), FastLogMLELVMGradientProcedure5.class);
    Assertions.assertEquals(MLELVMGradientProcedureFactory.create(y0, f[4], fl).getClass(), FastLogMLELVMGradientProcedure4.class);
    Assertions.assertEquals(MLELVMGradientProcedureFactory.create(y0, f[1], fl).getClass(), FastLogMLELVMGradientProcedure.class);
    Assertions.assertEquals(MLELVMGradientProcedureFactory.create(y1, f[6], fl).getClass(), FastLogMLELVMGradientProcedureX6.class);
    Assertions.assertEquals(MLELVMGradientProcedureFactory.create(y1, f[5], fl).getClass(), FastLogMLELVMGradientProcedureX5.class);
    Assertions.assertEquals(MLELVMGradientProcedureFactory.create(y1, f[4], fl).getClass(), FastLogMLELVMGradientProcedureX4.class);
    Assertions.assertEquals(MLELVMGradientProcedureFactory.create(y1, f[1], fl).getClass(), FastLogMLELVMGradientProcedureX.class);

    //@formatter:on
  }

  @SeededTest
  public void gradientProcedureLSQComputesSameAsGradientCalculator(RandomSeed seed) {
    gradientProcedureComputesSameAsGradientCalculator(seed, Type.LSQ);
  }

  @SeededTest
  public void gradientProcedureMLEComputesSameAsGradientCalculator(RandomSeed seed) {
    gradientProcedureComputesSameAsGradientCalculator(seed, Type.MLE);
  }

  @SeededTest
  public void gradientProcedureFastLogMLEComputesSameAsGradientCalculator(RandomSeed seed) {
    gradientProcedureComputesSameAsGradientCalculator(seed, Type.FastLogMLE, 1e-5);
  }

  private void gradientProcedureComputesSameAsGradientCalculator(RandomSeed seed, Type type) {
    gradientProcedureComputesSameAsGradientCalculator(seed, type, 0);
  }

  private void gradientProcedureComputesSameAsGradientCalculator(RandomSeed seed, Type type,
      double error) {
    gradientProcedureComputesSameAsGradientCalculator(seed, 4, type, error);
    gradientProcedureComputesSameAsGradientCalculator(seed, 5, type, error);
    gradientProcedureComputesSameAsGradientCalculator(seed, 6, type, error);
    gradientProcedureComputesSameAsGradientCalculator(seed, 11, type, error);
    gradientProcedureComputesSameAsGradientCalculator(seed, 21, type, error);
  }

  @SpeedTag
  @SeededTest
  public void gradientProcedureLSQIsNotSlowerThanGradientCalculator(RandomSeed seed) {
    gradientProcedureIsNotSlowerThanGradientCalculator(seed, Type.LSQ);
  }

  @SpeedTag
  @SeededTest
  public void gradientProcedureMLEIsNotSlowerThanGradientCalculator(RandomSeed seed) {
    gradientProcedureIsNotSlowerThanGradientCalculator(seed, Type.MLE);
  }

  @SpeedTag
  @SeededTest
  public void gradientProcedureFastLogMLEIsNotSlowerThanGradientCalculator(RandomSeed seed) {
    gradientProcedureIsNotSlowerThanGradientCalculator(seed, Type.FastLogMLE);
  }

  private void gradientProcedureIsNotSlowerThanGradientCalculator(RandomSeed seed, Type type) {
    gradientProcedureIsNotSlowerThanGradientCalculator(seed, 4, type);
    gradientProcedureIsNotSlowerThanGradientCalculator(seed, 5, type);
    gradientProcedureIsNotSlowerThanGradientCalculator(seed, 6, type);
    // 2 peaks
    gradientProcedureIsNotSlowerThanGradientCalculator(seed, 11, type);
    // 4 peaks
    gradientProcedureIsNotSlowerThanGradientCalculator(seed, 21, type);
  }

  private void gradientProcedureComputesSameAsGradientCalculator(RandomSeed seed, int nparams,
      Type type, double error) {
    final int iter = 10;

    final double[][] alpha = new double[nparams][nparams];
    final double[] beta = new double[nparams];

    final ArrayList<double[]> paramsList = new ArrayList<>(iter);
    final ArrayList<double[]> yList = new ArrayList<>(iter);

    final int[] x =
        createFakeData(RngUtils.create(seed.getSeedAsLong()), nparams, iter, paramsList, yList);
    final int n = x.length;
    final FakeGradientFunction func = new FakeGradientFunction(blockWidth, nparams);

    final boolean mle = type != Type.LSQ;
    final FastLog fastLog = (type == Type.FastLogMLE) ? getFastLog() : null;
    final GradientCalculator calc = GradientCalculatorFactory.newCalculator(nparams, mle);

    final String name = String.format("[%d] %b", nparams, mle);
    // Create messages
    final IndexSupplier msgR = new IndexSupplier(1, name + "Result: Not same ", null);
    final IndexSupplier msgOB = new IndexSupplier(1, name + "Observations: Not same beta ", null);
    final IndexSupplier msgOAl =
        new IndexSupplier(1, name + "Observations: Not same alpha linear ", null);
    final IndexSupplier msgOAm =
        new IndexSupplier(1, name + "Observations: Not same alpha matrix ", null);

    final DoubleDoubleBiPredicate predicate =
        (error == 0) ? TestHelper.doublesEqual() : TestHelper.doublesAreClose(error, 0);

    for (int i = 0; i < paramsList.size(); i++) {
      // Reference implementation
      final double s = calc.findLinearised(n, yList.get(i), paramsList.get(i), alpha, beta, func);
      // Procedure
      final LVMGradientProcedure p =
          LVMGradientProcedureFactory.create(yList.get(i), func, type, fastLog);
      p.gradient(paramsList.get(i));
      final double s2 = p.value;
      // Value may be different depending on log implementation
      msgR.set(0, i);
      TestAssertions.assertTest(s, s2, predicate, msgR);

      // Exactly the same ...
      Assertions.assertArrayEquals(p.beta, beta, msgOB.set(0, i));

      final double[] al = p.getAlphaLinear();
      Assertions.assertArrayEquals(al, new DenseMatrix64F(alpha).data, msgOAl.set(0, i));

      final double[][] am = p.getAlphaMatrix();
      Assertions.assertArrayEquals(am, alpha, msgOAm.set(0, i));
    }
  }

  private abstract class Timer {
    private int loops;
    int min;

    Timer() {}

    Timer(int min) {
      this.min = min;
    }

    long getTime() {
      // Run till stable timing
      long t1 = time();
      for (int i = 0; i < 10; i++) {
        final long t2 = t1;
        t1 = time();
        if (loops >= min && DoubleEquality.relativeError(t1, t2) < 0.02) {
          break;
        }
      }
      return t1;
    }

    long time() {
      loops++;
      long t = System.nanoTime();
      run();
      t = System.nanoTime() - t;
      // logger.fine(FunctionUtils.getSupplier("[%d] Time = %d", loops, t);
      return t;
    }

    abstract void run();
  }

  private void gradientProcedureIsNotSlowerThanGradientCalculator(RandomSeed seed,
      final int nparams, final Type type) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final int iter = 1000;
    final double[][] alpha = new double[nparams][nparams];
    final double[] beta = new double[nparams];

    final ArrayList<double[]> paramsList = new ArrayList<>(iter);
    final ArrayList<double[]> yList = new ArrayList<>(iter);

    final int[] x =
        createFakeData(RngUtils.create(seed.getSeedAsLong()), nparams, iter, paramsList, yList);
    final int n = x.length;
    final FakeGradientFunction func = new FakeGradientFunction(blockWidth, nparams);

    final boolean mle = type != Type.LSQ;
    final FastLog fastLog = (type == Type.FastLogMLE) ? getFastLog() : null;
    final GradientCalculator calc = GradientCalculatorFactory.newCalculator(nparams, mle);

    for (int i = 0; i < paramsList.size(); i++) {
      calc.findLinearised(n, yList.get(i), paramsList.get(i), alpha, beta, func);
    }

    for (int i = 0; i < paramsList.size(); i++) {
      final LVMGradientProcedure p =
          LVMGradientProcedureFactory.create(yList.get(i), func, type, fastLog);
      p.gradient(paramsList.get(i));
    }

    // Realistic loops for an optimisation
    final int loops = 15;

    // Run till stable timing
    final Timer t1 = new Timer() {
      @Override
      void run() {
        for (int i = 0, k = 0; i < iter; i++) {
          final GradientCalculator calc = GradientCalculatorFactory.newCalculator(nparams, mle);
          for (int j = loops; j-- > 0;) {
            calc.findLinearised(n, yList.get(i), paramsList.get(k++ % iter), alpha, beta, func);
          }
        }
      }
    };
    final long time1 = t1.getTime();

    final Timer t2 = new Timer(t1.loops) {
      @Override
      void run() {
        for (int i = 0, k = 0; i < iter; i++) {
          final LVMGradientProcedure p =
              LVMGradientProcedureFactory.create(yList.get(i), func, type, fastLog);
          for (int j = loops; j-- > 0;) {
            p.gradient(paramsList.get(k++ % iter));
          }
        }
      }
    };
    final long time2 = t2.getTime();

    logger.log(TestLogUtils.getTimingRecord(new TimingResult("GradientCalculator", time1),
        new TimingResult(() -> String.format("LVMGradientProcedure %d %s", nparams, type), time2)));
  }

  @SeededTest
  public void gradientProcedureLSQUnrolledComputesSameAsGradientProcedure(RandomSeed seed) {
    gradientProcedureUnrolledComputesSameAsGradientProcedure(seed, Type.LSQ, false);
  }

  @SeededTest
  public void gradientProcedureMLEUnrolledComputesSameAsGradientProcedure(RandomSeed seed) {
    gradientProcedureUnrolledComputesSameAsGradientProcedure(seed, Type.MLE, false);
  }

  @SeededTest
  public void gradientProcedureFastLogMLEUnrolledComputesSameAsGradientProcedure(RandomSeed seed) {
    gradientProcedureUnrolledComputesSameAsGradientProcedure(seed, Type.FastLogMLE, false);
  }

  @SeededTest
  public void gradientProcedureWLSQUnrolledComputesSameAsGradientProcedure(RandomSeed seed) {
    gradientProcedureUnrolledComputesSameAsGradientProcedure(seed, Type.WLSQ, false);
  }

  @SeededTest
  public void
      gradientProcedureLSQUnrolledComputesSameAsGradientProcedureWithPrecomputed(RandomSeed seed) {
    gradientProcedureUnrolledComputesSameAsGradientProcedure(seed, Type.LSQ, true);
  }

  @SeededTest
  public void
      gradientProcedureMLEUnrolledComputesSameAsGradientProcedureWithPrecomputed(RandomSeed seed) {
    gradientProcedureUnrolledComputesSameAsGradientProcedure(seed, Type.MLE, true);
  }

  @SeededTest
  public void gradientProcedureFastLogMLEUnrolledComputesSameAsGradientProcedureWithPrecomputed(
      RandomSeed seed) {
    gradientProcedureUnrolledComputesSameAsGradientProcedure(seed, Type.FastLogMLE, true);
  }

  @SeededTest
  public void
      gradientProcedureWLSQUnrolledComputesSameAsGradientProcedureWithPrecomputed(RandomSeed seed) {
    gradientProcedureUnrolledComputesSameAsGradientProcedure(seed, Type.WLSQ, true);
  }

  private void gradientProcedureUnrolledComputesSameAsGradientProcedure(RandomSeed seed, Type type,
      boolean precomputed) {
    gradientProcedureUnrolledComputesSameAsGradientProcedure(seed, 4, type, precomputed);
    gradientProcedureUnrolledComputesSameAsGradientProcedure(seed, 5, type, precomputed);
    gradientProcedureUnrolledComputesSameAsGradientProcedure(seed, 6, type, precomputed);
  }

  private void gradientProcedureUnrolledComputesSameAsGradientProcedure(RandomSeed seed,
      int nparams, Type type, boolean precomputed) {
    final int iter = 10;

    final ArrayList<double[]> paramsList = new ArrayList<>(iter);
    final ArrayList<double[]> yList = new ArrayList<>(iter);

    createFakeData(RngUtils.create(seed.getSeedAsLong()), nparams, iter, paramsList, yList);
    Gradient1Function func = new FakeGradientFunction(blockWidth, nparams);

    if (precomputed) {
      final double[] b = SimpleArrayUtils.newArray(func.size(), 0.1, 1.3);
      func = OffsetGradient1Function.wrapGradient1Function(func, b);
    }

    final FastLog fastLog = type == Type.FastLogMLE ? getFastLog() : null;

    final String name = String.format("[%d] %b", nparams, type);
    // Create messages
    final IndexSupplier msgR = new IndexSupplier(1, name + "Result: Not same ", null);
    final IndexSupplier msgOB = new IndexSupplier(1, name + "Observations: Not same beta ", null);
    final IndexSupplier msgOAl =
        new IndexSupplier(1, name + "Observations: Not same alpha linear ", null);
    final IndexSupplier msgOAm =
        new IndexSupplier(1, name + "Observations: Not same alpha matrix ", null);

    for (int i = 0; i < paramsList.size(); i++) {
      final LVMGradientProcedure p1 = createProcedure(type, yList.get(i), func, fastLog);
      p1.gradient(paramsList.get(i));

      final LVMGradientProcedure p2 =
          LVMGradientProcedureFactory.create(yList.get(i), func, type, fastLog);
      p2.gradient(paramsList.get(i));

      // Exactly the same ...
      Assertions.assertEquals(p1.value, p2.value, msgR.set(0, i));

      Assertions.assertArrayEquals(p1.beta, p2.beta, msgOB.set(0, i));

      Assertions.assertArrayEquals(p1.getAlphaLinear(), p2.getAlphaLinear(), msgOAl.set(0, i));

      final double[][] am1 = p1.getAlphaMatrix();
      final double[][] am2 = p2.getAlphaMatrix();
      Assertions.assertArrayEquals(am1, am2, msgOAm.set(0, i));
    }
  }

  private static LVMGradientProcedure createProcedure(Type type, double[] y, Gradient1Function func,
      FastLog fastLog) {
    switch (type) {
      case FastLogMLE:
        return new FastLogMLELVMGradientProcedure(y, func, fastLog);
      case MLE:
        return new MLELVMGradientProcedure(y, func);
      case WLSQ:
        return new WLSQLVMGradientProcedure(y, null, func);
      default:
        return new LSQLVMGradientProcedure(y, func);
    }
  }

  @SpeedTag
  @SeededTest
  public void gradientProcedureLSQIsFasterUnrolledThanGradientProcedure(RandomSeed seed) {
    gradientProcedureIsFasterUnrolledThanGradientProcedure(seed, Type.LSQ, false);
  }

  @SpeedTag
  @SeededTest
  public void gradientProcedureMLEIsFasterUnrolledThanGradientProcedure(RandomSeed seed) {
    gradientProcedureIsFasterUnrolledThanGradientProcedure(seed, Type.MLE, false);
  }

  @SpeedTag
  @SeededTest
  public void gradientProcedureFastLogMLEIsFasterUnrolledThanGradientProcedure(RandomSeed seed) {
    gradientProcedureIsFasterUnrolledThanGradientProcedure(seed, Type.FastLogMLE, false);
  }

  @SpeedTag
  @SeededTest
  public void gradientProcedureWLSQIsFasterUnrolledThanGradientProcedure(RandomSeed seed) {
    gradientProcedureIsFasterUnrolledThanGradientProcedure(seed, Type.WLSQ, false);
  }

  @SpeedTag
  @SeededTest
  public void
      gradientProcedureLSQIsFasterUnrolledThanGradientProcedureWithPrecomputed(RandomSeed seed) {
    gradientProcedureIsFasterUnrolledThanGradientProcedure(seed, Type.LSQ, true);
  }

  @SpeedTag
  @SeededTest
  public void
      gradientProcedureMLEIsFasterUnrolledThanGradientProcedureWithPrecomputed(RandomSeed seed) {
    gradientProcedureIsFasterUnrolledThanGradientProcedure(seed, Type.MLE, true);
  }

  @SpeedTag
  @SeededTest
  public void gradientProcedureFastLogMLEIsFasterUnrolledThanGradientProcedureWithPrecomputed(
      RandomSeed seed) {
    gradientProcedureIsFasterUnrolledThanGradientProcedure(seed, Type.FastLogMLE, true);
  }

  @SpeedTag
  @SeededTest
  public void
      gradientProcedureWLSQIsFasterUnrolledThanGradientProcedureWithPrecomputed(RandomSeed seed) {
    gradientProcedureIsFasterUnrolledThanGradientProcedure(seed, Type.WLSQ, true);
  }

  private void gradientProcedureIsFasterUnrolledThanGradientProcedure(RandomSeed seed, Type type,
      boolean precomputed) {
    gradientProcedureIsFasterUnrolledThanGradientProcedure(seed, 4, type, precomputed);
    gradientProcedureIsFasterUnrolledThanGradientProcedure(seed, 5, type, precomputed);
    gradientProcedureIsFasterUnrolledThanGradientProcedure(seed, 6, type, precomputed);
  }

  private void gradientProcedureIsFasterUnrolledThanGradientProcedure(RandomSeed seed,
      final int nparams, final Type type, final boolean precomputed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final int iter = 100;

    final ArrayList<double[]> paramsList = new ArrayList<>(iter);
    final ArrayList<double[]> yList = new ArrayList<>(iter);

    createData(RngUtils.create(seed.getSeedAsLong()), 1, iter, paramsList, yList);

    // Remove the timing of the function call by creating a dummy function
    final FakeGradientFunction fgf = new FakeGradientFunction(blockWidth, nparams);
    final Gradient1Function func;
    if (precomputed) {
      final double[] b = SimpleArrayUtils.newArray(fgf.size(), 0.1, 1.3);
      func = OffsetGradient1Function.wrapGradient1Function(fgf, b);
    } else {
      func = fgf;
    }

    final FastLog fastLog = type == Type.FastLogMLE ? getFastLog() : null;

    final IntArrayFormatSupplier msgA = new IntArrayFormatSupplier("A [%d]", 1);
    final IntArrayFormatSupplier msgB = new IntArrayFormatSupplier("B [%d]", 1);

    for (int i = 0; i < paramsList.size(); i++) {
      final LVMGradientProcedure p1 = createProcedure(type, yList.get(i), func, fastLog);
      p1.gradient(paramsList.get(i));
      p1.gradient(paramsList.get(i));

      final LVMGradientProcedure p2 =
          LVMGradientProcedureFactory.create(yList.get(i), func, type, fastLog);
      p2.gradient(paramsList.get(i));
      p2.gradient(paramsList.get(i));

      // Check they are the same
      Assertions.assertArrayEquals(p1.getAlphaLinear(), p2.getAlphaLinear(), msgA.set(0, i));
      Assertions.assertArrayEquals(p1.beta, p2.beta, msgB.set(0, i));
    }

    // Realistic loops for an optimisation
    final int loops = 15;

    // Run till stable timing
    final Timer t1 = new Timer() {
      @Override
      void run() {
        for (int i = 0, k = 0; i < paramsList.size(); i++) {
          final LVMGradientProcedure p1 = createProcedure(type, yList.get(i), func, fastLog);
          for (int j = loops; j-- > 0;) {
            p1.gradient(paramsList.get(k++ % iter));
          }
        }
      }
    };
    final long time1 = t1.getTime();

    final Timer t2 = new Timer(t1.loops) {
      @Override
      void run() {
        for (int i = 0, k = 0; i < paramsList.size(); i++) {
          final LVMGradientProcedure p2 =
              LVMGradientProcedureFactory.create(yList.get(i), func, type, fastLog);
          for (int j = loops; j-- > 0;) {
            p2.gradient(paramsList.get(k++ % iter));
          }
        }
      }
    };
    final long time2 = t2.getTime();

    logger
        .log(TestLogUtils.getTimingRecord(
            new TimingResult(
                () -> String.format("%s, Precomputed=%b : Standard", type, precomputed), time1),
            new TimingResult(() -> String.format("Unrolled %d", nparams), time2)));
  }

  @SeededTest
  public void gradientProcedureLSQComputesGradient(RandomSeed seed) {
    gradientProcedureComputesGradient(seed,
        new SingleFreeCircularErfGaussian2DFunction(blockWidth, blockWidth), Type.LSQ, false);
  }

  @SeededTest
  public void gradientProcedureMLEComputesGradient(RandomSeed seed) {
    gradientProcedureComputesGradient(seed,
        new SingleFreeCircularErfGaussian2DFunction(blockWidth, blockWidth), Type.MLE, false);
  }

  @SeededTest
  public void gradientProcedureFastLogMLECannotComputeGradient(RandomSeed seed) {
    Assertions.assertThrows(AssertionFailedError.class, () -> {
      gradientProcedureComputesGradient(seed,
          new SingleFreeCircularErfGaussian2DFunction(blockWidth, blockWidth), Type.FastLogMLE,
          false);
    });
  }

  @Disabled("This test now passes as the tolerance for computing the gradient has been lowered "
      + " so that the tests pass under a stress test using many different random seeds.")
  @SeededTest
  public void gradientProcedureFastLogMLECannotComputeGradientWithHighPrecision(RandomSeed seed) {
    // Try different precision
    for (int n = FastLog.N; n < 23; n++) {
      try {
        // logger.fine(FunctionUtils.getSupplier("Precision n=%d", n);
        fastLog = new TurboLog2(n);
        gradientProcedureComputesGradient(seed,
            new SingleFreeCircularErfGaussian2DFunction(blockWidth, blockWidth), Type.FastLogMLE,
            false);
      } catch (final AssertionError ex) {
        continue;
      } finally {
        // Reset
        fastLog = null;
      }
      return;
    }
    Assertions.fail();
  }

  @SeededTest
  public void gradientProcedureWLSQComputesGradient(RandomSeed seed) {
    gradientProcedureComputesGradient(seed,
        new SingleFreeCircularErfGaussian2DFunction(blockWidth, blockWidth), Type.WLSQ, false);
  }

  @SeededTest
  public void gradientProcedureLSQComputesGradientWithPrecomputed(RandomSeed seed) {
    gradientProcedureComputesGradient(seed,
        new SingleFreeCircularErfGaussian2DFunction(blockWidth, blockWidth), Type.LSQ, true);
  }

  @SeededTest
  public void gradientProcedureMLEComputesGradientWithPrecomputed(RandomSeed seed) {
    gradientProcedureComputesGradient(seed,
        new SingleFreeCircularErfGaussian2DFunction(blockWidth, blockWidth), Type.MLE, true);
  }

  @SeededTest
  public void gradientProcedureFastLogMLECannotComputeGradientWithPrecomputed(RandomSeed seed) {
    Assertions.assertThrows(AssertionFailedError.class, () -> {
      gradientProcedureComputesGradient(seed,
          new SingleFreeCircularErfGaussian2DFunction(blockWidth, blockWidth), Type.FastLogMLE,
          true);
    });
  }

  @SeededTest
  public void gradientProcedureWLSQComputesGradientWithPrecomputed(RandomSeed seed) {
    gradientProcedureComputesGradient(seed,
        new SingleFreeCircularErfGaussian2DFunction(blockWidth, blockWidth), Type.WLSQ, true);
  }

  @SuppressWarnings("null")
  private void gradientProcedureComputesGradient(RandomSeed seed, ErfGaussian2DFunction func,
      Type type, boolean precomputed) {
    final int nparams = func.getNumberOfGradients();
    final int[] indices = func.gradientIndices();

    final int iter = 100;

    final ArrayList<double[]> paramsList = new ArrayList<>(iter);
    final ArrayList<double[]> yList = new ArrayList<>(iter);

    createData(RngUtils.create(seed.getSeedAsLong()), 1, iter, paramsList, yList, true);

    // for the gradients
    final double delta = 1e-4;
    final DoubleEquality eq = new DoubleEquality(5e-2, 1e-16);

    final double[] b = (precomputed) ? new double[func.size()] : null;

    final FastLog fastLog = type == Type.FastLogMLE ? getFastLog() : null;

    // Must compute most of the time
    final int failureLimit = TestCounter.computeFailureLimit(iter, 0.1);
    final TestCounter failCounter = new TestCounter(failureLimit, nparams);

    for (int i = 0; i < paramsList.size(); i++) {
      final int ii = i;
      final double[] y = yList.get(i);
      final double[] a = paramsList.get(i);
      final double[] a2 = a.clone();

      LVMGradientProcedure p;
      if (precomputed) {
        // Mock fitting part of the function already
        for (int j = 0; j < b.length; j++) {
          b[j] = y[j] * 0.5;
        }
        p = LVMGradientProcedureFactory.create(y,
            OffsetGradient1Function.wrapGradient1Function(func, b), type, fastLog);
      } else {
        p = LVMGradientProcedureFactory.create(y, func, type, fastLog);
      }
      p.gradient(a);
      // double s = p.value;
      final double[] beta = p.beta.clone();
      for (int j = 0; j < nparams; j++) {
        final int jj = j;
        final int k = indices[j];
        // double d = Precision.representableDelta(a[k], (a[k] == 0) ? 1e-3 : a[k] * delta);
        final double d = Precision.representableDelta(a[k], delta);
        a2[k] = a[k] + d;
        p.value(a2);
        final double s1 = p.value;
        a2[k] = a[k] - d;
        p.value(a2);
        final double s2 = p.value;
        a2[k] = a[k];

        // Apply a factor of -2 to compute the actual gradients:
        // See Numerical Recipes in C++, 2nd Ed. Equation 15.5.6 for Nonlinear Models
        beta[j] *= -2;

        final double gradient = (s1 - s2) / (2 * d);
        // logger.fine(FunctionUtils.getSupplier("[%d,%d] %f (%s %f+/-%f) %f ?= %f", i, k, s,
        // Gaussian2DFunction.getName(k),
        // a[k], d, beta[j], gradient);

        failCounter.run(j, () -> eq.almostEqualRelativeOrAbsolute(beta[jj], gradient), () -> {
          Assertions.fail(() -> String.format("Not same gradient @ %d,%d: %s != %s (error=%s)", ii,
              jj, beta[jj], gradient, DoubleEquality.relativeError(beta[jj], gradient)));
        });
      }
    }
  }

  @SeededTest
  public void gradientProcedureLSQSupportsPrecomputed(RandomSeed seed) {
    gradientProcedureSupportsPrecomputed(seed, Type.LSQ);
  }

  @SeededTest
  public void gradientProcedureMLESupportsPrecomputed(RandomSeed seed) {
    gradientProcedureSupportsPrecomputed(seed, Type.MLE);
  }

  @SeededTest
  public void gradientProcedureFastLogMLESupportsPrecomputed(RandomSeed seed) {
    gradientProcedureSupportsPrecomputed(seed, Type.FastLogMLE, false);
  }

  @SeededTest
  public void gradientProcedureFastLogMLECannotSupportPrecomputedWithGradients(RandomSeed seed) {
    Assertions.assertThrows(AssertionFailedError.class, () -> {
      gradientProcedureSupportsPrecomputed(seed, Type.FastLogMLE);
    });
  }

  @SeededTest
  public void gradientProcedureWLSQSupportsPrecomputed(RandomSeed seed) {
    gradientProcedureSupportsPrecomputed(seed, Type.WLSQ);
  }

  private void gradientProcedureSupportsPrecomputed(RandomSeed seed, final Type type) {
    gradientProcedureSupportsPrecomputed(seed, type, true);
  }

  private void gradientProcedureSupportsPrecomputed(RandomSeed seed, final Type type,
      boolean checkGradients) {
    final int iter = 10;
    final UniformRandomProvider r = RngUtils.create(seed.getSeedAsLong());
    final GaussianSampler gs = GaussianSamplerUtils.createGaussianSampler(r, 0, noise);

    final ArrayList<double[]> paramsList = new ArrayList<>(iter);
    final ArrayList<double[]> yList = new ArrayList<>(iter);

    // 3 peaks
    createData(r, 3, iter, paramsList, yList, true);

    for (int i = 0; i < paramsList.size(); i++) {
      final double[] y = yList.get(i);
      // Add Gaussian read noise so we have negatives
      final double min = MathUtils.min(y);
      for (int j = 0; j < y.length; j++) {
        y[j] = y[i] - min + gs.sample();
      }
    }

    // We want to know that:
    // y|peak1+peak2+peak3 == y|peak1+peak2+peak3(precomputed)
    // We want to know when:
    // y|peak1+peak2+peak3 != y-peak3|peak1+peak2
    // i.e. we cannot subtract a precomputed peak from the data, it must be included in the fit
    // E.G. LSQ - subtraction is OK, MLE/WLSQ - subtraction is not allowed

    final Gaussian2DFunction f123 = GaussianFunctionFactory.create2D(3, blockWidth, blockWidth,
        GaussianFunctionFactory.FIT_ERF_FREE_CIRCLE, null);
    final Gaussian2DFunction f12 = GaussianFunctionFactory.create2D(2, blockWidth, blockWidth,
        GaussianFunctionFactory.FIT_ERF_FREE_CIRCLE, null);
    final Gaussian2DFunction f3 = GaussianFunctionFactory.create2D(1, blockWidth, blockWidth,
        GaussianFunctionFactory.FIT_ERF_FREE_CIRCLE, null);

    final FastLog fastLog = type == Type.FastLogMLE ? getFastLog() : null;

    final int nparams = f12.getNumberOfGradients();
    final int[] indices = f12.gradientIndices();
    final double[] b = new double[f12.size()];

    final DoubleEquality eq = new DoubleEquality(1e-8, 1e-16); // for checking strict equivalence

    // for the gradients
    final double delta = 1e-4;
    final DoubleEquality eq2 = new DoubleEquality(5e-2, 1e-16);

    final double[] a1peaks = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
    final double[] y_b = new double[b.length];

    // Count the number of failures for each gradient
    final int failureLimit = TestCounter.computeFailureLimit(iter, 0.1);
    final TestCounter failCounter = new TestCounter(failureLimit, nparams * 2);

    for (int i = 0; i < paramsList.size(); i++) {
      final int ii = i;
      final double[] y = yList.get(i);
      final double[] a3peaks = paramsList.get(i);
      // logger.fine(FunctionUtils.getSupplier("[%d] a=%s", i, Arrays.toString(a3peaks));

      final double[] a2peaks =
          Arrays.copyOf(a3peaks, 1 + 2 * Gaussian2DFunction.PARAMETERS_PER_PEAK);
      final double[] a2peaks2 = a2peaks.clone();
      for (int j = 1; j < a1peaks.length; j++) {
        a1peaks[j] = a3peaks[j + 2 * Gaussian2DFunction.PARAMETERS_PER_PEAK];
      }

      // Evaluate peak 3 to get the background and subtract it from the data to get the new data
      f3.initialise0(a1peaks);
      f3.forEach(new ValueProcedure() {
        int k = 0;

        @Override
        public void execute(double value) {
          b[k] = value;
          // Remove negatives for MLE
          if (type.isMLE()) {
            y[k] = Math.max(0, y[k]);
            y_b[k] = Math.max(0, y[k] - value);
          } else {
            y_b[k] = y[k] - value;
          }
          k++;
        }
      });

      final LVMGradientProcedure p123 = LVMGradientProcedureFactory.create(y, f123, type, fastLog);

      /////////////////////////////////////
      // These should be the same
      /////////////////////////////////////
      final LVMGradientProcedure p12b3 = LVMGradientProcedureFactory.create(y,
          OffsetGradient1Function.wrapGradient1Function(f12, b), type, fastLog);

      // Check they are the same
      p123.gradient(a3peaks);
      final double[][] m123 = p123.getAlphaMatrix();

      p12b3.gradient(a2peaks);
      double s = p12b3.value;
      final double[] beta = p12b3.beta.clone();
      double[][] alpha = p12b3.getAlphaMatrix();

      // logger.fine(FunctionUtils.getSupplier("MLE=%b [%d] p12b3 %f %f", type.isMLE(), i,
      // p123.value, s);

      if (!eq.almostEqualRelativeOrAbsolute(p123.value, s)) {
        Assertions.fail(FunctionUtils.getSupplier("p12b3 Not same value @ %d (error=%s) : %s == %s",
            i, DoubleEquality.relativeError(p123.value, s), p123.value, s));
      }
      if (!eq.almostEqualRelativeOrAbsolute(beta, p123.beta)) {
        Assertions
            .fail(FunctionUtils.getSupplier("p12b3 Not same gradient @ %d (error=%s) : %s vs %s", i,
                DoubleEquality.relativeError(beta, p123.beta), Arrays.toString(beta),
                Arrays.toString(p123.beta)));
      }
      for (int j = 0; j < alpha.length; j++) {
        // logger.fine(FunctionUtils.getSupplier("%s !=\n%s", Arrays.toString(alpha[j]),
        // Arrays.toString(m123[j]));
        if (!eq.almostEqualRelativeOrAbsolute(alpha[j], m123[j])) {
          Assertions
              .fail(FunctionUtils.getSupplier("p12b3 Not same alpha @ %d,%d (error=%s) : %s vs %s",
                  i, j, DoubleEquality.relativeError(alpha[j], m123[j]), Arrays.toString(alpha[j]),
                  Arrays.toString(m123[j])));
        }
      }

      // Check actual gradients are correct
      if (checkGradients) {
        for (int j = 0; j < nparams; j++) {
          final int jj = j;
          final int k = indices[j];
          // double d = Precision.representableDelta(a2peaks[k], (a2peaks[k] == 0) ? 1e-3 :
          // a2peaks[k] * delta);
          final double d = Precision.representableDelta(a2peaks[k], delta);
          a2peaks2[k] = a2peaks[k] + d;
          p12b3.value(a2peaks2);
          final double s1 = p12b3.value;
          a2peaks2[k] = a2peaks[k] - d;
          p12b3.value(a2peaks2);
          final double s2 = p12b3.value;
          a2peaks2[k] = a2peaks[k];

          // Apply a factor of -2 to compute the actual gradients:
          // See Numerical Recipes in C++, 2nd Ed. Equation 15.5.6 for Nonlinear Models
          beta[j] *= -2;

          final double gradient = (s1 - s2) / (2 * d);
          // logger.fine(FunctionUtils.getSupplier("[%d,%d] %f (%s %f+/-%f) %f ?= %f (%f)", i, k, s,
          // Gaussian2DFunction.getName(k), a2peaks[k], d, beta[j], gradient,
          // DoubleEquality.relativeError(gradient, beta[j]));
          failCounter.run(j, () -> eq2.almostEqualRelativeOrAbsolute(beta[jj], gradient), () -> {
            Assertions.fail(() -> String.format("Not same gradient @ %d,%d: %s != %s (error=%s)",
                ii, jj, beta[jj], gradient, DoubleEquality.relativeError(beta[jj], gradient)));
          });
        }
      }

      /////////////////////////////////////
      // This may be different
      /////////////////////////////////////
      final LVMGradientProcedure p12m3 =
          LVMGradientProcedureFactory.create(y_b, f12, type, fastLog);

      // Check these may be different.
      // Sometimes they are not different.

      p12m3.gradient(a2peaks);
      s = p12m3.value;
      System.arraycopy(p12m3.beta, 0, beta, 0, p12m3.beta.length);
      alpha = p12m3.getAlphaMatrix();

      // logger.fine(FunctionUtils.getSupplier("%s [%d] p12m3 %f %f", type, i, p123.value, s);

      // The test for different or equal is not robust to different random seeds.
      // ExtraAssertions.fail has been changed for TestLog.logFailure

      if (type != Type.LSQ) {
        if (eq.almostEqualRelativeOrAbsolute(p123.value, s)) {
          logger.log(TestLogUtils.getFailRecord("p12b3 Same value @ %d (error=%s) : %s == %s", i,
              DoubleEquality.relativeError(p123.value, s), p123.value, s));
        }
        if (eq.almostEqualRelativeOrAbsolute(beta, p123.beta)) {
          logger.log(TestLogUtils.getFailRecord("p12b3 Same gradient @ %d (error=%s) : %s vs %s", i,
              DoubleEquality.relativeError(beta, p123.beta), Arrays.toString(beta),
              Arrays.toString(p123.beta)));
        }

        // Note: Test the matrix is different by finding 1 different column
        int dj = -1;
        for (int j = 0; j < alpha.length; j++) {
          // logger.fine(FunctionUtils.getSupplier("%s !=\n%s\n", Arrays.toString(alpha[j]),
          // Arrays.toString(m123[j]));
          if (!eq.almostEqualRelativeOrAbsolute(alpha[j], m123[j])) {
            dj = j; // Different column
            break;
          }
        }
        if (dj == -1) {
          // Find biggest error for reporting. This helps set the test tolerance.
          double error = 0;
          dj = -1;
          for (int j = 0; j < alpha.length; j++) {
            final double e = DoubleEquality.relativeError(alpha[j], m123[j]);
            if (error <= e) {
              error = e;
              dj = j;
            }
          }
          logger.log(TestLogUtils.getFailRecord("p12b3 Same alpha @ %d,%d (error=%s) : %s vs %s", i,
              dj, error, Arrays.toString(alpha[dj]), Arrays.toString(m123[dj])));
        }
      } else {
        if (!eq.almostEqualRelativeOrAbsolute(p123.value, s)) {
          logger.log(TestLogUtils.getFailRecord("p12b3 Not same value @ %d (error=%s) : %s == %s",
              i, DoubleEquality.relativeError(p123.value, s), p123.value, s));
        }
        if (!eq.almostEqualRelativeOrAbsolute(beta, p123.beta)) {
          logger
              .log(TestLogUtils.getFailRecord("p12b3 Not same gradient @ %d (error=%s) : %s vs %s",
                  i, DoubleEquality.relativeError(beta, p123.beta), Arrays.toString(beta),
                  Arrays.toString(p123.beta)));
        }
        for (int j = 0; j < alpha.length; j++) {
          // logger.fine(FunctionUtils.getSupplier("%s !=\n%s\n", Arrays.toString(alpha[j]),
          // Arrays.toString(m123[j]));
          if (!eq.almostEqualRelativeOrAbsolute(alpha[j], m123[j])) {
            logger.log(
                TestLogUtils.getFailRecord("p12b3 Not same alpha @ %d,%d (error=%s) : %s vs %s", i,
                    j, DoubleEquality.relativeError(alpha[j], m123[j]), Arrays.toString(alpha[j]),
                    Arrays.toString(m123[j])));
          }
        }
      }

      // Check actual gradients are correct
      if (!checkGradients) {
        continue;
      }

      for (int j = 0; j < nparams; j++) {
        final int jj = j;
        final int k = indices[j];
        // double d = Precision.representableDelta(a2peaks[k], (a2peaks[k] == 0) ? 1e-3 : a2peaks[k]
        // * delta);
        final double d = Precision.representableDelta(a2peaks[k], delta);
        a2peaks2[k] = a2peaks[k] + d;
        p12m3.value(a2peaks2);
        final double s1 = p12m3.value;
        a2peaks2[k] = a2peaks[k] - d;
        p12m3.value(a2peaks2);
        final double s2 = p12m3.value;
        a2peaks2[k] = a2peaks[k];

        // Apply a factor of -2 to compute the actual gradients:
        // See Numerical Recipes in C++, 2nd Ed. Equation 15.5.6 for Nonlinear Models
        beta[j] *= -2;

        final double gradient = (s1 - s2) / (2 * d);
        // logger.fine(FunctionUtils.getSupplier("[%d,%d] %f (%s %f+/-%f) %f ?= %f (%f)", i, k, s,
        // Gaussian2DFunction.getName(k), a2peaks[k], d, beta[j], gradient,
        // DoubleEquality.relativeError(gradient, beta[j]));
        failCounter.run(nparams + j, () -> eq2.almostEqualRelativeOrAbsolute(beta[jj], gradient),
            () -> {
              Assertions.fail(() -> String.format("Not same gradient @ %d,%d: %s != %s (error=%s)",
                  ii, jj, beta[jj], gradient, DoubleEquality.relativeError(beta[jj], gradient)));
            });
      }
    }
  }

  /**
   * Create random elliptical Gaussian data an returns the data plus an estimate of the parameters.
   * Only the chosen parameters are randomised and returned for a maximum of (background, amplitude,
   * angle, xpos, ypos, xwidth, ywidth }
   *
   * @param npeaks the npeaks
   * @param params set on output
   * @param randomiseParams Set to true to randomise the params
   * @return the double[]
   */
  private double[] doubleCreateGaussianData(UniformRandomProvider r, int npeaks, double[] params,
      boolean randomiseParams) {
    final int n = blockWidth * blockWidth;

    // Generate a 2D Gaussian
    final ErfGaussian2DFunction func =
        (ErfGaussian2DFunction) GaussianFunctionFactory.create2D(npeaks, blockWidth, blockWidth,
            GaussianFunctionFactory.FIT_ERF_FREE_CIRCLE, null);
    params[0] = random(r, background);
    for (int i = 0, j = 0; i < npeaks; i++, j += Gaussian2DFunction.PARAMETERS_PER_PEAK) {
      params[j + Gaussian2DFunction.SIGNAL] = random(r, signal);
      params[j + Gaussian2DFunction.X_POSITION] = random(r, xpos);
      params[j + Gaussian2DFunction.Y_POSITION] = random(r, ypos);
      params[j + Gaussian2DFunction.X_SD] = random(r, xwidth);
      params[j + Gaussian2DFunction.Y_SD] = random(r, ywidth);
    }

    if (npeaks > 1) {
      // Move the peaks around so they do not overlap
      final double[] shift = SimpleArrayUtils.newArray(npeaks, -2, 4.0 / (npeaks - 1));
      RandomUtils.shuffle(shift, r);
      for (int i = 0, j = 0; i < npeaks; i++, j += Gaussian2DFunction.PARAMETERS_PER_PEAK) {
        params[j + Gaussian2DFunction.X_POSITION] += shift[i];
      }
      RandomUtils.shuffle(shift, r);
      for (int i = 0, j = 0; i < npeaks; i++, j += Gaussian2DFunction.PARAMETERS_PER_PEAK) {
        params[j + Gaussian2DFunction.Y_POSITION] += shift[i];
      }
    }

    final double[] y = new double[n];
    final CustomPoissonDistribution pd =
        new CustomPoissonDistribution(new RandomGeneratorAdapter(r), 1);
    func.initialise(params);
    for (int i = 0; i < y.length; i++) {
      // Add random Poisson noise
      pd.setMeanUnsafe(func.eval(i));
      y[i] = pd.sample();
    }

    if (randomiseParams) {
      params[0] = random(r, params[0]);
      for (int i = 0, j = 0; i < npeaks; i++, j += Gaussian2DFunction.PARAMETERS_PER_PEAK) {
        params[j + Gaussian2DFunction.SIGNAL] = random(r, params[j + Gaussian2DFunction.SIGNAL]);
        params[j + Gaussian2DFunction.X_POSITION] =
            random(r, params[j + Gaussian2DFunction.X_POSITION]);
        params[j + Gaussian2DFunction.Y_POSITION] =
            random(r, params[j + Gaussian2DFunction.Y_POSITION]);
        params[j + Gaussian2DFunction.X_SD] = random(r, params[j + Gaussian2DFunction.X_SD]);
        params[j + Gaussian2DFunction.Y_SD] = random(r, params[j + Gaussian2DFunction.Y_SD]);
      }
    }

    return y;
  }

  protected int[] createData(UniformRandomProvider r, int npeaks, int iter,
      ArrayList<double[]> paramsList, ArrayList<double[]> yList) {
    return createData(r, npeaks, iter, paramsList, yList, true);
  }

  protected int[] createData(UniformRandomProvider r, int npeaks, int iter,
      ArrayList<double[]> paramsList, ArrayList<double[]> yList, boolean randomiseParams) {
    final int[] x = new int[blockWidth * blockWidth];
    for (int i = 0; i < x.length; i++) {
      x[i] = i;
    }
    for (int i = 0; i < iter; i++) {
      final double[] params = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK * npeaks];
      final double[] y = doubleCreateGaussianData(r, npeaks, params, randomiseParams);
      paramsList.add(params);
      yList.add(y);
    }
    return x;
  }

  protected int[] createFakeData(UniformRandomProvider r, int nparams, int iter,
      ArrayList<double[]> paramsList, ArrayList<double[]> yList) {
    final int[] x = new int[blockWidth * blockWidth];
    for (int i = 0; i < x.length; i++) {
      x[i] = i;
    }
    for (int i = 0; i < iter; i++) {
      final double[] params = new double[nparams];
      final double[] y = createFakeData(r, params);
      paramsList.add(params);
      yList.add(y);
    }
    return x;
  }

  private double[] createFakeData(UniformRandomProvider r, double[] params) {
    final int n = blockWidth * blockWidth;

    for (int i = 0; i < params.length; i++) {
      params[i] = r.nextDouble();
    }

    final double[] y = new double[n];
    for (int i = 0; i < y.length; i++) {
      y[i] = r.nextDouble();
    }

    return y;
  }

  protected ArrayList<double[]> copyList(ArrayList<double[]> paramsList) {
    final ArrayList<double[]> params2List = new ArrayList<>(paramsList.size());
    for (int i = 0; i < paramsList.size(); i++) {
      params2List.add(paramsList.get(i).clone());
    }
    return params2List;
  }
}
