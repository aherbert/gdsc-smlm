package uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient;

import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.smlm.fitting.FisherInformationMatrix;
import uk.ac.sussex.gdsc.smlm.function.DummyGradientFunction;
import uk.ac.sussex.gdsc.smlm.function.FakeGradientFunction;
import uk.ac.sussex.gdsc.smlm.function.Gradient1Function;
import uk.ac.sussex.gdsc.smlm.function.OffsetGradient1Function;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import uk.ac.sussex.gdsc.smlm.function.gaussian.erf.ErfGaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultHelper;
import uk.ac.sussex.gdsc.test.api.TestAssertions;
import uk.ac.sussex.gdsc.test.api.TestHelper;
import uk.ac.sussex.gdsc.test.api.function.DoubleDoubleBiPredicate;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.junit5.SpeedTag;
import uk.ac.sussex.gdsc.test.rng.RngUtils;
import uk.ac.sussex.gdsc.test.utils.TestComplexity;
import uk.ac.sussex.gdsc.test.utils.TestLogUtils;
import uk.ac.sussex.gdsc.test.utils.TestSettings;
import uk.ac.sussex.gdsc.test.utils.functions.IntArrayFormatSupplier;

import org.apache.commons.rng.UniformRandomProvider;
import org.ejml.data.DenseMatrix64F;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.BeforeAll;

import java.util.ArrayList;
import java.util.logging.Logger;

@SuppressWarnings({"javadoc"})
public class PoissonGradientProcedureTest {
  private static Logger logger;

  @BeforeAll
  public static void beforeAll() {
    logger = Logger.getLogger(PoissonGradientProcedureTest.class.getName());
  }

  @AfterAll
  public static void afterAll() {
    logger = null;
  }

  DoubleEquality eq = new DoubleEquality(1e-6, 1e-16);

  int MAX_ITER = 20000;
  int blockWidth = 10;
  double background = 0.5;
  double signal = 100;
  double angle = Math.PI;
  double xpos = 5;
  double ypos = 5;
  double xwidth = 1.2;
  double ywidth = 1.2;

  private static double nextUniform(UniformRandomProvider r, double min, double max) {
    return min + r.nextDouble() * (max - min);
  }

  @SeededTest
  public void gradientProcedureFactoryCreatesOptimisedProcedures() {
    Assertions.assertEquals(
        PoissonGradientProcedureFactory.create(new DummyGradientFunction(6)).getClass(),
        PoissonGradientProcedure6.class);
    Assertions.assertEquals(
        PoissonGradientProcedureFactory.create(new DummyGradientFunction(5)).getClass(),
        PoissonGradientProcedure5.class);
    Assertions.assertEquals(
        PoissonGradientProcedureFactory.create(new DummyGradientFunction(4)).getClass(),
        PoissonGradientProcedure4.class);
  }

  @SeededTest
  public void gradientProcedureComputesSameAsGradientCalculator(RandomSeed seed) {
    gradientProcedureComputesSameAsGradientCalculator(seed, 4);
    gradientProcedureComputesSameAsGradientCalculator(seed, 5);
    gradientProcedureComputesSameAsGradientCalculator(seed, 6);
    gradientProcedureComputesSameAsGradientCalculator(seed, 11);
    gradientProcedureComputesSameAsGradientCalculator(seed, 21);
  }

  @SpeedTag
  @SeededTest
  public void gradientProcedureIsNotSlowerThanGradientCalculator(RandomSeed seed) {
    gradientProcedureIsNotSlowerThanGradientCalculator(seed, 4);
    gradientProcedureIsNotSlowerThanGradientCalculator(seed, 5);
    gradientProcedureIsNotSlowerThanGradientCalculator(seed, 6);
    // 2 peaks
    gradientProcedureIsNotSlowerThanGradientCalculator(seed, 11);
    // 4 peaks
    gradientProcedureIsNotSlowerThanGradientCalculator(seed, 21);
  }

  private void gradientProcedureComputesSameAsGradientCalculator(RandomSeed seed, int nparams) {
    final int iter = 10;
    final ArrayList<double[]> paramsList = new ArrayList<>(iter);

    createFakeParams(RngUtils.create(seed.getSeedAsLong()), nparams, iter, paramsList);
    final int n = blockWidth * blockWidth;
    final FakeGradientFunction func = new FakeGradientFunction(blockWidth, nparams);

    final GradientCalculator calc = GradientCalculatorFactory.newCalculator(nparams, false);

    // Create messages
    final IntArrayFormatSupplier msgOAl =
        getMessage(nparams, "[%d] Observations: Not same alpha linear @ %d");
    final IntArrayFormatSupplier msgOAm =
        getMessage(nparams, "[%d] Observations: Not same alpha matrix @ %d");

    final DoubleDoubleBiPredicate predicate = TestHelper.doublesAreClose(1e-10, 0);

    for (int i = 0; i < paramsList.size(); i++) {
      final PoissonGradientProcedure p = PoissonGradientProcedureFactory.create(func);
      p.computeFisherInformation(paramsList.get(i));
      final double[][] m = calc.fisherInformationMatrix(n, paramsList.get(i), func);
      // Not exactly the same ...
      final double[] al = p.getLinear();
      TestAssertions.assertArrayTest(al, new DenseMatrix64F(m).data, predicate, msgOAl.set(1, i));

      final double[][] am = p.getMatrix();
      TestAssertions.assertArrayTest(am, m, predicate, msgOAm.set(1, i));
    }
  }

  private static IntArrayFormatSupplier getMessage(int nparams, String format) {
    final IntArrayFormatSupplier msg = new IntArrayFormatSupplier(format, 2);
    msg.set(0, nparams);
    return msg;
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
      final int nparams) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final int iter = 1000;
    final ArrayList<double[]> paramsList = new ArrayList<>(iter);

    createFakeParams(RngUtils.create(seed.getSeedAsLong()), nparams, iter, paramsList);
    final int n = blockWidth * blockWidth;
    final FakeGradientFunction func = new FakeGradientFunction(blockWidth, nparams);

    final GradientCalculator calc = GradientCalculatorFactory.newCalculator(nparams, false);

    for (int i = 0; i < paramsList.size(); i++) {
      calc.fisherInformationMatrix(n, paramsList.get(i), func);
    }

    for (int i = 0; i < paramsList.size(); i++) {
      final PoissonGradientProcedure p = PoissonGradientProcedureFactory.create(func);
      p.computeFisherInformation(paramsList.get(i));
    }

    // Realistic loops for an optimisation
    final int loops = 15;

    // Run till stable timing
    final Timer t1 = new Timer() {
      @Override
      void run() {
        for (int i = 0, k = 0; i < iter; i++) {
          final GradientCalculator calc = GradientCalculatorFactory.newCalculator(nparams, false);
          for (int j = loops; j-- > 0;) {
            calc.fisherInformationMatrix(n, paramsList.get(k++ % iter), func);
          }
        }
      }
    };
    final long time1 = t1.getTime();

    final Timer t2 = new Timer(t1.loops) {
      @Override
      void run() {
        for (int i = 0, k = 0; i < iter; i++) {
          final PoissonGradientProcedure p = PoissonGradientProcedureFactory.create(func);
          for (int j = loops; j-- > 0;) {
            p.computeFisherInformation(paramsList.get(k++ % iter));
          }
        }
      }
    };
    final long time2 = t2.getTime();

    logger.log(TestLogUtils.getTimingRecord("GradientCalculator " + nparams, time1,
        "PoissonGradientProcedure", time2));
  }

  @SeededTest
  public void gradientProcedureUnrolledComputesSameAsGradientProcedure(RandomSeed seed) {
    gradientProcedureUnrolledComputesSameAsGradientProcedure(seed, 4, false);
    gradientProcedureUnrolledComputesSameAsGradientProcedure(seed, 5, false);
    gradientProcedureUnrolledComputesSameAsGradientProcedure(seed, 6, false);
  }

  @SeededTest
  public void gradientProcedureUnrolledComputesSameAsGradientProcedureWithPrecomputed(
      RandomSeed seed) {
    gradientProcedureUnrolledComputesSameAsGradientProcedure(seed, 4, true);
    gradientProcedureUnrolledComputesSameAsGradientProcedure(seed, 5, true);
    gradientProcedureUnrolledComputesSameAsGradientProcedure(seed, 6, true);
  }

  private void gradientProcedureUnrolledComputesSameAsGradientProcedure(RandomSeed seed,
      int nparams, boolean precomputed) {
    final int iter = 10;
    final ArrayList<double[]> paramsList = new ArrayList<>(iter);

    createFakeParams(RngUtils.create(seed.getSeedAsLong()), nparams, iter, paramsList);
    Gradient1Function func = new FakeGradientFunction(blockWidth, nparams);

    if (precomputed) {
      func = OffsetGradient1Function.wrapGradient1Function(func,
          SimpleArrayUtils.newArray(func.size(), 0.1, 1.3));
    }

    // Create messages
    final IntArrayFormatSupplier msgOAl =
        getMessage(nparams, "[%d] Observations: Not same alpha linear @ %d");
    final IntArrayFormatSupplier msgOAm =
        getMessage(nparams, "[%d] Observations: Not same alpha matrix @ %d");

    for (int i = 0; i < paramsList.size(); i++) {
      final PoissonGradientProcedure p1 = new PoissonGradientProcedure(func);
      p1.computeFisherInformation(paramsList.get(i));

      final PoissonGradientProcedure p2 = PoissonGradientProcedureFactory.create(func);
      p2.computeFisherInformation(paramsList.get(i));

      // Exactly the same ...
      Assertions.assertArrayEquals(p1.getLinear(), p2.getLinear(), msgOAl.set(1, i));

      final double[][] am1 = p1.getMatrix();
      final double[][] am2 = p2.getMatrix();
      Assertions.assertArrayEquals(am1, am2, msgOAm.set(1, i));
    }
  }

  @SpeedTag
  @SeededTest
  public void gradientProcedureIsFasterUnrolledThanGradientProcedure(RandomSeed seed) {
    gradientProcedureIsFasterUnrolledThanGradientProcedure(seed, 4, false);
    gradientProcedureIsFasterUnrolledThanGradientProcedure(seed, 5, false);
    gradientProcedureIsFasterUnrolledThanGradientProcedure(seed, 6, false);
  }

  @SpeedTag
  @SeededTest
  public void gradientProcedureIsFasterUnrolledThanGradientProcedureWithPrecomputed(
      RandomSeed seed) {
    gradientProcedureIsFasterUnrolledThanGradientProcedure(seed, 4, true);
    gradientProcedureIsFasterUnrolledThanGradientProcedure(seed, 5, true);
    gradientProcedureIsFasterUnrolledThanGradientProcedure(seed, 6, true);
  }

  private void gradientProcedureIsFasterUnrolledThanGradientProcedure(RandomSeed seed,
      final int nparams, final boolean precomputed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final int iter = 100;
    final ArrayList<double[]> paramsList = new ArrayList<>(iter);

    createFakeParams(RngUtils.create(seed.getSeedAsLong()), nparams, iter, paramsList);

    // Remove the timing of the function call by creating a dummy function
    final FakeGradientFunction f = new FakeGradientFunction(blockWidth, nparams);
    final Gradient1Function func =
        (precomputed)
            ? OffsetGradient1Function.wrapGradient1Function(f,
                SimpleArrayUtils.newArray(f.size(), 0.1, 1.3))
            : f;

    final IntArrayFormatSupplier msg = new IntArrayFormatSupplier("M [%d]", 1);
    for (int i = 0; i < paramsList.size(); i++) {
      final PoissonGradientProcedure p1 = new PoissonGradientProcedure(func);
      p1.computeFisherInformation(paramsList.get(i));
      p1.computeFisherInformation(paramsList.get(i));

      final PoissonGradientProcedure p2 = PoissonGradientProcedureFactory.create(func);
      p2.computeFisherInformation(paramsList.get(i));
      p2.computeFisherInformation(paramsList.get(i));

      // Check they are the same
      Assertions.assertArrayEquals(p1.getLinear(), p2.getLinear(), msg.set(0, i));
    }

    // Realistic loops for an optimisation
    final int loops = 15;

    // Run till stable timing
    final Timer t1 = new Timer() {
      @Override
      void run() {
        for (int i = 0, k = 0; i < paramsList.size(); i++) {
          final PoissonGradientProcedure p1 = new PoissonGradientProcedure(func);
          for (int j = loops; j-- > 0;) {
            p1.computeFisherInformation(paramsList.get(k++ % iter));
          }
        }
      }
    };
    final long time1 = t1.getTime();

    final Timer t2 = new Timer(t1.loops) {
      @Override
      void run() {
        for (int i = 0, k = 0; i < paramsList.size(); i++) {
          final PoissonGradientProcedure p2 = PoissonGradientProcedureFactory.create(func);
          for (int j = loops; j-- > 0;) {
            p2.computeFisherInformation(paramsList.get(k++ % iter));
          }
        }
      }
    };
    final long time2 = t2.getTime();

    logger.log(TestLogUtils.getTimingRecord("precomputed=" + precomputed + " Standard " + nparams,
        time1, "Unrolled", time2));
    // Assertions.assertTrue(time2 < time1);
  }

  @SeededTest
  public void crlbIsHigherWithPrecomputed(RandomSeed seed) {
    final int iter = 10;
    final UniformRandomProvider r = RngUtils.create(seed.getSeedAsLong());

    final ErfGaussian2DFunction func = (ErfGaussian2DFunction) GaussianFunctionFactory.create2D(1,
        10, 10, GaussianFunctionFactory.FIT_ERF_FREE_CIRCLE, null);

    final double[] a = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
    final int n = func.getNumberOfGradients();

    // Get a background
    final double[] b = new double[func.size()];
    for (int i = 0; i < b.length; i++) {
      b[i] = nextUniform(r, 1, 2);
    }

    for (int i = 0; i < iter; i++) {
      a[Gaussian2DFunction.BACKGROUND] = nextUniform(r, 0.1, 0.3);
      a[Gaussian2DFunction.SIGNAL] = nextUniform(r, 100, 300);
      a[Gaussian2DFunction.X_POSITION] = nextUniform(r, 4, 6);
      a[Gaussian2DFunction.Y_POSITION] = nextUniform(r, 4, 6);
      a[Gaussian2DFunction.X_SD] = nextUniform(r, 1, 1.3);
      a[Gaussian2DFunction.Y_SD] = nextUniform(r, 1, 1.3);

      final PoissonGradientProcedure p1 = PoissonGradientProcedureFactory.create(func);
      p1.computeFisherInformation(a);

      final PoissonGradientProcedure p2 = PoissonGradientProcedureFactory
          .create(OffsetGradient1Function.wrapGradient1Function(func, b));
      p2.computeFisherInformation(a);

      final FisherInformationMatrix m1 = new FisherInformationMatrix(p1.getLinear(), n);
      final FisherInformationMatrix m2 = new FisherInformationMatrix(p2.getLinear(), n);
      final double[] crlb1 = m1.crlb();
      final double[] crlb2 = m2.crlb();
      Assertions.assertNotNull(crlb1);
      Assertions.assertNotNull(crlb2);
      // logger.fine(FunctionUtils.getSupplier("%s : %s", Arrays.toString(crlb1),
      // Arrays.toString(crlb2));
      for (int j = 0; j < n; j++) {
        Assertions.assertTrue(crlb1[j] < crlb2[j]);
      }
    }
  }

  @SeededTest
  public void varianceMatchesFormula() {
    // Assumptions.assumeTrue(false);

    final double[] N_ = new double[] {20, 50, 100, 500};
    final double[] b2_ = new double[] {0, 1, 2, 4};
    final double[] s_ = new double[] {1, 1.2, 1.5};
    final double[] x_ = new double[] {4.8, 5, 5.5};
    final double a = 100;
    final int size = 10;
    final Gaussian2DFunction f = GaussianFunctionFactory.create2D(1, size, size,
        GaussianFunctionFactory.FIT_ERF_CIRCLE, null);
    final PoissonGradientProcedure p = PoissonGradientProcedureFactory.create(f);
    final int ix = f.findGradientIndex(Gaussian2DFunction.X_POSITION);
    final int iy = f.findGradientIndex(Gaussian2DFunction.Y_POSITION);
    final double[] params = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
    for (final double N : N_) {
      params[Gaussian2DFunction.SIGNAL] = N;
      for (final double b2 : b2_) {
        params[Gaussian2DFunction.BACKGROUND] = b2;
        for (final double s : s_) {
          final double ss = s * a;
          params[Gaussian2DFunction.X_SD] = s;
          for (final double x : x_) {
            params[Gaussian2DFunction.X_POSITION] = x;
            for (final double y : x_) {
              params[Gaussian2DFunction.Y_POSITION] = y;
              p.computeFisherInformation(params);
              final FisherInformationMatrix m1 = new FisherInformationMatrix(p.getLinear(), p.n);
              final double[] crlb = m1.crlb();
              if (crlb == null) {
                Assertions.fail("No variance");
              }
              @SuppressWarnings("null")
              final double o1 = Math.sqrt(crlb[ix]) * a;
              final double o2 = Math.sqrt(crlb[iy]) * a;
              final double e = Gaussian2DPeakResultHelper.getMLPrecisionX(a, ss, N, b2, false);
              // logger.fine(FunctionUtils.getSupplier("e = %f : o = %f %f", e, o1, o2);
              Assertions.assertEquals(e, o1, e * 5e-2);
              Assertions.assertEquals(e, o2, e * 5e-2);
            }
          }
        }
      }
    }
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

  protected void createFakeParams(UniformRandomProvider r, int nparams, int iter,
      ArrayList<double[]> paramsList) {
    for (int i = 0; i < iter; i++) {
      final double[] params = new double[nparams];
      createFakeParams(r, params);
      paramsList.add(params);
    }
  }

  private static void createFakeParams(UniformRandomProvider r, double[] params) {
    for (int i = 0; i < params.length; i++) {
      params[i] = r.nextDouble();
    }
  }

  protected ArrayList<double[]> copyList(ArrayList<double[]> paramsList) {
    final ArrayList<double[]> params2List = new ArrayList<>(paramsList.size());
    for (int i = 0; i < paramsList.size(); i++) {
      params2List.add(paramsList.get(i).clone());
    }
    return params2List;
  }
}
