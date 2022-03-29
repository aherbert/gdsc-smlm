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

package uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient;

import java.util.ArrayList;
import java.util.logging.Logger;
import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.BeforeAll;
import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.smlm.function.DummyGradientFunction;
import uk.ac.sussex.gdsc.smlm.function.FakeGradientFunction;
import uk.ac.sussex.gdsc.smlm.function.Gradient1Function;
import uk.ac.sussex.gdsc.smlm.function.OffsetGradient1Function;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import uk.ac.sussex.gdsc.smlm.function.gaussian.erf.ErfGaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultHelper;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.junit5.SpeedTag;
import uk.ac.sussex.gdsc.test.rng.RngFactory;
import uk.ac.sussex.gdsc.test.utils.RandomSeed;
import uk.ac.sussex.gdsc.test.utils.TestComplexity;
import uk.ac.sussex.gdsc.test.utils.TestLogging;
import uk.ac.sussex.gdsc.test.utils.TestSettings;
import uk.ac.sussex.gdsc.test.utils.functions.IndexSupplier;
import uk.ac.sussex.gdsc.test.utils.functions.IntArrayFormatSupplier;

@SuppressWarnings({"javadoc"})
class LsqVarianceGradientProcedureTest {
  private static Logger logger;

  @BeforeAll
  public static void beforeAll() {
    logger = Logger.getLogger(LsqVarianceGradientProcedureTest.class.getName());
  }

  @AfterAll
  public static void afterAll() {
    logger = null;
  }

  DoubleEquality eq = new DoubleEquality(1e-6, 1e-16);

  int maxIter = 20000;
  int blockWidth = 10;
  double background = 0.5;
  double signal = 100;
  double angle = Math.PI;
  double xpos = 5;
  double ypos = 5;
  double xwidth = 1.2;
  double ywidth = 1.2;

  private static double nextUniform(UniformRandomProvider rng, double min, double max) {
    return min + rng.nextDouble() * (max - min);
  }

  @SeededTest
  void gradientProcedureFactoryCreatesOptimisedProcedures() {
    Assertions.assertEquals(
        LsqVarianceGradientProcedureUtils.create(new DummyGradientFunction(6)).getClass(),
        LsqVarianceGradientProcedure6.class);
    Assertions.assertEquals(
        LsqVarianceGradientProcedureUtils.create(new DummyGradientFunction(5)).getClass(),
        LsqVarianceGradientProcedure5.class);
    Assertions.assertEquals(
        LsqVarianceGradientProcedureUtils.create(new DummyGradientFunction(4)).getClass(),
        LsqVarianceGradientProcedure4.class);
  }

  @SeededTest
  void gradientProcedureComputesSameAsGradientCalculator(RandomSeed seed) {
    gradientProcedureComputesSameAsGradientCalculator(seed, 4);
    gradientProcedureComputesSameAsGradientCalculator(seed, 5);
    gradientProcedureComputesSameAsGradientCalculator(seed, 6);
    gradientProcedureComputesSameAsGradientCalculator(seed, 11);
    gradientProcedureComputesSameAsGradientCalculator(seed, 21);
  }

  private void gradientProcedureComputesSameAsGradientCalculator(RandomSeed seed, int nparams) {
    final int iter = 10;

    final ArrayList<double[]> paramsList = new ArrayList<>(iter);

    createFakeParams(RngFactory.create(seed.get()), nparams, iter, paramsList);
    final int n = blockWidth * blockWidth;
    final FakeGradientFunction func = new FakeGradientFunction(blockWidth, nparams);

    final GradientCalculator calc = GradientCalculatorUtils.newCalculator(nparams, false);

    final IntArrayFormatSupplier msg =
        new IntArrayFormatSupplier("[%d] Observations: Not same variance @ %d", 2);
    msg.set(0, nparams);

    for (int i = 0; i < paramsList.size(); i++) {
      final LsqVarianceGradientProcedure p = LsqVarianceGradientProcedureUtils.create(func);
      p.variance(paramsList.get(i));
      final double[] e = calc.variance(n, paramsList.get(i), func);
      Assertions.assertArrayEquals(e, p.variance, msg.set(1, i));
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
      long time = System.nanoTime();
      run();
      time = System.nanoTime() - time;
      // logger.fine(FormatSupplier.getSupplier("[%d] Time = %d", loops, t);
      return time;
    }

    abstract void run();
  }

  @SeededTest
  void gradientProcedureIsNotSlowerThanGradientCalculator(RandomSeed seed) {
    gradientProcedureIsNotSlowerThanGradientCalculator(seed, 4);
    gradientProcedureIsNotSlowerThanGradientCalculator(seed, 5);
    gradientProcedureIsNotSlowerThanGradientCalculator(seed, 6);
    // 2 peaks
    gradientProcedureIsNotSlowerThanGradientCalculator(seed, 11);
    // 4 peaks
    gradientProcedureIsNotSlowerThanGradientCalculator(seed, 21);
  }

  private void gradientProcedureIsNotSlowerThanGradientCalculator(RandomSeed seed,
      final int nparams) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final int iter = 1000;

    final ArrayList<double[]> paramsList = new ArrayList<>(iter);

    createFakeParams(RngFactory.create(seed.get()), nparams, iter, paramsList);
    final int n = blockWidth * blockWidth;
    final FakeGradientFunction func = new FakeGradientFunction(blockWidth, nparams);

    final GradientCalculator calc = GradientCalculatorUtils.newCalculator(nparams, false);

    for (int i = 0; i < paramsList.size(); i++) {
      calc.variance(n, paramsList.get(i), func);
    }

    for (int i = 0; i < paramsList.size(); i++) {
      final LsqVarianceGradientProcedure p = LsqVarianceGradientProcedureUtils.create(func);
      p.variance(paramsList.get(i));
    }

    // Realistic loops for an optimisation
    final int loops = 15;

    // Run till stable timing
    final Timer t1 = new Timer() {
      @Override
      void run() {
        for (int i = 0, k = 0; i < iter; i++) {
          final GradientCalculator calc = GradientCalculatorUtils.newCalculator(nparams, false);
          for (int j = loops; j-- > 0;) {
            calc.variance(n, paramsList.get(k++ % iter), func);
          }
        }
      }
    };
    final long time1 = t1.getTime();

    final Timer t2 = new Timer(t1.loops) {
      @Override
      void run() {
        for (int i = 0, k = 0; i < iter; i++) {
          final LsqVarianceGradientProcedure p = LsqVarianceGradientProcedureUtils.create(func);
          for (int j = loops; j-- > 0;) {
            p.variance(paramsList.get(k++ % iter));
          }
        }
      }
    };
    final long time2 = t2.getTime();

    logger.log(TestLogging.getTimingRecord("GradientCalculator " + nparams, time1,
        "LSQVarianceGradientProcedure", time2));
  }

  @SeededTest
  public void
      gradientProcedureUnrolledComputesSameAsGradientProcedureWithoutPrecomputed(RandomSeed seed) {
    gradientProcedureUnrolledComputesSameAsGradientProcedure(seed, 4, false);
    gradientProcedureUnrolledComputesSameAsGradientProcedure(seed, 5, false);
    gradientProcedureUnrolledComputesSameAsGradientProcedure(seed, 6, false);
  }

  @SeededTest
  public void
      gradientProcedureUnrolledComputesSameAsGradientProcedureWithPrecomputed(RandomSeed seed) {
    gradientProcedureUnrolledComputesSameAsGradientProcedure(seed, 4, true);
    gradientProcedureUnrolledComputesSameAsGradientProcedure(seed, 5, true);
    gradientProcedureUnrolledComputesSameAsGradientProcedure(seed, 6, true);
  }

  private void gradientProcedureUnrolledComputesSameAsGradientProcedure(RandomSeed seed,
      int nparams, boolean precomputed) {
    final int iter = 10;

    final ArrayList<double[]> paramsList = new ArrayList<>(iter);

    createFakeParams(RngFactory.create(seed.get()), nparams, iter, paramsList);
    Gradient1Function func = new FakeGradientFunction(blockWidth, nparams);

    if (precomputed) {
      func = OffsetGradient1Function.wrapGradient1Function(func,
          SimpleArrayUtils.newArray(func.size(), 0.1, 1.3));
    }

    final IntArrayFormatSupplier msg =
        new IntArrayFormatSupplier("[%d] Observations: Not same variance @ %d", 2);
    msg.set(0, nparams);

    for (int i = 0; i < paramsList.size(); i++) {
      final LsqVarianceGradientProcedure p1 = new LsqVarianceGradientProcedure(func);
      p1.variance(paramsList.get(i));

      final LsqVarianceGradientProcedure p2 = LsqVarianceGradientProcedureUtils.create(func);
      p2.variance(paramsList.get(i));

      // Exactly the same ...
      Assertions.assertArrayEquals(p1.variance, p2.variance, msg.set(1, i));
    }
  }

  @SpeedTag
  @SeededTest
  public void
      gradientProcedureIsFasterUnrolledThanGradientProcedureWithoutPrecomputed(RandomSeed seed) {
    gradientProcedureIsFasterUnrolledThanGradientProcedure(seed, 4, false);
    gradientProcedureIsFasterUnrolledThanGradientProcedure(seed, 5, false);
    gradientProcedureIsFasterUnrolledThanGradientProcedure(seed, 6, false);
  }

  @SpeedTag
  @SeededTest
  public void
      gradientProcedureIsFasterUnrolledThanGradientProcedureWithPrecomputed(RandomSeed seed) {
    gradientProcedureIsFasterUnrolledThanGradientProcedure(seed, 4, true);
    gradientProcedureIsFasterUnrolledThanGradientProcedure(seed, 5, true);
    gradientProcedureIsFasterUnrolledThanGradientProcedure(seed, 6, true);
  }

  private void gradientProcedureIsFasterUnrolledThanGradientProcedure(RandomSeed seed,
      final int nparams, final boolean precomputed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final int iter = 100;

    final ArrayList<double[]> paramsList = new ArrayList<>(iter);

    createFakeParams(RngFactory.create(seed.get()), nparams, iter, paramsList);

    // Remove the timing of the function call by creating a dummy function
    final FakeGradientFunction f = new FakeGradientFunction(blockWidth, nparams);
    final Gradient1Function func =
        (precomputed)
            ? OffsetGradient1Function.wrapGradient1Function(f,
                SimpleArrayUtils.newArray(f.size(), 0.1, 1.3))
            : f;

    final IndexSupplier msg = new IndexSupplier(1, "M ", null);
    for (int i = 0; i < paramsList.size(); i++) {
      final LsqVarianceGradientProcedure p1 = new LsqVarianceGradientProcedure(func);
      p1.variance(paramsList.get(i));
      p1.variance(paramsList.get(i));

      final LsqVarianceGradientProcedure p2 = LsqVarianceGradientProcedureUtils.create(func);
      p2.variance(paramsList.get(i));
      p2.variance(paramsList.get(i));

      // Check they are the same
      Assertions.assertArrayEquals(p1.variance, p2.variance, msg.set(0, i));
    }

    // Realistic loops for an optimisation
    final int loops = 15;

    // Run till stable timing
    final Timer t1 = new Timer() {
      @Override
      void run() {
        for (int i = 0, k = 0; i < paramsList.size(); i++) {
          final LsqVarianceGradientProcedure p1 = new LsqVarianceGradientProcedure(func);
          for (int j = loops; j-- > 0;) {
            p1.variance(paramsList.get(k++ % iter));
          }
        }
      }
    };
    final long time1 = t1.getTime();

    final Timer t2 = new Timer(t1.loops) {
      @Override
      void run() {
        for (int i = 0, k = 0; i < paramsList.size(); i++) {
          final LsqVarianceGradientProcedure p2 = LsqVarianceGradientProcedureUtils.create(func);
          for (int j = loops; j-- > 0;) {
            p2.variance(paramsList.get(k++ % iter));
          }
        }
      }
    };
    final long time2 = t2.getTime();

    logger.log(TestLogging.getTimingRecord("precomputed=" + precomputed + " Standard " + nparams,
        time1, "Unrolled", time2));
  }

  @SeededTest
  void crlbIsHigherWithPrecomputed(RandomSeed seed) {
    final int iter = 10;
    final UniformRandomProvider rng = RngFactory.create(seed.get());

    final ErfGaussian2DFunction func = (ErfGaussian2DFunction) GaussianFunctionFactory.create2D(1,
        10, 10, GaussianFunctionFactory.FIT_ERF_FREE_CIRCLE, null);

    final double[] a = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
    final int n = func.getNumberOfGradients();

    // Get a background
    final double[] b = new double[func.size()];
    for (int i = 0; i < b.length; i++) {
      b[i] = nextUniform(rng, 1, 2);
    }

    for (int i = 0; i < iter; i++) {
      a[Gaussian2DFunction.BACKGROUND] = nextUniform(rng, 0.1, 0.3);
      a[Gaussian2DFunction.SIGNAL] = nextUniform(rng, 100, 300);
      a[Gaussian2DFunction.X_POSITION] = nextUniform(rng, 4, 6);
      a[Gaussian2DFunction.Y_POSITION] = nextUniform(rng, 4, 6);
      a[Gaussian2DFunction.X_SD] = nextUniform(rng, 1, 1.3);
      a[Gaussian2DFunction.Y_SD] = nextUniform(rng, 1, 1.3);

      final LsqVarianceGradientProcedure p1 = LsqVarianceGradientProcedureUtils.create(func);
      p1.variance(a);

      final LsqVarianceGradientProcedure p2 = LsqVarianceGradientProcedureUtils
          .create(OffsetGradient1Function.wrapGradient1Function(func, b));
      p2.variance(a);

      final double[] crlb1 = p1.variance;
      final double[] crlb2 = p2.variance;
      Assertions.assertNotNull(crlb1);
      Assertions.assertNotNull(crlb2);
      // logger.fine(FormatSupplier.getSupplier("%s : %s", Arrays.toString(crlb1),
      // Arrays.toString(crlb2));
      for (int j = 0; j < n; j++) {
        Assertions.assertTrue(crlb1[j] < crlb2[j]);
      }
    }
  }

  @SeededTest
  void varianceMatchesFormula() {
    // Assumptions.assumeTrue(false);

    final double[] n_ = new double[] {20, 50, 100, 500};
    final double[] b2_ = new double[] {0, 1, 2, 4};
    final double[] s_ = new double[] {1, 1.2, 1.5};
    final double[] x_ = new double[] {4.8, 5, 5.5};
    final double a = 100;
    final int size = 10;
    final Gaussian2DFunction f = GaussianFunctionFactory.create2D(1, size, size,
        GaussianFunctionFactory.FIT_ERF_CIRCLE, null);
    final LsqVarianceGradientProcedure p = LsqVarianceGradientProcedureUtils.create(f);
    final int ix = f.findGradientIndex(Gaussian2DFunction.X_POSITION);
    final int iy = f.findGradientIndex(Gaussian2DFunction.Y_POSITION);
    final double[] params = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
    for (final double n : n_) {
      params[Gaussian2DFunction.SIGNAL] = n;
      for (final double b2 : b2_) {
        params[Gaussian2DFunction.BACKGROUND] = b2;
        for (final double s : s_) {
          final double ss = s * a;
          params[Gaussian2DFunction.X_SD] = s;
          for (final double x : x_) {
            params[Gaussian2DFunction.X_POSITION] = x;
            for (final double y : x_) {
              params[Gaussian2DFunction.Y_POSITION] = y;
              if (p.variance(params) != LsqVarianceGradientProcedure.STATUS_OK) {
                Assertions.fail("No variance");
              }
              final double o1 = Math.sqrt(p.variance[ix]) * a;
              final double o2 = Math.sqrt(p.variance[iy]) * a;
              final double e = Gaussian2DPeakResultHelper.getPrecisionX(a, ss, n, b2, false);
              // logger.fine(FormatSupplier.getSupplier("e = %f : o = %f %f", e, o1, o2);
              Assertions.assertEquals(e, o1, e * 5e-2);
              Assertions.assertEquals(e, o2, e * 5e-2);
            }
          }
        }
      }
    }
  }

  protected int[] createFakeData(final UniformRandomProvider rng, int nparams, int iter,
      ArrayList<double[]> paramsList, ArrayList<double[]> yList) {
    final int[] x = new int[blockWidth * blockWidth];
    for (int i = 0; i < x.length; i++) {
      x[i] = i;
    }
    for (int i = 0; i < iter; i++) {
      final double[] params = new double[nparams];
      final double[] y = createFakeData(rng, params);
      paramsList.add(params);
      yList.add(y);
    }
    return x;
  }

  private double[] createFakeData(final UniformRandomProvider rng, double[] params) {
    final int n = blockWidth * blockWidth;

    for (int i = 0; i < params.length; i++) {
      params[i] = rng.nextDouble();
    }

    final double[] y = new double[n];
    for (int i = 0; i < y.length; i++) {
      y[i] = rng.nextDouble();
    }

    return y;
  }

  protected void createFakeParams(final UniformRandomProvider rng, int nparams, int iter,
      ArrayList<double[]> paramsList) {
    for (int i = 0; i < iter; i++) {
      final double[] params = new double[nparams];
      createFakeParams(rng, params);
      paramsList.add(params);
    }
  }

  private static void createFakeParams(final UniformRandomProvider rng, double[] params) {
    for (int i = 0; i < params.length; i++) {
      params[i] = rng.nextDouble();
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
