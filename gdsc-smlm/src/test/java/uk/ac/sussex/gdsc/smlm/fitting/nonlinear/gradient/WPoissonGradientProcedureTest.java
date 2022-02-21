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
import java.util.concurrent.ConcurrentHashMap;
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
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.junit5.SpeedTag;
import uk.ac.sussex.gdsc.test.rng.RngUtils;
import uk.ac.sussex.gdsc.test.utils.RandomSeed;
import uk.ac.sussex.gdsc.test.utils.TestComplexity;
import uk.ac.sussex.gdsc.test.utils.TestLogUtils;
import uk.ac.sussex.gdsc.test.utils.TestSettings;
import uk.ac.sussex.gdsc.test.utils.functions.IntArrayFormatSupplier;

@SuppressWarnings({"javadoc"})
class WPoissonGradientProcedureTest {
  private static Logger logger;
  private static ConcurrentHashMap<RandomSeed, double[]> dataCache;

  @BeforeAll
  public static void beforeAll() {
    logger = Logger.getLogger(WPoissonGradientProcedureTest.class.getName());
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

  DoubleEquality eq = new DoubleEquality(1e-6, 1e-16);

  int maxIter = 20000;
  static int blockWidth = 10;
  double background = 0.5;
  double signal = 100;
  double angle = Math.PI;
  double xpos = 5;
  double ypos = 5;
  double xwidth = 1.2;
  double ywidth = 1.2;

  private static double[] createData(RandomSeed source) {
    int count = blockWidth * blockWidth;
    final double[] var = new double[count];
    final UniformRandomProvider rng = RngUtils.create(source.get());
    while (count-- > 0) {
      // Range 0.9 to 1.1
      var[count] = 0.9 + 0.2 * rng.nextDouble();
    }
    return var;
  }

  @SeededTest
  void gradientProcedureFactoryCreatesOptimisedProcedures(RandomSeed seed) {
    final double[] var = dataCache.computeIfAbsent(seed, WPoissonGradientProcedureTest::createData);
    final double[] y = SimpleArrayUtils.newDoubleArray(var.length, 1);
    Assertions.assertEquals(
        WPoissonGradientProcedureUtils.create(y, var, new DummyGradientFunction(6)).getClass(),
        WPoissonGradientProcedure6.class);
    Assertions.assertEquals(
        WPoissonGradientProcedureUtils.create(y, var, new DummyGradientFunction(5)).getClass(),
        WPoissonGradientProcedure5.class);
    Assertions.assertEquals(
        WPoissonGradientProcedureUtils.create(y, var, new DummyGradientFunction(4)).getClass(),
        WPoissonGradientProcedure4.class);
  }

  @SeededTest
  void poissonGradientProcedureComputesSameAsWLsqGradientProcedure(RandomSeed seed) {
    poissonGradientProcedureComputesSameAsWLsqGradientProcedure(seed, 4);
    poissonGradientProcedureComputesSameAsWLsqGradientProcedure(seed, 5);
    poissonGradientProcedureComputesSameAsWLsqGradientProcedure(seed, 6);
    poissonGradientProcedureComputesSameAsWLsqGradientProcedure(seed, 11);
    poissonGradientProcedureComputesSameAsWLsqGradientProcedure(seed, 21);
  }

  private void poissonGradientProcedureComputesSameAsWLsqGradientProcedure(RandomSeed seed,
      int nparams) {
    final double[] var = dataCache.computeIfAbsent(seed, WPoissonGradientProcedureTest::createData);

    final int iter = 10;

    final ArrayList<double[]> paramsList = new ArrayList<>(iter);

    final UniformRandomProvider rng = RngUtils.create(seed.get());
    createFakeParams(rng, nparams, iter, paramsList);
    final FakeGradientFunction func = new FakeGradientFunction(blockWidth, nparams);

    final IntArrayFormatSupplier msgOa =
        getMessage(nparams, "[%d] Observations: Not same alpha @ %d");
    final IntArrayFormatSupplier msgOal =
        getMessage(nparams, "[%d] Observations: Not same alpha linear @ %d");

    for (int i = 0; i < paramsList.size(); i++) {
      final double[] y = createFakeData(rng);
      final WPoissonGradientProcedure p1 = WPoissonGradientProcedureUtils.create(y, var, func);
      p1.computeFisherInformation(paramsList.get(i));
      final WLsqLvmGradientProcedure p2 = new WLsqLvmGradientProcedure(y, var, func);
      p2.gradient(paramsList.get(i));

      // Exactly the same ...
      Assertions.assertArrayEquals(p1.data, p2.alpha, msgOa.set(1, i));
      Assertions.assertArrayEquals(p1.getLinear(), p2.getAlphaLinear(), msgOal.set(1, i));
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
      long time = System.nanoTime();
      run();
      time = System.nanoTime() - time;
      // logger.fine(FunctionUtils.getSupplier("[%d] Time = %d", loops, t);
      return time;
    }

    abstract void run();
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

    final UniformRandomProvider rng = RngUtils.create(seed.get());
    createFakeParams(rng, nparams, iter, paramsList);
    final Gradient1Function func = new FakeGradientFunction(blockWidth, nparams);

    final double[] v =
        (precomputed) ? dataCache.computeIfAbsent(seed, WPoissonGradientProcedureTest::createData)
            : null;

    final IntArrayFormatSupplier msg =
        getMessage(nparams, "[%d] Observations: Not same linear @ %d");
    for (int i = 0; i < paramsList.size(); i++) {
      final double[] y = createFakeData(rng);
      final WPoissonGradientProcedure p1 = new WPoissonGradientProcedure(y, v, func);
      p1.computeFisherInformation(paramsList.get(i));

      final WPoissonGradientProcedure p2 = WPoissonGradientProcedureUtils.create(y, v, func);
      p2.computeFisherInformation(paramsList.get(i));

      // Exactly the same ...
      Assertions.assertArrayEquals(p1.getLinear(), p2.getLinear(), msg.set(1, i));
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
    final ArrayList<double[]> yList = new ArrayList<>(iter);

    createFakeData(RngUtils.create(seed.get()), nparams, iter, paramsList, yList);

    // Remove the timing of the function call by creating a dummy function
    final FakeGradientFunction func = new FakeGradientFunction(blockWidth, nparams);
    final double[] v =
        (precomputed) ? dataCache.computeIfAbsent(seed, WPoissonGradientProcedureTest::createData)
            : null;
    final IntArrayFormatSupplier msg = new IntArrayFormatSupplier("M [%d]", 1);

    for (int i = 0; i < paramsList.size(); i++) {
      final double[] y = yList.get(i);
      final WPoissonGradientProcedure p1 = new WPoissonGradientProcedure(y, v, func);
      p1.computeFisherInformation(paramsList.get(i));

      final WPoissonGradientProcedure p2 = WPoissonGradientProcedureUtils.create(y, v, func);
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
          final WPoissonGradientProcedure p1 = new WPoissonGradientProcedure(yList.get(i), v, func);
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
          final WPoissonGradientProcedure p2 =
              WPoissonGradientProcedureUtils.create(yList.get(i), v, func);
          for (int j = loops; j-- > 0;) {
            p2.computeFisherInformation(paramsList.get(k++ % iter));
          }
        }
      }
    };
    final long time2 = t2.getTime();

    logger.log(TestLogUtils.getTimingRecord("precomputed=" + precomputed + " Standard " + nparams,
        time1, "Unrolled", time2));
  }

  @SpeedTag
  @SeededTest
  void gradientProcedureIsFasterThanWLseGradientProcedure(RandomSeed seed) {
    gradientProcedureIsFasterThanWLseGradientProcedure(seed, 4);
    gradientProcedureIsFasterThanWLseGradientProcedure(seed, 5);
    gradientProcedureIsFasterThanWLseGradientProcedure(seed, 6);
    gradientProcedureIsFasterThanWLseGradientProcedure(seed, 11);
  }

  private void gradientProcedureIsFasterThanWLseGradientProcedure(RandomSeed seed,
      final int nparams) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final int iter = 100;

    final ArrayList<double[]> paramsList = new ArrayList<>(iter);
    final ArrayList<double[]> yList = new ArrayList<>(iter);

    final double[] var = dataCache.computeIfAbsent(seed, WPoissonGradientProcedureTest::createData);
    createFakeData(RngUtils.create(seed.get()), nparams, iter, paramsList, yList);

    // Remove the timing of the function call by creating a dummy function
    final FakeGradientFunction func = new FakeGradientFunction(blockWidth, nparams);
    final IntArrayFormatSupplier msg = new IntArrayFormatSupplier("M [%d]", 1);
    for (int i = 0; i < paramsList.size(); i++) {
      final double[] y = yList.get(i);
      final WLsqLvmGradientProcedure p1 = WLsqLvmGradientProcedureUtils.create(y, var, func);
      p1.gradient(paramsList.get(i));

      final WPoissonGradientProcedure p2 = WPoissonGradientProcedureUtils.create(y, var, func);
      p2.computeFisherInformation(paramsList.get(i));

      // Check they are the same
      Assertions.assertArrayEquals(p1.getAlphaLinear(), p2.getLinear(), msg.set(0, i));
    }

    // Realistic loops for an optimisation
    final int loops = 15;

    // Run till stable timing
    final Timer t1 = new Timer() {
      @Override
      void run() {
        for (int i = 0, k = 0; i < paramsList.size(); i++) {
          final WLsqLvmGradientProcedure p1 =
              WLsqLvmGradientProcedureUtils.create(yList.get(i), var, func);
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
          final WPoissonGradientProcedure p2 =
              WPoissonGradientProcedureUtils.create(yList.get(i), var, func);
          for (int j = loops; j-- > 0;) {
            p2.computeFisherInformation(paramsList.get(k++ % iter));
          }
        }
      }
    };
    final long time2 = t2.getTime();

    logger.log(TestLogUtils.getTimingRecord("WLsqLvmGradientProcedure " + nparams, time1,
        "WPoissonGradientProcedure", time2));
  }

  protected int[] createFakeData(UniformRandomProvider rng, int nparams, int iter,
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

  private static double[] createFakeData(UniformRandomProvider rng, double[] params) {
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

  private static double[] createFakeData(UniformRandomProvider rng) {
    final int n = blockWidth * blockWidth;

    final double[] y = new double[n];
    for (int i = 0; i < y.length; i++) {
      y[i] = rng.nextDouble();
    }

    return y;
  }

  protected void createFakeParams(UniformRandomProvider rng, int nparams, int iter,
      ArrayList<double[]> paramsList) {
    for (int i = 0; i < iter; i++) {
      final double[] params = new double[nparams];
      createFakeParams(rng, params);
      paramsList.add(params);
    }
  }

  private static void createFakeParams(UniformRandomProvider rng, double[] params) {
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
