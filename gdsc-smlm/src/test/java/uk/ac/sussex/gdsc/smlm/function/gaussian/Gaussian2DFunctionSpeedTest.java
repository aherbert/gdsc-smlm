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

package uk.ac.sussex.gdsc.smlm.function.gaussian;

import java.util.ArrayList;
import java.util.concurrent.ConcurrentHashMap;
import java.util.logging.Logger;
import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.BeforeAll;
import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.junit5.SpeedTag;
import uk.ac.sussex.gdsc.test.rng.RngUtils;
import uk.ac.sussex.gdsc.test.utils.RandomSeed;
import uk.ac.sussex.gdsc.test.utils.TestComplexity;
import uk.ac.sussex.gdsc.test.utils.TestLogUtils;
import uk.ac.sussex.gdsc.test.utils.TestSettings;

/**
 * Contains speed tests for the fastest method for calculating the Hessian and gradient vector from
 * a Gaussian 2D Function.
 */
@SuppressWarnings({"javadoc"})
class Gaussian2DFunctionSpeedTest {
  private static Logger logger;
  private static ConcurrentHashMap<RandomSeed, Object> dataCache;

  @BeforeAll
  public static void beforeAll() {
    logger = Logger.getLogger(Gaussian2DFunctionSpeedTest.class.getName());
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

  private final int single = 1;
  private final int multi = 2;

  private static int blockWidth = 10;
  private static double background = 20;
  private static double amplitude = 10;
  private static double xpos = 5;
  private static double ypos = 5;
  private static double xwidth = 5;

  private static class Gaussian2DFunctionSpeedTestData {
    ArrayList<double[]> paramsListSinglePeak = new ArrayList<>();
    ArrayList<double[]> yListSinglePeak = new ArrayList<>();
    ArrayList<double[]> paramsListMultiPeak = new ArrayList<>();
    ArrayList<double[]> yListMultiPeak = new ArrayList<>();
    final UniformRandomProvider rand;

    Gaussian2DFunctionSpeedTestData(UniformRandomProvider rand) {
      this.rand = rand;
    }
  }

  private static Object createData(RandomSeed source) {
    return new Gaussian2DFunctionSpeedTestData(RngUtils.create(source.get()));
  }

  // private static ArrayList<double[]> paramsListSinglePeak = new ArrayList<>();
  // private static ArrayList<double[]> yListSinglePeak = new ArrayList<>();
  // private static ArrayList<double[]> paramsListMultiPeak = new ArrayList<>();
  // private static ArrayList<double[]> yListMultiPeak = new ArrayList<>();

  private static int[] x;

  static {
    x = new int[blockWidth * blockWidth];
    for (int i = 0; i < x.length; i++) {
      x[i] = i;
    }
  }

  private static Gaussian2DFunctionSpeedTestData ensureDataSingle(RandomSeed seed, int size) {
    final Gaussian2DFunctionSpeedTestData data = (Gaussian2DFunctionSpeedTestData) dataCache
        .computeIfAbsent(seed, Gaussian2DFunctionSpeedTest::createData);
    if (data.paramsListSinglePeak.size() < size) {
      synchronized (data.paramsListSinglePeak) {
        if (data.paramsListSinglePeak.size() < size) {
          createGaussianData(data.rand, 1, size, data.paramsListSinglePeak, data.yListSinglePeak);
        }
      }
    }
    return data;
  }

  private static Gaussian2DFunctionSpeedTestData ensureDataMulti(RandomSeed seed, int size) {
    final Gaussian2DFunctionSpeedTestData data = (Gaussian2DFunctionSpeedTestData) dataCache
        .computeIfAbsent(seed, Gaussian2DFunctionSpeedTest::createData);
    if (data.paramsListMultiPeak.size() < size) {
      synchronized (data.paramsListMultiPeak) {
        if (data.paramsListMultiPeak.size() < size) {
          createGaussianData(data.rand, 2, size, data.paramsListMultiPeak, data.yListMultiPeak);
        }
      }
    }
    return data;
  }

  @SeededTest
  void freeCircularComputesSameAsEllipticalSinglePeak(RandomSeed seed) {
    f1ComputesSameAsf2(seed, single, GaussianFunctionFactory.FIT_FREE_CIRCLE,
        GaussianFunctionFactory.FIT_ELLIPTICAL);
  }

  @SpeedTag
  @SeededTest
  void freeCircularFasterThanEllipticalSinglePeak(RandomSeed seed) {
    f1FasterThanf2(seed, single, GaussianFunctionFactory.FIT_FREE_CIRCLE,
        GaussianFunctionFactory.FIT_ELLIPTICAL);
  }

  @SeededTest
  void circularComputesSameAsFreeCircularSinglePeak(RandomSeed seed) {
    f1ComputesSameAsf2(seed, single, GaussianFunctionFactory.FIT_CIRCLE,
        GaussianFunctionFactory.FIT_FREE_CIRCLE);
  }

  @SpeedTag
  @SeededTest
  void circularFasterThanFreeCircularSinglePeak(RandomSeed seed) {
    f1FasterThanf2(seed, single, GaussianFunctionFactory.FIT_CIRCLE,
        GaussianFunctionFactory.FIT_FREE_CIRCLE);
  }

  @SeededTest
  void fixedComputesSameAsFreeCircularSinglePeak(RandomSeed seed) {
    f1ComputesSameAsf2(seed, single, GaussianFunctionFactory.FIT_FIXED,
        GaussianFunctionFactory.FIT_FREE_CIRCLE);
  }

  @SpeedTag
  @SeededTest
  void fixedFasterThanFreeCircularSinglePeak(RandomSeed seed) {
    f1FasterThanf2(seed, single, GaussianFunctionFactory.FIT_FIXED,
        GaussianFunctionFactory.FIT_FREE_CIRCLE);
  }

  @SeededTest
  void freeCircularComputesSameAsEllipticalSinglePeakNb(RandomSeed seed) {
    f1ComputesSameAsf2(seed, single, GaussianFunctionFactory.FIT_SIMPLE_NB_FREE_CIRCLE,
        GaussianFunctionFactory.FIT_SIMPLE_NB_ELLIPTICAL);
  }

  @SpeedTag
  @SeededTest
  void freeCircularFasterThanEllipticalSinglePeakNb(RandomSeed seed) {
    f1FasterThanf2(seed, single, GaussianFunctionFactory.FIT_SIMPLE_NB_FREE_CIRCLE,
        GaussianFunctionFactory.FIT_SIMPLE_NB_ELLIPTICAL);
  }

  @SeededTest
  void circularComputesSameAsFreeCircularSinglePeakNb(RandomSeed seed) {
    f1ComputesSameAsf2(seed, single, GaussianFunctionFactory.FIT_SIMPLE_NB_CIRCLE,
        GaussianFunctionFactory.FIT_SIMPLE_NB_FREE_CIRCLE);
  }

  @SpeedTag
  @SeededTest
  void circularFasterThanFreeCircularSinglePeakNb(RandomSeed seed) {
    f1FasterThanf2(seed, single, GaussianFunctionFactory.FIT_SIMPLE_NB_CIRCLE,
        GaussianFunctionFactory.FIT_SIMPLE_NB_FREE_CIRCLE);
  }

  @SeededTest
  void fixedComputesSameAsFreeCircularSinglePeakNb(RandomSeed seed) {
    f1ComputesSameAsf2(seed, single, GaussianFunctionFactory.FIT_SIMPLE_NB_FIXED,
        GaussianFunctionFactory.FIT_SIMPLE_NB_FREE_CIRCLE);
  }

  @SpeedTag
  @SeededTest
  void fixedFasterThanFreeCircularSinglePeakNb(RandomSeed seed) {
    f1FasterThanf2(seed, single, GaussianFunctionFactory.FIT_SIMPLE_NB_FIXED,
        GaussianFunctionFactory.FIT_SIMPLE_NB_FREE_CIRCLE);
  }

  @SeededTest
  void freeCircularComputesSameAsEllipticalMultiPeak(RandomSeed seed) {
    f1ComputesSameAsf2(seed, multi, GaussianFunctionFactory.FIT_FREE_CIRCLE,
        GaussianFunctionFactory.FIT_ELLIPTICAL);
  }

  @SpeedTag
  @SeededTest
  void freeCircularFasterThanEllipticalMultiPeak(RandomSeed seed) {
    f1FasterThanf2(seed, multi, GaussianFunctionFactory.FIT_FREE_CIRCLE,
        GaussianFunctionFactory.FIT_ELLIPTICAL);
  }

  @SeededTest
  void circularComputesSameAsFreeCircularMultiPeak(RandomSeed seed) {
    f1ComputesSameAsf2(seed, multi, GaussianFunctionFactory.FIT_CIRCLE,
        GaussianFunctionFactory.FIT_FREE_CIRCLE);
  }

  @SpeedTag
  @SeededTest
  void circularFasterThanFreeCircularMultiPeak(RandomSeed seed) {
    f1FasterThanf2(seed, multi, GaussianFunctionFactory.FIT_CIRCLE,
        GaussianFunctionFactory.FIT_FREE_CIRCLE);
  }

  @SeededTest
  void fixedComputesSameAsFreeCircularMultiPeak(RandomSeed seed) {
    f1ComputesSameAsf2(seed, multi, GaussianFunctionFactory.FIT_FIXED,
        GaussianFunctionFactory.FIT_FREE_CIRCLE);
  }

  @SpeedTag
  @SeededTest
  void fixedFasterThanFreeCircularMultiPeak(RandomSeed seed) {
    f1FasterThanf2(seed, multi, GaussianFunctionFactory.FIT_FIXED,
        GaussianFunctionFactory.FIT_FREE_CIRCLE);
  }

  @SeededTest
  void freeCircularComputesSameAsEllipticalMultiPeakNb(RandomSeed seed) {
    f1ComputesSameAsf2(seed, multi, GaussianFunctionFactory.FIT_SIMPLE_NB_FREE_CIRCLE,
        GaussianFunctionFactory.FIT_SIMPLE_NB_ELLIPTICAL);
  }

  @SpeedTag
  @SeededTest
  void freeCircularFasterThanEllipticalMultiPeakNb(RandomSeed seed) {
    f1FasterThanf2(seed, multi, GaussianFunctionFactory.FIT_SIMPLE_NB_FREE_CIRCLE,
        GaussianFunctionFactory.FIT_SIMPLE_NB_ELLIPTICAL);
  }

  @SeededTest
  void circularComputesSameAsFreeCircularMultiPeakNb(RandomSeed seed) {
    f1ComputesSameAsf2(seed, multi, GaussianFunctionFactory.FIT_SIMPLE_NB_CIRCLE,
        GaussianFunctionFactory.FIT_SIMPLE_NB_FREE_CIRCLE);
  }

  @SpeedTag
  @SeededTest
  void circularFasterThanFreeCircularMultiPeakNb(RandomSeed seed) {
    f1FasterThanf2(seed, multi, GaussianFunctionFactory.FIT_SIMPLE_NB_CIRCLE,
        GaussianFunctionFactory.FIT_SIMPLE_NB_FREE_CIRCLE);
  }

  @SeededTest
  void fixedComputesSameAsFreeCircularMultiPeakNb(RandomSeed seed) {
    f1ComputesSameAsf2(seed, multi, GaussianFunctionFactory.FIT_SIMPLE_NB_FIXED,
        GaussianFunctionFactory.FIT_SIMPLE_NB_FREE_CIRCLE);
  }

  @SpeedTag
  @SeededTest
  void fixedFasterThanFreeCircularMultiPeakNb(RandomSeed seed) {
    f1FasterThanf2(seed, multi, GaussianFunctionFactory.FIT_SIMPLE_NB_FIXED,
        GaussianFunctionFactory.FIT_SIMPLE_NB_FREE_CIRCLE);
  }

  void f1ComputesSameAsf2(RandomSeed seed, int npeaks, int flags1, int flags2) {
    final DoubleEquality eq = new DoubleEquality(1e-2, 1e-10);
    final int iter = 50;
    ArrayList<double[]> paramsList2;
    if (npeaks == 1) {
      paramsList2 = copyList(ensureDataSingle(seed, iter).paramsListSinglePeak, iter);
    } else {
      paramsList2 = copyList(ensureDataMulti(seed, iter).paramsListMultiPeak, iter);
    }

    final Gaussian2DFunction f1 =
        GaussianFunctionFactory.create2D(1, blockWidth, blockWidth, flags1, null);
    final Gaussian2DFunction f2 =
        GaussianFunctionFactory.create2D(1, blockWidth, blockWidth, flags2, null);

    final double[] dyda1 = new double[1 + npeaks * Gaussian2DFunction.PARAMETERS_PER_PEAK];
    final double[] dyda2 = new double[1 + npeaks * Gaussian2DFunction.PARAMETERS_PER_PEAK];

    final int[] gradientIndices = f1.gradientIndices();
    final int[] g1 = new int[gradientIndices.length];
    final int[] g2 = new int[gradientIndices.length];
    int nparams = 0;
    for (int i = 0; i < gradientIndices.length; i++) {
      final int index1 = f1.findGradientIndex(g1[i]);
      final int index2 = f2.findGradientIndex(g2[i]);
      if (index1 >= 0 && index2 >= 0) {
        g1[nparams] = index1;
        g2[nparams] = index2;
        nparams++;
      }
    }

    for (int i = 0; i < paramsList2.size(); i++) {
      f1.initialise(paramsList2.get(i));
      f2.initialise(paramsList2.get(i));

      for (int j = 0; j < x.length; j++) {
        final double y1 = f1.eval(x[j], dyda1);
        final double y2 = f2.eval(x[j], dyda2);

        if (!eq.almostEqualRelativeOrAbsolute(y1, y2)) {
          Assertions.fail(String.format("Not same y[%d] @ %d : %g != %g", j, i, y1, y2));
        }

        for (int ii = 0; ii < nparams; ii++) {
          if (!eq.almostEqualRelativeOrAbsolute(dyda1[g1[ii]], dyda2[g2[ii]])) {
            Assertions.fail(String.format("Not same dyda[%d] @ %d : %g != %g", j,
                gradientIndices[g1[ii]], dyda1[g1[ii]], dyda2[g2[ii]]));
          }
        }
      }
    }
  }

  void f1FasterThanf2(RandomSeed seed, int npeaks, int flags1, int flags2) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final int iter = 10000;
    ArrayList<double[]> paramsList2;
    if (npeaks == 1) {
      paramsList2 = copyList(ensureDataSingle(seed, iter).paramsListSinglePeak, iter);
    } else {
      paramsList2 = copyList(ensureDataMulti(seed, iter).paramsListMultiPeak, iter);
    }

    // Use the full list of parameters to build the functions
    final Gaussian2DFunction f1 =
        GaussianFunctionFactory.create2D(npeaks, blockWidth, blockWidth, flags1, null);
    final Gaussian2DFunction f2 =
        GaussianFunctionFactory.create2D(npeaks, blockWidth, blockWidth, flags2, null);

    final double[] dyda = new double[1 + npeaks * 6];

    for (int i = 0; i < paramsList2.size(); i++) {
      f1.initialise(paramsList2.get(i));
      for (int j = 0; j < x.length; j++) {
        f1.eval(x[j], dyda);
      }
    }

    long start1 = System.nanoTime();
    for (int i = 0; i < paramsList2.size(); i++) {
      f1.initialise(paramsList2.get(i));
      for (int j = 0; j < x.length; j++) {
        f1.eval(x[j], dyda);
      }
    }
    start1 = System.nanoTime() - start1;

    for (int i = 0; i < paramsList2.size(); i++) {
      f2.initialise(paramsList2.get(i));
      for (int j = 0; j < x.length; j++) {
        f2.eval(x[j], dyda);
      }
    }

    long start2 = System.nanoTime();
    for (int i = 0; i < paramsList2.size(); i++) {
      f2.initialise(paramsList2.get(i));
      for (int j = 0; j < x.length; j++) {
        f2.eval(x[j], dyda);
      }
    }
    start2 = System.nanoTime() - start2;

    logger.log(TestLogUtils.getTimingRecord(f1.getClass().getName(), start1,
        f2.getClass().getName(), start2));
  }

  /**
   * Create random elliptical Gaussian data an returns the data plus an estimate of the parameters.
   * Only the chosen parameters are randomised and returned for a maximum of (background, amplitude,
   * angle, xpos, ypos, xwidth, ywidth }
   *
   * @param rand the rand
   * @param npeaks the npeaks
   * @param params set on output
   * @return the data
   */
  private static double[] doubleCreateGaussianData(UniformRandomProvider rand, int npeaks,
      double[] params) {
    final int n = blockWidth * blockWidth;

    // Generate a 2D Gaussian
    final EllipticalGaussian2DFunction func =
        new EllipticalGaussian2DFunction(npeaks, blockWidth, blockWidth);
    params[0] = background + rand.nextFloat() * 5f;
    for (int i = 0, j = 0; i < npeaks; i++, j += Gaussian2DFunction.PARAMETERS_PER_PEAK) {
      params[j] = amplitude + rand.nextFloat() * 5f;
      params[j + Gaussian2DFunction.X_POSITION] = xpos + rand.nextFloat() * 2f;
      params[j + Gaussian2DFunction.Y_POSITION] = ypos + rand.nextFloat() * 2f;
      params[j + Gaussian2DFunction.X_SD] = xwidth + rand.nextFloat() * 2f;
      params[j + Gaussian2DFunction.Y_SD] = params[j + 4];
      params[j + Gaussian2DFunction.ANGLE] = 0f; // (double) (Math.PI / 4.0); // Angle
    }

    final double[] dyDa = new double[params.length];
    final double[] y = new double[n];
    func.initialise(params);
    for (int i = 0; i < y.length; i++) {
      // Add random noise
      y[i] = func.eval(i, dyDa)
          + ((rand.nextFloat() < 0.5f) ? -rand.nextFloat() * 5f : rand.nextFloat() * 5f);
    }

    // Randomise only the necessary parameters (i.e. not angle and X & Y widths should be the same)
    params[0] += ((rand.nextFloat() < 0.5f) ? -rand.nextFloat() : rand.nextFloat());
    for (int i = 0, j = 0; i < npeaks; i++, j += Gaussian2DFunction.PARAMETERS_PER_PEAK) {
      params[j + Gaussian2DFunction.X_POSITION] +=
          ((rand.nextFloat() < 0.5f) ? -rand.nextFloat() : rand.nextFloat());
      params[j + Gaussian2DFunction.Y_POSITION] +=
          ((rand.nextFloat() < 0.5f) ? -rand.nextFloat() : rand.nextFloat());
      params[j + Gaussian2DFunction.X_SD] +=
          ((rand.nextFloat() < 0.5f) ? -rand.nextFloat() : rand.nextFloat());
      params[j + Gaussian2DFunction.Y_SD] = params[j + Gaussian2DFunction.X_SD];
    }

    return y;
  }

  protected static void createGaussianData(UniformRandomProvider rand, int npeaks, int iter,
      ArrayList<double[]> paramsList, ArrayList<double[]> yList) {
    while (paramsList.size() < iter) {
      final double[] params = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK * npeaks];
      final double[] y = doubleCreateGaussianData(rand, npeaks, params);
      paramsList.add(params);
      yList.add(y);
    }
  }

  protected ArrayList<double[]> copyList(ArrayList<double[]> paramsList, int iter) {
    iter = Math.min(iter, paramsList.size());

    final ArrayList<double[]> params2List = new ArrayList<>(iter);
    for (int i = 0; i < iter; i++) {
      params2List.add(paramsList.get(i));
    }
    return params2List;
  }
}
