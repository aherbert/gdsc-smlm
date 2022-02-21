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

package uk.ac.sussex.gdsc.smlm.fitting.linear;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.ConcurrentHashMap;
import java.util.logging.Logger;
import org.apache.commons.rng.UniformRandomProvider;
import org.ejml.data.DenseMatrix64F;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.BeforeAll;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.GradientCalculator;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.SingleFreeCircularGaussian2DFunction;
import uk.ac.sussex.gdsc.test.api.TestAssertions;
import uk.ac.sussex.gdsc.test.api.TestHelper;
import uk.ac.sussex.gdsc.test.api.function.DoubleDoubleBiPredicate;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.junit5.SpeedTag;
import uk.ac.sussex.gdsc.test.rng.RngUtils;
import uk.ac.sussex.gdsc.test.utils.AssertionErrorCounter;
import uk.ac.sussex.gdsc.test.utils.RandomSeed;
import uk.ac.sussex.gdsc.test.utils.TestComplexity;
import uk.ac.sussex.gdsc.test.utils.TestLogUtils;
import uk.ac.sussex.gdsc.test.utils.TestSettings;

@SuppressWarnings({"javadoc"})
class SolverSpeedTest {

  private static Logger logger;
  private static ConcurrentHashMap<RandomSeed, Object> dataCache;

  @BeforeAll
  public static void beforeAll() {
    logger = Logger.getLogger(SolverSpeedTest.class.getName());
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

  private static class SolverSpeedTestData {
    ArrayList<float[][]> adata = new ArrayList<>();
    ArrayList<float[]> bdata = new ArrayList<>();
    final UniformRandomProvider rng;

    SolverSpeedTestData(UniformRandomProvider rng) {
      this.rng = rng;
    }
  }

  private static SolverSpeedTestData ensureData(RandomSeed seed, int size) {
    final SolverSpeedTestData data =
        (SolverSpeedTestData) dataCache.computeIfAbsent(seed, SolverSpeedTest::createData);
    final ArrayList<float[][]> adata = data.adata;
    final ArrayList<float[]> bdata = data.bdata;
    if (adata.size() < size) {
      synchronized (adata) {
        while (adata.size() < size) {
          final float[][] a = new float[6][6];
          final float[] b = new float[6];
          if (createSolverData(data.rng, a, b, false)) {
            adata.add(a);
            bdata.add(b);
          }
        }
      }
    }
    return data;
  }

  private static Object createData(RandomSeed source) {
    return new SolverSpeedTestData(RngUtils.create(source.get()));
  }

  @SeededTest
  void solveLinearAndGaussJordanReturnSameSolutionAndInversionResult(RandomSeed seed) {
    final int iter = 100;
    final SolverSpeedTestData data = ensureData(seed, iter);
    final ArrayList<double[][]> adata = copyAdouble(data.adata, iter);
    final ArrayList<double[]> bdata = copyBdouble(data.bdata, iter);
    final ArrayList<double[][]> adata2 = copyAdouble(data.adata, iter);
    final ArrayList<double[]> bdata2 = copyBdouble(data.bdata, iter);

    final GaussJordan solver = new GaussJordan();
    final EjmlLinearSolver solver2 = new EjmlLinearSolver();

    final int failureLimit = AssertionErrorCounter.computeFailureLimit(iter, 0.1);
    final AssertionErrorCounter failCounter = new AssertionErrorCounter(failureLimit, 2);

    final DoubleDoubleBiPredicate predicate = TestHelper.doublesAreClose(1e-2, 0);

    int fail = 0;
    for (int i = 0; i < adata.size(); i++) {
      final double[][] a1 = adata.get(i);
      final double[] b1 = bdata.get(i);
      final double[][] a2 = adata2.get(i);
      final double[] b2 = bdata2.get(i);
      final boolean r1 = solver.solve(a1, b1);
      final boolean r2 = solver2.solveLinear(a2, b2);
      solver2.invertLastA(a2);
      // Assertions.assertTrue("Different solve result @ " + i, r1 == r2);
      if (r1 && r2) {
        failCounter.run(0, () -> {
          TestAssertions.assertArrayTest(b1, b2, predicate, "Different b result");
        });
        failCounter.run(0, () -> {
          TestAssertions.assertArrayTest(a1, a2, predicate, "Different a result");
        });
      } else {
        fail++;
      }
    }
    if (fail > iter / 2) {
      Assertions.fail(String.format("Failed to solve %d / %d", fail, iter));
    }
  }

  @SeededTest
  void solveLinearAndGaussJordanReturnSameSolutionResult(RandomSeed seed) {
    final int iter = 100;
    final SolverSpeedTestData data = ensureData(seed, iter);
    final ArrayList<double[][]> adata = copyAdouble(data.adata, iter);
    final ArrayList<double[]> bdata = copyBdouble(data.bdata, iter);
    final ArrayList<double[][]> adata2 = copyAdouble(data.adata, iter);
    final ArrayList<double[]> bdata2 = copyBdouble(data.bdata, iter);

    final GaussJordan solver = new GaussJordan();
    final EjmlLinearSolver solver2 = new EjmlLinearSolver();

    final int failureLimit = AssertionErrorCounter.computeFailureLimit(iter, 0.1);
    final AssertionErrorCounter failCounter = new AssertionErrorCounter(failureLimit);

    final DoubleDoubleBiPredicate predicate = TestHelper.doublesAreClose(1e-2, 0);

    int fail = 0;
    for (int i = 0; i < iter; i++) {
      final double[][] a1 = adata.get(i);
      final double[] b1 = bdata.get(i);
      final double[][] a2 = adata2.get(i);
      final double[] b2 = bdata2.get(i);
      final boolean r1 = solver.solve(a1, b1);
      final boolean r2 = solver2.solve(a2, b2);
      // Assertions.assertTrue("Different solve result @ " + i, r1 == r2);
      if (r1 && r2) {
        failCounter.run(() -> {
          TestAssertions.assertArrayTest(b1, b2, predicate, "Different b result");
        });
      } else {
        fail++;
      }
    }
    if (fail > iter / 2) {
      Assertions.fail(String.format("Failed to solve %d / %d", fail, iter));
    }
  }

  @SeededTest
  void gaussJordanFloatAndDoubleReturnSameSolutionAndInversionResult(RandomSeed seed) {
    final int iter = 100;
    final SolverSpeedTestData data = ensureData(seed, iter);
    final ArrayList<float[][]> adata = copyAfloat(data.adata, iter);
    final ArrayList<float[]> bdata = copyBfloat(data.bdata, iter);
    final ArrayList<double[][]> adata2 = copyAdouble(data.adata, iter);
    final ArrayList<double[]> bdata2 = copyBdouble(data.bdata, iter);

    final GaussJordan solver = new GaussJordan();

    final int failureLimit = AssertionErrorCounter.computeFailureLimit(iter, 0.1);
    final AssertionErrorCounter failCounter = new AssertionErrorCounter(failureLimit, 2);

    final DoubleDoubleBiPredicate predicate = TestHelper.doublesAreClose(1e-2, 0);

    int fail = 0;
    for (int i = 0; i < adata.size(); i++) {
      final float[][] a1 = adata.get(i);
      final float[] b1 = bdata.get(i);
      final double[][] a2 = adata2.get(i);
      final double[] b2 = bdata2.get(i);
      final boolean r1 = solver.solve(a1, b1);
      final boolean r2 = solver.solve(a2, b2);
      // Assertions.assertTrue("Different solve result @ " + i, r1 == r2);
      if (r1 && r2) {
        final double[] db1 = SimpleArrayUtils.toDouble(b1);
        final double[][] da1 = new double[a1.length][];
        for (int j = a1.length; j-- > 0;) {
          da1[j] = SimpleArrayUtils.toDouble(a1[j]);
        }
        failCounter.run(0, () -> {
          TestAssertions.assertArrayTest(db1, b2, predicate, "Different b result");
        });
        failCounter.run(1, () -> {
          TestAssertions.assertArrayTest(da1, a2, predicate, "Different a result");
        });
      } else {
        fail++;
      }
    }
    if (fail > iter / 2) {
      Assertions.fail(String.format("Failed to solve %d / %d", fail, iter));
    }
  }

  @SpeedTag
  @SeededTest
  void solveLinearWithInversionIsNotFasterThanGaussJordanFloat(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final int iter = 10000;
    final SolverSpeedTestData data = ensureData(seed, iter);

    final GaussJordan solver = new GaussJordan();
    final EjmlLinearSolver solver2 = new EjmlLinearSolver();

    long t1 = Long.MAX_VALUE;
    long t2 = Long.MAX_VALUE;

    for (int loops = 5; loops-- > 0;) {
      final ArrayList<float[][]> adata = copyAfloat(data.adata, iter);
      final ArrayList<float[]> bdata = copyBfloat(data.bdata, iter);
      final ArrayList<double[]> adata2 = copyA2double(data.adata, iter);
      final ArrayList<double[]> bdata2 = copyBdouble(data.bdata, iter);

      final long start1 = System.nanoTime();
      runFloat(adata, bdata, iter, solver);
      t1 = Math.min(t1, System.nanoTime() - start1);

      final long start2 = System.nanoTime();
      solveLinearWithInversion(adata2, bdata2, iter, solver2);
      t2 = Math.min(t2, System.nanoTime() - start2);
    }

    logger.log(TestLogUtils.getTimingRecord("GaussJordanFloat", t1,
        "LinearSolver.solveLinearWithInversion", t2));
  }

  @SpeedTag
  @SeededTest
  void solveLinearIsFasterThanGaussJordanFloat(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final int iter = 10000;
    final SolverSpeedTestData data = ensureData(seed, iter);

    final GaussJordan solver = new GaussJordan();
    final EjmlLinearSolver solver2 = new EjmlLinearSolver();

    long t1 = Long.MAX_VALUE;
    long t2 = Long.MAX_VALUE;

    for (int loops = 5; loops-- > 0;) {
      final ArrayList<float[][]> adata = copyAfloat(data.adata, iter);
      final ArrayList<float[]> bdata = copyBfloat(data.bdata, iter);
      final ArrayList<double[]> adata2 = copyA2double(data.adata, iter);
      final ArrayList<double[]> bdata2 = copyBdouble(data.bdata, iter);

      final long start1 = System.nanoTime();
      runFloat(adata, bdata, iter, solver);
      t1 = Math.min(t1, System.nanoTime() - start1);

      final long start2 = System.nanoTime();
      solveLinear(adata2, bdata2, iter, solver2);
      t2 = Math.min(t2, System.nanoTime() - start2);
    }

    logger
        .log(TestLogUtils.getTimingRecord("GaussJordanFloat", t1, "LinearSolver.solveLinear", t2));
  }

  protected void runFloat(ArrayList<float[][]> adata, ArrayList<float[]> bdata, int iter,
      GaussJordan solver) {
    for (int i = 0; i < iter; i++) {
      solver.solve(adata.get(i), bdata.get(i));
    }
  }

  @SpeedTag
  @SeededTest
  void solveLinearWithInversionIsNotFasterThanGaussJordanDouble(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final int iter = 10000;
    final SolverSpeedTestData data = ensureData(seed, iter);

    final GaussJordan solver = new GaussJordan();
    final EjmlLinearSolver solver2 = new EjmlLinearSolver();

    long t1 = Long.MAX_VALUE;
    long t2 = Long.MAX_VALUE;

    for (int loops = 5; loops-- > 0;) {
      final ArrayList<double[][]> adata = copyAdouble(data.adata, iter);
      final ArrayList<double[]> bdata = copyBdouble(data.bdata, iter);
      final ArrayList<double[]> adata2 = copyA2double(data.adata, iter);
      final ArrayList<double[]> bdata2 = copyBdouble(data.bdata, iter);

      final long start1 = System.nanoTime();
      solveGaussJordan(adata, bdata, iter, solver);
      t1 = Math.min(t1, System.nanoTime() - start1);

      final long start2 = System.nanoTime();
      solveLinearWithInversion(adata2, bdata2, iter, solver2);
      t2 = Math.min(t2, System.nanoTime() - start2);
    }

    logger.log(TestLogUtils.getTimingRecord("GaussJordanDouble", t1,
        "LinearSolver.solveLinearWithInversion", t2));
  }

  @SpeedTag
  @SeededTest
  void solveLinearIsFasterThanGaussJordanDouble(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final int iter = 10000;
    final SolverSpeedTestData data = ensureData(seed, iter);

    final GaussJordan solver = new GaussJordan();
    final EjmlLinearSolver solver2 = new EjmlLinearSolver();

    long t1 = Long.MAX_VALUE;
    long t2 = Long.MAX_VALUE;

    for (int loops = 5; loops-- > 0;) {
      final ArrayList<double[][]> adata = copyAdouble(data.adata, iter);
      final ArrayList<double[]> bdata = copyBdouble(data.bdata, iter);
      final ArrayList<double[]> adata2 = copyA2double(data.adata, iter);
      final ArrayList<double[]> bdata2 = copyBdouble(data.bdata, iter);

      final long start1 = System.nanoTime();
      solveGaussJordan(adata, bdata, iter, solver);
      t1 = Math.min(t1, System.nanoTime() - start1);

      final long start2 = System.nanoTime();
      solveLinear(adata2, bdata2, iter, solver2);
      t2 = Math.min(t2, System.nanoTime() - start2);
    }

    logger
        .log(TestLogUtils.getTimingRecord("GaussJordanDouble", t1, "LinearSolver.solveLinear", t2));
  }

  @SpeedTag
  @SeededTest
  void solveCholeskyIsFasterThanGaussJordanDouble(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final int iter = 10000;
    final SolverSpeedTestData data = ensureData(seed, iter);

    final GaussJordan solver = new GaussJordan();
    final EjmlLinearSolver solver2 = new EjmlLinearSolver();

    long t1 = Long.MAX_VALUE;
    long t2 = Long.MAX_VALUE;

    for (int loops = 5; loops-- > 0;) {
      final ArrayList<double[][]> adata = copyAdouble(data.adata, iter);
      final ArrayList<double[]> bdata = copyBdouble(data.bdata, iter);
      final ArrayList<double[]> adata2 = copyA2double(data.adata, iter);
      final ArrayList<double[]> bdata2 = copyBdouble(data.bdata, iter);

      final long start1 = System.nanoTime();
      solveGaussJordan(adata, bdata, iter, solver);
      t1 = Math.min(t1, System.nanoTime() - start1);

      final long start2 = System.nanoTime();
      solveCholesky(adata2, bdata2, iter, solver2);
      t2 = Math.min(t2, System.nanoTime() - start2);
    }

    logger.log(
        TestLogUtils.getTimingRecord("GaussJordanDouble", t1, "LinearSolver.solveCholesky", t2));
  }

  @SpeedTag
  @SeededTest
  void solveCholeskyLdlTIsFasterThanGaussJordanDouble(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final int iter = 10000;
    final SolverSpeedTestData data = ensureData(seed, iter);

    final GaussJordan solver = new GaussJordan();
    final EjmlLinearSolver solver2 = new EjmlLinearSolver();

    long t1 = Long.MAX_VALUE;
    long t2 = Long.MAX_VALUE;

    for (int loops = 5; loops-- > 0;) {
      final ArrayList<double[][]> adata = copyAdouble(data.adata, iter);
      final ArrayList<double[]> bdata = copyBdouble(data.bdata, iter);
      final ArrayList<double[]> adata2 = copyA2double(data.adata, iter);
      final ArrayList<double[]> bdata2 = copyBdouble(data.bdata, iter);

      final long start1 = System.nanoTime();
      solveGaussJordan(adata, bdata, iter, solver);
      t1 = Math.min(t1, System.nanoTime() - start1);

      final long start2 = System.nanoTime();
      solveCholeskyLdlT(adata2, bdata2, iter, solver2);
      t2 = Math.min(t2, System.nanoTime() - start2);
    }

    logger.log(TestLogUtils.getTimingRecord("GaussJordanDouble", t1,
        "LinearSolver.solveCholeskyLDLT", t2));
  }

  @SpeedTag
  @SeededTest
  void solveIsFasterThanGaussJordanDouble(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final int iter = 10000;
    final SolverSpeedTestData data = ensureData(seed, iter);

    final GaussJordan solver = new GaussJordan();
    final EjmlLinearSolver solver2 = new EjmlLinearSolver();

    long t1 = Long.MAX_VALUE;
    long t2 = Long.MAX_VALUE;

    for (int loops = 5; loops-- > 0;) {
      final ArrayList<double[][]> adata = copyAdouble(data.adata, iter);
      final ArrayList<double[]> bdata = copyBdouble(data.bdata, iter);
      final ArrayList<double[]> adata2 = copyA2double(data.adata, iter);
      final ArrayList<double[]> bdata2 = copyBdouble(data.bdata, iter);

      final long start1 = System.nanoTime();
      solveGaussJordan(adata, bdata, iter, solver);
      t1 = Math.min(t1, System.nanoTime() - start1);

      final long start2 = System.nanoTime();
      solve(adata2, bdata2, iter, solver2);
      t2 = Math.min(t2, System.nanoTime() - start2);
    }

    logger.log(TestLogUtils.getTimingRecord("GaussJordanDouble", t1, "LinearSolver.solve", t2));
  }

  private static boolean createSolverData(UniformRandomProvider rand, float[][] alpha, float[] beta,
      boolean positiveDifinite) {
    // Generate a 2D Gaussian
    final SingleFreeCircularGaussian2DFunction func =
        new SingleFreeCircularGaussian2DFunction(10, 10);
    final double[] params = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
    params[Gaussian2DFunction.BACKGROUND] = 2 + rand.nextDouble() * 2;
    params[Gaussian2DFunction.SIGNAL] = 100 + rand.nextDouble() * 5;
    params[Gaussian2DFunction.X_POSITION] = 4.5 + rand.nextDouble();
    params[Gaussian2DFunction.Y_POSITION] = 4.5 + rand.nextDouble();
    params[Gaussian2DFunction.X_SD] = 1 + rand.nextDouble();
    params[Gaussian2DFunction.Y_SD] = 1 + rand.nextDouble();
    params[Gaussian2DFunction.ANGLE] = rand.nextDouble();

    final int[] x = new int[100];
    final double[] y = new double[100];
    func.initialise(params);
    for (int i = 0; i < x.length; i++) {
      // Add random noise
      y[i] = func.eval(i)
          + ((rand.nextDouble() < 0.5) ? -rand.nextDouble() * 5 : rand.nextDouble() * 5);
    }

    // Randomise parameters
    for (int i = 0; i < params.length; i++) {
      params[i] += (rand.nextDouble() < 0.5) ? -rand.nextDouble() : rand.nextDouble();
    }

    // Compute the Hessian and parameter gradient vector
    final GradientCalculator calc = new GradientCalculator(6);
    final double[][] alpha2 = new double[6][6];
    final double[] beta2 = new double[6];
    calc.findLinearised(y.length, y, params, alpha2, beta2, func);

    // Update the Hessian using a lambda shift
    final double lambda = 1.001;
    for (int i = 0; i < alpha2.length; i++) {
      alpha2[i][i] *= lambda;
    }

    // Copy back
    for (int i = 0; i < beta.length; i++) {
      beta[i] = (float) beta2[i];
      for (int j = 0; j < beta.length; j++) {
        alpha[i][j] = (float) alpha2[i][j];
      }
    }

    // Check for a positive definite matrix
    if (positiveDifinite) {
      final EjmlLinearSolver solver = new EjmlLinearSolver();
      return solver.solveCholeskyLdlT(copydouble(alpha), copydouble(beta));
    }

    return true;
  }

  private static ArrayList<float[][]> copyAfloat(ArrayList<float[][]> a, int iter) {
    iter = Math.min(a.size(), iter);
    final ArrayList<float[][]> a2 = new ArrayList<>(iter);
    for (int i = 0; i < iter; i++) {
      a2.add(copyfloat(a.get(i)));
    }
    return a2;
  }

  private static float[][] copyfloat(float[][] array) {
    final float[][] d2 = new float[array.length][array.length];
    for (int i = 0; i < array.length; i++) {
      for (int j = 0; j < array.length; j++) {
        d2[i][j] = array[i][j];
      }
    }
    return d2;
  }

  private static ArrayList<float[]> copyBfloat(ArrayList<float[]> bdata, int iter) {
    iter = Math.min(bdata.size(), iter);
    final ArrayList<float[]> b2 = new ArrayList<>(iter);
    for (int i = 0; i < iter; i++) {
      b2.add(Arrays.copyOf(bdata.get(i), bdata.get(i).length));
    }
    return b2;
  }

  private static ArrayList<double[][]> copyAdouble(ArrayList<float[][]> adata, int iter) {
    iter = Math.min(adata.size(), iter);
    final ArrayList<double[][]> a2 = new ArrayList<>(iter);
    for (int i = 0; i < iter; i++) {
      a2.add(copydouble(adata.get(i)));
    }
    return a2;
  }

  private static ArrayList<double[]> copyA2double(ArrayList<float[][]> adata, int iter) {
    iter = Math.min(adata.size(), iter);
    final ArrayList<double[]> a2 = new ArrayList<>(iter);
    for (int i = 0; i < iter; i++) {
      a2.add(new DenseMatrix64F(copydouble(adata.get(i))).data);
    }
    return a2;
  }

  private static double[][] copydouble(float[][] array) {
    final double[][] d2 = new double[array.length][array.length];
    for (int i = 0; i < array.length; i++) {
      for (int j = 0; j < array.length; j++) {
        d2[i][j] = array[i][j];
      }
    }
    return d2;
  }

  private static double[] copydouble(float[] array) {
    return SimpleArrayUtils.toDouble(array);
  }

  private static ArrayList<double[]> copyBdouble(ArrayList<float[]> bdata, int iter) {
    iter = Math.min(bdata.size(), iter);
    final ArrayList<double[]> b2 = new ArrayList<>(iter);
    for (int i = 0; i < iter; i++) {
      b2.add(copydouble(bdata.get(i)));
    }
    return b2;
  }

  protected void solveGaussJordan(ArrayList<double[][]> adata, ArrayList<double[]> bdata, int iter,
      GaussJordan solver) {
    iter = Math.min(iter, adata.size());
    for (int i = 0; i < iter; i++) {
      solver.solve(adata.get(i), bdata.get(i));
    }
  }

  protected void solveLinearWithInversion(ArrayList<double[]> adata, ArrayList<double[]> bdata,
      int iter, EjmlLinearSolver solver) {
    iter = Math.min(iter, adata.size());
    for (int i = 0; i < iter; i++) {
      final double[] data = adata.get(i);
      solver.solveLinear(data, bdata.get(i));
      solver.invertLastA(data);
    }
  }

  protected void solveLinear(ArrayList<double[]> adata, ArrayList<double[]> bdata, int iter,
      EjmlLinearSolver solver) {
    iter = Math.min(iter, adata.size());
    for (int i = 0; i < iter; i++) {
      solver.solveLinear(adata.get(i), bdata.get(i));
    }
  }

  protected void solveCholesky(ArrayList<double[]> adata, ArrayList<double[]> bdata, int iter,
      EjmlLinearSolver solver) {
    iter = Math.min(iter, adata.size());
    for (int i = 0; i < iter; i++) {
      solver.solveCholesky(adata.get(i), bdata.get(i));
    }
  }

  protected void solveCholeskyLdlT(ArrayList<double[]> adata, ArrayList<double[]> bdata, int iter,
      EjmlLinearSolver solver) {
    iter = Math.min(iter, adata.size());
    for (int i = 0; i < iter; i++) {
      solver.solveCholeskyLdlT(adata.get(i), bdata.get(i));
    }
  }

  protected void solve(ArrayList<double[]> adata, ArrayList<double[]> bdata, int iter,
      EjmlLinearSolver solver) {
    iter = Math.min(iter, adata.size());
    for (int i = 0; i < iter; i++) {
      solver.solve(adata.get(i), bdata.get(i));
    }
  }
}
