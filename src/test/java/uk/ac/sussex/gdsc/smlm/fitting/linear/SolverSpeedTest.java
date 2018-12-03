package uk.ac.sussex.gdsc.smlm.fitting.linear;

import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.GradientCalculator;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.SingleFreeCircularGaussian2DFunction;
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

import org.apache.commons.math3.util.FastMath;
import org.apache.commons.rng.UniformRandomProvider;
import org.ejml.data.DenseMatrix64F;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.BeforeAll;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.ConcurrentHashMap;
import java.util.function.Function;
import java.util.logging.Logger;

@SuppressWarnings({"javadoc"})
public class SolverSpeedTest implements Function<RandomSeed, Object> {
  private static Logger logger;
  private static ConcurrentHashMap<RandomSeed, Object> ConcurrentHashMap;

  @BeforeAll
  public static void beforeAll() {
    logger = Logger.getLogger(SolverSpeedTest.class.getName());
    ConcurrentHashMap = new ConcurrentHashMap<>();
  }

  @AfterAll
  public static void afterAll() {
    ConcurrentHashMap.clear();
    ConcurrentHashMap = null;
    logger = null;
  }

  private class SolverSpeedTestData {
    ArrayList<float[][]> Adata = new ArrayList<>();
    ArrayList<float[]> Bdata = new ArrayList<>();
    final UniformRandomProvider r;

    SolverSpeedTestData(UniformRandomProvider r) {
      this.r = r;
    }
  }

  private SolverSpeedTestData ensureData(RandomSeed seed, int size) {
    final SolverSpeedTestData data =
        (SolverSpeedTestData) ConcurrentHashMap.computeIfAbsent(seed, this);
    final ArrayList<float[][]> Adata = data.Adata;
    final ArrayList<float[]> Bdata = data.Bdata;
    if (Adata.size() < size) {
      synchronized (Adata) {
        while (Adata.size() < size) {
          final float[][] a = new float[6][6];
          final float[] b = new float[6];
          if (createData(data.r, a, b, false)) {
            Adata.add(a);
            Bdata.add(b);
          }
        }
      }
    }
    return data;
  }

  @Override
  public Object apply(RandomSeed source) {
    return new SolverSpeedTestData(RngUtils.create(source.getSeed()));
  }

  @SeededTest
  public void solveLinearAndGaussJordanReturnSameSolutionAndInversionResult(RandomSeed seed) {
    final int ITER = 100;
    final SolverSpeedTestData data = ensureData(seed, ITER);
    final ArrayList<double[][]> A = copyAdouble(data.Adata, ITER);
    final ArrayList<double[]> B = copyBdouble(data.Bdata, ITER);
    final ArrayList<double[][]> A2 = copyAdouble(data.Adata, ITER);
    final ArrayList<double[]> B2 = copyBdouble(data.Bdata, ITER);

    final GaussJordan solver = new GaussJordan();
    final EJMLLinearSolver solver2 = new EJMLLinearSolver();

    final int failureLimit = TestCounter.computeFailureLimit(ITER, 0.1);
    final TestCounter failCounter = new TestCounter(failureLimit, 2);

    final DoubleDoubleBiPredicate predicate = TestHelper.doublesAreClose(1e-2, 0);

    int c = 0;
    for (int i = 0; i < A.size(); i++) {
      final double[][] a = A.get(i);
      final double[] b = B.get(i);
      final double[][] a2 = A2.get(i);
      final double[] b2 = B2.get(i);
      final boolean r1 = solver.solve(a, b);
      final boolean r2 = solver2.solveLinear(a2, b2);
      solver2.invertLastA(a2);
      // Assertions.assertTrue("Different solve result @ " + i, r1 == r2);
      if (r1 && r2) {
        failCounter.run(0, () -> {
          TestAssertions.assertArrayTest(b, b2, predicate, "Different b result");
        });
        failCounter.run(0, () -> {
          TestAssertions.assertArrayTest(a, a2, predicate, "Different a result");
        });
      } else {
        c++;
      }
    }
    if (c > ITER / 2) {
      Assertions.fail(String.format("Failed to solve %d / %d", c, ITER));
    }
  }

  @SeededTest
  public void solveLinearAndGaussJordanReturnSameSolutionResult(RandomSeed seed) {
    final int ITER = 100;
    final SolverSpeedTestData data = ensureData(seed, ITER);
    final ArrayList<double[][]> A = copyAdouble(data.Adata, ITER);
    final ArrayList<double[]> B = copyBdouble(data.Bdata, ITER);
    final ArrayList<double[][]> A2 = copyAdouble(data.Adata, ITER);
    final ArrayList<double[]> B2 = copyBdouble(data.Bdata, ITER);

    final GaussJordan solver = new GaussJordan();
    final EJMLLinearSolver solver2 = new EJMLLinearSolver();

    final int failureLimit = TestCounter.computeFailureLimit(ITER, 0.1);
    final TestCounter failCounter = new TestCounter(failureLimit);

    final DoubleDoubleBiPredicate predicate = TestHelper.doublesAreClose(1e-2, 0);

    int c = 0;
    for (int i = 0; i < ITER; i++) {
      final double[][] a = A.get(i);
      final double[] b = B.get(i);
      final double[][] a2 = A2.get(i);
      final double[] b2 = B2.get(i);
      final boolean r1 = solver.solve(a, b);
      final boolean r2 = solver2.solve(a2, b2);
      // Assertions.assertTrue("Different solve result @ " + i, r1 == r2);
      if (r1 && r2) {
        failCounter.run(() -> {
          TestAssertions.assertArrayTest(b, b2, predicate, "Different b result");
        });
      } else {
        c++;
      }
    }
    if (c > ITER / 2) {
      Assertions.fail(String.format("Failed to solve %d / %d", c, ITER));
    }
  }

  @SeededTest
  public void gaussJordanFloatAndDoubleReturnSameSolutionAndInversionResult(RandomSeed seed) {
    final int ITER = 100;
    final SolverSpeedTestData data = ensureData(seed, ITER);
    final ArrayList<float[][]> A = copyAfloat(data.Adata, ITER);
    final ArrayList<float[]> B = copyBfloat(data.Bdata, ITER);
    final ArrayList<double[][]> A2 = copyAdouble(data.Adata, ITER);
    final ArrayList<double[]> B2 = copyBdouble(data.Bdata, ITER);

    final GaussJordan solver = new GaussJordan();

    final int failureLimit = TestCounter.computeFailureLimit(ITER, 0.1);
    final TestCounter failCounter = new TestCounter(failureLimit, 2);

    final DoubleDoubleBiPredicate predicate = TestHelper.doublesAreClose(1e-2, 0);

    int c = 0;
    for (int i = 0; i < A.size(); i++) {
      final float[][] a = A.get(i);
      final float[] b = B.get(i);
      final double[][] a2 = A2.get(i);
      final double[] b2 = B2.get(i);
      final boolean r1 = solver.solve(a, b);
      final boolean r2 = solver.solve(a2, b2);
      // Assertions.assertTrue("Different solve result @ " + i, r1 == r2);
      if (r1 && r2) {
        final double[] b1 = SimpleArrayUtils.toDouble(b);
        final double[][] a1 = new double[a.length][];
        for (int j = a1.length; j-- > 0;) {
          a1[j] = SimpleArrayUtils.toDouble(a[j]);
        }
        failCounter.run(0, () -> {
          TestAssertions.assertArrayTest(b1, b2, predicate, "Different b result");
        });
        failCounter.run(1, () -> {
          TestAssertions.assertArrayTest(a1, a2, predicate, "Different a result");
        });
      } else {
        c++;
      }
    }
    if (c > ITER / 2) {
      Assertions.fail(String.format("Failed to solve %d / %d", c, ITER));
    }
  }

  @SpeedTag
  @SeededTest
  public void solveLinearWithInversionIsNotFasterThanGaussJordanFloat(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final int ITER = 10000;
    final SolverSpeedTestData data = ensureData(seed, ITER);

    final GaussJordan solver = new GaussJordan();
    final EJMLLinearSolver solver2 = new EJMLLinearSolver();

    long t1 = Long.MAX_VALUE, t2 = Long.MAX_VALUE;

    for (int loops = 5; loops-- > 0;) {
      final ArrayList<float[][]> A = copyAfloat(data.Adata, ITER);
      final ArrayList<float[]> B = copyBfloat(data.Bdata, ITER);
      final ArrayList<double[]> A2 = copyA2double(data.Adata, ITER);
      final ArrayList<double[]> B2 = copyBdouble(data.Bdata, ITER);

      final long start1 = System.nanoTime();
      runFloat(A, B, ITER, solver);
      t1 = Math.min(t1, System.nanoTime() - start1);

      final long start2 = System.nanoTime();
      solveLinearWithInversion(A2, B2, ITER, solver2);
      t2 = Math.min(t2, System.nanoTime() - start2);
    }

    logger.log(TestLogUtils.getTimingRecord("GaussJordanFloat", t1,
        "LinearSolver.solveLinearWithInversion", t2));
  }

  @SpeedTag
  @SeededTest
  public void solveLinearIsFasterThanGaussJordanFloat(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final int ITER = 10000;
    final SolverSpeedTestData data = ensureData(seed, ITER);

    final GaussJordan solver = new GaussJordan();
    final EJMLLinearSolver solver2 = new EJMLLinearSolver();

    long t1 = Long.MAX_VALUE, t2 = Long.MAX_VALUE;

    for (int loops = 5; loops-- > 0;) {
      final ArrayList<float[][]> A = copyAfloat(data.Adata, ITER);
      final ArrayList<float[]> B = copyBfloat(data.Bdata, ITER);
      final ArrayList<double[]> A2 = copyA2double(data.Adata, ITER);
      final ArrayList<double[]> B2 = copyBdouble(data.Bdata, ITER);

      final long start1 = System.nanoTime();
      runFloat(A, B, ITER, solver);
      t1 = Math.min(t1, System.nanoTime() - start1);

      final long start2 = System.nanoTime();
      solveLinear(A2, B2, ITER, solver2);
      t2 = Math.min(t2, System.nanoTime() - start2);
    }

    logger
        .log(TestLogUtils.getTimingRecord("GaussJordanFloat", t1, "LinearSolver.solveLinear", t2));
  }

  protected void runFloat(ArrayList<float[][]> A, ArrayList<float[]> B, int ITER,
      GaussJordan solver) {
    for (int i = 0; i < ITER; i++) {
      solver.solve(A.get(i), B.get(i));
    }
  }

  @SpeedTag
  @SeededTest
  public void solveLinearWithInversionIsNotFasterThanGaussJordanDouble(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final int ITER = 10000;
    final SolverSpeedTestData data = ensureData(seed, ITER);

    final GaussJordan solver = new GaussJordan();
    final EJMLLinearSolver solver2 = new EJMLLinearSolver();

    long t1 = Long.MAX_VALUE, t2 = Long.MAX_VALUE;

    for (int loops = 5; loops-- > 0;) {
      final ArrayList<double[][]> A = copyAdouble(data.Adata, ITER);
      final ArrayList<double[]> B = copyBdouble(data.Bdata, ITER);
      final ArrayList<double[]> A2 = copyA2double(data.Adata, ITER);
      final ArrayList<double[]> B2 = copyBdouble(data.Bdata, ITER);

      final long start1 = System.nanoTime();
      solveGaussJordan(A, B, ITER, solver);
      t1 = Math.min(t1, System.nanoTime() - start1);

      final long start2 = System.nanoTime();
      solveLinearWithInversion(A2, B2, ITER, solver2);
      t2 = Math.min(t2, System.nanoTime() - start2);
    }

    logger.log(TestLogUtils.getTimingRecord("GaussJordanDouble", t1,
        "LinearSolver.solveLinearWithInversion", t2));
  }

  @SpeedTag
  @SeededTest
  public void solveLinearIsFasterThanGaussJordanDouble(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final int ITER = 10000;
    final SolverSpeedTestData data = ensureData(seed, ITER);

    final GaussJordan solver = new GaussJordan();
    final EJMLLinearSolver solver2 = new EJMLLinearSolver();

    long t1 = Long.MAX_VALUE, t2 = Long.MAX_VALUE;

    for (int loops = 5; loops-- > 0;) {
      final ArrayList<double[][]> A = copyAdouble(data.Adata, ITER);
      final ArrayList<double[]> B = copyBdouble(data.Bdata, ITER);
      final ArrayList<double[]> A2 = copyA2double(data.Adata, ITER);
      final ArrayList<double[]> B2 = copyBdouble(data.Bdata, ITER);

      final long start1 = System.nanoTime();
      solveGaussJordan(A, B, ITER, solver);
      t1 = Math.min(t1, System.nanoTime() - start1);

      final long start2 = System.nanoTime();
      solveLinear(A2, B2, ITER, solver2);
      t2 = Math.min(t2, System.nanoTime() - start2);
    }

    logger
        .log(TestLogUtils.getTimingRecord("GaussJordanDouble", t1, "LinearSolver.solveLinear", t2));
  }

  @SpeedTag
  @SeededTest
  public void solveCholeskyIsFasterThanGaussJordanDouble(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final int ITER = 10000;
    final SolverSpeedTestData data = ensureData(seed, ITER);

    final GaussJordan solver = new GaussJordan();
    final EJMLLinearSolver solver2 = new EJMLLinearSolver();

    long t1 = Long.MAX_VALUE, t2 = Long.MAX_VALUE;

    for (int loops = 5; loops-- > 0;) {
      final ArrayList<double[][]> A = copyAdouble(data.Adata, ITER);
      final ArrayList<double[]> B = copyBdouble(data.Bdata, ITER);
      final ArrayList<double[]> A2 = copyA2double(data.Adata, ITER);
      final ArrayList<double[]> B2 = copyBdouble(data.Bdata, ITER);

      final long start1 = System.nanoTime();
      solveGaussJordan(A, B, ITER, solver);
      t1 = Math.min(t1, System.nanoTime() - start1);

      final long start2 = System.nanoTime();
      solveCholesky(A2, B2, ITER, solver2);
      t2 = Math.min(t2, System.nanoTime() - start2);
    }

    logger.log(
        TestLogUtils.getTimingRecord("GaussJordanDouble", t1, "LinearSolver.solveCholesky", t2));
  }

  @SpeedTag
  @SeededTest
  public void solveCholeskyLDLTIsFasterThanGaussJordanDouble(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final int ITER = 10000;
    final SolverSpeedTestData data = ensureData(seed, ITER);

    final GaussJordan solver = new GaussJordan();
    final EJMLLinearSolver solver2 = new EJMLLinearSolver();

    long t1 = Long.MAX_VALUE, t2 = Long.MAX_VALUE;

    for (int loops = 5; loops-- > 0;) {
      final ArrayList<double[][]> A = copyAdouble(data.Adata, ITER);
      final ArrayList<double[]> B = copyBdouble(data.Bdata, ITER);
      final ArrayList<double[]> A2 = copyA2double(data.Adata, ITER);
      final ArrayList<double[]> B2 = copyBdouble(data.Bdata, ITER);

      final long start1 = System.nanoTime();
      solveGaussJordan(A, B, ITER, solver);
      t1 = Math.min(t1, System.nanoTime() - start1);

      final long start2 = System.nanoTime();
      solveCholeskyLDLT(A2, B2, ITER, solver2);
      t2 = Math.min(t2, System.nanoTime() - start2);
    }

    logger.log(TestLogUtils.getTimingRecord("GaussJordanDouble", t1,
        "LinearSolver.solveCholeskyLDLT", t2));
  }

  @SpeedTag
  @SeededTest
  public void solveIsFasterThanGaussJordanDouble(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final int ITER = 10000;
    final SolverSpeedTestData data = ensureData(seed, ITER);

    final GaussJordan solver = new GaussJordan();
    final EJMLLinearSolver solver2 = new EJMLLinearSolver();

    long t1 = Long.MAX_VALUE, t2 = Long.MAX_VALUE;

    for (int loops = 5; loops-- > 0;) {
      final ArrayList<double[][]> A = copyAdouble(data.Adata, ITER);
      final ArrayList<double[]> B = copyBdouble(data.Bdata, ITER);
      final ArrayList<double[]> A2 = copyA2double(data.Adata, ITER);
      final ArrayList<double[]> B2 = copyBdouble(data.Bdata, ITER);

      final long start1 = System.nanoTime();
      solveGaussJordan(A, B, ITER, solver);
      t1 = Math.min(t1, System.nanoTime() - start1);

      final long start2 = System.nanoTime();
      solve(A2, B2, ITER, solver2);
      t2 = Math.min(t2, System.nanoTime() - start2);
    }

    logger.log(TestLogUtils.getTimingRecord("GaussJordanDouble", t1, "LinearSolver.solve", t2));
  }

  private static boolean createData(UniformRandomProvider rand, float[][] alpha, float[] beta,
      boolean positiveDifinite) {
    // Generate a 2D Gaussian
    final SingleFreeCircularGaussian2DFunction func =
        new SingleFreeCircularGaussian2DFunction(10, 10);
    final double[] a = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
    a[Gaussian2DFunction.BACKGROUND] = 2 + rand.nextDouble() * 2;
    a[Gaussian2DFunction.SIGNAL] = 100 + rand.nextDouble() * 5;
    a[Gaussian2DFunction.X_POSITION] = 4.5 + rand.nextDouble();
    a[Gaussian2DFunction.Y_POSITION] = 4.5 + rand.nextDouble();
    a[Gaussian2DFunction.X_SD] = 1 + rand.nextDouble();
    a[Gaussian2DFunction.Y_SD] = 1 + rand.nextDouble();
    a[Gaussian2DFunction.ANGLE] = rand.nextDouble();

    final int[] x = new int[100];
    final double[] y = new double[100];
    func.initialise(a);
    for (int i = 0; i < x.length; i++) {
      // Add random noise
      y[i] = func.eval(i)
          + ((rand.nextDouble() < 0.5) ? -rand.nextDouble() * 5 : rand.nextDouble() * 5);
    }

    // Randomise parameters
    for (int i = 0; i < a.length; i++) {
      a[i] += (rand.nextDouble() < 0.5) ? -rand.nextDouble() : rand.nextDouble();
    }

    // Compute the Hessian and parameter gradient vector
    final GradientCalculator calc = new GradientCalculator(6);
    final double[][] alpha2 = new double[6][6];
    final double[] beta2 = new double[6];
    calc.findLinearised(y.length, y, a, alpha2, beta2, func);

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
      final EJMLLinearSolver solver = new EJMLLinearSolver();
      return solver.solveCholeskyLDLT(copydouble(alpha), copydouble(beta));
    }

    return true;
  }

  private static ArrayList<float[][]> copyAfloat(ArrayList<float[][]> a, int iter) {
    iter = FastMath.min(a.size(), iter);
    final ArrayList<float[][]> a2 = new ArrayList<>(iter);
    for (int i = 0; i < iter; i++) {
      a2.add(copyfloat(a.get(i)));
    }
    return a2;
  }

  private static float[][] copyfloat(float[][] d) {
    final float[][] d2 = new float[d.length][d.length];
    for (int i = 0; i < d.length; i++) {
      for (int j = 0; j < d.length; j++) {
        d2[i][j] = d[i][j];
      }
    }
    return d2;
  }

  private static ArrayList<float[]> copyBfloat(ArrayList<float[]> b, int iter) {
    iter = FastMath.min(b.size(), iter);
    final ArrayList<float[]> b2 = new ArrayList<>(iter);
    for (int i = 0; i < iter; i++) {
      b2.add(Arrays.copyOf(b.get(i), b.get(i).length));
    }
    return b2;
  }

  private static ArrayList<double[][]> copyAdouble(ArrayList<float[][]> a, int iter) {
    iter = FastMath.min(a.size(), iter);
    final ArrayList<double[][]> a2 = new ArrayList<>(iter);
    for (int i = 0; i < iter; i++) {
      a2.add(copydouble(a.get(i)));
    }
    return a2;
  }

  private static ArrayList<double[]> copyA2double(ArrayList<float[][]> a, int iter) {
    iter = FastMath.min(a.size(), iter);
    final ArrayList<double[]> a2 = new ArrayList<>(iter);
    for (int i = 0; i < iter; i++) {
      a2.add(new DenseMatrix64F(copydouble(a.get(i))).data);
    }
    return a2;
  }

  private static double[][] copydouble(float[][] d) {
    final double[][] d2 = new double[d.length][d.length];
    for (int i = 0; i < d.length; i++) {
      for (int j = 0; j < d.length; j++) {
        d2[i][j] = d[i][j];
      }
    }
    return d2;
  }

  private static ArrayList<double[]> copyBdouble(ArrayList<float[]> b, int iter) {
    iter = FastMath.min(b.size(), iter);
    final ArrayList<double[]> b2 = new ArrayList<>(iter);
    for (int i = 0; i < iter; i++) {
      b2.add(copydouble(b.get(i)));
    }
    return b2;
  }

  private static double[] copydouble(float[] d) {
    final double[] d2 = new double[d.length];
    for (int i = 0; i < d.length; i++) {
      d2[i] = d[i];
    }
    return d2;
  }

  protected void solveGaussJordan(ArrayList<double[][]> A, ArrayList<double[]> B, int ITER,
      GaussJordan solver) {
    ITER = FastMath.min(ITER, A.size());
    for (int i = 0; i < ITER; i++) {
      solver.solve(A.get(i), B.get(i));
    }
  }

  protected void solveLinearWithInversion(ArrayList<double[]> A, ArrayList<double[]> B, int ITER,
      EJMLLinearSolver solver) {
    ITER = FastMath.min(ITER, A.size());
    for (int i = 0; i < ITER; i++) {
      final double[] a = A.get(i);
      solver.solveLinear(a, B.get(i));
      solver.invertLastA(a);
    }
  }

  protected void solveLinear(ArrayList<double[]> A, ArrayList<double[]> B, int ITER,
      EJMLLinearSolver solver) {
    ITER = FastMath.min(ITER, A.size());
    for (int i = 0; i < ITER; i++) {
      solver.solveLinear(A.get(i), B.get(i));
    }
  }

  protected void solveCholesky(ArrayList<double[]> A, ArrayList<double[]> B, int ITER,
      EJMLLinearSolver solver) {
    ITER = FastMath.min(ITER, A.size());
    for (int i = 0; i < ITER; i++) {
      solver.solveCholesky(A.get(i), B.get(i));
    }
  }

  protected void solveCholeskyLDLT(ArrayList<double[]> A, ArrayList<double[]> B, int ITER,
      EJMLLinearSolver solver) {
    ITER = FastMath.min(ITER, A.size());
    for (int i = 0; i < ITER; i++) {
      solver.solveCholeskyLDLT(A.get(i), B.get(i));
    }
  }

  protected void solve(ArrayList<double[]> A, ArrayList<double[]> B, int ITER,
      EJMLLinearSolver solver) {
    ITER = FastMath.min(ITER, A.size());
    for (int i = 0; i < ITER; i++) {
      solver.solve(A.get(i), B.get(i));
    }
  }
}
