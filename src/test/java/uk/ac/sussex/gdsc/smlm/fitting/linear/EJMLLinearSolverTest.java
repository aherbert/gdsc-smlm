/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2019 Alex Herbert
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

import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.RandomGeneratorAdapter;
import uk.ac.sussex.gdsc.core.utils.TurboList;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.GradientCalculator;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.GradientCalculatorUtils;
import uk.ac.sussex.gdsc.smlm.function.ValueProcedure;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import uk.ac.sussex.gdsc.smlm.math3.distribution.CustomPoissonDistribution;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.junit5.SpeedTag;
import uk.ac.sussex.gdsc.test.rng.RngUtils;
import uk.ac.sussex.gdsc.test.utils.BaseTimingTask;
import uk.ac.sussex.gdsc.test.utils.TestComplexity;
import uk.ac.sussex.gdsc.test.utils.TestLogUtils;
import uk.ac.sussex.gdsc.test.utils.TestSettings;
import uk.ac.sussex.gdsc.test.utils.TimingService;
import uk.ac.sussex.gdsc.test.utils.functions.FunctionUtils;

import org.ejml.data.DenseMatrix64F;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;

@SuppressWarnings({"javadoc"})
public class EJMLLinearSolverTest {
  private static Logger logger;

  @BeforeAll
  public static void beforeAll() {
    logger = Logger.getLogger(EJMLLinearSolverTest.class.getName());
  }

  @AfterAll
  public static void afterAll() {
    logger = null;
  }

  //@formatter:off
  @Test
  public void canSolveLinearEquation()
  {
    final EjmlLinearSolver solver = new EjmlLinearSolver(5e-3, 1e-6);

    // Solves (one) linear equation, a x = b, for x[n]

    // Taken from https://en.wikipedia.org/wiki/Positive-definite_matrix
    final double[][] a = new double[][] {
      new double[] { 2, -1, 0 },
      new double[] { -1, 2, -1 },
      new double[] { 0, -1, 2 } };
    final double[] b = new double[] { 3, 3, 4 };

    // Expected solution
    final double[] x = new double[] { 4.75, 6.5, 5.25 };
    final double[][] a_inv = new double[][] {
      new double[] { 0.75, 0.5, 0.25 },
      new double[] { 0.5, 1, 0.5 },
      new double[] { 0.25, 0.5, 0.75 } };

    final boolean result = solver.solve(a, b);
    solver.invertLastA(a);

    Assertions.assertTrue(result, "Failed to invert");
    Assertions.assertArrayEquals(x, b, 1e-4f, "Bad solution");

    if (logger.isLoggable(Level.INFO)) {
      logger.info(FunctionUtils.getSupplier("x = %s", Arrays.toString(b)));
    }
    for (int i = 0; i < b.length; i++)
    {
      if (logger.isLoggable(Level.INFO)) {
        logger.info(FunctionUtils.getSupplier("a[%d] = %s", i, Arrays.toString(a[i])));
      }
      Assertions.assertArrayEquals(a_inv[i], a[i], 1e-4f, "Bad inversion");
    }
  }

  @Test
  public void canSolveLinearEquationWithZeroInB()
  {
    final EjmlLinearSolver solver = new EjmlLinearSolver(5e-3, 1e-6);

    // Solves (one) linear equation, a x = b, for x[n]

    // Taken from https://en.wikipedia.org/wiki/Positive-definite_matrix
    final double[][] a = new double[][] {
      new double[] { 2, -1, 0 },
      new double[] { -1, 2, -1 },
      new double[] { 0, -1, 2 } };
    final double[] b = new double[] { 3, 0, 4 };

    // Expected solution
    final double[] x = new double[] { 3.25, 3.5, 3.75 };
    final double[][] a_inv = new double[][] {
      new double[] { 0.75, 0.5, 0.25 },
      new double[] { 0.5, 1, 0.5 },
      new double[] { 0.25, 0.5, 0.75 } };

    final boolean result = solver.solve(a, b);
    solver.invertLastA(a);

    Assertions.assertTrue(result, "Failed to invert");
    Assertions.assertArrayEquals(x, b, 1e-4f, "Bad solution");

    if (logger.isLoggable(Level.INFO)) {
      logger.info(FunctionUtils.getSupplier("x = %s", Arrays.toString(b)));
    }
    for (int i = 0; i < b.length; i++)
    {
      if (logger.isLoggable(Level.INFO)) {
        logger.info(FunctionUtils.getSupplier("a[%d] = %s", i, Arrays.toString(a[i])));
      }
      Assertions.assertArrayEquals(a_inv[i], a[i], 1e-4f, "Bad inversion");
    }
  }

  @Test
  public void canSolveLinearEquationWithZeroInA()
  {
    final EjmlLinearSolver solver = new EjmlLinearSolver(5e-3, 1e-6);

    // Solves (one) linear equation, a x = b, for x[n]

    final double[][] a = new double[][] {
      new double[] { 2, 0, -1, 0 },
      new double[] { 0, 0, 0, 0 },
      new double[] { -1, 0, 2, -1 },
      new double[] { 0, 0, -1, 2 } };
    final double[] b = new double[] { 3, 0, 3, 4 };

    // Expected solution
    final double[] x = new double[] { 4.75, 0, 6.5, 5.25 };
    final double[][] a_inv = new double[][] {
      new double[] { 0.75, 0, 0.5, 0.25 },
      new double[] { 0, 0, 0, 0 },
      new double[] { 0.5, 0, 1, 0.5 },
      new double[] { 0.25, 0, 0.5, 0.75 } };

    final boolean result = solver.solve(a, b);
    solver.invertLastA(a);

    Assertions.assertTrue(result, "Failed to invert");
    Assertions.assertArrayEquals(x, b, 1e-4f, "Bad solution");

    if (logger.isLoggable(Level.INFO)) {
      logger.info(FunctionUtils.getSupplier("x = %s", Arrays.toString(b)));
    }
    for (int i = 0; i < b.length; i++)
    {
      if (logger.isLoggable(Level.INFO)) {
        logger.info(FunctionUtils.getSupplier("a[%d] = %s", i, Arrays.toString(a[i])));
      }
      Assertions.assertArrayEquals(a_inv[i], a[i], 1e-4f, "Bad inversion");
    }
  }

  @Test
  public void canSolveLinearEquationWithZerosInA()
  {
    final EjmlLinearSolver solver = new EjmlLinearSolver();
    final DoubleEquality eq = new DoubleEquality(5e-3, 1e-16);
    solver.setEqual(eq);

    // Solves (one) linear equation, a x = b, for x[n]

    final double[][] a = new double[][] {
      new double[] { 2, 0, -1, 0, 0, 0 },
      new double[] { 0, 0, 0, 0, 0, 0 },
      new double[] { -1, 0, 2, 0, 0, -1 },
      new double[] { 0, 0, 0, 0, 0, 0 },
      new double[] { 0, 0, 0, 0, 0, 0 },
      new double[] { 0, 0, -1, 0, 0, 2 } };
    final double[] b = new double[] { 3, 0, 3, 0, 0, 4 };

    // Expected solution
    final double[] x = new double[] { 4.75, 0, 6.5, 0, 0, 5.25 };
    final double[][] a_inv = new double[][] {
      new double[] { 0.75, 0, 0.5, 0, 0, 0.25 },
      new double[] { 0, 0, 0, 0, 0, 0 },
      new double[] { 0.5, 0, 1, 0, 0, 0.5 },
      new double[] { 0, 0, 0, 0, 0, 0 },
      new double[] { 0, 0, 0, 0, 0, 0 },
      new double[] { 0.25, 0, 0.5, 0, 0, 0.75 } };

    final boolean result = solver.solve(a, b);
    solver.invertLastA(a);

    Assertions.assertTrue(result, "Failed to invert");
    Assertions.assertArrayEquals(x, b, 1e-4f, "Bad solution");

    if (logger.isLoggable(Level.INFO)) {
      logger.info(FunctionUtils.getSupplier("x = %s", Arrays.toString(b)));
    }
    for (int i = 0; i < b.length; i++)
    {
      if (logger.isLoggable(Level.INFO)) {
        logger.info(FunctionUtils.getSupplier("a[%d] = %s", i, Arrays.toString(a[i])));
      }
      Assertions.assertArrayEquals(a_inv[i], a[i], 1e-4f, "Bad inversion");
    }
  }

  @Test
  public void canInvert()
  {
    final EjmlLinearSolver solver = EjmlLinearSolver.createForInversion(1e-2);

    // Solves (one) linear equation, a x = b, for x[n]

    // Taken from https://en.wikipedia.org/wiki/Positive-definite_matrix
    final double[][] a = new double[][] {
      new double[] { 2, -1, 0 },
      new double[] { -1, 2, -1 },
      new double[] { 0, -1, 2 } };

    // Expected solution
    final double[][] a_inv = new double[][] {
      new double[] { 0.75, 0.5, 0.25 },
      new double[] { 0.5, 1, 0.5 },
      new double[] { 0.25, 0.5, 0.75 } };

    final boolean result = solver.invert(a);

    Assertions.assertTrue(result, "Failed to invert");

    for (int i = 0; i < a[0].length; i++)
    {
      if (logger.isLoggable(Level.INFO)) {
        logger.info(FunctionUtils.getSupplier("a[%d] = %s", i, Arrays.toString(a[i])));
      }
      Assertions.assertArrayEquals(a_inv[i], a[i], 1e-4f, "Bad inversion");
    }
  }

  @Test
  public void canInvertWithZeros()
  {
    final EjmlLinearSolver solver = EjmlLinearSolver.createForInversion(1e-2);

    // Solves (one) linear equation, a x = b, for x[n]

    // Taken from https://en.wikipedia.org/wiki/Positive-definite_matrix
    final double[][] a = new double[][] {
      new double[] { 2, 0, -1, 0, 0, 0 },
      new double[] { 0, 0, 0, 0, 0, 0 },
      new double[] { -1, 0, 2, 0, 0, -1 },
      new double[] { 0, 0, 0, 0, 0, 0 },
      new double[] { 0, 0, 0, 0, 0, 0 },
      new double[] { 0, 0, -1, 0, 0, 2 } };

    // Expected solution
    final double[][] a_inv = new double[][] {
      new double[] { 0.75, 0, 0.5, 0, 0, 0.25 },
      new double[] { 0, 0, 0, 0, 0, 0 },
      new double[] { 0.5, 0, 1, 0, 0, 0.5 },
      new double[] { 0, 0, 0, 0, 0, 0 },
      new double[] { 0, 0, 0, 0, 0, 0 },
      new double[] { 0.25, 0, 0.5, 0, 0, 0.75 } };

    final boolean result = solver.invert(a);

    Assertions.assertTrue(result, "Failed to invert");

    for (int i = 0; i < a[0].length; i++)
    {
      if (logger.isLoggable(Level.INFO)) {
        logger.info(FunctionUtils.getSupplier("a[%d] = %s", i, Arrays.toString(a[i])));
      }
      Assertions.assertArrayEquals(a_inv[i], a[i], 1e-4f, "Bad inversion");
    }
  }

  @Test
  public void canInvertDiagonal()
  {
    final EjmlLinearSolver solver = EjmlLinearSolver.createForInversion(1e-2);

    // Solves (one) linear equation, a x = b, for x[n]

    // Taken from https://en.wikipedia.org/wiki/Positive-definite_matrix
    final double[][] a = new double[][] {
      new double[] { 2, -1, 0 },
      new double[] { -1, 2, -1 },
      new double[] { 0, -1, 2 } };

    // Expected solution
    final double[] e = new double[] { 0.75, 1, 0.75 };

    final double[] o = solver.invertDiagonal(a);

    Assertions.assertNotNull(o, "Failed to invert");

    if (logger.isLoggable(Level.INFO)) {
      logger.info(FunctionUtils.getSupplier("a diagonal = %s", Arrays.toString(o)));
    }
    Assertions.assertArrayEquals(e, o, 1e-4, "Bad inversion");
  }

  @Test
  public void canInvertDiagonalWithZeros()
  {
    final EjmlLinearSolver solver = EjmlLinearSolver.createForInversion(1e-2);

    // Solves (one) linear equation, a x = b, for x[n]

    // Taken from https://en.wikipedia.org/wiki/Positive-definite_matrix
    final double[][] a = new double[][] {
      new double[] { 2, 0, -1, 0, 0, 0 },
      new double[] { 0, 0, 0, 0, 0, 0 },
      new double[] { -1, 0, 2, 0, 0, -1 },
      new double[] { 0, 0, 0, 0, 0, 0 },
      new double[] { 0, 0, 0, 0, 0, 0 },
      new double[] { 0, 0, -1, 0, 0, 2 } };

    // Expected solution
    final double[] e = new double[] { 0.75, 0, 1, 0, 0, 0.75 };

    final double[] o = solver.invertDiagonal(a);

    Assertions.assertNotNull(o, "Failed to invert");

    if (logger.isLoggable(Level.INFO)) {
      logger.info(FunctionUtils.getSupplier("a diagonal = %s", Arrays.toString(o)));
    }
    Assertions.assertArrayEquals(e, o, 1e-4, "Bad inversion");
  }
  //@formatter:on

  private abstract class SolverTimingTask extends BaseTimingTask {
    DenseMatrix64F[] a;
    DenseMatrix64F[] b;
    final boolean badSolver;
    // No validation for a pure speed test
    EjmlLinearSolver solver = new EjmlLinearSolver();

    public SolverTimingTask(String name, DenseMatrix64F[] a, DenseMatrix64F[] b) {
      super(name + " " + a[0].numCols);
      // Clone the data
      this.a = a;
      this.b = b;
      // Check the solver gets a good answer
      solver.setEqual(new DoubleEquality(5e-3, 1e-6));
      int fail = 0;
      final Object data = getData(0);
      a = (DenseMatrix64F[]) ((Object[]) data)[0];
      b = (DenseMatrix64F[]) ((Object[]) data)[1];
      for (int i = 0; i < a.length; i++) {
        if (!solve(a[i], b[i])) {
          fail++;
        }
      }
      if (fail > 0) {
        logger
            .info(FunctionUtils.getSupplier(getName() + " failed to invert %d/%d", fail, a.length));
        badSolver = true;
      } else {
        badSolver = false;
      }
      solver.setEqual(null);
    }

    @Override
    public int getSize() {
      return 1;
    }

    @Override
    public Object getData(int i) {
      // Clone
      int n = b.length;
      final DenseMatrix64F[] a = new DenseMatrix64F[n];
      final DenseMatrix64F[] b = new DenseMatrix64F[n];
      while (n-- > 0) {
        a[n] = this.a[n].copy();
        b[n] = this.b[n].copy();
      }
      return new Object[] {a, b};
    }

    @Override
    public Object run(Object data) {
      final DenseMatrix64F[] a = (DenseMatrix64F[]) ((Object[]) data)[0];
      final DenseMatrix64F[] b = (DenseMatrix64F[]) ((Object[]) data)[1];
      for (int i = 0; i < a.length; i++) {
        solve(a[i], b[i]);
      }
      return null;
    }

    abstract boolean solve(DenseMatrix64F a, DenseMatrix64F b);
  }

  private class LinearSolverTimingTask extends SolverTimingTask {
    public LinearSolverTimingTask(DenseMatrix64F[] a, DenseMatrix64F[] b) {
      super("Linear Solver", a, b);
    }

    @Override
    boolean solve(DenseMatrix64F a, DenseMatrix64F b) {
      return solver.solveLinear(a, b);
    }
  }

  private class CholeskySolverTimingTask extends SolverTimingTask {
    public CholeskySolverTimingTask(DenseMatrix64F[] a, DenseMatrix64F[] b) {
      super("Cholesky Solver", a, b);
    }

    @Override
    boolean solve(DenseMatrix64F a, DenseMatrix64F b) {
      return solver.solveCholesky(a, b);
    }
  }

  private class CholeskyLDLTSolverTimingTask extends SolverTimingTask {
    public CholeskyLDLTSolverTimingTask(DenseMatrix64F[] a, DenseMatrix64F[] b) {
      super("CholeskyLDLT Solver", a, b);
    }

    @Override
    boolean solve(DenseMatrix64F a, DenseMatrix64F b) {
      return solver.solveCholeskyLDLT(a, b);
    }
  }

  private class PseudoInverseSolverTimingTask extends SolverTimingTask {
    public PseudoInverseSolverTimingTask(DenseMatrix64F[] a, DenseMatrix64F[] b) {
      super("PseudoInverse Solver", a, b);
    }

    @Override
    boolean solve(DenseMatrix64F a, DenseMatrix64F b) {
      return solver.solvePseudoInverse(a, b);
    }
  }

  private class DirectInversionSolverTimingTask extends SolverTimingTask {
    public DirectInversionSolverTimingTask(DenseMatrix64F[] a, DenseMatrix64F[] b) {
      super("DirectInversion Solver", a, b);
    }

    @Override
    boolean solve(DenseMatrix64F a, DenseMatrix64F b) {
      return solver.solveDirectInversion(a, b);
    }
  }

  // Create a speed test of the different methods
  @SpeedTag
  @SeededTest
  public void runSolverSpeedTest6(RandomSeed seed) {
    runSolverSpeedTest(seed, GaussianFunctionFactory.FIT_ERF_FREE_CIRCLE);
  }

  @SpeedTag
  @SeededTest
  public void runSolverSpeedTest5(RandomSeed seed) {
    runSolverSpeedTest(seed, GaussianFunctionFactory.FIT_ERF_CIRCLE);
  }

  @SpeedTag
  @SeededTest
  public void runSolverSpeedTest4(RandomSeed seed) {
    runSolverSpeedTest(seed, GaussianFunctionFactory.FIT_ERF_FIXED);
  }

  @SpeedTag
  @SeededTest
  public void runSolverSpeedTest3(RandomSeed seed) {
    runSolverSpeedTest(seed, GaussianFunctionFactory.FIT_SIMPLE_NB_FIXED);
  }

  @SpeedTag
  @SeededTest
  public void runSolverSpeedTest2(RandomSeed seed) {
    runSolverSpeedTest(seed, GaussianFunctionFactory.FIT_SIMPLE_NS_NB_FIXED);
  }

  private void runSolverSpeedTest(RandomSeed seed, int flags) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final Gaussian2DFunction f0 = GaussianFunctionFactory.create2D(1, 10, 10, flags, null);
    final int n = f0.size();
    final double[] y = new double[n];
    final TurboList<DenseMatrix64F> aList = new TurboList<>();
    final TurboList<DenseMatrix64F> bList = new TurboList<>();
    final double[] testbackground = new double[] {0.2, 0.7};
    final double[] testsignal1 = new double[] {30, 100, 300};
    final double[] testcx1 = new double[] {4.9, 5.3};
    final double[] testcy1 = new double[] {4.8, 5.2};
    final double[] testw1 = new double[] {1.1, 1.2, 1.5};
    final int np = f0.getNumberOfGradients();
    final GradientCalculator calc = GradientCalculatorUtils.newCalculator(np);
    final CustomPoissonDistribution pd = new CustomPoissonDistribution(
        new RandomGeneratorAdapter(RngUtils.create(seed.getSeedAsLong())), 1);
    // double lambda = 10;
    for (final double background : testbackground) {
      // Peak 1
      for (final double signal1 : testsignal1) {
        for (final double cx1 : testcx1) {
          for (final double cy1 : testcy1) {
            for (final double w1 : testw1) {
              final double[] p = new double[] {background, signal1, 0, cx1, cy1, w1, w1};
              f0.initialise(p);
              f0.forEach(new ValueProcedure() {
                int i = 0;

                @Override
                public void execute(double value) {
                  // Poisson data
                  pd.setMeanUnsafe(value);
                  y[i++] = pd.sample();
                }
              });
              final double[][] alpha = new double[np][np];
              final double[] beta = new double[np];
              // double ss =
              calc.findLinearised(n, y, p, alpha, beta, f0);
              // TestLog.fine(logger,"SS = %f", ss);
              // As per the LVM algorithm
              // for (int i = 0; i < np; i++)
              // alpha[i][i] *= lambda;
              aList.add(EjmlLinearSolver.toA(alpha));
              bList.add(EjmlLinearSolver.toB(beta));
            }
          }
        }
      }
    }

    final DenseMatrix64F[] a = aList.toArray(new DenseMatrix64F[aList.size()]);
    final DenseMatrix64F[] b = bList.toArray(new DenseMatrix64F[bList.size()]);
    final int runs = 100000 / a.length;
    final TimingService ts = new TimingService(runs);
    final TurboList<SolverTimingTask> tasks = new TurboList<>();
    // Added in descending speed order
    tasks.add(new PseudoInverseSolverTimingTask(a, b));
    tasks.add(new LinearSolverTimingTask(a, b));
    tasks.add(new CholeskySolverTimingTask(a, b));
    tasks.add(new CholeskyLDLTSolverTimingTask(a, b));
    tasks.add(new DirectInversionSolverTimingTask(a, b));
    for (final SolverTimingTask task : tasks) {
      if (!task.badSolver) {
        ts.execute(task);
      }
    }
    final int size = ts.getSize();
    ts.repeat();
    if (logger.isLoggable(Level.INFO)) {
      logger.info(ts.getReport(size));
    }

    // Since the speed is very similar at sizes 2-5 there is nothing to reliably assert
    // about the fastest of Cholesky/CholeskyLDLT/Direct.
    // Just check the PseudoInverse is slowest
    for (int i = 1; i < size; i++) {
      logger.log(TestLogUtils.getTimingRecord(ts.get(-(size)), ts.get(-i)));
    }

    if (np > 2) {
      // The Direct solver may not be faster at size=5
      int i = (np == 5) ? 2 : 1;
      final int size_1 = size - 1;
      for (; i < size_1; i++) {
        logger.log(TestLogUtils.getTimingRecord(ts.get(-(size_1)), ts.get(-i)));
      }
    }
  }

  private abstract class InversionTimingTask extends BaseTimingTask {
    DenseMatrix64F[] a;
    boolean[] ignore;
    final boolean badSolver;
    // No validation for a pure speed test
    EjmlLinearSolver solver = new EjmlLinearSolver();

    public InversionTimingTask(String name, DenseMatrix64F[] a, boolean[] ignore,
        double[][] answer) {
      super(name + " " + a[0].numCols);
      // Clone the data
      this.a = a;
      this.ignore = ignore;
      // Check the inversion gets a good answer
      solver.setInversionTolerance(1e-2);
      int fail = 0;
      final boolean[] failed = new boolean[a.length];
      for (int i = 0; i < a.length; i++) {
        final double[] b = invert(a[i].copy());
        if (b == null) {
          failed[i] = true;
          fail++;
        } else if (answer[i] == null) {
          answer[i] = b;
        } else {
          // Check against the existing answer
          for (int j = 0; j < b.length; j++) {
            if (!DoubleEquality.almostEqualRelativeOrAbsolute(b[j], answer[i][j], 1e-3, 1e-4)) {
              failed[i] = true;
              fail++;
              break;
            }
          }
        }
      }
      if (fail > 0) {
        logger
            .info(FunctionUtils.getSupplier(getName() + " failed to invert %d/%d", fail, a.length));
      }
      solver.setInversionTolerance(0);
      if (fail == a.length) {
        // This solver cannot do the inversion
        badSolver = true;
      } else {
        badSolver = false;
        // Flag those we cannot do as bad
        for (int i = 0; i < a.length; i++) {
          if (failed[i]) {
            ignore[i] = true;
          }
        }
      }
    }

    @Override
    public int getSize() {
      return 1;
    }

    @Override
    public Object getData(int i) {
      // Clone
      int n = a.length;
      final DenseMatrix64F[] a = new DenseMatrix64F[n];
      while (n-- > 0) {
        if (!ignore[n]) {
          a[n] = this.a[n].copy();
        }
      }
      return a;
    }

    @Override
    public Object run(Object data) {
      final DenseMatrix64F[] a = (DenseMatrix64F[]) data;
      for (int i = 0; i < a.length; i++) {
        if (!ignore[i]) {
          invert(a[i]);
        }
      }
      return null;
    }

    abstract double[] invert(DenseMatrix64F a);

    double[] extract(DenseMatrix64F a) {
      final int n = a.numCols;
      final double[] b = new double[n];
      for (int i = 0, j = 0; i < n; i++, j += n + 1) {
        b[i] = a.data[j];
      }
      return b;
    }
  }

  private class LinearInversionTimingTask extends InversionTimingTask {
    public LinearInversionTimingTask(DenseMatrix64F[] a, boolean[] ignore, double[][] answer) {
      super("Linear Inversion", a, ignore, answer);
    }

    @Override
    double[] invert(DenseMatrix64F a) {
      if (solver.invertLinear(a)) {
        return extract(a);
      }
      return null;
    }
  }

  private class CholeskyInversionTimingTask extends InversionTimingTask {
    public CholeskyInversionTimingTask(DenseMatrix64F[] a, boolean[] ignore, double[][] answer) {
      super("Cholesky Inversion", a, ignore, answer);
    }

    @Override
    double[] invert(DenseMatrix64F a) {
      if (solver.invertCholesky(a)) {
        return extract(a);
      }
      return null;
    }
  }

  private class CholeskyLDLTInversionTimingTask extends InversionTimingTask {
    public CholeskyLDLTInversionTimingTask(DenseMatrix64F[] a, boolean[] ignore,
        double[][] answer) {
      super("CholeskyLDLT Inversion", a, ignore, answer);
    }

    @Override
    double[] invert(DenseMatrix64F a) {
      if (solver.invertCholeskyLdlT(a)) {
        return extract(a);
      }
      return null;
    }
  }

  private class PseudoInverseInversionTimingTask extends InversionTimingTask {
    public PseudoInverseInversionTimingTask(DenseMatrix64F[] a, boolean[] ignore,
        double[][] answer) {
      super("PseudoInverse Inversion", a, ignore, answer);
    }

    @Override
    double[] invert(DenseMatrix64F a) {
      if (solver.invertPseudoInverse(a)) {
        return extract(a);
      }
      return null;
    }
  }

  private class DirectInversionInversionTimingTask extends InversionTimingTask {
    public DirectInversionInversionTimingTask(DenseMatrix64F[] a, boolean[] ignore,
        double[][] answer) {
      super("DirectInversion Inversion", a, ignore, answer);
    }

    @Override
    double[] invert(DenseMatrix64F a) {
      if (solver.invertDirectInversion(a)) {
        return extract(a);
      }
      return null;
    }
  }

  private class DiagonalDirectInversionInversionTimingTask extends InversionTimingTask {
    public DiagonalDirectInversionInversionTimingTask(DenseMatrix64F[] a, boolean[] ignore,
        double[][] answer) {
      super("DiagonalDirectInversion Inversion", a, ignore, answer);
    }

    @Override
    double[] invert(DenseMatrix64F a) {
      return EjmlLinearSolver.invertDiagonalDirectInversion(a);
    }
  }

  // Create a speed test of the different methods
  @SpeedTag
  @SeededTest
  public void runInversionSpeedTest6(RandomSeed seed) {
    runInversionSpeedTest(seed, GaussianFunctionFactory.FIT_ERF_FREE_CIRCLE);
  }

  @SpeedTag
  @SeededTest
  public void runInversionSpeedTest5(RandomSeed seed) {
    runInversionSpeedTest(seed, GaussianFunctionFactory.FIT_ERF_CIRCLE);
  }

  @SpeedTag
  @SeededTest
  public void runInversionSpeedTest4(RandomSeed seed) {
    runInversionSpeedTest(seed, GaussianFunctionFactory.FIT_ERF_FIXED);
  }

  @SpeedTag
  @SeededTest
  public void runInversionSpeedTest3(RandomSeed seed) {
    runInversionSpeedTest(seed, GaussianFunctionFactory.FIT_SIMPLE_NB_FIXED);
  }

  @SpeedTag
  @SeededTest
  public void runInversionSpeedTest2(RandomSeed seed) {
    runInversionSpeedTest(seed, GaussianFunctionFactory.FIT_SIMPLE_NS_NB_FIXED);
  }

  private void runInversionSpeedTest(RandomSeed seed, int flags) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final Gaussian2DFunction f0 = GaussianFunctionFactory.create2D(1, 10, 10, flags, null);
    final int n = f0.size();
    final double[] y = new double[n];
    final TurboList<DenseMatrix64F> aList = new TurboList<>();
    final double[] testbackground = new double[] {0.2, 0.7};
    final double[] testsignal1 = new double[] {30, 100, 300};
    final double[] testcx1 = new double[] {4.9, 5.3};
    final double[] testcy1 = new double[] {4.8, 5.2};
    final double[] testw1 = new double[] {1.1, 1.2, 1.5};
    final int np = f0.getNumberOfGradients();
    final GradientCalculator calc = GradientCalculatorUtils.newCalculator(np);
    final CustomPoissonDistribution pd = new CustomPoissonDistribution(
        new RandomGeneratorAdapter(RngUtils.create(seed.getSeedAsLong())), 1);
    // double lambda = 10;
    for (final double background : testbackground) {
      // Peak 1
      for (final double signal1 : testsignal1) {
        for (final double cx1 : testcx1) {
          for (final double cy1 : testcy1) {
            for (final double w1 : testw1) {
              final double[] p = new double[] {background, signal1, 0, cx1, cy1, w1, w1};
              f0.initialise(p);
              f0.forEach(new ValueProcedure() {
                int i = 0;

                @Override
                public void execute(double value) {
                  // Poisson data
                  pd.setMeanUnsafe(value);
                  y[i++] = pd.sample();
                }
              });
              final double[][] alpha = new double[np][np];
              final double[] beta = new double[np];
              // double ss =
              calc.findLinearised(n, y, p, alpha, beta, f0);
              // TestLog.fine(logger,"SS = %f", ss);
              // As per the LVM algorithm
              // for (int i = 0; i < np; i++)
              // alpha[i][i] *= lambda;
              aList.add(EjmlLinearSolver.toA(alpha));
            }
          }
        }
      }
    }

    final DenseMatrix64F[] a = aList.toArray(new DenseMatrix64F[aList.size()]);
    final boolean[] ignore = new boolean[a.length];
    final double[][] answer = new double[a.length][];
    final int runs = 100000 / a.length;
    final TimingService ts = new TimingService(runs);
    final TurboList<InversionTimingTask> tasks = new TurboList<>();
    // Added in descending speed order
    tasks.add(new PseudoInverseInversionTimingTask(a, ignore, answer));
    tasks.add(new LinearInversionTimingTask(a, ignore, answer));
    tasks.add(new CholeskyLDLTInversionTimingTask(a, ignore, answer));
    tasks.add(new CholeskyInversionTimingTask(a, ignore, answer));
    tasks.add(new DirectInversionInversionTimingTask(a, ignore, answer));
    tasks.add(new DiagonalDirectInversionInversionTimingTask(a, ignore, answer));
    for (final InversionTimingTask task : tasks) {
      if (!task.badSolver) {
        ts.execute(task);
      }
    }
    final int size = ts.getSize();
    ts.repeat();
    if (logger.isLoggable(Level.INFO)) {
      logger.info(ts.getReport(size));
    }

    // When it is present the DiagonalDirect is fastest (n<=5)
    if (np <= 5) {
      for (int i = 2; i <= size; i++) {
        logger.log(TestLogUtils.getTimingRecord(ts.get(-i), ts.get(-1)));
      }

      if (np < 5) {
        // n < 5 Direct is fastest
        for (int i = 3; i <= size; i++) {
          logger.log(TestLogUtils.getTimingRecord(ts.get(-i), ts.get(-2)));
        }
      } else {
        // Cholesky should be fastest. It is marginal over CholeskyLDLT.
        // and may not be faster than Direct at n=5 so that comparison is ignored.
        for (int i = 4; i <= size; i++) {
          logger.log(TestLogUtils.getTimingRecord(ts.get(-i), ts.get(-3)));
        }
      }
    } else {
      // No Direct inversion possible.
      // Cholesky should be fastest.
      for (int i = 2; i <= size; i++) {
        logger.log(TestLogUtils.getTimingRecord(ts.get(-i), ts.get(-1)));
      }
    }
  }
}
