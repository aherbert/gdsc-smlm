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
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.rng.UniformRandomProvider;
import org.ejml.data.DenseMatrix64F;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.function.Executable;
import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.LocalList;
import uk.ac.sussex.gdsc.smlm.GdscSmlmTestUtils;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.GradientCalculator;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.GradientCalculatorUtils;
import uk.ac.sussex.gdsc.smlm.function.ValueProcedure;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import uk.ac.sussex.gdsc.test.api.TestHelper;
import uk.ac.sussex.gdsc.test.api.function.DoubleDoubleBiPredicate;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.junit5.SpeedTag;
import uk.ac.sussex.gdsc.test.rng.RngUtils;
import uk.ac.sussex.gdsc.test.utils.RandomSeed;
import uk.ac.sussex.gdsc.test.utils.TestLogUtils.TestLevel;
import uk.ac.sussex.gdsc.test.utils.functions.FunctionUtils;

@SuppressWarnings({"javadoc"})
class EjmlLinearSolverTest {
  /** Used to control test logging. Actual output messages are written at the INFO level. */
  private static final Level LOG_LEVEL = TestLevel.TEST_INFO;
  private static Logger logger;

  @BeforeAll
  public static void beforeAll() {
    logger = Logger.getLogger(EjmlLinearSolverTest.class.getName());
  }

  @AfterAll
  public static void afterAll() {
    logger = null;
  }

  //@formatter:off
  @Test
  void canSolveLinearEquation()
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

    if (logger.isLoggable(LOG_LEVEL)) {
      logger.log(LOG_LEVEL, FunctionUtils.getSupplier("x = %s", Arrays.toString(b)));
    }
    for (int i = 0; i < b.length; i++)
    {
      if (logger.isLoggable(LOG_LEVEL)) {
        logger.log(LOG_LEVEL, FunctionUtils.getSupplier("a[%d] = %s", i, Arrays.toString(a[i])));
      }
      Assertions.assertArrayEquals(a_inv[i], a[i], 1e-4f, "Bad inversion");
    }
  }

  @Test
  void canSolveLinearEquationWithZeroInB()
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

    if (logger.isLoggable(LOG_LEVEL)) {
      logger.log(LOG_LEVEL, FunctionUtils.getSupplier("x = %s", Arrays.toString(b)));
    }
    for (int i = 0; i < b.length; i++)
    {
      if (logger.isLoggable(LOG_LEVEL)) {
        logger.log(LOG_LEVEL, FunctionUtils.getSupplier("a[%d] = %s", i, Arrays.toString(a[i])));
      }
      Assertions.assertArrayEquals(a_inv[i], a[i], 1e-4f, "Bad inversion");
    }
  }

  @Test
  void canSolveLinearEquationWithZeroInA()
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

    if (logger.isLoggable(LOG_LEVEL)) {
      logger.log(LOG_LEVEL, FunctionUtils.getSupplier("x = %s", Arrays.toString(b)));
    }
    for (int i = 0; i < b.length; i++)
    {
      if (logger.isLoggable(LOG_LEVEL)) {
        logger.log(LOG_LEVEL, FunctionUtils.getSupplier("a[%d] = %s", i, Arrays.toString(a[i])));
      }
      Assertions.assertArrayEquals(a_inv[i], a[i], 1e-4f, "Bad inversion");
    }
  }

  @Test
  void canSolveLinearEquationWithZerosInA()
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

    if (logger.isLoggable(LOG_LEVEL)) {
      logger.log(LOG_LEVEL, FunctionUtils.getSupplier("x = %s", Arrays.toString(b)));
    }
    for (int i = 0; i < b.length; i++)
    {
      if (logger.isLoggable(LOG_LEVEL)) {
        logger.log(LOG_LEVEL, FunctionUtils.getSupplier("a[%d] = %s", i, Arrays.toString(a[i])));
      }
      Assertions.assertArrayEquals(a_inv[i], a[i], 1e-4f, "Bad inversion");
    }
  }

  @Test
  void canInvert()
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
      if (logger.isLoggable(LOG_LEVEL)) {
        logger.log(LOG_LEVEL, FunctionUtils.getSupplier("a[%d] = %s", i, Arrays.toString(a[i])));
      }
      Assertions.assertArrayEquals(a_inv[i], a[i], 1e-4f, "Bad inversion");
    }
  }

  @Test
  void canInvertWithZeros()
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
      if (logger.isLoggable(LOG_LEVEL)) {
        logger.log(LOG_LEVEL, FunctionUtils.getSupplier("a[%d] = %s", i, Arrays.toString(a[i])));
      }
      Assertions.assertArrayEquals(a_inv[i], a[i], 1e-4f, "Bad inversion");
    }
  }

  @Test
  void canInvertDiagonal()
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

    if (logger.isLoggable(LOG_LEVEL)) {
      logger.log(LOG_LEVEL, FunctionUtils.getSupplier("a diagonal = %s", Arrays.toString(o)));
    }
    Assertions.assertArrayEquals(e, o, 1e-4, "Bad inversion");
  }

  @Test
  void canInvertDiagonalWithZeros()
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

    if (logger.isLoggable(LOG_LEVEL)) {
      logger.log(LOG_LEVEL, FunctionUtils.getSupplier("a diagonal = %s", Arrays.toString(o)));
    }
    Assertions.assertArrayEquals(e, o, 1e-4, "Bad inversion");
  }
  //@formatter:on

  // Allow matrix name A and vector names b & x for equation A x = b
  // CHECKSTYLE.OFF: MemberName
  // CHECKSTYLE.OFF: ParameterName

  private abstract class SolverExecutable implements Executable {
    String name;
    DenseMatrix64F[] a;
    DenseMatrix64F[] b;

    SolverExecutable(String name, DenseMatrix64F[] a, DenseMatrix64F[] b) {
      this.name = name + " " + a[0].numCols;
      this.a = a;
      this.b = b;
    }

    @Override
    public void execute() {
      // Use the solver validation
      EjmlLinearSolver solver = new EjmlLinearSolver();
      solver.setEqual(new DoubleEquality(5e-3, 1e-6));
      int fail = 0;
      for (int i = 0; i < a.length; i++) {
        DenseMatrix64F aa = a[i].copy();
        DenseMatrix64F bb = b[i].copy();
        if (!solve(solver, aa, bb)) {
          fail++;
        }
      }
      Assertions.assertEquals(0, fail,
          FunctionUtils.getSupplier("%s failed to invert %d/%d", name, fail, a.length));
    }

    abstract boolean solve(EjmlLinearSolver solver, DenseMatrix64F a, DenseMatrix64F b);
  }

  private class LinearSolverExecutable extends SolverExecutable {
    public LinearSolverExecutable(DenseMatrix64F[] a, DenseMatrix64F[] b) {
      super("Linear Solver", a, b);
    }

    @Override
    boolean solve(EjmlLinearSolver solver, DenseMatrix64F a, DenseMatrix64F b) {
      return solver.solveLinear(a, b);
    }
  }

  private class CholeskySolverExecutable extends SolverExecutable {
    public CholeskySolverExecutable(DenseMatrix64F[] a, DenseMatrix64F[] b) {
      super("Cholesky Solver", a, b);
    }

    @Override
    boolean solve(EjmlLinearSolver solver, DenseMatrix64F a, DenseMatrix64F b) {
      return solver.solveCholesky(a, b);
    }
  }

  private class CholeskyLdltSolverExecutable extends SolverExecutable {
    public CholeskyLdltSolverExecutable(DenseMatrix64F[] a, DenseMatrix64F[] b) {
      super("CholeskyLDLT Solver", a, b);
    }

    @Override
    boolean solve(EjmlLinearSolver solver, DenseMatrix64F a, DenseMatrix64F b) {
      return solver.solveCholeskyLdlT(a, b);
    }
  }

  private class PseudoInverseSolverExecutable extends SolverExecutable {
    public PseudoInverseSolverExecutable(DenseMatrix64F[] a, DenseMatrix64F[] b) {
      super("PseudoInverse Solver", a, b);
    }

    @Override
    boolean solve(EjmlLinearSolver solver, DenseMatrix64F a, DenseMatrix64F b) {
      return solver.solvePseudoInverse(a, b);
    }
  }

  private class DirectInversionSolverExecutable extends SolverExecutable {
    public DirectInversionSolverExecutable(DenseMatrix64F[] a, DenseMatrix64F[] b) {
      super("DirectInversion Solver", a, b);
    }

    @Override
    boolean solve(EjmlLinearSolver solver, DenseMatrix64F a, DenseMatrix64F b) {
      return solver.solveDirectInversion(a, b);
    }
  }

  @SpeedTag
  @SeededTest
  void runSolverTest6(RandomSeed seed) {
    runSolverTest(seed, GaussianFunctionFactory.FIT_ERF_FREE_CIRCLE);
  }

  @SpeedTag
  @SeededTest
  void runSolverTest5(RandomSeed seed) {
    runSolverTest(seed, GaussianFunctionFactory.FIT_ERF_CIRCLE);
  }

  @SpeedTag
  @SeededTest
  void runSolverTest4(RandomSeed seed) {
    runSolverTest(seed, GaussianFunctionFactory.FIT_ERF_FIXED);
  }

  @SpeedTag
  @SeededTest
  void runSolverTest3(RandomSeed seed) {
    runSolverTest(seed, GaussianFunctionFactory.FIT_SIMPLE_NB_FIXED);
  }

  @SpeedTag
  @SeededTest
  void runSolverTest2(RandomSeed seed) {
    runSolverTest(seed, GaussianFunctionFactory.FIT_SIMPLE_NS_NB_FIXED);
  }

  private void runSolverTest(RandomSeed seed, int flags) {
    final Gaussian2DFunction f0 = GaussianFunctionFactory.create2D(1, 10, 10, flags, null);
    final int n = f0.size();
    final double[] y = new double[n];
    final LocalList<DenseMatrix64F> aList = new LocalList<>();
    final LocalList<DenseMatrix64F> bList = new LocalList<>();
    final double[] testbackground = new double[] {0.2, 0.7};
    final double[] testsignal1 = new double[] {30, 100, 300};
    final double[] testcx1 = new double[] {4.9, 5.3};
    final double[] testcy1 = new double[] {4.8, 5.2};
    final double[] testw1 = new double[] {1.1, 1.2, 1.5};
    final int np = f0.getNumberOfGradients();
    final GradientCalculator calc = GradientCalculatorUtils.newCalculator(np);
    final UniformRandomProvider rng = RngUtils.create(seed.get());
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
                int index = 0;

                @Override
                public void execute(double value) {
                  // Poisson data
                  y[index++] = GdscSmlmTestUtils.createPoissonSampler(rng, value).sample();
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

    final DenseMatrix64F[] a = aList.toArray(new DenseMatrix64F[0]);
    final DenseMatrix64F[] b = bList.toArray(new DenseMatrix64F[0]);
    final List<Executable> list = new ArrayList<>(
        Arrays.asList(new PseudoInverseSolverExecutable(a, b), new LinearSolverExecutable(a, b),
            new CholeskySolverExecutable(a, b), new CholeskyLdltSolverExecutable(a, b)));
    if (np <= 5) {
      list.add(new DirectInversionSolverExecutable(a, b));
    }
    Assertions.assertAll(list);
  }

  private abstract class InversionExecutable implements Executable {
    String name;
    DenseMatrix64F[] a;
    double[][] answer;

    public InversionExecutable(String name, DenseMatrix64F[] a, double[][] answer) {
      this.name = name + " " + a[0].numCols;
      this.a = a;
      this.answer = answer;
    }

    @Override
    public void execute() {
      // Use the solver validation
      EjmlLinearSolver solver = new EjmlLinearSolver();
      solver.setInversionTolerance(1e-2);
      DoubleDoubleBiPredicate test = TestHelper.doublesAreClose(1e-3, 1e-4);
      int fail = 0;
      for (int i = 0; i < a.length; i++) {
        DenseMatrix64F aa = a[i].copy();
        double[] b = invert(solver, aa);

        // Check against the answer
        if (answer[i] == null) {
          // Not initialised yet
          answer[i] = b;
        } else if (b == null) {
          fail++;
        } else {
          // Check against the existing answer
          for (int j = 0; j < b.length; j++) {
            if (!test.test(b[j], answer[i][j])) {
              fail++;
              break;
            }
          }
        }
      }
      Assertions.assertEquals(0, fail,
          FunctionUtils.getSupplier("%s failed to invert %d/%d", name, fail, a.length));
    }

    abstract double[] invert(EjmlLinearSolver solver, DenseMatrix64F a);

    final double[] extract(DenseMatrix64F a) {
      final int n = a.numCols;
      final double[] b = new double[n];
      for (int i = 0, j = 0; i < n; i++, j += n + 1) {
        b[i] = a.data[j];
      }
      return b;
    }
  }


  private class LinearInversionExecutable extends InversionExecutable {
    public LinearInversionExecutable(DenseMatrix64F[] a, double[][] answer) {
      super("Linear Inversion", a, answer);
    }

    @Override
    double[] invert(EjmlLinearSolver solver, DenseMatrix64F a) {
      if (solver.invertLinear(a)) {
        return extract(a);
      }
      return null;
    }
  }


  private class CholeskyInversionExecutable extends InversionExecutable {
    public CholeskyInversionExecutable(DenseMatrix64F[] a, double[][] answer) {
      super("Cholesky Inversion", a, answer);
    }

    @Override
    double[] invert(EjmlLinearSolver solver, DenseMatrix64F a) {
      if (solver.invertCholesky(a)) {
        return extract(a);
      }
      return null;
    }
  }


  private class CholeskyLdltInversionExecutable extends InversionExecutable {
    public CholeskyLdltInversionExecutable(DenseMatrix64F[] a, double[][] answer) {
      super("CholeskyLDLT Inversion", a, answer);
    }

    @Override
    double[] invert(EjmlLinearSolver solver, DenseMatrix64F a) {
      if (solver.invertCholeskyLdlT(a)) {
        return extract(a);
      }
      return null;
    }
  }


  private class PseudoInverseInversionExecutable extends InversionExecutable {
    public PseudoInverseInversionExecutable(DenseMatrix64F[] a, double[][] answer) {
      super("PseudoInverse Inversion", a, answer);
    }

    @Override
    double[] invert(EjmlLinearSolver solver, DenseMatrix64F a) {
      if (solver.invertPseudoInverse(a)) {
        return extract(a);
      }
      return null;
    }
  }


  private class DirectInversionInversionExecutable extends InversionExecutable {
    public DirectInversionInversionExecutable(DenseMatrix64F[] a, double[][] answer) {
      super("DirectInversion Inversion", a, answer);
    }

    @Override
    double[] invert(EjmlLinearSolver solver, DenseMatrix64F a) {
      if (solver.invertDirectInversion(a)) {
        return extract(a);
      }
      return null;
    }
  }


  private class DiagonalDirectInversionInversionExecutable extends InversionExecutable {
    public DiagonalDirectInversionInversionExecutable(DenseMatrix64F[] a, double[][] answer) {
      super("DiagonalDirectInversion Inversion", a, answer);
    }

    @Override
    double[] invert(EjmlLinearSolver solver, DenseMatrix64F a) {
      return EjmlLinearSolver.invertDiagonalDirectInversion(a);
    }

  }

  // Create a speed test of the different methods
  @SpeedTag
  @SeededTest
  void runInversionTest6(RandomSeed seed) {
    runInversionTest(seed, GaussianFunctionFactory.FIT_ERF_FREE_CIRCLE);
  }

  @SpeedTag
  @SeededTest
  void runInversionTest5(RandomSeed seed) {
    runInversionTest(seed, GaussianFunctionFactory.FIT_ERF_CIRCLE);
  }

  @SpeedTag
  @SeededTest
  void runInversionTest4(RandomSeed seed) {
    runInversionTest(seed, GaussianFunctionFactory.FIT_ERF_FIXED);
  }

  @SpeedTag
  @SeededTest
  void runInversionTest3(RandomSeed seed) {
    runInversionTest(seed, GaussianFunctionFactory.FIT_SIMPLE_NB_FIXED);
  }

  @SpeedTag
  @SeededTest
  void runInversionTest2(RandomSeed seed) {
    runInversionTest(seed, GaussianFunctionFactory.FIT_SIMPLE_NS_NB_FIXED);
  }

  private void runInversionTest(RandomSeed seed, int flags) {
    final Gaussian2DFunction f0 = GaussianFunctionFactory.create2D(1, 10, 10, flags, null);
    final int n = f0.size();
    final double[] y = new double[n];
    final LocalList<DenseMatrix64F> aList = new LocalList<>();
    final double[] testbackground = new double[] {0.2, 0.7};
    final double[] testsignal1 = new double[] {30, 100, 300};
    final double[] testcx1 = new double[] {4.9, 5.3};
    final double[] testcy1 = new double[] {4.8, 5.2};
    final double[] testw1 = new double[] {1.1, 1.2, 1.5};
    final int np = f0.getNumberOfGradients();
    final GradientCalculator calc = GradientCalculatorUtils.newCalculator(np);
    final UniformRandomProvider rng = RngUtils.create(seed.get());
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
                int index = 0;

                @Override
                public void execute(double value) {
                  // Poisson data
                  y[index++] = GdscSmlmTestUtils.createPoissonSampler(rng, value).sample();
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

    final DenseMatrix64F[] a = aList.toArray(new DenseMatrix64F[0]);
    final double[][] answer = new double[a.length][];

    // Get the actual answer
    new LinearInversionExecutable(a, answer).execute();
    Assertions.assertTrue(Arrays.stream(answer).noneMatch(x -> x == null),
        "Cannot solve all inversions");

    final List<Executable> list =
        new ArrayList<>(Arrays.asList(new PseudoInverseInversionExecutable(a, answer),
            new LinearInversionExecutable(a, answer),
            new CholeskyLdltInversionExecutable(a, answer),
            new CholeskyInversionExecutable(a, answer)));
    if (np <= 5) {
      list.add(new DirectInversionInversionExecutable(a, answer));
      list.add(new DiagonalDirectInversionInversionExecutable(a, answer));
    }
    Assertions.assertAll(list);
  }
}
