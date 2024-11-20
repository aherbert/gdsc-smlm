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

package uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient;

import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.LogRecord;
import java.util.logging.Logger;
import org.apache.commons.numbers.core.Precision;
import org.apache.commons.rng.UniformRandomProvider;
import org.ejml.data.DMatrixRMaj;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.BeforeAll;
import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.Statistics;
import uk.ac.sussex.gdsc.smlm.GdscSmlmTestUtils;
import uk.ac.sussex.gdsc.smlm.fitting.linear.EjmlLinearSolver;
import uk.ac.sussex.gdsc.smlm.function.DummyGradientFunction;
import uk.ac.sussex.gdsc.smlm.function.FakeGradientFunction;
import uk.ac.sussex.gdsc.smlm.function.Gradient1Function;
import uk.ac.sussex.gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import uk.ac.sussex.gdsc.smlm.function.gaussian.erf.ErfGaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.erf.SingleFreeCircularErfGaussian2DFunction;
import uk.ac.sussex.gdsc.test.api.Predicates;
import uk.ac.sussex.gdsc.test.api.TestAssertions;
import uk.ac.sussex.gdsc.test.api.function.DoubleDoubleBiPredicate;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.junit5.SpeedTag;
import uk.ac.sussex.gdsc.test.rng.RngFactory;
import uk.ac.sussex.gdsc.test.utils.AssertionErrorCounter;
import uk.ac.sussex.gdsc.test.utils.RandomSeed;
import uk.ac.sussex.gdsc.test.utils.TestComplexity;
import uk.ac.sussex.gdsc.test.utils.TestLogging;
import uk.ac.sussex.gdsc.test.utils.TestLogging.TestLevel;
import uk.ac.sussex.gdsc.test.utils.TestSettings;
import uk.ac.sussex.gdsc.test.utils.functions.IndexSupplier;

/**
 * Contains speed tests for the fastest method for calculating the Hessian and gradient vector for
 * use in the LVM algorithm.
 *
 * <p>Note: This class is a test-bed for implementation strategies. The fastest strategy can then be
 * used for other gradient procedures.
 */
@SuppressWarnings({"javadoc"})
class LsqLvmGradientProcedureTest {
  private static Logger logger;

  @BeforeAll
  public static void beforeAll() {
    logger = Logger.getLogger(LsqLvmGradientProcedureTest.class.getName());
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

  /**
   * Create a gradient procedure.
   */
  interface BaseLsqLvmGradientProcedureFactory {
    /**
     * Create a LSQ LVM procedure.
     *
     * <p>Instance methods for testing.
     *
     * @param y the y
     * @param func the function
     * @return the LSQ LVM gradient procedure
     */
    default BaseLsqLvmGradientProcedure createProcedure(final double[] y,
        final Gradient1Function func) {
      return createProcedure(y, null, func);
    }

    /**
     * Create a LSQ LVM procedure.
     *
     * <p>Instance methods for testing.
     *
     * @param y the y
     * @param baseline the baseline
     * @param func the function
     * @return the LSQ LVM gradient procedure
     */
    BaseLsqLvmGradientProcedure createProcedure(final double[] y, double[] baseline,
        final Gradient1Function func);
  }

  private static double random(UniformRandomProvider rng, double value) {
    return value - value * 0.1 + rng.nextDouble() * 0.2;
  }

  @SeededTest
  void gradientProcedureFactoryCreatesOptimisedProcedures() {
    final double[] y = new double[0];
    Assertions.assertEquals(
        LsqLvmGradientProcedureMatrixUtils.create(y, null, new DummyGradientFunction(6)).getClass(),
        LsqLvmGradientProcedureMatrix6.class);
    Assertions.assertEquals(
        LsqLvmGradientProcedureMatrixUtils.create(y, null, new DummyGradientFunction(5)).getClass(),
        LsqLvmGradientProcedureMatrix5.class);
    Assertions.assertEquals(
        LsqLvmGradientProcedureMatrixUtils.create(y, null, new DummyGradientFunction(4)).getClass(),
        LsqLvmGradientProcedureMatrix4.class);

    Assertions.assertEquals(
        LsqLvmGradientProcedureLinearUtils.create(y, null, new DummyGradientFunction(6)).getClass(),
        LsqLvmGradientProcedureLinear6.class);
    Assertions.assertEquals(
        LsqLvmGradientProcedureLinearUtils.create(y, null, new DummyGradientFunction(5)).getClass(),
        LsqLvmGradientProcedureLinear5.class);
    Assertions.assertEquals(
        LsqLvmGradientProcedureLinearUtils.create(y, null, new DummyGradientFunction(4)).getClass(),
        LsqLvmGradientProcedureLinear4.class);
  }

  @SeededTest
  void gradientProcedureLinearComputesSameAsGradientCalculator(RandomSeed seed) {
    gradientProcedureComputesSameAsGradientCalculator(seed,
        LsqLvmGradientProcedureLinearUtils::create);
  }

  @SeededTest
  void gradientProcedureMatrixComputesSameAsGradientCalculator(RandomSeed seed) {
    gradientProcedureComputesSameAsGradientCalculator(seed,
        LsqLvmGradientProcedureMatrixUtils::create);
  }

  @SeededTest
  void gradientProcedureComputesSameAsGradientCalculator(RandomSeed seed) {
    gradientProcedureComputesSameAsGradientCalculator(seed, LsqLvmGradientProcedureUtils::create);
  }

  private void gradientProcedureComputesSameAsGradientCalculator(RandomSeed seed,
      BaseLsqLvmGradientProcedureFactory factory) {
    gradientProcedureComputesSameAsGradientCalculator(seed, 4, factory);
    gradientProcedureComputesSameAsGradientCalculator(seed, 5, factory);
    gradientProcedureComputesSameAsGradientCalculator(seed, 6, factory);
    gradientProcedureComputesSameAsGradientCalculator(seed, 11, factory);
    gradientProcedureComputesSameAsGradientCalculator(seed, 21, factory);
  }

  private void gradientProcedureComputesSameAsGradientCalculator(RandomSeed seed, int nparams,
      BaseLsqLvmGradientProcedureFactory factory) {
    final int iter = 10;

    final double[][] alpha = new double[nparams][nparams];
    final double[] beta = new double[nparams];

    final ArrayList<double[]> paramsList = new ArrayList<>(iter);
    final ArrayList<double[]> yList = new ArrayList<>(iter);

    final int[] x = createFakeData(RngFactory.create(seed.get()), nparams, iter, paramsList, yList);
    final int n = x.length;
    final FakeGradientFunction func = new FakeGradientFunction(blockWidth, nparams);

    final GradientCalculator calc = GradientCalculatorUtils.newCalculator(nparams, false);

    final String name = factory.getClass().getSimpleName();
    // Create messages
    final IndexSupplier msgR = new IndexSupplier(1, name + "Result: Not same ", null);
    final IndexSupplier msgOb = new IndexSupplier(1, name + "Observations: Not same beta ", null);
    final IndexSupplier msgOal =
        new IndexSupplier(1, name + "Observations: Not same alpha linear ", null);
    final IndexSupplier msgOam =
        new IndexSupplier(1, name + "Observations: Not same alpha matrix ", null);

    for (int i = 0; i < paramsList.size(); i++) {
      final BaseLsqLvmGradientProcedure p = factory.createProcedure(yList.get(i), func);
      p.gradient(paramsList.get(i));
      final double s = p.value;
      final double s2 = calc.findLinearised(n, yList.get(i), paramsList.get(i), alpha, beta, func);
      // Exactly the same ...
      Assertions.assertEquals(s, s2, msgR.set(0, i));
      Assertions.assertArrayEquals(p.beta, beta, msgOb.set(0, i));

      final double[] al = p.getAlphaLinear();
      Assertions.assertArrayEquals(al, new DMatrixRMaj(alpha).data, msgOal.set(0, i));

      final double[][] am = p.getAlphaMatrix();
      Assertions.assertArrayEquals(am, alpha, msgOam.set(0, i));
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
      long nanos = System.nanoTime();
      run();
      nanos = System.nanoTime() - nanos;
      // logger.fine(FormatSupplier.getSupplier("[%d] Time = %d", loops, t);
      return nanos;
    }

    abstract void run();
  }

  @SpeedTag
  @SeededTest
  void gradientProcedureLinearIsNotSlowerThanGradientCalculator(RandomSeed seed) {
    gradientProcedureIsNotSlowerThanGradientCalculator(seed,
        LsqLvmGradientProcedureLinearUtils::create);
  }

  @SpeedTag
  @SeededTest
  void gradientProcedureMatrixIsNotSlowerThanGradientCalculator(RandomSeed seed) {
    gradientProcedureIsNotSlowerThanGradientCalculator(seed,
        LsqLvmGradientProcedureMatrixUtils::create);
  }

  @SpeedTag
  @SeededTest
  void gradientProcedureIsNotSlowerThanGradientCalculator(RandomSeed seed) {
    gradientProcedureIsNotSlowerThanGradientCalculator(seed, LsqLvmGradientProcedureUtils::create);
  }

  private void gradientProcedureIsNotSlowerThanGradientCalculator(RandomSeed seed,
      BaseLsqLvmGradientProcedureFactory factory) {
    gradientProcedureIsNotSlowerThanGradientCalculator(seed, 4, factory);
    gradientProcedureIsNotSlowerThanGradientCalculator(seed, 5, factory);
    gradientProcedureIsNotSlowerThanGradientCalculator(seed, 6, factory);
    // 2 peaks
    gradientProcedureIsNotSlowerThanGradientCalculator(seed, 11, factory);
    // 4 peaks
    gradientProcedureIsNotSlowerThanGradientCalculator(seed, 21, factory);
  }

  private void gradientProcedureIsNotSlowerThanGradientCalculator(RandomSeed seed,
      final int nparams, final BaseLsqLvmGradientProcedureFactory factory) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final int iter = 1000;
    final double[][] alpha = new double[nparams][nparams];
    final double[] beta = new double[nparams];

    final ArrayList<double[]> paramsList = new ArrayList<>(iter);
    final ArrayList<double[]> yList = new ArrayList<>(iter);

    final int[] x = createFakeData(RngFactory.create(seed.get()), nparams, iter, paramsList, yList);
    final int n = x.length;
    final FakeGradientFunction func = new FakeGradientFunction(blockWidth, nparams);

    final GradientCalculator calc = GradientCalculatorUtils.newCalculator(nparams, false);

    for (int i = 0; i < paramsList.size(); i++) {
      calc.findLinearised(n, yList.get(i), paramsList.get(i), alpha, beta, func);
    }

    for (int i = 0; i < paramsList.size(); i++) {
      final BaseLsqLvmGradientProcedure p = factory.createProcedure(yList.get(i), func);
      p.gradient(paramsList.get(i));
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
          final BaseLsqLvmGradientProcedure p = factory.createProcedure(yList.get(i), func);
          for (int j = loops; j-- > 0;) {
            p.gradient(paramsList.get(k++ % iter));
          }
        }
      }
    };
    final long time2 = t2.getTime();

    logger.log(TestLogging.getTimingRecord("GradientCalculator " + nparams, time1,
        factory.getClass().getSimpleName(), time2));
  }

  @SeededTest
  void gradientProcedureUnrolledComputesSameAsGradientProcedure(RandomSeed seed) {
    // Test the method that will be used for the standard and unrolled versions
    // for all other 'gradient procedures'
    gradientProcedureUnrolledComputesSameAsGradientProcedure(seed, 4);
    gradientProcedureUnrolledComputesSameAsGradientProcedure(seed, 5);
    gradientProcedureUnrolledComputesSameAsGradientProcedure(seed, 6);
  }

  private void gradientProcedureUnrolledComputesSameAsGradientProcedure(RandomSeed seed,
      int nparams) {
    final int iter = 10;

    final ArrayList<double[]> paramsList = new ArrayList<>(iter);
    final ArrayList<double[]> yList = new ArrayList<>(iter);

    createFakeData(RngFactory.create(seed.get()), nparams, iter, paramsList, yList);
    final FakeGradientFunction func = new FakeGradientFunction(blockWidth, nparams);

    final String name = GradientCalculator.class.getSimpleName();
    // Create messages
    final IndexSupplier msgR = new IndexSupplier(1, name + "Result: Not same ", null);
    final IndexSupplier msgOb = new IndexSupplier(1, name + "Observations: Not same beta ", null);
    final IndexSupplier msgOal =
        new IndexSupplier(1, name + "Observations: Not same alpha linear ", null);
    final IndexSupplier msgOam =
        new IndexSupplier(1, name + "Observations: Not same alpha matrix ", null);

    for (int i = 0; i < paramsList.size(); i++) {
      final BaseLsqLvmGradientProcedure p1 =
          LsqLvmGradientProcedureUtils.create(yList.get(i), func);
      p1.gradient(paramsList.get(i));

      final BaseLsqLvmGradientProcedure p2 = new LsqLvmGradientProcedure(yList.get(i), func);
      p2.gradient(paramsList.get(i));

      // Exactly the same ...
      Assertions.assertEquals(p1.value, p2.value, msgR.set(0, i));

      Assertions.assertArrayEquals(p1.beta, p2.beta, msgOb.set(0, i));

      Assertions.assertArrayEquals(p1.getAlphaLinear(), p2.getAlphaLinear(), msgOal.set(0, i));

      final double[][] am1 = p1.getAlphaMatrix();
      final double[][] am2 = p2.getAlphaMatrix();
      Assertions.assertArrayEquals(am1, am2, msgOam.set(0, i));
    }
  }

  @SpeedTag
  @SeededTest
  void gradientProcedureIsFasterUnrolledThanGradientProcedureMatrix(RandomSeed seed) {
    gradientProcedure2IsFasterUnrolledThanGradientProcedure1(seed,
        LsqLvmGradientProcedureMatrixUtils::create, LsqLvmGradientProcedureUtils::create);
  }

  @SpeedTag
  @SeededTest
  void gradientProcedureLinearIsFasterUnrolledThanGradientProcedureMatrix(RandomSeed seed) {
    gradientProcedure2IsFasterUnrolledThanGradientProcedure1(seed,
        LsqLvmGradientProcedureMatrixUtils::create, LsqLvmGradientProcedureLinearUtils::create);
  }

  private void gradientProcedure2IsFasterUnrolledThanGradientProcedure1(RandomSeed seed,
      BaseLsqLvmGradientProcedureFactory factory1, BaseLsqLvmGradientProcedureFactory factory2) {
    // Assert the unrolled versions
    gradientProcedureLinearIsFasterThanGradientProcedureMatrix(seed, 4, factory1, factory2, true);
    gradientProcedureLinearIsFasterThanGradientProcedureMatrix(seed, 5, factory1, factory2, true);
    gradientProcedureLinearIsFasterThanGradientProcedureMatrix(seed, 6, factory1, factory2, true);
    gradientProcedureLinearIsFasterThanGradientProcedureMatrix(seed, 11, factory1, factory2, false);
    gradientProcedureLinearIsFasterThanGradientProcedureMatrix(seed, 21, factory1, factory2, false);
  }

  private void gradientProcedureLinearIsFasterThanGradientProcedureMatrix(RandomSeed seed,
      final int nparams, final BaseLsqLvmGradientProcedureFactory factory1,
      final BaseLsqLvmGradientProcedureFactory factory2, boolean doAssert) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final int iter = 100;

    final ArrayList<double[]> paramsList = new ArrayList<>(iter);
    final ArrayList<double[]> yList = new ArrayList<>(iter);

    createData(RngFactory.create(seed.get()), 1, iter, paramsList, yList);

    // Remove the timing of the function call by creating a dummy function
    final Gradient1Function func = new FakeGradientFunction(blockWidth, nparams);

    // Create messages
    final IndexSupplier msgA = new IndexSupplier(1, "A ", null);
    final IndexSupplier msgB = new IndexSupplier(1, "B ", null);
    for (int i = 0; i < paramsList.size(); i++) {
      final BaseLsqLvmGradientProcedure p1 = factory1.createProcedure(yList.get(i), func);
      p1.gradient(paramsList.get(i));
      p1.gradient(paramsList.get(i));

      final BaseLsqLvmGradientProcedure p2 = factory2.createProcedure(yList.get(i), func);
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
          final BaseLsqLvmGradientProcedure p = factory1.createProcedure(yList.get(i), func);
          for (int j = loops; j-- > 0;) {
            p.gradient(paramsList.get(k++ % iter));
          }
        }
      }
    };
    final long time1 = t1.getTime();

    final Timer t2 = new Timer(t1.loops) {
      @Override
      void run() {
        for (int i = 0, k = 0; i < paramsList.size(); i++) {
          final BaseLsqLvmGradientProcedure p2 = factory2.createProcedure(yList.get(i), func);
          for (int j = loops; j-- > 0;) {
            p2.gradient(paramsList.get(k++ % iter));
          }
        }
      }
    };
    final long time2 = t2.getTime();

    final LogRecord record =
        TestLogging.getTimingRecord(factory1.getClass().getSimpleName() + nparams, time1,
            factory2.getClass().getSimpleName(), time2);
    if (!doAssert && record.getLevel() == TestLevel.TEST_FAILURE) {
      record.setLevel(TestLevel.TEST_WARNING);
    }
    logger.log(record);
  }

  @SeededTest
  void gradientProcedureComputesGradient(RandomSeed seed) {
    gradientProcedureComputesGradient(seed,
        new SingleFreeCircularErfGaussian2DFunction(blockWidth, blockWidth));
  }

  private void gradientProcedureComputesGradient(RandomSeed seed, ErfGaussian2DFunction func) {
    final int nparams = func.getNumberOfGradients();
    final int[] indices = func.gradientIndices();

    final int iter = 100;

    final ArrayList<double[]> paramsList = new ArrayList<>(iter);
    final ArrayList<double[]> yList = new ArrayList<>(iter);

    createData(RngFactory.create(seed.get()), 1, iter, paramsList, yList, true);

    // for the gradients
    final double delta = 1e-4;
    final DoubleDoubleBiPredicate eq = Predicates.doublesAreClose(5e-2, 1e-16);
    final IndexSupplier msg = new IndexSupplier(2).setMessagePrefix("Gradient ");

    // Must compute most of the time
    final int failureLimit = AssertionErrorCounter.computeFailureLimit(iter, 0.1);
    final AssertionErrorCounter failCounter = new AssertionErrorCounter(failureLimit, nparams);

    for (int i = 0; i < paramsList.size(); i++) {
      msg.set(0, i);
      final double[] y = yList.get(i);
      final double[] a = paramsList.get(i);
      final double[] a2 = a.clone();
      final BaseLsqLvmGradientProcedure p = LsqLvmGradientProcedureUtils.create(y, func);
      p.gradient(a);
      // double s = p.ssx;
      final double[] beta = p.beta.clone();
      for (int j = 0; j < nparams; j++) {
        final int jj = j;
        msg.set(1, j);
        final int k = indices[j];
        final double d = Precision.representableDelta(a[k], (a[k] == 0) ? 1e-3 : a[k] * delta);
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
        // logger.fine(FormatSupplier.getSupplier("[%d,%d] %f (%s %f+/-%f) %f ?= %f", i, k, s,
        // func.getName(k), a[k], d, beta[j],
        // gradient);
        failCounter.run(j, () -> TestAssertions.assertTest(gradient, beta[jj], eq, msg));
      }
    }
  }

  @SeededTest
  void gradientProcedureComputesSameOutputWithBias(RandomSeed seed) {
    final ErfGaussian2DFunction func =
        new SingleFreeCircularErfGaussian2DFunction(blockWidth, blockWidth);
    final int nparams = func.getNumberOfGradients();

    final int iter = 100;

    final Level logLevel = TestLevel.TEST_DEBUG;
    final boolean debug = logger.isLoggable(logLevel);

    final ArrayList<double[]> paramsList = new ArrayList<>(iter);
    final ArrayList<double[]> yList = new ArrayList<>(iter);

    final ArrayList<double[]> alphaList = new ArrayList<>(iter);
    final ArrayList<double[]> betaList = new ArrayList<>(iter);
    final ArrayList<double[]> xList = new ArrayList<>(iter);

    // Manipulate the background
    final double defaultBackground = background;
    try {
      background = 1e-2;
      createData(RngFactory.create(seed.get()), 1, iter, paramsList, yList, true);

      final EjmlLinearSolver solver = new EjmlLinearSolver(1e-5, 1e-6);

      for (int i = 0; i < paramsList.size(); i++) {
        final double[] y = yList.get(i);
        final double[] a = paramsList.get(i);
        final BaseLsqLvmGradientProcedure p = LsqLvmGradientProcedureUtils.create(y, func);
        p.gradient(a);
        final double[] beta = p.beta;
        alphaList.add(p.getAlphaLinear());
        betaList.add(beta.clone());
        for (int j = 0; j < nparams; j++) {
          if (Math.abs(beta[j]) < 1e-6) {
            logger.log(TestLogging.getRecord(TestLevel.TEST_INFO, "[%d] Tiny beta %s %g", i,
                func.getGradientParameterName(j), beta[j]));
          }
        }
        // Solve
        if (!solver.solve(p.getAlphaMatrix(), beta)) {
          throw new AssertionError();
        }
        xList.add(beta);
        // System.out.println(Arrays.toString(beta));
      }

      // for (int b = 1; b < 1000; b *= 2)
      for (final double b : new double[] {-500, -100, -10, -1, -0.1, 0, 0.1, 1, 10, 100, 500}) {
        final Statistics[] rel = new Statistics[nparams];
        final Statistics[] abs = new Statistics[nparams];
        if (debug) {
          for (int i = 0; i < nparams; i++) {
            rel[i] = new Statistics();
            abs[i] = new Statistics();
          }
        }

        for (int i = 0; i < paramsList.size(); i++) {
          final double[] y = add(yList.get(i), b);
          final double[] a = paramsList.get(i).clone();
          a[0] += b;
          final BaseLsqLvmGradientProcedure p = LsqLvmGradientProcedureUtils.create(y, func);
          p.gradient(a);
          final double[] beta = p.beta;
          final double[] alpha2 = alphaList.get(i);
          final double[] beta2 = betaList.get(i);
          final double[] x2 = xList.get(i);

          Assertions.assertArrayEquals(beta2, beta, 1e-10, "Beta");
          Assertions.assertArrayEquals(alpha2, p.getAlphaLinear(), 1e-10, "Alpha");

          // Solve
          solver.solve(p.getAlphaMatrix(), beta);
          Assertions.assertArrayEquals(x2, beta, 1e-10, "X");

          if (debug) {
            for (int j = 0; j < nparams; j++) {
              rel[j].add(DoubleEquality.relativeError(x2[j], beta[j]));
              abs[j].add(Math.abs(x2[j] - beta[j]));
            }
          }
        }

        if (debug) {
          for (int i = 0; i < nparams; i++) {
            logger.log(TestLogging.getRecord(logLevel,
                "Bias = %.2f : %s : Rel %g +/- %g: Abs %g +/- %g", b,
                func.getGradientParameterName(i), rel[i].getMean(), rel[i].getStandardDeviation(),
                abs[i].getMean(), abs[i].getStandardDeviation()));
          }
        }
      }
    } finally {
      background = defaultBackground;
    }
  }

  private static double[] add(double[] data, double value) {
    data = data.clone();
    for (int i = 0; i < data.length; i++) {
      data[i] += value;
    }
    return data;
  }

  /**
   * Create random elliptical Gaussian data an returns the data plus an estimate of the parameters.
   * Only the chosen parameters are randomised and returned for a maximum of (background, amplitude,
   * angle, xpos, ypos, xwidth, ywidth }
   *
   * @param rng the random
   * @param npeaks the npeaks
   * @param params set on output
   * @param randomiseParams Set to true to randomise the params
   * @return the double[]
   */
  private double[] doubleCreateGaussianData(UniformRandomProvider rng, int npeaks, double[] params,
      boolean randomiseParams) {
    final int n = blockWidth * blockWidth;

    // Generate a 2D Gaussian
    final ErfGaussian2DFunction func =
        (ErfGaussian2DFunction) GaussianFunctionFactory.create2D(npeaks, blockWidth, blockWidth,
            GaussianFunctionFactory.FIT_ERF_FREE_CIRCLE, null);
    params[0] = random(rng, background);
    for (int i = 0, j = 1; i < npeaks; i++, j += 6) {
      params[j] = random(rng, signal);
      params[j + 2] = random(rng, xpos);
      params[j + 3] = random(rng, ypos);
      params[j + 4] = random(rng, xwidth);
      params[j + 5] = random(rng, ywidth);
    }

    final double[] y = new double[n];
    func.initialise(params);
    for (int i = 0; i < y.length; i++) {
      // Add random Poisson noise
      final double u = func.eval(i);
      y[i] = GdscSmlmTestUtils.createPoissonSampler(rng, u).sample();
    }

    if (randomiseParams) {
      params[0] = random(rng, params[0]);
      for (int i = 0, j = 1; i < npeaks; i++, j += 6) {
        params[j] = random(rng, params[j]);
        params[j + 2] = random(rng, params[j + 2]);
        params[j + 3] = random(rng, params[j + 3]);
        params[j + 4] = random(rng, params[j + 4]);
        params[j + 5] = random(rng, params[j + 5]);
      }
    }

    return y;
  }

  protected int[] createData(final UniformRandomProvider rng, int npeaks, int iter,
      ArrayList<double[]> paramsList, ArrayList<double[]> yList) {
    return createData(rng, npeaks, iter, paramsList, yList, true);
  }

  protected int[] createData(final UniformRandomProvider rng, int npeaks, int iter,
      ArrayList<double[]> paramsList, ArrayList<double[]> yList, boolean randomiseParams) {
    final int[] x = new int[blockWidth * blockWidth];
    for (int i = 0; i < x.length; i++) {
      x[i] = i;
    }
    for (int i = 0; i < iter; i++) {
      final double[] params = new double[1 + 6 * npeaks];
      final double[] y = doubleCreateGaussianData(rng, npeaks, params, randomiseParams);
      paramsList.add(params);
      yList.add(y);
    }
    return x;
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

  protected ArrayList<double[]> copyList(ArrayList<double[]> paramsList) {
    final ArrayList<double[]> params2List = new ArrayList<>(paramsList.size());
    for (int i = 0; i < paramsList.size(); i++) {
      params2List.add(paramsList.get(i).clone());
    }
    return params2List;
  }
}
