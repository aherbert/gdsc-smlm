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
import java.util.Arrays;
import java.util.logging.Logger;
import org.apache.commons.math3.util.Precision;
import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.sampling.distribution.SharedStateContinuousSampler;
import org.ejml.data.DenseMatrix64F;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Disabled;
import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.rng.RandomUtils;
import uk.ac.sussex.gdsc.core.utils.rng.SamplerUtils;
import uk.ac.sussex.gdsc.smlm.GdscSmlmTestUtils;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.LvmGradientProcedureUtils.Type;
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
import uk.ac.sussex.gdsc.test.api.TestAssertions;
import uk.ac.sussex.gdsc.test.api.TestHelper;
import uk.ac.sussex.gdsc.test.api.function.DoubleDoubleBiPredicate;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngUtils;
import uk.ac.sussex.gdsc.test.utils.AssertionErrorCounter;
import uk.ac.sussex.gdsc.test.utils.RandomSeed;
import uk.ac.sussex.gdsc.test.utils.TestLogUtils;
import uk.ac.sussex.gdsc.test.utils.functions.FunctionUtils;
import uk.ac.sussex.gdsc.test.utils.functions.IndexSupplier;

/**
 * Contains speed tests for the methods for calculating the Hessian and gradient vector for use in
 * the LVM algorithm.
 */
@SuppressWarnings({"javadoc"})
class LvmGradientProcedureTest {
  private static Logger logger;

  @BeforeAll
  public static void beforeAll() {
    logger = Logger.getLogger(LvmGradientProcedureTest.class.getName());
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

  int maxIter = 20000;
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

  private static double random(UniformRandomProvider rng, double value) {
    return value - value * 0.1 + rng.nextDouble() * 0.2;
  }

  @SeededTest
  void gradientProcedureFactoryCreatesOptimisedProcedures() {
    final DummyGradientFunction[] f = new DummyGradientFunction[7];
    for (int i = 1; i < f.length; i++) {
      f[i] = new DummyGradientFunction(i);
    }

    final LvmGradientProcedureUtils.Type mle = LvmGradientProcedureUtils.Type.MLE;
    final LvmGradientProcedureUtils.Type wlsq = LvmGradientProcedureUtils.Type.WLSQ;
    final LvmGradientProcedureUtils.Type lsq = LvmGradientProcedureUtils.Type.LSQ;
    final LvmGradientProcedureUtils.Type fmle = LvmGradientProcedureUtils.Type.FAST_LOG_MLE;

    final FastLog fl = getFastLog();

    //@formatter:off

    // Generic factory
    final double[] y0 = new double[1];
    final double[] y1 = new double[]{ 1 };
    Assertions.assertEquals(LvmGradientProcedureUtils.create(y0, f[6], lsq, fl).getClass(), LsqLvmGradientProcedure6.class);
    Assertions.assertEquals(LvmGradientProcedureUtils.create(y0, f[5], lsq, fl).getClass(), LsqLvmGradientProcedure5.class);
    Assertions.assertEquals(LvmGradientProcedureUtils.create(y0, f[4], lsq, fl).getClass(), LsqLvmGradientProcedure4.class);
    Assertions.assertEquals(LvmGradientProcedureUtils.create(y0, f[1], lsq, fl).getClass(), LsqLvmGradientProcedure.class);
    Assertions.assertEquals(LvmGradientProcedureUtils.create(y0, f[6], mle, fl).getClass(), MleLvmGradientProcedure6.class);
    Assertions.assertEquals(LvmGradientProcedureUtils.create(y0, f[5], mle, fl).getClass(), MleLvmGradientProcedure5.class);
    Assertions.assertEquals(LvmGradientProcedureUtils.create(y0, f[4], mle, fl).getClass(), MleLvmGradientProcedure4.class);
    Assertions.assertEquals(LvmGradientProcedureUtils.create(y0, f[1], mle, fl).getClass(), MleLvmGradientProcedure.class);
    Assertions.assertEquals(LvmGradientProcedureUtils.create(y1, f[6], mle, fl).getClass(), MleLvmGradientProcedureX6.class);
    Assertions.assertEquals(LvmGradientProcedureUtils.create(y1, f[5], mle, fl).getClass(), MleLvmGradientProcedureX5.class);
    Assertions.assertEquals(LvmGradientProcedureUtils.create(y1, f[4], mle, fl).getClass(), MleLvmGradientProcedureX4.class);
    Assertions.assertEquals(LvmGradientProcedureUtils.create(y1, f[1], mle, fl).getClass(), MleLvmGradientProcedureX.class);
    Assertions.assertEquals(LvmGradientProcedureUtils.create(y0, f[6], wlsq, fl).getClass(), WLsqLvmGradientProcedure6.class);
    Assertions.assertEquals(LvmGradientProcedureUtils.create(y0, f[5], wlsq, fl).getClass(), WLsqLvmGradientProcedure5.class);
    Assertions.assertEquals(LvmGradientProcedureUtils.create(y0, f[4], wlsq, fl).getClass(), WLsqLvmGradientProcedure4.class);
    Assertions.assertEquals(LvmGradientProcedureUtils.create(y0, f[1], wlsq, fl).getClass(), WLsqLvmGradientProcedure.class);
    Assertions.assertEquals(LvmGradientProcedureUtils.create(y0, f[6], fmle, fl).getClass(), FastLogMleLvmGradientProcedure6.class);
    Assertions.assertEquals(LvmGradientProcedureUtils.create(y0, f[5], fmle, fl).getClass(), FastLogMleLvmGradientProcedure5.class);
    Assertions.assertEquals(LvmGradientProcedureUtils.create(y0, f[4], fmle, fl).getClass(), FastLogMleLvmGradientProcedure4.class);
    Assertions.assertEquals(LvmGradientProcedureUtils.create(y0, f[1], fmle, fl).getClass(), FastLogMleLvmGradientProcedure.class);
    Assertions.assertEquals(LvmGradientProcedureUtils.create(y1, f[6], fmle, fl).getClass(), FastLogMleLvmGradientProcedureX6.class);
    Assertions.assertEquals(LvmGradientProcedureUtils.create(y1, f[5], fmle, fl).getClass(), FastLogMleLvmGradientProcedureX5.class);
    Assertions.assertEquals(LvmGradientProcedureUtils.create(y1, f[4], fmle, fl).getClass(), FastLogMleLvmGradientProcedureX4.class);
    Assertions.assertEquals(LvmGradientProcedureUtils.create(y1, f[1], fmle, fl).getClass(), FastLogMleLvmGradientProcedureX.class);

    // Dedicated factories
    Assertions.assertEquals(LsqLvmGradientProcedureUtils.create(y0, f[6]).getClass(), LsqLvmGradientProcedure6.class);
    Assertions.assertEquals(LsqLvmGradientProcedureUtils.create(y0, f[5]).getClass(), LsqLvmGradientProcedure5.class);
    Assertions.assertEquals(LsqLvmGradientProcedureUtils.create(y0, f[4]).getClass(), LsqLvmGradientProcedure4.class);
    Assertions.assertEquals(LsqLvmGradientProcedureUtils.create(y0, f[1]).getClass(), LsqLvmGradientProcedure.class);
    Assertions.assertEquals(MleLvmGradientProcedureUtils.create(y0, f[6]).getClass(), MleLvmGradientProcedure6.class);
    Assertions.assertEquals(MleLvmGradientProcedureUtils.create(y0, f[5]).getClass(), MleLvmGradientProcedure5.class);
    Assertions.assertEquals(MleLvmGradientProcedureUtils.create(y0, f[4]).getClass(), MleLvmGradientProcedure4.class);
    Assertions.assertEquals(MleLvmGradientProcedureUtils.create(y0, f[1]).getClass(), MleLvmGradientProcedure.class);
    Assertions.assertEquals(MleLvmGradientProcedureUtils.create(y1, f[6]).getClass(), MleLvmGradientProcedureX6.class);
    Assertions.assertEquals(MleLvmGradientProcedureUtils.create(y1, f[5]).getClass(), MleLvmGradientProcedureX5.class);
    Assertions.assertEquals(MleLvmGradientProcedureUtils.create(y1, f[4]).getClass(), MleLvmGradientProcedureX4.class);
    Assertions.assertEquals(MleLvmGradientProcedureUtils.create(y1, f[1]).getClass(), MleLvmGradientProcedureX.class);
    Assertions.assertEquals(WLsqLvmGradientProcedureUtils.create(y0, null, f[6]).getClass(), WLsqLvmGradientProcedure6.class);
    Assertions.assertEquals(WLsqLvmGradientProcedureUtils.create(y0, null, f[5]).getClass(), WLsqLvmGradientProcedure5.class);
    Assertions.assertEquals(WLsqLvmGradientProcedureUtils.create(y0, null, f[4]).getClass(), WLsqLvmGradientProcedure4.class);
    Assertions.assertEquals(WLsqLvmGradientProcedureUtils.create(y0, null, f[1]).getClass(), WLsqLvmGradientProcedure.class);
    Assertions.assertEquals(MleLvmGradientProcedureUtils.create(y0, f[6], fl).getClass(), FastLogMleLvmGradientProcedure6.class);
    Assertions.assertEquals(MleLvmGradientProcedureUtils.create(y0, f[5], fl).getClass(), FastLogMleLvmGradientProcedure5.class);
    Assertions.assertEquals(MleLvmGradientProcedureUtils.create(y0, f[4], fl).getClass(), FastLogMleLvmGradientProcedure4.class);
    Assertions.assertEquals(MleLvmGradientProcedureUtils.create(y0, f[1], fl).getClass(), FastLogMleLvmGradientProcedure.class);
    Assertions.assertEquals(MleLvmGradientProcedureUtils.create(y1, f[6], fl).getClass(), FastLogMleLvmGradientProcedureX6.class);
    Assertions.assertEquals(MleLvmGradientProcedureUtils.create(y1, f[5], fl).getClass(), FastLogMleLvmGradientProcedureX5.class);
    Assertions.assertEquals(MleLvmGradientProcedureUtils.create(y1, f[4], fl).getClass(), FastLogMleLvmGradientProcedureX4.class);
    Assertions.assertEquals(MleLvmGradientProcedureUtils.create(y1, f[1], fl).getClass(), FastLogMleLvmGradientProcedureX.class);

    //@formatter:on
  }

  @SeededTest
  void gradientProcedureLsqComputesSameAsGradientCalculator(RandomSeed seed) {
    gradientProcedureComputesSameAsGradientCalculator(seed, Type.LSQ);
  }

  @SeededTest
  void gradientProcedureMleComputesSameAsGradientCalculator(RandomSeed seed) {
    gradientProcedureComputesSameAsGradientCalculator(seed, Type.MLE);
  }

  @SeededTest
  void gradientProcedureFastLogMleComputesSameAsGradientCalculator(RandomSeed seed) {
    gradientProcedureComputesSameAsGradientCalculator(seed, Type.FAST_LOG_MLE, 1e-3);
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

  private void gradientProcedureComputesSameAsGradientCalculator(RandomSeed seed, int nparams,
      Type type, double error) {
    final int iter = 10;

    final double[][] alpha = new double[nparams][nparams];
    final double[] beta = new double[nparams];

    final ArrayList<double[]> paramsList = new ArrayList<>(iter);
    final ArrayList<double[]> yList = new ArrayList<>(iter);

    final int[] x = createFakeData(RngUtils.create(seed.get()), nparams, iter, paramsList, yList);
    final int n = x.length;
    final FakeGradientFunction func = new FakeGradientFunction(blockWidth, nparams);

    final boolean mle = type != Type.LSQ;
    final FastLog fastLog = (type == Type.FAST_LOG_MLE) ? getFastLog() : null;
    final GradientCalculator calc = GradientCalculatorUtils.newCalculator(nparams, mle);

    final String name = String.format("[%d] %b", nparams, mle);
    // Create messages
    final IndexSupplier msgR = new IndexSupplier(1, name + "Result: Not same ", null);
    final IndexSupplier msgOb = new IndexSupplier(1, name + "Observations: Not same beta ", null);
    final IndexSupplier msgOal =
        new IndexSupplier(1, name + "Observations: Not same alpha linear ", null);
    final IndexSupplier msgOam =
        new IndexSupplier(1, name + "Observations: Not same alpha matrix ", null);

    final DoubleDoubleBiPredicate predicate =
        (error == 0) ? TestHelper.doublesEqual() : TestHelper.doublesAreClose(error, 0);

    for (int i = 0; i < paramsList.size(); i++) {
      // Reference implementation
      final double s = calc.findLinearised(n, yList.get(i), paramsList.get(i), alpha, beta, func);
      // Procedure
      final LvmGradientProcedure p =
          LvmGradientProcedureUtils.create(yList.get(i), func, type, fastLog);
      p.gradient(paramsList.get(i));
      final double s2 = p.value;
      // Value may be different depending on log implementation
      msgR.set(0, i);
      TestAssertions.assertTest(s, s2, predicate, msgR);

      // Exactly the same ...
      Assertions.assertArrayEquals(p.beta, beta, msgOb.set(0, i));

      final double[] al = p.getAlphaLinear();
      Assertions.assertArrayEquals(al, new DenseMatrix64F(alpha).data, msgOal.set(0, i));

      final double[][] am = p.getAlphaMatrix();
      Assertions.assertArrayEquals(am, alpha, msgOam.set(0, i));
    }
  }

  @SeededTest
  void gradientProcedureLsqUnrolledComputesSameAsGradientProcedure(RandomSeed seed) {
    gradientProcedureUnrolledComputesSameAsGradientProcedure(seed, Type.LSQ, false);
  }

  @SeededTest
  void gradientProcedureMleUnrolledComputesSameAsGradientProcedure(RandomSeed seed) {
    gradientProcedureUnrolledComputesSameAsGradientProcedure(seed, Type.MLE, false);
  }

  @SeededTest
  void gradientProcedureFastLogMleUnrolledComputesSameAsGradientProcedure(RandomSeed seed) {
    gradientProcedureUnrolledComputesSameAsGradientProcedure(seed, Type.FAST_LOG_MLE, false);
  }

  @SeededTest
  void gradientProcedureWLsqUnrolledComputesSameAsGradientProcedure(RandomSeed seed) {
    gradientProcedureUnrolledComputesSameAsGradientProcedure(seed, Type.WLSQ, false);
  }

  @SeededTest
  public void
      gradientProcedureLsqUnrolledComputesSameAsGradientProcedureWithPrecomputed(RandomSeed seed) {
    gradientProcedureUnrolledComputesSameAsGradientProcedure(seed, Type.LSQ, true);
  }

  @SeededTest
  public void
      gradientProcedureMleUnrolledComputesSameAsGradientProcedureWithPrecomputed(RandomSeed seed) {
    gradientProcedureUnrolledComputesSameAsGradientProcedure(seed, Type.MLE, true);
  }

  @SeededTest
  void gradientProcedureFastLogMleUnrolledComputesSameAsGradientProcedureWithPrecomputed(
      RandomSeed seed) {
    gradientProcedureUnrolledComputesSameAsGradientProcedure(seed, Type.FAST_LOG_MLE, true);
  }

  @SeededTest
  public void
      gradientProcedureWLsqUnrolledComputesSameAsGradientProcedureWithPrecomputed(RandomSeed seed) {
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

    createFakeData(RngUtils.create(seed.get()), nparams, iter, paramsList, yList);
    Gradient1Function func = new FakeGradientFunction(blockWidth, nparams);

    if (precomputed) {
      final double[] b = SimpleArrayUtils.newArray(func.size(), 0.1, 1.3);
      func = OffsetGradient1Function.wrapGradient1Function(func, b);
    }

    final FastLog fastLog = type == Type.FAST_LOG_MLE ? getFastLog() : null;

    final String name = String.format("[%d] %b", nparams, type);
    // Create messages
    final IndexSupplier msgR = new IndexSupplier(1, name + "Result: Not same ", null);
    final IndexSupplier msgOb = new IndexSupplier(1, name + "Observations: Not same beta ", null);
    final IndexSupplier msgOal =
        new IndexSupplier(1, name + "Observations: Not same alpha linear ", null);
    final IndexSupplier msgOam =
        new IndexSupplier(1, name + "Observations: Not same alpha matrix ", null);

    for (int i = 0; i < paramsList.size(); i++) {
      final LvmGradientProcedure p1 = createProcedure(type, yList.get(i), func, fastLog);
      p1.gradient(paramsList.get(i));

      final LvmGradientProcedure p2 =
          LvmGradientProcedureUtils.create(yList.get(i), func, type, fastLog);
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

  private static LvmGradientProcedure createProcedure(Type type, double[] y, Gradient1Function func,
      FastLog fastLog) {
    switch (type) {
      case FAST_LOG_MLE:
        return new FastLogMleLvmGradientProcedure(y, func, fastLog);
      case MLE:
        return new MleLvmGradientProcedure(y, func);
      case WLSQ:
        return new WLsqLvmGradientProcedure(y, null, func);
      default:
        return new LsqLvmGradientProcedure(y, func);
    }
  }

  @SeededTest
  void gradientProcedureLsqComputesGradient(RandomSeed seed) {
    gradientProcedureComputesGradient(seed,
        new SingleFreeCircularErfGaussian2DFunction(blockWidth, blockWidth), Type.LSQ, false);
  }

  @SeededTest
  void gradientProcedureMleComputesGradient(RandomSeed seed) {
    gradientProcedureComputesGradient(seed,
        new SingleFreeCircularErfGaussian2DFunction(blockWidth, blockWidth), Type.MLE, false);
  }

  @SeededTest
  void gradientProcedureFastLogMleCannotComputeGradient(RandomSeed seed) {
    Assertions.assertThrows(AssertionError.class, () -> {
      gradientProcedureComputesGradient(seed,
          new SingleFreeCircularErfGaussian2DFunction(blockWidth, blockWidth), Type.FAST_LOG_MLE,
          false);
    });
  }

  @Disabled("This test now passes as the tolerance for computing the gradient has been lowered "
      + " so that the tests pass under a stress test using many different random seeds.")
  @SeededTest
  void gradientProcedureFastLogMleCannotComputeGradientWithHighPrecision(RandomSeed seed) {
    // Try different precision
    for (int n = FastLog.N; n < 23; n++) {
      try {
        // logger.fine(FunctionUtils.getSupplier("Precision n=%d", n);
        fastLog = new TurboLog2(n);
        gradientProcedureComputesGradient(seed,
            new SingleFreeCircularErfGaussian2DFunction(blockWidth, blockWidth), Type.FAST_LOG_MLE,
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
  void gradientProcedureWLsqComputesGradient(RandomSeed seed) {
    gradientProcedureComputesGradient(seed,
        new SingleFreeCircularErfGaussian2DFunction(blockWidth, blockWidth), Type.WLSQ, false);
  }

  @SeededTest
  void gradientProcedureLsqComputesGradientWithPrecomputed(RandomSeed seed) {
    gradientProcedureComputesGradient(seed,
        new SingleFreeCircularErfGaussian2DFunction(blockWidth, blockWidth), Type.LSQ, true);
  }

  @SeededTest
  void gradientProcedureMleComputesGradientWithPrecomputed(RandomSeed seed) {
    gradientProcedureComputesGradient(seed,
        new SingleFreeCircularErfGaussian2DFunction(blockWidth, blockWidth), Type.MLE, true);
  }

  @SeededTest
  void gradientProcedureFastLogMleCannotComputeGradientWithPrecomputed(RandomSeed seed) {
    Assertions.assertThrows(AssertionError.class, () -> {
      gradientProcedureComputesGradient(seed,
          new SingleFreeCircularErfGaussian2DFunction(blockWidth, blockWidth), Type.FAST_LOG_MLE,
          true);
    });
  }

  @SeededTest
  void gradientProcedureWLsqComputesGradientWithPrecomputed(RandomSeed seed) {
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

    createData(RngUtils.create(seed.get()), 1, iter, paramsList, yList, true);

    // for the gradients
    final double delta = 1e-4;
    final DoubleDoubleBiPredicate eq = TestHelper.doublesAreClose(5e-2, 1e-16);
    final IndexSupplier msg = new IndexSupplier(2).setMessagePrefix("Gradient ");

    final double[] b = (precomputed) ? new double[func.size()] : null;

    final FastLog fastLog = type == Type.FAST_LOG_MLE ? getFastLog() : null;

    // Must compute most of the time
    final int failureLimit = AssertionErrorCounter.computeFailureLimit(iter, 0.1);
    final AssertionErrorCounter failCounter = new AssertionErrorCounter(failureLimit, nparams);

    for (int i = 0; i < paramsList.size(); i++) {
      msg.set(0, i);
      final double[] y = yList.get(i);
      final double[] a = paramsList.get(i);
      final double[] a2 = a.clone();

      LvmGradientProcedure gp;
      if (precomputed) {
        // Mock fitting part of the function already
        for (int j = 0; j < b.length; j++) {
          b[j] = y[j] * 0.5;
        }
        gp = LvmGradientProcedureUtils.create(y,
            OffsetGradient1Function.wrapGradient1Function(func, b), type, fastLog);
      } else {
        gp = LvmGradientProcedureUtils.create(y, func, type, fastLog);
      }
      gp.gradient(a);
      // double s = p.value;
      final double[] beta = gp.beta.clone();
      for (int j = 0; j < nparams; j++) {
        final int jj = j;
        msg.set(1, j);
        final int k = indices[j];
        // double d = Precision.representableDelta(a[k], (a[k] == 0) ? 1e-3 : a[k] * delta);
        final double d = Precision.representableDelta(a[k], delta);
        a2[k] = a[k] + d;
        gp.value(a2);
        final double s1 = gp.value;
        a2[k] = a[k] - d;
        gp.value(a2);
        final double s2 = gp.value;
        a2[k] = a[k];

        // Apply a factor of -2 to compute the actual gradients:
        // See Numerical Recipes in C++, 2nd Ed. Equation 15.5.6 for Nonlinear Models
        beta[j] *= -2;

        final double gradient = (s1 - s2) / (2 * d);
        // logger.fine(FunctionUtils.getSupplier("[%d,%d] %f (%s %f+/-%f) %f ?= %f", i, k, s,
        // Gaussian2DFunction.getName(k),
        // a[k], d, beta[j], gradient);
        failCounter.run(j, () -> TestAssertions.assertTest(gradient, beta[jj], eq, msg));
      }
    }
  }

  @SeededTest
  void gradientProcedureLsqSupportsPrecomputed(RandomSeed seed) {
    gradientProcedureSupportsPrecomputed(seed, Type.LSQ);
  }

  @SeededTest
  void gradientProcedureMleSupportsPrecomputed(RandomSeed seed) {
    gradientProcedureSupportsPrecomputed(seed, Type.MLE);
  }

  @SeededTest
  void gradientProcedureFastLogMleSupportsPrecomputed(RandomSeed seed) {
    gradientProcedureSupportsPrecomputed(seed, Type.FAST_LOG_MLE, false);
  }

  @SeededTest
  void gradientProcedureFastLogMleCannotSupportPrecomputedWithGradients(RandomSeed seed) {
    Assertions.assertThrows(AssertionError.class, () -> {
      gradientProcedureSupportsPrecomputed(seed, Type.FAST_LOG_MLE);
    });
  }

  @SeededTest
  void gradientProcedureWLsqSupportsPrecomputed(RandomSeed seed) {
    gradientProcedureSupportsPrecomputed(seed, Type.WLSQ);
  }

  private void gradientProcedureSupportsPrecomputed(RandomSeed seed, final Type type) {
    gradientProcedureSupportsPrecomputed(seed, type, true);
  }

  private void gradientProcedureSupportsPrecomputed(RandomSeed seed, final Type type,
      boolean checkGradients) {
    final int iter = 10;
    final UniformRandomProvider rng = RngUtils.create(seed.get());
    final SharedStateContinuousSampler gs = SamplerUtils.createGaussianSampler(rng, 0, noise);

    final ArrayList<double[]> paramsList = new ArrayList<>(iter);
    final ArrayList<double[]> yList = new ArrayList<>(iter);

    // 3 peaks
    createData(rng, 3, iter, paramsList, yList, true);

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

    final FastLog fastLog = type == Type.FAST_LOG_MLE ? getFastLog() : null;

    final int nparams = f12.getNumberOfGradients();
    final int[] indices = f12.gradientIndices();
    final double[] b = new double[f12.size()];

    final DoubleEquality eq = new DoubleEquality(1e-8, 1e-16); // for checking strict equivalence

    // for the gradients
    final double delta = 1e-4;
    final DoubleDoubleBiPredicate eq2 = TestHelper.doublesAreClose(5e-2, 1e-16);
    final IndexSupplier msg = new IndexSupplier(2).setMessagePrefix("Gradient ");

    final double[] a1peaks = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
    final double[] y_b = new double[b.length];

    // Count the number of failures for each gradient
    final int failureLimit = AssertionErrorCounter.computeFailureLimit(iter, 0.1);
    final AssertionErrorCounter failCounter = new AssertionErrorCounter(failureLimit, nparams * 2);

    for (int i = 0; i < paramsList.size(); i++) {
      msg.set(0, i);
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
        int index = 0;

        @Override
        public void execute(double value) {
          b[index] = value;
          // Remove negatives for MLE
          if (type.isMle()) {
            y[index] = Math.max(0, y[index]);
            y_b[index] = Math.max(0, y[index] - value);
          } else {
            y_b[index] = y[index] - value;
          }
          index++;
        }
      });

      final LvmGradientProcedure p123 = LvmGradientProcedureUtils.create(y, f123, type, fastLog);

      /////////////////////////////////////
      // These should be the same
      /////////////////////////////////////
      final LvmGradientProcedure p12b3 = LvmGradientProcedureUtils.create(y,
          OffsetGradient1Function.wrapGradient1Function(f12, b), type, fastLog);

      // Check they are the same
      p123.gradient(a3peaks);
      final double[][] m123 = p123.getAlphaMatrix();

      p12b3.gradient(a2peaks);
      double value = p12b3.value;
      final double[] beta = p12b3.beta.clone();
      double[][] alpha = p12b3.getAlphaMatrix();

      // logger.fine(FunctionUtils.getSupplier("MLE=%b [%d] p12b3 %f %f", type.isMLE(), i,
      // p123.value, s);

      if (!eq.almostEqualRelativeOrAbsolute(p123.value, value)) {
        Assertions.fail(FunctionUtils.getSupplier("p12b3 Not same value @ %d (error=%s) : %s == %s",
            i, DoubleEquality.relativeError(p123.value, value), p123.value, value));
      }
      if (!almostEqualRelativeOrAbsolute(eq, beta, p123.beta)) {
        Assertions
            .fail(FunctionUtils.getSupplier("p12b3 Not same gradient @ %d (error=%s) : %s vs %s", i,
                relativeError(beta, p123.beta), Arrays.toString(beta), Arrays.toString(p123.beta)));
      }
      for (int j = 0; j < alpha.length; j++) {
        // logger.fine(FunctionUtils.getSupplier("%s !=\n%s", Arrays.toString(alpha[j]),
        // Arrays.toString(m123[j]));
        if (!almostEqualRelativeOrAbsolute(eq, alpha[j], m123[j])) {
          Assertions
              .fail(FunctionUtils.getSupplier("p12b3 Not same alpha @ %d,%d (error=%s) : %s vs %s",
                  i, j, relativeError(alpha[j], m123[j]), Arrays.toString(alpha[j]),
                  Arrays.toString(m123[j])));
        }
      }

      // Check actual gradients are correct
      if (checkGradients) {
        for (int j = 0; j < nparams; j++) {
          final int jj = j;
          msg.set(1, j);
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
          failCounter.run(j, () -> TestAssertions.assertTest(gradient, beta[jj], eq2, msg));
        }
      }

      /////////////////////////////////////
      // This may be different
      /////////////////////////////////////
      final LvmGradientProcedure p12m3 = LvmGradientProcedureUtils.create(y_b, f12, type, fastLog);

      // Check these may be different.
      // Sometimes they are not different.

      p12m3.gradient(a2peaks);
      value = p12m3.value;
      System.arraycopy(p12m3.beta, 0, beta, 0, p12m3.beta.length);
      alpha = p12m3.getAlphaMatrix();

      // logger.fine(FunctionUtils.getSupplier("%s [%d] p12m3 %f %f", type, i, p123.value, s);

      // The test for different or equal is not robust to different random seeds.
      // ExtraAssertions.fail has been changed for TestLog.logFailure

      if (type != Type.LSQ) {
        if (eq.almostEqualRelativeOrAbsolute(p123.value, value)) {
          logger.log(TestLogUtils.getFailRecord("p12b3 Same value @ %d (error=%s) : %s == %s", i,
              DoubleEquality.relativeError(p123.value, value), p123.value, value));
        }
        if (almostEqualRelativeOrAbsolute(eq, beta, p123.beta)) {
          logger.log(TestLogUtils.getFailRecord("p12b3 Same gradient @ %d (error=%s) : %s vs %s", i,
              relativeError(beta, p123.beta), Arrays.toString(beta), Arrays.toString(p123.beta)));
        }

        // Note: Test the matrix is different by finding 1 different column
        int dj = -1;
        for (int j = 0; j < alpha.length; j++) {
          // logger.fine(FunctionUtils.getSupplier("%s !=\n%s\n", Arrays.toString(alpha[j]),
          // Arrays.toString(m123[j]));
          if (!almostEqualRelativeOrAbsolute(eq, alpha[j], m123[j])) {
            dj = j; // Different column
            break;
          }
        }
        if (dj == -1) {
          // Find biggest error for reporting. This helps set the test tolerance.
          double error = 0;
          dj = -1;
          for (int j = 0; j < alpha.length; j++) {
            final double e = relativeError(alpha[j], m123[j]);
            if (error <= e) {
              error = e;
              dj = j;
            }
          }
          logger.log(TestLogUtils.getFailRecord("p12b3 Same alpha @ %d,%d (error=%s) : %s vs %s", i,
              dj, error, Arrays.toString(alpha[dj]), Arrays.toString(m123[dj])));
        }
      } else {
        if (!eq.almostEqualRelativeOrAbsolute(p123.value, value)) {
          logger.log(TestLogUtils.getFailRecord("p12b3 Not same value @ %d (error=%s) : %s == %s",
              i, DoubleEquality.relativeError(p123.value, value), p123.value, value));
        }
        if (!almostEqualRelativeOrAbsolute(eq, beta, p123.beta)) {
          logger.log(TestLogUtils.getFailRecord(
              "p12b3 Not same gradient @ %d (error=%s) : %s vs %s", i,
              relativeError(beta, p123.beta), Arrays.toString(beta), Arrays.toString(p123.beta)));
        }
        for (int j = 0; j < alpha.length; j++) {
          // logger.fine(FunctionUtils.getSupplier("%s !=\n%s\n", Arrays.toString(alpha[j]),
          // Arrays.toString(m123[j]));
          if (!almostEqualRelativeOrAbsolute(eq, alpha[j], m123[j])) {
            logger.log(
                TestLogUtils.getFailRecord("p12b3 Not same alpha @ %d,%d (error=%s) : %s vs %s", i,
                    j, relativeError(alpha[j], m123[j]), Arrays.toString(alpha[j]),
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
        failCounter.run(nparams + j, () -> TestAssertions.assertTest(gradient, beta[jj], eq2, msg));
      }
    }
  }

  private static boolean almostEqualRelativeOrAbsolute(DoubleEquality eq, double[] a, double[] b) {
    for (int i = 0; i < a.length; i++) {
      if (!eq.almostEqualRelativeOrAbsolute(a[i], b[i])) {
        return false;
      }
    }
    return true;
  }

  private static double relativeError(double[] a, double[] b) {
    double error = 0;
    for (int i = 0; i < a.length; i++) {
      error = Math.max(error, DoubleEquality.relativeError(a[i], b[i]));
    }
    return error;
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
  private double[] doubleCreateGaussianData(UniformRandomProvider rng, int npeaks, double[] params,
      boolean randomiseParams) {
    final int n = blockWidth * blockWidth;

    // Generate a 2D Gaussian
    final ErfGaussian2DFunction func =
        (ErfGaussian2DFunction) GaussianFunctionFactory.create2D(npeaks, blockWidth, blockWidth,
            GaussianFunctionFactory.FIT_ERF_FREE_CIRCLE, null);
    params[0] = random(rng, background);
    for (int i = 0, j = 0; i < npeaks; i++, j += Gaussian2DFunction.PARAMETERS_PER_PEAK) {
      params[j + Gaussian2DFunction.SIGNAL] = random(rng, signal);
      params[j + Gaussian2DFunction.X_POSITION] = random(rng, xpos);
      params[j + Gaussian2DFunction.Y_POSITION] = random(rng, ypos);
      params[j + Gaussian2DFunction.X_SD] = random(rng, xwidth);
      params[j + Gaussian2DFunction.Y_SD] = random(rng, ywidth);
    }

    if (npeaks > 1) {
      // Move the peaks around so they do not overlap
      final double[] shift = SimpleArrayUtils.newArray(npeaks, -2, 4.0 / (npeaks - 1));
      RandomUtils.shuffle(shift, rng);
      for (int i = 0, j = 0; i < npeaks; i++, j += Gaussian2DFunction.PARAMETERS_PER_PEAK) {
        params[j + Gaussian2DFunction.X_POSITION] += shift[i];
      }
      RandomUtils.shuffle(shift, rng);
      for (int i = 0, j = 0; i < npeaks; i++, j += Gaussian2DFunction.PARAMETERS_PER_PEAK) {
        params[j + Gaussian2DFunction.Y_POSITION] += shift[i];
      }
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
      for (int i = 0, j = 0; i < npeaks; i++, j += Gaussian2DFunction.PARAMETERS_PER_PEAK) {
        params[j + Gaussian2DFunction.SIGNAL] = random(rng, params[j + Gaussian2DFunction.SIGNAL]);
        params[j + Gaussian2DFunction.X_POSITION] =
            random(rng, params[j + Gaussian2DFunction.X_POSITION]);
        params[j + Gaussian2DFunction.Y_POSITION] =
            random(rng, params[j + Gaussian2DFunction.Y_POSITION]);
        params[j + Gaussian2DFunction.X_SD] = random(rng, params[j + Gaussian2DFunction.X_SD]);
        params[j + Gaussian2DFunction.Y_SD] = random(rng, params[j + Gaussian2DFunction.Y_SD]);
      }
    }

    return y;
  }

  protected int[] createData(UniformRandomProvider rng, int npeaks, int iter,
      ArrayList<double[]> paramsList, ArrayList<double[]> yList) {
    return createData(rng, npeaks, iter, paramsList, yList, true);
  }

  protected int[] createData(UniformRandomProvider rng, int npeaks, int iter,
      ArrayList<double[]> paramsList, ArrayList<double[]> yList, boolean randomiseParams) {
    final int[] x = new int[blockWidth * blockWidth];
    for (int i = 0; i < x.length; i++) {
      x[i] = i;
    }
    for (int i = 0; i < iter; i++) {
      final double[] params = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK * npeaks];
      final double[] y = doubleCreateGaussianData(rng, npeaks, params, randomiseParams);
      paramsList.add(params);
      yList.add(y);
    }
    return x;
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

  private double[] createFakeData(UniformRandomProvider rng, double[] params) {
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
