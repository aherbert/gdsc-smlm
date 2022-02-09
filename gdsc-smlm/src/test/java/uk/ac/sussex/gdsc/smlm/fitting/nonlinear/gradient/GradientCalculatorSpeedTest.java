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
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.util.Precision;
import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.sampling.distribution.PoissonSampler;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.BeforeAll;
import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.Statistics;
import uk.ac.sussex.gdsc.smlm.GdscSmlmTestUtils;
import uk.ac.sussex.gdsc.smlm.fitting.linear.EjmlLinearSolver;
import uk.ac.sussex.gdsc.smlm.function.CameraNoiseModel;
import uk.ac.sussex.gdsc.smlm.function.NonLinearFunction;
import uk.ac.sussex.gdsc.smlm.function.PoissonCalculator;
import uk.ac.sussex.gdsc.smlm.function.gaussian.EllipticalGaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.SingleCircularGaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.SingleEllipticalGaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.SingleFixedGaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.SingleFreeCircularGaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.SingleNbFixedGaussian2DFunction;
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
import uk.ac.sussex.gdsc.test.utils.functions.FunctionUtils;
import uk.ac.sussex.gdsc.test.utils.functions.IntArrayFormatSupplier;

/**
 * Contains speed tests for the fastest method for calculating the Hessian and gradient vector for
 * use in non-linear fitting.
 */
@SuppressWarnings({"javadoc"})
class GradientCalculatorSpeedTest {
  private static Logger logger;

  @BeforeAll
  public static void beforeAll() {
    logger = Logger.getLogger(GradientCalculatorSpeedTest.class.getName());
  }

  @AfterAll
  public static void afterAll() {
    logger = null;
  }

  DoubleEquality eq = new DoubleEquality(1e-6, 1e-16);

  static final int maxIter = 20000;
  int blockWidth = 10;
  double background = 0.5;
  double amplitude = 100;
  double angle = Math.PI;
  double xpos = 5;
  double ypos = 5;
  double xwidth = 1.2;
  double ywidth = 1.2;

  private static double random(UniformRandomProvider rng, double value) {
    return value - value * 0.1 + rng.nextDouble() * 0.2;
  }

  @SeededTest
  void gradientCalculatorFactoryCreatesOptimisedCalculators() {
    Assertions.assertEquals(GradientCalculatorUtils.newCalculator(3).getClass(),
        GradientCalculator3.class);
    Assertions.assertEquals(GradientCalculatorUtils.newCalculator(4).getClass(),
        GradientCalculator4.class);
    Assertions.assertEquals(GradientCalculatorUtils.newCalculator(5).getClass(),
        GradientCalculator5.class);
    Assertions.assertEquals(GradientCalculatorUtils.newCalculator(6).getClass(),
        GradientCalculator6.class);
    Assertions.assertEquals(GradientCalculatorUtils.newCalculator(7).getClass(),
        GradientCalculator7.class);
    Assertions.assertEquals(GradientCalculatorUtils.newCalculator(13).getClass(),
        GradientCalculator.class);

    Assertions.assertEquals(GradientCalculatorUtils.newCalculator(3, true).getClass(),
        MleGradientCalculator3.class);
    Assertions.assertEquals(GradientCalculatorUtils.newCalculator(4, true).getClass(),
        MleGradientCalculator4.class);
    Assertions.assertEquals(GradientCalculatorUtils.newCalculator(5, true).getClass(),
        MleGradientCalculator5.class);
    Assertions.assertEquals(GradientCalculatorUtils.newCalculator(6, true).getClass(),
        MleGradientCalculator6.class);
    Assertions.assertEquals(GradientCalculatorUtils.newCalculator(7, true).getClass(),
        MleGradientCalculator7.class);
    Assertions.assertEquals(GradientCalculatorUtils.newCalculator(13, true).getClass(),
        MleGradientCalculator.class);
  }

  @SeededTest
  void mleGradientCalculator7ComputesSameAsGradientCalculator(RandomSeed seed) {
    gradientCalculatorNComputesSameAsGradientCalculator(seed,
        new SingleEllipticalGaussian2DFunction(blockWidth, blockWidth), 7, true);
  }

  @SpeedTag
  @SeededTest
  void mleGradientCalculator7IsFasterThanGradientCalculator(RandomSeed seed) {
    gradientCalculatorNIsFasterThanGradientCalculator(seed,
        new SingleEllipticalGaussian2DFunction(blockWidth, blockWidth), 7, true);
  }

  @SeededTest
  void mleGradientCalculator6ComputesSameAsGradientCalculator(RandomSeed seed) {
    gradientCalculatorNComputesSameAsGradientCalculator(seed,
        new SingleFreeCircularGaussian2DFunction(blockWidth, blockWidth), 6, true);
  }

  @SpeedTag
  @SeededTest
  void mleGradientCalculator6IsFasterThanGradientCalculator(RandomSeed seed) {
    gradientCalculatorNIsFasterThanGradientCalculator(seed,
        new SingleFreeCircularGaussian2DFunction(blockWidth, blockWidth), 6, true);
  }

  @SeededTest
  void mleGradientCalculator5ComputesSameAsGradientCalculator(RandomSeed seed) {
    gradientCalculatorNComputesSameAsGradientCalculator(seed,
        new SingleCircularGaussian2DFunction(blockWidth, blockWidth), 5, true);
  }

  @SpeedTag
  @SeededTest
  void mleGradientCalculator5IsFasterThanGradientCalculator(RandomSeed seed) {
    gradientCalculatorNIsFasterThanGradientCalculator(seed,
        new SingleCircularGaussian2DFunction(blockWidth, blockWidth), 5, true);
  }

  @SeededTest
  void mleGradientCalculator4ComputesSameAsGradientCalculator(RandomSeed seed) {
    gradientCalculatorNComputesSameAsGradientCalculator(seed,
        new SingleFixedGaussian2DFunction(blockWidth, blockWidth), 4, true);
  }

  @SpeedTag
  @SeededTest
  void mleGradientCalculator4IsFasterThanGradientCalculator(RandomSeed seed) {
    gradientCalculatorNIsFasterThanGradientCalculator(seed,
        new SingleFixedGaussian2DFunction(blockWidth, blockWidth), 4, true);
  }

  @SeededTest
  void mleGradientCalculator3ComputesSameAsGradientCalculator(RandomSeed seed) {
    gradientCalculatorNComputesSameAsGradientCalculator(seed,
        new SingleNbFixedGaussian2DFunction(blockWidth, blockWidth), 3, true);
  }

  @SpeedTag
  @SeededTest
  void mleGradientCalculator3IsFasterThanGradientCalculator(RandomSeed seed) {
    gradientCalculatorNIsFasterThanGradientCalculator(seed,
        new SingleNbFixedGaussian2DFunction(blockWidth, blockWidth), 3, true);
  }

  @SeededTest
  void gradientCalculator7ComputesSameAsGradientCalculator(RandomSeed seed) {
    gradientCalculatorNComputesSameAsGradientCalculator(seed,
        new SingleEllipticalGaussian2DFunction(blockWidth, blockWidth), 7, false);
  }

  @SpeedTag
  @SeededTest
  void gradientCalculator7IsFasterThanGradientCalculator(RandomSeed seed) {
    gradientCalculatorNIsFasterThanGradientCalculator(seed,
        new SingleEllipticalGaussian2DFunction(blockWidth, blockWidth), 7, false);
  }

  @SeededTest
  void gradientCalculator6ComputesSameAsGradientCalculator(RandomSeed seed) {
    gradientCalculatorNComputesSameAsGradientCalculator(seed,
        new SingleFreeCircularGaussian2DFunction(blockWidth, blockWidth), 6, false);
  }

  @SpeedTag
  @SeededTest
  void gradientCalculator6IsFasterThanGradientCalculator(RandomSeed seed) {
    gradientCalculatorNIsFasterThanGradientCalculator(seed,
        new SingleFreeCircularGaussian2DFunction(blockWidth, blockWidth), 6, false);
  }

  @SeededTest
  void gradientCalculator5ComputesSameAsGradientCalculator(RandomSeed seed) {
    gradientCalculatorNComputesSameAsGradientCalculator(seed,
        new SingleCircularGaussian2DFunction(blockWidth, blockWidth), 5, false);
  }

  @SpeedTag
  @SeededTest
  void gradientCalculator5IsFasterThanGradientCalculator(RandomSeed seed) {
    gradientCalculatorNIsFasterThanGradientCalculator(seed,
        new SingleCircularGaussian2DFunction(blockWidth, blockWidth), 5, false);
  }

  @SeededTest
  void gradientCalculator4ComputesSameAsGradientCalculator(RandomSeed seed) {
    gradientCalculatorNComputesSameAsGradientCalculator(seed,
        new SingleFixedGaussian2DFunction(blockWidth, blockWidth), 4, false);
  }

  @SpeedTag
  @SeededTest
  void gradientCalculator4IsFasterThanGradientCalculator(RandomSeed seed) {
    gradientCalculatorNIsFasterThanGradientCalculator(seed,
        new SingleFixedGaussian2DFunction(blockWidth, blockWidth), 4, false);
  }

  @SeededTest
  void gradientCalculator3ComputesSameAsGradientCalculator(RandomSeed seed) {
    gradientCalculatorNComputesSameAsGradientCalculator(seed,
        new SingleNbFixedGaussian2DFunction(blockWidth, blockWidth), 3, false);
  }

  @SpeedTag
  @SeededTest
  void gradientCalculator3IsFasterThanGradientCalculator(RandomSeed seed) {
    gradientCalculatorNIsFasterThanGradientCalculator(seed,
        new SingleNbFixedGaussian2DFunction(blockWidth, blockWidth), 3, false);
  }

  private void gradientCalculatorNComputesSameAsGradientCalculator(RandomSeed seed,
      Gaussian2DFunction func, int nparams, boolean mle) {
    // Check the function is the correct size
    Assertions.assertEquals(nparams, func.gradientIndices().length);

    final int iter = 50;

    final double[][] alpha = new double[nparams][nparams];
    final double[] beta = new double[nparams];
    final double[][] alpha2 = new double[nparams][nparams];
    final double[] beta2 = new double[nparams];

    final ArrayList<double[]> paramsList = new ArrayList<>(iter);
    final ArrayList<double[]> yList = new ArrayList<>(iter);

    final int[] x = createData(RngUtils.create(seed.getSeed()), 1, iter, paramsList, yList);

    final GradientCalculator calc =
        (mle) ? new MleGradientCalculator(beta.length) : new GradientCalculator(beta.length);
    final GradientCalculator calc2 = GradientCalculatorUtils.newCalculator(nparams, mle);

    // Create messages
    IntArrayFormatSupplier msgR = new IntArrayFormatSupplier("Result: Not same @ %d", 1);
    IntArrayFormatSupplier msgB =
        new IntArrayFormatSupplier("Observations: Not same beta @ %d [%d]", 2);
    IntArrayFormatSupplier msgA =
        new IntArrayFormatSupplier("Observations: Not same alpha @ %d [%d][%d]", 3);

    for (int i = 0; i < paramsList.size(); i++) {
      final double s = calc.findLinearised(x, yList.get(i), paramsList.get(i), alpha, beta, func);
      final double s2 =
          calc2.findLinearised(x, yList.get(i), paramsList.get(i), alpha2, beta2, func);
      Assertions.assertTrue(eq.almostEqualRelativeOrAbsolute(s, s2), msgR.set(0, i));
      msgB.set(0, i);
      for (int j = 0; j < beta.length; j++) {
        Assertions.assertTrue(eq.almostEqualRelativeOrAbsolute(beta[j], beta2[j]), msgB.set(1, j));
      }
      msgA.set(0, i);
      for (int j = 0; j < beta.length; j++) {
        msgA.set(1, j);
        for (int k = 0; k < alpha[j].length; k++) {
          Assertions.assertTrue(eq.almostEqualRelativeOrAbsolute(alpha[j][k], alpha2[j][k]),
              msgA.set(2, k));
        }
      }
    }

    msgR = new IntArrayFormatSupplier("N-Result: Not same @ %d", 1);
    msgB = new IntArrayFormatSupplier("N-Observations: Not same beta @ %d [%d]", 2);
    msgA = new IntArrayFormatSupplier("N-Observations: Not same alpha @ %d [%d][%d]", 3);

    for (int i = 0; i < paramsList.size(); i++) {
      final double s =
          calc.findLinearised(x.length, yList.get(i), paramsList.get(i), alpha, beta, func);
      final double s2 =
          calc2.findLinearised(x.length, yList.get(i), paramsList.get(i), alpha2, beta2, func);
      Assertions.assertTrue(eq.almostEqualRelativeOrAbsolute(s, s2), msgR.set(0, i));
      msgB.set(0, i);
      for (int j = 0; j < beta.length; j++) {
        Assertions.assertTrue(eq.almostEqualRelativeOrAbsolute(beta[j], beta2[j]), msgB.set(1, j));
      }
      msgA.set(0, i);
      for (int j = 0; j < beta.length; j++) {
        msgA.set(1, j);
        for (int k = 0; k < alpha[j].length; k++) {
          Assertions.assertTrue(eq.almostEqualRelativeOrAbsolute(alpha[j][k], alpha2[j][k]),
              msgA.set(2, k));
        }
      }
    }

    if (!mle) {
      func.setNoiseModel(CameraNoiseModel.createNoiseModel(10, 0, true));

      msgR = new IntArrayFormatSupplier("Result+Noise: Not same @ %d", 1);
      msgB = new IntArrayFormatSupplier("Observations+Noise: Not same beta @ %d [%d]", 2);
      msgA = new IntArrayFormatSupplier("Observations+Noise: Not same alpha @ %d [%d][%d]", 3);

      for (int i = 0; i < paramsList.size(); i++) {
        final double s = calc.findLinearised(x, yList.get(i), paramsList.get(i), alpha, beta, func);
        final double s2 =
            calc2.findLinearised(x, yList.get(i), paramsList.get(i), alpha2, beta2, func);
        Assertions.assertTrue(eq.almostEqualRelativeOrAbsolute(s, s2), msgR.set(0, i));
        msgB.set(0, i);
        for (int j = 0; j < beta.length; j++) {
          Assertions.assertTrue(eq.almostEqualRelativeOrAbsolute(beta[j], beta2[j]),
              msgB.set(1, j));
        }
        msgA.set(0, i);
        for (int j = 0; j < beta.length; j++) {
          msgA.set(1, j);
          for (int k = 0; k < alpha[j].length; k++) {
            Assertions.assertTrue(eq.almostEqualRelativeOrAbsolute(alpha[j][k], alpha2[j][k]),
                msgA.set(2, k));
          }
        }
      }

      msgR = new IntArrayFormatSupplier("N-Result+Noise: Not same @ %d", 1);
      msgB = new IntArrayFormatSupplier("N-Observations+Noise: Not same beta @ %d [%d]", 2);
      msgA = new IntArrayFormatSupplier("N-Observations+Noise: Not same alpha @ %d [%d][%d]", 3);

      for (int i = 0; i < paramsList.size(); i++) {
        final double s =
            calc.findLinearised(x.length, yList.get(i), paramsList.get(i), alpha, beta, func);
        final double s2 =
            calc2.findLinearised(x.length, yList.get(i), paramsList.get(i), alpha2, beta2, func);
        Assertions.assertTrue(eq.almostEqualRelativeOrAbsolute(s, s2), msgR.set(0, i));
        msgB.set(0, i);
        for (int j = 0; j < beta.length; j++) {
          Assertions.assertTrue(eq.almostEqualRelativeOrAbsolute(beta[j], beta2[j]),
              msgB.set(1, j));
        }
        msgA.set(0, i);
        for (int j = 0; j < beta.length; j++) {
          msgA.set(1, j);
          for (int k = 0; k < alpha[j].length; k++) {
            Assertions.assertTrue(eq.almostEqualRelativeOrAbsolute(alpha[j][k], alpha2[j][k]),
                msgA.set(2, k));
          }
        }
      }
    }
  }

  private void gradientCalculatorNIsFasterThanGradientCalculator(RandomSeed seed,
      Gaussian2DFunction func, int nparams, boolean mle) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    // Check the function is the correct size
    Assertions.assertEquals(nparams, func.gradientIndices().length);

    final int iter = 10000;
    final double[][] alpha = new double[nparams][nparams];
    final double[] beta = new double[nparams];

    final ArrayList<double[]> paramsList = new ArrayList<>(iter);
    final ArrayList<double[]> yList = new ArrayList<>(iter);

    final int[] x = createData(RngUtils.create(seed.getSeed()), 1, iter, paramsList, yList);

    final GradientCalculator calc =
        (mle) ? new MleGradientCalculator(beta.length) : new GradientCalculator(beta.length);
    final GradientCalculator calc2 = GradientCalculatorUtils.newCalculator(nparams, mle);

    for (int i = 0; i < paramsList.size(); i++) {
      calc.findLinearised(x, yList.get(i), paramsList.get(i), alpha, beta, func);
    }

    for (int i = 0; i < paramsList.size(); i++) {
      calc2.findLinearised(x, yList.get(i), paramsList.get(i), alpha, beta, func);
    }

    long start1 = System.nanoTime();
    for (int i = 0; i < paramsList.size(); i++) {
      calc.findLinearised(x, yList.get(i), paramsList.get(i), alpha, beta, func);
    }
    start1 = System.nanoTime() - start1;

    long start2 = System.nanoTime();
    for (int i = 0; i < paramsList.size(); i++) {
      calc2.findLinearised(x, yList.get(i), paramsList.get(i), alpha, beta, func);
    }
    start2 = System.nanoTime() - start2;

    logger.log(TestLogUtils.getTimingRecord(
        ((mle) ? "MLE " : "") + "Linearised GradientCalculator " + nparams, start1,
        "GradientCalculator", start2));
  }

  @SeededTest
  void gradientCalculatorAssumedXIsFasterThanGradientCalculator(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final int iter = 10000;

    final ArrayList<double[]> paramsList = new ArrayList<>(iter);
    final ArrayList<double[]> yList = new ArrayList<>(iter);

    final int[] x = createData(RngUtils.create(seed.getSeed()), 1, iter, paramsList, yList);

    final GradientCalculator calc = new GradientCalculator6();
    final GradientCalculator calc2 = new GradientCalculator6();
    final SingleFreeCircularGaussian2DFunction func =
        new SingleFreeCircularGaussian2DFunction(blockWidth, blockWidth);
    final int n = x.length;
    final int ng = func.getNumberOfGradients();
    final double[][] alpha = new double[ng][ng];
    final double[] beta = new double[ng];

    for (int i = 0; i < paramsList.size(); i++) {
      calc.findLinearised(x, yList.get(i), paramsList.get(i), alpha, beta, func);
    }

    for (int i = 0; i < paramsList.size(); i++) {
      calc2.findLinearised(n, yList.get(i), paramsList.get(i), alpha, beta, func);
    }

    long start1 = System.nanoTime();
    for (int i = 0; i < paramsList.size(); i++) {
      calc.findLinearised(x, yList.get(i), paramsList.get(i), alpha, beta, func);
    }
    start1 = System.nanoTime() - start1;

    long start2 = System.nanoTime();
    for (int i = 0; i < paramsList.size(); i++) {
      calc2.findLinearised(n, yList.get(i), paramsList.get(i), alpha, beta, func);
    }
    start2 = System.nanoTime() - start2;

    logger.log(TestLogUtils.getTimingRecord("GradientCalculator", start1,
        "GradientCalculatorAssumed", start2));
  }

  @SeededTest
  void lvmGradientCalculatorComputesGradient(RandomSeed seed) {
    gradientCalculatorComputesGradient(seed, new GradientCalculator(7));
  }

  @SeededTest
  void mleGradientCalculatorComputesGradient(RandomSeed seed) {
    gradientCalculatorComputesGradient(seed, new MleGradientCalculator(7));
  }

  private void gradientCalculatorComputesGradient(RandomSeed seed, GradientCalculator calc) {
    final int nparams = calc.nparams;
    final Gaussian2DFunction func = new SingleEllipticalGaussian2DFunction(blockWidth, blockWidth);
    // Check the function is the correct size
    final int[] indices = func.gradientIndices();
    Assertions.assertEquals(nparams, indices.length);

    final int iter = 50;

    final double[] beta = new double[nparams];
    final double[] beta2 = new double[nparams];

    final ArrayList<double[]> paramsList = new ArrayList<>(iter);
    final ArrayList<double[]> yList = new ArrayList<>(iter);

    final int[] x = createData(RngUtils.create(seed.getSeed()), 1, iter, paramsList, yList, true);

    final double delta = 1e-3;
    final DoubleEquality eq = new DoubleEquality(1e-3, 1e-3);

    final IntArrayFormatSupplier msg = new IntArrayFormatSupplier("[%d] Not same gradient @ %d", 2);

    for (int i = 0; i < paramsList.size(); i++) {
      msg.set(0, i);

      final double[] y = yList.get(i);
      final double[] a = paramsList.get(i);
      final double[] a2 = a.clone();
      // double s =
      calc.evaluate(x, y, a, beta, func);

      for (int k = 0; k < nparams; k++) {
        final int j = indices[k];
        final double d = Precision.representableDelta(a[j], (a[j] == 0) ? 1e-3 : a[j] * delta);
        a2[j] = a[j] + d;
        final double s1 = calc.evaluate(x, y, a2, beta2, func);
        a2[j] = a[j] - d;
        final double s2 = calc.evaluate(x, y, a2, beta2, func);
        a2[j] = a[j];

        final double gradient = (s1 - s2) / (2 * d);
        // logger.fine(FunctionUtils.getSupplier("[%d,%d] %f (%s %f+/-%f) %f ?= %f", i, j, s,
        // func.getName(j), a[j], d, beta[k],
        // gradient));
        Assertions.assertTrue(eq.almostEqualRelativeOrAbsolute(beta[k], gradient), msg.set(1, j));
      }
    }
  }

  @SeededTest
  void mleGradientCalculatorComputesLikelihood() {
    //@formatter:off
    final NonLinearFunction func = new NonLinearFunction(){
      double u;
      @Override
      public void initialise(double[] a) { u = a[0]; }
      @Override
      public int[] gradientIndices() { return null; }
      @Override
      public double eval(int x, double[] dyda)  { return 0; }
      @Override
      public double eval(int x) {
        return u;
      }
      @Override
      public double evalw(int x, double[] dyda, double[] w) { return 0; }
      @Override
      public double evalw(int x, double[] w) { return 0; }
      @Override
      public boolean canComputeWeights() { return false; }
      @Override
      public int getNumberOfGradients() {  return 0; }
    };
    //@formatter:on

    final DoubleDoubleBiPredicate predicate = TestHelper.doublesAreClose(1e-10, 0);

    final int[] xx = SimpleArrayUtils.natural(100);
    final double[] xxx = SimpleArrayUtils.newArray(100, 0, 1.0);
    for (final double u : new double[] {0.79, 2.5, 5.32}) {
      double ll = 0;
      double oll = 0;
      final PoissonDistribution pd = new PoissonDistribution(u);

      // The likelihood and logLikelihood computation are delegated to the PoissonCalculator.
      // However these are tested here to check they are correct for the observations.
      // The logLikelihood function for the entire array of observations is then asserted.
      for (final int x : xx) {
        double obs = PoissonCalculator.likelihood(u, x);
        double exp = pd.probability(x);
        TestAssertions.assertTest(exp, obs, predicate, "likelihood");

        obs = PoissonCalculator.logLikelihood(u, x);
        exp = pd.logProbability(x);
        TestAssertions.assertTest(exp, obs, predicate, "log likelihood");

        oll += obs;
        ll += exp;
      }

      final MleGradientCalculator gc = new MleGradientCalculator(1);
      final double o = gc.logLikelihood(xxx, new double[] {u}, func);

      Assertions.assertEquals(oll, o,
          "sum log likelihood should exactly match the PoissonCalculator");
      TestAssertions.assertTest(ll, o, predicate, "sum log likelihood");
    }
  }

  @SeededTest
  void gradientCalculatorComputesSameOutputWithBias(RandomSeed seed) {

    final Gaussian2DFunction func = new SingleEllipticalGaussian2DFunction(blockWidth, blockWidth);
    final int nparams = func.getNumberOfGradients();
    final GradientCalculator calc = new GradientCalculator(nparams);
    final int n = func.size();

    final int iter = 50;

    final ArrayList<double[]> paramsList = new ArrayList<>(iter);
    final ArrayList<double[]> yList = new ArrayList<>(iter);

    final ArrayList<double[][]> alphaList = new ArrayList<>(iter);
    final ArrayList<double[]> betaList = new ArrayList<>(iter);
    final ArrayList<double[]> xList = new ArrayList<>(iter);

    // Manipulate the background
    final double defaultBackground = background;
    final boolean report = logger.isLoggable(Level.INFO);
    try {
      background = 1e-2;
      createData(RngUtils.create(seed.getSeed()), 1, iter, paramsList, yList, true);

      final EjmlLinearSolver solver = new EjmlLinearSolver(1e-5, 1e-6);

      for (int i = 0; i < paramsList.size(); i++) {
        final double[] y = yList.get(i);
        final double[] a = paramsList.get(i);
        final double[][] alpha = new double[nparams][nparams];
        final double[] beta = new double[nparams];
        calc.findLinearised(n, y, a, alpha, beta, func);
        alphaList.add(alpha);
        betaList.add(beta.clone());
        for (int j = 0; j < nparams; j++) {
          if (Math.abs(beta[j]) < 1e-6) {
            logger.info(FunctionUtils.getSupplier("[%d] Tiny beta %s %g", i,
                func.getGradientParameterName(j), beta[j]));
          }
        }
        // Solve
        if (!solver.solve(alpha, beta)) {
          throw new AssertionError();
        }
        xList.add(beta);
        // System.out.println(Arrays.toString(beta));
      }

      final double[][] alpha = new double[nparams][nparams];
      final double[] beta = new double[nparams];
      final Statistics[] rel = new Statistics[nparams];
      final Statistics[] abs = new Statistics[nparams];
      for (int i = 0; i < nparams; i++) {
        rel[i] = new Statistics();
        abs[i] = new Statistics();
      }

      final DoubleDoubleBiPredicate predicate = TestHelper.doublesAreClose(1e-10, 0);

      // for (int b = 1; b < 1000; b *= 2)
      // for (double b : new double[] { -500, -100, -10, -1, -0.1, 0.1, 1, 10, 100, 500 })
      for (final double b : new double[] {-10, -1, -0.1, 0.1, 1, 10}) {
        if (report) {
          for (int i = 0; i < nparams; i++) {
            rel[i].reset();
            abs[i].reset();
          }
        }

        for (int i = 0; i < paramsList.size(); i++) {
          final double[] y = add(yList.get(i), b);
          final double[] a = paramsList.get(i).clone();
          a[0] += b;
          calc.findLinearised(n, y, a, alpha, beta, func);
          final double[][] alpha2 = alphaList.get(i);
          final double[] beta2 = betaList.get(i);
          final double[] x2 = xList.get(i);

          TestAssertions.assertArrayTest(beta2, beta, predicate, "Beta");
          TestAssertions.assertArrayTest(alpha2, alpha, predicate, "Alpha");

          // Solve
          solver.solve(alpha, beta);
          Assertions.assertArrayEquals(x2, beta, 1e-10, "X");

          if (report) {
            for (int j = 0; j < nparams; j++) {
              rel[j].add(DoubleEquality.relativeError(x2[j], beta[j]));
              abs[j].add(Math.abs(x2[j] - beta[j]));
            }
          }
        }

        if (report) {
          for (int i = 0; i < nparams; i++) {
            logger.info(FunctionUtils.getSupplier("Bias = %.2f : %s : Rel %g +/- %g: Abs %g +/- %g",
                b, func.getGradientParameterName(i), rel[i].getMean(),
                rel[i].getStandardDeviation(), abs[i].getMean(), abs[i].getStandardDeviation()));
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

  @SeededTest
  void mleCalculatorComputesLogLikelihoodRatio(RandomSeed seed) {
    final EllipticalGaussian2DFunction func =
        new EllipticalGaussian2DFunction(1, blockWidth, blockWidth);
    final int n = blockWidth * blockWidth;
    final double[] a = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
    final UniformRandomProvider rng = RngUtils.create(seed.getSeed());
    final DoubleDoubleBiPredicate predicate = TestHelper.doublesAreClose(1e-10, 0);
    for (int run = 5; run-- > 0;) {
      a[Gaussian2DFunction.BACKGROUND] = random(rng, background);
      a[Gaussian2DFunction.SIGNAL] = random(rng, amplitude);
      a[Gaussian2DFunction.ANGLE] = random(rng, angle);
      a[Gaussian2DFunction.X_POSITION] = random(rng, xpos);
      a[Gaussian2DFunction.Y_POSITION] = random(rng, ypos);
      a[Gaussian2DFunction.X_SD] = random(rng, xwidth);
      a[Gaussian2DFunction.Y_SD] = random(rng, ywidth);

      // Simulate Poisson process
      func.initialise(a);
      final double[] u = new double[n];
      final double[] x = new double[n];
      for (int i = 0; i < n; i++) {
        final double value = func.eval(i);
        u[i] = value;
        // Add random Poisson noise
        if (value > 0) {
          x[i] = new PoissonSampler(rng, value).sample();
        }
      }

      final int ng = func.getNumberOfGradients();
      final double[][] alpha = new double[ng][ng];
      final double[] beta = new double[ng];

      final GradientCalculator calc = GradientCalculatorUtils.newCalculator(ng, true);

      final double llr = PoissonCalculator.logLikelihoodRatio(u, x);
      final double llr2 = calc.findLinearised(n, x, a, alpha, beta, func);
      // logger.fine(FunctionUtils.getSupplier("llr=%f, llr2=%f", llr, llr2));
      TestAssertions.assertTest(llr, llr2, predicate, "Log-likelihood ratio");
    }
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
    final EllipticalGaussian2DFunction func =
        new EllipticalGaussian2DFunction(npeaks, blockWidth, blockWidth);
    params[0] = random(rng, background);
    for (int i = 0, j = 0; i < npeaks; i++, j += Gaussian2DFunction.PARAMETERS_PER_PEAK) {
      params[j + Gaussian2DFunction.SIGNAL] = random(rng, amplitude);
      params[j + Gaussian2DFunction.ANGLE] = random(rng, angle);
      params[j + Gaussian2DFunction.X_POSITION] = random(rng, xpos);
      params[j + Gaussian2DFunction.Y_POSITION] = random(rng, ypos);
      params[j + Gaussian2DFunction.X_SD] = random(rng, xwidth);
      params[j + Gaussian2DFunction.Y_SD] = random(rng, ywidth);
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
        params[j + Gaussian2DFunction.ANGLE] = random(rng, params[j + Gaussian2DFunction.ANGLE]);
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

  protected ArrayList<double[]> copyList(ArrayList<double[]> paramsList) {
    final ArrayList<double[]> params2List = new ArrayList<>(paramsList.size());
    for (int i = 0; i < paramsList.size(); i++) {
      params2List.add(paramsList.get(i).clone());
    }
    return params2List;
  }
}
