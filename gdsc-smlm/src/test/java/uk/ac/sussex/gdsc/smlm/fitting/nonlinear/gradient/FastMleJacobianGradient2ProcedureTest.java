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
import org.apache.commons.math3.util.Precision;
import org.ejml.data.DenseMatrix64F;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.Test;
import uk.ac.sussex.gdsc.smlm.function.FakeGradientFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.HoltzerAstigmatismZModel;
import uk.ac.sussex.gdsc.smlm.function.gaussian.erf.ErfGaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.erf.MultiFreeCircularErfGaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.erf.SingleAstigmatismErfGaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.erf.SingleFreeCircularErfGaussian2DFunction;
import uk.ac.sussex.gdsc.test.api.Predicates;
import uk.ac.sussex.gdsc.test.api.TestAssertions;
import uk.ac.sussex.gdsc.test.api.function.DoubleDoubleBiPredicate;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.junit5.SpeedTag;
import uk.ac.sussex.gdsc.test.rng.RngFactory;
import uk.ac.sussex.gdsc.test.utils.AssertionErrorCounter;
import uk.ac.sussex.gdsc.test.utils.RandomSeed;
import uk.ac.sussex.gdsc.test.utils.functions.IndexSupplier;

/**
 * Contains speed tests for the methods for calculating the Hessian and gradient vector for use in
 * the LVM algorithm.
 */
@SuppressWarnings({"javadoc"})
class FastMleJacobianGradient2ProcedureTest extends FastMleGradient2ProcedureTest {
  // Skip super-class tests ...
  @Override
  @Test
  void gradientProcedureFactoryCreatesOptimisedProcedures() {
    Assumptions.assumeTrue(false);
  }

  @Override
  @SeededTest
  void gradientProcedureComputesSameLogLikelihoodAsMleGradientCalculator(RandomSeed seed) {
    Assumptions.assumeTrue(false);
  }

  @Override
  @SpeedTag
  @SeededTest
  void gradientProcedureIsNotSlowerThanGradientCalculator(RandomSeed seed) {
    Assumptions.assumeTrue(false);
  }

  @Override
  @SeededTest
  void gradientProcedureComputesSameWithPrecomputed(RandomSeed seed) {
    Assumptions.assumeTrue(false);
  }

  @Override
  @SeededTest
  void gradientProcedureUnrolledComputesSameAsGradientProcedure(RandomSeed seed) {
    Assumptions.assumeTrue(false);
  }

  @Override
  @SpeedTag
  @SeededTest
  void gradientProcedureIsFasterUnrolledThanGradientProcedure(RandomSeed seed) {
    Assumptions.assumeTrue(false);
  }

  @SeededTest
  void gradientProcedureComputesSameAsBaseGradientProcedure(RandomSeed seed) {
    // Test the base functionality of computing the partial derivatives is the same
    final DoubleDoubleBiPredicate equality = Predicates.doublesAreClose(1e-5, 0);
    gradientProcedureComputesSameAsBaseGradientProcedure(seed, 4, equality);
    gradientProcedureComputesSameAsBaseGradientProcedure(seed, 5, equality);
    gradientProcedureComputesSameAsBaseGradientProcedure(seed, 6, equality);
    gradientProcedureComputesSameAsBaseGradientProcedure(seed, 11, equality);
    gradientProcedureComputesSameAsBaseGradientProcedure(seed, 21, equality);
  }

  private void gradientProcedureComputesSameAsBaseGradientProcedure(RandomSeed seed, int nparams,
      DoubleDoubleBiPredicate equality) {
    final int iter = 10;

    final ArrayList<double[]> paramsList = new ArrayList<>(iter);
    final ArrayList<double[]> yList = new ArrayList<>(iter);

    createFakeData(RngFactory.create(seed.get()), nparams, iter, paramsList, yList);
    final FakeGradientFunction func = new FakeGradientFunction(blockWidth, nparams);

    for (int i = 0; i < paramsList.size(); i++) {
      final FastMleGradient2Procedure p =
          FastMleGradient2ProcedureUtils.createUnrolled(yList.get(i), func);
      final FastMleJacobianGradient2Procedure p2 =
          new FastMleJacobianGradient2Procedure(yList.get(i), func);
      p.computeSecondDerivative(paramsList.get(i));
      p2.computeSecondDerivative(paramsList.get(i));
      // Virtually the same ...
      TestAssertions.assertArrayTest(p.d1, p2.d1, equality);
      TestAssertions.assertArrayTest(p.d2, p2.d2, equality);
    }
  }

  @Override
  @SeededTest
  void gradientCalculatorComputesGradient(RandomSeed seed) {
    gradientCalculatorComputesGradient(seed, 1,
        new SingleFreeCircularErfGaussian2DFunction(blockWidth, blockWidth));
    gradientCalculatorComputesGradient(seed, 2,
        new MultiFreeCircularErfGaussian2DFunction(2, blockWidth, blockWidth));

    // Use a reasonable z-depth function from the Smith, et al (2010) paper (page 377)
    final double sx = 1.08;
    final double sy = 1.01;
    final double gamma = 0.389;
    final double d = 0.531;
    final double Ax = -0.0708;
    final double Bx = -0.073;
    final double Ay = 0.164;
    final double By = 0.0417;
    final HoltzerAstigmatismZModel zModel =
        HoltzerAstigmatismZModel.create(sx, sy, gamma, d, Ax, Bx, Ay, By);

    gradientCalculatorComputesGradient(seed, 1,
        new SingleAstigmatismErfGaussian2DFunction(blockWidth, blockWidth, zModel));
  }

  private void gradientCalculatorComputesGradient(RandomSeed seed, int npeaks,
      ErfGaussian2DFunction func) {
    // Check the first and second derivatives
    final int nparams = func.getNumberOfGradients();
    final int[] indices = func.gradientIndices();

    final int iter = 100;

    final ArrayList<double[]> paramsList = new ArrayList<>(iter);
    final ArrayList<double[]> yList = new ArrayList<>(iter);

    createData(RngFactory.create(seed.get()), npeaks, iter, paramsList, yList, true);

    // for the gradients
    final double delta = 1e-4;
    final DoubleDoubleBiPredicate eq = Predicates.doublesAreClose(5e-2, 1e-16);
    final IndexSupplier msg1 = new IndexSupplier(2).setMessagePrefix("Gradient1 ");
    final IndexSupplier msg2 = new IndexSupplier(2).setMessagePrefix("Gradient2 ");
    final IndexSupplier msg3 = new IndexSupplier(3).setMessagePrefix("GradientJ ");

    // Must compute most of the time
    final int failureLimit = AssertionErrorCounter.computeFailureLimit(iter, 0.1);
    // failureLimit = 0;
    final AssertionErrorCounter failCounter = new AssertionErrorCounter(failureLimit, 2 * nparams);
    final AssertionErrorCounter failCounter2 =
        new AssertionErrorCounter(failureLimit, nparams * nparams);

    for (int i = 0; i < paramsList.size(); i++) {
      msg1.set(0, i);
      msg2.set(0, i);
      msg3.set(0, i);
      final double[] y = yList.get(i);
      final double[] a = paramsList.get(i);
      final double[] a2 = a.clone();
      final FastMleJacobianGradient2Procedure p = new FastMleJacobianGradient2Procedure(y, func);
      // double ll = p.computeLogLikelihood(a);
      p.computeJacobian(a);
      final double[] d1 = p.d1.clone();
      final double[] d2 = p.d2.clone();
      final DenseMatrix64F J = DenseMatrix64F.wrap(nparams, nparams, p.getJacobianLinear());
      for (int j = 0; j < nparams; j++) {
        final int j_ = j;
        msg1.set(1, j);
        msg2.set(1, j);
        msg3.set(1, j);
        final int k = indices[j];
        final double d = Precision.representableDelta(a[k], (a[k] == 0) ? delta : a[k] * delta);
        a2[k] = a[k] + d;
        final double llh = p.computeLogLikelihood(a2);
        p.computeFirstDerivative(a2);
        final double[] d1h = p.d1.clone();
        a2[k] = a[k] - d;
        final double lll = p.computeLogLikelihood(a2);
        p.computeFirstDerivative(a2);
        final double[] d1l = p.d1.clone();
        a2[k] = a[k];

        final double gradient1 = (llh - lll) / (2 * d);
        final double gradient2 = (d1h[j] - d1l[j]) / (2 * d);
        // logger.fine(FormatSupplier.getSupplier("[%d,%d] ll - %f (%s %f+/-%f) d1 %f ?= %f : d2 %f
        // ?= %f", i, k, ll, func.getName(k), a[k], d,
        // gradient1, d1[j], gradient2, d2[j]);
        failCounter.run(j, () -> TestAssertions.assertTest(gradient1, d1[j_], eq, msg1));
        failCounter.run(nparams + j, () -> TestAssertions.assertTest(gradient2, d2[j_], eq, msg2));

        // Test the Jacobian ...

        for (int jj = 0; jj < nparams; jj++) {
          if (j == jj) {
            // This is done above
            // Check it anyway to ensure the Jacobian is correct
            // continue;
          }

          final int jj_ = jj;
          msg3.set(2, jj);
          final int kk = indices[jj];
          final double dd =
              Precision.representableDelta(a[kk], (a[kk] == 0) ? delta : a[kk] * delta);
          a2[kk] = a[kk] + dd;
          p.computeFirstDerivative(a2);
          System.arraycopy(p.d1, 0, d1h, 0, d1h.length);
          a2[kk] = a[kk] - dd;
          p.computeFirstDerivative(a2);
          System.arraycopy(p.d1, 0, d1l, 0, d1l.length);
          a2[kk] = a[kk];

          // Use index j even though we adjusted index jj
          final double gradient3 = (d1h[j] - d1l[j]) / (2 * dd);
          // logger.fine(FormatSupplier.getSupplier("[%d,%d,%d] (%s %f %s %f+/-%f) J %f ?= %f %b",
          // i,
          // k, kk, func.getName(k),
          // a[k], func.getName(kk), a[kk], dd, gradient3, J.get(j, jj), ok);
          // if (!ok)
          // {
          // ExtraAssertions.fail("Not same gradientJ @ [%d,%d]", j, jj);
          // }
          failCounter2.run(nparams * j + jj,
              () -> TestAssertions.assertTest(gradient3, J.get(j_, jj_), eq, msg3));
        }
      }
    }
  }
}
