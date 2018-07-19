/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
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

import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.util.Precision;
import org.ejml.data.DenseMatrix64F;
import org.junit.Assume;
import org.junit.Test;

import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.smlm.function.FakeGradientFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.HoltzerAstigmatismZModel;
import uk.ac.sussex.gdsc.smlm.function.gaussian.erf.ErfGaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.erf.MultiFreeCircularErfGaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.erf.SingleAstigmatismErfGaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.erf.SingleFreeCircularErfGaussian2DFunction;
import uk.ac.sussex.gdsc.test.TestCounter;
import uk.ac.sussex.gdsc.test.TestSettings;
import uk.ac.sussex.gdsc.test.junit4.TestAssert;

/**
 * Contains speed tests for the methods for calculating the Hessian and gradient vector
 * for use in the LVM algorithm.
 */
@SuppressWarnings({ "javadoc" })
public class FastMLEJacobianGradient2ProcedureTest extends FastMLEGradient2ProcedureTest
{
	// Skip super-class tests ...
	@Override
	@Test
	public void gradientProcedureFactoryCreatesOptimisedProcedures()
	{
		Assume.assumeTrue(false);
	}

	@Override
	@Test
	public void gradientProcedureComputesSameLogLikelihoodAsMLEGradientCalculator()
	{
		Assume.assumeTrue(false);
	}

	@Override
	@Test
	public void gradientProcedureIsNotSlowerThanGradientCalculator()
	{
		Assume.assumeTrue(false);
	}

	@Override
	@Test
	public void gradientProcedureComputesSameWithPrecomputed()
	{
		Assume.assumeTrue(false);
	}

	@Override
	@Test
	public void gradientProcedureUnrolledComputesSameAsGradientProcedure()
	{
		Assume.assumeTrue(false);
	}

	@Override
	@Test
	public void gradientProcedureIsFasterUnrolledThanGradientProcedure()
	{
		Assume.assumeTrue(false);
	}

	@Test
	public void gradientProcedureComputesSameAsBaseGradientProcedure()
	{
		// Test the base functionality of computing the partial derivatives is the same
		gradientProcedureComputesSameAsBaseGradientProcedure(4);
		gradientProcedureComputesSameAsBaseGradientProcedure(5);
		gradientProcedureComputesSameAsBaseGradientProcedure(6);
		gradientProcedureComputesSameAsBaseGradientProcedure(11);
		gradientProcedureComputesSameAsBaseGradientProcedure(21);
	}

	private void gradientProcedureComputesSameAsBaseGradientProcedure(int nparams)
	{
		final int iter = 10;
		rdg = new RandomDataGenerator(TestSettings.getRandomGenerator());

		final ArrayList<double[]> paramsList = new ArrayList<>(iter);
		final ArrayList<double[]> yList = new ArrayList<>(iter);

		createFakeData(nparams, iter, paramsList, yList);
		final FakeGradientFunction func = new FakeGradientFunction(blockWidth, nparams);

		for (int i = 0; i < paramsList.size(); i++)
		{
			final FastMLEGradient2Procedure p = FastMLEGradient2ProcedureFactory.createUnrolled(yList.get(i), func);
			final FastMLEJacobianGradient2Procedure p2 = new FastMLEJacobianGradient2Procedure(yList.get(i), func);
			p.computeSecondDerivative(paramsList.get(i));
			p2.computeSecondDerivative(paramsList.get(i));
			// Virtually the same ...
			TestAssert.assertArrayEqualsRelative(p.d1, p2.d1, 1e-5);
			TestAssert.assertArrayEqualsRelative(p.d2, p2.d2, 1e-5);
		}
	}

	@Override
	@Test
	public void gradientCalculatorComputesGradient()
	{
		gradientCalculatorComputesGradient(1, new SingleFreeCircularErfGaussian2DFunction(blockWidth, blockWidth));
		gradientCalculatorComputesGradient(2, new MultiFreeCircularErfGaussian2DFunction(2, blockWidth, blockWidth));

		// Use a reasonable z-depth function from the Smith, et al (2010) paper (page 377)
		final double sx = 1.08;
		final double sy = 1.01;
		final double gamma = 0.389;
		final double d = 0.531;
		final double Ax = -0.0708;
		final double Bx = -0.073;
		final double Ay = 0.164;
		final double By = 0.0417;
		final HoltzerAstigmatismZModel zModel = HoltzerAstigmatismZModel.create(sx, sy, gamma, d, Ax, Bx, Ay, By);

		gradientCalculatorComputesGradient(1,
				new SingleAstigmatismErfGaussian2DFunction(blockWidth, blockWidth, zModel));
	}

	private void gradientCalculatorComputesGradient(int nPeaks, ErfGaussian2DFunction func)
	{
		// Check the first and second derivatives
		final int nparams = func.getNumberOfGradients();
		final int[] indices = func.gradientIndices();

		final int iter = 100;
		rdg = new RandomDataGenerator(TestSettings.getRandomGenerator());

		final ArrayList<double[]> paramsList = new ArrayList<>(iter);
		final ArrayList<double[]> yList = new ArrayList<>(iter);

		createData(nPeaks, iter, paramsList, yList, true);

		// for the gradients
		final double delta = 1e-4;
		final DoubleEquality eq = new DoubleEquality(5e-2, 1e-16);

		// Must compute most of the time
		final int failureLimit = TestCounter.computeFailureLimit(iter, 0.1);
		//failureLimit = 0;
		final TestCounter failCounter = new TestCounter(failureLimit, 2 * nparams);
		final TestCounter failCounter2 = new TestCounter(failureLimit, nparams * nparams);

		for (int i = 0; i < paramsList.size(); i++)
		{
			final int ii = i;
			final double[] y = yList.get(i);
			final double[] a = paramsList.get(i);
			final double[] a2 = a.clone();
			final FastMLEJacobianGradient2Procedure p = new FastMLEJacobianGradient2Procedure(y, func);
			//double ll = p.computeLogLikelihood(a);
			p.computeJacobian(a);
			final double[] d1 = p.d1.clone();
			final double[] d2 = p.d2.clone();
			final DenseMatrix64F J = DenseMatrix64F.wrap(nparams, nparams, p.getJacobianLinear());
			for (int j = 0; j < nparams; j++)
			{
				final int j_ = j;
				final int k = indices[j];
				final double d = Precision.representableDelta(a[k], (a[k] == 0) ? delta : a[k] * delta);
				a2[k] = a[k] + d;
				final double llh = p.computeLogLikelihood(a2);
				p.computeFirstDerivative(a2);
				double[] d1h = p.d1.clone();
				a2[k] = a[k] - d;
				final double lll = p.computeLogLikelihood(a2);
				p.computeFirstDerivative(a2);
				double[] d1l = p.d1.clone();
				a2[k] = a[k];

				final double gradient1 = (llh - lll) / (2 * d);
				final double gradient2 = (d1h[j] - d1l[j]) / (2 * d);
				//System.out.printf("[%d,%d] ll - %f  (%s %f+/-%f) d1 %f ?= %f : d2 %f ?= %f\n", i, k, ll, func.getName(k), a[k], d,
				//		gradient1, d1[j], gradient2, d2[j]);
				failCounter.run(j, () -> {
					return eq.almostEqualRelativeOrAbsolute(gradient1, d1[j_]);
				}, () -> {
					TestAssert.fail("Not same gradient1 @ %d,%d: %s != %s (error=%s)", ii, j_, gradient1, d1[j_],
							DoubleEquality.relativeError(gradient1, d1[j_]));
				});
				failCounter.run(nparams + j, () -> {
					return eq.almostEqualRelativeOrAbsolute(gradient2, d2[j_]);
				}, () -> {
					TestAssert.fail("Not same gradient2 @ %d,%d: %s != %s (error=%s)", ii, j_, gradient2, d2[j_],
							DoubleEquality.relativeError(gradient2, d2[j_]));
				});

				// Test the Jacobian ...

				for (int jj = 0; jj < nparams; jj++)
				{
					if (j == jj)
					{
						// This is done above
						// Check it anyway to ensure the Jacobian is correct
						//continue;
					}

					final int jj_ = jj;
					final int kk = indices[jj];
					final double dd = Precision.representableDelta(a[kk], (a[kk] == 0) ? delta : a[kk] * delta);
					a2[kk] = a[kk] + dd;
					p.computeFirstDerivative(a2);
					d1h = p.d1.clone();
					a2[kk] = a[kk] - dd;
					p.computeFirstDerivative(a2);
					d1l = p.d1.clone();
					a2[kk] = a[kk];

					// Use index j even though we adjusted index jj
					final double gradient3 = (d1h[j] - d1l[j]) / (2 * dd);
					final boolean ok = eq.almostEqualRelativeOrAbsolute(gradient3, J.get(j, jj));
					//System.out.printf("[%d,%d,%d] (%s %f  %s %f+/-%f) J %f ?= %f  %b\n", i, k, kk, func.getName(k),
					//		a[k], func.getName(kk), a[kk], dd, gradient3, J.get(j, jj), ok);
					//if (!ok)
					//{
					//	TestAssert.fail("Not same gradientJ @ [%d,%d]", j, jj);
					//}
					failCounter2.run(nparams * j_ + jj_, () -> {
						return ok;
					}, () -> {
						TestAssert.fail("Not same gradientJ @ %d [%d,%d]: %s != %s (error=%s)", ii, j_, jj_, gradient3,
								J.get(j_, jj_), DoubleEquality.relativeError(gradient3, J.get(j_, jj_)));
					});
				}
			}
		}
	}
}
