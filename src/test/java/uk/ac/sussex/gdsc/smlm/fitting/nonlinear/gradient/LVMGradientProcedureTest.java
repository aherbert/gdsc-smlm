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
import java.util.Arrays;

import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.util.Precision;
import org.ejml.data.DenseMatrix64F;
import org.junit.Assert;
import org.junit.Test;

import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.Maths;
import uk.ac.sussex.gdsc.core.utils.Random;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.LVMGradientProcedureFactory.Type;
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
import uk.ac.sussex.gdsc.test.TestCounter;
import uk.ac.sussex.gdsc.test.TestLog;
import uk.ac.sussex.gdsc.test.TestSettings;
import uk.ac.sussex.gdsc.test.junit4.TestAssert;
import uk.ac.sussex.gdsc.test.junit4.TestAssume;

/**
 * Contains speed tests for the methods for calculating the Hessian and gradient vector
 * for use in the LVM algorithm.
 */
@SuppressWarnings({ "javadoc" })
public class LVMGradientProcedureTest
{
	static FastLog fastLog = null;

	static FastLog getFastLog()
	{
		if (fastLog == null)
			// Default
			fastLog = FastLogFactory.getFastLog();
		return fastLog;
	}

	int MAX_ITER = 20000;
	int blockWidth = 10;
	double Noise = 0.3;
	double Background = 0.5;
	// High signal required to enabled the SupportsPrecomputed test to distinguish
	// LSQ and MLE/WLSQ. With a value of 100 it is possible to get the same gradient
	// with 3 peaks as with 2 peaks minus the third peak from the data.
	double Signal = 1000;
	double Angle = Math.PI;
	double Xpos = 5;
	double Ypos = 5;
	double Xwidth = 1.2;
	double Ywidth = 1.2;

	RandomDataGenerator rdg;

	@Test
	public void gradientProcedureFactoryCreatesOptimisedProcedures()
	{
		final DummyGradientFunction[] f = new DummyGradientFunction[7];
		for (int i = 1; i < f.length; i++)
			f[i] = new DummyGradientFunction(i);

		final LVMGradientProcedureFactory.Type MLE = LVMGradientProcedureFactory.Type.MLE;
		final LVMGradientProcedureFactory.Type WLSQ = LVMGradientProcedureFactory.Type.WLSQ;
		final LVMGradientProcedureFactory.Type LSQ = LVMGradientProcedureFactory.Type.LSQ;
		final LVMGradientProcedureFactory.Type FMLE = LVMGradientProcedureFactory.Type.FastLogMLE;

		final FastLog fl = getFastLog();

		//@formatter:off

		// Generic factory
		final double[] y0 = new double[1];
		final double[] y1 = new double[]{ 1 };
		Assert.assertEquals(LVMGradientProcedureFactory.create(y0, f[6], LSQ, fl).getClass(), LSQLVMGradientProcedure6.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y0, f[5], LSQ, fl).getClass(), LSQLVMGradientProcedure5.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y0, f[4], LSQ, fl).getClass(), LSQLVMGradientProcedure4.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y0, f[1], LSQ, fl).getClass(), LSQLVMGradientProcedure.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y0, f[6], MLE, fl).getClass(), MLELVMGradientProcedure6.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y0, f[5], MLE, fl).getClass(), MLELVMGradientProcedure5.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y0, f[4], MLE, fl).getClass(), MLELVMGradientProcedure4.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y0, f[1], MLE, fl).getClass(), MLELVMGradientProcedure.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y1, f[6], MLE, fl).getClass(), MLELVMGradientProcedureX6.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y1, f[5], MLE, fl).getClass(), MLELVMGradientProcedureX5.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y1, f[4], MLE, fl).getClass(), MLELVMGradientProcedureX4.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y1, f[1], MLE, fl).getClass(), MLELVMGradientProcedureX.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y0, f[6], WLSQ, fl).getClass(), WLSQLVMGradientProcedure6.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y0, f[5], WLSQ, fl).getClass(), WLSQLVMGradientProcedure5.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y0, f[4], WLSQ, fl).getClass(), WLSQLVMGradientProcedure4.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y0, f[1], WLSQ, fl).getClass(), WLSQLVMGradientProcedure.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y0, f[6], FMLE, fl).getClass(), FastLogMLELVMGradientProcedure6.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y0, f[5], FMLE, fl).getClass(), FastLogMLELVMGradientProcedure5.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y0, f[4], FMLE, fl).getClass(), FastLogMLELVMGradientProcedure4.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y0, f[1], FMLE, fl).getClass(), FastLogMLELVMGradientProcedure.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y1, f[6], FMLE, fl).getClass(), FastLogMLELVMGradientProcedureX6.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y1, f[5], FMLE, fl).getClass(), FastLogMLELVMGradientProcedureX5.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y1, f[4], FMLE, fl).getClass(), FastLogMLELVMGradientProcedureX4.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y1, f[1], FMLE, fl).getClass(), FastLogMLELVMGradientProcedureX.class);

		// Dedicated factories
		Assert.assertEquals(LSQLVMGradientProcedureFactory.create(y0, f[6]).getClass(), LSQLVMGradientProcedure6.class);
		Assert.assertEquals(LSQLVMGradientProcedureFactory.create(y0, f[5]).getClass(), LSQLVMGradientProcedure5.class);
		Assert.assertEquals(LSQLVMGradientProcedureFactory.create(y0, f[4]).getClass(), LSQLVMGradientProcedure4.class);
		Assert.assertEquals(LSQLVMGradientProcedureFactory.create(y0, f[1]).getClass(), LSQLVMGradientProcedure.class);
		Assert.assertEquals(MLELVMGradientProcedureFactory.create(y0, f[6]).getClass(), MLELVMGradientProcedure6.class);
		Assert.assertEquals(MLELVMGradientProcedureFactory.create(y0, f[5]).getClass(), MLELVMGradientProcedure5.class);
		Assert.assertEquals(MLELVMGradientProcedureFactory.create(y0, f[4]).getClass(), MLELVMGradientProcedure4.class);
		Assert.assertEquals(MLELVMGradientProcedureFactory.create(y0, f[1]).getClass(), MLELVMGradientProcedure.class);
		Assert.assertEquals(MLELVMGradientProcedureFactory.create(y1, f[6]).getClass(), MLELVMGradientProcedureX6.class);
		Assert.assertEquals(MLELVMGradientProcedureFactory.create(y1, f[5]).getClass(), MLELVMGradientProcedureX5.class);
		Assert.assertEquals(MLELVMGradientProcedureFactory.create(y1, f[4]).getClass(), MLELVMGradientProcedureX4.class);
		Assert.assertEquals(MLELVMGradientProcedureFactory.create(y1, f[1]).getClass(), MLELVMGradientProcedureX.class);
		Assert.assertEquals(WLSQLVMGradientProcedureFactory.create(y0, null, f[6]).getClass(), WLSQLVMGradientProcedure6.class);
		Assert.assertEquals(WLSQLVMGradientProcedureFactory.create(y0, null, f[5]).getClass(), WLSQLVMGradientProcedure5.class);
		Assert.assertEquals(WLSQLVMGradientProcedureFactory.create(y0, null, f[4]).getClass(), WLSQLVMGradientProcedure4.class);
		Assert.assertEquals(WLSQLVMGradientProcedureFactory.create(y0, null, f[1]).getClass(), WLSQLVMGradientProcedure.class);
		Assert.assertEquals(MLELVMGradientProcedureFactory.create(y0, f[6], fl).getClass(), FastLogMLELVMGradientProcedure6.class);
		Assert.assertEquals(MLELVMGradientProcedureFactory.create(y0, f[5], fl).getClass(), FastLogMLELVMGradientProcedure5.class);
		Assert.assertEquals(MLELVMGradientProcedureFactory.create(y0, f[4], fl).getClass(), FastLogMLELVMGradientProcedure4.class);
		Assert.assertEquals(MLELVMGradientProcedureFactory.create(y0, f[1], fl).getClass(), FastLogMLELVMGradientProcedure.class);
		Assert.assertEquals(MLELVMGradientProcedureFactory.create(y1, f[6], fl).getClass(), FastLogMLELVMGradientProcedureX6.class);
		Assert.assertEquals(MLELVMGradientProcedureFactory.create(y1, f[5], fl).getClass(), FastLogMLELVMGradientProcedureX5.class);
		Assert.assertEquals(MLELVMGradientProcedureFactory.create(y1, f[4], fl).getClass(), FastLogMLELVMGradientProcedureX4.class);
		Assert.assertEquals(MLELVMGradientProcedureFactory.create(y1, f[1], fl).getClass(), FastLogMLELVMGradientProcedureX.class);

		//@formatter:on
	}

	@Test
	public void gradientProcedureLSQComputesSameAsGradientCalculator()
	{
		gradientProcedureComputesSameAsGradientCalculator(Type.LSQ);
	}

	@Test
	public void gradientProcedureMLEComputesSameAsGradientCalculator()
	{
		gradientProcedureComputesSameAsGradientCalculator(Type.MLE);
	}

	@Test
	public void gradientProcedureFastLogMLEComputesSameAsGradientCalculator()
	{
		gradientProcedureComputesSameAsGradientCalculator(Type.FastLogMLE, 1e-5);
	}

	private void gradientProcedureComputesSameAsGradientCalculator(Type type)
	{
		gradientProcedureComputesSameAsGradientCalculator(type, 0);
	}

	private void gradientProcedureComputesSameAsGradientCalculator(Type type, double error)
	{
		gradientProcedureComputesSameAsGradientCalculator(4, type, error);
		gradientProcedureComputesSameAsGradientCalculator(5, type, error);
		gradientProcedureComputesSameAsGradientCalculator(6, type, error);
		gradientProcedureComputesSameAsGradientCalculator(11, type, error);
		gradientProcedureComputesSameAsGradientCalculator(21, type, error);
	}

	@Test
	public void gradientProcedureLSQIsNotSlowerThanGradientCalculator()
	{
		gradientProcedureIsNotSlowerThanGradientCalculator(Type.LSQ);
	}

	@Test
	public void gradientProcedureMLEIsNotSlowerThanGradientCalculator()
	{
		gradientProcedureIsNotSlowerThanGradientCalculator(Type.MLE);
	}

	@Test
	public void gradientProcedureFastLogMLEIsNotSlowerThanGradientCalculator()
	{
		gradientProcedureIsNotSlowerThanGradientCalculator(Type.FastLogMLE);
	}

	private void gradientProcedureIsNotSlowerThanGradientCalculator(Type type)
	{
		gradientProcedureIsNotSlowerThanGradientCalculator(4, type);
		gradientProcedureIsNotSlowerThanGradientCalculator(5, type);
		gradientProcedureIsNotSlowerThanGradientCalculator(6, type);
		// 2 peaks
		gradientProcedureIsNotSlowerThanGradientCalculator(11, type);
		// 4 peaks
		gradientProcedureIsNotSlowerThanGradientCalculator(21, type);
	}

	private void gradientProcedureComputesSameAsGradientCalculator(int nparams, Type type, double error)
	{
		final int iter = 10;
		rdg = new RandomDataGenerator(TestSettings.getRandomGenerator());

		final double[][] alpha = new double[nparams][nparams];
		final double[] beta = new double[nparams];

		final ArrayList<double[]> paramsList = new ArrayList<>(iter);
		final ArrayList<double[]> yList = new ArrayList<>(iter);

		final int[] x = createFakeData(nparams, iter, paramsList, yList);
		final int n = x.length;
		final FakeGradientFunction func = new FakeGradientFunction(blockWidth, nparams);

		final boolean mle = type != Type.LSQ;
		final FastLog fastLog = (type == Type.FastLogMLE) ? getFastLog() : null;
		final GradientCalculator calc = GradientCalculatorFactory.newCalculator(nparams, mle);

		final String name = String.format("[%d] %b", nparams, mle);

		for (int i = 0; i < paramsList.size(); i++)
		{
			// Reference implementation
			final double s = calc.findLinearised(n, yList.get(i), paramsList.get(i), alpha, beta, func);
			// Procedure
			final LVMGradientProcedure p = LVMGradientProcedureFactory.create(yList.get(i), func, type, fastLog);
			p.gradient(paramsList.get(i));
			final double s2 = p.value;
			// Value may be different depending on log implementation
			Assert.assertEquals(name + " Result: Not same @ " + i, s, s2, Math.abs(s2) * error);
			// Exactly the same ...
			Assert.assertArrayEquals(name + " Observations: Not same beta @ " + i, p.beta, beta, 0);

			final double[] al = p.getAlphaLinear();
			Assert.assertArrayEquals(name + " Observations: Not same alpha @ " + i, al, new DenseMatrix64F(alpha).data,
					0);

			final double[][] am = p.getAlphaMatrix();
			for (int j = 0; j < nparams; j++)
				Assert.assertArrayEquals(name + " Observations: Not same alpha @ " + i, am[j], alpha[j], 0);
		}
	}

	private abstract class Timer
	{
		private int loops;
		int min;

		Timer()
		{
		}

		Timer(int min)
		{
			this.min = min;
		}

		long getTime()
		{
			// Run till stable timing
			long t1 = time();
			for (int i = 0; i < 10; i++)
			{
				final long t2 = t1;
				t1 = time();
				if (loops >= min && DoubleEquality.relativeError(t1, t2) < 0.02) // 2% difference
					break;
			}
			return t1;
		}

		long time()
		{
			loops++;
			long t = System.nanoTime();
			run();
			t = System.nanoTime() - t;
			//System.out.printf("[%d] Time = %d\n", loops, t);
			return t;
		}

		abstract void run();
	}

	private void gradientProcedureIsNotSlowerThanGradientCalculator(final int nparams, final Type type)
	{
		TestAssume.assumeSpeedTest();

		final int iter = 1000;
		rdg = new RandomDataGenerator(TestSettings.getRandomGenerator());
		final double[][] alpha = new double[nparams][nparams];
		final double[] beta = new double[nparams];

		final ArrayList<double[]> paramsList = new ArrayList<>(iter);
		final ArrayList<double[]> yList = new ArrayList<>(iter);

		final int[] x = createFakeData(nparams, iter, paramsList, yList);
		final int n = x.length;
		final FakeGradientFunction func = new FakeGradientFunction(blockWidth, nparams);

		final boolean mle = type != Type.LSQ;
		final FastLog fastLog = (type == Type.FastLogMLE) ? getFastLog() : null;
		final GradientCalculator calc = GradientCalculatorFactory.newCalculator(nparams, mle);

		for (int i = 0; i < paramsList.size(); i++)
			calc.findLinearised(n, yList.get(i), paramsList.get(i), alpha, beta, func);

		for (int i = 0; i < paramsList.size(); i++)
		{
			final LVMGradientProcedure p = LVMGradientProcedureFactory.create(yList.get(i), func, type, fastLog);
			p.gradient(paramsList.get(i));
		}

		// Realistic loops for an optimisation
		final int loops = 15;

		// Run till stable timing
		final Timer t1 = new Timer()
		{
			@Override
			void run()
			{
				for (int i = 0, k = 0; i < iter; i++)
				{
					final GradientCalculator calc = GradientCalculatorFactory.newCalculator(nparams, mle);
					for (int j = loops; j-- > 0;)
						calc.findLinearised(n, yList.get(i), paramsList.get(k++ % iter), alpha, beta, func);
				}
			}
		};
		final long time1 = t1.getTime();

		final Timer t2 = new Timer(t1.loops)
		{
			@Override
			void run()
			{
				for (int i = 0, k = 0; i < iter; i++)
				{
					final LVMGradientProcedure p = LVMGradientProcedureFactory.create(yList.get(i), func, type, fastLog);
					for (int j = loops; j-- > 0;)
						p.gradient(paramsList.get(k++ % iter));
				}
			}
		};
		final long time2 = t2.getTime();

		TestLog.logSpeedTestResult(time2 < time1,
				"GradientCalculator = %d : LVMGradientProcedure %d %s = %d : %fx\n", time1, nparams, type, time2,
				(1.0 * time1) / time2);
	}

	@Test
	public void gradientProcedureLSQUnrolledComputesSameAsGradientProcedure()
	{
		gradientProcedureUnrolledComputesSameAsGradientProcedure(Type.LSQ, false);
	}

	@Test
	public void gradientProcedureMLEUnrolledComputesSameAsGradientProcedure()
	{
		gradientProcedureUnrolledComputesSameAsGradientProcedure(Type.MLE, false);
	}

	@Test
	public void gradientProcedureFastLogMLEUnrolledComputesSameAsGradientProcedure()
	{
		gradientProcedureUnrolledComputesSameAsGradientProcedure(Type.FastLogMLE, false);
	}

	@Test
	public void gradientProcedureWLSQUnrolledComputesSameAsGradientProcedure()
	{
		gradientProcedureUnrolledComputesSameAsGradientProcedure(Type.WLSQ, false);
	}

	@Test
	public void gradientProcedureLSQUnrolledComputesSameAsGradientProcedureWithPrecomputed()
	{
		gradientProcedureUnrolledComputesSameAsGradientProcedure(Type.LSQ, true);
	}

	@Test
	public void gradientProcedureMLEUnrolledComputesSameAsGradientProcedureWithPrecomputed()
	{
		gradientProcedureUnrolledComputesSameAsGradientProcedure(Type.MLE, true);
	}

	@Test
	public void gradientProcedureFastLogMLEUnrolledComputesSameAsGradientProcedureWithPrecomputed()
	{
		gradientProcedureUnrolledComputesSameAsGradientProcedure(Type.FastLogMLE, true);
	}

	@Test
	public void gradientProcedureWLSQUnrolledComputesSameAsGradientProcedureWithPrecomputed()
	{
		gradientProcedureUnrolledComputesSameAsGradientProcedure(Type.WLSQ, true);
	}

	private void gradientProcedureUnrolledComputesSameAsGradientProcedure(Type type, boolean precomputed)
	{
		gradientProcedureUnrolledComputesSameAsGradientProcedure(4, type, precomputed);
		gradientProcedureUnrolledComputesSameAsGradientProcedure(5, type, precomputed);
		gradientProcedureUnrolledComputesSameAsGradientProcedure(6, type, precomputed);
	}

	private void gradientProcedureUnrolledComputesSameAsGradientProcedure(int nparams, Type type, boolean precomputed)
	{
		final int iter = 10;
		rdg = new RandomDataGenerator(TestSettings.getRandomGenerator());

		final ArrayList<double[]> paramsList = new ArrayList<>(iter);
		final ArrayList<double[]> yList = new ArrayList<>(iter);

		createFakeData(nparams, iter, paramsList, yList);
		Gradient1Function func = new FakeGradientFunction(blockWidth, nparams);

		if (precomputed)
		{
			final double[] b = SimpleArrayUtils.newArray(func.size(), 0.1, 1.3);
			func = OffsetGradient1Function.wrapGradient1Function(func, b);
		}

		final FastLog fastLog = type == Type.FastLogMLE ? getFastLog() : null;

		final String name = String.format("[%d] %b", nparams, type);
		for (int i = 0; i < paramsList.size(); i++)
		{
			final LVMGradientProcedure p1 = createProcedure(type, yList.get(i), func, fastLog);
			p1.gradient(paramsList.get(i));

			final LVMGradientProcedure p2 = LVMGradientProcedureFactory.create(yList.get(i), func, type, fastLog);
			p2.gradient(paramsList.get(i));

			// Exactly the same ...
			Assert.assertEquals(name + " Result: Not same @ " + i, p1.value, p2.value, 0);
			Assert.assertArrayEquals(name + " Observations: Not same beta @ " + i, p1.beta, p2.beta, 0);

			Assert.assertArrayEquals(name + " Observations: Not same alpha @ " + i, p1.getAlphaLinear(),
					p2.getAlphaLinear(), 0);

			final double[][] am1 = p1.getAlphaMatrix();
			final double[][] am2 = p2.getAlphaMatrix();
			for (int j = 0; j < nparams; j++)
				Assert.assertArrayEquals(name + " Observations: Not same alpha @ " + i, am1[j], am2[j], 0);
		}
	}

	private static LVMGradientProcedure createProcedure(Type type, double[] y, Gradient1Function func, FastLog fastLog)
	{
		switch (type)
		{
			case FastLogMLE:
				return new FastLogMLELVMGradientProcedure(y, func, fastLog);
			case MLE:
				return new MLELVMGradientProcedure(y, func);
			case WLSQ:
				return new WLSQLVMGradientProcedure(y, null, func);
			default:
				return new LSQLVMGradientProcedure(y, func);
		}
	}

	@Test
	public void gradientProcedureLSQIsFasterUnrolledThanGradientProcedure()
	{
		gradientProcedureIsFasterUnrolledThanGradientProcedure(Type.LSQ, false);
	}

	@Test
	public void gradientProcedureMLEIsFasterUnrolledThanGradientProcedure()
	{
		gradientProcedureIsFasterUnrolledThanGradientProcedure(Type.MLE, false);
	}

	@Test
	public void gradientProcedureFastLogMLEIsFasterUnrolledThanGradientProcedure()
	{
		gradientProcedureIsFasterUnrolledThanGradientProcedure(Type.FastLogMLE, false);
	}

	@Test
	public void gradientProcedureWLSQIsFasterUnrolledThanGradientProcedure()
	{
		gradientProcedureIsFasterUnrolledThanGradientProcedure(Type.WLSQ, false);
	}

	@Test
	public void gradientProcedureLSQIsFasterUnrolledThanGradientProcedureWithPrecomputed()
	{
		gradientProcedureIsFasterUnrolledThanGradientProcedure(Type.LSQ, true);
	}

	@Test
	public void gradientProcedureMLEIsFasterUnrolledThanGradientProcedureWithPrecomputed()
	{
		gradientProcedureIsFasterUnrolledThanGradientProcedure(Type.MLE, true);
	}

	@Test
	public void gradientProcedureFastLogMLEIsFasterUnrolledThanGradientProcedureWithPrecomputed()
	{
		gradientProcedureIsFasterUnrolledThanGradientProcedure(Type.FastLogMLE, true);
	}

	@Test
	public void gradientProcedureWLSQIsFasterUnrolledThanGradientProcedureWithPrecomputed()
	{
		gradientProcedureIsFasterUnrolledThanGradientProcedure(Type.WLSQ, true);
	}

	private void gradientProcedureIsFasterUnrolledThanGradientProcedure(Type type, boolean precomputed)
	{
		gradientProcedureIsFasterUnrolledThanGradientProcedure(4, type, precomputed);
		gradientProcedureIsFasterUnrolledThanGradientProcedure(5, type, precomputed);
		gradientProcedureIsFasterUnrolledThanGradientProcedure(6, type, precomputed);
	}

	private void gradientProcedureIsFasterUnrolledThanGradientProcedure(final int nparams, final Type type,
			final boolean precomputed)
	{
		TestAssume.assumeSpeedTest();

		final int iter = 100;
		rdg = new RandomDataGenerator(TestSettings.getRandomGenerator());

		final ArrayList<double[]> paramsList = new ArrayList<>(iter);
		final ArrayList<double[]> yList = new ArrayList<>(iter);

		createData(1, iter, paramsList, yList);

		// Remove the timing of the function call by creating a dummy function
		final FakeGradientFunction fgf = new FakeGradientFunction(blockWidth, nparams);
		final Gradient1Function func;
		if (precomputed)
		{
			final double[] b = SimpleArrayUtils.newArray(fgf.size(), 0.1, 1.3);
			func = OffsetGradient1Function.wrapGradient1Function(fgf, b);
		}
		else
			func = fgf;

		final FastLog fastLog = type == Type.FastLogMLE ? getFastLog() : null;

		for (int i = 0; i < paramsList.size(); i++)
		{
			final LVMGradientProcedure p1 = createProcedure(type, yList.get(i), func, fastLog);
			p1.gradient(paramsList.get(i));
			p1.gradient(paramsList.get(i));

			final LVMGradientProcedure p2 = LVMGradientProcedureFactory.create(yList.get(i), func, type, fastLog);
			p2.gradient(paramsList.get(i));
			p2.gradient(paramsList.get(i));

			// Check they are the same
			Assert.assertArrayEquals("A " + i, p1.getAlphaLinear(), p2.getAlphaLinear(), 0);
			Assert.assertArrayEquals("B " + i, p1.beta, p2.beta, 0);
		}

		// Realistic loops for an optimisation
		final int loops = 15;

		// Run till stable timing
		final Timer t1 = new Timer()
		{
			@Override
			void run()
			{
				for (int i = 0, k = 0; i < paramsList.size(); i++)
				{
					final LVMGradientProcedure p1 = createProcedure(type, yList.get(i), func, fastLog);
					for (int j = loops; j-- > 0;)
						p1.gradient(paramsList.get(k++ % iter));
				}
			}
		};
		final long time1 = t1.getTime();

		final Timer t2 = new Timer(t1.loops)
		{
			@Override
			void run()
			{
				for (int i = 0, k = 0; i < paramsList.size(); i++)
				{
					final LVMGradientProcedure p2 = LVMGradientProcedureFactory.create(yList.get(i), func, type, fastLog);
					for (int j = loops; j-- > 0;)
						p2.gradient(paramsList.get(k++ % iter));
				}
			}
		};
		final long time2 = t2.getTime();

		TestLog.logSpeedTestResult(time2 < time1, "%s, Precomputed=%b : Standard = %d : Unrolled %d = %d : %fx\n",
				type, precomputed, time1, nparams, time2, (1.0 * time1) / time2);
	}

	@Test
	public void gradientProcedureLSQComputesGradient()
	{
		gradientProcedureComputesGradient(new SingleFreeCircularErfGaussian2DFunction(blockWidth, blockWidth), Type.LSQ,
				false);
	}

	@Test
	public void gradientProcedureMLEComputesGradient()
	{
		gradientProcedureComputesGradient(new SingleFreeCircularErfGaussian2DFunction(blockWidth, blockWidth), Type.MLE,
				false);
	}

	@Test(expected = AssertionError.class)
	public void gradientProcedureFastLogMLECannotComputeGradient()
	{
		gradientProcedureComputesGradient(new SingleFreeCircularErfGaussian2DFunction(blockWidth, blockWidth),
				Type.FastLogMLE, false);
	}

	// This test now passes as the tolerance for computing the gradient has been lowered
	// so that the test passes under a stress test using many different random seeds.
	//@Test(expected = AssertionError.class)
	public void gradientProcedureFastLogMLECannotComputeGradientWithHighPrecision()
	{
		// Try different precision
		for (int n = FastLog.N; n < 23; n++)
		{
			try
			{
				//System.out.printf("Precision n=%d\n", n);
				fastLog = new TurboLog2(n);
				gradientProcedureComputesGradient(new SingleFreeCircularErfGaussian2DFunction(blockWidth, blockWidth),
						Type.FastLogMLE, false);
			}
			catch (final AssertionError e)
			{
				continue;
			}
			finally
			{
				// Reset
				fastLog = null;
			}
			return;
		}
		Assert.fail();
	}

	@Test
	public void gradientProcedureWLSQComputesGradient()
	{
		gradientProcedureComputesGradient(new SingleFreeCircularErfGaussian2DFunction(blockWidth, blockWidth),
				Type.WLSQ, false);
	}

	@Test
	public void gradientProcedureLSQComputesGradientWithPrecomputed()
	{
		gradientProcedureComputesGradient(new SingleFreeCircularErfGaussian2DFunction(blockWidth, blockWidth), Type.LSQ,
				true);
	}

	@Test
	public void gradientProcedureMLEComputesGradientWithPrecomputed()
	{
		gradientProcedureComputesGradient(new SingleFreeCircularErfGaussian2DFunction(blockWidth, blockWidth), Type.MLE,
				true);
	}

	@Test(expected = AssertionError.class)
	public void gradientProcedureFastLogMLECannotComputeGradientWithPrecomputed()
	{
		gradientProcedureComputesGradient(new SingleFreeCircularErfGaussian2DFunction(blockWidth, blockWidth),
				Type.FastLogMLE, true);
	}

	@Test
	public void gradientProcedureWLSQComputesGradientWithPrecomputed()
	{
		gradientProcedureComputesGradient(new SingleFreeCircularErfGaussian2DFunction(blockWidth, blockWidth),
				Type.WLSQ, true);
	}

	@SuppressWarnings("null")
	private void gradientProcedureComputesGradient(ErfGaussian2DFunction func, Type type, boolean precomputed)
	{
		final int nparams = func.getNumberOfGradients();
		final int[] indices = func.gradientIndices();

		final int iter = 100;
		rdg = new RandomDataGenerator(TestSettings.getRandomGenerator());

		final ArrayList<double[]> paramsList = new ArrayList<>(iter);
		final ArrayList<double[]> yList = new ArrayList<>(iter);

		createData(1, iter, paramsList, yList, true);

		// for the gradients
		final double delta = 1e-4;
		final DoubleEquality eq = new DoubleEquality(5e-2, 1e-16);

		final double[] b = (precomputed) ? new double[func.size()] : null;

		final FastLog fastLog = type == Type.FastLogMLE ? getFastLog() : null;

		// Must compute most of the time
		final int failureLimit = TestCounter.computeFailureLimit(iter, 0.1);
		final TestCounter failCounter = new TestCounter(failureLimit, nparams);

		for (int i = 0; i < paramsList.size(); i++)
		{
			final int ii = i;
			final double[] y = yList.get(i);
			final double[] a = paramsList.get(i);
			final double[] a2 = a.clone();

			LVMGradientProcedure p;
			if (precomputed)
			{
				// Mock fitting part of the function already
				for (int j = 0; j < b.length; j++)
					b[j] = y[j] * 0.5;
				p = LVMGradientProcedureFactory.create(y, OffsetGradient1Function.wrapGradient1Function(func, b), type,
						fastLog);
			}
			else
				p = LVMGradientProcedureFactory.create(y, func, type, fastLog);
			p.gradient(a);
			//double s = p.value;
			final double[] beta = p.beta.clone();
			for (int j = 0; j < nparams; j++)
			{
				final int jj = j;
				final int k = indices[j];
				//double d = Precision.representableDelta(a[k], (a[k] == 0) ? 1e-3 : a[k] * delta);
				final double d = Precision.representableDelta(a[k], delta);
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
				//System.out.printf("[%d,%d] %f  (%s %f+/-%f)  %f  ?=  %f\n", i, k, s, Gaussian2DFunction.getName(k),
				//		a[k], d, beta[j], gradient);

				failCounter.run(j, () -> {
					return eq.almostEqualRelativeOrAbsolute(beta[jj], gradient);
				}, () -> {
					TestAssert.fail("Not same gradient @ %d,%d: %s != %s (error=%s)", ii, jj, beta[jj], gradient,
							DoubleEquality.relativeError(beta[jj], gradient));
				});
			}
		}
	}

	@Test
	public void gradientProcedureLSQSupportsPrecomputed()
	{
		gradientProcedureSupportsPrecomputed(Type.LSQ);
	}

	@Test
	public void gradientProcedureMLESupportsPrecomputed()
	{
		gradientProcedureSupportsPrecomputed(Type.MLE);
	}

	@Test
	public void gradientProcedureFastLogMLESupportsPrecomputed()
	{
		gradientProcedureSupportsPrecomputed(Type.FastLogMLE, false);
	}

	@Test(expected = AssertionError.class)
	public void gradientProcedureFastLogMLECannotSupportPrecomputedWithGradients()
	{
		gradientProcedureSupportsPrecomputed(Type.FastLogMLE);
	}

	@Test
	public void gradientProcedureWLSQSupportsPrecomputed()
	{
		gradientProcedureSupportsPrecomputed(Type.WLSQ);
	}

	private void gradientProcedureSupportsPrecomputed(final Type type)
	{
		gradientProcedureSupportsPrecomputed(type, true);
	}

	private void gradientProcedureSupportsPrecomputed(final Type type, boolean checkGradients)
	{
		final int iter = 10;
		rdg = new RandomDataGenerator(TestSettings.getRandomGenerator());

		final ArrayList<double[]> paramsList = new ArrayList<>(iter);
		final ArrayList<double[]> yList = new ArrayList<>(iter);

		// 3 peaks
		createData(3, iter, paramsList, yList, true);

		for (int i = 0; i < paramsList.size(); i++)
		{
			final double[] y = yList.get(i);
			// Add Gaussian read noise so we have negatives
			final double min = Maths.min(y);
			for (int j = 0; j < y.length; j++)
				y[j] = y[i] - min + rdg.nextGaussian(0, Noise);
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

		final FastLog fastLog = type == Type.FastLogMLE ? getFastLog() : null;

		final int nparams = f12.getNumberOfGradients();
		final int[] indices = f12.gradientIndices();
		final double[] b = new double[f12.size()];

		final DoubleEquality eq = new DoubleEquality(1e-8, 1e-16); // for checking strict equivalence

		// for the gradients
		final double delta = 1e-4;
		final DoubleEquality eq2 = new DoubleEquality(5e-2, 1e-16);

		final double[] a1peaks = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
		final double[] y_b = new double[b.length];

		// Count the number of failures for each gradient
		final int failureLimit = TestCounter.computeFailureLimit(iter, 0.1);
		final TestCounter failCounter = new TestCounter(failureLimit, nparams * 2);

		for (int i = 0; i < paramsList.size(); i++)
		{
			final int ii = i;
			final double[] y = yList.get(i);
			final double[] a3peaks = paramsList.get(i);
			//System.out.printf("[%d] a=%s\n", i, Arrays.toString(a3peaks));

			final double[] a2peaks = Arrays.copyOf(a3peaks, 1 + 2 * Gaussian2DFunction.PARAMETERS_PER_PEAK);
			final double[] a2peaks2 = a2peaks.clone();
			for (int j = 1; j < a1peaks.length; j++)
				a1peaks[j] = a3peaks[j + 2 * Gaussian2DFunction.PARAMETERS_PER_PEAK];

			// Evaluate peak 3 to get the background and subtract it from the data to get the new data
			f3.initialise0(a1peaks);
			f3.forEach(new ValueProcedure()
			{
				int k = 0;

				@Override
				public void execute(double value)
				{
					b[k] = value;
					// Remove negatives for MLE
					if (type.isMLE())
					{
						y[k] = Math.max(0, y[k]);
						y_b[k] = Math.max(0, y[k] - value);
					}
					else
						y_b[k] = y[k] - value;
					k++;
				}
			});

			final LVMGradientProcedure p123 = LVMGradientProcedureFactory.create(y, f123, type, fastLog);

			/////////////////////////////////////
			// These should be the same
			/////////////////////////////////////
			final LVMGradientProcedure p12b3 = LVMGradientProcedureFactory.create(y,
					OffsetGradient1Function.wrapGradient1Function(f12, b), type, fastLog);

			// Check they are the same
			p123.gradient(a3peaks);
			final double[][] m123 = p123.getAlphaMatrix();

			p12b3.gradient(a2peaks);
			double s = p12b3.value;
			final double[] beta = p12b3.beta.clone();
			double[][] alpha = p12b3.getAlphaMatrix();

			//System.out.printf("MLE=%b [%d] p12b3  %f  %f\n", type.isMLE(), i, p123.value, s);

			if (!eq.almostEqualRelativeOrAbsolute(p123.value, s))
				TestAssert.fail("p12b3 Not same value @ %d (error=%s) : %s == %s", i,
						DoubleEquality.relativeError(p123.value, s), p123.value, s);
			if (!eq.almostEqualRelativeOrAbsolute(beta, p123.beta))
				TestAssert.fail("p12b3 Not same gradient @ %d (error=%s) : %s vs %s", i,
						DoubleEquality.relativeError(beta, p123.beta), Arrays.toString(beta),
						Arrays.toString(p123.beta));
			for (int j = 0; j < alpha.length; j++)
				//System.out.printf("%s !=\n%s\n", Arrays.toString(alpha[j]), Arrays.toString(m123[j]));
				if (!eq.almostEqualRelativeOrAbsolute(alpha[j], m123[j]))
					TestAssert.fail("p12b3 Not same alpha @ %d,%d (error=%s) : %s vs %s", i, j,
							DoubleEquality.relativeError(alpha[j], m123[j]), Arrays.toString(alpha[j]),
							Arrays.toString(m123[j]));

			// Check actual gradients are correct
			if (checkGradients)
				for (int j = 0; j < nparams; j++)
				{
					final int jj = j;
					final int k = indices[j];
					//double d = Precision.representableDelta(a2peaks[k], (a2peaks[k] == 0) ? 1e-3 : a2peaks[k] * delta);
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
					//System.out.printf("[%d,%d] %f  (%s %f+/-%f)  %f  ?=  %f  (%f)\n", i, k, s,
					//		Gaussian2DFunction.getName(k), a2peaks[k], d, beta[j], gradient,
					//		DoubleEquality.relativeError(gradient, beta[j]));
					failCounter.run(j, () -> {
						return eq2.almostEqualRelativeOrAbsolute(beta[jj], gradient);
					}, () -> {
						TestAssert.fail("Not same gradient @ %d,%d: %s != %s (error=%s)", ii, jj, beta[jj], gradient,
								DoubleEquality.relativeError(beta[jj], gradient));
					});
				}

			/////////////////////////////////////
			// This may be different
			/////////////////////////////////////
			final LVMGradientProcedure p12m3 = LVMGradientProcedureFactory.create(y_b, f12, type, fastLog);

			// Check these may be different.
			// Sometimes they are not different.

			p12m3.gradient(a2peaks);
			s = p12m3.value;
			System.arraycopy(p12m3.beta, 0, beta, 0, p12m3.beta.length);
			alpha = p12m3.getAlphaMatrix();

			//System.out.printf("%s [%d] p12m3  %f  %f\n", type, i, p123.value, s);

			// The test for different or equal is not robust to different random seeds.
			// TestAssert.fail has been changed for TestLog.logFailure

			if (type != Type.LSQ)
			{
				if (eq.almostEqualRelativeOrAbsolute(p123.value, s))
					TestLog.logFailure("p12b3 Same value @ %d (error=%s) : %s == %s\n", i,
							DoubleEquality.relativeError(p123.value, s), p123.value, s);
				if (eq.almostEqualRelativeOrAbsolute(beta, p123.beta))
					TestLog.logFailure("p12b3 Same gradient @ %d (error=%s) : %s vs %s\n", i,
							DoubleEquality.relativeError(beta, p123.beta), Arrays.toString(beta),
							Arrays.toString(p123.beta));

				// Note: Test the matrix is different by finding 1 different column
				int dj = -1;
				for (int j = 0; j < alpha.length; j++)
					//System.out.printf("%s !=\n%s\n\n", Arrays.toString(alpha[j]), Arrays.toString(m123[j]));
					if (!eq.almostEqualRelativeOrAbsolute(alpha[j], m123[j]))
					{
						dj = j; // Different column
						break;
					}
				if (dj == -1)
				{
					// Find biggest error for reporting. This helps set the test tolerance.
					double error = 0;
					dj = -1;
					for (int j = 0; j < alpha.length; j++)
					{
						final double e = DoubleEquality.relativeError(alpha[j], m123[j]);
						if (error <= e)
						{
							error = e;
							dj = j;
						}
					}
					TestLog.logFailure("p12b3 Same alpha @ %d,%d (error=%s) : %s vs %s\n", i, dj, error,
							Arrays.toString(alpha[dj]), Arrays.toString(m123[dj]));
				}
			}
			else
			{
				if (!eq.almostEqualRelativeOrAbsolute(p123.value, s))
					TestLog.logFailure("p12b3 Not same value @ %d (error=%s) : %s == %s\n", i,
							DoubleEquality.relativeError(p123.value, s), p123.value, s);
				if (!eq.almostEqualRelativeOrAbsolute(beta, p123.beta))
					TestLog.logFailure("p12b3 Not same gradient @ %d (error=%s) : %s vs %s\n", i,
							DoubleEquality.relativeError(beta, p123.beta), Arrays.toString(beta),
							Arrays.toString(p123.beta));
				for (int j = 0; j < alpha.length; j++)
					//System.out.printf("%s !=\n%s\n\n", Arrays.toString(alpha[j]), Arrays.toString(m123[j]));
					if (!eq.almostEqualRelativeOrAbsolute(alpha[j], m123[j]))
						TestLog.logFailure("p12b3 Not same alpha @ %d,%d (error=%s) : %s vs %s\n", i, j,
								DoubleEquality.relativeError(alpha[j], m123[j]), Arrays.toString(alpha[j]),
								Arrays.toString(m123[j]));
			}

			// Check actual gradients are correct
			if (!checkGradients)
				continue;

			for (int j = 0; j < nparams; j++)
			{
				final int jj = j;
				final int k = indices[j];
				//double d = Precision.representableDelta(a2peaks[k], (a2peaks[k] == 0) ? 1e-3 : a2peaks[k] * delta);
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
				//System.out.printf("[%d,%d] %f  (%s %f+/-%f)  %f  ?=  %f  (%f)\n", i, k, s,
				//		Gaussian2DFunction.getName(k), a2peaks[k], d, beta[j], gradient,
				//		DoubleEquality.relativeError(gradient, beta[j]));
				failCounter.run(nparams + j, () -> {
					return eq2.almostEqualRelativeOrAbsolute(beta[jj], gradient);
				}, () -> {
					TestAssert.fail("Not same gradient @ %d,%d: %s != %s (error=%s)", ii, jj, beta[jj], gradient,
							DoubleEquality.relativeError(beta[jj], gradient));
				});
			}
		}
	}

	/**
	 * Create random elliptical Gaussian data an returns the data plus an estimate of the parameters.
	 * Only the chosen parameters are randomised and returned for a maximum of (background, amplitude, angle, xpos,
	 * ypos, xwidth, ywidth }
	 *
	 * @param npeaks
	 *            the npeaks
	 * @param params
	 *            set on output
	 * @param randomiseParams
	 *            Set to true to randomise the params
	 * @return the double[]
	 */
	private double[] doubleCreateGaussianData(int npeaks, double[] params, boolean randomiseParams)
	{
		final int n = blockWidth * blockWidth;

		// Generate a 2D Gaussian
		final ErfGaussian2DFunction func = (ErfGaussian2DFunction) GaussianFunctionFactory.create2D(npeaks, blockWidth,
				blockWidth, GaussianFunctionFactory.FIT_ERF_FREE_CIRCLE, null);
		params[0] = random(Background);
		for (int i = 0, j = 0; i < npeaks; i++, j += Gaussian2DFunction.PARAMETERS_PER_PEAK)
		{
			params[j + Gaussian2DFunction.SIGNAL] = random(Signal);
			params[j + Gaussian2DFunction.X_POSITION] = random(Xpos);
			params[j + Gaussian2DFunction.Y_POSITION] = random(Ypos);
			params[j + Gaussian2DFunction.X_SD] = random(Xwidth);
			params[j + Gaussian2DFunction.Y_SD] = random(Ywidth);
		}

		if (npeaks > 1)
		{
			// Move the peaks around so they do not overlap
			final double[] shift = SimpleArrayUtils.newArray(npeaks, -2, 4.0 / (npeaks - 1));
			Random.shuffle(shift, rdg.getRandomGenerator());
			for (int i = 0, j = 0; i < npeaks; i++, j += Gaussian2DFunction.PARAMETERS_PER_PEAK)
				params[j + Gaussian2DFunction.X_POSITION] += shift[i];
			Random.shuffle(shift, rdg.getRandomGenerator());
			for (int i = 0, j = 0; i < npeaks; i++, j += Gaussian2DFunction.PARAMETERS_PER_PEAK)
				params[j + Gaussian2DFunction.Y_POSITION] += shift[i];
		}

		final double[] y = new double[n];
		func.initialise(params);
		for (int i = 0; i < y.length; i++)
			// Add random Poisson noise
			y[i] = rdg.nextPoisson(func.eval(i));

		if (randomiseParams)
		{
			params[0] = random(params[0]);
			for (int i = 0, j = 0; i < npeaks; i++, j += Gaussian2DFunction.PARAMETERS_PER_PEAK)
			{
				params[j + Gaussian2DFunction.SIGNAL] = random(params[j + Gaussian2DFunction.SIGNAL]);
				params[j + Gaussian2DFunction.X_POSITION] = random(params[j + Gaussian2DFunction.X_POSITION]);
				params[j + Gaussian2DFunction.Y_POSITION] = random(params[j + Gaussian2DFunction.Y_POSITION]);
				params[j + Gaussian2DFunction.X_SD] = random(params[j + Gaussian2DFunction.X_SD]);
				params[j + Gaussian2DFunction.Y_SD] = random(params[j + Gaussian2DFunction.Y_SD]); //params[j + 4];
			}
		}

		return y;
	}

	private double random(double d)
	{
		return d + rdg.nextUniform(-d * 0.1, d * 0.1);
	}

	protected int[] createData(int npeaks, int iter, ArrayList<double[]> paramsList, ArrayList<double[]> yList)
	{
		return createData(npeaks, iter, paramsList, yList, true);
	}

	protected int[] createData(int npeaks, int iter, ArrayList<double[]> paramsList, ArrayList<double[]> yList,
			boolean randomiseParams)
	{
		final int[] x = new int[blockWidth * blockWidth];
		for (int i = 0; i < x.length; i++)
			x[i] = i;
		for (int i = 0; i < iter; i++)
		{
			final double[] params = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK * npeaks];
			final double[] y = doubleCreateGaussianData(npeaks, params, randomiseParams);
			paramsList.add(params);
			yList.add(y);
		}
		return x;
	}

	protected int[] createFakeData(int nparams, int iter, ArrayList<double[]> paramsList, ArrayList<double[]> yList)
	{
		final int[] x = new int[blockWidth * blockWidth];
		for (int i = 0; i < x.length; i++)
			x[i] = i;
		for (int i = 0; i < iter; i++)
		{
			final double[] params = new double[nparams];
			final double[] y = createFakeData(params);
			paramsList.add(params);
			yList.add(y);
		}
		return x;
	}

	private double[] createFakeData(double[] params)
	{
		final int n = blockWidth * blockWidth;
		final RandomGenerator r = rdg.getRandomGenerator();

		for (int i = 0; i < params.length; i++)
			params[i] = r.nextDouble();

		final double[] y = new double[n];
		for (int i = 0; i < y.length; i++)
			y[i] = r.nextDouble();

		return y;
	}

	protected ArrayList<double[]> copyList(ArrayList<double[]> paramsList)
	{
		final ArrayList<double[]> params2List = new ArrayList<>(paramsList.size());
		for (int i = 0; i < paramsList.size(); i++)
			params2List.add(paramsList.get(i).clone());
		return params2List;
	}
}
