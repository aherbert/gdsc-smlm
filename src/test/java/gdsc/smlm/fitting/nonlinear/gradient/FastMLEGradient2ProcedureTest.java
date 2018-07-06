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
package gdsc.smlm.fitting.nonlinear.gradient;

import java.util.ArrayList;
import java.util.Arrays;

import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.util.Precision;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.utils.DoubleEquality;
import gdsc.smlm.function.DummyGradientFunction;
import gdsc.smlm.function.FakeGradientFunction;
import gdsc.smlm.function.Gradient2Function;
import gdsc.smlm.function.OffsetGradient2Function;
import gdsc.smlm.function.ValueProcedure;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import gdsc.smlm.function.gaussian.HoltzerAstigmatismZModel;
import gdsc.smlm.function.gaussian.erf.ErfGaussian2DFunction;
import gdsc.smlm.function.gaussian.erf.SingleAstigmatismErfGaussian2DFunction;
import gdsc.smlm.function.gaussian.erf.SingleFreeCircularErfGaussian2DFunction;
import gdsc.test.TestAssert;
import gdsc.test.TestSettings;

/**
 * Contains speed tests for the methods for calculating the Hessian and gradient vector
 * for use in the LVM algorithm.
 */
public class FastMLEGradient2ProcedureTest
{
	DoubleEquality eq = new DoubleEquality(1e-6, 1e-16);

	int MAX_ITER = 20000;
	int blockWidth = 10;
	double Background = 0.5;
	double Signal = 100;
	double Angle = Math.PI;
	double Xpos = 5;
	double Ypos = 5;
	double Xwidth = 1.2;
	double Ywidth = 1.2;

	RandomDataGenerator rdg;

	@Test
	public void gradientProcedureFactoryCreatesOptimisedProcedures()
	{
		double[] y = new double[0];
		Assert.assertEquals(FastMLEGradient2ProcedureFactory.createUnrolled(y, new DummyGradientFunction(4)).getClass(),
				FastMLEGradient2Procedure4.class);
		Assert.assertEquals(FastMLEGradient2ProcedureFactory.createUnrolled(y, new DummyGradientFunction(5)).getClass(),
				FastMLEGradient2Procedure5.class);
		Assert.assertEquals(FastMLEGradient2ProcedureFactory.createUnrolled(y, new DummyGradientFunction(6)).getClass(),
				FastMLEGradient2Procedure6.class);
	}

	@Test
	public void gradientProcedureComputesSameLogLikelihoodAsMLEGradientCalculator()
	{
		gradientProcedureComputesSameLogLikelihoodAsMLEGradientCalculator(4);
		gradientProcedureComputesSameLogLikelihoodAsMLEGradientCalculator(5);
		gradientProcedureComputesSameLogLikelihoodAsMLEGradientCalculator(6);
		gradientProcedureComputesSameLogLikelihoodAsMLEGradientCalculator(11);
		gradientProcedureComputesSameLogLikelihoodAsMLEGradientCalculator(21);
	}

	@Test
	public void gradientProcedureIsNotSlowerThanGradientCalculator()
	{
		// Note: The procedure does not have a lot of work within loops. It is only a single loop
		// so unrolling does not produce performance gains. The JVM can optimise this. 

		gradientProcedureIsNotSlowerThanGradientCalculator(4);
		gradientProcedureIsNotSlowerThanGradientCalculator(5);
		gradientProcedureIsNotSlowerThanGradientCalculator(6);
		gradientProcedureIsNotSlowerThanGradientCalculator(11);
		gradientProcedureIsNotSlowerThanGradientCalculator(21);
	}

	private void gradientProcedureComputesSameLogLikelihoodAsMLEGradientCalculator(int nparams)
	{
		int iter = 10;
		rdg = new RandomDataGenerator(TestSettings.getRandomGenerator());

		ArrayList<double[]> paramsList = new ArrayList<double[]>(iter);
		ArrayList<double[]> yList = new ArrayList<double[]>(iter);

		createFakeData(nparams, iter, paramsList, yList);
		FakeGradientFunction func = new FakeGradientFunction(blockWidth, nparams);

		MLEGradientCalculator calc = (MLEGradientCalculator) GradientCalculatorFactory.newCalculator(nparams, true);

		for (int i = 0; i < paramsList.size(); i++)
		{
			FastMLEGradient2Procedure p = FastMLEGradient2ProcedureFactory.createUnrolled(yList.get(i), func);
			double s = p.computeLogLikelihood(paramsList.get(i));
			double s2 = calc.logLikelihood(yList.get(i), paramsList.get(i), func);
			// Virtually the same ...
			TestAssert.assertEqualsRelative(s, s2, 1e-5, "[%d] Result: Not same @ %d", nparams, i);
		}
	}

	@Test
	public void gradientProcedureComputesSameWithPrecomputed()
	{
		int iter = 10;
		rdg = new RandomDataGenerator(TestSettings.getRandomGenerator());

		ErfGaussian2DFunction f1 = (ErfGaussian2DFunction) GaussianFunctionFactory.create2D(1, 10, 10,
				GaussianFunctionFactory.FIT_ERF_FREE_CIRCLE, null);
		ErfGaussian2DFunction f2 = (ErfGaussian2DFunction) GaussianFunctionFactory.create2D(2, 10, 10,
				GaussianFunctionFactory.FIT_ERF_FREE_CIRCLE, null);

		double[] a1 = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
		double[] a2 = new double[1 + 2 * Gaussian2DFunction.PARAMETERS_PER_PEAK];

		final double[] x = new double[f1.size()];
		final double[] b = new double[f1.size()];

		for (int i = 0; i < iter; i++)
		{
			a2[Gaussian2DFunction.BACKGROUND] = rdg.nextUniform(0.1, 0.3);
			a2[Gaussian2DFunction.SIGNAL] = rdg.nextUniform(100, 300);
			a2[Gaussian2DFunction.X_POSITION] = rdg.nextUniform(3, 5);
			a2[Gaussian2DFunction.Y_POSITION] = rdg.nextUniform(3, 5);
			a2[Gaussian2DFunction.Z_POSITION] = rdg.nextUniform(-2, 2);
			a2[Gaussian2DFunction.X_SD] = rdg.nextUniform(1, 1.3);
			a2[Gaussian2DFunction.Y_SD] = rdg.nextUniform(1, 1.3);
			a2[Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.SIGNAL] = rdg.nextUniform(100, 300);
			a2[Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.X_POSITION] = rdg.nextUniform(5, 7);
			a2[Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.Y_POSITION] = rdg.nextUniform(5, 7);
			a2[Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.Z_POSITION] = rdg.nextUniform(-3, 1);
			a2[Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.X_SD] = rdg.nextUniform(1, 1.3);
			a2[Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.Y_SD] = rdg.nextUniform(1, 1.3);

			// Simulate Poisson data
			f2.initialise0(a2);
			f1.forEach(new ValueProcedure()
			{
				int k = 0;

				@Override
				public void execute(double value)
				{
					x[k++] = (value > 0) ? rdg.nextPoisson(value) : 0;
				}
			});

			// Precompute peak 2 (no background)
			a1[Gaussian2DFunction.BACKGROUND] = 0;
			for (int j = 1; j < 7; j++)
				a1[j] = a2[Gaussian2DFunction.PARAMETERS_PER_PEAK + j];
			f1.initialise0(a1);
			f1.forEach(new ValueProcedure()
			{
				int k = 0;

				@Override
				public void execute(double value)
				{
					b[k++] = value;
				}
			});

			// Reset to peak 1
			for (int j = 0; j < 7; j++)
				a1[j] = a2[j];

			// Compute peak 1+2
			FastMLEGradient2Procedure p12 = FastMLEGradient2ProcedureFactory.create(x, f2);
			p12.computeSecondDerivative(a2);
			double[] d11 = Arrays.copyOf(p12.d1, f1.getNumberOfGradients());
			double[] d21 = Arrays.copyOf(p12.d2, f1.getNumberOfGradients());

			// Compute peak 1+(precomputed 2)
			FastMLEGradient2Procedure p1b2 = FastMLEGradient2ProcedureFactory.create(x,
					OffsetGradient2Function.wrapGradient2Function(f1, b));
			p1b2.computeSecondDerivative(a1);
			double[] d12 = p1b2.d1;
			double[] d22 = p1b2.d2;

			Assert.assertArrayEquals(" Result: Not same @ " + i, p12.u, p1b2.u, 1e-10);
			Assert.assertArrayEquals(" D1: Not same @ " + i, d11, d12, 1e-10);
			Assert.assertArrayEquals(" D2: Not same @ " + i, d21, d22, 1e-10);

			double[] v1 = p12.computeValue(a2);
			double[] v2 = p1b2.computeValue(a1);
			Assert.assertArrayEquals(" Value: Not same @ " + i, v1, v2, 1e-10);

			double[] d1 = Arrays.copyOf(p12.computeFirstDerivative(a2), f1.getNumberOfGradients());
			double[] d2 = p1b2.computeFirstDerivative(a1);
			Assert.assertArrayEquals(" 1st derivative: Not same @ " + i, d1, d2, 1e-10);
		}
	}

	private abstract class Timer
	{
		private int loops;
		int min = 5;

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
				long t2 = t1;
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

	private void gradientProcedureIsNotSlowerThanGradientCalculator(final int nparams)
	{
		TestSettings.assumeSpeedTest();

		final int iter = 1000;
		rdg = new RandomDataGenerator(TestSettings.getRandomGenerator());

		final ArrayList<double[]> paramsList = new ArrayList<double[]>(iter);
		final ArrayList<double[]> yList = new ArrayList<double[]>(iter);

		createFakeData(nparams, iter, paramsList, yList);
		final FakeGradientFunction func = new FakeGradientFunction(blockWidth, nparams);

		MLEGradientCalculator calc = (MLEGradientCalculator) GradientCalculatorFactory.newCalculator(nparams, true);

		for (int i = 0; i < paramsList.size(); i++)
			calc.logLikelihood(yList.get(i), paramsList.get(i), func);

		for (int i = 0; i < paramsList.size(); i++)
		{
			FastMLEGradient2Procedure p = FastMLEGradient2ProcedureFactory.createUnrolled(yList.get(i), func);
			p.computeLogLikelihood(paramsList.get(i));
		}

		// Realistic loops for an optimisation
		final int loops = 15;

		// Run till stable timing
		Timer t1 = new Timer()
		{
			@Override
			void run()
			{
				for (int i = 0, k = 0; i < iter; i++)
				{
					MLEGradientCalculator calc = (MLEGradientCalculator) GradientCalculatorFactory
							.newCalculator(nparams, true);
					for (int j = loops; j-- > 0;)
						calc.logLikelihood(yList.get(i), paramsList.get(k++ % iter), func);
				}
			}
		};
		long time1 = t1.getTime();

		Timer t2 = new Timer(t1.loops)
		{
			@Override
			void run()
			{
				for (int i = 0, k = 0; i < iter; i++)
				{
					FastMLEGradient2Procedure p = FastMLEGradient2ProcedureFactory.createUnrolled(yList.get(i), func);
					for (int j = loops; j-- > 0;)
						p.computeLogLikelihood(paramsList.get(k++ % iter));
				}
			}
		};
		long time2 = t2.getTime();

		TestSettings.logSpeedTestResult(time2 < time1,
				"GradientCalculator = %d : FastMLEGradient2Procedure %d = %d : %fx\n", time1, nparams, time2,
				(1.0 * time1) / time2);
	}

	@Test
	public void gradientProcedureUnrolledComputesSameAsGradientProcedure()
	{
		gradientProcedureUnrolledComputesSameAsGradientProcedure(4);
		gradientProcedureUnrolledComputesSameAsGradientProcedure(5);
		gradientProcedureUnrolledComputesSameAsGradientProcedure(6);
	}

	private void gradientProcedureUnrolledComputesSameAsGradientProcedure(int nparams)
	{
		int iter = 10;
		rdg = new RandomDataGenerator(TestSettings.getRandomGenerator());

		ArrayList<double[]> paramsList = new ArrayList<double[]>(iter);
		ArrayList<double[]> yList = new ArrayList<double[]>(iter);

		createFakeData(nparams, iter, paramsList, yList);
		FakeGradientFunction func = new FakeGradientFunction(blockWidth, nparams);

		FastMLEGradient2Procedure p1, p2;
		for (int i = 0; i < paramsList.size(); i++)
		{
			p1 = new FastMLEGradient2Procedure(yList.get(i), func);
			p2 = FastMLEGradient2ProcedureFactory.createUnrolled(yList.get(i), func);
			double[] a = paramsList.get(i);

			double ll1 = p1.computeLogLikelihood(a);
			double ll2 = p2.computeLogLikelihood(a);
			TestAssert.assertEquals(ll1, ll2, 0, "[%d] LL: Not same @ %d", nparams, i);

			p1 = new FastMLEGradient2Procedure(yList.get(i), func);
			p2 = FastMLEGradient2ProcedureFactory.createUnrolled(yList.get(i), func);
			p1.computeFirstDerivative(a);
			p2.computeFirstDerivative(a);
			TestAssert.assertArrayEquals(p1.u, p2.u, 0, "[%d]  first derivative value: Not same @ %d", nparams, i);
			TestAssert.assertArrayEquals(p1.d1, p2.d1, 0, "[%d]  first derivative: Not same @ %d", nparams, i);

			p1 = new FastMLEGradient2Procedure(yList.get(i), func);
			p2 = FastMLEGradient2ProcedureFactory.createUnrolled(yList.get(i), func);
			p1.computeSecondDerivative(a);
			p2.computeSecondDerivative(a);
			TestAssert.assertArrayEquals(p1.u, p2.u, 0, "[%d]  update value: Not same @ %d", nparams, i);
			TestAssert.assertArrayEquals(p1.d1, p2.d1, 0, "[%d]  update: Not same d1 @ %d", nparams, i);
			TestAssert.assertArrayEquals(p1.d2, p2.d2, 0, "[%d]  update: Not same d2 @ %d", nparams, i);
		}
	}

	@Test
	public void gradientProcedureIsFasterUnrolledThanGradientProcedure()
	{
		gradientProcedureLinearIsFasterThanGradientProcedure(4);
		gradientProcedureLinearIsFasterThanGradientProcedure(5);
		gradientProcedureLinearIsFasterThanGradientProcedure(6);
	}

	private void gradientProcedureLinearIsFasterThanGradientProcedure(final int nparams)
	{
		TestSettings.assumeMediumComplexity();

		final int iter = 100;
		rdg = new RandomDataGenerator(TestSettings.getRandomGenerator());

		final ArrayList<double[]> paramsList = new ArrayList<double[]>(iter);
		final ArrayList<double[]> yList = new ArrayList<double[]>(iter);

		createData(1, iter, paramsList, yList);

		// Remove the timing of the function call by creating a dummy function
		final Gradient2Function func = new FakeGradientFunction(blockWidth, nparams);

		for (int i = 0; i < paramsList.size(); i++)
		{
			FastMLEGradient2Procedure p1 = new FastMLEGradient2Procedure(yList.get(i), func);
			p1.computeSecondDerivative(paramsList.get(i));
			p1.computeSecondDerivative(paramsList.get(i));

			FastMLEGradient2Procedure p2 = FastMLEGradient2ProcedureFactory.createUnrolled(yList.get(i), func);
			p2.computeSecondDerivative(paramsList.get(i));
			p2.computeSecondDerivative(paramsList.get(i));

			// Check they are the same
			Assert.assertArrayEquals("D1 " + i, p1.d1, p2.d1, 0);
			Assert.assertArrayEquals("D2 " + i, p1.d2, p2.d2, 0);
		}

		// Realistic loops for an optimisation
		final int loops = 15;

		// Run till stable timing
		Timer t1 = new Timer()
		{
			@Override
			void run()
			{
				for (int i = 0, k = 0; i < paramsList.size(); i++)
				{
					FastMLEGradient2Procedure p1 = new FastMLEGradient2Procedure(yList.get(i), func);
					for (int j = loops; j-- > 0;)
						p1.computeSecondDerivative(paramsList.get(k++ % iter));
				}
			}
		};
		long time1 = t1.getTime();

		Timer t2 = new Timer(t1.loops)
		{
			@Override
			void run()
			{
				for (int i = 0, k = 0; i < paramsList.size(); i++)
				{
					FastMLEGradient2Procedure p2 = FastMLEGradient2ProcedureFactory.createUnrolled(yList.get(i), func);
					for (int j = loops; j-- > 0;)
						p2.computeSecondDerivative(paramsList.get(k++ % iter));
				}
			}
		};
		long time2 = t2.getTime();

		TestSettings.info("Standard = %d : Unrolled %d = %d : %fx\n", time1, nparams, time2, (1.0 * time1) / time2);
		Assert.assertTrue(time2 < time1 * 1.5);
	}

	@Test
	public void gradientCalculatorComputesGradient()
	{
		gradientCalculatorComputesGradient(new SingleFreeCircularErfGaussian2DFunction(blockWidth, blockWidth));

		// Use a reasonable z-depth function from the Smith, et al (2010) paper (page 377)
		double sx = 1.08;
		double sy = 1.01;
		double gamma = 0.389;
		double d = 0.531;
		double Ax = -0.0708;
		double Bx = -0.073;
		double Ay = 0.164;
		double By = 0.0417;
		HoltzerAstigmatismZModel zModel = HoltzerAstigmatismZModel.create(sx, sy, gamma, d, Ax, Bx, Ay, By);
		gradientCalculatorComputesGradient(new SingleAstigmatismErfGaussian2DFunction(blockWidth, blockWidth, zModel));
	}

	private void gradientCalculatorComputesGradient(ErfGaussian2DFunction func)
	{
		// Check the first and second derivatives
		int nparams = func.getNumberOfGradients();
		int[] indices = func.gradientIndices();

		int iter = 100;
		rdg = new RandomDataGenerator(TestSettings.getRandomGenerator());

		ArrayList<double[]> paramsList = new ArrayList<double[]>(iter);
		ArrayList<double[]> yList = new ArrayList<double[]>(iter);

		createData(1, iter, paramsList, yList, true);

		// for the gradients
		double delta = 1e-4;
		DoubleEquality eq = new DoubleEquality(5e-2, 1e-16);

		for (int i = 0; i < paramsList.size(); i++)
		{
			double[] y = yList.get(i);
			double[] a = paramsList.get(i);
			double[] a2 = a.clone();
			FastMLEGradient2Procedure p = FastMLEGradient2ProcedureFactory.create(y, func);
			//double ll = p.computeLogLikelihood(a);
			p.computeSecondDerivative(a);
			double[] d1 = p.d1.clone();
			double[] d2 = p.d2.clone();
			for (int j = 0; j < nparams; j++)
			{
				int k = indices[j];
				double d = Precision.representableDelta(a[k], (a[k] == 0) ? delta : a[k] * delta);
				a2[k] = a[k] + d;
				double llh = p.computeLogLikelihood(a2);
				p.computeFirstDerivative(a2);
				double[] d1h = p.d1.clone();
				a2[k] = a[k] - d;
				double lll = p.computeLogLikelihood(a2);
				p.computeFirstDerivative(a2);
				double[] d1l = p.d1.clone();
				a2[k] = a[k];

				double gradient1 = (llh - lll) / (2 * d);
				double gradient2 = (d1h[j] - d1l[j]) / (2 * d);
				//System.out.printf("[%d,%d] ll - %f  (%s %f+/-%f) d1 %f ?= %f : d2 %f ?= %f\n", i, k, ll, func.getName(k), a[k], d, 
				//		gradient1, d1[j], gradient2, d2[j]);
				Assert.assertTrue("Not same gradient1 @ " + j, eq.almostEqualRelativeOrAbsolute(gradient1, d1[j]));
				if (!eq.almostEqualRelativeOrAbsolute(gradient2, d2[j]))
				{
					Assert.fail(
							"Not same gradient2 @ " + j + " error = " + DoubleEquality.relativeError(gradient2, d2[j]));
				}
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
		int n = blockWidth * blockWidth;

		// Generate a 2D Gaussian
		ErfGaussian2DFunction func = (ErfGaussian2DFunction) GaussianFunctionFactory.create2D(npeaks, blockWidth,
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

		double[] y = new double[n];
		func.initialise(params);
		for (int i = 0; i < y.length; i++)
		{
			// Add random Poisson noise
			y[i] = rdg.nextPoisson(func.eval(i));
		}

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
		int[] x = new int[blockWidth * blockWidth];
		for (int i = 0; i < x.length; i++)
			x[i] = i;
		for (int i = 0; i < iter; i++)
		{
			double[] params = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK * npeaks];
			double[] y = doubleCreateGaussianData(npeaks, params, randomiseParams);
			paramsList.add(params);
			yList.add(y);
		}
		return x;
	}

	protected int[] createFakeData(int nparams, int iter, ArrayList<double[]> paramsList, ArrayList<double[]> yList)
	{
		int[] x = new int[blockWidth * blockWidth];
		for (int i = 0; i < x.length; i++)
			x[i] = i;
		for (int i = 0; i < iter; i++)
		{
			double[] params = new double[nparams];
			double[] y = createFakeData(params);
			paramsList.add(params);
			yList.add(y);
		}
		return x;
	}

	private double[] createFakeData(double[] params)
	{
		int n = blockWidth * blockWidth;
		RandomGenerator r = rdg.getRandomGenerator();

		for (int i = 0; i < params.length; i++)
		{
			params[i] = r.nextDouble();
		}

		double[] y = new double[n];
		for (int i = 0; i < y.length; i++)
		{
			y[i] = r.nextDouble() * 10;
		}

		return y;
	}

	protected ArrayList<double[]> copyList(ArrayList<double[]> paramsList)
	{
		ArrayList<double[]> params2List = new ArrayList<double[]>(paramsList.size());
		for (int i = 0; i < paramsList.size(); i++)
		{
			params2List.add(copydouble(paramsList.get(i)));
		}
		return params2List;
	}

	private double[] copydouble(double[] d)
	{
		double[] d2 = new double[d.length];
		for (int i = 0; i < d.length; i++)
			d2[i] = d[i];
		return d2;
	}
}
