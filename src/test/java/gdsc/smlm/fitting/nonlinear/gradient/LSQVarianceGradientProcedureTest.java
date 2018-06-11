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

import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.utils.DoubleEquality;
import gdsc.core.utils.SimpleArrayUtils;
import gdsc.smlm.TestSettings;
import gdsc.smlm.function.DummyGradientFunction;
import gdsc.smlm.function.FakeGradientFunction;
import gdsc.smlm.function.Gradient1Function;
import gdsc.smlm.function.OffsetGradient1Function;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import gdsc.smlm.function.gaussian.erf.ErfGaussian2DFunction;
import gdsc.smlm.results.Gaussian2DPeakResultHelper;

public class LSQVarianceGradientProcedureTest
{
	boolean speedTests = true;
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
		Assert.assertEquals(LSQVarianceGradientProcedureFactory.create(new DummyGradientFunction(6)).getClass(),
				LSQVarianceGradientProcedure6.class);
		Assert.assertEquals(LSQVarianceGradientProcedureFactory.create(new DummyGradientFunction(5)).getClass(),
				LSQVarianceGradientProcedure5.class);
		Assert.assertEquals(LSQVarianceGradientProcedureFactory.create(new DummyGradientFunction(4)).getClass(),
				LSQVarianceGradientProcedure4.class);
	}

	@Test
	public void gradientProcedureComputesSameAsGradientCalculator()
	{
		gradientProcedureComputesSameAsGradientCalculator(4);
		gradientProcedureComputesSameAsGradientCalculator(5);
		gradientProcedureComputesSameAsGradientCalculator(6);
		gradientProcedureComputesSameAsGradientCalculator(11);
		gradientProcedureComputesSameAsGradientCalculator(21);
	}

	@Test
	public void gradientProcedureIsNotSlowerThanGradientCalculator()
	{
		gradientProcedureIsNotSlowerThanGradientCalculator(4);
		gradientProcedureIsNotSlowerThanGradientCalculator(5);
		gradientProcedureIsNotSlowerThanGradientCalculator(6);
		// 2 peaks
		gradientProcedureIsNotSlowerThanGradientCalculator(11);
		// 4 peaks
		gradientProcedureIsNotSlowerThanGradientCalculator(21);
	}

	private void gradientProcedureComputesSameAsGradientCalculator(int nparams)
	{
		int iter = 10;
		rdg = new RandomDataGenerator(new Well19937c(30051977));

		ArrayList<double[]> paramsList = new ArrayList<double[]>(iter);

		createFakeParams(nparams, iter, paramsList);
		int n = blockWidth * blockWidth;
		FakeGradientFunction func = new FakeGradientFunction(blockWidth, nparams);

		GradientCalculator calc = GradientCalculatorFactory.newCalculator(nparams, false);

		String name = String.format("[%d]", nparams);

		for (int i = 0; i < paramsList.size(); i++)
		{
			LSQVarianceGradientProcedure p = LSQVarianceGradientProcedureFactory.create(func);
			p.variance(paramsList.get(i));
			double[] e = calc.variance(n, paramsList.get(i), func);
			Assert.assertArrayEquals(name + " Observations: Not same @ " + i, e, p.variance, 0);
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
		org.junit.Assume.assumeTrue(speedTests || TestSettings.RUN_SPEED_TESTS);

		final int iter = 1000;
		rdg = new RandomDataGenerator(new Well19937c(30051977));

		final ArrayList<double[]> paramsList = new ArrayList<double[]>(iter);

		createFakeParams(nparams, iter, paramsList);
		final int n = blockWidth * blockWidth;
		final FakeGradientFunction func = new FakeGradientFunction(blockWidth, nparams);

		GradientCalculator calc = GradientCalculatorFactory.newCalculator(nparams, false);

		for (int i = 0; i < paramsList.size(); i++)
			calc.variance(n, paramsList.get(i), func);

		for (int i = 0; i < paramsList.size(); i++)
		{
			LSQVarianceGradientProcedure p = LSQVarianceGradientProcedureFactory.create(func);
			p.variance(paramsList.get(i));
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
					GradientCalculator calc = GradientCalculatorFactory.newCalculator(nparams, false);
					for (int j = loops; j-- > 0;)
						calc.variance(n, paramsList.get(k++ % iter), func);
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
					LSQVarianceGradientProcedure p = LSQVarianceGradientProcedureFactory.create(func);
					for (int j = loops; j-- > 0;)
						p.variance(paramsList.get(k++ % iter));
				}
			}
		};
		long time2 = t2.getTime();

		log("GradientCalculator = %d : LSQVarianceGradientProcedure %d = %d : %fx\n", time1, nparams, time2,
				(1.0 * time1) / time2);
		if (TestSettings.ASSERT_SPEED_TESTS)
		{
			// Add contingency
			Assert.assertTrue(time2 < time1 * 1.5);
		}
	}

	@Test
	public void gradientProcedureUnrolledComputesSameAsGradientProcedure()
	{
		gradientProcedureUnrolledComputesSameAsGradientProcedure(4, false);
		gradientProcedureUnrolledComputesSameAsGradientProcedure(5, false);
		gradientProcedureUnrolledComputesSameAsGradientProcedure(6, false);
	}

	@Test
	public void gradientProcedureUnrolledComputesSameAsGradientProcedureWithPrecomputed()
	{
		gradientProcedureUnrolledComputesSameAsGradientProcedure(4, true);
		gradientProcedureUnrolledComputesSameAsGradientProcedure(5, true);
		gradientProcedureUnrolledComputesSameAsGradientProcedure(6, true);
	}

	private void gradientProcedureUnrolledComputesSameAsGradientProcedure(int nparams, boolean precomputed)
	{
		int iter = 10;
		rdg = new RandomDataGenerator(new Well19937c(30051977));

		ArrayList<double[]> paramsList = new ArrayList<double[]>(iter);

		createFakeParams(nparams, iter, paramsList);
		Gradient1Function func = new FakeGradientFunction(blockWidth, nparams);

		if (precomputed)
			func = OffsetGradient1Function.wrapGradient1Function(func,
					SimpleArrayUtils.newArray(func.size(), 0.1, 1.3));

		String name = String.format("[%d]", nparams);
		for (int i = 0; i < paramsList.size(); i++)
		{
			LSQVarianceGradientProcedure p1 = new LSQVarianceGradientProcedure(func);
			p1.variance(paramsList.get(i));

			LSQVarianceGradientProcedure p2 = LSQVarianceGradientProcedureFactory.create(func);
			p2.variance(paramsList.get(i));

			// Exactly the same ...
			Assert.assertArrayEquals(name + " Observations: Not same alpha @ " + i, p1.variance, p2.variance, 0);
		}
	}

	@Test
	public void gradientProcedureIsFasterUnrolledThanGradientProcedure()
	{
		gradientProcedureIsFasterUnrolledThanGradientProcedure(4, false);
		gradientProcedureIsFasterUnrolledThanGradientProcedure(5, false);
		gradientProcedureIsFasterUnrolledThanGradientProcedure(6, false);
	}

	@Test
	public void gradientProcedureIsFasterUnrolledThanGradientProcedureWithPrecomputed()
	{
		gradientProcedureIsFasterUnrolledThanGradientProcedure(4, true);
		gradientProcedureIsFasterUnrolledThanGradientProcedure(5, true);
		gradientProcedureIsFasterUnrolledThanGradientProcedure(6, true);
	}

	private void gradientProcedureIsFasterUnrolledThanGradientProcedure(final int nparams, final boolean precomputed)
	{
		org.junit.Assume.assumeTrue(speedTests || TestSettings.RUN_SPEED_TESTS);

		final int iter = 100;
		rdg = new RandomDataGenerator(new Well19937c(30051977));

		final ArrayList<double[]> paramsList = new ArrayList<double[]>(iter);

		createFakeParams(nparams, iter, paramsList);

		// Remove the timing of the function call by creating a dummy function
		FakeGradientFunction f = new FakeGradientFunction(blockWidth, nparams);
		final Gradient1Function func = (precomputed)
				? OffsetGradient1Function.wrapGradient1Function(f, SimpleArrayUtils.newArray(f.size(), 0.1, 1.3)) : f;

		for (int i = 0; i < paramsList.size(); i++)
		{
			LSQVarianceGradientProcedure p1 = new LSQVarianceGradientProcedure(func);
			p1.variance(paramsList.get(i));
			p1.variance(paramsList.get(i));

			LSQVarianceGradientProcedure p2 = LSQVarianceGradientProcedureFactory.create(func);
			p2.variance(paramsList.get(i));
			p2.variance(paramsList.get(i));

			// Check they are the same
			Assert.assertArrayEquals("M " + i, p1.variance, p2.variance, 0);
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
					LSQVarianceGradientProcedure p1 = new LSQVarianceGradientProcedure(func);
					for (int j = loops; j-- > 0;)
						p1.variance(paramsList.get(k++ % iter));
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
					LSQVarianceGradientProcedure p2 = LSQVarianceGradientProcedureFactory.create(func);
					for (int j = loops; j-- > 0;)
						p2.variance(paramsList.get(k++ % iter));
				}
			}
		};
		long time2 = t2.getTime();

		log("Precomputed=%b : Standard %d : Unrolled %d = %d : %fx\n", precomputed, time1, nparams, time2,
				(1.0 * time1) / time2);
		Assert.assertTrue(time2 < time1 * 1.1); // Allow margin of error
	}

	@Test
	public void crlbIsHigherWithPrecomputed()
	{
		int iter = 10;
		rdg = new RandomDataGenerator(new Well19937c(30051977));

		ErfGaussian2DFunction func = (ErfGaussian2DFunction) GaussianFunctionFactory.create2D(1, 10, 10,
				GaussianFunctionFactory.FIT_ERF_FREE_CIRCLE, null);

		double[] a = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
		int n = func.getNumberOfGradients();

		// Get a background
		double[] b = new double[func.size()];
		for (int i = 0; i < b.length; i++)
			b[i] = rdg.nextUniform(1, 2);

		for (int i = 0; i < iter; i++)
		{
			a[Gaussian2DFunction.BACKGROUND] = rdg.nextUniform(0.1, 0.3);
			a[Gaussian2DFunction.SIGNAL] = rdg.nextUniform(100, 300);
			a[Gaussian2DFunction.X_POSITION] = rdg.nextUniform(4, 6);
			a[Gaussian2DFunction.Y_POSITION] = rdg.nextUniform(4, 6);
			a[Gaussian2DFunction.X_SD] = rdg.nextUniform(1, 1.3);
			a[Gaussian2DFunction.Y_SD] = rdg.nextUniform(1, 1.3);

			LSQVarianceGradientProcedure p1 = LSQVarianceGradientProcedureFactory.create(func);
			p1.variance(a);

			LSQVarianceGradientProcedure p2 = LSQVarianceGradientProcedureFactory
					.create(OffsetGradient1Function.wrapGradient1Function(func, b));
			p2.variance(a);

			double[] crlb1 = p1.variance;
			double[] crlb2 = p2.variance;
			Assert.assertNotNull(crlb1);
			Assert.assertNotNull(crlb2);
			//System.out.printf("%s : %s\n", Arrays.toString(crlb1), Arrays.toString(crlb2));
			for (int j = 0; j < n; j++)
				Assert.assertTrue(crlb1[j] < crlb2[j]);
		}
	}

	@Test
	public void varianceMatchesFormula()
	{
		//Assume.assumeTrue(false);

		double[] N_ = new double[] { 20, 50, 100, 500 };
		double[] b2_ = new double[] { 0, 1, 2, 4 };
		double[] s_ = new double[] { 1, 1.2, 1.5 };
		double[] x_ = new double[] { 4.8, 5, 5.5 };
		double a = 100;
		int size = 10;
		Gaussian2DFunction f = GaussianFunctionFactory.create2D(1, size, size, GaussianFunctionFactory.FIT_ERF_CIRCLE,
				null);
		LSQVarianceGradientProcedure p = LSQVarianceGradientProcedureFactory.create(f);
		int ix = f.findGradientIndex(Gaussian2DFunction.X_POSITION);
		int iy = f.findGradientIndex(Gaussian2DFunction.Y_POSITION);
		double[] params = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
		for (double N : N_)
		{
			params[Gaussian2DFunction.SIGNAL] = N;
			for (double b2 : b2_)
			{
				params[Gaussian2DFunction.BACKGROUND] = b2;
				for (double s : s_)
				{
					double ss = s * a;
					params[Gaussian2DFunction.X_SD] = s;
					for (double x : x_)
					{
						params[Gaussian2DFunction.X_POSITION] = x;
						for (double y : x_)
						{
							params[Gaussian2DFunction.Y_POSITION] = y;
							if (p.variance(params) != LSQVarianceGradientProcedure.STATUS_OK)
								Assert.fail("No variance");
							double o1 = Math.sqrt(p.variance[ix]) * a;
							double o2 = Math.sqrt(p.variance[iy]) * a;
							double e = Gaussian2DPeakResultHelper.getPrecisionX(a, ss, N, b2, false);
							//System.out.printf("e = %f  :  o  =   %f   %f\n", e, o1, o2);
							Assert.assertEquals(e, o1, e * 5e-2);
							Assert.assertEquals(e, o2, e * 5e-2);
						}
					}
				}
			}
		}
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
			y[i] = r.nextDouble();
		}

		return y;
	}

	protected void createFakeParams(int nparams, int iter, ArrayList<double[]> paramsList)
	{
		for (int i = 0; i < iter; i++)
		{
			double[] params = new double[nparams];
			createFakeParams(params);
			paramsList.add(params);
		}
	}

	private void createFakeParams(double[] params)
	{
		RandomGenerator r = rdg.getRandomGenerator();
		for (int i = 0; i < params.length; i++)
		{
			params[i] = r.nextDouble();
		}
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

	void log(String format, Object... args)
	{
		System.out.printf(format, args);
	}
}
