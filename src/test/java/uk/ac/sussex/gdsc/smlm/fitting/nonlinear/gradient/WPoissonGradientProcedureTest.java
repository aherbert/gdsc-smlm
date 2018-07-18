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
import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Assert;
import org.junit.Test;

import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.smlm.function.DummyGradientFunction;
import uk.ac.sussex.gdsc.smlm.function.FakeGradientFunction;
import uk.ac.sussex.gdsc.smlm.function.Gradient1Function;
import uk.ac.sussex.gdsc.test.TestLog;
import uk.ac.sussex.gdsc.test.TestSettings;

@SuppressWarnings({ "javadoc" })
public class WPoissonGradientProcedureTest
{
	DoubleEquality eq = new DoubleEquality(1e-6, 1e-16);

	int MAX_ITER = 20000;
	static int blockWidth = 10;
	double Background = 0.5;
	double Signal = 100;
	double Angle = Math.PI;
	double Xpos = 5;
	double Ypos = 5;
	double Xwidth = 1.2;
	double Ywidth = 1.2;

	RandomDataGenerator rdg;

	static double[] var;
	static
	{
		int n = blockWidth * blockWidth;
		var = new double[n];
		final RandomGenerator r = TestSettings.getRandomGenerator();
		while (n-- > 0)
			// Range 0.9 to 1.1
			var[n] = 0.9 + 0.2 * r.nextDouble();
	}

	@Test
	public void gradientProcedureFactoryCreatesOptimisedProcedures()
	{
		final double[] y = SimpleArrayUtils.newDoubleArray(var.length, 1);
		Assert.assertEquals(WPoissonGradientProcedureFactory.create(y, var, new DummyGradientFunction(6)).getClass(),
				WPoissonGradientProcedure6.class);
		Assert.assertEquals(WPoissonGradientProcedureFactory.create(y, var, new DummyGradientFunction(5)).getClass(),
				WPoissonGradientProcedure5.class);
		Assert.assertEquals(WPoissonGradientProcedureFactory.create(y, var, new DummyGradientFunction(4)).getClass(),
				WPoissonGradientProcedure4.class);
	}

	@Test
	public void poissonGradientProcedureComputesSameAsWLSQGradientProcedure()
	{
		poissonGradientProcedureComputesSameAsWLSQGradientProcedure(4);
		poissonGradientProcedureComputesSameAsWLSQGradientProcedure(5);
		poissonGradientProcedureComputesSameAsWLSQGradientProcedure(6);
		poissonGradientProcedureComputesSameAsWLSQGradientProcedure(11);
		poissonGradientProcedureComputesSameAsWLSQGradientProcedure(21);
	}

	private void poissonGradientProcedureComputesSameAsWLSQGradientProcedure(int nparams)
	{
		final int iter = 10;
		rdg = new RandomDataGenerator(TestSettings.getRandomGenerator());

		final ArrayList<double[]> paramsList = new ArrayList<>(iter);

		createFakeParams(nparams, iter, paramsList);
		final FakeGradientFunction func = new FakeGradientFunction(blockWidth, nparams);

		final String name = String.format("[%d]", nparams);

		for (int i = 0; i < paramsList.size(); i++)
		{
			final double[] y = createFakeData();
			final WPoissonGradientProcedure p = WPoissonGradientProcedureFactory.create(y, var, func);
			p.computeFisherInformation(paramsList.get(i));
			final WLSQLVMGradientProcedure p2 = new WLSQLVMGradientProcedure(y, var, func);
			p2.gradient(paramsList.get(i));

			// Exactly the same ...
			Assert.assertArrayEquals(name + " Observations: Not same alpha @ " + i, p.data, p2.alpha, 0);
			Assert.assertArrayEquals(name + " Observations: Not same alpha @ " + i, p.getLinear(), p2.getAlphaLinear(),
					0);
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
		final int iter = 10;
		rdg = new RandomDataGenerator(TestSettings.getRandomGenerator());

		final ArrayList<double[]> paramsList = new ArrayList<>(iter);

		createFakeParams(nparams, iter, paramsList);
		final Gradient1Function func = new FakeGradientFunction(blockWidth, nparams);

		final double[] v = (precomputed) ? var : null;

		final String name = String.format("[%d]", nparams);
		for (int i = 0; i < paramsList.size(); i++)
		{
			final double[] y = createFakeData();
			final WPoissonGradientProcedure p1 = new WPoissonGradientProcedure(y, v, func);
			p1.computeFisherInformation(paramsList.get(i));

			final WPoissonGradientProcedure p2 = WPoissonGradientProcedureFactory.create(y, v, func);
			p2.computeFisherInformation(paramsList.get(i));

			// Exactly the same ...
			Assert.assertArrayEquals(name + " Observations: Not same alpha @ " + i, p1.getLinear(), p2.getLinear(), 0);
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
		TestSettings.assumeSpeedTest();

		final int iter = 100;
		rdg = new RandomDataGenerator(TestSettings.getRandomGenerator());

		final ArrayList<double[]> paramsList = new ArrayList<>(iter);
		final ArrayList<double[]> yList = new ArrayList<>(iter);

		createFakeData(nparams, iter, paramsList, yList);

		// Remove the timing of the function call by creating a dummy function
		final FakeGradientFunction func = new FakeGradientFunction(blockWidth, nparams);
		final double[] v = (precomputed) ? var : null;

		for (int i = 0; i < paramsList.size(); i++)
		{
			final double[] y = yList.get(i);
			final WPoissonGradientProcedure p1 = new WPoissonGradientProcedure(y, v, func);
			p1.computeFisherInformation(paramsList.get(i));

			final WPoissonGradientProcedure p2 = WPoissonGradientProcedureFactory.create(y, v, func);
			p2.computeFisherInformation(paramsList.get(i));

			// Check they are the same
			Assert.assertArrayEquals("M " + i, p1.getLinear(), p2.getLinear(), 0);
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
					final WPoissonGradientProcedure p1 = new WPoissonGradientProcedure(yList.get(i), v, func);
					for (int j = loops; j-- > 0;)
						p1.computeFisherInformation(paramsList.get(k++ % iter));
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
					final WPoissonGradientProcedure p2 = WPoissonGradientProcedureFactory.create(yList.get(i), v, func);
					for (int j = loops; j-- > 0;)
						p2.computeFisherInformation(paramsList.get(k++ % iter));
				}
			}
		};
		final long time2 = t2.getTime();

		TestLog.logSpeedTestResult(time2 < time1, "Precomputed=%b : Standard %d : Unrolled %d = %d : %fx\n",
				precomputed, time1, nparams, time2, (1.0 * time1) / time2);
	}

	@Test
	public void gradientProcedureIsFasterThanWLSEGradientProcedure()
	{
		gradientProcedureIsFasterThanWLSEGradientProcedure(4);
		gradientProcedureIsFasterThanWLSEGradientProcedure(5);
		gradientProcedureIsFasterThanWLSEGradientProcedure(6);
		gradientProcedureIsFasterThanWLSEGradientProcedure(11);
	}

	private void gradientProcedureIsFasterThanWLSEGradientProcedure(final int nparams)
	{
		TestSettings.assumeSpeedTest();

		final int iter = 100;
		rdg = new RandomDataGenerator(TestSettings.getRandomGenerator());

		final ArrayList<double[]> paramsList = new ArrayList<>(iter);
		final ArrayList<double[]> yList = new ArrayList<>(iter);

		createFakeData(nparams, iter, paramsList, yList);

		// Remove the timing of the function call by creating a dummy function
		final FakeGradientFunction func = new FakeGradientFunction(blockWidth, nparams);
		for (int i = 0; i < paramsList.size(); i++)
		{
			final double[] y = yList.get(i);
			final WLSQLVMGradientProcedure p1 = WLSQLVMGradientProcedureFactory.create(y, var, func);
			p1.gradient(paramsList.get(i));

			final WPoissonGradientProcedure p2 = WPoissonGradientProcedureFactory.create(y, var, func);
			p2.computeFisherInformation(paramsList.get(i));

			// Check they are the same
			Assert.assertArrayEquals("M " + i, p1.getAlphaLinear(), p2.getLinear(), 0);
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
					final WLSQLVMGradientProcedure p1 = WLSQLVMGradientProcedureFactory.create(yList.get(i), var, func);
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
					final WPoissonGradientProcedure p2 = WPoissonGradientProcedureFactory.create(yList.get(i), var, func);
					for (int j = loops; j-- > 0;)
						p2.computeFisherInformation(paramsList.get(k++ % iter));
				}
			}
		};
		final long time2 = t2.getTime();

		TestLog.logSpeedTestResult(time2 < time1,
				"WLSQLVMGradientProcedure %d : WPoissonGradientProcedure %d = %d : %fx\n", time1, nparams, time2,
				(1.0 * time1) / time2);
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

	private double[] createFakeData()
	{
		final int n = blockWidth * blockWidth;
		final RandomGenerator r = rdg.getRandomGenerator();

		final double[] y = new double[n];
		for (int i = 0; i < y.length; i++)
			y[i] = r.nextDouble();

		return y;
	}

	protected void createFakeParams(int nparams, int iter, ArrayList<double[]> paramsList)
	{
		for (int i = 0; i < iter; i++)
		{
			final double[] params = new double[nparams];
			createFakeParams(params);
			paramsList.add(params);
		}
	}

	private void createFakeParams(double[] params)
	{
		final RandomGenerator r = rdg.getRandomGenerator();
		for (int i = 0; i < params.length; i++)
			params[i] = r.nextDouble();
	}

	protected ArrayList<double[]> copyList(ArrayList<double[]> paramsList)
	{
		final ArrayList<double[]> params2List = new ArrayList<>(paramsList.size());
		for (int i = 0; i < paramsList.size(); i++)
			params2List.add(copydouble(paramsList.get(i)));
		return params2List;
	}

	private static double[] copydouble(double[] d)
	{
		final double[] d2 = new double[d.length];
		for (int i = 0; i < d.length; i++)
			d2[i] = d[i];
		return d2;
	}
}