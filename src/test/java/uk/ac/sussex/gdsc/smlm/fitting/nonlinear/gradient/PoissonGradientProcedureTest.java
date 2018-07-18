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
import org.ejml.data.DenseMatrix64F;
import org.junit.Assert;
import org.junit.Test;

import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.smlm.fitting.FisherInformationMatrix;
import uk.ac.sussex.gdsc.smlm.function.DummyGradientFunction;
import uk.ac.sussex.gdsc.smlm.function.FakeGradientFunction;
import uk.ac.sussex.gdsc.smlm.function.Gradient1Function;
import uk.ac.sussex.gdsc.smlm.function.OffsetGradient1Function;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import uk.ac.sussex.gdsc.smlm.function.gaussian.erf.ErfGaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultHelper;
import uk.ac.sussex.gdsc.test.TestAssert;
import uk.ac.sussex.gdsc.test.TestLog;
import uk.ac.sussex.gdsc.test.TestSettings;

@SuppressWarnings({ "javadoc" })
public class PoissonGradientProcedureTest
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

	@Test
	public void gradientProcedureFactoryCreatesOptimisedProcedures()
	{
		Assert.assertEquals(PoissonGradientProcedureFactory.create(new DummyGradientFunction(6)).getClass(),
				PoissonGradientProcedure6.class);
		Assert.assertEquals(PoissonGradientProcedureFactory.create(new DummyGradientFunction(5)).getClass(),
				PoissonGradientProcedure5.class);
		Assert.assertEquals(PoissonGradientProcedureFactory.create(new DummyGradientFunction(4)).getClass(),
				PoissonGradientProcedure4.class);
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
		final int iter = 10;
		final ArrayList<double[]> paramsList = new ArrayList<>(iter);

		createFakeParams(TestSettings.getRandomGenerator(), nparams, iter, paramsList);
		final int n = blockWidth * blockWidth;
		final FakeGradientFunction func = new FakeGradientFunction(blockWidth, nparams);

		final GradientCalculator calc = GradientCalculatorFactory.newCalculator(nparams, false);

		for (int i = 0; i < paramsList.size(); i++)
		{
			final PoissonGradientProcedure p = PoissonGradientProcedureFactory.create(func);
			p.computeFisherInformation(paramsList.get(i));
			final double[][] m = calc.fisherInformationMatrix(n, paramsList.get(i), func);
			// Not exactly the same ...
			final double[] al = p.getLinear();
			TestAssert.assertArrayEqualsRelative(al, new DenseMatrix64F(m).data, 1e-10,
					"[%d] Observations: Not same alphaLinear @ %d", nparams, i);

			final double[][] am = p.getMatrix();
			for (int j = 0; j < nparams; j++)
				TestAssert.assertArrayEqualsRelative(am[j], m[j], 1e-10,
						"[%d] Observations: Not same alphaMatrix @ %d,%d", nparams, i, j);
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

	private void gradientProcedureIsNotSlowerThanGradientCalculator(final int nparams)
	{
		TestSettings.assumeSpeedTest();

		final int iter = 1000;
		final ArrayList<double[]> paramsList = new ArrayList<>(iter);

		createFakeParams(TestSettings.getRandomGenerator(), nparams, iter, paramsList);
		final int n = blockWidth * blockWidth;
		final FakeGradientFunction func = new FakeGradientFunction(blockWidth, nparams);

		final GradientCalculator calc = GradientCalculatorFactory.newCalculator(nparams, false);

		for (int i = 0; i < paramsList.size(); i++)
			calc.fisherInformationMatrix(n, paramsList.get(i), func);

		for (int i = 0; i < paramsList.size(); i++)
		{
			final PoissonGradientProcedure p = PoissonGradientProcedureFactory.create(func);
			p.computeFisherInformation(paramsList.get(i));
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
					final GradientCalculator calc = GradientCalculatorFactory.newCalculator(nparams, false);
					for (int j = loops; j-- > 0;)
						calc.fisherInformationMatrix(n, paramsList.get(k++ % iter), func);
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
					final PoissonGradientProcedure p = PoissonGradientProcedureFactory.create(func);
					for (int j = loops; j-- > 0;)
						p.computeFisherInformation(paramsList.get(k++ % iter));
				}
			}
		};
		final long time2 = t2.getTime();

		TestLog.logSpeedTestResult(time2 < time1,
				"GradientCalculator = %d : PoissonGradientProcedure %d = %d : %fx\n", time1, nparams, time2,
				(1.0 * time1) / time2);
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
		final ArrayList<double[]> paramsList = new ArrayList<>(iter);

		createFakeParams(TestSettings.getRandomGenerator(), nparams, iter, paramsList);
		Gradient1Function func = new FakeGradientFunction(blockWidth, nparams);

		if (precomputed)
			func = OffsetGradient1Function.wrapGradient1Function(func,
					SimpleArrayUtils.newArray(func.size(), 0.1, 1.3));

		final String name = String.format("[%d]", nparams);
		for (int i = 0; i < paramsList.size(); i++)
		{
			final PoissonGradientProcedure p1 = new PoissonGradientProcedure(func);
			p1.computeFisherInformation(paramsList.get(i));

			final PoissonGradientProcedure p2 = PoissonGradientProcedureFactory.create(func);
			p2.computeFisherInformation(paramsList.get(i));

			// Exactly the same ...
			Assert.assertArrayEquals(name + " Observations: Not same alpha @ " + i, p1.getLinear(), p2.getLinear(), 0);

			final double[][] am1 = p1.getMatrix();
			final double[][] am2 = p2.getMatrix();
			for (int j = 0; j < nparams; j++)
				Assert.assertArrayEquals(name + " Observations: Not same alpha @ " + i, am1[j], am2[j], 0);
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
		final ArrayList<double[]> paramsList = new ArrayList<>(iter);

		createFakeParams(TestSettings.getRandomGenerator(), nparams, iter, paramsList);

		// Remove the timing of the function call by creating a dummy function
		final FakeGradientFunction f = new FakeGradientFunction(blockWidth, nparams);
		final Gradient1Function func = (precomputed)
				? OffsetGradient1Function.wrapGradient1Function(f, SimpleArrayUtils.newArray(f.size(), 0.1, 1.3))
				: f;

		for (int i = 0; i < paramsList.size(); i++)
		{
			final PoissonGradientProcedure p1 = new PoissonGradientProcedure(func);
			p1.computeFisherInformation(paramsList.get(i));
			p1.computeFisherInformation(paramsList.get(i));

			final PoissonGradientProcedure p2 = PoissonGradientProcedureFactory.create(func);
			p2.computeFisherInformation(paramsList.get(i));
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
					final PoissonGradientProcedure p1 = new PoissonGradientProcedure(func);
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
					final PoissonGradientProcedure p2 = PoissonGradientProcedureFactory.create(func);
					for (int j = loops; j-- > 0;)
						p2.computeFisherInformation(paramsList.get(k++ % iter));
				}
			}
		};
		final long time2 = t2.getTime();

		TestLog.logSpeedTestResult(time2 < time1, "Precomputed=%b : Standard %d : Unrolled %d = %d : %fx\n",
				precomputed, time1, nparams, time2, (1.0 * time1) / time2);
		//Assert.assertTrue(time2 < time1);
	}

	@Test
	public void crlbIsHigherWithPrecomputed()
	{
		final int iter = 10;
		final RandomDataGenerator rdg = new RandomDataGenerator(TestSettings.getRandomGenerator());

		final ErfGaussian2DFunction func = (ErfGaussian2DFunction) GaussianFunctionFactory.create2D(1, 10, 10,
				GaussianFunctionFactory.FIT_ERF_FREE_CIRCLE, null);

		final double[] a = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
		final int n = func.getNumberOfGradients();

		// Get a background
		final double[] b = new double[func.size()];
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

			final PoissonGradientProcedure p1 = PoissonGradientProcedureFactory.create(func);
			p1.computeFisherInformation(a);

			final PoissonGradientProcedure p2 = PoissonGradientProcedureFactory
					.create(OffsetGradient1Function.wrapGradient1Function(func, b));
			p2.computeFisherInformation(a);

			final FisherInformationMatrix m1 = new FisherInformationMatrix(p1.getLinear(), n);
			final FisherInformationMatrix m2 = new FisherInformationMatrix(p2.getLinear(), n);
			final double[] crlb1 = m1.crlb();
			final double[] crlb2 = m2.crlb();
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

		final double[] N_ = new double[] { 20, 50, 100, 500 };
		final double[] b2_ = new double[] { 0, 1, 2, 4 };
		final double[] s_ = new double[] { 1, 1.2, 1.5 };
		final double[] x_ = new double[] { 4.8, 5, 5.5 };
		final double a = 100;
		final int size = 10;
		final Gaussian2DFunction f = GaussianFunctionFactory.create2D(1, size, size, GaussianFunctionFactory.FIT_ERF_CIRCLE,
				null);
		final PoissonGradientProcedure p = PoissonGradientProcedureFactory.create(f);
		final int ix = f.findGradientIndex(Gaussian2DFunction.X_POSITION);
		final int iy = f.findGradientIndex(Gaussian2DFunction.Y_POSITION);
		final double[] params = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
		for (final double N : N_)
		{
			params[Gaussian2DFunction.SIGNAL] = N;
			for (final double b2 : b2_)
			{
				params[Gaussian2DFunction.BACKGROUND] = b2;
				for (final double s : s_)
				{
					final double ss = s * a;
					params[Gaussian2DFunction.X_SD] = s;
					for (final double x : x_)
					{
						params[Gaussian2DFunction.X_POSITION] = x;
						for (final double y : x_)
						{
							params[Gaussian2DFunction.Y_POSITION] = y;
							p.computeFisherInformation(params);
							final FisherInformationMatrix m1 = new FisherInformationMatrix(p.getLinear(), p.n);
							final double[] crlb = m1.crlb();
							if (crlb == null)
								Assert.fail("No variance");
							@SuppressWarnings("null")
							final
							double o1 = Math.sqrt(crlb[ix]) * a;
							final double o2 = Math.sqrt(crlb[iy]) * a;
							final double e = Gaussian2DPeakResultHelper.getMLPrecisionX(a, ss, N, b2, false);
							//System.out.printf("e = %f  :  o  =   %f   %f\n", e, o1, o2);
							Assert.assertEquals(e, o1, e * 5e-2);
							Assert.assertEquals(e, o2, e * 5e-2);
						}
					}
				}
			}
		}
	}

	protected int[] createFakeData(RandomGenerator r, int nparams, int iter, ArrayList<double[]> paramsList,
			ArrayList<double[]> yList)
	{
		final int[] x = new int[blockWidth * blockWidth];
		for (int i = 0; i < x.length; i++)
			x[i] = i;
		for (int i = 0; i < iter; i++)
		{
			final double[] params = new double[nparams];
			final double[] y = createFakeData(r, params);
			paramsList.add(params);
			yList.add(y);
		}
		return x;
	}

	private double[] createFakeData(RandomGenerator r, double[] params)
	{
		final int n = blockWidth * blockWidth;

		for (int i = 0; i < params.length; i++)
			params[i] = r.nextDouble();

		final double[] y = new double[n];
		for (int i = 0; i < y.length; i++)
			y[i] = r.nextDouble();

		return y;
	}

	protected void createFakeParams(RandomGenerator r, int nparams, int iter, ArrayList<double[]> paramsList)
	{
		for (int i = 0; i < iter; i++)
		{
			final double[] params = new double[nparams];
			createFakeParams(r, params);
			paramsList.add(params);
		}
	}

	private static void createFakeParams(RandomGenerator r, double[] params)
	{
		for (int i = 0; i < params.length; i++)
			params[i] = r.nextDouble();
	}

	protected ArrayList<double[]> copyList(ArrayList<double[]> paramsList)
	{
		final ArrayList<double[]> params2List = new ArrayList<>(paramsList.size());
		for (int i = 0; i < paramsList.size(); i++)
			params2List.add(paramsList.get(i).clone());
		return params2List;
	}
}
