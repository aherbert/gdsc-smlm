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

import org.apache.commons.math3.util.Precision;
import org.apache.commons.rng.UniformRandomProvider;
import org.ejml.data.DenseMatrix64F;
import org.junit.jupiter.api.Assertions;

import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.RandomGeneratorAdapter;
import uk.ac.sussex.gdsc.core.utils.Statistics;
import uk.ac.sussex.gdsc.smlm.fitting.linear.EJMLLinearSolver;
import uk.ac.sussex.gdsc.smlm.function.DummyGradientFunction;
import uk.ac.sussex.gdsc.smlm.function.FakeGradientFunction;
import uk.ac.sussex.gdsc.smlm.function.Gradient1Function;
import uk.ac.sussex.gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import uk.ac.sussex.gdsc.smlm.function.gaussian.erf.ErfGaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.erf.SingleFreeCircularErfGaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.math3.distribution.CustomPoissonDistribution;
import uk.ac.sussex.gdsc.test.LogLevel;
import uk.ac.sussex.gdsc.test.TestCounter;
import uk.ac.sussex.gdsc.test.TestLog;
import uk.ac.sussex.gdsc.test.TestSettings;
import uk.ac.sussex.gdsc.test.junit5.ExtraAssertions;
import uk.ac.sussex.gdsc.test.junit5.ExtraAssumptions;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.junit5.SpeedTag;

/**
 * Contains speed tests for the fastest method for calculating the Hessian and gradient vector
 * for use in the LVM algorithm.
 * <p>
 * Note: This class is a test-bed for implementation strategies. The fastest strategy can then be used for other
 * gradient procedures.
 */
@SuppressWarnings({ "javadoc" })
public class LSQLVMGradientProcedureTest
{
	DoubleEquality eq = new DoubleEquality(1e-6, 1e-16);

	int MAX_ITER = 20000;
	int blockWidth = 10;
	double background = 0.5;
	double signal = 100;
	double angle = Math.PI;
	double xpos = 5;
	double ypos = 5;
	double xwidth = 1.2;
	double ywidth = 1.2;

	private static double random(UniformRandomProvider r, double d)
	{
		return d - d * 0.1 + r.nextDouble() * 0.2;
	}

	@SeededTest
	public void gradientProcedureFactoryCreatesOptimisedProcedures()
	{
		final double[] y = new double[0];
		Assertions.assertEquals(
				LSQLVMGradientProcedureMatrixFactory.create(y, null, new DummyGradientFunction(6)).getClass(),
				LSQLVMGradientProcedureMatrix6.class);
		Assertions.assertEquals(
				LSQLVMGradientProcedureMatrixFactory.create(y, null, new DummyGradientFunction(5)).getClass(),
				LSQLVMGradientProcedureMatrix5.class);
		Assertions.assertEquals(
				LSQLVMGradientProcedureMatrixFactory.create(y, null, new DummyGradientFunction(4)).getClass(),
				LSQLVMGradientProcedureMatrix4.class);

		Assertions.assertEquals(
				LSQLVMGradientProcedureLinearFactory.create(y, null, new DummyGradientFunction(6)).getClass(),
				LSQLVMGradientProcedureLinear6.class);
		Assertions.assertEquals(
				LSQLVMGradientProcedureLinearFactory.create(y, null, new DummyGradientFunction(5)).getClass(),
				LSQLVMGradientProcedureLinear5.class);
		Assertions.assertEquals(
				LSQLVMGradientProcedureLinearFactory.create(y, null, new DummyGradientFunction(4)).getClass(),
				LSQLVMGradientProcedureLinear4.class);
	}

	@SeededTest
	public void gradientProcedureLinearComputesSameAsGradientCalculator(RandomSeed seed)
	{
		gradientProcedureComputesSameAsGradientCalculator(seed, new LSQLVMGradientProcedureLinearFactory());
	}

	@SeededTest
	public void gradientProcedureMatrixComputesSameAsGradientCalculator(RandomSeed seed)
	{
		gradientProcedureComputesSameAsGradientCalculator(seed, new LSQLVMGradientProcedureMatrixFactory());
	}

	@SeededTest
	public void gradientProcedureComputesSameAsGradientCalculator(RandomSeed seed)
	{
		gradientProcedureComputesSameAsGradientCalculator(seed, new LSQLVMGradientProcedureFactory());
	}

	private void gradientProcedureComputesSameAsGradientCalculator(RandomSeed seed,
			BaseLSQLVMGradientProcedureFactory factory)
	{
		gradientProcedureComputesSameAsGradientCalculator(seed, 4, factory);
		gradientProcedureComputesSameAsGradientCalculator(seed, 5, factory);
		gradientProcedureComputesSameAsGradientCalculator(seed, 6, factory);
		gradientProcedureComputesSameAsGradientCalculator(seed, 11, factory);
		gradientProcedureComputesSameAsGradientCalculator(seed, 21, factory);
	}

	@SpeedTag
	@SeededTest
	public void gradientProcedureLinearIsNotSlowerThanGradientCalculator(RandomSeed seed)
	{
		gradientProcedureIsNotSlowerThanGradientCalculator(seed, new LSQLVMGradientProcedureLinearFactory());
	}

	@SpeedTag
	@SeededTest
	public void gradientProcedureMatrixIsNotSlowerThanGradientCalculator(RandomSeed seed)
	{
		gradientProcedureIsNotSlowerThanGradientCalculator(seed, new LSQLVMGradientProcedureMatrixFactory());
	}

	@SpeedTag
	@SeededTest
	public void gradientProcedureIsNotSlowerThanGradientCalculator(RandomSeed seed)
	{
		gradientProcedureIsNotSlowerThanGradientCalculator(seed, new LSQLVMGradientProcedureFactory());
	}

	private void gradientProcedureIsNotSlowerThanGradientCalculator(RandomSeed seed,
			BaseLSQLVMGradientProcedureFactory factory)
	{
		gradientProcedureIsNotSlowerThanGradientCalculator(seed, 4, factory);
		gradientProcedureIsNotSlowerThanGradientCalculator(seed, 5, factory);
		gradientProcedureIsNotSlowerThanGradientCalculator(seed, 6, factory);
		// 2 peaks
		gradientProcedureIsNotSlowerThanGradientCalculator(seed, 11, factory);
		// 4 peaks
		gradientProcedureIsNotSlowerThanGradientCalculator(seed, 21, factory);
	}

	private void gradientProcedureComputesSameAsGradientCalculator(RandomSeed seed, int nparams,
			BaseLSQLVMGradientProcedureFactory factory)
	{
		final int iter = 10;

		final double[][] alpha = new double[nparams][nparams];
		final double[] beta = new double[nparams];

		final ArrayList<double[]> paramsList = new ArrayList<>(iter);
		final ArrayList<double[]> yList = new ArrayList<>(iter);

		final int[] x = createFakeData(TestSettings.getRandomGenerator(seed.getSeed()), nparams, iter, paramsList,
				yList);
		final int n = x.length;
		final FakeGradientFunction func = new FakeGradientFunction(blockWidth, nparams);

		final GradientCalculator calc = GradientCalculatorFactory.newCalculator(nparams, false);

		final String name = factory.getClass().getSimpleName();
		for (int i = 0; i < paramsList.size(); i++)
		{
			final BaseLSQLVMGradientProcedure p = factory.createProcedure(yList.get(i), func);
			p.gradient(paramsList.get(i));
			final double s = p.value;
			final double s2 = calc.findLinearised(n, yList.get(i), paramsList.get(i), alpha, beta, func);
			// Exactly the same ...
			ExtraAssertions.assertEquals(s, s2, "%s Result: Not same @ %d", name, i);
			ExtraAssertions.assertArrayEquals(p.beta, beta, "%s Observations: Not same beta @ %d", name, i);

			final double[] al = p.getAlphaLinear();
			ExtraAssertions.assertArrayEquals(al, new DenseMatrix64F(alpha).data,
					"%s Observations: Not same alpha linear @ %d", name, i);

			final double[][] am = p.getAlphaMatrix();
			ExtraAssertions.assertArrayEquals(am, alpha, "%s Observations: Not same alpha matrix @ %d", name, i);
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

	private void gradientProcedureIsNotSlowerThanGradientCalculator(RandomSeed seed, final int nparams,
			final BaseLSQLVMGradientProcedureFactory factory)
	{
		ExtraAssumptions.assumeSpeedTest();

		final int iter = 1000;
		final double[][] alpha = new double[nparams][nparams];
		final double[] beta = new double[nparams];

		final ArrayList<double[]> paramsList = new ArrayList<>(iter);
		final ArrayList<double[]> yList = new ArrayList<>(iter);

		final int[] x = createFakeData(TestSettings.getRandomGenerator(seed.getSeed()), nparams, iter, paramsList,
				yList);
		final int n = x.length;
		final FakeGradientFunction func = new FakeGradientFunction(blockWidth, nparams);

		final GradientCalculator calc = GradientCalculatorFactory.newCalculator(nparams, false);

		for (int i = 0; i < paramsList.size(); i++)
			calc.findLinearised(n, yList.get(i), paramsList.get(i), alpha, beta, func);

		for (int i = 0; i < paramsList.size(); i++)
		{
			final BaseLSQLVMGradientProcedure p = factory.createProcedure(yList.get(i), func);
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
					final GradientCalculator calc = GradientCalculatorFactory.newCalculator(nparams, false);
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
					final BaseLSQLVMGradientProcedure p = factory.createProcedure(yList.get(i), func);
					for (int j = loops; j-- > 0;)
						p.gradient(paramsList.get(k++ % iter));
				}
			}
		};
		final long time2 = t2.getTime();

		TestLog.logSpeedTestResult(time2 < time1, "GradientCalculator = %d : %s %d = %d : %fx\n", time1,
				factory.getClass().getSimpleName(), nparams, time2, (1.0 * time1) / time2);
	}

	@SeededTest
	public void gradientProcedureUnrolledComputesSameAsGradientProcedure(RandomSeed seed)
	{
		// Test the method that will be used for the standard and unrolled versions
		// for all other 'gradient procedures'
		gradientProcedureUnrolledComputesSameAsGradientProcedure(seed, 4);
		gradientProcedureUnrolledComputesSameAsGradientProcedure(seed, 5);
		gradientProcedureUnrolledComputesSameAsGradientProcedure(seed, 6);
	}

	private void gradientProcedureUnrolledComputesSameAsGradientProcedure(RandomSeed seed, int nparams)
	{
		final int iter = 10;

		final ArrayList<double[]> paramsList = new ArrayList<>(iter);
		final ArrayList<double[]> yList = new ArrayList<>(iter);

		createFakeData(TestSettings.getRandomGenerator(seed.getSeed()), nparams, iter, paramsList, yList);
		final FakeGradientFunction func = new FakeGradientFunction(blockWidth, nparams);

		final String name = GradientCalculator.class.getSimpleName();
		for (int i = 0; i < paramsList.size(); i++)
		{
			final BaseLSQLVMGradientProcedure p1 = LSQLVMGradientProcedureFactory.create(yList.get(i), func);
			p1.gradient(paramsList.get(i));

			final BaseLSQLVMGradientProcedure p2 = new LSQLVMGradientProcedure(yList.get(i), func);
			p2.gradient(paramsList.get(i));

			// Exactly the same ...
			ExtraAssertions.assertEquals(p1.value, p2.value, "%s Result: Not same @ %d", name, i);

			ExtraAssertions.assertArrayEquals(p1.beta, p2.beta, "%s Observations: Not same beta @ %d", name, i);

			ExtraAssertions.assertArrayEquals(p1.getAlphaLinear(), p2.getAlphaLinear(),
					"%s Observations: Not same alpha linear @ %d", name, i);

			final double[][] am1 = p1.getAlphaMatrix();
			final double[][] am2 = p2.getAlphaMatrix();
			ExtraAssertions.assertArrayEquals(am1, am2, "%s Observations: Not same alpha matrix @ %d", name, i);
		}
	}

	@SpeedTag
	@SeededTest
	public void gradientProcedureIsFasterUnrolledThanGradientProcedureMatrix(RandomSeed seed)
	{
		gradientProcedure2IsFasterUnrolledThanGradientProcedure1(seed, new LSQLVMGradientProcedureMatrixFactory(),
				new LSQLVMGradientProcedureFactory());
	}

	@SpeedTag
	@SeededTest
	public void gradientProcedureLinearIsFasterUnrolledThanGradientProcedureMatrix(RandomSeed seed)
	{
		gradientProcedure2IsFasterUnrolledThanGradientProcedure1(seed, new LSQLVMGradientProcedureMatrixFactory(),
				new LSQLVMGradientProcedureLinearFactory());
	}

	private void gradientProcedure2IsFasterUnrolledThanGradientProcedure1(RandomSeed seed,
			BaseLSQLVMGradientProcedureFactory factory1, BaseLSQLVMGradientProcedureFactory factory2)
	{
		// Assert the unrolled versions
		gradientProcedureLinearIsFasterThanGradientProcedureMatrix(seed, 4, factory1, factory2, true);
		gradientProcedureLinearIsFasterThanGradientProcedureMatrix(seed, 5, factory1, factory2, true);
		gradientProcedureLinearIsFasterThanGradientProcedureMatrix(seed, 6, factory1, factory2, true);
		gradientProcedureLinearIsFasterThanGradientProcedureMatrix(seed, 11, factory1, factory2, false);
		gradientProcedureLinearIsFasterThanGradientProcedureMatrix(seed, 21, factory1, factory2, false);
	}

	private void gradientProcedureLinearIsFasterThanGradientProcedureMatrix(RandomSeed seed, final int nparams,
			final BaseLSQLVMGradientProcedureFactory factory1, final BaseLSQLVMGradientProcedureFactory factory2,
			boolean doAssert)
	{
		ExtraAssumptions.assumeSpeedTest();

		final int iter = 100;

		final ArrayList<double[]> paramsList = new ArrayList<>(iter);
		final ArrayList<double[]> yList = new ArrayList<>(iter);

		createData(TestSettings.getRandomGenerator(seed.getSeed()), 1, iter, paramsList, yList);

		// Remove the timing of the function call by creating a dummy function
		final Gradient1Function func = new FakeGradientFunction(blockWidth, nparams);

		for (int i = 0; i < paramsList.size(); i++)
		{
			final BaseLSQLVMGradientProcedure p1 = factory1.createProcedure(yList.get(i), func);
			p1.gradient(paramsList.get(i));
			p1.gradient(paramsList.get(i));

			final BaseLSQLVMGradientProcedure p2 = factory2.createProcedure(yList.get(i), func);
			p2.gradient(paramsList.get(i));
			p2.gradient(paramsList.get(i));

			// Check they are the same
			ExtraAssertions.assertArrayEquals(p1.getAlphaLinear(), p2.getAlphaLinear(), "A %d", i);
			ExtraAssertions.assertArrayEquals(p1.beta, p2.beta, "B %d", i);
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
					final BaseLSQLVMGradientProcedure p = factory1.createProcedure(yList.get(i), func);
					for (int j = loops; j-- > 0;)
						p.gradient(paramsList.get(k++ % iter));
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
					final BaseLSQLVMGradientProcedure p2 = factory2.createProcedure(yList.get(i), func);
					for (int j = loops; j-- > 0;)
						p2.gradient(paramsList.get(k++ % iter));
				}
			}
		};
		final long time2 = t2.getTime();

		TestLog.logSpeedTestResult(!doAssert || time2 < time1, "%s = %d : %s %d = %d : %fx\n",
				factory1.getClass().getSimpleName(), time1, factory2.getClass().getSimpleName(), nparams, time2,
				(1.0 * time1) / time2);
	}

	@SeededTest
	public void gradientProcedureComputesGradient(RandomSeed seed)
	{
		gradientProcedureComputesGradient(seed, new SingleFreeCircularErfGaussian2DFunction(blockWidth, blockWidth));
	}

	private void gradientProcedureComputesGradient(RandomSeed seed, ErfGaussian2DFunction func)
	{
		final int nparams = func.getNumberOfGradients();
		final int[] indices = func.gradientIndices();

		final int iter = 100;

		final ArrayList<double[]> paramsList = new ArrayList<>(iter);
		final ArrayList<double[]> yList = new ArrayList<>(iter);

		createData(TestSettings.getRandomGenerator(seed.getSeed()), 1, iter, paramsList, yList, true);

		// for the gradients
		final double delta = 1e-4;
		final DoubleEquality eq = new DoubleEquality(5e-2, 1e-16);

		// Must compute most of the time
		final int failureLimit = TestCounter.computeFailureLimit(iter, 0.1);
		final TestCounter failCounter = new TestCounter(failureLimit, nparams);

		for (int i = 0; i < paramsList.size(); i++)
		{
			final int ii = i;
			final double[] y = yList.get(i);
			final double[] a = paramsList.get(i);
			final double[] a2 = a.clone();
			final BaseLSQLVMGradientProcedure p = LSQLVMGradientProcedureFactory.create(y, func);
			p.gradient(a);
			//double s = p.ssx;
			final double[] beta = p.beta.clone();
			for (int j = 0; j < nparams; j++)
			{
				final int jj = j;
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
				//System.out.printf("[%d,%d] %f  (%s %f+/-%f)  %f  ?=  %f\n", i, k, s, func.getName(k), a[k], d, beta[j],
				//		gradient);
				failCounter.run(j, () -> {
					return eq.almostEqualRelativeOrAbsolute(beta[jj], gradient);
				}, () -> {
					ExtraAssertions.fail("Not same gradient @ %d,%d: %s != %s (error=%s)", ii, jj, beta[jj], gradient,
							DoubleEquality.relativeError(beta[jj], gradient));
				});
			}
		}
	}

	@SeededTest
	public void gradientProcedureComputesSameOutputWithBias(RandomSeed seed)
	{
		final ErfGaussian2DFunction func = new SingleFreeCircularErfGaussian2DFunction(blockWidth, blockWidth);
		final int nparams = func.getNumberOfGradients();

		final int iter = 100;

		final LogLevel logLevel = LogLevel.INFO;
		final boolean debug = TestSettings.allow(logLevel);

		final ArrayList<double[]> paramsList = new ArrayList<>(iter);
		final ArrayList<double[]> yList = new ArrayList<>(iter);

		final ArrayList<double[]> alphaList = new ArrayList<>(iter);
		final ArrayList<double[]> betaList = new ArrayList<>(iter);
		final ArrayList<double[]> xList = new ArrayList<>(iter);

		// Manipulate the background
		final double defaultBackground = background;
		try
		{
			background = 1e-2;
			createData(TestSettings.getRandomGenerator(seed.getSeed()), 1, iter, paramsList, yList, true);

			final EJMLLinearSolver solver = new EJMLLinearSolver(1e-5, 1e-6);

			for (int i = 0; i < paramsList.size(); i++)
			{
				final double[] y = yList.get(i);
				final double[] a = paramsList.get(i);
				final BaseLSQLVMGradientProcedure p = LSQLVMGradientProcedureFactory.create(y, func);
				p.gradient(a);
				final double[] beta = p.beta;
				alphaList.add(p.getAlphaLinear());
				betaList.add(beta.clone());
				for (int j = 0; j < nparams; j++)
					if (Math.abs(beta[j]) < 1e-6)
						TestLog.info("[%d] Tiny beta %s %g\n", i, func.getGradientParameterName(j), beta[j]);
				// Solve
				if (!solver.solve(p.getAlphaMatrix(), beta))
					throw new AssertionError();
				xList.add(beta);
				//System.out.println(Arrays.toString(beta));
			}

			//for (int b = 1; b < 1000; b *= 2)
			for (final double b : new double[] { -500, -100, -10, -1, -0.1, 0, 0.1, 1, 10, 100, 500 })
			{
				final Statistics[] rel = new Statistics[nparams];
				final Statistics[] abs = new Statistics[nparams];
				if (debug)
				{
					for (int i = 0; i < nparams; i++)
					{
						rel[i] = new Statistics();
						abs[i] = new Statistics();
					}
				}

				for (int i = 0; i < paramsList.size(); i++)
				{
					final double[] y = add(yList.get(i), b);
					final double[] a = paramsList.get(i).clone();
					a[0] += b;
					final BaseLSQLVMGradientProcedure p = LSQLVMGradientProcedureFactory.create(y, func);
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

					if (debug)
					{
						for (int j = 0; j < nparams; j++)
						{
							rel[j].add(DoubleEquality.relativeError(x2[j], beta[j]));
							abs[j].add(Math.abs(x2[j] - beta[j]));
						}
					}
				}

				if (debug)
					for (int i = 0; i < nparams; i++)
						TestLog.log(logLevel, "Bias = %.2f : %s : Rel %g +/- %g: Abs %g +/- %g\n", b,
								func.getGradientParameterName(i), rel[i].getMean(), rel[i].getStandardDeviation(),
								abs[i].getMean(), abs[i].getStandardDeviation());
			}
		}
		finally
		{
			background = defaultBackground;
		}
	}

	private static double[] add(double[] d, double b)
	{
		d = d.clone();
		for (int i = 0; i < d.length; i++)
			d[i] += b;
		return d;
	}

	/**
	 * Create random elliptical Gaussian data an returns the data plus an estimate of the parameters.
	 * Only the chosen parameters are randomised and returned for a maximum of (background, amplitude, angle, xpos,
	 * ypos, xwidth, ywidth }
	 *
	 * @param r
	 *            the random
	 * @param npeaks
	 *            the npeaks
	 * @param params
	 *            set on output
	 * @param randomiseParams
	 *            Set to true to randomise the params
	 * @return the double[]
	 */
	private double[] doubleCreateGaussianData(UniformRandomProvider r, int npeaks, double[] params,
			boolean randomiseParams)
	{
		final int n = blockWidth * blockWidth;

		// Generate a 2D Gaussian
		final ErfGaussian2DFunction func = (ErfGaussian2DFunction) GaussianFunctionFactory.create2D(npeaks, blockWidth,
				blockWidth, GaussianFunctionFactory.FIT_ERF_FREE_CIRCLE, null);
		params[0] = random(r, background);
		for (int i = 0, j = 1; i < npeaks; i++, j += 6)
		{
			params[j] = random(r, signal);
			params[j + 2] = random(r, xpos);
			params[j + 3] = random(r, ypos);
			params[j + 4] = random(r, xwidth);
			params[j + 5] = random(r, ywidth);
		}

		final double[] y = new double[n];
		func.initialise(params);
		CustomPoissonDistribution pd = new CustomPoissonDistribution(new RandomGeneratorAdapter(r), 1);
		for (int i = 0; i < y.length; i++)
		{
			// Add random Poisson noise
			final double u = func.eval(i);
			pd.setMean(u);
			y[i] = pd.sample();
		}

		if (randomiseParams)
		{
			params[0] = random(r, params[0]);
			for (int i = 0, j = 1; i < npeaks; i++, j += 6)
			{
				params[j] = random(r, params[j]);
				params[j + 2] = random(r, params[j + 2]);
				params[j + 3] = random(r, params[j + 3]);
				params[j + 4] = random(r, params[j + 4]);
				params[j + 5] = random(r, params[j + 5]);
			}
		}

		return y;
	}

	protected int[] createData(final UniformRandomProvider r, int npeaks, int iter, ArrayList<double[]> paramsList,
			ArrayList<double[]> yList)
	{
		return createData(r, npeaks, iter, paramsList, yList, true);
	}

	protected int[] createData(final UniformRandomProvider r, int npeaks, int iter, ArrayList<double[]> paramsList,
			ArrayList<double[]> yList, boolean randomiseParams)
	{
		final int[] x = new int[blockWidth * blockWidth];
		for (int i = 0; i < x.length; i++)
			x[i] = i;
		for (int i = 0; i < iter; i++)
		{
			final double[] params = new double[1 + 6 * npeaks];
			final double[] y = doubleCreateGaussianData(r, npeaks, params, randomiseParams);
			paramsList.add(params);
			yList.add(y);
		}
		return x;
	}

	protected int[] createFakeData(final UniformRandomProvider r, int nparams, int iter, ArrayList<double[]> paramsList,
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

	private double[] createFakeData(final UniformRandomProvider r, double[] params)
	{
		final int n = blockWidth * blockWidth;

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
