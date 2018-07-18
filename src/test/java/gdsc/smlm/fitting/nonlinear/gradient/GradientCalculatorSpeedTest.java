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

import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.util.Precision;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.utils.DoubleEquality;
import gdsc.core.utils.SimpleArrayUtils;
import gdsc.core.utils.Statistics;
import gdsc.smlm.fitting.linear.EJMLLinearSolver;
import gdsc.smlm.function.CameraNoiseModel;
import gdsc.smlm.function.NonLinearFunction;
import gdsc.smlm.function.PoissonCalculator;
import gdsc.smlm.function.gaussian.EllipticalGaussian2DFunction;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.function.gaussian.SingleCircularGaussian2DFunction;
import gdsc.smlm.function.gaussian.SingleEllipticalGaussian2DFunction;
import gdsc.smlm.function.gaussian.SingleFixedGaussian2DFunction;
import gdsc.smlm.function.gaussian.SingleFreeCircularGaussian2DFunction;
import gdsc.smlm.function.gaussian.SingleNBFixedGaussian2DFunction;
import gdsc.test.LogLevel;
import gdsc.test.TestLog;
import gdsc.test.TestSettings;

/**
 * Contains speed tests for the fastest method for calculating the Hessian and gradient vector
 * for use in NonLinearFit
 */
@SuppressWarnings({ "javadoc" })
public class GradientCalculatorSpeedTest
{
	DoubleEquality eq = new DoubleEquality(1e-6, 1e-16);

	int MAX_ITER = 20000;
	int blockWidth = 10;
	double Background = 0.5;
	double Amplitude = 100;
	double Angle = Math.PI;
	double Xpos = 5;
	double Ypos = 5;
	double Xwidth = 1.2;
	double Ywidth = 1.2;

	RandomDataGenerator rdg;

	@Test
	public void gradientCalculatorFactoryCreatesOptimisedCalculators()
	{
		Assert.assertEquals(GradientCalculatorFactory.newCalculator(3).getClass(), GradientCalculator3.class);
		Assert.assertEquals(GradientCalculatorFactory.newCalculator(4).getClass(), GradientCalculator4.class);
		Assert.assertEquals(GradientCalculatorFactory.newCalculator(5).getClass(), GradientCalculator5.class);
		Assert.assertEquals(GradientCalculatorFactory.newCalculator(6).getClass(), GradientCalculator6.class);
		Assert.assertEquals(GradientCalculatorFactory.newCalculator(7).getClass(), GradientCalculator7.class);
		Assert.assertEquals(GradientCalculatorFactory.newCalculator(13).getClass(), GradientCalculator.class);

		Assert.assertEquals(GradientCalculatorFactory.newCalculator(3, true).getClass(), MLEGradientCalculator3.class);
		Assert.assertEquals(GradientCalculatorFactory.newCalculator(4, true).getClass(), MLEGradientCalculator4.class);
		Assert.assertEquals(GradientCalculatorFactory.newCalculator(5, true).getClass(), MLEGradientCalculator5.class);
		Assert.assertEquals(GradientCalculatorFactory.newCalculator(6, true).getClass(), MLEGradientCalculator6.class);
		Assert.assertEquals(GradientCalculatorFactory.newCalculator(7, true).getClass(), MLEGradientCalculator7.class);
		Assert.assertEquals(GradientCalculatorFactory.newCalculator(13, true).getClass(), MLEGradientCalculator.class);
	}

	@Test
	public void mleGradientCalculator7ComputesSameAsGradientCalculator()
	{
		gradientCalculatorNComputesSameAsGradientCalculator(
				new SingleEllipticalGaussian2DFunction(blockWidth, blockWidth), 7, true);
	}

	@Test
	public void mleGradientCalculator7IsFasterThanGradientCalculator()
	{
		gradientCalculatorNIsFasterThanGradientCalculator(
				new SingleEllipticalGaussian2DFunction(blockWidth, blockWidth), 7, true);
	}

	@Test
	public void mleGradientCalculator6ComputesSameAsGradientCalculator()
	{
		gradientCalculatorNComputesSameAsGradientCalculator(
				new SingleFreeCircularGaussian2DFunction(blockWidth, blockWidth), 6, true);
	}

	@Test
	public void mleGradientCalculator6IsFasterThanGradientCalculator()
	{
		gradientCalculatorNIsFasterThanGradientCalculator(
				new SingleFreeCircularGaussian2DFunction(blockWidth, blockWidth), 6, true);
	}

	@Test
	public void mleGradientCalculator5ComputesSameAsGradientCalculator()
	{
		gradientCalculatorNComputesSameAsGradientCalculator(
				new SingleCircularGaussian2DFunction(blockWidth, blockWidth), 5, true);
	}

	@Test
	public void mleGradientCalculator5IsFasterThanGradientCalculator()
	{
		gradientCalculatorNIsFasterThanGradientCalculator(new SingleCircularGaussian2DFunction(blockWidth, blockWidth),
				5, true);
	}

	@Test
	public void mleGradientCalculator4ComputesSameAsGradientCalculator()
	{
		gradientCalculatorNComputesSameAsGradientCalculator(new SingleFixedGaussian2DFunction(blockWidth, blockWidth),
				4, true);
	}

	@Test
	public void mleGradientCalculator4IsFasterThanGradientCalculator()
	{
		gradientCalculatorNIsFasterThanGradientCalculator(new SingleFixedGaussian2DFunction(blockWidth, blockWidth), 4,
				true);
	}

	@Test
	public void mleGradientCalculator3ComputesSameAsGradientCalculator()
	{
		gradientCalculatorNComputesSameAsGradientCalculator(new SingleNBFixedGaussian2DFunction(blockWidth, blockWidth),
				3, true);
	}

	@Test
	public void mleGradientCalculator3IsFasterThanGradientCalculator()
	{
		gradientCalculatorNIsFasterThanGradientCalculator(new SingleNBFixedGaussian2DFunction(blockWidth, blockWidth),
				3, true);
	}

	@Test
	public void gradientCalculator7ComputesSameAsGradientCalculator()
	{
		gradientCalculatorNComputesSameAsGradientCalculator(
				new SingleEllipticalGaussian2DFunction(blockWidth, blockWidth), 7, false);
	}

	@Test
	public void gradientCalculator7IsFasterThanGradientCalculator()
	{
		gradientCalculatorNIsFasterThanGradientCalculator(
				new SingleEllipticalGaussian2DFunction(blockWidth, blockWidth), 7, false);
	}

	@Test
	public void gradientCalculator6ComputesSameAsGradientCalculator()
	{
		gradientCalculatorNComputesSameAsGradientCalculator(
				new SingleFreeCircularGaussian2DFunction(blockWidth, blockWidth), 6, false);
	}

	@Test
	public void gradientCalculator6IsFasterThanGradientCalculator()
	{
		gradientCalculatorNIsFasterThanGradientCalculator(
				new SingleFreeCircularGaussian2DFunction(blockWidth, blockWidth), 6, false);
	}

	@Test
	public void gradientCalculator5ComputesSameAsGradientCalculator()
	{
		gradientCalculatorNComputesSameAsGradientCalculator(
				new SingleCircularGaussian2DFunction(blockWidth, blockWidth), 5, false);
	}

	@Test
	public void gradientCalculator5IsFasterThanGradientCalculator()
	{
		gradientCalculatorNIsFasterThanGradientCalculator(new SingleCircularGaussian2DFunction(blockWidth, blockWidth),
				5, false);
	}

	@Test
	public void gradientCalculator4ComputesSameAsGradientCalculator()
	{
		gradientCalculatorNComputesSameAsGradientCalculator(new SingleFixedGaussian2DFunction(blockWidth, blockWidth),
				4, false);
	}

	@Test
	public void gradientCalculator4IsFasterThanGradientCalculator()
	{
		gradientCalculatorNIsFasterThanGradientCalculator(new SingleFixedGaussian2DFunction(blockWidth, blockWidth), 4,
				false);
	}

	@Test
	public void gradientCalculator3ComputesSameAsGradientCalculator()
	{
		gradientCalculatorNComputesSameAsGradientCalculator(new SingleNBFixedGaussian2DFunction(blockWidth, blockWidth),
				3, false);
	}

	@Test
	public void gradientCalculator3IsFasterThanGradientCalculator()
	{
		gradientCalculatorNIsFasterThanGradientCalculator(new SingleNBFixedGaussian2DFunction(blockWidth, blockWidth),
				3, false);
	}

	private void gradientCalculatorNComputesSameAsGradientCalculator(Gaussian2DFunction func, int nparams, boolean mle)
	{
		// Check the function is the correct size
		Assert.assertEquals(nparams, func.gradientIndices().length);

		final int iter = 50;
		rdg = new RandomDataGenerator(TestSettings.getRandomGenerator());

		final double[][] alpha = new double[nparams][nparams];
		final double[] beta = new double[nparams];
		final double[][] alpha2 = new double[nparams][nparams];
		final double[] beta2 = new double[nparams];

		final ArrayList<double[]> paramsList = new ArrayList<>(iter);
		final ArrayList<double[]> yList = new ArrayList<>(iter);

		final int[] x = createData(1, iter, paramsList, yList);

		final GradientCalculator calc = (mle) ? new MLEGradientCalculator(beta.length) : new GradientCalculator(beta.length);
		final GradientCalculator calc2 = GradientCalculatorFactory.newCalculator(nparams, mle);

		for (int i = 0; i < paramsList.size(); i++)
		{
			final double s = calc.findLinearised(x, yList.get(i), paramsList.get(i), alpha, beta, func);
			final double s2 = calc2.findLinearised(x, yList.get(i), paramsList.get(i), alpha2, beta2, func);
			Assert.assertTrue("Result: Not same @ " + i, eq.almostEqualRelativeOrAbsolute(s, s2));
			Assert.assertTrue("Observations: Not same beta @ " + i, eq.almostEqualRelativeOrAbsolute(beta, beta2));
			for (int j = 0; j < beta.length; j++)
				Assert.assertTrue("Observations: Not same alpha @ " + i,
						eq.almostEqualRelativeOrAbsolute(alpha[j], alpha2[j]));
		}

		for (int i = 0; i < paramsList.size(); i++)
		{
			final double s = calc.findLinearised(x.length, yList.get(i), paramsList.get(i), alpha, beta, func);
			final double s2 = calc2.findLinearised(x.length, yList.get(i), paramsList.get(i), alpha2, beta2, func);
			Assert.assertTrue("N-Result: Not same @ " + i, eq.almostEqualRelativeOrAbsolute(s, s2));
			Assert.assertTrue("N-observations: Not same beta @ " + i, eq.almostEqualRelativeOrAbsolute(beta, beta2));
			for (int j = 0; j < beta.length; j++)
				Assert.assertTrue("N-observations: Not same alpha @ " + i,
						eq.almostEqualRelativeOrAbsolute(alpha[j], alpha2[j]));
		}

		if (!mle)
		{
			func.setNoiseModel(CameraNoiseModel.createNoiseModel(10, 0, true));

			for (int i = 0; i < paramsList.size(); i++)
			{
				final double s = calc.findLinearised(x, yList.get(i), paramsList.get(i), alpha, beta, func);
				final double s2 = calc2.findLinearised(x, yList.get(i), paramsList.get(i), alpha2, beta2, func);
				Assert.assertTrue("Result+Noise: Not same @ " + i, eq.almostEqualRelativeOrAbsolute(s, s2));
				Assert.assertTrue("Observations+Noise: Not same beta @ " + i,
						eq.almostEqualRelativeOrAbsolute(beta, beta2));
				for (int j = 0; j < beta.length; j++)
					Assert.assertTrue("Observations+Noise: Not same alpha @ " + i,
							eq.almostEqualRelativeOrAbsolute(alpha[j], alpha2[j]));
			}

			for (int i = 0; i < paramsList.size(); i++)
			{
				final double s = calc.findLinearised(x.length, yList.get(i), paramsList.get(i), alpha, beta, func);
				final double s2 = calc2.findLinearised(x.length, yList.get(i), paramsList.get(i), alpha2, beta2, func);
				Assert.assertTrue("N-Result+Noise: Not same @ " + i, eq.almostEqualRelativeOrAbsolute(s, s2));
				Assert.assertTrue("N-Observations+Noise: Not same beta @ " + i,
						eq.almostEqualRelativeOrAbsolute(beta, beta2));
				for (int j = 0; j < beta.length; j++)
					Assert.assertTrue("N-Observations+Noise: Not same alpha @ " + i,
							eq.almostEqualRelativeOrAbsolute(alpha[j], alpha2[j]));
			}
		}
	}

	private void gradientCalculatorNIsFasterThanGradientCalculator(Gaussian2DFunction func, int nparams, boolean mle)
	{
		TestSettings.assumeSpeedTest();

		// Check the function is the correct size
		Assert.assertEquals(nparams, func.gradientIndices().length);

		final int iter = 10000;
		rdg = new RandomDataGenerator(TestSettings.getRandomGenerator());
		final double[][] alpha = new double[nparams][nparams];
		final double[] beta = new double[nparams];

		final ArrayList<double[]> paramsList = new ArrayList<>(iter);
		final ArrayList<double[]> yList = new ArrayList<>(iter);

		final int[] x = createData(1, iter, paramsList, yList);

		final GradientCalculator calc = (mle) ? new MLEGradientCalculator(beta.length) : new GradientCalculator(beta.length);
		final GradientCalculator calc2 = GradientCalculatorFactory.newCalculator(nparams, mle);

		for (int i = 0; i < paramsList.size(); i++)
			calc.findLinearised(x, yList.get(i), paramsList.get(i), alpha, beta, func);

		for (int i = 0; i < paramsList.size(); i++)
			calc2.findLinearised(x, yList.get(i), paramsList.get(i), alpha, beta, func);

		long start1 = System.nanoTime();
		for (int i = 0; i < paramsList.size(); i++)
			calc.findLinearised(x, yList.get(i), paramsList.get(i), alpha, beta, func);
		start1 = System.nanoTime() - start1;

		long start2 = System.nanoTime();
		for (int i = 0; i < paramsList.size(); i++)
			calc2.findLinearised(x, yList.get(i), paramsList.get(i), alpha, beta, func);
		start2 = System.nanoTime() - start2;

		TestLog.logSpeedTestResult(start2 < start1,
				"%sLinearised GradientCalculator = %d : GradientCalculator%d = %d : %fx\n", (mle) ? "MLE " : "", start1,
				nparams, start2, (1.0 * start1) / start2);
	}

	@Test
	public void gradientCalculatorAssumedXIsFasterThanGradientCalculator()
	{
		TestSettings.assumeSpeedTest();

		final int iter = 10000;
		rdg = new RandomDataGenerator(TestSettings.getRandomGenerator());

		final ArrayList<double[]> paramsList = new ArrayList<>(iter);
		final ArrayList<double[]> yList = new ArrayList<>(iter);

		final int[] x = createData(1, iter, paramsList, yList);

		final GradientCalculator calc = new GradientCalculator6();
		final GradientCalculator calc2 = new GradientCalculator6();
		final SingleFreeCircularGaussian2DFunction func = new SingleFreeCircularGaussian2DFunction(blockWidth, blockWidth);
		final int n = x.length;
		final int ng = func.getNumberOfGradients();
		final double[][] alpha = new double[ng][ng];
		final double[] beta = new double[ng];

		for (int i = 0; i < paramsList.size(); i++)
			calc.findLinearised(x, yList.get(i), paramsList.get(i), alpha, beta, func);

		for (int i = 0; i < paramsList.size(); i++)
			calc2.findLinearised(n, yList.get(i), paramsList.get(i), alpha, beta, func);

		long start1 = System.nanoTime();
		for (int i = 0; i < paramsList.size(); i++)
			calc.findLinearised(x, yList.get(i), paramsList.get(i), alpha, beta, func);
		start1 = System.nanoTime() - start1;

		long start2 = System.nanoTime();
		for (int i = 0; i < paramsList.size(); i++)
			calc2.findLinearised(n, yList.get(i), paramsList.get(i), alpha, beta, func);
		start2 = System.nanoTime() - start2;

		TestLog.logSpeedTestResult(start2 < start1,
				"GradientCalculator = %d : GradientCalculatorAssumed = %d : %fx\n", start1, start2,
				(1.0 * start1) / start2);
	}

	@Test
	public void gradientCalculatorComputesGradient()
	{
		gradientCalculatorComputesGradient(new GradientCalculator(7));
	}

	@Test
	public void mleGradientCalculatorComputesGradient()
	{
		gradientCalculatorComputesGradient(new MLEGradientCalculator(7));
	}

	private void gradientCalculatorComputesGradient(GradientCalculator calc)
	{
		final int nparams = calc.nparams;
		final Gaussian2DFunction func = new SingleEllipticalGaussian2DFunction(blockWidth, blockWidth);
		// Check the function is the correct size
		final int[] indices = func.gradientIndices();
		Assert.assertEquals(nparams, indices.length);

		final int iter = 50;
		rdg = new RandomDataGenerator(TestSettings.getRandomGenerator());

		final double[] beta = new double[nparams];
		final double[] beta2 = new double[nparams];

		final ArrayList<double[]> paramsList = new ArrayList<>(iter);
		final ArrayList<double[]> yList = new ArrayList<>(iter);

		final int[] x = createData(1, iter, paramsList, yList, true);

		final double delta = 1e-3;
		final DoubleEquality eq = new DoubleEquality(1e-3, 1e-3);

		for (int i = 0; i < paramsList.size(); i++)
		{
			final double[] y = yList.get(i);
			final double[] a = paramsList.get(i);
			final double[] a2 = a.clone();
			//double s =
			calc.evaluate(x, y, a, beta, func);

			for (int k = 0; k < nparams; k++)
			{
				final int j = indices[k];
				final double d = Precision.representableDelta(a[j], (a[j] == 0) ? 1e-3 : a[j] * delta);
				a2[j] = a[j] + d;
				final double s1 = calc.evaluate(x, y, a2, beta2, func);
				a2[j] = a[j] - d;
				final double s2 = calc.evaluate(x, y, a2, beta2, func);
				a2[j] = a[j];

				final double gradient = (s1 - s2) / (2 * d);
				//System.out.printf("[%d,%d] %f  (%s %f+/-%f)  %f  ?=  %f\n", i, j, s, func.getName(j), a[j], d, beta[k],
				//		gradient);
				Assert.assertTrue("Not same gradient @ " + j, eq.almostEqualRelativeOrAbsolute(beta[k], gradient));
			}
		}
	}

	@Test
	public void mleGradientCalculatorComputesLikelihood()
	{
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
			public double eval(int x, double[] dyda, double[] w) { return 0; }
			@Override
			public double evalw(int x, double[] w) { return 0; }
			@Override
			public boolean canComputeWeights() { return false; }
			@Override
			public int getNumberOfGradients() {	return 0; }
		};
		//@formatter:on

		final int[] xx = SimpleArrayUtils.newArray(100, 0, 1);
		final double[] xxx = SimpleArrayUtils.newArray(100, 0, 1.0);
		for (final double u : new double[] { 0.79, 2.5, 5.32 })
		{
			double ll = 0;
			final PoissonDistribution pd = new PoissonDistribution(u);
			for (final int x : xx)
			{
				double o = MLEGradientCalculator.likelihood(u, x);
				double e = pd.probability(x);
				Assert.assertEquals("likelihood", e, o, e * 1e-10);

				o = MLEGradientCalculator.logLikelihood(u, x);
				e = pd.logProbability(x);
				Assert.assertEquals("log likelihood", e, o, Math.abs(e) * 1e-10);

				ll += e;
			}

			final MLEGradientCalculator gc = new MLEGradientCalculator(1);
			final double o = gc.logLikelihood(xxx, new double[] { u }, func);

			Assert.assertEquals("sum log likelihood", ll, o, Math.abs(ll) * 1e-10);
		}
	}

	@Test
	public void gradientCalculatorComputesSameOutputWithBias()
	{

		final Gaussian2DFunction func = new SingleEllipticalGaussian2DFunction(blockWidth, blockWidth);
		final int nparams = func.getNumberOfGradients();
		final GradientCalculator calc = new GradientCalculator(nparams);
		final int n = func.size();

		final int iter = 50;
		rdg = new RandomDataGenerator(TestSettings.getRandomGenerator());

		final ArrayList<double[]> paramsList = new ArrayList<>(iter);
		final ArrayList<double[]> yList = new ArrayList<>(iter);

		final ArrayList<double[][]> alphaList = new ArrayList<>(iter);
		final ArrayList<double[]> betaList = new ArrayList<>(iter);
		final ArrayList<double[]> xList = new ArrayList<>(iter);

		// Manipulate the background
		final double defaultBackground = Background;
		final boolean report = TestSettings.allow(LogLevel.INFO);
		try
		{
			Background = 1e-2;
			createData(1, iter, paramsList, yList, true);

			final EJMLLinearSolver solver = new EJMLLinearSolver(1e-5, 1e-6);

			for (int i = 0; i < paramsList.size(); i++)
			{
				final double[] y = yList.get(i);
				final double[] a = paramsList.get(i);
				final double[][] alpha = new double[nparams][nparams];
				final double[] beta = new double[nparams];
				calc.findLinearised(n, y, a, alpha, beta, func);
				alphaList.add(alpha);
				betaList.add(beta.clone());
				for (int j = 0; j < nparams; j++)
					if (Math.abs(beta[j]) < 1e-6)
						System.out.printf("[%d] Tiny beta %s %g\n", i, func.getGradientParameterName(j), beta[j]);
				// Solve
				if (!solver.solve(alpha, beta))
					throw new AssertionError();
				xList.add(beta);
				//System.out.println(Arrays.toString(beta));
			}

			final double[][] alpha = new double[nparams][nparams];
			final double[] beta = new double[nparams];
			final Statistics[] rel = new Statistics[nparams];
			final Statistics[] abs = new Statistics[nparams];
			for (int i = 0; i < nparams; i++)
			{
				rel[i] = new Statistics();
				abs[i] = new Statistics();
			}

			//for (int b = 1; b < 1000; b *= 2)
			//for (double b : new double[] { -500, -100, -10, -1, -0.1, 0.1, 1, 10, 100, 500 })
			for (final double b : new double[] { -10, -1, -0.1, 0.1, 1, 10 })
			{
				if (report)
					for (int i = 0; i < nparams; i++)
					{
						rel[i].reset();
						abs[i].reset();
					}

				for (int i = 0; i < paramsList.size(); i++)
				{
					final double[] y = add(yList.get(i), b);
					final double[] a = paramsList.get(i).clone();
					a[0] += b;
					calc.findLinearised(n, y, a, alpha, beta, func);
					final double[][] alpha2 = alphaList.get(i);
					final double[] beta2 = betaList.get(i);
					final double[] x2 = xList.get(i);

					Assert.assertArrayEquals("Beta", beta2, beta, 1e-10);
					for (int j = 0; j < nparams; j++)
						Assert.assertArrayEquals("Alpha", alpha2[j], alpha[j], 1e-10);

					// Solve
					solver.solve(alpha, beta);
					Assert.assertArrayEquals("X", x2, beta, 1e-10);

					if (report)
						for (int j = 0; j < nparams; j++)
						{
							rel[j].add(DoubleEquality.relativeError(x2[j], beta[j]));
							abs[j].add(Math.abs(x2[j] - beta[j]));
						}
				}

				if (report)
					for (int i = 0; i < nparams; i++)
						System.out.printf("Bias = %.2f : %s : Rel %g +/- %g: Abs %g +/- %g\n", b,
								func.getGradientParameterName(i), rel[i].getMean(), rel[i].getStandardDeviation(),
								abs[i].getMean(), abs[i].getStandardDeviation());
			}
		}
		finally
		{
			Background = defaultBackground;
		}
	}

	private static double[] add(double[] d, double b)
	{
		d = d.clone();
		for (int i = 0; i < d.length; i++)
			d[i] += b;
		return d;
	}

	@Test
	public void mleCalculatorComputesLogLikelihoodRatio()
	{
		final EllipticalGaussian2DFunction func = new EllipticalGaussian2DFunction(1, blockWidth, blockWidth);
		final int n = blockWidth * blockWidth;
		final double[] a = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
		rdg = new RandomDataGenerator(TestSettings.getRandomGenerator());
		for (int run = 5; run-- > 0;)
		{
			a[Gaussian2DFunction.BACKGROUND] = random(Background);
			a[Gaussian2DFunction.SIGNAL] = random(Amplitude);
			a[Gaussian2DFunction.ANGLE] = random(Angle);
			a[Gaussian2DFunction.X_POSITION] = random(Xpos);
			a[Gaussian2DFunction.Y_POSITION] = random(Ypos);
			a[Gaussian2DFunction.X_SD] = random(Xwidth);
			a[Gaussian2DFunction.Y_SD] = random(Ywidth);

			// Simulate Poisson process
			func.initialise(a);
			final double[] x = SimpleArrayUtils.newArray(n, 0, 1.0);
			final double[] u = new double[x.length];
			for (int i = 0; i < n; i++)
			{
				u[i] = func.eval(i);
				if (u[i] > 0)
					x[i] = rdg.nextPoisson(u[i]);
			}

			final int ng = func.getNumberOfGradients();
			final double[][] alpha = new double[ng][ng];
			final double[] beta = new double[ng];

			final GradientCalculator calc = GradientCalculatorFactory.newCalculator(ng, true);

			final double llr = PoissonCalculator.logLikelihoodRatio(u, x);
			final double llr2 = calc.findLinearised(n, x, a, alpha, beta, func);
			//System.out.printf("llr=%f, llr2=%f\n", llr, llr2);
			Assert.assertEquals("Log-likelihood ratio", llr, llr2, llr * 1e-10);
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
		final EllipticalGaussian2DFunction func = new EllipticalGaussian2DFunction(npeaks, blockWidth, blockWidth);
		params[0] = random(Background);
		for (int i = 0, j = 0; i < npeaks; i++, j += Gaussian2DFunction.PARAMETERS_PER_PEAK)
		{
			params[j + Gaussian2DFunction.SIGNAL] = random(Amplitude);
			params[j + Gaussian2DFunction.ANGLE] = random(Angle);
			params[j + Gaussian2DFunction.X_POSITION] = random(Xpos);
			params[j + Gaussian2DFunction.Y_POSITION] = random(Ypos);
			params[j + Gaussian2DFunction.X_SD] = random(Xwidth);
			params[j + Gaussian2DFunction.Y_SD] = random(Ywidth);
		}

		final double[] dy_da = new double[params.length];
		final double[] y = new double[n];
		func.initialise(params);
		for (int i = 0; i < y.length; i++)
			// Add random Poisson noise
			y[i] = rdg.nextPoisson(func.eval(i, dy_da));

		if (randomiseParams)
		{
			params[0] = random(params[0]);
			for (int i = 0, j = 0; i < npeaks; i++, j += Gaussian2DFunction.PARAMETERS_PER_PEAK)
			{
				params[j + Gaussian2DFunction.SIGNAL] = random(params[j + Gaussian2DFunction.SIGNAL]);
				params[j + Gaussian2DFunction.ANGLE] = random(params[j + Gaussian2DFunction.ANGLE]);
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

	protected ArrayList<double[]> copyList(ArrayList<double[]> paramsList)
	{
		final ArrayList<double[]> params2List = new ArrayList<>(paramsList.size());
		for (int i = 0; i < paramsList.size(); i++)
			params2List.add(paramsList.get(i).clone());
		return params2List;
	}

	void log(String format, Object... args)
	{
		System.out.printf(format, args);
	}
}
