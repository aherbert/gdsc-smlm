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
package gdsc.smlm.function;

import java.util.Arrays;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.Precision;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.utils.DoubleEquality;
import gdsc.core.utils.SimpleArrayUtils;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import gdsc.test.TestAssert;
import gdsc.test.TestSettings;

public class PoissonLikelihoodWrapperTest
{
	double alpha = 1 / 40.0;
	double[] photons = { 0.25, 0.5, 1, 2, 4, 10, 100, 1000 };
	// Set this at the range output from cumulativeProbabilityIsOneWithIntegerData
	int[] maxRange = { 6, 7, 10, 13, 17, 29, 149, 1141 };

	DoubleEquality eqPerDatum = new DoubleEquality(5e-2, 0.01);
	DoubleEquality eq = new DoubleEquality(5e-3, 0.001);

	static String[] NAME;
	static
	{
		NAME = new String[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
		for (int i = 0; i < NAME.length; i++)
			NAME[i] = Gaussian2DFunction.getName(i);
	}

	// Compute as per Numerical Recipes 5.7.
	// Approximate error accuracy in single precision: Ef
	// Step size for derivatives:
	// h ~ (Ef)^(1/3) * xc
	// xc is the characteristic scale over which x changes, assumed to be 1 (not x as per NR since x is close to zero)
	final double h_ = 0.01; //(double) (Math.pow(1e-3f, 1.0 / 3));

	int[] testx = new int[] { 4, 5, 6 };
	int[] testy = new int[] { 4, 5, 6 };
	// Do not test zero background since this is an edge case for the likelihood function
	double[] testbackground_ = new double[] { 10, 400 };

	double[] testsignal1_ = new double[] { 15, 55, 105 };
	double[] testcx1_ = new double[] { 4.9, 5.3 };
	double[] testcy1_ = new double[] { 4.8, 5.2 };
	double[] testcz1_ = new double[] { -1.5, 1.0 };
	double[][] testw1_ = new double[][] { { 1.1, 1.4 }, { 1.1, 1.7 }, { 1.5, 1.2 }, { 1.3, 1.7 }, };
	double[] testangle1_ = new double[] { Math.PI / 5, Math.PI / 3 };

	double[] testbackground, testsignal1, testangle1, testcx1, testcy1, testcz1;
	double[][] testw1;

	int maxx = 10;
	double background = 50;
	double angle = 0;
	double width = 5;

	@Test
	public void fitFixedComputesGradientPerDatum()
	{
		TestSettings.assumeMediumComplexity();
		functionComputesGradientPerDatum(GaussianFunctionFactory.FIT_FIXED);
	}

	@Test
	public void fitCircleComputesGradientPerDatum()
	{
		functionComputesGradientPerDatum(GaussianFunctionFactory.FIT_CIRCLE);
	}

	@Test
	public void fitFreeCircleComputesGradientPerDatum()
	{
		functionComputesGradientPerDatum(GaussianFunctionFactory.FIT_FREE_CIRCLE);
	}

	@Test
	public void fitEllipticalComputesGradientPerDatum()
	{
		TestSettings.assumeMediumComplexity();
		functionComputesGradientPerDatum(GaussianFunctionFactory.FIT_ELLIPTICAL);
	}

	@Test
	public void fitNBFixedComputesGradientPerDatum()
	{
		TestSettings.assumeMediumComplexity();
		functionComputesGradientPerDatum(GaussianFunctionFactory.FIT_SIMPLE_NB_FIXED);
	}

	@Test
	public void fitNBCircleComputesGradientPerDatum()
	{
		TestSettings.assumeMediumComplexity();
		functionComputesGradientPerDatum(GaussianFunctionFactory.FIT_SIMPLE_NB_CIRCLE);
	}

	@Test
	public void fitNBFreeCircleComputesGradientPerDatum()
	{
		TestSettings.assumeMediumComplexity();
		functionComputesGradientPerDatum(GaussianFunctionFactory.FIT_SIMPLE_NB_FREE_CIRCLE);
	}

	@Test
	public void fitNBEllipticalComputesGradientPerDatum()
	{
		TestSettings.assumeMediumComplexity();
		functionComputesGradientPerDatum(GaussianFunctionFactory.FIT_SIMPLE_NB_ELLIPTICAL);
	}

	private void functionComputesGradientPerDatum(int flags)
	{
		Gaussian2DFunction f1 = GaussianFunctionFactory.create2D(1, maxx, maxx, flags, null);
		// Setup
		testbackground = testbackground_;
		testsignal1 = testsignal1_;
		testcx1 = testcx1_;
		testcy1 = testcy1_;
		testcz1 = testcz1_;
		testw1 = testw1_;
		testangle1 = testangle1_;
		if (!f1.evaluatesBackground())
		{
			testbackground = new double[] { testbackground[0] };
		}
		if (!f1.evaluatesSignal())
		{
			testsignal1 = new double[] { testsignal1[0] };
		}
		// Position is always evaluated

		if (!f1.evaluatesZ())
		{
			testcz1 = new double[] { 0 };
		}
		boolean noSecondWidth = false;
		if (!f1.evaluatesSD0())
		{
			// Just use 1 width
			testw1 = new double[][] { testw1[0] };
			// If no width 0 then assume we have no width 1 as well
			noSecondWidth = true;
		}
		else if (!f1.evaluatesSD1())
		{
			// No evaluation of second width needs only variation in width 0 so truncate 
			testw1 = Arrays.copyOf(testw1, 2);
			noSecondWidth = true;
		}
		if (noSecondWidth)
		{
			for (int i = 0; i < testw1.length; i++)
			{
				testw1[i][1] = testw1[i][0];
			}
		}
		if (!f1.evaluatesAngle())
		{
			testangle1 = new double[] { 0 };
		}

		if (f1.evaluatesBackground())
			functionComputesTargetGradientPerDatum(f1, Gaussian2DFunction.BACKGROUND);
		if (f1.evaluatesSignal())
			functionComputesTargetGradientPerDatum(f1, Gaussian2DFunction.SIGNAL);
		functionComputesTargetGradientPerDatum(f1, Gaussian2DFunction.X_POSITION);
		functionComputesTargetGradientPerDatum(f1, Gaussian2DFunction.Y_POSITION);
		if (f1.evaluatesZ())
			functionComputesTargetGradientPerDatum(f1, Gaussian2DFunction.Z_POSITION);
		if (f1.evaluatesSD0())
			functionComputesTargetGradientPerDatum(f1, Gaussian2DFunction.X_SD);
		if (f1.evaluatesSD1())
			functionComputesTargetGradientPerDatum(f1, Gaussian2DFunction.Y_SD);
		if (f1.evaluatesAngle())
			functionComputesTargetGradientPerDatum(f1, Gaussian2DFunction.ANGLE);
	}

	private void functionComputesTargetGradientPerDatum(Gaussian2DFunction f1, int targetParameter)
	{
		int[] indices = f1.gradientIndices();
		int gradientIndex = findGradientIndex(f1, targetParameter);
		double[] dyda = new double[indices.length];
		double[] a;

		PoissonLikelihoodWrapper ff1;

		int n = maxx * maxx;
		int count = 0, total = 0;

		for (double background : testbackground)
			for (double signal1 : testsignal1)
				for (double cx1 : testcx1)
					for (double cy1 : testcy1)
						for (double cz1 : testcz1)
							for (double[] w1 : testw1)
								for (double angle1 : testangle1)
								{
									a = createParameters(background, signal1, cx1, cy1, cz1, w1[0], w1[1], angle1);

									// Create y as a function we would want to move towards
									double[] a2 = a.clone();
									a2[targetParameter] *= 1.1;
									f1.initialise(a2);
									double[] data = new double[maxx * maxx];
									for (int i = 0; i < n; i++)
										data[i] = f1.eval(i);

									ff1 = new PoissonLikelihoodWrapper(f1, a, data, n, alpha);

									// Numerically solve gradient. 
									// Calculate the step size h to be an exact numerical representation
									final double xx = a[targetParameter];

									// Get h to minimise roundoff error
									double h = Precision.representableDelta(xx, h_);

									for (int x : testx)
										for (int y : testy)
										{
											int i = y * maxx + x;
											a[targetParameter] = xx;
											ff1.likelihood(getVariables(indices, a), dyda, i);

											// Evaluate at (x+h) and (x-h)
											a[targetParameter] = xx + h;
											double value2 = ff1.likelihood(getVariables(indices, a), i);

											a[targetParameter] = xx - h;
											double value3 = ff1.likelihood(getVariables(indices, a), i);

											double gradient = (value2 - value3) / (2 * h);
											boolean ok = Math.signum(gradient) == Math.signum(dyda[gradientIndex]) ||
													Math.abs(gradient - dyda[gradientIndex]) < 0.1;
											//logf("[%s-%s]/2*%g : %g == %g\n", "" + value2, "" + value3, h, gradient,
											//		dyda[gradientIndex]);
											if (!ok)
												Assert.assertTrue(NAME[targetParameter] + ": " + gradient + " != " +
														dyda[gradientIndex], ok);
											ok = eqPerDatum.almostEqualRelativeOrAbsolute(gradient,
													dyda[gradientIndex]);
											if (ok)
												count++;
											total++;
										}
								}
		double p = (100.0 * count) / total;
		TestSettings.info("Per Datum %s : %s = %d / %d (%.2f)\n", f1.getClass().getSimpleName(), NAME[targetParameter], count, total,
				p);
		Assert.assertTrue(NAME[targetParameter] + " fraction too low per datum: " + p, p > 90);
	}

	@Test
	public void fitFixedComputesGradient()
	{
		TestSettings.assumeMediumComplexity();
		functionComputesGradient(GaussianFunctionFactory.FIT_FIXED);
	}

	@Test
	public void fitCircleComputesGradient()
	{
		functionComputesGradient(GaussianFunctionFactory.FIT_CIRCLE);
	}

	@Test
	public void fitFreeCircleComputesGradient()
	{
		functionComputesGradient(GaussianFunctionFactory.FIT_FREE_CIRCLE);
	}

	@Test
	public void fitEllipticalComputesGradient()
	{
		TestSettings.assumeMediumComplexity();
		// The elliptical function gradient evaluation is worse 
		DoubleEquality tmp = eq;
		eq = eqPerDatum;
		functionComputesGradient(GaussianFunctionFactory.FIT_ELLIPTICAL);
		eq = tmp;
	}

	@Test
	public void fitNBFixedComputesGradient()
	{
		TestSettings.assumeMediumComplexity();
		functionComputesGradient(GaussianFunctionFactory.FIT_SIMPLE_NB_FIXED);
	}

	@Test
	public void fitNBCircleComputesGradient()
	{
		TestSettings.assumeMediumComplexity();
		functionComputesGradient(GaussianFunctionFactory.FIT_SIMPLE_NB_CIRCLE);
	}

	@Test
	public void fitNBFreeCircleComputesGradient()
	{
		TestSettings.assumeMediumComplexity();
		functionComputesGradient(GaussianFunctionFactory.FIT_SIMPLE_NB_FREE_CIRCLE);
	}

	@Test
	public void fitNBEllipticalComputesGradient()
	{
		TestSettings.assumeMediumComplexity();
		// The elliptical function gradient evaluation is worse 
		DoubleEquality tmp = eq;
		eq = eqPerDatum;
		functionComputesGradient(GaussianFunctionFactory.FIT_SIMPLE_NB_ELLIPTICAL);
		eq = tmp;
	}

	private void functionComputesGradient(int flags)
	{
		Gaussian2DFunction f1 = GaussianFunctionFactory.create2D(1, maxx, maxx, flags, null);
		testbackground = testbackground_;
		testsignal1 = testsignal1_;
		testcx1 = testcx1_;
		testcy1 = testcy1_;
		testcz1 = testcz1_;
		testw1 = testw1_;
		testangle1 = testangle1_;
		if (!f1.evaluatesBackground())
		{
			testbackground = new double[] { testbackground[0] };
		}
		if (!f1.evaluatesSignal())
		{
			testsignal1 = new double[] { testsignal1[0] };
		}
		// Position is always evaluated

		if (!f1.evaluatesZ())
		{
			testcz1 = new double[] { 0 };
		}
		boolean noSecondWidth = false;
		if (!f1.evaluatesSD0())
		{
			// Just use 1 width
			testw1 = new double[][] { testw1[0] };
			// If no width 0 then assume we have no width 1 as well
			noSecondWidth = true;
		}
		else if (!f1.evaluatesSD1())
		{
			// No evaluation of second width needs only variation in width 0 so truncate 
			testw1 = Arrays.copyOf(testw1, 2);
			noSecondWidth = true;
		}
		if (noSecondWidth)
		{
			for (int i = 0; i < testw1.length; i++)
			{
				testw1[i][1] = testw1[i][0];
			}
		}
		if (!f1.evaluatesAngle())
		{
			testangle1 = new double[] { 0 };
		}

		double fraction = 90;
		if (f1.evaluatesBackground())
			functionComputesTargetGradient(f1, Gaussian2DFunction.BACKGROUND, fraction);
		if (f1.evaluatesSignal())
			functionComputesTargetGradient(f1, Gaussian2DFunction.SIGNAL, fraction);
		functionComputesTargetGradient(f1, Gaussian2DFunction.X_POSITION, fraction);
		functionComputesTargetGradient(f1, Gaussian2DFunction.Y_POSITION, fraction);
		if (f1.evaluatesZ())
			functionComputesTargetGradient(f1, Gaussian2DFunction.Z_POSITION, fraction);
		if (f1.evaluatesSD0())
			functionComputesTargetGradient(f1, Gaussian2DFunction.X_SD, fraction);
		if (f1.evaluatesSD1())
			functionComputesTargetGradient(f1, Gaussian2DFunction.Y_SD, fraction);
		if (f1.evaluatesAngle())
			functionComputesTargetGradient(f1, Gaussian2DFunction.ANGLE, fraction);
	}

	private void functionComputesTargetGradient(Gaussian2DFunction f1, int targetParameter, double threshold)
	{
		int[] indices = f1.gradientIndices();
		int gradientIndex = findGradientIndex(f1, targetParameter);
		double[] dyda = new double[indices.length];
		double[] a;

		PoissonLikelihoodWrapper ff1;

		int n = maxx * maxx;
		int count = 0, total = 0;

		for (double background : testbackground)
			for (double signal1 : testsignal1)
				for (double cx1 : testcx1)
					for (double cy1 : testcy1)
						for (double cz1 : testcz1)
							for (double[] w1 : testw1)
								for (double angle1 : testangle1)
								{
									a = createParameters(background, signal1, cx1, cy1, cz1, w1[0], w1[1], angle1);

									// Create y as a function we would want to move towards
									double[] a2 = a.clone();
									a2[targetParameter] *= 1.3;
									f1.initialise(a2);
									double[] data = new double[maxx * maxx];
									for (int i = 0; i < n; i++)
										data[i] = f1.eval(i);

									ff1 = new PoissonLikelihoodWrapper(f1, a, data, n, alpha);

									// Numerically solve gradient. 
									// Calculate the step size h to be an exact numerical representation
									final double xx = a[targetParameter];

									// Get h to minimise roundoff error
									double h = Precision.representableDelta(xx, h_);

									ff1.likelihood(getVariables(indices, a), dyda);

									// Evaluate at (x+h) and (x-h)
									a[targetParameter] = xx + h;
									double value2 = ff1.likelihood(getVariables(indices, a));

									a[targetParameter] = xx - h;
									double value3 = ff1.likelihood(getVariables(indices, a));

									double gradient = (value2 - value3) / (2 * h);
									boolean ok = Math.signum(gradient) == Math.signum(dyda[gradientIndex]) ||
											Math.abs(gradient - dyda[gradientIndex]) < 0.1;
									//logf("[%s-%s]/2*%g : %g == %g\n", "" + value2, "" + value3, h, gradient,
									//		dyda[gradientIndex]);
									if (!ok)
										Assert.assertTrue(
												NAME[targetParameter] + ": " + gradient + " != " + dyda[gradientIndex],
												ok);
									ok = eq.almostEqualRelativeOrAbsolute(gradient, dyda[gradientIndex]);
									if (ok)
										count++;
									total++;

								}
		double p = (100.0 * count) / total;
		TestSettings.info("%s : %s = %d / %d (%.2f)\n", f1.getClass().getSimpleName(), NAME[targetParameter], count, total, p);
		Assert.assertTrue(NAME[targetParameter] + " fraction too low: " + p, p > threshold);
	}

	private double[] getVariables(int[] indices, double[] a)
	{
		double[] variables = new double[indices.length];
		for (int i = 0; i < indices.length; i++)
		{
			variables[i] = a[indices[i]];

		}
		return variables;
	}

	private int findGradientIndex(Gaussian2DFunction f, int targetParameter)
	{
		int i = f.findGradientIndex(targetParameter);
		Assert.assertTrue("Cannot find gradient index", i >= 0);
		return i;
	}

	double[] createParameters(double... args)
	{
		return args;
	}

	@Test
	public void cumulativeProbabilityIsOneWithIntegerData()
	{
		// Initialise for large observed count
		PoissonLikelihoodWrapper.likelihood(1, photons[photons.length - 1] * 2);

		for (double p : photons)
			cumulativeProbabilityIsOneWithIntegerData(p);
	}

	private void cumulativeProbabilityIsOneWithIntegerData(final double mu)
	{
		double p = 0;
		int x = 0;

		// Evaluate an initial range. 
		// Poisson will have mean mu with a variance mu. 
		// At large mu it is approximately normal so use 3 sqrt(mu) for the range added to the mean
		if (mu > 0)
		{
			int max = (int) Math.ceil(mu + 3 * Math.sqrt(mu));
			for (; x <= max; x++)
			{
				final double pp = PoissonLikelihoodWrapper.likelihood(mu, x);
				//System.out.printf("x=%d, p=%f\n", x, pp);
				p += pp;
			}
			if (p > 1.01)
				Assert.fail("P > 1: " + p);
		}

		// We have most of the probability density. 
		// Now keep evaluating up until no difference
		final double changeTolerance = 1e-6;
		for (;; x++)
		{
			final double pp = PoissonLikelihoodWrapper.likelihood(mu, x);
			//System.out.printf("x=%d, p=%f\n", x, pp);
			p += pp;
			if (pp / p < changeTolerance)
				break;
		}
		TestSettings.info("mu=%f, p=%f, max=%d\n", mu, p, x);
		TestAssert.assertEquals(1, p, 0.02, "mu=%f", mu);
	}

	@Test
	public void cumulativeProbabilityIsOneWithRealDataForCountAbove4()
	{
		//for (int i = photons.length; i-- > 0;)
		for (int i = 0; i < photons.length; i++)
			if (photons[i] >= 4)
				cumulativeProbabilityIsOneWithRealData(photons[i], maxRange[i] + 1);
	}

	private void cumulativeProbabilityIsOneWithRealData(final double mu, int max)
	{
		double p = 0;

		SimpsonIntegrator in = new SimpsonIntegrator();

		p = in.integrate(20000, new UnivariateFunction()
		{
			@Override
			public double value(double x)
			{
				return PoissonLikelihoodWrapper.likelihood(mu, x);
			}
		}, 0, max);

		TestSettings.info("mu=%f, p=%f\n", mu, p);
		TestAssert.assertEquals(1, p, 0.02, "mu=%f", mu);
	}

	@Test
	public void cumulativeProbabilityIsOneFromLikelihoodForCountAbove4()
	{
		for (int i = 0; i < photons.length; i++)
			if (photons[i] >= 4)
				cumulativeProbabilityIsOneFromLikelihood(photons[i]);
	}

	private void cumulativeProbabilityIsOneFromLikelihood(final double mu)
	{
		// Determine upper limit for a Poisson
		int limit = new PoissonDistribution(mu).inverseCumulativeProbability(0.999);

		// Expand to allow for the gain
		int n = (int) Math.ceil(limit / alpha);

		// Evaluate all values from zero to large n
		double[] k = SimpleArrayUtils.newArray(n, 0, 1.0);
		double[] a = new double[0];
		double[] g = new double[0];

		NonLinearFunction nlf = new NonLinearFunction()
		{
			@Override
			public void initialise(double[] a)
			{
			}

			@Override
			public int[] gradientIndices()
			{
				return new int[0];
			}

			@Override
			public double eval(int x, double[] dyda, double[] w)
			{
				return 0;
			}

			@Override
			public double eval(int x)
			{
				return mu;
			}

			@Override
			public double eval(int x, double[] dyda)
			{
				return mu;
			}

			@Override
			public boolean canComputeWeights()
			{
				return false;
			}

			@Override
			public double evalw(int x, double[] w)
			{
				return 0;
			}

			@Override
			public int getNumberOfGradients()
			{
				return 0;
			}
		};
		PoissonLikelihoodWrapper f = new PoissonLikelihoodWrapper(nlf, a, Arrays.copyOf(k, n), n, alpha);

		// Keep evaluating up until no difference
		final double changeTolerance = 1e-6;
		double total = 0;
		double p = 0;
		int i = 0;
		while (i < n)
		{
			double nll = f.computeLikelihood(i);
			double nll2 = f.computeLikelihood(g, i);
			i++;
			Assert.assertEquals("computeLikelihood(i)", nll, nll2, 1e-10);
			total += nll;
			double pp = FastMath.exp(-nll);
			//System.out.printf("mu=%f, o=%f, i=%d, pp=%f\n", mu, mu / alpha, i, pp);
			p += pp;
			if (p > 0.5 && pp / p < changeTolerance)
				break;
		}

		TestSettings.info("mu=%f, limit=%d, p=%f\n", mu, limit, p);
		TestAssert.assertEquals(1, p, 0.02, "mu=%f", mu);

		// Check the function can compute the same total
		f = new PoissonLikelihoodWrapper(nlf, a, k, i, alpha);
		double sum = f.computeLikelihood();
		double sum2 = f.computeLikelihood(g);
		double delta = total * 1e-10;
		Assert.assertEquals("computeLikelihood", total, sum, delta);
		Assert.assertEquals("computeLikelihood with gradient", total, sum2, delta);
	}
}
