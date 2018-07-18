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
package uk.ac.sussex.gdsc.smlm.function;

import java.math.BigDecimal;
import java.math.MathContext;
import java.util.Arrays;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.apache.commons.math3.analysis.integration.UnivariateIntegrator;
import org.apache.commons.math3.distribution.ExponentialDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.Precision;
import org.junit.Assert;
import org.junit.Test;

import gnu.trove.list.array.TDoubleArrayList;
import uk.ac.sussex.gdsc.core.data.DataException;
import uk.ac.sussex.gdsc.core.math.QuadraticUtils;
import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import uk.ac.sussex.gdsc.test.TestAssert;
import uk.ac.sussex.gdsc.test.TestLog;
import uk.ac.sussex.gdsc.test.TestSettings;

@SuppressWarnings({ "javadoc" })
public class SCMOSLikelihoodWrapperTest
{
	private final double[] photons = { 1, 1.5, 2, 2.5, 3, 4, 5, 7.5, 10, 100, 1000 };

	DoubleEquality eqPerDatum = new DoubleEquality(5e12, 0.01);
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
	private final double h_ = 0.01; //(double) (Math.pow(1e-3f, 1.0 / 3));

	private final int[] testx = new int[] { 4, 5, 6 };
	private final int[] testy = new int[] { 4, 5, 6 };
	// Do not test zero background since this is an edge case for the likelihood function
	private final double[] testbackground_ = new double[] { 0.1, 1, 10 };
	private final double[] testsignal1_ = new double[] { 15, 55, 105 };
	private final double[] testangle1_ = new double[] { Math.PI / 5, Math.PI / 3 };
	private final double[] testcx1_ = new double[] { 4.9, 5.3 };
	private final double[] testcy1_ = new double[] { 4.8, 5.2 };
	private final double[] testcz1_ = new double[] { -1.5, 1.0 };
	private final double[][] testw1_ = new double[][] { { 1.1, 1.4 }, { 1.1, 1.7 }, { 1.5, 1.2 }, { 1.3, 1.7 }, };

	private double[] testbackground, testsignal1, testangle1, testcx1, testcy1, testcz1;
	private double[][] testw1;

	private static int maxx = 10;

	// Simulate per pixel noise
	private static float VAR = 57.9f;
	private static float G = 2.2f;
	private static float G_SD = 0.2f;
	private static float O = 100f;

	private static float[] var;
	private static float[] g;
	private static float[] o;
	private static float[] sd;

	static
	{
		final int n = maxx * maxx;
		var = new float[n];
		g = new float[n];
		o = new float[n];
		sd = new float[n];
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		final PoissonDistribution pd = new PoissonDistribution(rg, O, PoissonDistribution.DEFAULT_EPSILON,
				PoissonDistribution.DEFAULT_MAX_ITERATIONS);
		final ExponentialDistribution ed = new ExponentialDistribution(rg, VAR,
				ExponentialDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY);
		for (int i = 0; i < n; i++)
		{
			o[i] = pd.sample();
			var[i] = (float) ed.sample();
			sd[i] = (float) Math.sqrt(var[i]);
			g[i] = (float) (G + rg.nextGaussian() * G_SD);
		}
	}

	@Test
	public void fitFixedComputesGradientPerDatum()
	{
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
		functionComputesGradientPerDatum(GaussianFunctionFactory.FIT_ELLIPTICAL);
	}

	@Test
	public void fitNBFixedComputesGradientPerDatum()
	{
		functionComputesGradientPerDatum(GaussianFunctionFactory.FIT_SIMPLE_NB_FIXED);
	}

	@Test
	public void fitNBCircleComputesGradientPerDatum()
	{
		functionComputesGradientPerDatum(GaussianFunctionFactory.FIT_SIMPLE_NB_CIRCLE);
	}

	@Test
	public void fitNBFreeCircleComputesGradientPerDatum()
	{
		functionComputesGradientPerDatum(GaussianFunctionFactory.FIT_SIMPLE_NB_FREE_CIRCLE);
	}

	@Test
	public void fitNBEllipticalComputesGradientPerDatum()
	{
		functionComputesGradientPerDatum(GaussianFunctionFactory.FIT_SIMPLE_NB_ELLIPTICAL);
	}

	private void functionComputesGradientPerDatum(int flags)
	{
		final Gaussian2DFunction f1 = GaussianFunctionFactory.create2D(1, maxx, maxx, flags, null);
		// Setup
		// Setup
		testbackground = testbackground_;
		testsignal1 = testsignal1_;
		testcx1 = testcx1_;
		testcy1 = testcy1_;
		testcz1 = testcz1_;
		testw1 = testw1_;
		testangle1 = testangle1_;
		if (!f1.evaluatesBackground())
			testbackground = new double[] { testbackground[0] };
		if (!f1.evaluatesSignal())
		 testsignal1 = new double[] { testsignal1[0] };

		if (!f1.evaluatesZ())
			testcz1 = new double[] { 0 };
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
			for (int i = 0; i < testw1.length; i++)
				testw1[i][1] = testw1[i][0];
		if (!f1.evaluatesAngle())
			testangle1 = new double[] { 0 };

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
		final int[] indices = f1.gradientIndices();
		final int gradientIndex = findGradientIndex(f1, targetParameter);
		final double[] dyda = new double[indices.length];
		double[] a;

		SCMOSLikelihoodWrapper ff1;

		final int n = maxx * maxx;
		int count = 0, total = 0;

		final RandomDataGenerator rdg = new RandomDataGenerator(TestSettings.getRandomGenerator());

		for (final double background : testbackground)
			for (final double signal1 : testsignal1)
				for (final double cx1 : testcx1)
					for (final double cy1 : testcy1)
						for (final double cz1 : testcz1)
							for (final double[] w1 : testw1)
								for (final double angle1 : testangle1)
								{
									a = createParameters(background, signal1, cx1, cy1, cz1, w1[0], w1[1], angle1);

									// Create y as a function we would want to move towards
									final double[] a2 = a.clone();
									a2[targetParameter] *= 1.1;
									f1.initialise(a2);
									final double[] data = new double[n];
									for (int i = 0; i < n; i++)
									{
										// Simulate sCMOS camera
										final double u = f1.eval(i);
										data[i] = rdg.nextPoisson(u) * g[i] + rdg.nextGaussian(o[i], sd[i]);
									}

									ff1 = new SCMOSLikelihoodWrapper(f1, a, data, n, var, g, o);

									// Numerically solve gradient.
									// Calculate the step size h to be an exact numerical representation
									final double xx = a[targetParameter];

									// Get h to minimise roundoff error
									final double h = Precision.representableDelta(xx, h_);

									for (final int x : testx)
										for (final int y : testy)
										{
											final int i = y * maxx + x;
											a[targetParameter] = xx;
											ff1.likelihood(getVariables(indices, a), dyda, i);

											// Evaluate at (x+h) and (x-h)
											a[targetParameter] = xx + h;
											final double value2 = ff1.likelihood(getVariables(indices, a), i);

											a[targetParameter] = xx - h;
											final double value3 = ff1.likelihood(getVariables(indices, a), i);

											final double gradient = (value2 - value3) / (2 * h);
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
		final double p = (100.0 * count) / total;
		TestLog.info("Per Datum %s : %s = %d / %d (%.2f)\n", f1.getClass().getSimpleName(), NAME[targetParameter],
				count, total, p);
		Assert.assertTrue(NAME[targetParameter] + " fraction too low per datum: " + p, p > 90);
	}

	@Test
	public void fitFixedComputesGradient()
	{
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
		// The elliptical function gradient evaluation is worse
		final DoubleEquality tmp = eq;
		eq = eqPerDatum;
		functionComputesGradient(GaussianFunctionFactory.FIT_ELLIPTICAL);
		eq = tmp;
	}

	@Test
	public void fitNBFixedComputesGradient()
	{
		functionComputesGradient(GaussianFunctionFactory.FIT_SIMPLE_NB_FIXED);
	}

	@Test
	public void fitNBCircleComputesGradient()
	{
		functionComputesGradient(GaussianFunctionFactory.FIT_SIMPLE_NB_CIRCLE);
	}

	@Test
	public void fitNBFreeCircleComputesGradient()
	{
		functionComputesGradient(GaussianFunctionFactory.FIT_SIMPLE_NB_FREE_CIRCLE);
	}

	@Test
	public void fitNBEllipticalComputesGradient()
	{
		// The elliptical function gradient evaluation is worse
		final DoubleEquality tmp = eq;
		eq = eqPerDatum;
		functionComputesGradient(GaussianFunctionFactory.FIT_SIMPLE_NB_ELLIPTICAL);
		eq = tmp;
	}

	private void functionComputesGradient(int flags)
	{
		final Gaussian2DFunction f1 = GaussianFunctionFactory.create2D(1, maxx, maxx, flags, null);
		// Setup
		testbackground = testbackground_;
		testsignal1 = testsignal1_;
		testcx1 = testcx1_;
		testcy1 = testcy1_;
		testcz1 = testcz1_;
		testw1 = testw1_;
		testangle1 = testangle1_;
		if (!f1.evaluatesBackground())
			testbackground = new double[] { testbackground[0] };
		if (!f1.evaluatesSignal())
		 testsignal1 = new double[] { testsignal1[0] };

		if (!f1.evaluatesZ())
			testcz1 = new double[] { 0 };
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
			for (int i = 0; i < testw1.length; i++)
				testw1[i][1] = testw1[i][0];
		if (!f1.evaluatesAngle())
			testangle1 = new double[] { 0 };

		final double fraction = 85;
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
		final int[] indices = f1.gradientIndices();
		final int gradientIndex = findGradientIndex(f1, targetParameter);
		final double[] dyda = new double[indices.length];
		double[] a;

		SCMOSLikelihoodWrapper ff1;

		final int n = maxx * maxx;
		int count = 0, total = 0;

		final RandomDataGenerator rdg = new RandomDataGenerator(TestSettings.getRandomGenerator());

		for (final double background : testbackground)
			for (final double signal1 : testsignal1)
				for (final double cx1 : testcx1)
					for (final double cy1 : testcy1)
						for (final double cz1 : testcz1)
							for (final double[] w1 : testw1)
								for (final double angle1 : testangle1)
								{
									a = createParameters(background, signal1, cx1, cy1, cz1, w1[0], w1[1], angle1);

									// Create y as a function we would want to move towards
									final double[] a2 = a.clone();
									a2[targetParameter] *= 1.3;
									f1.initialise(a2);
									final double[] data = new double[n];
									for (int i = 0; i < n; i++)
									{
										// Simulate sCMOS camera
										final double u = f1.eval(i);
										data[i] = rdg.nextPoisson(u) * g[i] + rdg.nextGaussian(o[i], sd[i]);
									}

									ff1 = new SCMOSLikelihoodWrapper(f1, a, data, n, var, g, o);

									// Numerically solve gradient.
									// Calculate the step size h to be an exact numerical representation
									final double xx = a[targetParameter];

									// Get h to minimise roundoff error
									final double h = Precision.representableDelta(xx, h_);

									ff1.likelihood(getVariables(indices, a), dyda);

									// Evaluate at (x+h) and (x-h)
									a[targetParameter] = xx + h;
									final double value2 = ff1.likelihood(getVariables(indices, a));

									a[targetParameter] = xx - h;
									final double value3 = ff1.likelihood(getVariables(indices, a));

									final double gradient = (value2 - value3) / (2 * h);
									boolean ok = Math.signum(gradient) == Math.signum(dyda[gradientIndex]) ||
											Math.abs(gradient - dyda[gradientIndex]) < 0.1;
									//logf("[%s-%s]/2*%g : %g == %g\n", "" + value2, "" + value3, h, gradient,
									//		dyda[gradientIndex]);
									if (!ok)
										TestAssert.fail(NAME[targetParameter] + ": " + gradient + " != " + dyda[gradientIndex]);
									ok = eq.almostEqualRelativeOrAbsolute(gradient, dyda[gradientIndex]);
									if (ok)
										count++;
									total++;

								}
		final double p = (100.0 * count) / total;
		TestLog.info("%s : %s = %d / %d (%.2f)\n", f1.getClass().getSimpleName(), NAME[targetParameter], count,
				total, p);
		TestAssert.assertTrue(p > threshold, "%s fraction too low: %s", NAME[targetParameter], p);
	}

	private static double[] getVariables(int[] indices, double[] a)
	{
		final double[] variables = new double[indices.length];
		for (int i = 0; i < indices.length; i++)
			variables[i] = a[indices[i]];
		return variables;
	}

	private static int findGradientIndex(Gaussian2DFunction f, int targetParameter)
	{
		final int i = f.findGradientIndex(targetParameter);
		Assert.assertTrue("Cannot find gradient index", i >= 0);
		return i;
	}

	double[] createParameters(double... args)
	{
		return args;
	}

	double P_LIMIT = 0.999999;

	@Test
	public void cumulativeProbabilityIsOneWithRealDataForCountAbove8()
	{
		for (final double mu : photons)
		{
			// Determine upper limit for a Poisson
			double max = new PoissonDistribution(mu).inverseCumulativeProbability(P_LIMIT);

			// Determine lower limit
			final double sd = Math.sqrt(mu);
			double min = (int) Math.max(0, mu - 4 * sd);

			// Map to observed values using the gain and offset
			max = max * G + O;
			min = min * G + O;

			cumulativeProbabilityIsOneWithRealData(mu, min, max, mu > 8);
		}
	}

	private void cumulativeProbabilityIsOneWithRealData(final double mu, double min, double max, boolean test)
	{
		double p = 0;

		// Test using a standard Poisson-Gaussian convolution
		//min = -max;
		//final PoissonGaussianFunction pgf = PoissonGaussianFunction.createWithVariance(1, 1, VAR);

		final UnivariateIntegrator in = new SimpsonIntegrator();

		p = in.integrate(20000, new UnivariateFunction()
		{
			@Override
			public double value(double x)
			{
				double v;
				v = SCMOSLikelihoodWrapper.likelihood(mu, VAR, G, O, x);
				//v = pgf.probability(x, mu);
				//TestLog.debug("x=%f, v=%f\n", x, v);
				return v;
			}
		}, min, max);

		//TestLog.debug("mu=%f, p=%f\n", mu, p);
		if (test)
			TestAssert.assertEquals(P_LIMIT, p, 0.02, "mu=%f", mu);
	}

	@Test
	public void instanceLikelihoodMatches()
	{
		for (final double mu : photons)
			instanceLikelihoodMatches(mu, mu > 8);
	}

	private void instanceLikelihoodMatches(final double mu, boolean test)
	{
		// Determine upper limit for a Poisson
		final int limit = new PoissonDistribution(mu).inverseCumulativeProbability(P_LIMIT);

		// Map to observed values using the gain and offset
		final double max = limit * G;

		final double step = 0.1;

		final int n = (int) Math.ceil(max / step);

		// Evaluate all values from (zero+offset) to large n
		final double[] k = SimpleArrayUtils.newArray(n, O, step);
		final double[] a = new double[0];
		final double[] gradient = new double[0];

		final float var[] = newArray(n, VAR);
		final float g[] = newArray(n, G);
		final float o[] = newArray(n, O);

		final NonLinearFunction nlf = new NonLinearFunction()
		{
			@Override
			public void initialise(double[] a)
				{
		// Ignore
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
		SCMOSLikelihoodWrapper f = new SCMOSLikelihoodWrapper(nlf, a, k, n, var, g, o);

		double total = 0, p = 0;
		double maxp = 0;
		int maxi = 0;
		for (int i = 0; i < n; i++)
		{
			final double nll = f.computeLikelihood(i);
			final double nll2 = f.computeLikelihood(gradient, i);
			final double nll3 = SCMOSLikelihoodWrapper.negativeLogLikelihood(mu, var[i], g[i], o[i], k[i]);
			total += nll;
			TestAssert.assertEqualsRelative("computeLikelihood @" + i, nll3, nll, 1e-10);
			TestAssert.assertEqualsRelative("computeLikelihood+gradient @" + i, nll3, nll2, 1e-10);
			final double pp = FastMath.exp(-nll);
			if (maxp < pp)
			{
				maxp = pp;
				maxi = i;
				//TestLog.debug("mu=%f, e=%f, k=%f, pp=%f\n", mu, mu * G + O, k[i], pp);
			}
			p += pp * step;
		}

		// Expected max of the distribution is the mode of the Poisson distribution.
		// This has two modes for integer input counts. We take the mean of those.
		// https://en.wikipedia.org/wiki/Poisson_distribution
		// Note that the shift of VAR/(G*G) is a constant applied to both the expected and
		// observed values and consequently cancels when predicting the max, i.e. we add
		// a constant count to the observed values and shift the distribution by the same
		// constant. We can thus compute the mode for the unshifted distribution.
		final double lambda = mu;
		final double mode1 = Math.floor(lambda);
		final double mode2 = Math.ceil(lambda) - 1;
		final double kmax = ((mode1 + mode2) * 0.5) * G + O; // Scale to observed values
		//TestLog.debug("mu=%f, p=%f, maxp=%f @ %f  (expected=%f  %f)\n", mu, p, maxp, k[maxi], kmax, kmax - k[maxi]);
		TestAssert.assertEqualsRelative("k-max", kmax, k[maxi], 1e-3);

		if (test)
			TestAssert.assertEquals(P_LIMIT, p, 0.02, "mu=%f", mu);

		// Check the function can compute the same total
		double sum, sum2;
		sum = f.computeLikelihood();
		sum2 = f.computeLikelihood(gradient);
		TestAssert.assertEqualsRelative("computeLikelihood", total, sum, 1e-10);
		TestAssert.assertEqualsRelative("computeLikelihood with gradient", total, sum2, 1e-10);

		// Check the function can compute the same total after duplication
		f = f.build(nlf, a);
		sum = f.computeLikelihood();
		sum2 = f.computeLikelihood(gradient);
		TestAssert.assertEqualsRelative("computeLikelihood", total, sum, 1e-10);
		TestAssert.assertEqualsRelative("computeLikelihood with gradient", total, sum2, 1e-10);
	}

	private static float[] newArray(int n, float val)
	{
		final float[] a = new float[n];
		Arrays.fill(a, val);
		return a;
	}

	private abstract class BaseNonLinearFunction implements NonLinearFunction
	{
		double[] a;
		String name;

		BaseNonLinearFunction(String name)
		{
			this.name = name;
		}

		@Override
		public void initialise(double[] a)
		{
			this.a = a;
		}

		@Override
		public int[] gradientIndices()
		{
			return new int[1];
		}

		@Override
		public double eval(int x, double[] dyda, double[] w)
		{
			return 0;
		}

		@Override
		public double eval(int x, double[] dyda)
		{
			return 0;
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
			return 1;
		}
	}

	@Test
	public void canComputePValue()
	{
		final double n2 = maxx * maxx * 0.5;
		//@formatter:off
		canComputePValue(new BaseNonLinearFunction("Linear")
		{
			@Override
			public double eval(int x) {	return a[0] * (x-n2); }
		});
		canComputePValue(new BaseNonLinearFunction("Quadratic")
		{
			@Override
			public double eval(int x) {	return a[0] * (x-n2) * (x-n2); }
		});
		canComputePValue(new BaseNonLinearFunction("Linear+C")
		{
			@Override
			public double eval(int x) {	return 10 * a[0] + (x-n2); }
		});
		canComputePValue(new BaseNonLinearFunction("Gaussian")
		{
			@Override
			public double eval(int x) {	return 100 * FastMath.exp(-0.5 * Math.pow(x - n2, 2) / (a[0] * a[0])); }
		});
		//@formatter:on
	}

	private static void canComputePValue(BaseNonLinearFunction nlf)
	{
		TestLog.infoln(nlf.name);

		final int n = maxx * maxx;

		final double[] a = new double[] { 1 };

		// Simulate sCMOS camera
		nlf.initialise(a);
		final RandomDataGenerator rdg = new RandomDataGenerator(TestSettings.getRandomGenerator());
		final double[] k = SimpleArrayUtils.newArray(n, 0, 1.0);
		for (int i = 0; i < n; i++)
		{
			double u = nlf.eval(i);
			if (u > 0)
				u = rdg.nextPoisson(u);
			k[i] = u * g[i] + rdg.nextGaussian(o[i], sd[i]);
		}

		final SCMOSLikelihoodWrapper f = new SCMOSLikelihoodWrapper(nlf, a, k, n, var, g, o);

		final double oll = f.computeObservedLikelihood();
		double oll2 = 0;
		final double[] op = new double[n];
		for (int j = 0; j < n; j++)
		{
			op[j] = SCMOSLikelihoodWrapper.likelihood((k[j] - o[j]) / g[j], var[j], g[j], o[j], k[j]);
			oll2 -= Math.log(op[j]);
		}
		TestLog.info("oll=%f, oll2=%f\n", oll, oll2);
		TestAssert.assertEqualsRelative("Observed Log-likelihood", oll2, oll, 1e-10);

		final TDoubleArrayList list = new TDoubleArrayList();
		final int imin = 5, imax = 15;
		for (int i = imin; i <= imax; i++)
		{
			a[0] = (double) i / 10;
			final double ll = f.likelihood(a);
			list.add(ll);
			final double llr = f.computeLogLikelihoodRatio(ll);
			BigDecimal product = new BigDecimal(1);
			double ll2 = 0;
			for (int j = 0; j < n; j++)
			{
				final double p1 = SCMOSLikelihoodWrapper.likelihood(nlf.eval(j), var[j], g[j], o[j], k[j]);
				ll2 -= Math.log(p1);
				final double ratio = p1 / op[j];
				product = product.multiply(new BigDecimal(ratio));
			}
			final double llr2 = -2 * Math.log(product.doubleValue());
			final double q = f.computeQValue(ll);
			TestLog.info("a=%f, ll=%f, ll2=%f, llr=%f, llr2=%f, product=%s, p=%f\n", a[0], ll, ll2, llr, llr2,
					product.round(new MathContext(4)).toString(), q);

			// Only value if the product could be computed. Low ratios cause it to becomes
			// too small to store in a double.
			if (product.doubleValue() > 0)
				TestAssert.assertEqualsRelative("Log-likelihood", llr, llr2, 1e-10);
		}

		// Find min using quadratic fit
		final double[] data = list.toArray();
		int i = SimpleArrayUtils.findMinIndex(data);
		final double mina = (double) (imin + i) / 10;
		double fita = mina;
		try
		{
			if (i == 0)
				i++;
			if (i == data.length - 1)
				i--;
			final int i1 = i - 1;
			final int i2 = i;
			final int i3 = i + 1;

			fita = QuadraticUtils.findMinMax((double) (imin + i1) / 10, data[i1], (double) (imin + i2) / 10, data[i2],
					(double) (imin + i3) / 10, data[i3]);
		}
		catch (final DataException e)
		{
			// Ignore
		}

		// Allow a tolerance as the random data may alter the p-value computation.
		// Should allow it to be less than 2 increment either side of the answer.
		TestLog.info("min fit = %g => %g\n", mina, fita);
		Assert.assertEquals("min", 1, fita, 0.199);
	}
}
