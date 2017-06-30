package gdsc.smlm.function;

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
import org.apache.commons.math3.random.Well19937c;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.Precision;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.ij.Utils;
import gdsc.core.utils.DoubleEquality;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.function.gaussian.GaussianFunctionFactory;

public class SCMOSLikelihoodWrapperTest
{
	private double[] photons = { 1, 1.5, 2, 2.5, 3, 4, 5, 7.5, 10, 100, 1000 };

	DoubleEquality eqPerDatum = new DoubleEquality(5e12, 0.01);
	DoubleEquality eq = new DoubleEquality(5e-3, 0.001);

	static String[] NAME;
	static
	{
		Gaussian2DFunction f1 = GaussianFunctionFactory.create2D(1, 1, 1, GaussianFunctionFactory.FIT_FIXED, null);
		NAME = new String[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
		for (int i = 0; i < 7; i++)
			NAME[i] = f1.getName(i);
	}

	// Compute as per Numerical Recipes 5.7.
	// Approximate error accuracy in single precision: Ef
	// Step size for derivatives:
	// h ~ (Ef)^(1/3) * xc
	// xc is the characteristic scale over which x changes, assumed to be 1 (not x as per NR since x is close to zero)
	private final double h_ = 0.01; //(double) (Math.pow(1e-3f, 1.0 / 3));

	private int[] testx = new int[] { 4, 5, 6 };
	private int[] testy = new int[] { 4, 5, 6 };
	// Do not test zero background since this is an edge case for the likelihood function
	private double[] testbackground_ = new double[] { 0.1, 1, 10 };
	private double[] testsignal1_ = new double[] { 15, 55, 105 };
	private double[] testangle1_ = new double[] { (double) (Math.PI / 5), (double) (Math.PI / 3) };
	private double[] testcx1_ = new double[] { 4.9, 5.3 };
	private double[] testcy1_ = new double[] { 4.8, 5.2 };
	private double[] testcz1_ = new double[] { -1.5, 1.0 };
	private double[][] testw1_ = new double[][] { { 1.1, 1.4 }, { 1.1, 1.7 }, { 1.5, 1.2 }, { 1.3, 1.7 }, };

	private double[] testbackground, testsignal1, testangle1, testcx1, testcy1, testcz1;
	private double[][] testw1;

	private int maxx = 10;

	// Simulate per pixel noise
	private float VAR = 57.9f;
	private float G = 2.2f;
	private float G_SD = 0.2f;
	private float O = 100f;

	private float[] var, g, o, sd;

	public SCMOSLikelihoodWrapperTest()
	{
		int n = maxx * maxx;
		var = new float[n];
		g = new float[n];
		o = new float[n];
		sd = new float[n];
		RandomGenerator rg = new Well19937c(30051977);
		PoissonDistribution pd = new PoissonDistribution(rg, O, PoissonDistribution.DEFAULT_EPSILON,
				PoissonDistribution.DEFAULT_MAX_ITERATIONS);
		ExponentialDistribution ed = new ExponentialDistribution(rg, VAR,
				ExponentialDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY);
		for (int i = 0; i < n; i++)
		{
			o[i] = (float) pd.sample();
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
		Gaussian2DFunction f1 = GaussianFunctionFactory.create2D(1, maxx, maxx, flags, null);
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

		SCMOSLikelihoodWrapper ff1;

		int n = maxx * maxx;
		int count = 0, total = 0;

		RandomDataGenerator rdg = new RandomDataGenerator(new Well19937c(30051977));

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
									double[] data = new double[n];
									for (int i = 0; i < n; i++)
									{
										// Simulate sCMOS camera
										double u = f1.eval(i);
										data[i] = rdg.nextPoisson(u) * g[i] + rdg.nextGaussian(o[i], sd[i]);
									}

									ff1 = new SCMOSLikelihoodWrapper(f1, a, data, n, var, g, o);

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
		logf("Per Datum %s : %s = %d / %d (%.2f)\n", f1.getClass().getSimpleName(), NAME[targetParameter], count, total,
				p);
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
		DoubleEquality tmp = eq;
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
		DoubleEquality tmp = eq;
		eq = eqPerDatum;
		functionComputesGradient(GaussianFunctionFactory.FIT_SIMPLE_NB_ELLIPTICAL);
		eq = tmp;
	}

	private void functionComputesGradient(int flags)
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

		SCMOSLikelihoodWrapper ff1;

		int n = maxx * maxx;
		int count = 0, total = 0;

		RandomDataGenerator rdg = new RandomDataGenerator(new Well19937c(30051977));

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
									double[] data = new double[n];
									for (int i = 0; i < n; i++)
									{
										// Simulate sCMOS camera
										double u = f1.eval(i);
										data[i] = rdg.nextPoisson(u) * g[i] + rdg.nextGaussian(o[i], sd[i]);
									}

									ff1 = new SCMOSLikelihoodWrapper(f1, a, data, n, var, g, o);

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
		logf("%s : %s = %d / %d (%.2f)\n", f1.getClass().getSimpleName(), NAME[targetParameter], count, total, p);
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

	void log(String message)
	{
		System.out.println(message);
	}

	void logf(String format, Object... args)
	{
		System.out.printf(format, args);
	}

	double P_LIMIT = 0.999999;

	@Test
	public void cumulativeProbabilityIsOneWithRealDataForCountAbove8()
	{
		for (double mu : photons)
		{
			// Determine upper limit for a Poisson
			double max = new PoissonDistribution(mu).inverseCumulativeProbability(P_LIMIT);

			// Determine lower limit
			double sd = Math.sqrt(mu);
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

		UnivariateIntegrator in = new SimpsonIntegrator();

		p = in.integrate(20000, new UnivariateFunction()
		{
			public double value(double x)
			{
				double v;
				v = SCMOSLikelihoodWrapper.likelihood(mu, VAR, G, O, x);
				//v = pgf.probability(x, mu);
				//System.out.printf("x=%f, v=%f\n", x, v);
				return v;
			}
		}, min, max);

		//System.out.printf("mu=%f, p=%f\n", mu, p);
		if (test)
		{
			Assert.assertEquals(String.format("mu=%f", mu), P_LIMIT, p, 0.02);
		}
	}

	@Test
	public void instanceLikelihoodMatches()
	{
		for (double mu : photons)
			instanceLikelihoodMatches(mu, mu > 8);
	}

	private void instanceLikelihoodMatches(final double mu, boolean test)
	{
		// Determine upper limit for a Poisson
		int limit = new PoissonDistribution(mu).inverseCumulativeProbability(P_LIMIT);

		// Map to observed values using the gain and offset
		double max = limit * G;

		double step = 0.1;

		int n = (int) Math.ceil(max / step);

		// Evaluate all values from (zero+offset) to large n
		double[] k = Utils.newArray(n, O, step);
		double[] a = new double[0];
		double[] gradient = new double[0];

		float var[] = newArray(n, VAR);
		float g[] = newArray(n, G);
		float o[] = newArray(n, O);

		NonLinearFunction nlf = new NonLinearFunction()
		{
			public void initialise(double[] a)
			{
			}

			public int[] gradientIndices()
			{
				return new int[0];
			}

			public double eval(int x, double[] dyda, double[] w)
			{
				return 0;
			}

			public double eval(int x)
			{
				return mu;
			}

			public double eval(int x, double[] dyda)
			{
				return mu;
			}

			public boolean canComputeWeights()
			{
				return false;
			}

			public double evalw(int x, double[] w)
			{
				return 0;
			}

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
			double nll = f.computeLikelihood(i);
			double nll2 = f.computeLikelihood(gradient, i);
			double nll3 = SCMOSLikelihoodWrapper.negativeLogLikelihood(mu, var[i], g[i], o[i], k[i]);
			total += nll;
			Assert.assertEquals("computeLikelihood @" + i, nll3, nll, nll * 1e-10);
			Assert.assertEquals("computeLikelihood+gradient @" + i, nll3, nll2, nll * 1e-10);
			double pp = FastMath.exp(-nll);
			if (maxp < pp)
			{
				maxp = pp;
				maxi = i;
				//System.out.printf("mu=%f, e=%f, k=%f, pp=%f\n", mu, mu * G + O, k[i], pp);
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
		double lambda = mu;
		double mode1 = Math.floor(lambda);
		double mode2 = Math.ceil(lambda) - 1;
		double kmax = ((mode1 + mode2) * 0.5) * G + O; // Scale to observed values
		//System.out.printf("mu=%f, p=%f, maxp=%f @ %f  (expected=%f  %f)\n", mu, p, maxp, k[maxi], kmax, kmax - k[maxi]);
		Assert.assertEquals("k-max", kmax, k[maxi], kmax * 1e-3);

		if (test)
		{
			Assert.assertEquals(String.format("mu=%f", mu), P_LIMIT, p, 0.02);
		}

		// Check the function can compute the same total
		double delta = total * 1e-10;
		double sum, sum2;
		sum = f.computeLikelihood();
		sum2 = f.computeLikelihood(gradient);
		Assert.assertEquals("computeLikelihood", total, sum, delta);
		Assert.assertEquals("computeLikelihood with gradient", total, sum2, delta);

		// Check the function can compute the same total after duplication
		f = f.build(nlf, a);
		sum = f.computeLikelihood();
		sum2 = f.computeLikelihood(gradient);
		Assert.assertEquals("computeLikelihood", total, sum, delta);
		Assert.assertEquals("computeLikelihood with gradient", total, sum2, delta);
	}

	private float[] newArray(int n, float val)
	{
		float[] a = new float[n];
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

		public void initialise(double[] a)
		{
			this.a = a;
		}

		public int[] gradientIndices()
		{
			return new int[1];
		}

		public double eval(int x, double[] dyda, double[] w)
		{
			return 0;
		}

		public double eval(int x, double[] dyda)
		{
			return 0;
		}

		public boolean canComputeWeights()
		{
			return false;
		}

		public double evalw(int x, double[] w)
		{
			return 0;
		}

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
			public double eval(int x) {	return a[0] * (x-n2); }
		});		
		canComputePValue(new BaseNonLinearFunction("Quadratic")
		{
			public double eval(int x) {	return a[0] * (x-n2) * (x-n2); }
		});		
		canComputePValue(new BaseNonLinearFunction("Linear+C")
		{
			public double eval(int x) {	return 10 * a[0] + (x-n2); }
		});		
		canComputePValue(new BaseNonLinearFunction("Gaussian")
		{
			public double eval(int x) {	return 100 * FastMath.exp(-0.5 * Math.pow(x - n2, 2) / (a[0] * a[0])); }
		});		
		//@formatter:on
	}

	private void canComputePValue(BaseNonLinearFunction nlf)
	{
		System.out.println(nlf.name);

		int n = maxx * maxx;

		double[] a = new double[] { 1 };

		// Simulate sCMOS camera
		nlf.initialise(a);
		RandomDataGenerator rdg = new RandomDataGenerator(new Well19937c(30051977));
		double[] k = Utils.newArray(n, 0, 1.0);
		for (int i = 0; i < n; i++)
		{
			double u = nlf.eval(i);
			if (u > 0)
				u = rdg.nextPoisson(u);
			k[i] = u * g[i] + rdg.nextGaussian(o[i], sd[i]);
		}

		SCMOSLikelihoodWrapper f = new SCMOSLikelihoodWrapper(nlf, a, k, n, var, g, o);

		double oll = f.computeObservedLikelihood();
		double oll2 = 0;
		double[] op = new double[n];
		for (int j = 0; j < n; j++)
		{
			op[j] = SCMOSLikelihoodWrapper.likelihood((k[j] - o[j]) / g[j], var[j], g[j], o[j], k[j]);
			oll2 -= Math.log(op[j]);
		}
		System.out.printf("oll=%f, oll2=%f\n", oll, oll2);
		Assert.assertEquals("Observed Log-likelihood", oll2, oll, oll2 * 1e-10);

		double min = Double.POSITIVE_INFINITY;
		double mina = 0;
		for (int i = 5; i <= 15; i++)
		{
			a[0] = (double) i / 10;
			double ll = f.likelihood(a);
			double llr = f.computeLogLikelihoodRatio(ll);
			BigDecimal product = new BigDecimal(1);
			double ll2 = 0;
			for (int j = 0; j < n; j++)
			{
				double p1 = SCMOSLikelihoodWrapper.likelihood(nlf.eval(j), var[j], g[j], o[j], k[j]);
				ll2 -= Math.log(p1);
				double ratio = p1 / op[j];
				product = product.multiply(new BigDecimal(ratio));
			}
			double llr2 = -2 * Math.log(product.doubleValue());
			double q = f.computeQValue(ll);
			System.out.printf("a=%f, ll=%f, ll2=%f, llr=%f, llr2=%f, product=%s, p=%f\n", a[0], ll, ll2, llr, llr2,
					product.round(new MathContext(4)).toString(), q);
			if (min > ll)
			{
				min = ll;
				mina = a[0];
			}

			// Only value if the product could be computed. Low ratios cause it to becomes 
			// too small to store in a double.
			if (product.doubleValue() > 0)
				Assert.assertEquals("Log-likelihood", llr, llr2, llr * 1e-10);
		}

		Assert.assertEquals("min", 1, mina, 0);
	}
}
