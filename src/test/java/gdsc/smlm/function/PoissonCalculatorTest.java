package gdsc.smlm.function;

import java.math.BigDecimal;
import java.math.MathContext;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.apache.commons.math3.analysis.integration.UnivariateIntegrator;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.apache.commons.math3.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.ij.Utils;
import gdsc.smlm.TestSettings;

public class PoissonCalculatorTest
{
	double[] photons = { 1, 1.5, 2, 2.5, 3, 4, 5, 7.5, 10, 100, 1000 };
	private int maxx = 10;

	double P_LIMIT = 0.999999;

	@Test
	public void canComputeLikelihoodForIntegerData()
	{
		for (double u : photons)
		{
			PoissonDistribution pd = new PoissonDistribution(u);
			for (int x = 0; x < 100; x++)
			{
				double e = pd.probability(x);
				double o = PoissonCalculator.likelihood(u, x);
				if (e > 1e-100)
					Assert.assertEquals(e, o, e * 1e-10);
				e = pd.logProbability(x);
				o = PoissonCalculator.logLikelihood(u, x);
				Assert.assertEquals(e, o, Math.abs(e) * 1e-10);
			}
		}
	}

	@Test
	public void cumulativeProbabilityIsOneWithRealDataForCountAbove4()
	{
		for (double mu : photons)
		{
			// Determine upper limit for a Poisson
			double max = new PoissonDistribution(mu).inverseCumulativeProbability(P_LIMIT);

			// Determine lower limit
			double sd = Math.sqrt(mu);
			double min = (int) Math.max(0, mu - 4 * sd);

			cumulativeProbabilityIsOneWithRealData(mu, min, max, mu >= 4);
		}
	}

	private void cumulativeProbabilityIsOneWithRealData(final double mu, double min, double max, boolean test)
	{
		double p = 0;

		UnivariateIntegrator in = new SimpsonIntegrator();

		p = in.integrate(20000, new UnivariateFunction()
		{
			public double value(double x)
			{
				double v;
				v = PoissonCalculator.likelihood(mu, x);
				//v = pgf.probability(x, mu);
				//System.out.printf("x=%f, v=%f\n", x, v);
				return v;
			}
		}, min, max);

		System.out.printf("mu=%f, p=%f\n", mu, p);
		if (test)
		{
			Assert.assertEquals(String.format("mu=%f", mu), P_LIMIT, p, 0.02);
		}
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
	}

	@Test
	public void canComputeLogLikelihoodRatio()
	{
		final double n2 = maxx * maxx * 0.5;
		// Functions must produce a strictly positive output so add background
		//@formatter:off
		canComputeLogLikelihoodRatio(new BaseNonLinearFunction("Quadratic")
		{
			public double eval(int x) {	return 0.1 + a[0] * (x-n2) * (x-n2); }
		});		
		canComputeLogLikelihoodRatio(new BaseNonLinearFunction("Gaussian")
		{
			public double eval(int x) {	return 0.1 + 100 * FastMath.exp(-0.5 * Math.pow(x - n2, 2) / (a[0] * a[0])); }
		});		
		//@formatter:on
	}

	private void canComputeLogLikelihoodRatio(BaseNonLinearFunction nlf)
	{
		System.out.println(nlf.name);

		int n = maxx * maxx;

		double[] a = new double[] { 1 };

		// Simulate Poisson process
		nlf.initialise(a);
		RandomDataGenerator rdg = new RandomDataGenerator(new Well19937c(30051977));
		double[] x = Utils.newArray(n, 0, 1.0);
		double[] u = new double[x.length];
		for (int i = 0; i < n; i++)
		{
			u[i] = nlf.eval(i);
			if (u[i] > 0)
				x[i] = rdg.nextPoisson(u[i]);
		}

		double ll = PoissonCalculator.logLikelihood(u, x);
		double mll = PoissonCalculator.maximumLogLikelihood(x);
		double llr = -2 * (ll - mll);
		double llr2 = PoissonCalculator.logLikelihoodRatio(u, x);
		System.out.printf("llr=%f, llr2=%f\n", llr, llr2);
		Assert.assertEquals("Log-likelihood ratio", llr, llr2, llr * 1e-10);

		double[] op = new double[x.length];
		for (int i = 0; i < n; i++)
			op[i] = PoissonCalculator.maximumLikelihood(x[i]);

		double max = Double.NEGATIVE_INFINITY;
		double maxa = 0;

		//TestSettings.setLogLevel(gdsc.smlm.TestSettings.LogLevel.DEBUG);

		for (int i = 5; i <= 15; i++)
		{
			a[0] = (double) i / 10;
			nlf.initialise(a);
			for (int j = 0; j < n; j++)
				u[j] = nlf.eval(j);

			ll = PoissonCalculator.logLikelihood(u, x);
			llr = PoissonCalculator.logLikelihoodRatio(u, x);
			BigDecimal product = new BigDecimal(1);
			double ll2 = 0;
			for (int j = 0; j < n; j++)
			{
				double p1 = PoissonCalculator.likelihood(u[j], x[j]);
				ll2 += Math.log(p1);
				double ratio = p1 / op[j];
				product = product.multiply(new BigDecimal(ratio));
			}
			llr2 = -2 * Math.log(product.doubleValue());
			double p = 1 - PoissonCalculator.computePValue(llr, n - 1);
			TestSettings.info("a=%f, ll=%f, ll2=%f, llr=%f, llr2=%f, product=%s, p=%f\n", a[0], ll, ll2, llr, llr2,
					product.round(new MathContext(4)).toString(), p);
			if (max < ll)
			{
				max = ll;
				maxa = a[0];
			}

			// Only value if the product could be computed. Low ratios cause it to becomes 
			// too small to store in a double.
			if (product.doubleValue() > 0)
			{
				Assert.assertEquals("Log-likelihood", ll, ll2, Math.abs(ll2) * 1e-10);
				Assert.assertEquals("Log-likelihood ratio", llr, llr2, Math.abs(llr) * 1e-10);
			}
		}

		Assert.assertEquals("max", 1, maxa, 0);
	}

	@Test
	public void cannotSubtractConstantBackgroundAndComputeLogLikelihoodRatio()
	{
		int n = maxx * maxx;
		final double n2 = n * 0.5;
		final double n3 = n * 0.33;
		final double n4 = n * 0.66;
		// Functions must produce a strictly positive output so add background
		//@formatter:off
		cannotSubtractConstantBackgroundAndComputeLogLikelihoodRatio(
		new BaseNonLinearFunction("Quadratic")
		{
			public double eval(int x) {	return 0.1 + a[0] * (x-n2) * (x-n2); }
		},
		new BaseNonLinearFunction("Quadratic")
		{
			public double eval(int x) {	return 0.2 + 0.5 * a[0] * (x-n3) * (x-n3); }
		},
		new BaseNonLinearFunction("Quadratic")
		{
			public double eval(int x) {	return 0.3 + 0.75 * a[0] * (x-n4) * (x-n4); }
		});		
		cannotSubtractConstantBackgroundAndComputeLogLikelihoodRatio(
		new BaseNonLinearFunction("Gaussian")
		{
			public double eval(int x) {	return 0.1 + 100 * FastMath.exp(-0.5 * Math.pow(x - n2, 2) / (a[0] * a[0])); }
		},
		new BaseNonLinearFunction("Gaussian")
		{
			public double eval(int x) {	return 0.2 + 50 * FastMath.exp(-0.5 * Math.pow(x - n3, 2) / (a[0] * a[0])); }
		},
		new BaseNonLinearFunction("Gaussian")
		{
			public double eval(int x) {	return 0.3 + 75 * FastMath.exp(-0.5 * Math.pow(x - n4, 2) / (a[0] * a[0])); }
		});		
		//@formatter:on
	}

	private void cannotSubtractConstantBackgroundAndComputeLogLikelihoodRatio(BaseNonLinearFunction nlf1,
			BaseNonLinearFunction nlf2, BaseNonLinearFunction nlf3)
	{
		System.out.println(nlf1.name);

		int n = maxx * maxx;

		double[] a = new double[] { 1 };

		// Simulate Poisson process of 3 combined functions
		nlf1.initialise(a);
		nlf2.initialise(a);
		nlf3.initialise(a);
		RandomDataGenerator rdg = new RandomDataGenerator(new Well19937c(30051977));
		double[] x = Utils.newArray(n, 0, 1.0);
		double[] u = new double[x.length];
		double[] b1 = new double[x.length];
		double[] b2 = new double[x.length];
		double[] b3 = new double[x.length];
		for (int i = 0; i < n; i++)
		{
			b1[i] = nlf1.eval(i);
			b2[i] = nlf2.eval(i);
			b3[i] = nlf3.eval(i);
			u[i] = b1[i] + b2[i] + b3[i];
			if (u[i] > 0)
				x[i] = rdg.nextPoisson(u[i]);
		}

		// x is the target data
		// b1 is the already computed background
		// b2 is the first function to add to the model
		// b3 is the second function to add to the model

		// Compute the LLR of adding b3 to b2 given we already have b1 modelling data x
		double[] b12 = add(b1, b2);
		double ll1a = PoissonCalculator.logLikelihood(b12, x);
		double ll2a = PoissonCalculator.logLikelihood(add(b12, b3), x);
		double llra = -2 * (ll1a - ll2a);
		System.out.printf("x|(a+b+c) ll1=%f, ll2=%f, llra=%f\n", ll1a, ll2a, llra);

		// Compute the LLR of adding b3 to b2 given we already have x minus b1
		x = subtract(x, b1);
		double ll1b = PoissonCalculator.logLikelihood(b2, x);
		double ll2b = PoissonCalculator.logLikelihood(add(b2, b3), x);
		double llrb = -2 * (ll1b - ll2b);
		System.out.printf("x-a|(b+c) : ll1=%f, ll2=%f, llrb=%f\n", ll1b, ll2b, llrb);

		System.out.printf("llr=%f (%g), llr2=%f (%g)\n", llra, PoissonCalculator.computePValue(llra, 1), llrb,
				PoissonCalculator.computePValue(llrb, 1));
		Assert.assertNotEquals("Log-likelihood ratio", llra, llrb, llra * 1e-10);

	}

	private double[] add(double[] a, double[] b)
	{
		double[] c = new double[a.length];
		for (int i = 0; i < a.length; i++)
			c[i] = a[i] + b[i];
		return c;
	}

	private double[] subtract(double[] a, double[] b)
	{
		double[] c = new double[a.length];
		for (int i = 0; i < a.length; i++)
			c[i] = a[i] - b[i];
		return c;
	}
}
