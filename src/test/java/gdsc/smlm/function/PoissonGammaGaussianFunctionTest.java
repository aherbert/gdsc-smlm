package gdsc.smlm.function;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.apache.commons.math3.analysis.integration.UnivariateIntegrator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.utils.DoubleEquality;
import gdsc.core.utils.StoredDataStatistics;
import gdsc.smlm.function.PoissonGammaGaussianFunction.ConvolutionMode;

public class PoissonGammaGaussianFunctionTest
{
	// Noise is in Counts and gain is total gain.
	// This makes more sense when testing as the 
	// PoissonGammaGaussianFunction accepts 1/gain and noise as parameters.

	// TODO Fix these test conditions

	// Poisson-Gamma convolution sums to above 1 at lower gain.
	// due to the Dirac delta function.
	// Store the expected sum at different gain below 10 for testing.
	static double[] pgSum = new double[11];
	static
	{
		//// Compute
		//MathContext mc = new MathContext(4, RoundingMode.UP);
		//for (int g = 1; g <= 10; g++)
		//{
		//	double sum = 0;
		//	for (int c = 0;; c++)
		//	{
		//		double p = PoissonGammaFunction.poissonGamma(c, 1, g);
		//		sum += p;
		//		if (p / sum < 1e-6)
		//			break;
		//	}
		//	pgSum[g] = sum;
		//	BigDecimal bd = new BigDecimal(sum);
		//	System.out.printf("pgSum[%d] = %.3f;\n", g, bd.round(mc).doubleValue());
		//}

		// Compute the sum at expected photons=1. This produces 
		// the highest sum as the contribution from the Poisson-Gamma to c=0
		// will be the greatest.
		// These are rounded up to 3 d.p. provide a safer upper bound.
		pgSum[1] = 1.200;
		pgSum[2] = 1.096;
		pgSum[3] = 1.064;
		pgSum[4] = 1.047;
		pgSum[5] = 1.038;
		pgSum[6] = 1.032;
		pgSum[7] = 1.027;
		pgSum[8] = 1.024;
		pgSum[9] = 1.021;
		pgSum[10] = 1.019;
	}

	double[] photons = { 0, 0.25, 0.5, 1, 2, 4, 10, 100, 1000 };
	double[] highPhotons = { 5000, 10000 };
	double[] lowPhotons = { 1e-2, 1e-4, 1e-6 };
	double[] noise = { 0.3, 1, 3, 10 }; // in counts
	double[] totalGain = { 6.5, 45 };

	@Test
	public void cumulativeGaussianProbabilityIsCorrect()
	{
		for (double s : noise)
			for (double g : totalGain)
				cumulativeGaussianProbabilityIsCorrect(s, g);
	}

	private void cumulativeGaussianProbabilityIsCorrect(double s, double g)
	{
		// Read noise should be in proportion to the camera gain
		final PoissonGammaGaussianFunction f = new PoissonGammaGaussianFunction(1 / g, s);
		double range = 5 * s;
		int upper = (int) Math.ceil(range);
		int lower = (int) Math.floor(-range);
		SimpsonIntegrator in = new SimpsonIntegrator(1e-4, 1e-8, 3, 32);
		UnivariateFunction uf = new UnivariateFunction()
		{
			public double value(double x)
			{
				return f.gaussianPDF(x);
			}
		};
		for (int u = lower; u <= upper; u++)
		{
			double ux = u + 0.5;
			double lx = u - 0.5;
			double e = in.integrate(20000, uf, lx, ux);
			double o = f.gaussianCDF(ux) - f.gaussianCDF(lx);
			double o2 = f.gaussianCDF(lx, ux);
			Assert.assertEquals(e, o, e * 0.1);
			Assert.assertEquals(o, o2, o * 1e-6);
		}
	}

	// The Poisson-Gamma has a delta function at c=0. This causes problems
	// if not correctly integrated.

	@Test
	public void cumulativeProbabilityIsOneWithDiscreteCDFIntegration()
	{
		for (double p : photons)
			for (double s : noise)
				for (double g : totalGain)
				{
					cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.DISCRETE_PMF);
				}
	}

	@Test
	public void cumulativeProbabilityIsOneWithDiscretePDFIntegration()
	{
		for (double p : photons)
			for (double s : noise)
				for (double g : totalGain)
				{
					cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.DISCRETE_PDF);
				}
	}

	@Test
	public void cumulativeProbabilityIsOneWithApproximationAtGainAbove10()
	{
		for (double p : photons)
			for (double s : noise)
				for (double g : totalGain)
				{
					if (g < 10)
						continue;
					cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.APPROXIMATION);
				}
	}

	@Test
	public void cumulativeProbabilityIsOneWithSimpsonIntegrationAtGainAbove10()
	{
		for (double p : photons)
			for (double s : noise)
				for (double g : totalGain)
				{
					cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.SIMPSON_PDF);
				}
	}

	@Test
	public void cumulativeProbabilityIsOneWithLegendreGaussIntegrationAtGainAbove10()
	{
		for (double p : photons)
			for (double s : noise)
				for (double g : totalGain)
				{
					cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.LEGENDRE_GAUSS_PDF);
				}
	}

	@Test
	public void cumulativeProbabilityIsOneWithDiscreteCDFIntegrationAtHighPhotons()
	{
		for (double p : highPhotons)
			for (double s : noise)
				for (double g : totalGain)
				{
					cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.DISCRETE_PMF);
				}
	}

	@Test
	public void cumulativeProbabilityIsOneWithDiscretePDFIntegrationAtHighPhotons()
	{
		for (double p : highPhotons)
			for (double s : noise)
				for (double g : totalGain)
				{
					cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.DISCRETE_PDF);
				}
	}

	// Since the Poisson-Gamma approximation is poor at gain < 10 the
	// following methods that use integration as an approximation 
	// (with the Poisson-Gamma PMF x Gaussian PDF) underestimate the total sum. 

	@Test
	public void cumulativeProbabilityIsOneWithApproximationAtHighPhotonsAtGainAbove10()
	{
		for (double p : photons)
			for (double s : noise)
				for (double g : totalGain)
				{
					if (g < 10)
						continue;
					cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.APPROXIMATION);
				}
	}

	@Test
	public void cumulativeProbabilityIsOneWithSimpsonIntegrationAtHighPhotonsAtGainAbove10()
	{
		for (double p : highPhotons)
			for (double s : noise)
				for (double g : totalGain)
				{
					if (g < 10)
						continue;
					cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.SIMPSON_PDF);
				}
	}

	@Test
	public void cumulativeProbabilityIsOneWithLegendreGaussIntegrationAtHighPhotonsAtGainAbove10()
	{
		for (double p : highPhotons)
			for (double s : noise)
				for (double g : totalGain)
				{
					if (g < 10)
						continue;
					cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.LEGENDRE_GAUSS_PDF);
				}
	}

	@Test
	public void cumulativeProbabilityIsOneWithDiscreteCDFIntegrationAtLowPhotons()
	{
		for (double p : lowPhotons)
			for (double s : noise)
				for (double g : totalGain)
				{
					cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.DISCRETE_PMF);
				}
	}

	@Test
	public void cumulativeProbabilityIsOneWithDiscretePDFIntegrationAtLowPhotons()
	{
		for (double p : lowPhotons)
			for (double s : noise)
				for (double g : totalGain)
				{
					cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.DISCRETE_PDF);
				}
	}

	// Since the Poisson-Gamma approximation is poor at gain < 10 the
	// following methods that use integration as an approximation 
	// (with the Poisson-Gamma PMF x Gaussian PDF) underestimate the total sum. 

	@Test
	public void cumulativeProbabilityIsOneWithApproximationAtLowPhotonsAtGainAbove10()
	{
		for (double p : photons)
			for (double s : noise)
				for (double g : totalGain)
				{
					if (g < 10)
						continue;
					cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.APPROXIMATION);
				}
	}

	@Test
	public void cumulativeProbabilityIsOneWithSimpsonIntegrationAtLowPhotonsAtGainAbove10()
	{
		for (double p : lowPhotons)
			for (double s : noise)
				for (double g : totalGain)
				{
					if (g < 10)
						continue;
					cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.SIMPSON_PDF);
				}
	}

	@Test
	public void cumulativeProbabilityIsOneWithLegendreGaussIntegrationAtLowPhotonsAtGainAbove10()
	{
		for (double p : lowPhotons)
			for (double s : noise)
				for (double g : totalGain)
				{
					if (g < 10)
						continue;
					cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.LEGENDRE_GAUSS_PDF);
				}
	}

	@Test
	public void discretePDFIntegrationCloselyMatchesDiscreteCDFIntegration()
	{
		double[] e = closelyMatchesDiscreteCDFIntegration(0.34, ConvolutionMode.DISCRETE_PDF);
		System.out.printf("Discrete integration max error : rel = %g : abs = %g\n", e[0], e[1]);
	}

	@Test
	public void approximationCloselyMatchesDiscreteCDFIntegration()
	{
		double[] e = closelyMatchesDiscreteCDFIntegration(0.15, ConvolutionMode.APPROXIMATION);
		System.out.printf("Approximation max error : rel = %g : abs = %g\n", e[0], e[1]);
	}

	@Test
	public void simpsonIntegrationMatchesDiscreteCDFIntegration()
	{
		double[] e = closelyMatchesDiscreteCDFIntegration(0.25, ConvolutionMode.SIMPSON_PDF);
		System.out.printf("Simpson integration max error : rel = %g : abs = %g\n", e[0], e[1]);
	}

	@Test
	public void legedreGaussIntegrationMatchesDiscreteCDFIntegration()
	{
		double[] e = closelyMatchesDiscreteCDFIntegration(0.15, ConvolutionMode.LEGENDRE_GAUSS_PDF);
		System.out.printf("Simpson integration max error : rel = %g : abs = %g\n", e[0], e[1]);
	}

	@Test
	public void approximationFasterThanDiscretePDFIntegration()
	{
		fasterThan(ConvolutionMode.DISCRETE_PDF, ConvolutionMode.APPROXIMATION);
	}

	@Test
	public void discretePDFIntegrationFasterThanDiscretePDFIntegration()
	{
		fasterThan(ConvolutionMode.DISCRETE_PMF, ConvolutionMode.DISCRETE_PDF);
	}

	@Test
	public void simpsonIntegrationFasterThanDiscreteCDFIntegration()
	{
		fasterThan(ConvolutionMode.DISCRETE_PMF, ConvolutionMode.SIMPSON_PDF);
	}

	@Test
	public void simpsonIntegrationFasterThanLegendreGaussIntegration()
	{
		fasterThan(ConvolutionMode.LEGENDRE_GAUSS_PDF, ConvolutionMode.SIMPSON_PDF);
	}

	private void cumulativeProbabilityIsOne(final double mu, final double s, final double g,
			ConvolutionMode convolutionMode)
	{
		double p = cumulativeProbability(mu, s, g, convolutionMode);
		System.out.printf("%s : mu=%f, s=%f, g=%f, p=%f\n", getName(convolutionMode), mu, s, g, p);

		// Poisson-Gamma convolution approximation does not sum to 1 at lower gain 
		// so account for this during the test.
		double delta = 0.02;
		double upper = 1 + delta;
		double lower = 1 - delta;
		if (g < 10)
		{
			lower = pgSum[(int) g] - delta;
		}

		if (p < lower || p > upper)
			Assert.fail(String.format("mu=%f, s=%f, g=%f, p=%g", mu, s, g, p));
	}

	private double cumulativeProbability(final double mu, final double s, final double g,
			ConvolutionMode convolutionMode)
	{
		final PoissonGammaGaussianFunction f = new PoissonGammaGaussianFunction(1 / g, s);
		f.setConvolutionMode(convolutionMode);
		f.setMinimumProbability(0);
		double p = 0;
		int min = 1;
		int max = 0;

		// Evaluate an initial range. 
		// Gaussian should have >99% within +/- 3s
		// Poisson will have mean mu with a variance mu. 
		// At large mu it is approximately normal so use 3 sqrt(mu) for the range added to the mean
		if (mu > 0)
		{
			int[] range = PoissonGaussianFunctionTest.getRange(g, mu, s);
			min = range[0];
			max = range[1];
			for (int x = min; x <= max; x++)
			{
				final double pp = f.likelihood(x, mu);
				//System.out.printf("x=%d, p=%g\n", x, pp);
				p += pp;
			}
			//if (p > 1.01)
			//	Assert.fail("P > 1: " + p);
		}

		// We have most of the probability density. 
		// Now keep evaluating up and down until no difference
		final double changeTolerance = 1e-6;
		for (int x = min - 1;; x--)
		{
			min = x;
			final double pp = f.likelihood(x, mu);
			//System.out.printf("x=%d, p=%g\n", x, pp);
			p += pp;
			if (pp / p < changeTolerance)
				break;
		}
		for (int x = max + 1;; x++)
		{
			max = x;
			final double pp = f.likelihood(x, mu);
			//System.out.printf("x=%d, p=%g\n", x, pp);
			p += pp;
			if (pp / p < changeTolerance)
				break;
		}

		// This is a simple integral. Compute the full integral if necessary.
		if (p < 0.98 || p > 1.02)
		{
			// Do a formal integration
			UnivariateIntegrator in = new SimpsonIntegrator(1e-6, 1e-6, 4,
					SimpsonIntegrator.SIMPSON_MAX_ITERATIONS_COUNT);
			double pp = in.integrate(Integer.MAX_VALUE, new UnivariateFunction()
			{
				public double value(double x)
				{
					return f.likelihood(x, mu);
				}
			}, min, max);
			System.out.printf("%s : mu=%f, rn=%f, cg=%f, s=%f, g=%f, p=%g => %g\n", getName(convolutionMode), mu, s, g,
					s, g, p, pp);
			p = pp;
		}

		return p;
	}

	private double[] closelyMatchesDiscreteCDFIntegration(double error, ConvolutionMode convolutionMode)
	{
		//DoubleEquality eq = new DoubleEquality(error, 1e-7);
		double[] maxError = new double[2];
		for (double s : noise)
			for (double g : totalGain)
			{
				if (g < 10)
					continue;

				PoissonGammaGaussianFunction f1 = new PoissonGammaGaussianFunction(1 / g, s);
				f1.setConvolutionMode(ConvolutionMode.DISCRETE_PMF);
				f1.setMinimumProbability(0);

				PoissonGammaGaussianFunction f2 = new PoissonGammaGaussianFunction(1 / g, s);
				f2.setConvolutionMode(convolutionMode);
				f2.setMinimumProbability(0);

				for (double p : photons)
				{
					double pg = p * g;
					double min = pg * 0.5 - 5 * s;
					double max = 2 * pg;
					for (double x = min; x < max; x += 1)
					{
						double p1 = f1.likelihood(x, p);
						double p2 = f2.likelihood(x, p);

						double relativeError = DoubleEquality.relativeError(p1, p2);
						double absError = Math.abs(p1 - p2);
						boolean equal = relativeError <= error; //eq.almostEqualRelativeOrAbsolute(p1, p2);
						if (!equal)
						{
							// Ignore small probabilities
							if (p1 < 1e-3)
								continue;

							Assert.fail(String.format("s=%g, g=%g, p=%g, x=%g: %g != %g (%g)", s, g, p, x, p1, p2,
									relativeError));
						}
						if (maxError[0] < relativeError)
							maxError[0] = relativeError;
						if (maxError[1] < absError)
							maxError[1] = absError;
					}
				}
			}
		return maxError;
	}

	private void fasterThan(ConvolutionMode slow, ConvolutionMode fast)
	{
		//org.junit.Assume.assumeTrue(false);

		// Realistic EM-CCD parameters for speed test
		double s = 7.16;
		double g = 39.1;

		PoissonGammaGaussianFunction f1 = new PoissonGammaGaussianFunction(1 / g, s);
		f1.setConvolutionMode(slow);
		if (!slow.validAtBoundary())
			f1.setBoundaryConvolutionMode(ConvolutionMode.DISCRETE_PMF);

		PoissonGammaGaussianFunction f2 = new PoissonGammaGaussianFunction(1 / g, s);
		f2.setConvolutionMode(fast);
		if (!fast.validAtBoundary())
			f2.setBoundaryConvolutionMode(ConvolutionMode.DISCRETE_PMF);

		// Generate realistic data from the probability mass function
		double[][] samples = new double[photons.length][];
		for (int j = 0; j < photons.length; j++)
		{
			int start = (int) (4 * -s);
			int u = start;
			StoredDataStatistics stats = new StoredDataStatistics();
			while (stats.getSum() < 0.995)
			{
				double p = f1.likelihood(u, photons[j]);
				stats.add(p);
				if (u > 10 && p / stats.getSum() < 1e-6)
					break;
				u++;
			}

			// Generate cumulative probability
			double[] data = stats.getValues();
			for (int i = 1; i < data.length; i++)
				data[i] += data[i - 1];

			// Sample
			RandomGenerator rand = new Well19937c();
			double[] sample = new double[1000];
			for (int i = 0; i < sample.length; i++)
			{
				final double p = rand.nextDouble();
				int x = 0;
				while (x < data.length && data[x] < p)
					x++;
				sample[i] = x;
			}
			samples[j] = sample;
		}

		// Warm-up
		run(f1, samples, photons);
		run(f2, samples, photons);

		long t1 = 0;
		for (int i = 0; i < 5; i++)
			t1 += run(f1, samples, photons);

		long t2 = 0;
		for (int i = 0; i < 5; i++)
			t2 += run(f2, samples, photons);

		System.out.printf("%s  %d -> %s  %d = %f x\n", getName(f1), t1, getName(f2), t2, (double) t1 / t2);
	}

	private long run(PoissonGammaGaussianFunction f, double[][] samples, double[] photons)
	{
		long start = System.nanoTime();
		for (int j = 0; j < photons.length; j++)
		{
			final double p = photons[j];
			for (double x : samples[j])
				f.likelihood(x, p);
		}
		return System.nanoTime() - start;
	}

	private String getName(PoissonGammaGaussianFunction f)
	{
		return getName(f.getConvolutionMode());
	}

	private String getName(ConvolutionMode convolutionMode)
	{
		return convolutionMode.toString();
	}
}
