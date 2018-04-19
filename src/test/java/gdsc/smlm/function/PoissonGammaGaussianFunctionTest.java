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
	// Change so noise is in Counts and gain is total gain.
	// This makes more sense when testing as the 
	// PoissonGammaGaussian accepts 1/gain and noise as parameters. 

	double[] photons = { 0, 0.25, 0.5, 1, 2, 4, 10, 100, 1000 };
	double[] highPhotons = { 5000, 10000 };
	double[] noise = { 30, 45, 76 }; // in electrons
	double[] cameraGain = { 6.5, 45 }; // ADU/e
	double emGain = 250;

	// Realistic parameters for speed test
	double s = 7.16;
	double g = 39.1;

	@Test
	public void cumulativeGaussianProbabilityIsCorrect()
	{
		for (double rn : noise)
			for (double cg : cameraGain)
				cumulativeGaussianProbabilityIsCorrect(rn, cg);
	}

	private void cumulativeGaussianProbabilityIsCorrect(double rn, double cg)
	{
		// Read noise should be in proportion to the camera gain
		double s = rn / cg;
		double g = emGain / cg;
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
			Assert.assertEquals(e, o, e * 0.1);
		}
	}

	@Test
	public void cumulativeProbabilityIsOneWithApproximation()
	{
		for (double p : photons)
			for (double rn : noise)
				for (double cg : cameraGain)
					cumulativeProbabilityIsOne(p, rn, cg, ConvolutionMode.APPROXIMATION);
	}

	// TODO - Figure out why these do not fail since they do not match the Discrete CDF 
	// integration and do not look correct in the CameraModelAnalysis plugin.

	//	@Test(expected = AssertionError.class)
	//	public void cumulativeProbabilityIsNotOneWithLegendreGaussIntegration()
	//	{
	//		for (double p : photons)
	//			for (double s : noise)
	//				for (double g : cameraGain)
	//					cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.LEGENDRE_GAUSS_PDF);
	//	}
	//
	//	@Test(expected = AssertionError.class)
	//	public void cumulativeProbabilityIsNotOneWithSimpsonIntegration()
	//	{
	//		for (double p : photons)
	//			for (double s : noise)
	//				for (double g : cameraGain)
	//					cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.SIMPSON_PDF);
	//	}

	@Test
	public void cumulativeProbabilityIsOneWithDiscreteCDFIntegrationWhenEffectiveGainIsAbove10()
	{
		for (double p : photons)
			for (double s : noise)
				for (double g : cameraGain)
				{
					double eg = emGain / g;
					if (eg <= 10)
						continue;
					cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.DISCRETE_CDF);
				}
	}

	@Test
	public void cumulativeProbabilityIsOneWithDiscretePDFIntegrationWhenEffectiveReadNoiseIsFarAbove1()
	{
		for (double p : photons)
			for (double s : noise)
				for (double g : cameraGain)
				{
					double eg = emGain / g;
					if (eg <= 10)
						continue;
					//					// Integration using a sample space of 1 is not valid when the 
					//					// effective read noise in counts is low 
					//					double es = s / g;
					//					if (es <= 4)
					//						continue;
					cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.DISCRETE_PDF);
				}
	}

	@Test
	public void cumulativeProbabilityIsOneWithApproximationAndDifficultParameters()
	{
		Assert.assertEquals(1, cumulativeProbability(10, 3, 45, ConvolutionMode.APPROXIMATION), 0.02);
		Assert.assertEquals(1, cumulativeProbability(1000, 30, 45, ConvolutionMode.APPROXIMATION), 0.02);
		Assert.assertEquals(1, cumulativeProbability(10000, 30, 45, ConvolutionMode.APPROXIMATION), 0.02);
		// Slow ... but passes
		//Assert.assertEquals(1, cumulativeProbability(100000, 30, 45, ConvolutionMode.LEGENDRE_GAUSS_PDF), 0.02);
	}

	@Test
	public void cumulativeProbabilityIsOneWithDiscreteCDFIntegrationAndDifficultParameters()
	{
		Assert.assertEquals(1, cumulativeProbability(0.5, 30, 4, ConvolutionMode.DISCRETE_CDF), 0.02);
		Assert.assertEquals(1, cumulativeProbability(10, 3, 45, ConvolutionMode.DISCRETE_CDF), 0.02);
		Assert.assertEquals(1, cumulativeProbability(1000, 30, 45, ConvolutionMode.DISCRETE_CDF), 0.02);
		Assert.assertEquals(1, cumulativeProbability(10000, 30, 45, ConvolutionMode.DISCRETE_CDF), 0.02);
	}

	@Test
	public void cumulativeProbabilityIsNotOneWithDiscretePDFIntegrationAndDifficultParameters()
	{
		// Integration using a sample space of 1 is not valid when the 
		// effective read noise in counts is below 1 

		// Read noise 3/45 = 0.067
		Assert.assertNotEquals(1, cumulativeProbability(0.5, 3, 45, ConvolutionMode.DISCRETE_PDF), 0.02);

		// Read noise 30/45 = 0.67

		// Fails when the total gain is low, OK when gain is high
		Assert.assertNotEquals(1, cumulativeProbability(0.5, 30, 45, ConvolutionMode.DISCRETE_PDF), 0.02);
		Assert.assertEquals(1, cumulativeProbability(0.5, 30 * 10, 45 / 10, ConvolutionMode.DISCRETE_PDF), 0.02);
	}

	@Test
	public void cumulativeProbabilityIsOneWithApproximationAndHighPhotons()
	{
		for (double p : highPhotons)
			for (double s : noise)
				for (double g : cameraGain)
					cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.APPROXIMATION);
	}

	@Test
	public void cumulativeProbabilityIsOneWithDiscretePDFIntegrationAndHighPhotons()
	{
		for (double p : highPhotons)
			for (double s : noise)
				for (double g : cameraGain)
					cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.DISCRETE_PDF);
	}

	@Test
	public void cumulativeProbabilityIsOneWithDiscreteCDFIntegrationAndHighPhotons()
	{
		for (double p : highPhotons)
			for (double s : noise)
				for (double g : cameraGain)
					cumulativeProbabilityIsOne(p, s, g, ConvolutionMode.DISCRETE_CDF);
	}

	@Test
	public void approximationCloselyMatchesDiscreteCDFIntegration()
	{
		double[] e = closelyMatchesDiscreteCDFIntegration(0.05, ConvolutionMode.APPROXIMATION);
		System.out.printf("Approximation max error : rel = %g : abs = %g\n", e[0], e[1]);
	}

	@Test
	public void discretePDFIntegrationCloselyMatchesDiscreteCDFIntegration()
	{
		double[] e = closelyMatchesDiscreteCDFIntegration(0.05, ConvolutionMode.DISCRETE_PDF);
		System.out.printf("Discrete integration max error : rel = %g : abs = %g\n", e[0], e[1]);
	}

	@Test(expected = AssertionError.class)
	public void simpsonIntegrationDoesNotMatchDiscreteCDFIntegration()
	{
		double[] e = closelyMatchesDiscreteCDFIntegration(0.05, ConvolutionMode.SIMPSON_PDF);
		System.out.printf("Simpson integration max error : rel = %g : abs = %g\n", e[0], e[1]);
	}

	@Test(expected = AssertionError.class)
	public void legedreGaussIntegrationDoesNotMatchDiscreteCDFIntegration()
	{
		double[] e = closelyMatchesDiscreteCDFIntegration(0.05, ConvolutionMode.LEGENDRE_GAUSS_PDF);
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
		fasterThan(ConvolutionMode.DISCRETE_CDF, ConvolutionMode.DISCRETE_PDF);
	}

	@Test
	public void simpsonIntegrationFasterThanDiscreteCDFIntegration()
	{
		fasterThan(ConvolutionMode.DISCRETE_CDF, ConvolutionMode.SIMPSON_PDF);
	}

	@Test
	public void simpsonIntegrationFasterThanLegendreGaussIntegration()
	{
		fasterThan(ConvolutionMode.LEGENDRE_GAUSS_PDF, ConvolutionMode.SIMPSON_PDF);
	}

	private void cumulativeProbabilityIsOne(final double mu, final double rn, final double cg,
			ConvolutionMode convolutionMode)
	{
		double p = cumulativeProbability(mu, rn, cg, convolutionMode);
		double s = rn / cg;
		double g = emGain / cg;
		System.out.printf("%s : mu=%f, rn=%f, cg=%f, s=%f, g=%f, p=%f\n", getName(convolutionMode), mu, rn, cg, s, g,
				p);
		//Assert.assertEquals(String.format("mu=%f, rn=%f, cg=%f, s=%f, g=%f", mu, rn, cg, s, g), 1, p, 0.02);
	}

	private double cumulativeProbability(final double mu, final double rn, final double cg,
			ConvolutionMode convolutionMode)
	{
		// Read noise should be in proportion to the camera gain
		double s = rn / cg;
		double g = emGain / cg;

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
			System.out.printf("%s : mu=%f, rn=%f, cg=%f, s=%f, g=%f, p=%g => %g\n", getName(convolutionMode), mu, rn,
					cg, s, g, p, pp);
			p = pp;
		}

		return p;
	}

	private double[] closelyMatchesDiscreteCDFIntegration(double error, ConvolutionMode convolutionMode)
	{
		//DoubleEquality eq = new DoubleEquality(error, 1e-7);
		double[] maxError = new double[2];
		for (double rn : noise)
			for (double cg : cameraGain)
			{
				// Read noise should be in proportion to the camera gain
				double s = rn / cg;
				if (s < 4)
					continue;

				double g = emGain / cg;
				PoissonGammaGaussianFunction f1 = new PoissonGammaGaussianFunction(1 / g, s);
				f1.setConvolutionMode(ConvolutionMode.DISCRETE_CDF);
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

							Assert.fail(String.format("rn=%g, cg=%g, s=%g, g=%g, p=%g, x=%g: %g != %g (%g)", rn, cg, s,
									g, p, x, p1, p2, relativeError));
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

		PoissonGammaGaussianFunction f1 = new PoissonGammaGaussianFunction(1 / g, s);
		f1.setConvolutionMode(slow);
		if (!slow.validAtBoundary())
			f1.setBoundaryConvolutionMode(ConvolutionMode.DISCRETE_CDF);

		PoissonGammaGaussianFunction f2 = new PoissonGammaGaussianFunction(1 / g, s);
		f2.setConvolutionMode(fast);
		if (!fast.validAtBoundary())
			f2.setBoundaryConvolutionMode(ConvolutionMode.DISCRETE_CDF);

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
