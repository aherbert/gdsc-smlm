package gdsc.smlm.function;

import gdsc.core.utils.DoubleEquality;
import gdsc.core.utils.StoredDataStatistics;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.junit.Assert;
import org.junit.Test;

public class PoissonGammaGaussianFunctionTest
{
	double[] photons = { 0, 0.25, 0.5, 1, 2, 4, 10, 100, 1000 };
	//double[] photons = { 100, 1000 };
	double[] noise = { 1.1, 2.3, 4.5, 8.7 };
	//double[] noise = { 7 };
	double[] gain = { 39.1 };

	// Realistic parmeters for speed test
	double s = 7.16;
	double g = 39.1;
	
	@Test
	public void cumulativeProbabilityIsOneWithApproximation()
	{
		for (double p : photons)
			for (double s : noise)
				for (double g : gain)
					cumulativeProbabilityIsOne(p, s, g, true, false);
	}

	@Test
	public void cumulativeProbabilityIsOneWithFullIntegration()
	{
		for (double p : photons)
			for (double s : noise)
				for (double g : gain)
					cumulativeProbabilityIsOne(p, s, g, false, false);
	}

	@Test
	public void cumulativeProbabilityIsOneWithSimpleIntegration()
	{
		for (double p : photons)
			for (double s : noise)
				for (double g : gain)
					cumulativeProbabilityIsOne(p, s, g, false, true);
	}

	@Test
	public void approximationCloselyMatchesFullIntegration()
	{
		double e = closelyMatchesFullIntegration(5e-2, true, false);
		System.out.println("Approximation max error = " + e);
	}

	@Test
	public void simpleIntegrationCloselyMatchesFullIntegration()
	{
		double e = closelyMatchesFullIntegration(1e-5, false, true);
		System.out.println("Simple integration max error = " + e);
	}

	@Test
	public void approximationFasterThanFullIntegration()
	{
		double g = 39.1;
		double s = 7.16;

		PoissonGammaGaussianFunction f1 = new PoissonGammaGaussianFunction(1 / g, s);
		f1.setUseApproximation(false);
		f1.setUseSimpleIntegration(false);

		PoissonGammaGaussianFunction f2 = new PoissonGammaGaussianFunction(1 / g, s);
		f2.setUseSimpleIntegration(false);
		f2.setUseApproximation(true);

		fasterThan(f1, f2);
	}

	@Test
	public void approximationFasterThanSimpleIntegration()
	{
		PoissonGammaGaussianFunction f1 = new PoissonGammaGaussianFunction(1 / g, s);
		f1.setUseApproximation(false);
		f1.setUseSimpleIntegration(true);

		PoissonGammaGaussianFunction f2 = new PoissonGammaGaussianFunction(1 / g, s);
		f2.setUseSimpleIntegration(false);
		f2.setUseApproximation(true);

		fasterThan(f1, f2);
	}

	@Test
	public void simpleIntegrationFasterThanFullIntegration()
	{
		PoissonGammaGaussianFunction f1 = new PoissonGammaGaussianFunction(1 / g, s);
		f1.setUseApproximation(false);
		f1.setUseSimpleIntegration(false);

		PoissonGammaGaussianFunction f2 = new PoissonGammaGaussianFunction(1 / g, s);
		f2.setUseSimpleIntegration(true);
		f2.setUseApproximation(false);

		fasterThan(f1, f2);
	}

	private void cumulativeProbabilityIsOne(final double mu, final double s, final double g, boolean useApproximation,
			boolean useSimpleIntegration)
	{
		PoissonGammaGaussianFunction f = new PoissonGammaGaussianFunction(1 / g, s);

		f.setUseApproximation(useApproximation);
		f.setUseSimpleIntegration(useSimpleIntegration);
		double p = 0;
		int min = 1;
		int max = 0;

		// Evaluate an initial range. 
		// Gaussian should have >99% within +/- 3s
		// Poisson will have mean mu with a variance mu. 
		// At large mu it is approximately normal so use 3 sqrt(mu) for the range added to the mean
		if (mu > 0)
		{
			min = (int) -Math.ceil(3 * s);
			max = (int) Math.ceil(mu + 3 * Math.sqrt(mu));
			for (int x = min; x <= max; x++)
			{
				final double pp = f.likelihood(x, mu);
				//System.out.printf("x=%d, p=%f\n", x, pp);
				p += pp;
			}
			if (p > 1.01)
				Assert.fail("P > 1: " + p);
		}

		// We have most of the probability density. 
		// Now keep evaluating up and down until no difference
		final double changeTolerance = 1e-6;
		for (int x = min - 1;; x--)
		{
			final double pp = f.likelihood(x, mu);
			//System.out.printf("x=%d, p=%f\n", x, pp);
			p += pp;
			if (pp / p < changeTolerance)
				break;
		}
		for (int x = max + 1;; x++)
		{
			final double pp = f.likelihood(x, mu);
			//System.out.printf("x=%d, p=%f\n", x, pp);
			p += pp;
			if (pp / p < changeTolerance)
				break;
		}
		System.out.printf("mu=%f, s=%f, p=%f\n", mu, s, p);
		Assert.assertEquals(String.format("mu=%f, s=%f", mu, s), 1, p, 0.02);
	}

	private double closelyMatchesFullIntegration(double error, boolean useApproximation, boolean useSimpleIntegration)
	{
		//DoubleEquality eq = new DoubleEquality(error, 1e-7);
		double maxError = 0;
		for (double s : noise)
			for (double g : gain)
			{
				PoissonGammaGaussianFunction f1 = new PoissonGammaGaussianFunction(1 / g, s);
				f1.setUseApproximation(false);
				f1.setUseSimpleIntegration(false);

				PoissonGammaGaussianFunction f2 = new PoissonGammaGaussianFunction(1 / g, s);
				f2.setUseApproximation(useApproximation);
				f2.setUseSimpleIntegration(useSimpleIntegration);

				for (double p : photons)
				{
					for (double x = p * 0.5 - 5 * s; x < 2 * p; x += 1)
					{
						double p1 = f1.likelihood(x, p);
						double p2 = f2.likelihood(x, p);

						double relativeError = DoubleEquality.relativeError(p1, p2);
						boolean equal = relativeError <= error; //eq.almostEqualRelativeOrAbsolute(p1, p2);
						if (!equal)
						{
							Assert.assertTrue(String.format("s=%f, g=%f, p=%f, x=%f: %f != %f (%f)", s, g, p, x, p1,
									p2, relativeError), equal);
						}
						if (maxError < relativeError)
							maxError = relativeError;
					}
				}
			}
		return maxError;
	}

	private void fasterThan(PoissonGammaGaussianFunction f1, PoissonGammaGaussianFunction f2)
	{
		// Generate realistic data from the probability mass function
		double[][] samples = new double[photons.length][];
		for (int j = 0; j < photons.length; j++)
		{
			int start = (int) (4 * -s);
			int u = start;
			StoredDataStatistics stats = new StoredDataStatistics();
			while (stats.getSum() < 0.995)
			{
				stats.add(f1.likelihood(u, photons[j]));
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
		if (f.isUseApproximation())
			return "Approximation";
		if (f.isUseSimpleIntegration())
			return "Simple integration";
		return "Full integration";
	}
}
