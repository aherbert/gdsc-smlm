package gdsc.smlm.function;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.apache.commons.math3.analysis.integration.UnivariateIntegrator;
import org.junit.Assert;
import org.junit.Test;

public class PoissonGaussianFunction2Test
{
	double[] gain = PoissonGaussianFunctionTest.gain;
	double[] photons = PoissonGaussianFunctionTest.photons;
	double[] noise = PoissonGaussianFunctionTest.noise;

	@Test
	public void cumulativeProbabilityIsOneWithPicard()
	{
		for (double g : gain)
			for (double p : photons)
				for (double s : noise)
					cumulativeProbabilityIsOne(g, p, s, true);
	}

	@Test
	public void cumulativeProbabilityIsOneWithPade()
	{
		for (double g : gain)
			for (double p : photons)
				for (double s : noise)
					cumulativeProbabilityIsOne(g, p, s, false);
	}

	@Test
	public void cumulativeProbabilityIsNotOneWhenMeanIsLowAndNoiseIsLow()
	{
		// The cumulative likelihood is poor for low mean and low noise. 
		// It can over-predict or under predict. The pattern of over/under is unknown. 
		// For example in the following:

		// OVER
		Assert.assertTrue(1.02 < cumulativeProbability(1.7, 0.25, 0.01, true));
		Assert.assertTrue(1.02 < cumulativeProbability(1.7, 0.25, 0.1, true));
		// OK
		Assert.assertEquals(1, cumulativeProbability(1.7, 0.25, 0.3, true), 0.02);
		// UNDER
		Assert.assertTrue(0.98 > cumulativeProbability(1.7, 0.25, 0.5, true));
		// OK
		Assert.assertEquals(1, cumulativeProbability(1.7, 0.25, 0.75, true), 0.02);

		// Fine with higher mean
		Assert.assertEquals(1, cumulativeProbability(1.7, 10, 0.01, true), 0.02);
		Assert.assertEquals(1, cumulativeProbability(1.7, 10, 0.1, true), 0.02);
		Assert.assertEquals(1, cumulativeProbability(1.7, 10, 0.3, true), 0.02);
		Assert.assertEquals(1, cumulativeProbability(1.7, 10, 0.5, true), 0.02);
		Assert.assertEquals(1, cumulativeProbability(1.7, 10, 0.75, true), 0.02);
	}

	@Test
	public void probabilityMatchesLogProbability()
	{
		for (double g : gain)
			for (double p : photons)
				for (double s : noise)
				{
					probabilityMatchesLogProbability(g, p, s, true);
					probabilityMatchesLogProbability(g, p, s, false);
				}
	}

	private void cumulativeProbabilityIsOne(final double gain, final double mu, final double s, final boolean usePicard)
	{
		double p2 = cumulativeProbability(gain, mu, s, usePicard);
		String msg = String.format("g=%f, mu=%f, s=%f", gain, mu, s);
		Assert.assertEquals(msg, 1, p2, 0.02);
	}

	private double cumulativeProbability(final double gain, final double mu, final double s, final boolean usePicard)
	{
		// Note: The input s parameter is pre-gain.
		final PoissonGaussianFunction2 f = PoissonGaussianFunction2.createWithStandardDeviation(1.0 / gain, s * gain);
		f.setUsePicardApproximation(usePicard);

		//final PoissonGaussianFunction f2 = PoissonGaussianFunction.createWithStandardDeviation(1.0 / gain, mu*gain, s * gain);
		//f2.setUsePicardApproximation(usePicard);

		double p = 0;
		int min = 1;
		int max = 0;

		// Note: The input mu parameter is pre-gain.
		final double e = mu * gain;

		// Evaluate an initial range. 
		// Gaussian should have >99% within +/- s
		// Poisson will have mean mu with a variance mu. 
		// At large mu it is approximately normal so use 3 sqrt(mu) for the range added to the mean
		if (mu > 0)
		{
			int[] range = PoissonGaussianFunctionTest.getRange(gain, mu, s);
			min = range[0];
			max = range[1];
			for (int x = min; x <= max; x++)
			{
				final double pp = f.likelihood(x, e);
				//System.out.printf("x=%d, p=%f   %f\n", x, pp);
				//System.out.printf("x=%d, p=%f   %f\n", x, pp, f2.probability(x));
				p += pp;
			}
			//if (p > 1.01)
			//	Assert.fail("P > 1: " + p);
		}

		// We have most of the likelihood density. 
		// Now keep evaluating up and down until no difference
		final double changeTolerance = 1e-6;
		for (int x = min - 1;; x--)
		{
			min = x;
			final double pp = f.likelihood(x, e);
			//System.out.printf("x=%d, p=%f\n", x, pp);
			p += pp;
			if (pp == 0 || pp / p < changeTolerance)
				break;
		}
		for (int x = max + 1;; x++)
		{
			max = x;
			final double pp = f.likelihood(x, e);
			//System.out.printf("x=%d, p=%f\n", x, pp);
			p += pp;
			if (pp == 0 || pp / p < changeTolerance)
				break;
		}

		// Do a formal integration
		double p2 = 0;
		UnivariateIntegrator in = new SimpsonIntegrator(1e-6, 1e-6, 4, SimpsonIntegrator.SIMPSON_MAX_ITERATIONS_COUNT);
		p2 = in.integrate(Integer.MAX_VALUE, new UnivariateFunction()
		{
			public double value(double x)
			{
				return f.likelihood(x, e);
			}
		}, min, max);

		if (p2 < 0.98 || p2 > 1.02)
			System.out.printf("g=%f, mu=%f, s=%f p=%f  %f\n", gain, mu, s, p, p2);

		return p2;
	}

	private void probabilityMatchesLogProbability(final double gain, double mu, final double s,
			final boolean usePicard)
	{
		// Note: The input s parameter is pre-gain.
		PoissonGaussianFunction2 f = PoissonGaussianFunction2.createWithStandardDeviation(1.0 / gain, s * gain);
		f.setUsePicardApproximation(usePicard);

		// Evaluate an initial range. 
		// Gaussian should have >99% within +/- s
		// Poisson will have mean mu with a variance mu. 
		// At large mu it is approximately normal so use 3 sqrt(mu) for the range added to the mean
		int[] range = PoissonGaussianFunctionTest.getRange(gain, mu, s);
		int min = range[0];
		int max = range[1];
		// Note: The input mu parameter is pre-gain.
		final double e = mu * gain;
		for (int x = min; x <= max; x++)
		{
			final double p = f.likelihood(x, e);
			if (p == 0)
				continue;
			final double logP = f.logLikelihood(x, e);
			Assert.assertEquals(String.format("g=%f, mu=%f, s=%f", gain, mu, s), Math.log(p), logP,
					1e-3 * Math.abs(logP));
		}
	}
}
