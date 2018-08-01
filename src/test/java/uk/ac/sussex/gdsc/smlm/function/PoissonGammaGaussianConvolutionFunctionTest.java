package uk.ac.sussex.gdsc.smlm.function;

import java.util.logging.Level;
import java.util.logging.Logger;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.BeforeAll;

import java.util.function.Supplier;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.apache.commons.math3.analysis.integration.UnivariateIntegrator;
import org.junit.jupiter.api.Test;

import uk.ac.sussex.gdsc.smlm.function.PoissonGammaGaussianFunction.ConvolutionMode;
import uk.ac.sussex.gdsc.test.TestLog;
import uk.ac.sussex.gdsc.test.junit5.ExtraAssertions;

@SuppressWarnings({ "javadoc" })
public class PoissonGammaGaussianConvolutionFunctionTest
{
    private static Logger logger;

    @BeforeAll
    public static void beforeAll()
    {
        logger = Logger.getLogger(PoissonGammaGaussianConvolutionFunctionTest.class.getName());
    }

    @AfterAll
    public static void afterAll()
    {
        logger = null;
    }

	static double[] gain = { 6, 30 }; // ADU/electron above 1
	static double[] photons = PoissonGaussianFunctionTest.photons;
	static double[] noise = { 1, 10 }; // ADUs

	@Test
	public void cumulativeProbabilityIsOne()
	{
		for (final double g : gain)
			for (final double p : photons)
				for (final double s : noise)
					cumulativeProbabilityIsOne(g, p, s);
	}

	@Test
	public void probabilityMatchesLogProbability()
	{
		for (final double g : gain)
			for (final double p : photons)
				for (final double s : noise)
					probabilityMatchesLogProbability(g, p, s);
	}

	private static void cumulativeProbabilityIsOne(final double gain, final double mu, final double s)
	{
		final double p2 = cumulativeProbability(gain, mu, s);
		// This only works when the mean is above 2 if the gain is low
		if (mu > 2 || gain > 20)
			ExtraAssertions.assertEquals(1, p2, 0.02, "g=%f, mu=%f, s=%f", gain, mu, s);
	}

	private static double cumulativeProbability(final double gain, final double mu, double s)
	{
		final PoissonGammaGaussianConvolutionFunction f = PoissonGammaGaussianConvolutionFunction
				.createWithStandardDeviation(1.0 / gain, s);

		final PoissonGammaGaussianFunction f2 = new PoissonGammaGaussianFunction(1.0 / gain, s);
		f2.setConvolutionMode(ConvolutionMode.DISCRETE_PDF);

		double p = 0;
		int min = 1;
		int max = 0;

		// Note: The input mu parameter is pre-gain.
		final double e = mu;

		final boolean debug = false;

		// Evaluate an initial range.
		// Gaussian should have >99% within +/- s
		// Poisson will have mean mu with a variance mu.
		// At large mu it is approximately normal so use 3 sqrt(mu) for the range added to the mean
		if (mu > 0)
		{
			// Note: The input s parameter is after-gain so adjust.
			final int[] range = PoissonGaussianFunctionTest.getRange(gain, mu, s / gain);
			min = range[0];
			max = range[1];
			for (int x = min; x <= max; x++)
			{
				final double pp = f.likelihood(x, e);
				//System.out.printf("x=%d, p=%g\n", x, pp);
				if (debug)
					System.out.printf("x=%d, p=%f   %f\n", x, pp, f2.likelihood(x, e));
				p += pp;
			}
			//if (p > 1.01)
			//	Assertions.fail("P > 1: " + p);
		}

		// We have most of the likelihood density.
		// Now keep evaluating up and down until no difference
		final double changeTolerance = 1e-6;
		for (int x = min - 1;; x--)
		{
			min = x;
			final double pp = f.likelihood(x, e);
			//System.out.printf("x=%d, p=%g\n", x, pp);
			if (debug)
				System.out.printf("x=%d, p=%f   %f\n", x, pp, f2.likelihood(x, e));
			p += pp;
			if (pp == 0 || pp / p < changeTolerance)
				break;
		}
		for (int x = max + 1;; x++)
		{
			max = x;
			final double pp = f.likelihood(x, e);
			//System.out.printf("x=%d, p=%g\n", x, pp);
			if (debug)
				System.out.printf("x=%d, p=%f   %f\n", x, pp, f2.likelihood(x, e));
			p += pp;
			if (pp == 0 || pp / p < changeTolerance)
				break;
		}

		if (p < 0.98 || p > 1.02)
			TestLog.fine(logger,"g=%f, mu=%f, s=%f p=%f\n", gain, mu, s, p);

		// Do a formal integration
		double p2 = 0;
		final UnivariateIntegrator in = new SimpsonIntegrator(1e-4, 1e-6, 4, SimpsonIntegrator.SIMPSON_MAX_ITERATIONS_COUNT);
		p2 = in.integrate(Integer.MAX_VALUE, new UnivariateFunction()
		{
			@Override
			public double value(double x)
			{
				return f.likelihood(x, e);
			}
		}, min, max);

		if (p2 < 0.98 || p2 > 1.02)
			TestLog.info(logger,"g=%f, mu=%f, s=%f p=%f  %f\n", gain, mu, s, p, p2);

		return p2;
	}

	private static void probabilityMatchesLogProbability(final double gain, double mu, double s)
	{
		final PoissonGammaGaussianConvolutionFunction f = PoissonGammaGaussianConvolutionFunction
				.createWithStandardDeviation(1.0 / gain, s);

		// Evaluate an initial range.
		// Gaussian should have >99% within +/- s
		// Poisson will have mean mu with a variance mu.
		// At large mu it is approximately normal so use 3 sqrt(mu) for the range added to the mean
		// Note: The input s parameter is after-gain so adjust.
		final int[] range = PoissonGaussianFunctionTest.getRange(gain, mu, s / gain);
		final int min = range[0];
		final int max = range[1];
		// Note: The input mu parameter is pre-gain.
		final double e = mu;
		final Supplier<String> msg = () -> String.format("g=%f, mu=%f, s=%f", gain, mu, s);
		for (int x = min; x <= max; x++)
		{
			final double p = f.likelihood(x, e);
			if (p == 0)
				continue;
			final double logP = f.logLikelihood(x, e);
			ExtraAssertions.assertEqualsRelative(Math.log(p), logP, 1e-3, msg);
		}
	}
}
