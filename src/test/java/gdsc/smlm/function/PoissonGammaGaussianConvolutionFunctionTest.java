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

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.apache.commons.math3.analysis.integration.UnivariateIntegrator;
import org.junit.Assert;
import org.junit.Test;

import gdsc.smlm.function.PoissonGammaGaussianFunction.ConvolutionMode;

public class PoissonGammaGaussianConvolutionFunctionTest
{
	static double[] gain = { 6, 16, 30 }; // ADU/electron above 1
	static double[] photons = PoissonGaussianFunctionTest.photons;
	static double[] noise = { 1, 2, 4, 8, 16 }; // ADUs

	@Test
	public void cumulativeProbabilityIsOne()
	{
		for (double g : gain)
			for (double p : photons)
				for (double s : noise)
					cumulativeProbabilityIsOne(g, p, s);
	}

	@Test
	public void probabilityMatchesLogProbability()
	{
		for (double g : gain)
			for (double p : photons)
				for (double s : noise)
					probabilityMatchesLogProbability(g, p, s);
	}

	private void cumulativeProbabilityIsOne(final double gain, final double mu, final double s)
	{
		double p2 = cumulativeProbability(gain, mu, s);
		String msg = String.format("g=%f, mu=%f, s=%f", gain, mu, s);
		// This only works when the mean is above 2 if the gain is low
		if (mu > 2 || gain > 20)
			Assert.assertEquals(msg, 1, p2, 0.02);
	}

	private double cumulativeProbability(final double gain, final double mu, double s)
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

		boolean debug = false;

		// Evaluate an initial range. 
		// Gaussian should have >99% within +/- s
		// Poisson will have mean mu with a variance mu. 
		// At large mu it is approximately normal so use 3 sqrt(mu) for the range added to the mean
		if (mu > 0)
		{
			// Note: The input s parameter is after-gain so adjust.
			int[] range = PoissonGaussianFunctionTest.getRange(gain, mu, s / gain);
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
			//	Assert.fail("P > 1: " + p);
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
			System.out.printf("g=%f, mu=%f, s=%f p=%f\n", gain, mu, s, p);

		// Do a formal integration
		double p2 = 0;
		UnivariateIntegrator in = new SimpsonIntegrator(1e-6, 1e-6, 4, SimpsonIntegrator.SIMPSON_MAX_ITERATIONS_COUNT);
		p2 = in.integrate(Integer.MAX_VALUE, new UnivariateFunction()
		{
			@Override
			public double value(double x)
			{
				return f.likelihood(x, e);
			}
		}, min, max);

		if (p2 < 0.98 || p2 > 1.02)
			System.out.printf("g=%f, mu=%f, s=%f p=%f  %f\n", gain, mu, s, p, p2);

		return p2;
	}

	private void probabilityMatchesLogProbability(final double gain, double mu, double s)
	{
		PoissonGammaGaussianConvolutionFunction f = PoissonGammaGaussianConvolutionFunction
				.createWithStandardDeviation(1.0 / gain, s);

		// Evaluate an initial range. 
		// Gaussian should have >99% within +/- s
		// Poisson will have mean mu with a variance mu. 
		// At large mu it is approximately normal so use 3 sqrt(mu) for the range added to the mean
		// Note: The input s parameter is after-gain so adjust.
		int[] range = PoissonGaussianFunctionTest.getRange(gain, mu, s / gain);
		int min = range[0];
		int max = range[1];
		// Note: The input mu parameter is pre-gain.
		final double e = mu;
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
