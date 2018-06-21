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
import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.utils.StoredDataStatistics;
import gdsc.test.TestSettings;
import gdsc.test.TestSettings.LogLevel;
import gdsc.test.TestSettings.TestComplexity;

public class PoissonGaussianConvolutionFunctionTest
{
	double[] gain = PoissonGaussianFunctionTest.gain;
	double[] photons = PoissonGaussianFunctionTest.photons;
	double[] noise = PoissonGaussianFunctionTest.noise;

	@Test
	public void cumulativeProbabilityIsOneWithPDF()
	{
		for (double g : gain)
			for (double p : photons)
				for (double s : noise)
					cumulativeProbabilityIsOne(g, p, s, false);
	}

	@Test
	public void cumulativeProbabilityIsOneWithPMF()
	{
		for (double g : gain)
			for (double p : photons)
				for (double s : noise)
					cumulativeProbabilityIsOne(g, p, s, true);
	}

	@Test
	public void probabilityMatchesLogProbabilityWithPDF()
	{
		for (double g : gain)
			for (double p : photons)
				for (double s : noise)
					probabilityMatchesLogProbability(g, p, s, false);
	}

	@Test
	public void probabilityMatchesLogProbabilityWithPMF()
	{
		for (double g : gain)
			for (double p : photons)
				for (double s : noise)
					probabilityMatchesLogProbability(g, p, s, true);
	}

	private void cumulativeProbabilityIsOne(final double gain, final double mu, final double s, boolean computePMF)
	{
		double p2 = cumulativeProbability(gain, mu, s, computePMF);
		String msg = String.format("g=%f, mu=%f, s=%f, erf=%b", gain, mu, s, computePMF);
		Assert.assertEquals(msg, 1, p2, 0.02);
	}

	private double cumulativeProbability(final double gain, final double mu, final double s, boolean computePMF)
	{
		// Note: The input s parameter is pre-gain.
		final PoissonGaussianConvolutionFunction f = PoissonGaussianConvolutionFunction
				.createWithStandardDeviation(1.0 / gain, s * gain);
		f.setComputePMF(computePMF);

		//final PoissonGaussianConvolutionFunction f2 = PoissonGaussianConvolutionFunction.createWithStandardDeviation(1.0 / gain, mu*gain, s * gain);
		//f2.setUsePicardApproximation(usePicard);

		double p = 0;
		int min = 1;
		int max = 0;

		// Note: The input mu parameter is pre-gain.
		final double e = mu;

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
				//System.out.printf("x=%d, p=%g\n", x, pp);
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
			//System.out.printf("x=%d, p=%g\n", x, pp);
			p += pp;
			if (pp == 0 || pp / p < changeTolerance)
				break;
		}
		for (int x = max + 1;; x++)
		{
			max = x;
			final double pp = f.likelihood(x, e);
			//System.out.printf("x=%d, p=%g\n", x, pp);
			p += pp;
			if (pp == 0 || pp / p < changeTolerance)
				break;
		}

		double p2 = p;
		if (!computePMF)
		{
			// Do a formal integration if the PDF
			//if (p < 0.98 || p > 1.02)
			//	System.out.printf("g=%f, mu=%f, s=%f p=%f\n", gain, mu, s, p);
			UnivariateIntegrator in = new SimpsonIntegrator(1e-4, 1e-6, 4,
					SimpsonIntegrator.SIMPSON_MAX_ITERATIONS_COUNT);
			p2 = in.integrate(Integer.MAX_VALUE, new UnivariateFunction()
			{
				@Override
				public double value(double x)
				{
					return f.likelihood(x, e);
				}
			}, min, max);
		}

		if (p2 < 0.98 || p2 > 1.02)
			TestSettings.info("g=%f, mu=%f, s=%f p=%f  %f\n", gain, mu, s, p, p2);

		return p2;
	}

	private void probabilityMatchesLogProbability(final double gain, double mu, final double s, boolean computePMF)
	{
		// Note: The input s parameter is pre-gain.
		PoissonGaussianConvolutionFunction f = PoissonGaussianConvolutionFunction
				.createWithStandardDeviation(1.0 / gain, s * gain);
		f.setComputePMF(computePMF);

		// Evaluate an initial range. 
		// Gaussian should have >99% within +/- s
		// Poisson will have mean mu with a variance mu. 
		// At large mu it is approximately normal so use 3 sqrt(mu) for the range added to the mean
		int[] range = PoissonGaussianFunctionTest.getRange(gain, mu, s);
		int min = range[0];
		int max = range[1];
		// Note: The input mu parameter is pre-gain.
		final double e = mu;
		String msg = String.format("g=%f, mu=%f, s=%f, erf=%b", gain, mu, s, computePMF);
		for (int x = min; x <= max; x++)
		{
			final double p = f.likelihood(x, e);
			if (p == 0)
				continue;
			final double logP = f.logLikelihood(x, e);
			Assert.assertEquals(msg, Math.log(p), logP, 1e-3 * Math.abs(logP));
		}
	}

	@Test
	public void pdfFasterThanPMF()
	{
		TestSettings.assume(LogLevel.WARN, TestComplexity.MEDIUM);

		// Realistic CCD parameters for speed test
		double s = 7.16;
		double g = 3.1;

		PoissonGaussianConvolutionFunction f1 = PoissonGaussianConvolutionFunction.createWithStandardDeviation(1 / g,
				s);
		f1.setComputePMF(true);

		PoissonGaussianConvolutionFunction f2 = PoissonGaussianConvolutionFunction.createWithStandardDeviation(1 / g,
				s);
		f2.setComputePMF(false);

		RandomGenerator rg = TestSettings.getRandomGenerator();

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
			double[] sample = new double[1000];
			for (int i = 0; i < sample.length; i++)
			{
				final double p = rg.nextDouble();
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

		TestSettings.logSpeedTestResult(t2 < t1, "cdf  %d -> pdf  %d = %f x\n", t1, t2, (double) t1 / t2);
	}

	private long run(PoissonGaussianConvolutionFunction f, double[][] samples, double[] photons)
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

}
