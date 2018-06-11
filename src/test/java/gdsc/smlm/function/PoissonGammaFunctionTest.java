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

import java.util.Arrays;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.apache.commons.math3.analysis.integration.UnivariateIntegrator;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.utils.DoubleEquality;
import gdsc.core.utils.Maths;
import gnu.trove.list.array.TDoubleArrayList;

public class PoissonGammaFunctionTest
{
	static double[] gain = { 6, 16, 30 }; // ADU/electron above 1
	static double[] photons = { 0.001, 0.1, 0.25, 0.5, 1, 2, 4, 10, 100, 1000 };

	@Test
	public void cumulativeProbabilityIsOneWithPMF()
	{
		for (double g : gain)
			for (double p : photons)
				cumulativeProbabilityIsOne(g, p, false);
	}

	@Test
	public void cumulativeProbabilityIsOneWithPDF()
	{
		for (double g : gain)
			for (double p : photons)
				cumulativeProbabilityIsOne(g, p, true);
	}

	@Test
	public void probabilityMatchesLogProbability()
	{
		for (double g : gain)
			for (double p : photons)
				probabilityMatchesLogProbability(g, p);
	}

	private void cumulativeProbabilityIsOne(final double gain, final double mu, boolean pdf)
	{
		double p2 = cumulativeProbability(gain, mu, pdf);
		String msg = String.format("g=%f, mu=%f, pdf=%b", gain, mu, pdf);

		if (pdf)
		{
			Assert.assertEquals(msg, 1, p2, 0.02);
		}
		else
		{
			// This is not actually a PMF but is a PDF so requires integration.
			// This only works when the mean is above 2 if the gain is low
			if (mu > 2 || gain > 20)
				Assert.assertEquals(msg, 1, p2, 0.02);
		}
	}

	private double cumulativeProbability(final double gain, final double mu, boolean pdf)
	{
		final PoissonGammaFunction f = PoissonGammaFunction.createWithAlpha(1.0 / gain);

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
			int[] range = PoissonGaussianFunctionTest.getRange(gain, mu, 0);
			min = range[0];
			max = range[1];
			for (int x = min; x <= max; x++)
			{
				final double pp = f.likelihood(x, e);
				//System.out.printf("x=%d, p=%g\n", x, pp);
				if (debug)
					System.out.printf("x=%d, p=%f\n", x, pp);
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
				System.out.printf("x=%d, p=%f\n", x, pp);
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
				System.out.printf("x=%d, p=%f\n", x, pp);
			p += pp;
			if (pp == 0 || pp / p < changeTolerance)
				break;
		}

		double p2 = p;
		if (pdf)
		{
			// Do a formal integration
			if (p < 0.98 || p > 1.02)
				System.out.printf("g=%f, mu=%f, p=%f\n", gain, mu, p);
			UnivariateIntegrator in = new SimpsonIntegrator(1e-6, 1e-6, 4,
					SimpsonIntegrator.SIMPSON_MAX_ITERATIONS_COUNT);
			p2 = in.integrate(Integer.MAX_VALUE, new UnivariateFunction()
			{
				@Override
				public double value(double x)
				{
					//return f.likelihood(x, e);
					return PoissonGammaFunction.poissonGammaN(x, mu, gain);
				}
			}, min, max);

			p2 += PoissonGammaFunction.dirac(mu);
		}

		//if (p2 < 0.98 || p2 > 1.02)
		System.out.printf("g=%f, mu=%f, p=%f  %f\n", gain, mu, p, p2);

		return p2;
	}

	private void probabilityMatchesLogProbability(final double gain, double mu)
	{
		PoissonGammaFunction f = PoissonGammaFunction.createWithAlpha(1.0 / gain);

		// Evaluate an initial range. 
		// Gaussian should have >99% within +/- s
		// Poisson will have mean mu with a variance mu. 
		// At large mu it is approximately normal so use 3 sqrt(mu) for the range added to the mean
		// Note: The input s parameter is after-gain so adjust.
		int[] range = PoissonGaussianFunctionTest.getRange(gain, mu, 0);
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
			Assert.assertEquals(String.format("g=%f, mu=%f", gain, mu), Math.log(p), logP, 1e-6 * Math.abs(logP));
		}
	}

	@Test
	public void canComputePoissonGammaGradientWithInteger()
	{
		for (int j = 0; j < gain.length; j++)
			for (int i = 0; i < photons.length; i++)
				canComputePoissonGammaGradient(gain[j], photons[i], false);
	}

	@Test
	public void canComputePoissonGammaGradientWithReal()
	{
		for (int j = 0; j < gain.length; j++)
			for (int i = 0; i < photons.length; i++)
				canComputePoissonGammaGradient(gain[j], photons[i], true);
	}

	@SuppressWarnings("unused")
	private void canComputePoissonGammaGradient(final double gain, final double mu, boolean nonInteger)
	{
		final double o = mu;
		double delta = 1e-3; // * o;
		double uo = o + delta;
		double lo = o - delta;
		double diff = uo - lo;

		// The numerical gradient is poor around the switch between the use of the 
		// Bessel function and the approximation. So just count the errors.
		int fail = 0;
		double sum = 0;

		int[] range = PoissonGaussianFunctionTest.getRange(gain, mu, 0);
		int min = Math.max(0, range[0]);
		int max = range[1];
		double[] dp_dt = new double[1];
		double[] dp_dt2 = new double[1];
		double step = (nonInteger) ? 0.5 : 1;

		// When using the approximation the gradients are not as accurate
		boolean approx = (2 * Math.sqrt(max * o / gain) > 709);
		double tol = approx ? 0.05 : 1e-3;

		TDoubleArrayList list = new TDoubleArrayList();
		if (min != 0)
			list.add(0);
		for (double x = min; x <= max; x += step)
			list.add(x);

		for (double x : list.toArray())
		{
			double p1 = PoissonGammaFunction.poissonGamma(x, o, gain);
			double p2 = PoissonGammaFunction.poissonGamma(x, o, gain, dp_dt);
			Assert.assertEquals(p1, p2, 0);

			// Check partial gradient matches
			double p3 = PoissonGammaFunction.poissonGammaPartial(x, o, gain, dp_dt2);
			Assert.assertEquals(p1, p3, 0);
			Assert.assertEquals(dp_dt[0] + p1, dp_dt2[0], Math.abs(dp_dt[0] + p1) * 1e-8);

			// Check no dirac gradient matches
			p3 = PoissonGammaFunction.poissonGammaN(x, o, gain, dp_dt2);
			if (x == 0)
			{
				double dirac = PoissonGammaFunction.dirac(o);
				// Add the dirac contribution
				p3 += dirac;
				dp_dt2[0] -= dirac;
				Assert.assertEquals(p1, p3, p1 * 1e-8);
				Assert.assertEquals(dp_dt[0], dp_dt2[0], Math.abs(dp_dt[0]) * 1e-8);
			}
			else
			{
				Assert.assertEquals(p1, p3, 0);
				Assert.assertEquals(dp_dt[0], dp_dt2[0], Math.abs(dp_dt[0]) * 1e-8);
			}

			double up = PoissonGammaFunction.poissonGamma(x, uo, gain);
			double lp = PoissonGammaFunction.poissonGamma(x, lo, gain);

			double eg = dp_dt[0];
			double g = (up - lp) / diff;
			double error = DoubleEquality.relativeError(g, eg);
			double ox = x / gain;
			//System.out.printf("g=%g, mu=%g, x=%g (ox=%g), p=%g  g=%g  %g  error=%g\n", gain, mu, x, ox, p1, g, eg,
			//		error);

			if (error > tol)
			{
				fail++;
				sum += error;
			}
		}

		double f = (double) fail / list.size();
		System.out.printf("g=%g, mu=%g, failures=%g, mean=%f\n", gain, mu, f, Maths.div0(sum, fail));
		if (approx)
		{
			Assert.assertTrue(f < 0.2);
		}
		else
		{
			Assert.assertTrue(f < 0.01);
		}
	}

	@Test
	public void canComputeC0delta()
	{
		TDoubleArrayList list = new TDoubleArrayList();
		for (int exp = -12; exp < 6; exp++)
		{
			list.add(Math.pow(10, exp * 0.5));
		}
		for (int x = 2; x < 10; x++)
		{
			list.add(x / 10.0);
		}
		for (double x = 11; x <= 20; x++)
		{
			list.add(x / 10.0);
		}
		double[] p = list.toArray();
		Arrays.sort(p);

		double m = 5;

		for (double x : p)
		{
			double dirac = PoissonGammaFunction.dirac(x);
			double p0 = PoissonGammaFunction.poissonGammaN(0, x, m);
			double p01 = PoissonGammaFunction.poissonGammaN(1e-10, x, m);

			System.out.printf("p=%g  Dirac=%s   p0=%s (dirac:p0=%s)   p01=%s  (p0:p01 = %s)\n", x, dirac, p0,
					dirac / p0,
					//gdsc.core.utils.DoubleEquality.relativeError(p0, dirac),
					p01, p0 / p01
			//gdsc.core.utils.DoubleEquality.relativeError(p0, p01)
			);
		}
	}
}
