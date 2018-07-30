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
package uk.ac.sussex.gdsc.smlm.function;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import gnu.trove.list.array.TDoubleArrayList;
import uk.ac.sussex.gdsc.smlm.math3.distribution.CustomPoissonDistribution;
import uk.ac.sussex.gdsc.test.TestLog;
import uk.ac.sussex.gdsc.test.junit5.ExtraAssertions;

@SuppressWarnings({ "unused", "javadoc" })
public class PoissonFunctionTest
{
	static double[] gain = { 0.25, 0.5, 0.7, 1, 1.5, 1.7, 2, 2.2, 4, 8, 16 };
	static double[] photons = { 0.25, 0.5, 1, 2, 4, 10, 100, 1000 };

	@Test
	public void cumulativeProbabilityIsOne()
	{
		for (int j = 0; j < gain.length; j++)
			for (int i = 0; i < photons.length; i++)
			{
				final int[] result = cumulativeProbabilityIsOne(gain[j], photons[i]);
				TestLog.debug("minRange[%d][%d] = %d;\n", j, i, result[0]);
				TestLog.debug("maxRange[%d][%d] = %d;\n", j, i, result[1]);
			}
	}

	private static int[] cumulativeProbabilityIsOne(final double gain, final double mu)
	{
		final double o = mu;

		final PoissonFunction f = new PoissonFunction(1.0 / gain);
		double p = 0;

		final TDoubleArrayList values = new TDoubleArrayList();

		double maxp = 0;
		int maxc = 0;

		// Evaluate an initial range.
		// Poisson will have mean mu with a variance mu.
		// At large mu it is approximately normal so use 3 sqrt(mu) for the range added to the mean

		final int[] range = getRange(gain, mu);
		int min = range[0];
		int max = range[1];
		for (int x = min; x <= max; x++)
		{
			final double pp = f.likelihood(x, o);
			//System.out.printf("x=%d, p=%f\n", x, pp);
			p += pp;
			values.add(pp);
			if (maxp < pp)
			{
				maxp = pp;
				maxc = x;
			}
		}
		if (p > 1.01)
			Assertions.fail("P > 1: " + p);

		// We have most of the probability density.
		// Now keep evaluating up and down until no difference
		final double changeTolerance = 1e-6;
		if (min > 0)
		{
			values.reverse();
			for (int x = min - 1; x >= 0; x--)
			{
				min = x;
				final double pp = f.likelihood(x, o);
				//System.out.printf("x=%d, p=%f\n", x, pp);
				p += pp;
				values.add(pp);
				if (maxp < pp)
				{
					maxp = pp;
					maxc = x;
				}
				if (pp == 0 || pp / p < changeTolerance)
					break;
			}
			values.reverse();
		}
		for (int x = max + 1;; x++)
		{
			max = x;
			final double pp = f.likelihood(x, o);
			//System.out.printf("x=%d, p=%f\n", x, pp);
			p += pp;
			values.add(pp);
			if (maxp < pp)
			{
				maxp = pp;
				maxc = x;
			}
			if (pp == 0 || pp / p < changeTolerance)
				break;
		}

		// Find the range for 99.5% of the sum
		final double[] h = values.toArray();
		// Find cumulative
		for (int i = 1; i < h.length; i++)
			h[i] += h[i - 1];
		int minx = 0, maxx = h.length - 1;
		while (h[minx + 1] < 0.0025)
			minx++;
		while (h[maxx - 1] > 0.9975)
			maxx--;

		minx += min;
		maxx += min;

		TestLog.info("g=%f, mu=%f, o=%f, p=%f, min=%d, %f @ %d, max=%d\n", gain, mu, o, p, minx, maxp, maxc, maxx);
		return new int[] { minx, maxx };
	}

	static int[] getRange(final double gain, final double mu)
	{
		// Evaluate an initial range.
		// Poisson will have mean mu with a variance mu.
		// At large mu it is approximately normal so use 3 sqrt(mu) for the range added to the mean
		final double range = Math.max(1, Math.sqrt(mu));
		final int min = Math.max(0, (int) Math.floor(gain * (mu - 3 * range)));
		final int max = (int) Math.ceil(gain * (mu + 3 * range));
		return new int[] { min, max };
	}

	@Test
	public void probabilityMatchesLogProbabilty()
	{
		for (int j = 0; j < gain.length; j++)
			for (int i = 0; i < photons.length; i++)
				probabilityMatchesLogProbabilty(gain[j], photons[i]);
	}

	private static void probabilityMatchesLogProbabilty(final double gain, final double mu)
	{
		final double o = mu;

		final PoissonFunction f = new PoissonFunction(1.0 / gain);
		final double p = 0;

		final int[] range = getRange(gain, mu);
		final int min = range[0];
		final int max = range[1];
		for (int x = min; x <= max; x++)
		{
			final double v1 = Math.log(f.likelihood(x, o));
			final double v2 = f.logLikelihood(x, o);

			ExtraAssertions.assertEqualsRelative(v1, v2, 1e-8, "g=%f, mu=%f, x=%d", gain, mu, x);
		}
	}

	@Test
	public void probabilityMatchesPoissonWithNoGain()
	{
		for (int i = 0; i < photons.length; i++)
			probabilityMatchesPoissonWithNoGain(photons[i]);
	}

	private static void probabilityMatchesPoissonWithNoGain(final double mu)
	{
		final double o = mu;

		final PoissonFunction f = new PoissonFunction(1.0);
		final CustomPoissonDistribution pd = new CustomPoissonDistribution(null, mu);

		final double p = 0;

		final int[] range = getRange(1, mu);
		final int min = range[0];
		final int max = range[1];
		for (int x = min; x <= max; x++)
		{
			final double v1 = f.likelihood(x, o);
			final double v2 = pd.probability(x);

			ExtraAssertions.assertEqualsRelative(v1, v2, 1e-8, "g=%f, mu=%f, x=%d", gain, mu, x);
		}
	}
}
