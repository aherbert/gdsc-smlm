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
package uk.ac.sussex.gdsc.smlm.fitting;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.sampling.distribution.DiscreteInverseCumulativeProbabilityFunction;
import org.apache.commons.rng.sampling.distribution.DiscreteSampler;
import org.apache.commons.rng.sampling.distribution.InverseTransformDiscreteSampler;
import org.junit.jupiter.api.Assertions;

import uk.ac.sussex.gdsc.test.TestComplexity;
import uk.ac.sussex.gdsc.test.TestLog;
import uk.ac.sussex.gdsc.test.TestSettings;
import uk.ac.sussex.gdsc.test.junit5.ExtraAssertions;
import uk.ac.sussex.gdsc.test.junit5.ExtraAssumptions;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;

@SuppressWarnings({ "javadoc" })
public class BinomialFitterTest
{
	int[] N = new int[] { 2, 3, 4, 6 };
	double[] P = new double[] { 0.3, 0.5, 0.7 };
	int TRIALS = 10;
	int FAILURES = (int) (0.3 * TRIALS);

	// Note: This test is slow so only one test is run by default.

	// This is the level for the tests for standard parameters
	TestComplexity optionalTestComplexity = TestComplexity.HIGH;
	// This is the default level for all the tests
	TestComplexity nonEssentialTestComplexity = TestComplexity.VERY_HIGH;

	@SeededTest
	public void canFitBinomialWithKnownNUsingLeastSquaresEstimator(RandomSeed seed)
	{
		ExtraAssumptions.assume(nonEssentialTestComplexity);
		final UniformRandomProvider rg = TestSettings.getRandomGenerator(seed.getSeed());
		final boolean zeroTruncated = false;
		final boolean maximumLikelihood = false;
		for (final int n : N)
			for (final double p : P)
				fitBinomial(rg, n, p, zeroTruncated, maximumLikelihood, n, n);
	}

	@SeededTest
	public void canFitBinomialWithKnownNUsingMaximumLikelihood(RandomSeed seed)
	{
		final UniformRandomProvider rg = TestSettings.getRandomGenerator(seed.getSeed());
		final boolean zeroTruncated = false;
		final boolean maximumLikelihood = true;
		if (TestSettings.allow(optionalTestComplexity))
			for (final int n : N)
				for (final double p : P)
					fitBinomial(rg, n, p, zeroTruncated, maximumLikelihood, n, n);
		else
		{
			// This is the default test
			final int n = 2;
			final double p = 0.5;
			fitBinomial(rg, n, p, zeroTruncated, maximumLikelihood, n, n);
		}
	}

	@SeededTest
	public void canFitBinomialWithUnknownNUsingLeastSquaresEstimator(RandomSeed seed)
	{
		ExtraAssumptions.assume(nonEssentialTestComplexity);
		final UniformRandomProvider rg = TestSettings.getRandomGenerator(seed.getSeed());
		final boolean zeroTruncated = false;
		final boolean maximumLikelihood = false;
		for (final int n : N)
			for (final double p : P)
				fitBinomial(rg, n, p, zeroTruncated, maximumLikelihood, 1, n);
	}

	@SeededTest
	public void canFitBinomialWithUnknownNUsingMaximumLikelihood(RandomSeed seed)
	{
		ExtraAssumptions.assume(optionalTestComplexity);
		final UniformRandomProvider rg = TestSettings.getRandomGenerator(seed.getSeed());
		final boolean zeroTruncated = false;
		final boolean maximumLikelihood = true;

		// TODO - Sort out how to fit unknown N using MLE.
		// The problem is that the model returns a p of zero when n>N and this results in a negative infinity likelihood

		for (final int n : N)
			for (final double p : P)
				fitBinomial(rg, n, p, zeroTruncated, maximumLikelihood, 1, n);
	}

	@SeededTest
	public void canFitZeroTruncatedBinomialWithKnownNUsingLeastSquaresEstimator(RandomSeed seed)
	{
		ExtraAssumptions.assume(nonEssentialTestComplexity);
		final UniformRandomProvider rg = TestSettings.getRandomGenerator(seed.getSeed());
		final boolean zeroTruncated = true;
		final boolean maximumLikelihood = false;
		for (final int n : N)
			for (final double p : P)
				fitBinomial(rg, n, p, zeroTruncated, maximumLikelihood, n, n);
	}

	@SeededTest
	public void canFitZeroTruncatedBinomialWithKnownNUsingMaximumLikelihood(RandomSeed seed)
	{
		ExtraAssumptions.assume(nonEssentialTestComplexity);
		final UniformRandomProvider rg = TestSettings.getRandomGenerator(seed.getSeed());
		final boolean zeroTruncated = true;
		final boolean maximumLikelihood = true;
		for (final int n : N)
			for (final double p : P)
				fitBinomial(rg, n, p, zeroTruncated, maximumLikelihood, n, n);
	}

	@SeededTest
	public void canFitZeroTruncatedBinomialWithUnknownNUsingLeastSquaresEstimator(RandomSeed seed)
	{
		ExtraAssumptions.assume(nonEssentialTestComplexity);
		final UniformRandomProvider rg = TestSettings.getRandomGenerator(seed.getSeed());
		final boolean zeroTruncated = true;
		final boolean maximumLikelihood = false;
		for (final int n : N)
			for (final double p : P)
				fitBinomial(rg, n, p, zeroTruncated, maximumLikelihood, 1, n);
	}

	@SeededTest
	public void canFitZeroTruncatedBinomialWithUnknownNUsingMaximumLikelihood(RandomSeed seed)
	{
		ExtraAssumptions.assume(nonEssentialTestComplexity);
		final UniformRandomProvider rg = TestSettings.getRandomGenerator(seed.getSeed());
		final boolean zeroTruncated = true;
		final boolean maximumLikelihood = false;
		for (final int n : N)
			for (final double p : P)
				fitBinomial(rg, n, p, zeroTruncated, maximumLikelihood, 1, n);
	}

	@SeededTest
	public void sameFitBinomialWithKnownNUsing_LSE_Or_MLE(RandomSeed seed)
	{
		ExtraAssumptions.assume(nonEssentialTestComplexity);
		final UniformRandomProvider rg = TestSettings.getRandomGenerator(seed.getSeed());
		final boolean zeroTruncated = false;
		for (final int n : N)
			for (final double p : P)
				fitBinomialUsing_LSE_Or_MLE(rg, n, p, zeroTruncated, n, n);
	}

	@SeededTest
	public void sameFitZeroTruncatedBinomialWithKnownNUsing_LSE_Or_MLE(RandomSeed seed)
	{
		ExtraAssumptions.assume(nonEssentialTestComplexity);
		final UniformRandomProvider rg = TestSettings.getRandomGenerator(seed.getSeed());
		final boolean zeroTruncated = true;
		for (final int n : N)
			for (final double p : P)
				fitBinomialUsing_LSE_Or_MLE(rg, n, p, zeroTruncated, n, n);
	}

	private void fitBinomial(UniformRandomProvider rg, int n, double p, boolean zeroTruncated,
			boolean maximumLikelihood, int minN, int maxN)
	{
		final BinomialFitter bf = new BinomialFitter(null);
		//BinomialFitter bf = new BinomialFitter(new ConsoleLogger());
		bf.setMaximumLikelihood(maximumLikelihood);

		log("Fitting (n=%d, p=%f)\n", n, p);
		int fail = 0;
		for (int i = 0; i < TRIALS; i++)
		{
			final int[] data = createData(rg, n, p, false);
			final double[] fit = bf.fitBinomial(data, minN, maxN, zeroTruncated);
			final int fittedN = (int) fit[0];
			final double fittedP = fit[1];
			log("  Fitted (n=%d, p=%f)\n", fittedN, fittedP);
			try
			{
				Assertions.assertEquals(n, fittedN, "Failed to fit n");
				Assertions.assertEquals(p, fittedP, 0.05, "Failed to fit p");
			}
			catch (final AssertionError e)
			{
				fail++;
				log("    " + e.getMessage() + "\n");
			}
		}
		ExtraAssertions.assertTrue(fail <= FAILURES, "Too many failures (n=%d, p=%f): %d", n, p, fail);
	}

	private void fitBinomialUsing_LSE_Or_MLE(UniformRandomProvider rg, int n, double p, boolean zeroTruncated, int minN,
			int maxN)
	{
		final BinomialFitter bf = new BinomialFitter(null);
		//BinomialFitter bf = new BinomialFitter(new ConsoleLogger());

		log("Fitting (n=%d, p=%f)\n", n, p);
		int fail = 0;
		int c1 = 0;
		for (int i = 0; i < TRIALS; i++)
		{
			final int[] data = createData(rg, n, p, false);
			bf.setMaximumLikelihood(false);
			final double[] fitLSE = bf.fitBinomial(data, minN, maxN, zeroTruncated);
			bf.setMaximumLikelihood(true);
			final double[] fitMLE = bf.fitBinomial(data, minN, maxN, zeroTruncated);

			final int n1 = (int) fitLSE[0];
			final double p1 = fitLSE[1];
			final int n2 = (int) fitMLE[0];
			final double p2 = fitMLE[1];

			log("  Fitted LSE (n=%d, p=%f) == MLE (n=%d, p=%f)\n", n1, p1, n2, p2);

			try
			{
				Assertions.assertEquals(n1, n2, "Failed to match n");
				Assertions.assertEquals(p1, p2, 0.05, "Failed to match p");
			}
			catch (final AssertionError e)
			{
				fail++;
				log("    " + e.getMessage() + "\n");
			}
			if (Math.abs(p1 - p) < Math.abs(p2 - p))
				c1++;
		}
		log("  Closest LSE %d, MLE %d\n", c1, TRIALS - c1);
		if (fail > FAILURES)
		{
			final String msg = String.format("Too many failures (n=%d, p=%f): %d", n, p, fail);
			Assertions.fail(msg);
		}
	}

	private static int[] createData(UniformRandomProvider rg, int n, double p, boolean zeroTruncated)
	{
		final BinomialDistribution bd = new BinomialDistribution(null, n, p);
		final DiscreteSampler sampler = new InverseTransformDiscreteSampler(rg,
				new DiscreteInverseCumulativeProbabilityFunction()
				{
					@Override
					public int inverseCumulativeProbability(double p)
					{
						return bd.inverseCumulativeProbability(p);
					}
				});
		
		final int[] data = new int[2000];
		if (zeroTruncated)
		{
			if (p <= 0)
				throw new RuntimeException("p must be positive");
			for (int i = 0; i < data.length; i++)
			{
				int count;
				do
					count = sampler.sample();
				while (count == 0);
				data[i] = count;
			}
		}
		else
			for (int i = 0; i < data.length; i++)
				data[i] = sampler.sample();
		return data;
	}

	void log(String format, Object... args)
	{
		TestLog.info(format, args);
	}
}
