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

import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Assert;
import org.junit.Test;

import uk.ac.sussex.gdsc.test.TestAssert;
import uk.ac.sussex.gdsc.test.TestComplexity;
import uk.ac.sussex.gdsc.test.TestLog;
import uk.ac.sussex.gdsc.test.TestSettings;

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

	@Test
	public void canFitBinomialWithKnownNUsingLeastSquaresEstimator()
	{
		TestSettings.assume(nonEssentialTestComplexity);
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		final boolean zeroTruncated = false;
		final boolean maximumLikelihood = false;
		for (final int n : N)
			for (final double p : P)
				fitBinomial(rg, n, p, zeroTruncated, maximumLikelihood, n, n);
	}

	@Test
	public void canFitBinomialWithKnownNUsingMaximumLikelihood()
	{
		final RandomGenerator rg = TestSettings.getRandomGenerator();
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

	@Test
	public void canFitBinomialWithUnknownNUsingLeastSquaresEstimator()
	{
		TestSettings.assume(nonEssentialTestComplexity);
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		final boolean zeroTruncated = false;
		final boolean maximumLikelihood = false;
		for (final int n : N)
			for (final double p : P)
				fitBinomial(rg, n, p, zeroTruncated, maximumLikelihood, 1, n);
	}

	@Test
	public void canFitBinomialWithUnknownNUsingMaximumLikelihood()
	{
		TestSettings.assume(optionalTestComplexity);
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		final boolean zeroTruncated = false;
		final boolean maximumLikelihood = true;

		// TODO - Sort out how to fit unknown N using MLE.
		// The problem is that the model returns a p of zero when n>N and this results in a negative infinity likelihood

		for (final int n : N)
			for (final double p : P)
				fitBinomial(rg, n, p, zeroTruncated, maximumLikelihood, 1, n);
	}

	@Test
	public void canFitZeroTruncatedBinomialWithKnownNUsingLeastSquaresEstimator()
	{
		TestSettings.assume(nonEssentialTestComplexity);
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		final boolean zeroTruncated = true;
		final boolean maximumLikelihood = false;
		for (final int n : N)
			for (final double p : P)
				fitBinomial(rg, n, p, zeroTruncated, maximumLikelihood, n, n);
	}

	@Test
	public void canFitZeroTruncatedBinomialWithKnownNUsingMaximumLikelihood()
	{
		TestSettings.assume(nonEssentialTestComplexity);
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		final boolean zeroTruncated = true;
		final boolean maximumLikelihood = true;
		for (final int n : N)
			for (final double p : P)
				fitBinomial(rg, n, p, zeroTruncated, maximumLikelihood, n, n);
	}

	@Test
	public void canFitZeroTruncatedBinomialWithUnknownNUsingLeastSquaresEstimator()
	{
		TestSettings.assume(nonEssentialTestComplexity);
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		final boolean zeroTruncated = true;
		final boolean maximumLikelihood = false;
		for (final int n : N)
			for (final double p : P)
				fitBinomial(rg, n, p, zeroTruncated, maximumLikelihood, 1, n);
	}

	@Test
	public void canFitZeroTruncatedBinomialWithUnknownNUsingMaximumLikelihood()
	{
		TestSettings.assume(nonEssentialTestComplexity);
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		final boolean zeroTruncated = true;
		final boolean maximumLikelihood = false;
		for (final int n : N)
			for (final double p : P)
				fitBinomial(rg, n, p, zeroTruncated, maximumLikelihood, 1, n);
	}

	@Test
	public void sameFitBinomialWithKnownNUsing_LSE_Or_MLE()
	{
		TestSettings.assume(nonEssentialTestComplexity);
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		final boolean zeroTruncated = false;
		for (final int n : N)
			for (final double p : P)
				fitBinomialUsing_LSE_Or_MLE(rg, n, p, zeroTruncated, n, n);
	}

	@Test
	public void sameFitZeroTruncatedBinomialWithKnownNUsing_LSE_Or_MLE()
	{
		TestSettings.assume(nonEssentialTestComplexity);
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		final boolean zeroTruncated = true;
		for (final int n : N)
			for (final double p : P)
				fitBinomialUsing_LSE_Or_MLE(rg, n, p, zeroTruncated, n, n);
	}

	private void fitBinomial(RandomGenerator rg, int n, double p, boolean zeroTruncated, boolean maximumLikelihood,
			int minN, int maxN)
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
				Assert.assertEquals("Failed to fit n", n, fittedN);
				Assert.assertEquals("Failed to fit p", p, fittedP, 0.05);
			}
			catch (final AssertionError e)
			{
				fail++;
				log("    " + e.getMessage() + "\n");
			}
		}
		TestAssert.assertTrue(fail <= FAILURES, "Too many failures (n=%d, p=%f): %d", n, p, fail);
	}

	private void fitBinomialUsing_LSE_Or_MLE(RandomGenerator rg, int n, double p, boolean zeroTruncated, int minN,
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
				Assert.assertEquals("Failed to match n", n1, n2);
				Assert.assertEquals("Failed to match p", p1, p2, 0.05);
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
			Assert.fail(msg);
		}
	}

	private static int[] createData(RandomGenerator rg, int n, double p, boolean zeroTruncated)
	{
		final RandomDataGenerator rdg = new RandomDataGenerator(rg);
		final int[] data = new int[2000];
		if (zeroTruncated)
		{
			if (p <= 0)
				throw new RuntimeException("p must be positive");
			for (int i = 0; i < data.length; i++)
			{
				int count;
				do
					count = rdg.nextBinomial(n, p);
				while (count == 0);
				data[i] = count;
			}
		}
		else
			for (int i = 0; i < data.length; i++)
				data[i] = rdg.nextBinomial(n, p);
		return data;
	}

	void log(String format, Object... args)
	{
		TestLog.info(format, args);
	}
}
