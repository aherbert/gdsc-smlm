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

import java.util.Arrays;

import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.ejml.data.DenseMatrix64F;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import org.opentest4j.AssertionFailedError;

import uk.ac.sussex.gdsc.core.utils.Maths;
import uk.ac.sussex.gdsc.core.utils.Random;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.GradientCalculator;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.GradientCalculatorFactory;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import uk.ac.sussex.gdsc.test.TestLog;
import uk.ac.sussex.gdsc.test.TestSettings;
import uk.ac.sussex.gdsc.test.junit5.ExtraAssertions;

@SuppressWarnings({ "javadoc" })
public class FisherInformationMatrixTest
{
	@Test
	public void canComputeCRLB()
	{
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		for (int n = 1; n < 10; n++)
			canComputeCRLB(rg, n, 0, true);
	}

	@Test
	public void canComputeCRLBWithZeros()
	{
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		for (int n = 2; n < 10; n++)
		{
			canComputeCRLB(rg, n, 1, true);
			canComputeCRLB(rg, n, n / 2, true);
		}
	}

	@Test
	public void canComputeCRLBWithReciprocal()
	{
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		for (int n = 1; n < 10; n++)
			canComputeCRLB(rg, n, 0, false);
	}

	@Test
	public void canComputeCRLBWithReciprocalWithZeros()
	{
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		for (int n = 2; n < 10; n++)
		{
			canComputeCRLB(rg, n, 1, false);
			canComputeCRLB(rg, n, n / 2, false);
		}
	}

	@Test
	public void inversionDoesNotMatchReciprocal()
	{
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		for (int n = 1; n < 10; n++)
		{
			final FisherInformationMatrix m = createFisherInformationMatrix(rg, n, 0);
			final double[] crlb = m.crlb();
			final double[] crlb2 = m.crlbReciprocal();
			// These increasingly do not match with increasing number of parameters.
			TestLog.info("%s =? %s\n", Arrays.toString(crlb), Arrays.toString(crlb2));
			if (n > 1)
				// Just do a sum so we have a test
				Assertions.assertThrows(AssertionFailedError.class, () -> {
					Assertions.assertEquals(Maths.sum(crlb), Maths.sum(crlb2));
				});
		}
	}

	private static double[] canComputeCRLB(RandomGenerator rg, int n, int k, boolean invert)
	{
		final FisherInformationMatrix m = createFisherInformationMatrix(rg, n, k);

		// Invert for CRLB
		final double[] crlb = (invert) ? m.crlb() : m.crlbReciprocal();
		TestLog.info("n=%d, k=%d : %s\n", n, k, Arrays.toString(crlb));
		ExtraAssertions.assertNotNull(crlb, "CRLB failed: n=%d, k=%d", n, k);
		return crlb;
	}

	private static FisherInformationMatrix createFisherInformationMatrix(RandomGenerator rg, int n, int k)
	{
		final int maxx = 10;
		final int size = maxx * maxx;
		final RandomDataGenerator rdg = new RandomDataGenerator(rg);

		// Use a real Gaussian function here to compute the Fisher information.
		// The matrix may be sensitive to the type of equation used.
		int npeaks = 1;
		Gaussian2DFunction f = createFunction(maxx, npeaks);
		while (f.getNumberOfGradients() < n)
		{
			npeaks++;
			f = createFunction(maxx, npeaks);
		}

		final double[] a = new double[1 + npeaks * Gaussian2DFunction.PARAMETERS_PER_PEAK];
		a[Gaussian2DFunction.BACKGROUND] = rdg.nextUniform(1, 5);
		for (int i = 0, j = 0; i < npeaks; i++, j += Gaussian2DFunction.PARAMETERS_PER_PEAK)
		{
			a[j + Gaussian2DFunction.SIGNAL] = rdg.nextUniform(100, 300);
			// Non-overlapping peaks otherwise the CRLB are poor
			a[j + Gaussian2DFunction.X_POSITION] = rdg.nextUniform(2 + i * 2, 4 + i * 2);
			a[j + Gaussian2DFunction.Y_POSITION] = rdg.nextUniform(2 + i * 2, 4 + i * 2);
			a[j + Gaussian2DFunction.X_SD] = rdg.nextUniform(1.5, 2);
			a[j + Gaussian2DFunction.Y_SD] = rdg.nextUniform(1.5, 2);
		}
		f.initialise(a);

		final GradientCalculator c = GradientCalculatorFactory.newCalculator(f.getNumberOfGradients());
		double[][] I = c.fisherInformationMatrix(size, a, f);

		//TestLog.debug("n=%d, k=%d, I=\n", n, k);
		//for (int i = 0; i < I.length; i++)
		//	TestLog.debugln(Arrays.toString(I[i]));

		// Reduce to the desired size
		I = Arrays.copyOf(I, n);
		for (int i = 0; i < n; i++)
			I[i] = Arrays.copyOf(I[i], n);

		// Zero selected columns
		if (k > 0)
		{
			final int[] zero = Random.sample(k, n, rg); // new RandomDataGenerator(randomGenerator).nextPermutation(n, k);
			for (final int i : zero)
				for (int j = 0; j < n; j++)
					I[i][j] = I[j][i] = 0;
		}

		//TestLog.debug("n=%d, k=%d\n", n, k);
		//for (int i = 0; i < n; i++)
		//	TestLog.debugln(Arrays.toString(I[i]));

		// Create matrix
		return new FisherInformationMatrix(I, 1e-3);
	}

	private static Gaussian2DFunction createFunction(int maxx, int npeaks)
	{
		final Gaussian2DFunction f = GaussianFunctionFactory.create2D(npeaks, maxx, maxx,
				GaussianFunctionFactory.FIT_ERF_FREE_CIRCLE, null);
		return f;
	}

	void log(String format, Object... args)
	{
		TestLog.info(format, args);
	}

	@Test
	public void canProduceSubset()
	{
		final int k = 5;
		final int n = 10;

		final RandomGenerator randomGenerator = TestSettings.getRandomGenerator();
		final FisherInformationMatrix m = createRandomMatrix(randomGenerator, n);
		final DenseMatrix64F e = m.getMatrix();
		TestLog.infoln(e);

		for (int run = 1; run < 10; run++)
		{
			final int[] indices = Random.sample(k, n, randomGenerator);
			Arrays.sort(indices);
			final DenseMatrix64F o = m.subset(indices).getMatrix();
			TestLog.infoln(Arrays.toString(indices));
			TestLog.infoln(o);
			for (int i = 0; i < indices.length; i++)
				for (int j = 0; j < indices.length; j++)
					Assertions.assertEquals(e.get(indices[i], indices[j]), o.get(i, j));
		}
	}

	private static FisherInformationMatrix createRandomMatrix(RandomGenerator randomGenerator, int n)
	{
		final double[] data = new double[n * n];
		for (int i = 0; i < data.length; i++)
			data[i] = randomGenerator.nextDouble();
		return new FisherInformationMatrix(data, n);
	}

	@Test
	public void computeWithSubsetReducesTheCRLB()
	{
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		final Gaussian2DFunction f = createFunction(10, 1);
		final int perPeak = f.getGradientParametersPerPeak();
		// Create a matrix with 2 peaks + background
		final FisherInformationMatrix m = createFisherInformationMatrix(rg, 1 + 2 * perPeak, 0);
		// Subset each peak
		final int[] indices = SimpleArrayUtils.newArray(1 + perPeak, 0, 1);
		final FisherInformationMatrix m1 = m.subset(indices);
		for (int i = 1; i < indices.length; i++)
			indices[i] += perPeak;
		final FisherInformationMatrix m2 = m.subset(indices);

		//TestLog.debugln(m.getMatrix());
		//TestLog.debugln(m1.getMatrix());
		//TestLog.debugln(m2.getMatrix());

		final double[] crlb = m.crlb();
		final double[] crlb1 = m1.crlb();
		final double[] crlb2 = m2.crlb();
		final double[] crlbB = Arrays.copyOf(crlb1, crlb.length);
		System.arraycopy(crlb2, 1, crlbB, crlb1.length, perPeak);

		//TestLog.debugln(Arrays.toString(crlb));
		//TestLog.debugln(Arrays.toString(crlb1));
		//TestLog.debugln(Arrays.toString(crlb2));
		//TestLog.debugln(Arrays.toString(crlbB));

		// Removing the interaction between fit parameters lowers the bounds
		for (int i = 0; i < crlb.length; i++)
			Assertions.assertTrue(crlbB[i] < crlb[i]);
	}
}
