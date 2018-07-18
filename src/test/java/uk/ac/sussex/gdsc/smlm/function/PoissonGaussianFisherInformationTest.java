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

import org.junit.Assert;
import org.junit.Test;

import uk.ac.sussex.gdsc.test.TestLog;

@SuppressWarnings({ "javadoc" })
public class PoissonGaussianFisherInformationTest
{
	@Test
	public void canComputeFisherInformation()
	{
		//org.junit.Assume.assumeTrue(false);

		//double u;
		////u = Math.pow(10, -300);
		////u = 1.0 / Math.nextDown(Double.MAX_VALUE); // Smallest p with a non-infinite Fisher information
		//u = 1e-100;
		//PoissonGaussianFisherInformation f = new PoissonGaussianFisherInformation(0.20);
		//f.setMeanThreshold(Double.MAX_VALUE);
		//double I = f.getPoissonGaussianI(u);
		//double lower = f.getPoissonGaussianApproximationI(u);
		//double upper = PoissonFisherInformation.getPoissonI(u);
		//System.out.printf("s=%g u=%g I=%s  (%s - %s) alpha=%s\n", f.s, u, I, lower, upper, I / upper);
		//if (true)
		//	return;

		for (int i = 0; i < 4; i++)
		{
			final double s = (1 << i) * 0.25;
			canComputeFisherInformation(s);
		}
	}

	private static void canComputeFisherInformation(double s)
	{
		canComputeFisherInformation(new PoissonGaussianFisherInformation(s));
	}

	private static void canComputeFisherInformation(PoissonGaussianFisherInformation f)
	{
		f.setMeanThreshold(Double.MAX_VALUE);
		// Do not evaluate at very high mean. The switch to the approximation will occur
		// and the approximation is good.
		for (int exp = -20; exp < 6; exp++)
			canComputeFisherInformation(f, Math.pow(10, exp * 0.5));
	}

	private static void canComputeFisherInformation(PoissonGaussianFisherInformation f, double u)
	{
		final double I = f.getPoissonGaussianI(u);
		double lower = f.getPoissonGaussianApproximationI(u);
		double upper = PoissonFisherInformation.getPoissonI(u);
		//System.out.printf("s=%g u=%g I=%s  (%s - %s) alpha=%s\n", f.s, u, I, lower, upper, I / upper);
		// Allow a tolerance on the approximation at high mean.
		// The function does not compute the sum to infinity and so can underestimate
		// the value.
		if (u >= 10)
			lower *= 0.99;
		if (u >= 10)
			upper *= 1.01;
		Assert.assertTrue("Not less than Poisson information", I <= upper);
		Assert.assertTrue("Not greater than Poisson-Gaussian approximation information", I >= lower);
	}

	@Test
	public void canComputeFisherInformationWithLowestPossibleMean()
	{
		//org.junit.Assume.assumeTrue(false);

		// Lowest value where the reciprocal is not infinity.
		double u = Double.longBitsToDouble(0x4000000000001L);

		// Binary search for the min value
		final boolean doSearch = false;
		if (doSearch)
		{
			// This is the full 52-bit mantissa of a double with zero for the unbiased exponent,
			// i.e. the largest sub-normal number.
			long upper = (1L << 52) - 1;
			Assert.assertTrue(1.0 / Double.longBitsToDouble(upper) != Double.POSITIVE_INFINITY);
			long lower = 1;
			while (lower + 1 < upper)
			{
				// 1/Upper is not infinty
				// Test mid-point
				final long mid = (upper + lower) / 2;
				u = Double.longBitsToDouble(mid);
				if (1 / u == Double.POSITIVE_INFINITY)
					lower = mid;
				else
					// Mid point
					upper = mid;
			}

			u = Double.longBitsToDouble(upper);
			System.out.printf("upper = 0x%s = %s\n", Long.toHexString(upper), u);
		}

		Assert.assertTrue(1.0 / u != Double.POSITIVE_INFINITY);
		Assert.assertTrue(1.0 / Math.nextDown(u) == Double.POSITIVE_INFINITY);

		for (int i = 0; i < 4; i++)
		{
			final double s = (1 << i) * 0.25;
			final PoissonGaussianFisherInformation f = new PoissonGaussianFisherInformation(s);
			final double I = f.getPoissonGaussianI(u);
			final double I2 = f.getPoissonGaussianI(1e-100);
			final double lower = f.getPoissonGaussianApproximationI(u);
			final double upper = PoissonFisherInformation.getPoissonI(u);
			final double alpha = I / upper;
			TestLog.info("s=%g u=%g I=%s I(1e-100)=%s (%s - %s) alpha=%s\n", f.s, u, I, I2, lower, upper, alpha);
			Assert.assertTrue(I > lower);
			Assert.assertTrue(I < upper);

			// Test convergence
			Assert.assertEquals(I, I2, 1e-50);
		}
	}
}
