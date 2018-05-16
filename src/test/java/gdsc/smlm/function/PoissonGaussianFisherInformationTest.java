package gdsc.smlm.function;

import org.junit.Assert;
import org.junit.Test;

public class PoissonGaussianFisherInformationTest
{
	@Test
	public void canComputeRealFisherInformation()
	{
		//org.junit.Assume.assumeTrue(false);

		//double u;
		////u = Math.pow(10, -300);
		////u = 1.0 / Math.nextAfter(Double.MAX_VALUE, -1); // Smallest p with a non-infinite Fisher information
		//u = Double.MIN_NORMAL;
		//PoissonGaussianFisherInformation f = new RealPoissonGaussianFisherInformation(0.25);
		//f.setMeanThreshold(Double.MAX_VALUE);
		//double I = f.getPoissonGaussianI(u);
		//double lower = f.getPoissonGaussianApproximationI(u);
		//double upper = PoissonFisherInformation.getPoissonI(u);
		//System.out.printf("s=%g u=%g I=%s  (%s - %s) alpha=%s\n", f.s, u, I, lower, upper, I / upper);
		//if (true)
		//	return;

		for (int i = 0; i < 4; i++)
		{
			double s = (1 << i) * 0.25;
			canComputeRealFisherInformation(s);
		}
	}

	private void canComputeRealFisherInformation(double s)
	{
		canComputeFisherInformation(new RealPoissonGaussianFisherInformation(s));
	}

	@Test(expected = AssertionError.class)
	public void cannotComputeDiscreteFisherInformation()
	{
		// TODO - The discrete version currently fails this test.
		// Do the full differentiation of the PMF that is being used 
		// to determine if it is computing a correct Fisher information. 

		//org.junit.Assume.assumeTrue(false);

		for (int i = 0; i < 4; i++)
		{
			double s = (1 << i) * 0.25;
			canComputeDiscreteFisherInformation(s);
		}
	}

	private void canComputeDiscreteFisherInformation(double s)
	{
		canComputeFisherInformation(new DiscretePoissonGaussianFisherInformation(s));
	}

	private void canComputeFisherInformation(PoissonGaussianFisherInformation f)
	{
		f.setMeanThreshold(Double.MAX_VALUE);
		// Do not evaluate at very high mean. The switch to the approximation will occur
		// and the approximation is good.
		for (int exp = -20; exp < 6; exp++)
		{
			canComputeFisherInformation(f, Math.pow(10, exp * 0.5));
		}
	}

	private void canComputeFisherInformation(PoissonGaussianFisherInformation f, double u)
	{
		double I = f.getPoissonGaussianI(u);
		double lower = f.getPoissonGaussianApproximationI(u);
		double upper = PoissonFisherInformation.getPoissonI(u);
		System.out.printf("s=%g u=%g I=%s  (%s - %s) alpha=%s\n", f.s, u, I, lower, upper, I / upper);
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
	public void canComputeRealFisherInformationWithLowestPossibleMean()
	{
		//org.junit.Assume.assumeTrue(false);

		// Lowest value where the reciprocal is not infinity.
		double u = Double.longBitsToDouble(0x4000000000001L);

		// Binary search for the min value
		boolean doSearch = false;
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
				long mid = (upper + lower) / 2;
				u = Double.longBitsToDouble(mid);
				if (1 / u == Double.POSITIVE_INFINITY)
				{
					lower = mid;
				}
				else
				{
					// Mid point
					upper = mid;
				}
			}

			u = Double.longBitsToDouble(upper);
			System.out.printf("upper = 0x%s = %s\n", Long.toHexString(upper), u);
		}

		Assert.assertTrue(1.0 / u != Double.POSITIVE_INFINITY);
		Assert.assertTrue(1.0 / Math.nextAfter(u, -1) == Double.POSITIVE_INFINITY);

		for (int i = 0; i < 4; i++)
		{
			double s = (1 << i) * 0.25;
			PoissonGaussianFisherInformation f = new RealPoissonGaussianFisherInformation(s);
			double I = f.getPoissonGaussianI(u);
			double lower = f.getPoissonGaussianApproximationI(u);
			double upper = PoissonFisherInformation.getPoissonI(u);
			double alpha = I / upper;
			System.out.printf("s=%g u=%g I=%s  (%s - %s) alpha=%s\n", f.s, u, I, lower, upper, alpha);
			Assert.assertTrue(I > lower);
			Assert.assertTrue(I < upper);
		}
	}
}