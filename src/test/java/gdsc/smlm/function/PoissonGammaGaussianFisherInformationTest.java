package gdsc.smlm.function;

import org.junit.Assert;
import org.junit.Test;

public class PoissonGammaGaussianFisherInformationTest
{
	@Test
	public void canComputeRealFisherInformation()
	{
		//org.junit.Assume.assumeTrue(false);

		//canComputeRealFisherInformation(250, 13);
		//double m1 = 500;
		//double u = 350;
		//double s1 = 13;
		//PoissonGammaGaussianFisherInformation f = new RealPoissonGammaGaussianFisherInformation(m1, s1);
		//f.setMeanThreshold(Double.POSITIVE_INFINITY);
		//double I = f.getPoissonGammaGaussianI(u);
		//double upper = PoissonFisherInformation.getPoissonI(u);
		//System.out.printf("s=%g u=%g I=%s PoissonI=%s alpha=%s\n", s1, u, I, upper, I / upper);
		//if (true)
		//	return;

		double[] M = { 20, 100, 500 };
		double[] S = { 3, 13 };

		for (double m : M)
			for (double s : S)
				canComputeRealFisherInformation(m, s);
	}

	private void canComputeRealFisherInformation(double m, double s)
	{
		canComputeFisherInformation(new RealPoissonGammaGaussianFisherInformation(m, s));
	}

	private void canComputeFisherInformation(PoissonGammaGaussianFisherInformation f)
	{
		f.setMeanThreshold(Double.POSITIVE_INFINITY);

		// Due to a limited convolution (for an infinite range) 
		// the class works up to mean of about 300. Above that the approximation using
		// half the Poisson Fisher information should be used instead. 
		// 10^2 == 100 (OK), 10^2.5 == 316 (Fail)
		for (int exp = -12; exp <= 4; exp++)
		{
			canComputeFisherInformation(f, Math.pow(10, exp * 0.5));
		}
	}

	private void canComputeFisherInformation(PoissonGammaGaussianFisherInformation f, double u)
	{
		double I = f.getPoissonGammaGaussianI(u);
		double upper = PoissonFisherInformation.getPoissonI(u);
		System.out.printf("m=%g s=%g u=%g I=%s PoissonI=%s alpha=%s\n", f.m, f.s, u, I, upper, I / upper);
		Assert.assertTrue("Not less than Poisson information", I <= upper);
		// This is true at higher mean
		if (u > 1)
			Assert.assertTrue("Not above half the Poisson information", I >= 0.4999 * upper);
	}

	@Test
	public void canComputeAlpha()
	{
		org.junit.Assume.assumeTrue(false);

		// Compute the alpha using a range of gain and standard deviation

		int minm = 100;
		int maxm = 200;

		// When Poisson mean is low s does not matter as the Dirac is insignificant
		// and the convolution is mute. This may not be true when m is low but at higher
		// m the function is smooth and convolution has little effect.
		int mins = 5;
		int maxs = 5;

		// The relative Fisher information plateaus at low photons. 
		// Output where this is.
		double p = 1e-100;
		double p2 = 1e-200; //Double.MIN_NORMAL; //Double.longBitsToDouble(0x4000000000001L);

		double upper = PoissonFisherInformation.getPoissonI(p);
		double upper2 = PoissonFisherInformation.getPoissonI(p2);
		double lastAlpha = 1;

		for (int s = mins; s <= maxs; s++)
			for (int m = minm; m <= maxm; m++)
			{
				PoissonGammaGaussianFisherInformation f = new RealPoissonGammaGaussianFisherInformation(m, s);
				f.setCumulativeProbability(1 - 1e-12);
				double I = f.getPoissonGammaGaussianI(p);
				double I2 = f.getPoissonGammaGaussianI(p2);
				double alpha = I / upper;
				double alpha2 = I2 / upper2;
				double change = lastAlpha / alpha;
				System.out.printf("p=%g  p2=%s   m=%s  s=%s   I=%s (%s)  alpha=%s (%s)  (delta=%s)\n", p, p2, m, s, I,
						I2, alpha, alpha2, change);
				lastAlpha = alpha;
			}
	}

	@Test(expected = AssertionError.class)
	public void cannotComputeRealFisherInformationWithLowestPossibleMean()
	{
		org.junit.Assume.assumeTrue(false);

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

		computeRealFisherInformationWithMean(u);
	}

	@Test
	public void canComputeRealFisherInformationWithLowMean()
	{
		org.junit.Assume.assumeTrue(false);
		double u = 1e-300;
		computeRealFisherInformationWithMean(u);
	}

	private void computeRealFisherInformationWithMean(double u)
	{
		double[] M = { 20, 100, 500 };
		double[] S = { 3, 13 };

		for (double m : M)
			for (double s : S)
			{
				PoissonGammaGaussianFisherInformation f = new RealPoissonGammaGaussianFisherInformation(m, s);
				double I = f.getPoissonGammaGaussianI(u);
				double upper = PoissonFisherInformation.getPoissonI(u);
				double alpha = I / upper;
				System.out.printf("m=%g s=%g u=%g I=%s PoissonI=%s alpha=%s\n", f.m, f.s, u, I, upper, alpha);
				Assert.assertTrue(I < upper);
				Assert.assertTrue(alpha > 0);
			}
	}

	@Test
	public void canApproximateRealFisherInformationWithLowMean()
	{
		double[] M = { 20, 100, 500 };
		double[] S = { 3, 13 };

		for (double m : M)
			for (double s : S)
				canApproximateRealFisherInformationWithLowMean(m, s);
	}

	private void canApproximateRealFisherInformationWithLowMean(double m, double s)
	{
		canApproximateRealFisherInformationWithLowMean(new RealPoissonGammaGaussianFisherInformation(m, s));
	}

	private void canApproximateRealFisherInformationWithLowMean(PoissonGammaGaussianFisherInformation f)
	{
		PoissonGammaGaussianFisherInformation fe = (PoissonGammaGaussianFisherInformation) f.clone();

		f.setLowerMeanThreshold(1e-20);
		fe.setLowerMeanThreshold(0);

		canApproximateRealFisherInformationWithLowMean(fe, f, 1e-21);
		canApproximateRealFisherInformationWithLowMean(fe, f, 1e-100);
		canApproximateRealFisherInformationWithLowMean(fe, f, 1e-200);
	}

	private void canApproximateRealFisherInformationWithLowMean(PoissonGammaGaussianFisherInformation fe,
			PoissonGammaGaussianFisherInformation fo, double u)
	{
		double e = fe.getPoissonGammaGaussianI(u);
		double o = fo.getPoissonGammaGaussianI(u);

		//System.out.printf("m=%g s=%g u=%g e=%s o=%s\n", fe.m, fe.s, u, e, o);
		Assert.assertEquals(e, o, e * 1e-2);
	}
}
