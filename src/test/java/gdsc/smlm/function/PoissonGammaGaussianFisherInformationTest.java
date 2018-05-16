package gdsc.smlm.function;

import org.apache.commons.math3.util.FastMath;
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
		//for (int exp = -590, i=0; i<2; exp++, i++)
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
	public void canComputeLowerLimit()
	{
		// The relative Fisher information plateaus at low photons. Output where this is.
		for (double p : new double[] { 1e-9, 1e-10, 1e-11 })
		{
			double exp_p = FastMath.exp(-p);
			double upper = PoissonFisherInformation.getPoissonI(p);
			double lastAlpha = 1;
			double lastChange = 1;
			double[] M = { 20, 100, 500 };
			double[] S = { 3, 13 };

			for (double m : M)
				for (double s : S)
				{
					PoissonGammaGaussianFisherInformation f = new RealPoissonGammaGaussianFisherInformation(m, s);
					f.setCumulativeProbability(1 - 1e-12);
					double I = f.getPoissonGammaGaussianI(p);
					double alpha = I / upper;
					double change = lastAlpha / alpha;
					System.out.printf("p=%g  exp(-p)=%s   s=%g   I=%s  alpha=%s  (%s)  %s\n", p, exp_p, s, I, alpha,
							change, lastChange / change);
					lastAlpha = alpha;
					lastChange = change;
				}
		}
	}

	@Test(expected = AssertionError.class)
	public void cannotComputeRealFisherInformationWithLowestPossibleMean()
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

		double[] M = { 20, 100, 500 };
		double[] S = { 3, 13 };

		for (double m : M)
			for (double s : S)
			{
				PoissonGammaGaussianFisherInformation f = new RealPoissonGammaGaussianFisherInformation(m, s);
				double I = f.getPoissonGammaGaussianI(u);
				double upper = PoissonFisherInformation.getPoissonI(u);
				System.out.printf("m=%g s=%g u=%g I=%s PoissonI=%s alpha=%s\n", f.m, f.s, u, I, upper, I / upper);
				//Assert.assertTrue(I < upper);
				//Assert.assertTrue(alpha > 0);
			}
	}
}
