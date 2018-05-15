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
}
