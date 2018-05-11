package gdsc.smlm.function;

import org.junit.Assert;
import org.junit.Test;

public class PoissonGammaGaussianFisherInformationTest
{
	@Test
	public void canComputeRealFisherInformation()
	{
		//org.junit.Assume.assumeTrue(false);

		double[] M = { 20 };
		double[] S = { 3 };

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
		for (int exp = -12; exp < 8; exp++)
		{
			canComputeFisherInformation(f, Math.pow(10, exp * 0.5));
		}
	}

	private void canComputeFisherInformation(PoissonGammaGaussianFisherInformation f, double u)
	{
		double I = f.getPoissonGammaGaussianI(u);
		double upper = PoissonFisherInformation.getPoissonI(u);
		System.out.printf("s=%g u=%g I=%s PoissonI=%s alpha=%s\n", f.s, u, I, upper, I / upper);
		// Allow a tolerance on the approximation at high mean.
		// The function does not compute the sum to infinity and so can underestimate
		// the value.
		if (u >= 100)
			upper *= 1.001;
		//Assert.assertTrue("Not less than Poisson information", I <= upper);
	}
}
