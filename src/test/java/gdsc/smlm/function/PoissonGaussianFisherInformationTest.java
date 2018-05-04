package gdsc.smlm.function;

import org.junit.Assert;
import org.junit.Test;

public class PoissonGaussianFisherInformationTest
{
	@Test
	public void canComputeFisherInformation()
	{
		//org.junit.Assume.assumeTrue(false);

		canComputeFisherInformation(0.25);
		canComputeFisherInformation(0.5);
		canComputeFisherInformation(1);
		canComputeFisherInformation(2);
	}

	private void canComputeFisherInformation(double s)
	{
		PoissonGaussianFisherInformation f = new PoissonGaussianFisherInformation(s);
		f.setMeanThreshold(Double.MAX_VALUE);
		for (int exp = -12; exp < 8; exp++)
		{
			canComputeFisherInformation(f, Math.pow(10, exp * 0.5));
		}
	}

	private void canComputeFisherInformation(PoissonGaussianFisherInformation f, double u)
	{
		double I = f.getPoissonGaussianI(u);
		//System.out.printf("s=%g u=%g I=%s\n", f.s, u, I);
		double lower = f.getPoissonGaussianApproximationI(u);
		double upper = f.getPoissonI(u);
		// Allow a tolerance on the approximation at high mean.
		// The function does not compute the sum to infinity and so can underestimate
		// the value.
		if (u >= 10)
			lower *= 0.99;
		if (u >= 100)
			upper *= 1.001;
		Assert.assertTrue("Not less than Poisson information", I <= upper);
		Assert.assertTrue("Not greater than Poisson-Gaussian approximation information", I >= lower);
	}
}
