package gdsc.smlm.function;

import org.junit.Assert;
import org.junit.Test;

public class PoissonGaussianFisherInformationTest
{
	@Test
	public void canComputeRealFisherInformation()
	{
		//org.junit.Assume.assumeTrue(false);

		canComputeRealFisherInformation(0.25);
		canComputeRealFisherInformation(0.5);
		canComputeRealFisherInformation(1);
		canComputeRealFisherInformation(2);
	}

	private void canComputeRealFisherInformation(double s)
	{
		canComputeFisherInformation(new RealPoissonGaussianFisherInformation(s));
	}

	@Test(expected=AssertionError.class)
	public void cannotComputeDiscreteFisherInformation()
	{
		// TODO - The discrete version currently fails this test.
		// Do the full differentiation of the PMF that is being used 
		// to determine if it is computing a correct Fisher information. 
		
		//org.junit.Assume.assumeTrue(false);
		
		canComputeDiscreteFisherInformation(0.25);
		canComputeDiscreteFisherInformation(0.5);
		canComputeDiscreteFisherInformation(1);
		canComputeDiscreteFisherInformation(2);
	}

	private void canComputeDiscreteFisherInformation(double s)
	{
		canComputeFisherInformation(new DiscretePoissonGaussianFisherInformation(s));
	}
	
	private void canComputeFisherInformation(PoissonGaussianFisherInformation f)
	{
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
		double upper = PoissonFisherInformation.getPoissonI(u);
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
