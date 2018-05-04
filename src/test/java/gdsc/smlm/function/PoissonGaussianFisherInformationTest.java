package gdsc.smlm.function;

import org.junit.Assume;
import org.junit.Test;

public class PoissonGaussianFisherInformationTest
{
	@SuppressWarnings("unused")
	@Test
	public void canComputeFisherInformation()
	{
		//Assume.assumeTrue(false);
		
		PoissonGaussianFisherInformation f = new PoissonGaussianFisherInformation(0.5, 5);
		f.setMeanThreshold(Double.MAX_VALUE);

		for (int exp = -12; exp < 8; exp++)
		{
			canComputeFisherInformation(f, Math.pow(10, exp*0.5));
		}
		
		if (true)
			return;
		
		for (int exp = -6; exp < 0; exp++)
		{
			canComputeFisherInformation(f, Math.pow(10, exp));
		}

		for (double u = 1; u <= 100; u++)
		{
			canComputeFisherInformation(f, u);
		}

		// Check the switch over to 2-sided cumulative is 'smooth'
		canComputeFisherInformation(f, 101);
		canComputeFisherInformation(f, 102);
		
		for (int exp = 3; exp < 6; exp++)
		{
			canComputeFisherInformation(f, Math.pow(10, exp));
		}
	}
	
	private void canComputeFisherInformation(PoissonGaussianFisherInformation f, double u)
	{
		double I = f.getFisherInformation(u);
		//System.out.printf("u=%g I=%s\n", u, I);
	}
}
