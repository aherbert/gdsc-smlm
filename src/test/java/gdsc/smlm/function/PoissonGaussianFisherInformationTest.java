package gdsc.smlm.function;

import org.junit.Assert;
import org.junit.Test;

public class PoissonGaussianFisherInformationTest
{
	@Test
	public void canComputeFisherInformation()
	{
		PoissonGaussianFisherInformation f = new PoissonGaussianFisherInformation(2, 5);

		for (double u = 1; u < 100; u++)
		{
			double I = f.getFisherInformation(u);
			//System.out.printf("u=%g I=%s\n", u, I);
		}
	}
}
