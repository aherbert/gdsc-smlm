package gdsc.smlm.function;

import org.junit.Assert;
import org.junit.Test;

public class PoissonGaussianFisherInformationTest
{
	@Test
	public void canComputeFisherInformation()
	{
		PoissonGaussianFisherInformation f = new PoissonGaussianFisherInformation(2, 5);

		// TODO - let the function store the last computed convolution.
		// Create a plugin to compute the alpha parameter across a range of means.
		// Draw the P-G convolution and the computed function.
		// Is it smooth enough for a Simpson sum?
		// Maybe the simpson sum is not working at low mean as the convolution/function
		// is not smooth.
		
		for (double u = 1; u < 100; u++)
		{
			double I = f.getFisherInformation(u);
			//System.out.printf("u=%g I=%s\n", u, I);
		}
	}
}
