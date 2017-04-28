package gdsc.smlm.function.gaussian.erf;

import gdsc.smlm.function.gaussian.GaussianFunctionFactory;

public class SingleFixedErfGaussian2DFunctionTest extends ErfGaussian2DFunctionTest
{
	protected void init()
	{
		flags = GaussianFunctionFactory.FIT_ERF_FIXED; 
		f1 = new SingleFixedErfGaussian2DFunction(maxx, maxx, 2);
	}
}
