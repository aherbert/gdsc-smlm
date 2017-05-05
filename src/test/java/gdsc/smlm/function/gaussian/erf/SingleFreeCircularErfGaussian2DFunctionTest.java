package gdsc.smlm.function.gaussian.erf;

import gdsc.smlm.function.gaussian.GaussianFunctionFactory;

public class SingleFreeCircularErfGaussian2DFunctionTest extends ErfGaussian2DFunctionTest
{
	protected void init()
	{
		flags = GaussianFunctionFactory.FIT_ERF_FREE_CIRCLE; 
		f1 = new SingleFreeCircularErfGaussian2DFunction(maxx, maxx);
	}
}
