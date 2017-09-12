package gdsc.smlm.function.gaussian.erf;

import gdsc.smlm.function.gaussian.GaussianFunctionFactory;

public class SingleNBCircularErfGaussian2DFunctionTest extends ErfGaussian2DFunctionTest
{
	protected void init()
	{
		flags = GaussianFunctionFactory.FIT_ERF_NB_CIRCLE; 
		f1 = new SingleNBCircularErfGaussian2DFunction(maxx, maxy);
	}
}
