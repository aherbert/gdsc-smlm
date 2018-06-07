package gdsc.smlm.function.gaussian.erf;

import gdsc.smlm.function.gaussian.GaussianFunctionFactory;

public class SingleNBFreeCircularErfGaussian2DFunctionTest extends SingleErfGaussian2DFunctionTest
{
	protected void init()
	{
		flags = GaussianFunctionFactory.FIT_ERF_NB_FREE_CIRCLE; 
		f1 = new SingleNBFreeCircularErfGaussian2DFunction(maxx, maxy);
	}
}
