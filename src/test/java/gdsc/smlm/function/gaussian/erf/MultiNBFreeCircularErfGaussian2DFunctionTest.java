package gdsc.smlm.function.gaussian.erf;

import gdsc.smlm.function.gaussian.GaussianFunctionFactory;

public class MultiNBFreeCircularErfGaussian2DFunctionTest extends ErfGaussian2DFunctionTest
{
	protected void init()
	{
		flags = GaussianFunctionFactory.FIT_ERF_NB_FREE_CIRCLE; 
		f1 = new MultiNBFreeCircularErfGaussian2DFunction(1, maxx, maxy);
		f2 = new MultiNBFreeCircularErfGaussian2DFunction(2, maxx, maxy);
	}
}
