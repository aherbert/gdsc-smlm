package gdsc.smlm.function.gaussian.erf;

import gdsc.smlm.function.gaussian.GaussianFunctionFactory;

public class MultiFreeCircularErfGaussian2DFunctionTest extends ErfGaussian2DFunctionTest
{
	protected void init()
	{
		flags = GaussianFunctionFactory.FIT_ERF_FREE_CIRCLE; 
		f1 = new MultiFreeCircularErfGaussian2DFunction(1, maxx, maxy);
		f2 = new MultiFreeCircularErfGaussian2DFunction(2, maxx, maxy);
	}
}
