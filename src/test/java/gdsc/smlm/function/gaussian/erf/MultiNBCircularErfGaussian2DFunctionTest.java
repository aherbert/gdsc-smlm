package gdsc.smlm.function.gaussian.erf;

import gdsc.smlm.function.gaussian.GaussianFunctionFactory;

public class MultiNBCircularErfGaussian2DFunctionTest extends ErfGaussian2DFunctionTest
{
	protected void init()
	{
		flags = GaussianFunctionFactory.FIT_ERF_NB_CIRCLE; 
		f1 = new MultiNBCircularErfGaussian2DFunction(1, maxx, maxy);
		f2 = new MultiNBCircularErfGaussian2DFunction(2, maxx, maxy);
	}
}
