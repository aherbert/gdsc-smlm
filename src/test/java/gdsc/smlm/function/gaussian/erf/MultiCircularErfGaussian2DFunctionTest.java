package gdsc.smlm.function.gaussian.erf;

import gdsc.smlm.function.gaussian.GaussianFunctionFactory;

public class MultiCircularErfGaussian2DFunctionTest extends ErfGaussian2DFunctionTest
{
	protected void init()
	{
		flags = GaussianFunctionFactory.FIT_ERF_CIRCLE; 
		f1 = new MultiCircularErfGaussian2DFunction(1, maxx, maxy);
		f2 = new MultiCircularErfGaussian2DFunction(2, maxx, maxy);
	}
}
