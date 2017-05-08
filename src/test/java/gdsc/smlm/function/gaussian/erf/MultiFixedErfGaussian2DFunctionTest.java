package gdsc.smlm.function.gaussian.erf;

import gdsc.smlm.function.gaussian.GaussianFunctionFactory;

public class MultiFixedErfGaussian2DFunctionTest extends ErfGaussian2DFunctionTest
{
	protected void init()
	{
		flags = GaussianFunctionFactory.FIT_ERF_FIXED;
		f1 = new MultiFixedErfGaussian2DFunction(1, maxx, maxy);
		f2 = new MultiFixedErfGaussian2DFunction(2, maxx, maxy);
	}
}
