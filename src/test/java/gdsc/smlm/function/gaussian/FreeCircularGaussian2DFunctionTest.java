package gdsc.smlm.function.gaussian;

import gdsc.smlm.function.gaussian.FreeCircularGaussian2DFunction;
import gdsc.smlm.function.gaussian.GaussianFunctionFactory;

public class FreeCircularGaussian2DFunctionTest extends Gaussian2DFunctionTest
{
	protected void init()
	{
		flags = GaussianFunctionFactory.FIT_FREE_CIRCLE; 
		f1 = new FreeCircularGaussian2DFunction(1, maxx, maxx);
		f2 = new FreeCircularGaussian2DFunction(2, maxx, maxx);
	}
}
