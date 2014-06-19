package gdsc.smlm.fitting.function.gaussian;

import gdsc.smlm.fitting.function.GaussianFunctionFactory;
import gdsc.smlm.fitting.function.gaussian.FreeCircularGaussian2DFunction;

public class FreeCircularGaussian2DFunctionTest extends Gaussian2DFunctionTest
{
	protected void init()
	{
		flags = GaussianFunctionFactory.FIT_FREE_CIRCLE; 
		f1 = new FreeCircularGaussian2DFunction(1, maxx);
		f2 = new FreeCircularGaussian2DFunction(2, maxx);
	}
}
