package gdsc.smlm.fitting.function.gaussian;

import gdsc.smlm.fitting.function.GaussianFunctionFactory;
import gdsc.smlm.fitting.function.gaussian.SingleFreeCircularGaussian2DFunction;

public class SingleFreeCircularGaussian2DFunctionTest extends Gaussian2DFunctionTest
{
	protected void init()
	{
		flags = GaussianFunctionFactory.FIT_FREE_CIRCLE; 
		f1 = new SingleFreeCircularGaussian2DFunction(maxx);
	}
}
