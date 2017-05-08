package gdsc.smlm.function.gaussian;

import gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import gdsc.smlm.function.gaussian.SingleNBFreeCircularGaussian2DFunction;

public class SingleNBFreeCircularGaussian2DFunctionTest extends Gaussian2DFunctionTest
{
	protected void init()
	{
		flags = GaussianFunctionFactory.FIT_SIMPLE_NB_FREE_CIRCLE; 
		f1 = new SingleNBFreeCircularGaussian2DFunction(maxx, maxx);
	}
}
