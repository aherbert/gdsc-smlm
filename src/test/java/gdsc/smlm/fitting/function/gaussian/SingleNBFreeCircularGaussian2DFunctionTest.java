package gdsc.smlm.fitting.function.gaussian;

import gdsc.smlm.fitting.function.GaussianFunctionFactory;
import gdsc.smlm.fitting.function.gaussian.SingleNBFreeCircularGaussian2DFunction;

public class SingleNBFreeCircularGaussian2DFunctionTest extends Gaussian2DFunctionTest
{
	protected void init()
	{
		flags = GaussianFunctionFactory.FIT_NB_FREE_CIRCLE; 
		f1 = new SingleNBFreeCircularGaussian2DFunction(maxx);
	}
}
