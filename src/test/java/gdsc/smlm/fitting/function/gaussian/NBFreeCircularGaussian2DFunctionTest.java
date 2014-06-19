package gdsc.smlm.fitting.function.gaussian;

import gdsc.smlm.fitting.function.GaussianFunctionFactory;
import gdsc.smlm.fitting.function.gaussian.NBFreeCircularGaussian2DFunction;

public class NBFreeCircularGaussian2DFunctionTest extends Gaussian2DFunctionTest
{
	protected void init()
	{
		flags = GaussianFunctionFactory.FIT_NB_FREE_CIRCLE; 
		f1 = new NBFreeCircularGaussian2DFunction(1, maxx);
		f2 = new NBFreeCircularGaussian2DFunction(2, maxx);
	}
}
