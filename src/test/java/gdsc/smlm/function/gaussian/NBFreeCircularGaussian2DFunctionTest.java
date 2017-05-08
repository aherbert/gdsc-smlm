package gdsc.smlm.function.gaussian;

import gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import gdsc.smlm.function.gaussian.NBFreeCircularGaussian2DFunction;

public class NBFreeCircularGaussian2DFunctionTest extends Gaussian2DFunctionTest
{
	protected void init()
	{
		flags = GaussianFunctionFactory.FIT_NB_FREE_CIRCLE; 
		f1 = new NBFreeCircularGaussian2DFunction(1, maxx, maxy);
		f2 = new NBFreeCircularGaussian2DFunction(2, maxx, maxy);
	}
}
