package gdsc.smlm.function.gaussian;

import gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import gdsc.smlm.function.gaussian.NBCircularGaussian2DFunction;

public class NBCircularGaussian2DFunctionTest extends Gaussian2DFunctionTest
{
	protected void init()
	{
		flags = GaussianFunctionFactory.FIT_NB_CIRCLE; 
		f1 = new NBCircularGaussian2DFunction(1, maxx, maxy);
		f2 = new NBCircularGaussian2DFunction(2, maxx, maxy);
	}
}
