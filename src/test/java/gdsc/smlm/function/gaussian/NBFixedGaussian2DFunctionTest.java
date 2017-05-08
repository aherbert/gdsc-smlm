package gdsc.smlm.function.gaussian;

import gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import gdsc.smlm.function.gaussian.NBFixedGaussian2DFunction;

public class NBFixedGaussian2DFunctionTest extends Gaussian2DFunctionTest
{
	protected void init()
	{
		flags = GaussianFunctionFactory.FIT_SIMPLE_NB_FIXED; 
		f1 = new NBFixedGaussian2DFunction(1, maxx, maxy);
		f2 = new NBFixedGaussian2DFunction(2, maxx, maxy);
	}
}
