package gdsc.smlm.function.gaussian;

import gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import gdsc.smlm.function.gaussian.SingleNBFixedGaussian2DFunction;

public class SingleNBFixedGaussian2DFunctionTest extends Gaussian2DFunctionTest
{
	protected void init()
	{
		flags = GaussianFunctionFactory.FIT_SIMPLE_NB_FIXED; 
		f1 = new SingleNBFixedGaussian2DFunction(maxx, maxx);
	}
}
