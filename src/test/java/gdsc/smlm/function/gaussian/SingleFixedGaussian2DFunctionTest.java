package gdsc.smlm.function.gaussian;

import gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import gdsc.smlm.function.gaussian.SingleFixedGaussian2DFunction;

public class SingleFixedGaussian2DFunctionTest extends Gaussian2DFunctionTest
{
	protected void init()
	{
		flags = GaussianFunctionFactory.FIT_FIXED; 
		f1 = new SingleFixedGaussian2DFunction(maxx, maxx);
	}
}
