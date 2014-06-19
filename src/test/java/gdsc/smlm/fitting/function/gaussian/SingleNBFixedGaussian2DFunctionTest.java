package gdsc.smlm.fitting.function.gaussian;

import gdsc.smlm.fitting.function.GaussianFunctionFactory;
import gdsc.smlm.fitting.function.gaussian.SingleNBFixedGaussian2DFunction;

public class SingleNBFixedGaussian2DFunctionTest extends Gaussian2DFunctionTest
{
	protected void init()
	{
		flags = GaussianFunctionFactory.FIT_NB_FIXED; 
		f1 = new SingleNBFixedGaussian2DFunction(maxx);
	}
}
