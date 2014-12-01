package gdsc.smlm.function.gaussian;

import gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import gdsc.smlm.function.gaussian.SingleNBCircularGaussian2DFunction;

public class SingleNBCircularGaussian2DFunctionTest extends Gaussian2DFunctionTest
{
	protected void init()
	{
		flags = GaussianFunctionFactory.FIT_NB_CIRCLE; 
		f1 = new SingleNBCircularGaussian2DFunction(maxx);
	}
}
