package gdsc.smlm.function.gaussian;

import gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import gdsc.smlm.function.gaussian.SingleNBEllipticalGaussian2DFunction;

public class SingleNBEllipticalGaussian2DFunctionTest extends Gaussian2DFunctionTest
{
	protected void init()
	{
		flags = GaussianFunctionFactory.FIT_SIMPLE_NB_ELLIPTICAL; 
		f1 = new SingleNBEllipticalGaussian2DFunction(maxx, maxx);
	}
}
