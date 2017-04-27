package gdsc.smlm.function.gaussian;

import gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import gdsc.smlm.function.gaussian.SingleEllipticalGaussian2DFunction;

public class SingleEllipticalGaussian2DFunctionTest extends Gaussian2DFunctionTest
{
	protected void init()
	{
		flags = GaussianFunctionFactory.FIT_ELLIPTICAL; 
		f1 = new SingleEllipticalGaussian2DFunction(maxx, maxx);
	}
}
