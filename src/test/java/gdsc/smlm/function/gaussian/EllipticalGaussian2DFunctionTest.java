package gdsc.smlm.function.gaussian;

import gdsc.smlm.function.gaussian.EllipticalGaussian2DFunction;
import gdsc.smlm.function.gaussian.GaussianFunctionFactory;

public class EllipticalGaussian2DFunctionTest extends Gaussian2DFunctionTest
{
	protected void init()
	{
		flags = GaussianFunctionFactory.FIT_SIMPLE_ELLIPTICAL; 
		f1 = new EllipticalGaussian2DFunction(1, maxx, maxy);
		f2 = new EllipticalGaussian2DFunction(2, maxx, maxy);
	}
}
