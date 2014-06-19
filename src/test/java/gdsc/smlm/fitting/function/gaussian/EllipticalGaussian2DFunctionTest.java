package gdsc.smlm.fitting.function.gaussian;

import gdsc.smlm.fitting.function.GaussianFunctionFactory;
import gdsc.smlm.fitting.function.gaussian.EllipticalGaussian2DFunction;

public class EllipticalGaussian2DFunctionTest extends Gaussian2DFunctionTest
{
	protected void init()
	{
		flags = GaussianFunctionFactory.FIT_ELLIPTICAL; 
		f1 = new EllipticalGaussian2DFunction(1, maxx);
		f2 = new EllipticalGaussian2DFunction(2, maxx);
	}
}
