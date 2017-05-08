package gdsc.smlm.function.gaussian;

import gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import gdsc.smlm.function.gaussian.NBEllipticalGaussian2DFunction;

public class NBEllipticalGaussian2DFunctionTest extends Gaussian2DFunctionTest
{
	protected void init()
	{
		flags = GaussianFunctionFactory.FIT_NB_ELLIPTICAL; 
		f1 = new NBEllipticalGaussian2DFunction(1, maxx, maxy);
		f2 = new NBEllipticalGaussian2DFunction(2, maxx, maxy);
	}
}
