package gdsc.smlm.fitting.function.gaussian;

import gdsc.smlm.fitting.function.GaussianFunctionFactory;
import gdsc.smlm.fitting.function.gaussian.NBEllipticalGaussian2DFunction;

public class NBEllipticalGaussian2DFunctionTest extends Gaussian2DFunctionTest
{
	protected void init()
	{
		flags = GaussianFunctionFactory.FIT_NB_ELLIPTICAL; 
		f1 = new NBEllipticalGaussian2DFunction(1, maxx);
		f2 = new NBEllipticalGaussian2DFunction(2, maxx);
	}
}
