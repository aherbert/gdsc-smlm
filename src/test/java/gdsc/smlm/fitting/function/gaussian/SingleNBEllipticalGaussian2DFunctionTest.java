package gdsc.smlm.fitting.function.gaussian;

import gdsc.smlm.fitting.function.GaussianFunctionFactory;
import gdsc.smlm.fitting.function.gaussian.SingleNBEllipticalGaussian2DFunction;

public class SingleNBEllipticalGaussian2DFunctionTest extends Gaussian2DFunctionTest
{
	protected void init()
	{
		flags = GaussianFunctionFactory.FIT_NB_ELLIPTICAL; 
		f1 = new SingleNBEllipticalGaussian2DFunction(maxx);
	}
}
