package gdsc.smlm.function.gaussian;

import gdsc.smlm.function.gaussian.FixedGaussian2DFunction;
import gdsc.smlm.function.gaussian.GaussianFunctionFactory;

public class FixedGaussian2DFunctionTest extends Gaussian2DFunctionTest
{
	protected void init()
	{
		flags = GaussianFunctionFactory.FIT_FIXED; 
		f1 = new FixedGaussian2DFunction(1, maxx, maxy);
		f2 = new FixedGaussian2DFunction(2, maxx, maxy);
	}
}
