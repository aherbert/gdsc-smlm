package gdsc.smlm.fitting.function.gaussian;

import gdsc.smlm.fitting.function.GaussianFunctionFactory;
import gdsc.smlm.fitting.function.gaussian.FixedGaussian2DFunction;

public class FixedGaussian2DFunctionTest extends Gaussian2DFunctionTest
{
	protected void init()
	{
		flags = GaussianFunctionFactory.FIT_FIXED; 
		f1 = new FixedGaussian2DFunction(1, maxx);
		f2 = new FixedGaussian2DFunction(2, maxx);
	}
}
