package gdsc.smlm.function.gaussian;

import gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import gdsc.smlm.function.gaussian.SingleCircularGaussian2DFunction;

public class SingleCircularGaussian2DFunctionTest extends Gaussian2DFunctionTest
{
	protected void init()
	{
		flags = GaussianFunctionFactory.FIT_CIRCLE; 
		f1 = new SingleCircularGaussian2DFunction(maxx, maxx);
	}
}
