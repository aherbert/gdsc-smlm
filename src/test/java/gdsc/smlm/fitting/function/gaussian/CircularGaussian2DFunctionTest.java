package gdsc.smlm.fitting.function.gaussian;

import gdsc.smlm.fitting.function.GaussianFunctionFactory;
import gdsc.smlm.fitting.function.gaussian.CircularGaussian2DFunction;

public class CircularGaussian2DFunctionTest extends Gaussian2DFunctionTest
{
	protected void init()
	{
		flags = GaussianFunctionFactory.FIT_CIRCLE; 
		f1 = new CircularGaussian2DFunction(1, maxx);
		f2 = new CircularGaussian2DFunction(2, maxx);
	}
}
