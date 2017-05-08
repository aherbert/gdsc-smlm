package gdsc.smlm.function.gaussian;

public class NSFixedGaussian2DFunctionTest extends Gaussian2DFunctionTest
{
	protected void init()
	{
		flags = GaussianFunctionFactory.FIT_SIMPLE_NS_FIXED;
		f1 = new NSFixedGaussian2DFunction(1, maxx, maxy);
		f2 = new NSFixedGaussian2DFunction(2, maxx, maxy);
	}
}
