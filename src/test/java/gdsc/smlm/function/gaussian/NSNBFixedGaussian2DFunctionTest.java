package gdsc.smlm.function.gaussian;

public class NSNBFixedGaussian2DFunctionTest extends Gaussian2DFunctionTest
{
	protected void init()
	{
		flags = GaussianFunctionFactory.FIT_NS_NB_FIXED;
		f1 = new NSNBFixedGaussian2DFunction(1, maxx, maxx);
		f2 = new NSNBFixedGaussian2DFunction(2, maxx, maxx);
	}
}
