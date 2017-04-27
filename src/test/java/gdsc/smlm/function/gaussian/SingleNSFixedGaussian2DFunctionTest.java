package gdsc.smlm.function.gaussian;

public class SingleNSFixedGaussian2DFunctionTest extends Gaussian2DFunctionTest
{
	protected void init()
	{
		flags = GaussianFunctionFactory.FIT_NS_FIXED;
		f1 = new SingleNSFixedGaussian2DFunction(maxx, maxx);
	}
}
