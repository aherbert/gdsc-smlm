package gdsc.smlm.function.gaussian;

public class SingleNSNBFixedGaussian2DFunctionTest extends Gaussian2DFunctionTest
{
	protected void init()
	{
		flags = GaussianFunctionFactory.FIT_NS_NB_FIXED;
		f1 = new SingleNSNBFixedGaussian2DFunction(maxx);
	}
}
