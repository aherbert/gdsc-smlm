package gdsc.smlm.function.gaussian.erf;

import gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import gdsc.smlm.function.gaussian.HoltzerAstimatismZModel;

public class SingleAstigmatismErfGaussian2DFunctionTest extends ErfGaussian2DFunctionTest
{
	protected void init()
	{
		flags = GaussianFunctionFactory.FIT_ERF_ASTIGMATISM;
		// Use a reasonable z-depth function from the Smith, et al (2010) paper (page 377)
		double gamma = 0.389;
		double d = 0.531;
		double Ax = -0.0708;
		double Bx = -0.073;
		double Ay = 0.164;
		double By = 0.0417;
		zModel = HoltzerAstimatismZModel.create(gamma, d, Ax, Bx, Ay, By);
		f1 = new SingleAstigmatismErfGaussian2DFunction(maxx, maxx, 2, zModel);
	}

	protected void postInit()
	{
		// Even though the function does not evaluate the widths it can use them
		// to construct independent widths.
		// Test with different X and Y SD 		
		testw1 = new double[][] { { 1.1, 1.1 }, { 1.1, 1.2 }, { 1.1, 1.5 } };
	};
}
