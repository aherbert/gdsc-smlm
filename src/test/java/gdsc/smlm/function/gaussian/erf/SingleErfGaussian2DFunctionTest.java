package gdsc.smlm.function.gaussian.erf;

import org.junit.Assert;
import org.junit.Test;

import gdsc.core.utils.Maths;

public abstract class SingleErfGaussian2DFunctionTest extends ErfGaussian2DFunctionTest
{
	@Test
	public void canComputeIntegral()
	{
		SingleErfGaussian2DFunction ef1 = (SingleErfGaussian2DFunction) f1;

		double[] a;
		for (double background : testbackground)
			// Peak 1
			for (double signal1 : testsignal1)
				for (double cx1 : testcx1)
					for (double cy1 : testcy1)
						for (double cz1 : testcz1)
							for (double[] w1 : testw1)
								for (double angle1 : testangle1)
								{
									a = createParameters(background, signal1, cx1, cy1, cz1, w1[0], w1[1], angle1);
									double e = Maths.sum(f1.computeValues(a));
									double o = ef1.integral(a);
									Assert.assertEquals(e, o, e * 1e-8);
								}
	}
}
