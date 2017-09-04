package gdsc.smlm.function.gaussian.erf;

import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.utils.DoubleEquality;
import gdsc.smlm.function.StandardValueProcedure;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import gdsc.smlm.model.GaussianPSFModel;

public class ErfGaussian2DFunctionVsPSFModelTest
{
	private int width = 10;
	private int height = 9;

	@Test
	public void computesSameAsPSFModel()
	{
		RandomDataGenerator r = new RandomDataGenerator(new Well19937c(30051977));
		for (int i = 0; i < 10; i++)
		{
			//@formatter:off
			computesSameAsPSFModel(
					r.nextUniform(50, 100),
					r.nextUniform((width-1)/2.0, (width+1)/2.0),
					r.nextUniform((height-1)/2.0, (height+1)/2.0),
					r.nextUniform(0.5, 2),
					r.nextUniform(0.5, 2));
			//@formatter:on
		}
	}

	private void computesSameAsPSFModel(double sum, double x0, double x1, double s0, double s1)
	{
		Gaussian2DFunction f = GaussianFunctionFactory.create2D(1, width, height,
				GaussianFunctionFactory.FIT_ERF_FREE_CIRCLE, null);
		double[] a = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
		a[Gaussian2DFunction.SIGNAL] = sum;
		a[Gaussian2DFunction.X_POSITION] = x0;
		a[Gaussian2DFunction.Y_POSITION] = x1;
		a[Gaussian2DFunction.X_SD] = s0;
		a[Gaussian2DFunction.Y_SD] = s1;
		double[] o = new StandardValueProcedure().getValues(f, a);

		GaussianPSFModel m = new GaussianPSFModel(s0, s1);
		double[] e = new double[o.length];
		// Note that the Gaussian2DFunction has 0,0 at the centre of a pixel.
		// The model has 0.5,0.5 at the centre so add an offset.
		m.create2D(e, width, height, sum, x0 + 0.5, x1 + 0.5, false);

		// Since the model only computes within +/- 5 sd only check for equality 
		// when the model is not zero (and there is a reasonable amount of signal)

		for (int i = 0; i < e.length; i++)
			if (e[i] > 1e-2) // Only check where there is a reasonable amount of signal
			{
				double error = DoubleEquality.relativeError(e[i], o[i]);
				// We expect a small error since the ErfGaussian2DFunction uses a 
				// fast approximation of the Erf(..) (the error function). The PSFModel
				// uses the Apache commons implementation.
				if (error > 5e-4)
					Assert.fail(String.format("[%d] %s != %s  error = %f\n", i, Double.toString(e[i]),
							Double.toString(o[i]), error));
			}
	}
}
