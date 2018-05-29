package gdsc.smlm.model;

import org.junit.Assert;
import org.junit.Test;

import gdsc.smlm.function.Gradient1Procedure;
import gdsc.smlm.function.gaussian.AstigmatismZModel;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.function.gaussian.HoltzerAstigmatismZModel;
import gdsc.smlm.function.gaussian.erf.SingleAstigmatismErfGaussian2DFunction;

public class PSFModelGradient1FunctionTest
{
	@Test
	public void canComputeValueAndGradient()
	{
		// Use a reasonable z-depth function from the Smith, et al (2010) paper (page 377)
		double sx = 1.08;
		double sy = 1.01;
		double gamma = 0.389;
		double d = 0.531;
		double Ax = -0.0708;
		double Bx = -0.073;
		double Ay = 0.164;
		double By = 0.0417;
		AstigmatismZModel zModel = HoltzerAstigmatismZModel.create(sx, sy, gamma, d, Ax, Bx, Ay, By);

		// Small size ensure the PSF model covers the entire image
		int maxx = 11, maxy = 11;
		final double[] ve = new double[maxx * maxy];
		final double[] vo = new double[maxx * maxy];
		final double[][] ge = new double[maxx * maxy][];
		final double[][] go = new double[maxx * maxy][];

		PSFModelGradient1Function psf = new PSFModelGradient1Function(new GaussianPSFModel(zModel), maxx, maxy);
		Gaussian2DFunction f = new SingleAstigmatismErfGaussian2DFunction(maxx, maxy, zModel);
		double[] a2 = new double[Gaussian2DFunction.PARAMETERS_PER_PEAK + 1];

		double c = maxx * 0.5;
		for (int i = -1; i <= 1; i++)
		{
			double x0 = c + i * 0.33;
			for (int j = -1; j <= 1; j++)
			{
				double x1 = c + j * 0.33;
				for (int k = -1; k <= 1; k++)
				{
					double x2 = k * 0.33;
					for (double in : new double[] { 23.2, 405.67 })
					{
						// Background is constant for gradients so just use 1 value
						double[] a = new double[] { 2.2, in, x0, x1, x2 };
						psf.initialise1(a);
						psf.forEach(new Gradient1Procedure()
						{
							int i = 0;

							public void execute(double value, double[] dy_da)
							{
								vo[i] = value;
								go[i] = dy_da.clone();
								i++;
							}
						});
						a2[Gaussian2DFunction.BACKGROUND] = a[0];
						a2[Gaussian2DFunction.SIGNAL] = a[1];
						a2[Gaussian2DFunction.X_POSITION] = a[2] - 0.5;
						a2[Gaussian2DFunction.Y_POSITION] = a[3] - 0.5;
						a2[Gaussian2DFunction.Z_POSITION] = a[4];
						f.initialise1(a2);
						f.forEach(new Gradient1Procedure()
						{
							int i = 0;

							public void execute(double value, double[] dy_da)
							{
								ve[i] = value;
								ge[i] = dy_da.clone();
								i++;
							}
						});

						for (int ii = 0; ii < ve.length; ii++)
						{
							Assert.assertEquals(ve[ii], vo[ii], ve[ii] * 1e-8);
							for (int l = 0; l < 5; l++)
							{
								Assert.assertEquals(ge[ii][l], go[ii][l], Math.abs(ge[ii][l] * 1e-8));
							}
						}
					}
				}
			}
		}
	}
}
