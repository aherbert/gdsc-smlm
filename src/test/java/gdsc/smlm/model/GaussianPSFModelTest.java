package gdsc.smlm.model;

import java.util.Arrays;

import org.junit.Assert;
import org.junit.Test;

import gdsc.core.utils.DoubleEquality;
import gdsc.core.utils.Maths;
import gdsc.smlm.function.gaussian.AstigmatismZModel;
import gdsc.smlm.function.gaussian.HoltzerAstigmatismZModel;

public class GaussianPSFModelTest
{
	@Test
	public void canComputeValue()
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

		int maxx = 21, maxy = 21;
		double[] e = new double[maxx * maxy];
		double[] o = new double[maxx * maxy];
		double[] o2 = new double[maxx * maxy];
		double[][] g = new double[maxx * maxy][];

		PSFModel psf = new GaussianPSFModel(zModel);

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
					Arrays.fill(e, 0);
					psf.create3D(e, maxx, maxy, 1, x0, x1, x2, false);
					psf.getValue(maxx, maxy, x0, x1, x2, o);
					Assert.assertEquals(1, Maths.sum(o), 1e-3);
					psf.getValueAndGradient(maxx, maxy, x0, x1, x2, o2, g);
					for (int ii = 0; ii < e.length; ii++)
					{
						if (e[ii] == 0)
						{
							// Same insertion into the blank data region
							Assert.assertTrue(o[ii] == 0);
						}
						else if (e[ii] > 1e-8) // Only check where there is a reasonable amount of signal
						{
							double error = DoubleEquality.relativeError(e[ii], o[ii]);
							//System.out.printf("[%d,%d]   %g == %g    %g\n", ii/maxx, ii%maxx, e[ii], o[ii], error);
							// We expect a small error since the ErfGaussian2DFunction uses a 
							// fast approximation of the Erf(..) (the error function). The PSFModel
							// uses the Apache commons implementation.
							if (error > 1e-8)
								Assert.fail(String.format("[%d] %s != %s  error = %f\n", ii, Double.toString(e[ii]),
										Double.toString(o[ii]), error));
						}
					}
				}
			}
		}
	}

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

		//zModel  = new NullAstigmatismZModel(1.2, 1.2);

		int maxx = 21, maxy = 21;
		double[] o = new double[maxx * maxy];
		double[] o2 = new double[maxx * maxy];
		double[] o3 = new double[maxx * maxy];
		double[][] g = new double[maxx * maxy][];
		double[][] g2 = new double[maxx * maxy][3];
		double[] dx = new double[] { 1e-4, 1e-4, 1e-4 };

		PSFModel psf = new GaussianPSFModel(zModel);

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
					psf.getValue(maxx, maxy, x0, x1, x2, o);
					psf.getValueAndGradient(maxx, maxy, x0, x1, x2, o2, g);
					Assert.assertArrayEquals(o, o2, 0);

					// Check verses default manual gradient
					psf.computeValueAndGradient(maxx, maxy, x0, x1, x2, o3, g2, dx);
					Assert.assertArrayEquals(o, o3, 0);

					for (int ii = 0; ii < o.length; ii++)
						if (o[ii] > 1e-4) // Only check where there is a reasonable amount of signal
						{
							for (int l = 0; l < 3; l++)
							{
								double error = DoubleEquality.relativeError(g[ii][l], g2[ii][l]);
								//System.out.printf("[%d,%d]   %g == %g    %g\n", ii, l, g[ii][l], g2[ii][l], error);
								if (error > 5e-3)
									Assert.fail(String.format("[%d] %s != %s  error = %f\n", ii,
											Double.toString(g[ii][l]), Double.toString(g2[ii][l]), error));

							}
						}
				}
			}
		}
	}
}
