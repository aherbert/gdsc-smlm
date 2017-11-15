package gdsc.smlm.ij.utils;

import org.junit.Test;

import gdsc.smlm.function.StandardFloatValueProcedure;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import gdsc.smlm.function.gaussian.QuadraticAstigmatismZModel;

public class StackAlignerTest
{
	// TODO - Make this test the StackAligner with sub-pixel accuracy and non power of 2 images

	// Check the alignment of the ImageStack and Image3D is the same with padding and windowing

	final static double gamma = 5;
	final static int zDepth = 10;
	protected QuadraticAstigmatismZModel zModel = new QuadraticAstigmatismZModel(gamma, zDepth);

	private FloatImage3D createData(int x, int y, int z, double cx, double cy, double cz)
	{
		Gaussian2DFunction f = GaussianFunctionFactory.create2D(1, x, y, GaussianFunctionFactory.FIT_ASTIGMATISM,
				zModel);
		int length = x * y;
		float[] data = new float[z * length];
		double[] a = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
		a[Gaussian2DFunction.SIGNAL] = 1;
		a[Gaussian2DFunction.X_POSITION] = cx;
		a[Gaussian2DFunction.Y_POSITION] = cy;
		a[Gaussian2DFunction.X_SD] = 1;
		a[Gaussian2DFunction.Y_SD] = 1;
		StandardFloatValueProcedure p = new StandardFloatValueProcedure();
		for (int zz = 0; zz < z; zz++)
		{
			double dz = zz - cz;
			if (zz == 0 || zz == z - 1)
				System.out.printf("%f  %f %f\n", dz, zModel.getSx(dz), zModel.getSy(dz));
			a[Gaussian2DFunction.Z_POSITION] = dz;
			p.getValues(f, a, data, zz * length);
		}
		return new FloatImage3D(x, y, z, data);
	}

	@Test
	public void canCorrelate()
	{
		int maxz = 63;
		canCorrelate(16, 16, maxz, 7.5, 7.5, (maxz - 1) / 2.0, 8, 8, (maxz +2) / 2.0, 0.1, 10, 1e-2);
	}

	private void canCorrelate(int maxx, int maxy, int maxz, double cx1, double cy1, double cz1, double cx2, double cy2,
			double cz2, double window, int refinements, double error)
	{
		Image3D reference = createData(maxx, maxy, maxz, cx1, cy1, cz1);
		Image3D target = createData(maxx, maxy, maxz, cx2, cy2, cz2);

		//Utils.display("Ref", reference.getImageStack());
		//Utils.display("Tar", target.getImageStack());

		StackAligner a = new StackAligner(window);
		a.setReference(reference);
		//a.setSearchMode(SearchMode.BINARY);
		//a.setRelativeThreshold(1e-6);

		double[] e = new double[] { cx1 - cx2, cy1 - cy2, cz1 - cz2 };
		int cx = maxx / 2 - (int) Math.round(e[0]);
		int cy = maxy / 2 - (int) Math.round(e[1]);
		int cz = maxz / 2 - (int) Math.round(e[2]);
		int index = target.getIndex(cx, cy, cz);

		for (int i = 0; i <= refinements; i++)
		{
			double[] result = a.align(target.copy(), i, error);
			Image3D c = a.getCorrelation();

			System.out.printf("e %s %g, o %s\n", java.util.Arrays.toString(e), c.get(index),
					java.util.Arrays.toString(result));
			//for (int i = 0; i < 3; i++)
			//	Assert.assertEquals(e[i], result[i], 0.25);
		}
	}
}
