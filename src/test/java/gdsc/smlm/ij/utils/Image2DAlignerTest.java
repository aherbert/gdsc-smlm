package gdsc.smlm.ij.utils;

import org.junit.Assert;
import org.junit.Test;

import gdsc.smlm.function.StandardFloatValueProcedure;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.function.gaussian.GaussianFunctionFactory;

public class Image2DAlignerTest
{
	// TODO - Make this test the StackAligner with sub-pixel accuracy and non power of 2 images

	// Check the alignment of the ImageProcessor and Image2D is the same with padding and windowing

	private FloatImage2D createData(int x, int y, double cx, double cy)
	{
		Gaussian2DFunction f = GaussianFunctionFactory.create2D(1, x, y, GaussianFunctionFactory.FIT_ERF_FREE_CIRCLE,
				null);
		double[] a = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
		a[Gaussian2DFunction.SIGNAL] = 1;
		a[Gaussian2DFunction.X_POSITION] = cx;
		a[Gaussian2DFunction.Y_POSITION] = cy;
		a[Gaussian2DFunction.X_SD] = 1.2;
		a[Gaussian2DFunction.Y_SD] = 1.1;
		StandardFloatValueProcedure p = new StandardFloatValueProcedure();
		return new FloatImage2D(x, y, p.getValues(f, a));
	}

	@Test
	public void canCorrelate()
	{
		int maxx = 16;
		int maxy = 16;
		double cx = (maxx - 1) / 2.0;
		double cy = (maxy - 1) / 2.0;

		double[] shift = new double[] { 0, 1, 1.5, 2 };

		Image2D reference = createData(maxx, maxy, cx, cy);
		Image2DAligner a = new Image2DAligner();
		a.setReference(reference);

		for (double sy : shift)
			for (double sx : shift)
			{
				canCorrelate(a, maxx, maxy, cx, cy, cx + sx, cy + sy, 0.25, 0, 1);
				canCorrelate(a, maxx, maxy, cx, cy, cx + sx, cy + sy, 0.25, 5, 0.25);
			}
	}

	private void canCorrelate(Image2DAligner a, int maxx, int maxy, double cx1, double cy1, double cx2, double cy2,
			double window, int refinements, double tolerance)
	{
		double[] result;
		Image2D c;
		Image2D target = createData(maxx, maxy, cx2, cy2);

		//Utils.display("Ref", reference.getImageProcessor());
		//Utils.display("Tar", target.getImageProcessor());

		a.setEdgeWindow(window);
		//a.setSearchMode(gdsc.smlm.ij.utils.StackAligner.SearchMode.BINARY);
		//a.setRelativeThreshold(1e-6);

		double[] e = new double[] { cx1 - cx2, cy1 - cy2 };
		int cx = maxx / 2 - (int) Math.round(e[0]);
		int cy = maxy / 2 - (int) Math.round(e[1]);
		int index = target.getIndex(cx, cy);

		// Debug the convergence
		//for (int i = 0; i <= refinements; i++)
		//{
		//	result = a.align(target.copy(), i, error);
		//	c = a.getCorrelation();
		//
		//	System.out.printf("e %s %g, o %s\n", java.util.Arrays.toString(e), c.get(index),
		//			java.util.Arrays.toString(result));
		//}

		// Test
		result = a.align(target, refinements);
		c = a.getCorrelation();
		System.out.printf("e %s %g, o %s\n", java.util.Arrays.toString(e), c.get(index),
				java.util.Arrays.toString(result));

		for (int i = 0; i < 2; i++)
			Assert.assertEquals(e[i], result[i], tolerance);
	}
}
