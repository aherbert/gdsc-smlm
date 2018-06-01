package gdsc.smlm.function;

import java.awt.Shape;
import java.awt.geom.Ellipse2D;
import java.util.Arrays;

import org.junit.Assert;
import org.junit.Test;

import gdsc.core.utils.Maths;
import gdsc.core.utils.SimpleArrayUtils;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import gdsc.smlm.results.Gaussian2DPeakResultHelper;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.ShapeRoi;
import ij.process.FloatProcessor;

public class FunctionHelperTest
{
	@Test
	public void canGetMeanValue()
	{
		int n = 10;
		double[] values = SimpleArrayUtils.newArray(n, 1.0, 1.0);
		Assert.assertEquals(10, FunctionHelper.getMeanValue(values.clone(), 0), 0);
		double total = sum(values, n);
		Assert.assertEquals(total / n, FunctionHelper.getMeanValue(values.clone(), 1), 0);
		for (int i = 1; i < n; i++)
		{
			double sum = sum(values, i);
			Assert.assertEquals(sum / i, FunctionHelper.getMeanValue(values.clone(), sum / total), 0);
		}
	}

	private double sum(double[] values, int top)
	{
		double sum = 0;
		int n = 0;
		for (int i = values.length; i-- > 0 && n < top; n++)
			sum += values[i];
		return sum;
	}

	@Test
	public void canGetFractionalMeanValue()
	{
		int n = 10;
		double[] values = SimpleArrayUtils.newDoubleArray(n, 1.0);
		Assert.assertEquals(1, FunctionHelper.getMeanValue(values.clone(), 0), 0);
		for (int i = 1; i < n; i++)
		{
			double f = (double) i / n;
			Assert.assertEquals(1, FunctionHelper.getMeanValue(values.clone(), f), 0);
			Assert.assertEquals(1, FunctionHelper.getMeanValue(values.clone(), f - 0.5), 0);
		}

		Arrays.fill(values, 5, n, 2);
		// sum = 5*1 + 5*2 = 15
		Assert.assertEquals(2, FunctionHelper.getMeanValue(values.clone(), 5.0 / 15), 0);
		Assert.assertEquals(2, FunctionHelper.getMeanValue(values.clone(), 10.0 / 15), 0);
		Assert.assertEquals(11.0 / 6, FunctionHelper.getMeanValue(values.clone(), 11.0 / 15), 0);
		Assert.assertEquals(11.5 / 6.5, FunctionHelper.getMeanValue(values.clone(), 11.5 / 15), 0);
	}

	@Test
	public void canGetXValue()
	{
		int n = 10;
		double[] values = SimpleArrayUtils.newArray(n, 1.0, 1.0);
		Assert.assertEquals(0, FunctionHelper.getXValue(values.clone(), 0), 0);
		Assert.assertEquals(n, FunctionHelper.getXValue(values.clone(), 1), 0);
		double total = sum(values, n);
		for (int i = 1; i < n; i++)
		{
			double sum = sum(values, i);
			Assert.assertEquals(i, FunctionHelper.getXValue(values.clone(), sum / total), 0);
		}
	}

	@Test
	public void canGetFractionalXValue()
	{
		int n = 10;
		double[] values = SimpleArrayUtils.newDoubleArray(n, 1.0);
		Assert.assertEquals(0, FunctionHelper.getXValue(values.clone(), 0), 0);
		Assert.assertEquals(n, FunctionHelper.getXValue(values.clone(), 1), 0);
		for (int i = 1; i < n; i++)
		{
			double f = (double) i / n;
			Assert.assertEquals(i, FunctionHelper.getXValue(values.clone(), f), 0);
			Assert.assertEquals(i - 0.5, FunctionHelper.getXValue(values.clone(), f - 0.05), 1e-8);
		}

		Arrays.fill(values, 5, n, 2);
		// sum = 5*1 + 5*2 = 15
		Assert.assertEquals(2.5, FunctionHelper.getXValue(values.clone(), 5.0 / 15), 0);
		Assert.assertEquals(5, FunctionHelper.getXValue(values.clone(), 10.0 / 15), 0);
		Assert.assertEquals(6, FunctionHelper.getXValue(values.clone(), 11.0 / 15), 0);
		Assert.assertEquals(6.5, FunctionHelper.getXValue(values.clone(), 11.5 / 15), 0);
	}

	@Test
	public void canGetMeanValueForGaussian()
	{
		// This does not work. 
		// The number of samples to achieve the target intensity is much 
		// greater than the area of the ellipse covered by the Gaussian.
		// Plot tbe drawn gaussian.
		// Overlay an ellipse for each range of SD.
		// Create a mask using all pixels required to achieve the sum.
		// Does it roughly fit the ellipse?

		// This shows that the Gaussian2DPeakResultHelper.getMeanSignal is not working
		// The 2D cumulative needs to be computed. I am using a 1D cumulative.
		// map the 2D back into a smaller 1D range: range = 1 / sqrt(range^2 * 2)
		for (int range = 1; range <= 3; range++)
		{
			double a = Gaussian2DPeakResultHelper.cumulative(range);
			double b = Math.PI * Maths.pow2(Gaussian2DPeakResultHelper.cumulative(range) / 2);

			System.out.printf("%d  %g  %g  %g\n", range, a, b, a / b);
		}

		float intensity = 100;
		float sx = 20f;
		float sy = 20f;
		int size = 1 + 2 * (int) Math.ceil(Math.max(sx, sy) * 4);
		float[] a = Gaussian2DPeakResultHelper.createParams(0.f, intensity, size / 2f, size / 2f, 0.f, sx, sy, 0);
		Gaussian2DFunction f = GaussianFunctionFactory.create2D(1, size, size, GaussianFunctionFactory.FIT_FREE_CIRCLE
		//| GaussianFunctionFactory.FIT_SIMPLE
				, null);
		double[] values = f.computeValues(SimpleArrayUtils.toDouble(a));
		ImagePlus imp = new ImagePlus("gauss", new FloatProcessor(size, size, values));
		double cx = size / 2.;
		Shape shape = new Ellipse2D.Double(cx - sx, cx - sy, 2 * sx, 2 * sy);
		imp.setRoi(new ShapeRoi(shape));
		//imp.setRoi(new EllipseRoi(cx - sx, cx - sy, cx + sx, cx + sy, 1));
		IJ.save(imp, "/Users/ah403/1.tif");
		//imp.close();
		double scale = Maths.sum(values) / intensity;
		for (int range = 1; range <= 3; range++)
		{
			double e = Gaussian2DPeakResultHelper.getMeanSignal(intensity, sx, sy, range);
			double o = FunctionHelper.getMeanValue(values.clone(),
					scale * Gaussian2DPeakResultHelper.cumulative(range));
			System.out.printf("%d  %g %g  %g\n", range, e, o, e / o);
			//Assert.assertEquals(e, o, e*1e-2);
		}
	}
}
