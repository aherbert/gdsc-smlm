package gdsc.smlm.function;

import java.util.Arrays;

import org.junit.Assert;
import org.junit.Test;

import gdsc.core.utils.DoubleEquality;
import gdsc.core.utils.Maths;
import gdsc.core.utils.SimpleArrayUtils;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import gdsc.smlm.results.Gaussian2DPeakResultHelper;

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
		float intensity = 100;
		// Realistic standard deviations.
		// Only test the highest
		float[] s = { 1, 1.2f, 1.5f, 2, 2.5f };
		int n_1 = s.length - 1;
		// Flag to indicate that all levels should be run and the difference reported
		boolean debug = false;
		for (int i = (debug) ? 0 : n_1; i < s.length; i++)
		{
			float sx = s[i];
			for (int j = i; j < s.length; j++)
			{
				float sy = s[j];
				int size = 1 + 2 * (int) Math.ceil(Math.max(sx, sy) * 4);
				float[] a = Gaussian2DPeakResultHelper.createParams(0.f, intensity, size / 2f, size / 2f, 0.f, sx, sy,
						0);
				Gaussian2DFunction f = GaussianFunctionFactory.create2D(1, size, size,
						GaussianFunctionFactory.FIT_FREE_CIRCLE
						//| GaussianFunctionFactory.FIT_SIMPLE
						, null);
				double[] values = f.computeValues(SimpleArrayUtils.toDouble(a));
				//ImagePlus imp = new ImagePlus("gauss", new FloatProcessor(size, size, values));
				//double cx = size / 2.;
				//Shape shape = new Ellipse2D.Double(cx - sx, cx - sy, 2 * sx, 2 * sy);
				//imp.setRoi(new ShapeRoi(shape));
				//IJ.save(imp, "/Users/ah403/1.tif");
				double scale = Maths.sum(values) / intensity;
				for (int range = 1; range <= 3; range++)
				{
					double e = Gaussian2DPeakResultHelper.getMeanSignalUsingR(intensity, sx, sy, range);
					double o = FunctionHelper.getMeanValue(values.clone(),
							scale * Gaussian2DPeakResultHelper.cumulative2D(range));
					if (debug)
						System.out.printf("%g,%g   %d  %g %g  %g\n", sx, sy, range, e, o,
								DoubleEquality.relativeError(e, o));
					// Only test the highest
					if (i == n_1)
						Assert.assertEquals(e, o, e * 0.025);
				}
			}
		}
	}
}