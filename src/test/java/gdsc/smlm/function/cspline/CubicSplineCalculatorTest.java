package gdsc.smlm.function.cspline;

import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.math.interpolation.CubicSplinePosition;
import gdsc.core.math.interpolation.CustomTricubicFunction;
import gdsc.core.math.interpolation.CustomTricubicInterpolatingFunction;
import gdsc.core.math.interpolation.CustomTricubicInterpolator;
import gdsc.core.utils.Maths;
import gdsc.core.utils.SimpleArrayUtils;

public class CubicSplineCalculatorTest
{
	@Test
	public void canComputeCoefficientsForDistanceFunction()
	{
		double[] e = new double[64];
		int c = 0;
		for (int k = 0; k < 4; k++)
			for (int j = 0; j < 4; j++)
				for (int i = 0; i < 4; i++)
					e[c++] = Math.sqrt(i * i + j * j + k * k);
		CustomTricubicFunction f = CustomTricubicFunction.create(e);
		CubicSplinePosition[] s = new CubicSplinePosition[4];
		for (int i = 0; i < 4; i++)
			s[i] = new CubicSplinePosition((double) i / 3);
		double[][][] value = new double[4][4][4];
		double[] b = new double[64];
		c = 0;
		for (int k = 0; k < 4; k++)
			for (int j = 0; j < 4; j++)
				for (int i = 0; i < 4; i++)
				{
					value[i][j][k] = f.value(s[i], s[j], s[k]);
					b[c++] = value[i][j][k]; 
				}

		CubicSplineCalculator calc = new CubicSplineCalculator();
		double[] o = calc.compute(value);
		Assert.assertArrayEquals(e, o, 1e-6);
		
		o = calc.compute(b);
		Assert.assertArrayEquals(e, o, 1e-6);
	}

	@Test
	public void canComputeCoefficientsForGaussianFunction()
	{
		int x = 4, y = 4, z = 4;
		double xscale = 1, yscale = 0.5, zscale = 2.0;
		double[] xval = SimpleArrayUtils.newArray(x, 0, xscale);
		double[] yval = SimpleArrayUtils.newArray(y, 0, yscale);
		double[] zval = SimpleArrayUtils.newArray(z, 0, zscale);
		double[][][] fval = createData(x, y, z, null);
		CustomTricubicInterpolatingFunction f1 = new CustomTricubicInterpolator().interpolate(xval, yval, zval, fval);

		double[] e = f1.getCoefficients(1, 1, 1);
		
		CustomTricubicFunction f = CustomTricubicFunction.create(e);
		CubicSplinePosition[] s = new CubicSplinePosition[4];
		for (int i = 0; i < 4; i++)
			s[i] = new CubicSplinePosition((double) i / 3);
		double[][][] value = new double[4][4][4];
		double[] b = new double[64];
		int c = 0;
		for (int k = 0; k < 4; k++)
			for (int j = 0; j < 4; j++)
				for (int i = 0; i < 4; i++)
				{
					value[i][j][k] = f.value(s[i], s[j], s[k]);
					b[c++] = value[i][j][k]; 
				}
		
		CubicSplineCalculator calc = new CubicSplineCalculator();
		double[] o = calc.compute(value);
		Assert.assertArrayEquals(e, o, 1e-6);
		
		o = calc.compute(b);
		Assert.assertArrayEquals(e, o, 1e-6);
	}

	double[][][] createData(int x, int y, int z, RandomGenerator r)
	{
		double[][][] fval = new double[x][y][z];
		// Create a 2D Gaussian
		double s = 1.0;
		double cx = x / 2.0;
		double cy = y / 2.0;
		if (r != null)
		{
			s += r.nextDouble() - 0.5;
			cx += r.nextDouble() - 0.5;
			cy += r.nextDouble() - 0.5;
		}
		double[] otherx = new double[x];
		for (int zz = 0; zz < z; zz++)
		{
			double s2 = 2 * s * s;
			for (int xx = 0; xx < x; xx++)
				otherx[xx] = Maths.pow2(xx - cx) / s2;
			for (int yy = 0; yy < y; yy++)
			{
				double othery = Maths.pow2(yy - cy) / s2;
				for (int xx = 0; xx < x; xx++)
				{
					fval[xx][yy][zz] = Math.exp(otherx[xx] + othery);
				}
			}
			// Move Gaussian
			s += 0.1;
			cx += 0.1;
			cy -= 0.05;
		}
		return fval;
	}
}
