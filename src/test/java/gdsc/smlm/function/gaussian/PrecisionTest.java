/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 * 
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */
package gdsc.smlm.function.gaussian;

import org.apache.commons.math3.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

/**
 * Contains tests for the Gaussian functions in single or double precision
 * <p>
 * The tests show that there is very little (if any) time penalty when using double precision for the calculations.
 * However the precision of the single-precision functions is 1e-4 when using reasonable Gaussian parameters. This could
 * effect the convergence of optimisers/fitters if using single precision math.
 */
public class PrecisionTest
{
	int Single = 1;
	int Double = 2;

	private int MAX_ITER = 200000;
	double SPEED_UP_FACTOR = 1.1;

	int maxx = 10;
	// Use realistic values for a camera with a bias of 500
	static double[] params2 = new double[] { 500.23, 300.12, 0, 5.12, 5.23, 1.11, 1.11 };
	static float[] params1 = toFloat(params2);

	// Stripped down Gaussian functions copied from the gdsc.smlm.fitting.function.gaussian package
	public abstract class Gaussian
	{
		public static final int BACKGROUND = 0;
		public static final int AMPLITUDE = 1;
		public static final int ANGLE = 2;
		public static final int X_POSITION = 3;
		public static final int Y_POSITION = 4;
		public static final int X_SD = 5;
		public static final int Y_SD = 6;

		int maxx;

		public Gaussian(int maxx)
		{
			this.maxx = maxx;
		}

		public void setMaxX(int maxx)
		{
			this.maxx = maxx;
		}
	}

	public interface DoublePrecision
	{
		public void setMaxX(int maxx);

		public void initialise(double[] a);

		public double eval(final int x, final double[] dyda);

		public double eval(final int x);
	}

	public interface SinglePrecision
	{
		public void setMaxX(int maxx);

		public void initialise(float[] a);

		public float eval(final int x, final float[] dyda);

		public float eval(final int x);
	}

	public class DoubleCircularGaussian extends Gaussian implements DoublePrecision
	{
		double background;
		double amplitude;
		double x0pos;
		double x1pos;

		double aa;
		double aa2;
		double ax;

		public DoubleCircularGaussian(int maxx)
		{
			super(maxx);
		}

		public void initialise(double[] a)
		{
			background = a[BACKGROUND];
			amplitude = a[AMPLITUDE];
			x0pos = a[X_POSITION];
			x1pos = a[Y_POSITION];

			final double sx = a[X_SD];
			final double sx2 = sx * sx;
			final double sx3 = sx2 * sx;

			aa = -0.5 / sx2;
			aa2 = -2.0 * aa;

			// For the x-width gradient
			ax = 1.0 / sx3;
		}

		public double eval(final int x, final double[] dyda)
		{
			dyda[0] = 1.0;

			final int x1 = x / maxx;
			final int x0 = x % maxx;

			return background + gaussian(x0, x1, dyda);
		}

		private double gaussian(final int x0, final int x1, final double[] dy_da)
		{
			final double h = amplitude;

			final double dx = x0 - x0pos;
			final double dy = x1 - x1pos;
			final double dx2dy2 = dx * dx + dy * dy;

			dy_da[1] = FastMath.exp(aa * (dx2dy2));
			final double y = h * dy_da[1];
			final double yaa2 = y * aa2;
			dy_da[2] = yaa2 * dx;
			dy_da[3] = yaa2 * dy;

			dy_da[4] = y * (ax * (dx2dy2));

			return y;
		}

		public double eval(final int x)
		{
			final int x1 = x / maxx;
			final int x0 = x % maxx;

			final double dx = x0 - x0pos;
			final double dy = x1 - x1pos;

			return background + amplitude * FastMath.exp(aa * (dx * dx + dy * dy));
		}
	}

	public class SingleCircularGaussian extends Gaussian implements SinglePrecision
	{
		float background;
		float amplitude;
		float x0pos;
		float x1pos;

		float aa;
		float aa2;
		float ax;

		public SingleCircularGaussian(int maxx)
		{
			super(maxx);
		}

		public void initialise(float[] a)
		{
			background = a[BACKGROUND];
			amplitude = a[AMPLITUDE];
			x0pos = a[X_POSITION];
			x1pos = a[Y_POSITION];

			final float sx = a[X_SD];
			final float sx2 = sx * sx;
			final float sx3 = sx2 * sx;

			aa = -0.5f / sx2;
			aa2 = -2.0f * aa;

			ax = 1.0f / sx3;
		}

		public float eval(final int x, final float[] dyda)
		{
			dyda[0] = 1.0f;

			final int x1 = x / maxx;
			final int x0 = x % maxx;

			return background + gaussian(x0, x1, dyda);
		}

		private float gaussian(final int x0, final int x1, final float[] dy_da)
		{
			final float h = amplitude;

			final float dx = x0 - x0pos;
			final float dy = x1 - x1pos;
			final float dx2dy2 = dx * dx + dy * dy;

			dy_da[1] = (float) FastMath.exp(aa * (dx2dy2));
			final float y = h * dy_da[1];
			final float yaa2 = y * aa2;
			dy_da[2] = yaa2 * dx;
			dy_da[3] = yaa2 * dy;

			dy_da[4] = y * (ax * (dx2dy2));

			return y;
		}

		public float eval(final int x)
		{
			final int x1 = x / maxx;
			final int x0 = x % maxx;

			final float dx = x0 - x0pos;
			final float dy = x1 - x1pos;

			return background + amplitude * (float) (FastMath.exp(aa * (dx * dx + dy * dy)));
		}
	}

	public class DoubleFixedGaussian extends Gaussian implements DoublePrecision
	{
		double width;

		double background;
		double amplitude;
		double x0pos;
		double x1pos;

		double aa;
		double aa2;

		public DoubleFixedGaussian(int maxx)
		{
			super(maxx);
		}

		public void initialise(double[] a)
		{
			background = a[BACKGROUND];
			amplitude = a[AMPLITUDE];
			x0pos = a[X_POSITION];
			x1pos = a[Y_POSITION];
			width = a[X_SD];

			final double sx = a[X_SD];
			final double sx2 = sx * sx;

			aa = -0.5 / sx2;
			aa2 = -2.0 * aa;
		}

		public double eval(final int x, final double[] dyda)
		{
			dyda[0] = 1.0;

			final int x1 = x / maxx;
			final int x0 = x % maxx;

			return background + gaussian(x0, x1, dyda);
		}

		private double gaussian(final int x0, final int x1, final double[] dy_da)
		{
			final double h = amplitude;

			final double dx = x0 - x0pos;
			final double dy = x1 - x1pos;

			dy_da[1] = FastMath.exp(aa * (dx * dx + dy * dy));
			final double y = h * dy_da[1];
			final double yaa2 = y * aa2;
			dy_da[2] = yaa2 * dx;
			dy_da[3] = yaa2 * dy;

			return y;
		}

		public double eval(final int x)
		{
			final int x1 = x / maxx;
			final int x0 = x % maxx;

			final double dx = x0 - x0pos;
			final double dy = x1 - x1pos;

			return background + amplitude * FastMath.exp(aa * (dx * dx + dy * dy));
		}
	}

	public class SingleFixedGaussian extends Gaussian implements SinglePrecision
	{
		float width;

		float background;
		float amplitude;
		float x0pos;
		float x1pos;

		float aa;
		float aa2;

		public SingleFixedGaussian(int maxx)
		{
			super(maxx);
		}

		public void initialise(float[] a)
		{
			background = a[BACKGROUND];
			amplitude = a[AMPLITUDE];
			x0pos = a[X_POSITION];
			x1pos = a[Y_POSITION];
			width = a[X_SD];

			final float sx = a[X_SD];
			final float sx2 = sx * sx;

			aa = -0.5f / sx2;
			aa2 = -2.0f * aa;
		}

		public float eval(final int x, final float[] dyda)
		{
			dyda[0] = 1.0f;

			final int x1 = x / maxx;
			final int x0 = x % maxx;

			return background + gaussian(x0, x1, dyda);
		}

		private float gaussian(final int x0, final int x1, final float[] dy_da)
		{
			final float h = amplitude;

			final float dx = x0 - x0pos;
			final float dy = x1 - x1pos;

			dy_da[1] = (float) (FastMath.exp(aa * (dx * dx + dy * dy)));
			final float y = h * dy_da[1];
			final float yaa2 = y * aa2;
			dy_da[2] = yaa2 * dx;
			dy_da[3] = yaa2 * dy;

			return y;
		}

		public float eval(final int x)
		{
			final int x1 = x / maxx;
			final int x0 = x % maxx;

			final float dx = x0 - x0pos;
			final float dy = x1 - x1pos;

			return background + amplitude * (float) (FastMath.exp(aa * (dx * dx + dy * dy)));
		}
	}

	@Test
	public void circularFunctionPrecisionIs3sf()
	{
		functionsComputeSameValue(maxx, new SingleCircularGaussian(maxx), new DoubleCircularGaussian(maxx), 1e-3);
	}

	@Test
	public void circularFunctionPrecisionIs4sf()
	{
		functionsComputeSameValue(maxx, new SingleCircularGaussian(maxx), new DoubleCircularGaussian(maxx), 1e-4);
	}

	@Test(expected = java.lang.AssertionError.class)
	public void circularFunctionPrecisionIsNot5sf()
	{
		functionsComputeSameValue(maxx, new SingleCircularGaussian(maxx), new DoubleCircularGaussian(maxx), 1e-5);
	}

	@Test(expected = java.lang.AssertionError.class)
	public void circularFunctionsPrecisionIsNot3sfAtLargeXY()
	{
		int maxx = this.maxx;
		try
		{
			for (;;)
			{
				maxx *= 2;
				System.out.printf("maxx = %d\n", maxx);
				functionsComputeSameValue(maxx, new SingleCircularGaussian(maxx), new DoubleCircularGaussian(maxx),
						1e-3);
			}
		}
		catch (AssertionError e)
		{
			System.out.println(e.getMessage());
			//e.printStackTrace();
			throw e;
		}
	}

	@Test(expected = java.lang.AssertionError.class)
	public void circularSinglePrecisionIsNotMuchFasterWithGradients()
	{
		singlePrecisionIsFasterWithGradients(maxx, new SingleCircularGaussian(maxx), new DoubleCircularGaussian(maxx),
				false);
	}

	@Test(expected = java.lang.AssertionError.class)
	public void circularSinglePrecisionIsNotMuchFaster()
	{
		singlePrecisionIsFaster(maxx, new SingleCircularGaussian(maxx), new DoubleCircularGaussian(maxx), false);
	}

	@Test(expected = java.lang.AssertionError.class)
	public void circularSinglePrecisionIsNotMuchFasterWithGradientsNoSum()
	{
		singlePrecisionIsFasterWithGradients(maxx, new SingleCircularGaussian(maxx), new DoubleCircularGaussian(maxx),
				true);
	}

	@Test(expected = java.lang.AssertionError.class)
	public void circularSinglePrecisionIsNotMuchFasterNoSum()
	{
		singlePrecisionIsFaster(maxx, new SingleCircularGaussian(maxx), new DoubleCircularGaussian(maxx), true);
	}

	@Test
	public void fixedFunctionPrecisionIs3sf()
	{
		functionsComputeSameValue(maxx, new SingleFixedGaussian(maxx), new DoubleFixedGaussian(maxx), 1e-3);
	}

	@Test
	public void fixedFunctionPrecisionIs4sf()
	{
		functionsComputeSameValue(maxx, new SingleFixedGaussian(maxx), new DoubleFixedGaussian(maxx), 1e-4);
	}

	@Test(expected = java.lang.AssertionError.class)
	public void fixedFunctionPrecisionIsNot5sf()
	{
		functionsComputeSameValue(maxx, new SingleFixedGaussian(maxx), new DoubleFixedGaussian(maxx), 1e-5);
	}

	@Test(expected = java.lang.AssertionError.class)
	public void fixedFunctionsPrecisionIsNot3sfAtLargeXY()
	{
		int maxx = this.maxx;
		try
		{
			for (;;)
			{
				maxx *= 2;
				System.out.printf("maxx = %d\n", maxx);
				functionsComputeSameValue(maxx, new SingleFixedGaussian(maxx), new DoubleFixedGaussian(maxx), 1e-3);
			}
		}
		catch (AssertionError e)
		{
			System.out.println(e.getMessage());
			//e.printStackTrace();
			throw e;
		}
	}

	@Test(expected = java.lang.AssertionError.class)
	public void fixedSinglePrecisionIsNotMuchFasterWithGradients()
	{
		singlePrecisionIsFasterWithGradients(maxx, new SingleFixedGaussian(maxx), new DoubleFixedGaussian(maxx), false);
	}

	@Test(expected = java.lang.AssertionError.class)
	public void fixedSinglePrecisionIsNotMuchFaster()
	{
		singlePrecisionIsFaster(maxx, new SingleFixedGaussian(maxx), new DoubleFixedGaussian(maxx), false);
	}

	@Test(expected = java.lang.AssertionError.class)
	public void fixedSinglePrecisionIsNotMuchFasterWithGradientsNoSum()
	{
		singlePrecisionIsFasterWithGradients(maxx, new SingleFixedGaussian(maxx), new DoubleFixedGaussian(maxx), true);
	}

	@Test(expected = java.lang.AssertionError.class)
	public void fixedSinglePrecisionIsNotMuchFasterNoSum()
	{
		singlePrecisionIsFaster(maxx, new SingleFixedGaussian(maxx), new DoubleFixedGaussian(maxx), true);
	}

	private void functionsComputeSameValue(int maxx, SinglePrecision f1, DoublePrecision f2, final double precision)
	{
		f1.setMaxX(maxx);
		f2.setMaxX(maxx);
		float[] p1 = params1.clone();
		double[] p2 = params2.clone();
		p1[Gaussian.X_POSITION] = (float) (p2[Gaussian.X_POSITION] = (float) (0.123 + maxx / 2));
		p1[Gaussian.Y_POSITION] = (float) (p2[Gaussian.Y_POSITION] = (float) (0.789 + maxx / 2));
		f1.initialise(p1);
		f2.initialise(p2);
		final int n = p1.length;
		float[] g1 = new float[n];
		double[] g2 = new double[n];

		double t1 = 0, t2 = 0;
		double[] tg1 = new double[n];
		double[] tg2 = new double[n];

		for (int i = 0; i < maxx; i++)
		{
			float v1 = f1.eval(i);
			t1 += v1;
			double v2 = f2.eval(i);
			t2 += v2;
			Assert.assertEquals("Different values", v2, v1, precision);
			float vv1 = f1.eval(i, g1);
			double vv2 = f2.eval(i, g2);
			Assert.assertEquals("Different f1 values", v1, vv1, precision);
			Assert.assertEquals("Different f2 values", v2, vv2, precision);
			for (int j = 0; j < n; j++)
			{
				tg1[j] += g1[j];
				tg2[j] += g2[j];
			}
			Assert.assertArrayEquals("Different gradients", g2, toDouble(g1), precision);
		}
		Assert.assertArrayEquals("Different total gradients", tg2, tg1, precision);
		Assert.assertEquals("Different totals", t2, t1, precision);
	}

	private void singlePrecisionIsFasterWithGradients(int maxx, SinglePrecision f1, DoublePrecision f2, boolean noSum)
	{
		f1.setMaxX(maxx);
		f2.setMaxX(maxx);
		float[] p1 = params1.clone();
		double[] p2 = params2.clone();
		p1[Gaussian.X_POSITION] = (float) (p2[Gaussian.X_POSITION] = (float) (0.123 + maxx / 2));
		p1[Gaussian.Y_POSITION] = (float) (p2[Gaussian.Y_POSITION] = (float) (0.789 + maxx / 2));

		long time1, time2;

		if (noSum)
		{
			time1 = runSingleWithGradientsNoSum(maxx, f1, p1);
			time1 = runSingleWithGradientsNoSum(maxx, f1, p1);
			time1 += runSingleWithGradientsNoSum(maxx, f1, p1);
			time2 = runDoubleWithGradientsNoSum(maxx, f2, p2);
			time2 = runDoubleWithGradientsNoSum(maxx, f2, p2);
			time2 += runDoubleWithGradientsNoSum(maxx, f2, p2);
		}
		else
		{
			time1 = runSingleWithGradients(maxx, f1, p1);
			time1 = runSingleWithGradients(maxx, f1, p1);
			time1 += runSingleWithGradients(maxx, f1, p1);
			time2 = runDoubleWithGradients(maxx, f2, p2);
			time2 = runDoubleWithGradients(maxx, f2, p2);
			time2 += runDoubleWithGradients(maxx, f2, p2);
		}

		System.out.printf("%sGradient %s = %d, %s = %d => (%f)\n", (noSum) ? "No sum " : "", f1.getClass()
				.getSimpleName(), time1, f2.getClass().getSimpleName(), time2, (double) time2 / time1);
		Assert.assertTrue(time1 * SPEED_UP_FACTOR < time2);
	}

	@SuppressWarnings("unused")
	private long runSingleWithGradients(int maxx, SinglePrecision f, float[] p)
	{
		f.initialise(p);
		final int n = params1.length;
		float[] g = new float[n];
		double[] tg = new double[n];

		// Warm up
		for (int j = 0; j < 10; j++)
		{
			for (int i = 0; i < maxx; i++)
			{
				f.eval(i, g);
			}
		}

		long time = System.nanoTime();
		double sum = 0;
		for (int j = 0; j < MAX_ITER; j++)
		{
			sum = 0;
			for (int i = 0; i < maxx; i++)
			{
				sum += f.eval(i, g);
				for (int k = 0; k < n; k++)
					tg[k] += g[k];
			}
		}
		return System.nanoTime() - time;
	}

	private long runSingleWithGradientsNoSum(int maxx, SinglePrecision f, float[] p)
	{
		f.initialise(p);
		float[] g = new float[params1.length];

		// Warm up
		for (int j = 0; j < 10; j++)
		{
			for (int i = 0; i < maxx; i++)
			{
				f.eval(i, g);
			}
		}

		long time = System.nanoTime();
		for (int j = 0; j < MAX_ITER; j++)
		{
			for (int i = 0; i < maxx; i++)
			{
				f.eval(i, g);
			}
		}
		return System.nanoTime() - time;
	}

	@SuppressWarnings("unused")
	private long runDoubleWithGradients(int maxx, DoublePrecision f, double[] p)
	{
		f.initialise(p);
		final int n = params1.length;
		double[] g = new double[n];
		double[] tg = new double[n];
		
		// Warm up
		for (int j = 0; j < 10; j++)
		{
			for (int i = 0; i < maxx; i++)
			{
				f.eval(i, g);
			}
		}

		long time = System.nanoTime();
		double sum = 0;
		for (int j = 0; j < MAX_ITER; j++)
		{
			sum = 0;
			for (int i = 0; i < maxx; i++)
			{
				sum += f.eval(i, g);
				for (int k = 0; k < n; k++)
					tg[k] += g[k];
			}
		}
		return System.nanoTime() - time;
	}

	private long runDoubleWithGradientsNoSum(int maxx, DoublePrecision f, double[] p)
	{
		f.initialise(p);
		double[] g = new double[params1.length];

		// Warm up
		for (int j = 0; j < 10; j++)
		{
			for (int i = 0; i < maxx; i++)
			{
				f.eval(i, g);
			}
		}

		long time = System.nanoTime();
		for (int j = 0; j < MAX_ITER; j++)
		{
			for (int i = 0; i < maxx; i++)
			{
				f.eval(i, g);
			}
		}
		return System.nanoTime() - time;
	}

	private void singlePrecisionIsFaster(int maxx, SinglePrecision f1, DoublePrecision f2, boolean noSum)
	{
		f1.setMaxX(maxx);
		f2.setMaxX(maxx);
		float[] p1 = params1.clone();
		double[] p2 = params2.clone();
		p1[Gaussian.X_POSITION] = (float) (p2[Gaussian.X_POSITION] = (float) (0.123 + maxx / 2));
		p1[Gaussian.Y_POSITION] = (float) (p2[Gaussian.Y_POSITION] = (float) (0.789 + maxx / 2));

		long time1, time2;
		if (noSum)
		{
			time1 = runSingleNoSum(maxx, f1, p1);
			time1 = runSingleNoSum(maxx, f1, p1);
			time1 += runSingleNoSum(maxx, f1, p1);
			time2 = runDoubleNoSum(maxx, f2, p2);
			time2 = runDoubleNoSum(maxx, f2, p2);
			time2 += runDoubleNoSum(maxx, f2, p2);
		}
		else
		{
			time1 = runSingle(maxx, f1, p1);
			time1 = runSingle(maxx, f1, p1);
			time1 += runSingle(maxx, f1, p1);
			time2 = runDouble(maxx, f2, p2);
			time2 = runDouble(maxx, f2, p2);
			time2 += runDouble(maxx, f2, p2);
		}

		System.out.printf("%s%s = %d, %s = %d => (%f)\n", (noSum) ? "No sum " : "", f1.getClass().getSimpleName(),
				time1, f2.getClass().getSimpleName(), time2, (double) time2 / time1);
		Assert.assertTrue(time1 * SPEED_UP_FACTOR < time2);
	}

	@SuppressWarnings("unused")
	private long runSingle(int maxx, SinglePrecision f, float[] p)
	{
		// Warm up
		for (int j = 0; j < 10; j++)
		{
			f.initialise(p);
			for (int i = 0; i < maxx; i++)
			{
				f.eval(i);
			}
		}

		long time = System.nanoTime();
		double sum = 0;
		for (int j = 0; j < MAX_ITER; j++)
		{
			sum = 0;
			f.initialise(p);
			for (int i = 0; i < maxx; i++)
			{
				sum += f.eval(i);
			}
		}
		return System.nanoTime() - time;
	}

	private long runSingleNoSum(int maxx, SinglePrecision f, float[] p)
	{
		// Warm up
		for (int j = 0; j < 10; j++)
		{
			f.initialise(p);
			for (int i = 0; i < maxx; i++)
			{
				f.eval(i);
			}
		}

		long time = System.nanoTime();
		for (int j = 0; j < MAX_ITER; j++)
		{
			f.initialise(p);
			for (int i = 0; i < maxx; i++)
			{
				f.eval(i);
			}
		}
		return System.nanoTime() - time;
	}

	@SuppressWarnings("unused")
	private long runDouble(int maxx, DoublePrecision f, double[] p)
	{
		// Warm up
		for (int j = 0; j < 10; j++)
		{
			f.initialise(p);
			for (int i = 0; i < maxx; i++)
			{
				f.eval(i);
			}
		}

		long time = System.nanoTime();
		double sum = 0;
		for (int j = 0; j < MAX_ITER; j++)
		{
			sum = 0;
			f.initialise(p);
			for (int i = 0; i < maxx; i++)
			{
				sum += f.eval(i);
			}
		}
		return System.nanoTime() - time;
	}

	private long runDoubleNoSum(int maxx, DoublePrecision f, double[] p)
	{
		// Warm up
		for (int j = 0; j < 10; j++)
		{
			f.initialise(p);
			for (int i = 0; i < maxx; i++)
			{
				f.eval(i);
			}
		}

		long time = System.nanoTime();
		for (int j = 0; j < MAX_ITER; j++)
		{
			f.initialise(p);
			for (int i = 0; i < maxx; i++)
			{
				f.eval(i);
			}
		}
		return System.nanoTime() - time;
	}

	private static float[] toFloat(double[] p)
	{
		float[] f = new float[p.length];
		for (int i = 0; i < f.length; i++)
			f[i] = (float) p[i];
		return f;
	}

	private static double[] toDouble(float[] p)
	{
		double[] f = new double[p.length];
		for (int i = 0; i < f.length; i++)
			f[i] = p[i];
		return f;
	}

	void log(String format, Object... args)
	{
		System.out.printf(format, args);
	}
}
