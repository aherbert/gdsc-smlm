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
package gdsc.smlm.filters;

import java.awt.Rectangle;

import org.junit.Assert;
import org.junit.Test;

import gdsc.core.test.BaseTimingTask;
import gdsc.core.test.TimingService;
import gdsc.core.utils.DoubleEquality;
import gdsc.core.utils.Maths;
import gdsc.core.utils.Random;
import ij.plugin.filter.Convolver;
import ij.process.FloatProcessor;

public class KernelFilterTest
{
	int size = 256;
	int[] borders = { 0, 1, 2, 3, 5, 10 };

	static float[] createKernel(int kw, int kh)
	{
		// Simple linear ramp
		float[] k = new float[kw * kh];
		int cx = kw / 2;
		int cy = kh / 2;
		for (int y = 0, i = 0; y < kh; y++)
		{
			int dy2 = Maths.pow2(cy - y);
			for (int x = 0; x < kw; x++)
			{
				int dx2 = Maths.pow2(cx - x);
				k[i++] = (float) Math.sqrt(dx2 + dy2);
			}
		}
		// Invert
		float max = k[0];
		for (int i = 0; i < k.length; i++)
			k[i] = max - k[i];
		return k;
	}

	private abstract class FilterWrapper
	{
		final float[] kernel;
		final int kw, kh;
		String name;

		FilterWrapper(String name, float[] kernel, int kw, int kh)
		{
			this.name = name + " " + kw + "x" + kh;
			this.kernel = kernel;
			this.kw = kw;
			this.kh = kh;
		}

		String getName()
		{
			return name;
		}

		abstract float[] filter(float[] d, int border);

		abstract void setWeights(float[] w);
	}

	private class ConvolverWrapper extends FilterWrapper
	{
		Convolver kf = new Convolver();

		ConvolverWrapper(float[] kernel, int kw, int kh)
		{
			super(Convolver.class.getSimpleName(), kernel, kw, kh);
		}

		@Override
		float[] filter(float[] d, int border)
		{
			FloatProcessor fp = new FloatProcessor(size, size, d);
			if (border > 0)
			{
				Rectangle roi = new Rectangle(border, border, size - 2 * border, size - 2 * border);
				fp.setRoi(roi);
			}
			kf.convolveFloat(fp, kernel, kw, kh);
			return d;
		}

		@Override
		void setWeights(float[] w)
		{
			// Ignored
		}
	}

	private class KernelFilterWrapper extends FilterWrapper
	{
		KernelFilter kf = new KernelFilter(kernel, kw, kh);

		KernelFilterWrapper(float[] kernel, int kw, int kh)
		{
			super(KernelFilterTest.class.getSimpleName(), kernel, kw, kh);
		}

		@Override
		float[] filter(float[] d, int border)
		{
			kf.convolve(d, size, size, border);
			return d;
		}

		@Override
		void setWeights(float[] w)
		{
			kf.setWeights(w, size, size);
		}
	}

	private class ZeroKernelFilterWrapper extends FilterWrapper
	{
		ZeroKernelFilter kf = new ZeroKernelFilter(kernel, kw, kh);

		ZeroKernelFilterWrapper(float[] kernel, int kw, int kh)
		{
			super(ZeroKernelFilterWrapper.class.getSimpleName(), kernel, kw, kh);
		}

		@Override
		float[] filter(float[] d, int border)
		{
			kf.convolve(d, size, size, border);
			return d;
		}

		@Override
		void setWeights(float[] w)
		{
			kf.setWeights(w, size, size);
		}
	}

	@Test
	public void canRotate180()
	{
		for (int kw = 1; kw < 3; kw++)
			for (int kh = 1; kh < 3; kh++)
			{
				float[] kernel = createKernel(kw, kh);
				FloatProcessor fp = new FloatProcessor(kw, kh, kernel.clone());
				fp.flipHorizontal();
				fp.flipVertical();
				KernelFilter.rotate180(kernel);
				Assert.assertArrayEquals((float[]) fp.getPixels(), kernel, 0);
			}
	}

	@Test
	public void kernelFilterIsSameAsIJFilter()
	{
		int kw = 5, kh = 5;
		float[] kernel = createKernel(kw, kh);
		filter1IsSameAsFilter2(new KernelFilterWrapper(kernel, kw, kh), new ConvolverWrapper(kernel, kw, kh), false,
				1e-2);
	}

	@Test
	public void zeroKernelFilterIsSameAsIJFilter()
	{
		int kw = 5, kh = 5;
		float[] kernel = createKernel(kw, kh);
		filter1IsSameAsFilter2(new ZeroKernelFilterWrapper(kernel, kw, kh), new ConvolverWrapper(kernel, kw, kh), true,
				1e-2);
	}

	private void filter1IsSameAsFilter2(FilterWrapper f1, FilterWrapper f2, boolean internal, double tolerance)
	{
		Random rand = new Random(-30051976);
		float[] data = createData(rand, size, size);

		int testBorder = (internal) ? f1.kw / 2 : 0;
		for (int border : borders)
		{
			filter1IsSameAsFilter2(f1, f2, data, border, testBorder, tolerance);
		}
	}

	private void filter1IsSameAsFilter2(FilterWrapper f1, FilterWrapper f2, float[] data, int border, int testBorder,
			double tolerance)
	{
		float[] e = data.clone();
		f2.filter(e, border);
		float[] o = data.clone();
		f1.filter(o, border);

		double max = 0;
		if (testBorder == 0)
		{
			for (int i = 0; i < e.length; i++)
			{
				double d = DoubleEquality.relativeError(e[i], o[i]);
				if (max < d)
					max = d;
			}
		}
		else
		{
			int limit = size - testBorder;
			for (int y = testBorder; y < limit; y++)
			{
				for (int x = testBorder, i = y * size + x; x < limit; x++, i++)
				{
					double d = DoubleEquality.relativeError(e[i], o[i]);
					if (max < d)
						max = d;
				}
			}
		}

		System.out.printf("%s vs %s @ %d = %g\n", f1.getName(), f2.getName(), border, max);
		Assert.assertTrue(max < tolerance);
	}

	private class MyTimingTask extends BaseTimingTask
	{
		FilterWrapper filter;
		float[][] data;
		int border;

		public MyTimingTask(FilterWrapper filter, float[][] data, int border)
		{
			super(filter.getName() + " " + border);
			this.filter = filter;
			this.data = data;
			this.border = border;
		}

		@Override
		public int getSize()
		{
			return data.length;
		}

		@Override
		public Object getData(int i)
		{
			return data[i].clone();
		}

		@Override
		public Object run(Object data)
		{
			float[] d = (float[]) data;
			return filter.filter(d, border);
		}
	}

	@Test
	public void floatFilterIsFasterThanIJFilter()
	{
		floatFilterIsFasterThanIJFilter(5);
		floatFilterIsFasterThanIJFilter(11);
	}

	private void floatFilterIsFasterThanIJFilter(int k)
	{
		Random rand = new Random(-30051976);

		float[][] data = new float[10][];
		for (int i = 0; i < data.length; i++)
			data[i] = createData(rand, size, size);

		float[] kernel = createKernel(k, k);
		for (int border : borders)
		{
			TimingService ts = new TimingService();
			ts.execute(new MyTimingTask(new KernelFilterWrapper(kernel, k, k), data, border));
			ts.execute(new MyTimingTask(new ZeroKernelFilterWrapper(kernel, k, k), data, border));
			ts.execute(new MyTimingTask(new ConvolverWrapper(kernel, k, k), data, border));
			int size = ts.getSize();
			ts.repeat();
			ts.report(size);
			Assert.assertTrue(ts.get(-2).getMean() < ts.get(-1).getMean() * 1.1);
		}
	}

	private static float[] createData(Random rand, int width, int height)
	{
		float[] data = new float[width * height];
		for (int i = data.length; i-- > 0;)
			data[i] = i;

		rand.shuffle(data);

		return data;
	}
}
