package gdsc.smlm.filters;

import java.awt.Rectangle;

import org.apache.commons.math3.distribution.ExponentialDistribution;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.test.BaseTimingTask;
import gdsc.core.test.TimingService;
import gdsc.core.utils.DoubleEquality;
import gdsc.core.utils.Random;
import ij.plugin.filter.GaussianBlur;
import ij.process.FloatProcessor;

public class GaussianFilterTest
{
	double[] sigmas = new double[] { 12.4, 9.3, 5, 3.2, 2.1, 0.5 };
	int size = 256;

	private abstract class GFilter
	{
		boolean internal;
		String name;

		GFilter(String name, boolean internal)
		{
			this.name = name;
			this.internal = internal;
		}

		String getName()
		{
			if (internal)
				return name + " internal";
			return name;
		}

		float[] run(float[] d, double sigma)
		{
			if (internal)
				return filterInternal(d, sigma);
			return filter(d, sigma);
		}

		abstract float[] filter(float[] d, double sigma);

		abstract float[] filterInternal(float[] d, double sigma);

		abstract void setWeights(float[] w);
	}

	private class IJFilter extends GFilter
	{
		GaussianBlur gf = new GaussianBlur();
		
		IJFilter(boolean internal)
		{
			super(GaussianBlur.class.getSimpleName(), internal);
		}

		@Override
		float[] filter(float[] d, double sigma)
		{
			FloatProcessor fp = new FloatProcessor(size, size, d);
			gf.blurGaussian(fp, sigma, sigma, GaussianFilter.DEFAULT_ACCURACY);
			return d;
		}

		@Override
		float[] filterInternal(float[] d, double sigma)
		{
			FloatProcessor fp = new FloatProcessor(size, size, d);
			final int border = GaussianFilter.getBorder(sigma);
			Rectangle roi = new Rectangle(border, border, size - 2 * border, size - 2 * border);
			fp.setRoi(roi);
			gf.blurGaussian(fp, sigma, sigma, GaussianFilter.DEFAULT_ACCURACY);
			return d;
		}

		@Override
		void setWeights(float[] w)
		{
			// Ignored
		}
	}	
	
	private class FloatFilter extends GFilter
	{
		GaussianFilter gf = new GaussianFilter();

		FloatFilter(boolean internal)
		{
			super(GaussianFilter.class.getSimpleName(), internal);
		}

		@Override
		float[] filter(float[] d, double sigma)
		{
			gf.convolve(d, size, size, sigma);
			return d;
		}

		@Override
		float[] filterInternal(float[] d, double sigma)
		{
			gf.convolveInternal(d, size, size, sigma);
			return d;
		}

		@Override
		void setWeights(float[] w)
		{
			gf.setWeights(w, size, size);
		}
	}

	private class DoubleFilter extends GFilter
	{
		DoubleGaussianFilter gf = new DoubleGaussianFilter();

		DoubleFilter(boolean internal)
		{
			super(DoubleGaussianFilter.class.getSimpleName(), internal);
		}

		@Override
		float[] filter(float[] d, double sigma)
		{
			gf.convolve(d, size, size, sigma);
			return d;
		}

		@Override
		float[] filterInternal(float[] d, double sigma)
		{
			gf.convolveInternal(d, size, size, sigma);
			return d;
		}

		@Override
		void setWeights(float[] w)
		{
			gf.setWeights(w, size, size);
		}
	}

	private class DPFilter extends GFilter
	{
		DPGaussianFilter gf = new DPGaussianFilter();

		DPFilter(boolean internal)
		{
			super(DPGaussianFilter.class.getSimpleName(), internal);
		}

		@Override
		float[] filter(float[] d, double sigma)
		{
			gf.convolve(d, size, size, sigma);
			return d;
		}

		@Override
		float[] filterInternal(float[] d, double sigma)
		{
			gf.convolveInternal(d, size, size, sigma);
			return d;
		}

		@Override
		void setWeights(float[] w)
		{
			gf.setWeights(w, size, size);
		}
	}

	@Test
	public void floatFilterIsSameAsIJFilter()
	{
		filter1IsSameAsFilter2(new FloatFilter(false), new IJFilter(false), false, 1e-2);
	}

	@Test
	public void floatFilterInternalIsSameAsIJFilter()
	{
		filter1IsSameAsFilter2(new FloatFilter(true), new IJFilter(true), false, 1e-2);
	}
	
	@Test
	public void floatFilterIsSameAsDoubleFilter()
	{
		filter1IsSameAsFilter2(new FloatFilter(false), new DoubleFilter(false), false, 1e-2);
	}

	@Test
	public void floatFilterIsSameAsDoubleFilterWeighted()
	{
		filter1IsSameAsFilter2(new FloatFilter(false), new DoubleFilter(false), true, 1e-2);
	}

	@Test
	public void dpFloatFilterIsSameAsDoubleFilter()
	{
		filter1IsSameAsFilter2(new DPFilter(false), new DoubleFilter(false), false, 1e-2);
	}

	@Test
	public void dpFloatFilterIsSameAsDoubleFilterWeighted()
	{
		filter1IsSameAsFilter2(new DPFilter(false), new DoubleFilter(false), true, 1e-2);
	}

	private void filter1IsSameAsFilter2(GFilter f1, GFilter f2, boolean weighted, double tolerance)
	{
		Random rand = new Random(-30051976);
		float[] data = createData(rand, size, size);
		float[] w = null;
		if (weighted)
		{
			ExponentialDistribution ed = new ExponentialDistribution(rand, 57,
					ExponentialDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY);

			w = new float[data.length];
			for (int i = 0; i < w.length; i++)
			{
				w[i] = (float) (1.0 / Math.max(0.01, ed.sample()));
				//w[i] = (float) (1.0 / Math.max(0.01, rand.nextGaussian() * 0.2 + 2));
				//w[i] = 0.5f;
			}
			f1.setWeights(w);
			f2.setWeights(w);
		}

		for (double sigma : sigmas)
		{
			float[] e = data.clone();
			f2.run(e, sigma);
			float[] o = data.clone();
			f1.run(o, sigma);

			double max = 0;
			for (int i = 0; i < e.length; i++)
			{
				double d = DoubleEquality.relativeError(e[i], o[i]);
				if (max < d)
					max = d;
			}

			System.out.printf("%s vs %s w=%b @ %.1f = %g\n", f1.getName(), f2.getName(), weighted, sigma, max);
			Assert.assertTrue(max < tolerance);
		}
	}

	private class MyTimingTask extends BaseTimingTask
	{
		GFilter filter;
		float[][] data;
		double sigma;

		public MyTimingTask(GFilter filter, float[][] data, double sigma)
		{
			super(filter.getName() + " " + sigma);
			this.filter = filter;
			this.data = data;
			this.sigma = sigma;
		}

		public int getSize()
		{
			return data.length;
		}

		public Object getData(int i)
		{
			return data[i].clone();
		}

		public Object run(Object data)
		{
			float[] d = (float[]) data;
			return filter.run(d, sigma);
		}
	}

	@Test
	public void floatFilterIsFasterThanDoubleFilter()
	{
		Random rand = new Random(-30051976);

		float[][] data = new float[10][];
		for (int i = 0; i < data.length; i++)
			data[i] = createData(rand, size, size);

		TimingService ts = new TimingService();
		for (double sigma : sigmas)
		{
			ts.execute(new MyTimingTask(new FloatFilter(false), data, sigma));
			ts.execute(new MyTimingTask(new DPFilter(false), data, sigma));
			ts.execute(new MyTimingTask(new DoubleFilter(false), data, sigma));
		}
		int size = ts.getSize();
		ts.repeat();
		ts.report(size);
		int n = size / sigmas.length;
		for (int i = 0, j = size; i < sigmas.length; i++, j += n)
		{
			for (int k = 1; k < n; k++)
				Assert.assertTrue(ts.get(j).getMean() < ts.get(j + k).getMean() * 1.1);
		}
	}

	//@Test
	public void floatFilterInternalIsFasterThanDoubleFilterInternal()
	{
		Random rand = new Random(-30051976);

		float[][] data = new float[10][];
		for (int i = 0; i < data.length; i++)
			data[i] = createData(rand, size, size);

		TimingService ts = new TimingService();
		for (double sigma : sigmas)
		{
			ts.execute(new MyTimingTask(new FloatFilter(true), data, sigma));
			ts.execute(new MyTimingTask(new DPFilter(false), data, sigma));
			ts.execute(new MyTimingTask(new DoubleFilter(true), data, sigma));
		}
		int size = ts.getSize();
		ts.repeat();
		ts.report(size);
		int n = size / sigmas.length;
		for (int i = 0, j = size; i < sigmas.length; i++, j += n)
		{
			for (int k = 1; k < n; k++)
				Assert.assertTrue(ts.get(j).getMean() < ts.get(j + k).getMean() * 1.1);
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
