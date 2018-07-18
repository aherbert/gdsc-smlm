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
package uk.ac.sussex.gdsc.smlm.filters;

import java.awt.Rectangle;

import org.apache.commons.math3.distribution.ExponentialDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Assert;
import org.junit.Test;

import ij.plugin.filter.GaussianBlur;
import ij.process.FloatProcessor;
import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.Random;
import uk.ac.sussex.gdsc.test.BaseTimingTask;
import uk.ac.sussex.gdsc.test.LogLevel;
import uk.ac.sussex.gdsc.test.TestLog;
import uk.ac.sussex.gdsc.test.TestSettings;
import uk.ac.sussex.gdsc.test.TimingService;

@SuppressWarnings({ "javadoc" })
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
			final FloatProcessor fp = new FloatProcessor(size, size, d);
			gf.blurGaussian(fp, sigma, sigma, GaussianFilter.DEFAULT_ACCURACY);
			return d;
		}

		@Override
		float[] filterInternal(float[] d, double sigma)
		{
			final FloatProcessor fp = new FloatProcessor(size, size, d);
			final int border = GaussianFilter.getBorder(sigma);
			final Rectangle roi = new Rectangle(border, border, size - 2 * border, size - 2 * border);
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
		final RandomGenerator rand = TestSettings.getRandomGenerator();
		final float[] data = createData(rand, size, size);
		float[] w = null;
		if (weighted)
		{
			final ExponentialDistribution ed = new ExponentialDistribution(rand, 57,
					ExponentialDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY);

			w = new float[data.length];
			for (int i = 0; i < w.length; i++)
				w[i] = (float) (1.0 / Math.max(0.01, ed.sample()));
				//w[i] = (float) (1.0 / Math.max(0.01, rand.nextGaussian() * 0.2 + 2));
				//w[i] = 0.5f;
			f1.setWeights(w);
			f2.setWeights(w);
		}

		for (final double sigma : sigmas)
		{
			final float[] e = data.clone();
			f2.run(e, sigma);
			final float[] o = data.clone();
			f1.run(o, sigma);

			double max = 0;
			for (int i = 0; i < e.length; i++)
			{
				final double d = DoubleEquality.relativeError(e[i], o[i]);
				if (max < d)
					max = d;
			}

			TestLog.info("%s vs %s w=%b @ %.1f = %g\n", f1.getName(), f2.getName(), weighted, sigma, max);
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
			final float[] d = (float[]) data;
			return filter.run(d, sigma);
		}
	}

	@Test
	public void floatFilterIsFasterThanDoubleFilter()
	{
		TestSettings.assumeSpeedTest();

		final RandomGenerator rg = TestSettings.getRandomGenerator();

		final float[][] data = new float[10][];
		for (int i = 0; i < data.length; i++)
			data[i] = createData(rg, size, size);

		final TimingService ts = new TimingService();
		for (final double sigma : sigmas)
		{
			ts.execute(new MyTimingTask(new FloatFilter(false), data, sigma));
			ts.execute(new MyTimingTask(new DPFilter(false), data, sigma));
			ts.execute(new MyTimingTask(new DoubleFilter(false), data, sigma));
		}
		final int size = ts.getSize();
		ts.repeat();
		if (TestSettings.allow(LogLevel.INFO))
			ts.report(size);
		final int n = size / sigmas.length;
		for (int i = 0, j = size; i < sigmas.length; i++, j += n)
			for (int k = 1; k < n; k++)
			{
				final double t1 = ts.get(j).getMean();
				final double t2 = ts.get(j + k).getMean();
				TestLog.logSpeedTestResult(t1 < t2, "%s %s => %s %s = %.2fx\n", ts.get(j + k).getTask().getName(),
						t2, ts.get(j).getTask().getName(), t1, t2 / t1);
			}
	}

	@Test
	public void floatFilterInternalIsFasterThanDoubleFilterInternal()
	{
		TestSettings.assumeHighComplexity();

		final RandomGenerator rg = TestSettings.getRandomGenerator();

		final float[][] data = new float[10][];
		for (int i = 0; i < data.length; i++)
			data[i] = createData(rg, size, size);

		final TimingService ts = new TimingService();
		for (final double sigma : sigmas)
		{
			ts.execute(new MyTimingTask(new FloatFilter(true), data, sigma));
			ts.execute(new MyTimingTask(new DPFilter(false), data, sigma));
			ts.execute(new MyTimingTask(new DoubleFilter(true), data, sigma));
		}
		final int size = ts.getSize();
		ts.repeat();
		if (TestSettings.allow(LogLevel.INFO))
			ts.report(size);
		final int n = size / sigmas.length;
		for (int i = 0, j = size; i < sigmas.length; i++, j += n)
			for (int k = 1; k < n; k++)
			{
				final double t1 = ts.get(j).getMean();
				final double t2 = ts.get(j + k).getMean();
				TestLog.logSpeedTestResult(t1 < t2, "%s %s => %s %s = %.2fx\n", ts.get(j + k).getTask().getName(),
						t2, ts.get(j).getTask().getName(), t1, t2 / t1);
			}
	}

	private static float[] createData(RandomGenerator rg, int width, int height)
	{
		final float[] data = new float[width * height];
		for (int i = data.length; i-- > 0;)
			data[i] = i;

		Random.shuffle(data, rg);

		return data;
	}
}
