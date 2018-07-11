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
package org.apache.commons.math3.distribution;

import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Assert;
import org.junit.Test;

import gdsc.test.BaseTimingTask;
import gdsc.test.TestSettings;
import gdsc.test.TestSettings.LogLevel;
import gdsc.test.TimingService;

@SuppressWarnings({ "javadoc" })
public class CustomPoissonDistributionTest
{
	private abstract class MyTimingTask extends BaseTimingTask
	{
		RandomGenerator r;
		double mean;
		double min;
		int n, m = 10;

		public MyTimingTask(String name, double min, double max)
		{
			super(String.format("%s %.1f - %.1f", name, min, max));
			r = TestSettings.getRandomGenerator();
			this.min = min;
			mean = min;
			n = 0;
			while (mean < max)
			{
				n++;
				mean += 1;
			}
		}

		@Override
		public int getSize()
		{
			return 1;
		}

		@Override
		public Object getData(int i)
		{
			r.setSeed(TestSettings.getSeed());
			mean = min;
			return null;
		}
	}

	private class StaticTimingTask extends MyTimingTask
	{
		RandomDataGenerator rdg;

		public StaticTimingTask(double min, double max)
		{
			super("RandomDataGenerator", min, max);
			rdg = new RandomDataGenerator(r);
		}

		@Override
		public Object run(Object data)
		{
			long[] e = new long[n * m];
			for (int i = 0, k = 0; i < n; i++)
			{
				for (int j = 0; j < m; j++, k++)
				{
					e[k] = rdg.nextPoisson(mean);
				}
				mean += 1;
			}
			return e;
		}
	}

	private class InstanceTimingTask extends MyTimingTask
	{
		CustomPoissonDistribution dist;

		public InstanceTimingTask(double min, double max)
		{
			super("Instance", min, max);
			dist = new CustomPoissonDistribution(r, 1);
		}

		@Override
		public Object run(Object data)
		{
			long[] e = new long[n * m];
			for (int i = 0, k = 0; i < n; i++)
			{
				dist.setMean(mean);
				for (int j = 0; j < m; j++, k++)
				{
					e[k] = dist.sample();
				}
				mean += 1;
			}
			return e;
		}
	}

	@Test
	public void canCreateSamples()
	{
		StaticTimingTask t1 = new StaticTimingTask(0.5, 60);
		t1.getData(0);
		long[] e = (long[]) t1.run(null);

		InstanceTimingTask t2 = new InstanceTimingTask(0.5, 60);
		t2.getData(0);
		long[] o = (long[]) t2.run(null);

		Assert.assertArrayEquals(e, o);
	}

	@Test
	public void customDistributionIsFasterWithTinyMean()
	{
		TimingService ts = new TimingService(5);
		ts.execute(new StaticTimingTask(0.5, 10));
		ts.execute(new InstanceTimingTask(0.5, 10));

		int size = ts.getSize();
		ts.repeat(size);
		if (TestSettings.allow(LogLevel.INFO))
			ts.report(size);

		//Assert.assertTrue(ts.get(-1).getMean() < ts.get(-2).getMean());
		double t1 = ts.get(-1).getMean();
		double t2 = ts.get(-2).getMean();
		TestSettings.logSpeedTestResult(t1 < t2, "RandomDataGenerator  %s  vs CustomPoissonDistribution  %s : %.2f", t2,
				t1, t2 / t1);
	}

	@Test
	public void customDistributionIsFasterWithSmallMean()
	{
		TimingService ts = new TimingService(5);
		ts.execute(new StaticTimingTask(10, 38));
		ts.execute(new InstanceTimingTask(10, 38));

		int size = ts.getSize();
		ts.repeat(size);
		if (TestSettings.allow(LogLevel.INFO))
			ts.report(size);

		//Assert.assertTrue(ts.get(-1).getMean() < ts.get(-2).getMean());
		double t1 = ts.get(-1).getMean();
		double t2 = ts.get(-2).getMean();
		TestSettings.logSpeedTestResult(t1 < t2, "RandomDataGenerator  %s  vs CustomPoissonDistribution  %s : %.2f", t2,
				t1, t2 / t1);
	}

	@Test
	public void customDistributionIsFasterWithBigMean()
	{
		// When the mean is above 40 the PoissonDistribution switches to a different
		// sampling method and this is so slow that the speed increase from using
		// the instance class is negligible. However test it is still faster. If this fails
		// then Apache commons may have changed their implementation and the custom
		// class should be updated.

		TimingService ts = new TimingService(5);
		ts.execute(new StaticTimingTask(40.5, 60));
		ts.execute(new InstanceTimingTask(40.5, 60));

		int size = ts.getSize();
		ts.repeat(size);
		if (TestSettings.allow(LogLevel.INFO))
			ts.report(size);

		//Assert.assertTrue(ts.get(-1).getMean() < ts.get(-2).getMean());
		double t1 = ts.get(-1).getMean();
		double t2 = ts.get(-2).getMean();
		TestSettings.logSpeedTestResult(t1 < t2, "RandomDataGenerator  %s  vs CustomPoissonDistribution  %s : %.2f", t2,
				t1, t2 / t1);
	}
}
