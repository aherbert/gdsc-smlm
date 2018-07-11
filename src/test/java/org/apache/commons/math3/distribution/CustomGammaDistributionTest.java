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
public class CustomGammaDistributionTest
{
	private abstract class MyTimingTask extends BaseTimingTask
	{
		RandomGenerator r;
		double shape = 0.5;
		double scale = 300;
		int n = 1000, m = 10;

		public MyTimingTask(String name)
		{
			super(name);
			r = TestSettings.getRandomGenerator();
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
			shape = 0.5;
			return null;
		}
	}

	private class StaticTimingTask extends MyTimingTask
	{
		RandomDataGenerator rdg;

		public StaticTimingTask()
		{
			super("RandomDataGenerator");
			rdg = new RandomDataGenerator(r);
		}

		@Override
		public Object run(Object data)
		{
			double[] e = new double[n * m];
			for (int i = 0, k = 0; i < n; i++)
			{
				for (int j = 0; j < m; j++, k++)
				{
					e[k] = rdg.nextGamma(shape, scale);
				}
				shape += 1;
			}
			return e;
		}
	}

	private class InstanceTimingTask extends MyTimingTask
	{
		CustomGammaDistribution dist;

		public InstanceTimingTask()
		{
			super("Instance");
			dist = new CustomGammaDistribution(r, 1, scale);
		}

		@Override
		public Object run(Object data)
		{
			double[] e = new double[n * m];
			for (int i = 0, k = 0; i < n; i++)
			{
				dist.setShape(shape);
				for (int j = 0; j < m; j++, k++)
				{
					e[k] = dist.sample();
				}
				shape += 1;
			}
			return e;
		}
	}

	@Test
	public void canCreateSamples()
	{
		StaticTimingTask t1 = new StaticTimingTask();
		t1.getData(0);
		double[] e = (double[]) t1.run(null);

		InstanceTimingTask t2 = new InstanceTimingTask();
		t2.getData(0);
		double[] o = (double[]) t2.run(null);

		Assert.assertArrayEquals(e, o, 0);
	}

	@Test
	public void customDistributionIsFaster()
	{
		TestSettings.assumeMediumComplexity();

		TimingService ts = new TimingService(5);
		ts.execute(new StaticTimingTask());
		ts.execute(new InstanceTimingTask());

		int size = ts.getSize();
		ts.repeat(size);
		if (TestSettings.allow(LogLevel.INFO))
			ts.report(size);

		Assert.assertTrue(ts.get(-1).getMean() < ts.get(-2).getMean());
	}
}
