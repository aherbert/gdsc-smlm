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
package uk.ac.sussex.gdsc.smlm.math3.distribution;

import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.Assertions;

import uk.ac.sussex.gdsc.core.utils.RandomGeneratorAdapter;
import uk.ac.sussex.gdsc.test.BaseTimingTask;
import uk.ac.sussex.gdsc.test.LogLevel;
import uk.ac.sussex.gdsc.test.TestLog;
import uk.ac.sussex.gdsc.test.TestSettings;
import uk.ac.sussex.gdsc.test.TimingService;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.junit5.SpeedTag;

@SuppressWarnings({ "javadoc" })
public class CustomPoissonDistributionTest
{
	private abstract class MyTimingTask extends BaseTimingTask
	{
		RandomSeed seed;
		UniformRandomProvider r;
		double mean;
		double min;
		int n, m = 10;

		public MyTimingTask(String name, RandomSeed seed, double min, double max)
		{
			super(String.format("%s %.1f - %.1f", name, min, max));
			this.seed = seed;
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
			r = TestSettings.getRandomGenerator(seed.getSeed());
			mean = min;
			return null;
		}
	}

	private class StaticTimingTask extends MyTimingTask
	{
		public StaticTimingTask(RandomSeed seed, double min, double max)
		{
			super("RandomDataGenerator", seed, min, max);
		}

		@Override
		public Object run(Object data)
		{
			final RandomDataGenerator rdg = new RandomDataGenerator(new RandomGeneratorAdapter(r));
			final long[] e = new long[n * m];
			for (int i = 0, k = 0; i < n; i++)
			{
				for (int j = 0; j < m; j++, k++)
					e[k] = rdg.nextPoisson(mean);
				mean += 1;
			}
			return e;
		}
	}

	private class InstanceTimingTask extends MyTimingTask
	{
		public InstanceTimingTask(RandomSeed seed, double min, double max)
		{
			super("Instance", seed, min, max);
		}

		@Override
		public Object run(Object data)
		{
			final CustomPoissonDistribution dist = new CustomPoissonDistribution(new RandomGeneratorAdapter(r), 1);
			final long[] e = new long[n * m];
			for (int i = 0, k = 0; i < n; i++)
			{
				dist.setMean(mean);
				for (int j = 0; j < m; j++, k++)
					e[k] = dist.sample();
				mean += 1;
			}
			return e;
		}
	}

	@SeededTest
	public void canCreateSamples(RandomSeed seed)
	{
		final StaticTimingTask t1 = new StaticTimingTask(seed, 0.5, 60);
		t1.getData(0);
		final long[] e = (long[]) t1.run(null);

		final InstanceTimingTask t2 = new InstanceTimingTask(seed, 0.5, 60);
		t2.getData(0);
		final long[] o = (long[]) t2.run(null);

		Assertions.assertArrayEquals(e, o);
	}

	@SpeedTag
	@SeededTest
	public void customDistributionIsFasterWithTinyMean(RandomSeed seed)
	{
		final TimingService ts = new TimingService(5);
		ts.execute(new StaticTimingTask(seed, 0.5, 10));
		ts.execute(new InstanceTimingTask(seed, 0.5, 10));

		final int size = ts.getSize();
		ts.repeat(size);
		if (TestSettings.allow(LogLevel.INFO))
			ts.report(size);

		//Assertions.assertTrue(ts.get(-1).getMean() < ts.get(-2).getMean());
		final double t1 = ts.get(-1).getMean();
		final double t2 = ts.get(-2).getMean();
		TestLog.logSpeedTestResult(t1 < t2, "RandomDataGenerator  %s  vs CustomPoissonDistribution  %s : %.2f", t2, t1,
				t2 / t1);
	}

	@SpeedTag
	@SeededTest
	public void customDistributionIsFasterWithSmallMean(RandomSeed seed)
	{
		final TimingService ts = new TimingService(5);
		ts.execute(new StaticTimingTask(seed, 10, 38));
		ts.execute(new InstanceTimingTask(seed, 10, 38));

		final int size = ts.getSize();
		ts.repeat(size);
		if (TestSettings.allow(LogLevel.INFO))
			ts.report(size);

		//Assertions.assertTrue(ts.get(-1).getMean() < ts.get(-2).getMean());
		final double t1 = ts.get(-1).getMean();
		final double t2 = ts.get(-2).getMean();
		TestLog.logSpeedTestResult(t1 < t2, "RandomDataGenerator  %s  vs CustomPoissonDistribution  %s : %.2f", t2, t1,
				t2 / t1);
	}

	@SpeedTag
	@SeededTest
	public void customDistributionIsFasterWithBigMean(RandomSeed seed)
	{
		// When the mean is above 40 the PoissonDistribution switches to a different
		// sampling method and this is so slow that the speed increase from using
		// the instance class is negligible. However test it is still faster. If this fails
		// then Apache commons may have changed their implementation and the custom
		// class should be updated.

		final TimingService ts = new TimingService(5);
		ts.execute(new StaticTimingTask(seed, 40.5, 60));
		ts.execute(new InstanceTimingTask(seed, 40.5, 60));

		final int size = ts.getSize();
		ts.repeat(size);
		if (TestSettings.allow(LogLevel.INFO))
			ts.report(size);

		//Assertions.assertTrue(ts.get(-1).getMean() < ts.get(-2).getMean());
		final double t1 = ts.get(-1).getMean();
		final double t2 = ts.get(-2).getMean();
		TestLog.logSpeedTestResult(t1 < t2, "RandomDataGenerator  %s  vs CustomPoissonDistribution  %s : %.2f", t2, t1,
				t2 / t1);
	}
}
