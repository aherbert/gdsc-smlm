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
import uk.ac.sussex.gdsc.test.TestSettings;
import uk.ac.sussex.gdsc.test.TimingService;
import uk.ac.sussex.gdsc.test.junit5.ExtraAssumptions;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.junit5.SpeedTag;

@SuppressWarnings({ "javadoc" })
public class CustomGammaDistributionTest
{
	private abstract class MyTimingTask extends BaseTimingTask
	{
		RandomSeed seed;
		UniformRandomProvider r;
		double shape = 0.5;
		double scale = 300;
		int n = 1000, m = 10;

		public MyTimingTask(String name, RandomSeed seed)
		{
			super(name);
			this.seed = seed;
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
			shape = 0.5;
			return null;
		}
	}

	private class StaticTimingTask extends MyTimingTask
	{
		public StaticTimingTask(RandomSeed seed)
		{
			super("RandomDataGenerator", seed);
		}

		@Override
		public Object run(Object data)
		{
			RandomDataGenerator rdg = new RandomDataGenerator(new RandomGeneratorAdapter(r));
			final double[] e = new double[n * m];
			for (int i = 0, k = 0; i < n; i++)
			{
				for (int j = 0; j < m; j++, k++)
					e[k] = rdg.nextGamma(shape, scale);
				shape += 1;
			}
			return e;
		}
	}

	private class InstanceTimingTask extends MyTimingTask
	{
		public InstanceTimingTask(RandomSeed seed)
		{
			super("Instance", seed);
		}

		@Override
		public Object run(Object data)
		{
			CustomGammaDistribution dist = new CustomGammaDistribution(new RandomGeneratorAdapter(r), 1, scale);
			final double[] e = new double[n * m];
			for (int i = 0, k = 0; i < n; i++)
			{
				dist.setShape(shape);
				for (int j = 0; j < m; j++, k++)
					e[k] = dist.sample();
				shape += 1;
			}
			return e;
		}
	}

	@SeededTest
	public void canCreateSamples(RandomSeed seed)
	{
		final StaticTimingTask t1 = new StaticTimingTask(seed);
		t1.getData(0);
		final double[] e = (double[]) t1.run(null);

		final InstanceTimingTask t2 = new InstanceTimingTask(seed);
		t2.getData(0);
		final double[] o = (double[]) t2.run(null);

		Assertions.assertArrayEquals(e, o);
	}

	@SpeedTag
	@SeededTest
	public void customDistributionIsFaster(RandomSeed seed)
	{
		ExtraAssumptions.assumeMediumComplexity();

		final TimingService ts = new TimingService(5);
		ts.execute(new StaticTimingTask(seed));
		ts.execute(new InstanceTimingTask(seed));

		final int size = ts.getSize();
		ts.repeat(size);
		if (TestSettings.allow(LogLevel.INFO))
			ts.report(size);

		Assertions.assertTrue(ts.get(-1).getMean() < ts.get(-2).getMean());
	}
}
