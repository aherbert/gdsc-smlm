package org.apache.commons.math3.distribution;

import org.apache.commons.math3.distribution.CustomGammaDistribution;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.test.BaseTimingTask;
import gdsc.core.test.TimingService;

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
			r = new Well19937c();
		}

		public int getSize()
		{
			return 1;
		}

		public Object getData(int i)
		{
			r.setSeed(30051977);
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
		TimingService ts = new TimingService(5);
		ts.execute(new StaticTimingTask());
		ts.execute(new InstanceTimingTask());

		int size = ts.getSize();
		ts.repeat(size);
		ts.report(size);
		
		Assert.assertTrue(ts.get(-1).getMean() < ts.get(-2).getMean());
	}
}
