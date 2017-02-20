package gdsc.smlm.ij.frc;

import org.apache.commons.math3.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.test.BaseTimingTask;
import gdsc.core.test.TimingService;
import gdsc.core.utils.DoubleEquality;

public class FRCTest
{
	@Test
	public void canComputeSine()
	{
		int steps = 1000;
		double delta = 2 * Math.PI / steps;
		for (int i = 0; i <= steps; i++)
		{
			double a = i * delta;
			double cosA = Math.cos(a);
			double e = Math.sin(a);
			double o = FRC.getSine(a, cosA);
			//System.out.printf("%f  %f ?= %f\n", a, e, o);
			Assert.assertTrue(DoubleEquality.almostEqualRelativeOrAbsolute(o, e, 1e-6, 1e-10));
		}
	}

	private abstract class MyTimingTask extends BaseTimingTask
	{
		public MyTimingTask(String name)
		{
			super(name);
		}

		public int getSize()
		{
			return 1;
		}

		public Object getData(int i)
		{
			return null;
		}
	}

	@Test
	public void computeSineIsFaster()
	{
		int steps = 100000;
		double delta = 2 * Math.PI / steps;
		final double[] a = new double[steps + 1];
		final double[] cosA = new double[steps + 1];
		for (int i = 0; i <= steps; i++)
		{
			a[i] = i * delta;
			cosA[i] = Math.cos(a[i]);
		}

		TimingService ts = new TimingService(100);
		ts.execute(new MyTimingTask("sin")
		{
			public Object run(Object data)
			{
				double d = 0;
				for (int i = 0; i < a.length; i++)
					d += Math.sin(a[i]);
				return d;
			}
		});
		ts.execute(new MyTimingTask("FastMath.sin")
		{
			public Object run(Object data)
			{
				double d = 0;
				for (int i = 0; i < a.length; i++)
					d += FastMath.sin(a[i]);
				return d;
			}
		});
		ts.execute(new MyTimingTask("getSine")
		{
			public Object run(Object data)
			{
				double d = 0;
				for (int i = 0; i < a.length; i++)
					d += FRC.getSine(a[i], cosA[i]);
				return d;
			}
		});

		ts.repeat(ts.getSize());
		ts.report();
	}
}
