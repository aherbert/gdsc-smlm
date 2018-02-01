package gdsc.smlm.function;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.test.BaseTimingTask;
import gdsc.core.test.TimingService;
import gdsc.core.utils.DoubleEquality;

@SuppressWarnings("unused")
public class FastLogTest
{
	FFastLog fLog = new FFastLog();
	DFastLog dLog = new DFastLog();

	//@formatter:off
	private abstract class BaseLog
	{
		String name;
		BaseLog() { this.name = this.getClass().getSimpleName(); }
		abstract float log(float x);
		abstract double log(double x);
		abstract boolean isDoublePrecision();
	}
	private class Math_log extends BaseLog
	{
		float log(float x) { return (float) Math.log(x); }
		double log(double x) { return Math.log(x); }
		boolean isDoublePrecision() { return true; }
	}
	private class FFastLog_log extends BaseLog
	{
		FFastLog f;
		FFastLog_log(FFastLog f) { this.f=f; }
		FFastLog_log() { this(fLog); }
		float log(float x) { return f.log(x); }
		double log(double x) { return f.log(x); }
		boolean isDoublePrecision() { return false; }
	}
	private class FFastLog_fastLog extends BaseLog
	{
		FFastLog f;
		FFastLog_fastLog(FFastLog f) { this.f=f; }
		FFastLog_fastLog() { this(fLog); }
		float log(float x) { return f.fastLog(x); }
		double log(double x) { return f.fastLog(x); }
		boolean isDoublePrecision() { return false; }
	}
	private class FFastLog_log2 extends BaseLog
	{
		FFastLog f;
		FFastLog_log2(FFastLog f) { this.f=f; }
		float log(float x) { return f.log2(x); }
		double log(double x) { return f.log2(x); }
		boolean isDoublePrecision() { return false; }
	}
	private class FFastLog_fastLog2 extends BaseLog
	{
		FFastLog f;
		FFastLog_fastLog2(FFastLog f) { this.f=f; }
		float log(float x) { return f.fastLog2(x); }
		double log(double x) { return f.fastLog2(x); }
		boolean isDoublePrecision() { return false; }
	}
	private class DFastLog_log extends BaseLog
	{
		DFastLog f;
		DFastLog_log(DFastLog f) { this.f=f; }
		DFastLog_log() { this(dLog); }
		float log(float x) { return (float) f.log(x); }
		double log(double x) { return f.log(x); }
		boolean isDoublePrecision() { return true; }
	}
	private class DFastLog_fastLog extends BaseLog
	{
		DFastLog f;
		DFastLog_fastLog(DFastLog f) { this.f=f; }
		DFastLog_fastLog() { this(dLog); }
		float log(float x) { return (float) f.fastLog(x); }
		double log(double x) { return f.fastLog(x); }
		boolean isDoublePrecision() { return true; }
	}
	private class DFastLog_log2 extends BaseLog
	{
		DFastLog f;
		DFastLog_log2(DFastLog f) { this.f=f; }
		float log(float x) { return (float) f.log2(x); }
		double log(double x) { return f.log2(x); }
		boolean isDoublePrecision() { return true; }
	}
	private class DFastLog_fastLog2 extends BaseLog
	{
		DFastLog f;
		DFastLog_fastLog2(DFastLog f) { this.f=f; }
		float log(float x) { return (float) f.fastLog2(x); }
		double log(double x) { return f.fastLog2(x); }
		boolean isDoublePrecision() { return true; }
	}	//@formatter:on

	@Test
	public void canComputeFastLog_log()
	{
		canComputeLog(new FFastLog_log(), true);
	}

	@Test
	public void canComputeFastLog_fastLog()
	{
		canComputeLog(new FFastLog_fastLog(), false);
	}

	@Test
	public void canComputeDFastLog_log()
	{
		canComputeLog(new DFastLog_log(), true);
	}

	@Test
	public void canComputeDFastLog_fastLog()
	{
		canComputeLog(new DFastLog_fastLog(), false);
	}

	private void canComputeLog(BaseLog f, boolean edgeCases)
	{
		testLog(f, Float.NaN, edgeCases);
		testLog(f, Float.NEGATIVE_INFINITY, edgeCases);
		testLog(f, -Float.MAX_VALUE, edgeCases);
		testLog(f, -Float.MIN_VALUE, edgeCases);
		testLog(f, -2, edgeCases);
		testLog(f, -1, edgeCases);
		testLog(f, -0, edgeCases);
		testLog(f, 0, true);
		testLog(f, Float.MIN_VALUE, false); // Not enough precision to test this
		testLog(f, 1e-10f, true);
		testLog(f, 1, true);
		testLog(f, 2, true);
		testLog(f, 2048, true);
		testLog(f, Float.MAX_VALUE, true);
		testLog(f, Float.POSITIVE_INFINITY, edgeCases);
	}

	private void testLog(BaseLog f, float v, boolean test)
	{
		double e = Math.log(v);
		double o = f.log(v);
		double error = DoubleEquality.relativeError(e, o);
		System.out.printf("%s v=%g : %f vs %s (%g)\n", f.name, v, e, o, error);
		if (test)
		{
			if (Double.isNaN(e) && Double.isNaN(o))
				return;
			if (e == o)
				return;
			Assert.assertTrue(error < 1e-4);
		}
	}

	@Test
	public void canComputeDoubleFastLog_log()
	{
		canComputeDoubleLog(new FFastLog_log(), true);
	}

	@Test
	public void canComputeDoubleFastLog_fastLog()
	{
		canComputeDoubleLog(new FFastLog_fastLog(), false);
	}

	@Test
	public void canComputeDoubleDFast_log()
	{
		canComputeDoubleLog(new DFastLog_log(), true);
	}

	@Test
	public void canComputeDoubleDFast_fastLog()
	{
		canComputeDoubleLog(new DFastLog_fastLog(), false);
	}

	private void canComputeDoubleLog(BaseLog f, boolean edgeCases)
	{
		testDoubleLog(f, Double.NaN, edgeCases);
		testDoubleLog(f, Double.NEGATIVE_INFINITY, edgeCases);
		testDoubleLog(f, -Double.MAX_VALUE, edgeCases);
		testDoubleLog(f, -Double.MIN_VALUE, edgeCases);
		testDoubleLog(f, -2, edgeCases);
		testDoubleLog(f, -1, edgeCases);
		testDoubleLog(f, -0, edgeCases);
		testDoubleLog(f, 0, true);
		testDoubleLog(f, Double.MIN_VALUE, false); // Not enough precision to test this
		testDoubleLog(f, 1e-10f, true);
		testDoubleLog(f, 1, true);
		testDoubleLog(f, 2, true);
		testDoubleLog(f, 2048, true);
		testDoubleLog(f, Double.MAX_VALUE / 2, f.isDoublePrecision());
		testDoubleLog(f, Double.MAX_VALUE, f.isDoublePrecision());
		testDoubleLog(f, Double.POSITIVE_INFINITY, edgeCases);
	}

	private void testDoubleLog(BaseLog f, double v, boolean test)
	{
		double e = Math.log(v);
		double o = f.log(v);
		double error = DoubleEquality.relativeError(e, o);
		System.out.printf("%s v=%g : %f vs %s (%g)\n", f.name, v, e, o, error);
		if (test)
		{
			if (Double.isNaN(e) && Double.isNaN(o))
				return;
			if (e == o)
				return;
			Assert.assertTrue(error < 1e-4);
		}
	}

	// A robust error test using a uniform random number, report min,av,max,sd of the error.
	// The error 'should' always be negative as the truncation rounds down. However the table
	// pre-computes using an exponent offset which can lead to rounding up.

	@Test
	public void canTestFloatError()
	{
		// All float values is a lot so we do a representative set
		RandomGenerator r = new Well19937c(30051977);
		double lower = Float.MIN_VALUE, upper = Float.MAX_VALUE;
		float[] d = new float[100000];
		float[] logD = new float[d.length];
		for (int i = 0; i < d.length; i++)
		{
			float v = (float) nextUniform(r, lower, upper);
			d[i] = v;
			logD[i] = (float) Math.log(v);
		}

		for (int n = 0; n <= 23; n++)
		{
			canTestFloatError(n, d, logD);
		}
	}

	private double nextUniform(RandomGenerator r, double lower, double upper)
	{
		double u = r.nextDouble();
		return u * upper + (1.0 - u) * lower;
	}

	private class FPair
	{
		int i = 0;
		float f = 0;
	}

	private class Stats
	{
		double min, max, s, ss;
		int n = 1;

		Stats(double e)
		{
			min = max = s = e;
			ss = e * e;
		}

		void add(double e)
		{
			if (min > e)
				min = e;
			else if (max < e)
				max = e;
			s += e;
			ss += e * e;
			n++;
		}

		double getMean()
		{
			return s / n;
		};

		double getSD()
		{
			double sd = ss - ((double) s * s) / n;
			if (sd > 0.0)
				return Math.sqrt(sd / (n - 1));
			else
				return 0.0;
		}

		String summary()
		{
			return String.format("min=%g, max=%g, mean=%g, sd=%g", min, max, getMean(), getSD());
		}
	}

	private void canTestFloatError(int precision, float[] d, float[] logD)
	{
		FFastLog f = new FFastLog(precision);
		FPair pair = new FPair();
		if (!next(f, pair, d))
			return;

		double delta = logD[pair.i - 1] - pair.f;
		delta = Math.abs(delta);
		Stats s1 = new Stats(delta);
		Stats s2 = new Stats(Math.abs(delta / logD[pair.i - 1]));
		while (next(f, pair, d))
		{
			delta = logD[pair.i - 1] - pair.f;
			delta = Math.abs(delta);
			s1.add(delta);
			s2.add(Math.abs(delta / logD[pair.i - 1]));
		}
		System.out.printf("FloatError n=%d, c=%d : %s : relative %s\n", precision, s1.n, s1.summary(), s2.summary());
	}

	private boolean next(FFastLog f, FPair pair, float[] d)
	{
		while (pair.i < d.length)
		{
			pair.f = f.log(d[pair.i++]);
			if (pair.f != Float.NEGATIVE_INFINITY)
				return true;
		}
		return false;
	}

	@Test
	public void canTestDoubleError()
	{
		// All float values is a lot so we do a representative set
		RandomGenerator r = new Well19937c(30051977);
		double lower = Double.MIN_VALUE, upper = Double.MAX_VALUE;
		double[] d = new double[100000];
		double[] logD = new double[d.length];
		for (int i = 0; i < d.length; i++)
		{
			double v = nextUniform(r, lower, upper);
			d[i] = v;
			logD[i] = Math.log(v);
		}

		for (int n = 0; n <= 23; n++)
		{
			canTestDoubleError(n, d, logD);
		}
	}

	private class DPair
	{
		int i = 0;
		double f = 0;
	}

	private void canTestDoubleError(int precision, double[] d, double[] logD)
	{
		DFastLog f = new DFastLog(precision);
		DPair pair = new DPair();
		if (!next(f, pair, d))
			return;

		double delta = logD[pair.i - 1] - pair.f;
		delta = Math.abs(delta);
		Stats s1 = new Stats(delta);
		Stats s2 = new Stats(Math.abs(delta / logD[pair.i - 1]));
		while (next(f, pair, d))
		{
			delta = logD[pair.i - 1] - pair.f;
			delta = Math.abs(delta);
			s1.add(delta);
			s2.add(Math.abs(delta / logD[pair.i - 1]));
		}
		System.out.printf("DoubleError n=%d, c=%d : %s : relative %s\n", precision, s1.n, s1.summary(), s2.summary());
	}

	private boolean next(DFastLog f, DPair pair, double[] d)
	{
		while (pair.i < d.length)
		{
			pair.f = f.log(d[pair.i++]);
			if (pair.f != Double.NEGATIVE_INFINITY)
				return true;
		}
		return false;
	}

	// Speed test of float/double version verses the Math.log.

	private class FloatTimingTask extends BaseTimingTask
	{
		BaseLog log;
		float[] x;

		public FloatTimingTask(BaseLog log, int q, float[] x)
		{
			super(log.name + " q=" + q);
			this.log = log;
			this.x = x;
		}

		public int getSize()
		{
			return 1;
		}

		public Object getData(int i)
		{
			return null;
		}

		public Object run(Object data)
		{
			for (int i = 0; i < x.length; i++)
				log.log(x[i]);
			return null;
		}
	}

	@Test
	public void canTestFloatSpeed()
	{
		RandomGenerator r = new Well19937c(30051977);
		double upper = Float.MAX_VALUE;
		float[] x = new float[1000000];
		for (int i = 0; i < x.length; i++)
		{
			x[i] = (float) (r.nextDouble() * upper);
		}

		TimingService ts = new TimingService(5);
		ts.execute(new FloatTimingTask(new Math_log(), 0, x));
		for (int q : new int[] { 0, 7, 8, 9, 10, 11, 12, 13 })
		{
			int n = 23 - q;
			FFastLog f = new FFastLog(n);
			//ts.execute(new FloatTimingTask(new FFastLog_log2(f), q, x));
			//ts.execute(new FloatTimingTask(new FFastLog_fastLog2(f), q, x));
			ts.execute(new FloatTimingTask(new FFastLog_log(f), q, x));
			ts.execute(new FloatTimingTask(new FFastLog_fastLog(f), q, x));
		}

		int size = ts.getSize();
		ts.repeat(size);
		ts.report(size);
	}

	private class DoubleTimingTask extends BaseTimingTask
	{
		BaseLog log;
		double[] x;

		public DoubleTimingTask(BaseLog log, int q, double[] x)
		{
			super(log.name + " q=" + q);
			this.log = log;
			this.x = x;
		}

		public int getSize()
		{
			return 1;
		}

		public Object getData(int i)
		{
			return null;
		}

		public Object run(Object data)
		{
			for (int i = 0; i < x.length; i++)
				log.log(x[i]);
			return null;
		}
	}

	@Test
	public void canTestDoubleSpeed()
	{
		RandomGenerator r = new Well19937c(30051977);
		double upper = Double.MAX_VALUE;
		double[] x = new double[1000000];
		for (int i = 0; i < x.length; i++)
		{
			x[i] = r.nextDouble() * upper;
		}

		TimingService ts = new TimingService(5);
		ts.execute(new DoubleTimingTask(new Math_log(), 0, x));
		for (int q : new int[] { 0, 7, 8, 9, 10, 11, 12, 13 })
		{
			int n = 23 - q;
			DFastLog f = new DFastLog(n);
			//ts.execute(new DoubleTimingTask(new DFastLog_log2(f), q, x));
			//ts.execute(new DoubleTimingTask(new DFastLog_fastLog2(f), q, x));
			ts.execute(new DoubleTimingTask(new DFastLog_log(f), q, x));
			ts.execute(new DoubleTimingTask(new DFastLog_fastLog(f), q, x));
		}

		int size = ts.getSize();
		ts.repeat(size);
		ts.report(size);
	}
	
	@Test
	public void canTestFloatVsDoubleSpeed()
	{
		RandomGenerator r = new Well19937c(30051977);
		double upper = Double.MAX_VALUE;
		double[] x = new double[1000000];
		for (int i = 0; i < x.length; i++)
		{
			x[i] = r.nextDouble() * upper;
		}

		TimingService ts = new TimingService(5);
		ts.execute(new DoubleTimingTask(new Math_log(), 0, x));
		for (int q : new int[] { 0, 7, 8, 9, 10, 11, 12, 13 })
		{
			int n = 23 - q;
			FFastLog f = new FFastLog(n);
			DFastLog f2 = new DFastLog(n);
			ts.execute(new DoubleTimingTask(new FFastLog_log(f), q, x));
			ts.execute(new DoubleTimingTask(new FFastLog_fastLog(f), q, x));
			ts.execute(new DoubleTimingTask(new DFastLog_log(f2), q, x));
			ts.execute(new DoubleTimingTask(new DFastLog_fastLog(f2), q, x));
		}

		int size = ts.getSize();
		ts.repeat(size);
		ts.report(size);
	}
}
