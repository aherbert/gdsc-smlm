package gdsc.smlm.function;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.apache.commons.math3.util.FastMath;
import org.junit.Assert;
import org.junit.Assume;
import org.junit.Test;

import gdsc.core.test.BaseTimingTask;
import gdsc.core.test.TimingService;
import gdsc.core.utils.BitFlags;
import gdsc.core.utils.DoubleEquality;
import gdsc.core.utils.FloatEquality;
import gdsc.core.utils.Maths;
import gdsc.core.utils.TurboList;
import gdsc.smlm.function.ICSIFastLog.DataType;

@SuppressWarnings("unused")
public class FastLogTest
{
	ICSIFastLog iLog = ICSIFastLog.create(DataType.BOTH);
	FFastLog fLog = new FFastLog();
	DFastLog dLog = new DFastLog();
	TurboLog tLog = new TurboLog();

	//@formatter:off
	private class MathLog extends FastLog
	{
		@Override
		public double getBase()	{ return Math.E;}
		@Override
		public double getScale() { return LN2; }
		@Override
		public int getN() { return 52; }
		
		@Override
		public float log(float x) {	return (float) Math.log(x);	}
		@Override
		public float log2(float x) { return (float) (log(x) / LN2);	}
		@Override
		public float log(double x) { return (float) Math.log(x); }
		@Override
		public float log2(double x) { return (float) (log(x) / LN2); }
		@Override
		public double logD(double x) { return Math.log(x); }
		@Override
		public double  log2D(double x) { return (log(x) / LN2); }
		
		@Override
		public float fastLog(float x) { return log(x); }
		@Override
		public float fastLog2(float x) { return log2(x); }
		@Override
		public float fastLog(double x) { return log(x);	}
		@Override
		public float fastLog2(double x) { return log2(x); }
		@Override
		public double fastLogD(double x) { return log(x);	}
		@Override
		public double fastLog2D(double x) { return log2(x); }
	}
	private abstract class BaseTestLog
	{
		FastLog fl;
		String name;
		BaseTestLog(FastLog fl) { 
			this.name = this.getClass().getSimpleName() + " " + fl.getClass().getSimpleName();
			this.fl=fl; }
		abstract float log(float x);
		abstract double log(double x);
		int getN() { return fl.getN(); };
	}
	private class TestLog extends BaseTestLog
	{
		TestLog(FastLog fl) { super(fl); }
		float log(float x) { return fl.log(x); }
		double log(double x) { return fl.logD(x); }
	}
	private class TestFastLog extends BaseTestLog
	{
		TestFastLog(FastLog fl) { super(fl); }
		float log(float x) { return fl.fastLog(x); }
		double log(double x) { return fl.fastLogD(x); }
	}
	// To test Math.log(1+x).
	// This is what is used in the MLE LVM gradient calculator
	private class Test1PLog extends BaseTestLog
	{
		Test1PLog(FastLog fl) { super(fl); }
		float log(float x) { return (float) Math.log(1+x); }
		double log(double x) { return Math.log(1+x); }
	}
	private class TestLog1P extends BaseTestLog
	{
		TestLog1P(FastLog fl) { super(fl); }
		float log(float x) { return (float) Math.log1p(x); }
		double log(double x) { return Math.log1p(x); }
	}
	private class TestLog1PApache extends BaseTestLog
	{
		TestLog1PApache(FastLog fl) { super(fl); }
		float log(float x) { return (float) FastMath.log1p(x); }
		double log(double x) { return FastMath.log1p(x); }
	}
	//@formatter:on

	@Test
	public void canComputeFFastLog_fastLog()
	{
		canComputeLog(new TestFastLog(fLog), false);
	}

	@Test
	public void canComputeFFastLog_log()
	{
		canComputeLog(new TestLog(fLog), true);
	}

	@Test
	public void canComputeDFastLog_fastLog()
	{
		canComputeLog(new TestFastLog(dLog), false);
	}

	@Test
	public void canComputeDFastLog_log()
	{
		canComputeLog(new TestLog(dLog), true);
	}

	@Test
	public void canComputeICSCFastLog_fastLog()
	{
		canComputeLog(new TestFastLog(iLog), false);
	}

	@Test
	public void canComputeICSCFastLog_log()
	{
		canComputeLog(new TestLog(iLog), true);
	}

	@Test
	public void canComputeTurboLog_fastLog()
	{
		canComputeLog(new TestFastLog(tLog), false);
	}

	@Test
	public void canComputeTurboLog_log()
	{
		canComputeLog(new TestLog(tLog), true);
	}

	private void canComputeLog(BaseTestLog f, boolean edgeCases)
	{
		testLog(f, Float.NaN, edgeCases);
		testLog(f, Float.NEGATIVE_INFINITY, edgeCases);
		testLog(f, -Float.MAX_VALUE, edgeCases);
		testLog(f, -Float.MIN_VALUE, edgeCases);
		testLog(f, -2f, edgeCases);
		testLog(f, -1f, edgeCases);
		testLog(f, -0f, edgeCases);
		testLog(f, 0f, edgeCases);
		testLog(f, Float.MIN_VALUE, false); // Not enough precision to test this
		testLog(f, 1e-10f, true);
		testLog(f, 1f, true);
		testLog(f, 2f, true);
		testLog(f, 2048f, true);
		testLog(f, Float.MAX_VALUE, true);
		testLog(f, Float.POSITIVE_INFINITY, edgeCases);
	}

	private void testLog(BaseTestLog f, float v, boolean test)
	{
		float e = (float) Math.log(v);
		float o = f.log(v);
		float error = FloatEquality.relativeError(e, o);
		System.out.printf("%s v=%g : %f vs %s (%g)\n", f.name, v, e, o, error);
		if (test)
		{
			if (Double.isNaN(e) && Double.isNaN(o))
				return;
			if (e == o)
				return;
			Assert.assertTrue(error < 1e-4f);
		}
	}

	@Test
	public void canComputeDoubleFFast_fastLog()
	{
		canComputeDoubleLog(new TestFastLog(fLog), false);
	}

	@Test
	public void canComputeDoubleFFastLog_log()
	{
		canComputeDoubleLog(new TestLog(fLog), true);
	}

	@Test
	public void canComputeDoubleDFast_fastLog()
	{
		canComputeDoubleLog(new TestFastLog(dLog), false);
	}

	@Test
	public void canComputeDoubleDFastLog_log()
	{
		canComputeDoubleLog(new TestLog(dLog), true);
	}

	@Test
	public void canComputeDoubleICSCFast_fastLog()
	{
		canComputeDoubleLog(new TestFastLog(iLog), false);
	}

	@Test
	public void canComputeDoubleICSCFastLog_log()
	{
		canComputeDoubleLog(new TestLog(iLog), true);
	}

	@Test
	public void canComputeDoubleTurbo_fastLog()
	{
		canComputeDoubleLog(new TestFastLog(tLog), false);
	}

	@Test
	public void canComputeDoubleTurboLog_log()
	{
		canComputeDoubleLog(new TestLog(tLog), true);
	}

	private void canComputeDoubleLog(BaseTestLog f, boolean edgeCases)
	{
		testDoubleLog(f, Double.NaN, edgeCases);
		testDoubleLog(f, Double.NEGATIVE_INFINITY, edgeCases);
		testDoubleLog(f, -Double.MAX_VALUE, edgeCases);
		testDoubleLog(f, -Double.MIN_VALUE, edgeCases);
		testDoubleLog(f, -2d, edgeCases);
		testDoubleLog(f, -1d, edgeCases);
		testDoubleLog(f, -0d, edgeCases);
		testDoubleLog(f, 0d, edgeCases);
		testDoubleLog(f, Double.MIN_VALUE, false); // Not enough precision to test this
		testDoubleLog(f, 1e-10, true);
		testDoubleLog(f, 1d, true);
		testDoubleLog(f, 2d, true);
		testDoubleLog(f, 2048d, true);
		testDoubleLog(f, Double.MAX_VALUE / 2, true);
		testDoubleLog(f, Double.MAX_VALUE, true);
		testDoubleLog(f, Double.POSITIVE_INFINITY, edgeCases);
	}

	private void testDoubleLog(BaseTestLog f, double v, boolean test)
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
		float[] d = generateRandomFloats(1000000);
		float[] logD = new float[d.length];
		for (int i = 0; i < d.length; i++)
		{
			logD[i] = (float) Math.log(d[i]);
		}

		//int min=0,max=23;
		int min = 13, max = 13;

		for (int n = min; n <= max; n++)
		{
			canTestFloatError(new TestFastLog(ICSIFastLog.create(n, DataType.FLOAT)), d, logD);
		}
		for (int n = min; n <= max; n++)
		{
			canTestFloatError(new TestFastLog(new FFastLog(n)), d, logD);
		}
		for (int n = min; n <= max; n++)
		{
			canTestFloatError(new TestFastLog(new DFastLog(n)), d, logD);
		}
		for (int n = min; n <= max; n++)
		{
			canTestFloatError(new TestFastLog(new TurboLog(n)), d, logD);
		}
	}

	private float[] generateRandomFloats(int n)
	{
		RandomGenerator r = new Well19937c(30051977);
		float[] d = new float[n];
		for (int i = 0; i < d.length; i++)
			d[i] = nextUniformFloat(r);
		return d;
	}

	private float nextUniformFloat(RandomGenerator r)
	{
		int u = r.nextInt();
		// Mask out sign and the last bit of the exponent (avoid infinity and NaN)
		u = BitFlags.unset(u, 0x80000000 | 0x00800000);
		//assert ((u >> 23) & 0xff) < 255;
		return Float.intBitsToFloat(u);
	}

	@Test
	public void canTestFloatErrorRange()
	{
		Assume.assumeTrue(true);

		TurboList<TestFastLog> test = new TurboList<TestFastLog>();
		int n = 13;
		test.add(new TestFastLog(ICSIFastLog.create(n, DataType.FLOAT)));
		test.add(new TestFastLog(new FFastLog(n)));
		test.add(new TestFastLog(new DFastLog(n)));
		test.add(new TestFastLog(new TurboLog(n)));

		// Full range in blocks.
		// Only when the number is around 1 or min value are there significant errors
		float[] d = null, logD = null;
		
		// All
		//testFloatErrorRange(test, n, d, logD, 0, 255, 0);
		
		// Only a problem around min value and x==1
		//testFloatErrorRange(test, n, d, logD, 0, 2, 0);
		testFloatErrorRange(test, n, d, logD, 125, 130, 0);
		//testFloatErrorRange(test, n, d, logD, 253, 255, 0);
	}

	private void testFloatErrorRange(TurboList<TestFastLog> test, int n, float[] d, float[] logD, int mine, int maxe, int ee)
	{
		for (int e = mine; e < maxe; e += ee + 1)
		{
			d = generateFloats(e, e + ee, d);
			if (logD == null || logD.length < n)
				logD = new float[d.length];
			for (int i = 0; i < d.length; i++)
			{
				logD[i] = (float) Math.log(d[i]);
			}
			System.out.printf("e=%d-%d\n", e, e + ee);
			for (TestFastLog f : test)
				canTestFloatError(f, d, logD);
		}
	}

	private float[] generateFloats(int mine, int maxe, float[] d)
	{
		// Mantissa = 23-bit, Exponent = 8-bit
		int mbits = 23;
		mine = Maths.clip(0, 255, mine);
		maxe = Maths.clip(0, 255, maxe);
		if (mine > maxe)
			throw new IllegalStateException();
		int mn = (1 << mbits);
		int n = mn * (maxe - mine + 1);
		if (d == null || d.length < n)
			d = new float[n];
		int i = 0;
		for (int m = 0; m < mn; m++)
		{
			for (int e = mine; e <= maxe; e++)
			{
				int bits = m | (e << 23);
				float v = Float.intBitsToFloat(bits);
				//System.out.printf("%g = %s\n", v, Integer.toBinaryString(bits));
				d[i++] = v;
			}
		}
		return d;
	}

	private class FPair
	{
		int i = 0;
		float f = 0;
	}

	private class Stats
	{
		double min, max, s, ss, minv, maxv;
		int n = 1;

		Stats(double e, double v)
		{
			min = max = s = e;
			minv = maxv = v;
			ss = e * e;
		}

		void add(double e, double v)
		{
			if (min > e)
			{
				min = e;
				minv = v;
			}
			else if (max < e)
			{
				max = e;
				maxv = v;
			}
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
			return String.format("min=%s (%s), max=%s (%s), mean=%s, sd=%g", min, minv, max, maxv, getMean(), getSD());
		}
	}

	private void canTestFloatError(BaseTestLog f, float[] d, float[] logD)
	{
		FPair pair = new FPair();
		if (!next(f, pair, d))
			return;

		float v = logD[pair.i - 1];
		double delta = v - pair.f;
		delta = Math.abs(delta);
		//		if (delta > 1)
		//		{
		//			System.out.printf("Big error: %f %f\n", v, d[pair.i-1]);
		//		}
		Stats s1 = new Stats(delta, d[pair.i - 1]);
		Stats s2 = (v != 0) ? new Stats(Math.abs(delta / v), d[pair.i - 1]) : new Stats(0, d[pair.i - 1]);
		while (next(f, pair, d))
		{
			v = logD[pair.i - 1];
			delta = v - pair.f;
			delta = Math.abs(delta);
			//if (delta > 5)
			//{
			//	System.out.printf("Big error: [%g] %f %f %f\n", d[pair.i - 1], v, pair.f,
			//			v));
			//}
			s1.add(delta, d[pair.i - 1]);
			if (v != 0)
				s2.add(Math.abs(delta / v), d[pair.i - 1]);
		}
		System.out.printf("%s, n=%d, c=%d : %s : relative %s\n", f.name, f.getN(), s1.n, s1.summary(), s2.summary());
	}

	private boolean next(BaseTestLog f, FPair pair, float[] d)
	{
		while (pair.i < d.length)
		{
			float x = d[pair.i++];
			if (x == 0)// Skip infinity
				continue;
			pair.f = f.log(x);
			if (pair.f != Float.NEGATIVE_INFINITY)
				return true;
			//System.out.printf("%g\n", d[pair.i - 1]);
		}
		return false;
	}

	@Test
	public void canTestDoubleError()
	{
		// All float values is a lot so we do a representative set
		RandomGenerator r = new Well19937c(30051977);
		double lower = Double.MIN_VALUE, upper = Double.MAX_VALUE;
		double[] d = new double[10000000];
		double[] logD = new double[d.length];
		for (int i = 0; i < d.length; i++)
		{
			double v = nextUniformDouble(r);
			d[i] = v;
			logD[i] = Math.log(v);
		}

		//int min = 0, max = 23;
		int min = 13, max = 13;

		//		for (int n = min; n <= max; n++)
		//		{
		//			canTestDoubleError(new TestFastLog(ICSIFastLog.create(n, DataType.DOUBLE)), d, logD);
		//		}
		//		for (int n = min; n <= max; n++)
		//		{
		//			canTestDoubleError(new TestFastLog(new DFastLog(n)), d, logD);
		//		}
		//		for (int n = min; n <= max; n++)
		//		{
		//			canTestDoubleError(new TestFastLog(new FFastLog(n)), d, logD);
		//		}
		for (int n = min; n <= max; n++)
		{
			canTestDoubleError(new TestFastLog(new TurboLog(n)), d, logD);
		}
	}

	@Test
	public void canTestDoubleErrorLog1P()
	{
		// All float values is a lot so we do a representative set
		RandomGenerator r = new Well19937c(30051977);
		double lower = Double.MIN_VALUE, upper = Double.MAX_VALUE;
		double[] d = new double[100000];
		double[] logD = new double[d.length];
		for (int i = 0; i < d.length; i++)
		{
			double v = nextUniformDouble(r);
			d[i] = v;
			logD[i] = Math.log1p(v);
		}

		canTestDoubleError(new Test1PLog(new MathLog()), d, logD);
		canTestDoubleError(new TestLog1P(new MathLog()), d, logD);
	}

	private double nextUniformDouble(RandomGenerator r)
	{
		long u = r.nextLong();
		// Mask out sign and the last bit of the exponent (avoid infinity and NaN)
		u &= ~(0x8000000000000000L | 0x0010000000000000L);
		//assert ((u >> 52) & 0x7ffL) < 2047;
		return Double.longBitsToDouble(u);
	}

	private class DPair
	{
		int i = 0;
		double f = 0;
	}

	private void canTestDoubleError(BaseTestLog f, double[] d, double[] logD)
	{
		DPair pair = new DPair();
		if (!next(f, pair, d))
			return;

		double v = logD[pair.i - 1];
		double delta = v - pair.f;
		delta = Math.abs(delta);
		Stats s1 = new Stats(delta, d[pair.i - 1]);
		Stats s2 = (v != 0) ? new Stats(Math.abs(delta / v), d[pair.i - 1]) : new Stats(0, d[pair.i - 1]);
		while (next(f, pair, d))
		{
			v = logD[pair.i - 1];
			//System.out.printf("%g vs %g\n", v, pair.f);
			delta = v - pair.f;
			delta = Math.abs(delta);
			s1.add(delta, d[pair.i - 1]);
			if (v != 0)
				s2.add(Math.abs(delta / v), d[pair.i - 1]);
		}
		System.out.printf("%s, n=%d, c=%d : %s : relative %s\n", f.name, f.getN(), s1.n, s1.summary(), s2.summary());
	}

	private boolean next(BaseTestLog f, DPair pair, double[] d)
	{
		while (pair.i < d.length)
		{
			double x = d[pair.i++];
			if (x == 0)// Skip infinity
				continue;
			pair.f = f.log(d[pair.i++]);
			if (pair.f != Double.NEGATIVE_INFINITY)
				return true;
			//System.out.printf("%g\n", d[pair.i - 1]);
		}
		return false;
	}

	// Speed test of float/double version verses the Math.log.
	private abstract class DummyTimingTask extends BaseTimingTask
	{
		BaseTestLog log;
		int q;

		public DummyTimingTask(String name, BaseTestLog log, int q)
		{
			super(name + " " + log.name + " q=" + q);
			this.log = log;
			this.q = q;
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

	private class FloatTimingTask extends DummyTimingTask
	{
		float[] x;

		public FloatTimingTask(BaseTestLog log, int q, float[] x)
		{
			super("log(float)", log, q);
			this.x = x;
		}

		public Object run(Object data)
		{
			float[] r = new float[x.length];
			for (int i = 0; i < x.length; i++)
				r[i] = log.log(x[i]);
			return r;
		}
	}

	@Test
	public void canTestFloatSpeed()
	{
		RandomGenerator r = new Well19937c(30051977);
		float[] x = new float[1000000];
		for (int i = 0; i < x.length; i++)
		{
			x[i] = nextUniformFloat(r);
		}

		TimingService ts = new TimingService(5);
		ts.execute(new FloatTimingTask(new TestLog(new MathLog()), 0, x));
		for (int q : new int[] { 11 })
		//for (int q : new int[] { 0, 7, 8, 9, 10, 11, 12, 13 })
		{
			int n = 23 - q;
			ICSIFastLog f = ICSIFastLog.create(n, DataType.FLOAT);
			ts.execute(new FloatTimingTask(new TestLog(f), q, x));
			ts.execute(new FloatTimingTask(new TestFastLog(f), q, x));
			FFastLog ff = new FFastLog(n);
			ts.execute(new FloatTimingTask(new TestLog(ff), q, x));
			ts.execute(new FloatTimingTask(new TestFastLog(ff), q, x));
			DFastLog df = new DFastLog(n);
			ts.execute(new FloatTimingTask(new TestLog(df), q, x));
			ts.execute(new FloatTimingTask(new TestFastLog(df), q, x));
			TurboLog tf = new TurboLog(n);
			ts.execute(new FloatTimingTask(new TestLog(tf), q, x));
			ts.execute(new FloatTimingTask(new TestFastLog(tf), q, x));
		}

		int size = ts.getSize();
		ts.repeat(size);
		ts.report(size);
	}

	private class DoubleTimingTask extends DummyTimingTask
	{
		double[] x;

		public DoubleTimingTask(BaseTestLog log, int q, double[] x)
		{
			super("log(double)", log, q);
			this.x = x;
		}

		public Object run(Object data)
		{
			double[] r = new double[x.length];
			for (int i = 0; i < x.length; i++)
				r[i] = log.log(x[i]);
			return r;
		}
	}

	private class DoubleToFloatTimingTask extends DummyTimingTask
	{
		double[] x;
		float[] xf;

		public DoubleToFloatTimingTask(BaseTestLog log, int q, double[] x, float[] xf)
		{
			super("log((float)double)", log, q);
			this.x = x;
			this.xf = xf;
		}

		public Object run(Object data)
		{
			double[] r = new double[x.length];
			for (int i = 0; i < x.length; i++)
				r[i] = log.log((float) x[i]);
			//r[i] = log.log(xf[i]);
			return r;
		}
	}

	@Test
	public void canTestDoubleSpeed()
	{
		RandomGenerator r = new Well19937c(30051977);
		double[] x = new double[1000000];
		for (int i = 0; i < x.length; i++)
		{
			x[i] = nextUniformDouble(r);
		}

		TimingService ts = new TimingService(5);
		ts.execute(new DoubleTimingTask(new TestLog(new MathLog()), 0, x));
		for (int q : new int[] { 11 })
		//for (int q : new int[] { 0, 7, 8, 9, 10, 11, 12, 13 })
		{
			int n = 23 - q;
			ICSIFastLog f = ICSIFastLog.create(n, DataType.DOUBLE);
			ts.execute(new DoubleTimingTask(new TestLog(f), q, x));
			ts.execute(new DoubleTimingTask(new TestFastLog(f), q, x));
			DFastLog df = new DFastLog(n);
			//ts.execute(new DoubleTimingTask(new DFastLog_log2(f), q, x));
			//ts.execute(new DoubleTimingTask(new DTestFastLog2(f), q, x));
			ts.execute(new DoubleTimingTask(new TestLog(df), q, x));
			ts.execute(new DoubleTimingTask(new TestFastLog(df), q, x));
			TurboLog tf = new TurboLog(n);
			ts.execute(new DoubleTimingTask(new TestLog(tf), q, x));
			ts.execute(new DoubleTimingTask(new TestFastLog(tf), q, x));
		}

		int size = ts.getSize();
		ts.repeat(size);
		ts.report(size);
	}

	@Test
	public void canTestDoubleSpeedLog1P()
	{
		RandomGenerator r = new Well19937c(30051977);
		double[] x = new double[1000000];
		for (int i = 0; i < x.length; i++)
		{
			x[i] = nextUniformDouble(r);
		}

		MathLog f = new MathLog();

		TimingService ts = new TimingService(5);
		//ts.execute(new DoubleTimingTask(new TestLog(f), 0, x));
		ts.execute(new DoubleTimingTask(new Test1PLog(f), 0, x));
		ts.execute(new DoubleTimingTask(new TestLog1P(f), 0, x));
		ts.execute(new DoubleTimingTask(new TestLog1PApache(f), 0, x));
		//ts.execute(new DoubleTimingTask(new TestLog(f), 0, x));
		ts.execute(new DoubleTimingTask(new Test1PLog(f), 0, x));
		ts.execute(new DoubleTimingTask(new TestLog1P(f), 0, x));
		ts.execute(new DoubleTimingTask(new TestLog1PApache(f), 0, x));

		int size = ts.getSize();
		ts.repeat(size);
		ts.report(size);
	}

	@Test
	public void canTestFloatVsDoubleSpeed()
	{
		RandomGenerator r = new Well19937c(30051977);
		double[] x = new double[1000000];
		float[] xf = new float[x.length];
		for (int i = 0; i < x.length; i++)
		{
			x[i] = nextUniformFloat(r);
			xf[i] = (float) x[i];
		}

		TimingService ts = new TimingService(5);
		ts.execute(new DoubleTimingTask(new TestLog(new MathLog()), 0, x));
		for (int q : new int[] { 11 })
		//for (int q : new int[] { 0, 7, 8, 9, 10, 11, 12, 13 })
		{
			int n = 23 - q;
			ICSIFastLog ff = ICSIFastLog.create(n, DataType.FLOAT);
			ICSIFastLog fd = ICSIFastLog.create(n, DataType.DOUBLE);
			ICSIFastLog ff2 = ICSIFastLog.create(n, DataType.FLOAT);
			ts.execute(new DoubleToFloatTimingTask(new TestLog(ff), q, x, xf));
			ts.execute(new DoubleToFloatTimingTask(new TestFastLog(ff), q, x, xf));
			ts.execute(new FloatTimingTask(new TestLog(ff2), q, xf));
			ts.execute(new FloatTimingTask(new TestFastLog(ff2), q, xf));
			ts.execute(new DoubleTimingTask(new TestLog(fd), q, x));
			ts.execute(new DoubleTimingTask(new TestFastLog(fd), q, x));

			//			ts.execute(new DoubleToFloatTimingTask(new TestLog(ff), q, x, xf));
			//			ts.execute(new DoubleToFloatTimingTask(new TestFastLog(ff), q, x, xf));
			//			ts.execute(new FloatTimingTask(new TestLog(ff2), q, xf));
			//			ts.execute(new FloatTimingTask(new TestFastLog(ff2), q, xf));
			//			ts.execute(new DoubleTimingTask(new TestLog(fd), q, x));
			//			ts.execute(new DoubleTimingTask(new TestFastLog(fd), q, x));

			TurboLog tf = new TurboLog(n);
			ts.execute(new DoubleToFloatTimingTask(new TestLog(tf), q, x, xf));
			ts.execute(new DoubleToFloatTimingTask(new TestFastLog(tf), q, x, xf));
			ts.execute(new FloatTimingTask(new TestLog(tf), q, xf));
			ts.execute(new FloatTimingTask(new TestFastLog(tf), q, xf));
			ts.execute(new DoubleTimingTask(new TestLog(tf), q, x));
			ts.execute(new DoubleTimingTask(new TestFastLog(tf), q, x));

			// Slower as the look-up table is bigger
			FFastLog f1 = new FFastLog(n);
			DFastLog f2 = new DFastLog(n);
			ts.execute(new FloatTimingTask(new TestLog(f1), q, xf));
			ts.execute(new FloatTimingTask(new TestFastLog(f1), q, xf));
			ts.execute(new DoubleTimingTask(new TestLog(f2), q, x));
			ts.execute(new DoubleTimingTask(new TestFastLog(f2), q, x));
		}

		int size = ts.getSize();
		ts.repeat(size);
		//ts.report();
		ts.report(size);
	}
}
