package gdsc.smlm.function;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.Precision;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.test.BaseTimingTask;
import gdsc.core.test.TimingService;
import gdsc.core.utils.DoubleEquality;

public class FastLogTest
{
	FFastLog fLog = new FFastLog(11);
	DFastLog dLog = new DFastLog(11);

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
	private class FastLog_log extends BaseLog
	{
		float log(float x) { return fLog.log(x); }
		double log(double x) { return fLog.log(x); }
		boolean isDoublePrecision() { return false; }
	}
	private class FastLog_fastLog extends BaseLog
	{
		float log(float x) { return fLog.fastLog(x); }
		double log(double x) { return fLog.fastLog(x); }
		boolean isDoublePrecision() { return false; }
	}
	private class DFastLog_log extends BaseLog
	{
		float log(float x) { return (float) dLog.log(x); }
		double log(double x) { return dLog.log(x); }
		boolean isDoublePrecision() { return true; }
	}
	private class DFastLog_fastLog extends BaseLog
	{
		float log(float x) { return (float) dLog.fastLog(x); }
		double log(double x) { return dLog.fastLog(x); }
		boolean isDoublePrecision() { return true; }
	}
	//@formatter:on

	@Test
	public void canComputeFastLog_log()
	{
		canComputeLog(new FastLog_log(), true);
	}

	@Test
	public void canComputeFastLog_fastLog()
	{
		canComputeLog(new FastLog_fastLog(), false);
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
		canComputeDoubleLog(new FastLog_log(), true);
	}

	@Test
	public void canComputeDoubleFastLog_fastLog()
	{
		canComputeDoubleLog(new FastLog_fastLog(), false);
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
	
	// TODO:
	// A robust error test using a uniform random number, report min,av,max,sd of the error.
	// Do this for different precisions.
	// Speed test of float/double version verses the Math.log. 
	
	
}
