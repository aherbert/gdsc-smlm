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
package uk.ac.sussex.gdsc.smlm.function;

import org.apache.commons.math3.util.FastMath;
import org.junit.jupiter.api.Test;import uk.ac.sussex.gdsc.test.junit5.SeededTest;import uk.ac.sussex.gdsc.test.junit5.RandomSeed;import uk.ac.sussex.gdsc.test.junit5.SpeedTag;

import uk.ac.sussex.gdsc.test.BaseTimingTask;
import uk.ac.sussex.gdsc.test.LogLevel;
import uk.ac.sussex.gdsc.test.TestLog;
import uk.ac.sussex.gdsc.test.TestSettings;
import uk.ac.sussex.gdsc.test.TimingService;
import uk.ac.sussex.gdsc.test.junit5.ExtraAssumptions;

@SuppressWarnings({ "javadoc" })
public class FastMathTest
{
	//@formatter:off
	private abstract class FunctionTimingTask extends BaseTimingTask
	{
		double[] data;
		FunctionTimingTask(String name, double[] data) { super(name); this.data = data; }
		@Override
		public int getSize() { return 1; }
		@Override
		public Object getData(int i) { return null;	}
		@Override
		public Object run(Object o)
		{
			for (int i=0; i<data.length; i++)
				value(data[i]);
			return null;
		}
		abstract double value(double d);
	}
	private class MathPow1_3 extends FunctionTimingTask
	{
		final double THIRD = 1.0/3.0;
		MathPow1_3(double[] data) { super("Math pow 1/3", data); }
		@Override
		double value(double d)
		{
			return Math.pow(d, THIRD);
		}
	}
	private class FastMathPow1_3 extends FunctionTimingTask
	{
		final double THIRD = 1.0/3.0;
		FastMathPow1_3(double[] data) { super("FastMath pow 1/3", data); }
		@Override
		double value(double d)
		{
			return FastMath.pow(d, THIRD);
		}
	}
	private class FastMathCbrt extends FunctionTimingTask
	{
		FastMathCbrt(double[] data) { super("FastMath cbrt", data); }
		@Override
		double value(double d)
		{
			return FastMath.cbrt(d);
		}
	}
	//@formatter:on

	@Test
	public void cbrtIsFaster()
	{
		ExtraAssumptions.assumeSpeedTest();

		// Q. What is a suitable range for this test?
		final int range = 5;
		final int steps = 10000;
		final double[] x = new double[steps];
		final double total = 2 * range;
		final double step = total / steps;
		for (int i = 0; i < steps; i++)
			x[i] = -range + i * step;

		final TimingService ts = new TimingService(5);
		ts.execute(new MathPow1_3(x));
		ts.execute(new FastMathPow1_3(x));
		ts.execute(new FastMathCbrt(x));

		final int size = ts.getSize();
		ts.repeat(size);
		if (TestSettings.allow(LogLevel.INFO))
			ts.report();

		for (int k = 2; k <= 3; k++)
		{
			final double t1 = ts.get(-1).getMean();
			final double t2 = ts.get(-k).getMean();
			TestLog.logSpeedTestResult(t1 < t2, "%s %s => %s %s = %.2fx\n", ts.get(-k).getTask().getName(), t2,
					ts.get(-1).getTask().getName(), t1, t2 / t1);
		}
	}
}
