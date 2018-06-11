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
package gdsc.smlm.fitting.nonlinear;

import org.junit.Assert;
import org.junit.Test;

import gdsc.core.utils.BitFlags;

/**
 * Test the ToleranceChecker can converge as expected.
 */
public class ToleranceCheckerTest
{
	double NONE = ToleranceChecker.IGNORE_TOLERANCE;
	int IGNORE = ToleranceChecker.IGNORE_MAX_ITERATIONS;

	@Test(expected = IllegalArgumentException.class)
	public void throwsIfCannotConverge()
	{
		canConverge(false, NONE, NONE, NONE, NONE, IGNORE, 0);
	}

	@Test
	public void canConvergeOnMaximumRelativeValue()
	{
		canConverge(false, 1e-2, NONE, NONE, NONE, IGNORE, ToleranceChecker.STATUS_VALUE);
	}

	@Test
	public void canConvergeOnMinimumRelativeValue()
	{
		canConverge(true, 1e-2, NONE, NONE, NONE, IGNORE, ToleranceChecker.STATUS_VALUE);
	}

	@Test
	public void canConvergeOnMaximumAbsoluteValue()
	{
		canConverge(false, NONE, 1e-2, NONE, NONE, IGNORE, ToleranceChecker.STATUS_VALUE);
	}

	@Test
	public void canConvergeOnMinimumAbsoluteValue()
	{
		canConverge(true, NONE, 1e-2, NONE, NONE, IGNORE, ToleranceChecker.STATUS_VALUE);
	}

	@Test(expected = AssertionError.class)
	public void cannotConvergeOnMaximumRelativeValueIfMinimising()
	{
		canConverge(false, -1, 1e-2, NONE, NONE, NONE, IGNORE, 0);
	}

	@Test(expected = AssertionError.class)
	public void cannotConvergeOnMaximumAbsoluteValueIfMinimising()
	{
		canConverge(false, -1, NONE, 1e-2, NONE, NONE, IGNORE, 0);
	}

	@Test(expected = AssertionError.class)
	public void cannotConvergeOnMinimumRelativeValueIfMaximising()
	{
		canConverge(true, 1, 1e-2, NONE, NONE, NONE, IGNORE, 0);
	}

	@Test(expected = AssertionError.class)
	public void cannotConvergeOnMinimumAbsoluteValueIfMaximising()
	{
		canConverge(true, 1, NONE, 1e-2, NONE, NONE, IGNORE, 0);
	}

	@Test
	public void canConvergeOnRelativeParameters()
	{
		canConverge(true, NONE, NONE, 1e-2, NONE, IGNORE, ToleranceChecker.STATUS_PARAMETERS);
		canConverge(false, NONE, NONE, 1e-2, NONE, IGNORE, ToleranceChecker.STATUS_PARAMETERS);
	}

	@Test
	public void canConvergeOnAbsoluteParameters()
	{
		canConverge(true, NONE, NONE, NONE, 1e-2, IGNORE, ToleranceChecker.STATUS_PARAMETERS);
		canConverge(false, NONE, NONE, NONE, 1e-2, IGNORE, ToleranceChecker.STATUS_PARAMETERS);
	}

	@Test
	public void canConvergeOnIterations()
	{
		canConverge(true, NONE, NONE, NONE, NONE, -10, ToleranceChecker.STATUS_TARGET_ITERATIONS);
		canConverge(false, NONE, NONE, NONE, NONE, -10, ToleranceChecker.STATUS_TARGET_ITERATIONS);
	}

	@Test
	public void canConvergeOnMaxIterations()
	{
		canConverge(true, NONE, NONE, NONE, NONE, 20, ToleranceChecker.STATUS_MAX_ITERATIONS);
		canConverge(false, NONE, NONE, NONE, NONE, 20, ToleranceChecker.STATUS_MAX_ITERATIONS);
	}

	private void canConverge(boolean minimiseValue, double relativeValue, double absoluteValue,
			double relativeParameters, double absoluteParameters, int maxIterations, int expected)
	{
		double dir = (minimiseValue) ? -1 : 1;
		canConverge(minimiseValue, dir, relativeValue, absoluteValue, relativeParameters, absoluteParameters,
				maxIterations, expected);
	}

	private void canConverge(boolean minimiseValue, double dir, double relativeValue, double absoluteValue,
			double relativeParameters, double absoluteParameters, int maxIterations, int expected)
	{
		ToleranceChecker tc = new ToleranceChecker(minimiseValue, relativeValue, absoluteValue, relativeParameters,
				absoluteParameters, maxIterations);

		double v = 10;
		double v2 = 1 + v;
		double p = 20;
		double[] p2 = new double[] { 1 + p };
		for (int i = 0; i < 20; i++)
		{
			double v1 = v2;
			double[] p1 = p2;
			v *= 0.5;
			p *= 0.5;
			v2 = v1 + dir * v;
			//System.out.printf("v2 = %f\n", v2);
			p2 = new double[] { p1[0] + dir * p };
			int observed = tc.converged(v1, p1, v2, p2);
			if (observed != 0)
			{
				Assert.assertEquals(expected, observed);
				if (BitFlags.areSet(expected, ToleranceChecker.STATUS_TARGET_ITERATIONS))
					Assert.assertEquals(-maxIterations, tc.getIterations());
				if (BitFlags.areSet(expected, ToleranceChecker.STATUS_MAX_ITERATIONS))
					Assert.assertEquals(maxIterations, tc.getIterations());
				return;
			}
		}

		Assert.fail();
	}

	@Test
	public void canConvergeOnImprovedValueIfMaximising()
	{
		double tolerance = 1e-2;
		ToleranceChecker tc = new ToleranceChecker(false, NONE, tolerance, NONE, NONE, 100);
		Assert.assertEquals(0, tc.converged(0, null, 1, null));
		Assert.assertEquals(0, tc.converged(0, null, -1, null));
		Assert.assertEquals(0, tc.converged(0, null, 2 * tolerance, null));
		Assert.assertEquals(0, tc.converged(0, null, -2 * tolerance, null));
		Assert.assertEquals(ToleranceChecker.STATUS_VALUE, tc.converged(0, null, tolerance, null));
		Assert.assertEquals(0, tc.converged(0, null, -tolerance, null));
		Assert.assertEquals(ToleranceChecker.STATUS_VALUE, tc.converged(0, null, 0, null));
	}

	@Test
	public void canConvergeOnImprovedValueIfMinimising()
	{
		double tolerance = 1e-2;
		ToleranceChecker tc = new ToleranceChecker(true, NONE, tolerance, NONE, NONE, 100);
		Assert.assertEquals(0, tc.converged(0, null, 1, null));
		Assert.assertEquals(0, tc.converged(0, null, -1, null));
		Assert.assertEquals(0, tc.converged(0, null, 2 * tolerance, null));
		Assert.assertEquals(0, tc.converged(0, null, -2 * tolerance, null));
		Assert.assertEquals(0, tc.converged(0, null, tolerance, null));
		Assert.assertEquals(ToleranceChecker.STATUS_VALUE, tc.converged(0, null, -tolerance, null));
		Assert.assertEquals(ToleranceChecker.STATUS_VALUE, tc.converged(0, null, 0, null));
	}

	@Test
	public void canConvergeOnValueUsingZeroTolerance()
	{
		double tolerance = 0;
		ToleranceChecker tc;

		tc = new ToleranceChecker(false, NONE, NONE, NONE, NONE, 100);
		Assert.assertFalse(tc.checkValue);
		Assert.assertFalse(tc.checkParameters);

		tc = new ToleranceChecker(true, tolerance, NONE, NONE, NONE, 100);
		Assert.assertTrue(tc.checkValue);
		Assert.assertFalse(tc.checkParameters);
		Assert.assertEquals(0, tc.converged(0, null, -Double.MIN_VALUE, null));
		Assert.assertEquals(0, tc.converged(0, null, Double.MIN_VALUE, null));
		Assert.assertEquals(ToleranceChecker.STATUS_VALUE, tc.converged(0, null, 0, null));

		tc = new ToleranceChecker(false, tolerance, NONE, NONE, NONE, 100);
		Assert.assertTrue(tc.checkValue);
		Assert.assertFalse(tc.checkParameters);
		Assert.assertEquals(0, tc.converged(0, null, -Double.MIN_VALUE, null));
		Assert.assertEquals(0, tc.converged(0, null, Double.MIN_VALUE, null));
		Assert.assertEquals(ToleranceChecker.STATUS_VALUE, tc.converged(0, null, 0, null));

		tc = new ToleranceChecker(true, NONE, tolerance, NONE, NONE, 100);
		Assert.assertTrue(tc.checkValue);
		Assert.assertFalse(tc.checkParameters);
		Assert.assertEquals(0, tc.converged(0, null, -Double.MIN_VALUE, null));
		Assert.assertEquals(0, tc.converged(0, null, Double.MIN_VALUE, null));
		Assert.assertEquals(ToleranceChecker.STATUS_VALUE, tc.converged(0, null, 0, null));

		tc = new ToleranceChecker(false, NONE, tolerance, NONE, NONE, 100);
		Assert.assertTrue(tc.checkValue);
		Assert.assertFalse(tc.checkParameters);
		Assert.assertEquals(0, tc.converged(0, null, -Double.MIN_VALUE, null));
		Assert.assertEquals(0, tc.converged(0, null, Double.MIN_VALUE, null));
		Assert.assertEquals(ToleranceChecker.STATUS_VALUE, tc.converged(0, null, 0, null));
	}

	@Test
	public void canConvergeOnParametersUsingZeroTolerance()
	{
		double tolerance = 0;
		ToleranceChecker tc;
		double[] p = new double[1];

		tc = new ToleranceChecker(false, NONE, NONE, NONE, NONE, 100);
		Assert.assertFalse(tc.checkValue);
		Assert.assertFalse(tc.checkParameters);

		tc = new ToleranceChecker(true, NONE, NONE, tolerance, NONE, 100);
		Assert.assertFalse(tc.checkValue);
		Assert.assertTrue(tc.checkParameters);
		Assert.assertEquals(0, tc.converged(0, p, 0, new double[] { -Double.MIN_VALUE }));
		Assert.assertEquals(0, tc.converged(0, p, 0, new double[] { Double.MIN_VALUE }));
		Assert.assertEquals(ToleranceChecker.STATUS_PARAMETERS, tc.converged(0, p, 0, p));

		tc = new ToleranceChecker(true, NONE, NONE, NONE, tolerance, 100);
		Assert.assertFalse(tc.checkValue);
		Assert.assertTrue(tc.checkParameters);
		Assert.assertEquals(0, tc.converged(0, p, 0, new double[] { -Double.MIN_VALUE }));
		Assert.assertEquals(0, tc.converged(0, p, 0, new double[] { Double.MIN_VALUE }));
		Assert.assertEquals(ToleranceChecker.STATUS_PARAMETERS, tc.converged(0, p, 0, p));
	}
}
