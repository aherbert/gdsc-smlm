package gdsc.smlm.fitting.nonlinear;

import org.junit.Assert;
import org.junit.Test;

import gdsc.core.utils.BitFlags;

/**
 * Test the ToleranceChecker can converge as expected.
 */
public class ToleranceCheckerTest
{
	double IGNORE = -1;

	@Test(expected = IllegalArgumentException.class)
	public void throwsIfCannotConverge()
	{
		canConverge(false, IGNORE, IGNORE, IGNORE, IGNORE, 0, 0);
	}

	@Test
	public void canConvergeOnMaximumRelativeValue()
	{
		canConverge(false, 1e-2, IGNORE, IGNORE, IGNORE, 0, ToleranceChecker.STATUS_VALUE);
	}

	@Test
	public void canConvergeOnMinimumRelativeValue()
	{
		canConverge(true, 1e-2, IGNORE, IGNORE, IGNORE, 0, ToleranceChecker.STATUS_VALUE);
	}

	@Test
	public void canConvergeOnMaximumAbsoluteValue()
	{
		canConverge(false, IGNORE, 1e-2, IGNORE, IGNORE, 0, ToleranceChecker.STATUS_VALUE);
	}

	@Test
	public void canConvergeOnMinimumAbsoluteValue()
	{
		canConverge(true, IGNORE, 1e-2, IGNORE, IGNORE, 0, ToleranceChecker.STATUS_VALUE);
	}

	@Test(expected = AssertionError.class)
	public void cannotConvergeOnMaximumRelativeValueIfMinimising()
	{
		canConverge(false, -1, 1e-2, IGNORE, IGNORE, IGNORE, 0, 0);
	}

	@Test(expected = AssertionError.class)
	public void cannotConvergeOnMaximumAbsoluteValueIfMinimising()
	{
		canConverge(false, -1, IGNORE, 1e-2, IGNORE, IGNORE, 0, 0);
	}

	@Test(expected = AssertionError.class)
	public void cannotConvergeOnMinimumRelativeValueIfMaximising()
	{
		canConverge(true, 1, 1e-2, IGNORE, IGNORE, IGNORE, 0, 0);
	}

	@Test(expected = AssertionError.class)
	public void cannotConvergeOnMinimumAbsoluteValueIfMaximising()
	{
		canConverge(true, 1, IGNORE, 1e-2, IGNORE, IGNORE, 0, 0);
	}

	@Test
	public void canConvergeOnRelativeParameters()
	{
		canConverge(true, IGNORE, IGNORE, 1e-2, IGNORE, 0, ToleranceChecker.STATUS_PARAMETERS);
		canConverge(false, IGNORE, IGNORE, 1e-2, IGNORE, 0, ToleranceChecker.STATUS_PARAMETERS);
	}

	@Test
	public void canConvergeOnAbsoluteParameters()
	{
		canConverge(true, IGNORE, IGNORE, IGNORE, 1e-2, 0, ToleranceChecker.STATUS_PARAMETERS);
		canConverge(false, IGNORE, IGNORE, IGNORE, 1e-2, 0, ToleranceChecker.STATUS_PARAMETERS);
	}

	@Test
	public void canConvergeOnIterations()
	{
		canConverge(true, IGNORE, IGNORE, IGNORE, IGNORE, -10, ToleranceChecker.STATUS_TARGET_ITERATIONS);
		canConverge(false, IGNORE, IGNORE, IGNORE, IGNORE, -10, ToleranceChecker.STATUS_TARGET_ITERATIONS);
	}

	@Test
	public void canConvergeOnMaxIterations()
	{
		canConverge(true, IGNORE, IGNORE, IGNORE, IGNORE, 20, ToleranceChecker.STATUS_MAX_ITERATIONS);
		canConverge(false, IGNORE, IGNORE, IGNORE, IGNORE, 20, ToleranceChecker.STATUS_MAX_ITERATIONS);
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
		ToleranceChecker tc = new ToleranceChecker(false, IGNORE, tolerance, IGNORE, IGNORE, 100);
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
		ToleranceChecker tc = new ToleranceChecker(true, IGNORE, tolerance, IGNORE, IGNORE, 100);
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
		
		tc = new ToleranceChecker(false, IGNORE, IGNORE, IGNORE, IGNORE, 100);
		Assert.assertFalse(tc.checkValue);
		Assert.assertFalse(tc.checkParameters);
		
		tc = new ToleranceChecker(true, tolerance, IGNORE, IGNORE, IGNORE, 100);
		Assert.assertTrue(tc.checkValue);
		Assert.assertFalse(tc.checkParameters);
		Assert.assertEquals(0, tc.converged(0, null, -Double.MIN_VALUE, null));
		Assert.assertEquals(0, tc.converged(0, null, Double.MIN_VALUE, null));
		Assert.assertEquals(ToleranceChecker.STATUS_VALUE, tc.converged(0, null, 0, null));
		
		tc = new ToleranceChecker(false, tolerance, IGNORE, IGNORE, IGNORE, 100);
		Assert.assertTrue(tc.checkValue);
		Assert.assertFalse(tc.checkParameters);
		Assert.assertEquals(0, tc.converged(0, null, -Double.MIN_VALUE, null));
		Assert.assertEquals(0, tc.converged(0, null, Double.MIN_VALUE, null));
		Assert.assertEquals(ToleranceChecker.STATUS_VALUE, tc.converged(0, null, 0, null));
		
		tc = new ToleranceChecker(true, IGNORE, tolerance, IGNORE, IGNORE, 100);
		Assert.assertTrue(tc.checkValue);
		Assert.assertFalse(tc.checkParameters);
		Assert.assertEquals(0, tc.converged(0, null, -Double.MIN_VALUE, null));
		Assert.assertEquals(0, tc.converged(0, null, Double.MIN_VALUE, null));
		Assert.assertEquals(ToleranceChecker.STATUS_VALUE, tc.converged(0, null, 0, null));
		
		tc = new ToleranceChecker(false, IGNORE, tolerance, IGNORE, IGNORE, 100);
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
		
		tc = new ToleranceChecker(false, IGNORE, IGNORE, IGNORE, IGNORE, 100);
		Assert.assertFalse(tc.checkValue);
		Assert.assertFalse(tc.checkParameters);
		
		tc = new ToleranceChecker(true, IGNORE, IGNORE, tolerance, IGNORE, 100);
		Assert.assertFalse(tc.checkValue);
		Assert.assertTrue(tc.checkParameters);
		Assert.assertEquals(0, tc.converged(0, p, 0, new double[]{ -Double.MIN_VALUE}));
		Assert.assertEquals(0, tc.converged(0, p, 0, new double[]{ Double.MIN_VALUE}));
		Assert.assertEquals(ToleranceChecker.STATUS_PARAMETERS, tc.converged(0, p, 0, p));
		
		tc = new ToleranceChecker(true, IGNORE, IGNORE, IGNORE, tolerance, 100);
		Assert.assertFalse(tc.checkValue);
		Assert.assertTrue(tc.checkParameters);
		Assert.assertEquals(0, tc.converged(0, p, 0, new double[]{ -Double.MIN_VALUE}));
		Assert.assertEquals(0, tc.converged(0, p, 0, new double[]{ Double.MIN_VALUE}));
		Assert.assertEquals(ToleranceChecker.STATUS_PARAMETERS, tc.converged(0, p, 0, p));
	}
}
