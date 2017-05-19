package gdsc.smlm.fitting.nonlinear;

import org.junit.Assert;
import org.junit.Test;

import gdsc.core.utils.BitFlags;

/**
 * Test the ToleranceChecker can converge as expected.
 */
public class ToleranceCheckerTest
{
	@Test(expected = IllegalArgumentException.class)
	public void throwsIfCannotConverge()
	{
		canConverge(false, 0, 0, 0, 0, 0, 0);
	}

	@Test
	public void canConvergeOnMaximumRelativeValue()
	{
		canConverge(false, 1e-2, 0, 0, 0, 0, ToleranceChecker.STATUS_VALUE);
	}
	
	@Test
	public void canConvergeOnMinimumRelativeValue()
	{
		canConverge(true, 1e-2, 0, 0, 0, 0, ToleranceChecker.STATUS_VALUE);
	}

	@Test
	public void canConvergeOnMaximumAbsoluteValue()
	{
		canConverge(false, 0, 1e-2, 0, 0, 0, ToleranceChecker.STATUS_VALUE);
	}

	@Test
	public void canConvergeOnMinimumAbsoluteValue()
	{
		canConverge(true, 0, 1e-2, 0, 0, 0, ToleranceChecker.STATUS_VALUE);
	}

	@Test(expected = AssertionError.class)
	public void cannotConvergeOnMaximumRelativeValueIfMinimising()
	{
		canConverge(false, -1, 1e-2, 0, 0, 0, 0, 0);
	}

	@Test(expected = AssertionError.class)
	public void cannotConvergeOnMaximumAbsoluteValueIfMinimising()
	{
		canConverge(false, -1, 0, 1e-2, 0, 0, 0, 0);
	}

	@Test(expected = AssertionError.class)
	public void cannotConvergeOnMinimumRelativeValueIfMaximising()
	{
		canConverge(true, 1, 1e-2, 0, 0, 0, 0, 0);
	}

	@Test(expected = AssertionError.class)
	public void cannotConvergeOnMinimumAbsoluteValueIfMaximising()
	{
		canConverge(true, 1, 0, 1e-2, 0, 0, 0, 0);
	}

	@Test
	public void canConvergeOnRelativeParameters()
	{
		canConverge(true, 0, 0, 1e-2, 0, 0, ToleranceChecker.STATUS_PARAMETERS);
		canConverge(false, 0, 0, 1e-2, 0, 0, ToleranceChecker.STATUS_PARAMETERS);
	}

	@Test
	public void canConvergeOnAbsoluteParameters()
	{
		canConverge(true, 0, 0, 0, 1e-2, 0, ToleranceChecker.STATUS_PARAMETERS);
		canConverge(false, 0, 0, 0, 1e-2, 0, ToleranceChecker.STATUS_PARAMETERS);
	}

	@Test
	public void canConvergeOnIterations()
	{
		canConverge(true, 0, 0, 0, 0, -10, ToleranceChecker.STATUS_TARGET_ITERATIONS);
		canConverge(false, 0, 0, 0, 0, -10, ToleranceChecker.STATUS_TARGET_ITERATIONS);
	}

	@Test
	public void canConvergeOnMaxIterations()
	{
		canConverge(true, 0, 0, 0, 0, 20, ToleranceChecker.STATUS_MAX_ITERATIONS);
		canConverge(false, 0, 0, 0, 0, 20, ToleranceChecker.STATUS_MAX_ITERATIONS);
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
}
