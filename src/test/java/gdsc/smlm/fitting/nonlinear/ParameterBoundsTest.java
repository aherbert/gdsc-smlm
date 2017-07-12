package gdsc.smlm.fitting.nonlinear;

import org.junit.Assert;
import org.junit.Test;

import gdsc.smlm.function.FakeGradientFunction;

public class ParameterBoundsTest
{
	@Test
	public void canUpperBoundParameter()
	{
		canBoundParameter(1);
	}

	@Test
	public void canLowerBoundParameter()
	{
		canBoundParameter(-1);
	}

	private void canBoundParameter(double s)
	{
		ParameterBounds b = new ParameterBounds(new FakeGradientFunction(1, 1, 1, 1, 1));
		if (s < 0)
		{
			b.setBounds(new double[] { 2 * s }, null);
		}
		else
		{
			b.setBounds(null, new double[] { 2 * s });
		}
		double[] a1 = new double[1];
		double[] a2 = new double[1];
		double[] step = new double[] { s };
		b.applyBounds(a1, step, a2);
		Assert.assertArrayEquals("Step 1", a2, new double[] { 1 * s }, 0);
		b.applyBounds(a2, step, a1);
		Assert.assertArrayEquals("Step 2", a1, new double[] { 2 * s }, 0);
		b.applyBounds(a1, step, a2);
		// Should be bounded
		Assert.assertArrayEquals("Step 3", a2, new double[] { 2 * s }, 0);
	}

	@Test
	public void canDoubleBoundParameter()
	{
		ParameterBounds b = new ParameterBounds(new FakeGradientFunction(1, 1, 1, 1, 1));
		double s = 2;
		b.setBounds(new double[] { -s }, new double[] { s });
		double[] a1 = new double[1];
		double[] a2 = new double[1];
		b.applyBounds(a1, new double[] { 10 }, a2);
		Assert.assertArrayEquals("Step 10", a2, new double[] { s }, 0);

		b.applyBounds(a1, new double[] { -10 }, a2);
		Assert.assertArrayEquals("Step -10", a2, new double[] { -s }, 0);
	}

	@Test
	public void canStepParameter()
	{
		canStepParameter(0.5);
		canStepParameter(1);
		canStepParameter(-0.5);
		canStepParameter(-1);
	}

	private void canStepParameter(double s)
	{
		ParameterBounds b = new ParameterBounds(new FakeGradientFunction(1, 1, 1, 1, 1));
		double[] a1 = new double[1];
		double[] a2 = new double[1];
		double[] tmp;
		double[] step = new double[] { s };
		for (int i = 1; i <= 10; i++)
		{
			b.applyBounds(a1, step, a2);
			Assert.assertArrayEquals("Step " + i, a2, new double[] { i * s }, 0);
			tmp = a1;
			a1 = a2;
			a2 = tmp;
		}
	}
}
