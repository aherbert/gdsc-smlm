package gdsc.smlm.function;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.apache.commons.math3.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.utils.DoubleEquality;

public class ErfTest
{
	@Test
	public void erfxHasLowError()
	{
		RandomGenerator rg = new Well19937c(30051977);
		int range = 8;
		double max = 0;

		for (int xi = -range; xi <= range; xi++)
		{
			for (int i = 0; i < 5; i++)
			{
				double x = xi + rg.nextDouble();
				double o = Erf.erf(x);
				double e = org.apache.commons.math3.special.Erf.erf(x);
				double absError = Math.abs(o - e);
				if (absError < 1e-10)
					continue;
				double error = DoubleEquality.relativeError(o, e);
				if (max < error)
					max = error;
				//System.out.printf("x=%f, e=%f, o=%f, error=%f\n", x, e, o, error);
				Assert.assertTrue(error < 1.3e-4);
			}
		}
		System.out.printf("erfx max error = %f\n", max);
	}

	@Test
	public void erfxxHasLowError()
	{
		RandomGenerator rg = new Well19937c(30051977);

		int range = 3;
		double max = 0;

		for (int xi = -range; xi <= range; xi++)
		{
			for (int xi2 = -range; xi2 <= range; xi2++)
			{
				for (int i = 0; i < 5; i++)
				{
					double x = xi + rg.nextDouble();
					for (int j = 0; j < 5; j++)
					{
						double x2 = xi2 + rg.nextDouble();

						double o = Erf.erf(x, x2);
						double e = org.apache.commons.math3.special.Erf.erf(x, x2);
						double absError = Math.abs(o - e);
						if (absError < 1e-5)
							continue;
						double error = DoubleEquality.relativeError(o, e);
						if (max < error)
							max = error;
						//System.out.printf("x=%f, x2=%f, e=%f, o=%f, error=%f\n", x, x2, e, o, error);
						Assert.assertTrue(error < 0.04);
					}
				}
			}
		}

		System.out.printf("erfxx max error = %f\n", max);
	}

	@Test
	public void erfxxHasLowErrorForUnitBlocks()
	{
		int range = 8;
		double max = 0;

		for (int xi = -range; xi <= range; xi++)
		{
			double x = xi;
			double x2 = xi + 1;
			double o = Erf.erf(x, x2);
			double e = org.apache.commons.math3.special.Erf.erf(x, x2);
			double absError = Math.abs(o - e);
			if (absError < 1e-5)
				continue;
			double error = DoubleEquality.relativeError(o, e);
			if (max < error)
				max = error;
			//System.out.printf("x=%f, x2=%f, e=%f, o=%f, error=%f\n", x, x2, e, o, error);
			Assert.assertTrue(error < 0.04);
		}

		System.out.printf("erfxx unit max error = %f\n", max);
	}

	@Test
	public void erfxxHasLowerErrorThanGaussianApproximationForUnitBlocks()
	{
		int range = 5;
		double max = 0, max2 = 0;

		// Standard deviation
		double s = 1.3;
		final double twos2 = 2 * s * s;
		double norm = 1 / (Math.PI * twos2);
		final double denom = 1.0 / (Math.sqrt(2.0) * s);

		double sum1 = 0, sum2 = 0, sum3 = 0;

		for (int x = -range; x <= range; x++)
		{
			double o1 = 0.5 * Erf.erf((x - 0.5) * denom, (x + 0.5) * denom);
			double e1 = 0.5 * org.apache.commons.math3.special.Erf.erf((x - 0.5) * denom, (x + 0.5) * denom);
			for (int y = -range; y <= range; y++)
			{
				double o2 = 0.5 * Erf.erf((y - 0.5) * denom, (y + 0.5) * denom);
				double e2 = 0.5 * org.apache.commons.math3.special.Erf.erf((y - 0.5) * denom, (y + 0.5) * denom);

				double o = o1 * o2;
				double e = e1 * e2;
				double oo = norm * FastMath.exp(-(x * x + y * y) / twos2);

				sum1 += e;
				sum2 += o;
				sum3 += oo;

				double absError = Math.abs(o - e);
				if (e < 1e-4 || absError < 1e-10)
					continue;
				double error = DoubleEquality.relativeError(o, e);
				double error2 = DoubleEquality.relativeError(oo, e);
				if (max < error)
					max = error;
				if (max2 < error2)
					max2 = error2;
				//System.out.printf("x=%d, y=%d, e=%g, o=%g, o2=%g, error=%f, error2=%f\n", x, y, e, o, oo, error, error2);
				Assert.assertTrue(error < error2);
			}
		}

		Assert.assertTrue("Gaussian 2D integral is not 1", sum1 > 0.999);
		Assert.assertTrue("Erf approx integral is incorrect", DoubleEquality.relativeError(sum1, sum2) < 1e-3);
		Assert.assertTrue("Gaussian approx integral is incorrect", DoubleEquality.relativeError(sum1, sum3) < 1e-3);

		System.out.printf("Erf approx pixel unit max error = %f\n", max);
		System.out.printf("Gaussian approx pixel unit max error = %f\n", max2);
	}
}
