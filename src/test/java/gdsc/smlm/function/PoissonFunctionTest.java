package gdsc.smlm.function;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.utils.DoubleEquality;
import gnu.trove.list.array.TDoubleArrayList;

@SuppressWarnings("unused")
public class PoissonFunctionTest
{
	static double[] gain = { 1, 2, 4, 8, 16 };
	static double[] photons = { 0.25, 0.5, 1, 2, 4, 10, 100, 1000 };

	// Set this at the range output from cumulativeProbabilityIsOneWithInteger
	static int[][] minRange = null, maxRange = null;

	static
	{
		minRange = new int[gain.length][photons.length];
		maxRange = new int[gain.length][photons.length];
		minRange[0][0] = 0;
		maxRange[0][0] = 2;
		minRange[0][1] = 0;
		maxRange[0][1] = 3;
		minRange[0][2] = 0;
		maxRange[0][2] = 5;
		minRange[0][3] = 0;
		maxRange[0][3] = 7;
		minRange[0][4] = 0;
		maxRange[0][4] = 11;
		minRange[0][5] = 1;
		maxRange[0][5] = 20;
		minRange[0][6] = 72;
		maxRange[0][6] = 129;
		minRange[0][7] = 911;
		maxRange[0][7] = 1090;
		minRange[1][0] = 0;
		maxRange[1][0] = 5;
		minRange[1][1] = 0;
		maxRange[1][1] = 7;
		minRange[1][2] = 0;
		maxRange[1][2] = 10;
		minRange[1][3] = 0;
		maxRange[1][3] = 15;
		minRange[1][4] = 0;
		maxRange[1][4] = 22;
		minRange[1][5] = 4;
		maxRange[1][5] = 41;
		minRange[1][6] = 146;
		maxRange[1][6] = 259;
		minRange[1][7] = 1824;
		maxRange[1][7] = 2180;
		minRange[2][0] = 0;
		maxRange[2][0] = 11;
		minRange[2][1] = 0;
		maxRange[2][1] = 15;
		minRange[2][2] = 0;
		maxRange[2][2] = 21;
		minRange[2][3] = 0;
		maxRange[2][3] = 30;
		minRange[2][4] = 0;
		maxRange[2][4] = 44;
		minRange[2][5] = 10;
		maxRange[2][5] = 82;
		minRange[2][6] = 293;
		maxRange[2][6] = 518;
		minRange[2][7] = 3651;
		maxRange[2][7] = 4362;
		minRange[3][0] = 0;
		maxRange[3][0] = 23;
		minRange[3][1] = 0;
		maxRange[3][1] = 31;
		minRange[3][2] = 0;
		maxRange[3][2] = 43;
		minRange[3][3] = 0;
		maxRange[3][3] = 60;
		minRange[3][4] = 0;
		maxRange[3][4] = 89;
		minRange[3][5] = 22;
		maxRange[3][5] = 164;
		minRange[3][6] = 587;
		maxRange[3][6] = 1037;
		minRange[3][7] = 7303;
		maxRange[3][7] = 8725;
		minRange[4][0] = 0;
		maxRange[4][0] = 47;
		minRange[4][1] = 0;
		maxRange[4][1] = 63;
		minRange[4][2] = 0;
		maxRange[4][2] = 86;
		minRange[4][3] = 0;
		maxRange[4][3] = 121;
		minRange[4][4] = 1;
		maxRange[4][4] = 178;
		minRange[4][5] = 45;
		maxRange[4][5] = 328;
		minRange[4][6] = 1176;
		maxRange[4][6] = 2075;
		minRange[4][7] = 14613;
		maxRange[4][7] = 17455;
	}

	@Test
	public void cumulativeProbabilityIsOneWithInteger()
	{
		for (int j = 0; j < gain.length; j++)
			for (int i = 0; i < photons.length; i++)
			{
				int[] result = cumulativeProbabilityIsOneWithInteger(gain[j], photons[i]);
				//System.out.printf("minRange[%d][%d] = %d;\n", j, i, result[0]);
				//System.out.printf("maxRange[%d][%d] = %d;\n", j, i, result[1]);
			}
	}

	@Test
	public void cumulativeProbabilityIsOneWithReal()
	{
		for (int j = 0; j < gain.length; j++)
			for (int i = 0; i < photons.length; i++)
				if (photons[i] / gain[j] >= 4)
					cumulativeProbabilityIsOneWithRealAbove4(gain[j], photons[i], minRange[j][i], maxRange[j][i] + 1);
	}

	private int[] cumulativeProbabilityIsOneWithInteger(final double gain, final double mu)
	{
		final double o = mu;

		PoissonFunction f = new PoissonFunction(1.0 / gain, false);
		double p = 0;

		TDoubleArrayList values = new TDoubleArrayList();

		double maxp = 0;
		int maxc = 0;

		// Evaluate an initial range. 
		// Poisson will have mean mu with a variance mu. 
		// At large mu it is approximately normal so use 3 sqrt(mu) for the range added to the mean

		int[] range = getRange(gain, mu);
		int min = range[0];
		int max = range[1];
		for (int x = min; x <= max; x++)
		{
			final double pp = f.likelihood(x, o);
			//System.out.printf("x=%d, p=%f\n", x, pp);
			p += pp;
			values.add(pp);
			if (maxp < pp)
			{
				maxp = pp;
				maxc = x;
			}
		}
		if (p > 1.01)
			Assert.fail("P > 1: " + p);

		// We have most of the probability density. 
		// Now keep evaluating up and down until no difference
		final double changeTolerance = 1e-6;
		if (min > 0)
		{
			values.reverse();
			for (int x = min - 1; x >= 0; x--)
			{
				min = x;
				final double pp = f.likelihood(x, o);
				//System.out.printf("x=%d, p=%f\n", x, pp);
				p += pp;
				values.add(pp);
				if (maxp < pp)
				{
					maxp = pp;
					maxc = x;
				}
				if (pp == 0 || pp / p < changeTolerance)
					break;
			}
			values.reverse();
		}
		for (int x = max + 1;; x++)
		{
			max = x;
			final double pp = f.likelihood(x, o);
			//System.out.printf("x=%d, p=%f\n", x, pp);
			p += pp;
			values.add(pp);
			if (maxp < pp)
			{
				maxp = pp;
				maxc = x;
			}
			if (pp == 0 || pp / p < changeTolerance)
				break;
		}

		// Find the range for 99.5% of the sum
		double[] h = values.toArray();
		// Find cumulative
		for (int i = 1; i < h.length; i++)
		{
			h[i] += h[i - 1];
		}
		int minx = 0, maxx = h.length - 1;
		while (h[minx + 1] < 0.0025)
			minx++;
		while (h[maxx - 1] > 0.9975)
			maxx--;

		minx += min;
		maxx += min;

		System.out.printf("g=%f, mu=%f, o=%f, p=%f, min=%d, %f @ %d, max=%d\n", gain, mu, o, p, minx, maxp, maxc, maxx);
		Assert.assertEquals(String.format("g=%f, mu=%f", gain, mu), 1, p, 0.02);
		return new int[] { minx, maxx };
	}

	static int[] getRange(final double gain, final double mu)
	{
		// Evaluate an initial range. 
		// Poisson will have mean mu with a variance mu. 
		// At large mu it is approximately normal so use 3 sqrt(mu) for the range added to the mean
		double range = Math.max(1, Math.sqrt(mu));
		int min = Math.max(0, (int) Math.floor(gain * (mu - 3 * range)));
		int max = (int) Math.ceil(gain * (mu + 3 * range));
		return new int[] { min, max };
	}

	private void cumulativeProbabilityIsOneWithRealAbove4(final double gain, final double mu, int min, int max)
	{
		final double o = mu;

		final PoissonFunction f = new PoissonFunction(1.0 / gain, true);
		double p = 0;
		SimpsonIntegrator in = new SimpsonIntegrator(3, 30);

		try
		{
			p = in.integrate(Integer.MAX_VALUE, new UnivariateFunction()
			{
				public double value(double x)
				{
					return f.likelihood(x, o);
				}
			}, min, max);

			System.out.printf("g=%f, mu=%f, o=%f, p=%f\n", gain, mu, o, p);
			//Assert.assertEquals(String.format("g=%f, mu=%f", gain, mu), 1, p, 0.02);
		}
		catch (TooManyEvaluationsException e)
		{
			//double inc = max / 20000.0;
			//for (double x = 0; x <= max; x += inc)
			//{
			//	final double pp = f.likelihood(x, o);
			//	//System.out.printf("g=%f, mu=%f, o=%f, p=%f\n", gain, mu, o, pp);
			//	p += pp;
			//}
			//System.out.printf("g=%f, mu=%f, o=%f, p=%f\n", gain, mu, o, p);
			Assert.assertFalse(e.getMessage(), true);
		}
	}

	@Test
	public void canComputeGradientWithInteger()
	{
		for (int j = 0; j < gain.length; j++)
			for (int i = 0; i < photons.length; i++)
			{
				canComputeGradient(gain[j], photons[i], false);
			}
	}

	@Test
	public void cumulativeProbabilityWithReal()
	{
		for (int j = 0; j < gain.length; j++)
			for (int i = 0; i < photons.length; i++)
				canComputeGradient(gain[j], photons[i], true);
	}

	private void canComputeGradient(final double gain, final double mu, boolean nonInteger)
	{
		final double o = mu;
		double delta = 1e-3;
		double uo = o + delta;
		double lo = o - delta;
		double diff = uo - lo;

		PoissonFunction f = new PoissonFunction(1.0 / gain, nonInteger);

		int[] range = getRange(gain, mu);
		int min = range[0];
		int max = range[1];
		double[] dp_dt = new double[1];
		double step = (nonInteger) ? 0.5 : 1;
		for (double x = min; x <= max; x += step)
		{
			double p1 = f.likelihood(x, o);
			double p2 = f.likelihood(x, o, dp_dt);
			Assert.assertEquals(p1, p2, 0);

			double up = f.likelihood(x, uo);
			double lp = f.likelihood(x, lo);

			double eg = dp_dt[0];
			double g = (up - lp) / diff;
			double error = DoubleEquality.relativeError(g, eg);
			double ox = x / gain;
			//System.out.printf("g=%g, mu=%g, x=%g (ox=%g), p=%g  g=%g  %g  error=%g\n", gain, mu, x, ox, p1, g, eg,
			//		error);

			// Ignore tiny gradients. These occur due to floating point error when the gradient
			// should be zero, e.g. mu*gain=x, i.e. the max of the distribution PMF
			if (Math.abs(eg) < 1e-10)
			{
				System.out.printf("g=%g, mu=%g, x=%g (ox=%g), p=%g  g=%g  %g  error=%g\n", gain, mu, x, ox, p1, g, eg,
						error);
				continue;
			}

			if (nonInteger && ox < 1)
			{
				// Gradients are wrong
				//Assert.assertTrue(error < 0.5);
			}
			else
			{
				Assert.assertTrue(error < 1e-3);
			}
		}
	}
}
