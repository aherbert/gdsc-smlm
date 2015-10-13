package gdsc.smlm.fitting;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.junit.Assert;
import org.junit.Test;

public class JumpDistanceAnalysisTest
{
	double delta = 1e-1;
	double[] D = new double[] { 0.1, 1, 10 };
	RandomGenerator random = new Well19937c(System.currentTimeMillis() + System.identityHashCode(this));

	@Test
	public void canFitSinglePopulation100Samples()
	{
		for (double d : D)
		{
			fit(100, 0, new double[] { d }, new double[] { 1 });
		}
	}

	@Test
	public void canFitSinglePopulation500Samples()
	{
		for (double d : D)
		{
			fit(500, 0, new double[] { d }, new double[] { 1 });
		}
	}

	@Test
	public void canFitSinglePopulation2000Samples()
	{
		for (double d : D)
		{
			fit(2000, 0, new double[] { d }, new double[] { 1 });
		}
	}

	@Test
	public void canFit11DualPopulation100Samples()
	{
		for (int i = 0; i < D.length; i++)
		{
			for (int j = i + 1; j < D.length; j++)
			{
				fit(100, 0, new double[] { D[i], D[j] }, new double[] { 1, 1 });
			}
		}
	}

	@Test
	public void canFit11DualPopulation500Samples()
	{
		for (int i = 0; i < D.length; i++)
		{
			for (int j = i + 1; j < D.length; j++)
			{
				fit(500, 0, new double[] { D[i], D[j] }, new double[] { 1, 1 });
			}
		}
	}

	@Test
	public void canFit11DualPopulation2000Samples()
	{
		for (int i = 0; i < D.length; i++)
		{
			for (int j = i + 1; j < D.length; j++)
			{
				fit(2000, 0, new double[] { D[i], D[j] }, new double[] { 1, 1 });
			}
		}
	}

	@Test
	public void canFit12DualPopulation100Samples()
	{
		for (int i = 0; i < D.length; i++)
		{
			for (int j = i + 1; j < D.length; j++)
			{
				fit(100, 0, new double[] { D[i], D[j] }, new double[] { 1, 2 });
				fit(100, 0, new double[] { D[i], D[j] }, new double[] { 2, 1 });
			}
		}
	}

	@Test
	public void canFit12DualPopulation500Samples()
	{
		for (int i = 0; i < D.length; i++)
		{
			for (int j = i + 1; j < D.length; j++)
			{
				fit(500, 0, new double[] { D[i], D[j] }, new double[] { 1, 2 });
				fit(500, 0, new double[] { D[i], D[j] }, new double[] { 2, 1 });
			}
		}
	}

	@Test
	public void canFit12DualPopulation2000Samples()
	{
		for (int i = 0; i < D.length; i++)
		{
			for (int j = i + 1; j < D.length; j++)
			{
				fit(2000, 0, new double[] { D[i], D[j] }, new double[] { 1, 2 });
				fit(2000, 0, new double[] { D[i], D[j] }, new double[] { 2, 1 });
			}
		}
	}
	
	private void fit(int samples, int n, double[] d, double[] f)
	{
		JumpDistanceAnalysis.sort(d, f);
		double[] jumpsDistances = createData(samples, d, f);
		JumpDistanceAnalysis jd = new JumpDistanceAnalysis();
		jd.setFitRestarts(5);
		jd.setMinFraction(0.1);
		jd.setMinDifference(2);
		jd.setN((n > 0) ? n : 10);
		double[][] fit = jd.fitJumpDistances(jumpsDistances);
		double[] fitD = fit[0];
		double[] fitF = fit[1];
		AssertionError error = null;
		try
		{
			Assert.assertEquals("Failed to fit n", d.length, fitD.length);
			for (int i = 0; i < d.length; i++)
			{
				Assert.assertEquals("Failed to fit d", d[i], fitD[i], delta * d[i]);
				Assert.assertEquals("Failed to fit f", f[i], fitF[i], delta * f[i]);
			}
		}
		catch (AssertionError e)
		{
			error = e;
		}
		finally
		{
			log("%s sample=%d, n<=%d : %s = %s : %s = %s\n", (error == null) ? "Pass" : "Fail", samples, n,
					toString(d), toString(fitD), toString(f), toString(fitF));
			if (error != null)
				throw error;
		}
	}

	private String toString(double[] d)
	{
		if (d.length == 0)
			return "";
		if (d.length == 1)
			return format(d[0]);
		StringBuilder sb = new StringBuilder();
		sb.append(format(d[0]));
		for (int i = 1; i < d.length; i++)
			sb.append(',').append(format(d[i]));
		return sb.toString();
	}

	private String format(double d)
	{
		return String.format("%.3f", d);
	}

	/**
	 * Create random jump distances
	 * 
	 * @param n
	 *            Number of jump distances
	 * @param d
	 *            Diffusion rate (should be ascending order of magnitude)
	 * @param f
	 *            Fraction of population (will be updated with the actual fractions, normalised to sum to 1)
	 * @return The jump distances
	 */
	private double[] createData(int n, double[] d, double[] f)
	{
		// Convert diffusion co-efficient into the standard deviation for the random move in each dimension
		// For 1D diffusion: sigma^2 = 2D
		//                   sigma = sqrt(2D)
		d = d.clone();
		double sum = 0;
		for (int i = 0; i < d.length; i++)
		{
			d[i] = Math.sqrt(2 * d[i]);
			sum += f[i];
		}
		// Normalise fractions to 1
		for (int i = 0; i < f.length; i++)
			f[i] /= sum;

		// Get cumulative fraction
		double[] c = f.clone();
		sum = 0;
		for (int i = 0; i < f.length; i++)
		{
			sum += f[i];
			c[i] = sum;
			// Reset to count the actual fractions
			f[i] = 0;
		}

		double[] data = new double[n];
		for (int i = 0; i < data.length; i++)
		{
			// Pick the population using the fraction
			final int j = pick(c, random.nextDouble());
			f[j]++;
			// Get the x/y shifts
			final double x = random.nextGaussian() * d[j];
			final double y = random.nextGaussian() * d[j];
			// Get the squared jump distance
			data[i] = x * x + y * y;
		}
		for (int i = 0; i < f.length; i++)
			f[i] /= data.length;
		return data;
	}

	private int pick(double[] f, double nextDouble)
	{
		for (int i = 0; i < f.length; i++)
			if (nextDouble < f[i])
				return i;
		return f.length - 1;
	}

	void log(String format, Object... args)
	{
		System.out.printf(format, args);
	}
}
