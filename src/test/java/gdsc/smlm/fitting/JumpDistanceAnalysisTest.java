package gdsc.smlm.fitting;

import gdsc.smlm.utils.logging.Logger;

import java.util.Arrays;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.junit.Assert;
import org.junit.Test;

public class JumpDistanceAnalysisTest
{
	/*
	 * Based on the paper: Weimann, L., Ganzinger, K.A., McColl, J., Irvine, K.L., Davis, S.J., Gay, N.J., Bryant, C.E.,
	 * Klenerman, D. (2013) A Quantitative Comparison of Single-Dye Tracking Analysis Tools Using Monte Carlo
	 * Simulations. PLoS One 8, Issue 5, e64287
	 * 
	 * This paper simulated tracks using 150 particle over 30 frames with a SNR of 6. The spots were fit and then MSD or
	 * Jump Distance analysis performed. This means that each image could have 150*29 jumps = 4350. They compute
	 * the fit using 750 trajectories (21750 jumps) and repeat this 10 times to get a mean fit. They indicate success if
	 * D is within 10% of the true value (the threshold for F is not explicitly stated but appears to be around 20% or
	 * less). 2 population sample used D = 0.1, 0.02. Good results were obtained when F(mobile) > 20%.
	 */

	// TODO - revise this so that fitting always works and then reduce the sample size down gradually
	// to see if there is a fail limit.
	
	double delta = 2e-1;
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
		Logger logger = null; //new ConsoleLogger();
		JumpDistanceAnalysis jd = new JumpDistanceAnalysis(logger);
		jd.setFitRestarts(5);
		jd.setMinFraction(0.1);
		jd.setMinDifference(4);
		jd.setN((n > 0) ? n : 10);
		double[][] fit = jd.fitJumpDistancesMLE(jumpsDistances);
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
			double[] e1 = getError(d, fitD);
			double[] e2 = getError(f, fitF);
			log("%s sample=%d, n<=%d : %s = %s @ %s : %s = %s @ %s\n", (error == null) ? "+++ Pass" : "--- Fail",
					samples, n, toString(d), toString(fitD), toString(e1), toString(f), toString(fitF), toString(e2));
			if (error != null)
				throw error;
		}
	}

	private double[] getError(double[] e, double[] o)
	{
		double[] error = new double[o.length];
		for (int i = 0; i < o.length; i++)
			//error[i] = DoubleEquality.relativeError(o[i], e[i]);
			error[i] = (o[i] - e[i]) / e[i];
		return error;
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
		// See: https://en.wikipedia.org/wiki/Brownian_motion#Einstein.27s_theory
		double[] s = new double[d.length];
		double sum = 0;
		for (int i = 0; i < s.length; i++)
		{
			// 1D simulation
			s[i] = Math.sqrt(2 * d[i]);
			sum += f[i];
		}
		// Normalise fractions to 1
		for (int i = 0; i < f.length; i++)
			f[i] /= sum;

		// Get cumulative fraction
		double[] c = new double[f.length];
		sum = 0;
		for (int i = 0; i < f.length; i++)
		{
			sum += f[i];
			c[i] = sum;
			// Reset to count the actual fractions
			f[i] = 0;
		}

		double[] data = new double[n];

		// Pick the population using the fraction.
		// Do this before sampling since the nextGaussian function computes random variables
		// in pairs so we want to process all the same sample together
		int[] next = new int[n];
		if (c.length > 1)
			for (int i = 0; i < data.length; i++)
				next[i] = pick(c, random.nextDouble());
		Arrays.sort(next);

		for (int i = 0; i < data.length; i++)
		{
			// Pick the population using the fraction
			final int j = next[i];
			f[j]++;

			// Get the x/y shifts
			final double x = random.nextGaussian() * s[j];
			final double y = random.nextGaussian() * s[j];
			// Get the squared jump distance
			data[i] = x * x + y * y;
		}
		for (int i = 0; i < f.length; i++)
			f[i] /= data.length;

		// Debug 
		//gdsc.smlm.utils.StoredDataStatistics stats = new gdsc.smlm.utils.StoredDataStatistics(data);
		//gdsc.smlm.ij.utils.Utils.showHistogram(
		//		"MSD",
		//		stats,
		//		"MSD",
		//		0,
		//		0,
		//		n / 10,
		//		true,
		//		String.format("%s : %s : u=%f, sd=%f", toString(d), toString(f), stats.getMean(),
		//				stats.getStandardDeviation()));

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
