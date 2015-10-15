package gdsc.smlm.fitting;

import gdsc.smlm.fitting.JumpDistanceAnalysis.JumpDistanceCumulFunction;
import gdsc.smlm.fitting.JumpDistanceAnalysis.JumpDistanceFunction;
import gdsc.smlm.fitting.JumpDistanceAnalysis.MixedJumpDistanceCumulFunction;
import gdsc.smlm.fitting.JumpDistanceAnalysis.MixedJumpDistanceFunction;
import gdsc.smlm.utils.logging.Logger;

import java.util.ArrayList;
import java.util.Arrays;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.junit.Assert;
import org.junit.Test;

public class JumpDistanceAnalysisTest
{

	//Based on the paper: Weimann, L., Ganzinger, K.A., McColl, J., Irvine, K.L., Davis, S.J., 
	//Gay, N.J., Bryant, C.E., Klenerman, D. (2013) A Quantitative Comparison of Single-Dye 
	//Tracking Analysis Tools Using Monte Carlo Simulations. PLoS One 8, Issue 5, e64287
	// 
	//This paper simulated tracks using 150 particle over 30 frames with a SNR of 6. The spots
	//were fit and then MSD or Jump Distance analysis performed. This means that each image 
	//could have 150*29 jumps = 4350. They compute the fit using 750 trajectories (21750 jumps)
	//and repeat this 10 times to get a mean fit. They indicate success if 
	//D is within 10% of the true value (the threshold for F is not explicitly stated
	//but appears to be around 20% or less). 2 population sample used D = 0.1, 0.02. 
	//Good results were obtained when F(mobile) > 20%.

	// TODO - revise this so that fitting always works and then reduce the sample size down gradually
	// to see if there is a fail limit.

	// Test MLE fitting and histogram fitting separately.
	// Is MLE fitting worth doing. Can it be made better?

	double deltaD = 0.1;
	double deltaF = 0.2;
	// Used for testing single populations
	// Used for testing dual populations: 15-fold, 5-fold, 3-fold difference between pairs
	double[] D = new double[] { 0.2, 3, 1 };
	RandomGenerator random = new Well19937c(System.currentTimeMillis() + System.identityHashCode(this));

	@Test
	public void canIntegrateProbabilityToCumulativeWithSinglePopulation()
	{
		JumpDistanceAnalysis jd = new JumpDistanceAnalysis();
		jd.setMinD(0);
		jd.setMinFraction(0);
		SimpsonIntegrator si = new SimpsonIntegrator(1e-3, 1e-8, 2, SimpsonIntegrator.SIMPSON_MAX_ITERATIONS_COUNT);
		for (double d : D)
		{
			final double[] params = new double[] { d };
			final JumpDistanceFunction fp = jd.new JumpDistanceFunction(null, d);
			JumpDistanceCumulFunction fc = jd.new JumpDistanceCumulFunction(null, null, d);
			double x = d / 8;
			UnivariateFunction func = new UnivariateFunction()
			{
				public double value(double x)
				{
					return fp.evaluate(x, params);
				}
			};
			for (int i = 1; i < 10; i++, x *= 2)
			{
				double e = fc.evaluate(x, params);
				// Integrate
				double o = si.integrate(10000, func, 0, x);
				//log("Integrate d=%.1f : x=%.1f, e=%f, o=%f, iter=%d, eval=%d\n", d, x, e, o, si.getIterations(),
				//		si.getEvaluations());
				Assert.assertEquals("Failed to integrate", e, o, e * 1e-2);
			}
		}
	}

	@Test
	public void canIntegrateProbabilityToCumulativeWithMixedPopulation()
	{
		JumpDistanceAnalysis jd = new JumpDistanceAnalysis();
		jd.setMinD(0);
		jd.setMinFraction(0);
		SimpsonIntegrator si = new SimpsonIntegrator(1e-3, 1e-8, 2, SimpsonIntegrator.SIMPSON_MAX_ITERATIONS_COUNT);
		for (double d : D)
		{
			for (double f : new double[] { 0, 0.1, 0.2, 0.4, 0.7, 0.9, 1 })
			{
				final double[] params = new double[] { f, d, 1 - f, d * 0.1 };
				final MixedJumpDistanceFunction fp = jd.new MixedJumpDistanceFunction(null, d, 2);
				MixedJumpDistanceCumulFunction fc = jd.new MixedJumpDistanceCumulFunction(null, null, d, 2);
				double x = d / 8;
				UnivariateFunction func = new UnivariateFunction()
				{
					public double value(double x)
					{
						return fp.evaluate(x, params);
					}
				};
				for (int i = 1; i < 10; i++, x *= 2)
				{
					double e = fc.evaluate(x, params);
					// Integrate
					double o = si.integrate(10000, func, 0, x);
					//log("Integrate d=%.1f, f=%.1f : x=%.1f, e=%f, o=%f, iter=%d, eval=%d\n", d, f, x, e, o,
					//		si.getIterations(), si.getEvaluations());
					Assert.assertEquals("Failed to integrate", e, o, e * 1e-2);
				}
			}
		}
	}

	@Test
	public void canFitSinglePopulationMLE()
	{
		fitSinglePopulation(true);
	}

	@Test
	public void canFitSinglePopulation()
	{
		fitSinglePopulation(false);
	}

	private void fitSinglePopulation(boolean mle)
	{
		String title = String.format("%s Single  ", (mle) ? "MLE" : "LSQ");
		AssertionError error = null;
		NEXT_D: for (double d : D)
		{
			for (int samples = 500, k=0; k<6; samples *= 2, k++)
			{
				try
				{
					fit(title, samples, 0, new double[] { d }, new double[] { 1 }, mle);
					error = null;
					continue NEXT_D;
				}
				catch (AssertionError e)
				{
					error = e;
				}
			}
			if (error != null)
				throw error;
		}
	}

	@Test
	public void canFitDual0_2PopulationMLE()
	{
		fitDualPopulation(true, 0.2);
	}

	@Test
	public void canFitDual0_2Population()
	{
		fitDualPopulation(false, 0.2);
	}

	@Test
	public void canFitDual0_5PopulationMLE()
	{
		fitDualPopulation(true, 0.5);
	}

	@Test
	public void canFitDual0_5Population()
	{
		fitDualPopulation(false, 0.5);
	}

	@Test
	public void canFitDual0_8PopulationMLE()
	{
		fitDualPopulation(true, 0.8);
	}

	@Test
	public void canFitDual0_8Population()
	{
		fitDualPopulation(false, 0.8);
	}

	private void fitDualPopulation(boolean mle, double fraction)
	{
		String title = String.format("%s Dual=%.1f", (mle) ? "MLE" : "LSQ", fraction);
		AssertionError error = null;
		for (int i = 0; i < D.length; i++)
		{
			NEXT_D: for (int j = i + 1; j < D.length; j++)
			{
				for (int samples = 500, k=0; k<6; samples *= 2, k++)
				{
					try
					{
						fit(title, samples, 0, new double[] { D[i], D[j] }, new double[] { fraction, 1 - fraction },
								mle);
						error = null;
						continue NEXT_D;
					}
					catch (AssertionError e)
					{
						error = e;
					}
				}
				if (error != null)
					throw error;
			}
		}
	}

	private void fit(String title, int samples, int n, double[] d, double[] f, boolean mle)
	{
		JumpDistanceAnalysis.sort(d, f);
		double[] jumpsDistances = createData(samples, d, f);
		Logger logger = null;
		logger = new gdsc.smlm.utils.logging.ConsoleLogger();
		JumpDistanceAnalysis jd = new JumpDistanceAnalysis(logger);
		jd.setFitRestarts(5);
		jd.setMinFraction(0);
		jd.setMinDifference(2);
		jd.setN((n > 0) ? n : 10);
		double[][] fit = (mle) ? jd.fitJumpDistancesMLE(jumpsDistances) : jd.fitJumpDistances(jumpsDistances);
		double[] fitD = fit[0];
		double[] fitF = fit[1];
		AssertionError error = null;
		try
		{
			Assert.assertEquals("Failed to fit n", d.length, fitD.length);
			for (int i = 0; i < d.length; i++)
			{
				Assert.assertEquals("Failed to fit d", d[i], fitD[i], deltaD * d[i]);
				Assert.assertEquals("Failed to fit f", f[i], fitF[i], deltaF * f[i]);
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
			log("%s %s N=%d sample=%d, n<=%d : %s = %s [%s] : %s = %s [%s]\n", (error == null) ? "+++ Pass"
					: "--- Fail", title, d.length, samples, n, toString(d), toString(fitD), toString(e1), toString(f),
					toString(fitF), toString(e2));
			if (error != null)
				throw error;
		}
	}

	private double[] getError(double[] e, double[] o)
	{
		double[] error = new double[Math.min(e.length, o.length)];
		for (int i = 0; i < error.length; i++)
			//error[i] = DoubleEquality.relativeError(o[i], e[i]);
			error[i] = 100.0 * (o[i] - e[i]) / e[i];
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

	class DataSample
	{
		double[] d, f, s;
		double[] data = null;
		int[] sample = null;

		DataSample(double[] d, double[] f)
		{
			this.d = d.clone();

			// Convert diffusion co-efficient into the standard deviation for the random move in each dimension
			// For 1D diffusion: sigma^2 = 2D
			//                   sigma = sqrt(2D)
			// See: https://en.wikipedia.org/wiki/Brownian_motion#Einstein.27s_theory
			s = new double[d.length];
			double sum = 0;
			for (int i = 0; i < s.length; i++)
			{
				s[i] = Math.sqrt(2 * d[i]);
				sum += f[i];
			}
			this.f = new double[f.length];
			for (int i = 0; i < f.length; i++)
				this.f[i] = f[i] / sum;
		}

		public boolean equals(Object obj)
		{
			if (!(obj instanceof DataSample))
				return super.equals(obj);
			DataSample that = (DataSample) obj;
			if (that.d.length != this.d.length)
				return false;
			for (int i = d.length; i-- > 0;)
			{
				if (that.d[i] != this.d[i])
					return false;
				if (that.f[i] != this.f[i])
					return false;
			}
			return true;
		}

		void add(double[] data2, int[] sample2)
		{
			if (data == null)
			{
				data = data2;
				sample = sample2;
				return;
			}
			int size = data.length;
			int newSize = size + data2.length;
			data = Arrays.copyOf(data, newSize);
			sample = Arrays.copyOf(sample, newSize);
			System.arraycopy(data2, 0, data, size, data2.length);
			System.arraycopy(sample2, 0, sample, size, sample2.length);
		}

		int getSize()
		{
			return (data == null) ? 0 : data.length;
		}

		double[][] getSample(int size)
		{
			if (size > getSize())
			{
				int extra = size - getSize();

				// Get cumulative fraction
				double[] c = new double[f.length];
				double sum = 0;
				for (int i = 0; i < f.length; i++)
				{
					sum += f[i];
					c[i] = sum;
				}

				double[] data = new double[extra];

				// Pick the population using the fraction.
				// Do this before sampling since the nextGaussian function computes random variables
				// in pairs so we want to process all the same sample together
				int[] sample = new int[extra];
				if (c.length > 1)
					for (int i = 0; i < data.length; i++)
						sample[i] = pick(c, random.nextDouble());
				Arrays.sort(sample);

				for (int i = 0; i < data.length; i++)
				{
					// Pick the population using the fraction
					final int j = sample[i];
					// Get the x/y shifts
					final double x = random.nextGaussian() * s[j];
					final double y = random.nextGaussian() * s[j];
					// Get the squared jump distance
					data[i] = x * x + y * y;
				}
				add(data, sample);
			}

			// Build the sample data and return the D and fractions
			double[] data = Arrays.copyOf(this.data, size);
			double[] d = new double[this.d.length];
			double[] f = new double[d.length];
			for (int i = 0; i < size; i++)
			{
				final int j = sample[i];
				d[j] += data[i];
				f[j]++;
			}
			for (int i = 0; i < d.length; i++)
			{
				// 4D = MSD
				// D = MSD / 4
				d[i] = (d[i] / f[i]) / 4;
				f[i] /= size;
			}

			return new double[][] { data, d, f };
		}
	}

	static ArrayList<DataSample> samples = new ArrayList<DataSample>();

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
		// Cache the data so that if we run a second test with 
		// the same d and f we use the same data
		DataSample sample = new DataSample(d, f);
		int index = samples.indexOf(sample);
		if (index != -1)
			sample = samples.get(index);
		else
			samples.add(sample);

		double[][] dataSample = sample.getSample(n);
		double[] data = dataSample[0];
		double[] d2 = dataSample[1];
		double[] f2 = dataSample[2];

		// Update with the real values
		for (int i = 0; i < d.length; i++)
		{
			d[i] = d2[i];
			f[i] = f2[i];
		}

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
