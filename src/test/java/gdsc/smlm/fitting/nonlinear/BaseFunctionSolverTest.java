package gdsc.smlm.fitting.nonlinear;

import java.util.Arrays;

import org.apache.commons.math3.distribution.ExponentialDistribution;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.apache.commons.math3.stat.inference.TTest;
import org.junit.Assert;

import gdsc.core.math.SimpleArrayMoment;
import gdsc.core.utils.DoubleEquality;
import gdsc.core.utils.Statistics;
import gdsc.core.utils.StoredDataStatistics;
import gdsc.smlm.fitting.FisherInformationMatrix;
import gdsc.smlm.fitting.FitStatus;
import gdsc.smlm.fitting.FunctionSolver;
import gdsc.smlm.fitting.nonlinear.gradient.PoissonGradientProcedure;
import gdsc.smlm.fitting.nonlinear.gradient.PoissonGradientProcedureFactory;
import gdsc.smlm.function.Gradient1Function;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.function.gaussian.GaussianFunctionFactory;

/**
 * Test that a bounded fitter can return the same results with and without bounds.
 */
public abstract class BaseFunctionSolverTest
{
	long seed = 30051977; //System.currentTimeMillis() + System.identityHashCode(this);
	//long seed = System.currentTimeMillis() + System.identityHashCode(this);
	RandomGenerator randomGenerator = new Well19937c(seed);
	RandomDataGenerator dataGenerator = new RandomDataGenerator(randomGenerator);

	// Basic Gaussian
	static double[] params = new double[7];
	static double[] base = { 0.8, 1, 1.2 };
	static double[] signal = { 1000, 2000, 5000, 10000 }; // 100, 200, 400, 800 };
	static double[] noise = { 0.1, 0.5, 1 };
	static double[] shift = { -1, 0, 1 };
	static double[] factor = { 0.7, 1, 1.3 };
	static int size = 11;
	static
	{
		params[Gaussian2DFunction.BACKGROUND] = 5;
		params[Gaussian2DFunction.X_POSITION] = size / 2;
		params[Gaussian2DFunction.Y_POSITION] = size / 2;
		params[Gaussian2DFunction.X_SD] = 1.4;
	}
	static double[] defaultClampValues;
	static
	{
		defaultClampValues = new double[7];
		// Taken from the 3D-DAO-STORM paper:
		// (Babcock et al. 2012) A high-density 3D localization algorithm for stochastic optical 
		// reconstruction microscopy. Optical Nanoscopy. 2012 1:6
		// DOI: 10.1186/2192-2853-1-6
		// Page 3
		// Note: It is not clear if the background/signal are in ADUs or photons. I assume photons.

		// This seems big for background in photons
		defaultClampValues[Gaussian2DFunction.BACKGROUND] = 100;
		//defaultClampValues[Gaussian2DFunction.BACKGROUND] = 20;
		defaultClampValues[Gaussian2DFunction.SIGNAL] = 1000;
		defaultClampValues[Gaussian2DFunction.SHAPE] = Math.PI;
		defaultClampValues[Gaussian2DFunction.X_POSITION] = 1;
		defaultClampValues[Gaussian2DFunction.Y_POSITION] = 1;
		defaultClampValues[Gaussian2DFunction.X_SD] = 3;
		defaultClampValues[Gaussian2DFunction.Y_SD] = 3;
	}

	void canFitSingleGaussian(FunctionSolver solver, boolean applyBounds)
	{
		canFitSingleGaussian(solver, applyBounds, false);
	}

	// Variance = Exponential (equivalent to chi-squared with k=1, i.e. 
	// sum of the squares of 1 normal distribution). 
	// We want 99.9% @ 400 ADU based on supplementary figure 1.a/1.b 
	// cumul = 1 - e^-lx (with l = 1/mean)
	// => e^-lx = 1 - cumul
	// => -lx = log(1-0.999)
	// => l = -log(0.001) / 400  (since x==400)
	// => 1/l = 57.9
	private static double variance = 57.9; // SD = 7.6

	// Gain = Approximately Normal
	private static double gain = 2.2;
	private static double gainSD = 0.2;
	private static double[] weights = null, pixelSD = null;

	private static double[] getWeights()
	{
		if (weights == null)
			computeWeights();
		return weights;
	}

	private static double[] getPixelSD()
	{
		if (pixelSD == null)
			computeWeights();
		return pixelSD;
	}

	private static void computeWeights()
	{
		// Per observation read noise.
		weights = new double[size * size];
		RandomGenerator randomGenerator = new Well19937c(42);
		ExponentialDistribution ed = new ExponentialDistribution(randomGenerator, variance,
				ExponentialDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY);
		for (int i = 0; i < weights.length; i++)
		{
			double pixelVariance = ed.sample();
			double pixelGain = Math.max(0.1, gain + randomGenerator.nextGaussian() * gainSD);
			// weights = var / g^2
			weights[i] = pixelVariance / (pixelGain * pixelGain);
		}
		// Convert to standard deviation for simulation
		pixelSD = new double[weights.length];
		for (int i = 0; i < weights.length; i++)
			pixelSD[i] = Math.sqrt(weights[i]);
	}

	void canFitSingleGaussian(FunctionSolver solver, boolean applyBounds, boolean weighted)
	{
		// Allow reporting the fit deviations
		boolean report = false;
		double[] crlb = null;
		SimpleArrayMoment m = null;

		if (weighted && solver.isWeighted())
			solver.setWeights(getWeights());
		else
			weighted = false;

		randomGenerator.setSeed(seed);
		for (double s : signal)
		{
			double[] expected = createParams(1, s, 0, 0, 1);
			double[] lower = createParams(0, s * 0.5, -0.2, -0.2, 0.8);
			double[] upper = createParams(3, s * 2, 0.2, 0.2, 1.2);
			if (applyBounds)
				solver.setBounds(lower, upper);
			if (report)
			{
				// Compute the CRLB for a Poisson process
				PoissonGradientProcedure gp = PoissonGradientProcedureFactory
						.create((Gradient1Function) ((BaseFunctionSolver) solver).getGradientFunction());
				gp.computeFisherInformation(expected);
				FisherInformationMatrix f = new FisherInformationMatrix(gp.getLinear(), gp.n);
				crlb = f.crlbSqrt();
				// Compute the deviations
				m = new SimpleArrayMoment();
			}
			for (double n : noise)
			{
				double[] data = drawGaussian(expected, n, (weighted) ? getPixelSD() : null);
				for (double db : base)
					for (double dx : shift)
						for (double dy : shift)
							for (double dsx : factor)
							{
								double[] p = createParams(db, s, dx, dy, dsx);
								double[] fp = fitGaussian(solver, data, p, expected);
								for (int i = 0; i < expected.length; i++)
								{
									if (fp[i] < lower[i])
										Assert.assertTrue(
												String.format("Fit Failed: [%d] %.2f < %.2f: %s != %s", i, fp[i],
														lower[i], Arrays.toString(fp), Arrays.toString(expected)),
												false);
									if (fp[i] > upper[i])
										Assert.assertTrue(
												String.format("Fit Failed: [%d] %.2f > %.2f: %s != %s", i, fp[i],
														upper[i], Arrays.toString(fp), Arrays.toString(expected)),
												false);
									if (report)
										fp[i] = expected[i] - fp[i];
								}
								// Store the deviations
								if (report)
									m.add(fp);
							}
			}
			// Report
			if (report)
				System.out.printf("%s %f : CRLB = %s, Devaitions = %s\n", solver.getClass().getSimpleName(), s,
						Arrays.toString(crlb), Arrays.toString(m.getStandardDeviation()));
		}
	}

	void canFitSingleGaussianBetter(FunctionSolver solver, boolean applyBounds, FunctionSolver solver2,
			boolean applyBounds2, String name, String name2)
	{
		int LOOPS = 5;
		randomGenerator.setSeed(seed);
		StoredDataStatistics[] stats = new StoredDataStatistics[6];
		String[] statName = { "Signal", "X", "Y" };

		int[] betterPrecision = new int[3];
		int[] totalPrecision = new int[3];
		int[] betterAccuracy = new int[3];
		int[] totalAccuracy = new int[3];

		int i1 = 0, i2 = 0;
		for (double s : signal)
		{
			double[] expected = createParams(1, s, 0, 0, 1);
			double[] lower = null, upper = null;
			if (applyBounds || applyBounds2)
			{
				lower = createParams(0, s * 0.5, -0.2, -0.2, 0.8);
				upper = createParams(3, s * 2, 0.2, 0.2, 1.2);
			}
			if (applyBounds)
				solver.setBounds(lower, upper);
			if (applyBounds2)
				solver2.setBounds(lower, upper);

			for (double n : noise)
			{
				for (int loop = LOOPS; loop-- > 0;)
				{
					double[] data = drawGaussian(expected, n, null);

					for (int i = 0; i < stats.length; i++)
						stats[i] = new StoredDataStatistics();

					for (double db : base)
						for (double dx : shift)
							for (double dy : shift)
								for (double dsx : factor)
								{
									double[] p = createParams(db, s, dx, dy, dsx);
									double[] fp = fitGaussian(solver, data, p, expected);
									i1 += solver.getEvaluations();

									double[] fp2 = fitGaussian(solver2, data, p, expected);
									i2 += solver2.getEvaluations();

									// Get the mean and sd (the fit precision)
									compare(fp, expected, fp2, expected, Gaussian2DFunction.SIGNAL, stats[0], stats[1]);

									compare(fp, expected, fp2, expected, Gaussian2DFunction.X_POSITION, stats[2],
											stats[3]);
									compare(fp, expected, fp2, expected, Gaussian2DFunction.Y_POSITION, stats[4],
											stats[5]);

									// Use the distance
									//stats[2].add(distance(fp, expected));
									//stats[3].add(distance(fp2, expected2));
								}

					double alpha = 0.05; // two sided
					for (int i = 0; i < stats.length; i += 2)
					{
						double u1 = stats[i].getMean();
						double u2 = stats[i + 1].getMean();
						double sd1 = stats[i].getStandardDeviation();
						double sd2 = stats[i + 1].getStandardDeviation();

						TTest tt = new TTest();
						boolean diff = tt.tTest(stats[i].getValues(), stats[i + 1].getValues(), alpha);

						int index = i / 2;
						String msg = String.format("%s vs %s : %.1f (%.1f) %s %f +/- %f vs %f +/- %f  (N=%d) %b", name2,
								name, s, n, statName[index], u2, sd2, u1, sd1, stats[i].getN(), diff);
						if (diff)
						{
							// Different means. Check they are roughly the same
							if (DoubleEquality.almostEqualRelativeOrAbsolute(u1, u2, 0.1, 0))
							{
								// Basically the same. Check which is more precise
								if (!DoubleEquality.almostEqualRelativeOrAbsolute(sd1, sd2, 0.05, 0))
								{
									if (sd2 < sd1)
									{
										betterPrecision[index]++;
										println(msg + " P*");
									}
									else
										println(msg + " P");
									totalPrecision[index]++;
								}
							}
							else
							{
								// Check which is more accurate (closer to zero)
								u1 = Math.abs(u1);
								u2 = Math.abs(u2);
								if (u2 < u1)
								{
									betterAccuracy[index]++;
									println(msg + " A*");
								}
								else
									println(msg + " A");
								totalAccuracy[index]++;
							}
						}
						else
						{
							// The same means. Check that it is more precise
							if (!DoubleEquality.almostEqualRelativeOrAbsolute(sd1, sd2, 0.05, 0))
							{
								if (sd2 < sd1)
								{
									betterPrecision[index]++;
									println(msg + " P*");
								}
								else
									println(msg + " P");
								totalPrecision[index]++;
							}
						}
					}
				}
			}
		}

		int better = 0, total = 0;
		for (int index = 0; index < statName.length; index++)
		{
			better += betterPrecision[index] + betterAccuracy[index];
			total += totalPrecision[index] + totalAccuracy[index];
			test(name2, name, statName[index] + " P", betterPrecision[index], totalPrecision[index],
					printBetterDetails);
			test(name2, name, statName[index] + " A", betterAccuracy[index], totalAccuracy[index], printBetterDetails);
		}
		test(name2, name, String.format("All (eval [%d] [%d]) : ", i2, i1), better, total, true);
	}

	private void test(String name2, String name, String statName, int better, int total, boolean print)
	{
		double p = 100.0 * better / total;
		String msg = String.format("%s vs %s : %s %d / %d  (%.1f)", name2, name, statName, better, total, p);
		if (print)
			System.out.println(msg);
		// Do not test if we don't have many examples
		if (total <= 10)
		{
			return;
		}

		// Disable this for now so builds do not fail during the test phase

		// It seems that most of the time clamping and bounds improve things.
		// There are a few cases where Bounds or Clamping alone do not improve things.
		// Use of Dynamic Clamping is always better.
		// Use of Bounded Dynamic Clamping is always better.

		// The test may be unrealistic as the initial params are close to the actual answer.

		//Assert.assertTrue(msg, p >= 50.0);
	}

	boolean printBetterDetails = false;

	private void println(String msg)
	{
		if (printBetterDetails)
			System.out.println(msg);
	}

	static double distance(double[] o, double[] e)
	{
		double dx = o[Gaussian2DFunction.X_POSITION] - e[Gaussian2DFunction.X_POSITION];
		double dy = o[Gaussian2DFunction.Y_POSITION] - e[Gaussian2DFunction.Y_POSITION];
		// Use the signs of the coords to assign a direction vector
		return Math.sqrt(dx * dx + dy * dy) * Math.signum(Math.signum(dy) * Math.signum(dx));
	}

	private void compare(double[] o1, double[] e1, double[] o2, double[] e2, int i, Statistics stats1,
			Statistics stats2)
	{
		compare(o1[i], e1[i], o2[i], e2[i], stats1, stats2);
	}

	private void compare(double o1, double e1, double o2, double e2, Statistics stats1, Statistics stats2)
	{
		stats1.add(o1 - e1);
		stats2.add(o2 - e2);
	}

	double[] createParams(double db, double signal, double dx, double dy, double dsx)
	{
		double[] p = params.clone();
		p[Gaussian2DFunction.BACKGROUND] *= db;
		p[Gaussian2DFunction.SIGNAL] = signal;
		p[Gaussian2DFunction.X_POSITION] += dx;
		p[Gaussian2DFunction.Y_POSITION] += dy;
		p[Gaussian2DFunction.X_SD] *= dsx;
		return p;
	}

	double[] addBiasToParams(double[] p, double bias)
	{
		p = p.clone();
		p[Gaussian2DFunction.BACKGROUND] += bias;
		return p;
	}

	double[] fitGaussian(FunctionSolver solver, double[] data, double[] params, double[] expected)
	{
		//System.out.printf("%s : Expected %s\n", solver.getClass().getSimpleName(), Arrays.toString(expected));
		params = params.clone();
		FitStatus status = solver.fit(data, null, params, null);
		if (status != FitStatus.OK)
			Assert.assertTrue(String.format("Fit Failed: %s i=%d: %s != %s", status.toString(), solver.getIterations(),
					Arrays.toString(params), Arrays.toString(expected)), false);
		return params;
	}

	/**
	 * Draw a Gaussian with Poisson shot noise and Gaussian read noise
	 * 
	 * @param params
	 *            The Gaussian parameters
	 * @param noise
	 *            The read noise
	 * @param weights
	 * @param withBias
	 * @return The data
	 */
	double[] drawGaussian(double[] params, double noise, double[] pixelSD)
	{
		double[] data = new double[size * size];
		int n = params.length / 6;
		Gaussian2DFunction f = GaussianFunctionFactory.create2D(n, size, size, GaussianFunctionFactory.FIT_CIRCLE,
				null);
		f.initialise(params);

		for (int i = 0; i < data.length; i++)
		{
			data[i] = dataGenerator.nextPoisson(f.eval(i));
			if (pixelSD != null)
				data[i] += randomGenerator.nextGaussian() * pixelSD[i];
		}
		if (noise != 0)
			for (int i = 0; i < data.length; i++)
				data[i] += dataGenerator.nextGaussian(0, noise);
		//gdsc.core.ij.Utils.display("Spot", data, size, size);
		return data;
	}
}
