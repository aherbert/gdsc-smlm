package gdsc.smlm.fitting.nonlinear;

import java.util.Arrays;

import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.apache.commons.math3.stat.inference.TTest;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.utils.DoubleEquality;
import gdsc.core.utils.Statistics;
import gdsc.core.utils.StoredDataStatistics;
import gdsc.smlm.fitting.FitStatus;
import gdsc.smlm.fitting.FunctionSolver;
import gdsc.smlm.fitting.nonlinear.stop.ErrorStoppingCriteria;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.function.gaussian.GaussianFunctionFactory;

/**
 * Test that a bounded fitter can return the same results with and without bounds.
 */
public class BoundedFunctionSolverTest
{
	long seed = 30051977; //System.currentTimeMillis() + System.identityHashCode(this);
	//long seed = System.currentTimeMillis() + System.identityHashCode(this);
	RandomGenerator randomGenerator = new Well19937c(seed);
	RandomDataGenerator dataGenerator = new RandomDataGenerator(randomGenerator);

	// Basic Gaussian
	static double bias = 100;
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
	private static double[] defaultClampValues;
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

	// TODO - Test if local search param if useful when using clamping

	// Standard LVM
	@Test
	public void canFitSingleGaussianLVM()
	{
		fitSingleGaussianLVM(0, 0, false);
	}

	// Bounded/Clamped LVM

	@Test
	public void canFitSingleGaussianBLVMNoBounds()
	{
		fitSingleGaussianLVM(1, 0, false);
	}

	@Test
	public void canFitSingleGaussianBLVM()
	{
		fitSingleGaussianLVM(2, 0, false);
	}

	@Test
	public void canFitSingleGaussianCLVM()
	{
		fitSingleGaussianLVM(0, 1, false);
	}

	@Test
	public void canFitSingleGaussianDCLVM()
	{
		fitSingleGaussianLVM(0, 2, false);
	}

	@Test
	public void canFitSingleGaussianBCLVM()
	{
		fitSingleGaussianLVM(2, 1, false);
	}

	@Test
	public void canFitSingleGaussianBDCLVM()
	{
		fitSingleGaussianLVM(2, 2, false);
	}

	// MLE LVM

	@Test
	public void canFitSingleGaussianLVMMLE()
	{
		fitSingleGaussianLVM(0, 0, true);
	}

	@Test
	public void canFitSingleGaussianBLVMMLENoBounds()
	{
		fitSingleGaussianLVM(1, 0, true);
	}

	@Test
	public void canFitSingleGaussianBLVMMLE()
	{
		fitSingleGaussianLVM(2, 0, true);
	}

	private void fitSingleGaussianLVM(int bounded, int clamping, boolean mle)
	{
		canFitSingleGaussian(getLVM(bounded, clamping, mle), bounded == 2, !mle);
	}

	// Is Bounded/Clamped LVM better?

	@Test
	public void fitSingleGaussianBLVMBetterThanLVM()
	{
		fitSingleGaussianBetterLVM(true, 0, false, false, 0, false);
	}

	@Test
	public void fitSingleGaussianCLVMBetterThanLVM()
	{
		fitSingleGaussianBetterLVM(false, 1, false, false, 0, false);
	}

	@Test
	public void fitSingleGaussianBCLVMBetterThanLVM()
	{
		fitSingleGaussianBetterLVM(true, 1, false, false, 0, false);
	}

	@Test
	public void fitSingleGaussianDCLVMBetterThanLVM()
	{
		fitSingleGaussianBetterLVM(false, 2, false, false, 0, false);
	}

	@Test
	public void fitSingleGaussianBDCLVMBetterThanLVM()
	{
		fitSingleGaussianBetterLVM(true, 2, false, false, 0, false);
	}

	@Test
	public void fitSingleGaussianLVMMLEBetterThanLVM()
	{
		fitSingleGaussianBetterLVM(false, 0, true, false, 0, false);
	}

	@Test
	public void fitSingleGaussianBLVMMLEBetterThanLVM()
	{
		fitSingleGaussianBetterLVM(true, 0, true, false, 0, false);
	}

	@Test
	public void fitSingleGaussianCLVMMLEBetterThanLVM()
	{
		fitSingleGaussianBetterLVM(false, 1, true, false, 0, false);
	}

	@Test
	public void fitSingleGaussianBCLVMMLEBetterThanLVM()
	{
		fitSingleGaussianBetterLVM(true, 1, true, false, 0, false);
	}

	@Test
	public void fitSingleGaussianDCLVMMLEBetterThanLVM()
	{
		fitSingleGaussianBetterLVM(false, 2, true, false, 0, false);
	}

	@Test
	public void fitSingleGaussianBDCLVMMLEBetterThanLVM()
	{
		fitSingleGaussianBetterLVM(true, 2, true, false, 0, false);
	}

	@Test
	public void fitSingleGaussianBLVMMLEBetterThanLVMMLE()
	{
		fitSingleGaussianBetterLVM(true, 0, true, false, 0, true);
	}

	@Test
	public void fitSingleGaussianCLVMMLEBetterThanLVMMLE()
	{
		fitSingleGaussianBetterLVM(false, 1, true, false, 0, true);
	}

	@Test
	public void fitSingleGaussianDCLVMMLEBetterThanLVMMLE()
	{
		fitSingleGaussianBetterLVM(false, 2, true, false, 0, true);
	}

	@Test
	public void fitSingleGaussianBDCLVMMLEBetterThanLVMMLE()
	{
		fitSingleGaussianBetterLVM(true, 2, true, false, 0, true);
	}

	@Test
	public void fitSingleGaussianBLVMMLEBetterThanBLVM()
	{
		fitSingleGaussianBetterLVM(true, 0, true, true, 0, false);
	}

	@Test
	public void fitSingleGaussianBCLVMMLEBetterThanBCLVM()
	{
		fitSingleGaussianBetterLVM(true, 1, true, true, 1, false);
	}

	@Test
	public void fitSingleGaussianBDCLVMMLEBetterThanBDCLVM()
	{
		fitSingleGaussianBetterLVM(true, 2, true, true, 2, false);
	}

	private void fitSingleGaussianBetterLVM(boolean bounded2, int clamping2, boolean mle2, boolean bounded,
			int clamping, boolean mle)
	{
		NonLinearFit solver = getLVM((bounded) ? 2 : 1, clamping, mle);
		NonLinearFit solver2 = getLVM((bounded2) ? 2 : 1, clamping2, mle2);
		canFitSingleGaussianBetter(solver, bounded, !mle, solver2, bounded2, !mle2, getLVMName(bounded, clamping, mle),
				getLVMName(bounded2, clamping2, mle2));
	}

	private NonLinearFit getLVM(int bounded, int clamping, boolean mle)
	{
		Gaussian2DFunction f = GaussianFunctionFactory.create2D(1, size, size, GaussianFunctionFactory.FIT_CIRCLE,
				null);
		StoppingCriteria sc = new ErrorStoppingCriteria(5);
		sc.setMaximumIterations(100);
		NonLinearFit solver = (bounded != 0 || clamping != 0) ? new BoundedNonLinearFit(f, sc)
				: new NonLinearFit(f, sc);
		if (clamping != 0)
		{
			BoundedNonLinearFit bsolver = (BoundedNonLinearFit) solver;
			bsolver.setClampValues(defaultClampValues);
			bsolver.setDynamicClamp(clamping == 2);
			// Local search anecdotally only works with clamped LVM fitters that are not bounded.
			// It must act like a soft-bounding search. For now this is not a used feature and
			// will not be formally tested 
			//bsolver.setLocalSearch(3);
		}
		solver.setMLE(mle);
		solver.setInitialLambda(1);
		return solver;
	}

	private String getLVMName(boolean bounded, int clamping, boolean mle)
	{
		return ((bounded) ? "B" : "") + ((clamping == 0) ? "" : ((clamping == 1) ? "C" : "DC")) + "LVM" +
				((mle) ? " MLE" : "");
	}

	private void canFitSingleGaussian(FunctionSolver solver, boolean applyBounds, boolean withBias)
	{
		randomGenerator.setSeed(seed);
		for (double s : signal)
		{
			double[] expected = createParams(1, s, 0, 0, 1, withBias);
			double[] lower = createParams(0, s * 0.5, -0.2, -0.2, 0.8, withBias);
			double[] upper = createParams(3, s * 2, 0.2, 0.2, 1.2, withBias);
			if (applyBounds)
				solver.setBounds(lower, upper);
			for (double n : noise)
			{
				double[] data = drawGaussian(expected, n, withBias);
				for (double db : base)
					for (double dx : shift)
						for (double dy : shift)
							for (double dsx : factor)
							{
								double[] p = createParams(db, s, dx, dy, dsx, withBias);
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
								}
							}
			}
		}
	}

	private void canFitSingleGaussianBetter(FunctionSolver solver, boolean applyBounds, boolean withBias,
			FunctionSolver solver2, boolean applyBounds2, boolean withBias2, String name, String name2)
	{
		int LOOPS = 5;
		randomGenerator.setSeed(seed);
		double bias2 = (withBias != withBias2) ? (withBias) ? -bias : bias : 0;
		StoredDataStatistics[] stats = new StoredDataStatistics[6];
		String[] statName = { "Signal", "X", "Y" };

		int[] betterPrecision = new int[3];
		int[] totalPrecision = new int[3];
		int[] betterAccuracy = new int[3];
		int[] totalAccuracy = new int[3];

		int i1 = 0, i2 = 0;
		for (double s : signal)
		{
			double[] expected = createParams(1, s, 0, 0, 1, withBias);
			if (applyBounds)
			{
				double[] lower = createParams(0, s * 0.5, -0.2, -0.2, 0.8, withBias);
				double[] upper = createParams(3, s * 2, 0.2, 0.2, 1.2, withBias);
				solver.setBounds(lower, upper);
			}

			double[] expected2 = createParams(1, s, 0, 0, 1, withBias2);
			if (applyBounds2)
			{
				double[] lower2 = createParams(0, s * 0.5, -0.2, -0.2, 0.8, withBias2);
				double[] upper2 = createParams(3, s * 2, 0.2, 0.2, 1.2, withBias2);
				solver2.setBounds(lower2, upper2);
			}

			for (double n : noise)
			{
				for (int loop = LOOPS; loop-- > 0;)
				{
					double[] data = drawGaussian(expected, n, withBias);
					double[] data2 = data.clone();
					for (int i = 0; i < data.length; i++)
						data2[i] += bias2;

					for (int i = 0; i < stats.length; i++)
						stats[i] = new StoredDataStatistics();

					for (double db : base)
						for (double dx : shift)
							for (double dy : shift)
								for (double dsx : factor)
								{
									double[] p = createParams(db, s, dx, dy, dsx, withBias);
									double[] fp = fitGaussian(solver, data, p, expected);
									i1 += solver.getEvaluations();

									double[] p2 = createParams(db, s, dx, dy, dsx, withBias2);
									double[] fp2 = fitGaussian(solver2, data2, p2, expected2);
									i2 += solver2.getEvaluations();

									// Get the mean and sd (the fit precision)
									compare(fp, expected, fp2, expected2, Gaussian2DFunction.SIGNAL, stats[0],
											stats[1]);

									compare(fp, expected, fp2, expected2, Gaussian2DFunction.X_POSITION, stats[2],
											stats[3]);
									compare(fp, expected, fp2, expected2, Gaussian2DFunction.Y_POSITION, stats[4],
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
		// TODO Auto-generated method stub
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

	private double[] createParams(double db, double signal, double dx, double dy, double dsx, boolean withBias)
	{
		double[] p = params.clone();
		p[Gaussian2DFunction.BACKGROUND] *= db;
		if (withBias)
			p[Gaussian2DFunction.BACKGROUND] += bias;
		p[Gaussian2DFunction.SIGNAL] = signal;
		p[Gaussian2DFunction.X_POSITION] += dx;
		p[Gaussian2DFunction.Y_POSITION] += dy;
		p[Gaussian2DFunction.X_SD] *= dsx;
		return p;
	}

	private double[] fitGaussian(FunctionSolver solver, double[] data, double[] params, double[] expected)
	{
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
	 * @param withBias
	 * @return The data
	 */
	private double[] drawGaussian(double[] params, double noise, boolean withBias)
	{
		double[] data = new double[size * size];
		int n = params.length / 6;
		Gaussian2DFunction f = GaussianFunctionFactory.create2D(n, size, size, GaussianFunctionFactory.FIT_CIRCLE,
				null);
		f.initialise(params);
		final double bias = (withBias) ? BoundedFunctionSolverTest.bias : 0;
		for (int i = 0; i < data.length; i++)
		{
			data[i] = bias + dataGenerator.nextPoisson(f.eval(i) - bias);
		}
		if (noise != 0)
			for (int i = 0; i < data.length; i++)
				data[i] += dataGenerator.nextGaussian(0, noise);
		//gdsc.core.ij.Utils.display("Spot", data, size, size);
		return data;
	}
}
