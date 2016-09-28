package gdsc.smlm.fitting.nonlinear;

import java.util.Arrays;

import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.junit.Assert;
import org.junit.Test;

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
	RandomGenerator randomGenerator = new Well19937c(seed);
	RandomDataGenerator dataGenerator = new RandomDataGenerator(randomGenerator);

	// Basic Gaussian
	static double bias = 100;
	static double[] params = new double[7];
	static double[] base = { 0.9, 1, 1.1 };
	static double[] signal = { 5000, 10000 }; // 100, 200, 400, 800 };
	static double[] noise = { 0.1, 0.5, 1 };
	static double[] shift = { -1, 0, 1 };
	static double[] factor = { 0.7, 1, 1.3 };
	static int size = 15;
	static
	{
		params[Gaussian2DFunction.BACKGROUND] = 0.4;
		params[Gaussian2DFunction.X_POSITION] = size / 2;
		params[Gaussian2DFunction.Y_POSITION] = size / 2;
		params[Gaussian2DFunction.X_SD] = 1.4;
	}

	@Test
	public void nonLinearFitCanFitSingleGaussian()
	{
		Gaussian2DFunction f = GaussianFunctionFactory.create2D(1, size, GaussianFunctionFactory.FIT_CIRCLE);
		StoppingCriteria sc = new ErrorStoppingCriteria();
		sc.setMaximumIterations(100);
		FunctionSolver solver = new NonLinearFit(f, sc);
		canFitSingleGaussian(solver, false);
	}

	@Test
	public void boundedNonLinearFitCanFitSingleGaussian()
	{
		Gaussian2DFunction f = GaussianFunctionFactory.create2D(1, size, GaussianFunctionFactory.FIT_CIRCLE);
		StoppingCriteria sc = new ErrorStoppingCriteria();
		sc.setMaximumIterations(100);
		FunctionSolver solver = new BoundedNonLinearFit(f, sc);
		canFitSingleGaussian(solver, false);
	}

	@Test
	public void boundedNonLinearFitCanFitSingleGaussianWithBounds()
	{
		Gaussian2DFunction f = GaussianFunctionFactory.create2D(1, size, GaussianFunctionFactory.FIT_CIRCLE);
		StoppingCriteria sc = new ErrorStoppingCriteria();
		sc.setMaximumIterations(100);
		FunctionSolver solver = new BoundedNonLinearFit(f, sc);
		canFitSingleGaussian(solver, true);
	}

	@Test
	public void boundedNonLinearFitCanFitSingleGaussianBetterWithBounds()
	{
		Gaussian2DFunction f = GaussianFunctionFactory.create2D(1, size, GaussianFunctionFactory.FIT_CIRCLE);
		StoppingCriteria sc = new ErrorStoppingCriteria();
		sc.setMaximumIterations(100);
		FunctionSolver solver = new BoundedNonLinearFit(f, sc);
		canFitSingleGaussianBetter(solver);
	}

	private void canFitSingleGaussian(FunctionSolver solver, boolean applyBounds)
	{
		randomGenerator.setSeed(seed);
		for (double s : signal)
		{
			double[] expected = createParams(1, s, 0, 0, 1);
			double[] lower = createParams(0, s * 0.5, -0.2, -0.2, 0.8);
			double[] upper = createParams(3, s * 2, 0.2, 0.2, 1.2);
			if (applyBounds)
				solver.setBounds(lower, upper);
			for (double n : noise)
			{
				double[] data = drawGaussian(expected, n);
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
								}
							}
			}
		}
	}

	private void canFitSingleGaussianBetter(FunctionSolver solver)
	{
		randomGenerator.setSeed(seed);
		int count = 0;
		int better = 0;
		for (double s : signal)
		{
			double[] expected = createParams(1, s, 0, 0, 1);
			double[] lower = createParams(0, s * 0.5, -0.2, -0.2, 0.8);
			double[] upper = createParams(3, s * 2, 0.2, 0.2, 1.2);
			for (double n : noise)
			{
				double[] data = drawGaussian(expected, n);
				for (double db : base)
					for (double dx : shift)
						for (double dy : shift)
							for (double dsx : factor)
							{
								double[] p = createParams(db, s, dx, dy, dsx);
								double[] fp = fitGaussian(solver, data, p, expected);

								solver.setBounds(lower, upper);
								double[] fp2 = fitGaussian(solver, data, p, expected);
								solver.setBounds(null, null);

								for (int i = 0; i < upper.length; i++)
								{ // Ignore parameters we are not fitting
									if (upper[i] == 0)
										continue;
									double d1 = Math.abs(fp[i] - expected[i]);
									double d2 = Math.abs(fp2[i] - expected[i]);
									count++;
									if (d2 <= d1)
										better++;
								}
							}
			}
		}
		String msg = String.format("Better %d / %d", better, count);
		System.out.println(msg);
		Assert.assertTrue(msg, better > count / 2);
	}

	private double[] createParams(double db, double signal, double dx, double dy, double dsx)
	{
		double[] p = params.clone();
		p[Gaussian2DFunction.BACKGROUND] *= db;
		p[Gaussian2DFunction.BACKGROUND] += bias;
		p[Gaussian2DFunction.SIGNAL] = signal;
		p[Gaussian2DFunction.X_POSITION] += dx;
		p[Gaussian2DFunction.Y_POSITION] += dy;
		p[Gaussian2DFunction.X_SD] *= dsx;
		return p;
	}

	private double[] fitGaussian(FunctionSolver solver, double[] data, double[] params, double[] expected)
	{
		double[] error = new double[1];
		params = params.clone();
		FitStatus status = solver.fit(data.length, data, null, params, null, error, 0);
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
	 * @return The data
	 */
	private double[] drawGaussian(double[] params, double noise)
	{
		double[] data = new double[size * size];
		int n = params.length / 6;
		Gaussian2DFunction f = GaussianFunctionFactory.create2D(n, size, GaussianFunctionFactory.FIT_CIRCLE);
		f.initialise(params);
		for (int i = 0; i < data.length; i++)
		{
			data[i] = bias + dataGenerator.nextPoisson(f.eval(i) - bias);
		}
		if (noise != 0)
			for (int i = 0; i < data.length; i++)
				data[i] += dataGenerator.nextGaussian(0, noise);
		gdsc.core.ij.Utils.display("Spot", data, size, size);
		return data;
	}
}
