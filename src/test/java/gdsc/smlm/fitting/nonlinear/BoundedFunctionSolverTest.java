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
	static double[] signal = { 1000, 5000, 10000 }; // 100, 200, 400, 800 };
	static double[] noise = { 0.1, 0.5, 1 };
	static double[] shift = { -1, 0, 1 };
	static double[] factor = { 0.7, 1, 1.3 };
	static int size = 15;
	static
	{
		params[Gaussian2DFunction.BACKGROUND] = 1;
		params[Gaussian2DFunction.X_POSITION] = size / 2;
		params[Gaussian2DFunction.Y_POSITION] = size / 2;
		params[Gaussian2DFunction.X_SD] = 1.4;
	}

	// TODO - test the Clamping of the LVM fitter
	
	
	
	@Test
	public void canFitSingleGaussianLVM()
	{
		fitSingleGaussianBetterLVM(false, false, false);
	}

	@Test
	public void canFitSingleGaussianBLVMNoBounds()
	{
		fitSingleGaussianBetterLVM(true, false, false);
	}

	@Test
	public void canFitSingleGaussianBLVM()
	{
		fitSingleGaussianBetterLVM(true, true, false);
	}

	@Test
	public void canFitSingleGaussianLVMMLE()
	{
		fitSingleGaussianBetterLVM(false, false, true);
	}

	@Test
	public void canFitSingleGaussianBLVMMLENoBounds()
	{
		fitSingleGaussianBetterLVM(true, false, true);
	}

	@Test
	public void canFitSingleGaussianBLVMMLE()
	{
		fitSingleGaussianBetterLVM(true, true, true);
	}

	private void fitSingleGaussianBetterLVM(boolean bounded, boolean applyBounds, boolean mle)
	{
		Gaussian2DFunction f = GaussianFunctionFactory.create2D(1, size, GaussianFunctionFactory.FIT_CIRCLE);
		StoppingCriteria sc = new ErrorStoppingCriteria();
		sc.setMaximumIterations(100);
		NonLinearFit solver = (bounded) ? new BoundedNonLinearFit(f, sc) : new NonLinearFit(f, sc);
		solver.setMLE(mle);
		canFitSingleGaussian(solver, applyBounds, !mle);
	}
	
	@Test
	public void fitSingleGaussianBLVMBetterThanLVM()
	{
		fitSingleGaussianBetterLVM(true, false, false, false);
	}

	@Test
	public void fitSingleGaussianLVMMLEBetterThanLVM()
	{
		fitSingleGaussianBetterLVM(false, true, false, false);
	}

	@Test
	public void fitSingleGaussianBLVMMLEBetterThanLVM()
	{
		fitSingleGaussianBetterLVM(true, true, false, false);
	}

	@Test
	public void fitSingleGaussianBLVMMLEBetterThanLVMMLE()
	{
		fitSingleGaussianBetterLVM(true, true, false, true);
	}

	@Test
	public void fitSingleGaussianBLVMMLEBetterThanBLVM()
	{
		fitSingleGaussianBetterLVM(true, true, true, false);
	}

	private void fitSingleGaussianBetterLVM(boolean bounded2, boolean mle2, boolean bounded, boolean mle)
	{
		Gaussian2DFunction f = GaussianFunctionFactory.create2D(1, size, GaussianFunctionFactory.FIT_CIRCLE);
		StoppingCriteria sc = new ErrorStoppingCriteria();
		sc.setMaximumIterations(100);
		BoundedNonLinearFit solver = new BoundedNonLinearFit(f, sc);
		solver.setMLE(mle);
		BoundedNonLinearFit solver2 = new BoundedNonLinearFit(f, sc);
		solver2.setMLE(mle2);
		canFitSingleGaussianBetter(solver, bounded, !mle, solver2, bounded2, !mle2, getLVMName(bounded, mle),
				getLVMName(bounded2, mle2));
	}

	private String getLVMName(boolean bounded, boolean mle)
	{
		return ((bounded) ? "B" : "") + "LVM" + ((mle) ? " MLE" : "");
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
		randomGenerator.setSeed(seed);
		int count = 0;
		int better = 0;
		double bias2 = (withBias != withBias2) ? (withBias) ? -bias : bias : 0;
		for (double s : signal)
		{
			double[] expected = createParams(1, s, 0, 0, 1, withBias);
			double[] lower = createParams(0, s * 0.5, -0.2, -0.2, 0.8, withBias);
			double[] upper = createParams(3, s * 2, 0.2, 0.2, 1.2, withBias);

			double[] expected2 = createParams(1, s, 0, 0, 1, withBias2);
			double[] lower2 = createParams(0, s * 0.5, -0.2, -0.2, 0.8, withBias2);
			double[] upper2 = createParams(3, s * 2, 0.2, 0.2, 1.2, withBias2);

			for (double n : noise)
			{
				double[] data = drawGaussian(expected, n, withBias);
				double[] data2 = data.clone();
				for (int i = 0; i < data.length; i++)
					data2[i] += bias2;

				for (double db : base)
					for (double dx : shift)
						for (double dy : shift)
							for (double dsx : factor)
							{
								double[] p = createParams(db, s, dx, dy, dsx, withBias);
								if (applyBounds)
									solver.setBounds(lower, upper);
								double[] fp = fitGaussian(solver, data, p, expected);
								if (applyBounds)
									solver.setBounds(null, null);

								double[] p2 = createParams(db, s, dx, dy, dsx, withBias2);
								if (applyBounds2)
									solver2.setBounds(lower2, upper2);
								double[] fp2 = fitGaussian(solver2, data2, p2, expected2);
								if (applyBounds2)
									solver2.setBounds(null, null);

								for (int i = 0; i < upper.length; i++)
								{ // Ignore parameters we are not fitting
									if (upper[i] == 0)
										continue;
									double d1 = Math.abs(fp[i] - expected[i]);
									double d2 = Math.abs(fp2[i] - expected2[i]);
									count++;
									if (d2 <= d1)
										better++;
								}
							}
			}
		}
		String msg = String.format("%s vs %s : Better %d / %d", name2, name, better, count);
		System.out.println(msg);
		Assert.assertTrue(msg, better > count / 2);
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
	 * @param withBias
	 * @return The data
	 */
	private double[] drawGaussian(double[] params, double noise, boolean withBias)
	{
		double[] data = new double[size * size];
		int n = params.length / 6;
		Gaussian2DFunction f = GaussianFunctionFactory.create2D(n, size, GaussianFunctionFactory.FIT_CIRCLE);
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
