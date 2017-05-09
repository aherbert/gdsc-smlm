package gdsc.smlm.fitting.nonlinear.gradient;

import java.util.ArrayList;

import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.ejml.data.DenseMatrix64F;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.utils.DoubleEquality;
import gdsc.smlm.TestSettings;
import gdsc.smlm.function.Gradient1Function;
import gdsc.smlm.function.gaussian.erf.ErfGaussian2DFunction;
import gdsc.smlm.function.gaussian.erf.SingleFreeCircularErfGaussian2DFunction;

/**
 * Contains speed tests for the methods for calculating the Hessian and gradient vector
 * for use in the LVM algorithm.
 */
public class LVMGradientProcedureTest
{
	boolean speedTests = true;
	DoubleEquality eq = new DoubleEquality(6, 1e-16);

	int MAX_ITER = 20000;
	int blockWidth = 10;
	double Background = 0.5;
	double Signal = 100;
	double Angle = Math.PI;
	double Xpos = 5;
	double Ypos = 5;
	double Xwidth = 1.2;
	double Ywidth = 1.2;

	RandomDataGenerator rdg;

	@Test
	public void gradientCalculatorFactoryCreatesOptimisedCalculators()
	{
		double[] y = new double[0];
		Assert.assertEquals(LVMGradientProcedureFactory.create(y, new DummyGradientFunction(6), false).getClass(),
				LSQLVMGradientProcedure6.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y, new DummyGradientFunction(5), false).getClass(),
				LSQLVMGradientProcedure5.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y, new DummyGradientFunction(4), false).getClass(),
				LSQLVMGradientProcedure4.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y, new DummyGradientFunction(6), true).getClass(),
				MLELVMGradientProcedure6.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y, new DummyGradientFunction(5), true).getClass(),
				MLELVMGradientProcedure5.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y, new DummyGradientFunction(4), true).getClass(),
				MLELVMGradientProcedure4.class);
	}

	@Test
	public void gradientProcedureLSQComputesSameAsGradientCalculator()
	{
		gradientProcedureComputesSameAsGradientCalculator(false);
	}

	@Test
	public void gradientProcedureMLEComputesSameAsGradientCalculator()
	{
		gradientProcedureComputesSameAsGradientCalculator(true);
	}

	private void gradientProcedureComputesSameAsGradientCalculator(boolean mle)
	{
		gradientProcedureComputesSameAsGradientCalculator(4, mle);
		gradientProcedureComputesSameAsGradientCalculator(5, mle);
		gradientProcedureComputesSameAsGradientCalculator(6, mle);
		gradientProcedureComputesSameAsGradientCalculator(11, mle);
		gradientProcedureComputesSameAsGradientCalculator(21, mle);
	}

	@Test
	public void gradientProcedureLSQIsNotSlowerThanGradientCalculator()
	{
		gradientProcedureIsNotSlowerThanGradientCalculator(false);
	}

	@Test
	public void gradientProcedureMLEIsNotSlowerThanGradientCalculator()
	{
		gradientProcedureIsNotSlowerThanGradientCalculator(true);
	}

	private void gradientProcedureIsNotSlowerThanGradientCalculator(boolean mle)
	{
		gradientProcedureIsNotSlowerThanGradientCalculator(4, mle);
		gradientProcedureIsNotSlowerThanGradientCalculator(5, mle);
		gradientProcedureIsNotSlowerThanGradientCalculator(6, mle);
		// 2 peaks
		gradientProcedureIsNotSlowerThanGradientCalculator(11, mle);
		// 4 peaks
		gradientProcedureIsNotSlowerThanGradientCalculator(21, mle);
	}

	private void gradientProcedureComputesSameAsGradientCalculator(int nparams, boolean mle)
	{
		int iter = 10;
		rdg = new RandomDataGenerator(new Well19937c(30051977));

		double[][] alpha = new double[nparams][nparams];
		double[] beta = new double[nparams];

		ArrayList<double[]> paramsList = new ArrayList<double[]>(iter);
		ArrayList<double[]> yList = new ArrayList<double[]>(iter);

		int[] x = createFakeData(nparams, iter, paramsList, yList);
		int n = x.length;
		FakeGradientFunction func = new FakeGradientFunction(blockWidth, nparams);

		GradientCalculator calc = GradientCalculatorFactory.newCalculator(nparams, mle);

		String name = String.format("[%d] %b", nparams, mle);

		for (int i = 0; i < paramsList.size(); i++)
		{
			LVMGradientProcedure p = LVMGradientProcedureFactory.create(yList.get(i), func, mle);
			p.gradient(paramsList.get(i));
			double s = p.value;
			double s2 = calc.findLinearised(n, yList.get(i), paramsList.get(i), alpha, beta, func);
			// Exactly the same ...
			Assert.assertEquals(name + " Result: Not same @ " + i, s, s2, 0);
			Assert.assertArrayEquals(name + " Observations: Not same beta @ " + i, p.beta, beta, 0);

			double[] al = p.getAlphaLinear();
			Assert.assertArrayEquals(name + " Observations: Not same alpha @ " + i, al, new DenseMatrix64F(alpha).data,
					0);

			double[][] am = p.getAlphaMatrix();
			for (int j = 0; j < nparams; j++)
				Assert.assertArrayEquals(name + " Observations: Not same alpha @ " + i, am[j], alpha[j], 0);
		}
	}

	private abstract class Timer
	{
		private int loops;
		int min;

		Timer()
		{
		}

		Timer(int min)
		{
			this.min = min;
		}

		long getTime()
		{
			// Run till stable timing
			long t1 = time();
			for (int i = 0; i < 10; i++)
			{
				long t2 = t1;
				t1 = time();
				if (loops >= min && DoubleEquality.relativeError(t1, t2) < 0.02) // 2% difference
					break;
			}
			return t1;
		}

		long time()
		{
			loops++;
			long t = System.nanoTime();
			run();
			t = System.nanoTime() - t;
			//System.out.printf("[%d] Time = %d\n", loops, t);
			return t;
		}

		abstract void run();
	}

	private void gradientProcedureIsNotSlowerThanGradientCalculator(final int nparams, final boolean mle)
	{
		org.junit.Assume.assumeTrue(speedTests || TestSettings.RUN_SPEED_TESTS);

		final int iter = 1000;
		rdg = new RandomDataGenerator(new Well19937c(30051977));
		final double[][] alpha = new double[nparams][nparams];
		final double[] beta = new double[nparams];

		final ArrayList<double[]> paramsList = new ArrayList<double[]>(iter);
		final ArrayList<double[]> yList = new ArrayList<double[]>(iter);

		int[] x = createFakeData(nparams, iter, paramsList, yList);
		final int n = x.length;
		final FakeGradientFunction func = new FakeGradientFunction(blockWidth, nparams);

		GradientCalculator calc = GradientCalculatorFactory.newCalculator(nparams, mle);

		for (int i = 0; i < paramsList.size(); i++)
			calc.findLinearised(n, yList.get(i), paramsList.get(i), alpha, beta, func);

		for (int i = 0; i < paramsList.size(); i++)
		{
			LVMGradientProcedure p = LVMGradientProcedureFactory.create(yList.get(i), func, mle);
			p.gradient(paramsList.get(i));
		}

		// Realistic loops for an optimisation
		final int loops = 15;

		// Run till stable timing
		Timer t1 = new Timer()
		{
			@Override
			void run()
			{
				for (int i = 0, k = 0; i < iter; i++)
				{
					GradientCalculator calc = GradientCalculatorFactory.newCalculator(nparams, mle);
					for (int j = loops; j-- > 0;)
						calc.findLinearised(n, yList.get(i), paramsList.get(k++ % iter), alpha, beta, func);
				}
			}
		};
		long time1 = t1.getTime();

		Timer t2 = new Timer(t1.loops)
		{
			@Override
			void run()
			{
				for (int i = 0, k = 0; i < iter; i++)
				{
					LVMGradientProcedure p = LVMGradientProcedureFactory.create(yList.get(i), func, mle);
					for (int j = loops; j-- > 0;)
						p.gradient(paramsList.get(k++ % iter));
				}
			}
		};
		long time2 = t2.getTime();

		log("GradientCalculator = %d : LVMGradientProcedure %d %b = %d : %fx\n", time1, nparams, mle, time2,
				(1.0 * time1) / time2);
		if (TestSettings.ASSERT_SPEED_TESTS)
		{
			// Add contingency
			Assert.assertTrue(time2 < time1 * 1.5);
		}
	}

	@Test
	public void gradientProcedureLSQUnrolledComputesSameAsGradientProcedure()
	{
		gradientProcedureUnrolledComputesSameAsGradientProcedure(false);
	}

	@Test
	public void gradientProcedureMLEUnrolledComputesSameAsGradientProcedure()
	{
		gradientProcedureUnrolledComputesSameAsGradientProcedure(true);
	}

	private void gradientProcedureUnrolledComputesSameAsGradientProcedure(boolean mle)
	{
		gradientProcedureUnrolledComputesSameAsGradientProcedure(4, mle);
		gradientProcedureUnrolledComputesSameAsGradientProcedure(5, mle);
		gradientProcedureUnrolledComputesSameAsGradientProcedure(6, mle);
	}

	private void gradientProcedureUnrolledComputesSameAsGradientProcedure(int nparams, boolean mle)
	{
		int iter = 10;
		rdg = new RandomDataGenerator(new Well19937c(30051977));

		ArrayList<double[]> paramsList = new ArrayList<double[]>(iter);
		ArrayList<double[]> yList = new ArrayList<double[]>(iter);

		createFakeData(nparams, iter, paramsList, yList);
		FakeGradientFunction func = new FakeGradientFunction(blockWidth, nparams);

		String name = String.format("[%d] %b", nparams, mle);
		for (int i = 0; i < paramsList.size(); i++)
		{
			LVMGradientProcedure p1 = (mle) ? new MLELVMGradientProcedure(yList.get(i), func)
					: new LSQLVMGradientProcedure(yList.get(i), func);
			p1.gradient(paramsList.get(i));

			LVMGradientProcedure p2 = LVMGradientProcedureFactory.create(yList.get(i), func, mle);
			p2.gradient(paramsList.get(i));

			// Exactly the same ...
			Assert.assertEquals(name + " Result: Not same @ " + i, p1.value, p2.value, 0);
			Assert.assertArrayEquals(name + " Observations: Not same beta @ " + i, p1.beta, p2.beta, 0);

			Assert.assertArrayEquals(name + " Observations: Not same alpha @ " + i, p1.getAlphaLinear(),
					p2.getAlphaLinear(), 0);

			double[][] am1 = p1.getAlphaMatrix();
			double[][] am2 = p2.getAlphaMatrix();
			for (int j = 0; j < nparams; j++)
				Assert.assertArrayEquals(name + " Observations: Not same alpha @ " + i, am1[j], am2[j], 0);
		}
	}

	@Test
	public void gradientProcedureLSQIsFasterUnrolledThanGradientProcedure()
	{
		gradientProcedureIsFasterUnrolledThanGradientProcedure(false);
	}

	@Test
	public void gradientProcedureMLEIsFasterUnrolledThanGradientProcedure()
	{
		gradientProcedureIsFasterUnrolledThanGradientProcedure(true);
	}

	private void gradientProcedureIsFasterUnrolledThanGradientProcedure(boolean mle)
	{
		gradientProcedureLinearIsFasterThanGradientProcedure(4, mle);
		gradientProcedureLinearIsFasterThanGradientProcedure(5, mle);
		gradientProcedureLinearIsFasterThanGradientProcedure(6, mle);
	}

	private void gradientProcedureLinearIsFasterThanGradientProcedure(final int nparams, final boolean mle)
	{
		org.junit.Assume.assumeTrue(speedTests || TestSettings.RUN_SPEED_TESTS);

		final int iter = 100;
		rdg = new RandomDataGenerator(new Well19937c(30051977));

		final ArrayList<double[]> paramsList = new ArrayList<double[]>(iter);
		final ArrayList<double[]> yList = new ArrayList<double[]>(iter);

		createData(1, iter, paramsList, yList);

		// Remove the timing of the function call by creating a dummy function
		final Gradient1Function func = new FakeGradientFunction(blockWidth, nparams);

		for (int i = 0; i < paramsList.size(); i++)
		{
			LVMGradientProcedure p1 = (mle) ? new MLELVMGradientProcedure(yList.get(i), func)
					: new LSQLVMGradientProcedure(yList.get(i), func);
			p1.gradient(paramsList.get(i));
			p1.gradient(paramsList.get(i));

			LVMGradientProcedure p2 = LVMGradientProcedureFactory.create(yList.get(i), func, mle);
			p2.gradient(paramsList.get(i));
			p2.gradient(paramsList.get(i));

			// Check they are the same
			Assert.assertArrayEquals("A " + i, p1.getAlphaLinear(), p2.getAlphaLinear(), 0);
			Assert.assertArrayEquals("B " + i, p1.beta, p2.beta, 0);
		}

		// Realistic loops for an optimisation
		final int loops = 15;

		// Run till stable timing
		Timer t1 = new Timer()
		{
			@Override
			void run()
			{
				for (int i = 0, k = 0; i < paramsList.size(); i++)
				{
					LVMGradientProcedure p1 = (mle) ? new MLELVMGradientProcedure(yList.get(i), func)
							: new LSQLVMGradientProcedure(yList.get(i), func);
					for (int j = loops; j-- > 0;)
						p1.gradient(paramsList.get(k++ % iter));
				}
			}
		};
		long time1 = t1.getTime();

		Timer t2 = new Timer(t1.loops)
		{
			@Override
			void run()
			{
				for (int i = 0, k = 0; i < paramsList.size(); i++)
				{
					LVMGradientProcedure p2 = LVMGradientProcedureFactory.create(yList.get(i), func, mle);
					for (int j = loops; j-- > 0;)
						p2.gradient(paramsList.get(k++ % iter));
				}
			}
		};
		long time2 = t2.getTime();

		log("%s = %d : %s %d = %d : %fx\n", mle, time1, mle, nparams, time2, (1.0 * time1) / time2);
		Assert.assertTrue(time2 < time1);
	}

	@Test
	public void gradientCalculatorLSQComputesGradient()
	{
		gradientCalculatorComputesGradient(new SingleFreeCircularErfGaussian2DFunction(blockWidth, blockWidth), false);
	}

	@Test
	public void gradientCalculatorMLEComputesGradient()
	{
		gradientCalculatorComputesGradient(new SingleFreeCircularErfGaussian2DFunction(blockWidth, blockWidth), true);
	}

	private void gradientCalculatorComputesGradient(ErfGaussian2DFunction func, boolean mle)
	{
		int nparams = func.getNumberOfGradients();
		int[] indices = func.gradientIndices();

		int iter = 100;
		rdg = new RandomDataGenerator(new Well19937c(30051977));

		ArrayList<double[]> paramsList = new ArrayList<double[]>(iter);
		ArrayList<double[]> yList = new ArrayList<double[]>(iter);

		createData(1, iter, paramsList, yList, true);

		double delta = 1e-3;
		DoubleEquality eq = new DoubleEquality(3, 1e-3);

		for (int i = 0; i < paramsList.size(); i++)
		{
			double[] y = yList.get(i);
			double[] a = paramsList.get(i);
			double[] a2 = a.clone();
			LVMGradientProcedure p = LVMGradientProcedureFactory.create(y, func, mle);
			p.gradient(a);
			//double s = p.value;
			double[] beta = p.beta.clone();
			for (int j = 0; j < nparams; j++)
			{
				int k = indices[j];
				double d = (a[k] == 0) ? 1e-3 : a[k] * delta;
				a2[k] = a[k] + d;
				p.value(a2);
				double s1 = p.value;
				a2[k] = a[k] - d;
				p.value(a2);
				double s2 = p.value;
				a2[k] = a[k];

				// Apply a factor of -2 to compute the actual gradients:
				// See Numerical Recipes in C++, 2nd Ed. Equation 15.5.6 for Nonlinear Models
				beta[j] *= -2;

				double gradient = (s1 - s2) / (2 * d);
				//System.out.printf("[%d,%d] %f  (%s %f+/-%f)  %f  ?=  %f\n", i, k, s, func.getName(k), a[k], d, beta[j],
				//		gradient);
				Assert.assertTrue("Not same gradient @ " + j, eq.almostEqualComplement(beta[j], gradient));
			}
		}
	}

	/**
	 * Create random elliptical Gaussian data an returns the data plus an estimate of the parameters.
	 * Only the chosen parameters are randomised and returned for a maximum of (background, amplitude, angle, xpos,
	 * ypos, xwidth, ywidth }
	 *
	 * @param npeaks
	 *            the npeaks
	 * @param params
	 *            set on output
	 * @param randomiseParams
	 *            Set to true to randomise the params
	 * @return the double[]
	 */
	private double[] doubleCreateGaussianData(int npeaks, double[] params, boolean randomiseParams)
	{
		int n = blockWidth * blockWidth;

		// Generate a 2D Gaussian
		SingleFreeCircularErfGaussian2DFunction func = new SingleFreeCircularErfGaussian2DFunction(blockWidth,
				blockWidth);
		params[0] = random(Background);
		for (int i = 0, j = 1; i < npeaks; i++, j += 6)
		{
			params[j] = random(Signal);
			params[j + 2] = random(Xpos);
			params[j + 3] = random(Ypos);
			params[j + 4] = random(Xwidth);
			params[j + 5] = random(Ywidth);
		}

		double[] y = new double[n];
		func.initialise(params);
		for (int i = 0; i < y.length; i++)
		{
			// Add random Poisson noise
			y[i] = rdg.nextPoisson(func.eval(i));
		}

		if (randomiseParams)
		{
			params[0] = random(params[0]);
			for (int i = 0, j = 1; i < npeaks; i++, j += 6)
			{
				params[j] = random(params[j]);
				params[j + 2] = random(params[j + 2]);
				params[j + 3] = random(params[j + 3]);
				params[j + 4] = random(params[j + 4]);
				params[j + 5] = random(params[j + 5]);
			}
		}

		return y;
	}

	private double random(double d)
	{
		return d + rdg.nextUniform(-d * 0.1, d * 0.1);
	}

	protected int[] createData(int npeaks, int iter, ArrayList<double[]> paramsList, ArrayList<double[]> yList)
	{
		return createData(npeaks, iter, paramsList, yList, true);
	}

	protected int[] createData(int npeaks, int iter, ArrayList<double[]> paramsList, ArrayList<double[]> yList,
			boolean randomiseParams)
	{
		int[] x = new int[blockWidth * blockWidth];
		for (int i = 0; i < x.length; i++)
			x[i] = i;
		for (int i = 0; i < iter; i++)
		{
			double[] params = new double[1 + 6 * npeaks];
			double[] y = doubleCreateGaussianData(npeaks, params, randomiseParams);
			paramsList.add(params);
			yList.add(y);
		}
		return x;
	}

	protected int[] createFakeData(int nparams, int iter, ArrayList<double[]> paramsList, ArrayList<double[]> yList)
	{
		int[] x = new int[blockWidth * blockWidth];
		for (int i = 0; i < x.length; i++)
			x[i] = i;
		for (int i = 0; i < iter; i++)
		{
			double[] params = new double[nparams];
			double[] y = createFakeData(params);
			paramsList.add(params);
			yList.add(y);
		}
		return x;
	}

	private double[] createFakeData(double[] params)
	{
		int n = blockWidth * blockWidth;
		RandomGenerator r = rdg.getRandomGenerator();

		for (int i = 0; i < params.length; i++)
		{
			params[i] = r.nextDouble();
		}

		double[] y = new double[n];
		for (int i = 0; i < y.length; i++)
		{
			y[i] = r.nextDouble();
		}

		return y;
	}

	protected ArrayList<double[]> copyList(ArrayList<double[]> paramsList)
	{
		ArrayList<double[]> params2List = new ArrayList<double[]>(paramsList.size());
		for (int i = 0; i < paramsList.size(); i++)
		{
			params2List.add(copydouble(paramsList.get(i)));
		}
		return params2List;
	}

	private double[] copydouble(double[] d)
	{
		double[] d2 = new double[d.length];
		for (int i = 0; i < d.length; i++)
			d2[i] = d[i];
		return d2;
	}

	void log(String format, Object... args)
	{
		System.out.printf(format, args);
	}
}
