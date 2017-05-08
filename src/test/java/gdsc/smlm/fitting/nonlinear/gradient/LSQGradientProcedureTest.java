package gdsc.smlm.fitting.nonlinear.gradient;

import java.util.ArrayList;

import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.ejml.data.DenseMatrix64F;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.utils.DoubleEquality;
import gdsc.core.utils.Statistics;
import gdsc.smlm.TestSettings;
import gdsc.smlm.fitting.linear.EJMLLinearSolver;
import gdsc.smlm.function.Gradient1Function;
import gdsc.smlm.function.gaussian.erf.ErfGaussian2DFunction;
import gdsc.smlm.function.gaussian.erf.SingleFreeCircularErfGaussian2DFunction;

/**
 * Contains speed tests for the fastest method for calculating the Hessian and gradient vector
 * for use in NonLinearFit
 * <p>
 * Note: This class is a test-bed for implementation strategies. The fastest strategy can then be used for other
 * gradient procedures.
 */
public class LSQGradientProcedureTest
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
		Assert.assertEquals(LSQGradientProcedureMatrixFactory.create(y, new DummyGradientFunction(6)).getClass(),
				LSQGradientProcedureMatrix6.class);
		Assert.assertEquals(LSQGradientProcedureMatrixFactory.create(y, new DummyGradientFunction(5)).getClass(),
				LSQGradientProcedureMatrix5.class);
		Assert.assertEquals(LSQGradientProcedureMatrixFactory.create(y, new DummyGradientFunction(4)).getClass(),
				LSQGradientProcedureMatrix4.class);

		Assert.assertEquals(LSQGradientProcedureLinearFactory.create(y, new DummyGradientFunction(6)).getClass(),
				LSQGradientProcedureLinear6.class);
		Assert.assertEquals(LSQGradientProcedureLinearFactory.create(y, new DummyGradientFunction(5)).getClass(),
				LSQGradientProcedureLinear5.class);
		Assert.assertEquals(LSQGradientProcedureLinearFactory.create(y, new DummyGradientFunction(4)).getClass(),
				LSQGradientProcedureLinear4.class);

		Assert.assertEquals(LSQGradientProcedureFactory.create(y, new DummyGradientFunction(6)).getClass(),
				LSQGradientProcedure6.class);
		Assert.assertEquals(LSQGradientProcedureFactory.create(y, new DummyGradientFunction(5)).getClass(),
				LSQGradientProcedure5.class);
		Assert.assertEquals(LSQGradientProcedureFactory.create(y, new DummyGradientFunction(4)).getClass(),
				LSQGradientProcedure4.class);
	}

	@Test
	public void gradientProcedureLinearComputesSameAsGradientCalculator()
	{
		gradientProcedureComputesSameAsGradientCalculator(new LSQGradientProcedureLinearFactory());
	}

	@Test
	public void gradientProcedureMatrixComputesSameAsGradientCalculator()
	{
		gradientProcedureComputesSameAsGradientCalculator(new LSQGradientProcedureMatrixFactory());
	}

	@Test
	public void gradientProcedureComputesSameAsGradientCalculator()
	{
		gradientProcedureComputesSameAsGradientCalculator(new LSQGradientProcedureFactory());
	}

	private void gradientProcedureComputesSameAsGradientCalculator(BaseLSQGradientProcedureFactory factory)
	{
		gradientProcedureComputesSameAsGradientCalculator(4, factory);
		gradientProcedureComputesSameAsGradientCalculator(5, factory);
		gradientProcedureComputesSameAsGradientCalculator(6, factory);
		gradientProcedureComputesSameAsGradientCalculator(11, factory);
		gradientProcedureComputesSameAsGradientCalculator(21, factory);
	}

	@Test
	public void gradientProcedureLinearIsNotSlowerThanGradientCalculator()
	{
		gradientProcedureIsNotSlowerThanGradientCalculator(new LSQGradientProcedureLinearFactory());
	}

	@Test
	public void gradientProcedureMatrixIsNotSlowerThanGradientCalculator()
	{
		gradientProcedureIsNotSlowerThanGradientCalculator(new LSQGradientProcedureMatrixFactory());
	}

	@Test
	public void gradientProcedureIsNotSlowerThanGradientCalculator()
	{
		gradientProcedureIsNotSlowerThanGradientCalculator(new LSQGradientProcedureFactory());
	}

	private void gradientProcedureIsNotSlowerThanGradientCalculator(BaseLSQGradientProcedureFactory factory)
	{
		gradientProcedureIsNotSlowerThanGradientCalculator(4, factory);
		gradientProcedureIsNotSlowerThanGradientCalculator(5, factory);
		gradientProcedureIsNotSlowerThanGradientCalculator(6, factory);
		// 2 peaks
		gradientProcedureIsNotSlowerThanGradientCalculator(11, factory);
		// 4 peaks
		gradientProcedureIsNotSlowerThanGradientCalculator(21, factory);
	}

	private void gradientProcedureComputesSameAsGradientCalculator(int nparams, BaseLSQGradientProcedureFactory factory)
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

		GradientCalculator calc = GradientCalculatorFactory.newCalculator(nparams, false);

		String name = factory.getClass().getSimpleName();
		for (int i = 0; i < paramsList.size(); i++)
		{
			BaseLSQGradientProcedure p = factory.createProcedure(yList.get(i), func);
			p.gradient(paramsList.get(i));
			double s = p.ssx;
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

	private void gradientProcedureIsNotSlowerThanGradientCalculator(final int nparams,
			final BaseLSQGradientProcedureFactory factory)
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

		GradientCalculator calc = GradientCalculatorFactory.newCalculator(nparams, false);

		for (int i = 0; i < paramsList.size(); i++)
			calc.findLinearised(n, yList.get(i), paramsList.get(i), alpha, beta, func);

		for (int i = 0; i < paramsList.size(); i++)
		{
			BaseLSQGradientProcedure p = factory.createProcedure(yList.get(i), func);
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
					GradientCalculator calc = GradientCalculatorFactory.newCalculator(nparams, false);
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
					BaseLSQGradientProcedure p = factory.createProcedure(yList.get(i), func);
					for (int j = loops; j-- > 0;)
						p.gradient(paramsList.get(k++ % iter));
				}
			}
		};
		long time2 = t2.getTime();

		log("GradientCalculator = %d : %s %d = %d : %fx\n", time1, factory.getClass().getSimpleName(), nparams, time2,
				(1.0 * time1) / time2);
		if (TestSettings.ASSERT_SPEED_TESTS)
		{
			// Add contingency
			Assert.assertTrue(time2 < time1 * 1.5);
		}
	}

	@Test
	public void gradientProcedureUnrolledComputesSameAsGradientProcedure()
	{
		// Test the method that will be used for the standard and unrolled versions
		// for all other 'gradient procedures'
		gradientProcedureUnrolledComputesSameAsGradientProcedure(4);
		gradientProcedureUnrolledComputesSameAsGradientProcedure(5);
		gradientProcedureUnrolledComputesSameAsGradientProcedure(6);
	}

	private void gradientProcedureUnrolledComputesSameAsGradientProcedure(int nparams)
	{
		int iter = 10;
		rdg = new RandomDataGenerator(new Well19937c(30051977));

		ArrayList<double[]> paramsList = new ArrayList<double[]>(iter);
		ArrayList<double[]> yList = new ArrayList<double[]>(iter);

		createFakeData(nparams, iter, paramsList, yList);
		FakeGradientFunction func = new FakeGradientFunction(blockWidth, nparams);

		String name = GradientCalculator.class.getSimpleName();
		for (int i = 0; i < paramsList.size(); i++)
		{
			BaseLSQGradientProcedure p1 = LSQGradientProcedureFactory.create(yList.get(i), func);
			p1.gradient(paramsList.get(i));

			BaseLSQGradientProcedure p2 = new LSQGradientProcedure(yList.get(i), func);
			p2.gradient(paramsList.get(i));

			// Exactly the same ...
			Assert.assertEquals(name + " Result: Not same @ " + i, p1.ssx, p2.ssx, 0);
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
	public void gradientProcedureIsFasterUnrolledThanGradientProcedureMatrix()
	{
		gradientProcedure2IsFasterUnrolledThanGradientProcedure1(new LSQGradientProcedureMatrixFactory(),
				new LSQGradientProcedureFactory());
	}

	@Test
	public void gradientProcedureLinearIsFasterUnrolledThanGradientProcedureMatrix()
	{
		gradientProcedure2IsFasterUnrolledThanGradientProcedure1(new LSQGradientProcedureMatrixFactory(),
				new LSQGradientProcedureLinearFactory());
	}

	private void gradientProcedure2IsFasterUnrolledThanGradientProcedure1(BaseLSQGradientProcedureFactory factory1,
			BaseLSQGradientProcedureFactory factory2)
	{
		// Assert the unrolled versions
		gradientProcedureLinearIsFasterThanGradientProcedureMatrix(4, factory1, factory2, true);
		gradientProcedureLinearIsFasterThanGradientProcedureMatrix(5, factory1, factory2, true);
		gradientProcedureLinearIsFasterThanGradientProcedureMatrix(6, factory1, factory2, true);
		gradientProcedureLinearIsFasterThanGradientProcedureMatrix(11, factory1, factory2, false);
		gradientProcedureLinearIsFasterThanGradientProcedureMatrix(21, factory1, factory2, false);
	}

	private void gradientProcedureLinearIsFasterThanGradientProcedureMatrix(final int nparams,
			final BaseLSQGradientProcedureFactory factory1, final BaseLSQGradientProcedureFactory factory2,
			boolean doAssert)
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
			BaseLSQGradientProcedure p = factory1.createProcedure(yList.get(i), func);
			p.gradient(paramsList.get(i));
			p.gradient(paramsList.get(i));

			BaseLSQGradientProcedure p2 = factory2.createProcedure(yList.get(i), func);
			p2.gradient(paramsList.get(i));
			p2.gradient(paramsList.get(i));

			// Check they are the same
			Assert.assertArrayEquals("A " + i, p.getAlphaLinear(), p2.getAlphaLinear(), 0);
			Assert.assertArrayEquals("B " + i, p.beta, p2.beta, 0);
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
					BaseLSQGradientProcedure p = factory1.createProcedure(yList.get(i), func);
					for (int j = loops; j-- > 0;)
						p.gradient(paramsList.get(k++ % iter));
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
					BaseLSQGradientProcedure p2 = factory2.createProcedure(yList.get(i), func);
					for (int j = loops; j-- > 0;)
						p2.gradient(paramsList.get(k++ % iter));
				}
			}
		};
		long time2 = t2.getTime();

		log("%s = %d : %s %d = %d : %fx\n", factory1.getClass().getSimpleName(), time1,
				factory2.getClass().getSimpleName(), nparams, time2, (1.0 * time1) / time2);
		if (doAssert)
			Assert.assertTrue(time2 < time1);
	}

	@Test
	public void gradientCalculatorComputesGradient()
	{
		gradientCalculatorComputesGradient(new SingleFreeCircularErfGaussian2DFunction(blockWidth, blockWidth));
	}

	private void gradientCalculatorComputesGradient(ErfGaussian2DFunction func)
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
			BaseLSQGradientProcedure p = LSQGradientProcedureFactory.create(y, func);
			p.gradient(a);
			//double s = p.ssx;
			double[] beta = p.beta.clone();
			for (int j = 0; j < nparams; j++)
			{
				int k = indices[j];
				double d = (a[k] == 0) ? 1e-3 : a[k] * delta;
				a2[k] = a[k] + d;
				p.value(a2);
				double s1 = p.ssx;
				a2[k] = a[k] - d;
				p.value(a2);
				double s2 = p.ssx;
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

	@Test
	public void gradientCalculatorComputesSameOutputWithBias()
	{
		ErfGaussian2DFunction func = new SingleFreeCircularErfGaussian2DFunction(blockWidth, blockWidth);
		int nparams = func.getNumberOfGradients();

		int iter = 100;
		rdg = new RandomDataGenerator(new Well19937c(30051977));

		ArrayList<double[]> paramsList = new ArrayList<double[]>(iter);
		ArrayList<double[]> yList = new ArrayList<double[]>(iter);

		ArrayList<double[]> alphaList = new ArrayList<double[]>(iter);
		ArrayList<double[]> betaList = new ArrayList<double[]>(iter);
		ArrayList<double[]> xList = new ArrayList<double[]>(iter);

		// Manipulate the background
		double defaultBackground = Background;
		try
		{
			Background = 1e-2;
			createData(1, iter, paramsList, yList, true);

			EJMLLinearSolver solver = new EJMLLinearSolver(5, 1e-6);

			for (int i = 0; i < paramsList.size(); i++)
			{
				double[] y = yList.get(i);
				double[] a = paramsList.get(i);
				BaseLSQGradientProcedure p = LSQGradientProcedureFactory.create(y, func);
				p.gradient(a);
				double[] beta = p.beta;
				alphaList.add(p.getAlphaLinear());
				betaList.add(beta.clone());
				for (int j = 0; j < nparams; j++)
				{
					if (Math.abs(beta[j]) < 1e-6)
						System.out.printf("[%d] Tiny beta %s %g\n", i, func.getName(j), beta[j]);
				}
				// Solve
				if (!solver.solve(p.getAlphaMatrix(), beta))
					throw new AssertionError();
				xList.add(beta);
				//System.out.println(Arrays.toString(beta));
			}

			//for (int b = 1; b < 1000; b *= 2)
			for (double b : new double[] { -500, -100, -10, -1, -0.1, 0, 0.1, 1, 10, 100, 500 })
			{
				Statistics[] rel = new Statistics[nparams];
				Statistics[] abs = new Statistics[nparams];
				for (int i = 0; i < nparams; i++)
				{
					rel[i] = new Statistics();
					abs[i] = new Statistics();
				}

				for (int i = 0; i < paramsList.size(); i++)
				{
					double[] y = add(yList.get(i), b);
					double[] a = paramsList.get(i).clone();
					a[0] += b;
					BaseLSQGradientProcedure p = LSQGradientProcedureFactory.create(y, func);
					p.gradient(a);
					double[] beta = p.beta;
					double[] alpha2 = alphaList.get(i);
					double[] beta2 = betaList.get(i);
					double[] x2 = xList.get(i);

					Assert.assertArrayEquals("Beta", beta2, beta, 1e-10);
					Assert.assertArrayEquals("Alpha", alpha2, p.getAlphaLinear(), 1e-10);

					// Solve
					solver.solve(p.getAlphaMatrix(), beta);
					Assert.assertArrayEquals("X", x2, beta, 1e-10);

					for (int j = 0; j < nparams; j++)
					{
						rel[j].add(DoubleEquality.relativeError(x2[j], beta[j]));
						abs[j].add(Math.abs(x2[j] - beta[j]));
					}
				}

				for (int i = 0; i < nparams; i++)
					System.out.printf("Bias = %.2f : %s : Rel %g +/- %g: Abs %g +/- %g\n", b, func.getName(i),
							rel[i].getMean(), rel[i].getStandardDeviation(), abs[i].getMean(),
							abs[i].getStandardDeviation());
			}
		}
		finally
		{
			Background = defaultBackground;
		}
	}

	private double[] add(double[] d, double b)
	{
		d = d.clone();
		for (int i = 0; i < d.length; i++)
			d[i] += b;
		return d;
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
