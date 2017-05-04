package gdsc.smlm.fitting.nonlinear.gradient;

import java.util.ArrayList;

import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.ejml.data.DenseMatrix64F;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.utils.DoubleEquality;
import gdsc.core.utils.PseudoRandomGenerator;
import gdsc.smlm.TestSettings;
import gdsc.smlm.function.Gradient1Procedure;
import gdsc.smlm.function.GradientFunction;
import gdsc.smlm.function.ValueProcedure;
import gdsc.smlm.function.gaussian.SingleEllipticalGaussian2DFunction;
import gdsc.smlm.function.gaussian.erf.ErfGaussian2DFunction;
import gdsc.smlm.function.gaussian.erf.SingleCircularErfGaussian2DFunction;
import gdsc.smlm.function.gaussian.erf.SingleFixedErfGaussian2DFunction;
import gdsc.smlm.function.gaussian.erf.SingleFreeCircularErfGaussian2DFunction;

/**
 * Contains speed tests for the fastest method for calculating the Hessian and gradient vector
 * for use in NonLinearFit
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

	private class DummyGradientFunction implements GradientFunction
	{
		int n;

		DummyGradientFunction(int n)
		{
			this.n = n;
		}

		public int size()
		{
			return 0;
		}

		public void initialise(double[] a)
		{
		}

		public int[] gradientIndices()
		{
			return null;
		}

		public int getNumberOfGradients()
		{
			return n;
		}

		public void forEach(ValueProcedure procedure)
		{
		}

		public void forEach(Gradient1Procedure procedure)
		{
		}
	}

	@Test
	public void gradientCalculatorFactoryCreatesOptimisedCalculators()
	{
		double[] y = new double[0];
		Assert.assertEquals(LSQGradientProcedureMatrixFactory.create(y, new DummyGradientFunction(6)).getClass(),
				LSQGradientProcedureMatrix6.class);
	}

	@Test
	public void gradientProcedureComputesSameAsGradientCalculator()
	{
		gradientProcedureComputesSameAsGradientCalculator(
				new SingleFreeCircularErfGaussian2DFunction(blockWidth, blockWidth, 1), 6);
	}

	@Test
	public void gradientProcedureIsFasterThanGradientCalculator()
	{
		gradientProcedureIsFasterThanGradientCalculator(
				new SingleFreeCircularErfGaussian2DFunction(blockWidth, blockWidth, 1), 6);
	}

	private void gradientProcedureComputesSameAsGradientCalculator(ErfGaussian2DFunction func, int nparams)
	{
		// Check the function is the correct size
		Assert.assertEquals(nparams, func.gradientIndices().length);

		int iter = 100;
		rdg = new RandomDataGenerator(new Well19937c(30051977));

		double[][] alpha = new double[nparams][nparams];
		double[] beta = new double[nparams];

		ArrayList<double[]> paramsList = new ArrayList<double[]>(iter);
		ArrayList<double[]> yList = new ArrayList<double[]>(iter);

		int[] x = createData(1, iter, paramsList, yList);

		GradientCalculator calc = GradientCalculatorFactory.newCalculator(nparams, false);

		for (int i = 0; i < paramsList.size(); i++)
		{
			LSQGradientProcedureMatrix p = LSQGradientProcedureMatrixFactory.create(yList.get(i), func);
			p.run(paramsList.get(i));
			double s = p.ssx;
			double s2 = calc.findLinearised(x, yList.get(i), paramsList.get(i), alpha, beta, func);
			// Exactly the same ...
			Assert.assertEquals("Result: Not same @ " + i, s, s2, 0);
			Assert.assertArrayEquals("Observations: Not same beta @ " + i, p.beta, beta, 0);
			for (int j = 0; j < nparams; j++)
				Assert.assertArrayEquals("Observations: Not same alpha @ " + i, p.alpha[j], alpha[j], 0);
		}
	}

	private void gradientProcedureIsFasterThanGradientCalculator(ErfGaussian2DFunction func, int nparams)
	{
		org.junit.Assume.assumeTrue(speedTests || TestSettings.RUN_SPEED_TESTS);

		// Check the function is the correct size
		Assert.assertEquals(nparams, func.gradientIndices().length);

		int iter = 10000;
		rdg = new RandomDataGenerator(new Well19937c(30051977));
		double[][] alpha = new double[nparams][nparams];
		double[] beta = new double[nparams];

		ArrayList<double[]> paramsList = new ArrayList<double[]>(iter);
		ArrayList<double[]> yList = new ArrayList<double[]>(iter);

		int[] x = createData(1, iter, paramsList, yList);

		GradientCalculator calc = GradientCalculatorFactory.newCalculator(nparams, false);

		for (int i = 0; i < paramsList.size(); i++)
			calc.findLinearised(x, yList.get(i), paramsList.get(i), alpha, beta, func);

		for (int i = 0; i < paramsList.size(); i++)
		{
			LSQGradientProcedureMatrix p = LSQGradientProcedureMatrixFactory.create(yList.get(i), func);
			p.run(paramsList.get(i));
		}

		long start1 = System.nanoTime();
		for (int i = 0; i < paramsList.size(); i++)
			calc.findLinearised(x, yList.get(i), paramsList.get(i), alpha, beta, func);
		start1 = System.nanoTime() - start1;

		long start2 = System.nanoTime();
		for (int i = 0; i < paramsList.size(); i++)
		{
			LSQGradientProcedureMatrix p = LSQGradientProcedureMatrixFactory.create(yList.get(i), func);
			p.run(paramsList.get(i));
		}
		start2 = System.nanoTime() - start2;

		log("Linearised GradientCalculator = %d : GradientProcedure%d = %d : %fx\n", start1, nparams, start2,
				(1.0 * start1) / start2);
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(start2 < start1);
	}

	@Test
	public void gradientProcedureLinear4IsFasterThanGradientProcedureMatrix()
	{
		gradientProcedureLinearIsFasterThanGradientProcedureMatrix(
				new SingleFixedErfGaussian2DFunction(blockWidth, blockWidth, 1), 4);
	}

	@Test
	public void gradientProcedureLinear5IsFasterThanGradientProcedureMatrix()
	{
		gradientProcedureLinearIsFasterThanGradientProcedureMatrix(
				new SingleCircularErfGaussian2DFunction(blockWidth, blockWidth, 1), 5);
	}
	
	@Test
	public void gradientProcedureLinear6IsFasterThanGradientProcedureMatrix()
	{
		gradientProcedureLinearIsFasterThanGradientProcedureMatrix(
				new SingleFreeCircularErfGaussian2DFunction(blockWidth, blockWidth, 1), 6);
	}
	
	@Test
	public void gradientProcedureLinear7IsFasterThanGradientProcedureMatrix()
	{
		gradientProcedureLinearIsFasterThanGradientProcedureMatrix(
				new SingleEllipticalGaussian2DFunction(blockWidth, blockWidth), 7);
	}

	private void gradientProcedureLinearIsFasterThanGradientProcedureMatrix(GradientFunction func, final int nparams)
	{
		org.junit.Assume.assumeTrue(speedTests || TestSettings.RUN_SPEED_TESTS);

		// Check the function is the correct size
		Assert.assertEquals(nparams, func.gradientIndices().length);

		// Remove the timing of the function call by creating a dummy function
		final PseudoRandomGenerator r = new PseudoRandomGenerator(10000, new Well19937c(30051977));
		final double[] dy_da = new double[nparams];
		GradientFunction gf = new GradientFunction()
		{
			public int size()
			{
				return 0;
			}

			public void initialise(double[] a)
			{
				int seed = 0;
				for (int i = a.length; i-- > 0;)
					seed += Double.hashCode(a[i]);
				r.setSeed(seed);
			}

			public int[] gradientIndices()
			{
				return null;
			}

			public int getNumberOfGradients()
			{
				return nparams;
			}

			public void forEach(ValueProcedure procedure)
			{

			}

			public void forEach(Gradient1Procedure procedure)
			{
				for (int i = nparams; i-- > 0;)
					dy_da[i] = r.nextDouble();
				procedure.execute(r.nextDouble(), dy_da);
				//procedure.execute(0, dy_da);
			}
		};

		//func = gf;

		int iter = 10000;
		rdg = new RandomDataGenerator(new Well19937c(30051977));

		ArrayList<double[]> paramsList = new ArrayList<double[]>(iter);
		ArrayList<double[]> yList = new ArrayList<double[]>(iter);

		createData(1, iter, paramsList, yList);

		for (int i = 0; i < paramsList.size(); i++)
		{
			LSQGradientProcedureMatrix p = LSQGradientProcedureMatrixFactory.create(yList.get(i), func);
			p.run(paramsList.get(i));
			p.run(paramsList.get(i));

			LSQGradientProcedureLinear p2 = LSQGradientProcedureLinearFactory.create(yList.get(i), func);
			p2.run(paramsList.get(i));
			p2.run(paramsList.get(i));

			// Check they are the same
			DenseMatrix64F m1 = new DenseMatrix64F(p.alpha);
			Assert.assertArrayEquals("A " + i, m1.data, p2.alpha, 0);
			Assert.assertArrayEquals("B " + i, p.beta, p2.beta, 0);
		}

		long start1 = System.nanoTime();
		for (int i = 0; i < paramsList.size(); i++)
		{
			LSQGradientProcedureMatrix p = LSQGradientProcedureMatrixFactory.create(yList.get(i), func);
			p.run(paramsList.get(i));
			p.run(paramsList.get(i));
			p.run(paramsList.get(i));
		}
		start1 = System.nanoTime() - start1;

		long start2 = System.nanoTime();
		for (int i = 0; i < paramsList.size(); i++)
		{
			LSQGradientProcedureLinear p2 = LSQGradientProcedureLinearFactory.create(yList.get(i), func);
			p2.run(paramsList.get(i));
			p2.run(paramsList.get(i));
			p2.run(paramsList.get(i));
		}
		start2 = System.nanoTime() - start2;

		func = gf;

		long start3 = System.nanoTime();
		for (int i = 0; i < paramsList.size(); i++)
		{
			LSQGradientProcedureMatrix p = LSQGradientProcedureMatrixFactory.create(yList.get(i), func);
			p.run(paramsList.get(i));
			p.run(paramsList.get(i));
			p.run(paramsList.get(i));
		}
		start3 = System.nanoTime() - start3;

		long start4 = System.nanoTime();
		for (int i = 0; i < paramsList.size(); i++)
		{
			LSQGradientProcedureLinear p2 = LSQGradientProcedureLinearFactory.create(yList.get(i), func);
			p2.run(paramsList.get(i));
			p2.run(paramsList.get(i));
			p2.run(paramsList.get(i));
		}
		start4 = System.nanoTime() - start4;

		log("GradientProcedure = %d : GradientProcedureLinear %d = %d : %fx\n", start1, nparams, start2,
				(1.0 * start1) / start2);
		log("GradientProcedure = %d : GradientProcedureLinear %d dummy = %d : %fx\n", start3, nparams, start4,
				(1.0 * start3) / start4);
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(start2 < start1);
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(start4 < start3);
	}

	//	@Test
	//	public void gradientCalculatorComputesGradient()
	//	{
	//		gradientCalculatorComputesGradient(new GradientCalculator(7));
	//	}
	//
	//	@Test
	//	public void mleGradientCalculatorComputesGradient()
	//	{
	//		gradientCalculatorComputesGradient(new MLEGradientCalculator(7));
	//	}
	//
	//	private void gradientCalculatorComputesGradient(GradientCalculator calc)
	//	{
	//		int nparams = calc.nparams;
	//		Gaussian2DFunction func = new SingleEllipticalGaussian2DFunction(blockWidth, blockWidth);
	//		// Check the function is the correct size
	//		Assert.assertEquals(nparams, func.gradientIndices().length);
	//
	//		int iter = 100;
	//		rdg = new RandomDataGenerator(new Well19937c(30051977));
	//
	//		double[] beta = new double[nparams];
	//		double[] beta2 = new double[nparams];
	//
	//		ArrayList<double[]> paramsList = new ArrayList<double[]>(iter);
	//		ArrayList<double[]> yList = new ArrayList<double[]>(iter);
	//
	//		int[] x = createData(1, iter, paramsList, yList, true);
	//
	//		double delta = 1e-3;
	//		DoubleEquality eq = new DoubleEquality(3, 1e-3);
	//
	//		for (int i = 0; i < paramsList.size(); i++)
	//		{
	//			double[] y = yList.get(i);
	//			double[] a = paramsList.get(i);
	//			double[] a2 = a.clone();
	//			//double s = 
	//			calc.evaluate(x, y, a, beta, func);
	//
	//			for (int j = 0; j < nparams; j++)
	//			{
	//				double d = (a[j] == 0) ? 1e-3 : a[j] * delta;
	//				a2[j] = a[j] + d;
	//				double s1 = calc.evaluate(x, y, a2, beta2, func);
	//				a2[j] = a[j] - d;
	//				double s2 = calc.evaluate(x, y, a2, beta2, func);
	//				a2[j] = a[j];
	//
	//				double gradient = (s1 - s2) / (2 * d);
	//				//System.out.printf("[%d,%d] %f  (%f+/-%f)  %f  ?=  %f\n", i, j, s, a[j], d, beta[j], gradient);
	//				Assert.assertTrue("Not same gradient @ " + j, eq.almostEqualComplement(beta[j], gradient));
	//			}
	//		}
	//	}
	//
	//	@Test
	//	public void mleGradientCalculatorComputesLikelihood()
	//	{
//		//@formatter:off
//		NonLinearFunction func = new NonLinearFunction(){
//			double u;
//			public void initialise(double[] a) { u = a[0]; }
//			public int[] gradientIndices() { return null; }
//			public double eval(int x, double[] dyda)  { return 0; }
//			public double eval(int x) {
//				return u;
//			}
//			public double eval(int x, double[] dyda, double[] w) { return 0; }
//			public double evalw(int x, double[] w) { return 0; }
//			public boolean canComputeWeights() { return false; }};
//		//@formatter:on
	//
	//		int[] xx = Utils.newArray(100, 0, 1);
	//		double[] xxx = Utils.newArray(100, 0, 1.0);
	//		for (double u : new double[] { 0.79, 2.5, 5.32 })
	//		{
	//			double ll = 0;
	//			PoissonDistribution pd = new PoissonDistribution(u);
	//			for (int x : xx)
	//			{
	//				double o = MLEGradientCalculator.likelihood(u, x);
	//				double e = pd.probability(x);
	//				Assert.assertEquals("likelihood", e, o, e * 1e-10);
	//
	//				o = MLEGradientCalculator.logLikelihood(u, x);
	//				e = pd.logProbability(x);
	//				Assert.assertEquals("log likelihood", e, o, Math.abs(e) * 1e-10);
	//
	//				ll += e;
	//			}
	//
	//			MLEGradientCalculator gc = new MLEGradientCalculator(1);
	//			double o = gc.logLikelihood(xxx, new double[] { u }, func);
	//
	//			Assert.assertEquals("sum log likelihood", ll, o, Math.abs(ll) * 1e-10);
	//		}
	//	}
	//
	//	@Test
	//	public void gradientCalculatorComputesSameOutputWithBias()
	//	{
	//		Gaussian2DFunction func = new SingleEllipticalGaussian2DFunction(blockWidth, blockWidth);
	//		int nparams = func.getNumberOfGradients();
	//		GradientCalculator calc = new GradientCalculator(nparams);
	//		int n = func.size();
	//
	//		int iter = 100;
	//		rdg = new RandomDataGenerator(new Well19937c(30051977));
	//
	//		ArrayList<double[]> paramsList = new ArrayList<double[]>(iter);
	//		ArrayList<double[]> yList = new ArrayList<double[]>(iter);
	//
	//		ArrayList<double[][]> alphaList = new ArrayList<double[][]>(iter);
	//		ArrayList<double[]> betaList = new ArrayList<double[]>(iter);
	//		ArrayList<double[]> xList = new ArrayList<double[]>(iter);
	//
	//		// Manipulate the background
	//		double defaultBackground = Background;
	//		try
	//		{
	//			Background = 1e-2;
	//			createData(1, iter, paramsList, yList, true);
	//
	//			EJMLLinearSolver solver = new EJMLLinearSolver(5, 1e-6);
	//
	//			for (int i = 0; i < paramsList.size(); i++)
	//			{
	//				double[] y = yList.get(i);
	//				double[] a = paramsList.get(i);
	//				double[][] alpha = new double[nparams][nparams];
	//				double[] beta = new double[nparams];
	//				calc.findLinearised(n, y, a, alpha, beta, func);
	//				alphaList.add(alpha);
	//				betaList.add(beta.clone());
	//				for (int j = 0; j < nparams; j++)
	//				{
	//					if (Math.abs(beta[j]) < 1e-6)
	//						System.out.printf("[%d] Tiny beta %s %g\n", i, func.getName(j), beta[j]);
	//				}
	//				// Solve
	//				if (!solver.solve(alpha, beta))
	//					throw new AssertionError();
	//				xList.add(beta);
	//				//System.out.println(Arrays.toString(beta));
	//			}
	//
	//			double[][] alpha = new double[nparams][nparams];
	//			double[] beta = new double[nparams];
	//
	//			//for (int b = 1; b < 1000; b *= 2)
	//			for (double b : new double[] { -500, -100, -10, -1, -0.1, 0, 0.1, 1, 10, 100, 500 })
	//			{
	//				Statistics[] rel = new Statistics[nparams];
	//				Statistics[] abs = new Statistics[nparams];
	//				for (int i = 0; i < nparams; i++)
	//				{
	//					rel[i] = new Statistics();
	//					abs[i] = new Statistics();
	//				}
	//
	//				for (int i = 0; i < paramsList.size(); i++)
	//				{
	//					double[] y = add(yList.get(i), b);
	//					double[] a = paramsList.get(i).clone();
	//					a[0] += b;
	//					calc.findLinearised(n, y, a, alpha, beta, func);
	//					double[][] alpha2 = alphaList.get(i);
	//					double[] beta2 = betaList.get(i);
	//					double[] x2 = xList.get(i);
	//
	//					Assert.assertArrayEquals("Beta", beta2, beta, 1e-10);
	//					for (int j = 0; j < nparams; j++)
	//					{
	//						Assert.assertArrayEquals("Alpha", alpha2[j], alpha[j], 1e-10);
	//					}
	//
	//					// Solve
	//					solver.solve(alpha, beta);
	//					Assert.assertArrayEquals("X", x2, beta, 1e-10);
	//
	//					for (int j = 0; j < nparams; j++)
	//					{
	//						rel[j].add(DoubleEquality.relativeError(x2[j], beta[j]));
	//						abs[j].add(Math.abs(x2[j] - beta[j]));
	//					}
	//				}
	//
	//				for (int i = 0; i < nparams; i++)
	//					System.out.printf("Bias = %.2f : %s : Rel %g +/- %g: Abs %g +/- %g\n", b, func.getName(i),
	//							rel[i].getMean(), rel[i].getStandardDeviation(), abs[i].getMean(),
	//							abs[i].getStandardDeviation());
	//			}
	//		}
	//		finally
	//		{
	//			Background = defaultBackground;
	//		}
	//	}
	//
	//	private double[] add(double[] d, double b)
	//	{
	//		d = d.clone();
	//		for (int i = 0; i < d.length; i++)
	//			d[i] += b;
	//		return d;
	//	}

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
				blockWidth, 1);
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
