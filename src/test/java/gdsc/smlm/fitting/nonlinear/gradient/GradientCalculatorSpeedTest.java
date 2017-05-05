package gdsc.smlm.fitting.nonlinear.gradient;

import java.util.ArrayList;

import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.ij.Utils;
import gdsc.core.utils.DoubleEquality;
import gdsc.core.utils.Statistics;
import gdsc.smlm.TestSettings;
import gdsc.smlm.fitting.linear.EJMLLinearSolver;
import gdsc.smlm.function.CameraNoiseModel;
import gdsc.smlm.function.NonLinearFunction;
import gdsc.smlm.function.gaussian.EllipticalGaussian2DFunction;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.function.gaussian.SingleCircularGaussian2DFunction;
import gdsc.smlm.function.gaussian.SingleEllipticalGaussian2DFunction;
import gdsc.smlm.function.gaussian.SingleFixedGaussian2DFunction;
import gdsc.smlm.function.gaussian.SingleFreeCircularGaussian2DFunction;
import gdsc.smlm.function.gaussian.SingleNBFixedGaussian2DFunction;

/**
 * Contains speed tests for the fastest method for calculating the Hessian and gradient vector
 * for use in NonLinearFit
 */
public class GradientCalculatorSpeedTest
{
	boolean speedTests = false;
	DoubleEquality eq = new DoubleEquality(6, 1e-16);

	int MAX_ITER = 20000;
	int blockWidth = 10;
	double Background = 0.5;
	double Amplitude = 100;
	double Angle = Math.PI;
	double Xpos = 5;
	double Ypos = 5;
	double Xwidth = 1.2;
	double Ywidth = 1.2;

	RandomDataGenerator rdg;

	@Test
	public void gradientCalculatorFactoryCreatesOptimisedCalculators()
	{
		Assert.assertEquals(GradientCalculatorFactory.newCalculator(3).getClass(), GradientCalculator3.class);
		Assert.assertEquals(GradientCalculatorFactory.newCalculator(4).getClass(), GradientCalculator4.class);
		Assert.assertEquals(GradientCalculatorFactory.newCalculator(5).getClass(), GradientCalculator5.class);
		Assert.assertEquals(GradientCalculatorFactory.newCalculator(6).getClass(), GradientCalculator6.class);
		Assert.assertEquals(GradientCalculatorFactory.newCalculator(7).getClass(), GradientCalculator7.class);
		Assert.assertEquals(GradientCalculatorFactory.newCalculator(13).getClass(), GradientCalculator.class);

		Assert.assertEquals(GradientCalculatorFactory.newCalculator(3, true).getClass(), MLEGradientCalculator3.class);
		Assert.assertEquals(GradientCalculatorFactory.newCalculator(4, true).getClass(), MLEGradientCalculator4.class);
		Assert.assertEquals(GradientCalculatorFactory.newCalculator(5, true).getClass(), MLEGradientCalculator5.class);
		Assert.assertEquals(GradientCalculatorFactory.newCalculator(6, true).getClass(), MLEGradientCalculator6.class);
		Assert.assertEquals(GradientCalculatorFactory.newCalculator(7, true).getClass(), MLEGradientCalculator7.class);
		Assert.assertEquals(GradientCalculatorFactory.newCalculator(13, true).getClass(), MLEGradientCalculator.class);
	}

	@Test
	public void mleGradientCalculator7ComputesSameAsGradientCalculator()
	{
		gradientCalculatorNComputesSameAsGradientCalculator(
				new SingleEllipticalGaussian2DFunction(blockWidth, blockWidth), 7, true);
	}

	@Test
	public void mleGradientCalculator7IsFasterThanGradientCalculator()
	{
		gradientCalculatorNIsFasterThanGradientCalculator(
				new SingleEllipticalGaussian2DFunction(blockWidth, blockWidth), 7, true);
	}

	@Test
	public void mleGradientCalculator6ComputesSameAsGradientCalculator()
	{
		gradientCalculatorNComputesSameAsGradientCalculator(
				new SingleFreeCircularGaussian2DFunction(blockWidth, blockWidth), 6, true);
	}

	@Test
	public void mleGradientCalculator6IsFasterThanGradientCalculator()
	{
		gradientCalculatorNIsFasterThanGradientCalculator(
				new SingleFreeCircularGaussian2DFunction(blockWidth, blockWidth), 6, true);
	}

	@Test
	public void mleGradientCalculator5ComputesSameAsGradientCalculator()
	{
		gradientCalculatorNComputesSameAsGradientCalculator(
				new SingleCircularGaussian2DFunction(blockWidth, blockWidth), 5, true);
	}

	@Test
	public void mleGradientCalculator5IsFasterThanGradientCalculator()
	{
		gradientCalculatorNIsFasterThanGradientCalculator(new SingleCircularGaussian2DFunction(blockWidth, blockWidth),
				5, true);
	}

	@Test
	public void mleGradientCalculator4ComputesSameAsGradientCalculator()
	{
		gradientCalculatorNComputesSameAsGradientCalculator(new SingleFixedGaussian2DFunction(blockWidth, blockWidth),
				4, true);
	}

	@Test
	public void mleGradientCalculator4IsFasterThanGradientCalculator()
	{
		gradientCalculatorNIsFasterThanGradientCalculator(new SingleFixedGaussian2DFunction(blockWidth, blockWidth), 4,
				true);
	}

	@Test
	public void mleGradientCalculator3ComputesSameAsGradientCalculator()
	{
		gradientCalculatorNComputesSameAsGradientCalculator(new SingleNBFixedGaussian2DFunction(blockWidth, blockWidth),
				3, true);
	}

	@Test
	public void mleGradientCalculator3IsFasterThanGradientCalculator()
	{
		gradientCalculatorNIsFasterThanGradientCalculator(new SingleNBFixedGaussian2DFunction(blockWidth, blockWidth),
				3, true);
	}

	@Test
	public void gradientCalculator7ComputesSameAsGradientCalculator()
	{
		gradientCalculatorNComputesSameAsGradientCalculator(
				new SingleEllipticalGaussian2DFunction(blockWidth, blockWidth), 7, false);
	}

	@Test
	public void gradientCalculator7IsFasterThanGradientCalculator()
	{
		gradientCalculatorNIsFasterThanGradientCalculator(
				new SingleEllipticalGaussian2DFunction(blockWidth, blockWidth), 7, false);
	}

	@Test
	public void gradientCalculator6ComputesSameAsGradientCalculator()
	{
		gradientCalculatorNComputesSameAsGradientCalculator(
				new SingleFreeCircularGaussian2DFunction(blockWidth, blockWidth), 6, false);
	}

	@Test
	public void gradientCalculator6IsFasterThanGradientCalculator()
	{
		gradientCalculatorNIsFasterThanGradientCalculator(
				new SingleFreeCircularGaussian2DFunction(blockWidth, blockWidth), 6, false);
	}

	@Test
	public void gradientCalculator5ComputesSameAsGradientCalculator()
	{
		gradientCalculatorNComputesSameAsGradientCalculator(
				new SingleCircularGaussian2DFunction(blockWidth, blockWidth), 5, false);
	}

	@Test
	public void gradientCalculator5IsFasterThanGradientCalculator()
	{
		gradientCalculatorNIsFasterThanGradientCalculator(new SingleCircularGaussian2DFunction(blockWidth, blockWidth),
				5, false);
	}

	@Test
	public void gradientCalculator4ComputesSameAsGradientCalculator()
	{
		gradientCalculatorNComputesSameAsGradientCalculator(new SingleFixedGaussian2DFunction(blockWidth, blockWidth),
				4, false);
	}

	@Test
	public void gradientCalculator4IsFasterThanGradientCalculator()
	{
		gradientCalculatorNIsFasterThanGradientCalculator(new SingleFixedGaussian2DFunction(blockWidth, blockWidth), 4,
				false);
	}

	@Test
	public void gradientCalculator3ComputesSameAsGradientCalculator()
	{
		gradientCalculatorNComputesSameAsGradientCalculator(new SingleNBFixedGaussian2DFunction(blockWidth, blockWidth),
				3, false);
	}

	@Test
	public void gradientCalculator3IsFasterThanGradientCalculator()
	{
		gradientCalculatorNIsFasterThanGradientCalculator(new SingleNBFixedGaussian2DFunction(blockWidth, blockWidth),
				3, false);
	}

	private void gradientCalculatorNComputesSameAsGradientCalculator(Gaussian2DFunction func, int nparams, boolean mle)
	{
		// Check the function is the correct size
		Assert.assertEquals(nparams, func.gradientIndices().length);

		int iter = 100;
		rdg = new RandomDataGenerator(new Well19937c(30051977));

		double[][] alpha = new double[nparams][nparams];
		double[] beta = new double[nparams];
		double[][] alpha2 = new double[nparams][nparams];
		double[] beta2 = new double[nparams];

		ArrayList<double[]> paramsList = new ArrayList<double[]>(iter);
		ArrayList<double[]> yList = new ArrayList<double[]>(iter);

		int[] x = createData(1, iter, paramsList, yList);

		GradientCalculator calc = (mle) ? new MLEGradientCalculator(beta.length) : new GradientCalculator(beta.length);
		GradientCalculator calc2 = GradientCalculatorFactory.newCalculator(nparams, mle);

		for (int i = 0; i < paramsList.size(); i++)
		{
			double s = calc.findLinearised(x, yList.get(i), paramsList.get(i), alpha, beta, func);
			double s2 = calc2.findLinearised(x, yList.get(i), paramsList.get(i), alpha2, beta2, func);
			Assert.assertTrue("Result: Not same @ " + i, eq.almostEqualComplement(s, s2));
			Assert.assertTrue("Observations: Not same beta @ " + i, eq.almostEqualComplement(beta, beta2));
			for (int j = 0; j < beta.length; j++)
				Assert.assertTrue("Observations: Not same alpha @ " + i, eq.almostEqualComplement(alpha[j], alpha2[j]));
		}

		for (int i = 0; i < paramsList.size(); i++)
		{
			double s = calc.findLinearised(x.length, yList.get(i), paramsList.get(i), alpha, beta, func);
			double s2 = calc2.findLinearised(x.length, yList.get(i), paramsList.get(i), alpha2, beta2, func);
			Assert.assertTrue("N-Result: Not same @ " + i, eq.almostEqualComplement(s, s2));
			Assert.assertTrue("N-observations: Not same beta @ " + i, eq.almostEqualComplement(beta, beta2));
			for (int j = 0; j < beta.length; j++)
				Assert.assertTrue("N-observations: Not same alpha @ " + i,
						eq.almostEqualComplement(alpha[j], alpha2[j]));
		}

		if (!mle)
		{
			func.setNoiseModel(CameraNoiseModel.createNoiseModel(10, 0, true));

			for (int i = 0; i < paramsList.size(); i++)
			{
				double s = calc.findLinearised(x, yList.get(i), paramsList.get(i), alpha, beta, func);
				double s2 = calc2.findLinearised(x, yList.get(i), paramsList.get(i), alpha2, beta2, func);
				Assert.assertTrue("Result+Noise: Not same @ " + i, eq.almostEqualComplement(s, s2));
				Assert.assertTrue("Observations+Noise: Not same beta @ " + i, eq.almostEqualComplement(beta, beta2));
				for (int j = 0; j < beta.length; j++)
					Assert.assertTrue("Observations+Noise: Not same alpha @ " + i,
							eq.almostEqualComplement(alpha[j], alpha2[j]));
			}

			for (int i = 0; i < paramsList.size(); i++)
			{
				double s = calc.findLinearised(x.length, yList.get(i), paramsList.get(i), alpha, beta, func);
				double s2 = calc2.findLinearised(x.length, yList.get(i), paramsList.get(i), alpha2, beta2, func);
				Assert.assertTrue("N-Result+Noise: Not same @ " + i, eq.almostEqualComplement(s, s2));
				Assert.assertTrue("N-Observations+Noise: Not same beta @ " + i, eq.almostEqualComplement(beta, beta2));
				for (int j = 0; j < beta.length; j++)
					Assert.assertTrue("N-Observations+Noise: Not same alpha @ " + i,
							eq.almostEqualComplement(alpha[j], alpha2[j]));
			}
		}

		// Only the diagonal Fisher Information has been unrolled into the other calculators
		for (int i = 0; i < paramsList.size(); i++)
		{
			beta = calc.fisherInformationDiagonal(x.length, paramsList.get(i), func);
			beta2 = calc.fisherInformationDiagonal(x.length, paramsList.get(i), func);
			Assert.assertTrue("Not same diagonal @ " + i, eq.almostEqualComplement(beta, beta2));
		}
	}

	private void gradientCalculatorNIsFasterThanGradientCalculator(Gaussian2DFunction func, int nparams, boolean mle)
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

		GradientCalculator calc = (mle) ? new MLEGradientCalculator(beta.length) : new GradientCalculator(beta.length);
		GradientCalculator calc2 = GradientCalculatorFactory.newCalculator(nparams, mle);

		for (int i = 0; i < paramsList.size(); i++)
			calc.findLinearised(x, yList.get(i), paramsList.get(i), alpha, beta, func);

		for (int i = 0; i < paramsList.size(); i++)
			calc2.findLinearised(x, yList.get(i), paramsList.get(i), alpha, beta, func);

		long start1 = System.nanoTime();
		for (int i = 0; i < paramsList.size(); i++)
			calc.findLinearised(x, yList.get(i), paramsList.get(i), alpha, beta, func);
		start1 = System.nanoTime() - start1;

		long start2 = System.nanoTime();
		for (int i = 0; i < paramsList.size(); i++)
			calc2.findLinearised(x, yList.get(i), paramsList.get(i), alpha, beta, func);
		start2 = System.nanoTime() - start2;

		log("%sLinearised GradientCalculator = %d : GradientCalculator%d = %d : %fx\n", (mle) ? "MLE " : "", start1,
				nparams, start2, (1.0 * start1) / start2);
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(start2 < start1);

		for (int i = 0; i < paramsList.size(); i++)
			calc.fisherInformationDiagonal(x.length, paramsList.get(i), func);

		for (int i = 0; i < paramsList.size(); i++)
			calc2.fisherInformationDiagonal(x.length, paramsList.get(i), func);

		start1 = System.nanoTime();
		for (int i = 0; i < paramsList.size(); i++)
			calc.fisherInformationDiagonal(x.length, paramsList.get(i), func);
		start1 = System.nanoTime() - start1;

		start2 = System.nanoTime();
		for (int i = 0; i < paramsList.size(); i++)
			calc2.fisherInformationDiagonal(x.length, paramsList.get(i), func);
		start2 = System.nanoTime() - start2;

		log("%sFisher Diagonal GradientCalculator = %d : GradientCalculator%d = %d : %fx\n", (mle) ? "MLE " : "",
				start1, nparams, start2, (1.0 * start1) / start2);
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(start2 < start1);
	}

	@Test
	public void gradientCalculatorAssumedXIsFasterThanGradientCalculator()
	{
		org.junit.Assume.assumeTrue(speedTests || TestSettings.RUN_SPEED_TESTS);

		int iter = 10000;
		rdg = new RandomDataGenerator(new Well19937c(30051977));
		double[][] alpha = new double[7][7];
		double[] beta = new double[7];

		ArrayList<double[]> paramsList = new ArrayList<double[]>(iter);
		ArrayList<double[]> yList = new ArrayList<double[]>(iter);

		int[] x = createData(1, iter, paramsList, yList);

		GradientCalculator calc = new GradientCalculator6();
		GradientCalculator calc2 = new GradientCalculator6();
		SingleFreeCircularGaussian2DFunction func = new SingleFreeCircularGaussian2DFunction(blockWidth, blockWidth);
		int n = x.length;

		for (int i = 0; i < paramsList.size(); i++)
			calc.findLinearised(x, yList.get(i), paramsList.get(i), alpha, beta, func);

		for (int i = 0; i < paramsList.size(); i++)
			calc2.findLinearised(n, yList.get(i), paramsList.get(i), alpha, beta, func);

		long start1 = System.nanoTime();
		for (int i = 0; i < paramsList.size(); i++)
			calc.findLinearised(x, yList.get(i), paramsList.get(i), alpha, beta, func);
		start1 = System.nanoTime() - start1;

		long start2 = System.nanoTime();
		for (int i = 0; i < paramsList.size(); i++)
			calc2.findLinearised(n, yList.get(i), paramsList.get(i), alpha, beta, func);
		start2 = System.nanoTime() - start2;

		log("GradientCalculator = %d : GradientCalculatorAssumed = %d : %fx\n", start1, start2,
				(1.0 * start1) / start2);
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(start2 < start1);
	}

	@Test
	public void gradientCalculatorComputesGradient()
	{
		gradientCalculatorComputesGradient(new GradientCalculator(7));
	}

	@Test
	public void mleGradientCalculatorComputesGradient()
	{
		gradientCalculatorComputesGradient(new MLEGradientCalculator(7));
	}

	private void gradientCalculatorComputesGradient(GradientCalculator calc)
	{
		int nparams = calc.nparams;
		Gaussian2DFunction func = new SingleEllipticalGaussian2DFunction(blockWidth, blockWidth);
		// Check the function is the correct size
		Assert.assertEquals(nparams, func.gradientIndices().length);

		int iter = 100;
		rdg = new RandomDataGenerator(new Well19937c(30051977));

		double[] beta = new double[nparams];
		double[] beta2 = new double[nparams];

		ArrayList<double[]> paramsList = new ArrayList<double[]>(iter);
		ArrayList<double[]> yList = new ArrayList<double[]>(iter);

		int[] x = createData(1, iter, paramsList, yList, true);

		double delta = 1e-3;
		DoubleEquality eq = new DoubleEquality(3, 1e-3);

		for (int i = 0; i < paramsList.size(); i++)
		{
			double[] y = yList.get(i);
			double[] a = paramsList.get(i);
			double[] a2 = a.clone();
			double s = calc.evaluate(x, y, a, beta, func);

			for (int j = 0; j < nparams; j++)
			{
				double d = (a[j] == 0) ? 1e-3 : a[j] * delta;
				a2[j] = a[j] + d;
				double s1 = calc.evaluate(x, y, a2, beta2, func);
				a2[j] = a[j] - d;
				double s2 = calc.evaluate(x, y, a2, beta2, func);
				a2[j] = a[j];

				double gradient = (s1 - s2) / (2 * d);
				System.out.printf("[%d,%d] %f  (%s %f+/-%f)  %f  ?=  %f\n", i, j, s, func.getName(j), a[j], d, beta[j],
						gradient);
				Assert.assertTrue("Not same gradient @ " + j, eq.almostEqualComplement(beta[j], gradient));
			}
			break;
		}
	}

	@Test
	public void mleGradientCalculatorComputesLikelihood()
	{
		//@formatter:off
		NonLinearFunction func = new NonLinearFunction(){
			double u;
			public void initialise(double[] a) { u = a[0]; }
			public int[] gradientIndices() { return null; }
			public double eval(int x, double[] dyda)  { return 0; }
			public double eval(int x) {
				return u;
			}
			public double eval(int x, double[] dyda, double[] w) { return 0; }
			public double evalw(int x, double[] w) { return 0; }
			public boolean canComputeWeights() { return false; }};
		//@formatter:on

		int[] xx = Utils.newArray(100, 0, 1);
		double[] xxx = Utils.newArray(100, 0, 1.0);
		for (double u : new double[] { 0.79, 2.5, 5.32 })
		{
			double ll = 0;
			PoissonDistribution pd = new PoissonDistribution(u);
			for (int x : xx)
			{
				double o = MLEGradientCalculator.likelihood(u, x);
				double e = pd.probability(x);
				Assert.assertEquals("likelihood", e, o, e * 1e-10);

				o = MLEGradientCalculator.logLikelihood(u, x);
				e = pd.logProbability(x);
				Assert.assertEquals("log likelihood", e, o, Math.abs(e) * 1e-10);

				ll += e;
			}

			MLEGradientCalculator gc = new MLEGradientCalculator(1);
			double o = gc.logLikelihood(xxx, new double[] { u }, func);

			Assert.assertEquals("sum log likelihood", ll, o, Math.abs(ll) * 1e-10);
		}
	}

	@Test
	public void gradientCalculatorComputesSameOutputWithBias()
	{
		Gaussian2DFunction func = new SingleEllipticalGaussian2DFunction(blockWidth, blockWidth);
		int nparams = func.getNumberOfGradients();
		GradientCalculator calc = new GradientCalculator(nparams);
		int n = func.size();

		int iter = 100;
		rdg = new RandomDataGenerator(new Well19937c(30051977));

		ArrayList<double[]> paramsList = new ArrayList<double[]>(iter);
		ArrayList<double[]> yList = new ArrayList<double[]>(iter);

		ArrayList<double[][]> alphaList = new ArrayList<double[][]>(iter);
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
				double[][] alpha = new double[nparams][nparams];
				double[] beta = new double[nparams];
				calc.findLinearised(n, y, a, alpha, beta, func);
				alphaList.add(alpha);
				betaList.add(beta.clone());
				for (int j = 0; j < nparams; j++)
				{
					if (Math.abs(beta[j]) < 1e-6)
						System.out.printf("[%d] Tiny beta %s %g\n", i, func.getName(j), beta[j]);
				}
				// Solve
				if (!solver.solve(alpha, beta))
					throw new AssertionError();
				xList.add(beta);
				//System.out.println(Arrays.toString(beta));
			}

			double[][] alpha = new double[nparams][nparams];
			double[] beta = new double[nparams];

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
					calc.findLinearised(n, y, a, alpha, beta, func);
					double[][] alpha2 = alphaList.get(i);
					double[] beta2 = betaList.get(i);
					double[] x2 = xList.get(i);

					Assert.assertArrayEquals("Beta", beta2, beta, 1e-10);
					for (int j = 0; j < nparams; j++)
					{
						Assert.assertArrayEquals("Alpha", alpha2[j], alpha[j], 1e-10);
					}

					// Solve
					solver.solve(alpha, beta);
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
		EllipticalGaussian2DFunction func = new EllipticalGaussian2DFunction(npeaks, blockWidth, blockWidth);
		params[0] = random(Background);
		for (int i = 0, j = 1; i < npeaks; i++, j += 6)
		{
			params[j] = random(Amplitude);
			params[j + 1] = random(Angle);
			params[j + 2] = random(Xpos);
			params[j + 3] = random(Ypos);
			params[j + 4] = random(Xwidth);
			params[j + 5] = random(Ywidth);
		}

		double[] dy_da = new double[params.length];
		double[] y = new double[n];
		func.initialise(params);
		for (int i = 0; i < y.length; i++)
		{
			// Add random Poisson noise
			y[i] = rdg.nextPoisson(func.eval(i, dy_da));
		}

		if (randomiseParams)
		{
			params[0] = random(params[0]);
			for (int i = 0, j = 1; i < npeaks; i++, j += 6)
			{
				params[j] = random(params[j]);
				params[j + 1] = random(params[j + 1]);
				params[j + 2] = random(params[j + 2]);
				params[j + 3] = random(params[j + 3]);
				params[j + 4] = random(params[j + 4]);
				params[j + 5] = random(params[j + 5]); //params[j + 4];
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
