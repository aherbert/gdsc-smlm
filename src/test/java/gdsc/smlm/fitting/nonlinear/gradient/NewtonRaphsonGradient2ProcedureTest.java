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
import gdsc.smlm.function.Gradient2Function;
import gdsc.smlm.function.gaussian.erf.ErfGaussian2DFunction;
import gdsc.smlm.function.gaussian.erf.SingleFreeCircularErfGaussian2DFunction;

/**
 * Contains speed tests for the methods for calculating the Hessian and gradient vector
 * for use in the LVM algorithm.
 */
public class NewtonRaphsonGradient2ProcedureTest
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
	public void gradientProcedureFactoryCreatesOptimisedProcedures()
	{
		double[] y = new double[0];
		Assert.assertEquals(NewtonRaphsonGradient2ProcedureFactory.create(y, new DummyGradientFunction(4)).getClass(),
				NewtonRaphsonGradient2Procedure4.class);
		Assert.assertEquals(NewtonRaphsonGradient2ProcedureFactory.create(y, new DummyGradientFunction(5)).getClass(),
				NewtonRaphsonGradient2Procedure5.class);
		Assert.assertEquals(NewtonRaphsonGradient2ProcedureFactory.create(y, new DummyGradientFunction(6)).getClass(),
				NewtonRaphsonGradient2Procedure6.class);
	}

	@Test
	public void gradientProcedureComputesSameLogLikelihoodAsMLEGradientCalculator()
	{
		gradientProcedureComputesSameLogLikelihoodAsMLEGradientCalculator(4);
		gradientProcedureComputesSameLogLikelihoodAsMLEGradientCalculator(5);
		gradientProcedureComputesSameLogLikelihoodAsMLEGradientCalculator(6);
		gradientProcedureComputesSameLogLikelihoodAsMLEGradientCalculator(11);
		gradientProcedureComputesSameLogLikelihoodAsMLEGradientCalculator(21);
	}

	@Test
	public void gradientProcedureIsNotSlowerThanGradientCalculator()
	{
		// Note: The procedure does not have a lot of work within loops. It is only a single loop
		// so unrolling does not produce performance gains. The JVM can optimise this. 
		
		gradientProcedureIsNotSlowerThanGradientCalculator(4);
		gradientProcedureIsNotSlowerThanGradientCalculator(5);
		gradientProcedureIsNotSlowerThanGradientCalculator(6);
		gradientProcedureIsNotSlowerThanGradientCalculator(11);
		gradientProcedureIsNotSlowerThanGradientCalculator(21);
	}

	private void gradientProcedureComputesSameLogLikelihoodAsMLEGradientCalculator(int nparams)
	{
		int iter = 10;
		rdg = new RandomDataGenerator(new Well19937c(30051977));

		ArrayList<double[]> paramsList = new ArrayList<double[]>(iter);
		ArrayList<double[]> yList = new ArrayList<double[]>(iter);

		createFakeData(nparams, iter, paramsList, yList);
		FakeGradientFunction func = new FakeGradientFunction(blockWidth, nparams);

		MLEGradientCalculator calc = (MLEGradientCalculator) GradientCalculatorFactory.newCalculator(nparams, true);

		String name = String.format("[%d]", nparams);

		for (int i = 0; i < paramsList.size(); i++)
		{
			NewtonRaphsonGradient2Procedure p = NewtonRaphsonGradient2ProcedureFactory.create(yList.get(i), func);
			double s = p.computeLogLikelihood(paramsList.get(i));
			double s2 = calc.logLikelihood(yList.get(i), paramsList.get(i), func);
			// Exactly the same ...
			Assert.assertEquals(name + " Result: Not same @ " + i, s, s2, 0);
		}
	}

	private abstract class Timer
	{
		private int loops;
		int min = 5;

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

	private void gradientProcedureIsNotSlowerThanGradientCalculator(final int nparams)
	{
		org.junit.Assume.assumeTrue(speedTests || TestSettings.RUN_SPEED_TESTS);

		final int iter = 1000;
		rdg = new RandomDataGenerator(new Well19937c(30051977));

		final ArrayList<double[]> paramsList = new ArrayList<double[]>(iter);
		final ArrayList<double[]> yList = new ArrayList<double[]>(iter);

		createFakeData(nparams, iter, paramsList, yList);
		final FakeGradientFunction func = new FakeGradientFunction(blockWidth, nparams);

		MLEGradientCalculator calc = (MLEGradientCalculator) GradientCalculatorFactory.newCalculator(nparams, true);

		for (int i = 0; i < paramsList.size(); i++)
			calc.logLikelihood(yList.get(i), paramsList.get(i), func);

		for (int i = 0; i < paramsList.size(); i++)
		{
			NewtonRaphsonGradient2Procedure p = NewtonRaphsonGradient2ProcedureFactory.create(yList.get(i), func);
			p.computeLogLikelihood(paramsList.get(i));
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
					MLEGradientCalculator calc = (MLEGradientCalculator) GradientCalculatorFactory
							.newCalculator(nparams, true);
					for (int j = loops; j-- > 0;)
						calc.logLikelihood(yList.get(i), paramsList.get(k++ % iter), func);
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
					NewtonRaphsonGradient2Procedure p = NewtonRaphsonGradient2ProcedureFactory.create(yList.get(i),
							func);
					for (int j = loops; j-- > 0;)
						p.computeLogLikelihood(paramsList.get(k++ % iter));
				}
			}
		};
		long time2 = t2.getTime();

		log("GradientCalculator = %d : NewtonRaphsonGradient2Procedure %d = %d : %fx\n", time1, nparams, time2,
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

		String name = String.format("[%d]", nparams);
		for (int i = 0; i < paramsList.size(); i++)
		{
			NewtonRaphsonGradient2Procedure p1 = new NewtonRaphsonGradient2Procedure(yList.get(i), func);
			NewtonRaphsonGradient2Procedure p2 = NewtonRaphsonGradient2ProcedureFactory.create(yList.get(i), func);
			double[] a = paramsList.get(i);

			double ll1 = p1.computeLogLikelihood(a);
			double ll2 = p2.computeLogLikelihood(a);
			Assert.assertEquals(name + " LL: Not same @ " + i, ll1, ll2, 0);

			p1.computeFirstDerivative(a);
			p2.computeFirstDerivative(a);
			Assert.assertArrayEquals(name + " first derivative: Not same @ " + i, p1.d1, p2.d1, 0);

			p1.computeUpdate(a);
			p2.computeUpdate(a);
			Assert.assertArrayEquals(name + " update: Not same d1 @ " + i, p1.d1, p2.d1, 0);
			Assert.assertArrayEquals(name + " update: Not same d2 @ " + i, p1.d2, p2.d2, 0);
			Assert.assertArrayEquals(name + " update: Not same update @ " + i, p1.getUpdate(), p2.getUpdate(), 0);
		}
	}

	@Test
	public void gradientProcedureIsFasterUnrolledThanGradientProcedure()
	{
		gradientProcedureLinearIsFasterThanGradientProcedure(4);
		gradientProcedureLinearIsFasterThanGradientProcedure(5);
		gradientProcedureLinearIsFasterThanGradientProcedure(6);
	}

	private void gradientProcedureLinearIsFasterThanGradientProcedure(final int nparams)
	{
		org.junit.Assume.assumeTrue(speedTests || TestSettings.RUN_SPEED_TESTS);

		final int iter = 100;
		rdg = new RandomDataGenerator(new Well19937c(30051977));

		final ArrayList<double[]> paramsList = new ArrayList<double[]>(iter);
		final ArrayList<double[]> yList = new ArrayList<double[]>(iter);

		createData(1, iter, paramsList, yList);

		// Remove the timing of the function call by creating a dummy function
		final Gradient2Function func = new FakeGradientFunction(blockWidth, nparams);

		for (int i = 0; i < paramsList.size(); i++)
		{
			NewtonRaphsonGradient2Procedure p1 = new NewtonRaphsonGradient2Procedure(yList.get(i), func);
			p1.computeUpdate(paramsList.get(i));
			p1.computeUpdate(paramsList.get(i));

			NewtonRaphsonGradient2Procedure p2 = NewtonRaphsonGradient2ProcedureFactory.create(yList.get(i), func);
			p2.computeUpdate(paramsList.get(i));
			p2.computeUpdate(paramsList.get(i));

			// Check they are the same
			Assert.assertArrayEquals("Update " + i, p1.getUpdate(), p2.getUpdate(), 0);
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
					NewtonRaphsonGradient2Procedure p1 = new NewtonRaphsonGradient2Procedure(yList.get(i), func);
					for (int j = loops; j-- > 0;)
						p1.computeUpdate(paramsList.get(k++ % iter));
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
					NewtonRaphsonGradient2Procedure p2 = NewtonRaphsonGradient2ProcedureFactory.create(yList.get(i),
							func);
					for (int j = loops; j-- > 0;)
						p2.computeUpdate(paramsList.get(k++ % iter));
				}
			}
		};
		long time2 = t2.getTime();

		log("Standard = %d : Unrolled %d = %d : %fx\n", time1, nparams, time2, (1.0 * time1) / time2);
		Assert.assertTrue(time2 < time1 * 1.5);
	}

	// TODO - Check the first and second derivatives

	//	@Test
	//	public void gradientCalculatorLSQComputesGradient()
	//	{
	//		gradientCalculatorComputesGradient(new SingleFreeCircularErfGaussian2DFunction(blockWidth, blockWidth), false);
	//	}
	//
	//	@Test
	//	public void gradientCalculatorMLEComputesGradient()
	//	{
	//		gradientCalculatorComputesGradient(new SingleFreeCircularErfGaussian2DFunction(blockWidth, blockWidth), true);
	//	}
	//
	//	private void gradientCalculatorComputesGradient(ErfGaussian2DFunction func, boolean mle)
	//	{
	//		int nparams = func.getNumberOfGradients();
	//		int[] indices = func.gradientIndices();
	//
	//		int iter = 100;
	//		rdg = new RandomDataGenerator(new Well19937c(30051977));
	//
	//		ArrayList<double[]> paramsList = new ArrayList<double[]>(iter);
	//		ArrayList<double[]> yList = new ArrayList<double[]>(iter);
	//
	//		createData(1, iter, paramsList, yList, true);
	//
	//		double delta = 1e-3;
	//		DoubleEquality eq = new DoubleEquality(3, 1e-3);
	//
	//		for (int i = 0; i < paramsList.size(); i++)
	//		{
	//			double[] y = yList.get(i);
	//			double[] a = paramsList.get(i);
	//			double[] a2 = a.clone();
	//			NewtonRaphsonGradient2Procedure p = NewtonRaphsonGradient2ProcedureFactory.create(y, func);
	//			p.gradient(a);
	//			//double s = p.value;
	//			double[] beta = p.beta.clone();
	//			for (int j = 0; j < nparams; j++)
	//			{
	//				int k = indices[j];
	//				double d = (a[k] == 0) ? 1e-3 : a[k] * delta;
	//				a2[k] = a[k] + d;
	//				p.value(a2);
	//				double s1 = p.value;
	//				a2[k] = a[k] - d;
	//				p.value(a2);
	//				double s2 = p.value;
	//				a2[k] = a[k];
	//
	//				// Apply a factor of -2 to compute the actual gradients:
	//				// See Numerical Recipes in C++, 2nd Ed. Equation 15.5.6 for Nonlinear Models
	//				beta[j] *= -2;
	//
	//				double gradient = (s1 - s2) / (2 * d);
	//				//System.out.printf("[%d,%d] %f  (%s %f+/-%f)  %f  ?=  %f\n", i, k, s, func.getName(k), a[k], d, beta[j],
	//				//		gradient);
	//				Assert.assertTrue("Not same gradient @ " + j, eq.almostEqualComplement(beta[j], gradient));
	//			}
	//		}
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
