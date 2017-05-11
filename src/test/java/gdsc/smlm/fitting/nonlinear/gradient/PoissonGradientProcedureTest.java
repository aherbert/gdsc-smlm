package gdsc.smlm.fitting.nonlinear.gradient;

import java.util.ArrayList;

import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.ejml.data.DenseMatrix64F;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.ij.Utils;
import gdsc.core.utils.DoubleEquality;
import gdsc.smlm.TestSettings;
import gdsc.smlm.fitting.FisherInformationMatrix;
import gdsc.smlm.function.Gradient1Function;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import gdsc.smlm.function.gaussian.erf.ErfGaussian2DFunction;

public class PoissonGradientProcedureTest
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
		Assert.assertEquals(PoissonGradientProcedureFactory.create(new DummyGradientFunction(6)).getClass(),
				PoissonGradientProcedure6.class);
		Assert.assertEquals(PoissonGradientProcedureFactory.create(new DummyGradientFunction(5)).getClass(),
				PoissonGradientProcedure5.class);
		Assert.assertEquals(PoissonGradientProcedureFactory.create(new DummyGradientFunction(4)).getClass(),
				PoissonGradientProcedure4.class);
		double[] b = null;
		Assert.assertEquals(PoissonGradientProcedureFactory.create(b, new DummyGradientFunction(6)).getClass(),
				PoissonGradientProcedure6.class);
		Assert.assertEquals(PoissonGradientProcedureFactory.create(b, new DummyGradientFunction(5)).getClass(),
				PoissonGradientProcedure5.class);
		Assert.assertEquals(PoissonGradientProcedureFactory.create(b, new DummyGradientFunction(4)).getClass(),
				PoissonGradientProcedure4.class);
		b = new double[new DummyGradientFunction(4).size()];
		Assert.assertEquals(PoissonGradientProcedureFactory.create(b, new DummyGradientFunction(6)).getClass(),
				PoissonGradientProcedureB6.class);
		Assert.assertEquals(PoissonGradientProcedureFactory.create(b, new DummyGradientFunction(5)).getClass(),
				PoissonGradientProcedureB5.class);
		Assert.assertEquals(PoissonGradientProcedureFactory.create(b, new DummyGradientFunction(4)).getClass(),
				PoissonGradientProcedureB4.class);
		Assert.assertEquals(PoissonGradientProcedureFactory.create(b, new DummyGradientFunction(1)).getClass(),
				PoissonGradientProcedureB.class);
	}

	@Test
	public void gradientProcedureComputesSameAsGradientCalculator()
	{
		gradientProcedureComputesSameAsGradientCalculator(4);
		gradientProcedureComputesSameAsGradientCalculator(5);
		gradientProcedureComputesSameAsGradientCalculator(6);
		gradientProcedureComputesSameAsGradientCalculator(11);
		gradientProcedureComputesSameAsGradientCalculator(21);
	}

	@Test
	public void gradientProcedureIsNotSlowerThanGradientCalculator()
	{
		gradientProcedureIsNotSlowerThanGradientCalculator(4);
		gradientProcedureIsNotSlowerThanGradientCalculator(5);
		gradientProcedureIsNotSlowerThanGradientCalculator(6);
		// 2 peaks
		gradientProcedureIsNotSlowerThanGradientCalculator(11);
		// 4 peaks
		gradientProcedureIsNotSlowerThanGradientCalculator(21);
	}

	private void gradientProcedureComputesSameAsGradientCalculator(int nparams)
	{
		int iter = 10;
		rdg = new RandomDataGenerator(new Well19937c(30051977));

		ArrayList<double[]> paramsList = new ArrayList<double[]>(iter);

		createFakeParams(nparams, iter, paramsList);
		int n = blockWidth * blockWidth;
		FakeGradientFunction func = new FakeGradientFunction(blockWidth, nparams);

		GradientCalculator calc = GradientCalculatorFactory.newCalculator(nparams, false);

		String name = String.format("[%d]", nparams);

		for (int i = 0; i < paramsList.size(); i++)
		{
			PoissonGradientProcedure p = PoissonGradientProcedureFactory.create(func);
			p.computeFisherInformation(paramsList.get(i));
			double[][] m = calc.fisherInformationMatrix(n, paramsList.get(i), func);
			// Exactly the same ...
			double[] al = p.getLinear();
			Assert.assertArrayEquals(name + " Observations: Not same alpha @ " + i, al, new DenseMatrix64F(m).data, 0);

			double[][] am = p.getMatrix();
			for (int j = 0; j < nparams; j++)
				Assert.assertArrayEquals(name + " Observations: Not same alpha @ " + i, am[j], m[j], 0);
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

	private void gradientProcedureIsNotSlowerThanGradientCalculator(final int nparams)
	{
		org.junit.Assume.assumeTrue(speedTests || TestSettings.RUN_SPEED_TESTS);

		final int iter = 1000;
		rdg = new RandomDataGenerator(new Well19937c(30051977));

		final ArrayList<double[]> paramsList = new ArrayList<double[]>(iter);

		createFakeParams(nparams, iter, paramsList);
		final int n = blockWidth * blockWidth;
		final FakeGradientFunction func = new FakeGradientFunction(blockWidth, nparams);

		GradientCalculator calc = GradientCalculatorFactory.newCalculator(nparams, false);

		for (int i = 0; i < paramsList.size(); i++)
			calc.fisherInformationMatrix(n, paramsList.get(i), func);

		for (int i = 0; i < paramsList.size(); i++)
		{
			PoissonGradientProcedure p = PoissonGradientProcedureFactory.create(func);
			p.computeFisherInformation(paramsList.get(i));
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
						calc.fisherInformationMatrix(n, paramsList.get(k++ % iter), func);
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
					PoissonGradientProcedure p = PoissonGradientProcedureFactory.create(func);
					for (int j = loops; j-- > 0;)
						p.computeFisherInformation(paramsList.get(k++ % iter));
				}
			}
		};
		long time2 = t2.getTime();

		log("GradientCalculator = %d : PoissonGradientProcedure %d = %d : %fx\n", time1, nparams, time2,
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
		gradientProcedureUnrolledComputesSameAsGradientProcedure(4, false);
		gradientProcedureUnrolledComputesSameAsGradientProcedure(5, false);
		gradientProcedureUnrolledComputesSameAsGradientProcedure(6, false);
	}

	@Test
	public void gradientProcedureUnrolledComputesSameAsGradientProcedureWithPrecomputed()
	{
		gradientProcedureUnrolledComputesSameAsGradientProcedure(4, true);
		gradientProcedureUnrolledComputesSameAsGradientProcedure(5, true);
		gradientProcedureUnrolledComputesSameAsGradientProcedure(6, true);
	}

	private void gradientProcedureUnrolledComputesSameAsGradientProcedure(int nparams, boolean precomputed)
	{
		int iter = 10;
		rdg = new RandomDataGenerator(new Well19937c(30051977));

		ArrayList<double[]> paramsList = new ArrayList<double[]>(iter);

		createFakeParams(nparams, iter, paramsList);
		FakeGradientFunction func = new FakeGradientFunction(blockWidth, nparams);

		double[] b = (precomputed) ? Utils.newArray(func.size(), 0.1, 1.3) : null;

		String name = String.format("[%d]", nparams);
		for (int i = 0; i < paramsList.size(); i++)
		{
			PoissonGradientProcedure p1 = (precomputed) ? new PoissonGradientProcedureB(b, func)
					: new PoissonGradientProcedure(func);
			p1.computeFisherInformation(paramsList.get(i));

			PoissonGradientProcedure p2 = PoissonGradientProcedureFactory.create(b, func);
			p2.computeFisherInformation(paramsList.get(i));

			// Exactly the same ...
			Assert.assertArrayEquals(name + " Observations: Not same alpha @ " + i, p1.getLinear(), p2.getLinear(), 0);

			double[][] am1 = p1.getMatrix();
			double[][] am2 = p2.getMatrix();
			for (int j = 0; j < nparams; j++)
				Assert.assertArrayEquals(name + " Observations: Not same alpha @ " + i, am1[j], am2[j], 0);
		}
	}

	@Test
	public void gradientProcedureIsFasterUnrolledThanGradientProcedure()
	{
		gradientProcedureIsFasterUnrolledThanGradientProcedure(4, false);
		gradientProcedureIsFasterUnrolledThanGradientProcedure(5, false);
		gradientProcedureIsFasterUnrolledThanGradientProcedure(6, false);
	}

	@Test
	public void gradientProcedureIsFasterUnrolledThanGradientProcedureWithPrecomputed()
	{
		gradientProcedureIsFasterUnrolledThanGradientProcedure(4, true);
		gradientProcedureIsFasterUnrolledThanGradientProcedure(5, true);
		gradientProcedureIsFasterUnrolledThanGradientProcedure(6, true);
	}

	private void gradientProcedureIsFasterUnrolledThanGradientProcedure(final int nparams, final boolean precomputed)
	{
		org.junit.Assume.assumeTrue(speedTests || TestSettings.RUN_SPEED_TESTS);

		final int iter = 100;
		rdg = new RandomDataGenerator(new Well19937c(30051977));

		final ArrayList<double[]> paramsList = new ArrayList<double[]>(iter);

		createFakeParams(nparams, iter, paramsList);

		// Remove the timing of the function call by creating a dummy function
		final Gradient1Function func = new FakeGradientFunction(blockWidth, nparams);
		final double[] b = (precomputed) ? Utils.newArray(func.size(), 0.1, 1.3) : null;

		for (int i = 0; i < paramsList.size(); i++)
		{
			PoissonGradientProcedure p1 = (precomputed) ? new PoissonGradientProcedureB(b, func)
					: new PoissonGradientProcedure(func);
			p1.computeFisherInformation(paramsList.get(i));
			p1.computeFisherInformation(paramsList.get(i));

			PoissonGradientProcedure p2 = PoissonGradientProcedureFactory.create(b, func);
			p2.computeFisherInformation(paramsList.get(i));
			p2.computeFisherInformation(paramsList.get(i));

			// Check they are the same
			Assert.assertArrayEquals("M " + i, p1.getLinear(), p2.getLinear(), 0);
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
					PoissonGradientProcedure p1 = (precomputed) ? new PoissonGradientProcedureB(b, func)
							: new PoissonGradientProcedure(func);
					for (int j = loops; j-- > 0;)
						p1.computeFisherInformation(paramsList.get(k++ % iter));
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
					PoissonGradientProcedure p2 = PoissonGradientProcedureFactory.create(b, func);
					for (int j = loops; j-- > 0;)
						p2.computeFisherInformation(paramsList.get(k++ % iter));
				}
			}
		};
		long time2 = t2.getTime();

		log("Precomputed=%b : Standard %d : Unrolled %d = %d : %fx\n", precomputed, time1, nparams, time2,
				(1.0 * time1) / time2);
		Assert.assertTrue(time2 < time1);
	}

	@Test
	public void crlbIsHigherWithPrecomputed()
	{
		int iter = 10;
		rdg = new RandomDataGenerator(new Well19937c(30051977));

		ErfGaussian2DFunction func = (ErfGaussian2DFunction) GaussianFunctionFactory.create2D(1, 10, 10,
				GaussianFunctionFactory.FIT_ERF_FREE_CIRCLE, null);

		double[] a = new double[7];
		int n = func.getNumberOfGradients();

		// Get a background
		double[] b = new double[func.size()];
		for (int i = 0; i < b.length; i++)
			b[i] = rdg.nextUniform(1, 2);

		for (int i = 0; i < iter; i++)
		{
			a[Gaussian2DFunction.BACKGROUND] = rdg.nextUniform(0.1, 0.3);
			a[Gaussian2DFunction.SIGNAL] = rdg.nextUniform(100, 300);
			a[Gaussian2DFunction.X_POSITION] = rdg.nextUniform(4, 6);
			a[Gaussian2DFunction.Y_POSITION] = rdg.nextUniform(4, 6);
			a[Gaussian2DFunction.X_SD] = rdg.nextUniform(1, 1.3);
			a[Gaussian2DFunction.Y_SD] = rdg.nextUniform(1, 1.3);

			PoissonGradientProcedure p1 = PoissonGradientProcedureFactory.create(func);
			p1.computeFisherInformation(a);

			PoissonGradientProcedure p2 = PoissonGradientProcedureFactory.create(b, func);
			p2.computeFisherInformation(a);

			FisherInformationMatrix m1 = new FisherInformationMatrix(p1.getLinear(), n);
			FisherInformationMatrix m2 = new FisherInformationMatrix(p2.getLinear(), n);
			double[] crlb1 = m1.crlb();
			double[] crlb2 = m2.crlb();
			Assert.assertNotNull(crlb1);
			Assert.assertNotNull(crlb2);
			//System.out.printf("%s : %s\n", Arrays.toString(crlb1), Arrays.toString(crlb2));
			for (int j = 0; j < n; j++)
				Assert.assertTrue(crlb1[j] < crlb2[j]);
		}
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

	protected void createFakeParams(int nparams, int iter, ArrayList<double[]> paramsList)
	{
		for (int i = 0; i < iter; i++)
		{
			double[] params = new double[nparams];
			createFakeParams(params);
			paramsList.add(params);
		}
	}

	private void createFakeParams(double[] params)
	{
		RandomGenerator r = rdg.getRandomGenerator();
		for (int i = 0; i < params.length; i++)
		{
			params[i] = r.nextDouble();
		}
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
