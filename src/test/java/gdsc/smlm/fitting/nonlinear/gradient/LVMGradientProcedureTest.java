package gdsc.smlm.fitting.nonlinear.gradient;

import java.util.ArrayList;
import java.util.Arrays;

import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.apache.commons.math3.util.Precision;
import org.ejml.data.DenseMatrix64F;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.ij.Utils;
import gdsc.core.utils.DoubleEquality;
import gdsc.core.utils.Maths;
import gdsc.smlm.TestSettings;
import gdsc.smlm.fitting.nonlinear.gradient.LVMGradientProcedureFactory.Type;
import gdsc.smlm.function.DummyGradientFunction;
import gdsc.smlm.function.FakeGradientFunction;
import gdsc.smlm.function.Gradient1Function;
import gdsc.smlm.function.PrecomputedGradient1Function;
import gdsc.smlm.function.ValueProcedure;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import gdsc.smlm.function.gaussian.erf.ErfGaussian2DFunction;
import gdsc.smlm.function.gaussian.erf.SingleFreeCircularErfGaussian2DFunction;

/**
 * Contains speed tests for the methods for calculating the Hessian and gradient vector
 * for use in the LVM algorithm.
 */
public class LVMGradientProcedureTest
{
	boolean speedTests = true;
	DoubleEquality eq = new DoubleEquality(1e-6, 1e-16);

	int MAX_ITER = 20000;
	int blockWidth = 10;
	double Noise = 0.3;
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
		DummyGradientFunction[] f = new DummyGradientFunction[7];
		for (int i = 1; i < f.length; i++)
			f[i] = new DummyGradientFunction(i);

		LVMGradientProcedureFactory.Type MLE = LVMGradientProcedureFactory.Type.MLE;
		LVMGradientProcedureFactory.Type WLSQ = LVMGradientProcedureFactory.Type.WLSQ;
		LVMGradientProcedureFactory.Type LSQ = LVMGradientProcedureFactory.Type.LSQ;

		//@formatter:off
		
		// Generic factory
		double[] y = new double[0];
		double[] b = null;
		Assert.assertEquals(LVMGradientProcedureFactory.create(y, f[6], LSQ).getClass(), LSQLVMGradientProcedure6.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y, f[5], LSQ).getClass(), LSQLVMGradientProcedure5.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y, f[4], LSQ).getClass(), LSQLVMGradientProcedure4.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y, f[1], LSQ).getClass(), LSQLVMGradientProcedure.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y, f[6], MLE).getClass(), MLELVMGradientProcedure6.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y, f[5], MLE).getClass(), MLELVMGradientProcedure5.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y, f[4], MLE).getClass(), MLELVMGradientProcedure4.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y, f[1], MLE).getClass(), MLELVMGradientProcedure.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y, f[6], WLSQ).getClass(), WLSQLVMGradientProcedure6.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y, f[5], WLSQ).getClass(), WLSQLVMGradientProcedure5.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y, f[4], WLSQ).getClass(), WLSQLVMGradientProcedure4.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y, f[1], WLSQ).getClass(), WLSQLVMGradientProcedure.class);
		
		// Handle null background
		b = null;
		Assert.assertEquals(LVMGradientProcedureFactory.create(y, b, f[6], LSQ).getClass(), LSQLVMGradientProcedure6.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y, b, f[5], LSQ).getClass(), LSQLVMGradientProcedure5.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y, b, f[4], LSQ).getClass(), LSQLVMGradientProcedure4.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y, b, f[1], LSQ).getClass(), LSQLVMGradientProcedure.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y, b, f[6], MLE).getClass(), MLELVMGradientProcedure6.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y, b, f[5], MLE).getClass(), MLELVMGradientProcedure5.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y, b, f[4], MLE).getClass(), MLELVMGradientProcedure4.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y, b, f[1], MLE).getClass(), MLELVMGradientProcedure.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y, b, f[6], WLSQ).getClass(), WLSQLVMGradientProcedure6.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y, b, f[5], WLSQ).getClass(), WLSQLVMGradientProcedure5.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y, b, f[4], WLSQ).getClass(), WLSQLVMGradientProcedure4.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y, b, f[1], WLSQ).getClass(), WLSQLVMGradientProcedure.class);
		
		// Handle the pre-computed background
		b = new double[f[6].size()];
		Assert.assertEquals(LVMGradientProcedureFactory.create(y, b, f[6], LSQ).getClass(), LSQLVMGradientProcedure6.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y, b, f[5], LSQ).getClass(), LSQLVMGradientProcedure5.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y, b, f[4], LSQ).getClass(), LSQLVMGradientProcedure4.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y, b, f[1], LSQ).getClass(), LSQLVMGradientProcedure.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y, b, f[6], MLE).getClass(), MLELVMGradientProcedure6.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y, b, f[5], MLE).getClass(), MLELVMGradientProcedure5.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y, b, f[4], MLE).getClass(), MLELVMGradientProcedure4.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y, b, f[1], MLE).getClass(), MLELVMGradientProcedure.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y, b, f[6], WLSQ).getClass(), WLSQLVMGradientProcedure6.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y, b, f[5], WLSQ).getClass(), WLSQLVMGradientProcedure5.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y, b, f[4], WLSQ).getClass(), WLSQLVMGradientProcedure4.class);
		Assert.assertEquals(LVMGradientProcedureFactory.create(y, b, f[1], WLSQ).getClass(), WLSQLVMGradientProcedure.class);

		// Dedicated factories
		Assert.assertEquals(LSQLVMGradientProcedureFactory.create(y, f[6]).getClass(), LSQLVMGradientProcedure6.class);
		Assert.assertEquals(LSQLVMGradientProcedureFactory.create(y, f[5]).getClass(), LSQLVMGradientProcedure5.class);
		Assert.assertEquals(LSQLVMGradientProcedureFactory.create(y, f[4]).getClass(), LSQLVMGradientProcedure4.class);
		Assert.assertEquals(LSQLVMGradientProcedureFactory.create(y, f[1]).getClass(), LSQLVMGradientProcedure.class);
		Assert.assertEquals(MLELVMGradientProcedureFactory.create(y, f[6]).getClass(), MLELVMGradientProcedure6.class);
		Assert.assertEquals(MLELVMGradientProcedureFactory.create(y, f[5]).getClass(), MLELVMGradientProcedure5.class);
		Assert.assertEquals(MLELVMGradientProcedureFactory.create(y, f[4]).getClass(), MLELVMGradientProcedure4.class);
		Assert.assertEquals(MLELVMGradientProcedureFactory.create(y, f[1]).getClass(), MLELVMGradientProcedure.class);
		Assert.assertEquals(WLSQLVMGradientProcedureFactory.create(y, null, f[6]).getClass(), WLSQLVMGradientProcedure6.class);
		Assert.assertEquals(WLSQLVMGradientProcedureFactory.create(y, null, f[5]).getClass(), WLSQLVMGradientProcedure5.class);
		Assert.assertEquals(WLSQLVMGradientProcedureFactory.create(y, null, f[4]).getClass(), WLSQLVMGradientProcedure4.class);
		Assert.assertEquals(WLSQLVMGradientProcedureFactory.create(y, null, f[1]).getClass(), WLSQLVMGradientProcedure.class);

		// Handle null background
		b = null;
		Assert.assertEquals(LSQLVMGradientProcedureFactory.create(y, b, f[6]).getClass(), LSQLVMGradientProcedure6.class);
		Assert.assertEquals(LSQLVMGradientProcedureFactory.create(y, b, f[5]).getClass(), LSQLVMGradientProcedure5.class);
		Assert.assertEquals(LSQLVMGradientProcedureFactory.create(y, b, f[4]).getClass(), LSQLVMGradientProcedure4.class);
		Assert.assertEquals(LSQLVMGradientProcedureFactory.create(y, b, f[1]).getClass(), LSQLVMGradientProcedure.class);
		Assert.assertEquals(MLELVMGradientProcedureFactory.create(y, b, f[6]).getClass(), MLELVMGradientProcedure6.class);
		Assert.assertEquals(MLELVMGradientProcedureFactory.create(y, b, f[5]).getClass(), MLELVMGradientProcedure5.class);
		Assert.assertEquals(MLELVMGradientProcedureFactory.create(y, b, f[4]).getClass(), MLELVMGradientProcedure4.class);
		Assert.assertEquals(MLELVMGradientProcedureFactory.create(y, b, f[1]).getClass(), MLELVMGradientProcedure.class);
		Assert.assertEquals(WLSQLVMGradientProcedureFactory.create(y, b, null, f[6]).getClass(), WLSQLVMGradientProcedure6.class);
		Assert.assertEquals(WLSQLVMGradientProcedureFactory.create(y, b, null, f[5]).getClass(), WLSQLVMGradientProcedure5.class);
		Assert.assertEquals(WLSQLVMGradientProcedureFactory.create(y, b, null, f[4]).getClass(), WLSQLVMGradientProcedure4.class);
		Assert.assertEquals(WLSQLVMGradientProcedureFactory.create(y, b, null, f[1]).getClass(), WLSQLVMGradientProcedure.class);

		// Handle the pre-computed background
		b = new double[f[6].size()];
		Assert.assertEquals(LSQLVMGradientProcedureFactory.create(y, b, f[6]).getClass(), LSQLVMGradientProcedure6.class);
		Assert.assertEquals(LSQLVMGradientProcedureFactory.create(y, b, f[5]).getClass(), LSQLVMGradientProcedure5.class);
		Assert.assertEquals(LSQLVMGradientProcedureFactory.create(y, b, f[4]).getClass(), LSQLVMGradientProcedure4.class);
		Assert.assertEquals(LSQLVMGradientProcedureFactory.create(y, b, f[1]).getClass(), LSQLVMGradientProcedure.class);
		Assert.assertEquals(MLELVMGradientProcedureFactory.create(y, b, f[6]).getClass(), MLELVMGradientProcedure6.class);
		Assert.assertEquals(MLELVMGradientProcedureFactory.create(y, b, f[5]).getClass(), MLELVMGradientProcedure5.class);
		Assert.assertEquals(MLELVMGradientProcedureFactory.create(y, b, f[4]).getClass(), MLELVMGradientProcedure4.class);
		Assert.assertEquals(MLELVMGradientProcedureFactory.create(y, b, f[1]).getClass(), MLELVMGradientProcedure.class);
		Assert.assertEquals(WLSQLVMGradientProcedureFactory.create(y, b, null, f[6]).getClass(), WLSQLVMGradientProcedure6.class);
		Assert.assertEquals(WLSQLVMGradientProcedureFactory.create(y, b, null, f[5]).getClass(), WLSQLVMGradientProcedure5.class);
		Assert.assertEquals(WLSQLVMGradientProcedureFactory.create(y, b, null, f[4]).getClass(), WLSQLVMGradientProcedure4.class);
		Assert.assertEquals(WLSQLVMGradientProcedureFactory.create(y, b, null, f[1]).getClass(), WLSQLVMGradientProcedure.class);
		
		//@formatter:on
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
			LVMGradientProcedure p = LVMGradientProcedureFactory.create(yList.get(i), func,
					(mle) ? Type.MLE : Type.LSQ);
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

		final Type type = (mle) ? Type.MLE : Type.LSQ;

		for (int i = 0; i < paramsList.size(); i++)
		{
			LVMGradientProcedure p = LVMGradientProcedureFactory.create(yList.get(i), func, type);
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
					LVMGradientProcedure p = LVMGradientProcedureFactory.create(yList.get(i), func, type);
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
		gradientProcedureUnrolledComputesSameAsGradientProcedure(Type.LSQ, false);
	}

	@Test
	public void gradientProcedureMLEUnrolledComputesSameAsGradientProcedure()
	{
		gradientProcedureUnrolledComputesSameAsGradientProcedure(Type.MLE, false);
	}

	@Test
	public void gradientProcedureWLSQUnrolledComputesSameAsGradientProcedure()
	{
		gradientProcedureUnrolledComputesSameAsGradientProcedure(Type.WLSQ, false);
	}

	@Test
	public void gradientProcedureLSQUnrolledComputesSameAsGradientProcedureWithPrecomputed()
	{
		gradientProcedureUnrolledComputesSameAsGradientProcedure(Type.LSQ, true);
	}

	@Test
	public void gradientProcedureMLEUnrolledComputesSameAsGradientProcedureWithPrecomputed()
	{
		gradientProcedureUnrolledComputesSameAsGradientProcedure(Type.MLE, true);
	}

	@Test
	public void gradientProcedureWLSQUnrolledComputesSameAsGradientProcedureWithPrecomputed()
	{
		gradientProcedureUnrolledComputesSameAsGradientProcedure(Type.WLSQ, true);
	}

	private void gradientProcedureUnrolledComputesSameAsGradientProcedure(Type type, boolean precomputed)
	{
		gradientProcedureUnrolledComputesSameAsGradientProcedure(4, type, precomputed);
		gradientProcedureUnrolledComputesSameAsGradientProcedure(5, type, precomputed);
		gradientProcedureUnrolledComputesSameAsGradientProcedure(6, type, precomputed);
	}

	private void gradientProcedureUnrolledComputesSameAsGradientProcedure(int nparams, Type type, boolean precomputed)
	{
		int iter = 10;
		rdg = new RandomDataGenerator(new Well19937c(30051977));

		ArrayList<double[]> paramsList = new ArrayList<double[]>(iter);
		ArrayList<double[]> yList = new ArrayList<double[]>(iter);

		createFakeData(nparams, iter, paramsList, yList);
		Gradient1Function func = new FakeGradientFunction(blockWidth, nparams);

		if (precomputed)
		{
			double[] b = Utils.newArray(func.size(), 0.1, 1.3);
			func = new PrecomputedGradient1Function(func, b);
		}

		String name = String.format("[%d] %b", nparams, type);
		for (int i = 0; i < paramsList.size(); i++)
		{
			LVMGradientProcedure p1 = createProcedure(type, yList.get(i), func);
			p1.gradient(paramsList.get(i));

			LVMGradientProcedure p2 = LVMGradientProcedureFactory.create(yList.get(i), func, type);
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

	private LVMGradientProcedure createProcedure(Type type, double[] y, Gradient1Function func)
	{
		switch (type)
		{
			case MLE:
				return new MLELVMGradientProcedure(y, func);
			case WLSQ:
				return new WLSQLVMGradientProcedure(y, null, func);
			default:
				return new LSQLVMGradientProcedure(y, func);
		}
	}

	@Test
	public void gradientProcedureLSQIsFasterUnrolledThanGradientProcedure()
	{
		gradientProcedureIsFasterUnrolledThanGradientProcedure(Type.LSQ, false);
	}

	@Test
	public void gradientProcedureMLEIsFasterUnrolledThanGradientProcedure()
	{
		gradientProcedureIsFasterUnrolledThanGradientProcedure(Type.MLE, false);
	}

	@Test
	public void gradientProcedureWLSQIsFasterUnrolledThanGradientProcedure()
	{
		gradientProcedureIsFasterUnrolledThanGradientProcedure(Type.WLSQ, false);
	}

	@Test
	public void gradientProcedureLSQIsFasterUnrolledThanGradientProcedureWithPrecomputed()
	{
		gradientProcedureIsFasterUnrolledThanGradientProcedure(Type.LSQ, true);
	}

	@Test
	public void gradientProcedureMLEIsFasterUnrolledThanGradientProcedureWithPrecomputed()
	{
		gradientProcedureIsFasterUnrolledThanGradientProcedure(Type.MLE, true);
	}

	@Test
	public void gradientProcedureWLSQIsFasterUnrolledThanGradientProcedureWithPrecomputed()
	{
		gradientProcedureIsFasterUnrolledThanGradientProcedure(Type.WLSQ, true);
	}

	private void gradientProcedureIsFasterUnrolledThanGradientProcedure(Type type, boolean precomputed)
	{
		gradientProcedureLinearIsFasterThanGradientProcedure(4, type, precomputed);
		gradientProcedureLinearIsFasterThanGradientProcedure(5, type, precomputed);
		gradientProcedureLinearIsFasterThanGradientProcedure(6, type, precomputed);
	}

	private void gradientProcedureLinearIsFasterThanGradientProcedure(final int nparams, final Type type,
			final boolean precomputed)
	{
		org.junit.Assume.assumeTrue(speedTests || TestSettings.RUN_SPEED_TESTS);

		final int iter = 100;
		rdg = new RandomDataGenerator(new Well19937c(30051977));

		final ArrayList<double[]> paramsList = new ArrayList<double[]>(iter);
		final ArrayList<double[]> yList = new ArrayList<double[]>(iter);

		createData(1, iter, paramsList, yList);

		// Remove the timing of the function call by creating a dummy function
		FakeGradientFunction fgf = new FakeGradientFunction(blockWidth, nparams);
		final Gradient1Function func;
		if (precomputed)
		{
			final double[] b = Utils.newArray(fgf.size(), 0.1, 1.3);
			func = new PrecomputedGradient1Function(fgf, b);
		}
		else
		{
			func = fgf;
		}

		for (int i = 0; i < paramsList.size(); i++)
		{
			LVMGradientProcedure p1 = createProcedure(type, yList.get(i), func);
			p1.gradient(paramsList.get(i));
			p1.gradient(paramsList.get(i));

			LVMGradientProcedure p2 = LVMGradientProcedureFactory.create(yList.get(i), func, type);
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
					LVMGradientProcedure p1 = createProcedure(type, yList.get(i), func);
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
					LVMGradientProcedure p2 = LVMGradientProcedureFactory.create(yList.get(i), func, type);
					for (int j = loops; j-- > 0;)
						p2.gradient(paramsList.get(k++ % iter));
				}
			}
		};
		long time2 = t2.getTime();

		log("%s, Precomputed=%b : Standard = %d : Unrolled %d = %d : %fx\n", type, precomputed, time1, nparams, time2,
				(1.0 * time1) / time2);
		Assert.assertTrue(time2 < time1);
	}

	@Test
	public void gradientProcedureLSQComputesGradient()
	{
		gradientProcedureComputesGradient(new SingleFreeCircularErfGaussian2DFunction(blockWidth, blockWidth), Type.LSQ,
				false);
	}

	@Test
	public void gradientProcedureMLEComputesGradient()
	{
		gradientProcedureComputesGradient(new SingleFreeCircularErfGaussian2DFunction(blockWidth, blockWidth), Type.MLE,
				false);
	}

	@Test
	public void gradientProcedureWLSQComputesGradient()
	{
		gradientProcedureComputesGradient(new SingleFreeCircularErfGaussian2DFunction(blockWidth, blockWidth),
				Type.WLSQ, false);
	}

	@Test
	public void gradientProcedureLSQComputesGradientWithPrecomputed()
	{
		gradientProcedureComputesGradient(new SingleFreeCircularErfGaussian2DFunction(blockWidth, blockWidth), Type.LSQ,
				true);
	}

	@Test
	public void gradientProcedureMLEComputesGradientWithPrecomputed()
	{
		gradientProcedureComputesGradient(new SingleFreeCircularErfGaussian2DFunction(blockWidth, blockWidth), Type.MLE,
				true);
	}

	@Test
	public void gradientProcedureWLSQComputesGradientWithPrecomputed()
	{
		gradientProcedureComputesGradient(new SingleFreeCircularErfGaussian2DFunction(blockWidth, blockWidth),
				Type.WLSQ, true);
	}

	private void gradientProcedureComputesGradient(ErfGaussian2DFunction func, Type type, boolean precomputed)
	{
		int nparams = func.getNumberOfGradients();
		int[] indices = func.gradientIndices();

		int iter = 100;
		rdg = new RandomDataGenerator(new Well19937c(30051977));

		ArrayList<double[]> paramsList = new ArrayList<double[]>(iter);
		ArrayList<double[]> yList = new ArrayList<double[]>(iter);

		createData(1, iter, paramsList, yList, true);

		double delta = 1e-3;
		DoubleEquality eq = new DoubleEquality(1e-3, 1e-3);
		final double[] b = (precomputed) ? new double[func.size()] : null;

		for (int i = 0; i < paramsList.size(); i++)
		{
			double[] y = yList.get(i);
			double[] a = paramsList.get(i);
			double[] a2 = a.clone();

			if (precomputed)
			{
				// Mock fitting part of the function already
				for (int j = 0; j < b.length; j++)
				{
					b[j] = y[j] * 0.5;
				}
			}

			LVMGradientProcedure p = LVMGradientProcedureFactory.create(y, b, func, type);
			p.gradient(a);
			//double s = p.value;
			double[] beta = p.beta.clone();
			for (int j = 0; j < nparams; j++)
			{
				int k = indices[j];
				double d = Precision.representableDelta(a[k], (a[k] == 0) ? 1e-3 : a[k] * delta);
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
				Assert.assertTrue("Not same gradient @ " + j, eq.almostEqualRelativeOrAbsolute(beta[j], gradient));
			}
		}
	}

	@Test
	public void gradientProcedureLSQSupportsPrecomputed()
	{
		gradientProcedureSupportsPrecomputed(Type.LSQ);
	}

	@Test
	public void gradientProcedureMLESupportsPrecomputed()
	{
		gradientProcedureSupportsPrecomputed(Type.MLE);
	}

	@Test
	public void gradientProcedureWLSQSupportsPrecomputed()
	{
		gradientProcedureSupportsPrecomputed(Type.WLSQ);
	}

	private void gradientProcedureSupportsPrecomputed(final Type type)
	{
		int iter = 10;
		rdg = new RandomDataGenerator(new Well19937c(30051977));

		ArrayList<double[]> paramsList = new ArrayList<double[]>(iter);
		ArrayList<double[]> yList = new ArrayList<double[]>(iter);

		// 3 peaks
		createData(3, iter, paramsList, yList, true);

		for (int i = 0; i < paramsList.size(); i++)
		{
			final double[] y = yList.get(i);
			// Add Gaussian read noise so we have negatives
			double min = Maths.min(y);
			for (int j = 0; j < y.length; j++)
				y[j] = y[i] - min + rdg.nextGaussian(0, Noise);
		}

		// We want to know that:
		// y|peak1+peak2+peak3 == y|peak1+peak2+peak3(precomputed)
		// We want to know when:
		// y|peak1+peak2+peak3 != y-peak3|peak1+peak2
		// i.e. we cannot subtract a precomputed peak from the data, it must be included in the fit
		// E.G. LSQ - subtraction is OK, MLE/WLSQ - subtraction is not allowed

		Gaussian2DFunction f123 = GaussianFunctionFactory.create2D(3, blockWidth, blockWidth,
				GaussianFunctionFactory.FIT_ERF_FREE_CIRCLE, null);
		Gaussian2DFunction f12 = GaussianFunctionFactory.create2D(2, blockWidth, blockWidth,
				GaussianFunctionFactory.FIT_ERF_FREE_CIRCLE, null);
		Gaussian2DFunction f3 = GaussianFunctionFactory.create2D(1, blockWidth, blockWidth,
				GaussianFunctionFactory.FIT_ERF_FREE_CIRCLE, null);

		int nparams = f12.getNumberOfGradients();
		int[] indices = f12.gradientIndices();
		final double[] b = new double[f12.size()];

		double delta = 1e-3;
		DoubleEquality eq = new DoubleEquality(5e-3, 1e-6);
		double[] a1peaks = new double[7];
		final double[] y_b = new double[b.length];

		for (int i = 0; i < paramsList.size(); i++)
		{
			final double[] y = yList.get(i);
			double[] a3peaks = paramsList.get(i);
			double[] a2peaks = Arrays.copyOf(a3peaks, 1 + 2 * 6);
			double[] a2peaks2 = a2peaks.clone();
			for (int j = 1; j < 7; j++)
				a1peaks[j] = a3peaks[j + 2 * 6];

			// Evaluate peak 3 to get the background and subtract it from the data to get the new data
			f3.initialise0(a1peaks);
			f3.forEach(new ValueProcedure()
			{
				int k = 0;

				public void execute(double value)
				{
					b[k] = value;
					// Remove negatives for MLE
					if (type == Type.MLE)
					{
						y[k] = Math.max(0, y[k]);
						y_b[k] = Math.max(0, y[k] - value);
					}
					else
					{
						y_b[k] = y[k] - value;
					}
					k++;
				}
			});

			// These should be the same
			LVMGradientProcedure p123 = LVMGradientProcedureFactory.create(y, f123, type);
			LVMGradientProcedure p12b3 = LVMGradientProcedureFactory.create(y, b, f12, type);
			// This may be different
			LVMGradientProcedure p12m3 = LVMGradientProcedureFactory.create(y_b, f12, type);

			// Check they are the same
			p123.gradient(a3peaks);
			double[][] m123 = p123.getAlphaMatrix();

			p12b3.gradient(a2peaks);
			double s = p12b3.value;
			double[] beta = p12b3.beta.clone();
			double[][] alpha = p12b3.getAlphaMatrix();

			System.out.printf("MLE=%b [%d] p12b3  %f  %f\n", type, i, p123.value, s);

			Assert.assertTrue("p12b3 Not same value @ " + i, eq.almostEqualRelativeOrAbsolute(p123.value, s));
			Assert.assertTrue("p12b3 Not same gradient @ " + i, eq.almostEqualRelativeOrAbsolute(beta, p123.beta));
			for (int j = 0; j < alpha.length; j++)
				Assert.assertTrue("p12b3 Not same alpha @ " + j, eq.almostEqualRelativeOrAbsolute(alpha[j], m123[j]));

			// Check actual gradients are correct
			for (int j = 0; j < nparams; j++)
			{
				int k = indices[j];
				double d = Precision.representableDelta(a2peaks[k], (a2peaks[k] == 0) ? 1e-3 : a2peaks[k] * delta);
				a2peaks2[k] = a2peaks[k] + d;
				p12b3.value(a2peaks2);
				double s1 = p12b3.value;
				a2peaks2[k] = a2peaks[k] - d;
				p12b3.value(a2peaks2);
				double s2 = p12b3.value;
				a2peaks2[k] = a2peaks[k];

				// Apply a factor of -2 to compute the actual gradients:
				// See Numerical Recipes in C++, 2nd Ed. Equation 15.5.6 for Nonlinear Models
				beta[j] *= -2;

				double gradient = (s1 - s2) / (2 * d);
				System.out.printf("[%d,%d] %f  (%s %f+/-%f)  %f  ?=  %f  (%f)\n", i, k, s, f12.getName(k), a2peaks[k],
						d, beta[j], gradient, DoubleEquality.relativeError(gradient, beta[j]));
				Assert.assertTrue("Not same gradient @ " + j, eq.almostEqualRelativeOrAbsolute(beta[j], gradient));
			}

			// Check these may be different
			p12m3.gradient(a2peaks);
			s = p12m3.value;
			beta = p12m3.beta.clone();
			alpha = p12m3.getAlphaMatrix();

			System.out.printf("%s [%d] p12m3  %f  %f\n", type, i, p123.value, s);

			if (type != Type.LSQ)
			{
				Assert.assertFalse("p12b3 Same value @ " + i, eq.almostEqualRelativeOrAbsolute(p123.value, s));
				Assert.assertFalse("p12b3 Same gradient @ " + i, eq.almostEqualRelativeOrAbsolute(beta, p123.beta));
				for (int j = 0; j < alpha.length; j++)
				{
					//System.out.printf("%s != %s\n", Arrays.toString(alpha[j]), Arrays.toString(m123[j]));
					Assert.assertFalse("p12b3 Same alpha @ " + j, eq.almostEqualRelativeOrAbsolute(alpha[j], m123[j]));
				}
			}
			else
			{
				Assert.assertTrue("p12b3 Not same value @ " + i, eq.almostEqualRelativeOrAbsolute(p123.value, s));
				Assert.assertTrue("p12b3 Not same gradient @ " + i, eq.almostEqualRelativeOrAbsolute(beta, p123.beta));
				for (int j = 0; j < alpha.length; j++)
					Assert.assertTrue("p12b3 Not same alpha @ " + j,
							eq.almostEqualRelativeOrAbsolute(alpha[j], m123[j]));
			}

			// Check actual gradients are correct
			for (int j = 0; j < nparams; j++)
			{
				int k = indices[j];
				double d = Precision.representableDelta(a2peaks[k], (a2peaks[k] == 0) ? 1e-3 : a2peaks[k] * delta);
				a2peaks2[k] = a2peaks[k] + d;
				p12m3.value(a2peaks2);
				double s1 = p12m3.value;
				a2peaks2[k] = a2peaks[k] - d;
				p12m3.value(a2peaks2);
				double s2 = p12m3.value;
				a2peaks2[k] = a2peaks[k];

				// Apply a factor of -2 to compute the actual gradients:
				// See Numerical Recipes in C++, 2nd Ed. Equation 15.5.6 for Nonlinear Models
				beta[j] *= -2;

				double gradient = (s1 - s2) / (2 * d);
				System.out.printf("[%d,%d] %f  (%s %f+/-%f)  %f  ?=  %f  (%f)\n", i, k, s, f12.getName(k), a2peaks[k],
						d, beta[j], gradient, DoubleEquality.relativeError(gradient, beta[j]));
				Assert.assertTrue("Not same gradient @ " + j, eq.almostEqualRelativeOrAbsolute(beta[j], gradient));
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
		ErfGaussian2DFunction func = (ErfGaussian2DFunction) GaussianFunctionFactory.create2D(npeaks, blockWidth,
				blockWidth, GaussianFunctionFactory.FIT_ERF_FREE_CIRCLE, null);
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
