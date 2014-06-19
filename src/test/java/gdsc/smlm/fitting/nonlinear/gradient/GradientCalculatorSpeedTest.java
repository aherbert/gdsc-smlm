package gdsc.smlm.fitting.nonlinear.gradient;

import gdsc.smlm.TestSettings;
import gdsc.smlm.fitting.function.CCDCameraNoiseModel;
import gdsc.smlm.fitting.function.Gaussian2DFunction;
import gdsc.smlm.fitting.function.gaussian.EllipticalGaussian2DFunction;
import gdsc.smlm.fitting.function.gaussian.SingleCircularGaussian2DFunction;
import gdsc.smlm.fitting.function.gaussian.SingleEllipticalGaussian2DFunction;
import gdsc.smlm.fitting.function.gaussian.SingleFixedGaussian2DFunction;
import gdsc.smlm.fitting.function.gaussian.SingleFreeCircularGaussian2DFunction;
import gdsc.smlm.fitting.function.gaussian.SingleNBFixedGaussian2DFunction;
import gdsc.smlm.fitting.nonlinear.gradient.GradientCalculator;
import gdsc.smlm.fitting.nonlinear.gradient.GradientCalculator4;
import gdsc.smlm.fitting.nonlinear.gradient.GradientCalculator6;
import gdsc.smlm.fitting.nonlinear.gradient.GradientCalculator7;
import gdsc.smlm.fitting.utils.DoubleEquality;

import java.util.Random;

import java.util.ArrayList;

import org.junit.Assert;
import org.junit.Test;

/**
 * Contains speed tests for the fastest method for calculating the Hessian and gradient vector
 * for use in NonLinearFit
 */
public class GradientCalculatorSpeedTest
{
	DoubleEquality eq = new DoubleEquality(6, 1e-16);

	int MAX_ITER = 20000;
	int blockWidth = 10;
	float Background = 20;
	float Amplitude = 10;
	float Angle = 0;
	float Xpos = 5;
	float Ypos = 5;
	float Xwidth = 5;
	float Ywidth = 5;

	Random rand;

	@Test
	public void gradientCalculatorFactoryCreatesOptimisedCalculators()
	{
		Assert.assertEquals(GradientCalculatorFactory.newCalculator(3).getClass(), GradientCalculator3.class);		
		Assert.assertEquals(GradientCalculatorFactory.newCalculator(4).getClass(), GradientCalculator4.class);		
		Assert.assertEquals(GradientCalculatorFactory.newCalculator(5).getClass(), GradientCalculator5.class);		
		Assert.assertEquals(GradientCalculatorFactory.newCalculator(6).getClass(), GradientCalculator6.class);		
		Assert.assertEquals(GradientCalculatorFactory.newCalculator(7).getClass(), GradientCalculator7.class);		
		Assert.assertEquals(GradientCalculatorFactory.newCalculator(13).getClass(), GradientCalculator.class);		
	}
	
	@Test
	public void gradientCalculator7ComputesSameAsGradientCalculator()
	{
		gradientCalculatorNComputesSameAsGradientCalculator(new SingleEllipticalGaussian2DFunction(blockWidth), 7);
	}

	@Test
	public void gradientCalculator7IsFasterThanGradientCalculator()
	{
		gradientCalculatorNIsFasterThanGradientCalculator(new SingleEllipticalGaussian2DFunction(blockWidth), 7);
	}
	
	@Test
	public void gradientCalculator6ComputesSameAsGradientCalculator()
	{
		gradientCalculatorNComputesSameAsGradientCalculator(new SingleFreeCircularGaussian2DFunction(blockWidth), 6);
	}

	@Test
	public void gradientCalculator6IsFasterThanGradientCalculator()
	{
		gradientCalculatorNIsFasterThanGradientCalculator(new SingleFreeCircularGaussian2DFunction(blockWidth), 6);
	}
	
	@Test
	public void gradientCalculator5ComputesSameAsGradientCalculator()
	{
		gradientCalculatorNComputesSameAsGradientCalculator(new SingleCircularGaussian2DFunction(blockWidth), 5);
	}

	@Test
	public void gradientCalculator5IsFasterThanGradientCalculator()
	{
		gradientCalculatorNIsFasterThanGradientCalculator(new SingleCircularGaussian2DFunction(blockWidth), 5);
	}
	
	@Test
	public void gradientCalculator4ComputesSameAsGradientCalculator()
	{
		gradientCalculatorNComputesSameAsGradientCalculator(new SingleFixedGaussian2DFunction(blockWidth), 4);
	}

	@Test
	public void gradientCalculator4IsFasterThanGradientCalculator()
	{
		gradientCalculatorNIsFasterThanGradientCalculator(new SingleFixedGaussian2DFunction(blockWidth), 4);
	}
	
	@Test
	public void gradientCalculator3ComputesSameAsGradientCalculator()
	{
		gradientCalculatorNComputesSameAsGradientCalculator(new SingleNBFixedGaussian2DFunction(blockWidth), 3);
	}

	@Test
	public void gradientCalculator3IsFasterThanGradientCalculator()
	{
		gradientCalculatorNIsFasterThanGradientCalculator(new SingleNBFixedGaussian2DFunction(blockWidth), 3);
	}

	private void gradientCalculatorNComputesSameAsGradientCalculator(Gaussian2DFunction func, int nparams)
	{
		// Check the function is the correct size
		Assert.assertEquals(nparams, func.gradientIndices().length);

		int iter = 100;
		rand = new Random(30051977);
		
		double[][] alpha = new double[nparams][nparams];
		double[] beta = new double[nparams];
		double[][] alpha2 = new double[nparams][nparams];
		double[] beta2 = new double[nparams];

		ArrayList<float[]> paramsList = new ArrayList<float[]>(iter);
		ArrayList<float[]> yList = new ArrayList<float[]>(iter);

		int[] x = createData(1, iter, paramsList, yList);

		GradientCalculator calc = new GradientCalculator(beta.length);
		GradientCalculator calc2 = GradientCalculatorFactory.newCalculator(nparams);

		for (int i = 0; i < paramsList.size(); i++)
		{
			calc.findLinearised(x, yList.get(i), paramsList.get(i), alpha, beta, func);
			calc2.findLinearised(x, yList.get(i), paramsList.get(i), alpha2, beta2, func);
			Assert.assertTrue("Observations: Not same beta @ " + i, eq.almostEqualComplement(beta, beta2));
			for (int j = 0; j < beta.length; j++)
				Assert.assertTrue("Observations: Not same alpha @ " + i, eq.almostEqualComplement(alpha[j], alpha2[j]));
		}
		
		for (int i = 0; i < paramsList.size(); i++)
		{
			calc.findLinearised(x.length, yList.get(i), paramsList.get(i), alpha, beta, func);
			calc2.findLinearised(x.length, yList.get(i), paramsList.get(i), alpha2, beta2, func);
			Assert.assertTrue("N-observations: Not same beta @ " + i, eq.almostEqualComplement(beta, beta2));
			for (int j = 0; j < beta.length; j++)
				Assert.assertTrue("N-observations: Not same alpha @ " + i, eq.almostEqualComplement(alpha[j], alpha2[j]));
		} 
		
		func.setNoiseModel(new CCDCameraNoiseModel(10, true));
		
		for (int i = 0; i < paramsList.size(); i++)
		{
			calc.findLinearised(x, yList.get(i), paramsList.get(i), alpha, beta, func);
			calc2.findLinearised(x, yList.get(i), paramsList.get(i), alpha2, beta2, func);
			Assert.assertTrue("Observations+Noise: Not same beta @ " + i, eq.almostEqualComplement(beta, beta2));
			for (int j = 0; j < beta.length; j++)
				Assert.assertTrue("Observations+Noise: Not same alpha @ " + i, eq.almostEqualComplement(alpha[j], alpha2[j]));
		}
		
		for (int i = 0; i < paramsList.size(); i++)
		{
			calc.findLinearised(x.length, yList.get(i), paramsList.get(i), alpha, beta, func);
			calc2.findLinearised(x.length, yList.get(i), paramsList.get(i), alpha2, beta2, func);
			Assert.assertTrue("Observations+Noise: Not same beta @ " + i, eq.almostEqualComplement(beta, beta2));
			for (int j = 0; j < beta.length; j++)
				Assert.assertTrue("Observations+Noise: Not same alpha @ " + i, eq.almostEqualComplement(alpha[j], alpha2[j]));
		}
	}
	
	private void gradientCalculatorNIsFasterThanGradientCalculator(Gaussian2DFunction func, int nparams)
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);
		
		// Check the function is the correct size
		Assert.assertEquals(nparams, func.gradientIndices().length);

		int iter = 10000;
		rand = new Random(30051977);
		double[][] alpha = new double[nparams][nparams];
		double[] beta = new double[nparams];

		ArrayList<float[]> paramsList = new ArrayList<float[]>(iter);
		ArrayList<float[]> yList = new ArrayList<float[]>(iter);

		int[] x = createData(1, iter, paramsList, yList);

		GradientCalculator calc = new GradientCalculator(beta.length);
		GradientCalculator calc2 = GradientCalculatorFactory.newCalculator(nparams);

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

		log("GradientCalculator = %d : GradientCalculator%d = %d : %fx\n", start1, nparams, start2, (1.0 * start1) / start2);
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(start2 < start1);
	}
	
	@Test
	public void gradientCalculatorAssumedXIsFasterThanGradientCalculator()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		int iter = 10000;
		rand = new Random(30051977);
		double[][] alpha = new double[7][7];
		double[] beta = new double[7];

		ArrayList<float[]> paramsList = new ArrayList<float[]>(iter);
		ArrayList<float[]> yList = new ArrayList<float[]>(iter);

		int[] x = createData(1, iter, paramsList, yList);

		GradientCalculator calc = new GradientCalculator6();
		GradientCalculator calc2 = new GradientCalculator6();
		SingleFreeCircularGaussian2DFunction func = new SingleFreeCircularGaussian2DFunction(blockWidth);
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

		log("GradientCalculator = %d : GradientCalculatorAssumed = %d : %fx\n", start1, start2, (1.0 * start1) / start2);
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(start2 < start1);
	}

	/**
	 * Create random elliptical Gaussian data an returns the data plus an estimate of the parameters.
	 * Only the chosen parameters are randomised and returned for a maximum of (background, amplitude, angle, xpos,
	 * ypos, xwidth, ywidth }
	 * 
	 * @param params
	 *            set on output
	 * @return
	 */
	private float[] floatCreateGaussianData(int npeaks, float[] params)
	{
		int n = blockWidth * blockWidth;

		// Generate a 2D Gaussian
		EllipticalGaussian2DFunction func = new EllipticalGaussian2DFunction(npeaks, blockWidth);
		params[0] = Background + rand.nextFloat() * 5f;
		for (int i = 0, j = 1; i < npeaks; i++, j += 6)
		{
			params[j] = Amplitude + rand.nextFloat() * 5f;
			params[j + 1] = 0f; //(float) (Math.PI / 4.0); // Angle
			params[j + 2] = Xpos + rand.nextFloat() * 2f;
			params[j + 3] = Ypos + rand.nextFloat() * 2f;
			params[j + 4] = Xwidth + rand.nextFloat() * 2f;
			params[j + 5] = params[j + 4];
		}

		float[] dy_da = new float[params.length];
		float[] y = new float[n];
		func.initialise(params);
		for (int i = 0; i < y.length; i++)
		{
			// Add random noise
			y[i] = func.eval(i, dy_da) + ((rand.nextFloat() < 0.5f) ? -rand.nextFloat() * 5f : rand.nextFloat() * 5f);
		}

		// Randomise only the necessary parameters (i.e. not angle and X & Y widths should be the same)
		params[0] += ((rand.nextFloat() < 0.5f) ? -rand.nextFloat() : rand.nextFloat());
		for (int i = 0, j = 1; i < npeaks; i++, j += 6)
		{
			params[j + 1] += ((rand.nextFloat() < 0.5f) ? -rand.nextFloat() : rand.nextFloat());
			params[j + 3] += ((rand.nextFloat() < 0.5f) ? -rand.nextFloat() : rand.nextFloat());
			params[j + 4] += ((rand.nextFloat() < 0.5f) ? -rand.nextFloat() : rand.nextFloat());
			params[j + 5] = params[j + 4];
		}

		return y;
	}

	protected int[] createData(int npeaks, int iter, ArrayList<float[]> paramsList, ArrayList<float[]> yList)
	{
		int[] x = new int[blockWidth * blockWidth];
		for (int i = 0; i < x.length; i++)
			x[i] = i;
		for (int i = 0; i < iter; i++)
		{
			float[] params = new float[1 + 6 * npeaks];
			float[] y = floatCreateGaussianData(npeaks, params);
			paramsList.add(params);
			yList.add(y);
		}
		return x;
	}

	protected ArrayList<double[]> copyList(ArrayList<float[]> paramsList)
	{
		ArrayList<double[]> params2List = new ArrayList<double[]>(paramsList.size());
		for (int i = 0; i < paramsList.size(); i++)
		{
			params2List.add(copydouble(paramsList.get(i)));
		}
		return params2List;
	}

	private double[] copydouble(float[] d)
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
