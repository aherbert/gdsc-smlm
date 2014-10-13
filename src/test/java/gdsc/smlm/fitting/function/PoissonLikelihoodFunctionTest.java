package gdsc.smlm.fitting.function;

import gdsc.smlm.fitting.utils.DoubleEquality;

import java.util.Arrays;

import org.junit.Assert;
import org.junit.Test;

public class PoissonLikelihoodFunctionTest
{
	DoubleEquality eqPerDatum = new DoubleEquality(2, 0.01);
	DoubleEquality eq = new DoubleEquality(3, 0.001);

	String[] NAME = { "BACKGROUND", "AMPLITUDE", "ANGLE", "X_POSITION", "Y_POSITION", "X_SD", "Y_SD" };

	// Compute as per Numerical Recipes 5.7.
	// Approximate error accuracy in single precision: Ef
	// Step size for deerivatives:
	// h ~ (Ef)^(1/3) * xc
	// xc is the characteristic scale over which x changes, assumed to be 1 (not x as per NR since x is close to zero)
	final float h_ = (float) (Math.pow(1e-3f, 1.0 / 3));

	int[] testx = new int[] { 4, 5, 6 };
	int[] testy = new int[] { 4, 5, 6 };
	// Do not test zero background since this is an edge case for the likelihood function
	float[] testbackground_ = new float[] { 10, 400 };

	float[] testamplitude1_ = new float[] { 15, 55, 105 };
	float[] testangle1_ = new float[] { (float) (Math.PI / 5), (float) (Math.PI / 3) };
	float[] testcx1_ = new float[] { 4.9f, 5.3f };
	float[] testcy1_ = new float[] { 4.8f, 5.1f };
	float[][] testw1_ = new float[][] { { 1.1f, 1.4f }, { 1.1f, 1.7f }, { 1.5f, 1.2f }, { 1.3f, 1.7f }, };

	float[] testbackground, testamplitude1, testangle1, testcx1, testcy1;
	float[][] testw1;

	int maxx = 10;
	float background = 50;
	float angle = 0;
	float width = 5;

	@Test
	public void fitFixedComputesGradientPerDatum()
	{
		functionComputesGradientPerDatum(GaussianFunctionFactory.FIT_FIXED);
	}

	@Test
	public void fitCircleComputesGradientPerDatum()
	{
		functionComputesGradientPerDatum(GaussianFunctionFactory.FIT_CIRCLE);
	}

	@Test
	public void fitFreeCircleComputesGradientPerDatum()
	{
		functionComputesGradientPerDatum(GaussianFunctionFactory.FIT_FREE_CIRCLE);
	}

	@Test
	public void fitEllipticalComputesGradientPerDatum()
	{
		functionComputesGradientPerDatum(GaussianFunctionFactory.FIT_ELLIPTICAL);
	}

	@Test
	public void fitNBFixedComputesGradientPerDatum()
	{
		functionComputesGradientPerDatum(GaussianFunctionFactory.FIT_NB_FIXED);
	}

	@Test
	public void fitNBCircleComputesGradientPerDatum()
	{
		functionComputesGradientPerDatum(GaussianFunctionFactory.FIT_NB_CIRCLE);
	}

	@Test
	public void fitNBFreeCircleComputesGradientPerDatum()
	{
		functionComputesGradientPerDatum(GaussianFunctionFactory.FIT_NB_FREE_CIRCLE);
	}

	@Test
	public void fitNBEllipticalComputesGradientPerDatum()
	{
		functionComputesGradientPerDatum(GaussianFunctionFactory.FIT_NB_ELLIPTICAL);
	}

	private void functionComputesGradientPerDatum(int flags)
	{
		Gaussian2DFunction f1 = GaussianFunctionFactory.create2D(1, maxx, flags);
		// Setup
		testbackground = testbackground_;
		testamplitude1 = testamplitude1_;
		testangle1 = testangle1_;
		testcx1 = testcx1_;
		testcy1 = testcy1_;
		testw1 = testw1_;
		if (!f1.evaluatesBackground())
		{
			testbackground = new float[] { testbackground[0] };
		}
		if (!f1.evaluatesAmplitude())
		{
			testamplitude1 = new float[] { testamplitude1[0] };
		}
		if (!f1.evaluatesAngle())
		{
			testangle1 = new float[] { 0 };
		}
		// Position is always evaluated

		boolean noSecondWidth = false;
		if (!f1.evaluatesSD0())
		{
			// Just use 1 width
			testw1 = new float[][] { testw1[0] };
			// If no width 0 then assume we have no width 1 as well
			noSecondWidth = true;
		}
		else if (!f1.evaluatesSD1())
		{
			// No evaluation of second width needs only variation in width 0 so truncate 
			testw1 = Arrays.copyOf(testw1, 2);
			noSecondWidth = true;
		}
		if (noSecondWidth)
		{
			for (int i = 0; i < testw1.length; i++)
			{
				testw1[i][1] = testw1[i][0];
			}
		}

		if (f1.evaluatesBackground())
			functionComputesTargetGradientPerDatum(f1, Gaussian2DFunction.BACKGROUND);
		if (f1.evaluatesAmplitude())
			functionComputesTargetGradientPerDatum(f1, Gaussian2DFunction.AMPLITUDE);
		if (f1.evaluatesAngle())
			functionComputesTargetGradientPerDatum(f1, Gaussian2DFunction.ANGLE);
		functionComputesTargetGradientPerDatum(f1, Gaussian2DFunction.X_POSITION);
		functionComputesTargetGradientPerDatum(f1, Gaussian2DFunction.Y_POSITION);
		if (f1.evaluatesSD0())
			functionComputesTargetGradientPerDatum(f1, Gaussian2DFunction.X_SD);
		if (f1.evaluatesSD1())
			functionComputesTargetGradientPerDatum(f1, Gaussian2DFunction.Y_SD);
	}

	private void functionComputesTargetGradientPerDatum(Gaussian2DFunction f1, int targetParameter)
	{
		int[] indices = f1.gradientIndices();
		int gradientIndex = findGradientIndex(f1, targetParameter);
		double[] dyda = new double[indices.length];
		double[] dyda2 = new double[indices.length];
		float[] a;

		PoissonLikelihoodFunction ff1;

		int n = maxx * maxx;
		int count = 0, total = 0;

		for (float background : testbackground)
			for (float amplitude1 : testamplitude1)
				for (float angle1 : testangle1)
					for (float cx1 : testcx1)
						for (float cy1 : testcy1)
							for (float[] w1 : testw1)
							{
								a = createParameters(background, amplitude1, angle1, cx1, cy1, w1[0], w1[1]);

								// Create y as a function we would want to move towards
								float[] a2 = createParameters(background, amplitude1, angle1, cx1, cy1, w1[0], w1[1]);
								a2[targetParameter] *= 1.3;
								f1.initialise(a2);
								float[] data = new float[maxx * maxx];
								for (int i = 0; i < n; i++)
									data[i] = f1.eval(i);

								ff1 = new PoissonLikelihoodFunction(f1, a, data, n);

								// Numerically solve gradient. 
								// Calculate the step size h to be an exact numerical representation
								final float xx = a[targetParameter];

								// Get h to minimise roundoff error
								float h = h_; // ((xx == 0) ? 1 : xx) * h_;
								final float temp = xx + h;
								doNothing(temp);
								h = temp - xx;

								for (int x : testx)
									for (int y : testy)
									{
										int i = y * maxx + x;
										a[targetParameter] = xx;
										ff1.value(getVariables(indices, a), dyda, i);

										// Evaluate at (x+h) and (x-h)
										a[targetParameter] = xx + h;
										double value2 = ff1.value(getVariables(indices, a), dyda2, i);

										a[targetParameter] = xx - h;
										double value3 = ff1.value(getVariables(indices, a), dyda2, i);

										double gradient = (value2 - value3) / (2 * h);
										boolean ok = Math.signum(gradient) == Math.signum(dyda[gradientIndex]) ||
												Math.abs(gradient - dyda[gradientIndex]) < 0.1;
										//logf("[%s-%s]/2*%g : %g == %g\n", "" + value2, "" + value3, h, gradient,
										//		dyda[gradientIndex]);
										if (!ok)
											Assert.assertTrue(NAME[targetParameter] + ": " + gradient + " != " +
													dyda[gradientIndex], ok);
										ok = eqPerDatum.almostEqualComplement(gradient, dyda[gradientIndex]);
										if (ok)
											count++;
										total++;
									}
							}
		double p = (100.0 * count) / total;
		logf("Per Datum %s : %s = %d / %d (%.2f)\n", f1.getClass().getSimpleName(), NAME[targetParameter], count,
				total, p);
		Assert.assertTrue(NAME[targetParameter] + " fraction too low per datum: " + p, p > 90);
	}

	@Test
	public void fitFixedComputesGradient()
	{
		functionComputesGradient(GaussianFunctionFactory.FIT_FIXED);
	}

	@Test
	public void fitCircleComputesGradient()
	{
		functionComputesGradient(GaussianFunctionFactory.FIT_CIRCLE);
	}

	@Test
	public void fitFreeCircleComputesGradient()
	{
		functionComputesGradient(GaussianFunctionFactory.FIT_FREE_CIRCLE);
	}

	@Test
	public void fitEllipticalComputesGradient()
	{
		// The elliptical function gradient evaluation is worse 
		DoubleEquality tmp = eq;
		eq = eqPerDatum;
		functionComputesGradient(GaussianFunctionFactory.FIT_ELLIPTICAL);
		eq =tmp;
	}

	@Test
	public void fitNBFixedComputesGradient()
	{
		functionComputesGradient(GaussianFunctionFactory.FIT_NB_FIXED);
	}

	@Test
	public void fitNBCircleComputesGradient()
	{
		functionComputesGradient(GaussianFunctionFactory.FIT_NB_CIRCLE);
	}

	@Test
	public void fitNBFreeCircleComputesGradient()
	{
		functionComputesGradient(GaussianFunctionFactory.FIT_NB_FREE_CIRCLE);
	}

	@Test
	public void fitNBEllipticalComputesGradient()
	{
		// The elliptical function gradient evaluation is worse 
		DoubleEquality tmp = eq;
		eq = eqPerDatum;
		functionComputesGradient(GaussianFunctionFactory.FIT_NB_ELLIPTICAL);
		eq = tmp;
	}

	private void functionComputesGradient(int flags)
	{
		Gaussian2DFunction f1 = GaussianFunctionFactory.create2D(1, maxx, flags);
		// Setup
		testbackground = testbackground_;
		testamplitude1 = testamplitude1_;
		testangle1 = testangle1_;
		testcx1 = testcx1_;
		testcy1 = testcy1_;
		testw1 = testw1_;
		if (!f1.evaluatesBackground())
		{
			testbackground = new float[] { testbackground[0] };
		}
		if (!f1.evaluatesAmplitude())
		{
			testamplitude1 = new float[] { testamplitude1[0] };
		}
		if (!f1.evaluatesAngle())
		{
			testangle1 = new float[] { 0 };
		}
		// Position is always evaluated

		boolean noSecondWidth = false;
		if (!f1.evaluatesSD0())
		{
			// Just use 1 width
			testw1 = new float[][] { testw1[0] };
			// If no width 0 then assume we have no width 1 as well
			noSecondWidth = true;
		}
		else if (!f1.evaluatesSD1())
		{
			// No evaluation of second width needs only variation in width 0 so truncate 
			testw1 = Arrays.copyOf(testw1, 2);
			noSecondWidth = true;
		}
		if (noSecondWidth)
		{
			for (int i = 0; i < testw1.length; i++)
			{
				testw1[i][1] = testw1[i][0];
			}
		}

		double fraction = 90;
		if (f1.evaluatesBackground())
			functionComputesTargetGradient(f1, Gaussian2DFunction.BACKGROUND, fraction);
		if (f1.evaluatesAmplitude())
			functionComputesTargetGradient(f1, Gaussian2DFunction.AMPLITUDE, fraction);
		if (f1.evaluatesAngle())
			functionComputesTargetGradient(f1, Gaussian2DFunction.ANGLE, fraction);
		functionComputesTargetGradient(f1, Gaussian2DFunction.X_POSITION, fraction);
		functionComputesTargetGradient(f1, Gaussian2DFunction.Y_POSITION, fraction);
		if (f1.evaluatesSD0())
			functionComputesTargetGradient(f1, Gaussian2DFunction.X_SD, fraction);
		if (f1.evaluatesSD1())
			functionComputesTargetGradient(f1, Gaussian2DFunction.Y_SD, fraction);
	}

	private void functionComputesTargetGradient(Gaussian2DFunction f1, int targetParameter, double threshold)
	{
		int[] indices = f1.gradientIndices();
		int gradientIndex = findGradientIndex(f1, targetParameter);
		double[] dyda = new double[indices.length];
		double[] dyda2 = new double[indices.length];
		float[] a;

		PoissonLikelihoodFunction ff1;

		int n = maxx * maxx;
		int count = 0, total = 0;

		for (float background : testbackground)
			for (float amplitude1 : testamplitude1)
				for (float angle1 : testangle1)
					for (float cx1 : testcx1)
						for (float cy1 : testcy1)
							for (float[] w1 : testw1)
							{
								a = createParameters(background, amplitude1, angle1, cx1, cy1, w1[0], w1[1]);

								// Create y as a function we would want to move towards
								float[] a2 = createParameters(background, amplitude1, angle1, cx1, cy1, w1[0], w1[1]);
								a2[targetParameter] *= 1.3;
								f1.initialise(a2);
								float[] data = new float[maxx * maxx];
								for (int i = 0; i < n; i++)
									data[i] = f1.eval(i);

								ff1 = new PoissonLikelihoodFunction(f1, a, data, n);

								// Numerically solve gradient. 
								// Calculate the step size h to be an exact numerical representation
								final float xx = a[targetParameter];

								// Get h to minimise roundoff error
								float h = h_; // ((xx == 0) ? 1 : xx) * h_;
								final float temp = xx + h;
								doNothing(temp);
								h = temp - xx;

								ff1.value(getVariables(indices, a), dyda);

								// Evaluate at (x+h) and (x-h)
								a[targetParameter] = xx + h;
								double value2 = ff1.value(getVariables(indices, a), dyda2);

								a[targetParameter] = xx - h;
								double value3 = ff1.value(getVariables(indices, a), dyda2);

								double gradient = (value2 - value3) / (2 * h);
								boolean ok = Math.signum(gradient) == Math.signum(dyda[gradientIndex]) ||
										Math.abs(gradient - dyda[gradientIndex]) < 0.1;
								//logf("[%s-%s]/2*%g : %g == %g\n", "" + value2, "" + value3, h, gradient,
								//		dyda[gradientIndex]);
								if (!ok)
									Assert.assertTrue(NAME[targetParameter] + ": " + gradient + " != " +
											dyda[gradientIndex], ok);
								ok = eq.almostEqualComplement(gradient, dyda[gradientIndex]);
								if (ok)
									count++;
								total++;

							}
		double p = (100.0 * count) / total;
		logf("%s : %s = %d / %d (%.2f)\n", f1.getClass().getSimpleName(), NAME[targetParameter], count, total, p);
		Assert.assertTrue(NAME[targetParameter] + " fraction too low: " + p, p > threshold);
	}

	private double[] getVariables(int[] indices, float[] a)
	{
		double[] variables = new double[indices.length];
		for (int i = 0; i < indices.length; i++)
		{
			variables[i] = a[indices[i]];

		}
		return variables;
	}

	private int findGradientIndex(Gaussian2DFunction f, int targetParameter)
	{
		int i = f.findGradientIndex(targetParameter);
		Assert.assertTrue("Cannot find gradient index", i >= 0);
		return i;
	}

	private void doNothing(float f)
	{

	}

	float[] createParameters(float... args)
	{
		return args;
	}

	void log(String message)
	{
		System.out.println(message);
	}

	void logf(String format, Object... args)
	{
		System.out.printf(format, args);
	}
}
