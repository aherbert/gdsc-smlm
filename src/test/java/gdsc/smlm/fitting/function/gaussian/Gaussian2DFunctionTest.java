package gdsc.smlm.fitting.function.gaussian;

import java.util.Arrays;

import gdsc.smlm.fitting.function.Gaussian2DFunction;
import gdsc.smlm.fitting.function.GaussianFunction;
import gdsc.smlm.fitting.function.GaussianFunctionFactory;
import gdsc.smlm.fitting.utils.FloatEquality;

import org.junit.Assert;
import org.junit.Test;

public abstract class Gaussian2DFunctionTest
{
	FloatEquality eq = new FloatEquality(2, 1e-3f);
	FloatEquality eq2 = new FloatEquality(5, 1e-8f);

	int[] testx = new int[] { 4, 5, 6 };
	int[] testy = new int[] { 4, 5, 6 };
	float[] testbackground = new float[] { 0, 400 };

	float[] testamplitude1 = new float[] { 15, 55, 105 };
	float[] testangle1 = new float[] { (float) (Math.PI / 5), (float) (Math.PI / 3) };
	float[] testcx1 = new float[] { 4.9f, 5.3f };
	float[] testcy1 = new float[] { 4.8f, 5.1f };
	float[][] testw1 = new float[][] { { 1.1f, 1.2f }, { 1.1f, 1.7f }, { 1.5f, 1.2f }, { 1.5f, 1.7f }, };

	float[] testamplitude2 = new float[] { 20, 50 };
	float[] testangle2 = new float[] { (float) (Math.PI / 7), (float) (Math.PI / 11) };
	float[] testcx2 = new float[] { 4.8f, 5.3f };
	float[] testcy2 = new float[] { 5.1f, 4.9f };
	float[][] testw2 = new float[][] { { 1.2f, 1.4f }, { 1.2f, 1.5f }, { 1.3f, 1.4f }, { 1.3f, 1.5f }, };

	int maxx = 10;
	float background = 50;
	float angle = 0;
	float width = 5;
	Gaussian2DFunction f1;
	Gaussian2DFunction f2 = null;
	int flags;

	public Gaussian2DFunctionTest()
	{
		init();

		// Setup Tests
		if (!f1.evaluatesBackground())
		{
			testbackground = new float[] { testbackground[0] };
		}
		if (!f1.evaluatesAmplitude())
		{
			testamplitude1 = new float[] { testamplitude1[0] };
			testamplitude2 = new float[] { testamplitude2[0] };
		}
		if (!f1.evaluatesAngle())
		{
			testangle1 = new float[] { 0 };
			testangle2 = new float[] { 0 };
		}
		// Position is always evaluated

		boolean noSecondWidth = false;
		if (!f1.evaluatesSD0())
		{
			// Just use 1 width
			testw1 = new float[][] { testw1[0] };
			testw2 = new float[][] { testw2[0] };
			// If no width 0 then assume we have no width 1 as well
			noSecondWidth = true;
		}
		else if (!f1.evaluatesSD1())
		{
			// No evaluation of second width needs only variation in width 0 so truncate 
			testw1 = Arrays.copyOf(testw1, 2);
			testw2 = Arrays.copyOf(testw2, 2);
			noSecondWidth = true;
		}
		if (noSecondWidth)
		{
			for (int i = 0; i < testw1.length; i++)
			{
				testw1[i][1] = testw1[i][0];
				testw2[i][1] = testw2[i][0];
			}
		}
	}

	/**
	 * Create the Gaussian2DFunction for 1 and 2 peaks. Creates the flags for the factory
	 */
	protected abstract void init();

	@Test
	public void functionCreatesCorrectGradientIndices()
	{
		checkGradientIndices(1, f1);
		checkGradientIndices(2, f2);
	}

	private void checkGradientIndices(int npeaks, GaussianFunction gf)
	{
		if (gf == null)
			return;

		int[] gradientIndices = gf.gradientIndices();
		logf("Function%d %s %s\n", npeaks, gf.getClass().getName(), Arrays.toString(gradientIndices));

		Assert.assertEquals("Incorrect number of peaks", gf.getNPeaks(), npeaks);

		int p = 0;
		if (gf.evaluatesBackground())
			Assert.assertEquals("Background", 0, gradientIndices[p++]);
		for (int peak = 1, i = 1; peak <= npeaks; peak++, i += 6)
		{
			if (gf.evaluatesAmplitude())
				Assert.assertEquals("Amplitude", i, gradientIndices[p++]);
			if (gf.evaluatesAngle())
				Assert.assertEquals("Angle", i + 1, gradientIndices[p++]);
			if (gf.evaluatesPosition())
			{
				Assert.assertEquals("Position0", i + 2, gradientIndices[p++]);
				Assert.assertEquals("Position1", i + 3, gradientIndices[p++]);
			}
			if (gf.evaluatesSD0())
				Assert.assertEquals("Width0", i + 4, gradientIndices[p++]);
			if (gf.evaluatesSD1())
				Assert.assertEquals("Width1", i + 5, gradientIndices[p++]);
		}
	}

	@Test
	public void factoryCreatesCorrectFunction()
	{
		GaussianFunction f;

		if (f2 != null)
		{
			f = GaussianFunctionFactory.create2D(2, maxx, flags);
			Assert.assertTrue("Incorrect function2", f.getClass() == f2.getClass());
		}
		else
		{
			f = GaussianFunctionFactory.create2D(1, maxx, flags);
			Assert.assertTrue("Incorrect function1", f.getClass() == f1.getClass());
		}
	}

	@Test
	public void functionComputesTargetWithAndWithoutGradient()
	{
		float[] dyda = new float[f1.gradientIndices().length];
		float[] a;

		for (int x : testx)
			for (int y : testy)
				for (float background : testbackground)
					// Peak 1
					for (float amplitude1 : testamplitude1)
						for (float angle1 : testangle1)
							for (float cx1 : testcx1)
								for (float cy1 : testcy1)
									for (float[] w1 : testw1)
									{
										a = createParameters(background, amplitude1, angle1, cx1, cy1, w1[0], w1[1]);

										f1.initialise(a);
										float y1 = f1.eval(y * maxx + x, dyda);
										float y2 = f1.eval(y * maxx + x);

										Assert.assertTrue(y1 + " != " + y2, eq2.almostEqualComplement(y1, y2));
									}
	}

	@Test
	public void functionComputesBackgroundGradient()
	{
		if (f1.evaluatesBackground())
			functionComputesTargetGradient(Gaussian2DFunction.BACKGROUND);
	}

	@Test
	public void functionComputesAmplitudeGradient()
	{
		if (f1.evaluatesAmplitude())
			functionComputesTargetGradient(Gaussian2DFunction.AMPLITUDE);
	}

	@Test
	public void functionComputesAngleGradient()
	{
		if (f1.evaluatesAngle())
			functionComputesTargetGradient(Gaussian2DFunction.ANGLE);
	}

	@Test
	public void functionComputesXGradient()
	{
		functionComputesTargetGradient(Gaussian2DFunction.X_POSITION);
	}

	@Test
	public void functionComputesYGradient()
	{
		functionComputesTargetGradient(Gaussian2DFunction.Y_POSITION);
	}

	@Test
	public void functionComputesXWidthGradient()
	{
		if (f1.evaluatesSD0())
			functionComputesTargetGradient(Gaussian2DFunction.X_SD);
	}

	@Test
	public void functionComputesYWidthGradient()
	{
		if (f1.evaluatesSD1())
			functionComputesTargetGradient(Gaussian2DFunction.Y_SD);
	}

	private void functionComputesTargetGradient(int targetParameter)
	{
		int gradientIndex = findGradientIndex(f1, targetParameter);
		float[] dyda = new float[f1.gradientIndices().length];
		float[] dyda2 = new float[dyda.length];
		float[] a;
		float delta = 0.05f;

		for (int x : testx)
			for (int y : testy)
				for (float background : testbackground)
					// Peak 1
					for (float amplitude1 : testamplitude1)
						for (float angle1 : testangle1)
							for (float cx1 : testcx1)
								for (float cy1 : testcy1)
									for (float[] w1 : testw1)
									{
										a = createParameters(background, amplitude1, angle1, cx1, cy1, w1[0], w1[1]);

										f1.initialise(a);
										f1.eval(y * maxx + x, dyda);

										// Numerically solve gradient
										a[targetParameter] += delta;

										f1.initialise(a);
										float value2 = f1.eval(y * maxx + x, dyda2);

										a[targetParameter] -= 2 * delta;

										f1.initialise(a);
										float value3 = f1.eval(y * maxx + x, dyda2);

										float gradient = (value2 - value3) / (2 * delta);
										Assert.assertTrue(gradient + " != " + dyda[gradientIndex],
												eq.almostEqualComplement(gradient, dyda[gradientIndex]));
									}
	}

	private int findGradientIndex(Gaussian2DFunction f, int targetParameter)
	{
		int i = f.findGradientIndex(targetParameter);
		Assert.assertTrue("Cannot find gradient index", i >= 0);
		return i;
	}

	@Test
	public void functionComputesTargetWithAndWithoutGradientWith2Peaks()
	{
		if (f2 == null)
			return;

		float[] dyda = new float[f2.gradientIndices().length];
		float[] a;

		for (int x : testx)
			for (int y : testy)
				for (float background : testbackground)
					// Peak 1
					for (float amplitude1 : testamplitude1)
						for (float angle1 : testangle1)
							for (float cx1 : testcx1)
								for (float cy1 : testcy1)
									for (float[] w1 : testw1)
										// Peak 2
										for (float amplitude2 : testamplitude2)
											for (float angle2 : testangle2)
												for (float cx2 : testcx2)
													for (float cy2 : testcy2)
														for (float[] w2 : testw2)
														{
															a = createParameters(background, amplitude1, angle1, cx1,
																	cy1, w1[0], w1[1], amplitude2, angle2, cx2, cy2,
																	w2[0], w2[1]);

															f2.initialise(a);
															float y1 = f2.eval(y * maxx + x, dyda);
															float y2 = f2.eval(y * maxx + x);

															Assert.assertTrue(y1 + " != " + y2,
																	eq2.almostEqualComplement(y1, y2));
														}
	}

	@Test
	public void functionComputesBackgroundGradientWith2Peaks()
	{
		if (f2 != null)
		{
			if (f2.evaluatesBackground())
				functionComputesTargetGradientWith2Peaks(Gaussian2DFunction.BACKGROUND);
		}
	}

	@Test
	public void functionComputesAmplitudeGradientWith2Peaks()
	{
		if (f2 != null)
		{
			if (f2.evaluatesAmplitude())
			{
				functionComputesTargetGradientWith2Peaks(Gaussian2DFunction.AMPLITUDE);
				functionComputesTargetGradientWith2Peaks(Gaussian2DFunction.AMPLITUDE + 6);
			}
		}
	}

	@Test
	public void functionComputesAngleGradientWith2Peaks()
	{
		if (f2 != null)
		{
			if (f2.evaluatesAngle())
			{
				functionComputesTargetGradientWith2Peaks(Gaussian2DFunction.ANGLE);
				functionComputesTargetGradientWith2Peaks(Gaussian2DFunction.ANGLE + 6);
			}
		}
	}

	@Test
	public void functionComputesXGradientWith2Peaks()
	{
		if (f2 != null)
		{
			functionComputesTargetGradientWith2Peaks(Gaussian2DFunction.X_POSITION);
			functionComputesTargetGradientWith2Peaks(Gaussian2DFunction.X_POSITION + 6);
		}
	}

	@Test
	public void functionComputesYGradientWith2Peaks()
	{
		if (f2 != null)
		{
			functionComputesTargetGradientWith2Peaks(Gaussian2DFunction.Y_POSITION);
			functionComputesTargetGradientWith2Peaks(Gaussian2DFunction.Y_POSITION + 6);
		}
	}

	@Test
	public void functionComputesXWidthGradientWith2Peaks()
	{
		if (f2 != null)
		{
			if (f2.evaluatesSD0())
			{
				functionComputesTargetGradientWith2Peaks(Gaussian2DFunction.X_SD);
				functionComputesTargetGradientWith2Peaks(Gaussian2DFunction.X_SD + 6);
			}
		}
	}

	@Test
	public void functionComputesYWidthGradientWith2Peaks()
	{
		if (f2 != null)
		{
			if (f2.evaluatesSD1())
			{
				functionComputesTargetGradientWith2Peaks(Gaussian2DFunction.Y_SD);
				functionComputesTargetGradientWith2Peaks(Gaussian2DFunction.Y_SD + 6);
			}
		}
	}

	private void functionComputesTargetGradientWith2Peaks(int targetParameter)
	{
		int gradientIndex = findGradientIndex(f2, targetParameter);
		float[] dyda = new float[f2.gradientIndices().length];
		float[] dyda2 = new float[dyda.length];
		float[] a;
		float delta = 0.05f;

		for (int x : testx)
			for (int y : testy)
				for (float background : testbackground)
					// Peak 1
					for (float amplitude1 : testamplitude1)
						for (float angle1 : testangle1)
							for (float cx1 : testcx1)
								for (float cy1 : testcy1)
									for (float[] w1 : testw1)
										// Peak 2
										for (float amplitude2 : testamplitude2)
											for (float angle2 : testangle2)
												for (float cx2 : testcx2)
													for (float cy2 : testcy2)
														for (float[] w2 : testw2)
														{
															a = createParameters(background, amplitude1, angle1, cx1,
																	cy1, w1[0], w1[1], amplitude2, angle2, cx2, cy2,
																	w2[0], w2[1]);

															f2.initialise(a);
															f2.eval(y * maxx + x, dyda);

															// Numerically solve gradient
															a[targetParameter] += delta;

															f2.initialise(a);
															float value2 = f2.eval(y * maxx + x, dyda2);

															a[targetParameter] -= 2 * delta;

															f2.initialise(a);
															float value3 = f2.eval(y * maxx + x, dyda2);

															float gradient = (value2 - value3) / (2 * delta);
															Assert.assertTrue(gradient + " != " + dyda[gradientIndex],
																	eq.almostEqualComplement(gradient,
																			dyda[gradientIndex]));
														}
	}

	@Test
	public void functionComputesGaussianIntegral()
	{
		float background = 0;
		int maxx = 30;

		Gaussian2DFunction f = f1;
		f.setMaxX(maxx);
		Gaussian2DFunction f2 = GaussianFunctionFactory.create2D(1, maxx, GaussianFunctionFactory.FIT_ELLIPTICAL);

		try
		{
			for (float amplitude1 : testamplitude1)
				for (float angle1 : testangle1)
					for (float cx1 : new float[] { maxx / 2 + 0.373f })
						for (float cy1 : new float[] { maxx / 2 + 0.876f })
							for (float[] w1 : testw1)
							{
								float[] a = createParameters(background, amplitude1, angle1, cx1, cy1, w1[0], w1[1]);

								f.initialise(a);
								f2.initialise(a);
								double sum = 0;
								for (int index = maxx * maxx; index-- > 0;)
								{
									float r1 = f.eval(index);
									float r2 = f2.eval(index);
									//System.out.printf("%d,%d r1=%f\n", index%maxx, index/maxx, r1);
									sum += r1;
									final boolean ok = eq2.almostEqualComplement(r1, r2);
									if (!ok)
										Assert.assertTrue(
												String.format("%g != %g @ [%d,%d]", r1, r2, index / maxx, index % maxx),
												ok);
								}

								float integral = (float) (amplitude1 * 2.0 * Math.PI * w1[0] * w1[1]);
								Assert.assertTrue(sum + " != " + integral,
										eq.almostEqualComplement((float) sum, integral));
							}
		}
		catch (AssertionError e)
		{
			throw e;
		}
		finally
		{
			// Reset the function width
			f.setMaxX(this.maxx);
		}
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
