package gdsc.smlm.function.gaussian;

import java.util.Arrays;

import gdsc.core.utils.DoubleEquality;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.function.gaussian.GaussianFunctionFactory;

import org.junit.Assert;
import org.junit.Test;

public abstract class Gaussian2DFunctionTest
{
	protected DoubleEquality eq = new DoubleEquality(2, 1e-3);
	protected DoubleEquality eq2 = new DoubleEquality(5, 1e-8);

	// Compute as per Numerical Recipes 5.7.
	// Approximate error accuracy in single precision: Ef
	// Step size for derivatives:
	// h ~ (Ef)^(1/3) * xc
	// xc is the characteristic scale over which x changes, assumed to be 1 (not x as per NR since x is close to zero)
	protected final double h_ = 0.01; //(double) (Math.pow(1e-3f, 1.0 / 3));

	protected int[] testx = new int[] { 4, 5, 6 };
	protected int[] testy = new int[] { 4, 5, 6 };
	protected double[] testbackground = new double[] { 0, 400 };

	protected double[] testamplitude1 = new double[] { 15, 55, 105 };
	protected double[] testshape1 = new double[] { (double) (Math.PI / 5), (double) (Math.PI / 3) };
	protected double[] testcx1 = new double[] { 4.9, 5.3 };
	protected double[] testcy1 = new double[] { 4.8, 5.2 };
	protected double[][] testw1 = new double[][] { { 1.1, 1.2 }, { 1.5, 1.2 }, { 1.1, 1.7 }, { 1.5, 1.7 }, };

	protected double[] testamplitude2 = new double[] { 20, 50 };
	protected double[] testshape2 = new double[] { (double) (Math.PI / 7), (double) (Math.PI / 11) };
	protected double[] testcx2 = new double[] { 4.8, 5.3 };
	protected double[] testcy2 = new double[] { 5.1, 4.9 };
	protected double[][] testw2 = new double[][] { { 1.2, 1.4 }, { 1.3, 1.4 }, { 1.2, 1.5 }, { 1.3, 1.5 }, };

	protected int maxx = 10;
	protected double background = 50;
	protected double shape = 0;
	protected double width = 5;
	protected Gaussian2DFunction f1;
	protected Gaussian2DFunction f2 = null;
	protected int flags;

	public Gaussian2DFunctionTest()
	{
		init();

		// Setup Tests
		if (!f1.evaluatesBackground())
		{
			testbackground = new double[] { testbackground[0] };
		}
		if (!f1.evaluatesSignal())
		{
			testamplitude1 = new double[] { testamplitude1[0] };
			testamplitude2 = new double[] { testamplitude2[0] };
		}
		if (!f1.evaluatesShape())
		{
			testshape1 = new double[] { 0 };
			testshape2 = new double[] { 0 };
		}
		// Position is always evaluated

		boolean noSecondWidth = false;
		if (!f1.evaluatesSD0())
		{
			// Just use 1 width
			testw1 = new double[][] { testw1[0] };
			testw2 = new double[][] { testw2[0] };
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

	private void checkGradientIndices(int npeaks, Gaussian2DFunction gf)
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
			if (gf.evaluatesSignal())
				Assert.assertEquals(gf.getName(i), i, gradientIndices[p++]);
			if (gf.evaluatesShape())
				Assert.assertEquals(gf.getName(i + 1), i + 1, gradientIndices[p++]);
			if (gf.evaluatesPosition())
			{
				Assert.assertEquals(gf.getName(i + 2), i + 2, gradientIndices[p++]);
				Assert.assertEquals(gf.getName(i + 3), i + 3, gradientIndices[p++]);
			}
			if (gf.evaluatesSD0())
				Assert.assertEquals(gf.getName(i + 4), i + 4, gradientIndices[p++]);
			if (gf.evaluatesSD1())
				Assert.assertEquals(gf.getName(i + 5), i + 5, gradientIndices[p++]);
		}
	}

	@Test
	public void factoryCreatesCorrectFunction()
	{
		Gaussian2DFunction f;

		if (f2 != null)
		{
			f = GaussianFunctionFactory.create2D(2, maxx, maxx, flags);
			Assert.assertTrue("Incorrect function2", f.getClass() == f2.getClass());
		}
		else
		{
			f = GaussianFunctionFactory.create2D(1, maxx, maxx, flags);
			Assert.assertTrue("Incorrect function1", f.getClass() == f1.getClass());
		}
	}

	@Test
	public void functionComputesTargetWithAndWithoutGradient()
	{
		double[] dyda = new double[f1.gradientIndices().length];
		double[] a;

		for (int x : testx)
			for (int y : testy)
				for (double background : testbackground)
					// Peak 1
					for (double amplitude1 : testamplitude1)
						for (double shape1 : testshape1)
							for (double cx1 : testcx1)
								for (double cy1 : testcy1)
									for (double[] w1 : testw1)
									{
										a = createParameters(background, amplitude1, shape1, cx1, cy1, w1[0], w1[1]);

										f1.initialise(a);
										double y1 = f1.eval(y * maxx + x, dyda);
										double y2 = f1.eval(y * maxx + x);

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
		if (f1.evaluatesSignal())
			functionComputesTargetGradient(Gaussian2DFunction.SIGNAL);
	}

	@Test
	public void functionComputesShapeGradient()
	{
		if (f1.evaluatesShape())
			functionComputesTargetGradient(Gaussian2DFunction.SHAPE);
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
		double[] dyda = new double[f1.gradientIndices().length];
		double[] dyda2 = new double[dyda.length];
		double[] a;

		Gaussian2DFunction f1a = GaussianFunctionFactory.create2D(1, maxx, maxx, flags);
		Gaussian2DFunction f1b = GaussianFunctionFactory.create2D(1, maxx, maxx, flags);

		for (double background : testbackground)
			// Peak 1
			for (double amplitude1 : testamplitude1)
				for (double shape1 : testshape1)
					for (double cx1 : testcx1)
						for (double cy1 : testcy1)
							for (double[] w1 : testw1)
							{
								a = createParameters(background, amplitude1, shape1, cx1, cy1, w1[0], w1[1]);
								f1.initialise(a);

								// Numerically solve gradient. 
								// Calculate the step size h to be an exact numerical representation
								final double xx = a[targetParameter];

								// Get h to minimise roundoff error
								double h = h_; //((xx == 0) ? 1 : xx) * h_;
								final double temp = xx + h;
								doNothing(temp);
								h = temp - xx;

								// Evaluate at (x+h) and (x-h)
								a = createParameters(background, amplitude1, shape1, cx1, cy1, w1[0], w1[1]);
								a[targetParameter] = xx + h;
								f1a.initialise(a);

								a = createParameters(background, amplitude1, shape1, cx1, cy1, w1[0], w1[1]);
								a[targetParameter] = xx - h;
								f1b.initialise(a);

								for (int x : testx)
									for (int y : testy)
									{
										int i = y * maxx + x;
										f1.eval(i, dyda);
										double value2 = f1a.eval(i, dyda2);
										double value3 = f1b.eval(i, dyda2);

										double gradient = (value2 - value3) / (2 * h);
										//System.out.printf("[%d,%d] %f == [%d] %f?\n", x, y, gradient, gradientIndex, dyda[gradientIndex]);
										Assert.assertTrue(gradient + " != " + dyda[gradientIndex],
												eq.almostEqualComplement(gradient, dyda[gradientIndex]));
									}
							}
	}

	protected int findGradientIndex(Gaussian2DFunction f, int targetParameter)
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

		double[] dyda = new double[f2.gradientIndices().length];
		double[] a;

		for (int x : testx)
			for (int y : testy)
				for (double background : testbackground)
					// Peak 1
					for (double amplitude1 : testamplitude1)
						for (double shape1 : testshape1)
							for (double cx1 : testcx1)
								for (double cy1 : testcy1)
									for (double[] w1 : testw1)
										// Peak 2
										for (double amplitude2 : testamplitude2)
											for (double shape2 : testshape2)
												for (double cx2 : testcx2)
													for (double cy2 : testcy2)
														for (double[] w2 : testw2)
														{
															a = createParameters(background, amplitude1, shape1, cx1,
																	cy1, w1[0], w1[1], amplitude2, shape2, cx2, cy2,
																	w2[0], w2[1]);

															f2.initialise(a);
															double y1 = f2.eval(y * maxx + x, dyda);
															double y2 = f2.eval(y * maxx + x);

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
			if (f2.evaluatesSignal())
			{
				functionComputesTargetGradientWith2Peaks(Gaussian2DFunction.SIGNAL);
				functionComputesTargetGradientWith2Peaks(Gaussian2DFunction.SIGNAL + 6);
			}
		}
	}

	@Test
	public void functionComputesShapeGradientWith2Peaks()
	{
		if (f2 != null)
		{
			if (f2.evaluatesShape())
			{
				functionComputesTargetGradientWith2Peaks(Gaussian2DFunction.SHAPE);
				functionComputesTargetGradientWith2Peaks(Gaussian2DFunction.SHAPE + 6);
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
		double[] dyda = new double[f2.gradientIndices().length];
		double[] dyda2 = new double[dyda.length];
		double[] a;

		Gaussian2DFunction f2a = GaussianFunctionFactory.create2D(2, maxx, maxx, flags);
		Gaussian2DFunction f2b = GaussianFunctionFactory.create2D(2, maxx, maxx, flags);

		for (double background : testbackground)
			// Peak 1
			for (double amplitude1 : testamplitude1)
				for (double shape1 : testshape1)
					for (double cx1 : testcx1)
						for (double cy1 : testcy1)
							for (double[] w1 : testw1)
								// Peak 2
								for (double amplitude2 : testamplitude2)
									for (double shape2 : testshape2)
										for (double cx2 : testcx2)
											for (double cy2 : testcy2)
												for (double[] w2 : testw2)
												{
													a = createParameters(background, amplitude1, shape1, cx1, cy1,
															w1[0], w1[1], amplitude2, shape2, cx2, cy2, w2[0], w2[1]);

													f2.initialise(a);

													// Numerically solve gradient. 
													// Calculate the step size h to be an exact numerical representation
													final double xx = a[targetParameter];

													// Get h to minimise roundoff error
													double h = h_; //((xx == 0) ? 1 : xx) * h_;
													double temp = xx + h;
													doNothing(temp);
													h = temp - xx;

													// Evaluate at (x+h) and (x-h)
													a = createParameters(background, amplitude1, shape1, cx1, cy1,
															w1[0], w1[1], amplitude2, shape2, cx2, cy2, w2[0], w2[1]);
													a[targetParameter] = xx + h;
													f2a.initialise(a);

													a = createParameters(background, amplitude1, shape1, cx1, cy1,
															w1[0], w1[1], amplitude2, shape2, cx2, cy2, w2[0], w2[1]);
													a[targetParameter] = xx - h;
													f2b.initialise(a);

													for (int x : testx)
														for (int y : testy)
														{
															int i= y * maxx + x;
															f2.eval(i, dyda);
															double value2 = f2a.eval(i, dyda2);
															double value3 = f2b.eval(i, dyda2);

															double gradient = (value2 - value3) / (2 * h);
															Assert.assertTrue(gradient + " != " + dyda[gradientIndex],
																	eq.almostEqualComplement(gradient,
																			dyda[gradientIndex]));
														}
												}
	}

	protected void doNothing(double f)
	{

	}

	@Test
	public void functionComputesGaussian()
	{
		double background = 0;
		int maxx = 30;

		Gaussian2DFunction f = GaussianFunctionFactory.create2D(1, maxx, maxx, flags);
		Gaussian2DFunction f2;
		if ((flags & GaussianFunctionFactory.FIT_ERF) == 0)
			f2 = GaussianFunctionFactory.create2D(1, maxx, maxx, GaussianFunctionFactory.FIT_ELLIPTICAL);
		else
			f2 = GaussianFunctionFactory.create2D(1, maxx, maxx, GaussianFunctionFactory.FIT_ERF_FREE_CIRCLE);

		for (double amplitude1 : testamplitude1)
			for (double shape1 : testshape1)
				for (double cx1 : new double[] { maxx / 2 + 0.373f })
					for (double cy1 : new double[] { maxx / 2 + 0.876f })
						for (double[] w1 : testw1)
						{
							double[] a = createParameters(background, amplitude1, shape1, cx1, cy1, w1[0], w1[1]);

							f.initialise(a);
							f2.initialise(a);
							double sum = 0;
							for (int index = maxx * maxx; index-- > 0;)
							{
								double r1 = f.eval(index);
								double r2 = f2.eval(index);
								//System.out.printf("%d,%d r1=%f\n", index%maxx, index/maxx, r1);
								sum += r1;
								final boolean ok = eq2.almostEqualComplement(r1, r2);
								if (!ok)
									Assert.assertTrue(
											String.format("%g != %g @ [%d,%d]", r1, r2, index / maxx, index % maxx),
											ok);
							}

							Assert.assertTrue(sum + " != " + amplitude1,
									eq.almostEqualComplement((double) sum, amplitude1));
						}
	}

	protected double[] createParameters(double... args)
	{
		return args;
	}

	protected void log(String message)
	{
		System.out.println(message);
	}

	protected void logf(String format, Object... args)
	{
		System.out.printf(format, args);
	}
}
