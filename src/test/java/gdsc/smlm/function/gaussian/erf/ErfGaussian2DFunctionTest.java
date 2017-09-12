package gdsc.smlm.function.gaussian.erf;

import org.apache.commons.math3.util.Precision;
import org.ejml.data.DenseMatrix64F;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.ij.Utils;
import gdsc.core.test.BaseTimingTask;
import gdsc.core.test.TimingService;
import gdsc.core.utils.BitFlags;
import gdsc.core.utils.DoubleEquality;
import gdsc.core.utils.Statistics;
import gdsc.core.utils.TurboList;
import gdsc.smlm.TestSettings;
import gdsc.smlm.function.ExtendedGradient2Procedure;
import gdsc.smlm.function.Gradient1Procedure;
import gdsc.smlm.function.Gradient2Procedure;
import gdsc.smlm.function.ValueProcedure;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.function.gaussian.Gaussian2DFunctionTest;
import gdsc.smlm.function.gaussian.GaussianFunctionFactory;

public abstract class ErfGaussian2DFunctionTest extends Gaussian2DFunctionTest
{
	public ErfGaussian2DFunctionTest()
	{
		super();
		// The derivative check can be tighter with the ERF since it is a true integration
		h_ = 0.0001;
		eq3 = new DoubleEquality(5e-3, 1e-3); // For the Gaussian integral
	}

	@Test
	public void factoryDefaultsToErfGaussian2DFunction()
	{
		Gaussian2DFunction f, f2;

		int flags2 = BitFlags.unset(flags, GaussianFunctionFactory.FIT_ERF);

		if (this.f2 != null)
		{
			f = GaussianFunctionFactory.create2D(2, maxx, maxy, flags, zModel);
			f2 = GaussianFunctionFactory.create2D(2, maxx, maxy, flags2, zModel);
			Assert.assertTrue("Incorrect function2", f.getClass() == f2.getClass());
		}
		else
		{
			f = GaussianFunctionFactory.create2D(1, maxx, maxy, flags, zModel);
			f2 = GaussianFunctionFactory.create2D(1, maxx, maxy, flags2, zModel);
			Assert.assertTrue("Incorrect function1", f.getClass() == f2.getClass());
		}
	}

	@Test
	public void functionComputesSecondBackgroundGradient()
	{
		if (f1.evaluatesBackground())
			functionComputesSecondTargetGradient(Gaussian2DFunction.BACKGROUND);
	}

	@Test
	public void functionComputesSecondSignalGradient()
	{
		if (f1.evaluatesSignal())
			functionComputesSecondTargetGradient(Gaussian2DFunction.SIGNAL);
	}

	@Test
	public void functionComputesSecondXGradient()
	{
		functionComputesSecondTargetGradient(Gaussian2DFunction.X_POSITION);
	}

	@Test
	public void functionComputesSecondYGradient()
	{
		functionComputesSecondTargetGradient(Gaussian2DFunction.Y_POSITION);
	}

	@Test
	public void functionComputesSecondZGradient()
	{
		if (f1.evaluatesZ())
			functionComputesSecondTargetGradient(Gaussian2DFunction.Z_POSITION);
	}

	@Test
	public void functionComputesSecondXWidthGradient()
	{
		if (f1.evaluatesSD0())
			functionComputesSecondTargetGradient(Gaussian2DFunction.X_SD);
	}

	@Test
	public void functionComputesSecondYWidthGradient()
	{
		if (f1.evaluatesSD1())
			functionComputesSecondTargetGradient(Gaussian2DFunction.Y_SD);
	}

	@Test
	public void functionComputesSecondAngleGradient()
	{
		if (f1.evaluatesAngle())
			functionComputesSecondTargetGradient(Gaussian2DFunction.ANGLE);
	}

	private void functionComputesSecondTargetGradient(int targetParameter)
	{
		ErfGaussian2DFunction f1 = (ErfGaussian2DFunction) this.f1;
		int gradientIndex = findGradientIndex(f1, targetParameter);
		double[] dyda = new double[f1.getNumberOfGradients()];
		double[] dyda2 = new double[dyda.length];
		double[] a;

		// Test fitting of second derivatives 
		ErfGaussian2DFunction f1a = (ErfGaussian2DFunction) GaussianFunctionFactory.create2D(1, maxx, maxy, flags,
				zModel);
		ErfGaussian2DFunction f1b = (ErfGaussian2DFunction) GaussianFunctionFactory.create2D(1, maxx, maxy, flags,
				zModel);
		Statistics s = new Statistics();

		for (double background : testbackground)
			// Peak 1
			for (double signal1 : testsignal1)
				for (double cx1 : testcx1)
					for (double cy1 : testcy1)
						for (double cz1 : testcz1)
							for (double[] w1 : testw1)
								for (double angle1 : testangle1)
								{
									a = createParameters(background, signal1, cx1, cy1, cz1, w1[0], w1[1], angle1);
									f1.initialise2(a);

									// Numerically solve gradient. 
									// Calculate the step size h to be an exact numerical representation
									final double xx = a[targetParameter];

									// Get h to minimise roundoff error
									double h = Precision.representableDelta(xx, h_);

									// Evaluate at (x+h) and (x-h)
									a[targetParameter] = xx + h;
									f1a.initialise1(a.clone());

									a[targetParameter] = xx - h;
									f1b.initialise1(a.clone());

									for (int x : testx)
										for (int y : testy)
										{
											int i = y * maxx + x;
											f1a.eval(i, dyda);
											double value2 = dyda[gradientIndex];
											f1b.eval(i, dyda);
											double value3 = dyda[gradientIndex];
											f1.eval(i, dyda, dyda2);

											double gradient = (value2 - value3) / (2 * h);
											double error = DoubleEquality.relativeError(gradient, dyda2[gradientIndex]);
											s.add(error);
											Assert.assertTrue(gradient + " sign != " + dyda2[gradientIndex],
													(gradient * dyda2[gradientIndex]) >= 0);
											//System.out.printf("[%d,%d] %f == [%d] %f? (%g)\n", x, y, gradient,
											//		gradientIndex, dyda2[gradientIndex], error);
											Assert.assertTrue(gradient + " != " + dyda2[gradientIndex],
													eq.almostEqualRelativeOrAbsolute(gradient, dyda2[gradientIndex]));

										}
								}
		System.out.printf("functionComputesSecondTargetGradient %s %s (error %s +/- %s)\n",
				f1.getClass().getSimpleName(), Gaussian2DFunction.getName(targetParameter), Utils.rounded(s.getMean()),
				Utils.rounded(s.getStandardDeviation()));
	}

	@Test
	public void functionComputesSecondBackgroundGradientWith2Peaks()
	{
		org.junit.Assume.assumeNotNull(f2);
		if (f2.evaluatesBackground())
			functionComputesSecondTargetGradientWith2Peaks(Gaussian2DFunction.BACKGROUND);
	}

	@Test
	public void functionComputesSecondSignalGradientWith2Peaks()
	{
		org.junit.Assume.assumeNotNull(f2);
		if (f2.evaluatesSignal())
		{
			functionComputesSecondTargetGradientWith2Peaks(Gaussian2DFunction.SIGNAL);
			functionComputesSecondTargetGradientWith2Peaks(
					Gaussian2DFunction.SIGNAL + Gaussian2DFunction.PARAMETERS_PER_PEAK);
		}
	}

	@Test
	public void functionComputesSecondXGradientWith2Peaks()
	{
		org.junit.Assume.assumeNotNull(f2);
		functionComputesSecondTargetGradientWith2Peaks(Gaussian2DFunction.X_POSITION);
		functionComputesSecondTargetGradientWith2Peaks(
				Gaussian2DFunction.Y_POSITION + Gaussian2DFunction.PARAMETERS_PER_PEAK);
	}

	@Test
	public void functionComputesSecondYGradientWith2Peaks()
	{
		org.junit.Assume.assumeNotNull(f2);
		functionComputesSecondTargetGradientWith2Peaks(Gaussian2DFunction.Y_POSITION);
		functionComputesSecondTargetGradientWith2Peaks(
				Gaussian2DFunction.Y_POSITION + Gaussian2DFunction.PARAMETERS_PER_PEAK);
	}

	@Test
	public void functionComputesSecondZGradientWith2Peaks()
	{
		org.junit.Assume.assumeNotNull(f2);
		if (f2.evaluatesZ())
		{
			functionComputesSecondTargetGradientWith2Peaks(Gaussian2DFunction.Z_POSITION);
			functionComputesSecondTargetGradientWith2Peaks(
					Gaussian2DFunction.Z_POSITION + Gaussian2DFunction.PARAMETERS_PER_PEAK);
		}
	}

	@Test
	public void functionComputesSecondXWidthGradientWith2Peaks()
	{
		org.junit.Assume.assumeNotNull(f2);
		if (f2.evaluatesSD0())
		{
			functionComputesSecondTargetGradientWith2Peaks(Gaussian2DFunction.X_SD);
			functionComputesSecondTargetGradientWith2Peaks(
					Gaussian2DFunction.X_SD + Gaussian2DFunction.PARAMETERS_PER_PEAK);
		}
	}

	@Test
	public void functionComputesSecondYWidthGradientWith2Peaks()
	{
		org.junit.Assume.assumeNotNull(f2);
		if (f2.evaluatesSD1())
		{
			functionComputesSecondTargetGradientWith2Peaks(Gaussian2DFunction.Y_SD);
			functionComputesSecondTargetGradientWith2Peaks(
					Gaussian2DFunction.Y_SD + Gaussian2DFunction.PARAMETERS_PER_PEAK);
		}
	}

	@Test
	public void functionComputesSecondAngleGradientWith2Peaks()
	{
		org.junit.Assume.assumeNotNull(f2);
		if (f2.evaluatesAngle())
		{
			functionComputesSecondTargetGradientWith2Peaks(Gaussian2DFunction.ANGLE);
			functionComputesSecondTargetGradientWith2Peaks(
					Gaussian2DFunction.ANGLE + Gaussian2DFunction.PARAMETERS_PER_PEAK);
		}
	}

	private void functionComputesSecondTargetGradientWith2Peaks(int targetParameter)
	{
		ErfGaussian2DFunction f2 = (ErfGaussian2DFunction) this.f2;
		int gradientIndex = findGradientIndex(f2, targetParameter);
		double[] dyda = new double[f2.getNumberOfGradients()];
		double[] dyda2 = new double[dyda.length];
		double[] a;

		// Test fitting of second derivatives 
		ErfGaussian2DFunction f2a = (ErfGaussian2DFunction) GaussianFunctionFactory.create2D(2, maxx, maxy, flags,
				zModel);
		ErfGaussian2DFunction f2b = (ErfGaussian2DFunction) GaussianFunctionFactory.create2D(2, maxx, maxy, flags,
				zModel);
		Statistics s = new Statistics();

		for (double background : testbackground)
			// Peak 1
			for (double signal1 : testsignal1)
				for (double cx1 : testcx1)
					for (double cy1 : testcy1)
						for (double cz1 : testcz1)
							for (double[] w1 : testw1)
								for (double angle1 : testangle1)
									// Peak 2
									for (double signal2 : testsignal2)
										for (double cx2 : testcx2)
											for (double cy2 : testcy2)
												for (double cz2 : testcz2)
													for (double[] w2 : testw2)
														for (double angle2 : testangle2)
														{
															a = createParameters(background, signal1, cx1, cy1, cz1,
																	w1[0], w1[1], angle1, signal2, cx2, cy2, cz2, w2[0],
																	w2[1], angle2);
															f2.initialise2(a);

															// Numerically solve gradient. 
															// Calculate the step size h to be an exact numerical representation
															final double xx = a[targetParameter];

															// Get h to minimise roundoff error
															double h = Precision.representableDelta(xx, h_);

															// Evaluate at (x+h) and (x-h)
															a[targetParameter] = xx + h;
															f2a.initialise1(a.clone());

															a[targetParameter] = xx - h;
															f2b.initialise1(a.clone());

															for (int x : testx)
																for (int y : testy)
																{
																	int i = y * maxx + x;
																	f2a.eval(i, dyda);
																	double value2 = dyda[gradientIndex];
																	f2b.eval(i, dyda);
																	double value3 = dyda[gradientIndex];
																	f2.eval(i, dyda, dyda2);

																	double gradient = (value2 - value3) / (2 * h);
																	double error = DoubleEquality.relativeError(
																			gradient, dyda2[gradientIndex]);
																	s.add(error);
																	Assert.assertTrue(
																			gradient + " sign != " +
																					dyda2[gradientIndex],
																			(gradient * dyda2[gradientIndex]) >= 0);
																	//System.out.printf("[%d,%d] %f == [%d] %f? (%g)\n", x, y, gradient,
																	//		gradientIndex, dyda2[gradientIndex], error);
																	Assert.assertTrue(
																			gradient + " != " + dyda2[gradientIndex],
																			eq.almostEqualRelativeOrAbsolute(gradient,
																					dyda2[gradientIndex]));
																}
														}
		System.out.printf("functionComputesSecondTargetGradient %s [%d] %s (error %s +/- %s)\n",
				f2.getClass().getSimpleName(), Gaussian2DFunction.getPeak(targetParameter),
				Gaussian2DFunction.getName(targetParameter), Utils.rounded(s.getMean()),
				Utils.rounded(s.getStandardDeviation()));
	}

	private class FunctionTimingTask extends BaseTimingTask
	{
		Gaussian2DFunction f;
		ErfGaussian2DFunction f2;
		double[][] x;
		int order;
		final double[] dyda, d2yda2;
		final int n = f1.size();

		public FunctionTimingTask(Gaussian2DFunction f, double[][] x, int order)
		{
			super(f.getClass().getSimpleName() + " " + order + " eval");
			this.f = f;
			if (order == 2)
				f2 = (ErfGaussian2DFunction) f;
			this.x = x;
			this.order = order;
			dyda = new double[f.getNumberOfGradients()];
			d2yda2 = new double[f.getNumberOfGradients()];
		}

		public int getSize()
		{
			return 1;
		}

		public Object getData(int i)
		{
			return null;
		}

		public Object run(Object data)
		{
			double s = 0;
			f = f.copy();
			if (order == 0)
			{
				for (int i = 0; i < x.length; i++)
				{
					f.initialise0(x[i]);
					for (int j = 0; j < n; j++)
						s += f.eval(j);
				}
			}
			else if (order == 1)
			{
				for (int i = 0; i < x.length; i++)
				{
					f.initialise1(x[i]);
					for (int j = 0; j < n; j++)
						s += f.eval(j, dyda);
				}
			}
			else
			{
				for (int i = 0; i < x.length; i++)
				{
					f2.initialise2(x[i]);
					for (int j = 0; j < n; j++)
						s += f2.eval(j, dyda, d2yda2);
				}
			}
			return s;
		}
	}

	// Speed test verses equivalent Gaussian2DFunction
	@Test
	public void functionIsFasterThanEquivalentGaussian2DFunction()
	{
		int flags = this.flags & ~GaussianFunctionFactory.FIT_ERF;
		final Gaussian2DFunction gf = GaussianFunctionFactory.create2D(1, maxx, maxy, flags, zModel);

		boolean zDepth = (flags & GaussianFunctionFactory.FIT_Z) != 0;

		final TurboList<double[]> params = new TurboList<double[]>();
		final TurboList<double[]> params2 = new TurboList<double[]>();
		for (double background : testbackground)
			// Peak 1
			for (double signal1 : testsignal1)
				for (double cx1 : testcx1)
					for (double cy1 : testcy1)
						for (double cz1 : testcz1)
							for (double[] w1 : testw1)
								for (double angle1 : testangle1)
								{
									double[] a = createParameters(background, signal1, cx1, cy1, cz1, w1[0], w1[1],
											angle1);
									params.add(a);
									if (zDepth)
									{
										// Change to a standard free circular function
										a = a.clone();
										a[Gaussian2DFunction.X_SD] *= zModel.getSx(a[Gaussian2DFunction.Z_POSITION]);
										a[Gaussian2DFunction.Y_SD] *= zModel.getSy(a[Gaussian2DFunction.Z_POSITION]);
										a[Gaussian2DFunction.Z_POSITION] = 0;
										params2.add(a);
									}
								}
		double[][] x = params.toArray(new double[params.size()][]);
		double[][] x2 = (zDepth) ? params2.toArray(new double[params2.size()][]) : x;

		int runs = 10000 / x.length;
		TimingService ts = new TimingService(runs);
		ts.execute(new FunctionTimingTask(gf, x2, 1));
		ts.execute(new FunctionTimingTask(gf, x2, 0));
		ts.execute(new FunctionTimingTask(f1, x, 2));
		ts.execute(new FunctionTimingTask(f1, x, 1));
		ts.execute(new FunctionTimingTask(f1, x, 0));

		int size = ts.getSize();
		ts.repeat(size);
		ts.report();

		// Sometimes this fails, probably due to JVM optimisations, so skip for now
		if (!TestSettings.ASSERT_SPEED_TESTS)
			return;

		int n = ts.getSize() - 1;
		Assert.assertTrue("0 order", ts.get(n).getMean() < ts.get(n - 3).getMean());
		n--;
		Assert.assertTrue("1 order", ts.get(n).getMean() < ts.get(n - 3).getMean());
	}

	@Test
	public void functionComputesGradientForEach()
	{
		final ErfGaussian2DFunction f1 = (ErfGaussian2DFunction) this.f1;

		final int n = f1.size();
		double[] du_da = new double[f1.getNumberOfGradients()];
		double[] du_db = new double[f1.getNumberOfGradients()];
		final double[] d2u_da2 = new double[f1.getNumberOfGradients()];

		final double[] values = new double[n];
		final double[][] jacobian = new double[n][];
		final double[][] jacobian2 = new double[n][];

		for (double background : testbackground)
			// Peak 1
			for (double signal1 : testsignal1)
				for (double cx1 : testcx1)
					for (double cy1 : testcy1)
						for (double cz1 : testcz1)
							for (double[] w1 : testw1)
								for (double angle1 : testangle1)
								{
									double[] a = createParameters(background, signal1, cx1, cy1, cz1, w1[0], w1[1],
											angle1);
									f1.initialiseExtended2(a);

									// Compute single
									for (int i = 0; i < n; i++)
									{
										double o1 = f1.eval(i, du_da);
										double o2 = f1.eval(i, du_db, d2u_da2);
										Assert.assertEquals("Value", o1, o2, 1e-10);
										Assert.assertArrayEquals("Jacobian!=Jacobian", du_da, du_db, 1e-10);
										values[i] = o1;
										jacobian[i] = du_da.clone();
										jacobian2[i] = d2u_da2.clone();
									}

									// Use procedures
									f1.forEach(new ValueProcedure()
									{
										int i = 0;

										public void execute(double value)
										{
											Assert.assertEquals("Value ValueProcedure", values[i], value, 1e-10);
											i++;
										}
									});

									f1.forEach(new Gradient1Procedure()
									{
										int i = 0;

										public void execute(double value, double[] dy_da)
										{
											Assert.assertEquals("Value Gradient1Procedure", values[i], value, 1e-10);
											Assert.assertArrayEquals("du_da Gradient1Procedure", jacobian[i], dy_da,
													1e-10);
											i++;
										}
									});

									f1.forEach(new Gradient2Procedure()
									{
										int i = 0;

										public void execute(double value, double[] dy_da, double[] d2y_da2)
										{
											Assert.assertEquals("Value Gradient2Procedure", values[i], value, 1e-10);
											Assert.assertArrayEquals("du_da Gradient2Procedure", jacobian[i], dy_da,
													1e-10);
											Assert.assertArrayEquals("d2u_da2 Gradient2Procedure", jacobian2[i],
													d2y_da2, 1e-10);
											i++;
										}
									});

									f1.forEach(new ExtendedGradient2Procedure()
									{
										int i = 0;

										public void executeExtended(double value, double[] dy_da, double[] d2y_dadb)
										{
											Assert.assertEquals("Value ExtendedGradient2Procedure", values[i], value,
													1e-10);
											Assert.assertArrayEquals("du_da ExtendedGradient2Procedure", jacobian[i],
													dy_da, 1e-10);
											for (int j = 0, k = 0; j < d2u_da2.length; j++, k += d2u_da2.length + 1)
												d2u_da2[j] = d2y_dadb[k];
											Assert.assertArrayEquals("d2u_da2 Gradient2Procedure", jacobian2[i],
													d2u_da2, 1e-10);
											i++;
										}
									});
								}
	}

	@Test
	public void functionComputesExtendedGradientForEach()
	{
		final ErfGaussian2DFunction f1 = (ErfGaussian2DFunction) this.f1;

		final int nparams = f1.getNumberOfGradients();
		final int[] gradientIndices = f1.gradientIndices();
		final double[] du_da = new double[f1.getNumberOfGradients()];
		final double[] du_db = new double[f1.getNumberOfGradients()];

		final ErfGaussian2DFunction[] fHigh = new ErfGaussian2DFunction[nparams];
		final ErfGaussian2DFunction[] fLow = new ErfGaussian2DFunction[nparams];
		final double[] delta = new double[nparams];
		for (int j = 0; j < nparams; j++)
		{
			fHigh[j] = f1.copy();
			fLow[j] = f1.copy();
		}

		for (double background : testbackground)
			// Peak 1
			for (double signal1 : testsignal1)
				for (double cx1 : testcx1)
					for (double cy1 : testcy1)
						for (double cz1 : testcz1)
							for (double[] w1 : testw1)
								for (double angle1 : testangle1)
								{
									double[] a = createParameters(background, signal1, cx1, cy1, cz1, w1[0], w1[1],
											angle1);
									f1.initialiseExtended2(a);

									// Create a set of functions initialised +/- delta in each parameter
									for (int j = 0; j < nparams; j++)
									{
										int targetParameter = gradientIndices[j];
										// Numerically solve gradient. 
										// Calculate the step size h to be an exact numerical representation
										final double xx = a[targetParameter];

										// Get h to minimise roundoff error
										double h = Precision.representableDelta(xx, h_);

										// Evaluate at (x+h) and (x-h)
										a[targetParameter] = xx + h;
										fHigh[j].initialise1(a.clone());

										a[targetParameter] = xx - h;
										fLow[j].initialise1(a.clone());

										a[targetParameter] = xx;
										delta[j] = 2 * h;
									}

									f1.forEach(new ExtendedGradient2Procedure()
									{
										int i = -1;

										public void executeExtended(double value, double[] dy_da, double[] d2y_dadb)
										{
											i++;
											DenseMatrix64F m = DenseMatrix64F.wrap(nparams, nparams, d2y_dadb);
											for (int j = 0; j < nparams; j++)
											{
												// Evaluate the function +/- delta for parameter j
												fHigh[j].eval(i, du_da);
												fLow[j].eval(i, du_db);
												// Check the gradient with respect to parameter k
												for (int k = 0; k < nparams; k++)
												{
													double gradient = (du_da[k] - du_db[k]) / delta[j];
													boolean ok = eq.almostEqualRelativeOrAbsolute(gradient,
															m.get(j, k));
													if (!ok)
													{
														System.out.printf("%d [%d,%d] %f ?= %f\n", i, j, k, gradient,
																m.get(j, k));
														Assert.fail(String.format("%d [%d,%d] %f != %f", i, j, k,
																gradient, m.get(j, k)));
													}
												}
											}
										}
									});
								}
	}

	@Test
	public void functionComputesGradientForEachWith2Peaks()
	{
		org.junit.Assume.assumeNotNull(f2);
		final ErfGaussian2DFunction f2 = (ErfGaussian2DFunction) this.f2;

		final int n = f2.size();
		double[] du_da = new double[f2.getNumberOfGradients()];
		double[] du_db = new double[f2.getNumberOfGradients()];
		final double[] d2u_da2 = new double[f2.getNumberOfGradients()];

		final double[] values = new double[n];
		final double[][] jacobian = new double[n][];
		final double[][] jacobian2 = new double[n][];

		for (double background : testbackground)
			// Peak 1
			for (double signal1 : testsignal1)
				for (double cx1 : testcx1)
					for (double cy1 : testcy1)
						for (double cz1 : testcz1)
							for (double[] w1 : testw1)
								for (double angle1 : testangle1)
									// Peak 2
									for (double signal2 : testsignal2)
										for (double cx2 : testcx2)
											for (double cy2 : testcy2)
												for (double cz2 : testcz2)
													for (double[] w2 : testw2)
														for (double angle2 : testangle2)
														{
															double[] a = createParameters(background, signal1, cx1, cy1,
																	cz1, w1[0], w1[1], angle1, signal2, cx2, cy2, cz2,
																	w2[0], w2[1], angle2);
															f2.initialiseExtended2(a);

															// Compute single
															for (int i = 0; i < n; i++)
															{
																double o1 = f2.eval(i, du_da);
																double o2 = f2.eval(i, du_db, d2u_da2);
																Assert.assertEquals("Value", o1, o2, 1e-10);
																Assert.assertArrayEquals("Jacobian!=Jacobian", du_da,
																		du_db, 1e-10);
																values[i] = o1;
																jacobian[i] = du_da.clone();
																jacobian2[i] = d2u_da2.clone();
															}

															// Use procedures
															f2.forEach(new ValueProcedure()
															{
																int i = 0;

																public void execute(double value)
																{
																	Assert.assertEquals("Value ValueProcedure",
																			values[i], value, 1e-10);
																	i++;
																}
															});

															f2.forEach(new Gradient1Procedure()
															{
																int i = 0;

																public void execute(double value, double[] dy_da)
																{
																	Assert.assertEquals("Value Gradient1Procedure",
																			values[i], value, 1e-10);
																	Assert.assertArrayEquals("du_da Gradient1Procedure",
																			jacobian[i], dy_da, 1e-10);
																	i++;
																}
															});

															f2.forEach(new Gradient2Procedure()
															{
																int i = 0;

																public void execute(double value, double[] dy_da,
																		double[] d2y_da2)
																{
																	Assert.assertEquals("Value Gradient2Procedure",
																			values[i], value, 1e-10);
																	Assert.assertArrayEquals("du_da Gradient2Procedure",
																			jacobian[i], dy_da, 1e-10);
																	Assert.assertArrayEquals(
																			"d2u_da2 Gradient2Procedure", jacobian2[i],
																			d2y_da2, 1e-10);
																	i++;
																}
															});

															f2.forEach(new ExtendedGradient2Procedure()
															{
																int i = 0;

																public void executeExtended(double value,
																		double[] dy_da, double[] d2y_dadb)
																{
																	Assert.assertEquals(
																			"Value ExtendedGradient2Procedure",
																			values[i], value, 1e-10);
																	Assert.assertArrayEquals(
																			"du_da ExtendedGradient2Procedure",
																			jacobian[i], dy_da, 1e-10);
																	for (int j = 0, k = 0; j < d2u_da2.length; j++, k += d2u_da2.length +
																			1)
																		d2u_da2[j] = d2y_dadb[k];
																	Assert.assertArrayEquals(
																			"d2u_da2 Gradient2Procedure", jacobian2[i],
																			d2u_da2, 1e-10);
																	i++;
																}
															});
														}
	}

	@Test
	public void functionComputesExtendedGradientForEachWith2Peaks()
	{
		org.junit.Assume.assumeNotNull(f2);
		final ErfGaussian2DFunction f2 = (ErfGaussian2DFunction) this.f2;

		final int nparams = f2.getNumberOfGradients();
		final int[] gradientIndices = f2.gradientIndices();
		final double[] du_da = new double[f2.getNumberOfGradients()];
		final double[] du_db = new double[f2.getNumberOfGradients()];

		final ErfGaussian2DFunction[] fHigh = new ErfGaussian2DFunction[nparams];
		final ErfGaussian2DFunction[] fLow = new ErfGaussian2DFunction[nparams];
		final double[] delta = new double[nparams];
		for (int j = 0; j < nparams; j++)
		{
			fHigh[j] = f2.copy();
			fLow[j] = f2.copy();
		}

		for (double background : testbackground)
			// Peak 1
			for (double signal1 : testsignal1)
				for (double cx1 : testcx1)
					for (double cy1 : testcy1)
						for (double cz1 : testcz1)
							for (double[] w1 : testw1)
								for (double angle1 : testangle1)
									// Peak 2
									for (double signal2 : testsignal2)
										for (double cx2 : testcx2)
											for (double cy2 : testcy2)
												for (double cz2 : testcz2)
													for (double[] w2 : testw2)
														for (double angle2 : testangle2)
														{
															double[] a = createParameters(background, signal1, cx1, cy1,
																	cz1, w1[0], w1[1], angle1, signal2, cx2, cy2, cz2,
																	w2[0], w2[1], angle2);

															f2.initialiseExtended2(a);

															// Create a set of functions initialised +/- delta in each parameter
															for (int j = 0; j < nparams; j++)
															{
																int targetParameter = gradientIndices[j];
																// Numerically solve gradient. 
																// Calculate the step size h to be an exact numerical representation
																final double xx = a[targetParameter];

																// Get h to minimise roundoff error
																double h = Precision.representableDelta(xx, h_);

																// Evaluate at (x+h) and (x-h)
																a[targetParameter] = xx + h;
																fHigh[j].initialise1(a.clone());

																a[targetParameter] = xx - h;
																fLow[j].initialise1(a.clone());

																a[targetParameter] = xx;
																delta[j] = 2 * h;
															}

															f2.forEach(new ExtendedGradient2Procedure()
															{
																int i = -1;

																public void executeExtended(double value,
																		double[] dy_da, double[] d2y_dadb)
																{
																	i++;
																	//if (i!=f2.size()/2) return;
																	DenseMatrix64F m = DenseMatrix64F.wrap(nparams,
																			nparams, d2y_dadb);
																	for (int j = 0; j < nparams; j++)
																	{
																		// Evaluate the function +/- delta for parameter j
																		fHigh[j].eval(i, du_da);
																		fLow[j].eval(i, du_db);
																		// Check the gradient with respect to parameter k
																		for (int k = 0; k < nparams; k++)
																		{
																			double gradient = (du_da[k] - du_db[k]) /
																					delta[j];
																			boolean ok = eq
																					.almostEqualRelativeOrAbsolute(
																							gradient, m.get(j, k));
																			if (!ok)
																			{
																				System.out.printf(
																						"%d [%d,%d] %f ?= %f\n", i, j,
																						k, gradient, m.get(j, k));
																				Assert.fail(String.format(
																						"%d [%d,%d] %f != %f", i, j, k,
																						gradient, m.get(j, k)));
																			}
																		}
																	}
																}
															});
															//return;
														}
	}

	abstract class SimpleProcedure
	{
		ErfGaussian2DFunction f;
		double s = 0;

		SimpleProcedure(ErfGaussian2DFunction f)
		{
			this.f = f;
		}

		void reset()
		{
			s = 0;
		}

		void run(double[] a)
		{
			f = f.copy();
			initialise(a);
			forEach();
		}

		abstract void initialise(double[] a);

		abstract void forEach();
	}

	class Procedure0 extends SimpleProcedure implements ValueProcedure
	{
		Procedure0(ErfGaussian2DFunction f)
		{
			super(f);
		}

		@Override
		void initialise(double[] a)
		{
			f.initialise0(a);
		}

		@Override
		void forEach()
		{
			f.forEach(this);
		}

		public void execute(double value)
		{
			s += value;
		}
	}

	class Procedure1 extends SimpleProcedure implements Gradient1Procedure
	{
		Procedure1(ErfGaussian2DFunction f)
		{
			super(f);
		}

		@Override
		void initialise(double[] a)
		{
			f.initialise1(a);
		}

		@Override
		void forEach()
		{
			f.forEach(this);
		}

		public void execute(double value, double[] dy_da)
		{
			s += value;
		}
	}

	class Procedure2 extends SimpleProcedure implements Gradient2Procedure
	{
		Procedure2(ErfGaussian2DFunction f)
		{
			super(f);
		}

		@Override
		void initialise(double[] a)
		{
			f.initialise2(a);
		}

		@Override
		void forEach()
		{
			f.forEach(this);
		}

		public void execute(double value, double[] dy_da, double[] d2y_da2)
		{
			s += value;
		}
	}

	private class ForEachTimingTask extends BaseTimingTask
	{
		double[][] x;
		SimpleProcedure p;

		public ForEachTimingTask(ErfGaussian2DFunction f, double[][] x, int order)
		{
			super(f.getClass().getSimpleName() + " " + order + " forEach");
			this.x = x;
			if (order == 0)
			{
				p = new Procedure0(f);
			}
			else if (order == 1)
			{
				p = new Procedure1(f);
			}
			else
			{
				p = new Procedure2(f);
			}
		}

		public int getSize()
		{
			return 1;
		}

		public Object getData(int i)
		{
			return null;
		}

		public Object run(Object data)
		{
			p.reset();
			for (int i = 0; i < x.length; i++)
			{
				p.run(x[i]);
			}
			return p.s;
		}
	}

	// Speed test forEach verses equivalent eval() function calls
	@Test
	public void functionIsFasterUsingForEach()
	{
		final ErfGaussian2DFunction f1 = (ErfGaussian2DFunction) this.f1;

		final TurboList<double[]> params = new TurboList<double[]>();
		for (double background : testbackground)
			// Peak 1
			for (double signal1 : testsignal1)
				for (double cx1 : testcx1)
					for (double cy1 : testcy1)
						for (double cz1 : testcz1)
							for (double[] w1 : testw1)
								for (double angle1 : testangle1)
								{
									double[] a = createParameters(background, signal1, cx1, cy1, cz1, w1[0], w1[1],
											angle1);
									params.add(a);
								}
		double[][] x = params.toArray(new double[params.size()][]);

		int runs = 10000 / x.length;
		TimingService ts = new TimingService(runs);
		ts.execute(new FunctionTimingTask(f1, x, 2));
		ts.execute(new FunctionTimingTask(f1, x, 1));
		ts.execute(new FunctionTimingTask(f1, x, 0));
		ts.execute(new ForEachTimingTask(f1, x, 2));
		ts.execute(new ForEachTimingTask(f1, x, 1));
		ts.execute(new ForEachTimingTask(f1, x, 0));

		int size = ts.getSize();
		ts.repeat(size);
		ts.report();

		// Sometimes this fails when running all the tests, probably due to JVM optimisations
		// in the more used eval(...) functions, so skip for now
		if (!TestSettings.ASSERT_SPEED_TESTS)
			return;

		int n = ts.getSize() - 1;
		Assert.assertTrue("0 order", ts.get(n).getMean() < ts.get(n - 3).getMean());
		n--;
		Assert.assertTrue("1 order", ts.get(n).getMean() < ts.get(n - 3).getMean());
		n--;
		Assert.assertTrue("2 order", ts.get(n).getMean() < ts.get(n - 3).getMean());
	}
}
