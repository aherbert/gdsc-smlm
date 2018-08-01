package uk.ac.sussex.gdsc.smlm.function.gaussian.erf;

import java.util.logging.Level;
import java.util.logging.Logger;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.BeforeAll;

import org.apache.commons.math3.util.Precision;
import org.ejml.data.DenseMatrix64F;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.Test;

import uk.ac.sussex.gdsc.core.ij.Utils;
import uk.ac.sussex.gdsc.core.utils.BitFlags;
import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.Statistics;
import uk.ac.sussex.gdsc.core.utils.TurboList;
import uk.ac.sussex.gdsc.smlm.function.ExtendedGradient2Procedure;
import uk.ac.sussex.gdsc.smlm.function.Gradient1Procedure;
import uk.ac.sussex.gdsc.smlm.function.Gradient2Procedure;
import uk.ac.sussex.gdsc.smlm.function.IntegralValueProcedure;
import uk.ac.sussex.gdsc.smlm.function.ValueProcedure;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunctionTest;
import uk.ac.sussex.gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import uk.ac.sussex.gdsc.test.BaseTimingTask;
import uk.ac.sussex.gdsc.test.LogLevel;
import uk.ac.sussex.gdsc.test.MessageProvider;
import uk.ac.sussex.gdsc.test.TestLog;
import uk.ac.sussex.gdsc.test.TestSettings;
import uk.ac.sussex.gdsc.test.TimingService;
import uk.ac.sussex.gdsc.test.junit5.ExtraAssertions;
import uk.ac.sussex.gdsc.test.junit5.ExtraAssumptions;
import uk.ac.sussex.gdsc.test.junit5.SpeedTag;

@SuppressWarnings({ "javadoc" })
public abstract class ErfGaussian2DFunctionTest extends Gaussian2DFunctionTest
{
	private static Logger logger;

	@BeforeAll
	public static void beforeAll()
	{
		logger = Logger.getLogger(ErfGaussian2DFunctionTest.class.getName());
	}

	@AfterAll
	public static void afterAll()
	{
		logger = null;
	}

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

		final int flags2 = BitFlags.unset(flags, GaussianFunctionFactory.FIT_ERF);

		if (this.f2 != null)
		{
			f = GaussianFunctionFactory.create2D(2, maxx, maxy, flags, zModel);
			f2 = GaussianFunctionFactory.create2D(2, maxx, maxy, flags2, zModel);
			Assertions.assertTrue(f.getClass() == f2.getClass(), "Incorrect function2");
		}
		else
		{
			f = GaussianFunctionFactory.create2D(1, maxx, maxy, flags, zModel);
			f2 = GaussianFunctionFactory.create2D(1, maxx, maxy, flags2, zModel);
			Assertions.assertTrue(f.getClass() == f2.getClass(), "Incorrect function1");
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
		final ErfGaussian2DFunction f1 = (ErfGaussian2DFunction) this.f1;
		final int gradientIndex = findGradientIndex(f1, targetParameter);
		final double[] dyda = new double[f1.getNumberOfGradients()];
		final double[] dyda2 = new double[dyda.length];
		double[] a;

		// Test fitting of second derivatives
		final ErfGaussian2DFunction f1a = (ErfGaussian2DFunction) GaussianFunctionFactory.create2D(1, maxx, maxy, flags,
				zModel);
		final ErfGaussian2DFunction f1b = (ErfGaussian2DFunction) GaussianFunctionFactory.create2D(1, maxx, maxy, flags,
				zModel);
		final Statistics s = new Statistics();

		for (final double background : testbackground)
			// Peak 1
			for (final double signal1 : testsignal1)
				for (final double cx1 : testcx1)
					for (final double cy1 : testcy1)
						for (final double cz1 : testcz1)
							for (final double[] w1 : testw1)
								for (final double angle1 : testangle1)
								{
									a = createParameters(background, signal1, cx1, cy1, cz1, w1[0], w1[1], angle1);
									f1.initialise2(a);

									// Numerically solve gradient.
									// Calculate the step size h to be an exact numerical representation
									final double xx = a[targetParameter];

									// Get h to minimise roundoff error
									final double h = Precision.representableDelta(xx, h_);

									// Evaluate at (x+h) and (x-h)
									a[targetParameter] = xx + h;
									f1a.initialise1(a.clone());

									a[targetParameter] = xx - h;
									f1b.initialise1(a.clone());

									for (final int x : testx)
										for (final int y : testy)
										{
											final int i = y * maxx + x;
											f1a.eval(i, dyda);
											final double value2 = dyda[gradientIndex];
											f1b.eval(i, dyda);
											final double value3 = dyda[gradientIndex];
											f1.eval(i, dyda, dyda2);

											final double gradient = (value2 - value3) / (2 * h);
											final double error = DoubleEquality.relativeError(gradient,
													dyda2[gradientIndex]);
											s.add(error);
											Assertions.assertTrue((gradient * dyda2[gradientIndex]) >= 0,
													() -> gradient + " sign != " + dyda2[gradientIndex]);
											//TestLog.fine(logger,"[%d,%d] %f == [%d] %f? (%g)\n", x, y, gradient,
											//		gradientIndex, dyda2[gradientIndex], error);
											Assertions.assertTrue(
													eq.almostEqualRelativeOrAbsolute(gradient, dyda2[gradientIndex]),
													() -> gradient + " != " + dyda2[gradientIndex]);

										}
								}
		logger.info(() -> {
			return String.format("functionComputesSecondTargetGradient %s %s (error %s +/- %s)\n",
					f1.getClass().getSimpleName(), Gaussian2DFunction.getName(targetParameter),
					Utils.rounded(s.getMean()), Utils.rounded(s.getStandardDeviation()));
		});
	}

	@Test
	public void functionComputesSecondBackgroundGradientWith2Peaks()
	{
		Assumptions.assumeTrue(null != f2);
		if (f2.evaluatesBackground())
			functionComputesSecondTargetGradientWith2Peaks(Gaussian2DFunction.BACKGROUND);
	}

	@Test
	public void functionComputesSecondSignalGradientWith2Peaks()
	{
		Assumptions.assumeTrue(null != f2);
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
		Assumptions.assumeTrue(null != f2);
		functionComputesSecondTargetGradientWith2Peaks(Gaussian2DFunction.X_POSITION);
		functionComputesSecondTargetGradientWith2Peaks(
				Gaussian2DFunction.Y_POSITION + Gaussian2DFunction.PARAMETERS_PER_PEAK);
	}

	@Test
	public void functionComputesSecondYGradientWith2Peaks()
	{
		Assumptions.assumeTrue(null != f2);
		functionComputesSecondTargetGradientWith2Peaks(Gaussian2DFunction.Y_POSITION);
		functionComputesSecondTargetGradientWith2Peaks(
				Gaussian2DFunction.Y_POSITION + Gaussian2DFunction.PARAMETERS_PER_PEAK);
	}

	@Test
	public void functionComputesSecondZGradientWith2Peaks()
	{
		Assumptions.assumeTrue(null != f2);
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
		Assumptions.assumeTrue(null != f2);
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
		Assumptions.assumeTrue(null != f2);
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
		Assumptions.assumeTrue(null != f2);
		if (f2.evaluatesAngle())
		{
			functionComputesSecondTargetGradientWith2Peaks(Gaussian2DFunction.ANGLE);
			functionComputesSecondTargetGradientWith2Peaks(
					Gaussian2DFunction.ANGLE + Gaussian2DFunction.PARAMETERS_PER_PEAK);
		}
	}

	private void functionComputesSecondTargetGradientWith2Peaks(int targetParameter)
	{
		final ErfGaussian2DFunction f2 = (ErfGaussian2DFunction) this.f2;
		final int gradientIndex = findGradientIndex(f2, targetParameter);
		final double[] dyda = new double[f2.getNumberOfGradients()];
		final double[] dyda2 = new double[dyda.length];
		double[] a;

		// Test fitting of second derivatives
		final ErfGaussian2DFunction f2a = (ErfGaussian2DFunction) GaussianFunctionFactory.create2D(2, maxx, maxy, flags,
				zModel);
		final ErfGaussian2DFunction f2b = (ErfGaussian2DFunction) GaussianFunctionFactory.create2D(2, maxx, maxy, flags,
				zModel);
		final Statistics s = new Statistics();

		for (final double background : testbackground)
			// Peak 1
			for (final double signal1 : testsignal1)
				for (final double cx1 : testcx1)
					for (final double cy1 : testcy1)
						for (final double cz1 : testcz1)
							for (final double[] w1 : testw1)
								for (final double angle1 : testangle1)
									// Peak 2
									for (final double signal2 : testsignal2)
										for (final double cx2 : testcx2)
											for (final double cy2 : testcy2)
												for (final double cz2 : testcz2)
													for (final double[] w2 : testw2)
														for (final double angle2 : testangle2)
														{
															a = createParameters(background, signal1, cx1, cy1, cz1,
																	w1[0], w1[1], angle1, signal2, cx2, cy2, cz2, w2[0],
																	w2[1], angle2);
															f2.initialise2(a);

															// Numerically solve gradient.
															// Calculate the step size h to be an exact numerical representation
															final double xx = a[targetParameter];

															// Get h to minimise roundoff error
															final double h = Precision.representableDelta(xx, h_);

															// Evaluate at (x+h) and (x-h)
															a[targetParameter] = xx + h;
															f2a.initialise1(a.clone());

															a[targetParameter] = xx - h;
															f2b.initialise1(a.clone());

															for (final int x : testx)
																for (final int y : testy)
																{
																	final int i = y * maxx + x;
																	f2a.eval(i, dyda);
																	final double value2 = dyda[gradientIndex];
																	f2b.eval(i, dyda);
																	final double value3 = dyda[gradientIndex];
																	f2.eval(i, dyda, dyda2);

																	final double gradient = (value2 - value3) / (2 * h);
																	final double error = DoubleEquality.relativeError(
																			gradient, dyda2[gradientIndex]);
																	s.add(error);
																	Assertions.assertTrue(

																			(gradient * dyda2[gradientIndex]) >= 0,
																			() -> gradient + " sign != " +
																					dyda2[gradientIndex]);
																	//TestLog.fine(logger,"[%d,%d] %f == [%d] %f? (%g)\n", x, y, gradient,
																	//		gradientIndex, dyda2[gradientIndex], error);
																	Assertions.assertTrue(

																			eq.almostEqualRelativeOrAbsolute(gradient,
																					dyda2[gradientIndex]),
																			() -> gradient + " != " +
																					dyda2[gradientIndex]);
																}
														}
		logger.info(() -> {
			return String.format("functionComputesSecondTargetGradient %s [%d] %s (error %s +/- %s)\n",
					f2.getClass().getSimpleName(), Gaussian2DFunction.getPeak(targetParameter),
					Gaussian2DFunction.getName(targetParameter), Utils.rounded(s.getMean()),
					Utils.rounded(s.getStandardDeviation()));
		});
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

		@Override
		public int getSize()
		{
			return 1;
		}

		@Override
		public Object getData(int i)
		{
			return null;
		}

		@Override
		public Object run(Object data)
		{
			double s = 0;
			f = f.copy();
			if (order == 0)
				for (int i = 0; i < x.length; i++)
				{
					f.initialise0(x[i]);
					for (int j = 0; j < n; j++)
						s += f.eval(j);
				}
			else if (order == 1)
				for (int i = 0; i < x.length; i++)
				{
					f.initialise1(x[i]);
					for (int j = 0; j < n; j++)
						s += f.eval(j, dyda);
				}
			else
				for (int i = 0; i < x.length; i++)
				{
					f2.initialise2(x[i]);
					for (int j = 0; j < n; j++)
						s += f2.eval(j, dyda, d2yda2);
				}
			return s;
		}
	}

	// Speed test verses equivalent Gaussian2DFunction
	@SpeedTag
	@Test
	public void functionIsFasterThanEquivalentGaussian2DFunction()
	{
		ExtraAssumptions.assumeSpeedTest();

		final int flags = this.flags & ~GaussianFunctionFactory.FIT_ERF;
		final Gaussian2DFunction gf = GaussianFunctionFactory.create2D(1, maxx, maxy, flags, zModel);

		final boolean zDepth = (flags & GaussianFunctionFactory.FIT_Z) != 0;

		final TurboList<double[]> params = new TurboList<>();
		final TurboList<double[]> params2 = new TurboList<>();
		for (final double background : testbackground)
			// Peak 1
			for (final double signal1 : testsignal1)
				for (final double cx1 : testcx1)
					for (final double cy1 : testcy1)
						for (final double cz1 : testcz1)
							for (final double[] w1 : testw1)
								for (final double angle1 : testangle1)
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
		final double[][] x = params.toArray(new double[params.size()][]);
		final double[][] x2 = (zDepth) ? params2.toArray(new double[params2.size()][]) : x;

		final int runs = 10000 / x.length;
		final TimingService ts = new TimingService(runs);
		ts.execute(new FunctionTimingTask(gf, x2, 1));
		ts.execute(new FunctionTimingTask(gf, x2, 0));
		ts.execute(new FunctionTimingTask(f1, x, 2));
		ts.execute(new FunctionTimingTask(f1, x, 1));
		ts.execute(new FunctionTimingTask(f1, x, 0));

		final int size = ts.getSize();
		ts.repeat(size);
		if (logger.isLoggable(Level.INFO))
			ts.report();

		for (int i = 1; i <= 2; i++)
		{
			final double t1 = ts.get(-i).getMean();
			final double t2 = ts.get(-i - 3).getMean();
			TestLog.logTestResult(logger, t1 < t2,
					"ERF function %d  %s  vs equivalent Gaussian2DFunction  %s : %.2fx\n", i - 1, t1, t2, t2 / t1);
		}
	}

	@Test
	public void functionComputesGradientForEach()
	{
		final ErfGaussian2DFunction f1 = (ErfGaussian2DFunction) this.f1;

		final int n = f1.size();
		final double[] du_da = new double[f1.getNumberOfGradients()];
		final double[] du_db = new double[f1.getNumberOfGradients()];
		final double[] d2u_da2 = new double[f1.getNumberOfGradients()];

		final double[] values = new double[n];
		final double[][] jacobian = new double[n][];
		final double[][] jacobian2 = new double[n][];

		for (final double background : testbackground)
			// Peak 1
			for (final double signal1 : testsignal1)
				for (final double cx1 : testcx1)
					for (final double cy1 : testcy1)
						for (final double cz1 : testcz1)
							for (final double[] w1 : testw1)
								for (final double angle1 : testangle1)
								{
									final double[] a = createParameters(background, signal1, cx1, cy1, cz1, w1[0],
											w1[1], angle1);
									f1.initialiseExtended2(a);

									// Compute single
									for (int i = 0; i < n; i++)
									{
										final double o1 = f1.eval(i, du_da);
										final double o2 = f1.eval(i, du_db, d2u_da2);
										Assertions.assertEquals(o1, o2, 1e-10, "Value");
										Assertions.assertArrayEquals(du_da, du_db, 1e-10, "Jacobian!=Jacobian");
										values[i] = o1;
										jacobian[i] = du_da.clone();
										jacobian2[i] = d2u_da2.clone();
									}

									// Use procedures
									f1.forEach(new ValueProcedure()
									{
										int i = 0;

										@Override
										public void execute(double value)
										{
											Assertions.assertEquals(values[i], value, 1e-10, "Value ValueProcedure");
											i++;
										}
									});

									f1.forEach(new Gradient1Procedure()
									{
										int i = 0;

										@Override
										public void execute(double value, double[] dy_da)
										{
											Assertions.assertEquals(values[i], value, 1e-10,
													"Value Gradient1Procedure");
											Assertions.assertArrayEquals(jacobian[i], dy_da, 1e-10,
													"du_da Gradient1Procedure");
											i++;
										}
									});

									f1.forEach(new Gradient2Procedure()
									{
										int i = 0;

										@Override
										public void execute(double value, double[] dy_da, double[] d2y_da2)
										{
											Assertions.assertEquals(values[i], value, 1e-10,
													"Value Gradient2Procedure");
											Assertions.assertArrayEquals(jacobian[i], dy_da, 1e-10,
													"du_da Gradient2Procedure");
											Assertions.assertArrayEquals(jacobian2[i], d2y_da2, 1e-10,
													"d2u_da2 Gradient2Procedure");
											i++;
										}
									});

									f1.forEach(new ExtendedGradient2Procedure()
									{
										int i = 0;

										@Override
										public void executeExtended(double value, double[] dy_da, double[] d2y_dadb)
										{
											Assertions.assertEquals(values[i], value, 1e-10,
													"Value ExtendedGradient2Procedure");
											Assertions.assertArrayEquals(jacobian[i], dy_da, 1e-10,
													"du_da ExtendedGradient2Procedure");
											for (int j = 0, k = 0; j < d2u_da2.length; j++, k += d2u_da2.length + 1)
												d2u_da2[j] = d2y_dadb[k];
											Assertions.assertArrayEquals(jacobian2[i], d2u_da2, 1e-10,
													"d2u_da2 Gradient2Procedure");
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

		final ErfGaussian2DFunction[] fHigh = new ErfGaussian2DFunction[nparams];
		final ErfGaussian2DFunction[] fLow = new ErfGaussian2DFunction[nparams];
		final double[] delta = new double[nparams];
		for (int j = 0; j < nparams; j++)
		{
			fHigh[j] = f1.copy();
			fLow[j] = f1.copy();
		}

		for (final double background : testbackground)
			// Peak 1
			for (final double signal1 : testsignal1)
				for (final double cx1 : testcx1)
					for (final double cy1 : testcy1)
						for (final double cz1 : testcz1)
							for (final double[] w1 : testw1)
								for (final double angle1 : testangle1)
								{
									final double[] a = createParameters(background, signal1, cx1, cy1, cz1, w1[0],
											w1[1], angle1);
									f1.initialiseExtended2(a);

									// Create a set of functions initialised +/- delta in each parameter
									for (int j = 0; j < nparams; j++)
									{
										final int targetParameter = gradientIndices[j];
										// Numerically solve gradient.
										// Calculate the step size h to be an exact numerical representation
										final double xx = a[targetParameter];

										// Get h to minimise roundoff error
										final double h = Precision.representableDelta(xx, h_);

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
										final double[] du_da = new double[f1.getNumberOfGradients()];
										final double[] du_db = new double[f1.getNumberOfGradients()];

										@Override
										public void executeExtended(double value, double[] dy_da, double[] d2y_dadb)
										{
											i++;
											final DenseMatrix64F m = DenseMatrix64F.wrap(nparams, nparams, d2y_dadb);
											for (int j = 0; j < nparams; j++)
											{
												// Evaluate the function +/- delta for parameter j
												fHigh[j].eval(i, du_da);
												fLow[j].eval(i, du_db);
												// Check the gradient with respect to parameter k
												for (int k = 0; k < nparams; k++)
												{
													final double gradient = (du_da[k] - du_db[k]) / delta[j];
													final boolean ok = eq.almostEqualRelativeOrAbsolute(gradient,
															m.get(j, k));
													if (!ok)
													{
														TestLog.info(logger, "%d [%d,%d] %f ?= %f\n", i, j, k, gradient,
																m.get(j, k));
														ExtraAssertions.fail("%d [%d,%d] %f != %f", i, j, k, gradient,
																m.get(j, k));
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
		Assumptions.assumeTrue(null != f2);
		final ErfGaussian2DFunction f2 = (ErfGaussian2DFunction) this.f2;

		final int n = f2.size();
		final double[] du_da = new double[f2.getNumberOfGradients()];
		final double[] du_db = new double[f2.getNumberOfGradients()];
		final double[] d2u_da2 = new double[f2.getNumberOfGradients()];

		final double[] values = new double[n];
		final double[][] jacobian = new double[n][];
		final double[][] jacobian2 = new double[n][];

		for (final double background : testbackground)
			// Peak 1
			for (final double signal1 : testsignal1)
				for (final double cx1 : testcx1)
					for (final double cy1 : testcy1)
						for (final double cz1 : testcz1)
							for (final double[] w1 : testw1)
								for (final double angle1 : testangle1)
									// Peak 2
									for (final double signal2 : testsignal2)
										for (final double cx2 : testcx2)
											for (final double cy2 : testcy2)
												for (final double cz2 : testcz2)
													for (final double[] w2 : testw2)
														for (final double angle2 : testangle2)
														{
															final double[] a = createParameters(background, signal1,
																	cx1, cy1, cz1, w1[0], w1[1], angle1, signal2, cx2,
																	cy2, cz2, w2[0], w2[1], angle2);
															f2.initialiseExtended2(a);

															// Compute single
															for (int i = 0; i < n; i++)
															{
																final double o1 = f2.eval(i, du_da);
																final double o2 = f2.eval(i, du_db, d2u_da2);
																Assertions.assertEquals(o1, o2, 1e-10, "Value");
																Assertions.assertArrayEquals(du_da, du_db, 1e-10,
																		"Jacobian!=Jacobian");
																values[i] = o1;
																jacobian[i] = du_da.clone();
																jacobian2[i] = d2u_da2.clone();
															}

															// Use procedures
															f2.forEach(new ValueProcedure()
															{
																int i = 0;

																@Override
																public void execute(double value)
																{
																	Assertions.assertEquals(values[i], value, 1e-10,
																			"Value ValueProcedure");
																	i++;
																}
															});

															f2.forEach(new Gradient1Procedure()
															{
																int i = 0;

																@Override
																public void execute(double value, double[] dy_da)
																{
																	Assertions.assertEquals(values[i], value, 1e-10,
																			"Value Gradient1Procedure");
																	Assertions.assertArrayEquals(jacobian[i], dy_da,
																			1e-10, "du_da Gradient1Procedure");
																	i++;
																}
															});

															f2.forEach(new Gradient2Procedure()
															{
																int i = 0;

																@Override
																public void execute(double value, double[] dy_da,
																		double[] d2y_da2)
																{
																	Assertions.assertEquals(values[i], value, 1e-10,
																			"Value Gradient2Procedure");
																	Assertions.assertArrayEquals(jacobian[i], dy_da,
																			1e-10, "du_da Gradient2Procedure");
																	Assertions.assertArrayEquals(jacobian2[i], d2y_da2,
																			1e-10, "d2u_da2 Gradient2Procedure");
																	i++;
																}
															});

															f2.forEach(new ExtendedGradient2Procedure()
															{
																int i = 0;

																@Override
																public void executeExtended(double value,
																		double[] dy_da, double[] d2y_dadb)
																{
																	Assertions.assertEquals(

																			values[i], value, 1e-10,
																			"Value ExtendedGradient2Procedure");
																	Assertions.assertArrayEquals(

																			jacobian[i], dy_da, 1e-10,
																			"du_da ExtendedGradient2Procedure");
																	for (int j = 0, k = 0; j < d2u_da2.length; j++, k += d2u_da2.length +
																			1)
																		d2u_da2[j] = d2y_dadb[k];
																	Assertions.assertArrayEquals(jacobian2[i], d2u_da2,
																			1e-10, "d2u_da2 Gradient2Procedure");
																	i++;
																}
															});
														}
	}

	@Test
	public void functionComputesExtendedGradientForEachWith2Peaks()
	{
		Assumptions.assumeTrue(null != f2);
		final ErfGaussian2DFunction f2 = (ErfGaussian2DFunction) this.f2;

		final int nparams = f2.getNumberOfGradients();
		final int[] gradientIndices = f2.gradientIndices();

		final ErfGaussian2DFunction[] fHigh = new ErfGaussian2DFunction[nparams];
		final ErfGaussian2DFunction[] fLow = new ErfGaussian2DFunction[nparams];
		final double[] delta = new double[nparams];
		for (int j = 0; j < nparams; j++)
		{
			fHigh[j] = f2.copy();
			fLow[j] = f2.copy();
		}

		for (final double background : testbackground)
			// Peak 1
			for (final double signal1 : testsignal1)
				for (final double cx1 : testcx1)
					for (final double cy1 : testcy1)
						for (final double cz1 : testcz1)
							for (final double[] w1 : testw1)
								for (final double angle1 : testangle1)
									// Peak 2
									for (final double signal2 : testsignal2)
										for (final double cx2 : testcx2)
											for (final double cy2 : testcy2)
												for (final double cz2 : testcz2)
													for (final double[] w2 : testw2)
														for (final double angle2 : testangle2)
														{
															final double[] a = createParameters(background, signal1,
																	cx1, cy1, cz1, w1[0], w1[1], angle1, signal2, cx2,
																	cy2, cz2, w2[0], w2[1], angle2);

															f2.initialiseExtended2(a);

															// Create a set of functions initialised +/- delta in each parameter
															for (int j = 0; j < nparams; j++)
															{
																final int targetParameter = gradientIndices[j];
																// Numerically solve gradient.
																// Calculate the step size h to be an exact numerical representation
																final double xx = a[targetParameter];

																// Get h to minimise roundoff error
																final double h = Precision.representableDelta(xx, h_);

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
																final double[] du_da = new double[f2
																		.getNumberOfGradients()];
																final double[] du_db = new double[f2
																		.getNumberOfGradients()];

																@Override
																public void executeExtended(double value,
																		double[] dy_da, double[] d2y_dadb)
																{
																	i++;
																	//if (i!=f2.size()/2) return;
																	final DenseMatrix64F m = DenseMatrix64F
																			.wrap(nparams, nparams, d2y_dadb);
																	for (int j = 0; j < nparams; j++)
																	{
																		// Evaluate the function +/- delta for parameter j
																		fHigh[j].eval(i, du_da);
																		fLow[j].eval(i, du_db);
																		// Check the gradient with respect to parameter k
																		for (int k = 0; k < nparams; k++)
																		{
																			final double gradient = (du_da[k] -
																					du_db[k]) / delta[j];
																			final boolean ok = eq
																					.almostEqualRelativeOrAbsolute(
																							gradient, m.get(j, k));
																			if (!ok)
																			{
																				TestLog.info(logger,
																						"%d [%d,%d] %f ?= %f\n", i, j,
																						k, gradient, m.get(j, k));
																				ExtraAssertions.fail(
																						"%d [%d,%d] %f != %f", i, j, k,
																						gradient, m.get(j, k));
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

		@Override
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

		@Override
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

		@Override
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
				p = new Procedure0(f);
			else if (order == 1)
				p = new Procedure1(f);
			else
				p = new Procedure2(f);
		}

		@Override
		public int getSize()
		{
			return 1;
		}

		@Override
		public Object getData(int i)
		{
			return null;
		}

		@Override
		public Object run(Object data)
		{
			p.reset();
			for (int i = 0; i < x.length; i++)
				p.run(x[i]);
			return p.s;
		}
	}

	// Speed test forEach verses equivalent eval() function calls
	@SpeedTag
	@Test
	public void functionIsFasterUsingForEach()
	{
		ExtraAssumptions.assumeSpeedTest();

		final ErfGaussian2DFunction f1 = (ErfGaussian2DFunction) this.f1;

		final TurboList<double[]> params = new TurboList<>();
		for (final double background : testbackground)
			// Peak 1
			for (final double signal1 : testsignal1)
				for (final double cx1 : testcx1)
					for (final double cy1 : testcy1)
						for (final double cz1 : testcz1)
							for (final double[] w1 : testw1)
								for (final double angle1 : testangle1)
								{
									final double[] a = createParameters(background, signal1, cx1, cy1, cz1, w1[0],
											w1[1], angle1);
									params.add(a);
								}
		final double[][] x = params.toArray(new double[params.size()][]);

		final int runs = 10000 / x.length;
		final TimingService ts = new TimingService(runs);
		ts.execute(new FunctionTimingTask(f1, x, 2));
		ts.execute(new FunctionTimingTask(f1, x, 1));
		ts.execute(new FunctionTimingTask(f1, x, 0));
		ts.execute(new ForEachTimingTask(f1, x, 2));
		ts.execute(new ForEachTimingTask(f1, x, 1));
		ts.execute(new ForEachTimingTask(f1, x, 0));

		final int size = ts.getSize();
		ts.repeat(size);
		if (logger.isLoggable(Level.INFO))
			ts.report();

		for (int i = 1; i <= 3; i++)
		{
			final double t1 = ts.get(-i).getMean();
			final double t2 = ts.get(-i - 3).getMean();
			TestLog.logTestResult(logger, t1 < t2, "forEach %d  order  %s  vs eval(int)  %s : %.2fx\n", i - 1, t1, t2,
					t2 / t1);
		}
	}

	@Test
	public void functionCanComputeIntegral()
	{
		double[] a;
		for (final double background : testbackground)
			// Peak 1
			for (final double signal1 : testsignal1)
				for (final double cx1 : testcx1)
					for (final double cy1 : testcy1)
						for (final double cz1 : testcz1)
							for (final double[] w1 : testw1)
								for (final double angle1 : testangle1)
								{
									a = createParameters(background, signal1, cx1, cy1, cz1, w1[0], w1[1], angle1);
									final double e = new IntegralValueProcedure().getIntegral(f1, a);
									final double o = f1.integral(a);
									ExtraAssertions.assertEqualsRelative(e, o, 1e-8);
								}
	}

	@Test
	public void computeIntegralIsFaster()
	{
		ExtraAssumptions.assumeMediumComplexity();

		final TurboList<double[]> p = new TurboList<>();
		for (final double background : testbackground)
			// Peak 1
			for (final double signal1 : testsignal1)
				for (final double cx1 : testcx1)
					for (final double cy1 : testcy1)
						for (final double cz1 : testcz1)
							for (final double[] w1 : testw1)
								for (final double angle1 : testangle1)
									p.add(createParameters(background, signal1, cx1, cy1, cz1, w1[0], w1[1], angle1));
		final int n = (int) Math.ceil(10000.0 / p.size());
		double s1 = 0, s2 = 0;
		long t1 = System.nanoTime();
		for (int i = n; i-- > 0;)
			for (int j = p.size(); j-- > 0;)
				s1 += new IntegralValueProcedure().getIntegral(f1, p.getf(j));
		long t2 = System.nanoTime();
		for (int i = n; i-- > 0;)
			for (int j = p.size(); j-- > 0;)
				s2 += f1.integral(p.getf(j));
		final long t3 = System.nanoTime();
		t1 = t2 - t1;
		t2 = t3 - t2;
		TestLog.info(logger, "computeIntegralIsFaster %s %d vs %d (%gx)\n", f1.getClass().getSimpleName(), t1, t2,
				(double) t1 / t2);
		ExtraAssertions.assertEqualsRelative(s1, s2, 1e-3);
		Assertions.assertTrue(t2 < t1);
	}

	@Test
	public void functionCanComputeIntegralWith2Peaks()
	{
		Assumptions.assumeTrue(null != f2);

		double[] a;
		for (final double background : testbackground)
			// Peak 1
			for (final double signal1 : testsignal1)
				for (final double cx1 : testcx1)
					for (final double cy1 : testcy1)
						for (final double cz1 : testcz1)
							for (final double[] w1 : testw1)
								for (final double angle1 : testangle1)
									// Peak 2
									for (final double signal2 : testsignal2)
										for (final double cx2 : testcx2)
											for (final double cy2 : testcy2)
												for (final double cz2 : testcz2)
													for (final double[] w2 : testw2)
														for (final double angle2 : testangle2)
														{
															a = createParameters(background, signal1, cx1, cy1, cz1,
																	w1[0], w1[1], angle1, signal2, cx2, cy2, cz2, w2[0],
																	w2[1], angle2);
															final double e = new IntegralValueProcedure()
																	.getIntegral(f2, a);
															final double o = f2.integral(a);
															ExtraAssertions.assertEqualsRelative(e, o, 1e-8);
														}
	}

	@Test
	public void computeIntegralIsFasterWith2Peaks()
	{
		ExtraAssumptions.assumeMediumComplexity();
		Assumptions.assumeTrue(null != f2);

		final TurboList<double[]> p = new TurboList<>();
		for (final double background : testbackground)
			// Peak 1
			for (final double signal1 : testsignal1)
				for (final double cx1 : testcx1)
					for (final double cy1 : testcy1)
						for (final double cz1 : testcz1)
							for (final double[] w1 : testw1)
								for (final double angle1 : testangle1)
									// Peak 2
									for (final double signal2 : testsignal2)
										for (final double cx2 : testcx2)
											for (final double cy2 : testcy2)
												for (final double cz2 : testcz2)
													for (final double[] w2 : testw2)
														for (final double angle2 : testangle2)
															p.add(createParameters(background, signal1, cx1, cy1, cz1,
																	w1[0], w1[1], angle1, signal2, cx2, cy2, cz2, w2[0],
																	w2[1], angle2));
		final int n = (int) Math.ceil(10000.0 / p.size());
		double s1 = 0, s2 = 0;
		long t1 = System.nanoTime();
		for (int i = n; i-- > 0;)
			for (int j = p.size(); j-- > 0;)
				s1 += new IntegralValueProcedure().getIntegral(f2, p.getf(j));
		long t2 = System.nanoTime();
		for (int i = n; i-- > 0;)
			for (int j = p.size(); j-- > 0;)
				s2 += f2.integral(p.getf(j));
		final long t3 = System.nanoTime();
		t1 = t2 - t1;
		t2 = t3 - t2;
		TestLog.info(logger, "computeIntegralIsFasterWith2Peaks %s %d vs %d (%gx)\n", f1.getClass().getSimpleName(), t1,
				t2, (double) t1 / t2);
		ExtraAssertions.assertEqualsRelative(s1, s2, 1e-3);
		Assertions.assertTrue(t2 < t1);
	}
}
