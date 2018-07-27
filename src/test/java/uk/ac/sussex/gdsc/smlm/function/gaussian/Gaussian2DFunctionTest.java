/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */
package uk.ac.sussex.gdsc.smlm.function.gaussian;

import java.util.Arrays;

import org.apache.commons.math3.util.Precision;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.Test;import uk.ac.sussex.gdsc.test.junit5.SeededTest;import uk.ac.sussex.gdsc.test.junit5.RandomSeed;import uk.ac.sussex.gdsc.test.junit5.SpeedTag;

import uk.ac.sussex.gdsc.core.ij.Utils;
import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.Statistics;
import uk.ac.sussex.gdsc.test.LogLevel;
import uk.ac.sussex.gdsc.test.MessageProvider;
import uk.ac.sussex.gdsc.test.TestLog;
import uk.ac.sussex.gdsc.test.TestSettings;
import uk.ac.sussex.gdsc.test.junit5.ExtraAssertions;

@SuppressWarnings({ "javadoc" })
public abstract class Gaussian2DFunctionTest
{
	protected DoubleEquality eq = new DoubleEquality(1e-2, 1e-3);
	protected DoubleEquality eq2 = new DoubleEquality(1e-5, 1e-8);
	protected DoubleEquality eq3 = new DoubleEquality(1e-1, 1e-3); // For the Gaussian integral

	// Compute as per Numerical Recipes 5.7.
	// Approximate error accuracy in single precision: Ef
	// Step size for derivatives:
	// h ~ (Ef)^(1/3) * xc
	// xc is the characteristic scale over which x changes, assumed to be 1 (not x as per NR since x is close to zero)
	protected double h_ = 0.0001; //(double) (Math.pow(1e-3f, 1.0 / 3));

	protected int[] testx = new int[] { 4, 5, 6 };
	protected int[] testy = new int[] { 4, 5, 6 };
	protected double[] testbackground = new double[] { 1, 400 };

	protected double[] testsignal1 = new double[] { 15, 105 };
	protected double[] testcx1 = new double[] { 4.9, 5.3 };
	protected double[] testcy1 = new double[] { 4.8, 5.2 };
	protected double[] testcz1 = new double[] { -1.5, 1.0 };
	protected double[][] testw1 = new double[][] { { 1.1, 1.2 }, { 1.5, 1.2 }, { 1.1, 1.7 }, { 1.5, 1.7 }, };
	protected double[] testangle1 = new double[] { Math.PI / 5, Math.PI / 3 };

	protected double[] testsignal2 = new double[] { 20, 50 };
	protected double[] testcx2 = new double[] { 4.8, 5.3 };
	protected double[] testcy2 = new double[] { 5.1, 4.9 };
	protected double[] testcz2 = new double[] { -1.9, 0.7 };
	protected double[][] testw2 = new double[][] { { 1.2, 1.4 }, { 1.3, 1.4 }, { 1.2, 1.5 }, { 1.3, 1.5 }, };
	protected double[] testangle2 = new double[] { Math.PI / 7, Math.PI / 11 };

	// Different widths to test for non-square function evaluation
	protected int maxx = 8, maxy = 9;
	protected double background = 50;
	protected double angle = 0;
	protected double width = 5;
	protected Gaussian2DFunction f1;
	protected Gaussian2DFunction f2 = null;
	protected int flags;
	protected AstigmatismZModel zModel = null;

	public Gaussian2DFunctionTest()
	{
		init();

		// Setup Tests
		if (!f1.evaluatesBackground())
			testbackground = new double[] { testbackground[0] };
		if (!f1.evaluatesSignal())
		{
			testsignal1 = new double[] { testsignal1[0] };
			testsignal2 = new double[] { testsignal2[0] };
		}
		// XY Position is always evaluated
		if (!f1.evaluatesZ())
		{
			testcz1 = new double[] { 0 };
			testcz2 = new double[] { 0 };
		}

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
			for (int i = 0; i < testw1.length; i++)
			{
				testw1[i][1] = testw1[i][0];
				testw2[i][1] = testw2[i][0];
			}
		if (!f1.evaluatesAngle())
		{
			testangle1 = new double[] { 0 };
			testangle2 = new double[] { 0 };
		}

		postInit();
	}

	/**
	 * Create the Gaussian2DFunction for 1 and 2 peaks. Creates the flags for the factory
	 */
	protected abstract void init();

	protected void postInit()
	{
		// To be overridden
	}

	@Test
	public void functionCreatesCorrectGradientIndices()
	{
		checkGradientIndices(1, f1);
		checkGradientIndices(2, f2);
	}

	private static void checkGradientIndices(int npeaks, Gaussian2DFunction gf)
	{
		if (gf == null)
			return;

		final int[] gradientIndices = gf.gradientIndices();
		if (TestSettings.allow(LogLevel.INFO))
			TestLog.info("Function%d %s %s\n", npeaks, gf.getClass().getName(), Arrays.toString(gradientIndices));

		Assertions.assertEquals(gf.getNPeaks(), npeaks, "Incorrect number of peaks");

		int p = 0;
		if (gf.evaluatesBackground())
			Assertions.assertEquals(0, gradientIndices[p++], "Background");
		for (int peak = 1, i = 1; peak <= npeaks; peak++, i += Gaussian2DFunction.PARAMETERS_PER_PEAK)
		{
			final int ii = i;
			if (gf.evaluatesSignal())
				Assertions.assertEquals(i, gradientIndices[p++], () -> Gaussian2DFunction.getName(ii));
			if (gf.evaluatesPosition())
			{
				Assertions.assertEquals(i + 1, gradientIndices[p++], () -> Gaussian2DFunction.getName(ii + 1));
				Assertions.assertEquals(i + 2, gradientIndices[p++], () -> Gaussian2DFunction.getName(ii + 2));
			}
			if (gf.evaluatesZ())
				Assertions.assertEquals(i + 3, gradientIndices[p++], () -> Gaussian2DFunction.getName(ii + 3));
			if (gf.evaluatesSD0())
				Assertions.assertEquals(i + 4, gradientIndices[p++], () -> Gaussian2DFunction.getName(ii + 4));
			if (gf.evaluatesSD1())
				Assertions.assertEquals(i + 5, gradientIndices[p++], () -> Gaussian2DFunction.getName(ii + 5));
			if (gf.evaluatesAngle())
				Assertions.assertEquals(i + 6, gradientIndices[p++], () -> Gaussian2DFunction.getName(ii + 6));
		}
	}

	@Test
	public void factoryCreatesCorrectFunction()
	{
		Gaussian2DFunction f;

		if (f2 != null)
		{
			f = GaussianFunctionFactory.create2D(2, maxx, maxy, flags, zModel);
			Assertions.assertTrue(f.getClass() == f2.getClass(), "Incorrect function2");
		}
		else
		{
			f = GaussianFunctionFactory.create2D(1, maxx, maxy, flags, zModel);
			Assertions.assertTrue(f.getClass() == f1.getClass(), "Incorrect function1");
		}
	}

	@Test
	public void functionComputesTargetWithAndWithoutGradient()
	{
		final double[] dyda = new double[f1.gradientIndices().length];
		double[] a;

		boolean record = TestSettings.allow(LogLevel.INFO);

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

									f1.initialise(a);

									// Test the frozen version
									final int flags = GaussianFunctionFactory.freeze(this.flags, zModel, a);
									final Gaussian2DFunction f = GaussianFunctionFactory.create2D(1, f1.getMaxX(),
											f1.getMaxY(), flags, zModel);
									f.initialise(a);
									if (record)
									{
										record = false;
										TestLog.info("%s %d frozen to %s\n", f1.getClass().getSimpleName(), 1,
												f.getClass().getSimpleName());
									}

									for (final int x : testx)
										for (final int y : testy)
										{
											final int xx = y * maxx + x;
											final double y1 = f1.eval(xx, dyda);
											final double y2 = f1.eval(xx);

											Assertions.assertTrue(eq2.almostEqualRelativeOrAbsolute(y1, y2),
													() -> y1 + " != " + y2);

											final double y3 = f.eval(xx);

											Assertions.assertTrue(eq2.almostEqualRelativeOrAbsolute(y1, y3),
													() -> y1 + " != frozen " + y3);
										}
								}
	}

	@Test
	public void functionComputesBackgroundGradient()
	{
		Assumptions.assumeTrue(f1.evaluatesBackground());
		functionComputesTargetGradient(Gaussian2DFunction.BACKGROUND);
	}

	@Test
	public void functionComputesSignalGradient()
	{
		Assumptions.assumeTrue(f1.evaluatesSignal());
		functionComputesTargetGradient(Gaussian2DFunction.SIGNAL);
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
	public void functionComputesZGradient()
	{
		Assumptions.assumeTrue(f1.evaluatesZ());
		functionComputesTargetGradient(Gaussian2DFunction.Z_POSITION);
	}

	@Test
	public void functionComputesXWidthGradient()
	{
		Assumptions.assumeTrue(f1.evaluatesSD0());
		functionComputesTargetGradient(Gaussian2DFunction.X_SD);
	}

	@Test
	public void functionComputesYWidthGradient()
	{
		Assumptions.assumeTrue(f1.evaluatesSD1());
		functionComputesTargetGradient(Gaussian2DFunction.Y_SD);
	}

	@Test
	public void functionComputesAngleGradient()
	{
		Assumptions.assumeTrue(f1.evaluatesAngle());
		functionComputesTargetGradient(Gaussian2DFunction.ANGLE);
	}

	private void functionComputesTargetGradient(int targetParameter)
	{
		final int gradientIndex = findGradientIndex(f1, targetParameter);
		final double[] dyda = new double[f1.gradientIndices().length];
		final double[] dyda2 = new double[dyda.length];
		double[] a;

		final Gaussian2DFunction f1a = GaussianFunctionFactory.create2D(1, maxx, maxy, flags, zModel);
		final Gaussian2DFunction f1b = GaussianFunctionFactory.create2D(1, maxx, maxy, flags, zModel);
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
									f1.initialise(a);

									// Numerically solve gradient.
									// Calculate the step size h to be an exact numerical representation
									final double xx = a[targetParameter];

									// Get h to minimise roundoff error
									final double h = Precision.representableDelta(xx, h_);

									// Evaluate at (x+h) and (x-h)
									a[targetParameter] = xx + h;
									f1a.initialise(a.clone());

									a[targetParameter] = xx - h;
									f1b.initialise(a.clone());

									for (final int x : testx)
										for (final int y : testy)
										{
											final int i = y * maxx + x;
											f1.eval(i, dyda);
											final double value2 = f1a.eval(i, dyda2);
											final double value3 = f1b.eval(i, dyda2);

											final double gradient = (value2 - value3) / (2 * h);
											final double error = DoubleEquality.relativeError(gradient,
													dyda2[gradientIndex]);
											s.add(error);
											ExtraAssertions.assertTrue((gradient * dyda2[gradientIndex]) >= 0,
													"%s sign != %s", gradient, dyda2[gradientIndex]);
											//System.out.printf("[%d,%d] %f == [%d] %f? (%g)\n", x, y, gradient,
											//		gradientIndex, dyda2[gradientIndex], error);
											//System.out.printf("[%d,%d] %f == [%d] %f?\n", x, y, gradient, gradientIndex, dyda[gradientIndex]);
											ExtraAssertions.assertTrue(
													eq.almostEqualRelativeOrAbsolute(gradient, dyda[gradientIndex]),
													"%s != %s", gradient, dyda[gradientIndex]);

										}
								}
		TestLog.info(new MessageProvider()
		{
			@Override
			public String getMessage()
			{
				return String.format("functionComputesTargetGradient %s %s (error %s +/- %s)\n",
						f1.getClass().getSimpleName(), Gaussian2DFunction.getName(targetParameter),
						Utils.rounded(s.getMean()), Utils.rounded(s.getStandardDeviation()));
			}
		});
	}

	protected int findGradientIndex(Gaussian2DFunction f, int targetParameter)
	{
		final int i = f.findGradientIndex(targetParameter);
		Assertions.assertTrue(i >= 0, "Cannot find gradient index");
		return i;
	}

	@Test
	public void functionComputesTargetWithAndWithoutGradientWith2Peaks()
	{
		if (f2 == null)
			return;

		final double[] dyda = new double[f2.gradientIndices().length];
		double[] a;

		boolean record = TestSettings.allow(LogLevel.INFO);

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

															f2.initialise(a);

															// Test the frozen version
															final int flags = GaussianFunctionFactory.freeze(this.flags,
																	zModel, a);
															final Gaussian2DFunction f = GaussianFunctionFactory
																	.create2D(2, f2.getMaxX(), f2.getMaxY(), flags,
																			zModel);
															f.initialise(a);
															if (record)
															{
																record = false;
																TestLog.info("%s %d frozen to %s\n",
																		f2.getClass().getSimpleName(), 2,
																		f.getClass().getSimpleName());
															}

															for (final int x : testx)
																for (final int y : testy)
																{
																	final int xx = y * maxx + x;
																	final double y1 = f2.eval(xx, dyda);
																	final double y2 = f2.eval(xx);

																	Assertions.assertTrue(
																			eq2.almostEqualRelativeOrAbsolute(y1, y2),
																			() -> y1 + " != " + y2);

																	final double y3 = f.eval(xx);

																	Assertions.assertTrue(
																			eq2.almostEqualRelativeOrAbsolute(y1, y3),
																			() -> y1 + " != frozen " + y3);
																}
														}
	}

	@Test
	public void functionComputesBackgroundGradientWith2Peaks()
	{
		Assumptions.assumeTrue(null != f2);
		Assumptions.assumeTrue(f2.evaluatesBackground());
		functionComputesTargetGradientWith2Peaks(Gaussian2DFunction.BACKGROUND);
	}

	@Test
	public void functionComputesSignalGradientWith2Peaks()
	{
		Assumptions.assumeTrue(null != f2);
		Assumptions.assumeTrue(f2.evaluatesSignal());
		functionComputesTargetGradientWith2Peaks(Gaussian2DFunction.SIGNAL);
		functionComputesTargetGradientWith2Peaks(Gaussian2DFunction.SIGNAL + Gaussian2DFunction.PARAMETERS_PER_PEAK);
	}

	@Test
	public void functionComputesXGradientWith2Peaks()
	{
		Assumptions.assumeTrue(null != f2);
		functionComputesTargetGradientWith2Peaks(Gaussian2DFunction.X_POSITION);
		functionComputesTargetGradientWith2Peaks(
				Gaussian2DFunction.X_POSITION + Gaussian2DFunction.PARAMETERS_PER_PEAK);
	}

	@Test
	public void functionComputesYGradientWith2Peaks()
	{
		Assumptions.assumeTrue(null != f2);
		functionComputesTargetGradientWith2Peaks(Gaussian2DFunction.Y_POSITION);
		functionComputesTargetGradientWith2Peaks(
				Gaussian2DFunction.Y_POSITION + Gaussian2DFunction.PARAMETERS_PER_PEAK);
	}

	@Test
	public void functionComputesZGradientWith2Peaks()
	{
		Assumptions.assumeTrue(null != f2);
		Assumptions.assumeTrue(f2.evaluatesZ());
		functionComputesTargetGradientWith2Peaks(Gaussian2DFunction.Z_POSITION);
		functionComputesTargetGradientWith2Peaks(
				Gaussian2DFunction.Z_POSITION + Gaussian2DFunction.PARAMETERS_PER_PEAK);
	}

	@Test
	public void functionComputesXWidthGradientWith2Peaks()
	{
		Assumptions.assumeTrue(null != f2);
		Assumptions.assumeTrue(f2.evaluatesSD0());
		functionComputesTargetGradientWith2Peaks(Gaussian2DFunction.X_SD);
		functionComputesTargetGradientWith2Peaks(Gaussian2DFunction.X_SD + Gaussian2DFunction.PARAMETERS_PER_PEAK);
	}

	@Test
	public void functionComputesYWidthGradientWith2Peaks()
	{
		Assumptions.assumeTrue(null != f2);
		Assumptions.assumeTrue(f2.evaluatesSD1());
		functionComputesTargetGradientWith2Peaks(Gaussian2DFunction.Y_SD);
		functionComputesTargetGradientWith2Peaks(Gaussian2DFunction.Y_SD + Gaussian2DFunction.PARAMETERS_PER_PEAK);
	}

	@Test
	public void functionComputesAngleGradientWith2Peaks()
	{
		Assumptions.assumeTrue(null != f2);
		Assumptions.assumeTrue(f2.evaluatesAngle());
		functionComputesTargetGradientWith2Peaks(Gaussian2DFunction.ANGLE);
		functionComputesTargetGradientWith2Peaks(Gaussian2DFunction.ANGLE + Gaussian2DFunction.PARAMETERS_PER_PEAK);
	}

	private void functionComputesTargetGradientWith2Peaks(int targetParameter)
	{
		final int gradientIndex = findGradientIndex(f2, targetParameter);
		final double[] dyda = new double[f2.gradientIndices().length];
		final double[] dyda2 = new double[dyda.length];
		double[] a;

		final Gaussian2DFunction f2a = GaussianFunctionFactory.create2D(2, maxx, maxy, flags, zModel);
		final Gaussian2DFunction f2b = GaussianFunctionFactory.create2D(2, maxx, maxy, flags, zModel);
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

															f2.initialise(a);

															// Numerically solve gradient.
															// Calculate the step size h to be an exact numerical representation
															final double xx = a[targetParameter];

															// Get h to minimise roundoff error
															final double h = Precision.representableDelta(xx, h_);

															// Evaluate at (x+h) and (x-h)
															a[targetParameter] = xx + h;
															f2a.initialise(a.clone());

															a[targetParameter] = xx - h;
															f2b.initialise(a.clone());

															for (final int x : testx)
																for (final int y : testy)
																{
																	final int i = y * maxx + x;
																	f2.eval(i, dyda);
																	final double value2 = f2a.eval(i, dyda2);
																	final double value3 = f2b.eval(i, dyda2);

																	final double gradient = (value2 - value3) / (2 * h);
																	final double error = DoubleEquality.relativeError(
																			gradient, dyda2[gradientIndex]);
																	s.add(error);
																	ExtraAssertions.assertTrue(
																			(gradient * dyda2[gradientIndex]) >= 0,
																			"%s sign != %s", gradient,
																			dyda2[gradientIndex]);
																	//System.out.printf("[%d,%d] %f == [%d] %f? (%g)\n", x, y, gradient,
																	//		gradientIndex, dyda2[gradientIndex], error);
																	//System.out.printf("[%d,%d] %f == [%d] %f?\n", x, y, gradient, gradientIndex, dyda[gradientIndex]);
																	ExtraAssertions.assertTrue(
																			eq.almostEqualRelativeOrAbsolute(gradient,
																					dyda[gradientIndex]),
																			"%s != %s", gradient, dyda[gradientIndex]);
																}
														}
		TestLog.info(new MessageProvider()
		{
			@Override
			public String getMessage()
			{
				return String.format("functionComputesTargetGradientWith2Peaks %s [%d] %s (error %s +/- %s)\n",
						f2.getClass().getSimpleName(), Gaussian2DFunction.getPeak(targetParameter),
						Gaussian2DFunction.getName(targetParameter), Utils.rounded(s.getMean()),
						Utils.rounded(s.getStandardDeviation()));
			}
		});
	}

	@Test
	public void functionComputesGaussian()
	{
		final double background = 0;
		final int maxx = 30;

		final Gaussian2DFunction f = GaussianFunctionFactory.create2D(1, maxx, maxx, flags, zModel);
		Gaussian2DFunction f2;
		if ((flags & GaussianFunctionFactory.FIT_ERF) == 0)
			f2 = GaussianFunctionFactory.create2D(1, maxx, maxx, GaussianFunctionFactory.FIT_ELLIPTICAL, zModel);
		else
			f2 = GaussianFunctionFactory.create2D(1, maxx, maxx, GaussianFunctionFactory.FIT_ERF_FREE_CIRCLE, zModel);

		final boolean zDepth = (flags & GaussianFunctionFactory.FIT_Z) != 0;

		for (final double signal1 : testsignal1)
			for (final double cx1 : new double[] { maxx / 2 + 0.373f })
				for (final double cy1 : new double[] { maxx / 2 + 0.876f })
					for (final double cz1 : testcz1)
						for (final double[] w1 : testw1)
							for (final double angle1 : testangle1)
							{
								final double[] a = createParameters(background, signal1, cx1, cy1, cz1, w1[0], w1[1],
										angle1);

								f.initialise(a);
								if (zDepth)
								{
									// Change to a standard free circular function
									a[Gaussian2DFunction.X_SD] = zModel.getSx(a[Gaussian2DFunction.Z_POSITION]);
									a[Gaussian2DFunction.Y_SD] = zModel.getSy(a[Gaussian2DFunction.Z_POSITION]);
									a[Gaussian2DFunction.Z_POSITION] = 0;
								}
								f2.initialise(a);
								double sum = 0;
								for (int index = maxx * maxx; index-- > 0;)
								{
									final double r1 = f.eval(index);
									final double r2 = f2.eval(index);
									//System.out.printf("%d,%d r1=%f\n", index%maxx, index/maxx, r1);
									sum += r1;
									final boolean ok = eq2.almostEqualRelativeOrAbsolute(r1, r2);
									ExtraAssertions.assertTrue(ok, "%g != %g @ [%d,%d]", r1, r2, index / maxx,
											index % maxx);
								}

								ExtraAssertions.assertTrue(eq3.almostEqualRelativeOrAbsolute(sum, signal1), "%s != %s",
										sum, signal1);
							}
	}

	protected double[] createParameters(double... args)
	{
		return args;
	}
}
