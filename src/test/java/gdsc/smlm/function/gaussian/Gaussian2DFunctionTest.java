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
package gdsc.smlm.function.gaussian;

import java.util.Arrays;

import org.apache.commons.math3.util.Precision;
import org.junit.Assert;
import org.junit.Assume;
import org.junit.Test;

import gdsc.core.ij.Utils;
import gdsc.core.utils.DoubleEquality;
import gdsc.core.utils.Statistics;

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

	protected double[] testsignal1 = new double[] { 15, 55, 105 };
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
	protected int maxx = 10, maxy = 9;
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
		{
			testbackground = new double[] { testbackground[0] };
		}
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
		{
			for (int i = 0; i < testw1.length; i++)
			{
				testw1[i][1] = testw1[i][0];
				testw2[i][1] = testw2[i][0];
			}
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
	}

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
		log("Function%d %s %s\n", npeaks, gf.getClass().getName(), Arrays.toString(gradientIndices));

		Assert.assertEquals("Incorrect number of peaks", gf.getNPeaks(), npeaks);

		int p = 0;
		if (gf.evaluatesBackground())
			Assert.assertEquals("Background", 0, gradientIndices[p++]);
		for (int peak = 1, i = 1; peak <= npeaks; peak++, i += Gaussian2DFunction.PARAMETERS_PER_PEAK)
		{
			if (gf.evaluatesSignal())
				Assert.assertEquals(Gaussian2DFunction.getName(i), i, gradientIndices[p++]);
			if (gf.evaluatesPosition())
			{
				Assert.assertEquals(Gaussian2DFunction.getName(i + 1), i + 1, gradientIndices[p++]);
				Assert.assertEquals(Gaussian2DFunction.getName(i + 2), i + 2, gradientIndices[p++]);
			}
			if (gf.evaluatesZ())
				Assert.assertEquals(Gaussian2DFunction.getName(i + 3), i + 3, gradientIndices[p++]);
			if (gf.evaluatesSD0())
				Assert.assertEquals(Gaussian2DFunction.getName(i + 4), i + 4, gradientIndices[p++]);
			if (gf.evaluatesSD1())
				Assert.assertEquals(Gaussian2DFunction.getName(i + 5), i + 5, gradientIndices[p++]);
			if (gf.evaluatesAngle())
				Assert.assertEquals(Gaussian2DFunction.getName(i + 6), i + 6, gradientIndices[p++]);
		}
	}

	@Test
	public void factoryCreatesCorrectFunction()
	{
		Gaussian2DFunction f;

		if (f2 != null)
		{
			f = GaussianFunctionFactory.create2D(2, maxx, maxy, flags, zModel);
			Assert.assertTrue("Incorrect function2", f.getClass() == f2.getClass());
		}
		else
		{
			f = GaussianFunctionFactory.create2D(1, maxx, maxy, flags, zModel);
			Assert.assertTrue("Incorrect function1", f.getClass() == f1.getClass());
		}
	}

	@Test
	public void functionComputesTargetWithAndWithoutGradient()
	{
		double[] dyda = new double[f1.gradientIndices().length];
		double[] a;

		boolean record = true;

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

									f1.initialise(a);

									// Test the frozen version
									int flags = GaussianFunctionFactory.freeze(this.flags, zModel, a);
									Gaussian2DFunction f = GaussianFunctionFactory.create2D(1, f1.getMaxX(),
											f1.getMaxY(), flags, zModel);
									f.initialise(a);
									if (record)
									{
										record = false;
										log("%s %d frozen to %s\n", f1.getClass().getSimpleName(), 1,
												f.getClass().getSimpleName());
									}

									for (int x : testx)
										for (int y : testy)
										{
											int xx = y * maxx + x;
											double y1 = f1.eval(xx, dyda);
											double y2 = f1.eval(xx);

											Assert.assertTrue(y1 + " != " + y2,
													eq2.almostEqualRelativeOrAbsolute(y1, y2));

											y2 = f.eval(xx);

											Assert.assertTrue(y1 + " != frozen " + y2,
													eq2.almostEqualRelativeOrAbsolute(y1, y2));
										}
								}
	}

	@Test
	public void functionComputesBackgroundGradient()
	{
		Assume.assumeTrue(f1.evaluatesBackground());
		functionComputesTargetGradient(Gaussian2DFunction.BACKGROUND);
	}

	@Test
	public void functionComputesSignalGradient()
	{
		Assume.assumeTrue(f1.evaluatesSignal());
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
		Assume.assumeTrue(f1.evaluatesZ());
		functionComputesTargetGradient(Gaussian2DFunction.Z_POSITION);
	}

	@Test
	public void functionComputesXWidthGradient()
	{
		Assume.assumeTrue(f1.evaluatesSD0());
		functionComputesTargetGradient(Gaussian2DFunction.X_SD);
	}

	@Test
	public void functionComputesYWidthGradient()
	{
		Assume.assumeTrue(f1.evaluatesSD1());
		functionComputesTargetGradient(Gaussian2DFunction.Y_SD);
	}

	@Test
	public void functionComputesAngleGradient()
	{
		Assume.assumeTrue(f1.evaluatesAngle());
		functionComputesTargetGradient(Gaussian2DFunction.ANGLE);
	}

	private void functionComputesTargetGradient(int targetParameter)
	{
		int gradientIndex = findGradientIndex(f1, targetParameter);
		double[] dyda = new double[f1.gradientIndices().length];
		double[] dyda2 = new double[dyda.length];
		double[] a;

		Gaussian2DFunction f1a = GaussianFunctionFactory.create2D(1, maxx, maxy, flags, zModel);
		Gaussian2DFunction f1b = GaussianFunctionFactory.create2D(1, maxx, maxy, flags, zModel);
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
									f1.initialise(a);

									// Numerically solve gradient. 
									// Calculate the step size h to be an exact numerical representation
									final double xx = a[targetParameter];

									// Get h to minimise roundoff error
									double h = Precision.representableDelta(xx, h_);

									// Evaluate at (x+h) and (x-h)
									a[targetParameter] = xx + h;
									f1a.initialise(a.clone());

									a[targetParameter] = xx - h;
									f1b.initialise(a.clone());

									for (int x : testx)
										for (int y : testy)
										{
											int i = y * maxx + x;
											f1.eval(i, dyda);
											double value2 = f1a.eval(i, dyda2);
											double value3 = f1b.eval(i, dyda2);

											double gradient = (value2 - value3) / (2 * h);
											double error = DoubleEquality.relativeError(gradient, dyda2[gradientIndex]);
											s.add(error);
											Assert.assertTrue(gradient + " sign != " + dyda2[gradientIndex],
													(gradient * dyda2[gradientIndex]) >= 0);
											//System.out.printf("[%d,%d] %f == [%d] %f? (%g)\n", x, y, gradient,
											//		gradientIndex, dyda2[gradientIndex], error);
											//System.out.printf("[%d,%d] %f == [%d] %f?\n", x, y, gradient, gradientIndex, dyda[gradientIndex]);
											Assert.assertTrue(gradient + " != " + dyda[gradientIndex],
													eq.almostEqualRelativeOrAbsolute(gradient, dyda[gradientIndex]));
										}
								}
		System.out.printf("functionComputesTargetGradient %s %s (error %s +/- %s)\n", f1.getClass().getSimpleName(),
				Gaussian2DFunction.getName(targetParameter), Utils.rounded(s.getMean()),
				Utils.rounded(s.getStandardDeviation()));
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

		boolean record = true;

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

															f2.initialise(a);

															// Test the frozen version
															int flags = GaussianFunctionFactory.freeze(this.flags,
																	zModel, a);
															Gaussian2DFunction f = GaussianFunctionFactory.create2D(2,
																	f2.getMaxX(), f2.getMaxY(), flags, zModel);
															f.initialise(a);
															if (record)
															{
																record = false;
																log("%s %d frozen to %s\n",
																		f2.getClass().getSimpleName(), 2,
																		f.getClass().getSimpleName());
															}

															for (int x : testx)
																for (int y : testy)
																{
																	int xx = y * maxx + x;
																	double y1 = f2.eval(xx, dyda);
																	double y2 = f2.eval(xx);

																	Assert.assertTrue(y1 + " != " + y2,
																			eq2.almostEqualRelativeOrAbsolute(y1, y2));

																	y2 = f.eval(xx);

																	Assert.assertTrue(y1 + " != frozen " + y2,
																			eq2.almostEqualRelativeOrAbsolute(y1, y2));
																}
														}
	}

	@Test
	public void functionComputesBackgroundGradientWith2Peaks()
	{
		Assume.assumeNotNull(f2);
		Assume.assumeTrue(f2.evaluatesBackground());
		functionComputesTargetGradientWith2Peaks(Gaussian2DFunction.BACKGROUND);
	}

	@Test
	public void functionComputesSignalGradientWith2Peaks()
	{
		Assume.assumeNotNull(f2);
		Assume.assumeTrue(f2.evaluatesSignal());
		functionComputesTargetGradientWith2Peaks(Gaussian2DFunction.SIGNAL);
		functionComputesTargetGradientWith2Peaks(Gaussian2DFunction.SIGNAL + Gaussian2DFunction.PARAMETERS_PER_PEAK);
	}

	@Test
	public void functionComputesXGradientWith2Peaks()
	{
		Assume.assumeNotNull(f2);
		functionComputesTargetGradientWith2Peaks(Gaussian2DFunction.X_POSITION);
		functionComputesTargetGradientWith2Peaks(
				Gaussian2DFunction.X_POSITION + Gaussian2DFunction.PARAMETERS_PER_PEAK);
	}

	@Test
	public void functionComputesYGradientWith2Peaks()
	{
		Assume.assumeNotNull(f2);
		functionComputesTargetGradientWith2Peaks(Gaussian2DFunction.Y_POSITION);
		functionComputesTargetGradientWith2Peaks(
				Gaussian2DFunction.Y_POSITION + Gaussian2DFunction.PARAMETERS_PER_PEAK);
	}

	@Test
	public void functionComputesZGradientWith2Peaks()
	{
		Assume.assumeNotNull(f2);
		Assume.assumeTrue(f2.evaluatesZ());
		functionComputesTargetGradientWith2Peaks(Gaussian2DFunction.Z_POSITION);
		functionComputesTargetGradientWith2Peaks(
				Gaussian2DFunction.Z_POSITION + Gaussian2DFunction.PARAMETERS_PER_PEAK);
	}

	@Test
	public void functionComputesXWidthGradientWith2Peaks()
	{
		Assume.assumeNotNull(f2);
		Assume.assumeTrue(f2.evaluatesSD0());
		functionComputesTargetGradientWith2Peaks(Gaussian2DFunction.X_SD);
		functionComputesTargetGradientWith2Peaks(Gaussian2DFunction.X_SD + Gaussian2DFunction.PARAMETERS_PER_PEAK);
	}

	@Test
	public void functionComputesYWidthGradientWith2Peaks()
	{
		Assume.assumeNotNull(f2);
		Assume.assumeTrue(f2.evaluatesSD1());
		functionComputesTargetGradientWith2Peaks(Gaussian2DFunction.Y_SD);
		functionComputesTargetGradientWith2Peaks(Gaussian2DFunction.Y_SD + Gaussian2DFunction.PARAMETERS_PER_PEAK);
	}

	@Test
	public void functionComputesAngleGradientWith2Peaks()
	{
		Assume.assumeNotNull(f2);
		Assume.assumeTrue(f2.evaluatesAngle());
		functionComputesTargetGradientWith2Peaks(Gaussian2DFunction.ANGLE);
		functionComputesTargetGradientWith2Peaks(Gaussian2DFunction.ANGLE + Gaussian2DFunction.PARAMETERS_PER_PEAK);
	}

	private void functionComputesTargetGradientWith2Peaks(int targetParameter)
	{
		int gradientIndex = findGradientIndex(f2, targetParameter);
		double[] dyda = new double[f2.gradientIndices().length];
		double[] dyda2 = new double[dyda.length];
		double[] a;

		Gaussian2DFunction f2a = GaussianFunctionFactory.create2D(2, maxx, maxy, flags, zModel);
		Gaussian2DFunction f2b = GaussianFunctionFactory.create2D(2, maxx, maxy, flags, zModel);
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

															f2.initialise(a);

															// Numerically solve gradient. 
															// Calculate the step size h to be an exact numerical representation
															final double xx = a[targetParameter];

															// Get h to minimise roundoff error
															double h = Precision.representableDelta(xx, h_);

															// Evaluate at (x+h) and (x-h)
															a[targetParameter] = xx + h;
															f2a.initialise(a.clone());

															a[targetParameter] = xx - h;
															f2b.initialise(a.clone());

															for (int x : testx)
																for (int y : testy)
																{
																	int i = y * maxx + x;
																	f2.eval(i, dyda);
																	double value2 = f2a.eval(i, dyda2);
																	double value3 = f2b.eval(i, dyda2);

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
																	//System.out.printf("[%d,%d] %f == [%d] %f?\n", x, y, gradient, gradientIndex, dyda[gradientIndex]);
																	Assert.assertTrue(
																			gradient + " != " + dyda[gradientIndex],
																			eq.almostEqualRelativeOrAbsolute(gradient,
																					dyda[gradientIndex]));
																}
														}
		System.out.printf("functionComputesTargetGradientWith2Peaks %s [%d] %s (error %s +/- %s)\n",
				f2.getClass().getSimpleName(), Gaussian2DFunction.getPeak(targetParameter),
				Gaussian2DFunction.getName(targetParameter), Utils.rounded(s.getMean()),
				Utils.rounded(s.getStandardDeviation()));
	}

	@Test
	public void functionComputesGaussian()
	{
		double background = 0;
		int maxx = 30;

		Gaussian2DFunction f = GaussianFunctionFactory.create2D(1, maxx, maxx, flags, zModel);
		Gaussian2DFunction f2;
		if ((flags & GaussianFunctionFactory.FIT_ERF) == 0)
			f2 = GaussianFunctionFactory.create2D(1, maxx, maxx, GaussianFunctionFactory.FIT_ELLIPTICAL, zModel);
		else
			f2 = GaussianFunctionFactory.create2D(1, maxx, maxx, GaussianFunctionFactory.FIT_ERF_FREE_CIRCLE, zModel);

		boolean zDepth = (flags & GaussianFunctionFactory.FIT_Z) != 0;

		for (double signal1 : testsignal1)
			for (double cx1 : new double[] { maxx / 2 + 0.373f })
				for (double cy1 : new double[] { maxx / 2 + 0.876f })
					for (double cz1 : testcz1)
						for (double[] w1 : testw1)
							for (double angle1 : testangle1)
							{
								double[] a = createParameters(background, signal1, cx1, cy1, cz1, w1[0], w1[1], angle1);

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
									double r1 = f.eval(index);
									double r2 = f2.eval(index);
									//System.out.printf("%d,%d r1=%f\n", index%maxx, index/maxx, r1);
									sum += r1;
									final boolean ok = eq2.almostEqualRelativeOrAbsolute(r1, r2);
									if (!ok)
										Assert.assertTrue(
												String.format("%g != %g @ [%d,%d]", r1, r2, index / maxx, index % maxx),
												ok);
								}

								Assert.assertTrue(sum + " != " + signal1,
										eq3.almostEqualRelativeOrAbsolute(sum, signal1));
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

	protected void log(String format, Object... args)
	{
		System.out.printf(format, args);
	}
}
