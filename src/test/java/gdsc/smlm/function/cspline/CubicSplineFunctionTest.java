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
package gdsc.smlm.function.cspline;

import java.util.Arrays;

import org.apache.commons.math3.util.Precision;
import org.junit.Assert;
import org.junit.Assume;
import org.junit.Test;

import gdsc.core.data.DoubleStackTrivalueProvider;
import gdsc.core.ij.Utils;
import gdsc.core.math.interpolation.CustomTricubicInterpolatingFunction;
import gdsc.core.math.interpolation.CustomTricubicInterpolator;
import gdsc.core.utils.DoubleEquality;
import gdsc.core.utils.Statistics;
import gdsc.core.utils.TurboList;
import gdsc.smlm.function.Gradient1Function;
import gdsc.smlm.function.Gradient1Procedure;
import gdsc.smlm.function.Gradient2Function;
import gdsc.smlm.function.Gradient2Procedure;
import gdsc.smlm.function.StandardGradient1Procedure;
import gdsc.smlm.function.StandardGradient2Procedure;
import gdsc.smlm.function.StandardValueProcedure;
import gdsc.smlm.function.ValueProcedure;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import gdsc.smlm.function.gaussian.QuadraticAstigmatismZModel;
import gdsc.smlm.function.gaussian.erf.ErfGaussian2DFunction;
import gdsc.test.BaseTimingTask;
import gdsc.test.TimingService;

public abstract class CubicSplineFunctionTest
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
	protected double[] testbackground = new double[] { 0, 400 };
	protected double[] testsignal1 = new double[] { 15, 55, 105 };
	// Pick some to fall on the node boundaries as the second order
	// numerical gradients evaluate poorly on the node boundaries.
	protected double[] testcx1 = new double[] { 5.3, 5.0 };
	protected double[] testcy1 = new double[] { 4.5, 5.2 };
	protected double[] testcz1 = new double[] { -1.5, 1.1 };
	protected double[] testsignal2 = new double[] { 20, 50 };
	protected double[] testcx2 = new double[] { 4.8, 5.3 };
	protected double[] testcy2 = new double[] { 5.1, 4.9 };
	protected double[] testcz2 = new double[] { -1.9, 0.7 };

	// Different widths to test for non-square function evaluation
	protected int maxx = 10, maxy = 9;
	protected double background = 50;
	protected CubicSplineFunction f1;
	protected CubicSplineFunction f1f;
	protected CubicSplineFunction f2 = null;
	protected CubicSplineFunction f2f = null;

	// Test Astigmatic Gaussian
	final static double gamma = 2;
	final static int zDepth = 5;
	protected QuadraticAstigmatismZModel zModel = new QuadraticAstigmatismZModel(gamma, zDepth);

	final static CubicSplineData splineData, splineDataFloat;
	final static double cx, cy, cz;
	final static int scale;
	static
	{
		// Create a Guassian PSF twice the size of the test Gaussian for interpolation
		scale = 2;
		QuadraticAstigmatismZModel zModel = new QuadraticAstigmatismZModel(scale * gamma, scale * zDepth);
		int size = 40;
		Gaussian2DFunction f = GaussianFunctionFactory.create2D(1, size, size, GaussianFunctionFactory.FIT_ASTIGMATISM,
				zModel);
		double[] a = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
		a[Gaussian2DFunction.SIGNAL] = 1;
		a[Gaussian2DFunction.X_POSITION] = size / scale;
		a[Gaussian2DFunction.Y_POSITION] = size / scale;
		a[Gaussian2DFunction.X_SD] = scale;
		a[Gaussian2DFunction.Y_SD] = scale;

		// Create the Gaussian data for different z-depths
		int minz = -scale * zDepth;
		int maxz = -minz;
		double[][] val = new double[maxz - minz + 1][];
		StandardValueProcedure p = new StandardValueProcedure();
		for (int z = minz, i = 0; z <= maxz; z++, i++)
		{
			a[CubicSplineFunction.Z_POSITION] = z;
			val[i] = p.getValues(f, a);
		}
		DoubleStackTrivalueProvider fval = new DoubleStackTrivalueProvider(val, size, size);

		//Utils.display("CubicSplineData", val, size, size);

		//@formatter:off
		CustomTricubicInterpolatingFunction function = new CustomTricubicInterpolator.Builder()
				// The axis value are ignored ...
				.setIntegerAxisValues(true)
				.setFValue(fval)
				.interpolate();
		//@formatter:on
		splineData = new CubicSplineData(function);
		cx = a[CubicSplineFunction.X_POSITION];
		cy = a[CubicSplineFunction.Y_POSITION];
		cz = splineData.getMaxZ() / scale;

		function.toSinglePrecision();
		splineDataFloat = new CubicSplineData(function);
	}

	public CubicSplineFunctionTest()
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

		postInit();
	}

	/**
	 * Create the CubicSplineFunction for 1 and 2 peaks. Creates the flags for the factory
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

	private void checkGradientIndices(int npeaks, CubicSplineFunction cf)
	{
		if (cf == null)
			return;

		int[] gradientIndices = cf.gradientIndices();
		log("Function%d %s %s\n", npeaks, cf.getClass().getName(), Arrays.toString(gradientIndices));

		Assert.assertEquals("Incorrect number of peaks", cf.getN(), npeaks);

		int p = 0;
		if (cf.evaluatesBackground())
			Assert.assertEquals("Background", 0, gradientIndices[p++]);
		for (int peak = 1, i = 1; peak <= npeaks; peak++, i += CubicSplineFunction.PARAMETERS_PER_PEAK)
		{
			if (cf.evaluatesSignal())
				Assert.assertEquals(CubicSplineFunction.getName(i), i, gradientIndices[p++]);
			if (cf.evaluatesPosition())
			{
				Assert.assertEquals(CubicSplineFunction.getName(i + 1), i + 1, gradientIndices[p++]);
				Assert.assertEquals(CubicSplineFunction.getName(i + 2), i + 2, gradientIndices[p++]);
			}
			if (cf.evaluatesZ())
				Assert.assertEquals(CubicSplineFunction.getName(i + 3), i + 3, gradientIndices[p++]);
		}
	}

	@Test
	public void factoryCreatesCorrectFunction()
	{
		CubicSplineFunction f;

		if (f2 != null)
		{
			f = CubicSplineFunctionFactory.createCubicSplineFunction(splineData, maxx, maxy, cx, cy, cz, 2, 2);
			Assert.assertTrue("Incorrect function2", f.getClass() == f2.getClass());
		}
		else
		{
			f = CubicSplineFunctionFactory.createCubicSplineFunction(splineData, maxx, maxy, cx, cy, cz, 2, 1);
			Assert.assertTrue("Incorrect function1", f.getClass() == f1.getClass());
		}
	}

	@Test
	public void functionComputesTargetWithAndWithoutGradient()
	{
		StandardValueProcedure p0 = new StandardValueProcedure();
		StandardGradient1Procedure p1 = new StandardGradient1Procedure();
		StandardGradient2Procedure p2 = new StandardGradient2Procedure();

		for (double background : testbackground)
			// Peak 1
			for (double signal1 : testsignal1)
				for (double cx1 : testcx1)
					for (double cy1 : testcy1)
						for (double cz1 : testcz1)
						{
							double[] a = createParameters(background, signal1, cx1, cy1, cz1);

							double[] e = p0.getValues(f1, a);
							double[] o1 = p1.getValues(f1, a);
							double[] o2 = p2.getValues(f1, a);

							Assert.assertArrayEquals(e, o1, 0);
							Assert.assertArrayEquals(e, o2, 0);
							for (int i = e.length; i-- > 0;)
								Assert.assertArrayEquals(p1.dyda[i], p2.dyda[i], 0);
						}
	}

	@Test
	public void functionComputesBackgroundGradient1()
	{
		Assume.assumeTrue(f1.evaluatesBackground());
		functionComputesTargetGradient1(CubicSplineFunction.BACKGROUND);
	}

	@Test
	public void functionComputesSignalGradient1()
	{
		Assume.assumeTrue(f1.evaluatesSignal());
		functionComputesTargetGradient1(CubicSplineFunction.SIGNAL);
	}

	@Test
	public void functionComputesXGradient1()
	{
		functionComputesTargetGradient1(CubicSplineFunction.X_POSITION);
	}

	@Test
	public void functionComputesYGradient1()
	{
		functionComputesTargetGradient1(CubicSplineFunction.Y_POSITION);
	}

	@Test
	public void functionComputesZGradient1()
	{
		Assume.assumeTrue(f1.evaluatesZ());
		functionComputesTargetGradient1(CubicSplineFunction.Z_POSITION);
	}

	private void functionComputesTargetGradient1(int targetParameter)
	{
		int gradientIndex = findGradientIndex(f1, targetParameter);

		Statistics s = new Statistics();

		StandardValueProcedure p1a = new StandardValueProcedure();
		StandardValueProcedure p1b = new StandardValueProcedure();
		StandardGradient1Procedure p2 = new StandardGradient1Procedure();

		for (double background : testbackground)
			// Peak 1
			for (double signal1 : testsignal1)
				for (double cx1 : testcx1)
					for (double cy1 : testcy1)
						for (double cz1 : testcz1)
						{
							double[] a = createParameters(background, signal1, cx1, cy1, cz1);

							//System.out.println(java.util.Arrays.toString(a));

							// Evaluate all gradients 
							p2.getValues(f1, a);

							// Numerically solve gradient. 
							// Calculate the step size h to be an exact numerical representation
							final double xx = a[targetParameter];

							// Get h to minimise roundoff error
							double h = Precision.representableDelta(xx, h_);

							// Evaluate at (x+h) and (x-h)
							a[targetParameter] = xx + h;
							p1a.getValues(f1, a);

							a[targetParameter] = xx - h;
							p1b.getValues(f1, a);

							// Only test close to the XY centre
							for (int x : testx)
								for (int y : testy)
								{
									int i = y * maxx + x;
									double high = p1a.values[i];
									double low = p1b.values[i];

									double gradient = (high - low) / (2 * h);
									double dyda = p2.dyda[i][gradientIndex];
									double error = DoubleEquality.relativeError(gradient, dyda);
									s.add(error);
									Assert.assertTrue(gradient + " sign != " + dyda, (gradient * dyda) >= 0);
									//System.out.printf("[%d,%d] %f == [%d] %f? (%g)\n", x, y, gradient, gradientIndex, dyda, error);
									Assert.assertTrue(gradient + " != " + dyda,
											eq.almostEqualRelativeOrAbsolute(gradient, dyda));
								}
						}
		System.out.printf("functionComputesTargetGradient1 %s %s (error %s +/- %s)\n", f1.getClass().getSimpleName(),
				CubicSplineFunction.getName(targetParameter), Utils.rounded(s.getMean()),
				Utils.rounded(s.getStandardDeviation()));
	}

	protected int findGradientIndex(CubicSplineFunction f, int targetParameter)
	{
		int i = f.findGradientIndex(targetParameter);
		Assert.assertTrue("Cannot find gradient index", i >= 0);
		return i;
	}

	@Test
	public void functionComputesBackgroundGradient2()
	{
		Assume.assumeTrue(f1.evaluatesBackground());
		functionComputesTargetGradient2(CubicSplineFunction.BACKGROUND);
	}

	@Test
	public void functionComputesSignalGradient2()
	{
		Assume.assumeTrue(f1.evaluatesSignal());
		functionComputesTargetGradient2(CubicSplineFunction.SIGNAL);
	}

	@Test
	public void functionComputesXGradient2()
	{
		functionComputesTargetGradient2(CubicSplineFunction.X_POSITION);
	}

	@Test
	public void functionComputesYGradient2()
	{
		functionComputesTargetGradient2(CubicSplineFunction.Y_POSITION);
	}

	@Test
	public void functionComputesZGradient2()
	{
		Assume.assumeTrue(f1.evaluatesZ());
		functionComputesTargetGradient2(CubicSplineFunction.Z_POSITION);
	}

	private void functionComputesTargetGradient2(int targetParameter)
	{
		int gradientIndex = findGradientIndex(f1, targetParameter);

		Statistics s = new Statistics();

		StandardGradient1Procedure p1a = new StandardGradient1Procedure();
		StandardGradient1Procedure p1b = new StandardGradient1Procedure();
		StandardGradient2Procedure p2 = new StandardGradient2Procedure();

		for (double background : testbackground)
			// Peak 1
			for (double signal1 : testsignal1)
				for (double cx1 : testcx1)
					for (double cy1 : testcy1)
						for (double cz1 : testcz1)
						{
							double[] a = createParameters(background, signal1, cx1, cy1, cz1);

							//System.out.println(java.util.Arrays.toString(a));

							f1.initialise2(a);
							boolean test = !f1.isNodeBoundary(gradientIndex);
							// Comment out when printing errors
							if (!test)
								continue;

							// Evaluate all gradients 
							p2.getValues(f1, a);

							// Numerically solve gradient. 
							// Calculate the step size h to be an exact numerical representation
							final double xx = a[targetParameter];

							// Get h to minimise roundoff error
							double h = Precision.representableDelta(xx, h_);

							// Evaluate at (x+h) and (x-h)
							a[targetParameter] = xx + h;
							p1a.getValues(f1, a);

							a[targetParameter] = xx - h;
							p1b.getValues(f1, a);

							// Only test close to the XY centre
							for (int x : testx)
								for (int y : testy)
								{
									int i = y * maxx + x;
									double high = p1a.dyda[i][gradientIndex];
									double low = p1b.dyda[i][gradientIndex];

									double gradient = (high - low) / (2 * h);
									double d2yda2 = p2.d2yda2[i][gradientIndex];
									double error = DoubleEquality.relativeError(gradient, d2yda2);
									//System.out.printf("[%d,%d] %f == [%d] %f? (%g)\n", x, y, gradient, gradientIndex,	d2yda2, error);
									if (test)
									{
										s.add(error);
										Assert.assertTrue(gradient + " sign != " + d2yda2, (gradient * d2yda2) >= 0);
										Assert.assertTrue(gradient + " != " + d2yda2,
												eq.almostEqualRelativeOrAbsolute(gradient, d2yda2));
									}
								}
						}
		System.out.printf("functionComputesTargetGradient2 %s %s (error %s +/- %s)\n", f1.getClass().getSimpleName(),
				CubicSplineFunction.getName(targetParameter), Utils.rounded(s.getMean()),
				Utils.rounded(s.getStandardDeviation()));
	}

	@Test
	public void functionComputesTargetWithAndWithoutGradientWith2Peaks()
	{
		if (f2 == null)
			return;

		StandardValueProcedure p0 = new StandardValueProcedure();
		StandardGradient1Procedure p1 = new StandardGradient1Procedure();
		StandardGradient2Procedure p2 = new StandardGradient2Procedure();

		for (double background : testbackground)
			// Peak 1
			for (double signal1 : testsignal1)
				for (double cx1 : testcx1)
					for (double cy1 : testcy1)
						for (double cz1 : testcz1)
							// Peak 2
							for (double signal2 : testsignal2)
								for (double cx2 : testcx2)
									for (double cy2 : testcy2)
										for (double cz2 : testcz2)
										{
											double[] a = createParameters(background, signal1, cx1, cy1, cz1, signal2,
													cx2, cy2, cz2);

											double[] e = p0.getValues(f1, a);
											double[] o1 = p1.getValues(f1, a);
											double[] o2 = p2.getValues(f1, a);

											Assert.assertArrayEquals(e, o1, 0);
											Assert.assertArrayEquals(e, o2, 0);
											for (int i = e.length; i-- > 0;)
												Assert.assertArrayEquals(p1.dyda[i], p2.dyda[i], 0);
										}
	}

	@Test
	public void functionComputesBackgroundGradient1With2Peaks()
	{
		Assume.assumeNotNull(f2);
		Assume.assumeTrue(f2.evaluatesBackground());
		functionComputesTargetGradient1With2Peaks(CubicSplineFunction.BACKGROUND);
		functionComputesTargetGradient1With2Peaks(
				CubicSplineFunction.BACKGROUND + CubicSplineFunction.PARAMETERS_PER_PEAK);
	}

	@Test
	public void functionComputesSignalGradient1With2Peaks()
	{
		Assume.assumeNotNull(f2);
		Assume.assumeTrue(f2.evaluatesSignal());
		functionComputesTargetGradient1With2Peaks(CubicSplineFunction.SIGNAL);
		functionComputesTargetGradient1With2Peaks(CubicSplineFunction.SIGNAL + CubicSplineFunction.PARAMETERS_PER_PEAK);
	}

	@Test
	public void functionComputesXGradient1With2Peaks()
	{
		Assume.assumeNotNull(f2);
		functionComputesTargetGradient1With2Peaks(CubicSplineFunction.X_POSITION);
		functionComputesTargetGradient1With2Peaks(
				CubicSplineFunction.X_POSITION + CubicSplineFunction.PARAMETERS_PER_PEAK);
	}

	@Test
	public void functionComputesYGradient1With2Peaks()
	{
		Assume.assumeNotNull(f2);
		functionComputesTargetGradient1With2Peaks(CubicSplineFunction.Y_POSITION);
		functionComputesTargetGradient1With2Peaks(
				CubicSplineFunction.Y_POSITION + CubicSplineFunction.PARAMETERS_PER_PEAK);
	}

	@Test
	public void functionComputesZGradient1With2Peaks()
	{
		Assume.assumeNotNull(f2);
		Assume.assumeTrue(f2.evaluatesZ());
		functionComputesTargetGradient1With2Peaks(CubicSplineFunction.Z_POSITION);
		functionComputesTargetGradient1With2Peaks(
				CubicSplineFunction.Z_POSITION + CubicSplineFunction.PARAMETERS_PER_PEAK);
	}

	private void functionComputesTargetGradient1With2Peaks(int targetParameter)
	{
		int gradientIndex = findGradientIndex(f2, targetParameter);

		Statistics s = new Statistics();

		StandardValueProcedure p1a = new StandardValueProcedure();
		StandardValueProcedure p1b = new StandardValueProcedure();
		StandardGradient1Procedure p2 = new StandardGradient1Procedure();

		for (double background : testbackground)
			// Peak 1
			for (double signal1 : testsignal1)
				for (double cx1 : testcx1)
					for (double cy1 : testcy1)
						for (double cz1 : testcz1)
							// Peak 2
							for (double signal2 : testsignal2)
								for (double cx2 : testcx2)
									for (double cy2 : testcy2)
										for (double cz2 : testcz2)
										{
											double[] a = createParameters(background, signal1, cx1, cy1, cz1, signal2,
													cx2, cy2, cz2);

											//System.out.println(java.util.Arrays.toString(a));

											// Evaluate all gradients 
											p2.getValues(f2, a);

											// Numerically solve gradient. 
											// Calculate the step size h to be an exact numerical representation
											final double xx = a[targetParameter];

											// Get h to minimise roundoff error
											double h = Precision.representableDelta(xx, h_);

											// Evaluate at (x+h) and (x-h)
											a[targetParameter] = xx + h;
											p1a.getValues(f2, a);

											a[targetParameter] = xx - h;
											p1b.getValues(f2, a);

											// Only test close to the XY centre
											for (int x : testx)
												for (int y : testy)
												{
													int i = y * maxx + x;
													double high = p1a.values[i];
													double low = p1b.values[i];

													double gradient = (high - low) / (2 * h);
													double dyda = p2.dyda[i][gradientIndex];
													double error = DoubleEquality.relativeError(gradient, dyda);
													s.add(error);
													Assert.assertTrue(gradient + " sign != " + dyda,
															(gradient * dyda) >= 0);
													//System.out.printf("[%d,%d] %f == [%d] %f? (%g)\n", x, y, gradient, gradientIndex, dyda, error);
													Assert.assertTrue(gradient + " != " + dyda,
															eq.almostEqualRelativeOrAbsolute(gradient, dyda));
												}
										}
		System.out.printf("functionComputesTargetGradient1With2Peaks %s %s (error %s +/- %s)\n",
				f1.getClass().getSimpleName(), CubicSplineFunction.getName(targetParameter), Utils.rounded(s.getMean()),
				Utils.rounded(s.getStandardDeviation()));
	}

	@Test
	public void functionComputesBackgroundGradient2With2Peaks()
	{
		Assume.assumeNotNull(f2);
		Assume.assumeTrue(f2.evaluatesBackground());
		functionComputesTargetGradient2With2Peaks(CubicSplineFunction.BACKGROUND);
		functionComputesTargetGradient2With2Peaks(
				CubicSplineFunction.BACKGROUND + CubicSplineFunction.PARAMETERS_PER_PEAK);
	}

	@Test
	public void functionComputesSignalGradient2With2Peaks()
	{
		Assume.assumeNotNull(f2);
		Assume.assumeTrue(f2.evaluatesSignal());
		functionComputesTargetGradient2With2Peaks(CubicSplineFunction.SIGNAL);
		functionComputesTargetGradient2With2Peaks(CubicSplineFunction.SIGNAL + CubicSplineFunction.PARAMETERS_PER_PEAK);
	}

	@Test
	public void functionComputesXGradient2With2Peaks()
	{
		Assume.assumeNotNull(f2);
		functionComputesTargetGradient2With2Peaks(CubicSplineFunction.X_POSITION);
		functionComputesTargetGradient2With2Peaks(
				CubicSplineFunction.X_POSITION + CubicSplineFunction.PARAMETERS_PER_PEAK);
	}

	@Test
	public void functionComputesYGradient2With2Peaks()
	{
		Assume.assumeNotNull(f2);
		functionComputesTargetGradient2With2Peaks(CubicSplineFunction.Y_POSITION);
		functionComputesTargetGradient2With2Peaks(
				CubicSplineFunction.Y_POSITION + CubicSplineFunction.PARAMETERS_PER_PEAK);
	}

	@Test
	public void functionComputesZGradient2With2Peaks()
	{
		Assume.assumeNotNull(f2);
		Assume.assumeTrue(f2.evaluatesZ());
		functionComputesTargetGradient1With2Peaks(CubicSplineFunction.Z_POSITION);
		functionComputesTargetGradient1With2Peaks(
				CubicSplineFunction.Z_POSITION + CubicSplineFunction.PARAMETERS_PER_PEAK);
	}

	private void functionComputesTargetGradient2With2Peaks(int targetParameter)
	{
		int gradientIndex = findGradientIndex(f2, targetParameter);

		Statistics s = new Statistics();

		StandardGradient1Procedure p1a = new StandardGradient1Procedure();
		StandardGradient1Procedure p1b = new StandardGradient1Procedure();
		StandardGradient2Procedure p2 = new StandardGradient2Procedure();

		for (double background : testbackground)
			// Peak 1
			for (double signal1 : testsignal1)
				for (double cx1 : testcx1)
					for (double cy1 : testcy1)
						for (double cz1 : testcz1)
							// Peak 2
							for (double signal2 : testsignal2)
								for (double cx2 : testcx2)
									for (double cy2 : testcy2)
										for (double cz2 : testcz2)
										{
											double[] a = createParameters(background, signal1, cx1, cy1, cz1, signal2,
													cx2, cy2, cz2);

											//System.out.println(java.util.Arrays.toString(a));

											f2.initialise2(a);
											boolean test = !f2.isNodeBoundary(gradientIndex);
											// Comment out when printing errors
											if (!test)
												continue;

											// Evaluate all gradients 
											p2.getValues(f2, a);

											// Numerically solve gradient. 
											// Calculate the step size h to be an exact numerical representation
											final double xx = a[targetParameter];

											// Get h to minimise roundoff error
											double h = Precision.representableDelta(xx, h_);

											// Evaluate at (x+h) and (x-h)
											a[targetParameter] = xx + h;
											p1a.getValues(f2, a);

											a[targetParameter] = xx - h;
											p1b.getValues(f2, a);

											// Only test close to the XY centre
											for (int x : testx)
												for (int y : testy)
												{
													int i = y * maxx + x;
													double high = p1a.dyda[i][gradientIndex];
													double low = p1b.dyda[i][gradientIndex];

													double gradient = (high - low) / (2 * h);
													double d2yda2 = p2.d2yda2[i][gradientIndex];
													double error = DoubleEquality.relativeError(gradient, d2yda2);
													//System.out.printf("[%d,%d] %f == [%d] %f? (%g)\n", x, y, gradient, gradientIndex, d2yda2, error);
													if (test)
													{
														s.add(error);
														Assert.assertTrue(gradient + " sign != " + d2yda2,
																(gradient * d2yda2) >= 0);
														Assert.assertTrue(gradient + " != " + d2yda2,
																eq.almostEqualRelativeOrAbsolute(gradient, d2yda2));
													}
												}
										}
		System.out.printf("functionComputesTargetGradient2With2Peaks %s %s (error %s +/- %s)\n",
				f1.getClass().getSimpleName(), CubicSplineFunction.getName(targetParameter), Utils.rounded(s.getMean()),
				Utils.rounded(s.getStandardDeviation()));
	}

	@Test
	public void runSpeedTestWith1Peak()
	{
		speedTest(1, 0);
		speedTest(1, 1);
		speedTest(1, 2);
	}

	@Test
	public void runSpeedTestWith2Peaks()
	{
		speedTest(2, 0);
		speedTest(2, 1);
		speedTest(2, 2);
	}

	private class FunctionTimingTask extends BaseTimingTask
			implements ValueProcedure, Gradient1Procedure, Gradient2Procedure
	{
		Gradient1Function f1;
		Gradient2Function f2;
		double[][] x;
		int order;
		double s;

		public FunctionTimingTask(Gradient1Function f, double[][] x, int order)
		{
			super(f.getClass().getSimpleName() + " " + order);
			this.f1 = f;
			if (order > 1)
				throw new IllegalArgumentException("Gradient1Function for order>1");
			this.x = x;
			this.order = order;
		}

		public FunctionTimingTask(Gradient2Function f, double[][] x, int order)
		{
			super(f.getClass().getSimpleName() + " " + order);
			this.f1 = f;
			this.f2 = f;
			this.x = x;
			this.order = order;
		}

		public FunctionTimingTask(Gradient2Function f, double[][] x, int order, String suffix)
		{
			super(f.getClass().getSimpleName() + " " + order + suffix);
			this.f1 = f;
			this.f2 = f;
			this.x = x;
			this.order = order;
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
			s = 0;
			//f = f.copy();
			if (order == 0)
			{
				for (int i = 0; i < x.length; i++)
				{
					f1.initialise0(x[i]);
					f1.forEach((ValueProcedure) this);
				}
			}
			else if (order == 1)
			{
				for (int i = 0; i < x.length; i++)
				{
					f1.initialise1(x[i]);
					f1.forEach((Gradient1Procedure) this);
				}
			}
			else
			{
				for (int i = 0; i < x.length; i++)
				{
					f2.initialise2(x[i]);
					f2.forEach((Gradient2Procedure) this);
				}
			}
			return s;
		}

		@Override
		public void execute(double value)
		{
			s += value;
		}

		@Override
		public void execute(double value, double[] dy_da)
		{
			s += value;
		}

		@Override
		public void execute(double value, double[] dy_da, double[] d2y_da2)
		{
			s += value;
		}
	}

	private void speedTest(int n, int order)
	{
		CubicSplineFunction cf = (n == 2) ? f2 : f1;
		Assume.assumeNotNull(cf);
		CubicSplineFunction cff = (n == 2) ? f2f : f1f;
		ErfGaussian2DFunction gf = (ErfGaussian2DFunction) GaussianFunctionFactory.create2D(n, maxx, maxy,
				GaussianFunctionFactory.FIT_ASTIGMATISM, zModel);
		Gaussian2DFunction gf2 = (order < 2) ? GaussianFunctionFactory.create2D(n, maxx, maxy,
				GaussianFunctionFactory.FIT_SIMPLE_FREE_CIRCLE, zModel) : null;
		TurboList<double[]> l1 = new TurboList<double[]>();
		TurboList<double[]> l2 = new TurboList<double[]>();
		TurboList<double[]> l3 = new TurboList<double[]>();
		double[] a = new double[1 + n * CubicSplineFunction.PARAMETERS_PER_PEAK];
		double[] b = new double[1 + n * Gaussian2DFunction.PARAMETERS_PER_PEAK];
		double[] bb = null;
		a[CubicSplineFunction.BACKGROUND] = 0.1;
		b[Gaussian2DFunction.BACKGROUND] = 0.1;
		for (int i = 0; i < n; i++)
		{
			a[i * CubicSplineFunction.PARAMETERS_PER_PEAK + CubicSplineFunction.SIGNAL] = 10;
			b[i * Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.SIGNAL] = 10;
		}
		if (n == 2)
		{
			// Fix second peak parameters
			a[CubicSplineFunction.PARAMETERS_PER_PEAK + CubicSplineFunction.X_POSITION] = testcx1[0];
			a[CubicSplineFunction.PARAMETERS_PER_PEAK + CubicSplineFunction.Y_POSITION] = testcy1[0];
			a[CubicSplineFunction.PARAMETERS_PER_PEAK + CubicSplineFunction.Z_POSITION] = testcz1[0];
			b[Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.X_POSITION] = testcx1[0];
			b[Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.Y_POSITION] = testcy1[0];
			b[Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.Z_POSITION] = testcz1[0];
		}
		if (gf2 != null)
		{
			bb = b.clone();
			if (n == 2)
			{
				// Fix second peak parameters
				bb[Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.X_SD] = zModel.getSx(testcz1[0]);
				bb[Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.Y_SD] = zModel.getSy(testcz1[0]);
			}
		}
		for (int x = 0; x <= maxx; x++)
		{
			a[CubicSplineFunction.X_POSITION] = x;
			b[Gaussian2DFunction.X_POSITION] = x;
			for (int y = 0; y <= maxy; y++)
			{
				a[CubicSplineFunction.Y_POSITION] = y;
				b[Gaussian2DFunction.Y_POSITION] = y;
				for (int z = -zDepth; z <= zDepth; z++)
				{
					a[CubicSplineFunction.Z_POSITION] = z;
					b[Gaussian2DFunction.Z_POSITION] = z;
					l1.add(a.clone());
					l2.add(b.clone());
					if (gf2 != null)
					{
						bb[Gaussian2DFunction.X_SD] = zModel.getSx(z);
						bb[Gaussian2DFunction.Y_SD] = zModel.getSy(z);
						l3.add(bb.clone());
					}
				}
			}
		}
		double[][] x1 = l1.toArray(new double[l1.size()][]);
		double[][] x2 = l2.toArray(new double[l2.size()][]);
		double[][] x3 = l3.toArray(new double[l3.size()][]);

		TimingService ts = new TimingService(5);
		ts.execute(new FunctionTimingTask(gf, x2, order));
		if (gf2 != null)
			ts.execute(new FunctionTimingTask(gf2, x3, order));
		ts.execute(new FunctionTimingTask(cf, x1, order));
		ts.execute(new FunctionTimingTask(cff, x1, order, " single-precision"));

		int size = ts.getSize();
		ts.repeat(size);
		ts.report(size);
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
