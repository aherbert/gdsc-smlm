package gdsc.smlm.function.cspline;

import java.util.Arrays;

import gdsc.core.data.DoubleStackTrivalueProvider;
import gdsc.core.ij.Utils;
import gdsc.core.math.interpolation.CustomTricubicInterpolatingFunction;
import gdsc.core.math.interpolation.CustomTricubicInterpolator;
import gdsc.core.utils.DoubleEquality;
import gdsc.core.utils.SimpleArrayUtils;
import gdsc.core.utils.Statistics;
import gdsc.smlm.function.StandardGradient1Procedure;
import gdsc.smlm.function.StandardGradient2Procedure;
import gdsc.smlm.function.StandardValueProcedure;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import gdsc.smlm.function.gaussian.QuadraticAstigmatismZModel;

import org.apache.commons.math3.util.Precision;
import org.junit.Assert;
import org.junit.Assume;
import org.junit.Test;

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
	protected double[] testcx1 = new double[] { 4.9, 5.3 };
	protected double[] testcy1 = new double[] { 4.8, 5.2 };
	protected double[] testcz1 = new double[] { -1.5, 1.0 };
	protected double[] testsignal2 = new double[] { 20, 50 };
	protected double[] testcx2 = new double[] { 4.8, 5.3 };
	protected double[] testcy2 = new double[] { 5.1, 4.9 };
	protected double[] testcz2 = new double[] { -1.9, 0.7 };

	// Different widths to test for non-square function evaluation
	protected int maxx = 10, maxy = 9;
	protected double background = 50;
	protected CubicSplineFunction f1;
	protected CubicSplineFunction f2 = null;

	// Test Astigmatic Gaussian
	final static double gamma = 2;
	final static int zDepth = 5;
	protected QuadraticAstigmatismZModel zModel = new QuadraticAstigmatismZModel(gamma, zDepth);

	final static CubicSplineData splineData;
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
			a[Gaussian2DFunction.Z_POSITION] = z;
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
		cx = a[Gaussian2DFunction.X_POSITION];
		cy = a[Gaussian2DFunction.Y_POSITION];
		cz = splineData.getMaxZ() / scale;
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
						}
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
