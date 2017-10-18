package gdsc.smlm.function.cspline;

import java.util.Arrays;

import gdsc.core.ij.Utils;
import gdsc.core.math.interpolation.CustomTricubicInterpolatingFunction;
import gdsc.core.utils.DoubleEquality;
import gdsc.core.utils.Statistics;
import gdsc.smlm.function.gaussian.GaussianFunctionFactory;

import org.apache.commons.math3.util.Precision;
import org.junit.Assert;
import org.junit.Assume;
import org.junit.Test;

public abstract class CubicSplineFunctionTest
{
	final static CubicSplineData splineData;
	static
	{
		// Create a function for interpolation
		
		
		CustomTricubicInterpolatingFunction function = null;
		splineData = new CubicSplineData(function);
	}
	
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
	protected int flags;

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

	private void checkGradientIndices(int npeaks, CubicSplineFunction gf)
	{
		if (gf == null)
			return;

		int[] gradientIndices = gf.gradientIndices();
		log("Function%d %s %s\n", npeaks, gf.getClass().getName(), Arrays.toString(gradientIndices));

		Assert.assertEquals("Incorrect number of peaks", gf.getN(), npeaks);

		int p = 0;
		if (gf.evaluatesBackground())
			Assert.assertEquals("Background", 0, gradientIndices[p++]);
		for (int peak = 1, i = 1; peak <= npeaks; peak++, i += CubicSplineFunction.PARAMETERS_PER_PEAK)
		{
			if (gf.evaluatesSignal())
				Assert.assertEquals(CubicSplineFunction.getName(i), i, gradientIndices[p++]);
			if (gf.evaluatesPosition())
			{
				Assert.assertEquals(CubicSplineFunction.getName(i + 1), i + 1, gradientIndices[p++]);
				Assert.assertEquals(CubicSplineFunction.getName(i + 2), i + 2, gradientIndices[p++]);
			}
			if (gf.evaluatesZ())
				Assert.assertEquals(CubicSplineFunction.getName(i + 3), i + 3, gradientIndices[p++]);
		}
	}

	@Test
	public void factoryCreatesCorrectFunction()
	{
		CubicSplineFunction f;

		if (f2 != null)
		{
			//f = GaussianFunctionFactory.create2D(2, maxx, maxy, flags, null);
			//Assert.assertTrue("Incorrect function2", f.getClass() == f2.getClass());
		}
		else
		{
			//f = GaussianFunctionFactory.create2D(1, maxx, maxy, flags, null);
			//Assert.assertTrue("Incorrect function1", f.getClass() == f1.getClass());
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
