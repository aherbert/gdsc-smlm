package gdsc.smlm.fitting.function.gaussian;

import gdsc.smlm.TestSettings;
import gdsc.smlm.fitting.function.Gaussian2DFunction;
import gdsc.smlm.fitting.function.GaussianFunction;
import gdsc.smlm.fitting.function.GaussianFunctionFactory;
import gdsc.smlm.fitting.function.gaussian.EllipticalGaussian2DFunction;
import gdsc.smlm.fitting.utils.FloatEquality;

import java.util.Random;
import java.util.ArrayList;

import org.apache.commons.math3.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

/**
 * Contains speed tests for the fastest method for calculating the Hessian and gradient vector
 * for use in NonLinearFit
 */
public class SpeedTest
{
	int Single = 1;
	int Double = 2;

	private int MAX_ITER = 20000;
	private int blockWidth = 10;
	private float Background = 20;
	private float Amplitude = 10;
	private float Xpos = 5;
	private float Ypos = 5;
	private float Xwidth = 5;
	private Random rand;
	private static ArrayList<float[]> paramsListSinglePeak = null;
	private static ArrayList<float[]> yListSinglePeak;
	private static int[] x;
	private static ArrayList<float[]> paramsListDoublePeak;
	private static ArrayList<float[]> yListDoublePeak;

	public SpeedTest()
	{
		// TODO - the data generation could be static in the Gaussian2DFunctionTest since 
		// all gaussian data generation should be similar

		rand = new Random(30051977);

		if (paramsListSinglePeak == null)
		{
			paramsListSinglePeak = new ArrayList<float[]>(MAX_ITER);
			yListSinglePeak = new ArrayList<float[]>(MAX_ITER);
			x = createData(1, MAX_ITER, paramsListSinglePeak, yListSinglePeak);
			paramsListDoublePeak = new ArrayList<float[]>(MAX_ITER);
			yListDoublePeak = new ArrayList<float[]>(MAX_ITER);
			x = createData(2, MAX_ITER, paramsListDoublePeak, yListDoublePeak);
		}
	}

	@Test
	public void freeCircularComputesSameAsEllipticalSinglePeak()
	{
		f1ComputesSameAsf2(Single, GaussianFunctionFactory.FIT_FREE_CIRCLE, GaussianFunctionFactory.FIT_ELLIPTICAL);
	}

	@Test
	public void freeCircularFasterThanEllipticalSinglePeak()
	{
		f1FasterThanf2(Single, GaussianFunctionFactory.FIT_FREE_CIRCLE, GaussianFunctionFactory.FIT_ELLIPTICAL);
	}

	@Test
	public void circularComputesSameAsFreeCircularSinglePeak()
	{
		f1ComputesSameAsf2(Single, GaussianFunctionFactory.FIT_CIRCLE, GaussianFunctionFactory.FIT_FREE_CIRCLE);
	}

	@Test
	public void circularFasterThanFreeCircularSinglePeak()
	{
		f1FasterThanf2(Single, GaussianFunctionFactory.FIT_CIRCLE, GaussianFunctionFactory.FIT_FREE_CIRCLE);
	}

	@Test
	public void fixedComputesSameAsFreeCircularSinglePeak()
	{
		f1ComputesSameAsf2(Single, GaussianFunctionFactory.FIT_FIXED, GaussianFunctionFactory.FIT_FREE_CIRCLE);
	}

	@Test
	public void fixedFasterThanFreeCircularSinglePeak()
	{
		f1FasterThanf2(Single, GaussianFunctionFactory.FIT_FIXED, GaussianFunctionFactory.FIT_FREE_CIRCLE);
	}

	@Test
	public void freeCircularComputesSameAsEllipticalSinglePeakNB()
	{
		f1ComputesSameAsf2(Single, GaussianFunctionFactory.FIT_NB_FREE_CIRCLE,
				GaussianFunctionFactory.FIT_NB_ELLIPTICAL);
	}

	@Test
	public void freeCircularFasterThanEllipticalSinglePeakNB()
	{
		f1FasterThanf2(Single, GaussianFunctionFactory.FIT_NB_FREE_CIRCLE, GaussianFunctionFactory.FIT_NB_ELLIPTICAL);
	}

	@Test
	public void circularComputesSameAsFreeCircularSinglePeakNB()
	{
		f1ComputesSameAsf2(Single, GaussianFunctionFactory.FIT_NB_CIRCLE, GaussianFunctionFactory.FIT_NB_FREE_CIRCLE);
	}

	@Test
	public void circularFasterThanFreeCircularSinglePeakNB()
	{
		f1FasterThanf2(Single, GaussianFunctionFactory.FIT_NB_CIRCLE, GaussianFunctionFactory.FIT_NB_FREE_CIRCLE);
	}

	@Test
	public void fixedComputesSameAsFreeCircularSinglePeakNB()
	{
		f1ComputesSameAsf2(Single, GaussianFunctionFactory.FIT_NB_FIXED, GaussianFunctionFactory.FIT_NB_FREE_CIRCLE);
	}

	@Test
	public void fixedFasterThanFreeCircularSinglePeakNB()
	{
		f1FasterThanf2(Single, GaussianFunctionFactory.FIT_NB_FIXED, GaussianFunctionFactory.FIT_NB_FREE_CIRCLE);
	}

	@Test
	public void freeCircularComputesSameAsEllipticalDoublePeak()
	{
		f1ComputesSameAsf2(Double, GaussianFunctionFactory.FIT_FREE_CIRCLE, GaussianFunctionFactory.FIT_ELLIPTICAL);
	}

	@Test
	public void freeCircularFasterThanEllipticalDoublePeak()
	{
		f1FasterThanf2(Double, GaussianFunctionFactory.FIT_FREE_CIRCLE, GaussianFunctionFactory.FIT_ELLIPTICAL);
	}

	@Test
	public void circularComputesSameAsFreeCircularDoublePeak()
	{
		f1ComputesSameAsf2(Double, GaussianFunctionFactory.FIT_CIRCLE, GaussianFunctionFactory.FIT_FREE_CIRCLE);
	}

	@Test
	public void circularFasterThanFreeCircularDoublePeak()
	{
		f1FasterThanf2(Double, GaussianFunctionFactory.FIT_CIRCLE, GaussianFunctionFactory.FIT_FREE_CIRCLE);
	}

	@Test
	public void fixedComputesSameAsFreeCircularDoublePeak()
	{
		f1ComputesSameAsf2(Double, GaussianFunctionFactory.FIT_FIXED, GaussianFunctionFactory.FIT_FREE_CIRCLE);
	}

	@Test
	public void fixedFasterThanFreeCircularDoublePeak()
	{
		f1FasterThanf2(Double, GaussianFunctionFactory.FIT_FIXED, GaussianFunctionFactory.FIT_FREE_CIRCLE);
	}

	@Test
	public void freeCircularComputesSameAsEllipticalDoublePeakNB()
	{
		f1ComputesSameAsf2(Double, GaussianFunctionFactory.FIT_NB_FREE_CIRCLE,
				GaussianFunctionFactory.FIT_NB_ELLIPTICAL);
	}

	@Test
	public void freeCircularFasterThanEllipticalDoublePeakNB()
	{
		f1FasterThanf2(Double, GaussianFunctionFactory.FIT_NB_FREE_CIRCLE, GaussianFunctionFactory.FIT_NB_ELLIPTICAL);
	}

	@Test
	public void circularComputesSameAsFreeCircularDoublePeakNB()
	{
		f1ComputesSameAsf2(Double, GaussianFunctionFactory.FIT_NB_CIRCLE, GaussianFunctionFactory.FIT_NB_FREE_CIRCLE);
	}

	@Test
	public void circularFasterThanFreeCircularDoublePeakNB()
	{
		f1FasterThanf2(Double, GaussianFunctionFactory.FIT_NB_CIRCLE, GaussianFunctionFactory.FIT_NB_FREE_CIRCLE);
	}

	@Test
	public void fixedComputesSameAsFreeCircularDoublePeakNB()
	{
		f1ComputesSameAsf2(Double, GaussianFunctionFactory.FIT_NB_FIXED, GaussianFunctionFactory.FIT_NB_FREE_CIRCLE);
	}

	@Test
	public void fixedFasterThanFreeCircularDoublePeakNB()
	{
		f1FasterThanf2(Double, GaussianFunctionFactory.FIT_NB_FIXED, GaussianFunctionFactory.FIT_NB_FREE_CIRCLE);
	}

	void f1ComputesSameAsf2(int npeaks, int flags1, int flags2)
	{
		FloatEquality eq = new FloatEquality(2, 1e-10f);
		int iter = 2000;
		ArrayList<float[]> paramsList2 = (npeaks == 1) ? copyList(paramsListSinglePeak, iter) : copyList(
				paramsListDoublePeak, iter);

		Gaussian2DFunction f1 = GaussianFunctionFactory.create2D(1, blockWidth, flags1);
		Gaussian2DFunction f2 = GaussianFunctionFactory.create2D(1, blockWidth, flags2);

		float[] dyda1 = new float[1 + npeaks * 6];
		float[] dyda2 = new float[1 + npeaks * 6];

		int[] gradientIndices = f1.gradientIndices();
		int[] g1 = new int[gradientIndices.length];
		int[] g2 = new int[gradientIndices.length];
		int nparams = 0;
		for (int i = 0; i < gradientIndices.length; i++)
		{
			int index1 = f1.findGradientIndex(g1[i]);
			int index2 = f2.findGradientIndex(g2[i]);
			if (index1 >= 0 && index2 >= 0)
			{
				g1[nparams] = index1;
				g2[nparams] = index2;
				nparams++;
			}
		}

		for (int i = 0; i < paramsList2.size(); i++)
		{
			f1.initialise(paramsList2.get(i));
			f2.initialise(paramsList2.get(i));

			for (int j = 0; j < x.length; j++)
			{
				float y1 = f1.eval(x[j], dyda1);
				float y2 = f2.eval(x[j], dyda2);

				Assert.assertTrue("Not same y[" + j + "] @ " + i + " " + y1 + " != " + y2,
						eq.almostEqualComplement(y1, y2));

				for (int ii = 0; ii < nparams; ii++)
					Assert.assertTrue("Not same dyda[" + j + "] @ " + gradientIndices[g1[ii]] + ": " + dyda1[g1[ii]] +
							" != " + dyda2[g2[ii]], eq.almostEqualComplement(dyda1[g1[ii]], dyda2[g2[ii]]));
			}
		}
	}

	void f1FasterThanf2(int npeaks, int flags1, int flags2)
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);
		
		ArrayList<float[]> paramsList2 = (npeaks == 1) ? paramsListSinglePeak : paramsListDoublePeak;

		// Use the full list of parameters to build the functions
		GaussianFunction f1 = GaussianFunctionFactory.create2D(npeaks, blockWidth, flags1);
		GaussianFunction f2 = GaussianFunctionFactory.create2D(npeaks, blockWidth, flags2);

		float[] dyda = new float[1 + npeaks * 6];

		for (int i = 0; i < paramsList2.size(); i++)
		{
			f1.initialise(paramsList2.get(i));
			for (int j = 0; j < x.length; j++)
				f1.eval(x[j], dyda);
		}

		long start1 = System.nanoTime();
		for (int i = 0; i < paramsList2.size(); i++)
		{
			f1.initialise(paramsList2.get(i));
			for (int j = 0; j < x.length; j++)
				f1.eval(x[j], dyda);
		}
		start1 = System.nanoTime() - start1;

		for (int i = 0; i < paramsList2.size(); i++)
		{
			f2.initialise(paramsList2.get(i));
			for (int j = 0; j < x.length; j++)
				f2.eval(x[j], dyda);
		}

		long start2 = System.nanoTime();
		for (int i = 0; i < paramsList2.size(); i++)
		{
			f2.initialise(paramsList2.get(i));
			for (int j = 0; j < x.length; j++)
				f2.eval(x[j], dyda);
		}
		start2 = System.nanoTime() - start2;

		log("%s = %d : %s = %d : %fx\n", f1.getClass().getName(), start1, f2.getClass().getName(), start2,
				(1.0 * start2) / start1);
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(start2 > start1);
	}

	/**
	 * Create random elliptical Gaussian data an returns the data plus an estimate of the parameters.
	 * Only the chosen parameters are randomised and returned for a maximum of (background, amplitude, angle, xpos,
	 * ypos, xwidth, ywidth }
	 * 
	 * @param params
	 *            set on output
	 * @return
	 */
	private float[] floatCreateGaussianData(int npeaks, float[] params)
	{
		int n = blockWidth * blockWidth;

		// Generate a 2D Gaussian
		EllipticalGaussian2DFunction func = new EllipticalGaussian2DFunction(npeaks, blockWidth);
		params[0] = Background + rand.nextFloat() * 5f;
		for (int i = 0, j = 1; i < npeaks; i++, j += 6)
		{
			params[j] = Amplitude + rand.nextFloat() * 5f;
			params[j + 1] = 0f; //(float) (Math.PI / 4.0); // Angle
			params[j + 2] = Xpos + rand.nextFloat() * 2f;
			params[j + 3] = Ypos + rand.nextFloat() * 2f;
			params[j + 4] = Xwidth + rand.nextFloat() * 2f;
			params[j + 5] = params[j + 4];
		}

		float[] dy_da = new float[params.length];
		float[] y = new float[n];
		func.initialise(params);
		for (int i = 0; i < y.length; i++)
		{
			// Add random noise
			y[i] = func.eval(i, dy_da) + ((rand.nextFloat() < 0.5f) ? -rand.nextFloat() * 5f : rand.nextFloat() * 5f);
		}

		// Randomise only the necessary parameters (i.e. not angle and X & Y widths should be the same)
		params[0] += ((rand.nextFloat() < 0.5f) ? -rand.nextFloat() : rand.nextFloat());
		for (int i = 0, j = 1; i < npeaks; i++, j += 6)
		{
			params[j + 1] += ((rand.nextFloat() < 0.5f) ? -rand.nextFloat() : rand.nextFloat());
			params[j + 3] += ((rand.nextFloat() < 0.5f) ? -rand.nextFloat() : rand.nextFloat());
			params[j + 4] += ((rand.nextFloat() < 0.5f) ? -rand.nextFloat() : rand.nextFloat());
			params[j + 5] = params[j + 4];
		}

		return y;
	}

	protected int[] createData(int npeaks, int iter, ArrayList<float[]> paramsList, ArrayList<float[]> yList)
	{
		int[] x = new int[blockWidth * blockWidth];
		for (int i = 0; i < x.length; i++)
			x[i] = i;
		for (int i = 0; i < iter; i++)
		{
			float[] params = new float[1 + 6 * npeaks];
			float[] y = floatCreateGaussianData(npeaks, params);
			paramsList.add(params);
			yList.add(y);
		}
		return x;
	}

	protected ArrayList<float[]> copyList(ArrayList<float[]> paramsList, int iter)
	{
		iter = FastMath.min(iter, paramsList.size());

		ArrayList<float[]> params2List = new ArrayList<float[]>(iter);
		for (int i = 0; i < iter; i++)
		{
			params2List.add(paramsList.get(i));
		}
		return params2List;
	}

	void log(String format, Object... args)
	{
		System.out.printf(format, args);
	}
}
