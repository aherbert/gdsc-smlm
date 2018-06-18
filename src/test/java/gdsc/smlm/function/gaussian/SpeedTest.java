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

import java.util.ArrayList;
import java.util.Random;

import org.apache.commons.math3.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.utils.DoubleEquality;
import gdsc.test.TestSettings;

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
	private double Background = 20;
	private double Amplitude = 10;
	private double Xpos = 5;
	private double Ypos = 5;
	private double Xwidth = 5;
	private Random rand;
	private static ArrayList<double[]> paramsListSinglePeak = null;
	private static ArrayList<double[]> yListSinglePeak;
	private static int[] x;
	private static ArrayList<double[]> paramsListDoublePeak;
	private static ArrayList<double[]> yListDoublePeak;

	public SpeedTest()
	{
		// TODO - the data generation could be static in the Gaussian2DFunctionTest since 
		// all gaussian data generation should be similar

		rand = new Random(30051977);

		if (paramsListSinglePeak == null)
		{
			paramsListSinglePeak = new ArrayList<double[]>(MAX_ITER);
			yListSinglePeak = new ArrayList<double[]>(MAX_ITER);
			x = createData(1, MAX_ITER, paramsListSinglePeak, yListSinglePeak);
			paramsListDoublePeak = new ArrayList<double[]>(MAX_ITER);
			yListDoublePeak = new ArrayList<double[]>(MAX_ITER);
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
		f1ComputesSameAsf2(Single, GaussianFunctionFactory.FIT_SIMPLE_NB_FREE_CIRCLE,
				GaussianFunctionFactory.FIT_SIMPLE_NB_ELLIPTICAL);
	}

	@Test
	public void freeCircularFasterThanEllipticalSinglePeakNB()
	{
		f1FasterThanf2(Single, GaussianFunctionFactory.FIT_SIMPLE_NB_FREE_CIRCLE,
				GaussianFunctionFactory.FIT_SIMPLE_NB_ELLIPTICAL);
	}

	@Test
	public void circularComputesSameAsFreeCircularSinglePeakNB()
	{
		f1ComputesSameAsf2(Single, GaussianFunctionFactory.FIT_SIMPLE_NB_CIRCLE,
				GaussianFunctionFactory.FIT_SIMPLE_NB_FREE_CIRCLE);
	}

	@Test
	public void circularFasterThanFreeCircularSinglePeakNB()
	{
		f1FasterThanf2(Single, GaussianFunctionFactory.FIT_SIMPLE_NB_CIRCLE,
				GaussianFunctionFactory.FIT_SIMPLE_NB_FREE_CIRCLE);
	}

	@Test
	public void fixedComputesSameAsFreeCircularSinglePeakNB()
	{
		f1ComputesSameAsf2(Single, GaussianFunctionFactory.FIT_SIMPLE_NB_FIXED,
				GaussianFunctionFactory.FIT_SIMPLE_NB_FREE_CIRCLE);
	}

	@Test
	public void fixedFasterThanFreeCircularSinglePeakNB()
	{
		f1FasterThanf2(Single, GaussianFunctionFactory.FIT_SIMPLE_NB_FIXED,
				GaussianFunctionFactory.FIT_SIMPLE_NB_FREE_CIRCLE);
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
		f1ComputesSameAsf2(Double, GaussianFunctionFactory.FIT_SIMPLE_NB_FREE_CIRCLE,
				GaussianFunctionFactory.FIT_SIMPLE_NB_ELLIPTICAL);
	}

	@Test
	public void freeCircularFasterThanEllipticalDoublePeakNB()
	{
		f1FasterThanf2(Double, GaussianFunctionFactory.FIT_SIMPLE_NB_FREE_CIRCLE,
				GaussianFunctionFactory.FIT_SIMPLE_NB_ELLIPTICAL);
	}

	@Test
	public void circularComputesSameAsFreeCircularDoublePeakNB()
	{
		f1ComputesSameAsf2(Double, GaussianFunctionFactory.FIT_SIMPLE_NB_CIRCLE,
				GaussianFunctionFactory.FIT_SIMPLE_NB_FREE_CIRCLE);
	}

	@Test
	public void circularFasterThanFreeCircularDoublePeakNB()
	{
		f1FasterThanf2(Double, GaussianFunctionFactory.FIT_SIMPLE_NB_CIRCLE,
				GaussianFunctionFactory.FIT_SIMPLE_NB_FREE_CIRCLE);
	}

	@Test
	public void fixedComputesSameAsFreeCircularDoublePeakNB()
	{
		f1ComputesSameAsf2(Double, GaussianFunctionFactory.FIT_SIMPLE_NB_FIXED,
				GaussianFunctionFactory.FIT_SIMPLE_NB_FREE_CIRCLE);
	}

	@Test
	public void fixedFasterThanFreeCircularDoublePeakNB()
	{
		f1FasterThanf2(Double, GaussianFunctionFactory.FIT_SIMPLE_NB_FIXED,
				GaussianFunctionFactory.FIT_SIMPLE_NB_FREE_CIRCLE);
	}

	void f1ComputesSameAsf2(int npeaks, int flags1, int flags2)
	{
		DoubleEquality eq = new DoubleEquality(1e-2, 1e-10);
		int iter = 2000;
		ArrayList<double[]> paramsList2 = (npeaks == 1) ? copyList(paramsListSinglePeak, iter)
				: copyList(paramsListDoublePeak, iter);

		Gaussian2DFunction f1 = GaussianFunctionFactory.create2D(1, blockWidth, blockWidth, flags1, null);
		Gaussian2DFunction f2 = GaussianFunctionFactory.create2D(1, blockWidth, blockWidth, flags2, null);

		double[] dyda1 = new double[1 + npeaks * Gaussian2DFunction.PARAMETERS_PER_PEAK];
		double[] dyda2 = new double[1 + npeaks * Gaussian2DFunction.PARAMETERS_PER_PEAK];

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
				double y1 = f1.eval(x[j], dyda1);
				double y2 = f2.eval(x[j], dyda2);

				Assert.assertTrue("Not same y[" + j + "] @ " + i + " " + y1 + " != " + y2,
						eq.almostEqualRelativeOrAbsolute(y1, y2));

				for (int ii = 0; ii < nparams; ii++)
					Assert.assertTrue("Not same dyda[" + j + "] @ " + gradientIndices[g1[ii]] + ": " + dyda1[g1[ii]] +
							" != " + dyda2[g2[ii]], eq.almostEqualRelativeOrAbsolute(dyda1[g1[ii]], dyda2[g2[ii]]));
			}
		}
	}

	void f1FasterThanf2(int npeaks, int flags1, int flags2)
	{
		TestSettings.assumeMediumComplexity();

		ArrayList<double[]> paramsList2 = (npeaks == 1) ? paramsListSinglePeak : paramsListDoublePeak;

		// Use the full list of parameters to build the functions
		Gaussian2DFunction f1 = GaussianFunctionFactory.create2D(npeaks, blockWidth, blockWidth, flags1, null);
		Gaussian2DFunction f2 = GaussianFunctionFactory.create2D(npeaks, blockWidth, blockWidth, flags2, null);

		double[] dyda = new double[1 + npeaks * 6];

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
	private double[] doubleCreateGaussianData(int npeaks, double[] params)
	{
		int n = blockWidth * blockWidth;

		// Generate a 2D Gaussian
		EllipticalGaussian2DFunction func = new EllipticalGaussian2DFunction(npeaks, blockWidth, blockWidth);
		params[0] = Background + rand.nextFloat() * 5f;
		for (int i = 0, j = 0; i < npeaks; i++, j += Gaussian2DFunction.PARAMETERS_PER_PEAK)
		{
			params[j] = Amplitude + rand.nextFloat() * 5f;
			params[j + Gaussian2DFunction.X_POSITION] = Xpos + rand.nextFloat() * 2f;
			params[j + Gaussian2DFunction.Y_POSITION] = Ypos + rand.nextFloat() * 2f;
			params[j + Gaussian2DFunction.X_SD] = Xwidth + rand.nextFloat() * 2f;
			params[j + Gaussian2DFunction.Y_SD] = params[j + 4];
			params[j + Gaussian2DFunction.ANGLE] = 0f; //(double) (Math.PI / 4.0); // Angle
		}

		double[] dy_da = new double[params.length];
		double[] y = new double[n];
		func.initialise(params);
		for (int i = 0; i < y.length; i++)
		{
			// Add random noise
			y[i] = func.eval(i, dy_da) + ((rand.nextFloat() < 0.5f) ? -rand.nextFloat() * 5f : rand.nextFloat() * 5f);
		}

		// Randomise only the necessary parameters (i.e. not angle and X & Y widths should be the same)
		params[0] += ((rand.nextFloat() < 0.5f) ? -rand.nextFloat() : rand.nextFloat());
		for (int i = 0, j = 0; i < npeaks; i++, j += Gaussian2DFunction.PARAMETERS_PER_PEAK)
		{
			params[j + Gaussian2DFunction.X_POSITION] += ((rand.nextFloat() < 0.5f) ? -rand.nextFloat()
					: rand.nextFloat());
			params[j + Gaussian2DFunction.Y_POSITION] += ((rand.nextFloat() < 0.5f) ? -rand.nextFloat()
					: rand.nextFloat());
			params[j + Gaussian2DFunction.X_SD] += ((rand.nextFloat() < 0.5f) ? -rand.nextFloat() : rand.nextFloat());
			params[j + Gaussian2DFunction.Y_SD] = params[j + Gaussian2DFunction.X_SD];
		}

		return y;
	}

	protected int[] createData(int npeaks, int iter, ArrayList<double[]> paramsList, ArrayList<double[]> yList)
	{
		int[] x = new int[blockWidth * blockWidth];
		for (int i = 0; i < x.length; i++)
			x[i] = i;
		for (int i = 0; i < iter; i++)
		{
			double[] params = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK * npeaks];
			double[] y = doubleCreateGaussianData(npeaks, params);
			paramsList.add(params);
			yList.add(y);
		}
		return x;
	}

	protected ArrayList<double[]> copyList(ArrayList<double[]> paramsList, int iter)
	{
		iter = FastMath.min(iter, paramsList.size());

		ArrayList<double[]> params2List = new ArrayList<double[]>(iter);
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
