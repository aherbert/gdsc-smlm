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

import java.util.ArrayList;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.test.TestLog;
import uk.ac.sussex.gdsc.test.TestSettings;

/**
 * Contains speed tests for the fastest method for calculating the Hessian and gradient vector
 * from a Gaussian 2D Function
 */
@SuppressWarnings({ "javadoc" })
public class SpeedTest
{
	private final int Single = 1;
	private final int Multi = 2;

	private static int blockWidth = 10;
	private static double Background = 20;
	private static double Amplitude = 10;
	private static double Xpos = 5;
	private static double Ypos = 5;
	private static double Xwidth = 5;

	private static RandomGenerator rand = TestSettings.getRandomGenerator();

	private static ArrayList<double[]> paramsListSinglePeak = new ArrayList<>();
	private static ArrayList<double[]> yListSinglePeak = new ArrayList<>();
	private static ArrayList<double[]> paramsListMultiPeak = new ArrayList<>();
	private static ArrayList<double[]> yListMultiPeak = new ArrayList<>();

	private static int[] x;
	static
	{
		x = new int[blockWidth * blockWidth];
		for (int i = 0; i < x.length; i++)
			x[i] = i;
	}

	private static synchronized void ensureDataSingle(int size)
	{
		if (paramsListSinglePeak.size() < size)
			createData(1, size, paramsListSinglePeak, yListSinglePeak);
	}

	private static synchronized void ensureDataMulti(int size)
	{
		if (paramsListMultiPeak.size() < size)
			createData(2, size, paramsListMultiPeak, yListMultiPeak);
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
	public void freeCircularComputesSameAsEllipticalMultiPeak()
	{
		f1ComputesSameAsf2(Multi, GaussianFunctionFactory.FIT_FREE_CIRCLE, GaussianFunctionFactory.FIT_ELLIPTICAL);
	}

	@Test
	public void freeCircularFasterThanEllipticalMultiPeak()
	{
		f1FasterThanf2(Multi, GaussianFunctionFactory.FIT_FREE_CIRCLE, GaussianFunctionFactory.FIT_ELLIPTICAL);
	}

	@Test
	public void circularComputesSameAsFreeCircularMultiPeak()
	{
		f1ComputesSameAsf2(Multi, GaussianFunctionFactory.FIT_CIRCLE, GaussianFunctionFactory.FIT_FREE_CIRCLE);
	}

	@Test
	public void circularFasterThanFreeCircularMultiPeak()
	{
		f1FasterThanf2(Multi, GaussianFunctionFactory.FIT_CIRCLE, GaussianFunctionFactory.FIT_FREE_CIRCLE);
	}

	@Test
	public void fixedComputesSameAsFreeCircularMultiPeak()
	{
		f1ComputesSameAsf2(Multi, GaussianFunctionFactory.FIT_FIXED, GaussianFunctionFactory.FIT_FREE_CIRCLE);
	}

	@Test
	public void fixedFasterThanFreeCircularMultiPeak()
	{
		f1FasterThanf2(Multi, GaussianFunctionFactory.FIT_FIXED, GaussianFunctionFactory.FIT_FREE_CIRCLE);
	}

	@Test
	public void freeCircularComputesSameAsEllipticalMultiPeakNB()
	{
		f1ComputesSameAsf2(Multi, GaussianFunctionFactory.FIT_SIMPLE_NB_FREE_CIRCLE,
				GaussianFunctionFactory.FIT_SIMPLE_NB_ELLIPTICAL);
	}

	@Test
	public void freeCircularFasterThanEllipticalMultiPeakNB()
	{
		f1FasterThanf2(Multi, GaussianFunctionFactory.FIT_SIMPLE_NB_FREE_CIRCLE,
				GaussianFunctionFactory.FIT_SIMPLE_NB_ELLIPTICAL);
	}

	@Test
	public void circularComputesSameAsFreeCircularMultiPeakNB()
	{
		f1ComputesSameAsf2(Multi, GaussianFunctionFactory.FIT_SIMPLE_NB_CIRCLE,
				GaussianFunctionFactory.FIT_SIMPLE_NB_FREE_CIRCLE);
	}

	@Test
	public void circularFasterThanFreeCircularMultiPeakNB()
	{
		f1FasterThanf2(Multi, GaussianFunctionFactory.FIT_SIMPLE_NB_CIRCLE,
				GaussianFunctionFactory.FIT_SIMPLE_NB_FREE_CIRCLE);
	}

	@Test
	public void fixedComputesSameAsFreeCircularMultiPeakNB()
	{
		f1ComputesSameAsf2(Multi, GaussianFunctionFactory.FIT_SIMPLE_NB_FIXED,
				GaussianFunctionFactory.FIT_SIMPLE_NB_FREE_CIRCLE);
	}

	@Test
	public void fixedFasterThanFreeCircularMultiPeakNB()
	{
		f1FasterThanf2(Multi, GaussianFunctionFactory.FIT_SIMPLE_NB_FIXED,
				GaussianFunctionFactory.FIT_SIMPLE_NB_FREE_CIRCLE);
	}

	void f1ComputesSameAsf2(int npeaks, int flags1, int flags2)
	{
		final DoubleEquality eq = new DoubleEquality(1e-2, 1e-10);
		final int iter = 50;
		ArrayList<double[]> paramsList2;
		if (npeaks == 1)
		{
			ensureDataSingle(iter);
			paramsList2 = copyList(paramsListSinglePeak, iter);
		}
		else
		{
			ensureDataMulti(iter);
			paramsList2 = copyList(paramsListMultiPeak, iter);
		}

		final Gaussian2DFunction f1 = GaussianFunctionFactory.create2D(1, blockWidth, blockWidth, flags1, null);
		final Gaussian2DFunction f2 = GaussianFunctionFactory.create2D(1, blockWidth, blockWidth, flags2, null);

		final double[] dyda1 = new double[1 + npeaks * Gaussian2DFunction.PARAMETERS_PER_PEAK];
		final double[] dyda2 = new double[1 + npeaks * Gaussian2DFunction.PARAMETERS_PER_PEAK];

		final int[] gradientIndices = f1.gradientIndices();
		final int[] g1 = new int[gradientIndices.length];
		final int[] g2 = new int[gradientIndices.length];
		int nparams = 0;
		for (int i = 0; i < gradientIndices.length; i++)
		{
			final int index1 = f1.findGradientIndex(g1[i]);
			final int index2 = f2.findGradientIndex(g2[i]);
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
				final double y1 = f1.eval(x[j], dyda1);
				final double y2 = f2.eval(x[j], dyda2);

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
		TestSettings.assumeSpeedTest();

		final int iter = 10000;
		ArrayList<double[]> paramsList2;
		if (npeaks == 1)
		{
			ensureDataSingle(iter);
			paramsList2 = paramsListSinglePeak;
		}
		else
		{
			ensureDataMulti(iter);
			paramsList2 = paramsListMultiPeak;
		}

		// Use the full list of parameters to build the functions
		final Gaussian2DFunction f1 = GaussianFunctionFactory.create2D(npeaks, blockWidth, blockWidth, flags1, null);
		final Gaussian2DFunction f2 = GaussianFunctionFactory.create2D(npeaks, blockWidth, blockWidth, flags2, null);

		final double[] dyda = new double[1 + npeaks * 6];

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

		TestLog.logSpeedTestResult(start2 > start1, "%s = %d : %s = %d : %fx\n", f1.getClass().getName(), start1,
				f2.getClass().getName(), start2, (1.0 * start2) / start1);
	}

	/**
	 * Create random elliptical Gaussian data an returns the data plus an estimate of the parameters.
	 * Only the chosen parameters are randomised and returned for a maximum of (background, amplitude, angle, xpos,
	 * ypos, xwidth, ywidth }
	 *
	 * @param npeaks
	 *            the npeaks
	 * @param params
	 *            set on output
	 * @return the data
	 */
	private static double[] doubleCreateGaussianData(int npeaks, double[] params)
	{
		final int n = blockWidth * blockWidth;

		// Generate a 2D Gaussian
		final EllipticalGaussian2DFunction func = new EllipticalGaussian2DFunction(npeaks, blockWidth, blockWidth);
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

		final double[] dy_da = new double[params.length];
		final double[] y = new double[n];
		func.initialise(params);
		for (int i = 0; i < y.length; i++)
			// Add random noise
			y[i] = func.eval(i, dy_da) + ((rand.nextFloat() < 0.5f) ? -rand.nextFloat() * 5f : rand.nextFloat() * 5f);

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

	protected static void createData(int npeaks, int iter, ArrayList<double[]> paramsList, ArrayList<double[]> yList)
	{
		for (int i = 0; i < iter; i++)
		{
			final double[] params = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK * npeaks];
			final double[] y = doubleCreateGaussianData(npeaks, params);
			paramsList.add(params);
			yList.add(y);
		}
	}

	protected ArrayList<double[]> copyList(ArrayList<double[]> paramsList, int iter)
	{
		iter = FastMath.min(iter, paramsList.size());

		final ArrayList<double[]> params2List = new ArrayList<>(iter);
		for (int i = 0; i < iter; i++)
			params2List.add(paramsList.get(i));
		return params2List;
	}

	void log(String format, Object... args)
	{
		System.out.printf(format, args);
	}
}
