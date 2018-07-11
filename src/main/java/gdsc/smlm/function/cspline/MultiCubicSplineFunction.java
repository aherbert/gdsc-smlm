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

import gdsc.core.utils.SimpleArrayUtils;
import gdsc.smlm.function.Gradient1Procedure;
import gdsc.smlm.function.Gradient2Procedure;
import gdsc.smlm.function.ValueProcedure;
import gdsc.smlm.results.PeakResult;

/**
 * Represent a cubic spline function for multiple points.
 */
public class MultiCubicSplineFunction extends CubicSplineFunction
{
	/** The number of splines to draw */
	private int n = 1;

	/** The gradient indices, This can be cached for the same n. */
	private int[] gradientIndices;

	/** The n target splines. This is cached to re-use memory */
	private TargetSpline[] t = new TargetSpline[0];

	/**
	 * The working splines for the current evaluation.
	 */
	private TargetSpline[] working;

	/**
	 * The number of working splines. This is <=n depending on whether the spline is within the target region for each
	 * point.
	 */
	private int w;

	/**
	 * The working splines for the current evaluation of the current y-index.
	 */
	private TargetSpline[] workingY;

	/**
	 * The number of working splines for the current y-index. This is <=w depending on whether the working spline is
	 * within the target region for Y..
	 */
	private int wY;

	/**
	 * Instantiates a new cubic spline function.
	 *
	 * @param splineData
	 *            the spline data
	 * @param maxx
	 *            The maximum x value of the 2-dimensional data
	 * @param maxy
	 *            The maximum y value of the 2-dimensional data
	 * @throws IllegalArgumentException
	 *             If the function does not have an integer grid spacing from the origin
	 */
	public MultiCubicSplineFunction(CubicSplineData splineData, int maxx, int maxy) throws IllegalArgumentException
	{
		super(splineData, maxx, maxy);
	}

	/**
	 * Instantiates a new cubic spline function.
	 *
	 * @param splineData
	 *            the spline data
	 * @param maxx
	 *            The maximum x value of the 2-dimensional data
	 * @param maxy
	 *            The maximum y value of the 2-dimensional data
	 * @param cx
	 *            the x centre of the spline data
	 * @param cy
	 *            the y centre of the spline data
	 * @param cz
	 *            the z centre of the spline data
	 * @param scale
	 *            the scale of the spline data
	 * @throws IllegalArgumentException
	 *             If the function does not have an integer grid spacing from the origin
	 */
	public MultiCubicSplineFunction(CubicSplineData splineData, int maxx, int maxy, double cx, double cy, double cz,
			int scale) throws IllegalArgumentException
	{
		super(splineData, maxx, maxy, cx, cy, cz, scale);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.function.cspline.CubicSplineFunction#getN()
	 */
	@Override
	public int getN()
	{
		return n;
	}

	/**
	 * Sets the number of splines to draw.
	 *
	 * @param n
	 *            the new number of splines to draw
	 * @throws IllegalArgumentException
	 *             If the number is not strictly positive
	 */
	public void setN(int n) throws IllegalArgumentException
	{
		if (n < 1)
			throw new IllegalArgumentException();
		if (n != this.n)
			gradientIndices = null;
		this.n = n;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.function.GradientFunction#gradientIndices()
	 */
	@Override
	public int[] gradientIndices()
	{
		if (gradientIndices == null)
			gradientIndices = SimpleArrayUtils.newArray(getNumberOfGradients(), 0, 1);
		return gradientIndices;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.function.GradientFunction#getNumberOfGradients()
	 */
	@Override
	public int getNumberOfGradients()
	{
		return 1 + 4 * n;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.function.cspline.CubicSplineFunction#initialise(double[], int)
	 */
	@Override
	protected void initialise(double[] a, int order)
	{
		tB = a[PeakResult.BACKGROUND];
		// Ensure we have enough room
		if (t.length < n)
		{
			int m = t.length;
			t = Arrays.copyOf(t, n); // Preserve memory space
			boolean sp = splines[0][0].isSinglePrecision();
			while (m < n)
				t[m++] = (sp) ? new FloatTargetSpline() : new DoubleTargetSpline();
			working = new TargetSpline[n];
			workingY = new TargetSpline[n];
		}
		// Convert the target parameters to spline offset tables
		w = 0;
		for (int i = 0, j = PeakResult.INTENSITY; i < n; i++)
		{
			double tI = a[j++];
			double tX = a[j++];
			double tY = a[j++];
			double tZ = a[j++];
			if (t[i].initialise(i, tI, tX, tY, tZ, order))
				working[w++] = t[i];
		}
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.function.ValueFunction#forEach(gdsc.smlm.function.ValueProcedure)
	 */
	@Override
	public void forEach(ValueProcedure procedure)
	{
		for (int n = 0; n < w; n++)
			working[n].reset();

		for (int y = 0; y < maxy; y++)
		{
			// Get the working targets for this Y
			wY = 0;
			for (int n = 0; n < w; n++)
			{
				if (working[n].isNextYActive())
					workingY[wY++] = working[n];
			}

			if (wY == 0)
			{
				for (int x = 0; x < maxx; x++)
					procedure.execute(tB);
			}
			else
			{
				for (int x = 0; x < maxx; x++)
				{
					double I = tB;
					for (int n = 0; n < wY; n++)
					{
						I += workingY[n].value(x);
					}
					procedure.execute(I);
				}
			}
		}
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.function.Gradient1Function#forEach(gdsc.smlm.function.Gradient1Procedure)
	 */
	@Override
	public void forEach(Gradient1Procedure procedure)
	{
		final double[] duda = new double[getNumberOfGradients()];
		duda[0] = 1.0;

		for (int n = 0; n < w; n++)
			working[n].reset();

		for (int y = 0; y < maxy; y++)
		{
			// Get the working targets for this Y
			wY = 0;
			for (int n = 0; n < w; n++)
			{
				if (working[n].isNextYActive(duda))
					workingY[wY++] = working[n];
			}

			if (wY == 0)
			{
				for (int x = 0; x < maxx; x++)
					procedure.execute(tB, duda);
			}
			else
			{
				for (int x = 0; x < maxx; x++)
				{
					double I = tB;
					for (int n = 0; n < wY; n++)
					{
						I += workingY[n].value(x, duda);
					}
					procedure.execute(I, duda);
				}
			}
		}
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.function.Gradient2Function#forEach(gdsc.smlm.function.Gradient2Procedure)
	 */
	@Override
	public void forEach(Gradient2Procedure procedure)
	{
		final double[] duda = new double[getNumberOfGradients()];
		final double[] d2uda2 = new double[getNumberOfGradients()];
		duda[0] = 1.0;

		for (int n = 0; n < w; n++)
			working[n].reset();

		for (int y = 0; y < maxy; y++)
		{
			// Get the working targets for this Y
			wY = 0;
			for (int n = 0; n < w; n++)
			{
				if (working[n].isNextYActive(duda, d2uda2))
					workingY[wY++] = working[n];
			}

			if (wY == 0)
			{
				for (int x = 0; x < maxx; x++)
					procedure.execute(tB, duda, d2uda2);
			}
			else
			{
				for (int x = 0; x < maxx; x++)
				{
					double I = tB;
					for (int n = 0; n < wY; n++)
					{
						I += workingY[n].value(x, duda, d2uda2);
					}
					procedure.execute(I, duda, d2uda2);
				}
			}
		}
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.function.cspline.CubicSplineFunction#isNodeBoundary(int)
	 */
	@Override
	public boolean isNodeBoundary(int gradientIndex)
	{
		int parameterIndex = gradientIndices()[gradientIndex];
		if (parameterIndex == BACKGROUND)
			return false;

		int dimension = (parameterIndex - 1) % PARAMETERS_PER_PEAK;
		if (dimension == 0)
			return false; // Signal

		int peak = getPeak(parameterIndex);
		for (int n = 0; n < w; n++)
			if (working[n].id == peak)
				return working[n].isNodeBoundary(dimension - 1);
		return false;
	}
}
