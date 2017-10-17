package gdsc.smlm.function.cspline;

import java.util.Arrays;

import gdsc.core.utils.SimpleArrayUtils;
import gdsc.smlm.function.Gradient1Procedure;
import gdsc.smlm.function.Gradient2Procedure;
import gdsc.smlm.function.ValueProcedure;
import gdsc.smlm.results.PeakResult;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2017 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Represent a cubic spline function for multiple points
 */
public class MultiCubicSplineFunction extends CubicSplineFunction
{
	private int n = 1;
	private int[] gradientIndices;

	private TargetSpline[] t = new TargetSpline[0];

	private TargetSpline[] working;
	private int w;

	private TargetSpline[] workingY;
	private int wy;

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

	public int[] gradientIndices()
	{
		if (gradientIndices == null)
			gradientIndices = SimpleArrayUtils.newIntArray(getNumberOfGradients(), 0);
		return gradientIndices;
	}

	public int getNumberOfGradients()
	{
		return 1 + 4 * n;
	}

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
			if (t[i].initialise(tI, tX, tY, tZ, order))
				working[w++] = t[i];
		}
	}

	public void forEach(ValueProcedure procedure)
	{
		for (int n = 0; n < w; n++)
			working[w].reset();

		for (int y = 0; y < maxy; y++)
		{
			// Get the working targets for this Y 
			wy = 0;
			for (int n = 0; n < w; n++)
			{
				if (working[w].isNextYActive())
					workingY[wy++] = working[w];
			}

			if (wy == 0)
			{
				for (int x = 0; x < maxx; x++)
					procedure.execute(tB);
			}
			else
			{
				for (int x = 0; x < maxx; x++)
				{
					double I = tB;
					for (int n = 0; n < wy; n++)
					{
						I += workingY[wy].value(x);
					}
					procedure.execute(I);
				}
			}
		}
	}

	public void forEach(Gradient1Procedure procedure)
	{

	}

	public void forEach(Gradient2Procedure procedure)
	{

	}
}