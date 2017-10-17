package gdsc.smlm.function.cspline;

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
 * Represent a cubic spline function for a single points
 */
public class SingleCubicSplineFunction extends CubicSplineFunction
{
	private static final int[] gradientIndices = new int[] { 0, 1, 2, 3, 4 };

	private TargetSpline t;

	private TargetSpline working;

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
	public SingleCubicSplineFunction(CubicSplineData splineData, int maxx, int maxy) throws IllegalArgumentException
	{
		super(splineData, maxx, maxy);
		t = (splines[0][0].isSinglePrecision()) ? new FloatTargetSpline() : new DoubleTargetSpline();
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
	public SingleCubicSplineFunction(CubicSplineData splineData, int maxx, int maxy, double cx, double cy, double cz,
			int scale) throws IllegalArgumentException
	{
		super(splineData, maxx, maxy, cx, cy, cz, scale);
		t = (splines[0][0].isSinglePrecision()) ? new FloatTargetSpline() : new DoubleTargetSpline();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.cspline.CubicSplineFunction#getN()
	 */
	@Override
	public int getN()
	{
		return 1;
	}

	public int[] gradientIndices()
	{
		return gradientIndices;
	}

	public int getNumberOfGradients()
	{
		return 5;
	}

	protected void initialise(double[] a, int order)
	{
		tB = a[PeakResult.BACKGROUND];
		working = (t.initialise(a[PeakResult.INTENSITY], a[PeakResult.X], a[PeakResult.Y], a[PeakResult.Z], order)) ? t
				: null;
	}

	public void forEach(ValueProcedure procedure)
	{
		if (working == null)
		{
			// Special case
			for (int i = 0, size = maxx * maxy; i < size; i++)
				procedure.execute(tB);
		}
		else
		{
			working.reset();
			for (int y = 0; y < maxy; y++)
			{
				if (working.isNextYActive())
				{
					for (int x = 0; x < maxx; x++)
						procedure.execute(tB);
				}
				else
				{
					for (int x = 0; x < maxx; x++)
					{
						procedure.execute(tB + working.value(x));
					}
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