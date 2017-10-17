package gdsc.smlm.function.cspline;

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
 * Create a cubic spline function.
 */
public class CubicSplineFunctionFactory
{
	/**
	 * Instantiates a new cubic spline function to model n points.
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
	 * @param n
	 *            the number of points
	 * @return the cubic spline function
	 * @throws IllegalArgumentException
	 *             If the function does not have an integer grid spacing from the origin
	 */
	public static CubicSplineFunction createCubicSplineFunction(CubicSplineData splineData, int maxx, int maxy,
			double cx, double cy, double cz, int scale, int n) throws IllegalArgumentException
	{
		if (n == 1)
			return new SingleCubicSplineFunction(splineData, maxx, maxy, cx, cy, cz, scale);
		MultiCubicSplineFunction f = new MultiCubicSplineFunction(splineData, maxx, maxy, cx, cy, cz, scale);
		f.setN(n);
		return f;
	}
}