package gdsc.smlm.function.cspline;

import gdsc.core.math.interpolation.CustomTricubicInterpolatingFunction;
import gdsc.core.utils.SimpleArrayUtils;
import gdsc.smlm.function.Gradient1Procedure;
import gdsc.smlm.function.Gradient2Function;
import gdsc.smlm.function.Gradient2Procedure;
import gdsc.smlm.function.ValueProcedure;

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
 * Represent a cubic spline function.
 */
public class CubicSplineFunction implements Gradient2Function
{
	/** The cubic spline function. */
	private final CustomTricubicInterpolatingFunction function;

	/** The scale. */
	private int scale = 1;

	/** The cx. */
	private double cx;

	/** The cy. */
	private double cy;

	/** The cz. */
	private double cz;

	private final int maxx, maxy;

	private int n = 1;

	private int[] gradientIndices;

	/**
	 * Instantiates a new cubic spline function.
	 *
	 * @param function
	 *            the function
	 * @param maxx
	 *            The maximum x value of the 2-dimensional data
	 * @param maxy
	 *            The maximum y value of the 2-dimensional data
	 * @throws IllegalArgumentException
	 *             If the function does not have an integer grid spacing
	 */
	public CubicSplineFunction(CustomTricubicInterpolatingFunction function, int maxx, int maxy)
			throws IllegalArgumentException
	{
		if (function == null || !function.isInteger())
			throw new IllegalArgumentException("Require integer grid spacing for the cubic spline");
		if (n < 1)
			throw new IllegalArgumentException("Number of splines must be strictly positive");
		this.function = function;
		this.maxx = (maxx < 1) ? 1 : maxx;
		this.maxy = (maxy < 1) ? 1 : maxy;
		// Centre in the middle
		cx = (function.getMaxX() - function.getMinX()) / 2.0;
		cy = (function.getMaxY() - function.getMinY()) / 2.0;
		cz = (function.getMaxZ() - function.getMinZ()) / 2.0;
	}

	/**
	 * Gets the scale.
	 *
	 * @return the scale
	 */
	public int getScale()
	{
		return scale;
	}

	/**
	 * Sets the scale to map the cubic spline function to the integer grid. E.g. set a scale of 2 to render the spline
	 * at half its size.
	 *
	 * @param scale
	 *            the new scale
	 * @throws IllegalArgumentException
	 *             If the scale is not strictly positive
	 */
	public void setScale(int scale) throws IllegalArgumentException
	{
		if (scale < 1)
			throw new IllegalArgumentException();
		this.scale = scale;
	}

	/**
	 * Gets the number of splines to draw.
	 *
	 * @return the number of splines to draw
	 */
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

	/**
	 * Gets the centre X.
	 *
	 * @return the centre X
	 */
	public double getCentreX()
	{
		return cx;
	}

	/**
	 * Sets the centre X.
	 *
	 * @param cx
	 *            the new centre X
	 */
	public void setCentreX(double cx)
	{
		this.cx = cx;
	}

	/**
	 * Gets the centre Y.
	 *
	 * @return the centre Y
	 */
	public double getCentreY()
	{
		return cy;
	}

	/**
	 * Sets the centre Y.
	 *
	 * @param cy
	 *            the new centre Y
	 */
	public void setCentreY(double cy)
	{
		this.cy = cy;
	}

	/**
	 * Gets the centre Z.
	 *
	 * @return the centre Z
	 */
	public double getCentreZ()
	{
		return cz;
	}

	/**
	 * Sets the centre Z.
	 *
	 * @param cz
	 *            the new centre Z
	 */
	public void setCentreZ(double cz)
	{
		this.cz = cz;
	}

	/**
	 * Gets the function.
	 *
	 * @return the function
	 */
	public CustomTricubicInterpolatingFunction getFunction()
	{
		return function;
	}

	public void initialise0(double[] a)
	{
		// TODO Auto-generated method stub

	}

	public void forEach(ValueProcedure procedure)
	{
		// TODO Auto-generated method stub

	}

	public void initialise1(double[] a)
	{
		// TODO Auto-generated method stub

	}

	public void forEach(Gradient1Procedure procedure)
	{
		// TODO Auto-generated method stub

	}

	public int size()
	{
		return maxx * maxy;
	}

	public void initialise(double[] a)
	{
		initialise1(a);
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

	public void initialise2(double[] a)
	{
		// TODO Auto-generated method stub

	}

	public void forEach(Gradient2Procedure procedure)
	{
		// TODO Auto-generated method stub

	}
}