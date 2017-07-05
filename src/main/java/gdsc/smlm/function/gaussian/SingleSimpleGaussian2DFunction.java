package gdsc.smlm.function.gaussian;

import org.apache.commons.math3.util.FastMath;

import gdsc.core.utils.NotImplementedException;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2013 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Evaluates an 2-dimensional Gaussian function for a single peak.
 * <p>
 * The single parameter x in the {@link #eval(int, double[])} function is assumed to be a linear index into
 * 2-dimensional
 * data. The dimensions of the data must be specified to allow unpacking to coordinates.
 * <p>
 * Data should be packed in descending dimension order, e.g. Y,X : Index for [x,y] = MaxX*y + x.
 * <p>
 * This function uses the signal, position and x-sd parameters. It ignores the background parameter, y-sd and angle
 * parameter. It can be used for fast evaluation of the sum of a Gaussian within a region defined by maxx. maxy is
 * defined by the number of calls to the eval() functions. It does not support gradient evaluation.
 */
public class SingleSimpleGaussian2DFunction extends Gaussian2DFunction
{
	private static final int[] gradientIndices = new int[0];

	protected double x0pos;
	protected double x1pos;

	protected double height;
	protected double aa;

	/**
	 * Constructor
	 * 
	 * @param maxx
	 *            The maximum x value of the 2-dimensional data (used to unpack a linear index into coordinates)
	 * @param maxy
	 *            The maximum y value of the 2-dimensional data (used to unpack a linear index into coordinates)
	 */
	public SingleSimpleGaussian2DFunction(int maxx, int maxy)
	{
		super(maxx, maxy);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.gaussian.Gaussian2DFunction#copy()
	 */
	@Override
	public Gaussian2DFunction copy()
	{
		return new SingleSimpleGaussian2DFunction(maxx, maxy);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.fitting.function.NonLinearFunction#initialise(double[])
	 */
	public void initialise(double[] a)
	{
		x0pos = a[X_POSITION];
		x1pos = a[Y_POSITION];

		final double sx = a[X_SD];
		final double sx2 = sx * sx;

		final double n = ONE_OVER_TWO_PI / sx2;
		height = a[SIGNAL] * n;

		// All prefactors are negated since the Gaussian uses the exponential to the negative:
		// A * exp( -( a(x-x0)^2 + 2b(x-x0)(y-y0) + c(y-y0)^2 ) )

		aa = -0.5 / sx2;
	}

	/**
	 * Not implemented.
	 * <p>
	 * {@inheritDoc}
	 * @see gdsc.smlm.function.gaussian.Gaussian2DFunction#eval(int, double[])
	 */
	public double eval(final int x, final double[] dyda)
	{
		throw new NotImplementedException();
	}

	/**
	 * Evaluates an 2-dimensional circular Gaussian function for a single peak.
	 * <p>
	 * {@inheritDoc}
	 * @see gdsc.smlm.function.gaussian.Gaussian2DFunction#eval(int, double[])
	 */
	public double eval(final int x)
	{
		// Unpack the predictor into the dimensions
		final int x1 = x / maxx;
		final int x0 = x % maxx;

		final double dx = x0 - x0pos;
		final double dy = x1 - x1pos;

		return height * FastMath.exp(aa * (dx * dx + dy * dy));
	}

	@Override
	public int getNPeaks()
	{
		return 1;
	}

	@Override
	public boolean evaluatesBackground()
	{
		return false;
	}

	@Override
	public boolean evaluatesSignal()
	{
		return false;
	}

	@Override
	public boolean evaluatesAngle()
	{
		return false;
	}

	@Override
	public boolean evaluatesPosition()
	{
		return false;
	}

	@Override
	public boolean evaluatesSD0()
	{
		return false;
	}

	@Override
	public boolean evaluatesSD1()
	{
		return false;
	}

	@Override
	public int getGradientParametersPerPeak()
	{
		return 0;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.fitting.function.NonLinearFunction#gradientIndices()
	 */
	public int[] gradientIndices()
	{
		return gradientIndices;
	}
}
