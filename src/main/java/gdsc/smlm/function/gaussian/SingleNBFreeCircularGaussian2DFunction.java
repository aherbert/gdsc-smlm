package gdsc.smlm.function.gaussian;

import org.apache.commons.math3.util.FastMath;

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
 */
public class SingleNBFreeCircularGaussian2DFunction extends SingleFreeCircularGaussian2DFunction
{
	private static final int[] gradientIndices;
	static
	{
		gradientIndices = createGradientIndices(1, new SingleNBFreeCircularGaussian2DFunction(1, 1));
	}

	/**
	 * Constructor
	 * 
	 * @param maxx
	 *            The maximum x value of the 2-dimensional data (used to unpack a linear index into coordinates)
	 * @param maxy
	 *            The maximum y value of the 2-dimensional data (used to unpack a linear index into coordinates)
	 */
	public SingleNBFreeCircularGaussian2DFunction(int maxx, int maxy)
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
		return new SingleNBFreeCircularGaussian2DFunction(maxx, maxy);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.function.gaussian.SingleFreeCircularGaussian2DFunction#eval(int, double[])
	 */
	public double eval(final int x, final double[] dyda)
	{
		// Unpack the predictor into the dimensions
		final int x1 = x / maxx;
		final int x0 = x % maxx;

		return background + gaussian(x0, x1, dyda);
	}

	private double gaussian(final int x0, final int x1, final double[] dy_da)
	{
		final double dx = x0 - x0pos;
		final double dy = x1 - x1pos;
		final double dx2 = dx * dx;
		final double dxy = dx * dy;
		final double dy2 = dy * dy;

		// Calculate gradients

		if (zeroAngle)
		{
			final double exp = FastMath.exp(aa * dx2 + cc * dy2);
			dy_da[0] = n * exp;
			final double y = height * exp;

			dy_da[1] = y * (-2.0 * aa * dx);
			dy_da[2] = y * (-2.0 * cc * dy);

			dy_da[3] = y * (nx + ax * dx2);
			dy_da[4] = y * (ny + cy * dy2);
			return y;
		}
		else
		{
			final double exp = FastMath.exp(aa * dx2 + bb * dxy + cc * dy2);
			dy_da[0] = n * exp;
			final double y = height * exp;

			dy_da[1] = y * (-2.0 * aa * dx - bb * dy);
			dy_da[2] = y * (-2.0 * cc * dy - bb * dx);

			dy_da[3] = y * (nx + ax * dx2 + bx * dxy + cx * dy2);
			dy_da[4] = y * (ny + ay * dx2 + by * dxy + cy * dy2);

			return y;
		}
	}

	@Override
	public boolean evaluatesBackground()
	{
		return false;
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
